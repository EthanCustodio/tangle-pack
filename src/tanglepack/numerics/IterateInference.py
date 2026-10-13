"""Filling the crossing iterate table, exactly.

A crossing found on two polylines is known only to within the two segments it
was found on: its unstable canonical distance lies in the BRACKET
``U(x) = [u(p0), u(p1)]`` of its unstable segment, its stable one in
``S(x) = [s(q0), s(q1)]`` of its stable segment. The machine builds every node by
mapping another one -- ``cdist(f(p)) = stretch_u * cdist(p)`` on an unstable
branch, ``cdist(f^-1(q)) = stretch_s * cdist(q)`` on a stable one, refined nodes
taking the mean of their neighbours, which commutes with the scaling -- so the
map carries the brackets EXACTLY::

    y = M(x)  iff  y lies on the advanced branches (FixedPoint.advance_key),
                   U(y) meets stretch_u * U(x)  and  S(y) meets S(x) / stretch_s,
                   and y is the only crossing that does;

and the real map then confirms the match: ``|f(x) - y| <= |f(p1) - f(p0)| +
|q1' - q0'|``, the image of x's unstable segment plus y's own, in ONE batched
map call per pass.

The rule runs in ``compute_intersections``, in every ``iterate_bridge`` (the
bridge's endpoints and the crossings just born) and at the end of every blast,
so the topology trusts the table and never matches identities by cdist itself.

Dev Notes:

* Why each manifold's own ``stretch_param`` and not
  :meth:`~.FixedPoint.FixedPoint.per_step_beta`: the stretch is the conjugacy the
  polylines were BUILT with, so the brackets are carried to float rounding.
  ``per_step_beta`` differs from it by 3e-7 (k=2.8) to 6e-6 (inversion)
  relative, more than a deep bracket is wide (median 3e-8 at nine k=2.8 blasts).
* Why the map is a gate and never the matcher: deep in a blast the lambda lemma
  puts distinct crossings within 3e-6 of each other in the plane (u = 1.1e3 and
  u = 1.6e9 at nine k=2.8 blasts); their parameter brackets keep them apart.
* The two brackets are compared INDIVIDUALLY. The cdist product is preserved
  along an iterate chain (area preservation), so it cannot tell two chains apart.
* Refined points cache only a preiterate (``ManifoldMachine``), so point-level
  ``next_iterate`` links cannot register crossings; the cdist parameter can.
* Each side advances on its OWN fixed point: on a heteroclinic crossing the
  unstable branch belongs to another orbit than the stable one.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import TYPE_CHECKING, Iterable, Optional

import numpy as np

from .Intersection import Intersection, ManifoldKey
from .IntersectionRegistry import IntersectionRegistry

if TYPE_CHECKING:  # pragma: no cover - import cycle guard
    from .TangleWorkbench import TangleWorkbench

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

#: Relative float slack of a bracket. Refined nodes are means of their
#: neighbours, not products, so a carried bracket meets its image only up to
#: rounding, and nodes tie to the ULP at folds.
BRACKET_PAD = 1e-12

Bracket = tuple[float, float]


def _brackets(crossing: Intersection) -> tuple[Bracket, Bracket]:
    """The crossing's unstable and stable brackets, padded by :data:`BRACKET_PAD`."""
    (p0, p1), (q0, q1) = crossing.unstable_segment, crossing.stable_segment
    unstable = sorted((float(p0.get_cdist("unstable")), float(p1.get_cdist("unstable"))))
    stable = sorted((float(q0.get_cdist("stable")), float(q1.get_cdist("stable"))))
    return tuple(
        (lo - BRACKET_PAD * abs(lo), hi + BRACKET_PAD * abs(hi)) for lo, hi in (unstable, stable)
    )


def _meets(a: Bracket, b: Bracket) -> bool:
    """Two closed intervals share a point."""
    return a[0] <= b[1] and b[0] <= a[1]


def infer_iterates(workbench: "TangleWorkbench", sources: Optional[Iterable[int]] = None) -> int:
    """
    Record the ``+1`` image of every source crossing the exact rule finds.

    An anchor's image is the anchor of the advanced branches. Any other
    crossing's image is the unique crossing whose brackets meet the source's
    carried one map step (module docstring), confirmed by the real map in one
    batched call. Idempotent: a source already linked is skipped.

    Only crossings of N are linked (B1): an unresolved detection is neither a
    source nor a candidate image, so no link starts or ends at one.

    Args:
        workbench: The workbench whose registry and iterate table this fills.
        sources: The crossings to find images for, searched against the
            crossings of N; None means every one of them.

    Returns:
        Number of new iterate relationships recorded.

    Note:
        A bracket match the map contradicts, or an image that already has a
        preimage, is a contradiction between parameter and geometry: it is
        logged at WARNING and not linked (never observed on the fixtures).
    """
    registry = workbench.intersection_registry
    table = registry.iterate_table
    resolved = registry.resolved_ids(workbench.resolution_ratio_threshold)
    ids = registry.all_ids() if sources is None else dict.fromkeys(sources)
    pending = [i for i in ids if i in resolved and (i, 1) not in table]

    recorded = sum(
        _register_anchor_iterate(registry, i, registry[i])
        for i in pending
        if registry[i].is_synthetic
    )
    detected = [i for i in pending if not registry[i].is_synthetic]
    if not detected:
        return recorded

    widest = _widest_stable_brackets(registry)
    matches = {
        i: j
        for i in detected
        if (j := _bracketed_image(workbench, registry, i, widest, resolved)) is not None
    }
    for i, j in _map_gated(workbench, registry, matches).items():
        if (j, -1) in table:
            logger.warning(
                "crossing %d: its image %d already has the preimage %d; not linked",
                i,
                j,
                table[j, -1],
            )
            continue
        registry.register_iterate(i, 1, j)
        recorded += 1
    return recorded


def _register_anchor_iterate(
    registry: IntersectionRegistry, src_id: int, source: Intersection
) -> bool:
    """
    Record the forward image of an ANCHOR -- a periodic point's own crossing.

    A periodic point sits at canonical distance ``(0, 0)`` on every branch
    attached to it, and the map carries it to the next point of its orbit,
    whose anchor sits at ``(0, 0)`` on the ADVANCED branches. So the image
    is fully determined by the branch bookkeeping and
    :meth:`IntersectionRegistry.find` resolves it exactly (at ``(0, 0)`` it
    matches the point, whichever of its branches the probe names). On a
    period-1 orbit the anchor is its own image.

    Args:
        registry: The registry to read and write.
        src_id: Registry id of the crossing being iterated.
        source: The crossing itself.

    Returns:
        True if the source is an anchor and its image was recorded, else
        False (including for an anchor whose image branch has no anchor yet,
        which happens while a period > 1 orbit is only partly initialized).
    """
    if source.unstable_cdist != 0.0 or source.stable_cdist != 0.0:
        return False
    a_key, b_key = (
        key[0].advance_key(key, 1) for key in (source.manifold_a_key, source.manifold_b_key)
    )
    target = registry.find(
        Intersection.synthetic(
            coords=source.coords,
            unstable_cdist=0.0,
            stable_cdist=0.0,
            manifold_a_key=a_key,
            manifold_b_key=b_key,
        )
    )
    if target is None:
        logger.debug(
            "Anchor %d has no image anchor on branches (%s, %s) yet", src_id, a_key, b_key
        )
        return False
    registry.register_iterate(src_id, 1, target)
    return True


def _widest_stable_brackets(registry: IntersectionRegistry) -> dict[ManifoldKey, float]:
    """Per stable branch, the widest stable bracket of a detected crossing on it.

    A candidate image's stable cdist lies inside its own bracket, so the
    candidates of a predicted stable bracket are those within this width of it.
    """
    widest: dict[ManifoldKey, float] = defaultdict(float)
    for _, crossing in registry:
        if not crossing.is_synthetic:
            lo, hi = _brackets(crossing)[1]
            widest[crossing.manifold_b_key] = max(widest[crossing.manifold_b_key], hi - lo)
    return widest


def _bracketed_image(
    workbench: "TangleWorkbench",
    registry: IntersectionRegistry,
    src_id: int,
    widest: dict[ManifoldKey, float],
    resolved: frozenset[int],
) -> Optional[int]:
    """The unique crossing of N whose brackets meet the source's carried one step, or None."""
    source = registry[src_id]
    a_key, b_key = (
        key[0].advance_key(key, 1) for key in (source.manifold_a_key, source.manifold_b_key)
    )
    stretch_u = source.unstable_manifold.stretch_param
    stretch_s = workbench.manifolds[source.manifold_b_key].stretch_param
    (u_lo, u_hi), (s_lo, s_hi) = _brackets(source)
    unstable, stable = (u_lo * stretch_u, u_hi * stretch_u), (s_lo / stretch_s, s_hi / stretch_s)

    candidates = registry.candidates_near_stable_cdist(
        b_key, 0.5 * (stable[0] + stable[1]), 0.5 * (stable[1] - stable[0]) + widest[b_key]
    )
    images = [
        j
        for j, y in candidates
        if j in resolved
        and not y.is_synthetic
        and y.manifold_a_key == a_key
        and _meets(unstable, _brackets(y)[0])
        and _meets(stable, _brackets(y)[1])
    ]
    if len(images) > 1:
        logger.info(
            "crossing %d: %d crossings %s meet its carried brackets (under-resolved); "
            "not linked",
            src_id,
            len(images),
            images,
        )
        return None
    return images[0] if images else None


def _map_gated(
    workbench: "TangleWorkbench", registry: IntersectionRegistry, matches: dict[int, int]
) -> dict[int, int]:
    """The bracket matches the real map confirms, in ONE batched map call.

    A match ``x -> y`` stands when ``|f(x) - y|`` is at most the image of x's
    unstable segment plus y's own unstable segment: the scale at which the two
    polylines resolve the crossings.
    """
    if not matches:
        return {}
    sources = [registry[i] for i in matches]
    points = np.array(
        [
            np.ravel(point)
            for x in sources
            for point in (x.coords, *(node.get_point() for node in x.unstable_segment))
        ],
        dtype=float,
    )
    images = workbench.dynamical_system.map_batch(points).reshape(len(sources), 3, 2)

    confirmed: dict[int, int] = {}
    for (i, j), (fx, fp0, fp1) in zip(matches.items(), images):
        y = registry[j]
        q0, q1 = (np.ravel(node.get_point()) for node in y.unstable_segment)
        gate = np.linalg.norm(fp1 - fp0) + np.linalg.norm(q1 - q0)
        miss = float(np.linalg.norm(fx - np.ravel(y.coords)))
        if miss <= gate:
            confirmed[i] = j
        else:
            logger.warning(
                "crossing %d: the brackets name %d as its image, but the map "
                "misses it by %.3g (gate %.3g); not linked",
                i,
                j,
                miss,
                gate,
            )
    return confirmed
