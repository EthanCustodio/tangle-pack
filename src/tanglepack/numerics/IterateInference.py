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

* The brackets are read and carried PRECISELY (:mod:`~.precise`): each end is
  its node's cdist with its residual, carried by ``dd_mul`` / ``dd_div``, and
  the meet test compares the precise values, so a deep bracket narrower than a
  float's ULP still meets only what it meets.
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

import bisect
import logging
from collections import defaultdict
from typing import TYPE_CHECKING, Iterable, Optional

import numpy as np
from numpy.typing import NDArray

from .Intersection import Intersection, ManifoldKey
from .IntersectionRegistry import IntersectionRegistry
from .precise import compare, dd_div, dd_mul, two_sum

if TYPE_CHECKING:  # pragma: no cover - import cycle guard
    from .TangleWorkbench import TangleWorkbench

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

#: Relative slack of a bracket: the error of the precise carry. A node's cdist
#: and its image's are one double-double computation apart (a product by the
#: stretch, a mean of neighbours), exact to about 1e-32 relative; 1e-24 is far
#: above that and far below a deep bracket's width. It was 1e-12 on plain
#: floats (2026-10-30: the identity probe is unchanged by the switch).
BRACKET_PAD = 1e-24

#: A bracket row as precise pairs: ``(hi, lo)`` arrays of ``(u_lo, u_hi, s_lo, s_hi)``.
PreciseBrackets = tuple[NDArray[np.float64], NDArray[np.float64]]


def _padded(hi: NDArray[np.float64], lo: NDArray[np.float64]) -> PreciseBrackets:
    """Precise bracket ends ``(u_lo, u_hi, s_lo, s_hi)`` (last axis), each side
    widened by :data:`BRACKET_PAD` of its ends, lower ends down and upper ends up."""
    sign = np.where(np.arange(4) % 2 == 0, -1.0, 1.0)
    s, e = two_sum(hi, sign * BRACKET_PAD * np.abs(hi))
    return s, e + lo


def _brackets(crossing: Intersection) -> PreciseBrackets:
    """The crossing's precise brackets ``(u_lo, u_hi, s_lo, s_hi)``, padded by :data:`BRACKET_PAD`."""
    return _padded(*crossing.precise_brackets)


def infer_iterates(workbench: "TangleWorkbench", sources: Optional[Iterable[int]] = None) -> int:
    """
    Record the ``+1`` image of every source crossing the exact rule finds.

    An anchor's image is the anchor of the advanced branches. Any other
    crossing's image is the unique crossing whose brackets meet the source's
    carried one map step (module docstring), confirmed by the real map in one
    batched call. Idempotent: a source already linked is skipped.

    Only crossings of N are linked: an unresolved detection is neither a
    source nor a candidate image, so no link starts or ends at one.

    What the brackets cannot decide is recorded, never guessed: a source
    whose brackets meet several images, and every claimant of an image named
    by two sources (this pass's, the one already linked to it, or one an
    earlier pass left contested), go to ``registry.unresolved_links`` and none
    of them is linked. The map is a bijection, so at most one claimant is
    right, and no order of the passes may pick it. Each pass re-decides its
    own sources' entries.

    Args:
        workbench: The workbench whose registry and iterate table this fills.
        sources: The crossings to find images for, searched against the
            crossings of N; None means every one of them.

    Returns:
        Number of new iterate relationships recorded.

    Note:
        A bracket match the map contradicts, or an image claimed by several
        crossings, is a contradiction between parameter and geometry: it is
        logged at WARNING and not linked (never observed on the fixtures).
    """
    registry = workbench.intersection_registry
    table = registry.iterate_table
    resolved = registry.resolved_ids(workbench.resolution_ratio_threshold)
    ids = registry.all_ids() if sources is None else dict.fromkeys(sources)
    before = dict(registry.unresolved_links)
    for i in ids:
        registry.unresolved_links.pop(i, None)
    pending = [i for i in ids if i in resolved and (i, 1) not in table]

    recorded = sum(
        _register_anchor_iterate(registry, i, registry[i])
        for i in pending
        if registry[i].is_synthetic
    )
    detected = [i for i in pending if not registry[i].is_synthetic]
    matches: dict[int, int] = {}
    if detected:
        widest = _widest_stable_brackets(registry)
        matches = {
            i: j
            for i in detected
            if (j := _bracketed_image(workbench, registry, i, widest, resolved)) is not None
        }

    claimed_by: dict[int, list[int]] = defaultdict(list)
    for i, j in _map_gated(workbench, registry, matches).items():
        claimed_by[j].append(i)
    for j, claimants in claimed_by.items():
        held = table[j, -1]
        earlier = [s for s, images in registry.unresolved_links.items() if images == (j,)]
        contested = [*claimants, *([] if held is None else [held]), *earlier]
        if len(contested) == 1:
            registry.register_iterate(contested[0], 1, j)
            recorded += 1
        else:
            logger.warning("image %d is claimed by crossings %s; none linked", j, contested)
            if held is not None:
                registry.remove_iterate(held, 1)
            registry.unresolved_links.update(dict.fromkeys(contested, (j,)))
    if registry.unresolved_links != before:
        registry._generation += 1
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
    for branch_key, (_cdists, _ids, hi, lo) in registry._stable_branch_buckets().items():
        stable = _padded(hi, lo)[0][:, 2:]
        # fmax skips the NaN rows of crossings without brackets; 0.0 when none has one.
        widest[branch_key] = float(np.fmax.reduce(stable[:, 1] - stable[:, 0], initial=0.0))
    return widest


def _bracketed_image(
    workbench: "TangleWorkbench",
    registry: IntersectionRegistry,
    src_id: int,
    widest: dict[ManifoldKey, float],
    resolved: frozenset[int],
) -> Optional[int]:
    """The unique crossing of N whose brackets meet the source's carried one step, or None.

    Several such crossings are recorded in ``registry.unresolved_links``.
    """
    source = registry[src_id]
    a_key, b_key = (
        key[0].advance_key(key, 1) for key in (source.manifold_a_key, source.manifold_b_key)
    )
    stretch_u = source.unstable_manifold.stretch_param
    stretch_s = workbench.manifolds[source.manifold_b_key].stretch_param
    hi, lo = _brackets(source)
    (u_hi, u_lo), (s_hi, s_lo) = dd_mul(hi[:2], lo[:2], stretch_u), dd_div(hi[2:], lo[2:], stretch_s)
    carried = np.concatenate((u_hi, s_hi)), np.concatenate((u_lo, s_lo))
    stable = s_hi

    # A candidate's stable cdist lies within the predicted bracket widened by the
    # branch's widest stable bracket.
    cdists, ids, near_hi, near_lo = registry.stable_brackets(b_key)
    centre, window = 0.5 * (stable[0] + stable[1]), 0.5 * (stable[1] - stable[0]) + widest[b_key]
    first = bisect.bisect_left(cdists, centre - window)
    last = bisect.bisect_right(cdists, centre + window)
    near = _padded(near_hi[first:last], near_lo[first:last])  # NaN rows (no brackets) meet nothing
    # Each lower end of one bracket at most the upper end of the other, precisely.
    meets = np.all(
        (compare(*(c[[0, 2]] for c in carried), *(n[:, [1, 3]] for n in near)) <= 0)
        & (compare(*(n[:, [0, 2]] for n in near), *(c[[1, 3]] for c in carried)) <= 0),
        axis=1,
    )
    images = [
        j
        for j in ids[first:last][meets].tolist()
        if j in resolved and registry[j].manifold_a_key == a_key
    ]
    if len(images) > 1:
        registry.unresolved_links[src_id] = tuple(images)
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
            for point in (x.coords, *x.segment_ends("unstable")[0])
        ],
        dtype=float,
    )
    images = workbench.dynamical_system.map_batch(points).reshape(len(sources), 3, 2)

    confirmed: dict[int, int] = {}
    for (i, j), (fx, fp0, fp1) in zip(matches.items(), images):
        y = registry[j]
        q0, q1 = y.segment_ends("unstable")[0]
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
