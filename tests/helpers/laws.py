"""The physical-law checks of the law tier, one function per law.

Every check has the signature ``check_<law>(case) -> int``: it asserts its law
on a built :class:`cases.Case` and returns HOW MANY ITEMS it checked, so the
law tier can insist that a law was not vacuous (a count of 0 is a failure
unless the ``(law, case)`` pair is listed in ``cases.NOT_APPLICABLE``). The
checks only read the case; they never mutate the session (the law tier shares
one build per module and fingerprints it around every test). The two laws
that must mutate a session (a recompute) share ONE fresh build per case, the
``recompute`` layer of :data:`MUTATING_LAYERS`:
:func:`check_recompute_preserves_ids` then :func:`check_recompute_is_idempotent`.

A check reads the session's cached products (``session.trellis()``,
``arrangement()``, ``bridge_classes()``, ``minimal_trellis()``,
``iterated_partition()``, ``dual_graph()``, ``symbolic_dynamics()``); the
law-tier fixture builds every one of them before the first test, so a check
never triggers a build.

The checks are grouped into :data:`LAYERS` (and :data:`MUTATING_LAYERS`); the
law tier runs ONE test per (layer, case). The law id of a check (the key of
``cases.KNOWN_ISSUES`` / ``cases.NOT_APPLICABLE``) is its name without the
``check_`` prefix (:func:`law_id`).

Dev Notes:

* Tolerances are stated against library constants where the library has one:
  the area of an iterate link against ``collision_rtol`` (the slack
  ``StrongPip`` / ``Pseudoneighbor`` identify iterates with), the cdist
  scaling of a link against ``infer_iterate_table``'s ``cdist_rtol`` (the
  slack the links were matched with), the area of a region image against
  ``Arrangement.image_of``'s ``rtol``, and every cdist comparison against the
  registry's ``cdist_tol``.
* The geometric oracles here (crossing sign from the manifold nodes, the row
  of a bridge end, a direct hole's side of its bridge, the side of a face) are
  measured from the curves with public geometry helpers and the ``_side_of``
  kernel, independently of the combinatorial rule the library uses.
"""

from __future__ import annotations

import inspect
from collections import defaultdict
from typing import TYPE_CHECKING, Callable, Optional

import numpy as np
import pytest

from helpers.invariants import (
    assert_cdist_monotonic,
    assert_iterate_relation,
    assert_no_geometric_spikes,
    assert_one_to_one,
)
from tanglepack.numerics.IterateInference import IterateInference
from tanglepack.numerics.geometry import oriented_bridge_polyline
from tanglepack.topology.Arrangement import Arrangement
from tanglepack.topology.BridgeClass import (
    anchor_outward_key,
    class_sort_key,
    oriented_class,
)
from tanglepack.topology.MinimalTrellis import image_chain
from tanglepack.topology.Pseudoneighbor import compute_pseudoneighbors
from tanglepack.topology.StablePartition import (
    _side_of,
    bridge_for_pair,
    bridge_row_violation,
    owns_cdist,
    row_of_end,
    rows_of_bridge,
)
if TYPE_CHECKING:  # pragma: no cover - typing only
    from cases import Case

#: Relative slack of the per-link area check: the library's own iterate
#: identification slack (``collision_rtol`` default of ``StrongPip`` and
#: ``Pseudoneighbor``). Measured worst link at baseline: 1.4e-3 (p3, nested),
#: 2.3e-4 (k28 two blasts), 1.8e-4 (k10).
AREA_LINK_RTOL: float = inspect.signature(compute_pseudoneighbors).parameters[
    "collision_rtol"
].default

#: Relative slack of a ``+1`` link's cdist ratio against ``per_step_beta``:
#: the ``cdist_rtol`` the iterate table's links are matched with.
LINK_SCALING_RTOL: float = inspect.signature(
    IterateInference.infer_iterate_table
).parameters["cdist_rtol"].default

#: Relative slack of a region image's area (``Arrangement.image_of``'s ``rtol``).
REGION_AREA_RTOL: float = inspect.signature(Arrangement.image_of).parameters["rtol"].default

#: A law check: asserts on a case and returns the number of items checked.
LawCheck = Callable[["Case"], int]


def law_id(check: LawCheck) -> str:
    """The law id of a check: its name without the ``check_`` prefix."""
    return check.__name__.removeprefix("check_")


# --------------------------------------------------------------------------- #
# Small readers
# --------------------------------------------------------------------------- #
def _anchor_ids(case: "Case") -> list[int]:
    """Registry ids of the synthetic anchors (both cdists exactly 0, no segments)."""
    return [
        intersection_id
        for intersection_id, ix in case.registry
        if ix.seg_ids is None
        and ix.unstable_cdist is not None
        and ix.stable_cdist is not None
        and float(ix.unstable_cdist) == 0.0
        and float(ix.stable_cdist) == 0.0
    ]


def _tol(case: "Case") -> float:
    """The registry's cdist tolerance."""
    return float(case.registry.cdist_tol)


def _full(case: "Case"):
    """The all-fixed-points trellis (cached by the session)."""
    return case.session.trellis()


def _per_fp_trellises(case: "Case") -> list:
    """The per-fixed-point trellises, outermost fixed point first."""
    return [case.session.trellis(fp) for fp in case.fixed_points]


def _holes(case: "Case") -> list:
    """``(trellis, hole)`` for every hole of every per-fixed-point trellis."""
    return [(trellis, hole) for trellis in _per_fp_trellises(case) for hole in trellis.holes]


def _partitions(case: "Case") -> list:
    """``(trellis, result)`` for every stable partition of every per-fixed-point trellis."""
    return [
        (trellis, result)
        for trellis in _per_fp_trellises(case)
        for result in trellis.stable_partitions
    ]


def _identified_bridges(case: "Case") -> list:
    """The workbench's bridges that carry a ``BridgeId``."""
    return [bridge for bridge in case.workbench.bridges if bridge.id is not None]


def _stable_nodes(manifold) -> tuple[np.ndarray, np.ndarray]:
    """``(points, cdists)`` of a stable manifold, root (anchor) first."""
    nodes = manifold.get_point_array(return_nodes=True)
    points = np.array([np.asarray(n.get_point(), dtype=float).ravel() for n in nodes])
    cdists = np.array(
        [np.nan if n.get_cdist("stable") is None else float(n.get_cdist("stable")) for n in nodes]
    )
    return points, cdists


def _local_direction(manifold, stability: str, cdist: float) -> np.ndarray:
    """The curve direction at ``cdist``, taken in increasing canonical distance."""
    nodes = manifold.get_point_array(return_nodes=True)
    cdists = [n.get_cdist(stability) for n in nodes]
    index = int(np.searchsorted(cdists, cdist))
    index = min(max(index, 1), len(nodes) - 1)
    return np.asarray(nodes[index].get_point(), dtype=float).ravel() - np.asarray(
        nodes[index - 1].get_point(), dtype=float
    ).ravel()


def _anchorward_look(case: "Case", intersection_id: int) -> Optional[np.ndarray]:
    """The stable dynamical direction (toward the anchor) at a crossing, from the nodes."""
    ix = case.registry[intersection_id]
    manifold = case.workbench.manifolds.get(ix.manifold_b_key)
    if manifold is None:
        return None
    points, cdists = _stable_nodes(manifold)
    here = np.asarray(ix.coords, dtype=float)
    tol = _tol(case)
    below = [i for i in range(len(cdists)) if cdists[i] < float(ix.stable_cdist) - tol]
    above = [i for i in range(len(cdists)) if cdists[i] > float(ix.stable_cdist) + tol]
    if below:
        return points[max(below, key=lambda i: cdists[i])] - here
    if above:
        return here - points[min(above, key=lambda i: cdists[i])]
    return None


def _bridge_end_displacement(
    case: "Case", bridge, intersection_id: int
) -> Optional[np.ndarray]:
    """From a crossing INTO the bridge: the second clearly distinct node at that end."""
    registry = case.registry
    first, second = bridge.id
    poly = oriented_bridge_polyline(
        bridge, float(registry[first].unstable_cdist), float(registry[second].unstable_cdist)
    )
    if poly is None or len(poly) < 3:
        return None
    here = np.asarray(registry[intersection_id].coords, dtype=float)
    chord = float(np.linalg.norm(poly[-1] - poly[0]))
    eps = 1e-9 * (1.0 + chord)
    order = range(1, len(poly)) if intersection_id == first else range(len(poly) - 2, -1, -1)
    found = [i for i in order if float(np.linalg.norm(poly[i] - here)) > eps][:2]
    if not found:
        return None
    return poly[found[-1]] - here


def _geometric_row(case: "Case", bridge, intersection_id: int) -> Optional[str]:
    """The stable side a bridge approaches a crossing from, measured on the curves."""
    look = _anchorward_look(case, intersection_id)
    displacement = _bridge_end_displacement(case, bridge, intersection_id)
    if look is None or displacement is None:
        return None
    return _side_of(look, displacement)


def _point_side_of_bridge(case: "Case", bridge, point: np.ndarray) -> Optional[str]:
    """The side of a bridge (in its unstable dynamical direction) a point sits on.

    Measured against the nearest segment of the oriented polyline.
    """
    registry = case.registry
    first, second = bridge.id
    poly = oriented_bridge_polyline(
        bridge, float(registry[first].unstable_cdist), float(registry[second].unstable_cdist)
    )
    if poly is None or len(poly) < 2:
        return None
    a, b = poly[:-1], poly[1:]
    ab = b - a
    length2 = np.maximum(np.einsum("ij,ij->i", ab, ab), 1e-300)
    t = np.clip(np.einsum("ij,ij->i", point - a, ab) / length2, 0.0, 1.0)
    projected = a + t[:, None] * ab
    nearest = int(np.argmin(np.linalg.norm(projected - point, axis=1)))
    return _side_of(ab[nearest], point - projected[nearest])


def _within(value: float, lo: float, hi: float, tol: float) -> bool:
    """``lo - tol <= value <= hi + tol``."""
    return lo - tol <= value <= hi + tol


# --------------------------------------------------------------------------- #
# Case sanity
# --------------------------------------------------------------------------- #
def check_case_builds(case: "Case") -> int:
    """The case runs the full pipeline: every product exists and is non-empty."""
    session = case.session
    assert len(case.fixed_points) == case.expect.n_tangles
    assert tuple(fp.period for fp in case.fixed_points) == case.expect.periods
    for trellis in _per_fp_trellises(case):
        assert trellis.stable_partitions, "every fixed point is partitioned"
    products = [
        session.arrangement(),
        session.bridge_classes(),
        session.minimal_trellis(),
        session.iterated_partition(),
        session.dual_graph(),
        session.symbolic_dynamics(),
    ]
    assert len(products[1]) > 0, "the case has bridge classes"
    assert products[5].classes, "the symbolic dynamics covers the classes"
    return len(products)


def check_k_value_matches_inversion(case: "Case") -> int:
    """``k_value = 2 * period`` exactly on an inversion point (both eigenvalues < 0).

    An inversion point here is orientation preserving: ``det J = lambda_u *
    lambda_s = +1``; it has two branches per stability.
    """
    for fp in case.fixed_points:
        lambda_u = float(np.asarray(fp.unstable_eigenvalues[0]).ravel()[0])
        lambda_s = float(np.asarray(fp.stable_eigenvalues[0]).ravel()[0])
        inversion = bool(fp.check_inversion())
        assert inversion == (lambda_u < 0 and lambda_s < 0)
        assert fp.k_value == (2 if inversion else 1) * fp.period
        assert fp.num_branches == (2 if inversion else 1)
        assert lambda_u * lambda_s == pytest.approx(1.0)
    assert any(fp.check_inversion() for fp in case.fixed_points) == case.expect.inversion
    return len(case.fixed_points)


def check_every_branch_is_built(case: "Case") -> int:
    """Every key of every branch cycle is a manifold of the workbench, keyed by itself."""
    count = 0
    for fp in case.fixed_points:
        for stability in ("unstable", "stable"):
            for key in fp.branch_cycle(stability):
                manifold = case.workbench.manifolds.get(key)
                assert manifold is not None, f"branch {key[1:]} was never built"
                assert manifold.manifold_key == key
                assert len(manifold.get_point_array()) > 5, f"branch {key[1:]} never grew"
                count += 1
    return count


def check_every_branch_crosses(case: "Case") -> int:
    """Every unstable and every stable branch of every cycle carries a transverse
    (non-anchor) crossing (both branches of an inversion point included)."""
    anchors = set(_anchor_ids(case))
    seen = set()
    for intersection_id, ix in case.registry:
        if intersection_id not in anchors:
            seen.update((ix.manifold_a_key, ix.manifold_b_key))
    count = 0
    for fp in case.fixed_points:
        for stability in ("unstable", "stable"):
            for key in fp.branch_cycle(stability):
                assert key in seen, f"{stability} branch {key[1:]} crosses nothing"
                count += 1
    return count


# --------------------------------------------------------------------------- #
# Manifolds
# --------------------------------------------------------------------------- #
def check_manifold_cdist_monotone(case: "Case") -> int:
    """cdist is non-decreasing along every manifold and every bridge (ties allowed)."""
    curves = list(case.workbench.manifolds.values()) + list(case.workbench.bridges)
    for curve in curves:
        assert_cdist_monotonic(curve, strict=False)
    return len(curves)


def check_manifold_no_spikes(case: "Case") -> int:
    """No manifold carries a geometric spike (a scrambled node)."""
    manifolds = list(case.workbench.manifolds.values())
    for manifold in manifolds:
        assert_no_geometric_spikes(manifold)
    return len(manifolds)


def check_manifold_iterate_law(case: "Case") -> int:
    """``c_iterate = stretch_param * c`` along every manifold's iterate chain."""
    manifolds = list(case.workbench.manifolds.values())
    for manifold in manifolds:
        assert_iterate_relation(manifold)
    return len(manifolds)


def check_manifold_one_to_one(case: "Case") -> int:
    """The geometric and iterate lists of every manifold and bridge are acyclic and consistent."""
    curves = list(case.workbench.manifolds.values()) + list(case.workbench.bridges)
    for curve in curves:
        assert_one_to_one(curve)
    return len(curves)


# --------------------------------------------------------------------------- #
# Crossings
# --------------------------------------------------------------------------- #
def check_crossings_unstable_by_stable(case: "Case") -> int:
    """Every registered crossing is one unstable key (a) times one stable key (b)."""
    count = 0
    for intersection_id, ix in case.registry:
        assert ix.manifold_a_key is not None and ix.manifold_b_key is not None, intersection_id
        assert ix.manifold_a_key[1] == "unstable", (intersection_id, ix.manifold_a_key)
        assert ix.manifold_b_key[1] == "stable", (intersection_id, ix.manifold_b_key)
        count += 1
    return count


def _transversal_hits(polylines: list[np.ndarray], lo: np.ndarray, hi: np.ndarray) -> list:
    """Transversal interior crossings among the segments of the polylines inside a box.

    Two segments of ONE polyline that share a vertex never count; grazing
    contacts and shared endpoints are excluded by an interior parameter margin.
    """
    starts, ends, owners, positions = [], [], [], []
    for owner, points in enumerate(polylines):
        if len(points) < 2:
            continue
        a, b = points[:-1], points[1:]
        keep = np.all(np.maximum(a, b) >= lo, axis=1) & np.all(np.minimum(a, b) <= hi, axis=1)
        index = np.nonzero(keep)[0]
        starts.append(a[index])
        ends.append(b[index])
        owners.append(np.full(len(index), owner))
        positions.append(index)
    if not starts:
        return []
    a = np.concatenate(starts)
    b = np.concatenate(ends)
    owner = np.concatenate(owners)
    position = np.concatenate(positions)
    cell = max(float((hi - lo).max()) / 256.0, 1e-12)
    cmin = np.floor((np.minimum(a, b) - lo) / cell).astype(int)
    cmax = np.floor((np.maximum(a, b) - lo) / cell).astype(int)
    buckets: dict[tuple[int, int], list[int]] = defaultdict(list)
    for s in range(len(a)):
        for cx in range(cmin[s, 0], cmax[s, 0] + 1):
            for cy in range(cmin[s, 1], cmax[s, 1] + 1):
                buckets[(cx, cy)].append(s)
    pairs = {
        (members[i], members[j])
        for members in buckets.values()
        for i in range(len(members))
        for j in range(i + 1, len(members))
    }
    if not pairs:
        return []
    pair_array = np.array(sorted(pairs))
    i, j = pair_array[:, 0], pair_array[:, 1]
    adjacent = (owner[i] == owner[j]) & (np.abs(position[i] - position[j]) < 2)
    i, j = i[~adjacent], j[~adjacent]
    d1, d2 = b[i] - a[i], b[j] - a[j]
    den = d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]
    r = a[j] - a[i]
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (r[:, 0] * d2[:, 1] - r[:, 1] * d2[:, 0]) / den
        u = (r[:, 0] * d1[:, 1] - r[:, 1] * d1[:, 0]) / den
    margin = 1e-9
    hit = (np.abs(den) > 1e-15) & (t > margin) & (t < 1 - margin) & (u > margin) & (u < 1 - margin)
    return list(zip(owner[i[hit]], position[i[hit]], owner[j[hit]], position[j[hit]]))


def _unstable_image_curves(case: "Case") -> list:
    """The blast's image bridges: unstable curves that are not pieces of a manifold.

    A bridge cut from a manifold shares its nodes with it; an iterated image
    bridge has nodes of its own.
    """
    manifold_nodes = {
        id(node)
        for key, manifold in case.workbench.manifolds.items()
        if key[1] == "unstable"
        for node in manifold.get_point_array(return_nodes=True)
    }
    return [
        bridge
        for bridge in case.workbench.bridges
        if not any(id(node) in manifold_nodes for node in bridge.get_point_array(return_nodes=True))
    ]


def check_no_same_stability_crossing(case: "Case") -> int:
    """No two unstable curves cross, and no two stable curves cross (CLAUDE.md law).

    Checked geometrically on the computed polylines -- every manifold of every
    fixed point plus the blast's image bridges -- one curve with itself, two
    branches, two fixed points, inside the box of the registered crossings
    padded by 10%: the tangle region. Outside it the binary horseshoe's
    escaping arm is under-resolved and says nothing about the dynamics.
    """
    points = np.array([ix.coords for _iid, ix in case.registry], dtype=float)
    lo, hi = points.min(axis=0), points.max(axis=0)
    pad = 0.1 * float((hi - lo).max())
    lo, hi = lo - pad, hi + pad
    count = 0
    for stability in ("unstable", "stable"):
        curves = [
            (str(key[1:]), manifold)
            for key, manifold in case.workbench.manifolds.items()
            if key[1] == stability
        ]
        if stability == "unstable":
            curves += [(f"image bridge {b.id}", b) for b in _unstable_image_curves(case)]
        polylines = [np.asarray(curve.get_point_array(), dtype=float) for _name, curve in curves]
        hits = _transversal_hits(polylines, lo, hi)
        assert not hits, (
            f"{len(hits)} {stability} x {stability} crossing(s), first between "
            f"{curves[hits[0][0]][0]} segment {hits[0][1]} and "
            f"{curves[hits[0][2]][0]} segment {hits[0][3]}"
        )
        count += len(curves)
    return count


def check_crossings_cdists_defined(case: "Case") -> int:
    """Every registered crossing has two non-negative canonical distances."""
    count = 0
    for intersection_id, ix in case.registry:
        assert ix.unstable_cdist is not None and float(ix.unstable_cdist) >= 0.0, intersection_id
        assert ix.stable_cdist is not None and float(ix.stable_cdist) >= 0.0, intersection_id
        count += 1
    return count


def check_crossing_cdist_bracketed(case: "Case") -> int:
    """A crossing's unstable cdist lies between the cdists of the segment it was found on."""
    tol = _tol(case)
    count = 0
    for intersection_id, ix in case.registry:
        if ix.seg_ids is None:
            continue  # a synthetic anchor
        assert ix.unstable_segment is not None, intersection_id
        p0, p1 = ix.unstable_segment
        lo, hi = sorted((float(p0.get_cdist("unstable")), float(p1.get_cdist("unstable"))))
        assert _within(float(ix.unstable_cdist), lo, hi, tol), (
            f"crossing {intersection_id} at unstable cdist {ix.unstable_cdist} "
            f"outside its segment [{lo}, {hi}]"
        )
        count += 1
    return count


def _unstable_curve_at(case: "Case", intersection_id: int):
    """The unstable curve a crossing lies on: a bridge ending at it, else its manifold.

    A blast registers crossings on ITERATED bridges, which are images, not
    pieces of the grown manifold, so the manifold nodes do not pass through
    them; the bridge that ends at the crossing does.
    """
    for bridge_id in case.workbench.bridges_at(intersection_id):
        return case.workbench.bridge(bridge_id)
    return case.workbench.manifolds[case.registry[intersection_id].manifold_a_key]


def check_crossing_sign_is_cross_product(case: "Case") -> int:
    """``crossing_sign = sign(cross(u+, s+))``, rebuilt from the curve nodes.

    The directions are taken in increasing canonical distance from the nodes
    flanking the crossing on the unstable curve it lies on and on its stable
    manifold. Anchors are excluded (their rays are the eigendirections; see
    the anchor law).
    """
    count = 0
    for intersection_id, ix in case.registry:
        if ix.seg_ids is None:
            continue
        u_dir = _local_direction(
            _unstable_curve_at(case, intersection_id), "unstable", ix.unstable_cdist
        )
        s_dir = _local_direction(
            case.workbench.manifolds[ix.manifold_b_key], "stable", ix.stable_cdist
        )
        expected = int(np.sign(u_dir[0] * s_dir[1] - u_dir[1] * s_dir[0]))
        assert ix.crossing_sign == expected, (
            f"crossing {intersection_id}: sign {ix.crossing_sign}, cross product {expected}"
        )
        count += 1
    return count


def check_area_along_iterate_links(case: "Case", *, rtol: float = AREA_LINK_RTOL) -> int:
    """``unstable_cdist * stable_cdist`` is preserved by every registered ``+1`` link.

    Compared LINK BY LINK, not against the head of the chain, so the error of
    a long chain (period 3 runs a dozen links) does not accumulate. The
    anchors (product 0) are skipped.
    """
    registry = case.registry
    anchors = set(_anchor_ids(case))
    count = 0
    for intersection_id, ix in registry:
        image_id = registry.iterate_table[intersection_id, 1]
        if image_id is None or intersection_id in anchors:
            continue
        image = registry[image_id]
        before = float(ix.unstable_cdist) * float(ix.stable_cdist)
        after = float(image.unstable_cdist) * float(image.stable_cdist)
        assert after == pytest.approx(before, rel=rtol), (
            f"area not preserved on the link {intersection_id} -> {image_id}: "
            f"{before!r} -> {after!r}"
        )
        count += 1
    return count


def check_iterate_link_scales_by_beta(case: "Case") -> int:
    """A ``+1`` link scales the unstable cdist by ``per_step_beta`` and the stable by its inverse,
    and lands on the branch ``advance_key`` names."""
    registry = case.registry
    anchors = set(_anchor_ids(case))
    count = 0
    for intersection_id, ix in registry:
        image_id = registry.iterate_table[intersection_id, 1]
        if image_id is None or intersection_id in anchors:
            continue
        image = registry[image_id]
        fp = ix.manifold_a_key[0]
        beta = fp.per_step_beta("unstable")
        assert image.manifold_a_key == fp.advance_key(ix.manifold_a_key, 1), intersection_id
        assert image.manifold_b_key == ix.manifold_b_key[0].advance_key(ix.manifold_b_key, 1)
        assert float(image.unstable_cdist) == pytest.approx(
            beta * float(ix.unstable_cdist), rel=LINK_SCALING_RTOL
        ), (intersection_id, image_id)
        assert float(image.stable_cdist) == pytest.approx(
            float(ix.stable_cdist) / beta, rel=LINK_SCALING_RTOL
        ), (intersection_id, image_id)
        count += 1
    return count


def _crossing_fingerprint(case: "Case") -> dict[int, tuple]:
    """``{id: (keys, cdists)}`` of every registered crossing."""
    return {
        intersection_id: (
            ix.manifold_a_key,
            ix.manifold_b_key,
            float(ix.unstable_cdist),
            float(ix.stable_cdist),
        )
        for intersection_id, ix in case.registry
    }


def _same_crossing(a: tuple, b: tuple, tol: float) -> bool:
    """Two crossing fingerprints name the same crossing: same keys, cdists within ``tol``."""
    return a[:2] == b[:2] and abs(a[2] - b[2]) <= tol and abs(a[3] - b[3]) <= tol


def check_recompute_preserves_ids(case: "Case") -> int:
    """MUTATES (fresh build only): ``compute_intersections(preserve_ids=True)``
    keeps the id of every crossing that reappears.

    Crossings registered on blast children (image bridges, not manifolds) are
    not recomputed and may drop out; every crossing that IS recomputed keeps
    its id, keys and cdists.
    """
    tol = _tol(case)
    before = _crossing_fingerprint(case)
    case.workbench.compute_intersections(case.fixed_points, preserve_ids=True)
    after = _crossing_fingerprint(case)
    count = 0
    for intersection_id, value in after.items():
        old = before.get(intersection_id)
        if old is not None and _same_crossing(old, value, tol):
            count += 1
            continue
        moved = [iid for iid, prior in before.items() if _same_crossing(prior, value, tol)]
        assert not moved, f"crossing {moved} reappeared under the new id {intersection_id}"
    return count


def check_recompute_is_idempotent(case: "Case") -> int:
    """MUTATES (fresh build only): recomputing the crossings twice gives the same
    crossings under the same ids (cdists to the registry's ``cdist_tol``)."""
    tol = _tol(case)
    case.workbench.compute_intersections(case.fixed_points, preserve_ids=True)
    first = _crossing_fingerprint(case)
    case.workbench.compute_intersections(case.fixed_points, preserve_ids=True)
    second = _crossing_fingerprint(case)
    assert set(second) == set(first)
    for intersection_id, value in first.items():
        assert _same_crossing(value, second[intersection_id], tol), intersection_id
    return len(first)


# --------------------------------------------------------------------------- #
# Anchors
# --------------------------------------------------------------------------- #
def check_one_anchor_per_unstable_branch(case: "Case") -> int:
    """Every unstable branch carries exactly one anchor (author rule, 2026-10-05)."""
    per_branch: dict = defaultdict(list)
    for intersection_id in _anchor_ids(case):
        per_branch[case.registry[intersection_id].manifold_a_key].append(intersection_id)
    count = 0
    for key in case.workbench.manifolds:
        if key[1] != "unstable":
            continue
        ids = per_branch.get(key, [])
        assert len(ids) == 1, f"{len(ids)} anchors on unstable branch {key[1:]}: {ids}"
        count += 1
    return count


def check_anchor_sign_matches_eigendirections(case: "Case") -> int:
    """An anchor's ``crossing_sign`` is the handedness of the two oriented first segments."""
    count = 0
    for intersection_id in _anchor_ids(case):
        anchor = case.registry[intersection_id]
        u_nodes = case.workbench.manifolds[anchor.manifold_a_key].get_point_array()
        s_nodes = case.workbench.manifolds[anchor.manifold_b_key].get_point_array()
        u_dir = np.asarray(u_nodes[1], dtype=float) - np.asarray(u_nodes[0], dtype=float)
        s_dir = np.asarray(s_nodes[1], dtype=float) - np.asarray(s_nodes[0], dtype=float)
        expected = int(np.sign(u_dir[0] * s_dir[1] - u_dir[1] * s_dir[0]))
        assert anchor.crossing_sign == expected, intersection_id
        count += 1
    return count


# --------------------------------------------------------------------------- #
# Bridges
# --------------------------------------------------------------------------- #
def check_bridge_partial_iff_no_id(case: "Case") -> int:
    """A bridge is partial exactly when it has no ``BridgeId``."""
    bridges = case.workbench.bridges
    for bridge in bridges:
        assert bridge.partial == (bridge.id is None), bridge
        assert bridge.manifold_key is not None
    return len(bridges)


def check_bridge_id_in_unstable_order(case: "Case") -> int:
    """``BridgeId = (first, second)`` is its endpoints in strictly increasing unstable cdist.

    The one exempt tie is a zero-length bridge between two anchors, which
    only the known two-anchors-per-branch issue creates (see the anchor law).
    """
    registry = case.registry
    anchors = set(_anchor_ids(case))
    count = 0
    for bridge in _identified_bridges(case):
        assert bridge.id == (bridge.first_intersection, bridge.second_intersection)
        first, second = (registry[i] for i in bridge.id)
        if set(bridge.id) <= anchors:
            continue
        assert float(first.unstable_cdist) < float(second.unstable_cdist), bridge.id
        count += 1
    return count


def check_bridge_endpoints_on_own_branch(case: "Case") -> int:
    """A bridge's two crossings lie on ITS unstable branch, inside its cdist span.

    The codebase-audit regression (2026-09): a registry-wide nearest-cdist
    lookup handed the first bridge of one period-3 branch another branch's
    anchor.
    """
    registry = case.registry
    tol = _tol(case)
    count = 0
    for bridge in _identified_bridges(case):
        first, second = (registry[i] for i in bridge.id)
        root_u = float(bridge.root.get_cdist("unstable"))
        tail_u = float(bridge.tail.get_cdist("unstable"))
        lo, hi = sorted((root_u, tail_u))
        for end in (first, second):
            assert end.manifold_a_key == bridge.manifold_key, (
                f"bridge {bridge.id} on {bridge.manifold_key[1:]} got an endpoint "
                f"from {end.manifold_a_key[1:]}"
            )
            assert _within(float(end.unstable_cdist), lo, hi, tol), bridge.id
        count += 2
    return count


def check_bridges_at_exact(case: "Case") -> int:
    """``bridges_at(x)`` lists exactly the identified bridges ending at ``x``, once each."""
    workbench = case.workbench
    expected: dict[int, set] = defaultdict(set)
    for bridge in _identified_bridges(case):
        for end in bridge.id:
            expected[end].add(bridge.id)
    count = 0
    for intersection_id, _ix in case.registry:
        at = list(workbench.bridges_at(intersection_id))
        assert len(at) == len(set(at)), f"duplicate ids in bridges_at({intersection_id})"
        assert set(at) == expected.get(intersection_id, set()), intersection_id
        count += 1
    return count


def check_bridge_single_copy(case: "Case") -> int:
    """Each ``BridgeId`` is registered once and ``bridge(id)`` returns that copy."""
    workbench = case.workbench
    identified = _identified_bridges(case)
    ids = [bridge.id for bridge in identified]
    assert len(ids) == len(set(ids)), "a BridgeId is registered twice"
    for bridge in identified:
        assert workbench.bridge(bridge.id) is bridge, bridge.id
    return len(identified)


def check_bridges_do_not_overlap(case: "Case") -> int:
    """Two bridges of one unstable branch only touch: their cdist spans never overlap."""
    registry = case.registry
    tol = _tol(case)
    spans: dict = defaultdict(list)
    for bridge in _identified_bridges(case):
        spans[bridge.manifold_key].append(
            tuple(float(registry[i].unstable_cdist) for i in bridge.id)
        )
    count = 0
    for key, branch_spans in spans.items():
        branch_spans.sort()
        for (_lo_a, hi_a), (lo_b, hi_b) in zip(branch_spans, branch_spans[1:]):
            assert lo_b >= hi_a - tol, f"overlapping bridges on {key[1:]}: {hi_a} > {lo_b}"
            count += 1
    return count


def check_bridge_image_round_trip(case: "Case") -> int:
    """A bridge's derived image bridges lie on the advanced branch, and wherever the
    preimage of an image bridge is registered it contains the bridge."""
    workbench = case.workbench
    count = 0
    for bridge in _identified_bridges(case):
        image = workbench.image_bridges(bridge.id)
        if not image:
            continue
        image_key = bridge.fixed_point.advance_key(bridge.manifold_key, 1)
        for child_id in image:
            assert workbench.bridge(child_id).manifold_key == image_key, (bridge.id, child_id)
            preimage = workbench.preimage_bridges(child_id)
            if preimage is None:
                continue
            assert bridge.id in preimage, f"{child_id} is an image of {bridge.id} but not back"
            count += 1
    return count


# --------------------------------------------------------------------------- #
# Partition and holes
# --------------------------------------------------------------------------- #
def check_partition_covers_branch(case: "Case") -> int:
    """Each side's partition runs anchor to branch end in contiguous, positional elements,
    and every open element is a punched hole region of that side."""
    count = 0
    for trellis, result in _partitions(case):
        tol = float(trellis.registry.cdist_tol)
        intervals = result.intervals
        branch = trellis.branch(result.branch_key)
        s_max = float(trellis.intersection(branch.ordered_ids()[-1]).stable_cdist)
        assert intervals[0].lo_cdist <= tol
        assert intervals[-1].hi_cdist == pytest.approx(s_max, abs=tol)
        assert [iv.element_id for iv in intervals] == list(range(len(intervals)))
        for interval in intervals:
            assert interval.branch_key == result.branch_key and interval.side == result.side
            assert result.element(interval.element_id) is interval
        for prev, nxt in zip(intervals, intervals[1:]):
            assert prev.hi_cdist == pytest.approx(nxt.lo_cdist, abs=tol)
        open_spans = {
            frozenset((iv.lo_id, iv.hi_id))
            for iv in intervals
            if not iv.closed_lo and not iv.closed_hi and iv.lo_id != iv.hi_id
        }
        hole_spans = {
            frozenset(hole.bounding_ids)
            for hole in trellis.holes
            if hole.bounding_ids is not None
            and any(row == result.side for _iid, _which, row in (hole.openings or []))
        }
        assert open_spans <= hole_spans, open_spans - hole_spans
        count += 1
    return count


def _owners(intervals, cdist: float, tol: float) -> list[int]:
    """Element ids covering a cdist, from the raw interval fields (an end only when closed)."""
    owners = []
    for interval in intervals:
        if cdist < interval.lo_cdist - tol or cdist > interval.hi_cdist + tol:
            continue
        if abs(cdist - interval.lo_cdist) <= tol:
            covered = interval.closed_lo
        elif abs(cdist - interval.hi_cdist) <= tol:
            covered = interval.closed_hi
        else:
            covered = True
        if covered:
            owners.append(interval.element_id)
    return owners


def check_partition_unique_owner(case: "Case") -> int:
    """Every crossing of a partitioned branch has exactly one owning element per side."""
    count = 0
    for trellis, result in _partitions(case):
        tol = float(trellis.registry.cdist_tol)
        branch = trellis.branch(result.branch_key)
        for intersection_id in branch.intersection_ids:
            cdist = float(trellis.intersection(intersection_id).stable_cdist)
            owners = _owners(result.intervals, cdist, tol)
            assert owners == [result.element_of_intersection[intersection_id]], (
                f"crossing {intersection_id} on {result.branch_key[1:]} ({result.side}) "
                f"covered by {owners}"
            )
            count += 1
        assert set(result.element_of_intersection) == set(branch.intersection_ids)
    return count


def check_partition_singletons(case: "Case") -> int:
    """A pinched singleton ``[x, x]`` is closed at both ends and owns exactly ``x``."""
    count = 0
    for _trellis, result in _partitions(case):
        for interval in result.intervals:
            if interval.lo_id is None or interval.lo_id != interval.hi_id:
                continue
            assert interval.closed_lo and interval.closed_hi
            assert interval.lo_cdist == interval.hi_cdist
            owned = [
                iid
                for iid, element_id in result.element_of_intersection.items()
                if element_id == interval.element_id
            ]
            assert owned == [interval.lo_id]
            count += 1
    return count


def check_holes_are_classified(case: "Case") -> int:
    """Every hole has a bridge side, its anchorward bound first, finite coordinates and an orbit."""
    count = 0
    for trellis, hole in _holes(case):
        assert hole.bridge_side in ("left", "right")
        near, far = hole.bounding_ids
        assert hole.near_intersection_id == near
        assert float(trellis.intersection(near).stable_cdist) <= float(
            trellis.intersection(far).stable_cdist
        ) or trellis.intersection(near).manifold_b_key != trellis.intersection(far).manifold_b_key
        assert np.isfinite(np.asarray(hole.coords, dtype=float)).all()
        assert hole.iterate is not None and hole.origin is not None
        count += 1
    return count


def check_backward_holes_keep_bridge_side(case: "Case") -> int:
    """I1 (author rule, 2026-10-05): a hole maps backward onto the same side of the
    BRIDGE it is punched in -- never "the same side of the stable manifold".

    Every hole of an origin carries the ``bridge_side`` of the origin's
    reference hole (iterate 0), parity-flipped under an orientation-reversing
    map; and, measured on the curves, the reference hole's coordinates carried
    back ``|iterate|`` steps by the real inverse map sit on that side of the
    propagated hole's own bridge. Nothing is asserted about the stable-manifold
    ROW of a backward hole: it may differ from its origin's, and that is
    correct.
    """
    inverse = case.workbench.dynamical_system.map_inv
    count = 0
    for trellis in _per_fp_trellises(case):
        preserving = trellis.orientation_preserving
        reference = {
            hole.origin: hole
            for hole in trellis.holes
            if hole.pair is not None and hole.iterate == 0
        }
        for hole in trellis.holes:
            if hole.iterate == 0 and hole.pair is not None:
                continue
            origin = reference.get(hole.origin)
            assert origin is not None, f"hole {hole.bounding_ids} has no reference hole"
            flip = not preserving and hole.iterate % 2 == 1
            expected = origin.bridge_side
            if flip:
                expected = "right" if expected == "left" else "left"
            assert hole.bridge_side == expected, (
                f"hole of origin {hole.origin} at iterate {hole.iterate}: bridge side "
                f"{hole.bridge_side}, origin {origin.bridge_side}"
            )
            if hole.pair is None:
                point = np.asarray(origin.coords, dtype=float)
                for _ in range(-hole.iterate):
                    point = np.asarray(inverse(point), dtype=float)
                bridge = trellis.bridge_between(*hole.bounding_ids)
                assert bridge is not None and bridge.id is not None, hole.bounding_ids
                measured = _point_side_of_bridge(case, bridge, point)
                assert measured == expected, (
                    f"origin {hole.origin}'s hole carried back {-hole.iterate} step(s) "
                    f"lies on the {measured} of bridge {bridge.id}, expected {expected}"
                )
            count += 1
    return count


def check_direct_hole_side_is_coordinate_side(case: "Case") -> int:
    """A direct hole's ``bridge_side`` is the side of its bridge its coordinates sit on.

    The combinatorial crossing-sign rule (2026-10-02) must reproduce the
    geometric measurement hole for hole.
    """
    count = 0
    for trellis, hole in _holes(case):
        if hole.pair is None:
            continue
        bridge = bridge_for_pair(trellis, hole.pair)
        assert bridge is not None, hole.pair
        measured = _point_side_of_bridge(case, bridge, np.asarray(hole.coords, dtype=float))
        assert hole.bridge_side == measured, (
            f"direct hole {hole.bounding_ids} (iterate {hole.iterate}): side "
            f"{hole.bridge_side}, coordinates on the {measured}"
        )
        count += 1
    return count


def check_direct_hole_opens_inward_pair(case: "Case") -> int:
    """Openings law (1): a direct hole opens its inward pair -- outward of the near
    bound, anchorward of the far bound -- on its own bridge's row."""
    count = 0
    for trellis, hole in _holes(case):
        if hole.pair is None:
            continue
        near, far = hole.bounding_ids
        bridge = bridge_for_pair(trellis, hole.pair)
        assert hole.openings, hole.bounding_ids
        for intersection_id, which, row in hole.openings:
            assert which == ("outward" if intersection_id == near else "anchorward")
            endpoint = "first" if intersection_id == bridge.id[0] else "second"
            assert row == row_of_end(trellis, bridge.id, endpoint), (hole.bounding_ids, row)
        count += 1
    return count


def check_openings_on_own_bridge_row(case: "Case") -> int:
    """Openings law (1'): every opening of every hole sits on the row ITS OWN bridge has
    at that bound. Rows are never compared along an orbit: a backward hole keeps
    its side of the bridge (I1), not its row of the stable manifold."""
    count = 0
    for trellis, hole in _holes(case):
        bridge = trellis.bridge_between(*hole.bounding_ids)
        assert bridge is not None and bridge.id is not None, hole.bounding_ids
        for intersection_id, _which, row in hole.openings or []:
            endpoint = "first" if intersection_id == bridge.id[0] else "second"
            assert row == row_of_end(trellis, bridge.id, endpoint), (hole.bounding_ids, row)
            count += 1
    return count


def check_openings_missing_only_at_anchor_or_tail(case: "Case") -> int:
    """Openings law (3): a hole bound opens nothing only at an anchor or a trimmed tail."""
    count = 0
    for trellis, hole in _holes(case):
        opened = {iid for iid, _w, _r in hole.openings or []}
        tol = float(trellis.registry.cdist_tol)
        for intersection_id in hole.bounding_ids:
            if intersection_id in opened:
                continue
            ix = trellis.intersection(intersection_id)
            branch = trellis.branch(ix.manifold_b_key)
            tail = branch.ordered_ids()[-1] if branch is not None and len(branch) else None
            assert float(ix.stable_cdist) <= tol or intersection_id == tail, (
                f"hole {hole.bounding_ids} (iterate {hole.iterate}) opens nothing at "
                f"{intersection_id}, which is neither an anchor nor the branch end"
            )
        count += 1
    return count


def check_propagated_holes_land_on_predicted_branch(case: "Case") -> int:
    """A hole at iterate ``n`` lands in a bridge on ``cycle[(pos_ref + n) % k]``, the
    unstable branch cycle position of its reference bridge advanced by ``n``."""
    count = 0
    for trellis, fp in zip(_per_fp_trellises(case), case.fixed_points):
        cycle = fp.branch_cycle("unstable")
        k = len(cycle)
        reference_position = {}
        for pair in trellis.pseudoneighbors:
            if not pair.is_reference:
                continue
            bridge = bridge_for_pair(trellis, pair)
            assert bridge is not None and bridge.manifold_key in cycle
            reference_position[pair.as_tuple()] = cycle.index(bridge.manifold_key)
        for hole in trellis.holes:
            if hole.pair is not None:
                continue
            expected = cycle[(reference_position[hole.origin] + hole.iterate) % k]
            bridge = trellis.bridge_between(*hole.bounding_ids)
            assert bridge.manifold_key == expected, (
                f"hole of origin {hole.origin} at iterate {hole.iterate} landed on "
                f"{bridge.manifold_key[1:]}, expected {expected[1:]}"
            )
            count += 1
    return count


def check_propagation_terminates(case: "Case") -> int:
    """No origin punches the same bounding region twice (propagation stops at periodicity)."""
    count = 0
    for trellis in _per_fp_trellises(case):
        propagated = [hole for hole in trellis.holes if hole.iterate not in (0, None)]
        for origin in {hole.origin for hole in propagated}:
            regions = [
                frozenset(hole.bounding_ids) for hole in propagated if hole.origin == origin
            ]
            assert len(regions) == len(set(regions)), origin
            count += 1
    return count


def check_no_direct_hole_beyond_fundamental(case: "Case") -> int:
    """The FIRM half of holes-backward-only: nothing at iterate ``>= k_value`` is punched.

    The ``+1 .. +(k_value - 1)`` exemption is provisional and not tested.
    """
    count = 0
    for trellis, fp in zip(_per_fp_trellises(case), case.fixed_points):
        k = fp.k_value
        for hole in trellis.holes:
            assert hole.iterate < k, (hole.bounding_ids, hole.iterate)
            count += 1
        for pair in trellis.pseudoneighbors:
            if pair.iterate is not None and pair.iterate >= k:
                assert pair.hole is None, pair.as_tuple()
                count += 1
    return count


def check_reference_pairs_valid(case: "Case") -> int:
    """A reference pair is consecutive on its stable AND its unstable branch, on the
    strong pip's own branch, inside the reference window ``[f^k(q0), q0]``."""
    count = 0
    for trellis, fp in zip(_per_fp_trellises(case), case.fixed_points):
        pip = trellis.strong_pip
        tol = float(trellis.registry.cdist_tol)
        pip_ix = trellis.intersection(pip)
        top = float(pip_ix.stable_cdist)
        _key, bottom, _from_table = trellis.image_cdist(pip, fp.k_value, "stable")
        for pair in trellis.pseudoneighbors:
            if not pair.is_reference:
                continue
            a, b = pair.as_tuple()
            assert pair.branch_key == pip_ix.manifold_b_key
            stable_ids = trellis.branch_containing(a, "stable").intersection_ids
            unstable_ids = trellis.branch_containing(a, "unstable").intersection_ids
            assert b in stable_ids and abs(stable_ids.index(a) - stable_ids.index(b)) == 1
            assert b in unstable_ids and abs(unstable_ids.index(a) - unstable_ids.index(b)) == 1
            for iid in (a, b):
                assert _within(float(trellis.intersection(iid).stable_cdist), float(bottom), top, tol)
            count += 1
    return count


# --------------------------------------------------------------------------- #
# Arrangement
# --------------------------------------------------------------------------- #
def check_arrangement_euler(case: "Case") -> int:
    """``V - E + F = 2`` per connected component, with one open face per component at least."""
    arrangement = case.session.arrangement()
    components = arrangement.component_count
    assert arrangement.euler_characteristic == 2 * components
    assert len(arrangement.open_faces) >= components
    return components


def check_arrangement_regions_sound(case: "Case") -> int:
    """Every region is a closed cycle of arcs between consecutive corners whose
    representative point survives an independent ray cast."""
    arrangement = case.session.arrangement()
    for region in arrangement.regions:
        assert region.is_closed
        assert len(region.corners) == len(region.arcs) >= 2
        for index, arc in enumerate(region.arcs):
            assert arc.tail_id == region.corners[index]
            assert arc.head_id == region.corners[(index + 1) % len(region.corners)]
        assert region.verify_representative_point(), region.corners
        assert region.contains(region.representative_point)
    return len(arrangement.regions)


def check_arrangement_regions_disjoint(case: "Case") -> int:
    """No region's representative point lies in any OTHER region; a containing face
    swallows more than one region."""
    arrangement = case.session.arrangement()
    count = 0
    for region in arrangement.regions:
        point = region.representative_point
        for other in arrangement.regions:
            if other is not region:
                assert not other.contains(point), (region.corners, other.corners)
        count += 1
    for face in arrangement.containing_faces:
        assert face.is_closed and not face.is_minimal and face not in arrangement.regions
        swallowed = [r for r in arrangement.regions if face.contains(r.representative_point)]
        assert len(swallowed) > 1
    return count


def check_arrangement_image_of(case: "Case") -> int:
    """A region's combinatorial image is where the map sends its representative point,
    with the source's area."""
    arrangement = case.session.arrangement()
    forward = case.workbench.dynamical_system.map
    count = 0
    for region in arrangement.regions:
        image = arrangement.image_of(region, 1)
        if image is None:
            continue
        mapped = np.asarray(forward(np.asarray(region.representative_point, dtype=float)))
        located = [other for other in arrangement.regions if other.contains(mapped)]
        assert located == [image], (region.corners, [r.corners for r in located])
        assert abs(image.area) == pytest.approx(abs(region.area), rel=REGION_AREA_RTOL)
        count += 1
    return count


def check_arrangement_preimage_inverts_image(case: "Case") -> int:
    """Where both directions are recorded, image and preimage undo each other."""
    arrangement = case.session.arrangement()
    count = 0
    for region in arrangement.regions:
        image = arrangement.image_of(region, 1)
        if image is None:
            continue
        back = arrangement.preimage_of(image, 1)
        if back is None:
            continue
        assert back is region, region.corners
        count += 1
    return count


# --------------------------------------------------------------------------- #
# Bridge classes
# --------------------------------------------------------------------------- #
def _element_ends(case: "Case", bridge_id) -> tuple:
    """The partition element at each end of a bridge, on that end's row."""
    trellis = _full(case)
    by_branch_side = {(r.branch_key, r.side): r for _t, r in _partitions(case)}
    refs = []
    for endpoint, intersection_id in zip(("first", "second"), bridge_id):
        row = row_of_end(trellis, bridge_id, endpoint)
        result = by_branch_side[(trellis.intersection(intersection_id).manifold_b_key, row)]
        refs.append(result.ref(result.element_of_intersection[intersection_id]))
    return tuple(refs)


def check_every_bridge_in_one_class(case: "Case") -> int:
    """Every identified bridge of the trellis lands in exactly one class."""
    trellis = _full(case)
    table = case.session.bridge_classes()
    expected = sorted(bridge.id for bridge in trellis.bridges if bridge.id is not None)
    placed = [bid for entry in table for bid in entry.bridge_ids]
    assert len(placed) == len(set(placed)), "a bridge landed in two classes"
    assert sorted(placed) == expected
    return len(placed)


def check_class_orientation_anchor_outward(case: "Case") -> int:
    """A class runs anchor outward (``anchor_outward_key``), and a member's direction is
    ``+1`` exactly when its first end is the class source; a loop sits in one element."""
    trellis = _full(case)
    fixed_points = trellis.fixed_points
    table = case.session.bridge_classes()
    count = 0
    for entry in table:
        cls = entry.bridge_class
        assert not cls.is_loop or not entry.active
        assert anchor_outward_key(cls.source, fixed_points) <= anchor_outward_key(
            cls.target, fixed_points
        )
        for member in entry.members:
            x, y = _element_ends(case, member.bridge_id)
            if member.is_loop:
                assert x == y == member.loop_element
            else:
                assert {x, y} == {cls.source, cls.target}
                assert member.direction == (+1 if x == cls.source else -1)
                assert oriented_class(x, y, fixed_points) == (cls, member.direction)
            count += 1
    return count


def check_row_of_end_is_geometry(case: "Case") -> int:
    """The combinatorial row (crossing sign alone) is the side the bridge geometrically
    approaches the crossing from, wherever the geometry has an answer."""
    trellis = _full(case)
    count = 0
    mismatches = []
    for bridge in trellis.bridges:
        if bridge.id is None:
            continue
        for endpoint, intersection_id in zip(("first", "second"), bridge.id):
            if trellis.intersection(intersection_id).seg_ids is None:
                continue  # an anchor: no stable node below it
            measured = _geometric_row(case, bridge, intersection_id)
            if measured is None:
                continue
            combinatorial = row_of_end(trellis, bridge.id, endpoint)
            if combinatorial != measured:
                mismatches.append((bridge.id, endpoint, combinatorial, measured))
            count += 1
    assert not mismatches, mismatches
    return count


def check_bridge_rows_consistent(case: "Case") -> int:
    """I2: a bridge whose two crossings share a stable branch has one row at both ends."""
    trellis = _full(case)
    count = 0
    for bridge in trellis.bridges:
        if bridge.id is None:
            continue
        assert bridge_row_violation(trellis, bridge) is None, bridge.id
        keys = [trellis.intersection(i).manifold_b_key for i in bridge.id]
        if keys[0] != keys[1]:
            continue
        first, second = rows_of_bridge(trellis, bridge.id)
        assert first == second, (bridge.id, first, second)
        count += 1
    return count


def check_anchor_bridge_class_leads_its_tangle(case: "Case") -> int:
    """A tangle's anchor bridge (leaving an anchor) is a ``+1`` member of an ACTIVE class,
    and that class comes first in its tangle's group (smallest member cdist: 0)."""
    trellis = _full(case)
    table = case.session.bridge_classes()
    entries = list(table)
    anchors = set(_anchor_ids(case))
    tol = _tol(case)
    count = 0
    for bridge in trellis.bridges:
        if bridge.id is None or bridge.id[0] not in anchors or bridge.id[1] in anchors:
            continue
        entry = table.entry_of(bridge.id)
        assert entry.active, bridge.id
        assert entry.member(bridge.id).direction == +1, bridge.id
        assert entry.min_unstable_cdist == pytest.approx(0.0, abs=tol)
        group = [other for other in entries if other.tangle == entry.tangle]
        assert group[0].min_unstable_cdist == pytest.approx(0.0, abs=tol)
        count += 1
    return count


def check_class_table_order(case: "Case") -> int:
    """Classes group by tangle (connecting classes last), then by smallest member
    first-crossing unstable cdist, then ``class_sort_key``; members by ``BridgeId``."""
    trellis = _full(case)
    fixed_points = trellis.fixed_points
    table = case.session.bridge_classes()
    entries = list(table)
    for entry in entries:
        assert entry.bridge_ids == sorted(entry.bridge_ids)
        assert entry.min_unstable_cdist == min(
            float(trellis.intersection(bid[0]).unstable_cdist) for bid in entry.bridge_ids
        )
    keys = [
        (
            len(fixed_points) if entry.tangle is None else entry.tangle,
            entry.min_unstable_cdist,
            class_sort_key(entry.bridge_class, fixed_points),
        )
        for entry in entries
    ]
    assert keys == sorted(keys)
    return len(entries)


def check_classes_do_not_mix_tangles(case: "Case") -> int:
    """A class names one tangle (``entry.tangle``) exactly when its elements and its
    members' unstable branches all belong to that fixed point; every tangle has classes."""
    trellis = _full(case)
    fixed_points = trellis.fixed_points
    table = case.session.bridge_classes()
    seen = set()
    for entry in table:
        cls = entry.bridge_class
        owners = {id(cls.source.fixed_point), id(cls.target.fixed_point)}
        for bridge_id in entry.bridge_ids:
            owners.add(id(trellis.intersection(bridge_id[0]).manifold_a_key[0]))
        if len(owners) == 1:
            (owner,) = owners
            assert entry.tangle is not None and id(fixed_points[entry.tangle]) == owner
            seen.add(owner)
        else:
            assert entry.tangle is None
    assert seen == {id(fp) for fp in case.fixed_points}, "every tangle contributes classes"
    return len(table)


def check_classes_use_every_orbit_branch(case: "Case") -> int:
    """The classes of a fixed point name elements on every stable branch of its cycle
    (all k orbit indices of a period-k orbit, both branches of an inversion point)."""
    table = case.session.bridge_classes()
    count = 0
    for fp in case.fixed_points:
        named = {
            ref.branch_key
            for cls in table.classes
            for ref in (cls.source, cls.target)
            if ref.fixed_point is fp
        }
        expected = set(fp.branch_cycle("stable"))
        assert named == expected, (sorted(k[2:] for k in named), sorted(k[2:] for k in expected))
        count += len(expected)
    return count


def check_only_active_classes_lettered(case: "Case") -> int:
    """The session letters exactly the active classes, each with its own letter."""
    table = case.session.bridge_classes()
    letters = []
    for entry in table:
        assert (entry.letter is not None) == entry.active, entry.bridge_class
        if entry.letter is not None:
            letters.append(entry.letter)
    assert len(letters) == len(set(letters))
    return len(table)


# --------------------------------------------------------------------------- #
# Iterated partition, element names and the minimal trellis
# --------------------------------------------------------------------------- #
def _flanks(result, boundary_id) -> tuple:
    """(closed_hi before the boundary, closed_lo after it, whether it is a singleton)."""
    before = [iv for iv in result.intervals if iv.hi_id == boundary_id and iv.lo_id != boundary_id]
    after = [iv for iv in result.intervals if iv.lo_id == boundary_id and iv.hi_id != boundary_id]
    singleton = any(iv.lo_id == iv.hi_id == boundary_id for iv in result.intervals)
    return (
        before[0].closed_hi if before else None,
        after[0].closed_lo if after else None,
        singleton,
    )


def _boundaries(result) -> set:
    """Every element boundary id of one partition result."""
    return {iv.lo_id for iv in result.intervals} | {iv.hi_id for iv in result.intervals}


def check_iterated_child_inside_parent(case: "Case") -> int:
    """Each iterated element lies inside the homotopy parent that owns its midpoint."""
    iterated = case.session.iterated_partition()
    homotopy = iterated.homotopy
    tol = homotopy.tol
    count = 0
    for key, result in iterated.results.items():
        parent_result = homotopy.result(*key)
        for interval in result.intervals:
            parent = parent_result.element(interval.parent_element_id)
            assert interval.lo_cdist >= parent.lo_cdist - tol
            assert interval.hi_cdist <= parent.hi_cdist + tol
            mid = 0.5 * (interval.lo_cdist + interval.hi_cdist)
            assert owns_cdist(parent, mid, tol)
            assert interval.ref.branch_key == key[0] and interval.ref.side == key[1]
            count += 1
        assert len(result.intervals) >= len(parent_result.intervals)
    return count


def check_iterated_keeps_homotopy_boundaries(case: "Case") -> int:
    """Every homotopy boundary survives the cut with the same closedness on both flanks,
    and every NEW boundary is an endpoint of an image bridge cutting that row."""
    iterated = case.session.iterated_partition()
    homotopy = iterated.homotopy
    minimal = case.session.minimal_trellis()
    full = minimal.trellis
    image_ends: dict = defaultdict(set)
    for bid in minimal.image_bridge_ids:
        for endpoint, which in zip(bid, ("first", "second")):
            key = (full.intersection(endpoint).manifold_b_key, row_of_end(full, bid, which))
            image_ends[key].add(endpoint)
    count = 0
    for key, result in iterated.results.items():
        parent_result = homotopy.result(*key)
        for boundary in _boundaries(parent_result) - {None}:
            assert _flanks(result, boundary) == _flanks(parent_result, boundary), (key, boundary)
            count += 1
        new = _boundaries(result) - _boundaries(parent_result)
        assert new <= image_ends.get(key, set()), (key, new)
    assert set(iterated.results) == set(homotopy.results)
    return count


def check_iterated_unique_owner(case: "Case") -> int:
    """Every crossing of a cut branch has exactly one iterated owner per side, and every
    stable arc midpoint of the minimal arrangement one owner per side."""
    iterated = case.session.iterated_partition()
    minimal = case.session.minimal_trellis()
    full = minimal.trellis
    tol = iterated.homotopy.tol
    count = 0
    for key, result in iterated.results.items():
        for iid in full.branch(key[0]).intersection_ids:
            cdist = float(full.intersection(iid).stable_cdist)
            owners = [iv for iv in result.intervals if owns_cdist(iv, cdist, tol)]
            assert len(owners) == 1, (key, iid, len(owners))
            assert result.element_of_intersection[iid] == owners[0].element_id
            count += 1
    for face in minimal.arrangement.faces:
        for arc in face.arcs:
            if arc.kind != "stable":
                continue
            mid = 0.5 * (
                float(full.intersection(arc.lo_id).stable_cdist)
                + float(full.intersection(arc.hi_id).stable_cdist)
            )
            for side in ("left", "right"):
                iterated.interval_at(arc.branch_key, side, mid)
    return count


def check_iterated_cut_provenance(case: "Case") -> int:
    """Every element under a same-branch, same-row image bridge records an image bridge
    covering it (``cut_by``)."""
    iterated = case.session.iterated_partition()
    minimal = case.session.minimal_trellis()
    full = minimal.trellis
    tol = iterated.homotopy.tol
    count = 0
    for bid in minimal.image_bridge_ids:
        a, b = bid
        branch = full.intersection(a).manifold_b_key
        if branch != full.intersection(b).manifold_b_key:
            continue
        row = row_of_end(full, bid, "first")
        if row != row_of_end(full, bid, "second"):
            continue
        result = iterated.result(branch, row)
        near_c, far_c = sorted(float(full.intersection(i).stable_cdist) for i in bid)
        under = [
            iv for iv in result.intervals
            if iv.lo_cdist >= near_c - tol and iv.hi_cdist <= far_c + tol
        ]
        assert under, bid
        for interval in under:
            assert interval.cut_by is not None and minimal.is_image_bridge(interval.cut_by)
            lo_c, hi_c = sorted(float(full.intersection(i).stable_cdist) for i in interval.cut_by)
            assert lo_c <= interval.lo_cdist + tol and hi_c >= interval.hi_cdist - tol
            count += 1
    return count


def check_names_agree_with_structure(case: "Case") -> int:
    """Each iterated element's name is ``<side letter>_<parent id + 1>``, round-trips
    through ``ref_of``/``lookup``, and its children lie within their parent."""
    dual = case.session.dual_graph()
    naming = dual.naming
    iterated = case.session.iterated_partition()
    homotopy = iterated.homotopy
    count = 0
    for result in iterated.as_list():
        for interval in result.intervals:
            ref = result.ref(interval.element_id)
            name = naming.name(ref)
            parent = naming.parent_of(ref)
            assert parent.element_id == interval.parent_element_id
            assert parent.branch_key == ref.branch_key and parent.side == ref.side
            assert name.subscript == interval.parent_element_id + 1
            assert name.side_letter == result.side[0].upper()
            assert naming.ref_of(name) == ref
            assert naming.lookup(name.text) == ref
            assert naming.homotopy_name(parent).subscript == name.subscript
            assert ref in naming.children_of(parent)
            count += 1
    for parent in naming.homotopy_refs:
        span = homotopy.element(parent)
        for child in naming.children_of(parent):
            child_interval = iterated.element(child)
            assert child_interval.lo_cdist >= span.lo_cdist - homotopy.tol
            assert child_interval.hi_cdist <= span.hi_cdist + homotopy.tol
    return count


def check_minimal_trellis_bridges(case: "Case") -> int:
    """The minimal trellis keeps every hole bridge, maps exactly the active ones along
    their registered image chains (``image_chain``, skipped pairs excepted), and
    keeps nothing else."""
    session = case.session
    minimal = session.minimal_trellis()
    full = minimal.trellis
    table = session.bridge_classes()
    hole_set = {frozenset(b) for b in minimal.hole_bridge_ids}
    for _trellis, hole in _holes(case):
        bridge = full.bridge_between(*hole.bounding_ids)
        if bridge is not None and bridge.id is not None:
            assert frozenset(bridge.id) in hole_set, hole.bounding_ids
    assert len(hole_set) == len(minimal.hole_bridge_ids)
    for bid in minimal.inert_hole_bridge_ids + minimal.unclassed_hole_bridge_ids:
        assert minimal.is_hole_bridge(bid) and bid not in minimal.image_chains
    for bid in minimal.inert_hole_bridge_ids:
        assert not table.entry_of(bid).active
    for bid, chain in minimal.image_chains.items():
        assert table.entry_of(bid).active
        expected = image_chain(full, bid)
        assert expected is not None, bid
        skipped = {frozenset(p) for parent, p in minimal.skipped_pairs if parent == bid}
        assert [frozenset(p) for p in chain] == [
            frozenset(p) for p in expected if frozenset(p) not in skipped
        ], bid
    for bid in minimal.unmapped:
        assert image_chain(full, bid) is None, bid
    image_set = {frozenset(b) for b in minimal.image_bridge_ids}
    assert not (image_set & hole_set)
    for bid in minimal.image_bridge_ids:
        parents = minimal.parents_of(bid)
        assert parents and all(minimal.is_hole_bridge(p) for p in parents)
        assert minimal.is_image_bridge(bid) and minimal.is_kept(bid)
    for _parent, pair in minimal.skipped_pairs:
        assert full.bridge_between(*pair) is None
    kept_ids = {frozenset(b.id) for b in minimal.kept}
    assert kept_ids == hole_set | image_set
    assert all(not b.partial for b in minimal.kept)
    all_ids = {frozenset(b.id) for b in full.bridges if b.id is not None}
    assert {frozenset(b) for b in minimal.dropped_bridge_ids} == all_ids - kept_ids
    return len(kept_ids)


def check_minimal_trellis_nodes_and_arrangement(case: "Case") -> int:
    """The kept nodes of a stable branch are the kept bridges' endpoints, the anchors and
    the branch end; the sparse arrangement carries exactly the kept bridges and closes up."""
    minimal = case.session.minimal_trellis()
    full = minimal.trellis
    tol = float(full.registry.cdist_tol)
    endpoints = {e for b in minimal.kept for e in b.id}
    count = 0
    for branch in full.stable_branches:
        kept = minimal.kept_ids(branch.key)
        expected = {i for i in endpoints if full.intersection(i).manifold_b_key == branch.key}
        expected |= {
            i for i in branch.intersection_ids if float(full.intersection(i).stable_cdist) <= tol
        }
        if len(branch):
            expected.add(branch.intersection_ids[-1])
        assert set(kept) == expected, branch.key[1:]
        cdists = [float(full.intersection(i).stable_cdist) for i in kept]
        assert cdists == sorted(cdists)
        count += len(kept)
    assert minimal.sparse.registry is full.registry
    arrangement = minimal.arrangement
    assert arrangement.sparse
    unstable_pairs = {
        frozenset((arc.lo_id, arc.hi_id))
        for face in arrangement.faces
        for arc in face.arcs
        if arc.kind == "unstable"
    }
    assert unstable_pairs == {frozenset(b.id) for b in minimal.kept}
    assert arrangement.euler_characteristic == 2 * arrangement.component_count
    for face in arrangement.faces:
        for arc in face.arcs:
            if arc.kind == "stable":
                kept = minimal.kept_ids(arc.branch_key)
                assert arc.lo_id in kept and arc.hi_id in kept
    return count


# --------------------------------------------------------------------------- #
# Dual graph
# --------------------------------------------------------------------------- #
def _stable_edge_keys(dual) -> set:
    """Every stable edge of the dual graph's arrangement."""
    return {
        arc.edge_key
        for face in dual.arrangement.faces
        for arc in face.arcs
        if arc.kind == "stable"
    }


def check_dual_node_structure(case: "Case") -> int:
    """Each stable edge has two open WALL side nodes, or one solid UNIFIED node on the
    fundamental segment; a wall cannot be crossed; faces and nodes point at each
    other consistently, and a face node carries its member regions' corners and bridges."""
    dual = case.session.dual_graph()
    count = 0
    for edge_key in _stable_edge_keys(dual):
        by_side = dual.nodes_of_edge(edge_key)
        left, right = by_side["left"], by_side["right"]
        if left is right:
            assert left.is_unified and left.traversable and left.node_side == "both"
            assert left.key == (edge_key, "both") and set(left.faces) == {"left", "right"}
            assert left.fundamental_span is not None
            face = left.face_on("left")
            assert left.other_face(face) is left.face_on("right")
        else:
            for node, side in ((left, "left"), (right, "right")):
                other = "right" if side == "left" else "left"
                assert node.sides == (side,) and not node.traversable
                assert node.key == (edge_key, side) and list(node.faces) == [side]
                assert node.fundamental_span is None
                assert node.element_on(side) == node.elements[side]
                # A wall cannot be crossed and faces one side only.
                with pytest.raises(ValueError):
                    node.other_face(node.faces[side])
                with pytest.raises(ValueError):
                    node.element_on(other)
        count += 1
    for node in dual.stable_nodes.values():
        for side in node.sides:
            face = node.face_on(side)
            assert (node, side) in face.stable_nodes
            assert face.side_of(node) == side or node.is_unified
        assert dual.stable_nodes[node.key] is node
    for face in dual.face_nodes:
        for node, side in face.stable_nodes:
            assert node.faces[side] is face
        for node in face.exits:
            assert node.is_unified
        assert set(face.corners) == {
            corner for region in face.faces for corner in region.corners if corner >= 0
        }
        assert set(face.bridge_ids) == {bid for region in face.faces for bid in region.bridge_ids}
    assert set(dual.unified_nodes) | set(dual.wall_nodes) == set(dual.stable_nodes.values())
    assert len(dual.stable_nodes) == 2 * len(_stable_edge_keys(dual)) - len(dual.unified_nodes)
    return count


def check_dual_face_nodes(case: "Case") -> int:
    """One face node per merged face: a single unbounded outer node, every region its own
    region node, and a merged node only where an inner component's outer face is glued;
    a face's lazy image is a face node of the graph, and only a region has one."""
    dual = case.session.dual_graph()
    assert dual.unbounded.is_unbounded and dual.unbounded.kind == "outer"
    assert sum(1 for face in dual.face_nodes if face.is_unbounded) == 1
    assert sum(len(face.faces) for face in dual.face_nodes) == len(dual.arrangement.faces)
    for index, face in enumerate(dual.face_nodes):
        assert face.index == index
        if len(face.faces) > 1:
            assert face.kind == "outer"
    for region in dual.arrangement.regions:
        assert dual.face_of(region).is_region
    for face in dual.face_nodes:
        image = dual.image_face(face)
        assert image is None or image in dual.face_nodes
        if not face.is_region:
            assert image is None
    return len(dual.face_nodes)


def _anchorward_tangent(node, trellis) -> Optional[np.ndarray]:
    """The unit direction toward the anchor at a stable edge's midpoint vertex."""
    poly = node.arc.polyline(trellis)
    if len(poly) < 2:
        return None
    middle = len(poly) // 2
    tangent = np.asarray(poly[max(middle - 1, 0)], dtype=float) - np.asarray(
        poly[min(middle + 1, len(poly) - 1)], dtype=float
    )
    norm = float(np.linalg.norm(tangent))
    return tangent / norm if norm > 0.0 else None


def check_dual_face_side_is_geometry(case: "Case") -> int:
    """A face attaches to the side node of the stable edge it geometrically lies on."""
    dual = case.session.dual_graph()
    trellis = dual.arrangement.trellis
    count = 0
    checkable = 0
    for region in dual.arrangement.regions:
        face = dual.face_of(region)
        inside = region.representative_point
        checkable += len(region.stable_arcs)
        for arc in region.stable_arcs:
            side = "right" if arc.reverse else "left"
            node = dual.nodes_of_edge(arc.edge_key)[side]
            look = _anchorward_tangent(node, trellis)
            if look is None:
                continue
            poly = node.arc.polyline(trellis)
            midpoint = poly[len(poly) // 2]
            assert _side_of(look, np.asarray(inside) - midpoint) == node.side_of(face), (
                region.corners, arc.lo_id, arc.hi_id,
            )
            assert node.faces[side] is face
            count += 1
    assert count == checkable, "every (region, stable edge) pair is checkable"
    return count


def check_dual_unified_is_pip_segment(case: "Case") -> int:
    """The unified edges are exactly the edges of each strong pip's OWN branch meeting
    ``(cdist(f^k(q0)), cdist(q0)]``; no other branch unifies."""
    dual = case.session.dual_graph()
    full = dual.trellis
    tol = float(full.registry.cdist_tol)
    expected = set()
    branches = set()
    for trellis, fp in zip(_per_fp_trellises(case), case.fixed_points):
        pip = trellis.strong_pip
        branch_key = full.intersection(pip).manifold_b_key
        branches.add(branch_key)
        image = full.iterate(pip, fp.k_value)
        if image is not None:
            bottom = float(full.intersection(image).stable_cdist)
        else:
            _key, bottom, _from_table = full.image_cdist(pip, fp.k_value, "stable")
        top = float(full.intersection(pip).stable_cdist)
        assert dual.fundamental_segments[branch_key] == (
            pytest.approx(bottom, abs=tol),
            pytest.approx(top, abs=tol),
        )
        for edge_key in _stable_edge_keys(dual):
            node = dual.nodes_of_edge(edge_key)["left"]
            if (
                node.branch_key == branch_key
                and node.hi_cdist > bottom + tol
                and node.lo_cdist < top - tol
            ):
                expected.add(edge_key)
    assert set(dual.fundamental_segments) == branches
    assert {node.edge_key for node in dual.unified_nodes} == expected
    for node in dual.unified_nodes:
        assert node.fundamental_span == dual.fundamental_segments[node.branch_key]
    return len(expected)


def check_dual_bipartite_degree(case: "Case") -> int:
    """The graph view is bipartite face/stable; a wall has degree 1, a unified node 2."""
    import networkx as nx

    dual = case.session.dual_graph()
    graph = dual.graph()
    faces = {n for n, d in graph.nodes(data=True) if d["bipartite"] == 0}
    stable = {n for n, d in graph.nodes(data=True) if d["bipartite"] == 1}
    assert nx.is_bipartite(graph)
    assert len(faces) == len(dual.face_nodes) and len(stable) == len(dual.stable_nodes)
    for a, b, data in graph.edges(data=True):
        assert {a[0], b[0]} == {"face", "stable"} and data["side"] in ("left", "right")
    for node in dual.stable_nodes.values():
        assert graph.nodes[("stable", node.key)]["traversable"] == node.traversable
        degree = graph.degree(("stable", node.key))
        assert degree == 1 or (node.is_unified and degree == 2)
    return len(dual.stable_nodes)


# --------------------------------------------------------------------------- #
# Symbolic dynamics
# --------------------------------------------------------------------------- #
def _resolved(dyn) -> list:
    """The class records that have an itinerary."""
    return [cd for cd in dyn.classes.values() if cd.itinerary is not None]


def check_itineraries_even(case: "Case") -> int:
    """Every itinerary has even length (it splits into disjoint consecutive pairs)."""
    dyn = case.session.symbolic_dynamics()
    count = 0
    for cd in _resolved(dyn):
        assert len(cd.itinerary) % 2 == 0, cd.letter
        count += 1
    return count


def check_itinerary_pairs_same_side(case: "Case") -> int:
    """Every itinerary pair is two elements of ONE side (a bridge connects same-side
    elements, possibly on different branches); no token is flagged cross-side."""
    dyn = case.session.symbolic_dynamics()
    count = 0
    for cd in _resolved(dyn):
        for x, y in zip(cd.itinerary[::2], cd.itinerary[1::2]):
            assert x.side == y.side, (cd.letter, x.label, y.label)
            count += 1
        assert not any(symbol.cross_side for symbol in cd.symbols), cd.letter
    return count


def check_itinerary_pairs_are_classes(case: "Case") -> int:
    """Each itinerary pair, read through its homotopy parents, is a loop (same parent) or
    one class or its inverse -- exactly the word's tokens, in order."""
    dyn = case.session.symbolic_dynamics()
    naming = dyn.naming
    fixed_points = dyn.fixed_points
    count = 0
    for cd in _resolved(dyn):
        expected = []
        for x, y in zip(cd.itinerary[::2], cd.itinerary[1::2]):
            px, py = naming.parent_of(x), naming.parent_of(y)
            if px == py:
                continue
            expected.append(oriented_class(px, py, fixed_points))
        assert [(s.bridge_class, s.direction) for s in cd.symbols] == expected, cd.letter
        count += len(expected) + 1
    return count


def check_landings_contained(case: "Case") -> int:
    """Each homotopy element of a class maps one step into ONE iterated element."""
    dyn = case.session.symbolic_dynamics()
    count = 0
    for cd in dyn.classes.values():
        for landing in cd.landings:
            assert landing is not None and landing.target is not None, cd.letter
            assert landing.contained, (cd.letter, landing)
            count += 1
    return count


def check_refined_children_inherit_word(case: "Case") -> int:
    """A split class's children keep the parent ``BridgeClass`` and each carries the
    parent's word, rewritten in refined symbols."""
    dyn = case.session.symbolic_dynamics()
    count = 0
    for cls, children in dyn.refined.items():
        cd = dyn.classes[cls]
        words = []
        for child in children:
            assert child.parent == cls
            assert child.name in dyn.refined_rules
            words.append([(s.bridge_class, s.direction) for s in dyn.refined_rules[child.name]])
        assert all(word == [(s.bridge_class, s.direction) for s in cd.symbols] for word in words)
        assert len(children) >= 2
        count += len(children)
    return count


def check_member_matching_consistent(case: "Case") -> int:
    """A matched member is listed under exactly the child it was matched to; every
    non-loop member of a split class is matched or listed unmatched."""
    dyn = case.session.symbolic_dynamics()
    count = 0
    for cls, children in dyn.refined.items():
        for member in dyn.classes[cls].entry.members:
            if member.is_loop:
                continue
            bid = member.bridge_id
            assert bid in dyn.member_refinement or bid in dyn.unmatched_members, bid
            child = dyn.member_refinement.get(bid)
            for other in children:
                assert (bid in other.members) == (other is child), bid
            count += 1
    return count


def check_inert_classes_outside_transitions(case: "Case") -> int:
    """No inert class is a node of the transition graph or a row of the matrix."""
    dyn = case.session.symbolic_dynamics()
    inert = {cls for cls, cd in dyn.classes.items() if cd.kind == "inert"}
    graph = dyn.transition_graph()
    for _name, data in graph.nodes(data=True):
        assert data["bridge_class"] not in inert
    names, _matrix = dyn.transition_matrix()
    assert set(names) == set(graph.nodes)
    return len(dyn.classes)


def check_matrix_is_token_counts(case: "Case") -> int:
    """Entry ``(i, j)`` counts the tokens of symbol ``j`` (inert tokens skipped) in the
    refined word of active symbol ``i``; virtual symbols are sink rows."""
    dyn = case.session.symbolic_dynamics()
    names, matrix = dyn.transition_matrix()
    index = {name: position for position, name in enumerate(names)}
    expected = np.zeros_like(matrix)
    nodes = {node.name: node for node in dyn.symbol_nodes()}
    for name, symbols in dyn.refined_rules.items():
        node = nodes.get(name)
        if node is None or node.kind != "active":
            continue
        for symbol in symbols:
            if symbol.kind == "inert" or symbol.base not in index:
                continue
            expected[index[name], index[symbol.base]] += 1
    assert matrix.tolist() == expected.tolist()
    for name, node in nodes.items():
        if node.kind == "virtual":
            assert not matrix[index[name]].any()
    return int(matrix.sum()) + len(names)


#: The read-only law checks by layer: ONE law-tier test per (layer, case) runs
#: every check of its layer (``helpers.law_tier.run_layer``).
LAYERS: dict[str, list[LawCheck]] = {
    "case_sanity": [
        check_case_builds,
        check_k_value_matches_inversion,
        check_every_branch_is_built,
        check_every_branch_crosses,
    ],
    "manifolds": [
        check_manifold_cdist_monotone,
        check_manifold_no_spikes,
        check_manifold_iterate_law,
        check_manifold_one_to_one,
    ],
    "crossings": [
        check_crossings_unstable_by_stable,
        check_no_same_stability_crossing,
        check_crossings_cdists_defined,
        check_crossing_cdist_bracketed,
        check_crossing_sign_is_cross_product,
        check_area_along_iterate_links,
        check_iterate_link_scales_by_beta,
    ],
    "anchors": [check_one_anchor_per_unstable_branch, check_anchor_sign_matches_eigendirections],
    "bridges": [
        check_bridge_partial_iff_no_id,
        check_bridge_id_in_unstable_order,
        check_bridge_endpoints_on_own_branch,
        check_bridges_at_exact,
        check_bridge_single_copy,
        check_bridges_do_not_overlap,
        check_bridge_image_round_trip,
    ],
    "partition": [
        check_partition_covers_branch,
        check_partition_unique_owner,
        check_partition_singletons,
        check_holes_are_classified,
        check_backward_holes_keep_bridge_side,
        check_direct_hole_side_is_coordinate_side,
        check_direct_hole_opens_inward_pair,
        check_openings_on_own_bridge_row,
        check_openings_missing_only_at_anchor_or_tail,
        check_propagated_holes_land_on_predicted_branch,
        check_propagation_terminates,
        check_no_direct_hole_beyond_fundamental,
        check_reference_pairs_valid,
    ],
    "arrangement": [
        check_arrangement_euler,
        check_arrangement_regions_sound,
        check_arrangement_regions_disjoint,
        check_arrangement_image_of,
        check_arrangement_preimage_inverts_image,
    ],
    "classes": [
        check_every_bridge_in_one_class,
        check_class_orientation_anchor_outward,
        check_row_of_end_is_geometry,
        check_bridge_rows_consistent,
        check_anchor_bridge_class_leads_its_tangle,
        check_class_table_order,
        check_classes_do_not_mix_tangles,
        check_classes_use_every_orbit_branch,
        check_only_active_classes_lettered,
    ],
    "iterated": [
        check_iterated_child_inside_parent,
        check_iterated_keeps_homotopy_boundaries,
        check_iterated_unique_owner,
        check_iterated_cut_provenance,
        check_names_agree_with_structure,
        check_minimal_trellis_bridges,
        check_minimal_trellis_nodes_and_arrangement,
    ],
    "dual_graph": [
        check_dual_node_structure,
        check_dual_face_nodes,
        check_dual_face_side_is_geometry,
        check_dual_unified_is_pip_segment,
        check_dual_bipartite_degree,
    ],
    "symbolic": [
        check_itineraries_even,
        check_itinerary_pairs_same_side,
        check_itinerary_pairs_are_classes,
        check_landings_contained,
        check_refined_children_inherit_word,
        check_member_matching_consistent,
        check_inert_classes_outside_transitions,
        check_matrix_is_token_counts,
    ],
}

#: The layers whose laws mutate their session: run on ONE fresh build per case,
#: never the shared one (the checks run in order on that build).
MUTATING_LAYERS: dict[str, list[LawCheck]] = {
    "recompute": [check_recompute_preserves_ids, check_recompute_is_idempotent],
}

#: The two laws that mutate their session.
MUTATING: tuple[LawCheck, ...] = tuple(MUTATING_LAYERS["recompute"])


