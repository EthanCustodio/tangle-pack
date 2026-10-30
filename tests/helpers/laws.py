"""The physical-law checks of the law tier, one function per law.

Every check has the signature ``check_<law>(case) -> int``: it asserts its law
on a built :class:`cases.Case` and returns HOW MANY ITEMS it checked, so the
law tier can insist that a law was not vacuous (a count of 0 is a failure
unless the ``(law, case)`` pair is listed in ``cases.NOT_APPLICABLE``). The
checks only read the case; they never mutate the session (the law tier shares
one build per case per run and fingerprints it around every test). The two laws
that must mutate a session (a recompute) share ONE fresh build per case, the
``recompute`` layer of :data:`MUTATING_LAYERS`:
:func:`check_recompute_keeps_the_topology`, :func:`check_recompute_preserves_ids`
then :func:`check_recompute_is_idempotent`.

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
  scaling of a link against ``Trellis.SCALING_RTOL``, both widened by the
  linked crossings' own BRACKETS (the node cdists either side of each crossing
  on both curves: a crossing's cdist is known only to within them, WP2), the
  area of a region image against ``Arrangement.image_of``'s ``rtol``, and every
  cdist comparison against the registry's ``cdist_tol``.
* The iterate-table laws (WP2) recompute every bracket from the curve nodes
  and use the REAL map as the oracle of a link: ``|f(x) - y|`` within the image
  of x's unstable segment plus y's own.
* The geometric oracles here (crossing sign from the manifold nodes, the row
  of a bridge end, a direct hole's side of its bridge, the side of a face) are
  measured from the curves with public geometry helpers and the ``_side_of``
  kernel, independently of the combinatorial rule the library uses.
"""

from __future__ import annotations

import inspect
from collections import defaultdict
from types import SimpleNamespace
from typing import TYPE_CHECKING, Callable, Optional

import numpy as np
import pytest

from helpers.invariants import (
    assert_cdist_monotonic,
    assert_no_geometric_spikes,
    assert_one_to_one,
)
from tanglepack.numerics.BranchPoint import BranchPoint
from tanglepack.numerics.geometry import arc_polyline, oriented_bridge_polyline, point_in_polygon
from tanglepack.topology.Arrangement import Arrangement
from tanglepack.topology.BridgeClass import (
    anchor_outward_key,
    oriented_class,
)
from tanglepack.topology.DualGraph import face_side
from tanglepack.topology.DualWalk import trellis_itinerary, walk_chords
from tanglepack.topology.MinimalTrellis import image_chain
from tanglepack.topology.PartitionFamily import IteratedHomotopyPartition, _hole_stretches
from tanglepack.topology.Pseudoneighbor import compute_pseudoneighbors
from tanglepack.topology.StablePartition import (
    bridge_for_pair,
    owns_cdist,
    row_of_end,
    rows_of_bridge,
)
from tanglepack.topology.SymbolicDynamics import translate_itinerary
from tanglepack.topology.diagnostics import _side_of, bridge_row_violation
from tanglepack.topology.plotting import hole_point
from tanglepack.topology.Trellis import SCALING_RTOL
if TYPE_CHECKING:  # pragma: no cover - typing only
    from cases import Case
    from tanglepack import TangleSession
    from tanglepack.numerics.Bridge import BridgeId

#: Relative slack of the per-link area check: the library's own iterate
#: identification slack (``collision_rtol`` default of ``StrongPip`` and
#: ``Pseudoneighbor``). Measured worst link at baseline: 1.4e-3 (p3, nested),
#: 2.3e-4 (k28 two blasts), 1.8e-4 (k10).
AREA_LINK_RTOL: float = inspect.signature(compute_pseudoneighbors).parameters[
    "collision_rtol"
].default

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
def _is_anchor(ix) -> bool:
    """A synthetic crossing at canonical distance exactly (0, 0): a periodic point."""
    return (
        ix.seg_ids is None
        and ix.unstable_cdist is not None
        and ix.stable_cdist is not None
        and float(ix.unstable_cdist) == 0.0
        and float(ix.stable_cdist) == 0.0
    )


def _anchor_ids(case: "Case") -> list[int]:
    """Registry ids of the synthetic anchors (both cdists exactly 0, no segments)."""
    return [intersection_id for intersection_id, ix in case.registry if _is_anchor(ix)]


def _stable_branch_at(case: "Case", bridge_id, intersection_id: int):
    """The stable branch a bridge meets at one end, written independently of the
    library: the crossing's own, or at an anchor (one point on every branch, D4)
    ``s(o, b)`` of the bridge's ``u(o, b)`` (C7), read off the bridge's other end."""
    ix = case.registry[intersection_id]
    if not _is_anchor(ix):
        return ix.manifold_b_key
    other = bridge_id[1] if bridge_id[0] == intersection_id else bridge_id[0]
    fp, _, orbit_index, branch_index = case.registry[other].manifold_a_key
    return (fp, "stable", orbit_index, branch_index)


def _anchor_row(case: "Case", bridge, stable_key) -> str:
    """The row of a bridge at its anchor on one stable branch, from the signs alone:
    the anchor carries the sign of its home pair ``(u(o, 0), s(o, 0))``, and each
    branch index off that pair (the bridge's own, the stable one) runs the other way."""
    anchor = case.registry[bridge.id[0]]
    flips = (bridge.manifold_key[3] != anchor.manifold_a_key[3]) + (
        stable_key[3] != anchor.manifold_b_key[3]
    )
    return "left" if anchor.crossing_sign * (-1) ** flips > 0 else "right"


def _tol(case: "Case") -> float:
    """The registry's cdist tolerance."""
    return float(case.registry.cdist_tol)


def _full(case: "Case"):
    """The all-fixed-points trellis (cached by the session)."""
    return case.session.trellis()


def _per_fp_trellises(case: "Case") -> list:
    """The per-fixed-point trellises, outermost fixed point first."""
    return [case.session.trellis(fp) for fp in case.fixed_points]


def _staged(case: "Case"):
    """The trellis stages 4-6 ran on: the topology's, the selection's (WP9)."""
    return case.session.topology().trellis


def _holes(case: "Case") -> list:
    """``(trellis, hole)`` for every hole of the topology."""
    trellis = _staged(case)
    return [(trellis, hole) for hole in trellis.holes]


def _partitions(case: "Case") -> list:
    """``(trellis, result)`` for every stable partition of the topology."""
    trellis = _staged(case)
    return [(trellis, result) for result in trellis.stable_partitions]


def _identified_bridges(case: "Case") -> list:
    """The workbench's bridges that carry a ``BridgeId``."""
    return [bridge for bridge in case.workbench.bridges if bridge.id is not None]


def _stable_nodes(manifold) -> tuple[np.ndarray, np.ndarray]:
    """``(points, cdists)`` of a stable manifold, root (anchor) first."""
    points = np.asarray(manifold.points, dtype=float).reshape(-1, 2)
    cdists = np.array([np.nan if c is None else float(c) for c in manifold.cdists])
    return points, cdists


def _local_direction(manifold, stability: str, cdist: float) -> np.ndarray:
    """The curve direction at ``cdist``, taken in increasing canonical distance."""
    points = np.asarray(manifold.points, dtype=float).reshape(-1, 2)
    cdists = list(manifold.cdists)
    index = int(np.searchsorted(cdists, cdist))
    index = min(max(index, 1), len(points) - 1)
    return points[index] - points[index - 1]


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
    partitioned = {result.branch_key[0] for _trellis, result in _partitions(case)}
    assert partitioned == set(case.fixed_points), "every fixed point is partitioned"
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


def check_case_area_cutoff(case: "Case") -> int:
    """The build refined at its frozen cutoff, never at the library default.

    The workbench and its initializer share one machine, so one value covers
    every refinement of the case, the initial fundamental segments included.
    """
    assert case.session.area_cutoff == case.expect.area_cutoff
    return 1


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


def check_fixed_point_minimal_period(case: "Case") -> int:
    """Every fixed point's ``period`` is the minimal period of its coordinates.

    A seed that slid to a lower-period orbit would repeat its points and carry a
    wrong ``k_value``; ``construct_fixed_point`` refuses it, this re-checks the
    fixtures independently of that guard.
    """
    for fp in case.fixed_points:
        orbit = np.asarray(fp.coordinates, dtype=float)
        assert orbit.shape == (fp.period, 2)
        scale = max(1.0, float(np.max(np.abs(orbit))))
        for d in range(1, fp.period):
            if fp.period % d == 0:
                gap = float(np.max(np.abs(np.roll(orbit, -d, axis=0) - orbit)))
                assert gap > 1e-6 * scale, f"{fp.label}: period {fp.period} repeats every {d}"
    return len(case.fixed_points)


def check_map_is_area_preserving(case: "Case") -> int:
    """``|det J| = 1`` at every orbit point, so ``construct_fixed_point`` never warned.

    Measured by central differences of the map itself (exact up to rounding for
    a quadratic map), since not every case hands the workbench a Jacobian; the
    library's own guard must agree wherever it has one to read.
    """
    system = case.workbench.dynamical_system
    points = [np.asarray(p, float) for fp in case.fixed_points for p in fp.coordinates]
    h = 1e-6
    for point in points:
        columns = [
            (np.asarray(system.map(point + h * e)) - np.asarray(system.map(point - h * e)))
            / (2 * h)
            for e in np.eye(2)
        ]
        determinant = np.linalg.det(np.column_stack(columns))
        assert abs(abs(determinant) - 1) < 1e-6, f"|det J| = {abs(determinant)} at {point}"
        assert system.check_area_preserving(point) is not False
    return len(points)


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


def check_case_heteroclinic(case: "Case") -> int:
    """The selection trellis holds a crossing of one fixed point's unstable branch and
    another's stable branch exactly when the case says it is heteroclinic."""
    trellis = _full(case)
    heteroclinic = any(
        ix.manifold_a_key[0] is not ix.manifold_b_key[0]
        for iid, ix in trellis.registry
        if iid in trellis and not ix.is_anchor
    )
    assert heteroclinic == case.expect.heteroclinic
    return 1


# --------------------------------------------------------------------------- #
# Manifolds
# --------------------------------------------------------------------------- #
def check_manifold_cdist_monotone(case: "Case") -> int:
    """cdist is non-decreasing along every manifold and every bridge (ties allowed),
    and where two neighbours tie as plain floats their residuals put them in order.

    Not "the precise cdist never decreases": the merge orders a chain by the
    plain floats, and the k=2.8 and inversion cases hold neighbours one ULP
    apart whose residuals order them the other way (2026-10-30).
    """
    from tanglepack.numerics.precise import compare

    curves = list(case.workbench.manifolds.values()) + list(case.workbench.bridges)
    for curve in curves:
        assert_cdist_monotonic(curve, strict=False)
        _head, store, rows = curve._rows(stop_at_final=False)
        hi, lo = store.precise_cdist(rows)
        tied = np.flatnonzero(store.cdist[rows][1:] == store.cdist[rows][:-1])
        assert (compare(hi[tied + 1], lo[tied + 1], hi[tied], lo[tied]) > 0).all(), (
            curve.manifold_key, rows[tied].tolist())
    return len(curves)


def check_manifold_no_spikes(case: "Case") -> int:
    """No manifold carries a geometric spike (a scrambled node)."""
    manifolds = list(case.workbench.manifolds.values())
    for manifold in manifolds:
        assert_no_geometric_spikes(manifold)
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
    unstable = [m for key, m in case.workbench.manifolds.items() if key[1] == "unstable"]
    manifold_rows = set(np.concatenate([m.ids for m in unstable]).tolist()) if unstable else set()
    heads = {id(m.root) for m in unstable if isinstance(m.root, BranchPoint)}
    return [
        bridge
        for bridge in case.workbench.bridges
        if not (manifold_rows.intersection(bridge.ids.tolist())
                or (isinstance(bridge.root, BranchPoint) and id(bridge.root) in heads))
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
        lo, hi = sorted(float(c) for c in ix.segment_ends("unstable")[1])
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


def _brackets(ix) -> tuple[tuple[float, float], tuple[float, float]]:
    """A detected crossing's unstable and stable brackets, from the curve nodes.

    The crossing was found between two adjacent nodes of each curve, so each of
    its canonical distances is known only to within those nodes' cdists.
    """
    unstable, stable = (
        sorted(float(c) for c in ix.segment_ends(stability)[1]) for stability in ("unstable", "stable")
    )
    return (unstable[0], unstable[1]), (stable[0], stable[1])


def _meets(lo: float, hi: float, other_lo: float, other_hi: float) -> bool:
    """Two closed intervals meet, up to float rounding (1e-12 relative)."""
    pad = 1e-12 * max(abs(lo), abs(hi), abs(other_lo), abs(other_hi))
    return other_lo <= hi + pad and lo <= other_hi + pad


def _image_brackets(case: "Case", ix) -> tuple[tuple, tuple, tuple[float, float], tuple[float, float]]:
    """Where the exact rule puts a crossing's image: the advanced keys and its
    brackets carried one map step by the two manifolds' own ``stretch_param``."""
    (u_lo, u_hi), (s_lo, s_hi) = _brackets(ix)
    stretch_u = ix.unstable_manifold.stretch_param
    stretch_s = case.workbench.manifolds[ix.manifold_b_key].stretch_param
    a_key, b_key = (key[0].advance_key(key, 1) for key in (ix.manifold_a_key, ix.manifold_b_key))
    return a_key, b_key, (u_lo * stretch_u, u_hi * stretch_u), (s_lo / stretch_s, s_hi / stretch_s)


def _passes_map_gate(case: "Case", source, image) -> bool:
    """``|f(x) - y|`` is at most the image of x's unstable segment plus y's own segment."""
    def mapped(point) -> np.ndarray:
        return np.ravel(case.workbench.dynamical_system.map(np.asarray(point, dtype=float)))

    p0, p1 = source.segment_ends("unstable")[0]
    q0, q1 = image.segment_ends("unstable")[0]
    gate = np.linalg.norm(mapped(p1) - mapped(p0)) + np.linalg.norm(np.ravel(q1) - np.ravel(q0))
    return float(np.linalg.norm(mapped(source.coords) - np.ravel(image.coords))) <= gate


def _detected_links(case: "Case") -> list[tuple[int, int]]:
    """``(x, f(x))`` for every registered ``+1`` link of a detected (non-anchor) crossing."""
    registry = case.registry
    return [
        (source, target)
        for source, n, target in registry.iterate_table.items()
        if n == 1 and not registry[source].is_synthetic
    ]


def check_area_along_iterate_links(case: "Case", *, rtol: float = AREA_LINK_RTOL) -> int:
    """``unstable_cdist * stable_cdist`` is preserved by every registered homoclinic
    ``+1`` link, and scaled by ``beta_A / beta_B`` on a heteroclinic
    ``W^u(A) x W^s(B)`` link: each side moves on its own saddle (NUM-20).

    Compared LINK BY LINK, not against the head of the chain, so the error of
    a long chain (period 3 runs a dozen links) does not accumulate. The
    anchors (product 0) are skipped. The slack is ``rtol`` plus the relative
    widths of the four brackets the two products are read from.
    """
    registry = case.registry
    count = 0
    for source, target in _detected_links(case):
        ix, image = registry[source], registry[target]
        widths = sum(
            (hi - lo) / float(cdist)
            for crossing in (ix, image)
            for (lo, hi), cdist in zip(
                _brackets(crossing), (crossing.unstable_cdist, crossing.stable_cdist)
            )
        )
        factor = ix.manifold_a_key[0].per_step_beta("unstable") * ix.manifold_b_key[0].per_step_beta("stable")
        before = float(ix.unstable_cdist) * float(ix.stable_cdist) * factor
        after = float(image.unstable_cdist) * float(image.stable_cdist)
        assert after == pytest.approx(before, rel=rtol + widths), (
            f"area not preserved on the link {source} -> {target}: {before!r} -> {after!r}"
        )
        count += 1
    return count


def check_iterate_link_scales_by_beta(case: "Case") -> int:
    """A ``+1`` link scales the unstable cdist by the ``per_step_beta`` of its unstable
    saddle and the stable cdist by that of its stable saddle (NUM-20: on a
    heteroclinic link they differ), and lands on the branches ``advance_key`` names.

    The slack is ``SCALING_RTOL`` of the prediction plus the two crossings'
    bracket widths (the source's carried one step).
    """
    registry = case.registry
    count = 0
    for source, target in _detected_links(case):
        ix, image = registry[source], registry[target]
        a, b = ix.manifold_a_key[0], ix.manifold_b_key[0]
        beta_u, beta_s = a.per_step_beta("unstable"), b.per_step_beta("stable")
        assert image.manifold_a_key == a.advance_key(ix.manifold_a_key, 1), source
        assert image.manifold_b_key == b.advance_key(ix.manifold_b_key, 1)
        (u_lo, u_hi), (s_lo, s_hi) = _brackets(ix)
        (iu_lo, iu_hi), (is_lo, is_hi) = _brackets(image)
        predicted_u = beta_u * float(ix.unstable_cdist)
        predicted_s = beta_s * float(ix.stable_cdist)
        assert abs(float(image.unstable_cdist) - predicted_u) <= (
            SCALING_RTOL * predicted_u + beta_u * (u_hi - u_lo) + (iu_hi - iu_lo)
        ), (source, target)
        assert abs(float(image.stable_cdist) - predicted_s) <= (
            SCALING_RTOL * predicted_s + beta_s * (s_hi - s_lo) + (is_hi - is_lo)
        ), (source, target)
        count += 1
    return count


def check_iterate_table_partial_bijection(case: "Case") -> int:
    """The map is injective: no crossing is the ``+n`` image of two, and every
    forward link reads back as the backward one, ``table[table[x, n], -n] == x``."""
    table = case.registry.iterate_table
    preimage: dict[tuple[int, int], int] = {}
    count = 0
    for source, n, target in table.items():
        assert table[target, -n] == source, (source, n, target)
        assert preimage.setdefault((target, n), source) == source, (
            f"crossing {target} is the +{n} image of both {preimage[target, n]} and {source}"
        )
        count += 1
    return count


def check_iterate_links_pass_the_map_gate(case: "Case") -> int:
    """The real map confirms every ``+1`` link of a detected crossing: ``f(x)``
    lies within the image of x's unstable segment plus y's own segment of y."""
    registry = case.registry
    count = 0
    for source, target in _detected_links(case):
        assert _passes_map_gate(case, registry[source], registry[target]), (
            f"the link {source} -> {target} fails the map gate"
        )
        count += 1
    return count


def check_iterate_table_complete(case: "Case") -> int:
    """No unlinked detected crossing has an image the exact rule would link.

    A crossing y is the image of x when it sits on the advanced branches, both
    of its brackets meet x's carried one step, and the real map puts it within
    the gate. The table must already hold every such link (inference runs
    after every detection and every blast, BCG-3), except for the sources the
    rule cannot decide (``registry.unresolved_links``, B8): those hold at least
    one candidate, none is linked, and an image contested alone by them has
    no preimage either -- no claimant is chosen by the order of the passes.
    """
    registry = case.registry
    table = registry.iterate_table
    unresolved = registry.unresolved_links
    for source, candidates in unresolved.items():
        assert candidates, f"crossing {source} is unresolved with no candidate image"
        assert table[source, 1] is None, f"crossing {source} is unresolved but linked"
        if len(candidates) == 1:
            assert table[candidates[0], -1] is None, (
                f"image {candidates[0]} is contested by {source} but linked to "
                f"{table[candidates[0], -1]}"
            )
    on_branches: dict[tuple, list[int]] = defaultdict(list)
    for intersection_id, ix in registry:
        if not ix.is_synthetic:
            on_branches[ix.manifold_a_key, ix.manifold_b_key].append(intersection_id)
    count = 0
    for intersection_id, ix in registry:
        if ix.is_synthetic:
            continue
        count += 1
        if table[intersection_id, 1] is not None or intersection_id in unresolved:
            continue
        a_key, b_key, unstable, stable = _image_brackets(case, ix)
        images = [
            candidate
            for candidate in on_branches[a_key, b_key]
            if _meets(*unstable, *_brackets(registry[candidate])[0])
            and _meets(*stable, *_brackets(registry[candidate])[1])
            and _passes_map_gate(case, ix, registry[candidate])
        ]
        assert not images, f"crossing {intersection_id} is unlinked, but {images} is its image"
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


def _topology_snapshot(case: "Case") -> tuple:
    """What a recompute must not move: the bridges and their ``iterated`` flags, the
    number of partial pieces, the holes, the class table and the dynamics."""
    session, workbench = case.session, case.workbench
    bridges = sorted((bridge.id, bridge.iterated) for bridge in workbench.bridges if not bridge.partial)
    partials = sum(bridge.partial for bridge in workbench.bridges)
    holes = sorted((hole.bounding_ids, hole.iterate, hole.bridge_side) for hole in session.topology().holes)
    return (
        bridges,
        partials,
        holes,
        session.bridge_classes().describe(),
        session.symbolic_dynamics().describe(),
    )


def check_recompute_keeps_the_topology(case: "Case") -> int:
    """MUTATES (fresh build only), D1: a recompute of the crossings and a rebuild of
    the bridges, with nothing grown or blasted in between, change nothing.

    The blast's image curves are computed manifold (NUM-12): their crossings,
    bridges, partial pieces and ``iterated`` flags survive, so the bridges, holes,
    class table and words are exactly what they were.
    """
    before = _topology_snapshot(case)
    case.workbench.compute_intersections(case.fixed_points, preserve_ids=True)
    case.workbench.rebuild_bridges()
    after = _topology_snapshot(case)
    for name, old, new in zip(("bridges", "partial pieces", "holes", "classes", "dynamics"), before, after):
        assert new == old, f"a recompute moved the {name}"
    return len(before[0])


def check_recompute_preserves_ids(case: "Case") -> int:
    """MUTATES (fresh build only): ``compute_intersections(preserve_ids=True)``
    gives every crossing back under its own id, keys and cdists, the crossings
    on a blast's image curves included (NUM-12)."""
    tol = _tol(case)
    before = _crossing_fingerprint(case)
    case.workbench.compute_intersections(case.fixed_points, preserve_ids=True)
    after = _crossing_fingerprint(case)
    assert set(after) == set(before), (
        f"lost {sorted(set(before) - set(after))}, new {sorted(set(after) - set(before))}"
    )
    for intersection_id, value in after.items():
        assert _same_crossing(before[intersection_id], value, tol), intersection_id
    return len(after)


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
def check_one_anchor_per_periodic_point(case: "Case") -> int:
    """D4: each periodic point is ONE anchor, and every branch attached to it (both
    stabilities) starts there (replaces ``one_anchor_per_unstable_branch``)."""
    at_point: dict = defaultdict(list)
    for intersection_id in _anchor_ids(case):
        fp, _, orbit_index, _ = case.registry[intersection_id].manifold_a_key
        at_point[id(fp), orbit_index].append(intersection_id)
    count = 0
    for key, branch in _full(case).branches.items():
        ids = at_point.get((id(key[0]), key[2]), [])
        assert len(ids) == 1, f"{len(ids)} anchors at the point of branch {key[1:]}: {ids}"
        assert branch.intersection_ids[0] == ids[0], (key[1:], branch.intersection_ids[:2])
        count += 1
    return count


def check_anchor_rays_are_its_branches(case: "Case") -> int:
    """D4: the real arcs of the arrangement leaving an anchor are exactly one per
    trellis branch of its point that holds a crossing past it (one degree-4 node at
    an inversion point, never four coincident ones)."""
    trellis = _full(case)
    leaving: dict = defaultdict(set)
    for face in case.session.arrangement().faces:
        for arc in face.arcs:
            if arc is None:
                continue
            for end in (arc.lo_id, arc.hi_id):
                leaving[end].add((arc.kind, arc.branch_key, frozenset((arc.lo_id, arc.hi_id))))
    count = 0
    for intersection_id in _anchor_ids(case):
        fp, _, orbit_index, _ = case.registry[intersection_id].manifold_a_key
        expected = {
            key
            for key, branch in trellis.branches.items()
            if key[0] is fp and key[2] == orbit_index and len(branch) >= 2
        }
        rays = leaving[intersection_id]
        assert {branch_key for _kind, branch_key, _ends in rays} == expected, intersection_id
        assert len(rays) == len(expected), (intersection_id, sorted(r[1][1:] for r in rays))
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
def check_bridge_id_in_unstable_order(case: "Case") -> int:
    """``BridgeId = (first, second)`` is its endpoints in strictly increasing unstable cdist."""
    registry = case.registry
    count = 0
    for bridge in _identified_bridges(case):
        assert bridge.id == (bridge.first_intersection, bridge.second_intersection)
        first, second = (registry[i] for i in bridge.id)
        assert float(first.unstable_cdist) < float(second.unstable_cdist), bridge.id
        count += 1
    return count


def check_bridge_endpoints_on_own_branch(case: "Case") -> int:
    """A bridge's two crossings lie on ITS unstable branch, inside its cdist span; an
    anchor end is the anchor of the bridge's own periodic point (one point on every
    branch, D4).

    The codebase-audit regression (2026-09): a registry-wide nearest-cdist
    lookup handed the first bridge of one period-3 branch another branch's
    anchor.
    """
    registry = case.registry
    tol = _tol(case)
    count = 0
    for bridge in _identified_bridges(case):
        first, second = (registry[i] for i in bridge.id)
        root_u, tail_u = (
            float(end.get_cdist("unstable") if isinstance(end, BranchPoint) else end.cdist)
            for end in (bridge.root, bridge.tail)
        )
        lo, hi = sorted((root_u, tail_u))
        for end in (first, second):
            if _is_anchor(end):
                fp, _, orbit_index, _ = end.manifold_a_key
                assert (fp, orbit_index) == (bridge.manifold_key[0], bridge.manifold_key[2])
                assert end is first, bridge.id
            else:
                assert end.manifold_a_key == bridge.manifold_key, (
                    f"bridge {bridge.id} on {bridge.manifold_key[1:]} got an endpoint "
                    f"from {end.manifold_a_key[1:]}"
                )
            assert _within(float(end.unstable_cdist), lo, hi, tol), bridge.id
            count += 1
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
def _sign_openings(case: "Case", trellis, hole) -> set:
    """The ``(intersection_id, half, row, branch_key)`` halves a hole opens, derived from
    its bounds, its ``bridge_side`` and the registry signs alone (D6, never from
    ``hole.openings``): at a crossing the outward half iff ``(sign > 0) == (side ==
    "left")``, on the crossing's branch and the bridge's row there; at an anchor the
    outward half of the branch on which the bridge's row is its side, if computed."""
    bridge = trellis.bridge_between(*hole.bounding_ids)
    assert bridge is not None and bridge.id is not None, hole.bounding_ids
    fp, _, orbit_index, _ = bridge.manifold_key
    opened = set()
    for endpoint, intersection_id in zip(("first", "second"), bridge.id):
        ix = case.registry[intersection_id]
        if _is_anchor(ix):
            for branch_index in (0, 1):
                key = (fp, "stable", orbit_index, branch_index)
                if key in trellis.branches and _anchor_row(case, bridge, key) == hole.bridge_side:
                    opened.add((intersection_id, "outward", hole.bridge_side, key))
            continue
        positive = ix.crossing_sign > 0
        half = "outward" if positive == (hole.bridge_side == "left") else "anchorward"
        row = "left" if positive == (endpoint == "first") else "right"
        opened.add((intersection_id, half, row, ix.manifold_b_key))
    return opened


def check_partition_covers_branch(case: "Case") -> int:
    """Each side's partition runs anchor to branch end in contiguous, positional elements,
    and a non-singleton element is open at an end iff a hole opens that half there
    (stage 6). The opened halves are derived from each hole's bounds, bridge side and
    crossing signs (:func:`_sign_openings`), never read from ``hole.openings``, which the
    partition is built from. WP10: an element open at both ends need not be ONE hole's
    span -- between consecutive U-arcs it is a stretch of the face holding a hole of the
    other bridge, as on p3 at 16 unstable steps."""
    count = 0
    staged = _staged(case)
    opened = {half for hole in staged.holes for half in _sign_openings(case, staged, hole)}
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
        for iv in intervals:
            if iv.lo_id == iv.hi_id:
                continue  # a pinched singleton is closed at its one crossing
            for iid, closed, half in ((iv.lo_id, iv.closed_lo, "outward"), (iv.hi_id, iv.closed_hi, "anchorward")):
                if iid is not None:
                    is_open = (iid, half, result.side, result.branch_key) in opened
                    assert (not closed) == is_open, (result.branch_key[1:], result.side, iid, half)
        count += 1
    return count


def _owners(intervals, intersection_id: int, cdist: float) -> list[int]:
    """Element ids owning a crossing, from the raw interval fields and with no tolerance:
    the elements closed at its id when it is a boundary, else those covering its exact
    cdist (an end only when closed)."""
    if any(intersection_id in (iv.lo_id, iv.hi_id) for iv in intervals):
        return sorted(
            {
                iv.element_id
                for iv in intervals
                if (iv.lo_id == intersection_id and iv.closed_lo)
                or (iv.hi_id == intersection_id and iv.closed_hi)
            }
        )
    return [
        iv.element_id
        for iv in intervals
        if iv.lo_cdist < cdist < iv.hi_cdist
        or (cdist == iv.lo_cdist and iv.closed_lo)
        or (cdist == iv.hi_cdist and iv.closed_hi)
    ]


def check_partition_unique_owner(case: "Case") -> int:
    """Every crossing of a partitioned branch has exactly one owning element per side,
    by id at a boundary and by exact cdist elsewhere."""
    count = 0
    for trellis, result in _partitions(case):
        branch = trellis.branch(result.branch_key)
        for intersection_id in branch.intersection_ids:
            cdist = float(trellis.intersection(intersection_id).stable_cdist)
            owners = _owners(result.intervals, intersection_id, cdist)
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
    """Every hole has a bridge side, its anchorward bound first, a finite drawing point and an orbit.

    A direct hole's side comes from its crossing signs, a propagated hole's is
    carried along its orbit (author, 2026-10-07).
    """
    count = 0
    for trellis, hole in _holes(case):
        assert hole.bridge_side in ("left", "right")
        near, far = hole.bounding_ids
        assert hole.near_intersection_id == near
        assert float(trellis.intersection(near).stable_cdist) <= float(
            trellis.intersection(far).stable_cdist
        ) or trellis.intersection(near).manifold_b_key != trellis.intersection(far).manifold_b_key
        assert np.isfinite(hole_point(trellis, hole)).all()
        assert hole.iterate is not None and hole.origin is not None
        count += 1
    return count


def check_backward_holes_keep_bridge_side(case: "Case") -> int:
    """I1 (author rule, 2026-10-05): a hole maps backward onto the same side of the
    BRIDGE it is punched in -- never "the same side of the stable manifold".

    Every hole of an origin carries the ``bridge_side`` of the origin's
    reference hole (iterate 0), parity-flipped under an orientation-reversing
    map; and, measured on the curves, the reference hole's point carried
    back ``|iterate|`` steps by the real inverse map sit on that side of the
    propagated hole's own bridge. Nothing is asserted about the stable-manifold
    ROW of a backward hole: it may differ from its origin's, and that is
    correct.
    """
    inverse = case.workbench.dynamical_system.map_inv
    trellis = _staged(case)
    preserving = trellis.orientation_preserving
    reference = {
        hole.origin: hole
        for hole in trellis.holes
        if hole.pair is not None and hole.iterate == 0
    }
    count = 0
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
            point = hole_point(trellis, origin)
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
    """A direct hole's ``bridge_side`` is the side of its bridge its drawing point sits on.

    The combinatorial crossing-sign rule (2026-10-02) must reproduce the
    geometric measurement hole for hole.
    """
    count = 0
    for trellis, hole in _holes(case):
        if hole.pair is None:
            continue
        bridge = bridge_for_pair(trellis, hole.pair)
        assert bridge is not None, hole.pair
        measured = _point_side_of_bridge(case, bridge, hole_point(trellis, hole))
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
        for intersection_id, which, row, key in hole.openings:
            assert which == ("outward" if intersection_id == near else "anchorward")
            assert row == _opening_row(case, trellis, bridge, intersection_id, key), (
                hole.bounding_ids,
                row,
            )
        count += 1
    return count


def _opening_row(case: "Case", trellis, bridge, intersection_id: int, key) -> str:
    """The row an opening at one bound must sit on: its own bridge's row there, on the
    branch the record opens -- the crossing's own branch, or at an anchor the one the
    signs pick (:func:`_anchor_row`)."""
    if _is_anchor(case.registry[intersection_id]):
        return _anchor_row(case, bridge, key)
    assert key == case.registry[intersection_id].manifold_b_key, (bridge.id, intersection_id)
    endpoint = "first" if intersection_id == bridge.id[0] else "second"
    return row_of_end(trellis, bridge.id, endpoint)


def check_anchor_opening_faces_its_hole(case: "Case") -> int:
    """D6 at the anchor: a hole opening there opens the stable branch whose first node
    past the anchor lies on the hole's side of its own oriented bridge (geometry as
    the independent check of the sign rule)."""
    count = 0
    for trellis, hole in _holes(case):
        bridge = trellis.bridge_between(*hole.bounding_ids)
        for intersection_id, which, _row, key in hole.openings or []:
            if not _is_anchor(case.registry[intersection_id]):
                continue
            assert which == "outward", hole.bounding_ids
            first_point = np.asarray(case.workbench.manifolds[key].get_point_array()[1], dtype=float)
            assert _point_side_of_bridge(case, bridge, first_point) == hole.bridge_side, (
                hole.bounding_ids,
                key[1:],
            )
            count += 1
    return count


def check_propagated_holes_land_on_predicted_branch(case: "Case") -> int:
    """A hole at iterate ``n`` lands in a bridge on ``cycle[(pos_ref + n) % k]``, the
    unstable branch cycle position of its reference bridge advanced by ``n``, the cycle
    being the reference BRIDGE's fixed point's (A's for a heteroclinic pair on B, D7)."""
    trellis = _staged(case)
    reference_position = {}
    for pair in trellis.reference_pseudoneighbors:
        bridge = bridge_for_pair(trellis, pair)
        assert bridge is not None
        cycle = bridge.fixed_point.branch_cycle("unstable")
        reference_position[pair.as_tuple()] = (cycle, cycle.index(bridge.manifold_key))
    count = 0
    for hole in trellis.holes:
        if hole.pair is not None:
            continue
        cycle, position = reference_position[hole.origin]
        expected = cycle[(position + hole.iterate) % len(cycle)]
        bridge = trellis.bridge_between(*hole.bounding_ids)
        assert bridge.manifold_key == expected, (
            f"hole of origin {hole.origin} at iterate {hole.iterate} landed on "
            f"{bridge.manifold_key[1:]}, expected {expected[1:]}"
        )
        count += 1
    return count


def check_propagation_terminates(case: "Case") -> int:
    """No origin punches the same bounding region twice (propagation stops at periodicity)."""
    propagated = [hole for hole in _staged(case).holes if hole.iterate not in (0, None)]
    count = 0
    for origin in {hole.origin for hole in propagated}:
        regions = [frozenset(hole.bounding_ids) for hole in propagated if hole.origin == origin]
        assert len(regions) == len(set(regions)), origin
        count += 1
    return count


def check_no_forward_hole(case: "Case") -> int:
    """Holes propagate backward only: no hole sits at a forward iterate.

    On a period-k orbit too (author, 2026-10-09, reversing D5): the
    ``+1 .. +(k_value - 1)`` appearances on the other branches are recorded
    pairs, never holes.
    """
    trellis = _staged(case)
    count = 0
    for hole in trellis.holes:
        assert hole.iterate <= 0, (hole.bounding_ids, hole.iterate)
        count += 1
    for pair in trellis.pseudoneighbors:
        if pair.iterate is not None and pair.iterate > 0:
            assert pair.hole is None, pair.as_tuple()
            count += 1
    return count


def pip_orbit_cuts(trellis, pip: int) -> dict:
    """
    The cuts of the trellis of the pip ``T(q0) = (W^u[P, f^m(q0)], W^s[P, q0])`` (D10, D12),
    read off the iterate table's ``+1`` links and the registry alone.

    Stable branch ``j`` ends at ``f^j(q0)``, ``j < k``; the unstable branch of
    ``f^i(q0)`` ends there for the ``k`` deepest ``i``, ``m`` the last ``+1`` link of
    the pip's chain. A point the table does not hold is scaled from the pip by
    ``per_step_beta`` and flagged as placed.

    Args:
        trellis: A trellis holding the pip's registry.
        pip: The strong pip's id.

    Returns:
        ``{branch key: (cdist, placed)}``.
    """
    orbit = [pip]
    while (image := trellis.iterate(orbit[-1], 1)) is not None:
        orbit.append(image)
    q0 = trellis.intersection(pip)
    fp = q0.manifold_b_key[0]
    k, m = fp.k_value, len(orbit) - 1

    def point(i: int, stability: str) -> tuple:
        key = q0.manifold_b_key if stability == "stable" else q0.manifold_a_key
        if 0 <= i <= m:
            ix = trellis.intersection(orbit[i])
            own = ix.manifold_b_key if stability == "stable" else ix.manifold_a_key
            cdist = ix.stable_cdist if stability == "stable" else ix.unstable_cdist
            return own, (float(cdist), False)
        cdist = q0.stable_cdist if stability == "stable" else q0.unstable_cdist
        return fp.advance_key(key, i), (float(cdist) * fp.per_step_beta(stability) ** i, True)

    return dict(
        [point(j, "stable") for j in range(k)]
        + [point(i, "unstable") for i in range(m - k + 1, m + 1)]
    )


def check_trellis_crossings_belong_to_its_selection(case: "Case") -> int:
    """A trellis holds a crossing only when BOTH its branches are selected (an anchor
    lies on its own branches), and each of its bridges ends at crossings it holds: a
    tangle alone holds no heteroclinic crossing, its arcs run through them (WP9)."""
    count = 0
    for trellis in [*_per_fp_trellises(case), _full(case)]:
        for iid in trellis.own_intersection_ids:
            ix = trellis.intersection(iid)
            if not ix.is_anchor:
                assert {ix.manifold_a_key, ix.manifold_b_key} <= set(trellis.branches), iid
            count += 1
        for bridge in trellis.bridges:
            ends = [i for i in (bridge.first_intersection, bridge.second_intersection)
                    if i is not None]
            assert all(i in trellis for i in ends), (bridge.id, ends)
            count += 1
    return count


def check_trellis_is_cut_at_the_pip_orbit(case: "Case") -> int:
    """D10, D12: once a pip is chosen, the trellis holds exactly the crossings of ``T(q0)``
    in ``N``: none past a pip-orbit cut on its stable or its unstable branch, and every
    resolved crossing of its branches inside both cuts. The selection trellis is cut at
    every pip. A detection the polylines cannot resolve is in no trellis (Q13 a, B1)."""
    count = 0
    trellises = _per_fp_trellises(case)
    pips = [trellis.strong_pip for trellis in trellises]
    for trellis, chosen in [*((t, [t.strong_pip]) for t in trellises), (_full(case), pips)]:
        cuts = {
            key: cut
            for pip in chosen
            if pip is not None
            for key, cut in pip_orbit_cuts(trellis, pip).items()
        }
        tol = float(trellis.registry.cdist_tol)
        held = set(trellis.own_intersection_ids)
        for iid, ix in trellis.registry:
            if (
                ix.is_anchor
                or iid not in trellis.resolved_ids
                or not {ix.manifold_a_key, ix.manifold_b_key} <= set(trellis.branches)
            ):
                continue
            inside = all(
                float(cdist) <= cuts[key][0] + (tol if cuts[key][1] else 0.0)
                for key, cdist in (
                    (ix.manifold_b_key, ix.stable_cdist),
                    (ix.manifold_a_key, ix.unstable_cdist),
                )
                if key in cuts
            )
            assert (iid in held) == inside, (
                iid, iid in held, float(ix.unstable_cdist), float(ix.stable_cdist))
            count += 1
    return count


def check_reference_pairs_valid(case: "Case") -> int:
    """A reference pair is consecutive on its stable AND its unstable branch, on the
    strong pip's own branch, inside the reference window ``[f^k(q0), q0]``."""
    trellis = _staged(case)
    tol = float(trellis.registry.cdist_tol)
    count = 0
    for fp in case.fixed_points:
        pip = trellis.pips[fp]
        pip_ix = trellis.intersection(pip)
        top = float(pip_ix.stable_cdist)
        _key, bottom, _from_table = trellis.image_cdist(pip, fp.k_value, "stable")
        for pair in trellis.reference_pseudoneighbors:
            if pair.branch_key[0] is not fp:
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


def check_chosen_pip_is_strong(case: "Case") -> int:
    """The chosen pip passes the box test (A.4), re-derived from the registry and the
    fixed point alone: no crossing of the pip's fixed point, carried onto the pip's
    stable branch by the fewest backward steps, lies in the open box
    ``(0, c_s(q0)) x (0, c_u(q0))``. The whole registry is read, since classification
    un-chooses the pip first; q0's own orbit is recognised by a 2-D collision within
    ``SCALING_RTOL``, the one identity assumption."""
    trellis, count = _staged(case), 0
    for fp in case.fixed_points:
        q0 = case.registry[trellis.pips[fp]]
        cycle = fp.branch_cycle("stable")
        s0, u0 = float(q0.stable_cdist), float(q0.unstable_cdist)
        for rid, r in case.registry:
            if r.is_anchor or r.manifold_a_key[0] is not fp or r.manifold_b_key not in cycle:
                continue  # homoclinic crossings of fp only: a heteroclinic u-cdist is another saddle's
            steps = (cycle.index(r.manifold_b_key) - cycle.index(q0.manifold_b_key)) % len(cycle)
            s = float(r.stable_cdist) / fp.per_step_beta("stable") ** steps  # backward: stable grows
            u = float(r.unstable_cdist) / fp.per_step_beta("unstable") ** steps  # unstable shrinks
            count += 1
            if abs(s - s0) <= SCALING_RTOL * s0 and abs(u - u0) <= SCALING_RTOL * u0:
                continue
            assert not (0 < s < s0 and 0 < u < u0), (rid, steps, s, s0, u, u0)
    return count


def _key(ix, stability: str):
    """A crossing's manifold key of one stability."""
    return ix.manifold_a_key if stability == "unstable" else ix.manifold_b_key


def _cdist(ix, stability: str) -> float:
    """A crossing's canonical distance along its manifold of one stability."""
    return float(ix.unstable_cdist if stability == "unstable" else ix.stable_cdist)


def _scaled(ix, n: int) -> dict:
    """``M^n`` of a crossing placed by scaling alone: per stability, the advanced key and
    the cdist times ``per_step_beta ** n`` of THAT side's fixed point (D7)."""
    image = {}
    for stability in ("unstable", "stable"):
        key = _key(ix, stability)
        fp = key[0]
        image[stability] = (fp.advance_key(key, n), _cdist(ix, stability) * fp.per_step_beta(stability) ** n)
    return image


def _laps(cdist: float, lo: float, hi: float, beta: float, n0: int, k: int) -> list[int]:
    """Every ``n = n0 (mod k)`` with ``lo < cdist * beta**n < hi``: finitely many, ``beta != 1``."""
    a, b = sorted(np.log([lo / cdist, hi / cdist]) / np.log(beta))
    first = n0 + k * int(np.floor((a - n0) / k))
    return [n for n in range(first, int(np.ceil(b)) + 1, k) if lo < cdist * beta**n < hi]


def _table_image(registry, intersection_id: int, n: int) -> Optional[int]:
    """``M^n`` of a crossing read off the registry's iterate table: the direct entry,
    else the composition of its ``+-1`` links; None where the table stops."""
    if (direct := registry.iterate_table[intersection_id, n]) is not None:
        return direct
    step = 1 if n > 0 else -1
    for _ in range(abs(n)):
        intersection_id = registry.iterate_table[intersection_id, step]
        if intersection_id is None:
            return None
    return intersection_id


def check_reference_pairs_clear_of_X(case: "Case") -> int:
    """Both open intervals of every reference pair hold no point of X (the pseudoneighbor
    definition), by scaling alone: every crossing ``r`` of the cut N on the interval's
    branch cycle, carried ``n`` steps by ``per_step_beta`` (the cycle and beta of the
    interval's own fixed point, D7), for every ``n`` landing on the interval's branch
    strictly inside it. A landing the iterate table links to an end of the pair, or one
    colliding with an end on both cdists within ``SCALING_RTOL``, is that end. Law cases
    only: scaling alone is wrong at depth (A.5)."""
    trellis, registry, count = _staged(case), case.registry, 0
    crossings = [(rid, registry[rid]) for rid in trellis.own_intersection_ids
                 if not registry[rid].is_anchor]
    for pair in trellis.reference_pseudoneighbors:
        ends = pair.as_tuple()
        x0, x1 = (registry[i] for i in ends)
        for stability in ("unstable", "stable"):
            key = _key(x0, stability)
            cycle = key[0].branch_cycle(stability)
            lo, hi = sorted((_cdist(x0, stability), _cdist(x1, stability)))
            for rid, r in crossings:
                if _key(r, stability) not in cycle:
                    continue
                n0 = (cycle.index(key) - cycle.index(_key(r, stability))) % len(cycle)
                beta = key[0].per_step_beta(stability)
                for n in _laps(_cdist(r, stability), lo, hi, beta, n0, len(cycle)):
                    if _table_image(registry, rid, n) in ends or any(
                        all(
                            image_key == _key(x, s)
                            and abs(c - _cdist(x, s)) <= SCALING_RTOL * _cdist(x, s)
                            for s, (image_key, c) in _scaled(r, n).items()
                        )
                        for x in (x0, x1)
                    ):
                        continue  # the pair's own end
                    raise AssertionError(
                        f"pair {ends}: M^{n}({rid}) lies inside its open {stability} interval"
                    )
            count += 1
    return count


def check_propagated_hole_inside_container(case: "Case") -> int:
    """Every propagated hole's container holds the reference bridge's ends carried back
    by the REAL inverse map (the test may call the map; the topology never does):
    each carried end, projected onto the container's polyline, lies strictly between
    the container's two crossings -- or, when the iterate table links that end to a
    crossing of the container, on that shared end (the real map's overshoot of an
    exact crossing is noise)."""
    trellis, registry, count = _staged(case), case.registry, 0
    map_inv = case.workbench.dynamical_system.map_inv_batch
    for hole in trellis.holes:
        if hole.iterate is None or hole.iterate >= 0:
            continue
        container = trellis.bridge_between(*hole.bounding_ids)
        points = np.asarray(container.points, dtype=float).reshape(-1, 2)
        cdists = np.asarray(container.cdists, dtype=float)
        lo, hi = sorted(float(registry[i].unstable_cdist) for i in hole.bounding_ids)
        for end in hole.origin:
            if _table_image(registry, end, hole.iterate) in hole.bounding_ids:
                count += 1
                continue
            point = np.asarray(registry[end].coords, dtype=float).reshape(1, 2)
            for _ in range(-hole.iterate):
                point = map_inv(point)
            a, ab = points[:-1], points[1:] - points[:-1]
            t = np.clip(np.einsum("ij,ij->i", point - a, ab)
                        / np.maximum(np.einsum("ij,ij->i", ab, ab), 1e-300), 0.0, 1.0)
            nearest = int(np.argmin(np.linalg.norm(a + t[:, None] * ab - point, axis=1)))
            cdist = cdists[nearest] + t[nearest] * (cdists[nearest + 1] - cdists[nearest])
            assert lo < cdist < hi, (hole.origin, end, hole.iterate, hole.bounding_ids, cdist, lo, hi)
            count += 1
    return count


def reference_witness(trellis, pair: tuple[int, int]) -> Optional[str]:
    """
    Why a reference pair of an earlier recompute is no pair of ``trellis`` (BCG-4).

    The definition (``Pseudoneighbor_Algorithm.pdf``) kills a pair when its open
    stable or unstable interval holds a point of X = every iterate of a
    crossing. Read off the snapshot's orderings and the iterate table alone (no
    scaling, no tolerance): a crossing ``w`` strictly between the two on either
    branch, or strictly between their registered images ``M^n(x0)``,
    ``M^n(x1)`` on either branch (the image of an interval is an interval).
    The map preserves order along a branch, so a witness must hold in EVERY
    frame where the table registers ``M^m(w)`` beside the pair's images: one
    whose own chain falls outside the pair there is an order flip made by
    cdist noise, not a point of X between them, and is rejected.

    Args:
        trellis: The later recompute's per-fixed-point trellis.
        pair: The two crossing ids.

    Returns:
        The witness, described, or None when the table names none.
    """
    frames = {0: pair}
    for step in (1, -1):
        x, y, n = *pair, 0
        while (x := trellis.iterate(x, step)) is not None and (
            y := trellis.iterate(y, step)
        ) is not None:
            n += step
            frames[n] = (x, y)

    def inside(w: int, x: int, y: int, stability: str) -> Optional[bool]:
        """Whether w lies strictly between x and y on their branch; None off it."""
        branch = trellis.branch_containing(x, stability)
        if branch is None or not {y, w} <= set(branch.intersection_ids):
            return None
        i, j, k = (branch.intersection_ids.index(z) for z in (x, y, w))
        return min(i, j) < k < max(i, j)

    def consistent(w: int, n: int, stability: str) -> bool:
        """No other frame orders ``M^(m-n)(w)`` outside the pair's images."""
        for m, (x, y) in frames.items():
            image = trellis.iterate(w, m - n) if m != n else None
            if image is not None and inside(image, x, y, stability) is False:
                return False
        return True

    for n, (x, y) in sorted(frames.items(), key=lambda item: abs(item[0])):
        for stability in ("stable", "unstable"):
            branch = trellis.branch_containing(x, stability)
            if branch is None or y not in branch.intersection_ids:
                continue
            ids = branch.intersection_ids
            i, j = sorted((ids.index(x), ids.index(y)))
            for w in ids[i + 1:j]:
                if consistent(w, n, stability):
                    return f"{w} between {x} and {y} on their {stability} branch (frame {n})"
    return None


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
        result = by_branch_side[(_stable_branch_at(case, bridge_id, intersection_id), row)]
        refs.append(result.ref(result.element_of_intersection[intersection_id]))
    return tuple(refs)


def _evidence(case: "Case"):
    """Every bridge's class: the EVIDENCE table the alphabet is drawn from (D2)."""
    return case.session.minimal_trellis().classes


def check_every_bridge_in_one_class(case: "Case") -> int:
    """Every identified bridge of the trellis lands in exactly one class of the evidence
    table, or is an orphan loop (both ends in one element) that the table leaves
    unclassed (D2)."""
    trellis = _full(case)
    table = _evidence(case)
    expected = sorted(bridge.id for bridge in trellis.bridges if bridge.id is not None)
    placed = [bid for entry in table for bid in entry.bridge_ids] + list(table.orphans)
    assert len(placed) == len(set(placed)), "a bridge landed in two classes"
    assert sorted(placed) == expected
    for bid in table.orphans:
        x, y = _element_ends(case, bid)
        assert x == y, bid
    return len(placed)


def check_class_orientation_anchor_outward(case: "Case") -> int:
    """A class runs anchor outward (``anchor_outward_key``), and a member's direction is
    ``+1`` exactly when its first end is the class source; a loop member sits in one
    element, and no class is a loop (D2: a loop folds into a class, never stands alone)."""
    trellis = _full(case)
    fixed_points = trellis.fixed_points
    table = _evidence(case)
    count = 0
    for entry in table:
        cls = entry.bridge_class
        assert not cls.is_loop, cls.label
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


def check_loops_fold_forward(case: "Case") -> int:
    """D2: a loop folds FORWARD. Up its image-chain parents (each pair's parent is the one
    bridge whose ``image_chain`` holds it), through loop bridges only, the first bridge
    that connects two elements is the loop's ``folded_from``, a non-loop member of the
    loop's own class. An orphan is a loop whose climb never reaches such a bridge."""
    full = _full(case)
    table = _evidence(case)
    parent: dict = {}
    for bridge in full.bridges:
        if bridge.id is None:
            continue
        for pair in image_chain(full, bridge.id) or ():
            assert parent.setdefault(frozenset(pair), bridge.id) == bridge.id, pair
    count = len(table)
    for entry in table:
        for loop in entry.loops:
            assert loop.folded_from is not None, loop.bridge_id
            assert not entry.member(loop.folded_from).is_loop, loop.bridge_id
            current, seen = parent[frozenset(loop.bridge_id)], set()
            while frozenset(current) != frozenset(loop.folded_from):
                assert frozenset(current) not in seen, (loop.bridge_id, current)
                seen.add(frozenset(current))
                x, y = _element_ends(case, current)
                assert x == y, (loop.bridge_id, current)
                current = parent[frozenset(current)]
            count += 1
    for orphan in table.orphans:
        current, seen = orphan, set()
        while current is not None and frozenset(current) not in seen:
            seen.add(frozenset(current))
            x, y = _element_ends(case, current)
            assert x == y, f"orphan {orphan} has the non-loop ancestor {current}"
            current = parent.get(frozenset(current))
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
        keys = [_stable_branch_at(case, bridge.id, i) for i in bridge.id]
        if keys[0] != keys[1]:
            continue
        first, second = rows_of_bridge(trellis, bridge.id)
        assert first == second, (bridge.id, first, second)
        count += 1
    return count


def check_anchor_bridge_class_leads_its_tangle(case: "Case") -> int:
    """A tangle's anchor bridge (leaving an anchor) is a ``+1`` member of an ACTIVE class,
    and that class comes first in its tangle's group (smallest member cdist: 0).

    Only for a tangle whose every stable branch carries a hole: a branch whose
    reference images were never computed has none (D9), and the anchor bridge can then
    map onto a loop, an inert class."""
    trellis = _full(case)
    table = _evidence(case)
    dyn = case.session.symbolic_dynamics()
    entries = list(table)
    anchors = set(_anchor_ids(case))
    tol = _tol(case)
    holed = {opening[3] for _trellis, hole in _holes(case) for opening in hole.openings or ()}
    complete = {
        index
        for index, fp in enumerate(trellis.fixed_points)
        if all(branch.key in holed for branch in trellis.branches_of(fp, "stable"))
    }
    count = 0
    for bridge in trellis.bridges:
        if bridge.id is None or bridge.id[0] not in anchors or bridge.id[1] in anchors:
            continue
        entry = table.entry_of(bridge.id)
        if entry.tangle not in complete:
            continue
        assert dyn.classes[entry.bridge_class].kind == "active", bridge.id
        assert entry.member(bridge.id).direction == +1, bridge.id
        assert entry.min_unstable_cdist == pytest.approx(0.0, abs=tol)
        group = [other for other in entries if other.tangle == entry.tangle]
        assert group[0].min_unstable_cdist == pytest.approx(0.0, abs=tol)
        count += 1
    return count


def check_tangles_in_case_order(case: "Case") -> int:
    """NEST-6: the session's tangles are the case's fixed points in order (outermost
    first), so ``entry.tangle = t`` names ``case.fixed_points[t]``: the class's elements
    and its members' unstable branches all belong to that fixed point."""
    trellis = _full(case)
    assert [id(fp) for fp in trellis.fixed_points] == [id(fp) for fp in case.fixed_points]
    count = 0
    for entry in _evidence(case):
        if entry.tangle is None:
            continue
        owner = case.fixed_points[entry.tangle]
        cls = entry.bridge_class
        assert cls.source.fixed_point is owner and cls.target.fixed_point is owner, cls.label
        for bridge_id in entry.bridge_ids:
            assert trellis.intersection(bridge_id[0]).manifold_a_key[0] is owner, bridge_id
        count += 1
    return count


def check_classes_use_every_orbit_branch(case: "Case") -> int:
    """The classes of a fixed point name elements on every stable branch of its cycle
    (all k orbit indices of a period-k orbit, both branches of an inversion point)."""
    table = _evidence(case)
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


def check_alphabet_inside_minimal_trellis(case: "Case") -> int:
    """D2: the alphabet is the minimal trellis's bridges: a class of the evidence table
    is in the session's lettered table exactly when it has a member the minimal trellis
    keeps, and one outside it is never lettered (blast children are image evidence only)."""
    minimal = case.session.minimal_trellis()
    kept = set(minimal.kept_bridge_ids)
    table = case.session.bridge_classes()
    for entry in table:
        assert any(m.bridge_id in kept for m in entry.members), entry.bridge_class.label
    alphabet = set(table.classes)
    for entry in minimal.classes:
        has_kept = any(m.bridge_id in kept for m in entry.members)
        assert has_kept == (entry.bridge_class in alphabet), entry.bridge_class.label
        assert has_kept or entry.letter is None, entry.bridge_class.label
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
    """Each iterated element lies EXACTLY inside its homotopy parent (no tolerance: every
    boundary is a registry cdist), and the parent owns the element's midpoint exactly."""
    iterated = case.session.iterated_partition()
    homotopy = iterated.homotopy
    count = 0
    for key, result in iterated.results.items():
        parent_result = homotopy.result(*key)
        for interval in result.intervals:
            parent = parent_result.element(interval.parent_element_id)
            assert parent.lo_cdist <= interval.lo_cdist <= interval.hi_cdist <= parent.hi_cdist
            mid = 0.5 * (interval.lo_cdist + interval.hi_cdist)
            assert owns_cdist(parent, mid, 0.0)
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
            key = (_stable_branch_at(case, bid, endpoint), row_of_end(full, bid, which))
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
    """Every crossing of a cut branch has exactly one iterated owner per side (by id at a
    boundary, by exact cdist elsewhere), and every stable arc midpoint of the minimal
    arrangement one exact owner per side."""
    iterated = case.session.iterated_partition()
    minimal = case.session.minimal_trellis()
    full = minimal.trellis
    count = 0
    for key, result in iterated.results.items():
        for iid in full.branch(key[0]).intersection_ids:
            cdist = float(full.intersection(iid).stable_cdist)
            owners = _owners(result.intervals, iid, cdist)
            assert owners == [result.element_of_intersection[iid]], (key, iid, owners)
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
    count = 0
    for bid in minimal.image_bridge_ids:
        branch = _stable_branch_at(case, bid, bid[0])
        if branch != _stable_branch_at(case, bid, bid[1]):
            continue
        row = row_of_end(full, bid, "first")
        if row != row_of_end(full, bid, "second"):
            continue
        result = iterated.result(branch, row)
        near_c, far_c = sorted(float(full.intersection(i).stable_cdist) for i in bid)
        under = [
            iv for iv in result.intervals
            if near_c <= iv.lo_cdist and iv.hi_cdist <= far_c
        ]
        assert under, bid
        for interval in under:
            assert interval.cut_by is not None and minimal.is_image_bridge(interval.cut_by)
            lo_c, hi_c = sorted(float(full.intersection(i).stable_cdist) for i in interval.cut_by)
            assert lo_c <= interval.lo_cdist and interval.hi_cdist <= hi_c
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
            assert span.lo_cdist <= child_interval.lo_cdist
            assert child_interval.hi_cdist <= span.hi_cdist
    return count


def check_minimal_trellis_bridges(case: "Case") -> int:
    """The minimal trellis keeps every hole bridge, maps every one without a zeroth hole
    along its registered image chain (``image_chain``, skipped pairs excepted) or lists
    it unmapped, never maps a zeroth-hole bridge (author, 2026-10-14), and keeps
    nothing else."""
    session = case.session
    minimal = session.minimal_trellis()
    full = minimal.trellis
    hole_set = {frozenset(b) for b in minimal.hole_bridge_ids}
    zeroth = set()
    for _trellis, hole in _holes(case):
        bridge = full.bridge_between(*hole.bounding_ids)
        if bridge is not None and bridge.id is not None:
            assert frozenset(bridge.id) in hole_set, hole.bounding_ids
            if hole.pair is not None:
                zeroth.add(frozenset(bridge.id))
    assert len(hole_set) == len(minimal.hole_bridge_ids)
    for bid in minimal.inert_hole_bridge_ids + minimal.unclassed_hole_bridge_ids:
        assert minimal.is_hole_bridge(bid) and bid not in minimal.image_chains
    assert {frozenset(b) for b in minimal.inert_hole_bridge_ids} == zeroth
    for bid in [*minimal.image_chains, *minimal.unmapped]:
        assert frozenset(bid) not in zeroth, bid
    unclassed = {frozenset(b) for b in minimal.unclassed_hole_bridge_ids}
    mapped = {frozenset(b) for b in [*minimal.image_chains, *minimal.unmapped]}
    assert hole_set - zeroth == unclassed | mapped
    for bid, chain in minimal.image_chains.items():
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


def _distance_to_segments(point: np.ndarray, starts: np.ndarray, ends: np.ndarray) -> float:
    """The distance from a point to the nearest of some segments."""
    along = ends - starts
    t = np.clip(((point - starts) * along).sum(1) / np.maximum((along * along).sum(1), 1e-300), 0.0, 1.0)
    return float(np.linalg.norm(starts + t[:, None] * along - point, axis=1).min())


def check_dual_face_side_is_geometry(case: "Case") -> int:
    """A face attaches to the side node of the stable edge it geometrically lies on: a
    probe just off the edge's middle segment, on the node's side, lies in the face.

    The probe is LOCAL (WP10): pushed off the middle segment by a quarter of its
    distance to every OTHER segment of the face's boundary (the edge's own arc
    included, which can fold back as a hairpin), so it cannot leave the face through
    another wall. A half-plane test against the face's representative point failed on
    slivers and non-convex faces (p3 at 16 unstable steps).
    """
    dual = case.session.dual_graph()
    trellis = dual.arrangement.trellis
    count = 0
    checkable = 0
    for region in dual.arrangement.regions:
        face = dual.face_of(region)
        checkable += len(region.stable_arcs)
        boundary = np.asarray(region.boundary_points, dtype=float)
        walls = np.stack([boundary, np.roll(boundary, -1, axis=0)], axis=1)
        for arc in region.stable_arcs:
            side = "right" if arc.reverse else "left"
            node = dual.nodes_of_edge(arc.edge_key)[side]
            poly = np.asarray(node.arc.polyline(trellis), dtype=float)
            if len(poly) < 2:
                continue
            middle = max(1, len(poly) // 2)
            below, above = poly[middle - 1], poly[middle]  # the polyline runs anchor outward
            look = (below - above) / np.linalg.norm(below - above)
            left = np.array([-look[1], look[0]])
            midpoint = 0.5 * (below + above)
            starts, ends = walls[:, 0], walls[:, 1]
            itself = ((starts == below).all(1) & (ends == above).all(1)) | (
                (starts == above).all(1) & (ends == below).all(1)
            )
            others = walls[~itself]
            clearance = _distance_to_segments(midpoint, others[:, 0], others[:, 1]) if len(others) else np.inf
            push = 0.25 * min(clearance, float(np.linalg.norm(above - below)))
            probe = midpoint + (1.0 if node.side_of(face) == "left" else -1.0) * push * left
            assert region.contains(probe, tol=0.0), (region.corners, arc.lo_id, arc.hi_id)
            assert node.faces[side] is face
            count += 1
    assert count == checkable, "every (region, stable edge) pair is checkable"
    return count


def misplaced_unified_crossings(dual) -> list:
    """``[(stable branch key, crossing id)]`` for every crossing of a pip's fixed point
    whose place disagrees with the unified span. ``S \\ f(S)`` is read off the iterate
    table, independently of the dual graph's cdist arithmetic: a crossing lies in it
    exactly when it has no preimage IN THE TRELLIS (its preimage would sit beyond the
    end of the trellis's stable branch, which is cut at the pip orbit, D10), so it must
    be unified exactly then."""
    trellis = dual.trellis
    tol = float(trellis.registry.cdist_tol)
    owners = {key[0] for key in dual.fundamental_segments}
    misplaced = []
    for branch in trellis.stable_branches:
        if branch.key[0] not in owners:
            continue
        span = dual.fundamental_segments.get(branch.key)
        for crossing_id in branch.intersection_ids:
            cdist = float(trellis.intersection(crossing_id).stable_cdist)
            if cdist == 0.0:  # the anchor is its own preimage
                continue
            unified = span is not None and span[0] + tol < cdist <= span[1] + tol
            if unified != (trellis.iterate(crossing_id, -1) not in trellis):
                misplaced.append((branch.key, crossing_id))
    return misplaced


def check_dual_unified_is_where_f_of_s_stops(case: "Case") -> int:
    """The precondition of the walk (PK-5, D10): on every stable branch of a pip's
    fixed point the unified span is exactly ``S \\ f(S)``, where the image of a bridge
    may cross the stable manifold away from the images of its ends: a crossing is
    unified iff it has no preimage in the trellis."""
    dual = case.session.dual_graph()
    assert misplaced_unified_crossings(dual) == []
    return sum(len(branch) for branch in dual.trellis.stable_branches)


def check_stable_nodes_never_own_a_singleton(case: "Case") -> int:
    """A stable node's element on each faced side is never a singleton ``[x, x]`` (a
    singleton owns no stable edge), spans the node's edge, and owns its midpoint exactly."""
    dual = case.session.dual_graph()
    count = 0
    for node in dual.stable_nodes.values():
        for side in node.sides:
            interval = node.intervals[side]
            assert interval.lo_id != interval.hi_id, (node.key, side, interval.lo_id)
            assert interval.lo_cdist <= node.lo_cdist and node.hi_cdist <= interval.hi_cdist
            assert owns_cdist(interval, node.mid_cdist, 0.0), (node.key, side)
            count += 1
    return count


# --------------------------------------------------------------------------- #
# Markers and enclosure: the non-crossing rule (D13)
# --------------------------------------------------------------------------- #
def _unstable_edge_keys(dual) -> dict:
    """``{frozenset(bridge id): edge key}`` of every unstable arc of the dual graph's arrangement."""
    return {
        frozenset(arc.bridge_id): arc.edge_key
        for face in dual.arrangement.faces
        for arc in face.arcs
        if arc.kind == "unstable"
    }


def _closed_curve(case: "Case", dual, pairs: list, edges: dict) -> Optional[set]:
    """The edges of the closed curve made of the consecutive unstable arcs ``pairs`` and
    the stable segment between their two outer ends, or None when those ends lie on two
    stable branches. An anchor lies on every branch of its point (D4), so an anchor end
    closes along the other end's branch: any stable path between the ends makes a
    Jordan curve, whatever row the arc's end sits against there (C7)."""
    ends = (pairs[0][0], pairs[-1][1])
    branches = {case.registry[end].manifold_b_key for end in ends if not _is_anchor(case.registry[end])}
    if len(branches) != 1:
        return None
    lo, hi = sorted(float(case.registry[end].stable_cdist) for end in ends)
    stable = {
        node.edge_key
        for node in dual.stable_nodes_on(branches.pop())
        if lo <= node.lo_cdist and node.hi_cdist <= hi
    }
    return stable | {edges[frozenset(pair)] for pair in pairs}


def check_lobe_enclosure_is_geometry(case: "Case") -> int:
    """D13 decides enclosure combinatorially, so it must agree with the curves: the lobe
    of every kept bridge whose ends lie on one stable branch (the bridge and the stable
    segment between its ends) holds a periodic point by the dual graph's flood exactly
    when point-in-polygon on the lobe's polylines says so. A point within 1e-6 of the
    lobe is skipped."""
    dual = case.session.dual_graph()
    trellis = dual.trellis
    edges = _unstable_edge_keys(dual)
    arcs = {
        arc.edge_key: arc.reversed() if arc.reverse else arc
        for face in dual.arrangement.faces
        for arc in face.arcs
    }
    count = 0
    for bridge_id in dual.minimal.kept_bridge_ids:
        curve = _closed_curve(case, dual, [bridge_id], edges)
        if curve is None:
            continue
        (held,) = dual.enclosed(curve)
        bridge = arcs[edges[frozenset(bridge_id)]].polyline(trellis)
        stable = sorted(
            (arcs[key] for key in curve if key[0] == "stable"),
            key=lambda arc: float(case.registry[arc.lo_id].stable_cdist),
        )
        segment = np.vstack([arc.polyline(trellis) for arc in stable])
        if np.linalg.norm(segment[0] - bridge[-1]) > np.linalg.norm(segment[-1] - bridge[-1]):
            segment = segment[::-1]
        lobe = np.vstack([bridge, segment])
        for position, fixed_point in enumerate(trellis.fixed_points):
            for orbit_index in range(fixed_point.period):
                point = np.ravel(fixed_point.coordinates[orbit_index])
                if np.min(np.linalg.norm(lobe - point, axis=1)) < 1e-6:
                    continue
                inside = held[("orbit", position, orbit_index)]
                assert inside == point_in_polygon(point, lobe), (bridge_id, position, orbit_index)
                count += 1
    return count


def check_registered_images_keep_markers(case: "Case") -> int:
    """The physics of D13: the map is a homeomorphism, so where a kept bridge's whole
    registered image chain is kept, the closed curve of the chain and the stable segment
    between its ends holds the image of every marker the bridge's lobe holds, and of no
    marker the lobe leaves out (wherever both curves decide)."""
    dual = case.session.dual_graph()
    edges = _unstable_edge_keys(dual)
    count = 0
    for bridge_id in dual.minimal.kept_bridge_ids:
        chain = image_chain(dual.trellis, bridge_id)
        if not chain or any(frozenset(pair) not in edges for pair in chain):
            continue
        lobe = _closed_curve(case, dual, [bridge_id], edges)
        image = _closed_curve(case, dual, chain, edges)
        if lobe is None or image is None:
            continue
        (held,) = dual.enclosed(lobe)
        (image_held,) = dual.enclosed(image)
        for marker, inside in held.items():
            image_marker = dual.image_marker(marker)
            if inside is None or image_held.get(image_marker) is None:
                continue
            assert image_held[image_marker] == inside, (bridge_id, marker, inside)
            count += 1
    return count


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


def check_cross_side_only_within_a_branch(case: "Case") -> int:
    """D4: an itinerary pair on ONE stable branch joins two elements of one side, and
    only such a pair may be flagged cross-side; a pair across two branches never is (at
    an inversion point L of 0.0 is R of 0.1)."""
    dyn = case.session.symbolic_dynamics()
    count = 0
    for cd in _resolved(dyn):
        for x, y in zip(cd.itinerary[::2], cd.itinerary[1::2]):
            if x.branch_key == y.branch_key:
                assert x.side == y.side, (cd.letter, x.label, y.label)
            count += 1
        for symbol in cd.symbols:
            x, y = symbol.occurrence
            assert not (symbol.cross_side and x.branch_key != y.branch_key), (
                cd.letter, x.label, y.label)
    return count


def check_every_active_class_walks(case: "Case") -> int:
    """D3: the shortest dual-graph walk IS the word, so it exists: every active class
    whose landings are not singletons has a ``resolved`` walk (a singleton
    takes ``f(bridge)``, the trellis itinerary)."""
    dyn = case.session.symbolic_dynamics()
    count = 0
    for cd in dyn.classes.values():
        if cd.kind != "active" or any(l is not None and l.singleton for l in cd.landings):
            continue
        assert cd.search is not None and cd.search.status == "resolved", (
            cd.letter, cd.status, cd.unresolved_reason)
        count += 1
    return count


def check_walks_agree_with_registered_images(case: "Case") -> int:
    """D3 + D13: the walk is the word and a member's registered image its cross-check, so
    every registered image that gives a verdict agrees with the walk at word level and no
    class is contradicted by its evidence (C11: incomplete ``"stable end"`` evidence gives
    no verdict)."""
    dyn = case.session.symbolic_dynamics()
    count = 0
    for cd in dyn.classes.values():
        for evidence in cd.evidence:
            if evidence.word_agrees is None:
                continue
            assert evidence.word_agrees, (cd.letter, evidence.bridge_id, evidence.word, cd.word)
            count += 1
        assert cd.verified is not False, cd.letter
    return count


def check_member_lobe_classes_have_one_walk(case: "Case") -> int:
    """D13: a walked class with a member lobe -- a non-loop member the minimal trellis
    keeps whose two ends lie on one stable branch -- has exactly one walk (a ``resolved``
    search): the non-crossing rule picks its side, so only a class without one
    may stay ambiguous, and the image of the member is a walk holding the right
    markers, so one is found (no ``no_admissible_walk``, no ``truncated`` reason)."""
    dyn = case.session.symbolic_dynamics()
    kept = set(case.session.minimal_trellis().kept_bridge_ids)
    count = 0
    for cd in dyn.classes.values():
        lobes = [
            member
            for member in cd.entry.members
            if not member.is_loop
            and member.bridge_id in kept
            and len({_stable_branch_at(case, member.bridge_id, end) for end in member.bridge_id}) == 1
        ]
        if cd.search is None or not lobes:
            continue
        assert cd.search.status == "resolved", (
            cd.letter, cd.search.status, [walk.itinerary for walk in cd.search.walks])
        count += 1
    return count


def _span_meets(node, landing) -> bool:
    """A stable node's edge overlaps a landing span by more than the span's error (nothing
    at a table end, ``SCALING_RTOL`` relative at a scaled one)."""
    lo, hi = landing.span
    error = max(
        [SCALING_RTOL * abs(c) for c, exact in zip(landing.span, landing.from_table) if not exact]
        + [0.0]
    )
    return min(hi, node.hi_cdist) - max(lo, node.lo_cdist) > error


def check_walk_starts_on_the_landing(case: "Case") -> int:
    """D3, D13: the walk runs from the image of one element to the image of the other, so
    it leaves from a face beside a stable edge the source's landing span meets and
    arrives at one beside an edge the target's span meets (when some edge of the element
    meets it)."""
    dyn = case.session.symbolic_dynamics()
    dual = case.session.dual_graph()
    count = 0
    for cd in dyn.classes.values():
        if cd.source != "walk":
            continue
        walk = cd.search.walk
        for landing, face in ((cd.landings[0], walk.faces[0]), (cd.landings[1], walk.faces[-1])):
            target = landing.target
            meeting = [node for node in dual.stable_nodes_of(target) if _span_meets(node, landing)]
            if meeting:
                assert face in [node.face_on(target.side) for node in meeting], (cd.letter, target.label)
            count += 1
    return count


def check_alphabet_closed_under_words(case: "Case") -> int:
    """D3, D11: every symbol of every word names a class of the table. The seed classes
    are the session alphabet; every other class was added by the closure and is named
    in the word of the class that named it. A class without an itinerary makes the
    dynamics unreliable. (That a closure class is spanned by a registered bridge exactly
    when it is not flagged "no computed bridge yet" holds by construction: its entry IS
    the evidence entry when there is one, WP11.)"""
    dyn = case.session.symbolic_dynamics()
    seed = case.session.bridge_classes()
    count = 0
    for letter, word in dyn.rules.items():
        for symbol in word:
            assert symbol.bridge_class in dyn.table, (letter, symbol.text)
            count += 1
    assert list(dyn.classes) == dyn.table.classes
    for cls, cd in dyn.classes.items():
        assert (cd.named_by is None) == (cls in seed), cd.letter
        if cd.named_by is not None:
            namer = dyn.classes[cd.named_by]
            assert cls in {symbol.bridge_class for symbol in namer.symbols}, cd.letter
        if cd.itinerary is None:
            assert not dyn.is_reliable, cd.letter
    # D14: a class connects two tangles only where a computed bridge does.
    connecting = {
        frozenset((id(cls.source.fixed_point), id(cls.target.fixed_point)))
        for cls in dyn.classes
        if cls.source.fixed_point is not cls.target.fixed_point
    }
    spanned = {
        frozenset(id(_stable_branch_at(case, bridge.id, end)[0]) for end in bridge.id)
        for bridge in _staged(case).bridges
        if bridge.id is not None
    }
    for pair in connecting:
        assert pair in spanned, "a class connects two tangles no computed bridge connects"
    return count + len(dyn.classes) + len(connecting)


def _arc_branches(case: "Case", dyn) -> dict:
    """The unstable branches each class's arcs lie on: its members', or, for a class no
    computed bridge spans yet, one map step past its namer's (D11; a namer comes first)."""
    trellis = case.session.dual_graph().trellis
    branches: dict = {}
    for cls, cd in dyn.classes.items():
        if cd.entry.members:
            branches[cls] = {trellis.unstable_key_of(member.bridge_id) for member in cd.entry.members}
        else:
            branches[cls] = {key[0].advance_key(key, 1) for key in branches[cd.named_by]}
    return branches


def _chords_interleave(face, chord, other) -> bool:
    """Four distinct chord ends on one ring of the face that alternate around it."""
    ends = (chord[:2], chord[2:], other[:2], other[2:])
    for region in face.faces:
        at = {(arc.edge_key, face_side(arc)): i for i, arc in enumerate(region.arcs)}
        if all(end in at for end in ends):
            a, b, c, d = sorted(at[end] for end in ends[:2]) + [at[end] for end in ends[2:]]
            return len({a, b, c, d}) == 4 and (a < c < b) != (a < d < b)
    return False


def check_walks_avoid_disjoint_images(case: "Case") -> int:
    """D13, D14: the map is a homeomorphism and two unstable arcs never cross, so the
    images of arcs on disjoint sets of unstable branches are disjoint. Two classes whose
    arcs lie on disjoint branch sets (``_arc_branches``) have decided walks (``resolved``,
    of multiplicity 1: a walk standing for several face paths names
    none of them) whose chords never interleave on a ring of a face both visit. Counts
    the class pairs judged (a pair sharing no face cannot cross: a deep run may judge
    only such pairs)."""
    dyn = case.session.symbolic_dynamics()
    dual = case.session.dual_graph()
    branches = _arc_branches(case, dyn)
    decided = [
        (cls, cd, walk_chords(cd.search.walk))
        for cls, cd in dyn.classes.items()
        if cd.search is not None
        and cd.search.status == "resolved"
        and cd.search.walk.multiplicity == 1
    ]
    count = 0
    for i, (x, cx, chords_x) in enumerate(decided):
        for y, cy, chords_y in decided[i + 1:]:
            if branches[x] & branches[y]:
                continue
            count += 1
            for face, chord in chords_x:
                for other_face, other in chords_y:
                    if other_face is not face:
                        continue
                    crossing = _chords_interleave(face, chord, other)
                    assert dual.chords_cross(face, chord, other) == crossing
                    assert not crossing, (cx.letter, cy.letter, face.index)
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


def p1_itinerary_disagreements(
    session: "TangleSession",
) -> tuple[int, tuple["BridgeId", ...], tuple[str, ...]]:
    """
    Theorem T read off the registered images (A.8d), lookup only.

    Every non-loop member of every class whose image chain is registered is
    read at the ``P1`` level: ``trellis_itinerary`` on the topology's trellis and
    the iterated partition, then ``translate_itinerary``, as ``(class,
    direction, occurrence)`` tokens, freely reduced (a token next to its inverse
    on the same occurrence cancels). A chain with two consecutive crossings on
    one stable branch with one crossing sign passed the end of the computed
    stable branch (the "stable end" of C11, recomputed here from the chain) and
    a chain the iterated partition does not cover is unreadable: both are
    incomplete evidence and skipped.

    Args:
        session: A session with its symbolic dynamics built.

    Returns:
        ``(tested, skipped, disagreements)``: how many classes had two or more
        readable members, the skipped members, and one text per class whose
        members' reduced itineraries differ (its letter and label, then each
        distinct itinerary with the members reading it).
    """
    dyn = session.symbolic_dynamics()
    trellis = session.topology().trellis
    iterated = session.iterated_partition()
    letters = {cls: cd.letter for cls, cd in dyn.classes.items()}

    def stable_end(chain) -> bool:
        crossings = [chain[0][0], *(second for _first, second in chain)]
        marks = [
            (trellis.intersection(c).manifold_b_key, trellis.intersection(c).crossing_sign)
            for c in crossings
        ]
        return any(a == b for a, b in zip(marks, marks[1:]))

    def reduced(symbols) -> tuple:
        word: list = []
        for symbol in symbols:
            token = (symbol.bridge_class, symbol.direction, symbol.occurrence)
            if word and word[-1] == (token[0], -token[1], token[2]):
                word.pop()
            else:
                word.append(token)
        return tuple(word)

    def text(word) -> str:
        return " ".join(
            f"{letters.get(cls, cls.label)}{'^-1' if direction < 0 else ''}"
            f"({dyn.naming.name(lo).text}, {dyn.naming.name(hi).text})"
            for cls, direction, (lo, hi) in word
        )

    tested, skipped, disagreements = 0, [], []
    for cls, cd in dyn.classes.items():
        readings: dict[tuple, list] = {}
        for member in cd.entry.members:
            chain = None if member.is_loop else image_chain(trellis, member.bridge_id)
            if not chain:
                continue
            if stable_end(chain):
                skipped.append(member.bridge_id)
                continue
            try:
                itinerary = trellis_itinerary(trellis, iterated, chain, member.direction)
            except ValueError:
                skipped.append(member.bridge_id)
                continue
            symbols, _loops = translate_itinerary(
                itinerary, dyn.naming, dyn.table, letters, dyn.fixed_points
            )
            readings.setdefault(reduced(symbols), []).append(member.bridge_id)
        if sum(map(len, readings.values())) < 2:
            continue
        tested += 1
        if len(readings) > 1:
            disagreements.append(
                f"class {cd.letter} {cls.label}: "
                + " | ".join(f"{members}: {text(word)}" for word, members in readings.items())
            )
    return tested, tuple(skipped), tuple(disagreements)


def _lobe_image(trellis, workbench, chain: list, base_key) -> np.ndarray:
    """The polygon of ``f(H)``: the image chain's arcs, closed along the base."""
    def arc(a: int, b: int) -> np.ndarray:
        bridge = trellis.bridge_between(a, b)
        assert bridge is not None, (a, b)
        ends = [np.ravel(trellis.intersection(c).coords).astype(float) for c in (a, b)]
        lo, hi = sorted(trellis.intersection(c).unstable_cdist for c in (a, b))
        return np.vstack([ends[0], arc_polyline(bridge, lo, hi), ends[1]])

    first, last = chain[0][0], chain[-1][1]
    lo, hi = sorted(trellis.intersection(c).stable_cdist for c in (first, last))
    base = arc_polyline(workbench.manifolds[base_key], lo, hi)
    if trellis.intersection(last).stable_cdist > trellis.intersection(first).stable_cdist:
        base = base[::-1]
    return np.vstack([*(arc(a, b) for a, b in chain), base])


def _arclength_midpoint(points: np.ndarray) -> np.ndarray:
    """The point halfway along a polyline by arc length."""
    lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    along = np.concatenate([[0.0], np.cumsum(lengths)])
    i = min(int(np.searchsorted(along, along[-1] / 2, side="right")) - 1, len(lengths) - 1)
    t = (along[-1] / 2 - along[i]) / lengths[i] if lengths[i] else 0.0
    return points[i] + t * (points[i + 1] - points[i])


def base_chord_misreadings(session: "TangleSession") -> tuple[int, tuple[str, ...]]:
    """
    The base-branch chords of every lobe image, read by point-in-polygon (B7b).

    A lobe image ``f(H)`` holds no anchor, so walking its BASE branch from the
    anchor enters and leaves it in turn, the base counting as inside (Jordan).
    For every hole with a stretch whose consecutive image chain crosses its base
    branch at two or more interior crossings, the branch is cut at the chain's
    crossings on it (and the anchor and the branch tail at its ends); each gap,
    closed at its two crossings' own coordinates, is read at its arc-length
    midpoint against :func:`_lobe_image`, ``tol=0``. Polylines numerics already
    computed only; the map is never called. Checked:

    * (a) the base is a single gap;
    * (b) the readings alternate, the first gap outside unless it is the base
      (an anchor-cornered lobe, the inversion's);
    * (c) every inside gap but the base is an empty stretch of the iterated cut;
    * (d) every base or chord the cut declares for THIS lobe on its base branch
      (``_empty_stretches`` over this hole bridge alone) is one inside gap.

    Args:
        session: A session with its minimal trellis and homotopy partition.

    Returns:
        ``(lobes checked, one text per misreading)``.
    """
    minimal = session.minimal_trellis()
    homotopy = session.homotopy_partition()
    trellis, workbench = minimal.trellis, session.workbench
    at = trellis.intersection
    empty = IteratedHomotopyPartition._empty_stretches(minimal, homotopy)
    stretched = _hole_stretches(homotopy)
    checked, misreadings = 0, []
    for hole_bridge in minimal.hole_bridge_ids:
        if frozenset(hole_bridge) not in stretched:
            continue
        chain = image_chain(trellis, hole_bridge)
        if not chain or any(pair[1] != after[0] for pair, after in zip(chain, chain[1:])):
            continue
        crossings = [chain[0][0], *(second for _first, second in chain)]
        (base_key,) = {
            at(c).manifold_b_key for c in (crossings[0], crossings[-1]) if not at(c).is_anchor
        }
        if sum(at(c).manifold_b_key == base_key for c in crossings[1:-1]) < 2:
            continue
        checked += 1
        manifold = workbench.manifolds[base_key]
        ring = _lobe_image(trellis, workbench, chain, base_key)
        on_branch = sorted(
            {c for c in crossings if at(c).manifold_b_key == base_key or at(c).is_anchor},
            key=lambda c: at(c).stable_cdist,
        )
        ends = [
            ("anchor", 0.0, np.array([manifold.root.x, manifold.root.y], float)),
            *((c, at(c).stable_cdist, np.ravel(at(c).coords).astype(float)) for c in on_branch),
            ("tail", manifold.tail.cdist, np.array([manifold.tail.x, manifold.tail.y], float)),
        ]
        base = frozenset((crossings[0], crossings[-1]))
        gaps = {}  # (lo, hi) -> inside, anchor outward
        for (lo, lo_cdist, lo_point), (hi, hi_cdist, hi_point) in zip(ends, ends[1:]):
            if hi_cdist <= lo_cdist:
                continue  # the anchor is a chain crossing: no gap before it
            nodes = arc_polyline(manifold, lo_cdist, hi_cdist)
            gap = np.vstack([lo_point, nodes.reshape(-1, 2), hi_point])
            gaps[lo, hi] = frozenset((lo, hi)) == base or point_in_polygon(
                _arclength_midpoint(gap), ring, tol=0.0
            )
        readings = " ".join(
            f"({lo},{hi}){'B' if frozenset((lo, hi)) == base else 'in' if inside else 'out'}"
            for (lo, hi), inside in gaps.items()
        )

        def misread(clause: str) -> None:
            misreadings.append(f"hole bridge {hole_bridge} on {base_key[1:]}: {clause}; {readings}")

        if sum(frozenset(gap) == base for gap in gaps) != 1:
            misread("(a) the base is not one gap")
        inside = list(gaps.values())
        first_is_base = frozenset(next(iter(gaps))) == base
        if inside[0] != first_is_base or any(a == b for a, b in zip(inside, inside[1:])):
            misread("(b) the gaps do not alternate")
        for (lo, hi), is_inside in gaps.items():
            if is_inside and frozenset((lo, hi)) != base and not any(
                other == hi for other, _kind in empty.get(lo, ())
            ):
                misread(f"(c) the inside gap ({lo},{hi}) is not empty")
        own = IteratedHomotopyPartition._empty_stretches(
            SimpleNamespace(trellis=trellis, hole_bridge_ids=[hole_bridge]), homotopy
        )
        inside_gaps = {frozenset(gap) for gap, is_inside in gaps.items() if is_inside}
        for lo, stretches in own.items():
            for hi, kind in stretches:
                on_base = {at(c).manifold_b_key for c in (lo, hi) if not at(c).is_anchor}
                if (
                    lo < hi
                    and kind in ("base", "chord")
                    and on_base == {base_key}
                    and frozenset((lo, hi)) not in inside_gaps
                ):
                    misread(f"(d) its {kind} ({lo},{hi}) is not one inside gap")
    return checked, tuple(misreadings)


def check_members_share_the_P1_itinerary(case: "Case") -> int:
    """Theorem T (A.8d): every member of a class whose image is registered has the same
    ``P1``-level freely reduced itinerary; stable-end chains are skipped."""
    tested, _skipped, disagreements = p1_itinerary_disagreements(case.session)
    assert not disagreements, "\n".join(disagreements)
    return tested


def check_refined_matrix_is_an_in_splitting(case: "Case") -> int:
    """``M_1`` is the in-splitting of ``M`` (A.8d): sibling rows are equal, every row's
    children columns sum to its parent's entry of ``M``, and the nodes are the
    children of the split active classes and the letters of the unsplit ones."""
    dyn = case.session.symbolic_dynamics()
    names, class_matrix = dyn.transition_matrix(refined=False)
    refined_names, matrix = dyn.transition_matrix(refined=True)
    parent: dict[str, str] = {}
    for cls, cd in dyn.classes.items():
        if cd.kind == "active":
            children = dyn.refined.get(cls)
            parent.update({child.name: cd.letter for child in children or ()})
            if not children:
                parent[cd.letter] = cd.letter
    assert sorted(refined_names) == sorted(parent), (refined_names, sorted(parent))
    assert sorted(set(parent.values())) == sorted(names), (names, sorted(set(parent.values())))
    rows = {name: position for position, name in enumerate(refined_names)}
    columns = {letter: position for position, letter in enumerate(names)}
    children_of: dict[str, list[int]] = defaultdict(list)
    for name, letter in parent.items():
        children_of[letter].append(rows[name])
    for letter, siblings in children_of.items():
        for row in siblings[1:]:
            assert matrix[row].tolist() == matrix[siblings[0]].tolist(), (
                f"sibling rows of {letter} differ: {refined_names[siblings[0]]} "
                f"{matrix[siblings[0]].tolist()} vs {refined_names[row]} {matrix[row].tolist()}"
            )
    for name, row in rows.items():
        for letter, siblings in children_of.items():
            summed = int(matrix[row, siblings].sum())
            expected = int(class_matrix[columns[parent[name]], columns[letter]])
            assert summed == expected, (
                f"row {name}: the columns of {letter}'s children sum to {summed}, "
                f"M[{parent[name]}][{letter}] = {expected}"
            )
    return len(refined_names)


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


def check_inert_is_transient(case: "Case") -> int:
    """Q3 (I-a): a class of the closed table is inert iff, in the graph of its
    class-level words (loops dropped), neither a cycle nor an unresolved class is
    reachable from it. Reachability by its own search, never by ``inert_classes``."""
    dyn = case.session.symbolic_dynamics()
    edges = {cls: {symbol.bridge_class for symbol in cd.symbols} for cls, cd in dyn.classes.items()}
    unresolved = {cls for cls, cd in dyn.classes.items() if cd.itinerary is None}

    def reachable(start) -> set:
        seen, stack = set(), list(edges[start])
        while stack:
            node = stack.pop()
            if node not in seen:
                seen.add(node)
                stack.extend(edges[node])
        return seen

    reach = {cls: reachable(cls) for cls in edges}
    for cls, cd in dyn.classes.items():
        persistent = any(
            node in unresolved or node in reach[node] for node in reach[cls] | {cls}
        )
        assert (cd.kind == "inert") is not persistent, (cd.letter, cd.kind, cls.label)
        assert cd.entry.inert is (cd.kind == "inert"), (cd.letter, cls.label)
    return len(dyn.classes)


def check_zeroth_hole_classes_have_empty_words(case: "Case") -> int:
    """The author's rule (2026-10-14): the class of a bridge carrying a zeroth hole
    (punched on a pseudoneighbor pair) walks, its word is empty -- the lobe's one
    marker maps to no marker -- and so it is inert."""
    dyn = case.session.symbolic_dynamics()
    table = case.session.minimal_trellis().classes
    count = 0
    for trellis, hole in _holes(case):
        if hole.pair is None:
            continue
        cls = table.entry_of(trellis.bridge_between(*hole.bounding_ids).id).bridge_class
        cd = dyn.classes[cls]
        assert cd.itinerary is not None, (cd.letter, cls.label, cd.unresolved_reason)
        assert not cd.symbols, (cd.letter, cls.label, cd.word)
        assert cd.kind == "inert", (cd.letter, cls.label)
        count += 1
    return count


def check_matrix_is_token_counts(case: "Case") -> int:
    """Entry ``(i, j)`` counts the tokens of symbol ``j`` (inert tokens skipped) in the
    refined word of active symbol ``i``; every token of an active word is a node."""
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
            if symbol.kind == "inert":
                continue
            expected[index[name], index[symbol.base]] += 1
    assert matrix.tolist() == expected.tolist()
    return int(matrix.sum()) + len(names)



def check_k10_is_the_full_shift(case: "Case") -> int:
    """On k=10, and only there, the dynamics is the full 2-shift (author Q11, the one
    textbook pin): one active class whose word is two occurrences of itself in opposite
    directions, class matrix ``[[2]]``, refined matrix its in-splitting
    ``[[1, 1], [1, 1]]`` (A.8d: two occurrences, two children). Letter-free; every
    other case checks nothing (``cases.NOT_APPLICABLE``)."""
    if case.name != "k10":
        return 0
    dyn = case.session.symbolic_dynamics()
    active = [cls for cls, cd in dyn.classes.items() if cd.kind == "active"]
    assert len(active) == 1, active
    tokens = [symbol for symbol in dyn.classes[active[0]].symbols if symbol.kind == "active"]
    assert [symbol.bridge_class for symbol in tokens] == active * 2
    assert sorted(symbol.direction for symbol in tokens) == [-1, 1]
    assert len({symbol.occurrence for symbol in tokens}) == 2
    assert dyn.transition_matrix(refined=False)[1].tolist() == [[2]]
    assert dyn.transition_matrix()[1].tolist() == [[1, 1], [1, 1]]
    return len(tokens)

#: The read-only law checks by layer: ONE law-tier test per (layer, case) runs
#: every check of its layer (``helpers.law_tier.run_layer``).
LAYERS: dict[str, list[LawCheck]] = {
    "case_sanity": [
        check_case_builds,
        check_case_area_cutoff,
        check_k_value_matches_inversion,
        check_fixed_point_minimal_period,
        check_map_is_area_preserving,
        check_every_branch_is_built,
        check_every_branch_crosses,
        check_case_heteroclinic,
    ],
    "manifolds": [
        check_manifold_cdist_monotone,
        check_manifold_no_spikes,
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
        check_iterate_table_partial_bijection,
        check_iterate_links_pass_the_map_gate,
        check_iterate_table_complete,
    ],
    "anchors": [
        check_one_anchor_per_periodic_point,
        check_anchor_rays_are_its_branches,
        check_anchor_sign_matches_eigendirections,
    ],
    "bridges": [
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
        check_anchor_opening_faces_its_hole,
        check_propagated_holes_land_on_predicted_branch,
        check_propagation_terminates,
        check_no_forward_hole,
        check_trellis_crossings_belong_to_its_selection,
        check_trellis_is_cut_at_the_pip_orbit,
        check_reference_pairs_valid,
        check_chosen_pip_is_strong,
        check_reference_pairs_clear_of_X,
        check_propagated_hole_inside_container,
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
        check_loops_fold_forward,
        check_row_of_end_is_geometry,
        check_bridge_rows_consistent,
        check_anchor_bridge_class_leads_its_tangle,
        check_tangles_in_case_order,
        check_classes_use_every_orbit_branch,
        check_only_active_classes_lettered,
        check_alphabet_inside_minimal_trellis,
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
        check_dual_unified_is_where_f_of_s_stops,
        check_stable_nodes_never_own_a_singleton,
        check_dual_bipartite_degree,
        check_lobe_enclosure_is_geometry,
        check_registered_images_keep_markers,
    ],
    "symbolic": [
        check_itineraries_even,
        check_cross_side_only_within_a_branch,
        check_every_active_class_walks,
        check_walks_agree_with_registered_images,
        check_member_lobe_classes_have_one_walk,
        check_walks_avoid_disjoint_images,
        check_walk_starts_on_the_landing,
        check_alphabet_closed_under_words,
        check_landings_contained,
        check_members_share_the_P1_itinerary,
        check_refined_matrix_is_an_in_splitting,
        check_member_matching_consistent,
        check_inert_classes_outside_transitions,
        check_inert_is_transient,
        check_zeroth_hole_classes_have_empty_words,
        check_matrix_is_token_counts,
        check_k10_is_the_full_shift,
    ],
}

#: The laws that run on ``cases.LAW_CASES`` only, never on a deep build: the
#: scaling-only pair oracle (scaling alone is wrong at depth, A.5), the k=10
#: textbook pin (plan B0, section 6) and Theorem T, which has no class with two
#: chained members to test on p3_15, p3_6_blasts and inversion_period2_10 (B7);
#: its depth is the blast series' (``test_members_share_the_P1_itinerary_under_blasting``).
LAW_CASES_ONLY: tuple[str, ...] = (
    "reference_pairs_clear_of_X",
    "k10_is_the_full_shift",
    "members_share_the_P1_itinerary",
)

#: The layers whose laws mutate their session: run on ONE fresh build per case,
#: never the shared one (the checks run in order on that build).
MUTATING_LAYERS: dict[str, list[LawCheck]] = {
    "recompute": [
        check_recompute_keeps_the_topology,
        check_recompute_preserves_ids,
        check_recompute_is_idempotent,
    ],
}

#: The laws that mutate their session.
MUTATING: tuple[LawCheck, ...] = tuple(MUTATING_LAYERS["recompute"])


