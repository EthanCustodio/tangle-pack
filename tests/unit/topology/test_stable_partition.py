"""Hole punching, backward propagation, and the stable-manifold partition.

The synthetic tests fabricate holes directly on a hand-built trellis to pin the
interval logic (open/closed ends, singletons, the both-sides pinch) without
any manifold numerics; the Hénon tests run the real punch → propagate →
partition pipeline on a computed tangle. The hole and partition laws (coverage,
classification, inward pairs, the openings law, propagation landing and
termination, no hole at iterate >= k_value) run on every law case in
``tests/invariants/test_law_partition.py``.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

from cases import build_k28, build_period3
from helpers.logs import assert_logged
from helpers.fakes import bare_fixed_point
from tanglepack.numerics.FixedPoint import FixedPoint
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry
from tanglepack.topology.StablePartition import partition_stable_manifold
from tanglepack.topology.diagnostics import _bridge_side_of, _side_of
from tanglepack.topology.TopologyResults import Hole, PseudoneighborPair
from tanglepack.topology.Trellis import Trellis
from tanglepack.topology.TrellisBranch import TrellisBranch


def _fixed_point(period: int, lambda_u: float) -> FixedPoint:
    """A bare no-inversion fixed point whose full-cycle unstable eigenvalue is ``lambda_u``."""
    return bare_fixed_point(period, beta=lambda_u ** (-1.0 / period))


@pytest.fixture
def stable_line():
    """A trellis holding one stable branch with points at stable cdist 1..6.

    Returns ``(trellis, ids)`` where ``ids[j]`` is the intersection at stable
    cdist ``j + 1``. No holes are punched yet — each test fabricates its own.
    """
    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
    stable = (fp, "stable", 0, 0)
    unstable = (fp, "unstable", 0, 0)
    ids = [
        reg.add_synthetic(
            (float(s), 0.0), unstable_cdist=float(7 - s), stable_cdist=float(s),
            manifold_a_key=unstable, manifold_b_key=stable,
        )
        for s in range(1, 7)
    ]
    branch = TrellisBranch(
        key=stable, fixed_point=fp, stability="stable",
        orbit_index=0, branch_index=0, intersection_ids=list(ids),
    )
    trellis = Trellis(
        fixed_points=[fp], registry=reg, branches={stable: branch}, bridges=[]
    )
    return trellis, ids


def _hole(
    trellis,
    near: int,
    far: int,
    row: str,
    *,
    near_which: str = "outward",
    far_which: str = "anchorward",
    bridge_side: str = "left",
) -> Hole:
    """A fabricated hole for the region bounded by (near, far).

    Openings default to the inward pair (outward of the near bound,
    anchorward of the far bound) on ``row`` — the direct-hole case; pass
    ``near_which``/``far_which`` to fabricate a hole on the other side of
    its bridge. Each opening names its crossing's own stable branch.
    """
    branch = trellis.intersection(near).manifold_b_key
    return Hole(
        coords=(0.0, 0.0), near_intersection_id=near,
        bridge_side=bridge_side, bounding_ids=(near, far),
        openings=[(near, near_which, row, branch), (far, far_which, row, branch)],
    )


def _spans(intervals):
    """Intervals as comparable (lo, hi, closed_lo, closed_hi) tuples."""
    return [
        (iv.lo_cdist, iv.hi_cdist, iv.closed_lo, iv.closed_hi) for iv in intervals
    ]


def test_side_of_cross_product():
    """Left/right follow the standard orientation: standing on the manifold
    facing along the tangent (toward the anchor), positive cross = left hand
    (the author's dynamical-direction convention, 2026-07-09)."""
    tangent = np.array([0.0, -1.0])  # toward the anchor, pointing down
    assert _side_of(tangent, np.array([1.0, 0.0])) == "left"   # facing south, east is left
    assert _side_of(tangent, np.array([-1.0, 0.0])) == "right"
    assert _side_of(tangent, np.array([0.0, -2.0])) is None  # collinear


class _StubBridge:
    """A minimal bridge stand-in: endpoint registry ids and a fixed polyline."""

    def __init__(self, first: int, second: int, points, manifold_key=None) -> None:
        self.first_intersection = first
        self.second_intersection = second
        self.manifold_key = manifold_key
        self._points = np.asarray(points, dtype=np.float64)

    @property
    def partial(self) -> bool:
        return self.first_intersection is None or self.second_intersection is None

    def get_point_array(self, return_nodes: bool = False):
        return self._points


def test_bridge_side_orientation_normalization(stable_line):
    """A hole's bridge side is signed in the DYNAMICAL orientation (increasing
    unstable cdist), independent of the polyline's storage order: the same
    geometry stored root→tail either way classifies a query point identically.
    Here the dynamical direction runs -x (ids[0] at x=1 has the larger
    unstable cdist), so a point above the arc is on the RIGHT."""
    trellis, ids = stable_line  # ids[0]: (1,0), u=6; ids[1]: (2,0), u=5

    arc = [(2.0, 0.0), (1.5, 0.0), (1.0, 0.0)]
    stored_forward = _StubBridge(ids[1], ids[0], arc)         # u_first < u_second
    stored_reversed = _StubBridge(ids[0], ids[1], arc[::-1])  # u_first > u_second

    above, below = np.array([1.5, 1.0]), np.array([1.5, -1.0])
    assert _bridge_side_of(trellis, stored_forward, above) == "right"
    assert _bridge_side_of(trellis, stored_reversed, above) == "right"
    assert _bridge_side_of(trellis, stored_forward, below) == "left"
    assert _bridge_side_of(trellis, stored_reversed, below) == "left"


def test_partition_opens_hole_intervals(stable_line):
    """An interval bounding a hole is open; every other interval is closed."""
    trellis, ids = stable_line
    trellis.holes.append(_hole(trellis, ids[0], ids[1], "right"))

    result = partition_stable_manifold(trellis, (trellis.fixed_points[0], "stable", 0, 0), "right")

    assert _spans(result.intervals) == [
        (0.0, 1.0, True, True),
        (1.0, 2.0, False, False),
        (2.0, 6.0, True, True),
    ]


def test_partition_side_without_holes_is_one_closed_interval(stable_line):
    """The side with no holes partitions into a single closed interval."""
    trellis, ids = stable_line
    trellis.holes.append(_hole(trellis, ids[0], ids[1], "right"))

    result = partition_stable_manifold(trellis, (trellis.fixed_points[0], "stable", 0, 0), "left")

    assert _spans(result.intervals) == [(0.0, 6.0, True, True)]


def test_point_between_two_holes_is_closed_singleton(stable_line):
    """A boundary flanked by hole regions on both sides becomes [x, x]."""
    trellis, ids = stable_line
    trellis.holes.append(_hole(trellis, ids[0], ids[1], "right"))
    trellis.holes.append(_hole(trellis, ids[1], ids[2], "right"))

    result = partition_stable_manifold(trellis, (trellis.fixed_points[0], "stable", 0, 0), "right")

    assert _spans(result.intervals) == [
        (0.0, 1.0, True, True),
        (1.0, 2.0, False, False),
        (2.0, 2.0, True, True),
        (2.0, 3.0, False, False),
        (3.0, 6.0, True, True),
    ]


def test_partition_of_a_cut_trellis_ends_at_the_pip():
    """The partition reads the trellis cut at the chosen pip's orbit (D10): it ends
    at the pip, and the crossings past it are no boundaries. (Until 2026-10-08 the
    pip only moved the unified set and the partition ran past it.)"""
    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
    stable, unstable = (fp, "stable", 0, 0), (fp, "unstable", 0, 0)

    def crossing(s: float, u: float) -> int:
        return reg.add_synthetic(
            (s, u), unstable_cdist=u, stable_cdist=s, manifold_a_key=unstable, manifold_b_key=stable
        )

    ids = [crossing(float(s), float(7 - s)) for s in range(1, 7)]
    pip = crossing(4.5, 1.5)
    reg.register_iterate(pip, 1, ids[0])  # f(q0): the unstable cut lies past every crossing
    branch = TrellisBranch(
        key=stable, fixed_point=fp, stability="stable", orbit_index=0, branch_index=0,
        intersection_ids=sorted([*ids, pip], key=lambda i: reg[i].stable_cdist),
    )
    trellis = Trellis(fixed_points=[fp], registry=reg, branches={stable: branch}, bridges=[])
    trellis.strong_pip = pip

    trellis.holes.append(_hole(trellis, ids[0], ids[1], "right"))
    trellis.holes.append(_hole(trellis, ids[2], ids[3], "right"))

    result = partition_stable_manifold(trellis, stable, "right")

    assert _spans(result.intervals) == [
        (0.0, 1.0, True, True),
        (1.0, 2.0, False, False),
        (2.0, 3.0, True, True),
        (3.0, 4.0, False, False),
        (4.0, 4.5, True, True),
    ]


def test_hole_at_branch_end_pinches_endpoint_singleton(stable_line):
    """A hole abutting the outermost intersection opens the interior interval
    and leaves the branch-end point as a closed singleton."""
    trellis, ids = stable_line
    branch_key = (trellis.fixed_points[0], "stable", 0, 0)
    trellis.holes.append(_hole(trellis, ids[4], ids[5], "left"))

    result = partition_stable_manifold(trellis, branch_key, "left")

    assert _spans(result.intervals) == [
        (0.0, 5.0, True, True),
        (5.0, 6.0, False, False),
        (6.0, 6.0, True, True),
    ]


def test_all_holes_participate_regardless_of_flags(stable_line):
    """Every hole on a side shapes the partition identically — the iterate
    label is descriptive only (a zone-membership gate split otherwise
    identical holes by nudge luck and was removed)."""
    trellis, ids = stable_line
    branch_key = (trellis.fixed_points[0], "stable", 0, 0)
    deep = _hole(trellis, ids[0], ids[1], "left")
    deep.iterate = -9
    trellis.holes.append(deep)

    result = partition_stable_manifold(trellis, branch_key, "left")

    assert _spans(result.intervals) == [
        (0.0, 1.0, True, True),
        (1.0, 2.0, False, False),
        (2.0, 6.0, True, True),
    ]


def test_partition_warns_without_pseudoneighbors(stable_line, caplog):
    """Punching holes or partitioning a trellis whose fixed point has no
    recorded pseudoneighbors flags the missing compute_pseudoneighbors call."""
    trellis, _ids = stable_line
    branch_key = (trellis.fixed_points[0], "stable", 0, 0)

    with caplog.at_level("WARNING", logger="tanglepack.topology.Trellis"):
        trellis.punch_holes()
    assert_logged(caplog, logging.WARNING, "tanglepack.topology.Trellis")
    caplog.clear()
    with caplog.at_level("WARNING", logger="tanglepack.topology.Trellis"):
        trellis.partition_stable_manifold(branch_key)
    assert_logged(caplog, logging.WARNING, "tanglepack.topology.Trellis")


def test_shared_boundary_never_closed_on_both_sides(stable_line):
    """A partition assigns each point to exactly one piece: adjacent intervals
    must never both include their shared boundary ("]["). Interleaved hole
    regions (1,4) and (2,5) exercise the worst case: each piece is open
    exactly at the hole-facing sides of the regions' bounding intersections,
    and every shared point is owned by exactly one piece."""
    trellis, ids = stable_line
    branch_key = (trellis.fixed_points[0], "stable", 0, 0)
    trellis.holes.append(_hole(trellis, ids[0], ids[3], "left"))
    trellis.holes.append(_hole(trellis, ids[1], ids[4], "left"))

    result = partition_stable_manifold(trellis, branch_key, "left")

    for previous, current in zip(result.intervals, result.intervals[1:]):
        assert not (previous.closed_hi and current.closed_lo)
    assert _spans(result.intervals) == [
        (0.0, 1.0, True, True),
        (1.0, 2.0, False, True),   # inside (1,4) only at its lower end
        (2.0, 4.0, False, False),  # inside BOTH regions: open piece
        (4.0, 5.0, True, False),   # inside (2,5) only at its upper end
        (5.0, 6.0, True, True),
    ]


def test_two_holes_flanking_one_bridge_pinch_singletons(stable_line):
    """Holes on BOTH sides of one bridge (the author's "singleton bridge")
    emit complementary openings at both defining intersections, pinching each
    into a closed singleton with open intervals on both of its sides."""
    trellis, ids = stable_line
    branch_key = (trellis.fixed_points[0], "stable", 0, 0)
    trellis.holes.append(_hole(trellis, ids[1], ids[2], "left", bridge_side="left"))
    trellis.holes.append(
        _hole(
            trellis, ids[1], ids[2], "left",
            near_which="anchorward", far_which="outward", bridge_side="right",
        )
    )

    result = partition_stable_manifold(trellis, branch_key, "left")

    assert _spans(result.intervals) == [
        (0.0, 2.0, True, False),
        (2.0, 2.0, True, True),
        (2.0, 3.0, False, False),
        (3.0, 3.0, True, True),
        (3.0, 6.0, False, True),
    ]


def test_partition_requires_a_stable_branch(stable_line):
    """A key that is not a stable branch of the trellis raises."""
    trellis, _ids = stable_line
    with pytest.raises(ValueError):
        partition_stable_manifold(
            trellis, (trellis.fixed_points[0], "unstable", 0, 0), "left"
        )


def test_trellis_partition_of_one_branch_is_stored(stable_line):
    """``Trellis.partition_stable_manifold(branch)`` returns one result per side and keeps them."""
    trellis, ids = stable_line
    branch_key = (trellis.fixed_points[0], "stable", 0, 0)
    trellis.holes.append(_hole(trellis, ids[0], ids[1], "right"))
    results = trellis.partition_stable_manifold(branch_key)

    assert sorted(result.side for result in results) == ["left", "right"]
    assert trellis.stable_partitions == results


# --------------------------------------------------------------------------- #
# Real-tangle pipeline
# --------------------------------------------------------------------------- #
@pytest.fixture
def henon_with_holes(henon_tangle_with_bridges):
    """``(trellis, references)`` with pseudoneighbors computed and all holes
    (direct + backward-propagated) punched on the real tangle."""
    workbench, fp = henon_tangle_with_bridges
    trellis = Trellis.from_workbench(workbench, fp)
    trellis.classify_strong_pips()  # the pip's cut starts the reference window
    references = trellis.compute_pseudoneighbors()
    trellis.punch_holes()
    return trellis, references


def test_henon_propagated_hole_side_matches_coords(henon_with_holes):
    """A propagated hole's coordinates are nudged toward the carried point's
    side of the containing bridge, so re-classifying them reproduces the
    stored bridge_side."""
    trellis, _references = henon_with_holes

    checked = 0
    for hole in trellis.holes:
        if hole.pair is not None or hole.bounding_ids is None:
            continue  # not a propagated hole
        wanted = set(hole.bounding_ids)
        bridge = next(
            b
            for b in trellis.bridges
            if {b.first_intersection, b.second_intersection} == wanted
        )
        recomputed = _bridge_side_of(trellis, bridge, np.asarray(hole.coords))
        assert recomputed == hole.bridge_side
        checked += 1
    assert checked > 0, "propagation should have punched holes"


def test_describe_reports(henon_with_holes):
    """The verbose reports are plain strings built from the trellis state."""
    trellis, references = henon_with_holes
    trellis.partition_stable_manifold()

    assert references
    assert trellis.describe_pseudoneighbors()
    assert trellis.describe_holes()
    assert trellis.describe_stable_partitions()


# --------------------------------------------------------------------------- #
# Openings at the two ends of a branch (plan row 1.14)
# --------------------------------------------------------------------------- #
def _end_hole(trellis, iid: int, which: str, row: str) -> Hole:
    """A hole recording a single opening, at one end of the branch."""
    return Hole(
        coords=(0.0, 0.0), near_intersection_id=iid,
        bridge_side="left", bounding_ids=(iid, iid),
        openings=[(iid, which, row, trellis.intersection(iid).manifold_b_key)],
    )


def test_outward_opening_at_the_last_boundary_is_inert(stable_line):
    """A hole hanging off the removed tail leaves the branch end closed.

    Dev Notes: "a hole marked outward of the branch's outermost intersection
    or anchorward at the anchor-artifact contributes its opening to that pinch
    test even though no interval lies there." A boundary is owned by an
    interval only when that interval is closed at it, and the outward side of
    the outermost intersection holds no interval at all — so an opening there
    excludes nothing that exists, and the point stays inside the one interval
    it does have. The pinch at this end is driven by the ANCHORWARD opening
    (see test_hole_at_branch_end_pinches_endpoint_singleton); this pins that
    the outward one cannot pinch it on its own.
    """
    trellis, ids = stable_line
    branch_key = (trellis.fixed_points[0], "stable", 0, 0)
    trellis.holes.append(_end_hole(trellis, ids[5], "outward", "left"))

    result = partition_stable_manifold(trellis, branch_key, "left")

    assert _spans(result.intervals) == [(0.0, 6.0, True, True)]


def test_anchorward_opening_at_the_first_boundary_is_inert(stable_line):
    """A hole abutting the anchor artifact from below leaves it closed.

    The mirror of the branch-end case: nothing lies anchorward of the anchor
    artifact, so an opening there excludes nothing that exists and the point
    stays inside the first real interval. The anchor pinches into a singleton
    only when a hole opens OUTWARD of it.
    """
    trellis, ids = stable_line
    fp = trellis.fixed_points[0]
    branch_key = (fp, "stable", 0, 0)
    # The anchor artifact: the crossing the manifolds make at the fixed point.
    anchor = trellis.registry.add_synthetic(
        (0.0, 0.0), unstable_cdist=0.0, stable_cdist=0.0,
        manifold_a_key=(fp, "unstable", 0, 0), manifold_b_key=branch_key,
    )
    trellis.branch(branch_key).intersection_ids.insert(0, anchor)
    trellis.holes.append(_end_hole(trellis, anchor, "anchorward", "left"))

    result = partition_stable_manifold(trellis, branch_key, "left")

    assert _spans(result.intervals) == [(0.0, 6.0, True, True)]

    # ...whereas an outward opening at the same point does pinch it.
    trellis.holes[:] = [_end_hole(trellis, anchor, "outward", "left")]
    pinched = partition_stable_manifold(trellis, branch_key, "left")
    assert _spans(pinched.intervals) == [
        (0.0, 0.0, True, True),
        (0.0, 6.0, False, True),
    ]


# --------------------------------------------------------------------------- #
# Propagation region identity (plan row 1.15)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "builder", [lambda: build_k28(blasts=1), build_period3], ids=["k28_one_blast", "p3"]
)
def test_both_holes_of_a_singleton_bridge_survive_propagation(builder):
    """A bridge can bound a propagated hole on EACH of its sides (the author's
    "singleton bridge with a hole on either side"): propagation keeps both.

    A bridge bounds a region on each side, so the two defining crossings alone
    do not name a region; a de-duplication keyed on them would drop the second
    of the two. On these cases two orbits land on opposite sides of one bridge.
    """
    case = builder()
    found = 0
    for fp in case.fixed_points:
        by_bridge: dict[tuple, list] = {}
        for hole in case.session.trellis(fp).holes:
            if hole.bounding_ids is not None:
                by_bridge.setdefault(tuple(sorted(hole.bounding_ids)), []).append(hole)
        for holes in by_bridge.values():
            if {hole.bridge_side for hole in holes} == {"left", "right"}:
                found += 1
                assert all(hole.iterate is not None and hole.iterate < 0 for hole in holes)
    assert found, "the case must hold a bridge with a propagated hole on either side"


# --------------------------------------------------------------------------- #
# Cross-branch bounds (plan row 1.20)
# --------------------------------------------------------------------------- #
class _StubNode:
    """A manifold node: a phase-space point at a canonical distance."""

    def __init__(self, point, cdist: float) -> None:
        self._point = np.asarray(point, dtype=np.float64)
        self.cdist = float(cdist)

    def get_point(self):
        return self._point


class _StubManifold:
    """A stable branch running far from the chord of its own crossings."""

    def __init__(self, nodes) -> None:
        self._nodes = list(nodes)

    def get_point_array(self, return_nodes: bool = False):
        if return_nodes:
            return self._nodes
        return np.array([n.get_point() for n in self._nodes])


@pytest.fixture
def stable_line_with_manifold(stable_line):
    """``(trellis, ids, foreign_id)``: ``stable_line`` plus a real stable curve
    that bows far off the crossings' chord, and one crossing attributed to a
    SECOND stable branch (stable cdist 5, unstable cdist 0.5 — so the stable
    and unstable orderings against ``ids[0]`` disagree)."""
    trellis, ids = stable_line
    fp = trellis.fixed_points[0]
    branch_key = (fp, "stable", 0, 0)
    bowed = _StubManifold(
        [_StubNode((float(x), 100.0), float(x)) for x in range(1, 7)]
    )
    trellis.manifolds[branch_key] = bowed
    trellis.manifolds[(fp, "stable", 1, 0)] = bowed  # both branches have nodes
    foreign = trellis.registry.add_synthetic(
        (10.0, -10.0), unstable_cdist=0.5, stable_cdist=5.0,
        manifold_a_key=(fp, "unstable", 0, 0),
        manifold_b_key=(fp, "stable", 1, 0),
    )
    return trellis, ids, foreign


def test_near_far_orders_same_branch_bounds_by_stable_cdist(stable_line_with_manifold):
    """On one branch the toward-anchor member is the smaller stable cdist."""
    from tanglepack.topology.StablePartition import _near_far

    trellis, ids, _foreign = stable_line_with_manifold
    assert _near_far(trellis, ids[3], ids[1]) == (ids[1], ids[3])


def test_near_far_falls_back_to_unstable_order_across_branches(
    stable_line_with_manifold,
):
    """Stable cdists on two different branches are not comparable, so the
    ordering falls back to the unstable dynamical direction — which here
    reverses the (meaningless) stable answer."""
    from tanglepack.topology.StablePartition import _near_far

    trellis, ids, foreign = stable_line_with_manifold
    # Stable cdists say ids[0] (1.0) < foreign (5.0); unstable cdists say
    # foreign (0.5) < ids[0] (6.0). The unstable order is the one used.
    assert _near_far(trellis, ids[0], foreign) == (foreign, ids[0])
    assert _near_far(trellis, foreign, ids[0]) == (foreign, ids[0])


# --------------------------------------------------------------------------- #
# Map orientation (plan row 0.3 review finding)
# --------------------------------------------------------------------------- #
def test_orientation_preserved_for_a_positive_determinant(stable_line):
    """det J > 0 preserves orientation (the Hénon fixtures' b = 1 case)."""
    trellis, _ids = stable_line
    trellis.fixed_points[0].partial_jacobians[0] = np.array([[2.0, 1.0], [-1.0, 0.0]])  # det +1

    assert trellis.orientation_preserving is True


def test_orientation_reversed_for_a_negative_determinant(stable_line):
    """det J < 0 reverses orientation, so a hole's side flips every step."""
    trellis, _ids = stable_line
    trellis.fixed_points[0].partial_jacobians[0] = np.array([[0.0, 1.0], [1.0, 0.0]])  # det -1

    assert trellis.orientation_preserving is False


def test_orientation_is_read_off_the_fixed_point_not_the_map(stable_line):
    """The trellis never evaluates the map's Jacobian: a dynamical system whose
    jacobian contradicts the fixed point's single-step Jacobian is ignored."""
    trellis, _ids = stable_line
    trellis.fixed_points[0].partial_jacobians[0] = np.array([[0.0, 1.0], [1.0, 0.0]])  # det -1
    trellis.dynamical_system = type("System", (), {"jacobian": staticmethod(lambda p: np.eye(2))})()

    assert trellis.orientation_preserving is False


def test_orientation_defaults_to_preserving_without_evidence(stable_line):
    """An unsolved fixed point (no single-step Jacobian) defaults to preserving."""
    trellis, _ids = stable_line
    assert trellis.fixed_points[0].preserves_orientation is None
    assert trellis.orientation_preserving is True


# --------------------------------------------------------------------------- #
# Re-punching (plan row 1.16)
# --------------------------------------------------------------------------- #
def test_repunching_a_narrower_scope_clears_the_other_pairs_holes(henon_with_holes):
    """Punching again drops every hole outside the new scope.

    ``Trellis.punch_holes`` clears ``trellis.holes`` wholesale, so a pair
    outside the new scope must not keep pointing at a hole the trellis no
    longer carries: ``propagate_reference_holes`` seeds its already-punched
    orbits and regions from ``pair.hole``, and a stale one makes it propagate
    an orbit that was never punched. This tangle has a single pseudoneighbor
    pair, so the empty scope is its only narrower one.
    """
    trellis, _references = henon_with_holes
    punched = [p for p in trellis.pseudoneighbors if p.hole is not None]
    assert punched and trellis.holes

    trellis.punch_holes(pairs=[])

    assert all(pair.hole is None for pair in trellis.pseudoneighbors)
    assert trellis.holes == []

    # ...and the full scope punches them again.
    trellis.punch_holes()
    assert [p for p in trellis.pseudoneighbors if p.hole is not None] == punched
    assert trellis.holes


# --------------------------------------------------------------------------- #
# The containing bridge is a bisect over the ray's crossing order (PERF-4)
# --------------------------------------------------------------------------- #
def test_containing_bridge_brackets_by_bisect():
    """Bridges tile an unstable branch between consecutive crossings, so the
    container of a backward image is the bridge whose two crossings bracket its
    span, on the image's own branch; a span straddling a crossing has none."""
    from tanglepack.topology.StablePartition import _containing_bridge

    fp = _fixed_point(3, 4.0)
    reg = IntersectionRegistry()
    stable = (fp, "stable", 0, 0)
    rays = {(fp, "unstable", 0, 0): (0.5, 1.0, 2.0, 3.0), (fp, "unstable", 1, 0): (0.5, 1.5)}
    branches, bridges = {}, {}
    for key, cdists in rays.items():
        ids = [
            reg.add_synthetic(
                (float(key[2]), u), unstable_cdist=u, stable_cdist=u,
                manifold_a_key=key, manifold_b_key=stable,
            )
            for u in cdists
        ]
        branches[key] = TrellisBranch(
            key=key, fixed_point=fp, stability="unstable",
            orbit_index=key[2], branch_index=0, intersection_ids=ids,
        )
        for lo, hi in zip(ids, ids[1:]):
            bridges[key, lo] = _StubBridge(lo, hi, [(0.0, 0.0), (1.0, 0.0)], manifold_key=key)
    trellis = Trellis(
        fixed_points=[fp], registry=reg, branches=branches, bridges=list(bridges.values())
    )
    ray0, ray1 = rays
    first0 = branches[ray0].intersection_ids

    orders: dict = {}
    assert _containing_bridge(trellis, (1.2, 1.8), ray0, orders) is bridges[ray0, first0[1]]
    assert _containing_bridge(trellis, (1.0, 1.5), ray0, orders) is bridges[ray0, first0[1]]
    assert _containing_bridge(trellis, (1.5, 2.5), ray0, orders) is None  # straddles 2.0
    assert _containing_bridge(trellis, (3.2, 3.5), ray0, orders) is None  # past the last
    assert _containing_bridge(trellis, (0.2, 0.4), ray0, orders) is None  # before the first
    first1 = branches[ray1].intersection_ids[0]
    assert _containing_bridge(trellis, (0.6, 0.9), ray1, orders) is bridges[ray1, first1]


# --------------------------------------------------------------------------- #
# The measured side reads the nearest segment of the image sub-arc (INT-5)
# --------------------------------------------------------------------------- #
class _Node:
    """A bridge node: a point at an unstable canonical distance."""

    def __init__(self, point, cdist: float) -> None:
        self._point = np.asarray(point, dtype=np.float64)
        self.cdist = float(cdist)

    def get_point(self):
        return self._point


class _NodeBridge(_StubBridge):
    """A stub bridge that also hands out its nodes (the sub-arc slice reads them)."""

    def __init__(self, first: int, second: int, nodes) -> None:
        super().__init__(first, second, [node.get_point() for node in nodes])
        self._nodes = list(nodes)

    def get_point_array(self, return_nodes: bool = False):
        return self._nodes if return_nodes else self._points


def test_side_at_nearest_segment_reads_two_nodes(stable_line):
    """A propagated hole's sub-arc can lie inside one segment of its bridge: the
    slice then holds two nodes, and the side is still read there (its nearest
    segment), not left undecided. The dynamical direction runs -x here (ids[1]
    at x=2 has the smaller unstable cdist), so a point above the arc is RIGHT."""
    trellis, ids = stable_line
    bridge = _NodeBridge(ids[1], ids[0], [_Node((2.0, 0.0), 5.0), _Node((1.0, 0.0), 6.0)])
    span = (5.2, 5.8)
    assert _bridge_side_of(trellis, bridge, np.array([1.5, 1.0]), span) == "right"
    assert _bridge_side_of(trellis, bridge, np.array([1.5, -1.0]), span) == "left"


# --------------------------------------------------------------------------- #
# Holes propagate backward only (author's rule, 2026-09-16)
# --------------------------------------------------------------------------- #
class _FakePair:
    """The two attributes _is_forward_iterate reads off a pair."""

    def __init__(self, iterate, fixed_point):
        self.iterate = iterate
        self.branch_key = (fixed_point, "stable", 0, 0)


@pytest.mark.parametrize(
    ("period", "iterate", "dropped"),
    [
        (1, 0, False),   # a reference
        (1, -1, False),  # not forward: a backward pair is asserted away instead
        (1, -7, False),
        (1, 1, True),
        (1, 2, True),
        (3, -4, False),
        (3, 1, True),    # +1 .. +(k-1) on the other branches: no hole either
        (3, 2, True),
        (3, 3, True),
        (3, None, False),
    ],
)
def test_forward_pairs_get_no_hole(period, iterate, dropped):
    """Holes propagate backward only: every forward iterate is the image of a
    hole still attached to the stable manifold and is never punched, the
    ``+1 .. +(k_value - 1)`` appearances on the other branches of a period-k
    orbit included (author, 2026-10-09); references always are."""
    from tanglepack.topology.StablePartition import _is_forward_iterate

    fp = _fixed_point(period, 4.0)
    assert fp.k_value == period
    pair = _FakePair(iterate, fp)

    assert _is_forward_iterate(pair) is dropped




def test_a_backward_pair_is_never_punched(stable_line):
    """No backward pair is ever recorded (a reference lies in the window
    ``S \\ f(S)``, so the backward walk leaves T at once), and ``punch_holes``
    asserts it: the backward punch it no longer has is falsifiable."""
    from tanglepack.topology.StablePartition import punch_holes

    trellis, ids = stable_line
    pair = PseudoneighborPair(
        ids[0], ids[1], branch_key=(trellis.fixed_points[0], "stable", 0, 0),
        iterate=-1, origin=(ids[2], ids[3]),
    )

    with pytest.raises(AssertionError):
        punch_holes(trellis, [pair])


# --------------------------------------------------------------------------- #
# A scaled span end next to a crossing is flagged (NUM-5)
# --------------------------------------------------------------------------- #
class _IdBridge(_NodeBridge):
    """A node bridge on one unstable branch, with its BridgeId."""

    def __init__(self, first: int, second: int, nodes, manifold_key) -> None:
        super().__init__(first, second, nodes)
        self.manifold_key = manifold_key
        self.id = (first, second)


@pytest.mark.parametrize("from_table", [False, True], ids=["scaled", "table"])
def test_scaled_span_end_beside_a_crossing_warns(from_table, caplog):
    """A backward image's span end that came from the scaling law is good to
    ``SCALING_RTOL`` only, so one that close to a crossing of its branch may
    have picked the neighbour container: a WARNING. The same end read off the
    iterate table is exact (it is then the container's own end): no WARNING."""
    from tanglepack.topology.StablePartition import propagate_reference_holes
    from tanglepack.topology.Trellis import SCALING_RTOL

    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
    stable, unstable = (fp, "stable", 0, 0), (fp, "unstable", 0, 0)
    cdists = (1.0, 2.0, 4.0, 8.0)
    ids = [
        reg.add_synthetic(
            (u, 0.0), unstable_cdist=u, stable_cdist=u,
            manifold_a_key=unstable, manifold_b_key=stable,
        )
        for u in cdists
    ]
    for i, iid in enumerate(ids):
        reg[iid].crossing_sign = 1 if i % 2 == 0 else -1
    branches = {
        unstable: TrellisBranch(
            key=unstable, fixed_point=fp, stability="unstable",
            orbit_index=0, branch_index=0, intersection_ids=list(ids),
        )
    }
    bridges = [
        _IdBridge(lo, hi, [_Node((u, 1.0), u), _Node((v, 1.0), v)], unstable)
        for (lo, u), (hi, v) in zip(zip(ids, cdists), zip(ids[1:], cdists[1:]))
    ]
    trellis = Trellis(fixed_points=[fp], registry=reg, branches=branches, bridges=bridges)
    reference = PseudoneighborPair(
        ids[2], ids[3], branch_key=stable, is_reference=True, iterate=0, origin=(ids[2], ids[3])
    )
    reference.hole = Hole(
        coords=(6.0, 0.5), near_intersection_id=ids[2], pair=reference,
        bridge_side="left", bounding_ids=(ids[2], ids[3]), iterate=0, origin=reference.origin,
    )
    # One step back the reference bridge lands inside (ids[0], ids[1]), its far end
    # half a SCALING_RTOL short of ids[1]; the next step leaves the branch.
    landing = {ids[2]: (1.2, True), ids[3]: (2.0 * (1.0 - 0.5 * SCALING_RTOL), from_table)}
    trellis.image_cdist = lambda iid, n, stability, branch=None: (
        (unstable, *landing[iid]) if n == -1 else (unstable, 0.1, True)
    )
    caplog.set_level(logging.WARNING, logger="tanglepack.topology.StablePartition")

    holes = propagate_reference_holes(trellis, [reference])

    assert [hole.bounding_ids for hole in holes] == [(ids[0], ids[1])]
    assert_logged(
        caplog, "WARNING", "tanglepack.topology.StablePartition", count=0 if from_table else 1
    )
