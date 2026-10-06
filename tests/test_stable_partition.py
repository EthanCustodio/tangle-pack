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

import matplotlib

matplotlib.use("Agg")  # headless: exercise the plot helpers without a display
import matplotlib.pyplot as plt
import numpy as np
import pytest

from cases import build_k28, build_period3
from helpers.logs import assert_logged
from helpers.fakes import bare_fixed_point
from tanglepack.numerics.FixedPoint import FixedPoint
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry
from tanglepack.topology.plotting import plot_stable_partition
from tanglepack.topology.StablePartition import (
    _bridge_side_of,
    _side_of,
    partition_stable_manifold,
)
from tanglepack.topology.TopologyResults import Hole
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
    its bridge.
    """
    return Hole(
        coords=(0.0, 0.0), near_intersection_id=near,
        bridge_side=bridge_side, bounding_ids=(near, far),
        openings=[(near, near_which, row), (far, far_which, row)],
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
    trellis.holes.append(_hole(ids[0], ids[1], "right"))

    result = partition_stable_manifold(trellis, (trellis.fixed_points[0], "stable", 0, 0), "right")

    assert _spans(result.intervals) == [
        (0.0, 1.0, True, True),
        (1.0, 2.0, False, False),
        (2.0, 6.0, True, True),
    ]


def test_partition_side_without_holes_is_one_closed_interval(stable_line):
    """The side with no holes partitions into a single closed interval."""
    trellis, ids = stable_line
    trellis.holes.append(_hole(ids[0], ids[1], "right"))

    result = partition_stable_manifold(trellis, (trellis.fixed_points[0], "stable", 0, 0), "left")

    assert _spans(result.intervals) == [(0.0, 6.0, True, True)]


def test_point_between_two_holes_is_closed_singleton(stable_line):
    """A boundary flanked by hole regions on both sides becomes [x, x]."""
    trellis, ids = stable_line
    trellis.holes.append(_hole(ids[0], ids[1], "right"))
    trellis.holes.append(_hole(ids[1], ids[2], "right"))

    result = partition_stable_manifold(trellis, (trellis.fixed_points[0], "stable", 0, 0), "right")

    assert _spans(result.intervals) == [
        (0.0, 1.0, True, True),
        (1.0, 2.0, False, False),
        (2.0, 2.0, True, True),
        (2.0, 3.0, False, False),
        (3.0, 6.0, True, True),
    ]


def test_partition_ignores_strong_pip_cut(stable_line):
    """The strong-pip cut no longer splits the partition — the intervals run
    contiguously past it, and every hole opens its interval (the earlier
    interior/exterior reporting split was removed with the zone association,
    2026-07-16)."""
    trellis, ids = stable_line
    fp = trellis.fixed_points[0]
    branch_key = (fp, "stable", 0, 0)
    # A strong pip between the hole regions, at stable cdist 2.5.
    cut = trellis.registry.add_synthetic(
        (2.5, 0.0), unstable_cdist=4.5, stable_cdist=2.5,
        manifold_a_key=(fp, "unstable", 0, 0), manifold_b_key=branch_key,
    )
    trellis.branch(branch_key).intersection_ids.insert(2, cut)
    trellis.strong_pip = cut

    trellis.holes.append(_hole(ids[0], ids[1], "right"))
    trellis.holes.append(_hole(ids[2], ids[3], "right"))
    trellis.holes.append(_hole(ids[4], ids[5], "right"))  # outermost

    result = partition_stable_manifold(trellis, branch_key, "right")

    # No boundary at the cut; the region (5, 6) ends at the branch end, so
    # the outermost point owns itself as a singleton.
    assert 2.5 not in {iv.lo_cdist for iv in result.intervals}
    assert _spans(result.intervals) == [
        (0.0, 1.0, True, True),
        (1.0, 2.0, False, False),
        (2.0, 3.0, True, True),
        (3.0, 4.0, False, False),
        (4.0, 5.0, True, True),
        (5.0, 6.0, False, False),
        (6.0, 6.0, True, True),
    ]


def test_hole_at_branch_end_pinches_endpoint_singleton(stable_line):
    """A hole abutting the outermost intersection opens the interior interval
    and leaves the branch-end point as a closed singleton."""
    trellis, ids = stable_line
    branch_key = (trellis.fixed_points[0], "stable", 0, 0)
    trellis.holes.append(_hole(ids[4], ids[5], "left"))

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
    deep = _hole(ids[0], ids[1], "left")
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
    trellis.holes.append(_hole(ids[0], ids[3], "left"))
    trellis.holes.append(_hole(ids[1], ids[4], "left"))

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
    trellis.holes.append(_hole(ids[1], ids[2], "left", bridge_side="left"))
    trellis.holes.append(
        _hole(
            ids[1], ids[2], "left",
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


def test_plot_stable_partition_smoke(stable_line):
    """The number-line plot draws one row per partition and returns the axes."""
    trellis, ids = stable_line
    branch_key = (trellis.fixed_points[0], "stable", 0, 0)
    trellis.holes.append(_hole(ids[0], ids[1], "right"))
    results = trellis.partition_stable_manifold(branch_key)

    fig, ax = plt.subplots()
    try:
        drawn = plot_stable_partition(results, ax=ax)
        assert drawn is ax
        assert len(ax.get_yticklabels()) == 2  # left and right rows
    finally:
        plt.close(fig)

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


def test_henon_plot_helpers_smoke(henon_with_holes):
    """plot_pseudoneighbors and plot_holes draw on a real tangle."""
    trellis, _references = henon_with_holes

    fig, ax = plt.subplots()
    try:
        assert trellis.plot_pseudoneighbors(ax=ax) is not None
        handles = trellis.plot_holes(ax=ax)
        assert handles
    finally:
        plt.close(fig)


# --------------------------------------------------------------------------- #
# Openings at the two ends of a branch (plan row 1.14)
# --------------------------------------------------------------------------- #
def _end_hole(iid: int, which: str, row: str) -> Hole:
    """A hole recording a single opening, at one end of the branch."""
    return Hole(
        coords=(0.0, 0.0), near_intersection_id=iid,
        bridge_side="left", bounding_ids=(iid, iid),
        openings=[(iid, which, row)],
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
    trellis.holes.append(_end_hole(ids[5], "outward", "left"))

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
    trellis.holes.append(_end_hole(anchor, "anchorward", "left"))

    result = partition_stable_manifold(trellis, branch_key, "left")

    assert _spans(result.intervals) == [(0.0, 6.0, True, True)]

    # ...whereas an outward opening at the same point does pinch it.
    trellis.holes[:] = [_end_hole(anchor, "outward", "left")]
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
class _StubSystem:
    """A dynamical system stand-in carrying only a constant jacobian."""

    def __init__(self, jacobian) -> None:
        self._jacobian = np.asarray(jacobian, dtype=np.float64)
        self.jacobian = lambda point: self._jacobian


def test_orientation_preserved_for_a_positive_determinant(stable_line):
    """det J > 0 preserves orientation (the Hénon fixtures' b = 1 case)."""
    trellis, _ids = stable_line
    trellis.fixed_points[0].coordinates[0] = np.zeros((2, 1))
    trellis.dynamical_system = _StubSystem([[2.0, 1.0], [-1.0, 0.0]])  # det +1

    assert trellis.orientation_preserving is True


def test_orientation_reversed_for_a_negative_determinant(stable_line):
    """det J < 0 reverses orientation, so a hole's side flips every step."""
    trellis, _ids = stable_line
    trellis.fixed_points[0].coordinates[0] = np.zeros((2, 1))
    trellis.dynamical_system = _StubSystem([[0.0, 1.0], [1.0, 0.0]])  # det -1

    assert trellis.orientation_preserving is False


def test_orientation_from_eigenvalues_without_a_jacobian(stable_line):
    """Without a jacobian an odd-period saddle still decides it: the product of
    its eigenvalues is the determinant of the p-fold Jacobian, whose sign is
    the single-step sign when p is odd."""
    trellis, _ids = stable_line
    fp = trellis.fixed_points[0]
    assert fp.period == 1
    fp.stable_eigenvalues = [-0.25]
    assert trellis.orientation_preserving is False


def test_orientation_defaults_to_preserving_without_evidence(stable_line):
    """No jacobian and no usable eigenvalue pair defaults to preserving."""
    trellis, _ids = stable_line
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
# A bridge's unstable branch comes from its own key (plan row 1.1)
# --------------------------------------------------------------------------- #
@pytest.fixture
def two_branch_bridges():
    """``(trellis, cycle, bridges)``: two bridges on two unstable branches of a
    period-3 orbit, sharing an anchor artifact at unstable cdist 0.

    The shared anchor's ``manifold_a_key`` names branch 1, so the branch-1
    bridge is correctly labelled by its endpoints but the branch-0 bridge is
    NOT — exactly the situation the after-the-fact endpoint assignment creates
    (all three anchors sit at unstable cdist 0). The branch-0 bridge is also
    the tighter of the two, so an endpoint-key filter picks it as the branch-1
    container.
    """
    fp = _fixed_point(3, 4.0)
    cycle = [(fp, "unstable", i, 0) for i in range(3)]
    reg = IntersectionRegistry()
    stable = (fp, "stable", 0, 0)
    anchor = reg.add_synthetic(
        (0.0, 0.0), unstable_cdist=0.0, stable_cdist=0.0,
        manifold_a_key=cycle[1], manifold_b_key=stable,
    )
    end_0 = reg.add_synthetic(
        (1.0, 0.0), unstable_cdist=1.0, stable_cdist=1.0,
        manifold_a_key=cycle[0], manifold_b_key=stable,
    )
    end_1 = reg.add_synthetic(
        (2.0, 0.0), unstable_cdist=2.0, stable_cdist=2.0,
        manifold_a_key=cycle[1], manifold_b_key=stable,
    )
    bridges = {
        0: _StubBridge(anchor, end_0, [(0.0, 0.0), (1.0, 0.0)], manifold_key=cycle[0]),
        1: _StubBridge(anchor, end_1, [(0.0, 0.0), (2.0, 0.0)], manifold_key=cycle[1]),
    }
    trellis = Trellis(
        fixed_points=[fp], registry=reg, branches={}, bridges=list(bridges.values())
    )
    return trellis, cycle, bridges


def test_bridge_span_reads_the_branch_from_the_bridges_own_key(two_branch_bridges):
    """A bridge's cycle position comes from ``Bridge.manifold_key``, not from
    whichever endpoint happens to carry an unstable key."""
    from tanglepack.topology.StablePartition import _bridge_unstable_span

    trellis, cycle, bridges = two_branch_bridges
    span, pos = _bridge_unstable_span(trellis, bridges[0], cycle)

    assert span == (0.0, 1.0)
    assert pos == 0  # the shared anchor's own key would say 1


def test_containing_bridge_filters_on_the_bridges_own_key(two_branch_bridges):
    """The container of a backward image is looked up on the image's branch.

    The branch-0 bridge spans the query too and is tighter, so it wins on cdist
    evidence; only its own ``manifold_key`` rules it out of a branch-1 lookup.
    """
    from tanglepack.topology.StablePartition import _containing_bridge

    trellis, cycle, bridges = two_branch_bridges

    assert _containing_bridge(trellis, (0.1, 0.5), cycle[1]) is bridges[1]
    assert _containing_bridge(trellis, (0.1, 0.5), cycle[0]) is bridges[0]
    assert _containing_bridge(trellis, (0.1, 0.5), cycle[2]) is None



# --------------------------------------------------------------------------- #
# Backward images of a bridge span: iterate-table lookup first, scaling after
# --------------------------------------------------------------------------- #
def test_backward_endpoint_prefers_the_table_then_scales():
    """A bridge endpoint's preimage is read from the iterate table at its
    REGISTERED unstable cdist (5% off the scaled estimate here, which a
    scaled span would have carried); once the chain ends, scaling by 1/beta
    takes over from the last registered value and the id is dropped."""
    from tanglepack.topology.StablePartition import _backward_endpoint

    fp = _fixed_point(1, 4.0)
    beta = fp.per_step_beta("unstable")
    reg = IntersectionRegistry()
    unstable, stable = (fp, "unstable", 0, 0), (fp, "stable", 0, 0)
    root = reg.add_synthetic(
        (0.0, 0.0), unstable_cdist=2.1, stable_cdist=4.0,
        manifold_a_key=unstable, manifold_b_key=stable,
    )
    end = reg.add_synthetic(
        (0.0, 0.0), unstable_cdist=8.0, stable_cdist=1.0,
        manifold_a_key=unstable, manifold_b_key=stable,
    )
    reg.register_iterate(root, 1, end)
    trellis = Trellis(fixed_points=[fp], registry=reg, branches={}, bridges=[])

    linked = _backward_endpoint(trellis, end, 8.0, beta)
    assert linked == (root, 2.1)  # the registered value, not 8 / 4 = 2

    beyond = _backward_endpoint(trellis, *linked, beta)
    assert beyond == (None, 2.1 / beta)  # chain ended: scaled from 2.1

    assert _backward_endpoint(trellis, *beyond, beta) == (None, 2.1 / beta**2)

    reg2 = IntersectionRegistry()
    unlinked = reg2.add_synthetic(
        (0.0, 0.0), unstable_cdist=8.0, stable_cdist=1.0,
        manifold_a_key=unstable, manifold_b_key=stable,
    )
    trellis2 = Trellis(fixed_points=[fp], registry=reg2, branches={}, bridges=[])
    assert _backward_endpoint(trellis2, unlinked, 8.0, beta) == (None, 2.0)


# --------------------------------------------------------------------------- #
# Holes propagate backward only (author's rule, 2026-09-16)
# --------------------------------------------------------------------------- #
class _FakePair:
    """The two attributes _is_forward_beyond_fundamental reads off a pair."""

    def __init__(self, iterate, fixed_point):
        self.iterate = iterate
        self.branch_key = (fixed_point, "stable", 0, 0)
        self.intersection_a = 0


@pytest.mark.parametrize(
    ("period", "iterate", "dropped"),
    [
        (1, 0, False),   # a reference
        (1, -1, False),  # backward iterates always punch
        (1, -7, False),
        (1, 1, True),    # k_value = 1: no forward hole at all
        (1, 2, True),
        (3, 3, True),    # k_value = 3: a full branch return wraps past the fundamental segment
        (3, 4, True),
        (3, None, False),
    ],
)
def test_forward_pairs_beyond_the_fundamental_segment_get_no_hole(period, iterate, dropped):
    """The FIRM half of the backward-only rule: forward iterates ``>= k_value``
    are the images of holes still attached to the stable manifold and are never
    punched; references and backward iterates always are.

    Note:
        The ``+1 .. +(k_value - 1)`` exemption is PROVISIONAL (author,
        2026-09-16) and deliberately untested here.
    """
    from tanglepack.topology.StablePartition import _is_forward_beyond_fundamental

    fp = _fixed_point(period, 4.0)
    assert fp.k_value == period
    pair = _FakePair(iterate, fp)

    assert _is_forward_beyond_fundamental(None, pair) is dropped


