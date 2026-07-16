"""Hole punching, backward propagation, and the stable-manifold partition.

The synthetic tests fabricate holes directly on a hand-built trellis to pin the
interval logic (open/closed ends, singletons, the both-sides pinch) without
any manifold numerics; the Hénon tests run the real punch → propagate →
partition pipeline on a computed tangle.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless: exercise the plot helpers without a display
import matplotlib.pyplot as plt
import numpy as np
import pytest

from tanglepack.numerics.FixedPoint import FixedPoint
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry
from tanglepack.topology.StablePartition import (
    _bridge_side_of,
    _side_of,
    partition_stable_manifold,
    plot_stable_partition,
)
from tanglepack.topology.TopologyResults import Hole
from tanglepack.topology.Trellis import Trellis
from tanglepack.topology.TrellisBranch import TrellisBranch


def _fixed_point(period: int, lambda_u: float) -> FixedPoint:
    """A minimal no-inversion fixed point with a positive unstable eigenvalue."""
    fp = FixedPoint(period, 1)
    fp.unstable_eigenvalues = [lambda_u] * period
    fp.set_k_value()
    return fp


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

    def __init__(self, first: int, second: int, points) -> None:
        self.first_intersection = first
        self.second_intersection = second
        self._points = np.asarray(points, dtype=np.float64)

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
        trellis.partition_stable_manifold(branch_key)

    warnings = [r for r in caplog.records if "compute_pseudoneighbors" in r.message]
    assert len(warnings) == 2


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


def test_henon_holes_are_classified(henon_with_holes):
    """Every punched hole gets a bridge side, its toward-anchor member as
    near_intersection_id, bounding intersections on one stable branch, and its
    orbit identity (origin + iterate). No bridge on this tangle is degenerate,
    so no side may be None."""
    trellis, references = henon_with_holes

    assert len(trellis.holes) >= len(references)
    for hole in trellis.holes:
        assert hole.bridge_side in ("left", "right")
        near, far = hole.bounding_ids
        assert hole.near_intersection_id == near
        s_near = trellis.intersection(near).stable_cdist
        s_far = trellis.intersection(far).stable_cdist
        assert s_near <= s_far
        assert np.isfinite(hole.coords).all()
        assert hole.iterate is not None
        assert hole.origin is not None


def test_henon_reference_holes_hug_the_stable_manifold(henon_with_holes):
    """A reference hole plots near the pair's stable chord; an iterated pair's
    hole plots between the two neighbors on the unstable manifold (the bridge
    midpoint)."""
    from tanglepack.topology.StablePartition import (
        bridge_for_pair,
        _nearest_arc_point,
        _stable_arc_midpoint,
    )

    trellis, _references = henon_with_holes
    checked = 0
    for pair in trellis.pseudoneighbors:
        if pair.hole is None:
            continue
        bridge = bridge_for_pair(trellis, pair)
        a = trellis.intersection(pair.intersection_a)
        b = trellis.intersection(pair.intersection_b)
        near, far = (a, b) if a.stable_cdist <= b.stable_cdist else (b, a)
        stable_mid, _tangent = _stable_arc_midpoint(trellis, near, far)
        arc_point = _nearest_arc_point(bridge, stable_mid)
        coords = np.asarray(pair.hole.coords)
        to_chord = np.linalg.norm(coords - stable_mid)
        to_bridge = np.linalg.norm(coords - arc_point)
        if pair.is_reference:
            assert to_chord <= to_bridge
        else:
            assert to_bridge <= to_chord
        checked += 1
    assert checked > 0


def test_henon_partition_covers_branch(henon_with_holes):
    """Both side partitions run from the anchor to the branch end with
    contiguous intervals, open exactly at the punched hole regions."""
    trellis, _references = henon_with_holes
    fp = trellis.fixed_points[0]
    branch_key = (fp, "stable", 0, 0)

    results = trellis.partition_stable_manifold(branch_key)
    assert [r.side for r in results] == ["left", "right"]

    branch = trellis.branch(branch_key)
    s_max = trellis.intersection(branch.ordered_ids()[-1]).stable_cdist
    for result in results:
        intervals = result.intervals
        assert intervals[0].lo_cdist <= trellis.registry.cdist_tol
        assert intervals[-1].hi_cdist == pytest.approx(s_max)
        for prev, nxt in zip(intervals, intervals[1:]):
            assert prev.hi_cdist == pytest.approx(nxt.lo_cdist)
        open_spans = {
            (iv.lo_id, iv.hi_id)
            for iv in intervals
            if not iv.closed_lo and not iv.closed_hi
        }
        # A hole acts on this row exactly when one of its openings does.
        hole_spans = {
            tuple(sorted(h.bounding_ids, key=lambda i: trellis.intersection(i).stable_cdist))
            for h in trellis.holes
            if h.bounding_ids is not None
            and any(row == result.side for _iid, _which, row in (h.openings or []))
        }
        # Every open interval is a punched hole region on this side.
        assert open_spans <= hole_spans


def test_henon_direct_hole_side_reproduces_inward_pair(henon_with_holes):
    """A direct hole sits on the side of its bridge facing the pair's own
    stable segment, so its bridge-side-derived openings must reproduce the
    inward pair: outward of the near bound, anchorward of the far bound
    (modulo an anchor-artifact end, which opens nothing)."""
    trellis, _references = henon_with_holes

    checked = 0
    for pair in trellis.pseudoneighbors:
        hole = pair.hole
        if hole is None:
            continue
        assert hole.bridge_side in ("left", "right")
        near, far = hole.bounding_ids
        assert hole.openings
        for iid, which, _row in hole.openings:
            assert which == ("outward" if iid == near else "anchorward")
        checked += 1
    assert checked > 0


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

    pairs_report = trellis.describe_pseudoneighbors()
    assert f"{len(references)} reference pseudoneighbor pair(s):" in pairs_report
    assert "hole(s) on the left side" in trellis.describe_holes()
    partition_report = trellis.describe_stable_partitions()
    assert "left partition" in partition_report
    assert "right partition" in partition_report


def test_p3_propagation_terminates_at_periodicity(henon_p3_session):
    """Period-3: one backward step moves a bridge to the previous branch of
    the cycle, so termination must compare bridges at the same cycle residue.
    The regression signature of the broken (consecutive-step) rule is the
    backward orbit cycling through the same bridges at ever deeper iterates —
    i.e. one origin punching the same bounding region more than once."""
    session, fp3, _fp1, _zone = henon_p3_session
    trellis = session.trellis(fp3)
    trellis.classify_strong_pips()
    references = trellis.compute_pseudoneighbors()
    assert references

    trellis.punch_holes()

    propagated = [h for h in trellis.holes if h.iterate not in (0, None)]
    assert propagated, "backward propagation should punch holes"
    for origin in {h.origin for h in propagated}:
        regions = [h.bounding_ids for h in propagated if h.origin == origin]
        assert len(regions) == len(set(regions))


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
