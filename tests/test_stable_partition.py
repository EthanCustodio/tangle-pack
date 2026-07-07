"""Hole punching, backward propagation, and the stable-manifold partition.

The synthetic tests fabricate holes directly on a hand-built trellis to pin the
interval logic (open/closed ends, singletons, the resonance-zone cut, the
first-exterior-hole rule) without any manifold numerics; the Hénon tests run
the real punch → propagate → partition pipeline on a computed tangle.
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


def _hole(near: int, far: int, side: str, interior: bool = None) -> Hole:
    """A fabricated hole for the region bounded by (near, far)."""
    return Hole(
        coords=(0.0, 0.0), near_intersection_id=near,
        side=side, bounding_ids=(near, far), interior=interior,
    )


def _spans(intervals):
    """Intervals as comparable (lo, hi, closed_lo, closed_hi) tuples."""
    return [
        (iv.lo_cdist, iv.hi_cdist, iv.closed_lo, iv.closed_hi) for iv in intervals
    ]


def test_side_of_cross_product():
    """Left/right follow the sign of cross(tangent, displacement): negative
    cross is left (the author's convention)."""
    tangent = np.array([0.0, -1.0])  # toward the anchor, pointing down
    assert _side_of(tangent, np.array([-1.0, 0.0])) == "left"
    assert _side_of(tangent, np.array([1.0, 0.0])) == "right"
    assert _side_of(tangent, np.array([0.0, -2.0])) is None  # collinear


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
    assert result.cut_cdist is None
    assert result.exterior_intervals == []


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


def test_cut_splits_interval_and_keeps_first_exterior_hole(stable_line):
    """The branch splits at the resonance-zone cut (the cut point belongs to
    the interior side) and only the outermost exterior hole participates."""
    trellis, ids = stable_line
    fp = trellis.fixed_points[0]
    branch_key = (fp, "stable", 0, 0)
    # The cut sits at an intersection between the hole regions: stable cdist 2.5.
    cut = trellis.registry.add_synthetic(
        (2.5, 0.0), unstable_cdist=4.5, stable_cdist=2.5,
        manifold_a_key=(fp, "unstable", 0, 0), manifold_b_key=branch_key,
    )
    trellis.branch(branch_key).intersection_ids.insert(2, cut)
    trellis.strong_pip = cut

    trellis.holes.append(_hole(ids[0], ids[1], "right"))  # interior
    trellis.holes.append(_hole(ids[2], ids[3], "right"))  # exterior, inner: dropped
    trellis.holes.append(_hole(ids[4], ids[5], "right"))  # exterior, outermost: kept

    result = partition_stable_manifold(trellis, branch_key, "right")

    assert result.cut_cdist == 2.5
    assert _spans(result.interior_intervals) == [
        (0.0, 1.0, True, True),
        (1.0, 2.0, False, False),
        (2.0, 2.5, True, True),
    ]
    assert _spans(result.exterior_intervals) == [
        (2.5, 5.0, False, True),
        (5.0, 6.0, False, False),
    ]


def test_hole_in_another_zone_is_excluded(stable_line):
    """An out-of-zone hole (interior=False) at a non-reference iterate belongs
    to a different resonance zone and must not open any interval here, while
    an out-of-zone REFERENCE hole (iterate 0) is an exterior hole and does."""
    trellis, ids = stable_line
    branch_key = (trellis.fixed_points[0], "stable", 0, 0)
    other_zone = _hole(ids[0], ids[1], "left", interior=False)
    other_zone.iterate = -1
    trellis.holes.append(other_zone)

    result = partition_stable_manifold(trellis, branch_key, "left")
    assert _spans(result.intervals) == [(0.0, 6.0, True, True)]

    exterior_ref = _hole(ids[2], ids[3], "left", interior=False)
    exterior_ref.iterate = 0
    trellis.holes.append(exterior_ref)

    result = partition_stable_manifold(trellis, branch_key, "left")
    assert _spans(result.intervals) == [
        (0.0, 3.0, True, True),
        (3.0, 4.0, False, False),
        (4.0, 6.0, True, True),
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
    """Every punched hole gets a side, its toward-anchor member as
    near_intersection_id, bounding intersections on one stable branch, and its
    orbit identity (origin + iterate)."""
    trellis, references = henon_with_holes

    assert len(trellis.holes) >= len(references)
    for hole in trellis.holes:
        assert hole.side in ("left", "right")
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
        hole_spans = {
            tuple(sorted(h.bounding_ids, key=lambda i: trellis.intersection(i).stable_cdist))
            for h in trellis.holes
            if h.side == result.side and h.bounding_ids is not None
        }
        # Every open interval is a punched hole region on this side.
        assert open_spans <= hole_spans


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
