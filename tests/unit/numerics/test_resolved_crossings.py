"""Which detections are crossings of N (A.2, Q13 (a), B1, 2026-10-12).

Every detection is registered. It is a crossing of N when its resolution
ratio is at least the threshold (clause 1) and no neighbouring detection of the
same two branches sits closer along the unstable polyline than the larger of
their two chord deviations (clause 2: a run is out whole, median included).
An unresolved detection lies inside a bridge, is never linked and is named in
one WARNING. Hand-placed polylines on a real workbench (identity map, never
called by the topology): rules only, no case's numbers. The law cases are
resolved, so the rule moves nothing on them: N is their whole registry.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

from cases import BUILDERS, LAW_CASES
from helpers.fakes import bare_fixed_point
from helpers.logs import assert_logged
from helpers.polylines import polyline_workbench
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry
from tanglepack.topology.StrongPip import is_strong_pip
from tanglepack.topology.Trellis import Trellis

WORKBENCH = "tanglepack.numerics.TangleWorkbench"
THRESHOLD = 3.0

#: The x-axis, in segments longer than anything that crosses it.
STABLE = [(-9.0, 0.0), (-3.0, 0.0), (3.0, 0.0), (9.0, 0.0)]

#: Two clean crossings at x = -4.5 and 4.5, and between them a zig-zag of two
#: short segments crossing the axis three times within their chord deviation.
ZIG_ZAG = [
    (-8, 2), (-6, 1), (-5, 0.5), (-4, -0.5), (-3, -1),
    (-1, -1), (0, 1e-3), (1e-3, -1e-3), (1, 1),
    (3, 1), (4, 0.5), (5, -0.5), (6, -1), (8, -2),
]


def _workbench(unstable):
    """A workbench holding one unstable polyline and the x-axis as its stable curve."""
    fp = bare_fixed_point(label="A")
    workbench = polyline_workbench(
        fp, {(fp, "unstable", 0, 0): unstable, (fp, "stable", 0, 0): STABLE}
    )
    return workbench, fp


def _detections(registry, near: float = np.inf) -> list[int]:
    """The detected (non-synthetic) crossing ids within ``near`` of the origin in x."""
    return [iid for iid, ix in registry if not ix.is_synthetic and abs(ix.coords[0]) < near]


def test_a_zig_zag_across_the_stable_curve_is_registered_and_out_of_n(caplog):
    workbench, fp = _workbench(ZIG_ZAG)
    caplog.set_level(logging.WARNING, logger=WORKBENCH)
    registry = workbench.compute_intersections([fp])
    zig_zag = _detections(registry, near=1.0)
    clean = [iid for iid in _detections(registry) if iid not in zig_zag]

    assert len(zig_zag) == 3 and len(clean) == 2

    workbench.create_bridges()
    assert [bridge.id for bridge in workbench.bridges] == [tuple(clean)]

    table = registry.iterate_table
    assert all(table[iid, n] is None for iid in zig_zag for n in (1, -1))
    assert all(table[iid, 1] == iid for iid in clean), "the identity map links the clean ones"
    assert not {table[iid, n] for iid in clean for n in (1, -1)} & set(zig_zag)

    resolved = registry.resolved_ids(THRESHOLD)
    assert not set(zig_zag) & resolved
    assert set(clean) <= resolved
    assert sorted(registry.unresolved_crossings(THRESHOLD)) == sorted(zig_zag)

    assert_logged(caplog, logging.WARNING, WORKBENCH, count=1)
    caplog.clear()
    workbench.compute_intersections([fp], preserve_ids=True)
    assert_logged(caplog, logging.WARNING, WORKBENCH, count=0)


def test_a_pair_is_judged_by_the_larger_of_its_two_deviations():
    # One detection on a long segment (large chord deviation), its neighbour on
    # the short segment after the shared node (small deviation). They sit between
    # the two deviations apart: compared to its own deviation alone the short one
    # would read resolved; as a pair, both are out.
    unstable = [(-4, 2), (-2, 1), (0, -1e-3), (4e-3, 1e-3), (2, 1), (4, 2)]
    workbench, fp = _workbench(unstable)
    registry = workbench.compute_intersections([fp], diagnose=False)
    long_one, short_one = sorted(
        _detections(registry), key=lambda iid: registry[iid].chord_deviation, reverse=True
    )
    small, large = registry[short_one].chord_deviation, registry[long_one].chord_deviation
    separation = registry.run_separation(short_one)

    assert separation == registry.run_separation(long_one)
    assert small < separation < large
    assert registry[short_one].resolution_ratio >= THRESHOLD, "clause 1 must not decide it"
    assert not {long_one, short_one} & registry.resolved_ids(THRESHOLD)
    assert (registry.run_partner(short_one), registry.run_partner(long_one)) == (long_one, short_one)


def test_a_clean_transversal_crossing_is_in_n():
    workbench, fp = _workbench([(-1.0, -2.0), (-0.5, -1.0), (0.5, 1.0), (1.0, 2.0)])
    registry = workbench.compute_intersections([fp], diagnose=False)
    (crossing,) = _detections(registry)

    assert registry[crossing].resolution_ratio >= THRESHOLD
    assert registry.run_separation(crossing) is None
    assert crossing in registry.resolved_ids(THRESHOLD)
    assert registry.unresolved_crossings(THRESHOLD) == []


def test_an_unresolved_detection_is_on_no_trellis_branch():
    workbench, fp = _workbench(ZIG_ZAG)
    registry = workbench.compute_intersections([fp], diagnose=False)
    workbench.create_bridges()
    trellis = Trellis.from_workbench(workbench, [fp])
    zig_zag = set(_detections(registry, near=1.0))

    assert trellis.resolved_ids == registry.resolved_ids(THRESHOLD)
    assert all(not zig_zag & set(branch.intersection_ids) for branch in trellis.branches.values())
    assert all(set(branch.intersection_ids) <= trellis.resolved_ids for branch in trellis.branches.values())


@pytest.mark.parametrize("name", LAW_CASES)
def test_every_law_case_crossing_is_in_n(name: str) -> None:
    """The law cases are resolved, so N is the whole registry and the rule moves nothing there."""
    case = BUILDERS[name]()
    assert case.registry.unresolved_crossings(case.workbench.resolution_ratio_threshold) == []


@pytest.mark.parametrize("name", LAW_CASES)
def test_chord_deviation_is_none_exactly_on_synthetic_crossings(name: str) -> None:
    registry = BUILDERS[name]().registry
    assert all((ix.chord_deviation is None) == ix.is_synthetic for _, ix in registry)


def test_an_unresolved_detection_never_sets_the_stable_trim():
    # One clean crossing at x = -4.5, then a zig-zag at x = 6, past it along the
    # stable curve: the trim stops at the outermost crossing of N.
    fp = bare_fixed_point(label="A")
    unstable = [(-8, 2), (-6, 1), (-5, 0.5), (-4, -0.5), (-3, -1), (3, -1), (5, -1),
                (6, 1e-3), (6.001, -1e-3), (7, 1), (8, 2)]
    stable = [(float(x), 0.0) for x in range(-9, 10)]
    workbench = polyline_workbench(
        fp, {(fp, "unstable", 0, 0): unstable, (fp, "stable", 0, 0): stable}
    )
    registry = workbench.compute_intersections([fp], diagnose=False)
    zig_zag = registry.unresolved_crossings(THRESHOLD)
    (clean,) = [iid for iid in _detections(registry) if iid not in zig_zag]

    workbench.trim_stable_manifolds(fp)
    end = workbench.manifolds[fp, "stable", 0, 0].tail.cdist

    assert len(zig_zag) == 3
    assert registry[clean].stable_cdist <= end < min(registry[iid].stable_cdist for iid in zig_zag)


def _grazing_workbench():
    """A coarse parabola grazing the x-axis: its near-parallel chord gives a detection
    below the threshold, the next chord a clean crossing."""
    fp = bare_fixed_point(label="A")
    xs = np.arange(-3, 4) + 0.45
    grazing = np.column_stack([xs, 0.05 * xs**2 - 0.0125])
    line = np.column_stack([np.arange(8) - 3.5, np.zeros(8)])
    workbench = polyline_workbench(
        fp, {(fp, "unstable", 0, 0): line, (fp, "stable", 0, 0): grazing}
    )
    return workbench, fp


def test_a_detection_below_the_threshold_is_neither_source_nor_image():
    # Under the identity map the clean crossing is its own image, uniquely once the
    # grazing detection is no candidate.
    workbench, fp = _grazing_workbench()
    registry = workbench.compute_intersections([fp], diagnose=False)
    (unresolved,) = registry.unresolved_crossings(THRESHOLD)
    (clean,) = [iid for iid in _detections(registry) if iid != unresolved]
    table = registry.iterate_table

    assert registry[unresolved].resolution_ratio < THRESHOLD
    assert table[unresolved, 1] is None and table[unresolved, -1] is None
    assert table[clean, 1] == clean


def test_closing_all_iterates_never_waits_on_an_unresolved_detection():
    # Every crossing of N is linked already, so "all" is closed without a growth
    # round; the unresolved detection, never linked, would hold it to the cap.
    workbench, fp = _grazing_workbench()
    workbench.compute_intersections([fp], diagnose=False)
    assert workbench.grow_until_iterates_closed(fp, "all") == 0


def test_an_unresolved_detection_never_disqualifies_a_strong_pip():
    # q0 and one crossing of its own tangle inside its open box (0, 1) x (0, 1):
    # in N it blocks q0; out of N it is no crossing and blocks nothing.
    fp = bare_fixed_point(label="A")
    registry = IntersectionRegistry()
    keys = {"manifold_a_key": (fp, "unstable", 0, 0), "manifold_b_key": (fp, "stable", 0, 0)}
    q0 = registry.add_synthetic((0.0, 0.0), unstable_cdist=1.0, stable_cdist=1.0, **keys)
    inside = registry.add_synthetic((0.1, 0.1), unstable_cdist=0.5, stable_cdist=0.5, **keys)

    def strong(resolved):
        trellis = Trellis(fixed_points=[fp], registry=registry, branches={}, bridges=[],
                          resolved_ids=resolved)
        return is_strong_pip(trellis, q0)

    assert strong(None).blocking_intersection_id == inside
    assert strong(frozenset({q0})).is_strong_pip
