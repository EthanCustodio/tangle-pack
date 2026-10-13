"""Which detections are crossings of N (A.2, Q13 (a), B1, 2026-10-12).

Every detection is registered. It is a crossing of N when its resolution
ratio is at least the threshold (clause 1) and no neighbouring detection of the
same two branches sits closer along the unstable polyline than the larger of
their two chord deviations (clause 2: a run is out whole, median included).
An unresolved detection lies inside a bridge, is never linked and is named in
one WARNING. Hand-placed polylines on a real workbench (identity map, never
called by the topology): rules only, no case's numbers.
"""

from __future__ import annotations

import logging

import numpy as np

from helpers.fakes import bare_fixed_point
from helpers.logs import assert_logged
from helpers.polylines import polyline_workbench

WORKBENCH = "tanglepack.numerics.TangleWorkbench"
THRESHOLD = 3.0

#: The x-axis, in segments longer than anything that crosses it.
STABLE = [(-9.0, 0.0), (-3.0, 0.0), (3.0, 0.0), (9.0, 0.0)]


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
    # Two clean crossings at x = -4.5 and 4.5, and between them a zig-zag of two
    # short segments crossing the axis three times within their chord deviation.
    unstable = [
        (-8, 2), (-6, 1), (-5, 0.5), (-4, -0.5), (-3, -1),
        (-1, -1), (0, 1e-3), (1e-3, -1e-3), (1, 1),
        (3, 1), (4, 0.5), (5, -0.5), (6, -1), (8, -2),
    ]
    workbench, fp = _workbench(unstable)
    caplog.set_level(logging.WARNING, logger=WORKBENCH)
    registry = workbench.compute_intersections([fp])
    zig_zag = _detections(registry, near=1.0)
    clean = [iid for iid in _detections(registry) if iid not in zig_zag]

    assert len(zig_zag) == 3 and len(clean) == 2
    resolved = registry.resolved_ids(THRESHOLD)
    assert not set(zig_zag) & resolved
    assert set(clean) <= resolved
    assert sorted(registry.unresolved_crossings(THRESHOLD)) == sorted(zig_zag)

    workbench.create_bridges()
    assert [bridge.id for bridge in workbench.bridges] == [tuple(clean)]

    table = registry.iterate_table
    assert all(table[iid, n] is None for iid in zig_zag for n in (1, -1))
    assert all(table[iid, 1] == iid for iid in clean), "the identity map links the clean ones"
    assert not {table[iid, n] for iid in clean for n in (1, -1)} & set(zig_zag)

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
