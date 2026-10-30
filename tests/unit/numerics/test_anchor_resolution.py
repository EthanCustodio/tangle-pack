"""Clause 3 of N, the near-anchor guard (design note §A.2, 2026-10-29).

A detection within ``anchor_resolution_factor`` chord deviations of the orbit
point of one of its branches is registered but no crossing of N: it is listed
unresolved, named in the WARNING with the remedy "refine", and never linked.
Hand-placed polylines on a real workbench (identity map, never called by the
topology): the orbit point is placed by the test, and every assertion compares
a clearance with the factor, never a number of the geometry.
"""

from __future__ import annotations

import logging

from helpers.fakes import bare_fixed_point
from helpers.polylines import polyline_workbench
from tanglepack.numerics.IntersectionRegistry import ANCHOR_RESOLUTION_FACTOR

WORKBENCH = "tanglepack.numerics.TangleWorkbench"
THRESHOLD = 3.0
FACTOR = ANCHOR_RESOLUTION_FACTOR

#: The x-axis, in segments longer than anything that crosses it.
STABLE = [(-9.0, 0.0), (-3.0, 0.0), (3.0, 0.0), (9.0, 0.0)]

#: A bent curve crossing the axis twice, at x = -4.5 and 4.5.
BENT = [(-8, 2), (-6, 1), (-5, 0.5), (-4, -0.5), (-3, -1), (0, -1.5),
        (3, -1), (4, -0.5), (5, 0.5), (6, 1), (8, 2)]


def _workbench(unstable, orbit_point, factor=FACTOR):
    """A workbench of one unstable polyline over the x-axis, its fixed point at ``orbit_point``."""
    fp = bare_fixed_point(label="A", coordinates=[orbit_point])
    workbench = polyline_workbench(
        fp,
        {(fp, "unstable", 0, 0): unstable, (fp, "stable", 0, 0): STABLE},
        anchor_resolution_factor=factor,
    )
    return workbench, fp


def _detections(registry) -> list[int]:
    """The detected crossings, in unstable cdist order."""
    return [iid for iid in registry.by_unstable_cdist if not registry[iid].is_synthetic]


def test_a_detection_nearer_its_orbit_point_than_c_h_is_out_of_n(caplog):
    workbench, fp = _workbench(BENT, (-4.5, 0.05))
    caplog.set_level(logging.WARNING, logger=WORKBENCH)
    registry = workbench.compute_intersections([fp])
    near, far = _detections(registry)

    assert registry.orbit_clearance(near) < FACTOR <= registry.orbit_clearance(far)
    assert min(registry[iid].resolution_ratio for iid in (near, far)) >= THRESHOLD, "clause 1 must not decide"
    assert registry.run_partner(near) is None, "clause 2 must not decide"
    assert near not in registry.resolved_ids(THRESHOLD)
    assert far in registry.resolved_ids(THRESHOLD)
    assert registry.unresolved_crossings(THRESHOLD) == [near]
    (record,) = [r for r in caplog.records if r.levelno == logging.WARNING]
    assert "clause 3" in record.getMessage() and "refine" in record.getMessage()


def test_none_turns_clause_3_off():
    workbench, fp = _workbench(BENT, (-4.5, 0.05), factor=None)
    registry = workbench.compute_intersections([fp], diagnose=False)

    assert set(_detections(registry)) <= registry.resolved_ids(THRESHOLD)
    assert registry.unresolved_crossings(THRESHOLD) == []
    workbench.anchor_resolution_factor = FACTOR
    assert registry.unresolved_crossings(THRESHOLD), "the factor is read on every call"


def test_an_unresolved_image_is_never_linked():
    # The identity map links every crossing of N to itself; one out by clause 3 to nothing.
    workbench, fp = _workbench(BENT, (-4.5, 0.05))
    registry = workbench.compute_intersections([fp], diagnose=False)
    near, far = _detections(registry)
    table = registry.iterate_table

    assert table[far, 1] == far
    assert table[near, 1] is None and table[near, -1] is None


def test_anchors_are_exempt():
    workbench, fp = _workbench(BENT, (-4.5, 0.05))
    registry = workbench.compute_intersections([fp], diagnose=False)
    anchors = [iid for iid, ix in registry if ix.is_synthetic]

    assert all(registry.orbit_clearance(iid) is None for iid in anchors)
    assert set(anchors) <= registry.resolved_ids(THRESHOLD)


def test_the_factor_survives_a_registry_swap():
    workbench, fp = _workbench(BENT, (-4.5, 0.05), factor=7.5)
    workbench.compute_intersections([fp], diagnose=False)
    assert workbench.intersection_registry.anchor_resolution_factor == 7.5


def test_the_nearer_orbit_point_of_either_branch_decides():
    # A heteroclinic-style detection: its unstable branch belongs to a far
    # fixed point, its stable branch to one beside it.
    far, near_fp = (bare_fixed_point(label=label, coordinates=[xy]) for label, xy in
                    (("A", (50.0, 50.0)), ("B", (-4.5, 0.05))))
    workbench = polyline_workbench(
        far, {(far, "unstable", 0, 0): BENT, (near_fp, "stable", 0, 0): STABLE},
        anchor_resolution_factor=FACTOR,
    )
    workbench.fixed_points.append(near_fp)
    registry = workbench.compute_intersections([far, near_fp], diagnose=False)
    near, _far = _detections(registry)
    assert registry.orbit_clearance(near) < FACTOR
    assert near not in registry.resolved_ids(THRESHOLD)
