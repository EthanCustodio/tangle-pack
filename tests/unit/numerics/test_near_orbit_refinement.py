"""The opt-in graded mesh near the periodic orbit (``near_orbit_radius``, 2026-10-29).

Within ``r0`` of an orbit point a pair is refined when its curvature area is
at least ``area_cutoff * (r / r0) ** 2`` or its chord is longer than ``r``
(``r`` how near the chord passes the orbit point), down to a spacing floor; the
seed's segment at the orbit point is graded once at initialisation, carried
inward along the seed's ray. ``None`` is today's refinement exactly, and the
point budget still guards. Checked on the ``k = 2.5`` period-1 saddle of the
scale series at a coarse cutoff: rules only, no count of the build.
"""

from __future__ import annotations

import numpy as np
import pytest

from tanglepack import TangleWorkbench
from tanglepack.examples.henon import henon_jacobian, henon_map, henon_map_inverse
from tanglepack.numerics.budget import PointBudgetExceeded
from tanglepack.numerics.geometry import curvature_area_batch
from tanglepack.numerics.ManifoldMachine import SPACING_FLOOR

K, SEED = 2.5, (4.0, -4.0)
AREA_CUTOFF = 1e-5
R0 = 1e-3


def _workbench(near_orbit_radius="omit", **kwargs):
    """The scale saddle, oriented and seeded; ``near_orbit_radius`` passed unless omitted."""
    if near_orbit_radius != "omit":
        kwargs["near_orbit_radius"] = near_orbit_radius
    workbench = TangleWorkbench(henon_map(K, 1), henon_map_inverse(K, 1), henon_jacobian(K, 1),
                                area_cutoff=AREA_CUTOFF, **kwargs)
    fp = workbench.construct_fixed_point(list(SEED))
    workbench.orient_eigenvectors(fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])})
    workbench.initialize_both_manifolds(fp)
    return workbench, fp


def _floor(fp) -> float:
    return 64 * np.finfo(float).eps * max(1.0, float(np.abs(np.asarray(fp.coordinates)).max()))


def _distances(workbench, fp, stability):
    """Each point's distance to the orbit point, root first."""
    points = workbench.manifolds[(fp, stability, 0, 0)].get_point_array()
    return np.hypot(*(points - fp.coordinates[0]).T)


@pytest.mark.parametrize("stability", ["stable", "unstable"])
def test_the_seed_segment_is_graded_toward_the_orbit_point(stability):
    (workbench, fp), (plain, plain_fp) = _workbench(R0), _workbench(None)
    graded, seed = _distances(workbench, fp, stability), _distances(plain, plain_fp, stability)

    assert len(graded) > len(seed)
    assert np.all(np.diff(graded) > 0), "the points approach the orbit point monotonically"
    assert graded[1] < 1e-6 * seed[1], "the grading reaches far inside the seed's first point"
    assert graded[1] >= _floor(fp), "no point inside the spacing floor"
    inner = graded[1:]
    assert np.all(np.diff(inner) <= inner[:-1] + 1e-12 * inner[:-1]), "no chord longer than its distance"


def test_the_chord_deviation_falls_like_r_squared_near_the_orbit_point():
    workbench, fp = _workbench(R0)
    workbench.grow_n_times(fp, "unstable", num_iterations=4)
    points = workbench.manifolds[(fp, "unstable", 0, 0)].get_point_array()
    p0, p1 = points[1:-2], points[2:-1]
    areas = curvature_area_batch(p0, p1, points[:-3], points[3:])
    r = np.minimum(np.hypot(*(p0 - fp.coordinates[0]).T), np.hypot(*(p1 - fp.coordinates[0]).T))
    chords = np.hypot(*(p1 - p0).T)
    near = (r < R0) & (chords >= np.maximum(SPACING_FLOOR * r / R0, _floor(fp)))
    assert near.sum() > 10
    assert np.all(areas[near] < AREA_CUTOFF * (r[near] / R0) ** 2)


def test_none_is_todays_refinement_byte_for_byte():
    rows = []
    for radius in ("omit", None):
        workbench, fp = _workbench(radius)
        workbench.grow_n_times(fp, "unstable", num_iterations=4)
        workbench.grow_n_times(fp, "stable", num_iterations=4)
        store = workbench.point_store
        rows.append([np.asarray(column).tobytes() for column in (store.xy, store.cdist, store.image, store.preimage)])
    assert rows[0] == rows[1]


def test_the_point_budget_still_raises():
    workbench, fp = _workbench(R0, max_points_per_step=200)
    with pytest.raises(PointBudgetExceeded):
        workbench.grow_n_times(fp, "unstable", num_iterations=12)


def test_the_seed_grading_is_budgeted():
    with pytest.raises(PointBudgetExceeded):
        _workbench(R0, max_points_per_step=5)
