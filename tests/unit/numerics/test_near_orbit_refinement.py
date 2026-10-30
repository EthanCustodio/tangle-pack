"""The opt-in graded mesh near the periodic orbit (``near_orbit_radius``, 2026-10-29).

Within ``r0`` of an orbit point a pair is refined when its curvature area is
at least ``area_cutoff * (r / r0) ** 2`` or its chord is longer than ``r``
(``r`` how near the chord passes the orbit point, so a straight chord between
far nodes passing it is graded too), down to a spacing floor; the
seed's segment at the orbit point is graded once at initialisation, carried
inward along the seed's ray. ``None`` is today's refinement exactly, and the
point budget still guards. Checked on the ``k = 2.5`` period-1 saddle of the
scale series at a coarse cutoff, and the seed's walk along its chain on the
period-3 orbit and the ``k = 4`` inversion saddle (a branch flip on every
return): rules only, no count of the build.
"""

from __future__ import annotations

import numpy as np
import pytest

from cases import HENON_K4, HENON_P3, K4_INVERSION_SEED, P3_ORBIT_SEED, PERIOD3_ORIENTATION
from helpers.polylines import polyline_manifold
from tanglepack import TangleWorkbench
from tanglepack.examples.henon import henon_jacobian, henon_map, henon_map_inverse
from tanglepack.numerics.budget import PointBudgetExceeded
from tanglepack.numerics.geometry import curvature_area_batch
from tanglepack.numerics.ManifoldMachine import SPACING_FLOOR, ManifoldMachine

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


def _chord_distance(p0, p1, point):
    """How near each chord ``p0 -> p1`` passes ``point``."""
    chord, offset = p1 - p0, point - p0
    along = np.clip(np.sum(offset * chord, axis=1) / np.sum(chord * chord, axis=1), 0.0, 1.0)
    return np.hypot(*(offset - along[:, None] * chord).T)


def test_the_chord_deviation_falls_like_r_squared_near_the_orbit_point():
    workbench, fp = _workbench(R0)
    workbench.grow_n_times(fp, "unstable", num_iterations=4)
    points = workbench.manifolds[(fp, "unstable", 0, 0)].get_point_array()
    p0, p1 = points[1:-2], points[2:-1]  # the pair at the anchor is grade_seed's
    areas = curvature_area_batch(p0, p1, points[:-3], points[3:])
    r = _chord_distance(p0, p1, fp.coordinates[0])
    chords = np.hypot(*(p1 - p0).T)
    near = (r < R0) & (chords >= np.maximum(SPACING_FLOOR * r / R0, _floor(fp)))
    assert near.sum() > 10
    assert np.all(areas[near] < AREA_CUTOFF * (r[near] / R0) ** 2)
    assert np.all(chords[near] <= r[near]), "no chord within r0 longer than its distance"


def test_a_chord_passing_the_orbit_point_between_far_nodes_is_graded():
    # A linear saddle at the origin and a straight curve y = DELTA whose two
    # nodes lie 2 r0 out: both nodes are beyond r0 and the curve is straight,
    # so only the chord's own distance to the orbit point can refine it.
    delta = 1e-6
    stretch = np.array([2.0, 0.5])
    workbench = TangleWorkbench(lambda p: np.asarray(p) * stretch, lambda p: np.asarray(p) / stretch,
                                lambda p: np.diag(stretch), area_cutoff=AREA_CUTOFF, near_orbit_radius=R0)
    fp = workbench.construct_fixed_point([0.0, 0.0])
    store = workbench.point_store
    curve = polyline_manifold([(-2 * R0, delta), (2 * R0, delta)], (fp, "unstable", 0, 0), store=store)
    store.link(store.append([(-R0, 2 * delta), (R0, 2 * delta)], None), curve.ids)
    workbench._man_machine.refine_manifold(curve)

    points = curve.get_point_array()
    p0, p1 = points[:-1], points[1:]
    r, chords = _chord_distance(p0, p1, fp.coordinates[0]), np.hypot(*(p1 - p0).T)
    assert len(points) > 2 and np.all(points[:, 1] == delta)
    assert np.all(chords <= np.maximum(r, R0)), "within r0 no chord longer than its distance"
    assert chords.min() < 2 * delta, "the grading reaches the chord's distance"


def _grown_rows(near_orbit_radius):
    workbench, fp = _workbench(near_orbit_radius)
    workbench.grow_n_times(fp, "unstable", num_iterations=10)  # deep enough that refinement acts
    workbench.grow_n_times(fp, "stable", num_iterations=10)
    store = workbench.point_store
    return [np.asarray(column).tobytes() for column in (store.xy, store.cdist, store.image, store.preimage)]


def test_none_is_the_graded_rule_at_a_vanishing_radius(monkeypatch):
    # The reference: as r0 -> 0 the graded rule IS the plain one (no pair
    # passes within r0, so the area cutoff and spacing floor are today's), and
    # the seed is left ungraded. Two code paths, one refinement.
    plain = _grown_rows(None)
    monkeypatch.setattr(ManifoldMachine, "grade_seed", lambda self, *args: None)
    assert _grown_rows(1e-300) == plain


def test_the_point_budget_still_raises():
    workbench, fp = _workbench(R0, max_points_per_step=200)
    with pytest.raises(PointBudgetExceeded):
        workbench.grow_n_times(fp, "unstable", num_iterations=12)


def test_the_seed_grading_is_budgeted():
    with pytest.raises(PointBudgetExceeded):
        _workbench(R0, max_points_per_step=5)


SEEDS = {
    "period3": (HENON_P3, [list(p) for p in P3_ORBIT_SEED], PERIOD3_ORIENTATION),
    "inversion": (HENON_K4, list(K4_INVERSION_SEED), None),
}


def _seeded(name, near_orbit_radius):
    params, seed, orientation = SEEDS[name]
    workbench = TangleWorkbench(henon_map(*params), henon_map_inverse(*params), henon_jacobian(*params),
                                area_cutoff=AREA_CUTOFF, near_orbit_radius=near_orbit_radius)
    fp = workbench.construct_fixed_point(seed)
    if orientation is not None:
        workbench.orient_eigenvectors(fp, {s: np.array(v) for s, v in orientation.items()})
    workbench.initialize_both_manifolds(fp)
    return workbench, fp


@pytest.mark.parametrize("stability", ["stable", "unstable"])
@pytest.mark.parametrize("name", list(SEEDS))
def test_every_piece_of_the_seed_chain_is_graded_and_linked_along_it(name, stability):
    (workbench, fp), (plain, plain_fp) = _seeded(name, R0), _seeded(name, None)
    assert fp.k_value > 1
    store = workbench.point_store
    # The source of a carried point: its preimage on a stable seed, its image on an unstable one.
    source = store.preimage if stability == "stable" else store.image
    for key in fp.branch_cycle(stability):
        chain, orbit_point = fp.slot_chain(stability, key[2], key[3]), fp.coordinates[key[2]]
        rows = np.asarray(chain.ids)
        offsets = store.xy[rows] - orbit_point
        distances = np.hypot(*offsets.T)
        seed = plain_fp.slot_chain(stability, key[2], key[3])
        seed_ray = plain.point_store.xy[seed.ids[0]] - plain_fp.coordinates[key[2]]
        seed_innermost = np.hypot(*seed_ray)
        assert np.all(np.diff(distances) > 0), key
        assert _floor(fp) <= distances[0] < 1e-6 * seed_innermost, key
        # Inside the seed's innermost point, on its ray (where xy rounding leaves a direction).
        inside = (distances > 1e-9) & (distances < seed_innermost)
        cosines = offsets[inside] @ seed_ray / distances[inside] / seed_innermost
        assert inside.sum() > 5 and np.all(cosines > 1 - 1e-12), f"{key}: the carried points leave the seed's ray"

        # Each carried point comes from the piece one step outward, scaled by the stretch.
        outward = fp.advance_key(key, -1 if stability == "stable" else 1)
        expected = fp.slot_chain(stability, outward[2], outward[3]).index
        linked = rows[(source[rows] >= 0) & (store._chain[np.maximum(source[rows], 0)] >= 0)]
        assert len(linked) > len(seed), key
        assert np.all(store._chain[source[linked]] == expected), key
        np.testing.assert_allclose(store.cdist[source[linked]], store.cdist[linked] * chain.stretch_param,
                                   rtol=1e-9, err_msg=str(key))
