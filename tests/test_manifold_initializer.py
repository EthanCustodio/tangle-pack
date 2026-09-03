import pytest
import numpy as np

from tanglepack import ManifoldInitializer
from tanglepack import DynamicalSystem
from tanglepack import FixedPointSolver


def henon_map(point):
    """defines the henon map for binary horshoe parameters to test basic functionality"""

    k, b = (10, 1)

    x = point[0]
    y = point[1]

    return np.array([y - k + x ** 2, -b * x])


def henon_map_inverse(point):
    """defines the inverse henon map for"""

    k, b = (10, 1)

    x = point[0]
    y = point[1]

    return np.array([-y / b, x + k - (y ** 2) / (b ** 2)])


def test_initialization_unstable():

    henon = DynamicalSystem(henon_map, henon_map_inverse)

    initial_guess = [4, -4]

    fp_solver = FixedPointSolver(henon)

    fixed_point = fp_solver.construct_fixed_point(initial_guess)

    man_maker = ManifoldInitializer(henon)

    initial_segment = man_maker.get_initial_fundamental_segment(fixed_point, 0, 0, 'unstable')
    initial_points = initial_segment.get_point_array()

    assert len(initial_points) == 3

    assert np.linalg.norm(initial_points[1] - initial_points[0]) < np.linalg.norm(initial_points[2] - initial_points[0])


def test_initialization_stable():

    henon = DynamicalSystem(henon_map, henon_map_inverse)

    initial_guess = [4, -4]

    fp_solver = FixedPointSolver(henon)

    fixed_point = fp_solver.construct_fixed_point(initial_guess)

    man_maker = ManifoldInitializer(henon)

    initial_segment = man_maker.get_initial_fundamental_segment(fixed_point, 0, 0, 'stable')
    initial_points = initial_segment.get_point_array()

    assert len(initial_points) == 3

    assert np.linalg.norm(initial_points[1] - initial_points[0]) < np.linalg.norm(initial_points[2] - initial_points[0])




# --------------------------------------------------------------------------- #
# Plan 1.6 -- ManifoldInitializer called a nonexistent FixedPoint.has_inversion()
# --------------------------------------------------------------------------- #
def _k10_fixed_point():
    """The k=10 binary-horseshoe saddle at [4, -4], oriented, with k_value set."""
    henon = DynamicalSystem(henon_map, henon_map_inverse)
    fp_solver = FixedPointSolver(henon)
    fixed_point = fp_solver.construct_fixed_point([4, -4])
    man_maker = ManifoldInitializer(henon)
    man_maker.orient_manifolds(
        fixed_point,
        {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])},
    )
    return man_maker, fixed_point


def test_two_branch_initialization_grows_opposite_directions():
    """Branch 1 is the eigenvector direction negated, so the two first points
    straddle the fixed point. Previously an AttributeError (``has_inversion``)."""
    man_maker, fixed_point = _k10_fixed_point()

    segments = man_maker.construct_kevin_way(fixed_point, "unstable", num_branches=2)

    assert set(segments) == {(0, 0), (0, 1)}

    origin = np.asarray(fixed_point.coordinates[0]).reshape(-1)
    first_0 = np.asarray(
        man_maker.get_first_point(fixed_point, 0, 0, "unstable")
    ).reshape(-1)
    first_1 = np.asarray(
        man_maker.get_first_point(fixed_point, 0, 1, "unstable")
    ).reshape(-1)

    assert np.allclose(first_0 - origin, -(first_1 - origin))


# --------------------------------------------------------------------------- #
# Plan 1.7 -- construct_kevin_way walks the chain with FixedPoint.advance_key.
# The dict it returns must be unchanged on period 1 and period 3.
# --------------------------------------------------------------------------- #
def _p3_map(point):
    k, b = 2, 1
    x, y = point
    return np.stack([y - k + x**2, -b * x], axis=0)


def _p3_map_inverse(point):
    k, b = 2, 1
    x, y = point
    return np.stack([-y / b, x + k - (y**2) / (b**2)], axis=0)


def _p3_jacobian(point):
    k, b = 2, 1
    x, y = point
    return np.array([[2 * x, 1], [-b, 0]])


def _p3_fixed_point():
    system = DynamicalSystem(_p3_map, _p3_map_inverse, _p3_jacobian)
    fp_solver = FixedPointSolver(system)
    fixed_point = fp_solver.construct_fixed_point([[0, 1], [-1, 0], [-1, 1]])
    man_maker = ManifoldInitializer(system)
    man_maker.orient_manifolds(
        fixed_point,
        {"unstable": np.array([0, -1]), "stable": np.array([-1, -1])},
    )
    return man_maker, fixed_point


@pytest.mark.parametrize("stability", ["unstable", "stable"])
def test_kevin_way_period_one_keys_and_points(stability):
    man_maker, fixed_point = _k10_fixed_point()

    segments = man_maker.construct_kevin_way(fixed_point, stability)

    assert list(segments) == [(0, 0)]
    assert segments[(0, 0)].root is fixed_point.branch_points[0]
    assert len(segments[(0, 0)].get_point_array()) == 3


@pytest.mark.parametrize("stability", ["unstable", "stable"])
def test_kevin_way_period_three_orbit_chain(stability):
    """Chain order: unstable walks 0, 1, 2; stable walks 0, 2, 1 (the inverse map
    grows the stable fundamental segment)."""
    man_maker, fixed_point = _p3_fixed_point()

    segments = man_maker.construct_kevin_way(fixed_point, stability)

    assert set(segments) == {(0, 0), (1, 0), (2, 0)}
    for (orbit_index, _branch_index), manifold in segments.items():
        assert manifold.root is fixed_point.branch_points[orbit_index]
        # every segment sits on the orbit point it is keyed to
        first = manifold.walk_fwd(None, manifold.root, 0)
        assert first is not None
        assert np.linalg.norm(
            np.asarray(first.get_point()).reshape(-1)
            - np.asarray(fixed_point.coordinates[orbit_index]).reshape(-1)
        ) < np.linalg.norm(
            np.asarray(first.get_point()).reshape(-1)
            - np.asarray(
                fixed_point.coordinates[(orbit_index + 1) % fixed_point.period]
            ).reshape(-1)
        )
