import pytest
import numpy as np

from tanglepack.ManifoldInitializer import ManifoldInitializer
from tanglepack.DynamicalSystem import DynamicalSystem
from tanglepack.FixedPointSolver import FixedPointSolver


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

    fixed_point = fp_solver.construct_fixed_point(initial_guess, 1)

    man_maker = ManifoldInitializer(henon)

    initial_segment = man_maker.get_initial_fundamental_segment(fixed_point, 0, 0, 'unstable')
    initial_points = initial_segment.get_point_array(branch_index=0)

    assert len(initial_points) == 3

    assert np.linalg.norm(initial_points[1] - initial_points[0]) < np.linalg.norm(initial_points[2] - initial_points[0])


def test_initialization_stable():

    henon = DynamicalSystem(henon_map, henon_map_inverse)

    initial_guess = [4, -4]

    fp_solver = FixedPointSolver(henon)

    fixed_point = fp_solver.construct_fixed_point(initial_guess, 1)

    man_maker = ManifoldInitializer(henon)

    initial_segment = man_maker.get_initial_fundamental_segment(fixed_point, 0, 0, 'stable')
    initial_points = initial_segment.get_point_array(branch_index=0)

    assert len(initial_points) == 3

    assert np.linalg.norm(initial_points[1] - initial_points[0]) < np.linalg.norm(initial_points[2] - initial_points[0])


