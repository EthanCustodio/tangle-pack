import numpy as np
from tanglepack.FixedPointSolver import FixedPointSolver
from tanglepack.BranchPoint import BranchPoint
from tanglepack.DynamicalSystem import DynamicalSystem


def henon(point):
    """defines the henon map for binary horshoe parameters to test basic functionality"""

    k, b = (10, 1)

    x = point[0]
    y = point[1]

    return np.array([y - k + x ** 2, -b * x])


def henon_inverse(point):
    """defines the inverse henon map for"""

    k, b = (10, 1)

    x = point[0]
    y = point[1]

    return np.array([-y / b, x + k - (y ** 2) / (b ** 2)])


def henon_jacobian(point):
    """defines the jacobian for the henon map"""

    k, b = (10, 1)

    x = point[0]
    y = point[1]

    return np.array([[2*x, 1], [-b, 0]])


def test_initialization():

    system = DynamicalSystem(henon, henon_inverse, henon_jacobian)

    solver = FixedPointSolver(system)

    assert solver.dynamical_map is not None
    assert solver.dynamical_map_inverse is not None
    assert solver.jacobian_function is not None


def test_find_period_one_fixed_point():

    tolerance = 1e-5

    system = DynamicalSystem(henon, henon_inverse, henon_jacobian)

    solver = FixedPointSolver(system)

    initial_guess = [4, -4]

    fixed_point = solver.compute_fixed_point(initial_guess)

    difference = np.abs(solver.multipoint_shoot(fixed_point) - fixed_point)

    assert np.all(np.abs(difference) < tolerance)


def test_find_period_two_fixed_point():

    tolerance = 1e-5

    system = DynamicalSystem(henon, henon_inverse, henon_jacobian)

    solver = FixedPointSolver(system)

    initial_guess = [[4, -4], [4, -4]]

    fixed_point = solver.compute_fixed_point(initial_guess)

    difference = np.abs(solver.multipoint_shoot(fixed_point) - fixed_point)

    assert np.all(np.abs(difference) < tolerance)


def test_compute_single_jacobian():
    """test to check if the jacobian for henon map fixed point can be computed"""

    tolerance = 1e-8

    system = DynamicalSystem(henon, henon_inverse, henon_jacobian)

    solver = FixedPointSolver(system)

    initial_guess = np.array([4, -4])

    fixed_point = solver.compute_fixed_point(initial_guess)

    _, jacobian = solver.compute_jacobian(fixed_point)

    assert np.allclose(jacobian, henon_jacobian(fixed_point[0]), atol=tolerance)


def test_compute_two_jacobians():
    """test to check if a jacobian is cmoputed for each point in the orbit"""

    tolerance = 1e-8

    system = DynamicalSystem(henon, henon_inverse, henon_jacobian)

    solver = FixedPointSolver(system)

    initial_guess = np.array([[4, -4], [4, -4]])

    fixed_point = solver.compute_fixed_point(initial_guess)

    # in this case of a repeated fixed point the overal jacobian is incorrect
    jacobians, _ = solver.compute_jacobian(fixed_point)

    assert np.allclose(jacobians[0], henon_jacobian(fixed_point[0]), atol=tolerance)
    assert np.allclose(jacobians[1], henon_jacobian(fixed_point[0]), atol=tolerance)


def test_compute_eigenvectors():

    tolerance = 1e-6

    system = DynamicalSystem(henon, henon_inverse, henon_jacobian)

    solver = FixedPointSolver(system)

    initial_guess = np.array([4, -4])

    fixed_point = solver.compute_fixed_point(initial_guess)

    eigenvalues, eigenvectors = solver.compute_eigenvectors(fixed_point)

    _, jacobian = solver.compute_jacobian(fixed_point)

    unstable_vector = eigenvectors[0][0]
    stable_vector = eigenvectors[0][1]

    unstable_value = eigenvalues[0][0]
    stable_value = eigenvalues[0][1]

    Jv_unstable = jacobian @ unstable_vector
    Jv_stable = jacobian @ stable_vector

    lv_unstable = unstable_value * unstable_vector
    lv_stable = stable_value * stable_vector

    # checks if J v = lambda v 
    assert np.allclose(Jv_unstable, lv_unstable, atol=tolerance)
    assert np.allclose(Jv_stable, lv_stable, atol=tolerance)

    # checks if the eigenvalues multiply to 1 for the area preserving map
    assert np.isclose(eigenvalues[0][0] * eigenvalues[0][1], 1.0)


def test_build_fixed_point_object():

    system = DynamicalSystem(henon, henon_inverse, henon_jacobian)

    solver = FixedPointSolver(system)
    fixed_point = solver.construct_fixed_point([4, -4], num_branches=1)

    assert all(isinstance(x, BranchPoint) for x in fixed_point.branch_points)
    assert np.allclose(fixed_point.total_jacobian, henon_jacobian(fixed_point.coordinates[0]))

