import pytest
import numpy as np
from tanglepack.ManifoldMachine import ManifoldMachine
from tanglepack.DynamicalSystem import DynamicalSystem


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


def test_curvature_area():
    """tests the curvature area method of Manifold"""

    # area for these points should be 3.75 from analytic computatiaon
    analytic_area = 3.75

    points = np.array([[1, 1], [2, 4], [4, 5]])

    area = ManifoldMachine.curvature_area(points)

    assert abs(area - analytic_area) < 1e-8


def test_linear_fit_simple():
    """tests if the fit betwee nthe outer two points has the right center"""

    # Line through (0,0) and (2,2) should be y = x
    points = np.array([[0, 0], [1, 1], [2, 2]])
    line = ManifoldMachine.linear_fit(points)

    assert np.isclose(line(0), 0)
    assert np.isclose(line(1), 1)
    assert np.isclose(line(2), 2)


def test_linear_fit_ignores_middle_point():
    """tests if it only fits to outer two points"""

    points_a = np.array([[0, 0], [1, 999], [2, 2]])  # middle is garbage
    points_b = np.array([[0, 0], [1, 1], [2, 2]])

    line_a = ManifoldMachine.linear_fit(points_a)
    line_b = ManifoldMachine.linear_fit(points_b)

    # Coefficients should match
    assert np.allclose(line_a.coefficients, line_b.coefficients)


def test_parabolic_fit_quadratic():

    # y = x^2 through three points
    points = np.array([[-1, 1], [0, 0], [1, 1]])
    poly = ManifoldMachine.parabolic_fit(points)

    # Should exactly interpolate the input points
    for x, y in points:
        assert np.isclose(poly(x), y)

    # Should be close to y = x^2
    assert np.allclose(poly.coefficients, [1, 0, 0], atol=1e-8)


def test_parabolic_fit_general():
    
    # Fit y = 2x^2 + 3x + 1
    x = np.array([-1, 0, 1])
    y = 2*x**2 + 3*x + 1
    points = np.column_stack([x, y])

    poly = ManifoldMachine.parabolic_fit(points)

    for x_val, y_val in points:
        assert np.isclose(poly(x_val), y_val)


def test_machine_initialization():

    system = DynamicalSystem(henon_map, henon_map_inverse)

    machine = ManifoldMachine(system)

    assert machine.system is system
