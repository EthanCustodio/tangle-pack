import pytest
import numpy as np
import tanglepack



def henon(point):
    """defines the henon map for classical parameters to test basic functionality"""

    a, b = (1.4, 0.3)

    x = point[0]
    y = point[1]

    return [1 - a * x ** 2 + y, b * x]


def test_compute_fixed_point():
    """test to check if the fixed point for henon map can be computed"""

    tolerance = 1e-8

    guess = np.array([0.6, 0.2])
    # guess = np.array([[0.6, 0.2], [0.6, 0.2]])

    point = tanglepack.FixedPoint(henon, guess)

    fixed_point = point.fixed_point

    difference = point.accuracy

    print(f'\n Fixed Point: \n {fixed_point}')
    print(f'\n Difference: \n {difference}')

    assert np.allclose(difference, 0.0, atol=tolerance)


def test_compute_jacobian():
    """test to check if the jacobian for henon map fixed point can be computed"""

    tolerance = 1e-8

    guess = np.array([0.6, 0.2])
    # guess = np.array([[0.6, 0.2], [0.6, 0.2]])


    point = tanglepack.FixedPoint(henon, guess)

    fixed_point = point.fixed_point

    difference = point.accuracy

    point.compute_jacobian()

    assert np.allclose(difference, 0.0, atol=tolerance)


def test_compute_eigenvectors():
    """test to check if the eigenvectors for henon map fixed point can be computed"""

    tolerance = 1e-8

    guess = np.array([0.6, 0.2])
    # guess = np.array([[0.6, 0.2], [0.6, 0.2]])


    point = tanglepack.FixedPoint(henon, guess)

    fixed_point = point.fixed_point

    difference = point.accuracy

    point.compute_eigenvectors()

    print(f'eigens {point.eigenvectors}')

    assert np.allclose(difference, 0.0, atol=tolerance)

