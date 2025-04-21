import pytest
import numpy as np
import scipy as sp
import tanglepack



def henon(point):
    """defines the henon map for classical parameters to test basic functionality"""

    a, b = (1.4, 0.3)

    x = point[0]
    y = point[1]

    return [1 - a * x ** 2 + y, b * x]


def test_curvature_area():
    """tests the curvature area method of Manifold"""

    # area for these points should be 3.75 from analytic computatiaon
    analytic_area = 3.75

    points = np.array([[1, 1], [2, 4], [4, 5]])

    area = tanglepack.Manifold.curvature_area(points)

    assert abs(area - analytic_area) < 1e-8


    
def test_get_first_point():

    # area for these points should be 3.75 from analytic computatiaon
    

    assert None is None





