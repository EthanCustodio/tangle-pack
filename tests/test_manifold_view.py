import pytest
import numpy as np
from tanglepack.ManifoldView import ManifoldView
from tanglepack.DynamicalSystem import DynamicalSystem
from tanglepack.BaseManifold import BaseManifold
from tanglepack.Point import Point


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


def test_initialize_unstable():

    point = Point()

    manifold = BaseManifold(point, stability='unstable', stretch_param=1.0)
    system = DynamicalSystem(henon, henon_inverse)

    viewer = ManifoldView(manifold, system)

    assert viewer.stability == 'unstable'
    assert (viewer.map_fwd, viewer.map_back) == (system.map, system.map_inv)


def test_initialize_stable():

    point = Point()

    manifold = BaseManifold(point, stability='stable', stretch_param=1.0)
    system = DynamicalSystem(henon, henon_inverse)

    viewer = ManifoldView(manifold, system)

    assert viewer.stability == 'stable'
    assert (viewer.map_fwd, viewer.map_back) == (system.map_inv, system.map)


