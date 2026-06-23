import numpy as np

from tanglepack import ManifoldMachine
from tanglepack import DynamicalSystem


def henon_map(point):
    """Hénon map, binary-horseshoe parameters."""
    k, b = (10, 1)
    x, y = point[0], point[1]
    return np.array([y - k + x**2, -b * x])


def henon_map_inverse(point):
    """Inverse Hénon map."""
    k, b = (10, 1)
    x, y = point[0], point[1]
    return np.array([-y / b, x + k - (y**2) / (b**2)])


def test_machine_initialization():
    system = DynamicalSystem(henon_map, henon_map_inverse)
    machine = ManifoldMachine(system)
    assert machine.system is system
