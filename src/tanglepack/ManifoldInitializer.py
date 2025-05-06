from typing import Literal
from .FixedPoint import FixedPoint
from .DynamicalSystem import DynamicalSystem
from .ManifoldMachine import ManifoldMachine
from .BaseManifold import BaseManifold
from .Point import Point
import numpy as np


class ManifoldInitializer():


    def __init__(self, system: DynamicalSystem):

        self.system = system


    def get_first_point(self, fixed_point: FixedPoint, orbit_index, branch_index, stability: Literal["stable", "unstable"]):
        """
        Computes the first point from the fixed point based on a linear interpolation
        """

        step = fixed_point.accuracy

        if stability == 'unstable':
            direction_from_fixed_point = fixed_point.unstable_eigenvectors[orbit_index]
        else:
            direction_from_fixed_point = fixed_point.stable_eigenvectors[orbit_index]

        direction_from_fixed_point = np.asarray(direction_from_fixed_point).flatten()

        first_point = fixed_point.coordinates[orbit_index] + (step) * direction_from_fixed_point

        return np.array(first_point, dtype=np.float64).reshape(-1)
    

    def get_first_point_preiterate(self, fixed_point: FixedPoint, orbit_index, branch_index, stability: Literal["stable", "unstable"]):
        """
        Get the iterate of the first point
        """

        first_point = self.get_first_point(fixed_point, orbit_index, branch_index, stability)

        if stability == "unstable":
            first_preiterate = self.system.map_inv(first_point)
        
        else:  # stable branch
            first_preiterate = self.system.map(first_point)

        return first_preiterate
    

    def get_initial_fundamental_segment(self, fixed_point: FixedPoint, orbit_index, branch_index, stability: Literal["stable", "unstable"]):
        """
        Computes the initial fundamental segment from iterating the first point
        """

        first_point = self.get_first_point(fixed_point, orbit_index, branch_index, stability)
        distance_first = np.linalg.norm(first_point - fixed_point.coordinates[orbit_index])

        first_preiterate = self.get_first_point_preiterate(fixed_point, orbit_index, branch_index, stability)
        distance_prev = np.linalg.norm(first_preiterate - fixed_point.coordinates[orbit_index])

        alpha = distance_first / distance_prev

        first_point = Point(first_point[0], first_point[1], cdist=distance_first, edist=distance_first, stretch_param=alpha)
        first_preiterate = Point(first_preiterate[0], first_preiterate[1], cdist=distance_prev, edist=distance_prev, stretch_param=alpha)

        if stability == "unstable":

            fixed_point.branch_points[orbit_index].insert_point_forward(first_preiterate, branch_index)
            fixed_point.branch_points[orbit_index].insert_point_forward(first_point, branch_index)

        else:  # stable

            fixed_point.branch_points[orbit_index].insert_point_backward(first_preiterate, branch_index)
            fixed_point.branch_points[orbit_index].insert_point_backward(first_point, branch_index)

        return BaseManifold(fixed_point.branch_points[orbit_index], stability, alpha, tail=first_point)

