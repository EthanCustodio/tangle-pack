from typing import Literal
from .FixedPoint import FixedPoint
from .BranchPoint import BranchPoint
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
    

    def get_first_point_back(self, fixed_point: FixedPoint, orbit_index, branch_index, stability: Literal["stable", "unstable"]):
        """
        Get the iterate of the first point
        """

        first_point = self.get_first_point(fixed_point, orbit_index, branch_index, stability)

        if stability == "unstable":
            first_back = self.system.map_inv(first_point)
        
        else:  # stable branch
            first_back = self.system.map(first_point)

        return first_back
    

    def get_initial_fundamental_segment(self, fixed_point: FixedPoint, orbit_index, branch_index, stability: Literal["stable", "unstable"]):
        """
        Computes the initial fundamental segment from iterating the first point
        """

        first_point = self.get_first_point(fixed_point, orbit_index, branch_index, stability)
        distance_first = np.linalg.norm(first_point - fixed_point.coordinates[orbit_index])

        first_back = self.get_first_point_back(fixed_point, orbit_index, branch_index, stability)
        distance_prev = np.linalg.norm(first_back - fixed_point.coordinates[orbit_index])

        alpha = distance_first / distance_prev

        first_point = Point(first_point[0], first_point[1], cdist=distance_first, edist=distance_first, stretch_param=alpha)
        first_back = Point(first_back[0], first_back[1], cdist=distance_prev, edist=distance_prev, stretch_param=alpha)

        first_back.insert_next_iterate(first_point)

        if stability == "unstable":

            fixed_point.branch_points[orbit_index].insert_point_forward(first_point, branch_index)
            fixed_point.branch_points[orbit_index].insert_point_forward(first_back, branch_index)

        else:  # stable

            fixed_point.branch_points[orbit_index].insert_point_backward(first_point, branch_index)
            fixed_point.branch_points[orbit_index].insert_point_backward(first_back, branch_index)

        return BaseManifold(fixed_point.branch_points[orbit_index], stability, alpha, tail=first_point, branch_index=branch_index)


    def construct_manifold_from_point_list(self, points: list[Point], stability: Literal["stable", "unstable"], stretch_param, branch_index=None):
        """
        Constructs a manifold from a list of Point objects
        The given list is assumed to be in cdist ordering
        """

        manifold = BaseManifold(points[0], stability, stretch_param)

        current_point = manifold.root

        for i, point in enumerate(points):

            # start inserting the second point in the list
            if i == 0:
                continue

            self._insert_point_geometrically(current_point, point, manifold, branch_index)
                
            current_point = point

        manifold.tail = points[-1]
        return manifold


    def _insert_point_geometrically(self, p0: Point, new_point: Point, manifold: BaseManifold, branch_index=None):
        """helper function to insert points smartly
            based on stability and brach_point'ness"""

        if isinstance(p0, BranchPoint):
            if manifold.stability == "unstable":
                p0.insert_point_forward(new_point, branch_index=branch_index)
            else:
                p0.insert_point_backward(new_point, branch_index=branch_index)

        else:
            if manifold.stability == "unstable":
                p0.insert_point_forward(new_point)
            else:
                p0.insert_point_backward(new_point)