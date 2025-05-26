from tanglepack.Point import Point
import numpy as np
from collections import deque
import scipy.integrate as spi
from .DynamicalSystem import DynamicalSystem
from .BaseManifold import BaseManifold
from .ManifoldView import ManifoldView
from .BranchPoint import BranchPoint
from numpy.linalg import LinAlgError
import time

class ManifoldMachine():


    def __init__(self, system: DynamicalSystem):

        self.system = system
        self.area_cutoff = 1e-2

    
    def grow_manifold(self, manifold: BaseManifold):

        # this line meant to deal with roots that are fixed points
        # if it is an intersection behavior is unknown
        temp_root = manifold.root
        if isinstance(manifold.root, BranchPoint):
            manifold.root = manifold.walk_fwd(None, temp_root)

        # print(f"current manifold: {manifold.get_point_array()}")

        iterated_manifold = self.iterate_manifold(manifold)

        grown_manifold = self.merge_manifolds(manifold, iterated_manifold)

        grown_manifold.root = temp_root

        return grown_manifold


    def iterate_manifold(self, manifold: BaseManifold):
        """
        Iterates all uniterated points in a manifold
        """
        from .ManifoldInitializer import ManifoldInitializer
        initalizer = ManifoldInitializer(self.system)
        viewer = ManifoldView(manifold, self.system)

        # print(f'Iterating these points: \n{manifold.get_point_array()}')

        non_iterated_coords = manifold.get_non_iterated_point_array()
        non_iterated_cdists = manifold.get_non_iterated_cdist_array()
        non_iterated_points = manifold.get_non_iterated_point_array(return_nodes=True)

        if len(non_iterated_coords):
            # iterated_points = viewer.map_fwd(non_iterated_coords)
            iterated_points = np.vstack(
                [viewer.map_fwd(p) for p in non_iterated_coords]
            )
            distances = manifold.stretch_param * non_iterated_cdists

            xvals = iterated_points[:, 0]
            yvals = iterated_points[:, 1]

            new_points = [
                Point(x, y, cdist, stretch_param=manifold.stretch_param) 
                for x, y, cdist in zip(xvals, yvals, distances)
                ]
            
            for i, point in enumerate(non_iterated_points):
                point.insert_next_iterate(new_points[i])

            new_iterated_points = initalizer.construct_manifold_from_point_list(
                new_points, manifold.stability, manifold.stretch_param, manifold.branch_index
                )
                    
        old_points = manifold.get_iterated_point_array(return_nodes=True)

        old_iterated_points = BaseManifold(
            old_points[0], manifold.stability, manifold.stretch_param, 
            tail=old_points[-1], branch_index=manifold.branch_index
            )

        # old_iterated_points = initalizer.construct_manifold_from_point_list(
        #     old_points, manifold.stability, manifold.stretch_param, manifold.branch_index
        # )

        # if there were no points that needed to be mapped
        if not len(non_iterated_coords) == 0:
            # print(f"_____ \n old iterated: \n {old_iterated_points.get_point_array()}")
            # print(f"new iterated: \n {new_iterated_points.get_point_array()}")

            mapped_manifold = self.merge_manifolds(old_iterated_points, new_iterated_points)
            self.refine_manifold(mapped_manifold)

            # print(f'merged iterated: \n{mapped_manifold.get_point_array()} \n _________')
            assert (sorted(mapped_manifold.get_cdist_array()) == mapped_manifold.get_cdist_array()).all()

            return mapped_manifold
        
        else:
            # print(f'iterated manifold: \n{old_iterated_points.get_point_array()}')
            return old_iterated_points


    def merge_manifolds(self, manifold_1: BaseManifold, manifold_2: BaseManifold):

        head_1 = manifold_1.root
        head_2 = manifold_2.root

        # If there are negative cdists this merge method will fail
        if head_1.cdist <= head_2.cdist:
            current_point = head_1
            head_1 = manifold_1.walk_fwd(None, head_1)
            output_manifold = manifold_1
        else:
            current_point = head_2
            head_2 = manifold_2.walk_fwd(None, head_2)
            output_manifold = manifold_2

        # the termination point right after the manifold tail
        over_one_1 = manifold_1.walk_fwd(None, manifold_1.tail)
        over_one_2 = manifold_2.walk_fwd(None, manifold_2.tail)

        while head_1 is not over_one_1 and head_2 is not over_one_2:
            
            if head_1.cdist < head_2.cdist:
                self._insert_point_geometrically(current_point, head_1, manifold_1)
                head_1 = manifold_1.walk_fwd(None, head_1)
            
            elif head_2.cdist < head_1.cdist:
                self._insert_point_geometrically(current_point, head_2, manifold_2)
                head_2 = manifold_2.walk_fwd(None, head_2)

            # they are the same node
            else:
                if head_1 is head_2:
                    self._insert_point_geometrically(current_point, head_1, manifold_1)
                    head_1 = manifold_1.walk_fwd(None, head_1)
                    head_2 = manifold_2.walk_fwd(None, head_2)
            
            current_point = manifold_1.walk_fwd(None, current_point)

        # manifold_2 emptied first
        if head_1 is not over_one_1:
            self._insert_point_geometrically(current_point, head_1, manifold_1)
            output_manifold.tail = manifold_1.tail
        else:
            self._insert_point_geometrically(current_point, head_2, manifold_2)
            output_manifold.tail = manifold_2.tail


        return output_manifold
        

    # def merge_manifolds(self, manifold_1: BaseManifold, manifold_2: BaseManifold):
    #     """
    #     Merges the second manifold into the first manifold
    #     """

    #     prev_point_1 = manifold_1.root
    #     curr_point_1 = manifold_1.walk_fwd(None, prev_point_1)

    #     prev_point_2 = None
    #     curr_point_2 = manifold_2.root

    #     while curr_point_1 is not manifold_1.tail and curr_point_2 is not manifold_2.tail:

    #         next_point_2 = manifold_2.walk_fwd(prev_point_2, curr_point_2)

    #         # if you encounter the same point
    #         # the rest of the strip must be the same
    #         if curr_point_1 is curr_point_2:
                
    #             prev_point_2, curr_point_2 = curr_point_2, next_point_2
    #             # if manifold_2.tail.cdist > manifold_1.tail.cdist:
    #             #     manifold_1.tail = manifold_2.tail

    #             # return manifold_1

    #         elif prev_point_1.cdist < curr_point_2.cdist < curr_point_1.cdist:
                
    #             temp_point_2 = curr_point_2

    #             prev_point_2, curr_point_2 = None, next_point_2

    #             self._insert_point_geometrically(
    #                 prev_point_1, temp_point_2, manifold_1, manifold_1.branch_index)
                    
    #             # next_point_1 = manifold_1.walk_fwd(temp_point_2, curr_point_1)
    #             # prev_point_1, curr_point_1 = temp_point_2, next_point_1
    #             next_point_1 = manifold_1.walk_fwd(prev_point_1, temp_point_2)
    #             prev_point_1 = temp_point_2
    #             curr_point_1 = next_point_1
    #             assert prev_point_1.cdist < temp_point_2.cdist, "merge order broke"

    #             # next_point_1 = manifold_1.walk_fwd(prev_point_1, curr_point_1)
    #             # prev_point_1, curr_point_1 = curr_point_1, next_point_1
                

    #         else:
    #             next_point_1 = manifold_1.walk_fwd(prev_point_1, curr_point_1)
    #             prev_point_1, curr_point_1 = curr_point_1, next_point_1
                
        
    #     # if we reach the end of the first manifold first
    #     # append just attach the curr_point_2 to the end and update the tail
    #     if curr_point_1 is manifold_1.tail:
            
    #         if curr_point_1 is not curr_point_2:
    #             self._insert_point_geometrically(
    #                 curr_point_1, curr_point_2, manifold_1, manifold_1.branch_index
    #             )
    #         manifold_1.tail = manifold_2.tail

    #         return manifold_1
        
    #     # if we reach the end of the second manifold first
    #     # add in the last node and do not update the tail
    #     if curr_point_2 is manifold_2.tail:

    #         while curr_point_1 is not None:

    #             if prev_point_1.cdist < curr_point_2.cdist < curr_point_1.cdist:
    #                 self._insert_point_geometrically(
    #                     prev_point_1, curr_point_2, manifold_1, manifold_1.branch_index
    #                     )
    #                 return manifold_1

    #             next_point_1 = manifold_1.walk_fwd(prev_point_1, curr_point_1)
    #             prev_point_1, curr_point_1 = curr_point_1, next_point_1
    #             if prev_point_1 is manifold_1.tail:
    #                 break

    #         if prev_point_1 is not curr_point_2:
    #             self._insert_point_geometrically(
    #                 prev_point_1, curr_point_2, manifold_1, manifold_1.branch_index
    #                 )

    #         if curr_point_1 is None:
    #             manifold_1.tail = curr_point_2
            
    #         return manifold_1


    def refine_manifold(self, manifold: BaseManifold, branch_index=None, final_node=None):
        """
        Adds additional points in areas of the manifold with high curvature

        Checks every consecutive set of three points in the manifold.
        Performs a linear and a parabolic fit between them.
        If the area bounded by those curves is less than self.area_cutoff
        then add additional points.
        Iterates through the manifold until max_passes is reached.        
        """
        # print(f"refining this manifold: \n {manifold.get_point_array()}")
        # print(f"cdists {manifold.get_cdist_array()}")
        final_node = manifold.tail

        num_initial_points = len(manifold.get_point_array())

        previous_point = manifold.root
        current_point = manifold.walk_fwd(None, previous_point)

        while current_point is not None:

            # print(f'Checking these points: {(previous_point.get_point(), current_point.get_point())}')
            self.refine_two_points((previous_point, current_point), manifold, branch_index)

            if current_point is final_node:
                break

            next_point = manifold.walk_fwd(previous_point, current_point)

            previous_point, current_point = current_point, next_point


        num_final_points = len(manifold.get_point_array())

        print(f"{num_final_points - num_initial_points} POINTS ADDED DURING REFINEMENT")
        print(f'{num_final_points} Number of Current Points')


    def refine_two_points(self, points: tuple[Point, Point], manifold: BaseManifold, branch_index=None):
        """
        Iteratively refines the two given points.
        Uses a modified stack structure to iterate through the manifold 
        and continuously refine the points that it creates as needed.
        
        """
        viewer = ManifoldView(manifold, self.system)
        
        left_point, right_point = points

        pair_queue = deque([(left_point, right_point)])

        while pair_queue:

            p0, p1 = pair_queue.pop()

            # check if points are so close together to cause numerical instability
            if np.linalg.norm(p1.get_point() - p0.get_point()) < 1e-8:
                continue

            p2 = manifold.walk_fwd(p0, p1)

            # if we are at the last two points in a manifold bail 
            if p2 is None:
                break

            triplet = np.vstack((p0.get_point(), p1.get_point(), p2.get_point()))

            try:
                curvature_area = self.curvature_area(triplet)

            except LinAlgError: # singular Vandermond matrix
                continue

            if abs(curvature_area) < self.area_cutoff:
                continue

            p0_preiterate = self._get_preiterate(p0, manifold.stability)
            p1_preiterate = self._get_preiterate(p1, manifold.stability)

            new_point_coords_back = 0.5 * (p1_preiterate.get_point() + p0_preiterate.get_point())

            new_distance = 0.5 * (p0.cdist + p1.cdist)

            new_point_coords = viewer.map_fwd(new_point_coords_back)

            x = new_point_coords[0]
            y = new_point_coords[1]
            new_point = Point(x, y, new_distance, stretch_param=p0.stretch_param)

            self._cache_preiterate(new_point, new_point_coords_back, manifold.stability)

            self._insert_point_geometrically(p0, new_point, manifold, branch_index)

            pair_queue.append((p0, new_point))
            pair_queue.append((new_point, p1))
            pair_queue.append((p1, p2))


    def _get_preiterate(self, point: Point, stability: str):
        """helper function to get preiterate based on stability"""

        if stability == "unstable":
            return point.prev_iterate
        else:
            return point.next_iterate


    def _cache_preiterate(self, point: Point, preiterate_coords, stability: str):
        """helper function to cache preiterate based on stability"""

        # phoney cached preiterate point without cdist
        cached_preiterate = Point(preiterate_coords[0], preiterate_coords[1])

        if stability == "unstable":
            point.prev_iterate = cached_preiterate
        else:
            point.next_iterate = cached_preiterate


    def _insert_point_geometrically(self, p0: Point, new_point: Point, manifold: BaseManifold, branch_index=None):
        """helper function to insert points smartly
            based on stability and brach_point'ness"""

        if isinstance(p0, BranchPoint):
            if manifold.stability == "unstable":
                if p0.forward_branches[branch_index] is new_point:
                    return
                p0.insert_point_forward(new_point, branch_index=branch_index)
            else:
                if p0.backward_branches[branch_index] is new_point:
                    return
                p0.insert_point_backward(new_point, branch_index=branch_index)

        else:
            if manifold.stability == "unstable":
                if p0.forward is new_point:
                    return
                p0.insert_point_forward(new_point)
            else:
                if p0.backward is new_point:
                    return
                p0.insert_point_backward(new_point)


    @staticmethod
    def linear_fit(points):
        """
        Takes in three points and gives the linear fit between the first and last
        
        Parameters:
            points: list of three points
        """

        point_one = points[0]
        point_two = points[-1]

        m = (point_two[1] - point_one[1]) / (point_two[0] - point_one[0])
        b = point_one[1] - m * point_one[0]

        return np.poly1d([m, b])


    @staticmethod
    def parabolic_fit(points):
        """
        Takes in three points and gives the parabolic fit between them
        
        Parameters:
            points: list of three points
        """

        x_vals = points[:, 0]
        y_vals = points[:, 1]

        A = np.vstack([x_vals**2, x_vals, np.ones_like(x_vals)]).T
        # Solve for the coefficients [a, b, c].
        coefficient = np.linalg.solve(A, y_vals)

        poly = np.poly1d(coefficient)

        return poly


    @staticmethod
    def curvature_area(points):
            """
            Takes in three points and gives
            the area of the shape made by a parabolic and linear fit
            
            Parameters:
                points: list of three points
            """

            x_min = points[0, 0]
            x_max = points[-1, 0]

            parabola = ManifoldMachine.parabolic_fit(points)
            line = ManifoldMachine.linear_fit(points)

            difference = lambda x: np.abs(parabola(x) - line(x))

            area, _ = spi.quad(difference, x_min, x_max)

            return area




