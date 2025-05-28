from collections import deque
from typing import Literal
import logging

import numpy as np
import scipy.integrate as spi
from numpy.linalg import LinAlgError

from .Point import Point
from .DynamicalSystem import DynamicalSystem
from .BaseManifold import BaseManifold
from .ManifoldView import ManifoldView
from .BranchPoint import BranchPoint

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class ManifoldMachine:

    def __init__(self, system: DynamicalSystem):

        self.system = system
        self.area_cutoff = 1e-2

    def grow_manifold(self, manifold: BaseManifold):
        """
        Takes a manifold and iterates all uniterated points
        then adds those to the original manifold and returns
        the result.

        Parameters:
            manifold: manifold to grow
        """

        # this line meant to deal with roots that are fixed points
        # if it is an intersection behavior is unknown
        temp_root = manifold.root
        if isinstance(manifold.root, BranchPoint):
            manifold.root = manifold.walk_fwd(None, temp_root)

        iterated_manifold = self.iterate_manifold(manifold)

        grown_manifold = self.merge_manifolds(manifold, iterated_manifold)

        grown_manifold.root = temp_root

        return grown_manifold

    def iterate_manifold(self, manifold: BaseManifold):
        """
        Iterates all uniterated points in a manifold
        Returns a new manifold pointing only to the iterate

        Parameters:
            manifold: manifold to iterate
        """
        from .ManifoldInitializer import ManifoldInitializer

        initalizer = ManifoldInitializer(self.system)
        viewer = ManifoldView(manifold, self.system)

        non_iterated_coords = manifold.get_non_iterated_point_array()
        non_iterated_cdists = manifold.get_non_iterated_cdist_array()
        non_iterated_points = manifold.get_non_iterated_point_array(return_nodes=True)

        if len(non_iterated_coords):

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
                if manifold.stability == "unstable":
                    point.insert_next_iterate(new_points[i])
                else:
                    point.insert_prev_iterate(new_points[i])

            new_iterated_points = initalizer.construct_manifold_from_point_list(
                new_points,
                manifold.stability,
                manifold.stretch_param,
                manifold.branch_index,
            )

        old_points = manifold.get_iterated_point_array(return_nodes=True)

        old_iterated_points = BaseManifold(
            old_points[0],
            manifold.stability,
            manifold.stretch_param,
            tail=old_points[-1],
            branch_index=manifold.branch_index,
        )

        # if there were no points that needed to be mapped
        if not len(non_iterated_coords) == 0:

            mapped_manifold = self.merge_manifolds(
                old_iterated_points, new_iterated_points
            )
            self.refine_manifold(mapped_manifold)

            assert (
                sorted(mapped_manifold.get_cdist_array())
                == mapped_manifold.get_cdist_array()
            ).all()

            return mapped_manifold

        else:
            return old_iterated_points

    def merge_manifolds(self, manifold_1: BaseManifold, manifold_2: BaseManifold):
        """
        O(n) inplace method for merging two linked lists (manifold_1 and manifold_2)
        Merges two manifolds based on their cdists.
        Does not return a 3rd manifold, returns the whichever manifold has the smaller
        starting cdist with the new tail set

        Parameters:
            manifold_1: manifold to merge
            manifold_2: manifold to merge
        """

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

    def refine_manifold(
        self, manifold: BaseManifold, branch_index=None, final_node=None
    ):
        """
        Adds additional points in areas of the manifold with high curvature

        Checks every consecutive set of three points in the manifold.
        Performs a linear and a parabolic fit between them.
        If the area bounded by those curves is less than self.area_cutoff
        then add additional points.
        Iterates through the manifold until max_passes is reached.

        Parameters:
            manifold (BaseManifold): current manifold
            branch_index: Optional if starting at a fixed point
            final_node: Optional otherwise set as manifold.tail

        Refactor:
            branch_index and final_node can probably be eliminated from the input
        """
        logger.debug("Refining manifold: %r", manifold.get_point_array())

        final_node = manifold.tail

        num_initial_points = len(manifold.get_point_array())

        previous_point = manifold.root
        current_point = manifold.walk_fwd(None, previous_point)

        while current_point is not None:

            logger.debug(
                "Checking these points:",
                (previous_point.get_point(), current_point.get_point()),
            )

            self.refine_two_points(
                (previous_point, current_point), manifold, branch_index
            )

            if current_point is final_node:
                break

            next_point = manifold.walk_fwd(previous_point, current_point)

            previous_point, current_point = current_point, next_point

        num_final_points = len(manifold.get_point_array())

        logger.info(
            "%d Points added during refinement", num_final_points - num_initial_points
        )
        logger.info("%d NUMBER OF CURRENT POINTS", num_final_points)

    def refine_two_points(
        self,
        points: tuple[Point | BranchPoint, Point | BranchPoint],
        manifold: BaseManifold,
        branch_index=None,
    ):
        """
        Iteratively refines the two given points.
        Uses a modified stack structure to iterate through the manifold
        and continuously refine the points that it creates as needed.

        Parameters:
            points (tuple): set of two points to check for refinement
            manifold (BaseManifold): manifold the two points are on
            branch_index: Optional input if starting at a fixed point

        Refactor:
            branch_index is probably unnecessary since it was
            put into manifold
            adaptively change area_cuttoff refine as the manifolds increase in length
                or we are at smaller scales
        """
        viewer = ManifoldView(manifold, self.system)

        left_point, right_point = points
        pair_queue = deque([(left_point, right_point)])

        while pair_queue:

            p0, p1 = pair_queue.pop()

            # Check if points are so close together to cause numerical instability
            if np.linalg.norm(p1.get_point() - p0.get_point()) < 1e-8:
                continue

            p2 = manifold.walk_fwd(p0, p1)

            # If we are at the last two points in a manifold bail
            if p2 is None:
                break

            triplet = np.vstack((p0.get_point(), p1.get_point(), p2.get_point()))
            try:
                curvature_area = self._curvature_area(triplet)
            except LinAlgError:
                logger.debug(
                    "Singular Vandermond matrix at cdist %.3g–%.3g; skipping. "
                    "Points are likely near vertical",
                    p0.cdist,
                    p1.cdist,
                )
                continue

            if abs(curvature_area) < self.area_cutoff:
                continue

            new_point = self._get_refined_point(p0, p1, viewer, manifold.stability)

            self._insert_point_geometrically(p0, new_point, manifold, branch_index)

            pair_queue.append((p0, new_point))
            pair_queue.append((new_point, p1))
            pair_queue.append((p1, p2))

    def _get_refined_point(
        self,
        p0: Point | BranchPoint,
        p1: Point | BranchPoint,
        viewer: ManifoldView,
        stability: Literal["unstable", "stable"],
    ):
        """ """
        p0_preiterate = self._get_preiterate(p0, stability)
        p1_preiterate = self._get_preiterate(p1, stability)

        new_point_coords_back = 0.5 * (
            p1_preiterate.get_point() + p0_preiterate.get_point()
        )

        new_distance = 0.5 * (p0.cdist + p1.cdist)

        new_point_coords = viewer.map_fwd(new_point_coords_back)

        x = new_point_coords[0]
        y = new_point_coords[1]
        new_point = Point(x, y, new_distance, stretch_param=p0.stretch_param)

        self._cache_preiterate(new_point, new_point_coords_back, stability)

        return new_point

    @staticmethod
    def _get_preiterate(
        point: Point | BranchPoint, stability: Literal["unstable", "stable"]
    ):
        """
        Helper function to get a point's preiterate based on stability

        Parameters:
            point (Point or BranchPoint): point to get preiterate
            stability (str): stability of the manifold
        """

        if stability == "unstable":
            return point.prev_iterate
        else:
            return point.next_iterate

    @staticmethod
    def _cache_preiterate(
        point: Point | BranchPoint,
        preiterate_coords,
        stability: Literal["unstable", "stable"],
    ):
        """
        Helper function to cache preiterate based on stability.
        Only stores the coordinates to be accessed for refinement.

        Parameters:
            point (Point or BranchPoint): Current Points
            preiterate_coords: Preiterate coordinates to cache
            stability (str): stability of the manifold
        """

        # phoney cached preiterate point without cdist
        cached_preiterate = Point(preiterate_coords[0], preiterate_coords[1])

        if stability == "unstable":
            point.prev_iterate = cached_preiterate
        else:
            point.next_iterate = cached_preiterate

    @staticmethod
    def _insert_point_geometrically(
        p0: Point | BranchPoint,
        new_point: Point | BranchPoint,
        manifold: BaseManifold,
        branch_index: int = None,
    ):
        """
        Takes in a point (p0) on the manifold and inserts new_point after it.
        In this case we are inserting away from the fixed point. On the
        unstable manifold that is 'forward'. On the stable manifold that
        is 'backward'.

        Parameters:
            p0 (Point): Anchor point we wil insert after
            new_point (Point): Point to insert into the manifold
            manifold (BaseManifold): Current manifold
            branch_index (int): Optional parameter if inserting after fixed point

        Note:
            Use this method to smoothly handle insertion regardless of stability
            and BranchPoints

        REFACTOR:
            branch_index may not be necessary if that information is already in manifold
        """

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
    def _linear_fit(points):
        """
        Takes in three points and gives the linear fit between the first and last

        Parameters:
            points: list of three points
        """

        point_one = points[0]
        point_two = points[-1]

        # y = mx + b
        m = (point_two[1] - point_one[1]) / (point_two[0] - point_one[0])
        b = point_one[1] - m * point_one[0]

        return np.poly1d([m, b])

    @staticmethod
    def _parabolic_fit(points):
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
    def _curvature_area(points):
        """
        Takes in three points and gives
        the area of the shape made by a parabolic and linear fit

        Parameters:
            points: list of three points
        """

        x_min = points[0, 0]
        x_max = points[-1, 0]

        parabola = ManifoldMachine._parabolic_fit(points)
        line = ManifoldMachine._linear_fit(points)

        difference = lambda x: np.abs(parabola(x) - line(x))

        area, _ = spi.quad(difference, x_min, x_max)

        return area
