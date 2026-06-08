from collections import deque
from typing import Literal, Tuple, Optional
import logging

import numpy as np
import scipy.integrate as spi
from numpy.linalg import LinAlgError

from .Point import Point
from .FixedPoint import FixedPoint
from .DynamicalSystem import DynamicalSystem
from .BaseManifold import BaseManifold
from .ManifoldView import ManifoldView
from .BranchPoint import BranchPoint
from .Bridge import Bridge

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.setLevel(logging.INFO)


class ManifoldMachine:

    def __init__(self, system: DynamicalSystem):

        self.system = system
        self.area_cutoff = 1e-4

    def grow_manifold(
        self, fixed_point: FixedPoint, stability: Literal["unstable", "stable"]
    ):
        """
        Grows all the manifolds of a given stability from the fixed_point.
        This will grow all the manifolds attached to each iterate of the point.

        Parameters:
            fixed_point: point to grow the manifolds from
            stability: stability of the manifolds you want to grow

        Note:
            Currently you have to reset the tail after using this
        """

        branches = []
        iterated_branches = []
        temp_roots = []
        num_iterations = fixed_point.k_value

        # iterate all manifolds
        for i in range(num_iterations):

            branch_point = fixed_point.branch_points[i]

            # TODO fix this to account for inversion/reflection
            branch_num = 0  # if stability == "unstable" else 1

            manifold = BaseManifold(
                branch_point,
                stability,
                stretch_param=1,
                fixed_point=fixed_point,
                branch_index=branch_num,
            )
            branches.append(manifold)

            # this line meant to deal with roots that are fixed points
            # if it is an intersection behavior is unknown
            temp_root = manifold.root
            temp_roots.append(temp_root)
            if isinstance(manifold.root, BranchPoint):
                manifold.root = manifold.walk_fwd(None, temp_root)

            # you must walk forward before adding the stretch param
            # because the fixed point does not have one
            manifold.stretch_param = manifold.root.stretch_param

            manifold_iterate = self._iterate_without_refine(manifold)
            iterated_branches.append(manifold_iterate)

        # merge all manifolds
        iterated_branches = self._shift_list(iterated_branches)
        merged_manifolds = []
        for i in range(num_iterations):
            merged = self.merge_manifolds(branches[i], iterated_branches[i])
            merged_manifolds.append(merged)

        # refine all manifolds
        for i in range(num_iterations):
            self.refine_manifold(merged_manifolds[i])

        # reset all the roots
        for i in range(num_iterations):
            merged_manifolds[i].root = temp_roots[i]

        if num_iterations == 1:
            return merged_manifolds[0]

    # def new_grow_manifold(
    #     self, fixed_point: FixedPoint, stability: Literal["unstable", "stable"]
    # ):

    #     # construct a list of orbit indices based off the stability
    #     # this starts at the most recently iterated index and iterates
    #     # either clockwise or counterclockwise based on stability
    #     orbit_indices = fixed_point.get_iterable_array(stability, shift=1)

    #     current_manifold = BaseManifold(
    #         fixed_point.branch_points[orbit_indices[0]],
    #         stability,
    #         stretch_param=1,
    #         fixed_point=fixed_point,
    #         branch_index=0,
    #     )

    #     # this line meant to deal with roots that are fixed points
    #     # if it is an intersection behavior is unknown
    #     temp_root = current_manifold.root
    #     if isinstance(current_manifold.root, BranchPoint):
    #         current_manifold.root = current_manifold.walk_fwd(None, temp_root)

    #     # you must walk forward before adding the stretch param
    #     # because the fixed point does not have one
    #     current_manifold.stretch_param = current_manifold.root.stretch_param

    #     for i in range(fixed_point.period):

    #         for branch_index in fixed_point.get_branch_array():

    #             iterated_manifold = self.iterate_manifold(current_manifold)

    #             next_index = (i + 1) % fixed_point.period
    #             next_orbit_idx = orbit_indices[next_index]
    #             # self.merge_manifolds(uniterated_manifolds[i], iterated_manifold)

    #             next_manifold = BaseManifold(
    #                 root=fixed_point.branch_points[next_orbit_idx],
    #                 stability=stability,
    #                 stretch_param=current_manifold.stretch_param,
    #                 fixed_point=fixed_point,
    #                 branch_index=branch_index,
    #             )

    #             temp_root = next_manifold.root
    #             if isinstance(next_manifold.root, BranchPoint):
    #                 next_manifold.root = next_manifold.walk_fwd(None, temp_root)

    #             current_manifold = next_manifold
    #             # current_manifold = self.merge_manifolds(
    #             #     next_manifold, iterated_manifold
    #             # )
    #             # current_manifold._find_tail()

    def new_grow_manifold(
        self,
        fixed_point: FixedPoint,
        stability: Literal["unstable", "stable"],
        branch_index: Optional[int] = None,
    ):
        orbit_indices = fixed_point.get_iterable_array(stability, shift=1)

        if fixed_point.check_inversion():
            if branch_index is not None:
                logger.warning(
                    "branch_index is ignored for inversion points; "
                    "both branches are grown together."
                )
            branches_to_grow = fixed_point.get_branch_array()
        else:
            branches_to_grow = (
                list(range(fixed_point.num_branches))
                if branch_index is None
                else [branch_index]
            )

        for b in branches_to_grow:
            current_manifold = BaseManifold(
                fixed_point.branch_points[orbit_indices[0]],
                stability,
                stretch_param=1,
                fixed_point=fixed_point,
                branch_index=b,
            )

            temp_root = current_manifold.root
            if isinstance(current_manifold.root, BranchPoint):
                current_manifold.root = current_manifold.walk_fwd(None, temp_root)

            if current_manifold.root is None:
                continue  # branch not initialized

            current_manifold.stretch_param = current_manifold.root.stretch_param

            for i in range(fixed_point.period):

                for bi in (
                    fixed_point.get_branch_array()
                    if fixed_point.check_inversion()
                    else [b]
                ):

                    iterated_manifold = self.iterate_manifold(current_manifold)

                    next_index = (i + 1) % fixed_point.period
                    next_orbit_idx = orbit_indices[next_index]

                    next_manifold = BaseManifold(
                        root=fixed_point.branch_points[next_orbit_idx],
                        stability=stability,
                        stretch_param=current_manifold.stretch_param,
                        fixed_point=fixed_point,
                        branch_index=bi,
                    )

                    temp_root = next_manifold.root
                    if isinstance(next_manifold.root, BranchPoint):
                        next_manifold.root = next_manifold.walk_fwd(None, temp_root)

                    current_manifold = next_manifold

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

        # TODO include input for num_iterations
        # I don't think we need to do that anymore actually
        non_iterated_coords = manifold.get_non_iterated_point_array()
        non_iterated_cdists = manifold.get_non_iterated_cdist_array()
        non_iterated_points = manifold.get_non_iterated_point_array(return_nodes=True)

        old_points = manifold.get_iterated_point_array(return_nodes=True)

        if len(non_iterated_coords):

            iterated_points = np.vstack(
                [viewer.map_fwd(p) for p in non_iterated_coords]
            )
            distances = (
                (manifold.stretch_param * non_iterated_cdists).astype(float).ravel()
            )

            xvals = iterated_points[:, 0]
            yvals = iterated_points[:, 1]

            new_points = [
                Point(x, y, float(cdist), stretch_param=manifold.stretch_param)
                for x, y, cdist in zip(xvals, yvals, distances)
            ]

            # TODO include num_iterates
            # I don't think we need to do that anymore actually
            for i, point in enumerate(non_iterated_points):
                if manifold.stability == "unstable":
                    point.insert_next_iterate(new_points[i])
                else:
                    point.insert_prev_iterate(new_points[i])

            new_iterated_points = initalizer.construct_manifold_from_point_list(
                new_points,
                manifold.stability,
                manifold.stretch_param,
                manifold.fixed_point,
                manifold.branch_index,
            )

        # if there were points that needed to be mapped
        if not len(non_iterated_coords) == 0:

            # Only create old_iterated_points if there are some
            if len(old_points) > 0:
                old_iterated_points = BaseManifold(
                    old_points[0],
                    manifold.stability,
                    manifold.stretch_param,
                    manifold.fixed_point,
                    tail=old_points[-1],
                    branch_index=manifold.branch_index,
                )

                mapped_manifold = self.merge_manifolds(
                    old_iterated_points, new_iterated_points
                )
                self.refine_manifold(mapped_manifold)
            else:
                # No old iterated points - this is the first iteration of this bridge
                # Just refine the newly iterated points
                self.refine_manifold(new_iterated_points)
                mapped_manifold = new_iterated_points

            assert (
                sorted(mapped_manifold.get_cdist_array())
                == mapped_manifold.get_cdist_array()
            ).all()

            return mapped_manifold

        else:
            # All points already have iterates - collect the iterated points
            # This happens with bridges cut from the middle of a manifold
            if len(old_points) == 0:
                raise ValueError("No points in manifold to iterate")

            # The old_points already have their iterates computed
            # Just need to create a new manifold/bridge from those iterates
            # The iterated points are already geometrically connected in the iterate chain

            # Get the first and last iterated points
            first_point = old_points[0]
            last_point = old_points[-1]

            if manifold.stability == "unstable":
                first_iterate = first_point.next_iterate
                last_iterate = last_point.next_iterate
            else:
                first_iterate = first_point.prev_iterate
                last_iterate = last_point.prev_iterate

            if first_iterate is None or last_iterate is None:
                raise ValueError(
                    f"Points claim to have iterates but next_iterate/prev_iterate is None"
                )

            # Check if the input was a Bridge - if so, return a Bridge
            # I might want to do this for the other case as well, not just
            # if none of the points had already been iterated
            if isinstance(manifold, Bridge):
                iterated_bridge = Bridge(
                    root=first_iterate,
                    stability=manifold.stability,
                    stretch_param=manifold.stretch_param,
                    fixed_point=manifold.fixed_point,
                    tail=last_iterate,
                    name=manifold.name,
                    branch_index=manifold.branch_index,
                )
                return iterated_bridge
            else:
                # Construct a new BaseManifold
                iterated_manifold = BaseManifold(
                    root=first_iterate,
                    stability=manifold.stability,
                    stretch_param=manifold.stretch_param,
                    fixed_point=manifold.fixed_point,
                    tail=last_iterate,
                    branch_index=manifold.branch_index,
                )

                return iterated_manifold

    def iterate_bridge(self, manifold: Bridge):
        """
        Iterates a bridge forward and returns another bridge.

        Args:
            manifold (Bridge): Bridge to titerate forward.
        """

        # TODO consider how this method is handling bridge classes
        iterated_manifold = self.iterate_manifold(manifold)

        # we want to check if the resulting manifold conforms to our bridge standards
        # That could happen in Bridge if we want it to

        return iterated_manifold

    def cut_manifold(self, manifold: BaseManifold) -> list[Bridge]:
        """
        Takes a manifold and cuts it into Bridges that connect the intersection points.

        Args:
            manifold (BaseManifold): Manifold to be cut up.
        """

        final_node = manifold.tail

        previous_point = manifold.root
        current_point = manifold.walk_fwd(None, previous_point)

        bridges = []

        left_intersection = None
        right_intersection = None

        if isinstance(previous_point, BranchPoint):
            left_intersection = previous_point
            forming_bridge = True
        else:
            forming_bridge = False

        while current_point is not None:

            if isinstance(current_point, BranchPoint):
                if forming_bridge:
                    right_intersection = current_point
                else:
                    left_intersection = current_point
                    cached_previous_point = previous_point
                    forming_bridge = True

                if self._check_bridge_readiness(
                    left_intersection, right_intersection, forming_bridge
                ):

                    new_root = cached_previous_point
                    new_tail = manifold.walk_fwd(previous_point, right_intersection)

                    new_bridge = Bridge(
                        new_root,
                        manifold.stability,
                        manifold.stretch_param,
                        manifold.fixed_point,
                        manifold.name,
                        new_tail,
                        manifold.branch_index,
                    )

                    bridges.append(new_bridge)
                    left_intersection = right_intersection
                    cached_previous_point = previous_point
                    right_intersection = None
                    forming_bridge = True

            if current_point is final_node:
                break

            next_point = manifold.walk_fwd(previous_point, current_point)

            previous_point, current_point = current_point, next_point

        return bridges

    def _check_bridge_readiness(self, left, right, toggle) -> bool:
        """
        Helper function which returns True if a bridge is ready to be formed.

        Args:
            left (BranchPoint): First point which may form a bridge.
            right (BranchPoint): Second point which may form a bridge.
            toggle (bool): Flag telling if we are currently forming a bridge.

        Returns:
            bool: True if a bridge can be formed from left and right.
        """

        first_term = left is not None and right is not None

        return first_term and toggle

    def iterate_x_times(self, manifold: BaseManifold, num_times=1):
        """
        Iterates the manifold x times and returns a new manifold

        Parameters:
            manifold: manifold to iterate
            num_times: number of times to iterate
        """

        current_iterate = manifold

        for _ in range(num_times):
            current_iterate = self.iterate_manifold(current_iterate)

        return current_iterate

    def _iterate_without_refine(self, manifold: BaseManifold):
        """
        Iterates a manifold and returns the result without refinement
        Meant to be used internally for growing manifolds.

        Parameters:
            manifold: manifold to be iterated
        """

        from .ManifoldInitializer import ManifoldInitializer

        initalizer = ManifoldInitializer(self.system)
        viewer = ManifoldView(manifold, self.system)

        # TODO include input for num_iterations
        non_iterated_coords = manifold.get_non_iterated_point_array()
        non_iterated_cdists = manifold.get_non_iterated_cdist_array()
        non_iterated_points = manifold.get_non_iterated_point_array(return_nodes=True)

        number = [
            False if point.next_iterate is None else True
            for point in non_iterated_points
        ].count(True)
        print(f"Num Incorrectly labeled points: {number}")

        old_points = manifold.get_iterated_point_array(return_nodes=True)

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

            # TODO include num_iterates
            for i, point in enumerate(non_iterated_points):
                if manifold.stability == "unstable":
                    point.insert_next_iterate(new_points[i])
                else:
                    point.insert_prev_iterate(new_points[i])

            new_iterated_points = initalizer.construct_manifold_from_point_list(
                new_points,
                manifold.stability,
                manifold.stretch_param,
                manifold.fixed_point,
                manifold.branch_index,
            )

        # TODO include num_iterates
        # old_points = manifold.get_iterated_point_array(return_nodes=True)

        number = [True if point is None else False for point in old_points].count(True)
        print(f"Num Incorrectly labeled points: {number}")

        old_iterated_points = BaseManifold(
            old_points[0],
            manifold.stability,
            manifold.stretch_param,
            manifold.fixed_point,
            tail=old_points[-1],
            branch_index=manifold.branch_index,
        )

        # if there were points that needed to be mapped
        if not len(non_iterated_coords) == 0:

            mapped_manifold = self.merge_manifolds(
                old_iterated_points, new_iterated_points
            )

            assert (
                sorted(mapped_manifold.get_cdist_array())
                == mapped_manifold.get_cdist_array()
            ).all()

            return mapped_manifold

        else:
            return old_iterated_points

    def grow_x_times(
        self,
        fixed_point: FixedPoint,
        stability: Literal["unstable", "stable"],
        num_times=1,
        branch_index: Optional[int] = None,
    ):
        """
        Grows all the manifolds of the given stability from the fixed point
        by iterating them num_times and merging them back together

        Parameters:
            fixed_point: fixed point to grow the manifolds from
            stability: stability of the manifolds to grow
            num_times: number of times iterating the manifolds
        """

        for _ in range(num_times):
            self.new_grow_manifold(fixed_point, stability, branch_index)

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

        # NOTE If there are negative cdists this merge method will fail
        # Choose which manifold to start with based on the lower cdist
        if head_1.cdist <= head_2.cdist:
            current_point = head_1
            head_1 = manifold_1.walk_fwd(None, head_1)
            output_manifold = manifold_1
        else:
            current_point = head_2
            head_2 = manifold_2.walk_fwd(None, head_2)
            output_manifold = manifold_2

        # set the termination point right after the manifold tail
        over_one_1 = manifold_1.walk_fwd(None, manifold_1.tail)
        over_one_2 = manifold_2.walk_fwd(None, manifold_2.tail)

        while head_1 is not over_one_1 and head_2 is not over_one_2:

            if head_1.cdist < head_2.cdist:
                next_head = manifold_1.walk_fwd(None, head_1)
                if manifold_1.walk_fwd(None, current_point) is not head_1:
                    self._insert_point_geometrically(current_point, head_1, manifold_1)
                current_point = head_1
                head_1 = next_head
                # head_1 = manifold_1.walk_fwd(None, head_1)

            elif head_2.cdist < head_1.cdist:
                next_head = manifold_2.walk_fwd(None, head_2)
                if manifold_2.walk_fwd(None, current_point) is not head_2:
                    self._insert_point_geometrically(current_point, head_2, manifold_2)
                current_point = head_2
                head_2 = next_head
                # head_2 = manifold_2.walk_fwd(None, head_2)

            else:  # they are the same node
                if head_1 is head_2:
                    if manifold_1.walk_fwd(None, current_point) is not head_1:
                        self._insert_point_geometrically(
                            current_point, head_1, manifold_1
                        )
                    head_1 = manifold_1.walk_fwd(None, head_1)
                    head_2 = manifold_2.walk_fwd(None, head_2)
                    current_point = manifold_1.walk_fwd(None, current_point)
                else:
                    # Two different points share the same cdist.  One is
                    # newly generated (no iterate links); the other is the
                    # keeper.  Determine the keeper first, then advance the
                    # appropriate branch and skip the duplicate — no rewiring.
                    logger.debug(
                        "Duplicate cdist collision at cdist=%s; coords %s vs %s",
                        head_1.cdist,
                        head_1.get_point(),
                        head_2.get_point(),
                    )
                    head_1_has_iter = (
                        head_1.next_iterate is not None
                        or head_1.prev_iterate is not None
                    )
                    head_2_has_iter = (
                        head_2.next_iterate is not None
                        or head_2.prev_iterate is not None
                    )
                    if head_1_has_iter or not head_2_has_iter:
                        # head_1 is the keeper (or neither has iterates — pick head_1)
                        next_head_1 = manifold_1.walk_fwd(None, head_1)
                        if manifold_1.walk_fwd(None, current_point) is not head_1:
                            self._insert_point_geometrically(
                                current_point, head_1, manifold_1
                            )
                        current_point = head_1
                        head_1 = next_head_1
                        head_2 = manifold_2.walk_fwd(None, head_2)
                    else:
                        # head_2 is the keeper
                        next_head_2 = manifold_2.walk_fwd(None, head_2)
                        if manifold_2.walk_fwd(None, current_point) is not head_2:
                            self._insert_point_geometrically(
                                current_point, head_2, manifold_2
                            )
                        current_point = head_2
                        head_2 = next_head_2
                        head_1 = manifold_1.walk_fwd(None, head_1)

            # current_point = manifold_1.walk_fwd(None, current_point)

        # manifold_2 emptied first
        # WARNING: This may fail if we are merging manifold segments that
        # continue on past the tail node. This scheme doesn't insert the
        # rest of the remaining manifold, it only inserts the next node (head_1)
        # the fix is to not use this the _insert_point_geometrically method
        # probably just set
        # current_point.forward = head_1 and head_1.backward = current_point
        # but smartly based on stability of course
        # it is a bit trickier than that actually because we need the tail of
        # the non_exausted manifold to link to the over_one of the other manifold
        # NOTE This should be fixed now :)
        # manifold_2 emptied first
        if head_1 is not over_one_1:
            self._insert_point_geometrically(
                current_point, head_1, manifold_1, only_forward=True
            )
            output_manifold.tail = manifold_1.tail
            # self._insert_point_geometrically(
            #     output_manifold.tail, over_one_2, output_manifold, only_forward=True
            # )
        else:
            self._insert_point_geometrically(
                current_point, head_2, manifold_2, only_forward=True
            )
            output_manifold.tail = manifold_2.tail
            # self._insert_point_geometrically(
            #     output_manifold.tail, over_one_1, output_manifold, only_forward=True
            # )

        return output_manifold

    def refine_manifold(
        self, manifold: BaseManifold, branch_index=None, final_node=None
    ) -> set[Point]:
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

        num_initial_points = (
            len(manifold.get_point_array())
            if logger.isEnabledFor(logging.INFO)
            else None
        )

        previous_point = manifold.root
        current_point = manifold.walk_fwd(None, previous_point)

        modified_points = set()
        while current_point is not None:

            logger.debug(
                "Checking these points: %s",
                (previous_point.get_point(), current_point.get_point()),
            )

            added_points = self.refine_two_points(
                (previous_point, current_point), manifold, branch_index
            )
            modified_points.update(added_points)

            if current_point is final_node:
                break

            next_point = manifold.walk_fwd(previous_point, current_point)

            previous_point, current_point = current_point, next_point

        # logging logic
        if num_initial_points is not None:
            num_final_points = len(manifold.get_point_array())
            logger.info(
                "%d Points added during refinement",
                num_final_points - num_initial_points,
            )
            logger.info("%d NUMBER OF CURRENT POINTS", num_final_points)

        return modified_points

    def refine_two_points(
        self,
        points: tuple[Point | BranchPoint, Point | BranchPoint],
        manifold: BaseManifold,
        branch_index=None,
    ) -> set[Point]:
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

        modified_points = set()
        while pair_queue:

            p0, p1 = pair_queue.pop()

            # if one of the points is past the manifold continue
            if p0 is None or p1 is None:
                continue

            # Check if points are so close together to cause numerical instability
            if np.linalg.norm(p1.get_point() - p0.get_point()) < 1e-8:
                continue

            try:
                # curvature_area = self._curvature_area(triplet, viewer)
                curvature_area = self._curvature_area((p0, p1), viewer)
            except LinAlgError:
                logger.debug(
                    "Singular Vandermond matrix at cdist %.3gâ€“%.3g; skipping. "
                    "Points are likely near vertical",
                    p0.cdist,
                    p1.cdist,
                )
                continue

            if abs(curvature_area) < self.area_cutoff:
                continue

            # Skip refinement if either point's preiterate is missing.
            # This happens when bridge boundary points (root/tail inserted by
            # create_bridges) are mixed into a manifold with already-iterated
            # interior points and the iterate chain is incomplete for some members
            # of the merged manifold.
            if (ManifoldMachine._get_preiterate(p0, manifold.stability, 1) is None
                    or ManifoldMachine._get_preiterate(p1, manifold.stability, 1) is None):
                logger.debug(
                    "Skipping refinement for pair at cdist %.3g–%.3g: "
                    "preiterate missing on one or both points.",
                    p0.cdist,
                    p1.cdist,
                )
                continue

            new_point = self._get_refined_point(p0, p1, viewer, manifold.stability)

            self._insert_point_geometrically(p0, new_point, manifold, branch_index)
            # modified_points.update((p0, p1, p2, new_point))
            modified_points.update((p0, p1, new_point))

            pair_queue.append((p0, new_point))
            pair_queue.append((new_point, p1))

        return modified_points

    def _get_refined_point(
        self,
        p0: Point | BranchPoint,
        p1: Point | BranchPoint,
        viewer: ManifoldView,
        stability: Literal["unstable", "stable"],
    ):
        """ """

        # num_iterates = viewer.manifold.fixed_point.k_value
        num_iterates = 1

        p0_preiterate = self._get_preiterate(p0, stability, num_iterates)
        p1_preiterate = self._get_preiterate(p1, stability, num_iterates)

        new_point_coords_back = 0.5 * (
            p1_preiterate.get_point() + p0_preiterate.get_point()
        )

        new_distance = 0.5 * (float(p0.cdist) + float(p1.cdist))

        # new_point_coords = viewer.map_fwd(new_point_coords_back)
        new_point_coords = self._get_iterate(
            new_point_coords_back, viewer, num_iterates
        )

        x = new_point_coords[0]
        y = new_point_coords[1]
        new_point = Point(x, y, float(new_distance), stretch_param=p0.stretch_param)

        self._cache_preiterate(new_point, new_point_coords_back, stability)

        return new_point

    @staticmethod
    def _get_iterate(
        point: Point | BranchPoint, viewer: ManifoldView, num_iterates: int
    ):
        # for _ in range(num_iterates):
        point = viewer.map_fwd(point)

        return point

    @staticmethod
    def _shift_list(to_shift: list):
        """Shifts the list to the left one"""

        d = deque(to_shift)
        d.rotate(1)
        return list(d)

    @staticmethod
    def _get_preiterate(
        point: Point | BranchPoint,
        stability: Literal["unstable", "stable"],
        num_iterates: int = 1,
    ):
        """
        Helper function to get a point's preiterate based on stability

        Parameters:
            point (Point or BranchPoint): point to get preiterate
            stability (str): stability of the manifold
        """

        if stability == "unstable":
            return point.get_prev_iterate(num_iterates)
        else:
            return point.get_next_iterate(num_iterates)

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
        only_forward: bool = False,
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
                p0.insert_point_forward(
                    new_point, branch_index=branch_index, only_forward=only_forward
                )
            else:
                if p0.backward_branches[branch_index] is new_point:
                    return
                p0.insert_point_backward(
                    new_point, branch_index=branch_index, only_forward=only_forward
                )

        else:
            if manifold.stability == "unstable":
                if p0.forward is new_point:
                    return
                p0.insert_point_forward(new_point, only_forward=only_forward)
            else:
                if p0.backward is new_point:
                    return
                p0.insert_point_backward(new_point, only_forward=only_forward)

    @staticmethod
    def _linear_fit(points) -> Tuple[float, float]:
        """
        Takes in two points and gives the linear fit between the first and last

        Parameters:
            points: list of two points
        """

        if len(points) != 2:
            raise ValueError("Linear fit takes in two points")

        point_one = points[0]
        point_two = points[-1]

        # y = mx + b
        m = (point_two[1] - point_one[1]) / (point_two[0] - point_one[0])
        b = point_one[1] - m * point_one[0]

        return m, b

    @staticmethod
    def _parabolic_fit(points) -> Tuple[float, float, float]:
        """
        Takes in three points and gives the parabolic fit between them

        Parameters:
            points: list of three points
        """

        x_vals = points[:, 0]
        y_vals = points[:, 1]

        x0, x1, x2 = x_vals[0], x_vals[1], x_vals[2]

        # construct Vandermond matrix
        A = np.array([[x0**2, x0, 1], [x1**2, x1, 1], [x2**2, x2, 1]], dtype=float)

        Ainv = np.linalg.inv(A)

        # Solve for the coefficients [a, b, c].
        return tuple(Ainv @ y_vals)

    @staticmethod
    def _curvature_area(
        points: Tuple[Point | BranchPoint, Point | BranchPoint], viewer: ManifoldView
    ):
        """
        Takes in two points and computes the curvature between them.

        Fits 2 parabolas, a forward and a backward one
        Fits 1 line between the two points

        Computes 2 integrals from the difference tbetween the curves
        Returns the larger area

        Parameters:
            points: list of two points
        """

        if len(points) != 2:
            raise ValueError("Curvature Area takes in two points")

        p0 = points[0]
        p1 = points[1]

        point_vals = np.vstack((p0.get_point(), p1.get_point()))

        # compute the linear fit
        m, d = ManifoldMachine._linear_fit([p0.get_point(), p1.get_point()])

        # compute the left quadratic fit and area
        left = viewer.manifold.walk_back(p1, p0)
        if left is not None:

            left_points = np.vstack((left.get_point(), p0.get_point(), p1.get_point()))

            a_l, b_l, c_l = ManifoldMachine._parabolic_fit(left_points)

            left_area = abs(
                ManifoldMachine._compute_single_area(
                    point_vals, (m, d), (a_l, b_l, c_l)
                )
            )
        else:
            left_area = 0

        # compute the right quadratic fit and area
        right = viewer.manifold.walk_fwd(p0, p1)
        if right is not None:

            right_points = np.vstack(
                (p0.get_point(), p1.get_point(), right.get_point())
            )

            a_r, b_r, c_r = ManifoldMachine._parabolic_fit(right_points)

            right_area = abs(
                ManifoldMachine._compute_single_area(
                    point_vals, (m, d), (a_r, b_r, c_r)
                )
            )
        else:
            right_area = 0

        # return whichever area is larger
        area = left_area if left_area > right_area else right_area

        return area

    @staticmethod
    def _compute_single_area(
        points: list,
        linear_coef: tuple[float, float],
        quad_coef: tuple[float, float, float],
    ):
        """
        Computes the area between a line and parabola
        given the coefficients of the two curves and
        the endpoints

        Parameters:
            points: two endpoints
            linear_coef: m, d the linear coefficiencts
            quad_coef: a, b, c the quadratic coefficients
        """

        # unpack coefficients
        m, d = linear_coef
        a, b, c = quad_coef

        x_min = points[0, 0]
        x_max = points[-1, 0]

        # coefficients for the analytic solution to the integral
        a = a / 3
        b = 0.5 * (b - m)
        c = c - d

        x_3 = x_max**3 - x_min**3
        x_2 = x_max**2 - x_min**2
        x_1 = x_max - x_min

        area = a * x_3 + b * x_2 + c * x_1

        return abs(area)
