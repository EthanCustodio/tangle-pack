from collections import deque
from typing import Literal, Tuple, Optional
import logging

import numpy as np
import scipy.integrate as spi

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

"""
Dev Notes:

NOT A BUG -- refinement "explosion" at high development is just exponential
manifold growth into the map's escape region, and we never need to grow that far.

Symptom: at the k=2.1 inner period-3 fp ``construct_fixed_point([[0,1],[-1,0],
[-1,1]])`` with a tight area_cutoff (1e-7), a manifold arm grows gradually for a
while and then the point count and runtime blow up in a single iteration (in one
seen case to ~2e5 points; pushing one further hit millions and tens of GB).

Diagnosis (instrumented, 2026-06): the blow-up is NOT a tight-fold refinement
artifact (the earlier guess in this note). It is the arm ESCAPING to infinity. Two
signals settle it, traced per growth iteration on both the stable [1,1] and the
unstable [0,1] arm of that fp:
  * The minimum segment length never shrinks -- it sits at the ~1.4e-6 seed
    spacing the whole time. A tight fold would drive it toward zero; there is no
    tight fold.
  * The point count tracks coordinate divergence exactly. For the unstable [0,1]
    arm: iter 8 -> 478 pts, max|coord| ~ 3; iter 9 -> ~9.9k pts, max|coord| ~ 2e3;
    iter 10 -> ~1.5e6 pts, max|coord| ~ 2e26. The arm flies off to infinity, and
    refinement faithfully (and pointlessly) tries to resolve segments that span an
    ever-larger, diverging curve -- each inserted preiterate-midpoint maps even
    farther out, a divergence feedback, not fold resolution.

So the manifold is simply growing exponentially and escaping, which is the
expected behaviour of these maps for arms that leave the trapped region (cf.
``test_high_stretch_period3_growth``, which stops at iter 4 noting "the arm escapes
by the fifth iterate"). There is nothing to fix in ``refine_manifold`` /
``_refine_layer``: a fully developed tangle needs only a handful of iterations, and
arms are meant to be grown with a stop condition (``grow_until_turnaround`` /
``grow_until_arclength`` / ``grown_until_intersection``), NOT a large fixed
iteration count that walks an escaping arm out to infinity. The ``< 1e-8`` spatial
guard in ``_refine_layer`` stays as-is; it is fine for the regime we actually use.
See project memory ``solver-rootcause-and-seed-step``.
"""


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

            iterated_points = viewer.map_fwd_batch(non_iterated_coords)
            distances = np.ravel(manifold.stretch_param * non_iterated_cdists)

            xvals = iterated_points[:, 0]
            yvals = iterated_points[:, 1]

            new_points = [
                Point(x, y, float(cdist), stretch_param=manifold.stretch_param)
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
                manifold.fixed_point,
                manifold.branch_index,
            )
            new_iterated_points.manifold_key = manifold.manifold_key

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

            # cdist is the manifold's ordering key, so the geometric list must be
            # cdist-sorted. It is only *non-decreasing*, not strictly increasing:
            # where a high-stretch fold has collapsed adjacent cdists to within a
            # float ULP the refiner bridges the gap with equal-cdist points (a tie
            # is harmless -- such points are spliced geometrically, never re-sorted).
            assert (np.diff(mapped_manifold.get_cdist_array(), axis=0) >= 0).all()

            return mapped_manifold

        else:
            # All points already have iterates - collect the iterated points
            # This happens with bridges cut from the middle of a manifold
            if len(old_points) == 0:
                raise ValueError("No points in manifold to iterate")

            # The old_points already have their iterates computed
            # Just need to create a new manifold/bridge from those iterates
            # The iterated points are already geometrically connected in the iterate chain

            # ``old_points`` are already the iterates (forward images) of this
            # manifold's points -- ``get_iterated_point_array`` returned each point's
            # image. The iterated manifold is therefore bracketed by those images
            # directly; do NOT step one more iterate forward (that image does not
            # exist yet and was the source of the spurious "next_iterate is None"
            # failure when re-iterating an already-iterated child bridge).
            first_iterate = old_points[0]
            last_iterate = old_points[-1]

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
        iterated_manifold.manifold_key = (
            manifold.manifold_key
        )  # propagate for Tangle lookup

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
        logger.debug("Num incorrectly labeled points: %d", number)

        old_points = manifold.get_iterated_point_array(return_nodes=True)

        if len(non_iterated_coords):

            iterated_points = viewer.map_fwd_batch(non_iterated_coords)
            distances = np.ravel(manifold.stretch_param * non_iterated_cdists)

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
                manifold.fixed_point,
                manifold.branch_index,
            )

        number = [True if point is None else False for point in old_points].count(True)
        logger.debug("Num incorrectly labeled points: %d", number)

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

            # cdist is the manifold's ordering key, so the geometric list must be
            # cdist-sorted. It is only *non-decreasing*, not strictly increasing:
            # where a high-stretch fold has collapsed adjacent cdists to within a
            # float ULP the refiner bridges the gap with equal-cdist points (a tie
            # is harmless -- such points are spliced geometrically, never re-sorted).
            assert (np.diff(mapped_manifold.get_cdist_array(), axis=0) >= 0).all()

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

            elif head_2.cdist < head_1.cdist:
                next_head = manifold_2.walk_fwd(None, head_2)
                if manifold_2.walk_fwd(None, current_point) is not head_2:
                    self._insert_point_geometrically(current_point, head_2, manifold_2)
                current_point = head_2
                head_2 = next_head

            else:  # they share a cdist
                if head_1 is head_2:
                    # The same physical node appears in both lists.
                    if manifold_1.walk_fwd(None, current_point) is not head_1:
                        self._insert_point_geometrically(
                            current_point, head_1, manifold_1
                        )
                    head_1 = manifold_1.walk_fwd(None, head_1)
                    head_2 = manifold_2.walk_fwd(None, head_2)
                    current_point = manifold_1.walk_fwd(None, current_point)
                else:
                    # Two different points collided onto the same cdist. cdist is
                    # the only ordering key, so once it saturates (one ULP at large
                    # magnitude exceeds the true separation near a fold) it can no
                    # longer order the pair -- keeping both and ordering by the tied
                    # key scrambles them into a spike. The newly generated point has
                    # no iterate links yet; keep the established (iterate-linked)
                    # point and skip the duplicate.
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
        Adds additional points in areas of the manifold with high curvature.

        Refines breadth-first: every consecutive pair of points is a candidate.
        For each pair a left and right parabola plus a chord are fit and the
        bounded area is measured; pairs whose area exceeds ``self.area_cutoff``
        get a new point inserted at the (forward image of the) midpoint of their
        preiterates. A pair that is refined spawns its two child pairs, which are
        checked in the next layer. All of a layer's midpoints are mapped forward
        in a single batched call, so the dynamical map is evaluated once per
        layer rather than once per inserted point.

        Parameters:
            manifold (BaseManifold): current manifold
            branch_index: Optional if starting at a fixed point
            final_node: Unused; kept for backwards-compatible signature.

        Note:
            The per-pair decision depends only on local geometry and preiterates,
            so the breadth-first traversal reaches the same ``area < cutoff``
            stopping condition as the old depth-first scheme.
        """
        logger.debug("Refining manifold: %r", manifold.get_point_array())

        num_initial_points = (
            len(manifold.get_point_array())
            if logger.isEnabledFor(logging.INFO)
            else None
        )

        viewer = ManifoldView(manifold, self.system)

        # Layer 0: every consecutive geometric pair from root to tail.
        layer: list[tuple] = []
        previous_point = manifold.root
        current_point = manifold.walk_fwd(None, previous_point)
        while current_point is not None:
            layer.append((previous_point, current_point))
            if current_point is manifold.tail:
                break
            next_point = manifold.walk_fwd(previous_point, current_point)
            previous_point, current_point = current_point, next_point

        modified_points: set = set()
        while layer:
            layer = self._refine_layer(
                layer, manifold, viewer, branch_index, modified_points
            )

        if num_initial_points is not None:
            num_final_points = len(manifold.get_point_array())
            logger.debug(
                "%d points added during refinement (%d total)",
                num_final_points - num_initial_points,
                num_final_points,
            )

        return modified_points

    def _refine_layer(
        self,
        layer: list[tuple],
        manifold: BaseManifold,
        viewer: ManifoldView,
        branch_index,
        modified_points: set,
    ) -> list[tuple]:
        """
        Refine one breadth-first layer of pairs and return the next layer.

        Computes curvature areas for the whole layer at once, inserts a new
        point into each pair that needs refining (mapping all midpoints forward
        in a single batched call), and returns the child pairs of the refined
        pairs as the next layer.

        Args:
            layer: list of ``(p0, p1)`` node pairs to consider.
            manifold: the manifold being refined.
            viewer: view binding the manifold to the dynamical system.
            branch_index: branch index for inserting after a BranchPoint.
            modified_points: set accumulating every touched/created point.

        Returns:
            The next layer of ``(p0, p1)`` pairs (empty when refinement is done).
        """
        stability = manifold.stability

        # Gather node references and per-pair geometry, dropping pairs that are
        # degenerate or missing the data needed to refine (mirrors the old
        # per-pair guards: dead endpoints, near-coincident points, and missing
        # preiterates on bridge boundary points).
        p0_nodes, p1_nodes = [], []
        p0_xy, p1_xy, left_xy, right_xy = [], [], [], []
        pre_mid, c0_list, c1_list = [], [], []

        nan_row = np.array([np.nan, np.nan])

        for p0, p1 in layer:
            if p0 is None or p1 is None:
                continue

            a = p0.get_point()
            b = p1.get_point()
            if np.hypot(b[0] - a[0], b[1] - a[1]) < 1e-8:
                continue

            pre0 = ManifoldMachine._get_preiterate(p0, stability, 1)
            pre1 = ManifoldMachine._get_preiterate(p1, stability, 1)
            if pre0 is None or pre1 is None:
                logger.debug(
                    "Skipping refinement for pair at cdist %.3g-%.3g: "
                    "preiterate missing on one or both points.",
                    p0.cdist,
                    p1.cdist,
                )
                continue

            left = manifold.walk_back(p1, p0)
            right = manifold.walk_fwd(p0, p1)

            p0_nodes.append(p0)
            p1_nodes.append(p1)
            p0_xy.append(a)
            p1_xy.append(b)
            left_xy.append(left.get_point() if left is not None else nan_row)
            right_xy.append(right.get_point() if right is not None else nan_row)
            pre_mid.append(0.5 * (pre0.get_point() + pre1.get_point()))
            c0_list.append(float(p0.cdist))
            c1_list.append(float(p1.cdist))

        if not p0_nodes:
            return []

        areas = self._curvature_area_batch(
            np.array(p0_xy),
            np.array(p1_xy),
            np.array(left_xy),
            np.array(right_xy),
        )

        # Refine every pair whose curvature area exceeds the cutoff. The new point
        # takes the cdist midpoint of its neighbours; at a high-stretch fold the
        # curve has expanded so much that adjacent points are far apart in space
        # while their cdists have collapsed to within a float ULP, so that midpoint
        # may round exactly onto an endpoint. That tie is harmless here -- the point
        # is spliced geometrically between its neighbours, not re-sorted by cdist --
        # and bridging the gap is exactly what keeps the curve smooth instead of
        # leaving a coarse spike. The ``< 1e-8`` spacing guard above guarantees the
        # subdivision still terminates once the two points effectively coincide.
        c0_arr = np.array(c0_list)
        c1_arr = np.array(c1_list)
        mid_cd = 0.5 * (c0_arr + c1_arr)

        flagged = np.flatnonzero(areas >= self.area_cutoff)
        if flagged.size == 0:
            return []

        pre_mid = np.array(pre_mid)
        new_coords = viewer.map_fwd_batch(pre_mid[flagged])
        new_cdists = mid_cd[flagged]

        next_layer: list[tuple] = []
        for out_idx, pair_idx in enumerate(flagged):
            p0 = p0_nodes[pair_idx]
            p1 = p1_nodes[pair_idx]
            x, y = new_coords[out_idx]
            new_point = Point(
                x, y, float(new_cdists[out_idx]), stretch_param=manifold.stretch_param
            )
            # Cache the preiterate so the child pairs can be refined next layer.
            self._cache_preiterate(new_point, pre_mid[pair_idx], stability)

            self._insert_point_geometrically(p0, new_point, manifold, branch_index)
            modified_points.update((p0, p1, new_point))

            next_layer.append((p0, new_point))
            next_layer.append((new_point, p1))

        return next_layer

    def _get_refined_point(
        self,
        p0: Point | BranchPoint,
        p1: Point | BranchPoint,
        viewer: ManifoldView,
        stability: Literal["unstable", "stable"],
    ):
        """
        Reference single-pair refinement: the midpoint of the two preiterates
        mapped one step forward. The batched _refine_layer is the live path;
        this scalar version is kept as the pinned reference implementation.
        """

        p0_preiterate = self._get_preiterate(p0, stability)
        p1_preiterate = self._get_preiterate(p1, stability)

        new_point_coords_back = 0.5 * (
            p1_preiterate.get_point() + p0_preiterate.get_point()
        )

        new_distance = 0.5 * (float(p0.cdist) + float(p1.cdist))

        new_point_coords = self._get_iterate(new_point_coords_back, viewer)

        x = new_point_coords[0]
        y = new_point_coords[1]
        new_point = Point(x, y, float(new_distance), stretch_param=p0.stretch_param)

        self._cache_preiterate(new_point, new_point_coords_back, stability)

        return new_point

    @staticmethod
    def _get_iterate(point: Point | BranchPoint, viewer: ManifoldView):
        """Map a point one step forward in the manifold's stability direction."""
        return viewer.map_fwd(point)

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
        Takes in three points and gives the parabolic fit ``y = a x^2 + b x + c``
        between them.

        Uses the closed-form divided-difference (Newton) solution rather than
        building and inverting the 3x3 Vandermonde matrix. This is both faster
        (no ``np.linalg.inv``) and the basis for the vectorized curvature path.

        Parameters:
            points: list of three points

        Raises:
            ZeroDivisionError: if two of the three x-values coincide (the
                divided differences are then singular). Callers in the chord
                frame avoid this for non-degenerate triplets.
        """

        x_vals = points[:, 0]
        y_vals = points[:, 1]

        x0, x1, x2 = x_vals[0], x_vals[1], x_vals[2]
        y0, y1, y2 = y_vals[0], y_vals[1], y_vals[2]

        s01 = (y1 - y0) / (x1 - x0)
        s12 = (y2 - y1) / (x2 - x1)
        a = (s12 - s01) / (x2 - x0)
        b = s01 - a * (x0 + x1)
        c = y0 - a * x0**2 - b * x0

        return a, b, c

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

        # Work in a frame where the chord p0->p1 lies on the x-axis. The curvature
        # area is a geometric (rotation-invariant) quantity, but the line/parabola
        # fits are expressed as y-of-x, so a near-vertical segment would otherwise
        # give a singular Vandermonde matrix (and the chord slope would blow up).
        # Rotating the local points onto the chord keeps every fit well-conditioned
        # regardless of the segment's orientation.
        rotate = ManifoldMachine._chord_frame(p0.get_point(), p1.get_point())

        point_vals = np.vstack((rotate(p0.get_point()), rotate(p1.get_point())))

        # chord is the x-axis in this frame: m ~ 0, d ~ 0
        m, d = ManifoldMachine._linear_fit([point_vals[0], point_vals[1]])

        # compute the left quadratic fit and area
        left = viewer.manifold.walk_back(p1, p0)
        if left is not None:

            left_points = np.vstack(
                (rotate(left.get_point()), point_vals[0], point_vals[1])
            )

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
                (point_vals[0], point_vals[1], rotate(right.get_point()))
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
    def _curvature_area_batch(
        p0: np.ndarray,
        p1: np.ndarray,
        left: np.ndarray,
        right: np.ndarray,
    ) -> np.ndarray:
        """
        Vectorized form of :meth:`_curvature_area` over many pairs at once.

        Mirrors the scalar routine exactly: for each row the chord ``p0->p1`` is
        rotated onto the x-axis, a left parabola through ``(left, p0, p1)`` and a
        right parabola through ``(p0, p1, right)`` are fit with the closed-form
        divided differences, the closed-form integral between each parabola and
        the chord is taken, and the larger of the two areas is returned.

        Args:
            p0: ``(M, 2)`` first endpoints.
            p1: ``(M, 2)`` second endpoints.
            left: ``(M, 2)`` left neighbours; rows that are ``NaN`` (no
                neighbour) contribute ``0`` on the left side.
            right: ``(M, 2)`` right neighbours; ``NaN`` rows contribute ``0`` on
                the right side.

        Returns:
            ``(M,)`` array of curvature areas. Degenerate rows (zero-length
            chord, coincident x-values after rotation) yield ``0`` so they are
            never flagged for refinement.
        """
        p0 = np.asarray(p0, dtype=float)
        p1 = np.asarray(p1, dtype=float)
        left = np.asarray(left, dtype=float)
        right = np.asarray(right, dtype=float)

        chord = p1 - p0
        length = np.hypot(chord[:, 0], chord[:, 1])
        safe_len = np.where(length == 0.0, 1.0, length)
        ux = chord[:, 0] / safe_len
        uy = chord[:, 1] / safe_len

        def to_frame(pt):
            rel = pt - p0
            xf = ux * rel[:, 0] + uy * rel[:, 1]
            yf = -uy * rel[:, 0] + ux * rel[:, 1]
            return xf, yf

        x0, y0 = to_frame(p0)  # ~ (0, 0)
        x1, y1 = to_frame(p1)  # ~ (length, 0)

        # line through the two endpoints (the chord, ~ the x-axis here)
        with np.errstate(divide="ignore", invalid="ignore"):
            m = (y1 - y0) / (x1 - x0)
        d = y0 - m * x0

        def side_area(third):
            xt, yt = to_frame(third)
            with np.errstate(divide="ignore", invalid="ignore"):
                # parabola through (xt, yt), (x0, y0), (x1, y1)
                s_ab = (y0 - yt) / (x0 - xt)
                s_bc = (y1 - y0) / (x1 - x0)
                a = (s_bc - s_ab) / (x1 - xt)
                b = s_ab - a * (xt + x0)
                c = yt - a * xt**2 - b * xt

                # closed-form integral between parabola and chord (see
                # _compute_single_area): coefficients of (quad - line)
                aa = a / 3.0
                bb = 0.5 * (b - m)
                cc = c - d
                area = np.abs(
                    aa * (x1**3 - x0**3)
                    + bb * (x1**2 - x0**2)
                    + cc * (x1 - x0)
                )
            return np.where(np.isfinite(area), area, 0.0)

        left_area = side_area(left)
        right_area = side_area(right)

        area = np.maximum(left_area, right_area)
        area = np.where(length == 0.0, 0.0, area)
        return area

    @staticmethod
    def _chord_frame(p0_xy: np.ndarray, p1_xy: np.ndarray):
        """
        Return a function rotating points into the frame whose x-axis is the
        chord ``p0 -> p1`` (origin at ``p0``).

        Used so the curvature fits stay well-conditioned for segments of any
        orientation, including near-vertical ones. Falls back to the identity
        (translation only) for a degenerate zero-length chord.

        Args:
            p0_xy: Coordinates of the first endpoint.
            p1_xy: Coordinates of the second endpoint.

        Returns:
            Callable mapping an (x, y) point into the chord frame.
        """
        p0_xy = np.asarray(p0_xy, dtype=float).ravel()
        p1_xy = np.asarray(p1_xy, dtype=float).ravel()

        chord = p1_xy - p0_xy
        length = float(np.hypot(chord[0], chord[1]))
        if length == 0.0:
            return lambda pt: np.asarray(pt, dtype=float).ravel() - p0_xy

        ux, uy = chord[0] / length, chord[1] / length
        # rotation by -theta maps the chord direction onto +x
        rot = np.array([[ux, uy], [-uy, ux]], dtype=float)

        def to_frame(pt: np.ndarray) -> np.ndarray:
            return rot @ (np.asarray(pt, dtype=float).ravel() - p0_xy)

        return to_frame

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
