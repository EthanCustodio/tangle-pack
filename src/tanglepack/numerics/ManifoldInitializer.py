"""
Seeding the initial fundamental segment of a manifold from eigendata.

:class:`ManifoldInitializer` turns a solved :class:`~.FixedPoint.FixedPoint`
into the short starting curve that :class:`~.ManifoldMachine.ManifoldMachine`
then grows: it steps off the fixed point along the eigenvector, maps that seed
around the orbit, and orients the eigenvectors consistently.

Dev Notes:

One construction: one fundamental segment is seeded per chain and mapped
around it (``construct_kevin_way``).

branch_index: a branch index names an
EIGENDIRECTION, and each one seeds one chain of ``k_value`` pieces walked by
``FixedPoint.advance_key``. An inversion chain covers both branches, so asking
for branch 1 there adds nothing; without inversion branch 1 is its own chain.

Open requests (features, for the API plan): orient_manifolds only accepts
approximate directions for the 0th iterate of the cycle, and it requires both a
stable and an unstable direction. Allowing an arbitrary iterate, and allowing
one of the two directions to be given on its own, are both open.

Seed step. get_first_point seeds the fundamental segment a distance
``max(fixed_point.accuracy, _MIN_SEED_STEP)`` along the eigenvector. The dynamic
``accuracy`` term (orbit residual^(1/3)) scales the step with how well the orbit is
located and keeps it comfortably inside the linear regime around the fixed point --
this is the intended behaviour. The floor only guards the degenerate case: once the
solver (fsolve) locates an orbit exactly (the rational period-3 orbit at Henon k=2),
the residual is 0 and the bare step would be 0, giving a 0/0 nan. Note the side
effect of an accurate solver: for a generic k the residual is ~machine eps so the
seed is ~5e-6 and tangles develop slowly per iteration (many
iterations are needed for a full tangle). The point-count blow-up of a late step
(at k=2.1 fp3 stable: ~628 points at iter 6, ~240k at iter 7)
is NOT an open refinement bug -- see the ManifoldMachine Dev Notes: the arm is
escaping to infinity, and the answer is to grow with a stop condition rather than a
large fixed iteration count.
"""

from __future__ import annotations

from typing import Dict, Tuple, Optional
from numpy.typing import NDArray
from .FixedPoint import FixedPoint
from .DynamicalSystem import DynamicalSystem
from .ManifoldMachine import ManifoldMachine
from .BaseManifold import BaseManifold
from .Intersection import Stability
from .Point import Point
from .PointStore import PointStore
import numpy as np

import logging

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Floor for the seed-step length used to start a fundamental segment. The seed step
# itself is dynamic -- the located orbit's residual^(1/3) (``fixed_point.accuracy``)
# -- which scales the step with how well the orbit is pinned down and so keeps the
# seed comfortably inside the linear regime around the fixed point. The only failure
# mode is an exactly-located orbit (the rational period-3 orbit at Henon k=2), where
# the residual is 0 and the step would collapse to 0 (a 0/0 nan downstream). This
# floor is the "minimally small but numerically robust" length used in that case;
# it sits at the same scale (~cbrt machine eps) that a machine-precision orbit
# already produces, so k=2 stays continuous with its neighbours. See the seed-step
# dev note below.
_MIN_SEED_STEP = 5e-6


class ManifoldInitializer:
    """
    Toolbox getting the initial manifold segments from a fixed point.

    Attributes:
        system (DynamicalSystem): All the functions related to the dynamical system.
        machine (ManifoldMachine): The machine that maps and refines the initial
            fundamental segments of a ``k_value > 1`` orbit.
        store (PointStore): The store every seeded point is appended to.
    """

    def __init__(
        self,
        system: DynamicalSystem,
        machine: Optional[ManifoldMachine] = None,
        store: Optional[PointStore] = None,
    ) -> None:
        """
        Initializes the toolbox with the dynamical system and machine machinery

        Args:
            system (DynamicalSystem): All the functions related to the dynamical system.
            machine (Optional[ManifoldMachine]): The machine to refine with. The
                workbench passes its own, so the initial segments refine at ITS
                ``area_cutoff`` and under its point budget; a fresh machine at the
                library defaults when omitted.
            store (Optional[PointStore]): The store the seeded points go into
                (the workbench passes its own); a fresh one when omitted.
        """

        self.system = system
        self.machine = machine if machine is not None else ManifoldMachine(system)
        self.store = store if store is not None else PointStore()

    def get_first_point(
        self,
        fixed_point: FixedPoint,
        orbit_index: int,
        branch_index: int,
        stability: Stability,
    ) -> NDArray[np.float64]:
        """
        Computes the first point from the fixed point based on a linear interpolation.

        Note:
            This point will end up getting iterated towards the fixed point and be the
            second in the first fundamental segment.

        Args:
            fixed_point (FixedPoint): The fixed point to start from.
            orbit_index (int): The index corresponding to the iterate of the fixed
                point to start from.
            branch_index (int): The index corresponding to which branch to grow out.
                Currently not really being used.
            stability (Stability): Stability of the manifold
                to get the first point of.
        """

        # Dynamic seed step (orbit residual^(1/3)) keeps us in the linear regime;
        # floor it so an exactly-located orbit (residual 0, e.g. k=2) does not
        # collapse the step to 0 and produce a 0/0 nan downstream.
        step = max(fixed_point.accuracy, _MIN_SEED_STEP)

        if stability == "unstable":
            direction_from_fixed_point = fixed_point.unstable_eigenvectors[orbit_index]
        else:
            direction_from_fixed_point = fixed_point.stable_eigenvectors[orbit_index]

        direction_from_fixed_point = direction_from_fixed_point / np.linalg.norm(
            direction_from_fixed_point
        )

        if branch_index == 1 and not fixed_point.check_inversion():
            direction_from_fixed_point = -direction_from_fixed_point

        return fixed_point.coordinates[orbit_index] + step * direction_from_fixed_point

    def get_first_point_back(
        self,
        fixed_point: FixedPoint,
        orbit_index: int,
        branch_index: int,
        stability: Stability,
        first_point: Point,
    ) -> NDArray[np.float64]:
        """
        Get the iterate/preiterate (based on stability) of the first point.

        Args:
            fixed_point (FixedPoint): The fixed point we are attached to.
            orbit_index (int): The index of the iterate of the fixed point we are at.
            branch_index (int): The index of the branch of the fixed point we are at.
            stability (Stability): The stability of the manifold.
            first_point (Point): The point generated by taking a small step from
                the fixed_point.

        Return:
            np.ndarray: The (pre)iterate of the first point.
        """
        num_iterations = fixed_point.k_value

        first_back = first_point.get_point()
        curr_point = first_point

        if stability == "unstable":
            # this scheme establishes a 'fictitious' list of
            # iterates so we can insert the point we are looking for
            # at the proper iterate in the cycle
            # these fictitous points are named so because they are
            # not part of any geometrical manifold object
            # when the initial manifold segments are mande these points
            # are given cdists and made full dudes. Just not inserted
            # geometrically yet
            for _ in range(num_iterations):
                temp_back = self.system.map_inv(first_back)
                temp_point = Point(temp_back[0], temp_back[1], store=self.store)
                curr_point.insert_prev_iterate(temp_point)
                curr_point = curr_point.prev_iterate
                first_back = temp_back
            curr_point.next_iterate.prev_iterate = None

        else:  # stable branch
            for _ in range(num_iterations):
                temp_back = self.system.map(first_back)
                temp_point = Point(temp_back[0], temp_back[1], store=self.store)
                curr_point.insert_next_iterate(temp_point)
                curr_point = curr_point.next_iterate
                first_back = temp_back
            curr_point.prev_iterate.next_iterate = None

        return first_back

    def get_initial_fundamental_segment(
        self,
        fixed_point: FixedPoint,
        orbit_index: int,
        branch_index: int,
        stability: Stability,
    ) -> BaseManifold:
        """
        Computes the initial fundamental segment.
        Takes a small step from the fixed point then takes the
        iterate/preiterate (depends on stability) to get another point closer
        to the fixed point. Those two points make up the initial fundamental segment.

        Args:
            fixed_point (FixedPoint): The fixed point the segment is attached to.
            orbit_index (int): The iterate of the fixed point to get the segment from.
            branch_index (int): The branch of the fixed point the segment is from.
            stability (Stability): The stability of the manifold.

        Returns:
            BaseManifold: The resulting fundamental segment.
        """

        first_point = self.get_first_point(
            fixed_point, orbit_index, branch_index, stability
        )
        distance_first = np.linalg.norm(
            first_point - fixed_point.coordinates[orbit_index]
        )

        first_point = Point(
            first_point[0],
            first_point[1],
            cdist=distance_first,
            store=self.store,
        )

        first_back = self.get_first_point_back(
            fixed_point, orbit_index, branch_index, stability, first_point
        )
        distance_prev = np.linalg.norm(
            first_back - fixed_point.coordinates[orbit_index]
        )

        # distance_first / distance_prev is the growth over one full return to
        # this branch (k_value map steps); the fixed point converts it to the
        # per-map-step factor the points carry as their stretch parameter.
        alpha = fixed_point.per_step_factor(distance_first / distance_prev)

        first_point.stretch_param = alpha

        first_back = Point(
            first_back[0],
            first_back[1],
            cdist=distance_prev,
            stretch_param=alpha,
            store=self.store,
        )

        # The seed's (pre)iterate one full return back is first_back: the
        # unstable branch walks back along prev_iterate (the map's inverse), the
        # stable one along next_iterate, and the geometric list grows forward on
        # an unstable branch, backward on a stable one.
        unstable = stability == "unstable"
        if unstable:
            first_point.insert_prev_iterate(first_back, fixed_point.k_value)
        else:
            first_point.insert_next_iterate(first_back, fixed_point.k_value)

        # fill the cdist and stretch param of the 'fictitious' points added
        # when finding the (pre)iterate
        current_point = first_point
        for _ in range(fixed_point.k_value - 1):
            previous_point = current_point
            current_point = current_point.prev_iterate if unstable else current_point.next_iterate
            current_point.cdist = previous_point.cdist / alpha
            current_point.stretch_param = alpha

        branch_point = fixed_point.branch_points[orbit_index]
        insert = branch_point.insert_point_forward if unstable else branch_point.insert_point_backward
        insert(first_point, branch_index)
        insert(first_back, branch_index)

        return BaseManifold(
            fixed_point.branch_points[orbit_index],
            stability,
            alpha,
            fixed_point,
            tail=first_point,
            branch_index=branch_index,
            manifold_key=(fixed_point, stability, orbit_index, branch_index),
        )

    @staticmethod
    def orient_manifolds(
        fixed_point: FixedPoint,
        approx_dirs: Optional[dict[str, np.ndarray]] = None,
    ) -> None:
        """
        Allows the user to provide an approximate direction for the manifolds.
        This routine goes through the list of eigenvectors and checks if the
        approximate direction matches. Give the approximate direction for the
        0th iterate of the cycle and it will be mapped around. On an inversion
        point this chooses which half of each manifold is branch 0 (branch 1
        runs the other way), and it is propagated along the orbit the same way.

        Args:
            fixed_point (FixedPoint): Fixed point to orient the manifolds from.
            approx_dirs (dict): User supplied approximate directions of the
                manifolds to grow. It is a dictionary that looks like:
                    {"unstable": [ux, uy], "stable": [sx, sy]}.
        """

        unstable_approx = approx_dirs["unstable"]
        stable_approx = approx_dirs["stable"]

        unstable_dir = fixed_point.unstable_eigenvectors[0]
        stable_dir = fixed_point.stable_eigenvectors[0]

        # check if the eigenvectors and approximate directions align
        unstable_test = True if np.dot(unstable_approx, unstable_dir) > 0 else False
        stable_test = True if np.dot(stable_approx, stable_dir) > 0 else False

        # if the vectors do not align flip them
        if not unstable_test:
            fixed_point.unstable_eigenvectors[0] *= -1

        if not stable_test:
            fixed_point.stable_eigenvectors[0] *= -1

        period = fixed_point.period

        for i in range(period - 1):

            partial_jacobian = fixed_point.partial_jacobians[i]

            unstable_approx = partial_jacobian @ fixed_point.unstable_eigenvectors[i]
            stable_approx = partial_jacobian @ fixed_point.stable_eigenvectors[i]

            unstable_dir = fixed_point.unstable_eigenvectors[i + 1]
            stable_dir = fixed_point.stable_eigenvectors[i + 1]

            # check if the eigenvectors and approximate directions align
            unstable_test = True if np.dot(unstable_approx, unstable_dir) > 0 else False
            stable_test = True if np.dot(stable_approx, stable_dir) > 0 else False

            # if the vectors do not align flip them
            if not unstable_test:
                fixed_point.unstable_eigenvectors[i + 1] *= -1

            if not stable_test:
                fixed_point.stable_eigenvectors[i + 1] *= -1

    def construct_kevin_way(
        self,
        fixed_point: FixedPoint,
        stability: Stability,
        num_branches: int = 1,
    ) -> Dict[Tuple[int, int], BaseManifold]:
        """
        Constructs the initial segments by creating a single fundamental
        segment and mapping it all the way around the fixed point.

        Args:
            fixed_point (FixedPoint): Fixed point to grow the manifolds from.
            stability (Stability): Stability of the manifold.
            num_branches (int): Number of EIGENDIRECTIONS to seed. Must be 1 or 2.
                Each direction seeds one chain; a direction an earlier chain
                already covers (branch 1 of an inversion point) is skipped. For
                a point without inversion the two directions are two independent
                chains, so 2 is a legitimate request even though
                ``fixed_point.num_branches`` is 1 (that property counts the
                branches of ONE chain).

        Returns:
            Dict[Tuple[int, int], BaseManifold]: All initial fundamental segments.
                The tuple structure is (orbit_index, branch_index).

        Raises:
            ValueError: If ``num_branches`` is not 1 or 2.

        Note:
            ONE ordering rule: each chain is ``k_value`` pieces walked with
            :meth:`FixedPoint.advance_key`, ``+1`` for the unstable manifold and
            ``-1`` for the stable one (grown by the INVERSE map, so its orbit
            indices run ``0, p-1, ..., 1``). An inversion chain runs
            ``(0, 0), (1, 0), ..., (p-1, 0), (0, 1), ..., (p-1, 1)`` (reversed for
            the stable manifold) and so covers both branches.
        """

        if num_branches not in (1, 2):
            raise ValueError(f"num_branches must be 1 or 2, got {num_branches}")

        step = 1 if stability == "unstable" else -1
        initial_segments: Dict[Tuple[int, int], BaseManifold] = {}
        for branch_index in range(num_branches):
            key = (fixed_point, stability, 0, branch_index)
            if key[2:] in initial_segments:
                continue  # an inversion chain already covers both branches
            segment = self.get_initial_fundamental_segment(fixed_point, 0, branch_index, stability)
            segment.root = segment.first_node()
            initial_segments[key[2:]] = segment
            for _ in range(fixed_point.k_value - 1):
                key = fixed_point.advance_key(key, step)
                segment = self.machine.iterate_manifold(segment)
                # iterate_manifold builds the image without a branch identity (it
                # cannot advance the key itself); the chain walk knows it exactly.
                segment.manifold_key = key
                segment.branch_index = key[3]
                initial_segments[key[2:]] = segment
                self.machine._insert_point_geometrically(
                    fixed_point.branch_points[key[2]], segment.root, segment, key[3]
                )

        for (orbit_index, _branch), segment in initial_segments.items():
            segment.root = fixed_point.branch_points[orbit_index]

        return initial_segments
