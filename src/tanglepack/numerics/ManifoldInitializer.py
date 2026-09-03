from typing import Literal, Dict, Tuple, Optional
from numpy.typing import NDArray
from .FixedPoint import FixedPoint
from .BranchPoint import BranchPoint
from .DynamicalSystem import DynamicalSystem
from .ManifoldMachine import ManifoldMachine
from .BaseManifold import BaseManifold
from .Intersection import ManifoldKey
from .Point import Point
import numpy as np

import logging

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.setLevel(logging.INFO)

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

"""
Dev Notes:

Fully commit to one way of doing the manifold creation. 
Either Ethan way or Kevin way. The suspicion is that what was
wront with Ethan way was the mapping of the stable manifolds backwards 
in the right way like I did for the Kevin way. So if I fix that then 
I can just fully commit to that way. To do Ethan's way we will take 
the things I learned from Kevin's way like the output dict structure and such.
    There are two advantages and one disadvantage to doing Ethan way
    +1 It will be symmetrical, all manifolds will grow at the same rate
    +1 This class will have no dependency on ManifoldMachine which is awkward
    -1 p**2 unnessesary points are added to the manifold.

Fully sus out what to do about branch_index. It is unused in most of 
these methods, but I think it is given as an input elsewhere in the codebase.
Decide if that is necessary and make its inclusion in these methods more clear.
Potentially it is necessary to do things the Ethan way, in the Kevin way
I don't think it matters.

Look at the TODO in orient_manifold

Look at insert_point_geometrically. We have the same function in manifold machine.
We only need to have that functionality in one place. If the function here 
does something different then it needs a different name. If they both need this helper
then it needs to go down to a lower level. Potentially in BaseManifold. 

The get_initial_fundamental segment code is quite long.
Make it more concise or split it into multiple functions.

Seed step. get_first_point seeds the fundamental segment a distance
``max(fixed_point.accuracy, _MIN_SEED_STEP)`` along the eigenvector. The dynamic
``accuracy`` term (orbit residual^(1/3)) scales the step with how well the orbit is
located and keeps it comfortably inside the linear regime around the fixed point --
this is the intended behaviour. The floor only guards the degenerate case: once the
solver (fsolve) locates an orbit exactly (the rational period-3 orbit at Henon k=2),
the residual is 0 and the bare step would be 0, giving a 0/0 nan. Note the side
effect of an accurate solver: for a generic k the residual is ~machine eps so the
seed is ~5e-6 and tangles develop more slowly per iteration than they used to (more
iterations are needed for a full tangle). That interacts with the OPEN refinement
explosion -- see ManifoldMachine and the project memory: at k=2.1 fp3 stable the
point count is fine through iter 6 (~628) then blows up ~380x at iter 7 (~240k).
"""


class ManifoldInitializer:
    """
    Toolbox capable of getting the initial manifold segments from a fixed point
    and initializing manifolds from lists of points.

    Attributes:
        system (DynamicalSystem): All the functions related to the dynamical system.
        machine (ManifoldMachine): Toolbox containing methods for manifold
    """

    def __init__(self, system: DynamicalSystem):
        """
        Initializes the toolbox with the dynamical system and machine machinery

        Args:
            system (DynamicalSystem): All the functions related to the dynamical system.
        """

        self.system = system
        self.machine = ManifoldMachine(system)

    def get_first_point(
        self,
        fixed_point: FixedPoint,
        orbit_index: int,
        branch_index: int,
        stability: Literal["stable", "unstable"],
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
            stability (Literal["unstable", "stable]): Stability of the manifold
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

        direction_from_fixed_point = np.asarray(direction_from_fixed_point).flatten()
        direction_from_fixed_point /= np.linalg.norm(direction_from_fixed_point)

        if branch_index == 1 and not fixed_point.check_inversion():
            direction_from_fixed_point = -direction_from_fixed_point

        first_point = (
            fixed_point.coordinates[orbit_index] + (step) * direction_from_fixed_point
        )

        return np.array(first_point, dtype=np.float64).reshape(-1)

    def get_first_point_back(
        self,
        fixed_point: FixedPoint,
        orbit_index: int,
        branch_index: int,
        stability: Literal["stable", "unstable"],
        first_point: Point,
    ) -> NDArray[np.float64]:
        """
        Get the iterate/preiterate (based on stability) of the first point.

        Args:
            fixed_point (FixedPoint): The fixed point we are attached to.
            orbit_index (int): The index of the iterate of the fixed point we are at.
            branch_index (int): The index of the branch of the fixed point we are at.
            stability (Literal["unstable", "stable"]): The stability of the manifold.
            first_point (Point): The point generated by taking a small step from
                the fixed_point.

        Return:
            np.ndarray: The (pre)iterate of the first point.
        """
        num_iterations = fixed_point.k_value

        # first_back = self.get_first_point(
        #     fixed_point, orbit_index, branch_index, stability
        # )
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
                temp_point = Point(temp_back[0], temp_back[1])
                curr_point.insert_prev_iterate(temp_point)
                curr_point = curr_point.prev_iterate
                first_back = temp_back
            curr_point.next_iterate.prev_iterate = None

        else:  # stable branch
            for _ in range(num_iterations):
                temp_back = self.system.map(first_back)
                temp_point = Point(temp_back[0], temp_back[1])
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
        stability: Literal["stable", "unstable"],
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
            stability (Literal["unstable", "stable"]): The stability of the manifold.

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
            edist=distance_first,
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
            edist=distance_prev,
            stretch_param=alpha,
        )

        if stability == "unstable":
            # add in the new point to the iterate list
            first_point.insert_prev_iterate(first_back, fixed_point.k_value)

            # fill the cdist and stretch param info of the
            # 'fictitious' points added when finding the preiterate
            previous_point = None
            current_point = first_point
            for _ in range(fixed_point.k_value - 1):

                previous_point = current_point
                current_point = current_point.prev_iterate

                current_point.cdist = previous_point.cdist / alpha
                current_point.stretch_param = alpha

            # insert the new point in geometrically
            fixed_point.branch_points[orbit_index].insert_point_forward(
                first_point, branch_index
            )
            fixed_point.branch_points[orbit_index].insert_point_forward(
                first_back, branch_index
            )

        else:
            # add in the new point to the iterate list
            first_point.insert_next_iterate(first_back, fixed_point.k_value)

            # fill the cdist and stretch param info of the
            # 'fictitious' points added when finding the preiterate
            previous_point = None
            current_point = first_point
            for _ in range(fixed_point.k_value - 1):

                previous_point = current_point
                current_point = current_point.next_iterate

                current_point.cdist = previous_point.cdist / alpha
                current_point.stretch_param = alpha

            # insert the new point in geometrically
            fixed_point.branch_points[orbit_index].insert_point_backward(
                first_point, branch_index
            )
            fixed_point.branch_points[orbit_index].insert_point_backward(
                first_back, branch_index
            )

        return BaseManifold(
            fixed_point.branch_points[orbit_index],
            stability,
            alpha,
            fixed_point,
            tail=first_point,
            branch_index=branch_index,
            manifold_key=(fixed_point, stability, orbit_index, branch_index),
        )

    def construct_manifold_from_point_list(
        self,
        points: list[Point],
        stability: Literal["stable", "unstable"],
        stretch_param: float,
        fixed_point: FixedPoint,
        branch_index=None,
        *,
        manifold_key: Optional[ManifoldKey],
    ) -> BaseManifold:
        """
        Constructs a manifold from a list of Point objects.

        Note:
            The given list is assumed to be in cdist ordering.

        Args:
            points (list[Point]): List of points to turn into a manifold.
            stability (Literal["unstable", "stable"]): Stability of the resulting
                manifold.
            stretch_param (float): The stretch parameter for the new manifold.
            fixed_point (FixedPoint): The fixed point the new manifold emanates from.
            branch_index (int, optional): The branch index of the fixed point this
                manifold emanates from.
            manifold_key (Optional[ManifoldKey]): Keyword-only and required -- the
                branch the new manifold is (``None`` only for a transient forward
                image whose branch the caller has yet to advance).

        Returns:
            BaseManifold: The resulting manifold made from the list.
        """

        manifold = BaseManifold(
            points[0],
            stability,
            stretch_param,
            fixed_point,
            manifold_key=manifold_key,
        )

        current_point = manifold.root

        for i, point in enumerate(points):

            # start inserting the second point in the list
            if i == 0:
                continue

            self._insert_point_geometrically(
                current_point, point, manifold, branch_index
            )

            current_point = point

        manifold.tail = points[-1]
        return manifold

    def _insert_point_geometrically(
        self, p0: Point, new_point: Point, manifold: BaseManifold, branch_index=None
    ):
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

    @staticmethod
    def orient_manifolds(
        fixed_point: FixedPoint,
        approx_dirs: dict[str, np.ndarray] | None = None,
    ):
        """
        Allows the user to provide an approximate direction for the manifolds.
        This routine goes through the list of eigenvectors and checks if the
        approximate direction matches. Give the approximate direction for the
        0th iterate of the cycle and it will be mapped around.

        TODO make it so you can provide the approx_dirs for an arbitrary iterate
        TODO make it so you can pass in either stable or unstable dir

        Args:
            fixed_point (FixedPoint): Fixed point to orient the manifolds from.
            approx_dirs (dict): User supplied approximate directions of the
                manifolds to grow. It is a dictionary that looks like:
                    {"unstable": [ux, uy], "stable": [sx, sy]}.
        """

        # if there is inversion there is no need to flip
        if fixed_point.check_inversion():
            return None

        unstable_approx = approx_dirs["unstable"]
        stable_approx = approx_dirs["stable"]

        unstable_dir = fixed_point.unstable_eigenvectors[0].ravel()
        stable_dir = fixed_point.stable_eigenvectors[0].ravel()

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

            unstable_approx = (
                partial_jacobian @ fixed_point.unstable_eigenvectors[i]
            ).ravel()
            stable_approx = (
                partial_jacobian @ fixed_point.stable_eigenvectors[i]
            ).ravel()

            unstable_dir = fixed_point.unstable_eigenvectors[i + 1].ravel()
            stable_dir = fixed_point.stable_eigenvectors[i + 1].ravel()

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
        stability: Literal["unstable", "stable"],
        num_branches: int = 1,
    ) -> Dict[Tuple[int, int], BaseManifold]:
        """
        Constructs the initial segments by creating a single fundamental
        segment and mapping it all the way around the fixed point.

        Args:
            fixed_point (FixedPoint): Fixed point to grow the manifolds from.
            stability (Literal["unstable", "stable"]): Stability of the manifold.
            num_branches (int): Number of EIGENDIRECTIONS to seed. Must be 1 or 2.
                Ignored for inversion points, whose single chain always covers
                both branches. For a point without inversion the two directions
                are two independent chains, so 2 is a legitimate request even
                though ``fixed_point.num_branches`` is 1 (that property counts
                the branches of ONE chain -- plan 2.8).

        Returns:
            Dict[Tuple[int, int], BaseManifold]: All initial fundamental segments.
                The tuple structure is (orbit_index, branch_index).

        Raises:
            ValueError: If ``num_branches`` is not 1 or 2.

        Note:
            The chain is walked with :meth:`FixedPoint.advance_key`: ``+1`` for the
            unstable manifold and ``-1`` for the stable one. The stable fundamental
            segment is grown by the INVERSE map, which is why the orbit indices it
            visits run ``0, p-1, ..., 1`` (the same order ``get_iterable_array
            ("stable")`` produces).

        Note:
            For an inversion point the chain is now
            ``(0, 0), (1, 0), ..., (p-1, 0), (0, 1), ..., (p-1, 1)``. The previous
            implementation nested ``for orbit_index: for branch_index:`` and so
            visited ``(0, 1)`` immediately after ``(0, 0)``, which is physically
            wrong: the piece one map step after ``(0, 0)`` is ``(1, 0)``, and
            ``(0, 1)`` only comes a full orbit later, after ``(p-1, 0)``.
            Non-inversion points are unaffected -- the branch index never changes
            there.
        """

        if num_branches not in (1, 2):
            raise ValueError(f"num_branches must be 1 or 2, got {num_branches}")

        orbit_indices = fixed_point.get_iterable_array(stability)
        step = 1 if stability == "unstable" else -1

        if fixed_point.check_inversion():
            # num_branches is not consulted here: the single chain of k_value
            # pieces visits both branches whatever the caller asked for, and the
            # only value that differs from that outcome is the default 1 -- so
            # warning about it would fire on every ordinary call.
            branch_indices = fixed_point.get_branch_array()  # always [0, 1]

            initial_segments = {
                (orbit_index, branch_index): None
                for branch_index in branch_indices
                for orbit_index in orbit_indices
            }

            segment = self.get_initial_fundamental_segment(fixed_point, 0, 0, stability)
            initial_segments[(0, 0)] = segment
            segment.root = segment.walk_fwd(None, segment.root)

            key = (fixed_point, stability, 0, 0)
            for _ in range(fixed_point.k_value - 1):
                key = fixed_point.advance_key(key, step)
                orbit_index, branch_index = key[2], key[3]
                segment = self.machine.iterate_manifold(segment)
                # iterate_manifold builds the image without a branch identity (it
                # cannot advance the key itself); the chain walk knows it exactly.
                segment.manifold_key = key
                segment.branch_index = branch_index
                initial_segments[(orbit_index, branch_index)] = segment
                self.machine._insert_point_geometrically(
                    fixed_point.branch_points[orbit_index],
                    segment.root,
                    segment,
                    branch_index,
                )

        else:
            branch_indices = list(range(num_branches))

            initial_segments = {
                (orbit_index, branch_index): None
                for branch_index in branch_indices
                for orbit_index in orbit_indices
            }

            for b in branch_indices:
                segment = self.get_initial_fundamental_segment(
                    fixed_point, 0, b, stability
                )
                initial_segments[(0, b)] = segment
                segment.root = segment.walk_fwd(None, segment.root)

                key = (fixed_point, stability, 0, b)
                for _ in range(fixed_point.period - 1):
                    key = fixed_point.advance_key(key, step)
                    orbit_index, branch_index = key[2], key[3]
                    segment = self.machine.iterate_manifold(segment)
                    segment.manifold_key = key
                    segment.branch_index = branch_index
                    initial_segments[(orbit_index, branch_index)] = segment
                    self.machine._insert_point_geometrically(
                        fixed_point.branch_points[orbit_index],
                        segment.root,
                        segment,
                        branch_index,
                    )

        for dict_index in initial_segments:
            orbit_index = dict_index[0]
            initial_segments[dict_index].root = fixed_point.branch_points[orbit_index]

        return initial_segments
