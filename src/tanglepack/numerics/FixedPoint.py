from __future__ import annotations

from typing import Literal, Optional, Tuple
from collections import deque

import numpy as np
from .BranchPoint import BranchPoint

ManifoldKey = Tuple["FixedPoint", Literal["unstable", "stable"], int, int]

"""
Dev Notes:

Ensure that the shape of the coordinates array is consisntent everywhere
in the package. In this class it is (2, 1). I know in BasePoint it is (2,).

Check how the unstable eigenvectors is implemented when it comes to points 
with inversion. We may or may not need to store two vectors, the method
we use elsewhere I think just takes the initial vector and maps it around
k times so we don't need to be able to walk in both eigendirections likely.

Specify why we store the partial_jacobians. What are they used for?
"""


class FixedPoint:
    """
    Implements the data structure to store all information about a fixed point.

    Note:
        Does not include the methods to compute nor construct fixed points.
        Those are located in FixedPointSolver.py.

    Attributes:

        period (int): Period of the fixed point.
        num_branches (int): Number of branches attached to the fixed point for each
            manifold. One if the point has no inversion, two if it has inversion.
        branch_points (List[BranchPoint]): List of the individual BranchPoints that
            make up the fixed point.
        coordinates (List[np.ndarray]): Array storing the coordinates of the
            fixed point. The shape is (2, 1).

        accuracy (float): Computed as the norm of the difference between the
            fixed point and it's iterate. {latex}
        k_value (int): Number of iterations it takes to return to the segment
            of manifold you started from.

        unstable_eigenvectors (List[np.ndarray]): List of arrays storing the
            eigenvectors corresponding to the direction of the unstable manifold.
        unstable_eigenvalues (List[float]): List of arrays storing the eigenvalues
            corresponding to the unstable manifold. If an eigenvalue is < 0 the
            fixed point has inversion.
        stable_eigenvectors (List[np.ndarray]): List of arrays storing the eigenvectors
            corresponding to the direction fo the stable manifold.
        stable_eigenvalues (List[float]): List of arrays storing the eigenvalues
            corresponding to the stable manifold. If an eigenvalue is < 0 the
            fixed point has inversion.

        jacobians (List[np.ndarray]): List of the full cycle Jacobians at each
            iterate of the fixed point. These Jacobians are used to compute the
            eigenstuffs.
        partial_jacobians (List[np.ndarray]): List of the single step Jacobians
            at each iterate of the fixed point. These Jacobians are used to
    """

    def __init__(self, period: int, num_branches: int) -> None:
        """
        Allocates the data structures which store fixed point information.

        Note:
            Does not include the methods to compute nor construct fixed points.
            Those are located in FixedPointSolver.py.

        Args:
            period (int): Period of the fixed point.
            num_branches (int): Number of branches attached to the fixed point for each
                manifold. One if the point has no inversion, two if it has inversion.
        """

        self.period = period
        self.num_branches = num_branches

        self.branch_points = [
            BranchPoint(num_branches, (0.0, 0.0)) for _ in range(period)
        ]
        self.coordinates = [np.empty((2, 1)) for _ in range(period)]

        self.unstable_eigenvectors = [np.empty((2, 1)) for _ in range(period)]
        self.unstable_eigenvalues = [0.0] * period

        self.stable_eigenvectors = [np.empty((2, 1)) for _ in range(period)]
        self.stable_eigenvalues = [0.0] * period

        self.accuracy = 0.0

        self.jacobians = [np.empty((2, 2)) for _ in range(period)]
        self.partial_jacobians = [np.empty((2, 2)) for _ in range(period)]

        # Filled in by set_k_value() once the eigenvalues are known.
        self.k_value: Optional[int] = None

    def check_inversion(self) -> bool:
        """
        Returns True if the Fixed Point has inversion
        False if it does not.

        Raises:
            ValueError: If called before set_k_value().
        """

        if self.k_value is None:
            raise ValueError(
                "k_value has not been set yet; call set_k_value() after the "
                "eigenvalues are computed."
            )

        return self.period != self.k_value

    def set_k_value(self):
        """
        Sets the value 'k' which describes how many iterations it takes a nearby
        point to get back to that neighborhoods. If the fixed point has inversion
        this is double the period.
        """

        multiplier = 2 if any(x < 0 for x in self.unstable_eigenvalues) else 1

        self.k_value = self.period * multiplier

    def advance_key(self, key: ManifoldKey, n: int = 1) -> ManifoldKey:
        """
        Advance a manifold key by ``n`` applications of the map.

        This is the single source of truth for "which manifold piece does one map
        step send this one to". One application of M sends the manifold piece
        attached to orbit point ``i`` onto the piece at orbit point ``(i + 1) mod
        period``, for BOTH stabilities: M maps W^s(x_i) onto W^s(x_{i+1}) exactly
        as it maps W^u(x_i) onto W^u(x_{i+1}). ``stability`` is therefore carried
        through untouched -- it is the caller who decides whether it wants the
        forward step (``n = +1``, how the unstable fundamental segment is grown)
        or the backward one (``n = -1``, how the stable one is).

        The branch index labels the eigenvector direction at orbit point 0; the
        pieces at the other orbit points inherit it by propagation. Without
        inversion the branch never changes. With inversion (negative eigenvalue,
        ``k_value == 2 * period``) M^period returns to orbit point 0 on the
        OPPOSITE branch, so the branch flips exactly when the orbit index wraps
        (period-1 -> 0 going forward, 0 -> period-1 going backward) and two full
        orbits -- ``k_value`` steps -- return to the start.

        The map is a group action on the ``k_value`` pieces: ``advance_key`` is
        additive in ``n`` and has period ``k_value``.

        Args:
            key: A manifold key ``(fixed_point, stability, orbit_index,
                branch_index)``. Its fixed point must be this one.
            n: Number of map steps. Negative means applications of M^-1.
                Defaults to 1.

        Returns:
            ManifoldKey: The key ``n`` map steps along the chain.

        Raises:
            ValueError: If the key belongs to another fixed point, if
                ``set_k_value()`` has not been called, or if the branch index is
                outside the fixed point's branch range.

        Note:
            The ordering this induces on an inversion point is
            ``(0, 0), (1, 0), ..., (p-1, 0), (0, 1), ..., (p-1, 1)``: the piece
            after ``(p-1, 0)`` is ``(0, 1)``. The pre-2026-09 initializer instead
            walked ``for orbit_index: for branch_index:``, visiting ``(0, 1)``
            before ``(1, 0)``, which is not the physical chain.
        """

        fp, stability, orbit_index, branch_index = key

        if fp is not self:
            raise ValueError(
                "advance_key was given a key belonging to a different fixed point."
            )

        if self.k_value is None:
            raise ValueError(
                "k_value has not been set yet; call set_k_value() after the "
                "eigenvalues are computed."
            )

        if not 0 <= branch_index < self.num_branches:
            raise ValueError(
                f"branch_index {branch_index} is outside this fixed point's branch "
                f"range (num_branches={self.num_branches})."
            )

        if self.check_inversion():
            assert self.num_branches == 2, (
                "an inversion point must have two branches for the chain of "
                f"k_value pieces to close; got num_branches={self.num_branches}"
            )
            # The k_value pieces are laid out as orbit_index + period * branch_index,
            # so one map step is +1 on that single cyclic index and the branch flip
            # falls out of the wrap.
            position = (orbit_index + self.period * branch_index + n) % self.k_value
            return (fp, stability, position % self.period, position // self.period)

        return (fp, stability, (orbit_index + n) % self.period, branch_index)

    def get_iterable_array(
        self, stability: Literal["unstable", "stable"], shift: int = 0
    ) -> list[int]:
        """
        Construct a list of orbit indices based off the stability.
        Starting at the zeroth index this is the order a manifold will move
        between the indices of the fixed points

        Ex:
            unstable [0, 1, 2, 3]
            stable [0, 3, 2, 1]

        Parameters:
            stability: The stability of the manifold.
            shift: Amount to cyclically shift the resulting array
                    "takes the last element and moves it to the front".
                    This is equivalent to starting at another point in the cycle.
        """

        orbit_indices = [i for i in range(self.period)]
        if stability == "stable":
            temp = orbit_indices[1:]
            orbit_indices = [0] + temp[::-1]

        orbit_indices = deque(orbit_indices)
        orbit_indices.rotate(shift)

        return list(orbit_indices)

    def get_branch_array(self) -> list[int]:
        """
        Get an array to iterate over the branch indices.

        Returns:
            list[int]: A list of branch indices [0] or [0, 1],
                depending on whether the fixed point has inversion.
        """

        # determine if there is inversion
        if self.check_inversion():
            branch_indices = [0, 1]
        else:
            branch_indices = [0]

        return branch_indices
