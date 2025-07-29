from typing import Literal
from collections import deque

import numpy as np
from .BranchPoint import BranchPoint

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

    def check_inversion(self) -> bool:
        """
        Returns True if the Fixed Point has inversion
        False if it does not.
        """

        return True if self.period is not self.k_value else False

    def reset_accuracy(self):
        """
        Resets the accuracy based on the k value to account
        for mapping around the whole orbit.
        """

        self.accuracy = self.accuracy

    def set_k_value(self):
        """
        Sets the value 'k' which describes how many iterations it takes a nearby
        point to get back to that neighborhoods. If the fixed point has inversion
        this is double the period.
        """

        multiplier = 2 if any(x < 0 for x in self.unstable_eigenvalues) else 1

        self.k_value = self.period * multiplier

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
