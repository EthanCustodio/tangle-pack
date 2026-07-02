import logging
from typing import Optional

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import fsolve
from scipy.differentiate import jacobian as jacob
from .FixedPoint import FixedPoint
from .DynamicalSystem import DynamicalSystem

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

"""
Dev Notes:

num_branches. There should be a way to remove this entirely from the 
contruction method. That information is contained within the k_value
Though we may want to keep it if someone wants to compute both sides. 
We would need to think through that case more explicitly, it might
be detrimental to have it in like this actually.
"""


class FixedPointSolver:
    """
    Toolbox for computing fixed points of any period and initializing
    FixedPoint objects.

    Attributes:
        dynamical_map (func): Function for the map of the system.
        dynamical_map_inverse (func): Function for the inverse map of the system.
        jacobian_function (func, optional): Function to compute the Jacobian at
            any point in the system. Without an explicit Jacobian function finite
            difference will be used to compute the Jacobians.
    """

    def __init__(self, system: DynamicalSystem) -> None:
        """
        Initializes the functions associated with the dynamical system.

        Args:
            system (DynamicalSystem): Object storing all the system maps.
        """

        self.dynamical_map = system.map
        self.dynamical_map_inverse = system.map_inv
        self.jacobian_function = system.jacobian

    def construct_fixed_point(
        self, initial_guess: NDArray[np.float64], num_branches: int
    ) -> FixedPoint:
        """
        Computes the fixed point from an initial guess using a multipoint
        shooting Newton's method and initializes a FixedPoint object containing
        all the information.

        Args:
            initial_guess (np.ndarray): A (period, 2) array.
                Each row is an initial guess for one iterate.
            num_branches (int): Number of branches the fixed point has.
                Based on inversion.

        Returns:
            FixedPoint: The fully constructed fixed point.
        """

        period, _ = np.shape(np.atleast_2d(initial_guess))

        fixed_point = self.compute_fixed_point(initial_guess)

        difference = np.abs(self.multipoint_shoot(fixed_point) - fixed_point)
        accuracy = (np.average(np.linalg.norm(difference, axis=1))) ** (1 / 3)

        jacobians = self.compute_jacobian(fixed_point)
        partial_jacobians = self.compute_partial_jacobians(fixed_point)

        eigenvalues, eigenvectors = self.compute_eigenvectors(fixed_point, jacobians)

        point = FixedPoint(period, num_branches)

        point.coordinates = fixed_point
        point.accuracy = accuracy

        for i in range(period):

            p_x, p_y = point.coordinates[i]
            point.branch_points[i].set_x(p_x)
            point.branch_points[i].set_y(p_y)

            point.branch_points[i].cdist = 0.0

            point.unstable_eigenvectors[i] = eigenvectors[i][0]
            point.stable_eigenvectors[i] = eigenvectors[i][1]

            point.unstable_eigenvalues[i] = eigenvalues[i][0]
            point.stable_eigenvalues[i] = eigenvalues[i][1]

            point.jacobians[i] = jacobians[i]
            point.partial_jacobians[i] = partial_jacobians[i]

            point.branch_points[i].next_iterate = point.branch_points[(i + 1) % period]
            point.branch_points[i].prev_iterate = point.branch_points[(i - 1) % period]

        point.set_k_value()

        return point

    def compute_fixed_point(
        self, initial_guess: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """
        Computes the coordinates of a fixed point based on the initial guess. Uses
        a multipoint shooting Newton's method.

        Args:
            initial_guess (np.ndarray): A (period, 2) array representing an
                initial guess for the orbit.

        Returns:
            np.ndarray: A (period, 2) array representing the converged fixed point
                coordinates.
        """

        initial_guess_flattened = self.flatten_trajectory(initial_guess)

        # Multipoint shooting couples every coordinate (residual_i = map(x_{i-1}) -
        # x_i), so this is a multivariate root-find. ``scipy.optimize.newton`` on an
        # array solves each coordinate as an INDEPENDENT scalar problem, which only
        # "converges" when the guess is already the answer (e.g. the exact orbit at
        # k=2) and silently returns a non-orbit otherwise -- the source of the
        # k-sensitive, eigenvector-flip-requiring behaviour. ``fsolve`` (MINPACK
        # hybrd) solves the coupled system properly and, staying near the guess,
        # also keeps the orbit labelling consistent with ``initial_guess``.
        shape = np.shape(initial_guess_flattened)
        fixed_point_flattened = fsolve(
            lambda x: np.ravel(
                self.multipoint_shoot_flattened_difference(np.reshape(x, shape))
            ),
            np.ravel(initial_guess_flattened),
            xtol=1e-13,
            maxfev=10000,
        ).reshape(shape)

        fixed_point_full = self.unflatten_trajectory(fixed_point_flattened)

        return np.array(fixed_point_full)

    def compute_eigenvectors(
        self,
        fixed_point: NDArray[np.float64],
        jacobians: Optional[list[NDArray[np.float64]]] = None,
    ) -> tuple[list[NDArray[np.float64]], list[list[NDArray[np.float64]]]]:
        """
        Computes the eigenvectors for each iterate of the fixed point.

        Args:
            fixed_point (np.ndarray): (period, 2) array of orbit coordinates.
            jacobians (list of np.ndarray, optional): Full-cycle Jacobians at
                each iterate. Automatically computes the Jacobians if not provided.

        Returns:
            tuple: (eigenvalues, eigenvectors), where:
                eigenvalues: A list of (2, 1) arrays;
                    [unstable, stable] for each iterate.
                eigenvectors: A list of two (2, 1) arrays;
                    [unstable, stable] for each iterate.
        """

        if jacobians is None:
            jacobians = self.compute_jacobian(fixed_point)

        period, _ = np.shape(np.atleast_2d(fixed_point))

        eigenvector_list = [[np.empty((2, 1)), np.empty((2, 1))] for _ in range(period)]
        eigenvalue_list = [np.empty((2, 1)) for _ in range(period)]

        for i in range(period):

            jacobian = jacobians[i]

            eigenvalues, eigenvectors = np.linalg.eig(jacobian)

            unstable_index = np.argmax(np.abs(eigenvalues))

            # WARNING POORLY UNDERSTOOD why we multiply by -1. This is henon specific
            # TODO find a way to automatically choose the proper directions
            eigenvector_list[i][0] = eigenvectors[:, unstable_index].reshape(2, 1)
            eigenvector_list[i][1] = eigenvectors[:, 1 - unstable_index].reshape(2, 1)

            logger.debug("eigenvalues at orbit index %d: %s", i, eigenvalues)
            eigenvalue_list[i][0] = eigenvalues[unstable_index]
            eigenvalue_list[i][1] = eigenvalues[1 - unstable_index]

        return eigenvalue_list, eigenvector_list

    def compute_partial_jacobians(
        self, fixed_point: NDArray[np.float64]
    ) -> list[NDArray[np.float64]]:
        """
        Computes the partial step jacobians for each step of the
        fixed point.

        Args:
            fixed_point (np.ndarray): A (period, 2) array representing the fixed point.

        Returns:
            list[np.ndarray]: List of (2, 2) Jacobians for single steps at each iterate.
        """

        period, _ = np.shape(np.atleast_2d(fixed_point))

        difference = np.abs(self.multipoint_shoot(fixed_point) - fixed_point)
        initial_step = (np.average(np.linalg.norm(difference, axis=1))) ** (1 / 3)

        partial_jacobians = [np.empty((2, 2)) for i in range(period)]

        for i in range(period):

            # compute the single step jacobian at each iterate in the fixed point
            x_i = fixed_point[i]

            if self.jacobian_function is None:
                # if a jacobian function isn't provided, compute numerically. The
                # finite-difference step reuses the orbit residual^(1/3); when the
                # orbit is located exactly (residual 0) that is 0, which makes the
                # difference stencil degenerate (nan). Fall back to scipy's own
                # default step in that case.
                if initial_step > 0:
                    jacobian = jacob(self.dynamical_map, x_i, initial_step=initial_step)
                else:
                    jacobian = jacob(self.dynamical_map, x_i)
                jacobian = jacobian.df

            else:
                jacobian = self.jacobian_function(x_i)

            partial_jacobians[i] = jacobian

        return partial_jacobians

    def compute_jacobian(
        self, fixed_point: NDArray[np.float64]
    ) -> list[NDArray[np.float64]]:
        """
        Computes the Jacobian for the fixed point.
        Computes p matrices for a period p fixed point.
        This routine computes the Jacobian for each step and then
        returns the product of those matrices for the full cycle jacobians.
        Will use a jacobian function if provided
        otherwise it will switch to a finite difference.

        Args:
            fixed_point (np.ndarray): A (period, 2) array representing the fixed point.

        Returns:
            list[np.ndarray]: List of (2, 2) full-cycle Jacobians, one per iterate.
        """

        period, _ = np.shape(np.atleast_2d(fixed_point))

        partial_jacobians = self.compute_partial_jacobians(fixed_point)

        true_jacobians = []
        # compute each cyclic permutation to get the full cycle jacobians
        for shift in range(period):

            factors = [partial_jacobians[(shift + i) % period] for i in range(period)]
            factors = factors[::-1]

            if period == 1:
                product = factors[0].copy()
            else:
                product = np.linalg.multi_dot(factors)
            true_jacobians.append(product)

        return true_jacobians

    def multipoint_shoot_flattened_difference(
        self, trajectory: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """
        Takes the difference between an iterate and the current trajectory.

        Args:
            trajectory (np.ndarray): (2 * period, 1) array of points to map forward.

        Returns:
            np.ndarray: (2 * period, 1) array.
        """

        shoot = self.multipoint_shoot_flattened(trajectory)

        difference = shoot - trajectory

        return difference

    def multipoint_shoot_flattened(
        self, trajectory: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """
        Takes a flattened trajectory and iterates it forward.

        Args:
            trajectory (np.ndarray): (2 * period, 1) array of points to shoot forward.

        Returns:
            np.ndarray: (2 * period, 1) array of the mapped trajectory.
        """

        trajectory_full = self.unflatten_trajectory(trajectory)

        trajectory_full_mapped = self.multipoint_shoot(trajectory_full)

        trajectory_flattened_mapped = self.flatten_trajectory(trajectory_full_mapped)

        return trajectory_flattened_mapped

    @staticmethod
    def flatten_trajectory(trajectory: NDArray[np.float64]) -> NDArray[np.float64]:
        """
        Takes a trajectory and makes it a 2n x 1 vector where n is
        the number of iterates.

        Args:
            trajectory (np.ndarray): (period, 2) array of points representing
                a trajectory.

        Returns
            np.ndarray: (2 * period, 1) array of points representing a
                flattened trajectory.
        """

        number_iterates, _ = np.shape(np.atleast_2d(trajectory))

        trajectory_reshaped = np.reshape(trajectory, (2 * number_iterates, 1))

        return trajectory_reshaped

    @staticmethod
    def unflatten_trajectory(trajectory: NDArray[np.float64]) -> NDArray[np.float64]:
        """
        Takes a flattened trajectory and reformats it back into a p x 2 matrix.

        Args:
            trajectory (np.ndarray): (2 * period, 1) array of points representing a
                flattened trajectory.

        Returns:
            np.ndarray: (period, 2) array of points representing an
                unflattened trajectory.
        """

        number_iterates_doubled, _ = np.shape(np.atleast_2d(trajectory))

        number_iterates = number_iterates_doubled // 2

        trajectory_reshaped = np.reshape(trajectory, (number_iterates, 2))

        return trajectory_reshaped

    def multipoint_shoot(self, trajectory: NDArray[np.float64]) -> NDArray[np.float64]:
        """
        Takes a trajectory and dynamical_maps it forward once.

        Args:
            trajectory (np.ndarray): (period, 2) array of points representing
                a trajectory.

        Returns:
            np.ndarray: (period, 2) array of the resulting trajectory.
        """

        period, _ = np.shape(np.atleast_2d(trajectory))

        maped_trajectory = np.array(
            [self.dynamical_map(trajectory[i, :]) for i in range(period)]
        )

        maped_trajectory = np.roll(maped_trajectory, 1, axis=0)

        return maped_trajectory
