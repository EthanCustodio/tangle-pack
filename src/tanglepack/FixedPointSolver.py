import numpy as np
from scipy.optimize import newton as newton_method
from scipy.differentiate import jacobian as jacob
from .FixedPoint import FixedPoint
from .DynamicalSystem import DynamicalSystem


class FixedPointSolver:

    def __init__(self, system: DynamicalSystem):

        self.dynamical_map = system.map
        self.dynamical_map_inverse = system.map_inv
        self.jacobian_function = system.jacobian

    def construct_fixed_point(self, initial_guess, num_branches):

        period, _ = np.shape(np.atleast_2d(initial_guess))

        fixed_point = self.compute_fixed_point(initial_guess)

        # TODO make this more correct for orbits with p > 1
        # shoot once for each point in the orbit
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
        point.reset_accuracy()

        return point

    def compute_fixed_point(self, initial_guess):
        """
        Computes a fixed point based on the initial guess
        """

        initial_guess_flattened = self.flatten_trajectory(initial_guess)

        fixed_point_flattened = newton_method(
            self.multipoint_shoot_flattened_difference, initial_guess_flattened
        )

        fixed_point_full = self.unflatten_trajectory(fixed_point_flattened)

        return np.array(fixed_point_full)

    def compute_eigenvectors(self, fixed_point, jacobians=None):
        """
        Computes the eigenvectors for each iterate of the fixed point
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

            print(f"eigenvalues: {eigenvalues}")
            eigenvalue_list[i][0] = eigenvalues[unstable_index]
            eigenvalue_list[i][1] = eigenvalues[1 - unstable_index]

        return eigenvalue_list, eigenvector_list

    def compute_partial_jacobians(self, fixed_point):
        """
        computes the partial step jacobians for each step of the
        fixed point
        """

        period, _ = np.shape(np.atleast_2d(fixed_point))

        difference = np.abs(self.multipoint_shoot(fixed_point) - fixed_point)
        initial_step = (np.average(np.linalg.norm(difference, axis=1))) ** (1 / 3)

        partial_jacobians = [np.empty((2, 2)) for i in range(period)]

        for i in range(period):

            # compute the single step jacobian at each iterate in the fixed point
            x_i = fixed_point[i]

            if self.jacobian_function is None:
                # if a jacobian function isn't provided, compute numerically
                jacobian = jacob(self.dynamical_map, x_i, initial_step=initial_step)
                jacobian = jacobian.df

            else:
                jacobian = self.jacobian_function(x_i)

            partial_jacobians[i] = jacobian

        return partial_jacobians

    def compute_jacobian(self, fixed_point):
        """
        Computes the Jacobian for the fixed point
        Computes p matrices for a period p fixed point
        This routine computes the Jacobian for each step and then
        returns the product of those matrices
        Will use a jacobian function if provided
        otherwise it will switch to a finite difference
        """

        # NOTE
        # this should now be computing the correct
        # full cycle jacobians at each point

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

    def multipoint_shoot_flattened_difference(self, trajectory):
        """
        Takes the difference between an iterate and the current trajectory

        Parameters:
            trajectory: list of points to map forward
        """

        shoot = self.multipoint_shoot_flattened(trajectory)

        difference = shoot - trajectory

        return difference

    def multipoint_shoot_flattened(self, trajectory):
        """
        Takes a flattened trajectory and iterates it forward

        Parameters:
            trajectory: list of points to map forward
        """

        trajectory_full = self.unflatten_trajectory(trajectory)

        trajectory_full_maped = self.multipoint_shoot(trajectory_full)

        trajectory_flattened_maped = self.flatten_trajectory(trajectory_full_maped)

        return trajectory_flattened_maped

    @staticmethod
    def flatten_trajectory(trajectory):
        """
        Takes a trajectory and makes it a 2n x 1 vector where n is the number of iterates

        Parameters:
            trajectory: list of points
        """

        number_iterates, _ = np.shape(np.atleast_2d(trajectory))

        trajectory_reshaped = np.reshape(trajectory, (2 * number_iterates, 1))

        return trajectory_reshaped

    @staticmethod
    def unflatten_trajectory(trajectory):
        """
        Takes a flattened trajectory and reformats it back into a n x d matrix

        Parameters:
            trajectory: list of points
        """

        number_iterates_doubled, _ = np.shape(np.atleast_2d(trajectory))

        number_iterates = number_iterates_doubled // 2

        trajectory_reshaped = np.reshape(trajectory, (number_iterates, 2))

        return trajectory_reshaped

    def multipoint_shoot(self, trajectory):
        """
        Takes a trajectory and dynamical_maps it forward once

        Parameters:
            trajectory: list of points
        """

        period, _ = np.shape(np.atleast_2d(trajectory))

        maped_trajectory = np.array(
            [self.dynamical_map(trajectory[i, :]) for i in range(period)]
        )

        maped_trajectory = np.roll(maped_trajectory, 1, axis=0)

        return maped_trajectory
