import numpy as np
from scipy.optimize import newton as newton_method
from scipy.differentiate import jacobian as jacob
from .FixedPoint import FixedPoint


class FixedPointSolver():


    def __init__(self, dynamical_map, dynamical_map_inverse=None, jacobian_function=None):

        
        self.dynamical_map = dynamical_map
        self.dynamical_map_inverse = dynamical_map_inverse
        self.jacobian_function = jacobian_function


    def construct_fixed_point(self, initial_guess, num_branches):

        period, _ = np.shape(np.atleast_2d(initial_guess))

        fixed_point = self.compute_fixed_point(initial_guess)

        difference = np.abs(self.multipoint_shoot(fixed_point) - fixed_point)
        accuracy = (np.average(np.linalg.norm(difference, axis=1))) ** (1/3)

        eigenvalues, eigenvectors = self.compute_eigenvectors(fixed_point)
        jacobians, total_jacobian = self.compute_jacobian(fixed_point)

        point = FixedPoint(period, num_branches)

        point.coordinates = fixed_point
        point.accuracy = accuracy
        point.total_jacobian = total_jacobian

        for i in range(period):

            point.unstable_eigenvectors[i] = eigenvectors[i][0]
            point.stable_eigenvectors[i] = eigenvectors[i][1]
            
            point.unstable_eigenvalues[i] = eigenvalues[i][0]
            point.stable_eigenvalues[i] = eigenvalues[i][1]

            point.jacobians[i] = jacobians[i]

        return point


    def compute_fixed_point(self, initial_guess):
        """
        Computes a fixed point based on the initial guess
        """

        initial_guess_flattened = self.flatten_trajectory(initial_guess)

        fixed_point_flattened = newton_method(self.multipoint_shoot_flattened_difference, initial_guess_flattened)

        fixed_point_full = self.unflatten_trajectory(fixed_point_flattened)

        print(f"This is the fixed point {fixed_point_full}")

        return np.array(fixed_point_full)
    
    
    def compute_eigenvectors(self, fixed_point):
        """
        Computes the eigenvectors for each iterate of the fixed point
        """

        period, _ = np.shape(np.atleast_2d(fixed_point))

        eigenvector_list = [[np.empty((2, 1)), np.empty((2, 1))] for i in range(period)]
        eigenvalue_list = [np.empty((2, 1)) for i in range(period)]

        difference = np.abs(self.multipoint_shoot(fixed_point) - fixed_point)
        initial_step = (np.average(np.linalg.norm(difference, axis=1))) ** (1/3)

        for i in range(period):

            x_i = fixed_point[i]

            if self.jacobian_function is None:
                jacobian = jacob(self.dynamical_map, x_i, initial_step=initial_step)
                jacobian = jacobian.df

            else:
                jacobian = self.jacobian_function(x_i)

            eigenvalues, eigenvectors = np.linalg.eig(jacobian)

            unstable_index = np.argmax(np.abs(eigenvalues))

            # WARNING POORLY UNDERSTOOD why we multiply by -1. This is henon specific
            # TODO find a way to automatically choose the proper directions
            eigenvector_list[i][0] = -1*eigenvectors[:, unstable_index].reshape(2, 1)
            eigenvector_list[i][1] = eigenvectors[:, 1 - unstable_index].reshape(2, 1)

            eigenvalue_list[i][0] = eigenvalues[unstable_index]
            eigenvalue_list[i][1] = eigenvalues[1 - unstable_index]

        return eigenvalue_list, eigenvector_list


    def compute_jacobian(self, fixed_point):
        """
        Computes the Jacobian for the fixed point
        Computes p matrices for a period p fixed point
        This routine computes the Jacobian for each step and then
        returns the product of those matrices
        Will use a jacobian function if provided
        otherwise it will switch to a finite difference
        """

        period, _ = np.shape(np.atleast_2d(fixed_point))

        difference = np.abs(self.multipoint_shoot(fixed_point) - fixed_point)
        initial_step = (np.average(np.linalg.norm(difference, axis=1))) ** (1/3)

        jacobians = [np.empty((2, 2)) for i in range(period)]
        total_jacobian = np.identity(2)

        for i in range(period):
            # compute the jacobian at each iterate in the fixed point

            x_i = fixed_point[i] 

            if self.jacobian_function is None:
                # if a jacobian function isn't provided, compute numerically
                jacobian = jacob(self.dynamical_map, x_i, initial_step=initial_step)
                jacobian = jacobian.df

            else:
                jacobian = self.jacobian_function(x_i)

            jacobians[i] = jacobian
            total_jacobian = jacobian @ total_jacobian

        return jacobians, total_jacobian


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

        maped_trajectory = np.array([self.dynamical_map(trajectory[i, :]) for i in range(period)])

        return maped_trajectory
    