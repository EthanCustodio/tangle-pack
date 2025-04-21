import numpy as np
from scipy.optimize import newton as newton_method
from functools import partial
from scipy.differentiate import jacobian as jacob


class FixedPoint():

    def __init__(self, dynamical_map, initial_guess, dynamical_map_inverse=None, jacobian_function=None):
        #dynamical_map
        """
        Class to store and compute fixed points for a given dynamical_map.

        Parameters:
            dynamical_map: function that defines the dynamical map
            jacobian_function: a user defined function with the jacobian of the dynamical map
            initial_guess: initial trajectory guess (n x 2 matrix)
            fixed_point: the computed fixed point
            accuracy: the difference between the fixed point and it's iterate
        """

        self.dynamical_map = dynamical_map
        self.dynamical_map_inverse = dynamical_map_inverse
        self.jacobian_function = jacobian_function
        self.initial_guess = initial_guess
        self.period, _ = np.shape(np.atleast_2d(self.initial_guess))

        self.fixed_point = self.compute_fixed_point()

        self.difference = np.abs(self.multipoint_shoot(self.fixed_point) - self.fixed_point)
        self.accuracy = (np.average(np.linalg.norm(self.difference, axis=1))) ** (1/3)
        print(f"Fixed Point Accurarcy: {self.accuracy}")

        # eigenvectors are ordered: [unstable, stable]
        self.jacobian = self.compute_jacobian()
        self.eigenvectors = [[np.empty((2, 1)), np.empty((2, 1))] for i in range(self.period)]
        self.eigenvalues = [np.empty((2, 1)) for i in range(self.period)]
        self.compute_eigenvectors()
        self.stable_direction = None
        self.unstable_direction = None


    def compute_fixed_point(self):
        """
        Computes a fixed point based on the initial guess
        """

        initial_guess_flattened = self.flatten_trajectory(self.initial_guess)

        fixed_point_flattened = newton_method(self.multipoint_shoot_flattened_difference, initial_guess_flattened)

        fixed_point_full = self.unflatten_trajectory(fixed_point_flattened)

        print(f"This is the fixed point {fixed_point_full}")

        return np.array(fixed_point_full)
    
    
    def compute_eigenvectors(self):
        """
        Comptues the eigenvectors for each iterate of the fixed point
        """

        print("entered eigenvalues")

        initial_step = self.accuracy

        jacobian = self.jacobian

        eigenvalues, eigenvectors = np.linalg.eig(jacobian)

        # largest eigenvalue is unstable direction
        unstable_index = np.argmax(np.abs(eigenvalues))


        # WARNING POORLY UNDERSTOOD why we multiply by -1. This is henon specific
        # TODO find a way to automatically choose the proper directions
        self.eigenvectors[0][0] = -1*eigenvectors[:, unstable_index]
        self.eigenvalues[0][0] = eigenvalues[unstable_index]
        self.eigenvectors[0][1] = eigenvectors[:, 1 - unstable_index]
        self.eigenvalues[0][1] = eigenvalues[1 - unstable_index]


        for i in range(1, self.period):

            x_i = self.fixed_point[i]

            if self.jacobian_function is None:
                jacobian = jacob(self.dynamical_map, x_i, initial_step=initial_step)
                jacobian = jacobian.df

            else:
                jacobian = self.jacobian_function(x_i)

            self.eigenvectors[i][0] = jacobian @ eigenvectors[:, unstable_index]
            self.eigenvectors[i][1] = jacobian @ eigenvectors[:, 1 - unstable_index]


    def compute_jacobian(self):
        """
        Computes the Jacobian for the fixed point
        Computes p matrices for a period p fixed point
        This routine computes the Jacobian for each step and then
        returns the product of those matrices
        Will use a jacobian function if provided
        otherwise it will switch to a finite difference
        """

        # initial_step = sum(self.accuracy) / len(self.accuracy)
        initial_step = self.accuracy

        jacobians = [np.empty((2, 2)) for i in range(self.period)]
        total_jacobian = np.identity(2)

        for i in range(self.period):
            # compute the jacobian at each iterate in the fixed point

            x_i = self.fixed_point[i] 

            if self.jacobian_function is None:
                # if a jacobian function isn't provided, compute numerically
                jacobian = jacob(self.dynamical_map, x_i, initial_step=initial_step)
                jacobian = jacobian.df

            else:
                jacobian = self.jacobian_function(x_i)

            jacobians[i] = jacobian
            total_jacobian = jacobian @ total_jacobian

        return total_jacobian


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