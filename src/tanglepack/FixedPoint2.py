import numpy as np
from scipy.optimize import newton as newton_method
from functools import partial
from scipy.differentiate import jacobian as jacob
from .Point import Point
from .MultiPoint import MultiPoint


class FixedPoint2(Point):

    def __init__(self, dynamical_map, initial_guess, dynamical_map_inverse=None, jacobian_function=None):
        """
        Class to store and compute fixed points for a given dynamical_map.

        Parameters:
            dynamical_map: function that defines the dynamical map
            jacobian_function: a user defined function with the jacobian of the dynamical map
            initial_guess: initial trajectory guess (n x 2 matrix)
            fixed_point: the computed fixed point
            accuracy: the difference between the fixed point and it's iterate
        """
        super().__init__(cdist=0.0, edist=0.0)

        self.dynamical_map = dynamical_map
        self.dynamical_map_inverse = dynamical_map_inverse
        self.jacobian_function = jacobian_function
        self.initial_guess = initial_guess

        self.period, _ = np.shape(np.atleast_2d(self.initial_guess))

        self.fixed_point = self.compute_fixed_point()

        print(f'wow fixed point: {self.fixed_point}')

        # TODO comment on how the accuracy is computed
        self.difference = np.abs(self.multipoint_shoot(self.fixed_point) - self.fixed_point)
        self.accuracy = (np.average(np.linalg.norm(self.difference, axis=1))) ** (1/3)
        print(f"Fixed Point Accurarcy: {self.accuracy}")

        self.jacobians = [np.empty((2, 2)) for i in range(self.period)]
        self.jacobian = self.compute_jacobian()

        # eigenvectors are ordered: [unstable, stable]
        self.eigenvectors = [[np.empty((2, 1)), np.empty((2, 1))] for i in range(self.period)]
        self.eigenvalues = [np.empty((2, 1)) for i in range(self.period)]
        self.compute_eigenvectors()

        self.inversion_factor = self.compute_inversion_factor()

        self.unstable_directions = [np.empty((2, 1)) for i in range(self.period * self.inversion_factor)]
        self.stable_directions = [np.empty((2, 1)) for i in range(self.period * self.inversion_factor)]
        self.compute_directions()

        # declare all the MultiPoint objects for the ring structure
        self.set_points()

        self.get_first_fundamental_segments()


    def get_second_point(self, period, stability):
        """
        Computes the second point from the fixed point based on a linear interpolation
        """

        step = self.accuracy

        if stability == 'unstable':
            second_point = self.fixed_point[period] + (step) * self.unstable_directions[period]

        if stability == 'stable':
            second_point = self.fixed_point[period] + (step) * self.stable_directions[period]

        return np.array(second_point, dtype=np.float64).reshape(-1)
    
    
    def get_first_point(self, point, stability):
        """
        Computes the first point from the fixed point based on preiterating the second
        """

        first_point = np.array([point.x, point.y], dtype=np.float64)

        for i in range(self.period * self.inversion_factor):

            if stability == 'unstable':
                first_point = self.dynamical_map_inverse(first_point)

            if stability == 'stable':
                first_point = self.dynamical_map(first_point)

        return np.array(first_point, dtype=np.float64).reshape(-1)


    def get_first_fundamental_segments(self):
        """
        Computes the initial fundamental segments made of two points
        """

        for period in range(self.period):
            for branch in range(self.inversion_factor):
            # TODO watch out for this branch logic. It is likely not in iterate order

                second_point_unstable = self.get_second_point(period, 'unstable')
                second_point_unstable = Point(second_point_unstable[0], second_point_unstable[1])

                second_point_stable = self.get_second_point(period, 'stable')
                second_point_stable = Point(second_point_stable[0], second_point_stable[1])

                first_point_unstable = self.get_first_point(second_point_unstable, 'unstable')
                first_point_unstable = Point(first_point_unstable[0], first_point_unstable[1])

                first_point_stable = self.get_first_point(second_point_stable, 'stable')
                first_point_stable = Point(first_point_stable[0], first_point_stable[1])

                # next_manifolds is the unstable directions
                second_point_unstable.insert_iterate_after(first_point_unstable)
                second_point_unstable.insert_manifold_after(first_point_unstable)
                self.points[period].next_manifolds[branch] = first_point_unstable

                # prev_manifolds is the stable direction
                second_point_stable.insert_iterate_before(first_point_stable)
                second_point_stable.insert_manifold_after(first_point_stable)
                self.points[period].prev_manifolds[branch] = first_point_stable

                # update the stretch factor for each step in the orbit
                self.points[period].unstable_stretch_params[branch] = (
                        self.compute_canonical_distance(first_point_unstable, second_point_unstable, period, 'unstable')
                    )
                self.points[period].stable_stretch_params[branch] = (
                        self.compute_canonical_distance(first_point_stable, second_point_stable, period, 'stable')
                    )

                
    def compute_edist(self, point, period):
        """
        Computes the euclidean distance from the fixed point to the point
        """

        first_point_x = self.points[period].x
        first_point_y = self.points[period].y

        second_point_x = point.x
        second_point_y = point.y

        x_diff = np.abs(first_point_x - second_point_x)
        y_diff = np.abs(first_point_y - second_point_y)

        distance = np.sqrt(x_diff ** 2 + y_diff ** 2)

        point.edist = distance

        return distance


    def compute_canonical_distance(self, point1:Point, point2:Point, period:int, stability:str):
        """
        Computes the canonical distance for two iterates in the linear regime.
        Updates the alpha values for the fixed point in the process
        """

        edist1 = self.compute_edist(point1, period)
        edist2 = self.compute_edist(point2, period)

        if edist2 > edist1 and stability == "unstable":
            stretch_param = edist2 / edist1
        else:
            stretch_param = edist1 / edist2

        point1.cdist = edist1
        point2.cdist = edist2

        return stretch_param


    def compute_directions(self):
        """
        Computes the directions of all manifolds attached to the fixed point
        """

        # TODO make this generalize to k > 1 inversion factor
        # The problem is that the other directions we need aren't listed in the eigenvector object
        for i in range(self.period * self.inversion_factor):
            
            self.unstable_directions[i] = self.eigenvectors[i][0]
            self.stable_directions[i] = self.eigenvectors[i][1]


    def compute_inversion_factor(self):
        """
        Computes the inversion factor, k, for the fixed point
        """
        
        return 1


    def set_points(self):
        """
        Takes the computed fixed point and creates a list of 
            point objects for them
        This also creates a loop in the graph/linked list structure
        """
        x = self.fixed_point[:, 0]
        y = self.fixed_point[:, 1]

        # MultiPoint has next for unstable and prev for stable
        manifolds_per_point = self.inversion_factor

        self.points = [MultiPoint(x=x[i], y=y[i], num_manifolds=manifolds_per_point, cdist=0.0, edist=0.0) for i in range(self.period)]

        points_shifted_right = [self.points[-1]] + self.points[:-1]
        points_shifted_left = self.points[1:] + [self.points[0]]

        for point in points_shifted_right: 
            point.next_iterate = point
    
        for point in points_shifted_left:
            point.prev_iterate = point


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
        Computes the eigenvectors for each iterate of the fixed point
        """

        print("entered eigenvalues")

        initial_step = self.accuracy

        jacobian = self.jacobians[0]

        eigenvalues, eigenvectors = np.linalg.eig(jacobian)

        # largest eigenvalue is unstable direction
        unstable_index = np.argmax(np.abs(eigenvalues))

        # WARNING POORLY UNDERSTOOD why we multiply by -1. This is henon specific
        # TODO find a way to automatically choose the proper directions
        self.eigenvectors[0][0] = -1*eigenvectors[:, unstable_index]
        self.eigenvalues[0][0] = eigenvalues[unstable_index]
        self.eigenvectors[0][1] = eigenvectors[:, 1 - unstable_index]
        self.eigenvalues[0][1] = eigenvalues[1 - unstable_index]

        for i in range(0, self.period - 1):

            jacobian = self.jacobians[i]

            # TODO think if we need the eigenvalues too
            self.eigenvectors[i + 1][0] = jacobian @ self.eigenvectors[i][0]
            self.eigenvectors[i + 1][1] = jacobian @ self.eigenvectors[i][1]


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
            self.jacobians[i] = jacobian
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