import numpy as np
import scipy.integrate as spi
from .Points import Points
from .FixedPoint import FixedPoint
from numpy.linalg import LinAlgError


class Manifold():

    def __init__(self, fixed_point: FixedPoint = None):
        """
        Low level object to store manifold data and routines

        Parameters:
            fixed_point (FixedPoint): fixed_point that the manifold is attached to
        """

        self.fixed_point = fixed_point
        self.dynamical_map = fixed_point.dynamical_map
        self.dynamical_map_inverse = fixed_point.dynamical_map_inverse
        self.orientation = None

        self.points = Points()

        self.points.insert_point(0, self.fixed_point.fixed_point[0])

        self.maximum_spacing = 1e-4

        self.area_cutoff = 1e-1


    def refine_manifold(self, max_passes=20):
        """"
        Adds additional points in areas of the manifold with high curvature

        Checks every consecutive set of three points in the manifold.
        Performs a linear and a parabolic fit between them.
        If the area bounded by those curves is less than self.area_cutoff
        then add additional points.
        Iterates through the manifold until max_passes is reached.

        Parameters:
            max_passes: The max number of passes over the manifold.
            
        Note:
            If the routine terminates from reaching max_passes then there
            are regions above the area cutoff
        """

        pass_counter = 0
        
        while pass_counter < max_passes:
            refine_counter = 0
            old_length = len(self.points)
            
            i = 0
            print("---------")
            print(f"number of points: {old_length}")
            print(f"{pass_counter} passes")
            while i < old_length - 2:
                index = i + refine_counter
                three_points = np.array(self.points.points[index : index + 3])

                x_vals = three_points[:, 0]

                # check if points are so close together to cause numerical instability
                if abs(x_vals[1] - x_vals[0]) < 1e-8 and abs(x_vals[2] - x_vals[1]) < 1e-8:
                    i += 1
                    continue

                try:
                    area = Manifold.curvature_area(three_points)
                except LinAlgError:                         # singular Vandermonde
                    i += 1
                    continue

                if np.abs(area) > self.area_cutoff:
                    self.refine_three_points(three_points, left_index=index)
                    refine_counter += 2
                i += 1

            # If we didn't refine, we are done.
            if refine_counter == 0:
                break
            
            # Otherwise, re-check from the top, possibly with a new length.
            pass_counter += 1
        
        if pass_counter == max_passes:
            print("Warning: refine_manifold reached max_passes without finishing.")


    def refine_three_points(self, points, left_index):
        """
        Method that takes a set of three points and adds new points between
        
        Parameters:
            points: a list of three points to add two new points to
            left_index (int): the index corresponding to the first point
        """

        first_two = points[0:2]
        second_two = points[1:3]

        self.refine_two_points(first_two, left_index)
        self.refine_two_points(second_two, left_index + 2)

        
    def refine_two_points(self, points, left_index):
        """
        Takes two points and adds a new point between them.
        Maps the two points backwards, then does a linear fit between them.
        Interpolates along the linear fit and then maps those points forward.

        Parameters:
            points: list of two points to add a point between
            left_index: the index corresponding to the first point
        """

        right_index = left_index + 1

        first_point = points[0, :]
        second_point = points[1, :]

        self.check_no_inverse()

        # map the two points backwards
        first_point_back = self.dynamical_map_inverse(first_point)
        second_point_back = self.dynamical_map_inverse(second_point)

        # make a point between the pre-iterates
        point_back = (first_point_back + second_point_back) / 2

        point = self.dynamical_map(point_back)

        self.points.insert_point(right_index, point)


    def check_no_inverse(self):
        """
        Raises an error if no inverse map was specified
        """

        if self.dynamical_map_inverse is None:
            raise ValueError("You must include an inverse map to refine manifolds")


    def iterate_manifold(self):
        '''
        Iterates all points in the manifold that have not already been iterated
        '''

        uniterated_indices = [i for i, flag in enumerate(self.points.iterated_flags) if not flag]
        points_to_iterate = [self.points.points[i] for i in uniterated_indices]

        iterated_points = list(map(self.dynamical_map, points_to_iterate))

        self.points.append_points(iterated_points)

        self.refine_manifold()
    
    
    @staticmethod
    def linear_fit(points):
        """
        Takes in three points and gives the linear fit between the first and last
        
        Parameters:
            points: list of three points
        """

        point_one = points[0]
        point_two = points[-1]

        m = (point_two[1] - point_one[1]) / (point_two[0] - point_one[0])
        b = point_one[1] - m * point_one[0]

        return np.poly1d([m, b])


    @staticmethod
    def parabolic_fit(points):
        """
        Takes in three points and gives the parabolic fit between them
        
        Parameters:
            points: list of three points
        """

        x_vals = points[:, 0]
        y_vals = points[:, 1]

        A = np.vstack([x_vals**2, x_vals, np.ones_like(x_vals)]).T
        # Solve for the coefficients [a, b, c].
        coefficient = np.linalg.solve(A, y_vals)

        poly = np.poly1d(coefficient)

        return poly


    @staticmethod
    def curvature_area(points):
            """
            Takes in three points and gives
            the area of the shape made by a parabolic and linear fit
            
            Parameters:
                points: list of three points
            """

            x_min = points[0, 0]
            x_max = points[-1, 0]

            parabola = Manifold.parabolic_fit(points)
            line = Manifold.linear_fit(points)

            difference = lambda x: np.abs(parabola(x) - line(x))

            area, _ = spi.quad(difference, x_min, x_max)

            return area


class Bridge(Manifold):
    '''class for a bridge'''

    def __init__(self):
        super().__init__()


class FundamentalSegment(Manifold):
    '''class for a bridge'''

    def __init__(self):
        super().__init__()

