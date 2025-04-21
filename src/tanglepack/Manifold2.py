import numpy as np
import scipy.integrate as spi
from .Point import Point
from .FixedPoint2 import FixedPoint2
from .MultiPoint import MultiPoint
from numpy.linalg import LinAlgError


class Manifold2():


    def __init__(self, fixed_point:FixedPoint2, root_node:Point, stability:str, branch_index=0):

        self.fixed_point = fixed_point

        self.root_node = root_node

        self.maximum_spacing = 1e-4

        self.area_cutoff = 1e-1

        self.branch_index = branch_index
        self.is_multi_point = isinstance(self.root_node, MultiPoint)
        self.check_branch_index()

        self.stability = stability
        self.strech_param = None
        self.stability_housekeeping()

        self.point_list = None
        self.fill_point_list()


    def fill_point_list(self):
        """
        Iterates through all points in the manifold and stores their x, y values in a list
        """

        num_points = 0
        current_node = self.root_node

        while current_node.next_manifold != None:
            
            num_points += 1
            current_node = self.walk_manifold(current_node)


        self.point_list = np.array([np.zeros((2, 1)) for i in range(num_points)])
        current_node = self.root_node
        index = 0

        while current_node.next_manifold != None:

            self.point_list[index, :, 0] = current_node.get_point()

            index += 1
            current_node = self.walk_manifold(current_node)


    def stability_housekeeping(self):

        if self.stability == 'unstable':

            if self.is_multi_point:
                self.root_node.next_manifold = self.root_node.next_manifolds[self.branch_index]
                self.strech_param = self.root_node.unstable_stretch_params[self.branch_index]

        elif self.stability == 'stable':

            if self.is_multi_point:
                self.root_node.prev_manifold = self.root_node.prev_manifolds[self.branch_index]
                self.strech_param = self.root_node.stable_stretch_params[self.branch_index]

        else:
            raise ValueError("Must specify the stability of the manifold")


    def walk_manifold(self, point:Point):
        """
        Gives the next point along the manifold in the stability direction

        Parameters:
            point (Point): point to walk from
        """

        if self.stability == 'unstable':
            return point.next_manifold
        
        elif self.stability == 'stable':
            return point.prev_manifold
 

    def walk_manifold_inverse(self, point:Point):
        """
        Gives the previous point along the manifold in the stability direction

        Parameters:
            point (Point): point to walk from
        """

        if self.stability == 'unstable':
            return point.prev_manifold
        
        elif self.stability == 'stable':
            return point.next_manifold


    def map_manifold(self, point):

        if self.stability == 'unstable':
            return self.fixed_point.dynamical_map(point)
        
        elif self.stability == 'stable':
            return self.fixed_point.dynamical_map_inverse(point)


    def map_manifold_inverse(self, point):

        if self.stability == 'unstable':
            return self.fixed_point.dynamical_map_inverse(point)
        
        elif self.stability == 'stable':
            return self.fixed_point.dynamical_map(point)
        

    def check_branch_index(self):
        """
        Checks if the root node is a multipoint if a branch_index was given
        """

        if self.is_multi_point and self.branch_index is None:
            raise ValueError("Must specify branch_index when root_node has multiple branches")


    def refine_manifold(self, final_node=None, max_passes=20):
        # so the refinement will go til there are no more nodes
        # if a final_node is specified, like in a parent Segment class
        # then it will only refine that segment
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
            
            i = 0

            print("__________________")
            print(f"Refining! #{pass_counter}")
            print("------------------")

            if self.branch_index is None:
                current_point = self.root_node
            else:
                # I think the only time a MultiPoint will be the root is at a fixed point
                # at an intersection we'll have another point on the other side as the root
                current_point = self.root_node.next_manifolds[self.branch_index]

            while (current_point and current_point.next_manifold) != final_node:

                point1 = current_point.get_point()
                point2 = current_point.next_manifold.get_point()
                point3 = current_point.next_manifold.next_manifold.get_point()

                three_points = np.array([point1] + [point2] + [point3])

                x_vals = three_points[:, 0]

                # check if points are so close together to cause numerical instability
                if abs(x_vals[1] - x_vals[0]) < 1e-8 and abs(x_vals[2] - x_vals[1]) < 1e-8:
                    i += 1
                    current_point = current_point.next_manifold
                    continue

                try:
                    area = Manifold2.curvature_area(three_points)
                except LinAlgError:                         # singular Vandermonde
                    i += 1
                    current_point = current_point.next_manifold
                    continue

                if np.abs(area) > self.area_cutoff:
                    self.refine_three_points(three_points, current_point)
                    current_point = current_point.next_manifold.next_manifold
                    refine_counter += 1
                i += 1

            # If we didn't refine, we are done.
            if refine_counter == 0:
                break
            
            # Otherwise, re-check from the top, possibly with a new length.
            pass_counter += 1
        
        if pass_counter == max_passes:
            print("Warning: refine_manifold reached max_passes without finishing.")


    def iterate_manifold(self, final_node=None):

        print("__________________")
        print("Iterating Manifold")
        print("------------------")

        current_node = self.root_node

        # we may need to do a special handling of the first point
        # when we go to insert it into the ordered manifold list
        # there will possibly be no iterate of it to insert off of
        # this is primarily for bridges and other segments not attached
        # to a fixed point directly
        while current_node != final_node:

            if current_node.next_iterate != None:
                current_node = self.walk_manifold(current_node)
                continue
                
            else:
                new_point = self.map_manifold(current_node.get_point())
                new_dist = current_node.cdist * self.strech_param
                new_point = Point(x=new_point[0], y=new_point[1], cdist=new_dist, stretch_param=self.strech_param)

                new_point.insert_iterate_after(current_node)
                # this does not work if the first point has not already been iterated
                # if it has then it would be skipped above
                new_point.insert_manifold_after(current_node.prev_manifold.next_iterate)

            current_node = self.walk_manifold(current_node)

        self.refine_manifold()
        self.fill_point_list()


    def plot(self):

        pass

    
    def refine_three_points(self, points, first_point:Point):
        """
        Method that takes a set of three points and adds new points between
        
        Parameters:
            points: a list of three points to add two new points to
            left_index (int): the index corresponding to the first point
        """

        first_two = points[0:2]
        second_two = points[1:3]

        # self.refine_two_points(second_two, first_point)
        self.refine_two_points(first_two, first_point)

        
    def refine_two_points(self, points, left_point:Point):
        """
        Takes two points and adds a new point between them.
        Maps the two points backwards, then does a linear fit between them.
        Interpolates along the linear fit and then maps those points forward.

        Parameters:
            points: list of two points to add a point between
            left_index: the index corresponding to the first point
        """

        first_point = points[0, :]
        second_point = points[1, :]

        self.check_no_inverse()

        # map the two points backwards
        first_point_back = self.map_manifold_inverse(first_point)
        second_point_back = self.map_manifold_inverse(second_point)

        # make a point between the pre-iterates
        point_back = (first_point_back + second_point_back) / 2

        # compute the cdist of the new point
        alpha = self.strech_param
        distance = alpha * (((left_point.next_manifold.cdist / alpha) + (left_point.cdist / alpha)) / 2)

        # map the point back forward
        point = self.map_manifold(point_back)

        new_point = Point(x=point[0], y=point[1], cdist=distance)

        new_point.insert_manifold_after(left_point)


    def check_no_inverse(self):
        """
        Raises an error if no inverse map was specified
        """

        if self.fixed_point.dynamical_map_inverse is None:
            raise ValueError("You must include an inverse map to refine manifolds")
 
    
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

            parabola = Manifold2.parabolic_fit(points)
            line = Manifold2.linear_fit(points)

            difference = lambda x: np.abs(parabola(x) - line(x))

            area, _ = spi.quad(difference, x_min, x_max)

            return area

