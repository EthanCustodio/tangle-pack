from tanglepack.Point import Point
import numpy as np
import scipy.integrate as spi
from .DynamicalSystem import DynamicalSystem
from .BaseManifold import BaseManifold
from .ManifoldView import ManifoldView
from .BranchPoint import BranchPoint
from numpy.linalg import LinAlgError
import time

class ManifoldMachine():


    def __init__(self, system: DynamicalSystem):

        self.system = system
        self.area_cutoff = 1e-2


    def map_manifold(self, manifold: BaseManifold, final_node=None, branch_index=None):

        manifold = ManifoldView(manifold, self.system)

        current_point = manifold.root
        previous_point = None
        count = 0
        while current_point is not None:
            count +=1
            next_point = manifold.walk_fwd(previous_point, current_point, branch_index)
            previous_point, current_point = current_point, next_point

            print(f'The COUNT: {count}')
            time.sleep(0.1)

        current_point = manifold.root
        previous_point = None

        new_manifold = BaseManifold(None, manifold.stability, manifold.stretch_param, manifold.name, None)
        previous_new_manifold = None
        current_new_manifold = None

        while current_point is not None:
            print('ping')
            iterated = current_point.next_iterate
            if iterated is None:

                new_point = manifold.map_fwd(current_point.get_point())
                print(f'new point: {new_point}')
                print(f'current point: {current_point.get_point()}')
                time.sleep(0.1)
                new_dist = current_point.cdist * manifold.stretch_param
                print(f"New distance: {new_dist}")
                iterated = Point(x=new_point[0], y=new_point[1], cdist=new_dist, stretch_param=manifold.stretch_param)
                
                current_point.insert_next_iterate(iterated)
                
            if current_new_manifold is None:
                print('first')
                new_manifold.root = iterated
                current_new_manifold = new_manifold.root

            else:
                print('next')
                if isinstance(current_new_manifold, BranchPoint):
                    if new_manifold.stability == "unstable":
                        current_new_manifold.insert_point_forward(iterated, branch_index)
                    else:  # stable
                        current_new_manifold.insert_point_backward(iterated, branch_index)
                else:
                    # ordinary point → just splice after / before it
                    if new_manifold.stability == "unstable":
                        current_new_manifold.insert_point_forward(iterated)
                    else:
                        current_new_manifold.insert_point_backward(iterated)

                next_point_new = new_manifold.walk_fwd(previous_new_manifold, current_new_manifold, branch_index)
                previous_new_manifold, current_new_manifold = current_new_manifold, next_point_new


                if new_manifold.stability == "unstable":
                    iterated.forward = None
                else:
                    iterated.backward = None


                # previous_new_manifold = current_new_manifold          # always update tail pointer
                # if next_point_new is not None:                        # <-- guard!
                #     current_new_manifold = next_point_new
                        # if next_point_new is not None:
                    # previous_new_manifold, current_new_manifold = current_new_manifold, next_point_new

            next_point = manifold.walk_fwd(previous_point, current_point, branch_index)
            previous_point, current_point = current_point, next_point

            if current_point is final_node:
                break

        if new_manifold.tail is None:
            new_manifold.tail = previous_new_manifold

        return new_manifold
    

    def merge_manifolds(self, manifold_1: BaseManifold, manifold_2: BaseManifold, branch_index=None):

        manifold_1_original = manifold_1

        manifold_1 = ManifoldView(manifold_1, self.system)
        manifold_2 = ManifoldView(manifold_2, self.system)

        print(f'Manifold1: {manifold_1.get_point_array(branch_index=branch_index)}')
        print(f'Manifold1: {manifold_1.get_cdist_array(branch_index=branch_index)}')

        print(f'Manifold2: {manifold_2.get_point_array(branch_index=branch_index)}')
        print(f'Manifold2: {manifold_2.get_cdist_array(branch_index=branch_index)}')


        previous_point_1 = manifold_1.root
        current_point_1 = manifold_1.walk_fwd(None, previous_point_1, branch_index=branch_index)

        previous_point_2 = None
        current_point_2 = manifold_2.root

        while current_point_1 is not None:
            
            left_point = previous_point_1.cdist
            right_point = current_point_1.cdist

            test_point = current_point_2.cdist

            print(left_point, "<–", right_point, "  | inserting", test_point)
            time.sleep(0.1)

            if previous_point_1 is current_point_2:

                print('One')

                next_point_2 = manifold_2.walk_fwd(previous_point_2, current_point_2, branch_index)
                previous_point_2, current_point_2 = current_point_2, next_point_2

            elif left_point < test_point < right_point:

                print('inserted point')

                if isinstance(previous_point_1, BranchPoint):
                    # anchor is the branch leg we’re walking on
                    if manifold_1.stability == "unstable":
                        previous_point_1.insert_point_forward(current_point_2, branch_index)
                    else:  # stable
                        previous_point_1.insert_point_backward(current_point_2, branch_index)
                else:
                    # ordinary point → just splice after / before it
                    if manifold_1.stability == "unstable":
                        previous_point_1.insert_point_forward(current_point_2)
                    else:
                        previous_point_1.insert_point_backward(current_point_2)
                        

                next_point_2 = manifold_2.walk_fwd(previous_point_2, current_point_2, branch_index)
                previous_point_2, current_point_2 = current_point_2, next_point_2

                next_point_1 = manifold_1.walk_fwd(previous_point_1, current_point_1, branch_index)
                previous_point_1, current_point_1 = current_point_1, next_point_1

            else:

                print("THREE")
                next_point_1 = manifold_1.walk_fwd(previous_point_1, current_point_1, branch_index)
                previous_point_1, current_point_1 = current_point_1, next_point_1

            if current_point_1 is None:
                previous_point_1.insert_point_forward(current_point_2)
                manifold_1_original.tail = manifold_2.tail
                break

            if current_point_2 is None:
                break

        return manifold_1_original


    # def grow_manifold(self, manifold: BaseManifold, branch_index=None):

    #     mapped_manifold = self.map_manifold(manifold, manifold.tail, branch_index)

    #     grown_manifold = self.merge_manifolds(manifold, mapped_manifold, branch_index)
    #     # self.merge_manifolds(manifold, mapped_manifold, branch_index)

    #     grown_manifold = self.refine_manifold(grown_manifold, grown_manifold.tail, branch_index=branch_index)
    #     # self.refine_manifold(manifold, manifold.tail, branch_index=branch_index)

    #     return grown_manifold


    def grow_manifold(self, manifold1: BaseManifold, branch_index=None):

        original_manifold = manifold1

        print(f"GROWTH INPUT: {original_manifold.get_point_array(branch_index=branch_index)}")

        manifold = ManifoldView(manifold1, self.system)

        current_point = manifold.root
        previous_point = None

        while current_point is not None:
            stop = False
            # chooses the forward or backward iterate based on stability
            if manifold.stability == "unstable":
                iterate = current_point.next_iterate
            else:
                iterate = current_point.prev_iterate

            if iterate is not None:
                print('iterated already')
                next_point = manifold.walk_fwd(previous_point, current_point, branch_index)
                previous_point, current_point = current_point, next_point

            else:

                new_point = manifold.map_fwd(current_point.get_point())
                print(f'new point: {new_point}')
                print(f'current point: {current_point.get_point()}')
                # time.sleep(0.1)
                new_dist = current_point.cdist * manifold.stretch_param
                print(f"Current Distance: {current_point.cdist}")
                print(f"New distance: {new_dist}")
                iterate = Point(x=new_point[0], y=new_point[1], cdist=new_dist, stretch_param=manifold.stretch_param)

                current_point.insert_next_iterate(iterate)

                next_point = manifold.walk_fwd(previous_point, current_point, branch_index)
                
                
                if next_point == None or next_point is manifold.tail:
                    stop = True

                # temp_current_point = current_point
                # current_point = next_point

                # previous_point, current_point = current_point, next_point

                # insert the new point
                if isinstance(current_point, BranchPoint):
                    if manifold.stability == "unstable":
                        previous_point.next_iterate.insert_point_forward(iterate, branch_index)
                    else:  # stable
                        current_point.prev_iterate.insert_point_backward(iterate, branch_index)
                else:
                    # ordinary point → just splice after / before it
                    if manifold.stability == "unstable":
                        previous_point.next_iterate.insert_point_forward(iterate)
                    else:
                        current_point.prev_iterate.insert_point_backward(iterate)

                next_point = manifold.walk_fwd(previous_point, current_point, branch_index)

                previous_point, current_point = current_point, next_point

                # previous_point = temp_current_point

            # stopping condition
            # if current_point is None:
            if stop == True:
                print(f'stopping growth. Previous_point: {previous_point.get_point()}')
                # original_manifold.tail = previous_point
                # original_manifold.tail = current_point

                temp_current = current_point
                temp_previous = previous_point
                while temp_current is not None:
                    temp = manifold.walk_fwd(temp_previous, temp_current, branch_index)
                    temp_previous, temp_current = temp_current, temp

                original_manifold.tail = temp_previous

                break

        original_manifold = self.refine_manifold(original_manifold, original_manifold.tail, branch_index=branch_index)

        return original_manifold

                
    def iterate_manifold(self, manifold: BaseManifold, branch_index=None):
        """Maps a manifold forward, refines the manifold and returns the mapped manifold"""

        mapped_manifold = self.map_manifold(manifold, branch_index=branch_index)
        self.refine_manifold(manifold=mapped_manifold, final_node=mapped_manifold.tail, branch_index=branch_index)

        return mapped_manifold

    
    #TODO make it so upon iterating it sets the new final node in both iteration methods


    # def iterate_manifold(self, manifold: BaseManifold, final_node=None, branch_index=None):

    #     print("__________________")
    #     print("Iterating Manifold")
    #     print("------------------")

    #     # object to unify the manifold and dynamical system to access mapping
    #     manifold = ManifoldView(manifold, self.system)

    #     current_point = manifold.root
    #     previous_point = None

    #     # we may need to do a special handling of the first point
    #     # when we go to insert it into the ordered manifold list
    #     # there will possibly be no iterate of it to insert off of
    #     # this is primarily for bridges and other segments not attached
    #     # to a fixed point directly
    #     while current_point is not final_node:

    #         # if the current node has been iterated already move on
    #         if current_point.next_iterate != None:

    #             next_point = manifold.walk_fwd(previous_point, current_point, branch_index)
    #             previous_point, current_point = current_point, next_point
    #             continue
            
    #         # otherwise iterate
    #         else:

    #             new_point = manifold.map_fwd(current_point.get_point())
    #             new_dist = current_point.cdist * manifold.stretch_param
    #             new_point = Point(x=new_point[0], y=new_point[1], cdist=new_dist, stretch_param=manifold.stretch_param)
                
    #             current_point.insert_next_iterate(new_point)
            
    #         next_point = manifold.walk_fwd(previous_point, current_point, branch_index)

    #         # --- insert new_pt into the manifold list -----------------
    #         if isinstance(current_point, BranchPoint):
    #             # anchor is the branch leg we’re walking on
    #             if manifold.stability == "unstable":
    #                 current_point.insert_point_forward(new_point, branch_index)
    #             else:  # stable
    #                 current_point.insert_point_backward(new_point, branch_index)
    #         else:
    #             # ordinary point → just splice after / before it
    #             if manifold.stability == "unstable":
    #                 current_point.insert_point_forward(new_point)
    #             else:
    #                 current_point.insert_point_backward(new_point)


    #         previous_point, current_point = current_point, next_point

    #     self.refine_manifold(manifold=manifold, final_node=final_node, branch_index=branch_index)
    #     # self.fill_point_list()

    # TODO do not recheck parts that have aready been checked
    # recursively refine the points
    def refine_manifold(self, manifold1: BaseManifold, final_node=None, max_passes=20, branch_index=None):
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
        original_manifold = manifold1

        print(f'final node: {final_node.get_point()}')

        manifold = ManifoldView(manifold1, self.system)

        pass_counter = 0
        
        while pass_counter < max_passes:

            refine_counter = 0
            i = 0

            print("__________________")
            print(f"Refining! #{pass_counter}")
            print("------------------")

            print(f'Refinement input {original_manifold.get_point_array(branch_index=branch_index)}')

            current_point = manifold.root
            previous_point = None

            # while (current_point or manifold.walk_fwd(previous_point, current_point, branch_index)) != final_node:
            while current_point is not final_node:
                # print(f'Refinement input {original_manifold.get_point_array(branch_index=branch_index)}')
                # time.sleep(0.1)
            # while current_point is not final_node:

                next_point = manifold.walk_fwd(previous_point, current_point, branch_index)

                # print(f'next point: {next_point.get_point()}')

                if next_point is final_node:
                    break

                next_next_point = manifold.walk_fwd(current_point, next_point, branch_index)

                # print(f"next next point: {next_next_point.get_point()}")

                if next_next_point is None:
                    break

                point1 = current_point.get_point()
                point2 = next_point.get_point()
                point3 = next_next_point.get_point()

                three_points = np.array([point1] + [point2] + [point3])
                three_points_object = [current_point, next_point, next_next_point]

                print(f'currently checking: {three_points}')

                x_vals = three_points[:, 0]

                # check if points are so close together to cause numerical instability
                if abs(x_vals[1] - x_vals[0]) < 1e-8 and abs(x_vals[2] - x_vals[1]) < 1e-8:
                    i += 1
                    previous_point, current_point = current_point, next_point
                    continue

                try:
                    area = self.curvature_area(three_points)
                except LinAlgError:                         # singular Vandermonde
                    i += 1
                    previous_point, current_point = current_point, next_point
                    continue

                if np.abs(area) > self.area_cutoff:
                    print("Adding a point")
                    # self.refine_three_points(three_points, current_point)
                    self.refine_three_points(three_points_object, manifold.stability)
                    previous_point, current_point = current_point, next_point
                    refine_counter += 1
                    i += 1
                else:
                    previous_point, current_point = current_point, next_point
                    i += 1

                print(f'loop number: {i}')
                print(f'refinement counter: {refine_counter}')

            # If we didn't refine, we are done.
            if refine_counter == 0:
                print(f'Refinement output {original_manifold.get_point_array(branch_index=branch_index)}')
                print(f'refined cdists: {original_manifold.get_cdist_array(branch_index=branch_index)}')
                break
            
            # Otherwise, re-check from the top, possibly with a new length.
            pass_counter += 1
        
        if pass_counter == max_passes:
            print("Warning: refine_manifold reached max_passes without finishing.")

        return original_manifold


    # def refine_manifold(self, manifold1: BaseManifold, final_node=None, max_passes=20, branch_index=None):
    #     # so the refinement will go til there are no more nodes
    #     # if a final_node is specified, like in a parent Segment class
    #     # then it will only refine that segment
    #     """"
    #     Adds additional points in areas of the manifold with high curvature

    #     Checks every consecutive set of three points in the manifold.
    #     Performs a linear and a parabolic fit between them.
    #     If the area bounded by those curves is less than self.area_cutoff
    #     then add additional points.
    #     Iterates through the manifold until max_passes is reached.

    #     Parameters:
    #         max_passes: The max number of passes over the manifold.
            
    #     Note:
    #         If the routine terminates from reaching max_passes then there
    #         are regions above the area cutoff
    #     """
    #     original_manifold = manifold1

    #     print(f'final node: {final_node.get_point()}')

    #     manifold = ManifoldView(manifold1, self.system)

    #     refine_counter = 0
    #     i = 0

    #     print("__________________")
    #     print(f"Refining! #{pass_counter}")
    #     print("------------------")

    #     print(f'Refinement input {original_manifold.get_point_array(branch_index=branch_index)}')

    #     current_point = manifold.root
    #     previous_point = None

    #     # while (current_point or manifold.walk_fwd(previous_point, current_point, branch_index)) != final_node:
    #     while current_point is not final_node:
    #         # print(f'Refinement input {original_manifold.get_point_array(branch_index=branch_index)}')
    #         # time.sleep(0.1)
    #     # while current_point is not final_node:

    #         next_point = manifold.walk_fwd(previous_point, current_point, branch_index)

    #         # print(f'next point: {next_point.get_point()}')

    #         if next_point is final_node:
    #             break

    #         next_next_point = manifold.walk_fwd(current_point, next_point, branch_index)

    #         # print(f"next next point: {next_next_point.get_point()}")

    #         if next_next_point is None:
    #             break

    #         point1 = current_point.get_point()
    #         point2 = next_point.get_point()
    #         point3 = next_next_point.get_point()

    #         three_points = np.array([point1] + [point2] + [point3])
    #         three_points_object = [current_point, next_point, next_next_point]

    #         print(f'currently checking: {three_points}')

    #         x_vals = three_points[:, 0]

    #         # check if points are so close together to cause numerical instability
    #         if abs(x_vals[1] - x_vals[0]) < 1e-8 and abs(x_vals[2] - x_vals[1]) < 1e-8:
    #             i += 1
    #             previous_point, current_point = current_point, next_point
    #             continue

    #         try:
    #             area = self.curvature_area(three_points)
    #         except LinAlgError:                         # singular Vandermonde
    #             i += 1
    #             previous_point, current_point = current_point, next_point
    #             continue

    #         if np.abs(area) > self.area_cutoff:
    #             print("Adding a point")
    #             # self.refine_three_points(three_points, current_point)
    #             self.refine_three_points(three_points_object, manifold.stability)
    #             previous_point, current_point = current_point, next_point
    #             refine_counter += 1
    #             i += 1
    #         else:
    #             previous_point, current_point = current_point, next_point
    #             i += 1

    #         print(f'loop number: {i}')
    #         print(f'refinement counter: {refine_counter}')

    #     # If we didn't refine, we are done.
    #     if refine_counter == 0:
    #         print(f'Refinement output {original_manifold.get_point_array(branch_index=branch_index)}')
    #         print(f'refined cdists: {original_manifold.get_cdist_array(branch_index=branch_index)}')
    #         break
        
    #     # Otherwise, re-check from the top, possibly with a new length.
    #     pass_counter += 1
    
    #     if pass_counter == max_passes:
    #         print("Warning: refine_manifold reached max_passes without finishing.")

    #     return original_manifold




    # def refine_three_points(self, points, first_point:Point):
    #     """
    #     Method that takes a set of three points and adds new points between
        
    #     Parameters:
    #         points: a list of three points to add two new points to
    #         left_index (int): the index corresponding to the first point
    #     """
        
    #     first_two = points[0:2]
    #     # second_two = points[1:3]

    #     # self.refine_two_points(second_two, first_point)
    #     self.refine_two_points(first_two, first_point)


    def refine_three_points(self, points, stability):
        """
        Method that takes a set of three points and adds new points between
        
        Parameters:
            points: a list of three points to add two new points to
            left_index (int): the index corresponding to the first point
        """
        
        p0 = points[0]
        p1 = points[1]
        p2 = points[2]

        # p0 -> p1 -> p2
        p0_p1_coords = np.array([p0.get_point(), p1.get_point()])
        # self.refine_two_points(p0_p1_coords, stability)
        self.refine_two_points([p0, p1], stability)

        # p0 -> q0 -> p1 -> p2
        p1_p2_coords = np.array([p1.get_point(), p2.get_point()])
        # self.refine_two_points(p1_p2_coords, stability)
        self.refine_two_points([p1, p2], stability)

        # p0 -> q0 -> p1 -> q1 -> p2


    def refine_two_points(self, points: list[Point], stability):
        """
        Takes two points and adds a new point between them.
        Maps the two points backwards, then does a linear fit between them.
        Interpolates along the linear fit and then maps those points forward.

        Parameters:
            points: list of two points to add a point between
            left_index: the index corresponding to the first point
        """

        p0 = points[0]
        p1 = points[1]

        p0_coords = p0.get_point()
        p1_coords = p1.get_point()

        # map the two points backwards
        # TODO Take their preiterates
        if stability == "unstable":
            p0_back = p0.prev_iterate
            p1_back = p1.prev_iterate

            first_point_back = p0_back.get_point()
            second_point_back = p1_back.get_point()

            # if p0_back.cdist is None:
            #     p0.prev_iterate = None
            #     del p0_back

            # if p1_back.cdist is None:
            #     p1.prev_iterate = None
            #     del p1_back

        else:
            p0_back = p0.next_iterate
            p1_back = p1.next_iterate

            first_point_back = p0_back.get_point()
            second_point_back = p1_back.get_point()

            # if p0_back.cdist is None:
            #     p0.next_iterate = None
            #     del p0_back
            # if p1_back.cdist is None:
            #     p1.next_iterate = None
            #     del p1_back


        # make a point between the pre-iterates
        point_back = (first_point_back + second_point_back) / 2

        # compute the cdist of the new point
        alpha = p0.stretch_param
        distance = ((p0.forward.cdist) + (p0.cdist)) / 2

        # map the point back forward
        point = self.system.map(point_back)

        new_point = Point(x=point[0], y=point[1], cdist=distance)

        cached_preiterate = Point(x=point_back[0], y=point_back[1])
        
        if stability == "unstable":
            new_point.insert_prev_iterate(cached_preiterate)
        else:
            new_point.insert_next_iterate(cached_preiterate)

        p0.insert_point_forward(new_point)
 

    # def refine_two_points(self, points, left_point:Point):
    #     """
    #     Takes two points and adds a new point between them.
    #     Maps the two points backwards, then does a linear fit between them.
    #     Interpolates along the linear fit and then maps those points forward.

    #     Parameters:
    #         points: list of two points to add a point between
    #         left_index: the index corresponding to the first point
    #     """

    #     first_point = points[0, :]
    #     second_point = points[1, :]

    #     # map the two points backwards
    #     # TODO Take their preiterates
    #     first_point_back = self.system.map_inv(first_point)
    #     second_point_back = self.system.map_inv(second_point)

    #     # make a point between the pre-iterates
    #     point_back = (first_point_back + second_point_back) / 2

    #     # compute the cdist of the new point
    #     alpha = left_point.stretch_param
    #     distance = ((left_point.forward.cdist) + (left_point.cdist)) / 2

    #     # map the point back forward
    #     point = self.system.map(point_back)

    #     new_point = Point(x=point[0], y=point[1], cdist=distance)

    #     left_point.insert_point_forward(new_point)
 
    
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

            parabola = ManifoldMachine.parabolic_fit(points)
            line = ManifoldMachine.linear_fit(points)

            difference = lambda x: np.abs(parabola(x) - line(x))

            area, _ = spi.quad(difference, x_min, x_max)

            return area




