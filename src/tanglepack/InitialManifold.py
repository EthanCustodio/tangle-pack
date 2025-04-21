import numpy as np
from .Manifold import Manifold

class InitialManifold(Manifold):
    """Class for a manifold attached directly to a fixed point"""

    def __init__(self, fixed_point, stability='unstable'):
        """
        Potentially a class that doesn't need to exist. Should be in the fixed point class
        """

        super().__init__(fixed_point)

        self.stability = stability

        if self.stability == 'unstable':
            index = 0
        else:
            index = 1

        self.direction_from_fixed_point = self.fixed_point.eigenvectors[0][index]
        self.direction_from_fixed_point = self.direction_from_fixed_point.flatten()

        self.get_initial_fundamental_segment()


    def get_first_point(self):
        """
        Computes the first point from the fixed point based on a linear interpolation
        """

        step = self.fixed_point.accuracy

        first_point = self.fixed_point.fixed_point + (step) * self.direction_from_fixed_point

        print(f"eigen {self.direction_from_fixed_point}")

        return np.array(first_point, dtype=np.float64).reshape(-1)
    

    def get_first_point_preiterate(self):
        """
        Get the iterate of the first point
        """

        first_point = self.get_first_point()

        first_preiterate = self.fixed_point.dynamical_map_inverse(first_point)

        return first_preiterate
    

    def get_initial_fundamental_segment(self):
        """
        Computes the initial fundamental segment from iterating the first point
        """

        first_point = self.get_first_point()
        self.points.insert_point(1, first_point)

        first_preiterate = self.get_first_point_preiterate()
        self.points.insert_point(1, first_preiterate)

        self.refine_manifold()


    def map_fundamental_segment(self):
        """
        Maps the fundamental segment forward until it reaches the desired initial length
        """

        mapped_points = [self.dynamical_map(np.array(p, dtype=np.float64).reshape(2)) for p in self.points.points]

        for i in range(len(self.points)):

            self.points.insert_point(-1, mapped_points[i])

        self.refine_manifold()
        
        


        

        
