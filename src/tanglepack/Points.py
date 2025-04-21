import numpy as np

class Points():

    def __init__(self, points=None, edists=None):
        """
        Low level object to store points and routine to add them
        
        Parameters:
            points: list of points to initialize the object with
            edists: list of distances to initialize with
        """
        
        # edists is measured from the periodic point measured ALONG the manifold
        # for iterating one segment you cannot define edists in this same way
        # cdists is ultimately more important than edists 
        self.points = [] if points is None else points

        self.iterated_flags = [False for _ in self.points]

        if edists is None:
            self.edists = [None for _ in range(len(self.points))]
        else:
            self.edists = edists

        self.find_edists()
        #TODO: implement cdists


    # def __str__(self):
    #     """Make print(points) display like a NumPy array"""
    #     return str(self.as_array())
    """implement cannonical distance
        that will fix the iterate manifold problem
        we need to be inserting points into 
        the proper part of the manifold
    """

    def __len__(self):
        return len(self.points)


    def __getitem__(self, index):
        return self.points[index]


    def as_array(self):
        """
        Returns an array representation for numerical operations.
        """

        return np.array([[p[0], p[1]] for p in self.points])


    def find_edists(self):
        """
        Finds the euclidean distance for all points without an edist
        """

        pass


    def insert_point(self, index, point):
        """
        Inserts a point at the given index

        Parameters:
            index (int): index to insert the point
            points (array-like): 2D point to insert

        Note:
            This method allows you to insert a point anywhere in the manifold.
            The user must ensure it is the proper geometrical ordering.    
        """

        point = np.array(point, dtype=np.float64).reshape(-1)

        self.points.insert(index, point)

        # update the iterated flag list to mark the new point as un-iterated
        self.iterated_flags.insert(index, False)

        self.insert_edist(index)
        self.find_edists()


    def append_points(self, points):
        """
        Appends a list of points to the end of the current list of points
        
        Parameters:
            points (array-like): list of points to append
        """
    
        self.points.extend([np.array(pt, dtype=np.float64).reshape(-1) for pt in points])

        # update the iterated flag list to mark the new points as un-iterated
        self.iterated_flags.extend([False] * len(points))

        self.find_edists()


    def insert_edist(self, index, edist=None):
        """
        Inserts an euclidean distance
        """

        self.edists.insert(index, edist)