from __future__ import annotations
import numpy as np

class Point():

    def __init__(self, x=None, y=None, cdist=None, edist=None, stretch_param=None):
        """
        Basic functionality of a single point in a manifold.
        Contains the methods for two doubly linked lists.

        Manifold linked list tracks the ordering of points 
        walking along the manifold the point is on.

        Iterate linked list tracks the ordering of how points
        map to each other.

        Parameters:
            x (float): x-coordinate of the point
            y (float): y-coordinate of the point
            cdist (float): canonical distance from the fixed point
        """

        self.x = x
        self.y = y

        self.cdist = cdist
        self.edist = edist
        self.stretch_param = stretch_param

        self.next_manifold = None
        
        self.prev_manifold = None

        self.next_iterate = None
        self.prev_iterate = None


    def get_point(self):
        
        return np.array([self.x, self.y])


    def insert_manifold_after(self, node: Point):
        """
        Inserts this object after another point node in the linked list

        Paramters:
            node (Point): point to be inserted after
        """
        self.stretch_param = node.stretch_param

        self.next_manifold = node.next_manifold
        self.prev_manifold = node
        
        node.next_manifold = self


    def insert_manifold_before(self, node: Point):
        """
        Inserts this object before another point node in the linked list

        Paramters:
            node (Point): point to be inserted before
        """
        self.stretch_param = node.stretch_param
        
        self.prev_manifold = node.prev_manifold
        self.next_manifold = node
        
        node.prev_manifold = self


    def insert_iterate_after(self, node: Point):
        """
        Inserts this object after another point node in the linked list

        Paramters:
            node (Point): point to be inserted after
        """
        self.stretch_param = node.stretch_param


        self.next_iterate = node.next_iterate
        self.prev_iterate = node
        
        node.next_iterate = self


    def insert_iterate_before(self, node: Point):
        """
        Inserts this object before another point node in the linked list

        Paramters:
            node (Point): point to be inserted before
        """
        self.stretch_param = node.stretch_param

        self.prev_iterate = node.prev_iterate
        self.next_iterate = node
        
        node.prev_iterate = self