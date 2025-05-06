from __future__ import annotations
import numpy as np
from .BasePoint import BasePoint

class Point(BasePoint):

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

        super().__init__(x=x, y=y, cdist=cdist, edist=edist)

        self.stretch_param = stretch_param

        self.forward = None
        self.backward = None

        self.next_iterate = None
        self.prev_iterate = None


    def insert_point_forward(self, node: Point):
        """
        Inserts a point (node) object after this point in the linked list

        Paramters:
            node (Point): point to be inserted after
        """
        self.stretch_param = node.stretch_param

        if self.forward is not None:
            self.forward.backward = node
            node.forward = self.forward

        self.forward = node
        node.backward = self


    def insert_point_backward(self, node: Point):
        """
        Inserts a point (node) object before this point in the linked list

        Paramters:
            node (Point): point to be inserted before
        """
        self.stretch_param = node.stretch_param
        
        if self.backward is not None:
            self.backward.forward = node
            node.backward = self.backward

        self.backward = node
        node.forward = self


    def insert_next_iterate(self, node: Point):
        """
        Inserts this object after another point node in the linked list

        Paramters:
            node (Point): point to be inserted after
        """
        self.stretch_param = node.stretch_param

        if self.next_iterate is not None:
            self.next_iterate.prev_iterate = node
            node.next_iterate = self.next_iterate

        self.next_iterate = node
        node.prev_iterate = self


    def insert_prev_iterate(self, node: Point):
        """
        Inserts this object before another point node in the linked list

        Paramters:
            node (Point): point to be inserted before
        """
        self.stretch_param = node.stretch_param

        if self.prev_iterate is not None:
            self.prev_iterate.next_iterate = node
            node.prev_iterate = self.prev_iterate

        self.prev_iterate = node
        node.next_iterate = self

