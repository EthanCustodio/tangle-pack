from __future__ import annotations
import numpy as np
from .BasePoint import BasePoint

"""
Dev Notes:

This class is quite tight. Maybe descibe how the only_forward flag 
work a bit more.
"""


class Point(BasePoint):
    """
    Basic functionality of a single point in a manifold.
    Contains the methods for two doubly linked lists.

    BasePoint implements a linked list to track the iterates.
        next_iterate, prev_iterate
    This class implements a linked list to track the geometric ordering.
        forward, backward

    Attributes:
        stretch_param (float): Amount by which two points on the manifold separate by
            upon a single iteration of the map. Usually denoted 'alpha'.

        forward (Point or BasePoint):
        backward (Point or BasePoint):

        next_iterate (Point or BasePoint):
        prev_iterate (Point or BasePoint):
    """

    def __init__(self, x=None, y=None, cdist=None, edist=None, stretch_param=None):
        """
        Basic functionality of a single point in a manifold.
        Contains the methods for two doubly linked lists.

        'Geometrical' linked list tracks the ordering of points
        walking along the manifold the point is on.

        'Iterate' linked list tracks the ordering of how points
        map to each other.

        Parameters:
            x (float): x-coordinate of the point.
            y (float): y-coordinate of the point.
            cdist (float): Canonical distance from the fixed point.
        """

        super().__init__(x=x, y=y, cdist=cdist, edist=edist)

        self.stretch_param = stretch_param

        self.forward = None
        self.backward = None

    def insert_point_forward(self, node: Point, only_forward: bool = False):
        """
        Inserts a point (node) object after this point in the geometric linked list.

        Paramters:
            node (Point): Point to be inserted.
            only_forward (bool): Flag which switches whether or not we cut off the new
                point from its current manifold.
        """

        if self.forward is node:
            return

        # A freshly created point (e.g. a crossing separator) may not carry a
        # stretch parameter yet; let it inherit ours. Never the other way round.
        if node.stretch_param is None:
            node.stretch_param = self.stretch_param

        if self.forward is not None:
            if not only_forward:
                self.forward.backward = node
                node.forward = self.forward

        self.forward = node
        node.backward = self

    def insert_point_backward(self, node: Point, only_forward: bool = False):
        """
        Inserts a point (node) object before this point in the geometric linked list.

        Paramters:
            node (Point): Point to be inserted.
            only_forward (bool): Flag which switches whether or not we cut off the new
                point from its current manifold.
        """

        if self.backward is node:
            return

        # See insert_point_forward: the inserted node inherits a missing
        # stretch parameter; the host point is never overwritten.
        if node.stretch_param is None:
            node.stretch_param = self.stretch_param

        if self.backward is not None:
            if not only_forward:
                self.backward.forward = node
                node.backward = self.backward

        self.backward = node
        node.forward = self
