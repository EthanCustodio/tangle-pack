import numpy as np
from scipy.optimize import newton as newton_method
from functools import partial
from scipy.differentiate import jacobian as jacob
from .Point import Point


class MultiPoint(Point):

    def __init__(self, num_manifolds, x=None, y=None, cdist=None, edist=None):

        super().__init__(x=x, y=y, cdist=cdist, edist=edist)

        self.cdists = np.zeros(num_manifolds)
        self.edists = np.zeros(num_manifolds)

        self.next_manifolds = [None for i in range(num_manifolds)]
        self.unstable_stretch_params = [None for i in range(num_manifolds)]

        self.prev_manifolds = [None for i in range(num_manifolds)]
        self.stable_stretch_params = [None for i in range(num_manifolds)]

        self.next_iterate = None
        self.prev_iterate = None


    # these methods might not make much sense since a
    def insert_manifold_after(self, node: Point, index):
        """
        Inserts this object after another point node in the linked list

        Paramters:
            node (Point): point to be inserted after
        """

        self.next_manifolds[index] = node.next_manifold
        self.prev_manifolds[index] = node
        
        node.next_manifold[index] = self


    def insert_manifold_before(self, node: Point, index):
        """
        Inserts this object before another point node in the linked list

        Paramters:
            node (Point): point to be inserted before
        """
        
        self.prev_manifolds[index] = node.prev_manifold
        self.next_manifolds[index] = node
        
        node.prev_manifold[index] = self
