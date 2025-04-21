import numpy as np
from .Point import Point

class BranchPoint(Point):

    def __init__(self, num_branches, x=None, y=None, cdist=None, edist=None):

        super().__init__(x=x, y=y, cdist=cdist, edist=edist)

        self.cdists = np.zeros(num_branches)
        self.edists = np.zeros(num_branches)

        self.forward_branches = [None] * num_branches
        self.forward_stretch_params = [None] * num_branches

        self.backward_branches = [None] * num_branches
        self.backward_stretch_params = [None] * num_branches

        self.next_iterate = None
        self.prev_iterate = None


    # these methods might not make much sense since a
    def insert_point_forward(self, node: Point, index):
        """
        Inserts this object after another point node in the linked list

        Paramters:
            node (Point): point to be inserted after
        """

        self.forward_branches[index] = node.forward
        self.backward_branches[index] = node
        
        node.forward[index] = self


    def insert_point_backward(self, node: Point, index):
        """
        Inserts this object before another point node in the linked list

        Paramters:
            node (Point): point to be inserted before
        """
        
        self.backward_branches[index] = node.backward
        self.forward_branches[index] = node
        
        node.backward[index] = self
