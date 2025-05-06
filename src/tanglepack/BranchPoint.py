import numpy as np
from .Point import Point
from .BasePoint import BasePoint

class BranchPoint(BasePoint):

    def __init__(self, num_branches, x=None, y=None, cdist=None, edist=None):

        super().__init__(x=x, y=y, cdist=cdist, edist=edist)

        self.num_branches = num_branches

        self.cdists = np.zeros(num_branches)
        self.edists = np.zeros(num_branches)

        self.forward_branches = [None] * num_branches
        self.forward_stretch_params = [None] * num_branches

        self.backward_branches = [None] * num_branches
        self.backward_stretch_params = [None] * num_branches

        self.next_iterate = None
        self.prev_iterate = None


    def insert_point_forward(self, node: Point, branch_index):
        """
        Inserts a point (node) after this point node in the linked list at the given branch index

        Paramters:
            node (Point): point to be inserted after
        """

        if self.forward_branches[branch_index] is not None:
            self.forward_branches[branch_index].backward = node
            node.forward = self.forward_branches[branch_index]

        self.forward_branches[branch_index] = node        
        node.backward = self


    def insert_point_backward(self, node: Point, branch_index):
        """
        Inserts this object before another point node in the linked list

        Paramters:
            node (Point): point to be inserted before
        """

        if self.backward_branches[branch_index] is not None:
            self.backward_branches[branch_index].forward = node
            node.backward = self.backward_branches[branch_index]

        self.backward_branches[branch_index] = node
        node.forward = self


    def insert_next_iterate(self, node: Point):
        """
        Inserts this object after another point node in the linked list

        Paramters:
            node (Point): point to be inserted after
        """

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

        if self.prev_iterate is not None:
            self.prev_iterate.next_iterate = node
            node.prev_iterate = self.prev_iterate

        self.prev_iterate = node
        node.next_iterate = self
