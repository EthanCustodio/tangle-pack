from typing import Literal

import numpy as np
from .Point import Point
from .BasePoint import BasePoint


class BranchPoint(BasePoint):

    def __init__(
        self, num_branches, cdists: tuple[float, float] = None, x=None, y=None
    ):
        """

        Parameters:
            cdists: index 0 stores unstable canonical distance
                    index 1 stores stable canonical distance
        """

        super().__init__(x=x, y=y)

        self.num_branches = num_branches

        self.cdists = cdists

        self.edists = np.zeros(num_branches)

        self.forward_branches = [None] * num_branches
        self.forward_stretch_params = [None] * num_branches

        self.backward_branches = [None] * num_branches
        self.backward_stretch_params = [None] * num_branches

    def get_cdist(self, stability: Literal["unstable", "stable"]) -> float:
        """
        Get the canonical distance for the manifold of your choice
        """

        index = 0 if stability == "unstable" else 1

        return self.cdists[index]

    def insert_point_forward(
        self, node: Point, branch_index, only_forward: bool = False
    ):
        """
        Inserts a point (node) after this point node in the linked list at the given branch index

        Paramters:
            node (Point): point to be inserted after\
            only_forward (bool): flag which switches whether we cut off node
                from its current manifold
        """
        if self.forward_branches[branch_index] is node:
            pass

        else:
            if self.forward_branches[branch_index] is not None:
                if not only_forward:
                    self.forward_branches[branch_index].backward = node
                    node.forward = self.forward_branches[branch_index]

            self.forward_branches[branch_index] = node
            node.backward = self

    def insert_point_backward(
        self, node: Point, branch_index, only_forward: bool = False
    ):
        """
        Inserts this object before another point node in the linked list

        Paramters:
            node (Point): point to be inserted before
            only_forward (bool): flag which switches whether we cut off node
                from its current manifold
        """
        if self.backward_branches[branch_index] is node:
            pass

        else:
            if self.backward_branches[branch_index] is not None:
                if not only_forward:
                    self.backward_branches[branch_index].forward = node
                    node.backward = self.backward_branches[branch_index]

            self.backward_branches[branch_index] = node
            node.forward = self

    # def insert_next_iterate(self, node: Point):
    #     """
    #     Inserts this object after another point node in the linked list

    #     Paramters:
    #         node (Point): point to be inserted after
    #     """

    #     if self.next_iterate is not None:
    #         raise ValueError("next iterate already exists")

    #     self.next_iterate = node
    #     node.prev_iterate = self

    # def insert_prev_iterate(self, node: Point):
    #     """
    #     Inserts this object before another point node in the linked list

    #     Paramters:
    #         node (Point): point to be inserted before
    #     """

    #     if self.prev_iterate is not None:
    #         raise ValueError("previous iterate already exists")

    #     self.prev_iterate = node
    #     node.next_iterate = self
