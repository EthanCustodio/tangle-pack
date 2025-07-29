from typing import Literal

import numpy as np
from .Point import Point
from .BasePoint import BasePoint

"""
Dev Notes:

Potentially clarify the parameter language in insert_point_backwards()
so it is not only_forward, but only_backward. Or change them both to
be consistent.
"""


class BranchPoint(BasePoint):
    """
    Basic functionality for an intersection point in a manifold.
    Contains methods for two doubly linked lists.

    This class implements a doubly linked list to track geometrical point ordering.
    This class differs from Point in that there are up to two incoming and outgoing
    branches of both a stable and unstable manifold. Both intersection points and
    fixed points are represented by BranchPoints.
        EX: forward_branches, backward_branches
    The BasePoint class implements a linked list to track the iterates of the point.
    There is only one branch of this linked list.
        EX: next_iterate, prev_iterate

    Attributes:
        num_branches (int): _description_
        cdists (List[float]): Canonical distance from the fixed point for both the
            stable and unstable manifolds. Index 0 is unstable, index 1 is stable.
        edists (List[float]): Arclength distance from the fixed point for both the
            stable and unstable manifolds. Index 0 in unstable, index 1 is stable.

        forward_branches (List[Point]): List of the outgoing branches from the
            BranchPoint. In the case of a fixed point the forward branches will
            strictly be the unstable manifolds. In the case of an intersection
            there will be one unstable and one stable.
        forward_stretch_params (List[float]): Geometric stretch associated with
            the outgoing manifolds.

        backward_branches (List[Point]): List of the incoming branches from the
            BranchPoint. In the case of a fixed point the backward branches will
            strictly be the stable manifolds. In the case of an intersection
            there will be one unstable and one stable.
        backward_stretch_params (List[float]): Geometric stretch associated with
            the incoming manifolds.
    """

    def __init__(
        self,
        num_branches: Literal[1, 2],
        cdists: tuple[float, float] = None,
        x=None,
        y=None,
    ):
        """
        Initializes a BranchPoint.

        Parameters:
            x (float): x-coordinate of the point.
            y (float): y-coordinate of the point.
            num_branches (int 1 or 2): Number of branches per manifold attached to the
                branch.
                EX:
                    A fixed point without inversion will have 1 branch. An intersection
                    point will have 2 branches.
            cdists (float, float): index 0 stores unstable canonical distance
                    index 1 stores stable canonical distance.
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
        Get the canonical distance for the manifold of your choice.

        Args:
            stability (Literal["unstable", "stable"]): Stability of the manifold you
                want to access the cdist of. BranchPoints have both a unstable and
                stable cdist.

        Return:
            float: The canonical distance along the given manifold.
        """

        index = 0 if stability == "unstable" else 1

        return self.cdists[index]

    def insert_point_forward(
        self, node: Point, branch_index: int, only_forward: bool = False
    ):
        """
        Inserts a point (node) after this point node in the linked list
        at the given branch index.

        Args:
            node (Point): Point to be inserted into the manifold.
            branch_index (int): The branch of the BranchPoint you wish to insert
                the point into.
            only_forward (bool): Flag which switches whether we cut off node
                from its current manifold. Choosing True will not separate the new
                point from the points that are already forwardly connected to it.
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
        at the given branch index.

        Paramters:
            node (Point): point to be inserted before
            branch_index (int): The branch of the BranchPoint you wish to insert
                the point into.
            only_forward (bool): Flag which switches whether we cut off node
                from its current manifold. Choosing True will not separate the new
                point from the points that are already backwardly connected to it.
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
