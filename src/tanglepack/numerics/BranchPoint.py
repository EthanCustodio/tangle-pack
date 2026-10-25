"""
Branch points: the nodes where several manifold branches meet.

A :class:`BranchPoint` is a fixed point or a crossing, so it can carry up to
two incoming and two outgoing branches of each stability instead of the single
``forward`` / ``backward`` pair a plain :class:`~.Point.Point` has.

Dev Notes:

Potentially clarify the parameter language in insert_point_backwards()
so it is not only_forward, but only_backward. Or change them both to
be consistent.
"""

from __future__ import annotations

from typing import Literal, Optional

import numpy as np

from .Point import Point
from .BasePoint import BasePoint
from .Intersection import Stability


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
        num_branches (int): The eigendirection slots at this point (two: one
            per direction along each eigenvector).
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
        cdists: Optional[tuple[float, float]] = None,
        x: Optional[float] = None,
        y: Optional[float] = None,
    ) -> None:
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

    def get_cdist(self, stability: Stability) -> float:
        """
        Get the canonical distance for the manifold of your choice.

        Args:
            stability (Stability): Stability of the manifold you
                want to access the cdist of. BranchPoints have both a unstable and
                stable cdist.

        Return:
            float: The canonical distance along the given manifold.
        """

        index = 0 if stability == "unstable" else 1

        return self.cdists[index]

    def insert_point_forward(self, node: Point, branch_index: int) -> None:
        """
        Hangs ``node``'s chain off this point at the given branch index: the
        branch's slot names ``node`` and the chain names this point as its anchor.

        Args:
            node (Point): The first point of the branch, on a chain of its store.
            branch_index (int): The branch of the BranchPoint you wish to insert
                the point into.

        Raises:
            ValueError: The slot already holds a point of another chain.
        """
        chain = node.store.chain_of(node.row)
        held = self.forward_branches[branch_index]
        if held is not None and held.store.chain_of(held.row) is not chain:
            raise ValueError(f"branch {branch_index} already holds a point of another chain")
        self.forward_branches[branch_index] = node
        chain.anchor = self

    def insert_point_backward(self, node: Point, branch_index: int) -> None:
        """
        Hangs ``node``'s chain off this point at the given backward branch
        index: the slot names ``node`` and the chain names this point as its
        anchor.

        Args:
            node (Point): The first point of the branch, on a chain of its store.
            branch_index (int): The branch of the BranchPoint you wish to insert
                the point into.

        Raises:
            ValueError: The slot already holds a point of another chain.
        """
        chain = node.store.chain_of(node.row)
        held = self.backward_branches[branch_index]
        if held is not None and held.store.chain_of(held.row) is not chain:
            raise ValueError(f"branch {branch_index} already holds a point of another chain")
        self.backward_branches[branch_index] = node
        chain.anchor = self
