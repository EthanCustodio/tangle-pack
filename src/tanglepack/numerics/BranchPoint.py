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

from typing import TYPE_CHECKING, Literal, Optional

import numpy as np

from .Intersection import Stability

if TYPE_CHECKING:
    from .PointStore import Chain


class BranchPoint:
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

        forward_branches (List[Optional[Chain]]): The chains of the outgoing
            branches from the BranchPoint. In the case of a fixed point the
            forward branches will strictly be the unstable manifolds. In the
            case of an intersection there will be one unstable and one stable.
        forward_stretch_params (List[float]): Geometric stretch associated with
            the outgoing manifolds.

        backward_branches (List[Optional[Chain]]): The chains of the incoming
            branches from the BranchPoint. In the case of a fixed point the
            backward branches will strictly be the stable manifolds. In the
            case of an intersection there will be one unstable and one stable.
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

        self.x = x
        self.y = y

        self.cdist = None
        self.edist = None

        self.next_iterate = None
        self.prev_iterate = None

        self._set_coords()

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

    def insert_point_forward(self, chain: Chain, branch_index: int) -> None:
        """
        Hangs a chain off this point at the given branch index: the branch's
        slot holds the chain and the chain names this point as its anchor.

        Args:
            chain (Chain or Point): The branch's chain. A ``Point`` stands for
                the chain it lies on (the N7a adapter for hand-built curves,
                which N7b removes).
            branch_index (int): The branch of the BranchPoint you wish to insert
                the chain into.

        Raises:
            ValueError: The slot already holds another chain.
        """
        held = self.forward_branches[branch_index]
        if held is not None and held is not chain:
            raise ValueError(f"branch {branch_index} already holds another chain")
        self.forward_branches[branch_index] = chain
        chain.anchor = self

    def insert_point_backward(self, chain: Chain, branch_index: int) -> None:
        """
        Hangs a chain off this point at the given backward branch index: the
        slot holds the chain and the chain names this point as its anchor.

        Args:
            chain (Chain or Point): The branch's chain. A ``Point`` stands for
                the chain it lies on (the N7a adapter for hand-built curves,
                which N7b removes).
            branch_index (int): The branch of the BranchPoint you wish to insert
                the chain into.

        Raises:
            ValueError: The slot already holds another chain.
        """
        held = self.backward_branches[branch_index]
        if held is not None and held is not chain:
            raise ValueError(f"branch {branch_index} already holds another chain")
        self.backward_branches[branch_index] = chain
        chain.anchor = self

    def get_point(self) -> np.ndarray:
        """
        Get the coordinates of the point.

        Returns:
            np.ndarray: Array of shape (2, ) containing the coordinates of the point.
        """

        return self._coords

    def set_x(self, x: float) -> None:
        """
        Sets the x-coordinate of the point. Also resets the interal _coords.

        Args:
            x (float): The x-coordinate of the point.
        """

        self.x = x
        self._set_coords()

    def set_y(self, y: float) -> None:
        """
        Sets the y-coordinate of the point. Also resets the internal _coords.

        Args:
            y (float): The y-coordinate of the point.
        """

        self.y = y
        self._set_coords()

    def insert_next_iterate(self, node: BranchPoint, num_iterates: int = 1) -> None:
        """
        Inserts 'node' after this object 'num_iterates' forward
        in the iterate linked list.

        Args:
            node (BranchPoint): Iterate of the current point to insert into
                the linked list.
            num_iterations (int): Number of iterates forward to insert the point.

        Raises:
            ValueError: Intermediate iterates does not exist.
            ValueError: The desired iterate does not exist.
        """

        current_node = self

        for _ in range(num_iterates - 1):

            if current_node.next_iterate is None:
                raise ValueError("Intermediate iterates do not exist")

            current_node = current_node.next_iterate

        if current_node.next_iterate is not None:
            raise ValueError(f"Iterate {num_iterates} already exists")

        current_node.next_iterate = node
        node.prev_iterate = current_node

    def insert_prev_iterate(self, node: BranchPoint, num_iterates: int = 1) -> None:
        """
        Inserts node before this object 'num_iterates' backwards
        in the iterate linked list.

        Args:
            node (BranchPoint): Iterate of the current point to insert into
                the linked list.
            num_iterations (int): Number of iterates backward to insert the point.

        Raises:
            ValueError: Intermediate iterates does not exist.
            ValueError: The desired iterate does not exist.
        """

        current_node = self

        for i in range(num_iterates - 1):

            if current_node.prev_iterate is None:
                raise ValueError("Intermediate iterates do not exist")

            current_node = current_node.prev_iterate

        if current_node.prev_iterate is not None:
            raise ValueError(f"Pre-iterate {num_iterates} already exists")

        current_node.prev_iterate = node
        node.next_iterate = current_node

    def get_next_iterate(self, num_iterates: int = 1) -> BranchPoint:
        """
        Gets the iterate 'num_iterates' forward.

        Args:
            num_iterates (int): Number of forward iterates.

        Return:
            BranchPoint: The iterate 'num_iterates' forward.
        """

        current_node = self

        for _ in range(num_iterates):

            if current_node is None:
                return None

            current_node = current_node.next_iterate

        return current_node

    def get_prev_iterate(self, num_iterates: int = 1) -> BranchPoint:
        """
        Gets the pre-iterate 'num_iterates' backwards.

        Arge:
            num_iterates (int): Number of backward iterates.

        Return:
            BranchPoint: The iterate 'num_iterates' backward.
        """

        current_node = self

        for _ in range(num_iterates):

            if current_node is None:
                return None

            current_node = current_node.prev_iterate

        return current_node

    def has_iterate_link(self) -> bool:
        """
        Whether the point has an image or a preimage.

        The merge's tie rule reads it: of two points on one cdist the one with
        an iterate link is the established point and is kept.

        Returns:
            bool: True if either iterate link is set.
        """

        return self.next_iterate is not None or self.prev_iterate is not None

    def exists_next_iterate(self, num_iterates: int = 1) -> bool:
        """
        Checks if the iterate 'num_iterates' forward exists.

        Args:
            num_iterates (int): Number of iterates forward to check.
        """

        if self.get_next_iterate(num_iterates) is None:
            return False
        else:
            return True

    def exists_prev_iterate(self, num_iterates: int = 1) -> bool:
        """
        Checks if the pre-iterate 'num_iterates' backwards exists.

        Args:
            num_iterations (int): Number of iterates backwards to check.
        """

        if self.get_prev_iterate(num_iterates) is None:
            return False
        else:
            return True

    # ---------- internal helpers ----------
    def _set_coords(self) -> None:
        """
        Sets the coordinate array based on the x and y value of the point.
        """

        self._coords = np.array([self.x, self.y])
