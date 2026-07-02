from __future__ import annotations

from typing import Literal, Optional
import numpy as np

"""
Dev Notes:

Consider the data type of _coords more closely. Do I want (2,) or (2,1)?
"""


class BasePoint:
    """
    Lowest level data structure for points on a manifold. Doubly linked list connecting
    points to their next iterate and pre iterate.

    Contains methods for inserting, geting, and checking iterates. ie manipulating the
    linked list.

    Attributes:
        x (float): x-coordinate of the point.
        y (float): y-coordinate of the point.
        _coords (ndarray [2,]): Array coordinates of the point.
        cdist (float): Canonical distance from the fixed point. Defined as {latex here}
        edist (float): Arc-length distance from the fixed point.
        next_iterate (BasePoint): The iterate of the point. Default None.
        prev_iterate (BasePoint): The pre-iterate of the point. Default None.
    """

    def __init__(
        self, x: float = None, y: float = None, cdist: float = None, edist: float = None
    ):
        """
        Initializes a point.

        Args:
            x (float): x-coordinate of the point. Defaults to None.
            y (float): y-coordinate of the point. Defaults to None.
            cdist (float, optional): canonical distance from the fixed point. Defaults to None.
            edist (float, optional): the pre-iterate of the point. Defaults to None.
        """

        self.x = x
        self.y = y

        self.cdist = cdist
        self.edist = edist

        self.next_iterate = None
        self.prev_iterate = None

        self._set_coords()

    def get_point(self) -> np.ndarray:
        """
        Get the coordinates of the point.

        Returns:
            np.ndarray: Array of shape (2, ) containing the coordinates of the point.
        """

        return self._coords

    def get_cdist(
        self, stability: Optional[Literal["unstable", "stable"]] = None
    ) -> float:
        """
        Gets the cdist of the point.

        Args:
            stability (Literal["unstable", "stable"]): Included here because child class
                BranchPoint uses this method to distinguish the cdist on from the
                stable and unstable manifolds. Default to None.

        Returns:
            float: The canonical distance of the point.
        """

        return self.cdist

    def set_x(self, x: float):
        """
        Sets the x-coordinate of the point. Also resets the interal _coords.

        Args:
            x (float): The x-coordinate of the point.
        """

        self.x = x
        self._set_coords()

    def set_y(self, y: float):
        """
        Sets the y-coordinate of the point. Also resets the internal _coords.

        Args:
            y (float): The y-coordinate of the point.
        """

        self.y = y
        self._set_coords()

    def insert_next_iterate(self, node: BasePoint, num_iterates: int = 1):
        """
        Inserts 'node' after this object 'num_iterates' forward
        in the iterate linked list.

        Args:
            node (BasePoint): Iterate of the current point to insert into
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

    def insert_prev_iterate(self, node: BasePoint, num_iterates: int = 1):
        """
        Inserts node before this object 'num_iterates' backwards
        in the iterate linked list.

        Args:
            node (BasePoint): Iterate of the current point to insert into
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

    def get_next_iterate(self, num_iterates: int = 1) -> BasePoint:
        """
        Gets the iterate 'num_iterates' forward.

        Args:
            num_iterates (int): Number of forward iterates.

        Return:
            BasePoint: The iterate 'num_iterates' forward.
        """

        current_node = self

        for _ in range(num_iterates):

            if current_node is None:
                return None

            current_node = current_node.next_iterate

        return current_node

    def get_prev_iterate(self, num_iterates: int = 1) -> BasePoint:
        """
        Gets the pre-iterate 'num_iterates' backwards.

        Arge:
            num_iterates (int): Number of backward iterates.

        Return:
            BasePoint: The iterate 'num_iterates' backward.
        """

        current_node = self

        for _ in range(num_iterates):

            if current_node is None:
                return None

            current_node = current_node.prev_iterate

        return current_node

    def exists_next_iterate(self, num_iterates: int = 1):
        """
        Checks if the iterate 'num_iterates' forward exists.

        Args:
            num_iterates (int): Number of iterates forward to check.
        """

        if self.get_next_iterate(num_iterates) is None:
            return False
        else:
            return True

    def exists_prev_iterate(self, num_iterates: int = 1):
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
    def _set_coords(self):
        """
        Sets the coordinate array based on the x and y value of the point.
        """

        self._coords = np.array([self.x, self.y])
