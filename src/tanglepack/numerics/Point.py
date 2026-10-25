"""
A manifold point: geometric ordering on top of the iterate list.

:class:`Point` adds the second doubly linked list -- ``forward`` /
``backward``, the ordering of points *along the curve* -- to the iterate list
it inherits from :class:`~.BasePoint.BasePoint`.

A point is a slotted object over one row of a
:class:`~.PointStore.PointStore`: its coordinates, cdist and iterate links live
in the store's columns, and only the geometric links and the stretch parameter
live on the object. The store keeps one ``Point`` per row for the store's life,
so ``id()`` of a point and ``is`` between points stay meaningful.

Dev Notes:

This class is quite tight. Maybe descibe how the only_forward flag 
work a bit more.

``forward`` / ``backward`` are still object links: they are the adapter the
order arrays of N4 replace (numerics arrays plan).
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from .BasePoint import BasePoint
from .PointStore import PointStore, loose_store


class Point(BasePoint):
    """
    Basic functionality of a single point in a manifold.
    Contains the methods for two doubly linked lists.

    BasePoint implements a linked list to track the iterates.
        next_iterate, prev_iterate
    This class implements a linked list to track the geometric ordering.
        forward, backward

    Attributes:
        store (PointStore): The store holding the point's row.
        row (int): The point's row in ``store``.
        stretch_param (float): Amount by which two points on the manifold separate by
            upon a single iteration of the map. Usually denoted 'alpha'.

        forward (Point or BasePoint):
        backward (Point or BasePoint):

        next_iterate (Point): The image, through the store's image column.
        prev_iterate (Point): The preimage, through the store's preimage column.
    """

    __slots__ = ("store", "row", "forward", "backward", "stretch_param")

    def __init__(
        self,
        x: Optional[float] = None,
        y: Optional[float] = None,
        cdist: Optional[float] = None,
        stretch_param: Optional[float] = None,
        *,
        store: Optional[PointStore] = None,
    ) -> None:
        """
        Appends a point to a store.

        Args:
            x (Optional[float]): x-coordinate of the point.
            y (Optional[float]): y-coordinate of the point.
            cdist (Optional[float]): Canonical distance from the fixed point.
            stretch_param (Optional[float]): The manifold's stretch parameter.
            store (Optional[PointStore]): The store to append to; the
                process-wide :func:`~.PointStore.loose_store` when omitted.
        """
        store = store if store is not None else loose_store()
        self._attach(store, store.append_one(x, y, cdist), stretch_param)

    @classmethod
    def _over(cls, store: PointStore, row: int, stretch_param: Optional[float]) -> Point:
        """The point of a row already appended (in a batch)."""
        point = cls.__new__(cls)
        point._attach(store, row, stretch_param)
        return point

    def _attach(self, store: PointStore, row: int, stretch_param: Optional[float]) -> None:
        """Bind the object to its row and register it as the row's point."""
        self.store = store
        self.row = row
        self.stretch_param = stretch_param
        self.forward = None
        self.backward = None
        store.nodes[row] = self

    @property
    def x(self) -> Optional[np.float64]:
        """The x-coordinate (None if the point was built without one)."""
        value = self.store._xy[self.row, 0]
        return None if value != value else value

    @property
    def y(self) -> Optional[np.float64]:
        """The y-coordinate (None if the point was built without one)."""
        value = self.store._xy[self.row, 1]
        return None if value != value else value

    @property
    def cdist(self) -> Optional[float]:
        """The canonical distance, None for a phantom."""
        value = self.store._cdist[self.row]
        return None if value != value else float(value)

    @cdist.setter
    def cdist(self, value: Optional[float]) -> None:
        self.store._cdist[self.row] = np.nan if value is None else value

    def get_point(self) -> np.ndarray:
        """
        Get the coordinates of the point.

        Returns:
            np.ndarray: A read-only ``(2,)`` view of the point's row. A row's
            coordinates never change, so the view stays right after the store
            grows.
        """
        return self.store.xy[self.row]

    def has_iterate_link(self) -> bool:
        """
        Whether the point has an image or a preimage row.

        A phantom preimage counts, exactly as a phantom ``Point`` did before
        the store: the merge's tie rule (:meth:`BasePoint.has_iterate_link`)
        must read the same.

        Returns:
            bool: True if either link column names a row.
        """
        store = self.store
        return bool(store._image[self.row] >= 0 or store._preimage[self.row] >= 0)

    @property
    def next_iterate(self) -> Optional[Point]:
        """The image, or None."""
        row = self.store._image[self.row]
        return None if row < 0 else self.store.node(row)

    @next_iterate.setter
    def next_iterate(self, node: Optional[Point]) -> None:
        self.store._image[self.row] = self._row_of(node)

    @property
    def prev_iterate(self) -> Optional[Point]:
        """The preimage, or None."""
        row = self.store._preimage[self.row]
        return None if row < 0 else self.store.node(row)

    @prev_iterate.setter
    def prev_iterate(self, node: Optional[Point]) -> None:
        self.store._preimage[self.row] = self._row_of(node)

    def _row_of(self, node: Optional[Point]) -> int:
        """The row an iterate link to ``node`` stores (-1 for None).

        Raises:
            TypeError: ``node`` is not a ``Point`` of this point's store.
        """
        if node is None:
            return -1
        if not isinstance(node, Point) or node.store is not self.store:
            raise TypeError(
                f"an iterate link must name a Point of the same store, got {type(node).__name__}"
            )
        return node.row

    def insert_point_forward(self, node: Point, only_forward: bool = False) -> None:
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

    def insert_point_backward(self, node: Point, only_forward: bool = False) -> None:
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
