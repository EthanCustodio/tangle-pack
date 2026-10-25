"""
A manifold point: one row of a :class:`~.PointStore.PointStore`.

A point is a slotted object over one row of a store: its coordinates, cdist and
iterate links live in the store's columns, and its place along a curve is its
row's position on a :class:`~.PointStore.Chain`. The object holds nothing but
its store and row. The store keeps one ``Point`` per row that has one (a
phantom row gets one only when :meth:`PointStore.node` reads it) for the
store's life, so ``id()`` of a point and ``is`` between points stay meaningful.

Dev Notes:

The object itself is the adapter N7 of the numerics arrays plan deletes; until
then every reader that wants nodes gets them from ``store.nodes``.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from .BasePoint import BasePoint
from .PointStore import PointStore, loose_store


class Point(BasePoint):
    """
    A single point of a manifold, read through its store row.

    BasePoint's iterate list (``next_iterate`` / ``prev_iterate``) is the
    store's image and preimage columns; the order along the curve is the
    row's chain.

    Attributes:
        store (PointStore): The store holding the point's row.
        row (int): The point's row in ``store``.
        stretch_param (Optional[float]): Read-only: the stretch parameter of
            the point's chain, None for a row on no chain.
        next_iterate (Point): The image, through the store's image column.
        prev_iterate (Point): The preimage, through the store's preimage column.
    """

    __slots__ = ("store", "row")

    def __init__(
        self,
        x: Optional[float] = None,
        y: Optional[float] = None,
        cdist: Optional[float] = None,
        *,
        store: Optional[PointStore] = None,
    ) -> None:
        """
        Appends a point to a store.

        Args:
            x (Optional[float]): x-coordinate of the point.
            y (Optional[float]): y-coordinate of the point.
            cdist (Optional[float]): Canonical distance from the fixed point.
            store (Optional[PointStore]): The store to append to; the
                process-wide :func:`~.PointStore.loose_store` when omitted.
        """
        store = store if store is not None else loose_store()
        self._attach(store, store.append_one(x, y, cdist))

    @classmethod
    def _over(cls, store: PointStore, row: int) -> Point:
        """The point of a row already appended (in a batch)."""
        point = cls.__new__(cls)
        point._attach(store, row)
        return point

    def _attach(self, store: PointStore, row: int) -> None:
        """Bind the object to its row and register it as the row's point."""
        self.store = store
        self.row = row
        store.nodes[row] = self

    @property
    def stretch_param(self) -> Optional[float]:
        """The stretch parameter of the point's chain (None on no chain)."""
        chain = self.store.chain_of(self.row)
        return None if chain is None else chain.stretch_param

    @property
    def x(self) -> Optional[np.float64]:
        """The x-coordinate (None if the point was built without one)."""
        value = self.store.xy[self.row, 0]
        return None if value != value else value

    @property
    def y(self) -> Optional[np.float64]:
        """The y-coordinate (None if the point was built without one)."""
        value = self.store.xy[self.row, 1]
        return None if value != value else value

    @property
    def cdist(self) -> Optional[float]:
        """The canonical distance, None for a phantom."""
        value = self.store.cdist[self.row]
        return None if value != value else float(value)

    @cdist.setter
    def cdist(self, value: Optional[float]) -> None:
        self.store.set_cdist(self.row, value)

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

        A phantom preiterate counts, exactly as a phantom ``Point`` did before
        the store: the merge's tie rule (:meth:`BasePoint.has_iterate_link`)
        must read the same.

        Returns:
            bool: True if either link column names a row.
        """
        store = self.store
        return bool(store.image[self.row] >= 0 or store.preimage[self.row] >= 0)

    @property
    def next_iterate(self) -> Optional[Point]:
        """The image, or None."""
        row = self.store.image[self.row]
        return None if row < 0 else self.store.node(row)

    @next_iterate.setter
    def next_iterate(self, node: Optional[Point]) -> None:
        self.store._image[self.row] = self._row_of(node)

    @property
    def prev_iterate(self) -> Optional[Point]:
        """The preimage, or None."""
        row = self.store.preimage[self.row]
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
