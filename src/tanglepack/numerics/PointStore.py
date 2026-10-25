"""
The point store: every manifold point of a workbench as one row of arrays.

A row holds a point's coordinates, its canonical distance and its two iterate
links (the row of its image and of its preimage). A :class:`~.Point.Point` is a
small object over its row: it keeps only the geometric links ``forward`` /
``backward`` and its stretch parameter, and reads everything else here.

Rows are appended and never removed or reused, so a row id names one point for
the life of its store, and the coordinates of a row never change once written.

Phantoms -- the coordinate-only preiterates refinement caches for its next
layer (cdist NaN, on no curve, named by one point's preiterate link -- its
preimage on an unstable curve, its image on a stable one -- and naming nothing
themselves) -- are rows with no :class:`~.Point.Point` at all, until something
reads one through the adapter (:meth:`PointStore.node`).

Dev Notes:

Append-only until N10 (numerics arrays plan, Q7): the only rows the base would
have freed are the initializer's detached fictitious tails, one per seeded
chain; merge-dropped duplicates are off-curve but named by their preimage's
link. Compaction belongs to save/load.

Who owns a store (Q5): the workbench (``TangleWorkbench.point_store``), which
hands it to its initializer; the machine stays stateless and appends to the
store of the curve it grows. A point built without a store goes into the one
process-wide :func:`loose_store`, which is never freed: it is for points built
by hand (tests, notebooks), not for library paths.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

import numpy as np
from numpy.typing import ArrayLike, NDArray

if TYPE_CHECKING:
    from .Point import Point

#: Capacity grows by this factor when an append does not fit.
_GROWTH = 1.5


class PointStore:
    """
    Append-only columns of points: xy, cdist, image row and preimage row.

    Attributes:
        nodes (list[Optional[Point]]): Row -> its ``Point``; None for a phantom
            no reader has asked for.
        xy (np.ndarray): A read-only alias of the coordinate column, rebuilt on
            every resize; ``Point.get_point`` hands out views of it.
    """

    def __init__(self, capacity: int = 1024) -> None:
        """
        Creates an empty store.

        Args:
            capacity (int): Rows allocated up front. Defaults to 1024.
        """
        self._xy = np.empty((capacity, 2))
        self._cdist = np.empty(capacity)  # NaN = no cdist (a phantom)
        self._image = np.empty(capacity, dtype=np.int64)  # -1 = no link
        self._preimage = np.empty(capacity, dtype=np.int64)
        self._n = 0
        self.nodes: list[Optional[Point]] = []
        self._alias_xy()

    def __len__(self) -> int:
        """The number of rows written."""
        return self._n

    def append(self, xy: ArrayLike, cdist: Optional[ArrayLike]) -> NDArray[np.int64]:
        """
        Appends one row per coordinate pair, with no iterate links.

        Args:
            xy (ArrayLike): The ``(m, 2)`` coordinates.
            cdist (Optional[ArrayLike]): Their ``m`` canonical distances, or None
                for rows without one (phantoms).

        Returns:
            np.ndarray: The ids of the new rows, in order.
        """
        xy = np.asarray(xy, dtype=np.float64).reshape(-1, 2)
        lo, hi = self._n, self._n + len(xy)
        self._reserve(hi)
        self._xy[lo:hi] = xy
        self._cdist[lo:hi] = np.nan if cdist is None else cdist
        self._image[lo:hi] = -1
        self._preimage[lo:hi] = -1
        self._n = hi
        self.nodes.extend([None] * (hi - lo))
        return np.arange(lo, hi)

    def append_one(
        self, x: Optional[float], y: Optional[float], cdist: Optional[float]
    ) -> int:
        """
        Appends one row (a missing coordinate or cdist is stored as NaN).

        Args:
            x (Optional[float]): x-coordinate.
            y (Optional[float]): y-coordinate.
            cdist (Optional[float]): Canonical distance.

        Returns:
            int: The id of the new row.
        """
        row = self._n
        self._reserve(row + 1)
        self._xy[row] = (np.nan if x is None else x, np.nan if y is None else y)
        self._cdist[row] = np.nan if cdist is None else cdist
        self._image[row] = -1
        self._preimage[row] = -1
        self._n = row + 1
        self.nodes.append(None)
        return row

    def node(self, row: int) -> Point:
        """
        The ``Point`` of a row, created (once) for a phantom read through the adapter.

        Args:
            row (int): A row id of this store.

        Returns:
            Point: The row's point.
        """
        point = self.nodes[row]
        if point is None:
            from .Point import Point  # Point imports this module

            point = Point._over(self, int(row), None)
        return point

    @property
    def nbytes(self) -> int:
        """The bytes of the four columns, capacity included."""
        return self._xy.nbytes + self._cdist.nbytes + self._image.nbytes + self._preimage.nbytes

    def _reserve(self, rows: int) -> None:
        """Grow the columns so that ``rows`` rows fit."""
        capacity = len(self._cdist)
        if rows <= capacity:
            return
        capacity = max(rows, int(capacity * _GROWTH) + 1)
        for name in ("_xy", "_cdist", "_image", "_preimage"):
            old = getattr(self, name)
            new = np.empty((capacity,) + old.shape[1:], dtype=old.dtype)
            new[: self._n] = old[: self._n]
            setattr(self, name, new)
        self._alias_xy()

    def _alias_xy(self) -> None:
        """Point ``xy`` at the current coordinate column, read-only."""
        self.xy = self._xy.view()
        self.xy.flags.writeable = False


_LOOSE: Optional[PointStore] = None


def loose_store() -> PointStore:
    """The one process-wide store of points built without a store (never freed)."""
    global _LOOSE
    if _LOOSE is None:
        _LOOSE = PointStore()
    return _LOOSE
