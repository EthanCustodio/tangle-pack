"""
The point store: every manifold point of a workbench as one row of arrays.

A row holds a point's coordinates, its canonical distance and its two iterate
links (the row of its image and of its preimage). A :class:`~.Point.Point` is a
small object over its row and holds nothing of its own.

The order of points ALONG a curve is a :class:`Chain`: the rows of one
physical curve, in curve order (away from its periodic point). A row lies on
at most one chain, and two store columns say which and where (``_chain``,
``_position``), so both lookups are O(1) and nothing is ever sorted. A curve
(a branch, a bridge, a partial piece, a blast image curve) is a view of a
chain between its root and its tail (:class:`~.BaseManifold.BaseManifold`).

Rows are appended and never removed or reused, so a row id names one point for
the life of its store, and the coordinates of a row never change once written.

"No value" is NaN in every float array and -1 in every row column, and None on
a ``Point`` (until N7 of the numerics arrays plan). No curve row is ever NaN:
the NaN rows are phantoms and the initializer's fictitious tails, all off-curve.

The public columns (``xy``, ``cdist``, ``image``, ``preimage``) are fresh
read-only views per access, so no caller holds a buffer a resize replaces; a
``Point`` reads the private columns.

The engine makes no :class:`~.Point.Point` (N7a): it appends and links rows,
and a row gets its ``Point`` only when something reads it through the adapter
(:meth:`PointStore.node`), as a curve's ``root`` / ``tail`` does. Phantoms are
the coordinate-only preiterates refinement caches for its next layer (cdist
NaN, on no curve, named by one point's preiterate link -- its preimage on an
unstable curve, its image on a stable one -- and naming nothing themselves).

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

A chain never holds a ``BranchPoint`` (D5: its storage is N8's): a chain that
hangs off one names it as its ``anchor``, the left neighbour of position 0,
and the ``BranchPoint``'s slot for that branch holds the chain.
A chain a merge absorbs is retired (no rows, no anchor) and keeps its index.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple, Optional

import numpy as np
from numpy.typing import ArrayLike, NDArray

if TYPE_CHECKING:
    from .BranchPoint import BranchPoint
    from .Intersection import Stability


class End(NamedTuple):
    """A curve end as a VALUE: its store row, coordinates and cdist, no links.

    Compare two ends by ``id``: a row with no cdist reads NaN, and NaN != NaN.
    """

    id: int
    x: float
    y: float
    cdist: float

#: Capacity grows by this factor when an append does not fit.
_GROWTH = 1.5

#: ``_position`` of a row a merge left off its chain (a dropped duplicate).
LEFT_OFF = -2


class Chain:
    """
    The rows of one physical curve, in curve order (away from its periodic point).

    A BranchPoint is never in ``ids`` (D5); when the chain hangs off one, it is the
    chain's ``anchor`` and reads as the left neighbour of position 0.

    Attributes:
        store (PointStore): The store whose rows these are.
        index (int): The chain's index in ``store.chains`` (its ``_chain`` value).
        ids (np.ndarray): The rows, in curve order. Replaced, never written in
            place, by every insert and absorb.
        stretch_param (float): The per-map-step cdist factor of every point on
            the chain (usually ``alpha``).
        anchor (Optional[BranchPoint]): The periodic point the chain hangs off.
    """

    __slots__ = ("store", "index", "ids", "stretch_param", "anchor")

    def __init__(
        self, store: PointStore, index: int, ids: NDArray[np.int64], stretch_param: float
    ) -> None:
        """
        A chain over rows on no other chain (use :meth:`PointStore.new_chain`).

        Args:
            store (PointStore): The store of the rows.
            index (int): Its index in ``store.chains``.
            ids (np.ndarray): The rows in curve order.
            stretch_param (float): The chain's stretch parameter.
        """
        self.store = store
        self.index = index
        self.ids = np.empty(0, dtype=np.int64)
        self.stretch_param = stretch_param
        self.anchor: Optional[BranchPoint] = None
        self.insert(np.zeros(len(ids), dtype=np.int64), ids)

    def __len__(self) -> int:
        """The number of rows on the chain."""
        return len(self.ids)

    def xy_at(self, position: int) -> Optional[NDArray[np.float64]]:
        """
        The coordinates at a position: the anchor's at -1, None off either end.

        Args:
            position (int): A position on the chain.

        Returns:
            Optional[np.ndarray]: A read-only ``(2,)`` view of the row's
            coordinates, the anchor's point for -1, or None.
        """
        if position == -1:
            return None if self.anchor is None else self.anchor.get_point()
        if 0 <= position < len(self.ids):
            return self.store._xy_read[self.ids[position]]
        return None

    def insert(self, positions: ArrayLike, rows: ArrayLike) -> None:
        """
        Inserts rows on no chain before the given positions, in one ``np.insert``.

        Args:
            positions (ArrayLike): Insertion positions in the current ``ids``
                (``np.insert`` semantics: ascending, equal ones keep row order).
            rows (ArrayLike): The rows to insert.
        """
        positions = np.asarray(positions, dtype=np.int64)
        rows = np.asarray(rows, dtype=np.int64)
        if not len(rows):
            return
        assert (self.store._chain[rows] == -1).all(), "a row lies on one chain at most"
        self.ids = np.insert(self.ids, positions, rows)
        self._renumber(int(positions.min()))

    def absorb(self, other: Chain, ids: NDArray[np.int64], start: int = 0) -> None:
        """
        Takes the merged order of this chain and ``other``; ``other`` is retired.

        Every row of the two chains not in ``ids`` leaves its chain (a merge
        left it off, see :data:`LEFT_OFF`); ``other``'s anchor is carried over.

        Args:
            other (Chain): The chain merged into this one.
            ids (np.ndarray): The merged rows, in curve order.
            start (int): ``ids[:start]`` is this chain's unchanged prefix.
        """
        assert other is not self
        assert other.stretch_param == self.stretch_param, (
            f"merging chains of stretch {self.stretch_param!r} and {other.stretch_param!r}"
        )
        assert self.anchor is None or other.anchor is None, "two anchors on one chain"
        store = self.store
        before = np.concatenate((self.ids[start:], other.ids))
        store._chain[before] = -1
        self.ids = np.asarray(ids, dtype=np.int64)
        self.anchor = self.anchor if self.anchor is not None else other.anchor
        other.ids, other.anchor = np.empty(0, dtype=np.int64), None
        self._renumber(start)
        store._position[before[store._chain[before] == -1]] = LEFT_OFF

    def _renumber(self, start: int) -> None:
        """Rewrite the chain and position columns of ``ids[start:]``."""
        tail = self.ids[start:]
        self.store._chain[tail] = self.index
        self.store._position[tail] = np.arange(start, len(self.ids))


class PointStore:
    """
    Append-only columns of points: xy, cdist, image row and preimage row, and
    the chain and position of every row on a curve.

    Attributes:
        nodes (list[Optional[Point]]): Row -> its ``Point``; None for a row no
            reader has asked :meth:`node` for (every row the engine appends).
        chains (list[Chain]): Every chain of the store, by index (retired ones
            included).
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
        self._chain = np.empty(capacity, dtype=np.int32)  # -1 = on no chain
        self._position = np.empty(capacity, dtype=np.int32)  # -1, or LEFT_OFF
        self._n = 0
        self.chains: list[Chain] = []
        self._alias_xy()

    def __len__(self) -> int:
        """The number of rows written."""
        return self._n

    # ── the columns, read-only ─────────────────────────────────────────────

    @property
    def xy(self) -> NDArray[np.float64]:
        """The ``(n, 2)`` coordinates, a fresh read-only view per access."""
        return self._read_only(self._xy)

    @property
    def cdist(self) -> NDArray[np.float64]:
        """The canonical distances (NaN: none), a fresh read-only view."""
        return self._read_only(self._cdist)

    @property
    def image(self) -> NDArray[np.int64]:
        """The row of each row's image (-1: none), a fresh read-only view."""
        return self._read_only(self._image)

    @property
    def preimage(self) -> NDArray[np.int64]:
        """The row of each row's preimage (-1: none), a fresh read-only view."""
        return self._read_only(self._preimage)

    def _read_only(self, column: np.ndarray) -> np.ndarray:
        """The written rows of a column, as a view nobody can write through."""
        view = column[: self._n]
        view.flags.writeable = False
        return view

    # ── writers ────────────────────────────────────────────────────────────

    def append(self, xy: ArrayLike, cdist: Optional[ArrayLike]) -> NDArray[np.int64]:
        """
        Appends one row per coordinate pair, with no iterate links, on no chain.

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
        for column in (self._image, self._preimage, self._chain, self._position):
            column[lo:hi] = -1
        self._n = hi
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
        for column in (self._image, self._preimage, self._chain, self._position):
            column[row] = -1
        self._n = row + 1
        return row

    def link(self, preimages: ArrayLike, images: ArrayLike) -> None:
        """
        Links each preimage row to its image row, both ways.

        Args:
            preimages (ArrayLike): Rows whose image is set.
            images (ArrayLike): Their images, in the same order.
        """
        self._image[preimages] = images
        self._preimage[images] = preimages

    def link_phantoms(self, rows: ArrayLike, phantoms: ArrayLike, stability: Stability) -> None:
        """
        Names each row's phantom preiterate, one way (a phantom names nothing).

        Args:
            rows (ArrayLike): Curve rows.
            phantoms (ArrayLike): Their preiterates: preimages on an unstable
                curve, images on a stable one.
            stability (Stability): The curve's stability.
        """
        (self._preimage if stability == "unstable" else self._image)[rows] = phantoms

    def set_cdist(self, row: int, value: Optional[float]) -> None:
        """
        Writes one row's canonical distance (None is stored as NaN).

        Args:
            row (int): The row.
            value (Optional[float]): Its canonical distance.
        """
        self._cdist[row] = np.nan if value is None else value

    def new_chain(self, rows: ArrayLike, stretch_param: float) -> Chain:
        """
        A chain over rows on no chain yet, in the order given.

        Args:
            rows (ArrayLike): The rows, in curve order.
            stretch_param (float): The chain's stretch parameter.

        Returns:
            Chain: The new chain.
        """
        chain = Chain(self, len(self.chains), np.asarray(rows, dtype=np.int64), stretch_param)
        self.chains.append(chain)
        return chain

    # ── lookups ────────────────────────────────────────────────────────────

    def chain_of(self, row: int) -> Optional[Chain]:
        """The chain a row lies on, or None."""
        index = self._chain[row]
        return None if index < 0 else self.chains[index]

    def position_of(self, row: int) -> int:
        """The row's position on its chain (negative off every chain)."""
        return int(self._position[row])

    def end(self, row: int) -> End:
        """The :class:`End` record of a row (a fresh value per read)."""
        x, y = self._xy[row]
        return End(int(row), float(x), float(y), float(self._cdist[row]))

    @property
    def nbytes(self) -> int:
        """The bytes of the six columns, capacity included."""
        return sum(
            column.nbytes
            for column in (self._xy, self._cdist, self._image, self._preimage, self._chain, self._position)
        )

    def _reserve(self, rows: int) -> None:
        """Grow the columns so that ``rows`` rows fit."""
        capacity = len(self._cdist)
        if rows <= capacity:
            return
        capacity = max(rows, int(capacity * _GROWTH) + 1)
        for name in ("_xy", "_cdist", "_image", "_preimage", "_chain", "_position"):
            old = getattr(self, name)
            new = np.empty((capacity,) + old.shape[1:], dtype=old.dtype)
            new[: self._n] = old[: self._n]
            setattr(self, name, new)
        self._alias_xy()

    def _alias_xy(self) -> None:
        """A read-only alias of the whole coordinate column, rebuilt on every
        resize: ``Point.get_point`` hands out row views of it (a fresh view of
        the column per point read costs five times as much)."""
        self._xy_read = self._xy.view()
        self._xy_read.flags.writeable = False
