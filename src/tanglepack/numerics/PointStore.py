"""
The point store: every manifold point of a workbench as one row of arrays.

A row holds a point's coordinates, its canonical distance (a float and the
residual its rounding dropped, :mod:`~.precise`) and its two iterate links (the row of its image and of its preimage). A point is its row: no
object stands for it (N7b); a curve's end reads as an :class:`End`, a value
built from the row's columns.

The order of points ALONG a curve is a :class:`Chain`: the rows of one
physical curve, in curve order (away from its periodic point). A row lies on
at most one chain, and two store columns say which and where (``_chain``,
``_position``), so both lookups are O(1) and nothing is ever sorted. A curve
(a branch, a bridge, a partial piece, a blast image curve) is a view of a
chain between its root and its tail (:class:`~.BaseManifold.BaseManifold`).

Rows are appended and never reused, so a row id names one point until the
store is compacted, and the coordinates of a row never change once written.
Only a lean blast removes rows (:meth:`PointStore.compact`, called by
``TangleWorkbench.release_curves``): it renumbers the survivors in order and
hands back the old-to-new map every holder of a row is remapped through.

"No value" is NaN in every float array and -1 in every row column. No curve
row is ever NaN: the NaN rows are phantoms and the initializer's fictitious
tails, all off-curve.

The public columns (``xy``, ``cdist``, ``image``, ``preimage``) are fresh
read-only views per access, so no caller holds a buffer a resize replaces.

The engine appends and links rows (N7a). Phantoms are the coordinate-only
preiterates refinement caches for its next layer (cdist NaN, on no curve,
named by one point's preiterate link -- its preimage on an unstable curve, its
image on a stable one -- and naming nothing themselves).

Dev Notes:

Append-only but for the lean blast (2026-10-29): ``compact`` frees the rows a
released bridge held, and its phantoms. Rows nothing frees stay, the
initializer's detached fictitious tails and the merge-dropped duplicates
among them.

Who owns a store (Q5): the workbench (``TangleWorkbench.point_store``), which
hands it to its initializer; the machine stays stateless and appends to the
store of the curve it grows. A curve built by hand (tests, notebooks) names
its store (``store=``).

A chain never holds a ``BranchPoint`` (D5: it is a view of its ``FixedPoint``,
on no row): a chain that hangs off one names it as its ``anchor``, the left
neighbour of position 0, and the fixed point's slot table holds the chain.
A chain a merge absorbs is retired (no rows, no anchor) and keeps its index.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple, Optional

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .precise import normalise

if TYPE_CHECKING:
    from .BranchPoint import BranchPoint
    from .Intersection import Stability


class End(NamedTuple):
    """A curve end as a VALUE: its store row, coordinates and cdist, no links.

    ``cdist`` is the plain float; ``cdist + cdist_lo`` the precise value.
    Compare two ends by ``id``: a row with no cdist reads NaN, and NaN != NaN.
    """

    id: int
    x: float
    y: float
    cdist: float
    cdist_lo: float = 0.0


#: Capacity grows by this factor when an append does not fit.
_GROWTH = 1.5

#: ``_position`` of a row a merge left off its chain (a dropped duplicate).
LEFT_OFF = -2

#: Every column of the store, in the order a resize or a compaction copies them.
_COLUMNS = ("_xy", "_cdist", "_cdist_lo", "_image", "_preimage", "_chain", "_position")


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
    Append-only columns of points: xy, cdist and its residual, image row and
    preimage row, and the chain and position of every row on a curve.

    Attributes:
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
        self._cdist_lo = np.empty(capacity)  # the residual: 0.0 for a seed row, NaN with no cdist
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

    def append(
        self, xy: ArrayLike, cdist: Optional[ArrayLike], cdist_lo: ArrayLike = 0.0
    ) -> NDArray[np.int64]:
        """
        Appends one row per coordinate pair, with no iterate links, on no chain.

        Args:
            xy (ArrayLike): The ``(m, 2)`` coordinates.
            cdist (Optional[ArrayLike]): Their ``m`` canonical distances, or None
                for rows without one (phantoms).
            cdist_lo (ArrayLike): The residuals of the cdists (:mod:`~.precise`).
                Defaults to 0.0; ignored for rows without a cdist.

        Returns:
            np.ndarray: The ids of the new rows, in order.
        """
        xy = np.asarray(xy, dtype=np.float64).reshape(-1, 2)
        lo, hi = self._n, self._n + len(xy)
        self._reserve(hi)
        self._xy[lo:hi] = xy
        self._cdist[lo:hi] = np.nan if cdist is None else cdist
        self._cdist_lo[lo:hi] = np.nan if cdist is None else cdist_lo
        for column in (self._image, self._preimage, self._chain, self._position):
            column[lo:hi] = -1
        self._n = hi
        return np.arange(lo, hi)

    def append_one(
        self, x: Optional[float], y: Optional[float], cdist: Optional[float], cdist_lo: float = 0.0
    ) -> int:
        """
        Appends one row (a missing coordinate or cdist is stored as NaN).

        Args:
            x (Optional[float]): x-coordinate.
            y (Optional[float]): y-coordinate.
            cdist (Optional[float]): Canonical distance.
            cdist_lo (float): Its residual. Defaults to 0.0.

        Returns:
            int: The id of the new row.
        """
        row = self._n
        self._reserve(row + 1)
        self._xy[row] = (np.nan if x is None else x, np.nan if y is None else y)
        self._cdist[row] = np.nan if cdist is None else cdist
        self._cdist_lo[row] = np.nan if cdist is None else cdist_lo
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

    def set_cdist(self, row: int, value: Optional[float], lo: float = 0.0) -> None:
        """
        Writes one row's canonical distance (None is stored as NaN).

        Args:
            row (int): The row.
            value (Optional[float]): Its canonical distance.
            lo (float): Its residual. Defaults to 0.0.
        """
        self._cdist[row] = np.nan if value is None else value
        self._cdist_lo[row] = np.nan if value is None else lo

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
        return End(int(row), float(x), float(y), float(self._cdist[row]), float(self._cdist_lo[row]))

    def precise_cdist(self, rows: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """
        The precise canonical distances of rows, as normalised pairs.

        The one reader of the residual column: an order or a comparison that
        must hold where plain cdists tie reads this (:func:`.precise.rank`,
        :func:`.precise.compare`). ``hi`` may differ from the plain cdist by
        an ULP where the residual is that large.

        Args:
            rows (ArrayLike): Store rows.

        Returns:
            tuple: ``(hi, lo)``, each of the shape of ``rows``.
        """
        rows = np.asarray(rows, dtype=np.int64)
        return normalise(self._cdist[rows], self._cdist_lo[rows])

    @property
    def nbytes(self) -> int:
        """The bytes of the seven columns, capacity included."""
        return sum(getattr(self, name).nbytes for name in _COLUMNS)

    def compact(self, keep: NDArray[np.bool_]) -> NDArray[np.int64]:
        """
        Drops every row ``keep`` does not hold, renumbering the rest in order.

        Each surviving row keeps its coordinates, cdist (with its residual) and the links into rows
        that survive; a link into a dropped row becomes -1. Every chain keeps
        its surviving rows in their order, renumbered from position 0 (a
        retired chain stays empty). The capacity is not shrunk below the rows
        kept, so the arrays are reallocated once.

        Args:
            keep: One bool per written row.

        Returns:
            np.ndarray: The new row of each old row, -1 for a dropped one. Every
            other holder of a row (a curve end, a crossing's bracket, an edge
            key) is remapped through it by the caller.
        """
        keep = np.asarray(keep, dtype=bool)
        assert len(keep) == self._n, "one keep flag per written row"
        old_to_new = np.full(self._n, -1, dtype=np.int64)
        old_to_new[keep] = np.arange(int(keep.sum()))
        n = int(keep.sum())
        for name in _COLUMNS:
            old = getattr(self, name)
            new = np.empty((max(n, 1),) + old.shape[1:], dtype=old.dtype)
            new[:n] = old[: self._n][keep]
            setattr(self, name, new)
        for column in (self._image, self._preimage):
            linked = column[:n] >= 0
            column[:n][linked] = old_to_new[column[:n][linked]]
        self._n = n
        for chain in self.chains:
            rows = old_to_new[chain.ids]
            chain.ids = rows[rows >= 0]
            chain._renumber(0)
        self._alias_xy()
        return old_to_new

    def _reserve(self, rows: int) -> None:
        """Grow the columns so that ``rows`` rows fit."""
        capacity = len(self._cdist)
        if rows <= capacity:
            return
        capacity = max(rows, int(capacity * _GROWTH) + 1)
        for name in _COLUMNS:
            old = getattr(self, name)
            new = np.empty((capacity,) + old.shape[1:], dtype=old.dtype)
            new[: self._n] = old[: self._n]
            setattr(self, name, new)
        self._alias_xy()

    def _alias_xy(self) -> None:
        """A read-only alias of the whole coordinate column, rebuilt on every
        resize: :meth:`Chain.xy_at` and the refiner hand out row views of it (a
        fresh view of the column per point read costs five times as much)."""
        self._xy_read = self._xy.view()
        self._xy_read.flags.writeable = False
