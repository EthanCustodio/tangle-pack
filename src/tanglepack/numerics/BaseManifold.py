"""
A manifold curve as a view of a chain: walking, collecting and plotting one curve.

:class:`BaseManifold` stores its ``root`` and ``tail`` as rows of its
``store`` (the root may be the periodic point's
:class:`~.BranchPoint.BranchPoint`) and reads the rows between them off the
:class:`~.PointStore.Chain` they lie on: the chain is found from the root's
row (from the root's slot, which holds the chain, when the root is a
``BranchPoint``) on EVERY read and never cached, so a view stays right through
every insert, append and absorb of its chain.
It exposes the array / plot helpers every other layer reads a curve through.
:class:`~.Bridge.Bridge` is the subclass for a truncated arc.

Dev Notes:

Include return type hints including hints like Union[np.ndarray, list[Point]]

The five array getters are one gather: :meth:`BaseManifold._collect` takes the
row filter (all / iterated / non-iterated), what to read off each row (point
or cdist), and whether to stop at the tail; every public getter is a thin
wrapper over it. A ``BranchPoint`` root is the one node read one at a time
(``_has_iterate`` / ``_iterate_of``): it is no row (D5; its storage is N8's).

Ends are rows inside (N7a); ``root`` / ``tail`` read a row as its store's
``Point``, N3's adapter, made on first read and the same object on every read.
The node adapters (``_walk_nodes``, ``first_node``, ``return_nodes=True``)
have no library caller; N7b deletes them with ``Point``.

A curve reads past its own ends where its chain goes on: ``get_cdist_array``
runs to the chain's end (past a trimmed tail, into blast images growth merged
onto a branch), and so does ``_walk_nodes(stop_at_final=False)``. Only
removing an end from its chain breaks a view, which only the merge's keeper
branch (b) does to the existing point, and no case reaches it.

Versioning: every manifold carries a ``_version`` bumped by
:meth:`bump_version` and by the ``root``/``tail`` setters, and
:attr:`TangleWorkbench.generation` sums those versions, so a retrim or a regrow
invalidates every derived cache; :meth:`Trellis.cached_bridge_polyline` also
keys a bridge's polyline on its version. Points are inserted into the chain,
which the manifold does not observe, so the workbench bumps every bridge's
version on the two operations that can lay new points inside an already-cut
bridge (growth and a re-cut). Code that drives ``ManifoldMachine`` directly,
bypassing the workbench, must call ``bump_version()`` itself. No read is
memoised: every getter gathers its rows afresh (N6, 2026-10-26, removed the
bridge memo: it held about 33 MiB at scale and never served an array a fresh
read would not have).
"""

from __future__ import annotations

from typing import Iterator, Literal, Optional

import numpy as np
from numpy.typing import NDArray
import matplotlib.pyplot as plt

from .FixedPoint import FixedPoint
from .BranchPoint import BranchPoint
from .Point import Point
from .PointStore import Chain, PointStore
from .Intersection import ManifoldKey, Stability


class BaseManifold:
    """
    Implements a structure for keeping track of manifold characteristics: a view
    of the chain of points from ``root`` to ``tail``, with methods for
    extracting information from it and plotting.

    Attributes:
        root (Point or BranchPoint): First point in the manifold, stored as a
            row of ``store`` (or the periodic point) and read as its ``Point``.
        tail (Point or BranchPoint): Final point in the manifold, stored and
            read as ``root``. If not specified the tail is the last point of
            the root's chain.
        store (PointStore): The store of the curve's rows.
        ids (np.ndarray): The rows from the root to the tail.
        stability (Stability): Stability of the manifold.
        stretch_param (float): Amount by which two points on the manifold separate by
            upon a single iteration of the map. Usually denoted 'alpha'.
        fixed_point (FixedPoint): Fixed point that the manifold originates from.
        name (string, optional): Name of the manifold
        branch_index (int): If the manifold is attached to a fixed point with inversion
            this attribute specifies which branch the manifold eminates from.
        manifold_key (ManifoldKey): ``(fixed_point, stability, orbit_index,
            branch_index)`` -- the branch this manifold is. Required at construction.
        _version (int): Monotone counter bumped by :meth:`bump_version` and by
            any reassignment of ``root`` / ``tail``. It is what
            :attr:`TangleWorkbench.generation` reads to notice that manifold
            geometry moved, and what :meth:`Trellis.cached_bridge_polyline` keys
            a bridge's polyline on.
    """

    def __init__(
        self,
        root: Point | BranchPoint | int,
        stability: Stability,
        stretch_param: float,
        fixed_point: FixedPoint,
        name: str = "unnamed",
        tail: Optional[Point | BranchPoint | int] = None,
        branch_index: Optional[int] = None,
        *,
        manifold_key: Optional[ManifoldKey],
        store: Optional[PointStore] = None,
    ) -> None:
        """
        Initializes the manifold.

        Args:
            root (Point, BranchPoint or int): First point in the manifold: the
                periodic point, a row of ``store`` or a ``Point`` (kept as its row).
            stability (Stability): Stability of the manifold.
            stretch_param (float): Amount by which two points on the manifold separate by
            upon a single iteration of the map. Usually denoted 'alpha'.
            fixed_point (FixedPoint): Fixed point that the manifold originates from.
            name (str, optional): Name of the manifold. Defaults to "unnamed".
            tail (Optional[Point, BranchPoint or int], optional): Final point in
                the manifold, as ``root``. If not specified the tail is the last
                point of the root's chain.
            branch_index (Optional[int], optional): If the manifold is attached to a
                fixed point with inversion this attribute specifies which branch the
                manifold eminates from. Defaults to None.
            manifold_key (Optional[ManifoldKey]): Keyword-only and REQUIRED --
                ``(fixed_point, stability, orbit_index, branch_index)``, the single
                source of truth for which branch this manifold is. It is the key
                :attr:`TangleWorkbench.manifolds` stores the manifold under and the
                key every crossing detected on it records. Pass ``None`` only for a
                transient forward image whose branch the caller has yet to advance
                (``ManifoldMachine.iterate_manifold``); the caller must set it before
                the object is indexed or registered.
            store (Optional[PointStore]): The store of the curve's rows.
                Required for a row root; a ``Point`` end brings its own, and a
                periodic-point root takes its slot chain's.
        """

        self._version: int = 0
        self._store = store
        self._root: int | BranchPoint | None = None
        self._tail: int | BranchPoint | None = None
        self.stability = stability
        self.stretch_param = stretch_param
        self.fixed_point = fixed_point
        self.name = name
        self.branch_index = branch_index
        self.manifold_key = manifold_key
        self.root = root
        self.tail = tail
        if self._tail is None:
            self._find_tail()

    # ── versioning ──────────────────────────────────────────────────────────

    @property
    def version(self) -> int:
        """The manifold's mutation counter (see the module Dev Notes)."""
        return self._version

    def bump_version(self) -> None:
        """
        Record that this manifold's geometry changed.

        Called by the ``root`` / ``tail`` setters and by
        :class:`TangleWorkbench` around anything that can insert points into an
        existing curve (growth, a re-cut). Code that mutates a manifold's chain
        directly must call it too.
        """
        self._version += 1

    @property
    def store(self) -> Optional[PointStore]:
        """The store of the curve's rows: the one given, else (for a
        periodic-point root) its slot chain's, adopted on first read; None for
        a periodic point whose slot is empty."""
        if self._store is None and isinstance(self._root, BranchPoint):
            chain = self._slot_chain()
            if chain is not None:
                self._store = chain.store
        return self._store

    @property
    def root(self) -> Optional[Point | BranchPoint]:
        """First point in the manifold: the periodic point, or the root row as its
        store's ``Point`` (N3's adapter, made on first read and the same object on
        every read; N7b makes ends records). Reassigning it bumps the version."""
        return self._end(self._root)

    @root.setter
    def root(self, end: Optional[Point | BranchPoint | int]) -> None:
        self._root = self._stored(end)
        self.bump_version()

    @property
    def tail(self) -> Optional[Point | BranchPoint]:
        """Final point in the manifold, read as :attr:`root` is. Reassigning it
        bumps the version."""
        return self._end(self._tail)

    @tail.setter
    def tail(self, end: Optional[Point | BranchPoint | int]) -> None:
        # Deliberately unconditional -- do NOT short-circuit on `end ==
        # self._tail`. `_find_tail()` is what every growth pass ends with, and a
        # pass that only REFINED the curve re-assigns the same tail row; that
        # assignment is the only signal the manifold gets that its points moved.
        self._tail = self._stored(end)
        self.bump_version()

    @property
    def ids(self) -> NDArray[np.int64]:
        """The rows from the root to the tail, in curve order, as a fresh array
        (a periodic-point root has no row and is not in it)."""
        return np.array(self._rows()[2], dtype=np.int64)

    def _stored(self, end: Optional[Point | BranchPoint | int]) -> int | BranchPoint | None:
        """An end as the curve stores it: a ``Point`` becomes its row (its store
        adopted), an int stays a row of :attr:`store`; a BranchPoint or None is
        kept as given."""
        if end is None or isinstance(end, BranchPoint):
            return end
        if isinstance(end, Point):
            row, store = end.row, end.store
            if self._store is None:
                self._store = store
        else:
            row, store = int(end), self.store
        assert store is not None and store is self._store, "an end is a row of the curve's store"
        return row

    def _end(self, end: int | BranchPoint | None) -> Optional[Point | BranchPoint]:
        """A stored end as the API reads it (the row's ``Point``, the adapter)."""
        return end if end is None or isinstance(end, BranchPoint) else self._store.node(end)

    def _find_tail(self) -> None:
        """Set the tail to the last point of the root's chain."""
        _head, chain, _lo, _hi = self._span(stop_at_final=False)
        self.tail = self._root if chain is None else int(chain.ids[-1])

    def _span(
        self, final_node: Optional[Point | int] = None, stop_at_final: bool = True
    ) -> tuple[Optional[BranchPoint], Optional[Chain], int, int]:
        """
        Where the curve lies: its BranchPoint root (if any), its chain and the
        slice ``[lo, hi)`` of the chain from the root to ``final_node``.

        A BranchPoint root starts at position 0 of the chain in its slot.

        Args:
            final_node (Optional[Point | int]): Last node (inclusive), a row or
                its ``Point``. None (the default) reads to the tail, or to the
                chain's end when the tail is None; a node not on the chain at or
                after the root reads to the chain's end, as a walk that never
                meets it did.
            stop_at_final (bool): If False, ignore final_node and read to the
                chain's end.

        Returns:
            The BranchPoint root or None, the chain (None when the curve has no
            chain: a lone point, or a BranchPoint whose branch was never
            initialized) and the slice bounds.

        Raises:
            ValueError: The root is a BranchPoint and :attr:`branch_index` is None.
        """
        root = self._root
        head = root if isinstance(root, BranchPoint) else None
        if head is not None:
            chain, lo = self._slot_chain(), 0
        elif root is None:
            return None, None, 0, 0
        else:
            chain, lo = self._store.chain_of(root), self._store.position_of(root)
        if chain is None:
            return head, None, 0, 0
        hi = len(chain)
        if stop_at_final:
            final = self._tail if final_node is None else final_node
            if isinstance(final, Point):
                final = final.row
            if head is not None and final is head:
                hi = lo
            elif isinstance(final, (int, np.integer)) and chain.store.chain_of(final) is chain:
                position = chain.store.position_of(final)
                hi = position + 1 if position >= lo else hi
        return head, chain, lo, hi

    def _rows(
        self, final_node: Optional[Point | int] = None, stop_at_final: bool = True
    ) -> tuple[Optional[BranchPoint], Optional[PointStore], NDArray[np.int64]]:
        """The BranchPoint root (if any), the store and the rows of the curve after
        the root, in curve order (a lone point off every chain is its own row)."""
        head, chain, lo, hi = self._span(final_node, stop_at_final)
        if chain is not None:
            return head, chain.store, chain.ids[lo:hi]
        root = self._root
        if root is None or head is not None:
            return head, None, np.empty(0, dtype=np.int64)
        return None, self._store, np.array([root], dtype=np.int64)

    def _walk_nodes(
        self, final_node: Optional[Point] = None, stop_at_final: bool = True
    ) -> Iterator[Point]:
        """
        Yield every node from the root along the stability direction.

        An adapter N7b deletes (no library path calls it): each row is read
        as its store's ``Point``; read :attr:`ids` instead.

        Args:
            final_node (Optional[Point]): Last node to yield (inclusive).
                Defaults to the manifold's tail (None, the default, reads to
                the tail, or to the chain's end when the tail is None).
            stop_at_final (bool): If False, ignore final_node and read to the
                end of the root's chain.
        """
        head, store, rows = self._rows(final_node, stop_at_final)
        nodes = [] if store is None else [store.node(row) for row in rows.tolist()]
        return iter(nodes if head is None else [head, *nodes])

    def first_node(
        self, branch_index: Optional[int] = None
    ) -> Optional[Point | BranchPoint]:
        """
        Return the first ordinary point of this manifold, stepping past the root
        :class:`BranchPoint` when the manifold is anchored to one.

        A manifold rooted at a fixed point (or at a crossing) starts on a
        :class:`BranchPoint`, which carries up to two outgoing branches per
        stability; the point the curve actually begins at is position 0 of the
        chain its slot holds, on the branch this manifold belongs to. An
        adapter N7b deletes (no library path calls it): the engine reads the
        slot's chain (:meth:`_slot_chain`).

        Args:
            branch_index (Optional[int]): Which branch of the root BranchPoint to
                leave on. Defaults to :attr:`branch_index`.

        Returns:
            Optional[Point | BranchPoint]: The root itself when it is an ordinary
            point, otherwise the first point on the requested branch, as its
            store's ``Point``. ``None`` if the manifold has no root, or the
            branch was never initialized.

        Raises:
            ValueError: The root is a BranchPoint and neither ``branch_index``
                nor :attr:`branch_index` says which branch to leave on.
        """
        root = self.root
        if not isinstance(root, BranchPoint):
            return root
        chain = self._slot_chain(branch_index)
        return None if chain is None else chain.store.node(chain.ids[0])

    def _slot_chain(self, branch_index: Optional[int] = None) -> Optional[Chain]:
        """
        The chain hanging in the root BranchPoint's slot for this curve's
        stability and branch (None for an empty slot).

        Args:
            branch_index (Optional[int]): The slot's branch. Defaults to
                :attr:`branch_index`.

        Raises:
            ValueError: Neither ``branch_index`` nor :attr:`branch_index` says
                which branch to leave on.
        """
        if branch_index is None:
            branch_index = self.branch_index
        if branch_index is None:
            raise ValueError("Must supply branch_index when starting walk from root BranchPoint")
        root = self._root
        slots = root.forward_branches if self.stability == "unstable" else root.backward_branches
        return slots[branch_index]

    def _has_iterate(self, node: Point, num_iterates: int) -> bool:
        """Whether ``node`` has its ``num_iterates``-step image along this stability."""
        if self.stability == "unstable":
            return node.exists_next_iterate(num_iterates)
        return node.exists_prev_iterate(num_iterates)

    def _iterate_of(self, node: Point, num_iterates: int) -> Point:
        """The ``num_iterates``-step image of ``node`` along this stability."""
        if self.stability == "unstable":
            return node.get_next_iterate(num_iterates)
        return node.get_prev_iterate(num_iterates)

    def _collect(
        self,
        kind: Literal["all", "iterated", "non_iterated"] = "all",
        *,
        value: Literal["point", "cdist"] = "point",
        num_iterates: int = 1,
        final_node: Optional[Point] = None,
        return_nodes: bool = False,
        stop_at_final: bool = True,
    ) -> list[Point] | NDArray[np.float64]:
        """
        The single gather behind every array getter on this class.

        Args:
            kind (Literal["all", "iterated", "non_iterated"]): Which nodes to
                keep. ``"iterated"`` keeps only nodes that already have a
                ``num_iterates``-step image (and yields that IMAGE, not the node);
                ``"non_iterated"`` keeps only the nodes that do not.
            value (Literal["point", "cdist"]): What to read off each kept node
                when ``return_nodes`` is False. Ignored when it is True.
            num_iterates (int): Step count for the ``kind`` filter.
            final_node (Optional[Point]): Last node to read to (inclusive);
                defaults to the manifold's tail.
            return_nodes (bool): Return the nodes themselves instead of an
                array (the adapter N7b deletes: each row read as its store's
                ``Point``; no library path asks for it).
            stop_at_final (bool): False reads to the end of the root's chain,
                ignoring ``final_node`` and the tail.

        Returns:
            list[Point] or np.ndarray: The nodes, or an ``(N, 2)`` coordinate
            array / ``(N, 1)`` cdist array. An empty read gives ``[]`` or an
            empty array.

        Raises:
            ValueError: The BranchPoint root claims a ``num_iterates`` iterate
                that is None.
        """
        head, store, rows = self._select(kind, num_iterates, final_node, stop_at_final)
        if return_nodes:
            nodes = [] if store is None else [store.node(row) for row in rows.tolist()]
            return nodes if head is None else [head, *nodes]
        if head is None and not len(rows):
            return np.array([])
        if value == "cdist":
            cdists = store.cdist[rows] if len(rows) else np.empty(0)
            if head is not None:
                cdists = np.concatenate((np.array([head.cdist]), cdists))
            return cdists.reshape(-1, 1)
        points = store.xy[rows] if len(rows) else np.empty((0, 2))
        if head is not None:
            points = np.vstack((head.get_point(), points))
        return points

    def _select(
        self,
        kind: Literal["all", "iterated", "non_iterated"],
        num_iterates: int = 1,
        final_node: Optional[Point | int] = None,
        stop_at_final: bool = True,
    ) -> tuple[Optional[BranchPoint], Optional[PointStore], NDArray[np.int64]]:
        """
        The head, store and rows :meth:`_collect` gathers, in curve order.

        For ``"iterated"`` they are the ``num_iterates``-step images of the rows
        that have one, and the BranchPoint head's image; for ``"non_iterated"``
        the rows with none, and the head when it has none.

        Raises:
            ValueError: The BranchPoint root claims a ``num_iterates`` iterate
                that is None.
        """
        head, store, rows = self._rows(final_node, stop_at_final)
        if kind == "all":
            return head, store, rows
        if head is not None and self._has_iterate(head, num_iterates) != (kind == "iterated"):
            head = None
        if len(rows):
            iterates = self._iterate_rows(store, rows, num_iterates)
            rows = iterates[iterates >= 0] if kind == "iterated" else rows[iterates < 0]
        if head is not None and kind == "iterated":
            head = self._iterate_of(head, num_iterates)
            if head is None:
                raise ValueError("Iterate computed incorrectly, NoneType added")
        return head, store, rows

    def _iterate_rows(
        self, store: PointStore, rows: NDArray[np.int64], num_iterates: int
    ) -> NDArray[np.int64]:
        """The row of each row's ``num_iterates``-step image along this
        stability (-1 where there is none)."""
        link = store.image if self.stability == "unstable" else store.preimage
        for _ in range(num_iterates):
            rows = np.where(rows >= 0, link[np.maximum(rows, 0)], -1)
        return rows

    def get_point_array(
        self, final_node: Optional[Point] = None, return_nodes: bool = False
    ) -> list[Point] | NDArray[np.float64]:
        """
        Reads the manifold in curve order, root first, and returns either
        a list of Point objects or an array of (x, y) coordinates.

        Parameters:
            final_node (Optional[Point]): If specified, stop reading at this node.
            return_nodes (bool): If True, return a list of Point objects;
                else return an np.ndarray of [x, y]. Defaults to False.

        Returns:
            list[Point] or np.ndarray of shape (N, 2)

        Note:
            Every call gathers afresh: the array (or node list) is the caller's
            own, safe to modify.
        """
        return self._collect("all", final_node=final_node, return_nodes=return_nodes)

    @property
    def points(self) -> NDArray[np.float64]:
        """The ``(N, 2)`` coordinates from the root to the TAIL (the periodic point
        first when the root is one), a fresh array per read."""
        return self._collect("all")

    @property
    def cdists(self) -> NDArray[np.float64]:
        """The ``(N,)`` canonical distances from the root to the TAIL, in step with
        :attr:`points` (unlike :meth:`get_cdist_array`, which reads to the chain's end)."""
        return self._collect("all", value="cdist").reshape(-1)

    def get_cdist_array(
        self, final_node: Optional[Point] = None, return_nodes: bool = False
    ) -> list[Point] | NDArray[np.float64]:
        """
        Reads the manifold in curve order, root first, and returns either
        a list of Point objects or an array of their canonical distances.

        Note:
            This read has always covered the whole curve to the end of its
            chain; final_node is accepted for signature symmetry with the other
            getters but is not honored (a trimmed tail does not shorten the
            returned array).

        Returns:
            list[Point] or np.ndarray of shape (N, 1)
        """
        return self._collect(
            "all",
            value="cdist",
            final_node=final_node,
            return_nodes=return_nodes,
            stop_at_final=False,
        )

    def get_non_iterated_point_array(
        self,
        num_iterates: int = 1,
        final_node: Optional[Point] = None,
        return_nodes: bool = False,
    ) -> list[Point] | NDArray[np.float64]:
        """
        Returns the points that do not have a num_iterates-step iterate yet.

        Args:
            num_iterates (int, optional): Number of iterates. Defaults to 1.
            final_node (Optional[Point]): If specified, stop reading at this node.
            return_nodes (bool, optional): If True, return a list of Point objects;
                else return an np.ndarray of [x, y]. Defaults to False.

        Returns:
            list[Point] or np.ndarray of shape (N, 2)
        """
        return self._collect(
            "non_iterated",
            num_iterates=num_iterates,
            final_node=final_node,
            return_nodes=return_nodes,
        )

    def get_non_iterated_cdist_array(
        self, num_iterates: int = 1, final_node: Optional[Point] = None
    ) -> NDArray[np.float64]:
        """
        Returns the canonical distances of the points that do not have a
        num_iterates-step iterate yet.

        Args:
            num_iterates (int, optional): Number of iterates. Defaults to 1.
            final_node (Optional[Point]): If specified, stop reading at this node.

        Returns:
            np.ndarray of shape (N, 1)
        """
        return self._collect(
            "non_iterated",
            value="cdist",
            num_iterates=num_iterates,
            final_node=final_node,
        )

    def get_iterated_point_array(
        self,
        num_iterates: int = 1,
        final_node: Optional[Point] = None,
        return_nodes: bool = False,
    ) -> list[Point] | NDArray[np.float64]:
        """
        Reads the manifold in curve order, root first, and returns either
        a list of Point objects or an array of (x, y) coordinates corresponding to
        the num_iterates-step iterates of the points that have them.

        Parameters:
            final_node (Optional[Point]): If specified, stop reading at this node.
            return_nodes (bool): If True, return list of Point objects;
                else return np.ndarray of [x, y]. Defaults to False.

        Returns:
            list[Point] or np.ndarray of shape (N, 2)

        Raises:
            ValueError: A point claims to have an iterate that is actually None.
        """
        return self._collect(
            "iterated",
            num_iterates=num_iterates,
            final_node=final_node,
            return_nodes=return_nodes,
        )

    def plot(
        self,
        color: str = "blue",
        branch_index: Optional[int] = None,
        show_points: bool = False,
        **kwargs,
    ) -> None:
        """
        Plots the manifold points.

        Parameters:
            color (str): Color of the manifold line.
            branch_index (int, optional): Branch index if starting from a branch point.
            show_points (bool): Whether to show individual points.
            **kwargs: Additional kwargs for plt.plot().
        """
        points = self.get_point_array()

        if points.size == 0:
            raise ValueError("No points available to plot!")

        plt.plot(points[:, 0], points[:, 1], color=color, **kwargs)

        if show_points:
            plt.scatter(points[:, 0], points[:, 1], color=color, s=10, alpha=0.6)

        plt.title(f"Manifold Plot ({self.stability.capitalize()})")
        plt.axis("equal")
