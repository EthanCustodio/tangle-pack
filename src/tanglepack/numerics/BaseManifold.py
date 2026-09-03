"""
The linked-list manifold: walking, collecting and plotting one curve.

:class:`BaseManifold` owns a ``root`` and a ``tail`` into the geometric list
of :class:`~.Point.Point` objects that make up one branch of one manifold, and
exposes the walk / array / plot helpers every other layer reads a curve
through. :class:`~.Bridge.Bridge` is the subclass for a truncated arc.

Dev Notes:

Include return type hints including hints like Union[np.ndarray, list[Point]]

The five array getters are one walk: :meth:`BaseManifold._collect` takes the
node filter (all / iterated / non-iterated), what to read off each node (point
or cdist), and whether to stop at the tail; every public getter is a thin
wrapper over it. The old string-dispatched ``_iter_method`` /
``getattr(node, "exists_next_iterate")`` trick is gone -- ``_has_iterate`` and
``_iterate_of`` branch on the stability directly.

Versioning and the point-array memo (Phase 4): every manifold carries a
``_version`` bumped by :meth:`bump_version` and by the ``root``/``tail``
setters, and :attr:`TangleWorkbench.generation` sums those versions, so a
retrim or a regrow invalidates every derived cache. The point-array memo built
on that version is enabled ONLY on :class:`Bridge` (``_memoise_walks``): points
are inserted through ``Point.insert_point_forward/backward``, which the manifold
cannot observe, so a growable manifold has no airtight invalidation and a stale
array would be worse than no cache. A bridge is a fixed arc between two
crossings and the workbench bumps every bridge's version on the two operations
that can lay new points inside such an arc (growth and a re-cut). Code that
drives ``ManifoldMachine`` directly, bypassing the workbench, must call
``bump_version()`` itself.
"""

from __future__ import annotations

from typing import Iterator, Literal, Optional

import numpy as np
from numpy.typing import NDArray
import matplotlib.pyplot as plt

from .FixedPoint import FixedPoint
from .BranchPoint import BranchPoint
from .Point import Point
from .Intersection import ManifoldKey, Stability


class BaseManifold:
    """
    Implements a structure for keeping track of manifold characteristics. Contains
    methods for walking along manifolds, extracting information from manfiolds, and
    plotting.

    Attributes:
        root (Point or BranchPoint): First point in the manifold.
        tail (Point or BranchPoint): Final point in the manifold. If not specified the
            tail will be set my walking from the root until None is reached.
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
            geometry moved, and what invalidates the point-array memo.
    """

    #: Whether :meth:`get_point_array` is served from a memo. False here (a
    #: growable manifold gains points through ``Point`` without the manifold
    #: seeing it), True on :class:`Bridge`. See the module Dev Notes.
    _memoise_walks: bool = False

    def __init__(
        self,
        root: Point | BranchPoint,
        stability: Stability,
        stretch_param: float,
        fixed_point: FixedPoint,
        name: str = "unnamed",
        tail: Optional[Point | BranchPoint] = None,
        branch_index: Optional[int] = None,
        *,
        manifold_key: Optional[ManifoldKey],
    ) -> None:
        """
        Initializes the manifold.

        Args:
            root (Point or BranchPoint): First point in the manifold.
            stability (Stability): Stability of the manifold.
            stretch_param (float): Amount by which two points on the manifold separate by
            upon a single iteration of the map. Usually denoted 'alpha'.
            fixed_point (FixedPoint): Fixed point that the manifold originates from.
            name (str, optional): Name of the manifold. Defaults to "unnamed".
            tail (Optional[Point or BranchPoint], optional): Final point in the manifold.
                If not specified the tail will be set my walking from the root until None is reached.
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
        """

        self._version: int = 0
        self._walk_cache: dict = {}
        self.root = root
        self.tail = tail
        self.stability = stability
        self.stretch_param = stretch_param
        self.fixed_point = fixed_point
        self.name = name
        self.branch_index = branch_index
        self.manifold_key = manifold_key
        if self.tail is None:
            self._find_tail()

    # ── versioning ──────────────────────────────────────────────────────────

    @property
    def version(self) -> int:
        """The manifold's mutation counter (see the module Dev Notes)."""
        return self._version

    def bump_version(self) -> None:
        """
        Record that this manifold's geometry changed and drop its memoised walks.

        Called by the ``root`` / ``tail`` setters and by
        :class:`TangleWorkbench` around anything that can insert points into an
        existing curve (growth, a re-cut). Code that mutates a manifold's linked
        list directly must call it too.
        """
        self._version += 1
        self._walk_cache.clear()

    @property
    def root(self) -> Optional[Point | BranchPoint]:
        """First point in the manifold. Reassigning it bumps the version."""
        return self._root

    @root.setter
    def root(self, point: Optional[Point | BranchPoint]) -> None:
        self._root = point
        self.bump_version()

    @property
    def tail(self) -> Optional[Point | BranchPoint]:
        """Final point in the manifold. Reassigning it bumps the version."""
        return self._tail

    @tail.setter
    def tail(self, point: Optional[Point | BranchPoint]) -> None:
        # Deliberately unconditional -- do NOT short-circuit on `point is
        # self._tail`. `_find_tail()` is what every growth pass ends with, and a
        # pass that only REFINED the curve re-assigns the same tail object; that
        # assignment is the only signal the manifold gets that its points moved.
        self._tail = point
        self.bump_version()

    def _find_tail(self) -> None:
        """Walks until None is reached and set the tail"""

        previous_point = None
        current_point = self.root

        while current_point is not None:
            next_point = self.walk_fwd(previous_point, current_point)
            previous_point, current_point = current_point, next_point

        self.tail = previous_point

    def _walk_nodes(
        self, final_node: Optional[Point] = None, stop_at_final: bool = True
    ) -> Iterator[Point]:
        """
        Yield every node from the root along the stability direction.

        The single traversal all the array getters share.

        Args:
            final_node (Optional[Point]): Last node to yield (inclusive).
                Defaults to the manifold's tail; None walks the physical end
                of the linked list.
            stop_at_final (bool): If False, ignore final_node and walk the
                whole linked list.
        """
        if final_node is None:
            final_node = self.tail

        prev = None
        current = self.root

        while current is not None:
            yield current

            if stop_at_final and current is final_node:
                break

            prev, current = current, self.walk_fwd(
                prev, current, branch_index=self.branch_index
            )

    @staticmethod
    def _stack_points(
        points: list, return_nodes: bool
    ) -> list[Point] | NDArray[np.float64]:
        """Return the collected nodes as-is, or stacked into an array."""
        if not points:
            return [] if return_nodes else np.array([])
        return points if return_nodes else np.vstack(points)

    def first_node(
        self, branch_index: Optional[int] = None
    ) -> Optional[Point | BranchPoint]:
        """
        Return the first ordinary point of this manifold, stepping past the root
        :class:`BranchPoint` when the manifold is anchored to one.

        A manifold rooted at a fixed point (or at a crossing) starts on a
        :class:`BranchPoint`, which carries up to two outgoing branches per
        stability; the point the curve actually begins at is one ``walk_fwd``
        step past it, on the branch this manifold belongs to. This is the single
        implementation of that step -- ``ManifoldMachine._branch_view``,
        ``ManifoldInitializer.construct_kevin_way`` and the growth drivers all go
        through it.

        Args:
            branch_index (Optional[int]): Which branch of the root BranchPoint to
                leave on. Defaults to :attr:`branch_index`.

        Returns:
            Optional[Point | BranchPoint]: The root itself when it is an ordinary
            point, otherwise the first point on the requested branch. ``None`` if
            the manifold has no root, or the branch was never initialized.

        Raises:
            ValueError: The root is a BranchPoint and neither ``branch_index``
                nor :attr:`branch_index` says which branch to leave on.
        """
        root = self.root
        if not isinstance(root, BranchPoint):
            return root
        return self.walk_fwd(None, root, branch_index)

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
        The single walk behind every array getter on this class.

        Args:
            kind (Literal["all", "iterated", "non_iterated"]): Which nodes to
                keep. ``"iterated"`` keeps only nodes that already have a
                ``num_iterates``-step image (and yields that IMAGE, not the node);
                ``"non_iterated"`` keeps only the nodes that do not.
            value (Literal["point", "cdist"]): What to read off each kept node
                when ``return_nodes`` is False. Ignored when it is True.
            num_iterates (int): Step count for the ``kind`` filter.
            final_node (Optional[Point]): Last node to walk to (inclusive);
                defaults to the manifold's tail.
            return_nodes (bool): Return the nodes themselves instead of an array.
            stop_at_final (bool): False walks the whole linked list, ignoring
                ``final_node`` and the tail.

        Returns:
            list[Point] or np.ndarray: The nodes, or an ``(N, 2)`` coordinate
            array / ``(N, 1)`` cdist array. An empty walk gives ``[]`` or an
            empty array.

        Raises:
            ValueError: A node claims a ``num_iterates`` iterate that is None.
        """
        collected = []

        for node in self._walk_nodes(final_node, stop_at_final=stop_at_final):
            if kind != "all":
                if self._has_iterate(node, num_iterates) != (kind == "iterated"):
                    continue

                if kind == "iterated":
                    node = self._iterate_of(node, num_iterates)
                    if node is None:
                        raise ValueError("Iterate computed incorrectly, NoneType added")

            if return_nodes:
                collected.append(node)
            elif value == "cdist":
                collected.append(node.cdist)
            else:
                collected.append(node.get_point())

        return self._stack_points(collected, return_nodes)

    def get_point_array(
        self, final_node: Optional[Point] = None, return_nodes: bool = False
    ) -> list[Point] | NDArray[np.float64]:
        """
        Walks along the manifold in the stability direction and returns either
        a list of Point objects or an array of (x, y) coordinates.

        Parameters:
            final_node (Optional[Point]): If specified, stop walking at this node.
            return_nodes (bool): If True, return a list of Point objects;
                else return an np.ndarray of [x, y]. Defaults to False.

        Returns:
            list[Point] or np.ndarray of shape (N, 2)

        Note:
            On a class that memoises walks (:class:`Bridge`) the array is
            cached against the manifold's version and handed back READ-ONLY --
            it is shared with every other caller, so it must not be written to;
            take a copy to modify it. The node list is copied per call, so it is
            safe to mutate. See the module Dev Notes for the invalidation rules.
        """
        if not self._memoise_walks:
            return self._collect(
                "all", final_node=final_node, return_nodes=return_nodes
            )

        # The node itself is the key (Points hash by identity), so the memo
        # holds it alive and no recycled id can alias two different walks.
        cache_key = (final_node, bool(return_nodes))
        cached = self._walk_cache.get(cache_key)
        if cached is None:
            cached = self._collect(
                "all", final_node=final_node, return_nodes=return_nodes
            )
            if isinstance(cached, np.ndarray):
                cached.flags.writeable = False
            self._walk_cache[cache_key] = cached
        return list(cached) if return_nodes else cached

    def get_cdist_array(
        self, final_node: Optional[Point] = None, return_nodes: bool = False
    ) -> list[Point] | NDArray[np.float64]:
        """
        Walks along the manifold in the stability direction and returns either
        a list of Point objects or an array of their canonical distances.

        Note:
            This walk has always covered the whole linked list; final_node is
            accepted for signature symmetry with the other getters but is not
            honored (a trimmed tail does not shorten the returned array).

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
            final_node (Optional[Point]): If specified, stop walking at this node.
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
            final_node (Optional[Point]): If specified, stop walking at this node.

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
        Walks along the manifold in the stability direction and returns either
        a list of Point objects or an array of (x, y) coordinates corresponding to
        the num_iterates-step iterates of the points that have them.

        Parameters:
            final_node (Optional[Point]): If specified, stop walking at this node.
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

    def walk_fwd(
        self, prev: Optional[Point], node: Point, branch_index: Optional[int] = None
    ) -> Optional[Point]:
        """
        Return the next point along the manifold walking away from the fixed point.
        If `node` is a BranchPoint, we exit on the other branch of the
        same stability type we entered on.
        `prev` is the point we just came from (None at the root).

        Args:
            branch_index : int, optional
            Which branch of a BranchPoint to follow.  If omitted (None),
            `self.branch_index` is used.

        Returns:
            Point or BranchPoint
        """
        if branch_index is None:
            branch_index = self.branch_index

        if isinstance(node, BranchPoint):
            return self._branch_forward(prev, node, branch_index)

        # Ordinary point: follow whichever pointer is "forward"
        return node.forward if self.stability == "unstable" else node.backward

    def walk_back(
        self, nxt: Optional[Point], node: Point, branch_index: Optional[int] = None
    ) -> Optional[Point]:
        """
        The inverse of `walk_fwd`: step one link *backward* toward the fixed point.
        `nxt` is the point we are coming from.

        Args:
            branch_index : int, optional
            Which branch of a BranchPoint to follow.  If omitted (None),
            `self.branch_index` is used.

        Returns:
            Point or BranchPoint
        """
        if branch_index is None:
            branch_index = self.branch_index

        if isinstance(node, BranchPoint):
            return self._branch_backward(nxt, node, branch_index)

        return node.backward if self.stability == "unstable" else node.forward

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

    # ---------- internal helpers ----------
    def _branch_forward(
        self, prev: Point, bp: BranchPoint, branch_index: Optional[int] = None
    ) -> Point:
        """
        Choose the correct outgoing branch at a BranchPoint when moving
        'forward' along the manifold.

        Args:
            prev (Point): Point we walked into the BranchPoint from
            bp (BranchPoint): BranchPoint we are walking through
            branch_index (Optional[int], optional): If there is inversion we may
                potentially need this, currently unsused I think. Defaults to None.

        Raises:
            ValueError: Must supply branch_index when starting walk from
                root BranchPoint
            ValueError: Prev node is not connected to this BranchPoint

        Returns:
            Point
        """

        if self.stability == "stable":
            branches_out = bp.backward_branches
            branches_in = bp.forward_branches
        else:
            branches_out = bp.forward_branches
            branches_in = bp.backward_branches

        if prev is None:
            if branch_index is None:
                raise ValueError(
                    "Must supply branch_index when starting walk from root BranchPoint"
                )
            return branches_out[branch_index]

        for i, point in enumerate(branches_in):
            if point is prev:
                return branches_out[i]  # toggle branch
        raise ValueError("Prev node is not connected to this BranchPoint")

    def _branch_backward(
        self, nxt: Point, bp: BranchPoint, branch_index: Optional[int] = None
    ) -> Point:
        """
        Choose the correct outgoing branch when walking *backward*.
        Symmetric to _branch_forward.

        Args:
            nxt (Point): Point we walked into the BranchPoint from
            bp (BranchPoint): BranchPoint we are walking through
            branch_index (Optional[int], optional): If there is inversion we may
                potentially need this, currently unsused I think. Defaults to None.

        Raises:
            ValueError: Must supply branch_index when starting walk from
                root BranchPoint
            ValueError: Prev node is not connected to this BranchPoint

        Returns:
            Point
        """

        if self.stability == "stable":
            branches_out = bp.backward_branches
            branches_in = bp.forward_branches
        else:
            branches_out = bp.forward_branches
            branches_in = bp.backward_branches

        if nxt is None:
            if branch_index is None:
                raise ValueError(
                    "Must supply branch_index when starting walk from root BranchPoint"
                )
            return branches_in[branch_index]

        for i, point in enumerate(branches_out):
            if point is nxt:
                return branches_in[i]  # toggle branch
        raise ValueError("Prev node is not connected to this BranchPoint")
