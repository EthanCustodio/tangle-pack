from typing import Literal, Optional
from .FixedPoint import FixedPoint
from .BranchPoint import BranchPoint
from .Point import Point
from .Intersection import ManifoldKey
import numpy as np
import matplotlib.pyplot as plt

"""
Dev Notes:

Include return type hints including hints like Union[np.ndarray, list[Point]]

Implement _traverse() or something similar to reduce the redundancy in 
the array getting operations.
"""


class BaseManifold:
    """
    Implements a structure for keeping track of manifold characteristics. Contains
    methods for walking along manifolds, extracting information from manfiolds, and
    plotting.

    Attributes:
        root (Point or BranchPoint): First point in the manifold.
        tail (Point or BranchPoint): Final point in the manifold. If not specified the
            tail will be set my walking from the root until None is reached.
        stability (string ["unstable", "stable"]): Stability of the manifold.
        stretch_param (float): Amount by which two points on the manifold separate by
            upon a single iteration of the map. Usually denoted 'alpha'.
        fixed_point (FixedPoint): Fixed point that the manifold originates from.
        name (string, optional): Name of the manifold
        branch_index (int): If the manifold is attached to a fixed point with inversion
            this attribute specifies which branch the manifold eminates from.
    """

    def __init__(
        self,
        root: Point | BranchPoint,
        stability: Literal["stable", "unstable"],
        stretch_param: float,
        fixed_point: FixedPoint,
        name="unnamed",
        tail: Optional[Point | BranchPoint] = None,
        branch_index: Optional[int] = None,
        manifold_key: Optional[ManifoldKey] = None,
    ):
        """
        Initializes the manifold.

        Args:
            root (Point or BranchPoint): First point in the manifold.
            stability (Literal[stable, unstable]): Stability of the manifold.
            stretch_param (float): Amount by which two points on the manifold separate by
            upon a single iteration of the map. Usually denoted 'alpha'.
            fixed_point (FixedPoint): Fixed point that the manifold originates from.
            name (str, optional): Name of the manifold. Defaults to "unnamed".
            tail (Optional[Point or BranchPoint], optional): Final point in the manifold.
                If not specified the tail will be set my walking from the root until None is reached.
            branch_index (Optional[int], optional): If the manifold is attached to a
                fixed point with inversion this attribute specifies which branch the
                manifold eminates from. Defaults to None.
        """

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

    def _find_tail(self):
        """Walks until None is reached and set the tail"""

        previous_point = None
        current_point = self.root

        while current_point is not None:
            next_point = self.walk_fwd(previous_point, current_point)
            previous_point, current_point = current_point, next_point

        self.tail = previous_point

    def _walk_nodes(self, final_node=None, stop_at_final=True):
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
    def _stack_points(points, return_nodes):
        """Return the collected nodes as-is, or stacked into an array."""
        if not points:
            return [] if return_nodes else np.array([])
        return points if return_nodes else np.vstack(points)

    def get_point_array(self, final_node=None, return_nodes=False):
        """
        Walks along the manifold in the stability direction and returns either
        a list of Point objects or an array of (x, y) coordinates.

        Parameters:
            final_node (Optional[Point]): If specified, stop walking at this node.
            return_nodes (bool): If True, return a list of Point objects;
                else return an np.ndarray of [x, y]. Defaults to False.

        Returns:
            list[Point] or np.ndarray of shape (N, 2)
        """
        points = [
            node if return_nodes else node.get_point()
            for node in self._walk_nodes(final_node)
        ]
        return self._stack_points(points, return_nodes)

    def get_cdist_array(self, final_node=None, return_nodes=False):
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
        points = [
            node if return_nodes else node.cdist
            for node in self._walk_nodes(final_node, stop_at_final=False)
        ]
        return self._stack_points(points, return_nodes)

    def get_non_iterated_point_array(
        self, num_iterates: int = 1, final_node=None, return_nodes=False
    ):
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
        points = []
        for node in self._walk_nodes(final_node):
            has_iterate = getattr(node, self._iter_method("exists"))
            if not has_iterate(num_iterates):
                points.append(node if return_nodes else node.get_point().ravel())

        return self._stack_points(points, return_nodes)

    def get_non_iterated_cdist_array(self, num_iterates: int = 1, final_node=None):
        """
        Returns the canonical distances of the points that do not have a
        num_iterates-step iterate yet.

        Args:
            num_iterates (int, optional): Number of iterates. Defaults to 1.
            final_node (Optional[Point]): If specified, stop walking at this node.

        Returns:
            np.ndarray of shape (N, 1)
        """
        cdists = []
        for node in self._walk_nodes(final_node):
            has_iterate = getattr(node, self._iter_method("exists"))
            if not has_iterate(num_iterates):
                cdists.append(node.cdist)

        return self._stack_points(cdists, return_nodes=False)

    def get_iterated_point_array(
        self, num_iterates: int = 1, final_node=None, return_nodes=False
    ):
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
        points = []
        for node in self._walk_nodes(final_node):
            has_iterate = getattr(node, self._iter_method("exists"))
            if not has_iterate(num_iterates):
                continue

            iterate = getattr(node, self._iter_method("get"))(num_iterates)
            if iterate is None:
                raise ValueError("Iterate computed incorrectly, NoneType added")

            points.append(iterate if return_nodes else iterate.get_point().ravel())

        return self._stack_points(points, return_nodes)

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

    def plot(self, color="blue", branch_index=None, show_points=False, **kwargs):
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

    def plot_colormap(self):

        pts = self.get_point_array()  # (N, 2) ndarray
        if pts.size == 0:
            raise ValueError("No points available to plot!")

        idx = np.arange(len(pts))
        norm = plt.Normalize(idx.min(), idx.max())
        cmap = "coolwarm"

        # grey polyline so the geometry is clear, then color-by-index scatter
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.plot(pts[:, 0], pts[:, 1], color="0.7", lw=1, zorder=1)

        sc = ax.scatter(
            pts[:, 0],
            pts[:, 1],
            c=idx,
            cmap=cmap,
            norm=norm,
            s=40,
            edgecolor="k",
            zorder=3,
        )

        ax.set_aspect("equal")
        ax.set_title(f"Manifold Plot ({self.stability.capitalize()})")
        plt.tight_layout()
        return ax

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

    def _iter_attr(self) -> str:
        """Return the Point attribute that stores the iterate for this stability."""

        return "next_iterate" if self.stability == "unstable" else "prev_iterate"

    def _iter_method(self, prefix: str):
        """
        Returns the correct function name based on the manifold stablity

        Example:
            get_next_iterate
            check_prev_iterate
        """

        stability = "next" if self.stability == "unstable" else "prev"
        return f"{prefix}_{stability}_iterate"
