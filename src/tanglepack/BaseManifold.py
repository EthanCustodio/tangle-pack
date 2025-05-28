from typing import Literal, Optional
from .FixedPoint import FixedPoint
from .BranchPoint import BranchPoint
from .Point import Point
import numpy as np
import matplotlib.pyplot as plt


class BaseManifold:

    def __init__(
        self,
        root: Point,
        stability: Literal["stable", "unstable"],
        stretch_param: float,
        name="unnamed",
        tail: Optional[Point] = None,
        branch_index: Optional[int] = None,
    ):

        self.root = root
        self.tail = tail
        self.stability = stability
        self.stretch_param = stretch_param
        self.name = name
        self.branch_index = branch_index

    def get_point_array(self, final_node=None, return_nodes=False):
        """
        Walks along the manifold in the stability direction and returns either
        a list of Point objects or an array of (x, y) coordinates.

        Parameters:
            final_node (Optional[Point]): If specified, stop walking before this node.
            return_nodes (bool): If True, return list of Point objects; else return np.ndarray of [x, y].

        Returns:
            list[Point] or np.ndarray of shape (N, 2)
        """
        if final_node is None:
            final_node = self.tail

        branch_index = self.branch_index

        points = []
        prev = None
        current = self.root

        while current is not None:
            points.append(current if return_nodes else current.get_point())

            if current is final_node:
                break

            next_node = self.walk_fwd(prev, current, branch_index=branch_index)

            prev, current = current, next_node

        if not points:
            return np.array([]) if not return_nodes else []

        return points if return_nodes else np.vstack(points)

    def get_cdist_array(self, final_node=None, return_nodes=False):
        """
        Walks along the manifold in the stability direction and returns either
        a list of Point objects or an array of (x, y) coordinates.

        Parameters:
            final_node (Optional[Point]): If specified, stop walking before this node.
            return_nodes (bool): If True, return list of Point objects; else return np.ndarray of [x, y].

        Returns:
            list[Point] or np.ndarray of shape (N, 2)
        """
        if final_node is None:
            final_node = self.tail

        branch_index = self.branch_index

        points = []
        prev = None
        current = self.root

        while current is not None:
            points.append(current if return_nodes else current.cdist)

            next_node = self.walk_fwd(prev, current, branch_index=branch_index)

            prev, current = current, next_node

        if not points:
            return np.array([]) if not return_nodes else []

        return points if return_nodes else np.vstack(points)

    def _iter_attr(self) -> str:
        """Return the Point attribute that stores the iterate for this stability."""
        return "next_iterate" if self.stability == "unstable" else "prev_iterate"

    def get_non_iterated_point_array(self, final_node=None, return_nodes=False):
        """Returns a list of the iterate points if exists otherwise None"""

        if final_node is None:
            final_node = self.tail

        branch_index = self.branch_index

        points = []
        prev = None
        current = self.root

        while current is not None and current is not final_node:

            if getattr(current, self._iter_attr()) is None:
                points.append(current if return_nodes else current.get_point().ravel())

            next_node = self.walk_fwd(prev, current, branch_index=branch_index)

            prev, current = current, next_node

        if current is final_node:
            points.append(current if return_nodes else current.get_point().ravel())

        if not points:
            return np.array([]) if not return_nodes else []

        return points if return_nodes else np.vstack(points)

    def get_non_iterated_cdist_array(self, final_node=None):
        """Returns a list of the iterate points if exists otherwise None"""

        if final_node is None:
            final_node = self.tail

        branch_index = self.branch_index

        points = []
        prev = None
        current = self.root

        while current is not None and current is not final_node:

            if getattr(current, self._iter_attr()) is None:
                points.append(current.cdist)

            next_node = self.walk_fwd(prev, current, branch_index=branch_index)

            prev, current = current, next_node

        if current is final_node:
            points.append(current.cdist)

        if not points:
            return np.array([])

        return np.vstack(points)

    def get_iterated_point_array(self, final_node=None, return_nodes=False):
        """
        Walks along the manifold in the stability direction and returns either
        a list of Point objects or an array of (x, y) coordinates.

        Parameters:
            final_node (Optional[Point]): If specified, stop walking before this node.
            return_nodes (bool): If True, return list of Point objects; else return np.ndarray of [x, y].

        Returns:
            list[Point] or np.ndarray of shape (N, 2)
        """
        if final_node is None:
            final_node = self.tail

        branch_index = self.branch_index

        points = []
        prev = None
        current = self.root

        while current is not None and current is not final_node:

            if getattr(current, self._iter_attr()) is not None:
                if self.stability == "unstable":
                    points.append(
                        current.next_iterate
                        if return_nodes
                        else current.next_iterate.get_point()
                    )
                else:
                    points.append(
                        current.prev_iterate
                        if return_nodes
                        else current.prev_iterate.get_point()
                    )

            next_node = self.walk_fwd(prev, current, branch_index=branch_index)

            prev, current = current, next_node

        if current is final_node:
            points.append(current if return_nodes else current.get_point().ravel())

        if not points:
            return np.array([]) if not return_nodes else []

        return points if return_nodes else np.vstack(points)

    def walk_fwd(
        self, prev: Optional[Point], node: Point, branch_index: Optional[int] = None
    ) -> Optional[Point]:
        """
        Return the next point along the manifold walking away from the fixed point.
        If `node` is a BranchPoint, we exit on the other branch of the
        same stability type we entered on.
        `prev` is the point we just came from (None at the root).

        Parameters
        ----------
        branch_index : int, optional
        Which branch of a BranchPoint to follow.  If omitted (None),
        `self.branch_index` is used.
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

        Parameters
        ----------
        branch_index : int, optional
        Which branch of a BranchPoint to follow.  If omitted (None),
        `self.branch_index` is used.
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
        # plt.show()

    def plot_colormap(self):

        from matplotlib.collections import LineCollection

        pts = self.get_point_array()  # (N, 2) ndarray
        if pts.size == 0:
            raise ValueError("No points available to plot!")

        idx = np.arange(len(pts))  # 0 … N‑1
        norm = plt.Normalize(idx.min(), idx.max())
        cmap = "coolwarm"  # pick any Matplotlib cmap

        # 1⃣  draw the grey polyline so the geometry is clear
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.plot(pts[:, 0], pts[:, 1], color="0.7", lw=1, zorder=1)

        # 2⃣  colour‑by‑index scatter
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

        # 3⃣  optional: colour‑by‑index segments instead of scatter (uncomment)
        # segs  = np.stack([pts[:-1], pts[1:]], axis=1)
        # lc    = LineCollection(segs, cmap=cmap, norm=norm, linewidth=2)
        # lc.set_array(idx[:-1])
        # ax.add_collection(lc)

        # 4⃣  cosmetics
        # plt.colorbar(sc, ax=ax, pad=0.02, label="Point order (0 → last)")
        ax.set_aspect("equal")
        ax.set_title(f"Manifold Plot ({self.stability.capitalize()})")
        # ax.set_aspect("equal", adjustable="box")
        plt.tight_layout()
        return ax

    # ---------- internal helpers -------------------------------------
    def _branch_forward(
        self, prev: Point, bp: BranchPoint, branch_index: Optional[int] = None
    ) -> Point:
        """
        Choose the correct outgoing branch at a BranchPoint when moving
        'forward' along the manifold.
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
