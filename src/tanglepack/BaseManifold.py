from typing import Literal, Optional
from .FixedPoint import FixedPoint
from .BranchPoint import BranchPoint
from .Point import Point
import numpy as np
import matplotlib.pyplot as plt


class BaseManifold():

    def __init__(self, root: Point,
                 stability: Literal["stable", "unstable"],
                 stretch_param: float,
                 name = "unnamed",
                 tail: Optional[Point] = None):
        
        self.root = root
        self.tail = tail
        self.stability = stability
        self.stretch_param = stretch_param
        self.name = name


    def get_point_array(self, final_node=None, return_nodes=False, branch_index=None):
        """
        Walks along the manifold in the stability direction and returns either
        a list of Point objects or an array of (x, y) coordinates.

        Parameters:
            final_node (Optional[Point]): If specified, stop walking before this node.
            return_nodes (bool): If True, return list of Point objects; else return np.ndarray of [x, y].

        Returns:
            list[Point] or np.ndarray of shape (N, 2)
        """

        #TODO implement caching in this method

        points = []
        prev = None
        current = self.root

        while current is not None and current != final_node:
            points.append(current if return_nodes else current.get_point())

            next_node = self.walk_fwd(prev, current, branch_index=branch_index)

            prev, current = current, next_node

        if not points:
            return np.array([]) if not return_nodes else []

        return points if return_nodes else np.vstack(points)


    def walk_fwd(self, prev: Optional[Point], node: Point, branch_index: Optional[int] = None) -> Optional[Point]:  
        """
        Return the next point along the manifold walking away from the fixed point.
        If `node` is a BranchPoint, we exit on the other branch of the
        same stability type we entered on.
        `prev` is the point we just came from (None at the root).
        """
        if isinstance(node, BranchPoint):
            return self._branch_forward(prev, node, branch_index)

        # Ordinary point: follow whichever pointer is "forward"
        return node.forward if self.stability == "unstable" else node.backward


    def walk_back(self, nxt: Optional[Point], node: Point, branch_index: Optional[int] = None) -> Optional[Point]:
        """
        The inverse of `walk_fwd`: step one link *backward* toward the fixed point.
        `nxt` is the point we are coming from.
        """
        if isinstance(node, BranchPoint):
            return self._branch_backward(nxt, node, branch_index)

        return node.backward if self.stability == "unstable" else node.forward


    def plot(self, color='blue', branch_index=None, show_points=False, **kwargs):
        """
        Plots the manifold points.

        Parameters:
            color (str): Color of the manifold line.
            branch_index (int, optional): Branch index if starting from a branch point.
            show_points (bool): Whether to show individual points.
            **kwargs: Additional kwargs for plt.plot().
        """
        points = self.get_point_array(branch_index=branch_index)
        
        if points.size == 0:
            raise ValueError("No points available to plot!")

        plt.plot(points[:, 0], points[:, 1], color=color, **kwargs)

        if show_points:
            plt.scatter(points[:, 0], points[:, 1], color=color, s=10, alpha=0.6)

        plt.title(f'Manifold Plot ({self.stability.capitalize()})')
        plt.axis('equal')
        plt.show()


    # ---------- internal helpers -------------------------------------
    def _branch_forward(self, prev: Point, bp: BranchPoint, branch_index: Optional[int] = None) -> Point:
        """
        Choose the correct outgoing branch at a BranchPoint when moving
        'forward' along the manifold.
        """

        if self.stability == "stable":
            branches_out = bp.backward_branches
            branches_in  = bp.forward_branches
        else:
            branches_out = bp.forward_branches
            branches_in  = bp.backward_branches

        if prev is None:
            if branch_index is None:
                raise ValueError("Must supply branch_index when starting walk from root BranchPoint")
            return branches_out[branch_index]

        for i, point in enumerate(branches_in):
            if point is prev:
                return branches_out[i]  # toggle branch
        raise ValueError("Prev node is not connected to this BranchPoint")


    def _branch_backward(self, nxt: Point, bp: BranchPoint, branch_index: Optional[int] = None) -> Point:
        """
        Choose the correct outgoing branch when walking *backward*.
        Symmetric to _branch_forward.
        """

        if self.stability == "stable":
            branches_out = bp.backward_branches
            branches_in  = bp.forward_branches
        else:
            branches_out = bp.forward_branches
            branches_in  = bp.backward_branches

        if nxt is None:
            if branch_index is None:
                raise ValueError("Must supply branch_index when starting walk from root BranchPoint")
            return branches_in[branch_index]

        for i, point in enumerate(branches_out):
            if point is nxt:
                return branches_in[i]  # toggle branch
        raise ValueError("Prev node is not connected to this BranchPoint")
    
    