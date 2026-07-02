from __future__ import annotations
from typing import Optional, Literal

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .Intersection import Intersection

from numpy.typing import NDArray

import numpy as np

from .FixedPoint import FixedPoint
from .BranchPoint import BranchPoint
from .DynamicalSystem import DynamicalSystem
from .BaseManifold import BaseManifold
from .Point import Point


class Bridge(BaseManifold):
    """
    A bridge is defined as a segment of manifold that connects two intersection points
    together.

    Note:
        The root and tail points will not be intersection points so that when a bridge
        is mapped forward it is easy to compute where it intersects

    Args:
        BaseManifold (_type_): _description_
    """

    def __init__(
        self,
        root: Point,
        stability: Literal["stable", "unstable"],
        stretch_param: float,
        fixed_point: FixedPoint,
        tail: Point,
        name="unnamed",
        branch_index: Optional[int] = None,
    ):
        """

        Note:
            The root and tail cannot be BranchPoints since we want every bridge
            to have a point on either side of the intersection point so that when
            they are mapped forward we can find the new intersection.

            We want to set the tail first so that the BaseManifold __init__ does
            not get the tail via walk.

        Args:
            root (Point): The root of the bridge.
            tail (Point): The tail of the bridge.
            name (str, optional): Name of the bridge. Defaults to 'None'.

        Raises:
            ValueError: A tail must be specified to construct a bridge.
        """
        self._check_input_types(root, tail)

        super().__init__(
            root, stability, stretch_param, fixed_point, name, tail, branch_index
        )

        self.iterated: bool = False
        self.parent: Optional[Bridge] = None
        self.children: list[Bridge] = []
        self.next_bridge: Optional[Bridge] = None
        self.prev_bridge: Optional[Bridge] = None
        self.first_intersection: Optional[Intersection] = None
        self.second_intersection: Optional[Intersection] = None

        # we likely in here want to have these bridges have a quick reference
        # to the two intersection points that define it.
        # bridges are uniquely defined based on the two intersection points
        # that it connects.
        # A dictionary that has tuples as a key where the tuple is the two intersection
        # points that the bridge connects.
        # I think that storage mechanism should exist inside of the Tangle object

    def _check_input_types(self, root: Point, tail: Point):
        """
        Checks that the tail used to construct the bridge is a Point object
        rather than a BranchPoint.

        Args:
            root (Point): Root of the bridge. Adjacent to the first intersection point.
            tail (Point): Tail of the bridge. Adjacent to the second intersection point.

        Raises:
            ValueError: A tail must be specified to construct a bridge.
            TypeError: Tail must be a Point.
        """

        if tail is None:
            raise ValueError("A tail must be specified to construct a bridge.")

        if not isinstance(tail, Point):
            raise TypeError(f"Tail must be a Point, not {type(tail).__name__}")

    def map_forward(self):
        """
        Not implemented here; iteration is handled by TangleWorkbench.iterate_bridge().
        ManifoldMachine owns the map logic and BranchPoint insertion.
        """
        # if we are going to include this then we have to think about
        # the dependecy hell of combining this with ManifoldMachine.
        # We should probably just use manifold machine to map these forward
        # Although you can imagine that mapping a bridge we might have special things
        # Maybe we just create a method map_bridge which just calles iterate_manifold
        # and does some more behind the scenes work.

        pass
