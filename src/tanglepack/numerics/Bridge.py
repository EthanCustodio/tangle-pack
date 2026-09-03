"""
The arc of unstable manifold between two consecutive crossings.

A :class:`Bridge` is a :class:`~.BaseManifold.BaseManifold` cut out of one
unstable branch, carrying the identity of the cut with it: the branch it lives
on and the registry ids of the two crossings that bound it (its
:data:`BridgeId`).
"""

from __future__ import annotations
from typing import Optional

from .FixedPoint import FixedPoint
from .BaseManifold import BaseManifold
from .Intersection import ManifoldKey, Stability
from .Point import Point

#: Topological identity of a bridge: the registry ids of the two crossings it
#: connects, ordered by increasing unstable canonical distance (i.e. in the
#: unstable dynamical direction). A bridge IS the piece of unstable manifold
#: between two consecutive crossings, so this pair identifies it completely --
#: there is exactly one bridge per pair, on exactly one unstable branch.
BridgeId = tuple[int, int]


class Bridge(BaseManifold):
    """
    A bridge is defined as a segment of manifold that connects two intersection points
    together.

    A bridge is cut out of one unstable manifold between two consecutive crossings,
    so it carries its identity from the cut itself -- never inferred afterwards:

    * ``manifold_key`` -- the unstable branch it lives on, inherited from the curve
      it was cut out of. Every crossing later detected on the bridge (or on its
      forward image) records this key, so canonical distances are always compared
      within one branch.
    * ``first_intersection`` / ``second_intersection`` -- the registry ids of the two
      crossings the cut was made at, ordered by unstable canonical distance (i.e. in
      the unstable dynamical direction).

    Note:
        The root and tail points will not be intersection points so that when a bridge
        is mapped forward it is easy to compute where it intersects.

    Note:
        A leading or trailing piece of an iterated bridge can be bounded by fewer
        than two crossings -- it starts or ends mid-arc, before the first (or after
        the last) crossing on the image. Such a piece is not a bridge by definition,
        but it is still a real stretch of unstable manifold whose dynamics the blast
        frontier must carry forward, so it is kept with ``None`` for the missing
        endpoint(s) and reports :attr:`partial` as True.

    Note:
        A bridge is a FIXED arc once cut -- its two endpoints are crossings that
        do not move -- so its point walk is memoised (``_memoise_walks``). The
        two operations that can still lay new points inside such an arc (growing
        the parent manifold, re-cutting it against fresh crossings) bump every
        registered bridge's version through the workbench; see the
        :mod:`BaseManifold` Dev Notes.
    """

    _memoise_walks: bool = True

    def __init__(
        self,
        root: Point,
        stability: Stability,
        stretch_param: float,
        fixed_point: FixedPoint,
        tail: Point,
        name: str = "unnamed",
        branch_index: Optional[int] = None,
        *,
        manifold_key: Optional[ManifoldKey],
        first_intersection: Optional[int] = None,
        second_intersection: Optional[int] = None,
    ) -> None:
        """

        Note:
            The root and tail cannot be BranchPoints since we want every bridge
            to have a point on either side of the intersection point so that when
            they are mapped forward we can find the new intersection.

            We want to set the tail first so that the BaseManifold __init__ does
            not get the tail via walk.

        Args:
            root (Point): The root of the bridge.
            stability (Stability): Stability of the parent curve.
            stretch_param (float): The parent's stretch parameter.
            fixed_point (FixedPoint): The fixed point the parent emanates from.
            tail (Point): The tail of the bridge.
            name (str, optional): Name of the bridge. Defaults to 'unnamed'.
            branch_index (Optional[int]): Branch of an inversion fixed point.
            manifold_key (Optional[ManifoldKey]): Keyword-only and REQUIRED -- the
                unstable branch this bridge lives on (see the class docstring).
            first_intersection (Optional[int]): Registry id of the crossing at the
                bridge's low-cdist end. ``None`` only for a partial piece.
            second_intersection (Optional[int]): Registry id of the crossing at the
                bridge's high-cdist end. ``None`` only for a partial piece.

        Raises:
            ValueError: A tail must be specified to construct a bridge.
        """
        self._check_input_types(root, tail)

        super().__init__(
            root,
            stability,
            stretch_param,
            fixed_point,
            name,
            tail,
            branch_index,
            manifold_key=manifold_key,
        )

        self.iterated: bool = False
        self.first_intersection: Optional[int] = first_intersection
        self.second_intersection: Optional[int] = second_intersection

    @property
    def id(self) -> Optional[BridgeId]:
        """
        This bridge's :data:`BridgeId`, or ``None`` when it is :attr:`partial`.

        Derived from the endpoint ids set by the cut, never stored separately, so
        it can never disagree with them. Genealogy is derived from this id (see
        :meth:`TangleWorkbench.image_bridges`) rather than stored as object
        links, which is why a bridge carries no parent/children of its own.
        """
        if self.first_intersection is None or self.second_intersection is None:
            return None
        return (self.first_intersection, self.second_intersection)

    @property
    def partial(self) -> bool:
        """
        Whether an endpoint of this piece is not a crossing.

        Derived from the endpoint ids rather than stored, so it can never disagree
        with them. True for the leading/trailing piece of an iterated bridge that
        starts or ends mid-arc (see the class docstring); False for every bridge in
        the topological sense.
        """
        return self.first_intersection is None or self.second_intersection is None

    def _check_input_types(self, root: Point, tail: Point) -> None:
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

