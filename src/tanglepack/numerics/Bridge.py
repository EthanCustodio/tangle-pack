"""
The arc of unstable manifold between two consecutive crossings.

A :class:`Bridge` is a :class:`~.BaseManifold.BaseManifold` cut out of one
unstable branch, carrying the identity of the cut with it: the branch it lives
on and the registry ids of the two crossings that bound it (its
:data:`BridgeId`).
"""

from __future__ import annotations
from typing import Optional, Sequence

import numpy as np
from numpy.typing import NDArray

from .FixedPoint import FixedPoint
from .BaseManifold import BaseManifold
from .BranchPoint import BranchPoint
from .Intersection import ManifoldKey, Stability
from .PointStore import Chain, End, PointStore
from .geometry import polyline_midpoint

#: Topological identity of a bridge: the registry ids of the two crossings it
#: connects, ordered by increasing unstable canonical distance (i.e. in the
#: unstable dynamical direction). A bridge IS the piece of unstable manifold
#: between two consecutive crossings, so this pair identifies it completely --
#: there is exactly one bridge per pair, on exactly one unstable branch.
BridgeId = tuple[int, int]


class CurveReleased(RuntimeError):
    """A read of curve a lean blast freed (``TangleWorkbench.release_curves``).

    Raised by every row read of a released bridge, and by every workbench
    operation that would re-read freed curve once a lean blast has run. The
    remedy is to rebuild the session and blast it without ``lean``.
    """


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
        do not move. The two operations that can still lay new points inside such
        an arc (growing the parent manifold, re-cutting it against fresh
        crossings) bump every registered bridge's version through the workbench,
        which :meth:`Trellis.cached_bridge_polyline` keys on; see the
        :mod:`BaseManifold` Dev Notes.
    """

    def __init__(
        self,
        root: End | BranchPoint | int,
        stability: Stability,
        stretch_param: float,
        fixed_point: FixedPoint,
        tail: End | int,
        name: str = "unnamed",
        branch_index: Optional[int] = None,
        *,
        manifold_key: Optional[ManifoldKey],
        first_intersection: Optional[int] = None,
        second_intersection: Optional[int] = None,
        store: Optional[PointStore] = None,
    ) -> None:
        """

        Note:
            The root and tail cannot be BranchPoints since we want every bridge
            to have a point on either side of the intersection point so that when
            they are mapped forward we can find the new intersection.

            We want to set the tail first so that the BaseManifold __init__ does
            not get the tail via walk.

        Args:
            root (End, BranchPoint or int): The root of the bridge: a row of
                ``store`` (or its ``End``), or the periodic point of an anchor
                bridge.
            stability (Stability): Stability of the parent curve.
            stretch_param (float): The parent's stretch parameter.
            fixed_point (FixedPoint): The fixed point the parent emanates from.
            tail (End or int): The tail of the bridge, a row (or its ``End``).
            name (str, optional): Name of the bridge. Defaults to 'unnamed'.
            branch_index (Optional[int]): Branch of an inversion fixed point.
            manifold_key (Optional[ManifoldKey]): Keyword-only and REQUIRED -- the
                unstable branch this bridge lives on (see the class docstring).
            first_intersection (Optional[int]): Registry id of the crossing at the
                bridge's low-cdist end. ``None`` only for a partial piece.
            second_intersection (Optional[int]): Registry id of the crossing at the
                bridge's high-cdist end. ``None`` only for a partial piece.
            store (Optional[PointStore]): The store of the bridge's rows
                (see :class:`~.BaseManifold.BaseManifold`).

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
            store=store,
        )

        self.iterated: bool = False
        # The min_separation a blast dropped this bridge at (loom/Blast.py);
        # -1.0 for a drop by the width guard.
        self.dropped_at: Optional[float] = None
        # The lobe width the width guard carries down the image chain
        # (loom/Blast.py); None until a width blast seeds or inherits it.
        # The saddle whose stable factor contracts it: the fixed point of the
        # stable branch the seed ring's base lies on.
        self.width: Optional[float] = None
        self.width_saddle: Optional[FixedPoint] = None
        self.first_intersection: Optional[int] = first_intersection
        self.second_intersection: Optional[int] = second_intersection

    # What a lean blast leaves of a released bridge (release): class defaults,
    # set on the instance only when it is released.
    released: bool = False
    _test_point: Optional[NDArray[np.float64]] = None
    _point_count: int = 0

    def release(self) -> None:
        """
        Mark the bridge's interior as freed, keeping what readers still need.

        Its ends, id, key, ``iterated``, ``dropped_at``, ``width`` and
        ``width_saddle`` stay; its midpoint node (:attr:`test_point`) and point
        count are stored first.
        Every later row read raises :class:`CurveReleased`. Called by
        ``TangleWorkbench.release_curves`` before the rows go.
        """
        points = self.get_point_array()
        self._test_point = polyline_midpoint(points)
        self._point_count = len(points)
        self.released = True

    @property
    def test_point(self) -> Optional[NDArray[np.float64]]:
        """The bridge's midpoint node, the point zone membership is tested at
        (``TangleSession.classify_bridge``, a blast's ``in_zone``).

        The geometric middle NODE, so a bridge that is a zone's own unstable
        boundary arc lands exactly on that boundary and, containment being
        boundary-inclusive, is attributed to it. None for an empty bridge.
        Gathered afresh on each read (no read is memoised since N6); stored
        when the bridge is released.
        """
        return self._test_point if self.released else polyline_midpoint(self.get_point_array())

    @property
    def point_count(self) -> int:
        """The points of the bridge's polyline, a periodic-point root included;
        stored when it was released."""
        return self._point_count if self.released else len(self.get_point_array())

    def _span(
        self, final_node: Optional[End | int] = None, stop_at_final: bool = True
    ) -> tuple[Optional[BranchPoint], Optional[Chain], int, int]:
        """:meth:`BaseManifold._span`, refused once the bridge is released.

        Raises:
            CurveReleased: The bridge was released by a lean blast.
        """
        if self.released:
            raise CurveReleased(
                f"bridge {self.id} was released by a lean blast and its points are "
                "freed; rebuild without lean to read it"
            )
        return super()._span(final_node, stop_at_final)

    @classmethod
    def joined(cls, pieces: Sequence["Bridge"]) -> "Bridge":
        """
        Consecutive bridges of one curve as ONE arc, from the first piece's root to
        the last piece's tail.

        The trellis cut at a strong pip's orbit drops the crossings past it, and
        the unstable arc across a dropped crossing is one arc of the cut trellis.
        No map call and no registry write: the workbench never holds the
        joined arc.

        Args:
            pieces: Bridges in unstable order, each starting where the previous
                one ends, all on one branch.

        Returns:
            The joined bridge, keyed and bounded by the outer ends of the run.
        """
        first, last = pieces[0], pieces[-1]
        assert all(
            a.second_intersection == b.first_intersection is not None
            and a.manifold_key == b.manifold_key
            for a, b in zip(pieces, pieces[1:])
        ), "joined bridges must be consecutive pieces of one curve"
        return cls(
            first._root,
            first.stability,
            first.stretch_param,
            first.fixed_point,
            last._tail,
            store=first.store,
            branch_index=first.branch_index,
            manifold_key=first.manifold_key,
            first_intersection=first.first_intersection,
            second_intersection=last.second_intersection,
        )

    @property
    def head(self) -> Optional[End | BranchPoint]:
        """The bridge's first end, a read-only alias of :attr:`root`."""
        return self.root

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

    def _check_input_types(self, root: End | BranchPoint | int, tail: End | int) -> None:
        """
        Checks that the tail used to construct the bridge is a point of the
        curve (a row or its ``End``) rather than a BranchPoint.

        Args:
            root (End, BranchPoint or int): Root of the bridge. Adjacent to the
                first intersection point.
            tail (End or int): Tail of the bridge. Adjacent to the second
                intersection point.

        Raises:
            ValueError: A tail must be specified to construct a bridge.
            TypeError: The tail is a BranchPoint, not a row of the curve.
        """

        if tail is None:
            raise ValueError("A tail must be specified to construct a bridge.")

        if isinstance(tail, BranchPoint):
            raise TypeError(f"Tail must be a row of the curve, not {type(tail).__name__}")

