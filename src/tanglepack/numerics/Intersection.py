"""
A single stable x unstable crossing, and the aliases that name a manifold.

:class:`Intersection` records one crossing -- its coordinates, its canonical
distance along each of the two curves, and the geometry needed to cut a bridge
there. This module is also the single home of the :data:`Stability` and
:data:`ManifoldKey` aliases used across the package.
"""

from __future__ import annotations

from functools import cached_property
from typing import Optional, Literal, TYPE_CHECKING, Union

import numpy as np
from numpy.typing import NDArray

from .precise import two_sum

if TYPE_CHECKING:
    from .BaseManifold import BaseManifold
    from .BranchPoint import BranchPoint
    from .FixedPoint import FixedPoint
    from .PointStore import PointStore

#: The single spelling of a manifold's stability. Every ``stability`` parameter
#: in :mod:`tanglepack.numerics` is annotated with this alias.
Stability = Literal["unstable", "stable"]

# ManifoldKey = (fixed_point, stability, orbit_index, branch_index)
# Identical to the key type used in TangleWorkbench.manifolds.
ManifoldKey = tuple["FixedPoint", Stability, int, int]

#: One end of a bracket: a store row, or the periodic point's BranchPoint, which
#: has no row (D5: it is a view of its FixedPoint). It reads its orbit point and
#: cdist 0.0, and fails loudly as an index, where a -1 sentinel would silently
#: read the last row.
SegmentEnd = Union[int, "BranchPoint"]

#: A ``(lo, hi)`` interval of canonical distance.
Bracket = tuple[float, float]

#: A crossing's ORDER key, ``(cdist, residual)`` as Python floats: a tuple
#: compares the plain float first, so the residual only breaks a plain tie (where
#: it is the precise order). Where plain cdists differ they decide, as they decide
#: the merge's chain order: deep chains hold neighbours one ULP apart whose
#: residuals order them the other way (ManifoldMachine's Dev Notes).
CdistKey = tuple[float, float]


class Intersection:
    """
    Represents a single crossing between the stable and unstable manifolds.

    Two instances at the same coordinates are distinct objects (eq=False keeps
    identity-based equality and hash, so Intersections can live in sets/dicts).

    Attributes:
        coords: Geometric (x, y) of the crossing.
        unstable_cdist: Position along the unstable manifold at this crossing.
        stable_cdist: Position along the stable manifold at this crossing.
        unstable_cdist_lo: The residual of ``unstable_cdist`` (:mod:`~.precise`):
            the precise value is the sum. 0.0 for an anchor or a hand-built one.
        stable_cdist_lo: The residual of ``stable_cdist``.
        seg_ids: The pair of R-tree segment IDs that produced this crossing.
            None for synthetic intersections.
        label: Optional human-readable name.
        unstable_manifold: The unstable curve the crossing was detected on --
            an indexed manifold or an iterated bridge. ``Tangle.create_bridges``
            groups crossings by it, because a bridge and the manifold it was cut
            out of share a ``manifold_key`` but are distinct polylines.
        unstable_segment: The store rows of the two adjacent nodes of
            ``unstable_manifold`` that bracket the crossing, captured at resolve
            time (a :data:`SegmentEnd`: the periodic point's ``BranchPoint`` in
            place of a row). These are the head and tail candidates a cut at
            this crossing uses; a separator spliced into the same segment later
            adds a new row and moves neither of them. None for a synthetic
            crossing.
        stable_segment: The rows of the two adjacent nodes of the stable
            manifold the crossing was found between. With ``unstable_segment``
            they are the crossing's BRACKETS: each canonical distance is known
            only to within its two nodes' cdists, and the map carries those
            brackets exactly, which is how the iterate table is filled
            (:mod:`~.IterateInference`). None for a hand-built synthetic one.
        store: The :class:`~.PointStore.PointStore` the bracket rows belong to
            (rows are store-relative). None for a synthetic crossing.
        crossing_sign: Handedness of the crossing: ``+1`` when the stable
            direction is counter-clockwise from the unstable one, ``-1`` when it
            is clockwise, with both directions taken in INCREASING canonical
            distance (away from the anchor). It is the sign of
            ``cross(unstable_direction, stable_direction)``, so it fixes the
            cyclic order of the four manifold rays leaving this crossing --
            ``(u+, s+, u-, s-)`` counter-clockwise when positive and
            ``(u+, s-, u-, s+)`` when negative -- which is all the planar
            arrangement needs to build faces without sorting any angles (see
            :class:`~tanglepack.topology.Arrangement.Arrangement`). ``0`` marks a
            crossing whose sign was never computed (a hand-built synthetic one);
            the arrangement treats that as positive and logs it.
        resolution_ratio: How well the two polylines resolve this crossing:
            ``min(l_u, l_s) * |sin(theta)| / max(h_u, h_s)``, the separation the
            two segments gain across their shorter length over the larger chord
            deviation of the curves they stand for (see
            ``Tangle._resolution_ratio``). Below a few, the crossing cannot be
            told apart from polygon noise. None for an anchor or a synthetic
            crossing. When a later detection collides with a registered
            crossing, the stored (first) detection's ratio is the one kept.
        chord_deviation: ``max(h_u, h_s)``, the larger chord deviation of the
            two segments, i.e. how far the curves may stray from the polylines
            here. Two detections of one branch pair closer along the unstable
            polyline than the larger of their two deviations are one unresolved
            run, out of N (:mod:`~.IntersectionRegistry`).
            None exactly when the crossing is synthetic.
    """

    def __init__(
        self,
        coords: Optional[tuple[float, float]] = None,
        unstable_cdist: Optional[float] = None,
        stable_cdist: Optional[float] = None,
        seg_ids: Optional[frozenset[int]] = None,
        id: Optional[int] = None,
        label: Optional[str] = None,
        manifold_a_key: Optional[ManifoldKey] = None,
        manifold_b_key: Optional[ManifoldKey] = None,
        unstable_manifold: Optional["BaseManifold"] = None,
        unstable_segment: Optional[tuple[SegmentEnd, SegmentEnd]] = None,
        crossing_sign: int = 0,
        resolution_ratio: Optional[float] = None,
        stable_segment: Optional[tuple[SegmentEnd, SegmentEnd]] = None,
        chord_deviation: Optional[float] = None,
        store: Optional["PointStore"] = None,
        unstable_cdist_lo: float = 0.0,
        stable_cdist_lo: float = 0.0,
    ) -> None:
        """
        Build one crossing. Every field is optional so that a synthetic
        intersection (a periodic-point anchor, a test fixture) can be made
        without geometry; see the class docstring for what each one means.
        """
        self.coords = coords
        self.unstable_cdist = unstable_cdist
        self.stable_cdist = stable_cdist
        self.unstable_cdist_lo = unstable_cdist_lo
        self.stable_cdist_lo = stable_cdist_lo
        self.seg_ids = seg_ids
        self.id = id
        self.label = label
        self.manifold_a_key = manifold_a_key
        self.manifold_b_key = manifold_b_key
        self.unstable_manifold = unstable_manifold
        self.unstable_segment = unstable_segment
        self.stable_segment = stable_segment
        self.store = store
        self.crossing_sign = int(crossing_sign)
        self.resolution_ratio = resolution_ratio
        self.chord_deviation = chord_deviation

    @classmethod
    def from_segments(
        cls,
        coords: tuple[float, float],
        unstable_cdist: float,
        stable_cdist: float,
        seg1_id: int,
        seg2_id: int,
        manifold_a_key: Optional[ManifoldKey] = None,
        manifold_b_key: Optional[ManifoldKey] = None,
        label: Optional[str] = None,
        unstable_manifold: Optional["BaseManifold"] = None,
        unstable_segment: Optional[tuple[SegmentEnd, SegmentEnd]] = None,
        crossing_sign: int = 0,
        resolution_ratio: Optional[float] = None,
        stable_segment: Optional[tuple[SegmentEnd, SegmentEnd]] = None,
        chord_deviation: Optional[float] = None,
        store: Optional["PointStore"] = None,
        unstable_cdist_lo: float = 0.0,
        stable_cdist_lo: float = 0.0,
    ) -> Intersection:
        """Create an Intersection backed by two R-tree segment IDs."""
        return cls(
            coords=coords,
            unstable_cdist=unstable_cdist,
            stable_cdist=stable_cdist,
            seg_ids=frozenset({seg1_id, seg2_id}),
            label=label,
            manifold_a_key=manifold_a_key,
            manifold_b_key=manifold_b_key,
            unstable_manifold=unstable_manifold,
            unstable_segment=unstable_segment,
            crossing_sign=crossing_sign,
            resolution_ratio=resolution_ratio,
            stable_segment=stable_segment,
            chord_deviation=chord_deviation,
            store=store,
            unstable_cdist_lo=unstable_cdist_lo,
            stable_cdist_lo=stable_cdist_lo,
        )

    def cdist_key(self, stability: Stability) -> CdistKey:
        """
        The canonical distance on one side as an order key (:data:`CdistKey`).

        Every order or exact comparison of crossings by cdist reads this: two
        crossings whose plain cdists tie are told apart by their residuals, and
        any other two keep their plain order, the chain's.

        Args:
            stability: Which side.

        Returns:
            ``(cdist, residual)`` Python floats.
        """
        if stability == "unstable":
            return float(self.unstable_cdist), float(self.unstable_cdist_lo)
        return float(self.stable_cdist), float(self.stable_cdist_lo)

    def segment_ends(self, stability: Stability) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """
        The coordinates and canonical distances of one bracket's two ends.

        Args:
            stability: Which bracket: the unstable or the stable segment.

        Returns:
            The ``(2, 2)`` coordinates and the ``(2,)`` cdists of the two ends,
            in captured order. A ``BranchPoint`` end reads its own point and
            cdist 0.0.

        Note:
            Each store column is viewed once per call, never per end.
        """
        ends = self.unstable_segment if stability == "unstable" else self.stable_segment
        xy, cdist = self.store.xy, self.store.cdist
        coords = np.array(
            [xy[end] if isinstance(end, int) else end.get_point() for end in ends],
            dtype=np.float64,
        )
        cdists = np.array(
            [cdist[end] if isinstance(end, int) else 0.0 for end in ends], dtype=np.float64
        )
        return coords, cdists

    @cached_property
    def brackets(self) -> Optional[tuple[Bracket, Bracket]]:
        """
        The sorted unstable and stable cdist brackets, unpadded.

        Read once: a bracket's rows and their cdists are written before the
        crossing is detected and never change after.

        Returns:
            ``((u_lo, u_hi), (s_lo, s_hi))``, or None for a crossing without
            brackets (a synthetic one).
        """
        if self.unstable_segment is None:
            return None
        unstable, stable = (
            sorted(float(c) for c in self.segment_ends(stability)[1])
            for stability in ("unstable", "stable")
        )
        return (unstable[0], unstable[1]), (stable[0], stable[1])

    @cached_property
    def precise_brackets(self) -> Optional[tuple[tuple[float, ...], tuple[float, ...]]]:
        """
        :attr:`brackets` to the precision of the residuals: what the iterate
        inference compares.

        Returns:
            ``(hi, lo)``, two tuples of four Python floats: the normalised
            pairs of ``(u_lo, u_hi, s_lo, s_hi)``, each side's two ends in
            increasing precise order (a ``BranchPoint`` end at 0.0), or None
            for a crossing without brackets. Tuples, so a registry bucket of
            many crossings is one ``np.array`` call.
        """
        if self.unstable_segment is None:
            return None
        ends = [
            sorted(
                (0.0, 0.0) if not isinstance(end, int)
                else two_sum(float(self.store._cdist[end]), float(self.store._cdist_lo[end]))
                for end in segment
            )
            for segment in (self.unstable_segment, self.stable_segment)
        ]
        return tuple(zip(*(end for side in ends for end in side)))

    @property
    def fixed_points(self) -> tuple["FixedPoint", ...]:
        """
        Return the distinct FixedPoint objects involved in this intersection.

        A detected crossing carries both keys; a hand-built synthetic one may
        carry either or neither.

        Returns:
            The distinct fixed points in (a, b) order: one element for a
            homoclinic crossing, two for a heteroclinic one, and the empty
            tuple when no manifold key is set.
        """
        fps: list["FixedPoint"] = []
        for key in (self.manifold_a_key, self.manifold_b_key):
            if key is None:
                continue
            if not any(key[0] is seen for seen in fps):
                fps.append(key[0])
        return tuple(fps)

    @classmethod
    def synthetic(
        cls,
        coords: tuple[float, float],
        unstable_cdist: float,
        stable_cdist: float,
        label: Optional[str] = None,
        manifold_a_key: Optional[ManifoldKey] = None,
        manifold_b_key: Optional[ManifoldKey] = None,
        crossing_sign: int = 0,
    ) -> Intersection:
        """
        Create an Intersection not backed by a detected segment crossing.

        Use this for:
        - The periodic point itself (the anchor at cdist (0, 0), one crossing of
          every branch attached to it)
        - Manually specified turning points
        - Any crossing you want to declare programmatically

        Args:
            coords: Geometric (x, y) of the crossing.
            unstable_cdist: Position along the unstable manifold.
            stable_cdist: Position along the stable manifold.
            label: Optional human-readable name.
            manifold_a_key: Key of the unstable branch this crossing sits on.
            manifold_b_key: Key of the stable branch this crossing sits on.
            crossing_sign: Handedness of the crossing (see the class docstring).
                An anchor takes it from the oriented eigenvectors.

        Returns:
            An Intersection with ``seg_ids`` and ``id`` both None; the id is
            assigned by :meth:`IntersectionRegistry.add`.

        Note:
            The keys matter for anchors: the anchor of every point of a
            periodic orbit sits at cdist (0, 0), so the fixed point and orbit
            index of its keys are what keep them apart in the registry's
            collision test. An anchor is stored under its HOME pair
            ``(u(o, 0), s(o, 0))`` with that pair's crossing sign; branch 1 of
            either stability leaves it the opposite way.
        """
        return cls(
            coords=coords,
            unstable_cdist=unstable_cdist,
            stable_cdist=stable_cdist,
            seg_ids=None,
            label=label,
            manifold_a_key=manifold_a_key,
            manifold_b_key=manifold_b_key,
            crossing_sign=crossing_sign,
        )

    # --- helpers ---
    @property
    def is_anchor(self) -> bool:
        """A periodic point's own crossing: synthetic, at canonical distance (0, 0)."""
        return self.seg_ids is None and self.unstable_cdist == 0.0 and self.stable_cdist == 0.0

    @property
    def is_synthetic(self) -> bool:
        """True if this intersection was not detected from crossing segments."""
        return self.seg_ids is None

    def get_point(self) -> NDArray[np.float64]:
        """Return coords as a (2,) array, consistent with BranchPoint.get_point()."""
        return np.array(self.coords, dtype=np.float64)
