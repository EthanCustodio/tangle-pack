"""
The segment index: where two manifolds cross, and where a bridge is cut.

:class:`Tangle` registers every adjacent pair of points on every manifold as a
segment in an R-tree, finds candidate crossings by bounding-box query, resolves
them exactly, and cuts :class:`~.Bridge.Bridge` objects at the crossings it is
handed.

Dev Notes:

A periodic point IS a crossing here: ``TangleWorkbench._register_anchors`` puts
ONE anchor at canonical distance ``(0, 0)`` on every branch attached to it
(two stable and two unstable branches at an inversion point). It carries no
curve of its own, so ``create_bridges`` cuts every unstable branch rooted at the
periodic point at that branch's own root.

The Tangle holds no resolved crossings: every crossing lives in the workbench's
IntersectionRegistry. Cutting a bridge needs the two real manifold nodes
bracketing the crossing, and each Intersection carries them itself, as
``unstable_manifold`` (the curve it was detected on) and ``unstable_segment`` (the
store rows of the crossing's two bracketing nodes, captured at resolve time and
never re-read off the mutable segment afterwards). So a registered crossing is
self-sufficient: it pins its own geometry, and ``create_bridges`` takes the
crossings to cut at as an argument rather than reaching back into the index.

``create_bridges`` groups by the unstable manifold OBJECT, not by manifold_key: a
bridge and the manifold it was cut out of legitimately share a key while being
distinct polylines, so key-grouping would let one cut span two curves; the object
is the only unambiguous grouping.

Query-only edges (N1, 2026-10-24). An edge a query-only curve has queried stays
in ``_query_edges`` for the index's life, so on a ``reset=False`` recompute it is
never queried again: neither when its image curve is re-added (before N1 the
purge re-queried the edges only that curve owned, and the registry collapsed the
re-detections) nor when a later indexed branch contains it (the branch neither
indexes nor queries it, as before N1; the base also claimed it, for no effect).
Open, and older than N1: the accumulating ``reset=False`` pass misses real
crossings a fresh pass finds (k28 two blasts, grown one step on each stability:
16 crossings against 30; nested, stable grown one step: 48 against 54), and an
image curve laid over a newly grown stable span is not re-queried against it. A
later slice should decide whether ``reset=False`` re-queries or is retired.

Crossing resolution. Below a
``resolution_ratio`` of a few, a crossing of two near-tangent polylines is
undecidable at the current ``area_cutoff``: the polylines may straddle where the
curves do not cross, or miss a pair of crossings the curves do have. Every
detection is registered with its ratio and chord deviation, and the registry
decides which are crossings of N (design note §A.2: ratio at least the threshold, and in no
run of detections closer than their chord deviation); an unresolved one is reported with the remedy "refine" and read by no path. No
parity vote decides a run, so a recompute keeps N. Deciding them
would need a validated method, and interval arithmetic alone is not enough: the
splitting angle near ``k = -1`` reaches ~1e-31, below double precision, so it
would need multiprecision (and a parameterisation of the curves) as well. That is
research, not engineering, and is accepted as out of reach for now.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable, Optional

import numpy as np
from numpy.typing import NDArray
from rtree import index
from rtree.core import RTreeError

from .Intersection import Intersection, ManifoldKey, Stability
from .FixedPoint import FixedPoint
from .BaseManifold import BaseManifold
from .PointStore import PointStore
from .BranchPoint import BranchPoint
from .Bridge import Bridge
from .geometry import chord_deviation_batch

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Relative threshold shared by the two geometric predicates below. Both compare an
# area-like quantity (a cross product, a 2x2 determinant) that scales as the product
# of the two lengths involved, so dividing that product out makes the threshold a
# bound on sin(angle) between the two directions and the predicates scale invariant.
_EPS_REL = 1e-12

# The deviation floor of the resolution ratio, RELATIVE to the shorter segment: two
# straight segments get a large finite ratio instead of a division by zero, and the
# ratio stays invariant under translation, rotation and scaling.
_DEVIATION_FLOOR_REL = 1e-8


def _edge(end0: int | BranchPoint, end1: int | BranchPoint) -> int:
    """One int per edge, whichever way it is walked, from its two store rows; a
    BranchPoint end is -1 through the +1 offset, never an index."""
    a = -1 if isinstance(end0, BranchPoint) else end0
    b = -1 if isinstance(end1, BranchPoint) else end1
    lo, hi = (a, b) if a <= b else (b, a)
    return ((lo + 1) << 32) | (hi + 1)


def _edge_keys(r0: NDArray[np.int64], r1: NDArray[np.int64]) -> NDArray[np.int64]:
    """:func:`_edge` of many edges at once (-1 = a BranchPoint end)."""
    assert max(r0.max(initial=-1), r1.max(initial=-1)) < 2**31, "an edge key holds a row in 31 bits"
    return ((np.minimum(r0, r1) + 1) << 32) | (np.maximum(r0, r1) + 1)


def _curve_ends(
    manifold: BaseManifold,
) -> tuple[Optional[PointStore], NDArray[np.int64], NDArray[np.int64], Optional[BranchPoint]]:
    """A curve's store, its edges' rows ``(r0, r1)`` in curve order (-1 = the
    root BranchPoint), and that BranchPoint (None when the curve has none)."""
    head, store, rows = manifold._rows()
    ends = rows if head is None else np.concatenate((np.array([-1], dtype=np.int64), rows))
    return store, ends[:-1], ends[1:], head


def _segment_index(stream: Optional[list] = None) -> index.Index:
    """A 2-D rtree, bulk (stream) loaded from ``(id, bounds, None)`` entries if any."""
    properties = index.Property()
    properties.dimension = 2
    if stream:
        return index.Index(stream, properties=properties)
    return index.Index(properties=properties)


@dataclass(slots=True)
class _Segment:
    """
    One edge of a curve: two consecutive store rows, the first of them the
    root's BranchPoint where the curve starts at the periodic point (which has
    no row).
    """

    id: int
    manifold: BaseManifold
    store: PointStore
    row0: int | BranchPoint
    row1: int
    # False for segments that were only queried against the rtree, never
    # inserted into it (iterated bridges) -- _remove_segment must not try to
    # delete those from the rtree.
    in_rtree: bool = True

    @property
    def bounds(self) -> tuple[float, float, float, float]:
        """The bounding box ``(minx, miny, maxx, maxy)`` of the rows ``row0`` and ``row1``."""
        (x0, y0), (x1, y1) = self.coords()
        return (min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))

    def coords(self, xy: Optional[NDArray[np.float64]] = None) -> NDArray[np.float64]:
        """The ``(2, 2)`` coordinates of the two ends: read off the caller's view
        of ``store.xy`` (a batch passes one for all its reads), or one taken here."""
        xy = self.store.xy if xy is None else xy
        return np.array([end.get_point() if isinstance(end, BranchPoint) else xy[end]
                         for end in (self.row0, self.row1)], dtype=np.float64)

    def cdists(self, stability: Stability) -> list[float]:
        """The two ends' cdists as Python floats (a BranchPoint reads its own)."""
        cdist = self.store.cdist
        return [end.get_cdist(stability) if isinstance(end, BranchPoint) else float(cdist[end])
                for end in (self.row0, self.row1)]

    def intersects(self, other: _Segment) -> bool:
        """Whether this segment and ``other`` cross (endpoints included)."""
        p0, p1 = map(tuple, self.coords())
        q0, q1 = map(tuple, other.coords())
        return Tangle._do_segments_intersect((p0, p1), (q0, q1))


class Tangle:
    """
    Spatial index of manifold segments that detects and resolves crossings.

    Two kinds of curve meet the index. An INDEXED curve (a branch) has every
    edge (two consecutive store rows) stored as a _Segment in the rtree of its
    stability, with its owners. A QUERY-ONLY curve (a blast image curve, the
    bridges cut from it, a partial piece) is only queried: each new edge is one
    key in ``_query_edges``, and its _Segment is kept only when its query
    found a crossing pair, which resolution reads.

    Candidate crossings are found by querying a segment's bounding box against
    the OTHER stability's rtree, so only unstable x stable pairs are ever
    formed (CLAUDE.md's fundamental invariant: two curves of one stability
    never cross). They are resolved to true intersection points and cut into
    bridges.

    Attributes:
        _rtrees: one spatial rtree per stability (``"unstable"``, ``"stable"``),
            storing the indexed segments by their bounding boxes.
        _seg_lookup: The _Segment of every sid still read: every indexed
            segment, and every query-only segment that found a pair.
        _manifold_segs: The sids each indexed curve (and each bridge cut from
            one) owns.
        _edge_to_sid: The sid of every indexed edge, by its edge key.
        _query_edges: The edge key of every edge a query-only curve has
            queried, so an image laid over an earlier one is not queried twice.
        _intersecting_segments: Set of frozenset pairs of segment ids that
            cross.
        _store: The one PointStore whose rows every key and record name.
        _next_sid: The next segment id this index draws (ids are per index).

    An edge key (:func:`_edge`) is one int built from the edge's two store rows,
    the same whichever way the edge is walked; both key maps use it. Rows are
    store-relative, so a Tangle indexes the curves of one store (asserted).

    Note:
        The Tangle owns INDEX state only. A resolved crossing lives in the
        workbench's :class:`IntersectionRegistry` and nowhere else:
        :meth:`_resolve_crossing_pair` returns the :class:`Intersection` and the
        caller registers it. Each crossing carries the unstable manifold it was
        found on and the rows of its two bracketing nodes, which is everything
        :meth:`create_bridges` needs to cut at it.
    """

    def __init__(self) -> None:
        """
        Initializes an empty index.

        Attributes:
            _rtrees: one spatial rtree per stability, storing the indexed
                segments by their bounding boxes.
            _edge_to_sid: Maps an indexed edge (its :func:`_edge` key, built
                from its two rows) to its segment id, so a bridge cut from an
                indexed branch (which reuses the branch's rows) claims the
                EXISTING segment rather than being silently dropped.
            _seg_lookup: Dictionary keyed by segment id: every indexed segment,
                and every query-only segment whose query found a pair.
            _manifold_segs: The segment ids each indexed curve owns, a bridge
                cut from one included.
            _seg_manifolds: Reverse of _manifold_segs -- every curve that claims
                a given segment id. A segment stays indexed until its last owner
                is removed.
            _query_edges: One key per edge a query-only curve has queried
                (:func:`_edge`, the same key as ``_edge_to_sid``'s).
                It is never forgotten: a query-only curve holds no record, so
                there is nothing to release.
            _intersecting_segments: set of frozensets containing pairs of
                segment ids that intersect
            _processed_pairs: pairs already handed to _resolve_crossing_pair, so a
                crossing is resolved (and registered) exactly once
            _store: The PointStore of the first curve added, None until then.
            _next_sid: The next segment id this index draws. Ids are drawn per
                index, in order, so two builds in one process number their
                segments alike; a pair's operand order follows its two ids.
        """

        self._next_sid = 0
        self._rtrees = {"unstable": _segment_index(), "stable": _segment_index()}
        self._edge_to_sid: dict[int, int] = {}
        self._store: Optional[PointStore] = None
        self._seg_lookup: dict[int, _Segment] = {}
        self._manifold_segs: defaultdict[BaseManifold, set[int]] = defaultdict(set)
        self._seg_manifolds: defaultdict[int, set[BaseManifold]] = defaultdict(set)
        self._query_edges: set[int] = set()

        self._intersecting_segments: set[frozenset[int]] = set()
        self._processed_pairs: set[frozenset[int]] = set()

    @staticmethod
    def _key_of(seg: _Segment) -> ManifoldKey:
        """The ManifoldKey of the segment's manifold (required on every curve)."""
        key = seg.manifold.manifold_key
        assert key is not None, (
            "every curve indexed in a Tangle -- manifold or bridge -- must carry "
            "its manifold_key, so the crossings found on it record which two "
            "branches they sit on"
        )
        return key

    def clear_all(self) -> None:
        """
        Completely clear all Tangle state, useful when recomputing everything from scratch.

        The segment-id counter is kept: a recompute draws ids after the old ones.
        """
        self._rtrees = {"unstable": _segment_index(), "stable": _segment_index()}
        self._edge_to_sid.clear()
        self._seg_lookup.clear()
        self._manifold_segs.clear()
        self._seg_manifolds.clear()
        self._query_edges.clear()
        self._intersecting_segments.clear()
        self._processed_pairs.clear()
        self._store = None

    def resolve_crossings(self) -> list[Intersection]:
        """
        Resolve every detected crossing pair into an Intersection.

        Each crossing is found from its segment PAIR (frozenset of the two segment
        ids), so a single segment may participate in many crossings without any
        being overwritten. Every pair is unstable x stable: the index never pairs
        two segments of one stability (see the class docstring).

        Returns:
            The freshly resolved crossings, for the caller to register. The Tangle
            keeps no copy: the registry is the single source of truth.
        """
        resolved = (self._resolve_crossing_pair(pair) for pair in list(self._intersecting_segments))
        return [intersection for intersection in resolved if intersection is not None]

    def _resolve_crossing_pair(
        self, seg_id_pair: frozenset[int]
    ) -> Optional[Intersection]:
        """
        Resolve one detected crossing pair into an Intersection.

        Handles all the shared bookkeeping: the true crossing point, the cdists
        interpolated at that point (not the segment midpoint), the two branch
        keys, and the unstable manifold and the rows of the bracketing nodes a
        later cut needs. The result is returned, not stored -- the caller adds it
        to the registry.

        Returns:
            The new Intersection, or None when the pair was already processed,
            references a removed segment (dropped as stale), or is near-parallel
            and so has no well-conditioned crossing point.
        """
        if seg_id_pair in self._processed_pairs:
            return None

        seg1_id, seg2_id = tuple(seg_id_pair)

        if seg1_id not in self._seg_lookup or seg2_id not in self._seg_lookup:
            # drop stale pair and skip
            self._intersecting_segments.discard(seg_id_pair)
            return None

        seg_1, seg_2 = self._seg_lookup[seg1_id], self._seg_lookup[seg2_id]

        assert seg_1.manifold.stability != seg_2.manifold.stability, seg_id_pair

        found = self._find_true_intersection(seg_1, seg_2)
        if found is None:
            # near-parallel pair (already logged): discard it rather than let it
            # abort the whole intersection pass
            self._intersecting_segments.discard(seg_id_pair)
            self._processed_pairs.add(seg_id_pair)
            return None
        point, sin_theta = found

        u_seg, s_seg = (
            (seg_1, seg_2) if seg_1.manifold.stability == "unstable" else (seg_2, seg_1)
        )

        ratio, deviation = self._resolution_ratio(u_seg, s_seg, sin_theta)
        assert u_seg.store is s_seg.store, f"the bracket ends of {seg_id_pair} span two stores"
        # cdist interpolated at the true crossing (not the segment midpoint)
        intersection = Intersection.from_segments(
            coords=tuple(point),
            unstable_cdist=self._cdist_at_point(u_seg, point),
            stable_cdist=self._cdist_at_point(s_seg, point),
            seg1_id=seg1_id,
            seg2_id=seg2_id,
            manifold_a_key=Tangle._key_of(u_seg),
            manifold_b_key=Tangle._key_of(s_seg),
            unstable_manifold=u_seg.manifold,
            # The ORIGINAL endpoints' rows, captured now: a cut at this crossing
            # brackets it with these two rows, read off the crossing, so a
            # separator spliced in between them later does not move the bracket.
            unstable_segment=(u_seg.row0, u_seg.row1),
            stable_segment=(s_seg.row0, s_seg.row1),
            store=u_seg.store,
            crossing_sign=self._crossing_sign(u_seg, s_seg),
            resolution_ratio=ratio,
            chord_deviation=deviation,
        )

        self._processed_pairs.add(seg_id_pair)

        return intersection

    @staticmethod
    def _crossing_sign(u_seg: _Segment, s_seg: _Segment) -> int:
        """
        Handedness of one unstable x stable crossing.

        The sign of ``cross(unstable_direction, stable_direction)`` with both
        directions taken in INCREASING canonical distance -- i.e. pointing away
        from their anchoring periodic point. That is exactly the datum the planar
        arrangement needs: it fixes the counter-clockwise cyclic order of the four
        manifold rays leaving the crossing without any angle sorting (see
        :attr:`~.Intersection.Intersection.crossing_sign`).

        Segment endpoints are stored root to tail, which is increasing cdist, but
        the cdists are read rather than assumed so a reversed or hand-built
        segment cannot silently flip the sign of a whole tangle.

        Args:
            u_seg: The unstable segment of the crossing pair.
            s_seg: The stable segment.

        Returns:
            ``+1`` or ``-1``; ``0`` only for a degenerate (zero-length or exactly
            parallel) pair, which ``_find_true_intersection`` has already
            discarded by the time this is called.
        """
        def direction(seg: _Segment, stability: Stability) -> NDArray[np.float64]:
            a, b = seg.coords()
            ca, cb = seg.cdists(stability)
            vec = b - a
            return -vec if cb < ca else vec

        u_dir = direction(u_seg, "unstable")
        s_dir = direction(s_seg, "stable")
        cross = float(u_dir[0] * s_dir[1] - u_dir[1] * s_dir[0])
        return int(np.sign(cross))

    @staticmethod
    def _resolution_ratio(
        u_seg: _Segment, s_seg: _Segment, sin_theta: float
    ) -> tuple[float, float]:
        """
        How far a crossing stands above the polygon noise of its two polylines.

        ``min(l_u, l_s) * |sin(theta)|`` is how far the two segments separate
        across the shorter one; the chord deviation ``h = 1.5 * area / l`` of each
        (the refiner's curvature area, with the curve's neighbouring points, one
        side only where there is none) is how far the true curve may stray from
        its segment. Their ratio is ``min(l) * |sin| / max(h_u, h_s)``: well above
        one the crossing is a fact about the curves, near or below one it may be
        two near-tangent polylines straddling each other. The deviation is
        floored at ``_DEVIATION_FLOOR_REL * min(l)`` so straight segments stay
        finite and the ratio stays similarity invariant.

        Args:
            u_seg: The unstable segment of the crossing pair.
            s_seg: The stable segment.
            sin_theta: ``|sin|`` of the crossing angle.

        Returns:
            The resolution ratio, and ``max(h_u, h_s)`` unfloored: the chord
            deviation a neighbouring detection is measured against.
        """
        segments = (u_seg, s_seg)
        ends = np.array([seg.coords() for seg in segments])
        p0, p1 = ends[:, 0], ends[:, 1]
        left, right = zip(*(Tangle._outer_neighbours(seg) for seg in segments))
        deviation = float(chord_deviation_batch(p0, p1, np.array(left), np.array(right)).max())
        shorter = float(np.hypot(*(p1 - p0).T).min())
        floor = _DEVIATION_FLOOR_REL * shorter
        return shorter * sin_theta / max(deviation, floor), deviation

    @staticmethod
    def _outer_neighbours(
        seg: _Segment,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """The curve points just before and just after a segment, NaN for none.

        They are the chain's neighbours, read past the curve's own ends; the
        chain's anchor is the point before its first. A BranchPoint end (the
        periodic point) counts as having no neighbour: the curve on its far side
        is another branch.
        """
        store = seg.store

        def neighbour(end: int | BranchPoint, step: int) -> Optional[NDArray[np.float64]]:
            if isinstance(end, BranchPoint):
                return None
            chain = store.chain_of(end)
            return None if chain is None else chain.xy_at(store.position_of(end) + step)

        missing = np.array([np.nan, np.nan])
        return tuple(
            missing if xy is None else xy
            for xy in (neighbour(seg.row0, -1), neighbour(seg.row1, 1))
        )

    def _find_true_intersection(
        self, seg1: _Segment, seg2: _Segment
    ) -> Optional[tuple[NDArray[np.float64], float]]:
        """
        Solve two crossing segments for their true point of intersection.

        The 2x2 determinant scales like ``|AB| * |CD|``, so the near-parallel
        test divides that product out: with an absolute threshold a real
        crossing of two 1e-7-long segments looked degenerate and aborted the
        whole intersection pass with an exception.

        Args:
            seg1: One of the two segments.
            seg2: The other segment.

        Returns:
            The crossing point and ``|sin|`` of the crossing angle, or None if the
            two directions are parallel to within ``_EPS_REL`` (the caller
            discards such a pair).
        """
        a, b = seg1.coords()
        c, d = seg2.coords()
        # 2x2 matrix
        M = np.vstack((b - a, c - d)).T  # [[Bx-Ax, Cx-Dx], [By-Ay, Cy-Dy]]
        rhs = c - a

        det = float(np.linalg.det(M))
        scale = float(np.linalg.norm(b - a) * np.linalg.norm(d - c))
        if abs(det) <= _EPS_REL * scale:
            logger.warning(
                "Discarding near-parallel segment pair %s (%s) x %s (%s): "
                "relative det %.3g",
                seg1.id,
                Tangle._key_of(seg1),
                seg2.id,
                Tangle._key_of(seg2),
                det / scale if scale > 0.0 else float("inf"),
            )
            return None

        t, s = np.linalg.solve(M, rhs)  # t on AB, s on CD
        return a + t * (b - a), abs(det) / scale

    def add_manifold(
        self,
        manifold: BaseManifold,
        index_segments: bool = True,
        detect_crossings: bool = True,
    ) -> set[int]:
        """
        Registers a manifold's new edges and detects their crossings.

        Args:
            manifold (BaseManifold): The manifold to register.
            index_segments (bool): If True (default) each segment is recorded
                with its owner and inserted into the rtree of its stability so
                later additions can find it. Pass False for a blast image curve:
                it is unstable manifold, so the only curves it can cross are
                stable ones, and those are already in the index. It is queried
                only, and holds no record of its own (see the Note).
            detect_crossings (bool): If False, genuinely new segments are
                recorded without being queried against the rtree. Used when
                cutting bridges out of an indexed branch: every crossing on that
                stretch of curve was already found on the parent's segments, so
                querying again would only duplicate it.

        Returns:
            The ids of the segments the manifold added: empty when every edge of
            it was already registered, by an indexed curve or by an earlier
            query-only one.

        Note:
            An indexed edge that is already registered (a bridge reuses its
            parent's points) is mapped to the existing segment id and that id is
            recorded under ``manifold`` as well, so one segment may have several
            owners.

            A query-only curve (``index_segments=False`` with crossings
            detected) writes one key per new edge into ``_query_edges``, and a
            _Segment into ``_seg_lookup`` only when its query found a pair: the
            pair's resolution reads it. It owns nothing, so it has nothing to
            release.
        """
        if not index_segments and detect_crossings:
            return self._query_curve(manifold)

        if manifold in self._manifold_segs:
            # copy, because _remove_segment mutates the same set
            for sid in list(self._manifold_segs[manifold]):
                self._remove_segment(sid, manifold)

        new_sids = (
            self._insert_segment(segment, index_in_rtree=index_segments,
                                 detect_crossings=detect_crossings)
            for segment in self._segments_of(manifold)
        )
        return {sid for sid in new_sids if sid is not None}

    def _adopt_store(self, store: Optional[PointStore]) -> None:
        """Record the store on first sight; every later curve must share it."""
        if store is None:
            return
        if self._store is None:
            self._store = store
        assert store is self._store, "one Tangle indexes the rows of one store"

    def _query_curve(self, manifold: BaseManifold) -> set[int]:
        """
        The query-only add, on the curve's rows: one ``intersection_v`` for all its new edges.

        Each new edge takes, in curve order and whether it pairs or not, the
        next sid of this index (:meth:`_new_sids`): a pair's operand order comes
        from its two sid values, so they decide float bits, and a curve whose add
        returns no sid is dropped on a fresh recompute. A _Segment is kept only
        for an edge whose query found a pair.
        """
        store, r0, r1, head = _curve_ends(manifold)
        if store is None:
            return set()
        self._adopt_store(store)
        new = []
        for k, edge in enumerate(_edge_keys(r0, r1).tolist()):
            if edge in self._query_edges:
                continue
            indexed_sid = self._edge_to_sid.get(edge)
            if indexed_sid is not None:
                assert self._seg_lookup[indexed_sid].manifold.stability == manifold.stability, (
                    "an edge is shared only between a curve and a piece cut out of it")
                continue
            self._query_edges.add(edge)
            new.append(k)
        if not new:
            return set()

        sids = list(self._new_sids(len(new)))
        e0, e1 = r0[new], r1[new]
        xy = store.xy
        a, b = xy[np.maximum(e0, 0)], xy[e1]
        if head is not None:
            a[e0 < 0] = head.get_point()
        rtree = self._rtrees["stable" if manifold.stability == "unstable" else "unstable"]
        candidates, counts = rtree.intersection_v(np.minimum(a, b), np.maximum(a, b))
        candidates, stops = candidates.tolist(), np.cumsum(counts).tolist()
        for k in np.flatnonzero(counts).tolist():
            edge = (tuple(a[k]), tuple(b[k]))
            pairs = {
                frozenset((sids[k], cand_id))
                for cand_id in candidates[stops[k] - int(counts[k]):stops[k]]
                if Tangle._do_segments_intersect(
                    edge, tuple(map(tuple, self._seg_lookup[cand_id].coords(xy))))
            }
            if pairs:
                row0 = head if e0[k] < 0 else int(e0[k])
                self._seg_lookup[sids[k]] = _Segment(
                    sids[k], manifold, store, row0, int(e1[k]), in_rtree=False)
                self._intersecting_segments |= pairs
        return set(sids)

    def add_manifolds(self, manifolds: list[BaseManifold]) -> None:
        """
        Index every segment of several manifolds at once.

        On a freshly cleared Tangle the segments are handed to the rtree as one
        bulk (stream) load, which is far faster than inserting them one at a
        time — the per-segment inserts dominated the runtime of every full
        recompute. On a Tangle that already holds segments this falls back to
        ``add_manifold`` for each manifold.

        Args:
            manifolds (list[BaseManifold]): The manifolds to index.
        """
        if self._seg_lookup:
            for manifold in manifolds:
                self.add_manifold(manifold)
            return

        new_segments = []
        for manifold in manifolds:
            for segment in self._segments_of(manifold):
                edge_key = _edge(segment.row0, segment.row1)
                existing_sid = self._edge_to_sid.get(edge_key)
                if existing_sid is not None:
                    # Same edge reached from a second manifold: claim the segment
                    # that is already indexed rather than dropping it.
                    self._claim_segment(existing_sid, manifold)
                    continue

                segment.id = self._new_sids(1)[0]
                self._edge_to_sid[edge_key] = segment.id
                self._seg_lookup[segment.id] = segment
                self._claim_segment(segment.id, segment.manifold)
                new_segments.append(segment)

        if not new_segments:
            return

        for stability in self._rtrees:
            self._rtrees[stability] = _segment_index(
                [(seg.id, seg.bounds, None) for seg in new_segments if seg.manifold.stability == stability]
            )

        # Each unstable segment is queried against the stable index, so every
        # u x s pair is tested exactly once.
        for segment in new_segments:
            if segment.manifold.stability != "unstable":
                continue
            for candidate_id in self._rtrees["stable"].intersection(segment.bounds):
                if segment.intersects(self._seg_lookup[candidate_id]):
                    self._intersecting_segments.add(frozenset((segment.id, candidate_id)))

    def populate_intersections(self, new_sids: set[int]) -> list[Intersection]:
        """
        Resolve only the crossing pairs that involve one of ``new_sids``.

        Use this after add_manifold() when incrementally adding a blast image
        curve rather than rebuilding the full tangle from scratch.

        Args:
            new_sids: The segment ids that add_manifold() just returned.

        Returns:
            Every new detection, the unresolved ones included (which of them are
            crossings of N is the registry's decision, over all its detections),
            in unstable cdist order: the caller registers them in this order,
            never in that of a set of segment ids, so ids are reproducible.
        """
        resolved = (
            self._resolve_crossing_pair(pair)
            for pair in list(self._intersecting_segments)
            if pair & new_sids
        )
        return sorted(
            (intersection for intersection in resolved if intersection is not None),
            key=lambda ix: (ix.unstable_cdist, ix.stable_cdist),
        )

    # ------------- internal helpers -----------------
    def _new_sids(self, n: int) -> range:
        """The next ``n`` segment ids of this index, in order."""
        sids = range(self._next_sid, self._next_sid + n)
        self._next_sid += n
        return sids

    def _claim_segment(self, sid: int, manifold: BaseManifold) -> None:
        """Record ``manifold`` as an owner of segment ``sid`` (idempotent)."""
        self._manifold_segs[manifold].add(sid)
        self._seg_manifolds[sid].add(manifold)

    def _insert_segment(
        self,
        seg: _Segment,
        index_in_rtree: bool = True,
        detect_crossings: bool = True,
    ) -> Optional[int]:
        """
        Registers a segment, detects its crossings, and (optionally) inserts
        it into the rtree.

        An edge that is already registered is NOT duplicated: the existing
        segment id is claimed by ``seg.manifold``. A bridge is cut
        out of a live manifold and so reuses that manifold's points; dropping
        those edges would leave ``_manifold_segs[bridge]`` empty and make every
        ``for_manifold=`` filter and segment lookup blind to it.

        Parameters:
            seg (_Segment): segment to be inserted
            index_in_rtree (bool): If False the segment is registered and its
                crossings against the indexed segments are detected, but it is
                not itself made findable. See add_manifold for when this is
                the right call.
            detect_crossings (bool): If False the rtree is not queried for
                crossings of this segment. See add_manifold.

        Returns:
            The freshly allocated segment id, or None when the edge was already
            registered (and so only claimed).
        """
        edge_key = _edge(seg.row0, seg.row1)
        if edge_key in self._query_edges:
            # a query-only curve already queried it: as before N1, neither indexed
            # nor queried (the base also claimed it, for no effect; the
            # reset=False path, see Dev Notes)
            return None
        existing_sid = self._edge_to_sid.get(edge_key)
        if existing_sid is not None:
            # already indexed -> do NOT duplicate, just add an owner
            existing = self._seg_lookup[existing_sid]
            assert existing.manifold.stability == seg.manifold.stability, (
                "an edge is shared only between a manifold and a piece cut out of "
                "it, which have the same stability"
            )
            self._claim_segment(existing_sid, seg.manifold)
            return None

        # choose a new id there for a new segment
        sid = self._new_sids(1)[0] if seg.id is None else seg.id

        seg.id = sid  # assign the id to the segment

        # insert segment into the rtree and dictionary
        if index_in_rtree:
            self._rtrees[seg.manifold.stability].insert(sid, seg.bounds)
        seg.in_rtree = index_in_rtree
        self._edge_to_sid[edge_key] = sid
        self._seg_lookup[sid] = seg
        self._claim_segment(sid, seg.manifold)

        if not detect_crossings:
            return sid

        other_stability = "stable" if seg.manifold.stability == "unstable" else "unstable"
        for cand_id in self._rtrees[other_stability].intersection(seg.bounds):
            if seg.intersects(self._seg_lookup[cand_id]):
                self._intersecting_segments.add(frozenset((sid, cand_id)))

        return sid

    def release_manifold(self, manifold: BaseManifold) -> None:
        """
        Drop every claim ``manifold`` has on the index.

        A manifold that is discarded (most often a bridge dropped by
        :meth:`TangleWorkbench.clear_bridges`) must not stay an owner of its
        segments: it would keep them alive in ``_seg_lookup`` for the lifetime of
        the Tangle and keep the dead object reachable from ``_seg_manifolds``.
        Segments it shares with a still-live owner -- a bridge's segments belong
        to the unstable manifold it was cut out of as well -- are kept indexed for
        that owner; only segments whose last owner this was are dropped.

        Args:
            manifold: The manifold giving up its claims. Unknown manifolds are a
                no-op.

        Note:
            A manifold whose owner ran :meth:`clear_all` first (the usual case --
            every production ``rebuild_bridges`` follows a recompute, which clears
            the whole index) is already unknown here, so this is then a no-op. It
            does the real work on the direct ``clear_bridges`` path.
        """
        for sid in list(self._manifold_segs.get(manifold, ())):
            self._remove_segment(sid, manifold)
        self._manifold_segs.pop(manifold, None)

    def _remove_segment(self, sid: int, manifold: BaseManifold):
        """
        Drop one owner's claim on a segment, and the segment itself once the
        last owner is gone.

        A segment id may belong to several manifolds (a bridge shares its
        parent's points), so dropping one owner must leave the segment indexed
        for the others.

        Parameters:
            sid (int): segment id
            manifold (BaseManifold): The owner giving up the segment.
        """
        # Check if segment exists before trying to remove it
        if sid not in self._seg_lookup:
            return

        owners = self._seg_manifolds.get(sid, set())
        owners.discard(manifold)
        self._manifold_segs[manifold].discard(sid)
        if owners:
            return  # still referenced -> keep it indexed

        self._seg_manifolds.pop(sid, None)
        seg = self._seg_lookup.pop(sid)

        # Safely remove from rtree (query-only segments were never in it)
        if seg.in_rtree:
            try:
                self._rtrees[seg.manifold.stability].delete(sid, seg.bounds)
            except RTreeError:
                # Segment may have already been deleted from the rtree; anything
                # else deserves to surface rather than be swallowed.
                logger.debug(
                    "rtree delete failed for segment %s", sid, exc_info=True
                )

        self._edge_to_sid.pop(_edge(seg.row0, seg.row1), None)

    def _segments_of(self, manifold: BaseManifold):
        """
        Yield the segments of a manifold, from the root up to and including the
        segment that ends at ``manifold.tail``.

        ``tail`` bounds the indexed extent: for an untrimmed manifold it is the true
        end (kept current by ``BaseManifold._find_tail`` after every growth step), so
        the whole manifold is indexed; for a trimmed manifold (``tail`` moved to an
        interior point, e.g. by ``trim_stable_manifolds`` or the loom resonance-zone
        trim) only the segments up to the trim are indexed. ``tail is None`` falls
        back to the end of the root's chain.
        """
        store, r0, r1, head = _curve_ends(manifold)
        self._adopt_store(store)
        for row0, row1 in zip(r0.tolist(), r1.tolist()):
            yield _Segment(None, manifold, store, head if row0 < 0 else row0, row1)

    @staticmethod
    def _orientation(
        a: tuple[float, float],
        b: tuple[float, float],
        c: tuple[float, float],
        eps_rel: float = _EPS_REL,
    ) -> int:
        """
        Orientation of the ordered triple (a, b, c).

        ``val`` is the cross product of ``b - a`` and ``c - b``, i.e.
        ``|b - a| * |c - b| * sin(angle)``. Comparing it against an ABSOLUTE
        epsilon made the verdict depend on the size of the triangle: on a
        refined manifold (leg lengths ~1e-7) a perfectly real shallow crossing
        produces ``val ~ 1e-17`` and was reported collinear, so the crossing was
        missed. Dividing out the two leg lengths turns the threshold into one on
        ``sin(angle)``, which is scale invariant.

        Args:
            a: First point as an (x, y) pair.
            b: Second point as an (x, y) pair.
            c: Third point as an (x, y) pair.
            eps_rel: Threshold on ``sin(angle)`` below which the triple counts
                as collinear.

        Returns:
            0 if a, b, c are collinear, 1 if clockwise, 2 if counterclockwise.
        """
        ux, uy = b[0] - a[0], b[1] - a[1]
        vx, vy = c[0] - b[0], c[1] - b[1]
        val = uy * vx - ux * vy
        scale = np.hypot(ux, uy) * np.hypot(vx, vy)
        if abs(val) <= eps_rel * scale:
            # covers the exact-zero case, including a degenerate (zero-length) leg
            return 0
        return 1 if val > 0 else 2

    @staticmethod
    def _do_segments_intersect(
        segA: tuple[tuple[float, float], tuple[float, float]],
        segB: tuple[tuple[float, float], tuple[float, float]],
    ) -> bool:
        """
        Using orientation tests to see if segA and segB intersect in a proper point.
        segA, segB = ((x1,y1), (x2,y2)), ((x3,y3), (x4,y4))
        Returns True if there's a proper intersection, False otherwise.
        Ignores collinearity for simplicity (treat collinear as not intersecting).
        """
        p0_seg1, p1_seg1 = segA
        p0_seg2, p1_seg2 = segB

        o1 = Tangle._orientation(p0_seg1, p1_seg1, p0_seg2)
        o2 = Tangle._orientation(p0_seg1, p1_seg1, p1_seg2)
        o3 = Tangle._orientation(p0_seg2, p1_seg2, p0_seg1)
        o4 = Tangle._orientation(p0_seg2, p1_seg2, p1_seg1)

        # General case:
        if o1 != o2 and o3 != o4:
            return True

        return False  # no intersection if collinear or not straddling

    def create_bridges(
        self,
        crossings: Iterable[Intersection],
        for_manifold: Optional[BaseManifold] = None,
        fixed_point: Optional[FixedPoint] = None,
        *,
        cdist_tol: float = 0.0,
    ) -> list[Bridge]:
        """
        Cut unstable manifolds into bridges at the given crossings.

        The crossings to cut at are an INPUT: the Tangle holds no resolved
        crossings of its own (see the class docstring), so the caller passes the
        crossings of N -- or, when only one image is being cut, just those born
        on it. An unresolved detection is never passed: it lies inside a bridge. Each crossing carries the unstable manifold it was found on and the
        rows of the two nodes bracketing it, which is everything a cut needs.

        Crossings are grouped by their parent unstable manifold OBJECT and sorted by
        unstable canonical distance, so bridges never span two curves and a single
        segment that hosts several crossings is handled correctly (each crossing is
        a distinct cut, not deduped by segment id). The grouping is by object rather
        than by ``manifold_key`` because a bridge and the manifold it was cut out of
        legitimately share a key while being distinct polylines.

        Args:
            crossings: The resolved crossings to cut at, each with an assigned
                registry ``id``. An anchor cuts every unstable manifold rooted at
                its periodic point; any other synthetic crossing (no unstable
                segment) is skipped.
            for_manifold: If given, only cut crossings detected on this specific
                curve.
            fixed_point: If given, only build bridges on unstable manifolds that
                emanate from this fixed point. Each Bridge still records its own
                fixed_point, so a global call (fixed_point=None) followed by
                filtering on bridge.fixed_point is equivalent to per-fixed-point
                calls.
            cdist_tol: Unstable canonical-distance span below which the bridge
                between two consecutive crossings is shorter than the curve's
                resolution: it is still cut, with a WARNING (see the Note). Pass
                the registry's ``cdist_tol``.

        Returns:
            List of Bridge objects in unstable canonical-distance order per parent
            manifold, each carrying the registry ids of the two crossings it was
            cut at (its :data:`~.Bridge.BridgeId`). Adjacency is not stored on the
            bridges: consecutive bridges share an endpoint id, so the caller reads
            it off the ids (see ``TangleWorkbench.bridges_at``).

        Note:
            An anchor (a periodic point's own crossing) carries no curve: it is
            on every branch of its point, so it cuts each unstable manifold
            ROOTED there, at the manifold's first segment.

            Two crossings closer together than ``cdist_tol`` are still cut:
            however short, the manifold between them is real, so the bridge is
            cut and the situation is logged at WARNING. A real bridge must never
            disappear quietly.
        """
        crossings = list(crossings)
        # --- 1. Collect crossings grouped by their parent unstable manifold ---
        # entry: (unstable cdist, Intersection, orig_p0, orig_p1), the bracket
        # ends rows (or the root BranchPoint)
        manifold_crossings: dict[
            BaseManifold,
            list[tuple[float, Intersection, int | BranchPoint, int | BranchPoint]],
        ] = defaultdict(list)

        for intersection in crossings:
            manifold = intersection.unstable_manifold
            if manifold is None or intersection.unstable_segment is None:
                continue  # synthetic crossing: nothing to cut

            if for_manifold is not None and manifold is not for_manifold:
                continue

            if fixed_point is not None and manifold.fixed_point is not fixed_point:
                continue

            assert intersection.id is not None, (
                "a crossing must be registered (and so carry its id) before a "
                "bridge is cut at it -- the cut is what sets the bridge endpoints"
            )

            p0, p1 = intersection.unstable_segment
            manifold_crossings[manifold].append(
                (float(intersection.unstable_cdist), intersection, p0, p1)
            )

        # An anchor is ONE crossing on every branch of its periodic point:
        # it cuts each unstable branch rooted there at that branch's own root.
        anchors = {
            (ix.manifold_a_key[0], ix.manifold_a_key[2]): ix for ix in crossings if ix.is_anchor
        }
        for manifold, group in manifold_crossings.items():
            fp, _, orbit_index, _ = manifold.manifold_key
            anchor = anchors.get((fp, orbit_index))
            if anchor is not None and manifold._root is fp.branch_points[orbit_index]:
                head, _store, rows = manifold._rows()
                group.append((0.0, anchor, head, int(rows[0])))

        # --- 2. Build bridges as two points picked from the single manifold ---
        # A bridge is defined by a head and a tail point taken from its parent
        # unstable manifold: the real point just BELOW its first crossing and the
        # real point just ABOVE its last crossing. These are existing manifold points
        # that already carry the iterate links, so iterating a bridge follows the one
        # underlying manifold's links (the store's image column; a periodic-point head's
        # image is the next orbit point) instead of spawning a parallel,
        # slightly-offset set of points. Using fresh
        # 10%-offset "straddle" points was what laid a second polyline over existing
        # curve (the zig-zag), produced overlapping duplicate bridges, and -- because
        # each crossing got two distinct offset points that mapped to near-coincident
        # images -- tripped the unstable x unstable detector all along the tangle.
        cut: list[tuple[BaseManifold, Bridge]] = []  # (parent, bridge)
        for manifold, manifold_group in manifold_crossings.items():
            manifold_group.sort(key=lambda entry: entry[0])

            # The real points bracketing each crossing on the CURRENT polyline. A run
            # of crossings on one segment has no node between them, so each one after
            # the first gets a separator spliced right after the previous crossing's
            # lower bracket: (p0, S1), (S1, S2), ..., (Sn, p1), cdist increasing.
            brackets: list[tuple[int | BranchPoint, int]] = []
            for k, (_cdist, crossing, p0, p1) in enumerate(manifold_group):
                previous = manifold_group[k - 1] if k else None
                if previous is not None and previous[2] == p0 and previous[3] == p1:
                    below = brackets[-1][0]
                    separator = self._insert_crossing_separator(
                        below, p0, p1, previous[1].coords, crossing.coords, manifold
                    )
                    brackets[-1] = (below, separator)
                    brackets.append((separator, p1))
                else:
                    brackets.append((p0, p1))

            # Bridge (k, k+1) runs from the node below crossing k to the node above
            # crossing k+1. It inherits the parent's manifold_key, and its endpoints
            # are the crossings the cut was made at.
            for k, ((cdist_a, first, *_), (cdist_b, second, *_)) in enumerate(
                zip(manifold_group, manifold_group[1:])
            ):
                if cdist_b - cdist_a <= cdist_tol:
                    logger.warning(
                        "Crossings %s and %s are only %.3g apart in unstable cdist "
                        "(below the %.3g tolerance); "
                        "cutting the bridge between them anyway -- it is shorter "
                        "than the canonical-distance resolution and its geometry "
                        "should not be trusted",
                        first.id,
                        second.id,
                        cdist_b - cdist_a,
                        cdist_tol,
                    )
                cut.append((manifold, Bridge(
                    root=brackets[k][0],
                    store=manifold.store,
                    stability=manifold.stability,
                    stretch_param=manifold.stretch_param,
                    fixed_point=manifold.fixed_point,
                    tail=brackets[k + 1][1],
                    branch_index=manifold.branch_index,
                    manifold_key=manifold.manifold_key,
                    first_intersection=first.id,
                    second_intersection=second.id,
                )))

        # --- 3. Register each bridge's edges, by where it was cut ---
        # Done only after every cut, because _insert_crossing_separator may still
        # splice a point into a span an earlier bridge already covers. The edges
        # are the parent's own, so they are neither re-inserted into the rtree nor
        # re-queried: every crossing on them was found on the parent already. A
        # bridge cut from an indexed branch claims its segments; one cut from a
        # query-only curve only records its edges, the separators' included, so a
        # later image laid over them is not queried again.
        for parent, bridge in cut:
            if parent in self._manifold_segs:
                self.add_manifold(bridge, index_segments=False, detect_crossings=False)
            else:
                store, r0, r1, _ = _curve_ends(bridge)
                self._adopt_store(store)
                self._query_edges.update(_edge_keys(r0, r1).tolist())

        return [bridge for _, bridge in cut]

    def _insert_crossing_separator(
        self,
        after: int | BranchPoint,
        p0: int | BranchPoint,
        p1: int,
        coords1: tuple[float, float],
        coords2: tuple[float, float],
        manifold: BaseManifold,
    ) -> int:
        """
        Splice one real point between two crossings that share the ORIGINAL
        segment ``[p0, p1]``, right after ``after``, and return it.

        Two consecutive crossings on the same segment have no manifold point between
        them, so neighbouring bridges cannot bracket them with distinct real points.
        This refines the manifold by inserting a true point at the midpoint of the
        two crossings (its cdist interpolated against the original endpoints),
        restoring one-crossing-per-segment locally.

        Args:
            after: The row the separator follows (the periodic point at the
                chain's start): ``p0`` for the first separator on the segment,
                the previous separator for each later one.
            p0: Lower-cdist endpoint of the shared segment, a row or the
                periodic point.
            p1: Higher-cdist endpoint of the shared segment, a row.
            coords1: Coordinates of the lower crossing on the segment.
            coords2: Coordinates of the higher crossing on the segment.
            manifold: The manifold the segment belongs to.

        Returns:
            The row of the freshly inserted separator point. After the periodic
            point it becomes position 0 of the chain its slot holds.
        """
        mid_xy = 0.5 * (np.asarray(coords1) + np.asarray(coords2))
        store = manifold.store
        cdist = self._cdist_between(store, p0, p1, manifold.stability, mid_xy)
        separator = store.append_one(float(mid_xy[0]), float(mid_xy[1]), cdist)
        position = 0 if isinstance(after, BranchPoint) else store.position_of(after) + 1
        store.chain_of(p1).insert([position], [separator])
        return separator

    def _cdist_at_point(self, seg: _Segment, point: np.ndarray) -> float:
        """
        Interpolate the cdist of `seg`'s manifold at the true crossing `point`.

        Uses the fractional projection of `point` onto the segment, so two
        crossings that share one segment receive distinct cdists (a
        segment-midpoint value would give them one).
        """
        a, b = seg.coords()
        ca, cb = seg.cdists(seg.manifold.stability)
        return Tangle._interpolate(a, b, ca, cb, point)

    @staticmethod
    def _cdist_between(
        store: PointStore,
        p0: int | BranchPoint,
        p1: int | BranchPoint,
        stability: Stability,
        point: np.ndarray,
    ) -> float:
        """
        Interpolate cdist at `point` between the two explicit endpoints p0, p1.

        Endpoint-based variant of `_cdist_at_point`: takes the ends directly
        (rows of ``store``, or the periodic point, read through ``get_cdist``)
        rather than a `_Segment`, so `create_bridges` places a separator against
        the ORIGINAL segment endpoints captured before any point is spliced in.
        """
        (a, ca), (b, cb) = (
            (p.get_point(), p.get_cdist(stability)) if isinstance(p, BranchPoint)
            else (store.xy[p], float(store.cdist[p]))
            for p in (p0, p1)
        )
        return Tangle._interpolate(a, b, ca, cb, point)

    @staticmethod
    def _interpolate(
        a: NDArray[np.float64], b: NDArray[np.float64], ca: float, cb: float, point: np.ndarray
    ) -> float:
        """The cdist at ``point``, interpolated between ends ``a`` (``ca``) and ``b`` (``cb``)."""
        ab = b - a
        denom = float(ab @ ab)
        t = 0.0 if denom == 0.0 else float((np.asarray(point) - a) @ ab) / denom
        t = min(1.0, max(0.0, t))
        return (1 - t) * ca + t * cb
