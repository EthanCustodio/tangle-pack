"""
The segment index: where two manifolds cross, and where a bridge is cut.

:class:`Tangle` registers every adjacent pair of points on every manifold as a
segment in an R-tree, finds candidate crossings by bounding-box query, resolves
them exactly, and cuts :class:`~.Bridge.Bridge` objects at the crossings it is
handed.

Dev Notes:

A periodic point IS a crossing here: ``TangleWorkbench._register_anchors`` puts
ONE anchor at canonical distance ``(0, 0)`` on every branch attached to it (D4:
two stable and two unstable branches at an inversion point). It carries no
curve of its own, so ``create_bridges`` cuts every unstable branch rooted at the
periodic point at that branch's own root.

The Tangle holds no resolved crossings (plan row 2.1): every crossing lives in the
workbench's IntersectionRegistry. What used to make that impossible was geometry --
cutting a bridge needs the two real manifold nodes bracketing the crossing, which
only the segment index knew. Each Intersection now carries that itself, as
``unstable_manifold`` (the curve it was detected on) and ``unstable_segment`` (the
crossing's two bracketing points, captured at resolve time and never re-read off
the mutable segment afterwards). So a registered crossing is self-sufficient: it
pins its own geometry, and ``create_bridges`` takes the crossings to cut at as an
argument rather than reaching back into the index.

``create_bridges`` groups by the unstable manifold OBJECT, not by manifold_key: a
bridge and the manifold it was cut out of legitimately share a key while being
distinct polylines, so key-grouping would let one cut span two curves. Phase 3's
BridgeId may make a key-plus-generation grouping viable; until then the object is
the only unambiguous grouping.

Crossing resolution (2026-10-07; B1, Q13 (a), 2026-10-12). Below a
``resolution_ratio`` of a few, a crossing of two near-tangent polylines is
undecidable at the current ``area_cutoff``: the polylines may straddle where the
curves do not cross, or miss a pair of crossings the curves do have. Every
detection is registered with its ratio and chord deviation, and the registry
decides which are crossings of N (A.2: ratio at least the threshold, and in no
run of detections closer than their chord deviation); an unresolved one is reported with the remedy "refine" and read by no path. There
is no vote: the blast path once kept the median of a run of near-coincident
detections and nothing else did, so a recompute changed N. Deciding them
would need a validated method, and interval arithmetic alone is not enough: the
splitting angle near ``k = -1`` reaches ~1e-31, below double precision, so it
would need multiprecision (and a parameterisation of the curves) as well. That is
research, not engineering, and is accepted as out of reach for now.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass
from itertools import count
from typing import Iterable, Optional

import numpy as np
from numpy.typing import NDArray
from rtree import index
from rtree.core import RTreeError

from .Intersection import Intersection, ManifoldKey
from .FixedPoint import FixedPoint
from .BaseManifold import BaseManifold
from .Point import Point
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
    Helper class to store basic information about
    neighboring pairs of points.
    """

    id: int
    manifold: BaseManifold
    p0: Point | BranchPoint
    p0_seg1: Point | BranchPoint
    # False for segments that were only queried against the rtree, never
    # inserted into it (iterated bridges) -- _remove_segment must not try to
    # delete those from the rtree.
    in_rtree: bool = True

    @property
    def bounds(self) -> tuple[float, float, float, float]:
        """
        Returns the boundaries of the bounding box
        made from self.p0 and self.p0_seg1
        """

        x0, y0 = self.p0.get_point()
        x1, y1 = self.p0_seg1.get_point()

        return (min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))

    def intersects(self, other: "_Segment") -> bool:
        """Whether this segment and ``other`` cross (endpoints included)."""
        p0 = tuple(self.p0.get_point())  # (x, y)
        p0_seg1 = tuple(self.p0_seg1.get_point())
        q0 = tuple(other.p0.get_point())
        q1 = tuple(other.p0_seg1.get_point())

        return Tangle._do_segments_intersect((p0, p0_seg1), (q0, q1))


class Tangle:
    """
    Spatial index of manifold segments that detects and resolves crossings.

    Every adjacent pair of points on an indexed manifold is stored as a
    _Segment in the rtree of its stability. Candidate crossings are found by
    querying a segment's bounding box against the OTHER stability's rtree, so
    only unstable x stable pairs are ever formed (CLAUDE.md's fundamental
    invariant: two curves of one stability never cross). They are resolved to
    true intersection points and cut into bridges.

    Attributes:
        _rtrees: one spatial rtree per stability (``"unstable"``, ``"stable"``),
            storing its segments by their bounding boxes.
        _seg_lookup: Dictionary keyed by unique segment ids that stores
            segments.
        _manifold_segs: Dictionary keyed by manifolds which stores the set of
            all segment ids for that manifold.
        _intersecting_segments: Set of frozenset pairs of segment ids that
            cross.

    Note:
        The Tangle owns INDEX state only. A resolved crossing lives in the
        workbench's :class:`IntersectionRegistry` and nowhere else:
        :meth:`_resolve_crossing_pair` returns the :class:`Intersection` and the
        caller registers it. Each crossing carries the unstable manifold and the
        two bracketing points it was found on, which is everything
        :meth:`create_bridges` needs to cut at it.
    """

    _ids = count(0)  # global generator -> every segment is given a unique int key

    def __init__(self) -> None:
        """

        Initializes the Tangle object with an R-tree for spatial indexing
        of segments, a lookup dictionary for segments, and a mapping
        of manifolds to their segments.

        Attributes:
            _rtrees: one spatial rtree per stability, storing its segments by
                their bounding boxes
            _edge_to_sid: Maps an edge (the frozenset of its two point ids) to the
                id of the segment already registered for it, so a manifold that
                reuses another manifold's points (every bridge does) is given the
                EXISTING segment rather than being silently dropped.
            _seg_lookup: Dictionary keyed by unique segment ids that
                stores segments
            _manifold_segs Dictionary keyed by manifolds which stores
                a set of all segments
            _seg_manifolds: Reverse of _manifold_segs -- every manifold that
                claims a given segment id. A segment stays indexed until its last
                owner is removed.

            _intersecting_segments: set of frozensets containing pairs of
                segment ids that intersect
            _processed_pairs: pairs already handed to _resolve_crossing_pair, so a
                crossing is resolved (and registered) exactly once
        """

        self._rtrees = {"unstable": _segment_index(), "stable": _segment_index()}
        self._edge_to_sid: dict[frozenset[int], int] = {}
        self._seg_lookup: dict[int, _Segment] = {}
        self._manifold_segs: defaultdict[BaseManifold, set[int]] = defaultdict(set)
        self._seg_manifolds: defaultdict[int, set[BaseManifold]] = defaultdict(set)

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
        """
        self._rtrees = {"unstable": _segment_index(), "stable": _segment_index()}
        self._edge_to_sid.clear()
        self._seg_lookup.clear()
        self._manifold_segs.clear()
        self._seg_manifolds.clear()
        self._intersecting_segments.clear()
        self._processed_pairs.clear()

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
        keys, and the unstable manifold and bracketing points a later cut needs.
        The result is returned, not stored -- the caller adds it to the registry.

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
            # The ORIGINAL endpoints, captured now: a cut at this crossing brackets
            # it with these two real nodes, and they must not be re-read off the
            # segment later, after a separator point has been spliced in.
            unstable_segment=(u_seg.p0, u_seg.p0_seg1),
            stable_segment=(s_seg.p0, s_seg.p0_seg1),
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
        def direction(seg: _Segment, stability: str) -> NDArray[np.float64]:
            a, b = seg.p0, seg.p0_seg1
            vec = b.get_point() - a.get_point()
            return -vec if b.get_cdist(stability) < a.get_cdist(stability) else vec

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
        p0 = np.array([seg.p0.get_point() for seg in segments], dtype=float)
        p1 = np.array([seg.p0_seg1.get_point() for seg in segments], dtype=float)
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

        A BranchPoint end (the periodic point) counts as having no neighbour:
        the curve on its far side is another branch.
        """
        missing = np.array([np.nan, np.nan])
        before = (
            None
            if isinstance(seg.p0, BranchPoint)
            else seg.manifold.walk_back(seg.p0_seg1, seg.p0)
        )
        after = (
            None
            if isinstance(seg.p0_seg1, BranchPoint)
            else seg.manifold.walk_fwd(seg.p0, seg.p0_seg1)
        )
        return tuple(
            missing if node is None else node.get_point() for node in (before, after)
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
        p0_seg1, p1_seg1 = seg1.p0.get_point(), seg1.p0_seg1.get_point()
        p0_seg2, p1_seg2 = seg2.p0.get_point(), seg2.p0_seg1.get_point()

        a, b, c, d = map(np.asarray, (p0_seg1, p1_seg1, p0_seg2, p1_seg2))
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
    ) -> int:
        """
        Registers every segment of a manifold and detects its crossings.

        Args:
            manifold (BaseManifold): The manifold to register.
            index_segments (bool): If True (default) each segment is also
                inserted into the rtree of its stability so later additions can
                find it. Pass False for an iterated bridge: a bridge is unstable
                manifold, so the only curves it can cross are stable ones, and
                those are already in the index. Its crossings are found by
                querying; skipping the inserts is what keeps repeated bridge
                iteration fast.
            detect_crossings (bool): If False, genuinely new segments are
                registered without being queried against the rtree. Used when
                cutting bridges out of an already-indexed manifold: every
                crossing on that stretch of curve was already found on the
                parent's segments, so querying again would only duplicate it.

        Returns:
            The number of segments the manifold added: 0 when every edge of it was
            already registered by another owner.

        Note:
            An edge that is already registered (a bridge reuses its parent
            manifold's points) is mapped to the existing segment id and that id
            is recorded under ``manifold`` as well, so one segment may have
            several owners.
        """
        # ---------- A. purge old entries (if any) ----------
        if manifold in self._manifold_segs:
            # copy, because _remove_segment mutates the same set
            for sid in list(self._manifold_segs[manifold]):
                self._remove_segment(sid, manifold)

        # _insert_segment records each id in _manifold_segs itself.
        before = len(self._seg_lookup)
        for segment in self._segments_of(manifold):
            self._insert_segment(
                segment,
                index_in_rtree=index_segments,
                detect_crossings=detect_crossings,
            )
        return len(self._seg_lookup) - before

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
                edge_key = frozenset((id(segment.p0), id(segment.p0_seg1)))
                existing_sid = self._edge_to_sid.get(edge_key)
                if existing_sid is not None:
                    # Same edge reached from a second manifold: claim the segment
                    # that is already indexed rather than dropping it.
                    self._claim_segment(existing_sid, manifold)
                    continue

                segment.id = next(Tangle._ids)
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

    def populate_intersections_for_manifold(
        self, manifold: BaseManifold
    ) -> list[Intersection]:
        """
        Resolve only crossing pairs that involve a segment from `manifold`.

        Use this after add_manifold() when incrementally adding an iterated bridge
        rather than rebuilding the full tangle from scratch.

        Args:
            manifold: The newly added manifold to resolve intersections for.

        Returns:
            Every new detection, the unresolved ones included (which of them are
            crossings of N is the registry's decision, over all its detections),
            in unstable cdist order: the caller registers them in this order,
            never in that of a set of segment ids, so ids are reproducible (INT-7).
        """
        manifold_seg_ids = self._manifold_segs.get(manifold, set())
        resolved = (
            self._resolve_crossing_pair(pair)
            for pair in list(self._intersecting_segments)
            if pair & manifold_seg_ids
        )
        return sorted(
            (intersection for intersection in resolved if intersection is not None),
            key=lambda ix: (ix.unstable_cdist, ix.stable_cdist),
        )

    # ------------- internal helpers -----------------
    def _claim_segment(self, sid: int, manifold: BaseManifold) -> None:
        """Record ``manifold`` as an owner of segment ``sid`` (idempotent)."""
        self._manifold_segs[manifold].add(sid)
        self._seg_manifolds[sid].add(manifold)

    def _insert_segment(
        self,
        seg: _Segment,
        index_in_rtree: bool = True,
        detect_crossings: bool = True,
    ) -> int:
        """
        Registers a segment, detects its crossings, and (optionally) inserts
        it into the rtree.

        An edge that is already registered is NOT duplicated: the existing
        segment id is claimed by ``seg.manifold`` and returned. A bridge is cut
        out of a live manifold and so reuses that manifold's points; dropping
        those edges (the old behaviour) left ``_manifold_segs[bridge]`` empty and
        made every ``for_manifold=`` filter and segment lookup blind to it.

        Parameters:
            seg (_Segment): segment to be inserted
            index_in_rtree (bool): If False the segment is registered and its
                crossings against the indexed segments are detected, but it is
                not itself made findable. See add_manifold for when this is
                the right call.
            detect_crossings (bool): If False the rtree is not queried for
                crossings of this segment. See add_manifold.

        Returns:
            The segment id -- freshly allocated, or the id already registered
            for this edge.
        """
        # edge key defined by the id of two points
        edge_key = frozenset((id(seg.p0), id(seg.p0_seg1)))
        existing_sid = self._edge_to_sid.get(edge_key)
        if existing_sid is not None:
            # already indexed -> do NOT duplicate, just add an owner
            existing = self._seg_lookup[existing_sid]
            assert existing.manifold.stability == seg.manifold.stability, (
                "an edge is shared only between a manifold and a piece cut out of "
                "it, which have the same stability"
            )
            self._claim_segment(existing_sid, seg.manifold)
            return existing_sid

        # choose a new id there for a new segment
        sid = next(Tangle._ids) if seg.id is None else seg.id

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

        self._edge_to_sid.pop(frozenset((id(seg.p0), id(seg.p0_seg1))), None)

    def _segments_of(self, manifold: BaseManifold):
        """
        Yield the segments of a manifold, from the root up to and including the
        segment that ends at ``manifold.tail``.

        ``tail`` bounds the indexed extent: for an untrimmed manifold it is the true
        end (kept current by ``BaseManifold._find_tail`` after every growth step), so
        the whole manifold is indexed; for a trimmed manifold (``tail`` moved to an
        interior point, e.g. by ``trim_stable_manifolds`` or the loom resonance-zone
        trim) only the segments up to the trim are indexed. ``tail is None`` falls
        back to walking to the physical end of the linked list.
        """
        tail = manifold.tail
        prev_point = None
        curr_point = manifold.root

        while curr_point is not None:

            if curr_point is tail:
                break

            next_point = manifold.walk_fwd(prev_point, curr_point)
            if next_point is None:
                break

            yield _Segment(None, manifold, curr_point, next_point)

            prev_point, curr_point = curr_point, next_point

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
        two points bracketing it, which is everything a cut needs.

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
        # entry: (unstable cdist, Intersection, orig_p0, orig_p1)
        manifold_crossings: dict[
            BaseManifold,
            list[tuple[float, Intersection, Point, Point]],
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

        # An anchor is ONE crossing on every branch of its periodic point (D4):
        # it cuts each unstable branch rooted there at that branch's own root.
        anchors = {
            (ix.manifold_a_key[0], ix.manifold_a_key[2]): ix for ix in crossings if ix.is_anchor
        }
        for manifold, group in manifold_crossings.items():
            fp, _, orbit_index, _ = manifold.manifold_key
            anchor = anchors.get((fp, orbit_index))
            if anchor is not None and manifold.root is fp.branch_points[orbit_index]:
                nodes = manifold.get_point_array(return_nodes=True)
                group.append((0.0, anchor, nodes[0], nodes[1]))

        # --- 2. Build bridges as two points picked from the single manifold ---
        # A bridge is defined by a head and a tail point taken from its parent
        # unstable manifold: the real point just BELOW its first crossing and the
        # real point just ABOVE its last crossing. These are existing manifold points
        # that already carry the iterate links, so iterating a bridge follows the one
        # underlying manifold's iterate structure (head.next_iterate, tail.next_iterate)
        # instead of spawning a parallel, slightly-offset set of points. Using fresh
        # 10%-offset "straddle" points was what laid a second polyline over existing
        # curve (the zig-zag), produced overlapping duplicate bridges, and -- because
        # each crossing got two distinct offset points that mapped to near-coincident
        # images -- tripped the unstable x unstable detector all along the tangle.
        all_bridges: list[Bridge] = []
        for manifold, manifold_group in manifold_crossings.items():
            manifold_group.sort(key=lambda entry: entry[0])

            # The real points bracketing each crossing on the CURRENT polyline. A run
            # of crossings on one segment has no node between them, so each one after
            # the first gets a separator spliced right after the previous crossing's
            # lower bracket: (p0, S1), (S1, S2), ..., (Sn, p1), cdist increasing.
            brackets: list[tuple[Point, Point]] = []
            for k, (_cdist, crossing, p0, p1) in enumerate(manifold_group):
                previous = manifold_group[k - 1] if k else None
                if previous is not None and previous[2] is p0 and previous[3] is p1:
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
                all_bridges.append(Bridge(
                    root=brackets[k][0],
                    stability=manifold.stability,
                    stretch_param=manifold.stretch_param,
                    fixed_point=manifold.fixed_point,
                    tail=brackets[k + 1][1],
                    branch_index=manifold.branch_index,
                    manifold_key=manifold.manifold_key,
                    first_intersection=first.id,
                    second_intersection=second.id,
                ))

        # --- 3. Register each bridge's segments under the bridge itself ---
        # Done only after every cut, because _insert_crossing_separator may still
        # splice a point into a span an earlier bridge already covers. The segments
        # are the parent's own, so they are neither re-inserted into the rtree nor
        # re-queried: every crossing on them was found on the parent already.
        for bridge in all_bridges:
            self.add_manifold(bridge, index_segments=False, detect_crossings=False)

        return all_bridges

    def _insert_crossing_separator(
        self,
        after: Point,
        p0: Point,
        p1: Point,
        coords1: tuple[float, float],
        coords2: tuple[float, float],
        manifold: BaseManifold,
    ) -> Point:
        """
        Splice one real point between two crossings that share the ORIGINAL
        segment ``[p0, p1]``, right after ``after``, and return it.

        Two consecutive crossings on the same segment have no manifold point between
        them, so neighbouring bridges cannot bracket them with distinct real points.
        This refines the manifold by inserting a true point at the midpoint of the
        two crossings (its cdist interpolated against the original endpoints),
        restoring one-crossing-per-segment locally.

        Args:
            after: The node the separator follows: ``p0`` for the first separator
                on the segment, the previous separator for each later one.
            p0: Lower-cdist endpoint of the shared segment.
            p1: Higher-cdist endpoint of the shared segment.
            coords1: Coordinates of the lower crossing on the segment.
            coords2: Coordinates of the higher crossing on the segment.
            manifold: The manifold the segment belongs to.

        Returns:
            The freshly inserted separator point.
        """
        mid_xy = 0.5 * (np.asarray(coords1) + np.asarray(coords2))
        cdist = self._cdist_between(p0, p1, manifold.stability, mid_xy)
        separator = Point(float(mid_xy[0]), float(mid_xy[1]), cdist)
        if isinstance(after, BranchPoint):
            after.insert_point_forward(separator, manifold.branch_index)
        else:
            after.insert_point_forward(separator)
        return separator

    def _cdist_at_point(self, seg: _Segment, point: np.ndarray) -> float:
        """
        Interpolate the cdist of `seg`'s manifold at the true crossing `point`.

        Uses the fractional projection of `point` onto the segment, so two
        crossings that share one segment receive distinct cdists (unlike the old
        segment-midpoint value).
        """
        return self._cdist_between(
            seg.p0, seg.p0_seg1, seg.manifold.stability, point
        )

    @staticmethod
    def _cdist_between(
        p0: Point,
        p1: Point,
        stability: str,
        point: np.ndarray,
    ) -> float:
        """
        Interpolate cdist at `point` between the two explicit endpoints p0, p1.

        Endpoint-based variant of `_cdist_at_point`: takes the Point objects
        directly rather than reading them off a `_Segment`. This lets
        `create_bridges` evaluate cdists against the ORIGINAL segment endpoints
        captured before any boundary points are spliced in, so it never reads a
        mutated `seg.p0` / `seg.p0_seg1` (the source of the multi-crossing bridge
        corruption).
        """
        a = p0.get_point()
        b = p1.get_point()
        ab = b - a
        denom = float(ab @ ab)
        t = 0.0 if denom == 0.0 else float((np.asarray(point) - a) @ ab) / denom
        t = min(1.0, max(0.0, t))
        ca = p0.get_cdist(stability)
        cb = p1.get_cdist(stability)
        return (1 - t) * ca + t * cb
