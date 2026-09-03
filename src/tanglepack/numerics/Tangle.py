from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from rtree import index
from rtree.core import RTreeError
from collections import defaultdict
from itertools import count
from dataclasses import dataclass
from typing import Optional

from .Intersection import Intersection, ManifoldKey
from .BaseManifold import BaseManifold
from .Point import Point
from .BranchPoint import BranchPoint
from .Bridge import Bridge

import logging

logger = logging.getLogger(__name__)

# Relative threshold shared by the two geometric predicates below. Both compare an
# area-like quantity (a cross product, a 2x2 determinant) that scales as the product
# of the two lengths involved, so dividing that product out makes the threshold a
# bound on sin(angle) between the two directions and the predicates scale invariant.
_EPS_REL = 1e-12


"""
Dev Notes:

WARNING: this code must identify fixed points as intersection points. Consider 
    the complications of that 
"""


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
    def bounds(self):
        """
        Returns the boundaries of the bounding box
        made from self.p0 and self.p0_seg1
        """

        x0, y0 = self.p0.get_point()
        x1, y1 = self.p0_seg1.get_point()

        return (min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))

    def intersects(self, other) -> bool:
        p0 = tuple(self.p0.get_point())  # (x, y)
        p0_seg1 = tuple(self.p0_seg1.get_point())
        q0 = tuple(other.p0.get_point())
        q1 = tuple(other.p0_seg1.get_point())

        return Tangle._do_segments_intersect((p0, p0_seg1), (q0, q1))


class Tangle:
    """
    Spatial index of manifold segments that detects and resolves crossings.

    Every adjacent pair of points on an indexed manifold is stored as a
    _Segment in an rtree. Candidate crossings are found by bounding-box
    queries, resolved to true intersection points, and cut into bridges.

    Attributes:
        _rtree: spatial rtree storing all segments by their bounding boxes.
        _seg_lookup: Dictionary keyed by unique segment ids that stores
            segments.
        _manifold_segs: Dictionary keyed by manifolds which stores the set of
            all segment ids for that manifold.
        _intersecting_segments: Set of frozenset pairs of segment ids that
            cross.
        _intersections: List of resolved Intersection objects.
    """

    _ids = count(0)  # global generator -> every segment is given a unique int key

    def __init__(self):
        """

        Initializes the Tangle object with an R-tree for spatial indexing
        of segments, a lookup dictionary for segments, and a mapping
        of manifolds to their segments.

        Attributes:
            _rtree: spatial rtree to store all segments based on their bounding boxes
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
            _intersecting_coords: dictionary mapping segment ids to their
                intersection coordinates
            _intersecting_points: dictionary mapping segment ids to their
                corresponding BranchPoint objects
        """

        p = index.Property()
        p.dimension = 2
        self._rtree = index.Index(properties=p)  # R-tree
        self._edge_to_sid: dict[frozenset[int], int] = {}
        self._seg_lookup: dict[int, _Segment] = {}
        self._manifold_segs: defaultdict[BaseManifold, set[int]] = defaultdict(set)
        self._seg_manifolds: defaultdict[int, set[BaseManifold]] = defaultdict(set)

        self._intersecting_segments: set[frozenset[int]] = set()
        # keyed by the crossing's segment pair, so one segment can host many crossings
        self._intersecting_coords: dict[frozenset[int], tuple[float, float]] = {}
        self._intersecting_points: dict[frozenset[int], BranchPoint] = {}

        self._intersections: list[Intersection] = []
        self._intersection_by_seg: defaultdict[int, list[Intersection]] = defaultdict(
            list
        )
        self._processed_pairs: set[frozenset[int]] = set()

    @staticmethod
    def _key_of(seg: _Segment) -> Optional[ManifoldKey]:
        """Read the ManifoldKey stored on the segment's manifold, if set."""
        return getattr(seg.manifold, "manifold_key", None)

    def clear_all(self):
        """
        Completely clear all Tangle state, useful when recomputing everything from scratch.
        """
        # Reinitialize the rtree
        p = index.Property()
        p.dimension = 2
        self._rtree = index.Index(properties=p)

        # Clear all dictionaries and sets
        self._edge_to_sid.clear()
        self._seg_lookup.clear()
        self._manifold_segs.clear()
        self._seg_manifolds.clear()
        self._intersecting_segments.clear()
        self._intersecting_coords.clear()
        self._intersecting_points.clear()

        self._intersections.clear()
        self._intersection_by_seg.clear()
        self._processed_pairs.clear()

    def populate_intersection_dict(self):
        """
        Takes all intersection pairs in _intersection_segments and
        finds the true intersections and adds those to a list

        Resolve every detected crossing pair into an Intersection.

        Each crossing is keyed by its segment PAIR (frozenset of the two segment
        ids), so a single segment may participate in many crossings without any
        being overwritten. Only unstable x stable pairs are kept. A same-stability
        pair (u x u or s x s) is geometrically impossible (see CLAUDE.md's
        fundamental invariant) -- if one appears it is a polygonal/numerical
        artifact of two near-parallel manifolds straddling, so it is logged and
        discarded, never turned into an Intersection.
        """
        for seg_id_pair in list(self._intersecting_segments):
            self._resolve_crossing_pair(seg_id_pair)

    def _resolve_crossing_pair(self, seg_id_pair: frozenset[int]):
        """
        Resolve one detected crossing pair into an Intersection.

        Handles all the shared bookkeeping: the true crossing point, the cdists
        interpolated at that point (not the segment midpoint), the BranchPoint,
        and the intersection/coord/point registries.

        Returns:
            The new Intersection, or None when the pair was already processed,
            references a removed segment (dropped as stale), is a
            same-stability near-tangency (discarded per the fundamental
            invariant -- see populate_intersection_dict), or is near-parallel
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

        # A real crossing is always one unstable + one stable segment.
        if seg_1.manifold.stability == seg_2.manifold.stability:
            self._discard_same_stability(seg_id_pair, seg_1, seg_2)
            return None

        point = self._find_true_intersection(seg_1, seg_2)
        if point is None:
            # near-parallel pair (already logged): discard it rather than let it
            # abort the whole intersection pass
            self._intersecting_segments.discard(seg_id_pair)
            self._processed_pairs.add(seg_id_pair)
            return None

        # cdist interpolated at the true crossing (not the segment midpoint)
        seg_1_cdist = self._cdist_at_point(seg_1, point)
        seg_2_cdist = self._cdist_at_point(seg_2, point)

        unstable_cdist = (
            seg_1_cdist if seg_1.manifold.stability == "unstable" else seg_2_cdist
        )
        stable_cdist = (
            seg_1_cdist if seg_1.manifold.stability == "stable" else seg_2_cdist
        )

        branch_point = BranchPoint(
            2, (unstable_cdist, stable_cdist), point[0], point[1]
        )

        if seg_1.manifold.stability == "unstable":
            manifold_a_key, cdist_a = Tangle._key_of(seg_1), seg_1_cdist
            manifold_b_key, cdist_b = Tangle._key_of(seg_2), seg_2_cdist
        else:
            manifold_a_key, cdist_a = Tangle._key_of(seg_2), seg_2_cdist
            manifold_b_key, cdist_b = Tangle._key_of(seg_1), seg_1_cdist

        intersection = Intersection.from_segments(
            coords=tuple(point),
            unstable_cdist=cdist_a,
            stable_cdist=cdist_b,
            seg1_id=seg1_id,
            seg2_id=seg2_id,
            manifold_a_key=manifold_a_key,
            manifold_b_key=manifold_b_key,
        )

        self._intersections.append(intersection)
        self._intersection_by_seg[seg1_id].append(intersection)
        self._intersection_by_seg[seg2_id].append(intersection)

        self._intersecting_coords[seg_id_pair] = tuple(point)
        self._intersecting_points[seg_id_pair] = branch_point

        self._processed_pairs.add(seg_id_pair)

        return intersection

    def _find_true_intersection(
        self, seg1: _Segment, seg2: _Segment
    ) -> Optional[NDArray[np.float64]]:
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
            The crossing point, or None if the two directions are parallel to
            within ``_EPS_REL`` (the caller discards such a pair).
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
        # (Optional sanity: assert 0<=t<=1 and 0<=s<=1 if you didn't already know)
        return a + t * (b - a)

    def _discard_same_stability(
        self, seg_id_pair: "frozenset[int]", seg_1: _Segment, seg_2: _Segment
    ) -> None:
        """
        Drop a same-stability segment pair -- it can never be a real crossing.

        By the fundamental invariant two unstable manifolds (or two stable manifolds)
        never cross, whether they belong to the same fixed point or different ones. So
        every same-stability pair the spatial index turns up is a near-tangency: one
        curve folding near itself, two orbit branches sharing a resonance island, or
        two fixed points' manifolds running close. None is a malfunction to surface --
        it is expected geometry for a folded tangle -- so all are dropped (logged at
        debug, not warning). A genuine *transversal* self-crossing (a real numerics
        bug) is caught by the invariant test ``test_no_same_stability_crossing``, not
        by flagging every near-tangency here.
        """
        logger.debug(
            "Dropping same-stability (%s) segment pair %s as a near-tangency",
            seg_1.manifold.stability,
            tuple(seg_id_pair),
        )
        self._processed_pairs.add(seg_id_pair)

    def add_manifold(
        self,
        manifold: BaseManifold,
        index_segments: bool = True,
        detect_crossings: bool = True,
    ) -> None:
        """
        Registers every segment of a manifold and detects its crossings.

        Args:
            manifold (BaseManifold): The manifold to register.
            index_segments (bool): If True (default) each segment is also
                inserted into the rtree so later additions can find it. Pass
                False for an iterated bridge: a bridge is unstable manifold,
                so by the fundamental invariant the only curves it can
                legitimately cross are stable ones, and those are already in
                the index. Its crossings are found by querying; inserting its
                segments too would only let later bridges detect
                bridge-x-bridge (unstable x unstable) pairs, which are
                numerical artifacts that get discarded anyway. Skipping the
                inserts is what keeps repeated bridge iteration fast.
            detect_crossings (bool): If False, genuinely new segments are
                registered without being queried against the rtree. Used when
                cutting bridges out of an already-indexed manifold: every
                crossing on that stretch of curve was already found on the
                parent's segments, so querying again would only duplicate it.

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
        for segment in self._segments_of(manifold):
            self._insert_segment(
                segment,
                index_in_rtree=index_segments,
                detect_crossings=detect_crossings,
            )

    def add_manifolds(self, manifolds: list[BaseManifold]):
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

        properties = index.Property()
        properties.dimension = 2
        self._rtree = index.Index(
            ((seg.id, seg.bounds, None) for seg in new_segments),
            properties=properties,
        )

        # Every segment is queried, so each crossing pair comes up twice (once
        # from each side); keeping only the lower-id side does each exact
        # intersection test once.
        for segment in new_segments:
            for candidate_id in self._rtree.intersection(segment.bounds):
                if candidate_id >= segment.id:
                    continue

                other = self._seg_lookup[candidate_id]
                if other.manifold is segment.manifold:
                    continue

                if segment.intersects(other):
                    self._intersecting_segments.add(
                        frozenset((segment.id, candidate_id))
                    )

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
            List of newly created Intersection objects.
        """
        manifold_seg_ids = self._manifold_segs.get(manifold, set())
        new_intersections: list[Intersection] = []

        for seg_id_pair in list(self._intersecting_segments):
            # only process pairs where at least one segment belongs to manifold
            if not (seg_id_pair & manifold_seg_ids):
                continue

            intersection = self._resolve_crossing_pair(seg_id_pair)
            if intersection is not None:
                new_intersections.append(intersection)

        return self._collapse_noise_crossings(new_intersections)

    def _collapse_noise_crossings(
        self, intersections: list[Intersection], rtol: float = 1e-9
    ) -> list[Intersection]:
        """
        Collapse multiple detections of one transversal crossing into one.

        A forward-mapped (iterated-bridge) polyline carries coordinate noise
        amplified by the map; near a crossing it can zig-zag across the stable
        curve, so one physical crossing is detected several times within noise
        of a single point of the bridge — unstable cdists agreeing to ~1e-12
        RELATIVE, versus >= 1e-3 between genuine neighbouring crossings. Noise
        flips add detections in PAIRS while a genuine transversal crossing
        adds one, so within each noise-run (consecutive unstable cdists closer
        than ``rtol`` relative): an odd-length run keeps its median detection,
        an even-length run is a grazing artifact with no net crossing and is
        dropped entirely. Dropped crossings are purged from every tangle
        structure so bridge cutting never sees them.

        Args:
            intersections: Freshly resolved crossings of one query manifold.
            rtol: Relative unstable-cdist gap under which two detections are
                the same point of the manifold.

        Returns:
            The surviving crossings.
        """
        if len(intersections) < 2:
            return intersections

        ordered = sorted(intersections, key=lambda ix: ix.unstable_cdist)
        runs: list[list[Intersection]] = [[ordered[0]]]
        for previous, current in zip(ordered, ordered[1:]):
            gap = current.unstable_cdist - previous.unstable_cdist
            if gap <= rtol * max(abs(current.unstable_cdist), 1.0):
                runs[-1].append(current)
            else:
                runs.append([current])

        kept: list[Intersection] = []
        for run in runs:
            if len(run) == 1:
                kept.append(run[0])
                continue
            if len(run) % 2 == 1:
                keep = run[len(run) // 2]
                kept.append(keep)
            else:
                keep = None
            logger.info(
                "Collapsed %d noise detections of one crossing at unstable "
                "cdist %.9g (%s)",
                len(run),
                run[0].unstable_cdist,
                "kept median" if keep is not None else "grazing pair dropped",
            )
            for ix in run:
                if ix is not keep:
                    self._purge_crossing(ix)
        return kept

    def _purge_crossing(self, intersection: Intersection) -> None:
        """Remove one resolved crossing from every tangle structure."""
        pair = intersection.seg_ids
        if intersection in self._intersections:
            self._intersections.remove(intersection)
        if pair:
            for seg_id in pair:
                by_seg = self._intersection_by_seg.get(seg_id)
                if by_seg and intersection in by_seg:
                    by_seg.remove(intersection)
            self._intersecting_segments.discard(pair)
            self._intersecting_coords.pop(pair, None)
            self._intersecting_points.pop(pair, None)

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
            self._rtree.insert(sid, seg.bounds)
        seg.in_rtree = index_in_rtree
        self._edge_to_sid[edge_key] = sid
        self._seg_lookup[sid] = seg
        self._claim_segment(sid, seg.manifold)

        if not detect_crossings:
            return sid

        for cand_id in self._rtree.intersection(seg.bounds):

            if cand_id == sid:
                continue

            other = self._seg_lookup[cand_id]
            if other.manifold is seg.manifold:
                continue

            if seg.intersects(other):
                self._intersecting_segments.add(frozenset((sid, cand_id)))

        return sid

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
                self._rtree.delete(sid, seg.bounds)
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

    def _segments_touching(self, manifold, nodes):
        """
        Returns segment ids for all segments which contain a node
        """
        return {
            sid
            for sid in self._manifold_segs[manifold]
            if self._seg_lookup[sid].p0 in nodes
            or self._seg_lookup[sid].p0_seg1 in nodes
        }

    @staticmethod
    def _orientation(a, b, c, eps_rel: float = _EPS_REL) -> int:
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
    def _do_segments_intersect(segA, segB):
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

    def create_bridges(self, for_manifold=None, fixed_point=None):
        """
        Cut every indexed unstable manifold into bridges at its crossings with the
        stable manifold(s).

        Crossings are grouped by their parent unstable manifold and sorted by the
        true crossing cdist, so bridges never span two manifolds and a single
        segment that hosts several crossings is handled correctly (each crossing is
        a distinct cut, not deduped by segment id).

        Args:
            for_manifold: If given, only build bridges for crossings that involve a
                segment of this specific manifold.
            fixed_point: If given, only build bridges on unstable manifolds that
                emanate from this fixed point. Each Bridge still records its own
                fixed_point, so a global call (fixed_point=None) followed by
                filtering on bridge.fixed_point is equivalent to per-fixed-point
                calls.

        Returns:
            List of Bridge objects, doubly linked via next_bridge / prev_bridge.
        """
        from collections import defaultdict

        for_manifold_segs = (
            self._manifold_segs.get(for_manifold, set())
            if for_manifold is not None
            else None
        )

        # --- 1. Collect crossings grouped by their parent unstable manifold ---
        # Each crossing captures the ORIGINAL segment endpoints (p0, p1) now,
        # before any boundary point is spliced in. Nothing downstream reads
        # seg.p0 / seg.p0_seg1 again, so a segment that hosts several crossings is
        # never corrupted by a previous crossing's insertion (the line-696 bug).
        # entry: (cdist, crossing_coords, orig_p0, orig_p1)
        manifold_crossings: dict[
            BaseManifold,
            list[tuple[float, tuple[float, float], Point, Point]],
        ] = defaultdict(list)

        for sid_pair in self._intersecting_segments:
            if for_manifold_segs is not None and not (sid_pair & for_manifold_segs):
                continue
            if sid_pair not in self._intersecting_coords:
                continue  # not a resolved unstable x stable crossing

            sid1, sid2 = tuple(sid_pair)
            seg_1, seg_2 = self._seg_lookup[sid1], self._seg_lookup[sid2]

            if seg_1.manifold.stability == "unstable":
                u_seg = seg_1
            elif seg_2.manifold.stability == "unstable":
                u_seg = seg_2
            else:
                continue  # no unstable segment (shouldn't happen post-filter)

            if (
                fixed_point is not None
                and u_seg.manifold.fixed_point is not fixed_point
            ):
                continue

            coords = self._intersecting_coords[sid_pair]
            p0, p1 = u_seg.p0, u_seg.p0_seg1  # ORIGINAL endpoints
            cdist = self._cdist_between(
                p0, p1, u_seg.manifold.stability, np.asarray(coords)
            )
            manifold_crossings[u_seg.manifold].append((cdist, coords, p0, p1))

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
        for manifold, crossings in manifold_crossings.items():
            crossings.sort(key=lambda c: c[0])

            for i in range(len(crossings) - 1):
                _, coords1, p0_a, p1_a = crossings[i]
                _, coords2, p0_b, p1_b = crossings[i + 1]

                head = p0_a  # real point just below crossing i
                tail = p1_b  # real point just above crossing i+1

                # Two consecutive crossings on the SAME segment have no real point
                # between them; refine the manifold (splice one true point on the
                # segment, between the two crossings) so head/tail bracket distinct
                # crossings rather than spanning both.
                if p0_a is p0_b and p1_a is p1_b:
                    mid = self._insert_crossing_separator(
                        p0_a, p1_a, coords1, coords2, manifold
                    )
                    tail = mid

                bridge = Bridge(
                    root=head,
                    stability=manifold.stability,
                    stretch_param=manifold.stretch_param,
                    fixed_point=manifold.fixed_point,
                    tail=tail,
                    branch_index=manifold.branch_index,
                )
                # A bridge is a segment of its parent unstable manifold, so it lives
                # on the same branch. Carrying the parent's manifold_key means every
                # intersection later detected on this bridge (or on its forward
                # image) records which unstable branch it belongs to.
                bridge.manifold_key = manifold.manifold_key
                all_bridges.append(bridge)

        # --- 3. Register each bridge's segments under the bridge itself ---
        # Done only after every cut, because _insert_crossing_separator may still
        # splice a point into a span an earlier bridge already covers. The segments
        # are the parent's own, so they are neither re-inserted into the rtree nor
        # re-queried: every crossing on them was found on the parent already.
        for bridge in all_bridges:
            self.add_manifold(bridge, index_segments=False, detect_crossings=False)

        # --- 4. Wire next_bridge / prev_bridge doubly-linked list ---
        for i in range(len(all_bridges) - 1):
            all_bridges[i].next_bridge = all_bridges[i + 1]
            all_bridges[i + 1].prev_bridge = all_bridges[i]

        return all_bridges

    def _insert_crossing_separator(
        self,
        p0: Point,
        p1: Point,
        coords1: tuple[float, float],
        coords2: tuple[float, float],
        manifold: BaseManifold,
    ) -> Point:
        """
        Splice one real point into the adjacent segment ``[p0, p1]`` between two
        crossings that share it, and return it.

        Two consecutive crossings on the same segment have no manifold point between
        them, so neighbouring bridges cannot bracket them with distinct real points.
        This refines the manifold by inserting a true point on the segment at the
        midpoint of the two crossings (its cdist interpolated against the original
        endpoints), restoring one-crossing-per-segment locally.

        Args:
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
        if isinstance(p0, BranchPoint):
            p0.insert_point_forward(separator, manifold.branch_index)
        else:
            p0.insert_point_forward(separator)
        return separator

    def iter_intersection_coords(self) -> list[tuple[float, float]]:
        """Return the (x, y) of every detected crossing, exactly one per crossing."""
        return list(self._intersecting_coords.values())

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
