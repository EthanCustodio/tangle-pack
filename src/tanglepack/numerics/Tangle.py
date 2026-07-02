import numpy as np
from rtree import index
from rtree.core import RTreeError
from collections import defaultdict
from itertools import count
from dataclasses import dataclass
from typing import Literal, Optional

from .Intersection import Intersection, ManifoldKey
from .BaseManifold import BaseManifold
from .Point import Point
from .BranchPoint import BranchPoint
from .Bridge import Bridge

import logging

logger = logging.getLogger(__name__)


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
    """_summary_

    Attributes:
            _rtree:         spatial rtree to store all segments based
                            on their bounding boxes
            _seg_lookup:    Dictionary keyed by unique segment ids that
                            stores segments
            _manifold_segs: Dictionary keyed by manifolds which stores
                            a set of all segment ids for that manifold

    Raises:
        ValueError: _description_

    Returns:
        _type_: _description_

    Yields:
        _type_: _description_
    """

    _ids = count(0)  # global generator -> every segment is given a unique int key

    def __init__(self):
        """

        Initializes the Tangle object with an R-tree for spatial indexing
        of segments, a lookup dictionary for segments, and a mapping
        of manifolds to their segments.

        Attributes:
            _rtree: spatial rtree to store all segments based on their bounding boxes
            _edge_seen:
            _seg_lookup: Dictionary keyed by unique segment ids that
                stores segments
            _manifold_segs Dictionary keyed by manifolds which stores
                a set of all segments

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
        self._edge_seen: set[frozenset[int]] = set()
        self._seg_lookup: dict[int, _Segment] = {}
        self._manifold_segs: defaultdict[BaseManifold, set[int]] = defaultdict(set)

        self._intersecting_segments: set[frozenset[int]] = set()
        # keyed by the crossing's segment pair, so one segment can host many crossings
        self._intersecting_coords: dict[frozenset[int], tuple[float, float]] = {}
        self._intersecting_points: dict[frozenset[int], BranchPoint] = {}

        self._intersections: list[Intersection] = []
        self._intersection_by_seg: defaultdict[int, list[Intersection]] = defaultdict(
            list
        )
        self._processed_pairs: set[frozenset[int]] = set()

        # self.bridges = None

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
        self._edge_seen.clear()
        self._seg_lookup.clear()
        self._manifold_segs.clear()
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

            if seg_id_pair in self._processed_pairs:
                continue  # already processed this pair

            seg1_id, seg2_id = tuple(seg_id_pair)

            if seg1_id not in self._seg_lookup or seg2_id not in self._seg_lookup:
                # drop stale pair and skip
                self._intersecting_segments.discard(seg_id_pair)
                continue

            seg_1, seg_2 = self._seg_lookup[seg1_id], self._seg_lookup[seg2_id]

            # A real crossing is always one unstable + one stable segment.
            if seg_1.manifold.stability == seg_2.manifold.stability:
                self._discard_same_stability(seg_id_pair, seg_1, seg_2)
                continue

            point = self._find_true_intersection(seg_1, seg_2)

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

    def _find_true_intersection(self, seg1: _Segment, seg2: _Segment):
        """
        Takes in two segments that are intersecting and returns
        the true point of intersection between them
        """
        p0_seg1, p1_seg1 = seg1.p0.get_point(), seg1.p0_seg1.get_point()
        p0_seg2, p1_seg2 = seg2.p0.get_point(), seg2.p0_seg1.get_point()

        a, b, c, d = map(np.asarray, (p0_seg1, p1_seg1, p0_seg2, p1_seg2))
        # 2Ã—2 matrix
        M = np.vstack((b - a, c - d)).T  # [[Bx-Ax, Cx-Dx], [By-Ay, Cy-Dy]]
        rhs = c - a

        det = np.linalg.det(M)
        if abs(det) < 1e-16:  # nearly parallel; shouldn't happen here
            raise ValueError("Segments are parallel or degenerate")

        t, s = np.linalg.solve(M, rhs)  # t on AB, s on CD
        # (Optional sanity: assert 0<=t<=1 and 0<=s<=1 if you didn't already know)
        return a + t * (b - a)

    @staticmethod
    def _same_unstable_branch(seg_1: _Segment, seg_2: _Segment) -> bool:
        """
        Whether two segments lie on the same manifold branch -- one physical curve.

        A bridge and its forward images are all sections of the single unstable
        manifold of a given (fixed point, orbit index, branch), even though they are
        held as separate objects. Two segments of one curve can never cross (the
        manifold is simple), so the Tangle must treat the whole branch as "self" --
        the bare object-identity check is too narrow once a branch is split across
        bridge/image objects. Falls back to object identity when a key is missing.
        """
        k1 = getattr(seg_1.manifold, "manifold_key", None)
        k2 = getattr(seg_2.manifold, "manifold_key", None)
        if k1 is None or k2 is None:
            return seg_1.manifold is seg_2.manifold
        # key = (fixed_point, stability, orbit_index, branch_index)
        return k1[0] is k2[0] and k1[1] == k2[1] and k1[2] == k2[2] and k1[3] == k2[3]

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

    def add_manifold(self, manifold: BaseManifold):
        """
        Indexes every segment of a manifold in the rtree

        Note:
            Does not check if segments have already been inserted
        """
        # ---------- A. purge old entries (if any) ----------
        if manifold in self._manifold_segs:
            # copy â†’ because _remove_segment mutates the same set
            for sid in list(self._manifold_segs[manifold]):
                self._remove_segment(sid)

        # _insert_segment records each new id in _manifold_segs itself; collecting
        # its return values here would sweep the None it returns for already-seen
        # edges into the id set and poison later lookups.
        for segment in self._segments_of(manifold):
            self._insert_segment(segment)

    def intersections_for_segment(self, seg: _Segment):
        """
        Yield true intersections with segments already stores in the rtree.
        """
        for candidate_id in self._rtree.intersection(seg.bounds):

            other_segment = self._seg_lookup[candidate_id]
            if other_segment.manifold is seg.manifold:
                continue  # same manifold â€“ skip or apply rule

            if seg.intersects(other_segment):
                yield other_segment

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
            if seg_id_pair in self._processed_pairs:
                continue
            # only process pairs where at least one segment belongs to manifold
            if not (seg_id_pair & manifold_seg_ids):
                continue

            seg1_id, seg2_id = tuple(seg_id_pair)  # frozenset unpack

            if seg1_id not in self._seg_lookup or seg2_id not in self._seg_lookup:
                self._intersecting_segments.discard(seg_id_pair)
                continue

            seg_1 = self._seg_lookup[seg1_id]
            seg_2 = self._seg_lookup[seg2_id]

            if seg_1.manifold.stability == seg_2.manifold.stability:
                self._discard_same_stability(seg_id_pair, seg_1, seg_2)
                continue

            point = self._find_true_intersection(seg_1, seg_2)

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
            self._intersecting_coords[seg_id_pair] = tuple(point)
            self._intersecting_points[seg_id_pair] = branch_point

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
            new_intersections.append(intersection)

            self._processed_pairs.add(seg_id_pair)

        return new_intersections

    # ------------- internal helpers -----------------
    def _insert_segment(self, seg: _Segment) -> Optional[int]:
        """
        Inserts a segment into the rtree and local dictionary

        Parameters:
            seg (_Segment): segment to be inserted

        Returns:
            The segment id, or None if this edge was already indexed.
        """
        # edge key defined by the id of two points
        edge_key = frozenset((id(seg.p0), id(seg.p0_seg1)))
        if edge_key in self._edge_seen:
            return None  # already indexed -> do NOT duplicate

        # choose a new id there for a new segment
        sid = next(Tangle._ids) if seg.id is None else seg.id

        seg.id = sid  # assign the id to the segment

        # insert segment into the rtree and dictionary
        self._rtree.insert(sid, seg.bounds)
        self._edge_seen.add(edge_key)
        self._seg_lookup[sid] = seg
        self._manifold_segs[seg.manifold].add(sid)

        for cand_id in self._rtree.intersection(seg.bounds):

            if cand_id == sid:
                continue

            other = self._seg_lookup[cand_id]
            if other.manifold is seg.manifold:
                continue

            if seg.intersects(other):
                self._intersecting_segments.add(frozenset((sid, cand_id)))

        return sid

    def _remove_segment(self, sid: int):
        """
        Remove all references to the segment

        Parameters:
            sid (int): segment id
        """
        # Check if segment exists before trying to remove it
        if sid not in self._seg_lookup:
            return

        seg = self._seg_lookup.pop(sid)

        # Safely remove from rtree
        try:
            self._rtree.delete(sid, seg.bounds)
        except RTreeError:
            # Segment may have already been deleted from the rtree; anything
            # else deserves to surface rather than be swallowed.
            logger.debug("rtree delete failed for segment %s", sid, exc_info=True)

        # Safely remove from manifold's segment set
        if seg.manifold in self._manifold_segs:
            self._manifold_segs[seg.manifold].discard(sid)

        self._edge_seen.discard(frozenset((id(seg.p0), id(seg.p0_seg1))))

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

    def _orientation(a, b, c, eps=1e-15):
        """
        Returns:
        0 if a, b, c are collinear
        1 if clockwise
        2 if counterclockwise
        """
        val = (b[1] - a[1]) * (c[0] - b[0]) - (b[0] - a[0]) * (c[1] - b[1])
        if abs(val) < eps:
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

    def cut_manifold(self):
        pass

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

    def _boundary_point(
        self,
        p0: Point,
        p1: Point,
        stability: str,
        intersection_coords: tuple[float, float],
        side: Literal["root", "tail"],
    ) -> tuple[Point, float]:
        """
        Create a bridge boundary Point just outside the segment [p0, p1], offset
        10% from the crossing toward the endpoint on `side`, and return it together
        with its fractional position `t` along [p0, p1].

        The segment endpoints are passed in explicitly (not read from a live
        `_Segment`) so the offset and cdist are always computed against the
        ORIGINAL endpoints captured before any splicing. The returned `t` lets the
        caller splice several boundary points into one segment in correct
        geometric order.

        Args:
            p0: Lower-cdist endpoint of the original segment.
            p1: Higher-cdist endpoint of the original segment.
            stability: Manifold stability, for the cdist interpolation.
            intersection_coords: True (x, y) of the crossing.
            side: "root" offsets toward p0 (lower cdist); "tail" toward p1.

        Returns:
            (boundary_point, t) where t is the clamped fractional position of the
            boundary point along [p0, p1].
        """
        if side == "root":
            seg_point = p0
        elif side == "tail":
            seg_point = p1
        else:
            raise ValueError(f"Invalid side: {side}")

        new_point = self._linear_interpolation(
            intersection_coords, seg_point.get_point(), 0.1
        )

        # cdist evaluated at the boundary point's ACTUAL location (the 10% offset
        # point), not at the crossing: a Point's cdist must reflect where the point
        # physically sits, otherwise the error is baked in and scaled up every time
        # the bridge is iterated. Evaluated against the original endpoints p0, p1 so
        # it cannot be corrupted by earlier insertions.
        new_cdist = self._cdist_between(p0, p1, stability, np.asarray(new_point))

        boundary = Point(new_point[0], new_point[1], new_cdist)
        t = self._fractional_position(p0, p1, np.asarray(new_point))
        return boundary, t

    def _linear_interpolation(self, p0, p1, alpha):
        """
        Does a linear interpolation between two points p0 and p1 to get a
        point that is alpha away from the true intersection point.

        Args:
            p0 (Point): First Point
            p1 (Point): Second Point
            alpha (Float): percent distance to interpolate

        Returns:
            tuple(float, float): resulting point coordinates
        """
        x = (1 - alpha) * p0[0] + alpha * p1[0]
        y = (1 - alpha) * p0[1] + alpha * p1[1]
        return (x, y)

    def _cache_boundary_preiterate(
        self,
        boundary_point: Point,
        p0: Point,
        p1: Point,
        stability: str,
    ) -> None:
        """
        Approximate and cache `prev_iterate` (unstable) or `next_iterate` (stable)
        on a freshly created bridge boundary point.

        Uses the same weighted interpolation as `_linear_interpolation`: the boundary
        point sits alpha=0.1 away from the intersection toward seg_point, so the
        preiterate is interpolated with the same weight. The endpoints p0, p1 are
        passed explicitly (the ORIGINAL segment endpoints) so the cached preiterate
        is never read off a mutated `_Segment`.
        """
        if stability == "unstable":
            pre0 = p0.prev_iterate
            pre1 = p1.prev_iterate
        else:
            pre0 = p0.next_iterate
            pre1 = p1.next_iterate

        if pre0 is None or pre1 is None:
            return  # can't interpolate; leave None and rely on guard in refine

        # 0.9 * intersection_preiterate + 0.1 * seg_point_preiterate
        # (mirrors the 0.1 alpha used in _linear_interpolation for the point itself)
        coords = (
            0.9 * 0.5 * (pre0.get_point() + pre1.get_point()) + 0.1 * pre0.get_point()
        )

        cached = Point(coords[0], coords[1])
        if stability == "unstable":
            boundary_point.prev_iterate = cached
        else:
            boundary_point.next_iterate = cached

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

    @staticmethod
    def _fractional_position(p0: Point, p1: Point, point: np.ndarray) -> float:
        """Fractional projection t of `point` onto the segment [p0, p1], clamped to [0, 1]."""
        a = p0.get_point()
        b = p1.get_point()
        ab = b - a
        denom = float(ab @ ab)
        t = 0.0 if denom == 0.0 else float((np.asarray(point) - a) @ ab) / denom
        return min(1.0, max(0.0, t))
