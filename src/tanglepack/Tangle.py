import numpy as np
from rtree import index
from collections import defaultdict
from itertools import count
from dataclasses import dataclass

from .BaseManifold import BaseManifold
from .Point import Point
from .BranchPoint import BranchPoint


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

    _ids = count(0)  # global generator -> every segment is given a unique int key

    def __init__(self):
        """

        Parameters:
            _rtree:         spatial rtree to store all segments based
                            on their bounding boxes
            _seg_lookup:    Dictionary keyed by unique segment ids that
                            stores segments
            _manifold_segs: Dictionary keyed by manifolds which stores
                            a set of all segment ids for that manifold
        """

        p = index.Property()
        p.dimension = 2
        self._rtree = index.Index(properties=p)  # R-tree
        self._edge_seen: set[frozenset[int]] = set()
        self._seg_lookup: dict[int, _Segment] = {}
        self._manifold_segs: defaultdict[BaseManifold, set[int]] = defaultdict(set)

        self._intersecting_segments: set[frozenset[int]] = set()
        self._intersecting_coords: dict[int, tuple[float, float]] = {}
        self._intersecting_points: dict[int, BranchPoint] = {}

    def populate_intersection_dict(self):
        """
        Takes all intersection pairs in _intersection_segments and
        finds the true intersections and adds those to a list
        """

        for seg_id_pair in self._intersecting_segments:

            seg1_id, seg2_id = seg_id_pair
            seg_1, seg_2 = self._seg_lookup[seg1_id], self._seg_lookup[seg2_id]

            point = self._find_true_intersection(seg_1, seg_2)

            # Figure out how to compute the canonical distance
            # must compute differently for stable or unstable
            seg_1_cdist = 0.5 * (
                seg_1.p0.get_cdist(seg_1.manifold.stability)
                + seg_1.p0_seg1.get_cdist(seg_1.manifold.stability)
            )

            seg_2_cdist = 0.5 * (
                seg_2.p0.get_cdist(seg_2.manifold.stability)
                + seg_2.p0_seg1.get_cdist(seg_2.manifold.stability)
            )

            unstable_cdist = seg_1 if seg_1.manifold.stability == "unstable" else seg_2
            stable_cdist = seg_1 if seg_1.manifold.stability == "unstable" else seg_2

            branch_point = BranchPoint(2, point[0], point[1])

            self._intersecting_coords[seg1_id] = point
            self._intersecting_coords[seg2_id] = point
            self._intersecting_points[seg1_id] = branch_point
            self._intersecting_points[seg2_id] = branch_point

    def _find_true_intersection(self, seg1: _Segment, seg2: _Segment):
        """
        Takes in two segments that are intersecting and returns
        the true point of intersection between them
        """
        p0_seg1, p1_seg1 = seg1.p0.get_point(), seg1.p0_seg1.get_point()
        p0_seg2, p1_seg2 = seg2.p0.get_point(), seg2.p0_seg1.get_point()

        a, b, c, d = map(np.asarray, (p0_seg1, p1_seg1, p0_seg2, p1_seg2))
        # 2×2 matrix
        M = np.vstack((b - a, c - d)).T  # [[Bx-Ax, Cx-Dx], [By-Ay, Cy-Dy]]
        rhs = c - a

        det = np.linalg.det(M)
        if abs(det) < 1e-12:  # nearly parallel; shouldn't happen here
            raise ValueError("Segments are parallel or degenerate")

        t, s = np.linalg.solve(M, rhs)  # t on AB, s on CD
        # (Optional sanity: assert 0<=t<=1 and 0<=s<=1 if you didn't already know)
        return a + t * (b - a)

    def add_manifold(self, manifold: BaseManifold):
        """
        Indexes every segment of a manifold in the rtree

        Note:
            Does not check if segments have already been inserted
        """
        # ---------- A. purge old entries (if any) ----------
        if manifold in self._manifold_segs:
            # copy → because _remove_segment mutates the same set
            for sid in list(self._manifold_segs[manifold]):
                self._remove_segment(sid)

        seg_ids = {self._insert_segment(s) for s in self._segments_of(manifold)}
        self._manifold_segs[manifold].update(seg_ids)

    def update_manifold(self, manifold: BaseManifold, changed_points: set[Point]):
        """
        Re‑index only segments that touch the nodes in `changed_nodes`

        Note:
            new points from `ManifoldMachine.refine_manifold`.
        """
        # get which segments in manifold share a point with changed_points
        # and delete them from memory before we add them back
        dirty = self._segments_touching(manifold, changed_points)
        for sid in dirty:  # remove old versions
            self._remove_segment(sid)

        # Collects all the new segments which need to be added
        rebuilt_segments = set()
        for point in changed_points:
            prev_point = manifold.walk_back(None, point)
            next_point = manifold.walk_fwd(None, point)

            if prev_point is not None:
                rebuilt_segments.add(_Segment(manifold, prev_point, point))
            if next_point is not None:
                rebuilt_segments.add(_Segment(manifold, point, next_point))

        # new_ids = {self._insert_segment(s) for s in self._segments_of_nodes(dirty)}
        new_ids = {self._insert_segment(s) for s in rebuilt_segments}
        self._manifold_segs[manifold].update(new_ids)

    def intersections_for_segment(self, seg: _Segment):
        """
        Yield true intersections with segments already stores in the rtree.
        """
        for candidate_id in self._rtree.intersection(seg.bounds):

            other_segment = self._seg_lookup[candidate_id]
            if other_segment.manifold is seg.manifold:
                continue  # same manifold – skip or apply rule

            if seg.intersects(other_segment):
                yield other_segment

    # ------------- internal helpers -----------------
    def _insert_segment(self, seg: _Segment) -> int:
        """
        Inserts a segment into the rtree and local dictionary

        Parameters:
            seg (_Segment): segment to be inserted
        """
        edge_key = frozenset((id(seg.p0), id(seg.p0_seg1)))
        if edge_key in self._edge_seen:
            return None  # already indexed → do NOT duplicate

        # choose a new id there for a new segment
        sid = next(Tangle._ids) if seg.id is None else seg.id

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
        seg = self._seg_lookup.pop(sid)
        self._rtree.delete(sid, seg.bounds)
        self._manifold_segs[seg.manifold].remove(sid)
        self._edge_seen.discard(frozenset((id(seg.p0), id(seg.p0_seg1))))

    def _segments_of(self, manifold: BaseManifold):
        """
        Generator which returns the segments in a manifold
        """
        prev_point = None
        curr_point = manifold.root

        while curr_point is not None:

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

    # def _segments_of_nodes(self, old_ids):
    #     # rebuild the exact segments we just deleted
    #     for sid in old_ids:
    #         segment = self._seg_lookup.get(sid)
    #         if segment is None:  # already removed
    #             continue
    #         yield from (_Segment(None, segment.manifold, segment.p0, segment.p0_seg1))

    # def _promote_to_branchpoint(p_int, seg1, seg2):
    #     """Return BranchPoint + surrounding pointer surgery."""
    #     bp = BranchPoint(num_branches=2, x=p_int[0], y=p_int[1])
    #     # 1️⃣ split seg1 into (p0‑bp) and (bp‑p0_seg1)
    #     _splice_segment(seg1, bp)
    #     # 2️⃣ split seg2 similar
    #     _splice_segment(seg2, bp)
    #     # 3️⃣ create higher‑level Intersection wrapper if desired
    #     return bp

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
        (p0_seg1, p1_seg1) = segA
        (p0_seg2, p1_seg2) = segB

        o1 = Tangle._orientation(p0_seg1, p1_seg1, p0_seg2)
        o2 = Tangle._orientation(p0_seg1, p1_seg1, p1_seg2)
        o3 = Tangle._orientation(p0_seg2, p1_seg2, p0_seg1)
        o4 = Tangle._orientation(p0_seg2, p1_seg2, p1_seg1)

        # General case:
        if o1 != o2 and o3 != o4:
            return True

        return False  # no intersection if collinear or not straddling
