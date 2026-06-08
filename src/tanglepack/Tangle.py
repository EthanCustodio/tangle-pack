import numpy as np
from rtree import index
from collections import defaultdict
from itertools import count
from dataclasses import dataclass
from typing import Literal, Optional

from .Intersection import Intersection, ManifoldKey
from .BaseManifold import BaseManifold
from .Point import Point
from .BranchPoint import BranchPoint
from .Bridge import Bridge

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
        self._intersecting_coords: dict[int, tuple[float, float]] = {}
        self._intersecting_points: dict[int, BranchPoint] = {}

        self._intersections: list[Intersection] = []
        self._intersection_by_seg: dict[int, Intersection] = {}
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
        """

        for seg_id_pair in list(self._intersecting_segments):

            if seg_id_pair in self._processed_pairs:
                continue  # already processed this pair

            seg1_id, seg2_id = seg_id_pair

            if seg1_id not in self._seg_lookup or seg2_id not in self._seg_lookup:
                # drop stale pair and skip
                self._intersecting_segments.discard(seg_id_pair)
                continue

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

            unstable_cdist = (
                seg_1_cdist if seg_1.manifold.stability == "unstable" else seg_2_cdist
            )
            stable_cdist = (
                seg_1_cdist if seg_1.manifold.stability == "stable" else seg_2_cdist
            )

            branch_point = BranchPoint(
                2, (unstable_cdist, stable_cdist), point[0], point[1]
            )

            if (
                seg_1.manifold.stability == "unstable"
                or seg_2.manifold.stability != "unstable"
            ):
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
            self._intersection_by_seg[seg1_id] = intersection
            self._intersection_by_seg[seg2_id] = intersection

            self._intersecting_coords[seg1_id] = point
            self._intersecting_coords[seg2_id] = point
            self._intersecting_points[seg1_id] = branch_point
            self._intersecting_points[seg2_id] = branch_point

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
            # copy â†’ because _remove_segment mutates the same set
            for sid in list(self._manifold_segs[manifold]):
                self._remove_segment(sid)

        seg_ids = {self._insert_segment(s) for s in self._segments_of(manifold)}
        self._manifold_segs[manifold].update(seg_ids)

    def update_manifold(self, manifold: BaseManifold, changed_points: set[Point]):
        """
        Reâ€‘index only segments that touch the nodes in `changed_nodes`

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

            point = self._find_true_intersection(seg_1, seg_2)

            seg_1_cdist = 0.5 * (
                seg_1.p0.get_cdist(seg_1.manifold.stability)
                + seg_1.p0_seg1.get_cdist(seg_1.manifold.stability)
            )
            seg_2_cdist = 0.5 * (
                seg_2.p0.get_cdist(seg_2.manifold.stability)
                + seg_2.p0_seg1.get_cdist(seg_2.manifold.stability)
            )

            unstable_cdist = (
                seg_1_cdist if seg_1.manifold.stability == "unstable" else seg_2_cdist
            )
            stable_cdist = (
                seg_1_cdist if seg_1.manifold.stability == "stable" else seg_2_cdist
            )

            # keep existing BranchPoint creation for backward compat
            branch_point = BranchPoint(
                2, (unstable_cdist, stable_cdist), point[0], point[1]
            )
            self._intersecting_coords[seg1_id] = point
            self._intersecting_coords[seg2_id] = point
            self._intersecting_points[seg1_id] = branch_point
            self._intersecting_points[seg2_id] = branch_point

            if (
                seg_1.manifold.stability == "unstable"
                or seg_2.manifold.stability != "unstable"
            ):
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

            # new Intersection object
            # intersection = Intersection.from_segments(
            #     coords=tuple(point),
            #     unstable_cdist=unstable_cdist,
            #     stable_cdist=stable_cdist,
            #     seg1_id=seg1_id,
            #     seg2_id=seg2_id,
            # )
            self._intersections.append(intersection)
            self._intersection_by_seg[seg1_id] = intersection
            self._intersection_by_seg[seg2_id] = intersection
            new_intersections.append(intersection)

            self._processed_pairs.add(seg_id_pair)

        return new_intersections

    # ------------- internal helpers -----------------
    def _insert_segment(self, seg: _Segment) -> int:
        """
        Inserts a segment into the rtree and local dictionary

        Parameters:
            seg (_Segment): segment to be inserted
        """
        # edge key defined by the id of two points
        edge_key = frozenset((id(seg.p0), id(seg.p0_seg1)))
        if edge_key in self._edge_seen:
            return None  # already indexed â†’ do NOT duplicate

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
        except Exception:
            pass  # Segment may have already been deleted from rtree

        # Safely remove from manifold's segment set
        if seg.manifold in self._manifold_segs:
            self._manifold_segs[seg.manifold].discard(sid)

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

    def create_bridges(
        self, for_manifold: Optional[BaseManifold] = None
    ) -> list[Bridge]:

        bridges = {}

        intersecting_segments = []

        manifold_seg_ids = (
            self._manifold_segs.get(for_manifold, set())
            if for_manifold is not None
            else None
        )

        seen_ids = set()
        for sid_pair in self._intersecting_segments:

            if manifold_seg_ids is not None and not (sid_pair & manifold_seg_ids):
                continue  # skip pairs not involving the target manifold

            for sid in sid_pair:
                if sid in seen_ids:
                    continue
                seen_ids.add(sid)

                seg = self._seg_lookup[sid]
                if seg.manifold.stability != "unstable":
                    continue

                cdist = 0.5 * (
                    seg.p0.get_cdist("unstable") + seg.p0_seg1.get_cdist("unstable")
                )

                intersecting_segments.append((cdist, seg))

        intersecting_segments.sort(key=lambda x: x[0])  # sort by cdist

        for i in range(len(intersecting_segments) - 1):
            _, seg1 = intersecting_segments[i]
            _, seg2 = intersecting_segments[i + 1]

            root_point = self._get_nearby_point(seg1, "root")
            tail_point = self._get_nearby_point(seg2, "tail")

            # add in a point nearby the intersection so we do not have
            # dangling tails or bridges that are too short.
            seg1.p0.insert_point_forward(root_point, seg1.manifold.branch_index)
            seg1.p0 = root_point
            seg2.p0_seg1.insert_point_backward(tail_point, seg2.manifold.branch_index)
            seg2.p0_seg1 = tail_point

            bridge = Bridge(
                root=root_point,
                stability=seg1.manifold.stability,
                stretch_param=seg1.manifold.stretch_param,
                fixed_point=seg1.manifold.fixed_point,
                tail=tail_point,
            )

            bridges[(seg1.id, seg2.id)] = bridge

        # update the global bridge structure somehow here
        # if self.bridges is None:
        #     self.bridges = bridges
        # else:
        #     merged = self.bridges.copy()
        #     merged.update(bridges)
        #     self.bridges = merged

        # return list(self.bridges.values())
        bridge_list = list(bridges.values())
        for i in range(len(bridge_list) - 1):
            bridge_list[i].next_bridge = bridge_list[i + 1]
            bridge_list[i + 1].prev_bridge = bridge_list[i]

        return bridge_list

    def _get_nearby_point(self, seg: _Segment, side: Literal["root", "tail"]) -> Point:
        """
        Returns a nearby point to the segment's root point.
        This is used to create a bridge root point.
        """
        if side == "root":
            seg_point = seg.p0
        elif side == "tail":
            seg_point = seg.p0_seg1
        else:
            raise ValueError(f"Invalid side: {side}")

        intersection_coords = self._intersecting_coords[seg.id]

        new_point = self._linear_interpolation(
            intersection_coords, seg_point.get_point(), 0.1
        )

        # this is just an approximation, but it is guaranteed
        # to be between the two points
        new_cdist = (
            seg.p0.get_cdist(seg.manifold.stability)
            + seg.p0_seg1.get_cdist(seg.manifold.stability)
        ) / 2

        new_point = Point(new_point[0], new_point[1], new_cdist)

        return new_point

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
