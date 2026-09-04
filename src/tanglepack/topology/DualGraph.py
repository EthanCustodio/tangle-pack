"""
The dual graph of a trellis: faces on one side, stable arcs on the other.

The arrangement cuts the plane into faces; the stable arcs are the walls between
them. This module turns that picture into a bipartite graph — one
:class:`FaceNode` per face of the plane and one :class:`ArcNode` per stable arc,
with an edge whenever an arc bounds a face — and labels every arc node with the
partition element lying on each of its two sides. That labelling is what makes
the graph symbolic: a walk from face to face through arcs spells a word in the
tangle's bridge classes.

Three things happen on top of the raw arrangement:

* **Merging.** The face traversal walks each connected component's outer
  boundary as its own cycle, so a two-tangle arrangement reports two "outer"
  faces that are really one piece of plane (the nested one being a piece of the
  face of the outer tangle that swallows it). Union-find glues each component's
  outer face either onto the innermost face of another component that
  geometrically contains it, or onto one global unbounded node.
* **Elements.** The element on a side of an arc is the one whose span OWNS the
  arc's midpoint canonical distance, so a singleton element pinched at a
  crossing is never mistaken for the arc's element.
* **Fill.** An arc node is *filled* when an image bridge may cross the stable
  manifold there at a crossing this trellis has not computed yet — the outward
  stretch of each stable branch beyond the image of the previous branch's last
  crossing. A filled node is passable in a walk; a hollow one is a wall.

Dev Notes:

* The stored graph is deliberately bipartite. The natural object for a walk is
  the "through-face edge" ``(arc_in, face, arc_out)``, but storing those is
  quadratic in a face's boundary length; deriving them from the face hub is
  linear and answers the same question.
* A face node's side of an arc is combinatorial, not geometric: the arrangement
  traverses every face with the face on the RIGHT of each half-edge, and a
  stable arc traversed ``hi -> lo`` runs in the stable dynamical direction (the
  direction the partition's ``left``/``right`` is defined against). So
  ``reverse`` alone gives the side, with no angles and no polylines. The tests
  pin that against :func:`~.StablePartition._side_of` on every closed region.
* The fill rule reads the LAST crossing of each stable branch and the image of
  the last crossing of the branch one map step back. With every branch trimmed
  at the strong-pip orbit those two coincide on all but the pip's own branch,
  which is why the derived fill is ``(f^k(q0), q0]`` there and empty elsewhere;
  ``strong_pips=`` makes the constructor check exactly that and log it.

Open question: a merged outer node is a single node even though the piece of
plane it stands for is not simply connected (it wraps around the inner tangle).
Walks through it are therefore permitted between any two of its walls, which is
correct for the symbolic dynamics of one tangle but over-generous the moment two
tangles are genuinely linked by heteroclinic crossings.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Iterable, Literal, Optional, Sequence, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from ..numerics.geometry import point_in_polygon, signed_polygon_area
from .BridgeClass import _index_partitions
from .StablePartition import owns_cdist
from .TopologyResults import (
    Arc,
    ElementRef,
    OPPOSITE_SIDE,
    PartitionInterval,
    Region,
    Side,
    StablePartitionResult,
)

if TYPE_CHECKING:  # pragma: no cover - typing only
    import networkx

    from ..numerics.Bridge import BridgeId
    from ..numerics.FixedPoint import FixedPoint
    from ..numerics.Intersection import ManifoldKey
    from .Arrangement import Arrangement
    from .Trellis import Trellis
    from .TrellisBranch import TrellisBranch

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

#: The kinds of face a merged node can stand for.
FaceKind = Literal["region", "open", "outer"]


def _position_of(
    fixed_point: "FixedPoint", fixed_points: Sequence["FixedPoint"]
) -> int:
    """The index of a fixed point within a list, matched by IDENTITY."""
    for index, candidate in enumerate(fixed_points):
        if candidate is fixed_point:
            return index
    return len(fixed_points)


@dataclass
class ArcNode:
    """
    One stable arc of the arrangement, with the element on each of its sides.

    An arc runs between two consecutive crossings of one stable branch, so no
    partition boundary falls strictly inside it and each side carries exactly
    one element: the one owning the arc's midpoint canonical distance.

    Attributes:
        arc: The arc in its natural direction (``reverse=False``, lo to hi).
        key: The arc's :attr:`~.TopologyResults.Arc.edge_key` — this node's
            identity.
        branch_key: Manifold key of the stable branch the arc lies on.
        lo_cdist: Stable canonical distance of the arc's lower crossing.
        hi_cdist: Stable canonical distance of the arc's upper crossing.
        left: The element on the left of the branch's dynamical direction.
        right: The element on its right.
        filled: True when an image bridge may cross the stable manifold here
            (see :meth:`DualGraph._fill`); a filled node is passable in a walk.
        faces: The face node on each side, keyed by side. Both entries are
            filled once the graph is built; a cut arc whose two sides belong to
            one merged face node has the same node in both.
    """

    arc: Arc
    key: tuple
    branch_key: "ManifoldKey"
    lo_cdist: float
    hi_cdist: float
    left: ElementRef
    right: ElementRef
    filled: bool = False
    faces: dict[Side, "FaceNode"] = field(default_factory=dict)

    @property
    def mid_cdist(self) -> float:
        """The canonical distance halfway along the arc."""
        return 0.5 * (self.lo_cdist + self.hi_cdist)

    def element_on(self, side: Side) -> ElementRef:
        """
        The partition element on one side of this arc.

        Args:
            side: ``"left"`` or ``"right"`` of the stable dynamical direction.

        Returns:
            That side's :class:`~.TopologyResults.ElementRef`.

        Raises:
            ValueError: If ``side`` is neither "left" nor "right".
        """
        if side == "left":
            return self.left
        if side == "right":
            return self.right
        raise ValueError(f"side must be 'left' or 'right', got {side!r}")

    def face_on(self, side: Side) -> "FaceNode":
        """
        The face node lying on one side of this arc.

        Args:
            side: ``"left"`` or ``"right"``.

        Returns:
            The :class:`FaceNode` on that side.

        Raises:
            ValueError: If ``side`` is not a side, or the graph has not attached
                a face there (which the constructor asserts cannot happen).
        """
        if side not in ("left", "right"):
            raise ValueError(f"side must be 'left' or 'right', got {side!r}")
        face = self.faces.get(side)
        if face is None:
            raise ValueError(
                f"arc node {self.key} has no face on its {side} side; the dual "
                "graph was not fully built"
            )
        return face

    def side_of(self, face: "FaceNode") -> Side:
        """
        Which side of this arc a face lies on.

        Args:
            face: A face node adjacent to this arc.

        Returns:
            ``"left"`` or ``"right"``. When one merged node lies on BOTH sides
            (a cut arc inside one piece of plane) ``"left"`` is returned.

        Raises:
            ValueError: If ``face`` is not adjacent to this arc.
        """
        for side in ("left", "right"):
            if self.faces.get(side) is face:
                return side  # type: ignore[return-value]
        raise ValueError(
            f"face node {face.index} does not bound arc node {self.key}"
        )

    def other_face(self, face: "FaceNode") -> "FaceNode":
        """
        The face across this arc from one of its faces.

        Args:
            face: A face node adjacent to this arc.

        Returns:
            The face node on the other side — ``face`` itself when the same
            merged node lies on both sides.

        Raises:
            ValueError: If ``face`` is not adjacent to this arc.
        """
        side = self.side_of(face)
        return self.face_on("right" if side == "left" else "left")

    def midpoint(self, trellis: "Trellis") -> Optional[NDArray[np.float64]]:
        """
        The geometric midpoint of the arc.

        Args:
            trellis: The trellis resolving the arc's crossings and manifold.

        Returns:
            The ``(2,)`` midpoint of the arc's polyline (see
            :meth:`~.TopologyResults.Arc.midpoint`), or None for a degenerate
            arc.
        """
        return self.arc.midpoint(trellis)

    def __repr__(self) -> str:
        state = "filled" if self.filled else "hollow"
        return (
            f"ArcNode({self.arc.lo_id}-{self.arc.hi_id}, {state}, "
            f"L={self.left.label}, R={self.right.label})"
        )


@dataclass
class FaceNode:
    """
    One face of the plane: a face of the arrangement, or several merged into one.

    Attributes:
        index: Position of this node in :attr:`DualGraph.face_nodes`.
        faces: The arrangement faces this node stands for. More than one only
            for a merged outer node (see :meth:`DualGraph._merge_faces`).
        kind: ``"region"`` for a single closed minimal face, ``"outer"`` for any
            node that absorbed a component's outer face (or a containing face),
            ``"open"`` for anything else holding an open face.
        is_unbounded: True on the one node standing for the unbounded plane.
        arcs: The stable arcs on this node's boundary, as
            ``(arc node, the side of the ARC this face lies on)``.
    """

    index: int
    faces: list[Region]
    kind: FaceKind
    is_unbounded: bool = False
    arcs: list[tuple[ArcNode, Side]] = field(default_factory=list)

    @property
    def arc_nodes(self) -> list[ArcNode]:
        """The distinct arc nodes on this face's boundary, in traversal order."""
        seen: set[tuple] = set()
        nodes: list[ArcNode] = []
        for node, _side in self.arcs:
            if node.key in seen:
                continue
            seen.add(node.key)
            nodes.append(node)
        return nodes

    def side_of(self, arc: ArcNode) -> Side:
        """
        Which side of one of its boundary arcs this face lies on.

        Args:
            arc: An arc node on this face's boundary.

        Returns:
            ``"left"`` or ``"right"`` (``"left"`` when this face lies on both
            sides of the arc).

        Raises:
            ValueError: If ``arc`` does not bound this face.
        """
        return arc.side_of(self)

    @property
    def bridge_ids(self) -> list["BridgeId"]:
        """The bridge ids of every unstable arc on the member faces' boundaries."""
        found: list["BridgeId"] = []
        seen: set["BridgeId"] = set()
        for face in self.faces:
            for bridge_id in face.bridge_ids:
                if bridge_id in seen:
                    continue
                seen.add(bridge_id)
                found.append(bridge_id)
        return found

    @property
    def is_region(self) -> bool:
        """True when this node is a single closed minimal region of the tangle."""
        return self.kind == "region"

    def __repr__(self) -> str:
        unbounded = ", unbounded" if self.is_unbounded else ""
        return (
            f"FaceNode({self.index}, {self.kind}{unbounded}, "
            f"{len(self.faces)} face(s), {len(self.arcs)} stable wall(s))"
        )


class _UnionFind:
    """Union-find over face indices, with a synthetic root for the unbounded face."""

    def __init__(self) -> None:
        """Build an empty structure; items are added on first use."""
        self._parent: dict[int, int] = {}

    def find(self, item: int) -> int:
        """The representative of an item's class (adding the item if new)."""
        self._parent.setdefault(item, item)
        while self._parent[item] != item:
            self._parent[item] = self._parent[self._parent[item]]
            item = self._parent[item]
        return item

    def union(self, a: int, b: int) -> None:
        """Merge the classes of two items."""
        root_a, root_b = self.find(a), self.find(b)
        if root_a != root_b:
            self._parent[root_a] = root_b


class DualGraph:
    """
    The bipartite dual of a trellis's arrangement: faces, stable arcs, and fill.

    Build one over the ALL-fixed-points arrangement and the per-fixed-point
    stable partitions::

        dual = DualGraph(
            session.arrangement(),
            session.trellis(fp).stable_partitions,
            strong_pips=[session.trellis(fp).strong_pip],
        )

    Attributes:
        trellis: The trellis the arrangement was built over.
        arrangement: The arrangement whose faces and arcs this graph duals.
        partitions: The stable partitions, indexed by ``(branch_key, side)``.
        arc_nodes: One :class:`ArcNode` per stable arc, keyed by
            :attr:`~.TopologyResults.Arc.edge_key`, in branch-then-cdist order.
        face_nodes: The merged faces; ``face_nodes[i].index == i``.
        unbounded: The one face node standing for the unbounded plane.
        fill_segments: Per stable branch, the half-open canonical-distance span
            ``(c_img, cdist(T_j)]`` whose arcs are filled. A span with equal
            ends means nothing on that branch is filled.
    """

    def __init__(
        self,
        arrangement: "Arrangement",
        partitions: Iterable[StablePartitionResult],
        *,
        strong_pips: Optional[Iterable[int]] = None,
    ) -> None:
        """
        Args:
            arrangement: The arrangement to build over. Pass the
                ALL-fixed-points one: a face of a nested tangle is bounded by
                arcs of both fixed points.
            partitions: Every stable partition covering the arrangement's stable
                branches, on BOTH sides — normally each per-fixed-point
                trellis's ``stable_partitions``, concatenated.
            strong_pips: Registry ids of the strong pips, one per fixed point,
                used only to CHECK the derived fill against the expected
                ``(f^k(q0), q0]`` and log the comparison. A mismatch is a
                warning, never an error; None skips the check.

        Raises:
            ValueError: If two partitions cover the same ``(branch, side)``; if
                a stable arc's branch has no partition on a side; if no element
                (or more than one) owns an arc's midpoint; or if a connected
                component of the arrangement does not have exactly one
                positive-area face.

        Note:
            Nothing is cached and nothing is grown: the graph is a pure function
            of the arrangement, the partitions and the registry as they stand.
        """
        self.trellis: "Trellis" = arrangement.trellis
        self.arrangement = arrangement
        self.partitions: dict[tuple["ManifoldKey", Side], StablePartitionResult] = (
            _index_partitions(partitions)
        )
        self.arc_nodes: dict[tuple, ArcNode] = {}
        self.face_nodes: list[FaceNode] = []
        self.fill_segments: dict["ManifoldKey", tuple[float, float]] = {}
        self._face_node_of: dict[int, FaceNode] = {}

        self._build_arc_nodes()
        classes = self._merge_faces()
        self._build_face_nodes(classes)
        self._attach_arcs()
        self._fill(strong_pips)
        logger.debug("%s", self.summary())

    # ── B.1: arc nodes ──────────────────────────────────────────────────────

    def _build_arc_nodes(self) -> None:
        """One node per stable edge of the arrangement, ordered branch by branch."""
        found: dict[tuple, ArcNode] = {}
        for face in self.arrangement.faces:
            for arc in face.arcs:
                if arc.kind != "stable" or arc.edge_key in found:
                    continue
                found[arc.edge_key] = self._make_arc_node(arc)

        order = {
            branch.key: index
            for index, branch in enumerate(self.trellis.stable_branches)
        }
        for node in sorted(
            found.values(),
            key=lambda item: (order.get(item.branch_key, len(order)), item.lo_cdist),
        ):
            self.arc_nodes[node.key] = node

    def _make_arc_node(self, arc: Arc) -> ArcNode:
        """Build the node for one stable edge, with the element on each side."""
        natural = arc if not arc.reverse else arc.reversed()
        lo_cdist = float(self.trellis.intersection(arc.lo_id).stable_cdist)
        hi_cdist = float(self.trellis.intersection(arc.hi_id).stable_cdist)
        mid = 0.5 * (lo_cdist + hi_cdist)
        return ArcNode(
            arc=natural,
            key=natural.edge_key,
            branch_key=arc.branch_key,
            lo_cdist=lo_cdist,
            hi_cdist=hi_cdist,
            left=self.element_at(arc.branch_key, "left", mid),
            right=self.element_at(arc.branch_key, "right", mid),
        )

    # ── B.2: face nodes and the merge rule ──────────────────────────────────

    def _traversal_area(self, face: Region) -> float:
        """
        The signed area of a face's traversal polygon, open faces included.

        :attr:`~.TopologyResults.Region.area` refuses an open face — its
        boundary is not a closed curve — but the traversal itself always closes,
        reflecting off each dangling end and retracing that stub, which
        contributes exactly nothing to the shoelace sum. So the polygon is
        stitched here rather than through ``boundary_points``.
        """
        pieces = [arc.polyline(self.trellis) for arc in face.arcs]
        if not pieces:
            return 0.0
        ring = np.vstack(pieces)
        if len(ring) < 3:
            return 0.0
        return float(signed_polygon_area(ring))

    def _nodes_by_component(self) -> dict[int, list[int]]:
        """The real crossing ids of each connected component of the arrangement."""
        by_component: dict[int, list[int]] = {}
        for face in self.arrangement.faces:
            for corner in face.corners:
                if corner < 0:
                    continue
                component = self.arrangement.component_of(corner)
                if component is None:
                    continue
                ids = by_component.setdefault(component, [])
                if corner not in ids:
                    ids.append(corner)
        return by_component

    def _component_of_face(self, face: Region) -> Optional[int]:
        """The component a face belongs to, from its first real corner."""
        for corner in face.corners:
            if corner >= 0:
                return self.arrangement.component_of(corner)
        return None

    def _merge_faces(self) -> list[list[int]]:
        """
        Union-find over face indices; returns the classes as face-index lists.

        Each component's outer face — the single face of that component whose
        traversal polygon winds positively — is merged with the innermost face
        of ANOTHER component that geometrically contains it, or with a synthetic
        global root when nothing contains it.
        """
        faces = self.arrangement.faces
        root = len(faces)
        union = _UnionFind()
        union.find(root)
        for index in range(len(faces)):
            union.find(index)

        areas = [self._traversal_area(face) for face in faces]
        components = [self._component_of_face(face) for face in faces]
        nodes_by_component = self._nodes_by_component()
        containing = [
            (index, faces[index])
            for index in range(len(faces))
            if faces[index].is_closed and not faces[index].is_minimal
        ]

        for component, node_ids in nodes_by_component.items():
            outer = [
                index
                for index in range(len(faces))
                if components[index] == component and areas[index] > 0.0
            ]
            if len(outer) != 1:
                own_areas = [
                    round(areas[index], 6)
                    for index in range(len(faces))
                    if components[index] == component
                ]
                raise ValueError(
                    f"component {component} has {len(outer)} faces of positive "
                    f"traversal area (expected exactly one outer face); its "
                    f"face areas are {own_areas}"
                )
            host = self._innermost_container(
                outer[0], component, node_ids, containing, components, areas
            )
            union.union(outer[0], root if host is None else host)

        classes: dict[int, list[int]] = {}
        for index in range(len(faces)):
            classes.setdefault(union.find(index), []).append(index)
        classes.setdefault(union.find(root), []).append(root)
        return [
            members
            for _representative, members in sorted(
                classes.items(), key=lambda item: min(item[1])
            )
        ]

    def _innermost_container(
        self,
        outer_index: int,
        component: int,
        node_ids: list[int],
        containing: list[tuple[int, Region]],
        components: list[Optional[int]],
        areas: list[float],
    ) -> Optional[int]:
        """The smallest containing face of another component holding this one."""
        best: Optional[int] = None
        best_area = np.inf
        for index, face in containing:
            if components[index] == component:
                continue
            polygon = face.boundary_points
            if len(polygon) < 3:
                continue
            if not any(
                point_in_polygon(
                    self.trellis.intersection(node_id).get_point(), polygon
                )
                for node_id in node_ids
            ):
                continue
            magnitude = abs(areas[index])
            if magnitude < best_area:
                best, best_area = index, magnitude
        if best is not None:
            logger.debug(
                "merging outer face %s of component %s into containing face %s",
                self.arrangement.faces[outer_index].corners,
                component,
                self.arrangement.faces[best].corners,
            )
        return best

    def _build_face_nodes(self, classes: list[list[int]]) -> None:
        """Turn union-find classes into FaceNodes and index them by member face."""
        faces = self.arrangement.faces
        root = len(faces)
        for index, members in enumerate(classes):
            face_indices = [item for item in members if item != root]
            members_faces = [faces[item] for item in face_indices]
            is_unbounded = root in members
            node = FaceNode(
                index=index,
                faces=members_faces,
                kind=self._kind_of(members_faces, is_unbounded),
                is_unbounded=is_unbounded,
            )
            self.face_nodes.append(node)
            for face in members_faces:
                self._face_node_of[id(face)] = node

        unbounded = [node for node in self.face_nodes if node.is_unbounded]
        assert len(unbounded) == 1, (
            f"expected exactly one unbounded face node, found {len(unbounded)}"
        )
        self.unbounded = unbounded[0]

    def _kind_of(self, faces: list[Region], is_unbounded: bool) -> FaceKind:
        """Classify a merged face node (see :attr:`FaceNode.kind`)."""
        if is_unbounded or len(faces) > 1:
            return "outer"
        if any(face.is_closed and not face.is_minimal for face in faces):
            # A containing face is a level of nesting whether or not this
            # arrangement computed the component it swallows.
            return "outer"
        if any(not face.is_closed for face in faces):
            return "open"
        return "region"

    def _attach_arcs(self) -> None:
        """Record every (face, arc, side) slot on both endpoints of the edge."""
        for face in self.arrangement.faces:
            node = self._face_node_of[id(face)]
            for arc in face.arcs:
                if arc.kind != "stable":
                    continue
                arc_node = self.arc_nodes[arc.edge_key]
                side: Side = "right" if arc.reverse else "left"
                previous = arc_node.faces.get(side)
                assert previous is None or previous is node, (
                    f"two face nodes claim the {side} side of arc "
                    f"{arc_node.key}: {previous} and {node}"
                )
                arc_node.faces[side] = node
                node.arcs.append((arc_node, side))
        for arc_node in self.arc_nodes.values():
            assert "left" in arc_node.faces and "right" in arc_node.faces, (
                f"arc node {arc_node.key} has no face on "
                f"{'left' if 'left' not in arc_node.faces else 'right'}; every "
                "edge of a planar arrangement bounds a face on both sides"
            )

    # ── B.3: fill ───────────────────────────────────────────────────────────

    def _fill(self, strong_pips: Optional[Iterable[int]]) -> None:
        """Mark the arcs beyond the image of the previous branch's last crossing."""
        tol = self.trellis.registry.cdist_tol
        for branch in self.trellis.stable_branches:
            c_img, c_last = self._fill_segment(branch)
            self.fill_segments[branch.key] = (c_img, c_last)
            for node in self.arc_nodes_on(branch.key):
                if node.hi_cdist <= c_img + tol:
                    continue
                node.filled = True
                if node.lo_cdist < c_img - tol:
                    logger.info(
                        "arc %s of branch %s straddles the fill boundary "
                        "%.6g (span %.6g..%.6g); filling it",
                        node.key[1:3],
                        branch.key[1:],
                        c_img,
                        node.lo_cdist,
                        node.hi_cdist,
                    )
        self._check_fill_against_pips(strong_pips)

    def _fill_segment(self, branch: "TrellisBranch") -> tuple[float, float]:
        """The ``(c_img, cdist(T_j))`` fill span of one stable branch."""
        ordered = branch.ordered_ids()
        if len(ordered) <= 1:
            c_last = (
                float(self.trellis.intersection(ordered[-1]).stable_cdist)
                if ordered
                else 0.0
            )
            logger.info(
                "stable branch %s carries %d crossing(s); nothing to fill",
                branch.key[1:],
                len(ordered),
            )
            return c_last, c_last

        last_id = ordered[-1]
        c_last = float(self.trellis.intersection(last_id).stable_cdist)
        fixed_point = branch.key[0]
        previous = self.trellis.branch(fixed_point.advance_key(branch.key, -1))
        if previous is None or not previous.ordered_ids():
            logger.info(
                "stable branch %s has no computed predecessor branch; nothing "
                "to fill",
                branch.key[1:],
            )
            return c_last, c_last

        previous_ordered = previous.ordered_ids()
        if len(previous_ordered) <= 1:
            logger.info(
                "predecessor branch %s of %s has no crossing past its anchor, "
                "so the fill boundary falls at the anchor and the whole branch "
                "is conservatively filled",
                previous.key[1:],
                branch.key[1:],
            )

        previous_last = previous_ordered[-1]
        image_id = self.trellis.iterate(previous_last, 1)
        if image_id is not None:
            c_img = float(self.trellis.intersection(image_id).stable_cdist)
        else:
            c_img = float(
                self.trellis.intersection(previous_last).stable_cdist
                * fixed_point.per_step_beta("stable")
            )
            logger.info(
                "crossing %d (last of branch %s) has no registered image; "
                "scaling its stable cdist to %.6g for the fill boundary of %s",
                previous_last,
                previous.key[1:],
                c_img,
                branch.key[1:],
            )
        return c_img, c_last

    def _check_fill_against_pips(self, strong_pips: Optional[Iterable[int]]) -> None:
        """Log whether the derived fill matches ``(f^k(q0), q0]`` at each pip."""
        if strong_pips is None:
            return
        tol = self.trellis.registry.cdist_tol
        for pip_id in strong_pips:
            if pip_id is None:
                continue
            branch_key = self.trellis.intersection(pip_id).manifold_b_key
            if branch_key is None:
                logger.warning(
                    "strong pip %d carries no stable branch key; cannot check "
                    "the fill against it",
                    pip_id,
                )
                continue
            fixed_point = branch_key[0]
            _key, expected_lo, _from_table = self.trellis.image_cdist(
                pip_id, fixed_point.k_value, "stable"
            )
            expected = (
                float(expected_lo),
                float(self.trellis.intersection(pip_id).stable_cdist),
            )
            derived = self.fill_segments.get(branch_key)
            if derived is None:
                logger.warning(
                    "strong pip %d lies on stable branch %s, which has no fill "
                    "segment",
                    pip_id,
                    branch_key[1:],
                )
                continue
            if abs(derived[0] - expected[0]) <= tol and abs(
                derived[1] - expected[1]
            ) <= tol:
                logger.info(
                    "fill of branch %s is (f^%d(q0), q0] = (%.6g, %.6g] as expected",
                    branch_key[1:],
                    fixed_point.k_value,
                    *expected,
                )
            else:
                logger.warning(
                    "fill of branch %s is (%.6g, %.6g] but the strong pip %d "
                    "expects (f^%d(q0), q0] = (%.6g, %.6g]",
                    branch_key[1:],
                    derived[0],
                    derived[1],
                    pip_id,
                    fixed_point.k_value,
                    expected[0],
                    expected[1],
                )

    # ── lookups ─────────────────────────────────────────────────────────────

    def arc_nodes_on(self, branch_key: "ManifoldKey") -> list[ArcNode]:
        """
        Every arc node lying on one stable branch.

        Args:
            branch_key: Manifold key of the stable branch.

        Returns:
            The arc nodes, ordered by ``lo_cdist`` (anchor outward).
        """
        return sorted(
            (
                node
                for node in self.arc_nodes.values()
                if node.branch_key == branch_key
            ),
            key=lambda node: node.lo_cdist,
        )

    def arc_nodes_of(self, ref: ElementRef) -> list[ArcNode]:
        """
        Every arc node whose element on ``ref``'s side IS ``ref``.

        Args:
            ref: The element reference to look up.

        Returns:
            The arc nodes covered by that element, ordered by ``lo_cdist``. A
            singleton element pinched at a crossing owns no arc and answers
            with an empty list.
        """
        return [
            node
            for node in self.arc_nodes_on(ref.branch_key)
            if node.element_on(ref.side) == ref
        ]

    def element(self, ref: ElementRef) -> PartitionInterval:
        """
        The partition element one reference names.

        Args:
            ref: The element reference.

        Returns:
            The :class:`~.TopologyResults.PartitionInterval`.

        Raises:
            ValueError: If this graph holds no partition for the ref's
                ``(branch, side)``.
            IndexError: If the id is not an element of that partition.
        """
        result = self.partitions.get((ref.branch_key, ref.side))
        if result is None:
            raise ValueError(
                f"no {ref.side!r} partition for stable branch "
                f"{ref.branch_key[1:]}; this dual graph cannot resolve {ref}"
            )
        return result.element(ref.element_id)

    def element_at(
        self, branch_key: "ManifoldKey", side: Side, cdist: float
    ) -> ElementRef:
        """
        The element owning one canonical distance on one side of one branch.

        Args:
            branch_key: Manifold key of the stable branch.
            side: Which side's partition to read.
            cdist: The stable canonical distance to locate.

        Returns:
            The owning element's :class:`~.TopologyResults.ElementRef`.

        Raises:
            ValueError: If no partition covers that ``(branch, side)``, or if
                the distance lies outside the partitioned stretch, or (a broken
                partition) if more than one element owns it.
        """
        result = self.partitions.get((branch_key, side))
        if result is None:
            raise ValueError(
                f"no {side!r} partition covers stable branch {branch_key[1:]}; "
                "partition that branch (or pass its trellis's stable_partitions)"
            )
        tol = self.trellis.registry.cdist_tol
        owners = [
            interval
            for interval in result.intervals
            if owns_cdist(interval, cdist, tol)
        ]
        if len(owners) != 1:
            span = (
                f"{result.intervals[0].lo_cdist:.6g}.."
                f"{result.intervals[-1].hi_cdist:.6g}"
                if result.intervals
                else "empty"
            )
            raise ValueError(
                f"{len(owners)} elements of the {side!r} partition of branch "
                f"{branch_key[1:]} own stable cdist {cdist:.6g} (expected "
                f"exactly one); that partition spans {span}"
            )
        return result.ref(owners[0].element_id)

    def face_of(self, region: Region) -> FaceNode:
        """
        The node standing for one face of the arrangement.

        Args:
            region: A face of this graph's arrangement (a region, an open face
                or a containing face).

        Returns:
            The :class:`FaceNode` it was merged into (itself alone, usually).

        Raises:
            ValueError: If the face does not belong to this arrangement.
        """
        node = self._face_node_of.get(id(region))
        if node is None:
            raise ValueError(
                f"face {region.corners} is not a face of this dual graph's "
                "arrangement"
            )
        return node

    # ── C.3: element images ─────────────────────────────────────────────────

    def element_images(self, n: int = 1) -> dict[ElementRef, list[ElementRef]]:
        """
        Every element of every partition, mapped to the elements covering its image.

        The image span comes from :meth:`~.Trellis.Trellis.image_of_element` —
        the iterate table where it is populated and the canonical-distance
        scaling law where it is not — and the side is carried through the map's
        orientation, flipping once per step when the map reverses it.

        Args:
            n: Number of map steps; negative walks backward. Defaults to 1.

        Returns:
            A dict from each element's :class:`~.TopologyResults.ElementRef` to
            the refs of the elements covering its ``n``-th image, in increasing
            canonical distance. An element with an unbounded end (the anchor or
            the branch end) has no image arc and maps to the empty list.

        Raises:
            ValueError: If an image branch's partition is missing from this
                graph — partition every stable branch before asking for images.
        """
        images: dict[ElementRef, list[ElementRef]] = {}
        results = list(self.partitions.values())
        flip = bool(not self.trellis.orientation_preserving and n % 2)
        for (branch_key, side), result in self.partitions.items():
            fixed_point = branch_key[0]
            image_key = fixed_point.advance_key(branch_key, n)
            image_side: Side = OPPOSITE_SIDE[side] if flip else side
            for interval in result.intervals:
                ref = result.ref(interval.element_id)
                covering = self.trellis.image_of_element(
                    result,
                    interval.element_id,
                    n,
                    partitions=results,
                )
                if covering is None:
                    logger.debug(
                        "element %s has an unbounded end; no %d-image", ref, n
                    )
                    images[ref] = []
                    continue
                images[ref] = [
                    ElementRef(image_key, image_side, element_id)
                    for element_id in covering
                ]
        return images

    # ── views ───────────────────────────────────────────────────────────────

    def graph(self) -> "networkx.Graph":
        """
        This dual graph as a bipartite networkx graph.

        Nodes are ``("arc", edge_key)`` and ``("face", index)``, each carrying
        the :class:`ArcNode` / :class:`FaceNode` on its ``node`` attribute and a
        ``bipartite`` attribute (0 for faces, 1 for arcs). An edge joins a face
        to every stable arc on its boundary, with the side of the ARC the face
        lies on as the ``side`` attribute.

        Returns:
            The graph. It is simple, not a multigraph, so a cut arc with the
            same merged face on both sides contributes ONE edge whose ``side``
            records the first slot seen.
        """
        import networkx as nx

        graph = nx.Graph()
        for face in self.face_nodes:
            graph.add_node(("face", face.index), bipartite=0, node=face)
        for key, node in self.arc_nodes.items():
            graph.add_node(("arc", key), bipartite=1, node=node, filled=node.filled)
        for face in self.face_nodes:
            for arc_node, side in face.arcs:
                if graph.has_edge(("face", face.index), ("arc", arc_node.key)):
                    continue
                graph.add_edge(
                    ("face", face.index), ("arc", arc_node.key), side=side
                )
        return graph

    def summary(self) -> str:
        """A one-line description of the graph's size, kinds and fill."""
        kinds = {"region": 0, "open": 0, "outer": 0}
        for face in self.face_nodes:
            kinds[face.kind] += 1
        filled = sum(1 for node in self.arc_nodes.values() if node.filled)
        fixed_points = self.trellis.fixed_points
        segments = ", ".join(
            f"fp{_position_of(key[0], fixed_points)}:{key[2]}.{key[3]}:"
            f"({lo:.4g},{hi:.4g}]"
            for key, (lo, hi) in self.fill_segments.items()
        )
        return (
            f"DualGraph: {len(self.arc_nodes)} arc nodes ({filled} filled), "
            f"{len(self.face_nodes)} face nodes "
            f"({kinds['region']} region, {kinds['open']} open, "
            f"{kinds['outer']} outer); fill segments {segments}"
        )

    def __repr__(self) -> str:
        return f"<{self.summary()}>"
