"""
The dual graph of a minimal trellis: faces, and the stable edges seen from each side.

The sparse arrangement of a :class:`~tanglepack.topology.MinimalTrellis.MinimalTrellis`
cuts the plane into faces whose walls are the kept bridges and the stable
manifold. This module turns that picture into a graph built for walks:

* one CLOSED node per face of the plane (:class:`FaceNode`), after the same
  outer-face merging as before (union-find glues each component's outer face
  onto the innermost face of another component that geometrically contains it,
  or onto one global unbounded node);
* one OPEN node per (stable edge, side) — a :class:`StableNode` that is a wall:
  it knows the face on its side and the partition element on its side, and a
  walk cannot cross it;
* on the strong pip's stable fundamental segment, the two side nodes of a
  stable edge are UNIFIED into one SOLID node (``sides == ("left", "right")``,
  ``traversable``), which both faces attach to and a walk may cross.

A face node connects to the stable node on ITS side of each stable edge on its
boundary, and to nothing else: bridges are walls that shape faces but never
appear in the dual. So a walk alternates face -> stable node -> face, and only a
unified node leads anywhere.

Payload. A stable node carries, per side it faces, the element of the
partition it was built over (normally the iterated homotopy partition),
that element's interval, the homotopy element it refines and the image
bridge it lies under (``PartitionInterval.parent_element_id`` /
``cut_by``), the elements its element maps onto under one map step, and
the element's NAME (:class:`~.ElementNaming.ElementName`: side-tagged,
anchor outward, ``R_1^2`` for the second child of the first right-hand
homotopy element, plain ``R_3`` for an unsplit one). The names come from
one :class:`~.ElementNaming.ElementNaming` built over the partition's
homotopy family (:attr:`DualGraph.naming`); over a family with no homotopy
parent every element is its own parent and prints plain. A face node
carries its regions, corners and bridge ids and can resolve its image face
on request. Faces deliberately carry NO bridge class: a face may have
several bridges through it, and the class-to-face relation is not one to
one.

Dev Notes:

* The stored graph is bipartite (faces / stable nodes). Face-to-face adjacency
  through a unified node is derived (:meth:`StableNode.other_face`), not
  stored.
* A face's side of an edge is combinatorial: the arrangement traverses every
  face with the face on the RIGHT of each half-edge, and a stable arc traversed
  ``hi -> lo`` runs in the stable dynamical direction, so ``reverse`` alone
  gives the side (:func:`face_side`, which reads an unstable arc the same way
  in ITS dynamical direction, ``lo -> hi``).
* The fundamental segment is declared by the strong pip, exactly as the old
  fill was: on the pip's own stable branch, the arcs meeting
  ``(cdist(f^k(q0)), cdist(q0)]`` with ``k = k_value``, via
  :meth:`~tanglepack.topology.Trellis.Trellis.image_cdist`; nothing on any
  other branch. A straddling arc is unified and logged at INFO; no pip unifies
  nothing and warns.
* The trellis is cut at the pip orbit (``Trellis.cut_at_pip_orbit``, D10), so
  ``S \\ f(S)``, where the image of a bridge may cross the stable manifold
  away from the images of its ends, is exactly the pip window, and the other
  ``k - 1`` branches are walls as a consequence.
* Element images are read from the trellis (iterate table first, scaled cdist
  otherwise); the map is never called during the build. Face images go
  through :meth:`~tanglepack.topology.Arrangement.Arrangement.image_of`,
  which verifies with the real map, so they are resolved LAZILY by
  :meth:`DualGraph.image_face` and never at build time.
* A sparse component may have no closed face at all (a stable branch with no
  kept bridge). Its single open face is merged into the unbounded node with a
  WARNING rather than raising.
* A component's outer face is read off the rotation system, never an area
  (INT-10, WP10): it is the face whose sector at the component's EASTERNMOST
  polyline vertex holds due east, since everything else lies west of that
  vertex. The positive-traversal-area rule it replaces raised on the
  self-intersecting slivers of deep blasts, where no area is trustworthy.
  Nesting is then one point-in-polygon test per component
  (:meth:`DualGraph._innermost_container`), the chalkboard's one containment
  test per tangle.
* Markers and enclosure (D13, WP9b). The map is a homeomorphism, so a closed
  curve made of manifold pieces holds a marker -- a periodic point or a hole
  -- exactly when its image holds the marker's image. :attr:`DualGraph.markers`
  places them combinatorially: a periodic point at its anchor vertex (``f``
  moves it to the next orbit point), a hole BESIDE ITS OWN BRIDGE on its
  bridge side (``f`` moves hole ``(origin, n)`` to ``(origin, n + 1)``). The
  hole sits beside its own bridge, not in the trellis lobe between the bridge
  and the stable manifold, because a propagated hole lies beside the image
  SUB-ARC of its containing bridge (INT-2): the lobe-corner placement got 2 of
  83 member images wrong, all iterate ``-1`` holes. An unstable arc is never
  cut by a curve (curves cross the stable manifold only), so beside it is one
  piece of one face. A hole whose image is no hole (a reference, or the last
  ``+j`` of a period-k orbit) has that image placed beside the registered image
  chain of its bridge, at the image of the crossing it hugs (table only;
  nothing when the chain is not kept). :meth:`DualGraph.enclosed` floods the
  faces cut by a curve; a curve that crosses a face with several boundary
  components (a merged face, or the unbounded one with infinity) may pass each
  component on either side, so every placement is returned and the caller
  asks whether SOME routing holds the right markers: the routing is free
  because nothing but the markers and the walls is known inside a face.
"""

from __future__ import annotations

import itertools
import logging
from dataclasses import dataclass, field
from functools import cached_property
from typing import Hashable, Iterable, Literal, Optional, Sequence, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from ..numerics.geometry import point_in_polygon, signed_polygon_area
from .ElementNaming import ElementName, ElementNaming
from .MinimalTrellis import image_chain
from .TopologyResults import (
    Arc,
    ElementRef,
    OPPOSITE_SIDE,
    PartitionInterval,
    Region,
    Side,
    side_after,
)

if TYPE_CHECKING:  # pragma: no cover - typing only
    import networkx

    from ..numerics.Bridge import BridgeId
    from ..numerics.FixedPoint import FixedPoint
    from ..numerics.Intersection import ManifoldKey
    from .Arrangement import Arrangement
    from .MinimalTrellis import MinimalTrellis
    from .PartitionFamily import PartitionFamily
    from .Trellis import Trellis

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

#: The kinds of face a merged node can stand for.
FaceKind = Literal["region", "open", "outer"]

#: Which side(s) of its stable edge a stable node faces.
NodeSide = Literal["left", "right", "both"]

#: A point the map carries along an orbit (D13): ``("orbit", fixed-point
#: position, orbit index)`` for a periodic point, ``("hole", origin, iterate)``
#: for a hole.
Marker = tuple

#: Where a curve runs through a face: ``(edge key in, side in, edge key out,
#: side out)``, each end at the middle of a stable edge, on the side of the edge
#: the face lies on.
Chord = tuple


def face_side(arc: Arc) -> Side:
    """
    The side of an arc the face traversing it lies on, in the arc's dynamical direction.

    A face lies on the RIGHT of every half-edge it traverses, and the dynamical
    direction runs ``hi -> lo`` on a stable arc (toward the anchor) and
    ``lo -> hi`` on an unstable one (away from it).
    """
    return "right" if arc.reverse == (arc.kind == "stable") else "left"


def _position_of(
    fixed_point: "FixedPoint", fixed_points: Sequence["FixedPoint"]
) -> int:
    """The index of a fixed point within a list, matched by IDENTITY."""
    for index, candidate in enumerate(fixed_points):
        if candidate is fixed_point:
            return index
    return len(fixed_points)


@dataclass(eq=False)
class StableNode:
    """
    One stable edge of the arrangement seen from one side (or, unified, from both).

    Attributes:
        arc: The edge in its natural direction (``reverse=False``, lo to hi).
        edge_key: The edge's :attr:`~.TopologyResults.Arc.edge_key`.
        key: This node's identity: ``(edge_key, node_side)`` with ``node_side``
            ``"left"``, ``"right"`` or ``"both"``. What a walk records.
        branch_key: Manifold key of the stable branch the edge lies on.
        lo_cdist: Stable canonical distance of the edge's lower crossing.
        hi_cdist: Stable canonical distance of the edge's upper crossing.
        sides: The side(s) this node faces: one of ``("left",)``,
            ``("right",)`` or ``("left", "right")``.
        traversable: True iff the node is unified (on the fundamental
            segment): a walk may cross the stable manifold here.
        fundamental_span: The ``(lo, hi)`` fundamental-segment span the edge
            meets, or None off the segment.
        elements: Per faced side, the partition element owning the edge's
            midpoint on that side.
        intervals: Per faced side, that element's interval.
        parents: Per faced side, the homotopy element the element refines
            (a ref into the HOMOTOPY family), or None on an unrefined
            partition.
        cut_by: Per faced side, the image bridge the element lies under, or
            None.
        names: Per faced side, the element's
            :class:`~.ElementNaming.ElementName` (``R_1^2``), filled by
            :meth:`DualGraph._name_stable_nodes`.
        faces: The face node on each faced side, filled once the graph is
            built.
    """

    arc: Arc
    edge_key: tuple
    key: tuple
    branch_key: "ManifoldKey"
    lo_cdist: float
    hi_cdist: float
    sides: tuple[Side, ...]
    traversable: bool = False
    fundamental_span: Optional[tuple[float, float]] = None
    elements: dict[Side, ElementRef] = field(default_factory=dict)
    intervals: dict[Side, PartitionInterval] = field(default_factory=dict)
    parents: dict[Side, Optional[ElementRef]] = field(default_factory=dict)
    cut_by: dict[Side, Optional["BridgeId"]] = field(default_factory=dict)
    names: dict[Side, ElementName] = field(default_factory=dict)
    faces: dict[Side, "FaceNode"] = field(default_factory=dict)

    @property
    def label(self) -> str:
        """
        The node's element name(s), one per faced side, joined by ``" | "``.

        A unified node reads ``"R_1^2 | L_3"`` (left side first, matching
        :attr:`sides`); a side node reads its one name. A side with no name
        falls back to the element's :attr:`~.TopologyResults.ElementRef.label`.
        """
        return " | ".join(
            self.names[side].text if side in self.names else self.elements[side].label
            for side in self.sides
        )

    @property
    def node_side(self) -> NodeSide:
        """``"both"`` for a unified node, else the one side it faces."""
        return "both" if len(self.sides) == 2 else self.sides[0]

    @property
    def is_unified(self) -> bool:
        """True when this one node stands for both sides of its edge."""
        return len(self.sides) == 2

    @property
    def mid_cdist(self) -> float:
        """The canonical distance halfway along the edge."""
        return 0.5 * (self.lo_cdist + self.hi_cdist)

    def _require_side(self, side: Side) -> None:
        if side not in self.sides:
            raise ValueError(
                f"stable node {self.key} faces {self.sides}, not {side!r}"
            )

    def element_on(self, side: Side) -> ElementRef:
        """
        The partition element on one faced side.

        Args:
            side: ``"left"`` or ``"right"`` of the stable dynamical direction.

        Returns:
            That side's :class:`~.TopologyResults.ElementRef`.

        Raises:
            ValueError: If this node does not face ``side``.
        """
        self._require_side(side)
        return self.elements[side]

    def face_on(self, side: Side) -> "FaceNode":
        """
        The face node on one faced side.

        Args:
            side: ``"left"`` or ``"right"``.

        Returns:
            The :class:`FaceNode` there.

        Raises:
            ValueError: If this node does not face ``side`` or the graph has
                not attached a face there (asserted impossible after build).
        """
        self._require_side(side)
        face = self.faces.get(side)
        if face is None:
            raise ValueError(
                f"stable node {self.key} has no face on its {side} side; the dual "
                "graph was not fully built"
            )
        return face

    def side_of(self, face: "FaceNode") -> Side:
        """
        Which side of this edge a face lies on.

        Args:
            face: A face node adjacent to this node.

        Returns:
            ``"left"`` or ``"right"`` (``"left"`` when one merged node lies on
            both sides of a unified edge).

        Raises:
            ValueError: If ``face`` is not adjacent to this node.
        """
        for side in self.sides:
            if self.faces.get(side) is face:
                return side
        raise ValueError(f"face node {face.index} does not bound stable node {self.key}")

    def other_face(self, face: "FaceNode") -> "FaceNode":
        """
        The face across this edge from one of its faces — a step of a walk.

        Args:
            face: A face node adjacent to this node.

        Returns:
            The face on the other side (``face`` itself when one merged node
            lies on both sides).

        Raises:
            ValueError: If this node is a wall (not unified): there is no
                crossing here. Or if ``face`` is not adjacent.
        """
        if not self.is_unified:
            raise ValueError(
                f"stable node {self.key} is a wall; a walk cannot cross the stable "
                "manifold off the fundamental segment"
            )
        side = self.side_of(face)
        return self.face_on(OPPOSITE_SIDE[side])

    def midpoint(self, trellis: "Trellis") -> Optional[NDArray[np.float64]]:
        """
        The geometric midpoint of the edge.

        Args:
            trellis: The trellis resolving the edge's crossings and manifold
                (the minimal trellis's ``sparse`` snapshot or the full one).

        Returns:
            The ``(2,)`` midpoint of the polyline, or None for a degenerate arc.
        """
        return self.arc.midpoint(trellis)

    def __repr__(self) -> str:
        state = "solid" if self.is_unified else f"open/{self.sides[0]}"
        labels = ", ".join(f"{s[0].upper()}={ref.label}" for s, ref in self.elements.items())
        return f"StableNode({self.arc.lo_id}-{self.arc.hi_id}, {state}, {labels})"


@dataclass(eq=False)
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
        stable_nodes: The stable nodes on this face's boundary, as ``(node,
            the side of the EDGE this face lies on)``, in traversal order.
    """

    index: int
    faces: list[Region]
    kind: FaceKind
    is_unbounded: bool = False
    stable_nodes: list[tuple[StableNode, Side]] = field(default_factory=list)

    @property
    def stable_node_list(self) -> list[StableNode]:
        """The distinct stable nodes on this face's boundary, in traversal order."""
        seen: set[tuple] = set()
        nodes: list[StableNode] = []
        for node, _side in self.stable_nodes:
            if node.key in seen:
                continue
            seen.add(node.key)
            nodes.append(node)
        return nodes

    @property
    def corners(self) -> tuple[int, ...]:
        """The distinct real crossing ids at the corners of the member faces, sorted."""
        return tuple(
            sorted({corner for face in self.faces for corner in face.corners if corner >= 0})
        )

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

    @property
    def exits(self) -> list[StableNode]:
        """The unified (traversable) stable nodes on this face's boundary."""
        return [node for node in self.stable_node_list if node.is_unified]

    def side_of(self, node: StableNode) -> Side:
        """
        Which side of one of its boundary edges this face lies on.

        Args:
            node: A stable node on this face's boundary.

        Returns:
            ``"left"`` or ``"right"``.

        Raises:
            ValueError: If ``node`` does not bound this face.
        """
        return node.side_of(self)

    def __repr__(self) -> str:
        unbounded = ", unbounded" if self.is_unbounded else ""
        return (
            f"FaceNode({self.index}, {self.kind}{unbounded}, "
            f"{len(self.faces)} face(s), {len(self.stable_node_list)} stable node(s), "
            f"{len(self.exits)} exit(s))"
        )


class _UnionFind:
    """Union-find over hashable items (face indices, face pieces); items are added on first use."""

    def __init__(self) -> None:
        """Build an empty structure."""
        self._parent: dict[Hashable, Hashable] = {}

    def find(self, item: Hashable) -> Hashable:
        """The representative of an item's class (adding the item if new)."""
        self._parent.setdefault(item, item)
        while self._parent[item] != item:
            self._parent[item] = self._parent[self._parent[item]]
            item = self._parent[item]
        return item

    def union(self, a: Hashable, b: Hashable) -> None:
        """Merge the classes of two items."""
        root_a, root_b = self.find(a), self.find(b)
        if root_a != root_b:
            self._parent[root_a] = root_b


class DualGraph:
    """
    The dual graph of a minimal trellis over a partition family.

    Build one from the session's pieces::

        dual = DualGraph(
            session.minimal_trellis(),
            session.iterated_partition(),
            strong_pips=[session.trellis(fp).strong_pip],
        )

    Attributes:
        minimal: The minimal trellis dualled.
        trellis: Its full trellis (registry lookups, iterates, strong pips).
        arrangement: The minimal trellis's sparse arrangement.
        partition: The partition family the stable nodes' elements come from.
        naming: The :class:`~.ElementNaming.ElementNaming` the stable nodes'
            names come from: over ``partition.homotopy`` and ``partition``
            when the partition is an iterated homotopy partition with its
            homotopy family attached, else over ``partition`` alone (every
            element its own parent, superscript 1).
        stable_nodes: Every :class:`StableNode`, keyed by :attr:`StableNode.key`,
            in branch-then-cdist order (left before right before both).
        face_nodes: The merged faces; ``face_nodes[i].index == i``.
        unbounded: The one face node standing for the unbounded plane.
        fundamental_segments: Per strong-pip stable branch, the half-open span
            ``(cdist(f^k(q0)), cdist(q0)]`` whose edges are unified. Empty when
            no pip was supplied.
    """

    def __init__(
        self,
        minimal: "MinimalTrellis",
        partition: "PartitionFamily",
        *,
        strong_pips: Optional[Iterable[int]] = None,
    ) -> None:
        """
        Args:
            minimal: The minimal trellis (its ``arrangement`` is dualled).
            partition: The family resolving each stable node's element on each
                side; normally the session's iterated homotopy partition. It
                must cover every stable branch of the arrangement on both
                sides.
            strong_pips: Registry ids of the strong pips, one per fixed point
                (``None`` entries are skipped). Each declares the fundamental
                segment of its own stable branch. None or empty unifies
                nothing and logs a WARNING.

        Raises:
            ValueError: If a stable edge's branch has no partition on a side,
                if no element (or more than one) owns an edge's midpoint, or if a
                strong pip carries no stable branch key or two pips lie on one
                branch.
        """
        self.minimal = minimal
        self.trellis: "Trellis" = minimal.trellis
        self.arrangement: "Arrangement" = minimal.arrangement
        self.partition = partition
        self.naming: Optional[ElementNaming] = None
        self._naming_is_iterated = False
        self.stable_nodes: dict[tuple, StableNode] = {}
        self.face_nodes: list[FaceNode] = []
        self.fundamental_segments: dict["ManifoldKey", tuple[float, float]] = {}
        self._face_node_of: dict[int, FaceNode] = {}
        self._nodes_of_edge: dict[tuple, dict[Side, StableNode]] = {}

        self._fundamental_segments_from(strong_pips)
        self._build_stable_nodes()
        self._name_stable_nodes()
        classes = self._merge_faces()
        self._build_face_nodes(classes)
        self._attach()
        self._check()
        logger.debug("%s", self.summary())

    # ── the fundamental segment ─────────────────────────────────────────────

    def _fundamental_segments_from(self, strong_pips: Optional[Iterable[int]]) -> None:
        """Record ``(f^k(q0), q0]`` on each strong pip's own branch.

        The same pip gathered twice (from ``trellis(fp)`` and ``trellis([fp])``)
        counts once; two DIFFERENT pips on one branch raise.
        """
        pips = list(dict.fromkeys(pip for pip in (strong_pips or ()) if pip is not None))
        if not pips:
            logger.warning(
                "no strong pip supplied to the dual graph; nothing is unified and "
                "every stable node is a wall (TangleSession.dual_graph() gathers one "
                "pip per classified per-fixed-point trellis)"
            )
            return
        owner: dict["ManifoldKey", int] = {}
        for pip_id in pips:
            branch_key, c_lo, c_hi, from_table = self._pip_segment(pip_id)
            if branch_key in owner:
                raise ValueError(
                    f"strong pips {owner[branch_key]} and {pip_id} both lie on "
                    f"stable branch {branch_key[1:]}; a branch has one fundamental segment"
                )
            owner[branch_key] = pip_id
            self.fundamental_segments[branch_key] = (c_lo, c_hi)
            logger.info(
                "fundamental segment of branch %s is (f^%d(q0), q0] = (%.6g, %.6g] "
                "from strong pip %d (%s)",
                branch_key[1:],
                branch_key[0].k_value,
                c_lo,
                c_hi,
                pip_id,
                "iterate table" if from_table else "scaled cdist",
            )

    def _pip_segment(self, pip_id: int) -> tuple["ManifoldKey", float, float, bool]:
        """``(branch, cdist(f^k(q0)), cdist(q0), from_table)`` for one strong pip."""
        crossing = self.trellis.intersection(pip_id)
        branch_key = crossing.manifold_b_key
        if branch_key is None:
            raise ValueError(
                f"strong pip {pip_id} carries no stable branch key; it cannot "
                "declare a fundamental segment"
            )
        fixed_point = branch_key[0]
        image_key, c_lo, from_table = self.trellis.image_cdist(
            pip_id, fixed_point.k_value, "stable"
        )
        # k_value map steps bring a branch back to itself by definition
        # (FixedPoint.advance_key), so a different key means a corrupt table.
        assert image_key == branch_key, (
            f"the {fixed_point.k_value}-th image of strong pip {pip_id} landed on "
            f"branch {image_key[1:]} instead of its own {branch_key[1:]}"
        )
        if not from_table:
            logger.info(
                "strong pip %d has no registered %d-th iterate; its segment boundary "
                "%.6g is the scaled canonical distance",
                pip_id,
                fixed_point.k_value,
                c_lo,
            )
        return branch_key, float(c_lo), float(crossing.stable_cdist), from_table

    def _meets_fundamental_segment(
        self, branch_key: "ManifoldKey", lo: float, hi: float
    ) -> Optional[tuple[float, float]]:
        """The segment an edge meets, or None (straddlers are logged)."""
        segment = self.fundamental_segments.get(branch_key)
        if segment is None:
            return None
        c_lo, c_hi = segment
        if hi <= c_lo or lo >= c_hi:
            return None
        if lo < c_lo or hi > c_hi:
            logger.info(
                "stable edge %.6g..%.6g of branch %s straddles the fundamental "
                "segment (%.6g, %.6g]; unifying it",
                lo,
                hi,
                branch_key[1:],
                c_lo,
                c_hi,
            )
        return segment

    # ── stable nodes ────────────────────────────────────────────────────────

    def _build_stable_nodes(self) -> None:
        """Two open nodes per stable edge, or one solid node on the fundamental segment."""
        found: dict[tuple, Arc] = {}
        for face in self.arrangement.faces:
            for arc in face.arcs:
                if arc.kind == "stable" and arc.edge_key not in found:
                    found[arc.edge_key] = arc if not arc.reverse else arc.reversed()

        order = {branch.key: index for index, branch in enumerate(self.trellis.stable_branches)}
        rank = {"left": 0, "right": 1, "both": 2}
        nodes: list[StableNode] = []
        for arc in found.values():
            lo = float(self.trellis.intersection(arc.lo_id).stable_cdist)
            hi = float(self.trellis.intersection(arc.hi_id).stable_cdist)
            segment = self._meets_fundamental_segment(arc.branch_key, lo, hi)
            sides_groups: list[tuple[Side, ...]] = (
                [("left", "right")] if segment is not None else [("left",), ("right",)]
            )
            by_side: dict[Side, StableNode] = {}
            for sides in sides_groups:
                node_side: NodeSide = "both" if len(sides) == 2 else sides[0]
                node = StableNode(
                    arc=arc,
                    edge_key=arc.edge_key,
                    key=(arc.edge_key, node_side),
                    branch_key=arc.branch_key,
                    lo_cdist=lo,
                    hi_cdist=hi,
                    sides=sides,
                    traversable=segment is not None,
                    fundamental_span=segment,
                )
                # The midpoint is an average, not a registry cdist, but it lies
                # strictly inside a non-singleton element, so interval_at's exact
                # reading answers it.
                mid = node.mid_cdist
                for side in sides:
                    interval = self.partition.interval_at(arc.branch_key, side, mid)
                    node.intervals[side] = interval
                    node.elements[side] = self.partition.ref(
                        arc.branch_key, side, interval.element_id
                    )
                    node.parents[side] = (
                        None
                        if interval.parent_element_id is None
                        else ElementRef(arc.branch_key, side, interval.parent_element_id)
                    )
                    node.cut_by[side] = interval.cut_by
                    by_side[side] = node
                nodes.append(node)
            self._nodes_of_edge[arc.edge_key] = by_side

        for node in sorted(
            nodes,
            key=lambda n: (order.get(n.branch_key, len(order)), n.lo_cdist, rank[n.node_side]),
        ):
            self.stable_nodes[node.key] = node

    # ── payload: element names ──────────────────────────────────────────────

    def _build_naming(self) -> ElementNaming:
        """
        The naming the stable nodes' names come from (see :attr:`naming`).

        Over the partition's homotopy family and the partition itself when the
        partition is an :class:`~.PartitionFamily.IteratedHomotopyPartition`
        with its ``homotopy`` attached; over the partition alone otherwise, in
        which case every element is its own parent and prints plain.
        """
        homotopy = getattr(self.partition, "homotopy", None)
        if homotopy is not None:
            self._naming_is_iterated = True
            return ElementNaming(homotopy, self.partition)
        self._naming_is_iterated = False
        return ElementNaming(self.partition)

    def _name_stable_nodes(self) -> None:
        """Record, per node and side, the name of the element on that side."""
        self.naming = self._build_naming()
        for node in self.stable_nodes.values():
            for side in node.sides:
                ref = node.elements[side]
                node.names[side] = (
                    self.naming.name(ref)
                    if self._naming_is_iterated
                    else self.naming.homotopy_name(ref)
                )

    # ── face nodes and the merge rule ───────────────────────────────────────

    def _traversal_ring(self, face: Region) -> NDArray[np.float64]:
        """
        A face's traversal polygon, open faces included, without repeated vertices.

        :attr:`~.TopologyResults.Region.area` refuses an open face — its
        boundary is not a closed curve — but the traversal itself always closes,
        reflecting off each dangling end and retracing that stub. So the polygon
        is stitched here rather than through ``boundary_points``: consecutive
        arcs share their end vertex, kept once, and the closing vertex is dropped.
        """
        pieces = [np.asarray(arc.polyline(self.arrangement.trellis), dtype=float) for arc in face.arcs]
        if not pieces:
            return np.empty((0, 2))
        ring = np.vstack(pieces)
        ring = ring[np.r_[True, np.any(ring[1:] != ring[:-1], axis=1)]]
        if len(ring) > 1 and np.array_equal(ring[0], ring[-1]):
            ring = ring[:-1]
        return ring

    def _traversal_area(self, face: Region) -> float:
        """The signed area of a face's traversal polygon (a retraced stub adds nothing)."""
        ring = self._traversal_ring(face)
        return float(signed_polygon_area(ring)) if len(ring) >= 3 else 0.0

    def _outer_face(self, members: list[int]) -> int:
        """
        The outer face of one component, from the rotation system (INT-10).

        Every face is traversed with itself on the RIGHT, so at a vertex ``v``
        reached from ``prev`` and left towards ``next`` its sector runs clockwise
        from ``next - v`` to ``prev - v`` (the whole circle at a stub tip, where
        ``prev == next``). At the component's easternmost vertex everything lies
        west, so the one sector there holding due east is the outer face's.

        Args:
            members: The component's face indices.

        Returns:
            The outer face's index.
        """
        faces = self.arrangement.faces
        rings = {index: self._traversal_ring(faces[index]) for index in members}
        east = max((ring[np.argmax(ring[:, 0])] for ring in rings.values() if len(ring)), key=lambda v: v[0])

        def clockwise(start: float, end: float) -> float:
            return (start - end) % (2 * np.pi)

        holders = set()
        for index, ring in rings.items():
            for k in np.flatnonzero(np.all(ring == east, axis=1)):
                after, before = ring[(k + 1) % len(ring)] - east, ring[k - 1] - east
                out, back = np.arctan2(after[1], after[0]), np.arctan2(before[1], before[0])
                sector = clockwise(out, back)
                if sector == 0.0 or 0.0 < clockwise(out, 0.0) < sector:
                    holders.add(index)
        assert len(holders) == 1, (members, sorted(holders), tuple(east))
        return holders.pop()

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

        Each component's outer face (:meth:`_outer_face`, read off the rotation
        system) is merged with the innermost face of ANOTHER component that
        geometrically contains it, or with a synthetic global root when nothing
        contains it. A component with only open faces (a stable branch with no
        kept bridge) has all of its faces merged into the root, with a WARNING.
        """
        faces = self.arrangement.faces
        root = len(faces)
        union = _UnionFind()
        union.find(root)
        for index in range(len(faces)):
            union.find(index)

        components = [self._component_of_face(face) for face in faces]
        nodes_by_component = self._nodes_by_component()
        containing = [
            (index, faces[index])
            for index in range(len(faces))
            if faces[index].is_closed and not faces[index].is_minimal
        ]

        for component, node_ids in nodes_by_component.items():
            members = [index for index in range(len(faces)) if components[index] == component]
            if all(not faces[index].is_closed for index in members):
                logger.warning(
                    "component %s has no closed face (a stable branch with no kept "
                    "bridge?); its %d open face(s) join the unbounded node",
                    component,
                    len(members),
                )
                for index in members:
                    union.union(index, root)
                continue
            outer = self._outer_face(members)
            host = self._innermost_container(outer, component, node_ids, containing, components)
            union.union(outer, root if host is None else host)

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
                point_in_polygon(self.trellis.intersection(node_id).get_point(), polygon)
                for node_id in node_ids
            ):
                continue
            magnitude = abs(self._traversal_area(face))
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

    def _attach(self) -> None:
        """Connect every face to the stable node on ITS side of each boundary edge."""
        for face in self.arrangement.faces:
            node = self._face_node_of[id(face)]
            for arc in face.arcs:
                if arc.kind != "stable":
                    continue
                side = face_side(arc)
                stable_node = self._nodes_of_edge[arc.edge_key][side]
                previous = stable_node.faces.get(side)
                assert previous is None or previous is node, (
                    f"two face nodes claim the {side} side of stable edge "
                    f"{arc.edge_key}: {previous} and {node}"
                )
                stable_node.faces[side] = node
                node.stable_nodes.append((stable_node, side))

    # ── invariants ──────────────────────────────────────────────────────────

    def _check(self) -> None:
        """Assert the structure a walk relies on."""
        for node in self.stable_nodes.values():
            assert node.traversable == node.is_unified
            for side in node.sides:
                assert side in node.faces, (
                    f"stable node {node.key} has no face on its {side} side; every "
                    "edge of a planar arrangement bounds a face on both sides"
                )
        for edge_key, by_side in self._nodes_of_edge.items():
            assert set(by_side) == {"left", "right"}, (
                f"stable edge {edge_key} is missing a side node"
            )
            unified = by_side["left"] is by_side["right"]
            assert unified == (by_side["left"].fundamental_span is not None)
        for face in self.face_nodes:
            for node, side in face.stable_nodes:
                assert node.faces[side] is face

    # ── lookups ─────────────────────────────────────────────────────────────

    def nodes_of_edge(self, edge_key: tuple) -> dict[Side, StableNode]:
        """
        The stable node on each side of one stable edge.

        Args:
            edge_key: An :attr:`~.TopologyResults.Arc.edge_key` of a stable edge.

        Returns:
            ``{"left": node, "right": node}``; both entries are the same
            object for a unified edge.

        Raises:
            KeyError: If the edge is not in this graph.
        """
        return dict(self._nodes_of_edge[edge_key])

    def stable_nodes_on(self, branch_key: "ManifoldKey") -> list[StableNode]:
        """
        Every stable node on one stable branch.

        Args:
            branch_key: Manifold key of the stable branch.

        Returns:
            The nodes ordered by ``lo_cdist`` (anchor outward), left before
            right for a split edge.
        """
        return [
            node for node in self.stable_nodes.values() if node.branch_key == branch_key
        ]

    def stable_nodes_of(self, ref: ElementRef) -> list[StableNode]:
        """
        Every stable node whose element on ``ref``'s side IS ``ref``.

        Args:
            ref: The element reference (in this graph's partition family).

        Returns:
            The nodes, anchor outward. A singleton element pinched at a
            crossing owns no edge and answers with an empty list.
        """
        return [
            node
            for node in self.stable_nodes_on(ref.branch_key)
            if ref.side in node.sides and node.elements[ref.side] == ref
        ]

    @property
    def unified_nodes(self) -> list[StableNode]:
        """The solid (traversable) stable nodes, anchor outward per branch."""
        return [node for node in self.stable_nodes.values() if node.is_unified]

    @property
    def wall_nodes(self) -> list[StableNode]:
        """The open (non-traversable) stable nodes."""
        return [node for node in self.stable_nodes.values() if not node.is_unified]

    def element(self, ref: ElementRef) -> PartitionInterval:
        """
        The partition element one reference names (delegates to the family).

        Args:
            ref: The element reference.

        Returns:
            The :class:`~.TopologyResults.PartitionInterval`.
        """
        return self.partition.element(ref)

    def element_at(self, branch_key: "ManifoldKey", side: Side, cdist: float) -> ElementRef:
        """
        The element owning one canonical distance (delegates to the family).

        Args:
            branch_key: Manifold key of the stable branch.
            side: Which side's partition to read.
            cdist: The stable canonical distance to locate.

        Returns:
            The owning element's :class:`~.TopologyResults.ElementRef`.
        """
        return self.partition.element_at(branch_key, side, cdist)

    def face_of(self, region: Region) -> FaceNode:
        """
        The node standing for one face of the arrangement.

        Args:
            region: A face of this graph's arrangement.

        Returns:
            The :class:`FaceNode` it was merged into (itself alone, usually).

        Raises:
            ValueError: If the face does not belong to this arrangement.
        """
        node = self._face_node_of.get(id(region))
        if node is None:
            raise ValueError(
                f"face {region.corners} is not a face of this dual graph's arrangement"
            )
        return node

    def image_face(self, face: FaceNode, n: int = 1) -> Optional[FaceNode]:
        """
        The face node a region maps onto, resolved on request.

        Goes through :meth:`~tanglepack.topology.Arrangement.Arrangement.image_of`,
        which verifies the candidate against the real map, so it is never run
        during the build.

        Args:
            face: A face node of kind ``"region"``.
            n: Number of map steps (negative for the preimage).

        Returns:
            The image's face node, or None when the face is not a single
            region or the arrangement cannot resolve the image.
        """
        if not face.is_region or len(face.faces) != 1:
            return None
        image = self.arrangement.image_of(face.faces[0], n)
        return None if image is None else self._face_node_of.get(id(image))

    # ── markers and enclosure (D13) ─────────────────────────────────────────

    @cached_property
    def unstable_edges(self) -> dict["BridgeId", tuple]:
        """The edge key of each kept bridge's unstable arc, by :data:`BridgeId`."""
        return {
            tuple(arc.bridge_id): arc.edge_key
            for face in self.arrangement.faces
            for arc in face.arcs
            if arc.kind == "unstable"
        }

    @cached_property
    def markers(self) -> dict[Marker, tuple]:
        """
        Where each marker sits on the sparse arrangement (see the module Dev Notes).

        Returns:
            ``{marker: place}``, ``place`` being ``("vertex", crossing id)`` for a
            periodic point (its anchor) and ``("beside", unstable edge key, side)``
            for a hole (beside its own bridge on its bridge side), or for the
            image of a hole whose image is no hole (beside the kept image chain
            of its bridge, at the image of the crossing it hugs).
        """
        trellis = self.trellis
        places: dict[Marker, tuple] = {}
        for crossing_id, crossing in trellis.registry:
            if crossing.is_anchor and crossing_id in trellis:
                key = crossing.manifold_a_key
                orbit = ("orbit", _position_of(key[0], trellis.fixed_points), key[2])
                places[orbit] = ("vertex", crossing_id)
        holes = {(hole.origin, hole.iterate): hole for hole in self.minimal.holes}
        assert len(holes) == len(self.minimal.holes), "two holes share one orbit step"
        preserving = trellis.orientation_preserving
        for (origin, iterate), hole in holes.items():
            bridge = trellis.bridge_between(*hole.bounding_ids)
            edge = self.unstable_edges.get(bridge.id) if bridge is not None else None
            if edge is None:
                continue
            places[("hole", origin, iterate)] = ("beside", edge, hole.bridge_side)
            if (origin, iterate + 1) in holes:
                continue
            image = trellis.iterate(hole.near_intersection_id, 1)
            for pair in image_chain(trellis, bridge.id) or ():
                if image in pair and pair in self.unstable_edges:
                    side = side_after(hole.bridge_side, 1, preserving)
                    places[("hole", origin, iterate + 1)] = ("beside", self.unstable_edges[pair], side)
        return places

    def image_marker(self, marker: Marker) -> Optional[Marker]:
        """
        The image of a marker under one map step.

        Args:
            marker: A key of :attr:`markers`.

        Returns:
            The next point of a periodic orbit, or the hole one step forward
            on its orbit; None when that hole is neither a hole nor placed.
        """
        if marker[0] == "orbit":
            _kind, position, orbit_index = marker
            period = self.trellis.fixed_points[position].period
            return ("orbit", position, (orbit_index + 1) % period)
        image = ("hole", marker[1], marker[2] + 1)
        return image if image in self.markers else None

    def enclosed(
        self,
        along: Iterable[tuple],
        chords: Iterable[tuple[FaceNode, Chord]] = (),
        along_halves: Iterable[tuple] = (),
    ) -> Optional[list[dict[Marker, Optional[bool]]]]:
        """
        Which markers a closed curve holds.

        The curve runs along whole edges (``along``), along halves of the
        stable edges it ends a chord on (``along_halves``, each named
        ``(edge key, the end vertex of that half)``) and through faces
        (``chords``). The faces are flooded across every edge the curve does
        not run along; a chord cuts its face's ring at the middles of its two
        edges.

        Args:
            along: Edge keys the curve runs along.
            chords: ``(face node, chord)`` per face crossing.
            along_halves: Halves of split stable edges the curve runs along.

        Returns:
            One ``{marker: inside}`` per routing (``None`` for a marker on the
            curve); an EMPTY list when two chords of one face cross, so no
            routing keeps the curve simple; None when the curve is beyond this
            model: two chord ends on one side of one edge (one crossing per
            edge, at its middle), a chord whose ends lie on two boundary
            components, or chords on two components of one face.
        """
        along, halves = set(along), set(along_halves)
        by_face: dict[int, list[Chord]] = {}
        for face, chord in chords:
            by_face.setdefault(face.index, []).append(chord)
        ends = [end for chords_ in by_face.values() for c in chords_ for end in (c[:2], c[2:])]
        if len(set(ends)) < len(ends):
            return None
        split = {edge for edge, _side in ends}
        rings: dict[int, int] = {}
        for index, face_chords in by_face.items():
            ring = self._ring_of(self.face_nodes[index], face_chords)
            if ring is None:
                return None
            rings[index] = ring
        slots = [
            (index, component)
            for index, ring in rings.items()
            for component in self._components(self.face_nodes[index])
            if component != ring
        ]
        options = [range(len(by_face[index]) + 1) for index, _component in slots]
        routings = []
        for choice in itertools.product(*options):
            flood = self._flood(along, halves, by_face, rings, split, dict(zip(slots, choice)))
            if flood is None:
                return []
            union, outside, pieces_at = flood
            held: dict[Marker, Optional[bool]] = {}
            for marker, place in self.markers.items():
                roots = {union.find(piece) for piece in pieces_at.get(place, ())}
                held[marker] = roots.pop() != outside if len(roots) == 1 else None
            routings.append(held)
        return routings

    @staticmethod
    def _components(face: FaceNode) -> list:
        """A face node's boundary components: its rings, and infinity on the unbounded node."""
        return list(range(len(face.faces))) + (["infinity"] if face.is_unbounded else [])

    @staticmethod
    def _ring_of(face: FaceNode, chords: list[Chord]) -> Optional[int]:
        """The one ring of a face holding both ends of every chord, or None."""
        holding = {
            r
            for r, region in enumerate(face.faces)
            if all(
                any(arc.edge_key == edge and face_side(arc) == side for arc in region.arcs)
                for chord in chords
                for edge, side in (chord[:2], chord[2:])
            )
        }
        return holding.pop() if len(holding) == 1 else None

    def _flood(
        self,
        along: set,
        halves: set,
        by_face: dict[int, list[Chord]],
        rings: dict[int, int],
        split: set,
        placement: dict[tuple, int],
    ) -> Optional[tuple[_UnionFind, Hashable, dict[tuple, list]]]:
        """
        Flood the faces cut by one routing of a curve.

        A piece of a face ring is a whole arc, or the half of a split stable
        edge at one end vertex. Pieces of one face are joined unless a chord
        separates them; pieces on the two sides of an edge (or half) the curve
        does not run along are joined across it.

        Returns:
            ``(union-find, the root of infinity, pieces at each marker place)``,
            or None when two chords of one ring cross.
        """
        union = _UnionFind()

        def pieces(index: int, r: int, i: int, arc: Arc) -> list[tuple]:
            if arc.edge_key in split:
                return [(index, r, i, "tail"), (index, r, i, "head")]
            return [(index, r, i, "whole")]

        for face in self.face_nodes:
            components = self._components(face)
            ring = rings.get(face.index)
            for r, region in enumerate(face.faces):
                if r != ring:
                    for i, arc in enumerate(region.arcs):
                        for piece in pieces(face.index, r, i, arc):
                            union.union(piece, ("component", face.index, r))
            if ring is None:
                for component in components:
                    union.union(("component", face.index, component), ("component", face.index, components[0]))
                continue
            regions = self._cut_ring(union, face, ring, by_face[face.index], pieces)
            if regions is None:
                return None
            for component in components:
                if component != ring:
                    union.union(
                        ("component", face.index, component),
                        regions[placement[(face.index, component)]],
                    )

        pieces_at: dict[tuple, list] = {}
        for face in self.face_nodes:
            for r, region in enumerate(face.faces):
                for i, arc in enumerate(region.arcs):
                    for piece in pieces(face.index, r, i, arc):
                        part = piece[3]
                        vertices = {
                            "whole": (arc.tail_id, arc.head_id),
                            "tail": (arc.tail_id,),
                            "head": (arc.head_id,),
                        }[part]
                        for vertex in vertices:
                            pieces_at.setdefault(("vertex", vertex), []).append(piece)
                        if arc.kind == "unstable":
                            pieces_at.setdefault(("beside", arc.edge_key, face_side(arc)), []).append(piece)
                        half = (arc.edge_key, "whole" if part == "whole" else vertices[0])
                        if arc.edge_key not in along and half not in halves:
                            union.union(piece, ("edge",) + half)
        return union, union.find(("component", self.unbounded.index, "infinity")), pieces_at

    @staticmethod
    def _cut_ring(
        union: _UnionFind, face: FaceNode, ring: int, chords: list[Chord], pieces
    ) -> Optional[list[tuple]]:
        """
        Cut one ring by its chords: ``m`` chords that do not cross leave ``m + 1`` regions.

        The chord ends cut the ring into segments, each running from the head
        half of one cut edge to the tail half of the next. A chord from ``t``
        to ``u`` glues the segment ending at ``t`` to the one starting at
        ``u``, and the one starting at ``t`` to the one ending at ``u``.

        Returns:
            One representative per region, in ring order, or None when two
            chords cross.
        """
        arcs = face.faces[ring].arcs
        chord_at = {
            next(i for i, arc in enumerate(arcs) if arc.edge_key == edge and face_side(arc) == side): k
            for k, chord in enumerate(chords)
            for edge, side in (chord[:2], chord[2:])
        }
        cuts = sorted(chord_at)
        before: dict[int, tuple] = {}
        after: dict[int, tuple] = {}
        for j, (start, stop) in enumerate(zip(cuts, cuts[1:] + cuts[:1])):
            segment = ("segment", face.index, j)
            union.union((face.index, ring, start, "head"), segment)
            i = (start + 1) % len(arcs)
            while i != stop:
                for piece in pieces(face.index, ring, i, arcs[i]):
                    union.union(piece, segment)
                i = (i + 1) % len(arcs)
            union.union((face.index, ring, stop, "tail"), segment)
            after[start] = before[stop] = segment
        for k in range(len(chords)):
            t, u = (i for i in cuts if chord_at[i] == k)
            union.union(before[t], after[u])
            union.union(after[t], before[u])
        regions: dict[Hashable, tuple] = {}
        for start in cuts:
            regions.setdefault(union.find(after[start]), after[start])
        if len(regions) != len(chords) + 1:
            return None
        return list(regions.values())

    # ── views ───────────────────────────────────────────────────────────────

    def graph(self) -> "networkx.Graph":
        """
        This dual graph as a bipartite networkx graph.

        Nodes are ``("stable", key)`` (attributes ``node``, ``bipartite=1``,
        ``traversable``) and ``("face", index)`` (``node``, ``bipartite=0``).
        An edge joins a face to the stable node on its side of every stable
        edge on its boundary, with that side as the ``side`` attribute.

        Returns:
            The graph. It is simple, not a multigraph, so a unified edge with
            the same merged face on both sides contributes ONE edge whose
            ``side`` records the first slot seen.
        """
        import networkx as nx

        graph = nx.Graph()
        for face in self.face_nodes:
            graph.add_node(("face", face.index), bipartite=0, node=face)
        for key, node in self.stable_nodes.items():
            graph.add_node(("stable", key), bipartite=1, node=node, traversable=node.traversable)
        for face in self.face_nodes:
            for stable_node, side in face.stable_nodes:
                if graph.has_edge(("face", face.index), ("stable", stable_node.key)):
                    continue
                graph.add_edge(("face", face.index), ("stable", stable_node.key), side=side)
        return graph

    def summary(self) -> str:
        """A one-line description of the graph's size, kinds and fundamental segments."""
        kinds = {"region": 0, "open": 0, "outer": 0}
        for face in self.face_nodes:
            kinds[face.kind] += 1
        unified = len(self.unified_nodes)
        fixed_points = self.trellis.fixed_points
        segments = ", ".join(
            f"fp{_position_of(key[0], fixed_points)}:{key[2]}.{key[3]}:({lo:.4g},{hi:.4g}]"
            for key, (lo, hi) in self.fundamental_segments.items()
        )
        return (
            f"DualGraph: {len(self._nodes_of_edge)} stable edges -> "
            f"{len(self.stable_nodes)} stable nodes ({unified} solid, "
            f"{len(self.stable_nodes) - unified} open), "
            f"{len(self.face_nodes)} face nodes "
            f"({kinds['region']} region, {kinds['open']} open, {kinds['outer']} outer); "
            f"fundamental segments {segments}"
        )

    def __repr__(self) -> str:
        return f"<{self.summary()}>"
