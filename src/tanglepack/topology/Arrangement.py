"""The planar arrangement of a trellis: its nodes, its arcs, and its faces.

Dev Notes — why ``crossing_sign`` alone is enough

A planar subdivision is determined by its *rotation system*: the cyclic order, at
each vertex, of the edges leaving it. Building one normally means computing an
angle per edge and sorting. Here it is free, because the local picture at every
node of a tangle is the same picture.

Exactly four manifold rays leave a crossing, and the fundamental invariant says
which four: one unstable curve passes through it and one stable curve passes
through it, and nothing else can (two unstable manifolds never meet, nor two
stable ones, nor does a manifold meet itself). Name the four rays by the
direction the canonical distance moves along them:

* ``u+`` — along the unstable branch, away from its anchoring periodic point;
* ``u-`` — along the unstable branch, toward the anchor;
* ``s+`` — along the stable branch, away from its anchor;
* ``s-`` — along the stable branch, toward the anchor.

``u+`` and ``u-`` are opposite directions of one smooth curve, so they are
antipodal in angle; likewise ``s+`` and ``s-``. Put ``u+`` at angle 0. Then ``u-``
is at angle pi, and ``s+`` is at some angle theta with ``s-`` at ``theta + pi``.
Transversality (a tangency is not a crossing) rules out ``theta = 0`` and
``theta = pi``, so theta lies strictly in ``(0, pi)`` or strictly in ``(-pi, 0)``,
and there are exactly TWO possible rotation systems:

* ``theta in (0, pi)``, i.e. ``cross(u+, s+) > 0`` — counter-clockwise order
  ``(u+, s+, u-, s-)``;
* ``theta in (-pi, 0)``, i.e. ``cross(u+, s+) < 0`` — then ``s-`` sits at
  ``theta + pi in (0, pi)`` and the counter-clockwise order is
  ``(u+, s-, u-, s+)``.

The sign of ``cross(u+, s+)`` is exactly
:attr:`~tanglepack.numerics.Intersection.Intersection.crossing_sign`, so one
stored integer per crossing fixes the whole rotation system. No angles, no
sorting, and no dependence on how finely the polyline happens to be resolved near
the crossing.

Orientation convention
----------------------
Faces are traced with ``next(e) = the half-edge after twin(e) in the
counter-clockwise rotation at e's head``. That traversal walks a bounded face
CLOCKWISE (and the unbounded face counter-clockwise), so a closed
:class:`~.TopologyResults.Region` has NEGATIVE
:attr:`~.TopologyResults.Region.area`. Callers wanting a magnitude take ``abs``.

Dangling ends and open faces
----------------------------
A computed manifold stops somewhere. The slot in which it stops has no arc, so it
is filled with a VIRTUAL half-edge to a virtual node of degree one. A degree-one
node has a one-element rotation, so ``next`` of the edge into it is the edge back
out: the traversal reflects off the dangling end and carries on, exactly as a face
must wrap around a slit. Any face that picks up such a half-edge is bounded partly
by "we have not computed this far", not by manifold, so it is **open** and is not a
region. Growing the manifolds closes some of them; that is what the Phase 7 growth
driver is for.

Disconnected trellises and faces that are not regions
----------------------------------------------------
A trellis need not be connected. Two fixed points whose manifolds have no COMPUTED
heteroclinic crossing between them are two separate drawings on one plane, and the
face traversal — which only ever walks along arcs — cannot see one from the other.
It then reports each component's outer boundary as its own closed cycle, and one
of those cycles can geometrically SWALLOW the other component: in the nested
period-3 fixture the period-1 tangle's enclosing face is a clean four-crossing
cycle that contains the entire period-3 tangle and most of its faces.

Such a face is closed but it is not a region, because a region "contains no other
manifold piece". :meth:`Arrangement._classify_minimality` finds them — union-find
over the arcs for the components, then a ray cast of each other component's nodes
against each closed face's boundary — and moves them to
:attr:`Arrangement.containing_faces`, leaving :attr:`Arrangement.regions` pairwise
interior-disjoint, which is what a dual graph will need. Within a single component
the test is unnecessary (the faces of a connected planar subdivision are disjoint
by construction) and is skipped.

Inversion caveat
----------------
An inversion periodic point (``k_value == 2 * period``) has two branches of each
manifold, all four rays leaving ONE point, and
:meth:`~tanglepack.numerics.TangleWorkbench.TangleWorkbench._register_anchors`
registers one anchor per (unstable branch, stable branch) pair — four coincident
nodes where the topology wants a single degree-four node whose ``u-`` and ``s-``
slots are branch 1's ``u+`` and ``s+``. The arrangement is therefore not validated
on inversion fixtures; the region layer is built and tested on the non-inversion
ones.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Iterable, Literal, Optional, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from ..numerics.geometry import point_in_polygon
from .TopologyResults import Arc, Region, canonical_corners

if TYPE_CHECKING:
    from ..numerics.Bridge import BridgeId
    from .Trellis import Trellis

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def _is_cyclic_subsequence(
    sub: "list[int] | tuple[int, ...]", full: tuple[int, ...]
) -> bool:
    """Whether ``sub`` appears in ``full`` in order, under some rotation of ``full``."""
    if not sub:
        return True
    if len(sub) > len(full):
        return False
    for start in range(len(full)):
        rotated = full[start:] + full[:start]
        position = 0
        for item in rotated:
            if item == sub[position]:
                position += 1
                if position == len(sub):
                    return True
    return False


#: The four manifold rays leaving a crossing (see the module Dev Notes).
Slot = Literal["u+", "u-", "s+", "s-"]

_OPPOSITE_SLOT: dict[Slot, Slot] = {
    "u+": "u-",
    "u-": "u+",
    "s+": "s-",
    "s-": "s+",
}

#: Every slot, in the counter-clockwise order of the positive rotation.
_ALL_SLOTS: tuple[Slot, ...] = ("u+", "s+", "u-", "s-")

#: Counter-clockwise rotation at a node, by ``crossing_sign``.
_ROTATION: dict[int, tuple[Slot, ...]] = {
    1: _ALL_SLOTS,
    -1: ("u+", "s-", "u-", "s+"),
}


@dataclass
class _HalfEdge:
    """One directed side of one arc, or a dangling stub."""

    index: int
    tail: int
    head: int
    slot: Slot
    twin: int = -1
    arc: Optional[Arc] = None

    @property
    def is_virtual(self) -> bool:
        """True for a stub standing in for a manifold that has not been computed
        this far."""
        return self.arc is None


@dataclass
class _Node:
    """One node of the arrangement and the cyclic order of its rays."""

    id: int
    rotation: tuple[Slot, ...]
    slots: dict[Slot, int] = field(default_factory=dict)
    virtual: bool = False


class Arrangement:
    """
    The faces of a trellis, as a planar subdivision built from its crossings.

    Nodes are registry ids (anchors included — a periodic point is a crossing of
    its own two manifolds). Edges are the arcs between consecutive crossings: a
    stable arc for each consecutive pair on a stable branch, and a bridge for each
    consecutive pair on an unstable branch. Faces come from one linear traversal of
    the half-edges; the closed ones are :class:`~.TopologyResults.Region` objects.

    Build one with :meth:`from_trellis`, or take the session's cached one with
    :meth:`~tanglepack.loom.TangleSession.TangleSession.arrangement`.

    Attributes:
        trellis: The trellis this arrangement was built from. It resolves every id
            back to an :class:`~tanglepack.numerics.Intersection.Intersection`, a
            manifold, or a bridge.
        faces: Every face, in discovery order — open and non-minimal ones included.
        regions: The faces that are regions of the tangle: closed (no dangling
            end) AND minimal (enclosing no other component). See
            :attr:`open_faces` and :attr:`containing_faces` for the rest.
    """

    def __init__(self, trellis: "Trellis"):
        """
        Args:
            trellis: The trellis to build over. Prefer :meth:`from_trellis`, which
                also runs the build.
        """
        self.trellis = trellis
        self._nodes: dict[int, _Node] = {}
        self._half_edges: list[_HalfEdge] = []
        self._next_virtual_id: int = -1
        self.faces: list[Region] = []
        self.regions: list[Region] = []
        self._component_of: dict[int, int] = {}
        self._nodes_by_component: dict[int, list[int]] = {}
        self._region_by_corners: dict[tuple[int, ...], Region] = {}
        self._regions_at: dict[int, list[Region]] = {}
        self._regions_by_edge: dict[tuple, list[Region]] = {}
        self._face_of_half_edge: dict[int, Region] = {}

    # ── construction ────────────────────────────────────────────────────────

    @classmethod
    def from_trellis(cls, trellis: "Trellis") -> "Arrangement":
        """
        Build the arrangement of a trellis.

        Args:
            trellis: The trellis to build over. Pass the ALL-fixed-points trellis
                (``session.trellis(None)``) unless you deliberately want one
                tangle in isolation: a heteroclinic crossing belongs to two fixed
                points, and a single-fixed-point snapshot cuts the other side's
                arcs off, turning real faces into open ones.

        Returns:
            The built arrangement.

        Raises:
            AssertionError: If the unstable arcs read off the branch orderings are
                not exactly the trellis's bridges (see :meth:`_build_unstable_arcs`).
        """
        arrangement = cls(trellis)
        arrangement._build()
        return arrangement

    def _build(self) -> None:
        """Nodes, arcs, dangling stubs, faces, then which faces are regions."""
        self._build_nodes()
        self._build_stable_arcs()
        self._build_unstable_arcs()
        self._build_virtual_half_edges()
        self._build_faces()
        self._build_components()
        self._classify_minimality()
        self._index_regions()
        logger.debug(
            "Arrangement: %d nodes, %d half-edges, %d faces, %d regions "
            "(%d components)",
            len(self._nodes),
            len(self._half_edges),
            len(self.faces),
            len(self.regions),
            self.component_count,
        )

    def _build_nodes(self) -> None:
        """One node per crossing that lies on a branch of this trellis."""
        for branch in self.trellis.branches.values():
            for intersection_id in branch.intersection_ids:
                if intersection_id in self._nodes:
                    continue
                sign = self.trellis.intersection(intersection_id).crossing_sign
                if sign not in _ROTATION:
                    logger.warning(
                        "Crossing %d has no computed crossing_sign; assuming the "
                        "positive rotation (u+, s+, u-, s-)",
                        intersection_id,
                    )
                    sign = 1
                self._nodes[intersection_id] = _Node(
                    id=intersection_id, rotation=_ROTATION[sign]
                )

    def _build_stable_arcs(self) -> None:
        """One arc per consecutive pair on each stable branch."""
        for branch in self.trellis.stable_branches:
            ordered = branch.ordered_ids()
            for lo_id, hi_id in zip(ordered, ordered[1:]):
                self._add_edge(
                    Arc(
                        kind="stable",
                        lo_id=lo_id,
                        hi_id=hi_id,
                        branch_key=branch.key,
                        bridge_id=None,
                    )
                )

    def _build_unstable_arcs(self) -> None:
        """
        One arc per consecutive pair on each unstable branch; each is a bridge.

        "The arc between two consecutive crossings of an unstable branch" and "a
        bridge" are two names for one object, so the two constructions are
        cross-checked here rather than trusted. Only one direction of that check is
        an invariant, though:

        * **Every bridge is an arc.** A bridge is cut at two crossings of one
          unstable curve with nothing between them, so its endpoints have to be
          consecutive in that branch's ordering. A violation means the branch
          ordering and the bridge set describe different curves — asserted.
        * **Not every arc need be a bridge.** After a blast the branch also carries
          crossings born on ITERATED images of its bridges. The arcs between those
          were never cut as bridge objects (and a partial image piece has no id at
          all), so the pair is a real piece of unstable manifold with no Bridge
          behind it. That is expected, not a fault — logged, and the arc is built
          with ``bridge_id`` None so its geometry falls back to the parent branch
          rather than silently reading some other curve's points.

        Raises:
            AssertionError: If a bridge's endpoints are not consecutive on any
                unstable branch of this trellis.
        """
        pairs: set["BridgeId"] = set()
        cut: set["BridgeId"] = {
            bridge.id for bridge in self.trellis.bridges if bridge.id is not None
        }
        for branch in self.trellis.unstable_branches:
            ordered = branch.ordered_ids()
            for lo_id, hi_id in zip(ordered, ordered[1:]):
                pairs.add((lo_id, hi_id))
                self._add_edge(
                    Arc(
                        kind="unstable",
                        lo_id=lo_id,
                        hi_id=hi_id,
                        branch_key=branch.key,
                        bridge_id=(lo_id, hi_id) if (lo_id, hi_id) in cut else None,
                    )
                )

        assert not (cut - pairs), (
            "these bridges are not consecutive crossings of any unstable branch of "
            f"this trellis: {sorted(cut - pairs)} — the branch ordering and the "
            "bridge set describe different curves"
        )
        uncut = pairs - cut
        if uncut:
            logger.debug(
                "%d unstable arcs have no bridge object (iterated-image crossings "
                "or partial pieces): %s",
                len(uncut),
                sorted(uncut),
            )

    def _add_edge(self, arc: Arc) -> None:
        """Register both half-edges of one arc, in the two slots it occupies."""
        forward_slot: Slot = "s+" if arc.kind == "stable" else "u+"
        backward_slot: Slot = _OPPOSITE_SLOT[forward_slot]

        forward = self._new_half_edge(arc.lo_id, arc.hi_id, forward_slot, arc)
        backward = self._new_half_edge(
            arc.hi_id, arc.lo_id, backward_slot, arc.reversed()
        )
        forward.twin = backward.index
        backward.twin = forward.index

    def _new_half_edge(
        self, tail: int, head: int, slot: Slot, arc: Optional[Arc]
    ) -> _HalfEdge:
        """Allocate a half-edge and claim its slot at the tail node."""
        half_edge = _HalfEdge(
            index=len(self._half_edges), tail=tail, head=head, slot=slot, arc=arc
        )
        self._half_edges.append(half_edge)
        node = self._nodes[tail]
        assert slot not in node.slots, (
            f"node {tail} already has a {slot} ray; a crossing has exactly one "
            "unstable and one stable curve through it"
        )
        node.slots[slot] = half_edge.index
        return half_edge

    def _build_virtual_half_edges(self) -> None:
        """Fill every empty slot with a stub to a fresh degree-one node."""
        for node in list(self._nodes.values()):
            if node.virtual:
                continue
            for slot in _ALL_SLOTS:
                if slot in node.slots:
                    continue
                virtual_id = self._next_virtual_id
                self._next_virtual_id -= 1
                opposite = _OPPOSITE_SLOT[slot]
                self._nodes[virtual_id] = _Node(
                    id=virtual_id, rotation=(opposite,), virtual=True
                )
                out = self._new_half_edge(node.id, virtual_id, slot, None)
                back = self._new_half_edge(virtual_id, node.id, opposite, None)
                out.twin = back.index
                back.twin = out.index

    # ── face traversal ──────────────────────────────────────────────────────

    def _next(self, index: int) -> int:
        """The half-edge that follows ``index`` around its face."""
        twin = self._half_edges[self._half_edges[index].twin]
        node = self._nodes[twin.tail]
        rotation = node.rotation
        position = rotation.index(twin.slot)
        return node.slots[rotation[(position + 1) % len(rotation)]]

    def _build_faces(self) -> None:
        """Walk every ``next`` orbit once; each orbit is one face."""
        visited: set[int] = set()
        for start in range(len(self._half_edges)):
            if start in visited:
                continue
            cycle: list[int] = []
            current = start
            while current not in visited:
                visited.add(current)
                cycle.append(current)
                current = self._next(current)
            assert current == start, (
                "face traversal did not close: next() must be a permutation of "
                "the half-edges"
            )
            self._record_face(cycle)

    def _record_face(self, cycle: list[int]) -> None:
        """Turn one half-edge cycle into a Region and index its half-edges."""
        edges = [self._half_edges[index] for index in cycle]
        is_closed = not any(edge.is_virtual for edge in edges)

        raw_corners = [edge.tail for edge in edges]
        arcs = [edge.arc for edge in edges]
        if is_closed:
            # Rotate corners AND arcs together by the shift that canonicalises the
            # corner cycle, so arcs[i] still runs corners[i] -> corners[i + 1].
            corners, shift = min(
                (tuple(raw_corners[i:] + raw_corners[:i]), i)
                for i in range(len(raw_corners))
            )
            arcs = arcs[shift:] + arcs[:shift]
            assert corners == canonical_corners(raw_corners)
        else:
            corners = tuple(raw_corners)

        face = Region(
            corners=corners,
            arcs=[arc for arc in arcs if arc is not None],
            arrangement=self,
            is_closed=is_closed,
        )
        self.faces.append(face)
        for index in cycle:
            self._face_of_half_edge[index] = face

    # ── components and minimality ───────────────────────────────────────────

    def _build_components(self) -> None:
        """Union-find over the real nodes, joined by every real arc.

        A trellis need not be connected: two fixed points whose manifolds have no
        COMPUTED heteroclinic crossing between them are two separate drawings on
        one plane, and the face traversal — which only ever walks along arcs —
        cannot see one from the other.
        """
        parent: dict[int, int] = {
            node_id: node_id for node_id, node in self._nodes.items() if not node.virtual
        }

        def find(item: int) -> int:
            while parent[item] != item:
                parent[item] = parent[parent[item]]
                item = parent[item]
            return item

        for edge in self._half_edges:
            if edge.arc is None:
                continue
            a, b = find(edge.tail), find(edge.head)
            if a != b:
                parent[a] = b

        self._component_of = {node_id: find(node_id) for node_id in parent}
        self._nodes_by_component: dict[int, list[int]] = {}
        for node_id, root in self._component_of.items():
            self._nodes_by_component.setdefault(root, []).append(node_id)

    def _classify_minimality(self) -> None:
        """
        Mark every closed face that encloses a piece of ANOTHER component.

        A region of a tangle "contains no other manifold piece" (the plan's
        definition). Within one connected component that is automatic: the faces of
        a connected planar subdivision are interior-disjoint by construction. Across
        components it is NOT, and the failure is the big one — the outer face of the
        period-1 tangle in the nested period-3 fixture is a perfectly good closed
        cycle of four crossings that geometrically swallows the whole period-3
        tangle and 25 of its faces.

        So the containment test is run exactly where it can fail: for each closed
        face, against the nodes of the other components. One node inside is enough —
        a component is connected, so if any of its crossings is inside the face and
        none of its arcs crosses the face's boundary (they cannot; boundaries are
        manifold and manifolds do not cross except at crossings, which would have
        joined the components), all of it is.

        Faces that fail are kept, on :attr:`containing_faces`, with
        :attr:`~.TopologyResults.Region.is_minimal` False; they are still real
        cycles of the tangle, just not regions.
        """
        if len(self._nodes_by_component) <= 1:
            return
        coords = {
            node_id: self.trellis.intersection(node_id).get_point()
            for node_id in self._component_of
        }
        for face in self.faces:
            if not face.is_closed or not face.corners:
                continue
            own = self._component_of[face.corners[0]]
            polygon = face.boundary_points
            if len(polygon) < 3:
                continue
            for component, node_ids in self._nodes_by_component.items():
                if component == own:
                    continue
                if any(
                    point_in_polygon(coords[node_id], polygon) for node_id in node_ids
                ):
                    face.is_minimal = False
                    logger.debug(
                        "Face %s encloses another component's crossings; it is a "
                        "cycle of the tangle but not a region",
                        face.corners,
                    )
                    break

    def _index_regions(self) -> None:
        """Build the corner / node / edge lookups over the closed minimal faces."""
        for region in self.faces:
            if not (region.is_closed and region.is_minimal):
                continue
            self.regions.append(region)
            if region.corners in self._region_by_corners:
                logger.warning(
                    "Two closed faces share the corner cycle %s; keeping the first",
                    region.corners,
                )
            else:
                self._region_by_corners[region.corners] = region
            for corner in set(region.corners):
                self._regions_at.setdefault(corner, []).append(region)
            for arc in region.arcs:
                self._regions_by_edge.setdefault(arc.edge_key, []).append(region)

    # ── combinatorial lookups ───────────────────────────────────────────────

    def region(self, corners: "Iterable[int]") -> Optional[Region]:
        """
        The closed region with this corner cycle.

        Args:
            corners: The corner ids in cyclic order, in either direction and
                starting anywhere — the lookup canonicalises them.

        Returns:
            The region, or None if no closed face has that cycle.
        """
        items = tuple(corners)
        found = self._region_by_corners.get(canonical_corners(items))
        if found is not None:
            return found
        # A region's cycle is directed; an image under an orientation-REVERSING
        # map arrives wound the other way, so the reverse cycle names it too.
        return self._region_by_corners.get(canonical_corners(items[::-1]))

    def regions_at(self, intersection_id: int) -> list[Region]:
        """
        Every closed region having this crossing as a corner.

        Args:
            intersection_id: Registry id of the crossing.

        Returns:
            The regions, in discovery order (empty if the crossing is not a corner
            of any closed face).
        """
        return list(self._regions_at.get(intersection_id, []))

    def regions_bounded_by(self, arc: Arc) -> list[Region]:
        """
        The closed regions on either side of one arc.

        Direction is ignored: an arc and its reverse bound the same two faces.

        Args:
            arc: The arc, in either direction.

        Returns:
            Up to two regions; fewer when the face on a side is open.
        """
        return list(self._regions_by_edge.get(arc.edge_key, []))

    def face_at_stub(self, intersection_id: int, slot: Slot) -> Region:
        """
        The face a dangling manifold end sticks into.

        A slot with no computed arc carries a VIRTUAL half-edge to a degree-one
        node, and the traversal reflects off it: both the stub and its twin
        belong to the same face cycle. That face is therefore the one piece of
        plane the uncomputed continuation of the manifold starts in — which is
        exactly what a walk needs when it has to leave the computed picture at
        the outermost crossing of a branch.

        Args:
            intersection_id: Registry id of the crossing the ray leaves from.
            slot: Which of the four rays (``"u+"``, ``"u-"``, ``"s+"``,
                ``"s-"``) to follow.

        Returns:
            The (necessarily open) face containing that stub.

        Raises:
            ValueError: If the crossing is not a node of this arrangement, if
                the slot is not one of the four rays, or if the slot carries a
                real arc rather than a stub — a computed arc runs to another
                crossing and bounds two faces, so "the face at the stub" has no
                meaning there.
        """
        if slot not in _ALL_SLOTS:
            raise ValueError(
                f"slot must be one of {_ALL_SLOTS}, got {slot!r}"
            )
        node = self._nodes.get(intersection_id)
        if node is None:
            raise ValueError(
                f"crossing {intersection_id} is not a node of this arrangement"
            )
        index = node.slots.get(slot)
        if index is None:
            raise ValueError(
                f"node {intersection_id} has no {slot} ray at all; the "
                "arrangement was not fully built"
            )
        half_edge = self._half_edges[index]
        if not half_edge.is_virtual:
            raise ValueError(
                f"the {slot} ray of crossing {intersection_id} runs to crossing "
                f"{half_edge.head} along a computed arc, so it is not a "
                "dangling stub and bounds two faces rather than sticking into one"
            )
        face = self._face_of_half_edge.get(index)
        if face is None:  # pragma: no cover - the traversal covers every edge
            raise ValueError(
                f"the {slot} stub of crossing {intersection_id} belongs to no "
                "face; the face traversal did not visit every half-edge"
            )
        return face

    def image_of(
        self, region: Region, n: int = 1, *, rtol: float = 2e-2
    ) -> Optional[Region]:
        """
        The region a closed region maps to under ``n`` steps of the map.

        The map is a homeomorphism, so it carries crossings to crossings and the
        face they bound to the face their images bound: the corners go through the
        registry's iterate table and the resulting cycle is looked up.

        The combinatorics alone are not conclusive, though, and the way they fail
        is subtle. The image face often carries MORE corners than the source,
        because a trimmed stable manifold hides crossings and a crossing on the
        image boundary can have no detected preimage on the source boundary. But
        an extra arc between two mapped corners may equally be a chord across the
        image's INTERIOR, splitting it into several faces, and then the face
        carrying the mapped cycle is a proper SUB-face of the image, not the image.
        On the k=10 fixture that is exactly what happens: region ``(3, 7)`` of area
        43.02 maps onto ``(1, 4, 5, 6)`` of area 36.96, whose union with the little
        face ``(5, 6)`` recovers the 43.02.

        These maps are AREA-PRESERVING, which settles it without any guesswork: a
        candidate is accepted only if its area matches the source's (to ``rtol``,
        which is a bound on polyline resolution, not on the dynamics) AND the
        source's representative point, mapped forward for real, lands in it. A
        sub-face fails the first test; a same-area face elsewhere fails the second.

        Args:
            region: A closed region of this arrangement.
            n: Number of forward map steps (negative for backward — the iterate
               table stores both directions).
            rtol: Relative tolerance on the area match.

        Returns:
            The image region, or None when a corner has no recorded ``n``-iterate,
            when the face its images bound is open, when more than one closed face
            carries the mapped cycle, or when no single face passes the area and
            containment tests.

        Note:
            An image split across SEVERAL faces returns None rather than a union.
            The split is an artifact of what has been computed (the crossings that
            subdivide it have no preimage in the picture), so what to do with it —
            merge the faces, grow until the preimages appear, or carry a union as a
            first-class object — is a dual-graph decision, and the dual graph is
            deliberately out of scope here.
        """
        images: list[int] = []
        for corner in region.corners:
            image = self.trellis.iterate(corner, n)
            if image is None:
                return None
            images.append(image)

        candidates: list[Region] = []
        exact = self.region(images)
        if exact is not None:
            candidates.append(exact)
        else:
            candidates = [
                face
                for face in self.regions_at(images[0])
                if all(corner in face.corners for corner in images)
                and (
                    _is_cyclic_subsequence(images, face.corners)
                    or _is_cyclic_subsequence(images[::-1], face.corners)
                )
            ]
        if not candidates:
            return None
        if len(candidates) > 1:
            logger.debug(
                "Corner cycle %s is carried by %d closed faces; refusing to guess",
                images,
                len(candidates),
            )
            return None

        return self._verify_image(region, candidates[0], n, rtol=rtol)

    def _verify_image(
        self, region: Region, candidate: Region, n: int, *, rtol: float
    ) -> Optional[Region]:
        """Accept a candidate image only if the area and the dynamics agree."""
        source_area = abs(region.area)
        candidate_area = abs(candidate.area)
        if source_area <= 0.0 or abs(candidate_area - source_area) > rtol * source_area:
            logger.debug(
                "Face %s is not the image of %s: area %.6g vs %.6g (the map is "
                "area preserving, so the image is subdivided or elsewhere)",
                candidate.corners,
                region.corners,
                candidate_area,
                source_area,
            )
            return None

        point = region.representative_point
        if point is None:
            return None
        mapped = self._map_point(point, n)
        if not candidate.contains(mapped):
            logger.debug(
                "Face %s has the right area but %s's representative point does not "
                "map into it",
                candidate.corners,
                region.corners,
            )
            return None
        return candidate

    def _map_point(
        self, point: "NDArray[np.float64] | tuple[float, float]", n: int
    ) -> NDArray[np.float64]:
        """
        Apply the real map ``n`` times (its inverse when ``n`` is negative).

        Args:
            point: The ``(x, y)`` to carry forward (or backward).
            n: Number of map steps; negative uses the inverse map.

        Returns:
            The ``(2,)`` image point.

        Raises:
            ValueError: If the trellis carries no dynamical system.
        """
        system = self.trellis.dynamical_system
        if system is None:
            raise ValueError(
                "this arrangement's trellis carries no dynamical system, so an "
                "image cannot be verified against the map"
            )
        step = system.map if n >= 0 else system.map_inv
        current = np.asarray(point, dtype=float)
        for _ in range(abs(n)):
            current = np.asarray(step(current), dtype=float)
        return current

    def preimage_of(
        self, region: Region, n: int = 1, *, rtol: float = 2e-2
    ) -> Optional[Region]:
        """
        The region that maps ONTO a closed region in ``n`` steps.

        Args:
            region: A closed region of this arrangement.
            n: Number of map steps to go back.
            rtol: Relative tolerance on the area match.

        Returns:
            The preimage region, or None (see :meth:`image_of`).
        """
        return self.image_of(region, -n, rtol=rtol)

    # ── diagnostics ─────────────────────────────────────────────────────────

    @property
    def open_faces(self) -> list[Region]:
        """The faces with a dangling end — bounded partly by "not computed yet"."""
        return [face for face in self.faces if not face.is_closed]

    @property
    def containing_faces(self) -> list[Region]:
        """
        The closed faces that are not regions because they enclose another one.

        A closed cycle of arcs that geometrically swallows a whole other connected
        component of the trellis (see :meth:`_classify_minimality`) contains other
        manifold pieces, so by the definition of a region it is not one. It is kept
        here rather than discarded: it is a real cycle of the tangle, and the fact
        that one tangle nests inside a face of another is exactly what a nested
        fixture is about.
        """
        return [
            face for face in self.faces if face.is_closed and not face.is_minimal
        ]

    @property
    def component_count(self) -> int:
        """
        How many connected components the trellis's arcs form.

        More than one means two tangles have no COMPUTED crossing joining them —
        which is why :attr:`regions` needs the containment test at all.
        """
        return len(self._nodes_by_component)

    def component_of(self, intersection_id: int) -> Optional[int]:
        """
        A representative node id of the component containing this crossing.

        Args:
            intersection_id: Registry id of a crossing in this arrangement.

        Returns:
            The component's representative id, or None for an unknown crossing.
            The value is a label only; compare two of them for equality.
        """
        return self._component_of.get(intersection_id)

    @property
    def euler_characteristic(self) -> int:
        """
        ``V - E + F`` over the arrangement, virtual nodes and stubs included.

        Equals ``2 * component_count``, and that is a cheap end-to-end check that
        the face traversal found every face exactly once. The usual planar formula
        for a graph with C components is ``V - E + F = 1 + C``, which counts ONE
        shared unbounded face; this traversal instead walks each component's outer
        boundary as its own cycle, so it reports C outer faces rather than one and
        each component contributes its own ``V - E + F = 2``.
        """
        return (
            len(self._nodes) - len(self._half_edges) // 2 + len(self.faces)
        )

    def summary(self) -> str:
        """A one-line description of the arrangement's size."""
        real = sum(1 for node in self._nodes.values() if not node.virtual)
        virtual = len(self._nodes) - real
        return (
            f"Arrangement: {real} crossings (+{virtual} dangling ends), "
            f"{self.component_count} component(s), {len(self.faces)} faces, "
            f"{len(self.regions)} regions "
            f"(+{len(self.containing_faces)} containing, "
            f"{len(self.open_faces)} open)"
        )

    def __repr__(self) -> str:
        return f"<{self.summary()}>"
