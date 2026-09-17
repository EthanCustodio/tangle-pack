"""
The result schema of the topological layer.

Every dataclass the topological algorithms hand back --
:class:`PseudoneighborPair`, :class:`Hole`, :class:`PartitionInterval`,
:class:`Arc` and friends -- plus the :data:`Side` / :data:`Endpoint` labels
they are annotated with and :class:`ElementRef`, the global name of a partition
element, so that the algorithms, :class:`~.Trellis.Trellis` and the session
facade all name the same things.

Dev Notes:

These dataclasses are pure *result containers*. They hold the output of the
topological algorithms (Compute-Pseudoneighbors, Is-Strong-Pip, ...) so that a
single Trellis object can carry both the input trellis and every derived
classification. No algorithm logic lives here — only the schema of what the
algorithms produce.

Intersections are referenced by their integer registry ID rather than by object
so that a result survives a registry rebuild and is cheap to serialise.

Open question: a Hole is currently a single phase-space coordinate. We may
eventually want a Hole to also reference the bounding bridge(s) or the
enclosed region directly, rather than just via ``bounding_ids``.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Literal, Optional, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from ..numerics.geometry import (
    arc_polyline,
    point_in_polygon,
    polygon_interior_point,
    polyline_midpoint,
    signed_polygon_area,
)

if TYPE_CHECKING:
    from ..numerics.Bridge import BridgeId
    from ..numerics.FixedPoint import FixedPoint
    from ..numerics.Intersection import ManifoldKey
    from .Arrangement import Arrangement
    from .Trellis import Trellis

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# The two orientation labels the topological layer uses, defined here (with the
# result schema they annotate) so that StablePartition, Trellis, and the session
# facade all name the same thing.
Side = Literal["left", "right"]
Endpoint = Literal["first", "second"]

OPPOSITE_SIDE: dict[Side, Side] = {"left": "right", "right": "left"}


def endpoint_index(endpoint: Endpoint) -> int:
    """
    Position of a named bridge endpoint within its BridgeId.

    Args:
        endpoint: ``"first"`` or ``"second"``.

    Returns:
        0 for ``"first"``, 1 for ``"second"``.

    Raises:
        ValueError: If ``endpoint`` names neither end.
    """
    if endpoint == "first":
        return 0
    if endpoint == "second":
        return 1
    raise ValueError(f"endpoint must be 'first' or 'second', got {endpoint!r}")


@dataclass
class Hole:
    """
    A puncture placed in the region bounded by a pseudoneighbor pair.

    Per the Pseudoneighbor Algorithm, for each pseudoneighbor pair a hole is
    punched in the bounded region, placed infinitesimally close to one of the
    two pseudoneighbors.

    Attributes:
        coords: Phase-space (x, y) where the hole is placed.
        near_intersection_id: Registry ID of the pseudoneighbor the hole hugs.
        pair: The pseudoneighbor pair this hole belongs to (back-reference).
        bridge_side: Which side of the hole's own bridge the hole sits on,
            with the bridge oriented by the unstable dynamical direction —
            forward flow away from the fixed point, i.e. increasing unstable
            canonical distance. Standing on the bridge looking along that
            direction, positive cross(tangent, displacement) = left. Classified
            once at punch time: a direct hole from its own coordinates, a
            propagated hole from the backward-carried point. None when the
            geometry is degenerate (unorientable bridge or a point on the arc).
            Invariant: all holes sharing an ``origin`` share this side, since
            an orientation-preserving map carries a bridge's dynamical
            orientation to its image's (it alternates with the parity of
            ``iterate`` when the map reverses orientation). Checked by
            :func:`topology.StablePartition.check_holes_share_bridge_side`.
        bounding_ids: Registry IDs of the two intersections defining the bridge
            the hole belongs to — the pair's own bridge for a directly punched
            hole, the containing bridge for a propagated one. These are the
            partition boundaries of the Stable Manifold Partition Algorithm.
        openings: The partition intervals this hole opens, one record per
            defining intersection: ``(intersection_id, which, row)`` where
            ``which`` is ``"anchorward"``/``"outward"`` (which side of that
            intersection along the stable manifold the hole abuts — derived
            from ``bridge_side``) and ``row`` is the left/right partition the
            opening acts on. A hole with no openings does not participate in
            the partition.
        iterate: Position of this hole along its reference orbit — 0 for the
            reference hole itself, negative for backward iterates, positive
            for forward ones. None when unknown.
        origin: ``as_tuple()`` of the reference pseudoneighbor pair whose orbit
            this hole belongs to (holes with the same origin share a plot
            marker). None when unknown.
    """

    coords: tuple[float, float]
    near_intersection_id: int
    pair: Optional["PseudoneighborPair"] = None
    bridge_side: Optional[Side] = None
    bounding_ids: Optional[tuple[int, int]] = None
    openings: Optional[list[tuple[int, str, str]]] = None
    iterate: Optional[int] = None
    origin: Optional[tuple[int, int]] = None


@dataclass
class PseudoneighborPair:
    """
    A pair of intersection points classified as pseudoneighbors.

    Two intersections x, x' are a pseudoneighbor pair when the open stable and
    unstable intervals connecting them contain no trellis intersection nor any
    iterate of a trellis intersection (see Pseudoneighbor_Algorithm.pdf).

    Attributes:
        intersection_a: Registry ID of the first pseudoneighbor.
        intersection_b: Registry ID of the second pseudoneighbor.
        branch_key: Manifold key of the stable branch on which the pair was found.
        hole: The hole punched for this pair, if one has been placed.
        is_reference: True for a pair found on the reference window
            W^S(r_n, r_{n+p}) by Compute-Pseudoneighbors; False for a pair
            generated from a reference by iterating (see
            extend_pseudoneighbor_trajectories).
        iterate: Steps along the orbit from the reference pair — 0 for the
            reference itself, negative for backward iterates, positive for
            forward ones. None when unknown.
        origin: ``as_tuple()`` of the reference pair this pair was generated
            from (a reference is its own origin). None when unknown.
    """

    intersection_a: int
    intersection_b: int
    branch_key: Optional["ManifoldKey"] = None
    hole: Optional[Hole] = None
    is_reference: bool = False
    iterate: Optional[int] = None
    origin: Optional[tuple[int, int]] = None

    def as_tuple(self) -> tuple[int, int]:
        """Return the unordered ID pair as a sorted tuple (for set membership)."""
        return tuple(sorted((self.intersection_a, self.intersection_b)))


@dataclass
class StrongPipResult:
    """
    Classification of a single intersection as a strong pip (or not).

    A point q0 is a strong pip when no trellis intersection, mapped back onto
    the stable branch attached to z'_0, has a canonical distance less than that
    of q0 (see Strong_Pip_Algorithm.pdf).

    Attributes:
        intersection_id: Registry ID of the classified intersection.
        is_strong_pip: True if q0 is a strong pip.
        blocking_intersection_id: If not a strong pip, the registry ID of an
            intersection whose mapped-back canonical distance disqualifies q0.
            None when q0 is a strong pip.
    """

    intersection_id: int
    is_strong_pip: bool
    blocking_intersection_id: Optional[int] = None


@dataclass(frozen=True)
class ElementRef:
    """
    The global identity of one stable-partition element.

    An ``element_id`` alone only names an element *within* one
    :class:`StablePartitionResult`; the branch it lies on and the side whose
    holes cut it are what make the name global. This triple is that name — and
    being frozen it hashes, so it is usable as a dict key and as half of a
    :class:`~tanglepack.topology.BridgeClass.BridgeClass`.

    Equality and hashing are the dataclass defaults, which means the branch key
    is compared componentwise and its
    :class:`~tanglepack.numerics.FixedPoint.FixedPoint` by identity — two refs
    are equal exactly when they name the same element of the same live orbit.

    Attributes:
        branch_key: Manifold key of the stable branch the element lies on.
        side: Which side's partition the element belongs to.
        element_id: Index of the element within its result's ``intervals``,
            counted from the anchor outward.
    """

    branch_key: "ManifoldKey"
    side: Side
    element_id: int

    @property
    def fixed_point(self) -> "FixedPoint":
        """The periodic point the element's branch is anchored at."""
        return self.branch_key[0]

    @property
    def orbit_index(self) -> int:
        """Index of the branch's anchor within its periodic orbit."""
        return self.branch_key[2]

    @property
    def branch_index(self) -> int:
        """Index of the branch at that anchor (0, or 1 on an inversion point)."""
        return self.branch_key[3]

    @property
    def label(self) -> str:
        """
        A short, deterministic, address-free name for this element.

        Formatted ``p{period}@{orbit}.{branch}/{L|R}#{element_id}``, e.g.
        ``p3@1.0/L#2``. Two runs of the same trellis produce the same label, so
        it is safe to put in a report, a plot legend or a test expectation —
        unlike the default ``repr``, which prints the FixedPoint's address.

        Returns:
            The label string.
        """
        return (
            f"p{self.fixed_point.period}@{self.orbit_index}.{self.branch_index}"
            f"/{self.side[0].upper()}#{self.element_id}"
        )

    def __str__(self) -> str:
        return self.label


@dataclass
class PartitionInterval:
    """
    One interval of a stable-manifold partition, in canonical-distance units.

    Intervals run along a single stable branch from the anchor point outward.
    An interval bounding a punched hole is open at that end (the bounding
    intersection is excluded); all other ends are closed. A point bounding
    holes on both of its adjacent intervals appears as a degenerate closed
    singleton [x, x].

    An interval is also a partition *element*: the identifiable piece of stable
    branch the region layer will glue a face to. Its full identity is
    ``(branch_key, side, element_id)`` — the id alone only names it within one
    :class:`StablePartitionResult`.

    Attributes:
        lo_id: Registry ID of the lower (toward-anchor) boundary intersection,
            or None when the interval starts at the anchor point (cdist 0).
        hi_id: Registry ID of the upper (outward) boundary intersection, or
            None when the interval runs to the end of the computed branch.
        lo_cdist: Stable canonical distance of the lower boundary.
        hi_cdist: Stable canonical distance of the upper boundary.
        closed_lo: True if the lower boundary point belongs to the interval.
        closed_hi: True if the upper boundary point belongs to the interval.
        element_id: Index of this element within its result's ``intervals``
            (anchor outward). Stamped by
            :func:`topology.StablePartition.partition_stable_manifold`; -1 on an
            interval that has not been through it.
        branch_key: Manifold key of the stable branch this element lies on.
            None until stamped.
        side: Which side's partition this element belongs to. None until stamped.
        parent_element_id: For an element of a REFINED partition (see
            :class:`~tanglepack.topology.PartitionFamily.IteratedHomotopyPartition`),
            the ``element_id`` of the element of the partition it refines
            (same branch, same side) that owns this element's midpoint. None
            on an element of an unrefined partition.
        cut_by: For an element of a refined partition, the
            :data:`~tanglepack.numerics.Bridge.BridgeId` of the image bridge
            this element lies UNDER: the bridge lies on this side and the
            element sits inside the closed span between its two crossings.
            The innermost such bridge when several cover the element. None
            for an element under no image bridge and on an unrefined
            partition.
    """

    lo_id: Optional[int]
    hi_id: Optional[int]
    lo_cdist: float
    hi_cdist: float
    closed_lo: bool
    closed_hi: bool
    element_id: int = -1
    branch_key: Optional["ManifoldKey"] = None
    side: Optional[Side] = None
    parent_element_id: Optional[int] = None
    cut_by: Optional["BridgeId"] = None

    @property
    def ref(self) -> ElementRef:
        """
        This element's global identity.

        Returns:
            The :class:`ElementRef` naming this element.

        Raises:
            ValueError: If the interval has not been stamped by
                :func:`topology.StablePartition.partition_stable_manifold` —
                a raw interval carries no branch, no side and ``element_id``
                -1, so it has no global name to give.
        """
        if self.element_id < 0 or self.branch_key is None or self.side is None:
            raise ValueError(
                "this PartitionInterval has not been stamped with its identity "
                f"(element_id={self.element_id}, branch_key={self.branch_key}, "
                f"side={self.side}); run partition_stable_manifold first"
            )
        return ElementRef(self.branch_key, self.side, self.element_id)


@dataclass
class StablePartitionResult:
    """
    The partition of one stable branch induced by the holes on one side.

    The stable manifold is oriented by its dynamical direction — forward flow
    toward the fixed point, i.e. looking toward the anchor — and each of its
    two sides (left/right, positive cross = left) induces its own partition
    (see Stable_Manifold_Partition_Algorithm.pdf). An opening acts on the row
    named in the hole's ``openings`` records; which endpoints are open is
    derived from the hole's side of its bridge.

    The intervals double as the branch's partition *elements*; the two lookup
    tables built alongside them answer the two questions the region layer asks:
    "which element owns this crossing" and "which elements sit at this bridge's
    two ends".

    Attributes:
        branch_key: Manifold key of the partitioned stable branch.
        side: Which side's holes this partition is built from.
        intervals: The partition intervals, ordered from the anchor outward.
        element_of_intersection: Every crossing on this branch mapped to the
            single element that owns it — the element whose span contains it,
            with a shared endpoint owned by whichever neighbour is closed there
            (a pinched singleton owns its own point). Keyed by registry ID, so
            the table survives everything but a renumbering of the registry.
        elements_at_bridge: Convenience index over the bridges the building
            trellis actually held: every non-partial one of them with at least
            one endpoint on this branch, mapped to ``(element at first, element
            at second)`` in :data:`~tanglepack.numerics.Bridge.BridgeId` order,
            with ``None`` for an end that lies on a different stable branch.
            It is NOT the authority on bridge-to-element: a trellis holds a
            bridge under its UNSTABLE fixed point, so a heteroclinic bridge of
            W^u(fp1) whose endpoints sit on W^s(fp3) is absent from fp3's
            snapshot entirely and has no entry here even though fp3's partition
            owns both its endpoints. Because a BridgeId IS the endpoint id pair,
            :meth:`Trellis.element_for` falls back to
            ``element_of_intersection`` and answers anyway — go through it, or
            through
            :meth:`~tanglepack.loom.TangleSession.TangleSession.partition_element_for`
            for the cross-trellis case, rather than reading this table directly.

    Note:
        Both tables are keyed by registry ids and element indices of *this*
        result, so they are only meaningful for the trellis snapshot that built
        them. A trellis is dropped and rebuilt whenever the workbench generation
        moves (which is what a registry renumbering does), and its results go
        with it — nothing migrates ids across a reindex.
    """

    branch_key: "ManifoldKey"
    side: Side
    intervals: list[PartitionInterval] = field(default_factory=list)
    element_of_intersection: dict[int, int] = field(default_factory=dict)
    elements_at_bridge: dict["BridgeId", tuple[Optional[int], Optional[int]]] = field(
        default_factory=dict
    )

    def element(self, element_id: int) -> PartitionInterval:
        """
        Return the element with this id.

        Args:
            element_id: Index of the element within :attr:`intervals`.

        Returns:
            The PartitionInterval carrying that id.

        Raises:
            IndexError: If ``element_id`` is not an element of this result.
        """
        if not 0 <= element_id < len(self.intervals):
            raise IndexError(
                f"element {element_id} is not one of this partition's "
                f"{len(self.intervals)} elements"
            )
        return self.intervals[element_id]

    def ref(self, element_id: int) -> ElementRef:
        """
        The global identity of one of this result's elements.

        Args:
            element_id: Index of the element within :attr:`intervals`.

        Returns:
            An :class:`ElementRef` built from this result's own branch and side.

        Raises:
            IndexError: If ``element_id`` is not an element of this result.

        Note:
            Built from the result's ``branch_key``/``side`` rather than from the
            interval's own stamp, so it answers even for a result whose
            intervals were assembled by hand; the two agree on anything
            :func:`topology.StablePartition.partition_stable_manifold` produced.
        """
        self.element(element_id)
        return ElementRef(self.branch_key, self.side, element_id)


# --------------------------------------------------------------------------- #
# The region layer: the pieces of manifold a face is bounded by, and the face.
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Arc:
    """
    One edge of the arrangement: a piece of manifold between two crossings.

    An arc is the geometric edge; whatever labels the algorithms hang on it (the
    partition elements on a stable arc, the :data:`BridgeId` of an unstable one)
    are read off it, never stored twice. It is directed: ``reverse`` says whether
    a face traverses it from ``hi_id`` down to ``lo_id`` instead of the natural
    anchor-outward direction. Two arcs that differ only in ``reverse`` are the
    two half-edges of one edge and share an :attr:`edge_key`.

    Attributes:
        kind: ``"stable"`` or ``"unstable"`` — which manifold family the arc runs
            along. Per the fundamental invariant these are the only two kinds of
            manifold piece a face boundary can be made of.
        lo_id: Registry id of the crossing at the LOWER canonical distance.
        hi_id: Registry id of the crossing at the higher canonical distance.
        branch_key: Manifold key of the branch this arc lies on. Both crossings
            lie on it, and canonical distances are only comparable within it.
        bridge_id: For an unstable arc, the :data:`~tanglepack.numerics.Bridge.BridgeId`
            of the bridge it IS — an unstable arc and a bridge are the same
            object seen from the topological and the numerical side. ``None`` for
            a stable arc.
        reverse: True when the arc is traversed from ``hi_id`` to ``lo_id``.
    """

    kind: Literal["stable", "unstable"]
    lo_id: int
    hi_id: int
    branch_key: "ManifoldKey"
    bridge_id: Optional["BridgeId"] = None
    reverse: bool = False

    @property
    def edge_key(self) -> tuple:
        """The undirected identity of this edge: kind, endpoints, branch."""
        return (self.kind, self.lo_id, self.hi_id, self.branch_key)

    @property
    def tail_id(self) -> int:
        """The crossing this arc is traversed FROM."""
        return self.hi_id if self.reverse else self.lo_id

    @property
    def head_id(self) -> int:
        """The crossing this arc is traversed TO."""
        return self.lo_id if self.reverse else self.hi_id

    def reversed(self) -> "Arc":
        """The same edge traversed the other way."""
        return Arc(
            kind=self.kind,
            lo_id=self.lo_id,
            hi_id=self.hi_id,
            branch_key=self.branch_key,
            bridge_id=self.bridge_id,
            reverse=not self.reverse,
        )

    # ── geometry (lazy: nothing is built unless asked for) ──────────────────

    def polyline(self, trellis: "Trellis") -> NDArray[np.float64]:
        """
        The arc as an ``(N, 2)`` polyline, closed at its two crossings.

        The crossings themselves are not manifold nodes, so the real nodes strictly
        between them are taken from the underlying curve and the two exact crossing
        coordinates are put on the ends. That makes consecutive arcs of a face join
        exactly, with no gap and no overshoot past a corner.

        Args:
            trellis: The trellis resolving this arc's crossings, branch manifold
                and (for an unstable arc) bridge.

        Returns:
            The polyline in traversal order — ``lo`` to ``hi``, or the reverse when
            :attr:`reverse`. At minimum the two endpoint coordinates.
        """
        lo = trellis.intersection(self.lo_id)
        hi = trellis.intersection(self.hi_id)
        interior = self._interior_polyline(trellis, lo, hi)
        parts = [lo.get_point().reshape(1, 2)]
        if len(interior):
            parts.append(interior)
        parts.append(hi.get_point().reshape(1, 2))
        poly = np.vstack(parts)
        return poly[::-1] if self.reverse else poly

    def _interior_polyline(self, trellis, lo, hi) -> NDArray[np.float64]:
        """The real manifold nodes strictly inside this arc, in lo-to-hi order."""
        tol = trellis.registry.cdist_tol
        if self.kind == "stable":
            manifold = trellis.manifolds.get(self.branch_key)
            if manifold is None:
                return np.empty((0, 2))
            return arc_polyline(
                manifold,
                lo.stable_cdist,
                hi.stable_cdist,
                stability="stable",
                tol=-tol,
            )
        # A bridge's root and tail sit just PAST its crossings, so whichever curve
        # is used the polyline is clipped back to the arc's span rather than taken
        # whole -- and arc_polyline normalizes the storage direction, which is what
        # oriented_bridge_polyline does for the partition's own callers.
        curve = trellis.bridge_between(self.lo_id, self.hi_id)
        if curve is None and abs(lo.unstable_cdist) <= tol:
            # No bridge was cut here, but the arc starts AT the periodic point, so
            # it runs along the parent branch from its root -- which is exactly a
            # resonance-zone boundary arc (periodic point out to its pip, possibly
            # spanning several bridges). Anchored-ness is what makes the fallback
            # safe: an arc that starts at cdist 0 is on the parent curve by
            # construction, whereas an arc between two crossings born on an
            # ITERATED image is not, and reading the parent's nodes over its cdist
            # span would return an unrelated piece of curve.
            curve = trellis.manifolds.get(self.branch_key)
        if curve is None:
            # Nothing to walk: the polyline degrades to the chord between the two
            # crossings, which is all the caller's endpoints already give it.
            return np.empty((0, 2))
        return arc_polyline(
            curve,
            lo.unstable_cdist,
            hi.unstable_cdist,
            stability="unstable",
            tol=-tol,
        )

    def midpoint(self, trellis: "Trellis") -> Optional[NDArray[np.float64]]:
        """The middle node of this arc's polyline (see
        :func:`~tanglepack.numerics.geometry.polyline_midpoint`)."""
        return polyline_midpoint(self.polyline(trellis))

    # ── partition labels ────────────────────────────────────────────────────

    def elements(self, trellis: "Trellis") -> list["PartitionInterval"]:
        """
        The stable-partition elements lying on this arc.

        An arc runs between two CONSECUTIVE crossings of its branch, and partition
        boundaries are crossings, so an element either covers the whole arc or is a
        degenerate singleton pinched at one of its ends. Both are returned: every
        element of every :class:`StablePartitionResult` on this arc's branch whose
        canonical-distance span meets ``[lo, hi]``.

        Args:
            trellis: The trellis whose :attr:`~.Trellis.Trellis.stable_partitions`
                to read. Run the partition first; before that this is empty.

        Returns:
            The matching elements, each carrying its own ``side`` and
            ``element_id``. Always empty for an unstable arc — the partition is a
            partition of the STABLE manifold.
        """
        if self.kind != "stable":
            return []
        tol = trellis.registry.cdist_tol
        lo = trellis.intersection(self.lo_id).stable_cdist
        hi = trellis.intersection(self.hi_id).stable_cdist
        found: list["PartitionInterval"] = []
        for result in trellis.stable_partitions:
            if result.branch_key != self.branch_key:
                continue
            for interval in result.intervals:
                if interval.lo_cdist > hi + tol or interval.hi_cdist < lo - tol:
                    continue
                found.append(interval)
        return found


class Region:
    """
    A face of the arrangement: a minimal cycle of manifold arcs.

    A region is bounded strictly by pieces of manifold — stable arcs and bridges —
    and contains no other manifold piece. It is the unit the dual graph will be
    built over, so it is identified COMBINATORIALLY, by its corner crossings, and
    carries no geometry until asked: :attr:`boundary_points`, :attr:`area` and
    :attr:`representative_point` are each built once, on demand, and a caller that
    only walks adjacency never pays for a polygon.

    Attributes:
        corners: The registry ids of the crossings at the corners of the face, in
            traversal order, canonicalised to start at the smallest id (see
            :func:`canonical_corners`). This tuple IS the region's identity.
        arcs: The boundary arcs in traversal order, each directed
            (:attr:`Arc.reverse`) the way the face runs along it. On a closed
            region ``arcs[i]`` runs from ``corners[i]`` to ``corners[i + 1]``. On
            an OPEN face the correspondence does not hold: its ``corners`` carry
            negative placeholder ids where the boundary reflects off a dangling
            end, and ``arcs`` lists only the real manifold pieces.
        is_closed: True when the face is bounded entirely by computed manifold —
            i.e. no dangling end. An open face (a manifold that simply stops
            before the next crossing) is not a region of the tangle, only of what
            has been computed so far, and is excluded from
            :attr:`~.Arrangement.Arrangement.regions`.
        is_minimal: True when the face encloses no other piece of manifold. It is
            False only for a closed face that swallows a whole other connected
            component of the trellis — two tangles with no computed heteroclinic
            crossing between them are two drawings on one plane, and the face
            traversal, which walks only along arcs, cannot see one from the other.
            Such faces are collected on
            :attr:`~.Arrangement.Arrangement.containing_faces`. A face is a REGION
            iff it is both closed and minimal.
        arrangement: The arrangement this face belongs to, used to resolve
            neighbours and geometry.
    """

    def __init__(
        self,
        corners: tuple[int, ...],
        arcs: list[Arc],
        arrangement: "Arrangement",
        is_closed: bool = True,
        is_minimal: bool = True,
    ):
        """
        Args:
            corners: Corner ids in traversal order, already canonicalised.
            arcs: Boundary arcs in the same traversal order.
            arrangement: The owning arrangement.
            is_closed: Whether the face has no dangling end.
            is_minimal: Whether the face encloses no other manifold piece. Set by
                :meth:`~.Arrangement.Arrangement._classify_minimality` after the
                faces are built, since deciding it needs their geometry.
        """
        self.corners = corners
        self.arcs = arcs
        self.arrangement = arrangement
        self.is_closed = is_closed
        self.is_minimal = is_minimal
        self._boundary_points: Optional[NDArray[np.float64]] = None
        self._representative_point: Optional[NDArray[np.float64]] = None

    # ── combinatorial views ─────────────────────────────────────────────────

    @property
    def stable_arcs(self) -> list[Arc]:
        """The boundary arcs that run along a stable manifold."""
        return [arc for arc in self.arcs if arc.kind == "stable"]

    @property
    def bridge_ids(self) -> list["BridgeId"]:
        """The :data:`BridgeId` of every unstable arc on the boundary."""
        return [
            arc.bridge_id
            for arc in self.arcs
            if arc.kind == "unstable" and arc.bridge_id is not None
        ]

    def neighbor_across(self, arc: Arc) -> Optional["Region"]:
        """
        The region on the other side of one of this region's boundary arcs.

        Args:
            arc: A boundary arc of this region.

        Returns:
            The region sharing that arc, or None when the face across it is open
            (or when ``arc`` does not bound this region).
        """
        for other in self.arrangement.regions_bounded_by(arc):
            if other is not self:
                return other
        return None

    # ── geometry (lazy) ─────────────────────────────────────────────────────

    def _require_closed(self, what: str) -> None:
        """Refuse a geometric question an open face has no answer to."""
        if not self.is_closed:
            raise ValueError(
                f"an open face has no {what}: its boundary runs off the end of a "
                "computed manifold, so it encloses nothing. Grow the manifolds "
                f"until the face closes (see Arrangement.open_faces). Face: "
                f"{self.corners}"
            )

    @property
    def boundary_points(self) -> NDArray[np.float64]:
        """
        The closed boundary polygon as an ``(N, 2)`` array.

        Built once, by concatenating each arc's polyline in traversal order and
        dropping the duplicated corner shared by consecutive arcs. The ring is not
        explicitly closed — :func:`~tanglepack.numerics.geometry.signed_polygon_area`
        and :func:`~tanglepack.numerics.geometry.point_in_polygon` both close it
        themselves.

        Raises:
            ValueError: If this face is open. Its boundary is not a closed curve —
                part of it is "the manifold stops here" — so stitching its arcs
                into a ring would silently invent an edge that is not manifold.
        """
        self._require_closed("boundary polygon")
        if self._boundary_points is None:
            trellis = self.arrangement.trellis
            pieces: list[NDArray[np.float64]] = []
            for index, arc in enumerate(self.arcs):
                poly = arc.polyline(trellis)
                pieces.append(poly if index == 0 else poly[1:])
            ring = np.vstack(pieces) if pieces else np.empty((0, 2))
            if len(ring) > 1 and np.allclose(ring[0], ring[-1]):
                ring = ring[:-1]
            self._boundary_points = ring
        return self._boundary_points

    @property
    def area(self) -> float:
        """
        The signed area enclosed by :attr:`boundary_points`.

        Negative for a bounded face: the arrangement traverses its bounded faces
        clockwise (see the :class:`~.Arrangement.Arrangement` Dev Notes). Take
        ``abs`` for a magnitude.

        Raises:
            ValueError: If this face is open (see :attr:`boundary_points`).
        """
        return signed_polygon_area(self.boundary_points)

    @property
    def representative_point(self) -> Optional[NDArray[np.float64]]:
        """
        A point standing in for the region, built once and cached.

        The first candidate is the cheap one: the mean of the midpoints of the
        boundary arcs, which lands inside for the lens- and lobe-shaped faces a
        tangle mostly produces. It is only a heuristic, though -- a face bounded
        by one short stable arc and one long meandering bridge is a crescent whose
        arc-midpoint mean falls well outside it -- so the result is checked with
        the ray cast and, when it fails, replaced by the scanline construction of
        :func:`~tanglepack.numerics.geometry.polygon_interior_point`, which is
        guaranteed interior for a simple polygon.

        Returns:
            The ``(2,)`` point, or None for a closed face whose boundary is too
            degenerate to have an interior (both arcs collapsed to one chord).

        Raises:
            ValueError: If this face is open (see :attr:`boundary_points`).
        """
        self._require_closed("representative point")
        if self._representative_point is None:
            trellis = self.arrangement.trellis
            mids = [arc.midpoint(trellis) for arc in self.arcs]
            usable = [m for m in mids if m is not None]
            candidate: Optional[NDArray[np.float64]] = None
            if usable:
                candidate = np.mean(np.vstack(usable), axis=0).astype(np.float64)
                if not self.contains(candidate):
                    logger.debug(
                        "Region %s is not star-shaped about its arc-midpoint mean; "
                        "falling back to the scanline interior point",
                        self.corners,
                    )
                    candidate = None
            if candidate is None:
                candidate = polygon_interior_point(self.boundary_points)
            self._representative_point = candidate
        return self._representative_point

    def verify_representative_point(self, *, tol: float = 1e-9) -> bool:
        """
        Whether :attr:`representative_point` really lies in this region.

        Args:
            tol: Boundary tolerance handed to the ray cast.

        Returns:
            True if the representative point exists and is contained.
        """
        point = self.representative_point
        return point is not None and self.contains(point, tol=tol)

    def contains(
        self,
        point: "NDArray[np.float64] | tuple[float, float]",
        *,
        tol: float = 1e-9,
    ) -> bool:
        """
        Whether a point lies in this region, boundary included.

        Args:
            point: The ``(x, y)`` to test.
            tol: Distance below which the point counts as on the boundary.

        Returns:
            True if inside or on the boundary.

        Raises:
            ValueError: If this face is open (see :attr:`boundary_points`).
        """
        return point_in_polygon(point, self.boundary_points, tol=tol)

    def __len__(self) -> int:
        return len(self.corners)

    def __repr__(self) -> str:
        if not self.is_closed:
            state = "open"
        elif not self.is_minimal:
            state = "containing"
        else:
            state = "closed"
        return f"Region({state}, corners={self.corners})"


def canonical_corners(corners: "tuple[int, ...] | list[int]") -> tuple[int, ...]:
    """
    Rotate a face's corner cycle to a canonical starting point.

    A face is a CYCLE, so the same face can be written starting at any of its
    corners. Canonicalising to the rotation that starts at the smallest id (and,
    when that id repeats, to the lexicographically smallest such rotation) makes
    the tuple a usable dict key. The cyclic ORDER is preserved, never sorted: the
    order is the traversal, and reversing it would name the face on the other side.

    Args:
        corners: The corner ids in traversal order.

    Returns:
        The canonical rotation, or the empty tuple for an empty input.
    """
    items = tuple(corners)
    if not items:
        return items
    rotations = [items[i:] + items[:i] for i in range(len(items))]
    return min(rotations)
