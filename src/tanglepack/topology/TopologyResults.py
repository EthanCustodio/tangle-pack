from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from ..numerics.Bridge import BridgeId
    from ..numerics.Intersection import ManifoldKey

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

"""
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
