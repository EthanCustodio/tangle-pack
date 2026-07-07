from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from ..numerics.Intersection import ManifoldKey

"""
Dev Notes:

These dataclasses are pure *result containers*. They hold the output of the
topological algorithms (Compute-Pseudoneighbors, Is-Strong-Pip, ...) so that a
single Trellis object can carry both the input trellis and every derived
classification. No algorithm logic lives here — only the schema of what the
algorithms produce.

Intersections are referenced by their integer registry ID rather than by object
so that a result survives a registry rebuild and is cheap to serialise.

Open question: a Hole is currently a single phase-space coordinate. When the
region-punching / resonance-zone bookkeeping is implemented we may want a Hole
to also reference the bounding bridge(s) or the enclosed region.
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
        side: Which side of the oriented stable manifold (looking toward the
            anchor point) the hole lies on. None until classified.
        interior: True if the hole lies inside the resonance zone (its stable
            position is below the zone's cut on that branch), False if outside.
            None when no resonance-zone cut is known.
        bounding_ids: Registry IDs of the two intersections bounding the region
            the hole is punched in — the pair itself for a directly punched
            hole, the containing bridge's endpoints for a propagated one. These
            are the partition boundaries of the Stable Manifold Partition
            Algorithm.
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
    side: Optional[Literal["left", "right"]] = None
    interior: Optional[bool] = None
    bounding_ids: Optional[tuple[int, int]] = None
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

    Attributes:
        lo_id: Registry ID of the lower (toward-anchor) boundary intersection,
            or None when the interval starts at the anchor point (cdist 0).
        hi_id: Registry ID of the upper (outward) boundary intersection, or
            None when the interval runs to the end of the computed branch.
        lo_cdist: Stable canonical distance of the lower boundary.
        hi_cdist: Stable canonical distance of the upper boundary.
        closed_lo: True if the lower boundary point belongs to the interval.
        closed_hi: True if the upper boundary point belongs to the interval.
    """

    lo_id: Optional[int]
    hi_id: Optional[int]
    lo_cdist: float
    hi_cdist: float
    closed_lo: bool
    closed_hi: bool


@dataclass
class StablePartitionResult:
    """
    The partition of one stable branch induced by the holes on one side.

    The stable manifold is oriented looking toward the anchor point; holes
    punched between pseudoneighbor pairs fall on its left or right side, and
    each side induces its own partition (see
    Stable_Manifold_Partition_Algorithm.pdf). The branch splits at the
    resonance-zone cut into an interior part (inside the zone) and an exterior
    part, each partitioned separately.

    Attributes:
        branch_key: Manifold key of the partitioned stable branch.
        side: Which side's holes this partition is built from.
        interior_intervals: Intervals below the resonance-zone cut, ordered
            from the anchor outward.
        exterior_intervals: Intervals above the cut, ordered outward.
        cut_cdist: Stable canonical distance of the resonance-zone cut on this
            branch (the strong pip or its iterate), or None if no cut is known
            (then every interval is reported as interior).
    """

    branch_key: "ManifoldKey"
    side: Literal["left", "right"]
    interior_intervals: list[PartitionInterval] = field(default_factory=list)
    exterior_intervals: list[PartitionInterval] = field(default_factory=list)
    cut_cdist: Optional[float] = None

    @property
    def intervals(self) -> list[PartitionInterval]:
        """All intervals, interior then exterior, ordered from the anchor outward."""
        return self.interior_intervals + self.exterior_intervals
