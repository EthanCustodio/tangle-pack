from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from ..Intersection import ManifoldKey

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
    """

    coords: tuple[float, float]
    near_intersection_id: int
    pair: Optional["PseudoneighborPair"] = None


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
    """

    intersection_a: int
    intersection_b: int
    branch_key: Optional["ManifoldKey"] = None
    hole: Optional[Hole] = None

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
