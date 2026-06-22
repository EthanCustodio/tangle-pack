"""
Loom — the cross-layer ("meta") algorithms of tanglepack.

A loom weaves separate threads into a finished fabric. This subpackage weaves the
numerical *threads* (the manifolds — :mod:`tanglepack.numerics`) and the topological
*lattice* (the trellis — :mod:`tanglepack.topology`) into finished structures such as
resonance zones. Algorithms here read topological results (e.g. a chosen strong pip)
and act on the numerical layer (e.g. trimming a manifold and recomputing crossings).

Nothing in :mod:`tanglepack.numerics` or :mod:`tanglepack.topology` imports from here —
loom sits on top of both.

Entry points:
    * :class:`TangleSession` — the user-friendly facade tying both layers together.
    * :func:`define_resonance_zone` / :class:`ResonanceZone` — the first loom algorithm.
"""

from .ResonanceZone import (
    BoundaryArc,
    ResonanceZone,
    define_resonance_zone,
    trim_stable_at_intersection,
)
from .Blast import BlastResult, BlastStep, blast_zone
from .TangleSession import TangleSession

__all__ = [
    "TangleSession",
    "ResonanceZone",
    "BoundaryArc",
    "define_resonance_zone",
    "trim_stable_at_intersection",
    "blast_zone",
    "BlastResult",
    "BlastStep",
]
