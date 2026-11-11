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
    * :class:`Topology` — one recompute of the topology (stages 3-9) for a state and pip choice.
    * :class:`BridgeAlphabet` — the persistent letters of the session's bridge classes.
    * :func:`cartoon_zones` / :func:`plot_topology_cartoon` / :func:`plot_resonance_zones` — what the
      drawings need from the zones.
"""

from __future__ import annotations

from .ResonanceZone import (
    BoundaryArc,
    ResonanceZone,
    define_resonance_zone,
)
from .Blast import BlastResult, BlastStep, blast_zone
from .BridgeAlphabet import BridgeAlphabet
from .cartoon import cartoon_zones, plot_resonance_zones, plot_topology_cartoon
from .TangleSession import TangleSession
from .Topology import Topology

__all__ = [
    "TangleSession",
    "Topology",
    "ResonanceZone",
    "BoundaryArc",
    "define_resonance_zone",
    "blast_zone",
    "BlastResult",
    "BlastStep",
    "BridgeAlphabet",
    "cartoon_zones",
    "plot_resonance_zones",
    "plot_topology_cartoon",
]
