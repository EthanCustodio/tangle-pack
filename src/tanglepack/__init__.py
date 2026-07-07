"""Top-level tanglepack API.

The library is organized into three subpackages:

* :mod:`tanglepack.numerics` — the numerical engine (dynamical system, fixed
  points, manifold growth, intersection detection). Entry point: ``TangleWorkbench``.
* :mod:`tanglepack.topology` — the topological view of a computed tangle
  (``Trellis``, strong pips, pseudoneighbors).
* :mod:`tanglepack.loom` — cross-layer ("meta") algorithms that read topology
  results and act on the numerical layer (``TangleSession``, resonance zones).

The most-used names from every layer are re-exported here, so ``tanglepack.X``
keeps working regardless of which subpackage ``X`` now lives in.
"""

from . import numerics, topology, loom

from .numerics import (
    DynamicalSystem,
    BasePoint,
    Point,
    BranchPoint,
    FixedPointSolver,
    FixedPoint,
    BaseManifold,
    ManifoldView,
    ManifoldMachine,
    ManifoldInitializer,
    Bridge,
    Tangle,
    TangleWorkbench,
    IterateTable,
    IntersectionRegistry,
    Intersection,
    ManifoldKey,
    enable_gpu,
    disable_gpu,
)

from .topology import (
    Trellis,
    TrellisBranch,
    Hole,
    PartitionInterval,
    PseudoneighborPair,
    StablePartitionResult,
    StrongPipResult,
)

from .loom import (
    TangleSession,
    ResonanceZone,
    define_resonance_zone,
    trim_stable_at_intersection,
)
