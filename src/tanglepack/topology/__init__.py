"""
Topological layer of tanglepack.

This subpackage consumes the finished output of the numerical layer (manifolds,
intersections, bridges) via a TangleWorkbench and exposes it as a Trellis — the
object the topological algorithms read from and write their results to.

The numerical and topological layers are intentionally separable: nothing here
grows manifolds or computes intersections; it only interprets what the numerical
layer has already produced.
"""

from .Trellis import Trellis
from .TrellisBranch import TrellisBranch
from .TopologyResults import (
    Hole,
    PartitionInterval,
    PseudoneighborPair,
    StablePartitionResult,
    StrongPipResult,
)
from .StrongPip import (
    is_strong_pip,
    classify_strong_pips,
    forward_stable_branch_cycle,
)
from .Pseudoneighbor import (
    compute_pseudoneighbors,
    extend_pseudoneighbor_trajectories,
    forward_unstable_branch_cycle,
)
from .StablePartition import (
    bridge_for_pair,
    partition_stable_manifold,
    plot_stable_partition,
    propagate_reference_holes,
    punch_holes,
)

__all__ = [
    "Trellis",
    "TrellisBranch",
    "Hole",
    "PartitionInterval",
    "PseudoneighborPair",
    "StablePartitionResult",
    "StrongPipResult",
    "is_strong_pip",
    "classify_strong_pips",
    "forward_stable_branch_cycle",
    "compute_pseudoneighbors",
    "extend_pseudoneighbor_trajectories",
    "forward_unstable_branch_cycle",
    "bridge_for_pair",
    "partition_stable_manifold",
    "plot_stable_partition",
    "propagate_reference_holes",
    "punch_holes",
]
