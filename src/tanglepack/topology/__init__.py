"""
Topological layer of tanglepack.

This subpackage consumes the finished output of the numerical layer (manifolds,
intersections, bridges) via a TangleWorkbench and exposes it as a Trellis — the
object the topological algorithms read from and write their results to.

The numerical and topological layers are intentionally separable: nothing here
grows manifolds or computes intersections; it only interprets what the numerical
layer has already produced.
"""

from __future__ import annotations

from .Trellis import Trellis
from .TrellisBranch import TrellisBranch
from .Arrangement import Arrangement
from .TopologyResults import (
    Arc,
    ElementRef,
    Hole,
    PartitionInterval,
    PseudoneighborPair,
    Region,
    StablePartitionResult,
    StrongPipResult,
    canonical_corners,
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
from .plotting import (
    dual_graph_legend_handles,
    face_point,
    plot_dual_graph,
    plot_dual_graph_curved,
    plot_stable_partition,
)
from .BridgeClass import (
    BridgeClass,
    bridge_classes,
    class_sort_key,
    element_sort_key,
)
from .DualGraph import (
    ArcNode,
    DualGraph,
    FaceNode,
)
from .StablePartition import (
    bridge_for_pair,
    bridge_row_violation,
    bridge_side_violations,
    check_bridge_rows_consistent,
    check_holes_share_bridge_side,
    owns_cdist,
    partition_stable_manifold,
    propagate_reference_holes,
    punch_holes,
    row_of_end,
    rows_of_bridge,
    span_contains,
)

__all__ = [
    "Trellis",
    "TrellisBranch",
    "Arrangement",
    "Arc",
    "ArcNode",
    "BridgeClass",
    "DualGraph",
    "FaceNode",
    "ElementRef",
    "Region",
    "canonical_corners",
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
    "bridge_classes",
    "bridge_for_pair",
    "bridge_row_violation",
    "bridge_side_violations",
    "check_bridge_rows_consistent",
    "check_holes_share_bridge_side",
    "class_sort_key",
    "dual_graph_legend_handles",
    "element_sort_key",
    "face_point",
    "owns_cdist",
    "partition_stable_manifold",
    "plot_dual_graph",
    "plot_dual_graph_curved",
    "plot_stable_partition",
    "propagate_reference_holes",
    "punch_holes",
    "row_of_end",
    "rows_of_bridge",
    "span_contains",
]
