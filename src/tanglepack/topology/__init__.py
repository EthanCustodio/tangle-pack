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
    Hole,
    PartitionInterval,
    PseudoneighborPair,
    Region,
    StablePartitionResult,
    StrongPipResult,
    side_after,
)
from .plotting import (
    CartoonLayout,
    ZoneArcOrder,
    ZoneLayout,
    dual_graph_cartoon_layout,
    dual_graph_zone_layout,
    dual_graph_cartoon_legend_handles,
    dual_graph_legend_handles,
    face_point,
    minimal_trellis_legend_handles,
    plot_bridges_by_class,
    plot_dual_graph,
    plot_dual_graph_cartoon,
    plot_itinerary_table,
    plot_minimal_trellis,
    plot_stable_partition,
    plot_transition_graph,
    plot_walk,
    stable_node_point,
    walk_legend_handles,
)
from .BridgeClass import (
    BridgeClass,
    BridgeClassEntry,
    BridgeClassTable,
    bridge_classes,
    zone_label,
)
from .DualGraph import DualGraph
from .MinimalTrellis import (
    MinimalTrellis,
    minimal_trellis,
)
from .PartitionFamily import (
    HomotopyPartition,
    IteratedHomotopyPartition,
)
from .DualWalk import (
    MembersDisagree,
    Walk,
    WalkReason,
)
from .SymbolicDynamics import (
    ClassDynamics,
    SymbolicDynamics,
    free_reduce,
    inert_letters,
    letter,
    symbolic_dynamics,
)
from .StablePartition import check_holes_share_bridge_side
from .diagnostics import bridge_row_diagnosis

__all__ = [
    "Trellis",
    "TrellisBranch",
    "Arrangement",
    "Arc",
    "side_after",
    "MinimalTrellis",
    "minimal_trellis",
    "HomotopyPartition",
    "IteratedHomotopyPartition",
    "BridgeClass",
    "BridgeClassEntry",
    "BridgeClassTable",
    "DualGraph",
    "Walk",
    "WalkReason",
    "MembersDisagree",
    "ClassDynamics",
    "SymbolicDynamics",
    "inert_letters",
    "letter",
    "free_reduce",
    "symbolic_dynamics",
    "Region",
    "Hole",
    "PartitionInterval",
    "PseudoneighborPair",
    "StablePartitionResult",
    "StrongPipResult",
    "bridge_classes",
    "bridge_row_diagnosis",
    "check_holes_share_bridge_side",
    "dual_graph_legend_handles",
    "zone_label",
    "face_point",
    "plot_dual_graph",
    "plot_minimal_trellis",
    "minimal_trellis_legend_handles",
    "stable_node_point",
    "plot_stable_partition",
    "plot_walk",
    "plot_bridges_by_class",
    "plot_transition_graph",
    "plot_itinerary_table",
    "walk_legend_handles",
    "plot_dual_graph_cartoon",
    "dual_graph_cartoon_layout",
    "dual_graph_cartoon_legend_handles",
    "CartoonLayout",
    "ZoneArcOrder",
    "ZoneLayout",
    "dual_graph_zone_layout",
]
