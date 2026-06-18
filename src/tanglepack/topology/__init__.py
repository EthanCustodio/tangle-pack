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
from .TopologyResults import Hole, PseudoneighborPair, StrongPipResult
from .StrongPip import (
    is_strong_pip,
    classify_strong_pips,
    forward_stable_branch_cycle,
)

__all__ = [
    "Trellis",
    "TrellisBranch",
    "Hole",
    "PseudoneighborPair",
    "StrongPipResult",
    "is_strong_pip",
    "classify_strong_pips",
    "forward_stable_branch_cycle",
]
