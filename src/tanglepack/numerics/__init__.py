"""
Numerical layer of tanglepack.

This subpackage is the numerical engine: it defines the dynamical system, finds
fixed points, builds and grows the stable/unstable manifolds, and detects their
intersections — everything needed to produce a computed tangle. The topological
layer (:mod:`tanglepack.topology`) consumes its output; the cross-layer algorithms
(:mod:`tanglepack.loom`) sit on top of both. Nothing here depends on those layers.

The entry point for programmatic use is :class:`TangleWorkbench`.
"""

from __future__ import annotations

from .DynamicalSystem import DynamicalSystem

from .BasePoint import BasePoint
from .Point import Point
from .BranchPoint import BranchPoint

from .FixedPointSolver import FixedPointSolver
from .FixedPoint import FixedPoint

from .BaseManifold import BaseManifold
from .ManifoldView import ManifoldView
from .ManifoldMachine import ManifoldMachine
from .ManifoldInitializer import ManifoldInitializer

from .Bridge import Bridge
from .Tangle import Tangle
from .BridgeIterator import BridgeIterator
from .IterateInference import IterateInference
from .TangleWorkbench import TangleWorkbench
from . import graphviz

from .IterateTable import IterateTable
from .IntersectionRegistry import IntersectionRegistry
from .Intersection import Intersection, ManifoldKey, Stability

from .gpu import enable_gpu, disable_gpu

__all__ = [
    "DynamicalSystem",
    "BasePoint",
    "Point",
    "BranchPoint",
    "FixedPointSolver",
    "FixedPoint",
    "BaseManifold",
    "ManifoldView",
    "ManifoldMachine",
    "ManifoldInitializer",
    "Bridge",
    "Tangle",
    "BridgeIterator",
    "IterateInference",
    "TangleWorkbench",
    "graphviz",
    "IterateTable",
    "IntersectionRegistry",
    "Intersection",
    "ManifoldKey",
    "Stability",
    "enable_gpu",
    "disable_gpu",
]
