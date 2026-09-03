"""
Binds a manifold to a dynamical system without either owning the other.

:class:`ManifoldView` is a transient bundle of references that resolves
"forward" and "backward" to the right map for the manifold's stability, so the
numerical code can be written once for both stabilities.
"""

from __future__ import annotations

from .BaseManifold import BaseManifold
from .DynamicalSystem import DynamicalSystem


class ManifoldView:
    """
    *Brings* a manifold and a system together without either object
    owning the other.  Pure references, no new state.
    """

    def __init__(self, manifold: BaseManifold, system: DynamicalSystem) -> None:
        """
        Bind one manifold to one system.

        Args:
            manifold (BaseManifold): The curve being read or grown.
            system (DynamicalSystem): The system whose map drives it.
        """

        self.manifold = manifold
        self.system = system

        # cheap aliases so the numeric code looks neat
        self.root = manifold.root
        self.tail = manifold.tail
        self.stability = manifold.stability
        self.name = manifold.name
        self.stretch_param = manifold.stretch_param
        self.walk_fwd = manifold.walk_fwd
        self.walk_back = manifold.walk_back
        self.get_point_array = manifold.get_point_array
        self.get_cdist_array = manifold.get_cdist_array

        if manifold.stability == "unstable":
            self.map_fwd, self.map_back = system.map, system.map_inv
            self.map_fwd_batch = system.map_batch
            self.map_back_batch = system.map_inv_batch

        else:  # stable
            self.map_fwd, self.map_back = system.map_inv, system.map
            self.map_fwd_batch = system.map_inv_batch
            self.map_back_batch = system.map_batch
