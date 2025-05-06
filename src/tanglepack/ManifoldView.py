from .BaseManifold import BaseManifold
from .DynamicalSystem import DynamicalSystem

class ManifoldView:
    """
    *Brings* a manifold and a system together without either object
    owning the other.  Pure references, no new state.
    """
    def __init__(self, manifold: BaseManifold, system: DynamicalSystem):

        self.manifold = manifold
        self.system   = system

        # cheap aliases so the numeric code looks neat
        self.root = manifold.root
        self.tail = manifold.tail
        self.stability = manifold.stability
        self.name = manifold.name
        self.stretch_param = manifold.stretch_param
        self.walk_fwd  = manifold.walk_fwd
        self.walk_back = manifold.walk_back

        if system.map_inv is None:
            # you can choose to raise here or let numeric code handle it
            raise ValueError("This system has no inverse map.")
        
        if manifold.stability == "unstable":
            self.map_fwd, self.map_back = system.map, system.map_inv

        else:   # stable
            self.map_fwd, self.map_back = system.map_inv, system.map
