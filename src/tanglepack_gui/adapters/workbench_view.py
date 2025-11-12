from __future__ import annotations
import numpy as np


def manifold_arrays_for_fp(workbench, fp, stability=None):
    """Yield (key, Nx2 array) for manifolds of fp (optionally filter by stability)."""
    for (kfp, kstab, oi, bi), M in workbench.manifolds.items():
        if kfp is fp and (stability is None or kstab == stability):
            yield (kstab, oi, bi), M.get_point_array()  # Nx2
