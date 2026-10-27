"""
Regression (INV-5): the resonance zone of an inversion point is both of its loops.

The zone ring used to match a cut to the next by ORBIT INDEX alone. Both cuts of
the k=4 inversion saddle have orbit index 0, so one unstable branch was walked
twice and the other never: the polygon was neither loop. The ring now follows the
full branch key (:func:`~tanglepack.loom.ResonanceZone.ring_cycles`), and the
zone is the union of the two loops through the periodic point, each a cut's
unstable arc out to its pip and its stable arc back. The map carries one loop
onto the other, so the two enclose the same area.
"""

from __future__ import annotations

import numpy as np
import pytest

from cases import build_inversion
from tanglepack.numerics.geometry import arc_polyline, signed_polygon_area


def _loop_area(workbench, cut, tol: float) -> float:
    """The area of one cut's loop: unstable arc to the pip, stable arc back."""
    pip = np.asarray(cut.coords, dtype=float)
    out = arc_polyline(workbench.manifolds[cut.manifold_a_key], 0.0, cut.unstable_cdist, tol=tol)
    back = arc_polyline(workbench.manifolds[cut.manifold_b_key], 0.0, cut.stable_cdist, tol=tol)
    return abs(signed_polygon_area(np.vstack([*out, pip, *back[::-1]])))


def test_inversion_zone_is_the_union_of_its_loops():
    case = build_inversion(through="zones")
    (zone,) = case.zones
    workbench = case.workbench
    tol = float(workbench.intersection_registry.cdist_tol)
    loops = [_loop_area(workbench, cut, tol) for cut in zone.boundary_intersections]

    assert len(loops) == case.fixed_point.k_value
    assert loops[0] == pytest.approx(loops[1], rel=1e-3)
    assert zone.area == pytest.approx(sum(loops), rel=1e-6)
