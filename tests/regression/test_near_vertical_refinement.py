"""Regression: refinement curvature must be rotation invariant (near-vertical).

``_curvature_area`` fits a parabola and a line through manifold points to decide
whether a segment needs refining. The curvature area of a curved segment is a
geometric quantity and must not depend on the segment's orientation. The bug was
that the fit was expressed as ``y`` of ``x``: rotating a well-behaved curved
segment toward vertical made the x-values coincide, the Vandermonde matrix in
``_parabolic_fit`` went singular (``LinAlgError``) and the line fit divided by
``~0``. ``refine_two_points`` caught the error and skipped the pair, so
near-vertical stretches of a manifold silently lost resolution.

The fit now runs in the chord frame (``_chord_frame`` rotates the segment's chord
onto the x-axis before fitting), so the answer no longer depends on orientation.
This pins it: the curvature area of a curved segment equals that of the same
shape rotated 90 degrees.
"""

from __future__ import annotations

import numpy as np
import pytest

from tanglepack import BaseManifold, ManifoldMachine, ManifoldView, Point, DynamicalSystem


def _identity_system():
    return DynamicalSystem(lambda p: np.asarray(p, float), lambda p: np.asarray(p, float))


def _curved_manifold(coords):
    """Build a 4-point unstable manifold from (x, y) coords pA, p0, p1, pB."""
    pa, p0, p1, pb = (Point(x, y, cdist=float(i)) for i, (x, y) in enumerate(coords))
    pa.forward, p0.backward = p0, pa
    p0.forward, p1.backward = p1, p0
    p1.forward, pb.backward = pb, p1
    manifold = BaseManifold(
        pa, "unstable", 1.0, fixed_point=None, tail=pb, manifold_key=None
    )
    return manifold, (p0, p1)


def _area_for(coords):
    system = _identity_system()
    machine = ManifoldMachine(system)
    manifold, (p0, p1) = _curved_manifold(coords)
    viewer = ManifoldView(manifold, system)
    return machine._curvature_area((p0, p1), viewer)


@pytest.mark.regression
def test_curvature_area_is_rotation_invariant():
    horizontal = [(0.0, 0.0), (1.0, 1.0), (2.0, 1.0), (3.0, 0.0)]
    # rotate 90 degrees: (x, y) -> (-y, x); p0, p1 now share an x (vertical segment)
    vertical = [(-y, x) for (x, y) in horizontal]

    area_h = _area_for(horizontal)
    area_v = _area_for(vertical)

    assert np.isfinite(area_h) and np.isfinite(area_v)
    assert np.isclose(area_h, area_v, rtol=1e-6), (
        f"curvature area changed under rotation: {area_h} (horizontal) vs "
        f"{area_v} (vertical)"
    )
