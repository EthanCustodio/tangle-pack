"""Regression: refinement curvature is not rotation-invariant (near-vertical fails).

``_curvature_area`` fits a parabola ``y = a x^2 + b x + c`` and a line ``y = m x + b``
through manifold points. The curvature area of a curved segment is a geometric
quantity and must not depend on the segment's orientation. But because the fit is
expressed as ``y`` of ``x``, rotating a well-behaved curved segment toward vertical
makes the x-values coincide: the Vandermonde matrix in ``_parabolic_fit`` goes
singular (``LinAlgError``) and the line fit divides by ``~0``. ``refine_two_points``
then catches the error and skips the pair, so near-vertical stretches lose
resolution.

This test computes the curvature area of a curved segment, then of the *same shape
rotated 90 degrees*, and asserts they match. Today the rotated (vertical) case
raises / disagrees, so it xfails. The fix (fit in rotated/arclength coordinates)
makes the two equal.
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
    manifold = BaseManifold(pa, "unstable", 1.0, fixed_point=None, tail=pb)
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
