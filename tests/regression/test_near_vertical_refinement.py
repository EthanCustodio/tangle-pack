"""Regression: refinement curvature must be rotation invariant (near-vertical).

The curvature-area kernel fits a parabola and a line through manifold points to decide
whether a segment needs refining. The curvature area of a curved segment is a
geometric quantity and must not depend on the segment's orientation. The bug was
that the fit was expressed as ``y`` of ``x``: rotating a well-behaved curved
segment toward vertical made the x-values coincide, the Vandermonde matrix in
``_parabolic_fit`` went singular (``LinAlgError``) and the line fit divided by
``~0``. ``refine_two_points`` caught the error and skipped the pair, so
near-vertical stretches of a manifold silently lost resolution.

The fit now runs in the chord frame (the segment's chord is rotated onto the
x-axis before fitting), so the answer no longer depends on orientation. This pins
it on the vectorized kernel the refiner uses (``_curvature_area_batch``): the
curvature area of a curved segment equals that of the same shape rotated 90
degrees.
"""

from __future__ import annotations

import numpy as np

from tanglepack import ManifoldMachine


def _area_for(coords: list[tuple[float, float]]) -> float:
    """Curvature area of the middle segment ``p0 -> p1`` of ``pA, p0, p1, pB``."""
    pa, p0, p1, pb = (np.asarray([c], dtype=float) for c in coords)
    return float(ManifoldMachine._curvature_area_batch(p0, p1, pa, pb)[0])


def test_curvature_area_is_rotation_invariant():
    horizontal = [(0.0, 0.0), (1.0, 1.0), (2.0, 1.0), (3.0, 0.0)]
    # rotate 90 degrees: (x, y) -> (-y, x); p0, p1 now share an x (vertical segment)
    vertical = [(-y, x) for (x, y) in horizontal]

    area_h = _area_for(horizontal)
    area_v = _area_for(vertical)

    assert np.isfinite(area_h) and np.isfinite(area_v)
    assert area_h > 0, "a curved segment must have a positive curvature area"
    assert np.isclose(area_h, area_v, rtol=1e-6), (
        f"curvature area changed under rotation: {area_h} (horizontal) vs "
        f"{area_v} (vertical)"
    )
