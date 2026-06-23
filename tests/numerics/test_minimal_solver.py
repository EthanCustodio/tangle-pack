"""Minimal fixed-point solver sanity checks (replaces the stale solver tests)."""

from __future__ import annotations

import numpy as np

from tanglepack import DynamicalSystem, FixedPointSolver


def henon_map(point):
    k, b = 10, 1
    x, y = point[0], point[1]
    return np.array([y - k + x**2, -b * x])


def henon_map_inverse(point):
    k, b = 10, 1
    x, y = point[0], point[1]
    return np.array([-y / b, x + k - (y**2) / (b**2)])


def test_located_fixed_point_is_actually_fixed():
    system = DynamicalSystem(henon_map, henon_map_inverse)
    solver = FixedPointSolver(system)
    fp = solver.construct_fixed_point([4, -4], 1)

    coord = np.asarray(fp.coordinates[0], dtype=float).ravel()[:2]
    image = np.asarray(henon_map(coord), dtype=float).ravel()[:2]
    assert np.allclose(image, coord, atol=1e-9), (
        f"map({coord}) = {image} is not fixed"
    )


def test_eigenvalues_are_saddle_like():
    """An area-preserving saddle has eigenvalues lambda and 1/lambda."""
    system = DynamicalSystem(henon_map, henon_map_inverse)
    solver = FixedPointSolver(system)
    fp = solver.construct_fixed_point([4, -4], 1)

    u = float(np.abs(np.asarray(fp.unstable_eigenvalues, dtype=float).ravel()[0]))
    s = float(np.abs(np.asarray(fp.stable_eigenvalues, dtype=float).ravel()[0]))
    assert s < 1.0 < u, f"not a saddle: |stable|={s}, |unstable|={u}"
    assert np.isclose(u * s, 1.0, rtol=1e-6), (
        f"eigenvalue product {u * s} != 1 (area preservation)"
    )
