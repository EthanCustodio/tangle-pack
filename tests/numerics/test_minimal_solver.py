"""Minimal fixed-point solver sanity checks (replaces the stale solver tests)."""

from __future__ import annotations

import numpy as np
import pytest

from tanglepack import DynamicalSystem, FixedPointSolver
from tanglepack.examples import (
    HENON_K10,
    HENON_P3,
    henon_map as _henon_map_factory,
    henon_map_inverse as _henon_map_inverse_factory,
    saddle_guesses,
)

henon_map = _henon_map_factory(*HENON_K10)
henon_map_inverse = _henon_map_inverse_factory(*HENON_K10)


def test_located_fixed_point_is_actually_fixed():
    system = DynamicalSystem(henon_map, henon_map_inverse)
    solver = FixedPointSolver(system)
    fp = solver.construct_fixed_point(saddle_guesses(*HENON_K10)["saddle"])

    coord = np.asarray(fp.coordinates[0], dtype=float).ravel()[:2]
    image = np.asarray(henon_map(coord), dtype=float).ravel()[:2]
    assert np.allclose(image, coord, atol=1e-9), (
        f"map({coord}) = {image} is not fixed"
    )


def test_eigenvalues_are_saddle_like():
    """An area-preserving saddle has eigenvalues lambda and 1/lambda."""
    system = DynamicalSystem(henon_map, henon_map_inverse)
    solver = FixedPointSolver(system)
    fp = solver.construct_fixed_point(saddle_guesses(*HENON_K10)["saddle"])

    u = float(np.abs(np.asarray(fp.unstable_eigenvalues, dtype=float).ravel()[0]))
    s = float(np.abs(np.asarray(fp.stable_eigenvalues, dtype=float).ravel()[0]))
    assert s < 1.0 < u, f"not a saddle: |stable|={s}, |unstable|={u}"
    assert np.isclose(u * s, 1.0, rtol=1e-6), (
        f"eigenvalue product {u * s} != 1 (area preservation)"
    )


# ── 1.13 solver hardening ──────────────────────────────────────────────────


def constant_shift(point):
    """A map with no fixed point at all: f(p) = p + 1."""
    return np.asarray(point, dtype=float) + 1.0


def constant_shift_inverse(point):
    return np.asarray(point, dtype=float) - 1.0


def elliptic_map(point):
    """Rigid rotation by 0.4 rad about (1, 1): eigenvalues exp(+-0.4i).

    Written element-wise so it also accepts the batched (2, ...) arrays that
    scipy's finite-difference Jacobian feeds it.
    """
    theta = 0.4
    dx = point[0] - 1.0
    dy = point[1] - 1.0
    return np.array(
        [
            1.0 + np.cos(theta) * dx - np.sin(theta) * dy,
            1.0 + np.sin(theta) * dx + np.cos(theta) * dy,
        ]
    )


def elliptic_map_inverse(point):
    theta = -0.4
    dx = point[0] - 1.0
    dy = point[1] - 1.0
    return np.array(
        [
            1.0 + np.cos(theta) * dx - np.sin(theta) * dy,
            1.0 + np.sin(theta) * dx + np.cos(theta) * dy,
        ]
    )


def test_non_convergence_raises_with_the_minpack_message():
    """fsolve failing (ier != 1) must raise, not silently return the guess."""
    system = DynamicalSystem(constant_shift, constant_shift_inverse)
    solver = FixedPointSolver(system)

    with pytest.raises(ValueError, match="did not converge"):
        solver.compute_fixed_point(np.atleast_2d([1.0, 1.0]))


def test_non_saddle_fixed_point_raises():
    """A converged but elliptic fixed point is a user-input error, not a saddle."""
    system = DynamicalSystem(elliptic_map, elliptic_map_inverse)
    solver = FixedPointSolver(system)

    with pytest.raises(ValueError, match="saddle"):
        solver.construct_fixed_point([1.1, 0.9])


def test_orient_hook_is_honoured():
    """The optional orient hook may re-sign the eigenvectors it is handed."""
    system = DynamicalSystem(henon_map, henon_map_inverse)
    guess = saddle_guesses(*HENON_K10)["saddle"]
    baseline = FixedPointSolver(system).construct_fixed_point(guess)

    seen: list[int] = []

    def flip_unstable(orbit_index, unstable, stable):
        seen.append(orbit_index)
        return -unstable, stable

    oriented = FixedPointSolver(system, orient=flip_unstable).construct_fixed_point(
        guess
    )

    assert seen == [0]
    assert np.allclose(
        oriented.unstable_eigenvectors[0], -baseline.unstable_eigenvectors[0]
    )
    assert np.allclose(oriented.stable_eigenvectors[0], baseline.stable_eigenvectors[0])


def test_period_three_orbit_at_k_two_is_accepted():
    """The exact k=2 period-3 orbit (residual 0) must converge and be a saddle."""

    henon_k2 = _henon_map_factory(*HENON_P3)
    henon_k2_inverse = _henon_map_inverse_factory(*HENON_P3)

    solver = FixedPointSolver(DynamicalSystem(henon_k2, henon_k2_inverse))
    fp = solver.construct_fixed_point(saddle_guesses(*HENON_P3)["period_3"])

    for i in range(3):
        u = float(np.abs(np.asarray(fp.unstable_eigenvalues[i], dtype=float).ravel()[0]))
        s = float(np.abs(np.asarray(fp.stable_eigenvalues[i], dtype=float).ravel()[0]))
        assert s < 1.0 < u
