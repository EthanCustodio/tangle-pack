"""The forced Duffing example: an area-preserving, invertible, batch-capable map.

``examples/duffing.py`` is the expensive map the ``workers`` option is
measured on. Its leapfrog form must be what the library assumes of any map:
area-preserving (``|det J| = 1``), with an inverse that undoes it, and
batch-capable on axis 0. Its saddle must solve and its tangle cross at the
recipe's ``area_cutoff``, and the ``solve_ivp`` form must integrate the same
flow. Handed a CuPy array the leapfrog map runs as one device kernel, and
that kernel must give the numpy map's floats exactly (skipped without a GPU).
"""

from __future__ import annotations

import numpy as np
import pytest

import tanglepack
from tanglepack.examples.duffing import (
    duffing_map,
    duffing_map_inverse,
    duffing_map_ivp,
    duffing_map_ivp_inverse,
    saddle_guess,
)

#: Points near the saddle, on the loops and in both wells.
POINTS = np.array([[-0.06, 0.0], [0.5, 0.3], [-1.0, 0.2], [1.2, -0.4]])


def test_leapfrog_preserves_area():
    """The central-difference ``|det J|`` of the leapfrog map is 1 to 1e-8."""
    forward, h = duffing_map(), 1e-6
    for point in POINTS:
        jacobian = np.column_stack(
            [(forward(point + step) - forward(point - step)) / (2 * h) for step in np.eye(2) * h]
        )
        assert abs(abs(np.linalg.det(jacobian)) - 1) < 1e-8


def test_inverse_undoes_the_map():
    """``map_inv(map(x))`` and ``map(map_inv(x))`` return ``x`` to 1e-12."""
    forward, backward = duffing_map(), duffing_map_inverse()
    for point in POINTS:
        assert np.abs(backward(forward(point)) - point).max() < 1e-12
        assert np.abs(forward(backward(point)) - point).max() < 1e-12


def test_batch_equals_per_point():
    """A ``(2, N)`` batch maps to the per-point images exactly."""
    for fn in (duffing_map(), duffing_map_inverse()):
        assert np.array_equal(fn(POINTS.T).T, np.vstack([fn(point) for point in POINTS]))


def test_saddle_solves_and_the_tangle_crosses():
    """The saddle solves from ``saddle_guess`` and four steps of each manifold
    register a crossing of N besides the anchor at ``area_cutoff = 1e-5``."""
    session = tanglepack.TangleSession(duffing_map(), duffing_map_inverse(), area_cutoff=1e-5)
    fp = session.construct_fixed_point(saddle_guess())
    session.initialize_both_manifolds(fp)
    for stability in ("unstable", "stable"):
        session.grow_n_times(fp, stability, num_iterations=4)
    registry = session.compute_intersections([fp])
    assert any(not registry[iid].is_anchor for iid in registry.resolved_ids(3))


def test_ivp_form_integrates_the_same_flow():
    """At a step of ``T / 8000`` leapfrog agrees with ``solve_ivp`` to 1e-6,
    forward and backward, and the ``solve_ivp`` inverse undoes its map."""
    leapfrog, ivp = duffing_map(steps=8000), duffing_map_ivp(rtol=1e-11)
    leapfrog_inv, ivp_inv = duffing_map_inverse(steps=8000), duffing_map_ivp_inverse(rtol=1e-11)
    for point in POINTS:
        assert np.abs(leapfrog(point) - ivp(point)).max() < 1e-6
        assert np.abs(leapfrog_inv(point) - ivp_inv(point)).max() < 1e-6
        assert np.abs(ivp_inv(ivp(point)) - point).max() < 1e-8


def test_device_map_is_the_numpy_map():
    """On a CuPy batch the fused kernel's images equal the numpy loop's bit for bit,
    forward and backward, at the default step and at ten times it."""
    cp = pytest.importorskip("cupy")
    try:
        devices = cp.cuda.runtime.getDeviceCount()
    except Exception:  # noqa: BLE001 - any CUDA runtime failure means no device
        devices = 0
    if devices < 1:
        pytest.skip("CuPy is installed but no CUDA device is available")
    batch = np.random.default_rng(1).normal(0.0, 0.6, (2, 20_000))
    for fn in (duffing_map(), duffing_map_inverse(), duffing_map(steps=2000)):
        image = fn(cp.asarray(batch))
        assert isinstance(image, cp.ndarray)
        assert np.array_equal(cp.asnumpy(image), fn(batch))
