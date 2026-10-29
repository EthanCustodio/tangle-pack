"""The forced Duffing example: an area-preserving, invertible, batch-capable map.

``examples/duffing.py`` is the expensive map the ``workers`` option is
measured on. Its leapfrog form must be what the library assumes of any map:
area-preserving (``|det J| = 1``), with an inverse that undoes it, and
batch-capable on axis 0. Its saddle must solve and its tangle cross at the
recipe's ``area_cutoff``, and the ``solve_ivp`` form must integrate the same
flow.
"""

from __future__ import annotations

import numpy as np

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
