"""Regression: high-stretch period-3 growth must not scramble the manifold.

At the period-3 Hénon tangle for k=2.1 the per-step eigenvalue factor is ~4, so a
few growth iterations stretch an arm by orders of magnitude and pack many passes
near the saddle. Canonical distances there fall to floating-point spacing: two
geometrically distinct points end up with cdists that differ by a single ULP (or
collide exactly). cdist can no longer order them, and the manifold's linked list
darts out to a misplaced point and back -- the "super spiky inner tangle".

The real, observable invariant is geometric: the manifold is a smooth simple curve,
so no node may sit far off the chord of its neighbours (no over-long polyline
segment). An earlier attempt instead tried to force cdist to be strictly
injective (nudging every tie up by ``nextafter``); that made the cdist *values*
strictly increasing while leaving the *geometry* scrambled -- a cdist-only check
passed on a visibly broken manifold. This test therefore pins the geometry
directly. cdist is still required to be non-decreasing (ties are allowed: where a
fold collapses adjacent cdists below a ULP the refiner bridges the gap with
equal-cdist points, which are spliced geometrically and never re-sorted).
"""

from __future__ import annotations

import numpy as np
import pytest

import tanglepack
from invariants import assert_cdist_monotonic, assert_no_geometric_spikes


def _make(k: float):
    def henon_map(point):
        x, y = point
        return np.stack([y - k + x**2, -x], axis=0)

    def henon_map_inverse(point):
        x, y = point
        return np.stack([-y, x + k - y**2], axis=0)

    def henon_jacobian(point):
        x, y = point
        return np.array([[2 * x, 1], [-1, 0]])

    return henon_map, henon_map_inverse, henon_jacobian


@pytest.mark.regression
def test_period3_high_stretch_growth_is_not_scrambled():
    henon_map, henon_map_inverse, henon_jacobian = _make(2.1)
    wb = tanglepack.TangleWorkbench(henon_map, henon_map_inverse, henon_jacobian)
    wb._man_machine.area_cutoff = 1e-7

    fp3 = wb.construct_fixed_point([[0, 1], [-1, 0], [-1, 1]])
    wb.orient_eigenvectors(
        fp3, {"unstable": np.array([0, 1]), "stable": np.array([-1, -1])}
    )
    wb.initialize_both_manifolds(fp3)
    # Four unstable iterations reach the near-saddle accumulation where cdists fall
    # to float spacing -- enough to trigger the scramble (the buggy strictify build
    # left 11/7/37 spikes here). The arm escapes by the fifth iterate, so stop at
    # four. This shares the merge/refine code path that produced the user's stable
    # symptom (fp3 stable grown 5x), so it pins the same bug far more cheaply than
    # the heavy stable growth.
    wb.grow_n_times(fp3, "unstable", num_iterations=4)

    for orbit_index in range(fp3.period):
        manifold = wb.manifolds[(fp3, "unstable", orbit_index, 0)]
        # The geometric invariant -- this is what actually broke.
        assert_no_geometric_spikes(manifold)
        # cdist still orders the curve, only non-strictly (ties allowed).
        assert_cdist_monotonic(manifold, strict=False)
