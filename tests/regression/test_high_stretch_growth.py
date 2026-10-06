"""Regression: growth must not scramble a manifold's geometry (cdist-strict memory).

After a few growth iterations ``stretch_param * parent_cdist`` can collapse two
near-equal parent values onto the same float. cdist is the manifold's ordering
key, so a tie there used to let a later ``merge_manifolds`` swap two
geometrically distinct points and scramble the linked list into zig-zag spikes:
the "super spiky inner tangle", catastrophic at high stretch such as the
period-3 Hénon tangle at k=2.1 (per-step eigenvalue factor ~4), where cdists
near the saddle fall to floating-point spacing.

An earlier attempt forced cdist to be strictly injective (nudging every tie up
by ``nextafter``); that made the cdist VALUES strictly increasing while leaving
the GEOMETRY scrambled. The fix keeps the merge's tie handling (a collided point
is deduplicated, not re-sorted) and lets refinement bridge a genuine high-stretch
gap with equal-cdist points. The invariant that matters is therefore geometric:
the grown manifold stays a smooth simple curve with no node jutting off it, and
cdist stays non-decreasing (ties allowed).

Dev Notes:
    The period-3 run is VACUOUS on the current code (Phase 7a finding): four
    unstable steps leave 7 nodes per branch with adjacent cdist gaps of ~3e15
    ULP, so it never reaches the near-ULP regime it guards and a non-vacuity
    check would fail. Its parameters need re-deriving (author item); until then
    it stays as a cheap smoke of the same merge/refine path.
"""

from __future__ import annotations

import numpy as np

from cases import P3_ORBIT_SEED
from helpers.invariants import assert_cdist_monotonic, assert_no_geometric_spikes
from tanglepack import TangleWorkbench
from tanglepack.examples import henon_jacobian, henon_map, henon_map_inverse


def test_k10_growth_keeps_geometry_smooth(initialized) -> None:
    """Four unstable steps of the k=10 horseshoe: no spike, cdist non-decreasing."""
    workbench, fp = initialized
    workbench.grow_n_times(fp, "unstable", num_iterations=4)
    manifold = workbench.manifolds[(fp, "unstable", 0, 0)]
    assert_no_geometric_spikes(manifold)
    assert_cdist_monotonic(manifold, strict=False)


def test_period3_high_stretch_growth_is_not_scrambled() -> None:
    """Four unstable steps of the k=2.1 period-3 orbit (the buggy strictify build
    left 11/7/37 spikes here): every branch stays smooth and ordered."""
    workbench = TangleWorkbench(henon_map(2.1), henon_map_inverse(2.1), henon_jacobian(2.1))
    workbench._man_machine.area_cutoff = 1e-7
    fp3 = workbench.construct_fixed_point([list(point) for point in P3_ORBIT_SEED])
    workbench.orient_eigenvectors(
        fp3, {"unstable": np.array([0, 1]), "stable": np.array([-1, -1])}
    )
    workbench.initialize_both_manifolds(fp3)
    # The arm escapes by the fifth iterate, so stop at four.
    workbench.grow_n_times(fp3, "unstable", num_iterations=4)

    for orbit_index in range(fp3.period):
        manifold = workbench.manifolds[(fp3, "unstable", orbit_index, 0)]
        assert_no_geometric_spikes(manifold)
        assert_cdist_monotonic(manifold, strict=False)
