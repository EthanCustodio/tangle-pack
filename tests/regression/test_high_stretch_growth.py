"""Regression: growth must not scramble a manifold's geometry (cdist-strict memory).

cdist is the manifold's ordering key, and it can saturate: two geometrically
distinct points can end up with cdists equal, or one float ULP apart. A tie
there used to let a later ``merge_manifolds`` keep both points and order them by
a meaningless key, scrambling the linked list (the "super spiky inner tangle" of
the period-3 Hénon tangle at k=2.1).

An earlier attempt forced cdist to be strictly injective (nudging every tie up
by ``nextafter``) and stopped refinement wherever the cdist midpoint was not
representable; that made the cdist VALUES strictly increasing while leaving the
GEOMETRY scrambled. The fix keeps the merge's tie handling (a collided point is
deduplicated, not re-sorted) and lets refinement bridge a gap with equal-cdist
points. The invariant that matters is therefore geometric: the grown manifold
stays a smooth curve, with no node jutting off it and no hairpin, and cdist
stays non-decreasing (ties allowed).

Dev Notes:
    How the period-3 run reaches the near-ULP regime (re-derived 2026-10-05).
    With the default seed (``max(fixed_point.accuracy, 5e-6)``, the fsolve
    orbit making ``accuracy`` tiny) the k=2.1 arms never get there: adjacent
    relative cdist gaps stay above ~1e11 eps until the unstable arm escapes
    (step 10, ~1.5e6 points). Growth scales a gap and its cdist by the same
    stretch factor, so only refinement shrinks the relative gap, one halving per
    layer, and the seeded curve escapes long before ~52 halvings. A coarser seed
    (``accuracy = 5e-3``) does reach it in five stable steps: the fundamental
    segment is then curved enough that the images of its two ends land ~1e-5
    apart in space while their cdists agree to one ULP. That is the same
    collision the original bug mishandled (two distinct points the cdist can no
    longer order), produced at the seed seam instead of at a high-stretch fold.
    Evidence it has teeth: re-introducing the strictify + ``representable``
    guard (nudge tied image cdists and a merge collision apart by ``nextafter``,
    refine only where the cdist midpoint is strictly between its neighbours)
    leaves 20-30 hairpin reversals per branch here, while the current code
    leaves none; no segment is long enough to trip the spike check, which is
    why the reversal check is the one with teeth. See
    docs/test_suite_refactor_ledger.md (Author follow-up).
"""

from __future__ import annotations

import numpy as np

from cases import P3_ORBIT_SEED
from helpers.invariants import (
    assert_cdist_monotonic,
    assert_no_geometric_spikes,
    assert_no_reversals,
    manifold_cdists,
)
from tanglepack import TangleWorkbench
from tanglepack.examples import henon_jacobian, henon_map, henon_map_inverse

#: Seed step of the period-3 run: ``ManifoldInitializer`` seeds the fundamental
#: segment ``max(fixed_point.accuracy, 5e-6)`` along the eigenvector.
P3_SEED_STEP = 5e-3

#: Non-vacuity bound: some adjacent cdist pair must sit within this many float
#: epsilons of each other relative to its magnitude (one ULP is 0.5..1 eps).
NEAR_ULP_EPS = 4.0


def test_k10_growth_keeps_geometry_smooth(initialized) -> None:
    """Four unstable steps of the k=10 horseshoe: no spike, no hairpin, cdist
    non-decreasing."""
    workbench, fp = initialized
    workbench.grow_n_times(fp, "unstable", num_iterations=4)
    manifold = workbench.manifolds[(fp, "unstable", 0, 0)]
    assert_no_geometric_spikes(manifold)
    assert_no_reversals(manifold)
    assert_cdist_monotonic(manifold, strict=False)


def test_period3_near_ulp_growth_is_not_scrambled() -> None:
    """Five stable steps of the k=2.1 period-3 orbit from a coarse seed reach
    adjacent cdists one ULP apart on every branch, and every branch stays a
    smooth, ordered curve."""
    workbench = TangleWorkbench(henon_map(2.1), henon_map_inverse(2.1), henon_jacobian(2.1))
    workbench._man_machine.area_cutoff = 1e-7
    fp3 = workbench.construct_fixed_point([list(point) for point in P3_ORBIT_SEED])
    workbench.orient_eigenvectors(
        fp3, {"unstable": np.array([0, 1]), "stable": np.array([-1, -1])}
    )
    fp3.accuracy = P3_SEED_STEP
    workbench.initialize_both_manifolds(fp3)
    workbench.grow_n_times(fp3, "stable", num_iterations=5)

    eps = np.finfo(np.float64).eps
    for orbit_index in range(fp3.period):
        manifold = workbench.manifolds[(fp3, "stable", orbit_index, 0)]
        cdists = np.asarray(manifold_cdists(manifold, "stable"), dtype=float)
        gaps = np.diff(cdists)
        upper = cdists[1:]
        relative = gaps[upper > 0] / upper[upper > 0]
        assert relative.min() < NEAR_ULP_EPS * eps, (
            f"vacuous: branch {orbit_index}'s closest adjacent cdists are "
            f"{relative.min() / eps:.3g} eps apart, never near one ULP"
        )
        assert_no_geometric_spikes(manifold)
        assert_no_reversals(manifold)
        assert_cdist_monotonic(manifold, strict=False)
