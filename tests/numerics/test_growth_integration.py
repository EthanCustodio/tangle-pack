"""Invariant-driven integration: invariants hold after *every* growth step.

Grows the manifolds one iteration at a time and re-checks the canonical
invariants after each step, so a regression that only appears at higher iteration
counts is caught at the exact step it is introduced.
"""

from __future__ import annotations

import pytest

from invariants import (
    assert_cdist_monotonic,
    assert_iterate_relation,
    assert_one_to_one,
)


@pytest.mark.slow
def test_unstable_invariants_hold_at_every_step(initialized):
    workbench, fp = initialized
    for step in range(1, 6):
        workbench.grow_n_times(fp, "unstable", num_iterations=1)
        manifold = workbench.manifolds[(fp, "unstable", 0, 0)]
        assert_cdist_monotonic(manifold, strict=True)
        assert_iterate_relation(manifold, rtol=1e-6)
        assert_one_to_one(manifold)


@pytest.mark.slow
def test_bridges_individually_satisfy_invariants(small_tangle):
    workbench, fp = small_tangle
    workbench.trim_stable_manifolds(fp)
    bridges = workbench.create_bridges(fp)
    assert bridges, "fixture produced no bridges"

    for bridge in bridges:
        assert_cdist_monotonic(bridge, strict=False)
        assert_one_to_one(bridge)
