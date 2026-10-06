"""Invariant-driven integration: invariants hold after *every* growth step.

Grows the manifolds one iteration at a time and re-checks the canonical
invariants after each step, so a regression that only appears at higher iteration
counts is caught at the exact step it is introduced. (The same invariants on
every finished law case, bridges included, run in
``tests/invariants/test_law_manifolds.py``.)
"""

from __future__ import annotations


from helpers.invariants import (
    assert_cdist_monotonic,
    assert_iterate_relation,
    assert_no_geometric_spikes,
    assert_one_to_one,
)


def test_unstable_invariants_hold_at_every_step(initialized):
    workbench, fp = initialized
    for step in range(1, 6):
        workbench.grow_n_times(fp, "unstable", num_iterations=1)
        manifold = workbench.manifolds[(fp, "unstable", 0, 0)]
        assert_cdist_monotonic(manifold)
        assert_no_geometric_spikes(manifold)
        assert_iterate_relation(manifold, rtol=1e-6)
        assert_one_to_one(manifold)


