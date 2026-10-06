"""The invariant-driven heart of the suite.

Grow the unstable manifold of the binary-horseshoe Hénon saddle by N iterations
and assert all four geometric invariants hold on the result, for a range of N.
"""

from __future__ import annotations

import pytest

from invariants import (
    assert_cdist_monotonic,
    assert_iterate_relation,
    assert_no_geometric_spikes,
    assert_one_to_one,
)


@pytest.mark.parametrize("n_iter", [1, 2, 4, 5])
def test_unstable_growth_preserves_invariants(initialized, n_iter):
    workbench, fp = initialized
    workbench.grow_n_times(fp, "unstable", num_iterations=n_iter)
    manifold = workbench.manifolds[(fp, "unstable", 0, 0)]

    assert_cdist_monotonic(manifold)
    assert_no_geometric_spikes(manifold)
    assert_iterate_relation(manifold, rtol=1e-6)
    assert_one_to_one(manifold)


@pytest.mark.parametrize("n_iter", [1, 2, 4])
def test_stable_growth_preserves_invariants(initialized, n_iter):
    workbench, fp = initialized
    workbench.grow_n_times(fp, "stable", num_iterations=n_iter)
    manifold = workbench.manifolds[(fp, "stable", 0, 0)]

    assert_cdist_monotonic(manifold)
    assert_no_geometric_spikes(manifold)
    assert_iterate_relation(manifold, rtol=1e-6)
    assert_one_to_one(manifold)


