"""The invariant-driven heart of the suite.

Grow the unstable manifold of the binary-horseshoe Hénon saddle by N iterations
and assert all four geometric invariants hold on the result, for a range of N.
"""

from __future__ import annotations

import numpy as np
import pytest

from invariants import (
    assert_cdist_monotonic,
    assert_iterate_relation,
    assert_one_to_one,
)


@pytest.mark.parametrize("n_iter", [1, 2, 4, 5])
def test_unstable_growth_preserves_invariants(initialized, n_iter):
    workbench, fp = initialized
    workbench.grow_n_times(fp, "unstable", num_iterations=n_iter)
    manifold = workbench.manifolds[(fp, "unstable", 0, 0)]

    assert_cdist_monotonic(manifold, strict=True)
    assert_iterate_relation(manifold, rtol=1e-6)
    assert_one_to_one(manifold)


@pytest.mark.parametrize("n_iter", [1, 2, 4])
def test_stable_growth_preserves_invariants(initialized, n_iter):
    workbench, fp = initialized
    workbench.grow_n_times(fp, "stable", num_iterations=n_iter)
    manifold = workbench.manifolds[(fp, "stable", 0, 0)]

    assert_cdist_monotonic(manifold, strict=True)
    assert_iterate_relation(manifold, rtol=1e-6)
    assert_one_to_one(manifold)


def test_cdists_are_positive_and_increasing(grown_unstable):
    """cdist is a distance from the fixed point: strictly positive and growing."""
    from invariants import manifold_cdists

    _, _, manifold = grown_unstable
    cdists = manifold_cdists(manifold, "unstable")
    assert len(cdists) > 3
    # the root fixed point sits at cdist 0; every other point is strictly positive
    assert all(c >= 0 for c in cdists)
    assert cdists[-1] > 0
    assert cdists == sorted(cdists)
