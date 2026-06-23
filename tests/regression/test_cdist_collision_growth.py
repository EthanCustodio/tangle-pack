"""Regression: growth produces colliding cdists.

Canonical distance must be injective along a manifold -- two distinct points may
never share a cdist (a collision signals an unhandled edge case). After only a few
growth iterations of the binary-horseshoe unstable manifold, two geometrically
distinct, adjacent points acquire cdists equal to within machine epsilon. This is
the cdist-assignment degeneracy that also drives the runaway refinement at higher
iteration counts.

Expected fix: the refinement / iterate cdist assignment must keep distinct points'
cdists distinct. Until then this xfails.
"""

from __future__ import annotations

import pytest

from invariants import assert_no_cdist_collision


@pytest.mark.regression
@pytest.mark.xfail(strict=True, reason="adjacent points acquire identical cdist after ~4 iterations")
def test_growth_keeps_cdists_distinct(initialized):
    workbench, fp = initialized
    workbench.grow_n_times(fp, "unstable", num_iterations=4)
    manifold = workbench.manifolds[(fp, "unstable", 0, 0)]
    assert_no_cdist_collision(manifold, atol=1e-12)
