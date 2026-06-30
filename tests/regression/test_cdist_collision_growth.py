"""Regression: repeated growth must not scramble the manifold's geometry.

After a few growth iterations of the binary-horseshoe unstable manifold,
``stretch_param * parent_cdist`` can collapse two near-equal parent values onto
the same float. cdist is the manifold's ordering key, so a tie there used to let a
later ``merge_manifolds`` swap two geometrically distinct points and scramble the
linked list into zig-zag spikes (catastrophic at higher-stretch parameters such as
the period-3 Hénon tangle at k=2.1 -- see ``test_high_stretch_period3_growth``).

The fix keeps the merge's tie handling (a collided point is deduplicated, not
re-sorted) and lets refinement bridge any genuine high-stretch gap with
equal-cdist points. The invariant that matters is geometric: the grown manifold
stays a smooth simple curve with no node jutting off it. cdist remains
non-decreasing along the geometric ordering (ties allowed).
"""

from __future__ import annotations

import pytest

from invariants import assert_cdist_monotonic, assert_no_geometric_spikes


@pytest.mark.regression
def test_growth_keeps_geometry_smooth(initialized):
    workbench, fp = initialized
    workbench.grow_n_times(fp, "unstable", num_iterations=4)
    manifold = workbench.manifolds[(fp, "unstable", 0, 0)]
    assert_no_geometric_spikes(manifold)
    assert_cdist_monotonic(manifold, strict=False)
