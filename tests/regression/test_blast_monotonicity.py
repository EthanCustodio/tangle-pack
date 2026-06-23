"""Regression: blasting a resonance zone violates the cdist-monotonicity guard.

Blasting iterates the bridges that live inside a resonance zone. Iterating an
interior bridge forward should keep its image's cdists monotonic, but the wrong
cdist assigned to the bridge's straddle endpoints (see
``test_boundary_straddle_cdist``) plus the refinement degeneracy make the merged
image non-monotonic, tripping the ``assert sorted(cdists) == cdists`` guard in
``ManifoldMachine.iterate_manifold``. With ``strict=True`` that ``AssertionError``
propagates out of ``blast_zone`` instead of being swallowed.

This reproduces at the very first blast step. The fix is to assign straddle cdists
at their true location and keep refinement well-behaved; then the blast completes
with nothing skipped.
"""

from __future__ import annotations

import pytest


@pytest.mark.slow
@pytest.mark.regression
def test_blast_completes_without_monotonicity_failure(henon_p3_session):
    session, _fp3, fp1, inner_zone = henon_p3_session
    result = session.blast_zone(
        inner_zone, num_iterations=2, fixed_point=fp1, strict=True
    )
    assert result.skipped == 0
    assert result.completed_iterations == 2
