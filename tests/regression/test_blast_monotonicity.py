"""Regression: blasting a resonance zone tripped the cdist-monotonicity guard.

Blasting iterates the bridges that live inside a resonance zone. The bug was in
the bridge endpoints: a bridge used to be bounded by fresh "straddle" points
offset 10% from each crossing but carrying the cdist of the crossing itself (see
``test_boundary_straddle_cdist``). Iterating such a bridge scaled that wrong cdist
by ``stretch_param``, and the merged image came out non-monotonic, tripping the
``assert sorted(cdists) == cdists`` guard in ``ManifoldMachine.iterate_manifold``.
With ``strict=True`` that ``AssertionError`` propagated out of ``blast_zone``
instead of being swallowed, so the blast failed on its very first step.

Bridges now reuse the real manifold nodes that bracket their crossings, so every
endpoint cdist describes the endpoint's own location and survives iteration. This
pins that: the blast runs its full two iterations with nothing skipped.
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
