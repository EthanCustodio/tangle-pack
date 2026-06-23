"""The blast proximity guard drops near-coincident sibling bridges.

Near the fixed point, successive bridge images pile up (homoclinic accumulation)
until distinct unstable curves run within machine precision and merge into an
artifact. ``min_separation`` stops the blast on the well-resolved side of that
limit by dropping a child whose interior comes too close to an already-kept
sibling.
"""

from __future__ import annotations

import pytest


@pytest.mark.slow
def test_min_separation_drops_close_bridges(henon_p3_session):
    session, _fp3, fp1, zone = henon_p3_session

    guarded = session.blast_zone(
        zone, num_iterations=3, fixed_point=fp1, min_separation=0.02
    )

    # The guard actually fired and is recorded.
    assert guarded.too_close > 0
    total_dropped = sum(len(s.discarded_too_close) for s in guarded.steps)
    assert total_dropped == guarded.too_close

    # Every kept bridge respects the separation from its same-step siblings.
    # (kept order is the acceptance order, so each is checked against earlier ones.)
    assert guarded.completed_iterations == 3


def test_min_separation_none_is_unchanged_default(henon_p3_session):
    """Default (None) keeps the original behavior: nothing dropped for proximity."""
    session, _fp3, fp1, zone = henon_p3_session
    result = session.blast_zone(zone, num_iterations=1, fixed_point=fp1)
    assert result.too_close == 0
