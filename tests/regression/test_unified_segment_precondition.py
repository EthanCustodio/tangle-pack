"""Regression (PK-5 / DEC-20, WP7): the unified segment must be where f(S) stops.

The image of a bridge crosses the trellis's stable manifold only at the images
of its ends and in ``S \\ f(S)``; per stable branch ``S \\ f(S) = (c_s(f(end of
the previous branch)), end]``. The dual graph unifies the strong pip's window
instead, and the two agree exactly when the stable manifold is trimmed at the
pip orbit (every law case). Inversion x period 2 grown to 10 steps runs past the
pip orbit, and a walk there is unreachable; until the author settles DEC-20 the
dual graph states its precondition with one WARNING per branch where the two
differ. Rule-based: the branches are found independently (a crossing belongs to
``S \\ f(S)`` iff it has no registered preimage).
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

from cases import BUILDERS, build_inversion_period2
from helpers.laws import misplaced_unified_crossings
from helpers.logs import assert_logged

_LOGGER = "tanglepack.topology.DualGraph"


def _differing_branches(session) -> list:
    """The stable branches of a pip's fixed point where ``S \\ f(S)`` is not the unified
    span (read off the iterate table, not the dual graph's arithmetic)."""
    misplaced = misplaced_unified_crossings(session.dual_graph())
    return list(dict.fromkeys(key for key, _crossing in misplaced))


def test_a_stable_manifold_past_the_pip_orbit_is_warned(caplog):
    with np.errstate(over="ignore"):  # the 10th step's arms escape to infinity
        session = build_inversion_period2(unstable_steps=10, stable_steps=10).session
    differing = _differing_branches(session)
    assert differing, "the 10-step build no longer runs past the pip orbit"
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger=_LOGGER):
        session.dual_graph(rebuild=True)
    records = assert_logged(caplog, logging.WARNING, _LOGGER, count=len(differing))
    for key in differing:
        assert any(str(key[1:]) in record.getMessage() for record in records), key


@pytest.mark.parametrize("name", ["p3", "k10"])
def test_a_trellis_trimmed_at_the_pip_orbit_is_silent(caplog, name):
    session = BUILDERS[name]().session
    assert _differing_branches(session) == []
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger=_LOGGER):
        session.dual_graph(rebuild=True)
    assert_logged(caplog, logging.WARNING, _LOGGER, count=0)
