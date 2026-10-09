"""``side_after``: the one orientation rule for carrying a side along an orbit (C3, WP11)."""

from __future__ import annotations

import itertools

import pytest

from tanglepack.topology import side_after
from tanglepack.topology.TopologyResults import OPPOSITE_SIDE

STEPS = range(-4, 5)


@pytest.mark.parametrize("side", ["left", "right"])
def test_a_preserving_map_keeps_the_side(side):
    assert all(side_after(side, n, True) == side for n in STEPS)


@pytest.mark.parametrize("side", ["left", "right"])
def test_a_reversing_map_flips_the_side_every_step(side):
    for n in STEPS:
        assert side_after(side, n, False) == (side if n % 2 == 0 else OPPOSITE_SIDE[side])


@pytest.mark.parametrize("side, preserving", itertools.product(["left", "right"], [True, False]))
def test_steps_compose(side, preserving):
    for a, b in itertools.product(STEPS, STEPS):
        assert side_after(side_after(side, a, preserving), b, preserving) == side_after(
            side, a + b, preserving
        )
