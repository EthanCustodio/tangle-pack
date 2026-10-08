"""The solver refuses an orbit whose minimal period is a proper divisor of the request.

``FixedPointSolver.minimal_period`` is the smallest shift under which a cycle
repeats itself; ``construct_fixed_point`` raises ``ValueError`` when a solve
slides to such a lower-period orbit (the stress campaign's E8 found the
k=10 inversion saddle accepted as a "period-2 saddle" with a lost inversion).
Rules only: no message wording is pinned.
"""

from __future__ import annotations

import numpy as np
import pytest

import tanglepack
from tanglepack.examples import (
    HENON_K10,
    HENON_P3,
    henon_jacobian,
    henon_map,
    henon_map_inverse,
    saddle_guesses,
)
from tanglepack.numerics import FixedPointSolver

minimal_period = FixedPointSolver.minimal_period


def _workbench(params: tuple[float, float]) -> tanglepack.TangleWorkbench:
    return tanglepack.TangleWorkbench(
        henon_map(*params), henon_map_inverse(*params), henon_jacobian(*params)
    )


def _cycle(params: tuple[float, float], seed: list[float], period: int) -> np.ndarray:
    """``period`` points of the orbit of ``seed``, in map order."""
    points = [np.asarray(seed, dtype=float)]
    for _ in range(period - 1):
        points.append(np.ravel(henon_map(*params)(points[-1]))[:2])
    return np.array(points)


@pytest.mark.parametrize("base", [[[4.0, -4.0]], [[1.5, -0.5], [-2.0, 3.0]]])
@pytest.mark.parametrize("copies", [1, 2, 3, 6])
def test_tiled_orbit_has_its_base_period(base, copies):
    assert minimal_period(np.tile(base, (copies, 1))) == len(base)


def test_genuine_three_cycle_is_period_three():
    cycle = _cycle(HENON_P3, [0.0, 1.0], 3)
    assert minimal_period(cycle) == 3
    assert minimal_period(np.tile(cycle, (2, 1))) == 3


def test_round_off_jitter_keeps_the_base_period():
    rng = np.random.default_rng(0)
    tiled = np.tile(_cycle(HENON_P3, [0.0, 1.0], 3), (4, 1))
    assert minimal_period(tiled + 1e-13 * rng.standard_normal(tiled.shape)) == 3


def test_a_real_perturbation_is_full_length():
    tiled = np.tile(_cycle(HENON_P3, [0.0, 1.0], 3), (2, 1))
    tiled[4] += 1e-3
    assert minimal_period(tiled) == 6


@pytest.mark.parametrize("name", ["saddle", "inversion"])
@pytest.mark.parametrize("copies", [2, 3, 4, 6])
def test_period_one_saddle_tiled_is_refused(name, copies):
    """A period-1 saddle (with or without inversion) asked at a higher period."""
    workbench = _workbench(HENON_K10)
    point = workbench.construct_fixed_point(saddle_guesses(*HENON_K10)[name])
    with pytest.raises(ValueError):
        workbench.construct_fixed_point(np.tile(point.coordinates, (copies, 1)))
    assert len(workbench.fixed_points) == 1


def test_period_three_orbit_tiled_to_six_is_refused():
    """The slide that would hit the p3 fixture orbit."""
    workbench = _workbench(HENON_P3)
    orbit = workbench.construct_fixed_point(saddle_guesses(*HENON_P3)["period_3"])
    with pytest.raises(ValueError):
        workbench.construct_fixed_point(np.tile(orbit.coordinates, (2, 1)))


def test_inversion_saddle_requested_at_its_k_value_is_refused():
    """``k_value = 2`` is not a period: the inversion saddle is period 1."""
    workbench = _workbench(HENON_K10)
    point = workbench.construct_fixed_point(saddle_guesses(*HENON_K10)["inversion"])
    assert point.period == 1 and point.k_value == 2
    with pytest.raises(ValueError):
        workbench.construct_fixed_point(np.tile(point.coordinates, (point.k_value, 1)))
