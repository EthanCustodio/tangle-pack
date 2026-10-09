"""The workbench's API traps from the 2026-10-07 stress report, made loud.

* A ``grow_until_*`` driver at its cap raises :class:`MaxIterationsReached`, a
  ``ValueError`` carrying the cap, the stretch of one growth round and a
  suggested cap from ``1 / ln(stretch)`` -- and the suggestion is enough where
  the default cap of 10 is not (a weakly hyperbolic Hénon saddle).
* ``construct_fixed_point`` logs ONE WARNING when the Jacobian says the map
  does not preserve area (``|det J| != 1``), and none for ``b = 1``; the
  ``b = -1`` determinant passes the check because ``|det J|`` is compared.
* ``MaxIterationsReached`` and ``TURNAROUND_LOG_STRETCH`` are exported, and
  ``PointBudgetExceeded.blast_context`` always exists (WP11, API-3).

The stale-registry trap (``compute_intersections`` returns the registry it
installs) is tested with the other registry rules in
``test_single_source_of_truth.py``. Log tests check level and logger only.
"""

from __future__ import annotations

import logging
import math

import numpy as np
import pytest

from helpers.logs import assert_logged
from tanglepack.examples.henon import henon_jacobian, henon_map, henon_map_inverse
from tanglepack.numerics.DynamicalSystem import DynamicalSystem
from tanglepack.numerics.TangleWorkbench import (
    TURNAROUND_LOG_STRETCH,
    MaxIterationsReached,
    TangleWorkbench,
)

WORKBENCH_LOGGER = "tanglepack.numerics.TangleWorkbench"


def _henon_saddle(k: float, b: float = 1.0) -> tuple[TangleWorkbench, object]:
    """A workbench on Hénon ``(k, b)`` with its period-1 saddle, manifolds seeded."""
    workbench = TangleWorkbench(henon_map(k, b), henon_map_inverse(k, b), henon_jacobian(k, b))
    x = ((1 + b) + math.sqrt((1 + b) ** 2 + 4 * k)) / 2  # x**2 - (1+b) x - k = 0
    fp = workbench.construct_fixed_point([x, -b * x])
    workbench.orient_eigenvectors(fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])})
    workbench.initialize_both_manifolds(fp)
    return workbench, fp


# --------------------------------------------------------------------------- #
# the cap error
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    ("stability", "drive"),
    [
        ("stable", lambda wb, fp: wb.grow_until_turnaround(fp, "stable", max_iterations=1)),
        ("unstable", lambda wb, fp: wb.grow_until_intersection(fp, "unstable", max_iterations=1)),
        (
            "unstable",
            lambda wb, fp: wb.grow_until_arclength(fp, "unstable", length=1e9, max_iterations=2),
        ),
    ],
    ids=["turnaround", "intersection", "arclength"],
)
def test_the_cap_error_carries_a_suggested_cap(initialized, stability, drive):
    """The cap error is a ``ValueError`` whose suggestion follows the grown branch's stretch."""
    workbench, fp = initialized

    with pytest.raises(MaxIterationsReached) as info:
        drive(workbench, fp)

    error = info.value
    assert isinstance(error, ValueError)
    beta = fp.per_step_beta(stability)
    assert error.stretch_per_round == pytest.approx(beta if stability == "unstable" else 1 / beta)
    assert error.stretch_per_round > 1
    assert error.suggested_max_iterations > error.max_iterations


def test_the_suggested_cap_rescues_a_weakly_hyperbolic_saddle():
    """Where the default cap fails, a fresh run at the suggested cap turns around.

    Hénon at k = 0.3 stretches by about 4 per round, too little for a turnaround
    in 10 rounds (the stress report's small-k ``Max iterations`` class).
    """
    workbench, fp = _henon_saddle(0.3)
    with pytest.raises(MaxIterationsReached) as info:
        workbench.grow_until_turnaround(fp, "stable", max_iterations=10)
    suggested = info.value.suggested_max_iterations
    assert suggested == math.ceil(TURNAROUND_LOG_STRETCH / math.log(info.value.stretch_per_round))

    workbench, fp = _henon_saddle(0.3)
    workbench.grow_until_turnaround(fp, "stable", max_iterations=suggested)


def test_default_cap_scales_with_the_stretch():
    """NUM-15: with no cap given, the drivers take the rounds the stretch needs
    (never fewer than 10): the k = 0.3 saddle that fails at a cap of 10 turns around."""
    workbench, fp = _henon_saddle(0.3)
    stretch = 1.0 / fp.per_step_beta("stable")
    assert math.ceil(TURNAROUND_LOG_STRETCH / math.log(stretch)) > 10
    workbench.grow_until_turnaround(fp, "stable")


# --------------------------------------------------------------------------- #
# area preservation
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("b", [1.0, -1.0])
def test_an_area_preserving_jacobian_passes_in_either_orientation(b):
    """``|det J|`` is compared, so ``det J = -1`` (orientation reversing) passes."""
    system = DynamicalSystem(henon_map(4, b), henon_map_inverse(4, b), henon_jacobian(4, b))
    assert system.check_area_preserving(np.array([2.0, -2.0])) is True


def test_without_a_jacobian_there_is_nothing_to_check():
    """No Jacobian: the check answers None and the solver's fallback is not judged."""
    system = DynamicalSystem(henon_map(4, 0.5), henon_map_inverse(4, 0.5))
    assert system.check_area_preserving(np.array([2.0, -2.0])) is None


def test_a_dissipative_map_warns_once_at_construction(caplog):
    """A ``b = 0.5`` Hénon saddle is still constructed, with ONE WARNING."""
    caplog.set_level(logging.WARNING, logger=WORKBENCH_LOGGER)
    workbench, fp = _henon_saddle(10, b=0.5)
    assert fp in workbench.fixed_points
    assert_logged(caplog, "WARNING", WORKBENCH_LOGGER, count=1)


def test_an_area_preserving_map_constructs_silently(caplog):
    """``b = 1``: no WARNING from the workbench while constructing and seeding."""
    caplog.set_level(logging.WARNING, logger=WORKBENCH_LOGGER)
    _henon_saddle(10)
    assert not [r for r in caplog.records if r.name == WORKBENCH_LOGGER]


# --------------------------------------------------------------------------- #
# WP11 (API-3): the growth-cap symbols are exported, the blast context declared
# --------------------------------------------------------------------------- #
def test_the_growth_cap_symbols_are_exported():
    import tanglepack

    assert tanglepack.MaxIterationsReached is MaxIterationsReached
    assert tanglepack.TURNAROUND_LOG_STRETCH == TURNAROUND_LOG_STRETCH


def test_the_blast_context_is_declared():
    from tanglepack import PointBudgetExceeded

    assert PointBudgetExceeded(1, 2, "map", 0.0).blast_context is None
