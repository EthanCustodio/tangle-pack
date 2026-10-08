"""
Regression (2026-10-07, numerics stress campaign E7/E9): the resolution
WARNING fires exactly where an independent oracle finds a resolution problem.

Near a tangency the polylines register crossings that cannot be told apart
from polygon noise, silently until WP1. The rule, on any build: the
workbench logs a WARNING if and only if the oracle finds something -- a
same-stability crossing in the tangle box (the test-tier law
``helpers.laws.check_no_same_stability_crossing``, which does not use the
library diagnostic) or a registered crossing whose resolution ratio is below
the workbench's threshold. Parametrised over the near-tangency builder (a
positive case) and well-resolved law cases (k10, and p3, whose margin over the
threshold is the thinnest). No count is pinned.
"""

from __future__ import annotations

import logging

import pytest

from cases import BUILDERS
from helpers import laws

WORKBENCH = "tanglepack.numerics.TangleWorkbench"


def _oracle_finds_something(case) -> bool:
    """A same-stability crossing in the box, or a ratio below the threshold."""
    threshold = case.workbench.resolution_ratio_threshold
    if any(
        ix.resolution_ratio is not None and ix.resolution_ratio < threshold
        for _iid, ix in case.registry
    ):
        return True
    try:
        laws.check_no_same_stability_crossing(case)
    except AssertionError:
        return True
    return False


@pytest.mark.parametrize("name", ["near_tangency", "k10", "p3"])
def test_resolution_warning_iff_the_oracle_finds_a_problem(name: str, caplog):
    with caplog.at_level(logging.WARNING, logger=WORKBENCH):
        case = BUILDERS[name]()
    warned = any(
        record.name == WORKBENCH and record.levelno == logging.WARNING
        for record in caplog.records
    )
    assert warned == _oracle_finds_something(case)
