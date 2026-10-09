"""The opt-in point budget of ``blast_zone`` and its stop reasons.

``blast_zone(max_points=...)`` stops BETWEEN iterations, before one predicted to
take the blast over budget, on a state that can be blasted again; the machine's
per-call cap (``max_points_per_step``) is the backstop that raises mid-image.
``BlastResult.terminated_early`` is the ``"frontier_empty"`` stop reason, so the
two can never disagree. No point count, factor or iteration count of the fixture
is pinned: only the relations the budget promises.
"""

from __future__ import annotations

import logging

import pytest

import cases
from helpers import laws
from tanglepack import PointBudgetExceeded

BLAST = {"min_separation": 1e-5}


def _k28_zone():
    """The k=2.8 case at the zones stage, unblasted, and its one zone."""
    case = cases.build_k28(blasts=0, through="zones")
    return case, case.zones[0]


def _blast(case, zone, num_iterations, **kwargs):
    """Blast ``zone`` over the case's fixed points."""
    return case.session.blast_zone(
        zone, num_iterations, fixed_point=case.fixed_points, **BLAST, **kwargs
    )


def _assert_bridges_consistent(case):
    """The bridge laws still hold after a stopped blast."""
    laws.check_bridge_partial_iff_no_id(case)
    laws.check_bridge_endpoints_on_own_branch(case)
    laws.check_bridges_do_not_overlap(case)


def test_budget_stops_between_iterations_and_resumes():
    """A budget just above the first blast's cost stops a longer second blast
    cleanly, and the session can be blasted again afterwards."""
    case, zone = _k28_zone()
    first = _blast(case, zone, 1)
    assert first.completed_iterations == 1 and first.points_created > 0

    num_iterations = 4
    second = _blast(case, zone, num_iterations, max_points=first.points_created + 1)
    assert second.stop_reason == "point_budget"
    assert not second.terminated_early
    assert second.completed_iterations < num_iterations
    assert len(second.steps) == second.completed_iterations
    _assert_bridges_consistent(case)

    third = _blast(case, zone, 1)
    assert third.completed_iterations == 1 or third.stop_reason == "frontier_empty"
    _assert_bridges_consistent(case)


def test_stored_factor_predicts_the_next_call():
    """A zone blasted before is predicted from its stored factor (a one-point
    budget refuses the first iteration); a fresh zone has no factor, so its
    first iteration runs and the budget stops the one after."""
    case, zone = _k28_zone()
    fresh = _blast(case, zone, 3, max_points=1)
    assert fresh.completed_iterations == 1
    assert fresh.points_created > 0
    assert fresh.stop_reason == "point_budget"

    again = _blast(case, zone, 3, max_points=1)
    assert again.completed_iterations == 0
    assert again.stop_reason == "point_budget"


def test_without_prediction_the_budget_stops_once_reached():
    """``predict=False`` ignores the stored factor and stops once spent."""
    case, zone = _k28_zone()
    _blast(case, zone, 1)
    result = _blast(case, zone, 3, max_points=1, predict=False)
    assert result.completed_iterations == 1
    assert result.points_created >= 1
    assert result.stop_reason == "point_budget"


def test_stop_reason_without_budget():
    """No budget never stops on points, and ``terminated_early`` IS the
    ``"frontier_empty"`` reason."""
    case, zone = _k28_zone()
    for result in (_blast(case, zone, 2), case.session.blast_zone(zone, 2, fixed_point=[])):
        assert result.stop_reason in {"max_iterations", "frontier_empty"}
        assert result.terminated_early == (result.stop_reason == "frontier_empty")
    assert result.stop_reason == "frontier_empty"  # an empty fixed-point filter


def test_machine_budget_propagates_with_context(caplog):
    """A ``PointBudgetExceeded`` from one bridge image is never a skipped bridge:
    it propagates carrying the blast's context."""
    case, zone = _k28_zone()
    case.session.max_points_per_step = 1
    with caplog.at_level(logging.WARNING, logger="tanglepack.loom.Blast"):
        with pytest.raises(PointBudgetExceeded) as raised:
            _blast(case, zone, 1)
    exc = raised.value
    assert exc.points_created <= exc.budget == 1
    zone_key, iteration, _bridge_id = exc.blast_context
    assert zone_key == zone.key and iteration == 1
    assert not [r for r in caplog.records if r.name == "tanglepack.loom.Blast"]
