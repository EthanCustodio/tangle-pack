"""The per-step point budget, growth stats and the public refinement settings.

Rules pinned here (WP3 of the 2026-10-07 stress-report plan):

* a growth step (one ``k_value`` chain) or one bridge image that would create
  more points than ``max_points_per_step`` raises ``PointBudgetExceeded``
  before the overrunning map call (``points_created <= budget``, nothing
  evaluated or added past it), with the generation advanced, every bridge's
  version bumped and the partly grown curve still ordered; ``None``
  disables the budget;
* the count is per top-level call: consecutive steps or blasts that each fit
  never raise, however much they add up to;
* a bridge abandoned half-refined by an overrun refuses a second iterate
  instead of handing back its half-built image;
* every ``grow_until_*`` driver propagates the overrun;
* ``growth_stats`` (per branch), the prediction hint and ``escape_radius`` WARNINGs, and
  ``grow_until_bounded``;
* ``area_cutoff`` reaches the ONE shared machine from the workbench and the
  session, constructor and setter alike.

Budgets are derived from a probe build of the same case, never hard-coded
counts (tests/conftest.py policy), and the probe measures INDEPENDENTLY of the
counter under test: the change in the manifold's own point count for a growth
step, and the points the dynamical system evaluates for a bridge image (every
created point is one map evaluation). A broken counter therefore cannot skew
the budget a test sets.
"""

from __future__ import annotations

import logging
import time

import numpy as np
import pytest

from cases import (
    HENON_K28,
    K28_SADDLE_SEED,
    PERIOD1_ORIENTATION,
    _orientation,
    build_k10,
    build_k28,
)
from helpers.invariants import assert_cdist_monotonic, assert_one_to_one
from helpers.logs import assert_logged
from tanglepack import PointBudgetExceeded, TangleSession
from tanglepack.examples.henon import henon_jacobian, henon_map, henon_map_inverse

WORKBENCH_LOGGER = "tanglepack.numerics.TangleWorkbench"
STEPS = 8


def _seeded():
    """A k=10 workbench with its fundamental segments, and its saddle."""
    case = build_k10(through="seeded")
    return case.workbench, case.fixed_point


def _point_count(workbench, fp, stability: str = "unstable") -> int:
    """Points over every branch of one manifold."""
    return sum(
        len(manifold.get_point_array())
        for key, manifold in workbench.manifolds.items()
        if key[0] is fp and key[1] == stability
    )


def _awaiting_map(workbench, fp) -> int:
    """Unstable points without an image yet: what the next step's map stage maps."""
    return sum(
        len(manifold.get_non_iterated_point_array())
        for key, manifold in workbench.manifolds.items()
        if key[0] is fp and key[1] == "unstable"
    )


def _created_per_step(steps: int) -> list[int]:
    """Points each unstable growth step adds to a fresh k=10 build."""
    workbench, fp = _seeded()
    created = []
    for _ in range(steps):
        before = _point_count(workbench, fp)
        workbench.grow_n_times(fp, "unstable", 1)
        created.append(_point_count(workbench, fp) - before)
    return created


def _count_evaluations(workbench) -> list[int]:
    """Spy on the system's batched maps; returns the list of batch sizes."""
    batches: list[int] = []
    system = workbench.dynamical_system
    for name in ("map_batch", "map_inv_batch"):
        batched = getattr(system, name)

        def counted(coords, batched=batched):
            batches.append(len(coords))
            return batched(coords)

        setattr(system, name, counted)
    return batches


# --------------------------------------------------------------------------- #
# growth
# --------------------------------------------------------------------------- #
def test_overrunning_step_raises_in_refinement_with_state_consistent() -> None:
    """A budget that covers the last step's map stage but not its refinement
    raises in ``refine``, before evaluating or adding a point past the budget;
    the generation moves and the partly grown manifold is still an ordered,
    acyclic curve."""
    created = _created_per_step(STEPS)
    assert created[-1] > max(created[:-1]), "the probe needs a step bigger than all before"

    workbench, fp = _seeded()
    workbench.grow_n_times(fp, "unstable", STEPS - 1)
    budget = _awaiting_map(workbench, fp)
    assert max(created[:-1]) <= budget < created[-1], "the map stage must fit, the step not"
    workbench.max_points_per_step = budget
    evaluated = _count_evaluations(workbench)
    before = _point_count(workbench, fp)
    generation = workbench.generation

    with pytest.raises(PointBudgetExceeded) as raised:
        workbench.grow_n_times(fp, "unstable", 1)

    assert raised.value.budget == budget
    assert raised.value.stage == "refine"
    assert raised.value.points_created <= budget
    assert sum(evaluated) <= budget
    assert _point_count(workbench, fp) - before <= budget
    assert workbench.generation != generation
    for key, manifold in workbench.manifolds.items():
        if key[1] == "unstable":
            assert_cdist_monotonic(manifold)
            assert_one_to_one(manifold)


def test_map_stage_refuses_before_evaluating() -> None:
    """A budget below the map stage refuses the step before a single map
    evaluation."""
    workbench, fp = _seeded()
    workbench.max_points_per_step = 0
    evaluated = _count_evaluations(workbench)

    with pytest.raises(PointBudgetExceeded) as raised:
        workbench.grow_n_times(fp, "unstable", 1)

    assert raised.value.stage == "map"
    assert not evaluated


def test_overrun_drops_every_bridge_walk() -> None:
    """Refinement lays points inside already-cut bridges before it overruns,
    so the overrun still bumps every bridge's version."""
    case = build_k10(through="bridges")
    workbench, fp = case.workbench, case.fixed_point
    assert workbench.bridges, "the probe needs cut bridges"
    versions = [bridge.version for bridge in workbench.bridges]
    workbench.max_points_per_step = _awaiting_map(workbench, fp)

    with pytest.raises(PointBudgetExceeded) as raised:
        workbench.grow_n_times(fp, "unstable", 1)

    assert raised.value.stage == "refine"
    assert all(b.version > v for b, v in zip(workbench.bridges, versions))


def test_no_budget_never_raises() -> None:
    """``max_points_per_step = None`` grows what a tiny budget refuses."""
    workbench, fp = _seeded()
    workbench.max_points_per_step = None
    before = _point_count(workbench, fp)
    workbench.grow_n_times(fp, "unstable", STEPS)
    assert _point_count(workbench, fp) > before


def test_consecutive_steps_each_within_budget_never_raise() -> None:
    """The count restarts with every step: a budget of the largest single step
    survives all of them although their sum exceeds it."""
    created = _created_per_step(STEPS)
    assert sum(created) > max(created)

    workbench, fp = _seeded()
    workbench.max_points_per_step = max(created)
    workbench.grow_n_times(fp, "unstable", STEPS)


def _bounded(workbench, fp):
    workbench.grow_until_bounded(fp, "unstable", radius=1e300, max_iterations=2)


def _turnaround(workbench, fp):
    workbench.grow_until_turnaround(fp, "stable", max_iterations=2)


def _arclength(workbench, fp):
    workbench.grow_until_arclength(fp, "unstable", length=1e9, max_iterations=2)


def _intersection(workbench, fp):
    workbench.grow_until_intersection(fp, "unstable", max_iterations=2)


def _predicate(workbench, fp):
    workbench.grow_until(fp, lambda wb: False, max_iterations=2)


@pytest.mark.parametrize(
    "drive",
    [_bounded, _turnaround, _arclength, _intersection, _predicate],
    ids=["bounded", "turnaround", "arclength", "intersection", "grow_until"],
)
def test_every_driver_propagates_the_overrun(drive) -> None:
    """No ``grow_until_*`` driver swallows a ``PointBudgetExceeded``."""
    workbench, fp = _seeded()
    workbench.max_points_per_step = 0
    with pytest.raises(PointBudgetExceeded):
        drive(workbench, fp)


def test_growth_stats_record_the_last_step() -> None:
    """None before growth; afterwards the last step's own point count, a
    positive factor and a finite ``max|x|``; another branch has no stats."""
    workbench, fp = _seeded()
    assert workbench.growth_stats(fp, "unstable") is None
    workbench.grow_n_times(fp, "unstable", STEPS - 1)
    before = _point_count(workbench, fp)
    workbench.grow_n_times(fp, "unstable", 1)
    stats = workbench.growth_stats(fp, "unstable")
    assert stats.points_created == _point_count(workbench, fp) - before
    assert stats.growth_factor > 0
    assert np.isfinite(stats.max_abs_coord) and stats.max_abs_coord > 0
    assert workbench.growth_stats(fp, "unstable", branch_index=1) is None


def test_prediction_hint_warns_before_the_step(caplog) -> None:
    """A step whose predicted size exceeds the budget is announced at WARNING
    on the workbench logger before it runs (and here, overruns)."""
    workbench, fp = _seeded()
    workbench.grow_n_times(fp, "unstable", STEPS - 1)
    assert workbench.growth_stats(fp, "unstable").growth_factor is not None
    workbench.max_points_per_step = 0

    caplog.set_level(logging.WARNING, logger=WORKBENCH_LOGGER)
    with pytest.raises(PointBudgetExceeded):
        workbench.grow_n_times(fp, "unstable", 1)
    assert_logged(caplog, "WARNING", WORKBENCH_LOGGER, count=1)


def test_escape_radius_warns_only_when_set(caplog) -> None:
    """Off by default; set below the reached ``max|x|`` it warns."""
    caplog.set_level(logging.WARNING, logger=WORKBENCH_LOGGER)
    workbench, fp = _seeded()
    workbench.grow_n_times(fp, "unstable", 1)
    assert_logged(caplog, "WARNING", WORKBENCH_LOGGER, count=0)

    workbench.escape_radius = workbench.growth_stats(fp, "unstable").max_abs_coord / 2
    workbench.grow_n_times(fp, "unstable", 1)
    assert_logged(caplog, "WARNING", WORKBENCH_LOGGER, count=1)


def test_grow_until_bounded_stops_on_the_step_that_leaves() -> None:
    """It stops on the first step whose ``max|x|`` exceeds the radius."""
    workbench, fp = _seeded()
    reached = []
    for _ in range(STEPS):
        workbench.grow_n_times(fp, "unstable", 1)
        reached.append(workbench.growth_stats(fp, "unstable").max_abs_coord)
    radius = 0.5 * (reached[-2] + reached[-1])
    assert reached[-1] > radius > max(reached[:-1]), "the probe must leave on its last step"

    workbench, fp = _seeded()
    steps = workbench.grow_until_bounded(fp, "unstable", radius, max_iterations=STEPS + 2)

    assert steps == STEPS
    assert workbench.growth_stats(fp, "unstable").max_abs_coord > radius


def test_grow_until_bounded_raises_at_the_cap() -> None:
    """Like every ``grow_until_*`` sibling, the cap is a ``ValueError``."""
    workbench, fp = _seeded()
    with pytest.raises(ValueError):
        workbench.grow_until_bounded(fp, "unstable", radius=1e300, max_iterations=2)
    assert workbench.growth_stats(fp, "unstable") is not None


# --------------------------------------------------------------------------- #
# blasting
# --------------------------------------------------------------------------- #
def _k28_zone():
    """The k=2.8 session at its resonance zone, unblasted, and the zone."""
    case = build_k28(blasts=0, through="zones")
    return case.session, case.zones[0]


def _spy_bridge_images(workbench) -> tuple[list, list[int]]:
    """Record every bridge the workbench iterates and the map evaluations each
    image cost the machine (the last bridge is the one being iterated when a call
    raised). The iterate inference's map gate creates no point, so it is not counted."""
    evaluated = _count_evaluations(workbench)
    iterate_bridge = workbench.iterate_bridge
    machine = workbench._man_machine
    map_image = machine.iterate_bridge
    bridges: list = []
    costs: list[int] = []

    def spy(bridge):
        bridges.append(bridge)
        return iterate_bridge(bridge)

    def spy_image(bridge):
        evaluated.clear()
        try:
            return map_image(bridge)
        finally:
            costs.append(sum(evaluated))

    workbench.iterate_bridge = spy
    machine.iterate_bridge = spy_image
    return bridges, costs


def _largest_image(session, zone, blasts: int) -> int:
    """Points the largest single bridge image of ``blasts`` blasts creates."""
    _, costs = _spy_bridge_images(session.workbench)
    for _ in range(blasts):
        session.blast_zone(zone, 1)
    return max(costs, default=0)


def test_consecutive_blasts_each_within_budget_never_raise() -> None:
    """Blasts never accumulate a count across bridges or calls."""
    probe, probe_zone = _k28_zone()
    budget = _largest_image(probe, probe_zone, blasts=2)
    assert budget > 0

    session, zone = _k28_zone()
    session.max_points_per_step = budget
    session.blast_zone(zone, 1)
    session.blast_zone(zone, 1)


def test_overrunning_blast_raises_and_never_resumes_half_built() -> None:
    """An image over budget raises out of ``blast_zone`` (never skipped) with
    the generation advanced; the abandoned bridge then refuses a second
    iterate, stays un-iterated and has no registered image."""
    probe, probe_zone = _k28_zone()
    session, zone = _k28_zone()
    workbench = session.workbench
    budget = _largest_image(probe, probe_zone, blasts=1) - 1
    session.max_points_per_step = budget
    bridges, costs = _spy_bridge_images(workbench)
    generation = workbench.generation

    with pytest.raises(PointBudgetExceeded) as raised:
        session.blast_zone(zone, 1)
    assert raised.value.points_created <= budget
    assert costs[-1] <= budget
    assert workbench.generation != generation
    abandoned = bridges[-1]

    session.max_points_per_step = None
    with pytest.raises(RuntimeError) as again:
        session.blast_zone(zone, 1)
    assert not isinstance(again.value, PointBudgetExceeded)

    assert bridges[-1] is abandoned
    assert not abandoned.iterated
    if abandoned.id is not None:
        assert not workbench.image_bridges(abandoned.id)


# --------------------------------------------------------------------------- #
# area_cutoff: public and shared
# --------------------------------------------------------------------------- #
def test_area_cutoff_reaches_the_one_shared_machine() -> None:
    """Constructor and setter, workbench and session: one machine refines
    everything, the initializer's segments included."""
    session = TangleSession(
        henon_map(*HENON_K28), henon_map_inverse(*HENON_K28), area_cutoff=3e-6
    )
    workbench = session.workbench
    assert workbench.area_cutoff == session.area_cutoff == 3e-6
    assert workbench._man_maker.machine is workbench._man_machine

    session.area_cutoff = 2e-7
    assert workbench._man_machine.area_cutoff == 2e-7
    assert workbench._man_maker.machine.area_cutoff == 2e-7

    session.max_points_per_step = 123
    session.escape_radius = 5.0
    assert workbench._man_machine.max_points_per_step == session.max_points_per_step == 123
    assert workbench.escape_radius == session.escape_radius == 5.0
    assert not {"area_cutoff", "max_points_per_step", "escape_radius"} & set(vars(session))


def test_finer_cutoff_never_gives_fewer_points() -> None:
    """Refinement is monotone in the cutoff."""
    counts = []
    for cutoff in (1e-4, 1e-6):
        session = TangleSession(
            henon_map(*HENON_K28),
            henon_map_inverse(*HENON_K28),
            henon_jacobian(*HENON_K28),
            area_cutoff=cutoff,
        )
        fp = session.construct_fixed_point(list(K28_SADDLE_SEED))
        session.orient_eigenvectors(fp, _orientation(PERIOD1_ORIENTATION))
        session.initialize_both_manifolds(fp)
        session.grow_n_times(fp, "unstable", 6)
        counts.append(_point_count(session.workbench, fp))
    assert counts[1] >= counts[0]


# --------------------------------------------------------------------------- #
# perf: the killer step stops in bounded time
# --------------------------------------------------------------------------- #
@pytest.mark.perf
def test_escaping_growth_stops_in_bounded_time() -> None:
    """k=2.8 at 1e-5 grown 12+ steps with a small budget raises quickly
    (the stress campaign's E3 stalled ~190 s inside one step)."""
    session = TangleSession(
        henon_map(*HENON_K28),
        henon_map_inverse(*HENON_K28),
        henon_jacobian(*HENON_K28),
        area_cutoff=1e-5,
        max_points_per_step=200_000,
    )
    fp = session.construct_fixed_point(list(K28_SADDLE_SEED))
    session.orient_eigenvectors(fp, _orientation(PERIOD1_ORIENTATION))
    session.initialize_both_manifolds(fp)

    start = time.perf_counter()
    with pytest.raises(PointBudgetExceeded):
        session.grow_n_times(fp, "unstable", 16)
    assert time.perf_counter() - start < 60.0
