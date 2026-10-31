"""A representative blast iterates one bridge per refined class and holds the rest.

``blast_zone(representatives=...)`` reads the topology at the start of each
step: a bridge ending on the pip orbit, one the minimal trellis maps and one
the class table does not know are iterated; every other is bucketed by its
refined child (else its unmatched pair, else its class; a loop by class and
element) and one per bucket is iterated. The tests restate that rule on the
topology read before each call and compare the blast with itself on two
builds; no count or word of the case is pinned.
"""

from __future__ import annotations

from dataclasses import replace

import pytest

from cases import SERIES
from tanglepack.loom.Blast import _pip_orbits, _representatives

POLICIES = ["fewest_points", "smallest_id", "largest_id"]

#: Enough single blasts of the k28 series for the buckets to hold something.
BLASTS = 5


def _zone():
    """The k28 series' zone (k=2.8 at 1e-7, trimmed at ``f(q0)``), unblasted."""
    case = SERIES["k28_series"].build()
    zone = max(case.zones, key=lambda zone: zone.area)
    return case.session, zone.fixed_point, zone


def _blast(session, fp, zone, n, policy, **options):
    return session.blast_zone(
        zone, n, fixed_point=[fp], min_separation=SERIES["k28_series"].min_separation,
        diagnose=False, representatives=policy, **options,
    )


def _bucket(topology, bridge_id):
    """The bucket of a bridge the rule may hold, None for one it always iterates."""
    minimal, dynamics = topology.minimal, topology.dynamics
    if bridge_id in set(minimal.kept_bridge_ids) - set(minimal.inert_hole_bridge_ids):
        return None
    try:
        entry = minimal.classes.entry_of(bridge_id)
    except KeyError:
        return None
    member = next(m for m in entry.members if m.bridge_id == bridge_id)
    if member.is_loop:
        return entry.bridge_class, "loop", member.loop_element
    child = dynamics.member_refinement.get(bridge_id)
    if child is not None:
        return entry.bridge_class, child.name
    pair = dynamics.unmatched_members.get(bridge_id)
    return (entry.bridge_class, "unmatched", pair) if pair else (entry.bridge_class,)


def _steps(policy):
    """``BLASTS`` single representative blasts, each with the topology and point
    counts read before it; a held bridge is still un-iterated after its call."""
    session, fp, zone = _zone()
    steps = []
    for _ in range(BLASTS):
        topology = session.topology([fp])
        pip_ends = {x for orbit in _pip_orbits(session, [fp]).values() for x in orbit}
        points = {b.id: b.point_count for b in session.workbench.bridges}
        result = _blast(session, fp, zone, 1, policy)
        step = result.steps[0]
        assert not any(session.workbench.bridge(bid).iterated for bid in step.held)
        steps.append((topology, pip_ends, points, step))
    held = [step for *_, step in steps if step.held]
    assert held, "the build holds nothing: the rule is not exercised"
    return session, steps


@pytest.mark.parametrize("policy", POLICIES)
def test_every_iterated_bridge_is_one_the_rule_iterates(policy):
    """Each parent ends on the pip orbit, is mapped by the minimal trellis or
    unclassed, or is the only parent of its bucket."""
    _, steps = _steps(policy)
    for topology, pip_ends, _, step in steps:
        buckets = [
            _bucket(topology, b.id) for b in step.interior_parents
            if b.id is not None and not pip_ends.intersection(b.id)
        ]
        bucketed = [bucket for bucket in buckets if bucket is not None]
        assert len(bucketed) == len(set(bucketed)), "two representatives of one bucket"


@pytest.mark.parametrize("policy", POLICIES)
def test_a_held_bridge_shares_the_bucket_of_its_iterated_representative(policy):
    """Every held bridge's representative was iterated in its step, from the same bucket."""
    _, steps = _steps(policy)
    for topology, pip_ends, _, step in steps:
        parents = {b.id for b in step.interior_parents}
        for held, representative in step.held.items():
            assert representative in parents
            assert held not in parents
            assert not pip_ends.intersection(held)
            assert _bucket(topology, held) is not None
            assert _bucket(topology, held) == _bucket(topology, representative)


@pytest.mark.parametrize("policy", POLICIES)
def test_each_policy_picks_the_member_its_rule_names(policy):
    """``fewest_points`` (ties by id), ``smallest_id`` and ``largest_id`` against every held bridge."""
    _, steps = _steps(policy)
    for *_, points, step in steps:
        for held, representative in step.held.items():
            if policy == "fewest_points":
                assert (points[representative], representative) < (points[held], held)
            elif policy == "smallest_id":
                assert representative < held
            else:
                assert representative > held


def test_each_step_iterates_or_holds_every_candidate():
    """A step's candidates (its frontier and the bridges the step before held)
    are exactly its iterated parents and its held bridges."""
    session, fp, zone = _zone()
    result = _blast(session, fp, zone, BLASTS, "fewest_points")
    assert result.held > 0
    waiting = []
    for frontier, step in zip(result.interior_bridges_by_iteration, result.steps):
        candidates = [*frontier, *waiting]
        waiting = [b for b in candidates if b.id in step.held]
        assert len(waiting) == len(step.held)
        assert sorted(map(id, [*step.interior_parents, *waiting])) == sorted(map(id, candidates))


def test_an_unmatched_member_is_bucketed_by_its_pair():
    """No build has an unmatched member, so a held bridge and its representative
    are moved to the cached dynamics' ``unmatched_members``: with two pairs both
    are iterated, with one pair one of them is held."""
    session, fp, zone = _zone()
    for _ in range(3):
        _blast(session, fp, zone, 1, None)
    topology = session.topology([fp])
    candidates = [
        session.workbench.bridge(member.bridge_id)
        for entry in topology.minimal.classes for member in entry.members
    ]
    _, held = _representatives(session, [fp], candidates, "fewest_points")
    assert held, "the build holds nothing: the rule is not exercised"
    bridge_id, representative = next(iter(held.items()))
    bridge_class = topology.minimal.classes.entry_of(bridge_id).bridge_class
    dynamics = topology.dynamics
    for pairs in [
        ((bridge_class.source, bridge_class.target), (bridge_class.target, bridge_class.source)),
        ((bridge_class.source, bridge_class.target),) * 2,
    ]:
        for bid, pair in zip((bridge_id, representative), pairs):
            dynamics.member_refinement[bid] = None
            dynamics.unmatched_members[bid] = pair
        iterated, held = _representatives(session, [fp], candidates, "fewest_points")
        assert session.topology([fp]) is topology
        ids = {b.id for b in iterated}
        if pairs[0] == pairs[1]:
            assert len({bridge_id, representative} & ids) == 1
            assert {bridge_id, representative} & set(held)
        else:
            assert {bridge_id, representative} <= ids
            assert not {bridge_id, representative} & set(held)


def test_a_loop_member_is_bucketed_by_its_element():
    """No build holds a loop, so a held bridge and its representative are made
    loop members of their class entry: on two elements both are iterated, on
    one element one of them is held."""
    session, fp, zone = _zone()
    for _ in range(3):
        _blast(session, fp, zone, 1, None)
    topology = session.topology([fp])
    candidates = [
        session.workbench.bridge(member.bridge_id)
        for entry in topology.minimal.classes for member in entry.members
    ]
    _, held = _representatives(session, [fp], candidates, "fewest_points")
    assert held, "the build holds nothing: the rule is not exercised"
    bridge_id, representative = next(iter(held.items()))
    entry = topology.minimal.classes.entry_of(bridge_id)
    bridge_class = entry.bridge_class
    for elements in [(bridge_class.source, bridge_class.target), (bridge_class.source,) * 2]:
        loops = dict(zip((bridge_id, representative), elements))
        entry.members[:] = [
            replace(m, direction=0, loop_element=loops[m.bridge_id]) if m.bridge_id in loops else m
            for m in entry.members
        ]
        iterated, held = _representatives(session, [fp], candidates, "fewest_points")
        assert session.topology([fp]) is topology
        ids = {b.id for b in iterated}
        if elements[0] == elements[1]:
            assert len({bridge_id, representative} & ids) == 1
            assert {bridge_id, representative} & set(held)
        else:
            assert {bridge_id, representative} <= ids
            assert not {bridge_id, representative} & set(held)


def _state(session):
    """The registry, its links and the bridges with their ``iterated`` flags."""
    registry = session.workbench.intersection_registry
    return (
        {iid: (ix.unstable_cdist, ix.stable_cdist, ix.crossing_sign) for iid, ix in registry},
        {iid: registry.iterate_table[iid, 1] for iid, _ in registry},
        {b.id: b.iterated for b in session.workbench.bridges if b.id is not None},
    )


def test_single_representative_blasts_compose():
    """``n`` one-step calls leave the workbench one ``n``-step call leaves."""
    stepwise = _zone()
    for _ in range(BLASTS):
        _blast(*stepwise, 1, "fewest_points")
    at_once = _zone()
    result = _blast(*at_once, BLASTS, "fewest_points")
    assert result.held > 0
    assert _state(stepwise[0]) == _state(at_once[0])


def test_a_full_blast_iterates_every_held_bridge():
    """Nothing is remembered on a held bridge: the next blast without representatives takes it."""
    session, fp, zone = _zone()
    result = _blast(session, fp, zone, BLASTS, "fewest_points")
    held = set(result.steps[-1].held)
    assert len(held) == result.held > 0
    full = _blast(session, fp, zone, 1, None)
    assert held <= {b.id for b in full.steps[0].interior_parents}
    assert all(session.workbench.bridge(bid).iterated for bid in held)


def test_without_representatives_nothing_is_held():
    """``None`` iterates the whole frontier and records no hold."""
    session, fp, zone = _zone()
    result = _blast(session, fp, zone, BLASTS, None)
    assert result.held == 0
    assert all(step.held == {} for step in result.steps)
    for frontier, step in zip(result.interior_bridges_by_iteration, result.steps):
        assert step.interior_parents == frontier


def test_a_lean_representative_blast_frees_rows():
    """Lean and representatives compose: the width-guarded series frees rows."""
    session, fp, zone = _zone()
    series = SERIES["k28_series"]
    result = session.blast_zone(
        zone, series.blasts, fixed_point=[fp], min_separation="width", lean=True,
        representatives="fewest_points",
    )
    assert result.held > 0
    assert result.released_rows > 0
