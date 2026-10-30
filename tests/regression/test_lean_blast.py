"""A lean blast computes what a full blast computes, step by step.

A lean blast frees, after each iteration, every bridge interior neither the
topology nor a later blast reads (``TangleWorkbench.release_curves``). It must
not change anything computed: each series is blasted one iteration at a time,
lean and full with the same guard (none, or the width guard), and after every
blast the registry, the iterate table, the pip depth and the topology
(identities, words, the refined matrix, the inert classes) must agree. The
crossing ids are compared as they are: breadth-first order is unchanged, so a
lean build numbers its crossings as the full one does.

A series whose lean blasts free nothing compares nothing, and is a visible
skip (``inversion_series`` frees no row). The outer zone of a nested and of a
heteroclinic build is blasted lean and then the inner one, which iterates the
frontier the outer blasts had to keep; the keep rule's other clauses (the pip
orbit, the minimal trellis of every fixed point) are asserted directly after
each lean blast, and no series holds a joined arc over released curve, so the
joined-pieces rule is tested on its own
(``test_lean_blast_refusals.py::test_a_lean_blast_keeps_a_trellis_joined_arc``).
No count or word of a case is pinned; the two builds are compared with each
other.
"""

from __future__ import annotations

import pytest

from cases import BUILDERS, SERIES
from helpers.direct import view_of
from helpers.names import homotopy_pair
from tanglepack.loom.Blast import _lean_keep, _pieces, _pip_orbits
from tanglepack.numerics.Bridge import Bridge


def _canonical(value, fixed_points):
    """``value`` with fixed points as their indices and sets sorted: comparable across builds."""
    if any(value is fp for fp in fixed_points):
        return ("fp", next(i for i, fp in enumerate(fixed_points) if value is fp))
    if isinstance(value, (set, frozenset)):
        return ("set", tuple(sorted((_canonical(v, fixed_points) for v in value), key=repr)))
    if isinstance(value, dict):
        return tuple(sorted(
            ((_canonical(k, fixed_points), _canonical(v, fixed_points)) for k, v in value.items()),
            key=repr,
        ))
    if isinstance(value, (tuple, list)):
        return tuple(_canonical(v, fixed_points) for v in value)
    return value


def _step(session, result) -> dict:
    """What one blast left: registry, table, depth and topology, canonical."""
    workbench = session.workbench
    registry = workbench.intersection_registry
    topology = session.topology()
    dynamics = topology.dynamics
    view = {
        "registry": {
            iid: (ix.unstable_cdist, ix.stable_cdist, ix.manifold_a_key, ix.manifold_b_key, ix.crossing_sign)
            for iid, ix in registry
        },
        "table": {iid: registry.iterate_table[iid, 1] for iid, _ in registry},
        "pip_depth": result.pip_depth,
        "identities": {topology.identity(entry.bridge_class) for entry in topology.classes},
        "dynamics": view_of("symbolic_dynamics", dynamics),
        "inert": {homotopy_pair(dynamics, dynamics.classes[cls]) for cls in dynamics.inert_letters},
    }
    return _canonical(view, workbench.fixed_points)


def _blasts(
    build, blasts: int, separation, lean: bool, inner_blasts: int = 0
) -> tuple[list[dict], int, int]:
    """Blast a build's largest zone one iteration at a time, then its smallest.

    Returns:
        The steps, the rows the lean blasts freed, the store's final rows.
    """
    case = build()
    session = case.session
    zones = sorted(case.zones, key=lambda zone: zone.area)
    steps, released = [], 0
    for zone in [zones[-1]] * blasts + [zones[0]] * inner_blasts:
        result = session.blast_zone(
            zone, 1, fixed_point=[zone.fixed_point], min_separation=separation,
            diagnose=False, lean=lean,
        )
        released += result.released_rows
        steps.append(_step(session, result))
    return steps, released, len(session.workbench.point_store)


def _assert_equal(full: list[dict], lean: list[dict]) -> None:
    """Every step of the lean run equals the full run's, part by part."""
    assert len(full) == len(lean)
    for blast, (a, b) in enumerate(zip(full, lean), start=1):
        differing = [key for (key, value), (_, other) in zip(a, b) if value != other]
        assert not differing, f"blast {blast}: {differing} differ"


@pytest.mark.parametrize("separation", [None, "width"], ids=["unguarded", "width"])
@pytest.mark.parametrize("name", sorted(SERIES))
def test_a_lean_series_equals_the_full_one(name, separation):
    """Every step of a lean series equals the full series' step."""
    series = SERIES[name]
    lean, released, _ = _blasts(series.build, series.blasts, separation, lean=True)
    if not released:
        pytest.skip(f"a lean {name} frees no row: there is nothing to compare")
    full, _, _ = _blasts(series.build, series.blasts, separation, lean=False)
    _assert_equal(full, lean)


@pytest.mark.parametrize("name", ["nested_inner_blasts", "heteroclinic"])
def test_a_lean_outer_blast_equals_the_full_one(name):
    """The outer zone blasted lean, then the inner one: the inner frontier is
    another zone's while the outer zone is blasted, kept and iterated after.
    The heteroclinic tangle's per-fixed-point trellis joins arcs across a
    crossing of the other tangle."""
    full, _, _ = _blasts(BUILDERS[name], 3, "width", lean=False, inner_blasts=2)
    lean, released, _ = _blasts(BUILDERS[name], 3, "width", lean=True, inner_blasts=2)
    assert released > 0
    _assert_equal(full, lean)


def test_a_lean_blast_keeps_what_later_calls_read():
    """After each lean blast of the heteroclinic build's inner zone (the first
    of the suite's builds to release bridges there) nothing a later call reads
    is released: the outer zone's frontier, every bridge ending on the pip's
    orbit (``cut_at_pip_orbit`` cuts there), and the minimal trellis of the
    selection and of every fixed point (``topology()`` over all of them)."""
    case = BUILDERS["heteroclinic"]()
    session, workbench = case.session, case.session.workbench
    zone = min(case.zones, key=lambda zone: zone.area)
    for _ in range(5):
        session.blast_zone(zone, 1, fixed_point=[zone.fixed_point], min_separation="width", lean=True)
        frontiers = [
            b for b in workbench.uniiterated_bridges
            if b.dropped_at is None and session.classify_bridge(b) is not None
        ]
        orbit = [
            workbench.bridge(bid)
            for chain in _pip_orbits(session, [zone.fixed_point]).values()
            for x in chain for bid in workbench.bridges_at(x)
        ]
        minimal = [
            workbench.bridge(bid)
            for selection in ([zone.fixed_point], None)
            for bid in session.topology(selection).minimal.kept_bridge_ids
        ]
        for kind, bridges in (("frontier", frontiers), ("pip orbit", orbit), ("minimal", minimal)):
            assert bridges, kind
            released = [b.id for b in bridges if b.released]
            assert not released, f"{kind}: {released} released"
    assert any(b.released for b in workbench.bridges)


def test_a_lean_series_keeps_its_pip_orbit():
    """Through the width-guarded k28 series every bridge ending on the pip's
    orbit stays: k28 is the series whose lean blasts would release some."""
    series = SERIES["k28_series"]
    case = series.build()
    session, workbench = case.session, case.session.workbench
    zone = max(case.zones, key=lambda zone: zone.area)
    for _ in range(series.blasts):
        session.blast_zone(zone, 1, fixed_point=[zone.fixed_point], min_separation="width", lean=True)
        orbit = {
            bid for chain in _pip_orbits(session, [zone.fixed_point]).values()
            for x in chain for bid in workbench.bridges_at(x)
        }
        assert orbit
        assert not [bid for bid in orbit if workbench.bridge(bid).released]
    assert any(b.released for b in workbench.bridges)


def test_the_keep_rule_names_every_reader():
    """``_lean_keep`` of one fixed point's zone names the other zone's frontier,
    the minimal trellis of every fixed point, every bridge ending on the pip's
    orbit and the pieces of every joined arc a cached trellis holds."""
    case = BUILDERS["heteroclinic"]()
    session, workbench = case.session, case.session.workbench
    zone = min(case.zones, key=lambda zone: zone.area)
    session.blast_zone(zone, 2, fixed_point=[zone.fixed_point], min_separation="width", lean=True)
    keep = {
        k.id if isinstance(k, Bridge) else k
        for k in _lean_keep(session, [zone.fixed_point], lambda bridge: False)
    }
    readers = {
        "frontier": {
            b.id for b in workbench.uniiterated_bridges
            if b.dropped_at is None and session.classify_bridge(b) is not None
        },
        "minimal": set(session.topology(None).minimal.kept_bridge_ids),
        "pip orbit": {
            bid for chain in _pip_orbits(session, [zone.fixed_point]).values()
            for x in chain for bid in workbench.bridges_at(x)
        },
        "joined": {
            bid for trellis in session.cached_trellises
            for joined in trellis.joined_arcs for bid in _pieces(workbench, joined)
        },
    }
    for kind, bridge_ids in readers.items():
        assert bridge_ids, kind
        assert bridge_ids <= keep, f"{kind}: {sorted(bridge_ids - keep)} not kept"


def test_a_lean_blast_frees_rows():
    """The width-guarded k28 series frees rows: the store ends smaller than the full one's."""
    series = SERIES["k28_series"]
    _, _, full_rows = _blasts(series.build, series.blasts, "width", lean=False)
    _, released, lean_rows = _blasts(series.build, series.blasts, "width", lean=True)
    assert released > 0
    assert lean_rows < full_rows
