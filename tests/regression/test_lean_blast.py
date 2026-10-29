"""A lean blast computes what a full blast computes, step by step.

A lean blast frees, after each iteration, every bridge interior neither the
topology nor a later blast reads (``TangleWorkbench.release_curves``). It must
not change anything computed: each series is blasted one iteration at a time,
lean and full with the same guard (none, or the width guard), and after every
blast the registry, the iterate table, the pip depth and the topology
(identities, words, the refined matrix, the inert classes) must agree. The
crossing ids are compared as they are: breadth-first order is unchanged, so a
lean build numbers its crossings as the full one does.

The series cover the keep rule's cases: ``nested_outer_series`` keeps another
zone's frontier, ``inversion_series`` the pieces of joined arcs. No count or
word of a case is pinned; the two builds are compared with each other.
"""

from __future__ import annotations

import pytest

from cases import BUILDERS, SERIES
from helpers.direct import view_of
from helpers.names import homotopy_pair


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


def _blasts(build, blasts: int, separation, lean: bool) -> tuple[list[dict], int, int]:
    """Blast a build's largest zone one iteration at a time.

    Returns:
        The steps, the rows the lean blasts freed, the store's final rows.
    """
    case = build()
    session = case.session
    zone = max(case.zones, key=lambda zone: zone.area)
    steps, released = [], 0
    for _ in range(blasts):
        result = session.blast_zone(
            zone, 1, fixed_point=[zone.fixed_point], min_separation=separation,
            diagnose=False if not lean else None, lean=lean,
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
    full, _, _ = _blasts(series.build, series.blasts, separation, lean=False)
    lean, _, _ = _blasts(series.build, series.blasts, separation, lean=True)
    _assert_equal(full, lean)


@pytest.mark.parametrize("name", ["nested_inner_blasts", "heteroclinic"])
def test_a_lean_outer_blast_equals_the_full_one(name):
    """The outer zone blasted lean after the inner one was blasted (its frontier
    is another zone's, kept), and the heteroclinic tangle, whose per-fixed-point
    trellis joins arcs across a crossing of the other tangle."""
    full, _, _ = _blasts(BUILDERS[name], 3, "width", lean=False)
    lean, _, _ = _blasts(BUILDERS[name], 3, "width", lean=True)
    _assert_equal(full, lean)


def test_a_lean_blast_frees_rows():
    """The width-guarded k28 series frees rows: the store ends smaller than the full one's."""
    series = SERIES["k28_series"]
    _, _, full_rows = _blasts(series.build, series.blasts, "width", lean=False)
    _, released, lean_rows = _blasts(series.build, series.blasts, "width", lean=True)
    assert released > 0
    assert lean_rows < full_rows
