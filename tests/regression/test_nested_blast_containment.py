"""
Regression (WP9, NEST-4, C4): a blast iterates the bridges whose INNERMOST zone is
the blasted one.

With nested zones an outer blast used to take every bridge whose midpoint lay in
the outer polygon, the inner tangle's bridges inside the inner zone included, so
an unfiltered (``fixed_point=None``) blast depended on the order of the outer and
inner blasts (period 3 kept 4 holes outer-first, 42 inner-first). The nested
zones come from the tests-only :mod:`cases` builder.
"""

from __future__ import annotations

from collections import Counter

from cases import build_nested


def _zones(case):
    """``(inner, outer)`` zones of a nested build."""
    zones = sorted(case.zones, key=lambda zone: zone.area)
    return zones[0], zones[-1]


def _blast(case, zone) -> None:
    case.session.blast_zone(zone, num_iterations=1, fixed_point=None, min_separation=1e-4)


def test_an_unfiltered_outer_blast_iterates_only_its_innermost_bridges():
    case = build_nested(outer_blasts=0, through="zones")
    session = case.session
    inner, outer = _zones(case)
    inner_bridges = [
        bridge for bridge in session.workbench.uniiterated_bridges
        if session.classify_bridge(bridge) is inner
    ]
    assert inner_bridges
    result = session.blast_zone(outer, num_iterations=1, fixed_point=None, min_separation=1e-4)
    parents = result.steps[0].interior_parents
    assert parents and all(session.classify_bridge(bridge) is outer for bridge in parents)
    assert not [bridge.id for bridge in inner_bridges if bridge.iterated]


def _signature(case) -> tuple:
    """Holes per fixed point (by period), class count, and every class of each tangle
    spelled by its homotopy element names with its word."""
    session = case.session
    topology = session.topology()
    trellis = topology.trellis
    holes = Counter(
        trellis.intersection(hole.origin[0]).manifold_b_key[0].period for hole in topology.holes
    )
    dynamics = session.symbolic_dynamics()

    def spelled(bridge_class):
        names = (dynamics.naming.homotopy_name(ref) for ref in (bridge_class.source,
                                                                bridge_class.target))
        return tuple((name.orbit_code, name.short_text) for name in names)

    words = sorted(
        (cd.bridge_class.source.fixed_point.period, spelled(cd.bridge_class),
         tuple((spelled(s.bridge_class), s.direction) for s in cd.symbols))
        for cd in dynamics.classes.values()
    )
    return dict(holes), len(dynamics.classes), words


def test_an_unfiltered_nested_blast_does_not_depend_on_order():
    outer_first = build_nested(outer_blasts=0, through="zones")
    inner, outer = _zones(outer_first)
    for zone in [outer] * 4 + [inner] * 2:
        _blast(outer_first, zone)

    inner_first = build_nested(outer_blasts=0, through="zones")
    inner, outer = _zones(inner_first)
    for zone in [inner] * 2 + [outer] * 4:
        _blast(inner_first, zone)

    assert _signature(outer_first) == _signature(inner_first)
