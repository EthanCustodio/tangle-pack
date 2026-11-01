"""The width guard, ``min_separation="width"``: a lobe thinner than its polygon error is dropped.

A bridge's width is seeded from the ring it and its stable segment bound
(``|area| / L``) when it is first iterated, and contracts by the stable factor
down the image chain; a child is dropped when its width is below
``WIDTH_SEPARATION_FACTOR`` chord deviations ``1.5 * area_cutoff / L`` of its
own polyline (``L`` its median chord). Every test compares the guard with its
rule or with itself on two builds, never with a case number.
"""

from __future__ import annotations

import numpy as np
import pytest

from cases import BUILDERS, build_k28
from tanglepack.loom import Blast
from tanglepack.loom.Blast import WIDTH_DROP, WIDTH_SEPARATION_FACTOR, _pip_orbits


def _k28_zone():
    """The k=2.8 zone trimmed at ``f(q0)``, unblasted."""
    case = build_k28(blasts=0, through="zones")
    return case.session, case.fixed_points[0], case.zones[0]


def _blast(session, fp, zone, n):
    """``n`` width-guarded iterations of ``zone`` in one call."""
    return session.blast_zone(zone, n, fixed_point=[fp], min_separation="width", diagnose=False)


def _deviation(session, bridge) -> float:
    """The chord deviation of a bridge's polyline at the machine's cutoff."""
    chords = np.hypot(*np.diff(bridge.get_point_array(), axis=0).T)
    return 1.5 * session.workbench.area_cutoff / float(np.median(chords))


def test_a_child_is_dropped_iff_thinner_than_its_deviation():
    """Every width drop is below the bound; every other child is not."""
    session, fp, zone = _k28_zone()
    result = _blast(session, fp, zone, 8)
    dropped = [c for step in result.steps for c in step.discarded_too_close]
    kept = [c for step in result.steps for c in step.kept_interior]
    assert dropped and kept
    for child in dropped:
        assert child.dropped_at == WIDTH_DROP
        assert child.width < WIDTH_SEPARATION_FACTOR * _deviation(session, child)
    for child in kept:
        if child.width is not None:
            assert child.width >= WIDTH_SEPARATION_FACTOR * _deviation(session, child)


def _contractions(session, result) -> list[tuple["Bridge", "Bridge"]]:
    """``(parent, child)`` for every child a seeded parent of the first step
    gave its width, checking the child's width and saddle on the way."""
    pairs = []
    for parent in result.steps[0].interior_parents:
        if parent.width is None or parent.id is None:
            continue
        assert parent.width > 0.0
        contraction = parent.width_saddle.per_step_beta("stable")
        for child_id in session.workbench.image_bridges(parent.id) or ():
            child = session.workbench.bridge(child_id)
            if child.width_saddle is parent.width_saddle and child.width == pytest.approx(
                parent.width * contraction, rel=1e-12
            ):
                pairs.append((parent, child))
    return pairs


def test_a_width_contracts_by_the_stable_factor():
    """A child born of a seeded parent carries the parent's width times the stable factor."""
    session, fp, zone = _k28_zone()
    result = _blast(session, fp, zone, 1)
    pairs = _contractions(session, result)
    assert pairs
    assert all(parent.width_saddle is fp for parent, _ in pairs)


def test_a_heteroclinic_lobe_contracts_by_its_stable_saddle():
    """A lobe of ``W^u(B)`` based on ``W^s(A)`` contracts by A's stable factor
    (NUM-20): its saddle is the fixed point of the stable branch its ring lies on."""
    case = BUILDERS["heteroclinic_17u"]()
    session = case.session
    zone = max(case.zones, key=lambda zone: zone.area)
    result = session.blast_zone(zone, 1, fixed_point=None, min_separation="width", diagnose=False)
    registry = session.workbench.intersection_registry
    heteroclinic = [
        (parent, child) for parent, child in _contractions(session, result)
        if parent.width_saddle is not parent.fixed_point
    ]
    assert heteroclinic
    assert all(
        parent.width_saddle.per_step_beta("stable") != parent.fixed_point.per_step_beta("stable")
        for parent, _ in heteroclinic
    )
    for parent, _ in heteroclinic:
        base = {registry[i].manifold_b_key[0] for i in parent.id if not registry[i].is_anchor}
        assert base == {parent.width_saddle}


def test_the_pip_orbit_child_is_not_exempt(monkeypatch):
    """With a factor no lobe passes, no child with a width is kept, one ending
    on the pip's orbit included (the exemption went on 2026-10-31)."""
    monkeypatch.setattr(Blast, "WIDTH_SEPARATION_FACTOR", 1e30)
    session, fp, zone = _k28_zone()
    result = _blast(session, fp, zone, 3)
    orbit = {x for chain in _pip_orbits(session, [fp]).values() for x in chain}
    assert all(c.width is None for step in result.steps for c in step.kept_interior)
    assert any(set(c.id or ()) & orbit for step in result.steps for c in step.discarded_too_close)


def _state(session, fp):
    """``iterated``, ``dropped_at`` and ``width`` per bridge id, and the pip depth."""
    return (
        {b.id: (b.iterated, b.dropped_at, b.width) for b in session.workbench.bridges if b.id is not None},
        [len(orbit) - 1 for orbit in _pip_orbits(session, [fp]).values()],
    )


def test_single_width_blasts_compose():
    """``n`` single-iteration width blasts leave the workbench one ``n``-iteration blast leaves."""
    stepwise = _k28_zone()
    for _ in range(8):
        _blast(*stepwise, 1)
    at_once = _k28_zone()
    result = _blast(*at_once, 8)
    assert result.too_close > 0
    assert _state(*stepwise[:2]) == _state(*at_once[:2])


def test_a_width_drop_stands_until_none():
    """A later float or width blast skips a width drop; only None re-admits it."""
    session, fp, zone = _k28_zone()
    _blast(session, fp, zone, 8)
    dropped = [b for b in session.workbench.bridges if b.dropped_at == WIDTH_DROP]
    assert dropped
    for separation in ("width", 1e-5):
        result = session.blast_zone(zone, 1, fixed_point=[fp], min_separation=separation, diagnose=False)
        assert result.skipped_dropped >= len(dropped)
    result = session.blast_zone(zone, 1, fixed_point=[fp], min_separation=None, diagnose=False)
    frontier = {id(b) for b in result.interior_bridges_by_iteration[0]}
    assert all(id(b) in frontier for b in dropped)


@pytest.mark.parametrize(
    ("separation", "diagnose"),
    [(1e-5, None), ("resolution", None), (None, True)],
    ids=["float", "resolution", "diagnose"],
)
def test_a_lean_blast_refuses_what_reads_every_interior(separation, diagnose):
    """A geometric separation or a diagnose pass needs every interior: ValueError."""
    session, fp, zone = _k28_zone()
    with pytest.raises(ValueError):
        session.blast_zone(zone, 1, min_separation=separation, diagnose=diagnose, lean=True)
