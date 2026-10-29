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

from cases import build_k28
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
    """Every width drop is below the bound; every other child is not, the pip-orbit ones aside."""
    session, fp, zone = _k28_zone()
    result = _blast(session, fp, zone, 8)
    orbit = {x for chain in _pip_orbits(session, [fp]).values() for x in chain}
    dropped = [c for step in result.steps for c in step.discarded_too_close]
    kept = [c for step in result.steps for c in step.kept_interior]
    assert dropped and kept
    for child in dropped:
        assert child.dropped_at == WIDTH_DROP
        assert child.width < WIDTH_SEPARATION_FACTOR * _deviation(session, child)
    for child in kept:
        if not set(child.id or ()) & orbit and child.width is not None:
            assert child.width >= WIDTH_SEPARATION_FACTOR * _deviation(session, child)


def test_a_width_contracts_by_the_stable_factor():
    """A child born of a seeded parent carries the parent's width times the stable factor."""
    session, fp, zone = _k28_zone()
    frontier = [b for b in session.workbench.uniiterated_bridges if b.width is None]
    result = _blast(session, fp, zone, 1)
    contraction = fp.per_step_beta("stable")
    checked = 0
    for parent in result.steps[0].interior_parents:
        if parent.width is None or parent.id is None or parent not in frontier:
            continue
        assert parent.width > 0.0
        for child_id in session.workbench.image_bridges(parent.id) or ():
            child = session.workbench.bridge(child_id)
            if child in frontier:
                continue
            assert child.width == pytest.approx(parent.width * contraction, rel=1e-12)
            checked += 1
    assert checked


def test_the_pip_orbit_child_is_kept(monkeypatch):
    """With a factor no lobe passes, the only children kept end on the pip's orbit."""
    monkeypatch.setattr(Blast, "WIDTH_SEPARATION_FACTOR", 1e30)
    session, fp, zone = _k28_zone()
    result = _blast(session, fp, zone, 3)
    orbit = {x for chain in _pip_orbits(session, [fp]).values() for x in chain}
    kept = [c for step in result.steps for c in step.kept_interior]
    assert kept
    assert all(set(c.id or ()) & orbit or c.width is None for c in kept)
    assert any(set(c.id or ()) & orbit for c in kept)


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
