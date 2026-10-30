"""The blast's zone filter reads a pip-orbit bridge's pip end (2026-10-29).

An un-iterated bridge ARRIVING at ``f^j(q0)``, ``j >= 1`` (its second crossing
on the registered orbit of the blasted zone's own corner pip
``ResonanceZone.boundary_intersection``, whichever pip the session has chosen), is the zone's
whatever its midpoint says: near the periodic point the frozen zone polygon
cannot place a thin lobe's midpoint, and that filter, not the dynamics, set the
depth of ``T(q0)``. Another fixed point's pip-orbit bridges keep the
innermost-zone rule, so nested blasting stays order independent. Bridges are
found by rule (their ends on the composed orbit), never by id; a midpoint is
moved out of the zone by patching ``Bridge.test_point``.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from tanglepack.numerics.Bridge import Bridge

MIN_SEPARATION = 0.02
FAR = np.array([1e3, 1e3])


def _orbit(session, fp, pip) -> list[int]:
    """``[q0, f(q0), ..., f^m(q0)]`` of the pip ``q0``, composed through ``trellis.iterate``."""
    trellis = session.trellis(fp)
    orbit = [trellis.registry.find(pip)]
    while (image := trellis.iterate(orbit[-1], 1)) is not None:
        orbit.append(image)
    return orbit


def _arriving(session, zone, at=slice(1, None)) -> list[Bridge]:
    """The un-iterated bridges of the zone's fixed point whose second crossing is ``orbit[at]`` of its pip."""
    orbit = _orbit(session, zone.key, zone.boundary_intersection)[at]
    return [
        b for b in session.workbench.uniiterated_bridges
        if b.fixed_point is zone.key and b.id is not None and b.id[1] in orbit
    ]


def _iterated_first(result) -> set[int]:
    return {id(b) for b in result.interior_bridges_by_iteration[0]}


@pytest.fixture
def outside(monkeypatch):
    """Move the midpoints of chosen bridges far outside every zone."""
    moved: set[int] = set()
    original = Bridge.test_point.fget
    monkeypatch.setattr(
        Bridge, "test_point", property(lambda b: FAR if id(b) in moved else original(b))
    )
    return moved


def test_a_pip_orbit_bridge_outside_by_its_midpoint_is_iterated(henon_p3_session, outside):
    session, _fp3, fp1, zone = henon_p3_session
    session.blast_zone(zone, 1, fixed_point=fp1, min_separation=MIN_SEPARATION)
    arriving = _arriving(session, zone)
    assert arriving, "the blast must leave a bridge arriving at the deepest orbit point"
    outside.update(map(id, arriving))
    assert not any(zone.contains_point(b.test_point) for b in arriving)

    result = session.blast_zone(zone, 1, fixed_point=fp1, min_separation=MIN_SEPARATION)
    assert {id(b) for b in arriving} <= _iterated_first(result)
    assert all(b.iterated for b in arriving)


def test_the_zones_own_pip_decides_not_the_chosen_one(henon_p3_session, outside):
    session, _fp3, fp1, zone = henon_p3_session
    arriving = _arriving(session, zone)
    assert arriving
    outside.update(map(id, arriving))
    session.forget_strong_pip(fp1)

    result = session.blast_zone(zone, 1, fixed_point=fp1, min_separation=MIN_SEPARATION)
    assert {id(b) for b in arriving} <= _iterated_first(result)


def test_the_bridge_arriving_at_the_zones_corner_is_not_exempt(henon_p3_session, outside):
    # A zone cut at f(q0) while the chosen pip stays q0: the arc arriving at
    # f(q0) is that zone's boundary, so its midpoint decides.
    session, _fp3, fp1, zone = henon_p3_session
    orbit = _orbit(session, fp1, zone.boundary_intersection)
    registry = session.workbench.intersection_registry
    shifted = replace(zone, boundary_intersection=registry[orbit[1]])
    corner = _arriving(session, zone, at=slice(1, 2))
    assert corner and registry.find(session.chosen_strong_pip(fp1)) == orbit[0]
    outside.update(map(id, corner))

    result = session.blast_zone(shifted, 1, fixed_point=fp1, min_separation=MIN_SEPARATION)
    assert not {id(b) for b in corner} & _iterated_first(result)
    assert not any(b.iterated for b in corner)


def test_another_fixed_points_pip_orbit_bridge_keeps_the_innermost_zone_rule(henon_p3_session):
    session, fp3, _fp1, zone = henon_p3_session
    inner = session.resonance_zones[fp3]
    arriving = [b for b in _arriving(session, inner) if inner.contains_point(b.test_point)]
    assert arriving, "the inner tangle must hold bridges arriving on its own pip orbit"

    result = session.blast_zone(zone, 1, fixed_point=None, min_separation=MIN_SEPARATION)
    assert not {id(b) for b in arriving} & {id(b) for b in result.interior_bridges_by_iteration[0]}
    assert not any(b.iterated for b in arriving)


def test_a_bridge_arriving_at_an_orbit_point_registered_mid_blast_is_the_zones(henon_p3_session, monkeypatch):
    # Every bridge arriving on the zone's registered orbit, as it stands when
    # the midpoint is read, is moved out: a child arriving at an orbit point
    # the first iteration links is kept only by the orbit re-read.
    session, _fp3, fp1, zone = henon_p3_session
    registry = lambda: session.workbench.intersection_registry  # noqa: E731 - replaced by every recompute

    def registered_orbit() -> list[int]:
        orbit = [registry().find(zone.boundary_intersection)]
        while (image := registry().iterate_table[orbit[-1], 1]) is not None:
            orbit.append(image)
        return orbit

    original = Bridge.test_point.fget
    monkeypatch.setattr(Bridge, "test_point", property(
        lambda b: FAR if b.fixed_point is zone.key and b.id is not None and b.id[1] in registered_orbit()[1:]
        else original(b)
    ))
    before = registered_orbit()

    result = session.blast_zone(zone, 2, fixed_point=fp1, min_separation=MIN_SEPARATION)
    new_points = set(registered_orbit()) - set(before)
    arriving = [b for b in result.interior_bridges_by_iteration[1] if b.id is not None and b.id[1] in new_points]
    assert new_points, "the first iteration must link a deeper orbit point"
    assert arriving, "the bridge arriving at it is the second iteration's frontier"
    assert all(b.iterated for b in arriving)
