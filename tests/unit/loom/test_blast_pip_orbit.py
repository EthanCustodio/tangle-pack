"""The pip orbit under blasting: ``pip_depth`` and the proximity guard.

The depth ``m`` (``BlastResult.pip_depth``, a report since D12 was reversed
on 2026-10-31) is the deepest registered forward iterate of the pip. The
proximity guard no longer exempts a child ending on that orbit (the B3
exemption went on 2026-10-31): ``m`` cuts nothing, so the guard may set it.
The point ``f^{m+1}(q0)`` is born as an END of a child of the bridge ending at
``f^m(q0)``; that child is iterated in turn only inside a multi-iteration
call, so the tests blast three iterations at once. The child is found by rule
(an end on the composed orbit), never by id.
"""

from __future__ import annotations

MIN_SEPARATION = 0.02


def _pip_orbit(session, fp, chosen) -> list[int]:
    """``[q0, f(q0), ..., f^m(q0)]`` for the pip ``chosen``, composed through ``trellis.iterate``."""
    trellis = session.trellis(fp)
    orbit = [trellis.registry.find(chosen)]
    while (image := trellis.iterate(orbit[-1], 1)) is not None:
        orbit.append(image)
    return orbit


def _blast(session, zone, fp):
    """The fixture's three-iteration blast at the guard's separation."""
    return session.blast_zone(
        zone, num_iterations=3, fixed_point=fp, min_separation=MIN_SEPARATION
    )


def _ending_on(bridges, orbit) -> list:
    """The bridges with an end crossing on ``orbit``."""
    return [b for b in bridges if set(b.id or ()) & set(orbit)]


def test_pip_bearing_child_is_dropped(henon_p3_session):
    """With the pip chosen, a child ending on its orbit closer than the separation is dropped."""
    session, _fp3, fp1, zone = henon_p3_session
    result = _blast(session, zone, fp1)
    orbit = _pip_orbit(session, fp1, session.chosen_strong_pip(fp1))
    assert _ending_on([c for s in result.steps for c in s.discarded_too_close], orbit)


def test_pip_depth_is_the_registered_orbit(henon_p3_session):
    """``BlastResult.pip_depth`` is the depth ``trellis.iterate`` composes after the blast."""
    session, _fp3, fp1, zone = henon_p3_session
    result = _blast(session, zone, fp1)
    orbit = _pip_orbit(session, fp1, session.chosen_strong_pip(fp1))
    assert result.pip_depth == {fp1: len(orbit) - 1}
