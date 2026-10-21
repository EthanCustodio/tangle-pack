"""``n`` single-iteration blasts are one ``n``-iteration blast, drops included.

A child ``min_separation`` drops is remembered on the bridge (``Bridge.dropped_at``),
so a later call does not iterate it, and every step's proximity cloud is every
curve the workbench holds. Each test compares ``blast_zone`` with itself on two
builds of one zone, never with a case number.
"""

from __future__ import annotations

import pytest

from cases import ZONE_MIN_SEPARATION, build_k28, build_nested
from tanglepack.loom.Blast import _pip_orbits


def _p3_session_zone():
    """The ``henon_p3_session`` zone with the outer pip forgotten (nothing exempt)."""
    case = build_nested(outer_blasts=0, p1_area_cutoff=1e-4, through="zones")
    fp1 = case.fixed_points[0]
    case.session.forget_strong_pip(fp1)
    zone = max(case.session.resonance_zones.values(), key=lambda z: z.area)
    return case.session, fp1, zone


def _k28_zone():
    """The k=2.8 zone trimmed at ``f(q0)``, unblasted."""
    case = build_k28(blasts=0, through="zones")
    return case.session, case.fixed_points[0], case.zones[0]


def _blast(session, fp, zone, n, separation):
    """``n`` iterations of ``zone`` in one call."""
    return session.blast_zone(
        zone, n, fixed_point=[fp], min_separation=separation, diagnose=False
    )


def _state(session, fp):
    """``iterated`` per bridge id and the pip depth."""
    return (
        {b.id: b.iterated for b in session.workbench.bridges if b.id is not None},
        [len(orbit) - 1 for orbit in _pip_orbits(session, [fp]).values()],
    )


def _dropped(session):
    """The ids of the bridges a blast dropped."""
    return {b.id for b in session.workbench.bridges if b.dropped_at is not None}


@pytest.mark.parametrize(
    ("build", "separation", "n"),
    [(_p3_session_zone, 0.02, 4), (_k28_zone, ZONE_MIN_SEPARATION, 11)],
    ids=["p3_session_no_pip", "k28"],
)
def test_single_blasts_compose(build, separation, n):
    """``n`` calls of one iteration leave the workbench one ``n``-iteration call leaves.

    The p3 case catches a dropped child iterated by the next call; the k28 case,
    deep enough, also catches a call-start cloud that differs from the one a
    single call carries forward.
    """
    stepwise = build()
    for _ in range(n):
        _blast(*stepwise, 1, separation)
    at_once = build()
    result = _blast(*at_once, n, separation)
    assert result.too_close > 0
    assert _state(*stepwise[:2]) == _state(*at_once[:2])
    assert _dropped(stepwise[0]) == _dropped(at_once[0])


def test_drops_survive_a_recompute():
    """A ``preserve_ids`` recompute and ``rebuild_bridges`` keep every drop.

    Four single blasts, the recompute the session runs after growth or a zone
    (re)definition, one more blast: the same workbench as one five-iteration call.
    """
    stepwise = _k28_zone()
    for _ in range(4):
        _blast(*stepwise, 1, 0.02)
    session, fp, _zone = stepwise
    workbench = session.workbench
    workbench.compute_intersections(workbench.fixed_points, preserve_ids=True, diagnose=False)
    workbench.rebuild_bridges()
    last = _blast(*stepwise, 1, 0.02)
    at_once = _k28_zone()
    _blast(*at_once, 5, 0.02)
    assert _state(session, fp) == _state(*at_once[:2])
    assert last.skipped_dropped > 0
    assert _dropped(session) == _dropped(at_once[0])
