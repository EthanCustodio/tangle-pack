"""The session's strong-pip classification over fixed points and its workbench delegation.

``classify_strong_pips`` called with no fixed point classifies every fixed
point of the session and returns a ``{fixed_point: candidates}`` dict; called
with one fixed point it returns that trellis's own list. Shapes are checked on
the two-fixed-point nested session (stopped at its resonance zones) and on the
single k10 saddle; the candidates themselves are the strong-pip algorithm's
business, checked elsewhere.
"""

from __future__ import annotations

from cases import build_k10, build_nested
from tanglepack import TangleSession


def test_classify_strong_pips_covers_every_fixed_point() -> None:
    """With no fixed point, a dict over both tangles."""
    case = build_nested(outer_blasts=0, p1_area_cutoff=1e-4, through="zones")
    session = case.session
    fixed_points = set(case.fixed_points)
    assert len(fixed_points) == 2

    candidates = session.classify_strong_pips()
    assert set(candidates) == fixed_points
    for fp in fixed_points:
        assert candidates[fp], f"{fp} should have strong-pip candidates"
        assert session.strong_pip(fp) in candidates[fp]


def test_classify_strong_pips_of_one_fixed_point_is_its_list() -> None:
    """Passing one fixed point returns its own list, not a dict."""
    case = build_k10(through="bridges")
    session, fp = case.session, case.fixed_point

    candidates = session.classify_strong_pips(fp)
    assert isinstance(candidates, list)
    assert candidates == session.strong_pip_candidates(fp)


def test_session_delegates_unknown_attributes_to_the_workbench() -> None:
    """Workbench attributes and drivers are reachable on the session itself."""
    case = build_k10(through="bridges")
    session, workbench = case.session, case.workbench

    assert isinstance(session, TangleSession)
    assert session.generation == workbench.generation
    assert session.intersection_registry is workbench.intersection_registry
    for name in ("grow_until", "grow_until_iterates_closed", "grow_until_faces_closed"):
        assert getattr(session, name) == getattr(workbench, name)
