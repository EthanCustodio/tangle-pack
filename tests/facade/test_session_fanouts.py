"""The session's per-fixed-point fan-outs and its workbench delegation.

A fan-out called with no fixed point runs on every fixed point of the session
and returns a ``{fixed_point: result}`` dict; called with one fixed point it
returns that trellis's own list. The describe fan-outs concatenate the
per-trellis reports. Shapes are checked on the two-fixed-point nested session
(stopped at its resonance zones) and on the single k10 saddle; the contents
of each result are the per-trellis algorithms' business, checked elsewhere.
"""

from __future__ import annotations

from cases import build_k10, build_nested
from tanglepack import TangleSession


def test_no_argument_fanouts_return_one_entry_per_fixed_point() -> None:
    """classify / pseudoneighbors / holes / partition, each a dict over both tangles."""
    case = build_nested(outer_blasts=0, p1_area_cutoff=1e-4, through="zones")
    session = case.session
    fixed_points = set(case.fixed_points)
    assert len(fixed_points) == 2

    candidates = session.classify_strong_pips()
    assert set(candidates) == fixed_points
    for fp in fixed_points:
        assert candidates[fp], f"{fp} should have strong-pip candidates"
        assert session.strong_pip(fp) in candidates[fp]

    references = session.compute_pseudoneighbors()
    assert set(references) == fixed_points
    for fp in fixed_points:
        assert references[fp] == session.trellis(fp).reference_pseudoneighbors

    holes = session.punch_holes()
    assert set(holes) == fixed_points
    for fp in fixed_points:
        assert holes[fp] == session.trellis(fp).holes

    partitions = session.partition_stable_manifold()
    assert set(partitions) == fixed_points
    for fp in fixed_points:
        assert partitions[fp] == session.trellis(fp).stable_partitions
        assert partitions[fp], "a partitioned branch yields one result per side"


def test_single_fixed_point_fanouts_return_that_trellis_list() -> None:
    """Passing one fixed point returns its own list, not a dict."""
    case = build_k10(through="bridges")
    session, fp = case.session, case.fixed_point

    candidates = session.classify_strong_pips(fp)
    assert isinstance(candidates, list)
    assert candidates == session.strong_pip_candidates(fp)

    references = session.compute_pseudoneighbors(fp)
    assert isinstance(references, list) and references
    assert all(pair.is_reference for pair in references)

    holes = session.punch_holes(fp)
    assert isinstance(holes, list) and holes

    results = session.partition_stable_manifold(fp)
    assert isinstance(results, list)
    assert sorted(result.side for result in results) == ["left", "right"]


def test_describe_fanouts_mention_every_fixed_point() -> None:
    """The describe fan-outs hold every per-trellis report under its fixed point."""
    case = build_nested(outer_blasts=0, p1_area_cutoff=1e-4, through="zones")
    session = case.session
    session.classify_strong_pips()
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()

    holes_report = session.describe_holes()
    partition_report = session.describe_stable_partitions()
    for fp in case.fixed_points:
        assert repr(fp) in holes_report
        assert session.trellis(fp).describe_holes() in holes_report
        assert repr(fp) in partition_report
        assert session.trellis(fp).describe_stable_partitions() in partition_report


def test_session_delegates_unknown_attributes_to_the_workbench() -> None:
    """Workbench attributes and drivers are reachable on the session itself."""
    case = build_k10(through="bridges")
    session, workbench = case.session, case.workbench

    assert isinstance(session, TangleSession)
    assert session.generation == workbench.generation
    assert session.intersection_registry is workbench.intersection_registry
    for name in ("grow_until", "grow_until_iterates_closed", "grow_until_faces_closed"):
        assert getattr(session, name) == getattr(workbench, name)
