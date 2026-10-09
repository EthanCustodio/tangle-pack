"""The production invariant checks are WIRED into the pipeline, not re-run by tests.

``Trellis.punch_holes`` asserts the two topological invariants itself, once per
call: I2 (``check_bridge_rows_consistent``) on every bridge that carries a hole,
and I1 (``check_holes_share_bridge_side``) on all holes with the trellis's
orientation. Both read the crossing signs and the carried sides (D6); the
geometric readings (a bridge's polyline rows, a propagated hole's measured
side) are flags, one WARNING per kind. ``DualGraph`` runs its structural self-check on construction. These
tests spy on those call sites (both checks are imported inside ``punch_holes``,
so a module-attribute patch reaches them) instead of repeating the checks; the
invariants themselves run on every law case in ``tests/invariants/`` and the
checks are pinned to fire on hand-built violations in
``test_topology_invariants.py``.
"""

from __future__ import annotations

import importlib
import logging
from typing import Any, Callable

import pytest

from cases import Case, build_nested, build_period3
from helpers.logs import assert_logged

stable_partition = importlib.import_module("tanglepack.topology.StablePartition")
dual_graph_module = importlib.import_module("tanglepack.topology.DualGraph")
_LOGGER = "tanglepack.topology.StablePartition"


@pytest.fixture
def p3() -> Case:
    """The period-3 orbit alone, partitioned (12 holes in two orbits)."""
    return build_period3()


def _spy(
    monkeypatch: pytest.MonkeyPatch, owner: Any, name: str
) -> list[tuple[tuple, dict]]:
    """
    Wrap ``owner.name`` so every call is recorded and then forwarded.

    Args:
        monkeypatch: pytest's monkeypatch fixture (undoes the wrap).
        owner: The module or class holding the callable.
        name: The attribute to wrap.

    Returns:
        The call log, one ``(args, kwargs)`` per call.
    """
    real: Callable = getattr(owner, name)
    calls: list[tuple[tuple, dict]] = []

    def spy(*args: Any, **kwargs: Any) -> Any:
        calls.append((args, kwargs))
        return real(*args, **kwargs)

    monkeypatch.setattr(owner, name, spy)
    return calls


def test_punch_holes_checks_every_hole_bridge_once_and_every_orbit(
    p3: Case, monkeypatch: pytest.MonkeyPatch
) -> None:
    """One ``punch_holes`` call runs I2 once per hole-bearing bridge and I1 once on
    every hole with the trellis's orientation."""
    (fp,) = p3.fixed_points
    trellis = p3.session.trellis(fp)
    row_calls = _spy(monkeypatch, stable_partition, "check_bridge_rows_consistent")
    side_calls = _spy(monkeypatch, stable_partition, "check_holes_share_bridge_side")

    trellis.punch_holes()

    assert trellis.holes, "the period-3 case must punch holes"
    hole_bridges = {
        frozenset(hole.bounding_ids)
        for hole in trellis.holes
        if hole.bounding_ids is not None
    }
    checked = [
        frozenset((bridge.first_intersection, bridge.second_intersection))
        for (_trellis, bridge), _kwargs in row_calls
    ]
    assert all(args[0] is trellis for args, _kwargs in row_calls)
    assert len(checked) == len(set(checked)), "a hole bridge was checked twice"
    assert set(checked) == hole_bridges

    assert len(side_calls) == 1
    (holes,), kwargs = side_calls[0]
    assert list(holes) == list(trellis.holes)
    assert kwargs == {"orientation_preserving": trellis.orientation_preserving}


def test_punch_holes_raises_on_a_corrupted_hole_side(
    p3: Case, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A propagated hole whose side disagrees with its orbit stops ``punch_holes``."""
    (fp,) = p3.fixed_points
    trellis = p3.session.trellis(fp)
    real = stable_partition.propagate_reference_holes

    def corrupted(*args: Any, **kwargs: Any) -> list:
        holes = real(*args, **kwargs)
        assert holes, "the period-3 case propagates holes"
        holes[0].bridge_side = "left" if holes[0].bridge_side == "right" else "right"
        return holes

    monkeypatch.setattr(stable_partition, "propagate_reference_holes", corrupted)
    with pytest.raises(AssertionError):
        trellis.punch_holes()


def test_propagated_holes_carry_their_orbit_side(p3: Case) -> None:
    """D2 (author, 2026-10-07): a direct hole reads its side off its crossing
    signs, a propagated one carries its orbit's (parity-flipped per step under an
    orientation-reversing map) and keeps the curves' reading as
    ``measured_side``, which agrees wherever it is defined."""
    (fp,) = p3.fixed_points
    trellis = p3.session.trellis(fp)
    flips = not trellis.orientation_preserving
    direct = {hole.origin: hole for hole in trellis.holes if hole.pair is not None}
    for hole in trellis.holes:
        if hole.pair is not None:
            assert hole.side_source == "signs"
            continue
        assert hole.side_source == "carried"
        source = direct[hole.origin]
        expected = source.bridge_side
        if flips and (source.iterate - hole.iterate) % 2:
            expected = "left" if expected == "right" else "right"
        assert hole.bridge_side == expected
        assert hole.measured_side in (None, hole.bridge_side)


def test_punch_holes_flags_a_contradicting_measured_side(
    p3: Case, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """The measured side is the geometric cross-check of the carried one: a
    disagreement is one WARNING, never an error (D6)."""
    (fp,) = p3.fixed_points
    trellis = p3.session.trellis(fp)
    real = stable_partition.propagate_reference_holes

    def contradicted(*args: Any, **kwargs: Any) -> list:
        holes = real(*args, **kwargs)
        hole = next(hole for hole in holes if hole.measured_side is not None)
        hole.measured_side = "left" if hole.measured_side == "right" else "right"
        return holes

    monkeypatch.setattr(stable_partition, "propagate_reference_holes", contradicted)
    caplog.set_level(logging.WARNING, logger=_LOGGER)
    trellis.punch_holes()
    assert_logged(caplog, logging.WARNING, _LOGGER, count=1)


def test_punch_holes_raises_on_equal_signs(
    p3: Case, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A hole bridge on one stable branch whose two crossing signs are equal would
    cross that branch between consecutive crossings: the row assertion fires."""
    (fp,) = p3.fixed_points
    trellis = p3.session.trellis(fp)
    near, far = next(hole.bounding_ids for hole in trellis.holes if hole.pair is not None)
    monkeypatch.setattr(
        trellis.intersection(far), "crossing_sign", trellis.intersection(near).crossing_sign
    )
    with pytest.raises(AssertionError):
        trellis.punch_holes()


@pytest.mark.parametrize("build", [build_period3, build_nested], ids=["p3", "nested"])
@pytest.mark.parametrize("sign", [abs, lambda area: -area], ids=["all_positive", "all_flipped"])
def test_outer_faces_do_not_read_areas(build, sign, monkeypatch: pytest.MonkeyPatch) -> None:
    """INT-10: a component's outer face comes from the rotation system, so faking
    every traversal area (all positive, or every sign flipped) leaves the merged
    face nodes unchanged."""

    def face_nodes(case: Case) -> list[list[tuple[int, ...]]]:
        dual = case.session.dual_graph(rebuild=True)
        return sorted(sorted(face.corners for face in node.faces) for node in dual.face_nodes)

    case = build()
    expected = face_nodes(case)
    real = dual_graph_module.DualGraph._traversal_area
    monkeypatch.setattr(
        dual_graph_module.DualGraph, "_traversal_area", lambda self, face: sign(real(self, face))
    )
    assert face_nodes(case) == expected


def test_dual_graph_construction_runs_its_self_check(
    p3: Case, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Building the dual graph runs its structural self-check exactly once."""
    calls = _spy(monkeypatch, dual_graph_module.DualGraph, "_check")
    graph = p3.session.dual_graph()
    assert len(calls) == 1
    assert calls[0][0][0] is graph
