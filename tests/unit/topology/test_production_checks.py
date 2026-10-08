"""The production invariant checks are WIRED into the pipeline, not re-run by tests.

``Trellis.punch_holes`` asserts the two topological invariants itself, once per
call: I2 (``check_bridge_rows_consistent``) on every bridge that carries a hole,
and I1 (``check_holes_share_bridge_side``) on all holes with the trellis's
orientation -- since propagated holes carry their side (2026-10-07), I1's teeth
are each propagated hole's measured side. ``DualGraph`` runs its structural self-check on construction. These
tests spy on those call sites (both checks are imported inside ``punch_holes``,
so a module-attribute patch reaches them) instead of repeating the checks; the
invariants themselves run on every law case in ``tests/invariants/`` and the
checks are pinned to fire on hand-built violations in
``test_topology_invariants.py``.
"""

from __future__ import annotations

import importlib
from typing import Any, Callable

import pytest

from cases import Case, build_period3

stable_partition = importlib.import_module("tanglepack.topology.StablePartition")
dual_graph_module = importlib.import_module("tanglepack.topology.DualGraph")


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
    signs, a propagated one carries its orbit's and keeps the curves' reading as
    ``measured_side``, which agrees wherever it is defined."""
    (fp,) = p3.fixed_points
    trellis = p3.session.trellis(fp)
    assert trellis.orientation_preserving
    direct = {hole.origin: hole.bridge_side for hole in trellis.holes if hole.pair is not None}
    propagated = [hole for hole in trellis.holes if hole.pair is None]
    assert propagated, "the period-3 case propagates holes"
    assert {hole.side_source for hole in trellis.holes if hole.pair is not None} == {"signs"}
    for hole in propagated:
        assert hole.side_source == "carried"
        assert hole.bridge_side == direct[hole.origin]
        assert hole.measured_side in (None, hole.bridge_side)


def test_punch_holes_raises_on_a_contradicting_measured_side(
    p3: Case, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The carried side makes the orbit check hold by construction; the measured
    side is what keeps ``punch_holes`` a geometric guard."""
    (fp,) = p3.fixed_points
    trellis = p3.session.trellis(fp)
    real = stable_partition.propagate_reference_holes

    def contradicted(*args: Any, **kwargs: Any) -> list:
        holes = real(*args, **kwargs)
        hole = next(hole for hole in holes if hole.measured_side is not None)
        hole.measured_side = "left" if hole.measured_side == "right" else "right"
        return holes

    monkeypatch.setattr(stable_partition, "propagate_reference_holes", contradicted)
    with pytest.raises(AssertionError):
        trellis.punch_holes()


def test_two_outer_faces_raise_with_every_positive_face(
    p3: Case, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A component with more than one positively winding face keeps raising a
    ``ValueError``, now carrying every positive face (largest first)."""
    real = dual_graph_module.DualGraph._traversal_area
    monkeypatch.setattr(
        dual_graph_module.DualGraph,
        "_traversal_area",
        lambda self, face: abs(real(self, face)),
    )
    with pytest.raises(dual_graph_module.OuterFaceError) as caught:
        p3.session.dual_graph()
    assert isinstance(caught.value, ValueError)
    faces = caught.value.positive_faces
    assert len(faces) > 1
    areas = [area for _index, _corners, area in faces]
    assert areas == sorted(areas, reverse=True) and areas[-1] > 0.0
    assert all(corners for _index, corners, _area in faces)


def test_dual_graph_construction_runs_its_self_check(
    p3: Case, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Building the dual graph runs its structural self-check exactly once."""
    calls = _spy(monkeypatch, dual_graph_module.DualGraph, "_check")
    graph = p3.session.dual_graph()
    assert len(calls) == 1
    assert calls[0][0][0] is graph
