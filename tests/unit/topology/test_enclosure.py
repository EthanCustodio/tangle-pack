"""
The non-crossing rule's markers and enclosure (D13) on a hand-built lobe.

One stable branch with crossings ``0`` (its anchor) to ``5`` and one bridge
``(1, 5)`` above it make two faces: the lobe (the bridge and the stable edges
``1..5``) and the unbounded face. Real ``Arc`` / ``Region`` / ``FaceNode``
objects carry the combinatorics; nothing has geometry. The checks are the
rules themselves: the lobe holds what lies beside its bridge and leaves the
anchor out, two chords that interleave on one ring cross (no routing keeps the
curve simple), and a marker maps one map step along its orbit.
"""

from __future__ import annotations

import pytest

from tanglepack.topology.DualGraph import DualGraph, FaceNode, face_side
from tanglepack.topology.TopologyResults import Arc, Region

from helpers.fakes import bare_fixed_point, stable_key


def _dual(face_nodes: list[FaceNode], markers: dict, fixed_points=()) -> DualGraph:
    """A dual graph over hand-made faces and markers (the constructor needs a trellis)."""
    dual = object.__new__(DualGraph)
    dual.face_nodes = face_nodes
    dual.unbounded = next(face for face in face_nodes if face.is_unbounded)
    dual.trellis = type("Trellis", (), {"fixed_points": list(fixed_points)})()
    dual.__dict__["markers"] = markers
    return dual


@pytest.fixture
def lobe():
    """``(dual, lobe face, stable edges e[i] = (i, i+1), the bridge's edge, markers)``."""
    key = stable_key(bare_fixed_point())
    unstable_key = (key[0], "unstable", 0, 0)
    stable = [Arc("stable", i, i + 1, key) for i in range(5)]
    bridge = Arc("unstable", 1, 5, unstable_key, bridge_id=(1, 5))
    # Clockwise (the face on the right of every half-edge): along the bridge,
    # then back toward the anchor along the stable edges it spans.
    lobe_arcs = [bridge] + [arc.reversed() for arc in reversed(stable[1:])]
    outer_arcs = stable + [bridge.reversed(), stable[0].reversed()]
    lobe_face = FaceNode(0, [Region((1, 5, 4, 3, 2), lobe_arcs, None)], "region")
    outer_face = FaceNode(1, [Region((0, 1, 2, 3, 4, 5, 1), outer_arcs, None)], "outer", is_unbounded=True)
    markers = {
        ("orbit", 0, 0): ("vertex", 0),
        ("hole", 7, 0): ("beside", bridge.edge_key, face_side(bridge)),
    }
    return _dual([lobe_face, outer_face], markers), lobe_face, [arc.edge_key for arc in stable], bridge.edge_key


def test_the_lobe_holds_what_lies_beside_its_bridge_and_not_the_anchor(lobe):
    dual, _face, edges, bridge = lobe
    (held,) = dual.enclosed(edges[1:] + [bridge])
    assert held == {("orbit", 0, 0): False, ("hole", 7, 0): True}


def test_a_curve_along_nothing_holds_nothing(lobe):
    dual, _face, _edges, _bridge = lobe
    (held,) = dual.enclosed([])
    assert not any(held.values())


def _chord(face: FaceNode, a: tuple, b: tuple) -> tuple:
    side = {arc.edge_key: face_side(arc) for arc in face.faces[0].arcs}
    return (face, (a, side[a], b, side[b]))


def test_nested_chords_of_one_ring_keep_the_curve_simple(lobe):
    dual, face, edges, _bridge = lobe
    routings = dual.enclosed([], [_chord(face, edges[1], edges[4]), _chord(face, edges[2], edges[3])])
    assert routings is not None and len(routings) == 1


def test_interleaved_chords_of_one_ring_cross(lobe):
    """Chords whose ends alternate around the ring cross: no routing keeps ``f(b)`` simple."""
    dual, face, edges, _bridge = lobe
    assert dual.enclosed([], [_chord(face, edges[1], edges[3]), _chord(face, edges[2], edges[4])]) == []


@pytest.mark.parametrize(
    "ends, crossing",
    [((1, 3, 2, 4), True), ((1, 4, 2, 3), False), ((1, 3, 3, 4), False), ((1, 2, 3, 4), False)],
    ids=["interleaved", "nested", "shared_end_edge", "apart"],
)
def test_chords_cross_iff_their_ends_alternate_on_one_ring(lobe, ends, crossing):
    """Two disjoint simple arcs in a disk face never interleave; chords sharing an end
    edge can always be drawn apart (the order along that edge is free). The verdict is
    the enclosure's: crossing chords leave no routing that keeps the curve simple."""
    dual, face, edges, _bridge = lobe
    a, b, c, d = (edges[i] for i in ends)
    (_f, chord), (_g, other) = _chord(face, a, b), _chord(face, c, d)
    assert dual.chords_cross(face, chord, other) is crossing
    assert dual.chords_cross(face, other, chord) is crossing
    if len({a, b, c, d}) == 4:
        assert (dual.enclosed([], [(face, chord), (face, other)]) == []) is crossing


@pytest.mark.parametrize("period, inversion", [(1, False), (1, True), (3, False), (2, True)])
def test_a_periodic_point_marker_maps_along_its_orbit(period, inversion):
    """The marker of orbit point ``i`` maps to the point ``FixedPoint.advance_key`` names."""
    fixed_point = bare_fixed_point(period, inversion=inversion)
    dual = _dual([FaceNode(0, [], "outer", is_unbounded=True)], {}, [fixed_point])
    for orbit in range(period):
        image = fixed_point.advance_key(stable_key(fixed_point, orbit), 1)
        assert dual.image_marker(("orbit", 0, orbit)) == ("orbit", 0, image[2])


def test_a_hole_marker_maps_to_its_next_iterate_only_when_placed():
    markers = {("hole", 3, 0): ("vertex", 1), ("hole", 3, 1): ("vertex", 2)}
    dual = _dual([FaceNode(0, [], "outer", is_unbounded=True)], markers)
    assert dual.image_marker(("hole", 3, 0)) == ("hole", 3, 1)
    assert dual.image_marker(("hole", 3, 1)) is None
