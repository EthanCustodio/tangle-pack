"""
The dual graph of a minimal trellis over the iterated homotopy partition.

Pins, in the order they can break: composition of iterates (a Trellis
feature the fundamental segment relies on), the node structure (two open
side nodes per stable edge, one solid node on the fundamental segment),
handedness (each face attaches to the node on ITS side), the payload
(elements, parents, cutting bridges, images), and the fundamental segment.
"""

from __future__ import annotations

import dataclasses
import logging

import numpy as np
import pytest

from minimal_helpers import build_pieces
from tanglepack.topology.DualGraph import DualGraph, FaceNode, StableNode
from tanglepack.topology.StablePartition import _side_of
from tanglepack.topology.TopologyResults import ElementRef


# --------------------------------------------------------------------------- #
# Builders
# --------------------------------------------------------------------------- #
def _k10_graph(session, fp) -> tuple[DualGraph, object]:
    """The dual graph of the k=10 fixture, with its one strong pip."""
    pieces = build_pieces(session, [fp])
    return DualGraph(pieces.minimal, pieces.iterated, strong_pips=pieces.strong_pips), pieces


def _p3_graph(session, fp3, fp1) -> tuple[DualGraph, object]:
    """The dual graph of the nested period-3 fixture, with both pips."""
    pieces = build_pieces(session, [fp3, fp1])
    return DualGraph(pieces.minimal, pieces.iterated, strong_pips=pieces.strong_pips), pieces


def _anchorward_tangent(node: StableNode, trellis) -> "np.ndarray | None":
    """The unit direction toward the anchor at an edge's midpoint vertex."""
    poly = node.arc.polyline(trellis)
    if len(poly) < 2:
        return None
    middle = len(poly) // 2
    before = poly[max(middle - 1, 0)]
    after = poly[min(middle + 1, len(poly) - 1)]
    tangent = np.asarray(before, dtype=np.float64) - np.asarray(after, dtype=np.float64)
    norm = float(np.linalg.norm(tangent))
    return tangent / norm if norm > 0.0 else None


def _face_side_checks(dual: DualGraph) -> tuple[int, list[str]]:
    """Compare the combinatorial face side with the geometric one on every region."""
    trellis = dual.arrangement.trellis
    checks = 0
    failures: list[str] = []
    for region in dual.arrangement.regions:
        face = dual.face_of(region)
        inside = region.representative_point
        if inside is None:
            continue
        for arc in region.stable_arcs:
            side = "right" if arc.reverse else "left"
            node = dual.nodes_of_edge(arc.edge_key)[side]
            look = _anchorward_tangent(node, trellis)
            poly = node.arc.polyline(trellis)
            if look is None:
                continue
            midpoint = poly[len(poly) // 2]
            geometric = _side_of(look, np.asarray(inside) - midpoint)
            combinatorial = node.side_of(face)
            checks += 1
            if geometric != combinatorial or node.faces[side] is not face:
                failures.append(
                    f"edge {arc.lo_id}-{arc.hi_id} of region {region.corners}: "
                    f"geometry says {geometric}, combinatorics says {combinatorial}"
                )
    return checks, failures


def _checkable_pairs(dual: DualGraph) -> int:
    """Every (closed region, stable edge on its boundary) pair — none may be skipped."""
    return sum(
        len(region.stable_arcs)
        for region in dual.arrangement.regions
        if region.representative_point is not None
    )


def _expected_unified_keys(dual: DualGraph, trellis, pip=None) -> set[tuple]:
    """The edge keys the pip rule says must be unified on the pip's branch."""
    pip = trellis.strong_pip if pip is None else pip
    branch_key = trellis.intersection(pip).manifold_b_key
    fixed_point = branch_key[0]
    image = trellis.iterate(pip, fixed_point.k_value)
    assert image is not None, (
        "the fixture must have the pip's k-th iterate registered (link by link)"
    )
    boundary = float(dual.trellis.intersection(image).stable_cdist)
    top = float(dual.trellis.intersection(pip).stable_cdist)
    tol = dual.trellis.registry.cdist_tol
    return {
        edge_key
        for edge_key, by_side in dual._nodes_of_edge.items()
        if by_side["left"].branch_key == branch_key
        and by_side["left"].hi_cdist > boundary + tol
        and by_side["left"].lo_cdist < top - tol
    }


# --------------------------------------------------------------------------- #
# Pre-task: Trellis.iterate composes single steps
# --------------------------------------------------------------------------- #
def test_iterate_composes_single_steps_when_the_table_has_no_direct_entry(
    k10_partitioned,
):
    """``iterate(q0, 2)`` answers from two ``n = 1`` links, which is all inference
    records."""
    session, fp = k10_partitioned
    trellis = session.trellis(fp)
    table = trellis.registry.iterate_table
    q0 = trellis.strong_pip

    first = trellis.iterate(q0, 1)
    assert first is not None
    second = trellis.iterate(first, 1)
    assert second is not None, "the k=10 fixture must have a two-link chain"
    assert table[q0, 2] is None, (
        "this test is only meaningful while the table stores no direct 2-iterate"
    )
    assert trellis.iterate(q0, 2) == second
    assert trellis.iterate(q0, 0) == q0


def test_iterate_returns_none_as_soon_as_a_link_is_missing(k10_partitioned):
    session, fp = k10_partitioned
    trellis = session.trellis(fp)
    q0 = trellis.strong_pip

    walk = q0
    steps = 0
    while walk is not None:
        walk = trellis.iterate(walk, 1)
        steps += 1
        if steps > 20:
            pytest.skip("the chain does not terminate within 20 steps")
    assert trellis.iterate(q0, steps) is None


@pytest.mark.slow
def test_p3_iterate_composes_three_steps(p3_partitioned):
    """On the period-3 tangle ``iterate(q0, 3)`` walks the whole orbit and back."""
    session, fp3, _fp1 = p3_partitioned
    trellis = session.trellis(fp3)
    q0 = trellis.strong_pip

    walked = q0
    for _ in range(3):
        walked = trellis.iterate(walked, 1)
        assert walked is not None
    assert trellis.iterate(q0, 3) == walked

    orbit = trellis.registry.iterate_orbit(q0)
    assert orbit == [trellis.iterate(q0, n) for n in range(len(orbit))], (
        "the composed iterates must reproduce the registry's own forward orbit"
    )

    assert trellis.iterate(walked, -3) == q0, (
        "three backward steps must undo three forward ones"
    )
    if trellis.iterate(q0, -1) is None:
        assert trellis.iterate(q0, -3) is None, (
            "a missing first backward link must stop the composition"
        )


# --------------------------------------------------------------------------- #
# Node structure
# --------------------------------------------------------------------------- #
def _check_structure(dual: DualGraph, pieces) -> None:
    edges = {
        arc.edge_key
        for face in dual.arrangement.faces
        for arc in face.arcs
        if arc.kind == "stable"
    }
    assert set(dual._nodes_of_edge) == edges
    for edge_key, by_side in dual._nodes_of_edge.items():
        left, right = by_side["left"], by_side["right"]
        if left is right:
            assert left.is_unified and left.traversable and left.node_side == "both"
            assert left.key == (edge_key, "both") and set(left.faces) == {"left", "right"}
            assert left.fundamental_span is not None
        else:
            for node, side in ((left, "left"), (right, "right")):
                assert node.sides == (side,) and not node.traversable
                assert node.key == (edge_key, side) and list(node.faces) == [side]
                assert node.fundamental_span is None
                with pytest.raises(ValueError, match="wall"):
                    node.other_face(node.faces[side])
                with pytest.raises(ValueError, match="faces"):
                    node.element_on("right" if side == "left" else "left")
    for node in dual.stable_nodes.values():
        for side in node.sides:
            face = node.face_on(side)
            assert (node, side) in face.stable_nodes
            assert face.side_of(node) == side or node.is_unified
        assert dual.stable_nodes[node.key] is node
    for face in dual.face_nodes:
        for node, side in face.stable_nodes:
            assert node.faces[side] is face
        for node in face.exits:
            assert node.is_unified and node.other_face(face) is node.face_on(
                "right" if node.side_of(face) == "left" else "left"
            )
        assert set(face.corners) == {
            corner for region in face.faces for corner in region.corners if corner >= 0
        }
        assert set(face.bridge_ids) == {
            bid for region in face.faces for bid in region.bridge_ids
        }
    assert dual.unbounded.is_unbounded and dual.unbounded.kind == "outer"
    assert sum(1 for face in dual.face_nodes if face.is_unbounded) == 1
    assert set(dual.unified_nodes) | set(dual.wall_nodes) == set(dual.stable_nodes.values())
    assert dual.minimal is pieces.minimal and dual.partition is pieces.iterated
    assert dual.arrangement is pieces.minimal.arrangement


def test_k10_node_structure(k10_partitioned):
    session, fp = k10_partitioned
    dual, pieces = _k10_graph(session, fp)
    _check_structure(dual, pieces)
    assert dual.unified_nodes and dual.wall_nodes
    assert len(dual.stable_nodes) == 2 * len(dual._nodes_of_edge) - len(dual.unified_nodes)


@pytest.mark.slow
def test_p3_node_structure(p3_partitioned):
    session, fp3, fp1 = p3_partitioned
    dual, pieces = _p3_graph(session, fp3, fp1)
    _check_structure(dual, pieces)


def test_k10_face_nodes_are_one_per_face_with_a_single_outer(k10_partitioned):
    """The flat fixture has one component, so nothing merges: every face is a node."""
    session, fp = k10_partitioned
    dual, _ = _k10_graph(session, fp)
    arrangement = dual.arrangement
    assert len(dual.face_nodes) == len(arrangement.faces)
    outer = [face for face in dual.face_nodes if face.kind == "outer"]
    assert len(outer) == 1 and outer[0].is_unbounded
    for region in arrangement.regions:
        assert dual.face_of(region).is_region
    with pytest.raises(ValueError, match="not a face"):
        dual.face_of(session.arrangement().faces[0])


@pytest.mark.slow
def test_p3_merges_the_inner_outer_face_into_the_containing_face(p3_partitioned):
    """One outer node per level of nesting: the inner tangle's outer face is glued
    onto the face of the outer tangle that swallows it."""
    session, fp3, fp1 = p3_partitioned
    dual, _ = _p3_graph(session, fp3, fp1)
    merged = [face for face in dual.face_nodes if len(face.faces) > 1]
    assert dual.arrangement.component_count >= 2
    assert merged, "the nested fixture must merge at least one pair of faces"
    for face in merged:
        assert face.kind == "outer"


# --------------------------------------------------------------------------- #
# Handedness
# --------------------------------------------------------------------------- #
def test_k10_face_side_agrees_with_the_geometry(k10_partitioned):
    session, fp = k10_partitioned
    dual, _ = _k10_graph(session, fp)
    checks, failures = _face_side_checks(dual)
    assert not failures, "\n".join(failures)
    assert checks == _checkable_pairs(dual)


@pytest.mark.slow
def test_p3_face_side_agrees_with_the_geometry(p3_partitioned):
    session, fp3, fp1 = p3_partitioned
    dual, _ = _p3_graph(session, fp3, fp1)
    checks, failures = _face_side_checks(dual)
    assert not failures, "\n".join(failures)
    assert checks == _checkable_pairs(dual)


# --------------------------------------------------------------------------- #
# Payload
# --------------------------------------------------------------------------- #
def _check_payload(dual: DualGraph, pieces) -> None:
    iterated, homotopy, full = pieces.iterated, pieces.homotopy, pieces.full
    partitions = iterated.as_list()
    flip = not full.orientation_preserving
    for node in dual.stable_nodes.values():
        for side in node.sides:
            ref = node.element_on(side)
            interval = node.intervals[side]
            assert ref == iterated.element_at(node.branch_key, side, node.mid_cdist)
            assert interval is iterated.element(ref)
            assert node.parents[side] == ElementRef(
                node.branch_key, side, interval.parent_element_id
            )
            assert homotopy.element(node.parents[side]).lo_cdist <= node.lo_cdist + full.registry.cdist_tol
            assert node.cut_by[side] == interval.cut_by
            expected = full.image_of_element(
                iterated.result(node.branch_key, side), ref.element_id, 1,
                partitions=partitions,
            )
            image_key = node.branch_key[0].advance_key(node.branch_key, 1)
            image_side = ("right" if side == "left" else "left") if flip else side
            assert node.images[side] == [
                ElementRef(image_key, image_side, i) for i in (expected or [])
            ]
        assert dual.stable_nodes_of(node.elements[node.sides[0]]).count(node) == 1
    # No face carries a bridge class.
    for name in (f.name for f in dataclasses.fields(FaceNode)):
        assert "class" not in name and "letter" not in name
    # Face images resolve lazily, to a node of this graph or None.
    for face in dual.face_nodes:
        image = dual.image_face(face)
        assert image is None or image in dual.face_nodes
        if not face.is_region:
            assert image is None


def test_k10_payload(k10_partitioned):
    session, fp = k10_partitioned
    dual, pieces = _k10_graph(session, fp)
    _check_payload(dual, pieces)
    assert any(node.cut_by[s] is not None for node in dual.stable_nodes.values() for s in node.sides)
    assert any(node.images[s] for node in dual.stable_nodes.values() for s in node.sides)


@pytest.mark.slow
def test_p3_payload(p3_partitioned):
    session, fp3, fp1 = p3_partitioned
    dual, pieces = _p3_graph(session, fp3, fp1)
    _check_payload(dual, pieces)


# --------------------------------------------------------------------------- #
# Payload: element names
# --------------------------------------------------------------------------- #
#: The iterated-element names of the k=10 fixture (planning prototype, pinned
#: by the symbolic-dynamics plan): three children of R_1 and R_3, two of L_1.
K10_ITERATED_NAMES = {
    "L_1^1", "L_1^2", "L_2", "L_3",
    "R_1^1", "R_1^2", "R_1^3", "R_2", "R_3^1", "R_3^2", "R_3^3",
}


def test_k10_every_stable_node_side_is_named(k10_partitioned):
    """Each faced side carries an ElementName agreeing with the graph's naming."""
    from tanglepack.topology.ElementNaming import ElementName, ElementNaming

    session, fp = k10_partitioned
    dual, pieces = _k10_graph(session, fp)
    assert isinstance(dual.naming, ElementNaming)
    for node in dual.stable_nodes.values():
        assert set(node.names) == set(node.sides)
        for side in node.sides:
            name = node.names[side]
            assert isinstance(name, ElementName)
            assert name == dual.naming.name(node.elements[side])
            assert name.side_letter == side[0].upper()
            assert dual.naming.ref_of(name) == node.elements[side]


def test_k10_stable_node_names_are_the_expected_iterated_names(k10_partitioned):
    session, fp = k10_partitioned
    dual, _ = _k10_graph(session, fp)
    seen = {name.text for node in dual.stable_nodes.values() for name in node.names.values()}
    assert seen <= K10_ITERATED_NAMES
    # Every split element with an edge in the sparse arrangement shows up.
    assert any(name.is_split for node in dual.stable_nodes.values() for name in node.names.values())
    assert any("^" not in text for text in seen)


def test_k10_stable_node_label_joins_the_names(k10_partitioned):
    session, fp = k10_partitioned
    dual, _ = _k10_graph(session, fp)
    for node in dual.unified_nodes:
        assert node.label == f"{node.names['left'].text} | {node.names['right'].text}"
        assert " | " in node.label
    for node in dual.wall_nodes:
        assert node.label == node.names[node.sides[0]].text
        assert "|" not in node.label
    # With no name on a side the label falls back to the element label.
    node = dual.wall_nodes[0]
    saved = dict(node.names)
    try:
        node.names.clear()
        assert node.label == node.elements[node.sides[0]].label
    finally:
        node.names.update(saved)


def test_k10_homotopy_only_partition_names_plain(k10_partitioned):
    """Over a family with no homotopy parent every element is its own parent."""
    from tanglepack.topology.PartitionFamily import IteratedHomotopyPartition

    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    over_homotopy = DualGraph(pieces.minimal, pieces.homotopy, strong_pips=pieces.strong_pips)
    parentless = DualGraph(
        pieces.minimal,
        IteratedHomotopyPartition(pieces.iterated.as_list(), trellis=pieces.full),
        strong_pips=pieces.strong_pips,
    )
    for dual in (over_homotopy, parentless):
        for node in dual.stable_nodes.values():
            for side in node.sides:
                name = node.names[side]
                assert not name.is_split and "^" not in name.text
                assert name.subscript == node.elements[side].element_id + 1
                assert name == dual.naming.homotopy_name(node.elements[side])


# --------------------------------------------------------------------------- #
# The fundamental segment
# --------------------------------------------------------------------------- #
def test_k10_unified_nodes_are_the_edges_between_the_pip_and_its_image(k10_partitioned):
    session, fp = k10_partitioned
    dual, _ = _k10_graph(session, fp)
    trellis = session.trellis(fp)
    expected = _expected_unified_keys(dual, trellis)
    assert expected, "the k=10 fixture must unify something"
    assert {node.edge_key for node in dual.unified_nodes} == expected
    branch_key = trellis.intersection(trellis.strong_pip).manifold_b_key
    lo, hi = dual.fundamental_segments[branch_key]
    assert hi == pytest.approx(trellis.intersection(trellis.strong_pip).stable_cdist)
    for node in dual.unified_nodes:
        assert node.fundamental_span == (lo, hi)


@pytest.mark.slow
def test_p3_unifies_the_pip_segment_on_each_pips_own_branch(p3_partitioned):
    session, fp3, fp1 = p3_partitioned
    dual, _ = _p3_graph(session, fp3, fp1)
    expected: set[tuple] = set()
    branches = set()
    for fp in (fp3, fp1):
        trellis = session.trellis(fp)
        expected |= _expected_unified_keys(dual, trellis)
        branches.add(trellis.intersection(trellis.strong_pip).manifold_b_key)
    assert {node.edge_key for node in dual.unified_nodes} == expected
    assert set(dual.fundamental_segments) == branches


@pytest.mark.parametrize("pips", [None, [], [None]])
def test_no_pips_warns_and_unifies_nothing(k10_partitioned, caplog, pips):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.DualGraph"):
        dual = DualGraph(pieces.minimal, pieces.iterated, strong_pips=pips)
    assert any("no strong pip" in record.message for record in caplog.records)
    assert not dual.unified_nodes and not dual.fundamental_segments
    assert len(dual.stable_nodes) == 2 * len(dual._nodes_of_edge)


def test_k10_a_different_pip_moves_the_unified_set(k10_partitioned):
    session, fp = k10_partitioned
    trellis = session.trellis(fp)
    default = trellis.strong_pip
    alternatives = [c for c in trellis.strong_pip_candidates if c != default]
    if not alternatives:
        pytest.skip("the fixture has a single strong-pip candidate")
    pieces = build_pieces(session, [fp])
    before = DualGraph(pieces.minimal, pieces.iterated, strong_pips=[default])
    after = DualGraph(pieces.minimal, pieces.iterated, strong_pips=[alternatives[0]])
    assert {n.edge_key for n in after.unified_nodes} == _expected_unified_keys(
        after, trellis, pip=alternatives[0]
    )
    assert {n.edge_key for n in after.unified_nodes} != {
        n.edge_key for n in before.unified_nodes
    }


def test_two_pips_on_one_branch_are_rejected(k10_partitioned):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    pip = session.trellis(fp).strong_pip
    with pytest.raises(ValueError, match="one fundamental segment"):
        DualGraph(pieces.minimal, pieces.iterated, strong_pips=[pip, pip])


# --------------------------------------------------------------------------- #
# Views and errors
# --------------------------------------------------------------------------- #
def test_k10_graph_is_bipartite(k10_partitioned):
    import networkx as nx

    session, fp = k10_partitioned
    dual, _ = _k10_graph(session, fp)
    graph = dual.graph()
    faces = {n for n, d in graph.nodes(data=True) if d["bipartite"] == 0}
    stable = {n for n, d in graph.nodes(data=True) if d["bipartite"] == 1}
    assert nx.is_bipartite(graph)
    assert len(faces) == len(dual.face_nodes) and len(stable) == len(dual.stable_nodes)
    for a, b, data in graph.edges(data=True):
        assert {a[0], b[0]} == {"face", "stable"} and data["side"] in ("left", "right")
    for node in dual.stable_nodes.values():
        assert graph.nodes[("stable", node.key)]["traversable"] == node.traversable
        degree = graph.degree(("stable", node.key))
        assert degree == 1 or (node.is_unified and degree == 2)


def test_k10_summary_reports_counts(k10_partitioned):
    session, fp = k10_partitioned
    dual, _ = _k10_graph(session, fp)
    text = dual.summary()
    assert "stable edges" in text and "solid" in text and "fundamental segments" in text
    assert repr(dual) == f"<{text}>"
    assert repr(dual.unified_nodes[0]).startswith("StableNode(")
    assert "exit(s)" in repr(dual.face_nodes[0])


def test_a_missing_partition_side_is_rejected(k10_partitioned):
    from tanglepack.topology.PartitionFamily import PartitionFamily

    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    left_only = PartitionFamily(
        [r for r in pieces.iterated.as_list() if r.side == "left"], trellis=pieces.full
    )
    with pytest.raises(ValueError, match="no 'right' partition covers"):
        DualGraph(pieces.minimal, left_only, strong_pips=pieces.strong_pips)


def test_stable_node_side_of_rejects_a_foreign_face(k10_partitioned):
    session, fp = k10_partitioned
    dual, _ = _k10_graph(session, fp)
    node = dual.wall_nodes[0]
    foreign = next(f for f in dual.face_nodes if f is not node.faces[node.sides[0]])
    with pytest.raises(ValueError, match="does not bound"):
        node.side_of(foreign)
