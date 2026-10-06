"""
The dual graph of a minimal trellis over the iterated homotopy partition.

Pins composition of iterates (a Trellis feature the fundamental segment
relies on), face merging, the element names on the stable nodes, the
fundamental-segment API (no pip, a moved pip, a repeated pip, two pips on one
branch) and the constructor's rejections. The structural laws -- wall and
unified nodes, face sides against the geometry, the pip segment on the pip's
own branch, the bipartite degree rule -- run on every law case in
``tests/invariants/test_law_dual_graph.py``.
"""

from __future__ import annotations

import logging

import pytest

from helpers.logs import assert_logged
from tanglepack import TangleSession
from tanglepack.topology.DualGraph import DualGraph


# --------------------------------------------------------------------------- #
# Builders
# --------------------------------------------------------------------------- #
def _graph_with_pips(session: TangleSession, strong_pips) -> DualGraph:
    """A dual graph over the session's minimal trellis and iterated partition,
    with the strong pips given here instead of the session's own."""
    return DualGraph(
        session.minimal_trellis(), session.iterated_partition(), strong_pips=strong_pips
    )


def _edge_keys(dual: DualGraph) -> dict[tuple, list]:
    """Every stable edge's nodes, keyed by edge (public ``stable_nodes`` only)."""
    edges: dict[tuple, list] = {}
    for node in dual.stable_nodes.values():
        edges.setdefault(node.edge_key, []).append(node)
    return edges


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
        for edge_key, (node, *_rest) in _edge_keys(dual).items()
        if node.branch_key == branch_key
        and node.hi_cdist > boundary + tol
        and node.lo_cdist < top - tol
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
def test_k10_face_nodes_are_one_per_face_with_a_single_outer(k10_partitioned):
    """The flat fixture has one component, so nothing merges: every face is a node."""
    session, fp = k10_partitioned
    dual = session.dual_graph()
    arrangement = dual.arrangement
    assert len(dual.face_nodes) == len(arrangement.faces)
    outer = [face for face in dual.face_nodes if face.kind == "outer"]
    assert len(outer) == 1 and outer[0].is_unbounded
    for region in arrangement.regions:
        assert dual.face_of(region).is_region
    with pytest.raises(ValueError):
        dual.face_of(session.arrangement().faces[0])


def test_p3_merges_the_inner_outer_face_into_the_containing_face(p3_partitioned):
    """One outer node per level of nesting: the inner tangle's outer face is glued
    onto the face of the outer tangle that swallows it."""
    session, fp3, fp1 = p3_partitioned
    dual = session.dual_graph()
    merged = [face for face in dual.face_nodes if len(face.faces) > 1]
    assert dual.arrangement.component_count >= 2
    assert merged, "the nested fixture must merge at least one pair of faces"
    for face in merged:
        assert face.kind == "outer"


# --------------------------------------------------------------------------- #
# Payload: element names
# --------------------------------------------------------------------------- #
def test_k10_every_stable_node_side_is_named(k10_partitioned):
    """Each faced side carries an ElementName agreeing with the graph's naming."""
    from tanglepack.topology.ElementNaming import ElementName, ElementNaming

    session, _fp = k10_partitioned
    dual = session.dual_graph()
    assert isinstance(dual.naming, ElementNaming)
    for node in dual.stable_nodes.values():
        assert set(node.names) == set(node.sides)
        for side in node.sides:
            name = node.names[side]
            assert isinstance(name, ElementName)
            assert name == dual.naming.name(node.elements[side])
            assert name.side_letter == side[0].upper()
            assert dual.naming.ref_of(name) == node.elements[side]


def test_k10_stable_node_label_joins_the_names(k10_partitioned):
    session, fp = k10_partitioned
    dual = session.dual_graph()
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
    minimal = session.minimal_trellis()
    pips = [session.trellis(fp).strong_pip]
    over_homotopy = DualGraph(minimal, session.homotopy_partition(), strong_pips=pips)
    parentless = DualGraph(
        minimal,
        IteratedHomotopyPartition(
            session.iterated_partition().as_list(), trellis=session.trellis()
        ),
        strong_pips=pips,
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
@pytest.mark.parametrize("pips", [None, [], [None]])
def test_no_pips_warns_and_unifies_nothing(k10_partitioned, caplog, pips):
    session, _fp = k10_partitioned
    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.DualGraph"):
        dual = _graph_with_pips(session, pips)
    assert_logged(caplog, logging.WARNING, "tanglepack.topology.DualGraph")
    assert not dual.unified_nodes and not dual.fundamental_segments
    assert all(len(nodes) == 2 for nodes in _edge_keys(dual).values())


def test_k10_a_different_pip_moves_the_unified_set(k10_partitioned):
    session, fp = k10_partitioned
    trellis = session.trellis(fp)
    default = trellis.strong_pip
    # Registry ids are not reproducible between builds, so the alternative is
    # chosen by stable cdist among the candidates whose k-th iterate is
    # registered (the rule needs it), never by id.
    alternatives = sorted(
        (
            c
            for c in trellis.strong_pip_candidates
            if c != default and trellis.iterate(c, fp.k_value) is not None
        ),
        key=lambda c: trellis.intersection(c).stable_cdist,
    )
    if not alternatives:
        pytest.skip("the fixture has no other strong-pip candidate with a registered iterate")
    before = _graph_with_pips(session, [default])
    after = _graph_with_pips(session, [alternatives[0]])
    assert {n.edge_key for n in after.unified_nodes} == _expected_unified_keys(
        after, trellis, pip=alternatives[0]
    )
    assert {n.edge_key for n in after.unified_nodes} != {
        n.edge_key for n in before.unified_nodes
    }


def test_the_same_pip_twice_counts_once(k10_partitioned):
    """``trellis(fp)`` and ``trellis([fp])`` can both hand in the same pip."""
    session, fp = k10_partitioned
    pip = session.trellis(fp).strong_pip
    once = _graph_with_pips(session, [pip])
    twice = _graph_with_pips(session, [pip, pip])
    assert twice.fundamental_segments == once.fundamental_segments
    assert {n.edge_key for n in twice.unified_nodes} == {n.edge_key for n in once.unified_nodes}


def test_two_pips_on_one_branch_are_rejected(k10_partitioned):
    session, fp = k10_partitioned
    trellis = session.trellis(fp)
    pip = trellis.strong_pip
    branch = trellis.intersection(pip).manifold_b_key
    others = [
        c for c in trellis.strong_pip_candidates
        if c != pip and trellis.intersection(c).manifold_b_key == branch
    ]
    if not others:
        pytest.skip("the fixture has a single strong-pip candidate on the pip's branch")
    with pytest.raises(ValueError):
        _graph_with_pips(session, [pip, others[0]])


# --------------------------------------------------------------------------- #
# Views and errors
# --------------------------------------------------------------------------- #
def test_a_missing_partition_side_is_rejected(k10_partitioned):
    from tanglepack.topology.PartitionFamily import PartitionFamily

    session, fp = k10_partitioned
    left_only = PartitionFamily(
        [r for r in session.iterated_partition().as_list() if r.side == "left"],
        trellis=session.trellis(),
    )
    with pytest.raises(ValueError):
        DualGraph(
            session.minimal_trellis(), left_only, strong_pips=[session.trellis(fp).strong_pip]
        )


def test_stable_node_side_of_rejects_a_foreign_face(k10_partitioned):
    session, fp = k10_partitioned
    dual = session.dual_graph()
    node = dual.wall_nodes[0]
    foreign = next(f for f in dual.face_nodes if f is not node.faces[node.sides[0]])
    with pytest.raises(ValueError):
        node.side_of(foreign)
