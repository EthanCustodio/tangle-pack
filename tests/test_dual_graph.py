"""Phase B — the dual graph: arc nodes, merged face nodes, and the fill.

Four things are pinned here, in the order they can break.

**Composition of iterates** (the pre-task): ``Trellis.iterate`` composes single
steps when the table has no direct entry, so a chain that exists link by link is
answerable end to end.

**The graph's shape**: one arc node per stable edge of every face, a face on both
sides of every arc node, one face node per union-find class, and exactly one
``outer`` node per level of nesting — one on the flat k=10 fixture, two on the
nested period-3 one, where the inner tangle's outer face is merged into the outer
tangle's containing face.

**Handedness**: the face side of an arc is read off ``Arc.reverse`` alone. That is
only right if the arrangement's traversal orientation and the partition's
``left``/``right`` are the same orientation of the plane, so every closed region is
checked against :func:`~tanglepack.topology.StablePartition._side_of` applied to
the real geometry — the region's representative point against the arc's local
anchorward tangent.

**The fill**: the arcs an unseen crossing could land on. With every stable branch
trimmed at its tangle's strong-pip orbit the derived segment is ``(f^k(q0), q0]``
on each pip's own branch and empty on every other branch of that orbit.
"""

from __future__ import annotations

import numpy as np
import pytest

from tanglepack.topology.DualGraph import ArcNode, DualGraph, FaceNode
from tanglepack.topology.StablePartition import _side_of
from tanglepack.topology.TopologyResults import ElementRef


# --------------------------------------------------------------------------- #
# Builders
# --------------------------------------------------------------------------- #
def _k10_graph(session, fp) -> DualGraph:
    """The dual graph of the k=10 fixture, with its one strong pip."""
    trellis = session.trellis(fp)
    return DualGraph(
        session.arrangement(),
        trellis.stable_partitions,
        strong_pips=[trellis.strong_pip],
    )


def _p3_graph(session, fp3, fp1) -> DualGraph:
    """The dual graph of the nested fixture: both tangles' partitions and pips."""
    t3, t1 = session.trellis(fp3), session.trellis(fp1)
    return DualGraph(
        session.arrangement(),
        list(t3.stable_partitions) + list(t1.stable_partitions),
        strong_pips=[t3.strong_pip, t1.strong_pip],
    )


def _anchorward_tangent(node: ArcNode, trellis) -> "np.ndarray | None":
    """The unit direction toward the anchor at an arc's midpoint.

    Read from the two polyline vertices flanking the midpoint vertex rather than
    from ``_stable_frame`` at the arc's lo crossing: the frame has no anchorward
    direction at all at the anchor artifact (there are no manifold nodes below
    cdist 0), which would silently drop the arcs nearest the periodic point from
    the comparison. The polyline runs lo -> hi, i.e. outward, so the anchorward
    direction is the reverse.
    """
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
    trellis = dual.trellis
    checks = 0
    failures: list[str] = []
    for region in dual.arrangement.regions:
        face = dual.face_of(region)
        inside = region.representative_point
        if inside is None:
            continue
        for arc in region.stable_arcs:
            node = dual.arc_nodes[arc.edge_key]
            look = _anchorward_tangent(node, trellis)
            poly = node.arc.polyline(trellis)
            if look is None:
                continue
            midpoint = poly[len(poly) // 2]
            geometric = _side_of(look, np.asarray(inside) - midpoint)
            combinatorial = node.side_of(face)
            checks += 1
            if geometric != combinatorial:
                failures.append(
                    f"arc {arc.lo_id}-{arc.hi_id} of region {region.corners}: "
                    f"geometry says {geometric}, combinatorics says {combinatorial}"
                )
    return checks, failures


def _checkable_pairs(dual: DualGraph) -> int:
    """Every (closed region, stable arc on its boundary) pair — none may be skipped."""
    return sum(
        len(region.stable_arcs)
        for region in dual.arrangement.regions
        if region.representative_point is not None
    )


def _pip_branch_key(trellis):
    """The stable branch key the trellis's strong pip lies on."""
    return trellis.intersection(trellis.strong_pip).manifold_b_key


def _expected_fill_ids(dual: DualGraph, trellis) -> set[tuple]:
    """The arc-node keys the pip rule says must be filled on the pip's branch.

    ``(f^k(q0), q0]``: every arc of the pip's own branch lying at or above the
    canonical distance of the pip's ``k``-th image.
    """
    pip = trellis.strong_pip
    fixed_point = _pip_branch_key(trellis)[0]
    image = trellis.iterate(pip, fixed_point.k_value)
    assert image is not None, (
        "the fixture must have the pip's k-th iterate registered (link by link)"
    )
    boundary = float(dual.trellis.intersection(image).stable_cdist)
    tol = dual.trellis.registry.cdist_tol
    return {
        node.key
        for node in dual.arc_nodes_on(_pip_branch_key(trellis))
        if node.lo_cdist >= boundary - tol
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
# B.1 -- arc nodes
# --------------------------------------------------------------------------- #
def test_k10_every_stable_edge_has_exactly_one_arc_node(k10_partitioned):
    session, fp = k10_partitioned
    dual = _k10_graph(session, fp)

    edge_keys = {
        arc.edge_key
        for face in dual.arrangement.faces
        for arc in face.arcs
        if arc.kind == "stable"
    }
    assert set(dual.arc_nodes) == edge_keys
    assert len(dual.arc_nodes) == len(edge_keys)
    for node in dual.arc_nodes.values():
        assert node.arc.reverse is False
        assert node.faces.keys() == {"left", "right"}
        assert node.other_face(node.face_on("left")) is node.face_on("right")


def test_k10_arc_elements_own_the_midpoint(k10_partitioned):
    session, fp = k10_partitioned
    dual = _k10_graph(session, fp)

    for node in dual.arc_nodes.values():
        middle = 0.5 * (node.lo_cdist + node.hi_cdist)
        assert dual.element_at(node.branch_key, "left", middle) == node.left
        assert dual.element_at(node.branch_key, "right", middle) == node.right
        assert node.element_on("left") == node.left
        assert node in dual.arc_nodes_of(node.left)
        interval = dual.element(node.left)
        assert interval.lo_cdist <= middle <= interval.hi_cdist


def test_element_at_refuses_a_cdist_outside_the_partition(k10_partitioned):
    session, fp = k10_partitioned
    dual = _k10_graph(session, fp)
    (branch_key, side), result = next(iter(dual.partitions.items()))
    beyond = result.intervals[-1].hi_cdist + 10.0

    with pytest.raises(ValueError, match="own stable cdist"):
        dual.element_at(branch_key, side, beyond)


def test_arc_node_side_of_rejects_a_foreign_face(k10_partitioned):
    session, fp = k10_partitioned
    dual = _k10_graph(session, fp)
    node = next(iter(dual.arc_nodes.values()))
    stranger = FaceNode(index=999, faces=[], kind="open")

    with pytest.raises(ValueError, match="does not bound"):
        node.side_of(stranger)


# --------------------------------------------------------------------------- #
# B.2 -- face nodes and the merge rule
# --------------------------------------------------------------------------- #
def test_k10_face_nodes_are_one_per_face_with_a_single_outer(k10_partitioned):
    """The flat fixture has one component, so nothing merges: every face is a node
    and the one open outer face IS the unbounded node."""
    session, fp = k10_partitioned
    dual = _k10_graph(session, fp)
    arrangement = dual.arrangement

    assert len(dual.face_nodes) == len(arrangement.faces)
    assert [node.index for node in dual.face_nodes] == list(range(len(dual.face_nodes)))

    outer = [node for node in dual.face_nodes if node.kind == "outer"]
    assert len(outer) == 1
    assert outer[0] is dual.unbounded
    assert outer[0].is_unbounded
    assert sum(node.is_unbounded for node in dual.face_nodes) == 1
    assert [face.is_closed for face in outer[0].faces] == [False]

    for region in arrangement.regions:
        node = dual.face_of(region)
        assert node.is_region and node.kind == "region"
        assert node.faces == [region]
        assert node.bridge_ids == region.bridge_ids


def test_k10_every_arc_node_has_a_face_on_both_sides(k10_partitioned):
    session, fp = k10_partitioned
    dual = _k10_graph(session, fp)

    slots = 0
    for face in dual.face_nodes:
        for node, side in face.arcs:
            slots += 1
            assert node.face_on(side) is face
            assert face.side_of(node) in ("left", "right")
    assert slots == 2 * len(dual.arc_nodes)


@pytest.mark.slow
def test_p3_merges_the_inner_outer_face_into_the_containing_face(p3_partitioned):
    """One outer node per level of nesting.

    The period-1 tangle's own outer face is the unbounded plane; its CONTAINING
    face (the closed cycle that swallows the whole period-3 tangle) is the piece
    of plane the period-3 tangle's outer face is a part of, so the two merge.
    """
    session, fp3, fp1 = p3_partitioned
    dual = _p3_graph(session, fp3, fp1)
    arrangement = dual.arrangement

    assert len(dual.face_nodes) == len(arrangement.faces) - 1, (
        "exactly one pair of faces merges on the nested fixture"
    )
    outer = [node for node in dual.face_nodes if node.kind == "outer"]
    assert len(outer) == 2
    assert sum(node.is_unbounded for node in dual.face_nodes) == 1

    unbounded = dual.unbounded
    assert len(unbounded.faces) == 1 and not unbounded.faces[0].is_closed

    (merged,) = [node for node in outer if node is not unbounded]
    assert len(merged.faces) == 2
    containing = [face for face in merged.faces if face.is_closed]
    opened = [face for face in merged.faces if not face.is_closed]
    assert len(containing) == 1 and len(opened) == 1
    assert containing[0] in arrangement.containing_faces
    assert opened[0] in arrangement.open_faces
    assert arrangement.component_of(
        next(corner for corner in containing[0].corners if corner >= 0)
    ) != arrangement.component_of(
        next(corner for corner in opened[0].corners if corner >= 0)
    )

    for region in arrangement.regions:
        assert dual.face_of(region).is_region


# --------------------------------------------------------------------------- #
# Handedness: the combinatorial face side against the geometry
# --------------------------------------------------------------------------- #
def test_k10_face_side_agrees_with_the_geometry(k10_partitioned):
    session, fp = k10_partitioned
    dual = _k10_graph(session, fp)

    checks, failures = _face_side_checks(dual)
    assert not failures, "\n".join(failures)
    assert checks == _checkable_pairs(dual), (
        f"only {checks} of the regions' stable arcs were checked"
    )


@pytest.mark.slow
def test_p3_face_side_agrees_with_the_geometry(p3_partitioned):
    session, fp3, fp1 = p3_partitioned
    dual = _p3_graph(session, fp3, fp1)

    checks, failures = _face_side_checks(dual)
    assert not failures, "\n".join(failures)
    assert checks == _checkable_pairs(dual), (
        f"only {checks} of the regions' stable arcs were checked"
    )


# --------------------------------------------------------------------------- #
# B.3 -- the fill
# --------------------------------------------------------------------------- #
def test_k10_fill_is_the_arcs_between_the_pip_and_its_image(k10_partitioned):
    session, fp = k10_partitioned
    dual = _k10_graph(session, fp)
    trellis = session.trellis(fp)
    branch_key = _pip_branch_key(trellis)

    filled = {node.key for node in dual.arc_nodes.values() if node.filled}
    assert filled == _expected_fill_ids(dual, trellis)
    assert filled, "the k=10 fixture must fill something"
    for node in dual.arc_nodes.values():
        if node.key not in filled:
            assert not node.filled
        assert node.faces.keys() == {"left", "right"}

    image = trellis.iterate(trellis.strong_pip, fp.k_value)
    low, high = dual.fill_segments[branch_key]
    tol = trellis.registry.cdist_tol
    assert abs(low - trellis.intersection(image).stable_cdist) <= tol
    assert abs(high - trellis.intersection(trellis.strong_pip).stable_cdist) <= tol


@pytest.mark.slow
def test_p3_fill_is_the_pip_segment_on_each_pips_own_branch(p3_partitioned):
    """Each tangle fills ``(f^k(q0), q0]`` on its pip's branch and nothing else.

    The nested fixture has TWO strong pips — one per tangle — so two branches
    carry a fill: the period-3 pip's branch (a three-step segment) and the
    period-1 pip's branch (a one-step one). Every OTHER branch of the period-3
    orbit is empty, because its last crossing is exactly the image of the last
    crossing of the branch before it.
    """
    session, fp3, fp1 = p3_partitioned
    dual = _p3_graph(session, fp3, fp1)
    t3, t1 = session.trellis(fp3), session.trellis(fp1)
    tol = dual.trellis.registry.cdist_tol

    expected: set[tuple] = set()
    pip_branches = set()
    for trellis, fixed_point in ((t3, fp3), (t1, fp1)):
        branch_key = _pip_branch_key(trellis)
        pip_branches.add(branch_key)
        expected |= _expected_fill_ids(dual, trellis)

        image = trellis.iterate(trellis.strong_pip, fixed_point.k_value)
        low, high = dual.fill_segments[branch_key]
        assert abs(low - dual.trellis.intersection(image).stable_cdist) <= tol
        assert (
            abs(high - dual.trellis.intersection(trellis.strong_pip).stable_cdist)
            <= tol
        )
        assert high > low + tol

    filled = {node.key for node in dual.arc_nodes.values() if node.filled}
    assert filled == expected
    assert filled

    for branch in dual.trellis.stable_branches:
        if branch.key in pip_branches:
            continue
        low, high = dual.fill_segments[branch.key]
        assert abs(high - low) <= tol, (
            f"branch {branch.key[1:]} carries no pip, so its fill must be empty"
        )
        assert not any(node.filled for node in dual.arc_nodes_on(branch.key))


# --------------------------------------------------------------------------- #
# C.3 -- element images
# --------------------------------------------------------------------------- #
def _assert_element_images(dual: DualGraph, per_fixed_point) -> None:
    """Every bounded element has an image, and it is the trellis's own answer."""
    images = dual.element_images(1)
    assert len(images) == sum(
        len(result.intervals) for result in dual.partitions.values()
    )

    for (branch_key, side), result in dual.partitions.items():
        trellis = per_fixed_point[branch_key[0]]
        for interval in result.intervals:
            ref = result.ref(interval.element_id)
            covering = images[ref]
            if interval.lo_id is None or interval.hi_id is None:
                assert covering == []
                continue
            assert covering, f"bounded element {ref} has no image"
            own = trellis.image_of_element(result, interval.element_id, 1)
            assert own is not None
            assert [item.element_id for item in covering] == own
            for image_ref in covering:
                assert isinstance(image_ref, ElementRef)
                assert image_ref.branch_key == branch_key[0].advance_key(branch_key, 1)
                assert image_ref.side == side  # the fixtures preserve orientation


def test_k10_element_images_match_the_trellis(k10_partitioned):
    session, fp = k10_partitioned
    dual = _k10_graph(session, fp)
    assert dual.trellis.orientation_preserving
    _assert_element_images(dual, {fp: session.trellis(fp)})


@pytest.mark.slow
def test_p3_element_images_match_the_trellis(p3_partitioned):
    session, fp3, fp1 = p3_partitioned
    dual = _p3_graph(session, fp3, fp1)
    assert dual.trellis.orientation_preserving
    _assert_element_images(
        dual, {fp3: session.trellis(fp3), fp1: session.trellis(fp1)}
    )


# --------------------------------------------------------------------------- #
# Views
# --------------------------------------------------------------------------- #
def test_k10_graph_is_bipartite(k10_partitioned):
    import networkx as nx

    session, fp = k10_partitioned
    dual = _k10_graph(session, fp)
    graph = dual.graph()

    assert graph.number_of_nodes() == len(dual.arc_nodes) + len(dual.face_nodes)
    assert nx.is_bipartite(graph)
    # One edge per (face, arc) pair, NOT per slot: a cut arc with the same
    # merged face node on both sides is one edge of a simple graph.
    pairs = {
        (face.index, node.key)
        for face in dual.face_nodes
        for node, _side in face.arcs
    }
    assert graph.number_of_edges() == len(pairs)
    for kind_a, kind_b in ((a[0], b[0]) for a, b in graph.edges):
        assert {kind_a, kind_b} == {"arc", "face"}
    for face in dual.face_nodes:
        for node, side in face.arcs:
            assert graph.edges[("face", face.index), ("arc", node.key)]["side"] == side


def test_k10_summary_reports_counts(k10_partitioned):
    session, fp = k10_partitioned
    dual = _k10_graph(session, fp)
    text = dual.summary()
    assert f"{len(dual.arc_nodes)} arc nodes" in text
    assert f"{len(dual.face_nodes)} face nodes" in text
    assert "fill segments" in text


def test_duplicate_partitions_are_rejected(k10_partitioned):
    session, fp = k10_partitioned
    partitions = list(session.trellis(fp).stable_partitions)

    with pytest.raises(ValueError, match="two partitions"):
        DualGraph(session.arrangement(), partitions + partitions[:1])


def test_a_missing_partition_side_is_rejected(k10_partitioned):
    session, fp = k10_partitioned
    partitions = [
        result
        for result in session.trellis(fp).stable_partitions
        if result.side == "left"
    ]

    with pytest.raises(ValueError, match="partition covers stable branch"):
        DualGraph(session.arrangement(), partitions)
