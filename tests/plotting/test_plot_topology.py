"""Plotting tier, property half: what the drawings say about the topology.

Each test checks a property a reader of the figure relies on -- a face point
lies inside its region, a side node sits on its side, a walk threads face ->
node -> face, the cartoon has one node per element and keeps the cdist order,
its brackets show closedness, every kept bridge is drawn, a nested circle
sits inside its parent, the zone interior is inside its circle, no colour is
used twice. No style, colour value, angle, radius, font, z-order or legend
wording is pinned (author decision 5, 2026-10-05).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

from cases import Case, build_k10, build_nested, build_period3  # noqa: E402
from tanglepack.topology import plotting  # noqa: E402
from tanglepack.topology.DualWalk import Walk  # noqa: E402
from tanglepack.topology.StablePartition import _side_of  # noqa: E402


@pytest.fixture
def fig_ax():
    """A fresh ``(figure, axes)`` pair, closed after the test."""
    fig, ax = plt.subplots()
    yield fig, ax
    plt.close(fig)


@pytest.fixture
def nested() -> Case:
    """The nested period-1 + period-3 tangle, outer zone blasted twice."""
    return build_nested()


def _symbol_names(dynamics, *, refined: bool) -> set[str]:
    """Every symbol a class colour map is keyed by: refined children or the class letter."""
    names = set()
    for bridge_class, cd in dynamics.classes.items():
        children = dynamics.refined.get(bridge_class, []) if refined else []
        names |= {c.name for c in children} if children else {cd.letter}
    return names


def _active_walk(dynamics):
    """The walk of the first class resolved by a dual-graph walk."""
    return next(cd.search.walk for cd in dynamics.classes.values()
                if cd.search is not None and cd.search.walk is not None)


# --------------------------------------------------------------------------- #
# Dual graph over the plane
# --------------------------------------------------------------------------- #
def test_face_point_of_a_region_is_inside_it(k10_partitioned) -> None:
    session, _fp = k10_partitioned
    dg = session.dual_graph()
    regions = [f for f in dg.face_nodes if f.is_region]
    assert regions
    for face in regions:
        assert face.faces[0].contains(plotting.face_point(dg, face))


def test_face_point_of_a_region_returns_a_copy_not_the_cached_array(k10_partitioned) -> None:
    """Regression: mutating a returned face point must not move the region's cached point."""
    session, _fp = k10_partitioned
    dg = session.dual_graph()
    face = next(f for f in dg.face_nodes if f.is_region)
    point = plotting.face_point(dg, face)
    point += 1e6
    assert not np.allclose(point, face.faces[0].representative_point)


def test_face_point_of_unbounded_node_is_outside_the_edge_bbox(k10_partitioned) -> None:
    session, _fp = k10_partitioned
    dg = session.dual_graph()
    mids = np.vstack([n.midpoint(dg.arrangement.trellis) for n in dg.stable_nodes.values()])
    point = plotting.face_point(dg, dg.unbounded)
    assert (point < mids.min(axis=0)).any() or (point > mids.max(axis=0)).any()


def test_side_nodes_sit_on_their_side_and_unified_nodes_on_the_edge(k10_partitioned) -> None:
    """A wall node is drawn on its own side of its stable edge, a unified one on it.

    The side is read against the edge's anchorward direction at its midpoint,
    from the public arc polyline, with the partition's own ``_side_of`` kernel.
    """
    session, _fp = k10_partitioned
    dg = session.dual_graph()
    trellis = dg.arrangement.trellis
    checked = 0
    for node in dg.stable_nodes.values():
        point = plotting.stable_node_point(dg, node)
        mid = node.midpoint(trellis)
        if point is None or mid is None:
            continue
        mid = np.asarray(mid, dtype=float)
        if node.is_unified:
            assert np.allclose(point, mid)
            continue
        polyline = np.asarray(node.arc.polyline(trellis), dtype=float)
        nearest = int(np.argmin(np.linalg.norm(polyline[:-1] - mid, axis=1)
                                + np.linalg.norm(polyline[1:] - mid, axis=1)))
        anchorward = polyline[nearest] - polyline[nearest + 1]
        assert _side_of(anchorward, point - mid) == node.sides[0]
        checked += 1
    assert checked == len(dg.wall_nodes) > 0


def test_plot_dual_graph_draws_every_stable_node_at_its_point(k10_partitioned, fig_ax) -> None:
    session, _fp = k10_partitioned
    dg = session.dual_graph()
    _fig, ax = fig_ax
    plotting.plot_dual_graph(dg, ax=ax)
    drawn = np.vstack([c.get_offsets() for c in ax.collections if len(c.get_offsets())])
    for node in dg.stable_nodes.values():
        point = plotting.stable_node_point(dg, node)
        if point is None:
            continue
        assert np.isclose(drawn, point).all(axis=1).any(), node.key


def test_plot_walk_threads_face_node_face(k10_partitioned, fig_ax) -> None:
    """A walk's line runs start face, then node and face per step, and ends in an arrowhead."""
    session, _fp = k10_partitioned
    dg = session.dual_graph()
    walk = _active_walk(session.symbolic_dynamics())
    assert walk.steps
    fig, ax = fig_ax
    line = plotting.plot_walk(dg, walk, ax=ax)
    xy = line.get_xydata()
    expected = [plotting.face_point(dg, walk.faces[0])]
    for step in walk.steps:
        expected += [plotting.stable_node_point(dg, step.node),
                     plotting.face_point(dg, step.to_face)]
    assert np.allclose(xy, np.vstack(expected))
    arrows = [t for t in ax.texts if t.arrow_patch is not None]
    assert len(arrows) == 1 and np.allclose(arrows[0].xy, xy[-1])
    fig.canvas.draw()


def test_plot_walk_of_a_trivial_or_empty_walk(k10_partitioned, fig_ax) -> None:
    """A one-face walk is a single point with no arrowhead; a faceless walk draws nothing."""
    session, _fp = k10_partitioned
    dg = session.dual_graph()
    face = dg.face_nodes[0]
    _fig, ax = fig_ax
    line = plotting.plot_walk(dg, Walk(steps=[], itinerary=(), faces=[face]), ax=ax)
    assert line.get_xydata().shape == (1, 2)
    assert not [t for t in ax.texts if t.arrow_patch is not None]
    assert plotting.plot_walk(dg, Walk(steps=[], itinerary=(), faces=[]), ax=ax) is None


def test_plot_minimal_trellis_draws_every_kept_bridge(k10_partitioned, fig_ax) -> None:
    session, _fp = k10_partitioned
    minimal = session.minimal_trellis()
    _fig, ax = fig_ax
    plotting.plot_minimal_trellis(minimal, ax=ax)
    assert len(ax.lines) == len(minimal.kept_bridge_ids) + len(minimal.dropped_bridge_ids)
    ax.cla()
    plotting.plot_minimal_trellis(minimal, ax=ax, show_dropped=False)
    assert len(ax.lines) == len(minimal.kept_bridge_ids)


# --------------------------------------------------------------------------- #
# Symbolic-dynamics views
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("refined", [True, False])
def test_transition_graph_nodes_are_the_active_symbols(k10_partitioned, fig_ax, refined) -> None:
    """The drawn graph is the dynamics' own: same nodes and edges, every node active."""
    session, _fp = k10_partitioned
    dynamics = session.symbolic_dynamics()
    expected = dynamics.transition_graph(refined=refined)
    _fig, ax = fig_ax
    graph = plotting.plot_transition_graph(dynamics, ax=ax, refined=refined)
    assert set(graph.nodes) == set(expected.nodes)
    assert set(graph.edges) == set(expected.edges)
    assert all(graph.nodes[n]["kind"] == "active" for n in graph.nodes)


def test_plot_bridges_by_class_draws_every_classed_bridge(k10_partitioned, fig_ax) -> None:
    session, _fp = k10_partitioned
    trellis, dynamics = session.trellis(), session.symbolic_dynamics()

    def drawable(cds) -> list:
        """Member bridge ids of ``cds`` with a bridge object to draw."""
        return [m.bridge_id for cd in cds for m in cd.entry.members
                if trellis.bridge_between(*m.bridge_id) is not None]

    members = drawable(dynamics.classes.values())
    inert = drawable(cd for cd in dynamics.classes.values() if cd.kind != "active")
    _fig, ax = fig_ax
    colors = plotting.plot_bridges_by_class(trellis, dynamics, ax=ax)
    assert len(ax.lines) == len(members)
    assert set(colors) == _symbol_names(dynamics, refined=True)
    ax.cla()
    plotting.plot_bridges_by_class(trellis, dynamics, ax=ax, show_inert=False)
    assert len(ax.lines) == len(members) - len(inert)
    ax.cla()
    colors = plotting.plot_bridges_by_class(trellis, dynamics, ax=ax, refined=False)
    assert set(colors) == _symbol_names(dynamics, refined=False)


def test_an_unmatched_member_is_still_drawn_under_its_own_colour(k10_partitioned, fig_ax) -> None:
    """A member matched to no refined child is drawn, under its parent letter, in a colour of its own."""
    session, _fp = k10_partitioned
    trellis, dynamics = session.trellis(), session.symbolic_dynamics()
    cd = next(c for c in dynamics.classes.values() if c.kind == "active")
    children = dynamics.refined[cd.bridge_class]
    # Pretend one matched member matched nothing.
    bid = children[0].members[0]
    children[0].members.remove(bid)
    dynamics.member_refinement[bid] = None
    dynamics.unmatched_members[bid] = children[0].occurrence

    members = [m.bridge_id for c in dynamics.classes.values() for m in c.entry.members
               if trellis.bridge_between(*m.bridge_id) is not None]
    _fig, ax = fig_ax
    colors = plotting.plot_bridges_by_class(trellis, dynamics, ax=ax)
    assert len(ax.lines) == len(members)
    assert cd.letter in colors
    assert len(set(colors.values())) == len(colors)
    assert len(plotting.bridge_class_legend_handles(dynamics)) == len(colors)


def test_bridges_by_class_can_be_restricted_to_the_minimal_trellis(nested: Case, fig_ax) -> None:
    """``bridge_ids`` draws exactly the classed bridges it names (the figures' trellis panel)."""
    session = nested.session
    dynamics = session.symbolic_dynamics()
    members = {m.bridge_id for cd in dynamics.classes.values() for m in cd.entry.members}
    kept = set(session.minimal_trellis().kept_bridge_ids)
    assert members - kept, "the blasts add class members outside the minimal trellis"
    trellis = session.trellis()
    _fig, ax = fig_ax
    plotting.plot_bridges_by_class(trellis, dynamics, ax=ax, bridge_ids=kept)
    expected = sum(1 for bid in members & kept if trellis.bridge_between(*bid) is not None)
    assert len(ax.lines) == expected


@pytest.mark.parametrize("builder", [build_k10, build_nested], ids=["k10", "nested"])
@pytest.mark.parametrize("refined", [True, False])
def test_class_colors_never_repeat(builder, refined) -> None:
    """One colour per symbol, never reused, on one tangle and across two."""
    dynamics = builder().session.symbolic_dynamics()
    palette = plotting.class_colors(dynamics, refined=refined)
    assert set(palette) == _symbol_names(dynamics, refined=refined)
    assert len(set(palette.values())) == len(palette)


# --------------------------------------------------------------------------- #
# The line cartoon
# --------------------------------------------------------------------------- #
def _cartoon(k10_partitioned, ax):
    """``(dual graph, dynamics, layout)`` of the k=10 line cartoon drawn on ``ax``."""
    session, _fp = k10_partitioned
    dual, dynamics = session.dual_graph(), session.symbolic_dynamics()
    return dual, dynamics, plotting.plot_dual_graph_cartoon(dual, dynamics, ax=ax)


def test_cartoon_layout_has_one_node_per_element_per_side(k10_partitioned, fig_ax) -> None:
    _fig, ax = fig_ax
    dual, _dynamics, layout = _cartoon(k10_partitioned, ax)
    elements = sum(len(result.intervals) for result in dual.partition)
    assert len(layout.nodes) == elements == len(layout.segments) == len(layout.unified)
    assert len(layout.rows) == len(dual.partition.branch_keys)
    for ref, (x, y) in layout.nodes.items():
        row = layout.rows[ref.branch_key]
        assert (y > row) == (ref.side == "left")
        x_lo, x_hi, _y_bar = layout.segments[ref]
        assert min(x_lo, x_hi) <= x <= max(x_lo, x_hi)
        assert x_lo >= x_hi, "anchor on the right: the anchorward end is rightmost"
    anchor_refs = [ref for ref, (x_lo, _x_hi, _y) in layout.segments.items()
                   if x_lo == layout.width]
    assert anchor_refs and all(dual.partition.element(ref).lo_id is None for ref in anchor_refs)
    unified = {ref for ref in layout.nodes if layout.unified[ref]}
    assert unified == {
        ref for ref in layout.nodes if any(node.is_unified for node in dual.stable_nodes_of(ref))
    }
    assert unified, "the fundamental segment unifies something"


def test_cartoon_ordinal_coordinate_is_strictly_increasing(k10_partitioned, fig_ax) -> None:
    _fig, ax = fig_ax
    dual, _dynamics, layout = _cartoon(k10_partitioned, ax)
    for branch_key, reps in layout.ranks.items():
        assert reps == sorted(reps) and len(set(reps)) == len(reps)
        assert all(b - a > layout.tol for a, b in zip(reps, reps[1:]))
        assert [layout.rank(branch_key, c) for c in reps] == list(range(len(reps)))
        boundaries = {
            round(c, 9)
            for side in dual.partition.sides(branch_key)
            for iv in dual.partition.result(branch_key, side).intervals
            for c in (iv.lo_cdist, iv.hi_cdist)
        }
        assert len(reps) == len(boundaries)
    for ref, (x_lo, x_hi, _y) in layout.segments.items():
        iv = dual.partition.element(ref)
        if ref in layout.singletons:
            assert x_lo == x_hi and iv.lo_cdist == iv.hi_cdist
        else:
            assert x_lo - x_hi >= 1
    left = plotting.dual_graph_cartoon_layout(dual, anchor="left")
    for ref, (x_lo, x_hi, _y) in left.segments.items():
        assert (x_lo, x_hi) == tuple(layout.width - x for x in layout.segments[ref][:2])


def test_cartoon_brackets_match_closedness(k10_partitioned, fig_ax) -> None:
    _fig, ax = fig_ax
    dual, _dynamics, layout = _cartoon(k10_partitioned, ax)
    glyphs = {}
    for text in ax.texts:
        if text.get_text() in ("[", "(", "]", ")"):
            x, y = text.xy if hasattr(text, "xy") else text.get_position()
            glyphs[(round(x, 6), round(y, 6), text.get_text() in "[(")] = text.get_text()
    for ref, (x_lo, x_hi, y_bar) in layout.segments.items():
        if ref in layout.singletons:
            continue
        iv = dual.partition.element(ref)
        # Anchor on the right: the anchorward (lo) end closes with ] / ).
        assert glyphs[(round(x_lo, 6), round(y_bar, 6), False)] == ("]" if iv.closed_lo else ")")
        assert glyphs[(round(x_hi, 6), round(y_bar, 6), True)] == ("[" if iv.closed_hi else "(")


def test_cartoon_draws_every_kept_bridge_and_walks_element_to_element(k10_partitioned, fig_ax) -> None:
    """Every minimal-trellis bridge is drawn; an itinerary runs node to node, U-arc then S-curve."""
    _fig, ax = fig_ax
    dual, dynamics, layout = _cartoon(k10_partitioned, ax)
    bridges = [p for p in ax.patches if (p.get_gid() or "").startswith("bridge:")]
    assert len(bridges) == len(dual.minimal.kept_bridge_ids) == layout.bridges_drawn
    active = next(cd for cd in dynamics.classes.values() if cd.kind == "active")
    walk = [p for p in ax.patches if p.get_gid() == f"walk:{active.letter}"]
    itinerary = active.itinerary
    assert len(walk) == len(itinerary) - 1
    for index, patch in enumerate(walk):
        vertices = patch.get_path().vertices
        a, b = itinerary[index], itinerary[index + 1]
        assert tuple(vertices[0]) == tuple(layout.nodes[a])
        assert tuple(vertices[-1]) == tuple(layout.nodes[b])
        assert (a.side == b.side) == (index % 2 == 0)
    # A trivial (inert) itinerary [X', X'] draws no walk.
    assert layout.walks_drawn == sum(
        1 for cd in dynamics.classes.values()
        if cd.itinerary is not None and len(set(cd.itinerary)) > 1
    )


# --------------------------------------------------------------------------- #
# The circle cartoon
# --------------------------------------------------------------------------- #
def _assert_interior_inside(layout) -> int:
    """Every node of a branch's interior side is inside its circle, the other side outside."""
    checked = 0
    for ref, (x, y) in layout.nodes.items():
        cx, cy, radius = layout.circles[layout.arcs[ref.branch_key][0]]
        inside = np.hypot(x - cx, y - cy) < radius
        assert inside == (ref.side != layout.outside[ref.branch_key]), ref
        checked += 1
    return checked


def test_circle_interior_side_lies_inside_the_circle(k10_partitioned, fig_ax) -> None:
    """The zone-interior side (right by default, or the one asked for) is drawn inside."""
    session, _fp = k10_partitioned
    dual, dynamics = session.dual_graph(), session.symbolic_dynamics()
    (branch_key,) = dual.partition.branch_keys
    _fig, ax = fig_ax
    layout = plotting.plot_dual_graph_cartoon(dual, dynamics, ax=ax, shape="circle")
    assert layout.outside[branch_key] == "left"
    assert _assert_interior_inside(layout) == len(layout.nodes) > 0
    assert layout.bridges_drawn == len(dual.minimal.kept_bridge_ids)
    flipped = plotting.dual_graph_zone_layout(dual, interior_side={branch_key: "left"})
    assert flipped.outside[branch_key] == "right"
    assert _assert_interior_inside(flipped) > 0


def test_p3_circle_follows_the_zone_boundary(fig_ax) -> None:
    """The arcs go round the circle in the zone's ring order and direction, interiors inside.

    Whatever ring order and direction the zone boundary has, the circle keeps
    that cyclic order in that direction, and every zone-interior side is drawn
    inside the circle.
    """
    case = build_period3()
    session = case.session
    (fp,) = case.fixed_points
    interior, parents, arc_order = session.cartoon_zones()
    assert parents == {}
    assert set(interior) == set(fp.branch_cycle("stable"))
    assert set(interior.values()) <= {"left", "right"}
    order = arc_order[id(fp)]
    assert len(order.order) == len(set(order.order)) == fp.k_value
    assert set(order.order) == set(fp.branch_cycle("stable"))

    _fig, ax = fig_ax
    layout = session.plot_dual_graph_cartoon(ax=ax, shape="circle")
    anchors = {key: layout.arcs[key][1] % (2 * np.pi) for key in order.order}
    around = sorted(order.order, key=lambda key: anchors[key], reverse=order.clockwise)
    start = around.index(order.order[0])
    assert around[start:] + around[:start] == list(order.order)
    assert _assert_interior_inside(layout) == len(layout.nodes) > 0


def test_nested_inner_circle_sits_inside_the_outer_one(nested: Case, fig_ax) -> None:
    """The nested tangle's circle lies inside its parent's keep-out disc; lobes enclosing it are known."""
    session = nested.session
    outer, inner = nested.fixed_points
    _interior, parents, _arc_order = session.cartoon_zones()
    assert parents == {id(inner): id(outer)}
    enclosures = session.cartoon_enclosures(session.minimal_trellis().kept_bridge_ids, parents)
    assert enclosures, "the outer anchor lobes enclose the period-3 tangle"
    assert all(ids == frozenset({id(inner)}) for ids in enclosures.values())

    _fig, ax = fig_ax
    layout = session.plot_dual_graph_cartoon(ax=ax, shape="circle")
    assert layout.parents == {id(inner): id(outer)}
    ox, oy, o_radius = layout.circles[id(outer)]
    ix, iy, i_radius = layout.circles[id(inner)]
    gap = np.hypot(ix - ox, iy - oy)
    assert gap + i_radius < o_radius
    assert gap + i_radius <= layout.keepout[id(outer)] < o_radius
    for ref, (x, y) in layout.nodes.items():
        if ref.branch_key[0] is inner:
            assert np.hypot(x - ox, y - oy) < layout.keepout[id(outer)]
