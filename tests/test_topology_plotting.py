"""The topology plotting split (plan Phase 8).

Pins that every trellis plotter lives once in ``tanglepack.topology.plotting``
(the Trellis methods delegating to it), that ``StablePartition`` is pure
partition logic with no drawing code left, that the ``verbose`` kwarg reports
through ``logging`` rather than ``print``, and that the session's plot fan-outs
all run through the single ``_fanout_plot`` helper.
"""

from __future__ import annotations

import logging

import matplotlib

matplotlib.use("Agg")  # headless: exercise the plot helpers without a display
import matplotlib.pyplot as plt
import numpy as np
import pytest

from tanglepack import TangleSession
from tanglepack.loom.TangleSession import TangleSession as _SessionClass
from tanglepack.topology import StablePartition, plotting
from minimal_helpers import build_pieces
from tanglepack.topology.DualGraph import DualGraph
from tanglepack.topology.Trellis import Trellis


@pytest.fixture
def henon_session(henon_map, henon_map_inverse):
    """A single-saddle k=10 session with bridges cut."""
    session = TangleSession(henon_map, henon_map_inverse)
    fp = session.construct_fixed_point([4, -4])
    session.orient_eigenvectors(
        fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
    )
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=9)
    session.grow_until_turnaround(fp, "stable")
    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    return session, fp


@pytest.fixture
def populated_trellis(henon_session):
    """A trellis carrying strong pip, pseudoneighbors, holes and partitions."""
    session, fp = henon_session
    trellis = session.trellis(fp)
    trellis.classify_strong_pips()
    trellis.compute_pseudoneighbors()
    trellis.punch_holes()
    trellis.partition_stable_manifold()
    return trellis


PLOTTERS = (
    "plot_strong_pip_candidates",
    "plot_strong_pip",
    "plot_pseudoneighbors",
    "plot_holes",
    "plot_stable_partition",
)


# --------------------------------------------------------------------------- #
# the split itself
# --------------------------------------------------------------------------- #
def test_every_plotter_is_a_module_function_in_plotting():
    """All five drawing routines live in topology/plotting.py."""
    for name in PLOTTERS:
        assert callable(getattr(plotting, name)), name


def test_stable_partition_module_holds_no_drawing_code():
    """StablePartition.py is partition logic only — no plotters, no pyplot."""
    assert not [n for n in vars(StablePartition) if n.startswith("plot_")]
    assert not hasattr(StablePartition, "plt")


def test_trellis_plotters_delegate_to_plotting(populated_trellis, monkeypatch):
    """Each Trellis.plot_* method calls its plotting.py counterpart."""
    for name in PLOTTERS:
        calls = []
        monkeypatch.setattr(
            plotting, name, lambda *a, **kw: calls.append((a, kw)) or "sentinel"
        )
        assert getattr(populated_trellis, name)() == "sentinel", name
        assert len(calls) == 1, name


def test_hole_style_conventions_live_once():
    """The per-orbit hole marker/colour tables are plotting.py's, not Trellis's."""
    assert len(plotting.HOLE_MARKERS) == len(plotting.HOLE_COLORS)
    assert not hasattr(Trellis, "_HOLE_MARKERS")
    assert not hasattr(Trellis, "_HOLE_COLORS")


def test_plotters_draw_on_a_real_tangle(populated_trellis):
    """Behaviour is unchanged: every plotter still returns its handle(s)."""
    plt.figure()
    try:
        ax = plt.gca()
        assert populated_trellis.plot_strong_pip_candidates(ax=ax) is not None
        assert populated_trellis.plot_strong_pip(ax=ax) is not None
        assert populated_trellis.plot_pseudoneighbors(ax=ax) is not None
        assert populated_trellis.plot_holes(ax=ax)
        assert populated_trellis.plot_stable_partition(ax=ax) is not None
    finally:
        plt.close()


# --------------------------------------------------------------------------- #
# verbose → logging
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "method, describe",
    [
        ("compute_pseudoneighbors", "describe_pseudoneighbors"),
        ("punch_holes", "describe_holes"),
        ("partition_stable_manifold", "describe_stable_partitions"),
    ],
)
def test_verbose_logs_and_never_prints(
    populated_trellis, capsys, caplog, method, describe
):
    """verbose=True emits the describe_* report at INFO, not on stdout."""
    capsys.readouterr()
    with caplog.at_level(logging.INFO, logger="tanglepack.topology.Trellis"):
        getattr(populated_trellis, method)(verbose=True)

    assert capsys.readouterr().out == ""
    report = getattr(populated_trellis, describe)()
    assert any(report in record.getMessage() for record in caplog.records)


def test_verbose_prints_when_logging_is_unconfigured(populated_trellis, capsys):
    """With no logging handler at all, verbose=True still reaches stdout.

    ``logging.lastResort`` only emits at WARNING, so an INFO record from an
    unconfigured application would vanish — the report must be printed instead.
    """
    root = logging.getLogger()
    saved_handlers, saved_level = root.handlers[:], root.level
    root.handlers = []
    try:
        assert not root.hasHandlers(), "the fixture must leave logging unconfigured"
        capsys.readouterr()
        populated_trellis.punch_holes(verbose=True)
        out = capsys.readouterr().out
    finally:
        root.handlers = saved_handlers
        root.setLevel(saved_level)

    assert populated_trellis.describe_holes() in out


def test_verbose_survives_an_application_level_above_info(populated_trellis, caplog):
    """verbose=True forces the report through a WARNING-level configuration."""
    logger = logging.getLogger("tanglepack.topology.Trellis")
    previous = logger.level
    logger.setLevel(logging.WARNING)
    try:
        with caplog.at_level(logging.INFO):
            populated_trellis.punch_holes(verbose=True)
        assert any(
            record.levelno == logging.INFO
            and record.name == "tanglepack.topology.Trellis"
            for record in caplog.records
        )
    finally:
        logger.setLevel(previous)
    assert logger.level == previous, "the level must be restored after the call"


# --------------------------------------------------------------------------- #
# session fan-outs
# --------------------------------------------------------------------------- #
def test_session_exposes_one_fanout_helper_per_shape():
    """TangleSession routes its fan-outs through two shared helpers."""
    assert callable(_SessionClass._fanout_plot)
    assert callable(_SessionClass._fanout_call)


@pytest.mark.parametrize(
    "name",
    [
        "plot_strong_pip_candidates",
        "plot_strong_pip",
        "plot_pseudoneighbors",
        "plot_holes",
        "plot_stable_partition",
    ],
)
def test_session_plot_fanouts_go_through_fanout_plot(henon_session, monkeypatch, name):
    """Every session plot_* fan-out is implemented by _fanout_plot."""
    session, _ = henon_session
    calls = []
    monkeypatch.setattr(
        _SessionClass,
        "_fanout_plot",
        lambda self, *a, **kw: calls.append((a, kw)) or ["sentinel"],
    )
    assert getattr(session, name)() == ["sentinel"], name
    assert len(calls) == 1, name


@pytest.mark.parametrize(
    "name",
    [
        "classify_strong_pips",
        "compute_pseudoneighbors",
        "punch_holes",
        "partition_stable_manifold",
    ],
)
def test_session_call_fanouts_go_through_fanout_call(henon_session, monkeypatch, name):
    """Every session non-plot fan-out is implemented by _fanout_call."""
    session, _ = henon_session
    calls = []
    monkeypatch.setattr(
        _SessionClass,
        "_fanout_call",
        lambda self, *a, **kw: calls.append((a, kw)) or "sentinel",
    )
    assert getattr(session, name)() == "sentinel", name
    assert len(calls) == 1, name


# --------------------------------------------------------------------------- #
# dual graph
# --------------------------------------------------------------------------- #
def _k10_dual_graph(k10_partitioned) -> DualGraph:
    """The dual graph of the k10_partitioned fixture, with its strong pip."""
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    return DualGraph(pieces.minimal, pieces.iterated, strong_pips=pieces.strong_pips)


def test_plot_dual_graph_scatters_exactly_the_stable_nodes(k10_partitioned):
    """Collections 0 and 1 are the open and solid stable nodes, always in that order."""
    dg = _k10_dual_graph(k10_partitioned)
    fig, ax = plt.subplots()
    try:
        plotting.plot_dual_graph(dg, ax=ax)
        open_pts = ax.collections[0].get_offsets()
        solid_pts = ax.collections[1].get_offsets()
        assert len(open_pts) == len(dg.wall_nodes)
        assert len(solid_pts) == len(dg.unified_nodes)
        assert len(open_pts) + len(solid_pts) == len(dg.stable_nodes)
        assert ax.collections[0].get_facecolors().size == 0
    finally:
        plt.close(fig)


def test_side_nodes_sit_on_their_side_and_solid_nodes_on_the_edge(k10_partitioned):
    from tanglepack.topology.StablePartition import _side_of
    from tanglepack.topology.plotting import _anchorward_look

    dg = _k10_dual_graph(k10_partitioned)
    trellis = dg.arrangement.trellis
    checked = 0
    for node in dg.stable_nodes.values():
        point = plotting.stable_node_point(dg, node)
        mid = node.midpoint(trellis)
        if point is None or mid is None:
            continue
        if node.is_unified:
            assert np.allclose(point, mid)
            continue
        look = _anchorward_look(node, trellis)
        assert _side_of(look, point - mid) == node.sides[0]
        checked += 1
    assert checked == len(dg.wall_nodes)


def test_face_point_of_unbounded_node_is_outside_the_edge_bbox(k10_partitioned):
    dg = _k10_dual_graph(k10_partitioned)
    mids = np.vstack([n.midpoint(dg.arrangement.trellis) for n in dg.stable_nodes.values()])
    point = plotting.face_point(dg, dg.unbounded)
    assert (point < mids.min(axis=0)).any() or (point > mids.max(axis=0)).any()


def test_face_point_of_a_region_is_inside_it(k10_partitioned):
    dg = _k10_dual_graph(k10_partitioned)
    regions = [f for f in dg.face_nodes if f.is_region]
    assert regions
    for face in regions:
        assert face.faces[0].contains(plotting.face_point(dg, face))


def test_face_point_of_a_region_returns_a_copy_not_the_cached_array(k10_partitioned):
    dg = _k10_dual_graph(k10_partitioned)
    face = next(f for f in dg.face_nodes if f.is_region)
    point = plotting.face_point(dg, face)
    point += 1e6
    assert not np.allclose(point, face.faces[0].representative_point)


def test_scatter_kwargs_reject_facecolors_and_c(k10_partitioned):
    dg = _k10_dual_graph(k10_partitioned)
    fig, ax = plt.subplots()
    try:
        with pytest.raises(ValueError, match="facecolors"):
            plotting.plot_dual_graph(dg, ax=ax, facecolors="red")
        with pytest.raises(ValueError, match="'c'"):
            plotting.plot_dual_graph(dg, ax=ax, c="red")
        plotting.plot_dual_graph(dg, ax=ax, color="tab:green", s=30)
    finally:
        plt.close(fig)


def _clip_bbox(dg: DualGraph):
    trellis = dg.arrangement.trellis
    points = [n.midpoint(trellis) for n in dg.stable_nodes.values()]
    points = [p for p in points if p is not None] + [plotting.face_point(dg, dg.unbounded)]
    stacked = np.vstack(points)
    mins, maxs = stacked.min(axis=0), stacked.max(axis=0)
    pad = plotting.DUAL_GRAPH_PUSH_FRACTION * float(np.linalg.norm(maxs - mins))
    return mins - pad, maxs + pad


def test_clip_to_arcs_true_keeps_the_axes_within_the_padded_edge_bbox(k10_partitioned):
    dg = _k10_dual_graph(k10_partitioned)
    lo, hi = _clip_bbox(dg)
    fig, ax = plt.subplots()
    try:
        plotting.plot_dual_graph(dg, ax=ax)
        assert ax.get_xlim() == pytest.approx((lo[0], hi[0]))
        assert ax.get_ylim() == pytest.approx((lo[1], hi[1]))
    finally:
        plt.close(fig)


def test_clip_to_arcs_false_leaves_the_axes_to_autoscale(k10_partitioned):
    dg = _k10_dual_graph(k10_partitioned)
    lo, hi = _clip_bbox(dg)
    fig, ax = plt.subplots()
    try:
        plotting.plot_dual_graph(dg, ax=ax, clip_to_arcs=False)
        assert ax.get_xlim() != pytest.approx((lo[0], hi[0])) or ax.get_ylim() != pytest.approx(
            (lo[1], hi[1])
        )
    finally:
        plt.close(fig)


def test_show_labels_annotates_every_stable_and_face_node(k10_partitioned):
    """The default label style writes each stable node's element names as mathtext."""
    dg = _k10_dual_graph(k10_partitioned)
    fig, ax = plt.subplots()
    try:
        plotting.plot_dual_graph(dg, ax=ax, show_labels=True)
        texts = [t.get_text() for t in ax.texts]
        assert len(texts) == len(dg.stable_nodes) + len(dg.face_nodes)
        for node in dg.stable_nodes.values():
            expected = " | ".join(node.names[s].mathtext for s in node.sides)
            assert expected in texts
            assert plotting.stable_node_label(node) == expected
            assert expected.startswith("$") and expected.endswith("$")
        for face in dg.face_nodes:
            assert f"{face.index}:{face.kind}" in texts
        fig.canvas.draw()  # the mathtext must parse
    finally:
        plt.close(fig)


def test_show_labels_ref_style_keeps_the_element_labels(k10_partitioned):
    dg = _k10_dual_graph(k10_partitioned)
    fig, ax = plt.subplots()
    try:
        plotting.plot_dual_graph(dg, ax=ax, show_labels=True, label_style="ref")
        texts = [t.get_text() for t in ax.texts]
        assert len(texts) == len(dg.stable_nodes) + len(dg.face_nodes)
        for node in dg.stable_nodes.values():
            expected = "|".join(node.elements[s].label for s in node.sides)
            assert expected in texts
            assert plotting.stable_node_label(node, "ref") == expected
        with pytest.raises(ValueError, match="label_style"):
            plotting.plot_dual_graph(dg, ax=ax, show_labels=True, label_style="bogus")
    finally:
        plt.close(fig)


def test_stable_node_label_falls_back_to_the_element_label_without_a_name(k10_partitioned):
    dg = _k10_dual_graph(k10_partitioned)
    node = dg.unified_nodes[0]
    saved = dict(node.names)
    try:
        del node.names["left"]
        text = plotting.stable_node_label(node)
        assert text == f"{node.elements['left'].label} | {saved['right'].mathtext}"
    finally:
        node.names.update(saved)


@pytest.mark.parametrize(
    "text, expected",
    [
        ("R_1^2", "$R_{1}^{2}$"),
        ("R_3", "$R_{3}$"),
        ("a_2^-1", "$a_{2}^{-1}$"),
        ("u^-1", "$u^{-1}$"),
        ("u", "$u$"),
        ("new1", "$new1$"),
        ("p3@1.0:R_2^1", "p3@1.0:$R_{2}^{1}$"),
        ("$a_{1}$", "$a_{1}$"),
        ("", ""),
    ],
)
def test_name_mathtext(text, expected):
    assert plotting.name_mathtext(text) == expected


def test_name_mathtext_agrees_with_element_name_mathtext(k10_partitioned):
    dg = _k10_dual_graph(k10_partitioned)
    for node in dg.stable_nodes.values():
        for name in node.names.values():
            assert plotting.name_mathtext(name.text) == name.mathtext


def test_dual_graph_legend_handles_match_the_plotter():
    labels = [h.get_label() for h in plotting.dual_graph_legend_handles()]
    assert labels == [
        "stable node (open, wall)",
        "stable node (solid, traversable)",
        "face node",
        "face-node edge",
    ]


def test_plot_minimal_trellis_draws_every_kept_bridge(k10_partitioned):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    minimal = pieces.minimal
    fig, ax = plt.subplots()
    try:
        assert plotting.plot_minimal_trellis(minimal, ax=ax) is ax
        assert len(ax.lines) == len(minimal.kept_bridge_ids) + len(minimal.dropped_bridge_ids)
        ax.cla()
        plotting.plot_minimal_trellis(minimal, ax=ax, show_dropped=False)
        assert len(ax.lines) == len(minimal.kept_bridge_ids)
    finally:
        plt.close(fig)
    labels = [h.get_label() for h in plotting.minimal_trellis_legend_handles()]
    assert labels == ["hole bridge", "image bridge", "dropped bridge"]


def test_plot_stable_partition_accepts_row_labels(k10_partitioned):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    results = pieces.homotopy.as_list() + pieces.iterated.as_list()
    labels = [f"{r.side} {i}" for i, r in enumerate(results)]
    fig, ax = plt.subplots()
    try:
        plotting.plot_stable_partition(results, ax=ax, labels=labels)
        assert [t.get_text() for t in ax.get_yticklabels()] == labels
        with pytest.raises(ValueError, match="row labels"):
            plotting.plot_stable_partition(results, ax=ax, labels=labels[:1])
    finally:
        plt.close(fig)


def test_plot_stable_partition_draws_element_labels_at_midpoints(k10_partitioned):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    results = pieces.iterated.as_list()
    naming = _k10_dual_graph(k10_partitioned).naming

    def label(result, interval):
        ref = pieces.iterated.ref(result.branch_key, result.side, interval.element_id)
        return naming.name(ref).mathtext

    fig, ax = plt.subplots()
    try:
        plotting.plot_stable_partition(results, ax=ax)
        without = len(ax.texts)
        ax.cla()
        plotting.plot_stable_partition(results, ax=ax, element_labels=label)
        expected = {
            label(result, interval) for result in results for interval in result.intervals
        }
        drawn = {t.get_text() for t in ax.texts}
        assert expected <= drawn
        assert len(ax.texts) == without + sum(len(r.intervals) for r in results)
        # Each label sits at its interval's midpoint, above the row.
        by_text = {}
        for text in ax.texts:
            by_text.setdefault(text.get_text(), []).append(text)
        for row, result in enumerate(results):
            for interval in result.intervals:
                mid = 0.5 * (interval.lo_cdist + interval.hi_cdist)
                hits = [
                    t for t in by_text[label(result, interval)]
                    if t.xy[1] == row and np.isclose(t.xy[0], mid)
                ]
                assert hits and all(t.xyann[1] > 0 for t in hits)
        ax.cla()
        plotting.plot_stable_partition(results, ax=ax, element_labels=lambda r, i: "")
        assert len(ax.texts) == without
        fig.canvas.draw()
    finally:
        plt.close(fig)


# --------------------------------------------------------------------------- #
# walks
# --------------------------------------------------------------------------- #
class _FakeStep:
    def __init__(self, node, from_face, to_face):
        self.node = node
        self.entry_side = node.side_of(from_face)
        self.exit_side = "right" if self.entry_side == "left" else "left"
        self.from_face = from_face
        self.to_face = to_face


class _FakeWalk:
    """The shape ``plot_walk`` reads: ``faces`` and ``steps`` with node/to_face."""

    def __init__(self, faces, steps):
        self.faces = faces
        self.steps = steps


def _one_crossing_walk(dg: DualGraph) -> _FakeWalk:
    """A synthetic one-step walk through the first unified node of a graph."""
    node = dg.unified_nodes[0]
    start, goal = node.face_on("left"), node.face_on("right")
    return _FakeWalk([start, goal], [_FakeStep(node, start, goal)])


def test_plot_walk_threads_face_node_face(k10_partitioned):
    dg = _k10_dual_graph(k10_partitioned)
    walk = _one_crossing_walk(dg)
    fig, ax = plt.subplots()
    try:
        line = plotting.plot_walk(dg, walk, ax=ax)
        xy = line.get_xydata()
        assert xy.shape == (3, 2)
        assert np.allclose(xy[0], plotting.face_point(dg, walk.faces[0]))
        assert np.allclose(xy[1], plotting.stable_node_point(dg, walk.steps[0].node))
        assert np.allclose(xy[2], plotting.face_point(dg, walk.steps[0].to_face))
        assert line.get_color() == plotting.WALK_STYLE["color"]
        assert line.get_zorder() == plotting.WALK_STYLE["zorder"]
        # A start dot and an end arrowhead accompany the line.
        markers = [l for l in ax.lines if l is not line]
        assert len(markers) == 1 and np.allclose(markers[0].get_xydata()[0], xy[0])
        arrows = [t for t in ax.texts if t.arrow_patch is not None]
        assert len(arrows) == 1 and np.allclose(arrows[0].xy, xy[-1])
        fig.canvas.draw()
    finally:
        plt.close(fig)


def test_plot_walk_honours_overrides_and_trivial_walks(k10_partitioned):
    dg = _k10_dual_graph(k10_partitioned)
    walk = _one_crossing_walk(dg)
    fig, ax = plt.subplots()
    try:
        line = plotting.plot_walk(dg, walk, ax=ax, color="tab:purple", linewidth=4, zorder=3)
        assert line.get_color() == "tab:purple" and line.get_linewidth() == 4
        assert line.get_zorder() == 3
        arrows = [t for t in ax.texts if t.arrow_patch is not None]
        assert arrows[0].arrow_patch.get_edgecolor()[:3] == pytest.approx(
            plt.matplotlib.colors.to_rgb("tab:purple")
        )
        ax.cla()
        trivial = _FakeWalk([walk.faces[0]], [])
        line = plotting.plot_walk(dg, trivial, ax=ax)
        assert line.get_xydata().shape == (1, 2)
        assert not [t for t in ax.texts if t.arrow_patch is not None]
        ax.cla()
        assert plotting.plot_walk(dg, _FakeWalk([], []), ax=ax) is None
    finally:
        plt.close(fig)


def test_walk_legend_handles_match_the_plotter():
    (handle,) = plotting.walk_legend_handles()
    assert handle.get_label() == "walk"
    assert handle.get_color() == plotting.WALK_STYLE["color"]
    (handle,) = plotting.walk_legend_handles(color="tab:purple")
    assert handle.get_color() == "tab:purple"


def test_walk_zorder_sits_between_dual_graph_edges_and_nodes():
    """The documented ladder: walks above edges and face points, below nodes."""
    assert plotting.DUAL_GRAPH_FACE_STYLE["zorder"] < plotting.WALK_STYLE["zorder"]
    assert plotting.MINIMAL_TRELLIS_HOLE_STYLE["zorder"] < plotting.WALK_STYLE["zorder"]
    assert plotting.CLASS_BRIDGE_STYLE["zorder"] < plotting.WALK_STYLE["zorder"]
    assert plotting.WALK_STYLE["zorder"] < plotting.DUAL_GRAPH_ARC_STYLE["zorder"]


# --------------------------------------------------------------------------- #
# symbolic dynamics
# --------------------------------------------------------------------------- #
def _k10_dynamics(k10_partitioned):
    """``(pieces, dual graph, symbolic dynamics)`` of the k=10 fixture."""
    module = pytest.importorskip("tanglepack.topology.SymbolicDynamics")
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    if any(entry.letter is None for entry in pieces.table.active):
        pieces.table = session.bridge_classes([fp])
    dual = DualGraph(pieces.minimal, pieces.iterated, strong_pips=pieces.strong_pips)
    return pieces, dual, module.symbolic_dynamics(dual, pieces.table)


def _symbol_names(dynamics, *, refined: bool) -> set[str]:
    names = set()
    for bridge_class, cd in dynamics.classes.items():
        children = dynamics.refined.get(bridge_class, []) if refined else []
        names |= {c.name for c in children} if children else {cd.letter}
    return names


def test_class_colors_are_assigned_in_fixed_order(k10_partitioned):
    _pieces, _dual, dynamics = _k10_dynamics(k10_partitioned)
    refined = plotting.class_colors(dynamics, refined=True)
    plain = plotting.class_colors(dynamics, refined=False)
    assert set(refined) == _symbol_names(dynamics, refined=True)
    assert set(plain) == _symbol_names(dynamics, refined=False)
    assert list(refined.values()) == list(plotting.CLASS_COLORS[: len(refined)])
    assert len(set(refined.values())) == len(refined)
    labels = [h.get_label() for h in plotting.bridge_class_legend_handles(dynamics)]
    assert labels == [plotting.name_mathtext(name) for name in refined]


def test_plot_bridges_by_class_draws_every_classed_bridge(k10_partitioned):
    pieces, _dual, dynamics = _k10_dynamics(k10_partitioned)
    members = [
        m.bridge_id for cd in dynamics.classes.values() for m in cd.entry.members
        if pieces.full.bridge_between(*m.bridge_id) is not None
    ]
    inert_members = [
        m.bridge_id for cd in dynamics.classes.values() if cd.kind != "active"
        for m in cd.entry.members if pieces.full.bridge_between(*m.bridge_id) is not None
    ]
    fig, ax = plt.subplots()
    try:
        colors = plotting.plot_bridges_by_class(pieces.full, dynamics, ax=ax)
        assert len(ax.lines) == len(members)
        palette = plotting.class_colors(dynamics, refined=True)
        for name, color in colors.items():
            assert color == palette.get(name, plotting.CLASS_UNKNOWN_COLOR)
        assert set(colors) & _symbol_names(dynamics, refined=True)
        # An unmatched member of a split class is drawn under its parent letter.
        unmatched = {
            dynamics.classes[pieces.table.entry_of(bid).bridge_class].letter
            for bid, child in dynamics.member_refinement.items()
            if child is None
        }
        assert set(colors) == _symbol_names(dynamics, refined=True) | unmatched
        dashed = [l for l in ax.lines if l.get_linestyle() == "--"]
        assert len(dashed) == len(inert_members)
        assert ax.get_legend() is not None
        legend_texts = [t.get_text() for t in ax.get_legend().get_texts()]
        assert len(legend_texts) == len(colors)
        assert all(t.endswith(" (unmatched)") for t in legend_texts if "(unmatched)" in t)
        fig.canvas.draw()

        ax.cla()
        plotting.plot_bridges_by_class(pieces.full, dynamics, ax=ax, show_inert=False)
        assert len(ax.lines) == len(members) - len(inert_members)
        ax.cla()
        colors = plotting.plot_bridges_by_class(pieces.full, dynamics, ax=ax, refined=False)
        assert set(colors) == _symbol_names(dynamics, refined=False)
        with pytest.raises(ValueError, match="manages"):
            plotting.plot_bridges_by_class(pieces.full, dynamics, ax=ax, color="red")
        with pytest.raises(ValueError, match="manages"):
            plotting.plot_bridges_by_class(pieces.full, dynamics, ax=ax, linestyle=":")
    finally:
        plt.close(fig)


@pytest.mark.parametrize("refined", [True, False])
def test_plot_transition_graph_draws_every_symbol(k10_partitioned, refined):
    _pieces, _dual, dynamics = _k10_dynamics(k10_partitioned)
    expected = dynamics.transition_graph(refined=refined)
    fig, ax = plt.subplots()
    try:
        graph = plotting.plot_transition_graph(dynamics, ax=ax, refined=refined)
        assert set(graph.nodes) == set(expected.nodes)
        assert set(graph.edges) == set(expected.edges)
        assert "spectral radius" in ax.get_title()
        assert f"{dynamics.spectral_radius(refined=refined):.4g}" in ax.get_title()
        assert not ax.axison
        # One label per node, as mathtext.
        labels = {t.get_text() for t in ax.texts if t.get_text().startswith("$")}
        assert labels == {plotting.name_mathtext(node) for node in graph.nodes}
        # Sinks (inert/virtual) are drawn dashed, sources solid.
        sinks = [n for n in graph.nodes if graph.nodes[n]["kind"] != "active"]
        assert sinks, "the k=10 fixture has an inert class"
        node_collections = [
            c for c in ax.collections
            if c.get_offsets().shape[0] in (len(sinks), len(graph) - len(sinks))
        ]
        dashed = [c for c in node_collections if c.get_linestyle()[0][1] is not None]
        assert any(c.get_offsets().shape[0] == len(sinks) for c in dashed)
        # A self-loop exists and is drawn (as a ring off its node).
        loops = [e for e in graph.edges if e[0] == e[1]]
        assert loops, "the k=10 word a -> a u^-1 a^-1 has a self-loop"
        fig.canvas.draw()
    finally:
        plt.close(fig)


def test_plot_itinerary_table_lists_every_class(k10_partitioned):
    _pieces, _dual, dynamics = _k10_dynamics(k10_partitioned)
    fig, ax = plt.subplots()
    try:
        table = plotting.plot_itinerary_table(dynamics, ax=ax)
        assert table is not None and not ax.axison
        cells = table.get_celld()
        n_rows = 1 + max(row for row, _ in cells)
        n_cols = 1 + max(col for _, col in cells)
        assert n_rows == 1 + len(dynamics.classes)
        assert n_cols == len(plotting.ITINERARY_TABLE_COLUMNS)
        headers = [cells[(0, c)].get_text().get_text() for c in range(n_cols)]
        assert headers == list(plotting.ITINERARY_TABLE_COLUMNS)
        shown = plotting.itinerary_table_rows(dynamics)
        for r, row in enumerate(shown, start=1):
            assert [cells[(r, c)].get_text().get_text() for c in range(n_cols)] == row
        # The cells compile: every name and symbol is a mathtext run.
        rows = plotting.itinerary_table_rows(dynamics, mathtext=False)
        for shown_row, row, cd in zip(shown, rows, dynamics.classes.values()):
            assert shown_row[0].startswith(f"${cd.letter}$ = {{$")
            assert "$" not in "".join(row)
            if cd.itinerary is not None:
                assert shown_row[1].count("$") == 2 * len(cd.itinerary)
                assert shown_row[1].count(" | ") == row[1].count(" | ")
            assert shown_row[2] == " ".join(s.mathtext for s in cd.symbols)
            assert shown_row[-2:] == row[-2:]
        # Row content: names, not ids; status words; the k=10 active class is walked.
        for row, cd in zip(rows, dynamics.classes.values()):
            assert row[0].startswith(f"{cd.letter} = {{")
            assert "#" not in row[0] and "#" not in row[1]
            if cd.itinerary is not None:
                assert row[1] == dynamics.itinerary_text(cd.itinerary)
                assert row[1].count(" | ") == len(cd.itinerary) // 2 - 1
            assert row[2] == cd.word
            assert row[-1] in ("yes", "no", "-")
        widths = [cells[(0, c)].get_width() for c in range(n_cols)]
        statuses = {row[-2] for row in rows}
        assert statuses <= {"unique", "trivial", "ambiguous", "trellis"} | {
            s for s in statuses if s.startswith("unresolved: ")
        }
        active = next(cd for cd in dynamics.classes.values() if cd.kind == "active")
        active_row = rows[list(dynamics.classes.values()).index(active)]
        assert active_row[2] and active_row[3]  # word and refined word
        fig.canvas.draw()

        # The columns fill the axes exactly, widest text widest.
        assert all(w > 0.0 for w in widths) and abs(sum(widths) - 1.0) < 1e-9
        assert widths[1] == max(widths), "the itinerary column is the widest"

        ax.cla()
        plain_table = plotting.plot_itinerary_table(dynamics, ax=ax, mathtext=False)
        plain_cells = plain_table.get_celld()
        for r, row in enumerate(rows, start=1):
            assert [plain_cells[(r, c)].get_text().get_text() for c in range(n_cols)] == row
        # Widths are measured on the texts shown, so the (wider) plain texts
        # give the mathtext-heavy columns more room.
        plain_widths = [plain_cells[(0, c)].get_width() for c in range(n_cols)]
        assert abs(sum(plain_widths) - 1.0) < 1e-9 and plain_widths[1] > widths[1]

        ax.cla()
        table = plotting.plot_itinerary_table(dynamics, ax=ax, refined=False)
        n_cols = 1 + max(col for _, col in table.get_celld())
        assert n_cols == len(plotting.ITINERARY_TABLE_COLUMNS) - 1
    finally:
        plt.close(fig)


def test_plot_itinerary_table_columns_font_and_rows(k10_partitioned):
    """``columns`` picks a subset in the caller's order; the font grows to
    its cap on a wide axes and the rows grow with it; a narrow axes shrinks it."""
    _pieces, _dual, dynamics = _k10_dynamics(k10_partitioned)
    columns = ("word", "class", "refined word")
    rows = plotting.itinerary_table_rows(dynamics, mathtext=False, columns=columns)
    full = plotting.itinerary_table_rows(dynamics, mathtext=False)
    for row, full_row in zip(rows, full):
        assert row == [full_row[2], full_row[0], full_row[3]]
    with pytest.raises(ValueError, match="unknown itinerary table column"):
        plotting.itinerary_table_rows(dynamics, columns=("class", "bogus"))

    fig, ax = plt.subplots(figsize=(16, 6))
    try:
        table = plotting.plot_itinerary_table(
            dynamics, ax=ax, columns=columns, fontsize=20.0, row_height=3.0
        )
        cells = table.get_celld()
        n_cols = 1 + max(col for _, col in cells)
        assert n_cols == 3
        assert [cells[(0, c)].get_text().get_text() for c in range(n_cols)] == list(columns)
        assert cells[(1, 0)].get_fontsize() == 20.0
        # Rows are ``row_height`` ems of the font: 60 pt of a 6 in (432 pt) axes.
        axes_points = ax.get_window_extent().height / fig.dpi * 72.0
        assert cells[(1, 0)].get_height() * axes_points == pytest.approx(60.0, rel=0.05)
        with pytest.raises(ValueError):
            plotting.plot_itinerary_table(dynamics, ax=ax, columns=("class", "bogus"))
    finally:
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(2.5, 3))
    try:
        table = plotting.plot_itinerary_table(dynamics, ax=ax, fontsize=20.0)
        size = table.get_celld()[(1, 0)].get_fontsize()
        assert 4.0 <= size < 20.0, "the widest row does not fit at 20 pt on a 2.5 in axes"
    finally:
        plt.close(fig)


def test_unmatched_member_gets_its_parent_letter_slot(k10_partitioned):
    pieces, _dual, dynamics = _k10_dynamics(k10_partitioned)
    cd = next(c for c in dynamics.classes.values() if c.kind == "active")
    letter = cd.letter
    children = dynamics.refined[cd.bridge_class]
    assert dynamics.unmatched_members == {}, "the footprint rule matches every member"
    # Pretend one matched member matched nothing.
    bid = children[0].members[0]
    children[0].members.remove(bid)
    dynamics.member_refinement[bid] = None
    dynamics.unmatched_members[bid] = children[0].occurrence

    refined = plotting.class_colors(dynamics, refined=True)
    names = list(refined)
    assert names.index(letter) == names.index(f"{letter}_2") + 1
    assert refined[letter] == plotting.CLASS_COLORS[names.index(letter)]
    assert list(plotting.class_colors(dynamics, refined=False)) == [
        c.letter for c in dynamics.classes.values()
    ]
    labels = [h.get_label() for h in plotting.bridge_class_legend_handles(dynamics)]
    assert f"${letter}$ (unmatched)" in labels

    fig, ax = plt.subplots()
    try:
        colors = plotting.plot_bridges_by_class(pieces.full, dynamics, ax=ax)
        assert colors[letter] == refined[letter]
        bridge = pieces.full.bridge_between(*bid)
        points = bridge.get_point_array()
        drawn = [
            line for line in ax.lines
            if len(line.get_xdata()) == len(points) and line.get_color() == refined[letter]
        ]
        assert drawn, "the unmatched member is drawn in its parent letter's colour"
        assert plotting.CLASS_UNKNOWN_COLOR not in colors.values()
    finally:
        plt.close(fig)


# --------------------------------------------------------------------------- #
# dual-graph cartoon
# --------------------------------------------------------------------------- #
def _cartoon(k10_partitioned):
    """``(pieces, dual, dynamics, fig, ax, layout)`` of the k=10 cartoon."""
    pieces, dual, dynamics = _k10_dynamics(k10_partitioned)
    fig, ax = plt.subplots()
    layout = plotting.plot_dual_graph_cartoon(dual, dynamics, ax=ax)
    return pieces, dual, dynamics, fig, ax, layout


def test_cartoon_layout_has_one_node_per_element_per_side(k10_partitioned):
    _pieces, dual, _dynamics, fig, ax, layout = _cartoon(k10_partitioned)
    try:
        elements = sum(len(result.intervals) for result in dual.partition)
        assert len(layout.nodes) == elements == len(layout.segments) == len(layout.unified)
        assert len(layout.rows) == len(dual.partition.branch_keys)
        assert layout.anchor == "right"
        for ref, (x, y) in layout.nodes.items():
            row = layout.rows[ref.branch_key]
            assert (y > row) == (ref.side == "left")
            x_lo, x_hi, _y_bar = layout.segments[ref]
            assert min(x_lo, x_hi) <= x <= max(x_lo, x_hi)
            # Mirrored: the anchorward end is to the right of the outward end.
            assert x_lo >= x_hi
        anchor_refs = [ref for ref, (x_lo, _x_hi, _y) in layout.segments.items() if x_lo == layout.width]
        assert anchor_refs and all(dual.partition.element(ref).lo_id is None for ref in anchor_refs)
        # One circle per non-singleton element; singletons are dots on the bar.
        circles = [c for c in ax.collections if len(c.get_offsets()) > 0]
        assert len(circles) == 1
        assert len(circles[0].get_offsets()) == elements - len(layout.singletons)
        # Filled iff the element carries a unified stable node.
        unified = {ref for ref in layout.nodes if layout.unified[ref]}
        assert unified == {
            ref for ref in layout.nodes
            if any(node.is_unified for node in dual.stable_nodes_of(ref))
        }
        assert unified, "the k=10 fundamental segment unifies something"
    finally:
        plt.close(fig)


def test_cartoon_ordinal_coordinate_is_strictly_increasing(k10_partitioned):
    _pieces, dual, _dynamics, fig, _ax, layout = _cartoon(k10_partitioned)
    try:
        for branch_key, reps in layout.ranks.items():
            assert reps == sorted(reps) and len(set(reps)) == len(reps)
            assert all(b - a > layout.tol for a, b in zip(reps, reps[1:]))
            assert [layout.rank(branch_key, c) for c in reps] == list(range(len(reps)))
            assert [layout.position(branch_key, c) for c in reps] == [
                layout.width - r for r in range(len(reps))
            ]
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
        assert layout.width == max(len(r) - 1 for r in layout.ranks.values())
        left = plotting.dual_graph_cartoon_layout(dual, anchor="left")
        for ref, (x_lo, x_hi, _y) in left.segments.items():
            assert (x_lo, x_hi) == tuple(layout.width - x for x in layout.segments[ref][:2])
        with pytest.raises(ValueError):
            plotting.dual_graph_cartoon_layout(dual, anchor="top")
    finally:
        plt.close(fig)


def test_cartoon_brackets_match_closedness(k10_partitioned):
    _pieces, dual, _dynamics, fig, ax, layout = _cartoon(k10_partitioned)
    try:
        glyphs = {}
        for text in ax.texts:
            if text.get_text() in ("[", "(", "]", ")"):
                x, y = text.xy if hasattr(text, "xy") else text.get_position()
                glyphs[(round(x, 6), round(y, 6), text.get_text() in "[(")] = text.get_text()
        for ref, (x_lo, x_hi, y_bar) in layout.segments.items():
            iv = dual.partition.element(ref)
            if ref in layout.singletons:
                continue
            # Anchor on the right: the anchorward (lo) end closes with ] / ).
            assert glyphs[(round(x_lo, 6), round(y_bar, 6), False)] == ("]" if iv.closed_lo else ")")
            assert glyphs[(round(x_hi, 6), round(y_bar, 6), True)] == ("[" if iv.closed_hi else "(")
        # The k=10 right row reads [ ] ( ) [ ] ( ) [ ] ( ) [ ] under the footprint rule.
        (branch_key,) = layout.rows
        right = dual.partition.result(branch_key, "right").intervals
        assert [(iv.closed_lo, iv.closed_hi) for iv in right] == [
            (True, True), (False, False), (True, True), (False, False),
            (True, True), (False, False), (True, True),
        ]
    finally:
        plt.close(fig)


def test_cartoon_draws_every_kept_bridge_and_walks_element_to_element(k10_partitioned):
    _pieces, dual, dynamics, fig, ax, layout = _cartoon(k10_partitioned)
    try:
        gids = [patch.get_gid() for patch in ax.patches]
        bridge_gids = [g for g in gids if g and g.startswith("bridge:")]
        assert len(bridge_gids) == len(dual.minimal.kept_bridge_ids) == layout.bridges_drawn
        # Bridges and crossings are grey; only the walks carry class colours.
        from matplotlib.colors import to_hex
        for patch in ax.patches:
            if patch.get_gid() in ("unified",) or patch.get_gid().startswith("bridge:"):
                assert to_hex(patch.get_edgecolor()) in (
                    plotting.CLASS_UNKNOWN_COLOR, plotting.CARTOON_UNIFIED_STYLE["color"]
                )
        assert set(bridge_gids) == {f"bridge:{a}-{b}" for a, b in dual.minimal.kept_bridge_ids}
        assert gids.count("unified") == len({
            (n.elements["left"], n.elements["right"]) for n in dual.unified_nodes
        })
        active = next(cd for cd in dynamics.classes.values() if cd.kind == "active")
        walk = [p for p in ax.patches if p.get_gid() == f"walk:{active.letter}"]
        itinerary = active.itinerary
        assert len(walk) == len(itinerary) - 1 == 2 * (len(itinerary) // 2) - 1
        first = tuple(walk[0].get_path().vertices[0])
        last = tuple(walk[-1].get_path().vertices[-1])
        assert first == tuple(layout.nodes[itinerary[0]])
        assert last == tuple(layout.nodes[itinerary[-1]])
        # Same-side pairs are U-arcs leaving the line, crossings are S-curves.
        for index, patch in enumerate(walk):
            vertices = patch.get_path().vertices
            a, b = itinerary[index], itinerary[index + 1]
            assert tuple(vertices[0]) == tuple(layout.nodes[a])
            assert tuple(vertices[-1]) == tuple(layout.nodes[b])
            assert (a.side == b.side) == (index % 2 == 0)
        # The trivial inert itinerary [X', X'] draws no walk.
        inert = [cd for cd in dynamics.classes.values() if cd.kind != "active"]
        assert inert and all(
            not [p for p in ax.patches if p.get_gid() == f"walk:{cd.letter}"] for cd in inert
        )
        assert layout.walks_drawn == 1
        assert not ax.axison
    finally:
        plt.close(fig)


def test_cartoon_legend_handles_match_the_plotter(k10_partitioned):
    _pieces, _dual, dynamics, fig, _ax, _layout = _cartoon(k10_partitioned)
    try:
        fixed = [h.get_label() for h in plotting.dual_graph_cartoon_legend_handles()]
        assert fixed == [
            "hole bridge", "image bridge", "inert class", "unified crossing",
            "open node (wall)", "unified node", "trellis itinerary (no walk)",
        ]
        handles = plotting.dual_graph_cartoon_legend_handles(dynamics)
        labels = [h.get_label() for h in handles]
        active = next(cd for cd in dynamics.classes.values() if cd.kind == "active")
        assert labels == fixed + [f"walk of ${active.letter}$"]
        assert handles[-1].get_color() == plotting.class_colors(dynamics, refined=False)[active.letter]
        assert handles[-1].get_linewidth() == plotting.CARTOON_WALK_LINEWIDTHS[0]
        open_node, unified_node = handles[4], handles[5]
        assert open_node.get_markerfacecolor() == "white"
        assert unified_node.get_markerfacecolor() == "black"
    finally:
        plt.close(fig)


def test_cartoon_toggles_and_session_delegate(k10_partitioned):
    session, fp = k10_partitioned
    fig, ax = plt.subplots()
    try:
        layout = session.plot_dual_graph_cartoon(
            ax=ax, fixed_points=[fp], show_bridges=False, show_walks=False, show_unified=False,
        )
        assert isinstance(layout, plotting.CartoonLayout)
        assert layout.bridges_drawn == 0 and layout.walks_drawn == 0
        assert not ax.patches
    finally:
        plt.close(fig)


def test_cartoon_bridge_colour_labels_and_walkless_legend(k10_partitioned):
    """``bridge_color`` recolours every bridge (crossings stay grey), the names
    sit to the right of their nodes at ``label_fontsize`` (or outside the
    node on request) and the legend drops the walk entries with ``walks=False``."""
    from matplotlib.colors import to_hex

    _pieces, dual, dynamics, fig, ax, layout = _cartoon(k10_partitioned)
    plt.close(fig)
    purple = to_hex("tab:purple")
    fig, ax = plt.subplots()
    try:
        layout = plotting.plot_dual_graph_cartoon(
            dual, dynamics, ax=ax, show_walks=False, bridge_color="tab:purple",
            bridge_alpha=0.7, label_fontsize=22.0,
        )
        assert layout.walks_drawn == 0 and layout.bridges_drawn > 0
        bridges = [p for p in ax.patches if p.get_gid().startswith("bridge:")]
        assert bridges and all(to_hex(p.get_edgecolor()) == purple for p in bridges)
        assert all(p.get_alpha() == 0.7 for p in bridges)
        crossings = [p for p in ax.patches if p.get_gid() == "unified"]
        assert crossings and all(
            to_hex(p.get_edgecolor()) == plotting.CARTOON_UNIFIED_STYLE["color"] for p in crossings
        )
        assert not [p for p in ax.patches if p.get_gid().startswith("walk:")]
        # One name per element, anchored at its node, offset to the right, 22 pt.
        names = {
            plotting._cartoon_name(dual, ref): ref for ref in layout.nodes
        }
        labels = [t for t in ax.texts if t.get_text() in names]
        assert len(labels) == len(layout.nodes)
        for text in labels:
            ref = names[text.get_text()]
            assert tuple(text.xy) == tuple(layout.nodes[ref])
            assert text.xyann[0] > 0.0 and text.xyann[1] == 0.0
            assert text.get_ha() == "left" and text.get_va() == "center"
            assert text.get_fontsize() == 22.0

        ax.cla()
        plotting.plot_dual_graph_cartoon(dual, dynamics, ax=ax, label_position="outside")
        outside = [t for t in ax.texts if t.get_text() in names]
        assert len(outside) == len(layout.nodes)
        assert all(t.get_ha() == "center" and t.get_fontsize() == plotting.CARTOON_LABEL_FONTSIZE
                   for t in outside)
        with pytest.raises(ValueError, match="label_position"):
            plotting.plot_dual_graph_cartoon(dual, dynamics, ax=ax, label_position="above")

        handles = plotting.dual_graph_cartoon_legend_handles(
            dynamics, bridge_color="tab:purple", bridge_alpha=0.7, walks=False
        )
        assert [h.get_label() for h in handles] == [
            "hole bridge", "image bridge", "inert class", "unified crossing",
            "open node (wall)", "unified node",
        ]
        assert all(to_hex(h.get_color()) == purple and h.get_alpha() == 0.7 for h in handles[:3])

        # A uniform width: every bridge 1.2 wide, one "bridge" legend entry.
        ax.cla()
        layout = plotting.plot_dual_graph_cartoon(
            dual, dynamics, ax=ax, show_walks=False, bridge_linewidth=1.2
        )
        bridges = [p for p in ax.patches if p.get_gid().startswith("bridge:")]
        assert len(bridges) == layout.bridges_drawn
        assert all(p.get_linewidth() == 1.2 for p in bridges)
        uniform = plotting.dual_graph_cartoon_legend_handles(bridge_linewidth=1.2, walks=False)
        assert [h.get_label() for h in uniform] == [
            "bridge", "inert class", "unified crossing", "open node (wall)", "unified node",
        ]
        assert uniform[0].get_linewidth() == uniform[1].get_linewidth() == 1.2
    finally:
        plt.close(fig)
