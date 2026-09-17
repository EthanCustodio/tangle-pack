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
    dg = _k10_dual_graph(k10_partitioned)
    fig, ax = plt.subplots()
    try:
        plotting.plot_dual_graph(dg, ax=ax, show_labels=True)
        texts = [t.get_text() for t in ax.texts]
        assert len(texts) == len(dg.stable_nodes) + len(dg.face_nodes)
        for node in dg.stable_nodes.values():
            assert "|".join(node.elements[s].label for s in node.sides) in texts
        for face in dg.face_nodes:
            assert f"{face.index}:{face.kind}" in texts
    finally:
        plt.close(fig)


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
