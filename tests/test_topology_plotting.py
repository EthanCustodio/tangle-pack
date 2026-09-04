"""The topology plotting split (plan Phase 8).

Pins that every trellis plotter lives once in ``tanglepack.topology.plotting``
(the Trellis methods delegating to it), that ``StablePartition`` is pure
partition logic with no drawing code left, that the ``verbose`` kwarg reports
through ``logging`` rather than ``print``, and that the session's plot fan-outs
all run through the single ``_fanout_plot`` helper.
"""

from __future__ import annotations

import logging
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")  # headless: exercise the plot helpers without a display
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pytest

from tanglepack import TangleSession
from tanglepack.loom.TangleSession import TangleSession as _SessionClass
from tanglepack.topology import StablePartition, plotting
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
# dual graph (B.4)
# --------------------------------------------------------------------------- #
def _k10_dual_graph(k10_partitioned) -> DualGraph:
    """The dual graph of the k10_partitioned fixture, partitions only."""
    session, fp = k10_partitioned
    return DualGraph(session.arrangement(), session.trellis(fp).stable_partitions)


def test_plot_dual_graph_scatters_exactly_the_arc_nodes(k10_partitioned):
    """The two arc-node collections split hollow/filled and sum to arc_nodes.

    ``plot_dual_graph`` always scatters the hollow set first, the filled set
    second (see its Dev Notes), so ``ax.collections[0]``/``[1]`` are pinned to
    those two regardless of whether either set is empty.
    """
    dg = _k10_dual_graph(k10_partitioned)
    filled_expected = sum(1 for node in dg.arc_nodes.values() if node.filled)
    hollow_expected = len(dg.arc_nodes) - filled_expected

    plt.figure()
    try:
        ax = plt.gca()
        result = plotting.plot_dual_graph(dg, ax=ax)
        assert result is ax
        hollow_collection, filled_collection = ax.collections[0], ax.collections[1]
        hollow_count = len(hollow_collection.get_offsets())
        filled_count = len(filled_collection.get_offsets())
        assert hollow_count + filled_count == len(dg.arc_nodes)
        assert filled_count == filled_expected
        assert hollow_count == hollow_expected
    finally:
        plt.close()


def test_face_point_of_unbounded_node_is_outside_the_arc_bbox(k10_partitioned):
    """The unbounded node's face point is pushed clear of every arc midpoint."""
    dg = _k10_dual_graph(k10_partitioned)
    mids = [
        node.midpoint(dg.trellis) for node in dg.arc_nodes.values()
    ]
    mids = np.vstack([m for m in mids if m is not None])
    mins, maxs = mids.min(axis=0), mids.max(axis=0)

    point = plotting.face_point(dg, dg.unbounded)
    assert np.any(point < mins) or np.any(point > maxs)


def test_face_point_of_a_region_is_inside_it(k10_partitioned):
    """A region node's face point lies inside its own region's boundary."""
    dg = _k10_dual_graph(k10_partitioned)
    region_nodes = [node for node in dg.face_nodes if node.kind == "region"]
    assert region_nodes, "the k=10 fixture must have at least one region node"
    for node in region_nodes:
        point = plotting.face_point(dg, node)
        assert node.faces[0].contains(point)


def test_face_point_of_a_region_returns_a_copy_not_the_cached_array(k10_partitioned):
    """Mutating the returned point must not corrupt the region's own cache."""
    dg = _k10_dual_graph(k10_partitioned)
    region_node = next(node for node in dg.face_nodes if node.kind == "region")
    region = region_node.faces[0]
    cached_before = region.representative_point.copy()

    point = plotting.face_point(dg, region_node)
    point += 1000.0

    assert np.allclose(region.representative_point, cached_before)


def test_scatter_kwargs_reject_facecolors_and_c(k10_partitioned):
    """facecolors/c are internally managed and raise a clear ValueError."""
    dg = _k10_dual_graph(k10_partitioned)
    plt.figure()
    try:
        ax = plt.gca()
        with pytest.raises(ValueError):
            plotting.plot_dual_graph(dg, ax=ax, facecolors="red")
        with pytest.raises(ValueError):
            plotting.plot_dual_graph(dg, ax=ax, c="blue")
    finally:
        plt.close()


def _clip_bbox(dg: DualGraph):
    """The (mins, maxs, pad) triple plot_dual_graph's clip_to_arcs computes."""
    mids = [node.midpoint(dg.trellis) for node in dg.arc_nodes.values()]
    points = [m for m in mids if m is not None]
    points.append(plotting.face_point(dg, dg.unbounded))
    stacked = np.vstack(points)
    mins, maxs = stacked.min(axis=0), stacked.max(axis=0)
    pad = plotting.DUAL_GRAPH_PUSH_FRACTION * np.linalg.norm(maxs - mins)
    return mins, maxs, pad


def test_clip_to_arcs_true_keeps_the_axes_within_the_padded_arc_bbox(k10_partitioned):
    """clip_to_arcs=True limits the view to the padded arc/unbounded-point bbox."""
    dg = _k10_dual_graph(k10_partitioned)
    mins, maxs, pad = _clip_bbox(dg)

    plt.figure()
    try:
        ax = plt.gca()
        plotting.plot_dual_graph(dg, ax=ax, clip_to_arcs=True)
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        assert xlim[0] == pytest.approx(mins[0] - pad)
        assert xlim[1] == pytest.approx(maxs[0] + pad)
        assert ylim[0] == pytest.approx(mins[1] - pad)
        assert ylim[1] == pytest.approx(maxs[1] + pad)
    finally:
        plt.close()


def test_clip_to_arcs_false_lets_an_outlier_region_blow_out_the_axes(k10_partitioned):
    """clip_to_arcs=False leaves the axes autoscaled to everything drawn.

    On k=10 at least one bounded region's representative point sits well
    outside the tangle's own arc-midpoint span, so the unclipped view must
    reach past the padded bbox clip_to_arcs=True would have used.
    """
    dg = _k10_dual_graph(k10_partitioned)
    mins, maxs, pad = _clip_bbox(dg)

    plt.figure()
    try:
        ax = plt.gca()
        plotting.plot_dual_graph(dg, ax=ax, clip_to_arcs=False)
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        outside = (
            xlim[0] < mins[0] - pad or xlim[1] > maxs[0] + pad
            or ylim[0] < mins[1] - pad or ylim[1] > maxs[1] + pad
        )
        assert outside
    finally:
        plt.close()


def test_show_labels_annotates_every_arc_and_face_node(k10_partitioned):
    """show_labels=True adds one Text per plotted arc node and one per face node."""
    dg = _k10_dual_graph(k10_partitioned)
    expected_arc_labels = sum(
        1 for node in dg.arc_nodes.values() if node.midpoint(dg.trellis) is not None
    )

    plt.figure()
    try:
        ax = plt.gca()
        result = plotting.plot_dual_graph(dg, ax=ax, show_labels=True)
        assert result is ax
        assert len(ax.texts) == expected_arc_labels + len(dg.face_nodes)
    finally:
        plt.close()


# --------------------------------------------------------------------------- #
# transition graph (D.6)
# --------------------------------------------------------------------------- #
def _weighted_symbolic_dynamics() -> SimpleNamespace:
    """A SymbolicDynamics-like stand-in: only ``.transition_graph`` is read.

    One weight-2 edge (a -> b) and one weight-1 edge (b -> c), so the two
    cases plot_transition_graph's edge-label rule distinguishes both appear.
    """
    graph = nx.DiGraph()
    graph.add_node("a")
    graph.add_node("b")
    graph.add_node("c")
    graph.add_edge("a", "b", weight=2)
    graph.add_edge("b", "c", weight=1)
    return SimpleNamespace(transition_graph=graph)


def test_plot_transition_graph_labels_only_edges_with_weight_above_one():
    """A weight-2 edge gets a "2" label; the weight-1 edge gets none."""
    sd = _weighted_symbolic_dynamics()

    plt.figure()
    try:
        ax = plt.gca()
        result = plotting.plot_transition_graph(sd, ax=ax)
        assert result is ax
        texts = [t.get_text() for t in ax.texts]
        assert texts.count("2") == 1
        assert "1" not in texts
    finally:
        plt.close()


def test_plot_transition_graph_rejects_an_unknown_kwarg():
    """A kwarg that is neither a style key nor a draw_networkx_nodes param raises."""
    sd = _weighted_symbolic_dynamics()

    plt.figure()
    try:
        ax = plt.gca()
        with pytest.raises(TypeError):
            plotting.plot_transition_graph(sd, ax=ax, not_a_real_kwarg=123)
    finally:
        plt.close()


def test_plot_transition_graph_style_kwarg_overrides_the_default():
    """A known TRANSITION_GRAPH_STYLE name overrides the default (no raise)."""
    sd = _weighted_symbolic_dynamics()

    plt.figure()
    try:
        ax = plt.gca()
        result = plotting.plot_transition_graph(sd, ax=ax, node_size=1200)
        assert result is ax
    finally:
        plt.close()


def test_session_plot_transition_graph_matches_plotting_module(henon_session):
    """TangleSession.plot_transition_graph delegates to plotting.plot_transition_graph."""
    session, _fp = henon_session
    sd = _weighted_symbolic_dynamics()

    plt.figure()
    try:
        ax = plt.gca()
        result = session.plot_transition_graph(sd, ax=ax)
        assert result is ax
    finally:
        plt.close()
