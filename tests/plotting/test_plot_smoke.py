"""Plotting tier, smoke half: every plotter and session delegate runs headless.

Each test draws on the Agg backend and checks only that the call runs, hands
back its documented type and (where it is mathtext) that the canvas renders.
Nothing here pins a style, a colour, an angle, a font size, a z-order or a
legend wording (author decision 5, 2026-10-05): those are free to change.
The topological properties of the drawings are ``test_plot_topology.py``.
"""

from __future__ import annotations

from typing import Callable

import matplotlib

matplotlib.use("Agg")  # headless

import matplotlib.pyplot as plt  # noqa: E402
import networkx  # noqa: E402
import pytest  # noqa: E402
from matplotlib.artist import Artist  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.table import Table  # noqa: E402

from cases import build_k10  # noqa: E402
from tanglepack.topology import plotting  # noqa: E402

TRELLIS_PLOTTERS = (
    "plot_strong_pip_candidates",
    "plot_strong_pip",
    "plot_pseudoneighbors",
    "plot_holes",
    "plot_stable_partition",
)


@pytest.fixture
def fig_ax():
    """A fresh ``(figure, axes)`` pair, closed after the test."""
    fig, ax = plt.subplots()
    yield fig, ax
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Trellis plotters and the session fan-outs over them
# --------------------------------------------------------------------------- #
def test_trellis_plotters_draw_on_a_real_tangle(k10_partitioned, fig_ax) -> None:
    """Every ``Trellis.plot_*`` returns its handle(s) on the partitioned k=10 trellis."""
    session, fp = k10_partitioned
    trellis = session.trellis(fp)
    _fig, ax = fig_ax
    for name in TRELLIS_PLOTTERS:
        handle = getattr(trellis, name)(ax=ax)
        assert handle is not None, name
    assert trellis.plot_holes(ax=ax), "the k=10 trellis has holes to draw"


def test_session_plot_fanouts_draw_every_fixed_point(henon_p3_session, fig_ax) -> None:
    """Each session ``plot_*`` fan-out draws one handle per fixed point (two tangles)."""
    session, fp3, fp1, _zone = henon_p3_session
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()
    _fig, ax = fig_ax
    for name in ("plot_strong_pip_candidates", "plot_strong_pip",
                 "plot_pseudoneighbors", "plot_stable_partition"):
        handles = getattr(session, name)(ax=ax)
        assert len(handles) == 2, name
        assert len(getattr(session, name)(fp3, ax=ax)) == 1, name
    assert session.plot_holes(ax=ax), "both tangles have holes"


def test_session_plot_pseudoneighbors_computes_on_demand(k10_session, fig_ax) -> None:
    """``plot_pseudoneighbors`` computes the pairs first; ``plot_holes`` draws punched holes."""
    session, fp = k10_session
    _fig, ax = fig_ax
    assert len(session.plot_pseudoneighbors(ax=ax)) == 1
    session.trellis(fp).punch_holes()
    assert session.plot_holes(ax=ax)


def test_plot_stable_partition_rows_and_element_labels(k10_partitioned, fig_ax) -> None:
    """Homotopy and iterated rows on one axes, with row labels and one element label per interval."""
    session, _fp = k10_partitioned
    iterated = session.iterated_partition()
    results = session.homotopy_partition().as_list() + iterated.as_list()
    naming = session.dual_graph().naming
    fig, ax = fig_ax

    labels = [f"{r.side} {i}" for i, r in enumerate(results)]
    assert plotting.plot_stable_partition(results, ax=ax, labels=labels) is ax
    ax.cla()

    rows = iterated.as_list()

    def element_label(result, interval) -> str:
        """The iterated element's mathtext name."""
        ref = iterated.ref(result.branch_key, result.side, interval.element_id)
        return naming.name(ref).mathtext

    plotting.plot_stable_partition(rows, ax=ax)
    without = len(ax.texts)
    ax.cla()
    plotting.plot_stable_partition(rows, ax=ax, element_labels=element_label)
    assert len(ax.texts) == without + sum(len(r.intervals) for r in rows)
    fig.canvas.draw()
    ax.cla()
    plotting.plot_stable_partition(rows, ax=ax, element_labels=lambda r, i: "")
    assert len(ax.texts) == without


# --------------------------------------------------------------------------- #
# Dual graph, minimal trellis and walks
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"show_labels": True},
        {"show_labels": True, "label_style": "ref"},
        {"clip_to_arcs": False, "color": "tab:green", "s": 30},
    ],
    ids=["default", "names", "refs", "unclipped_restyled"],
)
def test_plot_dual_graph_draws(k10_partitioned, fig_ax, kwargs) -> None:
    """``plot_dual_graph`` returns the axes under every label / clip option; mathtext renders."""
    session, _fp = k10_partitioned
    fig, ax = fig_ax
    assert plotting.plot_dual_graph(session.dual_graph(), ax=ax, **kwargs) is ax
    fig.canvas.draw()


def test_plot_minimal_trellis_draws(k10_partitioned, fig_ax) -> None:
    """``plot_minimal_trellis`` returns the axes, with and without the dropped bridges."""
    session, _fp = k10_partitioned
    _fig, ax = fig_ax
    minimal = session.minimal_trellis()
    assert plotting.plot_minimal_trellis(minimal, ax=ax) is ax
    assert plotting.plot_minimal_trellis(minimal, ax=ax, show_dropped=False) is ax


def test_plot_walk_accepts_style_overrides(k10_partitioned, fig_ax) -> None:
    """``plot_walk`` returns a ``Line2D`` and takes line-style overrides."""
    session, _fp = k10_partitioned
    dynamics = session.symbolic_dynamics()
    walk = next(cd.search.walk for cd in dynamics.classes.values()
                if cd.search is not None and cd.search.walk is not None)
    fig, ax = fig_ax
    line = plotting.plot_walk(session.dual_graph(), walk, ax=ax,
                              color="tab:purple", linewidth=4, zorder=3)
    assert isinstance(line, Line2D)
    fig.canvas.draw()


# --------------------------------------------------------------------------- #
# Symbolic-dynamics views
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("refined", [True, False])
def test_plot_transition_graph_and_bridges_by_class_draw(k10_partitioned, fig_ax, refined) -> None:
    """The transition graph is a ``DiGraph``; the bridges-by-class plot returns its colour map."""
    session, _fp = k10_partitioned
    dynamics = session.symbolic_dynamics()
    fig, ax = fig_ax
    assert isinstance(plotting.plot_transition_graph(dynamics, ax=ax, refined=refined),
                      networkx.DiGraph)
    ax.cla()
    colors = plotting.plot_bridges_by_class(session.trellis(), dynamics, ax=ax, refined=refined)
    assert isinstance(colors, dict) and colors
    fig.canvas.draw()


@pytest.mark.parametrize(
    "figsize, kwargs",
    [
        ((6.4, 4.8), {}),
        ((6.4, 4.8), {"mathtext": False}),
        ((6.4, 4.8), {"refined": False}),
        ((16, 6), {"columns": ("word", "class", "refined word"), "fontsize": 20.0,
                   "row_height": 3.0}),
        ((2.5, 3), {"fontsize": 20.0}),
    ],
    ids=["default", "plain", "unrefined", "columns_wide", "narrow"],
)
def test_plot_itinerary_table_draws_one_row_per_class(k10_partitioned, figsize, kwargs) -> None:
    """``plot_itinerary_table`` returns a ``Table`` with a header plus one row per class."""
    session, _fp = k10_partitioned
    dynamics = session.symbolic_dynamics()
    fig, ax = plt.subplots(figsize=figsize)
    try:
        table = plotting.plot_itinerary_table(dynamics, ax=ax, **kwargs)
        assert isinstance(table, Table)
        cells = table.get_celld()
        assert 1 + max(row for row, _ in cells) == 1 + len(dynamics.classes)
        if "columns" in kwargs:
            assert 1 + max(col for _, col in cells) == len(kwargs["columns"])
        fig.canvas.draw()
    finally:
        plt.close(fig)


def test_itinerary_table_rows_follow_the_requested_columns(k10_partitioned) -> None:
    """``columns`` picks a subset of the full rows, in the caller's order."""
    session, _fp = k10_partitioned
    dynamics = session.symbolic_dynamics()
    full = plotting.itinerary_table_rows(dynamics, mathtext=False)
    index = {name: i for i, name in enumerate(plotting.ITINERARY_TABLE_COLUMNS)}
    columns = ("word", "class", "refined word")
    rows = plotting.itinerary_table_rows(dynamics, mathtext=False, columns=columns)
    assert len(rows) == len(full) == len(dynamics.classes)
    for row, full_row in zip(rows, full):
        assert row == [full_row[index[name]] for name in columns]
    for row in plotting.itinerary_table_rows(dynamics):
        assert all(isinstance(cell, str) for cell in row)


# --------------------------------------------------------------------------- #
# The dual-graph cartoon
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"anchor": "left"},
        {"show_walks": False, "bridge_color": "tab:purple", "bridge_alpha": 0.7,
         "label_fontsize": 22.0},
        {"label_position": "outside"},
        {"show_walks": False, "bridge_linewidth": 1.2},
        {"color_bridges": True},
        {"shape": "circle"},
        {"shape": "circle", "label_position": "outside", "color_bridges": True},
    ],
    ids=["default", "anchor_left", "recoloured", "outside_labels", "uniform_width",
         "class_coloured", "circle", "circle_outside_coloured"],
)
def test_plot_dual_graph_cartoon_draws(k10_partitioned, fig_ax, kwargs) -> None:
    """The cartoon returns its layout under every option, line or circle."""
    session, _fp = k10_partitioned
    fig, ax = fig_ax
    layout = plotting.plot_dual_graph_cartoon(
        session.dual_graph(), session.symbolic_dynamics(), ax=ax, **kwargs
    )
    expected = plotting.ZoneLayout if kwargs.get("shape") == "circle" else plotting.CartoonLayout
    assert isinstance(layout, expected)
    fig.canvas.draw()


# --------------------------------------------------------------------------- #
# Legend proxies
# --------------------------------------------------------------------------- #
def test_legend_handles_are_artists(k10_partitioned) -> None:
    """Every legend-handle helper returns a non-empty list of artists, under every option."""
    session, _fp = k10_partitioned
    dynamics = session.symbolic_dynamics()
    lists = [
        plotting.dual_graph_legend_handles(),
        plotting.minimal_trellis_legend_handles(),
        plotting.walk_legend_handles(),
        plotting.walk_legend_handles(color="tab:purple"),
        plotting.bridge_class_legend_handles(dynamics),
        plotting.bridge_class_legend_handles(dynamics, refined=False),
        plotting.dual_graph_cartoon_legend_handles(),
        plotting.dual_graph_cartoon_legend_handles(dynamics),
        plotting.dual_graph_cartoon_legend_handles(
            dynamics, bridge_color="tab:purple", bridge_alpha=0.7, walks=False
        ),
        plotting.dual_graph_cartoon_legend_handles(bridge_linewidth=1.2, walks=False),
        plotting.dual_graph_cartoon_legend_handles(walks=False, empty_side=True),
    ]
    for handles in lists:
        assert handles and all(isinstance(h, Artist) for h in handles)


# --------------------------------------------------------------------------- #
# Session plot delegates
# --------------------------------------------------------------------------- #
def test_session_plot_delegates_return_their_types(k10_partitioned) -> None:
    """The session's symbolic-dynamics and dual-graph delegates draw on the given axes.

    Forwarded keywords are checked by their effect, not by spying on the
    plotting module.
    """
    session, fp = k10_partitioned
    dyn = session.symbolic_dynamics([fp])
    walk = next(cd.search.walk for cd in dyn.classes.values()
                if cd.search is not None and cd.search.walk is not None)
    fig, (ax_graph, ax_table, ax_plane) = plt.subplots(1, 3)
    try:
        graph = session.plot_transition_graph(ax=ax_graph, fixed_points=[fp])
        assert isinstance(graph, networkx.DiGraph)
        unrefined = session.plot_transition_graph(dyn, ax=ax_graph, refined=False)
        assert set(unrefined.nodes) == set(dyn.transition_graph(refined=False).nodes)

        table = session.plot_itinerary_table(ax=ax_table, fixed_points=[fp])
        assert isinstance(table, Table) and table.axes is ax_table
        assert isinstance(session.plot_itinerary_table(dyn, ax=ax_table, refined=True), Table)

        colours = session.plot_bridges_by_class(ax=ax_plane, fixed_points=[fp])
        assert isinstance(colours, dict) and colours
        assert session.plot_bridges_by_class(dyn, ax=ax_plane, refined=False)
        assert session.plot_dual_graph(session.dual_graph([fp]), ax=ax_plane) is ax_plane
        assert session.plot_dual_graph(ax=ax_plane, show_labels=True) is ax_plane
        assert session.plot_minimal_trellis(ax=ax_plane) is ax_plane
        line = session.plot_walk(walk, ax=ax_plane, fixed_points=[fp], color="k")
        assert isinstance(line, Line2D) and line.axes is ax_plane
        assert ax_plane.lines, "the bridges and the walk are drawn on the plane axes"
    finally:
        plt.close(fig)


def test_session_cartoon_delegate_honours_its_toggles(k10_partitioned, fig_ax) -> None:
    """With every toggle off, the session cartoon draws no bridge, walk or crossing."""
    session, fp = k10_partitioned
    _fig, ax = fig_ax
    layout = session.plot_dual_graph_cartoon(
        ax=ax, fixed_points=[fp], show_bridges=False, show_walks=False, show_unified=False,
    )
    assert isinstance(layout, plotting.CartoonLayout)
    assert layout.bridges_drawn == 0 and layout.walks_drawn == 0
    assert not ax.patches


# --------------------------------------------------------------------------- #
# Rejected arguments (exception type only)
# --------------------------------------------------------------------------- #
def _dual(session):
    """The session's dual graph."""
    return session.dual_graph()


BAD_CALLS: dict[str, Callable] = {
    "dual_graph_facecolors": lambda s, ax: plotting.plot_dual_graph(_dual(s), ax=ax, facecolors="red"),
    "dual_graph_c": lambda s, ax: plotting.plot_dual_graph(_dual(s), ax=ax, c="red"),
    "dual_graph_label_style": lambda s, ax: plotting.plot_dual_graph(
        _dual(s), ax=ax, show_labels=True, label_style="bogus"),
    "bridges_by_class_color": lambda s, ax: plotting.plot_bridges_by_class(
        s.trellis(), s.symbolic_dynamics(), ax=ax, color="red"),
    "bridges_by_class_linestyle": lambda s, ax: plotting.plot_bridges_by_class(
        s.trellis(), s.symbolic_dynamics(), ax=ax, linestyle=":"),
    "table_rows_columns": lambda s, ax: plotting.itinerary_table_rows(
        s.symbolic_dynamics(), columns=("class", "bogus")),
    "table_columns": lambda s, ax: plotting.plot_itinerary_table(
        s.symbolic_dynamics(), ax=ax, columns=("class", "bogus")),
    "stable_partition_labels": lambda s, ax: plotting.plot_stable_partition(
        s.homotopy_partition().as_list(), ax=ax, labels=["only one"]),
    "cartoon_layout_anchor": lambda s, ax: plotting.dual_graph_cartoon_layout(_dual(s), anchor="top"),
    "cartoon_shape": lambda s, ax: plotting.plot_dual_graph_cartoon(
        _dual(s), s.symbolic_dynamics(), ax=ax, shape="hexagon"),
    "cartoon_label_position": lambda s, ax: plotting.plot_dual_graph_cartoon(
        _dual(s), s.symbolic_dynamics(), ax=ax, label_position="above"),
    "circle_label_position": lambda s, ax: plotting.plot_dual_graph_cartoon(
        _dual(s), s.symbolic_dynamics(), ax=ax, shape="circle", label_position="nowhere"),
}


@pytest.mark.parametrize("name", list(BAD_CALLS))
def test_plotters_reject_bad_arguments(name, fig_ax) -> None:
    """Reserved style keywords and unknown option values raise ``ValueError``."""
    session = build_k10().session
    _fig, ax = fig_ax
    with pytest.raises(ValueError):
        BAD_CALLS[name](session, ax)
