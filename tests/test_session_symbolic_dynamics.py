"""
TangleSession.symbolic_dynamics / describe_symbolic_dynamics and the plot
delegates.

Mirrors ``test_session_dual_graph.py``. Nothing here pins a registry id, a
bridge count or a fixture word: those are pinned once in ``tests/golden/``.
The cache contract and session = direct build are ``tests/facade/``.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from tanglepack.topology import plotting


# --------------------------------------------------------------------------- #
# k=10: the report
# --------------------------------------------------------------------------- #
def test_describe_symbolic_dynamics(k10_partitioned):
    session, fp = k10_partitioned
    text = session.describe_symbolic_dynamics([fp])
    assert text and isinstance(text, str)
    assert text == session.symbolic_dynamics([fp]).describe()
    for cd in session.symbolic_dynamics([fp]).classes.values():
        assert cd.letter in text


# --------------------------------------------------------------------------- #
# Plot delegates
# --------------------------------------------------------------------------- #
_PLOT_NAMES = (
    "plot_walk",
    "plot_bridges_by_class",
    "plot_transition_graph",
    "plot_itinerary_table",
    "walk_legend_handles",
)
_plots_missing = [name for name in _PLOT_NAMES if not hasattr(plotting, name)]
needs_plots = pytest.mark.skipif(
    bool(_plots_missing), reason=f"plotting is missing {_plots_missing}"
)


@needs_plots
def test_session_plot_delegates_draw_on_the_given_axes(k10_partitioned):
    """The delegates run on Agg and hand back what the plotting functions return."""
    import networkx
    from matplotlib.lines import Line2D
    from matplotlib.table import Table

    session, fp = k10_partitioned
    dyn = session.symbolic_dynamics([fp])
    active = next(cd for cd in dyn.classes.values() if cd.kind == "active")
    assert active.search is not None and active.search.walk is not None
    fig, (ax_graph, ax_table, ax_plane) = plt.subplots(1, 3)
    try:
        graph = session.plot_transition_graph(ax=ax_graph, fixed_points=[fp])
        assert isinstance(graph, networkx.DiGraph)
        assert isinstance(
            session.plot_transition_graph(dyn, ax=ax_graph, refined=False), networkx.DiGraph
        )
        table = session.plot_itinerary_table(ax=ax_table, fixed_points=[fp])
        assert isinstance(table, Table) and table.axes is ax_table
        assert isinstance(session.plot_itinerary_table(dyn, ax=ax_table), Table)
        colours = session.plot_bridges_by_class(ax=ax_plane, fixed_points=[fp])
        assert isinstance(colours, dict) and colours
        assert session.plot_bridges_by_class(dyn, ax=ax_plane, refined=False)
        assert session.plot_dual_graph(session.dual_graph([fp]), ax=ax_plane) is ax_plane
        line = session.plot_walk(active.search.walk, ax=ax_plane, fixed_points=[fp])
        assert isinstance(line, Line2D) and line.axes is ax_plane
        assert ax_plane.lines, "the bridges and the walk are drawn on the plane axes"
    finally:
        plt.close(fig)


@needs_plots
def test_session_plot_delegates_forward_to_plotting(k10_partitioned, monkeypatch):
    session, fp = k10_partitioned
    dyn = session.symbolic_dynamics()
    active = next(cd for cd in dyn.classes.values() if cd.kind == "active")
    walk = active.search.walk
    calls = []
    monkeypatch.setattr(
        plotting, "plot_transition_graph",
        lambda d, ax=None, **kw: calls.append(("graph", d, ax, kw)) or "graph",
    )
    monkeypatch.setattr(
        plotting, "plot_walk",
        lambda g, w, ax=None, **kw: calls.append(("walk", (g, w), ax, kw)) or "walk",
    )
    monkeypatch.setattr(
        plotting, "plot_bridges_by_class",
        lambda t, d, ax=None, **kw: calls.append(("bridges", (t, d), ax, kw)) or "bridges",
    )
    monkeypatch.setattr(
        plotting, "plot_itinerary_table",
        lambda d, ax=None, **kw: calls.append(("table", d, ax, kw)) or "table",
    )
    monkeypatch.setattr(
        plotting, "plot_dual_graph_cartoon",
        lambda g, d, ax=None, **kw: calls.append(("cartoon", (g, d), ax, kw)) or "cartoon",
    )
    assert session.plot_transition_graph(refined=False) == "graph"
    assert session.plot_walk(walk, color="k") == "walk"
    assert session.plot_bridges_by_class(show_inert=False) == "bridges"
    assert session.plot_itinerary_table(refined=True) == "table"
    assert calls[0][1] is dyn and calls[0][3] == {"refined": False}
    assert calls[1][1] == (session.dual_graph(), walk) and calls[1][3] == {"color": "k"}
    assert calls[2][1] == (session.trellis(), dyn) and calls[2][3] == {"show_inert": False}
    assert calls[3][1] is dyn and calls[3][3] == {"refined": True}
    assert session.plot_dual_graph_cartoon(show_walks=False) == "cartoon"
    assert calls[4][1] == (session.dual_graph(), dyn) and calls[4][3] == {"show_walks": False}
