"""
TangleSession.dual_graph / minimal_trellis / iterated_partition: describe and
the plot delegates.

The cache contract is ``tests/facade/test_session_caches.py``; session = direct
build is ``tests/facade/test_session_equivalence.py``.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from tanglepack.topology import plotting


def test_describe_iterated_partition(k10_partitioned):
    session, fp = k10_partitioned
    assert session.describe_iterated_partition() == session.iterated_partition().describe()


def test_session_plot_dual_graph_returns_axes(k10_partitioned):
    session, fp = k10_partitioned
    fig, ax = plt.subplots()
    try:
        assert session.plot_dual_graph(ax=ax) is ax
        assert session.plot_minimal_trellis(ax=ax) is ax
    finally:
        plt.close(fig)


def test_session_plot_delegates_to_plotting(k10_partitioned, monkeypatch):
    session, fp = k10_partitioned
    calls = []
    monkeypatch.setattr(
        plotting, "plot_dual_graph",
        lambda graph, ax=None, **kw: calls.append(("dual", graph, ax, kw)) or "dual",
    )
    monkeypatch.setattr(
        plotting, "plot_minimal_trellis",
        lambda minimal, ax=None, **kw: calls.append(("minimal", minimal, ax, kw)) or "minimal",
    )
    assert session.plot_dual_graph(show_labels=True) == "dual"
    assert session.plot_minimal_trellis(show_dropped=False) == "minimal"
    assert calls[0][1] is session.dual_graph() and calls[0][3] == {"show_labels": True}
    assert calls[1][1] is session.minimal_trellis() and calls[1][3] == {"show_dropped": False}
