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
