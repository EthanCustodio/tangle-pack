"""The ``verbose=`` policy: reports go through ``logging``, never a bare print.

A pipeline method called with ``verbose=True`` emits its ``describe_*``
report at INFO on its module logger. Only when the application has configured
no logging at all does the report fall back to stdout (``logging.lastResort``
would drop an INFO record). Assertions check the level, the logger and that
something was reported, never the wording (author decision 5, E6).
"""

from __future__ import annotations

import logging

import pytest

from helpers.logs import assert_logged

TRELLIS_LOGGER = "tanglepack.topology.Trellis"


@pytest.fixture
def trellis(k10_partitioned):
    """The partitioned k=10 trellis (pseudoneighbors, holes and partitions computed)."""
    session, fp = k10_partitioned
    return session.trellis(fp)


@pytest.mark.parametrize(
    "method, describe",
    [
        ("compute_pseudoneighbors", "describe_pseudoneighbors"),
        ("punch_holes", "describe_holes"),
        ("partition_stable_manifold", "describe_stable_partitions"),
    ],
)
def test_verbose_logs_and_never_prints(trellis, capsys, caplog, method, describe) -> None:
    """``verbose=True`` emits the ``describe_*`` report at INFO, not on stdout."""
    capsys.readouterr()
    with caplog.at_level(logging.INFO, logger=TRELLIS_LOGGER):
        getattr(trellis, method)(verbose=True)

    assert capsys.readouterr().out == ""
    assert getattr(trellis, describe)()
    assert_logged(caplog, logging.INFO, TRELLIS_LOGGER)


def test_verbose_prints_when_logging_is_unconfigured(trellis, capsys) -> None:
    """With no logging handler at all, ``verbose=True`` still reaches stdout.

    ``logging.lastResort`` only emits at WARNING, so an INFO record from an
    unconfigured application would vanish: the report is printed instead.
    """
    root = logging.getLogger()
    saved_handlers, saved_level = root.handlers[:], root.level
    root.handlers = []
    try:
        assert not root.hasHandlers(), "the test must leave logging unconfigured"
        capsys.readouterr()
        trellis.punch_holes(verbose=True)
        out = capsys.readouterr().out
    finally:
        root.handlers = saved_handlers
        root.setLevel(saved_level)

    assert out.strip()


def test_verbose_survives_an_application_level_above_info(trellis, caplog) -> None:
    """``verbose=True`` forces the report through a WARNING-level logger and restores the level."""
    logger = logging.getLogger(TRELLIS_LOGGER)
    previous = logger.level
    logger.setLevel(logging.WARNING)
    try:
        with caplog.at_level(logging.INFO):
            trellis.punch_holes(verbose=True)
        assert_logged(caplog, logging.INFO, TRELLIS_LOGGER)
        assert logger.level == logging.WARNING, "the call must restore the application's level"
    finally:
        logger.setLevel(previous)
