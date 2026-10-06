"""Log assertions by level and logger only (author decision 5: never pin wording)."""

from __future__ import annotations

import logging
from typing import Optional, Union

import pytest


def assert_logged(
    caplog: pytest.LogCaptureFixture,
    level: Union[int, str],
    logger: str,
    *,
    count: Optional[int] = None,
) -> list[logging.LogRecord]:
    """
    Assert that ``logger`` emitted at least one record at exactly ``level``.

    Args:
        caplog: pytest's log-capture fixture (set its level before the call
            under test, e.g. ``caplog.set_level(logging.WARNING, logger=...)``).
        level: The level, as a number or a name (``"WARNING"``).
        logger: The logger name, e.g. ``"tanglepack.topology.DualWalk"``.
            Records of child loggers (``logger + ".x"``) count too.
        count: When given, the exact number of matching records.

    Returns:
        The matching records, for further structural checks.

    Raises:
        AssertionError: No record (or not ``count`` records) matched.
    """
    levelno = logging.getLevelName(level) if isinstance(level, str) else level
    assert isinstance(levelno, int), f"unknown log level {level!r}"
    records = [
        record
        for record in caplog.records
        if record.levelno == levelno
        and (record.name == logger or record.name.startswith(logger + "."))
    ]
    if count is None:
        assert records, (
            f"no {logging.getLevelName(levelno)} record from {logger!r}; got "
            f"{[(r.name, r.levelname) for r in caplog.records]}"
        )
    else:
        assert len(records) == count, (
            f"{len(records)} {logging.getLevelName(levelno)} record(s) from "
            f"{logger!r}, expected {count}"
        )
    return records
