"""The D1 loop: blast, recompute the topology from scratch, compare.

The author's decisions of 2026-10-07: topology is a pure function of the
numerical state, recomputed from scratch after every blast (D1), and the class
count never decreases -- it plateaus and may jump again after much more
growth, so a decrease, or growth while nothing else changes, is a bug. The
alphabet is the minimal trellis's bridges only (D2).

Each series (:data:`cases.SERIES`) blasts one zone once at a time and records
a :class:`cases.SeriesStep` per recompute. The laws are rule-based, never a
pinned count:

* the TOTAL class count (active + inert) never decreases -- D1 says "the
  class count"; per kind it legitimately moves: a first blast registers image
  loops that turn active classes inert, and a reference-pair switch can merge
  inert ones (``docs/open_issues.md``, TEST-2);
* while the hole set stays the same between two recomputes, so does the class
  count (at least one such pair must exist, or the law checked nothing);
* every class has a member the minimal trellis keeps.

The per-kind counts are printed in every failure message as diagnostics.

``_series`` caches the immutable records of each series once per run; no
session survives the call, so no build state is shared between tests.
"""

from __future__ import annotations

from functools import cache

import pytest

from cases import SERIES, SeriesStep, blast_series, issue_marks


@cache
def _series(name: str) -> tuple[SeriesStep, ...]:
    """The records of one series, built once."""
    return blast_series(name)


def _table(steps: tuple[SeriesStep, ...]) -> str:
    """``blasts: holes, classes (active/inert), outside minimal`` per step."""
    return "; ".join(
        f"{s.blasts}: {len(s.holes)} holes, {s.classes} ({s.active}/{s.inert}), "
        f"{s.outside_minimal} outside"
        for s in steps
    )


def _params(law: str) -> list:
    """One param per series, carrying the law's known issue on it."""
    return [pytest.param(name, marks=issue_marks(law, name), id=name) for name in SERIES]


@pytest.mark.parametrize("name", _params("series_class_count_never_decreases"))
def test_class_count_never_decreases(name: str) -> None:
    """D1: one more blast never loses a class."""
    steps = _series(name)
    for before, after in zip(steps, steps[1:]):
        assert before.classes <= after.classes, _table(steps)


@pytest.mark.parametrize("name", _params("series_class_count_constant_while_holes_constant"))
def test_class_count_constant_while_holes_constant(name: str) -> None:
    """D1 + D2: the same hole set gives the same class count."""
    steps = _series(name)
    same = [(a, b) for a, b in zip(steps, steps[1:]) if a.holes == b.holes]
    assert same, f"no two consecutive recomputes share a hole set: {_table(steps)}"
    for before, after in same:
        assert before.classes == after.classes, _table(steps)


@pytest.mark.parametrize("name", _params("series_alphabet_inside_minimal_trellis"))
def test_alphabet_inside_minimal_trellis_under_blasting(name: str) -> None:
    """D2 after every blast: every class has a member the minimal trellis keeps."""
    steps = _series(name)
    assert all(s.outside_minimal == 0 for s in steps), _table(steps)
