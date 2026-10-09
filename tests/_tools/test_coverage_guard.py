"""The coverage guard matches lines by content, not by number (TEST-6, WP11)."""

from __future__ import annotations

import coverage_guard

PATH = "src/tanglepack/fake.py"

BASE = '''\
def f(x):
    if x:
        return 1
    return 2


def g():
    return 3
'''

# three lines inserted above f, nothing else changed
CANDIDATE = '''\
import math

TAU = 2 * math.pi
def f(x):
    if x:
        return 1
    return 2


def g():
    return 3

def h():
    return 4
'''


def _report(executed: list[int], missing: list[int]) -> dict:
    return {PATH: {"executed_lines": executed, "missing_lines": missing,
                   "summary": {"percent_covered": 100.0 * len(executed) / len(executed + missing)}}}


def _compare(base: dict, cand: dict) -> int:
    return coverage_guard.compare(base, cand, lambda _: BASE, lambda _: CANDIDATE)


def test_a_moved_covered_statement_is_not_a_miss(capsys) -> None:
    base = _report([1, 2, 3, 4, 7, 8], [])
    cand = _report([1, 3, 4, 5, 6, 7, 10, 11, 13], [14])
    assert _compare(base, cand) == 0
    assert "NEWLY MISSED" not in capsys.readouterr().out


def test_a_statement_newly_missed_is_reported_at_its_new_line(capsys) -> None:
    base = _report([1, 2, 3, 4, 7, 8], [])
    cand = _report([1, 3, 4, 5, 7, 10, 11, 13, 14], [6])
    assert _compare(base, cand) == 1
    assert "NEWLY MISSED: [6]" in capsys.readouterr().out


def test_a_new_uncovered_line_is_not_reported() -> None:
    base = _report([1, 2, 4, 7, 8], [3])
    cand = _report([1, 3, 4, 5, 7, 10, 11, 13], [6, 14])
    assert _compare(base, cand) == 0


def test_repeated_text_is_told_apart_by_occurrence() -> None:
    source = "def f():\n    pass\n    pass\n"
    keys = coverage_guard.line_keys(source, [2, 3])
    assert keys[2] != keys[3]
    assert keys[2][:2] == keys[3][:2] == ("f", "pass")
