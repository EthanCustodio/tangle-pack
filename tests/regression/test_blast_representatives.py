"""The D1 blast series under a representative blast.

Each series (:data:`cases.SERIES`) is blasted one iteration at a time to its
own depth with ``representatives="fewest_points"`` (one bridge per refined
class). The D1 series laws of ``test_class_count_under_blasting.py`` must
still hold on its records (:func:`cases.blast_series`): the class count never
decreases and is constant while the hole set is, so is the closed alphabet,
every class has a member the minimal trellis keeps, lost references are
witnessed, letters follow identity, the session trellis is history free, no
walk is ambiguous, a recompute keeps N and the cut, no bridge ends outside N,
the pip orbit never gets shallower and the trusted depth is within it. After
every blast the alphabet stays closed: every symbol of every word is a class.

Equivalence with the full series (``pip_depth`` included) is not asserted:
the representative blast holds bridges, so its X is thinner and its topology
may arrive later; ``scripts/blast_representatives.py`` measures it.
"""

from __future__ import annotations

from functools import cache

import pytest

from cases import SERIES, SeriesStep, blast_series

POLICY = "fewest_points"


@cache
def _series(name: str) -> tuple[SeriesStep, ...]:
    """The records of one representative series, built once."""
    return blast_series(name, representatives=POLICY)


def _violations(steps: tuple[SeriesStep, ...]) -> list[tuple]:
    """``(law, blasts)`` of every D1 series law a run of records breaks."""
    found = []
    for before, after in zip(steps, steps[1:]):
        if before.classes > after.classes:
            found.append(("class count decreases", after.blasts))
        if before.holes == after.holes and before.classes != after.classes:
            found.append(("class count moves under one hole set", after.blasts))
        if before.holes == after.holes and before.closed != after.closed:
            found.append(("closed alphabet moves under one hole set", after.blasts))
        if before.pip_depth > after.pip_depth:
            found.append(("pip orbit shallower", after.blasts))
    for step in steps:
        checks = {
            "class outside the minimal trellis": step.outside_minimal == 0,
            "lost reference unwitnessed": not step.unwitnessed,
            "session trellis depends on history": step.history_free,
            "ambiguous walk": step.ambiguous == 0,
            "recompute changes N or the cut": step.crossings == step.recomputed,
            "bridge ends outside N": all(
                set(bridge_id) <= crossings
                for crossings, bridges in (step.crossings, step.recomputed)
                for bridge_id in bridges
            ),
            "trusted depth beyond pip depth": step.blast_depths is None
            or step.blast_depths[1] <= step.blast_depths[0] == step.pip_depth,
        }
        found += [(law, step.blasts) for law, holds in checks.items() if not holds]
    letters: dict = {}
    for step in steps:
        for identity, letter in step.letters:
            if letters.setdefault(identity, letter) != letter:
                found.append(("letter changes", step.blasts))
    if len(set(letters.values())) != len(letters):
        found.append(("two identities share a letter", steps[-1].blasts))
    return found


@pytest.mark.parametrize("name", sorted(SERIES))
def test_the_series_laws_hold_under_representatives(name: str) -> None:
    """Every D1 law of the full series holds on the representative series."""
    steps = _series(name)
    assert len(steps) == SERIES[name].blasts + 1
    assert not _violations(steps)
    same = [(a, b) for a, b in zip(steps, steps[1:]) if a.holes == b.holes]
    assert same, "no two consecutive recomputes share a hole set: the hole laws checked nothing"


@pytest.mark.parametrize("name", sorted(SERIES))
def test_the_alphabet_stays_closed_under_representatives(name: str) -> None:
    """After every representative blast each word's symbols are classes of the table."""
    series = SERIES[name]
    case = series.build()
    session = case.session
    zone = max(case.zones, key=lambda zone: zone.area)
    for _ in range(series.blasts):
        result = session.blast_zone(
            zone, 1, fixed_point=[zone.fixed_point], min_separation=series.min_separation,
            diagnose=False, representatives=POLICY,
        )
        assert result.completed_iterations == 1 or result.stop_reason == "frontier_empty"
        dynamics = session.symbolic_dynamics()
        assert dynamics.rules
        for letter, word in dynamics.rules.items():
            assert all(symbol.bridge_class in dynamics.table for symbol in word), letter
