"""
Regression (BCG-1, WP6): blast children are image evidence, never letters.

Every blast adds bridges at the frontier with no registered image. They used to
be classed, lettered and fed to the transition matrix, so the alphabet grew
while the hole set stayed the same. D2: the alphabet is the classes of the
minimal trellis's bridges; every other bridge is evidence only, and evidence
never adds a symbol (a registered image naming a class outside the alphabet
used to leak a virtual ``new<n>`` into the graph); only a WORD may add a class,
by closing the alphabet under the word map (WP7). Rule-based: no count, letter
or word of the build is pinned.
"""

from __future__ import annotations

import pytest

from cases import SERIES, _blast, _partition

_BLASTS = 6


@pytest.fixture(scope="module")
def k28_blasted():
    """The k=2.8 zone blasted one iteration at a time, six times (~3 s)."""
    build, _ = SERIES["k28_series"]
    case = build()
    (zone,) = case.zones
    for _ in range(_BLASTS):
        _blast(case.session, zone, case.fixed_points, case.pips, diagnose=False)
    _partition(case.session)
    return case


def test_alphabet_is_the_minimal_trellis_classes(k28_blasted):
    session = k28_blasted.session
    minimal = session.minimal_trellis()
    kept = {frozenset(bid) for bid in minimal.kept_bridge_ids}
    table = session.bridge_classes()
    assert len(table) < len(minimal.classes), "the blasts left no evidence-only class"
    for entry in minimal.classes:
        has_kept = any(frozenset(m.bridge_id) in kept for m in entry.members)
        assert has_kept == (entry.bridge_class in table), entry.bridge_class.label


def test_words_and_matrix_stay_inside_the_alphabet(k28_blasted):
    session = k28_blasted.session
    dynamics = session.symbolic_dynamics()
    table = dynamics.table
    for cd in dynamics.classes.values():
        assert all(symbol.bridge_class in table for symbol in cd.symbols), cd.letter
        assert (cd.named_by is None) == (cd.bridge_class in session.bridge_classes())
    expected = [
        name
        for entry in table.active
        for name in (
            [child.name for child in dynamics.refined.get(entry.bridge_class, ())]
            or [dynamics.classes[entry.bridge_class].letter]
        )
    ]
    names, _matrix = dynamics.transition_matrix()
    assert sorted(names) == sorted(expected)
