"""Figure folders encode every non-default build parameter (API-7, WP11; DOC-4, WP12),
and the example inversion recipe is the law cases' (WP12)."""

from __future__ import annotations

import itertools

import pytest

import cases
from tanglepack.examples import henon_cases
from tanglepack.examples.henon_cases import _nested_slug, _period3_slug


def test_every_blast_count_pair_gets_its_own_folder():
    pairs = list(itertools.product(range(5), range(4)))
    assert len({_nested_slug(outer, inner) for outer, inner in pairs}) == len(pairs)


def test_the_heteroclinic_folder_differs_from_the_nested_one():
    for outer, inner in itertools.product(range(5), range(4)):
        nested = _nested_slug(outer, inner)
        assert nested.replace("nested", "heteroclinic", 1) != nested


def test_every_period3_build_gets_its_own_folder():
    grid = list(itertools.product((11, 13, 15, 16), (9, 10, 12), (0, 1, 8)))
    assert len({_period3_slug(*params) for params in grid}) == len(grid)


@pytest.mark.parametrize("blasts", [0, 1])
def test_the_example_inversion_recipe_is_the_law_cases(blasts):
    """``henon_cases.build_inversion`` (the figures) and ``cases.build_inversion``
    (the law cases) build the same topology: a drift guard, nothing pinned
    (registry ids are reproducible, WP8, so the reports compare exactly)."""
    example, law = henon_cases.build_inversion(blasts=blasts).session, cases.build_inversion(blasts=blasts).session
    assert example.symbolic_dynamics().describe() == law.symbolic_dynamics().describe()
    holes = [
        sorted((h.bounding_ids, h.bridge_side) for h in session.topology().trellis.holes)
        for session in (example, law)
    ]
    assert holes[0] == holes[1]


def test_every_inversion_build_gets_its_own_folder():
    slugs = {henon_cases.build_inversion(blasts=blasts).slug for blasts in range(3)}
    assert len(slugs) == 3
