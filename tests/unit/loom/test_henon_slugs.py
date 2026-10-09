"""Nested figure folders encode their blast counts (API-7, WP11)."""

from __future__ import annotations

import itertools

from tanglepack.examples.henon_cases import _nested_slug


def test_every_blast_count_pair_gets_its_own_folder():
    pairs = list(itertools.product(range(5), range(4)))
    assert len({_nested_slug(outer, inner) for outer, inner in pairs}) == len(pairs)


def test_the_heteroclinic_folder_differs_from_the_nested_one():
    for outer, inner in itertools.product(range(5), range(4)):
        nested = _nested_slug(outer, inner)
        assert nested.replace("nested", "heteroclinic", 1) != nested
