"""Plotting tier, rule half: the colour family of every class symbol.

CLAUDE.md (``class_colors``): no colour is used twice; one tangle draws from
:data:`~tanglepack.topology.plotting.CLASS_COLORS`; several tangles take one
:data:`~tanglepack.topology.plotting.TANGLE_COLOR_FAMILIES` family each
(distinct tangles, distinct families) and the classes connecting two tangles
(``BridgeClassEntry.tangle is None``) draw from
:data:`~tanglepack.topology.plotting.HETEROCLINIC_COLORS`. The author asked
for this as a rule (2026-10-05), not a style pin: the tests check FAMILY
MEMBERSHIP only, never a colour value or which slot of a family a symbol gets,
plus the rule that two families' colours never look alike (a CIE76 floor).

A family that runs out is extended from its colormap rather than repeating,
so a group larger than its family is only required to exhaust that family
first (and the extension colours to belong to no other family).
"""

from __future__ import annotations

import itertools
from types import SimpleNamespace
from typing import Optional

import numpy as np
from matplotlib.colors import to_hex, to_rgb
import pytest

from cases import build_k10, build_nested
from tanglepack.topology import plotting

_FAMILIES = [frozenset(to_hex(c) for c in family) for family in plotting.TANGLE_COLOR_FAMILIES]
_HETEROCLINIC = frozenset(to_hex(c) for c in plotting.HETEROCLINIC_COLORS)
_SINGLE = frozenset(to_hex(c) for c in plotting.CLASS_COLORS)


def _symbol_tangles(dynamics, *, refined: bool) -> dict[str, Optional[int]]:
    """``{symbol: tangle or None}``: refined children inherit their class's tangle."""
    tangles: dict[str, Optional[int]] = {}
    for bridge_class, cd in dynamics.classes.items():
        children = dynamics.refined.get(bridge_class, []) if refined else []
        for name in [c.name for c in children] or [cd.letter]:
            tangles[name] = cd.entry.tangle
    return tangles


def _assert_draws_from(colors: list[str], family: frozenset, others: frozenset) -> None:
    """``colors`` exhaust ``family`` before extending, and never borrow from ``others``."""
    inside = [c for c in colors if c in family]
    assert len(inside) == min(len(colors), len(family)), (colors, family)
    assert not set(colors) & others, set(colors) & others


def _assert_family_rule(dynamics, *, refined: bool) -> None:
    """Every tangle in its own family, the connecting classes heteroclinic, no repeats."""
    palette = {name: to_hex(color) for name, color in
               plotting.class_colors(dynamics, refined=refined).items()}
    tangles = _symbol_tangles(dynamics, refined=refined)
    assert set(palette) == set(tangles)
    assert len(set(palette.values())) == len(palette), "a colour is used twice"

    groups: dict[Optional[int], list[str]] = {}
    for name, tangle in tangles.items():
        groups.setdefault(tangle, []).append(palette[name])

    owner: dict[int, Optional[int]] = {}  # family index -> tangle
    for tangle, colors in groups.items():
        if tangle is None:
            _assert_draws_from(colors, _HETEROCLINIC, frozenset().union(*_FAMILIES))
            continue
        hits = {i for i, family in enumerate(_FAMILIES) if set(colors) & family}
        assert len(hits) == 1, f"tangle {tangle} draws from families {sorted(hits)}"
        (family_index,) = hits
        assert family_index not in owner, (
            f"tangles {owner[family_index]} and {tangle} share family {family_index}")
        owner[family_index] = tangle
        others = _HETEROCLINIC.union(*(f for i, f in enumerate(_FAMILIES) if i != family_index))
        _assert_draws_from(colors, _FAMILIES[family_index], others)


# --------------------------------------------------------------------------- #
# Hand-built class tables (connecting classes may not exist on real builds)
# --------------------------------------------------------------------------- #
def _dynamics(spec: list[tuple[str, Optional[int], int]]) -> SimpleNamespace:
    """
    A stand-in symbolic dynamics carrying only what ``class_colors`` reads.

    Args:
        spec: ``(letter, tangle, refined children)`` per class, in table order;
            ``tangle=None`` marks a class connecting two tangles.
    """
    classes, refined = {}, {}
    for position, (letter, tangle, children) in enumerate(spec):
        key = ("class", position)
        classes[key] = SimpleNamespace(
            letter=letter, kind="active",
            entry=SimpleNamespace(tangle=tangle, members=()))
        if children:
            refined[key] = [SimpleNamespace(name=f"{letter}_{i}", index=i)
                            for i in range(1, children + 1)]
    return SimpleNamespace(classes=classes, refined=refined, unmatched_members={})


_TWO_TANGLES_AND_CONNECTING = [
    ("a", 0, 2), ("b", 0, 0), ("c", 0, 0),
    ("d", 1, 0), ("e", 1, 3),
    ("f", None, 0), ("g", None, 2),
]


@pytest.mark.parametrize("refined", [True, False])
def test_hand_built_two_tangles_and_connecting_classes(refined) -> None:
    """Two tangles take two distinct families; connecting classes are heteroclinic."""
    _assert_family_rule(_dynamics(_TWO_TANGLES_AND_CONNECTING), refined=refined)


def test_hand_built_family_rule_ignores_table_order() -> None:
    """Interleaving the tangles in the table changes no symbol's family."""
    order = [5, 0, 3, 1, 6, 4, 2]
    _assert_family_rule(_dynamics([_TWO_TANGLES_AND_CONNECTING[i] for i in order]),
                        refined=True)


def test_hand_built_one_tangle_with_connecting_classes_still_splits_families() -> None:
    """A connecting class alongside one tangle: the tangle keeps a family of its own."""
    _assert_family_rule(_dynamics([("a", 0, 0), ("b", 0, 2), ("c", None, 0)]), refined=True)


def test_hand_built_family_overflow_extends_without_repeating() -> None:
    """A tangle with more symbols than its family exhausts it, then extends uniquely."""
    big = [(f"s{i}", 0, 0) for i in range(len(plotting.TANGLE_COLOR_FAMILIES[0]) + 3)]
    _assert_family_rule(_dynamics(big + [("t", 1, 0), ("h", None, 0)]), refined=True)


def test_hand_built_three_tangles_take_three_families() -> None:
    """Three tangles take three distinct families (TEST-5(c), WP9)."""
    spec = [("a", 0, 0), ("b", 0, 0), ("c", 1, 0), ("d", 1, 0), ("e", 2, 0), ("f", 2, 0),
            ("g", None, 0)]
    _assert_family_rule(_dynamics(spec), refined=True)


#: The smallest CIE76 colour difference allowed between two families' colours.
MIN_FAMILY_DELTA_E = 13.0


def _lab(color: str) -> np.ndarray:
    """CIE L*a*b* (D65) of a colour."""
    rgb = np.array(to_rgb(color))
    rgb = np.where(rgb > 0.04045, ((rgb + 0.055) / 1.055) ** 2.4, rgb / 12.92)
    xyz = np.array([[0.4124, 0.3576, 0.1805], [0.2126, 0.7152, 0.0722],
                    [0.0193, 0.1192, 0.9505]]) @ rgb / np.array([0.95047, 1.0, 1.08883])
    f = np.where(xyz > 0.008856, np.cbrt(xyz), 7.787 * xyz + 16 / 116)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])


def test_families_do_not_look_alike() -> None:
    """No colour of one family (or the connecting colours) is perceptually close
    to another's: "families never share a colour" means by eye, not by hex."""
    groups = [*plotting.TANGLE_COLOR_FAMILIES, plotting.HETEROCLINIC_COLORS]
    for one, other in itertools.combinations(groups, 2):
        for a, b in itertools.product(one, other):
            assert np.linalg.norm(_lab(a) - _lab(b)) >= MIN_FAMILY_DELTA_E, (a, b)


def test_hand_built_single_tangle_draws_from_class_colors() -> None:
    """One tangle and no connecting class: every colour from ``CLASS_COLORS``."""
    dynamics = _dynamics([("a", 0, 2), ("b", 0, 0), ("u", 0, 0)])
    palette = {n: to_hex(c) for n, c in plotting.class_colors(dynamics).items()}
    assert set(palette.values()) <= _SINGLE
    assert len(set(palette.values())) == len(palette)


# --------------------------------------------------------------------------- #
# Real builds
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("refined", [True, False])
def test_nested_tangles_each_draw_from_their_own_family(refined) -> None:
    """Period 1 + period 3: each tangle's symbols in one family, no family shared."""
    dynamics = build_nested().session.symbolic_dynamics()
    tangles = {cd.entry.tangle for cd in dynamics.classes.values()}
    assert len(tangles - {None}) >= 2, "the nested build should span two tangles"
    _assert_family_rule(dynamics, refined=refined)


def test_k10_single_tangle_draws_from_class_colors() -> None:
    """k=10 is one tangle: its symbols draw from ``CLASS_COLORS`` (none extended)."""
    dynamics = build_k10().session.symbolic_dynamics()
    palette = {n: to_hex(c) for n, c in plotting.class_colors(dynamics).items()}
    assert len(palette) <= len(_SINGLE)
    assert set(palette.values()) <= _SINGLE
