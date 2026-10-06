"""Golden facts of the period-3 orbit and the nested tangle.

AUTHOR-GATED: change only with the author's sign-off; there is no re-record
mechanism. These are the CLAUDE.md fixture facts of the higher-period cases
(``cases.build_period3()``, ``cases.build_nested()`` defaults), each pinned
ONCE in the suite (author decision 2, 2026-10-05).

Spelling: letters never appear. A class is its oriented homotopy element pair
``(source, target)``, anchor outward, in FULL names (the orbit-branch code is
the whole point here: ``R_(1.0;5)`` is element 5 on the right of orbit
branch ``1.0``); a word is the list of its tokens' pairs, an inverse token
reversed.

* Period 3 at 13 unstable steps is CLOSED: ``a -> b -> c -> a u^-1 w^-1``, no
  image bridges, and blasting changes nothing. Here ``a, b, c`` are
  ``X0, X1, X2`` below, each running from the anchor element of one orbit
  branch to element 5 of the next, and ``u``, ``w`` the inert classes ``U``,
  ``W``.
* Nested needs two outer blasts for every class to resolve.

Pending author sign-off (E8): the element names of ``X0``, ``X1``, ``X2``,
``U`` and ``W`` (recorded from the 2026-10-05 build; CLAUDE.md states the
words in letters only).
"""

from __future__ import annotations

import pytest

from cases import Case, build_nested, build_period3
from helpers.names import class_by_pair, homotopy_pair, word_in_names

pytestmark = pytest.mark.golden

#: The three active classes, one per orbit branch (E8: names pending sign-off).
X0 = ("R_(0.0;1)", "R_(1.0;5)")
X1 = ("R_(1.0;1)", "R_(2.0;5)")
X2 = ("R_(2.0;1)", "R_(0.0;5)")
#: The two inert classes in ``X2``'s word (E8: names pending sign-off).
U = ("L_(1.0;1)", "L_(1.0;3)")
W = ("R_(1.0;2)", "R_(1.0;4)")

#: ``X0 -> X1 -> X2 -> X0 U^-1 W^-1``.
P3_WORDS = {
    X0: [X1],
    X1: [X2],
    X2: [X0, (U[1], U[0]), (W[1], W[0])],
}


def _active_words(case: Case) -> dict:
    """``{active class pair: word}`` in full names."""
    dyn = case.session.symbolic_dynamics()
    return {
        homotopy_pair(dyn, cd): word_in_names(dyn, cd)
        for cd in dyn.classes.values()
        if cd.kind == "active"
    }


def test_p3_words_form_the_orbit_shift_chain() -> None:
    """Exactly three active classes, mapping one to the next; ``U``, ``W`` inert."""
    case = build_period3()
    assert _active_words(case) == P3_WORDS
    dyn = case.session.symbolic_dynamics()
    for pair in (U, W):
        assert class_by_pair(dyn, pair).kind == "inert"


def test_p3_is_closed() -> None:
    """The minimal trellis maps no hole bridge."""
    assert not build_period3().session.minimal_trellis().image_bridge_ids


def test_p3_words_survive_four_blasts() -> None:
    """Blasting the closed zone four times changes no word."""
    assert _active_words(build_period3(blasts=4)) == P3_WORDS


def test_p3_resolves_and_is_reliable() -> None:
    """Every class resolves (no unreachable walk) and is verified; the dynamics is reliable."""
    dyn = build_period3().session.symbolic_dynamics()
    for cd in dyn.classes.values():
        assert cd.itinerary is not None, cd.unresolved_reason
        assert cd.verified is True
    assert not dyn.virtual_classes
    assert dyn.is_reliable


def test_nested_resolves_everything_after_two_outer_blasts() -> None:
    """Every class resolves, with at least one unique walk and no unreachable one."""
    case = build_nested()
    assert case.expect.blasts == 2
    session = case.session
    dyn = session.symbolic_dynamics()
    for cd in dyn.classes.values():
        assert cd.itinerary is not None, cd.unresolved_reason
    statuses = {cd.search.status for cd in dyn.classes.values() if cd.search is not None}
    assert "unique" in statuses and "unreachable" not in statuses
    assert session.minimal_trellis().image_bridge_ids
    assert dyn.is_reliable
