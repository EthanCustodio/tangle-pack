"""Golden facts of the ``k=10`` horseshoe (``cases.build_k10()`` defaults).

AUTHOR-GATED: change only with the author's sign-off; there is no re-record
mechanism. These are the CLAUDE.md fixture facts of the k=10 case, each pinned
ONCE in the suite (author decision 2, 2026-10-05).

Spelling: letters (``a``, ``u``, ``a_1``) are presentation and never appear.
A class is its oriented homotopy element pair ``(source, target)``, anchor
outward; a word is the list of its tokens' pairs, an inverse token reversed;
a refined child is its oriented iterated pair. Names are in the SHORT form
CLAUDE.md writes the single-branch facts in (the branch code is ``0.0``,
checked below). In those terms CLAUDE.md's

* ``a = {R_1, R_3}`` (active) and ``u = {L_1, L_3}`` (inert);
* itinerary ``R_1^1 R_3^3 | L_3 L_1^2 | R_3^1 R_1^3``;
* ``a -> a u^-1 a^-1``;
* ``a_1 = (R_1^1, R_3^3)``, ``a_2 = (R_1^3, R_3^1)``,
  ``a_1 -> a_1 u^-1 a_2^-1`` and ``a_2 -> a_1 u^-1 a_2^-1``;
* refined matrix ``[[1,1],[1,1]]`` over ``(a_1, a_2)``

read as pinned below. Whether the itinerary came from a walk or the trellis
is NOT asserted (decision 2).

Pending author sign-off (E8): :func:`test_k10_inert_class_rests_on_its_folded_loop`.
"""

from __future__ import annotations

import pytest

from cases import Case, build_k10
from helpers.names import (
    brackets,
    class_by_pair,
    homotopy_pair,
    itinerary_pairs,
    matrix_by_classes,
    refined_children,
    refined_words_in_names,
    word_in_names,
)

pytestmark = pytest.mark.golden

#: The active class ``a`` and the inert class ``u``, oriented anchor outward.
X = ("R_1", "R_3")
U = ("L_1", "L_3")
#: The refined children of ``X``, anchor outward.
X1 = ("R_1^1", "R_3^3")
X2 = ("R_1^3", "R_3^1")


def _inverse(pair: tuple[str, str]) -> tuple[str, str]:
    """The pair of an inverse token."""
    return (pair[1], pair[0])


@pytest.fixture
def k10() -> Case:
    """A fresh k=10 build with its symbolic dynamics made."""
    case = build_k10()
    case.session.symbolic_dynamics()
    return case


def test_k10_iterated_rows(k10: Case) -> None:
    """Both iterated rows, anchor outward, with their closedness."""
    dyn = k10.session.symbolic_dynamics()
    naming = dyn.naming
    assert {name.branch_code for name in naming.names} == {"0.0"}
    assert sorted(name.short_text for name in naming.homotopy_names) == [
        "L_1", "L_2", "L_3", "R_1", "R_2", "R_3",
    ]
    rows = {
        side: brackets(naming, result, short=True)
        for (_key, side), result in k10.session.iterated_partition().results.items()
    }
    assert rows == {
        "right": [
            ("R_1^1", "[ ]"), ("R_1^2", "( )"), ("R_1^3", "[ ]"), ("R_2", "( )"),
            ("R_3^1", "[ ]"), ("R_3^2", "( )"), ("R_3^3", "[ ]"),
        ],
        "left": [("L_1^1", "[ )"), ("L_1^2", "[ ]"), ("L_2", "( )"), ("L_3", "[ ]")],
    }


def test_k10_classes(k10: Case) -> None:
    """Exactly two classes: ``X`` active, mapping over itself; ``U`` inert."""
    dyn = k10.session.symbolic_dynamics()
    kinds = {homotopy_pair(dyn, cd, short=True): cd.kind for cd in dyn.classes.values()}
    assert kinds == {X: "active", U: "inert"}
    active = class_by_pair(dyn, X, short=True)
    assert active.bridge_class in active.entry.image_classes


def test_k10_inert_class_rests_on_its_folded_loop(k10: Case) -> None:
    """``U`` is inert on its folded loop alone (E8: pending author sign-off).

    One member has no registered image (no evidence either way); the image
    evidence is exactly the folded loop.
    """
    inert = class_by_pair(k10.session.symbolic_dynamics(), U, short=True).entry
    assert inert.inert
    assert len(inert.unresolved) == 1
    assert not inert.image_classes
    assert inert.loops
    assert {pair for pair, _element in inert.image_loops} == {m.bridge_id for m in inert.loops}


def test_k10_itinerary_and_word(k10: Case) -> None:
    """``X``'s itinerary and word ``X U^-1 X^-1``; ``U`` is a trivial loop."""
    dyn = k10.session.symbolic_dynamics()
    active = class_by_pair(dyn, X, short=True)
    assert not active.ambiguous
    assert itinerary_pairs(dyn, active, short=True) == [
        ("R_1^1", "R_3^3"), ("L_3", "L_1^2"), ("R_3^1", "R_1^3"),
    ]
    assert word_in_names(dyn, active, short=True) == [X, _inverse(U), _inverse(X)]
    assert word_in_names(dyn, class_by_pair(dyn, U, short=True), short=True) == []


def test_k10_refinement(k10: Case) -> None:
    """``X`` splits into ``X1``, ``X2``, both with word ``X1 U^-1 X2^-1``; every member matches."""
    dyn = k10.session.symbolic_dynamics()
    active = class_by_pair(dyn, X, short=True)
    assert refined_children(dyn, active, short=True) == [X1, X2]
    word = [X1, _inverse(U), _inverse(X2)]
    assert refined_words_in_names(dyn, active, short=True) == {X1: word, X2: word}
    assert set(dyn.refined) == {active.bridge_class}
    assert dyn.unmatched_members == {}
    children = dyn.refined[active.bridge_class]
    assert {bid for child in children for bid in child.members} == {
        member.bridge_id for member in active.entry.members if not member.is_loop
    }


def test_k10_transition_matrix(k10: Case) -> None:
    """The refined matrix over ``(X1, X2)`` is all ones; ``U`` is in neither graph nor matrix."""
    dyn = k10.session.symbolic_dynamics()
    assert matrix_by_classes(dyn, short=True) == {
        X1: {X1: 1, X2: 1},
        X2: {X1: 1, X2: 1},
    }
    assert matrix_by_classes(dyn, refined=False, short=True) == {X: {X: 2}}


def test_k10_evidence(k10: Case) -> None:
    """``X``'s registered member images agree with its word; the dynamics is reliable."""
    dyn = k10.session.symbolic_dynamics()
    assert class_by_pair(dyn, X, short=True).verified is True
    assert dyn.is_reliable
