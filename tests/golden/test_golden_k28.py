"""Golden facts of the blasted ``k=2.8`` tangle (``cases.build_k28(blasts=1|2)``).

AUTHOR-GATED: change only with the author's sign-off; there is no re-record
mechanism. These are the CLAUDE.md fixture facts of the k=2.8 cases, each
pinned ONCE in the suite (author decision 2, 2026-10-05).

Spelling: letters never appear. A class is its oriented homotopy element pair
``(source, target)``, anchor outward; a word is the list of its tokens' pairs,
an inverse token reversed; a refined child is its oriented iterated pair.
Names are in the SHORT form CLAUDE.md writes the single-branch facts in.

* One blast: the single active class ``{R_1, R_5}`` is singleton to singleton
  with word ``X U^-1 V^-1``, ``U = {L_1, L_3}`` and ``V = {R_2, R_4}`` inert
  trivial loops; the exterior class ``U`` is inert through a VIRTUAL loop (an
  image pair no bridge spans). That the class is resolved "via the trellis
  path" is NOT asserted (decision 2: ``source`` is asserted nowhere); its two
  singleton landings are.
* Two blasts: ``A = {R_1, R_5}`` (the anchor class), ``B = {R_3, R_5}``,
  ``C = {R_1, R_3}``, ``U = {L_1, L_3}``; ``A -> A U^-1 B^-1``, ``B -> C``,
  ``C -> A U^-1 A^-1``; ``A`` refines into ``A1 = (R_1^1, R_5^3)`` and
  ``A2 = (R_1^3, R_5^1)``; the unrefined matrix over ``(A, B, C)`` is
  ``[[1,1,0],[0,0,1],[2,0,0]]``. CLAUDE.md's cut facts
  (``R_1^1=[0,10] R_1^2=(10,7) ...``) name registry ids of one build, which
  are not reproducible (Phase 0 probe P7), so they are pinned as the id-free
  relations P7 derived from them.

Pending author sign-off (E8): :func:`test_k28_one_blast_has_no_image_bridges`.
"""

from __future__ import annotations

import pytest

from cases import Case, build_k28
from helpers.names import (
    brackets,
    class_by_pair,
    homotopy_pair,
    matrix_by_classes,
    refined_children,
    refined_words_in_names,
    word_in_names,
)

pytestmark = pytest.mark.golden


def _inverse(pair: tuple[str, str]) -> tuple[str, str]:
    """The pair of an inverse token."""
    return (pair[1], pair[0])


# --------------------------------------------------------------------------- #
# One blast
# --------------------------------------------------------------------------- #
X = ("R_1", "R_5")
U = ("L_1", "L_3")
V = ("R_2", "R_4")


@pytest.fixture
def k28_one() -> Case:
    """A fresh one-blast k=2.8 build with its symbolic dynamics made."""
    case = build_k28(blasts=1)
    case.session.symbolic_dynamics()
    return case


def test_k28_one_blast_classes_and_word(k28_one: Case) -> None:
    """One active class, singleton to singleton, word ``X U^-1 V^-1``; ``U``, ``V`` trivial loops."""
    dyn = k28_one.session.symbolic_dynamics()
    kinds = {homotopy_pair(dyn, cd, short=True): cd.kind for cd in dyn.classes.values()}
    assert kinds == {X: "active", U: "inert", V: "inert"}
    active = class_by_pair(dyn, X, short=True)
    assert all(landing is not None and landing.singleton for landing in active.landings)
    assert word_in_names(dyn, active, short=True) == [X, _inverse(U), _inverse(V)]
    for pair in (U, V):
        inert = class_by_pair(dyn, pair, short=True)
        assert inert.itinerary is not None and inert.symbols == []
    assert active.verified is True
    assert dyn.is_reliable


def test_k28_one_blast_exterior_class_is_inert_through_a_virtual_loop(k28_one: Case) -> None:
    """``U``'s only image evidence is a loop whose pair no bridge object spans."""
    dyn = k28_one.session.symbolic_dynamics()
    exterior = class_by_pair(dyn, U, short=True).entry
    assert exterior.inert and not exterior.loops and not exterior.image_classes
    trellis = k28_one.session.trellis()
    ((pair, element),) = exterior.image_loops
    assert trellis.bridge_between(*pair) is None
    assert element.side == exterior.bridge_class.source.side


def test_k28_one_blast_has_no_image_bridges(k28_one: Case) -> None:
    """The minimal trellis maps no hole bridge (E8: pending author sign-off)."""
    assert not k28_one.session.minimal_trellis().image_bridge_ids


# --------------------------------------------------------------------------- #
# Two blasts
# --------------------------------------------------------------------------- #
A = ("R_1", "R_5")
B = ("R_3", "R_5")
C = ("R_1", "R_3")
A1 = ("R_1^1", "R_5^3")
A2 = ("R_1^3", "R_5^1")


@pytest.fixture
def k28_two() -> Case:
    """A fresh two-blast k=2.8 build with its symbolic dynamics made."""
    case = build_k28(blasts=2)
    case.session.symbolic_dynamics()
    return case


def test_k28_two_blasts_classes_in_table_order(k28_two: Case) -> None:
    """Three active classes, anchor class first (smallest member cdist), and one inert."""
    session = k28_two.session
    dyn = session.symbolic_dynamics()
    table = session.bridge_classes()
    ordered = [
        (homotopy_pair(dyn, dyn.classes[entry.bridge_class], short=True), entry.inert)
        for entry in table
    ]
    assert [pair for pair, inert in ordered if not inert] == [A, B, C]
    assert [pair for pair, inert in ordered if inert] == [U]
    assert table.entries[0].bridge_class == class_by_pair(dyn, A, short=True).bridge_class
    assert table.entries[0].min_unstable_cdist == pytest.approx(
        0.0, abs=k28_two.registry.cdist_tol
    )


def test_k28_two_blasts_words(k28_two: Case) -> None:
    """``A -> A U^-1 B^-1``, ``B -> C``, ``C -> A U^-1 A^-1``; ``U`` a trivial loop."""
    dyn = k28_two.session.symbolic_dynamics()
    words = {
        homotopy_pair(dyn, cd, short=True): word_in_names(dyn, cd, short=True)
        for cd in dyn.classes.values()
    }
    assert words == {
        A: [A, _inverse(U), _inverse(B)],
        B: [C],
        C: [A, _inverse(U), _inverse(A)],
        U: [],
    }


def test_k28_two_blasts_refinement(k28_two: Case) -> None:
    """Only ``A`` refines; ``A1``, ``A2`` share a word; ``C -> A1 U^-1 A2^-1``; all members match."""
    dyn = k28_two.session.symbolic_dynamics()
    a = class_by_pair(dyn, A, short=True)
    assert set(dyn.refined) == {a.bridge_class}
    assert refined_children(dyn, a, short=True) == [A1, A2]
    a_words = refined_words_in_names(dyn, a, short=True)
    assert set(a_words) == {A1, A2} and a_words[A1] == a_words[A2]
    c = class_by_pair(dyn, C, short=True)
    assert refined_words_in_names(dyn, c, short=True) == {C: [A1, _inverse(U), _inverse(A2)]}
    assert dyn.unmatched_members == {}
    children = dyn.refined[a.bridge_class]
    assert {bid for child in children for bid in child.members} == {
        member.bridge_id for member in a.entry.members if not member.is_loop
    }


def test_k28_two_blasts_transition_matrix(k28_two: Case) -> None:
    """The unrefined matrix over ``(A, B, C)``; ``U`` is absent."""
    dyn = k28_two.session.symbolic_dynamics()
    assert matrix_by_classes(dyn, refined=False, short=True) == {
        A: {A: 1, B: 1, C: 0},
        B: {A: 0, B: 0, C: 1},
        C: {A: 2, B: 0, C: 0},
    }


def test_k28_two_blasts_cut(k28_two: Case) -> None:
    """The iterated cut's closedness and P7's id-free relations (lobe base, chord, fold)."""
    session = k28_two.session
    dyn = session.symbolic_dynamics()
    naming = dyn.naming
    iterated = session.iterated_partition()
    shapes = {
        name: shape
        for result in iterated.results.values()
        for name, shape in brackets(naming, result, short=True)
    }
    assert {name: shapes[name] for name in ("R_1^1", "R_1^2", "R_1^3")} == {
        "R_1^1": "[ ]", "R_1^2": "( )", "R_1^3": "[ ]",
    }
    assert {name: shapes[name] for name in ("R_5^1", "R_5^2", "R_5^3")} == {
        "R_5^1": "[ ]", "R_5^2": "( )", "R_5^3": "[ ]",
    }
    assert {name: shapes[name] for name in ("L_1^1", "L_1^2")} == {
        "L_1^1": "[ )", "L_1^2": "[ ]",
    }

    def element(text: str):
        """The iterated interval printing ``text``."""
        return iterated.element(naming.lookup(text))

    def ends(text: str) -> set[int]:
        """The two bounding crossing ids of an element."""
        interval = element(text)
        return {interval.lo_id, interval.hi_id}

    r11, r12, r13 = element("R_1^1"), element("R_1^2"), element("R_1^3")
    r51, r52, r53 = element("R_5^1"), element("R_5^2"), element("R_5^3")
    l11, l12 = element("L_1^1"), element("L_1^2")
    assert r11.lo_cdist == pytest.approx(0.0, abs=k28_two.registry.cdist_tol)
    assert r11.hi_id == r12.lo_id and r12.hi_id == r13.lo_id
    assert r51.hi_id == r52.lo_id and r52.hi_id == r53.lo_id
    assert (l12.lo_id, l12.hi_id) == (r51.lo_id, r51.hi_id)
    assert l11.hi_id == l12.lo_id == r51.lo_id

    trellis = session.trellis(k28_two.fixed_point)
    holes = {hole.iterate: set(hole.bounding_ids) for hole in trellis.holes}
    assert set(holes) == {0, -1, -2, -3}
    assert ends("R_2") == holes[-1]
    # The hole (-1)'s image lobe: its base is R_1^2, its chord R_5^2 is the
    # reference hole's stretch (also L_2's), and the anchor lobe's fold
    # abuts the hole (-2) at the L_1^1 / L_1^2 boundary.
    assert ends("R_1^2") == {trellis.iterate(iid, 1) for iid in holes[-1]}
    assert ends("R_5^2") == holes[0] == ends("L_2")
    assert l11.hi_id in holes[-2]


def test_k28_two_blasts_evidence(k28_two: Case) -> None:
    """Every class resolves and is verified; the dynamics is reliable."""
    dyn = k28_two.session.symbolic_dynamics()
    for cd in dyn.classes.values():
        assert cd.itinerary is not None, cd.unresolved_reason
        assert cd.verified is True
    assert dyn.is_reliable
