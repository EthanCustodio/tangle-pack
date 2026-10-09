"""Tests for the symbolic dynamics built from bridge classes and dual-graph walks.

The synthetic tests use a fake naming and hand-built classes over hand-made
``ElementRef`` s; a word spelled here is the rule's output on that hand-made
input. No real build's word, matrix or letter is pinned (author, 2026-10-05).
"""

from __future__ import annotations

import importlib
import logging
from types import SimpleNamespace

import numpy as np
import pytest

from helpers.logs import assert_logged
from helpers.fakes import bare_fixed_point
from tanglepack.topology.BridgeClass import (
    BridgeClass,
    BridgeClassEntry,
    BridgeClassTable,
    BridgeMember,
    oriented_class,
)
from tanglepack.topology.TopologyResults import ElementRef

SD = importlib.import_module("tanglepack.topology.SymbolicDynamics")
Symbol = SD.Symbol
ClassDynamics = SD.ClassDynamics
SymbolicDynamics = SD.SymbolicDynamics
inert_letters = SD.inert_letters
translate_itinerary = SD.translate_itinerary
refine = SD.refine
symbolic_dynamics = SD.symbolic_dynamics


# --------------------------------------------------------------------------- #
# Synthetic scaffolding: the k=10 element layout, by hand
# --------------------------------------------------------------------------- #
class _FakeName:
    def __init__(self, text: str) -> None:
        self.text = text

    def __str__(self) -> str:
        return self.text


class _FakeNaming:
    """``parent_of`` / ``name`` / ``describe`` over explicit tables."""

    def __init__(self, parents: dict, names: dict) -> None:
        self._parents = parents
        self._names = names

    def parent_of(self, ref: ElementRef) -> ElementRef:
        return self._parents[ref]

    def name(self, ref: ElementRef) -> _FakeName:
        return _FakeName(self._names[ref])

    def describe(self) -> str:
        return "fake naming: " + " ".join(sorted(self._names.values()))


class _Layout:
    """The k=10 layout: homotopy L_1..L_3 / R_1..R_3, iterated as in the plan."""

    def __init__(self) -> None:
        self.fp = bare_fixed_point()
        self.key = (self.fp, "stable", 0, 0)
        self.fixed_points = [self.fp]
        # homotopy
        self.L = [ElementRef(self.key, "left", i) for i in range(3)]
        self.R = [ElementRef(self.key, "right", i) for i in range(3)]
        # iterated: L_1^1 L_1^2 L_2 L_3 ; R_1^1 R_1^2 R_1^3 R_2 R_3^1 R_3^2 R_3^3
        self.iL = [ElementRef(self.key, "left", i) for i in range(4)]
        self.iR = [ElementRef(self.key, "right", i) for i in range(7)]
        parents = {
            self.iL[0]: self.L[0], self.iL[1]: self.L[0],
            self.iL[2]: self.L[1], self.iL[3]: self.L[2],
            self.iR[0]: self.R[0], self.iR[1]: self.R[0], self.iR[2]: self.R[0],
            self.iR[3]: self.R[1],
            self.iR[4]: self.R[2], self.iR[5]: self.R[2], self.iR[6]: self.R[2],
        }
        names = {
            self.iL[0]: "L_1^1", self.iL[1]: "L_1^2", self.iL[2]: "L_2", self.iL[3]: "L_3",
            self.iR[0]: "R_1^1", self.iR[1]: "R_1^2", self.iR[2]: "R_1^3", self.iR[3]: "R_2",
            self.iR[4]: "R_3^1", self.iR[5]: "R_3^2", self.iR[6]: "R_3^3",
        }
        self.by_name = {text: ref for ref, text in names.items()}
        self.naming = _FakeNaming(parents, names)
        self.a = BridgeClass(self.R[0], self.R[2])      # active a = {R_1, R_3}
        self.u = BridgeClass(self.L[0], self.L[2])      # inert  u = {L_1, L_3}
        self.table = BridgeClassTable(
            [
                BridgeClassEntry(
                    self.u,
                    [BridgeMember((1, 12), +1), BridgeMember((9, 11), +1)],
                    inert=True,
                ),
                BridgeClassEntry(
                    self.a,
                    [BridgeMember((0, 5), +1), BridgeMember((8, 9), -1)],
                    inert=False,
                    letter="a",
                ),
            ]
        )
        self.letters = {self.a: "a", **inert_letters(self.table)}

    def refs(self, *texts: str) -> tuple[ElementRef, ...]:
        return tuple(self.by_name[text] for text in texts)

    @property
    def k10_itinerary(self) -> tuple[ElementRef, ...]:
        return self.refs("R_1^1", "R_3^3", "L_3", "L_1^2", "R_3^1", "R_1^3")

    def translate(self, itinerary):
        return translate_itinerary(
            itinerary, self.naming, self.table, self.letters, self.fixed_points
        )

    def dynamics(self, words: dict) -> SymbolicDynamics:
        """A SymbolicDynamics from ``{class: itinerary}`` over this layout."""
        records = {}
        for entry in self.table:
            cls = entry.bridge_class
            kind = "active" if entry.active else "inert"
            cd = ClassDynamics(entry, self.letters[cls], kind, (None, None))
            itinerary = words.get(cls)
            if itinerary is not None:
                cd.itinerary = tuple(itinerary)
                cd.source = "walk"
                cd.symbols, cd.loops = self.translate(itinerary)
            records[cls] = cd
        return SymbolicDynamics(self.naming, self.table, self.fixed_points, records)


@pytest.fixture
def layout() -> _Layout:
    return _Layout()


# --------------------------------------------------------------------------- #
# Symbols and letters
# --------------------------------------------------------------------------- #
def test_symbol_texts(layout):
    occ = layout.refs("R_1^1", "R_3^3")
    plain = Symbol(layout.a, "a", +1, "active", occ)
    assert plain.text == "a" and plain.base == "a" and plain.mathtext == "$a$"
    inverse = Symbol(layout.a, "a", -1, "active", occ)
    assert inverse.text == "a^-1" and inverse.mathtext == "$a^{-1}$"
    refined = Symbol(layout.a, "a", +1, "active", occ, refined_index=2)
    assert refined.text == "a_2" and refined.base == "a_2" and refined.mathtext == "$a_{2}$"
    refined_inverse = Symbol(layout.a, "a", -1, "active", occ, refined_index=2)
    assert refined_inverse.text == "a_2^-1" and refined_inverse.mathtext == "$a_{2}^{-1}$"
    inert = Symbol(layout.u, "u", -1, "inert", occ)
    assert inert.text == "u^-1" and str(inert) == "u^-1"


def test_inert_letter_series(layout):
    """Inert classes are lettered ``u, v, w, x, y, z, uu, ...`` in table order."""
    key = layout.key
    inert = [
        BridgeClass(ElementRef(key, "left", 0), ElementRef(key, "left", i))
        for i in range(1, 10)
    ]
    table = BridgeClassTable([BridgeClassEntry(cls, inert=True) for cls in inert])
    letters = inert_letters(table)
    assert [letters[cls] for cls in inert] == ["u", "v", "w", "x", "y", "z", "uu", "uv", "uw"]


def test_letter_series_are_disjoint():
    """One bijective numeral over two disjoint digit sets (SIMP-6): active letters use
    ``a..t``, inert ones ``u..z``, so the series never collide and nothing is skipped."""
    letter = SD.letter
    active = [letter(i) for i in range(2000)]
    inert = [letter(i, inert=True) for i in range(2000)]
    assert set("".join(active)) <= set("abcdefghijklmnopqrst")
    assert set("".join(inert)) <= set("uvwxyz")
    assert len(set(active)) == len(active) and len(set(inert)) == len(inert)
    assert letter(0) == "a" and letter(0, inert=True) == "u"
    for series in (active, inert):
        width = sum(len(text) == 1 for text in series)
        assert all(len(text) == 1 for text in series[:width])
        assert series[width] == series[0] * 2  # bijective: no "a0"-style gap
    with pytest.raises(ValueError):
        letter(-1)


# --------------------------------------------------------------------------- #
# Translation
# --------------------------------------------------------------------------- #
def test_k10_itinerary_translates_to_a_uinv_ainv(layout):
    symbols, loops = layout.translate(layout.k10_itinerary)
    assert [symbol.text for symbol in symbols] == ["a", "u^-1", "a^-1"]
    assert [symbol.kind for symbol in symbols] == ["active", "inert", "active"]
    assert loops == []
    # occurrences are class-source first: the inverse tokens are swapped
    assert symbols[0].occurrence == layout.refs("R_1^1", "R_3^3")
    assert symbols[1].occurrence == layout.refs("L_1^2", "L_3")
    assert symbols[2].occurrence == layout.refs("R_1^3", "R_3^1")
    assert not any(symbol.cross_side for symbol in symbols)


def test_odd_itinerary_rejected(layout):
    with pytest.raises(ValueError):
        layout.translate(layout.refs("R_1^1", "R_3^3", "L_3"))


def test_loops_dropped_from_word_but_recorded(layout):
    itinerary = layout.refs("R_1^1", "R_3^3", "L_1^1", "L_1^2", "R_3^1", "R_1^3")
    symbols, loops = layout.translate(itinerary)
    assert [symbol.text for symbol in symbols] == ["a", "a^-1"]
    assert loops == [layout.refs("L_1^1", "L_1^2")]


def test_a_class_outside_the_table_is_named_by_its_label(layout, caplog):
    """Only registered evidence can name a class outside the (closed) table: the token
    is named by the class label, counted active, and nothing is logged or lettered."""
    itinerary = layout.refs("R_1^1", "R_2", "R_2", "R_3^1")
    letters = dict(layout.letters)
    with caplog.at_level(logging.WARNING, logger=SD.__name__):
        symbols, _loops = translate_itinerary(
            itinerary, layout.naming, layout.table, letters, layout.fixed_points
        )
    assert [symbol.kind for symbol in symbols] == ["active", "active"]
    assert [symbol.letter for symbol in symbols] == [
        BridgeClass(layout.R[0], layout.R[1]).label,
        BridgeClass(layout.R[1], layout.R[2]).label,
    ]
    assert letters == layout.letters
    assert_logged(caplog, logging.WARNING, SD.__name__, count=0)


def test_cross_side_pair_flagged_with_warning(layout, caplog):
    itinerary = layout.refs("R_1^1", "L_3")
    with caplog.at_level(logging.WARNING, logger=SD.__name__):
        symbols, _loops = layout.translate(itinerary)
    assert len(symbols) == 1 and symbols[0].cross_side
    assert_logged(caplog, logging.WARNING, SD.__name__)


def test_cross_side_is_flagged_only_within_one_branch(caplog):
    """D4: the left of an inversion point's stable branch 0.0 is the right of 0.1, so a
    pair across two branches is never cross-side; on one branch it is, with a WARNING."""
    fp = bare_fixed_point(inversion=True)
    left_00 = ElementRef((fp, "stable", 0, 0), "left", 0)
    right_01 = ElementRef((fp, "stable", 0, 1), "right", 0)
    right_00 = ElementRef((fp, "stable", 0, 0), "right", 1)
    refs = (left_00, right_01, right_00)
    naming = _FakeNaming({ref: ref for ref in refs}, {ref: ref.label for ref in refs})
    across = oriented_class(left_00, right_01, [fp])[0]
    within = oriented_class(left_00, right_00, [fp])[0]
    table = BridgeClassTable(
        [
            BridgeClassEntry(across, [BridgeMember((0, 1), +1)], inert=False, letter="a"),
            BridgeClassEntry(within, [BridgeMember((2, 3), +1)], inert=False, letter="b"),
        ]
    )
    letters = {across: "a", within: "b"}
    for pair, flagged in (((left_00, right_01), False), ((left_00, right_00), True)):
        caplog.clear()
        with caplog.at_level(logging.WARNING, logger=SD.__name__):
            symbols, _loops = translate_itinerary(pair, naming, table, letters, [fp])
        assert [symbol.cross_side for symbol in symbols] == [flagged]
        assert [symbol.kind for symbol in symbols] == ["active"]
        assert_logged(caplog, logging.WARNING, SD.__name__, count=int(flagged))


# --------------------------------------------------------------------------- #
# Refinement
# --------------------------------------------------------------------------- #
def test_refinement_of_k10_word(layout):
    symbols, _ = layout.translate(layout.k10_itinerary)
    refined = refine({layout.a: symbols}, layout.fixed_points)
    assert set(refined) == {layout.a}
    children = refined[layout.a]
    assert [child.name for child in children] == ["a_1", "a_2"]
    assert children[0].occurrence == layout.refs("R_1^1", "R_3^3")
    assert children[1].occurrence == layout.refs("R_1^3", "R_3^1")
    assert all(child.parent == layout.a and child.letter == "a" for child in children)


def test_single_occurrence_stays_unrefined(layout):
    symbols, _ = layout.translate(layout.refs("R_1^1", "R_3^3", "L_3", "L_1^2"))
    assert refine({layout.a: symbols}, layout.fixed_points) == {}


def test_refined_words_of_k10(layout):
    dyn = layout.dynamics({layout.a: layout.k10_itinerary})
    assert dyn.word("a", refined=False) == "a u^-1 a^-1"
    assert dyn.word("a_1") == "a_1 u^-1 a_2^-1"
    assert dyn.word("a_2") == "a_1 u^-1 a_2^-1"
    assert set(dyn.refined_rules) == {"a_1", "a_2"}
    assert set(dyn.rules) == {"a"}
    with pytest.raises(KeyError):
        dyn.word("u")            # unresolved in this synthetic build: no rule
    with pytest.raises(KeyError):
        dyn.word("b")
    assert dyn.inert_letters == {layout.u: "u"}
    assert all(cd.named_by is None for cd in dyn.classes.values())


# --------------------------------------------------------------------------- #
# Transition structure
# --------------------------------------------------------------------------- #
def _k10_dynamics(layout) -> SymbolicDynamics:
    # the inert class maps to a trivial loop: (L_3, L_3) gives an empty word
    return layout.dynamics(
        {layout.a: layout.k10_itinerary, layout.u: layout.refs("L_3", "L_3")}
    )


def test_transition_matrix(layout):
    # The inert class u is in the words (a -> a u^-1 a^-1) but not in the
    # matrix: inert classes contribute nothing to the dynamics.
    dyn = _k10_dynamics(layout)
    names, matrix = dyn.transition_matrix()
    assert names == ["a_1", "a_2"]
    assert matrix.tolist() == [[1, 1], [1, 1]]
    names, matrix = dyn.transition_matrix(refined=False)
    assert names == ["a"]
    assert matrix.tolist() == [[2]]
    assert matrix.dtype == np.int64
    assert dyn.word("a", refined=False) == "a u^-1 a^-1", "the word still shows u"


def test_transition_graph_leaves_inert_classes_out(layout):
    networkx = pytest.importorskip("networkx")
    dyn = _k10_dynamics(layout)
    graph = dyn.transition_graph()
    assert isinstance(graph, networkx.DiGraph)
    assert set(graph.nodes) == {"a_1", "a_2"}
    assert "u" not in graph
    assert graph.nodes["a_1"]["parent"] == "a" and graph.nodes["a_1"]["kind"] == "active"
    assert graph.edges["a_1", "a_2"]["weight"] == 1
    assert graph.edges["a_1", "a_2"]["inverse"] == 1
    assert graph.edges["a_1", "a_1"]["inverse"] == 0
    assert set(graph.edges) == {("a_1", "a_1"), ("a_1", "a_2"), ("a_2", "a_1"), ("a_2", "a_2")}
    assert set(dyn.transition_graph(refined=False).nodes) == {"a"}


def test_inert_class_mapping_onto_inert_classes_is_silent(layout, caplog):
    """An inert class whose word names only inert symbols is inertness holding
    (on a period-k orbit it maps to the next branch's inert class, PK-8): it
    stays out of the matrix and nothing is logged."""
    dyn = layout.dynamics(
        {layout.a: layout.k10_itinerary, layout.u: layout.refs("L_1^1", "L_3")}
    )
    with caplog.at_level(logging.WARNING, logger=SD.__name__):
        names, matrix = dyn.transition_matrix()
    # u now occurs with two distinct pairs across the itineraries, so it refines
    # (u_1, u_2 exist as refined symbols) but stays out of the matrix.
    assert [child.name for child in dyn.refined[layout.u]] == ["u_1", "u_2"]
    assert names == ["a_1", "a_2"]
    assert matrix[0].tolist() == [1, 1]
    assert_logged(caplog, logging.WARNING, SD.__name__, count=0)


def test_inert_class_naming_an_active_class_warns_and_stays_out(layout, caplog):
    """An inert class whose word names an ACTIVE class contradicts its own
    inertness: one WARNING, and it still stays out of the matrix."""
    dyn = layout.dynamics(
        {layout.a: layout.k10_itinerary, layout.u: layout.refs("R_1^1", "R_3^3")}
    )
    with caplog.at_level(logging.WARNING, logger=SD.__name__):
        names, _matrix = dyn.transition_matrix()
    assert "u" not in names and not any(name.startswith("u_") for name in names)
    assert_logged(caplog, logging.WARNING, SD.__name__, count=1)


def test_describe_lists_every_class(layout):
    dyn = _k10_dynamics(layout)
    text = dyn.describe()
    assert text and repr(dyn)
    for cd in dyn.classes.values():
        assert cd.letter in text
    # a class with no itinerary is listed as unresolved, and still described
    partial = layout.dynamics({layout.a: layout.k10_itinerary})
    assert [cd.letter for cd in partial.unresolved] == ["u"]
    assert partial.describe()


def test_ambiguous_class_is_unreliable(layout):
    dyn = _k10_dynamics(layout)
    cd = dyn.classes[layout.a]
    assert cd.status == "resolved"
    cd.ambiguous = True
    assert cd.status == "ambiguous"


def test_unreadable_image_chain_is_skipped_as_evidence(layout, monkeypatch, caplog):
    """A registered chain the iterated family cannot own (``trellis_itinerary``
    raises ``ValueError``) is dropped from the evidence with a WARNING rather
    than propagating out of ``registered_evidence``."""
    entry = layout.table[layout.a]
    good = layout.k10_itinerary

    class _FakeTrellis:
        def __init__(self, fixed_points):
            self.fixed_points = fixed_points

    def fake_chain(_trellis, bridge_id):
        return [(0, 1), (1, 2), (2, 3)]

    def fake_itinerary(_trellis, _iterated, _chain, direction):
        if direction < 0:
            raise ValueError("crossing 2 of the image chain has no owner")
        return good

    monkeypatch.setattr(SD, "_image_chain", fake_chain)
    monkeypatch.setattr(SD, "trellis_itinerary", fake_itinerary)
    with caplog.at_level(logging.WARNING, logger=SD.__name__):
        evidence = SD.registered_evidence(
            _FakeTrellis(layout.fixed_points),
            object(),
            layout.naming,
            layout.table,
            dict(layout.letters),
            entry,
        )
    assert [e.bridge_id for e in evidence] == [(0, 5)]
    assert evidence[0].itinerary == good
    assert_logged(caplog, logging.WARNING, SD.__name__)


def test_registered_evidence_never_extends_the_letters(layout, monkeypatch):
    """D2: evidence names a class outside the alphabet by its label and never adds a
    symbol to the caller's letters (only a WORD closes the alphabet)."""
    entry = layout.table[layout.a]
    outside = layout.refs("R_1^1", "R_2")  # {R_1, R_2}: no class of the table

    class _FakeTrellis:
        def __init__(self, fixed_points):
            self.fixed_points = fixed_points

    monkeypatch.setattr(SD, "_image_chain", lambda _trellis, _bid: [(0, 1), (1, 2)])
    monkeypatch.setattr(SD, "trellis_itinerary", lambda *_args: outside)
    letters = dict(layout.letters)
    evidence = SD.registered_evidence(
        _FakeTrellis(layout.fixed_points), object(), layout.naming, layout.table, letters, entry
    )
    assert letters == layout.letters
    label = BridgeClass(layout.R[0], layout.R[1]).label
    assert evidence and all(
        (e.symbols[0].kind, e.symbols[0].letter) == ("active", label) for e in evidence
    )
    assert all(e.chain == ((0, 1), (1, 2)) for e in evidence)


def test_a_straddling_landing_is_unreliable(layout):
    """INT-8: a landing across a cut still resolves (the walk exists, D3), but a
    homotopy element's image lies in ONE iterated element, so a straddle counts
    against ``is_reliable``."""
    dyn = _k10_dynamics(layout)
    cd = dyn.classes[layout.a]

    def landing(contained: bool) -> SD.ElementLanding:
        covering = (layout.iR[0],) if contained else (layout.iR[0], layout.iR[1])
        return SD.ElementLanding(
            layout.R[0], layout.key, "right", (0.0, 1.0), (True, True), layout.iR[0],
            contained, False, covering, None,
        )

    cd.landings = (landing(True), landing(True))
    assert dyn.is_reliable
    cd.landings = (landing(True), landing(False))
    assert not dyn.is_reliable


# --------------------------------------------------------------------------- #
# The driver on the synthetic layout: closure and evidence diagnosis
# --------------------------------------------------------------------------- #
class _Crossings:
    """A trellis stand-in: fixed points, ``intersection(id)`` with key and sign, and
    every arc on the first fixed point's first unstable branch."""

    def __init__(self, fixed_points, crossings: dict) -> None:
        self.fixed_points = fixed_points
        self._crossings = crossings

    def intersection(self, crossing_id):
        key, sign = self._crossings[crossing_id]
        return SimpleNamespace(manifold_b_key=key, crossing_sign=sign)

    def unstable_key_of(self, _pair):
        return (self.fixed_points[0], "unstable", 0, 0)


def _drive(layout, monkeypatch, alphabet, evidence, walks, *, crossings=None,
           chains=None, registered=None, letter_for=None):
    """``symbolic_dynamics`` over the layout: every landing is contained, each class
    walks ``walks[class]``, ``chains`` / ``registered`` stand in for the iterate table."""
    def land(_trellis, _homotopy, _iterated, ref):
        return SD.ElementLanding(
            ref, ref.branch_key, ref.side, (0.0, 1.0), (True, True), ref, True, False,
            (ref,), None,
        )

    def resolve(_dual, _iterated, cd, _max_walks, _walls=()):
        cd.itinerary, cd.source = tuple(walks[cd.bridge_class]), "walk"

    def owner_pair(*_args):
        raise ValueError("no owners in the synthetic layout")

    monkeypatch.setattr(SD, "ElementNaming", lambda _homotopy, _iterated: layout.naming)
    monkeypatch.setattr(SD, "land_element", land)
    monkeypatch.setattr(SD, "_resolve_itinerary", resolve)
    monkeypatch.setattr(SD, "_member_pair", owner_pair)
    monkeypatch.setattr(SD, "_image_chain", lambda _t, bid: (chains or {}).get(bid))
    monkeypatch.setattr(SD, "trellis_itinerary", lambda _t, _i, chain, _d: registered[tuple(chain)])
    dual = SimpleNamespace(
        trellis=_Crossings(layout.fixed_points, crossings or {}),
        partition=SimpleNamespace(homotopy=object()),
        minimal=SimpleNamespace(classes=evidence),
    )
    return symbolic_dynamics(dual, alphabet, letter_for=letter_for)


def test_a_word_naming_a_spanned_class_closes_the_alphabet(layout, monkeypatch, caplog):
    """A word naming a class outside the alphabet that a registered bridge spans adds a
    copy of its evidence entry: walked in turn, lettered (an active one by
    ``letter_for``, an inert one in the inert series) and remembering the class that
    named it. No WARNING."""
    alphabet = BridgeClassTable([layout.table[layout.a]])
    walks = {layout.a: layout.k10_itinerary, layout.u: layout.refs("L_3", "L_3")}
    members = layout.table[layout.u].members
    asked = []

    def letter_for(cls):
        asked.append(cls)
        return "q"

    for inert, letter in ((True, "u"), (False, "q")):
        asked.clear()
        evidence = BridgeClassTable(
            [layout.table[layout.a], BridgeClassEntry(layout.u, members, inert=inert)]
        )
        caplog.clear()
        with caplog.at_level(logging.WARNING, logger=SD.__name__):
            dyn = _drive(layout, monkeypatch, alphabet, evidence, walks, letter_for=letter_for)
        assert_logged(caplog, logging.WARNING, SD.__name__, count=0)
        assert dyn.table.classes == [layout.a, layout.u]
        closed = dyn.classes[layout.u]
        assert closed.named_by == layout.a and dyn.classes[layout.a].named_by is None
        assert (closed.letter, closed.kind) == (letter, "inert" if inert else "active")
        assert asked == ([] if inert else [layout.u])
        assert closed.entry.members == members
        assert closed.itinerary == walks[layout.u]
        assert dyn.word("a", refined=False) == f"a {letter}^-1 a^-1"
        assert dyn.is_reliable


def test_a_pair_no_bridge_spans_is_a_flagged_class(layout, monkeypatch, caplog):
    """D11: words need no computed bridge. A pair of a walk that no registered bridge
    spans is a class all the same: added to the table, lettered by ``letter_for``,
    walked from its elements' landings, member-less and flagged ``spanned`` False
    ("no computed bridge yet"). Nothing is unresolved, nothing is warned, and the
    dynamics stays reliable."""
    alphabet = BridgeClassTable([layout.table[layout.a]])
    walks = {layout.a: layout.k10_itinerary, layout.u: layout.refs("L_3", "L_3")}
    asked = []

    def letter_for(cls):
        asked.append(cls)
        return "q"

    with caplog.at_level(logging.WARNING, logger=SD.__name__):
        dyn = _drive(layout, monkeypatch, alphabet, alphabet, walks, letter_for=letter_for)
    assert_logged(caplog, logging.WARNING, SD.__name__, count=0)
    assert dyn.table.classes == [layout.a, layout.u] and asked == [layout.u]
    bare = dyn.classes[layout.u]
    assert bare.entry.members == [] and bare.entry.spanned is False
    assert dyn.classes[layout.a].entry.spanned is True
    assert bare.itinerary == walks[layout.u] and bare.named_by == layout.a
    assert (bare.letter, bare.kind) == ("q", "active")
    assert dyn.unresolved == [] and dyn.is_reliable
    assert "no computed bridge yet" in dyn.describe()


@pytest.mark.parametrize(
    "signs, cause, run",
    [((+1, -1, -1, +1), "stable end", (11, 12)), ((+1, -1, +1, -1), "walk", ())],
)
def test_evidence_disagreement_names_a_stable_end(layout, monkeypatch, caplog, signs, cause, run):
    """INT-1: a registered word that parts from the walk carries a diagnosis: two
    consecutive chain crossings on one stable branch with one sign mean the image
    passed the end of the computed stable branch, and that evidence is INCOMPLETE:
    no verdict, INFO, the dynamics stays reliable (C11, WP10). Otherwise the walk is
    to blame: a WARNING, and the dynamics is not reliable."""
    chain = ((10, 11), (11, 12), (12, 13))
    registered = {chain: layout.refs("R_1^1", "R_3^3", "R_3^1", "R_1^3")}  # a a^-1
    crossings = {cid: (layout.key, sign) for cid, sign in zip((10, 11, 12, 13), signs)}
    walks = {layout.a: layout.k10_itinerary, layout.u: layout.refs("L_3", "L_3")}
    with caplog.at_level(logging.INFO, logger=SD.__name__):
        dyn = _drive(
            layout, monkeypatch, layout.table, layout.table, walks,
            crossings=crossings, chains={(0, 5): list(chain)}, registered=registered,
        )
    cd = dyn.classes[layout.a]
    (evidence,) = cd.evidence
    diagnosis = evidence.diagnosis
    assert (diagnosis.cause, diagnosis.crossings, diagnosis.position) == (cause, run, 1)
    assert diagnosis.walked == layout.refs("L_1^2", "L_3")
    assert diagnosis.registered == layout.refs("R_1^3", "R_3^1")
    incomplete = cause == "stable end"
    assert evidence.word_agrees is (None if incomplete else False)
    assert cd.verified is (None if incomplete else False)
    assert_logged(caplog, logging.WARNING, SD.__name__, count=0 if incomplete else 1)
    assert any("registered image" in r.getMessage() for r in caplog.records if r.levelno == logging.INFO) == incomplete
    assert dyn.is_reliable == incomplete


# --------------------------------------------------------------------------- #
# Real data: fallback letters and an unresolved class
# --------------------------------------------------------------------------- #
def test_unlettered_active_classes_get_fallback_letters(k10_partitioned):
    """A table built outside the session carries no letters: its active classes
    are lettered ``a, b, ...`` in table order, its inert ones ``u, v, ...``."""
    from tanglepack.topology import letter
    from tanglepack.topology.BridgeClass import bridge_classes

    session, _fp = k10_partitioned
    table = bridge_classes(session.trellis(), session.homotopy_partition().as_list())
    assert table.active and all(entry.letter is None for entry in table)

    dyn = symbolic_dynamics(session.dual_graph(), table)

    active = [dyn.classes[entry.bridge_class].letter for entry in table.active]
    assert active == [letter(i) for i in range(len(active))]
    inert = [dyn.classes[entry.bridge_class].letter for entry in table.inert]
    assert set(inert) == set(dyn.inert_letters.values())


def test_unresolved_class_carries_a_reason_and_warns(caplog):
    """The unblasted nested tangle is not grown far enough for every class: a
    class whose landing is a singleton and none of whose members has a
    registered image chain stays unresolved, says why, and warns; every other
    class still reads an even itinerary."""
    from cases import build_nested

    case = build_nested(outer_blasts=0, p1_area_cutoff=1e-4)
    with caplog.at_level(logging.WARNING, logger=SD.__name__):
        dyn = case.session.symbolic_dynamics()

    unresolved = [cd for cd in dyn.classes.values() if cd.itinerary is None]
    assert unresolved, "the unblasted nested tangle has an unresolved class"
    assert all(cd.unresolved_reason for cd in unresolved)
    assert any(landing.singleton for cd in unresolved for landing in cd.landings if landing)
    assert_logged(caplog, logging.WARNING, SD.__name__)
    for cd in dyn.classes.values():
        if cd.itinerary is not None:
            assert len(cd.itinerary) % 2 == 0


# --------------------------------------------------------------------------- #
# The unstable branches a class's arcs lie on (the walls, D13 / D14)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("period, inversion", [(3, False), (1, True), (2, True)])
def test_a_memberless_class_lies_one_step_past_its_namer(period, inversion):
    """A class no computed bridge spans yet (D11) is a piece of ``f(x)`` for ``x`` in
    the class that named it, so its arcs lie one map step past the namer's."""
    fp = bare_fixed_point(period, inversion=inversion)
    namer = {(fp, "unstable", 0, 0)}
    cls = BridgeClass(ElementRef((fp, "stable", 0, 0), "left", 0), ElementRef((fp, "stable", 0, 0), "left", 1))
    branches = SD._arc_branches(None, BridgeClassEntry(cls), frozenset(namer))
    assert branches == {fp.advance_key(key, 1) for key in namer}


def test_a_class_with_members_lies_on_its_members_branches():
    fp = bare_fixed_point(3)
    on = {(1, 2): (fp, "unstable", 2, 0), (3, 4): (fp, "unstable", 0, 0)}
    trellis = SimpleNamespace(unstable_key_of=on.__getitem__)
    cls = BridgeClass(ElementRef((fp, "stable", 0, 0), "left", 0), ElementRef((fp, "stable", 0, 0), "left", 1))
    entry = BridgeClassEntry(cls, [BridgeMember((1, 2), +1), BridgeMember((3, 4), -1)])
    assert SD._arc_branches(trellis, entry, frozenset({(fp, "unstable", 1, 0)})) == set(on.values())
