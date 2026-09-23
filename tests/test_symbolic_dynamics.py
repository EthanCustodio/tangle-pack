"""Tests for the symbolic dynamics built from bridge classes and dual-graph walks.

The synthetic tests use a fake naming and hand-built classes over hand-made
``ElementRef`` s; the fixture tests pin the author's k=10 example and the
k=2.8 one-blast singleton path (never by registry ids or bridge counts).
"""

from __future__ import annotations

import importlib
import logging

import numpy as np
import pytest

from tanglepack.topology.BridgeClass import (
    BridgeClass,
    BridgeClassEntry,
    BridgeClassTable,
    BridgeMember,
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
class _FakeFixedPoint:
    period = 1

    def __repr__(self) -> str:
        return "<fp>"


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
        self.fp = _FakeFixedPoint()
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


def test_inert_letter_series():
    series = [SD._inert_letter(i) for i in range(9)]
    assert series == ["u", "v", "w", "x", "y", "z", "uu", "uv", "uw"]
    assert SD._active_letter(0) == "a" and SD._active_letter(26) == "aa"


def test_inert_letters_skip_active_letters(layout):
    fp = layout.fp
    key = layout.key
    inert_a = BridgeClass(ElementRef(key, "left", 0), ElementRef(key, "left", 2))
    inert_b = BridgeClass(ElementRef(key, "left", 1), ElementRef(key, "left", 3))
    active = BridgeClass(ElementRef(key, "right", 0), ElementRef(key, "right", 2))
    table = BridgeClassTable(
        [
            BridgeClassEntry(inert_a, inert=True),
            BridgeClassEntry(inert_b, inert=True),
            BridgeClassEntry(active, inert=False, letter="u"),
        ]
    )
    letters = inert_letters(table)
    assert letters == {inert_a: "v", inert_b: "w"}
    letters = inert_letters(table, used=["v"])
    assert letters == {inert_a: "w", inert_b: "x"}


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


def test_virtual_class_named_new1_with_warning(layout, caplog):
    itinerary = layout.refs("R_1^1", "R_2", "R_2", "R_3^1")
    with caplog.at_level(logging.WARNING, logger=SD.__name__):
        symbols, _loops = layout.translate(itinerary)
    assert [symbol.kind for symbol in symbols] == ["virtual", "virtual"]
    assert [symbol.text for symbol in symbols] == ["new1", "new2"]
    assert any("virtual class" in record.message for record in caplog.records)
    # the shared letter map remembers the names
    assert layout.letters[BridgeClass(layout.R[0], layout.R[1])] == "new1"
    assert layout.letters[BridgeClass(layout.R[1], layout.R[2])] == "new2"
    # a second translation reuses them and warns no more
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger=SD.__name__):
        again, _ = layout.translate(itinerary)
    assert [symbol.text for symbol in again] == ["new1", "new2"]
    assert not any("virtual class" in record.message for record in caplog.records)


def test_cross_side_pair_flagged_with_warning(layout, caplog):
    itinerary = layout.refs("R_1^1", "L_3")
    with caplog.at_level(logging.WARNING, logger=SD.__name__):
        symbols, _loops = layout.translate(itinerary)
    assert len(symbols) == 1 and symbols[0].cross_side
    assert any("different sides" in record.message for record in caplog.records)


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
    assert dyn.virtual_classes == {}


# --------------------------------------------------------------------------- #
# Transition structure
# --------------------------------------------------------------------------- #
def _k10_dynamics(layout) -> SymbolicDynamics:
    # the inert class maps to a trivial loop: (L_3, L_3) gives an empty word
    return layout.dynamics(
        {layout.a: layout.k10_itinerary, layout.u: layout.refs("L_3", "L_3")}
    )


def test_transition_matrix_and_spectral_radius(layout):
    dyn = _k10_dynamics(layout)
    names, matrix = dyn.transition_matrix()
    assert names == ["a_1", "a_2", "u"]
    assert matrix.tolist() == [[1, 1, 1], [1, 1, 1], [0, 0, 0]]
    assert dyn.spectral_radius() == pytest.approx(2.0)
    names, matrix = dyn.transition_matrix(refined=False)
    assert names == ["a", "u"]
    assert matrix.tolist() == [[2, 1], [0, 0]]
    assert dyn.spectral_radius(refined=False) == pytest.approx(2.0)
    assert matrix.dtype == np.int64


def test_transition_graph_inert_symbols_are_sinks(layout):
    networkx = pytest.importorskip("networkx")
    dyn = _k10_dynamics(layout)
    graph = dyn.transition_graph()
    assert isinstance(graph, networkx.DiGraph)
    assert set(graph.nodes) == {"a_1", "a_2", "u"}
    assert graph.out_degree("u") == 0
    assert graph.nodes["u"]["inert"] is True and graph.nodes["u"]["kind"] == "inert"
    assert graph.nodes["a_1"]["parent"] == "a" and graph.nodes["a_1"]["kind"] == "active"
    assert graph.edges["a_1", "a_2"]["weight"] == 1
    assert graph.edges["a_1", "a_2"]["inverse"] == 1
    assert graph.edges["a_1", "u"]["inverse"] == 1
    assert graph.edges["a_1", "a_1"]["inverse"] == 0


def test_inert_class_with_word_warns_but_stays_sink(layout, caplog):
    dyn = layout.dynamics(
        {layout.a: layout.k10_itinerary, layout.u: layout.refs("L_1^1", "L_3")}
    )
    with caplog.at_level(logging.WARNING, logger=SD.__name__):
        names, matrix = dyn.transition_matrix()
    # u now occurs with two distinct pairs across the itineraries, so it refines
    assert names == ["a_1", "a_2", "u_1", "u_2"]
    assert matrix[2].tolist() == [0, 0, 0, 0] and matrix[3].tolist() == [0, 0, 0, 0]
    assert matrix[0].tolist() == [1, 1, 0, 1]
    assert any("sink" in record.message for record in caplog.records)
    assert dyn.spectral_radius() == pytest.approx(2.0)


def test_is_reliable_and_describe(layout):
    dyn = _k10_dynamics(layout)
    assert dyn.is_reliable
    text = dyn.describe()
    assert "a -> a u^-1 a^-1" in text
    assert "R_1^1 R_3^3 | L_3 L_1^2 | R_3^1 R_1^3" in text
    assert "a_1 = (R_1^1, R_3^3)" in text and "a_2 = (R_1^3, R_3^1)" in text
    assert "spectral radius: 2" in text
    assert repr(dyn).startswith("<SymbolicDynamics")
    # an unresolved class makes the result unreliable
    partial = layout.dynamics({layout.a: layout.k10_itinerary})
    assert not partial.is_reliable
    assert [cd.letter for cd in partial.unresolved] == ["u"]
    assert partial.describe()


def test_ambiguous_class_is_unreliable(layout):
    dyn = _k10_dynamics(layout)
    dyn.classes[layout.a].ambiguous = True
    assert not dyn.is_reliable
    dyn.classes[layout.a].ambiguous = False
    dyn.classes[layout.a].verified = False
    assert not dyn.is_reliable


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
    assert any("contributes no evidence" in r.message for r in caplog.records)


# --------------------------------------------------------------------------- #
# Fixture tests: k=10 (function-scoped) and k=2.8 one blast (session, slow)
# --------------------------------------------------------------------------- #
def _dual_and_table(session, fixed_points):
    from minimal_helpers import build_pieces
    from tanglepack.topology.DualGraph import DualGraph

    pieces = build_pieces(session, fixed_points)
    dual = DualGraph(pieces.minimal, pieces.iterated, strong_pips=pieces.strong_pips)
    return pieces, dual


def _names(dyn, refs):
    return [dyn.naming.name(ref).text for ref in refs]


@pytest.fixture
def k10_dynamics(k10_partitioned):
    session, fp = k10_partitioned
    pieces, dual = _dual_and_table(session, [fp])
    return pieces, dual, symbolic_dynamics(dual, pieces.table)


def test_k10_every_landing_contained(k10_dynamics):
    _pieces, _dual, dyn = k10_dynamics
    for cd in dyn.classes.values():
        for landing in cd.landings:
            assert landing is not None and landing.target is not None
            assert landing.contained
    # the anchor-side element lands in element 0
    for cd in dyn.classes.values():
        source, _target = cd.landings
        if source.source.element_id == 0:
            assert source.target.element_id == 0


def test_k10_active_class_word(k10_dynamics):
    _pieces, _dual, dyn = k10_dynamics
    active = [cd for cd in dyn.classes.values() if cd.kind == "active"]
    assert len(active) == 1
    cd = active[0]
    assert cd.source == "walk" and not cd.ambiguous and cd.unresolved_reason is None
    assert _names(dyn, cd.itinerary) == ["R_1^1", "R_3^3", "L_3", "L_1^2", "R_3^1", "R_1^3"]
    letter = cd.letter
    inert = [c.letter for c in dyn.classes.values() if c.kind == "inert"]
    assert inert == ["u"]
    assert cd.word == f"{letter} u^-1 {letter}^-1"
    assert dyn.word(letter, refined=False) == cd.word
    assert cd.verified is True
    assert cd.evidence, "the k=10 anchor class has registered member images"
    assert dyn.is_reliable


def test_k10_refinement_matches_every_member(k10_dynamics):
    _pieces, _dual, dyn = k10_dynamics
    cd = next(c for c in dyn.classes.values() if c.kind == "active")
    letter = cd.letter
    children = dyn.refined[cd.bridge_class]
    assert [child.name for child in children] == [f"{letter}_1", f"{letter}_2"]
    assert _names(dyn, children[0].occurrence) == ["R_1^1", "R_3^3"]
    assert _names(dyn, children[1].occurrence) == ["R_1^3", "R_3^1"]
    assert dyn.word(f"{letter}_1") == f"{letter}_1 u^-1 {letter}_2^-1"
    assert dyn.word(f"{letter}_2") == f"{letter}_1 u^-1 {letter}_2^-1"
    # Under the footprint cut rule every member's own endpoint pair is one of
    # the walked occurrences: nothing is left unmatched.
    assert dyn.unmatched_members == {}
    assert dyn.member_refinement and all(
        child is not None for child in dyn.member_refinement.values()
    )
    assert {bid for child in children for bid in child.members} == {
        member.bridge_id for member in cd.entry.members if not member.is_loop
    }
    assert all(
        bid in dyn.refined[cd.bridge_class][child.index - 1].members
        for bid, child in dyn.member_refinement.items()
    )


def test_k10_itineraries_even_and_matrix(k10_dynamics):
    _pieces, _dual, dyn = k10_dynamics
    for cd in dyn.classes.values():
        assert cd.itinerary is not None, cd.unresolved_reason
        assert len(cd.itinerary) % 2 == 0
        assert not any(symbol.cross_side for symbol in cd.symbols)
    assert dyn.spectral_radius() == pytest.approx(2.0)
    assert dyn.spectral_radius(refined=False) == pytest.approx(2.0)
    graph = dyn.transition_graph()
    assert graph.out_degree("u") == 0
    assert dyn.describe()


@pytest.mark.slow
def test_k28_one_blast_singleton_path(k28_partitioned):
    session, fp = k28_partitioned
    pieces, dual = _dual_and_table(session, [fp])
    dyn = symbolic_dynamics(dual, pieces.table)
    active = [cd for cd in dyn.classes.values() if cd.kind == "active"]
    assert len(active) == 1
    cd = active[0]
    assert cd.unresolved_reason is None
    assert any(landing.singleton for landing in cd.landings)
    assert cd.source == "trellis" and cd.search is None
    assert len(cd.itinerary) == 6
    assert len(cd.symbols) == 3
    assert cd.verified is True
    for other in dyn.classes.values():
        if other.kind == "inert":
            assert other.itinerary is not None, other.unresolved_reason
            assert other.symbols == [], other.word
    for other in dyn.classes.values():
        assert len(other.itinerary) % 2 == 0
    assert dyn.describe()
