"""Tests for the symbolic dynamics built from bridge classes and dual-graph walks.

The synthetic tests use a fake naming and hand-built classes over hand-made
``ElementRef`` s; a word spelled here is the rule's output on that hand-made
input. No real build's word, matrix or letter is pinned (author, 2026-10-05).
"""

from __future__ import annotations

import importlib
import logging

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


def test_virtual_class_named_new1_with_warning(layout, caplog):
    itinerary = layout.refs("R_1^1", "R_2", "R_2", "R_3^1")
    with caplog.at_level(logging.WARNING, logger=SD.__name__):
        symbols, _loops = layout.translate(itinerary)
    assert [symbol.kind for symbol in symbols] == ["virtual", "virtual"]
    assert [symbol.text for symbol in symbols] == ["new1", "new2"]
    assert_logged(caplog, logging.WARNING, SD.__name__)
    # the shared letter map remembers the names
    assert layout.letters[BridgeClass(layout.R[0], layout.R[1])] == "new1"
    assert layout.letters[BridgeClass(layout.R[1], layout.R[2])] == "new2"
    # a second translation reuses them and warns no more
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger=SD.__name__):
        again, _ = layout.translate(itinerary)
    assert [symbol.text for symbol in again] == ["new1", "new2"]
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
    assert dyn.virtual_classes == {}


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
    """D2: evidence names a class outside the alphabet locally and never adds a symbol
    to the caller's letters (only words make virtual classes)."""
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
    assert evidence and all(e.symbols[0].kind == "virtual" for e in evidence)


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
