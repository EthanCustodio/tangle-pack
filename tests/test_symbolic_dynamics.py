"""Phase D — walks, words, and the symbolic dynamics of a trellis.

The claim under test is a geometric one: the word this layer READS off a finite
picture is the word the map actually writes. Everything else here supports that.

**Geometric validation** (the important one). Every bridge's image is walked on
the trellis as it stands, then the unstable manifolds are grown one more step
with the registry ids preserved, and the image bridges that growth registers are
classed AGAINST THE ORIGINAL PARTITIONS. Those classes must be the walked word,
letter for letter. Classing the images against the original partitions rather
than against the grown session's own is what makes the comparison meaningful:
one more growth step re-partitions the stable branches (on k=10 the left side
goes from three elements to five), so element ids — and even an element's
``(lo_id, hi_id)`` — do not survive, while the original partition still owns
every crossing the grown picture puts on those branches.

How much of that is independent evidence has to be said plainly. For a case "i"
bridge it is near-tautological: the word is the classes of the bridges the
iterate table already names as the image, and the validation reads the SAME
table back through ``image_bridges``. What it does check there is that the row,
the element lookup and the ordering survive a re-partition. The independent
evidence is the case "ii" and "iii" bridges, whose words come from a walk the
table knows nothing about — and, on k=10, the fact that the two right-side
classes each hold one case "i" member and one case "iii" member and the two
spell the same word (that is D.4, and it is a real cross-check of the walk
against the table).

**The three cases.** A bridge whose two endpoint iterates are registered has its
image TILED by registered bridges (case "i", nothing walked); one with only the
first registered gets a registered prefix and a walked tail (case "ii"); one with
neither is walked end to end (case "iii"). Both fixtures exercise all three.

**New symbols.** A letter no computed bridge realises means the trellis is not
closed under the map at this extent. k=10 has none; the nested period-3 fixture
has exactly four, all of them the innermost self-class of a branch, and they are
pinned rather than asserted away.

**Uniqueness.** A face with two filled walls onto the same neighbour gives two
spellings of one image, which is an ``AmbiguousWalkError`` — built here by hand,
because neither fixture is tangled enough to produce one.
"""

from __future__ import annotations

import networkx as nx
import numpy as np
import pytest

from tanglepack import TangleSession
from tanglepack.examples import (
    HENON_K10,
    henon_map as _henon_map_factory,
    henon_map_inverse as _henon_map_inverse_factory,
)
from tanglepack.topology.BridgeClass import BridgeClass
from tanglepack.topology.DualGraph import (
    AmbiguousWalkError,
    ArcNode,
    DualGraph,
    FaceNode,
)
from tanglepack.topology.StablePartition import (
    owns_cdist,
    row_of_end,
    rows_of_bridge,
)
from tanglepack.topology.SymbolicDynamics import (
    SymbolicDynamics,
    symbol_label,
    symbolic_dynamics,
    word_of_bridge,
)
from tanglepack.topology.TopologyResults import Arc, ElementRef, OPPOSITE_SIDE


# --------------------------------------------------------------------------- #
# Builders
# --------------------------------------------------------------------------- #
def _k10_graph(session, fp) -> DualGraph:
    """The dual graph of the k=10 fixture, with its one strong pip."""
    trellis = session.trellis(fp)
    return DualGraph(
        session.arrangement(),
        trellis.stable_partitions,
        strong_pips=[trellis.strong_pip],
    )


def _p3_graph(session, fp3, fp1) -> DualGraph:
    """The dual graph of the nested fixture: both tangles' partitions and pips."""
    t3, t1 = session.trellis(fp3), session.trellis(fp1)
    return DualGraph(
        session.arrangement(),
        list(t3.stable_partitions) + list(t1.stable_partitions),
        strong_pips=[t3.strong_pip, t1.strong_pip],
    )


def _partition_index(partitions) -> dict:
    """Index stable partition results by ``(branch_key, side)``."""
    return {(result.branch_key, result.side): result for result in partitions}


def _fresh_k10_symbolic_dynamics() -> SymbolicDynamics:
    """Build a k=10 session from scratch and read its symbolic dynamics.

    Mirrors the ``k10_session`` / ``k10_partitioned`` fixtures. A second,
    independent build is what the determinism test compares against, and a
    function-scoped fixture can only be handed out once per test.
    """
    session = TangleSession(
        _henon_map_factory(*HENON_K10), _henon_map_inverse_factory(*HENON_K10)
    )
    fp = session.construct_fixed_point([4, -4])
    session.orient_eigenvectors(
        fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
    )
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=9)
    session.grow_until_turnaround(fp, "stable")
    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    session.infer_iterate_table()
    session.classify_strong_pips()
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()
    return symbolic_dynamics(_k10_graph(session, fp), session.bridge_classes())


# --------------------------------------------------------------------------- #
# Geometric validation: the walked word against what growth registers
# --------------------------------------------------------------------------- #
def _element_against(trellis, partitions, bridge_id, endpoint, tol) -> ElementRef:
    """The element of a FIXED partition that a bridge end sits against.

    Deliberately not ``element_of_intersection``: the crossing may be one growth
    has just added, which the original partition has never seen. Its stable
    canonical distance is all that is needed — the partition's intervals tile
    the branch, so exactly one of them owns it.
    """
    row = row_of_end(trellis, bridge_id, endpoint)
    intersection_id = bridge_id[0 if endpoint == "first" else 1]
    crossing = trellis.intersection(intersection_id)
    result = partitions[(crossing.manifold_b_key, row)]
    owners = [
        interval
        for interval in result.intervals
        if owns_cdist(interval, float(crossing.stable_cdist), tol)
    ]
    assert len(owners) == 1, (
        f"{len(owners)} elements of the {row!r} partition of "
        f"{crossing.manifold_b_key[1:]} own crossing {intersection_id} at "
        f"stable cdist {crossing.stable_cdist}"
    )
    return result.ref(owners[0].element_id)


def _class_against(trellis, partitions, bridge_id, tol) -> BridgeClass:
    """The class of a bridge measured against a partition it did not build."""
    return BridgeClass(
        _element_against(trellis, partitions, bridge_id, "first", tol),
        _element_against(trellis, partitions, bridge_id, "second", tol),
    )


def _grow_one_unstable_step(session, fixed_points, *, repartition: bool) -> None:
    """One more unstable step, ids preserved, bridges re-cut, table re-inferred."""
    for fp in fixed_points:
        session.grow_n_times(fp, "unstable", num_iterations=1)
    session.compute_intersections(list(fixed_points), preserve_ids=True)
    for fp in fixed_points:
        session.trim_stable_manifolds(fp)
        session.create_bridges(fp)
    session.infer_iterate_table()
    if repartition:
        session.classify_strong_pips()
        session.compute_pseudoneighbors()
        session.punch_holes()
        session.partition_stable_manifold()


def _coarsen_ref(ref, fine, coarse, tol) -> ElementRef:
    """The element of the ORIGINAL partition containing a finer element.

    One growth step splits elements rather than moving their boundaries, so a
    finer element lies inside exactly one coarse one — the one owning its
    midpoint canonical distance.
    """
    interval = fine[(ref.branch_key, ref.side)].element(ref.element_id)
    middle = 0.5 * (interval.lo_cdist + interval.hi_cdist)
    result = coarse[(ref.branch_key, ref.side)]
    owners = [
        candidate
        for candidate in result.intervals
        if owns_cdist(candidate, middle, tol)
    ]
    assert len(owners) == 1, (
        f"{len(owners)} coarse elements contain {ref.label} (midpoint "
        f"{middle}); one growth step should only ever refine a partition"
    )
    return result.ref(owners[0].element_id)


def _coarsen_class(bridge_class, fine, coarse, tol) -> BridgeClass:
    """A class of the finer partition, read on the coarser one."""
    return BridgeClass(
        _coarsen_ref(bridge_class.first, fine, coarse, tol),
        _coarsen_ref(bridge_class.second, fine, coarse, tol),
    )


def _validate_against_growth(session, fixed_points, dual_graph, symbolics):
    """Grow one more unstable step and compare every walked word with reality.

    Returns ``(validated, case counts)``; raises through the assertions on any
    letter that disagrees. See the module docstring for which of the three cases
    this is independent evidence for (the "ii" and "iii" ones) and which it is
    near-tautological for (case "i", where both sides read one iterate table).
    """
    partitions = _partition_index(dual_graph.partitions.values())
    tol = session.trellis().registry.cdist_tol

    _grow_one_unstable_step(session, fixed_points, repartition=False)
    grown = session.trellis()

    validated = 0
    cases = {"i": 0, "ii": 0, "iii": 0}
    for bridge_id, walk in symbolics.walks.items():
        images = session.image_bridges(bridge_id)
        assert images is not None, (
            f"bridge {bridge_id} still has no registered image after one more "
            "growth step; the validation has nothing to compare against"
        )
        actual = tuple(
            _class_against(grown, partitions, image, tol) for image in images
        )
        assert actual == walk.word, (
            f"bridge {bridge_id} (case {walk.case}): walked "
            f"{[letter.label for letter in walk.word]} but growth registered "
            f"{[letter.label for letter in actual]}"
        )
        cases[walk.case] += 1
        validated += 1
    return validated, cases


def test_k10_words_match_the_registered_images(k10_partitioned):
    """Every k=10 bridge's walked word is the word growth actually writes.

    Three of the seven bridges are case "i", where this is near-tautological.
    The load is carried by the other four — one case "ii" and three case "iii",
    all walked — and by the two right-side classes, each of which holds a case
    "i" bridge and a case "iii" bridge that spell the same word.
    """
    session, fp = k10_partitioned
    dual = _k10_graph(session, fp)
    symbolics = symbolic_dynamics(dual, session.bridge_classes())

    validated, cases = _validate_against_growth(session, [fp], dual, symbolics)

    assert validated == 7, "k=10 carries seven non-partial bridges"
    # Three images are fully registered, one has a registered prefix only, and
    # three are walked end to end -- all three cases are exercised.
    assert cases == {"i": 3, "ii": 1, "iii": 3}


def test_k10_survives_another_unstable_step(k10_partitioned):
    """Ten unstable steps still read cleanly — the regression that caught D.6.

    At ten steps the image crossings sit near unstable cdist 6000 with gaps of
    20 to 40 between them, so a window of 1% of the absolute canonical distance
    (which is what the scaling law's relative error would need) sweeps in the
    neighbouring crossing and the image span comes out wrong. The span's two
    ends are known by ID from the iterate table whenever they are known at all,
    so it is sliced out of ``ordered_ids()`` by index and no window exists to be
    mis-sized. Nine steps hid the bug; ten does not.
    """
    session, fp = k10_partitioned
    trellis = session.trellis(fp)
    coarse = _partition_index(trellis.stable_partitions)
    tol = trellis.registry.cdist_tol
    before = symbolic_dynamics(
        _k10_graph(session, fp), session.bridge_classes()
    )

    _grow_one_unstable_step(session, [fp], repartition=True)
    grown_trellis = session.trellis(fp)
    after = symbolic_dynamics(
        _k10_graph(session, fp), session.bridge_classes()
    )

    # 15 bridges now: seven images fully registered, one prefix-only, seven
    # walked; the left partition has refined from three elements to five, which
    # splits the old five symbols into seven.
    assert after.case_counts() == {"i": 7, "ii": 1, "iii": 7}
    assert len(after.symbols) == 7
    assert after.new_symbols == []
    assert after.rules == {
        "a": "a",
        "b": "a",
        "c": "a",
        "d": "b",
        "e": "c",
        "f": "feg",
        "g": "fdg",
    }

    # Every one of the finer symbols lies inside an original symbol, and read
    # at the original resolution its word is the original word.
    fine = _partition_index(grown_trellis.stable_partitions)
    compared = 0
    for bridge_class, word in after.words.items():
        original = before.words.get(
            _coarsen_class(bridge_class, fine, coarse, tol)
        )
        if original is None:
            continue
        assert (
            tuple(_coarsen_class(letter, fine, coarse, tol) for letter in word)
            == original
        ), f"{bridge_class.label} disagrees with its coarser self"
        compared += 1
    assert compared == len(after.words) == 7


@pytest.mark.slow
def test_p3_words_match_the_registered_images(p3_partitioned):
    """Every bridge of the nested fixture validates against one growth step.

    Twenty-six of the thirty are case "i" (near-tautological, see the module
    docstring); the four that are not are the ones this really tests.
    """
    session, fp3, fp1 = p3_partitioned
    dual = _p3_graph(session, fp3, fp1)
    symbolics = symbolic_dynamics(dual, session.bridge_classes())

    validated, cases = _validate_against_growth(
        session, [fp3, fp1], dual, symbolics
    )

    assert validated == 30, "the nested fixture carries thirty bridges"
    assert cases == {"i": 26, "ii": 2, "iii": 2}


# --------------------------------------------------------------------------- #
# The alphabet and its words
# --------------------------------------------------------------------------- #
def test_k10_alphabet_and_rules(k10_partitioned):
    """The k=10 tangle is a five-symbol substitution with no new symbols."""
    session, fp = k10_partitioned
    symbolics = symbolic_dynamics(
        _k10_graph(session, fp), session.bridge_classes()
    )

    assert list(symbolics.symbols) == ["a", "b", "c", "d", "e"]
    assert symbolics.new_symbols == []
    # The two-sided horseshoe: each outer symbol rewrites to three letters, one
    # of them itself, and the inner ones collapse onto the anchor element.
    assert symbolics.rules == {
        "a": "a",
        "b": "a",
        "c": "a",
        "d": "dce",
        "e": "dbe",
    }
    assert symbolics.case_counts() == {"i": 3, "ii": 1, "iii": 3}


def test_no_word_is_empty(k10_partitioned):
    """A bridge's image always crosses at least one face."""
    session, fp = k10_partitioned
    symbolics = symbolic_dynamics(
        _k10_graph(session, fp), session.bridge_classes()
    )
    assert symbolics.words
    for bridge_class, word in symbolics.words.items():
        assert word, f"{bridge_class.label} has an empty word"


def test_k10_every_letter_is_a_known_class(k10_partitioned):
    """k=10 is grown far enough to be closed under the map at this extent."""
    session, fp = k10_partitioned
    classes = session.bridge_classes()
    symbolics = symbolic_dynamics(_k10_graph(session, fp), classes)
    for word in symbolics.words.values():
        for letter in word:
            assert letter in classes, f"{letter.label} is a new symbol on k=10"


def test_anchor_class_word_begins_with_itself(k10_partitioned):
    """The image of the anchor bridge starts at the anchor and contains it."""
    session, fp = k10_partitioned
    classes = session.bridge_classes()
    symbolics = symbolic_dynamics(_k10_graph(session, fp), classes)

    trellis = session.trellis()
    anchors = [
        bridge_class
        for bridge_class, members in classes.items()
        if any(
            trellis.intersection(bridge_id[0]).is_synthetic
            for bridge_id in members
        )
    ]
    assert anchors, "k=10 has an anchor bridge, whose first crossing is synthetic"
    for bridge_class in anchors:
        word = symbolics.word_of(bridge_class)
        assert word[0] == bridge_class, (
            f"the anchor class {bridge_class.label} maps to "
            f"{[letter.label for letter in word]}, which does not begin with "
            "itself"
        )


def test_transition_matrix_rows_sum_to_word_lengths(k10_partitioned):
    """Each row of the substitution matrix counts the letters of one word."""
    session, fp = k10_partitioned
    symbolics = symbolic_dynamics(
        _k10_graph(session, fp), session.bridge_classes()
    )
    matrix = symbolics.transition_matrix()
    assert matrix.shape == (len(symbolics.symbols), len(symbolics.symbols))
    for row, bridge_class in enumerate(symbolics.symbols.values()):
        assert int(matrix[row].sum()) == len(symbolics.words.get(bridge_class, ()))


def test_word_of_bridge_matches_the_walk(k10_partitioned):
    """The convenience front returns exactly the walk's word."""
    session, fp = k10_partitioned
    classes = session.bridge_classes()
    dual = _k10_graph(session, fp)
    for members in classes.values():
        for bridge_id in members:
            assert word_of_bridge(dual, bridge_id, classes) == dual.walk_bridge(
                bridge_id, classes
            ).word


def test_symbol_labels_run_past_the_alphabet():
    """Labels are a..z, then a1..z1, then a2..."""
    assert [symbol_label(i) for i in (0, 25, 26, 51, 52)] == [
        "a",
        "z",
        "a1",
        "z1",
        "a2",
    ]
    with pytest.raises(ValueError, match="non-negative"):
        symbol_label(-1)


def test_describe_is_deterministic_across_builds(k10_partitioned):
    """Two independent builds of k=10 print byte-identical reports."""
    session, fp = k10_partitioned
    first = symbolic_dynamics(
        _k10_graph(session, fp), session.bridge_classes()
    ).describe()
    second = _fresh_k10_symbolic_dynamics().describe()
    assert first == second
    # The report names the fixed point by POSITION, never by address or period.
    assert "fp0@0.0/L#0 -> fp0@0.0/L#0" in first


@pytest.mark.slow
def test_p3_new_symbols_are_the_innermost_self_classes(p3_partitioned):
    """The nested fixture is NOT closed under the map: four letters are new.

    Each is the class of a bridge lying entirely inside a branch's innermost
    element -- exactly the bridge k=10 has (its class ``a``) and the period-3
    fixture has not grown far enough to cut.
    """
    session, fp3, fp1 = p3_partitioned
    classes = session.bridge_classes()
    symbolics = symbolic_dynamics(_p3_graph(session, fp3, fp1), classes)

    assert len(classes) == 30
    assert len(symbolics.new_symbols) == 4
    for symbol in symbolics.new_symbols:
        assert symbol not in classes
        assert symbol.first.side == symbol.second.side
        assert symbol.first.branch_key == symbol.second.branch_key
        assert symbol.first.element_id == symbol.second.element_id
        # The innermost element of a side: id 0 on the left, id 1 on the right
        # (the right partition's element 0 is the pinched anchor singleton).
        assert symbol.first.element_id in (0, 1)
        with pytest.raises(KeyError, match="new symbol"):
            symbolics.word_of(symbol)


# --------------------------------------------------------------------------- #
# The transition graph inside the resonance zone
# --------------------------------------------------------------------------- #
def test_k10_zone_classes_are_strongly_connected(k10_partitioned):
    """Inside the resonance zone the substitution mixes every symbol into every
    other.

    Outside it the symbols drain into the anchor class (``a -> a`` is a sink),
    so strong connectivity is a statement about the zone, not the whole tangle.
    """
    session, fp = k10_partitioned
    trellis = session.trellis(fp)
    session.resonance_zone(trellis.strong_pip)
    # Defining a zone trims and recomputes, which drops the cached trellis and
    # its partitions; rebuild them before reading the symbols off again.
    session.classify_strong_pips()
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()

    classes = session.bridge_classes()
    symbolics = symbolic_dynamics(_k10_graph(session, fp), classes)
    inside = {
        symbolics.label(bridge_class)
        for bridge_class, members in classes.items()
        if any(
            session.classify_bridge(session.bridge(bridge_id)) is not None
            for bridge_id in members
        )
    }

    subgraph = symbolics.transition_graph.subgraph(inside)
    # The zone holds the two outer classes d and e (the four bridges of the
    # resonance zone's own lobes); each rewrites to a word containing both.
    assert subgraph.number_of_nodes() == 2
    assert subgraph.number_of_edges() == 4
    assert nx.is_strongly_connected(subgraph)


# --------------------------------------------------------------------------- #
# Ambiguity, on a graph built by hand
# --------------------------------------------------------------------------- #
def _hand_arc(branch_key, lo_id, hi_id, lo, hi, element_id, filled) -> ArcNode:
    """One synthetic stable arc node carrying a distinct element on each side."""
    return ArcNode(
        arc=Arc(
            kind="stable",
            lo_id=lo_id,
            hi_id=hi_id,
            branch_key=branch_key,
            bridge_id=None,
        ),
        key=("stable", lo_id, hi_id, branch_key),
        branch_key=branch_key,
        lo_cdist=lo,
        hi_cdist=hi,
        left=ElementRef(branch_key, "left", element_id),
        right=ElementRef(branch_key, "right", element_id),
        filled=filled,
    )


def _walked_bridge(session, fp):
    """A k=10 bridge whose image is walked end to end, and where it lands.

    Registry ids are assigned in detection order and that order is not
    reproducible between builds (the crossing detector keys its segment edges on
    object identity), so no test may name a crossing by number. This picks a
    case "iii" bridge by its walk and reports everything the hand-built graphs
    below need: the two sides the image approaches its ends from, and the two
    stable canonical distances it lands at.
    """
    trellis = session.trellis(fp)
    classes = session.bridge_classes()
    dual = _k10_graph(session, fp)
    for members in classes.values():
        for bridge_id in members:
            walk = dual.walk_bridge(bridge_id, classes)
            if walk.case != "iii":
                continue
            rows = rows_of_bridge(trellis, bridge_id)
            ends = tuple(
                trellis.image_cdist(intersection_id, 1, "stable")[1]
                for intersection_id in bridge_id
            )
            if abs(ends[0] - ends[1]) > 1e-6:
                return bridge_id, rows, ends
    raise AssertionError("the k=10 fixture has a walked bridge with distinct ends")


def _hand_built_graph(session, fp, *, routes: int):
    """A synthetic dual graph with ``routes`` filled walls between two faces.

    The real trellis stays underneath — the walk still resolves the bridge's
    case and its two image canonical distances from the registry — but the arcs
    and faces it crosses are these, so the number of spellings is chosen here.
    """
    trellis = session.trellis(fp)
    bridge_id, rows, ends = _walked_bridge(session, fp)
    branch_key = trellis.stable_branches[0].key
    half = min(abs(ends[0] - ends[1]) / 4.0, 0.1)
    far_away = max(ends) + 10.0

    start = _hand_arc(
        branch_key, 900, 901, ends[0] - half, ends[0] + half, 0, filled=False
    )
    end = _hand_arc(
        branch_key, 902, 903, ends[1] - half, ends[1] + half, 1, filled=False
    )
    near = FaceNode(index=0, faces=[], kind="region")
    far = FaceNode(index=1, faces=[], kind="region")
    outside = FaceNode(index=2, faces=[], kind="outer", is_unbounded=True)

    other_start = OPPOSITE_SIDE[rows[0]]
    other_end = OPPOSITE_SIDE[rows[1]]
    start.faces = {rows[0]: near, other_start: outside}
    end.faces = {rows[1]: far, other_end: outside}
    near.arcs = [(start, rows[0])]
    far.arcs = [(end, rows[1])]
    outside.arcs = [(start, other_start), (end, other_end)]

    arcs = [start, end]
    for index in range(routes):
        wall = _hand_arc(
            branch_key,
            904 + 2 * index,
            905 + 2 * index,
            far_away + 2 * index,
            far_away + 2 * index + 1,
            2 + index,
            filled=True,
        )
        wall.faces = {"left": near, "right": far}
        near.arcs.append((wall, "left"))
        far.arcs.append((wall, "right"))
        arcs.append(wall)

    dual = DualGraph._from_parts(trellis, arcs, [near, far, outside])
    return dual, bridge_id


def test_two_filled_walls_on_one_face_are_ambiguous(k10_partitioned):
    """Two routes spelling two words raise rather than picking one.

    Neither fixture has a face with two filled walls onto the same neighbour, so
    the configuration is built by hand around the real k=10 trellis.
    """
    session, fp = k10_partitioned
    dual, bridge_id = _hand_built_graph(session, fp, routes=2)

    with pytest.raises(AmbiguousWalkError) as excinfo:
        dual.walk_bridge(bridge_id, {})
    error = excinfo.value
    assert error.bridge_id == bridge_id
    assert len(error.words) == 2
    assert len(error.example_paths) == 2
    # The two words differ only in which wall was crossed (elements 2 and 3).
    assert {word[0].second.element_id for word in error.words} == {2, 3}
    assert all(len(word) == 2 for word in error.words)
    assert "different ways" in str(error)


def test_one_filled_wall_spells_one_word(k10_partitioned):
    """The same hand-built graph with a single wall is unambiguous."""
    session, fp = k10_partitioned
    dual, bridge_id = _hand_built_graph(session, fp, routes=1)

    walk = dual.walk_bridge(bridge_id, {})
    assert walk.case == "iii"
    assert walk.prefix_len == 0
    assert [letter.second.element_id for letter in walk.word] == [2, 1]


def test_an_unreachable_end_is_an_error(k10_partitioned):
    """A wall between the two ends of an image is a named failure, not silence."""
    session, fp = k10_partitioned
    dual, bridge_id = _hand_built_graph(session, fp, routes=0)
    with pytest.raises(ValueError, match="no walk of the dual graph"):
        dual.walk_bridge(bridge_id, {})


# --------------------------------------------------------------------------- #
# The stub lookup the case "ii" walk starts from
# --------------------------------------------------------------------------- #
def test_face_at_stub_finds_the_open_face(k10_partitioned):
    """The outermost crossing's outward unstable ray sticks into an open face."""
    session, fp = k10_partitioned
    trellis = session.trellis(fp)
    arrangement = session.arrangement()
    last_id = trellis.unstable_branches[0].ordered_ids()[-1]

    face = arrangement.face_at_stub(last_id, "u+")
    assert not face.is_closed
    assert face in arrangement.open_faces

    with pytest.raises(ValueError, match="not a dangling stub"):
        arrangement.face_at_stub(last_id, "u-")
    with pytest.raises(ValueError, match="slot must be one of"):
        arrangement.face_at_stub(last_id, "x+")
    with pytest.raises(ValueError, match="not a node of this arrangement"):
        arrangement.face_at_stub(10_000, "u+")
