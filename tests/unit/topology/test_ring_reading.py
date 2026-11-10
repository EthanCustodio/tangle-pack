"""
A class read off the ring of its face (design note §A.8), on hand-built rings.

One disk face whose ring alternates four stable edges ``e[0..3]`` with four
kept bridges ``A, B, C, D`` (``B`` traversed against its dynamical direction).
Real ``Arc`` / ``Region`` / ``FaceNode`` objects carry the combinatorics;
nothing has geometry. The checks are the rules themselves: an arc through a
disk is the kept bridges on either side of it, read from its first end to its
second; it IS a kept bridge when it cuts that one alone off with no hole
between; and a tie of a class a word named is decided by the walks of those
bridges' classes, taken in turn.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from tanglepack.topology.BridgeClass import BridgeClass, BridgeClassEntry, BridgeClassTable, BridgeMember
from tanglepack.topology.DualGraph import DualGraph, FaceNode, face_side
from tanglepack.topology.DualWalk import Walk, WalkSearch, WalkStep, crossed_nodes, walk_itinerary
from tanglepack.topology.PartitionFamily import HomotopyPartition
from tanglepack.topology.SymbolicDynamics import _product_walk
from tanglepack.topology.TopologyResults import Arc, ElementRef, Region
from tanglepack.topology.diagnostics import classes_in_several_ring_positions

from helpers.fakes import bare_fixed_point, stable_key

A, B, C, D = (1, 2), (3, 4), (5, 6), (7, 0)


def _dual(face: FaceNode, markers: dict) -> DualGraph:
    """A dual graph over one hand-made face (the constructor needs a trellis)."""
    dual = object.__new__(DualGraph)
    dual.face_nodes = [face]
    dual.__dict__["markers"] = markers
    return dual


@pytest.fixture
def ring():
    """``(arcs, e)``: the ring ``e0 A e1 B^-1 e2 C e3 D`` and its stable edge keys."""
    key = stable_key(bare_fixed_point())
    unstable_key = (key[0], "unstable", 0, 0)
    stable = [Arc("stable", 2 * i, 2 * i + 1, key) for i in range(4)]
    bridges = [Arc("unstable", lo, hi, unstable_key, bridge_id=(lo, hi)) for lo, hi in (A, B, C, D)]
    bridges[1] = bridges[1].reversed()
    arcs = [arc for pair in zip(stable, bridges) for arc in pair]
    return arcs, [arc.edge_key for arc in stable]


def _face(arcs, **kind) -> FaceNode:
    return FaceNode(0, [Region(tuple(arc.tail_id for arc in arcs), arcs, None)], "region", **kind)


def _chord(face: FaceNode, a: tuple, b: tuple) -> tuple:
    side = {arc.edge_key: face_side(arc) for arc in face.faces[0].arcs}
    return (a, side[a], b, side[b])


def test_a_chord_is_the_kept_bridges_on_either_side_of_it(ring):
    """Read from the chord's first end to its second: along the ring on one side,
    against it on the other, each bridge signed by its own dynamical direction."""
    arcs, e = ring
    face = _face(arcs)
    dual = _dual(face, {})
    assert dual.ring_sides(face, _chord(face, e[0], e[2])) == (
        [(A, +1), (B, -1)],
        [(D, -1), (C, -1)],
    )
    # The same arc read the other way: the sides swap and every bridge turns round.
    assert dual.ring_sides(face, _chord(face, e[2], e[0])) == (
        [(C, +1), (D, +1)],
        [(B, +1), (A, -1)],
    )


def test_a_chord_cutting_one_kept_bridge_off_runs_alongside_it(ring):
    arcs, e = ring
    face = _face(arcs)
    dual = _dual(face, {})
    assert dual.alongside(face, _chord(face, e[0], e[1])) == [(A, +1)]
    assert dual.alongside(face, _chord(face, e[1], e[0])) == [(A, -1)]
    assert dual.alongside(face, _chord(face, e[1], e[2])) == [(B, -1)]
    # Two bridges on each side: a product of them, alongside neither.
    assert dual.alongside(face, _chord(face, e[0], e[2])) == []


def test_a_chord_through_a_strip_runs_alongside_both_its_bridges(ring):
    """A face that is a strip between two kept bridges: they are homotopic to each
    other, and a chord across the strip runs alongside both."""
    arcs, e = ring
    face = _face(arcs[:4])  # e0 A e1 B^-1
    dual = _dual(face, {})
    assert dual.alongside(face, _chord(face, e[0], e[1])) == [(A, +1), (B, +1)]


def test_a_hole_between_stops_the_reading_on_its_side(ring):
    """A hole beside a bridge, on this face's side of it, lies between the chord and
    the bridge: the chord is not homotopic to that side of the ring."""
    arcs, e = ring
    face = _face(arcs)
    bridge = arcs[1]
    dual = _dual(face, {("hole", 9, 0): ("beside", bridge.edge_key, face_side(bridge))})
    ahead, behind = dual.ring_sides(face, _chord(face, e[0], e[1]))
    assert ahead is None and [bridge_id for bridge_id, _along in behind] == [D, C, B]
    assert dual.alongside(face, _chord(face, e[0], e[1])) == []
    # A hole on the bridge's OTHER side lies in another face and stops nothing.
    other = "left" if face_side(bridge) == "right" else "right"
    dual = _dual(face, {("hole", 9, 0): ("beside", bridge.edge_key, other)})
    assert dual.alongside(face, _chord(face, e[0], e[1])) == [(A, +1)]


def test_the_ring_is_silent_off_a_disk(ring):
    """Only a disk determines an arc by its ends: the unbounded face, a merged face
    and a chord with both ends on one ring arc are not read."""
    arcs, e = ring
    chord = _chord(_face(arcs), e[0], e[1])
    unbounded = _face(arcs, is_unbounded=True)
    assert _dual(unbounded, {}).ring_sides(unbounded, chord) is None
    merged = _face(arcs)
    merged.faces.append(merged.faces[0])
    assert _dual(merged, {}).ring_sides(merged, chord) is None
    face = _face(arcs)
    assert _dual(face, {}).ring_sides(face, _chord(face, e[0], e[0])) is None


# ── the walk's itinerary ────────────────────────────────────────────────────


def _refs(count: int) -> list[ElementRef]:
    key = stable_key(bare_fixed_point())
    return [ElementRef(key, "right", index) for index in range(count)]


def _one_piece_walk(face, start_edge, goal_edge, itinerary) -> Walk:
    nodes = [SimpleNamespace(edge_key=edge) for edge in (start_edge, goal_edge)]
    return Walk([], tuple(itinerary), [face], start_node=nodes[0], goal_node=nodes[1])


def test_a_walk_piece_alongside_a_kept_bridge_is_read_at_its_ends(ring):
    """The piece joins the elements ``x, y``; the kept bridge it runs alongside sits
    in ``p, q``. The ring decides: the piece is read at the bridge's ends, turned
    round when it runs against the bridge. A piece already in the bridge's elements,
    and one alongside no single bridge, keep their own."""
    arcs, e = ring
    face = _face(arcs)
    dual = _dual(face, {})
    side = face_side(arcs[0])
    x, y, p, q = (ElementRef(ref.branch_key, side, ref.element_id) for ref in _refs(4))
    owners = {(A, "first"): p, (A, "second"): q}
    dual.trellis = object()
    dual.partition = SimpleNamespace(owner_at_end=lambda _t, bridge, end: owners[(bridge, end)])
    dual.naming = SimpleNamespace(parent_of=lambda ref: ref)

    def walk(start, goal, itinerary):
        return _one_piece_walk(face, start, goal, itinerary)

    assert walk_itinerary(dual, walk(e[0], e[1], (x, y))) == (p, q)
    assert walk_itinerary(dual, walk(e[1], e[0], (x, y))) == (q, p)
    assert walk_itinerary(dual, walk(e[0], e[1], (p, q))) == (p, q)
    assert walk_itinerary(dual, walk(e[0], e[2], (x, y))) == (x, y)


def test_in_a_strip_the_pair_says_which_of_the_two_bridges(ring):
    """Two kept bridges bounding a strip may be two classes of one homotopy class
    (an anchor bridge and its neighbour). A piece whose own elements name either
    keeps them; one naming neither is read at the first."""
    arcs, e = ring
    face = _face(arcs[:4])
    dual = _dual(face, {})
    side = face_side(arcs[0])
    x, y, p, q, r, s = (ElementRef(ref.branch_key, side, ref.element_id) for ref in _refs(6))
    owners = {(A, "first"): p, (A, "second"): q, (B, "first"): r, (B, "second"): s}
    dual.trellis = object()
    dual.partition = SimpleNamespace(owner_at_end=lambda _t, bridge, end: owners[(bridge, end)])
    dual.naming = SimpleNamespace(parent_of=lambda ref: ref)
    chord = _chord(face, e[0], e[1])
    assert dual.read_off_the_ring(face, chord, [p, q]) is None
    assert dual.read_off_the_ring(face, chord, [r, s]) is None
    assert dual.read_off_the_ring(face, chord, [x, y]) == (A, +1, [p, q])


def test_a_computed_bridge_alongside_a_kept_one_is_read_at_its_ends(p3_partitioned):
    """``Trellis.parallel_bridges`` is the one table ``owner_at_end`` consults: a
    listed bridge is read at the kept bridge's ends, turned round for ``-1``."""
    session, fp3, _fp1 = p3_partitioned
    trellis = session.trellis(fp3)
    family = HomotopyPartition.from_results(trellis.stable_partitions, trellis=trellis)
    ids = [tuple(b.id) for b in trellis.bridges if b.id is not None]
    ends = {bid: [family.owner_at_end(trellis, bid, end) for end in ("first", "second")] for bid in ids}
    mine, wall = next((a, b) for a in ids for b in ids if set(ends[a]).isdisjoint(ends[b]))
    try:
        trellis.parallel_bridges = {mine: (wall, +1)}
        assert [family.owner_at_end(trellis, mine, end) for end in ("first", "second")] == ends[wall]
        trellis.parallel_bridges = {mine: (wall, -1)}
        assert [family.owner_at_end(trellis, mine, end) for end in ("first", "second")] == ends[wall][::-1]
        assert [family.owner_at_end(trellis, wall, end) for end in ("first", "second")] == ends[wall]
    finally:
        trellis.parallel_bridges = {}


# ── the tie of a class a word named ─────────────────────────────────────────


def _walk(*crossed, multiplicity: int = 1) -> Walk:
    steps = [WalkStep(SimpleNamespace(key=key), "left", "right", None, None) for key in crossed]
    return Walk(steps, (), [], multiplicity=multiplicity)


def _record(cls, walks, *, status="resolved", named_by=None, itinerary=None):
    search = None if walks is None else WalkSearch(cls.source, cls.target, status, list(walks))
    return SimpleNamespace(
        bridge_class=cls, search=search, named_by=named_by, itinerary=itinerary,
        entry=BridgeClassEntry(cls),
    )


@pytest.fixture
def named(ring):
    """A namer whose one-piece walk is a chord ``e0 -> e2`` of the ring, naming the
    class ``x -> y``; the kept bridges ``A, B, C, D`` in four classes of their own.

    Returns ``(dual, records, the named class, the kept bridges' classes)``."""
    arcs, e = ring
    face = _face(arcs)
    dual = _dual(face, {})
    refs = _refs(12)
    x, y = refs[0], refs[1]
    classes = {bridge: BridgeClass(refs[2 + 2 * i], refs[3 + 2 * i]) for i, bridge in enumerate((A, B, C, D))}
    table = BridgeClassTable(
        [BridgeClassEntry(cls, [BridgeMember(bridge, +1)]) for bridge, cls in classes.items()]
    )
    dual.minimal = SimpleNamespace(classes=table)
    dual.trellis = SimpleNamespace(fixed_points=[x.branch_key[0]])
    namer_class = BridgeClass(refs[10], refs[11])
    namer_walk = _one_piece_walk(face, e[0], e[2], ())
    side = face_side(arcs[0])
    namer_walk.itinerary = (ElementRef(x.branch_key, side, 0), ElementRef(x.branch_key, side, 1))
    namer = _record(namer_class, [namer_walk], itinerary=(x, y))
    records = {namer_class: namer}
    return dual, records, BridgeClass(x, y), classes, namer_class


_NAMING = SimpleNamespace(parent_of=lambda ref: ref)


def test_a_tie_is_decided_by_the_product_of_the_walks_beside_the_naming_piece(named):
    """The naming chord has ``A B^-1`` on one side: the image crosses where ``A``'s
    walk does, then where ``B``'s does backwards. The other side has no walk to
    read (``D``'s class walked no dual graph) and is simply unavailable."""
    dual, records, cls, classes, namer_class = named
    records[classes[A]] = _record(classes[A], [_walk("n1", "n2")])
    records[classes[B]] = _record(classes[B], [_walk("n4", "n3")])
    records[classes[C]] = _record(classes[C], [_walk("n9")])
    records[classes[D]] = _record(classes[D], None)
    right, wrong = _walk("n1", "n2", "n3", "n4"), _walk("n1", "n2", "n4", "n3")
    cd = _record(cls, [wrong, right], status="ambiguous", named_by=namer_class)
    assert _product_walk(dual, cd, records, _NAMING) is right


def test_a_crossing_made_and_undone_cancels_in_the_product(named):
    dual, records, cls, classes, namer_class = named
    records[classes[A]] = _record(classes[A], [_walk("n1", "n2")])
    records[classes[B]] = _record(classes[B], [_walk("n3", "n2")])  # backwards: n2, n3
    right, wrong = _walk("n1", "n3"), _walk("n1", "n2")
    cd = _record(cls, [wrong, right], status="ambiguous", named_by=namer_class)
    assert _product_walk(dual, cd, records, _NAMING) is right


def test_both_readable_sides_must_agree_and_single_out_one_walk(named):
    dual, records, cls, classes, namer_class = named
    records[classes[A]] = _record(classes[A], [_walk("n1")])
    records[classes[B]] = _record(classes[B], [_walk()])
    records[classes[C]] = _record(classes[C], [_walk("n7")])
    records[classes[D]] = _record(classes[D], [_walk()])
    # One side says n1, the other n7 (read backwards along D, C): no walk is both.
    cd = _record(cls, [_walk("n1"), _walk("n7")], status="ambiguous", named_by=namer_class)
    assert _product_walk(dual, cd, records, _NAMING) is None


@pytest.mark.parametrize("spoil", ["no_namer", "namer_ambiguous", "side_unresolved", "no_match", "several_paths"])
def test_the_product_is_silent_where_it_cannot_be_read(named, spoil):
    """No namer, a namer with no single decided walk, a kept bridge whose class has
    none on every side, no tied walk equal to the product, or one standing for
    several face paths: the tie stays a tie."""
    dual, records, cls, classes, namer_class = named
    records[classes[A]] = _record(classes[A], [_walk("n1")])
    records[classes[B]] = _record(classes[B], [_walk()])
    walks = [_walk("n1"), _walk("n2")]
    if spoil == "no_namer":
        namer_class = None
    elif spoil == "namer_ambiguous":
        records[namer_class].search.status = "ambiguous"
    elif spoil == "side_unresolved":
        records[classes[B]].search.status = "ambiguous"
    elif spoil == "no_match":
        walks = [_walk("n5"), _walk("n2")]
    elif spoil == "several_paths":
        walks = [_walk("n1", multiplicity=2), _walk("n2")]
    cd = _record(cls, walks, status="ambiguous", named_by=namer_class)
    assert _product_walk(dual, cd, records, _NAMING) is None


def test_crossed_nodes_are_the_keys_of_the_nodes_a_walk_crosses():
    assert crossed_nodes(_walk("n1", "n2")) == ["n1", "n2"]
    assert crossed_nodes(_walk()) == []


# ── members in several ring positions ───────────────────────────────────────


def test_a_class_with_bridges_in_several_ring_positions_is_reported():
    """A member that is kept or runs alongside a kept bridge is in one position; a
    member cutting several bridges off is in the position those bridges name. Two
    positions in one class are reported; a member the ring does not read counts
    for nothing."""
    refs = _refs(6)
    kept, beside, product, other_product, unread = (1, 2), (3, 4), (5, 6), (7, 8), (9, 10)
    readings = {
        beside: [[(kept, +1)], [(A, +1), (B, +1)]],
        product: [[(A, +1), (B, +1)], [(C, +1), (D, +1)]],
        other_product: [[(A, +1), (C, +1)], [(B, +1), (D, +1)]],
    }

    def entry(index, *bridges):
        return BridgeClassEntry(
            BridgeClass(refs[index], refs[index + 1]), [BridgeMember(b, +1) for b in bridges]
        )

    one = entry(0, kept, beside, unread)
    products = entry(2, product, product[::-1])
    mixed = entry(4, kept, product)
    two_products = entry(1, product, other_product)
    dual = SimpleNamespace(
        unstable_edges={kept: None},
        minimal=SimpleNamespace(classes=BridgeClassTable([one, products, mixed, two_products])),
        chord_of=lambda bridge: ("face", bridge) if tuple(sorted(bridge)) in readings else None,
        ring_sides=lambda _face, bridge: tuple(readings[tuple(sorted(bridge))]),
        alongside=lambda _face, bridge: [(kept, +1)] if tuple(sorted(bridge)) == beside else [],
    )
    assert classes_in_several_ring_positions(dual) == [mixed, two_products]
