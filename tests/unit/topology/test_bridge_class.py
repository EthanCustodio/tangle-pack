"""Phase A — element identity, the combinatorial row, and bridge classes.

Three layers are pinned here, bottom up.

**ElementRef** (A.1) is the global name of a partition element. It has to hash,
compare and round-trip through both places that build it
(``PartitionInterval.ref`` and ``StablePartitionResult.ref``), and it has to
refuse to name an interval that was never stamped.

**The combinatorial row** (A.2) is the whole point of the phase: the side of the
stable branch a bridge approaches a crossing from, read off ``crossing_sign``
alone. Its agreement with the geometry on every bridge end (the independent
check on handedness) is the law ``row_of_end_is_geometry`` of
``tests/invariants/test_law_classes.py``; here it answers at an anchor and
rejects bad input.

**Bridge classes** (A.3) are the homotopy classes itineraries are written in:
the UNORDERED pair of elements a bridge connects, oriented anchor outward so a
member runs ``source -> target`` (``+1``) or back (``-1``). A loop (both ends
in one element) folds FORWARD into the class of the first non-loop bridge up
its image-chain parents (D2) and makes that class inert; a loop no image chain
reaches is an orphan, logged and left unclassed. The class laws on every law case (every bridge in one
class, anchor-outward orientation and member directions, the anchor bridge's
class, table order and tangle grouping, every orbit branch named, letters on
active classes only) live in ``tests/invariants/test_law_classes.py``; no
case-specific class fact is pinned. No law case
has a class connecting two tangles or a cross-branch tie, so the cross-branch
orientation and the connecting-classes-last grouping are pinned here on
hand-built references.
"""

from __future__ import annotations

import logging
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")  # headless: the session fixtures touch the plotting stack
import pytest

from helpers.logs import assert_logged
from helpers.fakes import bare_fixed_point, make_result
from tanglepack.topology.BridgeClass import (
    BridgeClass,
    BridgeClassEntry,
    bridge_classes,
    oriented_class,
)
from tanglepack.topology.MinimalTrellis import image_chain
from tanglepack.topology.StablePartition import (
    row_of_end,
    rows_of_bridge,
)
from tanglepack.topology.TopologyResults import ElementRef, PartitionInterval
from tanglepack.topology.Trellis import Trellis


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _classable_bridges(trellis):
    """Every bridge of the trellis that has a BridgeId (i.e. is not partial)."""
    return [bridge for bridge in trellis.bridges if bridge.id is not None]


def _all_partitions(session, fixed_points):
    """The stable partitions of the selected fixed points, gathered by the session."""
    selected = {id(fixed_point) for fixed_point in fixed_points}
    return [
        result
        for result in session.homotopy_partition(list(fixed_points)).as_list()
        if id(result.branch_key[0]) in selected
    ]


def _anchor_bridge_id(trellis):
    """The bridge running from a periodic point out to its first crossing."""
    for bridge in _classable_bridges(trellis):
        if trellis.intersection(bridge.id[0]).is_synthetic:
            return bridge.id
    return None


# --------------------------------------------------------------------------- #
# A.1 -- ElementRef
# --------------------------------------------------------------------------- #
def test_element_ref_hashes_compares_and_labels():
    fp = bare_fixed_point(3)
    key = (fp, "stable", 1, 0)
    a = ElementRef(key, "left", 2)
    b = ElementRef(key, "left", 2)
    other = ElementRef(key, "right", 2)

    assert a == b and hash(a) == hash(b)
    assert a != other
    assert len({a, b, other}) == 2
    assert a.fixed_point is fp
    assert a.orbit_index == 1 and a.branch_index == 0
    assert a.side == "left" and a.element_id == 2
    assert a.label == b.label and a.label != other.label
    assert str(a) == a.label


def test_element_ref_round_trips_through_both_builders(k10_partitioned):
    """``interval.ref`` and ``result.ref(id)`` must name the same element."""
    session, fp = k10_partitioned
    results = session.trellis(fp).stable_partitions
    assert results

    seen = set()
    for result in results:
        for interval in result.intervals:
            ref = interval.ref
            assert ref == result.ref(interval.element_id)
            assert ref.branch_key == result.branch_key
            assert ref.side == result.side
            assert ref.fixed_point is fp
            seen.add(ref)
    # An element's identity is global: no two elements of the branch share one.
    assert len(seen) == sum(len(r.intervals) for r in results)


def test_unstamped_interval_refuses_to_name_itself():
    raw = PartitionInterval(
        lo_id=None, hi_id=None, lo_cdist=0.0, hi_cdist=1.0,
        closed_lo=True, closed_hi=True,
    )
    with pytest.raises(ValueError):
        raw.ref

    half = PartitionInterval(
        lo_id=None, hi_id=None, lo_cdist=0.0, hi_cdist=1.0,
        closed_lo=True, closed_hi=True, element_id=0,
    )
    with pytest.raises(ValueError):
        half.ref  # stamped id but no branch/side


def test_result_ref_rejects_an_unknown_element_id(k10_partitioned):
    session, fp = k10_partitioned
    result = session.trellis(fp).stable_partitions[0]
    with pytest.raises(IndexError):
        result.ref(len(result.intervals))


# --------------------------------------------------------------------------- #
# A.2 -- the combinatorial row against the geometric one
# --------------------------------------------------------------------------- #
def test_row_of_end_answers_at_an_anchor(k10_partitioned):
    """The anchor is synthetic — no manifold nodes flank it — so the geometric
    row cannot always be read there, but the crossing sign always can."""
    session, _fp = k10_partitioned
    trellis = session.trellis()
    anchor_id = _anchor_bridge_id(trellis)
    assert anchor_id is not None, "the k=10 fixture must have an anchor bridge"
    assert trellis.intersection(anchor_id[0]).is_synthetic

    rows = rows_of_bridge(trellis, anchor_id)
    assert set(rows) <= {"left", "right"}
    assert rows == (
        row_of_end(trellis, anchor_id, "first"),
        row_of_end(trellis, anchor_id, "second"),
    )


def test_row_of_end_rejects_a_signless_crossing(k10_partitioned):
    session, _fp = k10_partitioned
    trellis = session.trellis()
    bridge_id = _classable_bridges(trellis)[0].id
    intersection = trellis.intersection(bridge_id[0])
    saved, intersection.crossing_sign = intersection.crossing_sign, 0
    try:
        with pytest.raises(ValueError):
            row_of_end(trellis, bridge_id, "first")
    finally:
        intersection.crossing_sign = saved


def test_row_of_end_rejects_an_unknown_endpoint_name(k10_partitioned):
    session, _fp = k10_partitioned
    trellis = session.trellis()
    bridge_id = _classable_bridges(trellis)[0].id
    with pytest.raises(ValueError):
        row_of_end(trellis, bridge_id, "middle")


# --------------------------------------------------------------------------- #
# A.3 -- bridge classes
# --------------------------------------------------------------------------- #
def _k10_table(session, fp):
    trellis = session.trellis()
    partitions = _all_partitions(session, [fp])
    return trellis, partitions, bridge_classes(trellis, partitions)


def test_loop_folds_into_the_class_whose_image_holds_it(k10_partitioned):
    """D2: the loop's ``folded_from`` is a non-loop member of the same class, and the
    image chain of that member holds the loop (directly, or through loop bridges)."""
    session, fp = k10_partitioned
    trellis, _partitions, table = _k10_table(session, fp)

    folded = [e for e in table.inert if e.loops]
    if not folded:
        pytest.skip("this build has no loop folded into an inert class")
    entry = folded[0]
    loop = entry.loops[0]
    ancestor = entry.member(loop.folded_from)
    assert not ancestor.is_loop

    parent = {
        frozenset(pair): bridge.id
        for bridge in _classable_bridges(trellis)
        for pair in image_chain(trellis, bridge.id) or ()
    }
    current = parent[frozenset(loop.bridge_id)]
    while current != ancestor.bridge_id:
        assert table.member_of(current).is_loop
        current = parent[frozenset(current)]


def test_orphan_loop_is_logged_and_unclassed(k10_partitioned, monkeypatch, caplog):
    """With no iterate registered no image chain holds a loop: it is an orphan, logged
    in a WARNING, recorded on the table and in no class (never lettered, D2)."""
    session, fp = k10_partitioned
    trellis = session.trellis()
    partitions = _all_partitions(session, [fp])
    before = bridge_classes(trellis, partitions)
    folded = [e for e in before.inert if e.loops]
    if not folded:
        pytest.skip("this build has no loop folded into an inert class")
    former = folded[0].bridge_class
    loop_id = folded[0].loops[0].bridge_id

    monkeypatch.setattr(trellis, "iterate", lambda intersection_id, n: None)
    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.BridgeClass"):
        table = bridge_classes(trellis, partitions)

    assert loop_id in table.orphans
    with pytest.raises(KeyError):
        table.entry_of(loop_id)
    assert not any(entry.bridge_class.is_loop for entry in table)
    assert_logged(caplog, logging.WARNING, "tanglepack.topology.BridgeClass")
    # The class the loop used to fold into is now active: no loop evidence.
    assert table[former].active


def test_restricted_to_keeps_the_entries_with_a_member_in_table_order(k10_partitioned):
    """The same entry objects, in table order, members kept whole; an id in either
    endpoint order matches and an unknown id is ignored."""
    session, fp = k10_partitioned
    _trellis, _partitions, table = _k10_table(session, fp)
    assert len(table) >= 2
    first, last = table.entries[0], table.entries[-1]
    wanted = [last.bridge_ids[0][::-1], first.bridge_ids[-1], (-1, -2)]

    restricted = table.restricted_to(wanted)
    assert [id(e) for e in restricted] == [id(first), id(last)]
    assert restricted.entries[0].members == first.members
    assert restricted.orphans == []
    assert len(table.restricted_to([])) == 0


def test_classes_are_keyed_by_element_pairs(k10_partitioned):
    session, fp = k10_partitioned
    trellis, partitions, table = _k10_table(session, fp)
    by_branch_side = {(r.branch_key, r.side): r for r in partitions}

    for entry in table:
        cls = entry.bridge_class
        assert isinstance(cls, BridgeClass)
        for member in entry.members:
            if member.is_loop:
                continue
            rows = rows_of_bridge(trellis, member.bridge_id)
            ends = (cls.source, cls.target) if member.direction > 0 else (cls.target, cls.source)
            for ref, row, intersection_id in zip(ends, rows, member.bridge_id):
                assert ref.side == row
                branch_key = trellis.intersection(intersection_id).manifold_b_key
                assert ref.branch_key == branch_key
                result = by_branch_side[(branch_key, row)]
                assert result.element_of_intersection[intersection_id] == ref.element_id


def test_bridge_ids_argument_restricts_the_classing(k10_partitioned):
    session, fp = k10_partitioned
    trellis = session.trellis()
    partitions = _all_partitions(session, [fp])
    chosen = sorted(bridge.id for bridge in _classable_bridges(trellis))[:3]

    table = bridge_classes(trellis, partitions, bridge_ids=chosen)
    assert table.bridge_ids == chosen


def test_table_lookups_symbols_and_report(k10_partitioned):
    session, fp = k10_partitioned
    _trellis, _partitions, table = _k10_table(session, fp)
    active, inert = table.active[0], table.inert[0]

    assert table.as_dict() == {e.bridge_class: e.bridge_ids for e in table}
    assert active.bridge_class in table and table[active.bridge_class] is active
    with pytest.raises(KeyError):
        table.entry_of((-1, -2))
    with pytest.raises(KeyError):
        active.member((-1, -2))

    # Unlettered: no symbol, but the report still prints.
    with pytest.raises(ValueError):
        active.symbol(active.bridge_ids[0])
    assert active.name == active.bridge_class.label

    active.letter = "a"
    forward = next(m for m in active.members if m.direction > 0)
    backward = next(m for m in active.members if m.direction < 0)
    assert active.symbol(forward.bridge_id) == "a"
    assert table.symbol(backward.bridge_id) == "a^-1"
    assert active.name == "a"

    report = table.describe()
    assert report
    for entry in table:
        assert entry.name in report


def test_a_missing_partition_names_the_branch(k10_partitioned):
    session, fp = k10_partitioned
    trellis = session.trellis()
    partitions = _all_partitions(session, [fp])

    bridge_id = _classable_bridges(trellis)[0].id
    row = row_of_end(trellis, bridge_id, "first")
    branch_key = trellis.intersection(bridge_id[0]).manifold_b_key
    kept = [
        result
        for result in partitions
        if (result.branch_key, result.side) != (branch_key, row)
    ]
    assert len(kept) < len(partitions)

    with pytest.raises(ValueError):
        bridge_classes(trellis, kept, bridge_ids=[bridge_id])
    # The same bridge classes once the missing (branch, row) partition is back.
    assert len(bridge_classes(trellis, partitions, bridge_ids=[bridge_id])) == 1


def test_duplicate_partitions_are_rejected(k10_partitioned):
    session, fp = k10_partitioned
    trellis = session.trellis()
    partitions = _all_partitions(session, [fp])
    with pytest.raises(ValueError):
        bridge_classes(trellis, partitions + partitions[:1])


# --------------------------------------------------------------------------- #
# A.4 -- inertness: a class is inert iff it maps to no active class
# --------------------------------------------------------------------------- #
from tanglepack.topology.BridgeClass import _resolve_inertness  # noqa: E402


_INERT_FP = bare_fixed_point(3)  # one fixed point, so refs compare by identity


def _ref(element_id: int, side: str = "left") -> ElementRef:
    return ElementRef((_INERT_FP, "stable", 0, 0), side, element_id)


def _class(a: int, b: int) -> BridgeClass:
    return BridgeClass(_ref(a), _ref(b))


def _entry(bridge_class, *, classes=(), loops=(), unresolved=()) -> BridgeClassEntry:
    return BridgeClassEntry(
        bridge_class,
        image_classes=list(classes),
        image_loops=[((0, 0), element) for element in loops],
        unresolved=list(unresolved),
    )


def test_resolve_inertness_fixed_point():
    """The rule on a hand-built table: a class whose only image is a loop is
    inert; one whose images are all inert is inert (propagated); a class that
    maps to itself stays active; mutual images stay active; a class with no
    resolvable image stays active for lack of evidence; and a class mapping to
    a virtual class with no entry (unknown activity) stays active."""
    loop_only = _class(0, 1)
    onto_inert = _class(1, 2)
    self_map = _class(2, 3)
    mutual_a, mutual_b = _class(3, 4), _class(4, 5)
    no_evidence = _class(5, 6)
    onto_unknown = _class(6, 7)
    unknown = _class(7, 8)  # no entry of its own

    entries = [
        _entry(loop_only, loops=[_ref(0)]),
        _entry(onto_inert, classes=[loop_only]),
        _entry(self_map, classes=[self_map, loop_only]),
        _entry(mutual_a, classes=[mutual_b]),
        _entry(mutual_b, classes=[mutual_a]),
        _entry(no_evidence, unresolved=[(1, 2)]),
        _entry(onto_unknown, classes=[unknown]),
    ]
    _resolve_inertness(entries)

    by_class = {entry.bridge_class: entry for entry in entries}
    assert by_class[loop_only].inert
    assert by_class[onto_inert].inert
    assert by_class[self_map].active
    assert by_class[mutual_a].active and by_class[mutual_b].active
    assert by_class[no_evidence].active and not by_class[no_evidence].has_image_evidence
    assert by_class[onto_unknown].active


# --------------------------------------------------------------------------- #
# Cross-branch orientation and tangle grouping (hand-built)
# --------------------------------------------------------------------------- #
def test_a_cross_branch_class_runs_from_the_anchor_nearer_element():
    """Across two orbit branches the SMALLER element id is the source, whichever
    branch it is on (anchor outward, 2026-09-30); a tie falls back to the
    run-stable branch order."""
    fp = bare_fixed_point(3)
    near = ElementRef((fp, "stable", 2, 0), "right", 0)  # later branch, element 0
    far = ElementRef((fp, "stable", 0, 0), "right", 1)  # earlier branch, element 1
    assert oriented_class(near, far, [fp]) == (BridgeClass(near, far), +1)
    assert oriented_class(far, near, [fp]) == (BridgeClass(near, far), -1)

    tie = ElementRef((fp, "stable", 0, 0), "right", 0)
    assert oriented_class(near, tie, [fp]) == (BridgeClass(tie, near), -1)


class _ClassTrellis:
    """The slice of a trellis ``bridge_classes`` reads: crossings, bridges, iterates.

    Every crossing is ``(unstable cdist, unstable key, stable key, sign)``;
    ``iterates`` maps ``(id, n)`` to the registered iterate (none by default, so
    every class is unresolved and active) and ``bridges`` lists the trellis's
    bridge ids.
    """

    #: The production anchor rules (no crossing here is an anchor).
    unstable_key_of = Trellis.unstable_key_of
    stable_key_at = Trellis.stable_key_at

    def __init__(
        self, crossings: dict, fixed_points: list, *, bridges=(), iterates=None
    ) -> None:
        self._crossings = crossings
        self.fixed_points = list(fixed_points)
        self.registry = {iid: self.intersection(iid) for iid in crossings}
        self.bridges = [SimpleNamespace(id=bid) for bid in bridges]
        self._iterates = dict(iterates or {})

    def intersection(self, iid: int) -> SimpleNamespace:
        """One crossing as the registry would hand it out."""
        cdist, unstable, stable, sign = self._crossings[iid]
        return SimpleNamespace(
            id=iid,
            unstable_cdist=cdist,
            stable_cdist=float(iid),
            manifold_a_key=unstable,
            manifold_b_key=stable,
            crossing_sign=sign,
            is_anchor=False,
        )

    def __contains__(self, iid) -> bool:
        """Every crossing of the fake is a crossing of its trellis (no pip cut)."""
        return iid in self._crossings

    def iterate(self, iid: int, n: int):
        """The registered iterate, or None."""
        return self._iterates.get((iid, n))

    def branch(self, key):
        """The crossings of one unstable branch, as ``_image_chain`` reads them."""
        on = [iid for iid, crossing in self._crossings.items() if crossing[1] == key]
        return SimpleNamespace(intersection_ids=sorted(on, key=lambda i: self._crossings[i][0]))


def test_classes_group_by_tangle_with_connecting_classes_last():
    """Tangle 0's classes, then tangle 1's, then the classes connecting two tangles
    (``tangle is None``) -- even when a connecting class has the smallest cdist."""
    outer, inner = bare_fixed_point(1, label="A"), bare_fixed_point(1, label="B")
    u_a, s_a = (outer, "unstable", 0, 0), (outer, "stable", 0, 0)
    u_b, s_b = (inner, "unstable", 0, 0), (inner, "stable", 0, 0)
    # sign +1 at the first end and -1 at the second puts both ends on the left row.
    crossings = {
        1: (5.0, u_a, s_a, +1), 2: (6.0, u_a, s_a, -1),  # tangle 0
        3: (1.0, u_b, s_b, +1), 4: (2.0, u_b, s_b, -1),  # tangle 1
        5: (0.5, u_a, s_a, +1), 6: (0.6, u_a, s_b, -1),  # connecting, smallest cdist
    }
    trellis = _ClassTrellis(crossings, [outer, inner])
    interval = (None, None, 0.0, 1.0, True, True)
    partitions = [
        make_result(s_a, "left", [interval] * 3, {1: 0, 2: 1, 5: 2}),
        make_result(s_b, "left", [interval] * 3, {3: 0, 4: 1, 6: 2}),
    ]
    table = bridge_classes(trellis, partitions, bridge_ids=[(1, 2), (3, 4), (5, 6)])

    assert [entry.tangle for entry in table] == [0, 1, None]
    assert [entry.bridge_ids for entry in table] == [[(1, 2)], [(3, 4)], [(5, 6)]]
    connecting = table.entries[-1]
    # Element 2 on both fixed points: the tie goes to the outer one (run-stable order).
    assert connecting.bridge_class.source.fixed_point is outer
    assert connecting.members[0].direction == +1


def _fold_born_trellis():
    """One period-1 tangle where a bridge's image holds a fold-born crossing.

    The parent ``(1, 2)`` joins left elements 0 and 2; its image runs from
    ``f(1) = 3`` to ``f(2) = 5`` across the fold-born crossing 4, which has no
    ``-1`` link. ``(3, 4)`` joins left elements 0 and 1; ``(4, 5)`` sits in right
    element 1 at both ends: a loop no preimage walk can reach.
    """
    fp = bare_fixed_point(1)
    u, s = (fp, "unstable", 0, 0), (fp, "stable", 0, 0)
    # Signs alternate along the unstable branch: +1 first / -1 second is the left row.
    crossings = {
        1: (1.0, u, s, +1), 2: (2.0, u, s, -1),
        3: (3.0, u, s, +1), 4: (4.0, u, s, -1), 5: (5.0, u, s, +1),
    }
    bridges = [(1, 2), (3, 4), (4, 5)]
    trellis = _ClassTrellis(
        crossings, [fp], bridges=bridges, iterates={(1, 1): 3, (2, 1): 5}
    )
    interval = (None, None, 0.0, 1.0, True, True)
    partitions = [
        make_result(s, "left", [interval] * 3, {1: 0, 2: 2, 3: 0, 4: 1, 5: 2}),
        make_result(s, "right", [interval] * 2, {1: 0, 2: 0, 3: 0, 4: 1, 5: 1}),
    ]
    return trellis, partitions, bridges


def test_a_fold_born_loop_folds_forward():
    """BCG-2's kernel: a loop whose crossing has no ``-1`` link folds into the class of
    the bridge whose image chain holds it, and is that class's loop evidence."""
    trellis, partitions, bridges = _fold_born_trellis()
    table = bridge_classes(trellis, partitions)

    assert table.orphans == []
    assert not any(entry.bridge_class.is_loop for entry in table)
    parent = table.entry_of((1, 2))
    loop = parent.member((4, 5))
    assert loop.is_loop and loop.folded_from == (1, 2)
    assert ((4, 5), loop.loop_element) in parent.image_loops
    assert sorted(table.bridge_ids) == sorted(bridges)


def test_an_image_pair_has_one_parent():
    """The map is a bijection, so the images of two disjoint arcs are disjoint: two
    bridges whose image chains share a pair break an invariant and raise."""
    trellis, partitions, _bridges = _fold_born_trellis()
    trellis.bridges.append(SimpleNamespace(id=(2, 3)))
    trellis._iterates[(3, 1)] = 4  # f(2, 3) = (5, 4): the pair (4, 5) again
    with pytest.raises(AssertionError):
        bridge_classes(trellis, partitions)
