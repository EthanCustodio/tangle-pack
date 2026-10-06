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
in one element) must fold into the class of the bridge it iterated from and
make that class inert. The class laws on every law case (every bridge in one
class, anchor-outward orientation and member directions, the anchor bridge's
class, table order and tangle grouping, every orbit branch named, letters on
active classes only) live in ``tests/invariants/test_law_classes.py``; the
k=10 and k=2.8 class facts are pinned once in ``tests/golden/``. No law case
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
    bridge_classes,
    oriented_class,
)
from tanglepack.topology.StablePartition import (
    row_of_end,
    rows_of_bridge,
)
from tanglepack.topology.TopologyResults import ElementRef, PartitionInterval


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


def test_loop_folds_into_its_preimage_class(k10_partitioned):
    """The loop's ancestor along the registry preimage chain is a member of the
    same class, and the chain really is the ``-1`` iterate of both ends."""
    session, fp = k10_partitioned
    trellis, _partitions, table = _k10_table(session, fp)

    loop = table.inert[0].loops[0]
    a, b = loop.bridge_id
    ancestor = table.inert[0].member(loop.folded_from)
    assert set(ancestor.bridge_id) == {trellis.iterate(a, -1), trellis.iterate(b, -1)}


def test_unresolved_loop_stands_alone_and_warns(k10_partitioned, monkeypatch, caplog):
    session, fp = k10_partitioned
    trellis = session.trellis()
    partitions = _all_partitions(session, [fp])
    before = bridge_classes(trellis, partitions)
    loop_id = before.inert[0].loops[0].bridge_id

    monkeypatch.setattr(trellis, "iterate", lambda intersection_id, n: None)
    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.BridgeClass"):
        table = bridge_classes(trellis, partitions)

    assert len(table) == 3
    entry = table.entry_of(loop_id)
    assert entry.bridge_class.is_loop and entry.inert
    assert entry.members[0].folded_from is None
    assert_logged(caplog, logging.WARNING, "tanglepack.topology.BridgeClass")
    # The class the loop used to fold into is now active: no loop evidence.
    former = before.inert[0].bridge_class
    assert table[former].active


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
from tanglepack.topology.BridgeClass import (  # noqa: E402  (section import)
    BridgeClassEntry,
    _resolve_inertness,
)


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


def test_an_unresolved_loop_class_is_inert_outright():
    entry = _entry(BridgeClass(_ref(2), _ref(2)))
    _resolve_inertness([entry])
    assert entry.bridge_class.is_loop and entry.inert




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
    """The slice of a trellis ``bridge_classes`` reads: crossings, no iterates.

    Every crossing is ``(unstable cdist, unstable key, stable key, sign)``; no
    crossing has a registered iterate, so every class is unresolved (active)
    and no loop occurs.
    """

    def __init__(self, crossings: dict, fixed_points: list) -> None:
        self._crossings = crossings
        self.fixed_points = list(fixed_points)

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
        )

    def iterate(self, iid: int, n: int) -> None:
        """No iterate is registered."""
        return None


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
