"""Phase A — element identity, the combinatorial row, and bridge classes.

Three layers are pinned here, bottom up.

**ElementRef** (A.1) is the global name of a partition element. It has to hash,
compare and round-trip through both places that build it
(``PartitionInterval.ref`` and ``StablePartitionResult.ref``), and it has to
refuse to name an interval that was never stamped.

**The combinatorial row** (A.2) is the whole point of the phase: the side of the
stable branch a bridge approaches a crossing from, read off ``crossing_sign``
alone. It is validated against the geometric ``_row_at`` — which walks real
manifold nodes and takes real cross products — on every bridge end of both
fixtures where the geometry has an answer. That comparison is also the
independent check on handedness: the partition's ``left`` and the crossing
sign's ``left`` must be the same orientation of the plane, and a flip in either
shows up as a wholesale disagreement here.

**Bridge classes** (A.3) are the homotopy classes itineraries are written in:
the UNORDERED pair of elements a bridge connects, oriented anchor outward so a
member runs ``source -> target`` (``+1``) or back (``-1``). Every non-partial
bridge must land in exactly one, the anchor bridge must be classed like any
other (it is the one whose geometry ``_row_at`` cannot read), a loop (both ends
in one element) must fold into the class of the bridge it iterated from and
make that class inert, and on the nested fixture — which has no computed
heteroclinic crossing — no class may mix the two tangles.
"""

from __future__ import annotations

import logging

import matplotlib

matplotlib.use("Agg")  # headless: the session fixtures touch the plotting stack
import pytest

from tanglepack.topology.BridgeClass import (
    BridgeClass,
    BridgeClassTable,
    bridge_classes,
    class_sort_key,
    element_sort_key,
    oriented_class,
)
from tanglepack.topology.StablePartition import (
    _row_at,
    _row_polyline,
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
    """The stable partitions of every per-fixed-point trellis, concatenated."""
    partitions = []
    for fixed_point in fixed_points:
        partitions.extend(session.trellis(fixed_point).stable_partitions)
    return partitions


def _row_agreements(trellis) -> tuple[int, list[str]]:
    """Compare row_of_end with the geometric _row_at on every readable end.

    Returns the number of ends where the geometry had an answer, together with a
    description of each disagreement.
    """
    compared = 0
    mismatches: list[str] = []
    for bridge in _classable_bridges(trellis):
        poly, _oriented = _row_polyline(trellis, bridge)
        if poly is None:
            continue
        for endpoint, intersection_id in zip(("first", "second"), bridge.id):
            geometric = _row_at(trellis, poly, intersection_id)
            if geometric is None:
                continue  # no manifold data / degenerate end: nothing to compare
            compared += 1
            combinatorial = row_of_end(trellis, bridge.id, endpoint)
            if combinatorial != geometric:
                mismatches.append(
                    f"bridge {bridge.id} {endpoint} end (crossing "
                    f"{intersection_id}, sign "
                    f"{trellis.intersection(intersection_id).crossing_sign}): "
                    f"combinatorial {combinatorial}, geometric {geometric}"
                )
    return compared, mismatches


def _anchor_bridge_id(trellis):
    """The bridge running from a periodic point out to its first crossing."""
    for bridge in _classable_bridges(trellis):
        if trellis.intersection(bridge.id[0]).is_synthetic:
            return bridge.id
    return None


# --------------------------------------------------------------------------- #
# A.1 -- ElementRef
# --------------------------------------------------------------------------- #
class _FakeFixedPoint:
    """The only thing an ElementRef asks of a fixed point is its period."""

    period = 3


def test_element_ref_hashes_compares_and_labels():
    fp = _FakeFixedPoint()
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
    assert a.label == "p3@1.0/L#2"
    assert str(a) == a.label
    assert other.label == "p3@1.0/R#2"


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
    with pytest.raises(ValueError, match="not been stamped"):
        raw.ref

    half = PartitionInterval(
        lo_id=None, hi_id=None, lo_cdist=0.0, hi_cdist=1.0,
        closed_lo=True, closed_hi=True, element_id=0,
    )
    with pytest.raises(ValueError, match="not been stamped"):
        half.ref  # stamped id but no branch/side


def test_result_ref_rejects_an_unknown_element_id(k10_partitioned):
    session, fp = k10_partitioned
    result = session.trellis(fp).stable_partitions[0]
    with pytest.raises(IndexError):
        result.ref(len(result.intervals))


# --------------------------------------------------------------------------- #
# A.2 -- the combinatorial row against the geometric one
# --------------------------------------------------------------------------- #
def test_row_of_end_agrees_with_the_geometry_on_k10(k10_partitioned):
    session, _fp = k10_partitioned
    compared, mismatches = _row_agreements(session.trellis())
    assert not mismatches, "\n".join(mismatches)
    assert compared >= 10, f"only {compared} bridge ends had a geometric row"


@pytest.mark.slow
def test_row_of_end_agrees_with_the_geometry_on_p3(p3_partitioned):
    session, _fp3, _fp1 = p3_partitioned
    compared, mismatches = _row_agreements(session.trellis())
    assert not mismatches, "\n".join(mismatches)
    assert compared >= 10, f"only {compared} bridge ends had a geometric row"


def test_same_branch_bridges_never_mismatch_on_k10(k10_partitioned):
    """Invariant I2, combinatorially: a bridge whose ends share a stable branch
    approaches both from the same side of it."""
    session, _fp = k10_partitioned
    _assert_no_same_branch_mismatch(session.trellis())


@pytest.mark.slow
def test_same_branch_bridges_never_mismatch_on_p3(p3_partitioned):
    session, _fp3, _fp1 = p3_partitioned
    _assert_no_same_branch_mismatch(session.trellis())


def _assert_no_same_branch_mismatch(trellis):
    checked = 0
    for bridge in _classable_bridges(trellis):
        keys = [trellis.intersection(i).manifold_b_key for i in bridge.id]
        if keys[0] != keys[1]:
            continue
        first, second = rows_of_bridge(trellis, bridge.id)
        assert first == second, (
            f"bridge {bridge.id} has both ends on stable branch {keys[0][1:]} "
            f"but rows {first} / {second}"
        )
        checked += 1
    assert checked, "the fixture must have at least one same-branch bridge"


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
        with pytest.raises(ValueError, match="no crossing sign"):
            row_of_end(trellis, bridge_id, "first")
    finally:
        intersection.crossing_sign = saved


def test_row_of_end_rejects_an_unknown_endpoint_name(k10_partitioned):
    session, _fp = k10_partitioned
    trellis = session.trellis()
    bridge_id = _classable_bridges(trellis)[0].id
    with pytest.raises(ValueError, match="first"):
        row_of_end(trellis, bridge_id, "middle")


# --------------------------------------------------------------------------- #
# A.3 -- bridge classes
# --------------------------------------------------------------------------- #
def _k10_table(session, fp):
    trellis = session.trellis()
    partitions = _all_partitions(session, [fp])
    return trellis, partitions, bridge_classes(trellis, partitions)


def _element_ends(trellis, partitions, bridge_id):
    """The element at each end of a bridge, on that end's row."""
    by_branch_side = {(r.branch_key, r.side): r for r in partitions}
    refs = []
    for endpoint, intersection_id in zip(("first", "second"), bridge_id):
        row = row_of_end(trellis, bridge_id, endpoint)
        branch_key = trellis.intersection(intersection_id).manifold_b_key
        result = by_branch_side[(branch_key, row)]
        refs.append(result.ref(result.element_of_intersection[intersection_id]))
    return tuple(refs)


def test_every_non_partial_bridge_lands_in_exactly_one_class(k10_partitioned):
    session, fp = k10_partitioned
    trellis, _partitions, table = _k10_table(session, fp)

    expected = sorted(bridge.id for bridge in _classable_bridges(trellis))
    placed = [bid for entry in table for bid in entry.bridge_ids]
    assert sorted(placed) == expected
    assert len(placed) == len(set(placed)), "a bridge landed in two classes"
    assert table.bridge_ids == expected
    assert isinstance(table, BridgeClassTable) and len(table) > 0


def test_k10_has_one_active_and_one_inert_class(k10_partitioned):
    """The k=10 fixture (one hole orbit, three elements per side) has exactly
    two homotopy classes: the zone-side pair with four bridges, two each way,
    and the exterior pair whose forward image is a loop pushed into the anchor
    element — folded in, making that class inert."""
    session, fp = k10_partitioned
    _trellis, _partitions, table = _k10_table(session, fp)

    assert len(table) == 2
    assert len(table.active) == 1 and len(table.inert) == 1
    active, inert = table.active[0], table.inert[0]

    directions = sorted(m.direction for m in active.members)
    assert directions == [-1, -1, 1, 1]
    assert not active.loops

    loops = inert.loops
    assert len(loops) == 1
    loop = loops[0]
    assert loop.folded_from is not None
    ancestor = inert.member(loop.folded_from)
    assert ancestor.direction == -1, "the loop's ancestor runs outer -> anchor element"
    assert loop.loop_element == inert.bridge_class.source, "the loop sits in the anchor element"
    assert sorted(m.direction for m in inert.members) == [-1, 0, 1]

    for entry in table:
        assert not entry.bridge_class.is_loop
        assert entry.bridge_class.source.side == entry.bridge_class.target.side
    assert active.bridge_class.source.side != inert.bridge_class.source.side


def test_k10_anchor_bridge_runs_source_to_target(k10_partitioned):
    """The anchor bridge leaves the periodic point (element with ``lo_id is
    None``) outward, so it is a ``+1`` member of the active class."""
    session, fp = k10_partitioned
    trellis, partitions, table = _k10_table(session, fp)

    anchor_id = _anchor_bridge_id(trellis)
    entry = table.entry_of(anchor_id)
    assert entry.active
    assert entry.member(anchor_id).direction == +1

    by_branch_side = {(r.branch_key, r.side): r for r in partitions}
    source = entry.bridge_class.source
    interval = by_branch_side[(source.branch_key, source.side)].element(source.element_id)
    assert interval.lo_id is None, "the class's source must start at the anchor"


def test_direction_matches_the_element_order_at_the_ends(k10_partitioned):
    session, fp = k10_partitioned
    trellis, partitions, table = _k10_table(session, fp)

    for entry in table:
        cls = entry.bridge_class
        assert element_sort_key(cls.source, trellis.fixed_points) <= element_sort_key(
            cls.target, trellis.fixed_points
        )
        for member in entry.members:
            x, y = _element_ends(trellis, partitions, member.bridge_id)
            if member.is_loop:
                assert x == y == member.loop_element
                continue
            assert {x, y} == {cls.source, cls.target}
            assert member.direction == (+1 if x == cls.source else -1)


def test_oriented_class_orders_anchor_outward(k10_partitioned):
    session, fp = k10_partitioned
    trellis, partitions, table = _k10_table(session, fp)
    cls = table.active[0].bridge_class
    fps = trellis.fixed_points

    assert oriented_class(cls.source, cls.target, fps) == (cls, +1)
    assert oriented_class(cls.target, cls.source, fps) == (cls, -1)
    loop_class, direction = oriented_class(cls.source, cls.source, fps)
    assert direction == 0 and loop_class.is_loop
    assert cls.source.element_id < cls.target.element_id


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
    assert "loop bridge" in caplog.text and str(loop_id) in caplog.text
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


def test_classes_and_members_come_out_in_the_documented_order(k10_partitioned):
    session, fp = k10_partitioned
    trellis, _partitions, table = _k10_table(session, fp)

    keys = table.classes
    by_class = {entry.bridge_class: entry for entry in table}
    assert keys == sorted(
        keys,
        key=lambda cls: (
            by_class[cls].min_unstable_cdist, class_sort_key(cls, trellis.fixed_points)
        ),
    )
    # The anchor bridge (unstable cdist 0) puts its class first.
    assert table.entries[0].min_unstable_cdist == 0.0
    for entry in table:
        assert entry.bridge_ids == sorted(entry.bridge_ids)
        assert entry.min_unstable_cdist == min(
            trellis.intersection(bid[0]).unstable_cdist for bid in entry.bridge_ids
        )

    # The key really is (fixed point, orbit, branch, side, element) twice over.
    for cls in keys:
        assert class_sort_key(cls, trellis.fixed_points) == (
            element_sort_key(cls.source, trellis.fixed_points)
            + element_sort_key(cls.target, trellis.fixed_points)
        )


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
    with pytest.raises(ValueError, match="no letter"):
        active.symbol(active.bridge_ids[0])
    assert active.name == active.bridge_class.label

    active.letter = "a"
    forward = next(m for m in active.members if m.direction > 0)
    backward = next(m for m in active.members if m.direction < 0)
    assert active.symbol(forward.bridge_id) == "a"
    assert table.symbol(backward.bridge_id) == "a^-1"
    assert active.name == "a"

    report = table.describe()
    assert "2 bridge class(es): 1 active, 1 inert" in report
    assert "a:" in report and "a^-1" in report
    assert "inert" in report and "folded from" in report
    assert inert.bridge_class.label in report


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

    with pytest.raises(ValueError) as excinfo:
        bridge_classes(trellis, kept, bridge_ids=[bridge_id])
    message = str(excinfo.value)
    assert str(branch_key[1:]) in message and repr(row) in message


def test_duplicate_partitions_are_rejected(k10_partitioned):
    session, fp = k10_partitioned
    trellis = session.trellis()
    partitions = _all_partitions(session, [fp])
    with pytest.raises(ValueError, match="two partitions"):
        bridge_classes(trellis, partitions + partitions[:1])


@pytest.mark.slow
def test_p3_classes_cover_both_tangles_and_mix_neither(p3_partitioned):
    """The nested fixture has no computed heteroclinic crossing, so every class
    must live entirely inside one tangle — and both tangles must appear."""
    session, fp3, fp1 = p3_partitioned
    trellis = session.trellis()
    table = bridge_classes(trellis, _all_partitions(session, [fp3, fp1]))

    expected = sorted(bridge.id for bridge in _classable_bridges(trellis))
    assert table.bridge_ids == expected

    seen = set()
    for cls in table.classes:
        owners = {id(cls.source.fixed_point), id(cls.target.fixed_point)}
        assert len(owners) == 1, f"class {cls} mixes two fixed points"
        seen |= owners
    assert seen == {id(fp3), id(fp1)}, "both tangles must contribute classes"


@pytest.mark.slow
def test_p3_period_three_classes_use_every_stable_branch(p3_partitioned):
    """A period-3 orbit anchors three stable branches; its bridges are spread
    over all of them, so the classes must name all three orbit indices."""
    session, fp3, fp1 = p3_partitioned
    trellis = session.trellis()
    table = bridge_classes(trellis, _all_partitions(session, [fp3, fp1]))

    orbits = {
        ref.orbit_index
        for cls in table.classes
        for ref in (cls.source, cls.target)
        if ref.fixed_point is fp3
    }
    assert orbits == {0, 1, 2}


# --------------------------------------------------------------------------- #
# A.4 -- inertness: a class is inert iff it maps to no active class
# --------------------------------------------------------------------------- #
from tanglepack.topology.BridgeClass import (  # noqa: E402  (section import)
    BridgeClassEntry,
    _resolve_inertness,
)


_INERT_FP = _FakeFixedPoint()  # one fixed point, so refs compare by identity


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


def test_k10_inert_class_rests_on_its_folded_loop_despite_an_unresolved_member(k10_partitioned):
    """k=10: the exterior class has one member with no registered image and one
    whose image is the folded loop; the loop is the evidence, the unresolved
    member is no counter-evidence. The active class maps over itself."""
    session, fp = k10_partitioned
    _trellis, _partitions, table = _k10_table(session, fp)
    active, inert = table.active[0], table.inert[0]

    assert inert.image_loops and not inert.image_classes
    assert inert.unresolved, "an exterior member's image is not grown yet"
    assert {pair for pair, _element in inert.image_loops} == {m.bridge_id for m in inert.loops}
    assert active.bridge_class in active.image_classes
    assert active.has_image_evidence


def test_k28_has_one_active_and_two_inert_classes(k28_partitioned):
    """The blasted k=2.8 tangle (holes backward only): the anchor bridge's class
    is the only active one; the exterior class is inert through a VIRTUAL loop
    (its image pair is registered but no bridge spans it); the interior class is
    inert through the folded blast-child loop."""
    session, fp = k28_partitioned
    trellis = session.trellis()
    partitions = _all_partitions(session, [fp])
    table = bridge_classes(trellis, partitions)

    assert len(table) == 3
    assert len(table.active) == 1 and len(table.inert) == 2
    active = table.active[0]
    anchor = next(iid for iid in trellis.own_intersection_ids
                  if trellis.intersection(iid).unstable_cdist == 0.0)
    assert [m.bridge_id[0] for m in active.members] == [anchor]
    assert active.bridge_class in active.image_classes
    assert {e.bridge_class for e in table.inert} <= set(active.image_classes)

    virtual, folded = sorted(table.inert, key=lambda e: len(e.loops))
    assert not virtual.loops and len(virtual.members) == 1
    assert len(virtual.image_loops) == 1 and not virtual.image_classes
    (pair, element), = virtual.image_loops
    assert trellis.bridge_between(*pair) is None, "the virtual loop has no bridge object"
    assert element.side == virtual.bridge_class.source.side

    assert len(folded.loops) == 1 and len(folded.members) == 2
    loop = folded.loops[0]
    assert loop.folded_from == next(m.bridge_id for m in folded.members if not m.is_loop)
    assert folded.image_loops == [(loop.bridge_id, loop.loop_element)]
    for entry in table:
        assert not entry.bridge_class.is_loop
