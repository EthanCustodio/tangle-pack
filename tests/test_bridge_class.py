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

**Bridge classes** (A.3) are the symbols the dual graph will be built over. Every
non-partial bridge must land in exactly one, the anchor bridge must be classed
like any other (it is the one whose geometry ``_row_at`` cannot read), and on the
nested fixture — which has no computed heteroclinic crossing — no class may mix
the two tangles.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless: the session fixtures touch the plotting stack
import pytest

from tanglepack.topology.BridgeClass import (
    BridgeClass,
    bridge_classes,
    class_sort_key,
    element_sort_key,
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
def test_every_non_partial_bridge_lands_in_exactly_one_class(k10_partitioned):
    session, fp = k10_partitioned
    trellis = session.trellis()
    classes = bridge_classes(trellis, _all_partitions(session, [fp]))

    expected = sorted(bridge.id for bridge in _classable_bridges(trellis))
    placed = sorted(bid for members in classes.values() for bid in members)
    assert placed == expected
    assert len(placed) == len(set(placed)), "a bridge landed in two classes"
    assert classes, "the fixture must produce at least one class"


def test_k10_anchor_bridge_is_classed_from_the_anchor_element(k10_partitioned):
    """The anchor bridge is a class like any other, and its first element is the
    one running from the periodic point (``lo_id is None``).

    The plan expected the anchor bridge to be a class of its own on this
    fixture; it is not. The k=10 partition is coarse — one hole orbit, three
    elements per side — so a second bridge whose ends land in the same two
    elements shares the class. That is a statement about the partition's
    fineness, not about the classing (it is exactly the D.4 diagnostic), so what
    is pinned here is the definitional part plus the fact that the sharing stays
    small.
    """
    session, fp = k10_partitioned
    trellis = session.trellis()
    partitions = _all_partitions(session, [fp])
    classes = bridge_classes(trellis, partitions)

    anchor_id = _anchor_bridge_id(trellis)
    owner = [cls for cls, members in classes.items() if anchor_id in members]
    assert len(owner) == 1, "the anchor bridge must land in exactly one class"
    assert len(classes[owner[0]]) <= 2

    by_branch_side = {(r.branch_key, r.side): r for r in partitions}
    first = owner[0].first
    interval = by_branch_side[(first.branch_key, first.side)].element(first.element_id)
    assert interval.lo_id is None, "the anchor's element must start at the anchor"


def test_classes_are_keyed_by_element_pairs(k10_partitioned):
    session, fp = k10_partitioned
    trellis = session.trellis()
    partitions = _all_partitions(session, [fp])
    classes = bridge_classes(trellis, partitions)
    by_branch_side = {(r.branch_key, r.side): r for r in partitions}

    for cls, members in classes.items():
        assert isinstance(cls, BridgeClass)
        for bridge_id in members:
            rows = rows_of_bridge(trellis, bridge_id)
            for ref, row, intersection_id in zip(
                (cls.first, cls.second), rows, bridge_id
            ):
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

    classes = bridge_classes(trellis, partitions, bridge_ids=chosen)
    placed = sorted(bid for members in classes.values() for bid in members)
    assert placed == chosen


def test_classes_and_members_come_out_in_the_documented_order(k10_partitioned):
    session, fp = k10_partitioned
    trellis = session.trellis()
    classes = bridge_classes(trellis, _all_partitions(session, [fp]))

    keys = list(classes)
    assert keys == sorted(
        keys, key=lambda cls: class_sort_key(cls, trellis.fixed_points)
    )
    for members in classes.values():
        assert members == sorted(members)

    # The key really is (fixed point, orbit, branch, side, element) twice over.
    for cls in keys:
        assert class_sort_key(cls, trellis.fixed_points) == (
            element_sort_key(cls.first, trellis.fixed_points)
            + element_sort_key(cls.second, trellis.fixed_points)
        )


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
    classes = bridge_classes(trellis, _all_partitions(session, [fp3, fp1]))

    expected = sorted(bridge.id for bridge in _classable_bridges(trellis))
    placed = sorted(bid for members in classes.values() for bid in members)
    assert placed == expected
    assert len(placed) == len(set(placed))

    seen = set()
    for cls in classes:
        owners = {id(cls.first.fixed_point), id(cls.second.fixed_point)}
        assert len(owners) == 1, f"class {cls} mixes two fixed points"
        seen |= owners
    assert seen == {id(fp3), id(fp1)}, "both tangles must contribute classes"


@pytest.mark.slow
def test_p3_period_three_classes_use_every_stable_branch(p3_partitioned):
    """A period-3 orbit anchors three stable branches; its bridges are spread
    over all of them, so the classes must name all three orbit indices."""
    session, fp3, fp1 = p3_partitioned
    trellis = session.trellis()
    classes = bridge_classes(trellis, _all_partitions(session, [fp3, fp1]))

    orbits = {
        ref.orbit_index
        for cls in classes
        for ref in (cls.first, cls.second)
        if ref.fixed_point is fp3
    }
    assert orbits == set(range(fp3.period))
