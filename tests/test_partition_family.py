"""
The partition-family base class and the homotopy family.

Pins: the homotopy family is exactly the trellis's own partition; its
signature is the session's cache signature; midpoint ownership is unique
across the full arrangement; the marks round trip reproduces every result
exactly (which validates the local interval builder against
``StablePartition._build_intervals`` without importing it).
"""

from __future__ import annotations

import pytest

from tanglepack.loom.TangleSession import TangleSession
from tanglepack.topology.PartitionFamily import (
    HomotopyPartition,
    IteratedHomotopyPartition,
    PartitionFamily,
    _intervals_from_marks,
    _marks_of,
)


def _boundaries(result):
    return [
        (iv.lo_id, iv.hi_id, iv.closed_lo, iv.closed_hi, iv.lo_cdist, iv.hi_cdist)
        for iv in result.intervals
    ]


def test_build_reproduces_the_trellis_partition(k10_partitioned):
    session, fp = k10_partitioned
    trellis = session.trellis(fp)
    family = HomotopyPartition.build(trellis)
    assert family.kind == "homotopy"
    assert len(family) == len(trellis.stable_partitions)
    for stored in trellis.stable_partitions:
        built = family.result(stored.branch_key, stored.side)
        assert _boundaries(built) == _boundaries(stored)
        assert built.element_of_intersection == stored.element_of_intersection


def test_from_results_signature_matches_the_session_signature(k10_partitioned):
    session, fp = k10_partitioned
    gathered = session._gathered_partitions()
    family = HomotopyPartition.from_results(gathered, trellis=session.trellis())
    assert family.signature() == TangleSession._partition_signature(gathered)
    assert family.as_list() == gathered


def test_element_at_owns_every_arc_midpoint_exactly_once(k10_partitioned):
    session, fp = k10_partitioned
    family = HomotopyPartition.from_results(
        session._gathered_partitions(), trellis=session.trellis()
    )
    arrangement = session.arrangement()
    seen = 0
    for face in arrangement.faces:
        for arc in face.arcs:
            if arc.kind != "stable":
                continue
            lo = session.trellis().intersection(arc.lo_id).stable_cdist
            hi = session.trellis().intersection(arc.hi_id).stable_cdist
            mid = 0.5 * (lo + hi)
            for side in ("left", "right"):
                ref = family.element_at(arc.branch_key, side, mid)
                assert family.element(ref) is family.interval_at(arc.branch_key, side, mid)
                assert ref.branch_key == arc.branch_key and ref.side == side
                seen += 1
    assert seen > 0


def test_element_at_refuses_a_cdist_outside_the_partition(k10_partitioned):
    session, fp = k10_partitioned
    family = HomotopyPartition.build(session.trellis(fp))
    key = family.branch_keys[0]
    end = family.result(key, "left").intervals[-1].hi_cdist
    with pytest.raises(ValueError, match="own stable cdist"):
        family.element_at(key, "left", end + 10.0)
    with pytest.raises(ValueError, match="own stable cdist"):
        family.element_at(key, "left", -1.0)


def test_duplicate_and_missing_results_raise(k10_partitioned):
    session, fp = k10_partitioned
    results = session.trellis(fp).stable_partitions
    with pytest.raises(ValueError, match="two partitions"):
        PartitionFamily(list(results) + [results[0]])
    left_only = [result for result in results if result.side == "left"]
    family = PartitionFamily(left_only, trellis=session.trellis())
    key = left_only[0].branch_key
    assert family.sides(key) == ["left"]
    with pytest.raises(ValueError, match="no 'right' partition covers"):
        family.result(key, "right")
    assert family.kind == "base"


def test_owner_of_intersection_matches_the_result_table(k10_partitioned):
    session, fp = k10_partitioned
    family = HomotopyPartition.build(session.trellis(fp))
    for result in family:
        for intersection_id, element_id in result.element_of_intersection.items():
            ref = family.owner_of_intersection(intersection_id, result.side)
            assert ref == result.ref(element_id)
    with pytest.raises(ValueError, match="on no"):
        family.owner_of_intersection(10**6, "left")


def test_describe_has_one_line_per_element(k10_partitioned):
    session, fp = k10_partitioned
    family = HomotopyPartition.build(session.trellis(fp))
    lines = family.describe().splitlines()
    elements = sum(len(result.intervals) for result in family)
    assert len(lines) == 1 + len(family) + elements
    assert "homotopy" in lines[0]
    assert repr(family).startswith("<HomotopyPartition")


def test_marks_round_trip_reproduces_every_result(k10_partitioned):
    """The local interval builder agrees with StablePartition's on real data."""
    session, fp = k10_partitioned
    for result in session.trellis(fp).stable_partitions:
        rebuilt = _intervals_from_marks(_marks_of(result))
        assert [
            (iv.lo_id, iv.hi_id, iv.closed_lo, iv.closed_hi)
            for iv in rebuilt
        ] == [
            (iv.lo_id, iv.hi_id, iv.closed_lo, iv.closed_hi)
            for iv in result.intervals
        ]
        for mine, theirs in zip(rebuilt, result.intervals):
            assert mine.lo_cdist == pytest.approx(theirs.lo_cdist)
            assert mine.hi_cdist == pytest.approx(theirs.hi_cdist)


@pytest.mark.slow
def test_marks_round_trip_on_the_nested_fixture(p3_partitioned):
    session, fp3, fp1 = p3_partitioned
    for fp in (fp3, fp1):
        for result in session.trellis(fp).stable_partitions:
            rebuilt = _intervals_from_marks(_marks_of(result))
            assert [
                (iv.lo_id, iv.hi_id, iv.closed_lo, iv.closed_hi) for iv in rebuilt
            ] == [
                (iv.lo_id, iv.hi_id, iv.closed_lo, iv.closed_hi)
                for iv in result.intervals
            ]


def test_kinds_are_distinct():
    assert PartitionFamily.kind == "base"
    assert HomotopyPartition.kind == "homotopy"
    assert IteratedHomotopyPartition.kind == "iterated_homotopy"
    assert issubclass(HomotopyPartition, PartitionFamily)
    assert issubclass(IteratedHomotopyPartition, PartitionFamily)
