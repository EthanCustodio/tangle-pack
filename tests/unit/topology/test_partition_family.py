"""
The partition-family base class and the homotopy family.

Pins: the homotopy family is exactly the trellis's own partition (the
session's family = a direct build is ``tests/facade/test_session_equivalence``);
midpoint ownership is unique
across the full arrangement; the marks round trip reproduces every result
exactly (the refinement rebuilds with the homotopy partition's own
``StablePartition.intervals_from_marks``; the round trip checks that the marks
of a result carry everything ``_build_intervals`` added); an image chain
crossing a stable branch an odd number of times is an incomplete footprint of
the iterated family, recorded and logged once at INFO, never a WARNING (B7).
"""

from __future__ import annotations

import importlib
import logging
from types import SimpleNamespace

import pytest

from helpers.logs import assert_logged
from tanglepack.topology.PartitionFamily import (
    HomotopyPartition,
    IteratedHomotopyPartition,
    PartitionFamily,
    _marks_of,
)
from tanglepack.topology.StablePartition import intervals_from_marks
from tanglepack.topology.TopologyResults import PartitionInterval

#: The module (the package attribute of the same name is the class).
family_module = importlib.import_module("tanglepack.topology.PartitionFamily")
LOGGER = "tanglepack.topology.PartitionFamily"


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


def test_element_at_owns_every_arc_midpoint_exactly_once(k10_partitioned):
    session, fp = k10_partitioned
    family = session.homotopy_partition()
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
    with pytest.raises(ValueError):
        family.element_at(key, "left", end + 10.0)
    with pytest.raises(ValueError):
        family.element_at(key, "left", -1.0)


def test_duplicate_and_missing_results_raise(k10_partitioned):
    session, fp = k10_partitioned
    results = session.trellis(fp).stable_partitions
    with pytest.raises(ValueError):
        PartitionFamily(list(results) + [results[0]])
    left_only = [result for result in results if result.side == "left"]
    family = PartitionFamily(left_only, trellis=session.trellis())
    key = left_only[0].branch_key
    assert family.sides(key) == ["left"]
    with pytest.raises(ValueError):
        family.result(key, "right")
    assert family.kind == "base"


def test_owner_of_intersection_matches_the_result_table(k10_partitioned):
    session, fp = k10_partitioned
    family = HomotopyPartition.build(session.trellis(fp))
    for result in family:
        for intersection_id, element_id in result.element_of_intersection.items():
            ref = family.owner_of_intersection(intersection_id, result.side, result.branch_key)
            assert ref == result.ref(element_id)
    with pytest.raises(ValueError):
        family.owner_of_intersection(10**6, "left", result.branch_key)


def test_describe_has_one_line_per_element(k10_partitioned):
    session, fp = k10_partitioned
    family = HomotopyPartition.build(session.trellis(fp))
    assert family.describe()
    assert repr(family)


def test_marks_round_trip_reproduces_every_result(k10_partitioned):
    """The local interval builder agrees with StablePartition's on real data."""
    session, fp = k10_partitioned
    for result in session.trellis(fp).stable_partitions:
        rebuilt = intervals_from_marks(_marks_of(result))
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


def test_marks_round_trip_on_the_nested_fixture(p3_partitioned):
    session, fp3, fp1 = p3_partitioned
    for fp in (fp3, fp1):
        for result in session.trellis(fp).stable_partitions:
            rebuilt = intervals_from_marks(_marks_of(result))
            assert [
                (iv.lo_id, iv.hi_id, iv.closed_lo, iv.closed_hi) for iv in rebuilt
            ] == [
                (iv.lo_id, iv.hi_id, iv.closed_lo, iv.closed_hi)
                for iv in result.intervals
            ]


def test_an_odd_fold_count_is_an_incomplete_footprint(monkeypatch, caplog):
    """The hole ``(1, 2)``'s image chain ``10 -> 11 -> 12 -> 13 -> 14`` has three
    interior crossings on its base branch: the last is unpaired, the footprint is
    recorded as ``(bridge, branch, 3)`` and logged once at INFO, never WARNING."""
    branch = ("fp", "stable", 0, 0)
    crossings = {
        iid: SimpleNamespace(manifold_b_key=branch, stable_cdist=float(iid), is_anchor=False)
        for iid in (1, 2, 10, 11, 12, 13, 14)
    }
    trellis = SimpleNamespace(intersection=crossings.__getitem__)
    minimal = SimpleNamespace(trellis=trellis, hole_bridge_ids=[(1, 2)])
    hole = PartitionInterval(
        lo_id=1, hi_id=2, lo_cdist=1.0, hi_cdist=2.0, closed_lo=False, closed_hi=False
    )
    homotopy = SimpleNamespace(results={(branch, "right"): SimpleNamespace(intervals=[hole])})
    chain = [(10, 11), (11, 12), (12, 13), (13, 14)]
    monkeypatch.setattr(family_module, "_image_chain", lambda trellis, bid: chain)

    with caplog.at_level(logging.INFO, logger=LOGGER):
        empty = IteratedHomotopyPartition._empty_stretches(minimal, homotopy)
        incomplete = IteratedHomotopyPartition._incomplete_footprints(minimal, homotopy)
    assert incomplete == [((1, 2), branch, 3)]
    assert (14, "base") in empty[10] and (12, "chord") in empty[11]
    assert 13 not in empty, "the unpaired crossing abuts no chord"
    assert_logged(caplog, logging.INFO, LOGGER, count=1)
    assert_logged(caplog, logging.WARNING, LOGGER, count=0)

    family = IteratedHomotopyPartition([], incomplete_footprints=incomplete)
    assert family.incomplete_footprints == incomplete
    assert "incomplete footprint" in family.describe()
