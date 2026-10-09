"""
Exact ownership of registry cdists (stress campaign E11, 2026-10-07).

Every partition boundary is a registry float, so a registered crossing is owned
with no tolerance, however close its neighbours sit, and since WP11 (C15) an
ESTIMATED cdist is located exactly too. Pure kernels over a hand-built registry: crossings closer than
the tolerance, a sliver of width ``1.5 * tol`` beside a closed singleton (the
shape that tripped the iterated child-inside-parent assertion at 9 blasts), and
the anchor detected by its exact cdist.
"""

from __future__ import annotations

import logging
from types import SimpleNamespace

import pytest

from tanglepack.topology.PartitionFamily import (
    PartitionFamily,
    _intervals_from_marks,
    _parent_interval,
)
from tanglepack.topology.StablePartition import _build_intervals, _element_of_intersection
from tanglepack.topology.TopologyResults import StablePartitionResult

TOL = 1e-6
BRANCH = ("stable", 0, 0)


class _Registry:
    """The registry slice the ownership kernels read: crossings by id, and the tolerance."""

    cdist_tol = TOL

    def __init__(self, cdists: dict[int, float]) -> None:
        self.cdists = cdists

    def intersection(self, iid: int) -> SimpleNamespace:
        return SimpleNamespace(stable_cdist=self.cdists[iid], manifold_b_key=BRANCH)


def _trellis(cdists: dict[int, float]) -> SimpleNamespace:
    registry = _Registry(cdists)
    return SimpleNamespace(registry=registry, intersection=registry.intersection)


class _Branch:
    """The branch slice the kernels read: its key and crossings, outermost first."""

    key = BRANCH

    def __init__(self, cdists: dict[int, float]) -> None:
        self.intersection_ids = list(cdists)
        self._outermost_first = sorted(cdists, key=cdists.get, reverse=True)

    def __len__(self) -> int:
        return len(self.intersection_ids)

    def ordered_ids(self, toward_anchor: bool = False) -> list[int]:
        return self._outermost_first if toward_anchor else self._outermost_first[::-1]


def _result(cdists, boundaries, open_outward=(), open_anchorward=()) -> StablePartitionResult:
    """A stamped result from marks, its ownership table built by the production kernel."""
    intervals = _intervals_from_marks(
        ({bid: cdists.get(bid, 0.0) for bid in boundaries}, set(open_outward), set(open_anchorward))
    )
    for element_id, interval in enumerate(intervals):
        interval.element_id, interval.branch_key, interval.side = element_id, BRANCH, "right"
    owners = _element_of_intersection(_trellis(cdists), _Branch(cdists), intervals)
    return StablePartitionResult(BRANCH, "right", intervals, owners)


def _sliver_beside_a_singleton():
    """``[0, 1) [1, 1] (1, s]`` with the sliver ``s - 1 = 1.5 tol``: its midpoint is
    within the tolerance of the closed singleton."""
    cdists = {1: 1.0, 2: 1.0 + 1.5 * TOL}
    return cdists, _result(cdists, [None, 1, 2], open_outward=[1], open_anchorward=[1])


def test_interval_at_gives_a_sliver_beside_a_singleton_its_own_midpoint():
    cdists, result = _sliver_beside_a_singleton()
    family = PartitionFamily([result], trellis=_trellis(cdists))
    sliver = result.intervals[-1]
    assert (sliver.lo_id, sliver.hi_id) == (1, 2) and not sliver.closed_lo
    mid = 0.5 * (sliver.lo_cdist + sliver.hi_cdist)
    assert family.interval_at(BRANCH, "right", mid) is sliver
    assert family.interval_at(BRANCH, "right", 1.0).lo_id == family.interval_at(
        BRANCH, "right", 1.0
    ).hi_id == 1


def test_interval_at_is_exact_for_an_estimated_cdist():
    """No absolute slack (C15, WP11): a cdist just past the closed outer end has no
    owner, however close it sits; the outer end itself is owned."""
    cdists, result = _sliver_beside_a_singleton()
    family = PartitionFamily([result], trellis=_trellis(cdists))
    assert family.interval_at(BRANCH, "right", cdists[2]) is result.intervals[-1]
    for estimated in (cdists[2] + 0.5 * TOL, cdists[2] + 10 * TOL):
        with pytest.raises(ValueError):
            family.interval_at(BRANCH, "right", estimated)


def test_crossings_closer_than_the_tolerance_each_have_one_owner():
    # [0, a] (a, b) [b, c]: a hole between two boundaries tol/2 apart, an
    # interior crossing d tol/5 inside the hole, the outermost crossing c.
    a, b = 1.0, 1.0 + 0.5 * TOL
    cdists = {1: a, 2: b, 3: a + 0.2 * TOL, 4: 2.0}
    result = _result(cdists, [None, 1, 2, 4], open_outward=[1], open_anchorward=[2])
    hole = next(iv for iv in result.intervals if (iv.lo_id, iv.hi_id) == (1, 2))
    by_id = {iv.element_id: iv for iv in result.intervals}
    owners = result.element_of_intersection
    assert set(owners) == set(cdists)
    assert by_id[owners[1]].hi_id == 1 and by_id[owners[1]].closed_hi
    assert by_id[owners[2]].lo_id == 2 and by_id[owners[2]].closed_lo
    assert owners[3] == hole.element_id
    assert by_id[owners[4]].hi_id == 4


def test_a_boundary_off_its_registry_cdist_is_refused():
    cdists = {1: 1.0, 2: 2.0}
    intervals = _intervals_from_marks(({None: 0.0, 1: 1.0 + TOL / 4, 2: 2.0}, {1}, set()))
    for element_id, interval in enumerate(intervals):
        interval.element_id = element_id
    with pytest.raises(AssertionError):
        _element_of_intersection(_trellis(cdists), _Branch(cdists), intervals)


def test_iterated_parent_is_found_by_exact_containment_and_singletons_by_id():
    # Homotopy [0, 1) [1, 1] (1, 2]; the iterated cut pinches a NEW crossing 3 a
    # sliver of 1.5 tol past the singleton: (1, 3) [3, 3] (3, 2].
    cdists = {1: 1.0, 2: 2.0, 3: 1.0 + 1.5 * TOL}
    homotopy = _result(cdists, [None, 1, 2], open_outward=[1], open_anchorward=[1])
    iterated = _result(
        cdists, [None, 1, 2, 3], open_outward=[1, 3], open_anchorward=[1, 3]
    )
    spans = [(iv.lo_id, iv.hi_id) for iv in iterated.intervals]
    assert spans == [(None, 1), (1, 1), (1, 3), (3, 3), (3, 2)]
    parents = [
        (_parent_interval(child, homotopy).lo_id, _parent_interval(child, homotopy).hi_id)
        for child in iterated.intervals
    ]
    assert parents == [(None, 1), (1, 1), (1, 2), (1, 2), (1, 2)]


def test_the_anchor_is_its_exact_cdist(caplog):
    """A boundary at exactly 0.0 is the anchor; any other gets a synthetic anchor
    below it, however close it sits (no absolute slack, C15, WP11)."""
    caplog.set_level(logging.INFO, logger="tanglepack.topology.StablePartition")
    anchored = {0: 0.0, 1: 1.0}
    intervals = _build_intervals(_trellis(anchored), _Branch(anchored), ({0: 0.0}, set(), set()))
    assert [iv.lo_id for iv in intervals] == [0]
    assert not caplog.records

    unanchored = {1: 1.0}
    intervals = _build_intervals(_trellis(unanchored), _Branch(unanchored), ({}, set(), set()))
    assert [(iv.lo_id, iv.lo_cdist) for iv in intervals] == [(None, 0.0)]

    near = {5: TOL / 2, 1: 1.0}
    intervals = _build_intervals(_trellis(near), _Branch(near), ({5: TOL / 2}, set(), set()))
    assert [(iv.lo_id, iv.lo_cdist) for iv in intervals] == [(None, 0.0), (5, TOL / 2)]
    assert not caplog.records
