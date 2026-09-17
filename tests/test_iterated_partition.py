"""
The iterated homotopy partition: the homotopy partition cut by image bridges.

The rule, in bridges: every image bridge of the minimal trellis cuts the side
of the stable manifold it lies on, the element under it closed at both its
crossings and the neighbours open; existing hole boundaries stand; the other
side is untouched at an image endpoint.
"""

from __future__ import annotations

import dataclasses
import importlib
import logging

import pytest

from k28_fixture import k28_partitioned  # noqa: F401
from minimal_helpers import build_pieces
from tanglepack.topology.PartitionFamily import IteratedHomotopyPartition
from tanglepack.topology.StablePartition import owns_cdist, row_of_end

#: The module (the package attribute of the same name is the class).
family_module = importlib.import_module("tanglepack.topology.PartitionFamily")


def _flanks(result, boundary_id):
    """(closed_hi of the interval ending at the boundary, closed_lo of the one starting there)."""
    before = [iv for iv in result.intervals if iv.hi_id == boundary_id and iv.lo_id != boundary_id]
    after = [iv for iv in result.intervals if iv.lo_id == boundary_id and iv.hi_id != boundary_id]
    singleton = any(iv.lo_id == iv.hi_id == boundary_id for iv in result.intervals)
    return (
        before[0].closed_hi if before else None,
        after[0].closed_lo if after else None,
        singleton,
    )


def _check_iterated(pieces):
    iterated, homotopy, minimal, full = (
        pieces.iterated, pieces.homotopy, pieces.minimal, pieces.full
    )
    assert iterated.kind == "iterated_homotopy"
    assert iterated.homotopy is homotopy and iterated.minimal is minimal
    assert set(iterated.results) == set(homotopy.results)
    tol = full.registry.cdist_tol
    image_ends_by_key: dict = {}
    for bid in minimal.image_bridge_ids:
        for endpoint, which in zip(bid, ("first", "second")):
            key = (full.intersection(endpoint).manifold_b_key, row_of_end(full, bid, which))
            image_ends_by_key.setdefault(key, set()).add(endpoint)

    for key, result in iterated.results.items():
        parent_result = homotopy.result(*key)
        # Refinement: each element sits inside the parent that owns its midpoint.
        for iv in result.intervals:
            parent = parent_result.element(iv.parent_element_id)
            assert iv.lo_cdist >= parent.lo_cdist - tol
            assert iv.hi_cdist <= parent.hi_cdist + tol
            mid = iv.lo_cdist if iv.lo_cdist == iv.hi_cdist else 0.5 * (iv.lo_cdist + iv.hi_cdist)
            assert owns_cdist(parent, mid, tol)
            assert iv.ref.branch_key == key[0] and iv.ref.side == key[1]
        # Every homotopy boundary survives with the same closedness on both flanks.
        for parent in parent_result.intervals:
            for boundary in (parent.lo_id, parent.hi_id):
                if boundary is None:
                    continue
                assert _flanks(result, boundary) == _flanks(parent_result, boundary)
        # New boundaries are exactly image endpoints on THIS side.
        old = {iv.lo_id for iv in parent_result.intervals} | {
            iv.hi_id for iv in parent_result.intervals
        }
        new = ({iv.lo_id for iv in result.intervals} | {iv.hi_id for iv in result.intervals}) - old
        assert new <= image_ends_by_key.get(key, set())
        # Exactly one owner for every crossing of the full branch.
        branch = full.branch(key[0])
        for iid in branch.intersection_ids:
            cdist = full.intersection(iid).stable_cdist
            owners = [iv for iv in result.intervals if owns_cdist(iv, cdist, tol)]
            assert len(owners) == 1
            assert result.element_of_intersection[iid] == owners[0].element_id
        assert len(result.intervals) >= len(parent_result.intervals)

    # Under each same-branch image bridge, on its side: closed at both ends
    # unless the end was already a boundary (then unchanged), open beyond.
    for bid in minimal.image_bridge_ids:
        a, b = bid
        branch_a = full.intersection(a).manifold_b_key
        if branch_a != full.intersection(b).manifold_b_key:
            continue
        row = row_of_end(full, bid, "first")
        if row != row_of_end(full, bid, "second"):
            continue
        result = iterated.result(branch_a, row)
        parent_result = homotopy.result(branch_a, row)
        old = {iv.lo_id for iv in parent_result.intervals} | {
            iv.hi_id for iv in parent_result.intervals
        }
        near, far = sorted(bid, key=lambda i: full.intersection(i).stable_cdist)
        before_near, after_near, _ = _flanks(result, near)
        before_far, after_far, _ = _flanks(result, far)
        if near not in old:
            assert after_near is True and before_near is False
        if far not in old:
            assert before_far is True and after_far is False
        # Every element under the bridge records an image bridge covering it.
        near_c = full.intersection(near).stable_cdist
        far_c = full.intersection(far).stable_cdist
        under = [
            iv for iv in result.intervals
            if iv.lo_cdist >= near_c - tol and iv.hi_cdist <= far_c + tol
        ]
        assert under
        for iv in under:
            assert iv.cut_by is not None and minimal.is_image_bridge(iv.cut_by)
            cover_near, cover_far = sorted(
                iv.cut_by, key=lambda i: full.intersection(i).stable_cdist
            )
            assert full.intersection(cover_near).stable_cdist <= iv.lo_cdist + tol
            assert full.intersection(cover_far).stable_cdist >= iv.hi_cdist - tol
        if under[0].lo_id == near and under[-1].hi_id == far:
            assert (near in old or under[0].closed_lo) and (far in old or under[-1].closed_hi)
        # The other side is not cut at an image endpoint by THIS bridge.
        other = iterated.result(branch_a, "right" if row == "left" else "left")
        other_parent = homotopy.result(branch_a, "right" if row == "left" else "left")
        other_old = {iv.lo_id for iv in other_parent.intervals} | {
            iv.hi_id for iv in other_parent.intervals
        }
        other_new = (
            {iv.lo_id for iv in other.intervals} | {iv.hi_id for iv in other.intervals}
        ) - other_old
        for endpoint in bid:
            if endpoint in other_new:
                # Only legitimate when another image bridge on that side ends there
                # (a crossing shared by two consecutive pieces of one image).
                assert endpoint in image_ends_by_key.get(
                    (branch_a, "right" if row == "left" else "left"), set()
                )

    # Every arc midpoint of the minimal arrangement has exactly one owner per side.
    for face in minimal.arrangement.faces:
        for arc in face.arcs:
            if arc.kind != "stable":
                continue
            mid = 0.5 * (
                full.intersection(arc.lo_id).stable_cdist
                + full.intersection(arc.hi_id).stable_cdist
            )
            for side in ("left", "right"):
                iterated.interval_at(arc.branch_key, side, mid)

    applied = [cut for cut in iterated.cuts if cut.opened is not None]
    if applied:
        assert iterated.signature() != homotopy.signature()
        assert sum(len(r.intervals) for r in iterated) > sum(len(r.intervals) for r in homotopy)


def test_k10_iterated_partition_invariants(k10_partitioned):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    _check_iterated(pieces)
    assert any(cut.opened for cut in pieces.iterated.cuts), "k=10 cuts something"
    assert any(iv.cut_by is not None for r in pieces.iterated for iv in r.intervals)


@pytest.mark.slow
def test_p3_iterated_partition_invariants(p3_partitioned):
    session, fp3, fp1 = p3_partitioned
    _check_iterated(build_pieces(session, [fp3, fp1]))


@pytest.mark.slow
def test_k28_iterated_partition_invariants(k28_partitioned):
    session, fp = k28_partitioned
    _check_iterated(build_pieces(session, [fp]))


def test_no_image_bridges_reproduces_the_homotopy_family(k10_partitioned):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    bare = dataclasses.replace(pieces.minimal, image_bridge_ids=[], image_chains={})
    iterated = IteratedHomotopyPartition.from_minimal(bare, pieces.homotopy)
    assert iterated.signature() == pieces.homotopy.signature()
    assert not iterated.cuts
    for result in iterated:
        for iv in result.intervals:
            assert iv.cut_by is None and iv.parent_element_id == iv.element_id


def test_a_row_invariant_violation_is_warned_and_skipped(
    k10_partitioned, monkeypatch, caplog
):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    victim = pieces.minimal.image_bridge_ids[0]
    original = family_module.row_of_end

    def twisted(trellis, bridge_id, endpoint):
        row = original(trellis, bridge_id, endpoint)
        if bridge_id == victim and endpoint == "second":
            return "right" if row == "left" else "left"
        return row

    monkeypatch.setattr(family_module, "row_of_end", twisted)
    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.PartitionFamily"):
        iterated = IteratedHomotopyPartition.from_minimal(pieces.minimal, pieces.homotopy)
    skipped = [cut for cut in iterated.cuts if cut.bridge_id == victim]
    assert skipped and all(cut.reason == "row invariant" for cut in skipped)
    assert any("row invariant" in record.message for record in caplog.records)
    assert sum(len(r.intervals) for r in iterated) < sum(
        len(r.intervals) for r in pieces.iterated
    )


def test_describe_mentions_parents_and_cuts(k10_partitioned):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    text = pieces.iterated.describe()
    assert "iterated_homotopy" in text and "parent #" in text and "cut by" in text
    assert repr(pieces.iterated.cuts[0]).startswith("Cut(")
