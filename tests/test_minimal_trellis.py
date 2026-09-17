"""
The minimal trellis: hole bridges, images of the active ones, nothing else.

Every test pins an invariant, never an id or a count of image bridges (the
inertness definition is changing on a concurrent branch, so the image set on
k=2.8 may move).
"""

from __future__ import annotations

import logging

import pytest

from k28_fixture import k28_partitioned  # noqa: F401
from minimal_helpers import build_pieces
from tanglepack.topology import BridgeClass
from tanglepack.topology.MinimalTrellis import image_chain, minimal_trellis


def _kept_pairs(arrangement):
    return {
        frozenset((arc.lo_id, arc.hi_id))
        for face in arrangement.faces
        for arc in face.arcs
        if arc.kind == "unstable"
    }


def _check_minimal(pieces):
    minimal, full, holes, table = pieces.minimal, pieces.full, pieces.holes, pieces.table
    hole_set = {frozenset(b) for b in minimal.hole_bridge_ids}
    # Every hole's bridge is kept as a hole bridge, and nothing else is one.
    for hole in holes:
        if hole.bounding_ids is None:
            continue
        bridge = full.bridge_between(*hole.bounding_ids)
        if bridge is not None and bridge.id is not None:
            assert frozenset(bridge.id) in hole_set
    assert len(hole_set) == len(minimal.hole_bridge_ids)
    # Inert and unclassed hole bridges are kept but never mapped.
    for bid in minimal.inert_hole_bridge_ids + minimal.unclassed_hole_bridge_ids:
        assert minimal.is_hole_bridge(bid)
        assert bid not in minimal.image_chains
    for bid in minimal.inert_hole_bridge_ids:
        assert not table.entry_of(bid).active
    # Every mapped hole bridge is active and its chain is image_chain's, kept
    # pair by kept pair (skipped ones excepted).
    for bid, chain in minimal.image_chains.items():
        assert table.entry_of(bid).active
        expected = image_chain(full, bid)
        assert expected is not None
        skipped = {frozenset(p) for parent, p in minimal.skipped_pairs if parent == bid}
        assert [frozenset(p) for p in chain] == [
            frozenset(p) for p in expected if frozenset(p) not in skipped
        ]
    for bid in minimal.unmapped:
        assert image_chain(full, bid) is None
    # Image bridges are chain pairs that are not hole bridges.
    image_set = {frozenset(b) for b in minimal.image_bridge_ids}
    assert not (image_set & hole_set)
    for bid in minimal.image_bridge_ids:
        parents = minimal.parents_of(bid)
        assert parents and all(minimal.is_hole_bridge(p) for p in parents)
        assert minimal.is_image_bridge(bid) and minimal.is_kept(bid)
    # Skipped pairs have no Bridge object.
    for _parent, pair in minimal.skipped_pairs:
        assert full.bridge_between(*pair) is None
    # Kept bridges are real, non-partial, and their endpoints are kept nodes.
    kept_ids = {frozenset(b.id) for b in minimal.kept}
    assert kept_ids == hole_set | image_set
    assert all(not b.partial for b in minimal.kept)
    dropped = {frozenset(b) for b in minimal.dropped_bridge_ids}
    all_ids = {frozenset(b.id) for b in full.bridges if b.id is not None}
    assert dropped == all_ids - kept_ids
    endpoints = {e for b in minimal.kept for e in b.id}
    tol = full.registry.cdist_tol
    for branch in full.stable_branches:
        kept = minimal.kept_ids(branch.key)
        expected = {
            i for i in endpoints if full.intersection(i).manifold_b_key == branch.key
        }
        expected |= {
            i for i in branch.intersection_ids
            if full.intersection(i).stable_cdist <= tol
        }
        if len(branch):
            expected.add(branch.intersection_ids[-1])
        assert set(kept) == expected
        assert set(kept) <= set(branch.intersection_ids)
        cdists = [full.intersection(i).stable_cdist for i in kept]
        assert cdists == sorted(cdists)
    # The sparse snapshot shares the parent's registry and manifolds.
    assert minimal.sparse.registry is full.registry
    assert minimal.sparse.manifolds is full.manifolds
    assert minimal.sparse._built_generation == full._built_generation
    # The sparse arrangement carries exactly the kept bridges as unstable arcs,
    # closes up (Euler), and stubs no unstable ray.
    arrangement = minimal.arrangement
    assert arrangement.sparse
    assert _kept_pairs(arrangement) == kept_ids
    assert arrangement.euler_characteristic == 2 * arrangement.component_count
    for node in arrangement._nodes.values():
        if node.virtual:
            continue
        for slot in ("u+", "u-"):
            if slot in node.slots:
                head = arrangement._half_edges[node.slots[slot]].head
                assert not arrangement._nodes[head].virtual
    for face in arrangement.faces:
        for arc in face.arcs:
            if arc.kind == "stable":
                kept = minimal.kept_ids(arc.branch_key)
                assert arc.lo_id in kept and arc.hi_id in kept
    assert minimal.arrangement is arrangement
    assert "MinimalTrellis" in minimal.summary() and minimal.summary() in repr(minimal)
    assert minimal.describe().count("\n") >= len(minimal.kept_bridge_ids)


def test_k10_minimal_trellis_invariants(k10_partitioned):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    _check_minimal(pieces)
    assert pieces.minimal.hole_bridge_ids
    assert pieces.minimal.image_bridge_ids, "k=10 maps at least one hole bridge forward"


def test_k10_minimal_trellis_drops_something(k10_partitioned):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    assert len(pieces.minimal.kept) < len(
        [b for b in pieces.full.bridges if b.id is not None]
    )


@pytest.mark.slow
def test_p3_minimal_trellis_invariants(p3_partitioned):
    session, fp3, fp1 = p3_partitioned
    _check_minimal(build_pieces(session, [fp3, fp1]))


@pytest.mark.slow
def test_k28_minimal_trellis_invariants(k28_partitioned):
    session, fp = k28_partitioned
    pieces = build_pieces(session, [fp])
    _check_minimal(pieces)
    assert pieces.minimal.image_bridge_ids or pieces.minimal.unmapped or pieces.minimal.skipped_pairs


def test_image_chain_matches_the_bridge_class_copy_when_present(k10_partitioned):
    theirs = getattr(BridgeClass, "_image_chain", None)
    if theirs is None:
        pytest.skip("BridgeClass._image_chain has not landed on this branch")
    session, fp = k10_partitioned
    full = session.trellis()
    for bridge in full.bridges:
        if bridge.id is not None:
            assert image_chain(full, bridge.id) == theirs(full, bridge.id)


def test_a_pair_with_no_bridge_object_is_skipped_with_a_warning(
    k10_partitioned, monkeypatch, caplog
):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    full = pieces.full
    victim = pieces.minimal.image_bridge_ids[0]
    original = full.bridge_between

    def hide(a, b):
        if frozenset((a, b)) == frozenset(victim):
            return None
        return original(a, b)

    monkeypatch.setattr(full, "bridge_between", hide)
    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.MinimalTrellis"):
        minimal = minimal_trellis(full, pieces.holes, pieces.table)
    assert any(frozenset(p) == frozenset(victim) for _parent, p in minimal.skipped_pairs)
    assert not minimal.is_kept(victim)
    assert any("no Bridge object" in record.message for record in caplog.records)


def test_an_unclassed_hole_bridge_is_kept_but_not_mapped(
    k10_partitioned, monkeypatch, caplog
):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    victim = next(iter(pieces.minimal.image_chains))
    original = pieces.table.entry_of

    def forget(bid):
        if frozenset(bid) == frozenset(victim):
            raise KeyError(bid)
        return original(bid)

    monkeypatch.setattr(pieces.table, "entry_of", forget)
    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.MinimalTrellis"):
        minimal = minimal_trellis(pieces.full, pieces.holes, pieces.table)
    assert victim in minimal.unclassed_hole_bridge_ids
    assert minimal.is_hole_bridge(victim) and victim not in minimal.image_chains
    assert any("no bridge class" in record.message for record in caplog.records)


def test_a_hole_with_no_spanning_bridge_contributes_nothing(k10_partitioned, caplog):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    from tanglepack.topology.TopologyResults import Hole

    ghost = Hole(coords=(0.0, 0.0), near_intersection_id=0, bounding_ids=(10**6, 10**6 + 1))
    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.MinimalTrellis"):
        minimal = minimal_trellis(pieces.full, pieces.holes + [ghost], pieces.table)
    assert [frozenset(b) for b in minimal.hole_bridge_ids] == [
        frozenset(b) for b in pieces.minimal.hole_bridge_ids
    ]
    assert any("no bridge of" in record.message for record in caplog.records)
