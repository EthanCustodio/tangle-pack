"""Bridge identity (``BridgeId``) and the derived genealogy built on top of it.

A bridge is the piece of unstable manifold between two consecutive crossings, so
its topological identity is exactly the ordered pair of those crossings' registry
ids, taken in the unstable dynamical direction. These tests pin that identity, the
reverse index built from it, the derived image/preimage genealogy that replaces the
stored ``parent``/``children`` links, and the segment-ownership release that
``clear_bridges`` now performs.
"""

from __future__ import annotations

import gc
import logging
import weakref

import pytest

from helpers.logs import assert_logged


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _full_bridges(workbench):
    """Every registered bridge that is bounded by two crossings (has an id)."""
    return [b for b in workbench.bridges if not b.partial]


def _iterable_interior_bridge(session, fp, zone):
    """An un-iterated interior bridge of ``fp`` inside ``zone``, or None."""
    for bridge in session.workbench.uniiterated_bridges:
        if bridge.fixed_point is not fp or bridge.partial:
            continue
        point = session._bridge_test_point(bridge)
        if point is not None and zone.contains_point(point):
            return bridge
    return None


# --------------------------------------------------------------------------- #
# Bridge.id
# --------------------------------------------------------------------------- #
def test_partial_pieces_are_held_outside_the_id_registry(henon_tangle_with_bridges):
    """A partial image piece is a bridge object but carries no identity."""
    workbench, _fp = henon_tangle_with_bridges

    partials = []
    for parent in list(workbench.uniiterated_bridges):
        if parent.iterated or parent.partial:
            continue
        partials.extend(c for c in workbench.iterate_bridge(parent) if c.partial)
        if partials:
            break
    assert partials, "iterating the k=10 bridges must produce a partial piece"

    registered = workbench.bridges
    identified = {id(workbench.bridge(b.id)) for b in registered if b.id is not None}
    for piece in partials:
        assert piece.id is None
        assert any(piece is b for b in registered), (
            "a partial is still a registered bridge"
        )
        assert id(piece) not in identified, "no id may look a partial up"

    # Every id-carrying bridge is the one copy its id looks up; only partials lack one.
    for bridge in registered:
        assert (bridge.id is None) == bridge.partial
        if bridge.id is not None:
            assert workbench.bridge(bridge.id) is bridge


# --------------------------------------------------------------------------- #
# registry: bridge(), bridges_at(), single copy
# --------------------------------------------------------------------------- #
def test_bridge_lookup_by_id(henon_tangle_with_bridges):
    """``workbench.bridge(bid)`` returns the one registered copy; KeyError else."""
    workbench, _fp = henon_tangle_with_bridges
    for bridge in _full_bridges(workbench):
        assert workbench.bridge(bridge.id) is bridge
    with pytest.raises(KeyError):
        workbench.bridge((-1, -2))


def test_standalone_collaborators_agree_with_the_workbench(henon_tangle_with_bridges):
    """A ``BridgeIterator`` / ``IterateInference`` built over a workbench reads
    that workbench: its genealogy matches the workbench's delegates, and a
    second inference pass over an already-inferred table adds nothing."""
    from tanglepack.numerics.BridgeIterator import BridgeIterator
    from tanglepack.numerics.IterateInference import IterateInference

    workbench, _fp = henon_tangle_with_bridges
    iterator = BridgeIterator(workbench)
    inference = IterateInference(workbench)
    assert iterator.workbench is workbench and inference.workbench is workbench

    bridges = _full_bridges(workbench)
    assert bridges
    for bridge in bridges:
        assert iterator.image_bridges(bridge.id) == workbench.image_bridges(bridge.id)
        assert iterator.preimage_bridges(bridge.id) == workbench.preimage_bridges(bridge.id)
    assert inference.infer_iterate_table() == 0


# --------------------------------------------------------------------------- #
# rebuild_bridges
# --------------------------------------------------------------------------- #
def test_bridge_ids_survive_rebuild_with_metadata(henon_tangle_with_bridges):
    """Ids and the ``iterated`` flag are carried across ``rebuild_bridges``."""
    workbench, fp = henon_tangle_with_bridges

    target = next(b for b in workbench.uniiterated_bridges if not b.partial)
    iterated_id = target.id
    workbench.iterate_bridge(target)
    assert workbench.bridge(iterated_id).iterated

    before = {b.id for b in _full_bridges(workbench)}
    workbench.rebuild_bridges(fp)
    after = {b.id for b in _full_bridges(workbench)}

    # Recutting the indexed manifolds reproduces every bridge cut from them; the
    # image bridges produced by iterate_bridge are not recut (they are not cut out
    # of an indexed manifold), so the recut set is a subset of what was held.
    assert after
    assert after <= before
    assert iterated_id in after
    assert workbench.bridge(iterated_id).iterated, (
        "rebuild_bridges must carry the iterated flag across by id"
    )


# --------------------------------------------------------------------------- #
# derived genealogy
# --------------------------------------------------------------------------- #
def test_image_bridges_matches_iterate_bridge(henon_p3_session):
    """``image_bridges`` reproduces the children ``iterate_bridge`` returned."""
    session, fp3, _fp1, zone = henon_p3_session
    workbench = session.workbench

    bridge = _iterable_interior_bridge(session, fp3, zone)
    assert bridge is not None, "the fixture must offer an interior bridge to iterate"
    bid = bridge.id

    children = workbench.iterate_bridge(bridge)
    child_ids = [c.id for c in children if not c.partial]
    assert child_ids, "the iterate must produce at least one full child bridge"

    assert workbench.image_bridges(bid) == child_ids
    for child_id in child_ids:
        # A crossing INTERIOR to the image has a preimage that is a crossing of the
        # stable manifold too far out to have been computed, so its n=-1 entry can
        # legitimately be missing; when it is there, it must point back.
        preimage = workbench.preimage_bridges(child_id)
        if preimage is not None:
            assert bid in preimage


def test_image_bridges_is_none_without_a_registered_iterate(
    henon_tangle_with_bridges,
):
    """No n=1 entry for an endpoint means no derivable image."""
    workbench, _fp = henon_tangle_with_bridges
    table = workbench.intersection_registry.iterate_table

    orphan = next(
        (
            b.id
            for b in _full_bridges(workbench)
            if table[b.first_intersection, 1] is None
            or table[b.second_intersection, 1] is None
        ),
        None,
    )
    assert orphan is not None, "the fixture must contain an outermost bridge"
    assert workbench.image_bridges(orphan) is None


def _walk_bridge_chain(workbench, start_id, end_id, manifold_key):
    """Consecutive bridges from ``start_id`` to ``end_id`` along one unstable branch.

    An independent derivation of a span's bridges: it follows the shared-endpoint
    adjacency (``bridges_at``) instead of filtering on canonical distance, so it
    catches an off-by-one in the interval bounds ``image_bridges`` uses.
    """
    chain: list[tuple[int, int]] = []
    current = start_id
    while current != end_id:
        onward = [
            bid
            for bid in workbench.bridges_at(current)
            if bid[0] == current and workbench.bridge(bid).manifold_key == manifold_key
        ]
        if not onward:
            return None
        assert len(onward) == 1, f"two bridges leave crossing {current} forward"
        chain.append(onward[0])
        current = onward[0][1]
    return chain


def test_preimage_bridges_equals_the_backward_chain(henon_p3_session):
    """``preimage_bridges`` is the bridge chain between the endpoints' preimages."""
    session, _fp3, _fp1, _zone = henon_p3_session
    workbench = session.workbench
    registry = workbench.intersection_registry
    table = registry.iterate_table

    checked = 0
    for bridge in _full_bridges(workbench):
        first, second = bridge.id
        back = (table[first, -1], table[second, -1])
        if any(b is None for b in back):
            assert workbench.preimage_bridges(bridge.id) is None
            continue

        # Order the two preimages by unstable cdist; the chain runs low to high.
        lo, hi = sorted(back, key=lambda i: registry[i].unstable_cdist)
        source_key = bridge.fixed_point.advance_key(bridge.manifold_key, -1)
        expected = _walk_bridge_chain(workbench, lo, hi, source_key)
        if expected is None:
            continue

        assert workbench.preimage_bridges(bridge.id) == expected
        assert workbench.image_bridges(bridge.id, -1) == expected
        checked += 1

    assert checked > 0, "the fixture must offer a walkable backward chain"


# --------------------------------------------------------------------------- #
# clear_bridges releases the segments it claimed
# --------------------------------------------------------------------------- #
def test_clear_bridges_releases_segment_ownership(grown_both):
    """``clear_bridges`` lets the discarded bridges go, and a recut reproduces them.

    A bridge shares its segments with the unstable manifold it was cut from, and
    the Tangle index records that claim. The bug was that ``clear_bridges`` left
    the claim in place, so every discarded bridge stayed reachable from the index
    and was never collected. Observed publicly: after the clear no weak reference
    to an old bridge survives garbage collection, and a fresh ``create_bridges``
    cuts the same id set again.
    """
    workbench, fp = grown_both
    workbench.compute_intersections([fp])
    workbench.trim_stable_manifolds(fp)
    bridges = workbench.create_bridges(fp)
    assert bridges
    ids_before = {b.id for b in bridges if b.id is not None}
    refs = [weakref.ref(b) for b in bridges]
    del bridges

    workbench.clear_bridges()
    gc.collect()

    assert not workbench.bridges
    leaked = [ref for ref in refs if ref() is not None]
    assert not leaked, f"{len(leaked)} discarded bridge(s) are still reachable"

    recut = workbench.create_bridges(fp)
    assert {b.id for b in recut if b.id is not None} == ids_before


# --------------------------------------------------------------------------- #
# session delegation
# --------------------------------------------------------------------------- #
def test_session_exposes_the_identity_api(henon_p3_session):
    """The facade reaches the workbench's bridge-identity methods."""
    session, _fp3, _fp1, _zone = henon_p3_session
    bid = next(b.id for b in session.workbench.bridges if not b.partial)
    assert session.bridge(bid) is session.workbench.bridge(bid)
    assert session.bridges_at(bid[0]) == session.workbench.bridges_at(bid[0])
    assert session.image_bridges(bid) == session.workbench.image_bridges(bid)
    assert session.preimage_bridges(bid) == session.workbench.preimage_bridges(bid)


# --------------------------------------------------------------------------- #
# metadata carry is gated on the registry the ids were captured against
# --------------------------------------------------------------------------- #
def test_rebuild_drops_metadata_after_an_unpreserved_recompute(
    henon_tangle_with_bridges, caplog
):
    """A renumbering recompute invalidates the ids, so no flag may be carried.

    ``BridgeId`` is a pair of registry ids, and ``compute_intersections`` without
    ``preserve_ids`` renumbers the registry from zero. The old pairs then name
    unrelated crossings, so carrying ``iterated`` across by id would silently mark a
    bridge that has never been iterated (or lose a flag that should have survived).
    Whether a stale pair happens to collide with a live one is luck, so the drop
    itself is what is pinned here, not just its absence of symptoms.
    """
    workbench, fp = henon_tangle_with_bridges

    target = next(b for b in workbench.uniiterated_bridges if not b.partial)
    workbench.iterate_bridge(target)
    assert any(b.iterated for b in workbench.bridges)

    workbench.compute_intersections([fp])  # preserve_ids=False: ids renumbered

    with caplog.at_level(
        logging.DEBUG, logger="tanglepack.numerics.TangleWorkbench"
    ):
        workbench.rebuild_bridges(fp)

    # rebuild_bridges reports the dropped carry (level and logger only)
    assert_logged(caplog, logging.DEBUG, "tanglepack.numerics.TangleWorkbench")

    wrongly_marked = [b.id for b in workbench.bridges if b.iterated]
    assert not wrongly_marked, (
        f"bridges {wrongly_marked} were marked iterated from stale ids"
    )


def test_rebuild_keeps_metadata_after_a_preserving_recompute(
    henon_tangle_with_bridges,
):
    """With ``preserve_ids=True`` the ids still mean the same crossings."""
    workbench, fp = henon_tangle_with_bridges

    target = next(b for b in workbench.uniiterated_bridges if not b.partial)
    iterated_id = target.id
    workbench.iterate_bridge(target)

    workbench.compute_intersections([fp], preserve_ids=True)
    workbench.rebuild_bridges(fp)

    assert workbench.bridge(iterated_id).iterated


# --------------------------------------------------------------------------- #
# inversion saddle (k_value == 2 * period)
# --------------------------------------------------------------------------- #
def test_iterating_an_inversion_bridge_lands_on_the_other_branch(henon_inversion):
    """``iterate_bridge`` on an inversion point advances the BRANCH index too.

    One map step sends a bridge of branch ``(0, 0)`` to ``(0, 1)`` and back
    (:meth:`FixedPoint.advance_key`), which is exactly what ``image_bridges``
    relies on to name the image branch: every child of an iterated bridge lies
    on the advanced key, and the DERIVED image contains every identified child
    (it can contain more: the image arc may also cover bridges already there).

    The id, ``bridges_at`` and single-copy laws on this saddle run in
    ``tests/invariants/test_law_bridges.py`` (case ``inversion``); a bridge
    between the saddle's two coincident anchors (the known two-anchors issue)
    is skipped here.
    """
    workbench, fp = henon_inversion
    bridges = workbench.create_bridges(fp)
    registry = workbench.intersection_registry
    proper = [
        b
        for b in bridges
        if not b.partial
        and registry[b.id[0]].unstable_cdist < registry[b.id[1]].unstable_cdist
    ]
    assert proper, "the inversion fixture must produce a non-degenerate bridge"

    exact = 0
    checked = 0
    for bridge in proper:
        if bridge.iterated:
            continue
        image_key = fp.advance_key(bridge.manifold_key, 1)
        assert image_key != bridge.manifold_key, "one map step must leave the branch"

        children = workbench.iterate_bridge(bridge)
        assert children
        for child in children:
            assert child.manifold_key == image_key

        derived = workbench.image_bridges(bridge.id)
        cut_ids = [c.id for c in children if not c.partial]
        if not cut_ids:
            continue
        assert derived is not None, (
            f"bridge {bridge.id} cut children {cut_ids} but has no derivable image"
        )
        for child_id in cut_ids:
            assert child_id in derived
        checked += 1
        if derived == cut_ids:
            exact += 1

    assert checked > 0, "no inversion bridge produced an identified child"
    assert exact > 0, "the derived image must reproduce at least one cut exactly"
