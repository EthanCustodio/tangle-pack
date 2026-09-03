"""Bridge identity (``BridgeId``) and the derived genealogy built on top of it.

A bridge is the piece of unstable manifold between two consecutive crossings, so
its topological identity is exactly the ordered pair of those crossings' registry
ids, taken in the unstable dynamical direction. These tests pin that identity, the
reverse index built from it, the derived image/preimage genealogy that replaces the
stored ``parent``/``children`` links, and the segment-ownership release that
``clear_bridges`` now performs.
"""

from __future__ import annotations

import logging

import pytest


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _full_bridges(workbench):
    """Every registered bridge that is bounded by two crossings (has an id)."""
    return [b for b in workbench.bridges if not b.partial]


def _check_id_ordering(workbench) -> int:
    """Assert every non-partial bridge's id is its endpoints in cdist order."""
    registry = workbench.intersection_registry
    checked = 0
    for bridge in _full_bridges(workbench):
        bid = bridge.id
        assert bid == (bridge.first_intersection, bridge.second_intersection)
        first, second = bid
        assert registry[first].unstable_cdist < registry[second].unstable_cdist, (
            f"bridge id {bid} is not ordered by increasing unstable cdist"
        )
        checked += 1
    return checked


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
def test_bridge_id_is_endpoints_in_unstable_order(henon_tangle_with_bridges):
    """``Bridge.id`` is ``(first, second)`` with strictly increasing cdist."""
    workbench, _fp = henon_tangle_with_bridges
    assert _check_id_ordering(workbench) > 0, "the fixture must produce bridges"


@pytest.mark.slow
def test_bridge_id_ordering_on_period_three(henon_p3_session):
    """The same ordering invariant on the nested period-3 tangle."""
    session, _fp3, _fp1, _zone = henon_p3_session
    assert _check_id_ordering(session.workbench) > 0


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

    identified = set(map(id, workbench._bridges.values()))
    for piece in partials:
        assert piece.id is None
        assert piece in workbench.bridges, "a partial is still a registered bridge"
        assert id(piece) not in identified, "a partial must stay out of _bridges"
        assert piece in workbench._partial_bridges

    # Every id-carrying bridge is in the dict and every dict entry has an id.
    for bridge in workbench.bridges:
        assert (bridge.id is None) == bridge.partial
    assert all(bid == b.id for bid, b in workbench._bridges.items())


def test_deleted_genealogy_attributes_are_gone(henon_tangle_with_bridges):
    """The stored parent/children/next/prev links no longer exist."""
    workbench, _fp = henon_tangle_with_bridges
    bridge = workbench.bridges[0]
    for name in ("parent", "children", "next_bridge", "prev_bridge"):
        assert not hasattr(bridge, name), f"Bridge.{name} must be deleted"


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


def test_bridges_at_indexes_exactly_the_two_endpoints(henon_tangle_with_bridges):
    """Every id appears in ``bridges_at`` of exactly its two endpoints."""
    workbench, _fp = henon_tangle_with_bridges
    ids = {b.id for b in _full_bridges(workbench)}
    assert ids

    registry = workbench.intersection_registry
    for iid in registry.all_ids():
        at = workbench.bridges_at(iid)
        assert len(at) == len(set(at)), f"duplicate ids in bridges_at({iid})"
        assert set(at) == {bid for bid in ids if iid in bid}


@pytest.mark.slow
def test_bridges_at_consistency_on_period_three(henon_p3_session):
    """The reverse index stays consistent after the period-3 build."""
    session, _fp3, _fp1, _zone = henon_p3_session
    workbench = session.workbench
    ids = {b.id for b in _full_bridges(workbench)}
    for bid in ids:
        for endpoint in bid:
            assert bid in workbench.bridges_at(endpoint)


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
@pytest.mark.slow
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


@pytest.mark.slow
def test_image_and_preimage_round_trip(henon_p3_session):
    """Wherever both directions are registered, image and preimage agree."""
    session, _fp3, _fp1, _zone = henon_p3_session
    workbench = session.workbench

    checked = 0
    for bridge in _full_bridges(workbench):
        image = workbench.image_bridges(bridge.id)
        if not image:
            continue
        for child_id in image:
            preimage = workbench.preimage_bridges(child_id)
            if preimage is None:
                continue
            assert bridge.id in preimage, (
                f"{child_id} is an image of {bridge.id} but not a preimage of it"
            )
            checked += 1

    assert checked > 0, "the fixture must exercise at least one round trip"


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


@pytest.mark.slow
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
    """No Bridge owns segments in the Tangle index after ``clear_bridges``."""
    from tanglepack.numerics.Bridge import Bridge

    workbench, fp = grown_both
    workbench.compute_intersections([fp])
    workbench.trim_stable_manifolds(fp)
    bridges = workbench.create_bridges(fp)
    assert bridges

    tangle = workbench.Tangle
    assert any(
        isinstance(m, Bridge)
        for owners in tangle._seg_manifolds.values()
        for m in owners
    ), "bridges must own segments before the clear"

    workbench.clear_bridges()

    leaked = {
        m
        for owners in tangle._seg_manifolds.values()
        for m in owners
        if isinstance(m, Bridge)
    }
    assert not leaked, f"{len(leaked)} discarded bridge(s) still own segments"
    assert not any(isinstance(m, Bridge) for m in tangle._manifold_segs)


# --------------------------------------------------------------------------- #
# session delegation
# --------------------------------------------------------------------------- #
@pytest.mark.slow
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

    assert any(
        "no longer name the same crossings" in record.getMessage()
        for record in caplog.records
    ), "rebuild_bridges must report that it dropped the stale metadata carry"

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
def test_bridge_identity_on_the_inversion_saddle(henon_inversion_initialized):
    """Ids, iteration and the derived image work on an inversion fixed point.

    The inversion path advances the BRANCH index as well as the orbit index
    (:meth:`FixedPoint.advance_key`), which is exactly what ``image_bridges`` relies
    on to name the image branch, so it is worth its own pass.

    It also exposes a pre-existing defect this phase does NOT fix: the periodic point
    is registered as one crossing per (unstable branch, stable branch) pair, i.e.
    FOUR crossings at canonical distance (0, 0), so each unstable branch carries two
    distinct anchor crossings and ``create_bridges`` cuts a degenerate zero-length
    bridge between them. The ordering check below therefore exempts a pair whose two
    endpoints are both anchors -- and asserts that this is the ONLY way a tie arises,
    so the exemption becomes vacuous once anchors are deduplicated (plan 6.3) rather
    than quietly hiding a new tie.
    """
    workbench, fp = henon_inversion_initialized
    workbench.grow_n_times(fp, "unstable", num_iterations=6)
    workbench.grow_n_times(fp, "stable", num_iterations=5)
    workbench.compute_intersections([fp])
    bridges = workbench.create_bridges(fp)
    registry = workbench.intersection_registry
    tol = registry.cdist_tol

    def is_anchor(iid):
        ix = registry[iid]
        return abs(ix.unstable_cdist) <= tol and abs(ix.stable_cdist) <= tol

    full = [b for b in bridges if not b.partial]
    assert full, "the inversion fixture must produce at least one full bridge"

    proper = []
    for bridge in full:
        first, second = bridge.id
        assert bridge.id == (bridge.first_intersection, bridge.second_intersection)
        if registry[first].unstable_cdist < registry[second].unstable_cdist:
            proper.append(bridge)
            continue
        assert is_anchor(first) and is_anchor(second), (
            f"bridge {bridge.id} ties on unstable cdist away from the anchor"
        )
    assert proper, "the inversion fixture must produce a non-degenerate bridge"

    # Every proper bridge iterates onto the branch one map step forward, and the
    # DERIVED image contains exactly the bridges the cut produced (it can contain
    # more: the image arc may also cover bridges that were already there).
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
