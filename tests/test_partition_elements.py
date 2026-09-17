"""Partition elements: identity, bridge lookup, and forward images (plan Phase 5).

A stable partition is not just a list of spans — each of its intervals is an
*element*, an identifiable piece of the stable branch that a region will later be
glued to. These tests pin the three lookups the region layer needs:

* every crossing on a branch is owned by exactly one element per side
  (``StablePartitionResult.element_of_intersection``);
* every non-partial bridge with an endpoint on the branch resolves to the element
  at that endpoint (``elements_at_bridge`` / :meth:`Trellis.element_for`), with a
  ``None`` for an endpoint that sits on a different branch;
* the cross-trellis fan-out (:meth:`TangleSession.partition_element_for`) finds an
  element for a bridge whose endpoint lands on a stable branch belonging to a
  *different* trellis of the same session — the gap a single trellis cannot close;
* an element maps forward to the element(s) covering its image arc
  (:meth:`Trellis.image_of_element`), and to ``None`` only when one of its ends
  is unbounded — an endpoint that is a crossing without a registered iterate is
  mapped by canonical-distance scaling instead (see ``tests/test_image_cdist.py``).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless: the session fixtures touch the plotting stack
import pytest

from tanglepack.topology.Trellis import Trellis


# --------------------------------------------------------------------------- #
# fixtures
# --------------------------------------------------------------------------- #
@pytest.fixture
def henon_partitioned(henon_tangle_with_bridges):
    """``trellis`` for the single k=10 saddle, punched and partitioned."""
    workbench, fp = henon_tangle_with_bridges
    trellis = Trellis.from_workbench(workbench, fp)
    trellis.classify_strong_pips()
    trellis.compute_pseudoneighbors()
    trellis.punch_holes()
    trellis.partition_stable_manifold()
    assert trellis.stable_partitions
    return trellis


def _wrapping_bridge(trellis):
    """A bridge of ``trellis`` whose two endpoints sit on different stable branches."""
    for bridge in trellis.bridges:
        if bridge.id is None:
            continue
        keys = [trellis.intersection(e).manifold_b_key for e in bridge.id]
        if keys[0] != keys[1]:
            return bridge
    return None


# --------------------------------------------------------------------------- #
# element identity
# --------------------------------------------------------------------------- #
def test_element_ids_are_positional_and_carry_their_branch(p3_partitioned):
    """Each interval knows its own index, branch, and side; ids are unique."""
    session, fp3, _fp1 = p3_partitioned
    trellis = session.trellis(fp3)

    for result in trellis.stable_partitions:
        ids = [iv.element_id for iv in result.intervals]
        assert ids == list(range(len(result.intervals)))
        assert len(set(ids)) == len(ids)
        for interval in result.intervals:
            assert interval.branch_key == result.branch_key
            assert interval.side == result.side
            assert result.element(interval.element_id) is interval


def _owners_from_interval_fields(intervals, cdist, tol):
    """Element ids covering a cdist, recomputed here from the raw interval fields.

    Deliberately independent of production code (and of
    ``element_of_intersection``): a point strictly inside a span is covered, a
    point on an end only when that end is closed.
    """
    owners = []
    for interval in intervals:
        if cdist < interval.lo_cdist - tol or cdist > interval.hi_cdist + tol:
            continue
        if abs(cdist - interval.lo_cdist) <= tol:
            covered = interval.closed_lo
        elif abs(cdist - interval.hi_cdist) <= tol:
            covered = interval.closed_hi
        else:
            covered = True
        if covered:
            owners.append(interval.element_id)
    return owners


def _check_ownership(trellis):
    """Every crossing on every partitioned branch has exactly one owning element."""
    tol = trellis.registry.cdist_tol
    checked = 0
    for result in trellis.stable_partitions:
        branch = trellis.branch(result.branch_key)
        for intersection_id in branch.intersection_ids:
            cdist = float(trellis.intersection(intersection_id).stable_cdist)
            owners = _owners_from_interval_fields(result.intervals, cdist, tol)
            assert owners == [result.element_of_intersection[intersection_id]], (
                f"crossing {intersection_id} at {cdist} on {result.branch_key[1:]} "
                f"({result.side}) is covered by {owners}, reported as "
                f"{result.element_of_intersection.get(intersection_id)}"
            )
            checked += 1
        # No stray keys: the table names crossings of this branch and nothing else.
        assert set(result.element_of_intersection) == set(branch.intersection_ids)
    assert checked, "the fixture should have crossings to own"


def test_element_of_intersection_covers_every_crossing_on_the_branch(p3_partitioned):
    """Every crossing lies in exactly one element's span, and that is its owner."""
    session, fp3, _fp1 = p3_partitioned
    _check_ownership(session.trellis(fp3))


def test_simple_tangle_elements_own_every_crossing_once(henon_partitioned):
    """k=10: the same ownership invariant on a single-saddle branch."""
    _check_ownership(henon_partitioned)


def test_singleton_elements_own_exactly_their_own_point(p3_partitioned):
    """A pinched singleton [x, x] owns its crossing and nothing else.

    Note:
        The plan asked for this on the k=10 tangle, but that tangle's partition
        has no pinched point at any growth the suite can afford (its holes never
        flank one boundary from both sides); the period-3 fixture pinches the
        anchor and branch-end crossings of every branch, so it is pinned there.
    """
    session, fp3, _fp1 = p3_partitioned
    trellis = session.trellis(fp3)

    singletons = 0
    for result in trellis.stable_partitions:
        for interval in result.intervals:
            if interval.lo_id is None or interval.lo_id != interval.hi_id:
                continue
            singletons += 1
            assert interval.closed_lo and interval.closed_hi
            assert interval.lo_cdist == interval.hi_cdist
            owned = [
                iid
                for iid, element_id in result.element_of_intersection.items()
                if element_id == interval.element_id
            ]
            assert owned == [interval.lo_id]
    assert singletons, "the period-3 partition should pinch singletons"


# --------------------------------------------------------------------------- #
# bridge -> element
# --------------------------------------------------------------------------- #
def test_every_workbench_bridge_touching_the_branch_resolves(p3_partitioned):
    """Bridges are enumerated from the WORKBENCH, not from the trellis snapshot.

    The expected set is derived independently: every non-partial bridge the
    workbench holds, keyed to a branch by its endpoints' ``manifold_b_key``. Each
    such endpoint must resolve through :meth:`Trellis.element_for` to the element
    ``element_of_intersection`` names — which is the real contract, and is
    stronger than what ``elements_at_bridge`` alone can carry (that table only
    covers the bridges the snapshot happens to hold; see
    ``test_element_for_resolves_a_bridge_missing_from_the_snapshot``).
    """
    session, fp3, _fp1 = p3_partitioned
    trellis = session.trellis(fp3)
    workbench_bridge_ids = [b.id for b in session.workbench.bridges if b.id is not None]
    assert workbench_bridge_ids

    checked = 0
    for result in trellis.stable_partitions:
        expected_keys = set()
        for bridge_id in workbench_bridge_ids:
            for index, endpoint in enumerate(bridge_id):
                if trellis.intersection(endpoint).manifold_b_key != result.branch_key:
                    continue
                expected_keys.add(bridge_id)
                name = "first" if index == 0 else "second"
                assert (
                    trellis.element_for(bridge_id, name, result.side)
                    == result.element_of_intersection[endpoint]
                )
                checked += 1
        # The table may not carry every one of them (a heteroclinic bridge is
        # filtered out of this trellis), but it must never invent a key.
        assert set(result.elements_at_bridge) <= expected_keys
    assert checked, "the period-3 bridges should touch the partitioned branches"


def test_elements_at_bridge_entries_agree_with_the_ownership_table(p3_partitioned):
    """Each stored entry names the owning element per end, None off-branch."""
    session, fp3, _fp1 = p3_partitioned
    trellis = session.trellis(fp3)

    for result in trellis.stable_partitions:
        for bridge_id, ends in result.elements_at_bridge.items():
            assert len(ends) == 2
            for endpoint, element_id in zip(bridge_id, ends):
                on_branch = (
                    trellis.intersection(endpoint).manifold_b_key == result.branch_key
                )
                if on_branch:
                    assert element_id == result.element_of_intersection[endpoint]
                else:
                    assert element_id is None


def test_element_for_resolves_both_ends_of_a_branch_wrapping_bridge(p3_partitioned):
    """A bridge spanning two stable branches resolves both ends in one trellis."""
    session, fp3, _fp1 = p3_partitioned
    trellis = session.trellis(fp3)

    bridge = _wrapping_bridge(trellis)
    assert bridge is not None, "the period-3 tangle has branch-wrapping bridges"

    for side in ("left", "right"):
        first = trellis.element_for(bridge.id, "first", side)
        second = trellis.element_for(bridge.id, "second", side)
        assert first is not None and second is not None
        # Each end resolves against its OWN branch's result.
        for endpoint, element_id in ((bridge.id[0], first), (bridge.id[1], second)):
            key = trellis.intersection(endpoint).manifold_b_key
            result = next(
                r
                for r in trellis.stable_partitions
                if r.branch_key == key and r.side == side
            )
            assert result.element_of_intersection[endpoint] == element_id


def test_element_for_resolves_a_bridge_missing_from_the_snapshot(henon_p3_session):
    """A bridge absent from ``trellis.bridges`` still resolves, via its BridgeId.

    This is the heteroclinic case in miniature, and the reason the lookup must
    not depend on the snapshot's bridge list. ``Trellis.from_workbench`` keeps a
    bridge under its UNSTABLE fixed point, so a bridge of W^u(fp1) cut at
    crossings on W^s(fp3) is filtered out of fp3's trellis — the very trellis
    whose partition owns both its endpoints. A BridgeId IS that endpoint pair, so
    the element is reachable regardless; here one bridge is dropped from the
    snapshot to reproduce exactly that state on the available fixture.
    """
    session, fp3, _fp1, _zone = henon_p3_session
    trellis = session.trellis(fp3)
    trellis.classify_strong_pips()
    trellis.compute_pseudoneighbors()
    trellis.punch_holes()

    stable_keys = {b.key for b in trellis.stable_branches}
    dropped = next(
        b
        for b in trellis.bridges
        if b.id is not None
        and all(trellis.intersection(e).manifold_b_key in stable_keys for e in b.id)
    )
    trellis.bridges = [b for b in trellis.bridges if b is not dropped]
    trellis.partition_stable_manifold()

    for side in ("left", "right"):
        for index, endpoint in enumerate(dropped.id):
            key = trellis.intersection(endpoint).manifold_b_key
            result = next(
                r
                for r in trellis.stable_partitions
                if r.branch_key == key and r.side == side
            )
            expected = result.element_of_intersection[endpoint]
            assert dropped.id not in result.elements_at_bridge

            name = "first" if index == 0 else "second"
            assert trellis.element_for(dropped.id, name, side) == expected
            assert session.partition_element_for(dropped.id, name, side) == expected


def test_session_closes_the_cross_trellis_gap(p3_partitioned):
    """A bridge of one tangle finds its element through the session, not one trellis.

    Covers the session's fan-out over the trellis CACHE, not the heteroclinic
    endpoint case (the nested fixture has no heteroclinic crossings — every
    crossing's two manifold keys share a fixed point; see
    ``test_element_for_resolves_a_bridge_missing_from_the_snapshot`` for that
    mechanism). The period-1 and period-3 tangles live in separate trellises of
    one session, so the period-3 trellis cannot answer for a period-1 bridge —
    neither its bridge list nor its partitions mention it — while the session
    scans every cached trellis and finds the answer in the period-1 one.
    """
    session, fp3, fp1 = p3_partitioned
    t3 = session.trellis(fp3)
    t1 = session.trellis(fp1)

    foreign = next(b for b in t1.bridges if b.id is not None)
    assert t3.element_for(foreign.id, "first", "left") is None
    assert t3.element_for(foreign.id, "second", "left") is None

    resolved = [
        session.partition_element_for(foreign.id, endpoint, "left")
        for endpoint in ("first", "second")
    ]
    assert any(element is not None for element in resolved)
    for endpoint, element_id in zip(("first", "second"), resolved):
        assert element_id == t1.element_for(foreign.id, endpoint, "left")


def test_element_for_rejects_an_unknown_endpoint_name(p3_partitioned):
    """Only "first" and "second" name a bridge endpoint."""
    session, fp3, _fp1 = p3_partitioned
    trellis = session.trellis(fp3)
    bridge = next(b for b in trellis.bridges if b.id is not None)

    with pytest.raises(ValueError):
        trellis.element_for(bridge.id, "third", "left")


# --------------------------------------------------------------------------- #
# element -> image element
# --------------------------------------------------------------------------- #
def test_image_of_element_lands_on_the_advanced_branch(p3_partitioned):
    """One forward step sends an element onto elements of ``advance_key(key, 1)``.

    The image elements tile the arc between the two endpoint images, so their
    combined cdist span contains both images' cdists.
    """
    session, fp3, _fp1 = p3_partitioned
    trellis = session.trellis(fp3)
    tol = trellis.registry.cdist_tol

    checked = 0
    for result in trellis.stable_partitions:
        image_key = fp3.advance_key(result.branch_key, 1)
        for interval in result.intervals:
            if interval.lo_id is None or interval.hi_id is None:
                continue
            images = [
                trellis.iterate(interval.lo_id, 1),
                trellis.iterate(interval.hi_id, 1),
            ]
            if any(image is None for image in images):
                continue

            element_ids = trellis.image_of_element(result, interval.element_id, 1)
            assert element_ids, "an image arc with both ends known is covered"

            image_result = next(
                r
                for r in trellis.stable_partitions
                if r.branch_key == image_key and r.side == result.side
            )
            covered = [image_result.element(e) for e in element_ids]
            assert all(iv.branch_key == image_key for iv in covered)

            span_lo = min(iv.lo_cdist for iv in covered)
            span_hi = max(iv.hi_cdist for iv in covered)
            for image in images:
                cdist = float(trellis.intersection(image).stable_cdist)
                assert span_lo - tol <= cdist <= span_hi + tol
            checked += 1
    assert checked, "the period-3 partition should have iterable elements"


def test_image_of_element_is_none_exactly_at_an_unbounded_end(p3_partitioned):
    """Only an unbounded end has no image arc; that is the sole None.

    An element bounded by the anchor (``lo_id is None``) or by the computed end of
    the branch (``hi_id is None``) has no crossing whose image would bound its
    image, so there is nothing to cover and the answer is None.
    """
    session, fp3, _fp1 = p3_partitioned
    trellis = session.trellis(fp3)

    unbounded = 0
    for result in trellis.stable_partitions:
        for interval in result.intervals:
            if interval.lo_id is None or interval.hi_id is None:
                assert trellis.image_of_element(result, interval.element_id, 1) is None
                unbounded += 1
            else:
                assert (
                    trellis.image_of_element(result, interval.element_id, 1)
                    is not None
                )
    assert unbounded, "the anchor/outermost elements are unbounded at one end"


def test_image_of_element_falls_back_when_an_iterate_is_missing(p3_partitioned):
    """A bounded element with no registered iterate is still mapped, not dropped.

    This is the Phase C.2 fallback: ``image_cdist`` scales the endpoint's stable
    canonical distance onto ``advance_key`` when the iterate table has no entry,
    so every element with two real endpoints has an image. Both tangles of the
    fixture are scanned: since holes propagate backward only (2026-09-16) the
    period-3 partition is coarse enough that every one of its bounded elements
    has registered endpoint iterates, and the unregistered ones live on the
    period-1 saddle's outermost elements.
    """
    session, fp3, fp1 = p3_partitioned

    fell_back = 0
    for fixed_point in (fp3, fp1):
        trellis = session.trellis(fixed_point)
        for result in trellis.stable_partitions:
            for interval in result.intervals:
                if interval.lo_id is None or interval.hi_id is None:
                    continue
                if (
                    trellis.iterate(interval.lo_id, 1) is not None
                    and trellis.iterate(interval.hi_id, 1) is not None
                ):
                    continue
                images = trellis.image_of_element(result, interval.element_id, 1)
                assert images is not None, "a bounded element always has an image arc"
                fell_back += 1
    assert fell_back, "the outermost elements have no registered forward iterate"


def test_image_of_element_rejects_an_unknown_element_id(p3_partitioned):
    """An element id outside the result is a caller error, not an empty answer."""
    session, fp3, _fp1 = p3_partitioned
    trellis = session.trellis(fp3)
    result = trellis.stable_partitions[0]

    with pytest.raises(IndexError):
        trellis.image_of_element(result, len(result.intervals), 1)
    with pytest.raises(IndexError):
        trellis.image_of_element(result, -1, 1)


def test_image_of_element_requires_the_image_branchs_partition(henon_p3_session):
    """Asking for an image before the image branch is partitioned is an error."""
    session, fp3, _fp1, _zone = henon_p3_session
    trellis = session.trellis(fp3)
    trellis.classify_strong_pips()
    trellis.compute_pseudoneighbors()
    trellis.punch_holes()

    branch_key = next(
        b.key for b in trellis.stable_branches if fp3.advance_key(b.key, 1) != b.key
    )
    result = trellis.partition_stable_manifold(branch_key)[0]
    element = next(
        iv
        for iv in result.intervals
        if iv.lo_id is not None
        and iv.hi_id is not None
        and trellis.iterate(iv.lo_id, 1) is not None
        and trellis.iterate(iv.hi_id, 1) is not None
    )

    # n = 0 never leaves the branch, so it answers without an image partition.
    assert trellis.image_of_element(result, element.element_id, 0) == [
        element.element_id
    ]
    with pytest.raises(ValueError):
        trellis.image_of_element(result, element.element_id, 1)
