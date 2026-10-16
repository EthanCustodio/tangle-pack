"""Partition elements: identity and bridge lookup (plan Phase 5).

A stable partition is not just a list of spans — each of its intervals is an
*element*, an identifiable piece of the stable branch that a region will later be
glued to. These tests pin the lookups the region layer needs (unique ownership
of every crossing, positional ids and singletons are laws of
``tests/invariants/test_law_partition.py``):

* every non-partial bridge with an endpoint on the branch resolves, on either
  side, to the element owning that endpoint
  (:meth:`PartitionFamily.owner_of_intersection` on the branch
  :meth:`Trellis.stable_key_at` names), read off the endpoint id its BridgeId
  carries, never off the trellis's bridge list;
* ``elements_at_bridge`` agrees with the ownership table, with a ``None`` for
  an endpoint that sits on a different branch.

An element's forward image is ``DualWalk.land_element`` (``test_dual_walk.py``).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless: the session fixtures touch the plotting stack
import pytest

from tanglepack.topology.PartitionFamily import HomotopyPartition


# --------------------------------------------------------------------------- #
# fixtures
# --------------------------------------------------------------------------- #
def _wrapping_bridge(trellis):
    """A bridge of ``trellis`` whose two endpoints sit on different stable branches."""
    for bridge in trellis.bridges:
        if bridge.id is None:
            continue
        keys = [trellis.intersection(e).manifold_b_key for e in bridge.id]
        if keys[0] != keys[1]:
            return bridge
    return None


def _owner(trellis, family, bridge_id, endpoint, side):
    """The element owning one end of ``bridge_id`` on ``side`` (any side, not just its row)."""
    branch_key = trellis.stable_key_at(endpoint, trellis.unstable_key_of(bridge_id))
    return family.owner_of_intersection(endpoint, side, branch_key)


def _expected(trellis, endpoint, side):
    """The owner read straight off the result of the endpoint's own stable branch."""
    key = trellis.intersection(endpoint).manifold_b_key
    result = next(
        r for r in trellis.stable_partitions if r.branch_key == key and r.side == side
    )
    return result.ref(result.element_of_intersection[endpoint])


# --------------------------------------------------------------------------- #
# bridge -> element
# --------------------------------------------------------------------------- #
def test_every_workbench_bridge_touching_the_branch_resolves(p3_partitioned):
    """Bridges are enumerated from the WORKBENCH, not from the trellis snapshot.

    The expected set is derived independently: every non-partial bridge the
    workbench holds, keyed to a branch by its endpoints' ``manifold_b_key``. Each
    such endpoint must resolve, on both sides, to the element
    ``element_of_intersection`` names — stronger than what ``elements_at_bridge``
    alone can carry (that table only covers the bridges the snapshot happens to
    hold; see ``test_a_bridge_missing_from_the_snapshot_resolves``).
    """
    session, fp3, _fp1 = p3_partitioned
    trellis = session.trellis(fp3)
    family = HomotopyPartition.from_results(trellis.stable_partitions, trellis=trellis)
    workbench_bridge_ids = [b.id for b in session.workbench.bridges if b.id is not None]
    assert workbench_bridge_ids

    checked = 0
    for result in trellis.stable_partitions:
        expected_keys = set()
        for bridge_id in workbench_bridge_ids:
            for endpoint in bridge_id:
                if trellis.intersection(endpoint).manifold_b_key != result.branch_key:
                    continue
                expected_keys.add(bridge_id)
                assert _owner(trellis, family, bridge_id, endpoint, result.side) == (
                    result.ref(result.element_of_intersection[endpoint])
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


def test_both_ends_of_a_branch_wrapping_bridge_resolve(p3_partitioned):
    """A bridge spanning two stable branches resolves both ends in one trellis."""
    session, fp3, _fp1 = p3_partitioned
    trellis = session.trellis(fp3)
    family = HomotopyPartition.from_results(trellis.stable_partitions, trellis=trellis)

    bridge = _wrapping_bridge(trellis)
    assert bridge is not None, "the period-3 tangle has branch-wrapping bridges"

    for side in ("left", "right"):
        first, second = (
            _owner(trellis, family, bridge.id, endpoint, side) for endpoint in bridge.id
        )
        # Each end resolves against its OWN branch's result.
        assert first.branch_key != second.branch_key
        assert first == _expected(trellis, bridge.id[0], side)
        assert second == _expected(trellis, bridge.id[1], side)


def test_a_bridge_missing_from_the_snapshot_resolves(henon_p3_session):
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
    family = HomotopyPartition.from_results(trellis.stable_partitions, trellis=trellis)

    for side in ("left", "right"):
        for endpoint in dropped.id:
            assert all(dropped.id not in r.elements_at_bridge for r in trellis.stable_partitions)
            assert _owner(trellis, family, dropped.id, endpoint, side) == _expected(
                trellis, endpoint, side
            )


def test_owner_at_end_rejects_an_unknown_endpoint_name(p3_partitioned):
    """Only "first" and "second" name a bridge endpoint."""
    session, fp3, _fp1 = p3_partitioned
    trellis = session.trellis(fp3)
    family = HomotopyPartition.from_results(trellis.stable_partitions, trellis=trellis)
    bridge = next(b for b in trellis.bridges if b.id is not None)

    with pytest.raises(ValueError):
        family.owner_at_end(trellis, bridge.id, "third")
