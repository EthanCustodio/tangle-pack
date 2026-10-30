"""The point store: every manifold point a row of arrays, phantoms rows only.

``test_merge_tie_keeps_a_point_whose_only_link_is_a_phantom`` pins the merge's
tie rule (no data tie depends on the phantom: in every tie of the corpus the
kept point also holds a real link). The other tests are of the store itself.
"""

from __future__ import annotations

import numpy as np
import pytest

from helpers.polylines import periodic_point
from tanglepack.examples.henon import henon_map, henon_map_inverse
from tanglepack.numerics import DynamicalSystem, ManifoldMachine
from tanglepack.numerics.PointStore import PointStore


def _curve(store: PointStore, rows: list[int]):
    """An unstable key-less curve through ``rows``."""
    return ManifoldMachine._manifold_from_rows(store, rows, "unstable", 1.0, None)


def _walk(manifold) -> list[int]:
    """The curve's rows from its root to the end of its chain."""
    store, root = manifold.store, int(manifold.ids[0])
    return store.chain_of(root).ids[store.position_of(root):].tolist()


def test_merge_tie_keeps_a_point_whose_only_link_is_a_phantom() -> None:
    """Of two points on one cdist the one with an iterate link is kept, even
    when that link is a coordinate-only phantom (cdist None) and the other
    point's link is a real one."""
    store = PointStore()
    head_1 = store.append_one(1.0, 1.0, 1.0)
    head_2 = store.append_one(1.0, 1.1, 1.0)
    store.link([store.append_one(0.5, 0.5, None)], [head_1])
    store.link([store.append_one(0.5, 0.6, 0.4)], [head_2])
    assert np.isnan(store.cdist[store.preimage[head_1]])
    manifold_1 = _curve(store, [store.append_one(0.0, 0.0, 0.0), head_1, store.append_one(2.0, 2.0, 2.0)])
    manifold_2 = _curve(store, [store.append_one(0.5, 0.5, 0.5), head_2, store.append_one(3.0, 3.0, 3.0)])
    machine = ManifoldMachine(DynamicalSystem(henon_map(10), henon_map_inverse(10)))

    merged = _walk(machine.merge_manifolds(manifold_1, manifold_2))

    assert head_1 in merged
    assert head_2 not in merged
    assert store.cdist[merged].tolist() == [0.0, 0.5, 1.0, 2.0, 3.0]


# ------------------------------------------------------------- the store

@pytest.fixture
def store():
    """A fresh store (the module exists only from the store's commit on)."""
    return pytest.importorskip("tanglepack.numerics.PointStore").PointStore(capacity=2)


def test_rows_are_append_only_across_a_resize(store) -> None:
    """Rows keep their ids and values when the columns grow; none is reused."""
    first = store.append_one(1.0, 2.0, 3.0)
    rows = store.append(np.arange(20.0).reshape(10, 2), np.arange(10.0))
    last = store.append_one(7.0, 8.0, None)

    assert (first, last) == (0, 11)
    np.testing.assert_array_equal(rows, np.arange(1, 11))
    assert len(store) == 12
    assert (*store.xy[first], store.cdist[first]) == (1.0, 2.0, 3.0)
    assert np.isnan(store.cdist[last])
    np.testing.assert_array_equal(store.xy[last], [7.0, 8.0])
    np.testing.assert_array_equal(store.xy[rows], np.arange(20.0).reshape(10, 2))


def test_a_view_taken_before_a_resize_still_reads_its_row(store) -> None:
    """A row's coordinates never change, so a view of the columns before they
    grew reads the same row after (Q9's read-only view)."""
    row = store.append_one(1.0, 2.0, None)
    view = store.xy[row]

    store.append(np.ones((50, 2)), np.ones(50))

    assert not np.shares_memory(view, store.xy)
    np.testing.assert_array_equal(view, [1.0, 2.0])
    np.testing.assert_array_equal(store.xy[row], [1.0, 2.0])


def test_an_end_is_a_value_of_its_row(store) -> None:
    """``end(row)`` is a record of the row's columns built per read: two reads
    are equal and are not one object, so ends are compared by value (by ``id``
    when the cdist may be NaN), never by identity."""
    row = store.append_one(1.0, 2.0, 3.0)

    first, second = store.end(row), store.end(row)

    assert first == second and first is not second
    assert (first.id, first.x, first.y, first.cdist) == (row, *store.xy[row], store.cdist[row])


def test_the_columns_are_fresh_read_only_views_of_the_written_rows(store) -> None:
    """``xy``, ``cdist``, ``image`` and ``preimage`` cover the rows written so
    far, refuse writes, and a view taken before a resize keeps its rows."""
    store.append_one(1.0, 2.0, 3.0)
    before = store.xy

    store.append(np.ones((5, 2)), np.ones(5))

    for column in (store.xy, store.cdist, store.image, store.preimage):
        assert len(column) == len(store) == 6
        with pytest.raises(ValueError):
            column[0] = 0
    assert store.xy is not store.xy and not np.shares_memory(before, store.xy)
    np.testing.assert_array_equal(before[0], [1.0, 2.0])
    assert np.isnan(store.cdist).sum() == 0 and (store.image == -1).all()


def test_the_writers_set_links_and_cdists(store) -> None:
    """``link`` writes both sides of each pair, ``link_phantoms`` one side (on the
    preiterate side of the stability), ``set_cdist`` one row (None as NaN)."""
    rows = store.append(np.zeros((4, 2)), np.arange(4.0))

    store.link(rows[:2], rows[2:])
    np.testing.assert_array_equal(store.image[rows], [2, 3, -1, -1])
    np.testing.assert_array_equal(store.preimage[rows], [-1, -1, 0, 1])

    (phantom,) = store.append([[9.0, 9.0]], None)
    store.link_phantoms([rows[0]], [phantom], "stable")
    assert store.image[rows[0]] == phantom and store.preimage[phantom] == -1

    store.set_cdist(rows[1], None)
    store.set_cdist(rows[2], 7.5)
    assert np.isnan(store.cdist[rows[1]]) and store.cdist[rows[2]] == 7.5


def test_a_chain_keeps_its_columns_through_insert_and_absorb(store) -> None:
    """Every chain write keeps the chain and position columns exact; an absorbed
    chain is retired and a row the merge left off is marked as such."""
    from tanglepack.numerics.PointStore import LEFT_OFF

    rows = store.append(np.zeros((6, 2)), np.arange(6.0))
    first = store.new_chain(rows[[0, 2]], 1.0)
    second = store.new_chain(rows[[3, 4]], 1.0)
    first.insert([1], rows[[1]])
    np.testing.assert_array_equal(first.ids, rows[:3])

    first.absorb(second, np.concatenate((first.ids, rows[[3]])), start=3)

    np.testing.assert_array_equal(first.ids, rows[:4])
    assert len(second) == 0 and second.anchor is None
    for position, row in enumerate(rows[:4]):
        assert store.chain_of(row) is first and store.position_of(row) == position
    assert store.chain_of(rows[4]) is None and store.position_of(rows[4]) == LEFT_OFF
    assert store.chain_of(rows[5]) is None and store.position_of(rows[5]) == -1


def test_an_absorbed_chain_hands_on_its_anchor_and_one_chain_holds_one(store) -> None:
    """The periodic point a chain hangs off survives the chain being absorbed;
    two chains each hanging off one cannot merge."""
    rows = store.append(np.zeros((4, 2)), np.arange(4.0))
    first, second = store.new_chain(rows[:2], 1.0), store.new_chain(rows[2:], 1.0)
    periodic = periodic_point()
    second.anchor = periodic

    first.absorb(second, rows)

    assert first.anchor is periodic and second.anchor is None and len(second) == 0
    other = store.new_chain(store.append(np.zeros((1, 2)), [9.0]), 1.0)
    other.anchor = periodic_point(1.0, 0.0)
    with pytest.raises(AssertionError, match="two anchors"):
        first.absorb(other, np.concatenate((first.ids, other.ids)))


def test_a_blasted_build_leaves_phantoms_as_rows(k28_partitioned) -> None:
    """Growth, a blast and the topology leave every row with no cdist either a
    phantom (no link of its own, named by a point's link) or an initializer's
    detached fictitious tail (it keeps its one link to the chain it was cut
    from, which names nothing back)."""
    pytest.importorskip("tanglepack.numerics.PointStore")
    session, _ = k28_partitioned
    store = session.workbench.point_store
    n = len(store)
    no_cdist = np.isnan(store.cdist)
    own_link = (store.image >= 0) | (store.preimage >= 0)
    named = np.zeros(n, bool)
    for column in (store.image, store.preimage):
        named[column[column >= 0]] = True

    phantoms, tails = no_cdist & ~own_link, no_cdist & own_link
    assert phantoms.any() and named[phantoms].all()
    assert tails.any() and not named[tails].any()


# ------------------------------------------------------------- compaction


def test_compact_keeps_every_surviving_row(store) -> None:
    """A compaction renumbers the kept rows in order: each keeps its point, its
    cdist and every link into a kept row (one into a dropped row reads -1), and
    every chain keeps its kept rows in order, renumbered from 0."""
    rows = store.append(np.arange(16.0).reshape(8, 2), np.arange(8.0))
    phantom = store.append_one(0.5, 0.5, None)
    chain = store.new_chain(rows[[0, 1, 2, 3, 4]], 1.0)
    other = store.new_chain(rows[[5, 6, 7]], 1.0)
    store.link(rows[[0, 1]], rows[[5, 6]])
    store.link_phantoms(rows[[3]], [phantom], "unstable")
    keep = np.ones(len(store), dtype=bool)
    keep[[rows[1], rows[2], phantom]] = False
    before_xy, before_cdist = store.xy.copy(), store.cdist.copy()

    old_to_new = store.compact(keep)

    assert len(store) == int(keep.sum())
    assert old_to_new.tolist() == [0, -1, -1, 1, 2, 3, 4, 5, -1]
    survivors = np.flatnonzero(keep)
    np.testing.assert_array_equal(store.xy, before_xy[survivors])
    np.testing.assert_array_equal(store.cdist, before_cdist[survivors])
    assert store.image[old_to_new[rows[0]]] == old_to_new[rows[5]]
    assert store.preimage[old_to_new[rows[5]]] == old_to_new[rows[0]]
    assert store.preimage[old_to_new[rows[6]]] == -1  # its preimage was dropped
    assert store.preimage[old_to_new[rows[3]]] == -1  # its phantom was dropped
    assert chain.ids.tolist() == old_to_new[rows[[0, 3, 4]]].tolist()
    assert other.ids.tolist() == old_to_new[rows[[5, 6, 7]]].tolist()
    for each in (chain, other):
        for position, row in enumerate(each.ids):
            assert store.chain_of(row) is each and store.position_of(row) == position


def _lean_k28():
    """The k28 zone blasted lean (width guard) until a bridge is released."""
    from cases import build_k28

    case = build_k28(blasts=0, through="zones")
    session, zone = case.session, case.zones[0]
    for _ in range(8):
        session.blast_zone(zone, 1, fixed_point=[zone.fixed_point], min_separation="width", lean=True)
        if any(bridge.released for bridge in session.workbench.bridges):
            return session
    raise AssertionError("no bridge was released")


def test_a_release_keeps_the_invariants() -> None:
    """After a release every curve still read is cdist non-decreasing and ends
    on rows of its chain, the iterate table is a partial bijection and every
    crossing is unstable x stable, its brackets rows of the store."""
    workbench = _lean_k28().workbench
    store = workbench.point_store
    curves = [*workbench.manifolds.values(), *(b for b in workbench.bridges if not b.released)]
    for curve in curves:
        cdists = curve.get_cdist_array().ravel()
        assert (np.diff(cdists) >= 0).all()
        for end in (curve.root, curve.tail):
            if hasattr(end, "id"):
                assert store.chain_of(end.id) is not None
    registry = workbench.intersection_registry
    images = [registry.iterate_table[iid, 1] for iid, _ in registry]
    linked = [image for image in images if image is not None]
    assert len(linked) == len(set(linked))
    for _iid, crossing in registry:
        assert crossing.manifold_a_key[1] == "unstable" and crossing.manifold_b_key[1] == "stable"
        for end in (*(crossing.unstable_segment or ()), *(crossing.stable_segment or ())):
            assert not isinstance(end, int) or 0 <= end < len(store)
        brackets = crossing.brackets if crossing.unstable_segment is not None else ()
        for stability, (lo, hi) in zip(("unstable", "stable"), brackets):
            assert sorted(crossing.segment_ends(stability)[1].tolist()) == [lo, hi]


def test_a_release_renumbers_stable_brackets_above_freed_rows() -> None:
    """Stable growth after a full blast puts stable bracket rows above blast
    rows a release frees: each bracket still reads its own two points."""
    from cases import build_k28

    case = build_k28(blasts=0, through="zones")
    session, zone = case.session, case.zones[0]
    fp = zone.fixed_point
    session.blast_zone(zone, 2, fixed_point=[fp], diagnose=False)
    session.grow_n_times(fp, "stable", 1)
    session.compute_intersections([fp], preserve_ids=True)
    session.create_bridges(fp)
    workbench = session.workbench
    registry = workbench.intersection_registry
    grown = len(workbench.point_store)

    def stable_ends():
        return {
            iid: tuple(map(np.ndarray.tolist, crossing.segment_ends("stable")))
            for iid, crossing in registry if crossing.stable_segment is not None
        }

    before = stable_ends()
    assert workbench.release_curves(set()) > 0
    assert len(workbench.point_store) < grown
    assert stable_ends() == before


def test_a_straddle_neighbours_images_survive_a_release() -> None:
    """A kept un-iterated bridge whose straddle rows an iterated neighbour already
    mapped keeps those images: released on its own, it iterates as on a full twin."""
    from cases import build_k28

    def frontier_with_images():
        case = build_k28(blasts=0, through="zones")
        session, zone = case.session, case.zones[0]
        session.blast_zone(zone, 1, fixed_point=[zone.fixed_point], diagnose=False)
        workbench = session.workbench
        for bridge in workbench.uniiterated_bridges:
            rows = bridge.ids
            if bridge.id is not None and (workbench.point_store.image[rows] >= 0).any():
                return workbench, bridge
        raise AssertionError("no frontier bridge shares an iterated neighbour's rows")

    full, full_bridge = frontier_with_images()
    lean, lean_bridge = frontier_with_images()
    assert full_bridge.id == lean_bridge.id
    assert lean.release_curves({lean_bridge.id}) > 0
    assert lean.lean and any(b.released for b in lean.bridges)

    full_children, lean_children = full.iterate_bridge(full_bridge), lean.iterate_bridge(lean_bridge)

    assert [c.id for c in lean_children] == [c.id for c in full_children]
    for a, b in zip(full_children, lean_children):
        np.testing.assert_array_equal(a.get_point_array(), b.get_point_array())


def test_a_release_leaves_a_bracket_end_at_the_periodic_point() -> None:
    """A bracket on the first segment from the anchor ends on the ``BranchPoint``,
    which is no row: a release keeps it and renumbers the other end."""
    from cases import build_k28
    from tanglepack.numerics.BranchPoint import BranchPoint

    case = build_k28(blasts=0, through="zones")
    session, zone = case.session, case.zones[0]
    session.blast_zone(zone, 2, fixed_point=[zone.fixed_point], diagnose=False)
    workbench = session.workbench
    store = workbench.point_store
    crossing = next(x for _iid, x in workbench.intersection_registry if x.stable_segment is not None)
    anchor = zone.fixed_point.branch_points[crossing.manifold_b_key[2]]
    row = crossing.stable_segment[1]
    point = store.xy[row].copy()
    crossing.stable_segment = (anchor, row)
    assert workbench.release_curves(set()) > 0
    lo, hi = crossing.stable_segment
    assert isinstance(lo, BranchPoint) and lo is anchor
    assert np.array_equal(store.xy[hi], point)


# ------------------------------------------------------------- the residual column


def test_the_residual_travels_with_its_row(store) -> None:
    """``cdist_lo`` is written by every writer (0.0 by default, NaN with no cdist),
    read by ``end`` and ``precise_cdist`` and kept by a resize and a compaction."""
    rows = store.append(np.zeros((3, 2)), [1.0, 1.0, 2.0], [1e-17, -1e-17, 0.0])
    plain = store.append_one(0.0, 0.0, 3.0)
    (phantom,) = store.append([[9.0, 9.0]], None)
    store.set_cdist(rows[2], 2.0, 2e-17)
    assert store.end(rows[0]).cdist_lo == 1e-17 and store.end(plain).cdist_lo == 0.0
    assert np.isnan(store._cdist_lo[phantom])

    store.append(np.zeros((5000, 2)), np.ones(5000))  # a resize
    hi, lo = store.precise_cdist(rows)
    assert hi.tolist() == [1.0, 1.0, 2.0] and lo.tolist() == [1e-17, -1e-17, 2e-17]

    keep = np.ones(len(store), dtype=bool)
    keep[rows[1]] = False
    old_to_new = store.compact(keep)
    assert store.precise_cdist(old_to_new[rows[[0, 2]]])[1].tolist() == [1e-17, 2e-17]
