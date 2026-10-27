"""The point store: every manifold point a row of arrays, phantoms rows only.

``test_merge_tie_keeps_a_point_whose_only_link_is_a_phantom`` pins the merge's
tie rule (no data tie depends on the phantom: in every tie of the corpus the
kept point also holds a real link). The other tests are of the store itself.
"""

from __future__ import annotations

import numpy as np
import pytest

from tanglepack import BranchPoint
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
    periodic = BranchPoint(1)
    second.anchor = periodic

    first.absorb(second, rows)

    assert first.anchor is periodic and second.anchor is None and len(second) == 0
    other = store.new_chain(store.append(np.zeros((1, 2)), [9.0]), 1.0)
    other.anchor = BranchPoint(1)
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
