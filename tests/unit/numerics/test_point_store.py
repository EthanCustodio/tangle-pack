"""The point store: every manifold point a row of arrays, phantoms rows only.

``test_merge_tie_keeps_a_point_whose_only_link_is_a_phantom`` pins the merge's
tie rule and uses only ``Point(x, y)`` and ``insert_prev_iterate``, so it runs
on the code before the store as well (no data tie depends on the phantom: in
every tie of the corpus the kept point also holds a real link). The other tests
are of the store itself.
"""

from __future__ import annotations

import numpy as np
import pytest

from tanglepack import Point
from tanglepack.examples.henon import henon_map, henon_map_inverse
from tanglepack.numerics import DynamicalSystem, ManifoldMachine


def _curve(points: list[Point]):
    """An unstable key-less curve through ``points``, linked forward."""
    return ManifoldMachine._manifold_from_points(points, "unstable", 1.0, None)


def _walk(manifold) -> list[Point]:
    """The curve's nodes from its root to the end of its chain."""
    return list(manifold._walk_nodes(stop_at_final=False))


def test_merge_tie_keeps_a_point_whose_only_link_is_a_phantom() -> None:
    """Of two points on one cdist the one with an iterate link is kept, even
    when that link is a coordinate-only phantom (cdist None) and the other
    point's link is a real one."""
    head_1 = Point(1.0, 1.0, 1.0)
    head_2 = Point(1.0, 1.1, 1.0)
    head_1.insert_prev_iterate(Point(0.5, 0.5))
    head_2.insert_prev_iterate(Point(0.5, 0.6, 0.4))
    assert head_1.prev_iterate.cdist is None
    manifold_1 = _curve([Point(0.0, 0.0, 0.0), head_1, Point(2.0, 2.0, 2.0)])
    manifold_2 = _curve([Point(0.5, 0.5, 0.5), head_2, Point(3.0, 3.0, 3.0)])
    machine = ManifoldMachine(DynamicalSystem(henon_map(10), henon_map_inverse(10)))

    merged = _walk(machine.merge_manifolds(manifold_1, manifold_2))

    assert any(node is head_1 for node in merged)
    assert not any(node is head_2 for node in merged)
    assert [node.cdist for node in merged] == [0.0, 0.5, 1.0, 2.0, 3.0]


# ------------------------------------------------------------- the store

@pytest.fixture
def store():
    """A fresh store (the module exists only from the store's commit on)."""
    return pytest.importorskip("tanglepack.numerics.PointStore").PointStore(capacity=2)


def test_rows_are_append_only_across_a_resize(store) -> None:
    """Rows keep their ids and values when the columns grow; none is reused."""
    first = Point(1.0, 2.0, 3.0, store=store)
    rows = store.append(np.arange(20.0).reshape(10, 2), np.arange(10.0))
    last = Point(7.0, 8.0, store=store)

    assert (first.row, last.row) == (0, 11)
    np.testing.assert_array_equal(rows, np.arange(1, 11))
    assert len(store) == len(store.nodes) == 12
    assert (first.x, first.y, first.cdist) == (1.0, 2.0, 3.0)
    assert last.cdist is None
    np.testing.assert_array_equal(store.xy[rows], np.arange(20.0).reshape(10, 2))


def test_a_view_taken_before_a_resize_still_reads_its_row(store) -> None:
    """A row's coordinates never change, so a ``get_point()`` view of the
    columns before they grew reads the same row after (Q9's read-only view)."""
    point = Point(1.0, 2.0, store=store)
    view = point.get_point()

    store.append(np.ones((50, 2)), np.ones(50))

    assert not np.shares_memory(view, store.xy)
    np.testing.assert_array_equal(view, [1.0, 2.0])
    np.testing.assert_array_equal(point.get_point(), [1.0, 2.0])


def test_link_setters_write_both_columns(store) -> None:
    """``insert_next_iterate`` writes the image row of one point and the
    preimage row of the other; clearing a link writes -1."""
    point, image = Point(0.0, 0.0, 1.0, store=store), Point(1.0, 1.0, 2.0, store=store)

    point.insert_next_iterate(image)

    assert (store.image[point.row], store.preimage[image.row]) == (image.row, point.row)
    assert point.next_iterate is image and image.prev_iterate is point
    assert point.has_iterate_link() and image.has_iterate_link()
    image.prev_iterate = None
    assert store.preimage[image.row] == -1 and not image.has_iterate_link()


def test_a_phantom_row_reads_as_a_point(store) -> None:
    """A row with no ``Point`` (a phantom) counts as a link and reads through
    ``node()`` as a ``Point`` with no cdist, created once."""
    point = Point(1.0, 1.0, 1.0, store=store)
    (phantom,) = store.append([[0.5, 0.5]], None)
    store.link_phantoms([point.row], [phantom], "unstable")

    assert store.nodes[phantom] is None
    assert point.has_iterate_link()
    preimage = point.prev_iterate
    assert isinstance(preimage, Point) and preimage.cdist is None
    assert preimage.row == phantom and store.node(phantom) is preimage
    np.testing.assert_array_equal(preimage.get_point(), [0.5, 0.5])


def test_a_link_across_stores_raises(store) -> None:
    """An iterate link names a ``Point`` of the same store."""
    other = type(store)()
    point = Point(0.0, 0.0, store=store)

    with pytest.raises(TypeError):
        point.next_iterate = Point(1.0, 1.0, store=other)
    with pytest.raises(TypeError):
        point.prev_iterate = object()


def test_the_columns_are_fresh_read_only_views_of_the_written_rows(store) -> None:
    """``xy``, ``cdist``, ``image`` and ``preimage`` cover the rows written so
    far, refuse writes, and a view taken before a resize keeps its rows."""
    Point(1.0, 2.0, 3.0, store=store)
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


def test_get_point_cannot_write_into_the_store(store) -> None:
    """``get_point()`` is a read-only view of the row."""
    point = Point(1.0, 2.0, store=store)

    with pytest.raises(ValueError):
        point.get_point()[0] = 5.0
    assert point.x == 1.0
    assert not hasattr(point, "__dict__")


def test_a_blasted_build_leaves_phantoms_as_rows(k28_partitioned) -> None:
    """Growth, a blast and the topology never give a phantom a ``Point``: a row
    with no cdist is either a phantom (no ``Point``, no link of its own, named
    by a point's link) or an initializer's detached fictitious tail (a
    ``Point`` that keeps its one link to the chain it was cut from)."""
    pytest.importorskip("tanglepack.numerics.PointStore")
    session, _ = k28_partitioned
    store = session.workbench.point_store
    n = len(store)
    no_cdist = np.isnan(store.cdist)
    has_point = np.array([node is not None for node in store.nodes])
    own_link = (store.image >= 0) | (store.preimage >= 0)
    named = np.zeros(n, bool)
    for column in (store.image, store.preimage):
        named[column[column >= 0]] = True

    phantoms = no_cdist & ~has_point
    assert phantoms.any()
    assert not (phantoms & own_link).any() and named[phantoms].all()
    assert own_link[no_cdist & has_point].all()
