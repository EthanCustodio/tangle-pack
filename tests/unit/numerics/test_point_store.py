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
    nodes, prev, node = [], None, manifold.root
    while node is not None:
        nodes.append(node)
        prev, node = node, manifold.walk_fwd(prev, node)
    return nodes


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


def test_link_setters_write_both_columns(store) -> None:
    """``insert_next_iterate`` writes the image row of one point and the
    preimage row of the other; clearing a link writes -1."""
    point, image = Point(0.0, 0.0, 1.0, store=store), Point(1.0, 1.0, 2.0, store=store)

    point.insert_next_iterate(image)

    assert (store._image[point.row], store._preimage[image.row]) == (image.row, point.row)
    assert point.next_iterate is image and image.prev_iterate is point
    assert point.has_iterate_link() and image.has_iterate_link()
    image.prev_iterate = None
    assert store._preimage[image.row] == -1 and not image.has_iterate_link()


def test_a_phantom_row_reads_as_a_point(store) -> None:
    """A row with no ``Point`` (a phantom) counts as a link and reads through
    ``node()`` as a ``Point`` with no cdist, created once."""
    point = Point(1.0, 1.0, 1.0, store=store)
    (phantom,) = store.append([[0.5, 0.5]], None)
    store._preimage[point.row] = phantom

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


def test_get_point_cannot_write_into_the_store(store) -> None:
    """``get_point()`` is a read-only view of the row."""
    point = Point(1.0, 2.0, store=store)

    with pytest.raises(ValueError):
        point.get_point()[0] = 5.0
    assert point.x == 1.0
    assert not hasattr(point, "__dict__")
