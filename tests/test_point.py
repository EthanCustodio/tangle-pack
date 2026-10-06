"""Linked-list behaviour of ``Point`` and ``BranchPoint``.

Both point kinds carry the iterate list (``next_iterate`` / ``prev_iterate``)
inherited from ``BasePoint``; the iterate tests are parametrized over the two
kinds. The geometric list differs: a ``Point`` has one ``forward`` /
``backward`` neighbour, a ``BranchPoint`` one per branch.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
import pytest

from tanglepack import BranchPoint, Point
from tanglepack.numerics.BasePoint import BasePoint

NUM_BRANCHES = 2

#: Factories for the two point kinds the iterate list is shared by.
POINT_KINDS: dict[str, Callable[[], BasePoint]] = {
    "point": lambda: Point(),
    "branch_point": lambda: BranchPoint(NUM_BRANCHES),
}


# ---------------------------------------------------------------- creation

def test_point_creation() -> None:
    """A fresh ``Point`` holds its coordinates and no geometric neighbours."""
    p = Point(1.0, 2.0)

    assert p.x == 1.0
    assert p.y == 2.0
    assert p.forward is None
    assert p.backward is None
    np.testing.assert_array_equal(p.get_point(), np.array([1.0, 2.0]))


def test_branch_point_creation() -> None:
    """A fresh ``BranchPoint`` has one empty slot per branch in each direction."""
    p = BranchPoint(NUM_BRANCHES, x=1.0, y=2.0)

    assert p.x == 1.0
    assert p.y == 2.0
    assert p.next_iterate is None
    assert p.prev_iterate is None
    assert len(p.forward_branches) == NUM_BRANCHES
    assert len(p.backward_branches) == NUM_BRANCHES
    np.testing.assert_array_equal(p.get_point(), np.array([1.0, 2.0]))


# ------------------------------------------------------- geometric (Point)

def test_insert_point_forward() -> None:
    """Forward insertion links both directions and leaves the ends open."""
    p1, p2 = Point(), Point()

    p1.insert_point_forward(p2)

    assert p1.forward is p2
    assert p2.backward is p1
    assert p1.backward is None
    assert p2.forward is None


def test_insert_point_backward() -> None:
    """Backward insertion links both directions and leaves the ends open."""
    p1, p2 = Point(), Point()

    p1.insert_point_backward(p2)

    assert p1.backward is p2
    assert p2.forward is p1
    assert p1.forward is None
    assert p2.backward is None


def test_insert_point_forward_connected() -> None:
    """Forward insertion splices into an existing list."""
    p1, p2, p3 = Point(), Point(), Point()

    p1.insert_point_forward(p3)
    p1.insert_point_forward(p2)

    assert (p1.forward, p2.forward, p3.forward) == (p2, p3, None)
    assert (p3.backward, p2.backward, p1.backward) == (p2, p1, None)


def test_insert_point_backward_connected() -> None:
    """Backward insertion splices into an existing list."""
    p1, p2, p3 = Point(), Point(), Point()

    p3.insert_point_backward(p1)
    p3.insert_point_backward(p2)

    assert (p1.forward, p2.forward, p3.forward) == (p2, p3, None)
    assert (p3.backward, p2.backward, p1.backward) == (p2, p1, None)


# ------------------------------------------------- geometric (BranchPoint)

def test_branch_point_insert_point_forward() -> None:
    """Each forward branch slot receives its own point."""
    p = BranchPoint(NUM_BRANCHES)
    p1, p2 = Point(), Point()

    p.insert_point_forward(p1, 0)
    p.insert_point_forward(p2, 1)

    assert p.forward_branches[0] is p1
    assert p.forward_branches[1] is p2
    assert p1.backward is p
    assert p2.backward is p


def test_branch_point_insert_point_backward() -> None:
    """Each backward branch slot receives its own point; forward slots stay empty."""
    p = BranchPoint(NUM_BRANCHES)
    p1, p2 = Point(), Point()

    p.insert_point_backward(p1, 0)
    p.insert_point_backward(p2, 1)

    assert p.backward_branches[0] is p1
    assert p.backward_branches[1] is p2
    assert list(p.forward_branches) == [None] * NUM_BRANCHES
    assert p1.forward is p and p2.forward is p
    assert p1.backward is None and p2.backward is None


def test_branch_point_insert_point_forward_connected() -> None:
    """Forward insertion on a branch splices between the root and its neighbour."""
    p = BranchPoint(NUM_BRANCHES)
    p1, p2 = Point(), Point()

    p.insert_point_forward(p2, 0)
    p.insert_point_forward(p1, 0)

    assert p.forward_branches[0] is p1
    assert (p1.forward, p2.forward) == (p2, None)
    assert (p2.backward, p1.backward) == (p1, p)
    assert p.backward_branches[0] is None


def test_branch_point_insert_point_backward_connected() -> None:
    """Backward insertion on a branch splices between the root and its neighbour."""
    p = BranchPoint(NUM_BRANCHES)
    p1, p2 = Point(), Point()

    p.insert_point_backward(p2, 0)
    p.insert_point_backward(p1, 0)

    assert p.backward_branches[0] is p1
    assert (p1.backward, p2.backward) == (p2, None)
    assert (p2.forward, p1.forward) == (p1, p)
    assert p.forward_branches[0] is None


# ------------------------------------------------ iterate list (both kinds)

@pytest.mark.parametrize("kind", sorted(POINT_KINDS))
def test_insert_next_iterate(kind: str) -> None:
    """``insert_next_iterate`` links both directions and leaves the ends open."""
    p1, p2 = POINT_KINDS[kind](), Point()

    p1.insert_next_iterate(p2)

    assert p1.next_iterate is p2
    assert p2.prev_iterate is p1
    assert p1.prev_iterate is None
    assert p2.next_iterate is None


@pytest.mark.parametrize("kind", sorted(POINT_KINDS))
def test_insert_prev_iterate(kind: str) -> None:
    """``insert_prev_iterate`` links both directions and leaves the ends open."""
    p1, p2 = POINT_KINDS[kind](), Point()

    p1.insert_prev_iterate(p2)

    assert p1.prev_iterate is p2
    assert p2.next_iterate is p1
    assert p1.next_iterate is None
    assert p2.prev_iterate is None


@pytest.mark.parametrize("kind", sorted(POINT_KINDS))
def test_insert_next_iterate_twice_raises(kind: str) -> None:
    """A point has at most one forward iterate."""
    p = POINT_KINDS[kind]()
    p.insert_next_iterate(Point())

    with pytest.raises(ValueError):
        p.insert_next_iterate(Point())


@pytest.mark.parametrize("kind", sorted(POINT_KINDS))
def test_insert_prev_iterate_twice_raises(kind: str) -> None:
    """A point has at most one backward iterate."""
    p = POINT_KINDS[kind]()
    p.insert_prev_iterate(Point())

    with pytest.raises(ValueError):
        p.insert_prev_iterate(Point())
