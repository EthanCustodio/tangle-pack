"""The merge's tie branches, each on a hand-built pair of curves.

``ManifoldMachine.merge_manifolds`` zips an existing run of points (the first
curve) with new images (the second) by canonical distance. Where two points
share a cdist it follows rules no case of the corpus reaches: (a) neither
point has an iterate link, (b) only the new one has, (c) the two ROOTS tie,
(d) one point is on both curves, (e) one curve ties with itself. Each test
pins the order the linked-list merge produced before the order arrays (N4),
and uses only ``Point(..., store=)``, the iterate setters,
``_manifold_from_points``, ``merge_manifolds`` and ``get_point_array``, so it
runs on both sides of that change.
"""

from __future__ import annotations

import pytest

from tanglepack import BaseManifold, Point
from tanglepack.examples.henon import henon_map, henon_map_inverse
from tanglepack.numerics import DynamicalSystem, ManifoldMachine
from tanglepack.numerics.PointStore import PointStore


@pytest.fixture
def machine() -> ManifoldMachine:
    """A machine; the merge never calls its map."""
    return ManifoldMachine(DynamicalSystem(henon_map(10), henon_map_inverse(10)))


@pytest.fixture
def store() -> PointStore:
    """A fresh store for each test's points."""
    return PointStore()


def _points(store: PointStore, cdists: list[float], y: float) -> list[Point]:
    """One point per cdist, on the line at height ``y``."""
    return [Point(c, y, c, store=store) for c in cdists]


def _curve(points: list[Point]) -> BaseManifold:
    """An unstable key-less curve through ``points``."""
    return ManifoldMachine._manifold_from_points(points, "unstable", 1.0, None)


def _merged(machine: ManifoldMachine, first: BaseManifold, second: BaseManifold) -> list:
    return machine.merge_manifolds(first, second).get_point_array(return_nodes=True)


def test_a_tie_without_links_keeps_the_existing_point(machine, store) -> None:
    """(a) Two points on one cdist, neither linked: the existing one stays."""
    old = _points(store, [0.0, 1.0, 2.0], 0.0)
    new = _points(store, [0.5, 1.0, 3.0], 1.0)

    merged = _merged(machine, _curve(old), _curve(new))

    assert merged == [old[0], new[0], old[1], old[2], new[2]]


def test_a_tie_where_only_the_new_point_is_linked_keeps_the_new_point(machine, store) -> None:
    """(b) Only the new point has an iterate link: it is kept and the existing
    point leaves the curve."""
    old = _points(store, [0.0, 1.0, 2.0], 0.0)
    new = _points(store, [0.5, 1.0, 3.0], 1.0)
    new[1].insert_prev_iterate(Point(0.0, 1.0, 0.25, store=store))

    merged = _merged(machine, _curve(old), _curve(new))

    assert merged == [old[0], new[0], new[1], old[2], new[2]]


def test_tied_roots_are_both_kept_the_existing_first(machine, store) -> None:
    """(c) The two roots share a cdist: the tie rule is not consulted, both
    stay, the existing root first."""
    old = _points(store, [1.0, 2.0], 0.0)
    new = _points(store, [1.0, 3.0], 1.0)

    merged = _merged(machine, _curve(old), _curve(new))

    assert merged == [old[0], new[0], old[1], new[1]]


def test_a_point_on_both_curves_appears_once(machine, store) -> None:
    """(d) Two curves over one list that share a point: the merge is their
    union, each point once."""
    points = _points(store, [0.0, 1.0, 2.0, 3.0, 4.0], 0.0)
    _curve(points)
    first = BaseManifold(points[0], "unstable", 1.0, None, tail=points[2], manifold_key=None)
    second = BaseManifold(points[2], "unstable", 1.0, None, tail=points[4], manifold_key=None)

    merged = _merged(machine, first, second)

    assert merged == points


def test_ties_within_one_curve_pair_one_to_one(machine, store) -> None:
    """(e) The new curve holds two points on a cdist the existing one holds
    once: the first pair ties (the existing point stays), the second new point
    is kept."""
    old = _points(store, [0.0, 1.0, 2.0], 0.0)
    new = _points(store, [1.0, 1.0, 3.0], 1.0)

    merged = _merged(machine, _curve(old), _curve(new))

    assert merged == [old[0], old[1], new[1], old[2], new[2]]
