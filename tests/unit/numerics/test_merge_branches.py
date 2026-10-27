"""The merge's tie branches, each on a hand-built pair of curves.

``ManifoldMachine.merge_manifolds`` zips an existing run of points (the first
curve) with new images (the second) by canonical distance. Where two points
share a cdist it follows rules no case of the corpus reaches: (a) neither
point has an iterate link, (b) only the new one has, (c) the two ROOTS tie,
(d) one point is on both curves, (e) one curve ties with itself. Each test
pins the order the linked-list merge produced before the order arrays (N4),
read as store rows (``PointStore.append`` / ``link``, ``_manifold_from_rows``,
``merge_manifolds`` and ``ids``). The last test merges onto a run headed
by a periodic point, the image of the bridge that leaves it, on a real tangle.
"""

from __future__ import annotations

import numpy as np
import pytest

from tanglepack import BaseManifold, BranchPoint
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


def _points(store: PointStore, cdists: list[float], y: float) -> list[int]:
    """One row per cdist, on the line at height ``y``."""
    return store.append([(c, y) for c in cdists], cdists).tolist()


def _curve(store: PointStore, rows: list[int]) -> BaseManifold:
    """An unstable key-less curve through ``rows``."""
    return ManifoldMachine._manifold_from_rows(store, rows, "unstable", 1.0, None)


def _merged(machine: ManifoldMachine, first: BaseManifold, second: BaseManifold) -> list[int]:
    return machine.merge_manifolds(first, second).ids.tolist()


def test_a_tie_without_links_keeps_the_existing_point(machine, store) -> None:
    """(a) Two points on one cdist, neither linked: the existing one stays."""
    old = _points(store, [0.0, 1.0, 2.0], 0.0)
    new = _points(store, [0.5, 1.0, 3.0], 1.0)

    merged = _merged(machine, _curve(store, old), _curve(store, new))

    assert merged == [old[0], new[0], old[1], old[2], new[2]]


def test_a_tie_where_only_the_new_point_is_linked_keeps_the_new_point(machine, store) -> None:
    """(b) Only the new point has an iterate link: it is kept and the existing
    point leaves the curve."""
    old = _points(store, [0.0, 1.0, 2.0], 0.0)
    new = _points(store, [0.5, 1.0, 3.0], 1.0)
    store.link([store.append_one(0.0, 1.0, 0.25)], [new[1]])

    merged = _merged(machine, _curve(store, old), _curve(store, new))

    assert merged == [old[0], new[0], new[1], old[2], new[2]]


def test_tied_roots_are_both_kept_the_existing_first(machine, store) -> None:
    """(c) The two roots share a cdist: the tie rule is not consulted, both
    stay, the existing root first."""
    old = _points(store, [1.0, 2.0], 0.0)
    new = _points(store, [1.0, 3.0], 1.0)

    merged = _merged(machine, _curve(store, old), _curve(store, new))

    assert merged == [old[0], new[0], old[1], new[1]]


def test_a_point_on_both_curves_appears_once(machine, store) -> None:
    """(d) Two curves over one list that share a point: the merge is their
    union, each point once."""
    points = _points(store, [0.0, 1.0, 2.0, 3.0, 4.0], 0.0)
    _curve(store, points)
    first = BaseManifold(points[0], "unstable", 1.0, None, tail=points[2], manifold_key=None, store=store)
    second = BaseManifold(points[2], "unstable", 1.0, None, tail=points[4], manifold_key=None, store=store)

    merged = _merged(machine, first, second)

    assert merged == points


def test_ties_within_one_curve_pair_one_to_one(machine, store) -> None:
    """(e) The new curve holds two points on a cdist the existing one holds
    once: the first pair ties (the existing point stays), the second new point
    is kept."""
    old = _points(store, [0.0, 1.0, 2.0], 0.0)
    new = _points(store, [1.0, 1.0, 3.0], 1.0)

    merged = _merged(machine, _curve(store, old), _curve(store, new))

    assert merged == [old[0], old[1], new[1], old[2], new[2]]


def test_an_anchor_bridge_image_is_its_branch_from_the_periodic_point(small_tangle) -> None:
    """The bridge leaving the periodic point maps onto the branch leaving it: its
    existing images (headed by the periodic point) and its new ones merge into
    the branch's own points, from the periodic point on (2026-10-25, N4: a path
    the order arrays had to keep, which no law case reaches)."""
    workbench, fp = small_tangle
    workbench.create_bridges(fp)
    (anchor_bridge,) = [bridge for bridge in workbench.bridges if isinstance(bridge.root, BranchPoint)]
    assert len(anchor_bridge.get_non_iterated_point_array()), "the merge must take new images"

    workbench.iterate_bridge(anchor_bridge)

    image = workbench._image_curves[-1]
    branch = workbench.manifolds[(fp, "unstable", 0, 0)]
    rows = image.ids.tolist()
    # The new images run on past the branch's tail, on the branch's own curve.
    branch_rows = branch.store.chain_of(branch.ids[0]).ids.tolist()
    assert image.root is fp.branch_points[0]
    assert rows == branch_rows[: len(rows)]
    assert len(rows) + 1 > len(branch.get_point_array())
    assert (np.diff(image.get_cdist_array(), axis=0) >= 0).all()
