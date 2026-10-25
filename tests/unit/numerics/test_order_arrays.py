"""The order arrays: points along a curve are the rows of a chain.

A :class:`~tanglepack.numerics.PointStore.Chain` holds the rows of one physical
curve in curve order, and two store columns say which chain a row is on and
where. The invariants are checked on a grown k=10 tangle, after a blast, and
after growth after a blast (which merges blast images into the branches):
cdist never decreases along a chain, a row lies on one chain at most and once,
the columns agree with every chain, every curve is a view of one chain from
root to tail, and the iterate links of the points with a cdist are a partial
bijection.

The merge states no case reaches and one chain cannot hold, (f)-(h), raise.
The iterate list both point kinds inherit is tested here too.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
import pytest

from tanglepack import BaseManifold, BranchPoint, Point
from tanglepack.examples.henon import henon_map, henon_map_inverse
from tanglepack.numerics import DynamicalSystem, ManifoldMachine
from tanglepack.numerics.BasePoint import BasePoint
from tanglepack.numerics.PointStore import PointStore

NUM_BRANCHES = 2

#: Factories for the two point kinds the iterate list is shared by.
POINT_KINDS: dict[str, Callable[[], BasePoint]] = {
    "point": lambda: Point(),
    "branch_point": lambda: BranchPoint(NUM_BRANCHES),
}


def assert_chains_hold(workbench) -> None:
    """Every invariant of the order arrays on a workbench's store and curves."""
    store = workbench.point_store
    n = len(store)
    chained = np.concatenate([chain.ids for chain in store.chains])
    assert len(np.unique(chained)) == len(chained), "a row lies on two chains, or twice on one"
    assert (store._chain[:n] >= 0).sum() == len(chained), "the chain column names an unchained row"
    for chain in store.chains:
        assert (store._chain[chain.ids] == chain.index).all()
        np.testing.assert_array_equal(store._position[chain.ids], np.arange(len(chain)))
        assert (np.diff(store.cdist[chain.ids]) >= 0).all(), "cdist decreases along a chain"
    # The initializer's detached fictitious tails (no cdist) keep a one-way link
    # to the seed they were cut from; every other link is one-to-one.
    has_cdist = ~np.isnan(store.cdist)
    for column in (store.image, store.preimage):
        linked = column[(column >= 0) & has_cdist]
        assert len(np.unique(linked)) == len(linked), "two rows share an image or a preimage"

    curves = [
        *workbench.manifolds.values(),
        *workbench.bridges,
        *workbench._image_curves,
    ]
    for curve in curves:
        _head, chain, lo, hi = curve._span()
        tail = curve.tail
        assert chain is not None and store.chain_of(tail.row) is chain
        assert lo <= store.position_of(tail.row) == hi - 1


@pytest.fixture
def blasted(small_tangle):
    """The 7-step k=10 tangle, cut and blasted once."""
    workbench, fp = small_tangle
    workbench.create_bridges(fp)
    workbench.iterate_all_bridges()
    return workbench, fp


def test_a_grown_tangle_holds_the_invariants(small_tangle) -> None:
    """Growth, refinement and the crossing cut keep the chains exact."""
    workbench, _fp = small_tangle
    assert_chains_hold(workbench)


def test_a_blast_holds_the_invariants(blasted) -> None:
    """The blast's image curves are chains (or runs of the branches) like any other."""
    workbench, _fp = blasted
    assert workbench._image_curves, "the fixture must blast"
    assert_chains_hold(workbench)


def test_growth_after_a_blast_holds_the_invariants(blasted) -> None:
    """Growth merges blast images into the branches: their chains are absorbed."""
    workbench, fp = blasted
    chains = sum(1 for chain in workbench.point_store.chains if len(chain))
    workbench.grow_n_times(fp, "unstable", num_iterations=1)
    workbench.compute_intersections([fp])
    workbench.rebuild_bridges()

    assert_chains_hold(workbench)
    assert sum(1 for chain in workbench.point_store.chains if len(chain)) <= chains


# --------------------------------------------------------- merges one chain cannot hold

@pytest.fixture
def machine() -> ManifoldMachine:
    """A machine; the merge never calls its map."""
    return ManifoldMachine(DynamicalSystem(henon_map(10), henon_map_inverse(10)))


def _curve(store: PointStore, cdists: list[float], y: float) -> BaseManifold:
    """An unstable key-less curve through one new point per cdist."""
    points = [Point(c, y, c, store=store) for c in cdists]
    return ManifoldMachine._manifold_from_points(points, "unstable", 1.0, None)


def _run(curve: BaseManifold, first: int, last: int) -> BaseManifold:
    """The run of ``curve``'s points from index ``first`` to ``last``."""
    nodes = curve.get_point_array(return_nodes=True)
    return BaseManifold(nodes[first], "unstable", 1.0, None, tail=nodes[last], manifold_key=None)


def test_a_run_that_runs_out_before_its_chain_ends_raises(machine) -> None:
    """(f) The existing run runs out first but its chain goes on past its tail."""
    store = PointStore()
    existing = _run(_curve(store, [0.0, 1.0, 2.0], 0.0), 0, 1)

    with pytest.raises(AssertionError, match="past its tail"):
        machine.merge_manifolds(existing, _curve(store, [1.5, 3.0], 1.0))


def test_a_lower_new_root_before_a_mid_chain_run_raises(machine) -> None:
    """(g) The new root comes first while the existing run starts mid-chain."""
    store = PointStore()
    existing = _run(_curve(store, [0.0, 1.0, 2.0], 0.0), 1, 2)

    with pytest.raises(AssertionError, match="beginning of its chain"):
        machine.merge_manifolds(existing, _curve(store, [0.5, 3.0], 1.0))


def test_runs_that_both_end_on_a_tie_raise(machine) -> None:
    """(h) The last points of both runs share a cdist."""
    store = PointStore()

    with pytest.raises(AssertionError, match="both runs end on a tie"):
        machine.merge_manifolds(_curve(store, [0.0, 2.0], 0.0), _curve(store, [1.0, 2.0], 1.0))


# ------------------------------------------------------------- points and branch points

def test_point_creation() -> None:
    """A fresh ``Point`` holds its coordinates and lies on no chain."""
    p = Point(1.0, 2.0)

    assert p.x == 1.0
    assert p.y == 2.0
    assert p.store.chain_of(p.row) is None and p.stretch_param is None
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


@pytest.mark.parametrize("direction", ["forward", "backward"])
def test_a_branch_point_hangs_a_chain_on_a_slot(direction: str) -> None:
    """``insert_point_*`` fills one slot and makes the branch point the chain's
    anchor, the left neighbour of its first position; the other slots stay empty."""
    p = BranchPoint(NUM_BRANCHES)
    stability = "unstable" if direction == "forward" else "stable"
    curve = ManifoldMachine._manifold_from_points([Point(), Point()], stability, 2.5, None)
    first = curve.root

    getattr(p, f"insert_point_{direction}")(first, 1)

    slots = p.forward_branches if direction == "forward" else p.backward_branches
    others = p.backward_branches if direction == "forward" else p.forward_branches
    chain = first.store.chain_of(first.row)
    assert list(slots) == [None, first] and list(others) == [None] * NUM_BRANCHES
    assert chain.anchor is p and chain.node_at(-1) is p and first.stretch_param == 2.5


def test_a_slot_holding_another_chain_refuses_a_second() -> None:
    """A branch is one chain: hanging another chain's point on its slot raises."""
    p = BranchPoint(NUM_BRANCHES)
    first = ManifoldMachine._manifold_from_points([Point()], "unstable", 1.0, None).root
    other = ManifoldMachine._manifold_from_points([Point()], "unstable", 1.0, None).root
    p.insert_point_forward(first, 0)

    with pytest.raises(ValueError):
        p.insert_point_forward(other, 0)
    assert p.forward_branches[0] is first


@pytest.mark.parametrize("stability", ["unstable", "stable"])
def test_each_node_follows_its_chain_from_the_anchor(initialized, stability):
    """Read back along the chain, every node's predecessor is the node before it
    on the curve, and the first point's is the root BranchPoint, the chain's
    anchor."""
    workbench, fp = initialized
    manifold = workbench.manifolds[(fp, stability, 0, 0)]
    nodes = manifold.get_point_array(return_nodes=True)
    first = nodes[1]
    store = first.store
    chain = store.chain_of(first.row)

    assert nodes[0] is fp.branch_points[0] is chain.anchor
    for previous, node in zip(nodes, nodes[1:]):
        assert chain.node_at(store.position_of(node.row) - 1) is previous
    assert chain.node_at(len(chain)) is None


# ------------------------------------------------ iterate list (both kinds)

@pytest.mark.parametrize("kind", sorted(POINT_KINDS))
def test_insert_next_iterate(kind: str) -> None:
    """``insert_next_iterate`` links both directions and leaves the ends open."""
    p1, p2 = POINT_KINDS[kind](), POINT_KINDS[kind]()

    p1.insert_next_iterate(p2)

    assert p1.next_iterate is p2
    assert p2.prev_iterate is p1
    assert p1.prev_iterate is None
    assert p2.next_iterate is None


@pytest.mark.parametrize("kind", sorted(POINT_KINDS))
def test_insert_prev_iterate(kind: str) -> None:
    """``insert_prev_iterate`` links both directions and leaves the ends open."""
    p1, p2 = POINT_KINDS[kind](), POINT_KINDS[kind]()

    p1.insert_prev_iterate(p2)

    assert p1.prev_iterate is p2
    assert p2.next_iterate is p1
    assert p1.next_iterate is None
    assert p2.prev_iterate is None


@pytest.mark.parametrize("kind", sorted(POINT_KINDS))
def test_insert_next_iterate_twice_raises(kind: str) -> None:
    """A point has at most one forward iterate."""
    p = POINT_KINDS[kind]()
    p.insert_next_iterate(POINT_KINDS[kind]())

    with pytest.raises(ValueError):
        p.insert_next_iterate(POINT_KINDS[kind]())


@pytest.mark.parametrize("kind", sorted(POINT_KINDS))
def test_insert_prev_iterate_twice_raises(kind: str) -> None:
    """A point has at most one backward iterate."""
    p = POINT_KINDS[kind]()
    p.insert_prev_iterate(POINT_KINDS[kind]())

    with pytest.raises(ValueError):
        p.insert_prev_iterate(POINT_KINDS[kind]())
