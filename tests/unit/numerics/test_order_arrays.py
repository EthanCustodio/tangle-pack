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
A point laid right after a periodic point (a crossing separator, a refined
pair at the anchor) becomes position 0 of the chain its slot holds.
The periodic point, a view of its fixed point's orbit, is tested here too.
"""

from __future__ import annotations

import numpy as np
import pytest

from helpers.fakes import bare_fixed_point
from helpers.polylines import periodic_point, polyline_manifold
from tanglepack import BaseManifold, BranchPoint
from tanglepack.examples.henon import henon_map, henon_map_inverse
from tanglepack.numerics import DynamicalSystem, ManifoldMachine
from tanglepack.numerics.PointStore import PointStore
from tanglepack.numerics.Tangle import Tangle

NUM_BRANCHES = 2


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
        tail = curve._tail
        assert chain is not None and store.chain_of(tail) is chain
        assert lo <= store.position_of(tail) == hi - 1


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


def test_ids_are_the_rows_from_root_to_tail(blasted) -> None:
    """``ids`` is a fresh array of the rows a curve reads from its root to its
    tail: cdist never decreases along it, no row is in it twice, and a
    periodic-point root, which has no row, is not in it."""
    workbench, _fp = blasted
    store = workbench.point_store
    curves = [*workbench.manifolds.values(), *workbench.bridges, *workbench._image_curves]
    assert any(isinstance(curve.root, BranchPoint) for curve in curves)
    for curve in curves:
        ids = curve.ids
        head, _store, rows = curve._rows()
        np.testing.assert_array_equal(ids, rows)
        assert ids.dtype == np.int64 and not np.shares_memory(ids, rows)
        assert len(np.unique(ids)) == len(ids)
        assert (np.diff(store.cdist[ids]) >= 0).all()
        assert len(ids) == len(curve.points) - (head is not None)


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
    """An unstable key-less curve through one new row per cdist."""
    rows = store.append([(c, y) for c in cdists], cdists)
    return ManifoldMachine._manifold_from_rows(store, rows, "unstable", 1.0, None)


def _run(curve: BaseManifold, first: int, last: int) -> BaseManifold:
    """The run of ``curve``'s points from index ``first`` to ``last``."""
    rows = curve.ids
    return BaseManifold(int(rows[first]), "unstable", 1.0, None, tail=int(rows[last]),
                        manifold_key=None, store=curve.store)


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


def test_a_new_run_that_runs_out_before_its_chain_ends_raises(machine) -> None:
    """(f) for the other run: an earlier merge's output, a run whose chain goes on
    past its tail, runs out first."""
    store = PointStore()
    earlier = _run(_curve(store, [1.0, 2.0, 3.0], 1.0), 0, 1)

    with pytest.raises(AssertionError, match="past its tail"):
        machine.merge_manifolds(_curve(store, [0.0, 5.0], 0.0), earlier)


def test_a_mid_chain_new_run_placed_second_raises(machine) -> None:
    """(g) for the other run: the existing root comes first while the new run
    starts in the middle of its chain."""
    store = PointStore()
    later = _run(_curve(store, [0.5, 1.0, 3.0], 1.0), 1, 2)

    with pytest.raises(AssertionError, match="beginning of its chain"):
        machine.merge_manifolds(_curve(store, [0.0, 2.0], 0.0), later)


def test_runs_that_both_end_on_a_tie_raise(machine) -> None:
    """(h) The last points of both runs share a cdist."""
    store = PointStore()

    with pytest.raises(AssertionError, match="both runs end on a tie"):
        machine.merge_manifolds(_curve(store, [0.0, 2.0], 0.0), _curve(store, [1.0, 2.0], 1.0))


# ------------------------------------------------------------- points and branch points

def test_branch_point_creation() -> None:
    """The periodic point of a hand-built period-1 fixed point is its own image,
    reads cdist 0.0 and has one empty slot per branch in each direction."""
    p = periodic_point(1.0, 2.0)

    assert p.x == 1.0
    assert p.y == 2.0
    assert p.next_iterate is p and p.prev_iterate is p
    assert p.get_cdist("unstable") == p.get_cdist("stable") == 0.0
    assert len(p.forward_branches) == NUM_BRANCHES
    assert len(p.backward_branches) == NUM_BRANCHES
    np.testing.assert_array_equal(p.get_point(), np.array([1.0, 2.0]))


@pytest.mark.parametrize("direction", ["forward", "backward"])
def test_a_branch_point_hangs_a_chain_on_a_slot(direction: str) -> None:
    """``attach`` puts the chain in one slot and makes the branch point
    its anchor, the left neighbour of its first position; the other slots stay
    empty."""
    p = periodic_point()
    stability = "unstable" if direction == "forward" else "stable"
    store = PointStore()
    curve = ManifoldMachine._manifold_from_rows(store, store.append([(0.0, 0.0)] * 2, None), stability, 2.5, None)
    chain = curve.store.chain_of(curve.ids[0])

    p.attach(chain, stability, 1)

    slots = p.forward_branches if direction == "forward" else p.backward_branches
    others = p.backward_branches if direction == "forward" else p.forward_branches
    assert list(slots) == [None, chain] and list(others) == [None] * NUM_BRANCHES
    assert chain.anchor is p and chain.stretch_param == 2.5
    np.testing.assert_array_equal(chain.xy_at(-1), p.get_point())


def test_a_slot_holding_another_chain_refuses_a_second() -> None:
    """A branch is one chain: hanging another chain on its slot raises."""
    p = periodic_point()
    store = PointStore()
    first, other = (
        ManifoldMachine._manifold_from_rows(store, store.append([(0.0, 0.0)], None), "unstable", 1.0, None)
        for _ in range(2)
    )
    first_chain = first.store.chain_of(first.ids[0])
    p.attach(first_chain, "unstable", 0)

    with pytest.raises(ValueError):
        p.attach(other.store.chain_of(other.ids[0]), "unstable", 0)
    assert p.forward_branches[0] is first_chain


def test_a_periodic_point_reads_its_orbit() -> None:
    """A periodic point is a view of its fixed point: its image and preimage are
    the next and previous orbit points (under inversion too: the branch flips,
    the point does not), its coordinates are the orbit point's and its cdist
    is 0.0 on both manifolds."""
    for period, inversion in ((3, False), (2, True)):
        orbit = [(float(o), 10.0 + o) for o in range(period)]
        fp = bare_fixed_point(period, inversion=inversion, coordinates=orbit)
        assert fp.num_branches == (2 if inversion else 1)

        for o, bp in enumerate(fp.branch_points):
            assert bp.next_iterate is fp.branch_points[(o + 1) % period]
            assert bp.prev_iterate is fp.branch_points[(o - 1) % period]
            assert (bp.x, bp.y) == orbit[o]
            np.testing.assert_array_equal(bp.get_point(), fp.coordinates[o])
            assert bp.get_cdist("unstable") == bp.get_cdist("stable") == 0.0


def test_a_periodic_point_cannot_be_written_through() -> None:
    """The point's coordinates are the fixed point's: ``get_point()`` is a
    read-only view, and a cdist is read on a named manifold only."""
    p = periodic_point(1.0, 2.0)

    with pytest.raises(ValueError):
        p.get_point()[0] = 5.0
    with pytest.raises(ValueError):
        p.get_cdist("neither")
    assert (p.x, p.y) == (1.0, 2.0)


@pytest.mark.parametrize("stability", ["unstable", "stable"])
def test_each_node_follows_its_chain_from_the_anchor(initialized, stability):
    """Read back along the chain, every row's predecessor is the row before it
    on the curve, and the first row's is the root BranchPoint, the chain's
    anchor."""
    workbench, fp = initialized
    manifold = workbench.manifolds[(fp, stability, 0, 0)]
    rows = manifold.ids
    store = manifold.store
    chain = store.chain_of(rows[0])

    assert manifold.root is fp.branch_points[0] is chain.anchor
    np.testing.assert_array_equal(store._position[rows], np.arange(len(rows)))
    np.testing.assert_array_equal(chain.ids[: len(rows)], rows)
    np.testing.assert_array_equal(chain.xy_at(-1), chain.anchor.get_point())
    assert chain.xy_at(len(chain)) is None


# ------------------------------------------------ a point right after the periodic point

def _assert_first_on_its_chain(periodic: BranchPoint, first: int, curve: BaseManifold) -> None:
    """Row ``first`` is position 0 of the chain ``periodic``'s slot holds, its
    left neighbour (position -1, which the resolution ratio's chord deviation
    reads) is the periodic point, and the curve reads it right after it."""
    chain = periodic.forward_branches[0]
    store = chain.store
    assert chain.anchor is periodic and store.chain_of(first) is chain
    assert store.position_of(first) == 0
    np.testing.assert_array_equal(chain.xy_at(-1), periodic.get_point())
    np.testing.assert_array_equal(store._position[chain.ids], np.arange(len(chain)))
    assert curve.root is periodic and curve.ids[0] == first


def test_a_separator_after_the_periodic_point_becomes_the_first_point() -> None:
    """A stable V crossing the segment that leaves the periodic point twice: the
    cut lays a separator between the two crossings, at position 0 of the chain."""
    fp = bare_fixed_point(label="A")
    periodic = periodic_point()
    tangle = Tangle()
    tangle.add_manifolds([polyline_manifold([(0.3, 1.0), (0.5, -1.0), (0.7, 1.0)], (fp, "stable", 0, 0))])
    curve = polyline_manifold([(1.0, 0.0), (2.0, 0.0)], (fp, "unstable", 0, 0), root=periodic)
    chain = periodic.forward_branches[0]
    seed = int(chain.ids[0])
    tangle.add_manifold(curve, index_segments=False)
    crossings = tangle.resolve_crossings()
    for i, crossing in enumerate(crossings):
        crossing.id = i

    (bridge,) = tangle.create_bridges(crossings, for_manifold=curve)

    separator = int(chain.ids[0])
    assert periodic.forward_branches[0] is chain
    assert separator != seed and chain.store.xy[separator, 0] == 0.5 and chain.ids[1] == seed
    assert bridge.root is periodic and bridge.ids.tolist() == [separator, seed]
    _assert_first_on_its_chain(periodic, separator, curve)


def test_refining_the_pair_at_the_periodic_point_lays_the_new_first_point() -> None:
    """The pair (periodic point, first point), refinable when the first point has
    a preimage: the new point goes to position 0 and the slot keeps its chain. Its
    left neighbour, the periodic point's, is no point (NaN); the base read the
    stable branch's first point there and raised before it could insert."""
    periodic = periodic_point()
    store = PointStore()
    xs = np.array([0.1, 0.2, 0.4, 0.8])
    rows = store.append(np.column_stack((xs, xs * xs)), xs)
    store.link(store.append(np.column_stack((xs / 2, xs * xs / 2)), xs / 2), rows)
    ManifoldMachine._manifold_from_rows(store, rows, "unstable", 2.0, None)
    periodic.attach(store.chain_of(rows[0]), "unstable", 0)
    curve = BaseManifold(periodic, "unstable", 2.0, None, tail=int(rows[-1]), branch_index=0,
                         manifold_key=None, store=store)
    double = lambda p: 2 * np.asarray(p)  # noqa: E731 - the map of the preimages above
    half = lambda p: 0.5 * np.asarray(p)  # noqa: E731
    machine = ManifoldMachine(DynamicalSystem(double, half), area_cutoff=1e-5)

    chain = periodic.forward_branches[0]

    machine.refine_manifold(curve)

    first = int(chain.ids[0])
    assert periodic.forward_branches[0] is chain and first != rows[0]
    assert 0.0 < chain.store.cdist[first] < xs[0]
    _assert_first_on_its_chain(periodic, first, curve)
    assert (np.diff(curve.get_cdist_array(), axis=0) >= 0).all()


@pytest.mark.parametrize("root", ["point", "empty_slot"])
def test_a_curve_of_one_node_reads_one_row(root: str) -> None:
    """A lone point, or a periodic point whose slot is empty, reads as one row in
    every getter. The lone point has no image; the periodic point always has
    one, the next orbit point."""
    store = PointStore()
    head = store.append_one(1.0, 2.0, 0.5) if root == "point" else periodic_point()
    curve = BaseManifold(head, "unstable", 1.0, None, tail=head, branch_index=0, manifold_key=None,
                         store=store)

    assert curve.get_point_array().shape == (1, 2)
    assert curve.get_cdist_array().shape == (1, 1)
    if root == "point":
        assert curve.get_non_iterated_point_array().shape == (1, 2)
        assert curve.get_iterated_point_array().size == 0
    else:  # the periodic point always has its image: the next orbit point
        assert curve.get_non_iterated_point_array().size == 0
        assert curve.get_iterated_point_array().shape == (1, 2)
    assert curve.ids.tolist() == ([head] if root == "point" else [])


# ------------------------------------------------ how a curve takes its ends

def test_a_periodic_point_root_takes_a_row_tail_without_store() -> None:
    """A curve rooted at a periodic point adopts its slot chain's store, so a
    row tail needs no ``store=``; the slot is read through ``branch_index``,
    which is set before the ends."""
    periodic = periodic_point()
    store = PointStore()
    rows = store.append([(0.1, 0.0), (0.2, 0.0), (0.3, 0.0)], [0.1, 0.2, 0.3])
    periodic.attach(store.new_chain(rows, 2.0), "unstable", 0)

    curve = BaseManifold(periodic, "unstable", 2.0, None, tail=int(rows[1]), branch_index=0,
                         manifold_key=None)

    assert curve.store is store
    assert curve.ids.tolist() == rows[:2].tolist()


def test_a_non_integer_end_raises() -> None:
    """An end is a row (an integer), an ``End`` or the periodic point: a float
    is no row and is not truncated to one."""
    store = PointStore()
    rows = store.append([(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)], [0.0, 1.0, 2.0])
    store.new_chain(rows, 1.0)

    with pytest.raises(TypeError):
        BaseManifold(int(rows[0]), "unstable", 1.0, None, tail=2.7, manifold_key=None, store=store)


def test_a_row_outside_the_store_raises() -> None:
    """A row end is a row of the curve's store: one past its last row, or the
    -1 sentinel, raises rather than reading spare capacity."""
    store = PointStore()
    rows = store.append([(0.0, 0.0), (1.0, 0.0)], [0.0, 1.0])
    store.new_chain(rows, 1.0)

    for row in (len(store), -1):
        with pytest.raises(IndexError):
            BaseManifold(int(rows[0]), "unstable", 1.0, None, tail=row, manifold_key=None, store=store)


def test_a_row_end_without_a_store_raises() -> None:
    """A row means nothing without its store: a row root with no ``store=`` (and
    no periodic point to adopt one from) raises."""
    with pytest.raises(ValueError):
        BaseManifold(0, "unstable", 1.0, None, tail=1, manifold_key=None)
