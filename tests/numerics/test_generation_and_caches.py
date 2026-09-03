"""Phase 4: the generation counter and the caches it invalidates.

Every derived view in the library (a Trellis snapshot, a memoised point array,
the registry's rank maps) is only valid for one *generation* of the tangle. This
module pins the three halves of that contract:

* :attr:`IntersectionRegistry.generation` and
  :attr:`TangleWorkbench.generation` advance on every mutation path and stand
  still for pure reads.
* the memoised point arrays agree with an un-memoised walk after growth.
* registry insertion stays near-linear (the sorted key arrays plus the collision
  prefilter replace the old O(N) rescan per insert).
"""

from __future__ import annotations

import gc
import time

import numpy as np
import pytest

from tanglepack.numerics.Intersection import Intersection
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _ix(u: float, s: float, key=("u",), b_key=("s",)) -> Intersection:
    """A synthetic crossing at the given canonical distances."""
    return Intersection.synthetic(
        coords=(u, s),
        unstable_cdist=u,
        stable_cdist=s,
        manifold_a_key=key,
        manifold_b_key=b_key,
    )


def _walked_points(manifold) -> np.ndarray:
    """The point array computed by walking, bypassing any memo."""
    return np.vstack([node.get_point() for node in manifold._walk_nodes()])


# --------------------------------------------------------------------------- #
# registry generation
# --------------------------------------------------------------------------- #
def test_registry_generation_bumps_on_add():
    registry = IntersectionRegistry()
    before = registry.generation

    registry.add(_ix(1.0, 2.0))

    assert registry.generation > before


def test_registry_generation_bumps_on_add_synthetic():
    registry = IntersectionRegistry()
    before = registry.generation

    registry.add_synthetic(coords=(0.0, 0.0), unstable_cdist=3.0, stable_cdist=4.0)

    assert registry.generation > before


def test_registry_generation_bumps_on_register_iterate():
    registry = IntersectionRegistry()
    a = registry.add(_ix(1.0, 2.0))
    b = registry.add(_ix(2.0, 1.0))
    before = registry.generation

    registry.register_iterate(a, 1, b)

    assert registry.generation > before


def test_registry_generation_bumps_on_reindex_from():
    old = IntersectionRegistry()
    old.add(_ix(1.0, 2.0))
    fresh = IntersectionRegistry()
    fresh.add(_ix(1.0, 2.0))
    before = fresh.generation

    fresh.reindex_from(old)

    assert fresh.generation > before


def test_registry_generation_stands_still_on_a_deduped_add():
    """A collision stores nothing, so it is not a mutation."""
    registry = IntersectionRegistry()
    registry.add(_ix(1.0, 2.0))
    before = registry.generation

    registry.add(_ix(1.0, 2.0))

    assert registry.generation == before


def test_registry_generation_stands_still_on_reads():
    registry = IntersectionRegistry()
    first = registry.add(_ix(1.0, 2.0))
    registry.add(_ix(2.0, 1.0))
    before = registry.generation

    registry.by_unstable_cdist
    registry.by_stable_cdist
    registry.all_ids()
    registry.unstable_rank(first)
    registry.stable_rank(first)
    registry.nearest_by_unstable_cdist(1.4)
    registry.on_cdist_range(0.0, 10.0)
    registry.graph()
    registry[first]

    assert registry.generation == before


# --------------------------------------------------------------------------- #
# registry ordering views stay correct through the sorted key arrays
# --------------------------------------------------------------------------- #
def test_registry_orderings_and_ranks_track_insertions():
    registry = IntersectionRegistry()
    ids = [registry.add(_ix(u, 10.0 - u)) for u in (3.0, 1.0, 2.0)]

    assert [registry[i].unstable_cdist for i in registry.by_unstable_cdist] == [
        1.0,
        2.0,
        3.0,
    ]
    assert [registry[i].stable_cdist for i in registry.by_stable_cdist] == [
        7.0,
        8.0,
        9.0,
    ]
    assert registry.unstable_rank(ids[1]) == 0
    assert registry.unstable_rank(ids[0]) == 2
    assert registry.stable_rank(ids[0]) == 0
    assert registry.nearest_by_unstable_cdist(2.9) == ids[0]


def test_registry_insert_is_near_linear():
    """20k inserts complete quickly and the per-insert cost does not blow up."""
    registry = IntersectionRegistry()

    def add_block(start: int, count: int) -> float:
        t0 = time.perf_counter()
        for i in range(start, start + count):
            registry.add_synthetic(
                coords=(float(i), 0.0),
                unstable_cdist=float(i),
                stable_cdist=float(20000 - i),
            )
        return time.perf_counter() - t0

    # Fix for a PRE-EXISTING flake (it failed roughly one full-suite run in three
    # at Phase 5's HEAD, while always passing standalone; more tests in the suite
    # made it worse, which is what surfaced it). The ratio below compares two
    # wall-clock blocks, so it is only meaningful against a settled heap: run
    # inside the full suite, with several hundred other tests' garbage still
    # uncollected, the FIRST block is measured on a fragmented allocator and comes
    # out artificially fast, and the ratio blows past any bound. Collecting first
    # makes the two blocks comparable. The bound itself is unchanged.
    gc.collect()

    total_start = time.perf_counter()
    first_block = add_block(0, 1000)
    add_block(1000, 18000)
    last_block = add_block(19000, 1000)
    total = time.perf_counter() - total_start

    assert len(registry) == 20000
    assert total < 2.0, f"20k inserts took {total:.2f}s"
    # The bisect is O(log N) but the list insert behind it is O(N), so the
    # per-insert cost does grow with N by construction: this bound only has to
    # catch a return to the O(N) RESCAN (which was ~50x across this range), and
    # it has to survive a loaded machine, so it is deliberately loose.
    assert last_block < 15 * first_block, (
        f"per-insert cost grew from {first_block:.4f}s/1k to {last_block:.4f}s/1k"
    )


# --------------------------------------------------------------------------- #
# workbench generation
# --------------------------------------------------------------------------- #
MUTATIONS = {
    "register_manifold": lambda wb, fp: wb.register_manifold(
        (fp, "unstable", 0, 0), wb.manifolds[(fp, "unstable", 0, 0)]
    ),
    "grow_n_times": lambda wb, fp: wb.grow_n_times(fp, "unstable", num_iterations=1),
    "compute_intersections": lambda wb, fp: wb.compute_intersections([fp]),
    "create_bridges": lambda wb, fp: wb.create_bridges(fp),
    "clear_bridges": lambda wb, fp: wb.clear_bridges(),
    "rebuild_bridges": lambda wb, fp: wb.rebuild_bridges(fp),
    "iterate_bridge": lambda wb, fp: wb.iterate_bridge(wb.uniiterated_bridges[0]),
    "trim_stable_manifolds": lambda wb, fp: wb.trim_stable_manifolds(fp),
}


@pytest.mark.parametrize("name", sorted(MUTATIONS))
def test_workbench_generation_bumps_on_every_mutation(
    name, henon_tangle_with_bridges
):
    workbench, fp = henon_tangle_with_bridges
    before = workbench.generation

    MUTATIONS[name](workbench, fp)

    assert workbench.generation > before, f"{name} did not bump the generation"


def test_workbench_generation_stands_still_on_reads(henon_tangle_with_bridges):
    workbench, fp = henon_tangle_with_bridges
    bridge = workbench.bridges[0]
    before = workbench.generation

    workbench.bridges
    workbench.uniiterated_bridges
    workbench.manifolds
    workbench.intersection_registry
    workbench.bridges_at(bridge.first_intersection)
    workbench.bridge(bridge.id)
    workbench.image_bridges(bridge.id)
    workbench.manifolds[(fp, "unstable", 0, 0)].get_point_array()
    bridge.get_point_array()
    workbench.generation

    assert workbench.generation == before


def test_workbench_generation_bumps_when_a_manifold_tail_moves(
    henon_tangle_with_bridges,
):
    """Manifold geometry is part of the generation, not just the registry."""
    workbench, fp = henon_tangle_with_bridges
    manifold = workbench.manifolds[(fp, "stable", 0, 0)]
    before = workbench.generation

    manifold.tail = manifold.tail

    assert workbench.generation > before


# --------------------------------------------------------------------------- #
# memoised point arrays
# --------------------------------------------------------------------------- #
def test_bridge_point_array_is_memoised(henon_tangle_with_bridges):
    workbench, _fp = henon_tangle_with_bridges
    bridge = workbench.bridges[0]

    first = bridge.get_point_array()
    second = bridge.get_point_array()

    assert second is first, "a bridge's point array should be served from the memo"
    assert np.array_equal(first, _walked_points(bridge))


def test_point_arrays_reflect_growth(henon_tangle_with_bridges):
    """grow, read, grow, read: the second read must see the new points.

    Every bridge is checked, not a sample: growth refines the parent curve, so
    the points land inside whichever arcs happen to be curved enough, and a
    bridge that gained none would not notice a broken invalidation.
    """
    workbench, fp = henon_tangle_with_bridges
    manifold = workbench.manifolds[(fp, "unstable", 0, 0)]
    bridges = list(workbench.bridges)

    before_manifold = manifold.get_point_array()
    before = {id(b): len(b.get_point_array()) for b in bridges}
    assert np.array_equal(before_manifold, _walked_points(manifold))

    workbench.grow_n_times(fp, "unstable", num_iterations=1)

    after_manifold = manifold.get_point_array()
    assert len(after_manifold) > len(before_manifold)
    assert np.array_equal(after_manifold, _walked_points(manifold))
    for bridge in bridges:
        assert np.array_equal(bridge.get_point_array(), _walked_points(bridge))
    grew = [b for b in bridges if len(b.get_point_array()) > before[id(b)]]
    assert grew, "growth should have refined at least one already-cut bridge"


def test_node_walk_is_not_shared_between_callers(henon_tangle_with_bridges):
    """``return_nodes=True`` hands out a private list, never the memo's own."""
    workbench, _fp = henon_tangle_with_bridges
    bridge = workbench.bridges[0]

    nodes = bridge.get_point_array(return_nodes=True)
    length = len(nodes)
    nodes.append(nodes[0])

    assert len(bridge.get_point_array(return_nodes=True)) == length


# --------------------------------------------------------------------------- #
# per-fixed-point branch memo
# --------------------------------------------------------------------------- #
def test_branch_position_map_is_memoised(fixed_point):
    _workbench, fp = fixed_point

    first, k_first = fp.branch_position_map("stable")
    second, k_second = fp.branch_position_map("stable")

    assert first == second and k_first == k_second == fp.k_value
    assert first is second, "the position map should be served from the memo"


def test_branch_memo_is_invalidated_by_set_k_value(fixed_point):
    """Recomputing k_value (inversion appears) must not serve the stale cycle."""
    _workbench, fp = fixed_point

    assert len(fp.branch_cycle("stable")) == fp.k_value
    _positions, k_before = fp.branch_position_map("stable")

    fp.unstable_eigenvalues = [-abs(v) for v in fp.unstable_eigenvalues]
    fp.stable_eigenvalues = [-abs(v) for v in fp.stable_eigenvalues]
    fp.set_k_value()

    _positions, k_after = fp.branch_position_map("stable")
    assert k_before == 1 and k_after == 2
    assert len(fp.branch_cycle("stable")) == 2


def test_bridge_point_arrays_survive_an_iterate(henon_tangle_with_bridges):
    """Iterating drops every registered bridge's memo, not just the child's.

    The image of a bridge whose points already have iterates is merged and
    refined over curve other bridges are cut from, and the cut of the image
    splices separator points into that same shared list, so an iterate can lay
    points inside an already-cut arc. It does not happen on this fixture (every
    image here is fresh curve), which is exactly why the contract is pinned
    directly -- every registered bridge's version must move, and every
    memoised array must still agree with a fresh walk.
    """
    workbench, _fp = henon_tangle_with_bridges
    registry = workbench.intersection_registry
    watched = list(workbench.bridges)
    versions = {id(b): b.version for b in watched}
    for bridge in watched:
        bridge.get_point_array()  # populate the memos

    # The OUTERMOST bridge: its image runs past the grown extent, so it is
    # really mapped forward and cut. An inner bridge's image is already-cut
    # curve, which iterate_bridge reuses without touching any geometry.
    outermost = max(
        (b for b in workbench.uniiterated_bridges if not b.partial),
        key=lambda b: registry[b.second_intersection].unstable_cdist,
    )
    workbench.iterate_bridge(outermost)

    for bridge in watched:
        assert bridge.version > versions[id(bridge)], (
            f"bridge {bridge.id} kept its memo across an iterate"
        )
        assert np.array_equal(bridge.get_point_array(), _walked_points(bridge))
