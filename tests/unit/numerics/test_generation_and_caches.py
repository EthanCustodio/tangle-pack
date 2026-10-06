"""Phase 4: the generation counter and the caches it invalidates.

Every derived view in the library (a Trellis snapshot, a memoised point array,
the registry's rank maps) is only valid for one *generation* of the tangle. This
module pins the three halves of that contract:

* :attr:`IntersectionRegistry.generation` and
  :attr:`TangleWorkbench.generation` advance on every mutation path and stand
  still for pure reads.
* the memoised point arrays agree with an un-memoised walk after growth.

Registry insertion staying near-linear is a wall-clock test, opt-in under the
``perf`` marker in ``tests/unit/numerics/test_registry_perf.py``.
"""

from __future__ import annotations


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
def _mutate_add(registry: IntersectionRegistry) -> None:
    """Store one new crossing."""
    registry.add(_ix(1.0, 2.0))


def _mutate_add_synthetic(registry: IntersectionRegistry) -> None:
    """Store one synthetic crossing."""
    registry.add_synthetic(coords=(0.0, 0.0), unstable_cdist=3.0, stable_cdist=4.0)


def _mutate_register_iterate(registry: IntersectionRegistry) -> None:
    """Link two stored crossings as one map step apart."""
    a = registry.add(_ix(1.0, 2.0))
    b = registry.add(_ix(2.0, 1.0))
    before = registry.generation
    registry.register_iterate(a, 1, b)
    assert registry.generation > before


def _mutate_reindex_from(registry: IntersectionRegistry) -> None:
    """Re-key a registry against an older one."""
    old = IntersectionRegistry()
    old.add(_ix(1.0, 2.0))
    registry.add(_ix(1.0, 2.0))
    before = registry.generation
    registry.reindex_from(old)
    assert registry.generation > before


@pytest.mark.parametrize(
    "mutate",
    [_mutate_add, _mutate_add_synthetic, _mutate_register_iterate, _mutate_reindex_from],
    ids=["add", "add_synthetic", "register_iterate", "reindex_from"],
)
def test_registry_generation_bumps_on_every_mutation(mutate) -> None:
    """Each registry mutation path advances ``generation``."""
    registry = IntersectionRegistry()
    before = registry.generation

    mutate(registry)

    assert registry.generation > before


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
def test_bridge_point_array_agrees_with_a_walk(henon_tangle_with_bridges):
    """Repeated reads of a bridge's point array agree with a fresh walk."""
    workbench, _fp = henon_tangle_with_bridges
    bridge = workbench.bridges[0]

    first = bridge.get_point_array()
    second = bridge.get_point_array()

    assert np.array_equal(first, second)
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
def test_branch_position_map_is_stable_across_calls(fixed_point):
    """Two reads of the position map agree and carry the point's ``k_value``."""
    _workbench, fp = fixed_point

    first, k_first = fp.branch_position_map("stable")
    second, k_second = fp.branch_position_map("stable")

    assert first == second and k_first == k_second == fp.k_value


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
