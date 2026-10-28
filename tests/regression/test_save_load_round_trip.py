"""
Regression (N10a): a saved workbench loads as the workbench it was, and goes on as it would.

Every law case is built once for this module, saved, loaded and re-saved: the
two saves are byte-equal (every array with its dtype, every JSON key), every
crossing keeps the type of each field and names a curve the loaded workbench
holds, and every slot holds the chain its branch runs on. A continuation
grows a fresh build and its loaded copy one step and recomputes both: the
crossings, the store and the bridges come out byte-equal. Segment ids are
drawn per index (``Tangle._next_sid``), so a load in a process that has built
other cases is as good as one in a fresh process.

No value is pinned: the original is compared with its loaded copy.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from cases import BUILDERS, LAW_CASES
from tanglepack import Intersection, TangleWorkbench


def _load_workbench(workbench, path):
    system = workbench.dynamical_system
    return TangleWorkbench.load(path, system.map, system.map_inv, system.jacobian)


def _assert_same_files(first, second):
    """Two saves are byte-equal: every array with its dtype, and the JSON."""
    with np.load(f"{first}.npz") as a, np.load(f"{second}.npz") as b:
        assert a.files == b.files
        for name in a.files:
            assert a[name].dtype == b[name].dtype and a[name].tobytes() == b[name].tobytes(), name
    with open(f"{first}.json") as a, open(f"{second}.json") as b:
        assert json.load(a) == json.load(b)


def _numbers(workbench) -> tuple:
    """The crossings' floats and ids, the iterate table, the bridges and the whole store."""
    registry = workbench.intersection_registry
    crossings = {
        crossing_id: (x.unstable_cdist, x.stable_cdist, x.coords and tuple(map(float, x.coords)),
                      x.resolution_ratio, x.crossing_sign, x.seg_ids and tuple(x.seg_ids))
        for crossing_id, x in registry
    }
    store = workbench.point_store
    columns = tuple(
        getattr(store, name)[: len(store)].tobytes()
        for name in ("_xy", "_cdist", "_image", "_preimage", "_chain", "_position")
    )
    return crossings, sorted(registry.iterate_table.items()), sorted(b.id for b in workbench.bridges), columns


@pytest.fixture(scope="module", params=LAW_CASES)
def saved_workbench(request, tmp_path_factory):
    """``(workbench, loaded copy, save path)`` of one law case."""
    workbench = BUILDERS[request.param]().workbench
    path = tmp_path_factory.mktemp(request.param) / "workbench"
    workbench.save(path)
    return workbench, _load_workbench(workbench, path), path


def test_a_loaded_workbench_resaves_byte_equal(saved_workbench):
    _, loaded, path = saved_workbench
    loaded.save(path.with_name("resaved"))
    _assert_same_files(path, path.with_name("resaved"))


def test_a_loaded_workbench_holds_the_same_numbers(saved_workbench):
    workbench, loaded, _ = saved_workbench
    assert _numbers(loaded) == _numbers(workbench)
    assert loaded.generation == workbench.generation


def test_a_loaded_crossing_keeps_its_types_and_its_curve(saved_workbench):
    workbench, loaded, _ = saved_workbench
    held = {id(curve) for curve in [*loaded.manifolds.values(), *loaded.bridges,
                                    *loaded._partial_bridges, *loaded._image_curves]}
    for (_, original), (_, copy) in zip(workbench.intersection_registry, loaded.intersection_registry):
        for field in vars(Intersection()):  # its fields, not the lazy brackets
            assert type(getattr(copy, field)) is type(getattr(original, field)), field
        if original.coords is not None:
            assert [type(c) for c in copy.coords] == [type(c) for c in original.coords]
        assert copy.unstable_manifold is None or id(copy.unstable_manifold) in held


def test_every_slot_holds_its_branch_chain(saved_workbench):
    _, loaded, _ = saved_workbench
    store = loaded.point_store
    for fixed_point in loaded.fixed_points:
        for (stability, orbit_index, branch), chain in fixed_point._slots.items():
            assert chain is store.chains[chain.index]
            assert chain.anchor is fixed_point.branch_points[orbit_index]
    for (fixed_point, stability, orbit_index, branch), manifold in loaded.manifolds.items():
        rows = manifold.ids
        if manifold.root is fixed_point.branch_points[orbit_index] and len(rows):
            assert store.chain_of(rows[0]) is fixed_point.slot_chain(stability, orbit_index, branch)


def _grow_and_recompute(workbench) -> tuple:
    for fixed_point in workbench.fixed_points:
        for stability in ("unstable", "stable"):
            workbench.grow_n_times(fixed_point, stability, 1)
    workbench.compute_intersections(workbench.fixed_points, preserve_ids=True, diagnose=True)
    workbench.rebuild_bridges()
    return _numbers(workbench)


@pytest.mark.parametrize("name", ["k10", "nested"])
def test_a_loaded_workbench_grows_as_the_original(name, tmp_path):
    workbench = BUILDERS[name]().workbench
    workbench.save(tmp_path / name)
    loaded = _load_workbench(workbench, tmp_path / name)
    assert _grow_and_recompute(loaded) == _grow_and_recompute(workbench)
