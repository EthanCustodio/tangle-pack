"""Saving and loading a workbench (N10a): the format, the file contents and the slot writer.

A save is plain numbers: every array is float64, int64 or int32, it is read
with ``allow_pickle=False``, and a save of another format is refused. The
segment index travels exactly: every record keeps its segment id, a crossing
pair keeps its iteration order (a pair's operand order), and the next id to
draw is the saved one.

A load maps a few linked rows forward and backward and refuses a map or an
inverse that does not give back their saved images (another Henon ``k``); a
save with no linked row is not checked, and the check maps a sample only.

Load writes every slot through :meth:`FixedPoint.attach_chain`, the one writer
of the slot table, and its input is now a file: a misspelt stability raises
instead of filling a slot no reader looks up.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from cases import build_k10
from helpers.fakes import bare_fixed_point
from tanglepack import TangleSession, TangleWorkbench
from tanglepack.examples import HENON_K10, henon_map, henon_map_inverse
from tanglepack.numerics import persist
from tanglepack.numerics.budget import PointBudgetExceeded
from tanglepack.numerics.PointStore import PointStore


def _chain():
    store = PointStore()
    return store.new_chain(store.append([(1.0, 0.0), (2.0, 0.0)], [1.0, 2.0]), 2.0)


def test_attach_chain_refuses_a_misspelt_stability():
    with pytest.raises(ValueError, match="'Unstable'"):
        bare_fixed_point().attach_chain(_chain(), "Unstable", 0, 0)


def test_attach_chain_refuses_a_third_branch():
    with pytest.raises(IndexError):
        bare_fixed_point().attach_chain(_chain(), "unstable", 0, 2)


def _load(workbench, path):
    system = workbench.dynamical_system
    return TangleWorkbench.load(path, system.map, system.map_inv, system.jacobian)


def test_a_save_of_another_format_is_refused(tmp_path):
    workbench = build_k10().workbench
    workbench.save(tmp_path / "k10")
    meta_path = tmp_path / "k10.json"
    meta = json.loads(meta_path.read_text())
    meta["format"] = 2
    meta_path.write_text(json.dumps(meta))
    with pytest.raises(ValueError, match="format 2.*format 1"):
        _load(workbench, tmp_path / "k10")


@pytest.mark.parametrize("version", [0, None])
def test_an_older_or_unversioned_save_is_refused(version, tmp_path):
    workbench = build_k10().workbench
    workbench.save(tmp_path / "k10")
    meta_path = tmp_path / "k10.json"
    meta = json.loads(meta_path.read_text())
    meta["format"] = version
    meta_path.write_text(json.dumps(meta))
    with pytest.raises(ValueError, match=f"format {version}.*format 1"):
        _load(workbench, tmp_path / "k10")


def test_the_state_the_law_cases_leave_at_rest_travels(tmp_path):
    # The law cases keep the default settings, abandon no bridge and leave no link
    # undecided, so their round trips cannot tell whether these travel.
    workbench = build_k10().workbench
    workbench.max_points_per_step = 10
    workbench.escape_radius = 7.5
    workbench.resolution_ratio_threshold = 2.5
    for bridge in workbench.bridges:
        try:
            workbench.iterate_bridge(bridge)
        except PointBudgetExceeded:
            break
    abandoned = {bridge.id for bridge in workbench._half_iterated_bridges}
    assert abandoned, "the budget must abandon a bridge for this test to mean anything"
    first, *rest = workbench.intersection_registry.resolved_ids(3.0)
    workbench.intersection_registry.unresolved_links = {first: tuple(rest[:2])}
    workbench.save(tmp_path / "k10")

    loaded = _load(workbench, tmp_path / "k10")
    settings = ("area_cutoff", "max_points_per_step", "escape_radius", "resolution_ratio_threshold")
    assert [getattr(loaded, name) for name in settings] == [getattr(workbench, name) for name in settings]
    assert loaded.intersection_registry.unresolved_links == workbench.intersection_registry.unresolved_links
    assert {bridge.id for bridge in loaded._half_iterated_bridges} == abandoned
    with pytest.raises(RuntimeError, match="half-refined"):
        loaded.iterate_bridge(loaded.bridge(next(iter(abandoned))))


def test_an_object_array_is_refused(tmp_path):
    workbench = build_k10().workbench
    workbench.save(tmp_path / "k10")
    with np.load(tmp_path / "k10.npz") as npz:
        arrays = dict(npz)
    arrays["store_cdist"] = np.array([{"cdist": 1.0}], dtype=object)
    np.savez(tmp_path / "k10.npz", **arrays)
    with pytest.raises(ValueError, match="allow_pickle"):
        _load(workbench, tmp_path / "k10")


def test_every_array_is_a_number_column_of_its_store_dtype(tmp_path):
    workbench = build_k10().workbench
    workbench.save(tmp_path / "k10")
    with np.load(tmp_path / "k10.npz", allow_pickle=False) as npz:
        dtypes = {name: npz[name].dtype for name in npz.files}
    assert set(dtypes.values()) <= {np.dtype(np.float64), np.dtype(np.int64), np.dtype(np.int32)}
    store = workbench.point_store
    assert all(dtypes[f"store{name}"] == getattr(store, name).dtype for name in persist._STORE_COLUMNS)


def test_the_segment_index_loads_with_its_ids_and_next_id(tmp_path):
    workbench = build_k10().workbench
    tangle = workbench.Tangle
    next_sid = tangle._next_sid
    workbench.save(tmp_path / "k10")
    assert tangle._next_sid == next_sid, "save draws no segment id"

    loaded = _load(workbench, tmp_path / "k10").Tangle
    assert loaded._next_sid == next_sid
    assert {sid: (s.row1, s.in_rtree) for sid, s in loaded._seg_lookup.items()} == {
        sid: (s.row1, s.in_rtree) for sid, s in tangle._seg_lookup.items()}
    assert loaded._edge_to_sid == tangle._edge_to_sid
    for pairs in ("_intersecting_segments", "_processed_pairs"):
        assert {tuple(p) for p in getattr(loaded, pairs)} == {tuple(p) for p in getattr(tangle, pairs)}


def test_a_pair_is_rebuilt_in_its_iteration_order():
    # Two ids that collide in a frozenset's table iterate in insertion order, so
    # both orders occur; the rebuild must give back the one that was saved.
    for a in range(64):
        for b in range(a):
            for saved in {tuple(frozenset((a, b))), tuple(frozenset((b, a)))}:
                assert tuple(persist._pair(*saved)) == saved


# ── the load-time map check ───────────────────────────────────────────────


@pytest.mark.parametrize("wrong", ["map", "map_inv"])
def test_a_load_with_another_map_is_refused(wrong, tmp_path):
    workbench = build_k10().workbench
    workbench.save(tmp_path / "k10")
    k, b = HENON_K10
    maps = {"map": henon_map(k, b), "map_inv": henon_map_inverse(k, b)}
    maps[wrong] = {"map": henon_map, "map_inv": henon_map_inverse}[wrong](k + 0.1, b)
    with pytest.raises(ValueError, match=f"the {wrong} given to load differs from the one the save was made with"):
        TangleWorkbench.load(tmp_path / "k10", maps["map"], maps["map_inv"])


def test_a_session_load_with_another_map_is_refused(tmp_path):
    session = build_k10().session
    session.save(tmp_path / "k10")
    k, b = HENON_K10
    with pytest.raises(ValueError, match="the map given to load differs"):
        TangleSession.load(tmp_path / "k10", henon_map(k + 0.1, b), henon_map_inverse(k, b))


def test_the_map_check_maps_a_sample_only(tmp_path):
    workbench = build_k10().workbench
    workbench.save(tmp_path / "k10")
    k, b = HENON_K10
    sizes = []

    def counted(inner):
        def apply(point):
            sizes.append(np.size(point) // 2)
            return inner(point)
        return apply

    TangleWorkbench.load(tmp_path / "k10", counted(henon_map(k, b)), counted(henon_map_inverse(k, b)))
    assert sizes and max(sizes) <= persist.MAP_CHECK_SAMPLE


def test_a_save_with_no_linked_row_loads_with_any_map(tmp_path):
    k, b = HENON_K10
    TangleWorkbench(henon_map(k, b), henon_map_inverse(k, b)).save(tmp_path / "empty")
    loaded = TangleWorkbench.load(tmp_path / "empty", henon_map(k + 1, b), henon_map_inverse(k + 1, b))
    assert not len(loaded.point_store)
