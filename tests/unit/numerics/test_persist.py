"""Saving and loading a workbench (N10a): the format, the file contents and the slot writer.

A save is plain numbers: every array is float64, int64 or int32, it is read
with ``allow_pickle=False``, and a save of another format is refused. The
segment index travels exactly: every record keeps its segment id, a crossing
pair keeps its iteration order (a pair's operand order), and the next id to
draw is the saved one.

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
from tanglepack import TangleWorkbench
from tanglepack.numerics import persist
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
