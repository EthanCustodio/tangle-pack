"""
Saving and loading a workbench: one ``.npz`` of arrays and one ``.json`` of structure.

``save(path)`` writes ``<path>.npz`` (plain ``np.savez``, uncompressed) and
``<path>.json``. Every float lives in the npz; the JSON holds structure, names,
ints and ``None`` -- no float is ever written as text. A loaded workbench is the
saved one: its re-save is byte-equal, and the growth, recomputes and blasts that
follow give the bits the original's would (N10a, 2026-10-27; measured on the law
cases, the blast series and the scale series).

What is saved is the workbench's whole state: lists and dicts in their insertion
orders, sets sorted (a set's iteration order is its hash table's, not state):

* the :class:`~.PointStore.PointStore` (written rows, capacity) and every chain,
  retired ones included;
* the fixed points with their eigendata and slot tables, each array field in
  the container it was set as (the solver leaves ``coordinates`` a
  ``(period, 2)`` ndarray, the other fields are lists of rows);
* every curve once, by object identity (a crossing's ``unstable_manifold`` is
  compared with ``is`` downstream); a bridge's ``width`` and ``width_saddle``
  (the width guard's) as columns written only when some bridge has a width,
  so a save without one is the save of format 1 before the guard existed;
* every :class:`~.Intersection.Intersection` once, the registry's first, in one
  table, with the registry's orders, iterate table and counters;
* the segment index EXACTLY: every record's segment id, the edge maps, the
  crossing pairs, each in its operand order (the sets themselves sorted), and
  the next id to draw. A pair's operand order follows its two ids and decides
  float bits, so an index rebuilt from scratch would not continue the same (the
  N10a analyst's ablation: 51 of 59 resumed series steps moved in the last bits);
* the bookkeeping later calls read: growth stats, reported findings, epochs,
  the generation, the machine's settings.

What is NOT saved is derived and rebuilt lazily: the registry's caches,
``Intersection.brackets``, the fixed points' memos, the system's batchability
probe (made on the map passed to ``load``) and the rtrees, which are
bulk-loaded from the saved records.

A lean workbench (``TangleWorkbench.lean``) is never saved: the curve a
lean blast freed cannot be restored, so the save raises ``CurveReleased``.

A load runs nothing it reads: the npz is opened with
``allow_pickle=False``, and the map, its inverse and the Jacobian are arguments
of ``load`` only. A file of another ``"format"`` is refused, and so is a map
that is not the save's: once the state is restored, :func:`check_map` maps a few
linked rows forward (and a few backward) and compares them with their saved
images (preimages). The check only reads the restored state.

Dev Notes:

Each saved number records its Python type (``float``, ``np.float64`` or ``int``): a
detected crossing carries ``np.float64`` coordinates and an anchor ``float``
ones, and texts quote their reprs. A periodic-point end (a
:class:`~.BranchPoint.BranchPoint`) is two ints of its own, ``(fixed point,
orbit index)``, beside the row column.

The session's remembered state (pips, letters, zones) is saved beside the
workbench by ``TangleSession.save``, through :class:`Encoder` /
:class:`Decoder` and :func:`write`; nothing here imports the layers above.
"""

from __future__ import annotations

import json
import os
from collections import defaultdict
from importlib import metadata
from typing import TYPE_CHECKING, Iterable, Optional

import numpy as np
from numpy.typing import NDArray

from .BaseManifold import BaseManifold
from .BranchPoint import BranchPoint
from .Bridge import Bridge
from .FixedPoint import FixedPoint
from .Intersection import Intersection, ManifoldKey
from .IntersectionRegistry import ANCHOR_RESOLUTION_FACTOR, IntersectionRegistry
from .PointStore import Chain
from .Tangle import _Segment, _segment_index
from .budget import GrowthStats

if TYPE_CHECKING:
    from .TangleWorkbench import TangleWorkbench

#: The one save format this tanglepack reads and writes. Bump it on ANY change
#: to what a save holds (a field added, dropped, renamed or re-encoded): a save
#: of another format is refused rather than misread. The law
#: ``test_a_load_restores_every_attribute`` fails when a class the save covers
#: gains an attribute this module does not save.
FORMAT_VERSION = 1
#: How many linked rows :func:`check_map` maps, per direction.
MAP_CHECK_SAMPLE = 8
#: The relative disagreement :func:`check_map` tolerates. A link the other map
#: made (a stable row's image is the point ``map_inv`` was applied to, an
#: unstable row's preimage the point ``map`` was applied to) re-maps only to
#: rounding: at most 9e-16 forward and 5e-15 backward on the law cases
#: (relative to the largest coordinate of the pair, at least 1).
MAP_CHECK_RTOL = 1e-10

STABILITIES = ("unstable", "stable")
_STORE_COLUMNS = ("_xy", "_cdist", "_image", "_preimage", "_chain", "_position")
_FIXED_POINT_ARRAYS = ("coordinates", "unstable_eigenvectors", "stable_eigenvectors", "jacobians")
_TYPE_TAGS = {float: "f", np.float64: "d", int: "i", type(None): "-"}
_NO_KEY = (-1, -1, -1, -1)
_NO_END = (-1, -1, -1)


# ── values ────────────────────────────────────────────────────────────────


def typed_floats(values: Iterable[Optional[float]]) -> tuple[NDArray[np.float64], str]:
    """Optional floats as a float64 array (NaN for None) and one type tag per value."""
    values = list(values)
    tags = "".join(_TYPE_TAGS[type(value)] for value in values)
    return np.array([np.nan if v is None else v for v in values], dtype=np.float64), tags


def floats_of(array: NDArray[np.float64], tags: str) -> list:
    """The values :func:`typed_floats` saved, each of its saved type."""
    cast = {"f": float, "d": np.float64, "i": int, "-": lambda _: None}
    return [cast[tag](value) for value, tag in zip(array.tolist(), tags)]


def _pair(a: int, b: int) -> frozenset[int]:
    """``frozenset((a, b))`` iterating as ``(a, b)``: two ids that collide in the
    set's table iterate in insertion order, and that order is a crossing's
    operand order."""
    pair = frozenset((a, b))
    return pair if tuple(pair) == (a, b) else frozenset((b, a))


def _ragged(groups: Iterable[Iterable[int]]) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    """Groups of ints as their lengths and their concatenation."""
    groups = [np.fromiter(group, dtype=np.int64) for group in groups]
    lengths = np.array([len(group) for group in groups], dtype=np.int64)
    return lengths, np.concatenate(groups) if groups else np.empty(0, dtype=np.int64)


def _unragged(lengths: NDArray[np.int64], flat: NDArray[np.int64]) -> list[list[int]]:
    """The groups :func:`_ragged` saved."""
    stops = np.cumsum(lengths).tolist()
    flat = flat.tolist()
    return [flat[stop - length:stop] for stop, length in zip(stops, lengths.tolist())]


def _rows(rows: list, width: int) -> NDArray[np.int64]:
    """A table of ints with ``width`` columns, empty or not."""
    return np.array(rows, dtype=np.int64).reshape(-1, width)


# ── references ────────────────────────────────────────────────────────────


class Encoder:
    """
    Numbers what a save refers to: fixed points by position, curves and
    crossings in first-seen order, each object once.

    Attributes:
        fixed_points: The workbench's fixed points.
        curves: Every curve the save names, in first-seen order.
        crossings: Every Intersection the save names, in first-seen order.
    """

    def __init__(self, workbench: TangleWorkbench) -> None:
        """
        Args:
            workbench: The workbench being saved.
        """
        self.fixed_points = list(workbench.fixed_points)
        self.curves: list[BaseManifold] = []
        self.crossings: list[Intersection] = []
        self._numbers: dict[int, int] = {}

    def fixed_point(self, fixed_point: FixedPoint) -> int:
        """The fixed point's position on the workbench."""
        return next(i for i, fp in enumerate(self.fixed_points) if fp is fixed_point)

    def key(self, key: Optional[ManifoldKey]) -> tuple[int, int, int, int]:
        """A manifold key as ``(fixed point, stability, orbit, branch)``; all -1 for None."""
        if key is None:
            return _NO_KEY
        fixed_point, stability, orbit_index, branch_index = key
        return self.fixed_point(fixed_point), STABILITIES.index(stability), orbit_index, branch_index

    def end(self, end: int | BranchPoint | None) -> tuple[int, int, int]:
        """A curve end as ``(row, fixed point, orbit index)``, -1 where it is no
        row or no periodic point; all -1 for None."""
        if end is None:
            return _NO_END
        if isinstance(end, BranchPoint):
            return -1, self.fixed_point(end.fixed_point), end.orbit_index
        return int(end), -1, -1

    def segment(self, segment: Optional[tuple]) -> tuple[int, ...]:
        """A crossing's bracket, its two nodes, as two :meth:`end` triples; all -1 for None."""
        return _NO_END * 2 if segment is None else self.end(segment[0]) + self.end(segment[1])

    def curve(self, curve: Optional[BaseManifold]) -> int:
        """The curve's number (-1 for None)."""
        return -1 if curve is None else self._number(curve, self.curves)

    def crossing(self, crossing: Optional[Intersection]) -> int:
        """The crossing's number (-1 for None)."""
        return -1 if crossing is None else self._number(crossing, self.crossings)

    def _number(self, obj: object, table: list) -> int:
        """The object's row in its table, appended on first sight."""
        number = self._numbers.get(id(obj))
        if number is None:
            number = self._numbers[id(obj)] = len(table)
            table.append(obj)
        return number


class Decoder:
    """
    The objects a load has rebuilt, by the numbers :class:`Encoder` gave them.

    Attributes:
        fixed_points: The loaded fixed points.
        store: The loaded point store.
        curves: The loaded curves, by number.
        crossings: The loaded Intersections, by number.
    """

    def __init__(self, fixed_points: list[FixedPoint], workbench: TangleWorkbench) -> None:
        """
        Args:
            fixed_points: The loaded fixed points.
            workbench: The workbench being loaded into.
        """
        self.fixed_points = fixed_points
        self.store = workbench.point_store
        self.curves: list[BaseManifold] = []
        self.crossings: list[Intersection] = []

    def key(self, saved: Iterable[int]) -> Optional[ManifoldKey]:
        """The manifold key :meth:`Encoder.key` saved."""
        fixed_point, stability, orbit_index, branch_index = map(int, saved)
        if fixed_point < 0:
            return None
        return self.fixed_points[fixed_point], STABILITIES[stability], orbit_index, branch_index

    def end(self, row: int, fixed_point: int, orbit_index: int) -> int | BranchPoint | None:
        """The curve end :meth:`Encoder.end` saved."""
        if row >= 0:
            return int(row)
        return None if fixed_point < 0 else self.fixed_points[fixed_point].branch_points[orbit_index]

    def segment(self, saved: list[int]) -> Optional[tuple]:
        """The bracket :meth:`Encoder.segment` saved."""
        return None if saved[0] < 0 and saved[1] < 0 else (self.end(*saved[:3]), self.end(*saved[3:]))

    def curve(self, number: int) -> Optional[BaseManifold]:
        """The curve of a number (None for -1)."""
        return None if number < 0 else self.curves[number]


# ── save ──────────────────────────────────────────────────────────────────


def save_workbench(workbench: TangleWorkbench, path: str | os.PathLike) -> None:
    """
    Write a workbench to ``<path>.npz`` and ``<path>.json``.

    Args:
        workbench: The workbench to save.
        path: The stem of the two files.
    """
    encoder = Encoder(workbench)
    meta, arrays = encode_workbench(workbench, encoder)
    write(path, "workbench", meta, arrays, encoder)


def encode_workbench(workbench: TangleWorkbench, encoder: Encoder) -> tuple[dict, dict]:
    """
    The workbench section of a save, as its JSON part and its arrays.

    The crossing and curve tables are left to :func:`write`, so a session can
    name more crossings first.

    Args:
        workbench: The workbench to save.
        encoder: The save's reference numbering.

    Returns:
        ``(meta, arrays)``.

    Raises:
        CurveReleased: The workbench is lean: a load could not rebuild the
            curve a lean blast freed.
    """
    workbench._refuse_if_lean("save")
    meta: dict = {}
    arrays: dict[str, np.ndarray] = {}
    _save_store(workbench, encoder, meta, arrays)
    _save_fixed_points(encoder, meta, arrays)
    _save_registry(workbench, encoder, meta, arrays)
    _save_bookkeeping(workbench, encoder, meta, arrays)  # numbers the listed curves first
    _save_index(workbench, encoder, meta, arrays)
    return meta, arrays


def write(path: str | os.PathLike, kind: str, meta: dict, arrays: dict, encoder: Encoder) -> None:
    """
    Append the crossing and curve tables and write the two files.

    Args:
        path: The stem of the two files.
        kind: ``"workbench"`` or ``"session"``.
        meta: The JSON part.
        arrays: The arrays.
        encoder: The save's reference numbering.
    """
    _save_crossings(encoder, meta, arrays)  # names more curves: before the curve table
    _save_curves(encoder, meta, arrays)
    try:
        version = metadata.version("tanglepack")
    except metadata.PackageNotFoundError:
        version = None
    path = os.fspath(path)
    np.savez(f"{path}.npz", **arrays)
    with open(f"{path}.json", "w") as handle:
        json.dump({"format": FORMAT_VERSION, "kind": kind, "tanglepack": version, **meta}, handle)


def _save_store(workbench: TangleWorkbench, encoder: Encoder, meta: dict, arrays: dict) -> None:
    """The store's written rows, its capacity, and every chain."""
    store = workbench.point_store
    for name in _STORE_COLUMNS:
        arrays[f"store{name}"] = getattr(store, name)[: len(store)]
    arrays["chain_lengths"], arrays["chain_rows"] = _ragged(chain.ids for chain in store.chains)
    arrays["chain_stretch"], stretch_types = typed_floats(chain.stretch_param for chain in store.chains)
    arrays["chain_anchor"] = _rows([encoder.end(chain.anchor)[1:] for chain in store.chains], 2)
    meta["store"] = {"capacity": len(store._cdist), "chain_stretch_types": stretch_types}


def _save_fixed_points(encoder: Encoder, meta: dict, arrays: dict) -> None:
    """Each fixed point's eigendata and its slot table, in insertion order."""
    meta["fixed_points"] = []
    for i, fixed_point in enumerate(encoder.fixed_points):
        for name in _FIXED_POINT_ARRAYS:
            arrays[f"fixed_point{i}_{name}"] = np.array(getattr(fixed_point, name), dtype=np.float64)
        arrays[f"fixed_point{i}_partial_jacobians"] = np.array(
            [np.full((2, 2), np.nan) if m is None else m for m in fixed_point.partial_jacobians],
            dtype=np.float64,
        )
        arrays[f"fixed_point{i}_scalars"], types = typed_floats(
            [*fixed_point.unstable_eigenvalues, *fixed_point.stable_eigenvalues, fixed_point.accuracy]
        )
        meta["fixed_points"].append({
            "period": fixed_point.period,
            "k_value": fixed_point.k_value,
            "label": fixed_point.label,
            "ndarray_fields": [name for name in _FIXED_POINT_ARRAYS
                               if isinstance(getattr(fixed_point, name), np.ndarray)],
            "scalar_types": types,
            "unsolved_partial_jacobians": [m is None for m in fixed_point.partial_jacobians],
            "slots": [
                [STABILITIES.index(stability), orbit_index, branch, chain.index]
                for (stability, orbit_index, branch), chain in fixed_point._slots.items()
            ],
        })


def _save_registry(workbench: TangleWorkbench, encoder: Encoder, meta: dict, arrays: dict) -> None:
    """The registry's crossings (first in the crossing table), orders, links and counters."""
    registry = workbench.intersection_registry
    for crossing_id, crossing in registry._store.items():
        assert crossing.id == crossing_id, "a registry crossing carries its own id"
        encoder.crossing(crossing)
    table = registry.iterate_table
    arrays["registry_unstable_order"] = np.array(registry._unstable_order, dtype=np.int64)
    arrays["registry_stable_order"] = np.array(registry._stable_order, dtype=np.int64)
    for name, links in (("forward", table._forward), ("backward", table._backward)):
        arrays[f"registry_{name}"] = _rows(
            [(source, n, target) for source, by_n in links.items() for n, target in by_n.items()], 3)
    arrays["registry_cdist_tol"], cdist_tol_type = typed_floats([registry.cdist_tol])
    meta["registry"] = {
        "cdist_tol_type": cdist_tol_type,
        "size": len(registry._store),
        "next_id": registry._next_id,
        "generation": registry._generation,
        "orders_version": registry._orders_version,
        "unresolved_links": [[source, list(images)] for source, images in registry.unresolved_links.items()],
    }


def _save_index(workbench: TangleWorkbench, encoder: Encoder, meta: dict, arrays: dict) -> None:
    """The segment index, exactly: records, edge maps, pairs (each in its operand order, sorted), owners."""
    tangle = workbench.Tangle
    arrays["index_segments"] = _rows([
        (segment.id, encoder.curve(segment.manifold), *encoder.end(segment.row0), segment.row1,
         segment.in_rtree)
        for segment in tangle._seg_lookup.values()
    ], 7)
    arrays["index_edges"] = _rows(list(tangle._edge_to_sid.items()), 2)
    arrays["index_query_edges"] = np.sort(np.fromiter(tangle._query_edges, dtype=np.int64))
    # each pair keeps its iteration order, its operand order; the sets are sorted
    for name, pairs in (("crossing", tangle._intersecting_segments), ("processed", tangle._processed_pairs)):
        arrays[f"index_{name}_pairs"] = _rows(sorted(tuple(pair) for pair in pairs), 2)
    arrays["index_owner_curves"] = np.array(
        [encoder.curve(curve) for curve in tangle._manifold_segs], dtype=np.int64)
    arrays["index_owned_counts"], arrays["index_owned"] = _ragged(
        sorted(sids) for sids in tangle._manifold_segs.values())
    arrays["index_claimed_sids"] = np.array(list(tangle._seg_manifolds), dtype=np.int64)
    arrays["index_claim_counts"], arrays["index_claims"] = _ragged(
        sorted(map(encoder.curve, owners)) for owners in tangle._seg_manifolds.values())
    meta["index"] = {"next_sid": tangle._next_sid, "has_store": tangle._store is not None}


def _save_bookkeeping(workbench: TangleWorkbench, encoder: Encoder, meta: dict, arrays: dict) -> None:
    """The workbench's curve lists, growth stats, reported findings, counters and settings."""
    assert all(bridge.id == bridge_id for bridge_id, bridge in workbench._bridges.items())
    machine = workbench._man_machine
    growth = list(workbench._growth_stats.items())
    artifacts = sorted(workbench._reported_artifacts)
    state = workbench._generation_state
    meta["workbench"] = {
        "manifolds": [[encoder.key(key), encoder.curve(m)] for key, m in workbench._manifolds.items()],
        "bridges": [encoder.curve(bridge) for bridge in workbench._bridges.values()],
        "partial_bridges": [encoder.curve(bridge) for bridge in workbench._partial_bridges],
        "image_curves": [encoder.curve(curve) for curve in workbench._image_curves],
        "half_iterated_bridges": sorted(map(encoder.curve, workbench._half_iterated_bridges)),
        "bridges_at": [[crossing_id, ids] for crossing_id, ids in workbench._bridges_at.items()],
        "growth_stats": [
            [encoder.fixed_point(fp), stability, branch, stats.points_created]
            for (fp, stability, branch), stats in growth
        ],
        "reported_flags": [workbench._reported_flags[0], sorted(workbench._reported_flags[1])],
        "reported_artifacts": [list(artifact[:3]) for artifact in artifacts],
        "registry_id_epoch": workbench._registry_id_epoch,
        "bridge_id_epoch": workbench._bridge_id_epoch,
        "mutations": workbench._mutations,
        "generation": workbench._generation,
        "generation_state": None if state is None else [state[0], state[1], list(state[2])],
        "points_created": machine.points_created,
        "max_points_per_step": machine.max_points_per_step,
    }
    arrays["growth_factors"], meta["workbench"]["growth_factor_types"] = typed_floats(
        stats.growth_factor for _, stats in growth)
    arrays["growth_max_abs_coords"], meta["workbench"]["growth_max_abs_coord_types"] = typed_floats(
        stats.max_abs_coord for _, stats in growth)
    arrays["reported_artifact_points"] = np.array(
        [artifact[3] for artifact in artifacts], dtype=np.float64).reshape(-1, 2)
    arrays["settings"], meta["workbench"]["settings_types"] = typed_floats([
        machine.area_cutoff, workbench.escape_radius, workbench.resolution_ratio_threshold,
        machine.max_abs_coord, workbench.anchor_resolution_factor, machine.near_orbit_radius,
    ])


def _save_crossings(encoder: Encoder, meta: dict, arrays: dict) -> None:
    """One row per Intersection: ints (brackets as ends), typed floats, labels."""
    crossings = encoder.crossings
    arrays["crossing_ints"] = _rows([
        (
            -1 if x.id is None else x.id,
            x.crossing_sign,
            *(tuple(x.seg_ids) if x.seg_ids is not None else (-1, -1)),
            encoder.curve(x.unstable_manifold),
            *encoder.segment(x.unstable_segment),
            *encoder.segment(x.stable_segment),
            x.store is not None,
            *encoder.key(x.manifold_a_key),
            *encoder.key(x.manifold_b_key),
        )
        for x in crossings
    ], 26)
    columns = {
        "x": [None if x.coords is None else x.coords[0] for x in crossings],
        "y": [None if x.coords is None else x.coords[1] for x in crossings],
        "unstable_cdist": [x.unstable_cdist for x in crossings],
        "stable_cdist": [x.stable_cdist for x in crossings],
        "resolution_ratio": [x.resolution_ratio for x in crossings],
        "chord_deviation": [x.chord_deviation for x in crossings],
    }
    meta["crossings"] = {"labels": [x.label for x in crossings], "types": {}}
    for name, values in columns.items():
        arrays[f"crossing_{name}"], meta["crossings"]["types"][name] = typed_floats(values)


def _save_curves(encoder: Encoder, meta: dict, arrays: dict) -> None:
    """One row per curve: ends, key, version and a bridge's flags; names in the JSON."""
    curves = encoder.curves
    assert all(type(curve) in (BaseManifold, Bridge) for curve in curves)
    arrays["curve_ints"] = _rows([
        (
            isinstance(curve, Bridge),
            *encoder.end(curve._root),
            *encoder.end(curve._tail),
            STABILITIES.index(curve.stability),
            encoder.fixed_point(curve.fixed_point),
            -1 if curve.branch_index is None else curve.branch_index,
            *encoder.key(curve.manifold_key),
            curve._version,
            getattr(curve, "iterated", False),
            *(-1 if crossing_id is None else crossing_id for crossing_id in (
                getattr(curve, "first_intersection", None), getattr(curve, "second_intersection", None))),
        )
        for curve in curves
    ], 18)
    arrays["curve_stretch"], stretch_types = typed_floats(curve.stretch_param for curve in curves)
    arrays["curve_dropped_at"], dropped_types = typed_floats(
        getattr(curve, "dropped_at", None) for curve in curves)
    meta["curves"] = {
        "names": [curve.name for curve in curves],
        "stretch_types": stretch_types,
        "dropped_at_types": dropped_types,
    }
    # Only a width blast sets widths: a save without one holds no width column,
    # byte-equal to a save before the width guard existed.
    widths = [getattr(curve, "width", None) for curve in curves]
    if any(width is not None for width in widths):
        arrays["curve_width"], meta["curves"]["width_types"] = typed_floats(widths)
        arrays["curve_width_saddle"] = np.array([
            -1 if (saddle := getattr(curve, "width_saddle", None)) is None
            else encoder.fixed_point(saddle)
            for curve in curves
        ], dtype=np.int64)


# ── load ──────────────────────────────────────────────────────────────────


def read(path: str | os.PathLike) -> tuple[dict, dict[str, np.ndarray]]:
    """
    The JSON part and the arrays of a save, refusing another format.

    Args:
        path: The stem of the two files.

    Returns:
        ``(meta, arrays)``.

    Raises:
        ValueError: The save is of another format, or the npz holds an object
            array (``allow_pickle=False``).
    """
    path = os.fspath(path)
    with open(f"{path}.json") as handle:
        meta = json.load(handle)
    if meta.get("format") != FORMAT_VERSION:
        raise ValueError(f"{path}.json is save format {meta.get('format')!r}; "
                         f"this tanglepack reads format {FORMAT_VERSION}")
    with np.load(f"{path}.npz", allow_pickle=False) as npz:
        arrays = {name: npz[name] for name in npz.files}
    return meta, arrays


def restore_workbench(workbench: TangleWorkbench, meta: dict, arrays: dict) -> Decoder:
    """
    Put a saved state into a freshly built workbench.

    Args:
        workbench: A workbench just built around the save's map.
        meta: The save's JSON part.
        arrays: The save's arrays.

    Returns:
        The rebuilt objects, for a session section to resolve its references.

    Raises:
        ValueError: The workbench's map disagrees with the saved points
            (:func:`check_map`).
    """
    assert not workbench.fixed_points and not len(workbench.point_store), "load into a fresh workbench"
    decoder = Decoder(_load_fixed_points(meta, arrays), workbench)
    workbench.fixed_points = decoder.fixed_points
    _load_store(decoder, meta, arrays)
    _load_curves(decoder, meta, arrays)
    _load_crossings(decoder, meta, arrays)
    _load_registry(workbench, decoder, meta, arrays)
    _load_index(workbench, decoder, meta, arrays)
    _load_bookkeeping(workbench, decoder, meta, arrays)
    check_map(workbench)
    return decoder


def check_map(workbench: TangleWorkbench) -> None:
    """
    Refuse a workbench whose map is not the one its points were made with.

    Maps up to :data:`MAP_CHECK_SAMPLE` rows with an image link, evenly spaced
    among them, through the workbench's batched map, and as many rows with a
    preimage link through the batched inverse, and compares each result with
    the saved xy of the linked row. A store with no linked row is not checked.

    Args:
        workbench: A workbench whose state is restored.

    Raises:
        ValueError: A mapped row disagrees with its saved image (preimage)
            beyond :data:`MAP_CHECK_RTOL`, relative to the larger coordinate of
            the two points and at least 1; the message names the largest
            disagreement.
    """
    store, system = workbench.point_store, workbench.dynamical_system
    n = len(store)
    for name, links, apply in (("map", store._image[:n], system.map_batch),
                               ("map_inv", store._preimage[:n], system.map_inv_batch)):
        linked = np.flatnonzero(links >= 0)
        if not len(linked):
            continue
        spaced = np.linspace(0, len(linked) - 1, min(MAP_CHECK_SAMPLE, len(linked))).astype(int)
        rows = np.unique(linked[spaced])
        source, saved = store._xy[rows], store._xy[links[rows]]
        mapped = apply(source)
        scale = np.maximum(1.0, np.maximum(np.abs(source).max(axis=1), np.abs(saved).max(axis=1)))
        error = np.abs(mapped - saved).max(axis=1) / scale
        worst = int(np.argmax(error))
        if not error[worst] <= MAP_CHECK_RTOL:  # NaN fails too
            raise ValueError(
                f"the {name} given to load differs from the one the save was made with: row "
                f"{rows[worst]} at {source[worst].tolist()} maps to {mapped[worst].tolist()}, its saved "
                f"{'image' if name == 'map' else 'preimage'} is {saved[worst].tolist()} "
                f"(relative error {error[worst]:.3g} > {MAP_CHECK_RTOL:g})")


def _load_fixed_points(meta: dict, arrays: dict) -> list[FixedPoint]:
    """The fixed points, their slot tables still empty."""
    fixed_points = []
    for i, saved in enumerate(meta["fixed_points"]):
        period = saved["period"]
        fixed_point = FixedPoint(period)
        for name in _FIXED_POINT_ARRAYS:
            saved_array = arrays[f"fixed_point{i}_{name}"]
            setattr(fixed_point, name, saved_array.copy() if name in saved["ndarray_fields"]
                    else [row.copy() for row in saved_array])
        partial_jacobians = arrays[f"fixed_point{i}_partial_jacobians"]
        fixed_point.partial_jacobians = [
            None if unsolved else row.copy()
            for row, unsolved in zip(partial_jacobians, saved["unsolved_partial_jacobians"])
        ]
        scalars = floats_of(arrays[f"fixed_point{i}_scalars"], saved["scalar_types"])
        fixed_point.unstable_eigenvalues = scalars[:period]
        fixed_point.stable_eigenvalues = scalars[period:2 * period]
        fixed_point.accuracy = scalars[-1]
        fixed_point.k_value = saved["k_value"]
        fixed_point.label = saved["label"]
        fixed_points.append(fixed_point)
    return fixed_points


def _load_store(decoder: Decoder, meta: dict, arrays: dict) -> None:
    """The store's columns at their saved capacity, every chain, and the slots."""
    store = decoder.store
    rows = len(arrays["store_cdist"])
    for name in _STORE_COLUMNS:
        saved, fresh = arrays[f"store{name}"], getattr(store, name)
        column = np.empty((meta["store"]["capacity"],) + fresh.shape[1:], dtype=fresh.dtype)
        column[:rows] = saved
        setattr(store, name, column)
    store._n = rows
    store._alias_xy()
    stretches = floats_of(arrays["chain_stretch"], meta["store"]["chain_stretch_types"])
    chain_rows = _unragged(arrays["chain_lengths"], arrays["chain_rows"])
    anchors = arrays["chain_anchor"].tolist()
    for index, (ids, stretch, anchor) in enumerate(zip(chain_rows, stretches, anchors)):
        chain = Chain(store, index, np.empty(0, dtype=np.int64), stretch)
        chain.ids = np.array(ids, dtype=np.int64)  # the columns already name it
        chain.anchor = decoder.end(-1, *anchor)
        store.chains.append(chain)
    for fixed_point, saved in zip(decoder.fixed_points, meta["fixed_points"]):
        for stability, orbit_index, branch, chain in saved["slots"]:
            fixed_point.attach_chain(store.chains[chain], STABILITIES[stability], orbit_index, branch)


def _load_curves(decoder: Decoder, meta: dict, arrays: dict) -> None:
    """Every curve, built by its constructor, then given its saved version and flags."""
    saved = meta["curves"]
    stretches = floats_of(arrays["curve_stretch"], saved["stretch_types"])
    dropped = floats_of(arrays["curve_dropped_at"], saved["dropped_at_types"])
    rows = arrays["curve_ints"].tolist()
    widths = (
        floats_of(arrays["curve_width"], saved["width_types"])
        if "curve_width" in arrays else [None] * len(rows)
    )
    saddles = [
        None if i < 0 else decoder.fixed_points[i]
        for i in arrays.get("curve_width_saddle", np.full(len(rows), -1)).tolist()
    ]
    for row, name, stretch, dropped_at, width, saddle in zip(
        rows, saved["names"], stretches, dropped, widths, saddles
    ):
        is_bridge, root, tail = row[0], decoder.end(*row[1:4]), decoder.end(*row[4:7])
        has_row = not (isinstance(root, BranchPoint) and isinstance(tail, BranchPoint))
        fields = dict(
            root=root, tail=tail, stability=STABILITIES[row[7]], stretch_param=stretch,
            fixed_point=decoder.fixed_points[row[8]], name=name,
            branch_index=None if row[9] < 0 else row[9], manifold_key=decoder.key(row[10:14]),
            store=decoder.store if has_row else None,
        )
        if is_bridge:
            first, second = (None if i < 0 else i for i in row[16:18])
            curve = Bridge(**fields, first_intersection=first, second_intersection=second)
            curve.iterated, curve.dropped_at = bool(row[15]), dropped_at
            curve.width, curve.width_saddle = width, saddle
        else:
            curve = BaseManifold(**fields)
        curve._version = row[14]
        decoder.curves.append(curve)


def _load_crossings(decoder: Decoder, meta: dict, arrays: dict) -> None:
    """Every Intersection, with its saved value types and object references."""
    saved = meta["crossings"]
    columns = {name: floats_of(arrays[f"crossing_{name}"], tags) for name, tags in saved["types"].items()}
    for i, (row, label) in enumerate(zip(arrays["crossing_ints"].tolist(), saved["labels"])):
        x, y = columns["x"][i], columns["y"][i]
        decoder.crossings.append(Intersection(
            coords=None if x is None else (x, y),
            unstable_cdist=columns["unstable_cdist"][i],
            stable_cdist=columns["stable_cdist"][i],
            seg_ids=None if row[2] < 0 else _pair(row[2], row[3]),
            id=None if row[0] < 0 else row[0],
            label=label,
            manifold_a_key=decoder.key(row[18:22]),
            manifold_b_key=decoder.key(row[22:26]),
            unstable_manifold=decoder.curve(row[4]),
            unstable_segment=decoder.segment(row[5:11]),
            stable_segment=decoder.segment(row[11:17]),
            store=decoder.store if row[17] else None,
            crossing_sign=row[1],
            resolution_ratio=columns["resolution_ratio"][i],
            chord_deviation=columns["chord_deviation"][i],
        ))


def _load_registry(workbench: TangleWorkbench, decoder: Decoder, meta: dict, arrays: dict) -> None:
    """The registry: its crossings in their order, orders, iterate table, links, counters."""
    saved = meta["registry"]
    (cdist_tol,) = floats_of(arrays["registry_cdist_tol"], saved["cdist_tol_type"])
    registry = IntersectionRegistry(cdist_tol)
    registry._store = {crossing.id: crossing for crossing in decoder.crossings[: saved["size"]]}
    registry._next_id = saved["next_id"]
    registry._generation = saved["generation"]
    registry._orders_version = saved["orders_version"]
    registry._unstable_order = arrays["registry_unstable_order"].tolist()
    registry._stable_order = arrays["registry_stable_order"].tolist()
    # the sorted keys are the crossings' own cdists, in the saved order
    registry._unstable_keys = [registry._store[i].unstable_cdist for i in registry._unstable_order]
    registry._stable_keys = [registry._store[i].stable_cdist for i in registry._stable_order]
    table = registry.iterate_table
    for name, links in (("forward", table._forward), ("backward", table._backward)):
        for source, n, target in arrays[f"registry_{name}"].tolist():
            links.setdefault(source, {})[n] = target
    registry.unresolved_links = {source: tuple(images) for source, images in saved["unresolved_links"]}
    workbench._intersection_registry = registry


def _load_index(workbench: TangleWorkbench, decoder: Decoder, meta: dict, arrays: dict) -> None:
    """The segment index as saved, its rtrees bulk-loaded from the records they held."""
    tangle = workbench.Tangle
    for sid, curve, row, fixed_point, orbit_index, row1, in_rtree in arrays["index_segments"].tolist():
        tangle._seg_lookup[sid] = _Segment(
            sid, decoder.curves[curve], decoder.store, decoder.end(row, fixed_point, orbit_index), row1,
            bool(in_rtree))
    tangle._edge_to_sid = dict(arrays["index_edges"].tolist())
    tangle._query_edges = set(arrays["index_query_edges"].tolist())
    tangle._intersecting_segments = {_pair(a, b) for a, b in arrays["index_crossing_pairs"].tolist()}
    tangle._processed_pairs = {_pair(a, b) for a, b in arrays["index_processed_pairs"].tolist()}
    owned = _unragged(arrays["index_owned_counts"], arrays["index_owned"])
    tangle._manifold_segs = defaultdict(set, {
        decoder.curves[curve]: set(sids)
        for curve, sids in zip(arrays["index_owner_curves"].tolist(), owned)
    })
    claims = _unragged(arrays["index_claim_counts"], arrays["index_claims"])
    tangle._seg_manifolds = defaultdict(set, {
        sid: {decoder.curves[curve] for curve in curves}
        for sid, curves in zip(arrays["index_claimed_sids"].tolist(), claims)
    })
    tangle._next_sid = meta["index"]["next_sid"]
    tangle._store = decoder.store if meta["index"]["has_store"] else None
    for stability in STABILITIES:
        tangle._rtrees[stability] = _segment_index([
            (segment.id, segment.bounds, None)
            for segment in tangle._seg_lookup.values()
            if segment.in_rtree and segment.manifold.stability == stability
        ])


def _load_bookkeeping(workbench: TangleWorkbench, decoder: Decoder, meta: dict, arrays: dict) -> None:
    """The curve lists, growth stats, reported findings, counters and settings."""
    saved = meta["workbench"]
    curve = decoder.curve
    workbench._manifolds = {decoder.key(key): curve(number) for key, number in saved["manifolds"]}
    workbench._bridges = {bridge.id: bridge for bridge in map(curve, saved["bridges"])}
    workbench._partial_bridges = [curve(number) for number in saved["partial_bridges"]]
    workbench._image_curves = [curve(number) for number in saved["image_curves"]]
    workbench._half_iterated_bridges = {curve(number) for number in saved["half_iterated_bridges"]}
    workbench._bridges_at = {
        crossing_id: [tuple(bridge_id) for bridge_id in ids] for crossing_id, ids in saved["bridges_at"]}

    factors = floats_of(arrays["growth_factors"], saved["growth_factor_types"])
    coords = floats_of(arrays["growth_max_abs_coords"], saved["growth_max_abs_coord_types"])
    workbench._growth_stats = {
        (decoder.fixed_points[fp], stability, branch): GrowthStats(points, factor, coord)
        for (fp, stability, branch, points), factor, coord in zip(saved["growth_stats"], factors, coords)
    }
    epoch, flagged = saved["reported_flags"]
    workbench._reported_flags = (epoch, frozenset(flagged))
    workbench._reported_artifacts = {
        (stability, curve_a, curve_b, (x, y))
        for (stability, curve_a, curve_b), (x, y) in zip(
            saved["reported_artifacts"], arrays["reported_artifact_points"].tolist())
    }
    workbench._registry_id_epoch = saved["registry_id_epoch"]
    workbench._bridge_id_epoch = saved["bridge_id_epoch"]
    workbench._mutations = saved["mutations"]
    workbench._generation = saved["generation"]
    state = saved["generation_state"]
    workbench._generation_state = None if state is None else (state[0], state[1], tuple(state[2]))

    machine = workbench._man_machine
    # The last two settings are optional (2026-10-29): an earlier format-1
    # file has four and loads the default guard and no graded mesh.
    area_cutoff, escape_radius, threshold, max_abs_coord, anchor_factor, near_orbit_radius = (
        floats_of(arrays["settings"], saved["settings_types"]) + [ANCHOR_RESOLUTION_FACTOR, None]
    )[:6]
    machine.area_cutoff, machine.max_points_per_step = area_cutoff, saved["max_points_per_step"]
    machine.points_created, machine.max_abs_coord = saved["points_created"], max_abs_coord
    workbench.escape_radius = escape_radius
    workbench.resolution_ratio_threshold = threshold
    workbench.anchor_resolution_factor = anchor_factor
    machine.near_orbit_radius = near_orbit_radius
