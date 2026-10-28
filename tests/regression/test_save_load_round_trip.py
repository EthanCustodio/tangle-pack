"""
Regression (N10a): a saved workbench or session loads as it was, and goes on as it would.

Every law case is built once for this module, saved, loaded and re-saved: the
two saves are byte-equal (every array with its dtype, every JSON key), every
crossing keeps the type of each field and names a curve the loaded workbench
holds, and every slot holds the chain its branch runs on. A loaded session
remembers the pips, the letters and the zones, so its first ``topology()``
reads as the original's. A continuation grows (k10, nested) or blasts (k28
two blasts, p3, heteroclinic) a fresh build and its loaded copy and recomputes
both: the crossings, the store, the bridges and the topology come out equal.
Segment ids are drawn per index (``Tangle._next_sid``), so a load in a process
that has built other cases is as good as one in a fresh process; one test
loads in a fresh interpreter all the same, against hidden process-global
state.

Every attribute of every object the save covers is either saved or listed as
derived, and a loaded object has the original's saved attributes with the same
value types: a class that gains an attribute fails here until a save holds it
or the attribute is listed as derived.

No value is pinned: the original is compared with its loaded copy.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from cases import BUILDERS, HENON_K28, LAW_CASES, ZONE_MIN_SEPARATION
from tanglepack import Intersection, TangleSession, TangleWorkbench
from tanglepack.loom.BridgeAlphabet import BridgeAlphabet
from tanglepack.loom.ResonanceZone import ResonanceZone
from tanglepack.numerics.BaseManifold import BaseManifold
from tanglepack.numerics.Bridge import Bridge
from tanglepack.numerics.FixedPoint import FixedPoint
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry
from tanglepack.numerics.IterateTable import IterateTable
from tanglepack.numerics.ManifoldMachine import ManifoldMachine
from tanglepack.numerics.PointStore import Chain, PointStore
from tanglepack.numerics.Tangle import Tangle


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


def test_a_loaded_index_keeps_each_crossing_pair_in_its_operand_order(saved_workbench):
    # A re-save cannot see this: a pair saved and rebuilt in one wrong order re-saves alike.
    workbench, loaded, _ = saved_workbench
    for pairs in ("_intersecting_segments", "_processed_pairs"):
        assert {tuple(pair) for pair in getattr(loaded.Tangle, pairs)} == {
            tuple(pair) for pair in getattr(workbench.Tangle, pairs)}


def test_a_loaded_fixed_point_keeps_its_array_containers(saved_workbench):
    # The solver leaves ``coordinates`` a (period, 2) ndarray; the other fields are lists.
    workbench, loaded, _ = saved_workbench
    for original, copy in zip(workbench.fixed_points, loaded.fixed_points):
        for name in ("coordinates", "unstable_eigenvectors", "stable_eigenvectors", "jacobians"):
            before, after = getattr(original, name), getattr(copy, name)
            assert type(after) is type(before), name
            assert np.shape(after) == np.shape(before), name


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


# ── sessions ──────────────────────────────────────────────────────────────


def _load_session(session, path):
    system = session.workbench.dynamical_system
    return TangleSession.load(path, system.map, system.map_inv, system.jacobian)


def _topology(session) -> tuple:
    """What a recompute of the topology reads: pips, classes, words and the refined matrix."""
    topology = session.topology()
    labels, matrix = topology.dynamics.transition_matrix(refined=True)
    return (topology.strong_pips, topology.classes.describe(), topology.dynamics.describe(),
            labels, matrix.tolist())


@pytest.fixture(scope="module", params=LAW_CASES)
def saved_session(request, tmp_path_factory):
    """``(case, loaded session, save path)`` of one law case."""
    case = BUILDERS[request.param]()
    case.session.topology().classes  # letters the classes: the alphabet is history
    path = tmp_path_factory.mktemp(request.param) / "session"
    case.session.save(path)
    return case, _load_session(case.session, path), path


def test_a_loaded_session_resaves_byte_equal(saved_session):
    _, loaded, path = saved_session
    loaded.save(path.with_name("resaved"))
    _assert_same_files(path, path.with_name("resaved"))


def test_a_loaded_session_recomputes_the_original_topology(saved_session):
    case, loaded, _ = saved_session
    assert _topology(loaded) == _topology(case.session)


def test_a_session_save_loads_as_its_workbench(saved_session):
    case, _, path = saved_session
    assert _numbers(_load_workbench(case.workbench, path)) == _numbers(case.workbench)


def test_a_workbench_save_is_no_session(tmp_path):
    workbench = BUILDERS["k10"]().workbench
    workbench.save(tmp_path / "k10")
    with pytest.raises(ValueError, match="TangleWorkbench.load"):
        _load_session(TangleSession(workbench.dynamical_system.map, workbench.dynamical_system.map_inv),
                      tmp_path / "k10")


def _blast_and_recompute(session) -> tuple:
    zone = max(session.resonance_zones.values(), key=lambda zone: zone.area)
    result = session.blast_zone(zone, num_iterations=1, fixed_point=[zone.fixed_point],
                                min_separation=ZONE_MIN_SEPARATION, diagnose=True)
    fixed_points = session.workbench.fixed_points
    outcome = (result.completed_iterations, result.points_created, result.too_close,
               {fixed_points.index(fp): depth for fp, depth in result.pip_depth.items()})
    return outcome, _numbers(session.workbench), _topology(session)


@pytest.mark.parametrize("name", ["k28_two_blasts", "p3", "heteroclinic"])
def test_a_loaded_session_blasts_as_the_original(name, tmp_path):
    session = BUILDERS[name]().session
    session.topology().classes
    session.save(tmp_path / name)
    loaded = _load_session(session, tmp_path / name)
    assert _blast_and_recompute(loaded) == _blast_and_recompute(session)


def test_a_loaded_session_remembers_its_dropped_children(tmp_path):
    case = BUILDERS["nested_unblasted"]()
    session, outer = case.session, case.fixed_points[0]
    session.forget_strong_pip(outer)  # no child is exempt from the proximity guard
    zone = max(session.resonance_zones.values(), key=lambda zone: zone.area)
    session.blast_zone(zone, 4, fixed_point=[outer], min_separation=0.02, diagnose=False)
    dropped = {b.id: b.dropped_at for b in session.workbench.bridges if b.dropped_at is not None}
    assert dropped, "the blast must drop a child for this test to mean anything"

    session.save(tmp_path / "session")
    loaded = _load_session(session, tmp_path / "session")
    assert {b.id: b.dropped_at for b in loaded.workbench.bridges if b.dropped_at is not None} == dropped
    for copy in (session, loaded):
        zone = max(copy.resonance_zones.values(), key=lambda zone: zone.area)
        copy.blast_zone(zone, 1, fixed_point=[copy.workbench.fixed_points[0]], min_separation=0.02,
                        diagnose=False)
    assert _numbers(loaded.workbench) == _numbers(session.workbench)


def test_a_zone_polygon_loads_as_frozen(tmp_path):
    # The polygon is history: blast_zone tests membership against the one captured
    # when the zone was defined. The reversed ring bounds the same region, and is
    # not what a capture from the loaded workbench would give.
    session = BUILDERS["k28_two_blasts"]().session
    zone = next(iter(session.resonance_zones.values()))
    zone.boundary_vertices = zone.boundary_vertices[::-1].copy()
    session.save(tmp_path / "session")
    (loaded,) = _load_session(session, tmp_path / "session").resonance_zones.values()
    assert loaded.boundary_vertices.tobytes() == zone.boundary_vertices.tobytes()


_FRESH_LOAD = """
import sys
from tanglepack import TangleSession
from tanglepack.examples import henon_jacobian, henon_map, henon_map_inverse
k, b = map(float, sys.argv[3:])
session = TangleSession.load(sys.argv[1], henon_map(k, b), henon_map_inverse(k, b), henon_jacobian(k, b))
session.save(sys.argv[2])
"""


def test_a_fresh_interpreter_loads_the_same_session(tmp_path):
    session = BUILDERS["k28_two_blasts"]().session
    session.save(tmp_path / "saved")
    source = Path(__file__).resolve().parents[2] / "src"
    arguments = [str(tmp_path / "saved"), str(tmp_path / "resaved"), *map(repr, HENON_K28)]
    subprocess.run([sys.executable, "-c", _FRESH_LOAD, *arguments],
                   check=True, env={"PYTHONPATH": str(source)}, timeout=120)
    _assert_same_files(tmp_path / "saved", tmp_path / "resaved")


# ── what a save covers ────────────────────────────────────────────────────

# Every attribute of a class the save covers is in exactly one of these two
# tables. A new attribute is in neither and fails the law below: save it (and
# bump ``persist.FORMAT_VERSION``) or, if it is derived, list it as not saved.
_SAVED: dict[type, frozenset[str]] = {
    TangleSession: frozenset({"workbench", "_strong_pips", "bridge_alphabet", "resonance_zones", "_blast_stats"}),
    TangleWorkbench: frozenset({
        "point_store", "fixed_points", "Tangle", "_man_machine", "_intersection_registry", "_manifolds",
        "_bridges", "_partial_bridges", "_image_curves", "_half_iterated_bridges", "_bridges_at",
        "_growth_stats", "_reported_flags", "_reported_artifacts", "_registry_id_epoch", "_bridge_id_epoch",
        "_mutations", "_generation", "_generation_state", "escape_radius", "resolution_ratio_threshold",
    }),
    PointStore: frozenset({"_xy", "_cdist", "_image", "_preimage", "_chain", "_position", "_n", "chains"}),
    Chain: frozenset({"store", "index", "ids", "stretch_param", "anchor"}),
    FixedPoint: frozenset({
        "period", "branch_points", "_slots", "coordinates", "unstable_eigenvectors", "unstable_eigenvalues",
        "stable_eigenvectors", "stable_eigenvalues", "accuracy", "jacobians", "partial_jacobians", "k_value",
        "label",
    }),
    BaseManifold: frozenset({
        "_root", "_tail", "_store", "_version", "stability", "stretch_param", "fixed_point", "name",
        "branch_index", "manifold_key",
    }),
    IntersectionRegistry: frozenset({
        "_store", "_next_id", "cdist_tol", "iterate_table", "unresolved_links", "_generation",
        "_orders_version", "_unstable_order", "_stable_order",
        "_unstable_keys", "_stable_keys",  # the saved crossings' cdists in the saved orders
    }),
    IterateTable: frozenset({"_forward", "_backward"}),
    Intersection: frozenset({
        "id", "label", "coords", "unstable_cdist", "stable_cdist", "seg_ids", "manifold_a_key",
        "manifold_b_key", "unstable_manifold", "unstable_segment", "stable_segment", "store",
        "crossing_sign", "resolution_ratio", "chord_deviation",
    }),
    Tangle: frozenset({
        "_next_sid", "_edge_to_sid", "_store", "_seg_lookup", "_manifold_segs", "_seg_manifolds",
        "_query_edges", "_intersecting_segments", "_processed_pairs",
    }),
    ManifoldMachine: frozenset({"area_cutoff", "max_points_per_step", "points_created", "max_abs_coord"}),
    ResonanceZone: frozenset({
        "fixed_point", "stable_branch_key", "boundary_intersection", "cut_intersections",
        "intersection_ids", "previous_tails", "boundary_arcs", "boundary_vertices",
    }),
    BridgeAlphabet: frozenset({"assigned"}),
}
_SAVED[Bridge] = _SAVED[BaseManifold] | {"first_intersection", "second_intersection", "iterated", "dropped_at"}

# Derived or lazy (``persist``'s "What is NOT saved"), or built by the load's
# constructors around the map it is given.
_NOT_SAVED: dict[type, frozenset[str]] = {
    TangleSession: frozenset({"_trellises", "_topologies"}),  # caches; topology() recomputes them
    TangleWorkbench: frozenset({"dynamical_system", "_fp_solver", "_man_maker", "_bridge_iterator"}),  # the map's
    PointStore: frozenset({"_xy_read"}),  # the read-only alias of _xy
    FixedPoint: frozenset({"_branch_cycle_cache", "_branch_position_cache"}),  # memos
    IntersectionRegistry: frozenset({  # caches, each behind its version counter
        "_rank_maps", "_rank_version", "_stable_branch_index", "_stable_branch_version", "_runs_cache",
        "_runs_version", "_resolved_cache", "_resolved_version",
    }),
    Intersection: frozenset({"brackets"}),  # lazy
    Tangle: frozenset({"_rtrees"}),  # bulk-loaded from the saved records
    ManifoldMachine: frozenset({"system", "_call_depth"}),  # the map's; zero between calls
}


def _attributes(obj) -> dict[str, type]:
    """An object's attribute names (its set slots for a slotted class) and their value types."""
    names = vars(obj) if hasattr(obj, "__dict__") else [n for n in obj.__slots__ if hasattr(obj, n)]
    return {name: type(getattr(obj, name)) for name in names}


def _covered_pairs(original: TangleSession, loaded: TangleSession):
    """Every object the save covers, with its loaded counterpart."""
    a, b = original.workbench, loaded.workbench
    yield original, loaded
    yield original.bridge_alphabet, loaded.bridge_alphabet
    yield from zip(original.resonance_zones.values(), loaded.resonance_zones.values(), strict=True)
    yield a, b
    yield a.point_store, b.point_store
    yield from zip(a.point_store.chains, b.point_store.chains, strict=True)
    yield from zip(a.fixed_points, b.fixed_points, strict=True)
    for curves in ("_manifolds", "_bridges"):
        yield from zip(getattr(a, curves).values(), getattr(b, curves).values(), strict=True)
    for curves in ("_partial_bridges", "_image_curves"):
        yield from zip(getattr(a, curves), getattr(b, curves), strict=True)
    yield a.intersection_registry, b.intersection_registry
    yield a.intersection_registry.iterate_table, b.intersection_registry.iterate_table
    yield from ((x, y) for (_, x), (_, y) in zip(a.intersection_registry, b.intersection_registry, strict=True))
    yield a.Tangle, b.Tangle
    yield a._man_machine, b._man_machine


def test_a_load_restores_every_attribute(saved_session):
    case, loaded, _ = saved_session
    for original, copy in _covered_pairs(case.session, loaded):
        cls = type(original)
        assert type(copy) is cls
        saved, not_saved = _SAVED[cls], _NOT_SAVED.get(cls, frozenset())
        before, after = _attributes(original), _attributes(copy)
        unlisted = before.keys() - saved - not_saved
        assert not unlisted, f"{cls.__name__}.{sorted(unlisted)}: save it (bump FORMAT_VERSION) or list it"
        assert saved <= before.keys(), f"{cls.__name__} no longer has {sorted(saved - before.keys())}"
        assert after.keys() - not_saved == before.keys() - not_saved, cls.__name__
        assert {n: after[n] for n in saved} == {n: before[n] for n in saved}, cls.__name__
