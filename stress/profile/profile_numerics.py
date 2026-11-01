"""
Memory and time profile of the numerics pipeline, by component and by stage.

The baseline of the numerics-arrays campaign (N0): the same script, rerun
unchanged on a later branch, draws the same figures, so the move from the
linked point lists to arrays is judged on measured bytes and seconds.

Usage (from the repository root)::

    env/bin/python stress/profile/profile_numerics.py            # every case, figures + json
    env/bin/python stress/profile/profile_numerics.py --cases k10 scale
    env/bin/python stress/profile/profile_numerics.py --plot-only --out DIR
    env/bin/python stress/profile/profile_numerics.py --compare BASE_DIR --out BRANCH_DIR
    env/bin/python stress/profile/profile_numerics.py --cases scale --modes time --workers 8
    env/bin/python stress/profile/profile_numerics.py --cases scale --lean --separation width

Output: ``figures/numerics_profile/<git-sha>/`` (``--out`` overrides): one
``<case>_<mode>.json`` per run and the figures (``stage_time.png``,
``memory_components.png``, ``scaling_bytes_per_point.png``,
``scaling_time_rss.png``, ``scaling_points.png``, ``stage_peak_rss.png``).

Each case runs in a FRESH subprocess per mode, so one case's memory never
leaks into another's:

* ``time`` -- the stages are timed by wrapping library functions (exclusive
  self time per call path, so ``grow > iterate > refine > map`` splits growth
  into map, refine, merge and walks), RSS is sampled every 20 ms and its peak
  recorded per top-level stage, and the memory census runs at the end (and,
  for the scaling case, after every blast; each scaling row also keeps the
  blast's ``same_stability_count``, caught by wrapping ``report_resolution``,
  so two runs compare the count blast by blast). The census is structural: it
  walks the objects the workbench and the session hold and sizes them with
  ``sys.getsizeof``.
* ``memory`` -- the same build under ``tracemalloc``: per top-level stage the
  Python-heap peak, at the end the traced total (the reference the census must
  sum to) and the top allocation sites.

``--lean`` blasts the scaling case lean (``blast_zone(lean=True)``, no
diagnose) and ``--separation`` sets its ``min_separation`` (a float,
``width`` or ``none``; ``1e-5`` by default, ``width`` with ``--lean``, which
refuses a float); each scaling row records the rows
the blast freed (``released_rows``) and its drops (``too_close``).

``--workers N`` builds every session with ``workers=N`` (the labels stay; the
json records ``workers``; at 1 a checkout without the option runs as before).
Only the main thread is timed: a split map batch counts under ``map``, the
sweep's threads under ``diagnose``.

Components of the census (bytes):

* ``points_curve`` -- point objects on a registered manifold (object,
  ``__dict__``, the ``_coords`` array, the x / y / cdist scalars; since the
  point store, a slotted object and its row ``int``); a point is told by its
  row (since N7a the engine makes none: these were the ends read through the
  adapter; since N7b a row has no object and this reads 0);
* ``points_image`` -- point objects whose row lies only on blast image curves
  (0 since N7b);
* ``points_phantom`` -- the coordinate-only preimages refinement caches
  (``cdist`` None, on no curve) while they were objects; since the point store
  they are rows with no object (``phantom_rows``), and this counts the objects
  with no cdist left (before N7a, the initializer's detached fictitious tails);
* ``points_other`` -- any other point object alive (separators, stale lists,
  merge-dropped duplicates; 0 since N7b);
* ``point_store`` -- the workbench's ``PointStore``: its columns at capacity
  (the chain and position columns too, since the order arrays) and, before
  N7b, its row -> point list. Its counts are of rows, so they read alike on
  every tree with a store: ``store_rows``, ``store_capacity``, ``curve_rows`` (on a
  walked manifold), ``image_rows`` (only on blast image curves),
  ``no_cdist_rows``, ``phantom_rows`` (no cdist, no link of their own) and
  ``off_curve_rows`` (a cdist, on no walked curve); ``n_points`` and the
  other ``n_*`` totals are these rows;
* ``curve_ids`` -- the store's ``Chain`` objects and their ``ids`` arrays, the
  order of points along every curve (count ``chains``, the live ones); 0 on a
  tree without them;
* ``branch_points`` -- the ``BranchPoint`` objects: the object, its
  ``__dict__`` and the lists or arrays in it (none since N8, where a
  ``BranchPoint`` is a view of its ``FixedPoint``);
* ``segment_index`` -- ``Tangle``'s Python side: ``_Segment`` objects, the
  edge / owner / lookup dicts and sets (the query-only edge set counted with
  ``edge_to_sid``), the crossing-pair sets;
* ``rtree`` -- the two rtrees' C memory (libspatialindex, invisible to
  tracemalloc): measured at the end by dropping them and reading RSS, and
  estimated per blast as entries x the calibrated bytes per entry;
* ``curve_objects`` -- manifold / bridge objects and their attributes;
* ``walk_memo`` -- left the census at N7a (N6 removed the bridges' walk
  memo); its key stays in ``profile_figures.MEMORY_GROUPS``, so an older
  record's memo bytes still stack;
* ``registry`` -- the intersection registry (crossings, orders, iterate table);
* ``topology`` -- the session's trellis and topology caches, alphabet, zones;
* ``workbench_other`` -- everything else the workbench holds, the fixed
  points' slot tables included (since N8); the walk stops at a ``Chain``,
  which ``curve_ids`` counts.

Dev Notes:
    The wrapped functions are listed in :data:`STAGES`; a renamed library
    function makes the run fail loudly (``AttributeError``) rather than drop a
    stage silently. When a later slice renames one on purpose, add the new name
    beside the old under the SAME label so the figures stay comparable.
"""

from __future__ import annotations

import argparse
import contextlib
import ctypes
import dataclasses
import functools
import gc
import importlib
import json
import logging
import os
import subprocess
import sys
import threading
import time
import tracemalloc
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable

import numpy as np
import psutil

REPO = Path(os.environ.get("TANGLEPACK_REPO", Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "tests"))

LAW_PROFILE_CASES = ("k10", "p3", "nested", "inversion", "heteroclinic")
ALL_CASES = LAW_PROFILE_CASES + ("scale",)
#: The scaling case: the k = 2.5 period-1 saddle (scripts/henon_k25_entropy.py's
#: ``p1``), blasted one single-iteration step at a time, topology after each.
SCALE_K = 2.5
SCALE_SEED = (4.0, -4.0)
SCALE_MIN_SEPARATION = 1e-5
#: Calibrated by building a stream-loaded 2-D rtree of 1e6 boxes (rtree 1.4).
RTREE_BYTES_PER_ENTRY = 68.0

#: Imported before the baseline (see ``run_worker``).
PRELOAD = ("tanglepack", "tanglepack.topology", "tanglepack.loom", "tanglepack.examples.henon",
           "tanglepack.examples.henon_cases", "cases", "networkx", "scipy.spatial",
           "matplotlib.pyplot", "tanglepack.topology.plotting")

# (module, attribute path, label). Top-level stages first; their labels head
# every call path. Inner labels split them.
TOP_STAGES = [
    ("tanglepack.numerics.TangleWorkbench", "TangleWorkbench.construct_fixed_point", "fixed_point"),
    ("tanglepack.numerics.TangleWorkbench", "TangleWorkbench.orient_eigenvectors", "fixed_point"),
    ("tanglepack.numerics.TangleWorkbench", "TangleWorkbench.initialize_both_manifolds", "initialise"),
    ("tanglepack.numerics.TangleWorkbench", "TangleWorkbench.grow_n_times", "grow"),
    ("tanglepack.numerics.TangleWorkbench", "TangleWorkbench.grow_until_turnaround", "grow"),
    ("tanglepack.numerics.TangleWorkbench", "TangleWorkbench.compute_intersections", "intersections"),
    ("tanglepack.numerics.TangleWorkbench", "TangleWorkbench.trim_stable_manifolds", "trim"),
    ("tanglepack.numerics.TangleWorkbench", "TangleWorkbench.create_bridges", "create_bridges"),
    ("tanglepack.loom.TangleSession", "TangleSession.classify_strong_pips", "pips"),
    ("tanglepack.loom.TangleSession", "TangleSession.resonance_zone", "zone"),
    ("tanglepack.loom.TangleSession", "TangleSession.add_resonance_zones", "zone"),
    ("tanglepack.loom.TangleSession", "TangleSession.blast_zone", "blast"),
    ("tanglepack.loom.TangleSession", "TangleSession.topology", "topology"),
]
INNER_STAGES = [
    ("tanglepack.numerics.DynamicalSystem", "DynamicalSystem.map_batch", "map"),
    ("tanglepack.numerics.DynamicalSystem", "DynamicalSystem.map_inv_batch", "map"),
    ("tanglepack.numerics.ManifoldMachine", "ManifoldMachine._iterate_manifold", "iterate"),
    ("tanglepack.numerics.ManifoldMachine", "ManifoldMachine.merge_manifolds", "merge"),
    ("tanglepack.numerics.ManifoldMachine", "ManifoldMachine.refine_manifold", "refine"),
    ("tanglepack.numerics.ManifoldMachine", "ManifoldMachine._image_chains", "image_chains"),
    ("tanglepack.numerics.BaseManifold", "BaseManifold._collect", "walk"),
    ("tanglepack.numerics.BaseManifold", "BaseManifold._find_tail", "walk"),
    ("tanglepack.numerics.Tangle", "Tangle.add_manifolds", "index"),
    ("tanglepack.numerics.Tangle", "Tangle.add_manifold", "index"),
    ("tanglepack.numerics.Tangle", "Tangle.resolve_crossings", "resolve"),
    ("tanglepack.numerics.Tangle", "Tangle.populate_intersections", "resolve"),
    ("tanglepack.numerics.Tangle", "Tangle.create_bridges", "cut"),
    ("tanglepack.numerics.IntersectionRegistry", "IntersectionRegistry.add", "register"),
    ("tanglepack.numerics.IntersectionRegistry", "IntersectionRegistry.resolved_ids", "resolved_ids"),
    ("tanglepack.numerics.TangleWorkbench", "TangleWorkbench._register_anchors", "register"),
    ("tanglepack.numerics.TangleWorkbench", "TangleWorkbench._register_bridge", "register_bridge"),
    ("tanglepack.numerics.TangleWorkbench", "TangleWorkbench.infer_iterates", "infer_iterates"),
    ("tanglepack.numerics.TangleWorkbench", "TangleWorkbench.report_resolution", "diagnose"),
    ("tanglepack.numerics.TangleWorkbench", "TangleWorkbench.iterate_bridge", "iterate_bridge"),
    ("tanglepack.numerics.TangleWorkbench", "TangleWorkbench._bump_generation", "generation"),
    ("tanglepack.topology.Trellis", "Trellis.from_workbench", "trellis_build"),
    ("tanglepack.loom.Blast", "_interior_points", "separation"),
    ("tanglepack.loom.Blast", "cKDTree", "separation"),
]


# --------------------------------------------------------------------------- #
# Stage timer: exclusive self time per call path, RSS peak per top-level stage
# --------------------------------------------------------------------------- #
class StageTimer:
    """Exclusive (self) time and call count per call path of stage labels."""

    def __init__(self, trace: bool) -> None:
        self.stack: list[list] = []  # [label, start, child_time]
        self.self_time: dict[tuple, float] = defaultdict(float)
        self.calls: dict[tuple, int] = defaultdict(int)
        self.top_events: list[dict] = []
        self.trace = trace
        self.proc = psutil.Process()
        self._window_peak = 0
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._sample, daemon=True)
        self.peak_rss = 0

    def _sample(self) -> None:
        while not self._stop.is_set():
            with contextlib.suppress(Exception):
                rss = self.proc.memory_info().rss
                self._window_peak = max(self._window_peak, rss)
                self.peak_rss = max(self.peak_rss, rss)
            self._stop.wait(0.02)

    def start(self) -> None:
        """Start the RSS sampler."""
        self._thread.start()

    def stop(self) -> None:
        """Stop the RSS sampler."""
        self._stop.set()
        self._thread.join(timeout=1)

    def enter(self, label: str) -> None:
        """Open a stage."""
        if not self.stack:
            self._window_peak = self.proc.memory_info().rss
            if self.trace:
                tracemalloc.reset_peak()
        self.stack.append([label, time.perf_counter(), 0.0])

    def exit(self) -> None:
        """Close the innermost stage."""
        end = time.perf_counter()
        path = tuple(frame[0] for frame in self.stack)
        label, start, child = self.stack.pop()
        elapsed = end - start
        self.self_time[path] += elapsed - child
        self.calls[path] += 1
        if self.stack:
            self.stack[-1][2] += elapsed
        else:
            rss_now = self.proc.memory_info().rss
            event = {"stage": label, "seconds": elapsed,
                     "peak_rss": max(self._window_peak, rss_now), "rss_after": rss_now}
            if self.trace:
                current, peak = tracemalloc.get_traced_memory()
                event.update(traced_after=current, traced_peak=peak)
            self.top_events.append(event)

    def wrap(self, fn: Callable, label: str) -> Callable:
        """``fn`` timed under ``label``."""
        timer = self

        @functools.wraps(fn)
        def timed(*args, **kwargs):
            timer.enter(label)
            try:
                return fn(*args, **kwargs)
            finally:
                timer.exit()

        return timed

    def summary(self) -> dict:
        """Self time per path, and per top-level stage x inner label."""
        by_stage: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
        for path, seconds in self.self_time.items():
            by_stage[path[0]][path[-1] if len(path) > 1 else "(self)"] += seconds
        return {
            "paths": [{"path": "/".join(p), "self_s": s, "calls": self.calls[p]}
                      for p, s in sorted(self.self_time.items(), key=lambda kv: -kv[1])],
            "by_stage": {k: dict(v) for k, v in by_stage.items()},
            "stage_total": {k: sum(v.values()) for k, v in by_stage.items()},
            "top_events": self.top_events,
        }


def install(timer: StageTimer) -> None:
    """Wrap every function of :data:`TOP_STAGES` and :data:`INNER_STAGES`."""
    for module_name, attr_path, label in TOP_STAGES + INNER_STAGES:
        module = importlib.import_module(module_name)
        owner = module
        *owners, name = attr_path.split(".")
        for part in owners:
            owner = getattr(owner, part)
        original = getattr(owner, name)  # AttributeError on a rename: loud
        if isinstance(original, (staticmethod, classmethod)):
            raise TypeError(f"{attr_path}: wrap the function, not a descriptor")
        if owners and isinstance(owner.__dict__.get(name), classmethod):
            fn = owner.__dict__[name].__func__
            setattr(owner, name, classmethod(timer.wrap(fn, label)))
        else:
            setattr(owner, name, timer.wrap(original, label))


def record_counts() -> list:
    """Wrap ``report_resolution`` (outside its timer) to collect every blast's
    ``same_stability_count``; ``blast_zone`` returns no report."""
    from tanglepack.numerics.TangleWorkbench import TangleWorkbench

    counts = []
    report_resolution = TangleWorkbench.report_resolution

    def recording(self, *args, **kwargs):
        report = report_resolution(self, *args, **kwargs)
        if report.same_stability_count is not None:
            counts.append(report.same_stability_count)
        return report

    TangleWorkbench.report_resolution = recording
    return counts


# --------------------------------------------------------------------------- #
# Memory census
# --------------------------------------------------------------------------- #
_libc = None


def malloc_trim() -> None:
    """Return freed heap pages to the OS (glibc), so RSS reads what is held."""
    global _libc
    with contextlib.suppress(Exception):
        if _libc is None:
            _libc = ctypes.CDLL("libc.so.6")
        _libc.malloc_trim(0)


def _walk_ids(manifold: Any) -> list[int]:
    """The store rows of a manifold, root to its PHYSICAL end (the end of its
    chain), the periodic point (no row) left out.

    Read off the chain (since the order arrays) rather than ``_collect``,
    which the timer wraps, and past a lean blast's read guard: a released
    bridge's rows are what is left of it."""
    from tanglepack.numerics.BaseManifold import BaseManifold

    _head, chain, lo, hi = BaseManifold._span(manifold, stop_at_final=False)
    return [] if chain is None else chain.ids[lo:hi].tolist()


def _deep_size(roots: list, seen: set, stop: Callable[[Any], bool]) -> int:
    """``sys.getsizeof`` over everything reachable from ``roots`` not stopped or seen."""
    import types

    skip_types = (type, types.ModuleType, types.FunctionType, types.BuiltinFunctionType,
                  types.MethodType, types.CodeType, logging.Logger, threading.Thread,
                  np.ufunc, np.dtype)
    total = 0
    stack = list(roots)
    while stack:
        obj = stack.pop()
        oid = id(obj)
        if oid in seen or isinstance(obj, skip_types) or stop(obj):
            continue
        seen.add(oid)
        total += sys.getsizeof(obj)
        if isinstance(obj, dict):
            stack.extend(obj.keys())
            stack.extend(obj.values())
        elif isinstance(obj, (list, tuple, set, frozenset)):
            stack.extend(obj)
        elif isinstance(obj, np.ndarray):
            if obj.base is not None:
                stack.append(obj.base)
        elif isinstance(obj, (str, bytes, int, float, complex, bool, np.generic)):
            pass
        else:
            d = getattr(obj, "__dict__", None)
            if isinstance(d, dict):
                stack.append(d)
            for cls in type(obj).__mro__:
                for slot in getattr(cls, "__slots__", ()):
                    if slot not in ("__dict__", "__weakref__"):
                        with contextlib.suppress(AttributeError):
                            stack.append(getattr(obj, slot))
    return total


def _tangle_bytes(tangle: Any) -> dict:
    """The segment index's Python side, sized structurally (no seen-set)."""
    out = {"segments": 0, "seg_lookup": 0, "edge_to_sid": 0, "owners": 0, "pairs": 0}
    lookup = tangle._seg_lookup
    out["seg_lookup"] = sys.getsizeof(lookup)
    for sid, seg in lookup.items():
        out["segments"] += sys.getsizeof(seg) + sys.getsizeof(sid)
    edges = tangle._edge_to_sid
    out["edge_to_sid"] = sys.getsizeof(edges) + sum(sys.getsizeof(k) for k in edges)
    query_edges = tangle._query_edges
    out["edge_to_sid"] += sys.getsizeof(query_edges) + sum(sys.getsizeof(k) for k in query_edges)
    owners = sys.getsizeof(tangle._manifold_segs) + sys.getsizeof(tangle._seg_manifolds)
    owners += sum(sys.getsizeof(s) for s in tangle._manifold_segs.values())
    owners += sum(sys.getsizeof(s) for s in tangle._seg_manifolds.values())
    out["owners"] = owners
    for pairs in (tangle._intersecting_segments, tangle._processed_pairs):
        out["pairs"] += sys.getsizeof(pairs) + sum(
            sys.getsizeof(fs) + sum(sys.getsizeof(i) for i in fs) for fs in pairs
        )
    out["rtree_entries"] = {s: int(idx.get_size()) for s, idx in tangle._rtrees.items()}
    out["n_segments"] = len(lookup)
    out["n_query_edges"] = len(query_edges)
    return out


def census(session: Any, *, measure_rtree: bool = False) -> dict:
    """Bytes by component of everything the session holds (see the module docstring)."""
    from tanglepack.numerics.BranchPoint import BranchPoint
    from tanglepack.numerics.BaseManifold import BaseManifold
    from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry
    from tanglepack.numerics.FixedPoint import FixedPoint
    from tanglepack.numerics.PointStore import Chain

    t0 = time.perf_counter()
    workbench = session.workbench
    store = getattr(workbench, "point_store", None)
    if store is None or not hasattr(BaseManifold, "_rows"):
        raise RuntimeError(
            f"the census reads a point store and curve rows (the order arrays, N4, on): "
            f"{REPO} has none; compare its committed census JSON instead (--compare)"
        )
    gc.collect()

    curves = list(workbench.manifolds.values())
    bridges = list(workbench._bridges.values()) + list(workbench._partial_bridges)
    images = list(workbench._image_curves)
    on_curve: set[int] = set()
    for manifold in curves:
        on_curve.update(_walk_ids(manifold))
    on_image: set[int] = set()
    for curve in images:
        on_image.update(_walk_ids(curve))

    counts = defaultdict(int)
    sizes = defaultdict(int)
    for obj in gc.get_objects():
        if isinstance(obj, BranchPoint):
            counts["branch_points"] += 1
            sizes["branch_points"] += sys.getsizeof(obj) + sys.getsizeof(obj.__dict__) + sum(
                sys.getsizeof(v) for v in obj.__dict__.values()
                if isinstance(v, (list, tuple, np.ndarray))
            )
    rows = len(store)
    on_a_curve = np.zeros(rows, bool)
    on_a_curve[list(on_curve | on_image)] = True
    has_cdist = ~np.isnan(store._cdist[:rows])
    own_link = (store._image[:rows] >= 0) | (store._preimage[:rows] >= 0)
    counts["store_rows"] = rows
    counts["store_capacity"] = len(store._cdist)
    counts["curve_rows"] = len(on_curve)
    counts["image_rows"] = len(on_image - on_curve)
    counts["no_cdist_rows"] = int((~has_cdist).sum())
    counts["phantom_rows"] = int((~has_cdist & ~own_link).sum())
    counts["off_curve_rows"] = int((has_cdist & ~on_a_curve).sum())
    nodes = getattr(store, "nodes", None)  # the row -> point list of a tree before N7b
    sizes["point_store"] = store.nbytes + (0 if nodes is None else sys.getsizeof(nodes))
    chains = getattr(store, "chains", None)
    if chains is not None:
        counts["chains"] = sum(1 for chain in chains if len(chain.ids))
        sizes["curve_ids"] = sys.getsizeof(chains) + sum(
            sys.getsizeof(chain) + sys.getsizeof(chain.ids) for chain in chains
        )
    del on_curve, on_image

    tangle = workbench.Tangle
    tangle_parts = _tangle_bytes(tangle)
    sizes["segment_index"] = sum(v for k, v in tangle_parts.items()
                                 if k not in ("rtree_entries", "n_segments", "n_query_edges"))
    rtree_entries = sum(tangle_parts["rtree_entries"].values())
    sizes["rtree"] = int(rtree_entries * RTREE_BYTES_PER_ENTRY)

    held_curves = {id(c) for c in curves + bridges + images}
    seen: set[int] = set()
    # curve objects (the manifold / bridge objects and their dicts)
    curve_bytes = 0
    for c in {id(c): c for c in curves + bridges + images}.values():
        curve_bytes += sys.getsizeof(c) + sys.getsizeof(c.__dict__)
        seen.update((id(c), id(c.__dict__)))
    sizes["curve_objects"] = curve_bytes

    def stop_numerics(o: Any) -> bool:
        return (isinstance(o, BranchPoint) or o is store or (isinstance(o, BaseManifold) and id(o) in held_curves)
                or o is workbench or o is tangle or isinstance(o, FixedPoint)
                or type(o).__name__ in ("Index", "Property", "DynamicalSystem",
                                        "ManifoldMachine", "ManifoldInitializer",
                                        "FixedPointSolver"))

    seen.add(id(tangle))
    registry = workbench.intersection_registry
    sizes["registry"] = _deep_size([registry], seen, stop_numerics)

    topo_roots = [session._trellises, session._topologies, session.bridge_alphabet,
                  session.resonance_zones, session._strong_pips]

    def stop_topology(o: Any) -> bool:
        return stop_numerics(o) or isinstance(o, (IntersectionRegistry,))

    sizes["topology"] = _deep_size(topo_roots, seen, stop_topology)
    seen.add(id(session.__dict__))
    # The fixed points' slot tables name chains, which ``curve_ids`` counts.
    sizes["workbench_other"] = _deep_size([workbench.__dict__], seen,
                                          lambda o: isinstance(o, (BranchPoint, Chain))
                                          or o is tangle or o is store)

    result = {
        "bytes": dict(sizes),
        "counts": dict(counts),
        "tangle_parts": {k: v for k, v in tangle_parts.items()},
        "n_points": counts["store_rows"],
        "n_curve_points": counts["curve_rows"],
        "n_image_points": counts["image_rows"],
        "n_phantom": counts["no_cdist_rows"],
        "n_bridges": len(bridges),
        "n_image_curves": len(images),
        "n_crossings": len(registry),
        "census_s": time.perf_counter() - t0,
    }
    if measure_rtree:
        proc = psutil.Process()
        gc.collect()
        malloc_trim()
        before = proc.memory_info().rss
        saved = tangle._rtrees
        from tanglepack.numerics.Tangle import _segment_index

        tangle._rtrees = {"unstable": _segment_index(), "stable": _segment_index()}
        del saved
        gc.collect()
        malloc_trim()
        after = proc.memory_info().rss
        result["rtree_measured_bytes"] = before - after
    return result


# --------------------------------------------------------------------------- #
# Cases
# --------------------------------------------------------------------------- #
def set_workers(workers: int) -> None:
    """Make ``workers`` the default of every ``TangleSession`` the run builds."""
    from tanglepack.loom.TangleSession import TangleSession

    TangleSession.__init__ = functools.partialmethod(TangleSession.__init__, workers=workers)


def build_law_case(name: str) -> Any:
    """Build one case through its frozen ``tests/cases.py`` builder."""
    import cases

    return cases.BUILDERS[name]().session


def force_topology(session: Any) -> None:
    """``session.topology()`` and its transition matrix, as the sweep scripts read it."""
    topology = session.topology()
    topology.dynamics.transition_matrix(refined=True)


def scale_build(area_cutoff: float) -> tuple[Any, Any, Any]:
    """The k = 2.5 period-1 recipe of ``scripts/henon_bridge_classes.py``, unblasted."""
    from tanglepack import TangleSession
    from tanglepack.examples.henon import henon_jacobian, henon_map, henon_map_inverse

    session = TangleSession(henon_map(SCALE_K, 1), henon_map_inverse(SCALE_K, 1),
                            henon_jacobian(SCALE_K, 1), area_cutoff=area_cutoff)
    fp = session.construct_fixed_point(list(SCALE_SEED))
    session.orient_eigenvectors(fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])})
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=10)
    session.grow_until_turnaround(fp, "stable")
    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    session.classify_strong_pips()
    trellis = session.trellis(fp)
    pip = trellis.iterate(trellis.strong_pip, 1)
    zone = session.resonance_zone(pip)
    return session, fp, zone


# --------------------------------------------------------------------------- #
# Worker
# --------------------------------------------------------------------------- #
def run_worker(case: str, mode: str, args: argparse.Namespace, out: Path) -> None:
    """One (case, mode) in this process; writes ``<case>_<mode>.json``."""
    logging.basicConfig(level=logging.ERROR)
    trace = mode == "memory"
    # Everything a build imports, imported BEFORE the baseline: the traced and
    # RSS growth then measure what the pipeline holds, not module code.
    for name in PRELOAD:
        importlib.import_module(name)
    if args.workers != 1:
        set_workers(args.workers)
    timer = StageTimer(trace)
    install(timer)
    gc.collect()
    malloc_trim()
    proc = psutil.Process()
    rss0 = proc.memory_info().rss
    if trace:
        tracemalloc.start(1)
    traced0 = tracemalloc.get_traced_memory()[0] if trace else None
    timer.start()
    record: dict = {"case": case, "mode": mode, "git_sha": git_sha(), "workers": args.workers,
                    "rss_baseline": rss0,
                    "traced_baseline": traced0, "series": [], "stopped": None}
    path = out / f"{case}_{mode}.json"
    t0 = time.perf_counter()

    def save() -> None:
        record["timer"] = timer.summary()
        record["wall_s"] = time.perf_counter() - t0
        record["peak_rss"] = timer.peak_rss
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps(record, indent=1, default=_jsonable))
        os.replace(tmp, path)

    try:
        if case != "scale":
            session = build_law_case(case)
            record["wall_build_s"] = time.perf_counter() - t0
        else:
            text = args.separation
            separation = None if text == "none" else text if text == "width" else float(text)
            record["scale"] = {"k": SCALE_K, "area_cutoff": args.scale_cutoff,
                               "min_separation": separation, "lean": args.lean,
                               "stop_gb": args.scale_stop_gb, "max_blasts": args.scale_max_blasts}
            counts = record_counts()
            session, fp, zone = scale_build(args.scale_cutoff)
            timer.enter("topology")
            try:
                force_topology(session)
            finally:
                timer.exit()
            max_blasts = args.scale_max_blasts if not trace else min(args.scale_max_blasts,
                                                                      args.scale_trace_blasts)
            flat = 0
            last_points = None
            for blast in range(max_blasts + 1):
                if blast:
                    t_blast = time.perf_counter()
                    result = session.blast_zone(zone, num_iterations=1, fixed_point=[fp],
                                                min_separation=separation,
                                                diagnose=not args.lean, lean=args.lean)
                    blast_s = time.perf_counter() - t_blast
                    timer.enter("topology")
                    try:
                        force_topology(session)
                    finally:
                        timer.exit()
                    step_s = time.perf_counter() - t_blast
                else:
                    result, blast_s, step_s = None, 0.0, 0.0
                events = [e for e in timer.top_events[-2:]] if blast else timer.top_events
                row = {
                    "blast": blast,
                    "blast_s": blast_s,
                    "step_s": step_s,
                    "wall_s": time.perf_counter() - t0,
                    "rss": proc.memory_info().rss,
                    "peak_rss_so_far": timer.peak_rss,
                    "step_peak_rss": max(e["peak_rss"] for e in events),
                    "pip_depth": ({str(k.label): v for k, v in result.pip_depth.items()}
                                  if result is not None else None),
                    "points_created": result.points_created if result is not None else 0,
                    "released_rows": result.released_rows if result is not None else 0,
                    "too_close": result.too_close if result is not None else 0,
                    "same_stability_count": (dataclasses.asdict(counts[-1])
                                             if blast and counts else None),
                }
                counts.clear()
                if trace:
                    row["traced"], row["traced_peak"] = tracemalloc.get_traced_memory()
                if not trace or blast == max_blasts:
                    row["census"] = census(session)
                record["series"].append(row)
                save()
                n_points = row.get("census", {}).get("n_points")
                print(f"[scale {args.scale_cutoff:g} {mode}] blast {blast}: "
                      f"{n_points} points, rss {row['rss'] / 2**30:.2f} GB, "
                      f"step peak {row['step_peak_rss'] / 2**30:.2f} GB, "
                      f"{row['step_s']:.1f} s", flush=True)
                if row["step_peak_rss"] >= args.scale_stop_gb * 2**30:
                    record["stopped"] = f"peak RSS crossed {args.scale_stop_gb} GB at blast {blast}"
                    break
                if n_points is not None:
                    flat = flat + 1 if n_points == last_points else 0
                    last_points = n_points
                    if flat >= 2:
                        record["stopped"] = f"plateau: no new points since blast {blast - 2}"
                        break
        if trace:
            record["traced_final"], record["traced_peak"] = tracemalloc.get_traced_memory()
            snapshot = tracemalloc.take_snapshot()
            record["top_files"] = [
                {"file": str(s.traceback[0].filename).replace(str(REPO) + "/", ""),
                 "bytes": s.size, "count": s.count}
                for s in snapshot.statistics("filename")[:15]
            ]
            noise = ("<frozen", "importlib", "linecache.py", "inspect.py", "tracemalloc.py")
            record["traced_import_noise"] = sum(
                s.size for s in snapshot.statistics("filename")
                if any(n in str(s.traceback[0].filename) for n in noise))
            record["top_lines"] = [
                {"site": f"{str(s.traceback[0].filename).replace(str(REPO) + '/', '')}:"
                         f"{s.traceback[0].lineno}", "bytes": s.size, "count": s.count}
                for s in snapshot.statistics("lineno")[:25]
            ]
            del snapshot
            record["census"] = census(session)
        else:
            record["rss_final"] = proc.memory_info().rss
            record["census"] = census(session, measure_rtree=True)
    except Exception as exc:  # noqa: BLE001 -- a failure is data
        import traceback

        record["stopped"] = f"{type(exc).__name__}: {exc}"
        record["traceback"] = traceback.format_exc()[-3000:]
    timer.stop()
    save()


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def git_sha() -> str:
    """The checkout's short sha (with ``-dirty`` when src/ has changes)."""
    with contextlib.suppress(Exception):
        sha = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=REPO,
                                      text=True).strip()
        dirty = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", "src"], cwd=REPO)
        return sha + ("-dirty" if dirty.returncode else "")
    return "unknown"


# --------------------------------------------------------------------------- #
# Parent
# --------------------------------------------------------------------------- #
def launch(case: str, mode: str, args: argparse.Namespace, out: Path) -> int:
    """Run one worker in a fresh interpreter, killing it past ``--cap-gb``."""
    cmd = [sys.executable, str(Path(__file__).resolve()), "--worker", case, mode,
           "--out", str(out), "--scale-cutoff", repr(args.scale_cutoff),
           "--scale-stop-gb", repr(args.scale_stop_gb),
           "--scale-max-blasts", str(args.scale_max_blasts),
           "--scale-trace-blasts", str(args.scale_trace_blasts), "--workers", str(args.workers),
           "--separation", args.separation, *(["--lean"] if args.lean else [])]
    child = subprocess.Popen(cmd, cwd=REPO)
    handle = psutil.Process(child.pid)
    while child.poll() is None:
        with contextlib.suppress(psutil.Error):
            if handle.memory_info().rss > args.cap_gb * 2**30:
                child.kill()
                note = out / f"{case}_{mode}.json"
                if note.exists():
                    data = json.loads(note.read_text())
                    data["stopped"] = f"killed by the {args.cap_gb} GB cap"
                    note.write_text(json.dumps(data, indent=1))
                break
        time.sleep(0.2)
    return child.wait()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--cases", nargs="+", default=list(ALL_CASES), choices=ALL_CASES)
    parser.add_argument("--modes", nargs="+", default=["time", "memory"],
                        choices=["time", "memory"])
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--scale-cutoff", type=float, default=1e-7)
    parser.add_argument("--scale-stop-gb", type=float, default=8.0,
                        help="stop the scaling series after the blast whose peak RSS crosses this")
    parser.add_argument("--scale-max-blasts", type=int, default=40)
    parser.add_argument("--scale-trace-blasts", type=int, default=8,
                        help="blasts of the scaling case run under tracemalloc (memory mode)")
    parser.add_argument("--cap-gb", type=float, default=12.0, help="kill a worker past this RSS")
    parser.add_argument("--workers", type=int, default=1, help="every session's workers")
    parser.add_argument("--lean", action="store_true", help="blast the scaling case lean")
    parser.add_argument("--separation", default=None,
                        help="the scaling blasts' min_separation: a float, width or none "
                             f"(default {SCALE_MIN_SEPARATION!r}, width with --lean)")
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--compare", type=Path, default=None,
                        help="a baseline output dir: figures draw it beside --out")
    parser.add_argument("--worker", nargs=2, metavar=("CASE", "MODE"), help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.separation is None:
        args.separation = "width" if args.lean else repr(SCALE_MIN_SEPARATION)
    if args.lean and args.separation not in ("width", "none"):
        parser.error("--lean takes --separation width or none: a float separation reads every interior")

    out = args.out or REPO / "figures" / "numerics_profile" / git_sha()
    out.mkdir(parents=True, exist_ok=True)
    if args.worker:
        run_worker(args.worker[0], args.worker[1], args, out)
        return
    if not args.plot_only:
        for case in args.cases:
            for mode in args.modes:
                t = time.perf_counter()
                code = launch(case, mode, args, out)
                print(f"{case} {mode}: exit {code}, {time.perf_counter() - t:.0f} s", flush=True)
    import profile_figures  # noqa: E402  (beside this file)

    profile_figures.draw(out, args.compare)


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
