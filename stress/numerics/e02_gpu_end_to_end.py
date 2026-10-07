"""E2: end-to-end GPU vs CPU on the real growth path.

Cases (unstable growth, then the stable side, crossings and bridges):

* ``k28``  -- k=2.8, ``area_cutoff = 1e-7``, 12 unstable steps (~4e5 points,
  ~6 s for the last step on the CPU; step 13 would be ~1e7 points), stable to
  turnaround. Past ~10 steps the unstable arm escapes to infinity
  (max|coord| ~1e6 at step 12), so most of the points are the escaping tip.
* ``p3``   -- the period-3 orbit of (2, 1), ``area_cutoff = 1e-7``, 16 unstable
  steps (~4e3 points; step 17 jumps to ~6e5 points and max|coord| ~1.6e6,
  the escape), 9 stable steps.
* ``k10``  -- k=10 at the default cutoff (1e-4), 9 steps: the cheap control.

Config kinds:

* ``timed``   -- (case, backend, repeat): per-step growth time and size series,
  then each later stage timed once.
* ``profile`` -- (case, backend): cProfile around the unstable growth; the
  top-40 functions by tottime plus the tottime of EVERY profiled function
  summed into categories (map / refinement / linked_list / intersections /
  other).
* ``mbp``     -- the ``min_batch_points`` sweep on k28 (1e9 = never on GPU).
* ``verify``  -- CPU and GPU sessions grown in lockstep in one process: point
  counts per step and the max coordinate deviation of the unstable manifolds.
"""

from __future__ import annotations

import cProfile
import gc
import pstats

import numpy as np

import harness
from maps import henon_session, max_abs_coord, total_points
from tanglepack.examples.henon_cases import PERIOD3_ORIENTATION
from maps import PERIOD1_ORIENTATION

EXPERIMENT = "e02_gpu_end_to_end"
TIMEOUT = 600
WORKERS = 1  # timings: never share the machine with another worker

#: case -> (k, b, area_cutoff, unstable steps full, unstable steps quick)
CASES = {
    "k28": (2.8, 1.0, 1e-7, 12, 10),
    "p3": (2.0, 1.0, 1e-7, 16, 13),
    "k10": (10.0, 1.0, None, 9, 9),
}
P3_SEED = [[0, 1], [-1, 0], [-1, 1]]
P3_STABLE_STEPS = 9
MIN_BATCH_SWEEP = [16, 64, 256, 1024, 4096, 10**9]


def configs(quick: bool) -> list[dict]:
    def steps(case):
        return CASES[case][4 if quick else 3]

    repeats = 1 if quick else 3
    out = [
        {"kind": "timed", "case": c, "backend": b, "repeat": r, "steps": steps(c)}
        for c in CASES for b in ("cpu", "gpu") for r in range(repeats)
    ]
    out += [{"kind": "profile", "case": c, "backend": b, "steps": steps(c)}
            for c in CASES for b in ("cpu", "gpu")]
    sweep = [64, 10**9] if quick else MIN_BATCH_SWEEP
    out += [{"kind": "mbp", "case": "k28", "backend": "gpu", "min_batch_points": m,
             "steps": steps("k28")} for m in sweep]
    out += [{"kind": "verify", "case": c, "steps": steps(c)} for c in CASES]
    return out


def build(case: str, gpu: bool, min_batch_points: int = 64):
    """A seeded session for one case, its fixed point, and a warmed-up map."""
    k, b, cutoff, *_ = CASES[case]
    session = henon_session(k, b, area_cutoff=cutoff, gpu=gpu, min_batch_points=min_batch_points)
    if case == "p3":
        fp = session.construct_fixed_point(P3_SEED)
        session.orient_eigenvectors(fp, PERIOD3_ORIENTATION)
    else:
        fp = session.construct_fixed_point([4.0, -4.0])
        session.orient_eigenvectors(fp, PERIOD1_ORIENTATION)
    session.initialize_both_manifolds(fp)
    # CUDA context, kernel compile and allocator warm-up stay out of the timings.
    system = session.workbench.dynamical_system
    warm = np.random.default_rng(0).uniform(-1, 1, size=(4096, 2))
    system.map_batch(warm)
    system.map_inv_batch(warm)
    return session, fp


def grow_stable(session, case: str, fp) -> None:
    if case == "p3":
        session.grow_n_times(fp, "stable", num_iterations=P3_STABLE_STEPS)
    else:
        session.grow_until_turnaround(fp, "stable")


def grow_unstable(session, fp, steps: int, rec, profiler=None) -> None:
    """Grow one step at a time, recording the per-step series.

    The profiler, when given, is enabled around the growth call only, so the
    stage's GPU sync and the ``gc.collect`` stay out of the profile.
    """
    for _ in range(steps):
        gc.collect()
        before = rec.stages.get("grow_unstable", 0.0)
        with rec.stage("grow_unstable"):
            if profiler:
                profiler.enable()
            session.grow_n_times(fp, "unstable", num_iterations=1)
            if profiler:
                profiler.disable()
        rec.append("step_s", rec.stages["grow_unstable"] - before)
        rec.append("points_unstable", total_points(session, "unstable"))
        rec.append("max_abs_coord", max_abs_coord(session))
        rec.checkpoint()


def downstream(session, case: str, fp, rec) -> None:
    """Time stable growth, crossings, trim + bridges and the iterate table once."""
    wb = session.workbench
    for name, call in [
        ("grow_stable", lambda: grow_stable(session, case, fp)),
        ("compute_intersections", lambda: session.compute_intersections([fp], infer_iterates=False)),
        ("trim", lambda: session.trim_stable_manifolds(fp)),
        ("create_bridges", lambda: session.create_bridges(fp)),
        ("infer_iterate_table", lambda: session.infer_iterate_table()),
    ]:
        gc.collect()
        with rec.stage(name):
            call()
        rec.checkpoint()
    rec.metric("points_stable", total_points(session, "stable"))
    rec.metric("segments", len(wb.Tangle._seg_lookup))
    rec.metric("crossings", len(wb.intersection_registry))
    rec.metric("bridges", len(wb.bridges))


def run(config: dict, rec) -> None:
    case, kind = config["case"], config["kind"]
    if kind == "verify":
        return verify(config, rec)
    gpu = config["backend"] == "gpu"
    session, fp = build(case, gpu, config.get("min_batch_points", 64))
    rec.metric("area_cutoff", session.workbench._man_machine.area_cutoff)
    rec.metric("lambda_u", float(np.abs(np.ravel(fp.unstable_eigenvalues[0])[0])))
    if kind == "profile":
        profiler = cProfile.Profile()
        grow_unstable(session, fp, config["steps"], rec, profiler)
        record_profile(profiler, rec)
        return
    grow_unstable(session, fp, config["steps"], rec)
    if kind == "timed":
        downstream(session, case, fp, rec)


# --------------------------------------------------------------------------- #
# Profiling
# --------------------------------------------------------------------------- #
def category(filename: str, function: str) -> str:
    """Which part of the growth machinery a profiled function belongs to."""
    if "cupy" in filename or "henon.py" in filename or function in {
        "map_batch", "map_inv_batch", "_apply_batch", "gpu_fn", "_detect_batchable"
    }:
        return "map"
    if "ManifoldMachine.py" in filename and any(
        word in function for word in ("refine", "curvature", "area", "fit", "chord", "frame")
    ):
        return "refinement"
    if any(name in filename for name in ("Point.py", "BaseManifold.py")) or function in {
        "_collect", "_get_iterate", "_get_preiterate", "_cache_preiterate"
    } or "insert" in function:
        return "linked_list"
    if "Tangle.py" in filename or "rtree" in filename or "Intersection" in filename:
        return "intersections"
    return "other"


def record_profile(profiler: cProfile.Profile, rec) -> None:
    stats = pstats.Stats(profiler).stats
    rows = [
        {"function": func, "file": file, "line": line, "ncalls": nc, "primitive_calls": cc,
         "tottime": tt, "cumtime": ct, "category": category(file, func)}
        for (file, line, func), (cc, nc, tt, ct, _callers) in stats.items()
    ]
    rows.sort(key=lambda r: r["tottime"], reverse=True)
    rec.metric("profile_top40", rows[:40])
    totals: dict[str, float] = {}
    for row in rows:
        totals[row["category"]] = totals.get(row["category"], 0.0) + row["tottime"]
    rec.metric("profile_categories_s", totals)
    rec.metric("profile_total_s", sum(totals.values()))


# --------------------------------------------------------------------------- #
# CPU / GPU agreement
# --------------------------------------------------------------------------- #
def unstable_arrays(session) -> dict:
    return {
        key: manifold.get_point_array()
        for key, manifold in session.workbench.manifolds.items()
        if key[1] == "unstable"
    }


def verify(config: dict, rec) -> None:
    """Grow a CPU and a GPU session in lockstep and compare their unstable manifolds."""
    cpu, fp_cpu = build(config["case"], gpu=False)
    gpu, fp_gpu = build(config["case"], gpu=True)
    for _ in range(config["steps"]):
        cpu.grow_n_times(fp_cpu, "unstable", num_iterations=1)
        gpu.grow_n_times(fp_gpu, "unstable", num_iterations=1)
        a, b = unstable_arrays(cpu), unstable_arrays(gpu)
        rec.append("points_cpu", sum(len(x) for x in a.values()))
        rec.append("points_gpu", sum(len(x) for x in b.values()))
        # Keys hold different FixedPoint objects: pair the branches by position.
        same_shape = all(x.shape == y.shape for x, y in zip(a.values(), b.values()))
        deviation = max(
            float(np.max(np.abs(x - y), initial=0.0)) for x, y in zip(a.values(), b.values())
        ) if same_shape else None
        rec.append("max_abs_deviation", deviation)
        rec.append("max_abs_coord", max_abs_coord(cpu))
        rec.checkpoint()
    rec.metric("final_counts_equal", rec.series["points_cpu"][-1] == rec.series["points_gpu"][-1])
    rec.metric("final_max_abs_deviation", rec.series["max_abs_deviation"][-1])
    first = next((i for i, d in enumerate(rec.series["max_abs_deviation"]) if d is None or d > 0), None)
    rec.metric("first_divergent_step", None if first is None else first + 1)


if __name__ == "__main__":
    harness.main(globals())
