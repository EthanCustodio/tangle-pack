"""E3: how the numerics scale with point count, on the CPU.

Two config kinds:

* ``growth`` -- one period-1 saddle per (k, area_cutoff), the unstable manifold
  grown one step at a time until it holds more than ``MAX_POINTS`` points, a
  step takes more than ``MAX_STEP_S``, ``MAX_STEPS`` steps are done, or the
  next step is predicted (points x last growth factor) to exceed
  ``MAX_PREDICTED_POINTS`` -- the guard against a single step that would only
  measure the RSS watchdog. Per step: sizes, time, memory, max|coord| and
  arclength, and the point growth factor to compare with ``lambda_u``.
* ``stages`` -- grow a fixed number of unstable steps (chosen from the growth
  curves so the unstable point count spans ~1e3 .. ~1e6), the stable side to
  turnaround, then time ``compute_intersections``, ``trim + create_bridges``
  and ``infer_iterate_table`` separately as functions of the point count.

Note: past ~10 steps the unstable arm of every period-1 Hénon saddle here
escapes to infinity (max|coord| ~1e3 at 1e4 points, ~1e6 at 4e5 points), and
the growth factor jumps far above lambda_u: most of the points are the
escaping tip, not the tangle.
"""

from __future__ import annotations

import gc
import time

import numpy as np
import psutil

import harness
from maps import max_abs_coord, saddle_session, total_points

EXPERIMENT = "e03_scaling"
TIMEOUT = 900
WORKERS = 1  # timings: never share the machine with another worker
# CPU only: stages skip the GPU sync, which would import CuPy (~250 MB RSS).

MAX_POINTS = 2_000_000
MAX_STEP_S = 90.0
MAX_STEPS = 25
MAX_PREDICTED_POINTS = 10_000_000

#: (k, area_cutoff) -> unstable step counts for the stage configs, read off
#: probe growth curves (points: k=2.8@1e-7 1.1e3 3.6e3 1.7e4 4.1e5;
#: k=2.8@1e-5 7.9e2 3.7e3 8.6e4; k=4 1.4e3 6.9e3 8.9e4;
#: k=6 1.2e3 4.5e3 3.6e4 2.2e6; k=10 2.2e3 1.3e4 5.2e5).
STAGE_STEPS = {
    (2.8, 1e-7): [9, 10, 11, 12],
    (2.8, 1e-5): [10, 11, 12],
    (4.0, 1e-7): [9, 10, 11],
    (6.0, 1e-7): [8, 9, 10, 11],
    (10.0, 1e-7): [8, 9, 10],
}


def configs(quick: bool) -> list[dict]:
    if quick:
        return [
            {"kind": "growth", "k": 2.8, "area_cutoff": 1e-7, "max_points": 20_000},
            {"kind": "growth", "k": 10.0, "area_cutoff": 1e-7, "max_points": 20_000},
            {"kind": "stages", "k": 2.8, "area_cutoff": 1e-7, "steps": 9},
            {"kind": "stages", "k": 10.0, "area_cutoff": 1e-7, "steps": 8},
        ]
    out = [{"kind": "growth", "k": k, "area_cutoff": cutoff, "max_points": MAX_POINTS}
           for k, cutoff in STAGE_STEPS]
    out += [{"kind": "stages", "k": k, "area_cutoff": cutoff, "steps": n}
            for (k, cutoff), steps in STAGE_STEPS.items() for n in steps]
    return out


def rss_mb() -> float:
    return psutil.Process().memory_info().rss / 2**20


def unstable_arclength(session) -> float:
    return sum(
        float(np.linalg.norm(np.diff(m.get_point_array(), axis=0), axis=1).sum())
        for key, m in session.workbench.manifolds.items()
        if key[1] == "unstable"
    )


def run(config: dict, rec) -> None:
    gc.collect()
    baseline = rss_mb()
    session, fp = saddle_session(config["k"], area_cutoff=config["area_cutoff"])
    lambda_u = float(np.abs(np.ravel(fp.unstable_eigenvalues[0])[0]))
    rec.metric("lambda_u", lambda_u)
    rec.metric("baseline_rss_mb", baseline)
    if config["kind"] == "growth":
        growth(session, fp, config, rec, baseline)
    else:
        stages(session, fp, config, rec, baseline)


def growth(session, fp, config: dict, rec, baseline: float) -> None:
    points = total_points(session, "unstable")
    factor = 1.0
    for step in range(1, MAX_STEPS + 1):
        if points * factor > MAX_PREDICTED_POINTS:
            rec.metric("stop_reason", f"predicted {points * factor:.3g} points next step")
            return
        gc.collect()
        t0 = time.perf_counter()
        session.grow_n_times(fp, "unstable", num_iterations=1)
        step_s = time.perf_counter() - t0
        previous, points = points, total_points(session, "unstable")
        factor = points / previous
        current = rss_mb()
        rec.append("step", step)
        rec.append("step_s", step_s)
        rec.append("points_unstable", points)
        rec.append("points_stable", total_points(session, "stable"))
        rec.append("growth_factor", factor)
        rec.append("rss_mb", current)
        rec.append("peak_rss_mb", rec.peak_rss_mb)
        rec.append("bytes_per_point", (current - baseline) * 2**20 / points)
        rec.append("max_abs_coord", max_abs_coord(session))
        rec.append("arclength_unstable", unstable_arclength(session))
        rec.checkpoint()
        if points > config["max_points"]:
            rec.metric("stop_reason", "max points")
            return
        if step_s > MAX_STEP_S:
            rec.metric("stop_reason", "step time")
            return
    rec.metric("stop_reason", "max steps")


def stages(session, fp, config: dict, rec, baseline: float) -> None:
    wb = session.workbench
    with rec.stage("grow_unstable", sync=False):
        session.grow_n_times(fp, "unstable", num_iterations=config["steps"])
    with rec.stage("grow_stable", sync=False):
        session.grow_until_turnaround(fp, "stable", max_iterations=15)
    points = total_points(session)
    rec.metric("points_unstable", total_points(session, "unstable"))
    rec.metric("points_stable", total_points(session, "stable"))
    rec.metric("max_abs_coord", max_abs_coord(session))
    rec.metric("rss_after_growth_mb", rss_mb())
    rec.metric("bytes_per_point", (rss_mb() - baseline) * 2**20 / points)
    rec.checkpoint()

    gc.collect()
    with rec.stage("compute_intersections", sync=False):
        session.compute_intersections([fp], infer_iterates=False)
    rec.metric("segments", len(wb.Tangle._seg_lookup))
    rec.metric("crossings", len(wb.intersection_registry))
    rec.metric("rss_after_intersections_mb", rss_mb())
    rec.checkpoint()

    gc.collect()
    with rec.stage("trim_and_bridges", sync=False):
        session.trim_stable_manifolds(fp)
        session.create_bridges(fp)
    rec.metric("bridges", len(wb.bridges))
    rec.checkpoint()

    gc.collect()
    with rec.stage("infer_iterate_table", sync=False):
        rec.metric("iterates_inferred", session.infer_iterate_table())
    rec.metric("rss_final_mb", rss_mb())


if __name__ == "__main__":
    harness.main(globals())
