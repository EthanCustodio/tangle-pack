"""Re-verification after the open-issues round (WP12): does every case run end to
end, and does the class count plateau?

The question of the round's last work package. Every case is a frozen
``tests/cases.py`` recipe (``e11_blast_depth.build``) and every blast is
followed by one ``session.topology()`` (D1: the topology is recomputed from
scratch), snapshotted by ``e11_blast_depth.record_topology``.

Studies (``study``)
-------------------
* ``series`` -- the D1 loop: blast the case's zone ``max_blasts`` times, one
  ``blast_zone`` call per step, ``num_iterations`` iterations per call
  (``diagnose=False``; ``max_points`` is the per-call budget, C6), at the
  blasting cutoff ``area_cutoff`` (every case is BUILT at its recipe's cutoff).
  Cases and depths: k28 and k32 to 13 blasts (the depths WP8b measured), p3 to
  20, nested_outer and nested_inner to 12, inversion (k=4) to 10, heteroclinic
  (its inner zone) to 8, each at 1e-7 and at the coarse 1e-5 (NUM-5), one
  iteration per call; plus the NUM-6 / NUM-7 cells at 1e-7 with TWO iterations
  per call (k28, k32, nested_outer to 6 calls, p3 to 10), to compare call ``n``
  with single blast ``2n``.
* ``depth`` -- growth instead of blasting: p3 at 11..18 unstable steps (9
  stable) and inversion (k=4) at 8..12 unstable steps (11 stable, 1e-5);
  the topology of the final state. Step 12 of the inversion arm escapes into
  the library's per-step point budget, recorded as ``stop_reason =
  "point_budget"`` (a datum, not a failure).

Not rerun (numerics untouched since the second edition): E1-E6, E8-E10, the
E7 escape and tangency studies and the full 58-cell E11 grid (separations
other than each recipe's, cutoffs 1e-6 / 1e-8, repeats).

A series stops early on its wall budget (``_budget`` = timeout - 120 s, the
next blast predicted at ~2x the last) or its memory budget (0.9 x
``RSS_CAP_GB``, peak RSS x the last growth factor), recorded as
``stop_reason`` ``budget`` / ``memory_budget``; ``max_blasts`` otherwise.

Record schema: e11's series (``n_classes`` = the seed table, ``n_closed`` =
the alphabet closed under the word map, ``n_flagged``, ``n_holes``,
``n_trellis_crossings`` = N, ``crossings`` = the registry, ``is_reliable``,
``n_ambiguous`` / ``n_disagree`` / ``n_straddle`` / ``n_unresolved``, the
lobe ratios, ``blast_s``, ``rss_mb``) plus ``blast_stop`` (each call's
``BlastResult.stop_reason``) and ``points_created``.
"""

from __future__ import annotations

import sys
import time

import harness
import e11_blast_depth as e11
from maps import curve_arrays, total_points

sys.path.insert(0, str(harness.REPO / "tests"))

EXPERIMENT = "e12_reverify"
TIMEOUT = 900
WORKERS = 8
RSS_CAP_GB = 12
MAX_POINTS = 4_000_000

#: Each recipe's blast separation and the series depth (WP8b's for k28 / k32).
SEPARATION = {"k28": 1e-5, "k32": 1e-5, "p3": 1e-5, "nested_outer": 1e-4,
              "nested_inner": 1e-4, "inversion": 1e-5, "heteroclinic": 1e-4}
DEPTH = {"k28": 13, "k32": 13, "p3": 20, "nested_outer": 12, "nested_inner": 12,
         "inversion": 10, "heteroclinic": 8}
TWO_ITERATIONS = {"k28": 6, "k32": 6, "nested_outer": 6, "p3": 10}


def configs(quick: bool) -> list[dict]:
    budgets = {"_timeout": TIMEOUT, "_budget": TIMEOUT - 120, "_rss_budget_gb": 0.9 * RSS_CAP_GB}
    if quick:
        return [
            {"study": "series", "case": "k28", "max_blasts": 3, "area_cutoff": 1e-7,
             "num_iterations": 1},
            {"study": "series", "case": "p3", "max_blasts": 2, "area_cutoff": 1e-7,
             "num_iterations": 1},
            {"study": "depth", "case": "p3", "unstable_steps": 11},
        ]
    series = [
        {"study": "series", "case": case, "max_blasts": depth, "area_cutoff": cutoff,
         "num_iterations": 1, **budgets}
        for case, depth in DEPTH.items() for cutoff in (1e-7, 1e-5)
    ] + [
        {"study": "series", "case": case, "max_blasts": calls, "area_cutoff": 1e-7,
         "num_iterations": 2, **budgets}
        for case, calls in TWO_ITERATIONS.items()
    ]
    depth = [{"study": "depth", "case": "p3", "unstable_steps": n, **budgets}
             for n in range(11, 19)]
    depth += [{"study": "depth", "case": "inversion", "unstable_steps": n, **budgets}
              for n in range(8, 13)]
    # The heavy period-1 series first, so the second-long cells fill in behind them.
    return sorted(series + depth, key=lambda c: (c["study"] != "series",
                                                 c["case"] not in ("k28", "k32", "nested_outer")))


def run(config: dict, rec: harness.Recorder) -> None:
    for name, value in config.items():
        if not name.startswith("_"):
            rec.metric(name, value)
    (run_series if config["study"] == "series" else run_depth)(config, rec)


# --------------------------------------------------------------------------- #
# series
# --------------------------------------------------------------------------- #
def run_series(config: dict, rec: harness.Recorder) -> None:
    t_start = time.perf_counter()
    with rec.stage("build"):
        session, fixed_points, pips, zone = e11.build(config["case"])
    rec.metric("build_s", rec.stages["build"])
    session.area_cutoff = config["area_cutoff"]
    wb = session.workbench
    known = {id(b) for b in wb.bridges}
    curves, curve_size = curve_arrays(session), total_points(session)

    def append_numerics(new_bridges, blast, blast_s, result) -> None:
        values = {"blast": blast, "blast_s": blast_s,
                  "blast_stop": result.stop_reason if result else None,
                  "points_created": result.points_created if result else 0,
                  **e11.numerics_snapshot(session, new_bridges),
                  **e11.law_checks(session, new_bridges, curves)}
        for name, value in values.items():
            rec.append(name, value)
        rec.append("rss_mb", rec.peak_rss_mb)

    previous, _ = e11.record_topology(rec, session, fixed_points, pips, None)
    append_numerics([], 0, 0.0, None)
    rec.checkpoint()
    stop_reason, last_cost = "max_blasts", 0.0
    budget = config.get("_budget", float("inf"))
    rss_cap_mb = config.get("_rss_budget_gb", float("inf")) * 1024
    for blast in range(1, config["max_blasts"] + 1):
        if time.perf_counter() - t_start + 2 * last_cost > budget:
            stop_reason = "budget"
            break
        peaks = rec.series["rss_mb"]
        growth = max(peaks[-1] / peaks[-2], 1.0) if len(peaks) > 1 and peaks[-2] else 1.0
        if peaks[-1] * growth > rss_cap_mb:
            stop_reason = "memory_budget"
            break
        t_blast = time.perf_counter()
        result = session.blast_zone(zone, num_iterations=config["num_iterations"],
                                    fixed_point=[zone.fixed_point],
                                    min_separation=SEPARATION[config["case"]],
                                    diagnose=False, max_points=MAX_POINTS)
        blast_s = time.perf_counter() - t_blast
        new_bridges = [b for b in wb.bridges if id(b) not in known]
        known |= {id(b) for b in new_bridges}
        if total_points(session) != curve_size:
            curves, curve_size = curve_arrays(session), total_points(session)
        append_numerics(new_bridges, blast, blast_s, result)
        snapshot, _delta = e11.record_topology(rec, session, fixed_points, pips, previous)
        previous = snapshot or previous
        rec.metric("blasts_done", blast)
        rec.checkpoint()
        last_cost = time.perf_counter() - t_blast
    rec.metric("stop_reason", stop_reason)


# --------------------------------------------------------------------------- #
# depth
# --------------------------------------------------------------------------- #
def run_depth(config: dict, rec: harness.Recorder) -> None:
    import cases
    from tanglepack.examples import henon_cases
    from tanglepack.numerics.budget import PointBudgetExceeded

    steps = config["unstable_steps"]
    with rec.stage("build"):
        try:
            if config["case"] == "p3":
                built = cases.build_period3(unstable_steps=steps)
            else:
                built = henon_cases.build_inversion(unstable_steps=steps)
        except PointBudgetExceeded as exc:
            rec.metric("stop_reason", "point_budget")
            rec.metric("point_budget_error", str(exc)[:300])
            return
    rec.metric("stop_reason", "done")
    session = built.session
    rec.metric("crossings", len(session.workbench.intersection_registry))
    rec.metric("unstable_points", total_points(session, "unstable"))
    snapshot, _ = e11.record_topology(rec, session, list(built.fixed_points), list(built.pips),
                                      None)
    for name in ("n_classes", "n_closed", "n_holes", "is_reliable", "n_flagged",
                 "n_ambiguous", "n_disagree", "n_straddle", "n_unresolved",
                 "n_trellis_crossings", "topology_outcome"):
        rec.metric(name, rec.series[name][-1])


if __name__ == "__main__":
    harness.main(globals())
