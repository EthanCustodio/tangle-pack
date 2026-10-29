"""
Wall time against ``workers``: the forced Duffing map and the scale blast series.

Usage (from the repository root)::

    env/bin/python stress/profile/bench_workers.py                    # everything
    env/bin/python stress/profile/bench_workers.py --workloads scale_1e-7 --workers 1 8
    env/bin/python stress/profile/bench_workers.py --table-only --out DIR

Every (workload, workers, repeat) runs in a FRESH interpreter pinned with
``taskset`` to the first ``workers`` cores of ``--cores`` (default 16-31: the
author's sweeps take the low cores; check ``uptime`` first). The table is the
median of ``--repeats`` (3) runs. Output: ``figures/numerics_profile/<sha>_workers/``,
one ``<workload>_w<N>_r<k>.json`` per run.

Workloads:

* ``duffing_leapfrog_thread`` / ``duffing_leapfrog_process`` -- the leapfrog
  Duffing map (``examples/duffing.py``, batch-capable) at ``area_cutoff = 1e-7``,
  7 steps of each manifold, then the crossings; the map split in a thread or a
  forked process pool;
* ``duffing_ivp_process`` -- the ``solve_ivp`` form (one point per call) at
  1e-5, 5 steps of each manifold, in a forked process pool;
* ``scale_1e-7`` / ``scale_1e-9`` -- ``profile_numerics.py``'s scale case: the
  k = 2.5 zone blasted one iteration at a time (``min_separation = 1e-5``,
  ``diagnose=True``) to a plateau (two blasts with no new point) or 40 blasts,
  no topology between them. A blast's time splits into ``diagnose``
  (``report_resolution``), ``separation`` (the interiors and the guard's
  KD-trees, built and queried) and the rest. ``sweep_serial`` is the part of ``diagnose``
  outside the sweep's pair blocks: the Amdahl ceiling of its threads (read at
  every worker count, meaningful at 1).

Only the main thread is timed; time in the pools counts where they were called.
"""

from __future__ import annotations

import argparse
import builtins
import json
import os
import statistics
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

WORKLOADS = ("duffing_leapfrog_thread", "duffing_leapfrog_process", "duffing_ivp_process",
             "scale_1e-7", "scale_1e-9")
MAX_BLASTS = 40


class Clock:
    """Seconds per label, the outermost call of a label counted once."""

    def __init__(self) -> None:
        self.seconds: dict[str, float] = defaultdict(float)
        self.depth: dict[str, int] = defaultdict(int)

    def wrap(self, fn, label: str):
        """``fn`` timed under ``label``."""

        def timed(*args, **kwargs):
            self.depth[label] += 1
            start = time.perf_counter()
            try:
                return fn(*args, **kwargs)
            finally:
                self.depth[label] -= 1
                if not self.depth[label]:
                    self.seconds[label] += time.perf_counter() - start

        return timed


def install(clock: Clock) -> None:
    """Time the diagnose, the separation and the sweep's pair blocks."""
    import importlib

    from tanglepack.numerics.TangleWorkbench import TangleWorkbench

    blast = importlib.import_module("tanglepack.loom.Blast")
    diagnostics = importlib.import_module("tanglepack.numerics.diagnostics")
    TangleWorkbench.report_resolution = clock.wrap(TangleWorkbench.report_resolution, "diagnose")
    blast._interior_points = clock.wrap(blast._interior_points, "separation")
    blast._Cloud.add = clock.wrap(blast._Cloud.add, "separation")
    blast._Cloud.near = clock.wrap(blast._Cloud.near, "separation")

    # The pair blocks run through ``map`` (one worker) or a thread pool's ``map``,
    # timed until every block's result is in.
    def blocks(fn, items):
        return list(builtins.map(fn, items))

    class TimedPool(diagnostics.ThreadPoolExecutor):
        """A thread pool whose ``map`` counts as the sweep's pair blocks."""

        def map(self, fn, *iterables, **kwargs):
            return clock.wrap(lambda: list(super(TimedPool, self).map(fn, *iterables, **kwargs)),
                              "sweep_blocks")()

    diagnostics.map = clock.wrap(blocks, "sweep_blocks")
    diagnostics.ThreadPoolExecutor = TimedPool


def duffing(form: str, executor: str, workers: int) -> dict:
    """Grow the Duffing tangle and find its crossings."""
    from tanglepack import TangleSession
    from tanglepack.examples.duffing import (
        duffing_map, duffing_map_inverse, duffing_map_ivp, duffing_map_ivp_inverse, saddle_guess)

    maps = ((duffing_map(), duffing_map_inverse()) if form == "leapfrog"
            else (duffing_map_ivp(), duffing_map_ivp_inverse()))
    cutoff, steps = (1e-7, 7) if form == "leapfrog" else (1e-5, 5)
    session = TangleSession(*maps, area_cutoff=cutoff, workers=workers, map_executor=executor)
    fp = session.construct_fixed_point(saddle_guess())
    session.initialize_both_manifolds(fp)
    start = time.perf_counter()
    for stability in ("unstable", "stable"):
        session.grow_n_times(fp, stability, num_iterations=steps)
    grown = time.perf_counter()
    registry = session.compute_intersections([fp])
    done = time.perf_counter()
    return {"growth_s": grown - start, "intersections_s": done - grown,
            "points": len(session.point_store), "crossings": len(registry.all_ids())}


def scale(cutoff: float, workers: int, clock: Clock) -> dict:
    """The scale case's blast series, no topology between the blasts."""
    import profile_numerics

    profile_numerics.set_workers(workers)
    start = time.perf_counter()
    session, fp, zone = profile_numerics.scale_build(cutoff)
    built = time.perf_counter()
    clock.seconds.clear()
    blasts, flat = [], 0
    for _ in range(MAX_BLASTS):
        t = time.perf_counter()
        result = session.blast_zone(zone, num_iterations=1, fixed_point=[fp],
                                    min_separation=profile_numerics.SCALE_MIN_SEPARATION,
                                    diagnose=True)
        blasts.append({"s": time.perf_counter() - t, "points_created": result.points_created})
        flat = flat + 1 if not result.points_created else 0
        if flat >= 2:
            break
    blast_s = sum(b["s"] for b in blasts)
    diagnose, separation = clock.seconds["diagnose"], clock.seconds["separation"]
    return {"build_s": built - start, "blast_s": blast_s, "diagnose_s": diagnose,
            "separation_s": separation, "rest_s": blast_s - diagnose - separation,
            "sweep_serial_s": diagnose - clock.seconds["sweep_blocks"],
            "blasts": blasts, "points": len(session.point_store)}


def run(workload: str, workers: int) -> dict:
    """One run in this interpreter."""
    import logging

    logging.basicConfig(level=logging.ERROR)
    clock = Clock()
    install(clock)
    record = {"workload": workload, "workers": workers, "affinity": sorted(os.sched_getaffinity(0))}
    if workload.startswith("duffing"):
        _, form, executor = workload.split("_")
        record |= duffing(form, executor, workers)
    else:
        record |= scale(float(workload.split("_")[1]), workers, clock)
    return record


def git_sha() -> str:
    """The checkout's short sha."""
    return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=REPO,
                                   text=True).strip()


def cores(spec: str) -> list[int]:
    """``"16-31"`` or ``"16,18,20"`` as a list of cores."""
    out: list[int] = []
    for part in spec.split(","):
        lo, _, hi = part.partition("-")
        out += range(int(lo), int(hi or lo) + 1)
    return out


#: Columns of the table per workload kind.
COLUMNS = {"duffing": ("growth_s", "intersections_s"),
           "scale": ("blast_s", "diagnose_s", "separation_s", "rest_s", "sweep_serial_s")}


def table(out: Path) -> str:
    """The median of every column per (workload, workers), and the speedups against 1."""
    runs: dict[tuple, list[dict]] = defaultdict(list)
    for path in sorted(out.glob("*_w*_r*.json")):
        record = json.loads(path.read_text())
        runs[record["workload"], record["workers"]].append(record)
    lines = []
    for workload in WORKLOADS:
        counts = sorted(w for name, w in runs if name == workload)
        if not counts:
            continue
        columns = COLUMNS[workload.split("_")[0]]
        median = {w: {c: statistics.median(r[c] for r in runs[workload, w]) for c in columns}
                  for w in counts}
        lines.append(f"{workload}  (median of {len(runs[workload, counts[0]])}; "
                     f"{runs[workload, counts[0]][0]['points']} points)")
        lines.append("  workers  " + "  ".join(f"{c:>16}" for c in columns))
        for w in counts:
            cells = []
            for c in columns:
                base = median[counts[0]][c]
                speedup = f" ({base / median[w][c]:.2f}x)" if median[w][c] > 0 else ""
                cells.append(f"{median[w][c]:8.2f}{speedup:>8}")
            lines.append(f"  {w:7d}  " + "  ".join(f"{cell:>16}" for cell in cells))
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--workloads", nargs="+", default=list(WORKLOADS), choices=WORKLOADS)
    parser.add_argument("--workers", nargs="+", type=int, default=[1, 4, 8, 16])
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--cores", default="16-31", help="the free cores to pin runs to")
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--table-only", action="store_true")
    parser.add_argument("--run", nargs=3, metavar=("WORKLOAD", "WORKERS", "PATH"),
                        help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.run:
        workload, workers, path = args.run
        Path(path).write_text(json.dumps(run(workload, int(workers)), indent=1))
        return

    out = args.out or REPO / "figures" / "numerics_profile" / f"{git_sha()}_workers"
    out.mkdir(parents=True, exist_ok=True)
    free = cores(args.cores)
    if not args.table_only:
        for workload in args.workloads:
            for repeat in range(args.repeats):
                for workers in args.workers:
                    if workers > len(free):
                        raise SystemExit(f"{workers} workers need more than the {len(free)} free cores")
                    pinned = ",".join(map(str, free[:workers]))
                    path = out / f"{workload}_w{workers}_r{repeat}.json"
                    t = time.perf_counter()
                    subprocess.run(["taskset", "-c", pinned, sys.executable, __file__, "--run",
                                    workload, str(workers), str(path)],
                                   cwd=REPO, check=True, env={**os.environ, "MPLBACKEND": "Agg"})
                    print(f"{workload} workers {workers} repeat {repeat}: "
                          f"{time.perf_counter() - t:.0f} s", flush=True)
    text = table(out)
    (out / "table.txt").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
