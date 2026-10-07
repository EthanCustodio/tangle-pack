"""Loading and derived quantities shared by make_figures.py and make_numbers.py.

Every number a figure shows and a macro quotes is computed here, once, so a
plot and the sentence citing it can never disagree. Nothing here imports
tanglepack or touches the GPU: it only reads ``results/*.jsonl``.

Records follow ``harness.py``: ``{experiment, index, config, outcome, error,
wall_s, stages, metrics, series, warnings, peak_rss_mb, ...}``. A record that
crashed or was killed may lack ``metrics``/``series``/``stages``; ``load``
fills them with empty dicts so callers only ever filter, never guard.
"""

from __future__ import annotations

import json
import math
import warnings
from collections import defaultdict
from functools import cache
from pathlib import Path
from typing import Any, Iterable, Optional

import numpy as np

RESULTS = Path(__file__).resolve().parent.parent / "results"

EXPERIMENTS = {
    "e01": "e01_map_throughput",
    "e02": "e02_gpu_end_to_end",
    "e03": "e03_scaling",
    "e04": "e04_area_cutoff",
    "e05": "e05_invariants",
    "e06": "e06_parameter_sweep",
    "e07": "e07_breaking",
    "e08": "e08_fixed_point_solver",
    "e09": "e09_intersections",
    "e10": "e10_determinism",
}
OUTCOMES = ("ok", "exception", "timeout", "oom", "invariant_assert", "crash")

#: Set by the scripts' ``--quick`` flag before the first ``load``.
QUICK = False
#: experiment prefix -> ("full" | "quick" | "missing", path, record count)
SOURCES: dict[str, tuple[str, str, int]] = {}


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #
def _read(path: Path) -> list[dict]:
    if not path.exists():
        return []
    records = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        try:
            records.append(json.loads(line))
        except json.JSONDecodeError:
            print(f"  [load] {path.name}: skipped a partial line (still being written?)")
    for r in records:
        for key in ("config", "metrics", "series", "stages", "warnings"):
            if not isinstance(r.get(key), dict):
                r[key] = {}
        r.setdefault("outcome", "crash")
    return sorted(records, key=lambda r: r.get("index", 0))


@cache
def load(prefix: str) -> list[dict]:
    """Every record of one experiment: full results unless QUICK, else quick.

    Without ``QUICK`` an experiment whose full file is missing or empty falls
    back to its quick file, with a printed warning.
    """
    name = EXPERIMENTS[prefix]
    full, quick = RESULTS / f"{name}.jsonl", RESULTS / f"{name}_quick.jsonl"
    if not QUICK:
        records = _read(full)
        if records:
            SOURCES[prefix] = ("full", full.name, len(records))
            return records
        print(f"  [load] WARNING {prefix}: {full.name} missing or empty, using {quick.name}")
    records = _read(quick)
    SOURCES[prefix] = ("quick" if records else "missing", quick.name, len(records))
    if not records:
        print(f"  [load] WARNING {prefix}: no records at all")
    return records


def ok(prefix: str, **config) -> list[dict]:
    """The ``ok`` records of an experiment whose config matches ``config``."""
    return [r for r in load(prefix) if r["outcome"] == "ok" and matches(r, config)]


def matches(record: dict, config: dict) -> bool:
    return all(record["config"].get(k) == v for k, v in config.items())


def outcome_counts(prefix: str) -> dict[str, int]:
    counts: dict[str, int] = defaultdict(int)
    for r in load(prefix):
        counts[r["outcome"]] += 1
    return dict(counts)


# --------------------------------------------------------------------------- #
# Small numerics
# --------------------------------------------------------------------------- #
def num(value: Any) -> float:
    """A JSON value as a float: None, 'nan' and non-numbers become nan."""
    if value is None or isinstance(value, bool):
        return math.nan if value is None else float(value)
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def arr(values: Iterable) -> np.ndarray:
    return np.array([num(v) for v in values], dtype=float)


def fit_loglog(x: Iterable, y: Iterable, x_min: float = 0.0) -> Optional[tuple[float, float]]:
    """Least-squares ``log y = slope log x + c`` over finite positive points >= x_min."""
    x, y = arr(x), arr(y)
    keep = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0) & (x >= x_min)
    if keep.sum() < 2 or len(np.unique(x[keep])) < 2:
        return None
    slope, intercept = np.polyfit(np.log10(x[keep]), np.log10(y[keep]), 1)
    return float(slope), float(intercept)


def finite(values: Iterable) -> np.ndarray:
    a = arr(values)
    return a[np.isfinite(a)]


def nanmax(values: Iterable) -> Optional[float]:
    a = finite(values)
    return float(a.max()) if len(a) else None


def nanmin(values: Iterable) -> Optional[float]:
    a = finite(values)
    return float(a.min()) if len(a) else None


def nanmedian(values: Iterable) -> Optional[float]:
    a = finite(values)
    return float(np.median(a)) if len(a) else None


def padded(rows: list[list]) -> np.ndarray:
    """Ragged per-repeat series as a (repeats, steps) array padded with nan."""
    width = max((len(r) for r in rows), default=0)
    out = np.full((len(rows), width), np.nan)
    for i, row in enumerate(rows):
        out[i, : len(row)] = arr(row)
    return out


# --------------------------------------------------------------------------- #
# E1 map throughput
# --------------------------------------------------------------------------- #
E1_BACKENDS = ("numpy", "cupy_e2e", "cupy_compute", "cupy_transfer", "map_batch_cpu", "map_batch_gpu")


def e01_sweeps() -> dict[str, dict[str, np.ndarray]]:
    """backend -> {n, median_s, q1_s, q3_s, points_per_s} (records that got at least one size)."""
    out = {}
    for r in load("e01"):
        if r["config"].get("kind") != "sweep" or not r["series"].get("n"):
            continue
        s = r["series"]
        out[r["config"]["backend"]] = {k: arr(s.get(k, [])) for k in ("n", "median_s", "q1_s", "q3_s", "points_per_s")}
    return out


def e01_crossovers() -> dict[str, dict]:
    """pair -> {n (None = no crossover), n_max, note, outcome} from the bisection configs."""
    out = {}
    for r in load("e01"):
        if r["config"].get("kind") != "crossover":
            continue
        m = r["metrics"]
        out[r["config"]["pair"]] = {
            "n": m.get("crossover_n"), "n_max": r["config"].get("hi"),
            "note": m.get("crossover_note"), "bracket": m.get("crossover_bracket"),
            "outcome": r["outcome"], "measured": "crossover_n" in m,
        }
    return out


def e01_sweep_crossover(cpu: dict, gpu: dict) -> Optional[float]:
    """The N where the two sweep curves cross (log-log interpolation), if they do."""
    common = np.intersect1d(cpu["n"], gpu["n"])
    if len(common) < 2:
        return None
    ratio = np.array([gpu["median_s"][gpu["n"] == n][0] / cpu["median_s"][cpu["n"] == n][0] for n in common])
    for i in range(len(common) - 1):
        a, b = np.log10(ratio[i]), np.log10(ratio[i + 1])
        if a >= 0 > b:
            t = a / (a - b)
            return float(10 ** (np.log10(common[i]) + t * np.log10(common[i + 1] / common[i])))
    return float(common[0]) if ratio[0] < 1 else None


def e01_split(sweeps: dict) -> Optional[dict[str, np.ndarray]]:
    """GPU end-to-end time split into transfer, compute and the rest, as fractions per N."""
    if not all(b in sweeps for b in ("cupy_e2e", "cupy_compute", "cupy_transfer")):
        return None
    e2e, comp, tran = (sweeps[b] for b in ("cupy_e2e", "cupy_compute", "cupy_transfer"))
    n = np.intersect1d(np.intersect1d(e2e["n"], comp["n"]), tran["n"])
    if not len(n):
        return None
    pick = lambda d: np.array([d["median_s"][d["n"] == x][0] for x in n])  # noqa: E731
    total, c, t = pick(e2e), pick(comp), pick(tran)
    rest = np.maximum(total - c - t, 0)
    norm = np.maximum(total, c + t)
    return {"n": n, "transfer": t / norm, "compute": c / norm, "other": rest / norm, "e2e_s": total}


# --------------------------------------------------------------------------- #
# E2 end-to-end GPU
# --------------------------------------------------------------------------- #
CASES = ("k28", "p3", "k10")
E2_STAGES = {  # display stage -> recorded stages summed
    "grow unstable": ("grow_unstable",),
    "grow stable": ("grow_stable",),
    "intersections": ("compute_intersections",),
    "trim + bridges": ("trim", "create_bridges"),
    "iterate table": ("infer_iterate_table",),
}


def e02_timed() -> dict[str, dict[str, list[dict]]]:
    out: dict = defaultdict(lambda: defaultdict(list))
    for r in ok("e02", kind="timed"):
        out[r["config"]["case"]][r["config"]["backend"]].append(r)
    return out


def e02_step_curves(timed: dict) -> dict:
    """case -> backend -> {points, median, q1, q3} per growth step over the repeats."""
    out: dict = defaultdict(dict)
    for case, by_backend in timed.items():
        for backend, records in by_backend.items():
            t = padded([r["series"].get("step_s", []) for r in records])
            p = padded([r["series"].get("points_unstable", []) for r in records])
            if not t.size:
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)  # all-nan columns of short repeats
                out[case][backend] = {
                    "points": np.nanmedian(p, axis=0), "median": np.nanmedian(t, axis=0),
                    "q1": np.nanpercentile(t, 25, axis=0), "q3": np.nanpercentile(t, 75, axis=0),
                    "repeats": len(records),
                }
    return out


def _stage_total(record: dict, names: Iterable[str]) -> float:
    values = [num(record["stages"].get(n)) for n in names]
    return math.nan if any(map(math.isnan, values)) else sum(values)


def e02_stage_speedups(timed: dict) -> dict[str, dict[str, Optional[float]]]:
    """case -> display stage (and 'total') -> median CPU time / median GPU time."""
    out = {}
    for case, by_backend in timed.items():
        if not by_backend.get("cpu") or not by_backend.get("gpu"):
            continue
        row = {}
        for stage, names in [*E2_STAGES.items(), ("total", sum(E2_STAGES.values(), ()))]:
            cpu = nanmedian(_stage_total(r, names) for r in by_backend["cpu"])
            gpu = nanmedian(_stage_total(r, names) for r in by_backend["gpu"])
            row[stage] = cpu / gpu if cpu and gpu else None
        out[case] = row
    return out


def e02_total_speedups(timed: dict) -> dict[str, list[float]]:
    """case -> per-repeat end-to-end speed-up, CPU total / GPU total paired by repeat."""
    names = sum(E2_STAGES.values(), ())
    out = {}
    for case, by_backend in timed.items():
        gpu = {r["config"].get("repeat"): _stage_total(r, names) for r in by_backend.get("gpu", [])}
        ratios = [_stage_total(r, names) / gpu[r["config"].get("repeat")]
                  for r in by_backend.get("cpu", []) if r["config"].get("repeat") in gpu]
        out[case] = [x for x in ratios if np.isfinite(x)]
    return out


def e02_best_step_speedup(curves: dict, min_cpu_s: float = 0.01) -> dict[str, Optional[dict]]:
    """case -> the largest per-step CPU/GPU median ratio over steps taking >= min_cpu_s on the CPU."""
    out = {}
    for case, by_backend in curves.items():
        cpu, gpu = by_backend.get("cpu"), by_backend.get("gpu")
        if not cpu or not gpu:
            continue
        m = min(len(cpu["median"]), len(gpu["median"]))
        ratio = cpu["median"][:m] / gpu["median"][:m]
        keep = np.isfinite(ratio) & (cpu["median"][:m] >= min_cpu_s)
        if not keep.any():
            out[case] = None
            continue
        i = int(np.argmax(np.where(keep, ratio, -np.inf)))
        out[case] = {"speedup": float(ratio[i]), "step": i + 1, "points": float(cpu["points"][i]),
                     "cpu_s": float(cpu["median"][i])}
    return out


PROFILE_CATEGORIES = ("map", "refinement", "linked_list", "intersections", "other")


def e02_profiles() -> dict[tuple[str, str], dict]:
    """(case, backend) -> category fractions of the profiled growth, f_map and the Amdahl bound."""
    out = {}
    for r in ok("e02", kind="profile"):
        cats = {k: num(v) for k, v in r["metrics"].get("profile_categories_s", {}).items()}
        total = sum(v for v in cats.values() if np.isfinite(v))
        if total <= 0:
            continue
        frac = {c: cats.get(c, 0.0) / total for c in PROFILE_CATEGORIES}
        f_map = frac["map"]
        out[(r["config"]["case"], r["config"]["backend"])] = {
            "fractions": frac, "total_s": total, "f_map": f_map,
            "amdahl": 1 / (1 - f_map) if f_map < 1 else math.inf,
        }
    return out


# --------------------------------------------------------------------------- #
# E3 scaling
# --------------------------------------------------------------------------- #
E3_STAGES = ("compute_intersections", "trim_and_bridges", "infer_iterate_table")


def e03_growth() -> dict[tuple[float, float], dict]:
    """(k, area_cutoff) -> growth series as arrays plus lambda_u and the fits.

    Failed runs keep their checkpointed series, so any record with steps counts.
    """
    out = {}
    for r in load("e03"):
        if r["config"].get("kind") != "growth" or not r["series"].get("step"):
            continue
        s = {k: arr(v) for k, v in r["series"].items()}
        baseline = num(r["metrics"].get("baseline_rss_mb"))
        big = s["points_unstable"] >= max(1e3, np.nanmax(s["points_unstable"]) / 100)
        step_fit = fit_loglog(s["points_unstable"], s["step_s"], x_min=1e3)
        mem = None
        if big.sum() >= 2 and len(np.unique(s["points_unstable"][big])) >= 2:
            mem = float(np.polyfit(s["points_unstable"][big], (s["rss_mb"][big] - baseline) * 2**20, 1)[0])
        out[(r["config"]["k"], r["config"]["area_cutoff"])] = {
            **s, "lambda_u": num(r["metrics"].get("lambda_u")), "baseline_rss_mb": baseline,
            "step_fit": step_fit, "bytes_per_point_fit": mem,
            "bytes_per_point_last": float(s["bytes_per_point"][-1]) if len(s.get("bytes_per_point", [])) else None,
            "stop_reason": r["metrics"].get("stop_reason"), "outcome": r["outcome"],
        }
    return out


def e03_stages() -> list[dict]:
    """One row per stages config: k, cutoff, segments, points and the stage times."""
    rows = []
    for r in load("e03"):
        if r["config"].get("kind") != "stages":
            continue
        m = r["metrics"]
        rows.append({
            "k": r["config"]["k"], "cutoff": r["config"]["area_cutoff"], "steps": r["config"]["steps"],
            "outcome": r["outcome"], "segments": num(m.get("segments")),
            "points": num(m.get("points_unstable")) + num(m.get("points_stable")),
            "crossings": num(m.get("crossings")), "bytes_per_point": num(m.get("bytes_per_point")),
            **{s: num(r["stages"].get(s)) for s in ("grow_unstable", "grow_stable", *E3_STAGES)},
        })
    return rows


def e03_stage_fits(rows: list[dict], x_min: float = 1e3) -> dict[str, Optional[tuple[float, float]]]:
    """Stage time against segment count, pooled over every k (segments >= x_min)."""
    return {s: fit_loglog([r["segments"] for r in rows], [r[s] for r in rows], x_min) for s in E3_STAGES}


# --------------------------------------------------------------------------- #
# E4 area cutoff
# --------------------------------------------------------------------------- #
MEASUREMENT_STAGES = {"invariance", "lobes"}


def e04_runs() -> dict[tuple[str, int], list[dict]]:
    """(case, unstable steps) -> records sorted by decreasing cutoff, with derived fields."""
    out: dict = defaultdict(list)
    for r in load("e04"):
        m, c = r["metrics"], r["config"]
        p50 = [num(v) for k, v in m.items() if k.startswith("inv_") and k.endswith("_p50")]
        p99 = [num(v) for k, v in m.items() if k.startswith("inv_") and k.endswith("_p99")]
        out[(c["case"], c["unstable_steps"])].append({
            "cutoff": c["area_cutoff"], "outcome": r["outcome"], "record": r,
            "time": sum(num(v) for k, v in r["stages"].items() if k not in MEASUREMENT_STAGES) if r["stages"] else math.nan,
            "inv_p50": nanmax(p50), "inv_p99": nanmax(p99),
            "points": num(m.get("points_unstable")) + num(m.get("points_stable")),
            "crossings": num(m.get("n_crossings")), "crossing_rows": m.get("crossings") or [],
        })
    return {key: sorted(rows, key=lambda row: -row["cutoff"]) for key, rows in sorted(out.items())}


def e04_convergence(rows: list[dict]) -> list[dict]:
    """Crossing deviation of each run against the finest-cutoff run of the same case and depth.

    Each reference crossing is matched to the nearest crossing of the run on
    the same (unstable, stable) branch pair; a reference crossing with no
    partner on its pair counts as unmatched.
    """
    from scipy.spatial import cKDTree

    usable = [row for row in rows if row["outcome"] == "ok" and row["crossing_rows"]]
    if len(usable) < 2:
        return []
    ref = min(usable, key=lambda row: row["cutoff"])

    def by_pair(crossing_rows):
        out = defaultdict(list)
        for row in crossing_rows:
            out[(row[5], row[6])].append((num(row[1]), num(row[2])))
        return {k: np.array(v) for k, v in out.items()}

    reference = by_pair(ref["crossing_rows"])
    result = []
    for row in usable:
        if row is ref:
            continue
        run = by_pair(row["crossing_rows"])
        devs, unmatched = [], 0
        for pair, points in reference.items():
            if pair not in run:
                unmatched += len(points)
                continue
            devs.extend(cKDTree(run[pair]).query(points)[0])
        devs = np.array(devs)
        result.append({"cutoff": row["cutoff"], "reference_cutoff": ref["cutoff"],
                       "median": float(np.median(devs)) if len(devs) else None,
                       "max": float(devs.max()) if len(devs) else None,
                       "n_reference": int(sum(map(len, reference.values()))),
                       "n_run": len(row["crossing_rows"]), "unmatched": unmatched})
    return result


# --------------------------------------------------------------------------- #
# E5 invariants
# --------------------------------------------------------------------------- #
def e05_cases() -> list[dict]:
    return [r for r in load("e05") if r["config"].get("case") != "gpu_k28"]


def e05_invariance_samples(record: dict) -> dict[str, np.ndarray]:
    """manifold tag ('A:unstable/0.0') -> sampled invariance errors."""
    out = {}
    for key, rows in record["series"].items():
        if key.startswith("inv_") and rows:
            out[key[4:]] = arr(row[1] for row in rows)
    return out


def e05_fp_rows(record: dict) -> list[dict]:
    return [orbit for fp in record["metrics"].get("fp") or [] for orbit in fp.get("orbit", [])]


def e05_lobe_ratios(record: dict) -> np.ndarray:
    rows = record["series"].get("lobes") or []
    return arr(abs(num(r[3]) / num(r[2]) - 1) for r in rows if num(r[2]) not in (0.0,) and np.isfinite(num(r[2])))


# --------------------------------------------------------------------------- #
# E6 parameter sweep
# --------------------------------------------------------------------------- #
E6_STAGES = ("none", "fixed_point", "eigen", "seeded", "grown_unstable", "grown",
             "intersected", "bridges", "trellis_pips")


def e06_family(config: dict) -> str:
    b, which = config.get("b"), config.get("which")
    if b == 1.0:
        return "b = 1, saddle" if which == "saddle" else "b = 1, other point"
    if b == -1.0:
        return "b = -1"
    return f"b = {b:g}"


def e06_rows() -> list[dict]:
    """One row per config: family, k, outcome, stage reached (index into E6_STAGES)."""
    rows = []
    for r in load("e06"):
        stage = r["metrics"].get("stage_reached", "none")
        rows.append({
            "family": e06_family(r["config"]), "k": num(r["config"].get("k")), "outcome": r["outcome"],
            "stage": stage, "stage_index": E6_STAGES.index(stage) if stage in E6_STAGES else 0,
            "error": r.get("error"), "metrics": r["metrics"], "stages": r["stages"], "wall_s": num(r.get("wall_s")),
        })
    return rows


# --------------------------------------------------------------------------- #
# E7 breaking
# --------------------------------------------------------------------------- #
def e07(study: str) -> list[dict]:
    return [r for r in load("e07") if r["config"].get("study") == study]


def e07_blowup_step(fixed: dict, threshold: float = 100.0) -> Optional[int]:
    """First step of the fixed-count escape run whose max|coord| exceeds threshold."""
    s = fixed["series"]
    for step, value in zip(s.get("step", []), s.get("max_abs", [])):
        if num(value) > threshold:
            return int(step)
    return None


# --------------------------------------------------------------------------- #
# E8 solver
# --------------------------------------------------------------------------- #
E8_CLASSES = ("saddle", "lower_period", "non_saddle", "no_convergence", "other")


def e08_class(solve: dict, period: int) -> str:
    """Outcome class of one solve; any converged orbit of smaller minimal period is 'lower_period'."""
    if "minimal_period" in solve and solve["minimal_period"] < period:
        return "lower_period"
    outcome = solve.get("outcome", "")
    if outcome == "saddle":
        return "saddle"
    if outcome in ("elliptic", "not_hyperbolic"):
        return "non_saddle"
    if outcome.startswith("no_convergence"):
        return "no_convergence"
    return "other"


def e08_map(config: dict) -> str:
    return f"k = {num(config['k']):g}"


def e08_table() -> dict[tuple[str, int], dict]:
    """(map, period) -> class counts, converged count, solve times, distinct orbits."""
    out = {}
    for r in load("e08"):
        c = r["config"]
        solves = r["series"].get("solves") or []
        counts = {k: 0 for k in E8_CLASSES}
        for s in solves:
            counts[e08_class(s, c["period"])] += 1
        converged = [s for s in solves if "residual" in s]
        out[(e08_map(c), c["period"])] = {
            "outcome": r["outcome"], "n": len(solves), "counts": counts,
            "converged": len(converged),
            "lower_period": sum(s["minimal_period"] < c["period"] for s in converged),
            "solve_s": arr(s.get("solve_s") for s in solves),
            "distinct": r["metrics"].get("distinct_orbit_list") or {},
            "horseshoe": r["metrics"].get("horseshoe_orbits"),
            "max_residual": num(r["metrics"].get("max_residual")),
        }
    return out


def e08_distinct_by_period(table: dict) -> dict[str, dict[int, int]]:
    """map -> minimal period -> distinct orbits found over every config of that map (any outcome)."""
    seen: dict = defaultdict(lambda: defaultdict(set))
    for (map_name, _), row in table.items():
        for label, orbits in row["distinct"].items():
            period = int(label.split(":")[0][1:])
            seen[map_name][period].update(orbits)
    return {m: {p: len(o) for p, o in sorted(d.items())} for m, d in seen.items()}


HORSESHOE_ORBITS = {1: 2, 2: 1, 3: 2, 4: 3, 5: 6, 6: 9}


# --------------------------------------------------------------------------- #
# E9 intersections
# --------------------------------------------------------------------------- #
E9_MISSED = ("near_parallel", "zigzag", "collision", "missed_unexplained")
E9_EXTRA = ("extra_touch", "extra_unexplained")


def e09_rows() -> list[dict]:
    rows = []
    for r in load("e09"):
        m = r["metrics"]
        rows.append({
            "case": r["config"].get("case"), "steps": r["config"].get("unstable_steps"),
            "cutoff": r["config"].get("area_cutoff"), "outcome": r["outcome"],
            "segments": num(m.get("unstable_segments")) + num(m.get("stable_segments")),
            "pairs": num(m.get("unstable_segments")) * num(m.get("stable_segments")),
            "library_s": num(m.get("library_s")), "brute_s": num(m.get("brute_s")),
            "library": num(m.get("library_crossings")) - num(m.get("anchors")),
            "brute": num(m.get("brute_crossings")), "matched": num(m.get("matched")),
            **{k: num(m.get(k)) for k in (*E9_MISSED, *E9_EXTRA)},
            "uu": m.get("uu_crossings", "missing"),
        })
    return rows


# --------------------------------------------------------------------------- #
# E10 determinism
# --------------------------------------------------------------------------- #
E10_STAGES = ("setup", "grow_unstable", "grow_stable", "intersections", "bridges")
E10_HASHES = ("unstable_hash", "stable_hash", "crossing_hash", "cdist_hash", "bridge_hash")


def e10_by_device() -> dict[str, list[dict]]:
    out: dict = defaultdict(list)
    for r in ok("e10"):
        out[r["config"]["device"]].append(r)
    return out


def cv(values: Iterable) -> Optional[float]:
    a = finite(values)
    return float(a.std(ddof=1) / a.mean()) if len(a) >= 2 and a.mean() > 0 else None


def e10_identity(by_device: dict) -> dict[str, dict[str, dict[str, Optional[bool]]]]:
    """hash -> kind ('exact'/'rounded') -> pair ('CPU-CPU', 'GPU-GPU', 'CPU-GPU') -> all equal?

    None when the pair has fewer than two runs to compare.
    """
    out: dict = {}
    for h in E10_HASHES:
        out[h] = {}
        for kind in ("exact", "rounded"):
            get = lambda rs: [(r["metrics"].get(h) or {}).get(kind) for r in rs]  # noqa: E731
            cpu, gpu = get(by_device.get("cpu", [])), get(by_device.get("gpu", []))
            same = lambda xs: None if len(xs) < 2 else len(set(xs)) == 1  # noqa: E731
            out[h][kind] = {"CPU-CPU": same(cpu), "GPU-GPU": same(gpu),
                            "CPU-GPU": None if not cpu or not gpu else len(set(cpu + gpu)) == 1}
    return out


def e10_crossing_deviation(by_device: dict) -> Optional[float]:
    """Max coordinate difference of the sorted crossing lists, every run against the first CPU run."""
    runs = [np.array(r["metrics"].get("crossing_coords") or [], dtype=float)
            for d in ("cpu", "gpu") for r in by_device.get(d, [])]
    if len(runs) < 2 or not runs[0].size:
        return None
    devs = [float(np.abs(x - runs[0]).max()) for x in runs[1:] if x.shape == runs[0].shape]
    return max(devs) if devs else None
