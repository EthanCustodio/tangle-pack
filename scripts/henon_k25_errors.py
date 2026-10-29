"""
Error analysis of the k=2.5 entropy sweeps: what limits the pip depth, as a function of ``area_cutoff``.

Each (case, cutoff) is rebuilt and blasted one single-iteration step at a time
exactly as ``henon_k25_entropy.py`` does (``min_separation`` 1e-5, until the
frontier is exhausted), but instead of the topology every step records the
numerical errors:

* the registered pip orbit ``f^j(q0)``: distance from its orbit point P, the
  resolution ratio and chord deviation of each crossing;
* the NEXT orbit point ``f^(m+1)(q0)``, the real map applied to ``f^m(q0)``:
  its distance from P, the distance from it to the nearest polyline node and
  segment of this blast's image curves and of the stable branch, the length and
  chord deviation of those two nearest segments, the ratio a crossing there
  would have, and the nearest registered detection off the orbit (resolved or
  not);
* the blast's new image curves: the chord deviation and segment length
  percentiles, the smallest relative cdist step (the cdist precision), their
  nearest node to P;
* the new detections' resolution ratios, the unresolved count, the
  same-stability count of the blast's diagnose;
* points, peak memory, time.

At the end of a sweep the whole state is histogrammed: chord deviation and
segment length per stability, and every detection's resolution ratio.

Records go to ``figures/henon_k25_entropy_post_array_refactor/errors/data/``
(one JSON per case and cutoff, written after every step); ``--plot-only``
draws ``errors.png`` (+ .svg, .fig.pickle) from them, together with the pruned
and unpruned sweep records of the entropy run.

Run with ``PYTHONPATH=src env/bin/python scripts/henon_k25_errors.py --cases p1
--cutoffs 1e-7 1e-8`` (``--max-rss-gb``, ``--max-seconds`` per sweep).
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import sys
import time
import traceback
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

sys.path.insert(0, str(Path(__file__).resolve().parent))

import henon_k25_entropy as entropy  # noqa: E402
from symbolic_figures import FIGURES_DIR, configure_logging, save_panel  # noqa: E402

from tanglepack.examples.henon import henon_map  # noqa: E402
from tanglepack.numerics.geometry import chord_deviation_batch  # noqa: E402

logger = logging.getLogger(__name__)

RUN_DIR = FIGURES_DIR / "henon_k25_entropy_post_array_refactor"
OUT_DIR = RUN_DIR / "errors"
DATA_DIR = OUT_DIR / "data"
#: Machine epsilon: the coordinate floor near P is about this times |P|.
EPS = float(np.finfo(float).eps)
#: log10 bin edges of the end-of-sweep histograms.
LOG_BINS = np.linspace(-20.0, 4.0, 241)


# ── geometry ────────────────────────────────────────────────────────────────


def deviations(points: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Per segment of a polyline: its length and its chord deviation (the refiner's kernel)."""
    if len(points) < 2:
        return np.empty(0), np.empty(0)
    nan = np.full((1, 2), np.nan)
    left = np.vstack([nan, points[:-2]])
    right = np.vstack([points[2:], nan])
    p0, p1 = points[:-1], points[1:]
    return np.hypot(*(p1 - p0).T), chord_deviation_batch(p0, p1, left, right)


def nearest_segment(points: NDArray[np.float64], target: NDArray[np.float64]) -> tuple[float, int]:
    """The distance from ``target`` to a polyline and the index of the nearest segment."""
    if len(points) < 2:
        return math.inf, -1
    p0, d = points[:-1], np.diff(points, axis=0)
    length2 = np.einsum("ij,ij->i", d, d)
    t = np.clip(np.einsum("ij,ij->i", target - p0, d) / np.where(length2 > 0, length2, 1.0), 0.0, 1.0)
    dist = np.hypot(*(p0 + t[:, None] * d - target).T)
    i = int(np.argmin(dist))
    return float(dist[i]), i


def segment_info(curves: list[NDArray[np.float64]], target: NDArray[np.float64]) -> dict:
    """Over several polylines, the segment nearest ``target``: distance, length, deviation, direction."""
    best = {"distance": math.inf}
    for points in curves:
        dist, i = nearest_segment(points, target)
        if i >= 0 and dist < best["distance"]:
            lengths, devs = deviations(points[max(i - 1, 0):i + 3])
            j = i - max(i - 1, 0)
            direction = points[i + 1] - points[i]
            best = {"distance": dist, "length": float(lengths[j]), "deviation": float(devs[j]),
                    "direction": (direction / max(np.hypot(*direction), 1e-300)).tolist(),
                    "node_distance": float(np.min(np.hypot(*(points - target).T)))}
    return best


def percentiles(values: NDArray[np.float64]) -> dict:
    """Min, 1st, 50th, 99th percentile and max of the positive values (empty: {})."""
    values = values[np.isfinite(values) & (values > 0)]
    if not values.size:
        return {}
    q = np.percentile(values, [0, 1, 50, 99, 100])
    return dict(zip(("min", "p1", "p50", "p99", "max"), map(float, q)))


def relative_steps(cdist: NDArray[np.float64]) -> NDArray[np.float64]:
    """The relative cdist step ``dc / c`` between consecutive nodes (positive cdists only)."""
    c = cdist[np.isfinite(cdist)]
    if len(c) < 2:
        return np.empty(0)
    step, base = np.diff(c), c[:-1]
    keep = base > 0
    return np.abs(step[keep]) / base[keep]


# ── per step ────────────────────────────────────────────────────────────────


def orbit_of(session, fp) -> list[int]:
    """The registered pip orbit ``[q0, f(q0), ..., f^m(q0)]`` (registry ids)."""
    registry = session.workbench.intersection_registry
    chosen = session.chosen_strong_pip(fp)
    x = registry.find(chosen) if chosen is not None else None
    orbit = []
    while x is not None:
        orbit.append(x)
        x = registry.iterate_table[x, 1]
    return orbit


def stable_key(ix):
    """The stable branch key of a crossing."""
    return ix.manifold_b_key if ix.manifold_b_key[1] == "stable" else ix.manifold_a_key


def orbit_rows(session, fp, mapping, new_curves: list, detections: dict) -> dict:
    """The registered pip orbit and the next orbit point's resolution, for one fixed point."""
    registry = session.workbench.intersection_registry
    orbit = orbit_of(session, fp)
    rows = []
    for j, iid in enumerate(orbit):
        ix = registry[iid]
        P = fp.coordinates[stable_key(ix)[2]]
        rows.append({"j": j, "id": iid, "distance": float(np.hypot(*(ix.coords - P))),
                     "ratio": ix.resolution_ratio, "deviation": ix.chord_deviation,
                     "stable_cdist": ix.stable_cdist, "unstable_cdist": ix.unstable_cdist})
    out = {"orbit": rows}
    if not orbit:
        return out
    last = registry[orbit[-1]]
    target = np.asarray(mapping(np.asarray(last.coords, dtype=float)), dtype=float)
    key = fp.advance_key(stable_key(last), 1)
    P = fp.coordinates[key[2]]
    stable_curves = [m.get_point_array() for k, m in session.manifolds.items() if k == key]
    unstable = segment_info(new_curves, target)
    stable = segment_info(stable_curves, target)
    nxt = {"distance": float(np.hypot(*(target - P))), "P_norm": float(np.hypot(*P)),
           "unstable": unstable, "stable": stable}
    if "length" in unstable and "length" in stable:
        sin = abs(np.cross(unstable["direction"], stable["direction"]))
        nxt["sin"] = float(sin)
        nxt["ratio_estimate"] = float(min(unstable["length"], stable["length"]) * sin
                                      / max(unstable["deviation"], stable["deviation"], 1e-300))
    others = np.array([i not in set(orbit) for i in detections["ids"]], dtype=bool)
    if others.any():
        dist = np.where(others, np.hypot(*(detections["coords"] - target).T), np.inf)
        i = int(np.argmin(dist))
        nxt["detection"] = {"distance": float(dist[i]), "ratio": detections["ratio"][i],
                            "resolved": bool(detections["resolved"][i])}
    out["next"] = nxt
    return out


def detection_table(workbench) -> dict:
    """Every registered detection: coordinates, ratio, resolved flag (anchors left out)."""
    registry = workbench.intersection_registry
    resolved = registry.resolved_ids(workbench.resolution_ratio_threshold)
    ids = [i for i in registry.all_ids() if registry[i].resolution_ratio is not None]
    return {"ids": ids,
            "coords": np.array([registry[i].coords for i in ids], dtype=float).reshape(-1, 2),
            "ratio": [float(registry[i].resolution_ratio) for i in ids],
            "resolved": [i in resolved for i in ids]}


def curve_stats(curves: list) -> dict:
    """Chord deviation, segment length and relative cdist step over a set of image curves."""
    lengths, devs, steps = [], [], []
    for points, cdist in curves:
        length, dev = deviations(points)
        lengths.append(length)
        devs.append(dev)
        steps.append(relative_steps(cdist))
    cat = (lambda parts: np.concatenate(parts) if parts else np.empty(0))
    step = cat(steps)
    return {"length": percentiles(cat(lengths)), "deviation": percentiles(cat(devs)),
            "relative_step": percentiles(step), "ties": int(np.sum(step == 0)),
            "points": int(sum(len(p) for p, _ in curves))}


def histograms(session) -> dict:
    """End of sweep: log10 histograms of chord deviation and segment length per stability, ratios."""
    workbench = session.workbench
    counts = {f"{kind}_{s}": np.zeros(len(LOG_BINS) - 1, dtype=np.int64)
              for kind in ("deviation", "length") for s in ("stable", "unstable", "image")}

    def add(points, stability):
        length, dev = deviations(points)
        for kind, values in (("length", length), ("deviation", dev)):
            values = values[values > 0]
            counts[f"{kind}_{stability}"] += np.histogram(np.log10(values), LOG_BINS)[0]

    for key, manifold in session.manifolds.items():
        add(manifold.get_point_array(), key[1])
    for curve in workbench._image_curves:  # the blast's computed curves (no public accessor)
        add(curve.get_point_array(), "image")
    table = detection_table(workbench)
    return {"bins": LOG_BINS.tolist(), **{k: v.tolist() for k, v in counts.items()},
            "ratios": table["ratio"], "resolved": table["resolved"]}


# ── the sweep ───────────────────────────────────────────────────────────────


def data_path(case: str, cutoff: float) -> Path:
    """``errors/data/<case>_a<cutoff>.json``."""
    return DATA_DIR / f"{case}_a{cutoff:.0e}.json"


def sweep(case: str, cutoff: float, max_blasts: int, max_seconds: Optional[float],
          max_rss_gb: Optional[float]) -> dict:
    """Build one case, blast it until its frontier is exhausted, and record the errors of every step."""
    path = data_path(case, cutoff)
    path.parent.mkdir(parents=True, exist_ok=True)
    record = {"case": case, "area_cutoff": cutoff, "min_separation": entropy.MIN_SEPARATION,
              "rows": [], "stopped": None}
    start = time.perf_counter()
    mapping = henon_map(entropy.K, entropy.B)
    try:
        session, fixed_points, zones = entropy.build_case(case, cutoff)
        workbench = session.workbench
        reports = []
        original = workbench.report_resolution

        def capture(*args, **kwargs):
            report = original(*args, **kwargs)
            reports.append(report)
            return report

        workbench.report_resolution = capture
        seen = set(detection_table(workbench)["ids"])
        previous_peak = entropy.resident_gb("VmHWM")
        for blasts in range(1, max_blasts + 1):
            if max_seconds is not None and time.perf_counter() - start > max_seconds:
                record["stopped"] = f"time limit ({max_seconds:g} s) before blast {blasts}"
                break
            peak = entropy.resident_gb("VmHWM")
            predicted = peak * peak / max(previous_peak, 1e-9)
            if max_rss_gb is not None and predicted > max_rss_gb:
                record["stopped"] = f"memory limit ({max_rss_gb:g} GB): blast {blasts} predicted at {predicted:.1f} GB"
                break
            previous_peak = peak
            row = {"blasts": blasts, "zones": {}}
            before = {f"p{fp.period}": len(orbit_of(session, fp)) for fp in fixed_points}
            new_curves = []
            reports.clear()
            for zone in zones:
                result = session.blast_zone(zone, num_iterations=1, fixed_point=[zone.fixed_point],
                                            min_separation=entropy.MIN_SEPARATION, diagnose=True)
                produced = [b for step in result.steps for b in step.produced_children]
                kept = [b for step in result.steps for b in step.kept_interior]
                curves = [(b.get_point_array(), b.get_cdist_array()) for b in produced]
                new_curves += [points for points, _ in curves]
                period = f"p{zone.fixed_point.period}"
                fp = zone.fixed_point
                near_P = [float(np.min(np.hypot(*(p - fp.coordinates[o]).T)))
                          for p, _ in curves if len(p) for o in range(fp.period)]
                row["zones"][period] = {
                    "pip_depth": {f"p{f.period}": d for f, d in result.pip_depth.items()},
                    "produced": len(produced), "kept": len(kept), "too_close": result.too_close,
                    "nearest_node_to_P": min(near_P) if near_P else None,
                    "curves": curve_stats(curves),
                }
            table = detection_table(workbench)
            new = [k for k, i in enumerate(table["ids"]) if i not in seen]
            seen.update(table["ids"])
            ratios = np.array([table["ratio"][k] for k in new], dtype=float)
            row["new_detections"] = len(new)
            row["new_ratio"] = percentiles(ratios)
            row["new_below_threshold"] = int(np.sum(ratios < workbench.resolution_ratio_threshold))
            row["unresolved"] = len(table["ids"]) - sum(table["resolved"])
            row["same_stability"] = sum(r.same_stability_count.in_region for r in reports
                                        if r.same_stability_count is not None)
            row["orbits"] = {f"p{fp.period}": orbit_rows(session, fp, mapping, new_curves, table)
                             for fp in fixed_points}
            for period, orbit in row["orbits"].items():
                orbit["advanced"] = len(orbit["orbit"]) > before[period]
            row["points"] = len(workbench.point_store)
            row["peak_gb"] = entropy.resident_gb("VmHWM")
            row["seconds"] = time.perf_counter() - start
            record["rows"].append(row)
            path.write_text(json.dumps(record, indent=1))
            depth = {p: z["pip_depth"] for p, z in row["zones"].items()}
            print(f"[{case} a={cutoff:g}] blast {blasts}: depth {depth}, unresolved {row['unresolved']}, "
                  f"{row['peak_gb']:.1f} GB, {row['seconds']:.0f} s", flush=True)
            if not any(z["kept"] for z in row["zones"].values()):
                record["stopped"] = f"frontier exhausted at blast {blasts}"
                break
        record["histograms"] = histograms(session)
    except Exception as exc:  # noqa: BLE001 -- a failing step ends the sweep, recorded
        record["stopped"] = f"{type(exc).__name__}: {exc}"
        logger.warning("%s a=%g stopped: %s\n%s", case, cutoff, exc, traceback.format_exc())
    path.write_text(json.dumps(record, indent=1))
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description="Error analysis of the k=2.5 entropy sweeps.")
    parser.add_argument("--cases", nargs="+", choices=entropy.CASES, default=list(entropy.CASES))
    parser.add_argument("--cutoffs", nargs="+", type=float, default=[1e-7])
    parser.add_argument("--max-blasts", type=int, default=100)
    parser.add_argument("--max-seconds", type=float, default=None)
    parser.add_argument("--max-rss-gb", type=float, default=64.0)
    parser.add_argument("--plot-only", action="store_true")
    args = parser.parse_args()
    configure_logging()
    plt.switch_backend("Agg")
    if not args.plot_only:
        for cutoff in args.cutoffs:
            for case in args.cases:
                sweep(case, cutoff, args.max_blasts, args.max_seconds, args.max_rss_gb)
    import henon_k25_errors_figures as figures  # noqa: PLC0415 -- the figure half
    figures.draw()


if __name__ == "__main__":
    main()
