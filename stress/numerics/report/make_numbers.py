"""Headline numbers, tables and the findings digest of the numerics report.

    env/bin/python stress/numerics/report/make_numbers.py [--quick]

Writes

* ``numbers.tex`` -- one ``\\newcommand{\\NumXxx}{value}`` per number the
  report quotes, each under a ``% source:`` line naming the experiment and how
  the value was computed; a missing value is ``n/a``;
* ``tables/*.tex`` -- booktabs tabular BODIES (``\\toprule`` .. ``\\bottomrule``;
  the column spec is suggested in the first comment line);
* ``summary.json`` -- per experiment: the results file used, outcome counts,
  every finding under its macro name (raw value) and structured details.

Every value comes from ``analysis``, which the figures read too.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
import time
from pathlib import Path
from typing import Any, Optional

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analysis as A  # noqa: E402

REPORT = Path(__file__).resolve().parent
TABLES = REPORT / "tables"

CASE_WORDS = {"k28": "KTwoEight", "p3": "PThree", "k10": "KTen", "nested": "Nested", "inversion": "Inversion"}
DIGITS = {1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five", 6: "Six", 7: "Seven", 8: "Eight", 9: "Nine",
          10: "Ten", 12: "Twelve"}
K_WORDS = {2.8: "KTwoEight", 4.0: "KFour", 6.0: "KSix", 10.0: "KTen"}

MACROS: list[tuple[str, str, str]] = []
SUMMARY: dict[str, Any] = {}


# --------------------------------------------------------------------------- #
# Formatting
# --------------------------------------------------------------------------- #
def sci_tex(x: float, digits: int = 2) -> str:
    mantissa, exponent = f"{x:.{digits - 1}e}".split("e")
    exponent = int(exponent)
    if math.isclose(abs(x), 10.0**exponent, rel_tol=1e-12):
        return f"{'-' if x < 0 else ''}10^{{{exponent}}}"
    return f"{mantissa}\\times10^{{{exponent}}}"


def tex(value: Any, kind: str = "auto") -> str:
    """A value as LaTeX: n/a, yes/no, escaped text or a number in the given style."""
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return "n/a"
    if isinstance(value, (bool, np.bool_)):
        return "yes" if value else "no"
    if isinstance(value, str):
        return value.replace("\\", "").replace("_", "\\_").replace("%", "\\%").replace("&", "\\&")
    if kind == "int" or (kind == "auto" and isinstance(value, (int, np.integer))):
        return f"{int(round(value)):,}".replace(",", "\\,")
    v = float(value)
    if v == 0:
        return "0"
    if kind == "pct":
        return f"{100 * v:.3g}\\%"
    if kind == "x":  # 3 decimals: an Amdahl bound of 1.004 must not print as 1
        return f"\\ensuremath{{{v:.3f}\\times}}" if v < 10 else f"\\ensuremath{{{v:.3g}\\times}}"
    if not 1e-2 <= abs(v) < 1e4 or kind == "sci" and abs(v) >= 1e3:
        return f"\\ensuremath{{{sci_tex(v)}}}"
    return f"{v:.3g}"


def put(exp: str, name: str, value: Any, source: str, kind: str = "auto", text: Optional[str] = None) -> None:
    """One finding: a ``\\Num<name>`` macro and its raw value in summary.json.

    ``text`` replaces the formatted macro body (LaTeX, not escaped).
    """
    assert re.fullmatch(r"[A-Za-z]+", name), name
    if isinstance(value, (np.floating, np.integer)):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        value = None
    MACROS.append((f"Num{name}", text or tex(value, kind), f"{A.EXPERIMENTS[exp]} {source}"))
    section(exp)["findings"][name] = value


def detail(exp: str, key: str, value: Any) -> None:
    section(exp)["details"][key] = value


def section(exp: str) -> dict:
    return SUMMARY.setdefault(exp, {"findings": {}, "details": {}, "notes": []})


def note(exp: str, text: str) -> None:
    section(exp)["notes"].append(text)


def write_table(name: str, colspec: str, header: list[str], rows: list[list[str]], caption: str) -> None:
    """A booktabs tabular body (no begin/end tabular) under a comment with its column spec."""
    lines = [f"% {caption}", f"% suggested column spec: {{{colspec}}}", "\\toprule",
             " & ".join(header) + " \\\\", "\\midrule"]
    lines += [" & ".join(row) + " \\\\" for row in rows] or ["\\multicolumn{%d}{c}{no data} \\\\" % len(header)]
    lines.append("\\bottomrule")
    (TABLES / f"{name}.tex").write_text("\n".join(lines) + "\n")
    print(f"  wrote tables/{name}.tex ({len(rows)} rows)")


# --------------------------------------------------------------------------- #
# E1
# --------------------------------------------------------------------------- #
def e01() -> None:
    sweeps, crossovers = A.e01_sweeps(), A.e01_crossovers()
    for pair, word in (("raw", "Raw"), ("map_batch", "MapBatch")):
        c = crossovers.get(pair)
        if c is None or not c["measured"]:
            put("e01", f"GpuCrossoverN{word}", None, f"crossover config '{pair}' missing or unfinished")
        elif c["n"] is None:
            put("e01", f"GpuCrossoverN{word}", f"none up to N = {c['n_max']:g}",
                f"crossover '{pair}': GPU median never below CPU median up to hi",
                text=f"none up to \\ensuremath{{N = {sci_tex(c['n_max'])}}}")
        else:
            put("e01", f"GpuCrossoverN{word}", c["n"], f"crossover '{pair}': bisection in log N of median GPU/CPU time",
                "sci")
        detail("e01", f"crossover_{pair}", c)
    for cpu, gpu, word in (("numpy", "cupy_e2e", "Raw"), ("map_batch_cpu", "map_batch_gpu", "MapBatch")):
        n = A.e01_sweep_crossover(sweeps[cpu], sweeps[gpu]) if cpu in sweeps and gpu in sweeps else None
        put("e01", f"GpuCrossoverNSweep{word}", n, f"sweep {cpu} vs {gpu}: log-log interpolated ratio = 1", "sci")
    for backend, word in (("numpy", "Numpy"), ("cupy_e2e", "CupyEndToEnd"), ("cupy_compute", "CupyCompute"),
                          ("map_batch_cpu", "MapBatchCpu"), ("map_batch_gpu", "MapBatchGpu")):
        d = sweeps.get(backend)
        put("e01", f"MapPeak{word}", None if d is None else A.nanmax(d["points_per_s"]),
            f"sweep {backend}: max over N of N / median call time (points/s)", "sci")
    split = A.e01_split(sweeps)
    if split is not None:
        put("e01", "GpuTransferFractionAtMaxN", split["transfer"][-1],
            f"cupy_transfer / cupy_e2e median at the largest N = {split['n'][-1]:.0f}", "pct")
        detail("e01", "split", {k: v.tolist() for k, v in split.items()})
    else:
        put("e01", "GpuTransferFractionAtMaxN", None, "needs cupy_e2e/compute/transfer sweeps")


# --------------------------------------------------------------------------- #
# E2
# --------------------------------------------------------------------------- #
def e02() -> None:
    timed = A.e02_timed()
    curves = A.e02_step_curves(timed)
    best = A.e02_best_step_speedup(curves)
    totals = A.e02_total_speedups(timed)
    stage_speedups = A.e02_stage_speedups(timed)
    profiles = A.e02_profiles()
    for case in A.CASES:
        w = CASE_WORDS[case]
        b = best.get(case)
        put("e02", f"GpuStepSpeedupBest{w}", b and b["speedup"],
            "timed: max over growth steps (CPU median >= 10 ms) of median CPU / median GPU step time", "x")
        v = totals.get(case) or []
        put("e02", f"GpuTotalSpeedupMedian{w}", A.nanmedian(v),
            "timed: median over repeats of CPU / GPU summed stage time (growth through iterate table)", "x")
        put("e02", f"GpuTotalSpeedupBest{w}", A.nanmax(v), "timed: max over repeats of CPU / GPU summed stage time", "x")
        p = profiles.get((case, "cpu"))
        put("e02", f"MapFraction{w}", p and p["f_map"], "profile, CPU: map-category tottime / all tottime", "pct")
        put("e02", f"AmdahlBound{w}", p and p["amdahl"], "profile, CPU: 1 / (1 - f_map)", "x")
    cpu_profiles = [p for (c, b), p in profiles.items() if b == "cpu"]
    put("e02", "MapFractionMax", A.nanmax(p["f_map"] for p in cpu_profiles), "profile, CPU: max f_map over cases", "pct")
    put("e02", "AmdahlBoundMax", A.nanmax(p["amdahl"] for p in cpu_profiles), "profile, CPU: max 1/(1-f_map)", "x")
    detail("e02", "best_step_speedup", best)
    detail("e02", "total_speedups", totals)
    detail("e02", "stage_speedups", stage_speedups)
    detail("e02", "profiles", {f"{c}/{b}": p for (c, b), p in profiles.items()})

    verify = A.ok("e02", kind="verify")
    put("e02", "GpuDeviationMax", A.nanmax(A.num(r["metrics"].get("final_max_abs_deviation")) for r in verify),
        "verify: max over cases of the final max |CPU - GPU| unstable coordinate", "sci")
    put("e02", "GpuCountsEqual", all(r["metrics"].get("final_counts_equal") for r in verify) if verify else None,
        "verify: CPU and GPU point counts equal at the last step in every case")
    detail("e02", "verify", {r["config"]["case"]: {k: r["metrics"].get(k) for k in
                             ("final_counts_equal", "final_max_abs_deviation", "first_divergent_step")} for r in verify})

    mbp = sorted(A.ok("e02", kind="mbp"), key=lambda r: A.num(r["stages"].get("grow_unstable")))
    if mbp:
        never = next((r for r in mbp if r["config"]["min_batch_points"] >= 1e9), None)
        best = mbp[0]["config"]["min_batch_points"]
        put("e02", "MinBatchBest", best, "mbp: min_batch_points with the smallest k28 unstable growth time "
            "(1e9 = never on GPU)", "int", text="never (CPU only)" if best >= 1e9 else None)
        put("e02", "MinBatchBestGain", never and A.num(never["stages"]["grow_unstable"]) /
            A.num(mbp[0]["stages"]["grow_unstable"]), "mbp: growth time at 'never' / at the best setting", "x")
        detail("e02", "min_batch", {r["config"]["min_batch_points"]: r["stages"].get("grow_unstable") for r in mbp})
    else:
        put("e02", "MinBatchBest", None, "mbp: no ok records")
        put("e02", "MinBatchBestGain", None, "mbp: no ok records")

    rows = []
    for case, row in stage_speedups.items():
        p = profiles.get((case, "cpu"))
        rows.append([case] + [tex(row.get(s), "x") for s in [*A.E2_STAGES, "total"]]
                    + [tex(p and p["f_map"], "pct"), tex(p and p["amdahl"], "x")])
    write_table("e02_speedup", "l" + "r" * (len(A.E2_STAGES) + 3),
                ["case", *A.E2_STAGES, "total", "$f_{\\rm map}$", "Amdahl"], rows,
                "E2: median CPU time / median GPU time per stage, and the CPU map fraction with its Amdahl bound")


# --------------------------------------------------------------------------- #
# E3
# --------------------------------------------------------------------------- #
def e03() -> None:
    growth = A.e03_growth()
    for (k, cutoff), g in sorted(growth.items()):
        w = K_WORDS.get(k, f"K{len(str(k))}") + ("" if cutoff == 1e-7 else "Coarse")
        put("e03", f"StepTimeExponent{w}", g["step_fit"] and g["step_fit"][0],
            f"growth k={k:g} cutoff {cutoff:g}: log-log slope of step time vs unstable points (>= 1e3 points)")
        put("e03", f"BytesPerPoint{w}", g["bytes_per_point_fit"],
            f"growth k={k:g} cutoff {cutoff:g}: slope of (RSS - baseline) vs points over the top two decades", "int")
        blowup = next((int(s) for s, c in zip(g["step"], g["max_abs_coord"]) if c > 100), None)
        put("e03", f"EscapeStep{w}", blowup, f"growth k={k:g} cutoff {cutoff:g}: first step with max|coord| > 100")
        detail("e03", f"growth_k{k:g}_c{cutoff:g}", {
            "steps": len(g["step"]), "max_points": A.nanmax(g["points_unstable"]), "stop_reason": g["stop_reason"],
            "outcome": g["outcome"], "lambda_u": g["lambda_u"], "last_growth_factor": float(g["growth_factor"][-1]),
            "step_fit": g["step_fit"], "bytes_per_point_fit": g["bytes_per_point_fit"]})
    exps = [g["step_fit"][0] for g in growth.values() if g["step_fit"]]
    put("e03", "StepTimeExponent", A.nanmedian(exps), "growth: median over (k, cutoff) of the step-time slope")
    mems = [g["bytes_per_point_fit"] for g in growth.values() if g["bytes_per_point_fit"]]
    put("e03", "BytesPerPoint", A.nanmedian(mems), "growth: median over (k, cutoff) of the RSS slope (bytes/point)", "int")
    put("e03", "BytesPerPointMax", A.nanmax(mems), "growth: max over (k, cutoff) of the RSS slope (bytes/point)", "int")
    rows = [r for r in A.e03_stages() if r["outcome"] == "ok"]
    fits = A.e03_stage_fits(rows)
    for stage, word in (("compute_intersections", "Intersections"), ("trim_and_bridges", "TrimBridges"),
                        ("infer_iterate_table", "IterateTable")):
        put("e03", f"{word}Exponent", fits[stage] and fits[stage][0],
            f"stages: log-log slope of {stage} time vs segments, pooled over k (>= 1e3 segments)")
    detail("e03", "stage_rows", rows)
    if any(r["outcome"] != "ok" for r in A.load("e03")):
        note("e03", f"non-ok records: {A.outcome_counts('e03')}")


# --------------------------------------------------------------------------- #
# E4
# --------------------------------------------------------------------------- #
def e04() -> None:
    runs = A.e04_runs()
    canonical = {}
    for case, steps in runs:
        canonical.setdefault(case, steps)
        canonical[case] = min(canonical[case], steps)
    for case in ("k28", "p3"):
        w = CASE_WORDS[case]
        rows = runs.get((case, canonical.get(case)), [])
        conv = {c["cutoff"]: c for c in A.e04_convergence(rows)}
        for cutoff in (1e-7, 1e-4):
            cw = f"Cutoff{DIGITS[round(-math.log10(cutoff))]}"
            r = next((r for r in rows if r["cutoff"] == cutoff and r["outcome"] == "ok"), None)
            where = f"{case} {canonical.get(case)} steps, cutoff {cutoff:g}"
            put("e04", f"InvPFifty{w}{cw}", r and r["inv_p50"], f"{where}: max over branches of the p50 invariance error", "sci")
            put("e04", f"InvPNinetyNine{w}{cw}", r and r["inv_p99"], f"{where}: max over branches of the p99 invariance error", "sci")
            put("e04", f"Time{w}{cw}", r and r["time"], f"{where}: summed pipeline stages (no measurement stages), s")
            put("e04", f"Points{w}{cw}", r and r["points"], f"{where}: unstable + stable points", "int")
            c = conv.get(cutoff)
            put("e04", f"CrossingDevMedian{w}{cw}", c and c["median"],
                f"{where}: median NN distance of the finest run's crossings to this run's", "sci")
        done = [r for r in rows if np.isfinite(r["points"])]
        fit = A.fit_loglog([r["cutoff"] for r in done], [r["points"] for r in done])
        put("e04", f"PointsCutoffExponent{w}", fit and fit[0], f"{case} canonical depth: log-log slope of points vs cutoff")
    detail("e04", "runs", {f"{c}/{s}": [{k: v for k, v in r.items() if k not in ("record", "crossing_rows")}
                                       for r in rows] for (c, s), rows in runs.items()})
    detail("e04", "convergence", {f"{c}/{s}": A.e04_convergence(rows) for (c, s), rows in runs.items()})

    table = []
    for (case, steps), rows in runs.items():
        conv = {c["cutoff"]: c for c in A.e04_convergence(rows)}
        for r in rows:
            c = conv.get(r["cutoff"])
            table.append([case, str(steps), tex(r["cutoff"], "sci"), r["outcome"].replace("_", "\\_"),
                          tex(r["points"], "int"), tex(r["time"]), tex(r["inv_p50"], "sci"), tex(r["inv_p99"], "sci"),
                          tex(r["crossings"], "int"), "ref." if c is None and r["outcome"] == "ok" and conv else
                          tex(c and c["median"], "sci")])
    write_table("e04_cutoff", "llrlrrrrrr",
                ["case", "steps", "cutoff", "outcome", "points", "time (s)", "inv.\\ p50", "inv.\\ p99",
                 "crossings", "crossing dev."], table,
                "E4: cost and accuracy per cutoff; crossing dev. = median distance to the finest run's crossings")


# --------------------------------------------------------------------------- #
# E5
# --------------------------------------------------------------------------- #
def e05() -> None:
    records = A.e05_cases()
    mets = [r["metrics"] for r in records]
    links = [A.num(row[3]) for r in records for row in r["series"].get("iter_links") or []]
    put("e05", "IterLinkErrMedian", A.nanmedian(links), "series iter_links (<= 2000 per case): median |f^n(x) - x_target|", "sci")
    put("e05", "IterLinkErrMax", A.nanmax(v for mm in mets for k, v in mm.items() if re.fullmatch(r"iter_n-?\d+_max", k)),
        "metrics iter_n<n>_max: max over cases and n", "sci")
    put("e05", "IterLinks", sum(A.num(mm.get("iter_links", 0)) for mm in mets), "metrics iter_links: total checked", "int")
    put("e05", "LobeRatioMax", A.nanmax(mm.get("lobe_ratio_max") for mm in mets), "metrics lobe_ratio_max: max |A_img/A - 1|", "sci")
    put("e05", "LobePairs", sum(A.num(mm.get("lobe_pairs", 0)) for mm in mets), "metrics lobe_pairs: total", "int")
    put("e05", "LobeSignFlips", sum(A.num(mm.get("lobe_sign_flips", 0)) for mm in mets), "metrics lobe_sign_flips: total", "int")
    orbits = [o for r in records for o in A.e05_fp_rows(r)]
    put("e05", "FpResidualMax", A.nanmax(o.get("fp_residual") for o in orbits), "metrics fp: max |f^p(x*) - x*|", "sci")
    put("e05", "DetDefectMax", A.nanmax(abs(A.num(o.get("det_minus_1"))) for o in orbits), "metrics fp: max ||lu ls| - 1|", "sci")
    put("e05", "EvecResidualMax", A.nanmax(max(A.num(o.get("evec_res_u")), A.num(o.get("evec_res_s"))) for o in orbits),
        "metrics fp: max |J v - lambda v| / |v|, both eigenvectors", "sci")
    put("e05", "JacobianDeviationMax", A.nanmax(o.get("jac_vs_analytic") for o in orbits),
        "metrics fp: max |stored - analytic| cycle Jacobian entry", "sci")
    inv_max = [v for mm in mets for k, v in mm.items() if k.startswith("inv_") and k.endswith("_max")]
    put("e05", "InvarianceMax", A.nanmax(inv_max), "metrics inv_*_max: largest invariance error over cases and branches", "sci")
    put("e05", "InvarianceBeyond", sum(A.num(v) for mm in mets for k, v in mm.items()
                                       if k.startswith("inv_") and k.endswith("_beyond")),
        "metrics inv_*_beyond: images past the image branch's tail (expected 0)", "int")
    for scope, word in (("cdist_manifolds", "Manifolds"), ("cdist_bridges", "Bridges")):
        for key, kw in (("decreasing", "Decreasing"), ("ties", "Ties"), ("steps", "Steps")):
            put("e05", f"Cdist{kw}{word}", sum(A.num((mm.get(scope) or {}).get(key, 0)) for mm in mets),
                f"metrics {scope}.{key}: total over cases", "int")
    put("e05", "ArrImageFound", sum(A.num(mm.get("arr_image_found", 0)) for mm in mets), "metrics arr_image_found: total", "int")
    put("e05", "ArrImageMissing", sum(A.num(mm.get("arr_image_missing", 0)) for mm in mets), "metrics arr_image_missing: total", "int")

    gpu = next((r for r in A.load("e05") if r["config"].get("case") == "gpu_k28"), None)
    gm = gpu["metrics"] if gpu else {}
    put("e05", "GpuHausdorffMax", A.nanmax((gm.get("hausdorff") or {}).values()),
        "gpu_k28: max polyline Hausdorff distance CPU vs GPU over manifolds", "sci")
    put("e05", "GpuSameLengths", all((gm.get("same_lengths") or {}).values()) if gm.get("same_lengths") else None,
        "gpu_k28: every manifold has the same point count on CPU and GPU")

    rows = []
    for r in records:
        mm, o = r["metrics"], A.e05_fp_rows(r)
        cm, cb = mm.get("cdist_manifolds") or {}, mm.get("cdist_bridges") or {}
        inv = [v for k, v in mm.items() if k.startswith("inv_") and k.endswith("_max")]
        rows.append([
            r["config"]["case"].replace("_", " "), r["outcome"].replace("_", "\\_"),
            tex(A.nanmax(x.get("fp_residual") for x in o), "sci"),
            tex(A.nanmax(abs(A.num(x.get("det_minus_1"))) for x in o), "sci"),
            tex(A.nanmax(max(A.num(x.get("evec_res_u")), A.num(x.get("evec_res_s"))) for x in o), "sci"),
            tex(A.nanmax(inv), "sci"), tex(mm.get("lobe_ratio_max"), "sci"),
            f"{tex(cm.get('decreasing'), 'int')}/{tex(cm.get('ties'), 'int')}",
            f"{tex(cb.get('decreasing'), 'int')}/{tex(cb.get('ties'), 'int')}",
            f"{tex(mm.get('arr_image_found'), 'int')}/{tex(mm.get('arr_image_missing'), 'int')}",
        ])
    write_table("e05_invariants", "llrrrrrccc",
                ["case", "outcome", "$|f^p(x^*)-x^*|$", "$||\\lambda_u\\lambda_s|-1|$", "eigvec res.",
                 "inv.\\ max", "lobe $|A'/A-1|$", "cdist dec./ties (man.)", "cdist dec./ties (bridges)",
                 "arr.\\ image found/missing"], rows,
                "E5: per-case invariants; maxima over the orbit points and the branches")
    detail("e05", "per_case", {r["config"]["case"]: {k: v for k, v in r["metrics"].items()
                                                    if not isinstance(v, (list, dict)) or k in ("cdist_manifolds", "cdist_bridges")}
                               for r in records})


# --------------------------------------------------------------------------- #
# E6
# --------------------------------------------------------------------------- #
STAGE_WORDS = {"none": "None", "fixed_point": "FixedPoint", "eigen": "Eigen", "seeded": "Seeded",
               "grown_unstable": "GrownUnstable", "grown": "Grown", "intersected": "Intersected",
               "bridges": "Bridges", "trellis_pips": "TrellisPips"}


def e06() -> None:
    rows = A.e06_rows()
    complete = [r for r in rows if r["outcome"] == "ok" and r["stage"] == "trellis_pips"]
    put("e06", "SweepConfigs", len(rows), "configs run", "int")
    put("e06", "SweepComplete", len(complete), "configs ok through trellis_pips", "int")
    put("e06", "SweepSuccessFraction", len(complete) / len(rows) if rows else None,
        "configs ok through trellis_pips / all configs", "pct")
    saddle = [r for r in rows if r["family"] == "b = 1, saddle"]
    saddle_ok = [r for r in saddle if r in complete]
    put("e06", "SaddleSuccessFraction", len(saddle_ok) / len(saddle) if saddle else None,
        "b = 1 saddle configs ok through trellis_pips / b = 1 saddle configs", "pct")
    put("e06", "SaddleMinKComplete", A.nanmin([r["k"] for r in saddle_ok]) if saddle_ok else None,
        "smallest k whose b = 1 saddle completes")
    for stage, word in STAGE_WORDS.items():
        put("e06", f"SweepStage{word}", sum(r["stage"] == stage for r in rows),
            f"configs whose last completed stage is {stage}", "int")
    for outcome in A.OUTCOMES:
        word = "".join(p.capitalize() for p in outcome.split("_"))
        put("e06", f"SweepOutcome{word}", sum(r["outcome"] == outcome for r in rows), f"configs with outcome {outcome}", "int")
    families = sorted({r["family"] for r in rows})
    table = []
    for fam in families:
        mine = [r for r in rows if r["family"] == fam]
        table.append([fam.replace("-", "$-$")] + [str(sum(r["stage"] == s for r in mine)) for s in A.E6_STAGES]
                     + [str(len(mine))])
    write_table("e06_stages", "l" + "r" * (len(A.E6_STAGES) + 1),
                ["family", *[s.replace("_", " ") for s in A.E6_STAGES], "total"], table,
                "E6: configs per family by last stage completed")
    errors: dict[str, int] = {}
    for r in rows:
        if r["outcome"] != "ok":
            key = f"{r['family']} | {(r['error'] or '')[:90]}"
            errors[key] = errors.get(key, 0) + 1
    detail("e06", "errors", errors)
    detail("e06", "rows", [{k: r[k] for k in ("family", "k", "outcome", "stage")} for r in rows])


# --------------------------------------------------------------------------- #
# E7
# --------------------------------------------------------------------------- #
def e07() -> None:
    escape = A.e07("escape")
    fixed = next((r for r in escape if r["config"].get("mode") == "fixed"), None)
    s = fixed["series"] if fixed else {}
    put("e07", "EscapeOutcome", fixed and fixed["outcome"], "escape fixed: record outcome")
    put("e07", "EscapeStepsDone", fixed and fixed["metrics"].get("steps_done"), "escape fixed: steps completed", "int")
    put("e07", "EscapeBlowupStep", fixed and A.e07_blowup_step(fixed), "escape fixed: first step with max|coord| > 100", "int")
    put("e07", "EscapePeakPoints", A.nanmax(s.get("points", [])), "escape fixed: max unstable points over steps", "int")
    put("e07", "EscapePeakRssMb", fixed and max(A.nanmax(s.get("rss_mb", [])) or 0, A.num(fixed.get("peak_rss_mb"))),
        "escape fixed: peak RSS (MB), max of per-step samples and the watchdog's", "int")
    put("e07", "EscapeMaxCoord", A.nanmax(s.get("max_abs", [])), "escape fixed: largest max|coord|", "sci")
    for r in escape:
        mode = r["config"].get("mode")
        if mode == "fixed":
            continue
        word = "Arclength" if mode == "arclength" else "Predicate"
        arg = r["config"].get("length", r["config"].get("radius"))
        suffix = "" if arg in (0.5, 2.5) else DIGITS.get(int(arg), "Other")
        put("e07", f"EscapeStop{word}{suffix}Points", r["metrics"].get("final_points"),
            f"escape {mode} {arg}: final unstable points ({r['outcome']})", "int")
    for r in A.e07("blasts"):
        case, mm = r["config"]["case"], r["metrics"]
        w = CASE_WORDS.get(case, case.capitalize())
        put("e07", f"MaxBlasts{w}", mm.get("blasts_done", 0), f"blasts {case}: blasts completed of {r['config']['max_blasts']}", "int")
        put("e07", f"BlastOutcome{w}", r["outcome"], f"blasts {case}: record outcome")
        topo = mm.get("topology_outcome")
        put("e07", f"BlastTopology{w}", None if topo is None else topo + ("" if mm.get("is_reliable", True) else " (unreliable)"),
            f"blasts {case}: topology stack on the final state")
        detail("e07", f"blasts_{case}", {"outcome": r["outcome"], "error": r.get("error"),
                                          **{k: mm.get(k) for k in ("blasts_done", "attempted_blast", "topology_outcome",
                                                                    "topology_error", "is_reliable", "n_classes",
                                                                    "n_unresolved", "n_virtual")}})
    depth = A.e07("depth")
    for case in ("p3", "inversion"):
        mine = [r for r in depth if r["config"].get("case") == case]
        w = CASE_WORDS[case]
        num_ok = [r["config"]["unstable_steps"] for r in mine if r["metrics"].get("numerics_done")]
        topo_ok = [r["config"]["unstable_steps"] for r in mine if r["metrics"].get("topology_outcome") == "ok"]
        put("e07", f"DepthNumericsMax{w}", max(num_ok) if num_ok else None, f"depth {case}: deepest run through bridges", "int")
        put("e07", f"DepthTopologyMax{w}", max(topo_ok) if topo_ok else None, f"depth {case}: deepest run with topology ok", "int")
        detail("e07", f"depth_{case}", {r["config"]["unstable_steps"]: {
            "outcome": r["outcome"], "topology": r["metrics"].get("topology_outcome"),
            "reliable": r["metrics"].get("is_reliable"), "error": r.get("error") or r["metrics"].get("topology_error")}
            for r in mine})
    tang = [r for r in A.e07("tangency") if "min_sin" in r["metrics"]]
    put("e07", "TangencyMinSin", A.nanmin(r["metrics"]["min_sin"] for r in tang), "tangency: smallest |sin| at any crossing", "sci")
    put("e07", "TangencySinBelowSix", sum(A.num(r["metrics"].get("n_sin_below_1e-6", 0)) for r in tang),
        "tangency: crossings with |sin| < 1e-6, total", "int")
    put("e07", "TangencyNearParallelWarnings",
        sum(c for r in A.e07("tangency") for k, c in r["warnings"].items() if "parallel" in k.lower()),
        "tangency: WARNING records mentioning 'parallel'", "int")


# --------------------------------------------------------------------------- #
# E8
# --------------------------------------------------------------------------- #
def e08() -> None:
    table = A.e08_table()
    n = sum(r["n"] for r in table.values())
    conv = sum(r["converged"] for r in table.values())
    lower = sum(r["lower_period"] for r in table.values())
    counts = {c: sum(r["counts"][c] for r in table.values()) for c in A.E8_CLASSES}
    put("e08", "SolverSolves", n, "solves (seeds x 2 variants x configs)", "int")
    put("e08", "SolverConvergedRate", conv / n if n else None, "solves where compute_fixed_point returned / all", "pct")
    put("e08", "SolverSaddleRate", counts["saddle"] / n if n else None, "saddle of the requested minimal period / all", "pct")
    put("e08", "SolverWrongPeriodRate", lower / conv if conv else None, "converged with minimal period < requested / converged", "pct")
    put("e08", "SolverNoConvergenceRate", counts["no_convergence"] / n if n else None, "fsolve did not converge / all", "pct")
    put("e08", "SolverMaxResidual", A.nanmax(r["max_residual"] for r in table.values()), "max |f(x_{i-1}) - x_i| of converged orbits", "sci")
    distinct = A.e08_distinct_by_period(table).get("k = 10", {})
    periods = sorted(A.HORSESHOE_ORBITS)
    matched = sum(distinct.get(p, 0) == A.HORSESHOE_ORBITS[p] for p in periods)
    put("e08", "HorseshoeMatched", f"{matched}/{len(periods)}" if distinct else None,
        "k = 10: periods whose distinct orbits found equal the full 2-shift count")
    put("e08", "HorseshoeFound", sum(distinct.values()) if distinct else None, "k = 10: distinct orbits found, periods 1-6", "int")
    put("e08", "HorseshoeExpected", sum(A.HORSESHOE_ORBITS.values()), "full 2-shift orbits of minimal period 1-6", "int")
    detail("e08", "class_counts", counts)
    detail("e08", "distinct_by_period", A.e08_distinct_by_period(table))
    rows = []
    for (map_name, p), r in sorted(table.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        found = A.e08_distinct_by_period(table).get(map_name, {}).get(p, 0)
        expected = A.HORSESHOE_ORBITS.get(p) if map_name == "k = 10" else None
        rows.append([f"${map_name}$", str(p), tex(r["converged"] / r["n"] if r["n"] else None, "pct"),
                     tex(r["counts"]["saddle"] / r["n"] if r["n"] else None, "pct"),
                     tex(r["lower_period"] / r["converged"] if r["converged"] else None, "pct"),
                     str(found) + ("" if expected is None else f"/{expected}"), tex(r["max_residual"], "sci")])
    write_table("e08_solver", "lrrrrrr",
                ["map", "period", "converged", "saddle", "lower period", "orbits found", "max residual"], rows,
                "E8: solver outcomes per map and requested period; orbits found = distinct of this minimal period over all configs (k = 10: / full 2-shift count)")


# --------------------------------------------------------------------------- #
# E9
# --------------------------------------------------------------------------- #
def e09() -> None:
    rows = [r for r in A.e09_rows() if r["outcome"] == "ok"]
    tot = lambda key: sum(r[key] for r in rows)  # noqa: E731
    put("e09", "IsectConfigs", len(rows), "ok configs", "int")
    put("e09", "IsectBruteTotal", tot("brute"), "brute-force proper crossings, total", "int")
    put("e09", "IsectMatchedTotal", tot("matched"), "brute crossings matched by the registry, total", "int")
    put("e09", "IsectMissedTotal", sum(tot(k) for k in A.E9_MISSED), "missed (all classes), total", "int")
    put("e09", "IsectMissedUnexplained", tot("missed_unexplained"), "missed, unexplained, total", "int")
    put("e09", "IsectExtraTotal", sum(tot(k) for k in A.E9_EXTRA), "extra registry crossings (no anchors), total", "int")
    put("e09", "IsectExtraUnexplained", tot("extra_unexplained"), "extra, unexplained, total", "int")
    uu = [r["uu"] for r in rows if isinstance(r["uu"], (int, float))]
    put("e09", "IsectUUTotal", sum(uu) if uu else None, "u x u proper crossings (non-adjacent segments), total", "int")
    put("e09", "IsectUUSkipped", sum(r["uu"] is None for r in rows), "configs above UU_MAX_SEGMENTS (u x u not checked)", "int")
    seg = [r["segments"] for r in rows]
    for key, word in (("library_s", "Library"), ("brute_s", "Brute")):
        fit = A.fit_loglog(seg, [r[key] for r in rows])
        put("e09", f"Isect{word}Exponent", fit and fit[0], f"log-log slope of {key} vs total segments")
    put("e09", "IsectBruteOverLibraryMax", A.nanmax(r["brute_s"] / r["library_s"] for r in rows),
        "max over configs of brute time / library time", "x")
    table = []
    for case in dict.fromkeys(r["case"] for r in rows):
        mine = [r for r in rows if r["case"] == case]
        s = lambda key: f"{sum(r[key] for r in mine):.0f}"  # noqa: E731
        uu_vals = [r["uu"] for r in mine]
        uu_text = "n/a" if all(v is None for v in uu_vals) else f"{sum(v for v in uu_vals if v is not None):.0f}" + (
            "$^*$" if None in uu_vals else "")
        table.append([case, str(len(mine)), s("brute"), s("matched"), s("near_parallel"), s("zigzag"),
                      s("collision"), s("missed_unexplained"), s("extra_touch"), s("extra_unexplained"), uu_text])
    write_table("e09_intersections", "lrrrrrrrrrr",
                ["case", "configs", "brute", "matched", "near-par.", "zigzag", "collision", "missed unexpl.",
                 "extra touch", "extra unexpl.", "$u\\times u$"], table,
                "E9: crossings vs brute force summed per case; * = some configs too large for the u x u check")


# --------------------------------------------------------------------------- #
# E10
# --------------------------------------------------------------------------- #
def verdict(values: list[Optional[bool]]) -> Optional[str]:
    values = [v for v in values if v is not None]
    if not values:
        return None
    return "identical" if all(values) else f"differs in {len(values) - sum(values)} of {len(values)} hashes"


def e10() -> None:
    by_device = A.e10_by_device()
    for device, word in (("cpu", "Cpu"), ("gpu", "Gpu")):
        records = by_device.get(device, [])
        cvs = {s: A.cv(r["stages"].get(s) for r in records) for s in A.E10_STAGES}
        put("e10", f"TimingCv{word}Max", A.nanmax(v for v in cvs.values() if v is not None) if any(cvs.values()) else None,
            f"{device}: max over stages of the stage-time CV over repeats", "pct")
        put("e10", f"TimingCv{word}Total", A.cv(sum(A.num(v) for v in r["stages"].values()) for r in records),
            f"{device}: CV of the summed stage time over repeats", "pct")
        put("e10", f"Repeats{word}", len(records), f"{device}: ok repeats", "int")
        detail("e10", f"cv_{device}", cvs)
    ident = A.e10_identity(by_device)
    for pair, word in (("CPU-CPU", "CpuCpu"), ("GPU-GPU", "GpuGpu"), ("CPU-GPU", "CpuGpu")):
        put("e10", f"BitIdent{word}", verdict([ident[h]["exact"][pair] for h in A.E10_HASHES]),
            f"exact sha256 of all five arrays agree {pair} (n/a with fewer than two runs)")
        put("e10", f"BitIdent{word}Rounded", verdict([ident[h]["rounded"][pair] for h in A.E10_HASHES]),
            f"sha256 of the arrays rounded to 1e-12 agree {pair}")
    put("e10", "CrossingDevCpuGpu", A.e10_crossing_deviation(by_device),
        "max |coordinate| difference of the sorted crossing lists, every run vs the first CPU run", "sci")
    detail("e10", "identity", ident)
    yes = lambda v: "n/a" if v is None else ("yes" if v else "\\textbf{no}")  # noqa: E731
    rows = [[h.replace("_hash", ""), *[yes(ident[h][kind][pair]) for kind in ("exact", "rounded")
                                        for pair in ("CPU-CPU", "GPU-GPU", "CPU-GPU")]] for h in A.E10_HASHES]
    write_table("e10_identity", "lcccccc",
                ["array", "exact CPU--CPU", "exact GPU--GPU", "exact CPU--GPU",
                 "$10^{-12}$ CPU--CPU", "$10^{-12}$ GPU--GPU", "$10^{-12}$ CPU--GPU"], rows,
                "E10: do the sha256 hashes agree across repeats (exact bytes, and rounded to 1e-12)?")


def clean(value: Any) -> Any:
    """JSON-safe: arrays to lists, numpy scalars to Python, nan/inf to None, tuple keys to strings."""
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean(v) for v in value]
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_numbers() -> None:
    lines = ["% Generated by stress/numerics/report/make_numbers.py -- do not edit.",
             f"% mode: {'quick' if A.QUICK else 'full'}; generated {time.strftime('%Y-%m-%d %H:%M:%S')}",
             "% sources: " + ", ".join(f"{k}={v[0]}({v[2]})" for k, v in sorted(A.SOURCES.items())), ""]
    for name, value, source in MACROS:
        lines += [f"% source: {source}", f"\\newcommand{{\\{name}}}{{{value}}}"]
    (REPORT / "numbers.tex").write_text("\n".join(lines) + "\n")
    print(f"  wrote numbers.tex ({len(MACROS)} macros)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--quick", action="store_true", help="use the *_quick.jsonl results")
    args = parser.parse_args()
    A.QUICK = args.quick
    TABLES.mkdir(parents=True, exist_ok=True)
    for exp, fn in (("e01", e01), ("e02", e02), ("e03", e03), ("e04", e04), ("e05", e05),
                    ("e06", e06), ("e07", e07), ("e08", e08), ("e09", e09), ("e10", e10)):
        print(exp)
        try:
            fn()
        except Exception as exc:  # noqa: BLE001 - one broken section never stops the rest
            import traceback

            traceback.print_exc()
            note(exp, f"numbers step failed: {type(exc).__name__}: {exc}")
        source = A.SOURCES.get(exp, ("missing", "", 0))
        section(exp).update({"source": source[0], "file": source[1], "records": source[2],
                             "outcomes": A.outcome_counts(exp)})
    write_numbers()
    digest = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"), "mode": "quick" if A.QUICK else "full",
              "experiments": {A.EXPERIMENTS[k]: SUMMARY[k] for k in sorted(SUMMARY)}}
    (REPORT / "summary.json").write_text(json.dumps(clean(digest), indent=1))
    print("  wrote summary.json")


if __name__ == "__main__":
    main()
