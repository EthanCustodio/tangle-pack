"""Extra report numbers computed straight from the full-run JSONL records.

    env/bin/python stress/numerics/report/make_numbers_extra.py

Writes ``numbers_extra.tex``: numbers the report quotes that ``numbers.tex``
(``make_numbers.py``) does not carry. Each macro sits under a ``% source:``
line naming the results file and the computation. Read-only on ``results/``.
"""

from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path
from typing import Any, Callable

import numpy as np

REPORT = Path(__file__).resolve().parent
RESULTS = REPORT.parent / "results"
OUT = REPORT / "numbers_extra.tex"

MACROS: list[tuple[str, str, str]] = []


def load(name: str) -> list[dict]:
    """All records of one experiment's full-run JSONL file."""
    with open(RESULTS / f"{name}.jsonl") as fh:
        return [json.loads(line) for line in fh if line.strip()]


def sci(x: float, digits: int = 2) -> str:
    """A float as a math-mode scientific number (at least two significant digits)."""
    digits = max(digits, 2)
    mantissa, exponent = f"{x:.{digits - 1}e}".split("e")
    exponent = int(exponent)
    if math.isclose(abs(x), 10.0**exponent, rel_tol=1e-9):
        return f"\\ensuremath{{10^{{{exponent}}}}}"
    return f"\\ensuremath{{{mantissa}\\times10^{{{exponent}}}}}"


def grp(n: float) -> str:
    """An integer with thin-space thousands groups."""
    return f"{int(round(n)):,}".replace(",", "\\,")


def pct(x: float, digits: int = 1) -> str:
    """A fraction as a percentage."""
    return f"{100 * x:.{digits}f}\\%"


def fixed(x: float, digits: int = 2) -> str:
    """A float with a fixed number of decimals."""
    return f"{x:.{digits}f}"


def times(x: float, digits: int = 2) -> str:
    """A ratio as a speed-up factor."""
    return f"\\ensuremath{{{x:.{digits}f}\\times}}"


def add(name: str, value: str, source: str) -> None:
    """Record one macro."""
    MACROS.append((name, value, source))


def slope(x: list[float], y: list[float]) -> float:
    """Log-log least-squares slope."""
    return float(np.polyfit(np.log10(x), np.log10(y), 1)[0])


# --------------------------------------------------------------------------- #
def e01() -> None:
    recs = load("e01_map_throughput")
    tp: dict[str, tuple[list, list]] = {}
    for r in recs:
        if r["config"]["kind"] != "sweep":
            continue
        s = r["series"]
        n = s["n"]
        tp[r["config"]["backend"]] = (n, [a / b for a, b in zip(n, s["median_s"])])
    n, t = tp["numpy"]
    i = int(np.argmax(t))
    add("NumMapPeakNumpyN", sci(n[i]), "e01 sweep numpy: batch size N of the peak throughput")
    add("NumMapPlateauNumpy", sci(t[-1]), "e01 sweep numpy: throughput at the largest N (points/s)")
    add("NumMapPeakOverPlateauNumpy", times(t[i] / t[-1], 1),
        "e01 sweep numpy: peak throughput / throughput at the largest N")
    nc, tc = tp["map_batch_cpu"]
    ng, tg = tp["map_batch_gpu"]
    add("NumMapBatchGpuOverCpuMaxN", times(tg[-1] / tc[-1]),
        "e01 sweep: map_batch GPU / map_batch CPU throughput at the largest N")
    nn, tn = tp["numpy"]
    ne, te = tp["cupy_e2e"]
    add("NumRawGpuOverCpuMaxN", times(te[-1] / tn[-1]),
        "e01 sweep: cupy end-to-end / numpy throughput at the largest N")


def e02() -> None:
    recs = load("e02_gpu_end_to_end")
    for r in recs:
        c = r["config"]
        if c["kind"] == "verify" and c["case"] == "k28":
            s = r["series"]
            add("NumGpuVerifyPointsCpuKTwoEight", grp(s["points_cpu"][-1]),
                "e02 verify k28: CPU unstable points at the last step (12)")
            add("NumGpuVerifyPointsGpuKTwoEight", grp(s["points_gpu"][-1]),
                "e02 verify k28: GPU unstable points at the last step (12)")
            add("NumGpuVerifyMaxCoordKTwoEight", sci(s["max_abs_coord"][-1], 1),
                "e02 verify k28: max |coord| at the last step")
            add("NumGpuFirstDivergentStepKTwoEight", str(r["metrics"]["first_divergent_step"]),
                "e02 verify k28: first step with a non-zero CPU-GPU deviation")
            dev, coord = s["max_abs_deviation"][-2], s["max_abs_coord"][-2]
            add("NumGpuVerifyDevStepElevenKTwoEight", sci(dev, 1),
                "e02 verify k28: max |CPU-GPU| coordinate at step 11")
            add("NumGpuVerifyRelDevStepElevenKTwoEight", sci(dev / coord, 1),
                "e02 verify k28: max |CPU-GPU| / max|coord| at step 11")
        if c["kind"] == "profile" and c["backend"] == "cpu":
            word = {"k28": "KTwoEight", "p3": "PThree", "k10": "KTen"}[c["case"]]
            cats = r["metrics"]["profile_categories_s"]
            tot = sum(cats.values())
            fr = {k: v / tot for k, v in cats.items()}
            add(f"NumProfRefinement{word}", pct(fr["refinement"], 0),
                f"e02 profile {c['case']} cpu: refinement share of tottime")
            add(f"NumProfLinkedList{word}", pct(fr["linked_list"], 0),
                f"e02 profile {c['case']} cpu: linked-list share of tottime")
            add(f"NumProfOther{word}", pct(fr["other"], 0),
                f"e02 profile {c['case']} cpu: other share of tottime")
            if c["case"] == "k28":
                total = r["metrics"]["profile_total_s"]
                top = r["metrics"]["profile_top40"]
                refine = next(t for t in top if t["function"] == "_refine_layer")
                point = next(t for t in top if t["function"] == "__init__" and t["file"].endswith("/Point.py"))
                arr = next(t for t in top if t["function"] == "<built-in method numpy.array>")
                add("NumProfTotalSKTwoEight", f"{total:.0f}", "e02 profile k28 cpu: total profiled seconds")
                add("NumProfRefineLayerShareKTwoEight", pct(refine["tottime"] / total, 0),
                    "e02 profile k28 cpu: ManifoldMachine._refine_layer tottime / total")
                add("NumProfPointInitCallsKTwoEight", sci(point["ncalls"], 2),
                    "e02 profile k28 cpu: Point.__init__ calls")
                add("NumProfPointInitShareKTwoEight", pct(point["cumtime"] / total, 0),
                    "e02 profile k28 cpu: Point.__init__ cumtime / total")
                add("NumProfNpArrayCallsKTwoEight", sci(arr["ncalls"], 2),
                    "e02 profile k28 cpu: numpy.array calls")
            add(f"NumProfRefPlusList{word}", pct(fr["refinement"] + fr["linked_list"], 0),
                f"e02 profile {c['case']} cpu: refinement + linked-list share of tottime")
    # best-step points for k28
    best = None
    for r in recs:
        c = r["config"]
        if c["kind"] == "timed" and c["case"] == "k28":
            best = max(best or 0, max(r["series"]["points_unstable"]))
    if best:
        add("NumGpuTimedMaxPointsKTwoEight", sci(best, 1),
            "e02 timed k28: unstable points after the last growth step")


def e03() -> None:
    recs = load("e03_scaling")
    words = {(2.8, 1e-7): "KTwoEight", (4.0, 1e-7): "KFour", (6.0, 1e-7): "KSix", (10.0, 1e-7): "KTen",
             (2.8, 1e-5): "KTwoEightCoarse"}
    kills = []
    for r in recs:
        c = r["config"]
        if c["kind"] != "growth":
            continue
        w = words[(c["k"], c["area_cutoff"])]
        s = r["series"]
        add(f"NumLastGrowthFactor{w}", fixed(s["growth_factor"][-1], 1),
            f"e03 growth k={c['k']} cutoff {c['area_cutoff']}: point growth factor of the last completed step")
        add(f"NumLambdaU{w}", fixed(r["metrics"]["lambda_u"], 2),
            f"e03 growth k={c['k']}: unstable eigenvalue lambda_u")
        add(f"NumLastMaxCoord{w}", sci(s["max_abs_coord"][-1], 1),
            f"e03 growth k={c['k']} cutoff {c['area_cutoff']}: max|coord| after the last completed step")
        add(f"NumLastPoints{w}", grp(s["points_unstable"][-1]),
            f"e03 growth k={c['k']} cutoff {c['area_cutoff']}: unstable points after the last completed step")
        add(f"NumGrowthSteps{w}", str(len(s["points_unstable"])),
            f"e03 growth k={c['k']} cutoff {c['area_cutoff']}: completed growth steps")
        if r["outcome"] == "oom":
            kills.append((s["points_unstable"][-1], r["process_wall_s"] - r["wall_s"]))
    add("NumScalingOomConfigs", str(len(kills)), "e03 growth: configs killed at the 24 GB RSS cap")
    add("NumScalingOomMinPointsHeld", grp(min(p for p, _ in kills)),
        "e03 growth oom: fewest unstable points held after the last completed step")
    add("NumScalingOomMaxPointsHeld", grp(max(p for p, _ in kills)),
        "e03 growth oom: most unstable points held after the last completed step")
    add("NumScalingOomKilledStepMinS", f"{min(t for _, t in kills):.0f}",
        "e03 growth oom: process wall minus checkpointed wall = time spent in the killed step (min, s)")
    add("NumScalingOomKilledStepMaxS", f"{max(t for _, t in kills):.0f}",
        "e03 growth oom: time spent in the killed step (max, s)")
    bpp = []
    for r in recs:
        if r["config"]["kind"] == "growth" and r["outcome"] == "ok":
            pass
    # ceiling estimate from the median bytes/point of numbers.tex (1058 B/point)
    add("NumPointCeilingTwentyFourGB", sci(24 * 2**30 / 1057.7246839312834, 1),
        "derived: 24 GiB / median RSS slope 1057.7 B/point (summary.json BytesPerPoint)")


def e04() -> None:
    recs = load("e04_area_cutoff")
    runs: dict[tuple, list] = {}
    for r in recs:
        c = r["config"]
        if r["outcome"] != "ok":
            continue
        m = r["metrics"]
        p50 = max(v for k, v in m.items() if k.startswith("inv_") and k.endswith("_p50"))
        runs.setdefault((c["case"], c["unstable_steps"]), []).append((c["area_cutoff"], p50, m["n_crossings"]))
    words = {("k28", 10): "KTwoEightTen", ("k28", 12): "KTwoEightTwelve", ("p3", 13): "PThreeThirteen",
             ("p3", 17): "PThreeSeventeen"}
    for key, v in runs.items():
        v.sort()
        fine = [row for row in v if row[0] <= 1e-6]
        add(f"NumInvSlope{words[key]}", fixed(slope([a for a, *_ in fine], [b for _, b, _ in fine]), 2),
            f"e04 {key}: log-log slope of the p50 invariance error vs cutoff, cutoffs 1e-9..1e-6")
        cross = {a: n for a, _, n in v}
        add(f"NumCrossings{words[key]}Ref", str(cross[1e-9]), f"e04 {key}: crossings at cutoff 1e-9")
        if 1e-3 in cross:
            add(f"NumCrossings{words[key]}CutoffThree", str(cross[1e-3]), f"e04 {key}: crossings at cutoff 1e-3")
        add(f"NumCrossings{words[key]}CutoffFour", str(cross[1e-4]), f"e04 {key}: crossings at cutoff 1e-4")
        add(f"NumCrossings{words[key]}CutoffFive", str(cross[1e-5]), f"e04 {key}: crossings at cutoff 1e-5")
    # crossing-position convergence slopes from the convergence dicts in summary.json
    summ = json.load(open(REPORT / "summary.json"))
    conv = summ["experiments"]["e04_area_cutoff"]["details"]["convergence"]
    for key, rows in conv.items():
        case, steps = key.split("/")
        w = words[(case, int(steps))]
        rows = [row for row in rows if row["median"] and row["cutoff"] <= 1e-5]
        add(f"NumCrossDevSlope{w}", fixed(slope([r["cutoff"] for r in rows], [r["median"] for r in rows]), 2),
            f"e04 {key}: log-log slope of the median crossing distance to the 1e-9 run vs cutoff, cutoffs <= 1e-5")


def e05() -> None:
    recs = load("e05_invariants")
    for r in recs:
        case = r["metrics"].get("case")
        if case == "k28_two_blasts":
            lobes = r["series"]["lobes"]
            a = [row for row in lobes if abs(abs(row[2]) - 7.0) < 0.1]
            add("NumLobeAreaKTwoEight", fixed(abs(a[0][2]), 3),
                "e05 k28 two blasts: |signed area| of the turnstile-sized lobes")
            add("NumLobeAreaKTwoEightCount", str(len(a)),
                "e05 k28 two blasts: lobes of |area| within 0.1 of 7")
            add("NumLobeAreaAbsErrKTwoEight", sci(max(abs(row[3] - row[2]) for row in a), 1),
                "e05 k28 two blasts: max |A_image - A| over those lobes")
        if case == "p3":
            lobes = r["series"]["lobes"]
            areas = [abs(row[2]) for row in lobes]
            add("NumLobeAreaPThreeMedian", sci(float(np.median(areas)), 2),
                "e05 p3: median |signed lobe area|")
            add("NumLobeAbsErrPThreeMax", sci(max(abs(row[3] - row[2]) for row in lobes), 1),
                "e05 p3: max |A_image - A| over the lobes")
            add("NumLobePairsPThree", str(len(lobes)), "e05 p3: lobe pairs compared")
        if case == "k10":
            fp = r["metrics"]["fp"][0]["orbit"][0]
            add("NumDetDefectKTen", sci(abs(fp["det_minus_1"]), 1), "e05 k10: ||lu ls| - 1|")
            add("NumJacDevKTen", sci(fp["jac_vs_analytic"], 1),
                "e05 k10: max |stored - analytic| Jacobian entry (the k10 recipe passes no analytic Jacobian)")
    others = [abs(o["det_minus_1"]) for r in recs if r["metrics"].get("fp") and r["metrics"]["case"] != "k10"
              for f in r["metrics"]["fp"] for o in f["orbit"]]
    add("NumDetDefectMaxAnalytic", sci(max(others), 1),
        "e05: max ||lu ls| - 1| over every case except k10 (analytic Jacobian)")


def e06() -> None:
    recs = load("e06_parameter_sweep")
    sad = [r for r in recs if r["config"]["which"] == "saddle" and r["config"]["b"] == 1.0]
    maxit = [r for r in sad if r["outcome"] == "exception" and "Max iterations" in (r["error"] or "")]
    add("NumSweepMaxIterFails", str(len(maxit)), "e06 b=1 saddle: 'Max iterations reached' exceptions")
    add("NumSweepMaxIterKMin", f"{min(r['config']['k'] for r in maxit):.2f}", "e06: smallest k of those")
    add("NumSweepMaxIterKMax", f"{max(r['config']['k'] for r in maxit):.3f}", "e06: largest k of those")
    ok = [r for r in sad if r["outcome"] == "ok"]
    zero = [r for r in ok if not r["metrics"].get("n_crossings")]
    add("NumSweepSaddleOk", str(len(ok)), "e06 b=1 saddle: configs ok through trellis_pips")
    add("NumSweepZeroCrossings", str(len(zero)), "e06 b=1 saddle ok: configs with 0 non-anchor crossings")
    add("NumSweepZeroCrossingsKMax", f"{max(r['config']['k'] for r in zero):.2f}",
        "e06 b=1 saddle ok: largest k with 0 crossings")
    add("NumSweepFirstCrossingK", f"{min(r['config']['k'] for r in ok if r['metrics'].get('n_crossings')):.2f}",
        "e06 b=1 saddle ok: smallest k with >= 1 crossing")
    add("NumSweepMaxCrossings", str(max(r["metrics"].get("n_crossings") or 0 for r in ok)),
        "e06 b=1 saddle ok: most crossings in any config")
    oth = [r for r in recs if r["config"]["which"] == "other"]
    ell = [r for r in oth if "not a saddle" in (r["error"] or "")]
    oom = [r for r in oth if r["outcome"] == "oom"]
    add("NumSweepOtherElliptic", str(len(ell)), "e06 b=1 other point: 'not a saddle' exceptions")
    add("NumSweepOtherEllipticKMax", f"{max(r['config']['k'] for r in ell):.2f}", "e06: largest k of those")
    add("NumSweepOtherOom", str(len(oom)), "e06 b=1 other point: oom at the 4 GB cap")
    add("NumSweepOtherOomKMin", f"{min(r['config']['k'] for r in oom):.1f}", "e06: smallest k of those")
    neg = [r for r in recs if r["config"]["b"] == -1.0]
    add("NumSweepBMinusOneFails", str(sum(1 for r in neg if r["outcome"] == "exception")),
        "e06 b=-1: exceptions")
    dis = [r for r in recs if r["config"]["b"] not in (1.0, -1.0)]
    add("NumSweepDissipativeOk", f"{sum(1 for r in dis if r['outcome'] == 'ok')}/{len(dis)}",
        "e06 b in {0.3,0.5,0.9}: ok / configs")


def e07() -> None:
    recs = load("e07_breaking")
    for r in recs:
        c = r["config"]
        if c["study"] == "escape" and c["mode"] == "fixed":
            add("NumEscapeKilledStepS", f"{r['process_wall_s'] - r['wall_s']:.0f}",
                "e07 escape fixed: process wall minus checkpointed wall (time inside the killed step, s)")
        if c["study"] == "depth" and c["case"] == "inversion" and c["unstable_steps"] == 7:
            add("NumInversionPointsSixSteps", grp(r["series"]["unstable_points"][-1]),
                "e07 depth inversion 7: unstable points after step 6, before the kill")
        if c["study"] == "depth" and c["case"] == "p3" and c["unstable_steps"] == 18:
            add("NumPThreePointsSeventeen", grp(r["series"]["unstable_points"][-1]),
                "e07 depth p3 18: unstable points after step 17, before the kill")
            add("NumPThreeStepSeventeenS", f"{r['series']['step_s'][-1]:.1f}",
                "e07 depth p3 18: time of growth step 17 (s)")
        if c["study"] == "blasts":
            w = {"k28": "KTwoEight", "p3": "PThree", "nested": "Nested"}[c["case"]]
            s, m = r["series"], r["metrics"]
            add(f"NumBlastCrossingsStart{w}", grp(s["crossings"][0]), f"e07 blasts {c['case']}: crossings before blasting")
            add(f"NumBlastCrossingsEnd{w}", grp(s["crossings"][-1]), f"e07 blasts {c['case']}: crossings after the last blast")
            add(f"NumBlastLastS{w}", f"{s['blast_s'][-1]:.2g}", f"e07 blasts {c['case']}: time of the last blast (s)")
            add(f"NumBlastClasses{w}", str(m["n_classes"]), f"e07 blasts {c['case']}: bridge classes after the last blast")
            add(f"NumBlastUnresolved{w}", str(m["n_unresolved"]), f"e07 blasts {c['case']}: unresolved classes")
            add(f"NumBlastVirtual{w}", str(m["n_virtual"]), f"e07 blasts {c['case']}: virtual classes")
    summ = json.load(open(REPORT / "summary.json"))["experiments"]["e07_breaking"]["details"]
    p3 = summ["depth_p3"]
    rel = [int(n) for n, d in p3.items() if d["reliable"]]
    unrel = [int(n) for n, d in p3.items() if d["reliable"] is False]
    add("NumDepthReliableMaxPThree", str(max(rel)), "e07 depth p3: deepest run with is_reliable = True")
    add("NumDepthUnreliableMinPThree", str(min(unrel)), "e07 depth p3: shallowest run with is_reliable = False")
    add("NumDepthOomPThree", " and ".join(str(n) for n in sorted(int(n) for n, d in p3.items() if d["outcome"] == "oom")),
        "e07 depth p3: unstable steps whose run hit the RSS cap")
    inv = summ["depth_inversion"]
    add("NumInversionAssertSteps", " and ".join(str(n) for n in sorted(int(n) for n, d in inv.items()
                                                                      if d["topology"] == "invariant_assert")),
        "e07 depth inversion: unstable steps whose topology stack raised an invariant AssertionError")
    add("NumInversionOomSteps", " and ".join(str(n) for n in sorted(int(n) for n, d in inv.items() if d["outcome"] == "oom")),
        "e07 depth inversion: unstable steps whose run hit the RSS cap")
    tan = [r for r in recs if r["config"]["study"] == "tangency"]
    add("NumTangencyConfigs", str(len(tan)), "e07 tangency: configs (24 k values x 3 cutoffs)")
    by = {}
    for r in tan:
        by[(r["config"]["area_cutoff"], round(r["config"]["k"] + 1, 6))] = r["metrics"]
    names = {1e-7: "Seven", 1e-10: "Ten", 1e-12: "Twelve"}
    for cut, w in names.items():
        add(f"NumTangencyCrossingsNearCutoff{w}", grp(by[(cut, 0.0001)]["non_anchor_crossings"]),
            f"e07 tangency k+1=1e-4 cutoff {cut}: non-anchor crossings")
        peak = max((m["non_anchor_crossings"], kk) for (cc, kk), m in by.items() if cc == cut)
        add(f"NumTangencyCrossingsPeakCutoff{w}", grp(peak[0]),
            f"e07 tangency cutoff {cut}: most non-anchor crossings over k")
    far = sorted({kk for (_, kk) in by if kk >= 0.158})
    same = all(len({by[(c, kk)]["non_anchor_crossings"] for c in names}) == 1 for kk in far)
    add("NumTangencyAgreeKPlusOne", "0.16" if same else "n/a",
        "e07 tangency: smallest k+1 from which all three cutoffs give the same crossing count")
    add("NumTangencyCrossingsFar", str(by[(1e-7, 1.0)]["non_anchor_crossings"]),
        "e07 tangency k+1 = 1: non-anchor crossings (every cutoff)")
    add("NumTangencyRefNear", sci(by[(1e-7, 0.0001)]["splitting_reference"], 1),
        "e07 tangency k+1=1e-4: heuristic splitting reference exp(-pi^2/ln lambda)")
    add("NumTangencyRefMid", sci(by[(1e-7, 0.025119)]["splitting_reference"], 1),
        "e07 tangency k+1=0.025: heuristic splitting reference")
    near = [m for (c, kk), m in by.items() if c == 1e-12 and kk <= 0.026]
    add("NumTangencyMinSinFineLo", sci(min(m["min_sin"] for m in near), 1),
        "e07 tangency cutoff 1e-12, k+1 <= 0.025: smallest per-config min |sin|")
    add("NumTangencyMinSinFineHi", sci(max(m["min_sin"] for m in near), 1),
        "e07 tangency cutoff 1e-12, k+1 <= 0.025: largest per-config min |sin|")
    add("NumTangencyMedianSinFine", sci(float(np.median([m["median_sin"] for m in near])), 1),
        "e07 tangency cutoff 1e-12, k+1 <= 0.025: median over configs of the median |sin|")
    add("NumTangencyWarningsTotal", str(sum(sum((r["warnings"] or {}).values()) for r in tan)),
        "e07 tangency: WARNING records of any kind, all configs")


def e08() -> None:
    recs = load("e08_fixed_point_solver")
    conv = none = slide = accepted = 0
    k10_slide = k10_conv = 0
    times_s = []
    for r in recs:
        p, k = r["config"]["period"], r["config"]["k"]
        for s in r["series"]["solves"]:
            times_s.append(s["solve_s"])
            mp = s.get("minimal_period")
            if mp is None:
                continue
            conv += 1
            if k == 10:
                k10_conv += 1
            if mp < p:
                if p == 2 and k < 3:
                    none += 1
                else:
                    slide += 1
                    if k == 10:
                        k10_slide += 1
                if s["outcome"] == "saddle":
                    accepted += 1
    add("NumSolverConverged", grp(conv), "e08: solves with a converged orbit (minimal period measured)")
    add("NumSolverLowerNoOrbit", grp(none), "e08: lower-period solves at requested period 2 with k < 3 (no period-2 orbit)")
    add("NumSolverLowerNoOrbitRate", pct(none / conv), "e08: those / converged")
    add("NumSolverLowerSlide", grp(slide), "e08: lower-period solves where the requested period has orbits")
    add("NumSolverLowerSlideRate", pct(slide / conv), "e08: those / converged")
    add("NumSolverLowerSlideRateKTen", pct(k10_slide / k10_conv), "e08 k = 10: lower-period / converged")
    add("NumSolverLowerAccepted", grp(accepted),
        "e08: lower-period solves that construct_fixed_point accepted as a saddle of the requested period")
    add("NumSolverMedianMs", f"{1e3 * float(np.median(times_s)):.2f}", "e08: median solve time (ms)")
    add("NumSolverPNinetyNineMs", f"{1e3 * float(np.percentile(times_s, 99)):.1f}", "e08: p99 solve time (ms)")


def e09() -> None:
    recs = load("e09_intersections")
    faster = sum(1 for r in recs if r["metrics"]["brute_s"] < r["metrics"]["library_s"])
    add("NumIsectBruteFasterConfigs", str(faster), "e09: configs where the brute force beat the library")
    big = max(recs, key=lambda r: r["metrics"]["unstable_segments"] + r["metrics"]["stable_segments"])
    m = big["metrics"]
    add("NumIsectLargestSegments", sci(m["unstable_segments"] + m["stable_segments"], 1),
        "e09: segments of the largest config")
    add("NumIsectLargestLibraryS", f"{m['library_s']:.1f}", "e09 largest config: library time (s)")
    add("NumIsectLargestBruteS", f"{m['brute_s']:.1f}", "e09 largest config: brute-force time (s)")
    uu = [r for r in recs if r["metrics"].get("uu_crossings")]
    add("NumIsectUUConfigs", str(len(uu)), "e09: configs with at least one u x u crossing")
    add("NumIsectUUCutoff", sci(uu[0]["config"]["area_cutoff"], 1) if all(
        r["config"].get("area_cutoff") == uu[0]["config"].get("area_cutoff") for r in uu) else "mixed",
        "e09: area_cutoff of every config with u x u crossings")


def e10() -> None:
    recs = load("e10_determinism")
    tot: dict[str, list] = {}
    for r in recs:
        tot.setdefault(r["config"]["device"], []).append(sum(r["stages"].values()))
    c, g = float(np.median(tot["cpu"])), float(np.median(tot["gpu"]))
    add("NumDetTotalCpuS", f"{c:.3f}", "e10: median summed stage time, CPU (s)")
    add("NumDetTotalGpuS", f"{g:.3f}", "e10: median summed stage time, GPU (s)")
    add("NumDetGpuSpeedup", times(c / g), "e10: CPU / GPU median summed stage time")


def campaign() -> None:
    names = ["e01_map_throughput", "e02_gpu_end_to_end", "e03_scaling", "e04_area_cutoff", "e05_invariants",
             "e06_parameter_sweep", "e07_breaking", "e08_fixed_point_solver", "e09_intersections",
             "e10_determinism", "e11_blast_depth"]
    out: Counter = Counter()
    n = 0
    for name in names:
        recs = load(name)
        n += len(recs)
        out.update(r["outcome"] for r in recs)
        digits = {"01": "One", "02": "Two", "03": "Three", "04": "Four", "05": "Five", "06": "Six",
                  "07": "Seven", "08": "Eight", "09": "Nine", "10": "Ten", "11": "Eleven"}[name[1:3]]
        add(f"NumConfigsE{digits}", str(len(recs)), f"{name}.jsonl: records")
    add("NumConfigsTotal", grp(n), "all full-run JSONL files (e01-e11): records")
    add("NumConfigsCampaign", grp(n - len(load("e11_blast_depth"))),
        "e01-e10 full-run JSONL files: records of the scripted campaign (e11 ran on its own)")
    import re
    done = re.findall(r"campaign done in ([0-9.]+) min", (RESULTS / "progress.log").read_text())
    add("NumCampaignQuickMin", done[-2], "progress.log: wall-clock of the last --quick campaign (min)")
    add("NumCampaignFullMin", done[-1], "progress.log: wall-clock of the full campaign e01-e10 (min; e11 ran on its own, see NumDeepWallMin)")
    for o in ("ok", "exception", "oom", "timeout", "invariant_assert"):
        add(f"NumOutcomeTotal{o.replace('_', ' ').title().replace(' ', '')}", str(out.get(o, 0)),
            f"all full-run JSONL files: outcome {o}")


CAPS = {  # (timeout s, RSS cap GB, workers): the scripts' TIMEOUT / RSS_CAP_GB / WORKERS constants
    "e01_map_throughput": (300, 24), "e02_gpu_end_to_end": (600, 24), "e03_scaling": (900, 24),
    "e04_area_cutoff": (900, 16), "e05_invariants": (600, 12), "e06_parameter_sweep": (240, 4),
    "e07_breaking": (600, 20), "e08_fixed_point_solver": (900, 24), "e09_intersections": (600, 16),
    "e10_determinism": (300, 24), "e11_blast_depth": (900, 12),
}


def counts_table() -> None:
    """tables/campaign_counts.tex: records and outcomes per experiment."""
    rows = []
    tot: Counter = Counter()
    for name, (timeout, cap) in CAPS.items():
        recs = load(name)
        out = Counter(r["outcome"] for r in recs)
        tot.update(out)
        tot["n"] += len(recs)
        label = name.replace("_", "\_")
        rows.append(f"{label} & {len(recs)} & {out.get('ok', 0)} & {out.get('exception', 0)} & "
                    f"{out.get('invariant_assert', 0)} & {out.get('oom', 0)} & {out.get('timeout', 0)} & "
                    f"{timeout} & {cap} \\\\")
    rows.append("\\midrule")
    rows.append(f"total & {tot['n']} & {tot['ok']} & {tot['exception']} & {tot['invariant_assert']} & "
                f"{tot['oom']} & {tot['timeout']} & & \\\\")
    body = ["% campaign: records and outcomes per experiment (full run); timeout / RSS cap = default per config",
            "% suggested column spec: {lrrrrrrrr}", "\\toprule",
            "experiment & configs & ok & exception & assert & oom & timeout & timeout (s) & cap (GB) \\\\",
            "\\midrule", *rows, "\\bottomrule"]
    (REPORT / "tables" / "campaign_counts.tex").write_text("\n".join(body) + "\n")


def main() -> None:
    counts_table()
    for fn in (e01, e02, e03, e04, e05, e06, e07, e08, e09, e10, campaign):
        fn()
    lines = ["% Generated by stress/numerics/report/make_numbers_extra.py -- do not edit.",
             "% Computed from the full-run results/*.jsonl (and summary.json for e04 convergence).", ""]
    seen = set()
    for name, value, source in MACROS:
        assert name not in seen, name
        assert name.isalpha(), name
        seen.add(name)
        lines.append(f"% source: {source}")
        lines.append(f"\\newcommand{{\\{name}}}{{{value}}}")
    OUT.write_text("\n".join(lines) + "\n")
    print(f"wrote {len(MACROS)} macros to {OUT}")


if __name__ == "__main__":
    main()
