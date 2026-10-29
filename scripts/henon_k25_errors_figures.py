"""
The figure half of ``henon_k25_errors.py``: one multi-panel error analysis against ``area_cutoff``.

Reads the error records (``errors/data/``) and the pruned and unpruned entropy
sweep records of ``figures/henon_k25_entropy_post_array_refactor/``, and draws
``errors/errors.png`` (+ .svg, .fig.pickle):

* (a, b) the registered pip orbit's distance from P against j, one line per
  cutoff, with the ``beta^-j`` slope, the float64 floor and the next point
  (hollow) the sweep stalled on;
* (c) the stall's anatomy against the cutoff: at the first blast that left the
  depth unchanged, the unstable chord deviation at ``f^(m+1)(q0)``, the stable
  segment there (times the crossing's sine) and the point's distance from P;
* (d) the stall's cause: the distance from ``f^(m+1)(q0)`` to the nearest new
  image segment in units of its chord deviation (near 1: the curve is there and
  the crossing failed the ratio, annotated; far above: the bridge carrying
  ``f^m(q0)`` was never iterated);
* (e) the depth reached against the cutoff, with the fitted levels per decade
  and the float64 depth, and the trustworthy depth (:func:`trustworthy_depth`);
* (f) chord deviation of the blasts' image curves against the cutoff, beside
  ``min_separation`` and the ``A^(2/3)`` and ``A`` slopes;
* (g) every detection's resolution ratio at the end of the sweep;
* (h) the smallest relative cdist step of each blast's image curves;
* (i) unresolved detections and same-stability polyline crossings per sweep;
* (j) the prune: active classes, pruned against unpruned, at matched blasts;
* (k) the cost: points and peak memory against the cutoff;
* (l) the registered orbit's contraction per step over beta (p3: per period
  over lambda): 1 while the orbit points are the real ones, scattered once they
  lie closer to P than the curves' chord deviation.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import henon_k25_entropy as entropy  # noqa: E402
import henon_k25_errors as errors  # noqa: E402
from symbolic_figures import save_panel  # noqa: E402

THRESHOLD = 3.0
BETA = {"p1": 5.561861466649672, "p3": 1.9132631930650705}


def load(directory: Path) -> dict[tuple[str, float], dict]:
    """Every record in a folder, keyed by ``(case, cutoff)``."""
    records = {}
    for path in sorted(directory.glob("*.json")):
        record = json.loads(path.read_text())
        if record["rows"]:
            records[(record["case"], record["area_cutoff"])] = record
    return records


def stall(record: dict, period: str):
    """The first row whose blast left ``period``'s depth unchanged after it had advanced, and its orbit block."""
    advanced = False
    for row in record["rows"]:
        orbit = row["orbits"].get(period)
        if not orbit or not orbit["orbit"]:
            continue
        if orbit["advanced"]:
            advanced = True
        elif advanced and "next" in orbit:
            return row, orbit
    return None, None


def depth_of(record: dict, period: str) -> int:
    """The final depth of ``period``'s pip orbit."""
    return max((len(row["orbits"][period]["orbit"]) - 1 for row in record["rows"] if period in row["orbits"]),
               default=0)


def trustworthy_depth(record: dict, period: str, linear: float = 1e-3, band: tuple = (0.8, 1.25)) -> int:
    """
    The deepest j up to which the registered orbit contracts like the real one.

    Inside the linear regime (``|f^j - P| < linear``) each step must divide the
    distance to P by beta (p3: each period by lambda) to within ``band``; the
    first step that does not ends the trustworthy orbit.
    """
    orbit = max((row["orbits"][period] for row in record["rows"]), key=lambda o: len(o["orbit"]))["orbit"]
    dist = [p["distance"] for p in orbit]
    step = 3 if period == "p3" else 1
    for j in range(step, len(dist)):
        if dist[j - step] < linear and not band[0] < dist[j - step] / dist[j] / BETA[period] ** step < band[1]:
            return j - 1
    return len(dist) - 1


def periods(case: str) -> list[str]:
    """The fixed points of a case."""
    return ["p1", "p3"] if case == "nested" else [case]


def draw() -> Path:
    """Draw the error figure from the records on disk."""
    records = load(errors.DATA_DIR)
    pruned = load(errors.RUN_DIR / "data")
    unpruned = load(errors.RUN_DIR / "components" / "unpruned" / "data")
    cutoffs = sorted({c for _, c in records} | {c for _, c in pruned})
    styles = entropy.cutoff_styles(cutoffs)
    fig, axes = plt.subplots(3, 4, figsize=(26, 17))
    a, b, c, d, e, f, g, h, i, j, k, l = axes.flat

    # (a, b) the pip orbit
    for ax, period in ((a, "p1"), (b, "p3")):
        for (case, cutoff), record in sorted(records.items()):
            if case != period:
                continue
            last = max((row["orbits"][period] for row in record["rows"]), key=lambda o: len(o["orbit"]))
            x = [p["j"] for p in last["orbit"]]
            y = [p["distance"] for p in last["orbit"]]
            st = styles[cutoff]
            ax.semilogy(x, y, color=st["color"], marker=st["marker"], markersize=4, label=f"cutoff {cutoff:g}")
            _, orbit = stall(record, period)
            if orbit is not None:
                ax.semilogy([len(orbit["orbit"])], [orbit["next"]["distance"]], marker=st["marker"],
                            markerfacecolor="white", color=st["color"], markersize=8, linestyle="")
        ref = [row for (case, _), rec in records.items() if case == period for row in rec["rows"]]
        if ref:
            P = max(o["next"]["P_norm"] for row in ref for o in [row["orbits"][period]] if "next" in o)
            ax.axhline(errors.EPS * P, color="0.3", linestyle="--", linewidth=1)
            ax.annotate("float64 floor  eps·|P|", (0, errors.EPS * P), xytext=(2, 3), textcoords="offset points",
                        fontsize=8, color="0.3")
            x0 = np.arange(0, 45)
            first = [p for rec in [rec for (cs, _), rec in records.items() if cs == period]
                     for p in max((r["orbits"][period] for r in rec["rows"]), key=lambda o: len(o["orbit"]))["orbit"]
                     if p["j"] == 3]
            if first:
                ax.semilogy(x0, first[0]["distance"] * BETA[period] ** (3.0 - x0), color="0.6", linestyle=":",
                            label=rf"$\beta^{{-j}}$, $\beta$ = {BETA[period]:.2f}")
        ax.set_ylim(1e-17, 10)
        ax.set_xlim(left=0)
        ax.set_title(f"({'a' if period == 'p1' else 'b'}) {period} alone: |f^j(q0) − P| (hollow: the next "
                     "point, never in N)", fontsize=10)
        ax.set_xlabel("j")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, ncol=2)

    # (c, d) the stall
    case_marker = {"p1": "o", "p3": "s", "nested": "^"}
    rows_cd = []
    for (case, cutoff), record in records.items():
        for period in periods(case):
            row, orbit = stall(record, period)
            if orbit is None or "length" not in orbit["next"]["unstable"]:
                continue
            nxt = orbit["next"]
            rows_cd.append((case, period, cutoff, nxt))
    for case, period, cutoff, nxt in rows_cd:
        color = entropy.CASE_COLORS[case] if case != "nested" else ("#2ca02c" if period == "p1" else "#98df8a")
        m = case_marker[case]
        c.loglog(cutoff, nxt["unstable"]["deviation"], marker=m, color=color, linestyle="")
        c.loglog(cutoff, nxt["stable"]["length"] * nxt.get("sin", 1.0), marker=m, markerfacecolor="white",
                 color=color, linestyle="")
        c.loglog(cutoff, nxt["distance"], marker="x", color=color, linestyle="")
        reach = nxt["unstable"]["distance"] / nxt["unstable"]["deviation"]
        made = reach < 10
        d.loglog(cutoff, reach, marker=m, color=color, linestyle="",
                 markerfacecolor=color if made else "white")
        if made:
            d.annotate(f"ratio {nxt.get('ratio_estimate', np.nan):.2g}", (cutoff, reach), xytext=(4, -3),
                       textcoords="offset points", fontsize=6, color=color)
    for case in entropy.CASES:
        c.plot([], [], marker=case_marker[case], color=entropy.CASE_COLORS[case], linestyle="", label=case)
    c.plot([], [], marker="o", color="0.3", linestyle="", label=r"unstable chord deviation $h_u$")
    c.plot([], [], marker="o", markerfacecolor="white", color="0.3", linestyle="",
           label=r"stable segment $l_s \cdot |\sin\theta|$")
    c.plot([], [], marker="x", color="0.3", linestyle="", label=r"$|f^{m+1}(q_0) - P|$")
    c.set_title("(c) the stall at f^(m+1)(q0): resolved iff $l_s|\\sin| / h_u \\geq 3$", fontsize=10)
    c.set_xlabel("area cutoff")
    c.invert_xaxis()
    c.grid(alpha=0.3)
    c.legend(fontsize=7)
    d.axhline(10, color="k", linestyle="--", linewidth=1)
    d.annotate("curve reaches f^(m+1): stalled on the ratio (< 3)", (cutoffs[0], 1.2), fontsize=7)
    d.annotate("no image curve near f^(m+1): its parent was never iterated", (cutoffs[0], 30), fontsize=7)
    d.set_title("(d) the stall's cause: nearest new image segment to f^(m+1)(q0), in units of $h_u$", fontsize=10)
    d.set_xlabel("area cutoff")
    d.invert_xaxis()
    d.grid(alpha=0.3)

    # (e) depth against cutoff
    for case in entropy.CASES:
        for period in periods(case):
            pts = sorted((cutoff, depth_of(rec, period)) for (cs, cutoff), rec in records.items() if cs == case)
            pts += sorted((cutoff, max(r["pip_depth"].get(period, 0) for r in rec["rows"]))
                          for (cs, cutoff), rec in pruned.items() if cs == case and (case, cutoff) not in records
                          and (rec["stopped"] or "").startswith("frontier exhausted"))
            pts.sort()
            if len(pts) < 2:
                continue
            x, y = np.array(pts).T
            style = "-" if case != "nested" else "--"
            e.semilogx(x, y, marker=case_marker[case], linestyle=style, color=entropy.CASE_COLORS[period],
                       label=f"{case}: {period}")
            fit = (x >= 1e-12) & (x <= 1e-7) if period == "p1" else (x <= 1e-9)
            if case != "nested" and fit.sum() >= 2:
                slope = np.polyfit(-np.log10(x[fit]), y[fit], 1)[0]
                alpha = slope * math.log(BETA[period]) / math.log(10)
                e.annotate(f"{slope:.2f} levels/decade (α = {alpha:.2f})", (x[fit][0], y[fit][0]),
                           xytext=(5, -12), textcoords="offset points", fontsize=8,
                           color=entropy.CASE_COLORS[period])
    for case in entropy.CASES:
        for period in periods(case):
            pts = sorted((cutoff, trustworthy_depth(rec, period)) for (cs, cutoff), rec in records.items() if cs == case)
            if pts:
                x, y = np.array(pts).T
                e.semilogx(x, y, marker=case_marker[case], markerfacecolor="white", linestyle=":",
                           color=entropy.CASE_COLORS[period], alpha=0.8,
                           label=f"{case}: {period} trustworthy" if case != "nested" else None)
    for period in ("p1", "p3"):
        ref = [o for (cs, _), rec in records.items() if cs == period for row in rec["rows"]
               for o in [row["orbits"][period]] if o["orbit"]]
        if ref:
            p = ref[0]["orbit"][min(3, len(ref[0]["orbit"]) - 1)]
            P = ref[0].get("next", {}).get("P_norm", 1.0)
            floor = p["j"] + math.log(p["distance"] / (errors.EPS * P)) / math.log(BETA[period])
            e.axhline(floor, color=entropy.CASE_COLORS[period], linestyle=":", linewidth=1)
            e.annotate(f"{period} float64 depth ≈ {floor:.0f}", (cutoffs[-1], floor), xytext=(0, 3),
                       textcoords="offset points", fontsize=8, color=entropy.CASE_COLORS[period])
    e.invert_xaxis()
    e.set_title("(e) pip depth m reached (filled) and trustworthy (hollow: the orbit still contracts by $\\beta$)", fontsize=10)
    e.set_xlabel("area cutoff")
    e.set_ylabel("m")
    e.grid(alpha=0.3)
    e.legend(fontsize=7)

    # (f) chord deviation of the image curves
    for case in entropy.CASES:
        for q, ls in (("p50", "-"), ("p99", "--"), ("max", ":")):
            pts = []
            for (cs, cutoff), rec in records.items():
                if cs == case:
                    vals = [z["curves"]["deviation"].get(q) for row in rec["rows"] for z in row["zones"].values()]
                    vals = [v for v in vals if v]
                    if vals:
                        pts.append((cutoff, np.median(vals) if q != "max" else max(vals)))
            if pts:
                x, y = np.array(sorted(pts)).T
                f.loglog(x, y, linestyle=ls, marker=case_marker[case], color=entropy.CASE_COLORS[case],
                         label=f"{case} {q}" if case == "p1" or q == "p50" else None)
    xs = np.array([min(cutoffs), max(cutoffs)])
    f.loglog(xs, 3e-6 * (xs / 1e-7) ** (2 / 3), color="0.6", linestyle=":", label=r"$\propto A^{2/3}$")
    f.loglog(xs, 3e-6 * (xs / 1e-7), color="0.8", linestyle="-.", label=r"$\propto A$")
    f.axhline(entropy.MIN_SEPARATION, color="k", linestyle="--", linewidth=1)
    f.annotate("min_separation", (xs[1], entropy.MIN_SEPARATION), xytext=(3, 3), textcoords="offset points",
               fontsize=8)
    f.invert_xaxis()
    f.set_title("(f) chord deviation of the blasts' image curves (median over blasts)", fontsize=10)
    f.set_xlabel("area cutoff")
    f.grid(alpha=0.3)
    f.legend(fontsize=7, ncol=2)

    # (g) resolution ratios
    for (case, cutoff), rec in sorted(records.items()):
        if case != "p1" or "histograms" not in rec:
            continue
        ratios = np.array(rec["histograms"]["ratios"], dtype=float)
        ratios = ratios[ratios > 0]
        if ratios.size:
            g.hist(np.log10(ratios), bins=np.linspace(-2, 8, 101), histtype="step", color=styles[cutoff]["color"],
                   label=f"{cutoff:g}: {np.sum(ratios < THRESHOLD)} below 3")
    g.axvline(math.log10(THRESHOLD), color="k", linestyle="--", linewidth=1)
    g.set_yscale("log")
    g.set_title("(g) p1: every detection's resolution ratio at the end of the sweep", fontsize=10)
    g.set_xlabel(r"$\log_{10}$ ratio")
    g.grid(alpha=0.3)
    g.legend(fontsize=7)

    # (h) cdist precision
    for (case, cutoff), rec in sorted(records.items()):
        if case not in ("p1", "p3"):
            continue
        x = [row["blasts"] for row in rec["rows"]]
        y = [min((z["curves"]["relative_step"].get("min", np.nan) for z in row["zones"].values()), default=np.nan)
             for row in rec["rows"]]
        h.semilogy(x, y, color=styles[cutoff]["color"], linestyle="-" if case == "p1" else "--",
                   label=f"{case} {cutoff:g}" if case == "p1" else None)
    h.axhline(1e-9, color="0.4", linestyle=":", linewidth=1)
    h.annotate("DEDUP_RTOL", (0, 1e-9), xytext=(3, 3), textcoords="offset points", fontsize=8)
    h.axhline(errors.EPS, color="k", linestyle="--", linewidth=1)
    h.annotate("float64 eps", (0, errors.EPS), xytext=(3, 3), textcoords="offset points", fontsize=8)
    h.set_title("(h) smallest relative cdist step dc/c of each blast's images (solid p1, dashed p3)", fontsize=10)
    h.set_xlabel("blasts")
    h.grid(alpha=0.3)
    h.legend(fontsize=7, ncol=2)

    # (i) unresolved and same-stability
    for case in entropy.CASES:
        pts = sorted((cutoff, rec["rows"][-1]["unresolved"], sum(r["same_stability"] for r in rec["rows"]))
                     for (cs, cutoff), rec in records.items() if cs == case)
        if pts:
            x, un, ss = np.array(pts).T
            i.semilogx(x, un, marker=case_marker[case], color=entropy.CASE_COLORS[case], label=f"{case} unresolved")
            i.semilogx(x, ss, marker=case_marker[case], color=entropy.CASE_COLORS[case], linestyle="--",
                       markerfacecolor="white", label=f"{case} same-stability (sum over blasts)")
    i.set_yscale("symlog", linthresh=1)
    i.invert_xaxis()
    i.set_title("(i) unresolved detections (end) and same-stability crossings in the tangle box", fontsize=10)
    i.set_xlabel("area cutoff")
    i.grid(alpha=0.3)
    i.legend(fontsize=7)

    # (j) the prune
    for (case, cutoff), rec in sorted(unpruned.items()):
        if (case, cutoff) not in pruned:
            continue
        mine = {r["blasts"]: r for r in pruned[(case, cutoff)]["rows"]}
        rows = [r for r in rec["rows"] if r["blasts"] in mine]
        x = [r["blasts"] for r in rows]
        color = styles[cutoff]["color"]
        last = rows[-1]
        h_text = (f"h {mine[last['blasts']]['entropy']:.7f} vs {last['entropy']:.7f}"
                  if abs(last["entropy"] - mine[last["blasts"]]["entropy"]) > 1e-12 else "h equal")
        j.plot(x, [mine[n]["active"] for n in x], color=color, marker=case_marker[case], markersize=4,
               label=f"{case} {cutoff:g}: {h_text} (blast {last['blasts']})")
        j.plot(x, [r["active"] for r in rows], color=color, linestyle="--", linewidth=2.5, alpha=0.5)
    j.set_title("(j) active classes, min_separation 1e-5 (solid) against none (dashed)", fontsize=10)
    j.set_xlabel("blasts")
    j.grid(alpha=0.3)
    j.legend(fontsize=7, title="h pruned vs unpruned at the last shared blast", title_fontsize=7)

    # (k) cost
    for case in entropy.CASES:
        pts = sorted((cutoff, rec["rows"][-1]["points"], max(r["peak_gb"] for r in rec["rows"]))
                     for (cs, cutoff), rec in pruned.items() if cs == case)
        if pts:
            x, pt, gb = np.array(pts).T
            k.loglog(x, pt / 1e6, marker=case_marker[case], color=entropy.CASE_COLORS[case], label=f"{case} Mpoints")
            k.loglog(x, gb, marker=case_marker[case], color=entropy.CASE_COLORS[case], linestyle="--",
                     markerfacecolor="white", label=f"{case} peak GB")
    k.invert_xaxis()
    k.set_title("(k) cost of the pruned sweep, at its last blast", fontsize=10)
    k.set_xlabel("area cutoff")
    k.grid(alpha=0.3)
    k.legend(fontsize=7)

    # (l) the orbit's own contraction: |f^j - P| / |f^(j+1) - P| per step, against the eigenvalue
    for (case, cutoff), rec in sorted(records.items()):
        if case not in ("p1", "p3"):
            continue
        orbit = max((row["orbits"][case] for row in rec["rows"]), key=lambda o: len(o["orbit"]))["orbit"]
        dist = np.array([p["distance"] for p in orbit])
        if case == "p3":  # one period: three steps, so compare with the full return lambda
            x = np.arange(3, len(dist))
            y = dist[x - 3] / dist[x] / BETA["p3"] ** 3
        else:
            x = np.arange(1, len(dist))
            y = dist[x - 1] / dist[x] / BETA["p1"]
        l.semilogy(x, y, color=styles[cutoff]["color"], marker=case_marker[case], markersize=4,
                   linestyle="-" if case == "p1" else "--", label=f"{case} {cutoff:g}")
    l.axhline(1.0, color="k", linewidth=1)
    l.axhspan(0.8, 1.25, color="0.85", zorder=0)
    l.set_title("(l) is the registered orbit the real one? contraction per step / $\\beta$ "
                "(p3: per period / $\\lambda$)", fontsize=10)
    l.set_xlabel("j")
    l.grid(alpha=0.3)
    l.legend(fontsize=6, ncol=3)

    fig.suptitle(f"Hénon k={entropy.K:g}: numerical errors of the entropy sweeps against the area cutoff "
                 f"(min_separation {entropy.MIN_SEPARATION:g})", fontsize=14)
    fig.tight_layout()
    return save_panel(fig, errors.OUT_DIR, "errors")


if __name__ == "__main__":
    plt.switch_backend("Agg")
    print(draw())
