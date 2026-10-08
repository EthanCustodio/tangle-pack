"""Figures of the numerics stress-test report.

    env/bin/python stress/numerics/report/make_figures.py [--quick] [--only e04,e07]

Writes one vector PDF per figure into ``report/figures/`` and a 110 dpi PNG
preview into ``report/figures/png/``. Every figure is one ``fig_eNN_<what>``
function reading ``analysis``; a figure with no usable data is skipped and a
panel with none says so, so a partial campaign still yields every figure it
can. Without ``--quick`` each experiment uses its full results and falls back
to the quick ones when the full file is missing or empty.
"""

from __future__ import annotations

import argparse
import logging
import math
import sys
import traceback
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analysis as A  # noqa: E402

REPORT = Path(__file__).resolve().parent
FIGURES = REPORT / "figures"
PNG = FIGURES / "png"

# The categorical palette, in its fixed order; never cycled.
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
INK, INK2, GRID, LIGHT = "#0b0b0b", "#52514e", "#e6e5e1", "#c9c7c1"
RAMP = LinearSegmentedColormap.from_list("blue_ramp", ["#cde2fb", "#0d366b"])

CPU, GPU = BLUE, ORANGE
CASE_COLORS = {"k28": BLUE, "p3": AQUA, "k10": VIOLET, "nested": MAGENTA, "inversion": YELLOW,
               "k32": ORANGE, "nested_outer": MAGENTA, "nested_inner": GREEN}
STATUS_COLORS = {"ok": GREEN, "exception": RED, "timeout": YELLOW, "oom": VIOLET,
                 "invariant_assert": MAGENTA, "crash": INK}
STATUS_LABELS = {"ok": "ok", "exception": "exception", "timeout": "timeout", "oom": "out of memory",
                 "invariant_assert": "invariant assert", "crash": "crash"}
# Pipeline stages; growth in the blue ramp, the later stages after the cases' colours.
STAGE_COLORS = {"grow_unstable": "#0d366b", "grow_stable": "#7fb0ee", "compute_intersections": GREEN,
                "trim_and_bridges": RED, "infer_iterate_table": VIOLET}
STAGE_LABELS = {"grow_unstable": "grow unstable", "grow_stable": "grow stable",
                "compute_intersections": "intersections", "trim_and_bridges": "trim + bridges",
                "infer_iterate_table": "iterate table"}
WIDE, TALL = (6.3, 2.6), (6.3, 4.2)
MARKERS = ("o", "s", "^", "D", "v", "P")

SKIPPED: list[str] = []


class Skip(Exception):
    """A figure with nothing to draw."""


def style() -> None:
    """Print rcParams: serif to match the LaTeX body, 9 pt, thin axes, light grid."""
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["Latin Modern Roman", "CMU Serif", "DejaVu Serif"],
        "mathtext.fontset": "cm", "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
        "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": INK2, "ytick.color": INK2, "xtick.labelcolor": INK, "ytick.labelcolor": INK,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4, "xtick.major.size": 3, "ytick.major.size": 3,
        "axes.grid": True, "axes.grid.which": "major", "grid.color": GRID, "grid.linewidth": 0.5,
        "axes.axisbelow": True, "lines.linewidth": 1.2, "lines.markersize": 4.5,
        "patch.linewidth": 0.5, "legend.frameon": False, "legend.handlelength": 1.8,
        "axes.prop_cycle": plt.cycler(color=[INK]), "axes.titlelocation": "left",
        "axes.unicode_minus": True, "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
        "pdf.fonttype": 42, "figure.dpi": 110, "figure.constrained_layout.use": True,
    })


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def note(message: str) -> None:
    SKIPPED.append(message)
    print(f"    - {message}")


def save(fig, name: str) -> None:
    fig.savefig(FIGURES / f"{name}.pdf")
    fig.savefig(PNG / f"{name}.png", dpi=110)
    plt.close(fig)
    print(f"  wrote {name}.pdf")


def empty(ax, message: str) -> None:
    """Mark a panel that has nothing to show."""
    ax.text(0.5, 0.5, f"no data\n{message}", transform=ax.transAxes, ha="center", va="center",
            color=INK2, fontsize=8)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    note(f"empty panel: {message}")


def ramp(n: int) -> list:
    """n colours of the one-hue blue ramp, light to dark, skipping the palest end."""
    return [RAMP(x) for x in (np.linspace(0.35, 1.0, n) if n > 1 else [0.8])]


def fit_line(ax, fit, x, label: str | None = None, **kw) -> None:
    """Draw a log-log fit ``(slope, intercept)`` over the range of x."""
    x = A.finite(x)
    x = x[x > 0]
    if fit is None or not len(x):
        return
    xs = np.geomspace(x.min(), x.max(), 20)
    ax.plot(xs, 10 ** fit[1] * xs ** fit[0], label=label, **{"color": INK2, "lw": 0.8, "ls": "--", **kw})


def slope(fit) -> str:
    return "n/a" if fit is None else f"{fit[0]:.2f}"


def status_handles(outcomes) -> list:
    return [Line2D([], [], ls="", marker="o", color=STATUS_COLORS[o], label=STATUS_LABELS[o])
            for o in A.OUTCOMES if o in set(outcomes)]


def sci(x: float) -> str:
    """10^e as mathtext for an exact power of ten, else a short float."""
    e = math.log10(x)
    return f"$10^{{{int(round(e))}}}$" if abs(e - round(e)) < 1e-9 else f"{x:.3g}"


def panel(ax, letter: str, title: str) -> None:
    ax.set_title(f"({letter}) {title}", fontsize=9)


# --------------------------------------------------------------------------- #
# E1 map throughput
# --------------------------------------------------------------------------- #
E1_STYLE = {  # backend -> colour, linestyle, marker, filled, label
    "numpy": (CPU, "-", "o", True, "NumPy (CPU)"),
    "map_batch_cpu": (CPU, "--", "D", False, "map_batch, CPU"),
    "cupy_e2e": (GPU, "-", "s", True, "CuPy end-to-end"),
    "cupy_compute": (GPU, ":", "^", False, "CuPy compute only"),
    "map_batch_gpu": (GPU, "--", "v", False, "map_batch, GPU"),
    "cupy_transfer": (INK2, "-.", "x", True, "transfer only"),
}


def fig_e01_map_throughput() -> None:
    sweeps = A.e01_sweeps()
    if not sweeps:
        raise Skip("no e01 sweep records")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=WIDE)
    for backend, (color, ls, marker, filled, label) in E1_STYLE.items():
        if backend not in sweeps:
            note(f"e01: backend {backend} missing")
            continue
        d = sweeps[backend]
        ax.plot(d["n"], d["points_per_s"], ls=ls, marker=marker, color=color, label=label,
                mfc=color if filled else "white")
        ax.fill_between(d["n"], d["n"] / d["q3_s"], d["n"] / d["q1_s"], color=color, alpha=0.15, lw=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    crossovers = A.e01_crossovers()
    none_text = []
    for pair, ls in (("raw", "-"), ("map_batch", "--")):
        c = crossovers.get(pair)
        if c is None:
            continue
        if c["n"]:
            ax.axvline(c["n"], color=INK, ls=ls, lw=0.8)
            ax.text(c["n"], 0.03, f" crossover\n ({pair.replace('_', ' ')})", transform=ax.get_xaxis_transform(),
                    fontsize=6.5, color=INK, va="bottom")
        elif c["measured"]:
            none_text.append(pair.replace("_", " "))
    if none_text:
        n_max = max(c["n_max"] for c in crossovers.values())
        ax.text(0.03, 0.97, f"no GPU crossover up to $N$ = {sci(n_max)}\n({', '.join(none_text)})",
                transform=ax.transAxes, fontsize=7, va="top", color=INK)
    ax.set_xlabel("batch size $N$ (points)")
    ax.set_ylabel("throughput (points/s)")
    panel(ax, "a", "map throughput")

    split = A.e01_split(sweeps)
    if split is None:
        empty(bx, "needs cupy_e2e, cupy_compute and cupy_transfer")
    else:
        bx.stackplot(split["n"], split["transfer"], split["compute"], split["other"],
                     colors=[LIGHT, GPU, YELLOW], labels=["transfer", "compute", "other (launch, alloc.)"],
                     alpha=0.9, lw=0)
        bx.set_xscale("log")
        bx.set_ylim(0, 1)
        bx.set_xlim(split["n"].min(), split["n"].max())
        bx.set_xlabel("batch size $N$ (points)")
        bx.set_ylabel("fraction of GPU end-to-end time")
        bx.legend(loc="upper left", fontsize=7, frameon=True, facecolor="white", edgecolor="none", framealpha=0.85)
        panel(bx, "b", "where the GPU time goes")
    fig.legend(*ax.get_legend_handles_labels(), loc="outside lower center", ncol=6, fontsize=7,
               handlelength=2.2, columnspacing=1.0)
    save(fig, "e01_map_throughput")


# --------------------------------------------------------------------------- #
# E2 end-to-end GPU
# --------------------------------------------------------------------------- #
def case_handles(cases) -> list:
    return [Line2D([], [], color=CASE_COLORS.get(c, INK), marker="o", label=c) for c in cases]


BACKEND_HANDLES = [Line2D([], [], color=INK, ls="-", marker="o", label="CPU"),
                   Line2D([], [], color=INK, ls="--", marker="o", mfc="white", label="GPU")]


def fig_e02_end_to_end() -> None:
    timed = A.e02_timed()
    if not timed:
        raise Skip("no ok e02 timed records")
    curves = A.e02_step_curves(timed)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=WIDE)
    for case in A.CASES:
        for backend, ls in (("cpu", "-"), ("gpu", "--")):
            d = curves.get(case, {}).get(backend)
            if d is None:
                note(f"e02: no timed {case}/{backend}")
                continue
            color = CASE_COLORS[case]
            ax.plot(d["points"], d["median"], ls=ls, marker="o", ms=3, color=color,
                    mfc=color if backend == "cpu" else "white")
            ax.fill_between(d["points"], d["q1"], d["q3"], color=color, alpha=0.15, lw=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("unstable points after the step")
    ax.set_ylabel("growth step time (s)")
    present = [c for c in A.CASES if c in curves]
    ax.legend(handles=case_handles(present) + BACKEND_HANDLES, loc="upper left", fontsize=7)
    reps = sorted({d["repeats"] for c in curves.values() for d in c.values()})
    panel(ax, "a", f"per-step growth (median of {'/'.join(map(str, reps))}, IQR band)")

    speedups = A.e02_stage_speedups(timed)
    if not speedups:
        empty(bx, "needs CPU and GPU timed runs of one case")
    else:
        stages = [*A.E2_STAGES, "total"]
        width = 0.8 / len(speedups)
        for i, (case, row) in enumerate(speedups.items()):
            x = np.arange(len(stages)) + (i - (len(speedups) - 1) / 2) * width
            y = [row.get(s) or np.nan for s in stages]
            bx.bar(x, y, width=width * 0.92, color=CASE_COLORS.get(case, INK), label=case)
        bx.axhline(1, color=INK, lw=0.7)
        bx.set_xticks(range(len(stages)), stages, rotation=25, ha="right")
        bx.set_ylabel("speed-up, CPU time / GPU time")
        bx.legend(loc="upper left", fontsize=7, ncol=len(speedups))
        bx.grid(axis="x", visible=False)
        top = np.nanmax([v or np.nan for r in speedups.values() for v in r.values()] + [1.0])
        bx.set_ylim(0, top * 1.25)
        panel(bx, "b", "speed-up per stage (median times)")
    save(fig, "e02_end_to_end")


PROFILE_COLORS = {"map": GPU, "refinement": AQUA, "linked_list": VIOLET, "intersections": MAGENTA, "other": LIGHT}


def fig_e02_profile() -> None:
    profiles = A.e02_profiles()
    if not profiles:
        raise Skip("no ok e02 profile records")
    rows = [(c, b) for c in A.CASES for b in ("cpu", "gpu") if (c, b) in profiles]
    fig, ax = plt.subplots(figsize=(6.3, 0.32 * len(rows) + 1.0))
    for i, key in enumerate(rows):
        left = 0.0
        for cat in A.PROFILE_CATEGORIES:
            w = profiles[key]["fractions"][cat]
            ax.barh(i, w, left=left, color=PROFILE_COLORS[cat], height=0.7,
                    label=cat.replace("_", " ") if i == 0 else None)
            left += w
        p = profiles[key]
        ax.text(1.01, i, f"$f_{{\\rm map}}$ = {100 * p['f_map']:.1f}%,  Amdahl bound {p['amdahl']:.3f}$\\times$"
                f"  ({p['total_s']:.2g} s)", va="center", fontsize=7, color=INK, transform=ax.get_yaxis_transform())
    ax.set_yticks(range(len(rows)), [f"{c} {b.upper()}" for c, b in rows])
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of profiled unstable-growth time (cProfile tottime)")
    ax.grid(axis="y", visible=False)
    fig.legend(*ax.get_legend_handles_labels(), loc="outside upper center", ncol=5, fontsize=7)
    save(fig, "e02_profile")


def fig_e02_min_batch() -> None:
    records = A.ok("e02", kind="mbp")
    if not records:
        raise Skip("no ok e02 min_batch_points records")
    records.sort(key=lambda r: r["config"]["min_batch_points"])
    labels = [("never" if r["config"]["min_batch_points"] >= 1e9 else f"{r['config']['min_batch_points']}")
              for r in records]
    fig, ax = plt.subplots(figsize=(4.2, 2.5))
    x = np.arange(len(records))
    ax.plot(x, [A.num(r["stages"].get("grow_unstable")) for r in records], marker="o", color=GPU,
            label="GPU session")
    cpu = [A.num(r["stages"].get("grow_unstable")) for r in A.e02_timed().get("k28", {}).get("cpu", [])]
    if cpu:
        ax.axhline(np.median(cpu), color=CPU, ls="--", lw=1, label="CPU session (timed runs, median)")
    ax.set_xticks(x, labels)
    ax.set_xlabel("min_batch_points (smallest batch sent to the GPU)")
    ax.set_ylabel("unstable growth time (s)")
    ax.set_ylim(bottom=0)
    steps = records[0]["config"].get("steps")
    ax.legend(loc="lower right", fontsize=7)
    ax.set_title(f"k28, {steps} unstable steps", fontsize=9)
    save(fig, "e02_min_batch")


def fig_e02_deviation() -> None:
    records = A.ok("e02", kind="verify")
    if not records:
        raise Skip("no ok e02 verify records")
    floor = 1e-17
    fig, (ax, bx) = plt.subplots(1, 2, figsize=WIDE, sharex=True)
    mismatch = False
    for r in records:
        case, s = r["config"]["case"], r["series"]
        dev = A.arr(s.get("max_abs_deviation", []))
        coord = A.arr(s.get("max_abs_coord", []))[: len(dev)]
        steps = np.arange(1, len(dev) + 1)
        color = CASE_COLORS.get(case, INK)
        for axis, y in ((ax, dev), (bx, dev / coord)):
            axis.plot(steps, np.where(np.isfinite(y), np.maximum(y, floor), np.nan), marker="o", ms=3,
                      color=color, label=case)
            bad = ~np.isfinite(dev)
            if bad.any():
                mismatch = True
                axis.plot(steps[bad], np.full(bad.sum(), 1.0), ls="", marker="x", color=color)
    finite_devs = [d for r in records for d in A.finite(r["series"].get("max_abs_deviation", []))]
    top = 10.0 ** math.ceil(math.log10(max(max(finite_devs, default=0), 1e-13)) + 1)
    if mismatch:
        top = 3.0
    for axis, letter, title, ylabel in ((ax, "a", "absolute", "max |CPU − GPU| coordinate"),
                                        (bx, "b", "relative to max |coord|", "max |CPU − GPU| / max |coord|")):
        axis.set_yscale("log")
        axis.set_ylim(floor / 3, top)
        ticks = [floor] + [10.0 ** e for e in range(-15, int(round(math.log10(top))) + 1, 3)]
        axis.set_yticks(ticks, ["0"] + [sci(t) for t in ticks[1:]])
        axis.minorticks_off()
        axis.axhline(np.finfo(float).eps, color=INK2, lw=0.6, ls=":")
        axis.text(1, np.finfo(float).eps * 1.6, r"machine $\epsilon$", fontsize=6.5, color=INK2)
        axis.set_xlabel("unstable growth step")
        axis.set_ylabel(ylabel)
        panel(axis, letter, title)
    ax.legend(loc="upper left", fontsize=7)
    if mismatch:
        ax.text(0.98, 0.97, "× = point counts differ", transform=ax.transAxes, ha="right", va="top", fontsize=7)
    save(fig, "e02_deviation")


# --------------------------------------------------------------------------- #
# E3 scaling
# --------------------------------------------------------------------------- #
E3_K_COLORS = {2.8: BLUE, 4.0: YELLOW, 6.0: GREEN, 10.0: VIOLET}


def e3_style(k: float, cutoff: float) -> dict:
    label = f"k = {k:g}" + ("" if cutoff == 1e-7 else f", cutoff {sci(cutoff)}")
    return {"color": E3_K_COLORS.get(k, INK), "ls": "-" if cutoff == 1e-7 else "--", "label": label}


def fig_e03_growth() -> None:
    growth = A.e03_growth()
    if not growth:
        raise Skip("no e03 growth records")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=WIDE, sharex=True)
    for (k, cutoff), g in sorted(growth.items()):
        st = e3_style(k, cutoff)
        ax.plot(g["step"], g["points_unstable"], marker="o", ms=3, **st)
        bx.plot(g["step"], g["max_abs_coord"], marker="o", ms=3, **st)
        anchor = int(np.argmax(g["points_unstable"] >= 100)) if (g["points_unstable"] >= 100).any() else None
        if anchor is not None and np.isfinite(g["lambda_u"]):
            n = g["step"][anchor:]
            ax.plot(n, g["points_unstable"][anchor] * abs(g["lambda_u"]) ** (n - n[0]),
                    color=st["color"], ls=":", lw=0.9)
        if g["outcome"] != "ok":
            ax.plot(g["step"][-1], g["points_unstable"][-1], marker="X", ms=7, color=STATUS_COLORS.get(g["outcome"], INK))
    ax.set_yscale("log")
    bx.set_yscale("log")
    handles, _ = ax.get_legend_handles_labels()
    handles.append(Line2D([], [], color=INK2, ls=":", label=r"$\propto |\lambda_u|^{n}$"))
    ax.legend(handles=handles, loc="upper left", fontsize=7)
    ax.set_xlabel("unstable step $n$")
    ax.set_ylabel("unstable points")
    bx.set_xlabel("unstable step $n$")
    bx.set_ylabel("max |coordinate|")
    panel(ax, "a", "point count")
    panel(bx, "b", "escape to infinity")
    save(fig, "e03_growth")


def fig_e03_cost() -> None:
    growth, rows = A.e03_growth(), [r for r in A.e03_stages() if r["outcome"] == "ok"]
    if not growth and not rows:
        raise Skip("no e03 records")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=WIDE)
    if growth:
        for (k, cutoff), g in sorted(growth.items()):
            st = e3_style(k, cutoff)
            st["label"] += f" (slope {slope(g['step_fit'])})"
            ax.plot(g["points_unstable"], g["step_s"], marker="o", ms=3, **st)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.legend(loc="upper left", fontsize=6.5)
        ax.set_xlabel("unstable points after the step")
        ax.set_ylabel("step time (s)")
        panel(ax, "a", r"growth step time (fit over $\geq 10^3$ points)")
    else:
        empty(ax, "e03 growth")
    if rows:
        fits = A.e03_stage_fits(rows)
        ks = sorted({r["k"] for r in rows})
        for stage in A.E3_STAGES:
            for r in rows:
                bx.plot(r["segments"], r[stage], ls="", marker=MARKERS[ks.index(r["k"]) % len(MARKERS)],
                        color=STAGE_COLORS[stage], mfc="white" if r["cutoff"] != 1e-7 else STAGE_COLORS[stage])
            fit_line(bx, fits[stage], [r["segments"] for r in rows], color=STAGE_COLORS[stage],
                     label=f"{STAGE_LABELS[stage]} (slope {slope(fits[stage])})")
        bx.set_xscale("log")
        bx.set_yscale("log")
        handles, _ = bx.get_legend_handles_labels()
        handles += [Line2D([], [], ls="", marker=MARKERS[i % len(MARKERS)], color=INK2, label=f"k = {k:g}")
                    for i, k in enumerate(ks)]
        fig.legend(handles=handles, loc="outside lower right", ncol=3, fontsize=6.5,
                   title="(b) open markers: cutoff $10^{-5}$", title_fontsize=6.5)
        bx.set_xlabel("segments (unstable + stable)")
        bx.set_ylabel("stage time (s)")
        panel(bx, "b", "downstream stages")
    else:
        empty(bx, "e03 stages")
    save(fig, "e03_cost")


def fig_e03_memory() -> None:
    growth = A.e03_growth()
    if not growth:
        raise Skip("no e03 growth records")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=WIDE, sharex=True)
    for (k, cutoff), g in sorted(growth.items()):
        st = e3_style(k, cutoff)
        fit = g["bytes_per_point_fit"]
        ax.plot(g["points_unstable"], np.maximum(g["rss_mb"] - g["baseline_rss_mb"], 1e-2), marker="o", ms=3,
                **{**st, "label": st["label"] + (f" ({fit:.0f} B/point)" if fit else "")})
        bx.plot(g["points_unstable"], g["bytes_per_point"], marker="o", ms=3, **st)
    rows = [r for r in A.e03_stages() if np.isfinite(r["bytes_per_point"])]
    for r in rows:
        bx.plot(r["points"], r["bytes_per_point"], ls="", marker="s", mfc="white",
                color=E3_K_COLORS.get(r["k"], INK))
    for axis in (ax, bx):
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel("unstable points")
    ax.set_ylabel("RSS above baseline (MB)")
    bx.set_ylabel("bytes per point")
    ax.legend(loc="upper left", fontsize=6.5)
    if rows:
        bx.legend(handles=[Line2D([], [], ls="", marker="s", mfc="white", color=INK2,
                                  label="stage runs (all points)")], loc="upper right", fontsize=6.5)
    panel(ax, "a", "resident memory (slope fit, top two decades)")
    panel(bx, "b", "bytes per point")
    save(fig, "e03_memory")


# --------------------------------------------------------------------------- #
# E4 area cutoff
# --------------------------------------------------------------------------- #
def e4_style(case: str, steps: int, runs: dict) -> dict:
    depths = sorted(s for c, s in runs if c == case)
    deep = steps != depths[0]
    return {"color": CASE_COLORS.get(case, INK), "marker": "s" if deep else "o",
            "label": f"{case}, {steps} steps"}


def fig_e04_pareto() -> None:
    runs = A.e04_runs()
    if not runs:
        raise Skip("no e04 records")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=WIDE)
    for (case, steps), rows in runs.items():
        st = e4_style(case, steps, runs)
        good = [r for r in rows if r["outcome"] == "ok" and r["inv_p50"] is not None]
        for r in rows:
            if r["outcome"] != "ok":
                note(f"e04: {case}/{steps} cutoff {r['cutoff']:g} is {r['outcome']}")
        if good:
            t = [r["time"] for r in good]
            ax.plot(t, [r["inv_p50"] for r in good], color=st["color"], marker=st["marker"], label=st["label"])
            ax.plot(t, [r["inv_p99"] for r in good], color=st["color"], marker=st["marker"], mfc="white", ls=":", lw=0.9)
            for r in good:
                ax.annotate(f"{math.log10(r['cutoff']):.0f}", (r["time"], r["inv_p50"]), xytext=(3, 3),
                            textcoords="offset points", fontsize=6, color=INK2)
        done = [r for r in rows if np.isfinite(r["crossings"])]
        bx.plot([r["cutoff"] for r in done], [r["crossings"] for r in done], color=st["color"],
                marker=st["marker"], label=st["label"])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("pipeline time (s, without the measurement stages)")
    ax.set_ylabel("invariance error (max over branches)")
    handles, _ = ax.get_legend_handles_labels()
    handles += [Line2D([], [], color=INK2, marker="o", label="p50"),
                Line2D([], [], color=INK2, marker="o", mfc="white", ls=":", label="p99")]
    ax.legend(handles=handles, loc="lower left", fontsize=6.5)
    panel(ax, "a", r"accuracy vs cost (labels: $\log_{10}$ cutoff)")
    bx.set_xscale("log")
    bx.invert_xaxis()
    bx.set_xlabel("area_cutoff (finer to the right)")
    bx.set_ylabel("registered crossings (no anchors)")
    bx.legend(loc="best", fontsize=6.5)
    panel(bx, "b", "crossings found")
    save(fig, "e04_pareto")


def fig_e04_convergence() -> None:
    runs = A.e04_runs()
    if not runs:
        raise Skip("no e04 records")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=WIDE)
    drew = False
    for (case, steps), rows in runs.items():
        st = e4_style(case, steps, runs)
        conv = [c for c in A.e04_convergence(rows) if c["median"] is not None]
        if conv:
            drew = True
            x = [c["cutoff"] for c in conv]
            fit = A.fit_loglog(x, [c["median"] for c in conv])
            ax.plot(x, [c["median"] for c in conv], color=st["color"], marker=st["marker"],
                    label=f"{st['label']} (slope {slope(fit)})")
            ax.plot(x, [c["max"] for c in conv], color=st["color"], marker=st["marker"], mfc="white", ls=":", lw=0.9)
        else:
            note(f"e04: {case}/{steps} has fewer than two ok runs with crossings, no convergence curve")
        done = [r for r in rows if np.isfinite(r["points"])]
        fit = A.fit_loglog([r["cutoff"] for r in done], [r["points"] for r in done])
        bx.plot([r["cutoff"] for r in done], [r["points"] for r in done], color=st["color"],
                marker=st["marker"], label=f"{st['label']} (slope {slope(fit)})")
    if drew:
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.invert_xaxis()
        handles, _ = ax.get_legend_handles_labels()
        handles += [Line2D([], [], color=INK2, marker="o", label="median"),
                    Line2D([], [], color=INK2, marker="o", mfc="white", ls=":", label="max")]
        ax.legend(handles=handles, loc="lower left", fontsize=6.5)
        ax.set_xlabel("area_cutoff (finer to the right)")
        ax.set_ylabel("distance to the finest run's crossing")
        panel(ax, "a", "crossing convergence")
    else:
        empty(ax, "needs two or more cutoffs per case")
    bx.set_xscale("log")
    bx.set_yscale("log")
    bx.invert_xaxis()
    bx.set_xlabel("area_cutoff (finer to the right)")
    bx.set_ylabel("points (unstable + stable)")
    bx.legend(loc="upper left", fontsize=6.5)
    panel(bx, "b", "cost in points")
    save(fig, "e04_convergence")


# --------------------------------------------------------------------------- #
# E5 invariants
# --------------------------------------------------------------------------- #
E5_STYLE = {  # case -> colour, linestyle, label
    "k10": (VIOLET, "-", "k10"), "k28_unblasted": (BLUE, "-", "k28"),
    "k28_one_blast": (BLUE, "--", "k28, 1 blast"), "k28_two_blasts": (BLUE, ":", "k28, 2 blasts"),
    "p3": (AQUA, "-", "p3"), "nested": (MAGENTA, "-", "nested"), "inversion": (YELLOW, "-", "inversion"),
}


def e5_style(case: str):
    return E5_STYLE.get(case, (INK, "-", case))


def fig_e05_invariance() -> None:
    records = [r for r in A.e05_cases() if r["series"]]
    if not records:
        raise Skip("no e05 case records with series")
    floor = 1e-17
    fig, axes = plt.subplots(1, 2, figsize=WIDE, sharey=True)
    for stability, ax, letter in (("unstable", axes[0], "a"), ("stable", axes[1], "b")):
        handles = []
        for r in records:
            color, ls, label = e5_style(r["config"]["case"])
            drawn = False
            for tag, err in A.e05_invariance_samples(r).items():
                if f":{stability}/" not in tag or not len(err):
                    continue
                e = np.sort(np.maximum(err, floor))
                ax.step(e, np.arange(1, len(e) + 1) / len(e), where="post", color=color, ls=ls, lw=1)
                drawn = True
            if drawn:
                handles.append(Line2D([], [], color=color, ls=ls, label=label))
        ax.set_xscale("log")
        ax.set_xlabel("invariance error (distance of the image to the curve)")
        panel(ax, letter, f"{stability} manifolds (one line per branch)")
        if handles:
            ax.legend(handles=handles, loc="upper left", fontsize=6.5)
        else:
            empty(ax, f"no {stability} samples")
    axes[0].set_ylabel("cumulative fraction of samples")
    save(fig, "e05_invariance")


def fig_e05_iterates_lobes() -> None:
    records = A.e05_cases()
    if not records:
        raise Skip("no e05 case records")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=WIDE)
    rng = np.random.default_rng(0)
    handles = []
    for r in records:
        rows = r["series"].get("iter_links") or []
        if not rows:
            continue
        color, ls, label = e5_style(r["config"]["case"])
        depth = A.arr(row[1] for row in rows)
        err = A.arr(row[3] for row in rows)
        ax.plot(depth + rng.uniform(-0.15, 0.15, len(depth)), np.maximum(err, 1e-17), ls="", marker="o",
                ms=3, color=color, alpha=0.7)
        handles.append(Line2D([], [], ls="", marker="o", color=color, label=label))
    if handles:
        ax.set_yscale("log")
        ax.set_xlabel("backward depth of the link's source in the iterate table")
        ax.set_ylabel(r"$|f^n(x_{\rm source}) - x_{\rm target}|$")
        ax.legend(handles=handles, loc="lower right", fontsize=6.5)
        panel(ax, "a", "iterate-link error")
    else:
        empty(ax, "no iterate links")

    cases = [r["config"]["case"] for r in records]
    any_ratio = False
    for i, r in enumerate(records):
        color, _, _ = e5_style(r["config"]["case"])
        ratio = A.e05_lobe_ratios(r)
        if len(ratio):
            any_ratio = True
            bx.plot(i + rng.uniform(-0.12, 0.12, len(ratio)) - 0.15, np.maximum(ratio, 1e-17), ls="",
                    marker="o", ms=3.5, color=color)
        arr_rows = [row for row in r["series"].get("arr_regions") or [] if row[1] is not None]
        if arr_rows:
            any_ratio = True
            ar = A.arr(abs(A.num(row[1]) / A.num(row[0]) - 1) for row in arr_rows)
            bx.plot(i + rng.uniform(-0.12, 0.12, len(ar)) + 0.15, np.maximum(ar, 1e-17), ls="", marker="s",
                    ms=3.5, color=color, mfc="white")
    if any_ratio:
        bx.set_yscale("log")
        bx.set_xticks(range(len(records)), [e5_style(c)[2] for c in cases], rotation=25, ha="right")
        bx.set_ylabel(r"$|A_{\rm image} / A - 1|$")
        bx.grid(axis="x", visible=False)
        bx.legend(handles=[Line2D([], [], ls="", marker="o", color=INK2, label="bridge lobes"),
                           Line2D([], [], ls="", marker="s", mfc="white", color=INK2, label="arrangement regions")],
                  loc="upper right", fontsize=6.5)
        panel(bx, "b", "area preservation")
    else:
        empty(bx, "no lobe or region pairs")
    save(fig, "e05_iterates_lobes")


# --------------------------------------------------------------------------- #
# E6 parameter sweep
# --------------------------------------------------------------------------- #
E6_FAMILIES = ("b = 1, saddle", "b = 1, other point", "b = -1", "b = 0.9", "b = 0.5", "b = 0.3")
E6_STAGE_LABELS = ("none", "fixed point", "eigenvectors", "seeded", "unstable grown", "both grown",
                   "intersected", "bridges", "trellis + pips")


def fig_e06_sweep_outcomes() -> None:
    rows = A.e06_rows()
    if not rows:
        raise Skip("no e06 records")
    families = [*E6_FAMILIES, *sorted({r["family"] for r in rows} - set(E6_FAMILIES))]
    ncols = 3
    nrows = math.ceil(len(families) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.3, 1.9 * nrows + 0.5), sharex=True, sharey=True)
    axes = np.atleast_2d(axes)
    for i, family in enumerate(families):
        ax = axes.flat[i]
        mine = sorted((r for r in rows if r["family"] == family), key=lambda r: r["k"])
        n_ok = sum(r["outcome"] == "ok" and r["stage"] == "trellis_pips" for r in mine)
        ax.set_title(family.replace("-", "−") + (f"   ({n_ok}/{len(mine)} complete)" if mine else ""), fontsize=8)
        ax.set_ylim(-0.6, len(A.E6_STAGES) - 0.4)
        if not mine:
            ax.text(0.5, 0.5, "no configs", transform=ax.transAxes, ha="center", va="center", color=INK2, fontsize=7.5)
            continue
        ax.plot([r["k"] for r in mine], [r["stage_index"] for r in mine], color=LIGHT, lw=0.7, zorder=1)
        for r in mine:
            ax.plot(r["k"], r["stage_index"], ls="", marker="o", ms=4.5, zorder=2,
                    color=STATUS_COLORS.get(r["outcome"], INK))
    for ax in axes.flat[len(families):]:
        ax.set_visible(False)
    for ax in axes.flat:
        ax.set_xscale("log")
        ax.set_yticks(range(len(A.E6_STAGES)), E6_STAGE_LABELS, fontsize=6.5)
        ax.grid(axis="x", visible=False)
    for ax in axes[-1]:
        ax.set_xlabel("$k$")
    fig.legend(handles=status_handles([r["outcome"] for r in rows]), loc="outside upper center",
               ncol=6, fontsize=7, title="outcome of the config (height = last stage completed)", title_fontsize=7)
    save(fig, "e06_sweep_outcomes")


def fig_e06_saddle_sweep() -> None:
    rows = sorted((r for r in A.e06_rows() if r["family"] == "b = 1, saddle"), key=lambda r: r["k"])
    if not rows:
        raise Skip("no b = 1 saddle configs in e06")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=WIDE, sharex=True)
    k = [r["k"] for r in rows]
    for key, color, marker, label in (("n_crossings", GREEN, "o", "crossings"),
                                      ("n_bridges", RED, "s", "bridges (with id)")):
        ax.plot(k, [A.num(r["metrics"].get(key)) for r in rows], color=color, marker=marker, ms=3.5,
                mfc="white" if marker == "s" else color, label=label)
    failed = [r for r in rows if r["outcome"] != "ok"]
    for r in failed:
        ax.plot(r["k"], 0, ls="", marker="X", ms=6, color=STATUS_COLORS.get(r["outcome"], INK), clip_on=False)
    ax.set_yscale("symlog", linthresh=1)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("$k$")
    ax.set_ylabel("count")
    handles, _ = ax.get_legend_handles_labels()
    ax.legend(handles=handles + [Line2D([], [], ls="", marker="X", color=STATUS_COLORS[o],
                                        label=f"failed: {STATUS_LABELS[o]}")
                                 for o in A.OUTCOMES if o in {r['outcome'] for r in failed}],
              loc="upper left", fontsize=6.5)
    panel(ax, "a", "what the recipe finds (b = 1 saddle)")
    for stage, key in (("grow_unstable", "grow_unstable"), ("grow_stable", "grow_stable"),
                       ("compute_intersections", "intersections")):
        bx.plot(k, [A.num(r["stages"].get(key)) for r in rows], marker="o", ms=3, color=STAGE_COLORS[stage],
                label=STAGE_LABELS[stage])
    bx.plot(k, [sum(A.num(v) for v in r["stages"].values()) or np.nan for r in rows], color=INK, lw=0.9,
            label="all stages")
    bx.set_xscale("log")
    bx.set_yscale("log")
    bx.set_xlabel("$k$")
    bx.set_ylabel("time (s)")
    bx.legend(loc="best", fontsize=6.5)
    panel(bx, "b", "cost")
    save(fig, "e06_saddle_sweep")


# --------------------------------------------------------------------------- #
# E7 breaking
# --------------------------------------------------------------------------- #
def fig_e07_escape() -> None:
    records = A.e07("escape")
    fixed = next((r for r in records if r["config"].get("mode") == "fixed" and r["series"].get("step")), None)
    if fixed is None:
        raise Skip("no e07 fixed-count escape run with steps")
    s = fixed["series"]
    step, points, max_abs = A.arr(s["step"]), A.arr(s["points"]), A.arr(s["max_abs"])
    fig, (ax, bx) = plt.subplots(1, 2, figsize=WIDE, sharex=True)
    ax.plot(step, points, color=INK, marker="o", ms=3, label="fixed step count")
    bx.plot(step, max_abs, color=INK, marker="o", ms=3, label="fixed step count")
    if fixed["outcome"] != "ok":
        for axis, y in ((ax, points), (bx, max_abs)):
            axis.plot(step[-1], y[-1], ls="", marker="X", ms=8, color=STATUS_COLORS.get(fixed["outcome"], INK),
                      label=f"run ended: {STATUS_LABELS.get(fixed['outcome'], fixed['outcome'])}")
    # Stop runs that end on the same step sit side by side, not on top of each other.
    stops = {"arclength": (VIOLET, "^", -0.15), "predicate": (AQUA, "D", 0.15)}
    for r in records:
        mode = r["config"].get("mode")
        if mode not in stops or r["outcome"] != "ok" or "final_points" not in r["metrics"]:
            if mode in stops:
                note(f"e07: escape {mode} run is {r['outcome']} without final metrics")
            continue
        m = r["metrics"]
        at = m.get("rounds")
        if at is None:
            hit = np.nonzero(points >= m["final_points"])[0]
            at = step[hit[0]] if len(hit) else step[-1]
        what = (f"arclength {r['config']['length']:g}" if mode == "arclength"
                else f"stop at max|x| > {r['config']['radius']:g}")
        color, marker, dx = stops[mode]
        ax.plot(at + dx, m["final_points"], ls="", marker=marker, ms=6, color=color, mfc="white", mew=1.2, label=what)
        bx.plot(at + dx, m["final_max_abs"], ls="", marker=marker, ms=6, color=color, mfc="white", mew=1.2, label=what)
    blowup = A.e07_blowup_step(fixed)
    for axis in (ax, bx):
        axis.set_yscale("log")
        axis.set_xlabel("unstable step")
        if blowup:
            axis.axvline(blowup, color=INK2, ls=":", lw=0.8)
    ax.set_ylabel("unstable points")
    bx.set_ylabel("max |coordinate|")
    ax.legend(loc="upper left", fontsize=6.5)
    panel(ax, "a", "period 3, k = 2.1: point count")
    panel(bx, "b", "escape" + (f" (max|x| > 100 at step {blowup})" if blowup else ""))
    save(fig, "e07_escape")


BLAST_CASES = ("k28", "p3", "nested")


def fig_e07_blasts() -> None:
    records = A.e07("blasts")
    if not any(r["series"].get("crossings") for r in records):
        raise Skip("no e07 blast series")
    fig, axes = plt.subplots(1, 3, figsize=(6.3, 2.4), sharex=True)
    for r in records:
        s, case = r["series"], r["config"].get("case")
        if not s.get("crossings"):
            note(f"e07: blasts {case} has no series ({r['outcome']})")
            continue
        color = CASE_COLORS.get(case, INK)
        idx = np.arange(len(s["crossings"]))
        m = r["metrics"]
        topo = m.get("topology_outcome")
        if topo is None:
            tag = "topology not run"
        elif topo != "ok":
            tag = f"topology {topo.replace('_', ' ')}"
        else:
            tag = f"{m.get('n_classes')} classes" + ("" if m.get("is_reliable", True) else
                                                      f", {m.get('n_unresolved')} unresolved, {m.get('n_virtual')} virtual")
        axes[0].plot(idx[1:], A.arr(s.get("blast_s", []))[1:], marker="o", ms=3, color=color, label=f"{case} ({tag})")
        axes[1].plot(idx, A.arr(s["crossings"]), marker="o", ms=3, color=color)
        axes[2].plot(idx, A.arr(s.get("bridges", [])), marker="o", ms=3, color=color)
        if r["outcome"] != "ok":
            for axis, key in ((axes[1], "crossings"), (axes[2], "bridges")):
                axis.plot(idx[-1], A.num(s[key][-1]), ls="", marker="X", ms=7, color=STATUS_COLORS.get(r["outcome"], INK))
    axes[0].set_yscale("log")
    axes[1].set_yscale("log")
    axes[2].set_yscale("log")
    axes[0].set_ylabel("blast time (s)")
    axes[1].set_ylabel("registered crossings")
    axes[2].set_ylabel("bridges")
    for axis, letter, title in zip(axes, "abc", ("time per blast", "crossings", "bridges")):
        axis.set_xlabel("blasts done")
        axis.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        panel(axis, letter, title)
    handles, _ = axes[0].get_legend_handles_labels()
    failed = {r["outcome"] for r in records if r["outcome"] != "ok"}
    handles += [Line2D([], [], ls="", marker="X", color=STATUS_COLORS[o], label=f"stopped: {STATUS_LABELS[o]}")
                for o in A.OUTCOMES if o in failed]
    fig.legend(handles=handles, loc="outside lower center", ncol=min(len(handles), 2), fontsize=7)
    save(fig, "e07_blasts")


def fig_e07_depth() -> None:
    records = A.e07("depth")
    if not records:
        raise Skip("no e07 depth records")
    cases = [c for c in ("p3", "inversion") if any(r["config"].get("case") == c for r in records)]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=WIDE)
    def pipeline_s(r):  # growth steps plus the timed stages; a killed run keeps its checkpoint
        return sum(A.finite(r["series"].get("step_s", []))) + sum(A.num(v) for v in r["stages"].values()) or np.nan

    for case in cases:
        mine = sorted((r for r in records if r["config"]["case"] == case), key=lambda r: r["config"]["unstable_steps"])
        steps = [r["config"]["unstable_steps"] for r in mine]
        ax.plot(steps, [pipeline_s(r) for r in mine], color=CASE_COLORS[case], marker="o", ms=3.5, label=case)
        for r in mine:
            if r["outcome"] != "ok":
                ax.plot(r["config"]["unstable_steps"], pipeline_s(r), ls="", marker="X", ms=7,
                        color=STATUS_COLORS.get(r["outcome"], INK))
    ax.set_yscale("log")
    ax.set_xlabel("unstable steps")
    ax.set_ylabel("pipeline time (s)")
    ax.legend(loc="upper left", fontsize=7)
    panel(ax, "a", "cost of depth")

    lanes, seen = [], set()
    for case in cases:
        lanes += [(case, "numerics"), (case, "topology")]
    for y, (case, layer) in enumerate(lanes):
        for r in (r for r in records if r["config"]["case"] == case):
            m = r["metrics"]
            if layer == "numerics":
                status = r["outcome"] if r["outcome"] != "ok" or m.get("numerics_done") else "exception"
                reliable = True
            else:
                status, reliable = m.get("topology_outcome"), m.get("is_reliable", True)
                if status is None:
                    continue
            seen.add(status)
            color = STATUS_COLORS.get(status, INK)
            bx.plot(r["config"]["unstable_steps"], y, ls="", marker="o", ms=6, color=color,
                    mfc=color if reliable else "white", mew=1.2)
    bx.set_yticks(range(len(lanes)), [f"{c} {layer}" for c, layer in lanes])
    bx.set_ylim(-0.6, len(lanes) - 0.4)
    bx.invert_yaxis()
    bx.set_xlabel("unstable steps")
    bx.grid(axis="y", visible=False)
    fig.legend(handles=status_handles(seen) + [Line2D([], [], ls="", marker="o", mfc="white", color=INK2,
                                                      label="symbolic result provisional (is_reliable=False)")],
               loc="outside lower center", ncol=4, fontsize=6.5)
    panel(bx, "b", "outcome per depth (numerics, topology)")
    save(fig, "e07_depth")


def fig_e07_tangency() -> None:
    records = [r for r in A.e07("tangency") if "non_anchor_crossings" in r["metrics"]]
    if not records:
        raise Skip("no e07 tangency records with crossings")
    cutoffs = sorted({r["config"]["area_cutoff"] for r in records}, reverse=True)
    colors = dict(zip(cutoffs, ramp(len(cutoffs))))
    fig, (ax, bx) = plt.subplots(1, 2, figsize=WIDE, sharex=True)
    for cutoff in cutoffs:
        mine = sorted((r for r in records if r["config"]["area_cutoff"] == cutoff), key=lambda r: r["config"]["k"])
        x = [r["config"]["k"] + 1 for r in mine]
        ax.plot(x, [A.num(r["metrics"]["non_anchor_crossings"]) for r in mine], marker="o", ms=3,
                color=colors[cutoff], label=f"cutoff {sci(cutoff)}")
        bx.plot(x, [A.num(r["metrics"].get("min_sin")) for r in mine], marker="o", ms=3, color=colors[cutoff],
                label=f"cutoff {sci(cutoff)}")
    ref = sorted({(r["config"]["k"] + 1, A.num(r["metrics"].get("splitting_reference"))) for r in records})
    bx.plot(*zip(*ref), color=INK, ls="--", lw=0.9, label=r"$e^{-\pi^2/\ln\lambda}$ (heuristic)")
    for axis in (ax, bx):
        axis.set_xscale("log")
        axis.set_xlabel("$k + 1$ (distance to the saddle-centre bifurcation)")
    ax.set_yscale("log")
    bx.set_yscale("log")
    ax.set_ylabel("crossings (no anchors)")
    bx.set_ylabel("smallest |sin| at a crossing")
    ax.legend(loc="upper right", fontsize=6.5)
    bx.legend(loc="lower right", fontsize=6.5)
    panel(ax, "a", "crossings near tangency")
    panel(bx, "b", "crossing angle vs splitting")
    save(fig, "e07_tangency")


# --------------------------------------------------------------------------- #
# E8 solver
# --------------------------------------------------------------------------- #
MAP_COLORS = {"k = 10": VIOLET, "k = 2.8": BLUE, "k = 2": AQUA, "k = 2.1": MAGENTA}
E8_COLORS = {"saddle": GREEN, "lower_period": YELLOW, "non_saddle": VIOLET, "no_convergence": RED, "other": INK}
E8_LABELS = {"saddle": "saddle, requested period", "lower_period": "converged to a lower period",
             "non_saddle": "elliptic / not hyperbolic", "no_convergence": "no convergence", "other": "other error"}


def e8_maps(table: dict) -> list[str]:
    return sorted({m for m, _ in table}, key=lambda m: [*MAP_COLORS].index(m) if m in MAP_COLORS else 99)


def fig_e08_solver_outcomes() -> None:
    table = A.e08_table()
    if not table:
        raise Skip("no e08 records")
    maps = e8_maps(table)
    fig, ax = plt.subplots(figsize=WIDE)
    ticks, labels = [], []
    for i, map_name in enumerate(maps):
        periods = sorted(p for m, p in table if m == map_name)
        for p in periods:
            row, x = table[(map_name, p)], i * (len(periods) + 1.5) + p
            bottom = 0.0
            for cls in A.E8_CLASSES:
                frac = row["counts"][cls] / row["n"] if row["n"] else 0
                ax.bar(x, frac, bottom=bottom, color=E8_COLORS[cls], width=0.82)
                bottom += frac
            if row["outcome"] != "ok":
                ax.plot(x, 1.04, marker="X", ls="", color=STATUS_COLORS.get(row["outcome"], INK), clip_on=False)
            ticks.append(x)
            labels.append(str(p))
        mid = i * (len(periods) + 1.5) + (min(periods) + max(periods)) / 2
        ax.text(mid, -0.15, map_name, ha="center", va="top", transform=ax.get_xaxis_transform(), fontsize=8)
    ax.set_xticks(ticks, labels, fontsize=7)
    ax.set_ylim(0, 1)
    ax.set_ylabel("fraction of solves")
    ax.grid(axis="x", visible=False)
    ax.set_xlabel("requested period (grouped by map)", labelpad=20)
    n = sorted({row["n"] for row in table.values()})
    fig.legend(handles=[Patch(color=E8_COLORS[c], label=E8_LABELS[c]) for c in A.E8_CLASSES],
              loc="outside upper center", ncol=3, fontsize=7,
              title=f"{'/'.join(map(str, n))} solves per bar (orbit + independent guesses)", title_fontsize=7)
    save(fig, "e08_solver_outcomes")


def fig_e08_orbits() -> None:
    table = A.e08_table()
    if not table:
        raise Skip("no e08 records")
    maps = e8_maps(table)
    distinct = A.e08_distinct_by_period(table)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=WIDE, sharex=True)
    periods = sorted(A.HORSESHOE_ORBITS)
    ax.plot(periods, [A.HORSESHOE_ORBITS[p] for p in periods], color=INK, ls="--", marker="x", lw=0.9,
            label="full 2-shift (k = 10 expected)")
    for map_name in maps:
        found = distinct.get(map_name, {})
        ax.plot(periods, [found.get(p, 0) for p in periods], marker="o", ms=3.5, color=MAP_COLORS.get(map_name, INK),
                label=map_name)
        rows = [(p, table[(map_name, p)]) for p in periods if (map_name, p) in table]
        bx.plot([p for p, _ in rows], [row["lower_period"] / row["converged"] if row["converged"] else np.nan
                                       for _, row in rows],
                marker="o", ms=3.5, color=MAP_COLORS.get(map_name, INK), label=map_name)
    ax.set_xlabel("minimal period")
    ax.set_ylabel("distinct orbits found (all configs)")
    ax.legend(loc="upper left", fontsize=6.5)
    panel(ax, "a", "orbits found vs the horseshoe count")
    bx.set_ylim(-0.03, 1.03)
    bx.set_xlabel("requested period")
    bx.set_ylabel("converged solves of lower minimal period")
    bx.legend(loc="upper right", fontsize=6.5)
    panel(bx, "b", "wrong-period rate")
    save(fig, "e08_orbits")


def fig_e08_solve_time() -> None:
    table = A.e08_table()
    if not table:
        raise Skip("no e08 records")
    maps = e8_maps(table)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=WIDE)
    for map_name in maps:
        color = MAP_COLORS.get(map_name, INK)
        times = np.sort(np.concatenate([A.finite(row["solve_s"]) for (m, _), row in table.items() if m == map_name]))
        if len(times):
            ax.step(times, np.arange(1, len(times) + 1) / len(times), where="post", color=color, label=map_name)
        periods = sorted(p for m, p in table if m == map_name)
        med, lo, hi = [], [], []
        for p in periods:
            t = A.finite(table[(map_name, p)]["solve_s"])
            q = np.percentile(t, [50, 25, 75]) if len(t) else [np.nan] * 3
            med.append(q[0]), lo.append(q[0] - q[1]), hi.append(q[2] - q[0])
        bx.errorbar(periods, med, yerr=[lo, hi], marker="o", ms=3.5, color=color, capsize=2, lw=1, label=map_name)
    ax.set_xscale("log")
    ax.set_xlabel("solve time (s)")
    ax.set_ylabel("cumulative fraction of solves")
    ax.legend(loc="upper left", fontsize=6.5)
    panel(ax, "a", "solve-time distribution")
    bx.set_yscale("log")
    bx.set_xlabel("requested period")
    bx.set_ylabel("solve time (s), median and IQR")
    panel(bx, "b", "by period")
    save(fig, "e08_solve_time")


# --------------------------------------------------------------------------- #
# E9 intersections
# --------------------------------------------------------------------------- #
def fig_e09_intersections() -> None:
    rows = [r for r in A.e09_rows() if r["outcome"] == "ok"]
    if not rows:
        raise Skip("no ok e09 records")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=WIDE)
    cases = [c for c in CASE_COLORS if any(r["case"] == c for r in rows)]
    for r in rows:
        color = CASE_COLORS.get(r["case"], INK)
        ax.plot(r["segments"], r["library_s"], ls="", marker="o", color=color)
        ax.plot(r["segments"], r["brute_s"], ls="", marker="s", color=color, mfc="white")
        bx.plot(r["brute"], r["library"], ls="", marker="o", color=color)
    seg = [r["segments"] for r in rows]
    lib_fit = A.fit_loglog(seg, [r["library_s"] for r in rows])
    brute_fit = A.fit_loglog(seg, [r["brute_s"] for r in rows])
    fit_line(ax, lib_fit, seg, ls="-")
    fit_line(ax, brute_fit, seg, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("segments (unstable + stable)")
    ax.set_ylabel("time (s)")
    methods = [Line2D([], [], ls="-", marker="o", color=INK2, label=f"library (slope {slope(lib_fit)})"),
               Line2D([], [], ls="--", marker="s", mfc="white", color=INK2,
                      label=f"brute force (slope {slope(brute_fit)})")]
    fig.legend(handles=methods + [Line2D([], [], ls="", marker="o", color=CASE_COLORS[c], label=c) for c in cases],
               loc="outside lower center", ncol=2 + len(cases), fontsize=7)
    panel(ax, "a", "crossing detection time")
    top = max(max(r["brute"], r["library"]) for r in rows) + 1
    bx.plot([0, top], [0, top], color=INK2, ls=":", lw=0.8)
    bx.set_xlabel("brute-force crossings")
    bx.set_ylabel("library crossings (no anchors)")
    missed = sum(sum(r[k] for k in A.E9_MISSED) for r in rows)
    extra = sum(sum(r[k] for k in A.E9_EXTRA) for r in rows)
    bx.text(0.04, 0.96, f"{len(rows)} configs: {missed:.0f} missed, {extra:.0f} extra",
            transform=bx.transAxes, va="top", fontsize=7)
    panel(bx, "b", "agreement")
    save(fig, "e09_intersections")


# --------------------------------------------------------------------------- #
# E10 determinism
# --------------------------------------------------------------------------- #
def fig_e10_determinism() -> None:
    by_device = A.e10_by_device()
    if not by_device:
        raise Skip("no ok e10 records")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=WIDE, width_ratios=[2.2, 1])
    rng = np.random.default_rng(0)
    for offset, device, color in ((-0.17, "cpu", CPU), (0.17, "gpu", GPU)):
        records = by_device.get(device, [])
        if not records:
            note(f"e10: no {device} runs")
            continue
        for i, stage in enumerate(A.E10_STAGES):
            t = A.arr(r["stages"].get(stage) for r in records)
            if len(t) >= 3:
                ax.boxplot(t, positions=[i + offset], widths=0.26, showfliers=False, manage_ticks=False,
                           medianprops={"color": color}, boxprops={"color": color, "lw": 0.7},
                           whiskerprops={"color": color, "lw": 0.7}, capprops={"color": color, "lw": 0.7})
            ax.plot(i + offset + rng.uniform(-0.06, 0.06, len(t)), t, ls="", marker="o", ms=3, color=color,
                    alpha=0.8, label=device.upper() if i == 0 else None)
            c = A.cv(t)
            if c is not None:
                ax.annotate(f"{100 * c:.1f}%", (i + offset, np.nanmax(t)), xytext=(0, 4),
                            textcoords="offset points", ha="center", fontsize=6)
        totals = [sum(A.num(v) for v in r["stages"].values()) for r in records]
        bx.plot([r["config"].get("repeat", j) for j, r in enumerate(records)], totals, marker="o", ms=3.5,
                color=color, label=device.upper())
    ax.set_yscale("log")
    ax.set_xticks(range(len(A.E10_STAGES)), [s.replace("_", " ") for s in A.E10_STAGES])
    ax.set_ylabel("stage time (s)")
    ax.grid(axis="x", visible=False)
    ax.legend(loc="upper left", fontsize=7)
    ax.margins(y=0.25)
    reps = max(len(v) for v in by_device.values())
    panel(ax, "a", "stage times over repeats" + (" (labels: CV)" if reps >= 2 else " (one run each: no CV)"))
    bx.set_xticks(range(reps))
    bx.set_xlim(-0.5, reps - 0.5)
    bx.set_xlabel("repeat (fresh process)")
    bx.set_ylabel("total time (s)")
    bx.legend(loc="best", fontsize=7)
    panel(bx, "b", "drift")
    save(fig, "e10_determinism")


# --------------------------------------------------------------------------- #
# E11 deep blasting
# --------------------------------------------------------------------------- #
# A changed class, word or matrix is what a blast is for (author, 2026-10-07), so
# change is drawn in the blue ramp; red is kept for genuine errors only.
EVENT_COLORS = {"not_run": "white", "start": LIGHT, "none": GRID, "change": "#7fb0ee",
                "structure": "#0d366b", "error": RED}
EVENT_LABELS = {"not_run": "not run (stopped by its budget)", "start": "before the first blast",
                "none": "no change", "change": "change flags set, canonical words unchanged",
                "structure": "canonical words changed (new structure)",
                "error": "topology error (invariant assert / exception)"}


def e11_cases(records) -> list[str]:
    return [c for c in A.E11_CASES if any(r["config"]["case"] == c for r in records)]


def provisional_handle(**kw) -> Line2D:
    return Line2D([], [], ls="", marker="o", mfc="white", color=INK2,
                  label="symbolic result provisional (is_reliable=False)", **kw)


def error_ticks(ax, record, color=RED) -> None:
    """Blasts whose topology raised, as crosses on the x axis."""
    s = record["series"]
    bad = [b for b, o in zip(s["blast"], s["topology_outcome"]) if o != "ok"]
    ax.plot(bad, np.full(len(bad), 0.03), ls="", marker="x", ms=4, mew=1.0, color=color,
            transform=ax.get_xaxis_transform())


def growth_label(case: str, record: dict, keys: dict) -> str:
    """'case (crossings x1.61, ...)': growth factor per blast over the last five blasts."""
    parts = [f"{what} $\\times${f:.2f}" for what, key in keys.items()
             if (f := A.e11_growth(record["series"][key]))]
    return f"{case.replace('_', ' ')} ({', '.join(parts)} per blast)"


def fig_e11_evolution() -> None:
    refs = A.e11_reference()
    if not refs:
        raise Skip("no e11 reference cells")
    fig, axes = plt.subplots(1, 3, figsize=(6.3, 2.7), sharex=True)
    for case, r in refs.items():
        s, color = r["series"], CASE_COLORS[case]
        blast, events = A.arr(s["blast"]), A.e11_events(r)
        for axis, key in ((axes[0], "n_active"), (axes[1], "n_inert")):
            y = A.arr(s[key])
            axis.plot(blast, y, color=color, lw=1.0)
            for marker, kinds, fill in (("o", ("structure",), color), ("o", ("change",), "white")):
                pick = np.array([e in kinds for e in events])
                axis.plot(blast[pick], y[pick], ls="", marker=marker, ms=3.5, color=color, mfc=fill)
            error_ticks(axis, r)
        axes[2].plot(blast, A.arr(s["crossings"]), color=color, marker="o", ms=2.5, label=case)
    axes[2].set_yscale("log")
    axes[0].set_ylabel("active bridge classes")
    axes[1].set_ylabel("inert bridge classes")
    axes[2].set_ylabel("registered crossings")
    for axis, letter, title in zip(axes, "abc", ("active classes", "inert classes", "crossings")):
        axis.set_xlabel("blasts done")
        axis.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        panel(axis, letter, title)
    handles = [Line2D([], [], color=CASE_COLORS[c], marker="o", label=growth_label(c, r, {"crossings": "crossings"}))
               for c, r in refs.items()] + [
        Line2D([], [], ls="", marker="o", color=INK2, label="canonical words changed"),
        Line2D([], [], ls="", marker="o", mfc="white", color=INK2, label="change flags set, canonical words unchanged"),
        Line2D([], [], ls="", marker="x", color=RED, label="topology error (no snapshot)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, fontsize=6.5)
    save(fig, "e11_evolution")


def fig_e11_events() -> None:
    records = A.e11_records()
    if not records:
        raise Skip("no e11 records")
    width = max(len(r["series"]["blast"]) for r in records)
    codes = np.zeros((len(records), width), dtype=int)
    for i, r in enumerate(records):
        events = A.e11_events(r)
        codes[i, : len(events)] = [A.E11_EVENTS.index(e) for e in events]
    fig, ax = plt.subplots(figsize=(6.3, 0.105 * len(records) + 1.25))
    cmap = matplotlib.colors.ListedColormap([EVENT_COLORS[e] for e in A.E11_EVENTS])
    ax.imshow(codes, cmap=cmap, vmin=-0.5, vmax=len(A.E11_EVENTS) - 0.5, aspect="auto", interpolation="none")
    for i, r in enumerate(records):
        for j, flag in enumerate(A.e11_provisional(r)):
            if flag:
                ax.plot(j, i, ls="", marker="o", ms=2.2, mfc="white", mec=INK, mew=0.5)
        s = r["series"]
        for j, outcome in enumerate(s["topology_outcome"]):
            if outcome == "exception":
                ax.plot(j, i, ls="", marker="x", ms=3, color=INK, mew=0.9)
    ax.set_xticks(np.arange(width) - 0.5, minor=True)
    ax.set_yticks(np.arange(len(records)) - 0.5, minor=True)
    ax.grid(which="minor", color="white", lw=0.6)
    ax.grid(which="major", visible=False)
    ax.tick_params(which="minor", length=0)
    ax.set_xticks(range(0, width, 2))
    ax.set_yticks(range(len(records)), [A.e11_label(r) for r in records], fontsize=5.5)
    ax.set_xlabel("blast (0 = before blasting)")
    for case in e11_cases(records):
        rows = [i for i, r in enumerate(records) if r["config"]["case"] == case]
        ax.axhline(rows[-1] + 0.5, color=INK, lw=0.8)
        ax.text(-0.31, (rows[0] + rows[-1]) / 2, case.replace("_", " "), transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=7.5, color=CASE_COLORS[case], fontweight="bold")
    handles = [Patch(facecolor=EVENT_COLORS[e], edgecolor=INK2, lw=0.4, label=EVENT_LABELS[e]) for e in A.E11_EVENTS]
    handles += [provisional_handle(ms=3.5),
                Line2D([], [], ls="", marker="x", color=INK, ms=4, label="the error is an exception (not an assert)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, fontsize=6.5)
    ax.set_title("what each blast changed, per config", fontsize=9)
    save(fig, "e11_events")


def fig_e11_cost() -> None:
    records = A.e11_records()
    refs = A.e11_reference()
    if not refs:
        raise Skip("no e11 reference cells")
    fig, axes = plt.subplots(1, 3, figsize=(6.3, 2.7), sharex=True)
    for case, r in refs.items():
        s, color = r["series"], CASE_COLORS[case]
        blast = A.arr(s["blast"])
        for axis, key, scale in ((axes[0], "bridge_points", 1), (axes[1], "blast_s", 1), (axes[2], "rss_mb", 1 / 1024)):
            keep = blast > 0 if key == "blast_s" else blast >= 0
            axis.plot(blast[keep], A.arr(s[key])[keep] * scale, color=color, marker="o", ms=2.5, label=case)
        stop = r["metrics"].get("stop_reason")
        if stop != "max_blasts":
            axes[2].plot(blast[-1], A.num(s["rss_mb"][-1]) / 1024, ls="", marker="s", ms=5, mfc="white",
                         color=color, mew=1.1)
    budgets = {r["config"].get("_rss_budget_gb") for r in records} - {None}
    for gb in budgets:
        axes[2].axhline(gb, color=INK2, ls="--", lw=0.8)
        axes[2].text(0.3, gb, f"memory budget {gb:g} GB", va="bottom", fontsize=6, color=INK2)
    for axis in axes[:2]:
        axis.set_yscale("log")
    axes[0].set_ylabel("points on all bridges")
    axes[1].set_ylabel("time of the blast (s)")
    axes[2].set_ylabel("peak RSS (GB)")
    for axis, letter, title in zip(axes, "abc", ("bridge points", "time per blast", "memory")):
        axis.set_xlabel("blasts done")
        axis.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        panel(axis, letter, title)
    handles = [Line2D([], [], color=CASE_COLORS[c], marker="o",
                      label=growth_label(c, r, {"points": "bridge_points", "time": "blast_s"})) for c, r in refs.items()]
    handles += [Line2D([], [], ls="", marker="s", mfc="white", color=INK2,
                       label="stopped by its budget (next blast would not fit)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, fontsize=6.5)
    save(fig, "e11_cost")


E11_COLUMNS = ([("separation", v) for v in (None, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)]
               + [("cutoff", v) for v in (1e-5, 1e-6, 1e-8)])


def e11_divergence_cell(record, reference) -> tuple[str, str]:
    """(text, colour) of one parameter cell: '=' identical, else first diverging blast of words / crossings."""
    if record is reference:
        return "ref", LIGHT
    if A.e11_identical(record, reference):
        return "=", GRID
    words, crossings = A.e11_divergence(record, reference), A.e11_divergence(record, reference, "crossings")
    text = f"{'–' if words is None else words} / {'–' if crossings is None else crossings}"
    return text, "#7fb0ee" if words is not None else "#cde2fb"


def fig_e11_parameters() -> None:
    records, refs = A.e11_records(), A.e11_reference()
    if not refs:
        raise Skip("no e11 reference cells")
    fig = plt.figure(figsize=(6.3, 4.6))
    top, bottom = fig.subfigures(2, 1, height_ratios=[1, 1.2])
    ax = top.subplots()
    bx, cx = bottom.subplots(1, 2)
    cases = [c for c in A.E11_CASES if c in refs]
    for i, case in enumerate(cases):
        ref = refs[case]
        for j, (kind, value) in enumerate(E11_COLUMNS):
            if kind == "separation" and value == A.E11_SEPARATION[case]:
                record = ref
            else:
                record = next((r for r in records if r["config"]["case"] == case and A.e11_kind(r) == (kind, value)), None)
            text, color = ("", "white") if record is None else e11_divergence_cell(record, ref)
            ax.add_patch(matplotlib.patches.Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor=color, edgecolor="white", lw=1))
            ax.text(j, i, text, ha="center", va="center", fontsize=6.5,
                    color="white" if color == "#0d366b" else INK)
    ax.set_xlim(-0.5, len(E11_COLUMNS) - 0.5)
    ax.set_ylim(len(cases) - 0.5, -0.5)
    ax.set_xticks(range(len(E11_COLUMNS)), ["none" if v is None else f"{v:g}" for _, v in E11_COLUMNS],
                  fontsize=6.5)
    ax.set_yticks(range(len(cases)), [c.replace("_", " ") for c in cases], fontsize=7)
    ax.axvline(5.5, color=INK, lw=0.8)
    ax.set_xlabel("min_separation at cutoff $10^{-7}$" + " " * 30 + "area_cutoff while blasting", fontsize=7.5)
    ax.grid(False)
    panel(ax, "a", "first blast whose canonical words (both completed) / crossings differ from the reference "
                   "(= identical at every blast, – never)")

    for axis, case, letter in ((bx, "k28", "b"), (cx, "p3", "c")):
        mine = [r for r in records if r["config"]["case"] == case and A.e11_kind(r)[0] in ("canonical", "cutoff")]
        if not mine:
            empty(axis, f"no {case} cutoff cells")
            continue
        cutoffs = sorted({r["config"]["area_cutoff"] for r in mine}, reverse=True)
        colors = dict(zip(cutoffs, ramp(len(cutoffs))))
        for r in sorted(mine, key=lambda r: -r["config"]["area_cutoff"]):
            cutoff = r["config"]["area_cutoff"]
            axis.plot(A.arr(r["series"]["blast"]), A.arr(r["series"]["n_classes"]), color=colors[cutoff],
                      marker="o", ms=2.2, lw=1.0, label=f"{sci(cutoff)}" + (" (ref.)" if cutoff == A.E11_CUTOFF else ""))
            error_ticks(axis, r, color=RED)
        axis.set_xlabel("blasts done")
        axis.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        axis.legend(title="area_cutoff", title_fontsize=6, fontsize=6, loc="upper left")
        panel(axis, letter, f"{case}: bridge classes")
    bx.set_ylabel("bridge classes")
    save(fig, "e11_parameters")


DOUBLE_STYLE = {"same": (INK, "o", INK), "different": (INK, "o", "white"),
                "error both": (RED, "x", RED), "error one": (RED, "o", "white")}


def fig_e11_agreement() -> None:
    records, refs = A.e11_records(), A.e11_reference()
    doubles = [r for r in records if A.e11_kind(r)[0] == "iterations" and r["config"]["case"] in refs]
    if not doubles:
        raise Skip("no e11 num_iterations=2 cells")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=WIDE, width_ratios=[1.5, 1])
    lanes = []
    for r in doubles:
        case = r["config"]["case"]
        for row in A.e11_double(r, refs[case]):
            y = len(lanes)
            ax.plot(row["n"], y, ls="", marker="o", ms=4, color=INK, mfc=INK if row["crossings"] else "white")
            color, marker, fill = DOUBLE_STYLE[row["words"]]
            ax.plot(row["n"], y + 1, ls="", marker=marker, ms=4, color=color, mfc=fill)
        lanes += [f"{case.replace('_', ' ')}: crossings", f"{case.replace('_', ' ')}: words"]
    ax.set_yticks(range(len(lanes)), lanes, fontsize=6)
    ax.set_ylim(len(lanes) - 0.5, -0.5)
    ax.set_xlabel("blast $n$ with 2 iterations  (vs single-iteration blast $2n$)")
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.grid(axis="y", visible=False)
    for y in range(2, len(lanes), 2):
        ax.axhline(y - 0.5, color=GRID, lw=0.6)
    panel(ax, "a", "two iterations per blast")

    rows = []
    for case in [c for c in A.E11_CASES if c in refs]:
        for kind, what in (("separation", "separations"), ("repeat", "repeats")):
            mine = [r for r in records if r["config"]["case"] == case and A.e11_kind(r)[0] == kind]
            if mine:
                rows.append((f"{case.replace('_', ' ')}, {what}", sum(A.e11_identical(r, refs[case]) for r in mine), len(mine)))
    for i, (label, same, total) in enumerate(rows):
        bx.barh(i, total, color=GRID, height=0.62)
        bx.barh(i, same, color=BLUE, height=0.62)
        bx.text(total + 0.1, i, f"{same}/{total}", va="center", fontsize=6.5)
    bx.set_yticks(range(len(rows)), [r[0] for r in rows], fontsize=6)
    bx.set_ylim(len(rows) - 0.5, -0.5)
    bx.set_xlim(0, max(r[2] for r in rows) + 1.2)
    bx.set_xlabel("cells identical to the reference")
    bx.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    bx.grid(axis="y", visible=False)
    panel(bx, "b", "determinism and min_separation")
    handles = [Line2D([], [], ls="", marker="o", color=INK, label="same"),
               Line2D([], [], ls="", marker="o", mfc="white", color=INK, label="different"),
               Line2D([], [], ls="", marker="x", color=RED, label="topology error in both")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, fontsize=6.5)
    save(fig, "e11_agreement")


def fig_e11_provisional() -> None:
    refs = A.e11_reference()
    if not refs:
        raise Skip("no e11 reference cells")
    fig, axes = plt.subplots(1, 3, figsize=(6.3, 2.6), sharex=True, sharey=True)
    for i, (case, r) in enumerate(refs.items()):
        s, color = r["series"], CASE_COLORS[case]
        blast, provisional = A.arr(s["blast"]), np.array([bool(p) for p in A.e11_provisional(r)])
        offset = 0.08 * (i - (len(refs) - 1) / 2)  # integer counts: cases side by side, not on top
        for axis, key in zip(axes, ("n_unresolved", "n_virtual", "n_ambiguous")):
            y = A.arr(s[key]) + offset
            axis.plot(blast, y, color=color, lw=1.0, label=case)
            axis.plot(blast[provisional], y[provisional], ls="", marker="o", ms=3, color=color, mfc="white")
            axis.plot(blast[~provisional], y[~provisional], ls="", marker="o", ms=3, color=color)
            error_ticks(axis, r)
    axes[0].set_ylabel("classes (cases offset by $\\pm$0.16)")
    for axis, letter, title in zip(axes, "abc", ("unresolved classes", "virtual classes (new1, ...)",
                                                 "ambiguous walks")):
        axis.set_xlabel("blasts done")
        axis.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        axis.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        panel(axis, letter, title)
    handles = case_handles(refs) + [Line2D([], [], ls="", marker="o", color=INK2, label="is_reliable=True"),
                                    provisional_handle(),
                                    Line2D([], [], ls="", marker="x", color=RED, label="topology error (no snapshot)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, fontsize=6.5)
    save(fig, "e11_provisional")


# --------------------------------------------------------------------------- #
# Overview
# --------------------------------------------------------------------------- #
def scorecard() -> list[tuple[str, float, float]]:
    """(check, passed, total) for the robustness experiments, each counted as stated."""
    out = []
    rows = [r for r in A.e06_rows() if r["family"] == "b = 1, saddle"]
    if rows:
        out.append(("E6 b = 1 saddles complete", sum(r["stage"] == "trellis_pips" for r in rows), len(rows)))
    table = A.e08_table()
    if table:
        n = sum(r["n"] for r in table.values())
        out.append(("E8 solves converge", sum(r["converged"] for r in table.values()), n))
        conv = sum(r["converged"] for r in table.values())
        out.append(("E8 converged at requested period", conv - sum(r["lower_period"] for r in table.values()), conv))
    e9 = [r for r in A.e09_rows() if r["outcome"] == "ok"]
    if e9:
        out.append(("E9 brute crossings matched", sum(r["matched"] for r in e9), sum(r["brute"] for r in e9)))
    e5 = A.e05_cases()
    if e5:
        steps = sum(A.num((r["metrics"].get("cdist_manifolds") or {}).get("steps", 0)) for r in e5)
        bad = sum(A.num((r["metrics"].get("cdist_manifolds") or {}).get("decreasing", 0)) for r in e5)
        if steps:
            out.append(("E5 cdist steps non-decreasing", steps - bad, steps))
    blasts = A.e07("blasts")
    if blasts:
        out.append(("E7 blasts completed", sum(A.num(r["metrics"].get("blasts_done", 0)) for r in blasts),
                    sum(r["config"].get("max_blasts", 0) for r in blasts)))
    deep = A.e11_records()
    if deep:
        blasts = sum(len(r["series"]["blast"]) - 1 for r in deep)
        laws = sum(s > 0 or c > 0 for r in deep for s, c in
                   zip(r["series"]["same_stability"][1:], r["series"]["cdist_decreasing"][1:]))
        out.append(("E11 blasts with no crossing/cdist violation", blasts - laws, blasts))
        snapshots = sum(len(r["series"]["blast"]) for r in deep)
        out.append(("E11 topology snapshots without error", sum(sum(A.e11_ok(r)) for r in deep), snapshots))
    ident = A.e10_identity(A.e10_by_device())
    pairs = [v["rounded"]["CPU-GPU"] for v in ident.values() if v["rounded"]["CPU-GPU"] is not None]
    if pairs:
        out.append(("E10 CPU = GPU hashes (1e-12)", sum(pairs), len(pairs)))
    return out


def fig_overview() -> None:
    fig, axes = plt.subplots(2, 2, figsize=(6.3, 4.6))
    (ax, bx), (cx, dx) = axes
    sweeps = A.e01_sweeps()
    if sweeps:
        for backend in ("numpy", "cupy_e2e", "map_batch_gpu"):
            if backend in sweeps:
                color, ls, marker, filled, label = E1_STYLE[backend]
                ax.plot(sweeps[backend]["n"], sweeps[backend]["points_per_s"], ls=ls, marker=marker, ms=3,
                        color=color, mfc=color if filled else "white", label=label)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("batch size $N$")
        ax.set_ylabel("points/s")
        ax.legend(loc="lower right", fontsize=6.5)
        panel(ax, "a", "map throughput (E1)")
    else:
        empty(ax, "E1")

    totals = A.e02_total_speedups(A.e02_timed())
    totals = {c: v for c, v in totals.items() if v}
    if totals:
        cases = [c for c in A.CASES if c in totals]
        for i, case in enumerate(cases):
            v = totals[case]
            bx.bar(i, np.median(v), color=CASE_COLORS[case], width=0.6)
            bx.plot(np.full(len(v), i), v, ls="", marker="o", ms=3, color=INK)
            bx.annotate(f"{np.median(v):.2f}$\\times$", (i, max(v)), xytext=(0, 4), textcoords="offset points",
                        ha="center", fontsize=7)
        bx.axhline(1, color=INK, lw=0.7)
        bx.set_ylim(0, 1.2 * max(max(v) for v in totals.values()))
        bx.set_xticks(range(len(cases)), cases)
        bx.set_ylabel("CPU time / GPU time")
        bx.grid(axis="x", visible=False)
        panel(bx, "b", "end-to-end GPU speed-up (E2)")
    else:
        empty(bx, "E2")

    runs = A.e04_runs()
    if runs:
        for (case, steps), rows in runs.items():
            st = e4_style(case, steps, runs)
            good = [r for r in rows if r["outcome"] == "ok" and r["inv_p50"] is not None]
            cx.plot([r["cutoff"] for r in good], [r["inv_p50"] for r in good], color=st["color"], marker=st["marker"],
                    ms=3, label=st["label"])
            cx.plot([r["cutoff"] for r in good], [r["inv_p99"] for r in good], color=st["color"], marker=st["marker"],
                    ms=3, mfc="white", ls=":", lw=0.9)
        cx.set_xscale("log")
        cx.set_yscale("log")
        cx.invert_xaxis()
        cx.set_xlabel("area_cutoff")
        cx.set_ylabel("invariance error")
        cx.legend(loc="upper right", fontsize=6, title="p50 solid, p99 dotted", title_fontsize=6)
        panel(cx, "c", "accuracy knob (E4)")
    else:
        empty(cx, "E4")

    card = scorecard()
    if card:
        for i, (label, passed, total) in enumerate(card):
            frac = passed / total if total else 0
            dx.barh(i, 1, color=GRID, height=0.62)
            dx.barh(i, frac, color=GREEN, height=0.62)
            dx.text(1.02, i, f"{passed:.0f}/{total:.0f}", va="center", fontsize=6.5,
                    transform=dx.get_yaxis_transform())
        dx.set_yticks(range(len(card)), [c[0] for c in card], fontsize=6.5)
        dx.invert_yaxis()
        dx.set_xlim(0, 1)
        dx.set_xlabel("fraction passing")
        dx.grid(axis="y", visible=False)
        panel(dx, "d", "robustness scorecard")
    else:
        empty(dx, "E5-E10")
    save(fig, "overview")


FIGURE_FUNCTIONS = [
    fig_e01_map_throughput,
    fig_e02_end_to_end, fig_e02_profile, fig_e02_min_batch, fig_e02_deviation,
    fig_e03_growth, fig_e03_cost, fig_e03_memory,
    fig_e04_pareto, fig_e04_convergence,
    fig_e05_invariance, fig_e05_iterates_lobes,
    fig_e06_sweep_outcomes, fig_e06_saddle_sweep,
    fig_e07_escape, fig_e07_blasts, fig_e07_depth, fig_e07_tangency,
    fig_e08_solver_outcomes, fig_e08_orbits, fig_e08_solve_time,
    fig_e09_intersections,
    fig_e10_determinism,
    fig_e11_evolution, fig_e11_events, fig_e11_cost, fig_e11_parameters, fig_e11_agreement,
    fig_e11_provisional,
    fig_overview,
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--quick", action="store_true", help="use the *_quick.jsonl results")
    parser.add_argument("--only", type=str, default=None, help="comma list of name fragments, e.g. e04,overview")
    args = parser.parse_args()
    A.QUICK = args.quick
    FIGURES.mkdir(parents=True, exist_ok=True)
    PNG.mkdir(parents=True, exist_ok=True)
    style()
    logging.getLogger("fontTools").setLevel(logging.ERROR)  # font subsetting chatter on PDF save
    failed = []
    for fn in FIGURE_FUNCTIONS:
        if args.only and not any(part in fn.__name__ for part in args.only.split(",")):
            continue
        print(f"{fn.__name__}")
        try:
            fn()
        except Skip as exc:
            note(f"{fn.__name__} skipped: {exc}")
        except Exception:  # noqa: BLE001 - one broken figure never stops the rest
            failed.append(fn.__name__)
            traceback.print_exc()
            plt.close("all")
    print("\nsources:", {k: v[0] + f" ({v[2]} records)" for k, v in sorted(A.SOURCES.items())})
    if SKIPPED:
        print(f"{len(SKIPPED)} skips/notes (listed above)")
    if failed:
        print("FAILED:", ", ".join(failed))
        sys.exit(1)


if __name__ == "__main__":
    main()
