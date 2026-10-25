"""
Figures of ``profile_numerics.py``: drawn from its json files alone.

``draw(out, compare=None)`` writes into ``out``; with ``compare`` (a baseline
output dir) every figure draws the baseline beside ``out``'s run (bars side by
side, baseline lines dashed), so the last slice of the campaign compares on one
axes.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

#: The eight categorical slots, fixed order (dataviz reference palette, light).
SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"

#: Census components folded into eight groups, in stacking order.
MEMORY_GROUPS = [
    ("Point on manifolds", ["points_curve", "branch_points"]),
    ("Point on blast images", ["points_image", "points_other"]),
    ("point store + phantoms", ["point_store", "points_phantom"]),
    ("segment index (py)", ["segment_index"]),
    ("rtree (C)", ["rtree"]),
    ("curve objects + ids + memo", ["curve_objects", "curve_ids", "walk_memo"]),
    ("registry", ["registry"]),
    ("topology caches", ["topology"]),
]
#: Inner stage labels folded into eight groups.
TIME_GROUPS = [
    ("map", ["map"]),
    ("refine", ["refine"]),
    ("merge + image chains", ["merge", "image_chains", "iterate"]),
    ("walks", ["walk"]),
    ("segment index + resolve", ["index", "resolve", "cut"]),
    ("register + iterates", ["register", "register_bridge", "infer_iterates", "resolved_ids"]),
    ("diagnose + separation", ["diagnose", "separation"]),
    ("other", None),
]
STAGE_ORDER = ["fixed_point", "initialise", "grow", "intersections", "trim", "create_bridges",
               "pips", "zone", "blast", "topology"]
CASE_ORDER = ["k10", "p3", "nested", "inversion", "heteroclinic", "scale"]


def _style(ax) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(axis="both", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def load(out: Path) -> dict:
    """``{(case, mode): record}`` of every json in ``out``."""
    runs = {}
    for path in sorted(out.glob("*_*.json")):
        data = json.loads(path.read_text())
        runs[(data["case"], data["mode"])] = data
    return runs


def memory_groups(census: dict) -> dict[str, float]:
    """Census bytes folded into :data:`MEMORY_GROUPS`."""
    b = census["bytes"]
    return {name: float(sum(b.get(k, 0) for k in keys)) for name, keys in MEMORY_GROUPS}


def time_groups(by_stage: dict[str, float]) -> dict[str, float]:
    """One stage's self times folded into :data:`TIME_GROUPS`."""
    known = {k for _, keys in TIME_GROUPS if keys for k in keys}
    out = {}
    for name, keys in TIME_GROUPS:
        if keys is None:
            out[name] = sum(v for k, v in by_stage.items() if k not in known)
        else:
            out[name] = sum(by_stage.get(k, 0.0) for k in keys)
    return out


def _legend(fig, names, colors) -> None:
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colors]
    fig.legend(handles, names, loc="lower center", ncol=4, frameon=False, fontsize=8,
               labelcolor=INK)


def figure_stage_time(runs: dict, base: dict, out: Path) -> Path:
    """Per case: time per top-level stage, stacked by inner label (self time)."""
    cases = [c for c in CASE_ORDER if (c, "time") in runs]
    fig, axes = plt.subplots(len(cases), 1, figsize=(9, 2.1 * len(cases) + 1.0), squeeze=False)
    for ax, case in zip(axes[:, 0], cases):
        sets = [("", runs[(case, "time")])]
        if (case, "time") in base:
            sets = [("base ", base[(case, "time")]), ("", runs[(case, "time")])]
        by_stage_all = [s[1]["timer"]["by_stage"] for s in sets]
        stages = [s for s in STAGE_ORDER if any(s in b for b in by_stage_all)]
        labels, rows = [], []
        for stage in stages:
            for prefix, rec in sets:
                labels.append(prefix + stage)
                rows.append(time_groups(rec["timer"]["by_stage"].get(stage, {})))
        y = np.arange(len(labels))[::-1]
        left = np.zeros(len(labels))
        for (name, _), color in zip(TIME_GROUPS, SLOTS):
            widths = np.array([r[name] for r in rows])
            ax.barh(y, widths, left=left, color=color, height=0.7, edgecolor="white",
                    linewidth=1)
            left += widths
        for yi, total in zip(y, left):
            ax.text(total, yi, f" {total:.2f} s", va="center", fontsize=7, color=MUTED)
        ax.set_yticks(y, labels, fontsize=8)
        wall = runs[(case, "time")]["wall_s"]
        ax.set_title(f"{case}: wall {wall:.1f} s (self time by inner stage)", fontsize=9,
                     loc="left", color=INK)
        ax.set_xlim(0, left.max() * 1.18 if left.max() else 1)
        _style(ax)
    axes[-1, 0].set_xlabel("seconds", fontsize=8, color=MUTED)
    _legend(fig, [n for n, _ in TIME_GROUPS], SLOTS)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    path = out / "stage_time.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def figure_memory(runs: dict, base: dict, out: Path) -> Path:
    """Per case: census components at the end of the build (stacked MB)."""
    cases = [c for c in CASE_ORDER if (c, "time") in runs and "census" in runs[(c, "time")]]
    labels, rows, refs = [], [], []
    for case in cases:
        sets = [("", runs)] if (case, "time") not in base else [("base ", base), ("", runs)]
        for prefix, src in sets:
            rec = src[(case, "time")]
            census = dict(rec["census"])
            census["bytes"] = dict(census["bytes"])
            if "rtree_measured_bytes" in census:
                census["bytes"]["rtree"] = max(census["rtree_measured_bytes"], 0)
            labels.append(f"{prefix}{case}\n{census['n_points']:,} pts")
            rows.append(memory_groups(census))
            mem = src.get((case, "memory"), {})
            # A memory run that stopped at another point count (the scaling case's
            # tracemalloc run is shorter by default) is no reference for this bar.
            if mem.get("census", {}).get("n_points") != census["n_points"]:
                mem = {}
            refs.append((
                (rec.get("rss_final", 0) - rec["rss_baseline"]) / 2**20,
                ((mem.get("traced_final") or 0) - (mem.get("traced_baseline") or 0)) / 2**20
                if mem else None,
            ))
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [3, 2]})
    x = np.arange(len(labels))
    bottom = np.zeros(len(labels))
    for (name, _), color in zip(MEMORY_GROUPS, SLOTS):
        h = np.array([r[name] for r in rows]) / 2**20
        ax.bar(x, h, bottom=bottom, color=color, width=0.7, edgecolor="white", linewidth=1)
        bottom += h
    for xi, (rss, traced) in zip(x, refs):
        ax.plot([xi - 0.38, xi + 0.38], [rss, rss], color=INK, linewidth=1.5)
        if traced is not None:
            ax.plot([xi - 0.38, xi + 0.38], [traced, traced], color=INK, linewidth=1.5,
                    linestyle=":")
    ax.set_xticks(x, labels, fontsize=7)
    ax.set_yscale("log")
    ax.set_ylabel("MB held at the end (log)", fontsize=8, color=MUTED)
    ax.set_title("memory by component; solid: RSS growth, dotted: tracemalloc growth",
                 fontsize=9, loc="left", color=INK)
    _style(ax)
    # shares
    shares = np.array([[r[n] for n, _ in MEMORY_GROUPS] for r in rows])
    shares = shares / shares.sum(axis=1, keepdims=True)
    left = np.zeros(len(labels))
    for j, color in enumerate(SLOTS):
        ax2.barh(x[::-1], shares[:, j], left=left, color=color, height=0.7, edgecolor="white",
                 linewidth=1)
        left += shares[:, j]
    ax2.set_yticks(x[::-1], [lb.split("\n")[0] for lb in labels], fontsize=7)
    ax2.set_xlim(0, 1)
    ax2.set_title("share of the census", fontsize=9, loc="left", color=INK)
    _style(ax2)
    _legend(fig, [n for n, _ in MEMORY_GROUPS], SLOTS)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    path = out / "memory_components.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def figure_stage_peaks(runs: dict, out: Path) -> Path:
    """Per case: the RSS peak inside each top-level stage call (transient peaks)."""
    cases = [c for c in CASE_ORDER if (c, "time") in runs and c != "scale"]
    fig, ax = plt.subplots(figsize=(9, 3.6))
    for case, color in zip(cases, SLOTS):
        rec = runs[(case, "time")]
        ev = rec["timer"]["top_events"]
        peaks = np.array([e["peak_rss"] - rec["rss_baseline"] for e in ev]) / 2**20
        after = np.array([e["rss_after"] - rec["rss_baseline"] for e in ev]) / 2**20
        ax.plot(np.arange(len(ev)), peaks, color=color, linewidth=2, label=case)
        ax.plot(np.arange(len(ev)), after, color=color, linewidth=1, linestyle=":")
    ax.set_xlabel("top-level stage call, in order", fontsize=8, color=MUTED)
    ax.set_ylabel("MB over baseline", fontsize=8, color=MUTED)
    ax.set_title("peak RSS inside each stage call (solid) and after it (dotted)", fontsize=9,
                 loc="left", color=INK)
    ax.legend(frameon=False, fontsize=8)
    _style(ax)
    fig.tight_layout()
    path = out / "stage_peak_rss.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def _series(rec: dict) -> list[dict]:
    return [r for r in rec.get("series", []) if "census" in r]


def figure_scaling(runs: dict, base: dict, out: Path) -> list[Path]:
    """The scaling series: bytes per point, time and RSS per blast, points per blast."""
    if ("scale", "time") not in runs:
        return []
    sets = [("this run", runs[("scale", "time")], "-")]
    if ("scale", "time") in base:
        sets.insert(0, ("base", base[("scale", "time")], "--"))
    paths = []

    # bytes per point
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 4.4))
    for label, rec, ls in sets:
        rows = _series(rec)
        pts = np.array([r["census"]["n_points"] for r in rows], float)
        rss = np.array([r["rss"] - rec["rss_baseline"] for r in rows], float)
        peak = np.array([r["step_peak_rss"] - rec["rss_baseline"] for r in rows], float)
        groups = np.array([[memory_groups(r["census"])[n] for n, _ in MEMORY_GROUPS] for r in rows])
        ax.plot(pts, rss / pts, color=INK, linestyle=ls, linewidth=2, label=f"RSS / point ({label})")
        ax.plot(pts, peak / pts, color=MUTED, linestyle=ls, linewidth=1.5,
                label=f"step peak RSS / point ({label})")
        ax.plot(pts, groups.sum(axis=1) / pts, color=SLOTS[0], linestyle=ls, linewidth=2,
                label=f"census / point ({label})")
        if ls == "-":
            ax2.stackplot(pts, (groups / pts[:, None]).T, colors=SLOTS, edgecolor="white",
                          linewidth=0.5)
    ax.set_xscale("log")
    ax.set_xlabel("points held", fontsize=8, color=MUTED)
    ax.set_ylabel("bytes per point", fontsize=8, color=MUTED)
    ax.set_title("bytes per point over the blast series", fontsize=9, loc="left", color=INK)
    ax.legend(frameon=False, fontsize=7)
    _style(ax)
    ax2.set_xscale("log")
    ax2.set_xlabel("points held", fontsize=8, color=MUTED)
    ax2.set_title("census bytes per point, by component", fontsize=9, loc="left", color=INK)
    _style(ax2)
    _legend(fig, [n for n, _ in MEMORY_GROUPS], SLOTS)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    paths.append(out / "scaling_bytes_per_point.png")
    fig.savefig(paths[-1], dpi=130)
    plt.close(fig)

    # time and RSS per blast (two charts, one axis each)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 4.0))
    for label, rec, ls in sets:
        rows = rec["series"]
        b = np.array([r["blast"] for r in rows])
        ax.plot(b, [r["blast_s"] for r in rows], color=SLOTS[0], linestyle=ls, linewidth=2,
                label=f"blast_zone ({label})")
        ax.plot(b, [r["step_s"] - r["blast_s"] for r in rows], color=SLOTS[1], linestyle=ls,
                linewidth=2, label=f"topology() ({label})")
        ax2.plot(b, np.array([r["step_peak_rss"] for r in rows]) / 2**30, color=SLOTS[0],
                 linestyle=ls, linewidth=2, label=f"peak RSS in the step ({label})")
        ax2.plot(b, np.array([r["rss"] for r in rows]) / 2**30, color=SLOTS[1], linestyle=ls,
                 linewidth=2, label=f"RSS after the step ({label})")
    ax.set_xlabel("blast", fontsize=8, color=MUTED)
    ax.set_ylabel("seconds per step", fontsize=8, color=MUTED)
    ax.set_yscale("log")
    ax.legend(frameon=False, fontsize=7)
    ax.set_title("wall time per blast step", fontsize=9, loc="left", color=INK)
    _style(ax)
    ax2.set_xlabel("blast", fontsize=8, color=MUTED)
    ax2.set_ylabel("GB", fontsize=8, color=MUTED)
    ax2.legend(frameon=False, fontsize=7)
    ax2.set_title("resident memory per blast step", fontsize=9, loc="left", color=INK)
    _style(ax2)
    fig.tight_layout()
    paths.append(out / "scaling_time_rss.png")
    fig.savefig(paths[-1], dpi=130)
    plt.close(fig)

    # points per blast
    fig, ax = plt.subplots(figsize=(8, 4.0))
    kinds = [("n_curve_points", "on manifolds"), ("n_image_points", "on blast images"),
             ("n_phantom", "phantom preimages")]
    for label, rec, ls in sets:
        rows = _series(rec)
        b = np.array([r["blast"] for r in rows])
        for (key, name), color in zip(kinds, SLOTS):
            ax.plot(b, [r["census"][key] for r in rows], color=color, linestyle=ls,
                    linewidth=2, label=f"{name} ({label})")
    ax.set_yscale("log")
    ax.set_xlabel("blast", fontsize=8, color=MUTED)
    ax.set_ylabel("Point objects", fontsize=8, color=MUTED)
    ax.legend(frameon=False, fontsize=7)
    ax.set_title("points per blast", fontsize=9, loc="left", color=INK)
    _style(ax)
    fig.tight_layout()
    paths.append(out / "scaling_points.png")
    fig.savefig(paths[-1], dpi=130)
    plt.close(fig)
    return paths


def draw(out: Path, compare: Optional[Path] = None) -> list[Path]:
    """Every figure of the runs in ``out`` (beside ``compare``'s when given)."""
    runs = load(out)
    base = load(compare) if compare else {}
    paths = []
    if any(m == "time" for _, m in runs):
        paths.append(figure_stage_time(runs, base, out))
        paths.append(figure_memory(runs, base, out))
        paths.append(figure_stage_peaks(runs, out))
        paths += figure_scaling(runs, base, out)
    for p in paths:
        print("wrote", p)
    return paths
