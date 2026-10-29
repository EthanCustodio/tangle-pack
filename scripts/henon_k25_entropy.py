"""
Topological entropy of the Hénon tangles at ``k = 2.5``, ``b = 1`` as a function
of the number of blasts (the maximum iterate of the zone's bridges).

Three cases, each built by its own k=2.5 pipeline driver with its blast count
set to zero, then blasted ONE single-iteration step at a time; the topology is
recomputed after every step (grow, blast and recompute commute, design note
§A.10):

* ``p1`` -- the period-1 saddle alone (``henon_k25_pipeline.py``: ten
  unstable steps, stable to turnaround, the zone cut at f(q0));
* ``p3`` -- the inner period-3 orbit alone (``henon_k25_period3_pipeline.py``:
  continued from k=2, 9 / 6 steps, the default strong pip);
* ``nested`` -- both together (``henon_k25_nested_pipeline.py``: outer pip
  f(q0), inner default); one step blasts the inner zone, then the outer.

At each step the entropy is ``h = ln rho(M)``, ``rho`` the spectral radius of
the transition matrix over the active symbols, both the refined matrix ``M_1``
(the library's default) and the class matrix ``M``. ``M_1`` is an in-splitting
of ``M``, so the two radii agree whenever Theorem T holds: a disagreement is a
flag. ``h`` is the entropy of the trellis's symbolic dynamics, a LOWER bound
on the map's topological entropy that rises as blasting exposes more of the
tangle; ``ln 2`` (the full horseshoe) bounds it from above. For the nested case
the entropy of each strongly connected component of the transition graph is
recorded too, keyed by the tangle(s) its symbols lie in: the whole matrix's
radius is the largest component's.

Each (case, cutoff) sweep writes one JSON file after every step (so a long or
interrupted run keeps what it computed) into
``figures/henon_k25_entropy/data/<case>_a<cutoff>.json``; the figures are drawn
from those files alone (``--plot-only`` redraws them without computing):

* ``entropy`` -- h against blasts, one panel per case, one line per cutoff
  (filled markers: ``is_reliable``; hollow: not), refined solid, class dashed;
* ``entropy_vs_depth`` -- h against the pip depth m, every cutoff and blast
  one point (the cutoff acts through the depth the pip orbit resolves to);
* ``entropy_overlay`` -- the three cases on one axes at the finest cutoff, the
  nested components drawn beside the cases alone;
* ``nested_components`` -- the nested run's per-component entropies against
  the p1-alone and p3-alone curves;
* ``alphabet`` -- active / inert / refined symbols, classes, minimal-trellis
  and trellis bridges, and the pip depth against blasts;
* ``health`` -- unresolved classes, ambiguous classes, unresolved detections
  (not crossings of N) and the M / M_1 radius gap against blasts;
* ``cost`` -- manifold points and wall time against blasts (log scale);
* ``convergence`` -- h at each blast count against the cutoff;
* ``length_growth`` -- the arc length of each blast's interior frontier and
  its per-blast growth rate beside h (a geometric cross-check: the length of
  the iterated unstable curve inside the zone grows like ``e^{h n}``);
* ``matrices`` -- the last step's refined transition matrix per case;
* ``spectra`` -- the last step's eigenvalues per case in the complex plane.

Every figure is written as PNG, SVG and a pickled figure (reopen a pickle with
``scripts/show_figure.py``).

Run with ``PYTHONPATH=src python3 scripts/henon_k25_entropy.py`` (``--cases p1
p3 nested``, ``--cutoffs 1e-7 1e-8``, ``--max-blasts 12``, ``--max-seconds``
per sweep, ``--plot-only``). The sweeps of one call run one after the other;
for parallel runs launch one process per case and cutoff and finish with
``--plot-only``.

Dev Notes:
    A step that raises (``PointBudgetExceeded``, a pip that stops being a
    candidate, an assertion) ends that sweep: the error is recorded in the
    JSON (``stopped``) and the figures show the steps before it.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import sys
import time
import traceback
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import networkx
import numpy as np
from matplotlib.lines import Line2D
from numpy.typing import NDArray

sys.path.insert(0, str(Path(__file__).resolve().parent))

import henon_k25_nested_pipeline as nested_driver  # noqa: E402
import henon_k25_period3_pipeline as p3_driver  # noqa: E402
from henon_bridge_classes import Case  # noqa: E402
from henon_bridge_classes import build as build_p1  # noqa: E402
from symbolic_figures import FIGURES_DIR, configure_logging, save_panel  # noqa: E402

from tanglepack import TangleSession  # noqa: E402
from tanglepack.topology.SymbolicDynamics import SymbolicDynamics  # noqa: E402

logger = logging.getLogger(__name__)

K = 2.5
B = 1
CASES = ("p1", "p3", "nested")
CASE_TITLES = {
    "p1": "period 1 (outer) alone",
    "p3": "period 3 (inner) alone",
    "nested": "nested: period 1 + period 3",
}
CASE_COLORS = {"p1": "#1f77b4", "p3": "#d62728", "nested": "#2ca02c"}
MIN_SEPARATION = 1e-5
#: The blast's same-stability count (``blast_zone(diagnose=)``); an unpruned sweep turns it off (~2.3x memory per step).
BLAST_DIAGNOSE = True
#: End a sweep after the first blast that keeps no new interior bridge in any zone (nothing left to iterate).
STOP_WHEN_EXHAUSTED = False
#: The folder the figures go in, and their data under it.
OUT_DIR = FIGURES_DIR / "henon_k25_entropy"
DATA_DIR = OUT_DIR / "data"


# ── building ────────────────────────────────────────────────────────────────


def build_case(case: str, area_cutoff: float) -> tuple[TangleSession, list, list]:
    """
    Build one case, unblasted, through its pipeline driver.

    Args:
        case: ``"p1"``, ``"p3"`` or ``"nested"``.
        area_cutoff: The refinement area cutoff.

    Returns:
        ``(session, fixed points outermost first, zones in blast order)``; a
        nested step blasts the inner zone first, as the nested driver does.
    """
    if case == "p1":
        session, fp, _ = build_p1(
            Case(k=K, unstable_steps=10, pip="image", area_cutoff=area_cutoff, blasts=0)
        )
        fixed_points = [fp]
    elif case == "p3":
        session, fp, _ = p3_driver.build(
            unstable_steps=13, stable_steps=9, blasts=0, scale=True, area_cutoff=area_cutoff
        )
        fixed_points = [fp]
    elif case == "nested":
        session, fixed_points, _ = nested_driver.build(
            outer_blasts=0, inner_blasts=0, area_cutoff=area_cutoff
        )
    else:
        raise ValueError(f"unknown case {case!r}")
    zones = sorted(session.resonance_zones.values(), key=lambda zone: zone.area)
    return session, fixed_points, zones


# ── entropy ─────────────────────────────────────────────────────────────────


def spectral_radius(matrix: NDArray[np.int64]) -> float:
    """The spectral radius of a square matrix (0 for an empty one)."""
    if matrix.size == 0:
        return 0.0
    return float(np.max(np.abs(np.linalg.eigvals(matrix.astype(float)))))


def entropy_of(radius: float) -> float:
    """``ln rho``, 0 when the radius is at most 1 (no exponential growth)."""
    return math.log(radius) if radius > 1.0 else 0.0


def components(dynamics: SymbolicDynamics, names: list[str], matrix: NDArray[np.int64]) -> list[dict]:
    """
    The strongly connected components of the refined transition graph that carry a cycle.

    Args:
        dynamics: The symbolic dynamics.
        names: The refined symbols, in matrix order.
        matrix: The refined transition matrix.

    Returns:
        One dict per component with its ``symbols``, the ``tangles`` its
        classes connect (the periods of the stable branches a class's two
        elements lie on: ``"p1"``, ``"p3"``, or ``"p1-p3"`` for a heteroclinic
        class) and its ``radius``
        and ``entropy``, largest radius first.
    """
    tangle_of = {}
    for node in dynamics.symbol_nodes(refined=True):
        ends = node.bridge_class.source, node.bridge_class.target
        tangle_of[node.name] = "-".join(sorted({f"p{ref.branch_key[0].period}" for ref in ends}))
    graph = networkx.DiGraph()
    graph.add_nodes_from(range(len(names)))
    graph.add_edges_from(zip(*np.nonzero(matrix)))
    rows = []
    for component in networkx.strongly_connected_components(graph):
        index = sorted(component)
        block = matrix[np.ix_(index, index)]
        if len(index) == 1 and block[0, 0] == 0:
            continue
        radius = spectral_radius(block)
        rows.append({
            "symbols": [names[i] for i in index],
            "tangles": sorted({tangle_of[names[i]] for i in index}),
            "radius": radius,
            "entropy": entropy_of(radius),
        })
    return sorted(rows, key=lambda row: -row["radius"])


def measure(session: TangleSession, fixed_points: list, blasts: int, seconds: float,
            pip_depth: dict, blast_points: int) -> dict:
    """
    Recompute the topology and read one step's row.

    Args:
        session: The session, blasted ``blasts`` times.
        fixed_points: Its fixed points.
        blasts: The blast count.
        seconds: Wall time spent so far in this sweep.
        pip_depth: The last blast's ``pip_depth``, by period.
        blast_points: Points the blasts have added (image curves are no
            registered manifold).

    Returns:
        The row: entropies, alphabet sizes, health counts, cost and the refined
        matrix.
    """
    topology = session.topology()
    dynamics = topology.dynamics
    names, refined = dynamics.transition_matrix(refined=True)
    letters, coarse = dynamics.transition_matrix(refined=False)
    rho_refined, rho_class = spectral_radius(refined), spectral_radius(coarse)
    workbench = session.workbench
    registry = workbench.intersection_registry
    trellis = topology.trellis
    eigenvalues = np.linalg.eigvals(refined.astype(float)) if refined.size else np.array([])
    kinds = [cd.kind for cd in dynamics.classes.values()]
    return {
        "blasts": blasts,
        "entropy": entropy_of(rho_refined),
        "entropy_class": entropy_of(rho_class),
        "radius": rho_refined,
        "radius_class": rho_class,
        "reliable": bool(dynamics.is_reliable),
        "classes": len(dynamics.classes),
        "active": kinds.count("active"),
        "inert": kinds.count("inert"),
        "refined_symbols": len(names),
        "unresolved_classes": len(dynamics.unresolved),
        "ambiguous_classes": sum(bool(cd.ambiguous) for cd in dynamics.classes.values()),
        "unspanned_classes": sum(not cd.entry.spanned for cd in dynamics.classes.values()),
        "minimal_bridges": len(topology.minimal.kept_bridge_ids),
        "trellis_bridges": len(trellis.bridges),
        "crossings": len(registry.resolved_ids(workbench.resolution_ratio_threshold)),
        "unresolved_detections": len(registry.unresolved_crossings(workbench.resolution_ratio_threshold)),
        "points": int(sum(len(m.get_point_array()) for m in session.manifolds.values())) + blast_points,
        "pip_depth": pip_depth,
        "seconds": seconds,
        "symbols": names,
        "matrix": refined.tolist(),
        "letters": letters,
        "matrix_class": coarse.tolist(),
        "eigenvalues": [[float(z.real), float(z.imag)] for z in eigenvalues],
        "components": components(dynamics, names, refined),
    }


def data_path(case: str, area_cutoff: float) -> Path:
    """``data/<case>_a<cutoff>.json``."""
    return DATA_DIR / f"{case}_a{area_cutoff:.0e}.json"


#: The columns of a sweep's summary CSV (``data/<case>_a<cutoff>.csv``, its last step).
SUMMARY_COLUMNS = ["case", "area_cutoff", "blasts", "entropy", "entropy_class", "reliable", "points",
                   "seconds", "peak_gb", "stopped"]


def write_summary_row(record: dict, path: Path) -> None:
    """One CSV row summing a sweep's record up at its last step (empty fields before the first)."""
    last = record["rows"][-1] if record["rows"] else {}
    row = {"case": record["case"], "area_cutoff": f"{record['area_cutoff']:g}", "stopped": record["stopped"] or ""}
    row.update({key: last.get(key, "") for key in SUMMARY_COLUMNS if key not in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS)
        writer.writeheader()
        writer.writerow(row)


def resident_gb(field: str = "VmRSS") -> float:
    """This process's resident memory in GB (Linux ``/proc``; 0 elsewhere); ``VmHWM`` is its peak."""
    try:
        for line in Path("/proc/self/status").read_text().splitlines():
            if line.startswith(f"{field}:"):
                return int(line.split()[1]) / 1e6
    except OSError:
        pass
    return 0.0


def arc_length(bridge) -> float:
    """The polyline length of a bridge."""
    points = bridge.get_point_array()
    return float(np.sum(np.hypot(*np.diff(points, axis=0).T))) if len(points) > 1 else 0.0


def sweep(case: str, area_cutoff: float, max_blasts: int, max_seconds: Optional[float],
          max_rss_gb: Optional[float] = None, workers: int = 1) -> dict:
    """
    Build one case and blast it step by step, measuring after every step.

    Args:
        case: The case.
        area_cutoff: The refinement area cutoff.
        max_blasts: The last blast count measured.
        max_seconds: Stop before a step once the sweep has run this long.
        max_rss_gb: Stop before a step whose peak resident memory, extrapolated
            by the last step's growth, would pass this many GB.
        workers: The session's ``workers`` (set once the case is built).

    Returns:
        The sweep's record (also written to :func:`data_path` after every
        step, with its summary row beside it as ``.csv``).
    """
    path = data_path(case, area_cutoff)
    path.parent.mkdir(parents=True, exist_ok=True)
    record = {"case": case, "k": K, "b": B, "area_cutoff": area_cutoff,
              "min_separation": MIN_SEPARATION, "rows": [], "stopped": None}
    start = time.perf_counter()

    def save() -> None:
        path.write_text(json.dumps(record, indent=1))
        write_summary_row(record, path.with_suffix(".csv"))

    try:
        session, fixed_points, zones = build_case(case, area_cutoff)
        session.workers = workers
        pip_depth: dict = {}
        blast_points = 0
        previous_peak = resident_gb("VmHWM")
        frontier: dict = {}
        for blasts in range(max_blasts + 1):
            if blasts:
                if max_seconds is not None and time.perf_counter() - start > max_seconds:
                    record["stopped"] = f"time limit ({max_seconds:g} s) before blast {blasts}"
                    break
                # A blast peaks far above the resting memory (about 10x on p1), and the
                # peak grows by the frontier's factor per blast: extrapolate the peak.
                peak = resident_gb("VmHWM")
                predicted = peak * peak / max(previous_peak, 1e-9)
                if max_rss_gb is not None and predicted > max_rss_gb:
                    record["stopped"] = (f"memory limit ({max_rss_gb:g} GB): blast {blasts} predicted "
                                         f"to peak at {predicted:.1f} GB")
                    break
                previous_peak = peak
                frontier = {}
                for zone in zones:
                    result = session.blast_zone(
                        zone, num_iterations=1, fixed_point=[zone.fixed_point],
                        min_separation=MIN_SEPARATION, diagnose=BLAST_DIAGNOSE,
                    )
                    blast_points += result.points_created
                    pip_depth.update({f"p{fp.period}": depth for fp, depth in result.pip_depth.items()})
                    kept = [bridge for step in result.steps for bridge in step.kept_interior]
                    frontier[f"p{zone.fixed_point.period}"] = {
                        "bridges": len(kept), "length": sum(arc_length(bridge) for bridge in kept),
                    }
            row = measure(session, fixed_points, blasts, time.perf_counter() - start, dict(pip_depth), blast_points)
            row["frontier"] = frontier
            row["rss_gb"] = resident_gb()
            row["peak_gb"] = resident_gb("VmHWM")
            record["rows"].append(row)
            save()
            print(
                f"[{case} a={area_cutoff:g}] blast {blasts:2d}: h = {row['entropy']:.6f} "
                f"(class {row['entropy_class']:.6f}), {row['active']} active / {row['classes']} classes, "
                f"{row['refined_symbols']} refined, reliable={row['reliable']}, "
                f"{row['points']} points, {row['seconds']:.0f} s",
                flush=True,
            )
            if STOP_WHEN_EXHAUSTED and blasts and not any(f["bridges"] for f in frontier.values()):
                record["stopped"] = f"frontier exhausted at blast {blasts}"
                break
    except Exception as exc:  # noqa: BLE001 -- a failing step ends the sweep, recorded
        record["stopped"] = f"{type(exc).__name__}: {exc}"
        logger.warning("%s a=%g stopped: %s\n%s", case, area_cutoff, exc, traceback.format_exc())
    save()
    return record


# ── figures ─────────────────────────────────────────────────────────────────


def load_records() -> dict[tuple[str, float], dict]:
    """Every sweep record on disk, keyed by ``(case, cutoff)``."""
    records = {}
    for path in sorted(DATA_DIR.glob("*.json")):
        record = json.loads(path.read_text())
        if record["rows"]:
            records[(record["case"], record["area_cutoff"])] = record
    return records


#: One colour per cutoff, coarsest first (the same in every row and figure).
CUTOFF_COLORS = {1e-7: "#e6a100", 1e-8: "#d1495b", 1e-9: "#2e86ab", 1e-10: "#6a4c93", 1e-11: "#3c9d5d",
                 1e-12: "#222222", 1e-13: "#a0522d", 1e-14: "#e377c2", 1e-15: "#17becf"}


def cutoff_styles(cutoffs: list[float]) -> dict[float, dict]:
    """One marker, colour and alpha per cutoff, coarsest faintest."""
    markers = "osD^v<>p*h"
    ordered = sorted(cutoffs, reverse=True)
    spare = plt.get_cmap("tab10")
    return {
        c: {"marker": markers[i % len(markers)], "alpha": 0.55 + 0.45 * (i + 1) / len(ordered),
            "color": CUTOFF_COLORS.get(c, spare(i % 10))}
        for i, c in enumerate(ordered)
    }


def _column(rows: list[dict], key: str) -> np.ndarray:
    return np.array([row[key] for row in rows], dtype=float)


def _plot_entropy_line(ax, rows, key, color, style, linestyle, label) -> None:
    """One entropy line, its markers filled where the step is reliable."""
    x, y = _column(rows, "blasts"), _column(rows, key)
    reliable = np.array([row["reliable"] for row in rows])
    ax.plot(x, y, color=color, linestyle=linestyle, alpha=style["alpha"], label=label)
    ax.scatter(x[reliable], y[reliable], color=color, marker=style["marker"], alpha=style["alpha"], zorder=3)
    ax.scatter(x[~reliable], y[~reliable], facecolors="white", edgecolors=color,
               marker=style["marker"], alpha=style["alpha"], zorder=3)


def _stop_note(ax, record: dict, color: str) -> None:
    """Mark where a sweep stopped on an error (a resource stop is no error)."""
    if record["stopped"] and not record["stopped"].startswith(("time limit", "memory limit", "killed", "frontier exhausted")):
        last = record["rows"][-1]
        ax.annotate("stopped", (last["blasts"], last["entropy"]), xytext=(6, -12),
                    textcoords="offset points", fontsize=8, color=color)


def _horseshoe(ax) -> None:
    ax.axhline(math.log(2), color="0.5", linewidth=0.8, linestyle=":")
    ax.text(0.01, math.log(2), " ln 2 (full horseshoe)", transform=ax.get_yaxis_transform(),
            va="bottom", fontsize=8, color="0.4")


def reliability_handles() -> list:
    """Legend entries for the marker fill and the line style."""
    return [
        Line2D([], [], color="0.3", marker="o", linestyle="", label="reliable"),
        Line2D([], [], color="0.3", marker="o", markerfacecolor="white", linestyle="", label="not reliable"),
        Line2D([], [], color="0.3", linestyle="-", label=r"refined $M_1$"),
        Line2D([], [], color="0.3", linestyle="--", label=r"class $M$"),
    ]


def figure_entropy(records: dict) -> Path:
    """h against blasts, one panel per case, one line per cutoff."""
    styles = cutoff_styles(sorted({c for _, c in records}))
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), sharey=True)
    for ax, case in zip(axes, CASES):
        for (name, cutoff), record in sorted(records.items()):
            if name != case:
                continue
            color = styles[cutoff]["color"]
            _plot_entropy_line(ax, record["rows"], "entropy", color, styles[cutoff], "-",
                               f"cutoff {cutoff:g}")
            _plot_entropy_line(ax, record["rows"], "entropy_class", color, styles[cutoff], "--", None)
            _stop_note(ax, record, color)
        _horseshoe(ax)
        ax.set_title(CASE_TITLES[case])
        ax.set_xlabel("blasts (maximum iterate)")
        ax.grid(alpha=0.3)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(handles=ax.get_legend_handles_labels()[0] + reliability_handles(), fontsize=8)
    axes[0].set_ylabel(r"$h = \ln \rho(M)$")
    fig.suptitle(f"Hénon k={K:g}, b={B}: entropy of the trellis's symbolic dynamics against blasts")
    fig.tight_layout()
    return save_panel(fig, OUT_DIR, "entropy")


def figure_entropy_vs_depth(records: dict) -> Path:
    """
    h against the pip depth ``m`` (the deepest registered orbit point of the
    pip, which cuts ``T(q0)``), every cutoff and blast count as one point.

    The cutoff and the blasts act on the entropy through the depth the pip
    orbit is resolved to: curves that differ against blasts collapse here.
    """
    styles = cutoff_styles(sorted({c for _, c in records}))
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), sharey=True)
    for ax, case in zip(axes, CASES):
        period = "p3" if case == "p3" else "p1"
        for (name, cutoff), record in sorted(records.items()):
            if name != case:
                continue
            rows = [row for row in record["rows"] if period in row["pip_depth"]]
            depth = np.array([row["pip_depth"][period] for row in rows], dtype=float)
            h = _column(rows, "entropy")
            ax.plot(depth, h, color=styles[cutoff]["color"], marker=styles[cutoff]["marker"], linestyle="",
                    alpha=styles[cutoff]["alpha"], label=f"cutoff {cutoff:g}")
        _horseshoe(ax)
        ax.set_title(CASE_TITLES[case])
        ax.set_xlabel(f"pip depth m of {period}")
        ax.grid(alpha=0.3)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=8)
    axes[0].set_ylabel(r"$h = \ln \rho(M_1)$")
    fig.suptitle(f"Hénon k={K:g}: entropy against the depth of the trellis T(q0)")
    fig.tight_layout()
    return save_panel(fig, OUT_DIR, "entropy_vs_depth")


def figure_convergence_vs_depth(records: dict) -> Path:
    """
    Where h settles against the pip depth, every cutoff one line, the y axis zoomed to the plateau.

    (a) p1 alone against its depth; (b) p3 alone; (c) nested against the p1
    depth; (d) nested against the p3 depth. Blasts that leave the depth
    unchanged collapse onto one point (the last h at that depth); the dotted
    lines are the other case's final h at its finest sweep that ran to completion.
    """
    styles = cutoff_styles(sorted({c for _, c in records}))
    complete = {key: record for key, record in records.items()
                if not record["stopped"] or record["stopped"].startswith("frontier exhausted")}
    best = {case: record["rows"][-1]["entropy"] for case, record in finest(complete).items()}
    panels = [("p1", "p1", "(a) period 1 alone", "nested"), ("p3", "p3", "(b) period 3 alone", None),
              ("nested", "p1", "(c) nested, against the period-1 depth", "p1"),
              ("nested", "p3", "(d) nested, against the period-3 depth", "p1")]
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    for ax, (case, period, title, other) in zip(axes.flat, panels):
        plateau = []
        for (name, cutoff), record in sorted(records.items()):
            if name != case:
                continue
            last = {row["pip_depth"][period]: row["entropy"] for row in record["rows"] if period in row["pip_depth"]}
            if not last:
                continue
            depth, h = zip(*sorted(last.items()))
            plateau += [v for v in h if v > 0]
            ax.plot(depth, h, color=styles[cutoff]["color"], marker=styles[cutoff]["marker"], markersize=5,
                    drawstyle="steps-post", alpha=styles[cutoff]["alpha"], label=f"cutoff {cutoff:g}")
        if other in best:
            ax.axhline(best[other], color=CASE_COLORS[other], linestyle=":", linewidth=1.2,
                       label=f"{other} finest final: {best[other]:.7f}")
        if plateau:
            top = max(plateau + ([best[other]] if other in best else []))
            near = [v for v in plateau + ([best[other]] if other in best else []) if v >= 0.75 * top]
            lo, hi = min(near), max(near)
            span = max(hi - lo, 1e-6)
            ax.set_ylim(lo - 0.25 * span, hi + 0.25 * span)
        ax.yaxis.set_major_formatter(plt.FormatStrFormatter("%.5f"))
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(f"pip depth m of {period}")
        ax.set_ylabel(r"$h = \ln \rho(M_1)$")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, loc="lower right")
    fig.suptitle(f"Hénon k={K:g}: where h settles against the depth of the trellis T(q0)")
    fig.tight_layout()
    return save_panel(fig, OUT_DIR, "convergence_vs_depth")


def finest(records: dict) -> dict[str, dict]:
    """Per case, the record at its finest cutoff."""
    best: dict[str, dict] = {}
    for (case, cutoff), record in records.items():
        if case not in best or cutoff < best[case]["area_cutoff"]:
            best[case] = record
    return best


def component_series(record: dict) -> dict[str, tuple[list, list]]:
    """Per tangle label (``p1``, ``p3``, ``p1-p3``, ``mixed (...)``), the largest component entropy by blast."""
    series: dict[str, tuple[list, list]] = {}
    for row in record["rows"]:
        best: dict[str, float] = {}
        for comp in row["components"]:
            tangles = comp["tangles"]
            label = tangles[0] if len(tangles) == 1 else "mixed (" + ", ".join(tangles) + ")"
            best[label] = max(best.get(label, 0.0), comp["entropy"])
        for label, h in best.items():
            series.setdefault(label, ([], []))
            series[label][0].append(row["blasts"])
            series[label][1].append(h)
    return series


def figure_overlay(records: dict) -> Path:
    """The three cases on one axes at each case's finest cutoff."""
    fig, ax = plt.subplots(figsize=(10, 6))
    style = {"marker": "o", "alpha": 1.0}
    for case, record in sorted(finest(records).items(), key=lambda item: CASES.index(item[0])):
        _plot_entropy_line(ax, record["rows"], "entropy", CASE_COLORS[case], style, "-",
                           f"{CASE_TITLES[case]} (cutoff {record['area_cutoff']:g})")
        _stop_note(ax, record, CASE_COLORS[case])
    _horseshoe(ax)
    ax.set_xlabel("blasts (maximum iterate)")
    ax.set_ylabel(r"$h = \ln \rho(M_1)$")
    ax.grid(alpha=0.3)
    ax.legend(handles=ax.get_legend_handles_labels()[0] + reliability_handles()[:2], fontsize=9)
    ax.set_title(f"Hénon k={K:g}, b={B}: entropy against blasts, the three cases")
    fig.tight_layout()
    return save_panel(fig, OUT_DIR, "entropy_overlay")


def figure_nested_components(records: dict) -> Optional[Path]:
    """The nested run's component entropies against the cases alone."""
    best = finest(records)
    if "nested" not in best:
        return None
    fig, ax = plt.subplots(figsize=(10, 6))
    for case in ("p1", "p3"):
        if case in best:
            rows = best[case]["rows"]
            ax.plot(_column(rows, "blasts"), _column(rows, "entropy"), color=CASE_COLORS[case],
                    linewidth=3, alpha=0.3, label=f"{case} alone (cutoff {best[case]['area_cutoff']:g})")
    palette = {"p1": CASE_COLORS["p1"], "p3": CASE_COLORS["p3"]}  # heteroclinic purple, mixed green, both dashed
    for label, (x, y) in sorted(component_series(best["nested"]).items()):
        color = palette.get(label, "#9467bd" if label == "p1-p3" else CASE_COLORS["nested"])
        ax.plot(x, y, marker="o", color=color, linestyle="-" if label in palette else "--",
                label=f"nested: component in {label}")
    rows = best["nested"]["rows"]
    ax.plot(_column(rows, "blasts"), _column(rows, "entropy"), color="black", linewidth=0.8,
            linestyle=":", label="nested: whole matrix")
    _horseshoe(ax)
    ax.set_xlabel("blasts of each zone")
    ax.set_ylabel("entropy of the strongly connected component")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    ax.set_title(f"Hénon k={K:g}: nested components against the tangles alone "
                 f"(cutoff {best['nested']['area_cutoff']:g})")
    fig.tight_layout()
    return save_panel(fig, OUT_DIR, "nested_components")


def _panel_grid(records: dict, panels: list[tuple[str, str, bool]], name: str, title: str) -> Path:
    """One row per case, one column per quantity, one line per cutoff."""
    styles = cutoff_styles(sorted({c for _, c in records}))
    fig, axes = plt.subplots(len(CASES), len(panels), figsize=(4.2 * len(panels), 3.4 * len(CASES)),
                             squeeze=False, sharex="row")
    for i, case in enumerate(CASES):
        for j, (key, label, log) in enumerate(panels):
            ax = axes[i, j]
            for (name_, cutoff), record in sorted(records.items()):
                if name_ != case:
                    continue
                rows = record["rows"]
                x = _column(rows, "blasts")
                if key == "pip_depth":
                    for period, linestyle in (("p1", "-"), ("p3", "--")):
                        y = [row["pip_depth"].get(period, np.nan) for row in rows]
                        if not np.all(np.isnan(np.array(y, dtype=float))):
                            ax.plot(x, y, color=styles[cutoff]["color"], linestyle=linestyle,
                                    marker=styles[cutoff]["marker"], alpha=styles[cutoff]["alpha"],
                                    label=f"{period} {cutoff:g}")
                    continue
                y = (np.abs(_column(rows, "radius") - _column(rows, "radius_class"))
                     if key == "radius_gap" else _column(rows, key))
                ax.plot(x, y, color=styles[cutoff]["color"], marker=styles[cutoff]["marker"],
                        alpha=styles[cutoff]["alpha"], label=f"{cutoff:g}")
            if log and key == "radius_gap":
                ax.set_yscale("symlog", linthresh=1e-10)  # below it: floating-point noise
            elif log:
                ax.set_yscale("symlog" if key == "unresolved_detections" else "log")
            ax.grid(alpha=0.3)
            if i == 0:
                ax.set_title(label, fontsize=10)
            if j == 0:
                ax.set_ylabel(CASE_TITLES[case], fontsize=9)
            if i == len(CASES) - 1:
                ax.set_xlabel("blasts")
            if ax.get_legend_handles_labels()[0] and j == len(panels) - 1:
                ax.legend(fontsize=7, title="cutoff", title_fontsize=7)
    fig.suptitle(title)
    fig.tight_layout()
    return save_panel(fig, OUT_DIR, name)


def figure_alphabet(records: dict) -> Path:
    """Alphabet and trellis sizes and the pip depth against blasts."""
    return _panel_grid(records, [
        ("active", "active classes", False),
        ("inert", "inert classes", False),
        ("refined_symbols", r"refined symbols (size of $M_1$)", False),
        ("minimal_bridges", "minimal-trellis bridges", False),
        ("trellis_bridges", "trellis bridges", True),
        ("pip_depth", "pip depth m (solid p1, dashed p3)", False),
    ], "alphabet", f"Hénon k={K:g}: alphabet, trellis and pip depth against blasts")


def figure_health(records: dict) -> Path:
    """The readings that say how far a step's matrix can be trusted."""
    return _panel_grid(records, [
        ("unresolved_classes", "unresolved classes (no itinerary)", False),
        ("ambiguous_classes", "ambiguous classes", False),
        ("unspanned_classes", "classes no computed bridge spans", False),
        ("unresolved_detections", "unresolved detections (not in N)", True),
        ("radius_gap", r"$|\rho(M_1) - \rho(M)|$ (in-splitting gap)", True),
    ], "health", f"Hénon k={K:g}: health of each step's symbolic dynamics")


def figure_cost(records: dict) -> Path:
    """Points held and wall time against blasts."""
    return _panel_grid(records, [
        ("points", "manifold points", True),
        ("crossings", "crossings of N", True),
        ("seconds", "cumulative wall time (s)", True),
    ], "cost", f"Hénon k={K:g}: cost of the sweep")


def figure_convergence(records: dict) -> Optional[Path]:
    """h at each blast count against the cutoff, one panel per case."""
    cutoffs = sorted({c for _, c in records})
    if len(cutoffs) < 2:
        return None
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), sharey=True)
    cmap = plt.get_cmap("viridis")
    for ax, case in zip(axes, CASES):
        by_blast: dict[int, list[tuple[float, float, bool]]] = {}
        for (name, cutoff), record in records.items():
            if name == case:
                for row in record["rows"]:
                    by_blast.setdefault(row["blasts"], []).append((cutoff, row["entropy"], row["reliable"]))
        top = max(by_blast, default=1)
        for blasts, points in sorted(by_blast.items()):
            points.sort()
            x, y, reliable = (np.array(v) for v in zip(*points))
            color = cmap(blasts / max(top, 1))
            ax.plot(x, y, color=color, marker="o", markersize=0, label=f"{blasts}")
            ax.scatter(x[reliable], y[reliable], color=color, zorder=3)
            ax.scatter(x[~reliable], y[~reliable], facecolors="white", edgecolors=color, zorder=3)
        ax.set_xscale("log")
        ax.invert_xaxis()
        _horseshoe(ax)
        ax.set_title(CASE_TITLES[case])
        ax.set_xlabel("area cutoff (finer to the right)")
        ax.grid(alpha=0.3)
        if by_blast:
            mappable = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, max(top, 1)))
            fig.colorbar(mappable, ax=ax, label="blasts")
    axes[0].set_ylabel(r"$h = \ln \rho(M_1)$")
    fig.suptitle(f"Hénon k={K:g}: does the entropy at a blast count depend on the cutoff?")
    fig.tight_layout()
    return save_panel(fig, OUT_DIR, "convergence")


def figure_length_growth(records: dict) -> Optional[Path]:
    """
    The geometric cross-check: the growth of the blast frontier's length.

    The interior frontier of a blast is the unstable curve the zone still holds
    after ``n`` iterations; its length grows like ``e^{h n}`` (Newhouse:
    length growth bounds the entropy from below, and is sharp for the curve's
    eventual growth). Top: ``ln L_n``; bottom: ``ln(L_n / L_{n-1})`` beside the
    symbolic ``h``.
    """
    best = finest(records)
    cases = [case for case in CASES if case in best and any("frontier" in row for row in best[case]["rows"])]
    if not cases:
        return None
    fig, axes = plt.subplots(2, len(cases), figsize=(6 * len(cases), 9), squeeze=False, sharex="col")
    for j, case in enumerate(cases):
        record = best[case]
        rows = [row for row in record["rows"] if row.get("frontier")]
        tangles = sorted({t for row in rows for t in row["frontier"]})
        for tangle, linestyle in zip(tangles, ("-", "--")):
            x = np.array([row["blasts"] for row in rows if tangle in row["frontier"]], dtype=float)
            length = np.array([row["frontier"][tangle]["length"] for row in rows if tangle in row["frontier"]])
            good = length > 0
            axes[0, j].plot(x[good], np.log(length[good]), color=CASE_COLORS[case], linestyle=linestyle,
                            marker="o", label=f"frontier of zone {tangle}")
            rate = np.log(length[1:] / length[:-1], where=(length[1:] > 0) & (length[:-1] > 0),
                          out=np.full(len(length) - 1, np.nan))
            axes[1, j].plot(x[1:], rate, color=CASE_COLORS[case], linestyle=linestyle, marker="o",
                            label=rf"$\ln(L_n/L_{{n-1}})$, zone {tangle}")
        h = _column(record["rows"], "entropy")
        axes[1, j].plot(_column(record["rows"], "blasts"), h, color="black", linewidth=1.2,
                        label=r"symbolic $h = \ln\rho(M_1)$")
        _horseshoe(axes[1, j])
        axes[0, j].set_title(f"{CASE_TITLES[case]} (cutoff {record['area_cutoff']:g})")
        axes[0, j].set_ylabel(r"$\ln L_n$ (frontier arc length)")
        axes[1, j].set_ylabel("growth rate per blast")
        axes[1, j].set_xlabel("blasts")
        for ax in axes[:, j]:
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8)
        axes[1, j].set_ylim(-0.5, 1.5)
    fig.suptitle(f"Hénon k={K:g}: length growth of the blasted unstable curve inside each zone")
    fig.tight_layout()
    return save_panel(fig, OUT_DIR, "length_growth")


def figure_matrices(records: dict) -> Path:
    """The last step's refined transition matrix per case (finest cutoff)."""
    best = finest(records)
    cases = [case for case in CASES if case in best]
    fig, axes = plt.subplots(1, len(cases), figsize=(7 * len(cases), 7), squeeze=False)
    for ax, case in zip(axes[0], cases):
        row = best[case]["rows"][-1]
        matrix = np.array(row["matrix"])
        names = row["symbols"]
        if matrix.size:
            if matrix.max() > 1:
                image = ax.imshow(np.where(matrix > 0, matrix, np.nan), cmap="viridis", interpolation="nearest")
                fig.colorbar(image, ax=ax, shrink=0.7, label="tokens")
            else:
                ax.imshow(matrix > 0, cmap="Greys", interpolation="nearest", vmin=0, vmax=1.6)
            ax.set_xticks(range(len(names)), [f"${n}$" for n in names], rotation=90, fontsize=6)
            ax.set_yticks(range(len(names)), [f"${n}$" for n in names], fontsize=6)
            if len(names) <= 20:
                for (r, c), value in np.ndenumerate(matrix):
                    if value:
                        ax.text(c, r, str(value), ha="center", va="center", fontsize=7,
                                color="white" if matrix.max() > 1 else "black")
        ax.set_title(f"{CASE_TITLES[case]}\n{row['blasts']} blasts, cutoff {best[case]['area_cutoff']:g}: "
                     f"{len(names)} symbols, h = {row['entropy']:.4f}", fontsize=10)
        ax.set_xlabel("symbol in the word")
        ax.set_ylabel("symbol")
    fig.suptitle(r"refined transition matrix $M_1$ at the last step (row $i$ counts the tokens of $j$ "
                 r"in $f(i)$)")
    fig.tight_layout()
    return save_panel(fig, OUT_DIR, "matrices")


def figure_spectra(records: dict) -> Path:
    """The last step's eigenvalues per case in the complex plane."""
    best = finest(records)
    fig, ax = plt.subplots(figsize=(8, 8))
    theta = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(theta), np.sin(theta), color="0.6", linewidth=0.8, linestyle=":")
    for case in CASES:
        if case not in best:
            continue
        row = best[case]["rows"][-1]
        z = np.array(row["eigenvalues"]).reshape(-1, 2)
        ax.scatter(z[:, 0], z[:, 1], color=CASE_COLORS[case], alpha=0.75,
                   label=f"{CASE_TITLES[case]} ({row['blasts']} blasts): rho = {row['radius']:.4f}")
        ax.plot(row["radius"] * np.cos(theta), row["radius"] * np.sin(theta),
                color=CASE_COLORS[case], linewidth=0.8, alpha=0.6)
    ax.set_aspect("equal")
    ax.axhline(0, color="0.8", linewidth=0.5)
    ax.axvline(0, color="0.8", linewidth=0.5)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    ax.set_title(r"eigenvalues of $M_1$ at the last step (circles: spectral radius; dotted: unit circle)")
    fig.tight_layout()
    return save_panel(fig, OUT_DIR, "spectra")


def write_summary(records: dict) -> Path:
    """A plain-text table of every row, beside the figures."""
    lines = []
    for (case, cutoff), record in sorted(records.items()):
        lines.append(f"{case}  cutoff {cutoff:g}  stopped: {record['stopped'] or '-'}")
        lines.append("  blasts  h(M1)     h(M)      reliable  active inert refined minimal trellis  "
                     "points     seconds")
        for row in record["rows"]:
            lines.append(
                f"  {row['blasts']:6d}  {row['entropy']:.6f}  {row['entropy_class']:.6f}  "
                f"{row['reliable']!s:8}  {row['active']:6d} {row['inert']:5d} {row['refined_symbols']:7d} "
                f"{row['minimal_bridges']:7d} {row['trellis_bridges']:7d}  {row['points']:9d}  "
                f"{row['seconds']:7.0f}"
            )
        lines.append("")
    path = OUT_DIR / "summary.txt"
    path.write_text("\n".join(lines))
    return path


def draw_figures() -> None:
    """Every figure from the data on disk."""
    records = load_records()
    if not records:
        print(f"no data in {DATA_DIR}")
        return
    for draw in (figure_entropy, figure_entropy_vs_depth, figure_convergence_vs_depth, figure_overlay, figure_nested_components, figure_alphabet,
                 figure_health, figure_cost, figure_convergence, figure_length_growth, figure_matrices, figure_spectra,
                 write_summary):
        path = draw(records)
        if path is not None:
            print(f"wrote {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Entropy against blasts for the Hénon k=2.5 tangles.")
    parser.add_argument("--cases", nargs="+", choices=CASES, default=list(CASES))
    parser.add_argument("--cutoffs", nargs="+", type=float, default=[1e-7, 1e-8],
                        help="refinement area cutoffs, one sweep per case and cutoff")
    parser.add_argument("--max-blasts", type=int, default=12, help="the last blast count measured")
    parser.add_argument("--max-seconds", type=float, default=None,
                        help="stop a sweep before the next blast once it has run this long")
    parser.add_argument("--max-rss-gb", type=float, default=20.0,
                        help="stop a sweep before a blast predicted to peak past this resident memory")
    parser.add_argument("--workers", type=int, default=1, help="the cores of each sweep's session (default 1)")
    parser.add_argument("--plot-only", action="store_true", help="redraw the figures from the data")
    args = parser.parse_args()
    configure_logging()
    plt.switch_backend("Agg")

    if not args.plot_only:
        for cutoff in args.cutoffs:
            for case in args.cases:
                sweep(case, cutoff, args.max_blasts, args.max_seconds, args.max_rss_gb, args.workers)
    draw_figures()


if __name__ == "__main__":
    main()
