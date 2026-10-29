"""
The computed tangle of each k=2.5 entropy case after 1..5 single-iteration blasts.

Each case (``p1``, ``p3``, ``nested``) is built and blasted exactly as the
entropy sweep does (``henon_k25_entropy.build_case``, one ``num_iterations=1``
blast of every zone per step, inner zone first, ``min_separation`` 1e-5), and
after each step one panel shows what the topology was computed from:

* the resonance zones, shaded;
* the stable manifolds (red) and the unstable manifolds as grown (light grey);
* every bridge the workbench holds (light blue: blast children included);
* the bridges of the selection's trellis ``T(q0)``, the trellis cut at the pip
  orbit that the symbolic dynamics reads (dark blue);
* the strong pip(s) (black circles) and the fixed points.

The panel title gives the step's entropy ``h = ln rho(M_1)``, the refined
alphabet size, the pip depth and the bridge counts. One view per case (the
outermost zone, padded), so the panels of a row compare directly.

Writes ``figures/henon_k25_entropy/blast_panels/<case>.{png,svg,fig.pickle}``
(one 1x5 row per case, large) and ``overview`` (the three rows together).

Run with ``PYTHONPATH=src python3 scripts/henon_k25_blast_panels.py``
(``--cutoff 1e-8``, ``--blasts 5``, ``--cases p1 p3 nested``).
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parent))

import henon_k25_entropy as entropy  # noqa: E402
from symbolic_figures import configure_logging, save_panel  # noqa: E402

from tanglepack import TangleSession  # noqa: E402
from tanglepack.loom import plot_resonance_zones  # noqa: E402

logger = logging.getLogger(__name__)

OUT_DIR = entropy.OUT_DIR / "blast_panels"
#: Fraction of the outermost zone's extent added on every side of the view.
VIEW_PAD = 0.15


def zone_view(session: TangleSession, zones: list) -> tuple[tuple[float, float], tuple[float, float]]:
    """The padded bounding box of the outermost zone's boundary."""
    outer = max(zones, key=lambda zone: zone.area)
    vertices = outer.boundary_vertices
    if vertices is None:
        vertices = outer.boundary_polygon(session.workbench)
    lo, hi = vertices.min(axis=0), vertices.max(axis=0)
    pad = VIEW_PAD * (hi - lo)
    return (lo[0] - pad[0], hi[0] + pad[0]), (lo[1] - pad[1], hi[1] + pad[1])


def snapshot(session: TangleSession, fixed_points: list, zones: list) -> dict:
    """
    The step's drawable state and headline numbers, read once after a blast.

    Returns:
        The polylines (stable, unstable, every bridge, trellis bridges), the
        pips and fixed points, and ``h``, the refined symbol count and the
        bridge counts.
    """
    topology = session.topology()
    names, refined = topology.dynamics.transition_matrix(refined=True)
    trellis = topology.trellis
    manifolds = session.manifolds
    return {
        "stable": [m.get_point_array() for key, m in manifolds.items() if key[1] == "stable"],
        "unstable": [m.get_point_array() for key, m in manifolds.items() if key[1] == "unstable"],
        "bridges": [b.get_point_array() for b in session.workbench.bridges],
        "trellis": [b.get_point_array() for b in trellis.bridges],
        "pip_points": _pip_points(session, topology),
        "fixed": np.vstack([np.asarray(fp.coordinates).reshape(-1, 2) for fp in fixed_points]),
        "entropy": entropy.entropy_of(entropy.spectral_radius(refined)),
        "symbols": len(names),
        "reliable": bool(topology.dynamics.is_reliable),
        "n_bridges": len(session.workbench.bridges),
        "n_trellis": len(trellis.bridges),
    }


def _pip_points(session: TangleSession, topology) -> np.ndarray:
    """The plane points of the chosen strong pips (registry ids; an empty array when none)."""
    registry = session.workbench.intersection_registry
    points = [registry[pip].coords for pip in topology.strong_pips if pip is not None]
    return np.asarray(points, dtype=float).reshape(-1, 2)


def draw_panel(ax, state: dict, zones: list, view: tuple, blasts: int, depth: dict) -> None:
    """One step's panel: zones, manifolds, every bridge, the trellis bridges, pips."""
    plot_resonance_zones(zones, ax=ax, alpha=0.18, zorder=0)
    for points in state["unstable"]:
        ax.plot(points[:, 0], points[:, 1], color="0.8", linewidth=0.5, zorder=1)
    for points in state["stable"]:
        ax.plot(points[:, 0], points[:, 1], color="tab:red", linewidth=0.8, zorder=2)
    for points in state["bridges"]:
        ax.plot(points[:, 0], points[:, 1], color="#9ecae1", linewidth=0.6, zorder=3)
    for points in state["trellis"]:
        ax.plot(points[:, 0], points[:, 1], color="#08519c", linewidth=0.8, zorder=4)
    if len(state["pip_points"]):
        ax.scatter(*state["pip_points"].T, s=60, facecolors="none", edgecolors="k", linewidths=1.2, zorder=6)
    ax.scatter(*state["fixed"].T, s=14, color="k", zorder=6)
    ax.set_xlim(*view[0])
    ax.set_ylim(*view[1])
    ax.set_aspect("equal")
    ax.tick_params(labelsize=7)
    depth_text = ", ".join(f"{period} m={m}" for period, m in sorted(depth.items())) or "-"
    ax.set_title(
        f"{blasts} blast{'s' if blasts != 1 else ''}: h = {state['entropy']:.4f}"
        f"{'' if state['reliable'] else ' (not reliable)'}\n"
        f"{state['symbols']} refined symbols, depth {depth_text}\n"
        f"{state['n_trellis']} trellis / {state['n_bridges']} computed bridges",
        fontsize=9,
    )


def legend_handles() -> list:
    """The legend shared by every row."""
    return [
        Patch(facecolor="#8c6bb1", alpha=0.3, label="resonance zone(s)"),
        Line2D([], [], color="tab:red", label="stable manifold"),
        Line2D([], [], color="0.8", label="unstable manifold (grown)"),
        Line2D([], [], color="#9ecae1", label="computed bridge (blast children included)"),
        Line2D([], [], color="#08519c", label=r"bridge of the trellis $T(q_0)$ (read by the topology)"),
        Line2D([], [], marker="o", linestyle="", markerfacecolor="none", markeredgecolor="k",
               label="strong pip"),
    ]


def run_case(case: str, cutoff: float, blasts: int) -> tuple[list, list, tuple]:
    """Build one case, blast it ``blasts`` times, snapshot after each step."""
    session, fixed_points, zones = entropy.build_case(case, cutoff)
    view = zone_view(session, zones)
    states, depths = [], []
    pip_depth: dict = {}
    for step in range(1, blasts + 1):
        for zone in zones:
            result = session.blast_zone(zone, num_iterations=1, fixed_point=[zone.fixed_point],
                                        min_separation=entropy.MIN_SEPARATION)
            pip_depth.update({f"p{fp.period}": depth for fp, depth in result.pip_depth.items()})
        states.append(snapshot(session, fixed_points, zones))
        depths.append(dict(pip_depth))
        logger.info("%s blast %d: h=%.6f, %d symbols", case, step, states[-1]["entropy"], states[-1]["symbols"])
    return states, depths, (zones, view)


def main() -> None:
    parser = argparse.ArgumentParser(description="The k=2.5 tangles after 1..N blasts, per case.")
    parser.add_argument("--cases", nargs="+", choices=entropy.CASES, default=list(entropy.CASES))
    parser.add_argument("--cutoff", type=float, default=1e-8)
    parser.add_argument("--blasts", type=int, default=5)
    args = parser.parse_args()
    configure_logging()
    plt.switch_backend("Agg")

    results = {case: run_case(case, args.cutoff, args.blasts) for case in args.cases}
    overview, grid = plt.subplots(len(results), args.blasts, figsize=(4.4 * args.blasts, 4.6 * len(results)),
                                  squeeze=False)
    for row, (case, (states, depths, (zones, view))) in enumerate(results.items()):
        fig, axes = plt.subplots(1, args.blasts, figsize=(7 * args.blasts, 7.6), squeeze=False)
        for col, (state, depth) in enumerate(zip(states, depths)):
            draw_panel(axes[0, col], state, zones, view, col + 1, depth)
            draw_panel(grid[row, col], state, zones, view, col + 1, depth)
        grid[row, 0].set_ylabel(entropy.CASE_TITLES[case], fontsize=10)
        fig.legend(handles=legend_handles(), loc="lower center", ncol=6, fontsize=10)
        fig.suptitle(f"Hénon k={entropy.K:g}: {entropy.CASE_TITLES[case]}, cutoff {args.cutoff:g}, "
                     f"blasts 1..{args.blasts}")
        fig.tight_layout(rect=(0, 0.05, 1, 1))
        save_panel(fig, OUT_DIR, case)
    overview.legend(handles=legend_handles(), loc="lower center", ncol=6, fontsize=9)
    overview.suptitle(f"Hénon k={entropy.K:g}: the computed tangle after 1..{args.blasts} blasts "
                      f"(cutoff {args.cutoff:g})")
    overview.tight_layout(rect=(0, 0.03, 1, 1))
    save_panel(overview, OUT_DIR, "overview")


if __name__ == "__main__":
    main()
