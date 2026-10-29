"""
Where the independent curve and the unpruned library blast differ, piece by piece.

``henon_k25_length_growth.py`` measures the unstable length kept in a zone two
ways, and they agree on the growth rate but not on the length (the curve holds
2.5x the blast's length in p1, 6.6x in p3). This script runs both to one step
``n`` of one case and sorts every segment of the curve's step-``n`` image by
where the blast has it:

* ``matched`` -- within ``MATCH_TOL`` of a bridge the blast kept at step ``n``
  (``kept_interior``, the length the black boxes plot);
* ``computed, not counted`` -- not one of those, but on unstable manifold the
  workbench holds (an earlier step's bridge, a child the blast found already
  known, the grown manifold);
* ``never computed`` -- on no curve the workbench holds;

and marks the blast's step-``n`` bridges the curve does not hold
(``blast only``). The lengths of each group are printed in the panel.

Writes ``figures/henon_k25_entropy/length_growth/compare_<case>_step<n>.{png,svg,fig.pickle}``.

Run with ``PYTHONPATH=src env/bin/python scripts/henon_k25_length_compare.py``
(``--cases p1 p3``, ``--steps 8 20`` one per case, ``--cutoff 1e-7``, ``--tol 1e-5``).
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.path import Path as PolygonPath
from numpy.typing import NDArray
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))

import henon_k25_entropy as entropy  # noqa: E402
import henon_k25_length_growth as growth  # noqa: E402
from symbolic_figures import configure_logging, save_panel  # noqa: E402

from tanglepack.examples.henon import henon_map  # noqa: E402
from tanglepack.loom import plot_resonance_zones  # noqa: E402

logger = logging.getLogger(__name__)

#: A curve segment within this distance of a blast polyline is the same manifold
#: (the curve's chord tolerance plus the blast's chord deviation at 1e-7, with margin).
MATCH_TOL = 5e-5
#: Nearest segment midpoints tested per query point.
NEIGHBOURS = 16

GROUPS = {
    "matched": ("0.65", "matched: a bridge the blast kept at this step"),
    "computed": ("#e6a100", "computed, not counted at this step"),
    "never": ("#d1495b", "never computed by the blast"),
    "blast_only": ("#2e86ab", "blast only: kept by the blast, not by the curve"),
}


class Segments:
    """The segments of a set of polylines, indexed by midpoint for nearest-distance queries."""

    def __init__(self, polylines: list[NDArray[np.float64]]) -> None:
        polylines = [p for p in polylines if len(p) > 1]
        self.a = np.vstack([p[:-1] for p in polylines]) if polylines else np.zeros((0, 2))
        self.b = np.vstack([p[1:] for p in polylines]) if polylines else np.zeros((0, 2))
        self.tree = cKDTree(0.5 * (self.a + self.b)) if len(self.a) else None

    def distance(self, probe: NDArray[np.float64]) -> NDArray[np.float64]:
        """Distance from each probe to the nearest of the ``NEIGHBOURS`` segments with the nearest midpoints."""
        if self.tree is None:
            return np.full(len(probe), np.inf)
        k = min(NEIGHBOURS, len(self.a))
        _, idx = self.tree.query(probe, k=k)
        idx = idx.reshape(len(probe), k)
        a, b = self.a[idx], self.b[idx]
        ab = b - a
        t = np.clip(((probe[:, None, :] - a) * ab).sum(axis=2) / np.maximum((ab * ab).sum(axis=2), 1e-300), 0, 1)
        return np.linalg.norm(probe[:, None, :] - (a + t[..., None] * ab), axis=2).min(axis=1)


def curve_segments(curve: growth.Curve) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """The curve's segments within one piece, as ``(starts, ends)``."""
    same = curve.ids[1:] == curve.ids[:-1]
    return curve.points[:-1][same], curve.points[1:][same]


def run_case(case: str, step: int, cutoff: float, tol: float) -> dict:
    """Both methods to ``step``; the curve's segments grouped, the blast's unmatched bridges."""
    session, _, zones = entropy.build_case(case, cutoff)
    zone = max(zones, key=lambda z: z.area) if case == "nested" else zones[0]
    inner = [z for z in session.resonance_zones.values() if z.key is not zone.key and z.area < zone.area]
    keep = PolygonPath(zone.boundary_vertices)
    exclude = [PolygonPath(z.boundary_vertices) for z in inner]
    graphs = [g for z in [zone, *inner] for arc in growth.stable_arcs(session, z) for g in growth.monotone_pieces(arc)]
    on_boundary = growth.unstable_arcs(session, zone)

    initial, known_length = None, 0.0
    for n in range(1, step + 1):
        result = session.blast_zone(zone, num_iterations=1, fixed_point=[zone.fixed_point], min_separation=None,
                                    diagnose=False)
        last = result.steps[-1]
        if initial is None:
            initial = [b.get_point_array().copy() for b in last.interior_parents]
        known_length = float(sum(entropy.arc_length(b) for b in last.already_known))
        logger.info("%s blast %d: %d kept, %d already known", case, n, len(last.kept_interior), len(last.already_known))
    kept = [b.get_point_array() for b in last.kept_interior]
    workbench = session.workbench
    held = [m.get_point_array() for key, m in session.manifolds.items() if key[1] == "unstable"]
    held += [b.get_point_array() for b in workbench.bridges]
    held += [np.asarray(c.get_point_array()) for c in getattr(workbench, "_image_curves", [])]

    curve = growth.Curve.from_polylines(initial)
    fmap = henon_map(entropy.K, entropy.B)
    for _ in range(step):
        curve, _ = growth.iterate_curve(fmap, curve, tol, graphs, keep, exclude, on_boundary)
    a, b = curve_segments(curve)
    mid = 0.5 * (a + b)
    lengths = np.hypot(*(b - a).T)
    near_kept = Segments(kept).distance(mid) < MATCH_TOL
    near_held = Segments(held).distance(mid) < MATCH_TOL
    group = np.where(near_kept, "matched", np.where(near_held, "computed", "never"))

    ka = np.vstack([p[:-1] for p in kept if len(p) > 1])
    kb = np.vstack([p[1:] for p in kept if len(p) > 1])
    curve_index = Segments([np.array([s, e]) for s, e in zip(a, b)]) if len(a) else Segments([])
    blast_only = curve_index.distance(0.5 * (ka + kb)) >= MATCH_TOL
    kept_lengths = np.hypot(*(kb - ka).T)
    return {
        "session": session, "zones": zones, "step": step,
        "segments": (a, b), "group": group, "lengths": lengths,
        "blast_only": (ka[blast_only], kb[blast_only]),
        "totals": {
            "curve": float(lengths.sum()),
            "blast": float(kept_lengths.sum()),
            "matched": float(lengths[group == "matched"].sum()),
            "computed": float(lengths[group == "computed"].sum()),
            "never": float(lengths[group == "never"].sum()),
            "blast_only": float(kept_lengths[blast_only].sum()),
            "already_known": known_length,
        },
    }


def _draw_segments(ax, a: NDArray[np.float64], b: NDArray[np.float64], color: str, width: float, z: int) -> None:
    """Segments as one NaN-separated line."""
    if not len(a):
        return
    xy = np.full((3 * len(a), 2), np.nan)
    xy[0::3], xy[1::3] = a, b
    ax.plot(xy[:, 0], xy[:, 1], color=color, linewidth=width, zorder=z)


def figure(case: str, data: dict) -> None:
    """The zone with the curve's segments coloured by group and the blast-only bridges."""
    fig, ax = plt.subplots(figsize=(11, 10))
    plot_resonance_zones(data["zones"], ax=ax, alpha=0.12, zorder=0)
    for key, m in data["session"].manifolds.items():
        if key[1] == "stable":
            p = m.get_point_array()
            ax.plot(p[:, 0], p[:, 1], color="0.85", linewidth=0.4, zorder=1)
    a, b = data["segments"]
    for name, width, z in (("matched", 0.6, 2), ("computed", 1.0, 3), ("never", 1.0, 4)):
        mask = data["group"] == name
        _draw_segments(ax, a[mask], b[mask], GROUPS[name][0], width, z)
    _draw_segments(ax, *data["blast_only"], GROUPS["blast_only"][0], 1.4, 5)
    view = data["zones"][-1].boundary_vertices
    lo, hi = view.min(axis=0), view.max(axis=0)
    pad = 0.08 * (hi - lo)
    ax.set_xlim(lo[0] - pad[0], hi[0] + pad[0])
    ax.set_ylim(lo[1] - pad[1], hi[1] + pad[1])
    ax.set_aspect("equal")
    t = data["totals"]
    rows = [
        f"curve length {t['curve']:.1f}   blast step-{data['step']} length {t['blast']:.1f}",
        f"matched {t['matched']:.1f} ({t['matched'] / t['curve']:.0%})",
        f"computed, not counted {t['computed']:.1f} ({t['computed'] / t['curve']:.0%})",
        f"never computed {t['never']:.1f} ({t['never'] / t['curve']:.0%})",
        f"blast only {t['blast_only']:.2f}   (step's already-known children: {t['already_known']:.1f})",
    ]
    ax.text(0.01, 0.01, "\n".join(rows), transform=ax.transAxes, fontsize=9, va="bottom",
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="0.7"))
    handles = [Line2D([], [], color=color, linewidth=2, label=label) for color, label in GROUPS.values()]
    ax.legend(handles=handles, loc="upper right", fontsize=9)
    ax.set_title(f"Hénon k={entropy.K:g}, {entropy.CASE_TITLES[case]}: the curve's step-{data['step']} image "
                 f"against the unpruned blast's", fontsize=11)
    fig.tight_layout()
    save_panel(fig, growth.OUT_DIR, f"compare_{case}_step{data['step']}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Where the independent curve and the library blast differ.")
    parser.add_argument("--cases", nargs="+", choices=entropy.CASES, default=["p1", "p3"])
    parser.add_argument("--steps", nargs="+", type=int, default=[8, 20], help="one step per case")
    parser.add_argument("--cutoff", type=float, default=1e-7)
    parser.add_argument("--tol", type=float, default=1e-5)
    args = parser.parse_args()
    configure_logging()
    plt.switch_backend("Agg")
    for case, step in zip(args.cases, args.steps):
        data = run_case(case, step, args.cutoff, args.tol)
        logger.info("%s step %d: %s", case, step, {k: round(v, 3) for k, v in data["totals"].items()})
        figure(case, data)


if __name__ == "__main__":
    main()
