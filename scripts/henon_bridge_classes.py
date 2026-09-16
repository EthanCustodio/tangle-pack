"""
Bridge classes on two Hénon tangles: k=10 and k=2.8.

One row per case, three panels each:

1. the tangle — stable manifold (red), unstable manifold and every bridge
   (blue), crossings with their registry ids, the chosen strong pip (green ring);
2. the tangle again with ONE bridge per class: the member with the smallest
   unstable canonical distance, coloured per class and labelled with the symbol
   it contributes (``a`` or ``a^-1``; an inert class has no letter and is
   drawn dashed);
3. the stable partition as number lines.

Each case's class table is printed to stdout. Both cases trim the stable
manifold at the chosen strong pip (a resonance zone) so the partition sees
that cut, and classes are computed once, after the final partition, so the
letters start at ``a``.

Run:  env/bin/python scripts/henon_bridge_classes.py

Dev Notes:
    The k=10 recipe is ``scripts/henon_k10_dual_graph.py`` (saddle at
    ``[4, -4]``, ten unstable steps, pip = f(q0)); the k=2.8 recipe is
    ``scripts/henon_pseudoneighbors.py`` (same saddle guess, ten steps,
    ``area_cutoff = 1e-7``, pip = f(q0), one blast of the zone). Registry ids
    are not reproducible between builds, so the pip is named by its relation to
    the default candidate (``"default"`` = q0, ``"image"`` = f(q0)), never by id.
    With ``PIP = "default"`` and one blast the k=2.8 case shows a folded loop
    (an inert class with three members) instead of four single-bridge classes.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Optional

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from tanglepack import TangleSession
from tanglepack.examples import henon_jacobian, henon_map, henon_map_inverse
from tanglepack.numerics.geometry import polyline_midpoint
from tanglepack.topology.BridgeClass import BridgeClassEntry, BridgeClassTable

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"
OUT_PATH = FIGURES_DIR / "henon_bridge_classes.png"

SADDLE_GUESS = [4, -4]

CLASS_COLORS = [
    "tab:blue",
    "tab:green",
    "tab:orange",
    "tab:purple",
    "tab:brown",
    "tab:pink",
    "tab:olive",
    "tab:cyan",
]


@dataclass(frozen=True)
class Case:
    """
    One tangle to build and draw.

    Attributes:
        k: The Hénon ``k`` parameter (``b = 1`` throughout).
        unstable_steps: How many map steps to grow the unstable manifold.
        pip: ``"default"`` for the primary crossing q0, ``"image"`` for f(q0).
        blasts: How many single-iteration blasts of the zone to run after
            trimming.
        area_cutoff: Optional refinement area cutoff for the manifold machine.
    """

    k: float
    unstable_steps: int
    pip: Literal["default", "image"]
    blasts: int = 0
    area_cutoff: Optional[float] = None

    @property
    def title(self) -> str:
        """``"Hénon k=<k>"``."""
        return f"Hénon k={self.k:g}"


CASES = [
    Case(k=10, unstable_steps=10, pip="image"),
    Case(k=2.8, unstable_steps=10, pip="image", blasts=1, area_cutoff=1e-7),
]


def _stable_points(session: TangleSession, fp) -> np.ndarray:
    """Every node of the fixed point's stable manifold, stacked."""
    return np.vstack(
        [
            manifold.get_point_array()
            for (fixed_point, stability, _o, _b), manifold in session.manifolds.items()
            if fixed_point is fp and stability == "stable"
        ]
    )


def _view(point_sets: list[np.ndarray], pad: float = 0.06) -> tuple:
    """
    Axis limits framing the given point sets, padded by ``pad`` of the diagonal.

    The unstable manifold runs far off (x ~ 2500 on k=10), so the view is
    framed on the stable manifold and the drawn bridges, never on the whole
    unstable manifold.

    Args:
        point_sets: ``(N, 2)`` arrays to frame.
        pad: Padding as a fraction of the bounding box diagonal.

    Returns:
        ``((xmin, xmax), (ymin, ymax))``.
    """
    points = np.vstack(point_sets)
    mins, maxs = points.min(axis=0), points.max(axis=0)
    margin = pad * float(np.linalg.norm(maxs - mins))
    return (mins[0] - margin, maxs[0] + margin), (mins[1] - margin, maxs[1] + margin)


def _apply_view(ax, view: tuple) -> None:
    """Fix the limits and keep the plane's aspect by resizing the box, not the data."""
    ax.set_xlim(view[0])
    ax.set_ylim(view[1])
    ax.set_aspect("equal", adjustable="box")


def build(case: Case) -> tuple[TangleSession, object, np.ndarray]:
    """
    Run the pipeline for one case up to the stable partition.

    Args:
        case: The case to build.

    Returns:
        ``(session, fixed_point, stable_points)`` with the stable manifold's
        nodes captured BEFORE the trim at the pip, for framing.
    """
    session = TangleSession(
        henon_map(case.k, 1), henon_map_inverse(case.k, 1), henon_jacobian(case.k, 1)
    )
    if case.area_cutoff is not None:
        session.workbench._man_machine.area_cutoff = case.area_cutoff
    fp = session.construct_fixed_point(SADDLE_GUESS)
    session.orient_eigenvectors(
        fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
    )
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=case.unstable_steps)
    session.grow_until_turnaround(fp, "stable")
    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    session.infer_iterate_table()
    stable_points = _stable_points(session, fp)

    session.classify_strong_pips()
    trellis = session.trellis(fp)
    pip = trellis.strong_pip
    if case.pip == "image":
        pip = trellis.iterate(pip, 1)
        if pip is None:
            raise ValueError(f"{case.title}: the default strong pip has no registered image")
    logger.info("%s: candidates %s, chosen pip %d", case.title, trellis.strong_pip_candidates, pip)

    # Trim at the pip and recompute; the recompute rebuilds the trellis, so
    # classify and set the pip again. Same after every blast.
    zone = session.resonance_zone(pip)
    session.classify_strong_pips()
    session.set_strong_pip(fp, pip)
    for _ in range(case.blasts):
        session.blast_zone(zone, num_iterations=1, fixed_point=[fp], min_separation=1e-5)
        session.classify_strong_pips()
        session.set_strong_pip(fp, pip)

    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()
    return session, fp, stable_points


def representative_bridge(session: TangleSession, entry: BridgeClassEntry):
    """
    The member of a class with the smallest unstable canonical distance.

    Args:
        session: The session whose registry holds the crossings.
        entry: The class entry.

    Returns:
        ``(bridge, member)`` for the member whose first crossing has the
        smallest unstable cdist.
    """
    registry = session.workbench.intersection_registry
    member = min(
        entry.members, key=lambda m: registry[m.bridge_id[0]].unstable_cdist
    )
    return session.workbench.bridge(member.bridge_id), member


def _label_point(points: np.ndarray, view: tuple) -> np.ndarray:
    """
    Where to put a bridge's label: its midpoint if in view, else the middle of
    the part of the polyline that is in view (an escaping lobe leaves the frame).
    """
    midpoint = polyline_midpoint(points)
    (xmin, xmax), (ymin, ymax) = view
    inside = (
        (points[:, 0] >= xmin)
        & (points[:, 0] <= xmax)
        & (points[:, 1] >= ymin)
        & (points[:, 1] <= ymax)
    )
    if xmin <= midpoint[0] <= xmax and ymin <= midpoint[1] <= ymax or not inside.any():
        return midpoint
    visible = points[inside]
    return visible[len(visible) // 2]


def _symbol_text(entry: BridgeClassEntry, member) -> str:
    """Mathtext for the symbol a member contributes (``$a$``, ``$a^{-1}$``, ``inert``)."""
    if entry.letter is None:
        return "inert"
    if member.direction < 0:
        return f"${entry.letter}^{{-1}}$"
    return f"${entry.letter}$"


def draw_tangle(ax, session: TangleSession, fp, view: tuple, title: str) -> None:
    """Panel 1: the whole tangle with crossings and the strong pip."""
    plt.sca(ax)
    session.plot_tangle(fp, "unstable", color="tab:blue", linewidth=0.5)
    # Blasted images live in bridge children, not the grown manifold.
    session.workbench.plot_all_bridges()
    session.plot_tangle(fp, "stable", color="tab:red", linewidth=0.8)
    session.workbench.plot_intersections(fp, show_ids=True, ax=ax, s=12, zorder=15)
    session.plot_strong_pip(
        fp, ax=ax, s=120, facecolors="none", linewidths=1.5, zorder=20
    )
    _apply_view(ax, view)
    ax.set_title(f"{title} -- tangle")


def draw_classes(
    ax, session: TangleSession, fp, view: tuple, table: BridgeClassTable, title: str
) -> None:
    """Panel 2: one bridge per class, the smallest-cdist member, labelled by symbol."""
    plt.sca(ax)
    session.plot_tangle(fp, "unstable", color="lightgray", linewidth=0.5)
    session.plot_tangle(fp, "stable", color="tab:red", linewidth=0.8)
    session.plot_strong_pip(
        fp, ax=ax, s=120, facecolors="none", linewidths=1.5, zorder=20
    )

    handles = []
    for index, entry in enumerate(table):
        bridge, member = representative_bridge(session, entry)
        color = CLASS_COLORS[index % len(CLASS_COLORS)]
        style = "--" if entry.inert else "-"
        points = bridge.get_point_array()
        ax.plot(points[:, 0], points[:, 1], style, color=color, linewidth=2.0, zorder=10)
        for endpoint in (points[0], points[-1]):
            ax.scatter(*endpoint, color=color, s=18, zorder=11)
        text = _symbol_text(entry, member)
        ax.annotate(
            text,
            xy=_label_point(points, view),
            xytext=(6, 6),
            textcoords="offset points",
            color=color,
            fontsize=11,
            fontweight="bold",
            zorder=12,
        )
        state = "inert" if entry.inert else "active"
        zone = "no zone" if entry.zone_key is None else "in zone"
        handles.append(
            Line2D(
                [0],
                [0],
                color=color,
                linestyle=style,
                linewidth=2,
                label=f"{entry.letter or 'inert'}: {entry.bridge_class.label} "
                f"({state}, {zone}) bridge {member.bridge_id}",
            )
        )
    ax.legend(handles=handles, loc="best", fontsize=7, framealpha=0.85)
    _apply_view(ax, view)
    ax.set_title(f"{title} -- one bridge per class (smallest cdist)")


def main() -> None:
    """Build both cases, print their class tables, save the figure."""
    fig, axes = plt.subplots(len(CASES), 3, figsize=(21, 7 * len(CASES)))
    for row, case in zip(axes, CASES):
        session, fp, stable_points = build(case)
        table = session.bridge_classes()
        # Frame on the stable manifold and the bridges inside the zone; an
        # escaping lobe (an inert class, or an exterior class not yet seen to
        # loop) runs far off and is allowed to leave the view.
        framed = [e for e in table if e.zone_key is not None] or list(table.active)
        view = _view(
            [stable_points]
            + [representative_bridge(session, e)[0].get_point_array() for e in framed]
        )

        print(case.title)
        print("=" * len(case.title))
        print(session.describe_bridge_classes())
        print()

        ax_tangle, ax_classes, ax_partition = row
        draw_tangle(ax_tangle, session, fp, view, case.title)
        draw_classes(ax_classes, session, fp, view, table, case.title)
        session.plot_stable_partition(fp, ax=ax_partition)
        ax_partition.set_title(f"{case.title} -- stable partition")

    fig.tight_layout()
    FIGURES_DIR.mkdir(exist_ok=True)
    fig.savefig(OUT_PATH, dpi=150)
    plt.close(fig)
    logger.info("saved %s", OUT_PATH)


if __name__ == "__main__":
    main()
