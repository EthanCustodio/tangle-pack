"""
The minimal trellis, the iterated homotopy partition and the new dual graph on
the Hénon cases of ``henon_bridge_classes.py``.

One row per case, three panels:

1. the tangle with the MINIMAL trellis highlighted — hole bridges solid,
   image bridges dashed, dropped bridges light gray — plus the holes and the
   chosen strong pip;
2. the homotopy partition (the rows the holes define) above the iterated
   homotopy partition (the same rows cut by the image bridges, ``[c1, c2]``
   closed under each image bridge on its own side), as number lines;
3. the tangle with the dual graph: an open circle just off each stable edge
   on each side (a wall), a solid circle on each edge of the fundamental
   segment (traversable), a dot per face joined to the nodes on its side.

The cases are imported from ``henon_bridge_classes.py`` (same recipes, same
``build``), so the two figures line up row for row. The third case is the
k=2.8 tangle blasted twice.

Run with ``PYTHONPATH=src python3 scripts/henon_minimal_dual_graph.py``; the
figure lands in ``figures/henon_minimal_dual_graph.png`` and each case's
minimal trellis, iterated partition and dual-graph summary are printed.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from henon_bridge_classes import CASES, Case, _apply_view, _view, build  # noqa: E402

from tanglepack import TangleSession  # noqa: E402
from tanglepack.topology import plotting  # noqa: E402

logger = logging.getLogger("henon_minimal_dual_graph")

FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"
OUT_PATH = FIGURES_DIR / "henon_minimal_dual_graph.png"


def _apply_plane_view(ax, view: tuple) -> None:
    """The shared view of a row's plane panels, WITHOUT the equal aspect.

    A k=10 bridge runs to x ~ 180 while the tangle sits within |y| < 20, so an
    equal-aspect panel is a thin strip in which nothing near the fixed point is
    legible; the dual graph is a topological picture and survives the stretch.
    """
    _apply_view(ax, view)
    ax.set_aspect("auto")


def _draw_stable(ax, session: TangleSession, fp) -> None:
    """The stable manifold, in red, as the backdrop of the plane panels."""
    plt.sca(ax)  # the workbench plotter draws on the current axes
    session.plot_tangle(fp, "stable", color="tab:red", linewidth=0.8)


def draw_minimal_trellis(ax, session: TangleSession, fp, view: tuple, title: str) -> None:
    """Panel 1: the tangle with the minimal trellis over it."""
    _draw_stable(ax, session, fp)
    minimal = session.minimal_trellis()
    plotting.plot_minimal_trellis(minimal, ax=ax)
    session.plot_holes(fp, ax=ax)
    session.plot_strong_pip(
        fp, ax=ax, s=120, facecolors="none", linewidths=1.5, zorder=20
    )
    ax.legend(handles=plotting.minimal_trellis_legend_handles(), loc="upper right", fontsize=8)
    _apply_plane_view(ax, view)
    ax.set_title(f"{title} -- minimal trellis", fontsize=10)


def draw_partitions(ax, session: TangleSession, title: str) -> None:
    """Panel 2: homotopy rows above iterated rows, one axes."""
    homotopy = session.homotopy_partition()
    iterated = session.iterated_partition()
    results = homotopy.as_list() + iterated.as_list()
    labels = [
        f"{family} {result.side} (orbit {result.branch_key[2]}.{result.branch_key[3]})"
        for family, family_results in (("homotopy", homotopy), ("iterated", iterated))
        for result in family_results
    ]
    plotting.plot_stable_partition(results, ax=ax, labels=labels)
    ax.set_title(f"{title} -- homotopy vs iterated homotopy partition", fontsize=9)


def draw_dual_graph(ax, session: TangleSession, fp, view: tuple, title: str) -> None:
    """Panel 3: the tangle with the dual graph over it."""
    _draw_stable(ax, session, fp)
    minimal = session.minimal_trellis()
    plotting.plot_minimal_trellis(
        minimal, ax=ax, show_dropped=False, color="lightgray", linewidth=0.8, linestyle="-"
    )
    dual = session.dual_graph()
    session.plot_dual_graph(dual, ax=ax, clip_to_arcs=False)
    ax.legend(handles=plotting.dual_graph_legend_handles(), loc="upper right", fontsize=8)
    _apply_plane_view(ax, view)
    ax.set_title(f"{title} -- dual graph", fontsize=10)


def report(session: TangleSession, title: str) -> None:
    """Print the three objects' reports for one case."""
    print(title)
    print("=" * len(title))
    print(session.minimal_trellis().describe())
    print(session.describe_iterated_partition())
    print(session.dual_graph().summary())
    print()


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    fig, axes = plt.subplots(
        len(CASES), 3, figsize=(26, 7 * len(CASES)),
        gridspec_kw={"width_ratios": [1.3, 1.0, 1.3]},
    )
    axes = np.atleast_2d(axes)
    for row, case in enumerate(CASES):
        session, fp, stable_points = build(case)
        minimal = session.minimal_trellis()
        view = _view(
            [stable_points]
            + [bridge.get_point_array() for bridge in minimal.kept]
        )
        report(session, case.title)
        draw_minimal_trellis(axes[row, 0], session, fp, view, case.title)
        draw_partitions(axes[row, 1], session, case.title)
        draw_dual_graph(axes[row, 2], session, fp, view, case.title)

    fig.tight_layout()
    FIGURES_DIR.mkdir(exist_ok=True)
    fig.savefig(OUT_PATH, dpi=150)
    logger.info("wrote %s", OUT_PATH)


if __name__ == "__main__":
    main()
