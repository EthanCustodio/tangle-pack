"""
The dual graph of the Hénon tangle: k=10 and the nested period-3 map.

Builds the two reference sessions used throughout the topology test suite (the
single k=10 saddle at ``[4, -4]``, and the nested period-3-inside-period-1
tangle), partitions each fully, and builds its
:class:`~tanglepack.topology.DualGraph.DualGraph` -- one node per face, one
per stable arc, the arcs between the strong pip and its ``k``-th iterate
filled -- without growing anything further. For each fixture this logs the
graph's summary and fill segments and saves a figure with the tangle and its
dual graph drawn twice: straight edges from each face point to its arc nodes
(:func:`~tanglepack.topology.plotting.plot_dual_graph`) and edges following
the unstable manifold inside each face
(:func:`~tanglepack.topology.plotting.plot_dual_graph_curved`).

Run:  env/bin/python scripts/henon_dual_graph.py

Dev Notes:
    ``build_k10``/``build_p3`` mirror the ``k10_partitioned``/``p3_partitioned``
    fixtures in ``tests/conftest.py`` exactly (down to the ``area_cutoff``
    overrides on the p3 fixed points) rather than importing them, since a
    script should not depend on the test tree. Recipe drift between the two
    is a risk worth watching if either changes.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from tanglepack import TangleSession
from tanglepack.examples import (
    HENON_K10,
    HENON_P3,
    henon_jacobian as _henon_jacobian_factory,
    henon_map as _henon_map_factory,
    henon_map_inverse as _henon_map_inverse_factory,
    saddle_guesses,
)
from tanglepack.numerics.FixedPoint import FixedPoint
from tanglepack.topology.DualGraph import DualGraph
from tanglepack.topology.plotting import dual_graph_legend_handles

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"


def build_k10() -> tuple[TangleSession, list[FixedPoint]]:
    """
    Build and fully partition the k=10 single-saddle session.

    Mirrors the ``k10_session`` / ``k10_partitioned`` fixtures of
    ``tests/conftest.py``: the saddle at ``[4, -4]``, unstable grown 9 steps,
    stable grown to turnaround, then intersections, trimming, bridges, the
    iterate table, and the four partition fan-outs.

    Returns:
        ``(session, [fp])``.
    """
    henon_map = _henon_map_factory(*HENON_K10)
    henon_map_inverse = _henon_map_inverse_factory(*HENON_K10)

    session = TangleSession(henon_map, henon_map_inverse)
    fp = session.construct_fixed_point([4, -4])
    session.orient_eigenvectors(
        fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
    )
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=9)
    session.grow_until_turnaround(fp, "stable")
    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    session.infer_iterate_table()

    session.classify_strong_pips()
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()
    return session, [fp]


def build_p3() -> tuple[TangleSession, list[FixedPoint]]:
    """
    Build and fully partition the nested period-3-inside-period-1 session.

    Mirrors the ``henon_p3_session`` / ``p3_partitioned`` fixtures of
    ``tests/conftest.py`` (equivalently ``scripts/henon_blast_period_3.py``),
    minus the resonance-zone step this script has no use for: the period-3
    orbit and the outer period-1 saddle are grown and crossed, then intersected,
    trimmed, bridged, the iterate table inferred, and the four partition
    fan-outs run over both fixed points.

    Returns:
        ``(session, [fp3, fp1])``.
    """
    henon_map = _henon_map_factory(*HENON_P3)
    henon_map_inverse = _henon_map_inverse_factory(*HENON_P3)
    henon_jacobian = _henon_jacobian_factory(*HENON_P3)

    session = TangleSession(henon_map, henon_map_inverse, henon_jacobian)

    session.workbench._man_machine.area_cutoff = 1e-7
    fp3 = session.construct_fixed_point(saddle_guesses(*HENON_P3)["period_3"])
    session.orient_eigenvectors(
        fp3, {"unstable": np.array([0, -1]), "stable": np.array([-1, -1])}
    )
    session.initialize_both_manifolds(fp3)
    session.grow_n_times(fp3, "unstable", num_iterations=13)
    session.grow_n_times(fp3, "stable", num_iterations=9)

    session.workbench._man_machine.area_cutoff = 1e-4
    fp1 = session.construct_fixed_point(saddle_guesses(*HENON_P3)["saddle"])
    session.orient_eigenvectors(
        fp1, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
    )
    session.initialize_both_manifolds(fp1)
    session.grow_n_times(fp1, "unstable", num_iterations=11)
    session.grow_until_turnaround(fp1, "stable")

    session.compute_intersections([fp3, fp1])
    session.trim_stable_manifolds(fp3)
    session.trim_stable_manifolds(fp1)
    session.create_bridges(fp3)
    session.create_bridges(fp1)
    session.infer_iterate_table()

    session.classify_strong_pips()
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()
    return session, [fp3, fp1]


def _manifold_bbox(
    session: TangleSession, fixed_point: FixedPoint
) -> tuple[float, float, float, float]:
    """
    The bounding box of every point on one fixed point's manifolds.

    Args:
        session: The session whose workbench holds the manifolds.
        fixed_point: The fixed point (matched by IDENTITY) to bound.

    Returns:
        ``(xmin, xmax, ymin, ymax)``.

    Raises:
        ValueError: If the fixed point has no manifolds in this session.
    """
    arrays = [
        manifold.get_point_array()
        for (fp, _stab, _orbit, _branch), manifold in
        session.workbench.manifolds.items()
        if fp is fixed_point
    ]
    if not arrays:
        raise ValueError(
            f"fixed point {fixed_point} has no manifolds in this session"
        )
    points = np.vstack(arrays)
    xmin, ymin = points.min(axis=0)
    xmax, ymax = points.max(axis=0)
    return float(xmin), float(xmax), float(ymin), float(ymax)


def _draw_tangle(session: TangleSession, fixed_points: list[FixedPoint], ax) -> None:
    """Every fixed point's manifolds on one axes: unstable blue, stable red."""
    plt.sca(ax)
    for fp in fixed_points:
        session.plot_tangle(fp, "unstable", color="tab:blue", linewidth=0.6)
        session.plot_tangle(fp, "stable", color="tab:red", linewidth=0.6)


def _zoom(ax, session: TangleSession, fixed_point: FixedPoint) -> None:
    """Limit ``ax`` to one fixed point's manifolds, padded 15%."""
    xmin, xmax, ymin, ymax = _manifold_bbox(session, fixed_point)
    pad_x = 0.15 * (xmax - xmin)
    pad_y = 0.15 * (ymax - ymin)
    ax.set_xlim(xmin - pad_x, xmax + pad_x)
    ax.set_ylim(ymin - pad_y, ymax + pad_y)


def report_and_plot(
    session: TangleSession,
    fixed_points: list[FixedPoint],
    title: str,
    out_path: Path,
    *,
    zoom_fixed_point: Optional[FixedPoint] = None,
) -> DualGraph:
    """
    Log a fixture's dual graph and save its figure.

    Reads ``session.dual_graph()`` (cached over every fixed point of the
    session, which is what ``fixed_points`` names), logs its one-line
    :meth:`~tanglepack.topology.DualGraph.DualGraph.summary` and prints the
    fill segment of every strong pip's branch.

    Without ``zoom_fixed_point`` the figure is a 1x2 layout: every fixed
    point's manifolds with the dual graph overlaid, straight edges on the
    left and curved edges on the right. With ``zoom_fixed_point`` (the nested
    p3 fixture) it is 2x2: the top row is the full tangle, the bottom row is
    zoomed on that fixed point's own manifolds (bounding box padded 15%, dual
    graph unclipped) -- otherwise the inner period-3 tangle is an illegible
    blob at the outer period-1 tangle's scale -- with straight edges on the
    left and curved on the right.

    Args:
        session: A fully partitioned session (see ``build_k10``/``build_p3``).
        fixed_points: The session's fixed points, for the manifold plot.
        title: Figure/report title.
        out_path: Where to save the PNG.
        zoom_fixed_point: When given, adds the zoomed bottom row around this
            fixed point (matched by identity).

    Returns:
        The :class:`~tanglepack.topology.DualGraph.DualGraph`.
    """
    dual_graph = session.dual_graph()
    logger.info("%s: %s", title, dual_graph.summary())

    print(f"{title}")
    print("=" * len(title))
    print(dual_graph.summary())
    for branch_key, (low, high) in dual_graph.fill_segments.items():
        print(
            f"  fill on stable branch {branch_key[1:]} (period "
            f"{branch_key[0].period}): ({low:.4g}, {high:.4g}]"
        )

    if zoom_fixed_point is None:
        fig, (ax_straight, ax_curved) = plt.subplots(1, 2, figsize=(16, 8))
        panels = [(ax_straight, False, False), (ax_curved, True, False)]
    else:
        fig, ((ax_straight, ax_curved), (ax_zs, ax_zc)) = plt.subplots(
            2, 2, figsize=(18, 16)
        )
        panels = [
            (ax_straight, False, False), (ax_curved, True, False),
            (ax_zs, False, True), (ax_zc, True, True),
        ]

    for ax, curved, zoomed in panels:
        _draw_tangle(session, fixed_points, ax)
        plotter = session.plot_dual_graph_curved if curved else session.plot_dual_graph
        plotter(dual_graph, ax=ax, clip_to_arcs=not zoomed)
        if zoomed:
            _zoom(ax, session, zoom_fixed_point)
        what = "curved edges" if curved else "straight edges"
        where = (
            f"zoom on period-{zoom_fixed_point.period} tangle" if zoomed else "tangle"
        )
        ax.set_title(f"{title} -- {where} + dual graph, {what}")
        ax.legend(
            handles=dual_graph_legend_handles(curved=curved),
            loc="best", fontsize=7, framealpha=0.8,
        )

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    logger.info("saved %s", out_path)

    return dual_graph


def main() -> None:
    """Build both fixtures, report their dual graphs, save both figures."""
    session_k10, fps_k10 = build_k10()
    report_and_plot(
        session_k10, fps_k10, "Hénon k=10",
        FIGURES_DIR / "henon_dual_graph_k10.png",
    )

    session_p3, fps_p3 = build_p3()
    fp3, _fp1 = fps_p3
    report_and_plot(
        session_p3, fps_p3, "Hénon nested period-3",
        FIGURES_DIR / "henon_dual_graph_p3.png",
        zoom_fixed_point=fp3,
    )


if __name__ == "__main__":
    main()
