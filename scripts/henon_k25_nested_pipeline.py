"""
The whole pipeline on the nested Hénon tangle at ``k = 2.5``, ``b = 1``: the
outer period-1 saddle and the inner period-3 orbit together.

The nested recipe (``henon_cases.build_nested``) is written for ``k = 2``;
this driver carries it to ``k = 2.5``:

* the period-3 orbit is continued from ``k = 2`` and grown as in
  ``henon_k25_period3_pipeline.py`` (13 / 9 steps scaled to the same stretch:
  9 / 6, each arm stopped early if it escapes);
* the period-1 saddle comes from the recipe's guess ``[4, -4]`` and
  orientation, its 11 unstable steps scaled the same way (by its own
  eigenvalue) and its stable side grown to turnaround;
* both are constructed outermost first (period 1 is ``A``, tangle 0), grown
  period 3 first, intersected, trimmed and bridged together; the inner zone is
  cut at its default strong pip, the outer at ``f(q0)``, the image of its
  default strong pip (the k=2.8 recipe's choice; ``--outer-pip default`` cuts
  it at ``q0``); the inner zone is blasted ``--inner-blasts`` times (default 6)
  and the outer ``--outer-blasts`` times (default 5; its blasts iterate the
  bridges of the annulus between the two tangles),
  each blast touching only its own tangle (``min_separation = 1e-5``, the
  period-3 runs'), and the topology is computed over both.

It prints the symbolic dynamics report and writes into
``figures/henon_k25_nested/nested_k25[_pipq0]_outerN_innerM/``:

* the six symbolic-itinerary panels of ``symbolic_figures.draw_panels``;
* ``trellis_original`` -- the trellis BEFORE minimization (both tangles);
* ``trellis_inner`` -- the minimal trellis zoomed on the inner tangle;
* ``combined`` -- the original and minimal trellis, the cartoon, the inner
  zoom (outlined on the minimal trellis), the itinerary table and the two
  transition graphs on one figure (the partitions stay their own panel).

The class figures (``henon_k25_period3_class_focus.py`` and
``henon_k25_period3_class_iterates.py``) take ``--case nested`` to draw this
run into the same folder.

Run with ``PYTHONPATH=src python3 scripts/henon_k25_nested_pipeline.py``
(``--outer-blasts 5``, ``--inner-blasts 6``, ``--outer-pip image``,
``--area-cutoff 1e-7``).
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent))

import henon_k25_period3_pipeline as p3  # noqa: E402
from henon_minimal_dual_graph import _apply_plane_view, _draw_stable  # noqa: E402
from henon_period3_k_sweep import grow_until_escape  # noqa: E402
from symbolic_figures import (  # noqa: E402
    FIGURES_DIR,
    _legend_below,
    configure_logging,
    draw_cartoon,
    draw_itinerary_table,
    draw_panels,
    draw_transition_graph,
    draw_trellis,
    plane_view,
    report,
    save_panel,
)

from tanglepack import TangleSession  # noqa: E402
from tanglepack.examples.henon import (  # noqa: E402
    HENON_P3,
    henon_jacobian,
    henon_map,
    henon_map_inverse,
    saddle_guesses,
)
from tanglepack.examples.henon_cases import PERIOD1_ORIENTATION  # noqa: E402
from tanglepack.topology import plotting  # noqa: E402
from tanglepack.topology.SymbolicDynamics import SymbolicDynamics  # noqa: E402

K = p3.K
B = p3.B
AREA_CUTOFF = p3.AREA_CUTOFF
MIN_SEPARATION = p3.MIN_SEPARATION
#: The k=2 recipe's growth: period 3 13 / 9, period 1 11 unstable (stable to turnaround).
P3_UNSTABLE, P3_STABLE, P1_UNSTABLE = 13, 9, 11
#: The default run: the outer zone cut at f(q0) and blasted 5 times, the inner 6 times.
OUTER_PIP = "image"
OUTER_BLASTS = 5
INNER_BLASTS = 6
#: The folder every run's panel folder goes in.
OUT_DIR = FIGURES_DIR / "henon_k25_nested"


def _period1_lambda(k: float) -> float:
    """The outer saddle's unstable eigenvalue at ``k`` (``x* = 1 + sqrt(1 + k)``)."""
    x = 1.0 + math.sqrt(1.0 + k)
    return x + math.sqrt(x * x - 1.0)


def build(
    *, outer_blasts: int, inner_blasts: int, area_cutoff: float = AREA_CUTOFF,
    outer_pip: str = OUTER_PIP,
) -> tuple[TangleSession, list, dict[str, int]]:
    """
    Run the pipeline on both tangles at ``K``.

    Args:
        outer_blasts: Single-iteration blasts of the period-1 zone.
        inner_blasts: Single-iteration blasts of the period-3 zone.
        area_cutoff: Refinement area cutoff.
        outer_pip: The period-1 zone's pip: ``"image"`` = f(q0) of the default
            strong pip q0 (the k=2.8 recipe), ``"default"`` = q0.

    Raises:
        ValueError: ``outer_pip == "image"`` and q0 has no registered image.

    Returns:
        ``(session, [fp1, fp3], steps grown)``, the topology computed.
    """
    guess3, orientation3, lambda_seed3 = p3.continue_orbit(K)
    session = TangleSession(
        henon_map(K, B), henon_map_inverse(K, B), henon_jacobian(K, B),
        area_cutoff=area_cutoff,
    )
    fp1 = session.construct_fixed_point(saddle_guesses(*HENON_P3)["saddle"])
    session.orient_eigenvectors(fp1, PERIOD1_ORIENTATION)
    fp3 = session.construct_fixed_point(guess3)
    session.orient_eigenvectors(fp3, orientation3)
    fixed_points = [fp1, fp3]

    grown: dict[str, int] = {}
    session.initialize_both_manifolds(fp3)
    factor3 = math.log(lambda_seed3) / math.log(float(fp3.unstable_eigenvalues[0]))
    radius = 2.0 * (1.0 + math.sqrt(1.0 + K))
    for stability, count in (("unstable", P3_UNSTABLE), ("stable", P3_STABLE)):
        steps, escaped = grow_until_escape(
            session, fp3, stability, radius=radius,
            max_steps=max(1, round(count * factor3)), step_budget=p3.STEP_BUDGET,
        )
        grown[f"p3 {stability}"] = steps
        print(f"period 3 {stability}: {steps} step(s){' (escaped)' if escaped else ''}")

    session.initialize_both_manifolds(fp1)
    factor1 = math.log(_period1_lambda(HENON_P3[0])) / math.log(float(fp1.unstable_eigenvalues[0]))
    grown["p1 unstable"] = max(1, round(P1_UNSTABLE * factor1))
    session.grow_n_times(fp1, "unstable", num_iterations=grown["p1 unstable"])
    session.grow_until_turnaround(fp1, "stable")
    print(f"period 1 unstable: {grown['p1 unstable']} step(s); stable to turnaround")

    session.compute_intersections(fixed_points)
    for fp in fixed_points:
        session.trim_stable_manifolds(fp)
    for fp in fixed_points:
        session.create_bridges(fp)
    session.classify_strong_pips()
    pips = [session.trellis(fp).strong_pip for fp in fixed_points]
    if outer_pip == "image":
        image = session.trellis(fp1).iterate(pips[0], 1)
        if image is None:
            raise ValueError("the outer default strong pip q0 has no registered image f(q0)")
        print(f"outer pip: f(q0) = crossing {image} (q0 = crossing {pips[0]})")
        pips[0] = image
    session.add_resonance_zones(pips, nesting={pips[1]: pips[0]})
    chosen = session.chosen_strong_pip(fp1)
    print(f"outer chosen strong pip after the trim: cdists (u, s) = ({chosen.unstable_cdist:.6g}, {chosen.stable_cdist:.6g})")
    zones = sorted(session.resonance_zones.values(), key=lambda zone: zone.area)
    inner, outer = zones[0], zones[-1]
    for name, zone, count in (("inner", inner, inner_blasts), ("outer", outer, outer_blasts)):
        for step in range(1, count + 1):
            result = session.blast_zone(
                zone, num_iterations=1, fixed_point=[zone.fixed_point],
                min_separation=MIN_SEPARATION,
            )
            print(
                f"{name} blast {step}: {len(result.all_interior_bridges())} interior bridge(s), "
                f"pip depth {result.pip_depth.get(zone.fixed_point)}"
            )
    session.topology()
    return session, fixed_points, grown


def run_slug(outer_blasts: int, inner_blasts: int, area_cutoff: float = AREA_CUTOFF,
             outer_pip: str = OUTER_PIP) -> str:
    """``nested_k25_outer5_inner6``, ``nested_k25_pipq0_outer2_inner6``, ...: the folder name of a run."""
    slug = "nested_k25"
    if outer_pip == "default":
        slug += "_pipq0"
    if area_cutoff != AREA_CUTOFF:
        slug += f"_a{area_cutoff:.0e}".replace("e-0", "e-")
    return slug + f"_outer{outer_blasts}_inner{inner_blasts}"


def run_title(grown: dict[str, int], outer_blasts: int, inner_blasts: int, area_cutoff: float,
              outer_pip: str = OUTER_PIP) -> str:
    """The run's title."""
    pip = "f(q0)" if outer_pip == "image" else "q0"
    return (
        f"Hénon k={K:g}, b={B}, period 1 + period 3 (p3 {grown['p3 unstable']}u/"
        f"{grown['p3 stable']}s, p1 {grown['p1 unstable']}u, cutoff {area_cutoff:g}, outer pip {pip}), "
        f"{inner_blasts} inner + {outer_blasts} outer blast(s)"
    )


def inner_view(session: TangleSession, fp, pad: float = 0.12) -> tuple:
    """A square view around one orbit's (trimmed) stable manifold and its orbit points."""
    points = [np.atleast_2d(np.asarray(c, dtype=float)) for c in fp.coordinates]
    for (owner, stability, _o, _b), manifold in session.manifolds.items():
        if owner is fp and stability == "stable":
            points.append(manifold.get_point_array())
    points = np.vstack(points)
    lo, hi = points.min(axis=0), points.max(axis=0)
    center, half = 0.5 * (lo + hi), 0.5 * float(np.max(hi - lo)) * (1 + 2 * pad)
    return ((center[0] - half, center[0] + half), (center[1] - half, center[1] + half))


def draw_inner_trellis(ax, session: TangleSession, fixed_points, dynamics: SymbolicDynamics,
                       title: str, *, outer_ax=None) -> None:
    """
    The minimal trellis zoomed on the inner period-3 tangle (its own plane view);
    with ``outer_ax`` the zoomed box is outlined on that (whole-trellis) panel.
    """
    view = inner_view(session, fixed_points[-1])
    draw_trellis(ax, session, fixed_points, dynamics, view, title)
    ax.set_title(f"{title} -- minimal trellis, zoomed on the inner period-3 tangle", fontsize=11)
    if outer_ax is not None:
        (x0, x1), (y0, y1) = view
        outer_ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, edgecolor="#555555",
                                     linewidth=1.2, linestyle="--", zorder=25))


def draw_original_trellis(ax, session: TangleSession, fixed_points, dynamics: SymbolicDynamics,
                          view, title) -> None:
    """The trellis before minimization (both tangles): every bridge, by refined class."""
    plt.sca(ax)
    minimal = session.minimal_trellis()
    trellis = session.trellis()
    for fp in fixed_points:
        plt.sca(ax)
        session.plot_tangle(fp, "unstable", color="lightgray", linewidth=0.6)
        _draw_stable(ax, session, fp)
    drawn = plotting.plot_bridges_by_class(trellis, dynamics, ax=ax, refined=True)
    trellis.plot_holes(ax=ax)
    for fp in fixed_points:
        session.trellis(fp).plot_strong_pip(ax=ax, s=120, facecolors="none", linewidths=1.5, zorder=20)
    legend = ax.get_legend()
    if legend is not None:
        legend.remove()
    shown = {plotting.name_mathtext(name) for name in drawn}
    handles = [
        handle for handle in plotting.bridge_class_legend_handles(dynamics, refined=True)
        if handle.get_label().split(" (")[0] in shown
    ]
    _legend_below(ax, handles)
    _apply_plane_view(ax, view)
    ax.set_title(
        f"{title} -- original trellis before minimization "
        f"({len(trellis.bridges)} bridges; minimal keeps {len(minimal.kept_bridge_ids)}), "
        "by refined class",
        fontsize=11,
    )


def draw_combined(session: TangleSession, fixed_points, dynamics: SymbolicDynamics, title: str,
                  folder: Path) -> Path:
    """All seven panels on one figure (the period-3 driver's layout)."""
    view = plane_view(session, fixed_points, session.minimal_trellis())
    full_title, title = title, f"nested k={K:g}"  # the suptitle carries the run
    table_height = 0.45 * len(dynamics.classes) + 2.0
    fig = plt.figure(figsize=(36, 13 + max(13.0, table_height)))
    grid = fig.add_gridspec(2, 6, height_ratios=[13, max(13.0, table_height)])
    draw_original_trellis(fig.add_subplot(grid[0, 0:2]), session, fixed_points, dynamics, view, title)
    whole = fig.add_subplot(grid[0, 2:4])
    draw_trellis(whole, session, fixed_points, dynamics, view, title)
    draw_cartoon(fig.add_subplot(grid[0, 4:6]), session, dynamics, title)
    draw_inner_trellis(fig.add_subplot(grid[1, 0:2]), session, fixed_points, dynamics, title,
                       outer_ax=whole)
    table_ax = fig.add_subplot(grid[1, 2:4])
    table_ax.set_axis_off()
    table_ax.set_title(f"{title} -- itineraries and words", fontsize=11)
    for column, refined, name in ((4, False, "class symbols"), (5, True, "refined symbols")):
        ax = fig.add_subplot(grid[1, column])
        draw_transition_graph(ax, dynamics, title, refined=refined)
        ax.set_title(f"transition graph ({name})", fontsize=11)
    fig.suptitle(full_title, fontsize=18)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    draw_itinerary_table(table_ax, dynamics)
    return save_panel(fig, folder, "combined")


def run(outer_blasts: int, inner_blasts: int, area_cutoff: float = AREA_CUTOFF,
        outer_pip: str = OUTER_PIP) -> Path:
    """Build one run, print its report, draw and save its panels; returns its folder."""
    session, fixed_points, grown = build(
        outer_blasts=outer_blasts, inner_blasts=inner_blasts, area_cutoff=area_cutoff,
        outer_pip=outer_pip,
    )
    dynamics = session.symbolic_dynamics()
    title = run_title(grown, outer_blasts, inner_blasts, area_cutoff, outer_pip)
    report(session, dynamics, title)
    p3.closure_summary(dynamics)

    folder = OUT_DIR / run_slug(outer_blasts, inner_blasts, area_cutoff, outer_pip)
    draw_panels(session, fixed_points, title, folder, dynamics=dynamics)
    fig, ax = plt.subplots(figsize=(12, 10))
    draw_inner_trellis(ax, session, fixed_points, dynamics, title)
    fig.tight_layout()
    save_panel(fig, folder, "trellis_inner")
    fig, ax = plt.subplots(figsize=(12, 10))
    draw_original_trellis(ax, session, fixed_points, dynamics,
                          plane_view(session, fixed_points, session.minimal_trellis()), title)
    fig.tight_layout()
    save_panel(fig, folder, "trellis_original")
    draw_combined(session, fixed_points, dynamics, title, folder)
    print(f"wrote {folder}")
    return folder


def add_run_arguments(parser: argparse.ArgumentParser) -> None:
    """The run's options (shared with the class figures' ``--case nested``)."""
    parser.add_argument("--outer-blasts", type=int, default=OUTER_BLASTS,
                        help="blasts of the period-1 zone (the annulus between the tangles)")
    parser.add_argument("--inner-blasts", type=int, default=INNER_BLASTS, help="blasts of the period-3 zone")
    parser.add_argument("--outer-pip", choices=("image", "default"), default=OUTER_PIP,
                        help="the period-1 zone's pip: f(q0) ('image') or q0 ('default')")
    parser.add_argument("--area-cutoff", type=float, default=AREA_CUTOFF)


def main() -> None:
    parser = argparse.ArgumentParser(description="The whole pipeline on the nested Hénon tangle at k=2.5.")
    add_run_arguments(parser)
    args = parser.parse_args()
    configure_logging()
    run(args.outer_blasts, args.inner_blasts, args.area_cutoff, args.outer_pip)


if __name__ == "__main__":
    main()
