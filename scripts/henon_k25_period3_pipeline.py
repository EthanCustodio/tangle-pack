"""
The whole pipeline on the period-3 Hénon orbit at ``k = 2.5``, ``b = 1``.

The period-3 recipe (``henon_cases.build_period3``) is written for ``k = 2``;
this driver carries it to ``k = 2.5``:

* the orbit is continued from ``k = 2`` in steps of 0.1, each solve seeded by
  the previous orbit and oriented by the previous ORIENTED eigenvectors (the
  ``k = 2`` seed is ``build_period3``'s guess and orientation), as in
  ``henon_period3_k_sweep.py``;
* the growth is the ``k = 2`` recipe's 13 / 9 steps scaled to the same total
  stretch, ``steps * ln lambda_u(2) / ln lambda_u(2.5)`` (9 / 6), each
  manifold stopped early if its arm escapes (``2 (1 + sqrt(1 + k))``) or a
  step is predicted over the point budget;
* ``area_cutoff = 1e-7``, crossings, stable trim, bridges, the default strong
  pip, the resonance zone at it, ``--blasts`` single-iteration blasts
  (``min_separation = 1e-5``) and the topology, as ``build_period3`` does.

It prints the symbolic dynamics report and writes into
``figures/henon_k25_period3/p3_k25[_uN][_sN][_blastN]/``:

* the six symbolic-itinerary panels of ``symbolic_figures.draw_panels``
  (``trellis`` -- the MINIMAL trellis --, ``partitions``, ``cartoon``,
  ``itineraries``, ``transitions``, ``transitions_refined``);
* ``trellis_original`` -- the trellis BEFORE minimization: the grown
  unstable manifold (grey), every bridge of the trellis, blast children
  included, coloured by refined class (inert dashed), the holes and the
  strong pips, framed like the minimal ``trellis`` panel;
* ``combined`` -- all seven panels on one figure.

Every file is written as PNG, SVG and a pickled figure (reopen a pickle with
``scripts/show_figure.py``).

Run with ``PYTHONPATH=src python3 scripts/henon_k25_period3_pipeline.py``
(``--blasts N [N ...]``, one folder per count, default 12 and 24; ``--unstable-steps`` / ``--stable-steps``, the
``k = 2`` counts before scaling, default 13 / 9; ``--no-scale`` grows them
verbatim).

Dev Notes:
    Measured 2026-10-19 at the default growth (9 / 6 steps, neither arm
    escapes): up to 5 blasts the trellis is closed with 5 classes and the
    three active symbols ``a, b, c`` (the ``k = 2`` regime); at 6 blasts the
    reference switches, as it does at ``k = 2``, and the alphabet holds at 28
    classes (27 refined symbols on 20 minimal-trellis bridges, reliable) to
    12 blasts while the full trellis grows (27 -> 46 bridges).
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from henon_minimal_dual_graph import _apply_plane_view, _draw_stable  # noqa: E402
from henon_period3_k_sweep import grow_until_escape  # noqa: E402
from symbolic_figures import (  # noqa: E402
    FIGURES_DIR,
    _legend_below,
    configure_logging,
    draw_cartoon,
    draw_itinerary_table,
    draw_panels,
    draw_partitions,
    draw_transition_graph,
    draw_trellis,
    plane_view,
    report,
    save_panel,
)

import tanglepack  # noqa: E402
from tanglepack import TangleSession  # noqa: E402
from tanglepack.examples.henon import (  # noqa: E402
    henon_jacobian,
    henon_map,
    henon_map_inverse,
    saddle_guesses,
)
from tanglepack.examples.henon_cases import PERIOD3_ORIENTATION  # noqa: E402
from tanglepack.topology import plotting  # noqa: E402
from tanglepack.topology.SymbolicDynamics import SymbolicDynamics  # noqa: E402

K = 2.5
B = 1
#: The recipe's own parameter and the continuation step towards ``K``.
K_SEED = 2.0
K_STEP = 0.1
AREA_CUTOFF = 1e-7
MIN_SEPARATION = 1e-5
STEP_BUDGET = 2_000_000
#: The folder every run's panel folder goes in.
OUT_DIR = FIGURES_DIR / "henon_k25_period3"


def continue_orbit(k: float) -> tuple[list, dict, float]:
    """
    Continue the period-3 orbit from ``K_SEED`` to ``k``.

    Args:
        k: The target parameter.

    Returns:
        ``(orbit guess, orientation hints, lambda_u at K_SEED)``: the orbit at
        ``k`` (as a multipoint guess), its oriented eigenvectors, and the seed's
        unstable eigenvalue (the step scaling's reference).
    """
    guess = saddle_guesses(K_SEED, B)["period_3"]
    orientation = dict(PERIOD3_ORIENTATION)
    lambda_seed = None
    for kk in np.round(np.arange(K_SEED, k + 1e-9, K_STEP), 10):
        session = tanglepack.TangleSession(
            henon_map(kk, B), henon_map_inverse(kk, B), henon_jacobian(kk, B)
        )
        fp = session.construct_fixed_point(guess)
        session.orient_eigenvectors(fp, orientation)
        if lambda_seed is None:
            lambda_seed = float(fp.unstable_eigenvalues[0])
        guess = [c.tolist() for c in fp.coordinates]
        orientation = {
            "unstable": np.array(fp.unstable_eigenvectors[0], dtype=float),
            "stable": np.array(fp.stable_eigenvectors[0], dtype=float),
        }
    return guess, orientation, lambda_seed


def build(
    *, unstable_steps: int, stable_steps: int, blasts: int, scale: bool,
    area_cutoff: float = AREA_CUTOFF,
) -> tuple[TangleSession, object, dict[str, int]]:
    """
    Run the pipeline on the period-3 orbit at ``K``.

    Args:
        unstable_steps: Unstable steps of the ``K_SEED`` recipe.
        stable_steps: Stable steps of the ``K_SEED`` recipe.
        blasts: Single-iteration blasts of the resonance zone.
        scale: Scale the step counts to the ``K_SEED`` recipe's stretch.

    Returns:
        ``(session, orbit, steps grown per stability)``, the topology computed.
    """
    guess, orientation, lambda_seed = continue_orbit(K)
    session = TangleSession(
        henon_map(K, B), henon_map_inverse(K, B), henon_jacobian(K, B),
        area_cutoff=area_cutoff,
    )
    fp = session.construct_fixed_point(guess)
    session.orient_eigenvectors(fp, orientation)
    session.initialize_both_manifolds(fp)
    factor = math.log(lambda_seed) / math.log(float(fp.unstable_eigenvalues[0])) if scale else 1.0
    radius = 2.0 * (1.0 + math.sqrt(1.0 + K))
    grown: dict[str, int] = {}
    for stability, count in (("unstable", unstable_steps), ("stable", stable_steps)):
        grown[stability], escaped = grow_until_escape(
            session, fp, stability, radius=radius,
            max_steps=max(1, round(count * factor)), step_budget=STEP_BUDGET,
        )
        print(f"{stability}: {grown[stability]} step(s){' (escaped)' if escaped else ''}")

    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    session.classify_strong_pips()
    zone = session.resonance_zone(session.trellis(fp).strong_pip)
    for step in range(1, blasts + 1):
        result = session.blast_zone(
            zone, num_iterations=1, fixed_point=[fp], min_separation=MIN_SEPARATION
        )
        print(
            f"blast {step}: {len(result.all_interior_bridges())} interior bridge(s), "
            f"pip depth {result.pip_depth.get(fp)}"
        )
    session.topology()
    return session, fp, grown


def draw_original_trellis(ax, session: TangleSession, fp, dynamics: SymbolicDynamics, view, title) -> None:
    """The trellis before minimization: every bridge, coloured by refined class."""
    plt.sca(ax)
    minimal = session.minimal_trellis()
    trellis = session.trellis(fp)
    session.plot_tangle(fp, "unstable", color="lightgray", linewidth=0.6)
    _draw_stable(ax, session, fp)
    drawn = plotting.plot_bridges_by_class(trellis, dynamics, ax=ax, refined=True)
    trellis.plot_holes(ax=ax)
    trellis.plot_strong_pip(ax=ax, s=120, facecolors="none", linewidths=1.5, zorder=20)
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


def draw_combined(session: TangleSession, fp, dynamics: SymbolicDynamics, title: str, folder: Path) -> Path:
    """
    All seven panels on one figure.

    Top row: original trellis, minimal trellis, cartoon. Bottom row:
    partitions, the itinerary table, the two transition graphs.
    """
    view = plane_view(session, [fp], session.minimal_trellis())
    table_height = 0.45 * len(dynamics.classes) + 2.0
    fig = plt.figure(figsize=(36, 13 + max(13.0, table_height)))
    grid = fig.add_gridspec(2, 6, height_ratios=[13, max(13.0, table_height)])
    draw_original_trellis(fig.add_subplot(grid[0, 0:2]), session, fp, dynamics, view, title)
    draw_trellis(fig.add_subplot(grid[0, 2:4]), session, [fp], dynamics, view, title)
    draw_cartoon(fig.add_subplot(grid[0, 4:6]), session, dynamics, title)
    draw_partitions(fig.add_subplot(grid[1, 0:2]), session, dynamics, title)
    table_ax = fig.add_subplot(grid[1, 2:4])
    table_ax.set_axis_off()
    table_ax.set_title(f"{title} -- itineraries and words", fontsize=11)
    for column, refined, name in ((4, False, "class symbols"), (5, True, "refined symbols")):
        ax = fig.add_subplot(grid[1, column])
        draw_transition_graph(ax, dynamics, title, refined=refined)
        ax.set_title(f"transition graph ({name})", fontsize=11)
    fig.suptitle(title, fontsize=18)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    draw_itinerary_table(table_ax, dynamics)
    return save_panel(fig, folder, "combined")


def run_slug(
    unstable_steps: int, stable_steps: int, blasts: int, scale: bool,
    area_cutoff: float = AREA_CUTOFF,
) -> str:
    """``p3_k25_blast6``, ``p3_k25_u16_noscale``, ...: the folder name of a run."""
    slug = "p3_k25"
    if unstable_steps != 13:
        slug += f"_u{unstable_steps}"
    if stable_steps != 9:
        slug += f"_s{stable_steps}"
    if not scale:
        slug += "_noscale"
    if area_cutoff != AREA_CUTOFF:
        slug += f"_a{area_cutoff:.0e}".replace("e-0", "e-")
    if blasts:
        slug += f"_blast{blasts}"
    return slug


def main() -> None:
    parser = argparse.ArgumentParser(description="The whole pipeline on the Hénon k=2.5 period-3 orbit.")
    parser.add_argument(
        "--blasts", type=int, nargs="+", default=[12, 24],
        help="single-iteration zone blasts; one run (and folder) per count",
    )
    parser.add_argument("--unstable-steps", type=int, default=13, help="unstable steps of the k=2 recipe")
    parser.add_argument("--stable-steps", type=int, default=9, help="stable steps of the k=2 recipe")
    parser.add_argument("--no-scale", action="store_true", help="grow the step counts verbatim")
    parser.add_argument(
        "--area-cutoff", type=float, default=AREA_CUTOFF,
        help="refinement area cutoff (default 1e-7; finer resolves deeper pip-orbit points)",
    )
    args = parser.parse_args()
    configure_logging()

    scale = not args.no_scale
    for blasts in args.blasts:
        run(args.unstable_steps, args.stable_steps, blasts, scale, args.area_cutoff)


def closure_summary(dynamics: SymbolicDynamics) -> None:
    """
    Print, per class, where it comes from and how its word is backed.

    A class is either EVIDENCE (a computed bridge spans it) or CLOSURE (the
    word of another class names it, ``named_by``); a closure class no computed
    bridge spans yet is "unspanned". ``verified`` compares the word with the
    registered images of its members (None: no member image to compare).
    """
    print("closure summary (class: source, spanned, named by, verified, evidence):")
    for cd in dynamics.classes.values():
        named_by = dynamics.classes[cd.named_by].letter if cd.named_by is not None else "-"
        print(
            f"  {cd.letter:>3}: {cd.source or 'none':7} spanned={cd.entry.spanned!s:5} "
            f"named_by={named_by:>3} verified={cd.verified!s:5} evidence={len(cd.evidence)}"
        )
    unspanned = [cd.letter for cd in dynamics.classes.values() if not cd.entry.spanned]
    print(f"  {len(dynamics.classes)} classes, {len(unspanned)} unspanned: {unspanned or 'none'}")
    print()


def run(
    unstable_steps: int, stable_steps: int, blasts: int, scale: bool,
    area_cutoff: float = AREA_CUTOFF,
) -> Path:
    """Build one run, print its report, draw and save its panels; returns its folder."""
    session, fp, grown = build(
        unstable_steps=unstable_steps, stable_steps=stable_steps, blasts=blasts, scale=scale,
        area_cutoff=area_cutoff,
    )
    dynamics = session.symbolic_dynamics()
    blasted = f", {blasts} blast(s)" if blasts else ""
    title = (
        f"Hénon k={K:g}, b={B}, period 3 ({grown['unstable']}u/{grown['stable']}s steps, "
        f"cutoff {area_cutoff:g}){blasted}"
    )
    report(session, dynamics, title)
    closure_summary(dynamics)

    folder = OUT_DIR / run_slug(unstable_steps, stable_steps, blasts, scale, area_cutoff)
    draw_panels(session, [fp], title, folder, dynamics=dynamics)
    fig, ax = plt.subplots(figsize=(12, 10))
    draw_original_trellis(ax, session, fp, dynamics, plane_view(session, [fp], session.minimal_trellis()), title)
    fig.tight_layout()
    save_panel(fig, folder, "trellis_original")
    draw_combined(session, fp, dynamics, title, folder)
    print(f"wrote {folder}")
    return folder

if __name__ == "__main__":
    main()
