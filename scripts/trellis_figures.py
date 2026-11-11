"""
The trellis figure package: a computed tangle and its trellis cartoon.

Writes four panels and a combined 2x2 figure of one case into
``figures/trellis_figures/<slug>/``, each as ``.png``, ``.svg`` and
``.fig.pickle`` (:func:`symbolic_figures.save_panel`):

1. ``tangle`` -- every computed manifold of every fixed point as fine lines
   (unstable blue, stable red), blast curves included;
2. ``pip_pseudoneighbors`` -- the trellis bridges in light grey over the
   stable manifold, the strong-pip candidates (small), the chosen pip and its
   orbit (large rings, ``q0, q1, ...``) and the pseudoneighbor pairs
   (numbered dot pairs joined by a thin segment);
3. ``minimal_trellis`` -- the minimal trellis's bridges by class, its holes,
   the pip orbit and the stable manifold;
4. ``cartoon`` -- the trellis cartoon
   (:func:`tanglepack.loom.plot_trellis_cartoon`): one representative per
   class, holes on;
5. ``combined`` -- the four on one figure.

Run with ``PYTHONPATH=src python3 scripts/trellis_figures.py --case k28
--blasts 2`` (cases ``k10``, ``k28``, ``p3``, ``nested``, ``inversion``);
``--show all``, ``--letters a b`` and ``--iterates 1`` pass on to the cartoon.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))

from symbolic_figures import (  # noqa: E402
    FIGURES_DIR,
    _legend_below,
    configure_logging,
    plane_view,
    save_panel,
)
from henon_minimal_dual_graph import _apply_plane_view, _draw_stable  # noqa: E402

from tanglepack import TangleSession  # noqa: E402
from tanglepack.loom import plot_trellis_cartoon  # noqa: E402
from tanglepack.topology import plotting  # noqa: E402
from tanglepack.topology.TrellisCartoon import pip_orbit  # noqa: E402

logger = logging.getLogger("trellis_figures")

#: The folder every case's figure folder goes in.
TRELLIS_DIR = FIGURES_DIR / "trellis_figures"

#: The panels, in drawing order (the combined figure is written after them).
PANELS = ("tangle", "pip_pseudoneighbors", "minimal_trellis", "cartoon")

#: Line width of the manifolds on the tangle panel.
TANGLE_LINEWIDTH = 0.4

#: Colours of the numbered pseudoneighbor pairs, cycled.
PAIR_COLORS = plt.get_cmap("tab10").colors


def _draw_pip_orbit(ax, trellis, *, labels: bool) -> None:
    """The chosen pips and their orbits as large open rings, ``q_j`` beside each."""
    for j, crossing in pip_orbit(trellis):
        x, y = trellis.intersection(crossing).get_point()
        ax.scatter(x, y, s=150, facecolors="none", edgecolors="green", linewidths=1.5, zorder=20)
        if labels:
            ax.annotate(f"$q_{{{j}}}$", (x, y), textcoords="offset points", xytext=(7, 7),
                        fontsize=11, color="green", zorder=21)


def draw_tangle(ax, session: TangleSession, topology, view: tuple, title: str) -> None:
    """Panel 1: every computed manifold, fine lines."""
    plt.sca(ax)  # the workbench plotter draws on the current axes
    for fp in session.fixed_points:
        session.plot_tangle(fp, "unstable", color="tab:blue", linewidth=TANGLE_LINEWIDTH)
        session.plot_tangle(fp, "stable", color="tab:red", linewidth=TANGLE_LINEWIDTH)
    for bridge in topology.trellis.bridges:  # blast image curves are bridges, not branches
        points = bridge.get_point_array()
        ax.plot(points[:, 0], points[:, 1], color="tab:blue", linewidth=TANGLE_LINEWIDTH)
    _apply_plane_view(ax, view)
    ax.set_title(f"{title} -- the computed tangle", fontsize=11)


def draw_pip_pseudoneighbors(ax, session: TangleSession, topology, view: tuple, title: str) -> None:
    """Panel 2: the trellis, the strong pip and its orbit, the pseudoneighbor pairs."""
    trellis = topology.trellis
    for fp in topology.fixed_points:
        _draw_stable(ax, session, fp)
        session.trellis(fp).plot_strong_pip_candidates(ax=ax, s=6)
    for bridge in trellis.bridges:
        points = bridge.get_point_array()
        ax.plot(points[:, 0], points[:, 1], color="lightgray", linewidth=0.8, zorder=1)
    _draw_pip_orbit(ax, trellis, labels=True)
    for number, pair in enumerate(trellis.pseudoneighbors, start=1):
        color = PAIR_COLORS[(number - 1) % len(PAIR_COLORS)]
        a, b = (trellis.intersection(i).get_point() for i in (pair.intersection_a, pair.intersection_b))
        ax.plot([a[0], b[0]], [a[1], b[1]], color=color, linewidth=0.8, zorder=13)
        ax.scatter([a[0], b[0]], [a[1], b[1]], color=color, s=14, zorder=14)
        ax.annotate(str(number), tuple(a), textcoords="offset points", xytext=(4, -10),
                    fontsize=9, color=color, zorder=15)
    _apply_plane_view(ax, view)
    ax.set_title(
        f"{title} -- strong pip, its orbit and {len(trellis.pseudoneighbors)} pseudoneighbor pair(s)",
        fontsize=11,
    )


def draw_minimal_trellis(ax, session: TangleSession, topology, view: tuple, title: str) -> None:
    """Panel 3: the minimal trellis by class, its holes and the pip orbit."""
    minimal, dynamics, trellis = topology.minimal, topology.dynamics, topology.trellis
    for fp in topology.fixed_points:
        _draw_stable(ax, session, fp)
    drawn = plotting.plot_bridges_by_class(
        trellis, dynamics, ax=ax, refined=False, bridge_ids=minimal.kept_bridge_ids,
    )
    plotting.plot_holes(trellis, ax=ax)
    _draw_pip_orbit(ax, trellis, labels=False)
    legend = ax.get_legend()
    if legend is not None:
        legend.remove()
    shown = {plotting.name_mathtext(name) for name in drawn}
    _legend_below(ax, [
        handle for handle in plotting.bridge_class_legend_handles(dynamics, refined=False)
        if handle.get_label().split(" (")[0] in shown
    ])
    _apply_plane_view(ax, view)
    ax.set_title(
        f"{title} -- minimal trellis ({len(minimal.hole_bridge_ids)} hole, "
        f"{len(minimal.image_bridge_ids)} image bridges), holes and classes",
        fontsize=11,
    )


def draw_cartoon(ax, session: TangleSession, topology, title: str, cartoon_kw: dict) -> None:
    """Panel 4: the trellis cartoon."""
    kwargs = {"holes": True, **cartoon_kw}
    layout = plot_trellis_cartoon(
        topology, session.resonance_zones.values(), session.workbench, ax=ax, **kwargs
    )
    ax.set_title(f"{title} -- trellis cartoon ({len(layout.bridges)} bridges)", fontsize=11)


def write_trellis_figures(
    session: TangleSession, folder: Path, title: str, *, cartoon_kw: Optional[dict] = None
) -> Path:
    """
    Draw and save the four panels and the combined figure of one built tangle.

    Args:
        session: The built, partitioned session (its zones lay the cartoon out).
        folder: Where the files go (created if missing).
        title: The case's title (prefixes every panel title).
        cartoon_kw: Passed on to the cartoon (``show``, ``letters``,
            ``iterates``, ...).

    Returns:
        ``folder``.
    """
    cartoon_kw = cartoon_kw or {}
    topology = session.topology()
    view = plane_view(session, topology.fixed_points, topology.minimal)
    plane = (draw_tangle, draw_pip_pseudoneighbors, draw_minimal_trellis)

    for name, draw in zip(PANELS, plane):
        fig, ax = plt.subplots(figsize=(11, 10))
        draw(ax, session, topology, view, title)
        fig.tight_layout()
        save_panel(fig, folder, name)
    fig, ax = plt.subplots(figsize=(11, 11))
    draw_cartoon(ax, session, topology, title, cartoon_kw)
    fig.tight_layout()
    save_panel(fig, folder, "cartoon")

    fig, axes = plt.subplots(2, 2, figsize=(20, 19))
    for ax, draw in zip(axes.flat, plane):
        draw(ax, session, topology, view, "")
    draw_cartoon(axes[1, 1], session, topology, "", cartoon_kw)
    for ax in axes.flat:
        ax.set_title(ax.get_title().removeprefix(" -- "), fontsize=12)
    fig.suptitle(title, fontsize=16)
    fig.tight_layout()
    save_panel(fig, folder, "combined")
    return folder


# ── the cases ───────────────────────────────────────────────────────────────


def build_case(case: str, blasts: Optional[int]) -> tuple[TangleSession, str, str]:
    """``(session, title, slug)`` of one named case, ``blasts`` overriding its default."""
    from tanglepack.examples import henon_cases

    if case in ("k10", "k28"):
        from henon_bridge_classes import Case, build

        k, cutoff, default = (10, 1e-4, 0) if case == "k10" else (2.8, 1e-7, 1)
        recipe = Case(k=k, unstable_steps=10, pip="image", area_cutoff=cutoff,
                      blasts=default if blasts is None else blasts)
        session, _fp, _points = build(recipe)
        slug = case + (f"_blast{recipe.blasts}" if recipe.blasts else "")
        return session, recipe.title, slug
    if case == "p3":
        built = henon_cases.build_period3(blasts=blasts or 0)
    elif case == "nested":
        built = henon_cases.build_nested(**({} if blasts is None else {"outer_blasts": blasts}))
    elif case == "inversion":
        built = henon_cases.build_inversion(blasts=blasts or 0)
    else:
        raise ValueError(f"unknown case {case!r}")
    return built.session, built.title, built.slug


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--case", choices=("k10", "k28", "p3", "nested", "inversion"), default="k28")
    parser.add_argument("--blasts", type=int, default=None, help="blasts (the case's default if unset)")
    parser.add_argument("--show", choices=("representatives", "all"), default="representatives")
    parser.add_argument("--letters", nargs="+", default=None, help="draw only these classes")
    parser.add_argument("--iterates", type=int, default=0, help="forward iterates to draw")
    args = parser.parse_args()
    configure_logging()
    session, title, slug = build_case(args.case, args.blasts)
    cartoon_kw = {"show": args.show, "letters": args.letters, "iterates": args.iterates}
    if args.show != "representatives":
        slug += f"_{args.show}"
    if args.letters:
        slug += "_" + "".join(args.letters)
    if args.iterates:
        slug += f"_it{args.iterates}"
    folder = write_trellis_figures(session, TRELLIS_DIR / slug, title, cartoon_kw=cartoon_kw)
    logger.info("wrote %s", folder)


if __name__ == "__main__":
    main()
