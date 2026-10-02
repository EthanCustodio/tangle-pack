"""
The symbolic-itinerary panels, one figure each, for any built tangle.

Shared by ``henon_symbolic_itineraries.py`` (the period-1 cases),
``henon_symbolic_itineraries_period3.py`` and
``henon_symbolic_itineraries_nested.py``. :func:`draw_panels` writes six
panels into one folder per case, each as ``<panel>.png``, ``<panel>.svg`` and
``<panel>.fig.pickle`` (reopen the pickle with ``scripts/show_figure.py`` to
pan and zoom it live):

1. ``trellis`` -- the MINIMAL trellis (its hole and image bridges, nothing
   else: blast children outside it are not drawn), each bridge coloured by
   its REFINED class (inert classes dashed), plus the holes and the chosen
   strong pips;
2. ``partitions`` -- the homotopy partition above the iterated homotopy
   partition as number lines, one row per (family, side, branch), each row
   labelled with its branch code (``A 1.0``) and every element with its name;
3. ``cartoon`` -- the dual graph as circular zones
   (``plot_dual_graph_cartoon(shape="circle")``): each fixed point a circle
   of ``2n`` arcs, stable arcs (labelled with their branch code) alternating
   with empty unstable ones, the zone-interior side of each branch inside
   the circle, nested zones as nested circles; bridges thin and black, no
   walk (the itineraries are the table's job);
4. ``itineraries`` -- the table: class, iterated itinerary (full element
   codes), word, refined word;
5. ``transitions`` -- the transition graph over the class symbols;
6. ``transitions_refined`` -- the transition graph over the refined symbols.

Dev Notes:
    The pickles need the same matplotlib version to load. Text keeps its
    point size when zooming, so zooming in on a crowded cartoon spreads the
    names apart -- that is what the pickles are for. ``draw_overview`` keeps
    the old 2x3 single-figure layout (line cartoon) for period-1 cases.
"""

from __future__ import annotations

import logging
import pickle
import sys
from collections import Counter
from pathlib import Path
from typing import Optional, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))

from henon_minimal_dual_graph import (  # noqa: E402
    _apply_plane_view,
    _draw_stable,
    _region_of_interest,
)

from tanglepack import TangleSession  # noqa: E402
from tanglepack.topology import plotting  # noqa: E402
from tanglepack.topology.SymbolicDynamics import SymbolicDynamics  # noqa: E402

logger = logging.getLogger("symbolic_figures")

FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"
#: The folder every case's panel folder goes in.
PANELS_DIR = FIGURES_DIR / "henon_symbolic_itineraries"

#: The panels, in drawing order.
PANELS = ("trellis", "partitions", "cartoon", "itineraries", "transitions", "transitions_refined")

#: Largest table font (points) on the itinerary panel.
TABLE_FONTSIZE = 22.0
#: Height of a table row in ems of its font.
TABLE_ROW_HEIGHT = 2.2
#: The table's columns (status and verified stay in the printed report).
TABLE_COLUMNS = ("class", "iterated itinerary", "word", "refined word")
#: Font size (points) of the cartoon's element names.
CARTOON_LABEL_FONTSIZE = 13.0
#: The cartoon's bridges: thin, black and opaque, all one width.
CARTOON_BRIDGE_COLOR = "black"
CARTOON_BRIDGE_ALPHA = 1.0
CARTOON_BRIDGE_LINEWIDTH = 1.2
#: Resolution of the PNGs.
DPI = 150


def save_panel(fig, folder: Path, name: str) -> Path:
    """
    Write one panel as PNG, SVG and a pickled figure, then close it.

    Args:
        fig: The panel's figure.
        folder: The case's folder (created if missing).
        name: The panel's file stem.

    Returns:
        The PNG's path.
    """
    folder.mkdir(parents=True, exist_ok=True)
    png = folder / f"{name}.png"
    fig.savefig(png, dpi=DPI, bbox_inches="tight")
    fig.savefig(folder / f"{name}.svg", bbox_inches="tight")
    _freeze_local_texts(fig)
    with open(folder / f"{name}.fig.pickle", "wb") as handle:
        pickle.dump(fig, handle)
    plt.close(fig)
    logger.info("wrote %s (+ .svg, .fig.pickle)", png)
    return png


def _freeze_local_texts(fig) -> None:
    """
    Replace texts of classes pickle cannot reach by plain copies.

    networkx draws curved edge labels as a ``Text`` subclass defined inside a
    function (``CurvedArrowText``); after the figure has been drawn once
    (the PNG) each one holds its final position, so a plain ``Text`` with the
    same content, placement and style replaces it.
    """
    from matplotlib.text import Text

    for ax in fig.axes:
        for text in list(ax.texts):
            if type(text) is Text or "<locals>" not in type(text).__qualname__:
                continue
            bbox = text.get_bbox_patch()
            copy = ax.text(
                *text.get_position(), text.get_text(),
                transform=text.get_transform(), rotation=text.get_rotation(),
                fontsize=text.get_fontsize(), color=text.get_color(),
                ha=text.get_horizontalalignment(), va=text.get_verticalalignment(),
                zorder=text.get_zorder(),
            )
            if bbox is not None:
                copy.set_bbox({
                    "boxstyle": bbox.get_boxstyle(), "facecolor": bbox.get_facecolor(),
                    "edgecolor": bbox.get_edgecolor(),
                })
            text.remove()


def _legend_below(ax, handles: list, ncol: int = 4) -> None:
    """A legend centred just below the axes."""
    ax.legend(
        handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.04),
        ncol=min(ncol, max(len(handles), 1)), fontsize=9, framealpha=0.9,
    )


def plane_view(session: TangleSession, fixed_points: Sequence, minimal) -> tuple:
    """The union of every fixed point's plane view (the folds near it)."""
    views = [_region_of_interest(session, fp, minimal) for fp in fixed_points]
    return (
        (min(v[0][0] for v in views), max(v[0][1] for v in views)),
        (min(v[1][0] for v in views), max(v[1][1] for v in views)),
    )


def _trellis(session: TangleSession, fixed_points: Sequence):
    """The trellis covering every fixed point (the single one's for one)."""
    return session.trellis(fixed_points[0]) if len(fixed_points) == 1 else session.trellis()


def draw_trellis(ax, session, fixed_points, dynamics, view, title) -> None:
    """Panel 1: the minimal trellis with every classed bridge coloured by refined class."""
    for fp in fixed_points:
        _draw_stable(ax, session, fp)
    minimal = session.minimal_trellis()
    plotting.plot_minimal_trellis(
        minimal, ax=ax, show_dropped=False, color="lightgray", linewidth=0.8, linestyle="-"
    )
    # Only the minimal trellis's bridges: every blast child is a class member
    # too, and drawing all members would add each blast's children here.
    drawn = plotting.plot_bridges_by_class(
        _trellis(session, fixed_points), dynamics, ax=ax, refined=True,
        bridge_ids=minimal.kept_bridge_ids,
    )
    for fp in fixed_points:
        session.plot_holes(fp, ax=ax)
        session.plot_strong_pip(fp, ax=ax, s=120, facecolors="none", linewidths=1.5, zorder=20)
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
        f"{title} -- minimal trellis ({len(minimal.kept_bridge_ids)} bridges: "
        f"{len(minimal.hole_bridge_ids)} hole, {len(minimal.image_bridge_ids)} image), "
        "by refined class",
        fontsize=11,
    )


def draw_partitions(ax, session, dynamics, title) -> None:
    """Panel 2: homotopy rows above iterated rows, rows by branch code, elements by name."""
    homotopy = session.homotopy_partition()
    iterated = session.iterated_partition()
    homotopy_results = homotopy.as_list()
    results = homotopy_results + iterated.as_list()
    naming = dynamics.naming
    labels = [
        f"{family} {result.side} ({naming.branch_code(result.branch_key)})"
        for family, family_results in (("homotopy", homotopy), ("iterated", iterated))
        for result in family_results
    ]

    def element_label(result, interval) -> str:
        ref = result.ref(interval.element_id)
        if any(result is candidate for candidate in homotopy_results):
            return naming.homotopy_name(ref).mathtext
        return naming.name(ref).mathtext

    plotting.plot_stable_partition(results, ax=ax, labels=labels, element_labels=element_label)
    ax.set_title(f"{title} -- homotopy vs iterated homotopy partition (element names)", fontsize=11)


def draw_cartoon(ax, session, dynamics, title, *, shape: str = "circle") -> None:
    """Panel 3: the dual graph as circular zones (or rows), thin black bridges, no walk."""
    layout = session.plot_dual_graph_cartoon(
        dynamics, ax=ax, shape=shape, show_walks=False, bridge_color=CARTOON_BRIDGE_COLOR,
        bridge_alpha=CARTOON_BRIDGE_ALPHA, bridge_linewidth=CARTOON_BRIDGE_LINEWIDTH,
        label_fontsize=CARTOON_LABEL_FONTSIZE,
    )
    handles = plotting.dual_graph_cartoon_legend_handles(
        bridge_color=CARTOON_BRIDGE_COLOR, bridge_alpha=CARTOON_BRIDGE_ALPHA,
        bridge_linewidth=CARTOON_BRIDGE_LINEWIDTH, walks=False, empty_side=shape == "circle",
    )
    _legend_below(ax, handles, ncol=3)
    if layout.bridges_drawn == 0:
        ax.text(
            0.02, 0.02, "no bridge to draw", transform=ax.transAxes, ha="left", va="bottom",
            fontsize=9, zorder=30,
            bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "#898781", "alpha": 0.9},
        )
    ax.set_title(f"{title} -- dual graph cartoon: elements and bridges", fontsize=11)


def draw_itinerary_table(ax, dynamics: SymbolicDynamics):
    """Panel 4's table, drawn AFTER the layout has settled the axes width."""
    return plotting.plot_itinerary_table(
        dynamics, ax=ax, refined=True, columns=TABLE_COLUMNS,
        fontsize=TABLE_FONTSIZE, row_height=TABLE_ROW_HEIGHT,
    )


def draw_transition_graph(ax, dynamics: SymbolicDynamics, title: str, *, refined: bool) -> None:
    """Panels 5 and 6: the transition graph over class or refined symbols."""
    plotting.plot_transition_graph(dynamics, ax=ax, refined=refined)
    ax.set_title(f"{title} -- {ax.get_title()}", fontsize=11)


def report(session: TangleSession, dynamics: SymbolicDynamics, title: str) -> None:
    """
    Print one case's symbolic dynamics report and a short health summary:
    how each class was resolved, the unresolved ones, the minimal trellis.
    """
    print(title)
    print("=" * len(title))
    print(dynamics.describe())
    minimal = session.minimal_trellis()
    sources = Counter(
        f"{cd.source or 'none'}/{cd.search.status if cd.search is not None else '-'}"
        for cd in dynamics.classes.values()
    )
    unresolved = [cd.letter for cd in dynamics.classes.values() if cd.itinerary is None]
    print()
    print("health:")
    print(f"  classes: {len(dynamics.classes)}; resolved by source/status: {dict(sources)}")
    print(f"  unresolved: {unresolved or 'none'}; reliable: {dynamics.is_reliable}")
    print(
        f"  minimal trellis: {len(minimal.hole_bridge_ids)} hole bridge(s), "
        f"{len(minimal.image_bridge_ids)} image bridge(s), "
        f"{len(minimal.skipped_pairs)} skipped pair(s)"
    )
    print()


def draw_panels(
    session: TangleSession,
    fixed_points: Sequence,
    title: str,
    folder: Path,
    *,
    dynamics: Optional[SymbolicDynamics] = None,
) -> list[Path]:
    """
    Draw and save the six panels of one case.

    Args:
        session: The built, partitioned session.
        fixed_points: Its fixed points, outermost first.
        title: The case's title (prefixes every panel title).
        folder: The case's folder.
        dynamics: The symbolic dynamics; ``session.symbolic_dynamics()`` by
            default.

    Returns:
        The PNG paths, in :data:`PANELS` order.
    """
    dynamics = dynamics if dynamics is not None else session.symbolic_dynamics()
    minimal = session.minimal_trellis()
    view = plane_view(session, fixed_points, minimal)
    naming = dynamics.naming
    branches = len(naming.homotopy.branch_keys)
    paths: list[Path] = []

    fig, ax = plt.subplots(figsize=(12, 10))
    draw_trellis(ax, session, fixed_points, dynamics, view, title)
    fig.tight_layout()
    paths.append(save_panel(fig, folder, "trellis"))

    rows = 2 * len(session.homotopy_partition().as_list())
    fig, ax = plt.subplots(figsize=(14, max(6.0, 0.9 * rows + 2.0)))
    draw_partitions(ax, session, dynamics, title)
    fig.tight_layout()
    paths.append(save_panel(fig, folder, "partitions"))

    side = 10.0 + 1.5 * branches
    fig, ax = plt.subplots(figsize=(side, side + 1.0))
    draw_cartoon(ax, session, dynamics, title)
    fig.tight_layout()
    paths.append(save_panel(fig, folder, "cartoon"))

    fig, ax = plt.subplots(figsize=(18, 0.75 * len(dynamics.classes) + 1.5))
    ax.set_axis_off()
    ax.set_title(f"{title} -- itineraries and words", fontsize=11)
    fig.tight_layout()
    draw_itinerary_table(ax, dynamics)
    paths.append(save_panel(fig, folder, "itineraries"))

    for refined, name in ((False, "transitions"), (True, "transitions_refined")):
        fig, ax = plt.subplots(figsize=(9, 9))
        draw_transition_graph(ax, dynamics, title, refined=refined)
        fig.tight_layout()
        paths.append(save_panel(fig, folder, name))
    return paths


def draw_overview(
    session: TangleSession, fixed_points: Sequence, title: str, path: Path,
    *, dynamics: Optional[SymbolicDynamics] = None,
) -> Path:
    """
    The old single 2x3 figure (line cartoon), for period-1 cases.

    Args:
        session: The built, partitioned session.
        fixed_points: Its fixed points.
        title: The case's title.
        path: Where to write the PNG.
        dynamics: The symbolic dynamics; ``session.symbolic_dynamics()`` by
            default.

    Returns:
        ``path``.
    """
    dynamics = dynamics if dynamics is not None else session.symbolic_dynamics()
    view = plane_view(session, fixed_points, session.minimal_trellis())
    fig, axes = plt.subplots(2, 3, figsize=(27, 15), gridspec_kw={"width_ratios": [1.3, 1.0, 1.3]})
    draw_trellis(axes[0, 0], session, fixed_points, dynamics, view, title)
    draw_partitions(axes[0, 1], session, dynamics, title)
    draw_cartoon(axes[0, 2], session, dynamics, title, shape="line")
    axes[1, 0].set_axis_off()
    axes[1, 0].set_title(f"{title} -- itineraries and words", fontsize=11)
    draw_transition_graph(axes[1, 1], dynamics, title, refined=False)
    draw_transition_graph(axes[1, 2], dynamics, title, refined=True)
    fig.tight_layout()
    draw_itinerary_table(axes[1, 0], dynamics)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    logger.info("wrote %s", path)
    return path


def configure_logging() -> None:
    """INFO for the scripts, WARNING for the library's chatter."""
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(name)s: %(message)s")
    logger.setLevel(logging.INFO)


__all__ = [
    "PANELS",
    "PANELS_DIR",
    "configure_logging",
    "draw_overview",
    "draw_panels",
    "report",
    "save_panel",
]
