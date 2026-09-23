"""
Symbolic itineraries of the bridge classes on the Hénon cases of
``henon_bridge_classes.py``: the words the dual-graph walks read off.

One figure per case (k=10; k=2.8 blasted once; k=2.8 blasted twice), six
panels in a 2x3 grid:

1. the tangle with the MINIMAL trellis underneath (light gray) and every
   classed bridge coloured by its REFINED class (``a_1``, ``a_2``, inert
   classes dashed), plus the holes and the chosen strong pip;
2. the homotopy partition above the iterated homotopy partition as number
   lines, every element labelled with its name (``L_1``, ``R_3^{2}``);
3. the dual graph as a CARTOON (``plot_dual_graph_cartoon``): each stable
   branch a line, anchor on the right, with its iterated elements as
   bracketed bars in an ordinal coordinate (left above, right below), one
   node per element (open = wall, filled = unified) named just to its
   right, the unified crossings as S-curves through the line and the
   minimal trellis's bridges as nested thin black U-arcs (inert classes
   dashed); no walk is drawn (the itineraries are the table's job);
4. the itinerary table: class, itinerary in element names, word and refined
   word (the status and verified columns are left to the printed report);
5. the transition graph over the class symbols;
6. the transition graph over the refined symbols.

The cases are imported from ``henon_bridge_classes.py`` (same recipes, same
``build``) and the plane panels are framed as in
``henon_minimal_dual_graph.py`` (the folds near the fixed point).

Run with ``PYTHONPATH=src python3 scripts/henon_symbolic_itineraries.py``;
the figures land in ``figures/henon_symbolic_itineraries_<case>.png``
(``k10``, ``k28_blast1``, ``k28_blast2``) and each case's symbolic dynamics
report is printed.

Dev Notes:
    The one-blast k=2.8 case has its single active class resolved through
    the singleton/trellis path (both of its elements are singletons and own
    no dual-graph node); the printed report shows ``source = trellis`` for
    it. The cartoon panel's text box appears only when no bridge was drawn.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))

from henon_bridge_classes import CASES, Case, build  # noqa: E402
from henon_minimal_dual_graph import (  # noqa: E402
    _apply_plane_view,
    _draw_stable,
    _region_of_interest,
)

from tanglepack import TangleSession  # noqa: E402
from tanglepack.topology import plotting  # noqa: E402
from tanglepack.topology.SymbolicDynamics import SymbolicDynamics  # noqa: E402

logger = logging.getLogger("henon_symbolic_itineraries")

FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"

#: Largest table font (points) on the itinerary panel; the plotter's own cap
#: (8 pt) is sized for a small axes, this panel is a third of a 27 in figure.
TABLE_FONTSIZE = 22.0
#: Height of a table row in ems of its font.
TABLE_ROW_HEIGHT = 2.2
#: The table's columns: the status and verified flags stay in the printed
#: report, which leaves the panel's width to the itineraries and words.
TABLE_COLUMNS = ("class", "itinerary", "word", "refined word")
#: Font size (points) of the cartoon's element names: between the table's
#: font and the plotter's default.
CARTOON_LABEL_FONTSIZE = 16.0
#: The cartoon's bridges: thin, black and opaque, all one width (no walk is
#: drawn over them, and the hole / image distinction is not this figure's
#: point).
CARTOON_BRIDGE_COLOR = "black"
CARTOON_BRIDGE_ALPHA = 1.0
CARTOON_BRIDGE_LINEWIDTH = 1.2


def case_slug(case: Case) -> str:
    """``k10``, ``k28_blast1``, ...: the file-name stem of a case."""
    k_text = f"{case.k:g}".replace(".", "")
    blasted = f"_blast{case.blasts}" if case.blasts else ""
    return f"k{k_text}{blasted}"


def out_path(case: Case) -> Path:
    """Where a case's figure is written."""
    return FIGURES_DIR / f"henon_symbolic_itineraries_{case_slug(case)}.png"


def _legend_below(ax, handles: list) -> None:
    """
    A plane panel's legend, centred just BELOW the axes.

    The plane panels keep the plane's aspect, which leaves room under the
    box; the folds' ends, the pip and the walks' nodes all sit in the upper
    right corner, where an inside legend would cover them.
    """
    ax.legend(
        handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.06),
        ncol=min(4, max(len(handles), 1)), fontsize=9, framealpha=0.9,
    )


def draw_bridges_by_class(
    ax, session: TangleSession, fp, dynamics: SymbolicDynamics, view: tuple, title: str
) -> None:
    """Panel 1: the minimal trellis with every classed bridge coloured by refined class."""
    _draw_stable(ax, session, fp)
    minimal = session.minimal_trellis()
    plotting.plot_minimal_trellis(
        minimal, ax=ax, show_dropped=False, color="lightgray", linewidth=0.8, linestyle="-"
    )
    plotting.plot_bridges_by_class(session.trellis(fp), dynamics, ax=ax, refined=True)
    session.plot_holes(fp, ax=ax)
    session.plot_strong_pip(
        fp, ax=ax, s=120, facecolors="none", linewidths=1.5, zorder=20
    )
    # plot_bridges_by_class placed its own legend inside; move it below.
    legend = ax.get_legend()
    if legend is not None:
        legend.remove()
    _legend_below(ax, plotting.bridge_class_legend_handles(dynamics, refined=True))
    _apply_plane_view(ax, view)
    ax.set_title(f"{title} -- minimal trellis, bridges by refined class", fontsize=11)


def draw_partitions(ax, session: TangleSession, dynamics: SymbolicDynamics, title: str) -> None:
    """Panel 2: homotopy rows above iterated rows, elements labelled by name."""
    homotopy = session.homotopy_partition()
    iterated = session.iterated_partition()
    homotopy_results = homotopy.as_list()
    results = homotopy_results + iterated.as_list()
    labels = [
        f"{family} {result.side} (orbit {result.branch_key[2]}.{result.branch_key[3]})"
        for family, family_results in (("homotopy", homotopy), ("iterated", iterated))
        for result in family_results
    ]
    naming = dynamics.naming

    def element_label(result, interval) -> str:
        ref = result.ref(interval.element_id)
        if any(result is candidate for candidate in homotopy_results):
            return naming.homotopy_name(ref).mathtext
        return naming.name(ref).mathtext

    plotting.plot_stable_partition(results, ax=ax, labels=labels, element_labels=element_label)
    ax.set_title(f"{title} -- homotopy vs iterated homotopy partition (element names)", fontsize=11)


def draw_dual_graph_cartoon(ax, session: TangleSession, dynamics: SymbolicDynamics, title: str) -> None:
    """Panel 3: the dual graph as a cartoon, names beside the nodes, thin black bridges, no walk."""
    layout = session.plot_dual_graph_cartoon(
        dynamics, ax=ax, show_walks=False, bridge_color=CARTOON_BRIDGE_COLOR,
        bridge_alpha=CARTOON_BRIDGE_ALPHA, bridge_linewidth=CARTOON_BRIDGE_LINEWIDTH,
        label_fontsize=CARTOON_LABEL_FONTSIZE,
    )
    handles = plotting.dual_graph_cartoon_legend_handles(
        bridge_color=CARTOON_BRIDGE_COLOR, bridge_alpha=CARTOON_BRIDGE_ALPHA,
        bridge_linewidth=CARTOON_BRIDGE_LINEWIDTH, walks=False,
    )
    ax.legend(
        handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.02),
        ncol=3, fontsize=9, framealpha=0.9,
    )
    if layout.bridges_drawn == 0:
        ax.text(
            0.02, 0.02, "no bridge to draw",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=9,
            bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "#898781", "alpha": 0.9},
            zorder=30,
        )
    ax.set_title(f"{title} -- dual graph cartoon: elements and bridges", fontsize=11)


def draw_itinerary_table(ax, dynamics: SymbolicDynamics, title: str):
    """
    Panel 4: the itinerary table.

    Called AFTER the figure's layout has settled the axes width: the plotter
    sizes its columns and font to fill that width (the axes is off, so the
    layout needs nothing but the title, set in :func:`title_itinerary_table`).

    Returns:
        The matplotlib ``Table`` (or None).
    """
    return plotting.plot_itinerary_table(
        dynamics, ax=ax, refined=True, columns=TABLE_COLUMNS,
        fontsize=TABLE_FONTSIZE, row_height=TABLE_ROW_HEIGHT,
    )


def title_itinerary_table(ax, title: str) -> None:
    """Panel 4's title, set before the layout so the layout leaves room for it."""
    ax.set_axis_off()
    ax.set_title(f"{title} -- itineraries and words", fontsize=11)


def draw_transition_graph(ax, dynamics: SymbolicDynamics, title: str, *, refined: bool) -> None:
    """Panels 5 and 6: the transition graph over class or refined symbols."""
    plotting.plot_transition_graph(dynamics, ax=ax, refined=refined)
    ax.set_title(f"{title} -- {ax.get_title()}", fontsize=11)


def report(dynamics: SymbolicDynamics, title: str) -> None:
    """Print one case's symbolic dynamics report."""
    print(title)
    print("=" * len(title))
    print(dynamics.describe())
    print()


def draw_case(case: Case) -> Path:
    """Build one case, print its report, draw and save its figure."""
    session, fp, _pre_trim_stable_points = build(case)
    dynamics = session.symbolic_dynamics()
    minimal = session.minimal_trellis()
    view = _region_of_interest(session, fp, minimal)
    report(dynamics, case.title)

    fig, axes = plt.subplots(
        2, 3, figsize=(27, 15), gridspec_kw={"width_ratios": [1.3, 1.0, 1.3]}
    )
    draw_bridges_by_class(axes[0, 0], session, fp, dynamics, view, case.title)
    draw_partitions(axes[0, 1], session, dynamics, case.title)
    draw_dual_graph_cartoon(axes[0, 2], session, dynamics, case.title)
    title_itinerary_table(axes[1, 0], case.title)
    draw_transition_graph(axes[1, 1], dynamics, case.title, refined=False)
    draw_transition_graph(axes[1, 2], dynamics, case.title, refined=True)

    fig.tight_layout()
    draw_itinerary_table(axes[1, 0], dynamics, case.title)
    FIGURES_DIR.mkdir(exist_ok=True)
    path = out_path(case)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    logger.info("wrote %s", path)
    return path


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    for case in CASES:
        draw_case(case)


if __name__ == "__main__":
    main()
