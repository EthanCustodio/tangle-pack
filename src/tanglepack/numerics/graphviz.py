"""Matplotlib / networkx rendering of the intersection graph.

The workbench owns the *data* -- the registry and the bridges -- and this module
owns the *picture*. Both entry points are kept out of
:mod:`~.TangleWorkbench` so the orchestrator carries no graph-plotting code;
:meth:`~.TangleWorkbench.TangleWorkbench.build_intersection_graph` and
:meth:`~.TangleWorkbench.TangleWorkbench.visualize_intersection_graph` remain as
thin forwarders, so no caller changes.

Dev Notes:

:func:`visualize_intersection_graph` deliberately takes a graph rather than a
workbench: the graph a caller wants drawn is often a decorated or filtered COPY of
:func:`build_intersection_graph`'s output, and nothing here needs to read the
workbench back.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Optional

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D
import networkx as nx

if TYPE_CHECKING:  # pragma: no cover - import cycle guard
    from .TangleWorkbench import TangleWorkbench

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def build_intersection_graph(workbench: "TangleWorkbench") -> nx.MultiDiGraph:
    """
    The registry's intersection graph, decorated for this workbench.

    The graph itself -- nodes, per-branch adjacency, iterate edges -- belongs
    to the :class:`IntersectionRegistry` (plan 2.7); this method supplies the
    one thing the registry cannot know, the bridges, and returns a COPY so
    callers may decorate or lay out the result without touching the live one.

    Unstable adjacency is exactly the registered bridges: a bridge IS the
    piece of unstable manifold between two consecutive crossings, so
    disconnected pieces of unstable manifold (an original bridge and its
    iterated children) stay separate paths instead of merging into one chain.
    Stable adjacency joins consecutive crossings ON EACH STABLE BRANCH;
    canonical distances on different branches are measured from different
    anchors and are not comparable.

    Returns:
        nx.MultiDiGraph: A copy of the registry graph.
    """
    bridges = list(workbench._bridges)
    return workbench._intersection_registry.graph(bridges=bridges).copy()


def visualize_intersection_graph(
    G: nx.MultiDiGraph,
    layout: str = "auto",
    figsize: tuple[int, int] = (12, 8),
    display_mode: str = "auto",
    compact_threshold: int = 20,
    node_size: Optional[int] = None,
    label_mode: str = "id",
    node_color_by: str = "none",
    show_iterate_edges: bool = True,
    save_path: Optional[str] = None,
) -> tuple:
    """
    Visualize the intersection graph with edges colored by type and stability.

    Args:
        G: The intersection graph from build_intersection_graph().
        layout: Layout algorithm. ``"auto"`` (default) picks
            ``"kamada_kawai"`` for ≤ 8 nodes and ``"stable_linear"``
            otherwise. Explicit options:

            * ``"stable_linear"`` — nodes sorted by stable arc-length on a
              horizontal line; bridges arch above and stable edges run flat.
              Usually the clearest layout for tangle graphs.
            * ``"unstable_linear"`` — same but sorted by unstable arc-length.
            * ``"cdist"`` — node position = (unstable_cdist, stable_cdist),
              revealing the full arc-length structure on labelled axes.
            * ``"kamada_kawai"``, ``"spring"``, ``"circular"``, ``"spectral"``
              — standard networkx force/geometric layouts.
        figsize: Figure size as (width, height).
        display_mode: One of ``"auto"``, ``"full"``, or ``"compact"``.
            ``"auto"`` switches to compact when the node count exceeds
            *compact_threshold*. ``"full"`` uses large, labeled nodes.
            ``"compact"`` uses small dots suitable for dense graphs.
        compact_threshold: Node count above which ``"auto"`` picks compact.
        node_size: Override the node area in points². Defaults to 800 in
            full mode and 80 in compact mode.
        label_mode: Content of node labels. One of:
            ``"id"`` — intersection ID only (default);
            ``"coords"`` — (x, y) phase-space coordinates;
            ``"cdist"`` — unstable and stable arc-lengths;
            ``"all"`` — ID, coordinates, and both cdists;
            ``"none"`` — no labels.
        node_color_by: How to color the nodes. One of:
            ``"none"`` — white (default);
            ``"unstable_cdist"`` — viridis colormap by unstable arc-length;
            ``"stable_cdist"`` — plasma colormap by stable arc-length;
            ``"fixed_point"`` — distinct color per originating fixed point.
        show_iterate_edges: Whether to draw iterate-type edges (rendered
            dashed in purple to distinguish them from adjacency edges).
        save_path: Optional file path to save the figure.

    Returns:
        (fig, ax) matplotlib Figure and Axes.
    """
    if G.number_of_nodes() == 0:
        logger.warning("Graph has no nodes to visualize")
        return None, None

    n_nodes = G.number_of_nodes()

    mode = (
        ("compact" if n_nodes > compact_threshold else "full")
        if display_mode == "auto"
        else display_mode
    )

    if node_size is None:
        node_size = 80 if mode == "compact" else 800

    # ── Layout ────────────────────────────────────────────────────────────
    if layout == "auto":
        layout = "kamada_kawai" if n_nodes <= 8 else "stable_linear"

    if layout == "stable_linear":
        sorted_nodes = sorted(
            G.nodes(), key=lambda nd: G.nodes[nd].get("stable_cdist") or 0.0
        )
        pos = {nd: (float(i), 0.0) for i, nd in enumerate(sorted_nodes)}
    elif layout == "unstable_linear":
        sorted_nodes = sorted(
            G.nodes(), key=lambda nd: G.nodes[nd].get("unstable_cdist") or 0.0
        )
        pos = {nd: (float(i), 0.0) for i, nd in enumerate(sorted_nodes)}
    elif layout == "cdist":
        pos = {
            node: (
                G.nodes[node].get("unstable_cdist") or 0.0,
                G.nodes[node].get("stable_cdist") or 0.0,
            )
            for node in G.nodes()
        }
    elif layout == "spring":
        pos = nx.spring_layout(G, k=1, iterations=50, seed=42)
    elif layout == "circular":
        pos = nx.circular_layout(G)
    elif layout == "kamada_kawai":
        pos = nx.kamada_kawai_layout(G)
    elif layout == "spectral":
        pos = nx.spectral_layout(G)
    else:
        pos = nx.spring_layout(G, seed=42)

    fig, ax = plt.subplots(figsize=figsize)

    # ── Node colors ───────────────────────────────────────────────────────
    cmap_obj = None
    norm_obj = None
    if node_color_by == "unstable_cdist":
        values = [G.nodes[nd].get("unstable_cdist") or 0.0 for nd in G.nodes()]
        cmap_obj = cm.viridis
        norm_obj = mcolors.Normalize(vmin=min(values), vmax=max(values))
        node_colors = [cmap_obj(norm_obj(v)) for v in values]
    elif node_color_by == "stable_cdist":
        values = [G.nodes[nd].get("stable_cdist") or 0.0 for nd in G.nodes()]
        cmap_obj = cm.plasma
        norm_obj = mcolors.Normalize(vmin=min(values), vmax=max(values))
        node_colors = [cmap_obj(norm_obj(v)) for v in values]
    elif node_color_by == "fixed_point":
        fps = list(
            dict.fromkeys(
                G.nodes[nd].get("manifold_a_key", (None,))[0] for nd in G.nodes()
            )
        )
        fp_idx = {fp: i for i, fp in enumerate(fps)}
        fp_cmap = cm.get_cmap("Set1", max(len(fps), 1))
        node_colors = [
            fp_cmap(fp_idx.get(G.nodes[nd].get("manifold_a_key", (None,))[0], 0))
            for nd in G.nodes()
        ]
    else:
        node_colors = ["white"] * n_nodes

    # ── Draw nodes ────────────────────────────────────────────────────────
    nx.draw_networkx_nodes(
        G,
        pos,
        node_color=node_colors,
        edgecolors="black",
        linewidths=1.5 if mode == "compact" else 2.0,
        node_size=node_size,
        ax=ax,
    )

    # ── Draw edges via FancyArrowPatch ────────────────────────────────────
    # Shrink endpoints so arrows touch the node boundary, not the centre.
    shrink = np.sqrt(node_size / np.pi)
    alpha = 0.55 if mode == "compact" else 0.80
    mutation = 10 if mode == "compact" else 18
    lw_adj = 1.5 if mode == "compact" else 2.0
    lw_iter = 1.2 if mode == "compact" else 1.8

    _EDGE_STYLE: dict[tuple[str, str], dict] = {
        ("adjacency", "unstable"): {
            "color": "#3b82f6",
            "lw": lw_adj,
            "base_rad": 0.20,
            "ls": "solid",
        },
        ("adjacency", "stable"): {
            "color": "#ef4444",
            "lw": lw_adj,
            "base_rad": -0.20,
            "ls": "solid",
        },
        ("iterate", "unstable"): {
            "color": "#a855f7",
            "lw": lw_iter,
            "base_rad": 0.38,
            "ls": "dashed",
        },
    }
    _FALLBACK = _EDGE_STYLE[("adjacency", "unstable")]

    # Track how many edges have been drawn for each (u, v) pair so that
    # parallel edges get staggered curvature and don't overlap.
    _pair_count: dict[tuple, int] = {}

    for u, v, _key, data in G.edges(keys=True, data=True):
        if u == v:
            continue
        edge_type = data.get("type", "adjacency")
        stability = data.get("stability", "unstable")

        if edge_type == "iterate" and not show_iterate_edges:
            continue

        style = _EDGE_STYLE.get((edge_type, stability), _FALLBACK)

        pair = (u, v)
        idx = _pair_count.get(pair, 0)
        _pair_count[pair] = idx + 1
        rad = style["base_rad"] + idx * 0.15 * np.sign(style["base_rad"] or 1)

        patch = FancyArrowPatch(
            posA=pos[u],
            posB=pos[v],
            connectionstyle=f"arc3,rad={rad}",
            arrowstyle="-|>",
            color=style["color"],
            linewidth=style["lw"],
            linestyle=style["ls"],
            alpha=alpha,
            mutation_scale=mutation,
            shrinkA=shrink,
            shrinkB=shrink,
            transform=ax.transData,
            zorder=2,
        )
        ax.add_patch(patch)

    # ── Labels ────────────────────────────────────────────────────────────
    if label_mode != "none":
        labels: dict[int, str] = {}
        for node in G.nodes():
            data = G.nodes[node]
            coords = data.get("coords")
            u_cd = data.get("unstable_cdist")
            s_cd = data.get("stable_cdist")
            parts: list[str] = []

            if label_mode in ("id", "all"):
                parts.append(str(node))
            if label_mode in ("coords", "all") and coords is not None:
                parts.append(f"({coords[0]:.2f},{coords[1]:.2f})")
            if label_mode in ("cdist", "all"):
                u_str = f"{u_cd:.2f}" if u_cd is not None else "?"
                s_str = f"{s_cd:.2f}" if s_cd is not None else "?"
                parts.append(f"u:{u_str}\ns:{s_str}")

            labels[node] = "\n".join(parts) if parts else str(node)

        font_size = 5 if mode == "compact" else 7
        nx.draw_networkx_labels(
            G, pos, labels, font_size=font_size, font_weight="bold", ax=ax
        )

    # ── Colorbar ──────────────────────────────────────────────────────────
    if cmap_obj is not None and norm_obj is not None:
        sm = cm.ScalarMappable(cmap=cmap_obj, norm=norm_obj)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
        cbar.set_label(
            "Unstable arc-length"
            if node_color_by == "unstable_cdist"
            else "Stable arc-length",
            fontsize=9,
        )

    # ── Axis appearance ───────────────────────────────────────────────────
    _linear_layouts = {"cdist", "stable_linear", "unstable_linear"}
    if layout in _linear_layouts:
        ax.set_axis_on()
        if layout == "stable_linear":
            ax.set_xlabel("Stable arc-length order", fontsize=10)
            ax.set_yticks([])
        elif layout == "unstable_linear":
            ax.set_xlabel("Unstable arc-length order", fontsize=10)
            ax.set_yticks([])
        else:
            ax.set_xlabel("Unstable arc-length (cdist)", fontsize=10)
            ax.set_ylabel("Stable arc-length (cdist)", fontsize=10)
            ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
        ax.margins(0.15)
    else:
        ax.axis("off")

    # ── Legend ────────────────────────────────────────────────────────────
    legend_handles = [
        Line2D([0], [0], color="#3b82f6", linewidth=2, label="Unstable adjacency"),
        Line2D([0], [0], color="#ef4444", linewidth=2, label="Stable adjacency"),
    ]
    if show_iterate_edges:
        legend_handles.append(
            Line2D(
                [0],
                [0],
                color="#a855f7",
                linewidth=1.5,
                linestyle="--",
                label="Iterate",
            )
        )
    ax.legend(handles=legend_handles, loc="upper right", fontsize=9)

    mode_label = f"{n_nodes} nodes · {mode} mode"
    ax.set_title(
        f"Intersection Graph  ({mode_label})", fontsize=12, fontweight="bold"
    )

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")

    plt.show()
    return fig, ax
