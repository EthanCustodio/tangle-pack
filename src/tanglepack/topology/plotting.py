"""
Drawing routines for the topological layer.

Every ``plot_*`` of the topological layer lives here so the marker, colour and
z-order conventions are defined once; :class:`~.Trellis.Trellis` and
:class:`~..loom.TangleSession.TangleSession` keep thin delegating wrappers.

Dev Notes — topological plotting.

The z-orders are a stack, not arbitrary numbers: black intersections (drawn by
the numerical layer) sit lowest, then the dual graph's own edges (z=1) and
face points (z=2) — kept low so a dual graph drawn over a trellis plot never
hides the numerics dots underneath it — then the minimal trellis's bridges
(dropped z=3, image z=4, hole z=5) and the bridges coloured by class (z=6),
then a highlighted dual-graph WALK (z=8: above every edge and bridge, below
the nodes it threads through), then the stable nodes (z=10), then the
magenta strong-pip CANDIDATES (z=11), then the green chosen pip (z=12) on top
of the set it was chosen from, then the orange pseudoneighbours (z=13), then
the holes (z=14) and their labels (z=15). Changing one means checking the
whole ladder.

Colours. The per-class palette (:data:`CLASS_COLORS`) is a fixed-order
categorical set: a class keeps its slot whatever else is drawn, and the same
slot is used for its bridges in the plane and its node in the transition
graph. The names drawn on stable nodes, number lines and graph nodes are
mathtext (``$R_1^{2}$``); :func:`name_mathtext` turns a plain symbol or
element name into it.

:func:`plot_stable_partition` is the odd one out: it takes the partition
RESULTS rather than a trellis, because a number-line figure is routinely built
from partitions gathered across several trellises (see
``scripts/henon_pseudoneighbors_nested.py``). The trellis method passes its own
``stable_partitions``.

:func:`plot_dual_graph` takes a :class:`~.DualGraph.DualGraph` directly, for the
same reason: a graph spans one arrangement but the partitions of several
trellises, so there is no single owning trellis to hang a method off.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from typing import Callable, Iterable, Literal, Optional, Sequence, TYPE_CHECKING, Union

import numpy as np
import matplotlib.pyplot as plt
from numpy.typing import NDArray

from .TopologyResults import PartitionInterval, StablePartitionResult

if TYPE_CHECKING:
    from .DualGraph import DualGraph, FaceNode, StableNode
    from .DualWalk import Walk
    from .MinimalTrellis import MinimalTrellis
    from .SymbolicDynamics import SymbolicDynamics
    from .Trellis import Trellis

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# ── shared style conventions ────────────────────────────────────────────────

#: Scatter defaults for the strong-pip candidate set.
STRONG_PIP_CANDIDATE_STYLE = {"color": "magenta", "s": 7, "zorder": 11}

#: Scatter defaults for the chosen strong pip and its cut points.
STRONG_PIP_STYLE = {"color": "green", "s": 7, "zorder": 12}

#: Scatter defaults for pseudoneighbor pair members.
PSEUDONEIGHBOR_STYLE = {"color": "darkorange", "s": 7, "zorder": 13}

#: Scatter defaults for punched holes (marker and colour come from the orbit).
HOLE_STYLE = {"s": 40, "zorder": 14}

#: One marker per reference-pseudoneighbor orbit: every hole descending from
#: the same reference (matching ``Hole.origin``) shares the symbol.
HOLE_MARKERS = ("x", "+", "*", "^", "s", "D", "v", "P")

#: One colour per reference-pseudoneighbor orbit, matching the markers, so an
#: orbit reads as a single consistent shape+colour across the figure.
HOLE_COLORS = (
    "purple", "teal", "darkgreen", "crimson",
    "chocolate", "navy", "olive", "deeppink",
)

#: Scatter defaults for dual-graph arc nodes (hollow = wall, filled = passable).
DUAL_GRAPH_ARC_STYLE = {"color": "black", "s": 20, "zorder": 10}

#: Scatter defaults for each dual-graph face node's representative point.
DUAL_GRAPH_FACE_STYLE = {"color": "dimgray", "s": 6, "zorder": 2}

#: Style for the thin lines from a face point to its boundary arc nodes.
DUAL_GRAPH_EDGE_STYLE = {"color": "lightgray", "linewidth": 0.5, "zorder": 1}

#: Fraction of a bounding box's diagonal used both to push the unbounded face
#: node's point outside the arc-node bbox and to pad ``plot_dual_graph``'s
#: ``clip_to_arcs`` axes limits.
DUAL_GRAPH_PUSH_FRACTION = 0.1

#: ``**scatter_kwargs`` names :func:`plot_dual_graph` manages itself (they
#: distinguish a hollow "wall" arc node from a filled "passable" one) and so
#: refuses to let a caller override.
_RESERVED_ARC_KWARGS = frozenset({"facecolors", "c"})


def _scatter_intersections(
    trellis: "Trellis",
    intersection_ids: Iterable[int],
    ax,
    style: dict,
    **scatter_kwargs,
):
    """Scatter the coordinates of intersection ids with a style's defaults."""
    coords = np.array([trellis.registry[iid].coords for iid in intersection_ids])
    for key, value in style.items():
        scatter_kwargs.setdefault(key, value)
    target = ax if ax is not None else plt
    return target.scatter(coords[:, 0], coords[:, 1], **scatter_kwargs)


# ── strong pips ─────────────────────────────────────────────────────────────


def plot_strong_pip_candidates(trellis: "Trellis", ax=None, **scatter_kwargs):
    """
    Scatter-plot a trellis's strong-pip candidates in magenta, on top of the tangle.

    Draw this before :func:`plot_strong_pip` so the single chosen pip (green) sits
    on top of the candidate set (magenta) it was chosen from.

    Args:
        trellis: The trellis whose candidates are drawn.
        ax: Optional matplotlib Axes to draw on. Defaults to the current axes
            (plt).
        **scatter_kwargs: Forwarded to scatter, overriding
            :data:`STRONG_PIP_CANDIDATE_STYLE` (e.g. pass s=30 for larger
            markers).

    Returns:
        The matplotlib PathCollection, or None if there are no candidates.
    """
    if not trellis.strong_pip_candidates:
        logger.info("No strong-pip candidates; call classify_strong_pips() first.")
        return None

    return _scatter_intersections(
        trellis,
        trellis.strong_pip_candidates,
        ax,
        STRONG_PIP_CANDIDATE_STYLE,
        **scatter_kwargs,
    )


def plot_strong_pip(trellis: "Trellis", ax=None, **scatter_kwargs):
    """
    Scatter-plot a trellis's strong-pip cut points in green, on top of the tangle.

    For a period-1 anchor this is the single chosen strong pip; for a period-k
    anchor it is the k points (strong pip + iterates) that bound the resonance
    zone — see :meth:`~tanglepack.topology.Trellis.Trellis.strong_pip_cut_points`.
    Call after plotting the tangle and intersections to highlight them.

    Args:
        trellis: The trellis whose strong pip is drawn.
        ax: Optional matplotlib Axes to draw on. Defaults to the current axes
            (plt).
        **scatter_kwargs: Forwarded to scatter, overriding
            :data:`STRONG_PIP_STYLE` (e.g. pass s=40 for a larger marker).

    Returns:
        The matplotlib PathCollection, or None if no strong pip has been chosen.
    """
    if trellis.strong_pip is None:
        logger.info(
            "No strong pip chosen; call classify_strong_pips() (or set_strong_pip()) first."
        )
        return None

    return _scatter_intersections(
        trellis,
        trellis.strong_pip_cut_points(),
        ax,
        STRONG_PIP_STYLE,
        **scatter_kwargs,
    )


# ── pseudoneighbors and holes ───────────────────────────────────────────────


def plot_pseudoneighbors(
    trellis: "Trellis",
    ax=None,
    *,
    include_trajectories: bool = True,
    **scatter_kwargs,
):
    """
    Scatter-plot a trellis's pseudoneighbor pair members, on top of the tangle.

    Every branch has pseudoneighbors — only the reference computation is
    restricted to the fundamental segment — so by default all recorded pairs
    (references and their iterated appearances) are drawn.

    Args:
        trellis: The trellis whose pairs are drawn.
        ax: Optional matplotlib Axes to draw on. Defaults to the current axes
            (plt).
        include_trajectories: If True (default) every recorded pair is drawn;
            pass False for the reference pairs only.
        **scatter_kwargs: Forwarded to scatter, overriding
            :data:`PSEUDONEIGHBOR_STYLE`.

    Returns:
        The matplotlib PathCollection, or None if there are no pairs.
    """
    pairs = (
        trellis.pseudoneighbors if include_trajectories
        else trellis.reference_pseudoneighbors
    )
    if not pairs:
        logger.info(
            "No pseudoneighbors to plot; call compute_pseudoneighbors() first."
        )
        return None

    ids = sorted({i for p in pairs for i in p.as_tuple()})
    return _scatter_intersections(
        trellis, ids, ax, PSEUDONEIGHBOR_STYLE, **scatter_kwargs
    )


def plot_holes(
    trellis: "Trellis", ax=None, *, show_iterates: bool = True, **scatter_kwargs
):
    """
    Scatter-plot a trellis's punched holes, on top of the tangle.

    Each reference pseudoneighbor's orbit gets its own marker symbol AND colour
    — the reference hole and every hole generated from it (forward, backward,
    or propagated) share them, from :data:`HOLE_MARKERS` / :data:`HOLE_COLORS`
    — and each hole is labelled with its iterate (0 = reference, negative =
    backward). Sides are reported by the partition, not encoded here.

    Args:
        trellis: The trellis whose holes are drawn.
        ax: Optional matplotlib Axes to draw on. Defaults to the current axes
            (plt).
        show_iterates: Annotate each hole with its iterate number (default True).
        **scatter_kwargs: Forwarded to scatter, overriding :data:`HOLE_STYLE`.
            A ``color=`` override disables the by-orbit colouring, a ``marker=``
            override the by-orbit markers.

    Returns:
        List of matplotlib PathCollections (one per orbit drawn), or None if
        there are no holes.
    """
    if not trellis.holes:
        logger.info("No holes to plot; call punch_holes() first.")
        return None

    for key, value in HOLE_STYLE.items():
        scatter_kwargs.setdefault(key, value)
    origins = sorted({h.origin for h in trellis.holes if h.origin is not None})
    style_of = {
        origin: (
            HOLE_MARKERS[i % len(HOLE_MARKERS)],
            HOLE_COLORS[i % len(HOLE_COLORS)],
        )
        for i, origin in enumerate(origins)
    }
    target = ax if ax is not None else plt

    handles = []
    groups = sorted(
        {h.origin for h in trellis.holes},
        key=lambda origin: (origin is None, origin or ()),
    )
    for origin in groups:
        batch = [h for h in trellis.holes if h.origin == origin]
        coords = np.array([h.coords for h in batch])
        marker, color = style_of.get(origin, ("x", "gray"))
        kwargs = dict(scatter_kwargs)
        kwargs.setdefault("marker", marker)
        kwargs.setdefault("color", color)
        handles.append(target.scatter(coords[:, 0], coords[:, 1], **kwargs))
        if show_iterates:
            axes = target if ax is not None else plt.gca()
            for hole in batch:
                if hole.iterate is None:
                    continue
                axes.annotate(
                    str(hole.iterate),
                    hole.coords,
                    textcoords="offset points",
                    xytext=(4, 4),
                    fontsize=8,
                    color=color,
                    zorder=15,
                )
    return handles


# ── the stable-manifold partition number line ───────────────────────────────

#: An interval narrower than this fraction of its row's span has its element
#: label staggered over two text rows (see ``plot_stable_partition``).
ELEMENT_LABEL_NARROW_FRACTION = 0.06


def plot_stable_partition(
    results: Union[StablePartitionResult, Iterable[StablePartitionResult]],
    ax=None,
    *,
    labels: Optional[Sequence[str]] = None,
    element_labels: Optional[
        Callable[[StablePartitionResult, PartitionInterval], str]
    ] = None,
    **line_kwargs,
):
    """
    Draw stable partitions as number lines (stable cdist on the x-axis).

    Each partition occupies one row; its intervals are horizontal segments
    with a bracket at each end — [ ] closed, ( ) open — every boundary
    labelled with its intersection id, and a circle marking each singleton
    piece [x, x] (a point pinched between two open intervals). Dotted
    vertical lines mark the beginning and end of the manifold.

    Args:
        results: One StablePartitionResult or an iterable of them.
        ax: Optional matplotlib Axes. Defaults to the current axes.
        labels: Optional y-axis label per row (one per result, in order),
            replacing the default ``"<side> (p<period>, orbit <n>)"`` — for
            instance to tell a homotopy row from an iterated row of the same
            branch and side on one axes.
        element_labels: Optional callable ``(result, interval) -> str``; when
            given, its text (mathtext allowed, e.g. ``"$R_1^{2}$"``) is drawn
            at each interval's midpoint just above the row. Return ``""`` to
            leave an interval unlabelled.
        **line_kwargs: Forwarded to the interval ``plot`` calls (e.g.
            ``linewidth``, ``color``).

    Raises:
        ValueError: If ``labels`` is given with a different length than
            ``results``.

    Returns:
        The Axes drawn on, or None if there was nothing to draw.
    """
    if isinstance(results, StablePartitionResult):
        results = [results]
    results = list(results)
    if not results:
        logger.info("No partitions to plot; call partition_stable_manifold() first.")
        return None

    target = ax if ax is not None else plt.gca()
    line_kwargs.setdefault("linewidth", 2)

    if labels is not None and len(labels) != len(results):
        raise ValueError(
            f"{len(labels)} row labels were given for {len(results)} partition(s)"
        )
    row_labels: list[str] = []
    extent_lo, extent_hi = np.inf, -np.inf
    for row, result in enumerate(results):
        color = "tab:blue" if result.side == "left" else "tab:red"
        row_labels.append(
            labels[row]
            if labels is not None
            else f"{result.side} (p{result.branch_key[0].period}, "
            f"orbit {result.branch_key[2]})"
        )
        labelled: set[int] = set()
        # Element labels of NARROW intervals (under this fraction of the
        # row's span) alternate between two text rows so a run of tightly
        # packed elements near the anchor stays legible.
        row_span = max(
            (i.hi_cdist for i in result.intervals), default=0.0
        ) - min((i.lo_cdist for i in result.intervals), default=0.0)
        narrow_count = 0
        for interval in result.intervals:
            extent_lo = min(extent_lo, interval.lo_cdist)
            extent_hi = max(extent_hi, interval.hi_cdist)
            target.plot(
                [interval.lo_cdist, interval.hi_cdist], [row, row],
                color=color, solid_capstyle="butt", **line_kwargs,
            )
            # Bracket per interval END, nudged inward so the two intervals
            # meeting at one intersection point stay individually readable:
            # e.g. "](" = closed on the anchor side, open on the outward side.
            if interval.lo_cdist == interval.hi_cdist:
                # A singleton piece [x, x] — a point pinched between two open
                # intervals. Brackets are illegible at zero width; a FILLED
                # circle marks the closed point, with the neighbours' open
                # parentheses either side.
                target.plot(
                    interval.lo_cdist, row, marker="o", markersize=9,
                    markerfacecolor=color, markeredgecolor=color,
                    zorder=7,
                )
            else:
                target.annotate(
                    "[" if interval.closed_lo else "(",
                    (interval.lo_cdist, row), textcoords="offset points",
                    xytext=(1, 0), ha="left", va="center",
                    fontsize=13, fontweight="bold", color=color, zorder=6,
                )
                target.annotate(
                    "]" if interval.closed_hi else ")",
                    (interval.hi_cdist, row), textcoords="offset points",
                    xytext=(-1, 0), ha="right", va="center",
                    fontsize=13, fontweight="bold", color=color, zorder=6,
                )
            if element_labels is not None:
                text = element_labels(result, interval)
                if text:
                    width = interval.hi_cdist - interval.lo_cdist
                    narrow = width < ELEMENT_LABEL_NARROW_FRACTION * row_span
                    lift = 7
                    if narrow:
                        lift = 7 if narrow_count % 2 == 0 else 18
                        narrow_count += 1
                    target.annotate(
                        text,
                        (0.5 * (interval.lo_cdist + interval.hi_cdist), row),
                        textcoords="offset points", xytext=(0, lift),
                        ha="center", va="bottom", fontsize=9, color="black",
                        zorder=8,
                    )
            # Label each boundary with the intersection's id in the tangle,
            # staggered over two text rows so tightly packed boundaries near
            # the anchor stay legible.
            for boundary_id, cdist in (
                (interval.lo_id, interval.lo_cdist),
                (interval.hi_id, interval.hi_cdist),
            ):
                if boundary_id is None or boundary_id in labelled:
                    continue
                drop = -14 if len(labelled) % 2 == 0 else -24
                labelled.add(boundary_id)
                target.annotate(
                    str(boundary_id), (cdist, row), textcoords="offset points",
                    xytext=(0, drop), ha="center", va="top", fontsize=8,
                    color="black",
                )

    # Beginning and end of the (trimmed) stable manifold.
    for cdist, name, align, nudge in (
        (extent_lo, "anchor", "left", 3),
        (extent_hi, "manifold end", "right", -3),
    ):
        target.axvline(cdist, color="gray", linestyle=":", alpha=0.6, zorder=0)
        target.annotate(
            name, (cdist, 1.0), xycoords=("data", "axes fraction"),
            textcoords="offset points", xytext=(nudge, -10), ha=align,
            fontsize=8, color="gray",
        )

    target.set_yticks(range(len(results)))
    target.set_yticklabels(row_labels)
    target.set_ylim(-0.6, len(results) - 0.4)
    target.set_xlabel("stable canonical distance")
    target.set_title(
        "Stable manifold partition — [ ] closed, ( ) open, ○ singleton; "
        "numbers = ids"
    )
    return target


# ── the minimal trellis ─────────────────────────────────────────────────────

#: Line defaults for the bridges a minimal trellis DROPPED.
MINIMAL_TRELLIS_DROPPED_STYLE = {"color": "lightgray", "linewidth": 0.6, "zorder": 3}
#: Line defaults for the hole bridges (the ones the homotopy partition is cut by).
MINIMAL_TRELLIS_HOLE_STYLE = {"color": "tab:purple", "linewidth": 1.8, "zorder": 5}
#: Line defaults for the image bridges (the ones the iterated partition adds).
MINIMAL_TRELLIS_IMAGE_STYLE = {
    "color": "tab:orange",
    "linewidth": 1.8,
    "linestyle": "--",
    "zorder": 4,
}


def plot_minimal_trellis(
    minimal: "MinimalTrellis",
    ax=None,
    *,
    show_dropped: bool = True,
    **line_kwargs,
):
    """
    Draw the bridges of a minimal trellis over the plane.

    Hole bridges are drawn solid (:data:`MINIMAL_TRELLIS_HOLE_STYLE`), image
    bridges dashed (:data:`MINIMAL_TRELLIS_IMAGE_STYLE`), and the bridges the
    reduction dropped in light gray (:data:`MINIMAL_TRELLIS_DROPPED_STYLE`)
    so the reader sees what was left out. The stable manifold is not drawn
    here; draw the tangle underneath first.

    Args:
        minimal: The minimal trellis to draw.
        ax: Optional matplotlib Axes. Defaults to the current axes.
        show_dropped: Draw the dropped bridges too (default True).
        **line_kwargs: Overrides applied to every bridge line (e.g.
            ``linewidth``).

    Returns:
        The Axes drawn on.
    """
    target = ax if ax is not None else plt.gca()
    groups = [
        (minimal.hole_bridge_ids, MINIMAL_TRELLIS_HOLE_STYLE),
        (minimal.image_bridge_ids, MINIMAL_TRELLIS_IMAGE_STYLE),
    ]
    if show_dropped:
        groups.insert(0, (minimal.dropped_bridge_ids, MINIMAL_TRELLIS_DROPPED_STYLE))
    for bridge_ids, base in groups:
        style = dict(base)
        style.update(line_kwargs)
        for bridge_id in bridge_ids:
            bridge = minimal.trellis.bridge_between(*bridge_id)
            if bridge is None:
                logger.debug("bridge %s has no object to draw; skipping it", bridge_id)
                continue
            points = bridge.get_point_array()
            if len(points) < 2:
                continue
            target.plot(points[:, 0], points[:, 1], **style)
    return target


def minimal_trellis_legend_handles() -> list:
    """
    Legend proxies matching :func:`plot_minimal_trellis`'s three bridge kinds.

    Returns:
        ``Line2D`` handles labelled "hole bridge", "image bridge", "dropped
        bridge", built from the same style tables the plotter uses.
    """
    from matplotlib.lines import Line2D

    return [
        Line2D([0], [0], label="hole bridge", **MINIMAL_TRELLIS_HOLE_STYLE),
        Line2D([0], [0], label="image bridge", **MINIMAL_TRELLIS_IMAGE_STYLE),
        Line2D([0], [0], label="dropped bridge", **MINIMAL_TRELLIS_DROPPED_STYLE),
    ]


# ── the dual graph ──────────────────────────────────────────────────────────

#: How far off its edge a side node is drawn: this fraction of the diagonal of
#: the bounding box of every edge midpoint in the graph.
DUAL_GRAPH_SIDE_OFFSET_FRACTION = 0.015


def _edge_midpoints(
    trellis: "Trellis", nodes: Iterable["StableNode"]
) -> list[NDArray[np.float64]]:
    """The usable (non-degenerate) midpoints of a run of stable nodes."""
    mids: list[NDArray[np.float64]] = []
    for node in nodes:
        mid = node.midpoint(trellis)
        if mid is None:
            logger.debug(
                "stable node %s has a degenerate midpoint; skipping it", node.key
            )
            continue
        mids.append(mid)
    return mids


def _midpoint_bbox(dual_graph: "DualGraph"):
    """``(mins, maxs, diagonal)`` of every edge midpoint, or None with none."""
    all_mids = _edge_midpoints(
        dual_graph.arrangement.trellis, dual_graph.stable_nodes.values()
    )
    if not all_mids:
        return None
    points = np.vstack(all_mids)
    mins = points.min(axis=0)
    maxs = points.max(axis=0)
    return mins, maxs, float(np.linalg.norm(maxs - mins))


def _anchorward_look(
    node: "StableNode", trellis: "Trellis"
) -> Optional[NDArray[np.float64]]:
    """
    The unit direction toward the anchor at the middle of a stable edge.

    Read from the two polyline vertices flanking the half-length point of the
    edge's natural (lo -> hi, anchor outward) polyline, negated; None for a
    degenerate polyline.
    """
    polyline = np.asarray(node.arc.polyline(trellis), dtype=float)
    if len(polyline) < 2:
        return None
    steps = np.linalg.norm(np.diff(polyline, axis=0), axis=1)
    total = float(steps.sum())
    if total <= 0.0:
        return None
    cumulative = np.concatenate([[0.0], np.cumsum(steps)])
    index = int(np.searchsorted(cumulative, 0.5 * total, side="right") - 1)
    index = min(max(index, 0), len(polyline) - 2)
    outward = polyline[index + 1] - polyline[index]
    norm = float(np.linalg.norm(outward))
    if norm <= 0.0:
        return None
    return -outward / norm


def stable_node_point(
    dual_graph: "DualGraph", node: "StableNode"
) -> Optional[NDArray[np.float64]]:
    """
    Where to draw one stable node.

    A unified (solid) node sits ON the midpoint of its edge. A side node sits
    just off it, on its side of the stable manifold: displaced along the
    left normal of the anchorward look (positive cross = left, the
    partition's own convention) by :data:`DUAL_GRAPH_SIDE_OFFSET_FRACTION` of
    the graph's midpoint bounding-box diagonal — to the left for a
    ``"left"`` node, to the right for a ``"right"`` one.

    Args:
        dual_graph: The graph ``node`` belongs to.
        node: The stable node to place.

    Returns:
        The ``(2,)`` point, or None when the edge's polyline is degenerate.
    """
    trellis = dual_graph.arrangement.trellis
    mid = node.midpoint(trellis)
    if mid is None:
        return None
    mid = np.asarray(mid, dtype=float)
    if node.is_unified:
        return mid
    look = _anchorward_look(node, trellis)
    bbox = _midpoint_bbox(dual_graph)
    if look is None or bbox is None:
        return mid
    normal = np.array([-look[1], look[0]])  # the left of the anchorward look
    offset = DUAL_GRAPH_SIDE_OFFSET_FRACTION * bbox[2]
    return mid + (normal if node.sides[0] == "left" else -normal) * offset


def _push_outside_bbox(
    dual_graph: "DualGraph", point: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Push a point outside the bounding box of every edge midpoint.

    Lands on the bbox corner in ``point``'s own quadrant (relative to the bbox
    centre), then a further :data:`DUAL_GRAPH_PUSH_FRACTION` of the bbox
    diagonal outward.
    """
    bbox = _midpoint_bbox(dual_graph)
    if bbox is None:
        return point
    mins, maxs, diagonal = bbox
    center = 0.5 * (mins + maxs)

    direction = point - center
    sign = np.sign(direction)
    sign[sign == 0.0] = 1.0
    corner = np.where(sign >= 0.0, maxs, mins)
    unit = sign / np.linalg.norm(sign)
    return corner + DUAL_GRAPH_PUSH_FRACTION * diagonal * unit


def face_point(
    dual_graph: "DualGraph", face_node: "FaceNode"
) -> NDArray[np.float64]:
    """
    A representative point for one face node of a dual graph.

    Args:
        dual_graph: The graph ``face_node`` belongs to.
        face_node: The face node to place.

    Returns:
        The ``(2,)`` point: a COPY of a region's own
        :attr:`~.TopologyResults.Region.representative_point` when the node
        stands for exactly one closed minimal face (never the region's own
        cached array); otherwise the mean of the midpoints of its boundary
        stable edges. Only the ONE node with
        :attr:`~.DualGraph.FaceNode.is_unbounded` set has that mean pushed
        outside the bounding box of every edge midpoint in the graph (see
        :func:`_push_outside_bbox`) — a bounded face's point is returned as
        computed, with no guarantee it falls outside any other face (a
        region's representative point is itself only a heuristic).

    Note:
        A region's ``representative_point`` can itself be None (a boundary too
        degenerate to have an interior); that, like an "open"/"outer" node,
        falls back to the edge-midpoint mean.
    """
    if face_node.kind == "region" and face_node.faces:
        point = face_node.faces[0].representative_point
        if point is not None:
            return point.copy()
        logger.debug(
            "face node %d (region) has no representative point; falling back "
            "to the mean of its edge midpoints",
            face_node.index,
        )

    mids = _edge_midpoints(dual_graph.arrangement.trellis, face_node.stable_node_list)
    if mids:
        mean = np.mean(np.vstack(mids), axis=0).astype(np.float64)
    else:
        mean = np.zeros(2, dtype=np.float64)
        logger.debug(
            "face node %d has no edge midpoints to average; defaulting to the "
            "origin before any outward push",
            face_node.index,
        )

    if not face_node.is_unbounded:
        return mean
    return _push_outside_bbox(dual_graph, mean)


def _clip_axes_to_arcs(ax, dual_graph: "DualGraph") -> None:
    """
    Limit an axes to the padded bbox of every edge midpoint and the unbounded
    face point.

    A bounded region's :attr:`~.TopologyResults.Region.representative_point`
    is only a heuristic and can, in a degenerate case, land far outside a
    tight tangle; left to matplotlib's autoscale that single outlier point
    blows the view out past anything useful, so :func:`plot_dual_graph`'s
    default view is clipped to what the edges (and the one unbounded point)
    actually span instead.
    """
    points = _edge_midpoints(
        dual_graph.arrangement.trellis, dual_graph.stable_nodes.values()
    )
    points = points + [face_point(dual_graph, dual_graph.unbounded)]
    if not points:
        return
    stacked = np.vstack(points)
    mins = stacked.min(axis=0)
    maxs = stacked.max(axis=0)
    pad = DUAL_GRAPH_PUSH_FRACTION * float(np.linalg.norm(maxs - mins))
    ax.set_xlim(mins[0] - pad, maxs[0] + pad)
    ax.set_ylim(mins[1] - pad, maxs[1] + pad)


#: ``_1``, ``^2`` and ``^-1`` in a plain symbol or element name.
_SCRIPT_PATTERN = re.compile(r"([_^])(-?[A-Za-z0-9]+)")


def name_mathtext(text: str) -> str:
    """
    A plain symbol or element name as matplotlib mathtext.

    ``R_1^2`` becomes ``$R_{1}^{2}$``, ``a_2^-1`` becomes ``$a_{2}^{-1}$``,
    ``u`` becomes ``$u$``; a branch-tagged element name keeps its tag as
    plain prefix text outside the math (``p3@1.0:R_2`` -> ``p3@1.0:$R_{2}$``,
    exactly as :attr:`~.ElementNaming.ElementName.mathtext` prints it).
    Already dollar-wrapped text is returned unchanged.

    Args:
        text: The name, as :attr:`~.ElementNaming.ElementName.text`,
            :attr:`~.SymbolicDynamics.Symbol.text` or a
            :attr:`~.SymbolicDynamics.RefinedClass.name` print it.

    Returns:
        The mathtext string.
    """
    if not text:
        return ""
    if text.startswith("$") and text.endswith("$"):
        return text
    tag, sep, body = text.rpartition(":")
    body = _SCRIPT_PATTERN.sub(lambda m: f"{m.group(1)}{{{m.group(2)}}}", body)
    return f"{tag}:${body}$" if sep else f"${body}$"


def stable_node_label(node: "StableNode", label_style: Literal["name", "ref"] = "name") -> str:
    """
    The annotation :func:`plot_dual_graph` draws on one stable node.

    Args:
        node: The stable node.
        label_style: ``"name"`` for the element names as mathtext, one per
            faced side joined by ``" | "`` (a side with no name falls back to
            its element label); ``"ref"`` for the element labels joined by
            ``"|"``.

    Returns:
        The text.
    """
    if label_style == "ref":
        return "|".join(node.elements[side].label for side in node.sides)
    return " | ".join(
        node.names[side].mathtext if side in node.names else node.elements[side].label
        for side in node.sides
    )


def plot_dual_graph(
    dual_graph: "DualGraph",
    ax=None,
    *,
    show_labels: bool = False,
    clip_to_arcs: bool = True,
    label_style: Literal["name", "ref"] = "name",
    **scatter_kwargs,
):
    """
    Draw a dual graph over the plane.

    Every stable node is a circle: an OPEN circle just off its edge on its own
    side for a wall (a side node), a SOLID circle on the edge for a unified,
    traversable node (:func:`stable_node_point` places both). Every face node
    is a small point (:func:`face_point`) joined by a thin line to the stable
    node on its side of each stable edge on its boundary. Bridges are not
    drawn; draw the tangle underneath first.

    Args:
        dual_graph: The graph to draw.
        ax: Optional matplotlib Axes. Defaults to the current axes.
        show_labels: Annotate each stable node with its element name(s) or
            label(s) (see ``label_style``) and each face node with
            ``index:kind``.
        clip_to_arcs: Limit the axes to the padded bounding box of the edge
            midpoints and the unbounded face point (see
            :func:`_clip_axes_to_arcs`). Pass False to keep matplotlib's
            autoscale, e.g. when drawing over a tangle whose view is set
            elsewhere.
        label_style: What a stable node's annotation says. ``"name"`` (the
            default) draws the element names as mathtext, one per faced side
            joined by ``" | "`` (``$R_1^{2}$ | $L_3$``), falling back to the
            element label on a side with no name; ``"ref"`` draws the
            :attr:`~.TopologyResults.ElementRef.label` texts joined by
            ``"|"`` (``p1@0.0/R#1|p1@0.0/L#2``).
        **scatter_kwargs: Overrides for the stable-node scatters on top of
            :data:`DUAL_GRAPH_ARC_STYLE` — ``color`` (the solid fill and,
            unless ``edgecolors`` is given, the open ring), ``edgecolors``,
            ``s``, ``zorder``, ...

    Returns:
        The Axes drawn on.

    Raises:
        ValueError: If ``facecolors`` or ``c`` is passed — open versus solid
            is what this plotter decides; use ``color`` / ``edgecolors``. Or
            if ``label_style`` is neither ``"name"`` nor ``"ref"``.

    Note:
        The two stable-node scatters are ALWAYS drawn, even when one set is
        empty, so they land at fixed positions 0 (open) and 1 (solid) in
        ``ax.collections`` — callers (tests included) can tell them apart
        without inspecting facecolors.
    """
    reserved = _RESERVED_ARC_KWARGS.intersection(scatter_kwargs)
    if reserved:
        raise ValueError(
            f"plot_dual_graph manages {sorted(reserved)} itself (open vs. "
            "solid stable nodes); pass 'color'/'edgecolors' instead"
        )
    if label_style not in ("name", "ref"):
        raise ValueError(
            f"label_style must be 'name' or 'ref', not {label_style!r}"
        )
    target = ax if ax is not None else plt.gca()

    style = dict(DUAL_GRAPH_ARC_STYLE)
    style.update(scatter_kwargs)
    color = style.pop("color")
    edgecolors = style.pop("edgecolors", color)

    placed: dict[tuple, NDArray[np.float64]] = {}
    open_points: list[NDArray[np.float64]] = []
    solid_points: list[NDArray[np.float64]] = []
    for node in dual_graph.stable_nodes.values():
        point = stable_node_point(dual_graph, node)
        if point is None:
            logger.debug("stable node %s has a degenerate midpoint; skipping it", node.key)
            continue
        placed[node.key] = point
        (solid_points if node.is_unified else open_points).append(point)

    open_coords = np.vstack(open_points) if open_points else np.empty((0, 2))
    solid_coords = np.vstack(solid_points) if solid_points else np.empty((0, 2))
    target.scatter(
        open_coords[:, 0], open_coords[:, 1],
        facecolors="none", edgecolors=edgecolors, **style,
    )
    target.scatter(solid_coords[:, 0], solid_coords[:, 1], color=color, **style)

    if show_labels:
        for node in dual_graph.stable_nodes.values():
            point = placed.get(node.key)
            if point is None:
                continue
            target.annotate(
                stable_node_label(node, label_style), point,
                textcoords="offset points", xytext=(3, 3), fontsize=7,
            )

    face_points: list[NDArray[np.float64]] = []
    for face_node in dual_graph.face_nodes:
        hub = face_point(dual_graph, face_node)
        face_points.append(hub)
        for stable_node, _side in face_node.stable_nodes:
            point = placed.get(stable_node.key)
            if point is None:
                continue
            target.plot(
                [hub[0], point[0]], [hub[1], point[1]], **DUAL_GRAPH_EDGE_STYLE,
            )
        if show_labels:
            target.annotate(
                f"{face_node.index}:{face_node.kind}", hub,
                textcoords="offset points", xytext=(3, -3), fontsize=7,
                color="dimgray",
            )
    if face_points:
        stacked = np.vstack(face_points)
        target.scatter(stacked[:, 0], stacked[:, 1], **DUAL_GRAPH_FACE_STYLE)

    if clip_to_arcs:
        _clip_axes_to_arcs(target, dual_graph)
    return target


def dual_graph_legend_handles() -> list:
    """
    Legend proxies matching :func:`plot_dual_graph`'s four marks.

    Returns:
        ``Line2D`` handles labelled "stable node (open, wall)", "stable node
        (solid, traversable)", "face node" and "face-node edge", built from the
        same style tables the plotter uses (marker sizes are
        ``sqrt(style["s"])`` in points).
    """
    from matplotlib.lines import Line2D

    arc_size = float(np.sqrt(DUAL_GRAPH_ARC_STYLE["s"]))
    face_size = float(np.sqrt(DUAL_GRAPH_FACE_STYLE["s"]))
    color = DUAL_GRAPH_ARC_STYLE["color"]
    return [
        Line2D(
            [0], [0], marker="o", linestyle="none", markersize=arc_size,
            markerfacecolor="none", markeredgecolor=color,
            label="stable node (open, wall)",
        ),
        Line2D(
            [0], [0], marker="o", linestyle="none", markersize=arc_size,
            markerfacecolor=color, markeredgecolor=color,
            label="stable node (solid, traversable)",
        ),
        Line2D(
            [0], [0], marker="o", linestyle="none", markersize=face_size,
            color=DUAL_GRAPH_FACE_STYLE["color"], label="face node",
        ),
        Line2D(
            [0], [0], color=DUAL_GRAPH_EDGE_STYLE["color"],
            linewidth=DUAL_GRAPH_EDGE_STYLE["linewidth"], label="face-node edge",
        ),
    ]


# ── a walk through the dual graph ───────────────────────────────────────────

#: Line defaults for a highlighted walk: above every edge and bridge (z=8),
#: below the stable nodes it threads through (z=10).
WALK_STYLE = {"color": "#1baf7a", "linewidth": 2.2, "zorder": 8, "solid_capstyle": "round"}

#: Marker size (points) of the walk's start dot and end arrowhead.
WALK_MARKER_SIZE = 7.0


def _walk_points(dual_graph: "DualGraph", walk: "Walk") -> list[NDArray[np.float64]]:
    """The polyline vertices of a walk: start face, then node and face per step."""
    if not walk.faces:
        return []
    points = [face_point(dual_graph, walk.faces[0])]
    for step in walk.steps:
        point = stable_node_point(dual_graph, step.node)
        if point is None:
            logger.debug(
                "walk step through stable node %s has no drawable point; the walk "
                "line skips it",
                step.node.key,
            )
        else:
            points.append(point)
        points.append(face_point(dual_graph, step.to_face))
    return points


def plot_walk(dual_graph: "DualGraph", walk: "Walk", ax=None, **line_kwargs):
    """
    Highlight one walk of the dual graph as a polyline over the plane.

    The line starts at the walk's first face (:func:`face_point`) and, per
    step, passes the stable node it crosses (:func:`stable_node_point`) and
    lands in the face it enters, so it threads through exactly the points
    :func:`plot_dual_graph` draws. A dot marks the start and an arrowhead the
    end. Draw the dual graph first (or after: the walk sits between the graph's
    edges and its nodes in z-order).

    Args:
        dual_graph: The graph the walk was found in.
        walk: A :class:`~.DualWalk.Walk` (anything with ``faces`` and
            ``steps`` whose steps carry ``node`` and ``to_face``).
        ax: Optional matplotlib Axes. Defaults to the current axes.
        **line_kwargs: Overrides on :data:`WALK_STYLE` (``color``,
            ``linewidth``, ``zorder``, ``linestyle``, ...); the end arrowhead
            and start dot follow the line's colour and z-order.

    Returns:
        The ``Line2D`` of the walk, or None when the walk visits no face.
    """
    target = ax if ax is not None else plt.gca()
    points = _walk_points(dual_graph, walk)
    if not points:
        logger.info("the walk visits no face; nothing to draw")
        return None
    style = dict(WALK_STYLE)
    style.update(line_kwargs)
    stacked = np.vstack(points)
    (line,) = target.plot(stacked[:, 0], stacked[:, 1], **style)
    color = line.get_color()
    zorder = line.get_zorder()
    target.plot(
        stacked[0, 0], stacked[0, 1], marker="o", markersize=WALK_MARKER_SIZE,
        color=color, linestyle="none", zorder=zorder,
    )
    if len(stacked) >= 2:
        target.annotate(
            "", xy=stacked[-1], xytext=stacked[-2],
            arrowprops={
                "arrowstyle": "-|>",
                "color": color,
                "linewidth": 0,
                "mutation_scale": 2.0 * WALK_MARKER_SIZE,
                "shrinkA": 0,
                "shrinkB": 0,
            },
            zorder=zorder,
        )
    return line


def walk_legend_handles(**line_kwargs) -> list:
    """
    Legend proxies matching :func:`plot_walk`.

    Args:
        **line_kwargs: The same overrides passed to :func:`plot_walk`, so the
            proxy matches the line drawn.

    Returns:
        One ``Line2D`` handle labelled "walk", with the walk's start dot as
        its marker.
    """
    from matplotlib.lines import Line2D

    style = dict(WALK_STYLE)
    style.update(line_kwargs)
    style.pop("zorder", None)
    style.pop("solid_capstyle", None)
    return [
        Line2D(
            [0], [0], marker="o", markersize=WALK_MARKER_SIZE, label="walk", **style
        ),
    ]


# ── bridges by class and the symbolic dynamics ──────────────────────────────

#: Fixed-order categorical colours for the (refined) bridge classes: a class
#: keeps its slot whatever else is on the axes, and the same slot colours its
#: node in :func:`plot_transition_graph`. Past eight, slots repeat.
CLASS_COLORS = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100",
    "#e87ba4", "#008300", "#4a3aa7", "#e34948",
)

#: Line defaults for a classed bridge (inert classes are additionally dashed).
CLASS_BRIDGE_STYLE = {"linewidth": 2.0, "zorder": 6, "solid_capstyle": "round"}

#: Colour of a symbol the palette does not know (a virtual class, say).
CLASS_UNKNOWN_COLOR = "#898781"

#: ``**line_kwargs`` names :func:`plot_bridges_by_class` manages itself.
_RESERVED_CLASS_KWARGS = frozenset({"color", "c", "linestyle", "ls", "label"})


def _class_slots(dynamics: "SymbolicDynamics", *, refined: bool) -> list[tuple[str, bool]]:
    """
    The symbol names of a dynamics in palette order, each with its inertness.

    Classes come in ``dynamics.classes`` order; with ``refined`` a class split
    into refined children contributes the children's names (``a_1, a_2``) in
    index order instead of its letter, followed by the letter itself when a
    member of the class matched no child (see :func:`_unmatched_parent_slots`).
    """
    slots: list[tuple[str, bool]] = []
    unmatched = _unmatched_parent_slots(dynamics, refined=refined)
    for bridge_class, class_dynamics in dynamics.classes.items():
        inert = class_dynamics.kind != "active"
        children = dynamics.refined.get(bridge_class, []) if refined else []
        if children:
            for child in sorted(children, key=lambda c: c.index):
                slots.append((child.name, inert))
            if class_dynamics.letter in unmatched:
                slots.append((class_dynamics.letter, inert))
        else:
            slots.append((class_dynamics.letter, inert))
    return slots


def _unmatched_parent_slots(dynamics: "SymbolicDynamics", *, refined: bool) -> set[str]:
    """
    The letters of split classes with a member that matched no refined child.

    Such a member is drawn under its parent letter, which therefore needs a
    palette slot of its own (after the children); empty unless ``refined``.
    """
    if not refined:
        return set()
    unmatched = getattr(dynamics, "unmatched_members", {})
    letters: set[str] = set()
    for bridge_class, class_dynamics in dynamics.classes.items():
        if not dynamics.refined.get(bridge_class):
            continue
        if any(member.bridge_id in unmatched for member in class_dynamics.entry.members):
            letters.add(class_dynamics.letter)
    return letters


def class_colors(dynamics: "SymbolicDynamics", *, refined: bool = True) -> dict[str, str]:
    """
    The colour of every (refined) class symbol of a symbolic dynamics.

    Args:
        dynamics: The :class:`~.SymbolicDynamics.SymbolicDynamics`.
        refined: Colour the refined children (``a_1``, ``a_2``) rather than
            their parent letter where a class was split.

    Returns:
        ``{symbol name: colour}`` with :data:`CLASS_COLORS` assigned in
        ``dynamics.classes`` order (refined children in index order).
    """
    return {
        name: CLASS_COLORS[index % len(CLASS_COLORS)]
        for index, (name, _inert) in enumerate(_class_slots(dynamics, refined=refined))
    }


def _member_name(
    dynamics: "SymbolicDynamics", bridge_id, letter: str, *, refined: bool
) -> str:
    """The symbol name a member bridge is drawn under."""
    if refined:
        child = dynamics.member_refinement.get(bridge_id)
        if child is not None:
            return child.name
    return letter


def plot_bridges_by_class(
    trellis: "Trellis",
    dynamics: "SymbolicDynamics",
    ax=None,
    *,
    refined: bool = True,
    show_inert: bool = True,
    **line_kwargs,
):
    """
    Draw every classed bridge, coloured by its (refined) class, with a legend.

    Each member bridge of each class is drawn as its polyline in the class's
    :data:`CLASS_COLORS` slot (:func:`class_colors`); with ``refined`` a
    member matched to a refined child takes the child's colour and name, an
    unmatched member its parent's. Inert classes are dashed. The legend lists
    the symbol names as mathtext. The stable manifold is not drawn; draw the
    tangle underneath first.

    Args:
        trellis: The trellis holding the bridges (``bridge_between``).
        dynamics: The :class:`~.SymbolicDynamics.SymbolicDynamics`.
        ax: Optional matplotlib Axes. Defaults to the current axes.
        refined: Colour and name by refined class where a class was split.
        show_inert: Draw the inert classes too (default True).
        **line_kwargs: Overrides on :data:`CLASS_BRIDGE_STYLE` applied to every
            bridge line (``linewidth``, ``alpha``, ``zorder``, ...).

    Returns:
        ``{symbol name: colour}`` for every symbol drawn, in legend order.

    Raises:
        ValueError: If ``color``, ``linestyle`` or ``label`` is passed: the
            class decides those.
    """
    reserved = _RESERVED_CLASS_KWARGS.intersection(line_kwargs)
    if reserved:
        raise ValueError(
            f"plot_bridges_by_class manages {sorted(reserved)} itself (one colour "
            "and line style per class); pass 'linewidth'/'alpha' instead"
        )
    from matplotlib.lines import Line2D

    target = ax if ax is not None else plt.gca()
    palette = class_colors(dynamics, refined=refined)
    inert_of = dict(_class_slots(dynamics, refined=refined))
    unmatched = _unmatched_parent_slots(dynamics, refined=refined)
    style = dict(CLASS_BRIDGE_STYLE)
    style.update(line_kwargs)

    drawn: dict[str, str] = {}
    for bridge_class, class_dynamics in dynamics.classes.items():
        inert = class_dynamics.kind != "active"
        if inert and not show_inert:
            continue
        for member in class_dynamics.entry.members:
            name = _member_name(dynamics, member.bridge_id, class_dynamics.letter, refined=refined)
            color = palette.get(name, CLASS_UNKNOWN_COLOR)
            bridge = trellis.bridge_between(*member.bridge_id)
            if bridge is None:
                logger.debug("bridge %s has no object to draw; skipping it", member.bridge_id)
                continue
            points = bridge.get_point_array()
            if len(points) < 2:
                continue
            target.plot(
                points[:, 0], points[:, 1], color=color,
                linestyle="--" if inert else "-", **style,
            )
            drawn.setdefault(name, color)
            if name not in palette:
                logger.debug("bridge %s is drawn under symbol %r, which has no palette slot", member.bridge_id, name)

    ordered = {name: drawn[name] for name, _ in _class_slots(dynamics, refined=refined) if name in drawn}
    ordered.update(drawn)
    handles = [
        Line2D(
            [0], [0], color=color,
            linestyle="--" if inert_of.get(name, False) else "-",
            linewidth=style["linewidth"],
            label=name_mathtext(name) + (" (unmatched)" if name in unmatched else ""),
        )
        for name, color in ordered.items()
    ]
    if handles:
        target.legend(handles=handles, loc="best", fontsize=8, framealpha=0.85)
    return ordered


def bridge_class_legend_handles(
    dynamics: "SymbolicDynamics", *, refined: bool = True
) -> list:
    """
    Legend proxies matching :func:`plot_bridges_by_class`, one per symbol.

    Args:
        dynamics: The :class:`~.SymbolicDynamics.SymbolicDynamics`.
        refined: Match the refined symbols where a class was split.

    Returns:
        ``Line2D`` handles labelled with each symbol's mathtext, in palette
        order, dashed for inert classes; a split class's own letter, present
        only when a member matched no child, is labelled ``(unmatched)``.
    """
    from matplotlib.lines import Line2D

    palette = class_colors(dynamics, refined=refined)
    unmatched = _unmatched_parent_slots(dynamics, refined=refined)
    return [
        Line2D(
            [0], [0], color=palette[name], linestyle="--" if inert else "-",
            linewidth=CLASS_BRIDGE_STYLE["linewidth"],
            label=name_mathtext(name) + (" (unmatched)" if name in unmatched else ""),
        )
        for name, inert in _class_slots(dynamics, refined=refined)
    ]


#: Node size (points squared) in :func:`plot_transition_graph`.
TRANSITION_NODE_SIZE = 900

#: Edge curvature in :func:`plot_transition_graph`: a pair of opposite edges
#: ``a -> b`` and ``b -> a`` bows apart, and a self-loop ``a -> a`` shows.
TRANSITION_CONNECTION_STYLE = "arc3,rad=0.15"


#: Radius (points) of a self-loop ring in :func:`plot_transition_graph`.
TRANSITION_LOOP_RADIUS = 11.0


def _draw_self_loop(
    ax, point: NDArray[np.float64], outward: NDArray[np.float64], *,
    color: str, linewidth: float, zorder: float,
) -> None:
    """
    A self-loop as a ring of :data:`TRANSITION_LOOP_RADIUS` points hanging off
    a node on its outward side, with an arrowhead on the ring.

    Drawn in POINT units through an offset transform so it keeps its size
    whatever the data scale; the node (drawn above, white-filled) hides the
    part of the ring that overlaps it.
    """
    from matplotlib.transforms import offset_copy

    node_radius = float(np.sqrt(TRANSITION_NODE_SIZE)) / 2.0
    ring = TRANSITION_LOOP_RADIUS
    distance = node_radius + 0.8 * ring
    norm = float(np.linalg.norm(outward))
    unit = outward / norm if norm > 0.0 else np.array([0.0, 1.0])
    centre = distance * unit
    figure = ax.get_figure()
    ax.scatter(
        [point[0]], [point[1]], s=(2.0 * ring) ** 2, facecolors="none",
        edgecolors=color, linewidths=linewidth, zorder=zorder,
        transform=offset_copy(ax.transData, figure, centre[0], centre[1], units="points"),
    )
    # Arrowhead on the ring at 90° counter-clockwise from the outward look,
    # pointing along the ring's clockwise tangent.
    theta = np.arctan2(unit[1], unit[0]) + 0.5 * np.pi
    tip = centre + ring * np.array([np.cos(theta), np.sin(theta)])
    heading = np.degrees(theta - 0.5 * np.pi)
    ax.scatter(
        [point[0]], [point[1]], s=60, marker=(3, 0, heading - 90.0), color=color,
        zorder=zorder,
        transform=offset_copy(ax.transData, figure, tip[0], tip[1], units="points"),
    )


def plot_transition_graph(
    dynamics: "SymbolicDynamics", ax=None, *, refined: bool = True, **kwargs
):
    """
    Draw the transition graph of a symbolic dynamics.

    Nodes are the (refined) symbols on a circle (``networkx.circular_layout``),
    coloured as in :func:`plot_bridges_by_class` and labelled with their
    mathtext; inert and virtual symbols are sinks drawn with a dashed
    outline. Each edge ``x -> y`` (``y`` occurs in the word of ``x``) is
    curved so that ``a -> b`` and ``b -> a`` bow apart and a self-loop shows,
    with width growing with its ``weight`` (multiplicity), which is also
    written on edges of weight above one. The title reports the spectral
    radius of the transition matrix.

    Args:
        dynamics: The :class:`~.SymbolicDynamics.SymbolicDynamics`.
        ax: Optional matplotlib Axes. Defaults to the current axes.
        refined: Use the refined symbols where a class was split.
        **kwargs: Overrides forwarded to ``networkx.draw_networkx_edges``
            (``edge_color``, ``arrowsize``, ``connectionstyle``, ...).

    Returns:
        The ``networkx.DiGraph`` drawn (as ``dynamics.transition_graph``
        returned it).
    """
    import networkx as nx

    target = ax if ax is not None else plt.gca()
    graph = dynamics.transition_graph(refined=refined)
    palette = class_colors(dynamics, refined=refined)
    radius = dynamics.spectral_radius(refined=refined)
    title = (
        f"transition graph ({'refined' if refined else 'class'} symbols); "
        f"spectral radius {radius:.4g}"
    )
    target.set_title(title, fontsize=10)
    target.set_axis_off()
    if graph.number_of_nodes() == 0:
        logger.info("the transition graph has no nodes; nothing to draw")
        return graph

    pos = nx.circular_layout(graph)
    nodes = list(graph.nodes)
    kinds = {node: graph.nodes[node].get("kind", "active") for node in nodes}
    sources = [node for node in nodes if kinds[node] == "active"]
    sinks = [node for node in nodes if kinds[node] != "active"]
    for group, linestyle in ((sources, "solid"), (sinks, "dashed")):
        if not group:
            continue
        collection = nx.draw_networkx_nodes(
            graph, pos, nodelist=group, ax=target, node_size=TRANSITION_NODE_SIZE,
            node_color="white", edgecolors=[palette.get(n, CLASS_UNKNOWN_COLOR) for n in group],
            linewidths=2.0,
        )
        collection.set_linestyle(linestyle)
        collection.set_zorder(3)
    nx.draw_networkx_labels(
        graph, pos, labels={node: name_mathtext(node) for node in nodes},
        ax=target, font_size=11,
    )

    edge_color = kwargs.pop("edge_color", "#52514e")
    edges = [edge for edge in graph.edges if edge[0] != edge[1]]
    loops = [edge for edge in graph.edges if edge[0] == edge[1]]
    weight_of = {edge: float(graph.edges[edge].get("weight", 1)) for edge in graph.edges}
    width_of = {edge: 0.8 + 0.8 * weight for edge, weight in weight_of.items()}
    if edges:
        edge_kwargs = {
            "edge_color": edge_color,
            "arrows": True,
            "arrowstyle": "-|>",
            "arrowsize": 14,
            "connectionstyle": TRANSITION_CONNECTION_STYLE,
            "node_size": TRANSITION_NODE_SIZE,
            "width": [width_of[edge] for edge in edges],
        }
        edge_kwargs.update(kwargs)
        nx.draw_networkx_edges(graph, pos, edgelist=edges, ax=target, **edge_kwargs)
    centre = np.mean(np.vstack([pos[node] for node in nodes]), axis=0)
    for edge in loops:
        point = np.asarray(pos[edge[0]], dtype=float)
        _draw_self_loop(
            target, point, point - centre, color=edge_color,
            linewidth=width_of[edge], zorder=2,
        )
    heavy = {
        edge: str(int(weight)) for edge, weight in weight_of.items()
        if weight > 1 and edge[0] != edge[1]
    }
    if heavy:
        nx.draw_networkx_edge_labels(
            graph, pos, edge_labels=heavy, ax=target, font_size=8,
            connectionstyle=TRANSITION_CONNECTION_STYLE, label_pos=0.35,
        )
    for edge in loops:
        if weight_of[edge] > 1:
            point = np.asarray(pos[edge[0]], dtype=float)
            outward = point - centre
            norm = float(np.linalg.norm(outward))
            unit = outward / norm if norm > 0.0 else np.array([0.0, 1.0])
            target.annotate(
                str(int(weight_of[edge])), point, textcoords="offset points",
                xytext=tuple(
                    (0.5 * np.sqrt(TRANSITION_NODE_SIZE) + 2.2 * TRANSITION_LOOP_RADIUS) * unit
                ),
                ha="center", va="center", fontsize=8, color=edge_color,
            )
    target.margins(0.15)
    return graph


#: Column headings of :func:`plot_itinerary_table`, in order.
ITINERARY_TABLE_COLUMNS = ("class", "itinerary", "word", "refined word", "status", "verified")

#: Largest font size (points) of :func:`plot_itinerary_table`; shrunk to fit.
ITINERARY_TABLE_FONTSIZE = 8.0

#: Height of a row of :func:`plot_itinerary_table`, in ems of its font.
ITINERARY_TABLE_ROW_HEIGHT = 2.0

#: Room a table cell adds around its text: a factor on the text's width (the
#: cell's own left padding is a tenth of its width) plus a flat gap in ems.
_TABLE_CELL_WIDTH_FACTOR = 1.15
_TABLE_CELL_GAP_EMS = 0.8


def _element_text(
    dynamics: "SymbolicDynamics", ref, *, homotopy: bool, mathtext: bool = False
) -> str:
    """An element's name text or mathtext (label when the dynamics carries no naming)."""
    naming = getattr(dynamics, "naming", None)
    if naming is None:
        return name_mathtext(ref.label) if mathtext else ref.label
    name = naming.homotopy_name(ref) if homotopy else naming.name(ref)
    return name.mathtext if mathtext else name.text


def _itinerary_text(dynamics: "SymbolicDynamics", itinerary, *, mathtext: bool = False) -> str:
    """``R_1^1 R_3^3 | L_3 L_1^2 | ...``: pairs joined by spaces, ``" | "`` between."""
    if itinerary is None:
        return ""
    names = [
        _element_text(dynamics, ref, homotopy=False, mathtext=mathtext) for ref in itinerary
    ]
    pairs = [" ".join(names[i : i + 2]) for i in range(0, len(names), 2)]
    return " | ".join(pairs)


def _word_text(symbols, *, mathtext: bool = False) -> str:
    """A word's symbols joined by spaces, each as text or mathtext; ``""`` when empty."""
    return " ".join(symbol.mathtext if mathtext else symbol.text for symbol in symbols)


def _refined_word_text(
    dynamics: "SymbolicDynamics", bridge_class, letter: str, *, mathtext: bool = False
) -> str:
    """
    The refined word of a class: every refined child of a class shares one
    word, so a split class's is read under its first child's name, an unsplit
    class's under its letter; ``""`` when the dynamics has none.
    """
    children = dynamics.refined.get(bridge_class, [])
    key = min(children, key=lambda c: c.index).name if children else letter
    rules = getattr(dynamics, "refined_rules", None)
    if rules is not None and key in rules:
        return _word_text(rules[key], mathtext=mathtext)
    try:
        return dynamics.word(key, refined=True)
    except KeyError:
        return ""


def _class_cell(
    dynamics: "SymbolicDynamics", class_dynamics, bridge_class, *, mathtext: bool = False
) -> str:
    """``a = {R_1, R_3}``, with the letter and the names as mathtext on request."""
    pair = ", ".join(
        _element_text(dynamics, ref, homotopy=True, mathtext=mathtext)
        for ref in (bridge_class.source, bridge_class.target)
    )
    letter = name_mathtext(class_dynamics.letter) if mathtext else class_dynamics.letter
    return f"{letter} = {{{pair}}}"


def _class_status(class_dynamics) -> str:
    """unique / trivial / ambiguous / unresolved: reason / trellis."""
    if class_dynamics.unresolved_reason:
        return f"unresolved: {class_dynamics.unresolved_reason}"
    if class_dynamics.ambiguous:
        return "ambiguous"
    if class_dynamics.source == "trellis":
        return "trellis"
    if class_dynamics.search is not None:
        return str(class_dynamics.search.status)
    return "-"


def _table_columns(refined: bool, columns: Optional[Sequence[str]]) -> list[str]:
    """
    The columns a table shows: ``columns`` as given (any subset of
    :data:`ITINERARY_TABLE_COLUMNS`, in the caller's order), else every
    column, without "refined word" when ``refined`` is False.

    Raises:
        ValueError: If ``columns`` names a heading that is not a column.
    """
    if columns is None:
        return [c for c in ITINERARY_TABLE_COLUMNS if refined or c != "refined word"]
    chosen = list(columns)
    unknown = [c for c in chosen if c not in ITINERARY_TABLE_COLUMNS]
    if unknown:
        raise ValueError(
            f"unknown itinerary table column(s) {unknown}; "
            f"choose from {list(ITINERARY_TABLE_COLUMNS)}"
        )
    return chosen


def itinerary_table_rows(
    dynamics: "SymbolicDynamics",
    *,
    refined: bool = True,
    mathtext: bool = True,
    columns: Optional[Sequence[str]] = None,
) -> list[list[str]]:
    """
    The rows :func:`plot_itinerary_table` draws, one per class.

    Args:
        dynamics: The :class:`~.SymbolicDynamics.SymbolicDynamics`.
        refined: Include the "refined word" column (when ``columns`` is
            None).
        mathtext: Render every name and symbol as mathtext (``$R_{1}^{1}$``,
            ``$u^{-1}$``) so sub- and superscripts compile; ``False`` gives
            the plain ``R_1^1`` / ``u^-1`` texts.
        columns: The columns to show, a subset of
            :data:`ITINERARY_TABLE_COLUMNS` in the caller's order; None
            shows them all (minus "refined word" when ``refined`` is False).

    Returns:
        Rows of text in the chosen column order: the class as
        ``letter = {X, Y}`` in homotopy names, the itinerary in iterated
        names (``" | "`` between disjoint pairs), the word, the refined
        word, the status and ``yes``/``no``/``-`` for ``verified``.

    Raises:
        ValueError: If ``columns`` names a heading that is not a column.
    """
    chosen = _table_columns(refined, columns)
    rows: list[list[str]] = []
    for bridge_class, class_dynamics in dynamics.classes.items():
        verified = class_dynamics.verified
        cells = {
            "class": _class_cell(dynamics, class_dynamics, bridge_class, mathtext=mathtext),
            "itinerary": _itinerary_text(dynamics, class_dynamics.itinerary, mathtext=mathtext),
            "word": _word_text(class_dynamics.symbols, mathtext=mathtext),
            "refined word": _refined_word_text(
                dynamics, bridge_class, class_dynamics.letter, mathtext=mathtext
            ),
            "status": _class_status(class_dynamics),
            "verified": "-" if verified is None else ("yes" if verified else "no"),
        }
        rows.append([cells[column] for column in chosen])
    return rows


def _text_width_ems(renderer, text: str, *, bold: bool = False) -> float:
    """
    The rendered width of ``text`` (mathtext compiled) in ems of its font,
    measured on ``renderer``; ``text``'s width grows linearly with the font
    size, so one measurement serves every size.
    """
    from matplotlib.cbook import is_math_text
    from matplotlib.font_manager import FontProperties

    reference = 10.0
    prop = FontProperties(size=reference, weight="bold" if bold else "normal")
    width, _height, _descent = renderer.get_text_width_height_descent(
        text, prop, ismath=is_math_text(text)
    )
    return width / renderer.points_to_pixels(reference)


def _table_column_ems(
    target, columns: Sequence[str], rows: Sequence[Sequence[str]]
) -> list[float]:
    """
    The width each column needs, in ems of the table's font: the widest of
    its header and its cells (rendered widths when the figure has a
    renderer, ~0.6 em per glyph otherwise), padded as a cell pads.
    """
    renderer = None
    try:
        renderer = target.get_figure().canvas.get_renderer()
    except (AttributeError, ValueError):  # a backend with no renderer
        pass
    widths: list[float] = []
    for index, column in enumerate(columns):
        texts = [(column, True)] + [(row[index], False) for row in rows]
        if renderer is not None:
            longest = max(_text_width_ems(renderer, text, bold=bold) for text, bold in texts)
        else:
            longest = 0.6 * max(len(text) for text, _bold in texts)
        widths.append(_TABLE_CELL_WIDTH_FACTOR * longest + _TABLE_CELL_GAP_EMS)
    return widths


def plot_itinerary_table(
    dynamics: "SymbolicDynamics",
    ax=None,
    *,
    refined: bool = True,
    mathtext: bool = True,
    columns: Optional[Sequence[str]] = None,
    fontsize: Optional[float] = None,
    row_height: float = ITINERARY_TABLE_ROW_HEIGHT,
):
    """
    Tabulate a symbolic dynamics: one row per class, axes turned off.

    Columns (:data:`ITINERARY_TABLE_COLUMNS`): the class (``a = {R_1, R_3}``),
    its itinerary in iterated-element names with ``" | "`` between the
    disjoint pairs, its word, its refined word (omitted when ``refined`` is
    False), its status (``unique`` / ``trivial`` / ``ambiguous`` /
    ``unresolved: <reason>`` / ``trellis``) and whether the registered
    images verified the word (``yes`` / ``no`` / ``-``); ``columns`` picks
    a subset.

    The table fills the axes' width: each column is as wide as its widest
    rendered text needs, and the font is ``fontsize`` unless the widest row
    would not fit at that size, in which case it is shrunk until it does.
    The rows grow with the font (``row_height`` ems each).

    Args:
        dynamics: The :class:`~.SymbolicDynamics.SymbolicDynamics`.
        ax: Optional matplotlib Axes. Defaults to the current axes.
        refined: Include the "refined word" column (when ``columns`` is None).
        mathtext: Compile the names and symbols (``$R_{1}^{1}$``, ``$u^{-1}$``)
            rather than printing ``R_1^1`` / ``u^-1``.
        columns: The columns to show, a subset of
            :data:`ITINERARY_TABLE_COLUMNS` in the caller's order; None
            shows them all.
        fontsize: The largest font size (points); defaults to
            :data:`ITINERARY_TABLE_FONTSIZE`.
        row_height: The height of a row in ems of the (final) font.

    Returns:
        The matplotlib ``Table``, or None when the dynamics has no class.

    Raises:
        ValueError: If ``columns`` names a heading that is not a column.
    """
    target = ax if ax is not None else plt.gca()
    target.set_axis_off()
    chosen = _table_columns(refined, columns)
    rows = itinerary_table_rows(dynamics, refined=refined, mathtext=mathtext, columns=chosen)
    if not rows:
        logger.info("the symbolic dynamics has no class; nothing to tabulate")
        return None
    # Column widths as fractions of the axes, proportional to the width each
    # column's widest text needs, so the table fills the axes exactly rather
    # than overflowing it as matplotlib's auto width does.
    column_ems = _table_column_ems(target, chosen, rows)
    needed_ems = float(sum(column_ems))
    widths = [ems / needed_ems for ems in column_ems]
    # Cells never wrap, so cap the font where the widest row fits the axes.
    size = ITINERARY_TABLE_FONTSIZE if fontsize is None else float(fontsize)
    try:
        figure = target.get_figure()
        axes_points = target.get_window_extent().width / figure.dpi * 72.0
        if needed_ems * size > axes_points > 0.0:
            size = max(4.0, axes_points / needed_ems)
    except (AttributeError, ValueError):  # a figure with no renderer yet
        pass
    table = target.table(
        cellText=rows, colLabels=chosen, colWidths=widths,
        loc="center", cellLoc="left", colLoc="left",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(size)
    # matplotlib sizes the rows for its default font (``Table.FONTSIZE``,
    # 1.2 em tall) whatever font the cells carry; rescale them to
    # ``row_height`` ems of the font actually used.
    table.scale(1.0, row_height * size / (1.2 * table.FONTSIZE))
    for (row, _col), cell in table.get_celld().items():
        cell.set_edgecolor("#e1e0d9")
        if row == 0:
            cell.set_text_props(fontweight="bold")
            cell.set_facecolor("#f0efec")
    return table


# ── dual-graph cartoon ──────────────────────────────────────────────────────

#: Offset (data units) of an element's bar and brackets from its stable line.
CARTOON_SEGMENT_OFFSET = 0.12

#: Offset of an element's node from its stable line (outside the bar).
CARTOON_NODE_OFFSET = 0.32

#: Offset of an element's name from its stable line (outside the node) when
#: the names are placed ``"outside"``.
CARTOON_LABEL_OFFSET = 0.62

#: Gap (points) between a node and its element's name placed ``"beside"`` it.
CARTOON_LABEL_GAP = 7.0

#: Font size (points) of the element names.
CARTOON_LABEL_FONTSIZE = 9.0

#: An arc between two nodes ``dx`` apart bulges ``BASE + SLOPE * |dx|`` away
#: from the line, so a bridge spanning more elements nests OUTSIDE a narrower
#: one, as on a chalkboard.
CARTOON_ARC_BASE = 0.45
CARTOON_ARC_SLOPE = 0.22

#: Half-width of the teardrop drawn for a loop (both ends in one element).
CARTOON_LOOP_HALFWIDTH = 0.3

#: Pull of the control points of the S-curve crossing the line at a unified node.
CARTOON_S_PULL = 0.5

#: Vertical gap between the rows of two partitioned stable branches.
CARTOON_ROW_GAP = 1.0

#: Marker area (points squared) of an element node.
CARTOON_NODE_SIZE = 60

#: Line widths of successive walks: later walks are drawn thinner ON TOP of
#: the earlier ones, since walks of different classes share arcs; all of them
#: are thinner than the translucent bridges they run along, so a bridge's
#: colour shows either side of the walk on it.
CARTOON_WALK_LINEWIDTHS = (2.6, 1.8, 1.2, 0.9)

#: Style of the stable lines, the unified crossings and the element bars.
CARTOON_LINE_STYLE = {"color": "#52514e", "linewidth": 1.5, "zorder": 1}
CARTOON_UNIFIED_STYLE = {"color": "#898781", "linewidth": 1.0, "zorder": 3}
CARTOON_SIDE_COLORS = {"left": "tab:blue", "right": "tab:red"}

#: Line widths of a hole bridge and of an image bridge in the cartoon, and
#: the bridges' opacity (a walk runs on top of them).
CARTOON_HOLE_LINEWIDTH = 6.0
CARTOON_IMAGE_LINEWIDTH = 3.5
CARTOON_BRIDGE_ALPHA = 0.45


@dataclass
class CartoonLayout:
    """
    Where :func:`plot_dual_graph_cartoon` put everything.

    Attributes:
        rows: ``{branch_key: y}`` of each partitioned stable branch's line;
            the first branch of the partition family is at 0, the next ones
            below it.
        ranks: ``{branch_key: [cdist, ...]}``, the representative boundary
            canonical distances of a branch in ascending order; the ordinal
            coordinate of a cdist is its index here.
        tol: The tolerance two boundary cdists are identified within.
        anchor: ``"right"`` (the default: the anchor at the right end of
            every line, elements running outward to the left, as on the
            chalkboard) or ``"left"``.
        segments: ``{element ref: (x_lo, x_hi, y)}`` of each element's bar:
            the positions of its anchorward and outward ends (``x_lo > x_hi``
            with the anchor on the right).
        nodes: ``{element ref: (x, y)}``, ONE node per element per side; a
            singleton's node is its dot on the bar.
        unified: ``{element ref: bool}``, whether any stable node of the
            element is unified (traversable).
        singletons: The elements ``[x, x]`` (drawn as a dot, no node circle:
            they own no stable edge, hence no dual-graph node).
        bridges_drawn: The number of minimal-trellis bridges drawn.
        walks_drawn: The number of class itineraries drawn.
    """

    rows: dict = field(default_factory=dict)
    ranks: dict = field(default_factory=dict)
    tol: float = 0.0
    anchor: str = "right"
    segments: dict = field(default_factory=dict)
    nodes: dict = field(default_factory=dict)
    unified: dict = field(default_factory=dict)
    singletons: set = field(default_factory=set)
    bridges_drawn: int = 0
    walks_drawn: int = 0

    def rank(self, branch_key, cdist: float) -> float:
        """The ordinal coordinate of a boundary cdist on a branch (see :attr:`ranks`)."""
        return _rank_of(self.ranks[branch_key], float(cdist), self.tol)

    def position(self, branch_key, cdist: float) -> float:
        """The x of a boundary cdist: its rank, mirrored (``width - rank``) with the anchor on the right."""
        rank = self.rank(branch_key, cdist)
        return rank if self.anchor == "left" else self.width - rank

    @property
    def width(self) -> float:
        """The largest ordinal coordinate over every row."""
        return max((len(reps) - 1 for reps in self.ranks.values()), default=0)

    def line_extent(self, branch_key) -> tuple[float, float]:
        """The x range of a branch's line, from its anchor to its outermost boundary."""
        span = max(len(self.ranks[branch_key]) - 1, 1)
        if self.anchor == "left":
            return 0.0, float(span)
        return self.width - span, self.width


def _rank_of(representatives: Sequence[float], cdist: float, tol: float) -> float:
    """The index of the representative within ``tol`` of ``cdist`` (else the nearest)."""
    import bisect

    index = bisect.bisect_left(representatives, cdist)
    candidates = [i for i in (index - 1, index) if 0 <= i < len(representatives)]
    best = min(candidates, key=lambda i: abs(representatives[i] - cdist))
    if abs(representatives[best] - cdist) > tol:
        logger.debug(
            "cdist %.6g is no boundary representative (nearest %.6g); using it anyway",
            cdist,
            representatives[best],
        )
    return float(best)


def _cluster(values: Iterable[float], tol: float) -> list[float]:
    """Sorted representatives of ``values``, values within ``tol`` merged."""
    representatives: list[float] = []
    for value in sorted(float(v) for v in values):
        if representatives and value - representatives[-1] <= tol:
            continue
        representatives.append(value)
    return representatives


def _side_sign(side) -> float:
    """``+1`` for the left side (drawn above the line), ``-1`` for the right."""
    return 1.0 if side == "left" else -1.0


def _row_label(branch_key) -> str:
    """``p3@1.0``: the fixed point's period, orbit index and branch index."""
    fixed_point, _stability, orbit, branch = branch_key
    return f"p{fixed_point.period}@{orbit}.{branch}"


def dual_graph_cartoon_layout(
    dual_graph: "DualGraph", *, anchor: Literal["left", "right"] = "right"
) -> CartoonLayout:
    """
    Lay out the cartoon of a dual graph: rows, ordinal coordinates, bars, nodes.

    Every boundary canonical distance of BOTH sides of a branch (the ends of
    every element of ``dual_graph.partition``, identified within the family's
    tolerance) gets consecutive ordinal coordinates, so a boundary shared by
    the two sides aligns and every non-singleton element is at least one unit
    wide whatever its real length. A singleton element ``[x, x]`` sits at its
    boundary's coordinate with zero width. With the anchor on the RIGHT (the
    default) the coordinate is mirrored, so every row's anchor sits at the
    same right end and the elements run outward to the left. The left side's
    bars, nodes and names go ABOVE the line, the right side's BELOW; each
    element gets ONE node at the middle of its bar.

    Args:
        dual_graph: The :class:`~.DualGraph.DualGraph` (its ``partition`` is
            the family drawn, normally the iterated homotopy partition).
        anchor: Which end of the lines the anchor is at.

    Returns:
        The :class:`CartoonLayout`.

    Raises:
        ValueError: For an ``anchor`` other than ``"left"`` / ``"right"``.
    """
    if anchor not in ("left", "right"):
        raise ValueError(f"anchor must be 'left' or 'right', got {anchor!r}")
    partition = dual_graph.partition
    layout = CartoonLayout(tol=max(float(partition.tol), 1e-12), anchor=anchor)
    for branch_key in partition.branch_keys:
        boundaries = [
            cdist
            for side in partition.sides(branch_key)
            for iv in partition.result(branch_key, side).intervals
            for cdist in (iv.lo_cdist, iv.hi_cdist)
        ]
        layout.ranks[branch_key] = _cluster(boundaries, layout.tol)
    y = 0.0
    for branch_key in partition.branch_keys:
        sides = partition.sides(branch_key)
        layout.rows[branch_key] = y
        for side in sides:
            sigma = _side_sign(side)
            result = partition.result(branch_key, side)
            for iv in result.intervals:
                ref = result.ref(iv.element_id)
                x_lo = layout.position(branch_key, iv.lo_cdist)
                x_hi = layout.position(branch_key, iv.hi_cdist)
                layout.segments[ref] = (x_lo, x_hi, y + sigma * CARTOON_SEGMENT_OFFSET)
                if x_lo == x_hi:
                    layout.singletons.add(ref)
                    layout.nodes[ref] = (x_lo, y + sigma * CARTOON_SEGMENT_OFFSET)
                else:
                    layout.nodes[ref] = (0.5 * (x_lo + x_hi), y + sigma * CARTOON_NODE_OFFSET)
                layout.unified[ref] = any(
                    node.is_unified for node in dual_graph.stable_nodes_of(ref)
                )
        # The next branch's row goes below the deepest arc this one can carry.
        deepest = CARTOON_ARC_BASE + CARTOON_ARC_SLOPE * max(len(layout.ranks[branch_key]) - 1, 1)
        y -= 2.0 * (CARTOON_NODE_OFFSET + deepest) + CARTOON_ROW_GAP
    return layout


def _u_arc(p, q, sigma_p: float, sigma_q: float):
    """
    A wide smooth arc from node ``p`` to node ``q``, leaving and arriving
    perpendicular to the line, bulging ``CARTOON_ARC_BASE + SLOPE * |dx|``
    away from it; a teardrop when ``p == q``.
    """
    from matplotlib.path import Path

    depth = CARTOON_ARC_BASE + CARTOON_ARC_SLOPE * abs(q[0] - p[0])
    height = 4.0 * depth / 3.0  # a cubic's apex is 3/4 of its control height
    if p == q:
        vertices = [
            p,
            (p[0] - CARTOON_LOOP_HALFWIDTH, p[1] + sigma_p * height),
            (p[0] + CARTOON_LOOP_HALFWIDTH, p[1] + sigma_p * height),
            p,
        ]
    else:
        vertices = [
            p,
            (p[0], p[1] + sigma_p * height),
            (q[0], q[1] + sigma_q * height),
            q,
        ]
    return Path(vertices, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])


def _s_curve(a, b):
    """A smooth S-curve from node ``a`` on one side of a line to node ``b`` on the other."""
    from matplotlib.path import Path

    toward = -1.0 if a[1] > b[1] else 1.0
    vertices = [
        a,
        (a[0], a[1] + toward * CARTOON_S_PULL),
        (b[0], b[1] - toward * CARTOON_S_PULL),
        b,
    ]
    return Path(vertices, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])


def _cubic_point(path, t: float):
    """The point at parameter ``t`` of a single-cubic path."""
    p0, p1, p2, p3 = (np.asarray(v, dtype=float) for v in path.vertices[:4])
    s = 1.0 - t
    return s**3 * p0 + 3 * s**2 * t * p1 + 3 * s * t**2 * p2 + t**3 * p3


def _add_curve(ax, path, *, gid: str, **style):
    """Draw a path as an unfilled patch; return the patch."""
    from matplotlib.patches import PathPatch

    patch = PathPatch(path, fill=False, gid=gid, capstyle="round", **style)
    ax.add_patch(patch)
    return patch


def _cartoon_name(dual_graph: "DualGraph", ref) -> str:
    """An element's name as mathtext (its label when the graph has no naming)."""
    naming = dual_graph.naming
    if naming is None:
        return name_mathtext(ref.label)
    if getattr(dual_graph.partition, "homotopy", None) is not None:
        return naming.name(ref).mathtext
    return naming.homotopy_name(ref).mathtext


def _bridge_symbol(dynamics: "SymbolicDynamics", bridge_id) -> tuple[Optional[str], bool]:
    """``(symbol name, inert)`` a bridge is coloured under; ``(None, False)`` when unclassed."""
    child = dynamics.member_refinement.get(bridge_id)
    if child is not None:
        class_dynamics = dynamics.classes.get(child.parent)
        return child.name, class_dynamics is not None and class_dynamics.kind != "active"
    try:
        entry = dynamics.table.entry_of(bridge_id)
    except KeyError:
        return None, False
    class_dynamics = dynamics.classes.get(entry.bridge_class)
    if class_dynamics is None:
        return None, False
    return class_dynamics.letter, class_dynamics.kind != "active"


def _bridge_end_elements(dual_graph: "DualGraph", bridge_id) -> tuple:
    """The iterated elements a bridge's two crossings are owned by on the bridge's rows."""
    from .StablePartition import row_of_end

    return tuple(
        dual_graph.partition.owner_of_intersection(
            bridge_id[index], row_of_end(dual_graph.trellis, bridge_id, endpoint)
        )
        for index, endpoint in ((0, "first"), (1, "second"))
    )


def _cartoon_itinerary(class_dynamics) -> Optional[tuple]:
    """The element-to-element itinerary a class is drawn from, or None when trivial."""
    itinerary = class_dynamics.itinerary
    if itinerary is None or len(itinerary) < 2 or len(itinerary) % 2:
        return None
    if len(itinerary) == 2 and itinerary[0] == itinerary[1]:
        return None
    return tuple(itinerary)


def plot_dual_graph_cartoon(
    dual_graph: "DualGraph",
    dynamics: "SymbolicDynamics",
    ax=None,
    *,
    anchor: Literal["left", "right"] = "right",
    show_bridges: bool = True,
    show_unified: bool = True,
    show_walks: bool = True,
    color_bridges: bool = False,
    bridge_color: str = CLASS_UNKNOWN_COLOR,
    bridge_alpha: float = CARTOON_BRIDGE_ALPHA,
    bridge_linewidth: Optional[float] = None,
    walk_linewidths: Sequence[float] = CARTOON_WALK_LINEWIDTHS,
    label_elements: bool = True,
    label_fontsize: float = CARTOON_LABEL_FONTSIZE,
    label_position: Literal["beside", "outside"] = "beside",
) -> CartoonLayout:
    """
    Draw the dual graph as a chalkboard cartoon rather than over the plane.

    Each partitioned stable branch is a horizontal line with its iterated
    elements as bars along it in an ORDINAL coordinate (see
    :func:`dual_graph_cartoon_layout`), the left side above, the right side
    below, bracketed like :func:`plot_stable_partition` (``[ ]`` closed,
    ``( )`` open, a dot for a singleton) and named. Each element has ONE
    node: an open circle when every stable node of the element is a wall, a
    filled one when the element carries a unified (traversable) node. Every
    unified stable node is a short S-curve crossing the line between its two
    side elements; every bridge of the minimal trellis is a wide U-arc on its
    row between the elements its crossings are owned by (deeper the wider,
    so bridges nest), in ``bridge_color`` (grey like the crossings by
    default; hole bridges thick, image bridges thin, inert classes dashed;
    ``color_bridges`` colours them by refined class like
    :func:`plot_bridges_by_class` instead); and every
    class itinerary is drawn element to element — a U-arc per same-side
    pair, an S-curve per crossing — in the class colour, from a dot at its
    landing element to an arrowhead at its target element (a
    trellis-resolved itinerary dotted). Anchor on the right by default. An
    element's name sits just to the RIGHT of its node by default
    (``label_position="beside"``): the bridges and walks leave a node away
    from the line, so nothing runs through a name there; ``"outside"`` puts
    it beyond the node, away from the line, where a nested arc may cross it.

    Args:
        dual_graph: The :class:`~.DualGraph.DualGraph`.
        dynamics: Its :class:`~.SymbolicDynamics.SymbolicDynamics` (classes,
            refinement, itineraries).
        ax: Optional matplotlib Axes. Defaults to the current axes.
        anchor: Which end of the lines the anchor is at (see
            :func:`dual_graph_cartoon_layout`).
        show_bridges: Draw the minimal trellis's bridges.
        show_unified: Draw the unified crossings.
        show_walks: Draw the class itineraries.
        color_bridges: Colour the bridges by refined class rather than
            ``bridge_color``.
        bridge_color: The colour of the bridges when ``color_bridges`` is
            False.
        bridge_alpha: The bridges' opacity (walks are drawn on top of them).
        bridge_linewidth: One line width for every bridge; None (the
            default) draws hole bridges :data:`CARTOON_HOLE_LINEWIDTH` and
            image bridges :data:`CARTOON_IMAGE_LINEWIDTH` wide.
        walk_linewidths: Line widths of successive itineraries (the last one
            repeats), later ones drawn thinner on top.
        label_elements: Annotate every element with its name.
        label_fontsize: Font size (points) of the element names.
        label_position: ``"beside"`` (to the right of the node) or
            ``"outside"`` (beyond the node, away from the line).

    Returns:
        The :class:`CartoonLayout`, with ``bridges_drawn`` / ``walks_drawn``
        filled in.
    """
    target = ax if ax is not None else plt.gca()
    layout = dual_graph_cartoon_layout(dual_graph, anchor=anchor)
    partition = dual_graph.partition
    several = len(layout.rows) > 1
    mirrored = layout.anchor == "right"
    y_min, y_max = 0.0, 0.0

    def track(*points) -> None:
        nonlocal y_min, y_max
        for _x, y in points:
            y_min, y_max = min(y_min, y), max(y_max, y)

    for branch_key, y in layout.rows.items():
        x_from, x_to = layout.line_extent(branch_key)
        target.plot(
            [x_from - 0.3, x_to + 0.3], [y, y], solid_capstyle="round", **CARTOON_LINE_STYLE
        )
        anchor_x = x_to + 0.3 if mirrored else x_from - 0.3
        target.annotate(
            "anchor" if not several else f"anchor  {_row_label(branch_key)}",
            (anchor_x, y), textcoords="offset points", xytext=(4 if mirrored else -4, 0),
            ha="left" if mirrored else "right", va="center", fontsize=8, color="gray",
        )
        track((0.0, y))
        for side in partition.sides(branch_key):
            color = CARTOON_SIDE_COLORS[side]
            sigma = _side_sign(side)
            result = partition.result(branch_key, side)
            for iv in result.intervals:
                ref = result.ref(iv.element_id)
                x_lo, x_hi, y_bar = layout.segments[ref]
                if x_lo == x_hi:
                    target.plot(
                        x_lo, y_bar, marker="o", markersize=7, markerfacecolor=color,
                        markeredgecolor=color, linestyle="none", zorder=7,
                    )
                else:
                    target.plot(
                        [x_lo, x_hi], [y_bar, y_bar], color=color, linewidth=2,
                        solid_capstyle="butt", zorder=5,
                    )
                    # A bracket opens toward the element's interior: the
                    # anchorward end is on the right when mirrored.
                    for x, closed, at_right in (
                        (x_lo, iv.closed_lo, mirrored), (x_hi, iv.closed_hi, not mirrored)
                    ):
                        if at_right:
                            glyph, nudge, align = ("]" if closed else ")"), -1, "right"
                        else:
                            glyph, nudge, align = ("[" if closed else "("), 1, "left"
                        target.annotate(
                            glyph, (x, y_bar), textcoords="offset points", xytext=(nudge, 0),
                            ha=align, va="center", fontsize=13, fontweight="bold",
                            color=color, zorder=6,
                        )
                if label_elements:
                    name = _cartoon_name(dual_graph, ref)
                    label_style = {"fontsize": label_fontsize, "color": "black", "zorder": 9}
                    if label_position == "beside":
                        # Nothing crosses a name beside its node (arcs leave
                        # the node away from the line), so no backing box,
                        # which would wash out the bar's brackets.
                        target.annotate(
                            name, layout.nodes[ref], textcoords="offset points",
                            xytext=(0.5 * np.sqrt(CARTOON_NODE_SIZE) + CARTOON_LABEL_GAP, 0),
                            ha="left", va="center", **label_style,
                        )
                    elif label_position == "outside":
                        target.annotate(
                            name, (0.5 * (x_lo + x_hi), y + sigma * CARTOON_LABEL_OFFSET),
                            ha="center", va="bottom" if sigma > 0 else "top",
                            bbox={"boxstyle": "round,pad=0.15", "facecolor": "white",
                                  "edgecolor": "none", "alpha": 0.85},
                            **label_style,
                        )
                    else:
                        raise ValueError(
                            f"label_position must be 'beside' or 'outside', not {label_position!r}"
                        )
                track((0.0, y + sigma * CARTOON_LABEL_OFFSET))

    # Nodes: one per element per side, open = wall, filled = unified.
    refs = [ref for ref in layout.nodes if ref not in layout.singletons]
    xy = np.array([layout.nodes[ref] for ref in refs]).reshape(-1, 2)
    filled = np.array([layout.unified[ref] for ref in refs], dtype=bool)
    target.scatter(
        xy[:, 0], xy[:, 1], s=CARTOON_NODE_SIZE,
        facecolors=["black" if f else "white" for f in filled],
        edgecolors="black", linewidths=1.2, zorder=8,
    )

    if show_unified:
        pairs = {
            (node.elements["left"], node.elements["right"])
            for node in dual_graph.unified_nodes
            if "left" in node.elements and "right" in node.elements
        }
        for left_ref, right_ref in sorted(pairs, key=lambda pair: layout.nodes[pair[0]]):
            _add_curve(
                target, _s_curve(layout.nodes[left_ref], layout.nodes[right_ref]),
                gid="unified", **CARTOON_UNIFIED_STYLE,
            )

    if show_bridges:
        palette = class_colors(dynamics, refined=True) if color_bridges else {}
        minimal = dual_graph.minimal
        for bridge_id in minimal.kept_bridge_ids:
            try:
                first, second = _bridge_end_elements(dual_graph, bridge_id)
            except (KeyError, ValueError) as error:
                logger.warning("bridge %s has no owned end elements (%s); not drawn", bridge_id, error)
                continue
            if first not in layout.nodes or second not in layout.nodes:
                logger.warning("bridge %s ends outside the drawn partition; not drawn", bridge_id)
                continue
            name, inert = _bridge_symbol(dynamics, bridge_id)
            color = palette.get(name, CLASS_UNKNOWN_COLOR) if name and color_bridges else bridge_color
            path = _u_arc(
                layout.nodes[first], layout.nodes[second],
                _side_sign(first.side), _side_sign(second.side),
            )
            _add_curve(
                target, path, gid=f"bridge:{bridge_id[0]}-{bridge_id[1]}",
                edgecolor=color, linestyle="--" if inert else "-",
                linewidth=(
                    bridge_linewidth if bridge_linewidth is not None
                    else CARTOON_HOLE_LINEWIDTH if minimal.is_hole_bridge(bridge_id)
                    else CARTOON_IMAGE_LINEWIDTH
                ),
                alpha=bridge_alpha,
                zorder=4,
            )
            track(tuple(_cubic_point(path, 0.5)))
            layout.bridges_drawn += 1

    if show_walks:
        palette = class_colors(dynamics, refined=False)
        drawn = 0
        for class_dynamics in dynamics.classes.values():
            itinerary = _cartoon_itinerary(class_dynamics)
            if itinerary is None or any(ref not in layout.nodes for ref in itinerary):
                continue
            letter = class_dynamics.letter
            color = palette.get(letter, CLASS_UNKNOWN_COLOR)
            width = walk_linewidths[min(drawn, len(walk_linewidths) - 1)]
            style = {
                "edgecolor": color, "linewidth": width, "zorder": 10 + drawn,
                "linestyle": ":" if class_dynamics.source == "trellis" else "-",
            }
            gid = f"walk:{letter}"
            last = None
            for index in range(0, len(itinerary), 2):
                a, b = itinerary[index], itinerary[index + 1]
                last = _u_arc(
                    layout.nodes[a], layout.nodes[b], _side_sign(a.side), _side_sign(b.side)
                )
                _add_curve(target, last, gid=gid, **style)
                track(tuple(_cubic_point(last, 0.5)))
                if index + 2 < len(itinerary):
                    c = itinerary[index + 2]
                    last = _s_curve(layout.nodes[b], layout.nodes[c])
                    _add_curve(target, last, gid=gid, **style)
            start = layout.nodes[itinerary[0]]
            target.plot(
                start[0], start[1], marker="o", markersize=WALK_MARKER_SIZE, color=color,
                linestyle="none", zorder=style["zorder"] + 0.5,
            )
            end = layout.nodes[itinerary[-1]]
            target.annotate(
                "", xy=end, xytext=tuple(_cubic_point(last, 0.92)),
                arrowprops={
                    "arrowstyle": "-|>", "color": color, "linewidth": 0,
                    "mutation_scale": 2.0 * WALK_MARKER_SIZE, "shrinkA": 0, "shrinkB": 0,
                },
                zorder=style["zorder"] + 0.5,
            )
            drawn += 1
        layout.walks_drawn = drawn

    target.set_aspect("auto")
    target.set_axis_off()
    if mirrored:
        target.set_xlim(-0.8, layout.width + 1.2)
    else:
        target.set_xlim(-1.2, layout.width + 0.8)
    target.set_ylim(y_min - 0.6, y_max + 0.6)
    return layout


def dual_graph_cartoon_legend_handles(
    dynamics: Optional["SymbolicDynamics"] = None,
    *,
    bridge_color: str = CLASS_UNKNOWN_COLOR,
    bridge_alpha: float = CARTOON_BRIDGE_ALPHA,
    bridge_linewidth: Optional[float] = None,
    walks: bool = True,
) -> list:
    """
    Legend proxies for :func:`plot_dual_graph_cartoon`.

    Args:
        dynamics: When given (and ``walks``), one ``walk of $x$`` handle per
            class whose itinerary the cartoon draws, in the class colour and
            with the successive :data:`CARTOON_WALK_LINEWIDTHS`, follows the
            fixed entries.
        bridge_color: The bridges' colour, as passed to the plotter.
        bridge_alpha: The bridges' opacity, as passed to the plotter.
        bridge_linewidth: The uniform bridge width, as passed to the
            plotter; when given, hole and image bridges look alike, so ONE
            ``bridge`` entry replaces the ``hole bridge`` / ``image bridge``
            pair.
        walks: Include the walk entries (the trellis-itinerary line style
            and the per-class walks); False for a cartoon drawn with
            ``show_walks=False``.

    Returns:
        ``Line2D`` handles: hole bridge, image bridge (or one ``bridge``
        entry under a uniform ``bridge_linewidth``), inert class, unified
        crossing, open node, unified node, then (with ``walks``) trellis
        itinerary (dotted) and one walk per drawn class (a start dot, dotted
        when trellis-resolved).
    """
    from matplotlib.lines import Line2D

    if bridge_linewidth is None:
        bridges = [
            Line2D([0], [0], color=bridge_color, alpha=bridge_alpha,
                   linewidth=CARTOON_HOLE_LINEWIDTH, label="hole bridge"),
            Line2D([0], [0], color=bridge_color, alpha=bridge_alpha,
                   linewidth=CARTOON_IMAGE_LINEWIDTH, label="image bridge"),
        ]
        inert_width = CARTOON_IMAGE_LINEWIDTH
    else:
        bridges = [
            Line2D([0], [0], color=bridge_color, alpha=bridge_alpha,
                   linewidth=bridge_linewidth, label="bridge"),
        ]
        inert_width = bridge_linewidth
    handles = bridges + [
        Line2D([0], [0], color=bridge_color, alpha=bridge_alpha,
               linewidth=inert_width, linestyle="--", label="inert class"),
        Line2D([0], [0], label="unified crossing", **{
            k: v for k, v in CARTOON_UNIFIED_STYLE.items() if k != "zorder"
        }),
        Line2D([0], [0], marker="o", markerfacecolor="white", markeredgecolor="black",
               linestyle="none", markersize=np.sqrt(CARTOON_NODE_SIZE), label="open node (wall)"),
        Line2D([0], [0], marker="o", markerfacecolor="black", markeredgecolor="black",
               linestyle="none", markersize=np.sqrt(CARTOON_NODE_SIZE), label="unified node"),
    ]
    if not walks:
        return handles
    handles.append(
        Line2D([0], [0], color="black", linewidth=CARTOON_WALK_LINEWIDTHS[0],
               linestyle=":", label="trellis itinerary (no walk)")
    )
    if dynamics is not None:
        palette = class_colors(dynamics, refined=False)
        drawn = 0
        for class_dynamics in dynamics.classes.values():
            if _cartoon_itinerary(class_dynamics) is None:
                continue
            letter = class_dynamics.letter
            width = CARTOON_WALK_LINEWIDTHS[min(drawn, len(CARTOON_WALK_LINEWIDTHS) - 1)]
            handles.append(
                Line2D(
                    [0], [0], color=palette.get(letter, CLASS_UNKNOWN_COLOR), linewidth=width,
                    linestyle=":" if class_dynamics.source == "trellis" else "-",
                    marker="o", markersize=WALK_MARKER_SIZE,
                    label=f"walk of {name_mathtext(letter)}",
                )
            )
            drawn += 1
    return handles
