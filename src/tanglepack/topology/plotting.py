"""
Drawing routines for the topological layer.

Every ``plot_*`` of the topological layer lives here so the marker, colour and
z-order conventions are defined once; :class:`~.Trellis.Trellis` and
:class:`~..loom.TangleSession.TangleSession` keep thin delegating wrappers.

Dev Notes — topological plotting.

The z-orders are a stack, not arbitrary numbers: black intersections (drawn by
the numerical layer) sit lowest, then the dual graph's own edges (z=1) and
face points (z=2) — kept low so a dual graph drawn over a trellis plot never
hides the numerics dots underneath it — then its stable nodes (z=10), then the
magenta strong-pip CANDIDATES (z=11), then the green chosen pip (z=12) on top
of the set it was chosen from, then the orange pseudoneighbours (z=13), then
the holes (z=14) and their labels (z=15). Changing one means checking the
whole ladder.

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
from typing import Iterable, Optional, Sequence, TYPE_CHECKING, Union

import numpy as np
import matplotlib.pyplot as plt
from numpy.typing import NDArray

from .TopologyResults import StablePartitionResult

if TYPE_CHECKING:
    from .DualGraph import DualGraph, FaceNode, StableNode
    from .MinimalTrellis import MinimalTrellis
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


def plot_stable_partition(
    results: Union[StablePartitionResult, Iterable[StablePartitionResult]],
    ax=None,
    *,
    labels: Optional[Sequence[str]] = None,
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


def plot_dual_graph(
    dual_graph: "DualGraph",
    ax=None,
    *,
    show_labels: bool = False,
    clip_to_arcs: bool = True,
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
        show_labels: Annotate each stable node with its element label(s) and
            each face node with ``index:kind``.
        clip_to_arcs: Limit the axes to the padded bounding box of the edge
            midpoints and the unbounded face point (see
            :func:`_clip_axes_to_arcs`). Pass False to keep matplotlib's
            autoscale, e.g. when drawing over a tangle whose view is set
            elsewhere.
        **scatter_kwargs: Overrides for the stable-node scatters on top of
            :data:`DUAL_GRAPH_ARC_STYLE` — ``color`` (the solid fill and,
            unless ``edgecolors`` is given, the open ring), ``edgecolors``,
            ``s``, ``zorder``, ...

    Returns:
        The Axes drawn on.

    Raises:
        ValueError: If ``facecolors`` or ``c`` is passed — open versus solid
            is what this plotter decides; use ``color`` / ``edgecolors``.

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
            text = "|".join(node.elements[side].label for side in node.sides)
            target.annotate(
                text, point, textcoords="offset points", xytext=(3, 3), fontsize=7,
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
