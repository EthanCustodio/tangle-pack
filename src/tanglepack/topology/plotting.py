"""
Drawing routines for the topological layer.

Every ``plot_*`` of the topological layer lives here so the marker, colour and
z-order conventions are defined once; :class:`~.Trellis.Trellis` and
:class:`~..loom.TangleSession.TangleSession` keep thin delegating wrappers.

Dev Notes — topological plotting.

The z-orders are a stack, not arbitrary numbers: black intersections (drawn by
the numerical layer) sit lowest, then the dual graph's own edges (z=1) and
face points (z=2) — kept low so a dual graph drawn over a trellis plot never
hides the numerics dots underneath it — then its arc nodes (z=10), then the
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
:func:`plot_transition_graph` is the same shape one level up: it takes a
:class:`~.SymbolicDynamics.SymbolicDynamics` directly and draws its
``networkx.DiGraph`` of symbols.
"""

from __future__ import annotations

import logging
from typing import Iterable, TYPE_CHECKING, Union

import numpy as np
import matplotlib.pyplot as plt
from numpy.typing import NDArray

from .TopologyResults import StablePartitionResult

if TYPE_CHECKING:
    from .DualGraph import ArcNode, DualGraph, FaceNode
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

#: Style defaults for :func:`plot_transition_graph`'s nodes, labels and edges.
#: ``node_color``/``node_size``/``edgecolors``/``font_size`` go to the node and
#: label draws; ``edge_color``/``arrows``/``arrowsize``/``connectionstyle`` go
#: to the edge draw.
TRANSITION_GRAPH_STYLE = {
    "node_color": "lightsteelblue",
    "node_size": 600,
    "edgecolors": "black",
    "font_size": 9,
    "edge_color": "gray",
    "arrows": True,
    "arrowsize": 12,
    "connectionstyle": "arc3,rad=0.1",
}

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
        **line_kwargs: Forwarded to the interval ``plot`` calls (e.g.
            ``linewidth``, ``color``).

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

    labels = []
    extent_lo, extent_hi = np.inf, -np.inf
    for row, result in enumerate(results):
        color = "tab:blue" if result.side == "left" else "tab:red"
        labels.append(
            f"{result.side} (p{result.branch_key[0].period}, "
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
    target.set_yticklabels(labels)
    target.set_ylim(-0.6, len(results) - 0.4)
    target.set_xlabel("stable canonical distance")
    target.set_title(
        "Stable manifold partition — [ ] closed, ( ) open, ○ singleton; "
        "numbers = ids"
    )
    return target


# ── the dual graph ──────────────────────────────────────────────────────────


def _arc_midpoints(
    trellis: "Trellis", nodes: Iterable["ArcNode"]
) -> list[NDArray[np.float64]]:
    """The usable (non-degenerate) midpoints of a run of arc nodes."""
    mids: list[NDArray[np.float64]] = []
    for node in nodes:
        mid = node.midpoint(trellis)
        if mid is None:
            logger.debug(
                "arc node %s has a degenerate midpoint; skipping it", node.key
            )
            continue
        mids.append(mid)
    return mids


def _push_outside_bbox(
    dual_graph: "DualGraph", point: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Push a point outside the bounding box of every arc node's midpoint.

    Lands on the bbox corner in ``point``'s own quadrant (relative to the bbox
    centre), then a further :data:`DUAL_GRAPH_PUSH_FRACTION` of the bbox
    diagonal outward.
    """
    all_mids = _arc_midpoints(dual_graph.trellis, dual_graph.arc_nodes.values())
    if not all_mids:
        return point
    points = np.vstack(all_mids)
    mins = points.min(axis=0)
    maxs = points.max(axis=0)
    center = 0.5 * (mins + maxs)
    diagonal = float(np.linalg.norm(maxs - mins))

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
        cached array); otherwise the mean of the midpoints of its boundary arc
        nodes. Only the ONE node with :attr:`~.DualGraph.FaceNode.is_unbounded`
        set has that mean pushed outside the bounding box of every arc node's
        midpoint in the graph (see :func:`_push_outside_bbox`) — a bounded
        face's point is returned as computed, with no guarantee it falls
        outside any other face (a region's representative point is itself only
        a heuristic; see its docstring for the degenerate case where it lands
        well outside a tight tangle).

    Note:
        A region's ``representative_point`` can itself be None (a boundary too
        degenerate to have an interior); that, like an "open"/"outer" node,
        falls back to the arc-midpoint mean.
    """
    if face_node.kind == "region" and face_node.faces:
        point = face_node.faces[0].representative_point
        if point is not None:
            return point.copy()
        logger.debug(
            "face node %d (region) has no representative point; falling back "
            "to the mean of its arc midpoints",
            face_node.index,
        )

    mids = _arc_midpoints(dual_graph.trellis, face_node.arc_nodes)
    if mids:
        mean = np.mean(np.vstack(mids), axis=0).astype(np.float64)
    else:
        mean = np.zeros(2, dtype=np.float64)
        logger.debug(
            "face node %d has no arc midpoints to average; defaulting to the "
            "origin before any outward push",
            face_node.index,
        )

    if not face_node.is_unbounded:
        return mean
    return _push_outside_bbox(dual_graph, mean)


def _clip_axes_to_arcs(ax, dual_graph: "DualGraph") -> None:
    """
    Limit an axes to the padded bbox of every arc midpoint and the unbounded
    face point.

    A bounded region's :attr:`~.TopologyResults.Region.representative_point`
    is only a heuristic and can, in a degenerate case, land far outside a
    tight tangle; left to matplotlib's autoscale that single outlier point
    blows the view out past anything useful, so :func:`plot_dual_graph`'s
    default view is clipped to what the arcs (and the one unbounded point)
    actually span instead.
    """
    points = _arc_midpoints(dual_graph.trellis, dual_graph.arc_nodes.values())
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
    Draw a dual graph: arc nodes, face points, and the edges between them.

    Arc nodes are drawn at their midpoint — a hollow circle when the arc is a
    wall, filled when it is passable (:attr:`~.DualGraph.ArcNode.filled`). Each
    face node gets a small dot at its :func:`face_point`, joined by a thin grey
    line to every arc node on its boundary.

    Args:
        dual_graph: The graph to draw.
        ax: Optional matplotlib Axes to draw on. Defaults to the current axes
            (plt).
        show_labels: If True, annotate each arc node with
            ``f"{left.label}|{right.label}"`` and each face node with its
            index and kind.
        clip_to_arcs: If True (default), set the axes limits to the bounding
            box of every arc midpoint plus the unbounded face node's point,
            padded by :data:`DUAL_GRAPH_PUSH_FRACTION` of that bbox's diagonal
            (see :func:`_clip_axes_to_arcs`) — a bounded region's
            representative point can otherwise sit far outside the tangle and
            blow out matplotlib's autoscale. Pass False to let the axes
            autoscale to everything drawn, outliers included.
        **scatter_kwargs: Forwarded to both the hollow and the filled arc-node
            scatter, overriding :data:`DUAL_GRAPH_ARC_STYLE`. ``facecolors``
            and ``c`` are managed internally (they are how a wall is told from
            a passable node) and raise ``ValueError`` if passed; use
            ``color`` (fill/edge colour) and ``edgecolors`` instead. Face
            points and the edges to them keep their own fixed style
            (:data:`DUAL_GRAPH_FACE_STYLE`, :data:`DUAL_GRAPH_EDGE_STYLE`).

    Returns:
        The Axes drawn on.

    Raises:
        ValueError: If ``scatter_kwargs`` tries to override ``facecolors`` or
            ``c``.
    """
    target = ax if ax is not None else plt.gca()

    reserved = _RESERVED_ARC_KWARGS.intersection(scatter_kwargs)
    if reserved:
        raise ValueError(
            f"plot_dual_graph manages {sorted(reserved)} itself (hollow vs. "
            "filled arc nodes); pass 'color'/'edgecolors' instead"
        )

    style = dict(DUAL_GRAPH_ARC_STYLE)
    style.update(scatter_kwargs)
    color = style.pop("color")
    edgecolors = style.pop("edgecolors", color)

    hollow: list[tuple] = []
    filled: list[tuple] = []
    for node in dual_graph.arc_nodes.values():
        mid = node.midpoint(dual_graph.trellis)
        if mid is None:
            logger.debug(
                "arc node %s has a degenerate midpoint; skipping it", node.key
            )
            continue
        (filled if node.filled else hollow).append((node, mid))

    # Always scattered, even with zero points, so the two arc-node collections
    # land at fixed positions (0, 1) in ``target.collections`` regardless of
    # whether either set is empty -- callers (tests included) can tell hollow
    # from filled without inspecting facecolors.
    hollow_coords = (
        np.vstack([mid for _node, mid in hollow]) if hollow else np.empty((0, 2))
    )
    filled_coords = (
        np.vstack([mid for _node, mid in filled]) if filled else np.empty((0, 2))
    )
    target.scatter(
        hollow_coords[:, 0], hollow_coords[:, 1],
        facecolors="none", edgecolors=edgecolors, **style,
    )
    target.scatter(filled_coords[:, 0], filled_coords[:, 1], color=color, **style)

    if show_labels:
        for node, mid in hollow + filled:
            target.annotate(
                f"{node.left.label}|{node.right.label}", mid,
                textcoords="offset points", xytext=(3, 3), fontsize=7,
            )

    face_points: dict[int, NDArray[np.float64]] = {}
    for face_node in dual_graph.face_nodes:
        point = face_point(dual_graph, face_node)
        face_points[face_node.index] = point
        for arc_node in face_node.arc_nodes:
            mid = arc_node.midpoint(dual_graph.trellis)
            if mid is None:
                continue
            target.plot(
                [point[0], mid[0]], [point[1], mid[1]], **DUAL_GRAPH_EDGE_STYLE,
            )
        if show_labels:
            target.annotate(
                f"{face_node.index}:{face_node.kind}", point,
                textcoords="offset points", xytext=(3, -3),
                fontsize=7, color=DUAL_GRAPH_FACE_STYLE["color"],
            )

    if face_points:
        coords = np.vstack(list(face_points.values()))
        target.scatter(coords[:, 0], coords[:, 1], **DUAL_GRAPH_FACE_STYLE)

    if clip_to_arcs:
        _clip_axes_to_arcs(target, dual_graph)

    return target


# ── the transition graph ────────────────────────────────────────────────────


def plot_transition_graph(
    symbolic_dynamics: "SymbolicDynamics",
    ax=None,
    **kwargs,
):
    """
    Draw a symbolic dynamics' transition graph: symbols as nodes, weighted edges.

    Nodes are laid out with ``networkx.spring_layout(seed=0)`` — deterministic
    across calls, unlike an unseeded spring layout — and labelled with the
    symbol itself. An edge ``u -> v`` is drawn whenever ``u``'s word contains
    ``v``, with an arrowhead; it is additionally labelled with its ``weight``
    (how many times ``v`` occurs in ``u``'s word) whenever that count exceeds 1,
    so a single occurrence keeps the picture uncluttered.

    Args:
        symbolic_dynamics: The
            :class:`~.SymbolicDynamics.SymbolicDynamics` whose
            :attr:`~.SymbolicDynamics.SymbolicDynamics.transition_graph` is
            drawn.
        ax: Optional matplotlib Axes to draw on. Defaults to the current axes
            (plt).
        **kwargs: Any name from :data:`TRANSITION_GRAPH_STYLE` overrides that
            style default. Anything else is forwarded straight through to
            ``nx.draw_networkx_nodes`` (e.g. ``alpha``, ``node_shape``,
            ``linewidths``) — so an unrecognized name raises the same
            ``TypeError`` calling that function directly would, rather than
            being silently dropped.

    Returns:
        The Axes drawn on.

    Raises:
        TypeError: If a keyword is neither a :data:`TRANSITION_GRAPH_STYLE`
            name nor a parameter of ``nx.draw_networkx_nodes``.
    """
    import networkx as nx

    target = ax if ax is not None else plt.gca()
    graph = symbolic_dynamics.transition_graph

    style = dict(TRANSITION_GRAPH_STYLE)
    extra_node_kwargs = {}
    for key, value in kwargs.items():
        if key in style:
            style[key] = value
        else:
            extra_node_kwargs[key] = value

    positions = nx.spring_layout(graph, seed=0)

    nx.draw_networkx_nodes(
        graph, positions, ax=target,
        node_color=style["node_color"], node_size=style["node_size"],
        edgecolors=style["edgecolors"], **extra_node_kwargs,
    )
    nx.draw_networkx_labels(graph, positions, ax=target, font_size=style["font_size"])
    nx.draw_networkx_edges(
        graph, positions, ax=target,
        edge_color=style["edge_color"], arrows=style["arrows"],
        arrowsize=style["arrowsize"], node_size=style["node_size"],
        connectionstyle=style["connectionstyle"],
    )

    edge_labels = {
        (u, v): data["weight"]
        for u, v, data in graph.edges(data=True)
        if data.get("weight", 1) > 1
    }
    if edge_labels:
        nx.draw_networkx_edge_labels(
            graph, positions, edge_labels=edge_labels, ax=target,
            font_size=style["font_size"],
        )

    target.set_axis_off()
    return target
