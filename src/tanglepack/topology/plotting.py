"""
Drawing routines for the topological layer.

Every ``plot_*`` of the topological layer lives here so the marker, colour and
z-order conventions are defined once; :class:`~.Trellis.Trellis` and
:class:`~..loom.TangleSession.TangleSession` keep thin delegating wrappers.

Dev Notes — topological plotting.

The z-orders are a stack, not arbitrary numbers: black intersections (drawn by
the numerical layer) sit lowest, then the magenta strong-pip CANDIDATES, then
the green chosen pip on top of the set it was chosen from, then the orange
pseudoneighbours, then the holes and their labels. Changing one means checking
the whole ladder.

:func:`plot_stable_partition` is the odd one out: it takes the partition
RESULTS rather than a trellis, because a number-line figure is routinely built
from partitions gathered across several trellises (see
``scripts/henon_pseudoneighbors_nested.py``). The trellis method passes its own
``stable_partitions``.
"""

from __future__ import annotations

import logging
from typing import Iterable, TYPE_CHECKING, Union

import numpy as np
import matplotlib.pyplot as plt

from .TopologyResults import StablePartitionResult

if TYPE_CHECKING:
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
