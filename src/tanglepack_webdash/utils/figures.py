# tanglepack_webdash/utils/figures.py
"""
Streamlined figure construction utilities for TanglePack web dashboard.
Follows backend patterns from BaseManifold and TangleWorkbench.
"""
from __future__ import annotations
import plotly.graph_objects as go
import numpy as np
from typing import Optional, List

# ============================================================================
# CONSTANTS
# ============================================================================

CLICK_SURFACE_NAME = "_click_surface"
DEFAULT_RANGE = [-15, 15]

# Color palettes
UNSTABLE_COLOR = "#3b82f6"  # blue
STABLE_COLOR = "#ef4444"  # red
FP_COLOR = "white"
INTERSECTION_COLOR = "white"

BRIDGE_PALETTE = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
    "#393b79",
    "#637939",
    "#8c6d31",
    "#843c39",
    "#7b4173",
    "#3182bd",
    "#e6550d",
    "#31a354",
    "#756bb1",
    "#636363",
]


# ============================================================================
# FIGURE INITIALIZATION
# ============================================================================


def blank_figure(xrange: tuple = None, yrange: tuple = None) -> go.Figure:
    """
    Create a blank figure with standard settings for dynamical systems plots.

    Args:
        xrange: Tuple (xmin, xmax) for initial x-axis range
        yrange: Tuple (ymin, ymax) for initial y-axis range

    Returns:
        Configured Plotly Figure with square aspect ratio and click surface
    """
    xrange = xrange or DEFAULT_RANGE
    yrange = yrange or DEFAULT_RANGE

    fig = go.Figure()
    fig.update_layout(
        clickmode="event+select",
        dragmode="pan",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=40, r=10, t=10, b=40),
        xaxis=dict(
            showgrid=True,
            gridwidth=0.5,
            gridcolor="rgba(255,255,255,0.06)",
            zeroline=False,
            range=list(xrange),
        ),
        yaxis=dict(
            showgrid=True,
            gridwidth=0.5,
            gridcolor="rgba(255,255,255,0.06)",
            zeroline=False,
            scaleanchor="x",
            scaleratio=1,
            range=list(yrange),
        ),
        uirevision="keep",  # preserve pan/zoom on redraws
        showlegend=False,
    )

    # Add invisible full-canvas trace to capture all clicks
    # _add_click_surface(fig, xrange, yrange)

    return fig


# ============================================================================
# CORE PLOTTING FUNCTIONS (mirrors backend patterns)
# ============================================================================


def plot_fixed_point(state, x: float, y: float, **kwargs) -> go.Figure:
    """
    Add fixed point marker to figure.

    Args:
        fig: Plotly figure
        x: x-coordinate
        y: y-coordinate
        **kwargs: Additional marker styling (size, color, etc.)

    Returns:
        Modified figure
    """

    marker_style = dict(
        size=kwargs.pop("size", 9),
        color=kwargs.pop("color", FP_COLOR),
        line=dict(width=1, color="black"),
    )
    marker_style.update(kwargs)

    state.fig.add_trace(
        go.Scattergl(
            x=[x],
            y=[y],
            mode="markers",
            marker=marker_style,
            name="FP",
            hovertemplate="FP: (%{x:.6f}, %{y:.6f})<extra></extra>",
        )
    )
    return state.fig


def plot_manifold(
    state, manifold, color: str = None, name: str = None, **kwargs
) -> go.Figure:
    """
    Add manifold curve to figure (mirrors BaseManifold.plot() pattern).

    Args:
        fig: Plotly figure
        manifold: BaseManifold object or (N,2) array of coordinates
        color: Line color (defaults based on stability if manifold object)
        name: Trace name
        **kwargs: Additional line styling

    Returns:
        Modified figure
    """
    # Handle both manifold objects and raw arrays
    if hasattr(manifold, "get_point_array"):
        arr = manifold.get_point_array()
        if color is None:
            color = UNSTABLE_COLOR if manifold.stability == "unstable" else STABLE_COLOR
        if name is None:
            name = f"{manifold.stability.capitalize()} manifold"
    else:
        arr = np.asarray(manifold)
        color = color or UNSTABLE_COLOR
        name = name or "Manifold"

    if arr.size == 0 or arr.shape[0] < 2:
        return state.fig

    line_style = dict(width=kwargs.pop("linewidth", 2), color=color)
    line_style.update(kwargs)

    state.fig.add_trace(
        go.Scatter(
            x=arr[:, 0],
            y=arr[:, 1],
            mode="lines",
            line=line_style,
            name=name,
            hoverinfo="none",
        )
    )
    return state.fig


def plot_intersections(
    state, workbench, fp, intersections: List = None, **kwargs
) -> go.Figure:
    """
    Add intersection points to figure (mirrors TangleWorkbench.plot_intersections()).

    Args:
        fig: Plotly figure
        workbench: TangleWorkbench object
        fp: FixedPoint object
        intersections: Optional pre-computed list of intersection coordinates
        **kwargs: Additional marker styling

    Returns:
        Modified figure
    """
    # Get intersection coordinates
    if intersections is None:
        arr = np.array(list(workbench.Tangle._intersecting_coords.values()))
        if arr.size == 0:
            return state.fig
    else:
        arr = np.array(list(intersections)) if intersections else np.array([])
        if arr.size == 0:
            return state.fig

    marker_style = dict(
        color=kwargs.pop("color", INTERSECTION_COLOR),
        size=kwargs.pop("size", 6),
        symbol=kwargs.pop("symbol", "circle"),
    )
    marker_style.update(kwargs)

    state.fig.add_trace(
        go.Scatter(
            x=arr[:, 0],
            y=arr[:, 1],
            mode="markers",
            name="Intersections",
            marker=marker_style,
            hovertemplate="x=%{x:.3f}<br>y=%{y:.3f}<extra>Intersection</extra>",
        )
    )
    return state.fig


def plot_bridges(
    state, bridges: List, palette: List[str] = None, **kwargs
) -> go.Figure:
    """
    Add bridge curves to figure (mirrors TangleWorkbench.plot_all_bridges()).

    Args:
        fig: Plotly figure
        bridges: List of Bridge objects
        palette: Color palette for bridges
        **kwargs: Additional line styling

    Returns:
        Modified figure
    """
    palette = palette or BRIDGE_PALETTE

    for i, bridge in enumerate(bridges):
        try:
            arr = bridge.get_point_array()
        except Exception:
            continue

        if arr.size == 0 or arr.shape[0] < 2:
            continue

        color = palette[i % len(palette)]
        state.fig.add_trace(
            go.Scatter(
                x=arr[:, 0],
                y=arr[:, 1],
                mode="lines",
                line=dict(width=kwargs.get("linewidth", 2), color=color),
                name=f"Bridge {i+1}",
                hoverinfo="skip",
            )
        )
    return state.fig


# ============================================================================
# TANGLE WORKBENCH PLOTTING METHODS (mirrors backend patterns)
# ============================================================================


def plot_tangle(state, workbench, fixed_point, stability: str, **kwargs) -> go.Figure:
    """
    Plot all manifolds of a given stability for a fixed point.
    Mirrors TangleWorkbench.plot_tangle() pattern.

    Args:
        fig: Plotly figure
        workbench: TangleWorkbench object
        fixed_point: FixedPoint object
        stability: "stable" or "unstable"
        **kwargs: Additional styling for manifold lines

    Returns:
        Modified figure
    """
    # Reset figure
    state.fig = blank_figure()

    # Plot all manifolds of this stability
    for (fp, stab, oi, bi), manifold in workbench.manifolds.items():
        if fp is fixed_point and stab == stability:
            try:
                state.fig = plot_manifold(
                    state, manifold, name=f"{stab}[o{oi},b{bi}]", **kwargs
                )
            except Exception:
                continue

    # Add fixed point markers
    for period in range(fixed_point.period):
        coords = fixed_point.coordinates[period]
        pt = coords if coords.ndim == 1 else coords.flatten()[:2]
        state.fig = plot_fixed_point(
            state,
            float(pt[0]),
            float(pt[1]),
        )

    return state.fig


def plot_all_bridges(
    state, bridges: List, colormap: str = "tab20", **kwargs
) -> go.Figure:
    """
    Plot all bridges with colormap.
    Mirrors TangleWorkbench.plot_all_bridges() pattern.

    Args:
        fig: Plotly figure
        bridges: List of Bridge objects
        colormap: Matplotlib colormap name (approximated)
        **kwargs: Additional line styling

    Returns:
        Modified figure
    """
    # Use our bridge palette (approximates tab20)
    return plot_bridges(state, bridges, palette=BRIDGE_PALETTE, **kwargs)


def plot_bridges_with_highlight(
    state, bridges: List, highlight_idx: int = None, palette: List[str] = None, **kwargs
) -> go.Figure:
    """
    Add bridge curves to figure with one highlighted.

    Args:
        state: WBState object
        bridges: List of Bridge objects
        highlight_idx: Index of bridge to highlight (None for no highlight)
        palette: Color palette for bridges
        **kwargs: Additional line styling

    Returns:
        Modified figure
    """
    palette = palette or BRIDGE_PALETTE

    for i, bridge in enumerate(bridges):
        try:
            arr = bridge.get_point_array()
        except Exception:
            continue

        if arr.size == 0 or arr.shape[0] < 2:
            continue

        if i == highlight_idx:
            # Highlighted bridge: thick gold line
            state.fig.add_trace(
                go.Scatter(
                    x=arr[:, 0],
                    y=arr[:, 1],
                    mode="lines",
                    line=dict(width=2, color="#FFD700"),
                    name=f"Bridge {i+1} (SELECTED)",
                    hoverinfo="skip",
                )
            )
        else:
            # Normal bridge
            color = palette[i % len(palette)]
            state.fig.add_trace(
                go.Scatter(
                    x=arr[:, 0],
                    y=arr[:, 1],
                    mode="lines",
                    line=dict(width=kwargs.get("linewidth", 2), color=color),
                    name=f"Bridge {i+1}",
                    hoverinfo="skip",
                )
            )
    return state.fig


# ============================================================================
# CONVENIENCE FUNCTIONS FOR BUILDING COMPLETE FIGURES
# ============================================================================


def build_figure_from_state(
    state, show_fp: bool = True, show_manifolds: bool = True
) -> go.Figure:
    """
    Build complete figure from session state (mirrors TangleWorkbench.plot_tangle()).

    Args:
        state: WBState object containing workbench and fp
        show_fp: Whether to include fixed point marker
        show_manifolds: Whether to include manifold curves

    Returns:
        Complete figure with all requested elements
    """
    state.fig = blank_figure()

    if state.fp is None:
        return state.fig

    # Fixed point
    if show_fp:
        coords = np.asarray(state.fp.coordinates)
        pt = (
            coords[0]
            if (coords.ndim == 2 and coords.shape[1] == 2)
            else coords.reshape(-1)[:2]
        )
        state.fig = plot_fixed_point(state, float(pt[0]), float(pt[1]))

    # Manifolds
    if show_manifolds and state.workbench is not None:
        for (kfp, stab, oi, bi), M in state.workbench.manifolds.items():
            if kfp is not state.fp:
                continue
            try:
                state.fig = plot_manifold(state, M, name=f"{stab}[o{oi},b{bi}]")
            except Exception:
                continue

    return state.fig


def build_bridge_figure(state) -> go.Figure:
    """
    Build figure showing bridges, intersections, and fixed point only.

    Args:
        state: WBState object with bridges attribute

    Returns:
        Figure with bridges, intersections, and FP
    """
    state.fig = blank_figure()

    if state.fp is None:
        return state.fig

    # Fixed point
    coords = np.asarray(state.fp.coordinates)
    pt = (
        coords[0]
        if (coords.ndim == 2 and coords.shape[1] == 2)
        else coords.reshape(-1)[:2]
    )
    plot_fixed_point(state.fig, float(pt[0]), float(pt[1]))

    # Intersections
    if state.workbench is not None:
        intersections = state.workbench.compute_intersections(state.fp)
        plot_intersections(state.fig, state.workbench, state.fp, intersections)

    # Bridges
    bridges = getattr(state, "bridges", None) or []
    if bridges:
        plot_bridges(state.fig, bridges)

    return state.fig


# ============================================================================
# LEGACY COMPATIBILITY (will be deprecated)
# ============================================================================


def add_fp_trace(state, x: float, y: float) -> go.Figure:
    """Legacy name for plot_fixed_point(). Use plot_fixed_point() instead."""
    return plot_fixed_point(state, x, y)


def add_manifold_line(state, arr: np.ndarray, *, color: str, name: str):
    """Legacy name for plot_manifold(). Use plot_manifold() instead."""
    return plot_manifold(state, arr, color=color, name=name)


def add_intersection_traces(state, workbench, fp, intersections, fig=None):
    """Legacy name. Use plot_intersections() instead."""
    state.fig = fig or go.Figure()
    return plot_intersections(state, workbench, fp, intersections)


def add_bridges_only(state, fig=None):
    """Legacy name. Use build_bridge_figure() instead."""
    return (
        build_bridge_figure(state)
        if state.fig is None
        else plot_bridges(state, getattr(state, "bridges", []))
    )


def set_square_aspect(state, enable: bool = True):
    """Set or remove square aspect ratio constraint."""
    if enable:
        state.fig.update_layout(yaxis=dict(scaleanchor="x", scaleratio=1))
    else:
        state.fig.update_layout(yaxis=dict(scaleanchor=None, scaleratio=None))
    return state.fig
