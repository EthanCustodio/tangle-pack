# tanglepack_webdash/utils/figures.py
from __future__ import annotations
import plotly.graph_objects as go
import numpy as np


def set_square_aspect(fig, enable: bool = True):
    if enable:
        fig.update_layout(yaxis=dict(scaleanchor="x", scaleratio=1))
    else:
        # remove anchor to go back to free aspect
        fig.update_layout(yaxis=dict(scaleanchor=None, scaleratio=None))
    return fig


def blank_figure():
    figure = go.Figure(
        layout=dict(
            paper_bgcolor="black",
            plot_bgcolor="black",
            xaxis=dict(color="white"),
            yaxis=dict(color="white"),
            showlegend=False,
        )
    )

    figure.update_layout(
        paper_bgcolor="#0e0f13",  # page background (softer black)
        plot_bgcolor="#1e1e22",  # plotting area (very dark gray)
    )
    set_square_aspect(figure, True)
    return figure


def add_fp_trace(fig: go.Figure, x: float, y: float) -> go.Figure:
    fig.add_trace(
        go.Scattergl(
            x=[x],
            y=[y],
            mode="markers",
            marker=dict(size=9, color="white", line=dict(width=1, color="black")),
            name="FP",
            hovertemplate="FP: (%{x:.6f}, %{y:.6f})<extra></extra>",
        )
    )
    set_square_aspect(fig, True)
    return fig


def add_manifold_line(fig: go.Figure, arr: np.ndarray, *, color: str, name: str):
    """arr is (N,2) coordinates."""
    if arr.size == 0:
        return
    fig.add_trace(
        go.Scattergl(
            x=arr[:, 0],
            y=arr[:, 1],
            mode="lines",
            line=dict(width=2, color=color),
            name=name,
            hoverinfo="skip",
        )
    )
