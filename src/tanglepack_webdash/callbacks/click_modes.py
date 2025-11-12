# src/tanglepack_webdash/callbacks/click_modes.py
from __future__ import annotations
from typing import List
import plotly.graph_objects as go
from dash import Dash, Input, Output, State, no_update

FP_GUESS_NAME = "FP guess"


def _remove_trace_by_name(fig_dict, name: str):
    fig_dict["data"] = [tr for tr in fig_dict.get("data", []) if tr.get("name") != name]


def _upsert_fp_guess_trace(fig_dict, pts: List[List[float]]):
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    idx = next(
        (
            i
            for i, tr in enumerate(fig_dict.get("data", []))
            if tr.get("name") == FP_GUESS_NAME
        ),
        None,
    )
    trace = dict(
        type="scatter",
        mode="markers",
        x=xs,
        y=ys,
        marker=dict(symbol="x", size=12, line=dict(width=2)),
        name=FP_GUESS_NAME,
        hovertemplate="FP guess: (%{x:.4f}, %{y:.4f})<extra></extra>",
    )
    if idx is None:
        fig_dict.setdefault("data", []).append(trace)
    else:
        fig_dict["data"][idx].update(x=xs, y=ys)


def register(app: Dash):
    # A) Compose overlay markers from the Store; DO NOT rebuild base figure here.
    # @app.callback(
    #     Output("plot", "figure", allow_duplicate=True),
    #     Input("fp-guess-points", "data"),
    #     State("plot", "figure"),
    #     prevent_initial_call=True,
    # )
    # def overlay_fp_guesses(pts, fig):
    #     if fig is None:
    #         return no_update
    #     pts = list(pts or [])
    #     _remove_trace_by_name(fig, FP_GUESS_NAME)
    #     if pts:
    #         _upsert_fp_guess_trace(fig, pts)
    #     return fig

    # B) Crosshair cursor class toggled by click-mode
    @app.callback(
        Output("plot", "className"),
        Input("click-mode", "value"),
        prevent_initial_call=False,
    )
    def _cursor(mode):
        return "clicking" if mode == "fp_guess" else ""
