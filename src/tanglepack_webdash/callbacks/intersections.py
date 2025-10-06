# tanglepack_webdash/callbacks/fixed_point.py
from __future__ import annotations
import numpy as np
import traceback
from dash import Dash, Input, Output, State as DashState, no_update
from tanglepack.TangleWorkbench import TangleWorkbench
from tanglepack_webdash.maps import PRESETS
from ..sessions import get_state
from ..utils.figures import (
    blank_figure,
    add_fp_trace,
    add_intersection_traces,
    add_manifold_traces,
)


def register(app: Dash):
    @app.callback(
        Output("status", "children", allow_duplicate=True),
        Input("btn-cmpt_intersections", "n_clicks"),
        DashState("sid", "data"),
        prevent_initial_call=True,
    )
    def compute_intersections(_n, sid):

        if not sid:
            return "❌ Missing session id."
        state = get_state(sid)

        if state.wb is None or state.fp is None:
            return no_update, "⚠️ Need a system + fixed point first.", ""

        try:
            intersections = state.wb.compute_intersections(state.fp)
            # fig = add_manifold_traces(state.wb, state.fp)
            # fig = add_intersection_traces(state.wb, state.fp, fig)
            # fig = add_manifold_traces(state.fp)
            # return fig, "✅ Grown until intersection.", ""
            return (
                (
                    f"✅ Grown until intersection. "
                    f"Number of intersections: {len(intersections)}",
                ),
                f"Intersections: {intersections}",
            )
        except Exception as e:
            return no_update, f"❌ Grow error: {e}", traceback.format_exc()
