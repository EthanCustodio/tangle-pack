# tanglepack_webdash/callbacks/fixed_point.py
from __future__ import annotations
import numpy as np
from dash import Dash, Input, Output, State as DashState, no_update
from tanglepack.TangleWorkbench import TangleWorkbench
from ..sessions import get_state
from ..parser import parse_map_text
from ..utils.figures import blank_figure, add_fp_trace


def register(app: Dash):
    @app.callback(
        Output("status", "children"),
        Input("btn-build", "n_clicks"),
        DashState("sid", "data"),  # 👈 grab the session id
        DashState("fx", "value"),
        DashState("fx_inv", "value"),
        prevent_initial_call=True,
    )
    def build_system(_, sid, fx_text, fx_inv_text):

        if not sid:
            return "❌ Missing session id."
        state = get_state(sid)

        try:
            f = parse_map_text(fx_text)
            f_inv = parse_map_text(fx_inv_text)
            # For FP-only MVP, we can pass f for both forward/inverse
            state.wb = TangleWorkbench(f, f_inv)
            state.fp = None
            return "✅ System built."
        except Exception as e:
            state.wb = None
            state.fp = None
            return f"❌ Build error: {e}"
