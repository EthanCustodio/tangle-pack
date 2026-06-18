# tanglepack_webdash/callbacks/fixed_point.py
from __future__ import annotations
import numpy as np
import traceback
from dash import Dash, Input, Output, State as DashState, no_update
from tanglepack import TangleWorkbench
from tanglepack_webdash.maps import PRESETS
from ..sessions import get_state
from ..utils.figures import (
    blank_figure,
    add_bridges_only,
    plot_tangle,
    plot_all_bridges,
    plot_intersections,
)


def register(app: Dash):

    @app.callback(
        Output("plot", "figure", allow_duplicate=True),
        Output("status", "children", allow_duplicate=True),
        Input("btn-create-bridges", "n_clicks"),
        DashState("sid", "data"),
        prevent_initial_call=True,
    )
    def on_make_bridges(n, sid):
        if not n:
            return no_update, no_update
        if not sid:
            return no_update, "❌ Missing session id."

        state = get_state(sid)
        if state.workbench is None or state.fp is None:
            return no_update, "⚠️ Build a system & fixed point first."

        try:
            # 1) Ensure intersections are current
            # state.workbench.compute_intersections(
            #     state.fp
            # )  # fills Tangle._intersecting_*  :contentReference[oaicite:3]{index=3}

            # 2) Build bridges from the current tangle view
            bridges = state.workbench.create_bridges(
                state.fp
            )  # returns list[Bridge]           :contentReference[oaicite:4]{index=4}

            # 3) Update bridges in session state
            # if state.bridges is None:
            #     state.bridges = bridges
            # else:
            #     merged = state.bridges.copy()
            #     merged.update(bridges)
            #     state.bridges = merged
            state.bridges = bridges

            # 4) Rebuild a "bridges-only" figure (FP + intersections + bridges)
            # state.fig = add_bridges_only(state, blank_figure())
            state.fig = plot_tangle(state, state.workbench, state.fp, "stable")
            state.fig = plot_all_bridges(state, bridges)
            state.fig = plot_intersections(state, state.workbench, state.fp)

            return state.fig, f"✅ Created {len(bridges)} bridge(s)."
        except Exception as e:
            return no_update, f"❌ Bridge error: {e}\n{traceback.format_exc()}"
