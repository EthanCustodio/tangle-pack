from __future__ import annotations
import numpy as np
import traceback
from dash import Dash, Input, Output, State as DashState, no_update, ctx
from ..sessions import get_state
from tanglepack_webdash.utils.figures import (
    blank_figure,
    add_fp_trace,
    add_manifold_line,
    build_figure_from_state,
    plot_tangle,
    plot_all_bridges,
    plot_intersections,
)


def register(app: Dash):
    @app.callback(
        Output("plot", "figure", allow_duplicate=True),
        Output("status", "children", allow_duplicate=True),
        Input("btn-reset", "n_clicks"),
        DashState("sid", "data"),
        prevent_initial_call=True,
    )
    def reset_manifolds(_n, sid):
        if not sid:
            return no_update, "❌ Missing session id."
        state = get_state(sid)
        if state.workbench is None or state.fp is None:
            return no_update, "ℹ️ Build + FP first."
        # Purge manifolds *for this FP only* so other FPs (if any) are untouched
        state.workbench.manifolds = {
            k: v for k, v in state.workbench.manifolds.items() if k[0] is not state.fp
        }
        state.fig = build_figure_from_state(state)
        return state.fig, "✅ Manifolds reset."

    # --- initialize both manifolds ---
    @app.callback(
        Output("plot", "figure", allow_duplicate=True),
        Output("status", "children", allow_duplicate=True),
        Input("btn-init", "n_clicks"),
        DashState("sid", "data"),
        prevent_initial_call=True,
    )
    def init_manifolds(_n, sid):
        if not sid:
            return no_update, "❌ Missing session id."
        state = get_state(sid)
        if state.workbench is None or state.fp is None:
            return no_update, "ℹ️ Build the system and find a fixed point first."

        try:
            # This should stamp branch indices internally; mirrors your script
            state.workbench.initialize_both_manifolds(state.fp)
        except Exception as e:
            return no_update, f"❌ Init error: {e}"

        state.fig = build_figure_from_state(state)
        return state.fig, "✅ Manifolds initialized."

    # --- grow one step (unstable or stable) ---
    @app.callback(
        Output("plot", "figure", allow_duplicate=True),
        Output("status", "children", allow_duplicate=True),
        Input("btn-grow-u", "n_clicks"),
        Input("btn-grow-s", "n_clicks"),
        DashState("sid", "data"),
        prevent_initial_call=True,
    )
    def grow_once(nu, ns, sid):
        if not sid:
            return no_update, "❌ Missing session id."
        state = get_state(sid)
        if state.workbench is None or state.fp is None:
            return no_update, "ℹ️ Build + FP + Init manifolds first."

        # Which button fired?
        which = ctx.triggered_id
        if which not in ("btn-grow-u", "btn-grow-s"):
            return no_update, no_update

        try:
            stab = "unstable" if which == "btn-grow-u" else "stable"
            state.workbench.grow_n_times(state.fp, stab, num_iterations=1)
        except Exception as e:
            return no_update, f"❌ Grow error: {e}"

        state.fig = build_figure_from_state(state)
        return state.fig, f"✅ Grew {stab} ×1."

    @app.callback(
        Output("plot", "figure", allow_duplicate=True),
        Output("status", "children", allow_duplicate=True),
        Output("debug", "children", allow_duplicate=True),
        Input("btn-grow-turn", "n_clicks"),
        DashState("sid", "data"),
        prevent_initial_call=True,
    )
    def grow_until_turnaround(_, sid):
        if not sid:
            return no_update, "❌ Missing session id.", ""
        state = get_state(sid)

        if state.workbench is None or state.fp is None:
            return no_update, "⚠️ Need a system + fixed point first.", ""

        try:
            state.workbench.grow_until_turnaround(state.fp, "stable")
            state.workbench.grow_until_turnaround(state.fp, "unstable")

            state.fig = build_figure_from_state(state)
            return state.fig, "✅ Grown until turnaround.", ""
        except Exception as e:
            return no_update, f"❌ Grow error: {e}", traceback.format_exc()

    @app.callback(
        Output("plot", "figure", allow_duplicate=True),
        Output("status", "children", allow_duplicate=True),
        Output("debug", "children", allow_duplicate=True),
        Input("btn-grow-inter", "n_clicks"),
        DashState("sid", "data"),
        prevent_initial_call=True,
    )
    def grow_until_intersection(_, sid):
        if not sid:
            return no_update, "❌ Missing session id.", ""
        state = get_state(sid)

        if state.workbench is None or state.fp is None:
            return no_update, "⚠️ Need a system + fixed point first.", ""

        try:
            state.workbench.grow_until_turnaround(state.fp, "stable")
            state.workbench.grown_until_intersection(state.fp, "unstable")

            state.fig = build_figure_from_state(state)
            return state.fig, "✅ Grown until intersection.", ""
        except Exception as e:
            return no_update, f"❌ Grow error: {e}", traceback.format_exc()

    @app.callback(
        Output("plot", "figure", allow_duplicate=True),
        Output("status", "children", allow_duplicate=True),
        Input("btn-trim_manifold", "n_clicks"),
        DashState("sid", "data"),
        prevent_initial_call=True,
    )
    def trim_stable_manifold(n, sid):
        if not sid:
            return no_update, "❌ Missing session id."
        state = get_state(sid)
        if state.workbench is None or state.fp is None:
            return no_update, "ℹ️ Build + FP + Init manifolds first."

        state.workbench.trim_stable_manifolds(state.fp)

        if state.bridges is None:
            state.fig = build_figure_from_state(state)
        else:
            state.fig = plot_tangle(state, state.workbench, state.fp, "stable")
            state.fig = plot_all_bridges(state, state.bridges)
            state.fig = plot_intersections(state, state.workbench, state.fp)

        # state.fig = build_figure_from_state(state)
        return state.fig, f"✅ Trimmed Stable Manifold."
