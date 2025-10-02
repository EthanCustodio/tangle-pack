from __future__ import annotations
import numpy as np
import traceback
from dash import Dash, Input, Output, State as DashState, no_update, ctx
from ..sessions import get_state
from tanglepack_webdash.utils.figures import (
    blank_figure,
    add_fp_trace,
    add_manifold_line,
    add_manifold_traces,
)


def _figure_from_state(st) -> "go.Figure":
    """Build a fresh figure with FP + all manifolds for the current FP."""
    fig = blank_figure()

    # Fixed point
    if st.fp is not None:
        coords = np.asarray(st.fp.coordinates)
        pt = (
            coords[0]
            if (coords.ndim == 2 and coords.shape[1] == 2)
            else coords.reshape(-1)[:2]
        )
        add_fp_trace(fig, float(pt[0]), float(pt[1]))

    # Manifolds for this fixed point
    if st.wb is not None and st.fp is not None:
        # Your TangleWorkbench.manifolds: { (fp, stability, oi, bi) : manifold }
        for (kfp, stab, oi, bi), M in st.wb.manifolds.items():
            if kfp is not st.fp:
                continue
            try:
                arr = M.get_point_array()  # (N,2)
            except Exception:
                # If a branch warning pops, skip this one
                continue
            color = "#3b82f6" if stab == "unstable" else "#ef4444"  # blue/red
            name = f"{stab}[o{oi},b{bi}]"
            add_manifold_line(fig, arr, color=color, name=name)

    return fig


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
        st = get_state(sid)
        if st.wb is None or st.fp is None:
            return no_update, "ℹ️ Build + FP first."
        # Purge manifolds *for this FP only* so other FPs (if any) are untouched
        st.wb.manifolds = {
            k: v for k, v in st.wb.manifolds.items() if k[0] is not st.fp
        }
        fig = _figure_from_state(st)
        return fig, "✅ Manifolds reset."

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
        st = get_state(sid)
        if st.wb is None or st.fp is None:
            return no_update, "ℹ️ Build the system and find a fixed point first."

        try:
            # This should stamp branch indices internally; mirrors your script
            st.wb.initialize_both_manifolds(st.fp)
        except Exception as e:
            return no_update, f"❌ Init error: {e}"

        fig = _figure_from_state(st)
        return fig, "✅ Manifolds initialized."

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
        st = get_state(sid)
        if st.wb is None or st.fp is None:
            return no_update, "ℹ️ Build + FP + Init manifolds first."

        # Which button fired?
        which = ctx.triggered_id
        if which not in ("btn-grow-u", "btn-grow-s"):
            return no_update, no_update

        try:
            stab = "unstable" if which == "btn-grow-u" else "stable"
            st.wb.grow_n_times(st.fp, stab, num_iterations=1)
        except Exception as e:
            return no_update, f"❌ Grow error: {e}"

        fig = _figure_from_state(st)
        return fig, f"✅ Grew {stab} ×1."

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

        if state.wb is None or state.fp is None:
            return no_update, "⚠️ Need a system + fixed point first.", ""

        try:
            state.wb.grow_until_turnaround(state.fp, "stable")
            state.wb.grow_until_turnaround(state.fp, "unstable")
            # fig = add_manifold_traces(state.fp)
            fig = _figure_from_state(state)
            return fig, "✅ Grown until turnaround.", ""
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

        if state.wb is None or state.fp is None:
            return no_update, "⚠️ Need a system + fixed point first.", ""

        try:
            # state.wb.grown_until_intersection(state.fp, "stable")
            state.wb.grown_until_intersection(state.fp, "unstable")
            fig = add_manifold_traces(state.fp)
            return fig, "✅ Grown until intersection.", ""
        except Exception as e:
            return no_update, f"❌ Grow error: {e}", traceback.format_exc()
