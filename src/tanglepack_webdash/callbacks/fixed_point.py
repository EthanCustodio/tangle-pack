# tanglepack_webdash/callbacks/fixed_point.py
from __future__ import annotations
import numpy as np
from dash import Dash, Input, Output, State, no_update
from tanglepack.TangleWorkbench import TangleWorkbench
from ..state import wb, fp  # module-level globals (simple local use)
from ..parser import parse_map_text
from ..utils.figures import blank_figure, add_fp_trace


def register(app: Dash):
    @app.callback(
        Output("status", "children"),
        Input("btn-build", "n_clicks"),
        State("fx", "value"),
        prevent_initial_call=True,
    )
    def build_system(_, fx_text):
        global wb, fp
        try:
            f = parse_map_text(fx_text)
            # For FP-only MVP, we can pass f for both forward/inverse
            wb = TangleWorkbench(f, f)
            fp = None
            return "✅ System built."
        except Exception as e:
            wb = None
            fp = None
            return f"❌ Build error: {e}"

    @app.callback(
        Output("plot", "figure", allow_duplicate=True),
        Output("status", "children", allow_duplicate=True),
        Input("btn-fp", "n_clicks"),
        State("x0y0", "value"),
        prevent_initial_call=True,
    )
    def find_fp(_, x0y0_text):
        global wb, fp
        if wb is None:
            return no_update, "ℹ️ Build the system first."
        try:
            x0, y0 = [float(s.strip()) for s in x0y0_text.split(",")]
        except Exception:
            return no_update, f"❌ Bad initial guess: {x0y0_text!r}"

        try:
            fp = wb.construct_fixed_point(np.array([x0, y0], float))
            coords = np.asarray(fp.coordinates)
            pt = (
                coords[0]
                if (coords.ndim == 2 and coords.shape[1] == 2)
                else coords.reshape(-1)[:2]
            )
        except Exception as e:
            return no_update, f"❌ Fixed-point error: {e}"

        fig = blank_figure()
        fig = add_fp_trace(fig, float(pt[0]), float(pt[1]))
        return fig, f"✅ FP: ({pt[0]:.6f}, {pt[1]:.6f})"
