# tanglepack_webdash/callbacks/fixed_point.py
from __future__ import annotations
import numpy as np
from dash import Dash, Input, Output, State as DashState, no_update
from tanglepack.TangleWorkbench import TangleWorkbench
from tanglepack_webdash.maps import PRESETS
from ..sessions import get_state
from ..parser import parse_map_text
from ..utils.wrappers import pointize
from ..utils.figures import blank_figure, add_fp_trace


def _check_inverse(f, finv, samples):
    """return max ||f(finv(p)) - p|| and ||finv(f(p)) - p|| over sample points"""
    errs = []
    for p in samples:
        p = np.asarray(p, float)
        e1 = np.linalg.norm(f(finv(p)) - p)
        e2 = np.linalg.norm(finv(f(p)) - p)
        errs.append(max(e1, e2))
    return float(np.max(errs))


def register(app: Dash):
    @app.callback(
        Output("status", "children", allow_duplicate=True),
        Input("btn-build", "n_clicks"),
        DashState("sid", "data"),
        DashState("map-preset", "value"),
        DashState("fx", "value"),
        DashState("finv", "value"),
        prevent_initial_call=True,
    )
    def build_system(_n, sid, preset, fx_text, finv_text):
        if not sid:
            return "❌ Missing session id."
        st = get_state(sid)

        try:
            if preset and preset != "custom":
                label, factory = PRESETS[preset]
                f, finv = factory()
                msg_src = f"preset: {label}"
            else:
                # only if you explicitly choose "custom"
                f = parse_map_text(fx_text)
                finv = parse_map_text(finv_text)
                msg_src = "custom text"
        except Exception as e:
            st.wb = None
            st.fp = None
            return f"❌ Build error: {e}"

        st.wb = TangleWorkbench(f, finv)
        st.fp = None
        # if you track orientation: st.did_orient = False
        return f"✅ System built from {msg_src}."
