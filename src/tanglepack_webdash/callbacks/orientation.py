from __future__ import annotations
import numpy as np
from dash import Dash, Input, Output, State as DashState
from tanglepack_webdash.sessions import get_state


def _parse_vec2(txt: str) -> np.ndarray:
    """Parse 'a, b' → unit vector [a, b]. Allows extra spaces. Zeros stay zero."""
    a, b = [float(s.strip()) for s in txt.split(",")]
    v = np.array([a, b], dtype=float)
    n = float(np.linalg.norm(v))
    return v if n == 0.0 else (v / n)


def register(app: Dash):
    @app.callback(
        Output("status", "children", allow_duplicate=True),
        Input("btn-orient", "n_clicks"),
        DashState("sid", "data"),
        DashState("u_dir", "value"),
        DashState("s_dir", "value"),
        prevent_initial_call=True,
    )
    def orient(_n, sid, u_txt, s_txt):
        if not sid:
            return "❌ Missing session id."
        state = get_state(sid)
        if state.workbench is None or state.fp is None:
            return "ℹ️ Build the system and find a fixed point first."

        try:
            u_dir = _parse_vec2(u_txt)
            s_dir = _parse_vec2(s_txt)
            state.workbench.orient_eigenvectors(
                state.fp, {"unstable": u_dir, "stable": s_dir}
            )
            return f"✅ Directions set. Unstable≈({u_dir[0]:.3f}, {u_dir[1]:.3f}), Stable≈({s_dir[0]:.3f}, {s_dir[1]:.3f})"
        except Exception as e:
            return f"❌ Orient error: {e}"
