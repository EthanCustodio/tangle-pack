# src/tanglepack_webdash/callbacks/session_init.py
from __future__ import annotations
from uuid import uuid4
from dash import Dash, Output, Input, State as DashState, no_update


def register(app: Dash):
    @app.callback(
        Output("sid", "data"),
        Input("url", "pathname"),  # fires on page load / navigation
        DashState("sid", "data"),
        prevent_initial_call=False,
    )
    def ensure_sid(_path, current_sid):
        return current_sid or str(uuid4())
