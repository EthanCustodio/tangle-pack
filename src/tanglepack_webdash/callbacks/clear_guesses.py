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
)


def register(app: Dash):
    @app.callback(
        Output("fp-guess-points", "data", allow_duplicate=True),
        Input("btn-clear-guesses", "n_clicks"),
        prevent_initial_call=True,
    )
    def _clear_points(n):
        return []
