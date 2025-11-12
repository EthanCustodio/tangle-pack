# tanglepack_webdash/callbacks/fixed_point.py
from __future__ import annotations
import traceback
import numpy as np
from dash import Dash, Input, Output, State as DashState, no_update
from tanglepack.TangleWorkbench import TangleWorkbench
from ..sessions import get_state
from ..utils.wrappers import pointize
from ..parser import parse_map_text
from ..utils.figures import blank_figure, add_fp_trace


def register(app: Dash):
    @app.callback(
        Output("plot", "figure", allow_duplicate=True),
        Output("status", "children", allow_duplicate=True),
        Output("debug", "children", allow_duplicate=True),
        Input("btn-fp", "n_clicks"),
        DashState("click-mode", "value"),  # add
        DashState("fp-guess-points", "data"),
        DashState("sid", "data"),
        DashState("x0y0", "value"),
        DashState("map-preset", "value"),
        prevent_initial_call=True,
    )
    def find_fp(_, click_mode, guess_points, sid, text_input, preset):
        """
        Finds the fixed point when the button "btn-fp" is clicked.
        Uses the text inputed guess if there are no points selected.
        Uses 'clicked' points for initial guess otherwise.

        Args:
            _ (_type_): _description_
            click_mode (_type_): _description_
            guess_points (_type_): _description_
            sid (_type_): Session id for the current instance of the program.
            text_input (_type_): Text input for the fixed point guess.
            preset (_type_): _description_

        Returns:
            _type_: _description_

        Todo:
            Include input functionality for higher period fixed points.
            Include a checkbox to toggle the text input?
            Cleanup design in general.
        """

        if not sid:
            return no_update, "❌ Missing session id."
        state = get_state(sid)

        if state.workbench is None:
            return no_update, "ℹ️ Build the system first."

        if click_mode != "fp_guess":
            try:
                x0, y0 = [float(s.strip()) for s in text_input.split(",")]
            except Exception:
                return no_update, f"❌ Bad initial guess: {text_input!r}"

        if click_mode == "fp_guess" and guess_points:
            x0, y0 = guess_points[-1]
        try:
            guess = [x0, y0]
            state.fp = state.workbench.construct_fixed_point(guess)
            # state.fp = state.workbench.construct_fixed_point(np.array([x0, y0], float))
            coords = np.asarray(state.fp.coordinates)
            pt = (
                coords[0]
                if (coords.ndim == 2 and coords.shape[1] == 2)
                else coords.reshape(-1)[:2]
            )
        except Exception as e:
            # return no_update, f"❌ Fixed-point error: {e}"
            return no_update, f"❌ Fixed-point error: {e}", traceback.format_exc()

        state.fig = blank_figure()
        state.fig = add_fp_trace(state, float(pt[0]), float(pt[1]))
        return state.fig, f"✅ FP: ({pt[0]:.6f}, {pt[1]:.6f})", ""
