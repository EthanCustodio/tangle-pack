# src/tanglepack_webdash/app.py
from __future__ import annotations
from pathlib import Path
from dash import Dash, Input, Output, State
from dash.dependencies import ClientsideFunction
from .layout.main_layout import build_layout
from .callbacks import (
    fixed_point,
    session_init,
    build_system,
    manifolds,
    orientation,
    intersections,
    bridges,
    click_modes,
    clear_guesses,
    bridge_ops,
)


def make_app() -> Dash:
    # Ensure Dash serves /assets (so clientside.js loads)
    pkg_dir = Path(__file__).parent
    assets_dir = pkg_dir / "assets"
    app = Dash(__name__, assets_folder=str(assets_dir))
    app.title = "Tangle Workbench"
    app.layout = build_layout()

    # ---- register your existing server-side callbacks
    session_init.register(app)
    build_system.register(app)
    fixed_point.register(app)
    manifolds.register(app)
    orientation.register(app)
    intersections.register(app)
    bridges.register(app)
    click_modes.register(app)
    clear_guesses.register(app)
    bridge_ops.register(app)

    app.clientside_callback(
        ClientsideFunction(namespace="fp", function_name="set_mode"),
        Output("mode-sync", "children"),
        Input("click-mode", "value"),
        prevent_initial_call=False,
    )

    app.clientside_callback(
        ClientsideFunction(namespace="fp", function_name="pull_points"),
        Output("fp-guess-points", "data"),
        Input("fp-sync", "n_intervals"),
        prevent_initial_call=False,
    )

    # B) Render list of points as "(x, y)" with 3 decimals
    app.clientside_callback(
        ClientsideFunction(namespace="fp", function_name="points_to_label"),
        Output("click-points", "children"),
        Input("fp-guess-points", "data"),
        prevent_initial_call=False,
    )

    # C) Keep a visible "Clicked Points" marker layer in sync with the store
    app.clientside_callback(
        ClientsideFunction(namespace="fp", function_name="sync_points_trace"),
        Output("points-sync", "children"),
        Input("fp-guess-points", "data"),
        prevent_initial_call=False,
    )

    # D) Live cursor readout (unchanged)
    app.clientside_callback(
        ClientsideFunction(namespace="fp", function_name="move_to_label"),
        Output("cursor-readout", "children"),
        Input("plot-events", "n_events"),
        State("plot-events", "event"),
        State("plot", "figure"),
        prevent_initial_call=False,
    )

    app.clientside_callback(
        ClientsideFunction(namespace="bridge", function_name="pull_selection"),
        Output("selected-bridge-idx", "data"),
        Input("fp-sync", "n_intervals"),  # Poll every 250ms
        prevent_initial_call=False,
    )

    return app


def main():
    app = make_app()

    @app.callback(
        Output("debug", "children"),
        Input("plot-events", "n_events"),
        State("plot-events", "event"),
        State("click-mode", "value"),
        State("fp-guess-points", "data"),
        prevent_initial_call=False,
    )
    def _dbg(n, evt, mode, pts):
        # keep this compact; helpful while testing
        return f"n_events={n} | mode={mode} | last.type={getattr(evt,'type',None)} | pts={len(pts) if pts else 0}"

    app.run(debug=True)


app = make_app()
server = app.server
