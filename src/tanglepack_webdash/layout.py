# tanglepack_webdash/layout.py
from __future__ import annotations
from dash import html, dcc
from .utils.figures import blank_figure


def build_layout():
    return html.Div(
        [
            dcc.Location(id="url"),  # allows initial-load callback
            dcc.Store(id="sid", storage_type="session"),  # holds session id per tab
            html.H3(
                "Tangle Workbench — Fixed Point",
                style={"color": "white", "margin": "10px"},
            ),
            html.Div(
                [
                    html.Label("Map and Inverse Map", style={"marginRight": "8px"}),
                    dcc.Input(
                        id="fx",
                        type="text",
                        value="y - 10 + x**2, -x",
                        style={"width": "420px"},
                    ),
                    dcc.Input(
                        id="fx_inv",
                        type="text",
                        value="-y, x + 10 - y**2",
                        style={"width": "420px"},
                    ),
                    html.Button(
                        "Build system",
                        id="btn-build",
                        n_clicks=0,
                        style={"marginLeft": "10px"},
                    ),
                ],
                className="controls",
            ),
            html.Div(
                [
                    html.Label("Initial guess (x0, y0)", style={"marginRight": "8px"}),
                    dcc.Input(
                        id="x0y0", type="text", value="4, -4", style={"width": "120px"}
                    ),
                    html.Button(
                        "Find fixed point",
                        id="btn-fp",
                        n_clicks=0,
                        style={"marginLeft": "10px"},
                    ),
                ],
                className="controls",
            ),
            dcc.Graph(id="plot", figure=blank_figure(), style={"height": "70vh"}),
            html.Div(id="status", style={"color": "#aaa", "padding": "8px"}),
        ],
        style={
            "background": "#000",
            "height": "100vh",
            "margin": 0,
            "padding": 0,
            "fontFamily": "system-ui",
        },
    )
