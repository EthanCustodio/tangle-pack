from dash import html, dcc

system_section = html.Div(
    className="controls",
    children=[
        html.Label("Preset:"),
        dcc.Dropdown(
            id="map-preset",
            options=[
                {"label": "Hénon (binary horseshoe)", "value": "henon_binary"},
                {"label": "Custom (use text boxes)", "value": "custom"},
            ],
            value="henon_binary",
            clearable=False,
            style={"width": "260px"},
        ),
        dcc.Input(
            id="fx", type="text", value="y - 10 + x**2, -x", style={"width": "360px"}
        ),
        dcc.Input(
            id="finv",
            type="text",
            value="-y, x + 10 - y**2",
            style={"width": "360px"},
        ),
        html.Button("Build system", id="btn-build", n_clicks=0),
    ],
)
