from dash import html, dcc

initializor = html.Div(
    className="controls",
    children=[
        dcc.Input(
            id="u_dir",
            type="text",
            value="-1, 0",
            style={"width": "120px"},
        ),
        dcc.Input(
            id="s_dir",
            type="text",
            value="0, 1",
            style={"width": "120px"},
        ),
        html.Button("Set directions", id="btn-orient", n_clicks=0),
        html.Button("Reset manifolds", id="btn-reset", n_clicks=0),
        html.Button("Init manifolds", id="btn-init", n_clicks=0),
    ],
)
