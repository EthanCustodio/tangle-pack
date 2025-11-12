from dash import html, dcc

manifold_section = html.Div(
    className="controls",
    children=[
        html.Button(
            "Grow ×1 (unstable)",
            id="btn-grow-u",
            n_clicks=0,
        ),
        html.Button("Grow ×1 (stable)", id="btn-grow-s", n_clicks=0),
        html.Button(
            "Grow until Turnaround",
            id="btn-grow-turn",
            n_clicks=0,
            className="btn",
        ),
        html.Button(
            "Grow until Intersection",
            id="btn-grow-inter",
            n_clicks=0,
            className="btn",
        ),
        html.Button(
            "Trim Stable Manifold",
            id="btn-trim_manifold",
            n_clicks=0,
            className="btn",
        ),
    ],
)
