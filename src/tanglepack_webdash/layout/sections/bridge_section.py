from dash import html, dcc

bridge_section = html.Div(
    className="controls",
    children=[
        html.Button(
            "Compute Intersections",
            id="btn-cmpt_intersections",
            n_clicks=0,
            className="btn",
        ),
        html.Button(
            "Create Bridges",
            id="btn-create-bridges",
            n_clicks=0,
            className="btn",
        ),
    ],
)
