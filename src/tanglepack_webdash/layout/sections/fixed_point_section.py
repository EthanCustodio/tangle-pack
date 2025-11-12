from dash import html, dcc

fixed_point_section = html.Div(
    className="controls",
    children=[
        html.Label("Initial Guess"),
        dcc.Input(id="x0y0", type="text", value="4, -4", style={"width": "120px"}),
        html.Button("Find Fixed Point", id="btn-fp", n_clicks=0),
    ],
)
