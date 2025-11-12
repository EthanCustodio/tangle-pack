from dash import html, dcc

click_mode = html.Div(
    className="controls",
    children=[
        html.Label("Click mode"),
        dcc.Dropdown(
            id="click-mode",
            options=[
                {
                    "label": "None (default Plotly)",
                    "value": "none",
                },
                {
                    "label": "Select fixed point guess",
                    "value": "fp_guess",
                },
                {
                    "label": "Surface of section seed",
                    "value": "sos_seed",
                },
                {
                    "label": "Select bridge",
                    "value": "select_bridge",
                },
            ],
            value="none",
            clearable=False,
            style={"width": "260px"},
        ),
        html.Button(
            "Clear guesses",
            id="btn-clear-guesses",
            n_clicks=0,
            className="btn",
        ),
    ],
)
