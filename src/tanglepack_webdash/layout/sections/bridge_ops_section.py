from dash import html, dcc

bridge_ops_section = html.Div(
    className="controls bridge-ops-section",
    children=[
        html.Label(
            "Bridge Operations", style={"fontWeight": "bold", "marginBottom": "8px"}
        ),
        html.Div(
            id="selected-bridge-info",
            children="No bridge selected",
            style={"fontSize": "12px", "color": "#9b9ca0", "marginBottom": "8px"},
        ),
        html.Button(
            "Iterate Bridge",
            id="btn-iterate-bridge",
            n_clicks=0,
            className="btn",
            disabled=True,  # Disabled until a bridge is selected
        ),
        html.Button(
            "Clear Selection",
            id="btn-clear-bridge-selection",
            n_clicks=0,
            className="btn",
        ),
    ],
)
