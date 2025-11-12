from dash import html, dcc
from dash_extensions import EventListener
from tanglepack_webdash.utils.figures import blank_figure

plot_area = html.Div(
    [
        html.Div(
            [
                # Live cursor readout bar
                html.Div(id="cursor-readout", className="cursor-readout"),
                # text list of clicked points
                html.Pre(
                    id="click-points",
                    style={
                        "margin": 0,
                        "padding": "6px 10px",
                        "whiteSpace": "pre",
                    },
                ),
                # hidden sink for marker sync
                html.Div(id="points-sync", style={"display": "none"}),
                dcc.Interval(id="fp-sync", interval=250, n_intervals=0),
                html.Div(id="mode-sync", style={"display": "none"}),
                # Wrap the Graph with EventListener; emit both mousemove + click
                EventListener(
                    id="plot-events",
                    logging=False,
                    useCapture=True,
                    # target_id="plot",
                    events=[
                        {
                            "event": "mousemove",
                            "props": ["clientX", "clientY", "type"],
                        },
                        {
                            "event": "click",
                            "props": ["clientX", "clientY", "type"],
                        },
                    ],
                    children=dcc.Graph(
                        id="plot",
                        figure=blank_figure(),
                        style={"height": "100%", "width": "100%"},
                        config={"displaylogo": False},
                    ),
                ),
            ],
            id="plot-wrap",
        ),
    ],
    className="main",
)
