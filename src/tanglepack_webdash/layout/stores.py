from dash import dcc, html

stores = [
    dcc.Store(id="sid", storage_type="session"),
    dcc.Store(id="fp-guess-points", storage_type="session", data=[]),
    dcc.Store(id="sos-seeds", storage_type="session", data=[]),
    dcc.Store(id="bridges-data", storage_type="session", data=None),
    dcc.Store(id="selected-bridge-idx", storage_type="session", data=None),
    html.Div(id="bridges-data-sync", style={"display": "none"}),
    html.Div(id="bridges-data-holder", style={"display": "none"}),
]
