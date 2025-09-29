# tanglepack_webdash/app.py
from __future__ import annotations
from dash import Dash
from .layout import build_layout
from .callbacks import fixed_point  # registers callbacks on import


def make_app() -> Dash:
    """
    Initializes the main app that contains everything

    Returns:
        Dash: the app itself!
    """

    app = Dash(__name__)
    app.title = "Tangle Workbench"
    app.layout = build_layout()
    # importing modules under .callbacks should register their callbacks with `app`
    fixed_point.register(app)
    return app


def main():
    """Program starting function"""

    app = make_app()
    app.run(debug=True)


# for gunicorn: "tanglepack_webdash.app:server"
app = make_app()
server = app.server

if __name__ == "__main__":
    main()
