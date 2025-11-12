from dash import html, dcc
from .stores import stores
from .sections.system_section import system_section
from .sections.fixed_point_section import fixed_point_section
from .sections.initialization import initializor
from .sections.manifold_section import manifold_section
from .sections.bridge_section import bridge_section
from .sections.click_mode import click_mode
from .components.plot_area import plot_area
from .sections.bridge_ops_section import bridge_ops_section


def build_layout():
    return html.Div(
        [
            dcc.Location(id="url"),
            *stores,
            html.Div(
                [
                    html.H3("Tangle Workbench"),
                    system_section,
                    fixed_point_section,
                    initializor,
                    manifold_section,
                    click_mode,
                    bridge_section,
                    bridge_ops_section,
                    html.Div(id="status", className="status"),
                    html.Pre(id="debug", className="debug"),
                ],
                className="sidebar",
            ),
            html.Div(plot_area),
        ],
        className="app",
    )
