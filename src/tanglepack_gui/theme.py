# src/tanglepack_gui/theme.py
from __future__ import annotations
from dataclasses import dataclass
import pyqtgraph as pg


@dataclass(frozen=True)
class PlotTheme:
    name: str
    background: str  # canvas bg
    foreground: str  # axes, text, grid color base
    fp_brush: str  # fixed-point fill
    fp_pen: str  # fixed-point outline


LIGHT = PlotTheme(
    name="light",
    background="w",
    foreground="k",
    fp_brush="k",
    fp_pen="k",
)

DARK = PlotTheme(
    name="dark",
    background="k",
    foreground="w",
    fp_brush="w",
    fp_pen="w",
)


def apply_pyqtgraph_theme(theme: PlotTheme) -> None:
    """Set global pg config for foreground/background so axes, grid, and defaults adapt."""
    pg.setConfigOption("background", theme.background)
    pg.setConfigOption("foreground", theme.foreground)
