"""
Reopen a pickled panel (``<panel>.fig.pickle``) in an interactive window.

The symbolic-itinerary scripts pickle every panel next to its PNG and SVG
(``figures/henon_symbolic_itineraries/<case>/<panel>.fig.pickle``). This
loads one or more of them under the Qt backend (PySide6, in ``env/``) and
shows them with the usual toolbar: pan, zoom, save. Text keeps its point
size while zooming, so zooming in on a crowded cartoon spreads the names
apart.

Run with ``env/bin/python scripts/show_figure.py path/to/cartoon.fig.pickle [...]``.

Note:
    A pickled figure loads only under the matplotlib version that wrote it.
"""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import matplotlib

matplotlib.use("QtAgg")

import matplotlib.pyplot as plt  # noqa: E402


def load(path: Path):
    """
    Unpickle one figure and hand it to the interactive backend.

    Args:
        path: The ``.fig.pickle`` file.

    Returns:
        The figure.
    """
    with open(path, "rb") as handle:
        fig = pickle.load(handle)
    # A figure pickled under Agg comes back with a non-interactive canvas;
    # give it a managed window.
    manager = plt.figure().canvas.manager
    plt.close(manager.canvas.figure)
    manager.canvas.figure = fig
    fig.set_canvas(manager.canvas)
    manager.set_window_title(f"{path.parent.name}/{path.name}")
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description="Show pickled matplotlib figures interactively.")
    parser.add_argument("paths", nargs="+", type=Path, help=".fig.pickle files")
    args = parser.parse_args()
    for path in args.paths:
        load(path)
    plt.show()


if __name__ == "__main__":
    main()
