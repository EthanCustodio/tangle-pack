"""
Symbolic itineraries of the bridge classes on the period-1 Hénon cases of
``henon_bridge_classes.py``: the words the dual-graph walks read off.

One folder per case (k=10; k=2.8 blasted once; k=2.8 blasted twice) under
``figures/henon_symbolic_itineraries/<case>/``, six panels each, every panel
written as PNG, SVG and a pickled figure (see ``symbolic_figures.py`` for
the panels; reopen a pickle with ``scripts/show_figure.py`` to pan and zoom).
The higher-period and nested cases have their own drivers,
``henon_symbolic_itineraries_period3.py`` and
``henon_symbolic_itineraries_nested.py``.

The cases are imported from ``henon_bridge_classes.py`` (same recipes, same
``build``) and the plane panels are framed as in
``henon_minimal_dual_graph.py`` (the folds near the fixed point).

Run with ``PYTHONPATH=src python3 scripts/henon_symbolic_itineraries.py``
(``--overview`` also writes the old single 2x3 figure per case, line cartoon,
to ``figures/henon_symbolic_itineraries_<case>.png``); each case's symbolic
dynamics report is printed.

Dev Notes:
    The one-blast k=2.8 case has its single active class resolved through
    the singleton/trellis path (both of its elements are singletons and own
    no dual-graph node); the printed report shows ``source = trellis`` for
    it. The cartoon panel's text box appears only when no bridge was drawn.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from henon_bridge_classes import CASES, Case, build  # noqa: E402
from symbolic_figures import (  # noqa: E402
    FIGURES_DIR,
    PANELS_DIR,
    configure_logging,
    draw_overview,
    draw_panels,
    report,
)


def case_slug(case: Case) -> str:
    """``k10``, ``k28_blast1``, ...: the folder name of a case."""
    k_text = f"{case.k:g}".replace(".", "")
    blasted = f"_blast{case.blasts}" if case.blasts else ""
    return f"k{k_text}{blasted}"


def draw_case(case: Case, *, overview: bool = False) -> Path:
    """Build one case, print its report, draw and save its panels."""
    session, fp, _pre_trim_stable_points = build(case)
    dynamics = session.symbolic_dynamics()
    report(session, dynamics, case.title)
    folder = PANELS_DIR / case_slug(case)
    draw_panels(session, [fp], case.title, folder, dynamics=dynamics)
    if overview:
        draw_overview(
            session, [fp], case.title,
            FIGURES_DIR / f"henon_symbolic_itineraries_{case_slug(case)}.png",
            dynamics=dynamics,
        )
    return folder


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument(
        "--overview", action="store_true",
        help="also write the old single 2x3 figure per case",
    )
    args = parser.parse_args()
    configure_logging()
    for case in CASES:
        draw_case(case, overview=args.overview)


if __name__ == "__main__":
    main()
