"""
Symbolic itineraries of a Hénon saddle with inversion.

Builds the ``k = 4``, ``b = 1`` Hénon map's saddle at ``(1 - sqrt5,
sqrt5 - 1)``, whose eigenvalues are both negative (``k_value = 2``: the two
branches of each manifold are one chain, swapped by the map), with
``tanglepack.examples.henon_cases.build_inversion`` -- the recipe of the
inversion law cases: 11 unstable and 11 stable steps at ``area_cutoff =
1e-5``, the zone trimmed at ``f(q0)``, optionally blasted -- and draws the six
symbolic-itinerary panels into
``figures/henon_symbolic_itineraries/inversion[_blastN]/`` (PNG, SVG and a
pickled figure each; see ``symbolic_figures.py``).

The cartoon panel is the figure eight of the inversion zone: one anchor with
two stable branches ``0.0`` and ``0.1`` (D4: each with its own two sides).

Run with ``PYTHONPATH=src python3 scripts/henon_symbolic_itineraries_inversion.py``
(``--blasts N``).

Dev Notes:
    Step 12 of either manifold escapes into the point budget, so the growth
    stays at 11. The class count rises 6, 7, 8, 9, 9, 9 with the blasts and
    drops at 6 blasts by a genuine reference switch (BCG-8, DEC-25, open).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from symbolic_figures import (
    PANELS_DIR,
    configure_logging,
    draw_panels,
    report,
)  # noqa: E402

from tanglepack.examples.henon_cases import build_inversion  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Symbolic itineraries of the k=4 Hénon saddle with inversion."
    )
    parser.add_argument(
        "--blasts", type=int, default=0, help="single-iteration zone blasts"
    )
    args = parser.parse_args()
    configure_logging()

    built = build_inversion(blasts=args.blasts)
    session = built.session
    dynamics = session.symbolic_dynamics()
    report(session, dynamics, built.title)
    if built.blast_sizes:
        print(f"blast sizes (interior bridges per blast): {built.blast_sizes}")
    draw_panels(
        session,
        built.fixed_points,
        built.title,
        PANELS_DIR / built.slug,
        dynamics=dynamics,
    )


if __name__ == "__main__":
    main()
