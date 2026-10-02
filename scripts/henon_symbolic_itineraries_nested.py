"""
Symbolic itineraries of the nested Hénon tangle: the outer period-1 saddle
and the inner period-3 orbit together.

Builds both fixed points of the nested Hénon map (``k = 2``, ``b = 1``) with
``tanglepack.examples.henon_cases.build_nested`` -- the recipe of
``henon_blast_period_3.py``: period 3 as in
``henon_symbolic_itineraries_period3.py``; period 1 from the guess
``[4, -4]``, orientation unstable ``[-1, 0]`` / stable ``[0, 1]``, 11
unstable steps and ``grow_until_turnaround`` on the stable side -- one
resonance zone per fixed point at its default strong pip, blasts the inner
zone (default none) and then the outer zone (default twice, the k=2.8
recipe), partitions everything, and draws the six symbolic-itinerary panels
into ``figures/henon_symbolic_itineraries/nested_p1_p3/`` (PNG, SVG and a
pickled figure each; see ``symbolic_figures.py``).

Element names carry the fixed point's letter (``A`` = the period-3 orbit,
constructed first; ``B`` = the period-1 saddle). The cartoon panel draws the
period-3 zone's circle (three stable arcs) INSIDE the period-1 zone's circle
(one stable arc), because the inner orbit lies in the outer zone.

Run with ``PYTHONPATH=src python3 scripts/henon_symbolic_itineraries_nested.py``
(``--outer-blasts N``, ``--inner-blasts N``, ``--p1-area-cutoff X``; the
``henon_p3_session`` fixture grows the period-1 manifolds at 1e-4).

Dev Notes:
    Without the outer blasts the outer class ``e`` is unresolved (a
    singleton landing with no registered image chain); after two it
    resolves and the outer classes read like k=2.8 two blasts. No
    heteroclinic crossing is computed by this recipe (the two tangles do
    not meet at this growth), so no walk or bridge joins the two circles.
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

from tanglepack.examples.henon_cases import build_nested  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Symbolic itineraries of the nested Hénon tangle."
    )
    parser.add_argument(
        "--outer-blasts", type=int, default=2, help="blasts of the outer zone"
    )
    parser.add_argument(
        "--inner-blasts", type=int, default=0, help="blasts of the inner zone"
    )
    parser.add_argument(
        "--p1-area-cutoff",
        type=float,
        default=1e-7,
        help="area cutoff while growing the period-1 manifolds",
    )
    args = parser.parse_args()
    configure_logging()

    built = build_nested(
        outer_blasts=args.outer_blasts,
        inner_blasts=args.inner_blasts,
        p1_area_cutoff=args.p1_area_cutoff,
    )
    session = built.session
    dynamics = session.symbolic_dynamics()
    report(session, dynamics, built.title)
    if built.blast_sizes:
        print(f"blast sizes (interior bridges per blast): {built.blast_sizes}")
    labels = ", ".join(f"{fp.label} = period {fp.period}" for fp in built.fixed_points)
    print(f"fixed points: {labels}")
    draw_panels(
        session,
        built.fixed_points,
        built.title,
        PANELS_DIR / built.slug,
        dynamics=dynamics,
    )


if __name__ == "__main__":
    main()
