"""
Symbolic itineraries of the inner period-3 Hénon orbit, alone.

Builds the period-3 saddle orbit of the nested Hénon map (``k = 2``,
``b = 1``) with ``tanglepack.examples.henon_cases.build_period3`` -- the
parameters of ``henon_blast_period_3.py`` and the ``henon_p3_session``
fixture: guess ``[[0, 1], [-1, 0], [-1, 1]]``, orientation unstable
``[0, -1]`` / stable ``[-1, -1]``, ``area_cutoff = 1e-7``, 13 unstable and 9
stable steps -- trims it at its default strong pip (a resonance zone with one
boundary pip per branch), optionally blasts the zone, partitions, and draws
the six symbolic-itinerary panels into
``figures/henon_symbolic_itineraries/p3[_blastN]/`` (PNG, SVG and a pickled
figure each; see ``symbolic_figures.py``).

The cartoon panel is the hexagon of the period-3 zone: three stable arcs
``0.0``, ``1.0``, ``2.0`` in forward-map order, alternating with three
empty (unstable) arcs.

Run with ``PYTHONPATH=src python3 scripts/henon_symbolic_itineraries_period3.py``
(``--blasts N``, ``--unstable-steps N``, ``--stable-steps N``).

Dev Notes:
    At the default growth the trellis is CLOSED: every active hole bridge's
    image is already a hole bridge on the next branch, so the words read
    ``a -> b``, ``b -> c``, ``c -> a u^-1 w^-1`` and blasting changes
    nothing topologically. Growing to 16 unstable steps, or blasting six
    times at the default growth (four and five are fine), currently trips
    the hole ``bridge_side`` invariant (see ``henon_cases.py``).
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

from tanglepack.examples.henon_cases import build_period3  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Symbolic itineraries of the period-3 Hénon orbit."
    )
    parser.add_argument(
        "--blasts", type=int, default=0, help="single-iteration zone blasts"
    )
    parser.add_argument(
        "--unstable-steps", type=int, default=13, help="unstable growth steps"
    )
    parser.add_argument(
        "--stable-steps", type=int, default=9, help="stable growth steps"
    )
    args = parser.parse_args()
    configure_logging()

    built = build_period3(
        unstable_steps=args.unstable_steps,
        stable_steps=args.stable_steps,
        blasts=args.blasts,
    )
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
