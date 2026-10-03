"""
Proof of ``hole-side-and-own-blast``: the six symbolic-itinerary panels per case.

Builds each case with ``henon_cases`` inside one source tree (``--tree``, the
branch by default) and draws ``symbolic_figures.draw_panels`` into
``figures/refactor_proof/hole-side-and-own-blast/panels/<label>/<case>/``,
one folder per case (``build_nested``'s slug ignores the blast counts). A
case that raises gets an ``ERROR.txt`` instead of panels.

Cases (``--cases``, all by default): ``p3_b4`` and ``p3_b6`` (period 3
alone), ``nested_o2_i0`` (the defaults), ``nested_o2_i4``,
``nested_o4_i4``, ``nested_o4_i4_outer_first`` (the inner blasts run after
the outer ones), ``nested_o6_i0``.

Run from the branch worktree::

    MPLBACKEND=Agg ../tangle-pack/env/bin/python scripts/refactor_proof/hole-side-and-own-blast/panels.py
    MPLBACKEND=Agg ../tangle-pack/env/bin/python scripts/refactor_proof/hole-side-and-own-blast/panels.py --tree ../tangle-pack-hole-side-and-own-blast-base --label base --cases nested_o6_i0 nested_o2_i4 p3_b6
"""

from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

BRANCH = Path(__file__).resolve().parents[3]
OUT = BRANCH / "figures" / "refactor_proof" / "hole-side-and-own-blast" / "panels"
CASES = ("p3_b4", "p3_b6", "nested_o2_i0", "nested_o2_i4", "nested_o4_i4", "nested_o4_i4_outer_first", "nested_o6_i0")


def build(hc, case: str):
    """One case, as that tree's recipes build it."""
    if case.startswith("p3_b"):
        return hc.build_period3(blasts=int(case[4:]))
    outer, inner = (int(part[1:]) for part in case.split("_")[1:3])
    if not case.endswith("outer_first"):
        return hc.build_nested(outer_blasts=outer, inner_blasts=inner)
    built = hc.build_nested(outer_blasts=outer)
    session = built.session
    zone = min(session.resonance_zones.values(), key=lambda z: z.area)
    for _ in range(inner):
        session.blast_zone(zone, 1, fixed_point=[zone.fixed_point], min_separation=1e-4)
        hc._repin(session, built.fixed_points, built.pips)
    hc._partition(session)
    built.title += f", then {inner} inner blast(s)"
    return built


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--tree", type=Path, default=BRANCH, help="source tree to run")
    parser.add_argument("--label", default="branch", help="output subfolder (base / branch)")
    parser.add_argument("--cases", nargs="+", default=CASES, choices=CASES)
    args = parser.parse_args()
    tree = args.tree.resolve()
    sys.path[:0] = [str(tree / "src"), str(tree / "scripts")]

    import tanglepack
    from symbolic_figures import configure_logging, draw_panels, report
    from tanglepack.examples import henon_cases as hc

    assert Path(tanglepack.__file__).is_relative_to(tree), tanglepack.__file__
    configure_logging()
    for case in args.cases:
        folder = OUT / args.label / case
        folder.mkdir(parents=True, exist_ok=True)
        try:
            built = build(hc, case)
            dynamics = built.session.symbolic_dynamics()
        except Exception:
            (folder / "ERROR.txt").write_text(traceback.format_exc())
            print(f"{case}: raised, see {folder / 'ERROR.txt'}")
            continue
        title = f"[{args.label}] {built.title}"
        report(built.session, dynamics, title)
        print(f"blast sizes: {built.blast_sizes}\n")
        draw_panels(built.session, built.fixed_points, title, folder, dynamics=dynamics)


if __name__ == "__main__":
    main()
