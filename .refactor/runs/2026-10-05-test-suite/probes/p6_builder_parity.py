"""P6: builder parity -- repoint conftest fixtures at the planned builders and see what fails.

Each variant rewrites tests/conftest.py TEMPORARILY (restored from the pristine
text in a ``finally``), runs the suite (no coverage) and records the failing
node ids and wall time. Variants:

  A   henon_tangle_with_bridges: 9 unstable steps in one go + infer_iterate_table
      (planned build_k10(through="bridges")) instead of 7 (+turnaround) + 2, no infer.
  B4  henon_p3_session: henon_cases.build_nested(outer_blasts=0, p1_area_cutoff=1e-4)
      (adds _repin + partition) instead of the hand recipe without _repin.
  B7  as B4 with p1_area_cutoff=1e-7 (the henon_cases default).
  C   k28_partitioned / k28_two_blasts_partitioned function-scoped.

Usage: env/bin/python p6_builder_parity.py [A B4 B7 C]
"""

from __future__ import annotations

import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
CONFTEST = ROOT / "tests" / "conftest.py"
OUT = Path(__file__).resolve().parent

A_OLD = '''    workbench, fp = grown_both
    workbench.grow_n_times(fp, "unstable", num_iterations=2)
    workbench.compute_intersections([fp])
    workbench.trim_stable_manifolds(fp)
    workbench.create_bridges(fp)
    return workbench, fp'''
A_NEW = '''    workbench, fp = initialized
    workbench.grow_n_times(fp, "unstable", num_iterations=9)
    workbench.grow_until_turnaround(fp, "stable")
    workbench.compute_intersections([fp])
    workbench.trim_stable_manifolds(fp)
    workbench.create_bridges(fp)
    workbench.infer_iterate_table()
    return workbench, fp'''


def variant_text(name: str, text: str) -> str:
    """Return the conftest text for one variant."""
    if name == "A":
        assert A_OLD in text
        text = text.replace("def henon_tangle_with_bridges(grown_both):", "def henon_tangle_with_bridges(initialized):")
        return text.replace(A_OLD, A_NEW)
    if name in ("B4", "B7"):
        cutoff = "1e-4" if name == "B4" else "1e-7"
        start = text.index("    session = TangleSession(_p3_map, _p3_map_inverse, _p3_jacobian)")
        end = text.index("    return session, fp3, fp1, inner_zone")
        body = (
            "    from tanglepack.examples.henon_cases import build_nested\n"
            f"    build = build_nested(outer_blasts=0, inner_blasts=0, p1_area_cutoff={cutoff})\n"
            "    session = build.session\n"
            "    fp1, fp3 = build.fixed_points\n"
            "    inner_zone = max(session.resonance_zones.values(), key=lambda z: z.area)\n"
        )
        return text[:start] + body + text[end:]
    if name == "C":
        return text.replace('@pytest.fixture(scope="session")\ndef k28', "@pytest.fixture\ndef k28")
    raise ValueError(name)


def run(name: str, pristine: str) -> None:
    """Run one variant and write its failures."""
    CONFTEST.write_text(variant_text(name, pristine))
    t0 = time.perf_counter()
    proc = subprocess.run(
        [str(ROOT / "env/bin/python"), "-m", "pytest", "-q", "-p", "no:cacheprovider", "-rfE", "--no-header"],
        cwd=ROOT, capture_output=True, text=True,
    )
    wall = time.perf_counter() - t0
    failed = sorted(set(re.findall(r"^(?:FAILED|ERROR) (\S+)", proc.stdout, re.M)))
    tail = proc.stdout.strip().splitlines()[-1] if proc.stdout.strip() else proc.stderr[-300:]
    with open(OUT / "p6_output.txt", "a") as fh:
        fh.write(f"\n== variant {name}: {tail} (wall {wall:.0f}s)\n")
        for nid in failed:
            fh.write(f"   {nid}\n")
    print(name, tail, f"{wall:.0f}s", len(failed), "failing")


if __name__ == "__main__":
    pristine = CONFTEST.read_text()
    try:
        for name in sys.argv[1:] or ["A", "B4", "B7", "C"]:
            run(name, pristine)
    finally:
        CONFTEST.write_text(pristine)
