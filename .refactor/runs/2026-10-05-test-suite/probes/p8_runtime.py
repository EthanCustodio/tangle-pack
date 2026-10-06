"""P8: runtime of each law case's full build (through the symbolic dynamics).

Times the partitioned build and then every downstream product, three repeats
per case, in one process (the first repeat carries import/JIT warm-up).
"""

from __future__ import annotations

import logging
import statistics
import sys
import time

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from _builds import build_k10, build_k28  # noqa: E402
from tanglepack import TangleSession  # noqa: E402
from tanglepack.examples import HENON_K10, henon_jacobian, henon_map, henon_map_inverse, saddle_guesses  # noqa: E402
from tanglepack.examples.henon_cases import build_nested, build_period3  # noqa: E402

logging.basicConfig(level=logging.ERROR)


def build_inversion():
    """The inversion saddle, 6 unstable / 5 stable steps, partitioned."""
    s = TangleSession(henon_map(*HENON_K10), henon_map_inverse(*HENON_K10), henon_jacobian(*HENON_K10))
    fp = s.construct_fixed_point(saddle_guesses(*HENON_K10)["inversion"])
    s.initialize_both_manifolds(fp)
    s.grow_n_times(fp, "unstable", num_iterations=6)
    s.grow_n_times(fp, "stable", num_iterations=5)
    s.compute_intersections([fp]); s.trim_stable_manifolds(fp); s.create_bridges(fp); s.infer_iterate_table()
    s.classify_strong_pips(); s.compute_pseudoneighbors(); s.punch_holes(); s.partition_stable_manifold()
    return s


CASES = {
    "k10": lambda: build_k10()[0],
    "k28_one_blast": lambda: build_k28(1)[0],
    "k28_two_blasts": lambda: build_k28(2)[0],
    "p3": lambda: build_period3().session,
    "nested": lambda: build_nested().session,
    "inversion": build_inversion,
}

if __name__ == "__main__":
    print(f"{'case':16s} {'build s (min/median of 3)':>26s} {'downstream s':>13s} {'total':>7s}")
    grand = 0.0
    for name, fn in CASES.items():
        builds, downs = [], []
        for _ in range(3):
            t0 = time.perf_counter(); s = fn(); t1 = time.perf_counter()
            s.arrangement(); s.bridge_classes(); s.minimal_trellis(); s.iterated_partition(); s.dual_graph(); s.symbolic_dynamics()
            t2 = time.perf_counter()
            builds.append(t1 - t0); downs.append(t2 - t1)
        b, d = statistics.median(builds), statistics.median(downs)
        grand += b + d
        print(f"{name:16s} {min(builds):10.2f} / {b:6.2f}         {d:13.2f} {b + d:7.2f}")
        sys.stdout.flush()
    print(f"sum of medians (one full build of every law case): {grand:.2f}s")
