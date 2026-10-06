"""Probe-local builders: the conftest recipes (k10, k28 one/two blasts, inversion).

Copied verbatim from tests/conftest.py at baseline so the probes need no pytest.
"""

from __future__ import annotations

import time

import numpy as np

from tanglepack import TangleSession
from tanglepack.examples import HENON_K10, henon_jacobian, henon_map, henon_map_inverse, saddle_guesses

HENON_K28 = (2.8, 1)


def build_k10(unstable_steps: int = 9) -> tuple[TangleSession, object]:
    """The conftest ``k10_partitioned`` recipe (9 unstable steps, infer on)."""
    session = TangleSession(henon_map(*HENON_K10), henon_map_inverse(*HENON_K10))
    fp = session.construct_fixed_point([4, -4])
    session.orient_eigenvectors(fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])})
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=unstable_steps)
    session.grow_until_turnaround(fp, "stable")
    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    session.infer_iterate_table()
    session.classify_strong_pips()
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()
    return session, fp


def build_k28(blasts: int) -> tuple[TangleSession, object]:
    """The conftest ``k28_partitioned`` / ``k28_two_blasts_partitioned`` recipe."""
    session = TangleSession(henon_map(*HENON_K28), henon_map_inverse(*HENON_K28), henon_jacobian(*HENON_K28))
    session.workbench._man_machine.area_cutoff = 1e-7
    fp = session.construct_fixed_point([4, -4])
    session.orient_eigenvectors(fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])})
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=10)
    session.grow_until_turnaround(fp, "stable")
    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    session.infer_iterate_table()
    session.classify_strong_pips()
    trellis = session.trellis(fp)
    pip = trellis.iterate(trellis.strong_pip, 1)
    zone = session.resonance_zone(pip)
    session.classify_strong_pips()
    session.set_strong_pip(fp, pip)
    for _ in range(blasts):
        session.blast_zone(zone, num_iterations=1, fixed_point=[fp], min_separation=1e-5)
        session.classify_strong_pips()
        session.set_strong_pip(fp, pip)
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()
    return session, fp


def timed(fn, *args, **kwargs):
    """Call ``fn`` and return ``(result, seconds)``."""
    t0 = time.perf_counter()
    out = fn(*args, **kwargs)
    return out, time.perf_counter() - t0
