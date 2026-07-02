"""Nested (period-1 outer + period-3 inner) Hénon strong-pip example.

Demonstrates the one-call session-level strong-pip API: instead of building a
Trellis per fixed point and repeating classify/plot for each, a single
``session.classify_strong_pips()`` classifies every tangle and a single
``session.plot_strong_pip_candidates()`` / ``session.plot_strong_pip()`` draws
them all. Pass a fixed point to either to scope it to one tangle.
"""

import logging

logging.basicConfig(
    level=logging.WARNING,
    format="%(asctime)s %(levelname)-8s %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)

import numpy as np
import matplotlib.pyplot as plt

from tanglepack import TangleSession


def henon_map(point):
    k, b = 2, 1
    x, y = point
    return np.stack([y - k + x**2, -b * x], axis=0)


def henon_map_inverse(point):
    k, b = 2, 1
    x, y = point
    return np.stack([-y / b, x + k - (y**2) / (b**2)], axis=0)


def henon_jacobian(point):
    k, b = 2, 1
    x, y = point
    return np.array([[2 * x, 1], [-b, 0]])


session = TangleSession(henon_map, henon_map_inverse, henon_jacobian)
session.workbench._man_machine.area_cutoff = 1e-7

# Period-3 inner fixed point.
fp3 = session.construct_fixed_point([[0, 1], [-1, 0], [-1, 1]])
session.orient_eigenvectors(
    fp3, {"unstable": np.array([0, -1]), "stable": np.array([-1, -1])}
)
session.initialize_both_manifolds(fp3)
session.grow_n_times(fp3, "unstable", num_iterations=14)
session.grow_n_times(fp3, "stable", num_iterations=9)

# Period-1 outer fixed point.
fp1 = session.construct_fixed_point([4, -4])
session.orient_eigenvectors(
    fp1, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
)
session.initialize_both_manifolds(fp1)
session.grow_n_times(fp1, "unstable", num_iterations=12)
session.grow_n_times(fp1, "stable", num_iterations=9)

# Nested tangles: co-index BOTH fixed points into one Tangle/registry so the
# inner (period-3) and outer (period-1) crossings coexist. Computing them
# separately would wipe the first tangle, since each compute_intersections()
# rebuilds the registry.
session.compute_intersections([fp3, fp1])
session.trim_stable_manifolds(fp3)
session.trim_stable_manifolds(fp1)
session.create_bridges(fp3)
session.create_bridges(fp1)
session.infer_iterate_table()

# One call classifies every tangle (inner and outer); returns {fixed_point: ids}.
candidates = session.classify_strong_pips()
for fp, ids in candidates.items():
    print(f"Strong-pip candidates for {fp}: {ids}")
    print(f"  chosen strong pip: {session.strong_pip(fp)}")

plt.figure()
for fp in (fp3, fp1):
    session.plot_tangle(fp, "stable", color="r", linewidth=1)
    session.plot_tangle(fp, "unstable", color="b", linewidth=1)
    session.plot_intersections(fp, show_ids=True)

# One call each draws candidates (magenta) and the chosen strong pips (green)
# for every tangle. Pass a fixed point to scope it, e.g.
# session.plot_strong_pip_candidates(fp3).
session.plot_strong_pip_candidates()
session.plot_strong_pip(s=12)

plt.xlim([-15, 15])
plt.ylim([-15, 15])
plt.show()
