"""
Nested resonance zones via the TangleSession facade.

Builds the nested Henon tangle (outer period-1 fixed point + inner period-3 orbit),
finds a strong pip on each, then defines and shades one resonance zone per fixed
point. The inner period-3 zone is a subset of the outer period-1 zone, drawn in a
different color.

Run:  PYTHONPATH=src python scripts/henon_resonance_zone.py
"""

import logging

logging.basicConfig(level=logging.WARNING)

import numpy as np
import matplotlib.pyplot as plt

from tanglepack import TangleSession


def henon_map(point):
    k, b = 2, 1
    x, y = point
    return np.array([y - k + x**2, -b * x])


def henon_map_inverse(point):
    k, b = 2, 1
    x, y = point
    return np.array([-y / b, x + k - (y**2) / (b**2)])


def henon_jacobian(point):
    k, b = 2, 1
    x, y = point
    return np.array([[2 * x, 1], [-b, 0]])


# One session ties the numerical engine (workbench) and topology (trellises) together.
session = TangleSession(henon_map, henon_map_inverse, henon_jacobian)

# --- inner period-3 orbit -------------------------------------------------------
session.workbench._man_machine.area_cutoff = 1e-7
fp3 = session.construct_fixed_point([[0, 1], [-1, 0], [-1, 1]])
session.orient_eigenvectors(
    fp3, {"unstable": np.array([0, -1]), "stable": np.array([-1, -1])}
)
session.initialize_both_manifolds(fp3)
session.grow_n_times(fp3, "unstable", num_iterations=10)
session.grow_n_times(fp3, "stable", num_iterations=6)

# --- outer period-1 fixed point -------------------------------------------------
session.workbench._man_machine.area_cutoff = 1e-4
fp1 = session.construct_fixed_point([4, -4])
session.orient_eigenvectors(
    fp1, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
)
session.initialize_both_manifolds(fp1)
session.grow_n_times(fp1, "unstable", num_iterations=7)
session.grow_until_turnaround(fp1, "stable")

# --- co-compute the nested tangle ----------------------------------------------
session.compute_intersections([fp3, fp1])
session.trim_stable_manifolds(fp3)
session.trim_stable_manifolds(fp1)
session.create_bridges(fp3)
session.create_bridges(fp1)
session.infer_iterate_table()

# --- pick a strong pip on each tangle ------------------------------------------
T1 = session.trellis(fp1)
T1.classify_strong_pips()
T3 = session.trellis(fp3)
T3.classify_strong_pips()
print("fp1 strong pip:", T1.strong_pip, "-> cut points:", T1.strong_pip_cut_points())
print("fp3 strong pip:", T3.strong_pip, "-> cut points:", T3.strong_pip_cut_points())
# The outer period-1 anchor cuts at a single pip; the inner period-3 anchor cuts at
# the strong pip AND its two iterates — one point on each of the three stable branches.

T1.set_strong_pip(7)

# --- define one resonance zone per fixed point (one final recompute) -----------
# Both pip ids come from the same registry; add_resonance_zones trims every stable
# branch (each cut point) then recomputes once, so the zones nest correctly: the
# inner period-3 zone is a subset of the outer period-1 zone.
session.add_resonance_zones([T1.strong_pip, T3.strong_pip])
for (fp, bi), rz in session.resonance_zones.items():
    name = "fp1 (period 1)" if fp is fp1 else "fp3 (period 3)"
    print(
        f"  zone {name}: {len(rz.cut_intersections)} cut point(s), "
        f"{len(rz.previous_tails)} trimmed branch(es)"
    )

# --- plot: shaded zones under the tangle ---------------------------------------
plt.figure(figsize=(8, 8))
session.plot_resonance_zones(alpha=0.3)

for fp in (fp3, fp1):
    session.workbench.plot_tangle(fp, "stable", color="r", linewidth=0.8)
    session.workbench.plot_tangle(fp, "unstable", color="b", linewidth=0.8)

# Rebuild trellises (the recompute rebuilt the registry) before plotting pips.
T1 = session.trellis(fp1, rebuild=True)
T1.classify_strong_pips()
T3 = session.trellis(fp3, rebuild=True)
T3.classify_strong_pips()
T1.plot_strong_pip(s=40)
T1.plot_strong_pip_candidates()
T3.plot_strong_pip(s=40)
T3.plot_strong_pip_candidates()

session.workbench.plot_intersections(fp1, show_ids=True)

plt.xlim([-8, 8])
plt.ylim([-8, 8])
plt.title("Nested resonance zones (outer period-1, inner period-3)")
plt.show()
