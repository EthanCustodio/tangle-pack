"""
Blast the inner period-3 resonance zone of the nested Hénon tangle.

Builds the same nested tangle as ``henon_bridge_classification.py`` (outer period-1
fixed point + inner period-3 orbit), defines one resonance zone per fixed point, then
"blasts" the inner period-3 zone: every un-iterated bridge inside the zone is iterated
forward, the children that stay inside the zone are kept, and the process repeats. This
exercises iterating *interior* bridges many times — the case that the cdist-precision
hardening makes stable.

The script prints the surviving interior-frontier size per generation (it stays bounded
rather than exploding) and shades the zone with every interior bridge produced.

Run:  PYTHONPATH=src python scripts/henon_blast_period_3.py
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


# --- build the nested tangle (identical to henon_bridge_classification.py) -------
session = TangleSession(henon_map, henon_map_inverse, henon_jacobian)

session.workbench._man_machine.area_cutoff = 1e-7
fp3 = session.construct_fixed_point([[0, 1], [-1, 0], [-1, 1]])
session.orient_eigenvectors(
    fp3, {"unstable": np.array([0, -1]), "stable": np.array([-1, -1])}
)
session.initialize_both_manifolds(fp3)
session.grow_n_times(fp3, "unstable", num_iterations=10)
session.grow_n_times(fp3, "stable", num_iterations=6)

session.workbench._man_machine.area_cutoff = 1e-4
fp1 = session.construct_fixed_point([4, -4])
session.orient_eigenvectors(
    fp1, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
)
session.initialize_both_manifolds(fp1)
session.grow_n_times(fp1, "unstable", num_iterations=7)
session.grow_until_turnaround(fp1, "stable")

session.compute_intersections([fp3, fp1])
session.trim_stable_manifolds(fp3)
session.trim_stable_manifolds(fp1)
session.create_bridges(fp3)
session.create_bridges(fp1)
session.infer_iterate_table()

T1 = session.trellis(fp1)
T1.classify_strong_pips()
T3 = session.trellis(fp3)
T3.classify_strong_pips()
T1.set_strong_pip(7)
session.add_resonance_zones([T1.strong_pip, T3.strong_pip])


# --- blast the inner period-3 zone ----------------------------------------------
inner_zone = max(session.resonance_zones.values(), key=lambda z: z.area)
result = session.blast_zone(inner_zone, num_iterations=1, fixed_point=fp1)

print("\nBlast of the inner period-3 resonance zone")
print("=" * 60)
print(f"requested iterations:      {result.num_iterations_requested}")
print(f"completed iterations:      {result.completed_iterations}")
print(f"terminated early:          {result.terminated_early}")
print(f"bridges skipped:           {result.skipped}")
print(f"max interior depth:        {result.max_depth_reached()}")
print(f"distinct interior bridges: {len(result.all_interior_bridges())}")
print("\ninterior frontier size by generation:")
print("  " + ", ".join(str(len(f)) for f in result.interior_bridges_by_iteration))


# --- plot -----------------------------------------------------------------------
plt.figure(figsize=(8, 8))
session.plot_resonance_zones(alpha=0.3)

for fp in (fp3, fp1):
    session.workbench.plot_tangle(fp, "stable", color="r", linewidth=0.8)
    session.workbench.plot_tangle(fp, "unstable", color="b", linewidth=0.6)

for bridge in result.all_interior_bridges():
    pts = bridge.get_point_array()
    if pts is not None and len(pts):
        pts = np.asarray(pts)
        plt.plot(pts[:, 0], pts[:, 1], color="orange", linewidth=1.0)

plt.xlim([-8, 8])
plt.ylim([-8, 8])
plt.title("Blasted inner period-3 resonance zone (interior bridge images)")
plt.show()
