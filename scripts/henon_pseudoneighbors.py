"""Pseudoneighbors and the stable-manifold partition on the k=10 horseshoe.

Runs the full pipeline of the Pseudoneighbor / Stable Manifold Partition
algorithms on the single-saddle binary-horseshoe Hénon tangle:

1. compute the reference pseudoneighbor pairs on W^S(r_n, r_{n+p}) and extend
   them into full trajectories through the iterate table;
2. punch a hole in each pair's bounded region and backward-propagate the
   reference bridges until they map onto themselves;
3. partition the stable manifold left/right (and interior/exterior at the
   strong-pip resonance-zone cut) by those holes.

Figure 1 shows the tangle with pseudoneighbors (orange) and holes (one
marker+color per orbit, labelled by iterate); figure 2 is the number-line
view of the partition.
"""

import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)

import numpy as np
import matplotlib.pyplot as plt

from tanglepack import TangleSession


def henon_map(point):
    k, b = 2.8, 1
    x, y = point
    return np.stack([y - k + x**2, -b * x], axis=0)


def henon_map_inverse(point):
    k, b = 2.8, 1
    x, y = point
    return np.stack([-y / b, x + k - (y**2) / (b**2)], axis=0)


def henon_jacobian(point):
    k, b = 2.8, 1
    x, y = point
    return np.array([[2 * x, 1], [-b, 0]])


session = TangleSession(henon_map, henon_map_inverse, henon_jacobian)
session.workbench._man_machine.area_cutoff = 1e-7

fp = session.construct_fixed_point([4, -4])
session.orient_eigenvectors(
    fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
)
session.initialize_both_manifolds(fp)
session.grow_n_times(fp, "unstable", num_iterations=10)
session.grow_until_turnaround(fp, "stable")
session.compute_intersections([fp])
session.trim_stable_manifolds(fp)
session.create_bridges(fp)

trellis = session.trellis(fp)

# The chosen strong pip anchors the whole topological analysis. Classify the
# candidates, choose one with set_strong_pip, then initialize the resonance
# zone at it: the zone trims the stable manifold at the pip and recomputes,
# so the crossings of the discarded tail leave the registry entirely. This is
# the necessary first step before partitioning — the zone is what defines
# left and right of the stable manifold.
trellis.classify_strong_pips()
print(f"Strong-pip candidates: {trellis.strong_pip_candidates}")
trellis.set_strong_pip(1)
pip = trellis.strong_pip
print(f"Chosen strong pip: {pip}")

session.add_resonance_zones([pip])  # trim at the pip + recompute (ids preserved)
zone = session.resonance_zones[(fp, 0)]
trellis = session.trellis(fp)  # fresh snapshot of the truncated trellis
trellis.classify_strong_pips(choose_default=False)
trellis.set_strong_pip(5)

# Blasting registers the children's new stable-manifold crossings, so the
# trellis snapshot must be refreshed (and the pip re-established) afterward.
session.blast_zone(zone, num_iterations=2, fixed_point=[fp], min_separation=1e-3)
trellis = session.trellis(fp)
trellis.classify_strong_pips(choose_default=False)
trellis.set_strong_pip(5)

# 1. Reference pseudoneighbors + full trajectories.
reference_pseudoneighbors = session.compute_pseudoneighbors(fp, verbose=True)

# 2. Holes: one per pair, plus the backward-propagated ones. The zone's
# containment test flags holes that belong to a different resonance zone
# (interior=False) so the partition can ignore them.
holes = trellis.punch_holes(in_zone=zone.contains_point, verbose=True)

# 3. The partition of each stable branch, both sides.
stable_partition = trellis.partition_stable_manifold(verbose=True)

# Figure 1: the tangle with pseudoneighbors and holes. plot_tangle draws the
# grown manifold linked-lists; the blasted images live in the bridge children,
# so plot_all_bridges draws the unstable manifold's full computed extent
# (same pattern as henon_blast_period_3.py).
plt.figure(figsize=(9, 9))
session.plot_tangle(fp, "stable", color="r", linewidth=1)
session.plot_tangle(fp, "unstable", color="b", linewidth=1)
session.workbench.plot_all_bridges()
session.plot_intersections(fp, show_ids=True)
session.plot_pseudoneighbors(fp, s=40)
session.plot_holes(fp)
plt.title("Pseudoneighbors (orange) and holes (marker+color = orbit, label = iterate)")

# Frame the tangle (the unstable manifold escapes to infinity).
points = np.array(
    [ix.coords for _, ix in trellis.registry] + [h.coords for h in trellis.holes]
)
pad = 0.25 * np.ptp(points, axis=0).max() + 1.0
plt.xlim(points[:, 0].min() - pad, points[:, 0].max() + pad)
plt.ylim(points[:, 1].min() - pad, points[:, 1].max() + pad)

# Figure 2: the stable-manifold partition as a number line.
fig, ax = plt.subplots(figsize=(9, 3))
trellis.plot_stable_partition(ax=ax)
fig.tight_layout()

plt.show()
