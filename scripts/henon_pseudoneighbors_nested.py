"""Pseudoneighbors and stable-manifold partitions on the nested period-3 tangle.

The nested (period-1 outer + period-3 inner) Hénon setup of
``henon_strong_pip_nested.py``, taken through the full pseudoneighbor
pipeline for BOTH tangles: strong pips → resonance zones (one batched
recompute) → reference pseudoneighbors → holes (zone-classified) →
left/right partitions of every stable branch.

Each fixed point runs on its own per-fixed-point trellis (the session fans
that out); within a trellis, hole punching and partitioning cover all of its
stable branches automatically — the period-3 trellis partitions all three.

Figure 1 shows both tangles with pseudoneighbors and holes; figure 2 stacks
the partition number lines of every branch and side.
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
from tanglepack.topology import plot_stable_partition


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
session.grow_n_times(fp3, "unstable", num_iterations=10)
session.grow_n_times(fp3, "stable", num_iterations=9)

# Period-1 outer fixed point.
session.workbench._man_machine.area_cutoff = 1e-4
fp1 = session.construct_fixed_point([4, -4])
session.orient_eigenvectors(
    fp1, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
)
session.initialize_both_manifolds(fp1)
session.grow_n_times(fp1, "unstable", num_iterations=11)
session.grow_until_turnaround(fp1, "stable")

# Co-index both tangles into one registry, then cut bridges.
session.compute_intersections([fp3, fp1])
session.trim_stable_manifolds(fp3)
session.trim_stable_manifolds(fp1)
session.create_bridges(fp3)
session.create_bridges(fp1)
session.infer_iterate_table()

# Strong pips (default choice per tangle), then both resonance zones in one
# batched recompute so every pip id stays valid.
session.classify_strong_pips()
pips = {fp: session.strong_pip(fp) for fp in (fp3, fp1)}
print(f"Strong pips: period-3 {pips[fp3]}, period-1 {pips[fp1]}")
session.add_resonance_zones([pips[fp1], pips[fp3]])

# The pseudoneighbor pipeline, per fixed point on the truncated trellises.
partitions = {}
for fp, name in ((fp3, "period-3"), (fp1, "period-1")):
    print(f"\n=== {name} tangle ===")
    trellis = session.trellis(fp)
    trellis.classify_strong_pips(choose_default=False)
    trellis.set_strong_pip(pips[fp])
    zone = session.resonance_zones[(fp, 0)]

    trellis.compute_pseudoneighbors(verbose=True)
    trellis.punch_holes(in_zone=zone.contains_point, verbose=True)
    partitions[name] = trellis.partition_stable_manifold(verbose=True)

# Figure 1: both tangles with pseudoneighbors and holes.
plt.figure(figsize=(9, 9))
for fp in (fp3, fp1):
    session.plot_tangle(fp, "stable", color="r", linewidth=1)
    session.plot_tangle(fp, "unstable", color="b", linewidth=1)
    session.plot_intersections(fp, show_ids=True)
session.plot_pseudoneighbors(compute=False, s=40)
session.plot_holes()
plt.title(
    "Nested tangles: pseudoneighbors (orange), holes (marker+color = orbit, "
    "label = iterate)"
)
plt.xlim([-15, 15])
plt.ylim([-15, 15])

# Figure 2: every branch's left/right partition as stacked number lines —
# one subplot per tangle, since their canonical-distance scales differ.
fig, axes = plt.subplots(
    len(partitions),
    1,
    figsize=(10, 1.4 * sum(len(p) for p in partitions.values())),
)
for axis, (name, results) in zip(np.atleast_1d(axes), partitions.items()):
    plot_stable_partition(results, ax=axis)
    axis.set_title(f"{name} stable manifold partition")
fig.tight_layout()

plt.show()
