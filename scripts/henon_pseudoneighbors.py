"""Pseudoneighbors and the stable-manifold partition on the k=10 horseshoe.

Runs the full pipeline of the Pseudoneighbor / Stable Manifold Partition
algorithms on the single-saddle binary-horseshoe Hénon tangle:

1. compute the reference pseudoneighbor pairs on W^S(r_n, r_{n+p}) and extend
   them into full trajectories through the iterate table;
2. punch a hole in each pair's bounded region and backward-propagate the
   reference bridges until they map onto themselves — each hole classified
   by its side of its bridge in the dynamical orientation;
3. partition the stable manifold left/right by those holes.

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
from tanglepack.examples import (
    henon_map as _henon_map_factory,
    henon_map_inverse as _henon_map_inverse_factory,
    henon_jacobian as _henon_jacobian_factory,
)


henon_map = _henon_map_factory(2.8, 1)
henon_map_inverse = _henon_map_inverse_factory(2.8, 1)
henon_jacobian = _henon_jacobian_factory(2.8, 1)


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
trellis.set_strong_pip(5)
pip = trellis.strong_pip
print(f"Chosen strong pip: {pip}")

session.add_resonance_zones([pip])  # trim at the pip + recompute (ids preserved)
zone = session.resonance_zones[(fp, 0)]
trellis = session.trellis(fp)  # fresh snapshot of the truncated trellis
trellis.classify_strong_pips(choose_default=False)
trellis.set_strong_pip(5)

# Blasting registers the children's new stable-manifold crossings, so the
# trellis snapshot must be refreshed (and the pip re-established) afterward.
session.blast_zone(zone, num_iterations=13, fixed_point=[fp], min_separation=1e-5)
trellis = session.trellis(fp)
trellis.classify_strong_pips(choose_default=False)
trellis.set_strong_pip(5)

# 1. Reference pseudoneighbors + full trajectories.
reference_pseudoneighbors = session.compute_pseudoneighbors(fp, verbose=True)

# 2. Holes: one per pair, plus the backward-propagated ones. Each hole is
# classified by its side of its bridge in the dynamical orientation
# (Hole.bridge_side); the zone plays no role here — it exists for blasting.
holes = trellis.punch_holes(verbose=True)

# 3. The partition of each stable branch, both sides.
stable_partition = trellis.partition_stable_manifold(verbose=True)

plt.show()
