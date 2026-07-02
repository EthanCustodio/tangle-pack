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

Set ``TANGLEPACK_GPU=1`` to run the batched map evaluations on the GPU (needs
CuPy; ``pip install tanglepack[gpu]``). Falls back to the CPU automatically if
CuPy is unavailable. The Hénon maps below are written batch-capable (coordinate
on axis 0, built with ``np.stack``) so a single vectorized call maps every point
of a layer/iteration at once -- on the CPU normally, on the GPU when enabled.
"""

import logging
import os

logging.basicConfig(level=logging.WARNING)

import numpy as np
import matplotlib.pyplot as plt

import tanglepack
from tanglepack import TangleSession

# Declare up-front whether this run should use the GPU. Enabling it only changes
# where the batched map is evaluated; the rest of the pipeline is untouched.
USE_GPU = os.environ.get("TANGLEPACK_GPU", "0") == "1"


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


# --- build the nested tangle (identical to henon_bridge_classification.py) -------
session = TangleSession(henon_map, henon_map_inverse, henon_jacobian)

if USE_GPU:
    try:
        tanglepack.enable_gpu(session)
        print("GPU acceleration enabled for the batched map evaluations.")
    except ImportError as exc:
        print(f"GPU requested but unavailable ({exc}); continuing on the CPU.")

session.workbench._man_machine.area_cutoff = 1e-7
fp3 = session.construct_fixed_point([[0, 1], [-1, 0], [-1, 1]])
session.orient_eigenvectors(
    fp3, {"unstable": np.array([0, -1]), "stable": np.array([-1, -1])}
)
session.initialize_both_manifolds(fp3)
session.grow_n_times(fp3, "unstable", num_iterations=13)
session.grow_n_times(fp3, "stable", num_iterations=9)

session.workbench._man_machine.area_cutoff = 1e-7
fp1 = session.construct_fixed_point([4, -4])
session.orient_eigenvectors(
    fp1, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
)
session.initialize_both_manifolds(fp1)
session.grow_n_times(fp1, "unstable", num_iterations=11)
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
# Pin a specific candidate for the outer zone rather than the default (the
# default, smallest unstable cdist, gives a different zone boundary). CAUTION:
# intersection ids are not stable across growth/refinement changes -- if this
# raises, pick again from the candidate list in the error message.
T1.set_strong_pip(10)

session.add_resonance_zones([T1.strong_pip, T3.strong_pip])


# --- blast both resonance zones, bridges from both fixed points -----------------
# The two zones overlap (the small period-3 zone sits inside the large period-1
# zone), and a bridge is iterated at most once, so blast the SMALLER, more specific
# zone first -- otherwise the large blast consumes the small zone's bridges and the
# small blast finds nothing left to iterate.
small_zone = min(session.resonance_zones.values(), key=lambda z: z.area)
large_zone = max(session.resonance_zones.values(), key=lambda z: z.area)

# min_separation stops iterating a bridge once its image folds back onto curve that
# already exists within that distance -- beyond it the folds are iterates-of-iterates
# whose accumulated error makes them cross ("zig-zag"). It keeps the blast on the
# well-resolved side of the precision limit; lower it to resolve finer folds.
result_small = session.blast_zone(
    small_zone, num_iterations=14, fixed_point=[fp1, fp3], min_separation=1e-4
)
result_large = session.blast_zone(
    large_zone, num_iterations=12, fixed_point=[fp1, fp3], min_separation=1e-3
)

for name, result in (("small", result_small), ("large", result_large)):
    print(f"\nBlast of the {name} resonance zone")
    print("=" * 60)
    print(f"completed iterations:      {result.completed_iterations}")
    print(f"terminated early:          {result.terminated_early}")
    print(f"distinct interior bridges: {len(result.all_interior_bridges())}")
    print("interior frontier size by generation:")
    print("  " + ", ".join(str(len(f)) for f in result.interior_bridges_by_iteration))


# --- plot -----------------------------------------------------------------------
# Every bridge the blast computed is registered in the workbench (`session.workbench
# .bridges`), so we plot the whole tangle straight from there -- no need to reach
# into the blast results. `plot_tangle` draws the grown manifold linked-lists (the
# red/blue curves); `plot_all_bridges` draws every registered bridge, original and
# iterated.
plt.figure(figsize=(8, 8))
session.plot_resonance_zones(alpha=0.3)

for fp in (fp3, fp1):
    session.workbench.plot_tangle(fp, "stable", color="r", linewidth=0.8)

session.workbench.plot_all_bridges()

print(f"\nbridges registered in the workbench: {len(session.workbench.bridges)}")

plt.xlim([-8, 8])
plt.ylim([-8, 8])
plt.title("Blasted resonance zones — all bridges from the workbench")
plt.show()
