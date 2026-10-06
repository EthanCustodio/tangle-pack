"""P2b: does deeper growth give the inversion saddle pseudoneighbors / holes?

Sweeps the unstable growth depth (stable 5 steps) and reports counts, timing,
anchors per unstable branch, virtual classes and is_reliable.
"""

from __future__ import annotations

import logging
import sys
import time
from collections import Counter

from tanglepack import TangleSession
from tanglepack.examples import HENON_K10, henon_jacobian, henon_map, henon_map_inverse, saddle_guesses

logging.basicConfig(level=logging.ERROR)

for steps in [int(a) for a in sys.argv[1:]] or [6, 7, 8]:
    t0 = time.perf_counter()
    s = TangleSession(henon_map(*HENON_K10), henon_map_inverse(*HENON_K10), henon_jacobian(*HENON_K10))
    fp = s.construct_fixed_point(saddle_guesses(*HENON_K10)["inversion"])
    s.initialize_both_manifolds(fp)
    s.grow_n_times(fp, "unstable", num_iterations=steps)
    s.grow_n_times(fp, "stable", num_iterations=5)
    s.compute_intersections([fp]); s.trim_stable_manifolds(fp); s.create_bridges(fp); s.infer_iterate_table()
    s.classify_strong_pips(); s.compute_pseudoneighbors(); s.punch_holes(); s.partition_stable_manifold()
    try:
        dyn = s.symbolic_dynamics(); rel = dyn.is_reliable; virt = list(dyn.virtual_classes.values()); ncls = len(dyn.classes)
    except Exception as exc:  # noqa: BLE001
        rel = virt = ncls = f"{type(exc).__name__}"
    T = s.trellis(fp)
    anchors = Counter(T.intersection(i).manifold_a_key[2:] for i in s.workbench.intersection_registry.all_ids() if T.intersection(i).label == "anchor")
    print(f"unstable={steps}: {time.perf_counter()-t0:.2f}s crossings={len(T.intersection_ids)} bridges={len(T.bridges)} "
          f"pn={len(T.pseudoneighbors)} holes={len(T.holes)} classes={ncls} virtual={virt} reliable={rel} anchors/branch={dict(anchors)}")
