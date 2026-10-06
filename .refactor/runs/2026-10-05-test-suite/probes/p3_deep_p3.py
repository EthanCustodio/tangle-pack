"""P3: the known-open deep period-3 runs (15 / 16 unstable steps, 6 blasts).

For each run: build time, I1 (check_holes_share_bridge_side), the classes
whose walk is unreachable / unresolved, the virtual ``new*`` symbols,
``is_reliable``, and the per-class itinerary status. These become the exact
assertions of the ``xfail(strict=True)`` deep-p3 tests.
"""

from __future__ import annotations

import logging
import sys
import time

from tanglepack.examples.henon_cases import build_period3
from tanglepack.topology import check_holes_share_bridge_side

logging.basicConfig(level=logging.ERROR)

RUNS = {"15_steps": {"unstable_steps": 15}, "16_steps": {"unstable_steps": 16}, "6_blasts": {"blasts": 6},
        "default_13": {}, "4_blasts": {"blasts": 4}, "5_blasts": {"blasts": 5}}

for name in sys.argv[1:] or list(RUNS):
    kwargs = RUNS[name]
    t0 = time.perf_counter()
    build = build_period3(**kwargs)
    t_build = time.perf_counter() - t0
    (fp3,) = build.fixed_points
    T = build.session.trellis(fp3)
    try:
        check_holes_share_bridge_side(T.holes, orientation_preserving=T.orientation_preserving)
        i1 = "holds"
    except AssertionError as exc:
        i1 = f"FAILS {str(exc)[:120]}"
    t1 = time.perf_counter()
    dyn = build.session.symbolic_dynamics()
    t_dyn = time.perf_counter() - t1
    statuses = {}
    unresolved = []
    for bc, cd in dyn.classes.items():
        status = cd.search.status if cd.search is not None else f"(source {cd.source})"
        statuses[str(cd.symbol if hasattr(cd, 'symbol') else bc)] = status
        if cd.unresolved_reason:
            unresolved.append((str(bc), cd.unresolved_reason[:80]))
    n_inert = sum(1 for e in build.session.bridge_classes().entries if e.inert)
    print(f"{name} {kwargs}: build {t_build:.1f}s, dyn {t_dyn:.2f}s, holes={len(T.holes)}, I1 {i1}")
    print(f"   classes={len(dyn.classes)} (inert {n_inert}), virtual={sorted(dyn.virtual_classes.values())}, is_reliable={dyn.is_reliable}")
    print(f"   walk statuses: {sorted(statuses.values())}")
    print(f"   unresolved: {unresolved}")
    print(f"   rules: {dyn.rules if hasattr(dyn, 'rules') else ''}")
    sys.stdout.flush()
