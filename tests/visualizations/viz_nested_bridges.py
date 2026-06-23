"""Visualize the create_bridges "nested cdist span" bug.

A manifold's bridges should partition its unstable canonical distance: consecutive
bridges share only their straddle endpoint, never a substantial interval. When an
iterated manifold is cut by ``Tangle.create_bridges`` during a blast, however, it
can emit two bridges whose cdist spans overlap by far more than the straddle margin
-- two bridge objects tracing essentially the same piece of curve. This is a
genuine cutting bug (distinct from the iterate-time duplicate the workbench already
canonicalizes away), so it is worth seeing.

This script builds the nested period-3 Hénon tangle (as in
``scripts/henon_blast_period_3.py``), blasts the inner zone, finds every pair of
bridges on one manifold whose cdist spans overlap by more than 10% of the smaller
span, and draws:

  * a full-tangle figure with every offending bridge highlighted, and
  * a zoom on the worst pair, the two bridges drawn in different colours so the
    geometric overlap is visible.

It adds nothing to the core library -- it only reads public bridge data and the
workbench's bridge-signature helpers.

Run:  MPLBACKEND=Agg PYTHONPATH=src python tests/visualizations/viz_nested_bridges.py
"""

from __future__ import annotations

import logging
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

from tanglepack import TangleSession

logging.basicConfig(level=logging.ERROR)

OVERLAP_FRACTION = 0.10  # span overlap above this fraction of the smaller span = bug
OUT_FULL = "/tmp/viz_nested_bridges_full.png"
OUT_ZOOM = "/tmp/viz_nested_bridges_zoom.png"


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


def build_and_blast():
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
    session.add_resonance_zones([T1.strong_pip, T3.strong_pip])

    zone = max(session.resonance_zones.values(), key=lambda z: z.area)
    session.blast_zone(zone, num_iterations=4, fixed_point=fp1)
    return session, (fp3, fp1)


def find_nested_pairs(workbench):
    """Pairs of same-manifold bridges whose cdist spans substantially overlap."""
    by_manifold = defaultdict(list)
    for bridge in workbench.bridges:
        sig = workbench._bridge_signature(bridge)
        if sig is not None:
            by_manifold[workbench._manifold_identity(bridge)].append((sig, bridge))

    pairs = []
    for items in by_manifold.values():
        items.sort(key=lambda it: it[0])
        for i in range(len(items)):
            (a0, a1), bi = items[i]
            for j in range(i + 1, len(items)):
                (b0, b1), bj = items[j]
                if b0 >= a1:  # sorted; no further overlaps with i
                    break
                overlap = min(a1, b1) - max(a0, b0)
                smaller = min(a1 - a0, b1 - b0)
                if smaller > 0 and overlap > OVERLAP_FRACTION * smaller:
                    pairs.append(((a0, a1), bi, (b0, b1), bj, overlap / smaller))
    pairs.sort(key=lambda p: -p[4])  # worst overlap first
    return pairs


def main():
    session, fps = build_and_blast()
    wb = session.workbench
    pairs = find_nested_pairs(wb)

    print(f"\nNested / overlapping bridge pairs found: {len(pairs)}")
    for (a0, a1), _bi, (b0, b1), _bj, frac in pairs[:12]:
        print(
            f"  cdist [{a0:.3f}, {a1:.3f}]  &  [{b0:.3f}, {b1:.3f}]   "
            f"overlap = {frac * 100:.0f}% of the smaller span"
        )
    if not pairs:
        print("  (none found)")
        return

    # --- full tangle with every offending bridge highlighted -----------------
    plt.figure(figsize=(9, 9))
    for fp in fps:
        wb.plot_tangle(fp, "unstable", color="0.75", linewidth=0.3)
    highlighted = set()
    for _sa, bi, _sb, bj, _f in pairs:
        for br in (bi, bj):
            if id(br) in highlighted:
                continue
            highlighted.add(id(br))
            pts = np.asarray(br.get_point_array())
            if len(pts):
                plt.plot(pts[:, 0], pts[:, 1], color="red", linewidth=1.4)
    plt.title(f"Overlapping-cdist bridge pairs ({len(pairs)} pairs, red)")
    plt.xlim([-8, 8])
    plt.ylim([-8, 8])
    plt.gca().set_aspect("equal")
    plt.tight_layout()
    plt.savefig(OUT_FULL, dpi=140)
    print(f"\nsaved full figure -> {OUT_FULL}")

    # --- zoom on the worst pair ---------------------------------------------
    (a0, a1), bi, (b0, b1), bj, frac = pairs[0]
    pi = np.asarray(bi.get_point_array())
    pj = np.asarray(bj.get_point_array())
    allpts = np.vstack([pi, pj])
    pad = 0.05 * max(np.ptp(allpts[:, 0]), np.ptp(allpts[:, 1]), 1e-3)

    plt.figure(figsize=(9, 9))
    plt.plot(pi[:, 0], pi[:, 1], "-o", color="tab:blue", ms=3, lw=1.2,
             label=f"bridge A  cdist [{a0:.2f}, {a1:.2f}]")
    plt.plot(pj[:, 0], pj[:, 1], "-s", color="tab:orange", ms=3, lw=1.2,
             label=f"bridge B  cdist [{b0:.2f}, {b1:.2f}]")
    plt.legend(loc="best", fontsize=9)
    plt.title(f"Worst overlapping pair — {frac * 100:.0f}% span overlap "
              "(two bridges, same curve)")
    plt.xlim(allpts[:, 0].min() - pad, allpts[:, 0].max() + pad)
    plt.ylim(allpts[:, 1].min() - pad, allpts[:, 1].max() + pad)
    plt.gca().set_aspect("equal")
    plt.tight_layout()
    plt.savefig(OUT_ZOOM, dpi=140)
    print(f"saved zoom figure -> {OUT_ZOOM}")

    plt.show()


if __name__ == "__main__":
    main()
