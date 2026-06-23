"""Overlay the unstable x unstable crossing artifacts on the blast figure.

Two unstable manifolds can never cross (fundamental invariant), so any
unstable x unstable segment pair the Tangle flags is a numerical artifact -- a
near-tangency or under-resolved fold -- which the Tangle logs and discards. There
are far fewer of these than before, but blasting still produces some. This script
reproduces the figure from ``scripts/henon_blast_period_3.py`` and scatters every
such artifact on top so they can be inspected.

The artifacts are captured with a logging handler on ``tanglepack.numerics.Tangle``
(the discard is logged as ``"Discarding impossible unstable x unstable segment
pair (id, id) ..."``). Each logged segment-id pair is looked up in the Tangle's
live ``_seg_lookup`` and the crossing point of the two segments is computed here --
so nothing is added to the core library.

Run:  MPLBACKEND=Agg PYTHONPATH=src python tests/visualizations/viz_unstable_artifacts.py
"""

from __future__ import annotations

import logging

import numpy as np
import matplotlib.pyplot as plt

from tanglepack import TangleSession

OUT = "/tmp/viz_unstable_artifacts.png"

# ── capture the discarded unstable x unstable pairs ─────────────────────────
_captured: list[tuple] = []


class _DiscardCapture(logging.Handler):
    def emit(self, record):
        msg = record.msg if isinstance(record.msg, str) else ""
        if "Discarding impossible" in msg and record.args:
            stab1, stab2, seg_pair = record.args
            if stab1 == "unstable" and stab2 == "unstable":
                _captured.append(tuple(seg_pair))


_tangle_logger = logging.getLogger("tanglepack.numerics.Tangle")
_tangle_logger.addHandler(_DiscardCapture())
_tangle_logger.setLevel(logging.WARNING)


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


def _segment_crossing(p1, p2, p3, p4):
    """Intersection point of segments p1-p2 and p3-p4, or their average if parallel."""
    d1 = p2 - p1
    d2 = p4 - p3
    den = d1[0] * d2[1] - d1[1] * d2[0]
    if abs(den) < 1e-15:
        return 0.5 * (p1 + p3)
    t = ((p3[0] - p1[0]) * d2[1] - (p3[1] - p1[1]) * d2[0]) / den
    return p1 + t * d1


def _artifact_points(tangle):
    """Resolve each captured segment-id pair to a crossing coordinate."""
    pts = []
    for seg_pair in _captured:
        ids = list(seg_pair)
        if not all(i in tangle._seg_lookup for i in ids):
            continue
        s0, s1 = tangle._seg_lookup[ids[0]], tangle._seg_lookup[ids[1]]
        pts.append(
            _segment_crossing(
                np.asarray(s0.p0.get_point()),
                np.asarray(s0.p0_seg1.get_point()),
                np.asarray(s1.p0.get_point()),
                np.asarray(s1.p0_seg1.get_point()),
            )
        )
    if not pts:
        return np.empty((0, 2))
    pts = np.array(pts)
    # de-duplicate coincident artifacts (same crossing re-detected each iterate)
    rounded = np.round(pts, 4)
    _, keep = np.unique(rounded, axis=0, return_index=True)
    return pts[np.sort(keep)]


def main():
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
    result = session.blast_zone(zone, num_iterations=4, fixed_point=fp1)

    artifacts = _artifact_points(session.workbench.Tangle)
    print(f"\nunstable x unstable artifacts logged: {len(_captured)}")
    print(f"distinct artifact locations: {len(artifacts)}")

    # --- the regular henon_blast_period_3 figure, plus the artifacts ---------
    plt.figure(figsize=(9, 9))
    session.plot_resonance_zones(alpha=0.3)
    for fp in (fp3, fp1):
        session.workbench.plot_tangle(fp, "stable", color="r", linewidth=0.8)
        session.workbench.plot_tangle(fp, "unstable", color="b", linewidth=0.6)
    for bridge in result.all_interior_bridges():
        pts = np.asarray(bridge.get_point_array())
        if len(pts):
            plt.plot(pts[:, 0], pts[:, 1], color="orange", linewidth=1.0)

    if len(artifacts):
        plt.scatter(
            artifacts[:, 0],
            artifacts[:, 1],
            marker="x",
            s=14,
            color="magenta",
            alpha=0.55,
            linewidths=0.8,
            zorder=20,
            label=f"unstable×unstable artifacts ({len(artifacts)})",
        )
        plt.legend(loc="upper right", fontsize=9)

    plt.xlim([-8, 8])
    plt.ylim([-8, 8])
    plt.gca().set_aspect("equal")
    plt.title("Blasted period-3 zone with unstable×unstable artifacts overlaid")
    plt.tight_layout()
    plt.savefig(OUT, dpi=140)
    print(f"saved figure -> {OUT}")

    plt.show()


if __name__ == "__main__":
    main()
