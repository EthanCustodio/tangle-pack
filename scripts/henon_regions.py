"""Plot the closed regions of a simple k=10 Hénon tangle with their representative points.

Builds the standard single-saddle tangle, computes crossings and bridges, builds the
planar arrangement, then fills every closed region with its own colour and marks the
representative point (the "centre" the region layer uses for containment and image
checks) with a labelled dot. Open faces and containing faces are not drawn.

Run:  env/bin/python scripts/henon_regions.py   (writes figures/henon_regions.png)
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

import tanglepack
from tanglepack.examples.henon import (
    henon_jacobian,
    henon_map,
    henon_map_inverse,
    saddle_guesses,
)

logging.basicConfig(level=logging.WARNING)

FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"
OUT_PATH = FIGURES_DIR / "henon_regions.png"

K, B = 10, 1
UNSTABLE_ITERATIONS = 10


def build_session() -> tuple[tanglepack.TangleSession, object]:
    """Grow the k=10 tangle far enough to enclose a handful of regions."""
    session = tanglepack.TangleSession(
        henon_map(K, B), henon_map_inverse(K, B), henon_jacobian(K, B)
    )
    fp = session.construct_fixed_point(saddle_guesses(K, B)["saddle"])
    session.orient_eigenvectors(
        fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
    )
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=UNSTABLE_ITERATIONS)
    session.grow_until_turnaround(fp, "stable")
    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    return session, fp


def main() -> None:
    session, fp = build_session()
    arrangement = session.arrangement()
    print(arrangement.summary())

    fig, ax = plt.subplots(figsize=(8, 8))
    session.workbench.plot_tangle(fp, "unstable", color="tab:blue", lw=0.8)
    session.workbench.plot_tangle(fp, "stable", color="tab:red", lw=0.8)

    colors = plt.get_cmap("tab20").resampled(max(len(arrangement.regions), 1))
    for index, region in enumerate(arrangement.regions):
        boundary = region.boundary_points
        ax.fill(
            boundary[:, 0],
            boundary[:, 1],
            color=colors(index),
            alpha=0.35,
            lw=0,
            label=f"region {region.corners} |A|={abs(region.area):.3g}",
        )
        centre = region.representative_point
        if centre is None:
            continue
        ax.plot(centre[0], centre[1], "o", color="k", ms=5, zorder=20)
        ax.annotate(
            str(index),
            centre,
            textcoords="offset points",
            xytext=(4, 4),
            fontsize=8,
            zorder=21,
        )
        print(
            f"region {index}: corners={region.corners} |area|={abs(region.area):.4g} "
            f"centre=({centre[0]:.4f}, {centre[1]:.4f}) "
            f"contains(centre)={region.contains(centre)}"
        )

    ax.plot(*fp.coordinates[0], "k*", ms=10, zorder=22)
    # Frame the view on the (trimmed) stable manifold: the primary lobes run
    # thousands of units out along the escaping unstable arm, so a bounding box
    # of the regions themselves would shrink the tangle to a dot.
    stable = np.vstack(
        [
            manifold.get_point_array()
            for (fixed_point, stability, _o, _b), manifold in session.manifolds.items()
            if fixed_point is fp and stability == "stable"
        ]
    )
    lo, hi = stable.min(axis=0), stable.max(axis=0)
    pad = 0.1 * float(np.linalg.norm(hi - lo))
    ax.set_xlim(lo[0] - pad, hi[0] + pad)
    ax.set_ylim(lo[1] - pad, hi[1] + pad)
    ax.set_aspect("equal")
    ax.set_title(
        f"k={K} Hénon tangle: {len(arrangement.regions)} closed regions and their centres"
    )
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=150)
    plt.close(fig)
    print(f"saved {OUT_PATH}")


if __name__ == "__main__":
    main()
