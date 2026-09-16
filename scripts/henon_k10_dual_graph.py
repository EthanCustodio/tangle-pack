"""
The whole pipeline on the single k=10 Hénon saddle, start to finish.

One fixed point, one tangle, every stage in order: map -> saddle -> manifolds
-> growth -> crossings -> bridges -> iterate table -> strong pip ->
pseudoneighbors -> holes -> stable partition -> bridge classes -> dual graph.
Logs the dual graph's summary and saves one 1x3 figure: the tangle with its
strong pip, pseudoneighbors and holes; the stable partition; and the tangle
with the dual graph overlaid.

This is also the demonstration of the pip-driven fill: ``PIP`` picks which
strong-pip candidate drives it, and the filled arc nodes move with it.

Run:  env/bin/python scripts/henon_k10_dual_graph.py

Dev Notes:
    The growth recipe is the ``k10_session`` fixture of ``tests/conftest.py``
    (saddle at ``[4, -4]``, stable grown to its turnaround) with the unstable
    manifold grown ``UNSTABLE_STEPS`` times: 9 in the fixture, 10 here, where
    the partition refines from 3 to 5 elements per side.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from tanglepack import TangleSession
from tanglepack.examples import HENON_K10, henon_map, henon_map_inverse
from tanglepack.topology.plotting import dual_graph_legend_handles

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"
OUT_PATH = FIGURES_DIR / "henon_k10_dual_graph.png"

SADDLE_GUESS = [4, -4]
UNSTABLE_STEPS = 10

# Which strong pip to cut at. "default" is the candidate with the smallest
# unstable cdist (the primary crossing q0, outermost on the stable branch);
# "image" is its forward image f(q0), one step in along the stable branch --
# the point at about (3.18, 4.05) on the right arm. An int is a raw registry
# id, which is NOT stable between builds (f(q0) is id 3 at 9 steps and id 5 at
# 10), so prefer the named choices.
#
# ``set_strong_pip`` records the choice on the trellis, and the dual-graph
# fill follows it directly: the arcs between the chosen pip and its k-th
# iterate are filled whether or not the manifold is trimmed there. The
# stable manifold itself stays trimmed at its outermost crossing, though,
# so the pseudoneighbors, holes and partition still see the default cut
# unless ``TRIM_AT_PIP`` trims the stable manifold at the chosen pip (a
# resonance zone) and recomputes. The recompute rebuilds the trellis, so
# the pip is classified and set again afterwards.
PIP: str | int = "image"
TRIM_AT_PIP = True


def _tangle_view(session: TangleSession, fp, pad: float = 0.08) -> tuple:
    """
    Axis limits framing the stable manifold, padded by ``pad`` of its diagonal.

    The unstable manifold runs off to x ~ 2500, so the stable manifold sets
    the view. Call this BEFORE trimming at the chosen pip: the trimmed
    manifold is far shorter and would crop the figure to a corner of the
    tangle.

    Args:
        session: The session holding the manifolds.
        fp: The fixed point whose stable manifold frames the view.
        pad: Padding as a fraction of the bounding box diagonal.

    Returns:
        ``((xmin, xmax), (ymin, ymax))``.
    """
    points = np.vstack(
        [
            manifold.get_point_array()
            for (fixed_point, stability, _o, _b), manifold in session.manifolds.items()
            if fixed_point is fp and stability == "stable"
        ]
    )
    mins, maxs = points.min(axis=0), points.max(axis=0)
    margin = pad * float(np.linalg.norm(maxs - mins))
    return (mins[0] - margin, maxs[0] + margin), (mins[1] - margin, maxs[1] + margin)


def _log_candidates(session: TangleSession, fp) -> None:
    """Log every strong-pip candidate with its point and canonical distances."""
    trellis = session.trellis(fp)
    for index, cid in enumerate(trellis.strong_pip_candidates):
        ix = trellis.intersection(cid)
        logger.info(
            "candidate %d: id %d at %s  stable cdist %.3f  unstable cdist %.1f%s",
            index,
            cid,
            np.round(ix.get_point(), 3),
            ix.stable_cdist,
            ix.unstable_cdist,
            "  (default)" if cid == trellis.strong_pip else "",
        )


def _choose_pip(session: TangleSession, fp) -> int:
    """
    Resolve ``PIP`` to a registry id among this trellis's strong-pip candidates.

    Args:
        session: The session, after ``classify_strong_pips``.
        fp: The fixed point whose trellis holds the candidates.

    Returns:
        The registry id to hand to ``set_strong_pip``.

    Raises:
        ValueError: If ``PIP`` names no candidate, or the default pip has no
            registered forward image.
    """
    trellis = session.trellis(fp)
    default = trellis.strong_pip
    if PIP == "default":
        return default
    if PIP == "image":
        image = trellis.iterate(default, 1)
        if image is None:
            raise ValueError(
                f"the default strong pip {default} has no registered image"
            )
        return image
    if PIP in trellis.strong_pip_candidates:
        return PIP
    raise ValueError(
        f"PIP={PIP!r} is not a strong-pip candidate; candidates are "
        f"{trellis.strong_pip_candidates}"
    )


def main() -> None:
    """Run the pipeline, log the dual graph, save the figure."""
    # 1. The map and the saddle.
    session = TangleSession(henon_map(*HENON_K10), henon_map_inverse(*HENON_K10))
    fp = session.construct_fixed_point(SADDLE_GUESS)
    session.orient_eigenvectors(
        fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
    )

    # 2. Manifolds: initialise both, grow the unstable a fixed number of steps
    #    and the stable until it turns around.
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=UNSTABLE_STEPS)
    session.grow_until_turnaround(fp, "stable")

    # 3. Crossings, trimming, bridges, iterate table.
    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    session.infer_iterate_table()
    logger.info("%s", session.trellis(fp).summary())

    # 4. The topological algorithms: strong pip, pseudoneighbors, holes, partition.
    session.classify_strong_pips()
    _log_candidates(session, fp)
    pip = session.set_strong_pip(fp, _choose_pip(session, fp))
    xlim, ylim = _tangle_view(session, fp)
    if TRIM_AT_PIP:
        # Trim the stable manifold at the pip and recompute (ids are preserved).
        # The recompute advances the workbench generation, which drops the
        # cached trellis and its strong-pip choice, so classify and set again.
        session.resonance_zone(pip)
        session.classify_strong_pips()
        session.set_strong_pip(fp, pip)
    logger.info("%s", session.trellis(fp).summary())

    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()
    pip = session.strong_pip(fp)

    logger.info(
        "strong pip: %s at %s",
        pip,
        np.round(session.trellis(fp).intersection(pip).get_point(), 3),
    )
    logger.info("%s", session.describe_stable_partitions(fp))

    # 5. Bridge classes and the dual graph.
    dual_graph = session.dual_graph()
    logger.info("%s", session.describe_bridge_classes())
    logger.info("%s", dual_graph.summary())

    title = f"Hénon k={HENON_K10[0]}"
    print(title)
    print("=" * len(title))
    print(session.describe_bridge_classes())
    print(dual_graph.summary())

    # 6. The figure.
    fig, (ax_tangle, ax_partition, ax_dual) = plt.subplots(1, 3, figsize=(21, 7))

    for ax in (ax_tangle, ax_dual):
        plt.sca(ax)
        session.plot_tangle(fp, "unstable", color="tab:blue", linewidth=0.6)
        session.plot_tangle(fp, "stable", color="tab:red", linewidth=0.6)

    # Both tangle panels share the pre-trim view (see _tangle_view) rather
    # than the dual-graph plotter's own clip to the surviving arcs.
    session.plot_dual_graph(dual_graph, ax=ax_dual, clip_to_arcs=False)
    ax_dual.legend(
        handles=dual_graph_legend_handles(), loc="best", fontsize=7, framealpha=0.8
    )
    ax_dual.set_title(f"{title} -- dual graph")

    session.plot_pseudoneighbors(fp, ax=ax_tangle)
    session.plot_holes(fp, ax=ax_tangle)
    session.workbench.plot_intersections(fp, show_ids=True, ax=ax_tangle)
    # The pip is also a pseudoneighbor and a crossing, so a filled dot the
    # same size vanishes under those markers: draw it last, as a green ring.
    session.plot_strong_pip(
        fp, ax=ax_tangle, s=120, facecolors="none", linewidths=1.5, zorder=20
    )
    for ax in (ax_tangle, ax_dual):
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    ax_tangle.set_title(f"{title} -- tangle, strong pip, pseudoneighbors, holes")

    # The partition plot lives on a canonical-distance axis of its own.
    session.plot_stable_partition(fp, ax=ax_partition)
    ax_partition.set_title(f"{title} -- stable partition")

    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=150)
    plt.close(fig)
    logger.info("saved %s", OUT_PATH)


if __name__ == "__main__":
    main()
