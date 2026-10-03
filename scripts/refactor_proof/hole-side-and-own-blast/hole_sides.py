"""
Proof of fix A: a propagated hole's side is read on its own image sub-arc.

Builds period 3 alone at 6 blasts (where 2315204 raised in
``check_holes_share_bridge_side``) and records every propagated hole as it
is punched: its containing bridge, the backward-carried reference point and
the unstable-cdist span of its image sub-arc. For each hole the side is read
twice, the old way (nearest vertex of the WHOLE bridge) and the new way (the
sub-arc only, ``_bridge_side_of(..., span)``). ``hole_sides.png`` draws one
panel per hole whose two readings differ: the whole bridge in grey, the
image sub-arc in black, the carried point, and the vertex each reading
measures against -- the old one on another fold.

Run from the branch worktree::

    MPLBACKEND=Agg ../tangle-pack/env/bin/python scripts/refactor_proof/hole-side-and-own-blast/hole_sides.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

BRANCH = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(BRANCH / "src"))
OUT = BRANCH / "figures" / "refactor_proof" / "hole-side-and-own-blast"

import tanglepack  # noqa: E402
from tanglepack.examples.henon_cases import build_period3  # noqa: E402
from tanglepack.topology import StablePartition as sp  # noqa: E402

assert Path(tanglepack.__file__).is_relative_to(BRANCH), tanglepack.__file__

BLASTS = 6


def record_punches() -> list[dict]:
    """Every propagated punch of a period-3 build: the hole and what its side was read from."""
    punches = []
    punch = sp._punch_in_bridge

    def recording(trellis, bridge, **kwargs):
        hole = punch(trellis, bridge, **kwargs)
        if hole is not None and kwargs.get("carried") is not None:
            punches.append({"trellis": trellis, "bridge": bridge, "hole": hole, **kwargs})
        return hole

    sp._punch_in_bridge = recording
    try:
        build = build_period3(blasts=BLASTS)
    finally:
        sp._punch_in_bridge = punch
    trellis = build.session.trellis(build.fixed_points[0])
    sp.check_holes_share_bridge_side(trellis.holes, orientation_preserving=trellis.orientation_preserving)
    return punches


def nearest(poly: np.ndarray, point: np.ndarray) -> int:
    """The index `_arc_side_of` measures at: the nearest vertex, kept off the ends."""
    i = int(np.argmin(np.linalg.norm(poly - point, axis=1)))
    return min(max(i, 1), len(poly) - 2)


def draw_hole(ax, punch: dict, reference_side: str) -> None:
    """One hole: whole bridge, image sub-arc, carried point and both readings."""
    trellis, bridge, carried, span = punch["trellis"], punch["bridge"], punch["carried"], punch["span"]
    poly = sp._oriented_bridge_polyline(trellis, bridge)
    cdists = sorted(node.cdist for node in bridge.get_point_array(return_nodes=True))
    lo = max(sp.bisect_left(cdists, span[0]) - 1, 0)
    hi = sp.bisect_right(cdists, span[1]) + 1
    sub = poly[lo:hi]
    old_side = sp._bridge_side_of(trellis, bridge, carried)
    new_side = sp._bridge_side_of(trellis, bridge, carried, span)
    old_i, new_i = nearest(poly, carried), lo + nearest(sub, carried)

    ax.plot(*poly.T, color="0.75", linewidth=1, label="containing bridge (whole)")
    ax.plot(*sub.T, color="black", linewidth=2.2, label="image sub-arc (span)")
    ax.annotate("", xy=poly[-1], xytext=poly[-2], arrowprops={"arrowstyle": "->", "color": "0.6"})
    ax.plot(*carried, "*", color="tab:purple", markersize=13, label="carried reference point")
    for i, colour, text in ((old_i, "tab:red", f"old: {old_side}"), (new_i, "tab:green", f"new: {new_side}")):
        ax.plot(*poly[i], "o", color=colour, markersize=9, markerfacecolor="none", markeredgewidth=2)
        ax.plot([carried[0], poly[i][0]], [carried[1], poly[i][1]], "--", color=colour, linewidth=1)
        ax.annotate(f"{text}\n(vertex at u-cdist {cdists[i]:.3g})", poly[i], color=colour, fontsize=8,
                    xytext=(6, 6), textcoords="offset points")
    focus = np.vstack([sub, carried, poly[old_i]])
    (x0, y0), (x1, y1) = focus.min(axis=0), focus.max(axis=0)
    pad = 0.15 * max(x1 - x0, y1 - y0)
    ax.set_xlim(x0 - pad, x1 + pad)
    ax.set_ylim(y0 - pad, y1 + pad)
    ax.set_aspect("equal")
    hole = punch["hole"]
    ax.set_title(
        f"origin {hole.origin}, iterate {hole.iterate}\n"
        f"sub-arc u-cdist {span[0]:.4g}..{span[1]:.4g}; iterate 0 is {reference_side}",
        fontsize=10,
    )


def main() -> None:
    punches = record_punches()
    changed = [
        p for p in punches
        if sp._bridge_side_of(p["trellis"], p["bridge"], p["carried"])
        != sp._bridge_side_of(p["trellis"], p["bridge"], p["carried"], p["span"])
    ]
    print(f"{len(punches)} propagated holes, {len(changed)} read differently on the sub-arc")
    trellis = punches[0]["trellis"]
    reference = {hole.origin: hole.bridge_side for hole in trellis.holes if hole.iterate == 0}
    fig, axes = plt.subplots(1, len(changed), figsize=(6.5 * len(changed), 6.5), squeeze=False)
    for ax, punch in zip(axes[0], changed):
        draw_hole(ax, punch, reference.get(punch["hole"].origin))
    axes[0][0].legend(fontsize=8, loc="best")
    fig.suptitle(
        f"Period 3 alone, {BLASTS} blasts: the propagated holes whose side the old reading got wrong "
        f"({len(changed)} of {len(punches)})",
        fontsize=13,
    )
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "hole_sides.png", dpi=150)
    print("wrote", OUT / "hole_sides.png")


if __name__ == "__main__":
    main()
