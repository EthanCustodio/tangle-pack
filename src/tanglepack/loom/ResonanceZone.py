from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Iterable, Optional, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from ..numerics.Intersection import Intersection, ManifoldKey

if TYPE_CHECKING:
    from ..numerics.BaseManifold import BaseManifold
    from ..numerics.BranchPoint import BranchPoint
    from ..numerics.FixedPoint import FixedPoint
    from ..numerics.Point import Point
    from ..numerics.TangleWorkbench import TangleWorkbench

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

"""
Dev Notes — Resonance zones from a chosen pip

A resonance zone is bounded by arcs of the stable and unstable manifolds that run
between the periodic point(s) and a primary intersection point (a pip). The loom
builds it by trimming the stable manifold at the pip and recomputing intersections
so the shortened stable arc takes effect.

Multi-period anchors: a period-k orbit's zone is bounded not by a single pip but by
the strong pip AND its k-1 forward iterates — one cut point on each of the k stable
branches (registry.iterate_orbit gathers them). Each branch is trimmed at its own
cut point. The closed boundary then alternates around the orbit: from a pip, in
along its stable branch to a periodic point, out along that point's unstable branch
to the next pip, and so on back to the start. For k = 1 this reduces to the single
loop z → (unstable) → pip → (stable) → z.

Why trim to ``segment.p0_seg1``: the stable manifold's canonical distance grows from
0 at the periodic point outward; the segment carrying a pip has its near endpoint
below the pip's stable cdist and its far endpoint (``p0_seg1``) above it, so setting
``manifold.tail = p0_seg1`` keeps the curve up to just past the pip — the convention
TangleWorkbench.trim_stable_manifolds uses, here generalized to chosen pips.

Registry-id caveat: recompute rebuilds the registry, so ``boundary_intersection_id``
is only valid against the registry that existed when the zone was defined. The cut
crossings are preserved as Intersection objects (``cut_intersections``), which carry
the canonical distances and branch keys and are stable across recomputes.

Inversion (k_value == 2*period) is not yet validated here (no inversion example in
the codebase) — the boundary traversal matches branches by orbit index, which is
exact only for the non-inversion case. Mirrors the StrongPip inversion caveat.
"""


@dataclass
class ResonanceZone:
    """
    A resonance zone defined by trimming stable manifold(s) at a pip (+ its iterates).

    Attributes:
        fixed_point: The periodic point whose zone this is.
        stable_branch_key: Stable branch of the primary pip (the strong pip).
        boundary_intersection: The primary pip (strong pip) as an Intersection.
        boundary_intersection_id: The primary pip's registry id at definition time.
            Valid only against the pre-recompute registry (see Dev Notes).
        cut_intersections: The pip and its iterates — one per stable branch — that
            bound the zone. Length 1 for period 1, k for a period-k orbit.
        intersection_ids: Registry ids present after the trim + recompute.
        previous_tails: Each trimmed stable branch's tail before trimming, keyed by
            manifold key, so the trim can be undone with :meth:`restore`.
    """

    fixed_point: "FixedPoint"
    stable_branch_key: ManifoldKey
    boundary_intersection: Intersection
    boundary_intersection_id: int
    cut_intersections: list[Intersection] = field(default_factory=list)
    intersection_ids: list[int] = field(default_factory=list)
    previous_tails: dict[ManifoldKey, "Point | BranchPoint"] = field(
        default_factory=dict
    )

    @property
    def branch_index(self) -> int:
        """The branch index of the primary stable branch (0, or 0/1 under inversion)."""
        return self.stable_branch_key[3]

    @property
    def key(self) -> tuple["FixedPoint", int]:
        """Storage key ``(fixed_point, branch_index)`` — one zone per branch.

        A non-inversion periodic point has a single branch (one zone); an inversion
        point (``k_value == 2 * period``) has two branches and so two zones.
        """
        return (self.fixed_point, self.branch_index)

    def boundary_polygon(
        self, workbench: "TangleWorkbench", *, close: bool = True
    ) -> NDArray[np.float64]:
        """
        The closed boundary of the resonance zone, as an (N, 2) array, for shading.

        Walks each cut point's unstable arc (periodic point → pip) and stable arc
        (periodic point → pip), then stitches them into one closed ring that
        alternates around the orbit: stable arc reversed (pip → periodic point),
        then the unstable arc out of that periodic point to the next pip, repeating
        until it returns to the start. For k = 1 this is the single z → pip → z loop.

        Args:
            workbench: The workbench whose manifolds back this zone.
            close: If True (default), repeat the first vertex at the end.

        Returns:
            (N, 2) array of (x, y) boundary vertices (empty if no cut points).
        """
        cuts = self.cut_intersections or [self.boundary_intersection]
        tol = workbench.intersection_registry.cdist_tol

        def arc(key: ManifoldKey, stability: str, cutoff: float, pip) -> list[NDArray]:
            """Manifold points from the root out to the pip (oriented root → pip)."""
            manifold = workbench.manifolds[key]
            nodes = manifold.get_point_array(return_nodes=True)
            pts = [n.get_point() for n in nodes if n.get_cdist(stability) <= cutoff + tol]
            pts.append(np.asarray(pip, dtype=float))
            return pts

        unstable_arcs, stable_arcs, u_orbit, b_orbit = [], [], [], []
        for ix in cuts:
            pip = ix.coords
            unstable_arcs.append(arc(ix.manifold_a_key, "unstable", ix.unstable_cdist, pip))
            stable_arcs.append(arc(ix.manifold_b_key, "stable", ix.stable_cdist, pip))
            u_orbit.append(ix.manifold_a_key[2])
            b_orbit.append(ix.manifold_b_key[2])

        # Index the cut whose unstable branch starts at each periodic point, so we can
        # hop pip -> periodic point (stable) -> next pip (unstable) around the orbit.
        by_u_orbit = {orbit: i for i, orbit in enumerate(u_orbit)}

        ring: list[NDArray] = []
        cur, start, n = 0, 0, len(cuts)
        for _ in range(n):
            ring.extend(stable_arcs[cur][::-1])  # pip -> periodic point
            nxt = by_u_orbit.get(b_orbit[cur])
            if nxt is None:
                # Iterates incomplete: fall back to closing this single pip's loop.
                ring.extend(unstable_arcs[cur][::-1])
                break
            ring.extend(unstable_arcs[nxt])  # periodic point -> next pip
            cur = nxt
            if cur == start:
                break

        if close and ring:
            ring = ring + [ring[0]]
        return np.vstack(ring) if ring else np.empty((0, 2))

    def restore(self, workbench: "TangleWorkbench", *, recompute: bool = True) -> None:
        """
        Undo the trim: put every trimmed stable branch's tail back where it was.

        Args:
            workbench: The workbench the zone was defined on.
            recompute: If True (default), recompute intersections over all fixed
                points so the registry reflects the restored manifolds.

        Raises:
            ValueError: If this zone captured no tails to restore.
        """
        if not self.previous_tails:
            raise ValueError("Resonance zone has no captured tails to restore.")
        for key, tail in self.previous_tails.items():
            workbench.manifolds[key].tail = tail
        if recompute:
            workbench.compute_intersections(list(workbench.fixed_points))


def _trim_stable_at(workbench: "TangleWorkbench", ix: Intersection) -> "BaseManifold":
    """Trim the stable branch of intersection ``ix`` to just past it. See
    :func:`trim_stable_at_intersection`."""
    key = ix.manifold_b_key
    if key is None:
        raise ValueError(
            "Intersection has no stable side (manifold_b_key); cannot trim a stable "
            "manifold at it."
        )
    manifold = workbench.manifolds[key]
    tangle = workbench.Tangle
    seg_ids = tangle._manifold_segs.get(manifold)
    if not seg_ids:
        raise ValueError(
            f"Stable manifold {key} has no indexed segments; run "
            "compute_intersections() before trimming."
        )

    target = ix.stable_cdist
    tol = workbench.intersection_registry.cdist_tol
    far = [(seg, tangle._seg_lookup[sid].p0_seg1.get_cdist("stable"))
           for sid, seg in ((s, tangle._seg_lookup[s]) for s in seg_ids)]
    beyond = [(seg, cdist) for seg, cdist in far if cdist >= target - tol]
    seg, _ = (min(beyond, key=lambda sc: sc[1]) if beyond
              else max(far, key=lambda sc: sc[1]))
    manifold.tail = seg.p0_seg1
    return manifold


def trim_stable_at_intersection(
    workbench: "TangleWorkbench", intersection_id: int
) -> "BaseManifold":
    """
    Truncate the stable manifold carrying an intersection to just past it.

    Resolves the stable branch from the intersection's ``manifold_b_key``, finds the
    segment whose far endpoint is the first one at or beyond the intersection's stable
    canonical distance, and sets ``manifold.tail`` to that endpoint — the chosen-pip
    generalization of :meth:`TangleWorkbench.trim_stable_manifolds`.

    The trim only moves the tail pointer; call ``workbench.compute_intersections(...)``
    afterward (as :func:`define_resonance_zone` does) to rebuild crossings on the
    shortened segment — ``Tangle._segments_of`` walks ``root → tail`` and respects it.

    Args:
        workbench: The workbench holding the manifold and registry.
        intersection_id: Registry id of the chosen crossing (e.g. a strong pip).

    Returns:
        The trimmed stable BaseManifold.

    Raises:
        ValueError: If the intersection has no stable side or its manifold is not indexed.
    """
    return _trim_stable_at(workbench, workbench.intersection_registry[intersection_id])


def define_resonance_zone(
    workbench: "TangleWorkbench",
    intersection_id: int,
    fixed_points: Optional[Iterable["FixedPoint"]] = None,
    *,
    recompute: bool = True,
) -> ResonanceZone:
    """
    Define a resonance zone by trimming the stable manifold(s) at a chosen pip.

    For a period-1 anchor this trims the single stable manifold at the pip. For a
    period-k anchor it gathers the pip and its k-1 iterates (one per stable branch,
    via ``registry.iterate_orbit``) and trims each branch at its own cut point, so the
    whole period-k zone is bounded consistently. By default it then recomputes
    intersections so the shortened stable arcs take effect; the recompute spans every
    fixed point (or ``fixed_points`` if given) so a co-indexed nested tangle survives.

    Args:
        workbench: The workbench to operate on.
        intersection_id: Registry id of the boundary pip (e.g. ``trellis.strong_pip``).
        fixed_points: Fixed points to re-index on recompute. Defaults to all.
        recompute: If True (default), recompute intersections after trimming.

    Returns:
        A :class:`ResonanceZone` recording the cut points, the trimmed branches, the
        recomputed crossings, and the pre-trim tails (for :meth:`ResonanceZone.restore`).
    """
    registry = workbench.intersection_registry
    primary = registry[intersection_id]
    key = primary.manifold_b_key
    if key is None:
        raise ValueError(
            f"Intersection {intersection_id} has no stable side; cannot define a "
            "resonance zone from it."
        )
    fixed_point = key[0]

    # Strong pip + its iterates: one cut point per stable branch (k for a period-k orbit).
    max_len = getattr(fixed_point, "k_value", None)
    cut_ids = registry.iterate_orbit(intersection_id, max_len=max_len)
    cut_intersections = [registry[c] for c in cut_ids]

    previous_tails: dict[ManifoldKey, "Point | BranchPoint"] = {}
    for ix in cut_intersections:
        bkey = ix.manifold_b_key
        if bkey is None or bkey in previous_tails:
            continue
        previous_tails[bkey] = workbench.manifolds[bkey].tail
        _trim_stable_at(workbench, ix)

    intersection_ids: list[int] = []
    if recompute:
        fps = (
            list(workbench.fixed_points)
            if fixed_points is None
            else list(fixed_points)
        )
        workbench.compute_intersections(fps, preserve_ids=True)
        intersection_ids = workbench.intersection_registry.all_ids()

    return ResonanceZone(
        fixed_point=fixed_point,
        stable_branch_key=key,
        boundary_intersection=primary,
        boundary_intersection_id=intersection_id,
        cut_intersections=cut_intersections,
        intersection_ids=intersection_ids,
        previous_tails=previous_tails,
    )
