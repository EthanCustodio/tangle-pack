from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Iterable, Literal, Optional, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from ..numerics.Intersection import Intersection, ManifoldKey
from ..numerics.geometry import arc_polyline, point_in_polygon, signed_polygon_area
from ..topology.TopologyResults import Arc

if TYPE_CHECKING:
    from ..numerics.BaseManifold import BaseManifold
    from ..numerics.IntersectionRegistry import IntersectionRegistry
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

Registry-id caveat: a recompute rebuilds the registry, so an id held across one is
only meaningful when the recompute preserved ids. The zone therefore stores no id
at all: the cut crossings are kept as Intersection objects (``cut_intersections``),
which carry the canonical distances and branch keys and are stable across
recomputes, and the primary pip's id is resolved against the live registry on
demand by ``resolve_boundary_intersection_id``.

Inversion (k_value == 2*period) is not yet validated here (no inversion example in
the codebase) — the boundary traversal matches branches by orbit index, which is
exact only for the non-inversion case. Mirrors the StrongPip inversion caveat.
"""


#: A resonance-zone boundary piece is described with the region layer's own
#: :class:`~tanglepack.topology.TopologyResults.Arc`, so a zone boundary and a
#: region boundary are made of the same thing and measure the same way. The name
#: is kept as an alias because the zone's arcs are a slightly wider notion: one
#: runs from a periodic point all the way OUT to its pip, which may span several
#: consecutive arcs of the arrangement rather than exactly one.
BoundaryArc = Arc


@dataclass
class ResonanceZone:
    """
    A resonance zone defined by trimming stable manifold(s) at a pip (+ its iterates).

    Attributes:
        fixed_point: The periodic point whose zone this is.
        stable_branch_key: Stable branch of the primary pip (the strong pip).
        boundary_intersection: The primary pip (strong pip) as an Intersection.
        cut_intersections: The pip and its iterates — one per stable branch — that
            bound the zone. Length 1 for period 1, k for a period-k orbit. Exposed as
            :attr:`boundary_intersections` (the preferred name).
        intersection_ids: Registry ids present after the trim + recompute.
        previous_tails: Each trimmed stable branch's tail before trimming, keyed by
            manifold key, so the trim can be undone with :meth:`restore`.
        boundary_arcs: The stable/unstable :class:`~tanglepack.topology.TopologyResults.Arc`
            objects forming the closed boundary, captured by
            :meth:`capture_boundary`. Each runs from a periodic point (canonical
            distance 0, the anchor crossing) out to a cut point. Empty until
            captured.
        boundary_vertices: Frozen (N, 2) snapshot of the closed boundary polygon taken
            at definition time, used for in/out-of-zone tests so classification is
            robust to later trims/recomputes. ``None`` until captured.
    """

    fixed_point: "FixedPoint"
    stable_branch_key: ManifoldKey
    boundary_intersection: Intersection
    cut_intersections: list[Intersection] = field(default_factory=list)
    intersection_ids: list[int] = field(default_factory=list)
    previous_tails: dict[ManifoldKey, "Point | BranchPoint"] = field(
        default_factory=dict
    )
    boundary_arcs: list[Arc] = field(default_factory=list)
    boundary_vertices: Optional[NDArray[np.float64]] = None

    @property
    def boundary_intersections(self) -> list[Intersection]:
        """The pips bounding this zone — the strong pip plus its per-branch iterates.

        For a period-1 anchor this is a single pip; for a period-k orbit it is the k
        cut points, one on each stable branch. Falls back to the primary
        :attr:`boundary_intersection` if no iterates were gathered.
        """
        return self.cut_intersections or [self.boundary_intersection]

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

        Rebuilt live from the workbench's current manifolds; for a snapshot frozen at
        definition time use :attr:`boundary_vertices` (set by :meth:`capture_boundary`).

        Args:
            workbench: The workbench whose manifolds back this zone.
            close: If True (default), repeat the first vertex at the end.

        Returns:
            (N, 2) array of (x, y) boundary vertices (empty if no cut points).
        """
        return self._build_boundary(workbench, close=close)[0]

    def _build_boundary(
        self, workbench: "TangleWorkbench", *, close: bool = True
    ) -> tuple[NDArray[np.float64], list[Arc]]:
        """Build the closed boundary polygon and the list of arcs that compose it.

        See :meth:`boundary_polygon` for the stitching order. Returns the (N, 2)
        vertex array and the composing :class:`~tanglepack.topology.TopologyResults.Arc`
        objects (one unstable and one stable arc per cut point), each running from a
        periodic point's anchor crossing out to that cut point.
        """
        cuts = self.cut_intersections or [self.boundary_intersection]
        registry = workbench.intersection_registry
        tol = registry.cdist_tol

        arcs: list[Arc] = []

        def arc(
            key: ManifoldKey,
            stability: Literal["unstable", "stable"],
            cutoff: float,
            pip: "NDArray[np.float64] | tuple[float, float]",
            pip_id: Optional[int],
        ) -> list[NDArray[np.float64]]:
            """Manifold points from the root out to the pip (oriented root → pip).

            Also records the piece as an :class:`Arc` when both of its endpoint
            crossings can be named: the anchor at the root of ``key`` and the pip
            at the far end. A missing id means the arc cannot be described in the
            shared vocabulary, so it is left out of ``arcs`` while its vertices
            still go into the ring.

            Args:
                key: The branch this arc runs along.
                stability: Which canonical distance ``cutoff`` is measured in.
                cutoff: The pip's canonical distance; the arc spans ``[0, cutoff]``.
                pip: The pip's coordinates, appended as the arc's far end.
                pip_id: The pip's registry id in the CURRENT registry, or None.

            Returns:
                The arc's vertices, root first.
            """
            manifold = workbench.manifolds[key]
            pts = list(arc_polyline(manifold, 0.0, cutoff, stability=stability, tol=tol))
            pts.append(np.asarray(pip, dtype=float))
            anchor_id = _anchor_id(registry, key)
            if anchor_id is None or pip_id is None:
                logger.debug(
                    "Not recording a boundary Arc on %s: anchor id %s, pip id %s",
                    key,
                    anchor_id,
                    pip_id,
                )
            else:
                arcs.append(
                    Arc(
                        kind=stability,
                        lo_id=anchor_id,
                        hi_id=pip_id,
                        branch_key=key,
                        bridge_id=None,
                    )
                )
            return pts

        unstable_arcs, stable_arcs, u_orbit, b_orbit = [], [], [], []
        for ix in cuts:
            pip = ix.coords
            pip_id = _live_id(registry, ix)
            unstable_arcs.append(
                arc(ix.manifold_a_key, "unstable", ix.unstable_cdist, pip, pip_id)
            )
            stable_arcs.append(
                arc(ix.manifold_b_key, "stable", ix.stable_cdist, pip, pip_id)
            )
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
        vertices = np.vstack(ring) if ring else np.empty((0, 2))
        return vertices, arcs

    def capture_boundary(self, workbench: "TangleWorkbench") -> None:
        """
        Freeze this zone's boundary geometry from the workbench's current manifolds.

        Stores the closed polygon as :attr:`boundary_vertices` and the composing
        :attr:`boundary_arcs`. Call this once the zone's own stable branch(es) have
        been trimmed (as :func:`define_resonance_zone` does), so the snapshot reflects
        the trimmed boundary and stays valid for in/out tests even after later zones
        trim or recompute the registry.

        Args:
            workbench: The workbench whose manifolds back this zone.
        """
        self.boundary_vertices, self.boundary_arcs = self._build_boundary(
            workbench, close=True
        )

    def resolve_boundary_intersection_id(
        self, workbench: "TangleWorkbench"
    ) -> Optional[int]:
        """
        This zone's primary pip, as an id in the workbench's CURRENT registry.

        The id is resolved on demand rather than stored, because a recompute that
        does not preserve ids would silently turn a stored one into a different
        crossing. The pip itself is kept as an :class:`Intersection`, whose
        canonical distances and branch keys survive any rebuild, and those are what
        the registry is asked for.

        Args:
            workbench: The workbench whose registry to resolve against.

        Returns:
            The pip's current registry id, or None if the crossing is no longer in
            the registry (it can be trimmed away by a later zone).
        """
        return _live_id(workbench.intersection_registry, self.boundary_intersection)

    @property
    def area(self) -> float:
        """Enclosed area of the captured boundary via the shoelace formula.

        Returns 0.0 if the boundary has not been captured or is degenerate. The value
        is the absolute area, so it is independent of vertex winding direction; the
        nested-zone classifier in :class:`TangleSession` orders zones by it.
        """
        return abs(signed_polygon_area(self.boundary_vertices))

    def contains_point(
        self,
        point: "NDArray[np.float64] | tuple[float, float]",
        *,
        tol: float = 1e-9,
    ) -> bool:
        """
        Whether ``point`` lies inside this zone, boundary included.

        A point on the boundary (e.g. the midpoint node of a bridge that forms this
        zone's own unstable boundary arc) counts as inside — such bridges belong to the
        innermost zone they bound. That inclusion is the shared policy of
        :func:`~tanglepack.numerics.geometry.point_in_polygon`, which this delegates
        to, so a zone and a region answer the same way.

        Args:
            point: An (x, y) coordinate.
            tol: Distance below which the point is taken to lie on the boundary.

        Returns:
            True if the point is inside or on the boundary; False otherwise (including
            when the boundary has not been captured).
        """
        verts = self.boundary_vertices
        if verts is None:
            return False
        return point_in_polygon(point, verts, tol=tol)

    def restore(self, workbench: "TangleWorkbench", *, recompute: bool = True) -> None:
        """
        Undo the trim: put every trimmed stable branch's tail back where it was.

        The mirror image of :func:`define_resonance_zone`, which trims, recomputes and
        recuts the bridges: restoring the tails lengthens the stable arcs again, so the
        crossing set changes and the existing bridges — cut against the *trimmed*
        manifolds — are stale. The recompute is therefore followed by
        :meth:`TangleWorkbench.rebuild_bridges`, which returns the bridge set to
        exactly what it was before the trim. Like the trim, the recompute preserves
        registry ids, so ids held across the round trip keep their meaning.

        Args:
            workbench: The workbench the zone was defined on.
            recompute: If True (default), recompute intersections over all fixed
                points and rebuild the bridges so both reflect the restored manifolds.

        Note:
            ``recompute=False`` deliberately leaves the registry *and* the bridges
            stale: it is the batched path, where several zones are restored in turn
            and the caller runs one ``compute_intersections`` + ``rebuild_bridges``
            at the end (mirroring :meth:`TangleSession.add_resonance_zones`).

        Raises:
            ValueError: If this zone captured no tails to restore.
        """
        if not self.previous_tails:
            raise ValueError("Resonance zone has no captured tails to restore.")
        for key, tail in self.previous_tails.items():
            workbench.manifolds[key].tail = tail
        if recompute:
            # preserve_ids mirrors define_resonance_zone: ids held across the trim
            # (a strong pip, a BridgeId) must still name the same crossings after the
            # restore, and rebuild_bridges carries its per-bridge metadata by id.
            workbench.compute_intersections(
                list(workbench.fixed_points), preserve_ids=True
            )
            workbench.rebuild_bridges()


def _anchor_id(
    registry: "IntersectionRegistry", key: ManifoldKey
) -> Optional[int]:
    """The registry id of the anchor crossing at the root of one branch.

    A branch is anchored at a periodic point, which is registered as a crossing of
    that point's unstable and stable branches at canonical distance ``(0, 0)`` (see
    ``TangleWorkbench._register_anchors``). The probe names the pair of branches
    exactly, so the registry's own collision test resolves it without any
    coordinate comparison.

    Args:
        registry: The registry to resolve against.
        key: The branch whose anchor is wanted; its orbit index and branch index
            pick out the partner branch of the other stability.

    Returns:
        The anchor's registry id, or None if that branch pair has no anchor.
    """
    fixed_point, stability, orbit_index, branch_index = key
    unstable = (
        key if stability == "unstable"
        else (fixed_point, "unstable", orbit_index, branch_index)
    )
    stable = (
        key if stability == "stable"
        else (fixed_point, "stable", orbit_index, branch_index)
    )
    return registry.find(
        Intersection.synthetic(
            coords=(0.0, 0.0),
            unstable_cdist=0.0,
            stable_cdist=0.0,
            manifold_a_key=unstable,
            manifold_b_key=stable,
        )
    )


def _live_id(
    registry: "IntersectionRegistry", intersection: Intersection
) -> Optional[int]:
    """The id of ``intersection`` in the CURRENT registry, matched not assumed.

    Args:
        registry: The registry to resolve against.
        intersection: The crossing to find, by canonical distances and branch keys.

    Returns:
        Its current registry id, or None if it is no longer registered.
    """
    return registry.find(intersection)


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
    intersections so the shortened stable arcs take effect and recuts the bridges
    against the new crossings; the recompute spans every fixed point (or
    ``fixed_points`` if given) so a co-indexed nested tangle survives.

    Args:
        workbench: The workbench to operate on.
        intersection_id: Registry id of the boundary pip (e.g. ``trellis.strong_pip``).
        fixed_points: Fixed points to re-index on recompute. Defaults to all.
        recompute: If True (default), recompute intersections after trimming.

    Returns:
        A :class:`ResonanceZone` recording the cut points, the trimmed branches, the
        recomputed crossings, the pre-trim tails (for :meth:`ResonanceZone.restore`),
        and a frozen boundary snapshot (for bridge classification).
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
        # The trimmed stable arc yields a different crossing set, so the old bridges
        # are stale — recut them against the new (shorter) stable manifold.
        workbench.rebuild_bridges()
        intersection_ids = workbench.intersection_registry.all_ids()

    zone = ResonanceZone(
        fixed_point=fixed_point,
        stable_branch_key=key,
        boundary_intersection=primary,
        cut_intersections=cut_intersections,
        intersection_ids=intersection_ids,
        previous_tails=previous_tails,
    )
    # Freeze the boundary now that this zone's stable branch(es) are trimmed, so the
    # snapshot used for bridge classification is robust to later trims/recomputes.
    zone.capture_boundary(workbench)
    return zone
