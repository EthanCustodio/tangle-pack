from __future__ import annotations

import logging
from typing import Iterable, Optional, TYPE_CHECKING

from ..numerics.TangleWorkbench import TangleWorkbench
from ..numerics.DynamicalSystem import MapFunc, JacFunc
from ..topology.Trellis import Trellis, _is_single_fixed_point
from .ResonanceZone import ResonanceZone, define_resonance_zone
from .Blast import BlastResult, blast_zone

if TYPE_CHECKING:
    from ..numerics.FixedPoint import FixedPoint

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class TangleSession:
    """
    User-friendly facade that weaves the numerical and topological layers together.

    A session owns one :class:`TangleWorkbench` (the numerical engine) and a cache of
    :class:`Trellis` objects (the topological view), and hosts the cross-layer
    ("loom") algorithms such as resonance-zone construction. The workbench is exposed
    directly as ``self.workbench`` and, for convenience, any attribute not found on
    the session is delegated to it — so a session is a drop-in superset of a
    workbench: ``session.construct_fixed_point(...)``, ``session.grow_n_times(...)``,
    ``session.compute_intersections(...)`` all work, while the session adds
    :meth:`trellis` and :meth:`resonance_zone` on top.

    Attributes:
        workbench: The underlying TangleWorkbench.

    Example:
        >>> session = TangleSession(my_map, my_map_inv, my_jac)
        >>> fp = session.construct_fixed_point([4, -4])          # delegated to workbench
        >>> session.initialize_both_manifolds(fp)
        >>> session.grow_n_times(fp, "unstable", num_iterations=7)
        >>> session.grow_until_turnaround(fp, "stable")
        >>> session.compute_intersections(fp)
        >>> T = session.trellis(fp)                               # built + cached
        >>> T.classify_strong_pips()
        >>> rz = session.resonance_zone(T.strong_pip)             # trim + recompute
        >>> T = session.trellis(fp, rebuild=True)                 # registry changed
    """

    def __init__(
        self,
        dynamical_map: MapFunc,
        dynamical_map_inverse: MapFunc,
        jacobian_function: Optional[JacFunc] = None,
    ):
        """
        Build a session around a fresh workbench for the given map.

        Args:
            dynamical_map: The forward map.
            dynamical_map_inverse: Its inverse.
            jacobian_function: Optional analytic Jacobian.
        """
        self.workbench = TangleWorkbench(
            dynamical_map, dynamical_map_inverse, jacobian_function
        )
        self._trellises: dict = {}
        # One resonance zone per (fixed_point, branch_index): a non-inversion point
        # has a single branch (one zone); an inversion point has two. Insertion order
        # is preserved so plotting/shading is deterministic.
        self.resonance_zones: dict[tuple["FixedPoint", int], ResonanceZone] = {}

    # ── attribute delegation ─────────────────────────────────────────────────

    def __getattr__(self, name: str):
        """Delegate unknown attributes to the workbench (drop-in superset)."""
        # __getattr__ runs only when normal lookup fails. Guard against recursion
        # during construction (before self.workbench exists).
        workbench = self.__dict__.get("workbench")
        if workbench is not None and hasattr(workbench, name):
            return getattr(workbench, name)
        raise AttributeError(
            f"{type(self).__name__!r} object (and its workbench) has no attribute {name!r}"
        )

    # ── topological view ─────────────────────────────────────────────────────

    @property
    def fixed_points(self) -> list["FixedPoint"]:
        """The workbench's fixed points."""
        return self.workbench.fixed_points

    def trellis(
        self,
        fixed_points: Optional["FixedPoint | Iterable[FixedPoint]"] = None,
        *,
        rebuild: bool = False,
    ) -> Trellis:
        """
        Build (and cache) the Trellis for one or more fixed points.

        Trellises are cached per fixed-point key so repeated calls are cheap. Because
        a Trellis is a snapshot, rebuild it (``rebuild=True``) after anything that
        changes the registry — growing manifolds, recomputing intersections, or
        defining a resonance zone. :meth:`resonance_zone` clears the cache for you.

        A cached Trellis whose registry is no longer the workbench's current one
        (every :meth:`compute_intersections` builds a fresh registry — e.g. inside
        a resonance-zone recompute) is transparently rebuilt, so a session never
        hands back a snapshot bound to a stale registry. The candidate/strong-pip
        slots of the rebuilt Trellis start empty; re-run classification (or use the
        session-level :meth:`classify_strong_pips`, which does so for you).

        Args:
            fixed_points: A single FixedPoint, an iterable of them, or None for all
                of the workbench's fixed points. The cache key is by identity.
            rebuild: Force a rebuild even if a cached Trellis exists.

        Returns:
            The Trellis for the requested fixed point(s).
        """
        cache_key = self._cache_key(fixed_points)
        cached = self._trellises.get(cache_key)
        registry = self.workbench.intersection_registry
        stale = cached is not None and (
            cached.registry is not registry
            # Blasting (iterate_bridge) adds crossings to the SAME registry
            # object; the snapshot's per-branch orderings miss them.
            or getattr(cached, "_built_registry_size", None) != len(registry)
        )
        if rebuild or cached is None or stale:
            self._trellises[cache_key] = Trellis.from_workbench(
                self.workbench, fixed_points
            )
        return self._trellises[cache_key]

    def invalidate_trellises(self) -> None:
        """Drop all cached Trellises (they are snapshots of a now-stale registry)."""
        self._trellises.clear()

    @staticmethod
    def _cache_key(fixed_points):
        """A hashable cache key for a trellis selection."""
        if fixed_points is None:
            return None
        if isinstance(fixed_points, (list, tuple, set, frozenset)):
            return frozenset(id(fp) for fp in fixed_points)
        return id(fixed_points)

    def _resolve_fixed_points(
        self, fixed_point: "Optional[FixedPoint | Iterable[FixedPoint]]"
    ) -> list["FixedPoint"]:
        """Normalize a fixed-point selector to a list (None → every fixed point)."""
        if fixed_point is None:
            return list(self.workbench.fixed_points)
        if _is_single_fixed_point(fixed_point):
            return [fixed_point]
        return list(fixed_point)

    # ── strong-pip convenience (per fixed point) ─────────────────────────────
    #
    # Each fixed point owns its own (single-fixed-point) Trellis, so a nested
    # session has one trellis per tangle. These helpers fan a single call out
    # across the fixed points and route every access through :meth:`trellis`, so
    # they always operate on a fresh (auto-rebuilt) snapshot rather than a stale
    # ``T1``/``T3`` variable kept across a resonance-zone recompute.

    def classify_strong_pips(
        self,
        fixed_point: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        **kwargs,
    ):
        """
        Classify strong-pip candidates for one, several, or every fixed point.

        For each selected fixed point this builds (or reuses) its per-fixed-point
        Trellis and runs :meth:`Trellis.classify_strong_pips`, so the inner and
        outer tangles of a nested session are each classified against their own
        crossings. Routing through :meth:`trellis` means a trellis left stale by a
        resonance-zone recompute is rebuilt and reclassified rather than silently
        reusing dropped candidates.

        Args:
            fixed_point: A single FixedPoint (→ returns that fixed point's candidate
                list), an iterable of them, or None for every fixed point on the
                workbench (→ returns a ``{fixed_point: candidates}`` dict).
            **kwargs: Forwarded to :meth:`Trellis.classify_strong_pips` (``tol``,
                ``collision_rtol``, ``choose_default``).

        Returns:
            A candidate-id list for a single fixed point, or a dict mapping each
            fixed point to its candidate-id list.
        """
        results = {
            fp: self.trellis(fp).classify_strong_pips(**kwargs)
            for fp in self._resolve_fixed_points(fixed_point)
        }
        if _is_single_fixed_point(fixed_point):
            return results[fixed_point]
        return results

    def set_strong_pip(self, fixed_point: "FixedPoint", intersection_id: int) -> int:
        """Choose ``intersection_id`` as the strong pip for ``fixed_point``'s trellis."""
        return self.trellis(fixed_point).set_strong_pip(intersection_id)

    def strong_pip(self, fixed_point: "FixedPoint") -> Optional[int]:
        """The chosen strong-pip id for ``fixed_point``'s trellis (or None)."""
        return self.trellis(fixed_point).strong_pip

    def strong_pip_candidates(self, fixed_point: "FixedPoint") -> list[int]:
        """The strong-pip candidate ids for ``fixed_point``'s trellis."""
        return self.trellis(fixed_point).strong_pip_candidates

    def plot_strong_pip_candidates(
        self,
        fixed_point: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        *,
        classify: bool = True,
        ax=None,
        **scatter_kwargs,
    ) -> list:
        """
        Scatter the strong-pip candidates of one, several, or every fixed point.

        Loops over the selected fixed points' trellises (default: all), so the
        candidates of both the inner and outer tangle of a nested session are
        drawn in a single call — no per-trellis bookkeeping. See
        :meth:`Trellis.plot_strong_pip_candidates`.

        Args:
            fixed_point: Fixed point selector; None (default) plots every one.
            classify: If a trellis has no candidates yet, classify it first
                (default True). Pass False to plot only already-classified trellises.
            ax: Optional matplotlib Axes (defaults to the current axes).
            **scatter_kwargs: Forwarded to :meth:`Trellis.plot_strong_pip_candidates`.

        Returns:
            List of the matplotlib handles drawn (one per fixed point that had
            candidates to plot).
        """
        handles = []
        for fp in self._resolve_fixed_points(fixed_point):
            trellis = self.trellis(fp)
            if classify and not trellis.strong_pip_candidates:
                trellis.classify_strong_pips()
            handle = trellis.plot_strong_pip_candidates(ax=ax, **scatter_kwargs)
            if handle is not None:
                handles.append(handle)
        return handles

    def plot_strong_pip(
        self,
        fixed_point: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        *,
        classify: bool = True,
        ax=None,
        **scatter_kwargs,
    ) -> list:
        """
        Scatter the chosen strong-pip cut points of one, several, or every fixed point.

        Companion to :meth:`plot_strong_pip_candidates`. For each selected fixed
        point this draws the strong pip (and, for a period-k anchor, its iterates
        bounding the resonance zone) via :meth:`Trellis.plot_strong_pip`.

        Args:
            fixed_point: Fixed point selector; None (default) plots every one.
            classify: If a trellis has no strong pip chosen yet, classify it first
                (which also picks the default strong pip). Pass False to plot only
                trellises that already have one.
            ax: Optional matplotlib Axes (defaults to the current axes).
            **scatter_kwargs: Forwarded to :meth:`Trellis.plot_strong_pip`.

        Returns:
            List of the matplotlib handles drawn (one per fixed point with a strong
            pip).
        """
        handles = []
        for fp in self._resolve_fixed_points(fixed_point):
            trellis = self.trellis(fp)
            if classify and trellis.strong_pip is None:
                trellis.classify_strong_pips()
            handle = trellis.plot_strong_pip(ax=ax, **scatter_kwargs)
            if handle is not None:
                handles.append(handle)
        return handles

    # ── pseudoneighbor convenience (per fixed point) ─────────────────────────

    def compute_pseudoneighbors(
        self,
        fixed_point: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        **kwargs,
    ):
        """
        Find reference pseudoneighbors for one, several, or every fixed point.

        For each selected fixed point this builds (or reuses) its
        per-fixed-point Trellis and runs
        :meth:`Trellis.compute_pseudoneighbors`, so each tangle of a nested
        session is checked against its own crossings.

        Args:
            fixed_point: A single FixedPoint (→ returns that fixed point's
                reference-pair list), an iterable of them, or None for every
                fixed point (→ returns a ``{fixed_point: references}`` dict).
            **kwargs: Forwarded to :meth:`Trellis.compute_pseudoneighbors`
                (``extend``, ``collision_rtol``, ``match_rtol``, ``tol``).

        Returns:
            A reference-pair list for a single fixed point, or a dict mapping
            each fixed point to its list.
        """
        results = {
            fp: self.trellis(fp).compute_pseudoneighbors(**kwargs)
            for fp in self._resolve_fixed_points(fixed_point)
        }
        if _is_single_fixed_point(fixed_point):
            return results[fixed_point]
        return results

    def plot_pseudoneighbors(
        self,
        fixed_point: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        *,
        compute: bool = True,
        ax=None,
        **scatter_kwargs,
    ) -> list:
        """
        Scatter the pseudoneighbors of one, several, or every fixed point.

        Args:
            fixed_point: Fixed point selector; None (default) plots every one.
            compute: If a trellis has no pseudoneighbors yet, compute them
                first (default True). Pass False to plot only already-computed
                trellises.
            ax: Optional matplotlib Axes (defaults to the current axes).
            **scatter_kwargs: Forwarded to :meth:`Trellis.plot_pseudoneighbors`
                (including ``include_trajectories``).

        Returns:
            List of the matplotlib handles drawn (one per fixed point with
            pairs to plot).
        """
        handles = []
        for fp in self._resolve_fixed_points(fixed_point):
            trellis = self.trellis(fp)
            if compute and not trellis.pseudoneighbors:
                trellis.compute_pseudoneighbors()
            handle = trellis.plot_pseudoneighbors(ax=ax, **scatter_kwargs)
            if handle is not None:
                handles.append(handle)
        return handles

    def plot_holes(
        self,
        fixed_point: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        *,
        ax=None,
        **scatter_kwargs,
    ) -> list:
        """
        Scatter the punched holes of one, several, or every fixed point.

        Holes must have been punched already (:meth:`Trellis.punch_holes`);
        trellises without holes are skipped.

        Args:
            fixed_point: Fixed point selector; None (default) plots every one.
            ax: Optional matplotlib Axes (defaults to the current axes).
            **scatter_kwargs: Forwarded to :meth:`Trellis.plot_holes`.

        Returns:
            List of the matplotlib handles drawn.
        """
        handles = []
        for fp in self._resolve_fixed_points(fixed_point):
            handle = self.trellis(fp).plot_holes(ax=ax, **scatter_kwargs)
            if handle is not None:
                handles.extend(handle)
        return handles

    # ── cross-layer ("loom") algorithms ──────────────────────────────────────

    def resonance_zone(
        self,
        intersection_id: int,
        fixed_points: Optional[Iterable["FixedPoint"]] = None,
        *,
        recompute: bool = True,
    ) -> ResonanceZone:
        """
        Define and store one resonance zone by trimming the stable manifold at a pip.

        The zone is stored in :attr:`resonance_zones` keyed by ``(fixed_point,
        branch_index)``, so each periodic point keeps its own zone (two under
        inversion) and a later call for a different fixed point does not overwrite it.
        Trims are cumulative (each sets a manifold's tail), so building zones for
        several fixed points nests correctly — the inner zone ends up a subset of the
        outer one. When recomputing, the Trellis cache is invalidated (the registry
        was rebuilt).

        Args:
            intersection_id: Registry id of the boundary pip (e.g. ``trellis.strong_pip``).
            fixed_points: Fixed points to re-index on recompute. Defaults to all.
            recompute: If True (default), recompute intersections after trimming. Pass
                False when batching several zones (see :meth:`add_resonance_zones`),
                so all pip ids stay valid against one registry until a single final
                recompute.

        Returns:
            The stored :class:`ResonanceZone`.
        """
        rz = define_resonance_zone(
            self.workbench, intersection_id, fixed_points, recompute=recompute
        )
        self.resonance_zones[rz.key] = rz
        if recompute:
            self.invalidate_trellises()
        return rz

    def add_resonance_zones(
        self,
        intersection_ids: Iterable[int],
        fixed_points: Optional[Iterable["FixedPoint"]] = None,
    ) -> dict[tuple["FixedPoint", int], ResonanceZone]:
        """
        Define resonance zones for several pips at once (one recompute at the end).

        All pips are trimmed first with the registry untouched, so every id in
        ``intersection_ids`` stays valid against the registry that exists now; a
        single recompute then re-indexes the whole co-indexed tangle with every trim
        in effect. Use this to build, e.g., the outer period-1 zone and the inner
        period-3 zone together.

        Args:
            intersection_ids: Boundary-pip registry ids, one per zone.
            fixed_points: Fixed points to re-index on the final recompute. Defaults to all.

        Returns:
            The :attr:`resonance_zones` mapping.
        """
        ids = list(intersection_ids)
        for iid in ids:
            # Trim only (no recompute) so the remaining ids keep their meaning.
            self.resonance_zone(iid, recompute=False)
        if ids:
            fps = (
                list(self.workbench.fixed_points)
                if fixed_points is None
                else list(fixed_points)
            )
            self.workbench.compute_intersections(fps, preserve_ids=True)
            # Stable manifolds are now trimmed; recut bridges against the new crossings.
            self.workbench.rebuild_bridges()
            self.invalidate_trellises()
            recomputed = self.workbench.intersection_registry.all_ids()
            for rz in self.resonance_zones.values():
                rz.intersection_ids = recomputed
        return self.resonance_zones

    # ── bridge ↔ resonance-zone classification ───────────────────────────────

    @staticmethod
    def _bridge_test_point(bridge) -> Optional["NDArray"]:
        """The midpoint node of a bridge, used as its representative point.

        Returns the geometric middle node so a bridge that forms a zone's own unstable
        boundary arc lands exactly on that zone's boundary (and, with boundary-inclusive
        containment, is attributed to it). ``None`` for an empty bridge.
        """
        pts = bridge.get_point_array()
        if pts is None or len(pts) == 0:
            return None
        return pts[len(pts) // 2]

    def classify_bridge(self, bridge) -> Optional[ResonanceZone]:
        """
        Determine which resonance zone a single bridge lies in.

        Tests the bridge's representative midpoint against every stored zone's frozen
        boundary (boundary included) and returns the innermost containing zone — the
        one of smallest area when zones nest. A bridge on a zone's own unstable
        boundary arc is counted as inside that (innermost) zone.

        Args:
            bridge: A :class:`~tanglepack.numerics.Bridge.Bridge`.

        Returns:
            The innermost :class:`ResonanceZone` containing the bridge, or ``None`` if
            it lies outside every zone (or no zones are defined).
        """
        point = self._bridge_test_point(bridge)
        if point is None:
            return None
        containing = [
            rz for rz in self.resonance_zones.values() if rz.contains_point(point)
        ]
        if not containing:
            return None
        return min(containing, key=lambda rz: rz.area)

    def classify_bridges(
        self, bridges: Optional[Iterable] = None
    ) -> dict:
        """
        Classify many bridges by resonance zone.

        Args:
            bridges: Bridges to classify. Defaults to every bridge on the workbench.

        Returns:
            A dict mapping each bridge to its innermost containing :class:`ResonanceZone`
            (or ``None``); see :meth:`classify_bridge`.
        """
        bridges = (
            list(self.workbench.bridges) if bridges is None else list(bridges)
        )
        return {bridge: self.classify_bridge(bridge) for bridge in bridges}

    # ── blasting ─────────────────────────────────────────────────────────────

    def blast_zone(
        self,
        zone: "ResonanceZone | tuple",
        num_iterations: int,
        *,
        fixed_point: "Optional[FixedPoint | list[FixedPoint]]" = None,
        strict: bool = False,
        min_separation: Optional[float] = None,
    ) -> BlastResult:
        """
        Repeatedly iterate the bridges inside a resonance zone.

        Each step maps every un-iterated interior bridge forward (cutting the image
        into child bridges at its new stable-manifold crossings) and keeps the children
        that land back inside the zone, recursing up to ``num_iterations`` times.
        Bridges that leave the zone are never re-iterated — iterating exterior bridges
        many times grows them exponentially. See :mod:`tanglepack.loom.Blast`.

        Args:
            zone: A :class:`ResonanceZone` or its ``(fixed_point, branch_index)`` key.
            num_iterations: Maximum number of blast steps.
            fixed_point: Restrict to bridges from this fixed point (default: all).
            strict: Re-raise a bridge's forward-map failure instead of skipping it.
            min_separation: Drop a child bridge whose interior comes within this
                distance of an already-kept sibling, to avoid precision-driven merges
                of near-coincident unstable curves (default: disabled).

        Returns:
            A :class:`~tanglepack.loom.Blast.BlastResult` genealogy of the blast.

        Note:
            Blasting registers the children's new stable-manifold crossings, so
            any cached Trellis is invalidated — take a fresh :meth:`trellis`
            (and re-establish its strong pip) before continuing topological
            work.
        """
        result = blast_zone(
            self,
            zone,
            num_iterations,
            fixed_point=fixed_point,
            strict=strict,
            min_separation=min_separation,
        )
        self.invalidate_trellises()
        return result

    def plot_resonance_zones(
        self,
        ax=None,
        colors: Optional[list] = None,
        alpha: float = 0.3,
        **fill_kwargs,
    ):
        """
        Shade every stored resonance zone as a filled region, each a distinct color.

        Fills are translucent (low ``alpha``) so the red stable / blue unstable
        manifold lines and any nested inner zone remain visible on top. Zones are
        drawn in insertion order.

        Args:
            ax: Optional matplotlib Axes. Defaults to the current axes (plt).
            colors: Optional list of fill colors cycled across zones. Defaults to a
                distinct qualitative palette.
            alpha: Fill transparency (default 0.3).
            **fill_kwargs: Forwarded to ``fill`` (e.g. ``zorder``, ``label``).

        Returns:
            List of the matplotlib polygon handles drawn (one per zone).
        """
        import matplotlib.pyplot as plt

        if not self.resonance_zones:
            logger.info("No resonance zones to plot; call resonance_zone() first.")
            return []

        palette = colors or [
            "#8c6bb1",  # purple
            "#41ab5d",  # green
            "#fd8d3c",  # orange
            "#6baed6",  # light blue
            "#df65b0",  # magenta
            "#a6761d",  # ochre
        ]
        target = ax if ax is not None else plt
        handles = []
        for i, rz in enumerate(self.resonance_zones.values()):
            poly = (
                rz.boundary_vertices
                if rz.boundary_vertices is not None
                else rz.boundary_polygon(self.workbench)
            )
            if poly.size == 0:
                continue
            color = palette[i % len(palette)]
            handles.append(
                target.fill(
                    poly[:, 0], poly[:, 1], color=color, alpha=alpha, **fill_kwargs
                )
            )
        return handles
