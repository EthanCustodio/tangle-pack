from __future__ import annotations

import logging
from typing import Iterable, Optional, TYPE_CHECKING

from ..numerics.TangleWorkbench import TangleWorkbench
from ..numerics.DynamicalSystem import MapFunc, JacFunc
from ..topology.Trellis import Trellis
from .ResonanceZone import ResonanceZone, define_resonance_zone

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

        Args:
            fixed_points: A single FixedPoint, an iterable of them, or None for all
                of the workbench's fixed points. The cache key is by identity.
            rebuild: Force a rebuild even if a cached Trellis exists.

        Returns:
            The Trellis for the requested fixed point(s).
        """
        cache_key = self._cache_key(fixed_points)
        if rebuild or cache_key not in self._trellises:
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
            self.invalidate_trellises()
            recomputed = self.workbench.intersection_registry.all_ids()
            for rz in self.resonance_zones.values():
                rz.intersection_ids = recomputed
        return self.resonance_zones

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
            poly = rz.boundary_polygon(self.workbench)
            if poly.size == 0:
                continue
            color = palette[i % len(palette)]
            handles.append(
                target.fill(
                    poly[:, 0], poly[:, 1], color=color, alpha=alpha, **fill_kwargs
                )
            )
        return handles
