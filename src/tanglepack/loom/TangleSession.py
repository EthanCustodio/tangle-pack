"""
The user-facing facade over the numerical and topological layers.

:class:`TangleSession` owns one
:class:`~..numerics.TangleWorkbench.TangleWorkbench` and a cache of
:class:`~..topology.Trellis.Trellis` objects, delegates any unknown attribute
to the workbench, and hosts the cross-layer ("loom") algorithms -- resonance
zones and blasting -- that need both halves at once.
"""

from __future__ import annotations

import logging
import warnings
from typing import Callable, Iterable, Optional, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from ..numerics.TangleWorkbench import TangleWorkbench
from ..numerics.DynamicalSystem import MapFunc, JacFunc
from ..numerics.geometry import polyline_midpoint
from ..topology import plotting
from ..topology.Arrangement import Arrangement
from ..topology.BridgeClass import BridgeClass, bridge_classes as _bridge_classes
from ..topology.DualGraph import DualGraph
from ..topology.TopologyResults import endpoint_index
from ..topology.Trellis import Trellis, _is_single_fixed_point
from .ResonanceZone import ResonanceZone, define_resonance_zone
from .Blast import BlastResult, blast_zone

if TYPE_CHECKING:
    from matplotlib.axes import Axes

    from ..numerics.Bridge import Bridge, BridgeId
    from ..numerics.FixedPoint import FixedPoint
    from ..topology.TopologyResults import Endpoint, Side, StablePartitionResult

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

    The predicate-driven growth drivers
    (:meth:`~tanglepack.numerics.TangleWorkbench.TangleWorkbench.grow_until` and
    its ``grow_until_iterates_closed`` / ``grow_until_faces_closed`` wrappers) are
    delegated to the workbench like any other growth call, and need no help from
    the session: every round they run moves the workbench generation, so the next
    :meth:`trellis` or :meth:`arrangement` call rebuilds itself.

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
        # (workbench.generation, partition signature, classes) per
        # trellis-selection cache key; see bridge_classes().
        self._bridge_classes: dict = {}
        # (workbench.generation, partition signature, strong pips, graph) per
        # trellis-selection cache key; see dual_graph().
        self._dual_graphs: dict = {}
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

        Trellises are cached per fixed-point key so repeated calls are cheap, and
        a stale one is rebuilt automatically; ``rebuild=True`` forces a rebuild
        even when the cached snapshot is still valid.

        Staleness is one comparison: a Trellis records the workbench
        :attr:`~tanglepack.numerics.TangleWorkbench.TangleWorkbench.generation`
        it was built at, and anything that changes the tangle — growth, a
        recompute (which swaps in a fresh registry), a blast adding crossings to
        the same registry, a re-cut, a trim, a resonance zone or its restore —
        advances that counter, so the cached snapshot is dropped and rebuilt on
        the next call. No caller has to invalidate by hand. The
        candidate/strong-pip slots of the rebuilt Trellis start empty; re-run
        classification (or use the session-level :meth:`classify_strong_pips`,
        which does so for you).

        Args:
            fixed_points: A single FixedPoint, an iterable of them, or None for all
                of the workbench's fixed points. The cache key is by identity.
            rebuild: Force a rebuild even if a cached Trellis exists.

        Returns:
            The Trellis for the requested fixed point(s).
        """
        cache_key = self._cache_key(fixed_points)
        cached = self._trellises.get(cache_key)
        stale = cached is not None and self._is_stale(cached)
        if rebuild or cached is None or stale:
            self._trellises[cache_key] = Trellis.from_workbench(
                self.workbench, fixed_points
            )
        return self._trellises[cache_key]

    def arrangement(
        self,
        fixed_points: Optional["FixedPoint | Iterable[FixedPoint]"] = None,
        *,
        rebuild: bool = False,
    ) -> Arrangement:
        """
        Build (and cache) the planar :class:`Arrangement` of a trellis.

        The default selection is EVERY fixed point, which is what the region layer
        wants: a heteroclinic crossing belongs to two tangles at once, and a
        single-fixed-point arrangement would cut the other side's arcs off and turn
        real faces into open ones.

        Cached ON the trellis, and so invalidated by exactly the same comparison:
        a Trellis is a snapshot of one workbench
        :attr:`~tanglepack.numerics.TangleWorkbench.TangleWorkbench.generation`, so
        growth, a recompute, a re-cut, a trim, a blast or a resonance zone drops the
        trellis and its arrangement together. There is no second cache to go stale.

        Args:
            fixed_points: A single FixedPoint, an iterable of them, or None (the
                default) for all of them.
            rebuild: Force a rebuild even if a valid arrangement is cached.

        Returns:
            The Arrangement for that selection.
        """
        trellis = self.trellis(fixed_points, rebuild=rebuild)
        if rebuild:
            trellis._arrangement = None
        return trellis.arrangement

    def bridge_classes(
        self,
        fixed_points: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        *,
        rebuild: bool = False,
    ) -> dict[BridgeClass, list["BridgeId"]]:
        """
        Group the bridges of a trellis by the pair of partition elements they connect.

        The classes are computed on the trellis selected by ``fixed_points`` (all
        fixed points by default, mirroring :meth:`arrangement`), but the stable
        partitions used to resolve each bridge end's element are gathered from
        EVERY cached, non-stale per-fixed-point trellis (see
        :meth:`_gathered_partitions`) rather than from the selected trellis's own
        (empty) :attr:`~tanglepack.topology.Trellis.Trellis.stable_partitions`. A
        heteroclinic bridge is a piece of one fixed point's unstable manifold, so
        it lives only in that fixed point's trellis, while the elements its ends
        sit against belong to the partitions of whichever (possibly different)
        fixed point owns the stable branch it crosses — the same split
        :meth:`partition_element_for` exists to bridge for a single lookup, done
        here for every bridge of the selection at once. See
        :func:`tanglepack.topology.BridgeClass.bridge_classes` for the class
        semantics themselves.

        Cached by :attr:`~tanglepack.numerics.TangleWorkbench.TangleWorkbench.generation`,
        exactly like :meth:`trellis`: any growth, recompute, re-cut, trim, blast
        or resonance-zone change advances the generation, so the next call
        recomputes rather than returning a stale grouping. The generation alone
        cannot see re-partitioning, though — punching holes and partitioning are
        Trellis-level mutations that do not touch the workbench — so the cache
        is also keyed by a structural :meth:`_partition_signature` of the
        gathered partitions; see the Note below.

        Args:
            fixed_points: A single FixedPoint, an iterable of them, or None (the
                default) for all of them — selects which trellis's bridges are
                classed.
            rebuild: Force a recompute even if a cached grouping exists for the
                current generation and partition signature.

        Returns:
            A dict from :class:`~tanglepack.topology.BridgeClass.BridgeClass` to
            the sorted list of that class's
            :data:`~tanglepack.numerics.Bridge.BridgeId` s (see
            :func:`tanglepack.topology.BridgeClass.bridge_classes`).

        Raises:
            ValueError: Propagated from
                :func:`tanglepack.topology.BridgeClass.bridge_classes` when a
                bridge end needs a partition that no gathered
                :class:`~tanglepack.topology.TopologyResults.StablePartitionResult`
                covers.

        Note:
            Call :meth:`partition_stable_manifold` (directly or via the session
            fan-out) on every fixed point of interest first — with nothing
            partitioned yet, every bridge end fails to resolve and this raises.
            A branch missing a partition on either side is logged as a WARNING,
            naming which side(s), before that happens.

            The cache hit test is ``(generation, signature)`` together: a
            :class:`~tanglepack.topology.Trellis.Trellis` can be re-partitioned
            (e.g. after :meth:`~tanglepack.topology.Trellis.Trellis.clear_results`
            and a fresh :meth:`partition_stable_manifold`) without the workbench
            generation moving at all, since the mutation happens on the Trellis,
            not the workbench. :meth:`_partition_signature` reads off exactly
            the interval boundaries :func:`tanglepack.topology.BridgeClass.bridge_classes`
            depends on, so a different partition of the same branch is detected
            even though the generation is unchanged.
        """
        cache_key = self._cache_key(fixed_points)
        trellis = self.trellis(fixed_points)
        partitions = self._gathered_partitions()
        signature = self._partition_signature(partitions)

        cached = self._bridge_classes.get(cache_key)
        if (
            not rebuild
            and cached is not None
            and cached[0] == self.workbench.generation
            and cached[1] == signature
        ):
            return cached[2]

        self._warn_unpartitioned_branches(trellis, partitions)
        classes = _bridge_classes(trellis, partitions)

        self._bridge_classes[cache_key] = (self.workbench.generation, signature, classes)
        return classes

    @staticmethod
    def _partition_signature(partitions: list["StablePartitionResult"]) -> tuple:
        """
        A structural signature of a gathered partition list, for cache keying.

        Two gathered-partition lists compare equal under this signature exactly
        when they name the same ``(branch_key, side)`` pairs with the same
        interval boundaries (endpoint ids and open/closed flags) in the same
        order — the fields :func:`tanglepack.topology.BridgeClass.bridge_classes`
        actually reads to resolve a bridge end to an element. Used by
        :meth:`bridge_classes` to detect a re-partition that leaves the
        workbench generation untouched (partitioning is a Trellis-level
        mutation); the dual-graph cache reuses it for the same reason, hence
        its generality (no dependence on anything bridge-class-specific).

        Args:
            partitions: The gathered
                :class:`~tanglepack.topology.TopologyResults.StablePartitionResult` s
                (see :meth:`_gathered_partitions`).

        Returns:
            A hashable, order-sensitive tuple: one
            ``(branch_key, side, (interval boundaries, ...))`` entry per result,
            each interval boundary itself a ``(lo_id, hi_id, closed_lo,
            closed_hi)`` tuple.
        """
        return tuple(
            (
                result.branch_key,
                result.side,
                tuple(
                    (interval.lo_id, interval.hi_id, interval.closed_lo, interval.closed_hi)
                    for interval in result.intervals
                ),
            )
            for result in partitions
        )

    def _gathered_partitions(self) -> list["StablePartitionResult"]:
        """
        Every stable partition of every cached, still-valid per-fixed-point trellis.

        Walks :attr:`_trellises` — the same scan :meth:`partition_element_for`
        does — skipping any snapshot :meth:`_is_stale` has flagged (its
        partitions are keyed by ids the workbench may have renumbered since),
        and concatenates the survivors' :attr:`~tanglepack.topology.Trellis.Trellis.stable_partitions`.
        Deduplicated by ``(branch_key, side)``, keeping the first result seen for
        a given key — cache iteration order is insertion order, so this favors
        whichever trellis was built (and cached) first.

        Note:
            General on purpose: the dual-graph gathering reuses this helper
            rather than re-implementing the scan.

        Returns:
            The deduplicated list of :class:`~tanglepack.topology.TopologyResults.StablePartitionResult`.
        """
        gathered: dict[tuple, "StablePartitionResult"] = {}
        for trellis in self._trellises.values():
            if self._is_stale(trellis):
                continue
            for result in trellis.stable_partitions:
                key = (result.branch_key, result.side)
                gathered.setdefault(key, result)
        return list(gathered.values())

    @staticmethod
    def _warn_unpartitioned_branches(
        trellis: Trellis, partitions: list["StablePartitionResult"]
    ) -> None:
        """Log a WARNING naming every stable branch of ``trellis`` missing a
        partition on either side, and which side(s) are missing."""
        covered = {(result.branch_key, result.side) for result in partitions}
        missing = []
        for branch in trellis.stable_branches:
            missing_sides = [
                side for side in ("left", "right") if (branch.key, side) not in covered
            ]
            if missing_sides:
                missing.append((branch.key, missing_sides))
        if missing:
            logger.warning(
                "%d stable branch(es) missing a partition on at least one side: %s",
                len(missing),
                ", ".join(f"{key[1:]} missing {sides}" for key, sides in missing),
            )

    def dual_graph(
        self,
        fixed_points: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        *,
        rebuild: bool = False,
    ) -> DualGraph:
        """
        Build (and cache) the bipartite :class:`~tanglepack.topology.DualGraph.DualGraph`
        of a trellis's arrangement.

        Mirrors :meth:`arrangement` and :meth:`bridge_classes`: the arrangement
        dualled is the selected trellis's own (all fixed points by default), and
        the stable partitions the graph's arc-node elements are resolved against
        are gathered from EVERY cached, non-stale per-fixed-point trellis (see
        :meth:`_gathered_partitions`), exactly as :meth:`bridge_classes` does.
        The strong pips passed to :class:`~tanglepack.topology.DualGraph.DualGraph`
        — they drive its fill, ``(f^k(q0), q0]`` on each pip's own branch —
        are gathered the same way, one per cached per-fixed-point trellis whose
        fixed point falls inside the selection (see :meth:`_gathered_strong_pips`).

        Cached on ``(workbench.generation, partition signature, gathered strong
        pips)``: the first two exactly as :meth:`bridge_classes` (the same
        :meth:`_partition_signature` catches a re-partition that leaves the
        workbench generation untouched), the third because a pip choice
        changes the graph's ``filled`` flags and ``fill_segments``.

        Args:
            fixed_points: A single FixedPoint, an iterable of them, or None (the
                default) for all of them — selects which trellis's arrangement
                is dualled.
            rebuild: Force a rebuild even if a cached graph exists for the
                current generation and partition signature.

        Returns:
            The :class:`~tanglepack.topology.DualGraph.DualGraph`.

        Raises:
            ValueError: Propagated from
                :class:`~tanglepack.topology.DualGraph.DualGraph` when a stable
                arc's branch has no partition on a side, or when an arrangement
                component fails the graph's own topological checks.

        Note:
            Like :meth:`bridge_classes`, this needs every fixed point of
            interest already partitioned (:meth:`partition_stable_manifold`),
            and only trellises still cached at call time contribute their
            partitions and strong pips.

            Calling :meth:`~tanglepack.topology.Trellis.Trellis.set_strong_pip`
            on a cached per-fixed-point trellis, even without re-partitioning
            or touching the workbench, changes the gathered pips and so the
            next call rebuilds the graph with the new fill. A selection with
            no classified pip builds a graph that fills nothing (the
            constructor warns).
        """
        cache_key = self._cache_key(fixed_points)
        trellis = self.trellis(fixed_points)
        partitions = self._gathered_partitions()
        signature = self._partition_signature(partitions)
        strong_pips = tuple(self._gathered_strong_pips(fixed_points))

        cached = self._dual_graphs.get(cache_key)
        if (
            not rebuild
            and cached is not None
            and cached[0] == self.workbench.generation
            and cached[1] == signature
            and cached[2] == strong_pips
        ):
            return cached[3]

        self._warn_unpartitioned_branches(trellis, partitions)
        dg = DualGraph(trellis.arrangement, partitions, strong_pips=strong_pips)

        self._dual_graphs[cache_key] = (
            self.workbench.generation,
            signature,
            strong_pips,
            dg,
        )
        return dg

    def _gathered_strong_pips(
        self, fixed_points: "Optional[FixedPoint | Iterable[FixedPoint]]"
    ) -> list[int]:
        """
        The chosen strong pip of every cached, non-stale per-fixed-point trellis
        whose fixed point falls inside ``fixed_points``.

        Walks :attr:`_trellises` like :meth:`_gathered_partitions`, but keeps
        only trellises built for exactly one fixed point (an all-fixed-points
        trellis carries no strong pip of its own) and skips a trellis with no
        strong pip chosen yet.
        """
        selected = {id(fp) for fp in self._resolve_fixed_points(fixed_points)}
        pips: list[int] = []
        for trellis in self._trellises.values():
            if self._is_stale(trellis):
                continue
            if len(trellis.fixed_points) != 1:
                continue
            if id(trellis.fixed_points[0]) not in selected:
                continue
            if trellis.strong_pip is not None:
                pips.append(trellis.strong_pip)
        return pips

    def plot_dual_graph(
        self,
        dual_graph: Optional[DualGraph] = None,
        ax: Optional["Axes"] = None,
        **kwargs,
    ) -> "Axes":
        """
        Draw a dual graph, defaulting to :meth:`dual_graph`'s own (all fixed
        points, cached) result.

        Thin delegate to
        :func:`~tanglepack.topology.plotting.plot_dual_graph`; see its
        docstring for the drawing itself.

        Args:
            dual_graph: The graph to draw. Defaults to ``self.dual_graph()``.
            ax: Optional matplotlib Axes to draw on. Defaults to the current
                axes (plt).
            **kwargs: Forwarded to
                :func:`~tanglepack.topology.plotting.plot_dual_graph` (e.g.
                ``show_labels``, ``clip_to_arcs``).

        Returns:
            The Axes drawn on.
        """
        graph = dual_graph if dual_graph is not None else self.dual_graph()
        return plotting.plot_dual_graph(graph, ax=ax, **kwargs)

    def invalidate_trellises(self) -> None:
        """
        No-op alias kept for one release.

        Cached Trellises are invalidated by the workbench generation counter
        (see :meth:`trellis`), so nothing has to be dropped by hand any more.
        Calls are harmless; new code should simply take a fresh
        :meth:`trellis`.

        Warns:
            DeprecationWarning: Always -- the behaviour changed from "drop every
                cached Trellis" to "do nothing", which an out-of-repo caller has
                to be told about while the alias still exists.
        """
        warnings.warn(
            "TangleSession.invalidate_trellises() is a no-op and will be removed: "
            "the Trellis cache is keyed by TangleWorkbench.generation, so a stale "
            "snapshot is rebuilt by the next trellis() call.",
            DeprecationWarning,
            stacklevel=2,
        )

    def _is_stale(self, trellis: Trellis) -> bool:
        """Whether a cached Trellis snapshot predates the workbench's current state."""
        return trellis._built_generation != self.workbench.generation

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

    def _fanout_call(
        self,
        name: str,
        fixed_point: "Optional[FixedPoint | Iterable[FixedPoint]]",
        **kwargs,
    ):
        """
        Run one Trellis method across a fixed-point selection.

        Every non-plot fan-out shares this shape: build (or reuse) each selected
        fixed point's trellis via :meth:`trellis` — so a snapshot left stale by a
        resonance-zone recompute is rebuilt rather than silently reused — call
        ``name`` on it, and collect the results.

        Args:
            name: The :class:`Trellis` method to call.
            fixed_point: A single FixedPoint, an iterable of them, or None for
                every fixed point on the workbench.
            **kwargs: Forwarded to the trellis method.

        Returns:
            The single fixed point's own result when one FixedPoint was passed,
            otherwise a ``{fixed_point: result}`` dict.
        """
        results = {
            fp: getattr(self.trellis(fp), name)(**kwargs)
            for fp in self._resolve_fixed_points(fixed_point)
        }
        if _is_single_fixed_point(fixed_point):
            return results[fixed_point]
        return results

    def _fanout_plot(
        self,
        name: str,
        fixed_point: "Optional[FixedPoint | Iterable[FixedPoint]]",
        ax,
        *,
        prepare: Optional[Callable[[Trellis], None]] = None,
        flatten: bool = False,
        **kwargs,
    ) -> list:
        """
        Draw one Trellis plotter across a fixed-point selection.

        Args:
            name: The :class:`Trellis` ``plot_*`` method to call.
            fixed_point: A single FixedPoint, an iterable of them, or None for
                every fixed point on the workbench.
            ax: Optional matplotlib Axes, forwarded to every plotter.
            prepare: Optional callable run on each trellis before drawing, for
                the plotters that compute their input on demand.
            flatten: True when the plotter returns a LIST of handles (extend the
                result) rather than a single handle (append it).
            **kwargs: Forwarded to the trellis plotter.

        Returns:
            The matplotlib handles drawn, skipping every trellis whose plotter
            returned None (nothing to draw).
        """
        handles: list = []
        for fp in self._resolve_fixed_points(fixed_point):
            trellis = self.trellis(fp)
            if prepare is not None:
                prepare(trellis)
            handle = getattr(trellis, name)(ax=ax, **kwargs)
            if handle is None:
                continue
            if flatten:
                handles.extend(handle)
            else:
                handles.append(handle)
        return handles

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
        return self._fanout_call("classify_strong_pips", fixed_point, **kwargs)

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

        def prepare(trellis: Trellis) -> None:
            if classify and not trellis.strong_pip_candidates:
                trellis.classify_strong_pips()

        return self._fanout_plot(
            "plot_strong_pip_candidates",
            fixed_point,
            ax,
            prepare=prepare,
            **scatter_kwargs,
        )

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

        def prepare(trellis: Trellis) -> None:
            if classify and trellis.strong_pip is None:
                trellis.classify_strong_pips()

        return self._fanout_plot(
            "plot_strong_pip", fixed_point, ax, prepare=prepare, **scatter_kwargs
        )

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
        return self._fanout_call("compute_pseudoneighbors", fixed_point, **kwargs)

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

        def prepare(trellis: Trellis) -> None:
            if compute and not trellis.pseudoneighbors:
                trellis.compute_pseudoneighbors()

        return self._fanout_plot(
            "plot_pseudoneighbors", fixed_point, ax, prepare=prepare, **scatter_kwargs
        )

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
        return self._fanout_plot(
            "plot_holes", fixed_point, ax, flatten=True, **scatter_kwargs
        )

    # ── hole / partition convenience (per fixed point) ───────────────────────

    def punch_holes(
        self,
        fixed_point: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        **kwargs,
    ):
        """
        Punch the holes of one, several, or every fixed point.

        For each selected fixed point this builds (or reuses) its per-fixed-point
        Trellis and runs :meth:`Trellis.punch_holes`, so each tangle of a nested
        session is punched from its own pseudoneighbor pairs.

        Args:
            fixed_point: A single FixedPoint (→ returns that fixed point's hole
                list), an iterable of them, or None for every fixed point (→
                returns a ``{fixed_point: holes}`` dict).
            **kwargs: Forwarded to :meth:`Trellis.punch_holes` (``pairs``,
                ``epsilon``, ``propagate``, ``verbose``).

        Returns:
            A hole list for a single fixed point, or a dict mapping each fixed
            point to its list.

        Raises:
            AssertionError: Propagated from :meth:`Trellis.punch_holes` when one
                of its two topological invariants fails.
        """
        return self._fanout_call("punch_holes", fixed_point, **kwargs)

    def partition_stable_manifold(
        self,
        fixed_point: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        **kwargs,
    ):
        """
        Partition the stable branches of one, several, or every fixed point.

        For each selected fixed point this builds (or reuses) its per-fixed-point
        Trellis and runs :meth:`Trellis.partition_stable_manifold`, which by
        default covers every stable branch of that trellis on both sides.

        Args:
            fixed_point: A single FixedPoint (→ returns that fixed point's
                result list), an iterable of them, or None for every fixed point
                (→ returns a ``{fixed_point: results}`` dict).
            **kwargs: Forwarded to :meth:`Trellis.partition_stable_manifold`
                (``branch_key``, ``verbose``).

        Returns:
            A StablePartitionResult list for a single fixed point, or a dict
            mapping each fixed point to its list.
        """
        return self._fanout_call(
            "partition_stable_manifold", fixed_point, **kwargs
        )

    def partition_element_for(
        self,
        bridge_id: "BridgeId",
        endpoint: "Endpoint",
        side: "Side",
    ) -> Optional[int]:
        """
        The partition element at one end of a bridge, searched across trellises.

        A session keeps one trellis per tangle, and
        :meth:`Trellis.element_for` only reads its own partitions, so no single
        trellis can answer for an endpoint whose stable branch belongs to a
        different one. This scans every cached (non-stale) trellis and returns
        the first element found. Bridge ids and registry ids are session-wide
        (one registry), so there is no ambiguity about which bridge is meant.

        This is what resolves a HETEROCLINIC bridge. A bridge of W^u(fp1) cut at
        crossings with W^s(fp3) belongs to fp1 (a bridge is a piece of unstable
        manifold), so it appears only in fp1's trellis — whose partitions cover
        fp1's stable branches, not the fp3 branches its endpoints actually sit
        on. The element lives in fp3's trellis, which does not hold the bridge.
        Both halves are covered because :meth:`Trellis.element_for` falls back to
        the endpoint id carried by the BridgeId itself rather than requiring the
        bridge object to be in the snapshot.

        Args:
            bridge_id: The bridge's :data:`~tanglepack.numerics.Bridge.BridgeId`.
            endpoint: ``"first"`` or ``"second"`` — which end of the id.
            side: Which side's partition to read.

        Returns:
            The element id in the trellis whose branch carries that endpoint, or
            None when no cached trellis has partitioned it.

        Raises:
            ValueError: If ``endpoint`` is neither "first" nor "second" — checked
                here, so an empty trellis cache cannot swallow a typo.

        Note:
            Only trellises already cached are searched, and only while they are
            still valid for the current workbench generation: a stale snapshot's
            partitions are keyed by ids the workbench may have renumbered, so
            answering from one would be worse than answering None. Call
            :meth:`partition_stable_manifold` first to populate the cache.
        """
        endpoint_index(endpoint)  # reject a bad name even with nothing cached
        for trellis in self._trellises.values():
            if self._is_stale(trellis):
                continue
            element_id = trellis.element_for(bridge_id, endpoint, side)
            if element_id is not None:
                return element_id
        return None

    def plot_stable_partition(
        self,
        fixed_point: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        *,
        ax=None,
        **line_kwargs,
    ) -> list:
        """
        Draw the stable partitions of one, several, or every fixed point.

        Partitions must have been built already
        (:meth:`partition_stable_manifold`); trellises without them are skipped.

        Args:
            fixed_point: Fixed point selector; None (default) plots every one.
            ax: Optional matplotlib Axes (defaults to the current axes).
            **line_kwargs: Forwarded to :meth:`Trellis.plot_stable_partition`.

        Returns:
            List of the Axes drawn on (one per fixed point with partitions).
        """
        return self._fanout_plot(
            "plot_stable_partition", fixed_point, ax, **line_kwargs
        )

    def describe_holes(
        self,
        fixed_point: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
    ) -> str:
        """
        Human-readable hole report for one, several, or every fixed point.

        Args:
            fixed_point: Fixed point selector; None (default) reports on every one.

        Returns:
            The per-trellis reports, each under a header naming its fixed point.
        """
        return self._describe("describe_holes", fixed_point)

    def describe_stable_partitions(
        self,
        fixed_point: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
    ) -> str:
        """
        Human-readable partition report for one, several, or every fixed point.

        Args:
            fixed_point: Fixed point selector; None (default) reports on every one.

        Returns:
            The per-trellis reports, each under a header naming its fixed point.
        """
        return self._describe("describe_stable_partitions", fixed_point)

    def _describe(
        self,
        method: str,
        fixed_point: "Optional[FixedPoint | Iterable[FixedPoint]]",
    ) -> str:
        """Concatenate one Trellis describe_* report per selected fixed point."""
        sections = []
        for fp in self._resolve_fixed_points(fixed_point):
            sections.append(f"{fp!r}:")
            sections.append(getattr(self.trellis(fp), method)())
        return "\n".join(sections)

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
        outer one. A recompute advances the workbench generation, so any cached
        Trellis is rebuilt on the next :meth:`trellis` call.

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
            recomputed = self.workbench.intersection_registry.all_ids()
            for rz in self.resonance_zones.values():
                rz.intersection_ids = recomputed
        return self.resonance_zones

    # ── bridge ↔ resonance-zone classification ───────────────────────────────

    @staticmethod
    def _bridge_test_point(bridge: "Bridge") -> Optional[NDArray[np.float64]]:
        """The midpoint node of a bridge, used as its representative point.

        Returns the geometric middle node so a bridge that forms a zone's own unstable
        boundary arc lands exactly on that zone's boundary (and, with boundary-inclusive
        containment, is attributed to it). ``None`` for an empty bridge. The walk
        behind it is memoised on the bridge, so classifying every bridge of a tangle
        costs one walk each, not one per zone.
        """
        return polyline_midpoint(bridge.get_point_array())

    def classify_bridge(self, bridge: "Bridge") -> Optional[ResonanceZone]:
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
            the workbench generation advances and any cached Trellis is rebuilt
            on the next :meth:`trellis` call — take a fresh one (and
            re-establish its strong pip) before continuing topological work.
        """
        return blast_zone(
            self,
            zone,
            num_iterations,
            fixed_point=fixed_point,
            strict=strict,
            min_separation=min_separation,
        )

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
