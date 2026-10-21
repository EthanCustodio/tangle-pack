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
from typing import Iterable, Literal, Optional, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from ..numerics.TangleWorkbench import TangleWorkbench
from ..numerics.DynamicalSystem import MapFunc, JacFunc
from ..numerics.ManifoldMachine import DEFAULT_AREA_CUTOFF
from ..numerics.budget import DEFAULT_MAX_POINTS_PER_STEP
from ..numerics.geometry import polyline_midpoint
from ..topology.Arrangement import Arrangement
from ..topology.BridgeClass import BridgeClassTable
from ..topology.DualGraph import DualGraph
from ..topology.MinimalTrellis import MinimalTrellis
from ..topology.PartitionFamily import HomotopyPartition, IteratedHomotopyPartition
from ..topology.SymbolicDynamics import SymbolicDynamics
from ..topology.Trellis import Trellis, _is_single_fixed_point
from .ResonanceZone import ResonanceZone, define_resonance_zone
from .Blast import BlastResult, blast_zone
from .BridgeAlphabet import BridgeAlphabet
from .Topology import Topology

if TYPE_CHECKING:
    from ..numerics.Bridge import Bridge
    from ..numerics.FixedPoint import FixedPoint
    from ..numerics.Intersection import Intersection

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

    ``compute_intersections`` REPLACES the workbench's intersection registry
    and returns the new one: re-fetch ``session.intersection_registry`` after
    every call (growth drivers and resonance-zone trims call it too) instead of
    holding a reference across one.

    Attributes:
        workbench: The underlying TangleWorkbench.

    Example:
        >>> session = TangleSession(my_map, my_map_inv, my_jac)
        >>> fp = session.construct_fixed_point([4, -4])          # delegated to workbench
        >>> session.initialize_both_manifolds(fp)
        >>> session.grow_n_times(fp, "unstable", num_iterations=7)
        >>> session.grow_until_turnaround(fp, "stable")
        >>> session.compute_intersections(fp)
        >>> session.classify_strong_pips()                        # default pip, remembered
        >>> rz = session.resonance_zone(session.strong_pip(fp))   # trim + recompute
        >>> topo = session.topology()                             # holes .. words, one pass
        >>> topo.classes, topo.dual_graph, topo.dynamics
    """

    def __init__(
        self,
        dynamical_map: MapFunc,
        dynamical_map_inverse: MapFunc,
        jacobian_function: Optional[JacFunc] = None,
        *,
        area_cutoff: float = DEFAULT_AREA_CUTOFF,
        max_points_per_step: Optional[int] = DEFAULT_MAX_POINTS_PER_STEP,
        escape_radius: Optional[float] = None,
    ):
        """
        Build a session around a fresh workbench for the given map.

        Args:
            dynamical_map: The forward map.
            dynamical_map_inverse: Its inverse.
            jacobian_function: Optional analytic Jacobian.
            area_cutoff: The workbench's refinement threshold (1e-7 by
                default, see :data:`~.numerics.ManifoldMachine.DEFAULT_AREA_CUTOFF`).
            max_points_per_step: The workbench's per-step point budget; ``None``
                disables it.
            escape_radius: ``max|coordinate|`` above which a growth step warns;
                ``None`` (default) never warns.
        """
        self.workbench = TangleWorkbench(
            dynamical_map,
            dynamical_map_inverse,
            jacobian_function,
            area_cutoff=area_cutoff,
            max_points_per_step=max_points_per_step,
            escape_radius=escape_radius,
        )
        # One Trellis snapshot and one Topology per fixed-point selection.
        self._trellises: dict[frozenset, Trellis] = {}
        self._topologies: dict[frozenset, Topology] = {}
        # The chosen strong pip of each fixed point, carried across generations
        # as a crossing (its cdists and keys), never as an id.
        self._strong_pips: dict["FixedPoint", "Intersection"] = {}
        # The letters of the bridge classes seen so far, keyed by class identity
        # (Topology.identity). Outlives every recompute so a class keeps its
        # letter across growth and blasts and a new class takes the next one.
        self.bridge_alphabet = BridgeAlphabet()
        # One resonance zone per fixed point: an inversion point's zone is the
        # union of its two loops. Insertion order is preserved so plotting/shading
        # is deterministic.
        self.resonance_zones: dict["FixedPoint", ResonanceZone] = {}
        # The last blast iteration's growth factor (new points per frontier
        # point) per zone key: blast_zone(max_points=...) predicts a call's
        # first iteration from it, since blasts usually come one iteration at
        # a time. Stale once the zone changes; the budget is opt-in.
        self._blast_stats: dict["FixedPoint", float] = {}

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

    # ── workbench settings ───────────────────────────────────────────────────
    #
    # __getattr__ delegates READS only: without these properties
    # ``session.area_cutoff = x`` would set a dead attribute on the session.

    @property
    def area_cutoff(self) -> float:
        """The workbench's refinement threshold (read and set through)."""
        return self.workbench.area_cutoff

    @area_cutoff.setter
    def area_cutoff(self, value: float) -> None:
        self.workbench.area_cutoff = value

    @property
    def max_points_per_step(self) -> Optional[int]:
        """The workbench's per-step point budget (read and set through)."""
        return self.workbench.max_points_per_step

    @max_points_per_step.setter
    def max_points_per_step(self, value: Optional[int]) -> None:
        self.workbench.max_points_per_step = value

    @property
    def escape_radius(self) -> Optional[float]:
        """The workbench's escape WARNING radius (read and set through)."""
        return self.workbench.escape_radius

    @escape_radius.setter
    def escape_radius(self, value: Optional[float]) -> None:
        self.workbench.escape_radius = value

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
        the next call. No caller has to invalidate by hand. A rebuilt
        single-fixed-point trellis gets its fixed point's chosen strong pip back:
        the remembered pip is found again by its cdists and keys and
        tested alone (the other candidates are classified only when read);
        a pip that is no longer a candidate is left unset with a WARNING. Once
        a pip is chosen the trellis is cut at its orbit
        (``Trellis.cut_at_pip_orbit``): a single fixed point's when
        the pip is set or restored, a selection of several at the pips of its
        fixed points on every read (a no-op while they are unchanged), so the
        snapshot is a function of the generation and the pips alone, never of
        which product ran first. A trellis holds a crossing only when both
        its branches are selected: a fixed point's own trellis holds no
        heteroclinic crossing (its arcs are joined across them), the
        selection's holds them.

        Args:
            fixed_points: A single FixedPoint, an iterable of them, or None for all
                of the workbench's fixed points. A selection is a SET (by
                identity): ``trellis(fp)`` and ``trellis([fp])`` are one
                snapshot, its fixed points in workbench order.
            rebuild: Force a rebuild even if a cached Trellis exists.

        Returns:
            The Trellis for the requested fixed point(s).
        """
        key, selected = self._selection(fixed_points)
        trellis = self._trellises.get(key)
        if rebuild or trellis is None or self._is_stale(trellis):
            trellis = self._trellises[key] = Trellis.from_workbench(self.workbench, selected)
            if len(selected) == 1 and selected[0] in self._strong_pips:
                self._restore_strong_pip(trellis)
        if len(selected) > 1:
            trellis.cut_at_pip_orbit(
                pip for fp in selected if (pip := self.trellis(fp).strong_pip) is not None
            )
        return trellis

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
        single-fixed-point arrangement holds neither it nor the other tangle's
        arcs.

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

    def topology(
        self,
        fixed_points: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        *,
        rebuild: bool = False,
    ) -> Topology:
        """
        The topology of a selection: stages 3-9 of the chalkboard, recomputed.

        One :class:`~tanglepack.loom.Topology.Topology` per selection, rebuilt
        when the workbench generation or the selection's strong pips change:
        each selected fixed point's pip is chosen on its own trellis (the
        default pip, remembered, when none was chosen), then pseudoneighbors,
        holes and the stable partition are recomputed from scratch on the
        SELECTION's trellis, and the bridge classes, minimal trellis, iterated
        partition, dual graph and words are its lazy stages. A manual stage
        call on a per-fixed-point trellis (``punch_holes(...)`` with custom
        arguments, say) fills that trellis for plotting but never feeds a
        product; after ``topology()`` of several fixed points the
        per-fixed-point trellises carry no stage output (``topology(fp)``
        fills one).

        Args:
            fixed_points: A single FixedPoint, an iterable of them, or None (the
                default) for all of them. Partitions, holes and pips come from
                the selection only, so a heteroclinic bridge needs both its
                fixed points selected, and so does a heteroclinic pair (a
                pair of ``W^u(A)`` crossings on B's stable branch, holed in B's
                partition, its backward chain holing A's).
            rebuild: Force a recompute at the same state.

        Returns:
            The topology.

        Raises:
            ValueError: A selected fixed point's chosen strong pip is no
                longer a candidate (choose another with :meth:`set_strong_pip` or drop
                it with :meth:`forget_strong_pip`).
        """
        key, selected = self._selection(fixed_points)
        topology = self._topologies.get(key)
        if (
            rebuild
            or topology is None
            or topology.generation != self.workbench.generation
            or topology.strong_pips != tuple(self.strong_pip(fp) for fp in selected)
        ):
            topology = self._topologies[key] = Topology(self, selected)
        return topology

    def bridge_classes(
        self,
        fixed_points: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        *,
        rebuild: bool = False,
    ) -> BridgeClassTable:
        """
        The ALPHABET of a selection: its minimal trellis's classes, lettered and zoned.

        A stage of :meth:`topology`. The alphabet is the classes of the
        minimal trellis's ``kept_bridge_ids``; ``minimal_trellis().classes``
        holds every bridge's class (the evidence table, see
        :func:`tanglepack.topology.BridgeClass.bridge_classes` for the class
        semantics). Every ACTIVE class gets its letter from
        :attr:`bridge_alphabet` by its identity (the bounding crossings of its
        two elements): a class seen before keeps its letter across
        recomputes, a new one takes the next unused, inert classes stay
        unlettered. Every class records the resonance zone its member bridges
        lie in (:meth:`classify_bridge`, unanimous over the members, else a
        WARNING and no zone).

        Args:
            fixed_points: The selection, as for :meth:`topology`.
            rebuild: Recompute the topology first.

        Returns:
            A :class:`~tanglepack.topology.BridgeClass.BridgeClassTable` in
            table order (grouped by tangle, smallest member cdist first).
        """
        return self.topology(fixed_points, rebuild=rebuild).classes

    def describe_bridge_classes(
        self,
        fixed_points: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
    ) -> str:
        """
        A multi-line report of :meth:`bridge_classes`, one line per class.

        Args:
            fixed_points: The trellis selection, as for :meth:`bridge_classes`.

        Returns:
            :meth:`~tanglepack.topology.BridgeClass.BridgeClassTable.describe`
            of the table.
        """
        return self.bridge_classes(fixed_points).describe()

    def homotopy_partition(
        self, fixed_points: "Optional[FixedPoint | Iterable[FixedPoint]]" = None
    ) -> HomotopyPartition:
        """
        The selection's stable partitions as a
        :class:`~tanglepack.topology.PartitionFamily.HomotopyPartition`; a
        stage of :meth:`topology`.

        Args:
            fixed_points: The selection, as for :meth:`topology`.

        Returns:
            The family.
        """
        return self.topology(fixed_points).homotopy

    def minimal_trellis(
        self,
        fixed_points: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        *,
        rebuild: bool = False,
    ) -> MinimalTrellis:
        """
        The :class:`~tanglepack.topology.MinimalTrellis.MinimalTrellis` of a
        selection: its hole bridges and the first images of those in an
        active class; a stage of :meth:`topology`. Activity is read from the
        evidence table (every bridge's class), kept as
        ``MinimalTrellis.classes``.

        Args:
            fixed_points: The selection, as for :meth:`topology`.
            rebuild: Recompute the topology first.

        Returns:
            The minimal trellis.
        """
        return self.topology(fixed_points, rebuild=rebuild).minimal

    def iterated_partition(
        self,
        fixed_points: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        *,
        rebuild: bool = False,
    ) -> IteratedHomotopyPartition:
        """
        The :class:`~tanglepack.topology.PartitionFamily.IteratedHomotopyPartition`:
        the homotopy partition cut by the image bridges of
        :meth:`minimal_trellis`; a stage of :meth:`topology`.

        Args:
            fixed_points: The selection, as for :meth:`topology`.
            rebuild: Recompute the topology first.

        Returns:
            The refined family.
        """
        return self.topology(fixed_points, rebuild=rebuild).iterated

    def describe_iterated_partition(
        self, fixed_points: "Optional[FixedPoint | Iterable[FixedPoint]]" = None
    ) -> str:
        """
        Human-readable report of :meth:`iterated_partition`, interval notation.

        Args:
            fixed_points: The trellis selection.

        Returns:
            The report, one line per element with its parent and cutting bridge.
        """
        return self.iterated_partition(fixed_points).describe()

    def dual_graph(
        self,
        fixed_points: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        *,
        rebuild: bool = False,
    ) -> DualGraph:
        """
        The :class:`~tanglepack.topology.DualGraph.DualGraph` of a selection;
        a stage of :meth:`topology`.

        Built over the minimal trellis's sparse arrangement and the iterated
        partition; each selected fixed point's strong pip declares the
        fundamental segment ``(f^k(q0), q0]`` of its own stable branch, where
        the two side nodes of a stable edge unify into one traversable node.

        Args:
            fixed_points: The selection, as for :meth:`topology`.
            rebuild: Recompute the topology first.

        Returns:
            The dual graph.

        Raises:
            ValueError: Propagated from the constructor when an arrangement
                component fails the graph's own topological checks.
        """
        return self.topology(fixed_points, rebuild=rebuild).dual_graph

    def symbolic_dynamics(
        self,
        fixed_points: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        *,
        rebuild: bool = False,
    ) -> SymbolicDynamics:
        """
        The :class:`~tanglepack.topology.SymbolicDynamics.SymbolicDynamics` of
        a selection; a stage of :meth:`topology`.

        Read off :meth:`dual_graph` and :meth:`bridge_classes`; see
        :func:`~tanglepack.topology.SymbolicDynamics.symbolic_dynamics` for the
        algorithm. An active class the closure of the alphabet adds takes its
        letter from :attr:`bridge_alphabet` by identity too, so it persists;
        :meth:`bridge_classes` stays the alphabet handed in (the dynamics'
        ``table`` is the closed one). A class whose image is not registered
        is left unresolved with a WARNING: grow or blast and call again.

        Args:
            fixed_points: The selection, as for :meth:`topology`.
            rebuild: Recompute the topology first.

        Returns:
            The symbolic dynamics.
        """
        return self.topology(fixed_points, rebuild=rebuild).dynamics

    def describe_symbolic_dynamics(
        self,
        fixed_points: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
    ) -> str:
        """
        A multi-line report of :meth:`symbolic_dynamics`.

        Args:
            fixed_points: The trellis selection, as for :meth:`symbolic_dynamics`.

        Returns:
            :meth:`~tanglepack.topology.SymbolicDynamics.SymbolicDynamics.describe`
            of the (cached) dynamics: the naming, every class with its
            itinerary and word, the refinement and the transition matrix.
        """
        return self.symbolic_dynamics(fixed_points).describe()

    def _is_stale(self, trellis: Trellis) -> bool:
        """Whether a cached Trellis snapshot predates the workbench's current state."""
        return trellis._built_generation != self.workbench.generation

    def _selection(self, fixed_points) -> tuple[frozenset, list["FixedPoint"]]:
        """A selection's cache key (its fixed-point identities, a SET) and its
        fixed points in workbench order (None: every fixed point)."""
        key = frozenset(id(fp) for fp in self._resolve_fixed_points(fixed_points))
        return key, [fp for fp in self.workbench.fixed_points if id(fp) in key]

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
        crossings. A fixed point whose pip was chosen before gets it back, or
        is left without one (one WARNING) when it is no longer a candidate:
        a chosen pip is never replaced by the default. Otherwise the default
        choice is remembered.

        Args:
            fixed_point: A single FixedPoint (→ returns that fixed point's candidate
                list), an iterable of them, or None for every fixed point on the
                workbench (→ returns a ``{fixed_point: candidates}`` dict).
            **kwargs: Forwarded to :meth:`Trellis.classify_strong_pips`
                (``collision_rtol``, ``choose_default``; the last only for a
                fixed point with no chosen pip).

        Returns:
            A candidate-id list for a single fixed point, or a dict mapping each
            fixed point to its candidate-id list.
        """
        results = {}
        for fp in self._resolve_fixed_points(fixed_point):
            trellis = self.trellis(fp)
            remembered = fp in self._strong_pips
            found = trellis.strong_pip is not None  # a lost pip was warned at the rebuild
            results[fp] = trellis.classify_strong_pips(
                **({**kwargs, "choose_default": False} if remembered else kwargs)
            )
            if remembered:
                self._restore_strong_pip(trellis, warn=found)
            if trellis.strong_pip is not None:
                self._strong_pips[fp] = trellis.strong_pip_intersection
        if _is_single_fixed_point(fixed_point):
            return results[fixed_point]
        return results

    def _restore_strong_pip(self, trellis: Trellis, *, warn: bool = True) -> None:
        """Choose the remembered pip on a single-fixed-point trellis, leaving it
        unset (with a WARNING when ``warn``) if it is no longer a candidate. Only
        that crossing is tested: the candidates stay unclassified until read."""
        (fixed_point,) = trellis.fixed_points
        chosen = self._strong_pips[fixed_point]
        pip = trellis.registry.find(chosen)
        if pip is not None and trellis._is_candidate(pip):
            trellis.strong_pip = pip  # tested just now: set_strong_pip would test it again
            return
        if warn:
            logger.warning(
                "the chosen strong pip of %r (stable cdist %.6g on %s) is no longer a "
                "strong-pip candidate; choose another with set_strong_pip() or drop it "
                "with forget_strong_pip()",
                fixed_point, chosen.stable_cdist, chosen.manifold_b_key[1:],
            )

    def set_strong_pip(self, fixed_point: "FixedPoint", intersection_id: int) -> int:
        """Choose ``intersection_id`` as the strong pip of ``fixed_point``, for
        every later generation too (it is remembered as a crossing)."""
        trellis = self.trellis(fixed_point)
        trellis.set_strong_pip(intersection_id)
        self._strong_pips[fixed_point] = trellis.strong_pip_intersection
        return intersection_id

    def forget_strong_pip(self, fixed_point: "FixedPoint") -> None:
        """Drop the pip chosen for ``fixed_point``: its trellis is uncut and the
        next :meth:`topology` chooses the default again. The remedy for a pip
        that is no longer a candidate."""
        self._strong_pips.pop(fixed_point, None)
        self.trellis(fixed_point).strong_pip = None

    def strong_pip(self, fixed_point: "FixedPoint") -> Optional[int]:
        """The chosen strong-pip id for ``fixed_point``'s trellis (or None)."""
        return self.trellis(fixed_point).strong_pip

    def chosen_strong_pip(self, fixed_point: "FixedPoint") -> Optional["Intersection"]:
        """The pip the session carries for ``fixed_point`` across generations, as a
        crossing (it may no longer be a candidate), or None when none was chosen."""
        return self._strong_pips.get(fixed_point)

    def strong_pip_candidates(self, fixed_point: "FixedPoint") -> list[int]:
        """The strong-pip candidate ids for ``fixed_point``'s trellis."""
        return self.trellis(fixed_point).strong_pip_candidates

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

        The zone is stored in :attr:`resonance_zones` keyed by its fixed point
        (one zone per fixed point, inversion included), so a later call for a
        different fixed point does not overwrite it. The zone's pip becomes
        that fixed point's chosen strong pip (the zone is cut at it) when it is
        a strong-pip candidate on the trimmed trellis; otherwise the zone is
        still stored, the previous choice stays and a WARNING says so.
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
        pip = self.workbench.intersection_registry[intersection_id]
        rz = define_resonance_zone(
            self.workbench, intersection_id, fixed_points, recompute=recompute
        )
        self.resonance_zones[rz.key] = rz
        # The trimmed trellis is built without the previous choice, which the
        # trim may have cut away: only the zone's pip is tested here.
        previous = self._strong_pips.pop(rz.key, None)
        trellis = self.trellis(rz.key)
        found = trellis.registry.find(pip)
        if found is not None and trellis._is_candidate(found):
            trellis.strong_pip = found  # tested just now, as _restore_strong_pip does
            self._strong_pips[rz.key] = pip
            return rz
        logger.warning(
            "the zone pip %d of %r is not a strong-pip candidate; the zone is "
            "stored, the chosen strong pip is unchanged",
            intersection_id, rz.key,
        )
        if previous is not None:
            self._strong_pips[rz.key] = previous
            self._restore_strong_pip(trellis)
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

    # ── blasting ─────────────────────────────────────────────────────────────

    def blast_zone(
        self,
        zone: "ResonanceZone | FixedPoint",
        num_iterations: int,
        *,
        fixed_point: "Optional[FixedPoint | list[FixedPoint]]" = None,
        min_separation: "Optional[float | Literal['resolution']]" = None,
        diagnose: bool = True,
        max_points: Optional[int] = None,
        predict: bool = True,
    ) -> BlastResult:
        """
        Repeatedly iterate the bridges inside a resonance zone.

        Each step maps every un-iterated bridge whose INNERMOST zone is this one
        (:meth:`classify_bridge`) forward (cutting the image into child bridges
        at its new stable-manifold crossings) and keeps the children that land back
        in it, recursing up to ``num_iterations`` times. An outer zone's blast never
        iterates a bridge inside a smaller stored zone of another fixed point (the
        blasted zone itself need not be stored), so nested blasting does not depend
        on the order of the zones; with ``fixed_point=None`` it still iterates an
        inner tangle's bridges in the annulus between them (they are the outer
        zone's). Bridges that leave the zone are never re-iterated — iterating
        exterior bridges many times grows them exponentially. See
        :mod:`tanglepack.loom.Blast`.

        Args:
            zone: A :class:`ResonanceZone` or its ``(fixed_point, branch_index)`` key.
            num_iterations: Maximum number of blast steps.
            fixed_point: Restrict to bridges from this fixed point (default: all).
            min_separation: Drop a child bridge whose interior comes within this
                distance of curve the workbench already holds or of an
                already-kept sibling, to avoid precision-driven merges of
                near-coincident unstable curves (default: disabled). A drop is
                remembered on the bridge and later calls skip it unless they
                pass a smaller separation or None, so single-iteration blasts
                compose. ``"resolution"`` derives it from the machine's
                ``area_cutoff`` (see :func:`~tanglepack.loom.Blast.blast_zone`).
            diagnose: Finish with one resolution report on the workbench logger
                (default; see :func:`~tanglepack.loom.Blast.blast_zone`).
            max_points: Opt-in point budget of the whole blast; it stops
                BETWEEN iterations on a resumable state with
                ``stop_reason == "point_budget"`` (default: no budget).
            predict: Refuse an iteration whose predicted cost (frontier points
                times the zone's last growth factor) would exceed
                ``max_points`` (default); False stops only once the budget is
                reached.

        Returns:
            A :class:`~tanglepack.loom.Blast.BlastResult` genealogy of the blast.

        Note:
            Blasting registers the children's new stable-manifold crossings, so
            the workbench generation advances and any cached Trellis is rebuilt
            on the next :meth:`trellis` call — take a fresh one (and
            re-establish its strong pip) before continuing topological work.

        Note:
            Any failure of one bridge's image propagates. A
            :class:`~tanglepack.numerics.budget.PointBudgetExceeded` leaves the
            session not resumable; rebuild it.
        """
        return blast_zone(
            self,
            zone,
            num_iterations,
            fixed_point=fixed_point,
            min_separation=min_separation,
            diagnose=diagnose,
            max_points=max_points,
            predict=predict,
        )
