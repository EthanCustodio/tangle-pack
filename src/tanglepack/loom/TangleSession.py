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
from ..numerics.ManifoldMachine import DEFAULT_AREA_CUTOFF
from ..numerics.budget import DEFAULT_MAX_POINTS_PER_STEP
from ..numerics.geometry import polyline_midpoint
from ..topology import plotting
from ..topology.Arrangement import Arrangement
from ..topology.BridgeClass import BridgeClassTable
from ..topology.DualGraph import DualGraph
from ..topology.MinimalTrellis import MinimalTrellis
from ..topology.PartitionFamily import HomotopyPartition, IteratedHomotopyPartition
from ..topology.SymbolicDynamics import SymbolicDynamics
from ..topology.TopologyResults import endpoint_index
from ..topology.Trellis import Trellis, _is_single_fixed_point
from .ResonanceZone import ResonanceZone, define_resonance_zone, ring_cycles
from .Blast import BlastResult, blast_zone
from .BridgeAlphabet import BridgeAlphabet
from .Topology import Topology

if TYPE_CHECKING:
    import networkx
    from matplotlib.axes import Axes
    from matplotlib.lines import Line2D
    from matplotlib.table import Table

    from ..numerics.Bridge import Bridge, BridgeId
    from ..numerics.FixedPoint import FixedPoint
    from ..numerics.Intersection import Intersection
    from ..topology.DualWalk import Walk
    from ..topology.plotting import CartoonLayout
    from ..topology.TopologyResults import Endpoint, Side

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
        # as a crossing (its cdists and keys), never as an id (API-1).
        self._strong_pips: dict["FixedPoint", "Intersection"] = {}
        # The letters of the bridge classes seen so far, keyed by class identity
        # (Topology.identity). Outlives every recompute so a class keeps its
        # letter across growth and blasts and a new class takes the next one.
        self.bridge_alphabet = BridgeAlphabet()
        # One resonance zone per fixed point (C8): an inversion point's zone is the
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
        single-fixed-point trellis gets its fixed point's chosen strong pip back
        (API-1): the remembered pip is found again by its cdists and keys and
        tested alone (the other candidates are classified only when read, WP10);
        a pip that is no longer a candidate is left unset with a WARNING. Once
        a pip is chosen the trellis is cut at its orbit
        (``Trellis.cut_at_pip_orbit``, D10, D12): a single fixed point's when
        the pip is set or restored, a selection of several at the pips of its
        fixed points on every read (a no-op while they are unchanged), so the
        snapshot is a function of the generation and the pips alone, never of
        which product ran first (D1). A trellis holds a crossing only when both
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
        The topology of a selection: stages 3-9 of the chalkboard, recomputed (D1).

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
                fixed points selected, and so does a heteroclinic pair (D7: a
                pair of ``W^u(A)`` crossings on B's stable branch, holed in B's
                partition, its backward chain holing A's).
            rebuild: Force a recompute at the same state.

        Returns:
            The topology.

        Raises:
            ValueError: A selected fixed point's chosen strong pip is no
                longer a candidate (choose another with :meth:`set_strong_pip`).
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

        A stage of :meth:`topology`. D2: the alphabet is the classes of the
        minimal trellis's ``kept_bridge_ids``; ``minimal_trellis().classes``
        holds every bridge's class (the evidence table, see
        :func:`tanglepack.topology.BridgeClass.bridge_classes` for the class
        semantics). Every ACTIVE class gets its letter from
        :attr:`bridge_alphabet` by its identity (the bounding crossings of its
        two elements, C2): a class seen before keeps its letter across
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

    def plot_transition_graph(
        self,
        dynamics: Optional[SymbolicDynamics] = None,
        ax: Optional["Axes"] = None,
        *,
        fixed_points: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        **kwargs,
    ) -> "networkx.DiGraph":
        """
        Draw a transition graph, defaulting to :meth:`symbolic_dynamics`'s own
        (cached) result for ``fixed_points``.

        Thin delegate to
        :func:`~tanglepack.topology.plotting.plot_transition_graph`.

        Args:
            dynamics: The dynamics to draw. Defaults to
                ``self.symbolic_dynamics(fixed_points)``.
            ax: Optional matplotlib Axes to draw on. Defaults to the current
                axes (plt).
            fixed_points: The trellis selection used when ``dynamics`` is
                omitted.
            **kwargs: Forwarded to
                :func:`~tanglepack.topology.plotting.plot_transition_graph`
                (e.g. ``refined``).

        Returns:
            The ``networkx.DiGraph`` drawn, as the plotting function returns it.
        """
        target = dynamics if dynamics is not None else self.symbolic_dynamics(fixed_points)
        return plotting.plot_transition_graph(target, ax=ax, **kwargs)

    def plot_walk(
        self,
        walk: "Walk",
        ax: Optional["Axes"] = None,
        *,
        fixed_points: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        **kwargs,
    ) -> "Optional[Line2D]":
        """
        Highlight one dual-graph walk over the plane, on
        :meth:`dual_graph`'s (cached) graph for ``fixed_points``.

        Thin delegate to :func:`~tanglepack.topology.plotting.plot_walk`. The
        walk is normally one of a class's
        :attr:`~tanglepack.topology.SymbolicDynamics.ClassDynamics.search`
        results from :meth:`symbolic_dynamics`.

        Args:
            walk: The :class:`~tanglepack.topology.DualWalk.Walk` to draw.
            ax: Optional matplotlib Axes to draw on. Defaults to the current
                axes (plt).
            fixed_points: The trellis selection whose dual graph the walk
                belongs to.
            **kwargs: Forwarded to
                :func:`~tanglepack.topology.plotting.plot_walk` (line
                overrides).

        Returns:
            The ``Line2D`` of the walk, or None when the walk visits no face,
            as the plotting function returns it.
        """
        graph = self.dual_graph(fixed_points)
        return plotting.plot_walk(graph, walk, ax=ax, **kwargs)

    def plot_bridges_by_class(
        self,
        dynamics: Optional[SymbolicDynamics] = None,
        ax: Optional["Axes"] = None,
        *,
        fixed_points: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        **kwargs,
    ) -> dict[str, object]:
        """
        Draw every classed bridge coloured by its (refined) class, defaulting
        to :meth:`symbolic_dynamics`'s own (cached) result for ``fixed_points``.

        Thin delegate to
        :func:`~tanglepack.topology.plotting.plot_bridges_by_class`, drawn over
        the (cached) :meth:`trellis` of the same selection.

        Args:
            dynamics: The dynamics whose classes colour the bridges. Defaults
                to ``self.symbolic_dynamics(fixed_points)``.
            ax: Optional matplotlib Axes to draw on. Defaults to the current
                axes (plt).
            fixed_points: The trellis selection.
            **kwargs: Forwarded to
                :func:`~tanglepack.topology.plotting.plot_bridges_by_class`
                (e.g. ``refined``, ``show_inert``).

        Returns:
            ``{symbol name: colour}`` for every symbol drawn, in legend order,
            as the plotting function returns it.
        """
        trellis = self.trellis(fixed_points)
        target = dynamics if dynamics is not None else self.symbolic_dynamics(fixed_points)
        return plotting.plot_bridges_by_class(trellis, target, ax=ax, **kwargs)

    def plot_itinerary_table(
        self,
        dynamics: Optional[SymbolicDynamics] = None,
        ax: Optional["Axes"] = None,
        *,
        fixed_points: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        **kwargs,
    ) -> "Optional[Table]":
        """
        Tabulate a symbolic dynamics (one row per class), defaulting to
        :meth:`symbolic_dynamics`'s own (cached) result for ``fixed_points``.

        Thin delegate to
        :func:`~tanglepack.topology.plotting.plot_itinerary_table`.

        Args:
            dynamics: The dynamics to tabulate. Defaults to
                ``self.symbolic_dynamics(fixed_points)``.
            ax: Optional matplotlib Axes to draw on. Defaults to the current
                axes (plt).
            fixed_points: The trellis selection used when ``dynamics`` is
                omitted.
            **kwargs: Forwarded to
                :func:`~tanglepack.topology.plotting.plot_itinerary_table`
                (e.g. ``refined``).

        Returns:
            The matplotlib ``Table``, or None when the dynamics has no class,
            as the plotting function returns it.
        """
        target = dynamics if dynamics is not None else self.symbolic_dynamics(fixed_points)
        return plotting.plot_itinerary_table(target, ax=ax, **kwargs)

    def plot_dual_graph_cartoon(
        self,
        dynamics: Optional[SymbolicDynamics] = None,
        ax: Optional["Axes"] = None,
        *,
        fixed_points: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        **kwargs,
    ) -> "CartoonLayout":
        """
        Draw the dual graph as a chalkboard cartoon (elements as bars along a
        line, one node each, bridges as arcs, itineraries element to element)
        for ``fixed_points``, on :meth:`dual_graph`'s (cached) graph and
        :meth:`symbolic_dynamics`'s (cached) dynamics.

        Thin delegate to
        :func:`~tanglepack.topology.plotting.plot_dual_graph_cartoon`.

        Args:
            dynamics: The dynamics to draw. Defaults to
                ``self.symbolic_dynamics(fixed_points)``.
            ax: Optional matplotlib Axes to draw on. Defaults to the current
                axes (plt).
            fixed_points: The trellis selection.
            **kwargs: Forwarded to
                :func:`~tanglepack.topology.plotting.plot_dual_graph_cartoon`
                (``show_bridges``, ``show_walks``, ``walk_linewidths``, ...).

        Returns:
            The :class:`~tanglepack.topology.plotting.CartoonLayout` the
            plotting function returns.
        """
        graph = self.dual_graph(fixed_points)
        target = dynamics if dynamics is not None else self.symbolic_dynamics(fixed_points)
        if kwargs.get("shape") == "circle" and not all(
            name in kwargs for name in ("interior_side", "parents", "arc_order")
        ):
            interior_side, parents, arc_order = self.cartoon_zones()
            kwargs.setdefault("interior_side", interior_side)
            kwargs.setdefault("parents", parents)
            kwargs.setdefault("arc_order", arc_order)
        if kwargs.get("shape") == "circle" and "enclosures" not in kwargs and kwargs.get("parents"):
            kwargs["enclosures"] = self.cartoon_enclosures(
                graph.trellis, graph.minimal.kept_bridge_ids, kwargs["parents"]
            )
        return plotting.plot_dual_graph_cartoon(graph, target, ax=ax, **kwargs)

    def cartoon_enclosures(self, trellis: Trellis, bridge_ids: Iterable, parents: dict) -> dict:
        """
        Which nested tangles each bridge's lobe encloses in the plane.

        A bridge between two crossings of ONE stable branch, closed by the
        stable segment between them, bounds a lobe; a fixed point nested in
        that branch's zone (``parents``) is enclosed when the centroid of its
        orbit lies inside the lobe. The circular cartoon routes such a
        bridge's arc around the nested circle's far side, and every other
        inward arc between the stable arc and the nested circle, so the
        cartoon keeps the plane's topology.

        Args:
            trellis: The trellis the bridges belong to (the dual graph's: a
                bridge of a pip-cut trellis may be a joined arc no other
                snapshot holds).
            bridge_ids: The bridges to test (the minimal trellis's).
            parents: ``{id(fixed_point): id(parent fixed point)}`` from
                :meth:`cartoon_zones`.

        Returns:
            ``{bridge_id: frozenset(id(nested fixed point), ...)}`` for the
            bridges enclosing at least one nested fixed point.
        """
        from ..numerics.geometry import arc_polyline, point_in_polygon

        workbench = self.workbench
        registry = trellis.registry
        tol = registry.cdist_tol
        by_id = {id(fp): fp for fp in workbench.fixed_points}
        children: dict = {}
        for child, parent in parents.items():
            if child in by_id:
                children.setdefault(parent, []).append(by_id[child])
        centroids = {
            id(fp): np.mean(
                [np.ravel(np.asarray(c, dtype=float))[:2] for c in fp.coordinates], axis=0
            )
            for group in children.values()
            for fp in group
        }
        enclosures: dict = {}
        for bridge_id in bridge_ids:
            first, second = registry[bridge_id[0]], registry[bridge_id[1]]
            key = first.manifold_b_key
            if key is None or key != second.manifold_b_key:
                continue
            nested = children.get(id(key[0]))
            if not nested or key not in workbench.manifolds:
                continue
            bridge = trellis.bridge_between(*bridge_id)
            assert bridge is not None, f"bridge {bridge_id} is not a bridge of the trellis"
            c_first, c_second = float(first.stable_cdist), float(second.stable_cdist)
            stable = np.asarray(
                arc_polyline(
                    workbench.manifolds[key], min(c_first, c_second), max(c_first, c_second),
                    stability="stable", tol=tol,
                ),
                dtype=float,
            ).reshape(-1, 2)
            if c_second > c_first:
                stable = stable[::-1]  # back from the second crossing to the first
            lobe = np.vstack([np.asarray(bridge.get_point_array(), dtype=float), stable])
            if len(lobe) < 3:
                continue
            inside = frozenset(
                id(fp) for fp in nested if point_in_polygon(centroids[id(fp)], lobe)
            )
            if inside:
                enclosures[bridge_id] = inside
        return enclosures

    def cartoon_zones(self) -> tuple[dict, dict, dict]:
        """
        What the circular dual-graph cartoon needs from the resonance zones.

        For every stable branch bounding a zone, which side of it faces the
        zone's interior: a few points along the branch's boundary arc (anchor
        to cut point) are pushed a short distance to the LEFT of the stable
        dynamical direction (toward the anchor) and tested against the zone's
        captured boundary polygon; the majority decides. And which zone nests
        in which: a fixed point's circle goes inside the smallest OTHER zone
        containing its first orbit point. And how each zone's stable arcs go
        around it: the ring order of its boundary (pip_j, stable arc to z_j,
        unstable side from z_j to pip_{j+1}, as ``ResonanceZone`` stitches
        it), its rotational direction (the sign of the boundary polygon's
        area) and the angle of every orbit point about the orbit's centroid.

        Returns:
            ``(interior_side, parents, arc_order)``: ``{branch_key: "left" |
            "right"}`` for every branch bounding a zone, ``{id(fixed_point):
            id(parent fixed point)}`` for the nested fixed points, and
            ``{id(fixed_point): ZoneArcOrder}``. All empty without zones (the
            cartoon then draws each branch's right side inside its circle, in
            ``branch_cycle`` order, and puts the circles side by side).
        """
        from ..numerics.geometry import arc_polyline

        workbench = self.workbench
        tol = workbench.intersection_registry.cdist_tol
        interior: dict = {}
        zones = sorted(self.resonance_zones.values(), key=lambda zone: zone.area)
        for zone in zones:
            vertices = zone.boundary_vertices
            if vertices is None or len(vertices) < 3:
                continue
            scale = float(np.linalg.norm(vertices.max(axis=0) - vertices.min(axis=0)))
            eps = 1e-3 * max(scale, 1e-9)
            for cut in zone.boundary_intersections:
                key = cut.manifold_b_key
                if key is None or key in interior or key not in workbench.manifolds:
                    continue
                points = np.asarray(
                    arc_polyline(
                        workbench.manifolds[key], 0.0, cut.stable_cdist,
                        stability="stable", tol=tol,
                    ),
                    dtype=float,
                )
                votes = {"left": 0, "right": 0}
                count = len(points)
                for i in sorted({int(count * f) for f in (0.3, 0.45, 0.6, 0.75)}):
                    if not 1 <= i < count:
                        continue
                    toward_anchor = points[i - 1] - points[i]
                    norm = float(np.linalg.norm(toward_anchor))
                    if norm == 0.0:
                        continue
                    left = np.array([-toward_anchor[1], toward_anchor[0]]) / norm
                    in_left = zone.contains_point(points[i] + eps * left)
                    in_right = zone.contains_point(points[i] - eps * left)
                    if in_left != in_right:
                        votes["left" if in_left else "right"] += 1
                if votes["left"] == votes["right"]:
                    logger.info(
                        "cartoon zones: no side of stable branch %s is clearly inside "
                        "its zone; the cartoon default is used",
                        key[1:],
                    )
                    continue
                interior[key] = "left" if votes["left"] > votes["right"] else "right"
        from ..numerics.geometry import signed_polygon_area
        from ..topology.plotting import ZoneArcOrder

        arc_order: dict = {}
        for zone in zones:
            order = _zone_ring_order(zone)
            vertices = zone.boundary_vertices
            if not order or vertices is None or len(vertices) < 3:
                continue
            fixed_point = zone.fixed_point
            points = [np.ravel(np.asarray(c, dtype=float))[:2] for c in fixed_point.coordinates]
            centroid = np.mean(points, axis=0)
            angles = {}
            if len(order) > 1:
                for key in order:
                    offset = points[key[2]] - centroid
                    angles[key] = float(np.arctan2(offset[1], offset[0]))
            arc_order.setdefault(
                id(fixed_point),
                ZoneArcOrder(
                    order=order,
                    clockwise=signed_polygon_area(vertices) < 0,
                    anchor_angles=angles,
                ),
            )
        parents: dict = {}
        for zone in zones:
            fixed_point = zone.fixed_point
            point = np.ravel(np.asarray(fixed_point.coordinates[0], dtype=float))[:2]
            containers = [
                other
                for other in zones
                if other.fixed_point is not fixed_point
                and other.area > zone.area
                and other.contains_point(point)
            ]
            if containers:
                parent = min(containers, key=lambda other: other.area).fixed_point
                parents[id(fixed_point)] = id(parent)
        return interior, parents, arc_order

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

    def plot_minimal_trellis(
        self,
        minimal: Optional[MinimalTrellis] = None,
        ax: Optional["Axes"] = None,
        **kwargs,
    ) -> "Axes":
        """
        Draw a minimal trellis, defaulting to :meth:`minimal_trellis`'s own
        (all fixed points, cached) result.

        Thin delegate to
        :func:`~tanglepack.topology.plotting.plot_minimal_trellis`.

        Args:
            minimal: The minimal trellis to draw. Defaults to
                ``self.minimal_trellis()``.
            ax: Optional matplotlib Axes to draw on. Defaults to the current
                axes (plt).
            **kwargs: Forwarded to
                :func:`~tanglepack.topology.plotting.plot_minimal_trellis`
                (e.g. ``show_dropped``).

        Returns:
            The Axes drawn on.
        """
        target = minimal if minimal is not None else self.minimal_trellis()
        return plotting.plot_minimal_trellis(target, ax=ax, **kwargs)

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
        crossings. A fixed point whose pip was chosen before gets it back, or
        is left without one (one WARNING) when it is no longer a candidate:
        a chosen pip is never replaced by the default. Otherwise the default
        choice is remembered (API-1).

        Args:
            fixed_point: A single FixedPoint (→ returns that fixed point's candidate
                list), an iterable of them, or None for every fixed point on the
                workbench (→ returns a ``{fixed_point: candidates}`` dict).
            **kwargs: Forwarded to :meth:`Trellis.classify_strong_pips` (``tol``,
                ``collision_rtol``, ``choose_default``; the last only for a
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
        that crossing is tested (F7): the candidates stay unclassified until read."""
        (fixed_point,) = trellis.fixed_points
        chosen = self._strong_pips[fixed_point]
        pip = trellis.registry.find(chosen)
        if pip is not None and trellis._is_candidate(pip):
            trellis.strong_pip = pip  # tested just now: set_strong_pip would test it again
            return
        if warn:
            logger.warning(
                "the chosen strong pip of %r (stable cdist %.6g on %s) is no longer a "
                "strong-pip candidate; choose another with set_strong_pip()",
                fixed_point, chosen.stable_cdist, chosen.manifold_b_key[1:],
            )

    def set_strong_pip(self, fixed_point: "FixedPoint", intersection_id: int) -> int:
        """Choose ``intersection_id`` as the strong pip of ``fixed_point``, for
        every later generation too (it is remembered as a crossing)."""
        trellis = self.trellis(fixed_point)
        trellis.set_strong_pip(intersection_id)
        self._strong_pips[fixed_point] = trellis.strong_pip_intersection
        return intersection_id

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
                self.classify_strong_pips(*trellis.fixed_points)

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
                through :meth:`classify_strong_pips` (the default pip, unless one
                was chosen before: a lost chosen pip stays unset and is not drawn).
                Pass False to plot only trellises that already have one.
            ax: Optional matplotlib Axes (defaults to the current axes).
            **scatter_kwargs: Forwarded to :meth:`Trellis.plot_strong_pip`.

        Returns:
            List of the matplotlib handles drawn (one per fixed point with a strong
            pip).
        """

        def prepare(trellis: Trellis) -> None:
            if classify and trellis.strong_pip is None:
                self.classify_strong_pips(*trellis.fixed_points)

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
                (``extend``, ``collision_rtol``, ``tol``).

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
                ``epsilon``, ``propagate``, ``diagnose``, ``verbose``).

        Returns:
            A hole list for a single fixed point, or a dict mapping each fixed
            point to its list.

        Raises:
            AssertionError: Propagated from :meth:`Trellis.punch_holes` when one
                of its two topological invariants fails (I1 on the carried
                sides, I2 on the crossing signs); a geometric disagreement is a
                WARNING, never an error (D6).
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

        A HETEROCLINIC bridge (a piece of W^u(fp1) between crossings with
        W^s(fp3)) exists only in a trellis holding both tangles (a tangle alone
        holds no heteroclinic crossing), so it is answered by the selection's
        trellis once ``topology()`` has partitioned it;
        :meth:`Trellis.element_for` reads the endpoint id carried by the
        BridgeId itself, never the bridge object.

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

        The zone is stored in :attr:`resonance_zones` keyed by its fixed point
        (one zone per fixed point, inversion included), so a later call for a
        different fixed point does not overwrite it. The zone's pip becomes
        that fixed point's chosen strong pip (the zone is cut at it).
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
        self._strong_pips[rz.key] = pip
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
        zone: "ResonanceZone | FixedPoint",
        num_iterations: int,
        *,
        fixed_point: "Optional[FixedPoint | list[FixedPoint]]" = None,
        min_separation: Optional[float] = None,
        diagnose: bool = True,
        max_points: Optional[int] = None,
        predict: bool = True,
    ) -> BlastResult:
        """
        Repeatedly iterate the bridges inside a resonance zone.

        Each step maps every un-iterated bridge whose INNERMOST zone is this one
        (:meth:`classify_bridge`, C4) forward (cutting the image into child bridges
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
                distance of an already-kept sibling, to avoid precision-driven merges
                of near-coincident unstable curves (default: disabled).
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
            Any failure of one bridge's image propagates (NUM-10). A
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


def _zone_ring_order(zone: ResonanceZone) -> list:
    """
    The stable branch keys of a zone's boundary in ring order (see
    :func:`~tanglepack.loom.ResonanceZone.ring_cycles`).

    Returns:
        The branch keys, cycle after cycle; empty when a cycle does not close
        (iterates missing).
    """
    cuts = zone.boundary_intersections
    cycles = ring_cycles(cuts)
    if any(cuts[cycle[-1]].manifold_b_key[2:] != cuts[cycle[0]].manifold_a_key[2:] for cycle in cycles):
        return []
    return [cuts[index].manifold_b_key for cycle in cycles for index in cycle]
