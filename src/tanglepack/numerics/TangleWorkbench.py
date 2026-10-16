"""The orchestrator of the numerical layer.

:class:`TangleWorkbench` owns the objects a computed tangle is made of -- the
dynamical system, the fixed points, the manifolds, the :class:`~.Tangle.Tangle`
index, the :class:`~.IntersectionRegistry.IntersectionRegistry` and the bridge
set -- and drives construction, growth, intersection detection, bridge cutting,
trimming and manifold plotting.

Two collaborators own the work that reads and extends that state without
being part of the orchestration:

* :mod:`~.BridgeIterator` -- mapping a bridge forward, and deriving the bridge
  genealogy from the crossing iterate table;
* :mod:`~.IterateInference` -- filling that iterate table.

Each is reached through a same-named method on the workbench, so no caller has
to know where the body lives. The Dev Notes for each of those areas live in the
module that owns it.

Dev Notes:

The ``grow_until_*`` drivers cap their rounds at ``max_iterations``; the default
(None) scales the cap with the stretch the stop condition needs,
``max(10, ceil(log_stretch / ln(stretch per round)) + 1)``, so a
small ``lambda_u`` gets the depth it needs and no cap falls below 10. ``grow_until_arclength`` stops on the polyline ARC
LENGTH of the grown branch (not the tail's canonical distance, which agrees
only in the linear regime).
"""

from __future__ import annotations

import logging
import math
from collections import Counter
from typing import Callable, Literal, Iterable, Optional, Sequence, Union

import numpy as np
import matplotlib.pyplot as plt

from .BridgeIterator import BridgeIterator
from .DynamicalSystem import DynamicalSystem, MapFunc, JacFunc
from .FixedPointSolver import FixedPointSolver
from .ManifoldInitializer import ManifoldInitializer
from .ManifoldMachine import DEFAULT_AREA_CUTOFF, ManifoldMachine
from .budget import DEFAULT_MAX_POINTS_PER_STEP, GrowthStats
from .Tangle import Tangle
from .FixedPoint import FixedPoint, position_letter
from .BaseManifold import BaseManifold
from .Bridge import Bridge, BridgeId
from .Intersection import Intersection, ManifoldKey, Stability
from .IntersectionRegistry import DEDUP_RTOL, IntersectionRegistry
from .IterateInference import infer_iterates
from .diagnostics import (
    ResolutionReport,
    SameStabilityCount,
    SameStabilityCrossing,
    count_same_stability_crossings,
    same_stability_crossings,
)


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

#: How much a stable segment must be stretched (log of the factor) before it
#: turns around: rounds * ln(stretch per round) measured 14.8-17.1 on the Hénon
#: saddle over k = 0.3 .. 10 from the default ~5e-6 seed, with margin.
TURNAROUND_LOG_STRETCH = 18.0


class MaxIterationsReached(ValueError):
    """
    A ``grow_until_*`` driver hit its iteration cap before its stop condition.

    One growth round stretches a branch by ``stretch_per_round`` (``lambda_u``
    per map step on an area-preserving map), so the depth a stop condition
    needs scales with ``1 / ln(stretch_per_round)``: a small ``lambda_u``
    (Hénon at small ``k``) needs more rounds than the default cap of 10.

    Attributes:
        max_iterations: The cap that was reached.
        stretch_per_round: The canonical-distance factor of one growth round on
            the grown manifold.
        suggested_max_iterations: ``ceil(log_stretch / ln(stretch_per_round))``
            for the stretch the condition is estimated to need, at least one
            more than the cap.
    """

    def __init__(
        self,
        max_iterations: int,
        stretch_per_round: float,
        log_stretch: float,
        what: str,
    ) -> None:
        """
        Build the error and its message.

        Args:
            max_iterations: The cap that was reached.
            stretch_per_round: The grown manifold's factor per growth round.
            log_stretch: Natural log of the total stretch the stop condition is
                estimated to need, counted from the start of the call.
            what: The stop condition that was not met, for the message.
        """
        self.max_iterations = max_iterations
        self.stretch_per_round = stretch_per_round
        self.suggested_max_iterations = max(
            math.ceil(log_stretch / math.log(stretch_per_round)), max_iterations + 1
        )
        super().__init__(
            f"Max iterations ({max_iterations}) reached without {what}. One growth "
            f"round stretches the branch by {stretch_per_round:.4g}, so the depth "
            f"needed scales with 1/ln({stretch_per_round:.4g}); try "
            f"max_iterations={self.suggested_max_iterations}, or stop on a "
            f"condition of your own with grow_until(predicate) or "
            f"grow_until_bounded(radius)."
        )


class TangleWorkbench:
    """
    The entry point for programmatic use of the numerical layer.

    A workbench holds one dynamical system and everything computed from it: the
    fixed points found so far, the manifolds grown from them (keyed by
    ``(fixed_point, stability, orbit_index, branch_index)``), the
    :class:`~.Tangle.Tangle` that indexes those manifolds and resolves their
    crossings, the :class:`~.IntersectionRegistry.IntersectionRegistry` that is
    the single source of truth for the crossings themselves, and the bridge set
    cut out of the unstable manifolds.

    The scientific workflow it drives is: construct a fixed point, orient its
    eigenvectors, initialize its manifolds, grow them, compute the
    intersections, cut the bridges, and iterate those bridges forward::

        wb = TangleWorkbench(my_map, my_map_inverse)
        fp = wb.construct_fixed_point([4, -4])
        wb.orient_eigenvectors(fp, {"unstable": [-1, 0], "stable": [0, 1]})
        wb.initialize_both_manifolds(fp)
        wb.grow_n_times(fp, "unstable", num_iterations=8)
        wb.grow_until_turnaround(fp, "stable")
        wb.compute_intersections(fp)
        wb.create_bridges(fp)

    Everything derived from a workbench -- a :class:`~.topology.Trellis.Trellis`,
    an arrangement, a resonance zone -- is only valid at one
    :attr:`generation`, and a cache checks staleness by comparing that one
    integer.

    Args:
        dynamical_map: The area-preserving map.
        dynamical_map_inverse: Its inverse.
        jacobian_function: Optional analytic Jacobian; a finite-difference one is
            used when it is omitted.
        area_cutoff: Refinement threshold of the one shared
            :class:`~.ManifoldMachine.ManifoldMachine` (see :attr:`area_cutoff`).
        max_points_per_step: Its per-step point budget, ``None`` to disable
            (see :attr:`max_points_per_step`).
        escape_radius: Optional ``max|coordinate|`` above which a growth step
            logs a WARNING; off by default.
        resolution_ratio_threshold: Part of the definition of N: a detection
            whose resolution ratio is below it is no crossing (see the
            constructor).

    Note:
        Bridge iteration and iterate-table inference are delegated to
        :class:`~.BridgeIterator.BridgeIterator` and
        :func:`~.IterateInference.infer_iterates`. The methods here forward to
        them and their signatures are the contract; the bodies (and their Dev
        Notes) live in those modules.
    """

    # Fixed palette of cool hues (blues, cyans, teals, greens, purples) used to
    # color bridges. Deliberately excludes any warm hue near red so a bridge is
    # never mistaken for the stable manifold (always plotted red). The ordering
    # interleaves hue families so consecutive bridges stay visually distinct.
    _COOL_BRIDGE_PALETTE = [
        "#1f77b4",  # blue
        "#2ca02c",  # green
        "#9467bd",  # purple
        "#17becf",  # cyan
        "#0050ef",  # strong blue
        "#008080",  # teal
        "#6a3d9a",  # deep purple
        "#33a02c",  # leaf green
        "#1b9e77",  # blue-green
        "#386cb0",  # slate blue
        "#7570b3",  # indigo
        "#66c2a5",  # mint
        "#5e4fa2",  # violet
        "#3690c0",  # ocean blue
        "#41ab5d",  # emerald
        "#54278f",  # royal purple
    ]

    def __init__(
        self,
        dynamical_map: MapFunc,
        dynamical_map_inverse: MapFunc,
        jacobian_function: Optional[JacFunc] = None,
        *,
        area_cutoff: float = DEFAULT_AREA_CUTOFF,
        max_points_per_step: Optional[int] = DEFAULT_MAX_POINTS_PER_STEP,
        escape_radius: Optional[float] = None,
        resolution_ratio_threshold: Optional[float] = 3.0,
    ) -> None:
        """
        Build an empty workbench around one dynamical system.

        Args:
            dynamical_map (MapFunc): The area-preserving map.
            dynamical_map_inverse (MapFunc): Its inverse.
            jacobian_function (Optional[JacFunc]): Analytic Jacobian; a
                finite-difference one is used when it is omitted.
            area_cutoff (float): Refinement threshold of the shared machine;
                :data:`~.ManifoldMachine.DEFAULT_AREA_CUTOFF` (1e-7, converged)
                by default. Pass a coarser one (1e-4: ~10x fewer points) only
                for a recipe whose topology is known not to depend on it.
            max_points_per_step (Optional[int]): Its per-step point budget;
                ``None`` disables it.
            escape_radius (Optional[float]): ``max|coordinate|`` above which a
                growth step warns; ``None`` (default) never warns.
            resolution_ratio_threshold (Optional[float]): Part of the
                definition of N (design note §A.2): a detection whose
                :attr:`~.Intersection.Intersection.resolution_ratio` is below it
                is registered but no crossing of N -- nothing is cut at it,
                linked to it or trimmed at it, and :meth:`report_resolution`
                names it with the remedy "refine". Derived bound 1: at a ratio
                of 1 the segments separate by no more than the polygon error,
                and the crossing cannot be told apart from noise. The default
                keeps a margin 3 above it; the corpus minimum, over the law
                cases, is 3.038 (heteroclinic). None turns this
                clause off; a run of detections closer than their chord
                deviation is out of N whatever the threshold. It is read by
                every detection pass and every trellis built afterwards: after
                changing it on a workbench that holds crossings, run
                :meth:`compute_intersections` to re-cut.
        """

        self.dynamical_system = DynamicalSystem(
            dynamical_map, dynamical_map_inverse, jacobian_function
        )

        self._fp_solver = FixedPointSolver(self.dynamical_system)
        # ONE machine refines everything, the initializer's fundamental segments
        # included, so area_cutoff and the point budget reach every refinement.
        self._man_machine = ManifoldMachine(
            self.dynamical_system, area_cutoff, max_points_per_step
        )
        self._man_maker = ManifoldInitializer(self.dynamical_system, self._man_machine)
        self.escape_radius = escape_radius
        self._growth_stats: dict[
            tuple[FixedPoint, Stability, Optional[int]], GrowthStats
        ] = {}

        self.Tangle = Tangle()

        self.fixed_points = []
        # manifolds are keyed like (fixed_point, stability, orbit_index, branch_index)
        self._manifolds: dict[tuple[FixedPoint, Stability, int, int], BaseManifold] = {}

        self._intersection_registry = IntersectionRegistry()
        self.resolution_ratio_threshold = resolution_ratio_threshold
        # What report_resolution has already logged, so a repeated pass (a zone
        # recompute, the next blast) only speaks about NEW findings: the flagged
        # crossing ids with the registry id epoch they belong to, and the
        # same-stability hits of a full listing (a blast only counts its own).
        self._reported_flags: tuple[int, frozenset[int]] = (0, frozenset())
        self._reported_artifacts: set[tuple] = set()
        # Bridges are keyed by their topological identity (the ordered endpoint
        # crossing pair), so there is exactly one object per computed bridge and
        # genealogy can be derived rather than stored. Partial pieces have no such
        # identity and are held separately; ``_bridges_at`` is the reverse index
        # from a crossing id to the bridges that end there.
        self._bridges: dict[BridgeId, Bridge] = {}
        self._partial_bridges: list[Bridge] = []
        self._bridges_at: dict[int, list[BridgeId]] = {}
        # Every blast image curve: computed unstable manifold the branches do not
        # hold, re-indexed by every recompute until growth covers it.
        self._image_curves: list[Bridge] = []
        # A BridgeId is a pair of REGISTRY ids, so it only means anything while the
        # registry still numbers its crossings the same way. This counter is bumped
        # every time a recompute renumbers them (see compute_intersections);
        # _bridge_id_epoch is its value when the current bridges were cut, and
        # rebuild_bridges refuses to carry metadata across a mismatch.
        self._registry_id_epoch: int = 0
        self._bridge_id_epoch: int = 0
        # Bridges whose image a PointBudgetExceeded abandoned half-refined; see
        # BridgeIterator.iterate_bridge, which refuses to iterate them again.
        self._half_iterated_bridges: set[Bridge] = set()

        # Generation bookkeeping (see the `generation` property). `_mutations`
        # counts the changes this object makes that nothing else can observe
        # (the bridge set, the registry swap); `_generation_state` is the last
        # composite state seen and `_generation` the monotone token handed out.
        self._mutations: int = 0
        self._generation_state: Optional[tuple] = None
        self._generation: int = 0

        # The collaborator that reads and extends this state (see the module
        # docstring). It holds no state of its own; the workbench forwards its
        # same-named methods to it.
        self._bridge_iterator = BridgeIterator(self)

    # ── refinement and the point budget ──────────────────────────────────────

    @property
    def area_cutoff(self) -> float:
        """
        Curvature area above which the shared machine refines a pair of points.

        Setting it takes effect on the next growth step, initial fundamental
        segments included (the initializer refines with the same machine).
        """
        return self._man_machine.area_cutoff

    @area_cutoff.setter
    def area_cutoff(self, value: float) -> None:
        self._man_machine.area_cutoff = value

    @property
    def max_points_per_step(self) -> Optional[int]:
        """
        Points one growth step (or one bridge image) may create before it raises
        :class:`~.budget.PointBudgetExceeded`; ``None`` disables the budget.
        """
        return self._man_machine.max_points_per_step

    @max_points_per_step.setter
    def max_points_per_step(self, value: Optional[int]) -> None:
        self._man_machine.max_points_per_step = value

    def growth_stats(
        self,
        fixed_point: FixedPoint,
        stability: Stability,
        branch_index: Optional[int] = 0,
    ) -> Optional[GrowthStats]:
        """
        What the last completed growth step of one manifold branch cost.

        Args:
            fixed_point: The fixed point whose manifold grew.
            stability: Which manifold.
            branch_index: The ``branch_index`` the step was grown with (None =
                every branch at once). Ignored on an inversion point, whose one
                chain covers both branches.

        Returns:
            The :class:`~.budget.GrowthStats` of the last step
            :meth:`grow_n_times` completed with that branch selection, or None
            before the first. Steps of different branches are never compared,
            so the growth factor always relates one branch to itself.
        """
        return self._growth_stats.get(
            self._growth_key(fixed_point, stability, branch_index)
        )

    @staticmethod
    def _growth_key(
        fixed_point: FixedPoint, stability: Stability, branch_index: Optional[int]
    ) -> tuple[FixedPoint, Stability, Optional[int]]:
        """The ``growth_stats`` key: an inversion point grows one chain whatever
        branch is asked for."""
        return (
            fixed_point,
            stability,
            None if fixed_point.check_inversion() else branch_index,
        )

    # ── generation ───────────────────────────────────────────────────────────

    @property
    def generation(self) -> int:
        """
        Monotone token that changes whenever anything derived from this
        workbench goes stale.

        One integer comparison is all a cache needs: a :class:`Trellis`, an
        arrangement, or any other snapshot records the generation it was built
        at and is stale as soon as the two differ. The token advances on:

        * any registry mutation -- a crossing added, a renumbering, an iterate
          registered (:attr:`IntersectionRegistry.generation`);
        * a registry swap, i.e. a resetting :meth:`compute_intersections`;
        * any bridge-set mutation -- :meth:`create_bridges`,
          :meth:`clear_bridges`, :meth:`rebuild_bridges`, :meth:`iterate_bridge`;
        * any geometry change a REGISTERED manifold can see -- a ``root``/``tail``
          reassignment on something in :attr:`manifolds` (growth ends in
          ``_find_tail``, a trim and a resonance-zone restore both set tails),
          or an explicit ``bump_version()`` on one of those manifolds. A bridge's
          version is deliberately NOT summed here: bridges come and go with the
          bridge-set mutations already listed, and each one's version is bumped
          in the same breath (see :meth:`_bump_generation`);
        * :meth:`register_manifold`.

        NOT tracked: mutations made behind the library's back -- points spliced
        into a linked list by code driving :class:`ManifoldMachine` directly
        without calling ``bump_version``, an ``Intersection``'s fields edited in
        place, or a ``FixedPoint``'s eigendata replaced. Reads never advance it,
        including the lazy rebuilds of derived caches (the registry's rank
        maps).

        Returns:
            int: The current generation. Compare for equality only; the value
            counts observed changes, not mutations.

        Note:
            This is orthogonal to ``_registry_id_epoch`` / ``_bridge_id_epoch``
            (Phase 3), which answer a different question: whether registry ids
            still NAME the same crossings. Adding a crossing changes the
            generation but preserves every existing id, so the epochs stay put;
            a renumbering recompute moves both. Staleness reads the generation,
            id validity reads the epochs.
        """
        state = (
            self._mutations,
            self._intersection_registry.generation,
            tuple(manifold.version for manifold in self._manifolds.values()),
        )
        if state != self._generation_state:
            self._generation_state = state
            self._generation += 1
        return self._generation

    def _bump_generation(self, *, invalidate_walks: bool = False) -> None:
        """Record a mutation only this object can see (see :attr:`generation`).

        Args:
            invalidate_walks: Also drop every registered bridge's memoised point
                array. Pass True for the operations that can lay new points
                INSIDE an already-cut arc -- growing a parent manifold, or
                re-cutting one against fresh crossings.
        """
        self._mutations += 1
        if invalidate_walks:
            for bridge in self.bridges:
                bridge.bump_version()

    @property
    def manifolds(self) -> dict[tuple[FixedPoint, Stability, int, int], BaseManifold]:
        """
        Every initialized manifold, keyed by
        ``(fixed_point, stability, orbit_index, branch_index)``.

        The key is not a second source of truth: it is exactly the manifold's own
        ``manifold_key``, and :meth:`register_manifold` is the only supported way in,
        precisely so the two can never drift apart.
        """
        return self._manifolds

    def register_manifold(
        self,
        key: tuple[FixedPoint, Stability, int, int],
        manifold: BaseManifold,
    ) -> None:
        """
        Store ``manifold`` under ``key``, asserting it agrees with its own key.

        Args:
            key: The manifold key to store under.
            manifold: The manifold, whose ``manifold_key`` must be ``key``.

        Raises:
            AssertionError: If the manifold's own key is not ``key``.
        """
        assert manifold.manifold_key == key, (
            f"manifold carries key {manifold.manifold_key} but is being registered "
            f"under {key}; the manifold's own key is the source of truth"
        )
        self._manifolds[key] = manifold
        self._bump_generation()

    def construct_fixed_point(
        self, initial_guess: Sequence | np.ndarray
    ) -> FixedPoint:
        """
        Constructs a fixed point for a given initial guess.
        Adds that fixed point to the class storage.

        Args:
            initial_guess: A ``(period, 2)`` array (or a single ``(2,)`` guess),
                one row per orbit point.

        Returns:
            FixedPoint: The constructed fixed point, already registered here and
                labelled ``A, B, ...`` in construction order. Its branch count is
                derived (``FixedPoint.num_branches``).

        Raises:
            ValueError: If the solve does not converge, if the orbit's minimal
                period is a proper divisor of the requested one (``period`` is
                the minimal period; a saddle with inversion is requested at its
                period, never at ``k_value``), if it is not a real saddle, or if
                its eigenvalues disagree in sign
                (:meth:`FixedPointSolver.construct_fixed_point`).

        Note:
            With a Jacobian, area preservation is checked once at the orbit's
            first point (:meth:`DynamicalSystem.check_area_preserving`) and a
            failure is logged as one WARNING: the construction still succeeds,
            but nothing downstream is validated for a dissipative map.
        """

        fixed_point = self._fp_solver.construct_fixed_point(initial_guess)
        orbit_point = fixed_point.coordinates[0]
        if self.dynamical_system.check_area_preserving(orbit_point) is False:
            logger.warning(
                "the map does not preserve area at %s (|det J| != 1); the "
                "topological layer assumes an area-preserving map",
                orbit_point,
            )
        position = len(self.fixed_points)
        fixed_point.label = position_letter(position)

        self.fixed_points.append(fixed_point)

        return fixed_point

    def orient_eigenvectors(
        self,
        fixed_point: FixedPoint,
        approx_dirs: Optional[dict[str, np.ndarray]] = None,
    ) -> None:
        """
        User supplies an approximate direction for the stable and unstable manifolds.
        If the eigenvectors align with that direction nothing happens, otherwise
        the eigenvectors are flipped.

        Input like:
            approx_dirs = {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}

        Note:
            The directions orient branch 0 at orbit point 0 and are propagated
            along the orbit. On an inversion point they choose which half of
            each manifold is branch 0; branch 1 runs the other way.
        """

        self._man_maker.orient_manifolds(fixed_point, approx_dirs)

    def initialize_manifold(
        self, fixed_point: FixedPoint, stability: Stability, num_branches: int = 1
    ) -> dict[tuple[int, int], BaseManifold]:
        """
        Build and register the fundamental segments of one manifold.

        Args:
            fixed_point: The fixed point to grow from.
            stability: Which manifold to initialize.
            num_branches: How many EIGENDIRECTIONS to seed, 1 or 2. This is a
                request, not a property of the fixed point: a simple saddle has
                two unstable directions and the caller may want both. On an
                inversion point the single chain of ``k_value`` pieces seeded
                on branch 0 already covers both branches
                (``FixedPoint.num_branches`` is 2 there), so 2 adds nothing.

        Returns:
            The initial segments, keyed by ``(orbit_index, branch_index)``.

        Raises:
            ValueError: If ``num_branches`` is neither 1 nor 2.
            PointBudgetExceeded: If refining the initial segments of a
                ``k_value > 1`` orbit (which maps them along the chain) would
                create more than :attr:`max_points_per_step` points. Nothing of
                this manifold is registered then.
        """

        initial_segments = self._man_maker.construct_kevin_way(
            fixed_point, stability, num_branches
        )

        for (orbit_index, branch_index), manifold in initial_segments.items():

            key = (fixed_point, stability, orbit_index, branch_index)
            self.register_manifold(key, manifold)

        return initial_segments

    def initialize_both_manifolds(
        self, fixed_point: FixedPoint, num_branches: int = 1
    ) -> tuple[
        dict[tuple[int, int], BaseManifold], dict[tuple[int, int], BaseManifold]
    ]:
        """
        Initialize the unstable and stable manifolds of one fixed point.

        Args:
            fixed_point: The fixed point to grow from.
            num_branches: Passed to :meth:`initialize_manifold` for both
                stabilities -- how many eigendirections to seed.

        Returns:
            ``(unstable_segments, stable_segments)``.

        Raises:
            PointBudgetExceeded: As :meth:`initialize_manifold`.
        """

        unstable_segments = self.initialize_manifold(
            fixed_point, "unstable", num_branches
        )
        stable_segments = self.initialize_manifold(fixed_point, "stable", num_branches)

        return (unstable_segments, stable_segments)

    def grow_n_times(
        self,
        fixed_point: FixedPoint,
        stability: Stability,
        num_iterations: int,
        branch_index: int = 0,
    ) -> None:
        """
        Grow one manifold a fixed number of iterations.

        Args:
            fixed_point (FixedPoint): The fixed point whose manifold is grown.
            stability (Stability): Which manifold to grow.
            num_iterations (int): How many times to apply the map.
            branch_index (int): Which eigenvector branch to grow. Defaults to 0.

        Raises:
            ValueError: If the manifold has not been initialized.
            PointBudgetExceeded: If one step would create more than
                :attr:`max_points_per_step` points. The tails and the generation
                are still updated, but the manifold is partly mapped and
                under-refined: the session is NOT resumable; rebuild it.

        Note:
            Growth refines the existing curve as well as extending it, so it
            bumps the generation and drops every memoised point array. It also
            re-finds the physical end of every branch of the grown stability, so
            stable growth undoes a stable trim (a resonance zone's included):
            the trim is a view, and the topology is cut at the pip orbit anyway.

        Note:
            A fixed count has no notion of escape: an arm leaving the trapped
            region multiplies its point count within a step or two, guarded only
            by the point budget. Before each step a WARNING hints when the last
            step's count times its growth factor exceeds the budget (a hint, not
            a refusal: past escape it underestimates by an order of magnitude),
            and with :attr:`escape_radius` set a step reaching it warns too.
            Prefer the stop-condition drivers (:meth:`grow_until_bounded`,
            :meth:`grow_until_turnaround`, :meth:`grow_until`). At the default
            :attr:`area_cutoff` (1e-7) a step makes ~10x the points it made at
            the old 1e-4, so a count tuned at 1e-4 can now run into the budget
            (k=10: 11 steps at 1e-7 is ~2e7 points).

            One iteration stretches a branch by ``per_step_beta("unstable")``
            (``lambda_u`` per map step), so the depth a tangle needs scales with
            ``1 / ln(lambda_u)``. A fixed count that works at one parameter can
            leave a weakly hyperbolic saddle (Hénon at small ``k``) with ZERO
            crossings: that is a depth problem, not a bug. Grow on a stop
            condition (``grow_until_*``, ``grow_until(predicate)``) instead.
        """

        self._require_manifold(fixed_point, stability, branch_index)

        machine = self._man_machine
        key = self._growth_key(fixed_point, stability, branch_index)
        try:
            for _ in range(num_iterations):
                last = self._growth_stats.get(key)
                budget = machine.max_points_per_step
                predicted = (
                    last.points_created * last.growth_factor
                    if last is not None and last.growth_factor is not None
                    else None
                )
                if predicted is not None and budget is not None and predicted > budget:
                    logger.warning(
                        "the next %s step of %s may create ~%.3g points, over the "
                        "%d point budget (last step %d, growth factor %.3g): the "
                        "arm is probably escaping; prefer a stop condition such "
                        "as grow_until_bounded",
                        stability,
                        fixed_point.label,
                        predicted,
                        budget,
                        last.points_created,
                        last.growth_factor,
                    )

                machine.new_grow_manifold(fixed_point, stability, branch_index)

                created = machine.points_created
                factor = (
                    created / last.points_created
                    if last is not None and last.points_created
                    else None
                )
                stats = GrowthStats(created, factor, machine.max_abs_coord)
                self._growth_stats[key] = stats
                if (
                    self.escape_radius is not None
                    and stats.max_abs_coord > self.escape_radius
                ):
                    logger.warning(
                        "the %s manifold of %s reached |coordinate| %.3g, beyond "
                        "escape_radius %.3g: the arm is escaping",
                        stability,
                        fixed_point.label,
                        stats.max_abs_coord,
                        self.escape_radius,
                    )
        finally:
            for (fp, stab, _orbit, _branch), manifold in self.manifolds.items():
                if fp is fixed_point and stab == stability:
                    manifold._find_tail()

            # Growth refines the existing curve as well as extending it, so points
            # can land inside an already-cut bridge: drop the memoised walks.
            self._bump_generation(invalidate_walks=True)

    def plot_tangle(
        self,
        fixed_point: FixedPoint,
        stability: Stability,
        **kwargs,
    ) -> None:
        """
        Plot every branch of one manifold of one fixed point, and the orbit.

        Args:
            fixed_point (FixedPoint): The fixed point to plot around.
            stability (Stability): Which manifold to plot.
            **kwargs: Forwarded to :meth:`BaseManifold.plot`.
        """

        for (fp, stab, _orbit_index, _branch_index), manifold in self.manifolds.items():

            if fp is fixed_point and stab == stability:

                manifold.plot(**kwargs)

        for period in range(fixed_point.period):
            plt.scatter(*fixed_point.coordinates[period], c="k", s=12)

    def grow_until_turnaround(
        self,
        fixed_point: FixedPoint,
        stability: Stability,
        max_iterations: Optional[int] = None,
        branch_index: int = 0,
    ) -> None:
        """
        Grows the manifold until a turnaround is detected or max_iterations is reached.

        Args:
            fixed_point (FixedPoint): The fixed point whose manifold is to be grown.
            stability (Stability): The stability type of the manifold ('stable' or 'unstable').
            max_iterations (Optional[int]): Maximum number of growth iterations
                to attempt before giving up. None (default) scales it with the
                stretch a turnaround needs (see the module Dev Notes), never
                below 10.
            branch_index (int, optional): Which eigenvector branch to grow.
                Defaults to 0.

        Raises:
            ValueError: If the manifold has not been initialized.
            MaxIterationsReached: A ``ValueError`` raised when the cap is reached
                before a turnaround is detected; it carries ``suggested_max_iterations``.
        """

        manifold = self._require_manifold(fixed_point, stability, branch_index)

        root = fixed_point.branch_points[0]
        first_point = manifold.first_node(branch_index)

        root_coord = root._coords
        first_point_coords = first_point._coords

        initial_direction = np.asarray(first_point_coords) - np.asarray(root_coord)

        if max_iterations is None:
            max_iterations = self._default_cap(fixed_point, stability, TURNAROUND_LOG_STRETCH)
        for _ in range(max_iterations):

            self.grow_n_times(
                fixed_point, stability, num_iterations=1, branch_index=branch_index
            )

            tail = self.manifolds.get((fixed_point, stability, 0, branch_index)).tail
            first_point = self.manifolds.get((fixed_point, stability, 0, branch_index))
            first_point = first_point.walk_back(None, tail)

            tail_coords = tail._coords
            first_point_coords = first_point._coords

            final_direction = tail_coords - first_point_coords

            dot_product = float(np.dot(initial_direction, final_direction))

            if dot_product < 0:
                return None

        raise self._cap_reached(
            fixed_point, stability, max_iterations, TURNAROUND_LOG_STRETCH,
            "a turnaround",
        )

    def grow_until_bounded(
        self,
        fixed_point: FixedPoint,
        stability: Stability,
        radius: float,
        max_iterations: Optional[int] = None,
        branch_index: int = 0,
    ) -> int:
        """
        Grow one step at a time until the manifold reaches ``|coordinate| > radius``.

        The escape stop: an arm leaving the trapped region is caught on the step
        it crosses ``radius``, before the steps whose point count explodes (the
        stress campaign found it the one stop condition that does). The
        driver to prefer at the default :attr:`area_cutoff` (1e-7), where a step
        past escape costs ~10x the points it costs at 1e-4.

        Args:
            fixed_point (FixedPoint): The fixed point whose manifold is grown.
            stability (Stability): Which manifold to grow.
            radius (float): The ``max|coordinate|`` that ends the growth.
            max_iterations (Optional[int]): Growth steps to attempt before giving
                up. None (default) scales it as :meth:`grow_until_turnaround`
                does.
            branch_index (int): Which eigenvector branch to grow. Defaults to 0;
                ignored on an inversion point, whose one chain covers both.

        Returns:
            The number of steps grown, the last being the one that crossed
            ``radius``.

        Raises:
            ValueError: If the manifold has not been initialized, or if the cap
                is reached with the manifold still inside ``radius``.
            PointBudgetExceeded: As :meth:`grow_n_times`.
        """
        self._require_manifold(fixed_point, stability, branch_index)

        if max_iterations is None:
            max_iterations = self._default_cap(fixed_point, stability, TURNAROUND_LOG_STRETCH)
        for step in range(1, max_iterations + 1):
            self.grow_n_times(fixed_point, stability, 1, branch_index)
            stats = self.growth_stats(fixed_point, stability, branch_index)
            if stats.max_abs_coord > radius:
                return step

        raise ValueError(
            f"Max iterations ({max_iterations}) reached with the {stability} "
            f"manifold still inside radius {radius}; choose a higher cap."
        )

    def grow_until_arclength(
        self,
        fixed_point: FixedPoint,
        stability: Stability,
        length: float,
        max_iterations: Optional[int] = None,
        branch_index: int = 0,
    ) -> None:
        """
        Grow the manifold until its polyline is at least ``length`` long.

        Args:
            fixed_point (FixedPoint): The fixed point whose manifold is grown.
            stability (Stability): Which manifold to grow.
            length (float): Target arc length of the branch's polyline, from the
                periodic point to the tail.
            max_iterations (Optional[int]): Growth iterations to attempt before
                giving up. None (default) scales it with ``length`` over the
                starting arc length (see the module Dev Notes), never below 10.
            branch_index (int): Which eigenvector branch to grow. Defaults to 0.

        Raises:
            ValueError: If the manifold has not been initialized.
            MaxIterationsReached: A ``ValueError`` raised when the cap is reached
                before the polyline reaches ``length``; it carries
                ``suggested_max_iterations``.
        """

        self._require_manifold(fixed_point, stability, branch_index)

        def arc_length() -> float:
            points = self.manifolds[(fixed_point, stability, 0, branch_index)].get_point_array()
            return float(np.linalg.norm(np.diff(points, axis=0), axis=1).sum())

        log_stretch = math.log(max(length / arc_length(), 1.0))
        if max_iterations is None:
            max_iterations = self._default_cap(fixed_point, stability, log_stretch)
        for _ in range(max_iterations):
            if arc_length() >= length:
                return None
            self.grow_n_times(
                fixed_point, stability, num_iterations=1, branch_index=branch_index
            )

        if arc_length() >= length:
            return None

        raise self._cap_reached(
            fixed_point, stability, max_iterations, log_stretch,
            f"the polyline reaching arc length {length:g}",
        )

    def grow_until_intersection(
        self,
        fixed_point: FixedPoint,
        stability: Stability,
        max_iterations: Optional[int] = None,
        branch_index: int = 0,
    ) -> None:
        """
        Grow the manifold until a new crossing appears.

        Args:
            fixed_point (FixedPoint): The fixed point whose manifold is grown.
            stability (Stability): Which manifold to grow.
            max_iterations (Optional[int]): Growth iterations to attempt before
                giving up. None (default) scales it as
                :meth:`grow_until_turnaround` does.
            branch_index (int): Which eigenvector branch to grow. Defaults to 0.

        Raises:
            ValueError: If the manifold has not been initialized.
            MaxIterationsReached: A ``ValueError`` raised when the cap is reached
                before a new crossing appears; it carries ``suggested_max_iterations``.
        """

        self._require_manifold(fixed_point, stability, branch_index)

        self.compute_intersections(fixed_point, infer_iterates=False, diagnose=False)

        num_initial_intersections = len(self._intersection_registry)
        logger.info("current intersections %d", num_initial_intersections)

        if max_iterations is None:
            max_iterations = self._default_cap(fixed_point, stability, TURNAROUND_LOG_STRETCH)
        for _ in range(max_iterations):

            self.grow_n_times(
                fixed_point, stability, num_iterations=1, branch_index=branch_index
            )

            self.compute_intersections(
                fixed_point, infer_iterates=False, diagnose=False
            )
            num_current_intersections = len(self._intersection_registry)

            if num_current_intersections > num_initial_intersections:
                return None

        raise self._cap_reached(
            fixed_point, stability, max_iterations, TURNAROUND_LOG_STRETCH,
            "a new crossing",
        )

    # ── predicate-driven growth ──────────────────────────────────────────────

    def grow_until(
        self,
        fixed_point: FixedPoint,
        predicate: Callable[["TangleWorkbench"], bool],
        *,
        grow: Sequence[Stability] = ("unstable", "stable"),
        max_iterations: int = 10,
        branch_index: int = 0,
    ) -> int:
        """
        Grow a tangle one iteration at a time until the caller's predicate holds.

        Each round grows ONE iteration on every stability named in ``grow``, then
        recomputes the crossings of every fixed point that has manifolds --
        together, in one :class:`~.Tangle.Tangle`, so heteroclinic crossings are
        seen -- with ``preserve_ids=True``, and re-infers the iterate table. The
        predicate is then asked whether the tangle is big enough.

        ``preserve_ids=True`` is what makes the loop usable at all: a predicate is
        written in terms of crossing ids the caller held BEFORE the loop started,
        and a renumbering recompute would hand those ids to unrelated crossings
        halfway through.

        The predicate is checked BEFORE the first growth, so a request that is
        already satisfied costs nothing.

        Args:
            fixed_point: The fixed point whose manifolds are grown. (Every fixed
                point's crossings are recomputed; only this one's manifolds move.)
            predicate: Called with this workbench after each recompute (and once
                before any growth). Return True to stop.
            grow: Which stabilities to advance each round. Defaults to both.
            max_iterations: Growth rounds to attempt before giving up. Defaults
                to 10, matching :meth:`grow_until_turnaround`.
            branch_index: Which eigenvector branch to grow. Defaults to 0. On an
                INVERSION point every branch index names the one chain, which
                covers both branches.

        Returns:
            The number of growth rounds performed: 0 when the predicate was
            already true.

        Raises:
            ValueError: If ``grow`` is empty, if a named manifold has not been
                initialized, or if the cap is reached with the predicate still
                false.

        Note:
            Only ONE branch grows. On a simple saddle whose two eigendirections
            were both seeded (``num_branches=2``) the other branch stays where it
            is, so a predicate about the whole tangle can sit false forever; drive
            the branches in separate calls, or grow the other one by hand first.

        Note:
            If the workbench already carries bridges they are recut each round
            (:meth:`rebuild_bridges`): growth inserts new crossings between the
            endpoints of existing bridges, so the old cut stops describing the
            curve. A blast's image curves are recomputed and recut with the
            branches, so the blast frontier, its partial pieces and the
            ``iterated`` flags survive the round.
        """
        stabilities = tuple(grow)
        if not stabilities:
            raise ValueError("grow must name at least one stability to advance.")

        for stability in stabilities:
            if self.manifolds.get((fixed_point, stability, 0, branch_index)) is None:
                raise ValueError(
                    f"Manifold for fixed point {fixed_point} with stability "
                    f"{stability} and branch_index {branch_index} has not been "
                    "initialized. Please run initialize_manifold first."
                )

        if predicate(self):
            return 0

        for round_index in range(1, max_iterations + 1):

            for stability in stabilities:
                self.grow_n_times(
                    fixed_point, stability, num_iterations=1, branch_index=branch_index
                )

            self.compute_intersections(
                self._fixed_points_with_manifolds(),
                infer_iterates=True,
                preserve_ids=True,
                diagnose=False,
            )
            if self._bridges or self._partial_bridges:
                self.rebuild_bridges()

            if predicate(self):
                self.report_resolution(self._fixed_points_with_manifolds())
                return round_index

        raise ValueError(
            f"Max iterations ({max_iterations}) reached without the growth "
            "predicate becoming true; choose a higher cap or a different "
            "stopping condition."
        )

    def _fixed_points_with_manifolds(self) -> list[FixedPoint]:
        """The fixed points that own at least one registered manifold, in insertion order."""
        selected: list[FixedPoint] = []
        seen: set[FixedPoint] = set()
        for fp, _stability, _orbit_index, _branch_index in self._manifolds:
            if fp not in seen:
                seen.add(fp)
                selected.append(fp)
        return selected

    def _frozen_ids(
        self,
        ids: "Union[str, Iterable[int]]",
        *,
        fixed_point: Optional[FixedPoint] = None,
        allow_all: bool = False,
    ) -> list[int]:
        """Resolve an id selector against the registry as it stands right now.

        The list is deliberately materialised here rather than re-read inside a
        predicate: the crossings born during a growth loop must not move the
        goalposts. ``"all"`` means every non-anchor crossing of N OF
        ``fixed_point`` -- scoping it to the whole registry would put another
        tangle's crossings in a set that only this fixed point's growth can ever
        close, and an unresolved detection is never linked, so the driver would
        run to its cap on it.
        """
        registry = self._intersection_registry
        if isinstance(ids, str):
            if not allow_all or ids != "all":
                raise ValueError(
                    f"Unknown id selector {ids!r}; pass an iterable of registry ids"
                    + (' or "all".' if allow_all else ".")
                )
            assert fixed_point is not None, '"all" needs a fixed point to scope to'
            resolved = registry.resolved_ids(self.resolution_ratio_threshold)
            selected = sorted(
                iid
                for iid, ix in registry
                if iid in resolved
                and not ix.is_synthetic
                and any(fp is fixed_point for fp in ix.fixed_points)
            )
        else:
            selected = [int(iid) for iid in ids]
            missing = [iid for iid in selected if iid not in registry]
            if missing:
                raise ValueError(f"Crossing ids {missing} are not in the registry.")

        if not selected:
            raise ValueError(
                "A growth driver needs at least one crossing id: an empty set is "
                "satisfied vacuously and the call would do nothing. Compute the "
                "intersections first, or name the crossings explicitly."
            )
        return selected

    def grow_until_iterates_closed(
        self,
        fixed_point: FixedPoint,
        ids: "Union[str, Iterable[int]]" = "all",
        direction: Literal["forward", "backward"] = "forward",
        *,
        max_iterations: int = 10,
        branch_index: int = 0,
    ) -> int:
        """
        Grow until every crossing in a set has its next iterate detected.

        The image of a crossing one map step forward sits at a LARGER unstable
        canonical distance (stretched by the unstable eigenvalue) and a smaller
        stable one, so detecting it needs the unstable manifold grown; a backward
        image is the mirror image of that statement and needs the stable manifold.
        The wrapper picks the stability from ``direction`` -- growing the other one
        would only cost time.

        The id set is frozen at call time. Growth manufactures new crossings every
        round, and re-reading "all" inside the predicate would add a fresh unclosed
        crossing for every one just closed, so the loop would run to the cap on
        every tangle.

        Args:
            fixed_point: The fixed point whose manifold is grown.
            ids: Registry ids to close, or ``"all"`` (the default) for every
                non-anchor crossing OF ``fixed_point`` in the registry at call
                time. Only this fixed point's manifolds grow, so a set scoped to
                the whole registry would hold another tangle's crossings hostage
                and run to the cap in any multi-saddle session. Anchors are left
                out because their iterates are declared, not detected (an anchor's
                forward image is the anchor of the NEXT point of the periodic
                orbit -- itself only when the period is 1).
            direction: ``"forward"`` for the ``n=+1`` entry, ``"backward"`` for
                ``n=-1``.
            max_iterations: Growth rounds to attempt. Defaults to 10.
            branch_index: Which eigenvector branch to grow. Defaults to 0.

        Returns:
            The number of growth rounds performed (0 if already closed).

        Raises:
            ValueError: If ``direction`` is not one of the two directions, if an
                id is not in the registry, if the id set is empty, or if the cap is
                reached with some iterate still undetected.
        """
        if direction not in ("forward", "backward"):
            raise ValueError(
                f"direction must be 'forward' or 'backward', not {direction!r}."
            )

        step = 1 if direction == "forward" else -1
        stability: Stability = "unstable" if direction == "forward" else "stable"
        targets = self._frozen_ids(ids, fixed_point=fixed_point, allow_all=True)

        def closed(workbench: "TangleWorkbench") -> bool:
            table = workbench.intersection_registry.iterate_table
            return all((iid, step) in table for iid in targets)

        return self.grow_until(
            fixed_point,
            closed,
            grow=(stability,),
            max_iterations=max_iterations,
            branch_index=branch_index,
        )

    def grow_until_faces_closed(
        self,
        fixed_point: FixedPoint,
        ids: "Iterable[int]",
        *,
        grow: Sequence[Stability] = ("unstable", "stable"),
        max_iterations: int = 10,
        branch_index: int = 0,
    ) -> int:
        """
        Grow until every face touching a set of crossings is closed.

        A face of the arrangement is open when its boundary runs off the end of a
        computed manifold, which means the region there is not yet determined: the
        answer to "what lies on this side of that crossing" is still "grow more".
        This driver grows until that answer exists for every crossing in ``ids``.

        The arrangement is rebuilt each round over the same fixed points the
        recompute used, rather than taken from a session cache -- every round moves
        the workbench generation, so a cache would miss every time anyway. It is
        deliberately not a single-fixed-point snapshot: that holds no heteroclinic
        crossing and none of the partner's arcs, turning real faces into open ones.

        Args:
            fixed_point: The fixed point whose manifolds are grown.
            ids: Registry ids whose incident faces must all close. Required and
                non-empty: an empty set closes vacuously and would make the call a
                silent no-op.
            grow: Which stabilities to advance. Defaults to both, because a face
                is bounded by arcs of both.
            max_iterations: Growth rounds to attempt. Defaults to 10.
            branch_index: Which eigenvector branch to grow. Defaults to 0.

        Returns:
            The number of growth rounds performed (0 if already closed).

        Raises:
            ValueError: If ``ids`` is empty, if an id is not in the registry, or
                if the cap is reached with some incident face still open.

        Note:
            In practice this is a predicate CHECK with a cap rather than a driver
            that converges. No growth-driven open-to-closed transition has been
            observed on the k=10 fixture: the real corners of the unbounded face
            stay on it however far the branches are grown (growing a branch pushes
            its tip past them but they remain on the outer envelope), and the
            crossings born by a round are born already interior. Expect either an
            immediate 0 or the cap.

        Note:
            A crossing on the outer envelope never closes, the inversion anchor
            included (one degree-4 node of its four branches): that is a
            fact about where it sits, not a missing growth step. A
            non-inversion anchor seeded on one eigendirection keeps virtual
            ``u-`` / ``s-`` stubs (the other half of each manifold is not
            computed), so its faces through them are open as well.
        """
        from ..topology.Trellis import Trellis

        targets = self._frozen_ids(ids, fixed_point=fixed_point, allow_all=False)

        def closed(workbench: "TangleWorkbench") -> bool:
            selection = workbench._fixed_points_with_manifolds()
            arrangement = Trellis.from_workbench(workbench, selection).arrangement
            corners: set[int] = set()
            open_corners: set[int] = set()
            for face in arrangement.faces:
                # Virtual tail nodes carry negative ids and are not crossings; a
                # target is always a real registry id, so they only add noise.
                for corner in (c for c in face.corners if c >= 0):
                    corners.add(corner)
                    if not face.is_closed:
                        open_corners.add(corner)
            # A crossing that is not a node of the arrangement at all has no
            # determined surroundings either, so it does not count as closed.
            return all(iid in corners and iid not in open_corners for iid in targets)

        return self.grow_until(
            fixed_point,
            closed,
            grow=grow,
            max_iterations=max_iterations,
            branch_index=branch_index,
        )

    def compute_intersections(
        self,
        fixed_points: FixedPoint | Iterable[FixedPoint],
        *,
        reset: bool = True,
        infer_iterates: bool = True,
        preserve_ids: bool = False,
        diagnose: bool = True,
    ) -> IntersectionRegistry:
        """
        Compute intersections among the manifolds of one or more fixed points.

        All manifolds of every supplied fixed point are indexed into the SAME Tangle
        before crossings are resolved, so homoclinic crossings (within one fixed
        point) and heteroclinic crossings (between two fixed points) are detected
        together.

        Args:
            fixed_points: A single FixedPoint or an iterable of FixedPoints whose
                manifolds should be co-indexed and intersected.
            reset: If True (default) the Tangle and registry are cleared first. Pass
                False to accumulate further manifolds into an existing computation.
            infer_iterates: If True (default) fill the iterate table from the freshly
                computed crossings via :meth:`infer_iterates` (the M^1 forward iterate
                of every intersection, by the exact bracket rule). Pass False to
                skip it — used by the growth loops that call this many times and do not
                need the table.
            preserve_ids: If True (and ``reset`` is True), re-align the rebuilt
                registry against the one being replaced so any crossing that reappears
                keeps its previous id (see :meth:`IntersectionRegistry.reindex_from`).
                Used by the resonance-zone recompute so a strong pip chosen as id N is
                still id N after the stable manifolds are trimmed. Defaults to False —
                the growth loops renumber freely.
            diagnose: If True (default) finish with :meth:`report_resolution`:
                one WARNING for new unresolved crossings, one for new
                same-stability polyline crossings in the tangle region. The
                growth loops pass False and report once at the end.

        Returns:
            The registry the workbench now holds. With ``reset`` (the default)
            this is a NEW object: a reference to :attr:`intersection_registry`
            taken before the call is stale (its crossings and ids belong to the
            old manifolds), so re-fetch the registry after every call or use
            this return value. Coordinates are ``[ix.coords for _, ix in
            registry]``.
        """
        fixed_points = (
            [fixed_points] if isinstance(fixed_points, FixedPoint) else list(fixed_points)
        )

        old_registry = self._intersection_registry
        if reset:
            if len(old_registry) > 0:
                logger.debug(
                    "compute_intersections replaces a registry of %d crossings; "
                    "references to the old one are stale",
                    len(old_registry),
                )
            self.Tangle.clear_all()
            self._intersection_registry = IntersectionRegistry()
            # The registry object itself changed; its own generation restarts at
            # zero, so the swap has to be recorded here or the two could cancel.
            self._bump_generation()

        # gather every manifold first so a fresh Tangle can bulk-load the whole
        # segment set into the rtree in one pass
        manifolds_to_index = []
        for fp in fixed_points:
            manifolds_to_index.extend(self._iter_manifolds(fp, "unstable"))
            manifolds_to_index.extend(self._iter_manifolds(fp, "stable"))
        self.Tangle.add_manifolds(manifolds_to_index)
        # A blast image is computed manifold too: re-indexed query-only,
        # as the blast indexed it. On a FRESH index a curve that adds no segment of
        # its own lies on a grown branch now (growth merged it in) and is dropped;
        # an accumulating call cannot tell (the bridges cut out of the curve still
        # own its edges), so it keeps every curve.
        kept = []
        for curve in self._image_curves:
            if (all(curve.fixed_point is not fp for fp in fixed_points)
                    or self.Tangle.add_manifold(curve, index_segments=False)
                    or not reset):
                kept.append(curve)
            else:
                self.Tangle.release_manifold(curve)
        self._image_curves = kept

        # Anchors FIRST: every periodic point is a crossing of its own manifolds,
        # and declaring it beats detecting it (see _register_anchors).
        self._register_anchors(fixed_points)

        # Registration order is the order along each unstable branch, never
        # the iteration order of a set of segment ids: ids are reproducible.
        position = {id(fp): i for i, fp in enumerate(self.fixed_points)}

        def along(ix: Intersection) -> tuple:
            (u_fp, _, u_orbit, u_branch), (s_fp, _, s_orbit, s_branch) = (
                ix.manifold_a_key, ix.manifold_b_key
            )
            return (position[id(u_fp)], u_orbit, u_branch, ix.unstable_cdist,
                    position[id(s_fp)], s_orbit, s_branch, ix.stable_cdist)

        for intersection in sorted(self.Tangle.resolve_crossings(), key=along):
            intersection.id = self._intersection_registry.add(intersection)

        if reset and preserve_ids and len(old_registry) > 0:
            self._intersection_registry.reindex_from(old_registry)
        elif reset:
            # The registry was rebuilt from zero: every id held by a caller (a strong
            # pip, a BridgeId) now names a different crossing, if anything at all.
            self._registry_id_epoch += 1

        if infer_iterates:
            self.infer_iterates()

        if diagnose:
            self.report_resolution(fixed_points)

        return self._intersection_registry

    # ── resolution diagnostics (see diagnostics) ─────────────────────────────

    def check_resolution(
        self, fixed_points: "Optional[FixedPoint | Iterable[FixedPoint]]" = None
    ) -> ResolutionReport:
        """
        What the polylines fail to resolve, on demand and without logging.

        Args:
            fixed_points: Whose crossings and curves to check; None checks every
                fixed point.

        Returns:
            The registered detections that are not crossings of N (below
            :attr:`resolution_ratio_threshold`, or in a run) and every
            same-stability polyline crossing
            (:func:`~.diagnostics.same_stability_crossings`). Neither feeds any
            reliability verdict; the unresolved detections are already out of N
            by definition, and the report only names them.
        """
        selected = self._selection(fixed_points)
        return ResolutionReport(
            threshold=self.resolution_ratio_threshold,
            unresolved_crossings=self._unresolved_crossings(selected),
            same_stability_crossings=tuple(same_stability_crossings(self, selected)),
        )

    def report_resolution(
        self,
        fixed_points: "Optional[FixedPoint | Iterable[FixedPoint]]" = None,
        *,
        new_bridges: Optional[Iterable[Bridge]] = None,
    ) -> ResolutionReport:
        """
        :meth:`check_resolution`, logging what no earlier report has said.

        One WARNING when unresolved detections appear that were not reported
        before, naming every unresolved one (id, coordinates, ratio, clause,
        and for a run its partner and separation) and ending with the remedy
        "refine: lower `area_cutoff` or grow the stretch" and the issue id; one
        WARNING when same-stability polyline crossings
        appear in the tangle region, one INFO when they appear only outside it
        (an escaping arm). Each finding is reported once, so the zone
        recomputes and the next blast stay quiet about what is already known.
        Everything goes to this module's logger.

        Args:
            fixed_points: As in :meth:`check_resolution`.
            new_bridges: A blast's new image bridges: the same-stability
                part is then only COUNTED, over the hits these bridges take part
                in (:func:`~.diagnostics.count_same_stability_crossings`), and
                nothing is remembered -- every such hit is new. A blast adds no
                stable curve and every hit it can add involves one of its images,
                so the full listing would re-find the earlier blasts' hits.

        Returns:
            The report; with ``new_bridges`` its ``same_stability_crossings``
            is empty and ``same_stability_count`` holds the count.
        """
        if new_bridges is None:
            report = self.check_resolution(fixed_points)
            self._report_artifacts(report.same_stability_crossings)
        else:
            count = count_same_stability_crossings(self, new_bridges)
            report = ResolutionReport(
                threshold=self.resolution_ratio_threshold,
                unresolved_crossings=self._unresolved_crossings(self._selection(fixed_points)),
                same_stability_crossings=(),
                same_stability_count=count,
            )
            self._report_count(count)
        self._report_unresolved(report.unresolved_crossings, report.threshold)
        return report

    def _selection(
        self, fixed_points: "Optional[FixedPoint | Iterable[FixedPoint]]"
    ) -> list[FixedPoint]:
        """The fixed points a resolution check covers (None: all of them)."""
        if fixed_points is None:
            return list(self.fixed_points)
        if isinstance(fixed_points, FixedPoint):
            return [fixed_points]
        return list(fixed_points)

    def _unresolved_crossings(self, selected: list[FixedPoint]) -> tuple[int, ...]:
        """The selection's detections that are not crossings of N."""
        registry = self._intersection_registry
        return tuple(
            iid
            for iid in registry.unresolved_crossings(self.resolution_ratio_threshold)
            if any(fp is chosen for fp in registry[iid].fixed_points for chosen in selected)
        )

    def pending_resolution_flags(self) -> list[int]:
        """The unresolved detections no report has logged yet, in unstable cdist order."""
        reported = self._reported_flag_ids()
        return [
            iid
            for iid in self._intersection_registry.unresolved_crossings(
                self.resolution_ratio_threshold
            )
            if iid not in reported
        ]

    def _reported_flag_ids(self) -> frozenset[int]:
        """The flagged ids already reported, while the registry still numbers them so."""
        epoch, ids = self._reported_flags
        return ids if epoch == self._registry_id_epoch else frozenset()

    def _report_unresolved(
        self, unresolved: tuple[int, ...], threshold: Optional[float]
    ) -> None:
        """One WARNING naming every unresolved detection, if any of them is new."""
        reported = self._reported_flag_ids()
        self._reported_flags = (self._registry_id_epoch, reported | frozenset(unresolved))
        if all(iid in reported for iid in unresolved):
            return
        registry = self._intersection_registry
        lines = []
        for iid in unresolved:
            ix = registry[iid]
            clauses = []
            if threshold is not None and ix.resolution_ratio is not None and ix.resolution_ratio < threshold:
                clauses.append(f"ratio below {threshold:.3g} (clause 1)")
            partner = registry.run_partner(iid)
            if partner is not None:
                clauses.append(
                    f"in a run with {partner}, the nearest detection "
                    f"{registry.run_separation(iid):.3g} away along the polyline (clause 2)"
                )
            lines.append(
                f"  {iid} at ({ix.coords[0]:.12g}, {ix.coords[1]:.12g}): ratio "
                f"{ix.resolution_ratio:.4g}, " + "; ".join(clauses)
            )
        logger.warning(
            "%d of %d detections are not crossings of N: the polylines do not resolve "
            "them at area_cutoff %.3g, so they are registered and read by no path "
            "(an unresolved one lies inside a bridge and is never linked):\n%s\n"
            "refine: lower `area_cutoff` or grow the stretch; NUM-3/4",
            len(unresolved),
            sum(not ix.is_synthetic for _iid, ix in registry),
            self._man_machine.area_cutoff,
            "\n".join(lines),
        )

    def _report_artifacts(self, hits: tuple[SameStabilityCrossing, ...]) -> None:
        """One WARNING for new same-stability hits in the tangle region, one INFO outside."""

        def key(hit: SameStabilityCrossing) -> tuple:
            return (hit.stability, hit.curve_a, hit.curve_b, hit.point)

        inside = [hit for hit in hits if hit.in_tangle_region]
        outside = [hit for hit in hits if not hit.in_tangle_region]
        if any(key(hit) not in self._reported_artifacts for hit in inside):
            pairs = sorted({f"{hit.curve_a} x {hit.curve_b}" for hit in inside})
            logger.warning(
                "%d same-stability polyline crossing(s) in the tangle region (%s): "
                "two curves of one stability never cross, so these are polygon "
                "artifacts of under-resolved polylines; an area_cutoff below %.3g "
                "is needed to refine them at all (current %.3g)",
                len(inside),
                "; ".join(pairs[:5]) + (" ..." if len(pairs) > 5 else ""),
                min(hit.suggested_area_cutoff for hit in inside),
                self._man_machine.area_cutoff,
            )
        if any(key(hit) not in self._reported_artifacts for hit in outside):
            logger.info(
                "%d same-stability polyline crossing(s) outside the tangle region: "
                "an escaping arm is under-resolved",
                len(outside),
            )
        self._reported_artifacts.update(key(hit) for hit in hits)

    def _report_count(self, count: SameStabilityCount) -> None:
        """One WARNING for a blast's same-stability hits in the tangle region, one INFO outside."""
        if count.in_region:
            pairs = count.curve_pairs
            logger.warning(
                "%d same-stability polyline crossing(s) of the new image bridges in "
                "the tangle region (%s): two curves of one stability never cross, so "
                "these are polygon artifacts of under-resolved polylines; an "
                "area_cutoff below %.3g is needed to refine them at all (current %.3g)",
                count.in_region,
                "; ".join(pairs[:5]) + (" ..." if len(pairs) > 5 else ""),
                count.suggested_area_cutoff,
                self._man_machine.area_cutoff,
            )
        if count.outside:
            logger.info(
                "%d same-stability polyline crossing(s) of the new image bridges "
                "outside the tangle region: an escaping arm is under-resolved",
                count.outside,
            )

    def _register_anchors(self, fixed_points: Iterable[FixedPoint]) -> list[int]:
        """
        Declare every periodic point as a crossing of its own manifolds.

        A periodic point sits at canonical distance 0 on every branch attached to
        it, so it IS a crossing -- the ``(0, 0)`` anchor the topological layer
        treats as the innermost node of every one of those branches. It is ONE
        crossing: at an inversion point the two stable and two unstable
        branches all meet there, and the arrangement draws it as one degree-4
        node.

        Anchors are registered BEFORE :meth:`Tangle.resolve_crossings`, so:

        * they take the lowest registry ids, deterministically ordered by fixed
          point and orbit index rather than by whatever order the rtree
          returned its candidate pairs;
        * the root segments detected crossing there collide with the anchor in
          :meth:`IntersectionRegistry.add` (exactly ``(0, 0)`` at the same
          periodic point) and are dropped in its favour -- one crossing, and the
          id belongs to the declared one, which ``preserve_ids`` relies on.

        The anchor is stored under its HOME pair ``(u(o, 0), s(o, 0))`` with
        that pair's crossing sign; branch 1 of either stability leaves it the
        opposite way. It carries no curve of its own:
        :meth:`Tangle.create_bridges` cuts every unstable branch rooted at the
        point at that branch's root.

        Args:
            fixed_points: The fixed points whose periodic points to anchor.

        Returns:
            The registry ids of the anchors, in registration order. A point
            whose branch-0 manifolds are missing or have fewer than two nodes
            gets none.
        """
        registry = self._intersection_registry
        anchor_ids: list[int] = []
        for fp in fixed_points:
            for orbit_index in range(fp.period):
                u_key = (fp, "unstable", orbit_index, 0)
                s_key = (fp, "stable", orbit_index, 0)
                if u_key not in self.manifolds or s_key not in self.manifolds:
                    continue
                u_nodes = self.manifolds[u_key].get_point_array(return_nodes=True)
                s_nodes = self.manifolds[s_key].get_point_array(return_nodes=True)
                if len(u_nodes) < 2 or len(s_nodes) < 2:
                    continue
                coords = u_nodes[0].get_point()
                anchor = Intersection.synthetic(
                    coords=(float(coords[0]), float(coords[1])),
                    unstable_cdist=0.0,
                    stable_cdist=0.0,
                    label="anchor",
                    manifold_a_key=u_key,
                    manifold_b_key=s_key,
                    crossing_sign=self._anchor_crossing_sign(u_nodes, s_nodes),
                )
                anchor_ids.append(registry.add(anchor))
        return anchor_ids

    @staticmethod
    def _anchor_crossing_sign(u_nodes: list, s_nodes: list) -> int:
        """
        Handedness of an anchor, from the two oriented eigendirections.

        The eigendirections are read off the fundamental segments rather than off
        ``FixedPoint``: the initializer lays each fundamental segment along the
        oriented eigenvector, so the first node-to-node step IS that eigenvector
        -- the one of branch 0, the anchor's home pair. Branch 1 of either
        stability runs the other way, so each branch index other than 0 flips
        the sign of a pair (``cross(-u, s) = -cross(u, s)``).

        Args:
            u_nodes: The unstable branch's nodes, root first.
            s_nodes: The stable branch's nodes, root first.

        Returns:
            The sign of ``cross(unstable eigendirection, stable eigendirection)``:
            ``+1``, ``-1``, or ``0`` if the two are exactly parallel (impossible
            for a real saddle, whose eigenvectors are independent).
        """
        u_dir = u_nodes[1].get_point() - u_nodes[0].get_point()
        s_dir = s_nodes[1].get_point() - s_nodes[0].get_point()
        return int(np.sign(float(u_dir[0] * s_dir[1] - u_dir[1] * s_dir[0])))

    def plot_intersections(
        self,
        fp: Optional[FixedPoint] = None,
        ax: Optional[plt.Axes] = None,
        show_ids: bool = False,
        id_fontsize: int = 8,
        **scatter_kwargs,
    ) -> None:
        """
        Scatter-plot computed intersections, optionally restricted to one fixed point.

        In a nested / multi-tangle session the registry holds the crossings of
        every fixed point at once. Pass ``fp`` to plot only that tangle's
        intersections; pass None to plot all of them. If nothing has been computed
        yet, the intersections for ``fp`` are computed first.

        Args:
            fp: Fixed point whose intersections to plot. None plots every computed
                intersection.
            ax: Optional matplotlib Axes. Defaults to the current axes (plt).
            show_ids: If True, label each intersection with its registry id
                (the same ids used by the Trellis / strong-pip API).
            id_fontsize: Font size of the id labels. Adjust to taste.
            **scatter_kwargs: Forwarded to scatter (defaults: s=7, zorder=10,
                color="k").
        """
        if len(self._intersection_registry) == 0 and fp is not None:
            self.compute_intersections(fp, infer_iterates=False, diagnose=False)

        # Pull coords (and ids, for labelling) from the registry so we can filter
        # by fixed point; an intersection belongs to fp if fp is among the fixed
        # points of its two manifold sides.
        items = [
            (iid, ix)
            for iid, ix in self._intersection_registry
            if fp is None or fp in ix.fixed_points
        ]
        if not items:
            logger.info("No intersections to plot.")
            return
        pts = np.array([ix.coords for _iid, ix in items])

        # sensible defaults; caller can override with kwargs
        scatter_kwargs.setdefault("s", 7)
        scatter_kwargs.setdefault("zorder", 10)
        scatter_kwargs.setdefault("color", "k")
        target = ax if ax is not None else plt
        target.scatter(pts[:, 0], pts[:, 1], **scatter_kwargs)

        if show_ids:
            for iid, intersection in items:
                x, y = intersection.coords
                target.annotate(
                    str(iid),
                    (x, y),
                    textcoords="offset points",
                    xytext=(4, 4),
                    fontsize=id_fontsize,
                    zorder=scatter_kwargs["zorder"] + 1,
                )

    def _register_bridge(self, bridge: Bridge) -> Bridge:
        """
        Store one freshly cut bridge, or return the copy already held.

        A bridge is uniquely identified by the two crossings it connects
        (:data:`~.Bridge.BridgeId`), so a freshly cut piece whose id is already
        registered IS that bridge -- the stored object is returned in its place and
        the duplicate gives up the segments the cut claimed for it. Partial pieces
        have no id and are simply appended.

        Every bridge enters here, so this is where its invariants are asserted:
        a bridge is an arc of ITS unstable branch between its two
        crossings. A non-anchor end sits on the bridge's own branch, inside its node
        span up to the registry's dedup tolerance; an anchor end comes first and is
        the bridge's own periodic point (one anchor on every branch there); and
        the id runs in unstable-cdist order.
        """
        bid = bridge.id
        if bid is None:
            self._partial_bridges.append(bridge)
            return bridge

        existing = self._bridges.get(bid)
        if existing is not None:
            if existing is not bridge:
                self.Tangle.release_manifold(bridge)
            return existing

        lo, hi = sorted(float(point.get_cdist("unstable")) for point in (bridge.root, bridge.tail))
        slack = DEDUP_RTOL * max(abs(hi), 1.0)
        fixed_point, _, orbit_index, _ = bridge.manifold_key
        first, second = (self._intersection_registry[end] for end in bid)
        assert first.is_anchor or first.unstable_cdist <= second.unstable_cdist, bid
        for end in bid:
            crossing = self._intersection_registry[end]
            if crossing.is_anchor:
                assert end == bid[0], bid
                assert crossing.manifold_a_key[0] is fixed_point and crossing.manifold_a_key[2] == orbit_index, bid
            else:
                assert crossing.manifold_a_key == bridge.manifold_key, (bid, end, crossing.manifold_a_key[1:])
                assert lo - slack <= crossing.unstable_cdist <= hi + slack, (bid, end, lo, crossing.unstable_cdist, hi)

        self._bridges[bid] = bridge
        for endpoint in bid:
            self._bridges_at.setdefault(endpoint, []).append(bid)
        return bridge

    def create_bridges(self, fixed_point: Optional[FixedPoint] = None) -> list[Bridge]:
        """
        Cut indexed unstable manifolds into bridges.

        The cut is at the crossings of N: an unresolved detection lies inside
        a bridge. A blast's image curves are indexed unstable manifold too:
        one crossed twice or more is cut like a branch, one crossed fewer
        than twice is held whole as a partial piece.

        Args:
            fixed_point: If given, only build bridges for that fixed point's
                unstable manifolds. If None, build bridges for every indexed
                unstable manifold at once. Each bridge keeps its own fixed_point
                linkage either way.

        Returns:
            The bridges cut, each the single registered copy for its
            :data:`~.Bridge.BridgeId`.
        """
        registry = self._intersection_registry
        resolved = registry.resolved_ids(self.resolution_ratio_threshold)
        crossings = [ix for iid, ix in registry if iid in resolved]
        bridges = self.Tangle.create_bridges(
            crossings,
            fixed_point=fixed_point,
            cdist_tol=self._intersection_registry.cdist_tol,
        )
        # Stamp the id epoch only when starting from nothing. Cutting ON TOP of an
        # existing set after a renumbering recompute leaves stale ids mixed in with
        # fresh ones, and keeping the older epoch is what makes the next
        # rebuild_bridges refuse to carry metadata across that mixture.
        if not self._bridges and not self._partial_bridges:
            self._bridge_id_epoch = self._registry_id_epoch
        # Cutting inserts a separator point at every crossing, which can land
        # inside an arc an earlier cut already claimed: drop the memoised walks.
        self._bump_generation(invalidate_walks=True)
        cut = [self._register_bridge(bridge) for bridge in bridges]
        # An image curve of an indexed fixed point crossed fewer than twice is one
        # partial piece, as at blast time.
        crossings_on = Counter(id(ix.unstable_manifold) for ix in crossings)
        indexed = {id(ix.manifold_a_key[0]) for ix in crossings if ix.is_anchor}
        held = {id(piece) for piece in self._partial_bridges}
        self._partial_bridges += [
            curve for curve in self._image_curves
            if id(curve.fixed_point) in indexed
            and (fixed_point is None or curve.fixed_point is fixed_point)
            and crossings_on[id(curve)] < 2
            and id(curve) not in held
        ]
        return cut

    def clear_bridges(self, fixed_point: Optional[FixedPoint] = None) -> None:
        """
        Discard the registered bridges (e.g. before recutting after a retrim).

        A bridge owns the segments it was cut over (it shares them with the
        unstable manifold it came from), so discarding it must also release those
        claims -- otherwise the dead object stays reachable from the Tangle's
        index and its segments are never collected. The parent's ownership is
        untouched: only a segment whose LAST owner was the bridge is dropped.

        Note:
            Every production call reaches here just after a recompute, which has
            already run ``Tangle.clear_all()`` -- so the release is usually a no-op
            on an empty index. It matters on the direct path (``clear_bridges`` or
            ``rebuild_bridges`` with no recompute in between), where the index
            would otherwise leak.

        Args:
            fixed_point: If given, discard only that fixed point's bridges (and
                partial pieces); the others keep their objects and ids.
        """
        def dropped(bridge: Bridge) -> bool:
            return fixed_point is None or bridge.fixed_point is fixed_point

        # A partial piece is an image curve: computed manifold, which stays indexed.
        for bridge in self._bridges.values():
            if dropped(bridge):
                self.Tangle.release_manifold(bridge)
        for bid in [bid for bid, bridge in self._bridges.items() if dropped(bridge)]:
            del self._bridges[bid]
        self._partial_bridges[:] = [b for b in self._partial_bridges if not dropped(b)]
        self._bridges_at.clear()
        for bid in self._bridges:
            for endpoint in bid:
                self._bridges_at.setdefault(endpoint, []).append(bid)
        self._bump_generation()

    def rebuild_bridges(
        self, fixed_point: Optional[FixedPoint] = None
    ) -> list[Bridge]:
        """
        Clear the existing bridges and recut them against the current crossings.

        Bridges are segments of an unstable manifold between consecutive crossings, so
        trimming a stable manifold (e.g. to define a resonance zone) and recomputing
        intersections invalidates them: the shorter stable arc produces fewer crossings
        and therefore different bridges. Call this after such a recompute so the stale
        bridges are dropped and recut against the new (shorter) stable manifold.

        Per-bridge metadata is carried across BY ID: a recut bridge with the same
        :data:`~.Bridge.BridgeId` is the same bridge, so its ``iterated`` flag
        survives -- but ONLY while those ids still name the same crossings. A
        :meth:`compute_intersections` without ``preserve_ids`` rebuilds the registry
        from zero, and an old pair then names an unrelated pair of crossings; that is
        tracked by the id epoch, and on a mismatch the carry is dropped (logged at
        debug) rather than silently marking a bridge that was never iterated. With
        ``preserve_ids=True`` the recompute runs
        :meth:`IntersectionRegistry.reindex_from`, every re-detected crossing gets its
        old id back, and the carry stands.

        A blast's image curves are recut too: they are computed manifold, indexed by
        every recompute, so a rebuild keeps every image bridge, every partial
        piece (an image curve crossed fewer than twice, whose ``iterated`` flag lives
        on the curve itself) and every flag carried by id.

        Args:
            fixed_point: If given, only rebuild bridges for that fixed point's unstable
                manifolds (the other fixed points' bridges stay, same objects and
                ids); otherwise rebuild for every indexed unstable manifold. After a
                renumbering recompute every fixed point is rebuilt: the surviving
                bridges' ids no longer name their crossings.

        Returns:
            The freshly created bridges.
        """
        carried: dict[BridgeId, bool] = {}
        if self._bridge_id_epoch == self._registry_id_epoch:
            carried = {bid: bridge.iterated for bid, bridge in self._bridges.items()}
        elif self._bridges:
            logger.debug(
                "rebuild_bridges: the intersection registry was renumbered since the "
                "%d registered bridge id(s) were cut (id epoch %d -> %d), so they no "
                "longer name the same crossings; rebuilding every fixed point's "
                "bridges and dropping the per-bridge metadata",
                len(self._bridges),
                self._bridge_id_epoch,
                self._registry_id_epoch,
            )
            fixed_point = None

        self.clear_bridges(fixed_point)
        rebuilt = self.create_bridges(fixed_point)
        for bid, bridge in self._bridges.items():
            if bid in carried:
                bridge.iterated = carried[bid]
        return rebuilt

    @property
    def bridges(self) -> list[Bridge]:
        """
        All bridges registered so far (originals and iterated children).

        Identified bridges come first, in the order they were cut, followed by the
        partial pieces (an image bounded by fewer than two crossings), which have
        no :data:`~.Bridge.BridgeId` of their own.
        """
        return list(self._bridges.values()) + list(self._partial_bridges)

    @property
    def uniiterated_bridges(self) -> list[Bridge]:
        """All bridges that have not yet been iterated forward."""
        return [b for b in self.bridges if not b.iterated]

    def bridge(self, bridge_id: BridgeId) -> Bridge:
        """
        The registered bridge with this :data:`~.Bridge.BridgeId`.

        Args:
            bridge_id: ``(first, second)`` crossing ids in unstable cdist order.

        Returns:
            Bridge: The single registered copy.

        Raises:
            KeyError: If no bridge with that id is registered.
        """
        return self._bridges[bridge_id]

    def bridges_at(self, intersection_id: int) -> list[BridgeId]:
        """
        The ids of the bridges that end at one crossing.

        A crossing interior to an unstable branch bounds two bridges (the one
        before and the one after it); an outermost crossing bounds one; a crossing
        on a stretch that has not been cut bounds none.

        Args:
            intersection_id: Registry id of the crossing.

        Returns:
            list[BridgeId]: The ids, in the order the bridges were cut.
        """
        return list(self._bridges_at.get(intersection_id, ()))

    # ── bridge iteration and genealogy (see BridgeIterator) ─────────────────

    def image_bridges(
        self, bridge_id: BridgeId, n: int = 1
    ) -> Optional[list[BridgeId]]:
        """
        The bridges tiling the ``n``-th forward image of one bridge.

        Delegates to :meth:`~.BridgeIterator.BridgeIterator.image_bridges`, which
        documents the derivation and the ``None`` case.
        """
        return self._bridge_iterator.image_bridges(bridge_id, n)

    def preimage_bridges(
        self, bridge_id: BridgeId, n: int = 1
    ) -> Optional[list[BridgeId]]:
        """
        The bridges tiling the ``n``-th BACKWARD image of one bridge.

        Delegates to :meth:`~.BridgeIterator.BridgeIterator.preimage_bridges`.
        """
        return self._bridge_iterator.preimage_bridges(bridge_id, n)

    @property
    def intersection_registry(self) -> IntersectionRegistry:
        """The single source of truth for this workbench's crossings.

        Note:
            :meth:`compute_intersections` replaces this object (unless
            ``reset=False``), so re-fetch it after every call rather than
            holding a reference across one.
        """
        return self._intersection_registry

    def iterate_bridge(self, bridge: Bridge) -> list[Bridge]:
        """
        Map a bridge forward one iterate and cut the image into bridges.

        Delegates to :meth:`~.BridgeIterator.BridgeIterator.iterate_bridge`, which
        documents the reuse path, the partial-bridge case and what is raised.
        """
        return self._bridge_iterator.iterate_bridge(bridge)

    def iterate_all_bridges(self) -> list[Bridge]:
        """
        Iterate every bridge that has not yet been mapped forward.

        Delegates to :meth:`~.BridgeIterator.BridgeIterator.iterate_all_bridges`.
        """
        return self._bridge_iterator.iterate_all_bridges()

    # ── iterate table inference (see IterateInference) ───────────────────────

    def infer_iterates(self, sources: Optional[Iterable[int]] = None) -> int:
        """
        Record the ``+1`` image of every source crossing, by the exact bracket rule.

        Delegates to :func:`~.IterateInference.infer_iterates`;
        :meth:`compute_intersections`, every :meth:`iterate_bridge` and the end
        of every blast call it, so the table is complete without a manual call.

        Args:
            sources: Registry ids to find images for; None means all of them.

        Returns:
            Number of new iterate relationships recorded.
        """
        return infer_iterates(self, sources)

    def plot_all_bridges(self, bridges: Optional[list[Bridge]] = None) -> None:
        """
        Plot a list of bridges. If no list is supplied, plots all registered bridges.

        Args:
            bridges: List of bridges to plot. Defaults to every registered bridge.
        """
        if bridges is None:
            bridges = self.bridges
        n = len(bridges)
        if n == 0:
            return
        # Cycle a fixed cool-color palette by index (mod len) rather than
        # resampling a colormap across all n bridges. Resampling makes adjacent
        # bridges nearly identical in hue, so a run of consecutive bridges reads as
        # a single colour spanning several intersections; cycling guarantees
        # neighbouring bridges are always visually distinct. The palette is
        # restricted to cool hues (blues, cyans, teals, greens, purples) so that
        # no bridge is ever confused with the stable manifold, which is always red.
        palette = self._COOL_BRIDGE_PALETTE
        for i, bridge in enumerate(bridges):
            bridge.plot(color=palette[i % len(palette)])

    def trim_stable_at(self, key: ManifoldKey, cdist: float) -> BaseManifold:
        """
        Cut a stable branch at the first node past its root at or beyond ``cdist``.

        The one stable trim: :meth:`trim_stable_manifolds` cuts at the outermost
        crossing, a resonance zone at its pip. The node is the one the crossing's
        own stable segment ends at, so the crossing stays on the branch. Only the
        tail pointer moves, so the next growth of that stability undoes the trim
        (:meth:`grow_n_times`); recompute the intersections afterwards.

        Args:
            key: The stable branch.
            cdist: The stable canonical distance to keep.

        Returns:
            The (trimmed) stable manifold. A branch with no node at or past
            ``cdist`` is left untrimmed, with a WARNING.
        """
        manifold = self.manifolds[key]
        new_tail = next(
            (
                node
                for node in manifold._walk_nodes(stop_at_final=False)
                if node is not manifold.root and node.get_cdist("stable") >= cdist
            ),
            None,
        )
        if new_tail is None:
            logger.warning(
                "trim_stable_at: %s has no node at or past stable cdist %r; leaving it "
                "untrimmed",
                key,
                cdist,
            )
        else:
            manifold.tail = new_tail
        return manifold

    def trim_stable_manifolds(self, fixed_point: FixedPoint) -> None:
        """
        Trim each stable branch of a fixed point just past its outermost crossing.

        Args:
            fixed_point (FixedPoint): fixed point manifolds will be
                trimmed from

        Note:
            A stable branch that nothing crosses yet -- a branch of a period > 1
            orbit whose crossings have not developed, or a second fixed point whose
            tangle has not reached this one -- is left untouched rather than
            trimmed. There is no outermost crossing to trim to.

        Note:
            The outermost crossing of N is read off the registry (grouped by the
            stable branch key each crossing records), the single source of truth
            for crossings; an unresolved detection never sets the trim. The cut
            is :meth:`trim_stable_at`.
        """
        registry = self._intersection_registry
        outermost: dict[ManifoldKey, float] = {}
        for iid in registry.resolved_ids(self.resolution_ratio_threshold):
            intersection = registry[iid]
            key = intersection.manifold_b_key
            outermost[key] = max(outermost.get(key, -np.inf), float(intersection.stable_cdist))

        for manifold in self._iter_manifolds(fixed_point, "stable"):
            target = outermost.get(manifold.manifold_key)
            # A branch whose only crossing is the anchor (stable cdist 0) has no
            # outermost crossing to trim to: trimming there would collapse the
            # manifold onto its own root and destroy the branch.
            if target is None or target <= 0.0:
                logger.debug(
                    "trim_stable_manifolds: no crossing past the anchor on %s, "
                    "leaving it untrimmed",
                    manifold.manifold_key,
                )
                continue
            self.trim_stable_at(manifold.manifold_key, target)

    @staticmethod
    def _stretch_per_round(fixed_point: FixedPoint, stability: Stability) -> float:
        """The canonical-distance factor of one growth round on a branch."""
        beta = fixed_point.per_step_beta(stability)
        # stable growth runs the inverse map, which stretches by 1 / beta_s
        return beta if stability == "unstable" else 1.0 / beta

    def _default_cap(self, fixed_point: FixedPoint, stability: Stability, log_stretch: float) -> int:
        """A driver's default round cap: the rounds the needed stretch takes, never below 10."""
        rounds = math.ceil(log_stretch / math.log(self._stretch_per_round(fixed_point, stability)))
        return max(10, rounds + 1)

    def _cap_reached(
        self,
        fixed_point: FixedPoint,
        stability: Stability,
        max_iterations: int,
        log_stretch: float,
        what: str,
    ) -> MaxIterationsReached:
        """The cap error of a ``grow_until_*`` driver, with its suggested cap."""
        return MaxIterationsReached(
            max_iterations, self._stretch_per_round(fixed_point, stability), log_stretch, what
        )

    def _require_manifold(
        self,
        fixed_point: FixedPoint,
        stability: Stability,
        branch_index: Optional[int] = 0,
    ) -> BaseManifold:
        """
        The manifold a growth driver is about to grow, or a ValueError.

        Args:
            fixed_point (FixedPoint): The manifold's fixed point.
            stability (Stability): Which manifold.
            branch_index (Optional[int]): Which eigenvector branch; ``None`` is
                read as branch 0.

        Returns:
            BaseManifold: The manifold registered under
            ``(fixed_point, stability, 0, branch_index)``.

        Raises:
            ValueError: No such manifold has been initialized.
        """
        key = (fixed_point, stability, 0, branch_index if branch_index is not None else 0)
        manifold = self.manifolds.get(key)
        if manifold is None:
            raise ValueError(
                f"Manifold for fixed point {fixed_point} with stability {stability} "
                f"and branch_index {branch_index} has not been initialized."
            )
        return manifold

    def _iter_manifolds(
        self, fp: FixedPoint, stability: Optional[Stability] = None
    ) -> Iterable[BaseManifold]:
        """Yield all manifolds for a fixed point (optionally filter by stability)."""
        for (kfp, kstab, _oi, _bi), M in self.manifolds.items():
            if kfp is fp and (stability is None or kstab == stability):
                yield M

