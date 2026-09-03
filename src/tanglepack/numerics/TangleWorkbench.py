"""The orchestrator of the numerical layer.

:class:`TangleWorkbench` owns the objects a computed tangle is made of -- the
dynamical system, the fixed points, the manifolds, the :class:`~.Tangle.Tangle`
index, the :class:`~.IntersectionRegistry.IntersectionRegistry` and the bridge
set -- and drives construction, growth, intersection detection, bridge cutting,
trimming and manifold plotting.

Three collaborators own the work that reads and extends that state without
being part of the orchestration:

* :mod:`~.BridgeIterator` -- mapping a bridge forward, and deriving the bridge
  genealogy from the crossing iterate table;
* :mod:`~.IterateInference` -- filling that iterate table;
* :mod:`~.graphviz` -- the networkx/matplotlib view of the intersection graph.

Each is reached through a same-named method on the workbench, so no caller has
to know where the body lives. The Dev Notes for each of those areas live in the
module that owns it.

Dev Notes:

``grow_until_arclength`` stops on the tail's CANONICAL distance, not on the true
arc length along the curve. The two agree up to the per-branch scaling only while
the curve is close to the linear regime, so the stop length is approximate; using
the accumulated ``edist`` would make it exact.
"""

from __future__ import annotations

import logging
from typing import Callable, Literal, Iterable, Optional, Sequence, Union

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

from . import graphviz
from .BridgeIterator import BridgeIterator
from .DynamicalSystem import DynamicalSystem, MapFunc, JacFunc
from .FixedPointSolver import FixedPointSolver
from .ManifoldInitializer import ManifoldInitializer
from .ManifoldMachine import ManifoldMachine
from .Tangle import Tangle
from .FixedPoint import FixedPoint
from .BaseManifold import BaseManifold
from .Bridge import Bridge, BridgeId
from .Intersection import Intersection, ManifoldKey, Stability
from .IntersectionRegistry import IntersectionRegistry
from .IterateInference import IterateInference


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


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

    Note:
        Bridge iteration, iterate-table inference and graph plotting are
        delegated to :class:`~.BridgeIterator.BridgeIterator`,
        :class:`~.IterateInference.IterateInference` and :mod:`~.graphviz`. The
        methods here forward to them and their signatures are the contract; the
        bodies (and their Dev Notes) live in those modules.
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
    ) -> None:
        """
        Build an empty workbench around one dynamical system.

        Args:
            dynamical_map (MapFunc): The area-preserving map.
            dynamical_map_inverse (MapFunc): Its inverse.
            jacobian_function (Optional[JacFunc]): Analytic Jacobian; a
                finite-difference one is used when it is omitted.
        """

        self.dynamical_system = DynamicalSystem(
            dynamical_map, dynamical_map_inverse, jacobian_function
        )

        self._fp_solver = FixedPointSolver(self.dynamical_system)
        self._man_maker = ManifoldInitializer(self.dynamical_system)
        self._man_machine = ManifoldMachine(self.dynamical_system)

        self.Tangle = Tangle()

        self.fixed_points = []
        # manifolds are keyed like (fixed_point, stability, orbit_index, branch_index)
        self._manifolds: dict[tuple[FixedPoint, Stability, int, int], BaseManifold] = {}

        self._intersection_registry = IntersectionRegistry()
        # Bridges are keyed by their topological identity (the ordered endpoint
        # crossing pair), so there is exactly one object per computed bridge and
        # genealogy can be derived rather than stored. Partial pieces have no such
        # identity and are held separately; ``_bridges_at`` is the reverse index
        # from a crossing id to the bridges that end there.
        self._bridges: dict[BridgeId, Bridge] = {}
        self._partial_bridges: list[Bridge] = []
        self._bridges_at: dict[int, list[BridgeId]] = {}
        # A BridgeId is a pair of REGISTRY ids, so it only means anything while the
        # registry still numbers its crossings the same way. This counter is bumped
        # every time a recompute renumbers them (see compute_intersections);
        # _bridge_id_epoch is its value when the current bridges were cut, and
        # rebuild_bridges refuses to carry metadata across a mismatch.
        self._registry_id_epoch: int = 0
        self._bridge_id_epoch: int = 0

        # Generation bookkeeping (see the `generation` property). `_mutations`
        # counts the changes this object makes that nothing else can observe
        # (the bridge set, the registry swap); `_generation_state` is the last
        # composite state seen and `_generation` the monotone token handed out.
        self._mutations: int = 0
        self._generation_state: Optional[tuple] = None
        self._generation: int = 0

        # The two collaborators that read and extend this state (see the module
        # docstring). They are constructed here, hold no state of their own, and
        # the workbench forwards its same-named methods to them.
        self._bridge_iterator = BridgeIterator(self)
        self._iterate_inference = IterateInference(self)

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
        including the lazy rebuilds of derived caches (the registry graph's
        adjacency edges, the rank maps).

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
            FixedPoint: The constructed fixed point, already registered here. Its
                branch count is derived (``FixedPoint.num_branches``).
        """

        fixed_point = self._fp_solver.construct_fixed_point(initial_guess)

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
                two unstable directions and the caller may want both. It is
                ignored on an inversion point, where the single chain of
                ``k_value`` pieces already covers both branches
                (``FixedPoint.num_branches`` is 2 there).

        Returns:
            The initial segments, keyed by ``(orbit_index, branch_index)``.

        Raises:
            ValueError: If ``num_branches`` is neither 1 nor 2.
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

        Note:
            Growth refines the existing curve as well as extending it, so it
            bumps the generation and drops every memoised point array.
        """

        self._require_manifold(fixed_point, stability, branch_index)

        self._man_machine.grow_x_times(
            fixed_point, stability, num_iterations, branch_index
        )

        for (fp, stab, _orbit_index, _branch_index), manifold in self.manifolds.items():

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
        max_iterations: int = 10,
        branch_index: int = 0,
    ) -> None:
        """
        Grows the manifold until a turnaround is detected or max_iterations is reached.

        Args:
            fixed_point (FixedPoint): The fixed point whose manifold is to be grown.
            stability (Stability): The stability type of the manifold ('stable' or 'unstable').
            max_iterations (int, optional): Maximum number of growth iterations
                to attempt before giving up. Defaults to 10.
            branch_index (int, optional): Which eigenvector branch to grow.
                Defaults to 0.

        Raises:
            ValueError: If the manifold has not been initialized, or if the cap
                is reached before a turnaround is detected.
        """

        manifold = self._require_manifold(fixed_point, stability, branch_index)

        root = fixed_point.branch_points[0]
        first_point = manifold.first_node(branch_index)

        root_coord = root._coords
        first_point_coords = first_point._coords

        initial_direction = np.asarray(first_point_coords) - np.asarray(root_coord)

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

        else:
            raise ValueError(
                "Max iterations reached, choose a higher cap or a different method."
            )

    def grow_until_arclength(
        self,
        fixed_point: FixedPoint,
        stability: Stability,
        length: float,
        max_iterations: int = 10,
        branch_index: int = 0,
    ) -> None:
        """
        Grow the manifold until its tail passes the requested canonical distance.

        Args:
            fixed_point (FixedPoint): The fixed point whose manifold is grown.
            stability (Stability): Which manifold to grow.
            length (float): Target canonical distance of the tail.
            max_iterations (int): Growth iterations to attempt before giving up.
                Defaults to 10, matching :meth:`grow_until_turnaround`.
            branch_index (int): Which eigenvector branch to grow. Defaults to 0.

        Raises:
            ValueError: If the manifold has not been initialized, or if the cap is
                reached before the tail passes ``length``.
        """

        self._require_manifold(fixed_point, stability, branch_index)

        def tail_distance() -> float:
            return self.manifolds[(fixed_point, stability, 0, branch_index)].tail.cdist

        for _ in range(max_iterations):

            if tail_distance() >= length:
                return None

            self.grow_n_times(
                fixed_point, stability, num_iterations=1, branch_index=branch_index
            )

        if tail_distance() >= length:
            return None

        raise ValueError(
            "Max iterations reached, choose a higher cap or a shorter length."
        )

    def grow_until_intersection(
        self,
        fixed_point: FixedPoint,
        stability: Stability,
        max_iterations: int = 10,
        branch_index: int = 0,
    ) -> None:
        """
        Grow the manifold until a new crossing appears.

        Args:
            fixed_point (FixedPoint): The fixed point whose manifold is grown.
            stability (Stability): Which manifold to grow.
            max_iterations (int): Growth iterations to attempt before giving up.
                Defaults to 10.
            branch_index (int): Which eigenvector branch to grow. Defaults to 0.

        Raises:
            ValueError: If the manifold has not been initialized, or if the cap is
                reached without a new crossing.
        """

        self._require_manifold(fixed_point, stability, branch_index)

        self.compute_intersections(fixed_point, infer_iterates=False)

        num_initial_intersections = len(self._intersection_registry)
        logger.info("current intersections %d", num_initial_intersections)

        for _ in range(max_iterations):

            self.grow_n_times(
                fixed_point, stability, num_iterations=1, branch_index=branch_index
            )

            self.compute_intersections(fixed_point, infer_iterates=False)
            num_current_intersections = len(self._intersection_registry)

            if num_current_intersections > num_initial_intersections:
                return None

        raise ValueError("Max iterations reached, no intersection found")

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
                INVERSION point it is ignored (one chain covers both branches) and
                ``None`` is forwarded to :meth:`grow_n_times` so the per-round
                "branch_index is ignored" warning is not emitted ten times.

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
            curve. That recut DROPS the image bridges and partial pieces produced
            by :meth:`iterate_bridge`, so a blast frontier does not survive a
            growth round -- see :meth:`_recut_bridges_after_growth`, which warns
            about what it cost and clears the ``iterated`` flag of every bridge
            whose image can no longer be derived.
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

        # An inversion point's two branches are one chain; passing a branch index
        # would only earn a warning from ManifoldMachine on every single round.
        growth_branch = None if fixed_point.check_inversion() else branch_index

        if predicate(self):
            return 0

        for round_index in range(1, max_iterations + 1):

            for stability in stabilities:
                self.grow_n_times(
                    fixed_point, stability, num_iterations=1, branch_index=growth_branch
                )

            self.compute_intersections(
                self._fixed_points_with_manifolds(),
                infer_iterates=True,
                preserve_ids=True,
            )
            if self._bridges or self._partial_bridges:
                self._recut_bridges_after_growth()

            if predicate(self):
                return round_index

        raise ValueError(
            f"Max iterations ({max_iterations}) reached without the growth "
            "predicate becoming true; choose a higher cap or a different "
            "stopping condition."
        )

    def _recut_bridges_after_growth(self) -> None:
        """
        Recut the bridge set against the grown curves, and repair the blast state.

        Two things go wrong if a growth round leaves the bridges alone or recuts
        them naively:

        * an uncut bridge is no longer a pair of CONSECUTIVE crossings on its
          branch (growing the stable manifold drops new crossings inside it), and
          the arrangement asserts exactly that;
        * :meth:`rebuild_bridges` carries the ``iterated`` flag across by id, but
          the recut DROPS the image bridges and partial pieces that flag was
          recording. A bridge left marked ``iterated`` whose
          :meth:`image_bridges` is ``None`` is invisible to
          :meth:`uniiterated_bridges` and can never be re-derived, so a later
          blast silently skips it.

        The flag is therefore cleared wherever the image is gone, and the whole
        repair is logged at WARNING: losing a frontier is a real cost to whoever
        was mid-blast, not an implementation detail.

        Note:
            Only a ``None`` image clears the flag. An EMPTY image list is a
            different statement -- the image span is known and simply carries no
            cut bridges yet -- and re-iterating that bridge would duplicate work
            rather than recover anything.
        """
        before_ids = set(self._bridges)
        partials_before = len(self._partial_bridges)

        self.rebuild_bridges()

        dropped = sorted(before_ids - set(self._bridges))
        cleared: list[BridgeId] = []
        for bid, bridge in self._bridges.items():
            if bridge.iterated and self.image_bridges(bid) is None:
                bridge.iterated = False
                cleared.append(bid)

        if dropped or partials_before or cleared:
            logger.warning(
                "growth round recut the bridge set: dropped %d identified bridge(s) "
                "%s and %d partial piece(s); cleared the 'iterated' flag on %d "
                "bridge(s) %s whose forward image can no longer be derived. A blast "
                "frontier does not survive a growth round.",
                len(dropped),
                dropped,
                partials_before,
                len(cleared),
                sorted(cleared),
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
        goalposts. ``"all"`` means every non-anchor crossing OF ``fixed_point`` --
        scoping it to the whole registry would put another tangle's crossings in a
        set that only this fixed point's growth can ever close.
        """
        registry = self._intersection_registry
        if isinstance(ids, str):
            if not allow_all or ids != "all":
                raise ValueError(
                    f"Unknown id selector {ids!r}; pass an iterable of registry ids"
                    + (' or "all".' if allow_all else ".")
                )
            assert fixed_point is not None, '"all" needs a fixed point to scope to'
            id_of = {id(ix): iid for iid, ix in registry}
            selected = sorted(
                id_of[id(ix)]
                for ix in registry.from_fixed_point(fixed_point)
                if not ix.is_synthetic
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
        deliberately not a single-fixed-point snapshot: that would cut a
        heteroclinic partner's arcs off and turn real faces into open ones.

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
            An ANCHOR never closes, but that is a limitation of the arrangement's
            anchor model rather than a fact about periodic points. A periodic point
            is one degree-four node whose ``u-`` and ``s-`` rays are the OTHER
            branch's ``u+`` and ``s+``;
            :meth:`_register_anchors` instead registers one anchor per (unstable
            branch, stable branch) pair, and the arrangement has no notion of a
            branch continuing THROUGH a node, so those two slots are filled with
            virtual stubs and every face through one of them is open. That is the
            deferred item recorded in the "Inversion caveat" Dev Note of
            :mod:`tanglepack.topology.Arrangement`; until it is closed, asking for
            an anchor's faces runs to the cap.
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
    ) -> list[tuple[float, float]]:
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
                of every intersection, by canonical-distance mapping). Pass False to
                skip it — used by the growth loops that call this many times and do not
                need the table.
            preserve_ids: If True (and ``reset`` is True), re-align the rebuilt
                registry against the one being replaced so any crossing that reappears
                keeps its previous id (see :meth:`IntersectionRegistry.reindex_from`).
                Used by the resonance-zone recompute so a strong pip chosen as id N is
                still id N after the stable manifolds are trimmed. Defaults to False —
                the growth loops renumber freely.

        Returns:
            List of (x, y) coordinates, one per detected crossing.
        """
        if isinstance(fixed_points, FixedPoint):
            fixed_points = [fixed_points]

        old_registry = self._intersection_registry
        if reset:
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

        # Anchors FIRST: every periodic point is a crossing of its own manifolds,
        # and declaring it beats detecting it (see _register_anchors).
        self._register_anchors(fixed_points)

        for intersection in self.Tangle.resolve_crossings():
            intersection.id = self._intersection_registry.add(intersection)

        if reset and preserve_ids and len(old_registry) > 0:
            self._intersection_registry.reindex_from(old_registry)
        elif reset:
            # The registry was rebuilt from zero: every id held by a caller (a strong
            # pip, a BridgeId) now names a different crossing, if anything at all.
            self._registry_id_epoch += 1

        if infer_iterates:
            self.infer_iterates()

        return [ix.coords for _iid, ix in self._intersection_registry]

    def _register_anchors(self, fixed_points: Iterable[FixedPoint]) -> list[int]:
        """
        Declare every periodic point as a crossing of its own two manifolds.

        A periodic point sits at canonical distance 0 on every branch attached to
        it, so it IS a crossing of each of its unstable branches with each of its
        stable branches -- the ``(0, 0)`` anchor the topological layer treats as
        the innermost node of a branch. Until now that crossing appeared only by
        accident, when the two fundamental segments happened to be detected as an
        rtree pair; declaring it makes it a fact of the model rather than a
        by-product of the index.

        Anchors are registered BEFORE :meth:`Tangle.resolve_crossings`, so:

        * they take the lowest registry ids, deterministically ordered by the
          manifold insertion order rather than by whatever order the rtree
          returned its candidate pairs;
        * the accidentally detected crossing at the same place collides with the
          anchor in :meth:`IntersectionRegistry.add` (same ``(0, 0)`` cdists,
          same two branch keys) and is dropped in its favour -- one crossing, not
          two, and the id belongs to the declared one.

        Purging the detected copy afterwards instead would work equally well
        geometrically but would leave the anchor's id dependent on detection
        order, which is exactly what ``preserve_ids`` has to be able to rely on.

        The anchor carries the same cut geometry a detected crossing carries --
        its unstable branch's root and the node after it -- so
        :meth:`Tangle.create_bridges` cuts at it exactly as before and the bridge
        set is unchanged.

        Args:
            fixed_points: The fixed points whose periodic points to anchor.

        Returns:
            The registry ids of the anchors, in registration order.

        Note:
            One anchor is registered per (unstable branch, stable branch) pair
            meeting at each periodic point. For a point without inversion that is
            the single pair (branch 0, branch 0). An inversion point has two
            branches of each manifold and so gets FOUR anchors at one coordinate,
            which keeps every branch anchored but is not the topologically right
            model -- the four rays leaving an inversion point form ONE degree-4
            node. See the Dev Notes in
            :mod:`tanglepack.topology.Arrangement`.
        """
        registry = self._intersection_registry
        anchor_ids: list[int] = []
        for fp in fixed_points:
            unstable = self._branch_items(fp, "unstable")
            stable = self._branch_items(fp, "stable")
            for u_key, u_manifold in unstable:
                u_nodes = u_manifold.get_point_array(return_nodes=True)
                if len(u_nodes) < 2:
                    logger.debug(
                        "No anchor on %s: the branch has fewer than two nodes",
                        u_key,
                    )
                    continue
                for s_key, s_manifold in stable:
                    if s_key[2] != u_key[2]:
                        continue  # a different point of the periodic orbit
                    s_nodes = s_manifold.get_point_array(return_nodes=True)
                    if len(s_nodes) < 2:
                        continue
                    coords = u_nodes[0].get_point()
                    anchor = Intersection.synthetic(
                        coords=(float(coords[0]), float(coords[1])),
                        unstable_cdist=0.0,
                        stable_cdist=0.0,
                        label="anchor",
                        manifold_a_key=u_key,
                        manifold_b_key=s_key,
                        unstable_manifold=u_manifold,
                        unstable_segment=(u_nodes[0], u_nodes[1]),
                        crossing_sign=self._anchor_crossing_sign(u_nodes, s_nodes),
                    )
                    anchor_ids.append(registry.add(anchor))
        return anchor_ids

    def _branch_items(
        self, fp: FixedPoint, stability: Stability
    ) -> list[tuple[ManifoldKey, BaseManifold]]:
        """The ``(key, manifold)`` pairs of one fixed point's branches."""
        return [
            (key, manifold)
            for key, manifold in self.manifolds.items()
            if key[0] is fp and key[1] == stability
        ]

    @staticmethod
    def _anchor_crossing_sign(u_nodes: list, s_nodes: list) -> int:
        """
        Handedness of an anchor, from the two oriented eigendirections.

        The eigendirections are read off the fundamental segments rather than off
        ``FixedPoint``: the initializer lays each fundamental segment along the
        oriented eigenvector, so the first node-to-node step IS that eigenvector
        -- and, unlike the stored eigenvector, it already carries the sign of the
        branch it belongs to, which is what an inversion point's second branch
        needs (there the two branches run along opposite eigendirections and one
        stored vector cannot describe both).

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
            self.compute_intersections(fp, infer_iterates=False)

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

        self._bridges[bid] = bridge
        for endpoint in bid:
            self._bridges_at.setdefault(endpoint, []).append(bid)
        return bridge

    def create_bridges(self, fixed_point: Optional[FixedPoint] = None) -> list[Bridge]:
        """
        Cut indexed unstable manifolds into bridges.

        Args:
            fixed_point: If given, only build bridges for that fixed point's
                unstable manifolds. If None, build bridges for every indexed
                unstable manifold at once. Each bridge keeps its own fixed_point
                linkage either way.

        Returns:
            The bridges cut, each the single registered copy for its
            :data:`~.Bridge.BridgeId`.
        """
        crossings = [ix for _iid, ix in self._intersection_registry]
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
        return [self._register_bridge(bridge) for bridge in bridges]

    def clear_bridges(self) -> None:
        """
        Discard every registered bridge (e.g. before recutting after a retrim).

        A bridge owns the segments it was cut over (it shares them with the
        unstable manifold it came from), so discarding it must also release those
        claims -- otherwise the dead object stays reachable from the Tangle's
        index and its segments are never collected. The parent's ownership is
        untouched: only a segment whose LAST owner was the bridge is dropped.

        Note:
            Every production call reaches here just after a recompute, which has
            already run ``Tangle.clear_all()`` -- so the release is usually a no-op
            on an empty index. It matters on the direct path (``clear_bridges`` or
            ``rebuild_bridges`` with no recompute in between), which is exactly the
            path that used to leak.
        """
        for bridge in self.bridges:
            self.Tangle.release_manifold(bridge)
        self._bridges.clear()
        self._partial_bridges.clear()
        self._bridges_at.clear()
        self._bump_generation()

    def rebuild_bridges(
        self, fixed_point: Optional[FixedPoint] = None
    ) -> list[Bridge]:
        """
        Clear all existing bridges and recut them against the current crossings.

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

        Only bridges cut out of an INDEXED unstable manifold are recut; the image
        bridges produced by :meth:`iterate_bridge` are not, so they are dropped by the
        rebuild as they always have been.

        Args:
            fixed_point: If given, only rebuild bridges for that fixed point's unstable
                manifolds; otherwise rebuild for every indexed unstable manifold.

        Returns:
            The freshly created bridges.

        Note:
            The CLEAR is global even when ``fixed_point`` is given: every registered
            bridge is discarded and only the named fixed point's are recut, so a
            per-fixed-point rebuild drops the other fixed points' bridges. That is
            pre-existing behaviour and every production caller rebuilds everything.
        """
        carried: dict[BridgeId, bool] = {}
        if self._bridge_id_epoch == self._registry_id_epoch:
            carried = {bid: bridge.iterated for bid, bridge in self._bridges.items()}
        elif self._bridges:
            logger.debug(
                "rebuild_bridges: the intersection registry was renumbered since the "
                "%d registered bridge id(s) were cut (id epoch %d -> %d), so they no "
                "longer name the same crossings; dropping the per-bridge metadata",
                len(self._bridges),
                self._bridge_id_epoch,
                self._registry_id_epoch,
            )

        self.clear_bridges()
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
        """The single source of truth for this workbench's crossings."""
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

    def infer_iterate_table(self, cdist_rtol: float = 0.05) -> int:
        """
        Record the n=1 forward iterate of every iterated bridge's endpoints.

        Delegates to
        :meth:`~.IterateInference.IterateInference.infer_iterate_table`.
        """
        return self._iterate_inference.infer_iterate_table(cdist_rtol)

    def infer_iterates(self, cdist_rtol: float = 0.05) -> int:
        """
        Fill the iterate table for *every* registered crossing.

        Delegates to :meth:`~.IterateInference.IterateInference.infer_iterates`;
        :meth:`compute_intersections` calls it automatically.
        """
        return self._iterate_inference.infer_iterates(cdist_rtol)

    # ── graph view (see graphviz) ────────────────────────────────────────────

    def build_intersection_graph(self) -> nx.MultiDiGraph:
        """
        The registry's intersection graph, decorated with this workbench's bridges.

        Delegates to :func:`~.graphviz.build_intersection_graph`.

        Returns:
            nx.MultiDiGraph: A copy of the registry graph.
        """
        return graphviz.build_intersection_graph(self)

    def visualize_intersection_graph(
        self,
        G: nx.MultiDiGraph,
        layout: str = "auto",
        figsize: tuple[int, int] = (12, 8),
        display_mode: str = "auto",
        compact_threshold: int = 20,
        node_size: Optional[int] = None,
        label_mode: str = "id",
        node_color_by: str = "none",
        show_iterate_edges: bool = True,
        save_path: Optional[str] = None,
    ) -> tuple:
        """
        Draw an intersection graph.

        Delegates to :func:`~.graphviz.visualize_intersection_graph`, which
        documents every option.

        Returns:
            (fig, ax) matplotlib Figure and Axes.
        """
        return graphviz.visualize_intersection_graph(
            G,
            layout=layout,
            figsize=figsize,
            display_mode=display_mode,
            compact_threshold=compact_threshold,
            node_size=node_size,
            label_mode=label_mode,
            node_color_by=node_color_by,
            show_iterate_edges=show_iterate_edges,
            save_path=save_path,
        )

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

    def trim_stable_manifolds(self, fixed_point: FixedPoint) -> None:
        """
        Trims the stable manifolds attached to the fixed point
        to just after the last intersection point.

        Args:
            fixed_point (FixedPoint): fixed point manifolds will be
                trimmed from

        Note:
            A stable branch that nothing crosses yet -- a branch of a period > 1
            orbit whose crossings have not developed, or a second fixed point whose
            tangle has not reached this one -- is left untouched rather than
            trimmed. There is no outermost crossing to trim to.

        Note:
            The outermost crossing is read off the registry (grouped by the stable
            branch key each crossing records), not off the segment index: the
            registry is the single source of truth for resolved crossings. The new
            tail is then the first node of the branch at or past that crossing --
            the same node the old segment-based lookup returned.
        """

        outermost: dict[tuple, float] = {}
        for _iid, intersection in self._intersection_registry:
            key = intersection.manifold_b_key
            if key is None:
                continue
            cdist = float(intersection.stable_cdist)
            if cdist > outermost.get(key, -np.inf):
                outermost[key] = cdist

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

            # the first node at or past the outermost crossing -- the node the
            # crossing's own segment ends at. get_cdist takes the stability because
            # the branch point carries one cdist per stability.
            new_tail = None
            for node in manifold._walk_nodes(stop_at_final=False):
                if node is manifold.root:
                    continue
                if node.get_cdist("stable") >= target:
                    new_tail = node
                    break

            if new_tail is None:
                logger.warning(
                    "trim_stable_manifolds: %s has no node at or past its outermost "
                    "crossing (stable cdist %r); leaving it untrimmed",
                    manifold.manifold_key,
                    target,
                )
                continue

            manifold.tail = new_tail

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

