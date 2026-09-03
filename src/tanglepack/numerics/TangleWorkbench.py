import logging
from typing import Annotated, Callable, Literal, Iterable, Optional
import numpy.typing as npt
from numpy.typing import NDArray

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import networkx as nx

from .DynamicalSystem import DynamicalSystem, MapFunc, JacFunc
from .FixedPointSolver import FixedPointSolver
from .ManifoldInitializer import ManifoldInitializer
from .ManifoldMachine import ManifoldMachine
from .BranchPoint import BranchPoint
from .Tangle import Tangle
from .FixedPoint import FixedPoint
from .BaseManifold import BaseManifold
from .Bridge import Bridge, BridgeId
from .Intersection import Intersection, ManifoldKey
from .IntersectionRegistry import IntersectionRegistry

Stability = Literal["unstable", "stable"]

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

"""
Dev Notes:

Cutting an iterated bridge (plan row 2.2) only ever cuts at the crossings just
born on that image -- the ones whose ``unstable_manifold`` IS the image. A crossing
found on an edge the image SHARES with an already-indexed curve is registered but
not cut at here, because its bracketing points lie on the other curve; the cut that
uses it is the one made on that curve. Widening the cut set to every crossing whose
cdist falls in the image's span would reintroduce the cross-curve bridges plan row
1.1 removed, so it must not be done without a per-curve span test.

Partial bridges (an image piece bounded by fewer than two crossings) have no
endpoint pair, hence no ``BridgeId``, hence no identity to dedupe on. They are
kept out of ``_bridges`` in a separate ``_partial_bridges`` list: they are real
stretches of unstable manifold the blast frontier must carry forward, but they are
not bridges in the topological sense and never appear in a genealogy answer. Two
consequences to keep in mind:

* repeated blasting of the same region can accumulate several Bridge objects over
  the same partial arc (nothing identifies them as the same arc), which is why the
  nested period-3 blast script registers 141 bridges where the pre-2.2 code
  reported 140;
* ``image_bridges`` can therefore lose a mid-arc leading/trailing piece of an
  image, so the derived frontier is the bridges of the image, not the whole image.

Both are by definition, not a gap to be papered over with a cdist-only signature --
that guesswork is exactly what plan row 2.2 removed.
"""


class TangleWorkbench:

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
        jacobian_function: JacFunc | None = None,
    ):

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

    def construct_fixed_point(self, initial_guess) -> FixedPoint:
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
        self, fixed_point: FixedPoint, approx_dirs: dict[str, np.ndarray] | None = None
    ):
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
    ):
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

    def initialize_both_manifolds(self, fixed_point: FixedPoint, num_branches: int = 1):
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

        key = (
            fixed_point,
            stability,
            0,
            branch_index if branch_index is not None else 0,
        )
        if self.manifolds.get(key) is None:
            raise ValueError(
                f"Manifold for fixed point {fixed_point} with stability {stability} "
                f"and branch_index {branch_index} has not been initialized."
            )

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
            max_iterations (int, optional): Maximum number of iterations to grow. Defaults to 50.
        """

        if self.manifolds.get((fixed_point, stability, 0, branch_index)) is None:
            raise ValueError(f"""Manifold for fixed point {fixed_point} 
                    with stability {stability} has not been initialized.
                    Please run initialize_manifold first.""")

        root = fixed_point.branch_points[0]
        first_point = self.manifolds.get((fixed_point, stability, 0, branch_index))
        first_point = first_point.walk_fwd(None, root, branch_index)

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

        if self.manifolds.get((fixed_point, stability, 0, branch_index)) is None:
            raise ValueError(f"""Manifold for fixed point {fixed_point}
                    with stability {stability} has not been initialized.
                    Please run initialize_manifold first.""")

        # TODO change this so it uses the actual arclength
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

        if self.manifolds.get((fixed_point, stability, 0, branch_index)) is None:
            raise ValueError(f"""Manifold for fixed point {fixed_point}
                    with stability {stability} has not been initialized.
                    Please run initialize_manifold first.""")

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

    def compute_intersections(
        self,
        fixed_points,
        *,
        reset: bool = True,
        infer_iterates: bool = True,
        preserve_ids: bool = False,
    ):
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

    def plot_intersections(
        self, fp=None, ax=None, show_ids=False, id_fontsize=8, **scatter_kwargs
    ):
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
        bridges = self.Tangle.create_bridges(crossings, fixed_point=fixed_point)
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

    def image_bridges(
        self, bridge_id: BridgeId, n: int = 1
    ) -> Optional[list[BridgeId]]:
        """
        The bridges tiling the ``n``-th forward image of one bridge.

        Genealogy is DERIVED, not stored: a bridge spans the arc between crossings
        ``a`` and ``b``, so its image spans the arc between ``f^n(a)`` and
        ``f^n(b)`` on the branch ``n`` map steps forward
        (:meth:`FixedPoint.advance_key`). The bridges of that image are exactly the
        registered bridges of the image branch whose two endpoints both lie inside
        that span. Both iterates are read from the registry's iterate table -- no
        cdist guesswork -- so the answer is only as complete as the table.

        Args:
            bridge_id: The bridge whose image is wanted.
            n: Number of map steps; negative walks backward.

        Returns:
            The image bridges' ids in increasing unstable canonical distance, or
            ``None`` when either endpoint has no registered ``n``-iterate. The list
            can be empty when the image arc is grown but not yet cut.

        Raises:
            KeyError: If ``bridge_id`` is not registered.

        Note:
            Partial pieces have no id and are never returned, so an image whose
            leading (or trailing) stretch runs past its last crossing loses that
            stretch here. That is by definition: a partial arc is not a bridge.
        """
        return self._image_bridge_ids(self.bridge(bridge_id), n)

    def preimage_bridges(
        self, bridge_id: BridgeId, n: int = 1
    ) -> Optional[list[BridgeId]]:
        """
        The bridges tiling the ``n``-th BACKWARD image of one bridge.

        Exactly :meth:`image_bridges` with the sign of ``n`` flipped; see there for
        the semantics and the ``None`` case.

        Args:
            bridge_id: The bridge whose preimage is wanted.
            n: Number of backward map steps (positive means ``f^-n``).

        Returns:
            The preimage bridges' ids in increasing unstable canonical distance, or
            ``None`` when either endpoint has no registered ``-n``-iterate.

        Raises:
            KeyError: If ``bridge_id`` is not registered.
        """
        return self._image_bridge_ids(self.bridge(bridge_id), -n)

    def _image_span(
        self, bridge: Bridge, n: int
    ) -> Optional[tuple[tuple[int, int], tuple[float, float], ManifoldKey]]:
        """The ``n``-iterates of a bridge's endpoints, their cdists, and the image key.

        The two ids come back ORDERED BY UNSTABLE CDIST (low, high), matching the
        cdists returned alongside them. Which of ``f^n(a)``, ``f^n(b)`` is the lower
        is not assumed: an inversion branch or a bookkeeping slip could swap them,
        and a caller comparing an ordered chain against them must not silently invert.
        """
        if bridge.id is None:
            return None
        table = self._intersection_registry.iterate_table
        images = [table[endpoint, n] for endpoint in bridge.id]
        if any(image is None for image in images):
            return None

        registry = self._intersection_registry
        images.sort(key=lambda image: float(registry[image].unstable_cdist))
        lo, hi = (float(registry[image].unstable_cdist) for image in images)
        image_key = bridge.fixed_point.advance_key(bridge.manifold_key, n)
        return (images[0], images[1]), (lo, hi), image_key

    def _image_bridge_ids(self, bridge: Bridge, n: int) -> Optional[list[BridgeId]]:
        """The ids of the registered bridges inside a bridge's n-th image span."""
        if n == 0:
            return [bridge.id] if bridge.id is not None else None

        span = self._image_span(bridge, n)
        if span is None:
            return None
        _images, (lo, hi), image_key = span

        registry = self._intersection_registry
        tol = registry.cdist_tol
        inside: list[tuple[float, BridgeId]] = []
        for bid, other in self._bridges.items():
            if other.manifold_key != image_key:
                continue
            first = float(registry[bid[0]].unstable_cdist)
            second = float(registry[bid[1]].unstable_cdist)
            if first >= lo - tol and second <= hi + tol:
                inside.append((first, bid))

        return [bid for _cdist, bid in sorted(inside)]

    @property
    def intersection_registry(self) -> IntersectionRegistry:
        return self._intersection_registry

    def iterate_bridge(self, bridge: Bridge) -> list[Bridge]:
        """
        Map a bridge forward one iterate, add the result to the tangle, detect new
        intersections with the stable manifold, cut the result into new bridges,
        and return those bridges.

        Marks the original bridge as iterated. The genealogy is not stored on the
        bridge: ``workbench.image_bridges(bridge.id)`` derives it from the crossing
        iterate table afterwards.

        Args:
            bridge: A bridge created by create_bridges() or a previous iterate_bridge().

        Returns:
            List of new Bridge objects from cutting the iterated result. If the
            image crosses the stable manifold fewer than twice, returns a
            single-element list holding the unsplit image, marked
            :attr:`Bridge.partial` (it is bounded by fewer than two crossings, so
            it is not a bridge in the topological sense).

        Raises:
            ValueError: If bridge has already been iterated.
            ValueError: If create_bridges() has not been called yet.
        """
        if bridge.iterated:
            raise ValueError(
                "This bridge has already been iterated. Its image is "
                "workbench.image_bridges(bridge.id)."
            )
        if not self._bridges and not self._partial_bridges:
            raise ValueError(
                "No bridges registered. Call create_bridges() before iterate_bridge()."
            )

        # 0. If this bridge's forward image is a stretch of manifold that has already
        #    been grown, that section is already cut into bridges -- reuse those
        #    existing objects instead of mapping the points forward again. Re-mapping
        #    would lay a second, slightly-different polyline over curve that already
        #    exists (the "zig-zag"), spawn overlapping duplicate bridges, and trip the
        #    same-stability (unstable x unstable) detector along the near-coincident
        #    pair. Only when the image runs past the grown extent is there genuinely
        #    new curve to compute (handled below).
        existing_image = self._existing_image_bridges(bridge)
        if existing_image is not None:
            bridge.iterated = True
            self._bump_generation()
            return existing_image

        # 1. map forward
        iterated = self._man_machine.iterate_bridge(bridge)

        # The iterated bridge is M(bridge): a segment of the unstable manifold one
        # orbit step forward of the parent's branch. Recording that key -- before the
        # image is indexed -- is what lets every crossing detected on it carry its
        # unstable branch identity (and lets its children advance the key again).
        image_key = bridge.fixed_point.advance_key(bridge.manifold_key, 1)
        iterated.manifold_key = image_key

        # 2. register with the tangle. The bridge is unstable manifold, so every
        #    curve it can legitimately cross is stable and already in the rtree
        #    from compute_intersections -- querying finds all its crossings, and
        #    skipping the rtree inserts keeps repeated blasting fast.
        self.Tangle.add_manifold(iterated, index_segments=False)

        # 3. resolve only new crossings involving the iterated bridge
        new_intersections = self.Tangle.populate_intersections_for_manifold(iterated)

        for ix in new_intersections:
            # A crossing that collides with one already registered IS that crossing;
            # take the registry's id for it so the cut below records the canonical id.
            ix.id = self._intersection_registry.add(ix)

        # 4. cut at crossings. Only the crossings just born on this image are cut at:
        #    they alone carry bracketing points that lie on the image's own polyline.
        new_bridges = (
            self.Tangle.create_bridges(new_intersections, for_manifold=iterated)
            if new_intersections
            else []
        )

        # A bridge is bounded by two crossings, so create_bridges yields nothing when
        # the image crosses the stable manifold fewer than twice -- either zero times
        # (no crossing) or exactly once (a single new homoclinic point cannot bound a
        # bridge on its own). Either way the forward image is still a real piece of
        # unstable manifold whose dynamics must be carried forward, so keep it as one
        # PARTIAL bridge (both endpoint ids stay None) rather than dropping it. (Any
        # new crossing was already registered above; it is picked up as a bridge
        # boundary on later iterates, once a second crossing appears alongside it.)
        # ManifoldMachine.iterate_bridge may return a BaseManifold, so wrap it as a
        # Bridge -- downstream consumers (uniiterated_bridges, genealogy) need the
        # Bridge attributes.
        if not new_bridges:
            if isinstance(iterated, Bridge):
                new_bridges = [iterated]
            else:
                new_bridges = [
                    Bridge(
                        root=iterated.root,
                        stability=iterated.stability,
                        stretch_param=iterated.stretch_param,
                        fixed_point=iterated.fixed_point,
                        tail=iterated.tail,
                        branch_index=iterated.branch_index,
                        manifold_key=image_key,
                    )
                ]

        # 5. record the forward iterate of the parent's two endpoint crossings
        self._register_endpoint_iterates(bridge, new_intersections)

        # 6. Single-copy invariant: a bridge is uniquely defined by the two
        #    intersections it connects, so a freshly cut child whose BridgeId is
        #    already registered IS that bridge -- _register_bridge hands back the
        #    stored object instead of a duplicate. The caller still sees every piece
        #    of the image; it just points at the one stored copy. Iterating the
        #    fixed-point bridge, for instance, re-traces curve already held, so its
        #    children all resolve to existing bridges and nothing is added.
        children = [self._register_bridge(child) for child in new_bridges]
        bridge.iterated = True
        # Mapping forward can refine (and cut at) curve that already exists --
        # the image of a bridge whose points all have iterates is laid over
        # grown manifold -- so every memoised walk is dropped, not just the
        # children's.
        self._bump_generation(invalidate_walks=True)

        return children

    def infer_iterate_table(self, cdist_rtol: float = 0.05) -> int:
        """
        Scan all iterated bridges and record the n=1 forward iterate relationship for
        each boundary intersection.

        Bridge topology identifies *which* intersections to process (only the two
        endpoints of each iterated bridge). The image f(i_src) is then identified by
        canonical distance and branch identity rather than phase-space coordinates,
        which drift under the nonlinear map and are only approximate. Under one
        application of M the image lies on the branches one orbit step forward of the
        source (see FixedPoint.advance_key), with the unstable cdist stretched and
        the stable cdist contracted by the per-map-step factor
        FixedPoint.per_step_beta.

        The stable branch (manifold_b_key) is the reliable discriminator: the stable
        manifold is always indexed, so every intersection carries it. The unstable
        branch (manifold_a_key) is missing on intersections born from an iterated
        bridge (those bridges are not indexed manifolds), so it is used only as an
        extra constraint when both source and candidate have it. Stable branch plus
        canonical distance resolves the arc-length ambiguity that motivated the old
        coordinate-based match.

        Args:
            cdist_rtol: Maximum relative canonical-distance error between the predicted
                image and a candidate (on the correct stable branch) for the match to
                be accepted. Defaults to 0.05, comfortably above per-step scaling noise.

        Returns:
            Number of new iterate relationships recorded.
        """
        registry = self._intersection_registry
        recorded = 0

        for bridge in self._bridges.values():
            if not bridge.iterated:
                continue
            for src_id in (bridge.first_intersection, bridge.second_intersection):
                if src_id is None:
                    continue
                if self._register_forward_iterate(src_id, registry, cdist_rtol):
                    recorded += 1

        return recorded

    def infer_iterates(self, cdist_rtol: float = 0.05) -> int:
        """
        Fill the iterate table for *every* intersection by canonical-distance mapping.

        For each registered intersection this records its n=1 forward image M(i) — the
        intersection one orbit step forward, found by predicting the image's branches
        (:meth:`FixedPoint.advance_key`) and canonical distances (each scaled by its
        own branch's :meth:`FixedPoint.per_step_beta`) and matching against the
        registry. This generalizes :meth:`infer_iterate_table`,
        which records the same relationship but only for bridge-boundary intersections.

        :meth:`compute_intersections` calls this automatically, so the table is dense
        as soon as intersections are computed; iterating bridges then keeps it current
        through the bridge machinery. Idempotent: it skips intersections that already
        have an n=1 entry, so repeated calls are cheap.

        Args:
            cdist_rtol: Maximum relative canonical-distance error for a match (default
                0.05, comfortably above per-step scaling noise).

        Returns:
            Number of new iterate relationships recorded.
        """
        registry = self._intersection_registry
        recorded = 0
        for src_id in registry.all_ids():
            if self._register_forward_iterate(src_id, registry, cdist_rtol):
                recorded += 1
        return recorded

    def _register_forward_iterate(
        self, src_id: int, registry: IntersectionRegistry, cdist_rtol: float
    ) -> bool:
        """
        Record the n=1 forward iterate of one intersection, if it can be identified.

        The image f(src) lies on the branches one map step forward (see
        :meth:`FixedPoint.advance_key`), with each canonical distance scaled by its
        own side's per-map-step factor (:meth:`FixedPoint.per_step_beta`): the
        unstable cdist stretched on the unstable branch's fixed point, the stable
        cdist contracted on the stable branch's. The two coincide on a homoclinic
        crossing; on a HETEROCLINIC one the two branches belong to different
        periodic orbits with different k_values, and using one point's factor for
        both cdists is simply wrong.

        The image is matched by stable branch (always present) plus canonical
        distance, with the unstable branch as an extra constraint when both source
        and candidate carry it (intersections born from iterated bridges have no
        unstable branch). Coordinates are deliberately not used — they drift under
        the nonlinear map, whereas canonical distances and branch keys are exact.

        Returns:
            True if a new iterate edge was recorded, else False.
        """
        if (src_id, 1) in registry.iterate_table:
            return False

        source = registry[src_id]
        if self._register_anchor_self_iterate(src_id, source, registry):
            return True

        # The image's stable branch and stable canonical distance are both
        # predicted, so the candidates are one bisected slice of that branch's
        # own cdist order rather than the whole registry. A crossing outside the
        # window fails the relative stable-cdist test by construction, so the
        # match this narrows to is the same one the full scan would return.
        prediction = self._image_prediction(source)
        if prediction is None:
            return False
        _a_key, b_key, _u_pred, s_pred = prediction
        # The same tolerance _match_image applies below, read once: the window
        # is exactly the set of candidates whose relative stable-cdist error can
        # still come in under cdist_rtol.
        cdist_tol = self._intersection_registry.cdist_tol
        window = cdist_rtol * (abs(s_pred) + cdist_tol)
        candidates = registry.candidates_near_stable_cdist(b_key, s_pred, window)

        best_id, _best_err = self._match_image(
            source, candidates, cdist_rtol, exclude={src_id}
        )
        if best_id is None:
            return False

        registry.register_iterate(src_id, 1, best_id)
        return True

    def _register_anchor_self_iterate(
        self, src_id: int, source: Intersection, registry: IntersectionRegistry
    ) -> bool:
        """
        Record ``f(anchor) = anchor`` for a PERIOD-1 orbit's own periodic point.

        The periodic point sits on both its manifolds at canonical distance
        ``(0, 0)`` and is registered as a crossing there. It is a fixed point of
        the map, so on a period-1 orbit -- and only there -- it is its own forward
        image. That is the one crossing :meth:`_match_image` can never find: the
        candidate it needs is the source, which the "an image is not its source"
        guard rightly excludes for every other crossing (at cdist ``c > 0`` the
        prediction is ``beta * c``, which is nowhere near ``c``).

        On a period > 1 orbit the anchor's image is the anchor of the NEXT orbit
        branch, a different crossing at its own ``(0, 0)``, which the normal match
        finds; the advanced-key test below is what tells the two cases apart.

        Returns:
            True if the self-iterate was recorded, else False.
        """
        tol = registry.cdist_tol
        if abs(source.unstable_cdist) > tol or abs(source.stable_cdist) > tol:
            return False

        prediction = self._image_prediction(source)
        if prediction is None:
            return False

        a_key, b_key, _u_pred, _s_pred = prediction
        if b_key != source.manifold_b_key:
            return False
        if a_key is not None and a_key != source.manifold_a_key:
            return False

        registry.register_iterate(src_id, 1, src_id)
        return True

    @staticmethod
    def _image_prediction(
        source: Intersection,
    ) -> Optional[tuple[Optional[tuple], tuple, float, float]]:
        """
        What the n=1 image of ``source`` must look like: its two branch keys and
        its two canonical distances.

        Each side advances on -- and is scaled by -- its OWN fixed point: on a
        heteroclinic crossing the unstable branch belongs to a different periodic
        orbit than the stable one, with a different ``k_value``, so using one
        point's per-step factor for both cdists is simply wrong.

        Returns:
            ``(advanced unstable key or None, advanced stable key, predicted
            unstable cdist, predicted stable cdist)``, or None when the crossing
            has no stable branch or its fixed point carries no eigendata.
        """
        if source.manifold_b_key is None:
            return None
        stable_fp = source.manifold_b_key[0]
        if (
            not getattr(stable_fp, "unstable_eigenvalues", None)
            or stable_fp.k_value is None
        ):
            return None

        b_key = stable_fp.advance_key(source.manifold_b_key, 1)
        a_key = (
            source.manifold_a_key[0].advance_key(source.manifold_a_key, 1)
            if source.manifold_a_key is not None
            else None
        )
        unstable_fp = (
            stable_fp if source.manifold_a_key is None else source.manifold_a_key[0]
        )
        u_pred = source.unstable_cdist * unstable_fp.per_step_beta("unstable")
        s_pred = source.stable_cdist * stable_fp.per_step_beta("stable")
        return a_key, b_key, u_pred, s_pred

    def _match_image(
        self,
        source: Intersection,
        candidates: Iterable[tuple[int, Intersection]],
        cdist_rtol: float,
        *,
        exclude: "set[int] | frozenset[int]" = frozenset(),
        accept: Optional[Callable[[Intersection], bool]] = None,
    ) -> tuple[Optional[int], float]:
        """
        Pick the crossing among ``candidates`` that is the n=1 image of ``source``.

        The single place the branch/beta prediction of an iterate is turned into a
        match, shared by the registry-wide heuristic
        (:meth:`_register_forward_iterate`) and the explicit bridge-endpoint
        registration (:meth:`_register_endpoint_iterates`). A candidate must sit on
        the advanced STABLE branch -- canonical distance restarts at 0 on every
        branch, so on a period > 1 or heteroclinic tangle the cdist window alone
        would happily pick a crossing on the wrong branch -- and, when both sides
        carry one, on the advanced unstable branch too. The two canonical distances
        are then compared INDIVIDUALLY; their product is never used, since two
        different iterate chains share it (CLAUDE.md, area preservation).

        Args:
            source: The crossing whose forward image is wanted.
            candidates: ``(id, Intersection)`` pairs to choose from.
            cdist_rtol: Maximum relative canonical-distance error for a match.
            exclude: Ids that may not be chosen (the source itself, and any image
                already claimed by another source in the same call).
            accept: Optional extra predicate a candidate must satisfy -- used to
                gate on the real map's image coordinates.

        Returns:
            ``(matched id or None, best relative error seen)``.
        """
        prediction = self._image_prediction(source)
        if prediction is None:
            return None, float("inf")

        a_key, b_key, u_pred, s_pred = prediction
        cdist_tol = self._intersection_registry.cdist_tol

        best_id, best_err = None, float("inf")
        for tgt_id, tgt in candidates:
            if tgt_id is None or tgt_id in exclude:
                continue
            if tgt.manifold_b_key != b_key:
                continue
            if (
                a_key is not None
                and tgt.manifold_a_key is not None
                and tgt.manifold_a_key != a_key
            ):
                continue
            if accept is not None and not accept(tgt):
                continue
            u_rel = abs(tgt.unstable_cdist - u_pred) / (abs(u_pred) + cdist_tol)
            s_rel = abs(tgt.stable_cdist - s_pred) / (abs(s_pred) + cdist_tol)
            err = max(u_rel, s_rel)
            if err < best_err:
                best_err, best_id = err, tgt_id

        if best_id is not None and best_err <= cdist_rtol:
            return best_id, best_err
        return None, best_err

    def _existing_image_bridges(self, bridge: Bridge) -> Optional[list[Bridge]]:
        """
        The already-computed bridges that tile this bridge's forward image, if the
        image is entirely covered by them.

        A bridge spanning crossings ``a -> b`` images onto the arc
        ``f(a) -> f(b)`` of the branch one map step forward
        (:meth:`FixedPoint.advance_key`). If that arc has already been grown and
        cut, those bridges ARE the forward image -- there is no need (and it is
        actively harmful) to re-map and re-cut it: re-mapping lays a second,
        slightly different polyline over curve that already exists (the "zig-zag"),
        spawns overlapping duplicate bridges, and trips the same-stability
        (unstable x unstable) detector along the near-coincident pair.

        The image is only reusable when the existing bridges TILE it exactly: the
        first starts at ``f(a)``, the last ends at ``f(b)``, and each shares an
        endpoint with the next. Anything less means part of the image is curve that
        has not been computed, which is precisely the case that must fall through
        to real iteration.

        Args:
            bridge: The bridge about to be iterated.

        Returns:
            The existing image bridges in unstable cdist order, or ``None`` to fall
            through to real iteration.
        """
        span = self._image_span(bridge, 1)
        if span is None:
            logger.debug(
                "bridge %s: no reusable image -- an endpoint has no registered n=1 "
                "iterate, so the image will be recomputed by mapping it forward",
                bridge.id,
            )
            return None
        (lo_image, hi_image), _cdists, _key = span

        image_ids = self._image_bridge_ids(bridge, 1)
        if not image_ids:
            logger.debug(
                "bridge %s: image arc %s -> %s carries no registered bridge; "
                "recomputing",
                bridge.id,
                lo_image,
                hi_image,
            )
            return None

        # The chain must run from f(a) to f(b) with no gap: the endpoint ids come
        # from _image_span already ordered by cdist, as image_ids is.
        if image_ids[0][0] != lo_image or image_ids[-1][1] != hi_image:
            logger.debug(
                "bridge %s: image arc %s -> %s is only partly cut (%s); recomputing",
                bridge.id,
                lo_image,
                hi_image,
                image_ids,
            )
            return None
        if any(lhs[1] != rhs[0] for lhs, rhs in zip(image_ids, image_ids[1:])):
            logger.debug(
                "bridge %s: image bridges %s do not tile the arc; recomputing",
                bridge.id,
                image_ids,
            )
            return None

        return [self._bridges[bid] for bid in image_ids]


    def _register_endpoint_iterates(
        self,
        bridge: Bridge,
        image_crossings: list[Intersection],
        cdist_rtol: float = 0.05,
    ) -> int:
        """
        Record the n=1 forward iterate of a bridge's two endpoint crossings.

        The image of a bridge's endpoint crossing is a crossing of the image bridge,
        so it is one of the crossings just born on that image -- no registry-wide
        search and no heuristic is needed. The match is made by
        :meth:`_match_image`: the advanced stable (and, when present, unstable)
        branch key, then both canonical distances individually, each scaled by its
        own branch's :meth:`FixedPoint.per_step_beta`.

        Two further guards apply here, where the registry-wide heuristic cannot use
        them:

        * The REAL map is a sanity gate. ``f(source.coords)`` is where the image
          must be, so a candidate further than ``1e-3 * max(1, |f|)`` from it is
          rejected outright. That both catches a cdist coincidence and stops a false
          registration when the true image lies beyond the computed manifold and so
          is not among the crossings at all.
        * The map is injective, so the two endpoints of one bridge cannot share an
          image: a candidate claimed by the first endpoint is off limits to the
          second.

        Args:
            bridge: The parent bridge that was just iterated.
            image_crossings: The crossings born on its forward image, ids assigned.
            cdist_rtol: Maximum relative canonical-distance error for a match.

        Returns:
            Number of new iterate relationships recorded.
        """
        if not image_crossings:
            return 0

        registry = self._intersection_registry
        candidates = [(ix.id, ix) for ix in image_crossings if ix.id is not None]
        claimed: set[int] = set()
        recorded = 0

        for src_id in (bridge.first_intersection, bridge.second_intersection):
            if src_id is None or src_id not in registry:
                continue
            if (src_id, 1) in registry.iterate_table:
                continue

            source = registry[src_id]
            if source.coords is None:
                continue

            image = np.asarray(
                self.dynamical_system.map(np.asarray(source.coords, dtype=float)),
                dtype=float,
            ).ravel()
            gate = 1e-3 * max(1.0, float(np.linalg.norm(image)))

            # image/gate bound as defaults: the predicate is consumed inside this
            # iteration, and binding makes that independent of the loop variable.
            def near_the_true_image(
                candidate: Intersection,
                image: NDArray[np.float64] = image,
                gate: float = gate,
            ) -> bool:
                if candidate.coords is None:
                    return False
                offset = np.asarray(candidate.coords, dtype=float) - image
                return float(np.linalg.norm(offset)) <= gate

            best_id, best_err = self._match_image(
                source,
                candidates,
                cdist_rtol,
                exclude=claimed | {src_id},
                accept=near_the_true_image,
            )

            if best_id is None:
                logger.debug(
                    "no image of crossing %s among the %d crossings born on the "
                    "image of its bridge (best relative cdist error %.3g, %d "
                    "candidate(s) already claimed)",
                    src_id,
                    len(candidates),
                    best_err,
                    len(claimed),
                )
                continue

            registry.register_iterate(src_id, 1, best_id)
            claimed.add(best_id)
            recorded += 1

        return recorded

    def build_intersection_graph(self) -> nx.MultiDiGraph:
        """
        The registry's intersection graph, decorated for this workbench.

        The graph itself -- nodes, per-branch adjacency, iterate edges -- belongs
        to the :class:`IntersectionRegistry` (plan 2.7); this method supplies the
        one thing the registry cannot know, the bridges, and returns a COPY so
        callers may decorate or lay out the result without touching the live one.

        Unstable adjacency is exactly the registered bridges: a bridge IS the
        piece of unstable manifold between two consecutive crossings, so
        disconnected pieces of unstable manifold (an original bridge and its
        iterated children) stay separate paths instead of merging into one chain.
        Stable adjacency joins consecutive crossings ON EACH STABLE BRANCH;
        canonical distances on different branches are measured from different
        anchors and are not comparable.

        Returns:
            nx.MultiDiGraph: A copy of the registry graph.
        """
        bridges = list(self._bridges)
        return self._intersection_registry.graph(bridges=bridges).copy()

    def iterate_all_bridges(self) -> list[Bridge]:
        """
        Iterate all bridges that have not yet been mapped forward.

        Returns:
            All new bridges produced across all iterations.
        """
        pending = list(
            self.uniiterated_bridges
        )  # snapshot before loop mutates _bridges

        all_new: list[Bridge] = []
        for bridge in pending:
            all_new.extend(self.iterate_bridge(bridge))

        return all_new

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
        Visualize the intersection graph with edges colored by type and stability.

        Args:
            G: The intersection graph from build_intersection_graph().
            layout: Layout algorithm. ``"auto"`` (default) picks
                ``"kamada_kawai"`` for ≤ 8 nodes and ``"stable_linear"``
                otherwise. Explicit options:

                * ``"stable_linear"`` — nodes sorted by stable arc-length on a
                  horizontal line; bridges arch above and stable edges run flat.
                  Usually the clearest layout for tangle graphs.
                * ``"unstable_linear"`` — same but sorted by unstable arc-length.
                * ``"cdist"`` — node position = (unstable_cdist, stable_cdist),
                  revealing the full arc-length structure on labelled axes.
                * ``"kamada_kawai"``, ``"spring"``, ``"circular"``, ``"spectral"``
                  — standard networkx force/geometric layouts.
            figsize: Figure size as (width, height).
            display_mode: One of ``"auto"``, ``"full"``, or ``"compact"``.
                ``"auto"`` switches to compact when the node count exceeds
                *compact_threshold*. ``"full"`` uses large, labeled nodes.
                ``"compact"`` uses small dots suitable for dense graphs.
            compact_threshold: Node count above which ``"auto"`` picks compact.
            node_size: Override the node area in points². Defaults to 800 in
                full mode and 80 in compact mode.
            label_mode: Content of node labels. One of:
                ``"id"`` — intersection ID only (default);
                ``"coords"`` — (x, y) phase-space coordinates;
                ``"cdist"`` — unstable and stable arc-lengths;
                ``"all"`` — ID, coordinates, and both cdists;
                ``"none"`` — no labels.
            node_color_by: How to color the nodes. One of:
                ``"none"`` — white (default);
                ``"unstable_cdist"`` — viridis colormap by unstable arc-length;
                ``"stable_cdist"`` — plasma colormap by stable arc-length;
                ``"fixed_point"`` — distinct color per originating fixed point.
            show_iterate_edges: Whether to draw iterate-type edges (rendered
                dashed in purple to distinguish them from adjacency edges).
            save_path: Optional file path to save the figure.

        Returns:
            (fig, ax) matplotlib Figure and Axes.
        """
        import matplotlib.colors as mcolors
        from matplotlib.patches import FancyArrowPatch
        from matplotlib.lines import Line2D

        if G.number_of_nodes() == 0:
            logger.warning("Graph has no nodes to visualize")
            return None, None

        n_nodes = G.number_of_nodes()

        mode = (
            ("compact" if n_nodes > compact_threshold else "full")
            if display_mode == "auto"
            else display_mode
        )

        if node_size is None:
            node_size = 80 if mode == "compact" else 800

        # ── Layout ────────────────────────────────────────────────────────────
        if layout == "auto":
            layout = "kamada_kawai" if n_nodes <= 8 else "stable_linear"

        if layout == "stable_linear":
            sorted_nodes = sorted(
                G.nodes(), key=lambda nd: G.nodes[nd].get("stable_cdist") or 0.0
            )
            pos = {nd: (float(i), 0.0) for i, nd in enumerate(sorted_nodes)}
        elif layout == "unstable_linear":
            sorted_nodes = sorted(
                G.nodes(), key=lambda nd: G.nodes[nd].get("unstable_cdist") or 0.0
            )
            pos = {nd: (float(i), 0.0) for i, nd in enumerate(sorted_nodes)}
        elif layout == "cdist":
            pos = {
                node: (
                    G.nodes[node].get("unstable_cdist") or 0.0,
                    G.nodes[node].get("stable_cdist") or 0.0,
                )
                for node in G.nodes()
            }
        elif layout == "spring":
            pos = nx.spring_layout(G, k=1, iterations=50, seed=42)
        elif layout == "circular":
            pos = nx.circular_layout(G)
        elif layout == "kamada_kawai":
            pos = nx.kamada_kawai_layout(G)
        elif layout == "spectral":
            pos = nx.spectral_layout(G)
        else:
            pos = nx.spring_layout(G, seed=42)

        fig, ax = plt.subplots(figsize=figsize)

        # ── Node colors ───────────────────────────────────────────────────────
        cmap_obj = None
        norm_obj = None
        if node_color_by == "unstable_cdist":
            values = [G.nodes[nd].get("unstable_cdist") or 0.0 for nd in G.nodes()]
            cmap_obj = cm.viridis
            norm_obj = mcolors.Normalize(vmin=min(values), vmax=max(values))
            node_colors = [cmap_obj(norm_obj(v)) for v in values]
        elif node_color_by == "stable_cdist":
            values = [G.nodes[nd].get("stable_cdist") or 0.0 for nd in G.nodes()]
            cmap_obj = cm.plasma
            norm_obj = mcolors.Normalize(vmin=min(values), vmax=max(values))
            node_colors = [cmap_obj(norm_obj(v)) for v in values]
        elif node_color_by == "fixed_point":
            fps = list(
                dict.fromkeys(
                    G.nodes[nd].get("manifold_a_key", (None,))[0] for nd in G.nodes()
                )
            )
            fp_idx = {fp: i for i, fp in enumerate(fps)}
            fp_cmap = cm.get_cmap("Set1", max(len(fps), 1))
            node_colors = [
                fp_cmap(fp_idx.get(G.nodes[nd].get("manifold_a_key", (None,))[0], 0))
                for nd in G.nodes()
            ]
        else:
            node_colors = ["white"] * n_nodes

        # ── Draw nodes ────────────────────────────────────────────────────────
        nx.draw_networkx_nodes(
            G,
            pos,
            node_color=node_colors,
            edgecolors="black",
            linewidths=1.5 if mode == "compact" else 2.0,
            node_size=node_size,
            ax=ax,
        )

        # ── Draw edges via FancyArrowPatch ────────────────────────────────────
        # Shrink endpoints so arrows touch the node boundary, not the centre.
        shrink = np.sqrt(node_size / np.pi)
        alpha = 0.55 if mode == "compact" else 0.80
        mutation = 10 if mode == "compact" else 18
        lw_adj = 1.5 if mode == "compact" else 2.0
        lw_iter = 1.2 if mode == "compact" else 1.8

        _EDGE_STYLE: dict[tuple[str, str], dict] = {
            ("adjacency", "unstable"): {
                "color": "#3b82f6",
                "lw": lw_adj,
                "base_rad": 0.20,
                "ls": "solid",
            },
            ("adjacency", "stable"): {
                "color": "#ef4444",
                "lw": lw_adj,
                "base_rad": -0.20,
                "ls": "solid",
            },
            ("iterate", "unstable"): {
                "color": "#a855f7",
                "lw": lw_iter,
                "base_rad": 0.38,
                "ls": "dashed",
            },
        }
        _FALLBACK = _EDGE_STYLE[("adjacency", "unstable")]

        # Track how many edges have been drawn for each (u, v) pair so that
        # parallel edges get staggered curvature and don't overlap.
        _pair_count: dict[tuple, int] = {}

        for u, v, _key, data in G.edges(keys=True, data=True):
            if u == v:
                continue
            edge_type = data.get("type", "adjacency")
            stability = data.get("stability", "unstable")

            if edge_type == "iterate" and not show_iterate_edges:
                continue

            style = _EDGE_STYLE.get((edge_type, stability), _FALLBACK)

            pair = (u, v)
            idx = _pair_count.get(pair, 0)
            _pair_count[pair] = idx + 1
            rad = style["base_rad"] + idx * 0.15 * np.sign(style["base_rad"] or 1)

            patch = FancyArrowPatch(
                posA=pos[u],
                posB=pos[v],
                connectionstyle=f"arc3,rad={rad}",
                arrowstyle="-|>",
                color=style["color"],
                linewidth=style["lw"],
                linestyle=style["ls"],
                alpha=alpha,
                mutation_scale=mutation,
                shrinkA=shrink,
                shrinkB=shrink,
                transform=ax.transData,
                zorder=2,
            )
            ax.add_patch(patch)

        # ── Labels ────────────────────────────────────────────────────────────
        if label_mode != "none":
            labels: dict[int, str] = {}
            for node in G.nodes():
                data = G.nodes[node]
                coords = data.get("coords")
                u_cd = data.get("unstable_cdist")
                s_cd = data.get("stable_cdist")
                parts: list[str] = []

                if label_mode in ("id", "all"):
                    parts.append(str(node))
                if label_mode in ("coords", "all") and coords is not None:
                    parts.append(f"({coords[0]:.2f},{coords[1]:.2f})")
                if label_mode in ("cdist", "all"):
                    u_str = f"{u_cd:.2f}" if u_cd is not None else "?"
                    s_str = f"{s_cd:.2f}" if s_cd is not None else "?"
                    parts.append(f"u:{u_str}\ns:{s_str}")

                labels[node] = "\n".join(parts) if parts else str(node)

            font_size = 5 if mode == "compact" else 7
            nx.draw_networkx_labels(
                G, pos, labels, font_size=font_size, font_weight="bold", ax=ax
            )

        # ── Colorbar ──────────────────────────────────────────────────────────
        if cmap_obj is not None and norm_obj is not None:
            sm = cm.ScalarMappable(cmap=cmap_obj, norm=norm_obj)
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
            cbar.set_label(
                "Unstable arc-length"
                if node_color_by == "unstable_cdist"
                else "Stable arc-length",
                fontsize=9,
            )

        # ── Axis appearance ───────────────────────────────────────────────────
        _linear_layouts = {"cdist", "stable_linear", "unstable_linear"}
        if layout in _linear_layouts:
            ax.set_axis_on()
            if layout == "stable_linear":
                ax.set_xlabel("Stable arc-length order", fontsize=10)
                ax.set_yticks([])
            elif layout == "unstable_linear":
                ax.set_xlabel("Unstable arc-length order", fontsize=10)
                ax.set_yticks([])
            else:
                ax.set_xlabel("Unstable arc-length (cdist)", fontsize=10)
                ax.set_ylabel("Stable arc-length (cdist)", fontsize=10)
                ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
            ax.margins(0.15)
        else:
            ax.axis("off")

        # ── Legend ────────────────────────────────────────────────────────────
        legend_handles = [
            Line2D([0], [0], color="#3b82f6", linewidth=2, label="Unstable adjacency"),
            Line2D([0], [0], color="#ef4444", linewidth=2, label="Stable adjacency"),
        ]
        if show_iterate_edges:
            legend_handles.append(
                Line2D(
                    [0],
                    [0],
                    color="#a855f7",
                    linewidth=1.5,
                    linestyle="--",
                    label="Iterate",
                )
            )
        ax.legend(handles=legend_handles, loc="upper right", fontsize=9)

        mode_label = f"{n_nodes} nodes · {mode} mode"
        ax.set_title(
            f"Intersection Graph  ({mode_label})", fontsize=12, fontweight="bold"
        )

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches="tight")

        plt.show()
        return fig, ax

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

    def trim_stable_manifolds(self, fixed_point: FixedPoint):
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

    def _iter_manifolds(self, fp, stability: Stability | None = None) -> Iterable:
        """Yield all manifolds for a fixed point (optionally filter by stability)."""
        for (kfp, kstab, _oi, _bi), M in self.manifolds.items():
            if kfp is fp and (stability is None or kstab == stability):
                yield M

