import logging
from typing import Callable, Literal, Iterable, Optional
import numpy.typing as npt
from typing_extensions import Annotated

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
from .Bridge import Bridge
from .Intersection import Intersection
from .IntersectionRegistry import IntersectionRegistry

Stability = Literal["unstable", "stable"]

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


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
        self.manifolds: dict[tuple[FixedPoint, Stability, int, int], BaseManifold] = {}

        self._intersection_registry = IntersectionRegistry()
        self._bridges: list[Bridge] = []

    def construct_fixed_point(self, initial_guess) -> FixedPoint:
        """
        Constructs a fixed point for a given initial guess.
        Adds that fixed point to the class storage.

        Args:
            initial_guess (_type_): _description_

        Returns:
            FixedPoint: _description_
        """

        # TODO go through the fixed point classes and eliminate this 2
        fixed_point = self._fp_solver.construct_fixed_point(initial_guess, 2)

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

        initial_segments = self._man_maker.construct_kevin_way(
            fixed_point, stability, num_branches
        )

        for (orbit_index, branch_index), manifold in initial_segments.items():

            key = (fixed_point, stability, orbit_index, branch_index)
            manifold.manifold_key = key  # ← new
            self.manifolds[key] = manifold

        return initial_segments

    def initialize_both_manifolds(self, fixed_point: FixedPoint, num_branches: int = 1):

        unstable_segments = self.initialize_manifold(
            fixed_point, "unstable", num_branches
        )
        stable_segments = self.initialize_manifold(fixed_point, "stable", num_branches)

        # for (orbit_index, branch_index), manifold in unstable_segments.items():

        #     self.manifolds[(fixed_point, "unstable", orbit_index, branch_index)] = (
        #         manifold
        #     )

        # for (orbit_index, branch_index), manifold in stable_segments.items():

        #     self.manifolds[(fixed_point, "stable", orbit_index, branch_index)] = (
        #         manifold
        #     )

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

        # self._close_iterate_gaps(fixed_point, stability)

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
        branch_index: int = 0,
    ):

        if self.manifolds.get((fixed_point, stability, 0, branch_index)) is None:
            raise ValueError(f"""Manifold for fixed point {fixed_point} 
                    with stability {stability} has not been initialized.
                    Please run initialize_manifold first.""")

        # TODO change this so it uses the actual arclength
        current_distance = self.manifolds.get(
            (fixed_point, stability, 0, branch_index)
        ).tail.cdist

        while current_distance < length:

            self.grow_n_times(
                fixed_point, stability, num_iterations=1, branch_index=branch_index
            )

            current_distance = self.manifolds.get(
                (fixed_point, stability, 0, branch_index)
            ).tail.cdist

        else:
            return None

    def grown_until_intersection(
        self, fixed_point: FixedPoint, stability: Stability, max_iterations: int = 10
    ):

        if self.manifolds.get((fixed_point, stability, 0, 0)) is None:
            raise ValueError(f"""Manifold for fixed point {fixed_point} 
                    with stability {stability} has not been initialized.
                    Please run initialize_manifold first.""")

        self.compute_intersections(fixed_point, infer_iterates=False)

        num_initial_intersections = len(self.Tangle._intersecting_segments)
        print(f"current intersections {num_initial_intersections}")

        for _ in range(max_iterations):

            self.grow_n_times(fixed_point, stability, num_iterations=1)

            self.compute_intersections(fixed_point, infer_iterates=False)
            num_current_intersections = len(self.Tangle._intersecting_segments)

            if num_current_intersections > num_initial_intersections:
                return None

        else:
            raise ValueError("Max iterations reached, no intersection found")

    def index_manifolds(self, fixed_point: FixedPoint, stability: Stability | None):

        count = 0
        for M in self._iter_manifolds(fixed_point, stability):
            self.Tangle.add_manifold(M)
            count += 1

        return self

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

        for fp in fixed_points:
            self.index_manifolds(fp, "unstable")
            self.index_manifolds(fp, "stable")

        self.Tangle.populate_intersection_dict()

        for intersection in self.Tangle._intersections:
            self._intersection_registry.add(intersection)

        if reset and preserve_ids and len(old_registry) > 0:
            self._intersection_registry.reindex_from(old_registry)

        if infer_iterates:
            self.infer_iterates()

        return self.Tangle.iter_intersection_coords()

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

    def create_bridges(self, fixed_point: Optional[FixedPoint] = None):
        """
        Cut indexed unstable manifolds into bridges.

        Args:
            fixed_point: If given, only build bridges for that fixed point's
                unstable manifolds. If None, build bridges for every indexed
                unstable manifold at once. Each bridge keeps its own fixed_point
                linkage either way.

        Returns:
            The newly created bridges.
        """
        bridges = self.Tangle.create_bridges(fixed_point=fixed_point)
        self._bridges.extend(bridges)
        self._assign_bridge_intersections(bridges)
        return bridges

    def create_resonance_zone(self, fixed_point: FixedPoint):

        # grow stable manifold until it turns around
        self.grow_until_turnaround(fixed_point, "stable")

        # grow unstable manifold until it intersects the stable
        self.grown_until_intersection(fixed_point, "unstable")

        # compute all intersections (possibly redundant)
        intersections = self.compute_intersections(fixed_point)

        # cut the unstable manifold into a bridge
        bridges = self.create_bridges(fixed_point)

        # trim stable manifold
        # self.trim_stable_at_first_intersection(fixed_point)

        self.plot_tangle(fixed_point, "stable", color="r")
        self.plot_all_bridges(bridges)
        self.plot_intersections(fixed_point)

    @property
    def uniiterated_bridges(self) -> list[Bridge]:
        """All bridges that have not yet been iterated forward."""
        return [b for b in self._bridges if not b.iterated]

    @property
    def intersection_registry(self) -> IntersectionRegistry:
        return self._intersection_registry

    def iterate_bridge(self, bridge: Bridge) -> list[Bridge]:
        """
        Map a bridge forward one iterate, add the result to the tangle, detect new
        intersections with the stable manifold, cut the result into new bridges,
        and return those bridges.

        Marks the original bridge as iterated and wires parent/child links.

        Args:
            bridge: A bridge created by create_bridges() or a previous iterate_bridge().

        Returns:
            List of new Bridge objects from cutting the iterated result.
            If the iterated bridge makes no new crossings, returns a single-element
            list containing the unsplit iterated bridge.

        Raises:
            ValueError: If bridge has already been iterated.
            ValueError: If create_bridges() has not been called yet.
        """
        if bridge.iterated:
            raise ValueError(
                "This bridge has already been iterated. Check bridge.children for the results."
            )
        if not self._bridges:
            raise ValueError(
                "No bridges registered. Call create_bridges() before iterate_bridge()."
            )

        # 1. map forward
        iterated = self._man_machine.iterate_bridge(bridge)

        # The iterated bridge is M(bridge): a segment of the unstable manifold one
        # orbit step forward of the parent's branch. Recording that key lets the
        # crossings detected on it carry their unstable branch identity (and lets
        # those bridges advance the key again if iterated further).
        if bridge.manifold_key is not None:
            iterated.manifold_key = self._advance_key_forward(
                bridge.manifold_key, bridge.fixed_point
            )

        # 2. add to tangle (stable manifold already indexed from compute_intersections)
        self.Tangle.add_manifold(iterated)

        # 3. resolve only new crossings involving the iterated bridge
        new_intersections = self.Tangle.populate_intersections_for_manifold(iterated)

        for ix in new_intersections:
            self._intersection_registry.add(ix)

        # 4. cut at crossings
        if new_intersections:
            new_bridges = self.Tangle.create_bridges(for_manifold=iterated)
        else:
            # No crossings: the iterated manifold is itself one unsplit bridge.
            # ManifoldMachine.iterate_bridge returns a BaseManifold, so wrap it as
            # a Bridge (carrying its manifold_key) — downstream consumers
            # (uniiterated_bridges, infer_iterate_table, genealogy) require Bridge
            # attributes such as .iterated / .children. This branch is only reached
            # when an iterated bridge happens not to cross the stable manifold.
            if isinstance(iterated, Bridge):
                new_bridges = [iterated]
            else:
                wrapped = Bridge(
                    root=iterated.root,
                    stability=iterated.stability,
                    stretch_param=iterated.stretch_param,
                    fixed_point=iterated.fixed_point,
                    tail=iterated.tail,
                    branch_index=iterated.branch_index,
                )
                wrapped.manifold_key = getattr(iterated, "manifold_key", None)
                new_bridges = [wrapped]

        # 5. wire genealogy
        bridge.iterated = True
        bridge.children = new_bridges
        for nb in new_bridges:
            nb.parent = bridge

        # 6. register
        self._bridges.extend(new_bridges)
        self._assign_bridge_intersections(new_bridges)

        return new_bridges

    def infer_iterate_table(self, cdist_rtol: float = 0.05) -> int:
        """
        Scan all iterated bridges and record the n=1 forward iterate relationship for
        each boundary intersection.

        Bridge topology identifies *which* intersections to process (only the two
        endpoints of each iterated bridge). The image f(i_src) is then identified by
        canonical distance and branch identity rather than phase-space coordinates,
        which drift under the nonlinear map and are only approximate. Under one
        application of M the image lies on the branches one orbit step forward of the
        source (see _advance_key_forward), with the unstable cdist stretched and the
        stable cdist contracted by the per-step factor lambda_u ** (1 / k_value).

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

        for bridge in self._bridges:
            if not bridge.iterated or not bridge.children:
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
        (``_advance_key_forward``) and canonical distances (unstable stretched by the
        per-step factor ``lambda_u ** (1 / k_value)``, stable contracted by it) and
        matching against the registry. This generalizes :meth:`infer_iterate_table`,
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

        The image f(src) lies on the branches one orbit step forward (see
        ``_advance_key_forward``), with the unstable canonical distance stretched and
        the stable one contracted by the per-step factor ``lambda_u ** (1 / k_value)``.
        The image is matched by stable branch (always present) plus canonical distance,
        with the unstable branch as an extra constraint when both source and candidate
        carry it (intersections born from iterated bridges have no unstable branch).
        Coordinates are deliberately not used — they drift under the nonlinear map,
        whereas canonical distances and branch keys are exact.

        Returns:
            True if a new iterate edge was recorded, else False.
        """
        if (src_id, 1) in registry.iterate_table:
            return False
        src = registry[src_id]
        if src.manifold_b_key is None:
            return False
        fp = src.manifold_b_key[0]
        evals = getattr(fp, "unstable_eigenvalues", None)
        if not evals:
            return False
        lambda_u = float(np.abs(np.asarray(evals[0]).ravel()[0]))
        beta = lambda_u ** (1.0 / fp.k_value)  # per single-map-step cdist factor

        b_key = self._advance_key_forward(src.manifold_b_key, fp)
        a_key = (
            self._advance_key_forward(src.manifold_a_key, fp)
            if src.manifold_a_key is not None
            else None
        )
        u_pred = src.unstable_cdist * beta
        s_pred = src.stable_cdist / beta

        best_id, best_err = None, float("inf")
        for tgt_id, tgt in registry:
            if tgt_id == src_id:
                continue
            if tgt.manifold_b_key != b_key:
                continue
            if (
                a_key is not None
                and tgt.manifold_a_key is not None
                and tgt.manifold_a_key != a_key
            ):
                continue
            u_rel = abs(tgt.unstable_cdist - u_pred) / (abs(u_pred) + registry.cdist_tol)
            s_rel = abs(tgt.stable_cdist - s_pred) / (abs(s_pred) + registry.cdist_tol)
            err = max(u_rel, s_rel)
            if err < best_err:
                best_err, best_id = err, tgt_id

        if best_id is not None and best_err <= cdist_rtol:
            registry.register_iterate(src_id, 1, best_id)
            return True
        return False

    def _advance_key_forward(
        self, key: tuple[FixedPoint, Stability, int, int], fixed_point: FixedPoint
    ) -> tuple[FixedPoint, Stability, int, int]:
        """
        Return the manifold key one forward application of the map advances `key` to.

        Under one application of M the orbit index advances by one (M cycles the
        periodic orbit). For a fixed point with inversion the branch index flips each
        time the orbit index wraps past the last orbit point (a full orbit returns to
        the same point on the opposite eigenvector branch; two full orbits = k_value
        steps return to the start).

        Note:
            The non-inversion path is exercised by the period-1 and period-3 tests;
            the inversion branch-flip is implemented from first principles but has not
            yet been validated against a computed inversion trellis.
        """
        fp, stability, orbit_index, branch_index = key
        period = fixed_point.period
        if orbit_index == period - 1:
            new_orbit_index = 0
            new_branch_index = (
                1 - branch_index if fixed_point.check_inversion() else branch_index
            )
        else:
            new_orbit_index = orbit_index + 1
            new_branch_index = branch_index
        return (fp, stability, new_orbit_index, new_branch_index)

        return recorded

    def populate_registry(self) -> IntersectionRegistry:
        """Rebuild the intersection registry from the current Tangle state."""
        self._intersection_registry = IntersectionRegistry()
        for intersection in self.Tangle._intersections:
            self._intersection_registry.add(intersection)
        return self._intersection_registry

    def _assign_bridge_intersections(self, bridges: list[Bridge]) -> None:
        """
        Populate first_intersection / second_intersection on each bridge by
        matching its root and tail cdists against the registry.

        Bridge root and tail points are created at the same cdist as the
        intersection they flank (both computed as the midpoint of the crossing
        segment), so a nearest-unstable-cdist lookup is exact up to float precision.

        Args:
            bridges: Freshly created Bridge objects whose endpoint fields are unset.
        """
        registry = self._intersection_registry
        for bridge in bridges:
            root_u = bridge.root.get_cdist("unstable")
            tail_u = bridge.tail.get_cdist("unstable")
            bridge.first_intersection = registry.nearest_by_unstable_cdist(root_u)
            bridge.second_intersection = registry.nearest_by_unstable_cdist(tail_u)

    def build_intersection_graph(self) -> nx.MultiDiGraph:
        """
        Build the intersection graph from bridge topology.

        Unstable edges are derived from registered bridges (one edge per bridge,
        connecting the two intersection points it spans), so disconnected pieces
        of the unstable manifold — e.g., an original bridge and its iterated
        children — appear as separate paths rather than one merged chain.

        Stable edges connect consecutive intersections in stable-cdist order,
        which is correct because the stable manifold is a single connected piece.

        Iterate edges come from the iterate table registered by infer_iterate_table().
        """
        registry = self._intersection_registry
        G = nx.MultiDiGraph()

        for ix_id, ix in registry:
            G.add_node(
                ix_id,
                coords=ix.coords,
                unstable_cdist=ix.unstable_cdist,
                stable_cdist=ix.stable_cdist,
                manifold_a_key=ix.manifold_a_key,
                manifold_b_key=ix.manifold_b_key,
                label=ix.label,
            )

        stable_ids = registry.by_stable_cdist
        for i in range(len(stable_ids) - 1):
            u, v = stable_ids[i], stable_ids[i + 1]
            G.add_edge(u, v, key=f"adj_s_{i}", type="adjacency", stability="stable")

        for bridge in self._bridges:
            a = bridge.first_intersection
            b = bridge.second_intersection
            if a is not None and b is not None and a != b and a in G and b in G:
                G.add_edge(a, b, type="adjacency", stability="unstable")

        for src_id, n_to_tgt in registry.iterate_table._forward.items():
            for n, tgt_id in n_to_tgt.items():
                if src_id in G and tgt_id in G:
                    G.add_edge(
                        src_id,
                        tgt_id,
                        key=f"iter_{src_id}_{n}",
                        type="iterate",
                        stability="unstable",
                        n=n,
                    )

        return G

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
            print("Graph has no nodes to visualize")
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
            bridges: List of bridges to plot. Defaults to self._bridges.
        """
        if bridges is None:
            bridges = self._bridges
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
        """

        intersecting_seg_ids = self.Tangle._intersecting_segments
        # get a list of all the intersecting segment ids
        intersecting_seg_ids = [n for pair in intersecting_seg_ids for n in pair]

        for manifold in self._iter_manifolds(fixed_point, "stable"):

            all_segs = list(self.Tangle._manifold_segs[manifold])

            # find all the intersecting segment ids that are on this manifold
            candidate_segs = list(set(all_segs) & set(intersecting_seg_ids))

            # get the segments
            segs = [
                self.Tangle._seg_lookup[candidate_segs[k]]
                for k in range(len(candidate_segs))
            ]

            # find the segment with the largest cdist
            max_seg = max(segs, key=lambda s: s.p0_seg1.get_cdist())

            # set the new truncated tail
            new_tail = max_seg.p0_seg1
            manifold.tail = new_tail

    def _iter_manifolds(self, fp, stability: Stability | None = None) -> Iterable:
        """Yield all manifolds for a fixed point (optionally filter by stability)."""
        for (kfp, kstab, _oi, _bi), M in self.manifolds.items():
            if kfp is fp and (stability is None or kstab == stability):
                yield M

    def _close_iterate_gaps(self, fixed_point: FixedPoint, stability: str) -> None:
        """One extra iterate pass on orbit 0 to wire up dangling tail points."""
        orbit_indices = fixed_point.get_iterable_array(stability, shift=1)
        current_manifold = BaseManifold(
            fixed_point.branch_points[orbit_indices[0]],
            stability,
            stretch_param=1,
            fixed_point=fixed_point,
            branch_index=0,
        )
        temp_root = current_manifold.root
        if isinstance(current_manifold.root, BranchPoint):
            current_manifold.root = current_manifold.walk_fwd(None, temp_root)
        if current_manifold.root is None:
            return
        current_manifold.stretch_param = current_manifold.root.stretch_param
        self._man_machine.iterate_manifold(
            current_manifold
        )  # side-effect: sets next_iterate
