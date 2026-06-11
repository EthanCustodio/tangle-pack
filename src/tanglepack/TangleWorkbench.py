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


class TangleWorkbench:

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

        self.compute_intersections(fixed_point)

        num_initial_intersections = len(self.Tangle._intersecting_segments)
        print(f"current intersections {num_initial_intersections}")

        for _ in range(max_iterations):

            self.grow_n_times(fixed_point, stability, num_iterations=1)

            self.compute_intersections(fixed_point)
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

    def compute_intersections(self, fixed_point):
        """This currently only computes homoclinic intersections, we will have
        to modify"""

        # Completely clear and rebuild the Tangle state to avoid stale references
        self.Tangle.clear_all()
        self._intersection_registry = IntersectionRegistry()

        self.index_manifolds(fixed_point, "unstable")
        self.index_manifolds(fixed_point, "stable")

        self.Tangle.populate_intersection_dict()

        for intersection in self.Tangle._intersections:
            self._intersection_registry.add(intersection)

        return list(self.Tangle._intersecting_coords.values())

    def plot_intersections(self, fp, ax=None, **scatter_kwargs):
        """
        Scatter-plot the last computed intersections for this session.
        If none computed yet for this fp, computes them first.
        """
        pts = np.array(list(self.Tangle._intersecting_coords.values()))
        if pts.size == 0:
            pts = self.compute_intersections(fp)
            if pts.size == 0:
                self.log.info("No intersections to plot.")
                return

        # sensible defaults; caller can override with kwargs
        scatter_kwargs.setdefault("s", 7)
        scatter_kwargs.setdefault("zorder", 10)
        scatter_kwargs.setdefault("color", "k")
        plt.scatter(pts[:, 0], pts[:, 1], **scatter_kwargs)

    def create_bridges(self, fixed_point: FixedPoint):

        bridges = self.Tangle.create_bridges()  # for_manifold=None → all intersections
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
            # no crossings — iterated bridge is already a valid unsplit bridge
            new_bridges = [iterated]

        # 5. wire genealogy
        bridge.iterated = True
        bridge.children = new_bridges
        for nb in new_bridges:
            nb.parent = bridge

        # 6. register
        self._bridges.extend(new_bridges)
        self._assign_bridge_intersections(new_bridges)

        return new_bridges

    def infer_iterate_table(self, coord_tol: float = 0.01) -> int:
        """
        Scan all iterated bridges and record the n=1 forward iterate relationship for
        each boundary intersection.

        Bridge topology identifies *which* intersections to process (only the two
        endpoints of each iterated bridge). The map is then applied directly to their
        coordinates and the closest registered intersection within ``coord_tol`` is
        recorded as the image. This avoids cdist ambiguity that arises when original
        and new intersections overlap in arc-length.

        Args:
            coord_tol: Maximum Euclidean distance (in phase-space units) between
                ``f(i_src)`` and a candidate intersection for the match to be accepted.
                Defaults to 0.01, which is ~10× the expected linear-intersection
                approximation error for well-resolved manifolds.

        Returns:
            Number of new iterate relationships recorded.
        """
        registry = self._intersection_registry
        f = self.dynamical_system.map
        all_items = [(ix_id, ix.get_point()) for ix_id, ix in registry]
        recorded = 0

        for bridge in self._bridges:
            if not bridge.iterated or not bridge.children:
                continue

            for src_id in (bridge.first_intersection, bridge.second_intersection):
                if src_id is None or (src_id, 1) in registry.iterate_table:
                    continue

                f_coords = f(registry[src_id].get_point())

                best_id, best_dist = None, float("inf")
                for tgt_id, tgt_coords in all_items:
                    d = float(np.linalg.norm(f_coords - tgt_coords))
                    if d < best_dist:
                        best_dist, best_id = d, tgt_id

                if best_id is not None and best_id != src_id and best_dist < coord_tol:
                    registry.register_iterate(src_id, 1, best_id)
                    recorded += 1

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
        G: nx.DiGraph,
        layout: str = "kamada_kawai",
        figsize: tuple = (12, 8),
        node_size_scale: float = 300,
        save_path: str = None,
        show_labels: bool = True,
        label_mode: str = "coords",
    ):
        """
        Visualize the intersection graph with edges colored by stability.

        Args:
            G: The intersection graph from build_intersection_graph()
            layout: Layout algorithm ("spring", "circular", "kamada_kawai", "spectral")
            figsize: Figure size as (width, height)
            node_size_scale: Scale factor for node sizes
            save_path: Optional path to save the figure
            show_labels: Whether to show node labels
            label_mode: Content of node labels. One of:
                ``"coords"`` — (x, y) coordinates (default)
                ``"cdist"``  — unstable and stable arc-length distances
                ``"id"``     — the intersection point ID
                ``"all"``    — coordinates, cdist, and ID together

        Returns:
            matplotlib figure and axis objects
        """
        import matplotlib.pyplot as plt

        if G.number_of_nodes() == 0:
            print("Graph has no nodes to visualize")
            return None, None

        fig, ax = plt.subplots(figsize=figsize)

        # Choose layout
        if layout == "spring":
            pos = nx.spring_layout(G, k=1, iterations=50)
        elif layout == "circular":
            pos = nx.circular_layout(G)
        elif layout == "kamada_kawai":
            pos = nx.kamada_kawai_layout(G)
        elif layout == "spectral":
            pos = nx.spectral_layout(G)
        else:
            pos = nx.spring_layout(G)

        mode_size_multiplier = {"coords": 1.0, "cdist": 1.0, "id": 1.0, "all": 2.5}
        effective_scale = node_size_scale * mode_size_multiplier.get(label_mode, 1.0)

        # Draw nodes (all same color - they're intersection points)
        nx.draw_networkx_nodes(
            G,
            pos,
            node_color="white",
            edgecolors="black",
            linewidths=2,
            node_size=[effective_scale * (1 + G.degree(n)) for n in G.nodes()],
            ax=ax,
        )

        def _arc3_mid_arrow_coords(
            p0_data: np.ndarray, p1_data: np.ndarray, rad: float
        ) -> tuple[np.ndarray, np.ndarray]:
            """Return (xy, xytext) in data coords for a midpoint arrow on an arc3 edge.

            Computed in display coordinates using matplotlib's exact Arc3 formula
            (ctrl = mid + rad*(dy, -dx)), then converted back to data coordinates.
            """
            trans = ax.transData
            inv = trans.inverted()
            p0 = trans.transform(p0_data)
            p1 = trans.transform(p1_data)
            mid = (p0 + p1) * 0.5
            dp = p1 - p0
            ctrl = np.array([mid[0] + rad * dp[1], mid[1] - rad * dp[0]])
            bm = 0.25 * p0 + 0.5 * ctrl + 0.25 * p1
            norm = np.linalg.norm(dp)
            direction = dp / norm if norm > 1e-6 else np.array([1.0, 0.0])
            eps = (
                direction * 4.0
            )  # 4 display-pixel offset gives the arrow its direction
            return inv.transform(bm + eps), inv.transform(bm - eps)

        # Separate edges by stability and draw them
        unstable_edges = [
            (u, v) for u, v, d in G.edges(data=True) if d["stability"] == "unstable"
        ]
        stable_edges = [
            (u, v) for u, v, d in G.edges(data=True) if d["stability"] == "stable"
        ]

        for edgelist, color, rad in [
            (unstable_edges, "#3b82f6", 0.2),
            (stable_edges, "#ef4444", -0.2),
        ]:
            nx.draw_networkx_edges(
                G,
                pos,
                edgelist=edgelist,
                edge_color=color,
                alpha=0.5,
                arrows=True,
                arrowsize=25,
                width=2.0,
                connectionstyle=f"arc3,rad={rad}",
                ax=ax,
            )
            for u, v in edgelist:
                p0, p1 = np.array(pos[u]), np.array(pos[v])
                xy, xytext = _arc3_mid_arrow_coords(p0, p1, rad)
                ax.annotate(
                    "",
                    xy=xy,
                    xytext=xytext,
                    arrowprops=dict(
                        arrowstyle="-|>",
                        color=color,
                        lw=1.5,
                        mutation_scale=18,
                    ),
                    alpha=0.85,
                )

        if show_labels:
            labels = {}
            for node in G.nodes():
                data = G.nodes[node]
                coords = data.get("coords")
                u_cd = data.get("unstable_cdist")
                s_cd = data.get("stable_cdist")

                parts: list[str] = []
                if label_mode in ("coords", "all"):
                    if coords is not None:
                        parts.append(f"({coords[0]:.2f}, {coords[1]:.2f})")
                    else:
                        parts.append("(?)")
                if label_mode in ("cdist", "all"):
                    u_str = f"{u_cd:.3f}" if u_cd is not None else "?"
                    s_str = f"{s_cd:.3f}" if s_cd is not None else "?"
                    parts.append(f"u:{u_str} s:{s_str}")
                if label_mode in ("id", "all"):
                    parts.append(f"id:{node}")

                labels[node] = "\n".join(parts) if parts else str(node)

            nx.draw_networkx_labels(
                G, pos, labels, font_size=7, font_weight="bold", ax=ax
            )

        ax.set_title(
            "Intersection Graph\n"
            + "(white nodes = intersections, blue edges = unstable flow, red edges = stable flow)",
            fontsize=12,
            fontweight="bold",
        )

        from matplotlib.lines import Line2D

        legend_handles = [
            Line2D([0], [0], color="#3b82f6", linewidth=2.5, label="Unstable"),
            Line2D([0], [0], color="#ef4444", linewidth=2.5, label="Stable"),
        ]
        ax.legend(handles=legend_handles, loc="upper right", fontsize=10)
        ax.axis("off")

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches="tight")
            print(f"Graph visualization saved to {save_path}")

        plt.show()
        return fig, ax

    # def plot_all_bridges(self, bridges):

    #     n = len(bridges)
    #     colors = cm.get_cmap("tab20", n)  # or 'tab10', 'nipy_spectral', etc.

    #     for i, bridge in enumerate(bridges):

    #         vibe = colors(i)
    #         bridge.plot(color=vibe)

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
        colors = cm.get_cmap("tab20", n)
        for i, bridge in enumerate(bridges):
            bridge.plot(color=colors(i))

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
