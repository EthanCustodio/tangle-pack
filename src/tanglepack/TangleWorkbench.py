from typing import Callable, Literal, Iterable
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
from .Tangle import Tangle
from .FixedPoint import FixedPoint
from .BaseManifold import BaseManifold

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

    def initialize_manifold(self, fixed_point: FixedPoint, stability: Stability):

        initial_segments = self._man_maker.construct_kevin_way(fixed_point, stability)

        for (orbit_index, branch_index), manifold in initial_segments.items():

            self.manifolds[(fixed_point, stability, orbit_index, branch_index)] = (
                manifold
            )

        return initial_segments

    def initialize_both_manifolds(self, fixed_point: FixedPoint):

        unstable_segments = self.initialize_manifold(fixed_point, "unstable")
        stable_segments = self.initialize_manifold(fixed_point, "stable")

        for (orbit_index, branch_index), manifold in unstable_segments.items():

            self.manifolds[(fixed_point, "unstable", orbit_index, branch_index)] = (
                manifold
            )

        for (orbit_index, branch_index), manifold in stable_segments.items():

            self.manifolds[(fixed_point, "stable", orbit_index, branch_index)] = (
                manifold
            )

        return (unstable_segments, stable_segments)

    def grow_n_times(
        self,
        fixed_point: FixedPoint,
        stability: Stability,
        num_iterations: int,
    ) -> None:

        if self.manifolds.get((fixed_point, stability, 0, 0)) is None:
            raise ValueError(f"""Manifold for fixed point {fixed_point} 
                    with stability {stability} has not been initialized.
                    Please run initialize_manifold first.""")

        self._man_machine.grow_x_times(fixed_point, stability, num_iterations)

        for (fp, stab, _orbit_index, _branch_index), manifold in self.manifolds.items():

            if fp is fixed_point and stab == stability:

                manifold._find_tail()

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
    ) -> None:
        """
        Grows the manifold until a turnaround is detected or max_iterations is reached.

        Args:
            fixed_point (FixedPoint): The fixed point whose manifold is to be grown.
            stability (Stability): The stability type of the manifold ('stable' or 'unstable').
            max_iterations (int, optional): Maximum number of iterations to grow. Defaults to 50.
        """

        if self.manifolds.get((fixed_point, stability, 0, 0)) is None:
            raise ValueError(f"""Manifold for fixed point {fixed_point} 
                    with stability {stability} has not been initialized.
                    Please run initialize_manifold first.""")

        root = fixed_point.branch_points[0]
        first_point = self.manifolds.get((fixed_point, stability, 0, 0))
        first_point = first_point.walk_fwd(None, root, 0)

        root_coord = root._coords
        first_point_coords = first_point._coords

        initial_direction = np.asarray(first_point_coords) - np.asarray(root_coord)

        for _ in range(max_iterations):

            self.grow_n_times(fixed_point, stability, num_iterations=1)

            tail = self.manifolds.get((fixed_point, stability, 0, 0)).tail
            first_point = self.manifolds.get((fixed_point, stability, 0, 0))
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
        self, fixed_point: FixedPoint, stability: Stability, length: float
    ):

        if self.manifolds.get((fixed_point, stability, 0, 0)) is None:
            raise ValueError(f"""Manifold for fixed point {fixed_point} 
                    with stability {stability} has not been initialized.
                    Please run initialize_manifold first.""")

        # TODO change this so it uses the actual arclength
        current_distance = self.manifolds.get((fixed_point, stability, 0, 0)).tail.cdist

        while current_distance < length:

            self.grow_n_times(fixed_point, stability, num_iterations=1)

            current_distance = self.manifolds.get(
                (fixed_point, stability, 0, 0)
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

        self.index_manifolds(fixed_point, "unstable")
        self.index_manifolds(fixed_point, "stable")

        self.Tangle.populate_intersection_dict()

        points = self.Tangle._intersecting_coords.values()

        return points

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

        bridges = self.Tangle.create_bridges()

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

    def visualize_intersection_graph(
        self,
        G: nx.DiGraph,
        layout: str = "spring",
        figsize: tuple = (12, 8),
        node_size_scale: float = 300,
        save_path: str = None,
        show_labels: bool = True,
    ):
        """
        Visualize the intersection graph with edges colored by stability.

        Args:
            G: The intersection graph from build_intersection_graph()
            layout: Layout algorithm ("spring", "circular", "kamada_kawai", "spectral")
            figsize: Figure size as (width, height)
            node_size_scale: Scale factor for node sizes
            save_path: Optional path to save the figure
            show_labels: Whether to show node labels (coordinates)

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

        # Draw nodes (all same color - they're intersection points)
        nx.draw_networkx_nodes(
            G,
            pos,
            node_color="white",
            edgecolors="black",
            linewidths=2,
            node_size=[node_size_scale * (1 + G.degree(n)) for n in G.nodes()],
            ax=ax,
        )

        # Separate edges by stability and draw them
        unstable_edges = [
            (u, v) for u, v, d in G.edges(data=True) if d["stability"] == "unstable"
        ]
        stable_edges = [
            (u, v) for u, v, d in G.edges(data=True) if d["stability"] == "stable"
        ]

        nx.draw_networkx_edges(
            G,
            pos,
            edgelist=unstable_edges,
            edge_color="#3b82f6",  # blue
            alpha=0.6,
            arrows=True,
            arrowsize=20,
            width=2.5,
            label="Unstable",
            ax=ax,
        )

        nx.draw_networkx_edges(
            G,
            pos,
            edgelist=stable_edges,
            edge_color="#ef4444",  # red
            alpha=0.6,
            arrows=True,
            arrowsize=20,
            width=2.5,
            label="Stable",
            ax=ax,
        )

        # Create labels showing coordinates
        if show_labels:
            labels = {}
            for node in G.nodes():
                coords = G.nodes[node].get("coords")
                if coords is not None:
                    labels[node] = f"({coords[0]:.2f}, {coords[1]:.2f})"
                else:
                    labels[node] = str(node)

            nx.draw_networkx_labels(
                G, pos, labels, font_size=7, font_weight="bold", ax=ax
            )

        ax.set_title(
            "Intersection Graph\n"
            + "(white nodes = intersections, blue edges = unstable flow, red edges = stable flow)",
            fontsize=12,
            fontweight="bold",
        )
        ax.legend(loc="upper right", fontsize=10)
        ax.axis("off")

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches="tight")
            print(f"Graph visualization saved to {save_path}")

        return fig, ax

    def plot_all_bridges(self, bridges):

        n = len(bridges)
        colors = cm.get_cmap("tab20", n)  # or 'tab10', 'nipy_spectral', etc.

        for i, bridge in enumerate(bridges):

            vibe = colors(i)
            bridge.plot(color=vibe)

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
