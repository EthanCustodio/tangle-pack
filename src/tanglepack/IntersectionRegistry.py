from __future__ import annotations

import bisect
from typing import Callable, Literal, Optional, TYPE_CHECKING

import networkx as nx
import numpy as np
from numpy.typing import NDArray

from .Intersection import Intersection, ManifoldKey
from .IterateTable import IterateTable

if TYPE_CHECKING:
    from .FixedPoint import FixedPoint


class IntersectionRegistry:
    """
    Master store of all intersection points in the tangle.

    IDs are contiguous integers (0, 1, 2, …) assigned in insertion order.
    Duplicate intersections (same cdists within tolerance) are detected and
    deduped — the existing ID is returned rather than creating a second entry.

    Primary interface:
        registry.add(intersection)               → int (assigned ID)
        registry[id]                             → Intersection
        registry.iterate_table[id, n]            → int or None
        registry.by_unstable_cdist               → list[int]  (sorted by u-cdist)
        registry.by_stable_cdist                 → list[int]  (sorted by s-cdist)
        registry.graph                           → nx.MultiDiGraph (live)

    Query interface (all return list[Intersection]):
        registry.on_interval(lo, hi)             → pre-images that map into [lo, hi]
        registry.on_cdist_range(lo, hi)          → intersections with cdist in [lo, hi]
        registry.from_fixed_point(fp)            → intersections involving fp
        registry.from_branch(branch_index)       → intersections on given branch
        registry.filter(predicate)               → arbitrary predicate

    Attributes:
        _store: dict[int, Intersection]
        _next_id: int
        cdist_tol: float
        iterate_table: IterateTable
        _unstable_order: list[int]  — IDs sorted ascending by unstable_cdist
        _stable_order: list[int]    — IDs sorted ascending by stable_cdist
        _cdist_index: dict[tuple[float, float], int]  — secondary collision index
        _graph: nx.MultiDiGraph
        _graph_adjacency_dirty: bool
    """

    def __init__(self, cdist_tol: float = 1e-6):
        self._store: dict[int, Intersection] = {}
        self._next_id: int = 0
        self.cdist_tol = cdist_tol
        self.iterate_table = IterateTable()

        self._unstable_order: list[int] = []
        self._stable_order: list[int] = []
        self._cdist_index: dict[tuple[float, float], int] = {}

        self._graph: nx.MultiDiGraph = nx.MultiDiGraph()
        self._graph_adjacency_dirty: bool = False

    # ── core insert / lookup ───────────────────────────────────────────────

    def add(self, intersection: Intersection) -> int:
        """
        Register an intersection, returning its unique ID.

        If a collision is detected (another intersection with cdists within
        self.cdist_tol), the existing ID is returned and no duplicate is stored.
        On a new insertion, the node is added to self.graph immediately.

        Args:
            intersection: The Intersection to register.

        Returns:
            The integer ID (new or existing on collision).
        """
        existing = self._find_collision(intersection)
        if existing is not None:
            return existing

        new_id = self._next_id
        self._next_id += 1
        intersection.id = new_id
        self._store[new_id] = intersection

        self._insert_into_unstable_order(new_id, intersection.unstable_cdist)
        self._insert_into_stable_order(new_id, intersection.stable_cdist)

        key = self._cdist_key(intersection)
        self._cdist_index[key] = new_id

        # Add node to the live graph.
        self._graph.add_node(
            new_id,
            coords=intersection.coords,
            unstable_cdist=intersection.unstable_cdist,
            stable_cdist=intersection.stable_cdist,
            manifold_a_key=intersection.manifold_a_key,
            manifold_b_key=intersection.manifold_b_key,
            label=intersection.label,
        )
        self._graph_adjacency_dirty = True

        return new_id

    def add_synthetic(
        self,
        coords: tuple[float, float],
        unstable_cdist: float,
        stable_cdist: float,
        label: Optional[str] = None,
        manifold_a_key: Optional[ManifoldKey] = None,
        manifold_b_key: Optional[ManifoldKey] = None,
    ) -> int:
        """
        Low-level utility: add an intersection not backed by a detected segment
        crossing (e.g., a manually placed reference point).

        For most workflows, prefer computing real intersections and querying them
        with on_interval() and related methods.
        """
        return self.add(
            Intersection(
                coords=coords,
                unstable_cdist=unstable_cdist,
                stable_cdist=stable_cdist,
                label=label,
                manifold_a_key=manifold_a_key,
                manifold_b_key=manifold_b_key,
            )
        )

    def __getitem__(self, id: int) -> Intersection:
        """registry[id] → Intersection. Raises KeyError if id is unknown."""
        return self._store[id]

    def __len__(self) -> int:
        return len(self._store)

    def __contains__(self, id: int) -> bool:
        return id in self._store

    def __iter__(self):
        """Iterate over all (id, Intersection) pairs in insertion order."""
        return iter(self._store.items())

    # ── ordering views ─────────────────────────────────────────────────────

    @property
    def by_unstable_cdist(self) -> list[int]:
        """IDs sorted ascending by unstable_cdist."""
        return list(self._unstable_order)

    @property
    def by_stable_cdist(self) -> list[int]:
        """IDs sorted ascending by stable_cdist."""
        return list(self._stable_order)

    def unstable_rank(self, id: int) -> int:
        """0-based position of intersection `id` in the W^u ordering."""
        return self._unstable_order.index(id)

    def stable_rank(self, id: int) -> int:
        """0-based position of intersection `id` in the W^s ordering."""
        return self._stable_order.index(id)

    def all_ids(self) -> list[int]:
        """All registered IDs in insertion order."""
        return list(self._store.keys())

    # ── live graph ─────────────────────────────────────────────────────────

    @property
    def graph(self) -> nx.MultiDiGraph:
        """
        The live intersection graph.

        Nodes are intersection IDs. Two edge types are maintained:
          - type="adjacency", stability="unstable"/"stable":
              connects consecutive intersections in the respective sorted ordering.
              These are rebuilt lazily whenever the sorted order changes.
          - type="iterate", n=<int>:
              directed edge (p → f^n(p)) added by register_iterate().

        The graph is always up to date w.r.t. nodes and iterate edges.
        Adjacency edges are rebuilt on first access after any add().
        """
        if self._graph_adjacency_dirty:
            self._rebuild_adjacency_edges()
        return self._graph

    def _rebuild_adjacency_edges(self):
        """Remove and rebuild all adjacency-type edges from the current sorted lists."""
        stale = [
            (u, v, k)
            for u, v, k, d in self._graph.edges(keys=True, data=True)
            if d.get("type") == "adjacency"
        ]
        for u, v, k in stale:
            self._graph.remove_edge(u, v, k)

        for i in range(len(self._unstable_order) - 1):
            u, v = self._unstable_order[i], self._unstable_order[i + 1]
            self._graph.add_edge(
                u,
                v,
                key=f"adj_u_{i}",
                type="adjacency",
                stability="unstable",
            )

        for i in range(len(self._stable_order) - 1):
            u, v = self._stable_order[i], self._stable_order[i + 1]
            self._graph.add_edge(
                u,
                v,
                key=f"adj_s_{i}",
                type="adjacency",
                stability="stable",
            )

        self._graph_adjacency_dirty = False

    # ── iterate table ──────────────────────────────────────────────────────

    def register_iterate(self, source_id: int, n: int, target_id: int):
        """
        Record f^n(source) = target.

        Delegates to iterate_table and wires the directed edge into the graph.
        """
        self.iterate_table.register_iterate(source_id, n, target_id)
        if source_id in self._graph and target_id in self._graph:
            self._graph.add_edge(
                source_id,
                target_id,
                key=f"iter_{source_id}_{n}",
                type="iterate",
                n=n,
            )

    def infer_iterates(
        self,
        max_depth: int = 10,
        tol_multiplier: float = 10.0,
    ) -> int:
        """
        Scan all registered intersections and infer iterate relationships by
        predicting cdists and checking for matches.

        Reads lambda_u per intersection from its manifold_a_key, so this works
        correctly with multiple fixed points in the same registry.

        Unlike v1, this method takes no lambda_u parameter.

        Args:
            max_depth: How many iterate levels to search forward and backward.
            tol_multiplier: Scale factor on cdist_tol for matching.

        Returns:
            Number of new iterate relationships recorded.
        """
        tol = self.cdist_tol * tol_multiplier
        recorded = 0

        for source_id, source in self._store.items():
            lambda_u = self._get_lambda_u(source)
            if lambda_u is None:
                continue

            for n in range(1, max_depth + 1):
                if (source_id, n) not in self.iterate_table:
                    pred_u = (lambda_u**n) * source.unstable_cdist
                    pred_s = source.stable_cdist / (lambda_u**n)
                    target_id = self.find_by_cdist(pred_u, pred_s, tol)
                    if target_id is not None:
                        self.register_iterate(source_id, n, target_id)
                        recorded += 1

                if (source_id, -n) not in self.iterate_table:
                    pred_u = source.unstable_cdist / (lambda_u**n)
                    pred_s = source.stable_cdist * (lambda_u**n)
                    target_id = self.find_by_cdist(pred_u, pred_s, tol)
                    if target_id is not None:
                        self.register_iterate(source_id, -n, target_id)
                        recorded += 1

        return recorded

    # ── query interface ────────────────────────────────────────────────────

    def filter(self, predicate: Callable[[Intersection], bool]) -> list[Intersection]:
        """
        Return all registered intersections that satisfy predicate.

        This is the extensibility core. All named query methods below are thin
        wrappers that build a predicate and delegate here.

        Args:
            predicate: A callable that takes an Intersection and returns bool.

        Returns:
            List of matching Intersection objects (no guaranteed ordering).

        Example:
            registry.filter(lambda ix: ix.unstable_cdist > 10.0)
        """
        return [ix for ix in self._store.values() if predicate(ix)]

    def on_interval(
        self,
        lo: float,
        hi: float,
        stability: Literal["unstable", "stable"] = "unstable",
        fixed_point: Optional[FixedPoint] = None,
        branch_index: Optional[int] = None,
    ) -> list[Intersection]:
        """
        Return all intersections p such that f(p) has cdist on `stability` in [lo, hi].

        Uses the eigenvalue formula — no additional manifold numerics required:
            f(p).unstable_cdist = lambda_u × p.unstable_cdist
            f(p).stable_cdist   = p.stable_cdist / lambda_u

        lambda_u is read per intersection from manifold_a_key (for stability="unstable")
        or manifold_b_key (for stability="stable"), so a registry containing intersections
        from multiple fixed points with different eigenvalues is handled correctly.

        Args:
            lo: Lower cdist bound for f(p) on the given manifold.
            hi: Upper cdist bound for f(p) on the given manifold.
            stability: Which manifold's cdist to apply the interval on.
            fixed_point: If given, only consider intersections involving this FP.
            branch_index: If given, only consider intersections on this branch.

        Returns:
            List of Intersection objects (the sources, not their images).
        """
        results = []
        for ix in self._store.values():
            # Optional provenance filters
            if fixed_point is not None or branch_index is not None:
                key = (
                    ix.manifold_a_key if stability == "unstable" else ix.manifold_b_key
                )
                if key is None:
                    continue
                if fixed_point is not None and key[0] is not fixed_point:
                    continue
                if branch_index is not None and key[3] != branch_index:
                    continue

            lambda_u = self._get_lambda_u(ix)
            if lambda_u is None:
                continue

            if stability == "unstable":
                image_cdist = lambda_u * ix.unstable_cdist
            else:
                image_cdist = ix.stable_cdist / lambda_u

            if lo <= image_cdist <= hi:
                results.append(ix)

        return results

    def on_cdist_range(
        self,
        lo: float,
        hi: float,
        stability: Literal["unstable", "stable"] = "unstable",
    ) -> list[Intersection]:
        """
        Return all intersections whose CURRENT cdist on `stability` is in [lo, hi].

        Unlike on_interval(), no iteration is applied — this filters by where the
        intersection already sits, not where it maps to.

        Args:
            lo: Lower bound of the cdist range.
            hi: Upper bound of the cdist range.
            stability: Which manifold cdist to filter on.

        Returns:
            List of matching Intersection objects sorted by the chosen cdist.
        """
        attr = "unstable_cdist" if stability == "unstable" else "stable_cdist"
        return self.filter(lambda ix: lo <= getattr(ix, attr) <= hi)

    def from_fixed_point(
        self,
        fp: FixedPoint,
        stability: Optional[Literal["unstable", "stable"]] = None,
    ) -> list[Intersection]:
        """
        Return all intersections that involve the given fixed point.

        Checks both manifold_a_key and manifold_b_key. Optionally restrict to
        intersections where fp appears on the specified stability side.

        Args:
            fp: The FixedPoint to filter by.
            stability: If given, only match on that side (manifold_a for "unstable",
                manifold_b for "stable").

        Returns:
            List of matching Intersection objects.
        """

        def pred(ix: Intersection) -> bool:
            a_match = ix.manifold_a_key is not None and ix.manifold_a_key[0] is fp
            b_match = ix.manifold_b_key is not None and ix.manifold_b_key[0] is fp
            if stability == "unstable":
                return a_match
            if stability == "stable":
                return b_match
            return a_match or b_match

        return self.filter(pred)

    def from_branch(
        self,
        branch_index: int,
        stability: Optional[Literal["unstable", "stable"]] = None,
    ) -> list[Intersection]:
        """
        Return all intersections on the given branch.

        Args:
            branch_index: 0 or 1.
            stability: If given, only check the corresponding manifold key side.

        Returns:
            List of matching Intersection objects.
        """

        def pred(ix: Intersection) -> bool:
            a_match = (
                ix.manifold_a_key is not None and ix.manifold_a_key[3] == branch_index
            )
            b_match = (
                ix.manifold_b_key is not None and ix.manifold_b_key[3] == branch_index
            )
            if stability == "unstable":
                return a_match
            if stability == "stable":
                return b_match
            return a_match or b_match

        return self.filter(pred)

    # ── cdist-based lookup ─────────────────────────────────────────────────

    def find_by_cdist(
        self,
        unstable_cdist: float,
        stable_cdist: float,
        tol: Optional[float] = None,
    ) -> Optional[int]:
        """
        Find the ID of an intersection whose cdists are within tol of the given values.

        Args:
            unstable_cdist: Target unstable cdist.
            stable_cdist: Target stable cdist.
            tol: Search tolerance. Defaults to self.cdist_tol.

        Returns:
            Matching ID, or None.
        """
        if tol is None:
            tol = self.cdist_tol
        for id, existing in self._store.items():
            if (
                abs(existing.unstable_cdist - unstable_cdist) < tol
                and abs(existing.stable_cdist - stable_cdist) < tol
            ):
                return id
        return None

    # ── array exports ──────────────────────────────────────────────────────

    def as_forward_array(self, max_depth: int = 5) -> NDArray[np.int64]:
        """
        Dense array A[i, d-1] = ID of f^d(ids[i]) where ids = all_ids() in order.
        Shape: (N, max_depth). -1 = unknown.
        """
        return self.iterate_table.as_forward_array(self.all_ids(), max_depth)

    def as_backward_array(self, max_depth: int = 5) -> NDArray[np.int64]:
        """Dense array B[i, d-1] = ID of f^{-d}(ids[i]). Shape (N, max_depth). -1 = unknown."""
        return self.iterate_table.as_backward_array(self.all_ids(), max_depth)

    def unstable_order_array(self) -> NDArray[np.float64]:
        """(N, 3) array [id, unstable_cdist, stable_cdist] sorted by unstable_cdist."""
        rows = [
            [iid, self._store[iid].unstable_cdist, self._store[iid].stable_cdist]
            for iid in self._unstable_order
        ]
        return np.array(rows, dtype=np.float64)

    def stable_order_array(self) -> NDArray[np.float64]:
        """(N, 3) array [id, unstable_cdist, stable_cdist] sorted by stable_cdist."""
        rows = [
            [iid, self._store[iid].unstable_cdist, self._store[iid].stable_cdist]
            for iid in self._stable_order
        ]
        return np.array(rows, dtype=np.float64)

    # ── internal helpers ───────────────────────────────────────────────────

    def _get_lambda_u(self, intersection: Intersection) -> Optional[float]:
        """
        Read the unstable eigenvalue magnitude from the intersection's manifold keys.

        For stability="unstable", uses manifold_a_key (the unstable manifold side).
        For stability="stable", uses manifold_b_key.
        If neither key is set, returns None and the intersection is skipped by callers.
        """
        key = intersection.manifold_a_key or intersection.manifold_b_key
        if key is None:
            return None
        fp = key[0]
        if not hasattr(fp, "unstable_eigenvalues") or not fp.unstable_eigenvalues:
            return None
        return abs(fp.unstable_eigenvalues[0])

    def _find_collision(self, intersection: Intersection) -> Optional[int]:
        """Linear scan for an existing intersection within cdist_tol."""
        for id, existing in self._store.items():
            if (
                abs(existing.unstable_cdist - intersection.unstable_cdist)
                < self.cdist_tol
                and abs(existing.stable_cdist - intersection.stable_cdist)
                < self.cdist_tol
            ):
                return id
        return None

    def _cdist_key(self, intersection: Intersection) -> tuple[float, float]:
        digits = max(0, -int(np.floor(np.log10(self.cdist_tol))) - 1)
        return (
            round(intersection.unstable_cdist, digits),
            round(intersection.stable_cdist, digits),
        )

    def _insert_into_unstable_order(self, id: int, unstable_cdist: float):
        keys = [self._store[i].unstable_cdist for i in self._unstable_order]
        pos = bisect.bisect_left(keys, unstable_cdist)
        self._unstable_order.insert(pos, id)

    def _insert_into_stable_order(self, id: int, stable_cdist: float):
        keys = [self._store[i].stable_cdist for i in self._stable_order]
        pos = bisect.bisect_left(keys, stable_cdist)
        self._stable_order.insert(pos, id)
