from __future__ import annotations

import bisect
import logging
from typing import Callable, Iterable, Literal, Optional, TYPE_CHECKING

import networkx as nx
import numpy as np
from numpy.typing import NDArray

from .Intersection import Intersection, ManifoldKey
from .IterateTable import IterateTable

if TYPE_CHECKING:
    from .FixedPoint import FixedPoint

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class IntersectionRegistry:
    """
    Master store of all intersection points in the tangle.

    IDs are contiguous integers (0, 1, 2, …) assigned in insertion order.
    An intersection is uniquely identified by (unstable_cdist, stable_cdist,
    unstable branch, stable branch) — the two cdists plus the two manifold keys it
    crosses on. Duplicates (a match on all four, cdists within tolerance) are
    detected and deduped — the existing ID is returned rather than creating a
    second entry. The branch keys are required because canonical distance is
    measured per branch: the periodic anchor of every branch sits at cdist (0, 0),
    so cdist alone would wrongly merge the distinct anchors of a period-p orbit
    into one. Branch keys are used rather than phase-space coordinates because the
    coordinates are only approximate (drift under the nonlinear map), while the
    cdists and branch identity are exact/stable.

    Primary interface:
        registry.add(intersection)               → int (assigned ID)
        registry[id]                             → Intersection
        registry.iterate_table[id, n]            → int or None
        registry.by_unstable_cdist               → list[int]  (sorted by u-cdist)
        registry.by_stable_cdist                 → list[int]  (sorted by s-cdist)
        registry.graph(bridges=None)             → nx.MultiDiGraph (live)

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
        _generation: int  — bumped by every mutation; read via :attr:`generation`
        _orders_version: int  — bumped by the mutations that change the stored
            crossings (insert, renumber); the ordering-derived views key on it
        _unstable_order: list[int]  — IDs sorted ascending by unstable_cdist
        _unstable_keys: list[float] — those IDs' unstable cdists, same order
        _stable_order: list[int]    — IDs sorted ascending by stable_cdist
        _stable_keys: list[float]   — those IDs' stable cdists, same order
        _graph: nx.MultiDiGraph
        _graph_adjacency_dirty: bool
        _graph_bridges: Optional[list[tuple[int, int]]]

    Note:
        The two sorted key arrays are maintained on insertion (bisect), so no
        query ever rebuilds them by scanning the store. Every derived view
        (rank maps, the per-stable-branch candidate index) is memoised against
        ``_orders_version`` and rebuilt lazily on the next read after an insert
        or a renumbering.
    """

    def __init__(self, cdist_tol: float = 1e-6):
        self._store: dict[int, Intersection] = {}
        self._next_id: int = 0
        self.cdist_tol = cdist_tol
        self.iterate_table = IterateTable()

        # Monotone mutation counter. Anything derived from this registry (a
        # Trellis snapshot, the workbench's generation) is valid only for the
        # value it was built at.
        self._generation: int = 0
        # Bumped by the subset of mutations that change the STORED CROSSINGS --
        # an insert or a renumbering, not an iterate relation. The views below
        # are functions of those alone, so keying them on this rather than on
        # the full generation keeps them alive across the thousands of
        # register_iterate calls one infer_iterates pass makes.
        self._orders_version: int = 0

        self._unstable_order: list[int] = []
        self._unstable_keys: list[float] = []
        self._stable_order: list[int] = []
        self._stable_keys: list[float] = []

        # Lazily rebuilt views, each stamped with the _orders_version it was
        # built at.
        self._rank_maps: Optional[tuple[dict[int, int], dict[int, int]]] = None
        self._rank_version: int = -1
        self._stable_branch_index: Optional[
            dict[Optional[ManifoldKey], tuple[list[float], list[int]]]
        ] = None
        self._stable_branch_version: int = -1

        self._graph: nx.MultiDiGraph = nx.MultiDiGraph()
        self._graph_adjacency_dirty: bool = False
        # The bridge endpoint pairs the current unstable adjacency was built
        # from (None = the per-branch cdist fallback).
        self._graph_bridges: Optional[list[tuple[int, int]]] = None

    # ── generation ─────────────────────────────────────────────────────────

    @property
    def generation(self) -> int:
        """
        Monotone counter bumped by every mutation of the registry's contents.

        Mutations are: a new intersection stored (:meth:`add`, and therefore
        :meth:`add_synthetic`), a renumbering (:meth:`reindex_from`), and a new
        iterate relation (:meth:`register_iterate`). An :meth:`add` that
        collides with a stored crossing stores nothing and is NOT a mutation --
        it is a lookup that happens to be spelled like an insert.

        Rebuilding a derived cache (the graph's adjacency edges, the rank maps)
        is not a mutation either: those are functions of the contents, so they
        never invalidate anything.

        Returns:
            int: The current generation.
        """
        return self._generation

    # ── core insert / lookup ───────────────────────────────────────────────

    def add(self, intersection: Intersection) -> int:
        """
        Register an intersection, returning its unique ID.

        If a collision is detected (another intersection with cdists within
        self.cdist_tol), the existing ID is returned and no duplicate is stored.
        On a new insertion, the node is added to the live graph immediately.

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
        self._generation += 1
        self._orders_version += 1

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
            Intersection.synthetic(
                coords=coords,
                unstable_cdist=unstable_cdist,
                stable_cdist=stable_cdist,
                label=label,
                manifold_a_key=manifold_a_key,
                manifold_b_key=manifold_b_key,
            )
        )

    def reindex_from(self, old: "IntersectionRegistry") -> dict[int, int]:
        """
        Re-number this registry so crossings carried over from ``old`` keep their id.

        IDs are normally insertion-order integers, so any rebuild of the registry
        (e.g. the recompute after trimming a stable manifold for a resonance zone)
        renumbers every crossing from zero. That invalidates ids the caller is
        holding — a strong pip chosen as id 7 becomes some unrelated integer. This
        re-aligns the freshly built registry against ``old``: every intersection
        that matches one in ``old`` (same canonical distances within ``cdist_tol``
        and the same two branch keys — the exact/stable identity used by
        :meth:`_find_collision`) is given back its old id. Crossings with no match
        in ``old`` (genuinely new geometry) are assigned fresh ids above every
        preserved one, so no id is ever reused for a different point.

        Call this after all intersections have been added but before the iterate
        table is filled (the table is keyed by id). The sorted orderings, cdist
        index, and graph nodes are rebuilt with the final ids.

        Args:
            old: The registry that existed before the rebuild, whose ids should be
                preserved wherever the same crossing reappears.

        Returns:
            Mapping ``{final_id: matched_old_id_or_None}`` for every intersection
            now stored — useful for logging/debugging the remap.
        """
        intersections = list(self._store.values())

        used: set[int] = set()
        final_id: dict[int, int] = {}  # id(intersection object) -> assigned id
        for ix in intersections:
            match = old._find_collision(ix)
            if match is not None and match not in used:
                final_id[id(ix)] = match
                used.add(match)

        counter = 0
        for ix in intersections:
            if id(ix) not in final_id:
                while counter in used:
                    counter += 1
                final_id[id(ix)] = counter
                used.add(counter)

        # Rebuild every id-keyed structure from scratch with the final ids.
        self._store = {}
        self._unstable_order = []
        self._unstable_keys = []
        self._stable_order = []
        self._stable_keys = []
        self._graph = nx.MultiDiGraph()
        self.iterate_table = IterateTable()
        self._next_id = (max(used) + 1) if used else 0
        self._generation += 1
        self._orders_version += 1

        remap: dict[int, int] = {}
        for ix in intersections:
            fid = final_id[id(ix)]
            ix.id = fid
            self._store[fid] = ix
            self._insert_into_unstable_order(fid, ix.unstable_cdist)
            self._insert_into_stable_order(fid, ix.stable_cdist)
            self._graph.add_node(
                fid,
                coords=ix.coords,
                unstable_cdist=ix.unstable_cdist,
                stable_cdist=ix.stable_cdist,
                manifold_a_key=ix.manifold_a_key,
                manifold_b_key=ix.manifold_b_key,
                label=ix.label,
            )
            remap[fid] = old._find_collision(ix)
        self._graph_adjacency_dirty = True
        self._graph_bridges = None

        return remap

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
        """0-based position of intersection `id` in the W^u ordering.

        Raises:
            KeyError: If ``id`` is not registered.
        """
        return self._ranks()[0][id]

    def stable_rank(self, id: int) -> int:
        """0-based position of intersection `id` in the W^s ordering.

        Raises:
            KeyError: If ``id`` is not registered.
        """
        return self._ranks()[1][id]

    def _ranks(self) -> tuple[dict[int, int], dict[int, int]]:
        """The (unstable, stable) id-to-rank maps, rebuilt when the orders change."""
        if self._rank_maps is None or self._rank_version != self._orders_version:
            self._rank_maps = (
                {iid: rank for rank, iid in enumerate(self._unstable_order)},
                {iid: rank for rank, iid in enumerate(self._stable_order)},
            )
            self._rank_version = self._orders_version
        return self._rank_maps

    def all_ids(self) -> list[int]:
        """All registered IDs in insertion order."""
        return list(self._store.keys())

    # ── live graph ─────────────────────────────────────────────────────────

    def graph(
        self, bridges: Optional[Iterable[tuple[int, int]]] = None
    ) -> nx.MultiDiGraph:
        """
        The live intersection graph.

        Nodes are intersection IDs. Two edge types are maintained:
          - type="adjacency", stability="unstable"/"stable":
              connects crossings that are NEIGHBOURS ON ONE BRANCH -- see
              :meth:`_rebuild_adjacency_edges`. Rebuilt lazily.
          - type="iterate", n=<int>:
              directed edge (p → f^n(p)) added by register_iterate().

        The graph is always up to date w.r.t. nodes and iterate edges.

        Args:
            bridges: Optional ``(first_id, second_id)`` endpoint pairs of the
                bridges cut out of the unstable manifold. When given, these ARE
                the unstable adjacency: a bridge is by definition the piece of
                unstable manifold between two consecutive crossings, and the
                registry has no other way to know where the manifold was cut
                (an iterated bridge lays new curve the per-branch cdist order
                cannot see). Without them the unstable side falls back to the
                per-branch canonical-distance order.

        Returns:
            nx.MultiDiGraph: The live graph. Callers that decorate it must work
            on a copy.

        Note:
            The returned object is owned by the registry and is rebuilt in
            place; a different ``bridges`` argument rebuilds the adjacency. In
            particular a later call WITHOUT ``bridges`` reverts the live unstable
            adjacency to the per-branch cdist fallback, so a caller that wants
            the bridge adjacency must pass it every time (or hold a copy).
        """
        bridge_pairs = None if bridges is None else [tuple(b) for b in bridges]
        if self._graph_adjacency_dirty or bridge_pairs != self._graph_bridges:
            self._rebuild_adjacency_edges(bridge_pairs)
        return self._graph

    def _rebuild_adjacency_edges(
        self, bridges: Optional[list[tuple[int, int]]] = None
    ) -> None:
        """
        Rebuild every adjacency edge from the PER-BRANCH orderings.

        Canonical distance is measured from each branch's own anchor, so a global
        cdist sort is meaningless across branches: on a period-p orbit it
        interleaves the p stable branches and joins crossings that are not
        neighbours on any curve. Adjacency is therefore built one branch at a
        time -- stable edges from the crossings sharing a ``manifold_b_key``, in
        stable-cdist order; unstable edges from the given bridges, or from the
        crossings sharing a ``manifold_a_key`` in unstable-cdist order when no
        bridges are supplied.
        """
        stale = [
            (u, v, k)
            for u, v, k, d in self._graph.edges(keys=True, data=True)
            if d.get("type") == "adjacency"
        ]
        for u, v, k in stale:
            self._graph.remove_edge(u, v, k)

        if bridges is None:
            unstable_pairs = self._consecutive_on_branch("unstable")
        else:
            unstable_pairs = [
                (u, v)
                for u, v in bridges
                if u is not None
                and v is not None
                and u != v
                and u in self._graph
                and v in self._graph
            ]

        for i, (u, v) in enumerate(unstable_pairs):
            self._graph.add_edge(
                u,
                v,
                key=f"adj_u_{i}",
                type="adjacency",
                stability="unstable",
            )

        for i, (u, v) in enumerate(self._consecutive_on_branch("stable")):
            self._graph.add_edge(
                u,
                v,
                key=f"adj_s_{i}",
                type="adjacency",
                stability="stable",
            )

        self._graph_bridges = bridges
        self._graph_adjacency_dirty = False

    def _consecutive_on_branch(
        self, stability: Literal["unstable", "stable"]
    ) -> list[tuple[int, int]]:
        """Neighbouring id pairs within each branch, in that branch's cdist order."""
        order = (
            self._unstable_order if stability == "unstable" else self._stable_order
        )
        per_branch: dict[ManifoldKey, list[int]] = {}
        for ix_id in order:  # already sorted by cdist, so each group inherits it
            ix = self._store[ix_id]
            key = ix.manifold_a_key if stability == "unstable" else ix.manifold_b_key
            if key is None:
                continue  # no branch identity: it belongs to no arc
            per_branch.setdefault(key, []).append(ix_id)

        return [
            (ids[i], ids[i + 1])
            for ids in per_branch.values()
            for i in range(len(ids) - 1)
        ]

    # ── iterate table ──────────────────────────────────────────────────────

    def register_iterate(self, source_id: int, n: int, target_id: int):
        """
        Record f^n(source) = target.

        Delegates to iterate_table and wires the directed edge into the graph.
        """
        self.iterate_table.register_iterate(source_id, n, target_id)
        self._generation += 1
        if source_id in self._graph and target_id in self._graph:
            self._graph.add_edge(
                source_id,
                target_id,
                key=f"iter_{source_id}_{n}",
                type="iterate",
                stability="unstable",
                n=n,
            )

    def iterate_orbit(
        self, start_id: int, max_len: Optional[int] = None
    ) -> list[int]:
        """
        Return an intersection's forward orbit — one cut point per stable branch.

        Chases the iterate table's n=1 links forward from ``start_id``
        (``start -> M(start) -> M^2(start) -> …``), stopping after one full lap of the
        orbit (an iterate landing on an already-visited stable branch), when the table
        has no further forward link, or at ``max_len`` points. For a period-p anchor
        this yields the start plus its p-1 iterates — the points that bound a period-p
        resonance zone.

        The table is populated by :meth:`TangleWorkbench.infer_iterates` (called
        automatically by ``compute_intersections``), so this just reads it; if the
        table has not been filled, only ``start_id`` is returned.

        Args:
            start_id: Registry id to start from (e.g. a strong pip).
            max_len: Hard cap on the number of points returned. Defaults to the
                anchor's ``k_value`` (one full lap of the orbit).

        Returns:
            Ordered list ``[start_id, M(start_id), …]``; always includes ``start_id``.
        """
        bkey = self._store[start_id].manifold_b_key
        if max_len is None:
            max_len = getattr(bkey[0], "k_value", 1) if bkey is not None else 1

        orbit = [start_id]
        seen_branches = {bkey}
        current = start_id
        while len(orbit) < max_len:
            nxt = self.iterate_table[current, 1]
            if nxt is None:
                break
            branch = self._store[nxt].manifold_b_key
            if branch in seen_branches:
                break
            orbit.append(nxt)
            seen_branches.add(branch)
            current = nxt
        return orbit

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

            lambda_u = self._get_lambda_u(ix, stability)
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

    def nearest_by_unstable_cdist(self, cdist: float) -> Optional[int]:
        """
        Return the ID of the intersection whose unstable_cdist is closest to cdist.

        Uses the sorted _unstable_order list for O(log n) positioning.

        Args:
            cdist: Target unstable cdist value.

        Returns:
            Matching ID, or None if the registry is empty.
        """
        if not self._unstable_order:
            return None
        pos = bisect.bisect_left(self._unstable_keys, cdist)
        candidates = []
        if pos < len(self._unstable_order):
            candidates.append(self._unstable_order[pos])
        if pos > 0:
            candidates.append(self._unstable_order[pos - 1])
        return min(candidates, key=lambda id_: abs(self._store[id_].unstable_cdist - cdist))

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

    def _get_lambda_u(
        self,
        intersection: Intersection,
        stability: Literal["unstable", "stable"] = "unstable",
    ) -> Optional[float]:
        """
        Read the unstable eigenvalue magnitude from the intersection's manifold keys.

        A crossing can involve two different fixed points (a heteroclinic
        crossing), whose eigenvalues differ, so the side matters: the image of a
        crossing along the unstable manifold is governed by the unstable
        branch's fixed point, and along the stable manifold by the stable
        branch's fixed point.

        Args:
            intersection: The crossing whose eigenvalue is wanted.
            stability: Which side to read: "unstable" uses ``manifold_a_key``,
                "stable" uses ``manifold_b_key``. Falls back to the other key
                (logging at debug) only when the requested one is None.

        Returns:
            The eigenvalue magnitude as a Python float, or None when no key is
            set or the fixed point has no eigenvalues yet (callers skip the
            intersection in that case).
        """
        preferred = (
            intersection.manifold_a_key
            if stability == "unstable"
            else intersection.manifold_b_key
        )
        key = preferred
        if key is None:
            key = (
                intersection.manifold_b_key
                if stability == "unstable"
                else intersection.manifold_a_key
            )
            if key is not None:
                logger.debug(
                    "no %s-side manifold key on intersection %s; falling back to "
                    "the other side's eigenvalue",
                    stability,
                    intersection.id,
                )
        if key is None:
            return None

        fp = key[0]
        eigenvalues = getattr(fp, "unstable_eigenvalues", None)
        if eigenvalues is None or len(eigenvalues) == 0:
            return None
        return float(abs(np.ravel(eigenvalues[0])[0]))

    def _find_collision(self, intersection: Intersection) -> Optional[int]:
        """
        The id of the stored intersection that is the same crossing, if any.

        A duplicate must match on both canonical distances (within cdist_tol) and
        on the two branches it crosses (manifold_a_key = unstable branch,
        manifold_b_key = stable branch). The branch test is what keeps the
        per-branch anchors apart: they all share cdist (0, 0) but lie on different
        branches, so a cdist-only test would collapse them into one. Branch keys
        are exact discrete labels, so this needs no positional tolerance and does
        not rely on the (approximate, drift-prone) phase-space coordinates.

        The unstable cdist is the prefilter: the sorted key array is bisected for
        the window ``[u - tol, u + tol]`` and only those candidates are tested in
        full, so an insert costs O(log N) plus the (normally empty) window rather
        than a scan of the whole store. When several stored crossings fall inside
        the window and match, the smallest id wins -- ids are handed out in
        increasing order and preserved by :meth:`reindex_from`, so that is the
        oldest of them and the deterministic choice.
        """
        tol = self.cdist_tol
        u = intersection.unstable_cdist
        lo = bisect.bisect_left(self._unstable_keys, u - tol)
        hi = bisect.bisect_right(self._unstable_keys, u + tol)

        best: Optional[int] = None
        for pos in range(lo, hi):
            candidate_id = self._unstable_order[pos]
            existing = self._store[candidate_id]
            if (
                abs(existing.unstable_cdist - u) < tol
                and abs(existing.stable_cdist - intersection.stable_cdist) < tol
                and existing.manifold_a_key == intersection.manifold_a_key
                and existing.manifold_b_key == intersection.manifold_b_key
            ):
                if best is None or candidate_id < best:
                    best = candidate_id
        return best

    def candidates_near_stable_cdist(
        self, branch_key: Optional[ManifoldKey], cdist: float, window: float
    ) -> list[tuple[int, Intersection]]:
        """
        The crossings on one stable branch whose stable cdist is within a window.

        The bucketed lookup behind the iterate search: the image of a crossing
        must sit on a KNOWN stable branch (the advanced key) at a KNOWN stable
        canonical distance (the source's, scaled by the per-step factor), so the
        candidate set is a bisected slice of that branch's own sorted cdists
        rather than the whole registry.

        Args:
            branch_key: The stable branch (``manifold_b_key``) to search.
            cdist: Centre of the stable-cdist window.
            window: Half-width of the window (inclusive on both ends).

        Returns:
            ``(id, Intersection)`` pairs in increasing stable cdist.
        """
        keys, ids = self._stable_branch_buckets().get(branch_key, ([], []))
        lo = bisect.bisect_left(keys, cdist - window)
        hi = bisect.bisect_right(keys, cdist + window)
        return [(ids[pos], self._store[ids[pos]]) for pos in range(lo, hi)]

    def _stable_branch_buckets(
        self,
    ) -> dict[Optional[ManifoldKey], tuple[list[float], list[int]]]:
        """Per-stable-branch ``(sorted cdists, ids)``, rebuilt when the orders change."""
        if (
            self._stable_branch_index is None
            or self._stable_branch_version != self._orders_version
        ):
            buckets: dict[
                Optional[ManifoldKey], tuple[list[float], list[int]]
            ] = {}
            for iid in self._stable_order:  # already in ascending stable cdist
                ix = self._store[iid]
                keys, ids = buckets.setdefault(ix.manifold_b_key, ([], []))
                keys.append(ix.stable_cdist)
                ids.append(iid)
            self._stable_branch_index = buckets
            self._stable_branch_version = self._orders_version
        return self._stable_branch_index

    def _insert_into_unstable_order(self, id: int, unstable_cdist: float):
        pos = bisect.bisect_left(self._unstable_keys, unstable_cdist)
        self._unstable_order.insert(pos, id)
        self._unstable_keys.insert(pos, unstable_cdist)

    def _insert_into_stable_order(self, id: int, stable_cdist: float):
        pos = bisect.bisect_left(self._stable_keys, stable_cdist)
        self._stable_order.insert(pos, id)
        self._stable_keys.insert(pos, stable_cdist)
