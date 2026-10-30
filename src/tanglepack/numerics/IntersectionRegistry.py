"""
The master store of crossings and their iterate table.

:class:`IntersectionRegistry` assigns every crossing a contiguous integer id,
dedupes on (unstable cdist, stable cdist, unstable branch, stable branch), and
keeps the ordered per-branch views plus the
:class:`~.IterateTable.IterateTable` that maps an id to its iterates.

It also decides which detections are crossings of N (design note §A.2): every
detection is registered, and :meth:`IntersectionRegistry.resolved_ids` keeps the
ones the polylines resolve. A run of detections closer along the unstable
polyline than their chord deviation is judged over the WHOLE registry, because
on the blast path one run can straddle two ``iterate_bridge`` calls.
"""

from __future__ import annotations

import bisect
import logging
from collections import defaultdict
from typing import Callable, Iterator, Optional, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from .Intersection import Intersection, ManifoldKey, Stability
from .IterateTable import IterateTable

if TYPE_CHECKING:
    from .FixedPoint import FixedPoint

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

#: Two detections are the same crossing when their branch keys agree and both
#: canonical distances agree to this RELATIVE tolerance. A
#: re-detection is bit-identical in practice (the same polylines); distinct
#: crossings on one branch pair stay >= 4e-3 apart even at 11 k=2.8 blasts, while
#: an absolute 1e-6 would merge deep ones. An anchor (cdists exactly 0.0) matches
#: only another exact 0.0 at the same periodic point.
DEDUP_RTOL = 1e-9


#: Per stable branch: the sorted stable cdists, the ids, their brackets.
_StableBucket = tuple[list[float], NDArray[np.int64], NDArray[np.float64]]
_EMPTY_BUCKET: _StableBucket = ([], np.empty(0, dtype=np.int64), np.empty((0, 4)))


def _close(a: float, b: float) -> bool:
    """``a`` and ``b`` agree to :data:`DEDUP_RTOL`, relative to the larger."""
    return abs(a - b) <= DEDUP_RTOL * max(abs(a), abs(b))


def _polyline_gap(a: Intersection, b: Intersection) -> Optional[float]:
    """Distance from ``a`` to ``b`` along the unstable polyline, None unless adjacent.

    On one stored segment it is the straight distance; on two segments sharing a
    node (row equality: a store keeps one point per row, and a ``BranchPoint``
    end equals only itself) it runs through that node. The segments are the
    ones captured at detection, before any separator was spliced in; a
    separator adds a new row and moves none.
    """
    assert a.store is b.store, (a.id, b.id)
    shared = [end for end in a.unstable_segment if end in b.unstable_segment]
    if not shared:
        return None
    a_xy, b_xy = np.asarray(a.coords, dtype=float), np.asarray(b.coords, dtype=float)
    if len(shared) == 2:
        return float(np.hypot(*(a_xy - b_xy)))
    node = a.segment_ends("unstable")[0][a.unstable_segment.index(shared[0])]
    return float(np.hypot(*(a_xy - node)) + np.hypot(*(node - b_xy)))


class IntersectionRegistry:
    """
    Master store of all intersection points in the tangle.

    IDs are contiguous integers (0, 1, 2, …) assigned in insertion order.
    An intersection is uniquely identified by (unstable_cdist, stable_cdist,
    unstable branch, stable branch) — the two cdists plus the two manifold keys it
    crosses on. Duplicates (a match on all four, cdists equal to
    :data:`DEDUP_RTOL`) are detected and deduped — the existing ID is returned rather than creating a
    second entry. The branch keys are required because canonical distance is
    measured per branch: the anchor of every periodic point sits at cdist (0, 0),
    so cdist alone would wrongly merge the distinct anchors of a period-p orbit
    into one. A periodic point is ONE anchor, on every branch attached to it,
    so at exactly (0, 0) two keys match when they name the same fixed
    point and orbit index. Branch keys are used rather than phase-space coordinates because the
    coordinates are only approximate (drift under the nonlinear map), while the
    cdists and branch identity are exact/stable.

    Primary interface:
        registry.add(intersection)               → int (assigned ID)
        registry[id]                             → Intersection
        registry.iterate_table[id, n]            → int or None
        registry.by_unstable_cdist               → list[int]  (sorted by u-cdist)
        registry.by_stable_cdist                 → list[int]  (sorted by s-cdist)
        registry.resolved_ids(threshold)         → frozenset[int]  (N)

    Query interface (all return list[Intersection]):
        registry.from_fixed_point(fp)            → intersections involving fp
        registry.filter(predicate)               → arbitrary predicate

    Attributes:
        _store: dict[int, Intersection]
        _next_id: int
        cdist_tol: float  — an ABSOLUTE canonical-distance slack for geometry
            only: a degenerate bridge cut (``Tangle.create_bridges``,
            ``BridgeIterator``) and clipping a polyline to an arc. Never
            identity, ownership or containment (the topology reads none),
            and NOT the dedup rule, which is :data:`DEDUP_RTOL`
        iterate_table: IterateTable  — a partial bijection: the map is injective
        unresolved_links: dict[int, tuple[int, ...]]  — source id -> the candidate
            images the exact rule cannot decide between (several bracketed
            images, or one image claimed by several sources); none of them is
            linked, so the set is disjoint from the table's ``+1`` sources
        _generation: int  — bumped by every mutation; read via :attr:`generation`
        _orders_version: int  — bumped by the mutations that change the stored
            crossings (insert, renumber); the ordering-derived views key on it
        _unstable_order: list[int]  — IDs sorted ascending by unstable_cdist
        _unstable_keys: list[float] — those IDs' unstable cdists, same order
        _stable_order: list[int]    — IDs sorted ascending by stable_cdist
        _stable_keys: list[float]   — those IDs' stable cdists, same order

    Note:
        The two sorted key arrays are maintained on insertion (bisect), so no
        query ever rebuilds them by scanning the store. Every derived view
        (rank maps, the per-stable-branch candidate index) is memoised against
        ``_orders_version`` and rebuilt lazily on the next read after an insert
        or a renumbering.
    """

    def __init__(self, cdist_tol: float = 1e-6) -> None:
        """
        Build an empty registry.

        Args:
            cdist_tol (float): The absolute canonical-distance slack of
                degenerate bridge cuts and polyline clipping (never identity,
                and never the dedup rule; see :data:`DEDUP_RTOL`).
        """
        self._store: dict[int, Intersection] = {}
        self._next_id: int = 0
        self.cdist_tol = cdist_tol
        self.iterate_table = IterateTable()
        self.unresolved_links: dict[int, tuple[int, ...]] = {}

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
        self._stable_branch_index: Optional[dict[Optional[ManifoldKey], _StableBucket]] = None
        self._stable_branch_version: int = -1
        # Each detection's nearest neighbour and each run member's partner, and N
        # per threshold.
        self._runs_cache: tuple[dict[int, tuple[float, int]], dict[int, int]] = ({}, {})
        self._runs_version: int = -1
        self._resolved_cache: dict[Optional[float], frozenset[int]] = {}
        self._resolved_version: int = -1

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

        Rebuilding a derived cache (the rank maps, the candidate index) is not
        a mutation either: those are functions of the contents, so they
        never invalidate anything.

        Returns:
            int: The current generation.
        """
        return self._generation

    # ── core insert / lookup ───────────────────────────────────────────────

    def add(self, intersection: Intersection) -> int:
        """
        Register an intersection, returning its unique ID.

        If a collision is detected (the same two branches and cdists equal to
        :data:`DEDUP_RTOL`), the existing ID is returned and no duplicate is stored.

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
        with filter() and related methods.
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

    def find(self, intersection: Intersection) -> Optional[int]:
        """
        The id of the stored crossing that IS this one, if it is already here.

        The public face of the collision test :meth:`add` applies: a match needs
        both canonical distances equal to :data:`DEDUP_RTOL` AND both branch keys
        equal (at ``(0, 0)``: the same periodic point).
        Use it to resolve a crossing you can describe but do not hold -- the
        image of an anchor, say, which is the anchor of the next orbit point and
        is fully specified by its ``(0, 0)`` cdists plus its two advanced keys.

        Args:
            intersection: A probe describing the crossing to find. It is not
                stored and needs no id.

        Returns:
            The registry id of the matching crossing, or None.
        """
        return self._find_collision(intersection)

    def remap_rows(self, old_to_new: NDArray[np.int64]) -> None:
        """
        Renumber every crossing's bracket rows after the store was compacted.

        The rows keep their points (``PointStore.compact``), so no cdist, order
        or cache changes; only the row numbers the brackets name do.

        Args:
            old_to_new: The compaction's map; every bracket row survives it.
        """
        def remap(segment: Optional[tuple]) -> Optional[tuple]:
            """A bracket's rows renumbered; its ``BranchPoint`` ends pass through."""
            if segment is None:
                return None
            ends = tuple(int(old_to_new[end]) if isinstance(end, int) else end for end in segment)
            assert all(not isinstance(end, int) or end >= 0 for end in ends), "a bracket row was freed"
            return ends

        for crossing in self._store.values():
            crossing.unstable_segment = remap(crossing.unstable_segment)
            crossing.stable_segment = remap(crossing.stable_segment)

    def reindex_from(self, old: "IntersectionRegistry") -> None:
        """
        Re-number this registry so crossings carried over from ``old`` keep their id.

        IDs are normally insertion-order integers, so any rebuild of the registry
        (e.g. the recompute after trimming a stable manifold for a resonance zone)
        renumbers every crossing from zero. That invalidates ids the caller is
        holding — a strong pip chosen as id 7 becomes some unrelated integer. This
        re-aligns the freshly built registry against ``old``: every intersection
        that matches one in ``old`` (same canonical distances to :data:`DEDUP_RTOL`
        and the same two branch keys — the exact/stable identity used by
        :meth:`_find_collision`) is given back its old id. Crossings with no match
        in ``old`` (genuinely new geometry) are assigned fresh ids above every
        preserved one, so no id is ever reused for a different point.

        Call this after all intersections have been added but before the iterate
        table is filled (the table is keyed by id). The sorted orderings and the
        cdist index are rebuilt with the final ids.

        Args:
            old: The registry that existed before the rebuild, whose ids should be
                preserved wherever the same crossing reappears.

        Note:
            Nothing migrates id-keyed derived state across this call. A Trellis
            (and every result hanging off it, partition elements included) is a
            snapshot dropped and rebuilt when the workbench generation moves,
            which a renumbering does.
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
        self.iterate_table = IterateTable()
        self.unresolved_links = {}
        self._next_id = (max(used) + 1) if used else 0
        self._generation += 1
        self._orders_version += 1

        for ix in intersections:
            fid = final_id[id(ix)]
            ix.id = fid
            self._store[fid] = ix
            self._insert_into_unstable_order(fid, ix.unstable_cdist)
            self._insert_into_stable_order(fid, ix.stable_cdist)

    def __getitem__(self, id: int) -> Intersection:
        """registry[id] → Intersection. Raises KeyError if id is unknown."""
        return self._store[id]

    def __len__(self) -> int:
        return len(self._store)

    def __contains__(self, id: int) -> bool:
        return id in self._store

    def __iter__(self) -> Iterator[tuple[int, Intersection]]:
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

    def run_separation(self, id: int) -> Optional[float]:
        """
        How far a detection's nearest same-curve neighbour lies along the polyline.

        Args:
            id: A registered crossing.

        Returns:
            The distance along the unstable polyline to the nearest other
            detection of the same ``(u_key, s_key)`` that is its neighbour in
            unstable cdist order and sits on the same or a node-sharing unstable
            segment; None when there is none (always for a synthetic crossing).
        """
        nearest = self._runs()[0].get(id)
        return None if nearest is None else nearest[0]

    def run_partner(self, id: int) -> Optional[int]:
        """The nearest neighbour a detection forms a run with (clause 2), or None when in no run."""
        return self._runs()[1].get(id)

    def resolved_ids(self, threshold: Optional[float]) -> frozenset[int]:
        """
        N: the registered crossings the polylines resolve (design note §A.2).

        A crossing is in N when it is synthetic, or when (clause 1) its
        resolution ratio is at least ``threshold`` and (clause 2) it is in no
        run. Two neighbouring detections (:meth:`run_separation`) are a run
        when they sit no farther apart than the LARGER of their two chord
        deviations: the pair is judged as a pair, so both are out, and a run is
        out whole, median included. Every path that reads the registry as the
        crossing set reads this set.

        Args:
            threshold: The smallest resolution ratio a crossing may have
                (``TangleWorkbench.resolution_ratio_threshold``); None turns
                clause 1 off, and clause 2 still applies.

        Returns:
            The ids of N, cached until the stored crossings change.
        """
        if self._resolved_version != self._orders_version:
            self._resolved_cache, self._resolved_version = {}, self._orders_version
        if threshold not in self._resolved_cache:
            runs = self._runs()[1]
            self._resolved_cache[threshold] = frozenset(
                iid
                for iid, ix in self._store.items()
                if iid not in runs
                and (threshold is None or ix.resolution_ratio is None or ix.resolution_ratio >= threshold)
            )
        return self._resolved_cache[threshold]

    def unresolved_crossings(self, threshold: Optional[float]) -> list[int]:
        """
        The registered detections that are not crossings of N.

        Args:
            threshold: As in :meth:`resolved_ids`.

        Returns:
            The complement of :meth:`resolved_ids`, in unstable cdist order:
            the detections below ``threshold`` (clause 1) and every member of
            every run (clause 2). Anchors and synthetic crossings are never
            listed.
        """
        resolved = self.resolved_ids(threshold)
        return [iid for iid in self._unstable_order if iid not in resolved]

    def _runs(self) -> tuple[dict[int, tuple[float, int]], dict[int, int]]:
        """Per detection its nearest neighbour ``(separation, id)``, and per run member its partner.

        Rebuilt when the stored crossings change (an insert or a renumbering),
        never by an iterate link.
        """
        if self._runs_version != self._orders_version:
            groups: dict[tuple, list[Intersection]] = defaultdict(list)
            for iid in self._unstable_order:
                ix = self._store[iid]
                if ix.unstable_segment is not None:
                    groups[ix.manifold_a_key, ix.manifold_b_key].append(ix)
            nearest: dict[int, tuple[float, int]] = {}
            in_run: dict[int, tuple[float, int]] = {}
            for group in groups.values():
                for a, b in zip(group, group[1:]):
                    gap = _polyline_gap(a, b)
                    if gap is None:
                        continue
                    # A pair is judged as a pair: by the larger of its two deviations.
                    run = gap <= max(a.chord_deviation, b.chord_deviation)
                    for ix, other in ((a, b), (b, a)):
                        nearest[ix.id] = min(nearest.get(ix.id, (np.inf, -1)), (gap, other.id))
                        if run:
                            in_run[ix.id] = min(in_run.get(ix.id, (np.inf, -1)), (gap, other.id))
            partners = {iid: partner for iid, (_gap, partner) in in_run.items()}
            self._runs_cache, self._runs_version = (nearest, partners), self._orders_version
        return self._runs_cache

    def all_ids(self) -> list[int]:
        """All registered IDs in insertion order."""
        return list(self._store.keys())

    # ── iterate table ──────────────────────────────────────────────────────

    def register_iterate(self, source_id: int, n: int, target_id: int) -> None:
        """
        Record f^n(source) = target.

        Delegates to iterate_table, which records the backward twin too.

        Raises:
            AssertionError: If the link would overwrite another: the map is
                injective, so the table is a partial bijection.
        """
        assert self.iterate_table[source_id, n] in (None, target_id), (source_id, n)
        assert self.iterate_table[target_id, -n] in (None, source_id), (target_id, -n)
        self.iterate_table.register_iterate(source_id, n, target_id)
        self._generation += 1

    def remove_iterate(self, source_id: int, n: int) -> None:
        """
        Forget f^n(source) and its backward twin.

        Raises:
            KeyError: If the table holds no such link.
        """
        self.iterate_table.remove(source_id, n)
        self._generation += 1

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

        The table is filled by :meth:`TangleWorkbench.infer_iterates`, which
        ``compute_intersections``, every ``iterate_bridge`` and the end of every
        blast run, so this just reads it; if the table has not been filled, only
        ``start_id`` is returned.

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

    def from_fixed_point(
        self,
        fp: FixedPoint,
        stability: Optional[Stability] = None,
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

    # ── internal helpers ───────────────────────────────────────────────────

    def _find_collision(self, intersection: Intersection) -> Optional[int]:
        """
        The id of the stored intersection that is the same crossing, if any.

        A duplicate must match on the two branches it crosses (manifold_a_key =
        unstable branch, manifold_b_key = stable branch) and on both canonical
        distances to :data:`DEDUP_RTOL`. The branch test is what keeps the
        anchors of a period-p orbit apart: they all sit at cdist (0, 0) but on
        the branches of different orbit points. At exactly (0, 0) only the fixed
        point and orbit index of the two keys are compared, because a periodic
        point is one crossing of all its branches: the root segments of
        branch pairs other than the anchor's own resolve to the anchor. Branch
        keys are exact discrete labels, and the (approximate, drift-prone)
        phase-space coordinates are never read.

        The unstable cdist is the prefilter: the sorted key array is bisected for
        the relative window around it, so an insert costs O(log N) plus the
        (normally empty) window. When several stored crossings match, the
        smallest id wins -- ids are handed out in increasing order and preserved
        by :meth:`reindex_from`, so that is the oldest and the deterministic
        choice.
        """
        u, s = intersection.unstable_cdist, intersection.stable_cdist
        # A periodic point is ONE crossing of every branch attached to it:
        # at exactly (0, 0) any two of its branches meet at its one anchor.
        at_anchor = u == 0.0 and s == 0.0

        def branches(crossing: Intersection) -> tuple:
            keys = (crossing.manifold_a_key, crossing.manifold_b_key)
            if at_anchor and None not in keys:
                return tuple((key[0], key[2]) for key in keys)
            return keys

        own = branches(intersection)
        lo = bisect.bisect_left(self._unstable_keys, u - DEDUP_RTOL * abs(u))
        hi = bisect.bisect_right(self._unstable_keys, u + DEDUP_RTOL * abs(u))
        matches = (
            candidate_id
            for candidate_id in self._unstable_order[lo:hi]
            if _close(self._store[candidate_id].unstable_cdist, u)
            and _close(self._store[candidate_id].stable_cdist, s)
            and branches(self._store[candidate_id]) == own
        )
        return min(matches, default=None)

    def stable_brackets(self, branch_key: Optional[ManifoldKey]) -> _StableBucket:
        """
        The crossings of one stable branch, in stable-cdist order, with their brackets.

        The bucketed lookup behind the iterate search: the image of a crossing
        must sit on a KNOWN stable branch (the advanced key) at a KNOWN stable
        canonical distance (the source's, scaled by the per-step factor), so the
        candidates are a bisected slice of that branch's own sorted cdists
        rather than the whole registry.

        Args:
            branch_key: The stable branch (``manifold_b_key``).

        Returns:
            The stable cdists (a list: the bisect keys), the ids (int64) and
            their unpadded brackets ``(u_lo, u_hi, s_lo, s_hi)`` per row (an
            ``(m, 4)`` float array, NaN where a crossing has none). Empty for a
            branch no crossing is on.

        Note:
            Rebuilt when the stored crossings change (an insert or a
            renumbering), never by an iterate link; each crossing's brackets
            are read once, into ``Intersection.brackets``.
        """
        return self._stable_branch_buckets().get(branch_key, _EMPTY_BUCKET)

    def _stable_branch_buckets(self) -> dict[Optional[ManifoldKey], _StableBucket]:
        """Per stable branch its :meth:`stable_brackets`, rebuilt when the orders change."""
        if (
            self._stable_branch_index is None
            or self._stable_branch_version != self._orders_version
        ):
            members: dict[Optional[ManifoldKey], list[int]] = defaultdict(list)
            for iid in self._stable_order:  # already in ascending stable cdist
                members[self._store[iid].manifold_b_key].append(iid)
            nan = ((np.nan, np.nan), (np.nan, np.nan))
            self._stable_branch_index = {
                branch_key: (
                    [self._store[iid].stable_cdist for iid in ids],
                    np.array(ids, dtype=np.int64),
                    np.array(
                        [self._store[iid].brackets or nan for iid in ids], dtype=np.float64
                    ).reshape(-1, 4),
                )
                for branch_key, ids in members.items()
            }
            self._stable_branch_version = self._orders_version
        return self._stable_branch_index

    def _insert_into_unstable_order(self, id: int, unstable_cdist: float) -> None:
        pos = bisect.bisect_left(self._unstable_keys, unstable_cdist)
        self._unstable_order.insert(pos, id)
        self._unstable_keys.insert(pos, unstable_cdist)

    def _insert_into_stable_order(self, id: int, stable_cdist: float) -> None:
        pos = bisect.bisect_left(self._stable_keys, stable_cdist)
        self._stable_order.insert(pos, id)
        self._stable_keys.insert(pos, stable_cdist)
