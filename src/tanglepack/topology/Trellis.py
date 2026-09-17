"""
The topological snapshot of a computed tangle.

:class:`Trellis` is the entry point for the topological half of the library:
built from a :class:`~..numerics.TangleWorkbench.TangleWorkbench` once the
numerical phase is complete, it exposes the per-branch ordering of crossings
and the algorithms (strong pip, pseudoneighbors, stable partition) that run on
it, without reaching back into the linked-list manifold machinery.

Dev Notes:

Trellis is the entry point for the topological half of the library. It is built
from a TangleWorkbench once the numerical phase is complete (manifolds grown,
intersections computed, bridges cut, iterate table inferred) and from then on
exposes everything the topological algorithms need without reaching back into the
linked-list manifold machinery.

Snapshot semantics: from_workbench() captures the *derived* topological
structure (the per-branch ordering of intersections) at build time, but holds
live references to the registry, bridges, and fixed points. Growing manifolds or
recomputing intersections on the workbench afterward invalidates a Trellis —
rebuild it. The snapshot records the workbench generation it was built at
(``_built_generation``), which is what TangleSession.trellis() compares to
decide that, so nothing has to be invalidated by hand. Everything a Trellis
memoises on top of the snapshot (the endpoint-to-bridge index, the oriented
bridge polylines) is a function of the snapshot alone and dies with it.

The algorithm results are id-keyed (partition elements included) and so are the
snapshot's own branch orderings. A registry renumbering (``reindex_from``) moves
the workbench generation, so the whole Trellis -- results and all -- is dropped
and rebuilt rather than having its id-keyed tables remapped: remapping only the
partition maps would leave the branch orderings and bridge list stale, which is
the more dangerous half. Nothing therefore consumes ``reindex_from``'s remap.

Open items on the partition-element lookups (Phase 5), both unexercised by any
available fixture:

* ``image_of_element`` flips the requested side once per map step when the map
  reverses orientation (det J < 0). Every fixture in the repo is
  orientation-PRESERVING (b = 1), so that branch has never run against real
  data; the flip mirrors the rule ``bridge_side_violations`` applies to a hole's
  backward chain, and should be validated together with it on a det J < 0 map.
* An element bounded by the anchor (``lo_id is None``) or by the end of the
  computed branch (``hi_id is None``) has no image arc at all -- there is no
  crossing whose image would bound it -- so ``image_of_element`` answers None
  there. That is now the ONLY None it returns: an endpoint that IS a crossing
  but has no registered iterate is mapped by ``image_cdist``, which falls back
  from the iterate table to ``advance_key`` plus ``per_step_beta``. The
  deliberate anchor registration of 6.3 would give the anchor-bounded case a
  real id and close what is left. ``image_of_element(..., use_table=False)``
  forces the scaling branch for both endpoints; it exists so a test can pin the
  two answers as equal wherever the table has both iterates.

  That equality is not free. A canonical distance is an arc length on a
  polyline, so ``c * per_step_beta ** n`` reproduces the registered image only
  to about 1e-3 RELATIVE -- three orders above the registry's ``cdist_tol`` --
  and the image branch's partition boundaries sit exactly ON the registered
  images. A scaled endpoint therefore lands just past a boundary and would drag
  a neighbouring element into the answer, so ``image_of_element`` snaps a scaled
  endpoint onto a boundary it lands within ``SNAP_RTOL`` of before computing the
  cover. Should a fixture ever put two boundaries inside one window, the nearest
  wins and the ambiguity is logged at DEBUG.

  What the snap is and is not evidence for. On BOTH repository fixtures the snap
  never fires on the production path: every endpoint ``use_table=True`` resolves
  by scaling (6 of them on p3) lands with no boundary inside its window, so the
  answer there is the raw scaled span. The snap is exercised only by
  ``use_table=False``, where scaling is forced for endpoints the table could
  have answered and the boundary genuinely is the right target. So the
  table-versus-scaling equality test validates the scaling LAW; it does not
  validate the true fallback, which has no ground truth to be checked against --
  an unregistered iterate is unregistered precisely because nothing knows where
  it went. The real fallback is pinned by coverage alone (every bounded element
  gets a non-empty cover on the advanced branch) and stays that way until
  iterate inference reaches deeper than the current +/-1.

  ``_snap_to_partition_boundary`` rebuilds and re-sorts the image branch's
  boundary set on every call. That is nothing at fixture size (tens of
  intervals, two calls per element); if the dual graph ever makes element images
  hot, memoise the sorted boundaries per ``StablePartitionResult`` and bisect.

This object deliberately contains NO algorithm logic. Compute-Pseudoneighbors,
Is-Strong-Pip, and friends will live in their own modules and read from / write
to a Trellis instance. The accessors here (branch lookup, orderings, cdist
scaling, next-intersection) are the shared primitives those algorithms compose.
"""

from __future__ import annotations

import logging
from typing import Iterable, Optional, Union, TYPE_CHECKING

import numpy as np

from ..numerics.Intersection import Intersection, ManifoldKey, Stability
from ..numerics.IntersectionRegistry import IntersectionRegistry
from . import plotting
from .TrellisBranch import TrellisBranch
from .TopologyResults import (
    Endpoint,
    Hole,
    OPPOSITE_SIDE,
    PseudoneighborPair,
    Side,
    StablePartitionResult,
    StrongPipResult,
    endpoint_index,
)

if TYPE_CHECKING:
    from .Arrangement import Arrangement
    from ..numerics.FixedPoint import FixedPoint
    from ..numerics.DynamicalSystem import DynamicalSystem
    from ..numerics.Bridge import Bridge, BridgeId
    from ..numerics.TangleWorkbench import TangleWorkbench

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

#: Relative accuracy of the canonical-distance scaling law against the iterate
#: table. A canonical distance is an arc length measured on a POLYLINE, and two
#: branches of one orbit are refined independently, so ``c * per_step_beta ** n``
#: reproduces the recorded image only to the discretisation error of the two
#: curves -- 1.35e-3 at worst across the repository's fixtures, not the
#: registry's 1e-6 ``cdist_tol``. This is the AGREEMENT bound: what a test
#: comparing a derived distance against a registered one may demand, set an
#: order of magnitude above the observed error.
SCALING_RTOL = 1e-2

#: Half-width, relative to the distance itself, of the window inside which a
#: SCALED canonical distance is taken to BE a partition boundary
#: (:func:`_snap_to_partition_boundary`). Deliberately tighter than
#: :data:`SCALING_RTOL`: this one decides an element's identity, so it wants the
#: smallest window that still swallows the scaling error (max displacement
#: actually taken on the fixtures: 1.35e-3) rather than the loosest one a test
#: would tolerate. At 2e-3 the p3 and k=10 fixtures snap every forced-scaling
#: endpoint onto the boundary the iterate table names -- 0 disagreements out of
#: 56 and 2 respectively -- while the nearest WRONG boundary stays some 15x
#: outside the window.
SNAP_RTOL = 2e-3


def _snap_to_partition_boundary(
    cdist: float, result: StablePartitionResult, tol: float
) -> float:
    """Pull a scaled cdist onto a partition boundary it lands within its error of."""
    window = max(tol, abs(cdist) * SNAP_RTOL)
    boundaries = sorted(
        {iv.lo_cdist for iv in result.intervals}
        | {iv.hi_cdist for iv in result.intervals}
    )
    if not boundaries:
        return cdist
    inside = [b for b in boundaries if abs(b - cdist) <= window]
    if not inside:
        return cdist
    nearest = min(inside, key=lambda b: abs(b - cdist))
    if len(inside) > 1:
        logger.debug(
            "scaled cdist %.12g has %d partition boundaries within %.3g; "
            "snapping to the nearest (%.12g)",
            cdist,
            len(inside),
            window,
            nearest,
        )
    return nearest


class Trellis:
    """
    The finite trellis T = T^U ∩ T^S: the topological view of a computed tangle.

    A Trellis gathers the intersections, branches, bridges, and eigen/orbit data
    of one or more fixed points into a single object keyed for topological work,
    and provides storage slots for the outputs of the topological algorithms
    (pseudoneighbors, strong-pip classifications, holes).

    Build one with :meth:`from_workbench`.

    Attributes:
        fixed_points: The fixed points whose manifolds make up this trellis.
        registry: The master IntersectionRegistry (shared with the workbench).
        branches: Mapping from manifold key to TrellisBranch.
        bridges: Bridges belonging to the trellis's fixed points.
        dynamical_system: The underlying map/inverse (for algorithms that need it;
            many steps can instead use canonical-distance scaling — see scale_cdist).
        pseudoneighbors: Output slot — pseudoneighbor pairs found by the algorithm.
        holes: Output slot — holes punched for pseudoneighbor pairs.
        stable_partitions: Output slot — stable-manifold partitions built from
            the holes (one per branch and side).
        strong_pip_candidates: Output slot — IDs of every intersection that
            qualifies as a strong pip.
        strong_pip: The single chosen strong pip (an intersection ID), or None.
            Chosen from strong_pip_candidates via set_strong_pip().
    """

    def __init__(
        self,
        fixed_points: list["FixedPoint"],
        registry: IntersectionRegistry,
        branches: dict[ManifoldKey, TrellisBranch],
        bridges: list["Bridge"],
        dynamical_system: Optional["DynamicalSystem"] = None,
        manifolds: Optional[dict] = None,
        generation: int = 0,
    ) -> None:
        """
        Build a trellis from an already-computed tangle.

        Prefer :meth:`from_workbench`, which derives every argument below.

        Args:
            fixed_points (list[FixedPoint]): The orbits this trellis covers.
            registry (IntersectionRegistry): The crossings, and their iterates.
            branches (dict[ManifoldKey, TrellisBranch]): The per-branch ordering
                of those crossings.
            bridges (list[Bridge]): The bridges cut out of the unstable branches.
            dynamical_system (Optional[DynamicalSystem]): The map, when a
                geometric predicate needs it.
            manifolds (Optional[dict]): The workbench's manifold dict, for the
                algorithms that walk a real curve.
            generation (int): The workbench generation this snapshot was taken
                at; a cache compares it to decide staleness.
        """
        self.fixed_points = fixed_points
        self.registry = registry
        self.branches = branches
        self.bridges = bridges
        self.dynamical_system = dynamical_system
        # Live reference to the workbench's manifolds (keyed like branches).
        # Purely geometric consumers (hole placement inside narrow lobes) walk
        # the actual curve nodes through this; None degrades to chord fallbacks.
        self.manifolds = manifolds or {}
        # The workbench generation this snapshot was built at. A Trellis is a
        # snapshot -- the per-branch orderings and the bridge list are frozen at
        # build time -- so it is valid exactly while the workbench still reports
        # this generation; TangleSession.trellis compares the two.
        self._built_generation = generation

        # Memos, all functions of the snapshot alone (see the properties that
        # build them). The bridge lookup is rebuilt when `bridges` is replaced
        # or appended to; the polylines are keyed by bridge version.
        self._bridge_by_endpoints: Optional[dict[frozenset, "Bridge"]] = None
        self._bridge_lookup_len: int = -1
        self._bridge_polylines: dict = {}
        self._arrangement = None

        # ── algorithm output slots (filled by topological algorithms) ────────
        self.pseudoneighbors: list[PseudoneighborPair] = []
        self.holes: list[Hole] = []
        self.stable_partitions: list[StablePartitionResult] = []

        # Every intersection that *qualifies* as a strong pip is a candidate.
        # A trellis has exactly one actual strong pip — a unique choice made from
        # the candidates (default: smallest unstable cdist; user-overridable via
        # set_strong_pip()).
        self.strong_pip_candidates: list[int] = []
        self.strong_pip: Optional[int] = None

        # Cache for the orientation_preserving property (see below).
        self._orientation_preserving: Optional[bool] = None

    # ── map properties ──────────────────────────────────────────────────────

    @property
    def orientation_preserving(self) -> bool:
        """
        True when the underlying map preserves orientation (det J > 0).

        Area-preserving maps have a constant Jacobian determinant sign, so one
        evaluation settles it: the determinant is read at the first fixed
        point's first orbit coordinate when the dynamical system carries a
        jacobian. Without one, an odd-period fixed point still decides it —
        the product of its two eigenvalues is the determinant of the p-fold
        Jacobian, whose sign equals the single-step sign for odd p. Failing
        both, the answer defaults to True with a debug log.

        Orientation decides how the hole-side invariant reads a backward chain:
        an orientation-preserving map carries a bridge's dynamical orientation
        to its image's unchanged, while a reversing one flips the side at
        every step (see
        :func:`topology.StablePartition.bridge_side_violations`).

        Returns:
            True if the map preserves orientation, False if it reverses it.
        """
        if self._orientation_preserving is None:
            self._orientation_preserving = self._compute_orientation_preserving()
        return self._orientation_preserving

    def _compute_orientation_preserving(self) -> bool:
        """Sign of det J from the map, else from an odd-period eigenvalue pair."""
        jacobian = getattr(self.dynamical_system, "jacobian", None)
        coords = self.fixed_points[0].coordinates[0] if self.fixed_points else None
        if jacobian is not None and coords is not None:
            point = np.asarray(coords, dtype=np.float64).ravel()[:2]
            if np.isfinite(point).all():
                determinant = float(np.linalg.det(np.asarray(jacobian(point))))
                if np.isfinite(determinant) and determinant != 0.0:
                    return determinant > 0.0

        for fixed_point in self.fixed_points:
            if fixed_point.period % 2 == 0:
                continue  # even p squares away the single-step sign
            lam_u = getattr(fixed_point, "unstable_eigenvalues", None)
            lam_s = getattr(fixed_point, "stable_eigenvalues", None)
            if not lam_u or not lam_s:
                continue
            product = float(np.asarray(lam_u[0]).ravel()[0]) * float(
                np.asarray(lam_s[0]).ravel()[0]
            )
            if product != 0.0:
                return product > 0.0

        logger.debug(
            "No jacobian and no odd-period eigenvalue pair on this trellis; "
            "assuming the map preserves orientation"
        )
        return True

    # ── construction ────────────────────────────────────────────────────────

    @classmethod
    def from_workbench(
        cls,
        workbench: "TangleWorkbench",
        fixed_points: Optional[Union["FixedPoint", Iterable["FixedPoint"]]] = None,
    ) -> "Trellis":
        """
        Build a Trellis from a workbench whose numerical phase is complete.

        The caller is responsible for having grown the manifolds, computed
        intersections, created bridges, and (if needed) inferred the iterate
        table before calling this — the same state a finished script leaves
        behind.

        Args:
            workbench: The TangleWorkbench to snapshot.
            fixed_points: A single FixedPoint, an iterable of them, or None to
                use every fixed point registered on the workbench.

        Returns:
            A Trellis whose branches are populated and ordered, with empty
            algorithm-output slots.
        """
        if fixed_points is None:
            selected = list(workbench.fixed_points)
        elif _is_single_fixed_point(fixed_points):
            selected = [fixed_points]
        else:
            selected = list(fixed_points)

        selected_set = set(selected)
        registry = workbench.intersection_registry

        # One branch per workbench manifold belonging to a selected fixed point.
        branches: dict[ManifoldKey, TrellisBranch] = {}
        for key, _manifold in workbench.manifolds.items():
            fp, stability, orbit_index, branch_index = key
            if fp not in selected_set:
                continue
            branches[key] = TrellisBranch(
                key=key,
                fixed_point=fp,
                stability=stability,
                orbit_index=orbit_index,
                branch_index=branch_index,
                intersection_ids=[],
            )

        # Bucket each intersection onto the branch(es) it lies on. An
        # intersection's unstable side is manifold_a_key, its stable side is
        # manifold_b_key.
        for ix_id, ix in registry:
            if ix.manifold_a_key in branches:
                branches[ix.manifold_a_key].intersection_ids.append(ix_id)
            if ix.manifold_b_key in branches:
                branches[ix.manifold_b_key].intersection_ids.append(ix_id)

        # Order each branch by canonical distance (anchor outward).
        for branch in branches.values():
            attr = (
                "unstable_cdist" if branch.stability == "unstable" else "stable_cdist"
            )
            branch.intersection_ids.sort(key=lambda i: getattr(registry[i], attr))

        bridges = [
            b
            for b in workbench.bridges
            if getattr(b, "fixed_point", None) in selected_set
        ]

        logger.debug(
            "Built Trellis: %d fixed points, %d branches, %d intersections, %d bridges",
            len(selected),
            len(branches),
            len(registry),
            len(bridges),
        )

        return cls(
            fixed_points=selected,
            registry=registry,
            branches=branches,
            bridges=bridges,
            dynamical_system=workbench.dynamical_system,
            manifolds=workbench.manifolds,
            generation=workbench.generation,
        )

    # ── snapshot memos ──────────────────────────────────────────────────────

    def bridge_between(self, first_id: int, second_id: int) -> Optional["Bridge"]:
        """
        The bridge whose two endpoint crossings are exactly this pair.

        A dict lookup on a per-trellis index built from :attr:`bridges`, so the
        callers that resolve one bridge per pseudoneighbor pair (there is one
        per pair, per iterate, per backward step) no longer rescan the bridge
        list each time. The index is rebuilt when the bridge list changes
        length -- a Trellis is a snapshot, so that is the only way its bridges
        move without a rebuild of the whole object.

        Args:
            first_id: One endpoint's registry id.
            second_id: The other endpoint's registry id.

        Returns:
            The matching Bridge, or None if no bridge spans the pair.
        """
        if (
            self._bridge_by_endpoints is None
            or self._bridge_lookup_len != len(self.bridges)
        ):
            index: dict[frozenset, "Bridge"] = {}
            for bridge in self.bridges:
                if bridge.first_intersection is None or bridge.second_intersection is None:
                    continue  # a partial arc has no endpoint pair to be found by
                index.setdefault(
                    frozenset((bridge.first_intersection, bridge.second_intersection)),
                    bridge,
                )
            self._bridge_by_endpoints = index
            self._bridge_lookup_len = len(self.bridges)
        return self._bridge_by_endpoints.get(frozenset((first_id, second_id)))

    @property
    def arrangement(self) -> "Arrangement":
        """
        The planar arrangement of this trellis: its arcs and the faces they bound.

        Built on first access and kept for the life of the trellis. That IS caching
        by generation: a Trellis is a snapshot valid only while the workbench still
        reports :attr:`_built_generation`, so the session drops the whole object —
        arrangement included — the moment anything moves. Nothing has to be
        invalidated by hand.

        Prefer the ALL-fixed-points trellis (``session.arrangement()``): a
        heteroclinic crossing belongs to two tangles, and a single-fixed-point
        snapshot cuts the other side's arcs off, turning real faces into open ones.

        Returns:
            The :class:`~tanglepack.topology.Arrangement.Arrangement`.
        """
        if self._arrangement is None:
            from .Arrangement import Arrangement

            self._arrangement = Arrangement.from_trellis(self)
        return self._arrangement

    def cached_bridge_polyline(self, bridge: "Bridge", build):
        """
        Memoise one oriented polyline per bridge, keyed by the bridge's version.

        The geometry callers (hole placement, row reading, side tests) ask for
        the same bridge's oriented polyline many times per partition; building
        it walks the bridge and orients it. The memo is keyed by the bridge's
        own mutation counter (:attr:`BaseManifold.version`), so it follows the
        same invalidation as the point array underneath it.

        Args:
            bridge: The bridge whose polyline is wanted.
            build: Zero-argument callable that computes it on a miss. Its
                result -- including ``None`` for an unorientable bridge -- is
                what gets cached.

        Returns:
            Whatever ``build`` returned, from the memo when it is still valid.

        Note:
            An object with no ``version`` (a hand-built stand-in in a test) is
            never cached: without the counter there is nothing to invalidate
            against, and a stale polyline is worse than an extra walk.
        """
        version = getattr(bridge, "version", None)
        if version is None:
            return build()

        key = id(bridge)
        entry = self._bridge_polylines.get(key)
        if entry is not None and entry[0] is bridge and entry[1] == version:
            return entry[2]
        value = build()
        self._bridge_polylines[key] = (bridge, version, value)
        return value

    # ── intersection access ─────────────────────────────────────────────────

    def intersection(self, intersection_id: int) -> Intersection:
        """Return the Intersection with the given registry ID."""
        return self.registry[intersection_id]

    @property
    def intersection_ids(self) -> list[int]:
        """All intersection IDs in the trellis (insertion order)."""
        return self.registry.all_ids()

    @property
    def own_intersection_ids(self) -> list[int]:
        """
        IDs of the intersections that actually belong to this trellis.

        The registry is shared across every fixed point on the workbench, so in a
        nested / multi-tangle session it holds intersections this trellis does not
        own. An intersection belongs here iff it lies on one of this trellis's
        branches (its unstable or stable side resolves to a branch we hold), which
        is precisely the set ``from_workbench`` bucketed onto the branches. The
        topological algorithms default to this set so that, e.g., the inner
        period-3 trellis is not classified against the outer fixed point's
        crossings.

        Returns:
            Sorted list of owned intersection IDs.
        """
        ids: set[int] = set()
        for branch in self.branches.values():
            ids.update(branch.intersection_ids)
        return sorted(ids)

    @property
    def by_stable_cdist(self) -> list[int]:
        """All intersection IDs sorted ascending by stable canonical distance."""
        return self.registry.by_stable_cdist

    @property
    def by_unstable_cdist(self) -> list[int]:
        """All intersection IDs sorted ascending by unstable canonical distance."""
        return self.registry.by_unstable_cdist

    def iterate(self, intersection_id: int, n: int) -> Optional[int]:
        """
        The ID of ``M^n(intersection_id)``, composing single steps if need be.

        The registry's iterate table stores only what inference actually
        recorded, and inference records the links it can see — on the
        repository's fixtures that is the ``n = +/-1`` links and little else.
        A chain ``q -> M(q) -> M^2(q) -> M^3(q)`` therefore exists link by link
        while ``table[q, 3]`` is empty. So the direct entry is preferred where
        the table has one and, failing that, ``|n|`` single steps are composed
        in the direction of ``sign(n)``.

        Args:
            intersection_id: Registry ID of the crossing.
            n: Number of MAP STEPS; negative walks backward. ``0`` is the
                identity.

        Returns:
            The registry ID of the image, or ``None`` as soon as a link of the
            chain is missing.

        Note:
            The composed answer is a genuine table answer, not an estimate:
            every step is a recorded iterate relation. For the estimate used
            when the table has no chain at all, see :meth:`image_cdist`.
        """
        table = self.registry.iterate_table
        direct = table[intersection_id, n]
        if direct is not None or n == 0:
            return direct
        step = 1 if n > 0 else -1
        current = intersection_id
        for _ in range(abs(n)):
            current = table[current, step]
            if current is None:
                return None
        return current

    # ── branch access ───────────────────────────────────────────────────────

    def branch(self, key: ManifoldKey) -> Optional[TrellisBranch]:
        """Return the TrellisBranch for a manifold key, or None."""
        return self.branches.get(key)

    @property
    def stable_branches(self) -> list[TrellisBranch]:
        """All stable branches in the trellis."""
        return [b for b in self.branches.values() if b.stability == "stable"]

    @property
    def unstable_branches(self) -> list[TrellisBranch]:
        """All unstable branches in the trellis."""
        return [b for b in self.branches.values() if b.stability == "unstable"]

    def branches_of(
        self,
        fixed_point: Optional["FixedPoint"] = None,
        stability: Optional[Stability] = None,
    ) -> list[TrellisBranch]:
        """
        Return branches filtered by fixed point and/or stability.

        Args:
            fixed_point: If given, only branches anchored to this fixed point.
            stability: If given, only branches of this stability.

        Returns:
            List of matching TrellisBranch objects.
        """
        result = []
        for b in self.branches.values():
            if fixed_point is not None and b.fixed_point is not fixed_point:
                continue
            if stability is not None and b.stability != stability:
                continue
            result.append(b)
        return result

    def branch_containing(
        self, intersection_id: int, stability: Stability
    ) -> Optional[TrellisBranch]:
        """
        Return the branch of the given stability that an intersection lies on.

        Reads the intersection's manifold_a_key (unstable) or manifold_b_key
        (stable) and looks up the corresponding branch.

        Args:
            intersection_id: Registry ID of the intersection.
            stability: Which side's branch to resolve.

        Returns:
            The TrellisBranch, or None if the side key is unset or absent.
        """
        ix = self.registry[intersection_id]
        key = ix.manifold_a_key if stability == "unstable" else ix.manifold_b_key
        if key is None:
            return None
        return self.branches.get(key)

    def next_intersection(
        self,
        intersection_id: int,
        stability: Stability,
        toward_anchor: bool = True,
    ) -> Optional[int]:
        """
        Return the neighbouring intersection along a manifold branch.

        This is the next-intersection primitive used by the topological
        algorithms: walking along the manifold of the given stability, in the
        direction toward (default) or away from the anchoring periodic point.

        Args:
            intersection_id: Registry ID of the reference intersection.
            stability: Which manifold to walk along.
            toward_anchor: Direction of travel.

        Returns:
            Registry ID of the neighbour, or None at the end of the branch.
        """
        branch = self.branch_containing(intersection_id, stability)
        if branch is None:
            return None
        return branch.neighbor(intersection_id, toward_anchor=toward_anchor)

    # ── canonical-distance helpers ──────────────────────────────────────────

    def lambda_u(self, fixed_point: "FixedPoint") -> Optional[float]:
        """
        Unstable eigenvalue magnitude for a fixed point, governing cdist scaling.

        Returns None if the fixed point has no eigenvalues set.
        """
        evals = getattr(fixed_point, "unstable_eigenvalues", None)
        if not evals:
            return None
        return float(abs(np.asarray(evals[0]).ravel()[0]))

    def scale_cdist(
        self,
        cdist: float,
        n: int,
        stability: Stability,
        fixed_point: "FixedPoint",
    ) -> Optional[float]:
        """
        Scale a canonical distance under n applications of the MAP.

        Uses the eigenvalue relation rather than the dynamical map itself:
            unstable: c_dist(M^n) = c_dist · beta^n
            stable:   c_dist(M^n) = c_dist / beta^n
        where ``beta = fixed_point.per_step_beta("unstable")`` is the factor of
        ONE map step. Pass a negative n for backward iterates (the Strong Pip
        Algorithm maps intersections back onto a stable branch this way).

        Args:
            cdist: The starting canonical distance.
            n: Number of MAP STEPS (negative for backward). One full return to
                the same branch is ``k_value`` steps, so a caller counting branch
                returns passes ``returns * fixed_point.k_value`` -- which scales
                by ``|lambda| ** num_branches`` per return, i.e. by the
                eigenvalue itself only when the point has no inversion.
            stability: Which manifold's cdist is being scaled.
            fixed_point: Fixed point supplying the per-step factor.

        Returns:
            The scaled canonical distance, or None if the eigenvalue is
            unavailable.

        Note:
            The unit is map steps, not branch returns: ``advance_key``,
            ``per_step_beta`` and this method all count the same thing, so a
            crossing's image after ``n`` steps is on ``advance_key(key, n)`` at
            ``scale_cdist(c, n, ...)``.
        """
        if self.lambda_u(fixed_point) is None:
            return None
        beta = fixed_point.per_step_beta("unstable")
        if stability == "unstable":
            return cdist * (beta ** n)
        return cdist / (beta ** n)

    def _keyed_cdist(
        self, intersection_id: int, stability: Stability
    ) -> tuple[ManifoldKey, float]:
        """The branch key and canonical distance of one crossing on one side."""
        crossing = self.registry[intersection_id]
        if stability == "unstable":
            key, cdist = crossing.manifold_a_key, crossing.unstable_cdist
        elif stability == "stable":
            key, cdist = crossing.manifold_b_key, crossing.stable_cdist
        else:
            raise ValueError(
                f"stability must be 'unstable' or 'stable', got {stability!r}"
            )
        if key is None or cdist is None:
            raise ValueError(
                f"intersection {intersection_id} has no {stability} branch key / "
                f"canonical distance; it cannot be mapped on that side"
            )
        return key, float(cdist)

    def _scaled_image_cdist(
        self, intersection_id: int, n: int, stability: Stability
    ) -> tuple[ManifoldKey, float]:
        """Image key and cdist of a crossing from advance_key/per_step_beta alone."""
        key, cdist = self._keyed_cdist(intersection_id, stability)
        fixed_point = key[0]
        return (
            fixed_point.advance_key(key, n),
            cdist * fixed_point.per_step_beta(stability) ** n,
        )

    def image_cdist(
        self,
        intersection_id: int,
        n: int = 1,
        stability: Stability = "stable",
    ) -> tuple[ManifoldKey, float, bool]:
        """
        Where one crossing lands after ``n`` map steps, on one side.

        A crossing sits on two branches at once, and this answers for the one of
        the requested stability: the branch its image lies on and the canonical
        distance it lies at. The registry's iterate table is the authority
        wherever it is populated; otherwise the answer is derived, from
        :meth:`FixedPoint.advance_key` for the branch and
        :meth:`FixedPoint.per_step_beta` for the distance --
        ``cdist * per_step_beta(stability) ** n``, which contracts a stable
        distance and expands an unstable one, forward or backward.

        The derived answer is exact for the true crossing -- canonical distance
        is the linearising coordinate the eigenvalue scales -- but the two cdists
        it is computed from are arc lengths on independently refined polylines,
        so it reproduces the table's value only to :data:`SCALING_RTOL`
        (about 1e-3 on the repository's fixtures), NOT to the registry's
        ``cdist_tol``. Callers that compare a derived distance against a
        registered one must budget for that; :meth:`image_of_element` does.

        Args:
            intersection_id: Registry ID of the crossing.
            n: Number of MAP STEPS; negative walks backward. ``0`` is the
                identity. Defaults to 1.
            stability: Which of the crossing's two branches to follow --
                ``"stable"`` (the default) reads ``manifold_b_key`` and the
                stable cdist, ``"unstable"`` reads ``manifold_a_key`` and the
                unstable cdist.

        Returns:
            ``(branch_key, cdist, from_table)``: the image's manifold key, its
            canonical distance on that branch, and whether the answer came from
            the iterate table (``True``) or from the scaling law (``False``).
            ``n = 0`` answers with the crossing's own key and cdist and
            ``from_table=True``.

        Raises:
            ValueError: If the crossing carries no branch key (or no canonical
                distance) on the requested side, or if ``stability`` is neither
                "unstable" nor "stable".

        Note:
            This is deliberately NOT :meth:`scale_cdist`, which divides by the
            UNSTABLE per-step factor for a stable distance. The two agree on an
            area-preserving map, where ``per_step_beta("stable")`` is defined as
            the reciprocal of the unstable one, and this method stays right when
            they are two independent estimates.

        Note:
            Never compare the product of a crossing's stable and unstable
            distances to decide that two crossings are iterates: two different
            chains generally share that product. Compare the two sides
            individually (this method, once per stability).
        """
        key, cdist = self._keyed_cdist(intersection_id, stability)
        if n == 0:
            return key, cdist, True

        image_id = self.iterate(intersection_id, n)
        if image_id is not None:
            image_crossing = self.registry[image_id]
            image_key = (
                image_crossing.manifold_a_key
                if stability == "unstable"
                else image_crossing.manifold_b_key
            )
            image_cdist_value = (
                image_crossing.unstable_cdist
                if stability == "unstable"
                else image_crossing.stable_cdist
            )
            if image_key is not None and image_cdist_value is not None:
                return image_key, float(image_cdist_value), True
            logger.debug(
                "iterate %d of intersection %d has no %s branch key; "
                "falling back to canonical-distance scaling",
                image_id,
                intersection_id,
                stability,
            )

        image_key, image_cdist_value = self._scaled_image_cdist(
            intersection_id, n, stability
        )
        return image_key, image_cdist_value, False

    # ── result storage ──────────────────────────────────────────────────────

    def add_pseudoneighbor(self, pair: PseudoneighborPair) -> None:
        """Record a pseudoneighbor pair (and its hole, if attached)."""
        self.pseudoneighbors.append(pair)
        if pair.hole is not None:
            self.holes.append(pair.hole)

    @property
    def reference_pseudoneighbors(self) -> list[PseudoneighborPair]:
        """The recorded pairs found on the reference window W^S(r_n, r_{n+p})."""
        return [p for p in self.pseudoneighbors if p.is_reference]

    @property
    def strong_pip_intersection(self) -> Optional[Intersection]:
        """The chosen strong pip as an Intersection, or None if unset."""
        return None if self.strong_pip is None else self.registry[self.strong_pip]

    # ── topological algorithms ──────────────────────────────────────────────

    def is_strong_pip(
        self,
        intersection_id: int,
        *,
        tol: Optional[float] = None,
        collision_rtol: float = 1e-2,
    ) -> StrongPipResult:
        """
        Test whether one intersection qualifies as a strong pip (pure query).

        Thin wrapper around :func:`topology.StrongPip.is_strong_pip`. Does not
        change trellis state — use classify_strong_pips() to populate candidates
        and set_strong_pip() to choose the actual strong pip.

        Args:
            intersection_id: Registry ID of the candidate.
            tol: Optional canonical-distance slack (defaults to the registry's
                cdist_tol).
            collision_rtol: Relative slack on the cdist collision test (default
                1e-2) so that q0's own orbit does not disqualify it.

        Returns:
            The StrongPipResult.
        """
        from .StrongPip import is_strong_pip as _is_strong_pip

        return _is_strong_pip(
            self, intersection_id, tol=tol, collision_rtol=collision_rtol
        )

    def classify_strong_pips(
        self,
        intersection_ids: Optional[Iterable[int]] = None,
        *,
        tol: Optional[float] = None,
        collision_rtol: float = 1e-2,
        choose_default: bool = True,
    ) -> list[int]:
        """
        Find every intersection that qualifies as a strong pip (the *candidates*).

        Many intersections can satisfy the strong-pip condition, but a trellis has
        exactly one actual strong pip — a unique choice. This method collects all
        qualifying candidates into ``self.strong_pip_candidates`` and, by default,
        chooses the candidate with the smallest unstable canonical distance as the
        actual strong pip (``self.strong_pip``). Override that choice with
        set_strong_pip(), or pass ``choose_default=False`` to leave it unset.

        Args:
            intersection_ids: Intersections to test; defaults to this trellis's own
                intersections (see :attr:`own_intersection_ids`), not the whole
                shared registry.
            tol: Optional canonical-distance slack.
            collision_rtol: Relative slack on the cdist collision test (default
                1e-2) so that a point's own orbit does not disqualify it.
            choose_default: If True (default), also select the default strong pip
                (smallest unstable cdist) from the candidates.

        Returns:
            The list of candidate intersection IDs.
        """
        from .StrongPip import classify_strong_pips as _classify_strong_pips

        # Default to *this* trellis's intersections, not the whole shared registry,
        # so a nested trellis is classified only against its own tangle's crossings.
        if intersection_ids is None:
            intersection_ids = self.own_intersection_ids

        results = _classify_strong_pips(
            self, intersection_ids, tol=tol, collision_rtol=collision_rtol
        )
        self.strong_pip_candidates = sorted(
            iid for iid, r in results.items() if r.is_strong_pip
        )
        self.strong_pip = None
        if choose_default and self.strong_pip_candidates:
            self.select_default_strong_pip()
        return list(self.strong_pip_candidates)

    def set_strong_pip(self, intersection_id: int) -> int:
        """
        Choose a specific candidate as the trellis's one actual strong pip.

        Args:
            intersection_id: A registry ID that must be in
                ``self.strong_pip_candidates``.

        Returns:
            The chosen intersection ID.

        Raises:
            ValueError: If the ID is not among the classified candidates. Run
                classify_strong_pips() first, or pick from strong_pip_candidates.
        """
        if intersection_id not in self.strong_pip_candidates:
            raise ValueError(
                f"Intersection {intersection_id} is not a strong-pip candidate. "
                f"Candidates: {self.strong_pip_candidates}. "
                "Run classify_strong_pips() first, or choose one of these."
            )
        self.strong_pip = intersection_id
        return intersection_id

    def select_default_strong_pip(self) -> Optional[int]:
        """
        Choose the candidate with the smallest unstable cdist as the strong pip.

        This is the default selection criterion. It is a thin policy built on
        set_strong_pip() — write your own selection (e.g. by stable cdist, or any
        other rule) the same way and call set_strong_pip() with the winner.

        Returns:
            The chosen strong-pip ID, or None if there are no candidates.
        """
        if not self.strong_pip_candidates:
            return None
        chosen = min(
            self.strong_pip_candidates,
            key=lambda iid: self.registry[iid].unstable_cdist,
        )
        return self.set_strong_pip(chosen)

    def compute_pseudoneighbors(
        self,
        *,
        extend: bool = True,
        collision_rtol: float = 1e-2,
        match_rtol: float = 0.05,
        tol: Optional[float] = None,
        verbose: bool = False,
    ) -> list[PseudoneighborPair]:
        """
        Find the reference pseudoneighbor pairs and record them on the trellis.

        Thin wrapper around :func:`topology.Pseudoneighbor.compute_pseudoneighbors`.
        Clears any previously recorded pseudoneighbors and holes first. With
        ``extend=True`` (default) the reference pairs are also mapped through
        the iterate table and the full trajectories recorded alongside them.

        The reference window starts at the chosen strong pip's cut point on
        each branch (run classify_strong_pips() / set_strong_pip() first);
        without one, the outermost intersection is used as a fallback.

        Args:
            extend: Also record the iterated (non-reference) pairs.
            collision_rtol: Relative slack of the endpoint-collision test (a
                landing skipped only when BOTH its cdists match an endpoint).
            match_rtol: Relative tolerance of the r_{n+p} cdist fallback.
            tol: Absolute canonical-distance slack (defaults to the registry's
                ``cdist_tol``).
            verbose: Report :meth:`describe_pseudoneighbors` when done —
                logged at INFO (this module's logger is raised to INFO for the
                call, so a higher application-level setting cannot swallow it),
                or printed when logging is unconfigured.

        Returns:
            The reference pairs (also available as
            :attr:`reference_pseudoneighbors`).
        """
        from .Pseudoneighbor import (
            compute_pseudoneighbors as _compute,
            extend_pseudoneighbor_trajectories as _extend,
        )

        self.pseudoneighbors.clear()
        self.holes.clear()
        references = _compute(
            self, collision_rtol=collision_rtol, match_rtol=match_rtol, tol=tol
        )
        for pair in references:
            self.add_pseudoneighbor(pair)
        if extend:
            for pair in _extend(self, references):
                self.add_pseudoneighbor(pair)
        if verbose:
            _log_report(self.describe_pseudoneighbors())
        return references

    def describe_pseudoneighbors(self) -> str:
        """Human-readable report of the recorded pseudoneighbor pairs."""
        references = self.reference_pseudoneighbors
        lines = [f"{len(references)} reference pseudoneighbor pair(s):"]
        for pair in references:
            a = self.registry[pair.intersection_a]
            b = self.registry[pair.intersection_b]
            lines.append(
                f"  ({pair.intersection_a}, {pair.intersection_b})  "
                f"stable cdists ({a.stable_cdist:.4g}, {b.stable_cdist:.4g})  "
                f"unstable cdists ({a.unstable_cdist:.4g}, {b.unstable_cdist:.4g})"
            )
        lines.append(
            f"{len(self.pseudoneighbors)} pair(s) total including trajectories"
        )
        return "\n".join(lines)

    def punch_holes(
        self,
        pairs: Optional[Iterable[PseudoneighborPair]] = None,
        *,
        epsilon: float = 0.05,
        propagate: bool = True,
        verbose: bool = False,
    ) -> list[Hole]:
        """
        Punch the holes for the recorded pseudoneighbor pairs.

        Thin wrapper around :func:`topology.StablePartition.punch_holes`; with
        ``propagate=True`` (default) the reference bridges are also mapped
        backward and their generated holes punched (see
        :func:`topology.StablePartition.propagate_reference_holes`). All holes
        are stored in :attr:`holes`. Holes propagate backward only: a recorded
        pair at a forward iterate ``>= k_value`` gets no hole (its region is
        the image of one still attached to the stable manifold).

        Args:
            pairs: Pairs to punch holes for; defaults to every recorded pair.
            epsilon: Inward nudge of the hole off the manifold it hugs.
            propagate: Also punch the backward-propagated holes.
            verbose: Report :meth:`describe_holes` when done —
                logged at INFO (this module's logger is raised to INFO for the
                call, so a higher application-level setting cannot swallow it),
                or printed when logging is unconfigured.

        Returns:
            The punched holes.

        Raises:
            AssertionError: If two holes of one ``origin`` disagree on their
                side of the bridge (I1), or if a bridge carrying a hole
                approaches its two defining crossings from opposite sides of
                their shared stable branch (I2).

        Note:
            Holes only exist after punching, so this is where the two
            topological invariants are checked — once per orbit and once per
            bridge, not once per hole. Only the bridges that actually carry a
            hole are checked; the rest are not part of this result.
        """
        from .StablePartition import (
            check_bridge_rows_consistent as _check_rows,
            check_holes_share_bridge_side as _check_sides,
            punch_holes as _punch,
            propagate_reference_holes as _propagate,
        )

        self._warn_missing_pseudoneighbors()
        self.holes.clear()
        # holes is cleared wholesale, so every pair's back-reference goes with
        # it: a pair outside this call's scope must not keep pointing at a hole
        # the trellis no longer carries (propagate_reference_holes seeds its
        # already-punched orbits and regions from pair.hole).
        for pair in self.pseudoneighbors:
            pair.hole = None
        holes = _punch(self, pairs, epsilon=epsilon)
        if propagate:
            holes += _propagate(self)
        self.holes.extend(holes)

        punched = {
            frozenset(hole.bounding_ids)
            for hole in self.holes
            if hole.bounding_ids is not None
        }
        for bridge in self.bridges:
            ends = frozenset((bridge.first_intersection, bridge.second_intersection))
            if ends not in punched:
                continue
            _check_rows(self, bridge)
        _check_sides(
            self.holes, orientation_preserving=self.orientation_preserving
        )

        if verbose:
            _log_report(self.describe_holes())
        return holes

    def describe_holes(self) -> str:
        """Human-readable report of the punched holes, by bridge side."""
        lines = []
        for side in ("left", "right"):
            count = sum(1 for h in self.holes if h.bridge_side == side)
            lines.append(f"{count} hole(s) on the {side} side of their bridge")
        unclassified = sum(1 for h in self.holes if h.bridge_side is None)
        if unclassified:
            lines.append(f"{unclassified} hole(s) with no classified side")
        return "\n".join(lines)

    def partition_stable_manifold(
        self,
        branch_key: Optional[ManifoldKey] = None,
        *,
        verbose: bool = False,
    ) -> list[StablePartitionResult]:
        """
        Partition stable branch(es) by the punched holes, one result per side.

        Thin wrapper around
        :func:`topology.StablePartition.partition_stable_manifold`, run for
        both sides of each selected branch. Results are stored in
        :attr:`stable_partitions` (replacing previous ones for those branches).

        Args:
            branch_key: A single stable branch to partition, or None (default)
                for every stable branch of the trellis.
            verbose: Report :meth:`describe_stable_partitions` when done —
                logged at INFO (this module's logger is raised to INFO for the
                call, so a higher application-level setting cannot swallow it),
                or printed when logging is unconfigured.

        Returns:
            The partition results (two per branch: left and right).
        """
        from .StablePartition import partition_stable_manifold as _partition

        self._warn_missing_pseudoneighbors()
        keys = (
            [branch_key]
            if branch_key is not None
            else [b.key for b in self.stable_branches]
        )
        results = [
            _partition(self, key, side) for key in keys for side in ("left", "right")
        ]
        self.stable_partitions = [
            p for p in self.stable_partitions if p.branch_key not in keys
        ] + results
        if verbose:
            _log_report(self.describe_stable_partitions())
        return results

    def _warn_missing_pseudoneighbors(self) -> None:
        """Warn for every fixed point with no recorded pseudoneighbor pairs.

        Hole punching and partitioning cover every fixed point held by this
        trellis automatically, but only from the pairs already recorded — a
        fixed point whose pseudoneighbors were never computed silently
        contributes nothing, so flag it.
        """
        covered = {
            pair.branch_key[0]
            for pair in self.pseudoneighbors
            if pair.branch_key is not None
        }
        for fp in self.fixed_points:
            if fp not in covered:
                logger.warning(
                    "No pseudoneighbors recorded for %r; run "
                    "compute_pseudoneighbors() first — holes and partitions "
                    "will be empty for it",
                    fp,
                )

    def describe_stable_partitions(self) -> str:
        """Human-readable report of the stored partitions, interval notation."""
        lines = []
        for result in self.stable_partitions:
            lines.append(
                f"{result.side} partition (p{result.branch_key[0].period}, "
                f"orbit {result.branch_key[2]}):"
            )
            for iv in result.intervals:
                lo = "[" if iv.closed_lo else "("
                hi = "]" if iv.closed_hi else ")"
                lo_name = "anchor" if iv.lo_id is None else f"id {iv.lo_id}"
                hi_name = "end" if iv.hi_id is None else f"id {iv.hi_id}"
                singleton = (
                    "  (singleton)" if iv.lo_cdist == iv.hi_cdist else ""
                )
                lines.append(
                    f"  {lo}{iv.lo_cdist:.4g}, {iv.hi_cdist:.4g}{hi}  "
                    f"({lo_name} → {hi_name}){singleton}"
                )
        return "\n".join(lines)

    # ── partition elements ──────────────────────────────────────────────────

    def element_for(
        self,
        bridge_id: "BridgeId",
        endpoint: Endpoint,
        side: Side,
    ) -> Optional[int]:
        """
        The partition element at one end of a bridge, on one side.

        Scans this trellis's :attr:`stable_partitions` for the result whose
        branch carries that endpoint. A bridge's two endpoints need not lie on
        the same stable branch (on a period-q orbit the bridge spanning an
        anchor has one end on each of two branches of the cycle), which is why
        this is a scan over results rather than a lookup in one of them.

        A :data:`~tanglepack.numerics.Bridge.BridgeId` IS the ordered pair of its
        endpoints' registry ids, so the answer never depends on this trellis
        holding the bridge object: when
        :attr:`StablePartitionResult.elements_at_bridge` has no entry for the id
        — the bridge belongs to another fixed point, or was simply not in the
        snapshot's bridge list — the endpoint id is looked up directly in
        :attr:`StablePartitionResult.element_of_intersection`. That is what makes
        a heteroclinic bridge resolvable: a bridge of W^u(fp1) crossing W^s(fp3)
        is filtered out of fp3's snapshot (``from_workbench`` keeps bridges by
        their UNSTABLE fixed point) while the element that owns its endpoint
        lives precisely in fp3's partition.

        Args:
            bridge_id: The bridge's :data:`~tanglepack.numerics.Bridge.BridgeId`.
            endpoint: ``"first"`` or ``"second"`` — which end of the id.
            side: Which side's partition to read.

        Returns:
            The element id within the result for that endpoint's branch, or
            ``None`` when no partition of this trellis covers that endpoint —
            the branch was not partitioned, or belongs to another trellis, which
            :meth:`~tanglepack.loom.TangleSession.TangleSession.partition_element_for`
            searches.

        Raises:
            ValueError: If ``endpoint`` is neither "first" nor "second".

        Note:
            The element id alone only names an element within one result; pair
            it with the branch key and side (or read
            :attr:`~.TopologyResults.PartitionInterval.branch_key`) to identify
            the element globally.
        """
        index = endpoint_index(endpoint)
        for result in self.stable_partitions:
            if result.side != side:
                continue
            ends = result.elements_at_bridge.get(bridge_id)
            if ends is not None and ends[index] is not None:
                return ends[index]
            # No entry: the bridge object is not in this snapshot. Its id still
            # names the two crossings, so ask the branch that owns one of them.
            direct = result.element_of_intersection.get(bridge_id[index])
            if direct is not None:
                return direct
        return None

    def image_of_element(
        self,
        result: StablePartitionResult,
        element_id: int,
        n: int = 1,
        *,
        use_table: bool = True,
        partitions: Optional[Iterable[StablePartitionResult]] = None,
    ) -> Optional[list[int]]:
        """
        The elements covering the ``n``-th image of one partition element.

        An element spans the stable arc between two crossings, so its image
        spans the arc between their ``n``-images on the branch ``n`` map steps
        forward (:meth:`FixedPoint.advance_key`). Each endpoint's image is
        located by :meth:`image_cdist`: the registry's iterate table where it is
        populated, and the canonical-distance scaling law where it is not, so an
        element whose iterates were never registered still has an image. The
        answer is the elements of the image branch's partition that meet that
        span. Side is carried through unchanged by an orientation-preserving map
        and flips once per step otherwise.

        The image partition is not a refinement of the image of this partition:
        the hole orbits are propagated backward for finitely many steps, so the
        image branch can be missing an innermost boundary that this branch has.
        The answer is therefore every element of the image branch that MEETS the
        image arc as a set (open ends included: two pieces touching at a single
        cdist meet only when both are closed there), not only the ones the arc
        contains. Their combined span always covers both endpoint images.

        Args:
            result: The partition result the element belongs to.
            element_id: The element's id within ``result``.
            n: Number of map steps; negative walks backward. Defaults to 1.
            use_table: When False, BOTH endpoint images are taken from the
                scaling branch of :meth:`image_cdist` even where the iterate
                table has an entry. This is a validation hook: the two answers
                must agree on every element whose iterates are registered, and
                a test that pins that is what licenses the fallback elsewhere.
                Defaults to True (table first, per endpoint).
            partitions: Where to look for the image branch's partition result.
                Defaults to this trellis's own :attr:`stable_partitions`; pass
                an explicit collection when the trellis holding the geometry is
                not the one holding the partitions (the all-fixed-points trellis
                carries no partitions of its own, so the dual graph hands it the
                per-fixed-point results it was built from).

        Returns:
            The image element ids in increasing stable canonical distance, or
            ``None`` when an end of the element is unbounded (the anchor or the
            branch end), which is the only case with no image arc to cover. An
            empty list means the image arc lies outside the partitioned stretch
            of the image branch (possible walking backward, which moves outward).

        Raises:
            IndexError: If ``element_id`` is not an element of ``result``.
            ValueError: If the image branch's partition (same side, advanced
                key) is not among ``partitions`` — partition every stable branch
                before asking for images.

        Note:
            An endpoint image that came from the SCALING law is snapped onto a
            partition boundary of the image branch when it lands within
            :data:`SNAP_RTOL` of one (:func:`_snap_to_partition_boundary`). The
            image branch's boundaries are themselves images of crossings, so a
            scaled endpoint a discretisation error away from one IS that
            boundary, and without the snap that ~1e-3 relative error would spill
            the arc into a neighbouring element.

            On both repository fixtures the snap never fires on the default
            ``use_table=True`` path: no endpoint that falls back to scaling there
            has a boundary inside its window, so those answers are the raw scaled
            span. It fires only under ``use_table=False``, where scaling is
            forced for endpoints the table could have answered — which is exactly
            the comparison that pins the scaling law, and is NOT evidence about
            the true fallback. An unregistered iterate has no ground truth to be
            checked against; the genuine fallback is pinned by coverage alone
            until iterate inference reaches past +/-1.
        """
        from .StablePartition import owns_cdist, span_contains

        interval = result.element(element_id)
        if n == 0:
            return [element_id]
        if interval.lo_id is None or interval.hi_id is None:
            return None

        image_ends: list[tuple[float, bool]] = []
        for end_id in (interval.lo_id, interval.hi_id):
            if use_table:
                _key, cdist, from_table = self.image_cdist(end_id, n, "stable")
                if not from_table:
                    logger.debug(
                        "element %d of %s: crossing %d has no registered %d-iterate; "
                        "scaling its stable cdist instead",
                        element_id,
                        result.branch_key[1:],
                        end_id,
                        n,
                    )
            else:
                _key, cdist = self._scaled_image_cdist(end_id, n, "stable")
                from_table = False
            image_ends.append((cdist, from_table))

        fixed_point = result.branch_key[0]
        image_key = fixed_point.advance_key(result.branch_key, n)
        side = result.side
        if not self.orientation_preserving and n % 2:
            side = OPPOSITE_SIDE[side]

        candidates = (
            self.stable_partitions if partitions is None else list(partitions)
        )
        image_result = next(
            (
                other
                for other in candidates
                if other.branch_key == image_key and other.side == side
            ),
            None,
        )
        if image_result is None:
            raise ValueError(
                f"no {side} partition stored for the image branch "
                f"{image_key[1:]}; partition it before asking for images"
            )

        tol = self.registry.cdist_tol
        image_cdists = [
            cdist
            if from_table
            else _snap_to_partition_boundary(cdist, image_result, tol)
            for cdist, from_table in image_ends
        ]
        ends = sorted(
            (
                (image_cdists[0], interval.closed_lo),
                (image_cdists[1], interval.closed_hi),
            ),
            key=lambda end: end[0],
        )
        (lo, closed_lo), (hi, closed_hi) = ends

        covering: list[int] = []
        for candidate in image_result.intervals:
            overlap_lo = max(lo, candidate.lo_cdist)
            overlap_hi = min(hi, candidate.hi_cdist)
            if overlap_lo > overlap_hi + tol:
                continue
            if overlap_hi - overlap_lo > tol:
                covering.append(candidate.element_id)
                continue
            # The two only touch: they meet iff both are closed at that point.
            point = 0.5 * (overlap_lo + overlap_hi)
            if span_contains(lo, hi, closed_lo, closed_hi, point, tol) and owns_cdist(
                candidate, point, tol
            ):
                covering.append(candidate.element_id)
        return covering

    def clear_results(self) -> None:
        """Empty all algorithm-output slots, leaving the input trellis intact."""
        self.pseudoneighbors.clear()
        self.holes.clear()
        self.stable_partitions.clear()
        self.strong_pip_candidates.clear()
        self.strong_pip = None

    def strong_pip_cut_points(self) -> list[int]:
        """
        The strong pip together with its iterates — one cut point per stable branch.

        For a period-1 anchor this is just the strong pip. For a period-k anchor the
        resonance zone is bounded by the strong pip *and its k-1 forward iterates*,
        which land one on each stable branch; this returns all of them (via
        ``registry.iterate_orbit``, capped at the fixed point's ``k_value``). Requires
        the iterate table to be inferred; otherwise only the strong pip is returned.

        Returns:
            Ordered list of intersection ids, or an empty list if no strong pip is set.
        """
        if self.strong_pip is None:
            return []
        key = self.registry[self.strong_pip].manifold_b_key
        max_len = getattr(key[0], "k_value", None) if key is not None else None
        return self.registry.iterate_orbit(self.strong_pip, max_len=max_len)

    # ── plotting ────────────────────────────────────────────────────────────
    #
    # Every drawing routine lives in :mod:`topology.plotting`, which owns the
    # shared marker/colour/z-order conventions; these methods are the trellis-
    # bound aliases of those functions.

    def plot_strong_pip_candidates(self, ax=None, **scatter_kwargs):
        """
        Scatter-plot the strong-pip candidates in magenta, on top of the tangle.

        Delegates to :func:`topology.plotting.plot_strong_pip_candidates`; draw
        this before :meth:`plot_strong_pip` so the single chosen pip (green)
        sits on top of the candidate set (magenta) it was chosen from.

        Args:
            ax: Optional matplotlib Axes to draw on. Defaults to the current axes
                (plt).
            **scatter_kwargs: Forwarded to scatter, overriding
                :data:`~topology.plotting.STRONG_PIP_CANDIDATE_STYLE`.

        Returns:
            The matplotlib PathCollection, or None if there are no candidates.
        """
        return plotting.plot_strong_pip_candidates(self, ax=ax, **scatter_kwargs)

    def plot_strong_pip(self, ax=None, **scatter_kwargs):
        """
        Scatter-plot the strong-pip cut points in green, on top of the tangle.

        Delegates to :func:`topology.plotting.plot_strong_pip`. For a period-1
        anchor this is the single chosen strong pip; for a period-k anchor it is
        the k points (strong pip + iterates) that bound the resonance zone — see
        :meth:`strong_pip_cut_points`.

        Args:
            ax: Optional matplotlib Axes to draw on. Defaults to the current axes
                (plt).
            **scatter_kwargs: Forwarded to scatter, overriding
                :data:`~topology.plotting.STRONG_PIP_STYLE`.

        Returns:
            The matplotlib PathCollection, or None if no strong pip has been chosen.
        """
        return plotting.plot_strong_pip(self, ax=ax, **scatter_kwargs)

    def plot_pseudoneighbors(
        self, ax=None, *, include_trajectories: bool = True, **scatter_kwargs
    ):
        """
        Scatter-plot the pseudoneighbor pair members, on top of the tangle.

        Delegates to :func:`topology.plotting.plot_pseudoneighbors`. Every branch
        has pseudoneighbors — only the reference computation is restricted to the
        fundamental segment — so by default all recorded pairs (references and
        their iterated appearances) are drawn.

        Args:
            ax: Optional matplotlib Axes to draw on. Defaults to the current
                axes (plt).
            include_trajectories: If True (default) every recorded pair is
                drawn; pass False for the reference pairs only.
            **scatter_kwargs: Forwarded to scatter, overriding
                :data:`~topology.plotting.PSEUDONEIGHBOR_STYLE`.

        Returns:
            The matplotlib PathCollection, or None if there are no pairs.
        """
        return plotting.plot_pseudoneighbors(
            self, ax=ax, include_trajectories=include_trajectories, **scatter_kwargs
        )

    def plot_holes(self, ax=None, *, show_iterates: bool = True, **scatter_kwargs):
        """
        Scatter-plot the punched holes, on top of the tangle.

        Delegates to :func:`topology.plotting.plot_holes`: each reference
        pseudoneighbor's orbit gets its own marker symbol AND color, and each
        hole is labelled with its iterate (0 = reference, negative = backward).
        Sides are reported by the partition, not encoded here.

        Args:
            ax: Optional matplotlib Axes to draw on. Defaults to the current
                axes (plt).
            show_iterates: Annotate each hole with its iterate number
                (default True).
            **scatter_kwargs: Forwarded to scatter, overriding
                :data:`~topology.plotting.HOLE_STYLE`.

        Returns:
            List of matplotlib PathCollections (one per orbit drawn), or None
            if there are no holes.
        """
        return plotting.plot_holes(
            self, ax=ax, show_iterates=show_iterates, **scatter_kwargs
        )

    def plot_stable_partition(self, ax=None, **line_kwargs):
        """
        Draw the stored stable partitions as number lines (one row each).

        Delegates to :func:`topology.plotting.plot_stable_partition`, passing
        :attr:`stable_partitions` — the x-axis is stable canonical distance,
        closed endpoints are filled markers and open endpoints hollow ones.

        Args:
            ax: Optional matplotlib Axes to draw on. Defaults to the current
                axes.
            **line_kwargs: Forwarded to the interval plot calls.

        Returns:
            The Axes drawn on, or None if no partitions are stored.
        """
        return plotting.plot_stable_partition(
            self.stable_partitions, ax=ax, **line_kwargs
        )

    # ── misc ────────────────────────────────────────────────────────────────

    def summary(self) -> str:
        """Return a one-line human-readable summary of the trellis contents."""
        return (
            f"Trellis: {len(self.fixed_points)} fixed point(s), "
            f"{len(self.branches)} branch(es), "
            f"{len(self.registry)} intersection(s), "
            f"{len(self.bridges)} bridge(s) | "
            f"{len(self.pseudoneighbors)} pseudoneighbor(s), "
            f"{len(self.strong_pip_candidates)} strong-pip candidate(s), "
            f"chosen strong pip: {self.strong_pip}"
        )

    def __repr__(self) -> str:
        return f"<{self.summary()}>"


def _is_single_fixed_point(obj) -> bool:
    """True if obj is a single FixedPoint rather than an iterable of them."""
    from ..numerics.FixedPoint import FixedPoint

    return isinstance(obj, FixedPoint)


def _log_report(report: str) -> None:
    """
    Emit a ``describe_*`` report for a ``verbose=True`` call.

    ``verbose`` is an alias for "report this call", so the text has to reach the
    user either way:

    * With logging CONFIGURED (the root logger has a handler), the report goes to
      the log stream at INFO — this module's logger is raised to INFO for the
      duration and restored in a ``finally``, so an application level above INFO
      does not swallow it.
    * With logging UNCONFIGURED, ``logging.lastResort`` only emits at WARNING, so
      an INFO record would vanish; the report is printed instead, which is what a
      plain REPL or notebook caller of ``verbose=True`` expects.

    Args:
        report: The already-formatted report text.
    """
    if not logging.getLogger().hasHandlers():
        print(report)
        return
    previous = logger.level
    logger.setLevel(logging.INFO)
    try:
        logger.info(report)
    finally:
        logger.setLevel(previous)
