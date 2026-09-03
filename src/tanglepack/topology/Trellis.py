from __future__ import annotations

import logging
from typing import Iterable, Literal, Optional, Union, TYPE_CHECKING

import numpy as np
import matplotlib.pyplot as plt

from ..numerics.Intersection import Intersection, ManifoldKey
from ..numerics.IntersectionRegistry import IntersectionRegistry
from .TrellisBranch import TrellisBranch
from .TopologyResults import (
    Hole,
    PseudoneighborPair,
    StablePartitionResult,
    StrongPipResult,
)

if TYPE_CHECKING:
    from ..numerics.FixedPoint import FixedPoint
    from ..numerics.DynamicalSystem import DynamicalSystem
    from ..numerics.Bridge import Bridge
    from ..numerics.TangleWorkbench import TangleWorkbench

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

Stability = Literal["unstable", "stable"]

"""
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

This object deliberately contains NO algorithm logic. Compute-Pseudoneighbors,
Is-Strong-Pip, and friends will live in their own modules and read from / write
to a Trellis instance. The accessors here (branch lookup, orderings, cdist
scaling, next-intersection) are the shared primitives those algorithms compose.
"""


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
    ):
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
        Return the ID of M^n(intersection_id), or None if not recorded.

        Delegates to the registry's iterate table. n may be negative for
        backward iterates.
        """
        return self.registry.iterate_table[intersection_id, n]

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
            verbose: Print :meth:`describe_pseudoneighbors` when done.

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
            print(self.describe_pseudoneighbors())
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
        are stored in :attr:`holes`.

        Args:
            pairs: Pairs to punch holes for; defaults to every recorded pair.
            epsilon: Inward nudge of the hole off the manifold it hugs.
            propagate: Also punch the backward-propagated holes.
            verbose: Print :meth:`describe_holes` when done.

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
            print(self.describe_holes())
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
            verbose: Print :meth:`describe_stable_partitions` when done.

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
            print(self.describe_stable_partitions())
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

    def clear_results(self) -> None:
        """Empty all algorithm-output slots, leaving the input trellis intact."""
        self.pseudoneighbors.clear()
        self.holes.clear()
        self.stable_partitions.clear()
        self.strong_pip_candidates.clear()
        self.strong_pip = None

    # ── plotting ────────────────────────────────────────────────────────────

    def plot_strong_pip_candidates(self, ax=None, **scatter_kwargs):
        """
        Scatter-plot the strong-pip candidates in magenta, on top of the tangle.

        Draw this before plot_strong_pip() so the single chosen pip (green) sits on
        top of the candidate set (magenta) it was chosen from.

        Args:
            ax: Optional matplotlib Axes to draw on. Defaults to the current axes
                (plt).
            **scatter_kwargs: Forwarded to scatter. Defaults are color="magenta",
                s=7 (matching plot_intersections), zorder=11 (above the black
                intersections, below the green strong pip); override any of them
                (e.g. pass s=30 for larger markers).

        Returns:
            The matplotlib PathCollection, or None if there are no candidates.
        """
        if not self.strong_pip_candidates:
            logger.info(
                "No strong-pip candidates; call classify_strong_pips() first."
            )
            return None

        coords = np.array(
            [self.registry[iid].coords for iid in self.strong_pip_candidates]
        )
        scatter_kwargs.setdefault("color", "magenta")
        scatter_kwargs.setdefault("s", 7)
        scatter_kwargs.setdefault("zorder", 11)
        target = ax if ax is not None else plt
        return target.scatter(coords[:, 0], coords[:, 1], **scatter_kwargs)

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

    def plot_strong_pip(self, ax=None, **scatter_kwargs):
        """
        Scatter-plot the strong-pip cut points in green, on top of the tangle.

        For a period-1 anchor this is the single chosen strong pip; for a period-k
        anchor it is the k points (strong pip + iterates) that bound the resonance
        zone — see :meth:`strong_pip_cut_points`. Call after plotting the tangle and
        intersections to highlight them.

        Args:
            ax: Optional matplotlib Axes to draw on. Defaults to the current axes
                (plt).
            **scatter_kwargs: Forwarded to scatter. Defaults are color="green",
                s=7 (matching plot_intersections), zorder=12 (above the black
                intersections); override any of them (e.g. pass s=40 for a larger
                marker).

        Returns:
            The matplotlib PathCollection, or None if no strong pip has been chosen.
        """
        if self.strong_pip is None:
            logger.info(
                "No strong pip chosen; call classify_strong_pips() (or set_strong_pip()) first."
            )
            return None

        ids = self.strong_pip_cut_points()
        coords = np.array([self.registry[i].coords for i in ids])
        scatter_kwargs.setdefault("color", "green")
        scatter_kwargs.setdefault("s", 7)
        scatter_kwargs.setdefault("zorder", 12)
        target = ax if ax is not None else plt
        return target.scatter(coords[:, 0], coords[:, 1], **scatter_kwargs)

    def plot_pseudoneighbors(
        self, ax=None, *, include_trajectories: bool = True, **scatter_kwargs
    ):
        """
        Scatter-plot the pseudoneighbor pair members, on top of the tangle.

        Every branch has pseudoneighbors — only the reference computation is
        restricted to the fundamental segment — so by default all recorded
        pairs (references and their iterated appearances) are drawn.

        Args:
            ax: Optional matplotlib Axes to draw on. Defaults to the current
                axes (plt).
            include_trajectories: If True (default) every recorded pair is
                drawn; pass False for the reference pairs only.
            **scatter_kwargs: Forwarded to scatter. Defaults are
                color="darkorange", s=7 (matching plot_intersections),
                zorder=13 (above the strong-pip markers); override any of them.

        Returns:
            The matplotlib PathCollection, or None if there are no pairs.
        """
        pairs = (
            self.pseudoneighbors if include_trajectories
            else self.reference_pseudoneighbors
        )
        if not pairs:
            logger.info(
                "No pseudoneighbors to plot; call compute_pseudoneighbors() first."
            )
            return None

        ids = sorted({i for p in pairs for i in p.as_tuple()})
        coords = np.array([self.registry[i].coords for i in ids])
        scatter_kwargs.setdefault("color", "darkorange")
        scatter_kwargs.setdefault("s", 7)
        scatter_kwargs.setdefault("zorder", 13)
        target = ax if ax is not None else plt
        return target.scatter(coords[:, 0], coords[:, 1], **scatter_kwargs)

    # One marker per reference-pseudoneighbor orbit: every hole descending
    # from the same reference (matching Hole.origin) shares the symbol.
    _HOLE_MARKERS = ("x", "+", "*", "^", "s", "D", "v", "P")

    # One color per reference-pseudoneighbor orbit, matching the markers, so
    # an orbit reads as a single consistent shape+color across the figure.
    _HOLE_COLORS = (
        "purple", "teal", "darkgreen", "crimson",
        "chocolate", "navy", "olive", "deeppink",
    )

    def plot_holes(self, ax=None, *, show_iterates: bool = True, **scatter_kwargs):
        """
        Scatter-plot the punched holes, on top of the tangle.

        Each reference pseudoneighbor's orbit gets its own marker symbol AND
        color — the reference hole and every hole generated from it (forward,
        backward, or propagated) share them — and each hole is labelled with
        its iterate (0 = reference, negative = backward). Sides are reported
        by the partition, not encoded here.

        Args:
            ax: Optional matplotlib Axes to draw on. Defaults to the current
                axes (plt).
            show_iterates: Annotate each hole with its iterate number
                (default True).
            **scatter_kwargs: Forwarded to scatter. Defaults are s=40,
                zorder=14; override any of them. A color= override disables
                the by-orbit coloring, a marker= override the by-orbit markers.

        Returns:
            List of matplotlib PathCollections (one per orbit drawn), or None
            if there are no holes.
        """
        if not self.holes:
            logger.info("No holes to plot; call punch_holes() first.")
            return None

        scatter_kwargs.setdefault("s", 40)
        scatter_kwargs.setdefault("zorder", 14)
        origins = sorted({h.origin for h in self.holes if h.origin is not None})
        style_of = {
            origin: (
                self._HOLE_MARKERS[i % len(self._HOLE_MARKERS)],
                self._HOLE_COLORS[i % len(self._HOLE_COLORS)],
            )
            for i, origin in enumerate(origins)
        }
        target = ax if ax is not None else plt

        handles = []
        groups = sorted(
            {h.origin for h in self.holes},
            key=lambda origin: (origin is None, origin or ()),
        )
        for origin in groups:
            batch = [h for h in self.holes if h.origin == origin]
            coords = np.array([h.coords for h in batch])
            marker, color = style_of.get(origin, ("x", "gray"))
            kwargs = dict(scatter_kwargs)
            kwargs.setdefault("marker", marker)
            kwargs.setdefault("color", color)
            handles.append(target.scatter(coords[:, 0], coords[:, 1], **kwargs))
            if show_iterates:
                axes = target if ax is not None else plt.gca()
                for hole in batch:
                    if hole.iterate is None:
                        continue
                    axes.annotate(
                        str(hole.iterate),
                        hole.coords,
                        textcoords="offset points",
                        xytext=(4, 4),
                        fontsize=8,
                        color=color,
                        zorder=15,
                    )
        return handles

    def plot_stable_partition(self, ax=None, **line_kwargs):
        """
        Draw the stored stable partitions as number lines (one row each).

        See :func:`topology.StablePartition.plot_stable_partition` — the x-axis
        is stable canonical distance, closed endpoints are filled markers and
        open endpoints hollow ones.

        Args:
            ax: Optional matplotlib Axes to draw on. Defaults to the current
                axes.
            **line_kwargs: Forwarded to the interval plot calls.

        Returns:
            The Axes drawn on, or None if no partitions are stored.
        """
        from .StablePartition import plot_stable_partition as _plot

        return _plot(self.stable_partitions, ax=ax, **line_kwargs)

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
