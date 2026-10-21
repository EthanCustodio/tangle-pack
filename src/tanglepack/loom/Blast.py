"""
Blasting — driving the dynamics inside a resonance zone.

"Blasting" a resonance zone repeatedly iterates the bridges that live *inside* it. One
blast step takes every un-iterated bridge whose INNERMOST zone is this one, maps it
forward one step under the dynamical map (which automatically cuts the image into child
bridges at its new crossings with the stable manifold), then keeps only the children
whose innermost zone is *still* this one. Those survivors are the next step's input, so
the process recurses up to ``num_iterations`` times.

Zones nest: a bridge belongs to the smallest zone containing its midpoint
(:meth:`TangleSession.classify_bridge`), so an outer zone's blast never iterates a
bridge inside a smaller zone of another fixed point stored on the session, and nested
blasting does not depend on the order of the zones. The blasted zone itself need not
be stored (:func:`define_resonance_zone`, or a replaced handle). With ``fixed_point=None`` an outer blast still iterates an inner tangle's bridges
lying in the annulus between the two zones: they are the outer zone's bridges. The
recipes pass the zone's own fixed point, which leaves them alone.

Bridges that leave the zone are deliberately dropped from the frontier and never
iterated again: an unstable manifold grows by the per-step eigenvalue factor each
iterate, so iterating an *exterior* bridge many times would blow up exponentially. The
zone interior is the only region where iterating many times stays bounded, and is the
intended use of this feature.

A bridge whose forward map fails aborts the blast: every such
failure is a bug.

Memory is budgeted at two levels. The machine's per-call cap
(``max_points_per_step``) is always on and RAISES
:class:`~tanglepack.numerics.budget.PointBudgetExceeded` mid-image, leaving the
session not resumable. ``blast_zone(max_points=...)`` is the opt-in budget of a
whole blast: it stops BETWEEN iterations, before one whose predicted cost would
take the blast over budget, so every completed iteration is whole and the session
can be blasted again (``stop_reason == "point_budget"``). The prediction is the
next frontier's point count times the last observed growth factor (points the
images added per frontier point), from this call or -- stored on the session per
zone -- from the previous one; blast factors are stable (~1.7 per step on k=2.8),
unlike growth past escape. The budget counts points, not RSS: blasting holds
roughly 1.7 KB per point, growth roughly 1.1 KB.

A blast ends with a registry-wide ``infer_iterates``, so the iterate table is
exact after it, whatever the stop reason (a raised ``PointBudgetExceeded``
excepted: that state is not resumable).

This is a loom-layer algorithm: it reads the topological notion of a zone and acts on
the numerical layer (:meth:`TangleWorkbench.iterate_bridge`). It is exposed on the
facade as :meth:`tanglepack.loom.TangleSession.blast_zone`.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal, Optional

import numpy as np
from scipy.spatial import cKDTree

from ..numerics.budget import PointBudgetExceeded
from .ResonanceZone import ResonanceZone

if TYPE_CHECKING:
    from ..numerics.Bridge import Bridge, BridgeId
    from ..numerics.FixedPoint import FixedPoint
    from .TangleSession import TangleSession

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def _seen_key(bridge: "Bridge") -> "BridgeId | int":
    """The cycle-guard key of a bridge in one blast.

    A bridge IS the two crossings it connects, so its ``BridgeId`` is the right
    identity: the workbench keeps exactly one object per id, and a child that
    re-traces an existing bridge comes back as that object. A partial piece
    (bounded by fewer than two crossings) has no id, so it falls back to object
    identity -- two partials over the same arc are not recognised as one, which is
    by definition, not an oversight.
    """
    return bridge.id if bridge.id is not None else id(bridge)


def _interior_points(bridge: "Bridge") -> np.ndarray:
    """The bridge's points with the outer 10% trimmed off each end.

    Adjacent bridges legitimately meet at their shared straddle endpoints, so the
    proximity test ignores the ends and only compares the *interiors*: two bridges
    are a precision problem when their bodies run nearly coincident, not when they
    merely touch at a shared crossing.
    """
    pts = np.asarray(bridge.get_point_array())
    if pts.ndim != 2 or len(pts) < 5:
        return pts if pts.ndim == 2 else np.empty((0, 2))
    margin = max(1, len(pts) // 10)
    return pts[margin:-margin]


def _min_interior_distance(child: "Bridge", reference: np.ndarray) -> float:
    """Smallest distance from ``child``'s interior to a reference point cloud."""
    pts = _interior_points(child)
    if len(pts) == 0 or len(reference) == 0:
        return float("inf")
    dists, _ = cKDTree(reference).query(pts)
    return float(np.min(dists))


def _pip_orbits(
    session: "TangleSession", fixed_points: list["FixedPoint"]
) -> dict["FixedPoint", list[int]]:
    """``[q0, f(q0), ..., f^m(q0)]`` per fixed point with a chosen pip: its registered orbit.

    Composed link by link, as ``Trellis._pip_orbit_cuts`` reads the depth ``m`` of
    ``T(q0)``: the iterate table holds ``+1`` links only, so
    ``IterateTable.forward_depth`` is 1 for any linked crossing.
    """
    registry = session.workbench.intersection_registry
    orbits = {}
    for fp in fixed_points:
        chosen = session.chosen_strong_pip(fp)
        x = registry.find(chosen) if chosen is not None else None
        if x is None:
            continue
        orbits[fp] = [x]
        while (x := registry.iterate_table[x, 1]) is not None:
            orbits[fp].append(x)
    return orbits


@dataclass
class BlastStep:
    """
    Record of one blast iteration.

    Attributes:
        iteration: 1-based index of this step within the blast.
        interior_parents: The un-iterated interior bridges fed into this step.
        produced_children: Every child bridge produced by iterating the parents.
        kept_interior: The children that landed back inside the zone — the next step's
            input.
        discarded_exterior: The children that fell outside the zone (dropped, never
            re-iterated).
    """

    iteration: int
    interior_parents: list["Bridge"]
    produced_children: list["Bridge"]
    kept_interior: list["Bridge"]
    discarded_exterior: list["Bridge"]
    discarded_too_close: list["Bridge"] = field(default_factory=list)
    already_known: list["Bridge"] = field(default_factory=list)


@dataclass
class BlastResult:
    """
    The genealogy produced by :func:`blast_zone`.

    Attributes:
        zone: The resonance zone that was blasted.
        fixed_points: The fixed points whose bridges took part (every fixed
            point of the workbench when the call passed none).
        num_iterations_requested: The ``num_iterations`` argument.
        steps: One :class:`BlastStep` per iteration actually run.
        interior_bridges_by_iteration: The surviving interior frontier after each
            generation. Index 0 is the initial frontier (before any iteration); index
            ``d`` is the frontier that survived ``d`` iterations. A non-empty entry at
            index ``d`` means interior bridges lasted ``d`` generations.
        completed_iterations: Number of steps actually run (≤ requested; fewer if the
            frontier emptied first or the point budget stopped it).
        stop_reason: Why the blast stopped: ``"max_iterations"`` (it ran them
            all), ``"frontier_empty"`` (no interior bridge was left to iterate) or
            ``"point_budget"`` (the next iteration was predicted to exceed
            ``max_points``).
        points_created: Points held by the bridges this blast added to the
            workbench (images that re-trace existing curve add none).
        too_close: Children dropped by ``min_separation``, over every step.
        skipped_dropped: Interior bridges left out of the initial frontier
            because an earlier call dropped them at a separation this call
            does not undercut (``Bridge.dropped_at``).
        already_known: Children that resolved to a bridge already enqueued in
            this blast (kept as the one copy, not re-iterated), over every step.
        pip_depth: ``m`` per selected fixed point with a chosen pip, read after
            the final ``infer_iterates``: the length of the pip's registered
            ``+1`` chain, ``f^m(q0)`` being its deepest registered iterate and
            so the depth ``T(q0)`` is cut at (data, never a stop condition).
    """

    zone: "ResonanceZone"
    fixed_points: list["FixedPoint"]
    num_iterations_requested: int
    steps: list[BlastStep] = field(default_factory=list)
    interior_bridges_by_iteration: list[list["Bridge"]] = field(default_factory=list)
    completed_iterations: int = 0
    stop_reason: Literal["max_iterations", "frontier_empty", "point_budget"] = (
        "max_iterations"
    )
    points_created: int = 0
    too_close: int = 0
    skipped_dropped: int = 0
    already_known: int = 0
    pip_depth: dict["FixedPoint", int] = field(default_factory=dict)

    def all_interior_bridges(self) -> list["Bridge"]:
        """Every distinct bridge that was ever an interior survivor, in order."""
        seen: dict[int, "Bridge"] = {}
        for frontier in self.interior_bridges_by_iteration:
            for b in frontier:
                seen.setdefault(id(b), b)
        return list(seen.values())


def blast_zone(
    session: "TangleSession",
    zone: "ResonanceZone | FixedPoint",
    num_iterations: int,
    *,
    fixed_point: "Optional[FixedPoint | list[FixedPoint]]" = None,
    min_separation: Optional[float] = None,
    diagnose: bool = True,
    max_points: Optional[int] = None,
    predict: bool = True,
) -> BlastResult:
    """
    Repeatedly iterate the bridges inside a resonance zone (see the module docstring).

    Args:
        session: The session owning the zone and its workbench.
        zone: A :class:`ResonanceZone`, or its key (the fixed point) to look up
            in ``session.resonance_zones``.
        num_iterations: Maximum number of blast steps.
        fixed_point: If given, only bridges emanating from this fixed point (or, if a
            list/iterable of fixed points is passed, any of them) participate;
            otherwise (``None``) every fixed point's interior bridges do.
        min_separation: A cost guard: if given, a freshly produced interior
            bridge is dropped from the frontier when its interior runs within
            ``min_separation`` of curve already in the tangle or of an
            already-kept sibling. Unstable curves can never truly cross, so a
            child this close re-traces existing curve and is about to merge into
            a numerical artifact; stopping there keeps the blast on the
            well-resolved side of the precision limit. A child ending on the
            orbit of a selected fixed point's chosen pip is exempt and logged
            at INFO, so the guard never decides the depth
            of ``T(q0)``: that is decided by whether ``f^{m+1}(q0)`` is detected
            and resolved (design note §A.4). A drop is remembered on the bridge
            (``Bridge.dropped_at``) and the zone's later calls skip it, unless
            they pass a smaller separation or None. The curve tested against is
            every curve the workbench holds (dropped and exterior children
            included), so ``n`` single-iteration blasts equal one
            ``n``-iteration blast. ``None`` (default) disables the guard.
        diagnose: If True (default) finish with ONE
            :meth:`~tanglepack.numerics.TangleWorkbench.TangleWorkbench.report_resolution`
            pass over the blasted fixed points: deep blasting is where
            under-resolved crossings and same-stability polyline crossings
            appear, and ``iterate_bridge`` never passes through
            ``compute_intersections``. Its same-stability part only COUNTS the
            hits of the bridges this blast added. It logs on the workbench's
            logger, never on this module's.
        max_points: Opt-in point budget of the whole blast (points its new
            bridges hold). Checked BETWEEN iterations: the blast stops, without
            an exception and on a consistent, resumable state, before an
            iteration that would take it over budget (``stop_reason ==
            "point_budget"``, logged at INFO). ``None`` (default) never stops.
        predict: If True (default) an iteration is refused when the points
            created so far plus its PREDICTED cost (frontier points times the
            last growth factor of this call, else the one stored on the session
            for this zone) exceed ``max_points``; with no factor yet (a zone's
            first iteration) it runs, backed by the machine's per-call cap. If
            False the blast only stops once it has already reached
            ``max_points``, so it can overshoot by one iteration.

    Note:
        Bridges that re-trace curve an existing bridge already holds -- e.g. the
        self-image produced when the fixed-point-anchored bridge is mapped forward --
        are dropped inside :meth:`TangleWorkbench.iterate_bridge` (a bridge is
        identified by its two intersections), so they never reach this loop.

    Returns:
        A :class:`BlastResult` genealogy.

    Raises:
        KeyError: If ``zone`` is a key not present in ``session.resonance_zones``.
        ValueError: If the resolved zone has no captured boundary (cannot test
            containment).
        PointBudgetExceeded: If one bridge's image would exceed the machine's
            per-call point budget. It always propagates, its message extended
            with the zone, iteration and bridge (also on ``exc.blast_context``
            as ``(zone.key, iteration, bridge.id)``), and leaves the session NOT
            resumable (that bridge's points are linked to a half-refined image);
            rebuild with fewer blasts or a larger budget. ``max_points`` is the
            budget that stops cleanly instead.
    """
    if not isinstance(zone, ResonanceZone):
        zone = session.resonance_zones[zone]
    if zone.boundary_vertices is None:
        raise ValueError(
            "Resonance zone has no captured boundary; define it via "
            "session.resonance_zone()/add_resonance_zones() before blasting."
        )

    # A bridge belongs to its INNERMOST zone. The blasted zone need not be
    # stored (define_resonance_zone, a replaced handle), so it is compared with
    # the stored zones of the OTHER fixed points rather than looked up.
    inner_zones = [
        other for other in session.resonance_zones.values()
        if other.key is not zone.key and other.area < zone.area
    ]

    def in_zone(bridge: "Bridge") -> bool:
        point = session._bridge_test_point(bridge)
        return (
            point is not None
            and zone.contains_point(point)
            and not any(inner.contains_point(point) for inner in inner_zones)
        )

    fixed_points = session._resolve_fixed_points(fixed_point)

    def fp_ok(bridge: "Bridge") -> bool:
        return any(bridge.fixed_point is fp for fp in fixed_points)

    workbench = session.workbench
    in_reach = [b for b in workbench.uniiterated_bridges if fp_ok(b) and in_zone(b)]
    # A drop stands until a call passes a smaller separation or None.
    frontier = [
        b for b in in_reach
        if b.dropped_at is None or min_separation is None or min_separation < b.dropped_at
    ]

    result = BlastResult(
        zone=zone,
        fixed_points=fixed_points,
        num_iterations_requested=num_iterations,
    )
    result.skipped_dropped = len(in_reach) - len(frontier)
    result.interior_bridges_by_iteration.append(list(frontier))
    logger.info(
        "blast: %d interior bridges in the initial frontier, %d skipped as dropped",
        len(frontier), result.skipped_dropped,
    )

    # Cycle guard: bridges already enqueued in this blast. A blast FOLLOWS the
    # dynamics -- when an interior bridge maps onto another bridge (existing or new)
    # that bridge is the next frontier and must be iterated in turn, which is the
    # whole point. We only refuse to enqueue a bridge we have already enqueued in
    # THIS blast, so the forward orbit cannot loop forever. (The single-copy
    # invariant lives in iterate_bridge, not here.) A skipped drop counts as
    # enqueued, so a child re-tracing it is already known, as within one call.
    seen: set["BridgeId | int"] = {_seen_key(b) for b in in_reach}

    # Accumulated point cloud of every bridge interior already in the tangle, used by
    # the proximity guard. A blasted bridge folds back through the same turnarounds at
    # exponentially growing cdist; once a new fold piles onto curve that already
    # exists (within ``min_separation``) the points there are iterates of iterates --
    # the floating-point error has grown until the near-coincident folds' polylines
    # cross (the "zig-zag"). Stopping a bridge whose image lands that close to ANY
    # existing fold (not just a same-generation sibling) keeps the blast on the
    # well-resolved side of that limit while still letting it extend into fresh space.
    reference_points: list[np.ndarray] = []
    if min_separation is not None:
        for b in workbench.bridges:
            interior = _interior_points(b)
            if len(interior):
                reference_points.append(interior)

    # The point budget's bookkeeping: the bridges the workbench already holds
    # (an image piece resolving to one of them adds no points) and the growth
    # factor of the zone's last blast iteration, this call's or a previous one's.
    held: set["BridgeId | int"] = {_seen_key(b) for b in workbench.bridges}
    added: list["Bridge"] = []
    factor: Optional[float] = session._blast_stats.get(zone.key)

    for iteration in range(1, num_iterations + 1):
        if not frontier:
            result.stop_reason = "frontier_empty"
            break

        frontier_points = sum(
            len(b.get_point_array()) for b in frontier if not b.iterated
        )
        if max_points is not None:
            predicted = factor * frontier_points if predict and factor else 0.0
            if (
                result.points_created + predicted > max_points
                if predict
                else result.points_created >= max_points
            ):
                logger.info(
                    "blast: point budget stops before step %d — %d points created, "
                    "%.0f predicted for %d frontier points, budget %d",
                    iteration,
                    result.points_created,
                    predicted,
                    frontier_points,
                    max_points,
                )
                result.stop_reason = "point_budget"
                break

        produced: list["Bridge"] = []
        kept: list["Bridge"] = []
        discarded: list["Bridge"] = []
        too_close: list["Bridge"] = []
        already_known: list["Bridge"] = []

        # Snapshot the accumulated manifold once per generation (one KDTree build,
        # not one per child). Bridges kept within this generation are still caught by
        # the same-generation ``kept_points`` cloud below.
        ref_tree = (
            cKDTree(np.vstack(reference_points))
            if (min_separation is not None and reference_points)
            else None
        )
        kept_points: list[np.ndarray] = []

        for parent in frontier:
            if parent.iterated:
                continue
            try:
                children = workbench.iterate_bridge(parent)
            except PointBudgetExceeded as exc:
                exc.blast_context = (zone.key, iteration, parent.id)
                exc.args = (
                    f"{exc.args[0]} [blasting zone {zone.key}, step {iteration}, "
                    f"bridge {parent.id}]",
                )
                raise

            produced.extend(children)

        # Iterating the bridge ending at f^m(q0) linked f^{m+1}(q0), so the
        # orbit read now marks every child ending on it.
        orbit = {
            x: (fp, n)
            for fp, chain in _pip_orbits(session, fixed_points).items()
            for n, x in enumerate(chain)
        }
        for child in produced:
            # A bridge already in the trellis (the fixed-point bridge's self-image
            # and any image piece that re-traces grown curve resolve to the
            # existing persistent copy): keep the single copy, do not re-iterate.
            if _seen_key(child) in seen:
                already_known.append(child)
                continue

            if not (not child.iterated and fp_ok(child) and in_zone(child)):
                seen.add(_seen_key(child))
                discarded.append(child)
                continue

            pip_end = next((x for x in child.id or () if x in orbit), None)
            if min_separation is not None and pip_end is not None:
                owner, n = orbit[pip_end]
                logger.info(
                    "blast: child %s ends at f^%d(q0) of %s (crossing %d), "
                    "kept regardless of separation",
                    child.id, n, owner.label, pip_end,
                )
            elif min_separation is not None:
                interior = _interior_points(child)
                too_near_existing = (
                    ref_tree is not None
                    and len(interior)
                    and float(ref_tree.query(interior)[0].min()) < min_separation
                )
                too_near_sibling = (
                    bool(kept_points)
                    and _min_interior_distance(child, np.vstack(kept_points))
                    < min_separation
                )
                if too_near_existing or too_near_sibling:
                    seen.add(_seen_key(child))
                    child.dropped_at = min_separation
                    too_close.append(child)
                    continue

            seen.add(_seen_key(child))
            kept.append(child)
            if min_separation is not None:
                interior = _interior_points(child)
                if len(interior):
                    reference_points.append(interior)
                    kept_points.append(interior)

        result.steps.append(
            BlastStep(
                iteration=iteration,
                interior_parents=list(frontier),
                produced_children=produced,
                kept_interior=kept,
                discarded_exterior=discarded,
                discarded_too_close=too_close,
                already_known=already_known,
            )
        )
        result.too_close += len(too_close)
        result.already_known += len(already_known)
        # Every child is curve the workbench now holds, dropped or exterior
        # alike: the next step's cloud is the one a later call would build.
        if min_separation is not None:
            reference_points += [
                p for c in discarded + too_close if len(p := _interior_points(c))
            ]

        new_points = 0
        for child in produced:
            if _seen_key(child) not in held:
                held.add(_seen_key(child))
                added.append(child)
                new_points += len(child.get_point_array())
        result.points_created += new_points
        if frontier_points > 0:
            factor = new_points / frontier_points
            session._blast_stats[zone.key] = factor

        result.completed_iterations = iteration
        result.interior_bridges_by_iteration.append(list(kept))
        logger.info(
            "blast: step %d — %d produced, %d kept interior, %d discarded, "
            "%d already known, %d too close",
            iteration,
            len(produced),
            len(kept),
            len(discarded),
            len(already_known),
            len(too_close),
        )
        frontier = kept

    # The table is exact after every blast: a crossing born early in
    # the blast may only now have its image registered.
    workbench.infer_iterates()
    result.pip_depth = {
        fp: len(orbit) - 1 for fp, orbit in _pip_orbits(session, fixed_points).items()
    }
    if diagnose:
        workbench.report_resolution(fixed_points, new_bridges=added)
    return result


