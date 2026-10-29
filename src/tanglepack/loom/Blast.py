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

A LEAN blast (``lean=True``, 2026-10-29) frees, after each iteration, the
curve nothing reads again: it links the iterate table, recomputes the topology
of the selection and keeps the bridges the topology and later blasts can read
(:func:`_lean_keep`); every other bridge's interior goes
(:meth:`~tanglepack.numerics.TangleWorkbench.TangleWorkbench.release_curves`).
The geometric guard reads every interior, so a lean blast takes the width
guard (``min_separation="width"``) or none. The width guard needs no point
cloud: a lobe's width (its area over its unstable length) contracts by the
stable factor per step, and a child thinner than ``WIDTH_SEPARATION_FACTOR``
chord deviations of its own polyline is below the resolution.

This is a loom-layer algorithm: it reads the topological notion of a zone and acts on
the numerical layer (:meth:`TangleWorkbench.iterate_bridge`). It is exposed on the
facade as :meth:`tanglepack.loom.TangleSession.blast_zone`.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Callable, Literal, Optional

import numpy as np
from scipy.spatial import cKDTree

from ..numerics.Bridge import CurveReleased
from ..numerics.ManifoldMachine import DEFAULT_AREA_CUTOFF
from ..numerics.budget import PointBudgetExceeded
from ..numerics.geometry import arc_polyline, signed_polygon_area
from .ResonanceZone import ResonanceZone

if TYPE_CHECKING:
    from ..numerics.Bridge import Bridge, BridgeId
    from ..numerics.FixedPoint import FixedPoint
    from ..numerics.TangleWorkbench import TangleWorkbench
    from .TangleSession import TangleSession

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

#: ``min_separation="resolution"`` at the default ``area_cutoff``: the recipes'
#: separation. It scales as the polylines' chord deviation, ``area_cutoff ** (2/3)``
#: (a segment of a curve with curvature ``k`` refined to area ``A`` is
#: ``L ~ (12 A / k) ** (1/3)`` long and strays ``h = 1.5 A / L`` from its chord), so
#: a finer cutoff resolves, and keeps, curves that run closer together. The author's
#: rule, 2026-10-21: the guard follows the resolution.
RESOLUTION_SEPARATION: float = 1e-5

#: ``min_separation="width"`` drops a child whose width is below this many
#: chord deviations ``h = 1.5 * area_cutoff / L`` of its own polyline (``L`` its
#: median chord): a lobe thinner than the polygon error is not resolved.
WIDTH_SEPARATION_FACTOR: float = 1.0

#: ``Bridge.dropped_at`` of a width drop: below every float separation, so only
#: ``min_separation=None`` re-admits it.
WIDTH_DROP: float = -1.0


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


class _Cloud:
    """The interiors the proximity guard tests against, as a few KD-trees.

    Bentley and Saxe's logarithmic method: the points sit in trees of distinct
    sizes, and adding a batch merges it with every tree no larger than it, so
    there are O(log n) trees and each point is rebuilt O(log n) times. A query
    reads every tree. The verdict is that of one tree over the union -- the
    distance to a union is the least over its parts -- without building one
    tree per child over everything kept so far. The queries run in the caller:
    a child's interior is too small for ``workers`` threads to pay (2026-10-29,
    the 1e-7 scale series: separation 4.0 s at one worker, 5.8 s at 16).
    """

    def __init__(self) -> None:
        self.trees: list[cKDTree] = []

    def add(self, interiors: list[np.ndarray]) -> None:
        """Hold the points of ``interiors`` too (one batch: one merge)."""
        interiors = [points for points in interiors if len(points)]
        if not interiors:
            return
        points = np.vstack(interiors)
        while self.trees and self.trees[-1].n <= len(points):
            points = np.vstack((self.trees.pop().data, points))
        # The sliding-midpoint tree builds twice as fast; its nearest points are exact too.
        self.trees.append(cKDTree(points, balanced_tree=False, compact_nodes=False))

    def near(self, points: np.ndarray, separation: float) -> bool:
        """Whether a point of ``points`` lies closer than ``separation`` to one held.

        The queries stop at ``separation``: a nearer point is found exactly, and
        a farther one cannot change the verdict.
        """
        return len(points) > 0 and any(
            float(tree.query(points, distance_upper_bound=separation)[0].min())
            < separation
            for tree in self.trees
        )


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


def _ring_width(session: "TangleSession", bridge: "Bridge") -> Optional[float]:
    """The width of a bridge's lobe: ``|area| / L`` of the ring the bridge and
    its stable segment bound, ``L`` the unstable arc length.

    None when the two crossings lie on two stable branches (no ring), or the
    bridge is a partial piece. At an anchor an arc of ``u(o, b)`` meets
    ``s(o, b)``.
    """
    if bridge.partial:
        return None
    registry = session.workbench.intersection_registry
    first, second = (registry[i] for i in bridge.id)
    fp, _, orbit, branch = bridge.manifold_key
    keys = [
        (fp, "stable", orbit, branch) if x.is_anchor else x.manifold_b_key for x in (first, second)
    ]
    if keys[0] != keys[1]:
        return None
    stable = session.workbench.manifolds[keys[0]]
    unstable = np.vstack((
        first.coords,
        arc_polyline(bridge, first.unstable_cdist, second.unstable_cdist),
        second.coords,
    ))
    back = arc_polyline(
        stable, first.stable_cdist, second.stable_cdist,
        reverse=second.stable_cdist > first.stable_cdist,
    )
    length = float(np.hypot(*np.diff(unstable, axis=0).T).sum())
    return abs(signed_polygon_area(np.vstack((unstable, back)))) / length


def _too_thin(session: "TangleSession", child: "Bridge") -> bool:
    """Whether a child's width is below ``WIDTH_SEPARATION_FACTOR`` chord deviations."""
    if child.width is None:
        return False
    chords = np.hypot(*np.diff(child.get_point_array(), axis=0).T)
    deviation = 1.5 * session.workbench.area_cutoff / float(np.median(chords))
    return child.width < WIDTH_SEPARATION_FACTOR * deviation


def _lean_keep(
    session: "TangleSession",
    fixed_points: list["FixedPoint"],
    in_zone: Callable[["Bridge"], bool],
) -> list["BridgeId | Bridge"]:
    """What a lean blast keeps: the bridges the topology and later blasts read.

    The future frontiers (every un-iterated, undropped bridge some stored zone,
    or the blasted one, holds), the minimal trellis of the selection (and of
    every fixed point, when the selection is not all of them), the pieces of
    every joined arc a trellis holds (``Bridge.joined`` reads their rows) and
    every bridge ending on a pip-orbit crossing. The straddle rows a frontier
    shares with an iterated neighbour are kept by ``release_curves`` itself.
    """
    workbench = session.workbench
    keep: list["BridgeId | Bridge"] = [
        b for b in workbench.uniiterated_bridges
        if b.dropped_at is None and (session.classify_bridge(b) is not None or in_zone(b))
    ]
    selections = [fixed_points]
    if set(map(id, fixed_points)) != set(map(id, workbench.fixed_points)):
        selections.append(None)
    trellises = []
    for selection in selections:
        topology = session.topology(selection)
        keep += topology.minimal.kept_bridge_ids
        trellises.append(topology.trellis)
    # Every cached trellis memoises its joined arcs (a re-cut hands them back).
    for trellis in {id(t): t for t in [*trellises, *session._trellises.values()]}.values():
        for joined in trellis._joined.values():
            keep += _pieces(workbench, joined)
    for orbit in _pip_orbits(session, fixed_points).values():
        for x in orbit:
            keep += workbench.bridges_at(x)
    return keep


def _pieces(workbench: "TangleWorkbench", joined: "Bridge") -> list["BridgeId"]:
    """The ids of the consecutive bridges a joined arc runs through."""
    pieces, at = [], joined.first_intersection
    while at != joined.second_intersection:
        bid = next(
            bid for bid in workbench.bridges_at(at)
            if bid[0] == at and workbench.bridge(bid).manifold_key == joined.manifold_key
        )
        pieces.append(bid)
        at = bid[1]
    return pieces


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
        released_rows: Store rows a lean blast freed, over every step (0 for
            a blast that is not lean).
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
    released_rows: int = 0

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
    min_separation: "Optional[float | Literal['resolution', 'width']]" = None,
    diagnose: Optional[bool] = None,
    max_points: Optional[int] = None,
    predict: bool = True,
    lean: bool = False,
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
            and resolved (design note §A.4). ``"resolution"`` ties it to the
            machine's ``area_cutoff``: ``RESOLUTION_SEPARATION * (area_cutoff /
            DEFAULT_AREA_CUTOFF) ** (2/3)``, read when the call starts. A drop
            is remembered on the bridge
            (``Bridge.dropped_at``) and the zone's later calls skip it, unless
            they pass a smaller separation or None. The curve tested against is
            every curve the workbench holds (dropped and exterior children
            included), so ``n`` single-iteration blasts equal one
            ``n``-iteration blast. ``"width"`` is the analytic guard instead
            (no point cloud): a child is dropped when its lobe width, carried
            down the image chain from the ring area of the first bridge
            iterated (``Bridge.width``), is below ``WIDTH_SEPARATION_FACTOR``
            chord deviations of its polyline; the pip-orbit exemption holds,
            and a width drop (``dropped_at == WIDTH_DROP``) is re-admitted by
            ``None`` only (never after a lean blast released the dropped
            bridge). ``None`` (default) disables the guard.
        diagnose: If True (the default unless ``lean``) finish with ONE
            :meth:`~tanglepack.numerics.TangleWorkbench.TangleWorkbench.report_resolution`
            pass over the blasted fixed points: deep blasting is where
            under-resolved crossings and same-stability polyline crossings
            appear, and ``iterate_bridge`` never passes through
            ``compute_intersections``. Its same-stability part only COUNTS the
            hits of the bridges this blast added. It logs on the workbench's
            logger, never on this module's. A lean blast cannot diagnose: the
            same-stability reading needs every curve.
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
        lean: If True, after each iteration link the iterate table, recompute
            the topology of the selection and free the curve no reader needs
            (:func:`_lean_keep`); ``BlastResult.released_rows`` counts the rows.
            The workbench is lean from then on and refuses whatever would
            re-read freed curve (growth, recomputes, zones, save, a blast
            that is not lean), raising
            :class:`~tanglepack.numerics.Bridge.CurveReleased`.

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
            containment), or ``lean`` is passed with a geometric separation
            (a float or ``"resolution"``) or with ``diagnose=True``.
        CurveReleased: A blast that is not lean on a lean workbench, or a
            ``min_separation=None`` blast that would re-admit a dropped bridge
            a lean blast released.
        PointBudgetExceeded: If one bridge's image would exceed the machine's
            per-call point budget. It always propagates, its message extended
            with the zone, iteration and bridge (also on ``exc.blast_context``
            as ``(zone.key, iteration, bridge.id)``), and leaves the session NOT
            resumable (that bridge's points are linked to a half-refined image);
            rebuild with fewer blasts or a larger budget. ``max_points`` is the
            budget that stops cleanly instead.
    """
    if lean and min_separation not in (None, "width"):
        raise ValueError(
            f"a lean blast frees the interiors min_separation={min_separation!r} reads; "
            "pass 'width' or None"
        )
    if lean and diagnose:
        raise ValueError("a lean blast cannot diagnose: the same-stability reading needs every curve")
    if not lean:
        session.workbench._refuse_if_lean("a blast that is not lean")
    diagnose = not lean if diagnose is None else diagnose
    width = min_separation == "width"
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
        point = bridge.test_point
        return (
            point is not None
            and zone.contains_point(point)
            and not any(inner.contains_point(point) for inner in inner_zones)
        )

    fixed_points = session._resolve_fixed_points(fixed_point)

    def fp_ok(bridge: "Bridge") -> bool:
        return any(bridge.fixed_point is fp for fp in fixed_points)

    workbench = session.workbench
    if min_separation == "resolution":
        min_separation = RESOLUTION_SEPARATION * (
            workbench.area_cutoff / DEFAULT_AREA_CUTOFF
        ) ** (2 / 3)
        logger.info("blast: min_separation %.3g from the resolution", min_separation)
    in_reach = [b for b in workbench.uniiterated_bridges if fp_ok(b) and in_zone(b)]
    # A drop stands until a call passes a smaller separation or None.
    frontier = [
        b for b in in_reach
        if b.dropped_at is None
        or min_separation is None
        or (not width and min_separation < b.dropped_at)
    ]
    # A lean blast released the bridges it dropped (_lean_keep): re-admitting
    # one would read freed curve mid-blast.
    released = [b.id for b in frontier if b.released]
    if released:
        raise CurveReleased(
            f"bridges {released} were dropped and then released by a lean blast, so "
            "min_separation=None cannot re-admit them; rebuild the session and blast "
            "it without lean"
        )

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
    cloud = _Cloud()
    if min_separation is not None and not width:
        cloud.add([_interior_points(b) for b in workbench.bridges])

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

        frontier_points = sum(b.point_count for b in frontier if not b.iterated)
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

        for parent in frontier:
            if parent.iterated:
                continue
            if width and parent.width is None:
                parent.width = _ring_width(session, parent)
            try:
                children = workbench.iterate_bridge(parent)
            except PointBudgetExceeded as exc:
                exc.blast_context = (zone.key, iteration, parent.id)
                exc.args = (
                    f"{exc.args[0]} [blasting zone {zone.key}, step {iteration}, "
                    f"bridge {parent.id}]",
                )
                raise

            if width and parent.width is not None:
                contraction = parent.fixed_point.per_step_beta("stable")
                for child in children:
                    if child.width is None:
                        child.width = parent.width * contraction
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
            elif width and _too_thin(session, child):
                seen.add(_seen_key(child))
                child.dropped_at = WIDTH_DROP
                too_close.append(child)
                continue
            elif (
                min_separation is not None
                and not width
                and cloud.near(_interior_points(child), min_separation)
            ):
                seen.add(_seen_key(child))
                child.dropped_at = min_separation
                too_close.append(child)
                continue

            seen.add(_seen_key(child))
            kept.append(child)
            # A kept child is curve every later child keeps clear of.
            if min_separation is not None and not width:
                cloud.add([_interior_points(child)])

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
        if min_separation is not None and not width:
            cloud.add([_interior_points(c) for c in discarded + too_close])

        new_points = 0
        for child in produced:
            if _seen_key(child) not in held:
                held.add(_seen_key(child))
                added.append(child)
                new_points += child.point_count
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
        if lean:
            # The table and the topology the keep rule reads are this
            # iteration's: a later one never needs what they do not name.
            workbench.infer_iterates()
            result.released_rows += workbench.release_curves(
                _lean_keep(session, fixed_points, in_zone)
            )

    # The table is exact after every blast: a crossing born early in
    # the blast may only now have its image registered.
    workbench.infer_iterates()
    result.pip_depth = {
        fp: len(orbit) - 1 for fp, orbit in _pip_orbits(session, fixed_points).items()
    }
    if diagnose:
        workbench.report_resolution(fixed_points, new_bridges=added)
    return result


