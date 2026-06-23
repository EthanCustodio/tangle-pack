"""
Blasting — driving the dynamics inside a resonance zone.

"Blasting" a resonance zone repeatedly iterates the bridges that live *inside* it. One
blast step takes every un-iterated bridge whose midpoint lies in the zone, maps it
forward one step under the dynamical map (which automatically cuts the image into child
bridges at its new crossings with the stable manifold), then keeps only the children
whose midpoint is *still* inside the zone. Those survivors are the next step's input, so
the process recurses up to ``num_iterations`` times.

Bridges that leave the zone are deliberately dropped from the frontier and never
iterated again: an unstable manifold grows by the per-step eigenvalue factor each
iterate, so iterating an *exterior* bridge many times would blow up exponentially. The
zone interior is the only region where iterating many times stays bounded, and is the
intended use of this feature.

The orchestration is intentionally tolerant: a bridge whose forward map fails (a rare
under-resolved iterate) is skipped and counted rather than aborting the whole blast, so
a long run always returns a :class:`BlastResult`. Pass ``strict=True`` to surface such
failures instead.

This is a loom-layer algorithm: it reads the topological notion of a zone and acts on
the numerical layer (:meth:`TangleWorkbench.iterate_bridge`). It is exposed on the
facade as :meth:`tanglepack.loom.TangleSession.blast_zone`.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Optional

import numpy as np

if TYPE_CHECKING:
    from ..numerics.Bridge import Bridge
    from ..numerics.FixedPoint import FixedPoint
    from .ResonanceZone import ResonanceZone
    from .TangleSession import TangleSession

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


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
    try:
        from scipy.spatial import cKDTree

        tree = cKDTree(reference)
        dists, _ = tree.query(pts)
        return float(np.min(dists))
    except Exception:  # pragma: no cover - scipy always present, brute-force fallback
        diffs = pts[:, None, :] - reference[None, :, :]
        return float(np.sqrt((diffs**2).sum(-1)).min())


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
        skipped: Parents whose forward map raised and were skipped (empty unless an
            under-resolved bridge was hit; see the module docstring).
    """

    iteration: int
    interior_parents: list["Bridge"]
    produced_children: list["Bridge"]
    kept_interior: list["Bridge"]
    discarded_exterior: list["Bridge"]
    skipped: list["Bridge"] = field(default_factory=list)
    discarded_too_close: list["Bridge"] = field(default_factory=list)
    already_known: list["Bridge"] = field(default_factory=list)


@dataclass
class BlastResult:
    """
    The genealogy produced by :func:`blast_zone`.

    Attributes:
        zone: The resonance zone that was blasted.
        fixed_point: The fixed-point filter applied (``None`` = all fixed points).
        num_iterations_requested: The ``num_iterations`` argument.
        steps: One :class:`BlastStep` per iteration actually run.
        interior_bridges_by_iteration: The surviving interior frontier after each
            generation. Index 0 is the initial frontier (before any iteration); index
            ``d`` is the frontier that survived ``d`` iterations. A non-empty entry at
            index ``d`` means interior bridges lasted ``d`` generations.
        completed_iterations: Number of steps actually run (≤ requested; fewer if the
            frontier emptied first).
        terminated_early: True if the blast stopped because no interior bridges
            remained before reaching ``num_iterations``.
        skipped: Total number of bridges skipped across all steps (forward-map
            failures).
    """

    zone: "ResonanceZone"
    fixed_point: Optional["FixedPoint"]
    num_iterations_requested: int
    steps: list[BlastStep] = field(default_factory=list)
    interior_bridges_by_iteration: list[list["Bridge"]] = field(default_factory=list)
    completed_iterations: int = 0
    terminated_early: bool = False
    skipped: int = 0
    too_close: int = 0
    already_known: int = 0

    def max_depth_reached(self) -> int:
        """The deepest generation that still held at least one interior bridge."""
        depth = 0
        for d, frontier in enumerate(self.interior_bridges_by_iteration):
            if frontier:
                depth = d
        return depth

    def all_interior_bridges(self) -> list["Bridge"]:
        """Every distinct bridge that was ever an interior survivor, in order."""
        seen: dict[int, "Bridge"] = {}
        for frontier in self.interior_bridges_by_iteration:
            for b in frontier:
                seen.setdefault(id(b), b)
        return list(seen.values())


def blast_zone(
    session: "TangleSession",
    zone: "ResonanceZone | tuple",
    num_iterations: int,
    *,
    fixed_point: "Optional[FixedPoint | list[FixedPoint]]" = None,
    strict: bool = False,
    min_separation: Optional[float] = None,
) -> BlastResult:
    """
    Repeatedly iterate the bridges inside a resonance zone (see the module docstring).

    Args:
        session: The session owning the zone and its workbench.
        zone: A :class:`ResonanceZone`, or its ``(fixed_point, branch_index)`` key to
            look up in ``session.resonance_zones``.
        num_iterations: Maximum number of blast steps.
        fixed_point: If given, only bridges emanating from this fixed point (or, if a
            list/iterable of fixed points is passed, any of them) participate;
            otherwise (``None``) every fixed point's interior bridges do.
        strict: If True, re-raise a bridge's forward-map failure instead of skipping it.
        min_separation: If given, a freshly produced interior bridge is dropped from
            the frontier when its interior runs within ``min_separation`` of an
            already-kept sibling. Unstable curves can never truly cross, so two
            bridge images that come this close are about to merge into a numerical
            artifact; stopping there keeps the blast on the well-resolved side of the
            precision limit. ``None`` (default) disables the guard.

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
    """
    zone = _resolve_zone(session, zone)
    if zone.boundary_vertices is None:
        raise ValueError(
            "Resonance zone has no captured boundary; define it via "
            "session.resonance_zone()/add_resonance_zones() before blasting."
        )

    def in_zone(bridge: "Bridge") -> bool:
        point = session._bridge_test_point(bridge)
        return point is not None and zone.contains_point(point)

    # Normalize the fixed-point filter to a set of allowed fixed points (or None for
    # "all"). Accepts a single FixedPoint or any iterable of them, so a caller can
    # blast bridges emanating from several fixed points in one pass.
    if fixed_point is None:
        allowed_fps = None
    elif isinstance(fixed_point, (list, tuple, set, frozenset)):
        allowed_fps = list(fixed_point)
    else:
        allowed_fps = [fixed_point]

    def fp_ok(bridge: "Bridge") -> bool:
        return allowed_fps is None or any(
            bridge.fixed_point is fp for fp in allowed_fps
        )

    workbench = session.workbench
    frontier = [
        b for b in workbench.uniiterated_bridges if fp_ok(b) and in_zone(b)
    ]

    result = BlastResult(
        zone=zone,
        fixed_point=fixed_point,
        num_iterations_requested=num_iterations,
    )
    result.interior_bridges_by_iteration.append(list(frontier))
    logger.info("blast: %d interior bridges in the initial frontier", len(frontier))

    # Cycle guard: bridges already enqueued in this blast. A blast FOLLOWS the
    # dynamics -- when an interior bridge maps onto another bridge (existing or new)
    # that bridge is the next frontier and must be iterated in turn, which is the
    # whole point. We only refuse to enqueue a bridge we have already enqueued in
    # THIS blast, so the forward orbit cannot loop forever. (The single-copy
    # invariant lives in iterate_bridge, not here.)
    seen: set[int] = {id(b) for b in frontier}

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

    for iteration in range(1, num_iterations + 1):
        if not frontier:
            result.terminated_early = True
            break

        produced: list["Bridge"] = []
        kept: list["Bridge"] = []
        discarded: list["Bridge"] = []
        skipped: list["Bridge"] = []
        too_close: list["Bridge"] = []
        already_known: list["Bridge"] = []

        # Snapshot the accumulated manifold once per generation (one KDTree build,
        # not one per child). Bridges kept within this generation are still caught by
        # the same-generation ``kept_points`` cloud below.
        from scipy.spatial import cKDTree

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
            except (AssertionError, AttributeError, ValueError) as exc:
                if strict:
                    raise
                logger.warning(
                    "blast: skipping bridge that failed to iterate (%s: %s)",
                    type(exc).__name__,
                    exc,
                )
                skipped.append(parent)
                continue

            produced.extend(children)
            for child in children:
                # A bridge already in the trellis (the fixed-point bridge's self-image
                # and any image piece that re-traces grown curve resolve to the
                # existing persistent copy): keep the single copy, do not re-iterate.
                if id(child) in seen:
                    already_known.append(child)
                    continue

                if not (not child.iterated and fp_ok(child) and in_zone(child)):
                    seen.add(id(child))
                    discarded.append(child)
                    continue

                if min_separation is not None:
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
                        seen.add(id(child))
                        too_close.append(child)
                        continue

                seen.add(id(child))
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
                skipped=skipped,
                discarded_too_close=too_close,
                already_known=already_known,
            )
        )
        result.skipped += len(skipped)
        result.too_close += len(too_close)
        result.already_known += len(already_known)
        result.completed_iterations = iteration
        result.interior_bridges_by_iteration.append(list(kept))
        logger.info(
            "blast: step %d — %d produced, %d kept interior, %d discarded, "
            "%d already known, %d too close, %d skipped",
            iteration,
            len(produced),
            len(kept),
            len(discarded),
            len(already_known),
            len(too_close),
            len(skipped),
        )
        frontier = kept

    return result


def _resolve_zone(
    session: "TangleSession", zone: "ResonanceZone | tuple"
) -> "ResonanceZone":
    """Resolve ``zone`` from a ResonanceZone or its ``(fixed_point, branch_index)`` key."""
    from .ResonanceZone import ResonanceZone

    if isinstance(zone, ResonanceZone):
        return zone
    return session.resonance_zones[zone]
