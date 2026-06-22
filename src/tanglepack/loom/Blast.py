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

if TYPE_CHECKING:
    from ..numerics.Bridge import Bridge
    from ..numerics.FixedPoint import FixedPoint
    from .ResonanceZone import ResonanceZone
    from .TangleSession import TangleSession

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


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
    fixed_point: Optional["FixedPoint"] = None,
    strict: bool = False,
) -> BlastResult:
    """
    Repeatedly iterate the bridges inside a resonance zone (see the module docstring).

    Args:
        session: The session owning the zone and its workbench.
        zone: A :class:`ResonanceZone`, or its ``(fixed_point, branch_index)`` key to
            look up in ``session.resonance_zones``.
        num_iterations: Maximum number of blast steps.
        fixed_point: If given, only bridges emanating from this fixed point
            participate; otherwise every fixed point's interior bridges do.
        strict: If True, re-raise a bridge's forward-map failure instead of skipping it.

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

    def fp_ok(bridge: "Bridge") -> bool:
        return fixed_point is None or bridge.fixed_point is fixed_point

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

    for iteration in range(1, num_iterations + 1):
        if not frontier:
            result.terminated_early = True
            break

        produced: list["Bridge"] = []
        kept: list["Bridge"] = []
        discarded: list["Bridge"] = []
        skipped: list["Bridge"] = []

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
                if not child.iterated and fp_ok(child) and in_zone(child):
                    kept.append(child)
                else:
                    discarded.append(child)

        result.steps.append(
            BlastStep(
                iteration=iteration,
                interior_parents=list(frontier),
                produced_children=produced,
                kept_interior=kept,
                discarded_exterior=discarded,
                skipped=skipped,
            )
        )
        result.skipped += len(skipped)
        result.completed_iterations = iteration
        result.interior_bridges_by_iteration.append(list(kept))
        logger.info(
            "blast: step %d — %d produced, %d kept interior, %d discarded, %d skipped",
            iteration,
            len(produced),
            len(kept),
            len(discarded),
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
