"""
The stable-manifold partition: holes, intervals and their propagation.

Places a hole in each region bounded by a pseudoneighbor pair and its bridge,
propagates the holes around the orbit, and cuts each stable branch into the
partition intervals the region layer reads.

Dev Notes — Stable Manifold Partition (Stable_Manifold_Partition_Algorithm.pdf)

The PDF is less formal than the Pseudoneighbor one; the following semantics were
fixed in conversation with the author (July 2026):

* Hole placement: a pseudoneighbor pair and its connecting unstable arc (a
  bridge) bound a closed region; the hole sits inside that region. A REFERENCE
  pair's hole is plotted near the stable manifold, an iterated/propagated hole
  between the two neighbors on the unstable manifold — per the author. Both
  step between the middles of the two BOUNDARY ARCS (the actual stable-arc
  node between the pair, walked via trellis.manifolds, and the bridge's middle
  node): in a narrow curved lobe the pair's chord midpoint can lie OUTSIDE the
  region, so chord-based placement mislocated holes (author-reported on the
  period-3 tangle). Every hole
  carries its orbit identity: ``origin`` (the reference pair it descends
  from, which selects the plot marker) and ``iterate`` (0 for the reference,
  negative backward, positive forward — the plot label).
  ``near_intersection_id`` is fixed by convention to the toward-anchor member
  (smaller stable cdist) — the choice is arbitrary but must be consistent.
* Partition elements (Phase 5): every interval is stamped with
  ``(branch_key, side, element_id)`` and two lookup tables are built beside it,
  ``element_of_intersection`` (crossing -> owning element, by id at a
  boundary and by EXACT cdist elsewhere, see "Exact ownership" below) and
  ``elements_at_bridge``. The latter is only a
  convenience index over the bridges the building trellis HELD: a trellis keeps a
  bridge under its unstable fixed point, so a heteroclinic bridge never appears
  in the trellis that owns its endpoints' stable branch. The authoritative
  lookup is therefore :meth:`Trellis.element_for`, which falls back to the
  endpoint ids carried by the BridgeId itself.
* Side (author, 2026-07-16): a hole's side is left/right of its own BRIDGE,
  i.e. of the unstable manifold, in the bridge's dynamical orientation —
  forward flow points AWAY from the fixed point, so the bridge is oriented by
  increasing unstable cdist; standing on it looking that way, positive
  cross(tangent, displacement) = left (``_side_of``). It is decided ONCE and
  stored as ``Hole.bridge_side``: a direct hole reads it off its crossing
  signs (2026-10-02), a propagated hole CARRIES it from the hole it was
  propagated from, parity-flipped under an orientation-reversing map
  (author, 2026-10-07; ``Hole.side_source``). Holes know nothing of the
  stable manifold's sides or of resonance zones. The stable partition keeps
  its own left/right, the ROW, defined by the STABLE dynamical direction
  (looking toward the anchor).
* Signs decide (D6, 2026-10-08; chalkboard stage 6). At a crossing ``x``
  with ``sign = crossing_sign(x) = sign(cross(u+, s+))`` the outward stable
  ray ``s+`` lies on the left of the bridge's ``u+`` exactly when the sign
  is positive, so: ``row(first) = left iff sign > 0``, ``row(second) = left
  iff sign < 0`` (:func:`row_of_end`, the bridge leaves ``first`` along
  ``u+`` and enters ``second`` along ``u-``), and the stable half on a
  hole's side of its bridge is ``outward iff (sign > 0) == (side ==
  "left")`` (:func:`_hole_openings`). No polyline decides anything. The
  geometry stays as a cross-check that FLAGS under-resolution
  (:func:`geometry_flags`, one WARNING per kind per ``Trellis.punch_holes``):
  :func:`bridge_row_diagnosis` reads each end's approach on the polyline,
  which on a 2-4 node bridge is the vector to the far node (NUM-1, the E11
  p3 / nested_inner 1e-5 cells from blast 6), and ``Hole.measured_side``
  reads the backward-carried reference point at the nearest segment of the
  image sub-arc, which escapes the lobe when the reference lobe is
  under-resolved (NUM-2, k28 1e-5/1e-6 blasts 15-17; inversion x period 2 at
  10 steps). Before WP4 the openings and the I2 assertion read the polyline
  (a nearest vertex of the WHOLE bridge for the halves, INT-4), which opened
  nothing, or the same half at both ends, at depth: the unreachable class of
  the deep p3 runs (PK-1) was those openings.
* Hole generation (chalkboard stage 5): each reference bridge is mapped
  backward step by step; every backward image lies within some existing
  bridge, and a hole is punched in that containing bridge's region.
  Termination is the author's "the bridge begins to map onto itself" rule,
  stated k-aware: one backward step moves the image to the PREVIOUS unstable
  branch of the cycle and a bridge lives on one branch, so meeting a visited
  bridge again means it maps onto itself after a multiple of k steps (for
  k = 1 the consecutive self-map condition). The image span is the start
  bridge's two crossing ids moved by :meth:`Trellis.image_cdist` (table
  first, else ``advance_key`` and beta-scaling; nothing else scales a cdist
  here), on the branch ``advance_key`` names, and the container is a bisect
  over that branch's crossing order (:func:`_containing_bridge`; PERF-4).
  The hole sits on the WHOLE containing bridge, by design (INT-2, checked
  2026-10-08 on every law case, the deep p3 runs and four deep E11 states:
  600+ lookups, every span a STRICT sub-arc of its container and at most one
  of its two ends a registered crossing; the other end is a crossing of the
  unstable branch with ``f^-1(W^s)``, beyond the trimmed stable manifold, so
  the trellis has no node to cut the sub-arc at). The reference hole's
  coordinates are carried backward with the real inverse map (the one use
  of the map here): the carried point places the hole and is the measured
  side's flag — two orbits landing in the same bridge differ exactly there,
  and BOTH are punched (per the author: no cross-orbit deduplication, even
  on the same bridge and side; only a same-orbit same-iterate duplicate of a
  recorded pair's hole is skipped).
* Superseded: resonance-zone association (removed 2026-07-16). Holes used
  to carry an ``interior`` flag (in_zone point-in-polygon test, else a
  cdist-vs-strong-pip-cut fallback) and the partition split its intervals at
  the zone cut for interior/exterior reporting. Both are gone: which region
  a hole belongs to is a dynamical-direction question (its side of the
  bridge), not a zone-membership one — the descriptive flag invited
  zone-gated filtering, which the author had already rejected ("all holes on
  manifolds should be treated equal"), and the deep backward sub-arcs hug
  the zone's own boundary bridge, so the point-in-polygon test split
  otherwise identical holes by the luck of a nudge direction. Resonance
  zones remain in the loom layer for blasting only.
* Partition intervals (semantics fixed with the author, July 2026 — this
  SUPERSEDES the earlier "only the outermost beyond-cut hole" rule and the
  piece-containment openness test): EVERY punched hole participates.
  Boundaries are the defining intersections of the holes' bridges
  (Hole.bounding_ids — the pair's own bridge for a direct hole, the
  containing bridge for a propagated one) plus the anchor and the branch's
  outermost intersection. The openness rule is LOCAL to each defining
  intersection: the arc crosses the stable manifold there, so each side of
  the bridge holds exactly one of the two stable half-intervals, and the
  hole opens the one on ITS side (Hole.openings, by the signs). A direct
  hole, facing its pair's own stable segment, thereby opens the inward pair
  (outward of the near bound, anchorward of the far bound). A hole whose
  region hangs off the removed tail opens the OUTWARD side of the outermost
  intersection — two holes flanking one bridge (the author's "singleton
  bridge with a hole on either side") open complementary sides of BOTH
  defining intersections and pinch each into a closed singleton. In general
  a boundary point excluded by both neighbouring intervals owns itself as a
  degenerate closed singleton [x, x]; a hole marked outward of the branch's
  outermost intersection or anchorward at the anchor contributes its
  opening to that pinch test even though no interval lies there. No point
  is special — the fixed point and the strong pip pinch into singletons by
  exactly the same rule. Each point belongs to exactly one piece ("](" and
  ")[" occur, "][" cannot).
* Invariants (2026-09, plan row 0.3). Two facts are physically exact and
  are asserted in ``Trellis.punch_holes``, once per orbit and once per hole
  bridge: (I1) all holes of one ``origin`` share ``bridge_side``, because a
  propagated hole is the backward image of the reference hole and an
  orientation-PRESERVING map (det J > 0, both Hénon fixtures have b = 1)
  carries the bridge's dynamical orientation to the image bridge's
  unchanged — under an orientation-REVERSING map the side alternates with
  the parity of ``Hole.iterate``
  (``check_holes_share_bridge_side(..., orientation_preserving=False)``);
  propagated holes carry their side, so I1 guards the direct holes of one
  orbit (the reference and its ``+1 .. +(k-1)`` appearances) and the
  carrying itself; and (I2) a bridge whose two crossings lie on the same
  stable branch has the same row at both ends, because its endpoints are
  consecutive crossings so the arc between them never crosses the stable
  manifold again — with rows read off the signs, I2 says the two signs
  DIFFER (equal signs still assert, D6). Failures name every reading:
  :class:`SideDiagnosis` per hole, :class:`RowDiagnosis` per bridge (crossing
  signs, combinatorial and geometric rows, approach angle, node count).
* BACKWARD ONLY (author's rule, 2026-09-16). Once a hole lands on the stable
  manifold its forward iterates get NO hole, even when the trellis has
  registered the image crossings (a blast does exactly that) and a bridge
  spans them: the image region is still attached to the stable manifold, and
  only the (uncomputed) image of the bridge on the far side of it could
  justify a boundary there. Before this rule ``punch_holes`` punched a direct
  hole for EVERY recorded pair, and ``extend_pseudoneighbor_trajectories``
  records forward pairs too, so the k=2.8 blasted tangle acquired a hole in
  the blast-child bridge (the +1 image of a reference pair) that split the
  partition and produced a fourth bridge class; k=10 escaped only because its
  reference pair's forward images were never registered (the manifold was not
  grown that far). PERIOD-k CHOICE, PROVISIONAL: forward iterates
  ``1 .. k_value - 1`` are still punched — they are how the other branches of
  the orbit get their fundamental-segment holes, since the reference window
  lives on the strong pip's branch alone (Pseudoneighbor Dev Notes) — and only
  ``iterate >= k_value`` is dropped. The author asked for this to be recorded
  as a choice to revisit, not a settled rule (``_is_forward_beyond_fundamental``
  is the one place it lives).
* Propagation start: a reference whose pair has no spanning Bridge object
  (after a blast there is no bridge between a parent-manifold crossing and
  a blast-child crossing) punches no direct hole, but its orbit still
  propagates — the backward chain continues from the orbit's most-backward
  directly-punched pair instead (author, 2026-07-09: the un-punchable pair
  itself is fine to skip, its deeper backward holes are not).
* Exact ownership (2026-10-07, stress campaign E11). Every partition boundary
  IS a registry float (asserted in :func:`_element_of_intersection`; the
  synthetic ``(None, 0.0)`` anchor aside), so a registered crossing is owned
  without any tolerance: a boundary crossing by the element closed at its id
  (the singleton if it has one, else the closed neighbour), any other by its
  exact stable cdist. The registry's absolute ``cdist_tol`` (1e-6) is a
  DEDUPLICATION tolerance, not an ownership one: deep blasts put real
  crossings closer than it (68% of consecutive k=2.8 gaps by blast 14), and
  a sliver ``(pip, outer]`` of width 1.1e-6 beyond the zone's pip had its
  midpoint swallowed by the closed pip singleton, which tripped the
  child-inside-parent assertion of the iterated partition at 9 blasts. The
  tolerance survives only where a cdist is ESTIMATED (``land_element``'s
  image spans, ``Trellis.image_of_element``): ``PartitionFamily.interval_at``
  tries exact first and falls back to it. The anchor is detected by its exact
  cdist 0.0 (a registered anchor is synthetic at exactly 0.0); a boundary
  within the tolerance of 0 but not on it is still taken as the anchor, at
  INFO.
* The anchor (WP3, 2026-10-08). A periodic point is ONE crossing on every
  branch attached to it (D4). Both stable branches leave it, so an anchor
  opening faces the outward half of the branch the signs pick
  (:func:`_anchor_branch`) and every opening record names the branch it
  opens; the anchor has no geometric row of its own. A stable branch with
  NO crossing past its anchor (an inversion saddle grown 4-5 steps) is ONE
  element ``[anchor, end]`` to the end of the computed branch (C1), logged
  at INFO on this module's logger; grow the unstable manifold to resolve it.

Caveats: inversion (k_value == 2*period) follows the same cycle bookkeeping as
StrongPip/Pseudoneighbor and is exercised by the inversion law cases (the k=4
saddle at 0/1/2 blasts and the k=4.1 period-2 orbit). Unstable cdists are only
comparable on one unstable branch, so every bridge carries its branch
(``Bridge.manifold_key``, required at construction since plan row 2.3: an
endpoint-derived branch once handed one branch's anchor to another's first
bridge). A pseudoneighbor pair with no spanning Bridge object (after blasting
there is no bridge between a parent-manifold crossing and a blast-child
crossing) punches no hole and therefore contributes nothing to the partition —
an open gap in the punching layer, see punch_holes' warning; a backward image
landing on such a pair has no container and ends its chain. A bridge whose
endpoints have equal unstable cdists (impossible for a real bridge — its
endpoints are consecutive distinct crossings — but reachable via numerical
collapse) cannot be oriented: it has no measured side and no geometric row.
"""

from __future__ import annotations

import logging
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from typing import Iterable, Optional, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from ..numerics.geometry import oriented_bridge_polyline
from .TopologyResults import (
    Endpoint,
    endpoint_index,
    Hole,
    OPPOSITE_SIDE,
    PartitionInterval,
    PseudoneighborPair,
    Side,
    SideSource,
    StablePartitionResult,
)

if TYPE_CHECKING:
    from .Trellis import Trellis
    from .TrellisBranch import TrellisBranch
    from ..numerics.Bridge import Bridge, BridgeId
    from ..numerics.FixedPoint import FixedPoint
    from ..numerics.Intersection import Intersection, ManifoldKey

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def bridge_for_pair(trellis: "Trellis", pair: PseudoneighborPair) -> Optional["Bridge"]:
    """
    Return the bridge whose endpoints are exactly this pseudoneighbor pair.

    A pseudoneighbor pair is consecutive along its unstable branch (its open
    unstable interval is empty), so after ``create_bridges`` a bridge normally
    spans it.

    Args:
        trellis: The Trellis whose bridges to search.
        pair: The pseudoneighbor pair.

    Returns:
        The matching Bridge, or None if no bridge spans the pair.

    Note:
        The search itself is the trellis's own endpoint index
        (:meth:`Trellis.bridge_between`), so repeated lookups over one snapshot
        cost a dict probe each.
    """
    first, second = pair.as_tuple()
    return trellis.bridge_between(first, second)


def punch_holes(
    trellis: "Trellis",
    pairs: Optional[Iterable[PseudoneighborPair]] = None,
    *,
    epsilon: float = 0.05,
) -> list[Hole]:
    """
    Punch a hole in the region bounded by each pseudoneighbor pair.

    Placement inside the bounded region depends on the pair's role: a
    REFERENCE pair's hole hugs the stable manifold (the pair's chord midpoint,
    nudged ``epsilon`` toward the bridge), while an iterated pair's hole sits
    between the two neighbors on the unstable manifold (the bridge midpoint,
    nudged ``epsilon`` toward the chord). Each hole's side of the pair's
    bridge, in the bridge's dynamical orientation (``Hole.bridge_side``), and
    its openings are read off the crossing signs (D6); it carries its orbit
    identity (``origin``/``iterate``) and is attached to its pair
    (``pair.hole``).

    Holes propagate BACKWARD only: a pair recorded at a forward iterate
    ``>= k_value`` of its fixed point is the image of a region still attached
    to the stable manifold and gets no direct hole (it is logged at DEBUG and
    left with ``pair.hole = None``). Forward iterates ``1 .. k_value - 1`` are
    the reference trajectory's appearances on the other branches' fundamental
    segments and are punched like references. See the module Dev Notes.

    Args:
        trellis: The Trellis carrying the pairs and bridges.
        pairs: Pairs to punch holes for. Defaults to every pair recorded on the
            trellis (``trellis.pseudoneighbors``).
        epsilon: Inward step off the manifold, as a fraction of the pair's
            chord length.

    Returns:
        The punched holes. Pairs with no spanning bridge are skipped with a
        warning; forward-iterate pairs beyond the fundamental segment are
        skipped silently (DEBUG).
    """
    if pairs is None:
        pairs = trellis.pseudoneighbors
    pairs = list(pairs)

    # Punching is idempotent per call: a pair that gets no hole this time must
    # not keep the one an earlier call attached, because propagate_reference_holes
    # reads pair.hole to seed the already-punched orbits and regions.
    for pair in pairs:
        pair.hole = None

    # Holes propagate BACKWARD only (author's rule, 2026-09-16): a forward
    # iterate of a punched hole is the image of a region still attached to the
    # stable manifold and gets no hole of its own. The forward iterates
    # +1 .. +(k_value - 1) are exempt: they are the reference trajectory's
    # appearances on the OTHER branches' fundamental segments (the reference
    # window lives on the strong pip's branch only). See the module Dev Notes.
    forward = [pair for pair in pairs if _is_forward_beyond_fundamental(trellis, pair)]
    if forward:
        logger.debug(
            "Not punching %d forward-iterate pair(s) (iterate >= k_value): "
            "holes propagate backward only; %s",
            len(forward),
            [(pair.as_tuple(), pair.iterate) for pair in forward],
        )
        dropped = {id(pair) for pair in forward}
        pairs = [pair for pair in pairs if id(pair) not in dropped]

    holes: list[Hole] = []
    skipped = 0
    for pair in pairs:
        bridge = bridge_for_pair(trellis, pair)
        if bridge is None:
            skipped += 1
            continue

        near_id, far_id = _near_far(trellis, pair.intersection_a, pair.intersection_b)
        near = trellis.intersection(near_id)
        far = trellis.intersection(far_id)
        # Placement anchors on the boundary-arc midpoints; the hole sits
        # strictly between the pair's stable arc and the bridge arc, so its
        # coordinates classify its side of the bridge.
        stable_mid, _stable_tangent = _stable_arc_midpoint(trellis, near, far)
        arc_point = _nearest_arc_point(bridge, stable_mid)
        if arc_point is None:
            skipped += 1
            continue

        chord_len = float(np.linalg.norm(far.get_point() - near.get_point()))
        if pair.is_reference:
            coords = _step_into(stable_mid, arc_point, epsilon, chord_len)
        else:
            coords = _step_into(arc_point, stable_mid, epsilon, chord_len)
        coords = (float(coords[0]), float(coords[1]))

        # A direct hole sits in the lobe between the bridge and the pair's own
        # stable segment, so its side is read off the crossing signs like
        # `row_of_end` (2026-10-02): the bridge leaves `first` along u+, the
        # segment runs from `first` along s+ when it heads outward, and s+ is
        # on the left of u+ iff crossing_sign > 0. Measuring the nearest
        # bridge vertex picked the wrong fold on thin lobes.
        first, second = sorted((near, far), key=lambda x: x.unstable_cdist)
        outward = 1 if second.stable_cdist > first.stable_cdist else -1
        bridge_side = "left" if first.crossing_sign * outward > 0 else "right"
        openings = _hole_openings(trellis, bridge, bridge_side)

        hole = Hole(
            coords=coords,
            near_intersection_id=near_id,
            pair=pair,
            bridge_side=bridge_side,
            bounding_ids=(near_id, far_id),
            openings=openings,
            iterate=pair.iterate,
            origin=pair.origin,
            side_source="signs",
        )
        pair.hole = hole
        holes.append(hole)

    if skipped:
        logger.warning("Skipped %d pair(s) with no spanning bridge", skipped)
    logger.debug("Punched %d hole(s)", len(holes))
    return holes


def propagate_reference_holes(
    trellis: "Trellis",
    references: Optional[Iterable[PseudoneighborPair]] = None,
    *,
    max_steps: int = 50,
) -> list[Hole]:
    """
    Punch the holes generated by mapping each reference bridge backward.

    Each reference pair's bridge is mapped backward one step at a time: its two
    crossings move by :meth:`Trellis.image_cdist` onto the branch
    ``advance_key`` names, and the bridge whose consecutive crossings bracket
    that span (:func:`_containing_bridge`) gets a hole on the WHOLE bridge
    (INT-2: the image sub-arc's far end is a crossing with ``f^-1(W^s)``, past
    the trimmed stable manifold, so the trellis has nothing to cut it at). The
    recursion stops when the backward orbit becomes periodic — for a period-k
    orbit one backward step moves the image to the PREVIOUS unstable branch
    of the cycle, and a bridge lives on one branch, so meeting a visited
    bridge again means it has mapped onto itself after a multiple of k steps
    (for k = 1 exactly the "image lands in the same bridge" rule). The
    recursion also stops when no containing bridge is found (the image left
    the computed trellis) or after ``max_steps`` (safety guard).

    Requires the reference holes to exist already (call :func:`punch_holes`
    first); references without a hole are skipped.

    Args:
        trellis: The Trellis carrying pairs, bridges, and eigen/orbit data.
        references: Reference pairs to propagate. Defaults to the trellis's
            recorded reference pairs.
        max_steps: Upper bound on backward steps per reference.

    Returns:
        The newly punched holes. Every reference orbit punches its own chain —
        two orbits landing in the same bridge (even on the same side) each get
        their own hole; only a hole whose orbit and iterate are already
        represented by a recorded pair's hole is skipped as a true duplicate.
    """
    if references is None:
        references = [p for p in trellis.pseudoneighbors if p.is_reference]
    references = list(references)

    new_holes: list[Hole] = []
    seen: set[tuple[Optional[tuple[int, int]], Optional[int]]] = {
        (pair.origin, pair.iterate)
        for pair in trellis.pseudoneighbors
        if pair.hole is not None
    }
    # Regions already punched, per origin: a propagated hole must not land a
    # second time on a region its orbit already punched directly (a reference
    # and its +1 .. +(k-1) appearances are punched before the chains run).
    punched_regions: set[tuple] = {
        _region_key(pair.origin, pair.hole)
        for pair in trellis.pseudoneighbors
        if pair.hole is not None and pair.hole.bounding_ids is not None
    }

    orders: dict = {}
    for ref in references:
        origin = ref.origin if ref.origin is not None else ref.as_tuple()
        start = ref
        bridge = bridge_for_pair(trellis, ref) if ref.hole is not None else None
        if bridge is None:
            # The reference itself punched nothing (typically: no spanning
            # Bridge object between a parent-manifold crossing and a
            # blast-child crossing). Its orbit still owns backward holes —
            # continue the chain from the most-backward directly punched
            # pair instead.
            fallbacks = sorted(
                (
                    p
                    for p in trellis.pseudoneighbors
                    if p.origin == origin
                    and p.hole is not None
                    and p.iterate is not None
                ),
                key=lambda p: p.iterate,
            )
            start, bridge = next(
                (
                    (p, b)
                    for p in fallbacks
                    if (b := bridge_for_pair(trellis, p)) is not None
                ),
                (None, None),
            )
            if start is None:
                logger.warning(
                    "Orbit %s has no punched pair with a spanning bridge; "
                    "cannot propagate (run punch_holes() first?)",
                    origin,
                )
                continue
            logger.info(
                "Reference %s has no punched bridge; propagating backward "
                "from its iterate %+d pair %s instead",
                ref.as_tuple(),
                start.iterate,
                start.as_tuple(),
            )

        # The image span is tracked through the start bridge's two CROSSINGS,
        # each moved by image_cdist (the iterate table, else advance_key and
        # beta-scaling), on the branch the bridge itself lies on.
        start_key, start_ids = bridge.manifold_key, bridge.id
        advance = start_key[0].advance_key
        carried = np.asarray(start.hole.coords, dtype=np.float64)
        start_iterate = start.iterate or 0
        start_side = start.hole.bridge_side
        flips = not trellis.orientation_preserving
        # A bridge lives on one branch, so meeting it again is meeting it at the
        # same branch-cycle residue: the backward orbit has become periodic.
        visited = {bridge.id}

        for step in range(1, max_steps + 1):
            iterate = start_iterate - step
            span = tuple(sorted(
                trellis.image_cdist(end_id, -step, "unstable", branch=start_key)[1]
                for end_id in start_ids
            ))
            carried = _map_backward(trellis, carried)

            containing = _containing_bridge(trellis, span, advance(start_key, -step), orders)
            if containing is None:
                logger.debug(
                    "Backward image of %s left the computed bridges", ref.as_tuple()
                )
                break
            if containing.id in visited:
                break  # the backward orbit has become periodic — terminate
            visited.add(containing.id)

            key = (origin, iterate)
            if key not in seen:
                # The side is CARRIED along the orbit (author, 2026-10-07): the
                # map carries the bridge's dynamical orientation to the image
                # bridge's, flipping it every step only when it reverses
                # orientation. The carried point is measured against it.
                side = start_side
                if flips and step % 2:
                    side = OPPOSITE_SIDE[side]
                hole = _punch_in_bridge(
                    trellis, containing, bridge_side=side,
                    iterate=iterate, origin=origin, carried=carried, span=span,
                )
                if hole is not None:
                    region = _region_key(origin, hole)
                    if region not in punched_regions:
                        punched_regions.add(region)
                        seen.add(key)
                        new_holes.append(hole)

    logger.debug("Propagated %d hole(s) from %d reference(s)", len(new_holes), len(references))
    return new_holes


def partition_stable_manifold(
    trellis: "Trellis",
    branch_key: "ManifoldKey",
    side: Side,
) -> StablePartitionResult:
    """
    Partition one stable branch by the holes punched on one of its sides.

    Every punched hole on this branch and side participates. Walking from the
    anchor point outward, the hole regions' bounding intersections delimit the
    intervals; each hole opens the interval on the hole's side of each of its
    two bounding intersections (an interval bounding a hole is open at that
    end), every other end is closed, and a boundary point excluded by both
    neighbouring intervals — hole regions on both of its sides, or a hole
    abutting the branch end or the anchor — becomes a degenerate closed
    singleton (see the module Dev Notes).

    Each interval is stamped with its element identity (``element_id``,
    ``branch_key``, ``side``) and the two element lookup tables are built on top
    of it: :attr:`StablePartitionResult.element_of_intersection` (which element
    owns each crossing on the branch) and
    :attr:`StablePartitionResult.elements_at_bridge` (which elements sit at each
    bridge's two endpoints).

    Args:
        trellis: The Trellis carrying the holes to partition by.
        branch_key: Key of the stable branch to partition.
        side: Which side's holes to use.

    Returns:
        The StablePartitionResult for this branch and side.

    Raises:
        ValueError: If branch_key does not name a stable branch of the trellis.
        AssertionError: If a crossing on the branch is owned by anything other
            than exactly one element (the partition would not be a partition).
    """
    branch = trellis.branch(branch_key)
    if branch is None or branch.stability != "stable":
        raise ValueError(f"{branch_key} is not a stable branch of this trellis")

    marks = _hole_marks_on_branch(trellis, branch_key, side)
    intervals = _build_intervals(trellis, branch, marks)
    for element_id, interval in enumerate(intervals):
        interval.element_id = element_id
        interval.branch_key = branch_key
        interval.side = side

    element_of_intersection = _element_of_intersection(trellis, branch, intervals)
    return StablePartitionResult(
        branch_key=branch_key,
        side=side,
        intervals=intervals,
        element_of_intersection=element_of_intersection,
        elements_at_bridge=_elements_at_bridge(
            trellis, branch_key, element_of_intersection
        ),
    )


# ── geometric / bookkeeping helpers ─────────────────────────────────────────


def _region_key(
    origin: Optional[tuple[int, int]], hole: Hole
) -> tuple[Optional[tuple[int, int]], tuple[int, ...], Optional[Side]]:
    """Identity of the region one orbit punched, for propagation de-duplication.

    The two defining crossings alone do not name a region: a bridge has a
    region on EACH of its sides, and one orbit can legitimately land in both
    (the author's "singleton bridge with a hole on either side"). The side is
    therefore part of the key.
    """
    bounding = tuple(sorted(hole.bounding_ids)) if hole.bounding_ids else ()
    return (origin, bounding, hole.bridge_side)


def _side_of(tangent: NDArray[np.float64], displacement: NDArray[np.float64]) -> Optional[Side]:
    """Side of the displacement relative to the tangent, standard orientation.

    Standing on the manifold facing along ``tangent`` (toward the anchor),
    positive cross = left hand — the author's dynamical-direction convention
    (2026-07-09), which supersedes the earlier flipped calibration.
    """
    cross = float(tangent[0] * displacement[1] - tangent[1] * displacement[0])
    if cross > 0.0:
        return "left"
    if cross < 0.0:
        return "right"
    return None


def _anchorward_look(
    trellis: "Trellis", ix: "Intersection"
) -> Optional[NDArray[np.float64]]:
    """Unit stable direction toward the anchor at a crossing, read off the manifold.

    Toward the nearest node below the crossing's stable cdist, else away from
    the nearest one above (node cdists are monotone along a branch, so both are
    a bisect); None without manifold nodes.
    """
    manifold = trellis.manifolds.get(ix.manifold_b_key)
    if manifold is None:
        return None
    nodes = manifold.get_point_array(return_nodes=True)
    cdists = [node.cdist for node in nodes]
    tol = trellis.registry.cdist_tol
    here = np.asarray(ix.get_point(), dtype=np.float64)
    below = bisect_left(cdists, ix.stable_cdist - tol) - 1
    above = bisect_right(cdists, ix.stable_cdist + tol)
    if below >= 0:
        look = np.asarray(nodes[below].get_point(), dtype=np.float64) - here
    elif above < len(nodes):
        look = here - np.asarray(nodes[above].get_point(), dtype=np.float64)
    else:
        return None
    norm = float(np.linalg.norm(look))
    return look / norm if norm > 0.0 else None


def _approach_displacement(
    poly: NDArray[np.float64], q: NDArray[np.float64]
) -> Optional[NDArray[np.float64]]:
    """How the arc approaches its end nearest ``q``: the vector from ``q`` into the bridge.

    The end node sits just past the crossing, so the vector runs to the second
    node clearly distinct from ``q`` (the first when there is only one). On a
    two-node polyline that is the FAR end, which is why this reading only
    flags (NUM-1).
    """
    eps = 1e-9 * (1.0 + float(np.linalg.norm(poly[-1] - poly[0])))
    if np.linalg.norm(poly[0] - q) > np.linalg.norm(poly[-1] - q):
        poly = poly[::-1]
    distinct = np.flatnonzero(np.linalg.norm(poly[1:] - q, axis=1) > eps)[:2]
    return poly[1 + distinct[-1]] - q if len(distinct) else None


def _oriented_bridge_polyline(
    trellis: "Trellis", bridge: "Bridge"
) -> Optional[NDArray[np.float64]]:
    """The bridge polyline oriented by the unstable dynamical direction.

    A bridge is unstable manifold, so its dynamical direction — the direction
    of forward flow — points away from the fixed point, i.e. along increasing
    unstable canonical distance. Storage order (root→tail) need not match;
    this is the single place that mismatch is normalized: the polyline is
    reversed when the first endpoint's unstable cdist exceeds the second's.
    Endpoint cdists always exist, even for keyless iterated-bridge children,
    so no manifold key is consulted.

    Args:
        trellis: The Trellis resolving the bridge's endpoint intersections.
        bridge: The bridge whose polyline to orient.

    Returns:
        The (N, 2) polyline in dynamical orientation, or None when the bridge
        has fewer than two points, an unresolved endpoint, or endpoints of
        equal unstable cdist (orientation undecidable — logged as a warning).

    Note:
        Memoised per trellis and per bridge version
        (:meth:`Trellis.cached_bridge_polyline`); the result is shared, so it
        must not be written to.
    """
    return trellis.cached_bridge_polyline(
        bridge, lambda: _build_oriented_bridge_polyline(trellis, bridge)
    )


def _build_oriented_bridge_polyline(
    trellis: "Trellis", bridge: "Bridge"
) -> Optional[NDArray[np.float64]]:
    """The uncached body of :func:`_oriented_bridge_polyline`.

    Resolves the bridge's two endpoint canonical distances against the trellis
    and hands the orientation itself to
    :func:`~tanglepack.numerics.geometry.oriented_bridge_polyline`, which is the
    single implementation shared with the region layer.
    """
    if bridge.first_intersection is None or bridge.second_intersection is None:
        return None
    return oriented_bridge_polyline(
        bridge,
        trellis.intersection(bridge.first_intersection).unstable_cdist,
        trellis.intersection(bridge.second_intersection).unstable_cdist,
        tol=trellis.registry.cdist_tol,
    )


def _bridge_side_of(
    trellis: "Trellis",
    bridge: "Bridge",
    point: NDArray[np.float64],
    span: Optional[tuple[float, float]] = None,
) -> Optional[Side]:
    """Which side of ``bridge`` the point sits on, in dynamical orientation.

    Standing on the bridge looking along the unstable dynamical direction
    (away from the fixed point, increasing unstable cdist), positive
    cross(tangent, displacement) = left — the same convention as
    :func:`_side_of` on the stable manifold. The side is read at the nearest
    SEGMENT (the foot of the perpendicular, clipped to the segment), so a
    sub-arc lying inside one segment of a coarse bridge still answers.

    Args:
        trellis: The Trellis resolving the bridge's endpoint intersections.
        bridge: The bridge to classify against.
        point: Phase-space point to classify.
        span: Optional unstable-cdist extent of the sub-arc to read the side
            against (a propagated hole's image sub-arc); the whole bridge
            when omitted.

    Returns:
        ``"left"`` or ``"right"``, or None when the bridge cannot be oriented
        or the point sits on the arc.
    """
    poly = _oriented_bridge_polyline(trellis, bridge)
    if poly is None:
        return None
    if span is not None:
        # The image sub-arc only: a far fold of the same bridge can sit nearer
        # the point than the sub-arc it belongs to (2026-10-02).
        nodes = bridge.get_point_array(return_nodes=True)
        cdists = sorted(node.cdist for node in nodes)
        lo = max(bisect_left(cdists, span[0]) - 1, 0)
        poly = poly[lo : bisect_right(cdists, span[1]) + 1]
    point = np.asarray(point, dtype=np.float64)
    start, step = poly[:-1], np.diff(poly, axis=0)
    along = np.einsum("ij,ij->i", point - start, step)
    length2 = np.einsum("ij,ij->i", step, step)
    t = np.clip(along / np.maximum(length2, np.finfo(np.float64).tiny), 0.0, 1.0)
    foot = start + t[:, None] * step
    nearest = int(np.argmin(np.linalg.norm(point - foot, axis=1)))
    return _side_of(step[nearest], point - foot[nearest])


def row_of_end(
    trellis: "Trellis", bridge_id: "BridgeId", endpoint: Endpoint
) -> Side:
    """
    The row of a bridge end, read off the crossing sign alone.

    The *row* is the side of the stable branch, in the stable dynamical
    orientation (looking toward the anchor), that the bridge approaches its
    crossing from. It DECIDES the row (D6); :func:`bridge_row_diagnosis`
    reads the same quantity off the polylines, as a resolution flag only.

    Four manifold rays leave a crossing ``p``, and
    :attr:`~tanglepack.numerics.Intersection.Intersection.crossing_sign` is
    ``s = sign(cross(u+, s+))``. The stable look direction is ``s-``, so
    ``cross(s-, u+) = cross(u+, s+) = s``: the outward-unstable ray ``u+`` lies
    on the ``left`` of the stable dynamical direction exactly when ``s > 0``,
    and ``u-`` on the opposite side. A
    :data:`~tanglepack.numerics.Bridge.BridgeId` is ordered by unstable
    canonical distance, so the bridge leaves its ``"first"`` end along ``u+``
    and enters its ``"second"`` end along ``u-``. Hence
    ``row(first) = left iff sign(first) > 0`` and
    ``row(second) = left iff sign(second) < 0``.

    It needs no manifold nodes, so it answers at an anchor and on a
    polyline of two nodes alike, and it never returns None. For a same-branch
    bridge the two ends share a row exactly when the two crossing signs
    differ (:func:`bridge_row_violation`).

    Args:
        trellis: The Trellis resolving the bridge's endpoint crossings.
        bridge_id: The bridge's :data:`~tanglepack.numerics.Bridge.BridgeId`.
        endpoint: ``"first"`` or ``"second"`` — which end of the id to read.

    Returns:
        ``"left"`` or ``"right"``.

    Raises:
        ValueError: If ``endpoint`` names neither end, or if the crossing's
            ``crossing_sign`` is 0 (undetermined handedness — a tangency or a
            crossing registered without one), in which case no side of the
            stable branch is defined.
    """
    intersection_id = bridge_id[endpoint_index(endpoint)]
    sign = trellis.intersection(intersection_id).crossing_sign
    if sign == 0:
        raise ValueError(
            f"crossing {intersection_id} of bridge {bridge_id} has no crossing "
            "sign, so the side of the stable branch its bridge approaches from "
            "is undetermined"
        )
    outward = sign > 0 if endpoint == "first" else sign < 0
    return "left" if outward else "right"


def rows_of_bridge(
    trellis: "Trellis", bridge_id: "BridgeId"
) -> tuple[Side, Side]:
    """
    The rows at both ends of a bridge (see :func:`row_of_end`).

    Args:
        trellis: The Trellis resolving the bridge's endpoint crossings.
        bridge_id: The bridge's :data:`~tanglepack.numerics.Bridge.BridgeId`.

    Returns:
        ``(row at "first", row at "second")``.

    Raises:
        ValueError: If either crossing has no crossing sign.
    """
    return (
        row_of_end(trellis, bridge_id, "first"),
        row_of_end(trellis, bridge_id, "second"),
    )


@dataclass(frozen=True)
class RowEnd:
    """
    One end of a bridge, read by its crossing sign and by its polyline.

    Attributes:
        intersection_id: Registry id of the crossing.
        crossing_sign: Its ``crossing_sign``.
        combinatorial_row: :func:`row_of_end` (None at a zero sign), the row
            that decides (D6).
        geometric_row: The side of the stable branch, looking toward the
            anchor, the polyline approaches the crossing from (None at an
            anchor, which has no stable frame of its own, or when undecidable).
        approach_angle: Degrees between the arc's approach displacement and
            the stable tangent, in ``[0, 90]`` (None when undecidable). A
            fraction of a degree is a near-tangent end.
        displacement: Length of that approach displacement (None when
            undecidable).
        unstable_cdist: The crossing's unstable canonical distance.
        stable_cdist: The crossing's stable canonical distance.
    """

    intersection_id: int
    crossing_sign: float
    combinatorial_row: Optional[Side]
    geometric_row: Optional[Side]
    approach_angle: Optional[float]
    displacement: Optional[float]
    unstable_cdist: float
    stable_cdist: float

    def describe(self) -> str:
        """One line naming every reading of this end."""
        angle = "-" if self.approach_angle is None else f"{self.approach_angle:.3g} deg"
        disp = "-" if self.displacement is None else f"{self.displacement:.3g}"
        return (
            f"crossing {self.intersection_id}: sign {self.crossing_sign:+g}, "
            f"combinatorial row {self.combinatorial_row}, geometric row "
            f"{self.geometric_row}, approach {angle} over |d| {disp}, cdists "
            f"(u {self.unstable_cdist!r}, s {self.stable_cdist!r})"
        )


@dataclass(frozen=True)
class RowDiagnosis:
    """
    Every reading of a bridge's two rows, for the row assertion and the flag.

    Attributes:
        bridge_id: The two crossing ids, ordered by unstable cdist (the
            :data:`~tanglepack.numerics.Bridge.BridgeId` order).
        branch_key: The stable branch holding the first crossing (None when
            it lies on no branch of the trellis).
        ends: The ``"first"`` and ``"second"`` ends, in ``bridge_id`` order.
        n_nodes: Nodes in the bridge polyline the geometric rows are read on.
        chord: Distance between the two crossings.
    """

    bridge_id: tuple[int, int]
    branch_key: Optional["ManifoldKey"]
    ends: tuple[RowEnd, RowEnd]
    n_nodes: int
    chord: float

    @property
    def signs_alternate(self) -> bool:
        """True when the two crossing signs differ (the combinatorial rows then agree)."""
        return self.ends[0].crossing_sign * self.ends[1].crossing_sign < 0

    @property
    def geometry_disagrees(self) -> bool:
        """Some end's polyline reads a row other than its crossing sign: under-resolved."""
        return any(
            end.geometric_row is not None and end.geometric_row != end.combinatorial_row
            for end in self.ends
        )

    @property
    def verdict(self) -> str:
        """Which reading to distrust."""
        if self.signs_alternate:
            angles = [e.approach_angle for e in self.ends if e.approach_angle is not None]
            angle = f"{min(angles):.3g} deg" if angles else "an undecidable angle"
            return (
                "signs alternate: the combinatorial rows agree; the geometric "
                f"reading on a {self.n_nodes}-node polyline at {angle} is likely "
                "below resolution"
            )
        return "signs equal: the combinatorial rows disagree as well"

    def describe(self) -> str:
        """The diagnosis as text: the bridge, both ends, and the verdict."""
        key = self.branch_key
        branch = "None" if key is None else f"({key[1]}, {key[2]}, {key[3]})"
        return "\n".join(
            [
                f"bridge {self.bridge_id} on stable branch {branch}: "
                f"{self.n_nodes} node(s), chord {self.chord:.3g}",
                *(f"  {end.describe()}" for end in self.ends),
                f"  verdict: {self.verdict}",
            ]
        )


def bridge_row_diagnosis(
    trellis: "Trellis", bridge: "Bridge"
) -> Optional[RowDiagnosis]:
    """
    Read both ends of a bridge by their crossing signs and by the polylines.

    The geometric row of an end is the side of the stable branch, looking
    toward the anchor (:func:`_anchorward_look`), of the arc's approach
    displacement there (:func:`_approach_displacement`). It is the
    cross-check of :func:`row_of_end`, never the decision (D6).

    Args:
        trellis: The Trellis carrying the intersections, branches and manifolds.
        bridge: The bridge to diagnose.

    Returns:
        The :class:`RowDiagnosis`, or None for a partial piece or a bridge
        that cannot be oriented.
    """
    poly = _oriented_bridge_polyline(trellis, bridge)
    if poly is None or bridge.id is None:
        return None
    ends = []
    for endpoint, iid in zip(("first", "second"), bridge.id):
        ix = trellis.intersection(iid)
        disp = look = None
        if not ix.is_anchor:
            disp = _approach_displacement(poly, np.asarray(ix.get_point(), dtype=np.float64))
            look = None if disp is None else _anchorward_look(trellis, ix)
        geometric_row = angle = displacement = None
        if look is not None:
            geometric_row = _side_of(look, disp)
            displacement = float(np.linalg.norm(disp))
            angle = float(np.degrees(np.arctan2(
                abs(float(look[0] * disp[1] - look[1] * disp[0])),
                abs(float(np.dot(look, disp))),
            )))
        ends.append(RowEnd(
            intersection_id=iid,
            crossing_sign=float(ix.crossing_sign),
            combinatorial_row=(
                row_of_end(trellis, bridge.id, endpoint) if ix.crossing_sign else None
            ),
            geometric_row=geometric_row,
            approach_angle=angle,
            displacement=displacement,
            unstable_cdist=ix.unstable_cdist,
            stable_cdist=ix.stable_cdist,
        ))
    branch = trellis.branch_containing(bridge.id[0], "stable")
    first, second = (
        np.asarray(trellis.intersection(i).get_point(), dtype=np.float64) for i in bridge.id
    )
    return RowDiagnosis(
        bridge_id=bridge.id,
        branch_key=None if branch is None else branch.key,
        ends=tuple(ends),
        n_nodes=len(poly),
        chord=float(np.linalg.norm(second - first)),
    )


def bridge_row_violation(
    trellis: "Trellis", bridge: "Bridge"
) -> Optional[str]:
    """
    Describe how a bridge breaks the row invariant, or None if it holds.

    A bridge whose two defining crossings lie on the SAME stable branch must
    approach both of them from the same side of that branch: the lobe lies
    entirely on one side of the stable manifold, because the endpoints are
    consecutive crossings and the arc between them therefore never crosses it
    again. The rows are the crossing signs' (:func:`rows_of_bridge`, D6), so
    the invariant fails exactly when the two signs are EQUAL. Bridges whose
    crossings sit on different stable branches (heteroclinic, or a period > 1
    orbit) are exempt, and so is a partial piece. A polyline that reads other
    rows while the signs alternate is under-resolved; :func:`geometry_flags`
    reports it.

    Args:
        trellis: The Trellis carrying the intersections, branches and manifolds.
        bridge: The bridge to check.

    Returns:
        A human-readable violation message followed by the bridge's
        :class:`RowDiagnosis` (built only on failure), or None when the
        invariant holds or the check does not apply.
    """
    if bridge.id is None:
        return None
    first, second = bridge.id
    unstable_key = trellis.unstable_key_of(bridge.id)
    branch = trellis.branch(trellis.stable_key_at(first, unstable_key))
    if branch is None or branch is not trellis.branch(trellis.stable_key_at(second, unstable_key)):
        return None  # different stable branches (heteroclinic / period > 1)
    row_first, row_second = rows_of_bridge(trellis, bridge.id)
    if row_first == row_second:
        return None
    diagnosis = bridge_row_diagnosis(trellis, bridge)
    return (
        f"Bridge {bridge.id} on stable branch (orbit {branch.orbit_index}, branch "
        f"{branch.branch_index}) approaches its two defining crossings from opposite "
        f"sides of that branch (equal crossing signs): {first} -> {row_first}, "
        f"{second} -> {row_second}. A bridge's endpoints are consecutive crossings, "
        f"so the arc between them cannot cross the stable manifold again."
        + ("" if diagnosis is None else f"\n{diagnosis.describe()}")
    )


def check_bridge_rows_consistent(trellis: "Trellis", bridge: "Bridge") -> None:
    """
    Assert that a bridge's two ends agree on their row.

    See :func:`bridge_row_violation` for the invariant and the exemptions.

    Args:
        trellis: The Trellis carrying the intersections, branches and manifolds.
        bridge: The bridge to check.

    Raises:
        AssertionError: If the bridge's two crossings share a stable branch and
            carry equal crossing signs (opposite rows).
    """
    message = bridge_row_violation(trellis, bridge)
    assert message is None, message


def _orbit_side(hole: Hole, orientation_preserving: bool) -> Optional[Side]:
    """A hole's bridge side normalized to iterate 0 of its orbit.

    Under an orientation-preserving map the side is already invariant along
    the orbit; under an orientation-reversing one it flips every step, so an
    odd iterate is normalized by flipping it back.
    """
    if hole.bridge_side is None:
        return None
    if orientation_preserving:
        return hole.bridge_side
    if hole.iterate is None:
        return None
    return hole.bridge_side if hole.iterate % 2 == 0 else OPPOSITE_SIDE[hole.bridge_side]


@dataclass(frozen=True)
class SideDiagnosis:
    """
    One hole that fails a hole-side check, with every reading of its side.

    Attributes:
        bounding_ids: The hole's bridge (its two defining crossings).
        origin: The hole's orbit (``Hole.origin``).
        iterate: The hole's place on it (``Hole.iterate``).
        side_source: Where its ``bridge_side`` came from (``Hole.side_source``).
        bridge_side: The stored side.
        measured_side: The side measured on the curves (``Hole.measured_side``).
        expected_side: The orbit's side at this iterate, from the first hole of
            the origin (parity-flipped under an orientation-reversing map);
            None when the orbit check does not apply.
        reference_iterate: The iterate of that first hole.
    """

    bounding_ids: Optional[tuple[int, int]]
    origin: Optional[tuple[int, int]]
    iterate: Optional[int]
    side_source: Optional[SideSource]
    bridge_side: Optional[Side]
    measured_side: Optional[Side]
    expected_side: Optional[Side]
    reference_iterate: Optional[int]

    @property
    def orbit_disagrees(self) -> bool:
        """The stored side breaks the orbit invariant."""
        return self.expected_side is not None and self.bridge_side != self.expected_side

    @property
    def measured_disagrees(self) -> bool:
        """The measured side contradicts the stored (carried) one."""
        return self.measured_side is not None and self.measured_side != self.bridge_side

    def describe(self) -> str:
        """One line naming the hole and each reading that disagrees."""
        parts = [
            f"hole {self.bounding_ids} of origin {self.origin} at iterate "
            f"{self.iterate} ({self.side_source}): bridge_side {self.bridge_side}"
        ]
        if self.orbit_disagrees:
            parts.append(
                f"its orbit expects {self.expected_side} (from iterate "
                f"{self.reference_iterate})"
            )
        if self.measured_disagrees:
            parts.append(f"measured {self.measured_side} on its image sub-arc")
        return ", ".join(parts)


def side_diagnoses(
    holes: Iterable[Hole], *, orientation_preserving: bool = True
) -> list[SideDiagnosis]:
    """
    Diagnose every hole that disagrees with its orbit or with its measured side.

    The ORBIT check: a propagated hole is the backward image of its orbit's
    reference hole, and the map carries a bridge's dynamical orientation to
    the image bridge's, so an orientation-PRESERVING map (det J > 0 — both
    Hénon fixtures have b = 1) leaves the left/right classification fixed
    along the whole backward chain. Under an orientation-REVERSING map the
    side alternates with the parity of :attr:`Hole.iterate`, which
    ``orientation_preserving=False`` encodes. Since propagated holes CARRY
    their side (author, 2026-10-07) this guards the direct holes of one
    orbit (the reference and its ``+1 .. +(k-1)`` appearances, each read off
    its own crossing signs) and the carrying itself.

    The MEASURED reading is the geometric cross-check: a hole whose
    ``measured_side`` (the reference hole carried back by the real inverse
    map, read on the image sub-arc) contradicts its stored side is
    under-resolved. It is a flag (:func:`geometry_flags`), never an error (D6).

    Holes with no ``origin`` or no comparable side skip the orbit check, and
    holes with no ``measured_side`` skip the measured one; both counts are
    logged at INFO.

    Args:
        holes: The holes to check, typically ``trellis.holes``.
        orientation_preserving: True when the map has det J > 0.

    Returns:
        One :class:`SideDiagnosis` per disagreeing hole, in input order.
    """
    holes = list(holes)
    orbit_side: dict[tuple[int, int], tuple[Side, Optional[int]]] = {}
    unplaced = 0
    for hole in holes:
        side = _orbit_side(hole, orientation_preserving)
        if hole.origin is None or side is None:
            unplaced += 1
        else:
            orbit_side.setdefault(hole.origin, (side, hole.iterate))
    unmeasured = sum(
        1 for hole in holes if hole.side_source == "carried" and hole.measured_side is None
    )
    if unplaced or unmeasured:
        logger.info(
            "Hole-side check skips %d hole(s) with no comparable orbit side and "
            "%d propagated hole(s) with no measured side", unplaced, unmeasured
        )

    diagnoses: list[SideDiagnosis] = []
    for hole in holes:
        expected = reference_iterate = None
        if hole.origin in orbit_side and _orbit_side(hole, orientation_preserving):
            expected, reference_iterate = orbit_side[hole.origin]
            if not orientation_preserving and hole.iterate % 2:
                expected = OPPOSITE_SIDE[expected]
        diagnosis = SideDiagnosis(
            bounding_ids=hole.bounding_ids,
            origin=hole.origin,
            iterate=hole.iterate,
            side_source=hole.side_source,
            bridge_side=hole.bridge_side,
            measured_side=hole.measured_side,
            expected_side=expected,
            reference_iterate=reference_iterate,
        )
        if diagnosis.orbit_disagrees or diagnosis.measured_disagrees:
            diagnoses.append(diagnosis)
    return diagnoses


def bridge_side_violations(
    holes: Iterable[Hole], *, orientation_preserving: bool = True
) -> list[str]:
    """
    Describe every hole whose side of its bridge breaks the orbit invariant (I1).

    See :func:`side_diagnoses` for the check and which holes it skips; a
    measured side that disagrees is a flag, not a violation.

    Args:
        holes: The holes to check, typically ``trellis.holes``.
        orientation_preserving: True when the map has det J > 0.

    Returns:
        One line per hole whose orbit disagrees (its :class:`SideDiagnosis`).
        Empty when I1 holds.
    """
    return [
        diagnosis.describe()
        for diagnosis in side_diagnoses(holes, orientation_preserving=orientation_preserving)
        if diagnosis.orbit_disagrees
    ]


def check_holes_share_bridge_side(
    holes: Iterable[Hole], *, orientation_preserving: bool = True
) -> None:
    """
    Assert that all holes of one origin share their side of the bridge (I1).

    See :func:`side_diagnoses` for the parity rule for orientation-reversing
    maps and which holes are skipped. Measured sides are not compared here:
    they are the geometric flag (:func:`geometry_flags`).

    Args:
        holes: The holes to check, typically ``trellis.holes``.
        orientation_preserving: True when the map has det J > 0.

    Raises:
        AssertionError: If two holes of the same origin sit on different sides
            of their bridges (up to the parity flip when the map reverses
            orientation). The message lists every :class:`SideDiagnosis`.
    """
    messages = bridge_side_violations(
        holes, orientation_preserving=orientation_preserving
    )
    assert not messages, "; ".join(messages)


def geometry_flags(
    trellis: "Trellis", bridges: Iterable["Bridge"], holes: Iterable[Hole]
) -> tuple[list[RowDiagnosis], list[SideDiagnosis]]:
    """
    The geometric readings that disagree with the sign rule (D6).

    The crossing signs decide rows and openings and the orbit carries hole
    sides; the polylines only cross-check them. A bridge whose polyline reads
    another row than its signs at some end, or a propagated hole whose
    carried point lies on the other side of its image sub-arc, marks an
    under-resolved polyline (or a carried point that left its lobe): a reason
    to refine, never an error.

    Args:
        trellis: The Trellis carrying the bridges and manifolds.
        bridges: The bridges to read, typically the hole bridges.
        holes: The holes whose measured sides to compare.

    Returns:
        ``(rows, sides)``: the :class:`RowDiagnosis` of every bridge whose
        ``geometry_disagrees``, and the :class:`SideDiagnosis` of every hole
        whose ``measured_disagrees``.
    """
    rows = [
        diagnosis
        for diagnosis in (bridge_row_diagnosis(trellis, bridge) for bridge in bridges)
        if diagnosis is not None and diagnosis.geometry_disagrees
    ]
    sides = [
        diagnosis
        for diagnosis in side_diagnoses(
            holes, orientation_preserving=trellis.orientation_preserving
        )
        if diagnosis.measured_disagrees
    ]
    return rows, sides


def _warn_geometry_flags(
    trellis: "Trellis", bridges: Iterable["Bridge"], holes: Iterable[Hole]
) -> None:
    """One WARNING per kind of :func:`geometry_flags` (a resolution flag, D6)."""
    rows, sides = geometry_flags(trellis, bridges, holes)
    if rows:
        logger.warning(
            "%d hole bridge(s) read a geometric row other than their crossing "
            "signs: under-resolved polylines, flag only (D6); %s\n%s",
            len(rows),
            [diagnosis.bridge_id for diagnosis in rows[:8]],
            rows[0].describe(),
        )
    if sides:
        logger.warning(
            "%d propagated hole(s) measured on the other side of their bridge than "
            "their orbit carries: under-resolved lobes, flag only (D6); %s",
            len(sides),
            "; ".join(diagnosis.describe() for diagnosis in sides[:5]),
        )


def _hole_openings(
    trellis: "Trellis", bridge: "Bridge", bridge_side: Side
) -> list[tuple[int, str, Side, "ManifoldKey"]]:
    """The stable halves a hole on ``bridge_side`` of ``bridge`` opens, one per bound (D6).

    The bridge passes each bound along ``u+`` (increasing unstable cdist), and
    the outward stable ray ``s+`` lies on the left of ``u+`` exactly when
    ``crossing_sign = sign(cross(u+, s+)) > 0``. The hole's side of the bridge
    therefore holds the outward half iff ``(sign > 0) == (side == "left")``, on
    the bridge's row there (:func:`row_of_end`). With alternating signs this is
    the inward pair of a direct hole (outward of the near bound, anchorward of
    the far one); a hole facing a trimmed tail opens its removed outward half,
    which still takes part in the singleton pinch. At an ANCHOR both stable
    branches leave the point: the hole faces the outward half of the one
    :func:`_anchor_branch` picks, and opens nothing there when that branch is
    not computed.

    Returns:
        The ``Hole.openings`` records ``(intersection_id, which, row,
        branch_key)``, each naming the stable branch it opens.
    """
    openings = []
    for endpoint, iid in zip(("first", "second"), bridge.id):
        ix = trellis.intersection(iid)
        if ix.is_anchor:
            key = _anchor_branch(trellis, bridge.id, bridge_side)
            if key is not None:
                openings.append((iid, "outward", bridge_side, key))
            continue
        outward = (ix.crossing_sign > 0) == (bridge_side == "left")
        openings.append((
            iid,
            "outward" if outward else "anchorward",
            row_of_end(trellis, bridge.id, endpoint),
            ix.manifold_b_key,
        ))
    return openings


def _anchor_branch(
    trellis: "Trellis", bridge_id: "BridgeId", row: Side
) -> Optional["ManifoldKey"]:
    """The stable branch at a bridge's anchor end on which the bridge's row is ``row``.

    The bridge of ``u(o, b)`` meets ``s(o, b)`` with the anchor's own sign
    and ``s(o, 1 - b)``, which runs the other way, with the opposite one;
    ``row(first)`` is left iff the sign is positive (:func:`row_of_end`).
    None when that branch is not in the trellis (a point seeded on one
    eigendirection: the hole faces the uncomputed half).
    """
    fp, _, orbit_index, branch_index = trellis.unstable_key_of(bridge_id)
    home_sign = trellis.intersection(bridge_id[0]).crossing_sign
    if (home_sign > 0) != (row == "left"):
        branch_index = 1 - branch_index
    key = (fp, "stable", orbit_index, branch_index)
    return key if key in trellis.branches else None


def _step_into(
    anchor: NDArray[np.float64],
    target: NDArray[np.float64],
    epsilon: float,
    chord_len: float,
) -> NDArray[np.float64]:
    """Step off a boundary-arc midpoint toward the opposite boundary's midpoint.

    The step is ``epsilon`` times the smaller of the pair chord and the
    arc-to-arc distance: the chord cap keeps a huge lobe (whose opposite
    boundary is far away, and whose interior curves off the straight line)
    from carrying the hole visibly away, the distance cap keeps a narrow lobe
    from being overshot.
    """
    direction = target - anchor
    dist = float(np.linalg.norm(direction))
    if dist == 0.0:
        return anchor
    return anchor + epsilon * min(chord_len, dist) * direction / dist


def _stable_arc_midpoint(
    trellis: "Trellis",
    near: "Intersection",
    far: "Intersection",
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Middle of the stable arc between two intersections, and its tangent.

    Walks the actual stable-manifold nodes between the pair's stable cdists
    (the trellis carries a live reference to the workbench manifolds), so the
    point lies ON the region's stable boundary even when the arc curves far
    from the pair's chord; the tangent is the local curve direction there,
    oriented toward the anchor (decreasing cdist). Falls back to the chord
    midpoint and chord direction when the manifold nodes are unavailable, and
    when the two bounds sit on DIFFERENT stable branches (heteroclinic, or a
    bridge of a period > 1 orbit): there is then no single arc between them,
    and walking one branch's nodes between two cdists measured on different
    branches would return an unrelated piece of curve.
    """
    key = near.manifold_b_key
    if key is None or key != far.manifold_b_key:
        logger.debug(
            "Stable bounds sit on different branches (%s, %s); using the chord "
            "midpoint",
            key,
            far.manifold_b_key,
        )
        return (
            (near.get_point() + far.get_point()) / 2.0,
            near.get_point() - far.get_point(),
        )
    manifold = trellis.manifolds.get(key)
    if manifold is not None:
        lo, hi = sorted((near.stable_cdist, far.stable_cdist))
        inside = [
            node.get_point()
            for node in manifold.get_point_array(return_nodes=True)
            if lo <= node.cdist <= hi
        ]
        if len(inside) >= 2:
            mid = len(inside) // 2
            prev_i, next_i = max(mid - 1, 0), min(mid + 1, len(inside) - 1)
            # Nodes run root -> tail = increasing cdist; toward the anchor
            # is the decreasing-cdist direction.
            tangent = np.asarray(inside[prev_i], dtype=np.float64) - np.asarray(
                inside[next_i], dtype=np.float64
            )
            return np.asarray(inside[mid], dtype=np.float64), tangent
        if inside:
            return (
                np.asarray(inside[0], dtype=np.float64),
                near.get_point() - far.get_point(),
            )
    return (
        (near.get_point() + far.get_point()) / 2.0,
        near.get_point() - far.get_point(),
    )


def _near_far(trellis: "Trellis", id_a: int, id_b: int) -> tuple[int, int]:
    """Order two intersection ids as (toward-anchor, outward).

    Stable canonical distances are only comparable along ONE stable branch, so
    the usual stable ordering applies only when both crossings sit on the same
    branch. When they sit on different branches (heteroclinic, or a bridge of a
    period > 1 orbit) the ordering falls back to unstable canonical distance —
    the unstable dynamical direction, which also runs away from the anchor and
    is comparable along the single unstable branch both crossings share, since
    the two are the defining crossings of one bridge.
    """
    a = trellis.intersection(id_a)
    b = trellis.intersection(id_b)
    if a.manifold_b_key is not None and a.manifold_b_key == b.manifold_b_key:
        return (id_a, id_b) if a.stable_cdist <= b.stable_cdist else (id_b, id_a)
    logger.debug(
        "Crossings %d and %d sit on different stable branches; ordering them "
        "by unstable cdist (the dynamical direction)",
        id_a,
        id_b,
    )
    return (id_a, id_b) if a.unstable_cdist <= b.unstable_cdist else (id_b, id_a)


def _nearest_arc_point(
    bridge: "Bridge", point: NDArray[np.float64]
) -> Optional[NDArray[np.float64]]:
    """The bridge node closest to ``point``, excluding the flanking end nodes.

    A bridge's root and tail sit just PAST its intersections (see Bridge), on
    the stable curve rather than across the lobe, so they are skipped when
    enough nodes exist.
    """
    points = bridge.get_point_array()
    if points is None or len(points) == 0:
        return None
    points = np.asarray(points, dtype=np.float64)
    if len(points) > 4:
        points = points[1:-1]
    distances = np.linalg.norm(points - point, axis=1)
    return points[int(np.argmin(distances))]


def _pair_fixed_point(trellis: "Trellis", pair: PseudoneighborPair) -> "FixedPoint":
    """The fixed point owning the pair's stable branch."""
    key = pair.branch_key or trellis.intersection(pair.intersection_a).manifold_b_key
    return key[0]


def _is_forward_beyond_fundamental(trellis: "Trellis", pair: PseudoneighborPair) -> bool:
    """True for a forward-iterated pair that has wrapped past the fundamental segment.

    A pair at ``iterate >= k_value`` is the image of one already on a fundamental
    segment (its own branch's, after a full branch return) and gets no direct
    hole; ``1 .. k_value - 1`` are the trajectory's first appearances on the
    other branches and do. References and backward iterates return False.
    """
    if pair.iterate is None or pair.iterate <= 0:
        return False
    k_value = getattr(_pair_fixed_point(trellis, pair), "k_value", 1)
    return pair.iterate >= k_value


def _containing_bridge(
    trellis: "Trellis",
    span: tuple[float, float],
    branch_key: "ManifoldKey",
    orders: dict,
) -> Optional["Bridge"]:
    """The bridge of ``branch_key`` whose two crossings bracket ``span``.

    Bridges tile an unstable branch between consecutive crossings, so this is
    a bisect over the branch's crossing order (anchor first), memoised per
    branch in ``orders`` for one propagation. A span that straddles a crossing,
    or lies past the last one, has no container; so has a pair of consecutive
    crossings no Bridge object spans.
    """
    if branch_key not in orders:
        ids = trellis.branch(branch_key).intersection_ids
        orders[branch_key] = (ids, [trellis.intersection(i).unstable_cdist for i in ids])
    ids, cdists = orders[branch_key]
    i = bisect_right(cdists, span[0]) - 1
    if i < 0 or i + 1 == len(ids) or span[1] > cdists[i + 1]:
        return None
    return trellis.bridge_between(ids[i], ids[i + 1])


def _map_backward(
    trellis: "Trellis", coords: NDArray[np.float64]
) -> Optional[NDArray[np.float64]]:
    """One backward map step of a phase-space point (None without a map)."""
    if coords is None or trellis.dynamical_system is None:
        return None
    return np.asarray(
        trellis.dynamical_system.map_inv(coords), dtype=np.float64
    ).ravel()


def _punch_in_bridge(
    trellis: "Trellis",
    bridge: "Bridge",
    *,
    bridge_side: Side,
    iterate: int,
    origin: Optional[tuple[int, int]],
    carried: Optional[NDArray[np.float64]],
    span: tuple[float, float],
) -> Optional[Hole]:
    """Punch a propagated hole between the two backward-imaged neighbors.

    ``span`` is the backward image's unstable-cdist extent inside the
    containing bridge — the two imaged pseudoneighbor points bound it, and the
    hole plots at the middle of that sub-arc (two orbits landing in the same
    bridge occupy different sub-arcs). ``bridge_side`` is the side carried
    from the orbit (``Hole.bridge_side``, ``side_source == "carried"``); the
    crossing signs turn it into the intervals the hole opens at the bridge's
    defining intersections (:func:`_hole_openings`). ``carried`` is the
    reference hole's coordinates mapped backward alongside — a point inside
    the true image region — which places the hole and is read on the image
    sub-arc as ``Hole.measured_side``, the geometric FLAG on the carried side.
    Without it (no inverse map) the containing bridge's own midpoint and
    chord place the hole.
    """
    near_id, far_id = _near_far(
        trellis, bridge.first_intersection, bridge.second_intersection
    )
    near = trellis.intersection(near_id)
    far = trellis.intersection(far_id)

    midpoint, u_tangent, arc_chord = _image_arc_midpoint(bridge, span)
    if midpoint is None:
        return None
    stable_mid, _stable_tangent = _stable_arc_midpoint(trellis, near, far)

    if carried is not None and u_tangent is not None:
        # Nudge off the unstable arc toward the carried point's side of it,
        # capped by the distance to the carried point so a narrow lobe is
        # not overshot.
        normal = np.array([-u_tangent[1], u_tangent[0]])
        norm = float(np.linalg.norm(normal))
        if norm > 0.0:
            normal /= norm
            if float(np.dot(normal, carried - midpoint)) < 0.0:
                normal = -normal
            step = min(
                0.05 * arc_chord,
                0.5 * float(np.linalg.norm(carried - midpoint)),
            )
            coords = midpoint + step * normal
        else:
            coords = _step_into(midpoint, stable_mid, 0.05, arc_chord)
    else:
        coords = _step_into(midpoint, stable_mid, 0.05, arc_chord)

    measured_side = (
        _bridge_side_of(trellis, bridge, carried, span) if carried is not None else None
    )
    return Hole(
        coords=(float(coords[0]), float(coords[1])),
        near_intersection_id=near_id,
        bridge_side=bridge_side,
        bounding_ids=(near_id, far_id),
        openings=_hole_openings(trellis, bridge, bridge_side),
        iterate=iterate,
        origin=origin,
        measured_side=measured_side,
        side_source="carried",
    )


def _image_arc_midpoint(
    bridge: "Bridge", span: Optional[tuple[float, float]]
) -> tuple[Optional[NDArray[np.float64]], Optional[NDArray[np.float64]], float]:
    """Middle of the image sub-arc of a bridge: (coords, local tangent, chord).

    ``span`` bounds the sub-arc in unstable cdist; the node closest to the
    span's middle is taken, with its local curve tangent and the straight-line
    chord length between the sub-arc's endpoints. Without a span the whole
    bridge is the arc.
    """
    nodes = bridge.get_point_array(return_nodes=True)
    if not nodes:
        return None, None, 0.0
    if span is None:
        lo_i, mid_i, hi_i = 0, len(nodes) // 2, len(nodes) - 1
    else:
        cdists = [node.cdist for node in nodes]
        u_lo, u_hi = span
        lo_i = min(range(len(nodes)), key=lambda i: abs(cdists[i] - u_lo))
        hi_i = min(range(len(nodes)), key=lambda i: abs(cdists[i] - u_hi))
        mid_i = min(
            range(len(nodes)), key=lambda i: abs(cdists[i] - (u_lo + u_hi) / 2.0)
        )
    midpoint = nodes[mid_i].get_point()
    prev_i, next_i = max(mid_i - 1, 0), min(mid_i + 1, len(nodes) - 1)
    tangent = (
        nodes[next_i].get_point() - nodes[prev_i].get_point()
        if next_i > prev_i
        else None
    )
    arc_chord = float(
        np.linalg.norm(nodes[hi_i].get_point() - nodes[lo_i].get_point())
    )
    return midpoint, tangent, arc_chord


def _hole_marks_on_branch(
    trellis: "Trellis",
    branch_key: "ManifoldKey",
    side: Side,
) -> tuple[dict[int, float], set[int], set[int]]:
    """
    Boundary intersections and open-interval marks on one branch and row.

    Every punched hole participates — a hole is a hole, with no exceptions by
    iterate or endpoint identity. Each opening record names the stable branch
    it opens (its intersection's own, or at the anchor the one the signs pick),
    so a bridge whose two defining intersections sit on different branches of
    a periodic orbit contributes to each one. A hole without openings does not
    participate.

    Returns:
        (boundaries, open_outward, open_anchorward): boundary intersection
        ids mapped to their stable cdists, and the ids whose outward /
        anchorward interval a hole opens.
    """
    boundaries: dict[int, float] = {}
    open_outward: set[int] = set()
    open_anchorward: set[int] = set()
    for hole in trellis.holes:
        if not hole.openings:
            logger.debug(
                "Hole at %s has no openings; it does not partition", hole.coords
            )
            continue
        for iid, which, row, opened in hole.openings:
            if (row, opened) != (side, branch_key):
                continue
            boundaries[iid] = trellis.intersection(iid).stable_cdist
            (open_outward if which == "outward" else open_anchorward).add(iid)
    return boundaries, open_outward, open_anchorward


def _build_intervals(
    trellis: "Trellis",
    branch: "TrellisBranch",
    marks: tuple[dict[int, float], set[int], set[int]],
) -> list[PartitionInterval]:
    """Assemble the ordered partition intervals from the hole marks.

    ``marks`` is :func:`_hole_marks_on_branch`'s output: the boundary
    intersections and, per boundary, whether a hole opens its outward and/or
    anchorward interval. A boundary point excluded by both neighbouring
    intervals owns itself as a closed singleton (see the module Dev Notes).
    """
    boundary_map, open_outward, open_anchorward = marks
    boundaries: list[tuple[Optional[int], float]] = list(boundary_map.items())
    outer_id = branch.ordered_ids(toward_anchor=True)[0] if len(branch) else None
    if outer_id is not None:
        boundaries.append((outer_id, trellis.intersection(outer_id).stable_cdist))

    # The partition starts at the anchor point (cdist 0). When a boundary
    # already sits there (the registered anchor, synthetic at exactly 0.0), it
    # plays that role itself; otherwise a synthetic anchor boundary is added.
    if not any(cdist == 0.0 for _, cdist in boundaries):
        tol = trellis.registry.cdist_tol
        near_anchor = [bid for bid, cdist in boundaries if cdist <= tol]
        if near_anchor:
            logger.info(
                "branch %s: no boundary at the anchor's exact cdist 0, but %s lie(s) "
                "within cdist_tol=%g of it; taken as the anchor",
                branch.key[1:],
                near_anchor,
                tol,
            )
        else:
            boundaries.append((None, 0.0))

    # Sort and deduplicate by id (the anchor and branch end keep their slots).
    unique: dict = {}
    for bid, cdist in boundaries:
        unique[bid if bid is not None else ("anchor", cdist)] = (bid, cdist)
    ordered = sorted(unique.values(), key=lambda b: b[1])
    if len(ordered) == 1:
        # C1: no crossing past the anchor yet -- one element to the end of the
        # computed branch.
        ((bid, cdist),) = ordered
        end = float(trellis.manifolds[branch.key].tail.get_cdist("stable"))
        logger.info(
            "branch %s has no crossing past its anchor; one element to its end "
            "(grow the unstable manifold)",
            branch.key[1:],
        )
        assert bid not in open_outward, (
            f"a hole opens outward of the anchor of {branch.key[1:]}, which has no "
            "crossing past it"
        )
        return [PartitionInterval(bid, None, cdist, end, True, True)]

    intervals: list[PartitionInterval] = []
    for (lo_id, lo_c), (hi_id, hi_c) in zip(ordered, ordered[1:]):
        intervals.append(
            PartitionInterval(
                lo_id, hi_id, lo_c, hi_c,
                closed_lo=lo_id not in open_outward,
                closed_hi=hi_id not in open_anchorward,
            )
        )

    # A boundary point excluded by both neighbouring intervals owns itself as
    # a closed singleton: an intersection with hole regions on both sides, the
    # outermost intersection when a hole abuts it from the interior, or the
    # anchor-artifact intersection when a hole abuts it from outside. Each
    # point thereby belongs to exactly one piece; an interior boundary always
    # bounds a hole, so at least one of its intervals is open there and "]["
    # cannot occur. Two holes flanking one bridge (one on each side) emit
    # complementary openings at both defining intersections, so this same
    # rule pinches each of them into a singleton with open intervals on both
    # sides — the "singleton bridge" case needs no extra code.
    with_singletons: list[PartitionInterval] = []
    for position, (bid, cdist) in enumerate(ordered):
        before = intervals[position - 1] if position > 0 else None
        after = intervals[position] if position < len(intervals) else None
        owned = (before is not None and before.closed_hi) or (
            after is not None and after.closed_lo
        )
        if not owned and (before is not None or after is not None):
            with_singletons.append(
                PartitionInterval(bid, bid, cdist, cdist, True, True)
            )
        if after is not None:
            with_singletons.append(after)
    return with_singletons


def span_contains(
    lo: float,
    hi: float,
    closed_lo: bool,
    closed_hi: bool,
    cdist: float,
    tol: float,
) -> bool:
    """
    Whether a (possibly half-open) canonical-distance span contains a point.

    Args:
        lo: Lower end of the span.
        hi: Upper end of the span.
        closed_lo: True if the lower end belongs to the span.
        closed_hi: True if the upper end belongs to the span.
        cdist: The canonical distance to test.
        tol: Comparison tolerance (the registry's ``cdist_tol``).

    Returns:
        True if the span contains the point.
    """
    if cdist < lo - tol or cdist > hi + tol:
        return False
    if abs(cdist - lo) <= tol:
        return closed_lo
    if abs(cdist - hi) <= tol:
        return closed_hi
    return True


def owns_cdist(interval: PartitionInterval, cdist: float, tol: float) -> bool:
    """
    Whether a partition element contains a point at this stable cdist.

    A point strictly inside the span belongs to the element; a point sitting on
    an end belongs to it only when that end is closed. Adjacent elements are
    never both closed at their shared end (see :func:`_build_intervals`), so
    exactly one element owns every point of the branch.

    Args:
        interval: The partition element (interval) to test.
        cdist: Stable canonical distance of the point.
        tol: Comparison tolerance (the registry's ``cdist_tol``).

    Returns:
        True if this element owns the point.
    """
    return span_contains(
        interval.lo_cdist,
        interval.hi_cdist,
        interval.closed_lo,
        interval.closed_hi,
        cdist,
        tol,
    )


def _element_of_intersection(
    trellis: "Trellis",
    branch: "TrellisBranch",
    intervals: list[PartitionInterval],
) -> dict[int, int]:
    """Map every crossing on the branch to the one element that owns it, exactly.

    A boundary crossing is owned by the element closed at its id (its singleton
    if it has one, else the closed neighbour); any other crossing by its exact
    stable cdist, no tolerance — every boundary is a registry float (asserted
    here), so two crossings closer than ``cdist_tol`` still fall on their own
    sides of a boundary (module Dev Notes, "Exact ownership"). The intervals
    tile the branch from the anchor to the outermost crossing, so every
    crossing lands in exactly one — asserted, since a crossing owned twice (or
    not at all) means the openings did not alternate and the result is not a
    partition.
    """
    closed_at: dict[int, list[int]] = {}
    for interval in intervals:
        for bid, cdist, closed in (
            (interval.lo_id, interval.lo_cdist, interval.closed_lo),
            (interval.hi_id, interval.hi_cdist, interval.closed_hi),
        ):
            if bid is None:
                continue  # the synthetic anchor at 0.0
            registered = float(trellis.intersection(bid).stable_cdist)
            assert cdist == registered, (
                f"boundary {bid} of element {interval.element_id} on {branch.key[1:]} "
                f"sits at stable cdist {cdist!r}, not its registry cdist {registered!r}"
            )
            owned = closed_at.setdefault(bid, [])
            if closed and interval.element_id not in owned:
                owned.append(interval.element_id)

    owners: dict[int, int] = {}
    for intersection_id in branch.intersection_ids:
        cdist = float(trellis.intersection(intersection_id).stable_cdist)
        matches = closed_at.get(intersection_id)
        if matches is None:
            matches = [
                interval.element_id
                for interval in intervals
                if owns_cdist(interval, cdist, 0.0)
            ]
        assert len(matches) == 1, (
            f"crossing {intersection_id} at stable cdist {cdist!r} on "
            f"{branch.key[1:]} is owned by {len(matches)} partition elements "
            f"{matches}, expected exactly one; candidates "
            + ", ".join(
                f"#{iv.element_id} {'[' if iv.closed_lo else '('}{iv.lo_cdist!r}, "
                f"{iv.hi_cdist!r}{']' if iv.closed_hi else ')'}"
                for iv in intervals
                if iv.lo_cdist <= cdist <= iv.hi_cdist
            )
        )
        owners[intersection_id] = matches[0]
    return owners


def _elements_at_bridge(
    trellis: "Trellis",
    branch_key: "ManifoldKey",
    element_of_intersection: dict[int, int],
) -> dict["BridgeId", tuple[Optional[int], Optional[int]]]:
    """Map each bridge touching this branch to the elements at its two endpoints.

    Only non-partial bridges have an id and appear here. An endpoint whose stable
    side lies on a different branch is ``None``: this result cannot name it, and
    the caller resolves it through the trellis (another branch of the same
    trellis) or the session (another trellis's branch).
    """
    at_bridge: dict["BridgeId", tuple[Optional[int], Optional[int]]] = {}
    for bridge in trellis.bridges:
        bridge_id = bridge.id
        if bridge_id is None:
            continue  # a partial arc is not a bridge and has no identity
        ends: list[Optional[int]] = []
        unstable_key = trellis.unstable_key_of(bridge_id)
        for endpoint in bridge_id:
            if trellis.stable_key_at(endpoint, unstable_key) != branch_key:
                ends.append(None)
                continue
            element_id = element_of_intersection.get(endpoint)
            assert element_id is not None, (
                f"crossing {endpoint} of bridge {bridge_id} lies on "
                f"{branch_key[1:]} but no partition element owns it"
            )
            ends.append(element_id)
        if ends[0] is None and ends[1] is None:
            continue  # neither end touches this branch
        at_bridge[bridge_id] = (ends[0], ends[1])
    return at_bridge
