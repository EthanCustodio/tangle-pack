"""
The stable-manifold partition: holes, intervals and their propagation.

Punches a hole in each region bounded by a reference pseudoneighbor pair and
its bridge, propagates the holes backward along the orbit, and cuts each
stable branch into the partition intervals the region layer reads.

Dev Notes — Stable Manifold Partition (Stable_Manifold_Partition_Algorithm.pdf)

The PDF is less formal than the Pseudoneighbor one; the following semantics were
fixed in conversation with the author (July 2026):

* Holes are combinatorial (B4b, 2026-10-13): a pseudoneighbor pair and its
  connecting unstable arc (a bridge) bound a closed region, and the hole IS
  ``(origin, iterate, bounding_ids, bridge_side, openings)`` -- the region's
  bridge, its side of it and the stable halves it opens. Nothing here places
  it in the plane or calls the map: where to draw it is
  ``plotting.hole_point``, a view. ``origin`` is the reference pair it
  descends from (the plot marker), ``iterate`` 0 for the reference and
  negative backward (the plot label).
  ``near_intersection_id`` is fixed by convention to the toward-anchor member
  (smaller stable cdist) — the choice is arbitrary but must be consistent.
* Partition elements (Phase 5): every interval is stamped with
  ``(branch_key, side, element_id)`` and one lookup table is built beside it,
  ``element_of_intersection`` (crossing -> owning element, by id at a
  boundary and by EXACT cdist elsewhere, see "Exact ownership" below). A
  bridge end's element is read from it by the endpoint id the BridgeId
  carries (``PartitionFamily.owner_at_end`` on the end's row,
  ``owner_of_intersection`` on any side), never from the bridge list.
* Side (author, 2026-07-16): a hole's side is left/right of its own BRIDGE,
  i.e. of the unstable manifold, in the bridge's dynamical orientation —
  forward flow points AWAY from the fixed point, so the bridge is oriented by
  increasing unstable cdist; standing on it looking that way, positive
  cross(tangent, displacement) = left (``diagnostics._side_of``). It is
  decided ONCE and stored as ``Hole.bridge_side``: a direct hole reads it off
  its crossing signs (2026-10-02), a propagated hole CARRIES it from the hole
  it was propagated from, parity-flipped under an orientation-reversing map
  (author, 2026-10-07). Holes know nothing of the
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
  geometry stays as a cross-check that FLAGS under-resolution, in
  ``topology/diagnostics.py`` (``geometry_flags``, one WARNING per kind per
  ``Trellis.punch_holes(diagnose=True)``): ``bridge_row_diagnosis`` reads
  each end's approach on the polyline, which on a 2-4 node bridge is the
  vector to the far node (NUM-1, the E11 p3 / nested_inner 1e-5 cells from
  blast 6), and the MEASURED side of a propagated hole reads the
  backward-carried reference point at the nearest segment of the image
  sub-arc, which escapes the lobe when the reference lobe is
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
  A SCALED span end (no table entry) is good to ``SCALING_RTOL`` only, so one
  that close to a crossing of the branch is one WARNING per branch: the
  container could be its neighbour (NUM-5). A table end is exact; it is then
  the container's own end crossing, and is not watched (B4b, 2026-10-13).
  The hole sits on the WHOLE containing bridge, by design (INT-2, checked
  2026-10-08 on every law case, the deep p3 runs and four deep E11 states:
  600+ lookups, every span a STRICT sub-arc of its container and at most one
  of its two ends a registered crossing; the other end is a crossing of the
  unstable branch with ``f^-1(W^s)``, beyond the trimmed stable manifold, so
  the trellis has no node to cut the sub-arc at). Two orbits landing in the
  same bridge differ only by their sub-arcs, and BOTH are punched (per the
  author: no cross-orbit deduplication, even on the same bridge and side);
  the reference point carried back by the real inverse map places a
  propagated hole in figures and is the diagnostics' measured-side flag,
  never part of the algorithm.
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
  the parity of ``Hole.iterate`` (``side_after``, the one rule, C3)
  (``check_holes_share_bridge_side(..., orientation_preserving=False)``);
  propagated holes carry their side, so I1 guards the direct holes of one
  orbit (its references) and the carrying itself; and (I2) a bridge whose two crossings lie on the same
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
  grown that far). PERIOD-k (author, 2026-10-09, reversing D5): the rule has
  no exemption. EVERY forward iterate is dropped (``_is_forward_iterate``),
  the ``+1 .. +(k_value - 1)`` appearances on the other branches of the orbit
  included; those branches get their holes from the backward chain, which
  walks the whole unstable branch cycle. Punching the ``+j`` appearances (D5,
  2026-10-07) added holes the period-3 figures showed as forward iterates
  and four inert classes at 0-5 blasts.
* Propagation start: every chain starts at its reference. No backward pair
  is ever recorded (a reference lies in the window ``S \\ f(S)``, so one of
  its preimages leaves T and the backward walk stops at once), which
  ``punch_holes`` asserts; the direct punch of backward pairs, its two
  de-duplications and the fallback that continued an orbit from its
  most-backward directly punched pair (author, 2026-07-09) never ran on any
  build and went in B4b (2026-10-13). A reference whose pair has no spanning
  Bridge object punches no hole and starts no chain (``punch_holes`` warns).
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
  child-inside-parent assertion of the iterated partition at 9 blasts.
  Since WP11 (C15) the topology reads no absolute cdist slack at all: an
  anchor is told by identity (``Intersection.is_anchor``, synthetic at
  exactly 0.0), ownership and containment are exact, and the one slack left
  is the RELATIVE error of a scaled cdist (``SCALING_RTOL`` at a scaled
  landing end, ``collision_rtol`` in ``Trellis.same_point``). Forcing the slack
  to 0 moved no output on any fixture. ``cdist_tol`` survives only as
  geometry -- clipping a polyline to its arc (``Arc._interior_polyline``,
  ``diagnostics._oriented_bridge_polyline``) -- never as identity or ownership.
* The anchor (WP3, 2026-10-08). A periodic point is ONE crossing on every
  branch attached to it (D4). Both stable branches leave it, so an anchor
  opening faces the outward half of the branch the signs pick
  (:func:`_hole_openings`) and every opening record names the branch it
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
from bisect import bisect_right
from itertools import count
from typing import Iterable, Optional, TYPE_CHECKING

from .Pseudoneighbor import _branch_text
from .Trellis import SCALING_RTOL
from .TopologyResults import (
    Endpoint,
    endpoint_index,
    Hole,
    PartitionInterval,
    PseudoneighborPair,
    Side,
    StablePartitionResult,
    side_after,
)

if TYPE_CHECKING:
    from .Trellis import Trellis
    from .TrellisBranch import TrellisBranch
    from ..numerics.Bridge import Bridge, BridgeId
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
) -> list[Hole]:
    """
    Punch a hole in the region bounded by each reference pseudoneighbor pair.

    The hole is the pair's bridge, its side of that bridge in the bridge's
    dynamical orientation (``Hole.bridge_side``) and the openings, both read
    off the crossing signs (D6); it carries its orbit identity
    (``origin``/``iterate``) and is attached to its pair (``pair.hole``).
    Where to draw it is :func:`plotting.hole_point`'s business.

    Holes propagate BACKWARD only: a pair recorded at any forward iterate is
    the image of a region still attached to the stable manifold and gets no
    direct hole (it is logged at DEBUG and left with ``pair.hole = None``),
    on a period-k orbit too, and no backward pair is ever recorded (the
    reference lies in the window ``S \\ f(S)``), so only the references
    (``iterate == 0``) are punched here; :func:`propagate_reference_holes`
    punches their backward chains. See the module Dev Notes.

    Args:
        trellis: The Trellis carrying the pairs and bridges.
        pairs: Pairs to punch holes for. Defaults to every pair recorded on the
            trellis (``trellis.pseudoneighbors``).

    Returns:
        The punched holes. Pairs with no spanning bridge are skipped with a
        warning; forward-iterate pairs are skipped silently (DEBUG).

    Raises:
        AssertionError: If a pair is recorded at a backward iterate.
    """
    if pairs is None:
        pairs = trellis.pseudoneighbors
    pairs = list(pairs)
    assert all(p.iterate >= 0 for p in pairs), "holes propagate backward only; no backward pair is recorded"

    # Punching is idempotent per call: a pair that gets no hole this time must
    # not keep the one an earlier call attached, because propagate_reference_holes
    # reads pair.hole to seed the already-punched orbits and regions.
    for pair in pairs:
        pair.hole = None

    # Holes propagate BACKWARD only (author's rule, 2026-09-16): a forward
    # iterate of a punched hole is the image of a region still attached to the
    # stable manifold and gets no hole of its own, on every branch of a
    # period-k orbit (author, 2026-10-09). See the module Dev Notes.
    forward = [pair for pair in pairs if _is_forward_iterate(pair)]
    pairs = [pair for pair in pairs if not _is_forward_iterate(pair)]
    if forward:
        logger.debug(
            "Not punching %d forward-iterate pair(s): "
            "holes propagate backward only; %s",
            len(forward),
            [(pair.as_tuple(), pair.iterate) for pair in forward],
        )

    holes: list[Hole] = []
    skipped: list[PseudoneighborPair] = []
    for pair in pairs:
        bridge = bridge_for_pair(trellis, pair)
        if bridge is None:
            skipped.append(pair)
            continue

        near_id, far_id = _near_far(trellis, pair.intersection_a, pair.intersection_b)
        near = trellis.intersection(near_id)
        far = trellis.intersection(far_id)

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
            near_intersection_id=near_id,
            pair=pair,
            bridge_side=bridge_side,
            bounding_ids=(near_id, far_id),
            openings=openings,
            iterate=pair.iterate,
            origin=pair.origin,
        )
        pair.hole = hole
        holes.append(hole)

    if skipped:
        named = [
            f"{pair.as_tuple()} (iterate {pair.iterate} of {pair.origin}) on "
            + _branch_text(pair.branch_key or trellis.intersection(pair.intersection_a).manifold_b_key)
            for pair in skipped
        ]
        logger.warning(
            "Skipped %d pair(s) with no spanning bridge: the arc between them is not "
            "computed (the image of a bridge no blast iterated, D9), so their stable "
            "branches get no hole from them: %s", len(skipped), "; ".join(named),
        )
    logger.debug("Punched %d hole(s)", len(holes))
    return holes


def propagate_reference_holes(
    trellis: "Trellis",
    references: Optional[Iterable[PseudoneighborPair]] = None,
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
    the computed trellis). It needs no step cap: every step either stops or
    visits a new bridge, so it ends within the trellis's finite bridge count
    (PK-9).

    Requires the reference holes to exist already (call :func:`punch_holes`
    first); a reference without a hole (no Bridge spans its pair, which
    :func:`punch_holes` warned about) starts no chain.

    Args:
        trellis: The Trellis carrying pairs, bridges, and eigen/orbit data.
        references: Reference pairs to propagate. Defaults to the trellis's
            recorded reference pairs.

    Returns:
        The newly punched holes. Every reference orbit punches its own chain —
        two orbits landing in the same bridge (even on the same side) each get
        their own hole.
    """
    if references is None:
        references = [p for p in trellis.pseudoneighbors if p.is_reference]
    references = list(references)

    new_holes: list[Hole] = []
    orders: dict = {}
    warned: set = set()
    for ref in references:
        if ref.hole is None:
            continue
        origin = ref.origin if ref.origin is not None else ref.as_tuple()
        bridge = bridge_for_pair(trellis, ref)
        # The image span is tracked through the start bridge's two CROSSINGS,
        # each moved by image_cdist (the iterate table, else advance_key and
        # beta-scaling), on the branch the bridge itself lies on.
        start_key, start_ids = bridge.manifold_key, bridge.id
        advance = start_key[0].advance_key
        start_side = ref.hole.bridge_side
        preserving = trellis.orientation_preserving
        # A bridge lives on one branch, so meeting it again is meeting it at the
        # same branch-cycle residue: the backward orbit has become periodic.
        visited = {bridge.id}

        for step in count(1):
            # The span's interior is never a crossing of N (it lies past the pip cut); one end
            # may be a registered crossing, which is then the container's own end.
            ends = [trellis.image_cdist(end_id, -step, "unstable", branch=start_key) for end_id in start_ids]
            span = tuple(sorted(c for _, c, _ in ends))

            branch_key = advance(start_key, -step)
            containing = _containing_bridge(trellis, span, branch_key, orders)
            if containing is None:
                logger.debug(
                    "Backward image of %s left the computed bridges", ref.as_tuple()
                )
                break
            if branch_key not in warned:
                # A scaled end is good to SCALING_RTOL only: that close to a crossing
                # the bisect may have picked the neighbour bridge (NUM-5).
                ids, cdists = orders[branch_key]
                near = next((
                    (c, iid) for _, c, from_table in ends if not from_table
                    for iid, x in zip(ids, cdists) if abs(c - x) <= SCALING_RTOL * abs(x)
                ), None)
                if near is not None:
                    warned.add(branch_key)
                    logger.warning(
                        "Backward image %+d of reference %s on %s: its scaled end %.9g lies "
                        "within SCALING_RTOL of crossing %d, so its container %s could be the "
                        "neighbour bridge (NUM-5); grow or blast to register the image",
                        -step, ref.as_tuple(), _branch_text(branch_key), *near, containing.id,
                    )
            if containing.id in visited:
                break  # the backward orbit has become periodic — terminate
            visited.add(containing.id)

            # The side is CARRIED along the orbit (author, 2026-10-07): the
            # map carries the bridge's dynamical orientation to the image
            # bridge's, flipping it every step only when it reverses
            # orientation.
            new_holes.append(_punch_in_bridge(
                trellis, containing, bridge_side=side_after(start_side, -step, preserving),
                iterate=-step, origin=origin,
            ))

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
    ``branch_key``, ``side``) and the element lookup table is built on top
    of it: :attr:`StablePartitionResult.element_of_intersection` (which element
    owns each crossing on the branch).

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
    )


# ── geometric / bookkeeping helpers ─────────────────────────────────────────


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
    from .diagnostics import bridge_row_violation

    message = bridge_row_violation(trellis, bridge)
    assert message is None, message


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
    from .diagnostics import bridge_side_violations

    messages = bridge_side_violations(
        holes, orientation_preserving=orientation_preserving
    )
    assert not messages, "; ".join(messages)


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
    which still takes part in the singleton pinch. At an ANCHOR the same rule
    reads on ``s(o, b)``, the branch the bridge of ``u(o, b)`` meets there
    (:meth:`Trellis.stable_key_at`), with one identification: the two stable
    branches of the point are one smooth curve through it, so the anchorward
    half of ``s(o, b)`` IS the outward half of ``s(o, 1 - b)``. Nothing opens
    when that branch was never computed (a point seeded on one
    eigendirection).

    Returns:
        The ``Hole.openings`` records ``(intersection_id, which, row,
        branch_key)``, each naming the stable branch it opens.
    """
    openings = []
    for endpoint, iid in zip(("first", "second"), bridge.id):
        ix = trellis.intersection(iid)
        outward = (ix.crossing_sign > 0) == (bridge_side == "left")
        if ix.is_anchor:
            key = trellis.stable_key_at(iid, trellis.unstable_key_of(bridge.id))
            if not outward:
                # s(o, b) and s(o, 1 - b) are one curve through the point: its
                # anchorward half here is the outward half of s(o, 1 - b).
                fp, stability, orbit_index, branch_index = key
                key = (fp, stability, orbit_index, 1 - branch_index)
            if key in trellis.branches:
                openings.append((iid, "outward", bridge_side, key))
            continue
        openings.append((
            iid,
            "outward" if outward else "anchorward",
            row_of_end(trellis, bridge.id, endpoint),
            ix.manifold_b_key,
        ))
    return openings


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


def _is_forward_iterate(pair: PseudoneighborPair) -> bool:
    """True for a forward-iterated pair: it gets no direct hole.

    Any ``iterate > 0`` is the image of a region still attached to the stable
    manifold, on whichever branch of the orbit it lands (author, 2026-10-09,
    reversing D5's ``+1 .. +(k_value - 1)`` exemption). References and
    backward iterates return False.
    """
    return pair.iterate is not None and pair.iterate > 0


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


def _punch_in_bridge(
    trellis: "Trellis",
    bridge: "Bridge",
    *,
    bridge_side: Side,
    iterate: int,
    origin: Optional[tuple[int, int]],
) -> Hole:
    """Punch a propagated hole on the whole containing bridge.

    ``bridge_side`` is the side carried from the orbit (``Hole.bridge_side``);
    the crossing signs turn it into the intervals the hole opens at the
    bridge's defining intersections (:func:`_hole_openings`).
    """
    near_id, far_id = _near_far(
        trellis, bridge.first_intersection, bridge.second_intersection
    )
    return Hole(
        near_intersection_id=near_id,
        bridge_side=bridge_side,
        bounding_ids=(near_id, far_id),
        openings=_hole_openings(trellis, bridge, bridge_side),
        iterate=iterate,
        origin=origin,
    )


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
                "Hole on %s has no openings; it does not partition", hole.bounding_ids
            )
            continue
        for iid, which, row, opened in hole.openings:
            if (row, opened) != (side, branch_key):
                continue
            boundaries[iid] = trellis.intersection(iid).stable_cdist
            (open_outward if which == "outward" else open_anchorward).add(iid)
    return boundaries, open_outward, open_anchorward


#: One boundary set of a partition on one (branch, side): the cdist of every
#: boundary id (None = the synthetic anchor), and the ids whose outward /
#: anchorward interval is open there.
Marks = tuple[dict[Optional[int], float], set[Optional[int]], set[Optional[int]]]


def _build_intervals(
    trellis: "Trellis",
    branch: "TrellisBranch",
    marks: Marks,
) -> list[PartitionInterval]:
    """Assemble the ordered partition intervals from the hole marks.

    ``marks`` is :func:`_hole_marks_on_branch`'s output: the boundary
    intersections and, per boundary, whether a hole opens its outward and/or
    anchorward interval. The branch adds its own two ends -- the outermost
    crossing and the anchor -- and :func:`intervals_from_marks` does the rest.
    """
    boundary_map, open_outward, open_anchorward = marks
    boundaries = dict(boundary_map)
    if len(branch):
        outer_id = branch.ordered_ids(toward_anchor=True)[0]
        boundaries[outer_id] = trellis.intersection(outer_id).stable_cdist

    # The partition starts at the anchor point (cdist 0). When a boundary
    # already sits there (the registered anchor, synthetic at exactly 0.0), it
    # plays that role itself; otherwise a synthetic anchor boundary is added.
    if 0.0 not in boundaries.values():
        boundaries[None] = 0.0

    if len(boundaries) == 1:
        # C1: no crossing past the anchor yet -- one element to the end of the
        # computed branch.
        ((bid, cdist),) = boundaries.items()
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
    return intervals_from_marks((boundaries, open_outward, open_anchorward))


def intervals_from_marks(marks: Marks) -> list[PartitionInterval]:
    """
    The ordered intervals of one branch and side, from its marks.

    Shared by the homotopy partition (:func:`_build_intervals`, after it adds
    the branch's two ends) and its refinement (``PartitionFamily``, which
    rebuilds from a result's own marks plus the cuts). Element ids are stamped
    by the caller.

    Args:
        marks: The boundary cdists by id and the ids opened outward /
            anchorward (:data:`Marks`).

    Returns:
        The intervals in increasing cdist, singletons included.

    Note:
        A boundary point excluded by both neighbouring intervals owns itself
        as a closed singleton: an intersection with hole regions on both
        sides, the outermost intersection when a hole abuts it from the
        interior, or the anchor-artifact intersection when a hole abuts it
        from outside. Each point thereby belongs to exactly one piece; an
        interior boundary always bounds a hole, so at least one of its
        intervals is open there and "][" cannot occur. Two holes flanking one
        bridge (one on each side) emit complementary openings at both defining
        intersections, so this same rule pinches each of them into a singleton
        with open intervals on both sides -- the "singleton bridge" case needs
        no extra code.
    """
    boundaries, open_outward, open_anchorward = marks
    ordered = sorted(boundaries.items(), key=lambda b: b[1])
    intervals = [
        PartitionInterval(
            lo_id, hi_id, lo_c, hi_c,
            closed_lo=lo_id not in open_outward,
            closed_hi=hi_id not in open_anchorward,
        )
        for (lo_id, lo_c), (hi_id, hi_c) in zip(ordered, ordered[1:])
    ]

    with_singletons: list[PartitionInterval] = []
    for position, (bid, cdist) in enumerate(ordered):
        before = intervals[position - 1] if position > 0 else None
        after = intervals[position] if position < len(intervals) else None
        owned = (before is not None and before.closed_hi) or (
            after is not None and after.closed_lo
        )
        if not owned and (before is not None or after is not None):
            with_singletons.append(PartitionInterval(bid, bid, cdist, cdist, True, True))
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
        tol: The slack at the ends: 0.0 (exact) for a registry cdist, the
            relative error of a scaled cdist at a scaled landing end (C15).

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
    never both closed at their shared end (see :func:`intervals_from_marks`),
    so exactly one element owns every point of the branch.

    Args:
        interval: The partition element (interval) to test.
        cdist: Stable canonical distance of the point.
        tol: The slack at the ends, as in :func:`span_contains`: 0.0 (exact)
            for a registry cdist (C15).

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
