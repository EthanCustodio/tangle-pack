from __future__ import annotations

import logging
from typing import Callable, Iterable, Literal, Optional, TYPE_CHECKING, Union

import numpy as np
import matplotlib.pyplot as plt
from numpy.typing import NDArray

from .Pseudoneighbor import forward_unstable_branch_cycle
from .TopologyResults import Hole, PartitionInterval, PseudoneighborPair, StablePartitionResult

if TYPE_CHECKING:
    from .Trellis import Trellis
    from .TrellisBranch import TrellisBranch
    from ..numerics.Bridge import Bridge
    from ..numerics.FixedPoint import FixedPoint
    from ..numerics.Intersection import Intersection, ManifoldKey

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

Side = Literal["left", "right"]

"""
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
* Side (REDONE per the author, 2026-07-09 — supersedes every earlier side
  test, including the k=10 sign calibration and the local-midpoint lobe
  test): left/right is defined by the dynamical direction — standing ON the
  stable manifold looking toward the anchor point, with the STANDARD
  orientation (positive cross(toward-anchor tangent, displacement) = left).
  A hole's side is found by following its bridge to a defining intersection
  and asking which side of the stable manifold the arc approaches that
  intersection from. This is well-defined because a bridge's endpoints are
  consecutive crossings: the arc between them never crosses the (trimmed)
  stable manifold, so the approach side is the side the whole arc — and the
  lobe it bounds — locally hugs. The row is computed per defining
  intersection (each end against its own local stable tangent; on multi-
  branch tangles the two ends may legitimately sit on different branches).
* Hole generation: each reference bridge is mapped backward step by step;
  every backward image lies within some existing bridge, and a hole is punched
  in that containing bridge's region. Termination is the author's "the bridge
  begins to map onto itself" rule, stated k-aware: one backward step moves the
  image to the PREVIOUS unstable branch of the cycle, so a bridge can only
  map onto itself after a multiple of k steps — the recursion stops when the
  containing bridge repeats one already visited at the same branch-cycle
  residue (step mod k). For k = 1 that is exactly the consecutive self-map
  condition; comparing consecutive steps for k > 1 can never fire (the bridges
  live on different branches) and would run the backward orbit far past the
  topologically significant holes. The containing bridge is found by
  canonical-distance bookkeeping (one backward step divides unstable cdists
  by beta = lambda_u^(1/k) and steps the branch cycle back by one). The
  reference hole's coordinates are carried backward
  with the real inverse map: the carried point lies inside the true image
  region, so it decides which side of the unstable arc the hole plots on and
  which side of the stable manifold it counts toward — two orbits landing in
  the same bridge differ exactly there, and BOTH are punched (per the author:
  no cross-orbit deduplication, even on the same bridge and side; only a
  same-orbit same-iterate duplicate of a recorded pair's hole is skipped).
* Interior/exterior: Hole.interior records whether the hole's position lies
  inside the resonance zone (via the in_zone test, else the cdist-vs-cut
  fallback). It is DESCRIPTIVE ONLY — the partition treats every hole on a
  side equally, with no exceptions by iterate, endpoint identity, or zone
  membership ("all holes on manifolds should be treated equal" — author.
  A zone-membership gate was tried and rejected: the deep backward sub-arcs
  hug the zone's own boundary bridge, so the point-in-polygon test split
  otherwise identical holes by the luck of a nudge direction).
* Partition intervals (semantics fixed with the author, July 2026 — this
  SUPERSEDES the earlier "only the outermost beyond-cut hole" rule and the
  piece-containment openness test): EVERY punched hole participates.
  Boundaries are the defining intersections of the holes' bridges
  (Hole.bounding_ids — the pair's own bridge for a direct hole, the
  containing bridge for a propagated one) plus the anchor and the branch's
  outermost intersection. The openness rule is LOCAL to each defining
  intersection and keyed to the hole's side of the BRIDGE ARC: at each
  defining intersection the arc crosses the stable manifold, so each side
  of the arc contains exactly one of the two stable half-intervals there —
  the hole opens the interval on ITS side (Hole.openings). A direct hole
  sits on the side facing the pair's own stable segment, which reduces to
  the inward pair (outward of the near bound, anchorward of the far
  bound). A propagated hole whose region hangs off the removed tail can
  instead open the OUTWARD side of a defining intersection — two holes
  flanking one bridge (the author's "singleton bridge with a hole on
  either side") thereby open complementary sides of BOTH defining
  intersections and pinch each into a closed singleton. In general a
  boundary point excluded by both neighbouring intervals owns itself as a
  degenerate closed singleton [x, x]; a hole marked outward of the
  branch's outermost intersection or anchorward at the anchor-artifact
  contributes its opening to that pinch test even though no interval lies
  there. No point is special — the fixed point and the strong pip pinch
  into singletons by exactly the same rule. Each point belongs to exactly
  one piece ("](" and ")[" occur, "][" cannot). A closed interval
  straddling the resonance-zone cut is split there, with the cut point
  assigned to the interior side; interior/exterior is a REPORTING split
  only and never drops a hole.
* Propagation start: a reference whose pair has no spanning Bridge object
  (after a blast there is no bridge between a parent-manifold crossing and
  a blast-child crossing) punches no direct hole, but its orbit still
  propagates — the backward chain continues from the orbit's most-backward
  directly-punched pair instead (author, 2026-07-09: the un-punchable pair
  itself is fine to skip, its deeper backward holes are not).

Caveats: inversion (k_value == 2*period) follows the same cycle bookkeeping as
StrongPip/Pseudoneighbor but is unvalidated. Containing-bridge lookup compares
unstable cdists, which are only comparable on the same unstable branch; bridges
whose endpoints lack a manifold_a_key (iterated-bridge children) are accepted
on cdist evidence alone. A pseudoneighbor pair with no spanning Bridge object
(after blasting there is no bridge between a parent-manifold crossing and a
blast-child crossing) punches no hole and therefore contributes nothing to the
partition — an open gap in the punching layer, see punch_holes' warning.
"""


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
    """
    wanted = set(pair.as_tuple())
    for bridge in trellis.bridges:
        if {bridge.first_intersection, bridge.second_intersection} == wanted:
            return bridge
    return None


def punch_holes(
    trellis: "Trellis",
    pairs: Optional[Iterable[PseudoneighborPair]] = None,
    *,
    epsilon: float = 0.05,
    in_zone: Optional[Callable[[tuple[float, float]], bool]] = None,
) -> list[Hole]:
    """
    Punch a hole in the region bounded by each pseudoneighbor pair.

    Placement inside the bounded region depends on the pair's role: a
    REFERENCE pair's hole hugs the stable manifold (the pair's chord midpoint,
    nudged ``epsilon`` toward the bridge), while an iterated pair's hole sits
    between the two neighbors on the unstable manifold (the bridge midpoint,
    nudged ``epsilon`` toward the chord). Each hole is classified by side
    (left/right of the oriented stable manifold) and as inside or outside the
    resonance zone, and carries its orbit identity (``origin``/``iterate``).
    The hole is attached to its pair (``pair.hole``).

    Args:
        trellis: The Trellis carrying the pairs, bridges, and cut points.
        pairs: Pairs to punch holes for. Defaults to every pair recorded on the
            trellis (``trellis.pseudoneighbors``).
        epsilon: Inward step off the manifold, as a fraction of the pair's
            chord length.
        in_zone: Optional membership test for the resonance zone (e.g.
            ``ResonanceZone.contains_point``). When given it decides
            ``Hole.interior`` from the hole's actual position — required to
            recognise a hole belonging to a different zone even though its
            stable interval is inside the cut. Without it, the cdist-vs-cut
            fallback is used.

    Returns:
        The punched holes. Pairs with no spanning bridge are skipped with a
        warning.
    """
    if pairs is None:
        pairs = trellis.pseudoneighbors
    pairs = list(pairs)

    cut_map = _cut_cdist_by_branch(trellis)
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
        # Placement anchors on the boundary-arc midpoints; side and openings
        # come from the bridge-approach rule (see _hole_openings).
        stable_mid, _stable_tangent = _stable_arc_midpoint(trellis, near, far)
        arc_point = _nearest_arc_point(bridge, stable_mid)
        if arc_point is None:
            skipped += 1
            continue
        openings, side = _hole_openings(trellis, bridge)

        chord_len = float(np.linalg.norm(far.get_point() - near.get_point()))
        if pair.is_reference:
            coords = _step_into(stable_mid, arc_point, epsilon, chord_len)
        else:
            coords = _step_into(arc_point, stable_mid, epsilon, chord_len)
        coords = (float(coords[0]), float(coords[1]))

        hole = Hole(
            coords=coords,
            near_intersection_id=near_id,
            pair=pair,
            side=side,
            interior=(
                bool(in_zone(coords))
                if in_zone is not None
                else _is_interior(
                    far.stable_cdist, cut_map.get(near.manifold_b_key), trellis
                )
            ),
            bounding_ids=(near_id, far_id),
            openings=openings,
            iterate=pair.iterate,
            origin=pair.origin,
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
    in_zone: Optional[Callable[[tuple[float, float]], bool]] = None,
) -> list[Hole]:
    """
    Punch the holes generated by mapping each reference bridge backward.

    Each reference pair's bridge is mapped backward one step at a time; every
    image lies within some existing bridge, and a hole is punched in that
    containing bridge's region. The recursion stops when the backward orbit
    becomes periodic — for a period-k orbit one backward step moves the image
    to the PREVIOUS unstable branch of the cycle, so a bridge can only map
    onto itself after a multiple of k steps; termination therefore compares
    the containing bridge against those already visited at the same
    branch-cycle residue (step count mod k), which for k = 1 reduces exactly
    to the "image lands in the same bridge" self-map rule. The recursion also
    stops when no containing bridge is found (the image left the computed
    trellis) or after ``max_steps`` (safety guard).

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

    cut_map = _cut_cdist_by_branch(trellis)
    new_holes: list[Hole] = []
    seen: set[tuple[Optional[tuple[int, int]], Optional[int]]] = {
        (pair.origin, pair.iterate)
        for pair in trellis.pseudoneighbors
        if pair.hole is not None
    }
    # Regions already punched, per origin. A keyless containing bridge can be
    # matched from more than one cycle residue (cdist evidence alone), so the
    # same region could otherwise be punched twice for one orbit.
    punched_regions: set[tuple] = {
        (pair.origin, tuple(sorted(pair.hole.bounding_ids)))
        for pair in trellis.pseudoneighbors
        if pair.hole is not None and pair.hole.bounding_ids is not None
    }

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

        fixed_point = _pair_fixed_point(trellis, ref)
        lambda_u = trellis.lambda_u(fixed_point)
        if lambda_u is None:
            logger.warning("No eigenvalue for %s; cannot propagate", ref.as_tuple())
            continue
        cycle = forward_unstable_branch_cycle(fixed_point)
        k = len(cycle)
        beta = lambda_u ** (1.0 / k)

        span, pos = _bridge_unstable_span(trellis, bridge, cycle)
        if pos is None and k > 1:
            # Without the branch identity the containing-bridge lookup would
            # compare cdists across branches, which is meaningless for k > 1.
            logger.warning(
                "Bridge of reference %s has no resolvable unstable branch; "
                "cannot propagate on a period-%d orbit",
                ref.as_tuple(),
                fixed_point.period,
            )
            continue
        carried = np.asarray(start.hole.coords, dtype=np.float64)
        start_iterate = start.iterate or 0

        # Containing bridges already visited, per branch-cycle residue. A
        # repeat at the same residue means the backward orbit is periodic
        # from here on — the k-aware generalization of "the bridge maps onto
        # itself" (see the docstring). Residues are indexed by orbit iterate
        # so a chain started from a fallback pair stays aligned.
        visited: dict[int, set[int]] = {(-start_iterate) % k: {id(bridge)}}

        for step in range(1, max_steps + 1):
            iterate = start_iterate - step
            span = (span[0] / beta, span[1] / beta)
            pos = (pos - 1) % k if pos is not None else None
            carried = _map_backward(trellis, carried)

            containing = _containing_bridge(trellis, span, cycle[pos] if pos is not None else None)
            if containing is None:
                logger.debug(
                    "Backward image of %s left the computed bridges", ref.as_tuple()
                )
                break
            residue = (-iterate) % k
            if id(containing) in visited.get(residue, set()):
                break  # the backward orbit has become periodic — terminate
            visited.setdefault(residue, set()).add(id(containing))

            key = (origin, iterate)
            if key not in seen:
                hole = _punch_in_bridge(
                    trellis, containing, cut_map,
                    iterate=iterate, origin=origin, carried=carried, span=span,
                    in_zone=in_zone,
                )
                if hole is not None:
                    region = (origin, tuple(sorted(hole.bounding_ids)))
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
    singleton (see the module Dev Notes). The branch splits into interior and
    exterior parts at the resonance-zone cut for reporting; the split never
    drops a hole.

    Args:
        trellis: The Trellis carrying the holes to partition by.
        branch_key: Key of the stable branch to partition.
        side: Which side's holes to use.

    Returns:
        The StablePartitionResult for this branch and side.

    Raises:
        ValueError: If branch_key does not name a stable branch of the trellis.
    """
    branch = trellis.branch(branch_key)
    if branch is None or branch.stability != "stable":
        raise ValueError(f"{branch_key} is not a stable branch of this trellis")

    cut_id, cut_cdist = _cuts_by_branch(trellis).get(branch_key, (None, None))
    marks = _hole_marks_on_branch(trellis, branch_key, side)
    intervals = _build_intervals(trellis, branch, marks, cut_cdist, cut_id)
    if cut_cdist is None:
        return StablePartitionResult(
            branch_key=branch_key, side=side,
            interior_intervals=intervals, cut_cdist=None,
        )

    # Same "at or below the cut" slack as _is_interior, so a hole and the
    # partition never disagree about which side of the zone it is on.
    threshold = cut_cdist + trellis.registry.cdist_tol
    return StablePartitionResult(
        branch_key=branch_key,
        side=side,
        interior_intervals=[iv for iv in intervals if iv.hi_cdist <= threshold],
        exterior_intervals=[iv for iv in intervals if iv.hi_cdist > threshold],
        cut_cdist=cut_cdist,
    )


def plot_stable_partition(
    results: Union[StablePartitionResult, Iterable[StablePartitionResult]],
    ax=None,
    **line_kwargs,
):
    """
    Draw stable partitions as number lines (stable cdist on the x-axis).

    Each partition occupies one row; its intervals are horizontal segments
    with a bracket at each end — [ ] closed, ( ) open — every boundary
    labelled with its intersection id, and a circle marking each singleton
    piece [x, x] (a point pinched between two open intervals). Dotted
    vertical lines mark the beginning and end of the manifold.

    Args:
        results: One StablePartitionResult or an iterable of them.
        ax: Optional matplotlib Axes. Defaults to the current axes.
        **line_kwargs: Forwarded to the interval ``plot`` calls (e.g.
            ``linewidth``, ``color``).

    Returns:
        The Axes drawn on, or None if there was nothing to draw.
    """
    if isinstance(results, StablePartitionResult):
        results = [results]
    results = list(results)
    if not results:
        logger.info("No partitions to plot; call partition_stable_manifold() first.")
        return None

    target = ax if ax is not None else plt.gca()
    line_kwargs.setdefault("linewidth", 2)

    labels = []
    extent_lo, extent_hi = np.inf, -np.inf
    for row, result in enumerate(results):
        color = "tab:blue" if result.side == "left" else "tab:red"
        labels.append(
            f"{result.side} (p{result.branch_key[0].period}, "
            f"orbit {result.branch_key[2]})"
        )
        labelled: set[int] = set()
        for interval in result.intervals:
            extent_lo = min(extent_lo, interval.lo_cdist)
            extent_hi = max(extent_hi, interval.hi_cdist)
            target.plot(
                [interval.lo_cdist, interval.hi_cdist], [row, row],
                color=color, solid_capstyle="butt", **line_kwargs,
            )
            # Bracket per interval END, nudged inward so the two intervals
            # meeting at one intersection point stay individually readable:
            # e.g. "](" = closed on the anchor side, open on the outward side.
            if interval.lo_cdist == interval.hi_cdist:
                # A singleton piece [x, x] — a point pinched between two open
                # intervals. Brackets are illegible at zero width; a FILLED
                # circle marks the closed point, with the neighbours' open
                # parentheses either side.
                target.plot(
                    interval.lo_cdist, row, marker="o", markersize=9,
                    markerfacecolor=color, markeredgecolor=color,
                    zorder=7,
                )
            else:
                target.annotate(
                    "[" if interval.closed_lo else "(",
                    (interval.lo_cdist, row), textcoords="offset points",
                    xytext=(1, 0), ha="left", va="center",
                    fontsize=13, fontweight="bold", color=color, zorder=6,
                )
                target.annotate(
                    "]" if interval.closed_hi else ")",
                    (interval.hi_cdist, row), textcoords="offset points",
                    xytext=(-1, 0), ha="right", va="center",
                    fontsize=13, fontweight="bold", color=color, zorder=6,
                )
            # Label each boundary with the intersection's id in the tangle,
            # staggered over two text rows so tightly packed boundaries near
            # the anchor stay legible.
            for boundary_id, cdist in (
                (interval.lo_id, interval.lo_cdist),
                (interval.hi_id, interval.hi_cdist),
            ):
                if boundary_id is None or boundary_id in labelled:
                    continue
                drop = -14 if len(labelled) % 2 == 0 else -24
                labelled.add(boundary_id)
                target.annotate(
                    str(boundary_id), (cdist, row), textcoords="offset points",
                    xytext=(0, drop), ha="center", va="top", fontsize=8,
                    color="black",
                )

    # Beginning and end of the (trimmed) stable manifold.
    for cdist, name, align, nudge in (
        (extent_lo, "anchor", "left", 3),
        (extent_hi, "manifold end", "right", -3),
    ):
        target.axvline(cdist, color="gray", linestyle=":", alpha=0.6, zorder=0)
        target.annotate(
            name, (cdist, 1.0), xycoords=("data", "axes fraction"),
            textcoords="offset points", xytext=(nudge, -10), ha=align,
            fontsize=8, color="gray",
        )

    target.set_yticks(range(len(results)))
    target.set_yticklabels(labels)
    target.set_ylim(-0.6, len(results) - 0.4)
    target.set_xlabel("stable canonical distance")
    target.set_title(
        "Stable manifold partition — [ ] closed, ( ) open, ○ singleton; "
        "numbers = ids"
    )
    return target


# ── geometric / bookkeeping helpers ─────────────────────────────────────────


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


def _cross_sign(a: NDArray[np.float64], b: NDArray[np.float64]) -> Optional[float]:
    """Sign of the 2D cross product a × b (None when collinear)."""
    cross = float(a[0] * b[1] - a[1] * b[0])
    if cross == 0.0:
        return None
    return float(np.sign(cross))


def _stable_frame(
    trellis: "Trellis", ix: "Intersection"
) -> tuple[Optional[NDArray[np.float64]], Optional[NDArray[np.float64]]]:
    """(anchorward, outward) unit directions of the stable branch at ``ix``.

    Walked from the live manifold nodes around the intersection's stable
    cdist. Either may be None: anchorward at the anchor artifact (no nodes
    below), outward at the trimmed branch end (no nodes above), both when
    the manifold nodes are unavailable.
    """
    key = ix.manifold_b_key
    manifold = trellis.manifolds.get(key) if key is not None else None
    if manifold is None:
        return None, None
    tol = trellis.registry.cdist_tol
    here = np.asarray(ix.get_point(), dtype=np.float64)
    below = above = None
    below_c, above_c = -np.inf, np.inf
    for node in manifold.get_point_array(return_nodes=True):
        point = np.asarray(node.get_point(), dtype=np.float64)
        if not np.linalg.norm(point - here) > 0.0:
            continue
        if node.cdist < ix.stable_cdist - tol and node.cdist > below_c:
            below, below_c = point, node.cdist
        elif node.cdist > ix.stable_cdist + tol and node.cdist < above_c:
            above, above_c = point, node.cdist
    anchorward = (below - here) / np.linalg.norm(below - here) if below is not None else None
    outward = (above - here) / np.linalg.norm(above - here) if above is not None else None
    return anchorward, outward


def _bridge_end_geometry(
    poly: NDArray[np.float64], q: NDArray[np.float64]
) -> Optional[tuple[NDArray[np.float64], NDArray[np.float64]]]:
    """Approach displacement and local tangent at the arc end nearest ``q``.

    The displacement points from the intersection INTO the bridge (the root
    and tail nodes sit just past the crossings, so the walk starts one node
    in and takes the second clearly-distinct node when one exists); the
    tangent is the local arc direction in root→tail order, matching the
    orientation `_arc_side_of` signs against.
    """
    if len(poly) < 2:
        return None
    chord = float(np.linalg.norm(poly[-1] - poly[0]))
    eps = 1e-9 * (1.0 + chord)
    if np.linalg.norm(poly[0] - q) <= np.linalg.norm(poly[-1] - q):
        indices = range(1, len(poly))
    else:
        indices = range(len(poly) - 2, -1, -1)
    picked = None
    found = 0
    for i in indices:
        if float(np.linalg.norm(poly[i] - q)) > eps:
            picked = i
            found += 1
            if found == 2:
                break
    if picked is None:
        return None
    tangent = poly[min(picked + 1, len(poly) - 1)] - poly[max(picked - 1, 0)]
    norm = float(np.linalg.norm(tangent))
    if norm == 0.0:
        return None
    return poly[picked] - q, tangent / norm


def _arc_side_of(
    poly: NDArray[np.float64], point: NDArray[np.float64]
) -> Optional[float]:
    """Sign of ``point``'s side of the arc polyline (root→tail orientation)."""
    if len(poly) < 3:
        return None
    i = int(np.argmin(np.linalg.norm(poly - point, axis=1)))
    i = min(max(i, 1), len(poly) - 2)
    return _cross_sign(poly[i + 1] - poly[i - 1], point - poly[i])


def _hole_openings(
    trellis: "Trellis",
    bridge: "Bridge",
    carried: Optional[NDArray[np.float64]] = None,
) -> tuple[list[tuple[int, str, Side]], Optional[Side]]:
    """The openings and display side of a hole attached to ``bridge``.

    Implements the author's rule (2026-07-09): follow the bridge arc to each
    of its two defining intersections. The ROW (left/right) is the side of
    the stable manifold the arc approaches that intersection from — standing
    on the manifold looking toward the anchor, standard orientation. The
    interval opened there (anchorward or outward of the intersection) is the
    stable half-interval lying on the HOLE'S side of the arc: the arc crosses
    the stable manifold at the intersection, so each side of the arc contains
    exactly one of the two halves. Without ``carried`` (a direct hole, which
    sits on the side of the arc facing the pair's own stable segment) the
    openings reduce to the inward pair — outward of the near bound,
    anchorward of the far bound.

    At a trimmed branch end whose outward half is removed, a hole facing that
    half still records an ``"outward"`` opening (it participates in the
    singleton pinch test); at the anchor artifact a hole facing the far half
    faces the OTHER stable branch and opens nothing here.

    Returns:
        (openings, side): the ``Hole.openings`` records and the row at the
        far (outermost) defining intersection — a genuine transversal
        crossing even when the near one is the anchor artifact.
    """
    near_id, far_id = _near_far(
        trellis, bridge.first_intersection, bridge.second_intersection
    )
    points = bridge.get_point_array()
    if points is None or len(points) < 2:
        return [], None
    poly = np.asarray(points, dtype=np.float64)

    arc_side = None
    if carried is not None:
        arc_side = _arc_side_of(poly, np.asarray(carried, dtype=np.float64))
        if arc_side is None:
            logger.warning(
                "Carried point of a propagated hole sits on the bridge arc; "
                "falling back to inward openings"
            )

    openings: list[tuple[int, str, Side]] = []
    display_side: Optional[Side] = None
    for iid in (near_id, far_id):
        ix = trellis.intersection(iid)
        q = np.asarray(ix.get_point(), dtype=np.float64)
        end = _bridge_end_geometry(poly, q)
        if end is None:
            continue
        disp, t_end = end
        anchorward, outward = _stable_frame(trellis, ix)
        look = anchorward if anchorward is not None else (
            -outward if outward is not None else None
        )
        if look is None:
            continue  # no manifold data at this end
        row = _side_of(look, disp)
        if row is None:
            continue
        if arc_side is None:
            which = "outward" if iid == near_id else "anchorward"
        else:
            side_anchorward = (
                _cross_sign(t_end, anchorward) if anchorward is not None else None
            )
            side_outward = (
                _cross_sign(t_end, outward) if outward is not None else None
            )
            if side_anchorward == arc_side:
                which = "anchorward"
            elif side_outward == arc_side:
                which = "outward"
            elif side_outward is None and side_anchorward is not None:
                # Trimmed branch end: the hole's side holds the removed tail.
                which = "outward"
            else:
                # Anchor artifact: the hole's side holds the OTHER branch.
                logger.debug(
                    "Hole faces the opposite stable branch at intersection %d",
                    iid,
                )
                continue
        openings.append((iid, which, row))
        if iid == far_id or display_side is None:
            display_side = row
    return openings, display_side


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
    midpoint and chord direction when the manifold nodes are unavailable.
    """
    key = near.manifold_b_key
    manifold = trellis.manifolds.get(key) if key is not None else None
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
    """Order two intersection ids as (toward-anchor, outward) by stable cdist."""
    if trellis.intersection(id_a).stable_cdist <= trellis.intersection(id_b).stable_cdist:
        return id_a, id_b
    return id_b, id_a


def _bridge_midpoint(bridge: "Bridge") -> Optional[NDArray[np.float64]]:
    """The bridge's geometric middle node, or None for an empty bridge."""
    points = bridge.get_point_array()
    if points is None or len(points) == 0:
        return None
    return np.asarray(points[len(points) // 2], dtype=np.float64)


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


def _cut_cdist_by_branch(trellis: "Trellis") -> dict["ManifoldKey", float]:
    """Map each stable branch key to its resonance-zone cut stable cdist."""
    return {key: cdist for key, (_iid, cdist) in _cuts_by_branch(trellis).items()}


def _cuts_by_branch(trellis: "Trellis") -> dict["ManifoldKey", tuple[int, float]]:
    """Map each stable branch key to its cut point's (intersection id, cdist)."""
    cuts: dict["ManifoldKey", tuple[int, float]] = {}
    for iid in trellis.strong_pip_cut_points():
        ix = trellis.intersection(iid)
        if ix.manifold_b_key is not None:
            cuts[ix.manifold_b_key] = (iid, ix.stable_cdist)
    return cuts


def _is_interior(outer_cdist: float, cut_cdist: Optional[float], trellis: "Trellis") -> Optional[bool]:
    """Interior iff the region's outer stable cdist is at or below the zone cut."""
    if cut_cdist is None:
        return None
    return outer_cdist <= cut_cdist + trellis.registry.cdist_tol


def _pair_fixed_point(trellis: "Trellis", pair: PseudoneighborPair) -> "FixedPoint":
    """The fixed point owning the pair's stable branch."""
    key = pair.branch_key or trellis.intersection(pair.intersection_a).manifold_b_key
    return key[0]


def _bridge_unstable_span(
    trellis: "Trellis", bridge: "Bridge", cycle: list["ManifoldKey"]
) -> tuple[tuple[float, float], Optional[int]]:
    """A bridge's (lo, hi) unstable-cdist span and its branch-cycle position."""
    a = trellis.intersection(bridge.first_intersection)
    b = trellis.intersection(bridge.second_intersection)
    lo, hi = sorted((a.unstable_cdist, b.unstable_cdist))
    pos = None
    for ix in (a, b):
        if ix.manifold_a_key is not None and ix.manifold_a_key in cycle:
            pos = cycle.index(ix.manifold_a_key)
            break
    return (lo, hi), pos


def _containing_bridge(
    trellis: "Trellis",
    span: tuple[float, float],
    branch_key: Optional["ManifoldKey"],
    rtol: float = 1e-3,
) -> Optional["Bridge"]:
    """
    The bridge whose unstable-cdist span contains ``span``.

    Prefers bridges whose endpoints are known to lie on ``branch_key``; a
    bridge with unknown endpoint branches (iterated children) is accepted on
    cdist evidence alone (see Dev Notes). Among multiple containers the
    tightest (smallest) span wins.
    """
    lo, hi = span
    slack = rtol * (hi - lo)
    best = None
    best_width = np.inf
    for bridge in trellis.bridges:
        if bridge.first_intersection is None or bridge.second_intersection is None:
            continue
        a = trellis.intersection(bridge.first_intersection)
        b = trellis.intersection(bridge.second_intersection)
        keys = {a.manifold_a_key, b.manifold_a_key} - {None}
        if branch_key is not None and keys and branch_key not in keys:
            continue
        b_lo, b_hi = sorted((a.unstable_cdist, b.unstable_cdist))
        if b_lo - slack <= lo and hi <= b_hi + slack:
            if (b_hi - b_lo) < best_width:
                best, best_width = bridge, b_hi - b_lo
    return best


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
    cut_map: dict["ManifoldKey", float],
    *,
    iterate: Optional[int] = None,
    origin: Optional[tuple[int, int]] = None,
    carried: Optional[NDArray[np.float64]] = None,
    span: Optional[tuple[float, float]] = None,
    in_zone: Optional[Callable[[tuple[float, float]], bool]] = None,
) -> Optional[Hole]:
    """Punch a propagated hole between the two backward-imaged neighbors.

    ``span`` is the backward image's unstable-cdist extent inside the
    containing bridge — the two imaged pseudoneighbor points bound it, and the
    hole plots at the middle of that sub-arc (two orbits landing in the same
    bridge occupy different sub-arcs). ``carried`` is the reference hole's
    coordinates mapped backward alongside — a point inside the true image
    region — and decides which side of the unstable arc the hole sits on,
    which in turn decides the intervals it opens at the containing bridge's
    defining intersections (see :func:`_hole_openings`). Without it the
    containing bridge's own midpoint and chord are the fallback and the
    openings default to the inward pair.
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

    openings, side = _hole_openings(trellis, bridge, carried=carried)
    coords = (float(coords[0]), float(coords[1]))
    return Hole(
        coords=coords,
        near_intersection_id=near_id,
        side=side,
        interior=(
            bool(in_zone(coords))
            if in_zone is not None
            else _is_interior(
                far.stable_cdist, cut_map.get(near.manifold_b_key), trellis
            )
        ),
        bounding_ids=(near_id, far_id),
        openings=openings,
        iterate=iterate,
        origin=origin,
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


def _openings_of(
    trellis: "Trellis", hole: Hole
) -> list[tuple[int, str, Side]]:
    """A hole's opening records, deriving the inward default when absent.

    Fabricated holes (tests, external callers) may carry only ``bounding_ids``
    and ``side``; they open the region between their bounds on that side —
    outward of the near bound, anchorward of the far bound.
    """
    if hole.openings is not None:
        return hole.openings
    if hole.bounding_ids is None or hole.side is None:
        return []
    near, far = _near_far(trellis, *hole.bounding_ids)
    return [(near, "outward", hole.side), (far, "anchorward", hole.side)]


def _hole_marks_on_branch(
    trellis: "Trellis",
    branch_key: "ManifoldKey",
    side: Side,
) -> tuple[dict[int, float], set[int], set[int]]:
    """
    Boundary intersections and open-interval marks on one branch and row.

    Every punched hole participates — a hole is a hole, with no exceptions by
    iterate, endpoint identity, or zone membership (``Hole.interior`` is
    descriptive only). Each opening record lands on the stable branch of its
    own intersection, so a bridge whose two defining intersections sit on
    different branches of a periodic orbit contributes to each one.

    Returns:
        (boundaries, open_outward, open_anchorward): boundary intersection
        ids mapped to their stable cdists, and the ids whose outward /
        anchorward interval a hole opens.
    """
    boundaries: dict[int, float] = {}
    open_outward: set[int] = set()
    open_anchorward: set[int] = set()
    for hole in trellis.holes:
        for iid, which, row in _openings_of(trellis, hole):
            if row != side:
                continue
            ix = trellis.intersection(iid)
            if ix.manifold_b_key != branch_key:
                continue
            boundaries[iid] = ix.stable_cdist
            (open_outward if which == "outward" else open_anchorward).add(iid)
    return boundaries, open_outward, open_anchorward


def _build_intervals(
    trellis: "Trellis",
    branch: "TrellisBranch",
    marks: tuple[dict[int, float], set[int], set[int]],
    cut_cdist: Optional[float],
    cut_id: Optional[int] = None,
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
    # already sits there (the anchor-artifact intersection), it plays that role
    # itself; otherwise a synthetic anchor boundary is added.
    tol = trellis.registry.cdist_tol
    if not any(cdist <= tol for _, cdist in boundaries):
        boundaries.append((None, 0.0))

    # Sort and deduplicate by id (the anchor and branch end keep their slots).
    unique: dict = {}
    for bid, cdist in boundaries:
        unique[bid if bid is not None else ("anchor", cdist)] = (bid, cdist)
    ordered = sorted(unique.values(), key=lambda b: b[1])

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
    # cannot occur.
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
    intervals = with_singletons

    if cut_cdist is not None:
        intervals = _split_at_cut(intervals, cut_cdist, cut_id)
    return intervals


def _split_at_cut(
    intervals: list[PartitionInterval],
    cut_cdist: float,
    cut_id: Optional[int] = None,
) -> list[PartitionInterval]:
    """Split a non-hole interval straddling the zone cut; the cut point is interior.

    ``cut_id`` is the cut intersection's registry id (the strong pip or its
    iterate on this branch), labelling the new shared boundary.
    """
    out: list[PartitionInterval] = []
    for iv in intervals:
        is_closed = iv.closed_lo and iv.closed_hi
        if is_closed and iv.lo_cdist < cut_cdist < iv.hi_cdist:
            out.append(
                PartitionInterval(iv.lo_id, cut_id, iv.lo_cdist, cut_cdist, iv.closed_lo, True)
            )
            out.append(
                PartitionInterval(cut_id, iv.hi_id, cut_cdist, iv.hi_cdist, False, iv.closed_hi)
            )
        else:
            out.append(iv)
    return out
