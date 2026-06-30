from __future__ import annotations

import logging
import math
from typing import Iterable, Optional, TYPE_CHECKING

from ..numerics.Intersection import ManifoldKey
from .TopologyResults import StrongPipResult

if TYPE_CHECKING:
    from .Trellis import Trellis
    from ..numerics.FixedPoint import FixedPoint

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

"""
Dev Notes — Is-Strong-Pip (Strong_Pip_Algorithm.pdf)

q0 is a strong pip iff, once every trellis intersection is mapped back onto the
open stable arc W^S(z'_0, q0), none of them has an UNSTABLE canonical distance less
than q0's. The unstable cdist is the disqualifier we compare; the stable cdist only
tells us which mapped points actually land on that arc.

Mapping back (the map-free shortcut):
  * Each intersection lives on a stable branch anchored at its periodic point. The
    candidate q0's branch B = the stable branch containing q0. To bring intersection
    r onto B we map it back i_r steps, where i_r = minimal backward steps in the
    forward branch cycle that carry r's stable branch onto B.
  * Canonical distance scales geometrically under the map, so we never need the real
    inverse map: one forward step of M scales the unstable cdist by lambda_u^(1/k_value)
    and the stable cdist by lambda_u^(-1/k_value) (k_value forward steps = one full
    return to the same branch, scaling by lambda_u^(±1), the full-cycle eigenvalue).
    The PDF writes this for the stable cdist purely as an *example* of the map-free
    scaling — the same holds for the unstable cdist independently; both are used.
    The per-step factor was confirmed empirically against the period-3 Henon trellis.
  * "Get it onto W^S(z'_0, q0)" (algorithm step 2) is a membership test: r counts
    only if its mapped STABLE cdist s_r·lambda_u^(i_r/k) lands on the open arc, i.e.
    is < s0. This is also what fixes the number of backward maps (the lift): mapping
    the minimal i_r steps always shrinks the unstable cdist, so without the stable
    arc test every candidate is trivially disqualified and nothing is a strong pip
    (verified empirically). For the survivors, the mapped UNSTABLE cdist
    u_r·lambda_u^(-i_r/k) is compared against q0's u0 (algorithm step 3).

A disqualifier must land strictly inside the open box (0, s0) x (0, u0). q0's own
orbit maps back onto q0 itself (a collision on BOTH canonical distances), so without
guarding against it q0's iterates falsely block it. This is exactly the period-3
failure that period-1 never hits: with k_value=1 the mapping is the identity (i_r is
always 0) and q0's iterates land outside the box, but for k_value>1 M(q0) maps back
onto q0 to within canonical-distance noise. We skip a point only when it collides
with q0 on BOTH cdists — a 2D coincidence that means the same point. We deliberately
do NOT use the invariant action product (stable_cdist * unstable_cdist) for this:
the map is area-preserving so the product is preserved along one chain, but two
*different* chains generally have *different* products, and two products being equal
does not make two points the same chain — so comparing the product alone would
wrongly merge distinct chains. Comparing each cdist individually is the safe test.

The fixed point itself is reported by the detector as an intersection at canonical
distance 0 on both manifolds; it is the anchor z'_0, not a transverse crossing, so
it is excluded from N (both as a candidate and as a disqualifier).

Same-periodic-point restriction: a strong pip is a property of a single periodic
point's homoclinic tangle, so only crossings between two manifolds of q0's own
fixed point may disqualify it. A nested session detects heteroclinic crossings too
(e.g. a period-1 outer tangle meeting a period-3 inner tangle); those are real and
kept for other purposes, but they must not enter N for strong-pip classification.
The disqualifier loop therefore skips any r whose stable branch is not on q0's
fixed point (it is absent from the branch cycle) or whose unstable branch is known
to be on a different fixed point. The whole orbit of a period-p point counts as one
"periodic point" here — the restriction is on the FixedPoint object, not on the
orbit index or branch index. The unstable key (manifold_a_key) is legitimately None
on points born from iterated bridges; those are kept (the stable key is the reliable
discriminator), so only a *known* foreign unstable key is rejected.

Branch cycle / inversion: forward_stable_branch_cycle() is validated for the
non-inversion case (k_value == period). The inversion ordering (k_value == 2*period,
branch_index flips after each full orbit) is implemented from first principles but
has not yet been validated against a computed inversion trellis — revisit when one
is available.
"""


def forward_stable_branch_cycle(fixed_point: "FixedPoint") -> list[ManifoldKey]:
    """
    Return the fixed point's stable branches in forward-iteration (M) order.

    Element j maps to element j+1 (mod len) under one application of M. The cycle
    has length ``fixed_point.k_value`` — one entry per stable branch (period
    branches without inversion, 2*period with inversion).

    The orbit cycling under M is taken from ``get_iterable_array("unstable")``
    (the M-forward orbit order, along which the unstable manifold grows); the
    stable manifold of z_i maps under M to the stable manifold of M(z_i), so the
    stable branches cycle in the same orbit order. With inversion, the branch
    index flips after each full pass around the orbit.

    Args:
        fixed_point: The fixed point whose stable branches to order.

    Returns:
        List of manifold keys (fixed_point, "stable", orbit_index, branch_index)
        in forward-iteration order.
    """
    orbit_order = fixed_point.get_iterable_array("unstable")
    branch_indices = fixed_point.get_branch_array()
    cycle: list[ManifoldKey] = []
    for branch_index in branch_indices:
        for orbit_index in orbit_order:
            cycle.append((fixed_point, "stable", orbit_index, branch_index))
    return cycle


def _branch_position_map(
    fixed_point: "FixedPoint",
) -> tuple[dict[ManifoldKey, int], int]:
    """Return (position-of-each-stable-branch, cycle length k_value)."""
    cycle = forward_stable_branch_cycle(fixed_point)
    return {key: i for i, key in enumerate(cycle)}, len(cycle)


def is_strong_pip(
    trellis: "Trellis",
    intersection_id: int,
    *,
    tol: Optional[float] = None,
    collision_rtol: float = 1e-2,
    _cache: Optional[dict] = None,
) -> StrongPipResult:
    """
    Classify whether a single intersection is a strong pip.

    Implements Is-Strong-Pip(q0, T): every trellis intersection is mapped back onto
    q0's stable branch using the canonical-distance scaling law, and q0 is a strong
    pip iff none of them lands strictly inside the open box (0, s0) x (0, u0) — on
    the stable arc W^S(z'_0, q0) and inside the unstable arc.

    q0's own orbit maps back onto q0 itself (a collision on both canonical
    distances), so without guarding against it q0's iterates falsely disqualify q0 —
    exactly the failure seen for period > 1, where M(q0) maps back onto q0. The guard
    is a direct cdist collision check (match on BOTH cdists), which only ever skips
    q0's genuine iterates; a distinct iterate chain that merely happens to be nearby
    (or shares the same invariant action product) is not mistaken for q0.

    Args:
        trellis: The Trellis to read intersections and eigen/orbit data from.
        intersection_id: Registry ID of the candidate q0.
        tol: Absolute canonical-distance slack. Defaults to the registry's
            ``cdist_tol``.
        collision_rtol: Relative slack for the cdist collision test that identifies
            q0's own orbit. A point is treated as q0 (and skipped) only if its mapped
            stable AND unstable cdists are both within ``collision_rtol`` of q0's.
            Defaults to 1e-2, which absorbs canonical-distance scaling noise while
            still distinguishing distinct intersections.
        _cache: Internal per-fixed-point cache used by classify_strong_pips to
            avoid rebuilding the branch cycle for every candidate.

    Returns:
        A StrongPipResult. When not a strong pip, ``blocking_intersection_id`` is
        the intersection whose mapped-back unstable cdist is smallest (the innermost
        disqualifier).

    Raises:
        ValueError: If q0 has no stable branch, its fixed point has no eigenvalue,
            or q0's branch is missing from the computed branch cycle.

    Note:
        If q0 is itself the anchoring fixed point (canonical distance ~0 on both
        manifolds), it has no interior and this returns True vacuously.
        classify_strong_pips() skips that artifact; a direct call does not.
    """
    q0 = trellis.intersection(intersection_id)
    branch_key = q0.manifold_b_key
    if branch_key is None:
        raise ValueError(
            f"Intersection {intersection_id} has no stable branch (manifold_b_key); "
            "cannot classify as a strong pip."
        )

    fixed_point = branch_key[0]
    s0 = q0.stable_cdist
    u0 = q0.unstable_cdist
    lambda_u = trellis.lambda_u(fixed_point)
    if lambda_u is None:
        raise ValueError(
            f"Fixed point of intersection {intersection_id} has no unstable "
            "eigenvalue; cannot scale canonical distances."
        )

    if _cache is not None and fixed_point in _cache:
        pos_map, k = _cache[fixed_point]
    else:
        pos_map, k = _branch_position_map(fixed_point)
        if _cache is not None:
            _cache[fixed_point] = (pos_map, k)

    pos_B = pos_map.get(branch_key)
    if pos_B is None:
        raise ValueError(
            f"Stable branch {branch_key} of intersection {intersection_id} is not "
            "in the forward branch cycle; the trellis may be inconsistent."
        )

    if tol is None:
        tol = trellis.registry.cdist_tol

    # Map every intersection back onto q0's stable branch by the minimal i_r
    # backward steps in the branch cycle, scaling stable cdist by lambda_u^(i_r/k)
    # and unstable cdist by lambda_u^(-i_r/k) (both map-free). A disqualifier is one
    # that lands strictly inside the open box (0, s0) x (0, u0): on the stable arc
    # W^S(z'_0, q0) (stable_rep < s0) and inside the unstable arc (unstable_rep < u0).
    # Among those, the innermost (smallest mapped unstable cdist) is reported.
    best_rep = math.inf
    best_id: Optional[int] = None
    for r_id, r in trellis.registry:
        if r_id == intersection_id:
            continue
        if abs(r.stable_cdist) <= tol and abs(r.unstable_cdist) <= tol:
            continue  # the anchoring fixed point, not a transverse intersection
        r_branch = r.manifold_b_key
        if r_branch is None:
            continue
        # A strong pip is a property of one periodic point's own tangle, so only
        # crossings between two manifolds of THIS fixed point may disqualify q0.
        # pos_r below already confirms r's stable branch belongs to fixed_point;
        # here we also reject the case where r's unstable branch is known to
        # belong to a *different* fixed point — a heteroclinic crossing (e.g. a
        # period-1 outer tangle meeting this period-3 inner tangle). Those
        # crossings are real and kept elsewhere, but must not be used to classify
        # a strong pip. A None unstable key is an iterated-bridge homoclinic point
        # on this tangle (manifold_b_key is the reliable discriminator throughout
        # the code), so it is kept.
        if r.manifold_a_key is not None and r.manifold_a_key[0] is not fixed_point:
            continue
        pos_r = pos_map.get(r_branch)
        if pos_r is None:
            # Stable branch of a different fixed point — cannot map onto B.
            continue
        backward_steps = (pos_r - pos_B) % k
        stable_rep = r.stable_cdist * (lambda_u ** (backward_steps / k))
        unstable_rep = r.unstable_cdist * (lambda_u ** (-backward_steps / k))

        # q0's own orbit maps back onto q0 itself — a collision on BOTH cdists. That
        # is the same point, not a disqualifier, so skip it. Requiring a collision on
        # both cdists (a 2D coincidence) means only q0's genuine iterates are skipped;
        # a distinct chain that merely shares q0's stable OR unstable cdist — or an
        # equal "action" product — is not mistaken for q0.
        collides_s = abs(stable_rep - s0) <= collision_rtol * s0 + tol
        collides_u = abs(unstable_rep - u0) <= collision_rtol * u0 + tol
        if collides_s and collides_u:
            continue

        if stable_rep < s0 and unstable_rep < u0 and unstable_rep < best_rep:
            best_rep = unstable_rep
            best_id = r_id

    if best_id is not None:
        logger.debug(
            "Intersection %s is not a strong pip: %s maps into the box at u=%.6g < %.6g",
            intersection_id,
            best_id,
            best_rep,
            u0,
        )
        return StrongPipResult(intersection_id, False, blocking_intersection_id=best_id)

    return StrongPipResult(intersection_id, True)


def classify_strong_pips(
    trellis: "Trellis",
    intersection_ids: Optional[Iterable[int]] = None,
    *,
    tol: Optional[float] = None,
    collision_rtol: float = 1e-2,
) -> dict[int, StrongPipResult]:
    """
    Classify many intersections, sharing the branch-cycle cache (pure query).

    Returns a result per intersection without mutating the trellis. The qualifying
    results (``is_strong_pip is True``) are the strong-pip *candidates*; selecting
    the single actual strong pip from them is the caller's job — see
    Trellis.classify_strong_pips() / set_strong_pip().

    Args:
        trellis: The Trellis to classify within.
        intersection_ids: Candidates to classify. Defaults to every intersection
            in the trellis.
        tol: Canonical-distance slack passed through to is_strong_pip.
        collision_rtol: Relative collision slack passed through to is_strong_pip.

    Returns:
        Mapping from intersection ID to its StrongPipResult.
    """
    if intersection_ids is None:
        intersection_ids = trellis.registry.all_ids()

    anchor_tol = tol if tol is not None else trellis.registry.cdist_tol

    cache: dict = {}
    results: dict[int, StrongPipResult] = {}
    for iid in intersection_ids:
        ix = trellis.intersection(iid)
        # The fixed point is reported by the detector as an intersection at
        # canonical distance 0 on both manifolds. It is the anchor z'_0, not a
        # transverse crossing, so its open interval W^S(z'_0, q0) is empty and
        # classifying it as a pip is meaningless — skip it. A genuine near-anchor
        # crossing has a large canonical distance on one manifold and is kept.
        if abs(ix.stable_cdist) <= anchor_tol and abs(ix.unstable_cdist) <= anchor_tol:
            logger.debug("Skipping anchor/fixed-point intersection %s", iid)
            continue
        results[iid] = is_strong_pip(
            trellis, iid, tol=tol, collision_rtol=collision_rtol, _cache=cache
        )

    logger.debug(
        "Classified %d intersection(s): %d candidate(s)",
        len(results),
        sum(1 for r in results.values() if r.is_strong_pip),
    )
    return results
