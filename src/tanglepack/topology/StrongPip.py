"""
Is-Strong-Pip: choosing the strong primary intersection point of a trellis.

Tests a candidate crossing against the definition (no trellis intersection
mapped back onto its stable arc undercuts it in unstable canonical distance)
and picks the strong pip a fixed point's topological algorithms start from.

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
    inverse map: one forward step of M multiplies the unstable cdist by
    beta = FixedPoint.per_step_beta("unstable") and divides the stable cdist by it.
    beta is |lambda_u|^(1/period) -- the period-th root, because the stored
    eigenvalue is that of the full-cycle Jacobian DM^period. Without inversion
    (k_value == period) a branch return is one orbit and costs lambda_u; with
    inversion it takes k_value = 2*period steps and costs lambda_u squared.
    The PDF writes this for the stable cdist purely as an *example* of the map-free
    scaling — the same holds for the unstable cdist independently; both are used.
    The per-step factor was confirmed empirically against the period-3 Henon trellis.
  * "Get it onto W^S(z'_0, q0)" (algorithm step 2) is a membership test: r counts
    only if its mapped STABLE cdist s_r·beta^(i_r) lands on the open arc, i.e.
    is < s0. This is also what fixes the number of backward maps (the lift): mapping
    the minimal i_r steps always shrinks the unstable cdist, so without the stable
    arc test every candidate is trivially disqualified and nothing is a strong pip
    (verified empirically). For the survivors, the mapped UNSTABLE cdist
    u_r·beta^(-i_r) is compared against q0's u0 (algorithm step 3).

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

Same-pair restriction: a strong pip is a property of ONE pair of manifolds, the
unstable manifold of q0's unstable fixed point and the stable branch cycle of its
stable fixed point, so only crossings of that pair may disqualify q0 (r's unstable
branch on q0's unstable fixed point, r's stable branch in q0's stable cycle). For
a homoclinic q0 that is its own tangle, as before; for a heteroclinic q0 in
``W^u(A) x W^s(B)`` (2026-11-10) it is the other ``W^u(A) x W^s(B)`` crossings,
mapped back in B's stable cycle, while B's homoclinic crossings are kept out. The
whole orbit of a period-p point counts as one "periodic point" here -- the
restriction is on the FixedPoint objects, not on the orbit or branch index. A
heteroclinic q0 needs A and B of the same period and the same ``k_value``
(:func:`same_branch_cycle`; ``ValueError`` otherwise): the pair's branch cycles
then advance together, and the crossings of one pair are
carried onto each other by the map. :func:`classify_strong_pips` skips such a
crossing (it is no candidate); a direct call raises. Every detected crossing
carries both keys (the Tangle asserts it), blast-born ones included.

Open question (not implemented): a ring of heteroclinic pips bounds one zone whose
boundary also runs along the OTHER pairs' arcs (``W^u[A, q0]`` and ``W^s[A, q0']``,
whose crossings are homoclinic). A simple closed boundary would also need the
ring's pips tested against those crossings; the box test reads only q0's own pair.

Lookup first: a crossing's backward representative is the iterate table's answer
when it has one, and otherwise Trellis.image_cdist -- the one crossing-level scaling
path (table first, then per_step_beta).

Branch cycle / inversion: FixedPoint.branch_position_map follows advance_key, so the
inversion ordering (k_value == 2*period, branch_index flips after each full orbit)
needs no special case; it is exercised by the inversion law cases (the k=4 saddle
at 0/1/2 blasts and the k=4.1 period-2 orbit).
"""

from __future__ import annotations

import logging
import math
from typing import Iterable, Optional, TYPE_CHECKING

from .TopologyResults import StrongPipResult
from .Trellis import same_point

if TYPE_CHECKING:
    from ..numerics.FixedPoint import FixedPoint
    from .Trellis import Trellis

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def same_branch_cycle(a: "FixedPoint", b: "FixedPoint") -> bool:
    """True when ``a`` and ``b`` have one period and one ``k_value``, so their branch
    cycles advance together (a heteroclinic pip or ring between them is defined)."""
    return (a.period, a.k_value) == (b.period, b.k_value)


def is_strong_pip(
    trellis: "Trellis",
    intersection_id: int,
    *,
    collision_rtol: float = 1e-2,
) -> StrongPipResult:
    """
    Classify whether a single intersection is a strong pip.

    Implements Is-Strong-Pip(q0, T): every crossing of q0's own pair (the unstable
    manifold of q0's unstable fixed point and the stable branch cycle of its
    stable one, module Dev Notes) is mapped back onto q0's stable branch using
    the canonical-distance scaling law (one map step is
    :meth:`FixedPoint.per_step_beta`), and q0 is a strong pip iff none of them lands
    strictly inside the open box (0, s0) x (0, u0) — on the stable arc
    W^S(z'_0, q0) and inside the unstable arc.

    q0's own orbit maps back onto q0 itself (a collision on both canonical
    distances), so without guarding against it q0's iterates falsely disqualify q0 —
    exactly the failure seen for period > 1, where M(q0) maps back onto q0. The guard
    is a direct cdist collision check (match on BOTH cdists), which only ever skips
    q0's genuine iterates; a distinct iterate chain that merely happens to be nearby
    (or shares the same invariant action product) is not mistaken for q0.

    Args:
        trellis: The Trellis to read intersections and eigen/orbit data from.
        intersection_id: Registry ID of the candidate q0.
        collision_rtol: Relative slack for the cdist collision test that identifies
            q0's own orbit among crossings the iterate table does NOT link (a
            linked crossing is mapped back by lookup and recognised by id). An
            unlinked point is treated as q0 (and skipped) only if its scaled
            stable AND unstable cdists are both within ``collision_rtol`` of q0's.
            Defaults to 1e-2, which absorbs canonical-distance scaling noise while
            still distinguishing distinct intersections.

    Returns:
        A StrongPipResult. When not a strong pip, ``blocking_intersection_id`` is
        the intersection whose mapped-back unstable cdist is smallest (the innermost
        disqualifier).

    Raises:
        ValueError: If q0 has no stable branch, q0's branch is missing from
            the computed branch cycle, or q0 is heteroclinic between fixed
            points of different periods or ``k_value`` (:func:`same_branch_cycle`).

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

    fixed_point, unstable_fixed_point = branch_key[0], q0.manifold_a_key[0]
    if not same_branch_cycle(unstable_fixed_point, fixed_point):
        raise ValueError(
            f"Intersection {intersection_id} joins fixed points of periods "
            f"{unstable_fixed_point.period} and {fixed_point.period} (k_value "
            f"{unstable_fixed_point.k_value} and {fixed_point.k_value}); a heteroclinic "
            "strong pip needs both of one period and one k_value"
        )
    s0, u0 = q0.cdist_key("stable"), q0.cdist_key("unstable")
    pos_map, k = fixed_point.branch_position_map("stable")

    pos_B = pos_map.get(branch_key)
    if pos_B is None:
        raise ValueError(
            f"Stable branch {branch_key} of intersection {intersection_id} is not "
            "in the forward branch cycle; the trellis may be inconsistent."
        )

    # Map every intersection back onto q0's stable branch by the minimal i_r
    # backward steps in the branch cycle, scaling stable cdist by beta^(i_r) and
    # unstable cdist by beta^(-i_r) (both map-free). A disqualifier is one
    # that lands strictly inside the open box (0, s0) x (0, u0): on the stable arc
    # W^S(z'_0, q0) (stable_rep < s0) and inside the unstable arc (unstable_rep < u0).
    # Among those, the innermost (smallest mapped unstable cdist) is reported.
    best_rep = (math.inf, 0.0)
    best_id: Optional[int] = None
    for r_id, r in trellis.registry:
        if r_id == intersection_id or r_id not in trellis.resolved_ids:
            continue  # an unresolved detection is no crossing and disqualifies nothing
        if r.is_anchor:
            continue  # the anchoring fixed point, not a transverse intersection
        # Only q0's own pair disqualifies it (module Dev Notes): r's unstable
        # branch on q0's unstable fixed point, and (pos_r below) r's stable branch
        # in q0's stable cycle.
        if r.manifold_a_key[0] is not unstable_fixed_point:
            continue
        pos_r = pos_map.get(r.manifold_b_key)
        if pos_r is None:
            # Stable branch of a different fixed point — cannot map onto B.
            continue
        backward_steps = (pos_r - pos_B) % k
        # Lookup first: when the iterate table knows M^-backward_steps(r), the
        # representative is a registered crossing with exact identity and
        # cdists. q0's own orbit then maps back onto q0 itself, by id.
        back = trellis.iterate(r_id, -backward_steps)
        if back is not None:
            if back == intersection_id:
                continue
            stable_rep, unstable_rep = (trellis.intersection(back).cdist_key(side) for side in ("stable", "unstable"))
        else:
            # Unlinked: backward_steps applications of M^-1 by the cdist
            # definition, through the one scaling path.
            _, stable_rep, _ = trellis.image_cdist_key(r_id, -backward_steps, "stable")
            _, unstable_rep, _ = trellis.image_cdist_key(r_id, -backward_steps, "unstable")

            # q0's own orbit maps back onto q0 itself — a collision on BOTH
            # cdists. That is the same point, not a disqualifier, so skip it.
            # Requiring a collision on both cdists (a 2D coincidence) means only
            # q0's genuine iterates are skipped; a distinct chain that merely
            # shares q0's stable OR unstable cdist — or an equal "action"
            # product — is not mistaken for q0.
            if same_point(unstable_rep[0], stable_rep[0], u0[0], s0[0], collision_rtol):
                continue

        # Cdist keys (Intersection.cdist_key; a scaled one is (c, 0.0)).
        if stable_rep < s0 and unstable_rep < u0 and unstable_rep < best_rep:
            best_rep = unstable_rep
            best_id = r_id

    if best_id is not None:
        logger.debug(
            "Intersection %s is not a strong pip: %s maps into the box at u=%.6g < %.6g",
            intersection_id,
            best_id,
            best_rep[0],
            u0[0],
        )
        return StrongPipResult(intersection_id, False, blocking_intersection_id=best_id)

    return StrongPipResult(intersection_id, True)


def classify_strong_pips(
    trellis: "Trellis",
    intersection_ids: Optional[Iterable[int]] = None,
    *,
    collision_rtol: float = 1e-2,
) -> dict[int, StrongPipResult]:
    """
    Classify many intersections, sharing the branch-cycle cache (pure query).

    Returns a result per intersection without mutating the trellis, skipping the
    anchors and the heteroclinic crossings whose fixed points differ in period
    or ``k_value`` (neither can be a strong pip, :func:`same_branch_cycle`). The qualifying
    results (``is_strong_pip is True``) are the strong-pip *candidates*; selecting
    the single actual strong pip from them is the caller's job — see
    Trellis.classify_strong_pips() / set_strong_pip().

    Args:
        trellis: The Trellis to classify within.
        intersection_ids: Candidates to classify. Defaults to every intersection
            in the trellis.
        collision_rtol: Relative collision slack passed through to is_strong_pip.

    Returns:
        Mapping from intersection ID to its StrongPipResult.
    """
    if intersection_ids is None:
        intersection_ids = trellis.registry.all_ids()

    results: dict[int, StrongPipResult] = {}
    for iid in intersection_ids:
        ix = trellis.intersection(iid)
        # The anchor z'_0 is not a transverse crossing: its open interval
        # W^S(z'_0, q0) is empty, so classifying it as a pip is meaningless.
        if ix.is_anchor:
            logger.debug("Skipping anchor/fixed-point intersection %s", iid)
            continue
        if not same_branch_cycle(ix.manifold_a_key[0], ix.manifold_b_key[0]):
            continue  # heteroclinic across branch cycles: no strong pip (is_strong_pip raises)
        results[iid] = is_strong_pip(trellis, iid, collision_rtol=collision_rtol)

    logger.debug(
        "Classified %d intersection(s): %d candidate(s)",
        len(results),
        sum(1 for r in results.values() if r.is_strong_pip),
    )
    return results
