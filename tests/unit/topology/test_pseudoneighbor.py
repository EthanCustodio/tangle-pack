"""Compute-Pseudoneighbors: reference pairs, the all-iterates interval check,
and trajectory extension.

The synthetic tests hand-build a registry (explicit canonical distances and
branch keys) plus the TrellisBranch ordering, so every geometric case is pinned
exactly. The structural validity of the reference pairs on real tangles is the
law ``reference_pairs_valid`` of ``tests/invariants/test_law_partition.py``.

The definition tests BOTH open intervals of a pair (B2), and every point of X
visits the window once a lap, so an extra crossing on the window's stable branch
lands strictly inside some pair's stable interval unless it is an iterate of a
window crossing. A crossing that probes the unstable interval alone therefore
sits on the saddle's OTHER stable branch, ``OTHER_STABLE`` (branch index 1, an
independent chain without inversion): its stable orbit never visits the window.
"""

from __future__ import annotations

import inspect
import logging

import pytest

from helpers.logs import assert_logged
from helpers.fakes import bare_fixed_point, synthetic_trellis as _trellis
from tanglepack.numerics.FixedPoint import FixedPoint
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry
from tanglepack.topology.Pseudoneighbor import (
    compute_pseudoneighbors,
    extend_pseudoneighbor_trajectories,
)
from tanglepack.topology.TopologyResults import PseudoneighborPair
from tanglepack.topology.Trellis import Trellis

#: The no-inversion saddle's second stable branch: never on the window's cycle.
OTHER_STABLE = 1

#: The default slack of the 2D endpoint-collision fallback; the table-linked
#: tests drift a landing by TWICE this, so the unlinked fallback must reject it.
COLLISION_RTOL = inspect.signature(compute_pseudoneighbors).parameters["collision_rtol"].default


def _fixed_point(period: int, lambda_u: float) -> FixedPoint:
    """A bare no-inversion fixed point whose full-cycle unstable eigenvalue is ``lambda_u``."""
    return bare_fixed_point(period, beta=lambda_u ** (-1.0 / period))


def _reference_window(fp: FixedPoint, registry: IntersectionRegistry):
    """The standard period-1 (lambda=4) reference window: r_n at stable cdist 4,
    one intermediate point, and r_end = M(r_n) at stable cdist 1, with the
    iterate link registered. Returns (r_n, x1, r_end) ids."""
    stable = (fp, "stable", 0, 0)
    unstable = (fp, "unstable", 0, 0)
    r_n = registry.add_synthetic(
        (4.0, 1.0), unstable_cdist=1.0, stable_cdist=4.0,
        manifold_a_key=unstable, manifold_b_key=stable,
    )
    x1 = registry.add_synthetic(
        (3.0, 2.0), unstable_cdist=2.0, stable_cdist=3.0,
        manifold_a_key=unstable, manifold_b_key=stable,
    )
    r_end = registry.add_synthetic(
        (1.0, 4.0), unstable_cdist=4.0, stable_cdist=1.0,
        manifold_a_key=unstable, manifold_b_key=stable,
    )
    registry.register_iterate(r_n, 1, r_end)
    return r_n, x1, r_end


def test_reference_pairs_found_on_clean_interval():
    """Consecutive window pairs with empty unstable intervals are references."""
    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
    r_n, x1, r_end = _reference_window(fp, reg)

    pairs = compute_pseudoneighbors(_trellis(reg, fp))

    assert [(p.intersection_a, p.intersection_b) for p in pairs] == [
        (r_n, x1),
        (x1, r_end),
    ]
    assert all(p.is_reference for p in pairs)
    assert all(p.branch_key == (fp, "stable", 0, 0) for p in pairs)


def test_direct_intersection_punctures_interval():
    """A crossing whose unstable cdist lies inside (u0, u1) rejects the pair."""
    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
    r_n, x1, r_end = _reference_window(fp, reg)
    # Inside the (1.0, 2.0) unstable interval of the (r_n, x1) pair at n=0;
    # outside the (2.0, 4.0) interval at every iterate (1.5 * 4^m).
    reg.add_synthetic(
        (0.5, 1.5), unstable_cdist=1.5, stable_cdist=0.5,
        manifold_a_key=(fp, "unstable", 0, 0), manifold_b_key=(fp, "stable", 0, OTHER_STABLE),
    )

    pairs = compute_pseudoneighbors(_trellis(reg, fp))

    assert [(p.intersection_a, p.intersection_b) for p in pairs] == [(x1, r_end)]


def test_iterate_punctures_interval():
    """A crossing outside the interval whose iterate lands inside rejects it.

    The candidate sits at unstable cdist 0.4; one forward map step (lambda = 4)
    carries it to 1.6, inside the (1.0, 2.0) interval — pins the log-ratio
    m-range arithmetic.
    """
    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
    r_n, x1, r_end = _reference_window(fp, reg)
    reg.add_synthetic(
        (0.5, 0.4), unstable_cdist=0.4, stable_cdist=0.5,
        manifold_a_key=(fp, "unstable", 0, 0), manifold_b_key=(fp, "stable", 0, OTHER_STABLE),
    )

    pairs = compute_pseudoneighbors(_trellis(reg, fp))

    assert [(p.intersection_a, p.intersection_b) for p in pairs] == [(x1, r_end)]


def test_endpoint_iterates_do_not_disqualify():
    """A crossing that is a (noisy) genuine iterate of a pair member — BOTH of
    its scaled cdists collide with the endpoint — does not reject the pair,
    while a half-matching crossing near the boundary does."""
    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
    r_n, x1, r_end = _reference_window(fp, reg)
    stable = (fp, "stable", 0, 0)
    unstable = (fp, "unstable", 0, 0)
    # M(x1) with cdist noise: backward landing (2.002, 3.002) collides with
    # x1 (u=2, s=3) on BOTH cdists — its own orbit, not a puncture.
    reg.add_synthetic(
        (0.75, 8.0), unstable_cdist=8.008, stable_cdist=0.7505,
        manifold_a_key=unstable, manifold_b_key=stable,
    )

    pairs = compute_pseudoneighbors(_trellis(reg, fp))
    assert len(pairs) == 2

    # Control: same unstable landing (2.002 in the open (2, 4) interval) but a
    # stable position far from x1's — a genuine distinct point, punctures.
    reg.add_synthetic(
        (0.2, 8.008), unstable_cdist=8.008, stable_cdist=0.2,
        manifold_a_key=unstable, manifold_b_key=stable,
    )
    pairs = compute_pseudoneighbors(_trellis(reg, fp))
    assert (x1, r_end) not in [(p.intersection_a, p.intersection_b) for p in pairs]


@pytest.mark.parametrize("middle_branch", [0, 1], ids=["same_branch", "other_branch"])
def test_pair_on_different_unstable_branches_rejected(middle_branch: int):
    """Consecutive stable-branch points on two different unstable branches are
    never a pair -- no single unstable arc connects them.

    The standard window, with the middle crossing's unstable side on the
    saddle's OTHER unstable branch (branch index 1, an independent chain of a
    point without inversion); ``same_branch`` is the control.
    """
    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
    stable = (fp, "stable", 0, 0)
    unstable = (fp, "unstable", 0, 0)
    r_n = reg.add_synthetic(
        (4.0, 1.0), unstable_cdist=1.0, stable_cdist=4.0,
        manifold_a_key=unstable, manifold_b_key=stable,
    )
    x1 = reg.add_synthetic(
        (3.0, 2.0), unstable_cdist=2.0, stable_cdist=3.0,
        manifold_a_key=(fp, "unstable", 0, middle_branch), manifold_b_key=stable,
    )
    r_end = reg.add_synthetic(
        (1.0, 4.0), unstable_cdist=4.0, stable_cdist=1.0,
        manifold_a_key=unstable, manifold_b_key=stable,
    )
    reg.register_iterate(r_n, 1, r_end)

    pairs = compute_pseudoneighbors(_trellis(reg, fp))
    found = [(p.intersection_a, p.intersection_b) for p in pairs]

    if middle_branch == 0:
        assert found == [(r_n, x1), (x1, r_end)]
    else:
        assert (r_n, x1) not in found and (x1, r_end) not in found


def test_branch_residue_gates_iterates():
    """Period-3: a candidate lands on the target unstable branch only for step
    counts matching its cycle residue; a cdist that would fall inside the
    interval on the WRONG branch does not disqualify."""
    fp = _fixed_point(3, 8.0)  # k = 3, beta = 2
    reg = IntersectionRegistry()
    stable = (fp, "stable", 0, 0)
    unstable0 = (fp, "unstable", 0, 0)
    r_n = reg.add_synthetic(
        (4.0, 1.0), unstable_cdist=1.0, stable_cdist=4.0,
        manifold_a_key=unstable0, manifold_b_key=stable,
    )
    x1 = reg.add_synthetic(
        (3.0, 2.0), unstable_cdist=2.0, stable_cdist=3.0,
        manifold_a_key=unstable0, manifold_b_key=stable,
    )
    r_end = reg.add_synthetic(
        (0.5, 8.0), unstable_cdist=8.0, stable_cdist=0.5,
        manifold_a_key=unstable0, manifold_b_key=stable,
    )
    # The table holds M^3(r_n) = r_end as one direct link, no single steps.
    reg.register_iterate(r_n, 3, r_end)

    # On unstable branch 1 (cycle position 1): reaches branch 0 only after
    # d = (0 - 1) % 3 = 2 steps, giving cdist 0.6 * 2^2 = 2.4 — outside
    # (1.0, 2.0) at every full return. Naively applying d = 1 would give 1.2,
    # inside the interval.
    reg.add_synthetic(
        (0.6, 0.6), unstable_cdist=0.6, stable_cdist=0.6,
        manifold_a_key=(fp, "unstable", 1, 0), manifold_b_key=stable,
    )

    pairs = compute_pseudoneighbors(_trellis(reg, fp))
    assert (r_n, x1) in [(p.intersection_a, p.intersection_b) for p in pairs]

    # Control: the same cdist on branch 0 itself (residue 0) disqualifies —
    # 0.15 * 8 = 1.2 lands inside (1.0, 2.0) after one full return.
    reg.add_synthetic(
        (0.7, 0.15), unstable_cdist=0.15, stable_cdist=0.7,
        manifold_a_key=unstable0, manifold_b_key=stable,
    )
    pairs = compute_pseudoneighbors(_trellis(reg, fp))
    assert (r_n, x1) not in [(p.intersection_a, p.intersection_b) for p in pairs]


def test_inversion_branch_return_scales_by_beta_to_the_k():
    """Inversion: a branch return is ``k_value = 2`` map steps and scales an
    unstable cdist by ``beta ** 2`` (= ``lambda_u ** 2``), not by ``lambda_u``.

    A candidate at u = 1.5 on the pair's own unstable branch returns to it at
    1.5 * 4 = 6, outside the pair's interval (2, 4); a ``lambda_u`` lap would put
    it at 1.5 * 2 = 3, inside, and puncture the pair."""
    fp = bare_fixed_point(1, inversion=True, beta=0.5)  # lambda_u = -2, beta = 2
    assert fp.k_value == 2
    reg = IntersectionRegistry()
    s0, u0 = (fp, "stable", 0, 0), (fp, "unstable", 0, 0)
    s1, u1 = (fp, "stable", 0, 1), (fp, "unstable", 0, 1)

    def _point(s, u, s_key, u_key):
        return reg.add_synthetic(
            (s, u), unstable_cdist=u, stable_cdist=s,
            manifold_a_key=u_key, manifold_b_key=s_key,
        )

    r_n = _point(4.0, 1.0, s0, u0)
    x1 = _point(3.0, 2.0, s0, u0)
    middle = _point(2.0, 2.0, s1, u1)  # M(r_n), on the other half-branch
    r_end = _point(1.0, 4.0, s0, u0)  # M^2(r_n): one full return
    reg.register_iterate(r_n, 1, middle)
    reg.register_iterate(middle, 1, r_end)
    _point(0.9, 1.5, s0, u0)

    pairs = compute_pseudoneighbors(_trellis(reg, fp))

    assert (x1, r_end) in [(p.intersection_a, p.intersection_b) for p in pairs]


def test_walk_stops_at_reference_window_end():
    """Pairs below r_{n+p} are outside the reference window and not returned."""
    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
    r_n, x1, r_end = _reference_window(fp, reg)
    # At M(r_end)'s cdists (unlinked): its orbit visits the window at r_end itself.
    below = reg.add_synthetic(
        (0.25, 16.0), unstable_cdist=16.0, stable_cdist=0.25,
        manifold_a_key=(fp, "unstable", 0, 0), manifold_b_key=(fp, "stable", 0, 0),
    )

    pairs = compute_pseudoneighbors(_trellis(reg, fp))

    members = {i for p in pairs for i in p.as_tuple()}
    assert below not in members
    assert len(pairs) == 2


def test_reference_end_is_a_table_iterate(caplog):
    """r_{n+p} = M^p(r_n) comes from the iterate table only: without the link
    there is no reference window (one WARNING), however well a crossing's
    cdists would match the scaled ones."""
    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
    stable = (fp, "stable", 0, 0)
    unstable = (fp, "unstable", 0, 0)
    for s_cdist, u_cdist in ((4.0, 1.0), (3.0, 2.0), (1.0, 4.0)):
        reg.add_synthetic(
            (s_cdist, u_cdist), unstable_cdist=u_cdist, stable_cdist=s_cdist,
            manifold_a_key=unstable, manifold_b_key=stable,
        )

    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.Pseudoneighbor"):
        pairs = compute_pseudoneighbors(_trellis(reg, fp))

    assert pairs == []
    assert_logged(caplog, logging.WARNING, "tanglepack.topology.Pseudoneighbor", count=1)


def test_heteroclinic_candidate_ignored():
    """A candidate whose unstable side belongs to a different fixed point can
    never iterate onto this tangle's unstable manifold, so it cannot puncture."""
    fp = _fixed_point(1, 4.0)
    other = _fixed_point(1, 6.0)
    reg = IntersectionRegistry()
    _reference_window(fp, reg)
    reg.add_synthetic(
        (0.5, 1.5), unstable_cdist=1.5, stable_cdist=0.5,
        manifold_a_key=(other, "unstable", 0, 0), manifold_b_key=(fp, "stable", 0, OTHER_STABLE),
    )

    pairs = compute_pseudoneighbors(_trellis(reg, fp, other))

    assert len(pairs) == 2


def test_reference_window_starts_at_chosen_strong_pip():
    """With a strong pip chosen, r_n is its cut point — not the outermost
    intersection — and the trellis cuts N there (D10): crossings beyond the cut
    are no crossings of the trellis, so they leave the check set."""
    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
    r_n, x1, r_end = _reference_window(fp, reg)
    stable = (fp, "stable", 0, 0)
    unstable = (fp, "unstable", 0, 0)
    # Beyond the pip: the manifold's (untrimmed) outermost intersection ...
    outer = reg.add_synthetic(
        (16.0, 0.25), unstable_cdist=0.25, stable_cdist=16.0,
        manifold_a_key=unstable, manifold_b_key=stable,
    )
    # ... and a crossing whose iterate (0.375 * 4 = 1.5) would puncture the
    # (1.0, 2.0) interval if it were not excluded as beyond the cut.
    reg.add_synthetic(
        (8.0, 0.375), unstable_cdist=0.375, stable_cdist=8.0,
        manifold_a_key=unstable, manifold_b_key=stable,
    )

    trellis = _trellis(reg, fp)
    trellis.strong_pip = r_n  # cut the manifold at the pip (stable cdist 4)

    pairs = compute_pseudoneighbors(trellis)

    assert [(p.intersection_a, p.intersection_b) for p in pairs] == [
        (r_n, x1),
        (x1, r_end),
    ]
    assert outer not in {i for p in pairs for i in p.as_tuple()}


def test_single_reference_window_per_fixed_point():
    """Period-3: the reference window lives ONLY on the branch of greatest
    stable cdist (the pip's branch); the trajectory's appearances on the other
    branches come from extension through the iterate table, never as
    independent references."""
    fp = _fixed_point(3, 8.0)  # k = 3, beta = 2
    reg = IntersectionRegistry()

    def _point(orbit, s, u):
        return reg.add_synthetic(
            (s, u), unstable_cdist=u, stable_cdist=s,
            manifold_a_key=(fp, "unstable", orbit, 0),
            manifold_b_key=(fp, "stable", orbit, 0),
        )

    # Branch 0 carries the globally outermost intersection: the window
    # [4 -> 0.5] with a middle point, closed by the direct link M^3(a0) = e0.
    a0, b0, e0 = _point(0, 4.0, 1.0), _point(0, 3.0, 2.0), _point(0, 0.5, 8.0)
    reg.register_iterate(a0, 3, e0)
    # Branch 1 holds the forward images (stable /beta, unstable *beta).
    a1, b1, e1 = _point(1, 2.0, 2.0), _point(1, 1.5, 4.0), _point(1, 0.25, 16.0)
    for source, image in ((a0, a1), (b0, b1), (e0, e1)):
        reg.register_iterate(source, 1, image)

    trellis = _trellis(reg, fp)
    references = compute_pseudoneighbors(trellis)

    assert [(p.intersection_a, p.intersection_b) for p in references] == [
        (a0, b0),
        (b0, e0),
    ]
    assert all(p.branch_key == (fp, "stable", 0, 0) for p in references)

    extended = extend_pseudoneighbor_trajectories(trellis, references)
    assert {(p.as_tuple(), p.iterate) for p in extended} == {
        ((min(a1, b1), max(a1, b1)), 1),
        ((min(b1, e1), max(b1, e1)), 1),
    }


def test_backward_iterate_punctures_interval():
    """The definition's X runs over ALL n in Z: a candidate's BACKWARD iterate
    landing inside the pair's open unstable interval rejects the pair, even
    though that iterate's stable position lies beyond the pip cut (on the
    removed tail) — e.g. the preimages of blasted crossings. The candidate itself
    is a crossing of T(q0): the pip's orbit is registered to M^2(q0), so the
    unstable branch is cut past it (D12)."""
    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
    r_n, x1, r_end = _reference_window(fp, reg)
    stable = (fp, "stable", 0, 0)
    unstable = (fp, "unstable", 0, 0)
    deepest = reg.add_synthetic(
        (0.25, 16.0), unstable_cdist=16.0, stable_cdist=0.25,
        manifold_a_key=unstable, manifold_b_key=stable,
    )
    reg.register_iterate(r_end, 1, deepest)
    # Backward iterate lands at u = 6 / 4 = 1.5 inside (1, 2); its stable
    # position 2 * 4 = 8 exceeds the cut, and it still punctures.
    reg.add_synthetic(
        (2.0, 6.0), unstable_cdist=6.0, stable_cdist=2.0,
        manifold_a_key=unstable, manifold_b_key=stable,
    )

    trellis = _trellis(reg, fp)
    trellis.strong_pip = r_n

    pairs = compute_pseudoneighbors(trellis)

    assert (r_n, x1) not in [(p.intersection_a, p.intersection_b) for p in pairs]


def _three_step_chain(fp, reg):
    """Two crossings and their first two forward images, the table linked."""
    stable = (fp, "stable", 0, 0)
    unstable = (fp, "unstable", 0, 0)

    def _point(u, s):
        return reg.add_synthetic(
            (s, u), unstable_cdist=u, stable_cdist=s,
            manifold_a_key=unstable, manifold_b_key=stable,
        )

    a0, b0 = _point(1.0, 4.0), _point(2.0, 3.0)
    a1, b1 = _point(4.0, 1.0), _point(8.0, 0.75)
    a2, b2 = _point(16.0, 0.25), _point(32.0, 0.1875)
    for src, dst in ((a0, a1), (a1, a2), (b0, b1), (b1, b2)):
        reg.register_iterate(src, 1, dst)
    return (a0, b0), (a1, b1), (a2, b2)


def test_trajectory_extension_forward_stops_at_the_end_of_the_chain():
    """A reference maps forward through the iterate table into non-reference
    pairs, stopping at the end of a chain."""
    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
    first, second, third = _three_step_chain(fp, reg)
    stable = (fp, "stable", 0, 0)

    trellis = _trellis(reg, fp)
    ref = PseudoneighborPair(*first, branch_key=stable, is_reference=True)
    out = extend_pseudoneighbor_trajectories(trellis, [ref])

    assert [(p.intersection_a, p.intersection_b) for p in out] == [second, third]
    assert all(not p.is_reference for p in out)
    assert all(p.branch_key == stable for p in out)
    assert [p.iterate for p in out] == [1, 2]


@pytest.mark.parametrize(
    "references, expected",
    [((0, 1), [2]), ((0, 2), [1])],
    ids=["against_a_reference", "against_another_trajectory"],
)
def test_trajectory_extension_deduplicates(references, expected):
    """Two references on ONE trajectory: the forward chain of one runs into the
    other reference (never re-emitted), and the gap between two references is
    emitted once although both walk into it."""
    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
    chain = _three_step_chain(fp, reg)
    stable = (fp, "stable", 0, 0)

    trellis = _trellis(reg, fp)
    refs = [
        PseudoneighborPair(*chain[i], branch_key=stable, is_reference=True)
        for i in references
    ]
    out = extend_pseudoneighbor_trajectories(trellis, refs)

    emitted = [(p.intersection_a, p.intersection_b) for p in out]
    assert emitted == [chain[i] for i in expected]
    assert len(set(emitted)) == len(emitted)
    assert not set(emitted) & {chain[i] for i in references}


def test_trellis_wrapper_populates_and_clears_slots():
    """Trellis.compute_pseudoneighbors stores references (+ trajectories) in
    the pseudoneighbors slot; clear_results empties it."""
    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
    _reference_window(fp, reg)
    trellis = _trellis(reg, fp)

    references = trellis.compute_pseudoneighbors()

    assert trellis.reference_pseudoneighbors == references
    assert len(references) == 2
    trellis.clear_results()
    assert trellis.pseudoneighbors == []


# --------------------------------------------------------------------------- #
# Strong-pip cut coverage (plan row 1.18)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "period, placed_warnings", [(1, 0), (3, 1)], ids=["period_1", "period_3_placed"]
)
def test_missing_pip_iterates_are_cut_at_their_predicted_cdist(caplog, period, placed_warnings):
    """A period-3 anchor is cut on every stable branch (PK-6, D5): with no
    ``+1``/``+2`` links for the pip, its cuts on the other two branches are
    placed by ``image_cdist`` (stable cdist ``4 * beta^j``) and flagged in one
    WARNING of the trellis, which cuts N when the pip is chosen (D10). A
    crossing beyond such a predicted cut is no crossing of the trellis, so it
    cannot puncture a window pair. A period-1 anchor has one stable branch,
    which the pip alone cuts: nothing to place, nothing to report.

    The window ``(4, 1) -> (3, 2) -> M^k = (4 / lambda, lambda)`` sits on orbit
    point 0 with the pip at its outer end and the table holds the one ``M^k``
    link. ``z`` (period 3 only) lies on ``s(1, 0)`` at stable cdist 3, past the
    predicted cut 4 * 0.5 = 2, at unstable cdist 1.5 inside the first pair's
    ``(1, 2)``: uncut, it would puncture that pair.
    """
    lambda_u = 4.0 if period == 1 else 8.0
    fp = _fixed_point(period, lambda_u)
    reg = IntersectionRegistry()

    def _point(stable_cdist, unstable_cdist, stable_orbit=0):
        return reg.add_synthetic(
            (stable_cdist, unstable_cdist),
            unstable_cdist=unstable_cdist, stable_cdist=stable_cdist,
            manifold_a_key=(fp, "unstable", 0, 0),
            manifold_b_key=(fp, "stable", stable_orbit, 0),
        )

    pip = _point(4.0, 1.0)
    x1 = _point(3.0, 2.0)
    r_end = _point(4.0 / lambda_u, lambda_u)
    reg.register_iterate(pip, period, r_end)
    if period == 3:
        _point(3.0, 1.5, stable_orbit=1)
    trellis = _trellis(reg, fp)
    assert fp.k_value == period

    with caplog.at_level(logging.WARNING):
        trellis.strong_pip = pip
        pairs = compute_pseudoneighbors(trellis)

    assert [p.as_tuple() for p in pairs] == [
        tuple(sorted((pip, x1))), tuple(sorted((x1, r_end)))
    ]
    assert_logged(caplog, logging.WARNING, "tanglepack.topology.Trellis", count=placed_warnings)
    assert_logged(caplog, logging.WARNING, "tanglepack.topology.Pseudoneighbor", count=0)


@pytest.mark.parametrize(
    "period, linked_steps, emitted",
    [(3, 1, [1]), (3, 0, []), (1, 0, [])],
    ids=["period_3_plus_2_missing", "period_3_unlinked", "period_1_needs_no_image"],
)
def test_missing_forward_appearance_is_silent(caplog, period, linked_steps, emitted):
    """A reference trajectory is extended forward only while the table links
    its members, and a missing ``+j`` appearance on another branch of a
    period-k orbit costs nothing: forward pairs are never holed (author,
    2026-10-09, reversing D5), so nothing is invented and nothing warns."""
    fp = _fixed_point(period, 8.0)
    reg = IntersectionRegistry()

    def _point(orbit, s, u):
        return reg.add_synthetic(
            (s, u), unstable_cdist=u, stable_cdist=s,
            manifold_a_key=(fp, "unstable", orbit, 0),
            manifold_b_key=(fp, "stable", orbit, 0),
        )

    a0, b0 = _point(0, 4.0, 1.0), _point(0, 3.0, 2.0)
    if linked_steps:
        reg.register_iterate(a0, 1, _point(1, 2.0, 2.0))
        reg.register_iterate(b0, 1, _point(1, 1.5, 4.0))
    trellis = _trellis(reg, fp)
    ref = PseudoneighborPair(a0, b0, branch_key=(fp, "stable", 0, 0), is_reference=True)

    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.Pseudoneighbor"):
        out = extend_pseudoneighbor_trajectories(trellis, [ref])

    assert [p.iterate for p in out] == emitted
    assert_logged(caplog, logging.WARNING, "tanglepack.topology.Pseudoneighbor", count=0)


def test_table_linked_deep_iterate_does_not_disqualify():
    """A pair member's own deep iterate is recognised by iterate-table LOOKUP,
    not by the cdist collision. Six forward links from x1 end at a crossing whose
    registered stable cdist has drifted 2 * collision_rtol — the k=2.8
    script at 8 blasts lost (10, 9) exactly this way — yet the pair survives
    because the table maps the landing back onto x1 by id. Without the links
    the collision fallback still rejects it (the unlinked behaviour is
    unchanged)."""
    fp = _fixed_point(1, 4.0)
    stable = (fp, "stable", 0, 0)
    unstable = (fp, "unstable", 0, 0)

    def chain(reg, x1, link: bool):
        prev = x1
        for n in range(1, 7):
            u = 2.0 * 4.0**n
            s = 3.0 / 4.0**n
            if n == 6:
                u *= 1.001  # scaled landing 2.002: strictly inside (2, 4)
                # Scaled stable 2.94 vs x1's 3: no collision, and inside the stable
                # interval of (x1, r_end) only, which the landing punctures anyway.
                s *= 1.0 - 2.0 * COLLISION_RTOL
            nxt = reg.add_synthetic(
                (s, u), unstable_cdist=u, stable_cdist=s,
                manifold_a_key=unstable, manifold_b_key=stable,
            )
            if link:
                reg.register_iterate(prev, 1, nxt)
            prev = nxt

    reg = IntersectionRegistry()
    r_n, x1, r_end = _reference_window(fp, reg)
    chain(reg, x1, link=True)
    pairs = compute_pseudoneighbors(_trellis(reg, fp))
    assert [(p.intersection_a, p.intersection_b) for p in pairs] == [
        (r_n, x1),
        (x1, r_end),
    ]

    reg = IntersectionRegistry()
    r_n, x1, r_end = _reference_window(fp, reg)
    chain(reg, x1, link=False)
    pairs = compute_pseudoneighbors(_trellis(reg, fp))
    assert [(p.intersection_a, p.intersection_b) for p in pairs] == [(r_n, x1)]


def test_table_linked_landing_uses_registered_cdist():
    """When the table knows a candidate's landing, the landing's REGISTERED
    unstable cdist decides, not the scaled estimate. Here the estimate (3.0)
    sits inside (2, 4) but the linked landing is registered at 4.5, outside, so
    the pair (x1, r_end) survives; unlinked, the same candidate punctures."""
    fp = _fixed_point(1, 4.0)
    unstable = (fp, "unstable", 0, 0)
    other_stable = (fp, "stable", 0, OTHER_STABLE)

    def candidate(reg, link: bool):
        z = reg.add_synthetic(
            (0.5, 4.5), unstable_cdist=4.5, stable_cdist=0.5,
            manifold_a_key=unstable, manifold_b_key=other_stable,
        )
        z1 = reg.add_synthetic(
            (0.125, 18.0), unstable_cdist=18.0, stable_cdist=0.125,
            manifold_a_key=unstable, manifold_b_key=other_stable,
        )
        r = reg.add_synthetic(
            (0.03125, 48.0), unstable_cdist=48.0, stable_cdist=0.03125,
            manifold_a_key=unstable, manifold_b_key=other_stable,
        )
        if link:
            reg.register_iterate(z, 1, z1)
            reg.register_iterate(z1, 1, r)

    reg = IntersectionRegistry()
    r_n, x1, r_end = _reference_window(fp, reg)
    candidate(reg, link=True)
    pairs = compute_pseudoneighbors(_trellis(reg, fp))
    assert (x1, r_end) in [(p.intersection_a, p.intersection_b) for p in pairs]

    reg = IntersectionRegistry()
    r_n, x1, r_end = _reference_window(fp, reg)
    candidate(reg, link=False)
    pairs = compute_pseudoneighbors(_trellis(reg, fp))
    assert (x1, r_end) not in [(p.intersection_a, p.intersection_b) for p in pairs]


@pytest.mark.parametrize(
    ("r_unstable", "survives"), [(128.2, True), (130.0, False)], ids=["below_c", "between_c_f"]
)
def test_unlinked_landing_is_read_against_the_pairs_images(r_unstable, survives):
    """An UNLINKED candidate is compared, on its own unstable branch, with the
    pair's REGISTERED images, not scaled back onto the pair (BCG-5).

    x1 and r_end are linked three steps forward to c (registered at 128.5, a
    drift of 0.5 off 2 * 4**3) and f (256). A crossing r at 128.2 with a far
    stable cdist scales back to 2.003, inside the pair's (2, 4), yet lies BELOW
    c on its own polyline, outside (128.5, 256): the image of W^U(x1, r_end)
    does not hold it, and the pair survives. At 130.0 it lies between c and f
    and punctures. Nested at 10 outer blasts lost its reference this way: the
    scaling noise at depth 8 exceeded the gap.
    """
    fp = _fixed_point(1, 4.0)
    stable = (fp, "stable", 0, 0)
    unstable = (fp, "unstable", 0, 0)
    reg = IntersectionRegistry()
    r_n, x1, r_end = _reference_window(fp, reg)
    for start, (u, s) in ((x1, (2.0, 3.0)), (r_end, (4.0, 1.0))):
        previous = start
        for n in (1, 2, 3):
            image_u = 128.5 if (start, n) == (x1, 3) else u * 4.0**n
            image = reg.add_synthetic(
                (s / 4.0**n, image_u), unstable_cdist=image_u, stable_cdist=s / 4.0**n,
                manifold_a_key=unstable, manifold_b_key=stable,
            )
            reg.register_iterate(previous, 1, image)
            previous = image
    reg.add_synthetic(
        (0.01, r_unstable), unstable_cdist=r_unstable, stable_cdist=0.01,
        manifold_a_key=unstable, manifold_b_key=(fp, "stable", 0, OTHER_STABLE),
    )

    pairs = [(p.intersection_a, p.intersection_b) for p in compute_pseudoneighbors(_trellis(reg, fp))]

    assert (r_n, x1) in pairs
    assert ((x1, r_end) in pairs) is survives


@pytest.mark.parametrize("noise", [1e-6, -1e-6], ids=["outward", "inward"])
def test_endpoint_own_unlinked_orbit_is_not_a_puncture(noise):
    """The third case of the frame reading above: r_end's images are LINKED,
    so the frame holds, but x1's own images are registered UNLINKED a noise
    hair off ``2 * 4**n``. In r's frame x1's end is placed by scaling, and its
    own orbit sitting a hair inside that placed end is x1 again (a 2D
    collision), not a puncture: both pairs survive at either sign."""
    fp = _fixed_point(1, 4.0)
    stable = (fp, "stable", 0, 0)
    unstable = (fp, "unstable", 0, 0)
    reg = IntersectionRegistry()
    r_n, x1, r_end = _reference_window(fp, reg)
    previous = r_end
    for n in (1, 2, 3):
        image = reg.add_synthetic(
            (1.0 / 4**n, 4.0 * 4**n), unstable_cdist=4.0 * 4**n, stable_cdist=1.0 / 4**n,
            manifold_a_key=unstable, manifold_b_key=stable,
        )
        reg.register_iterate(previous, 1, image)
        previous = image
    for n in (1, 2, 3):
        u = 2.0 * 4**n * (1 + noise)
        reg.add_synthetic(
            (3.0 / 4**n, u), unstable_cdist=u, stable_cdist=3.0 / 4**n,
            manifold_a_key=unstable, manifold_b_key=stable,
        )

    pairs = [(p.intersection_a, p.intersection_b) for p in compute_pseudoneighbors(_trellis(reg, fp))]

    assert (r_n, x1) in pairs
    assert (x1, r_end) in pairs


@pytest.mark.parametrize("frame", ["scaled", "own"])
def test_collision_on_another_stable_manifold_still_punctures(frame):
    """The 2D endpoint collision compares stable cdists, so it only names an
    endpoint's own orbit on the pair's OWN stable manifold (D7, WP9).

    A heteroclinic window of ``W^u(A)`` crossings on B's stable branch
    (A: lambda 4, B: lambda 6). An A x A crossing r lands one step back at
    unstable cdist 2.01, inside the pair (x1, r_end)'s (2, 4), and its stable
    cdist matches x1's on W^s(B) by numbers alone: read by scaling (no image
    of x1 or r_end registered) or in r's own frame (r_end's image registered,
    x1's placed), it still punctures, since r's orbit stays on W^s(A).
    """
    a, b = _fixed_point(1, 4.0), _fixed_point(1, 6.0)
    a_u, a_s, b_s = (a, "unstable", 0, 0), (a, "stable", 0, 0), (b, "stable", 0, 0)
    reg = IntersectionRegistry()

    def point(u, s, stable=b_s):
        return reg.add_synthetic((s, u), unstable_cdist=u, stable_cdist=s,
                                 manifold_a_key=a_u, manifold_b_key=stable)

    r_n, x1, r_end = point(1.0, 6.0), point(2.0, 3.0), point(4.0, 1.0)
    reg.register_iterate(r_n, 1, r_end)
    if frame == "own":
        reg.register_iterate(r_end, 1, point(16.0, 1.0 / 6.0))
    # x1's stable cdist one step forward on W^s(B) is 0.5, one step back from r 4 * s.
    point(8.04, 0.5 if frame == "own" else 0.75, stable=a_s)

    pairs = [(p.intersection_a, p.intersection_b)
             for p in compute_pseudoneighbors(_trellis(reg, a, b))]

    assert (r_n, x1) in pairs
    assert (x1, r_end) not in pairs


# --------------------------------------------------------------------------- #
# The open STABLE interval (B2): the same three readings as the unstable one
# --------------------------------------------------------------------------- #
def _pairs(trellis) -> list[tuple[int, int]]:
    return [(p.intersection_a, p.intersection_b) for p in compute_pseudoneighbors(trellis)]


def test_scaled_iterate_punctures_stable_interval():
    """A crossing below the window with no table links whose stable cdist, scaled
    one step back (0.5 * 4 = 2), lands inside the stable interval (1, 3) of
    (x1, r_end) kills that pair; neither end has a registered image, so only
    scaling can place it. Its unstable orbit (16 * 4^m) only touches the
    window's unstable ends, so the unstable interval stays clear: the stable
    interval decides."""
    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
    r_n, x1, r_end = _reference_window(fp, reg)
    reg.add_synthetic(
        (0.5, 16.0), unstable_cdist=16.0, stable_cdist=0.5,
        manifold_a_key=(fp, "unstable", 0, 0), manifold_b_key=(fp, "stable", 0, 0),
    )

    assert _pairs(_trellis(reg, fp)) == [(r_n, x1)]


def test_table_landing_punctures_stable_interval():
    """A table landing inside the stable interval kills the pair by id. The
    landing ``w`` (stable cdist 3.5 in (3, 4), unstable 8, outside both pairs'
    unstable intervals) is no crossing of the trellis, as past a cut, so x0 and
    x1 stay consecutive and ``w`` is no candidate; ``r = M(w)`` is one, and the
    table names ``M^-1(r) = w`` inside (r_n, x1)."""
    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
    r_n, x1, r_end = _reference_window(fp, reg)
    stable, unstable = (fp, "stable", 0, 0), (fp, "unstable", 0, 0)
    w = reg.add_synthetic(
        (3.5, 8.0), unstable_cdist=8.0, stable_cdist=3.5,
        manifold_a_key=unstable, manifold_b_key=stable,
    )
    r = reg.add_synthetic(
        (0.875, 32.0), unstable_cdist=32.0, stable_cdist=0.875,
        manifold_a_key=unstable, manifold_b_key=stable,
    )
    reg.register_iterate(w, 1, r)
    branches = _trellis(reg, fp).branches
    for branch in branches.values():
        branch.intersection_ids.remove(w)
    trellis = Trellis(fixed_points=[fp], registry=reg, branches=branches, bridges=[])
    assert w not in trellis and r in trellis

    assert _pairs(trellis) == [(x1, r_end)]


def test_own_orbit_does_not_puncture_stable_interval():
    """An unlinked crossing at M(x1)'s cdists up to noise (8.008, 0.7495) scales
    back to (2.002, 2.998): inside the stable interval (1, 3) of (x1, r_end), but
    it 2D-collides with x1, so it is x1's own orbit and the pair stays."""
    from tanglepack.topology.Pseudoneighbor import _interval_is_clear

    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
    r_n, x1, r_end = _reference_window(fp, reg)
    reg.add_synthetic(
        (0.7495, 8.008), unstable_cdist=8.008, stable_cdist=0.7495,
        manifold_a_key=(fp, "unstable", 0, 0), manifold_b_key=(fp, "stable", 0, 0),
    )
    trellis = _trellis(reg, fp)

    assert _interval_is_clear(
        trellis, x1, r_end, "stable", trellis.own_intersection_ids, COLLISION_RTOL
    )
    assert (x1, r_end) in _pairs(trellis)


@pytest.mark.parametrize(
    ("unstable_owner", "survives"), [("A", True), ("B", False)], ids=["own_orbit", "other_unstable"]
)
def test_collision_on_another_unstable_manifold_still_punctures_stable_interval(
    unstable_owner, survives
):
    """The stable reading's 2D collision compares UNSTABLE cdists, so it only names
    an endpoint's own orbit on the pair's own unstable manifold (D7 mirrored).

    A heteroclinic window of ``W^u(A)`` crossings on B's stable branch (A: lambda
    4, B: lambda 6), stable intervals (3, 6) and (1, 3). A crossing r on that
    branch at stable cdist 0.4985 lands one step back at 2.991, inside (1, 3),
    and its unstable cdist scaled one step back is 2.002, x1's 2 by numbers.
    On ``W^u(A)`` that is x1's own orbit up to noise and the pair stays; on
    ``W^u(B)`` (a B x B crossing) r's orbit never meets x1's, so it punctures."""
    a, b = _fixed_point(1, 4.0), _fixed_point(1, 6.0)
    a_u, b_u, b_s = (a, "unstable", 0, 0), (b, "unstable", 0, 0), (b, "stable", 0, 0)
    reg = IntersectionRegistry()

    def point(u, s, unstable=a_u):
        return reg.add_synthetic((s, u), unstable_cdist=u, stable_cdist=s,
                                 manifold_a_key=unstable, manifold_b_key=b_s)

    r_n, x1, r_end = point(1.0, 6.0), point(2.0, 3.0), point(4.0, 1.0)
    reg.register_iterate(r_n, 1, r_end)
    if unstable_owner == "A":
        point(2.002 * 4.0, 0.4985)
    else:
        point(2.002 * 6.0, 0.4985, unstable=b_u)

    pairs = _pairs(_trellis(reg, a, b))

    assert (r_n, x1) in pairs
    assert ((x1, r_end) in pairs) is survives
