"""Compute-Pseudoneighbors: reference pairs, the all-iterates interval check,
and trajectory extension.

The synthetic tests hand-build a registry (explicit canonical distances and
branch keys) plus the TrellisBranch ordering, so every geometric case is pinned
exactly. The structural validity of the reference pairs on real tangles is the
law ``reference_pairs_valid`` of ``tests/invariants/test_law_partition.py``.
"""

from __future__ import annotations

import inspect
import logging

import pytest

from helpers.logs import assert_logged
from helpers.fakes import bare_fixed_point
from tanglepack.numerics.FixedPoint import FixedPoint
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry
from tanglepack.topology.Pseudoneighbor import (
    compute_pseudoneighbors,
    extend_pseudoneighbor_trajectories,
)
from tanglepack.topology.TopologyResults import PseudoneighborPair
from tanglepack.topology.Trellis import Trellis
from tanglepack.topology.TrellisBranch import TrellisBranch

#: The default slack of the 2D endpoint-collision fallback; the table-linked
#: tests drift a landing by TWICE this, so the unlinked fallback must reject it.
COLLISION_RTOL = inspect.signature(compute_pseudoneighbors).parameters["collision_rtol"].default


def _fixed_point(period: int, lambda_u: float) -> FixedPoint:
    """A bare no-inversion fixed point whose full-cycle unstable eigenvalue is ``lambda_u``."""
    return bare_fixed_point(period, beta=lambda_u ** (-1.0 / period))


def _trellis(registry: IntersectionRegistry, *fixed_points: FixedPoint) -> Trellis:
    """A Trellis whose branches are bucketed and ordered from the registry,
    mirroring what Trellis.from_workbench derives from a workbench."""
    branches: dict = {}
    for ix_id, ix in registry:
        for key in (ix.manifold_a_key, ix.manifold_b_key):
            if key is None or key[0] not in fixed_points:
                continue
            branch = branches.setdefault(
                key,
                TrellisBranch(
                    key=key,
                    fixed_point=key[0],
                    stability=key[1],
                    orbit_index=key[2],
                    branch_index=key[3],
                    intersection_ids=[],
                ),
            )
            branch.intersection_ids.append(ix_id)
    for key, branch in branches.items():
        attr = "unstable_cdist" if key[1] == "unstable" else "stable_cdist"
        branch.intersection_ids.sort(key=lambda i: getattr(registry[i], attr))
    return Trellis(
        fixed_points=list(fixed_points),
        registry=registry,
        branches=branches,
        bridges=[],
    )


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
        manifold_a_key=(fp, "unstable", 0, 0), manifold_b_key=(fp, "stable", 0, 0),
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
        manifold_a_key=(fp, "unstable", 0, 0), manifold_b_key=(fp, "stable", 0, 0),
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
    # No full forward chain in the table: r_end is found by the cdist fallback.

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


def test_walk_stops_at_reference_window_end():
    """Pairs below r_{n+p} are outside the reference window and not returned."""
    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
    r_n, x1, r_end = _reference_window(fp, reg)
    below = reg.add_synthetic(
        (0.5, 16.0), unstable_cdist=16.0, stable_cdist=0.5,
        manifold_a_key=(fp, "unstable", 0, 0), manifold_b_key=(fp, "stable", 0, 0),
    )

    pairs = compute_pseudoneighbors(_trellis(reg, fp))

    members = {i for p in pairs for i in p.as_tuple()}
    assert below not in members
    assert len(pairs) == 2


def test_reference_end_fallback_by_cdist_match():
    """With no iterate links, r_{n+p} is located by matching BOTH scaled cdists."""
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
        manifold_a_key=unstable, manifold_b_key=stable,
    )
    r_end = reg.add_synthetic(
        (1.0, 4.0), unstable_cdist=4.0, stable_cdist=1.0,
        manifold_a_key=unstable, manifold_b_key=stable,
    )

    pairs = compute_pseudoneighbors(_trellis(reg, fp))

    assert [(p.intersection_a, p.intersection_b) for p in pairs] == [
        (r_n, x1),
        (x1, r_end),
    ]


def test_heteroclinic_candidate_ignored():
    """A candidate whose unstable side belongs to a different fixed point can
    never iterate onto this tangle's unstable manifold, so it cannot puncture."""
    fp = _fixed_point(1, 4.0)
    other = _fixed_point(1, 6.0)
    reg = IntersectionRegistry()
    _reference_window(fp, reg)
    reg.add_synthetic(
        (0.5, 1.5), unstable_cdist=1.5, stable_cdist=0.5,
        manifold_a_key=(other, "unstable", 0, 0), manifold_b_key=(fp, "stable", 0, 0),
    )

    pairs = compute_pseudoneighbors(_trellis(reg, fp, other))

    assert len(pairs) == 2


def test_unknown_unstable_branch_is_conservative():
    """A candidate with no unstable key (iterated-bridge point) is tested at
    every cycle residue and still punctures the interval."""
    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
    r_n, x1, r_end = _reference_window(fp, reg)
    reg.add_synthetic(
        (0.5, 1.5), unstable_cdist=1.5, stable_cdist=0.5,
        manifold_a_key=None, manifold_b_key=(fp, "stable", 0, 0),
    )

    pairs = compute_pseudoneighbors(_trellis(reg, fp))

    assert [(p.intersection_a, p.intersection_b) for p in pairs] == [(x1, r_end)]


def test_reference_window_starts_at_chosen_strong_pip():
    """With a strong pip chosen, r_n is its cut point — not the outermost
    intersection — and crossings beyond the cut leave the check set N."""
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
    branches come from extension — via the cdist fallback when the iterate
    table has no links — never as independent references."""
    fp = _fixed_point(3, 8.0)  # k = 3, beta = 2
    reg = IntersectionRegistry()

    def _point(orbit, s, u):
        return reg.add_synthetic(
            (s, u), unstable_cdist=u, stable_cdist=s,
            manifold_a_key=(fp, "unstable", orbit, 0),
            manifold_b_key=(fp, "stable", orbit, 0),
        )

    # Branch 0 carries the globally outermost intersection: the window
    # [4 -> 0.5] with a middle point. No iterate links anywhere.
    a0, b0, e0 = _point(0, 4.0, 1.0), _point(0, 3.0, 2.0), _point(0, 0.5, 8.0)
    # Branch 1 holds the forward images (stable /beta, unstable *beta).
    a1, b1, e1 = _point(1, 2.0, 2.0), _point(1, 1.5, 4.0), _point(1, 0.25, 16.0)

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
    removed tail) — e.g. the preimages of blasted crossings."""
    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
    r_n, x1, r_end = _reference_window(fp, reg)
    stable = (fp, "stable", 0, 0)
    unstable = (fp, "unstable", 0, 0)
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
    "period, uncut_warnings", [(1, 0), (3, 1)], ids=["period_1_cut", "period_3_uncut"]
)
def test_strong_pip_cuts_warn_when_branches_are_left_uncut(caplog, period, uncut_warnings):
    """A period-3 anchor needs one cut per stable branch; with no iterate table
    only the pip itself cuts, and the other two branches are left running to
    the manifold end -- which must be reported. A period-1 anchor has one
    stable branch, which the pip alone cuts: nothing to report.

    The window ``(4, 1) -> (3, 2) -> M^k = (4 / lambda, lambda)`` sits on orbit
    point 0 with the pip at its outer end, so the reference search itself
    succeeds silently in both cases.
    """
    lambda_u = 4.0 if period == 1 else 8.0
    fp = _fixed_point(period, lambda_u)
    reg = IntersectionRegistry()

    def _point(stable_cdist, unstable_cdist):
        return reg.add_synthetic(
            (stable_cdist, unstable_cdist),
            unstable_cdist=unstable_cdist, stable_cdist=stable_cdist,
            manifold_a_key=(fp, "unstable", 0, 0), manifold_b_key=(fp, "stable", 0, 0),
        )

    pip = _point(4.0, 1.0)
    _point(3.0, 2.0)
    _point(4.0 / lambda_u, lambda_u)
    trellis = _trellis(reg, fp)
    trellis.strong_pip = pip
    assert fp.k_value == period

    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.Pseudoneighbor"):
        pairs = compute_pseudoneighbors(trellis)

    assert len(pairs) == 2
    assert_logged(
        caplog, logging.WARNING, "tanglepack.topology.Pseudoneighbor", count=uncut_warnings
    )


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
                s *= 1.0 + 2.0 * COLLISION_RTOL  # scaled stable 3.06 vs x1's 3: outside it
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
    stable = (fp, "stable", 0, 0)
    unstable = (fp, "unstable", 0, 0)

    def candidate(reg, link: bool):
        z = reg.add_synthetic(
            (0.5, 4.5), unstable_cdist=4.5, stable_cdist=0.5,
            manifold_a_key=unstable, manifold_b_key=stable,
        )
        z1 = reg.add_synthetic(
            (0.125, 18.0), unstable_cdist=18.0, stable_cdist=0.125,
            manifold_a_key=unstable, manifold_b_key=stable,
        )
        r = reg.add_synthetic(
            (0.03125, 48.0), unstable_cdist=48.0, stable_cdist=0.03125,
            manifold_a_key=unstable, manifold_b_key=stable,
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
