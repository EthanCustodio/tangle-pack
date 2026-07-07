"""Compute-Pseudoneighbors: reference pairs, the all-iterates interval check,
and trajectory extension.

The synthetic tests hand-build a registry (explicit canonical distances and
branch keys) plus the TrellisBranch ordering, so every geometric case is pinned
exactly; the Hénon test checks structural invariants on a real tangle without
hard-coding registry ids.
"""

from __future__ import annotations

from tanglepack.numerics.FixedPoint import FixedPoint
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry
from tanglepack.topology.Pseudoneighbor import (
    compute_pseudoneighbors,
    extend_pseudoneighbor_trajectories,
    forward_unstable_branch_cycle,
)
from tanglepack.topology.TopologyResults import PseudoneighborPair
from tanglepack.topology.Trellis import Trellis
from tanglepack.topology.TrellisBranch import TrellisBranch


def _fixed_point(period: int, lambda_u: float) -> FixedPoint:
    """A minimal no-inversion fixed point with a positive unstable eigenvalue."""
    fp = FixedPoint(period, 1)
    fp.unstable_eigenvalues = [lambda_u] * period
    fp.set_k_value()
    return fp


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


def test_pair_on_different_unstable_branches_rejected():
    """Consecutive stable-branch points on two different unstable branches are
    never a pair — no single unstable arc connects them."""
    fp = _fixed_point(3, 8.0)
    reg = IntersectionRegistry()
    stable = (fp, "stable", 0, 0)
    r_n = reg.add_synthetic(
        (4.0, 1.0), unstable_cdist=1.0, stable_cdist=4.0,
        manifold_a_key=(fp, "unstable", 0, 0), manifold_b_key=stable,
    )
    x1 = reg.add_synthetic(
        (3.0, 2.0), unstable_cdist=2.0, stable_cdist=3.0,
        manifold_a_key=(fp, "unstable", 1, 0), manifold_b_key=stable,
    )
    r_end = reg.add_synthetic(
        (2.0, 4.0), unstable_cdist=2.0, stable_cdist=2.0,
        manifold_a_key=(fp, "unstable", 0, 0), manifold_b_key=stable,
    )
    reg.register_iterate(r_n, 1, reg.add_synthetic(
        (3.5, 1.5), unstable_cdist=1.26, stable_cdist=2.0 * (2.0 / 2.0),
        manifold_a_key=(fp, "unstable", 1, 0), manifold_b_key=(fp, "stable", 1, 0),
    ))
    # Chase M^3(r_n) back onto this branch: register the remaining two links.
    mid = reg.iterate_table[r_n, 1]
    mid2 = reg.add_synthetic(
        (3.2, 1.8), unstable_cdist=1.59, stable_cdist=2.5,
        manifold_a_key=(fp, "unstable", 2, 0), manifold_b_key=(fp, "stable", 2, 0),
    )
    reg.register_iterate(mid, 1, mid2)
    reg.register_iterate(mid2, 1, r_end)

    pairs = compute_pseudoneighbors(_trellis(reg, fp))

    assert (r_n, x1) not in [(p.intersection_a, p.intersection_b) for p in pairs]


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


def test_trajectory_extension_forward_with_dedup():
    """References map forward through the iterate table into non-reference
    pairs, stopping at the end of a chain, without duplicates."""
    fp = _fixed_point(1, 4.0)
    reg = IntersectionRegistry()
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

    trellis = _trellis(reg, fp)
    ref = PseudoneighborPair(a0, b0, branch_key=stable, is_reference=True)
    out = extend_pseudoneighbor_trajectories(trellis, [ref])

    assert [(p.intersection_a, p.intersection_b) for p in out] == [
        (a1, b1),
        (a2, b2),
    ]
    assert all(not p.is_reference for p in out)
    assert all(p.branch_key == stable for p in out)


def test_forward_unstable_branch_cycle_matches_orbit_order():
    """The unstable cycle runs the orbit in M-forward order, length k_value."""
    fp = _fixed_point(3, 8.0)
    cycle = forward_unstable_branch_cycle(fp)
    assert cycle == [
        (fp, "unstable", 0, 0),
        (fp, "unstable", 1, 0),
        (fp, "unstable", 2, 0),
    ]


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


def test_henon_reference_pairs_are_structurally_valid(henon_tangle_with_bridges):
    """On a real tangle every reference pair is a consecutive intersection pair
    on BOTH its stable and unstable branch, inside the reference window."""
    workbench, fp = henon_tangle_with_bridges
    trellis = Trellis.from_workbench(workbench, fp)

    references = trellis.compute_pseudoneighbors(extend=False)
    assert references, "the k=10 horseshoe tangle should contain reference pairs"

    lambda_u = trellis.lambda_u(fp)
    stable_branch = trellis.branch((fp, "stable", 0, 0))
    s_max = trellis.intersection(stable_branch.ordered_ids()[-1]).stable_cdist

    for pair in references:
        a, b = pair.as_tuple()
        s_branch = trellis.branch_containing(a, "stable")
        u_branch = trellis.branch_containing(a, "unstable")
        s_ids, u_ids = s_branch.intersection_ids, u_branch.intersection_ids
        assert abs(s_ids.index(a) - s_ids.index(b)) == 1
        assert abs(u_ids.index(a) - u_ids.index(b)) == 1
        for iid in (a, b):
            s = trellis.intersection(iid).stable_cdist
            assert s >= s_max / lambda_u * (1 - 1e-6)
            assert s <= s_max * (1 + 1e-6)
