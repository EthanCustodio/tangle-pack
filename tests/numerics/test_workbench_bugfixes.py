"""Regression tests for the Phase 1 TangleWorkbench bug fixes.

Covers plan rows 1.2 (bridge endpoint assignment must stay on the bridge's own
unstable branch), 1.11 (``trim_stable_manifolds`` on a branch with no crossings)
and 1.12 (growth-loop iteration caps and the ``grow_until_intersection`` rename).
"""

from __future__ import annotations

import pytest


# --------------------------------------------------------------------------- #
# 1.2 -- _assign_bridge_intersections must not steal endpoints from another branch
# --------------------------------------------------------------------------- #
def test_bridge_endpoints_stay_on_the_bridges_own_branch(henon_p3_session):
    """On a period-3 tangle every unstable branch has its own crossing at
    unstable cdist ~0, so a registry-wide nearest-cdist lookup hands the first
    bridge of one branch an endpoint belonging to another."""
    session, _fp3, _fp1, _zone = henon_p3_session
    workbench = session.workbench
    registry = workbench._intersection_registry

    keyed = [b for b in workbench.bridges if b.manifold_key is not None]
    assert keyed, "expected at least one keyed bridge on the period-3 tangle"

    for bridge in keyed:
        first_id = bridge.first_intersection
        second_id = bridge.second_intersection
        assert first_id is not None and second_id is not None
        assert first_id in registry and second_id in registry
        assert first_id != second_id

        first = registry[first_id]
        second = registry[second_id]

        for endpoint in (first, second):
            assert endpoint.manifold_a_key in (bridge.manifold_key, None), (
                f"bridge {bridge.manifold_key} got an endpoint from "
                f"{endpoint.manifold_a_key}"
            )

        root_u = bridge.root.get_cdist("unstable")
        tail_u = bridge.tail.get_cdist("unstable")
        assert (
            root_u <= first.unstable_cdist < second.unstable_cdist <= tail_u
        ), (
            f"bridge {bridge.manifold_key} spans [{root_u}, {tail_u}] but was "
            f"given endpoints at {first.unstable_cdist}, {second.unstable_cdist}"
        )


# --------------------------------------------------------------------------- #
# 1.11 -- trim_stable_manifolds
# --------------------------------------------------------------------------- #
def test_trim_stable_manifolds_without_crossings(fixed_point):
    """A stable branch that nothing crosses is skipped, not a ``max()`` on [].

    Note:
        The brief suggested "unstable grown twice, no crossings"; on the k=10
        saddle the unstable and stable manifolds always meet at the anchor
        (unstable cdist 0), so that configuration still has one crossing. A
        workbench with only the stable manifold initialized is the deterministic
        no-crossing case.
    """
    workbench, fp = fixed_point
    workbench.initialize_manifold(fp, "stable")
    workbench.grow_until_turnaround(fp, "stable")
    workbench.compute_intersections([fp], infer_iterates=False)

    assert len(workbench._intersection_registry) == 0

    manifolds = list(workbench._iter_manifolds(fp, "stable"))
    assert manifolds
    tails_before = [m.tail for m in manifolds]

    workbench.trim_stable_manifolds(fp)

    assert [m.tail for m in manifolds] == tails_before


def test_trim_stable_manifolds_reads_the_stable_cdist(small_tangle):
    """``BranchPoint.get_cdist`` has no default stability, so a crossing segment
    anchored at the branch point used to raise ``TypeError``."""
    workbench, fp = small_tangle

    crossing_seg_ids = {
        n for pair in workbench.Tangle._intersecting_segments for n in pair
    }
    manifold = workbench.manifolds[(fp, "stable", 0, 0)]
    on_this = sorted(crossing_seg_ids & workbench.Tangle._manifold_segs[manifold])
    assert len(on_this) >= 2

    innermost = min(
        on_this, key=lambda i: workbench.Tangle._seg_lookup[i].p0_seg1.cdist
    )
    expected = max(
        (workbench.Tangle._seg_lookup[i] for i in on_this if i != innermost),
        key=lambda s: s.p0_seg1.get_cdist("stable"),
    ).p0_seg1

    # anchor the innermost crossing segment at the fixed point's BranchPoint
    workbench.Tangle._seg_lookup[innermost].p0_seg1 = manifold.root

    workbench.trim_stable_manifolds(fp)

    assert manifold.tail is expected


def test_trim_stable_manifolds_cuts_just_past_the_outermost_crossing(small_tangle):
    """The new tail is the low node of the outermost crossing segment, and no
    crossing is left beyond it."""
    workbench, fp = small_tangle

    crossing_seg_ids = {
        n for pair in workbench.Tangle._intersecting_segments for n in pair
    }
    assert crossing_seg_ids

    for manifold in workbench._iter_manifolds(fp, "stable"):
        on_this = crossing_seg_ids & workbench.Tangle._manifold_segs[manifold]
        if not on_this:
            continue
        expected = max(
            (workbench.Tangle._seg_lookup[i] for i in on_this),
            key=lambda s: s.p0_seg1.get_cdist("stable"),
        ).p0_seg1

        workbench.trim_stable_manifolds(fp)

        assert manifold.tail is expected
        outermost = max(
            ix.stable_cdist for _id, ix in workbench._intersection_registry
        )
        assert manifold.tail.get_cdist("stable") >= outermost
        return

    pytest.fail("no stable manifold carried a crossing")


# --------------------------------------------------------------------------- #
# 1.12 -- growth loop caps and the grow_until_intersection rename
# --------------------------------------------------------------------------- #
def test_grow_until_intersection_stops_once_a_crossing_exists(initialized):
    workbench, fp = initialized
    workbench.grow_until_turnaround(fp, "stable")

    workbench.grow_until_intersection(fp, "unstable", max_iterations=10)

    assert len(workbench.Tangle._intersecting_segments) > 0


def test_grow_until_intersection_raises_at_the_cap(initialized):
    workbench, fp = initialized

    with pytest.raises(ValueError, match="Max iterations"):
        workbench.grow_until_intersection(fp, "unstable", max_iterations=1)


def test_grow_until_intersection_honours_branch_index(initialized):
    """branch_index must reach the manifold lookup: an uninitialized branch is
    reported as such rather than silently growing branch 0."""
    workbench, fp = initialized

    with pytest.raises(ValueError, match="has not been initialized"):
        workbench.grow_until_intersection(fp, "unstable", branch_index=1)


def test_grow_until_arclength_raises_at_the_cap(initialized):
    workbench, fp = initialized

    with pytest.raises(ValueError, match="Max iterations"):
        workbench.grow_until_arclength(
            fp, "unstable", length=1e9, max_iterations=2
        )


def test_grow_until_arclength_returns_when_long_enough(initialized):
    workbench, fp = initialized

    workbench.grow_until_arclength(fp, "unstable", length=1.0, max_iterations=10)

    tail = workbench.manifolds[(fp, "unstable", 0, 0)].tail
    assert tail.get_cdist("unstable") >= 1.0


def test_grown_until_intersection_is_gone(initialized):
    workbench, _fp = initialized
    assert not hasattr(workbench, "grown_until_intersection")


# --------------------------------------------------------------------------- #
# 1.7 -- the workbench key advance delegates to FixedPoint.advance_key
# --------------------------------------------------------------------------- #
def test_advance_key_forward_returns_the_next_orbit_point(henon_p3_session):
    """One map step advances the orbit index and leaves the branch alone (the
    period-3 orbit has no inversion, so ``k_value == period``)."""
    session, fp3, fp1, _zone = henon_p3_session
    workbench = session.workbench

    assert fp3.period == 3 and fp3.k_value == 3

    assert workbench._advance_key_forward((fp3, "unstable", 2, 0), fp3) == (
        fp3,
        "unstable",
        0,
        0,
    )
    assert workbench._advance_key_forward((fp3, "stable", 0, 0), fp3) == (
        fp3,
        "stable",
        1,
        0,
    )
    assert workbench._advance_key_forward((fp1, "unstable", 0, 0), fp1) == (
        fp1,
        "unstable",
        0,
        0,
    )
