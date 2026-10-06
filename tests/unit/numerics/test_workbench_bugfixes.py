"""Regression tests for the Phase 1 TangleWorkbench bug fixes.

Covers plan rows 1.11 (``trim_stable_manifolds`` on a branch with no crossings)
and 1.12 (growth-loop iteration caps and the ``grow_until_intersection`` rename).
Row 1.2 (a bridge's endpoints stay on its own unstable branch) is the law
``bridge_endpoints_on_own_branch`` of ``tests/invariants/test_law_bridges.py``,
run on the period-3 and nested cases where it once broke.
"""

from __future__ import annotations



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

    assert len(workbench.intersection_registry) == 0

    manifolds = [
        m for key, m in workbench.manifolds.items() if key[0] is fp and key[1] == "stable"
    ]
    assert manifolds
    tails_before = [m.tail for m in manifolds]

    workbench.trim_stable_manifolds(fp)

    assert [m.tail for m in manifolds] == tails_before


# --------------------------------------------------------------------------- #
# 1.12 -- growth loop caps and the grow_until_intersection rename
# --------------------------------------------------------------------------- #
def test_grow_until_intersection_stops_once_a_crossing_exists(initialized):
    """The driver returns as soon as the registry holds a crossing it did not before."""
    workbench, fp = initialized
    workbench.grow_until_turnaround(fp, "stable")

    workbench.grow_until_intersection(fp, "unstable", max_iterations=10)

    registry = workbench.intersection_registry
    assert any(not ix.is_synthetic for _iid, ix in registry), (
        "the driver returned without a detected (non-anchor) crossing"
    )


def test_grow_until_arclength_returns_when_long_enough(initialized):
    workbench, fp = initialized

    workbench.grow_until_arclength(fp, "unstable", length=1.0, max_iterations=10)

    tail = workbench.manifolds[(fp, "unstable", 0, 0)].tail
    assert tail.get_cdist("unstable") >= 1.0
