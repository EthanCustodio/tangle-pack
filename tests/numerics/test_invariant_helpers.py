"""Exercise the shared invariant helpers on real manifolds and a real registry.

Two of the checks in ``tests/invariants.py`` had no consumer, so nothing pinned
the invariants they encode:

* ``assert_no_cdist_collision`` — on low-stretch growth the canonical distance is
  injective, so no two distinct nodes may share a cdist. (Ties are legitimate only
  at a high-stretch fold, which the k=10 binary horseshoe never reaches at these
  depths; see ``tests/regression/test_high_stretch_period3_growth.py``.)
* ``assert_area_preserved_along_chain`` — the map is area preserving, so one
  forward step scales the unstable canonical distance up and the stable one down
  by the same per-step factor and the product ``unstable_cdist * stable_cdist`` is
  invariant *along one iterate chain*.

Note:
    Per CLAUDE.md the canonical-area product is only invariant *within* a chain;
    two different chains generally carry different products, so the product must
    never be used to decide chain membership. These tests therefore only walk
    links the iterate table already records.
"""

from __future__ import annotations

import pytest

from invariants import assert_area_preserved_along_chain, assert_no_cdist_collision


def _manifold(workbench, fp, stability: str, branch_index: int = 0):
    return workbench.manifolds[(fp, stability, 0, branch_index)]


def test_fundamental_segments_have_injective_cdist(initialized):
    """Both fundamental segments come out of the initializer collision-free."""
    workbench, fp = initialized
    for stability in ("unstable", "stable"):
        manifold = _manifold(workbench, fp, stability)
        nodes = manifold.get_point_array(return_nodes=True)
        assert len(nodes) >= 3, "fundamental segment is degenerate"
        assert_no_cdist_collision(manifold)


def test_low_stretch_growth_keeps_cdist_injective(initialized):
    """Two growth iterations of the k=10 horseshoe keep cdist injective.

    The k=10 binary horseshoe is low stretch: after two iterations the tightest
    cdist gap is still ~1e-7, orders of magnitude above the 1e-12 collision
    tolerance, so any tie here is a merge/refine bug, not a genuine fold.
    """
    workbench, fp = initialized
    for stability in ("unstable", "stable"):
        workbench.grow_n_times(fp, stability, num_iterations=2)

    for stability in ("unstable", "stable"):
        manifold = _manifold(workbench, fp, stability)
        nodes = manifold.get_point_array(return_nodes=True)
        assert len(nodes) > 3, "growth added no points"
        assert_no_cdist_collision(manifold)


def test_area_preserved_along_every_recorded_iterate_chain(
    henon_tangle_with_bridges,
):
    """Every n=1 link in the iterate table preserves the canonical-area product.

    ``compute_intersections`` fills the table via ``infer_iterates``, which matches
    an image by predicted branch keys plus canonical distances (unstable stretched
    by the per-step factor, stable contracted by it). The prediction preserves the
    product exactly, so the *matched* crossing must too, up to the tolerance of the
    match itself (observed worst case here ~2e-4).
    """
    workbench, _fp = henon_tangle_with_bridges
    registry = workbench.intersection_registry

    chain_starts = [
        iid for iid in registry.all_ids() if registry.iterate_table[iid, 1] is not None
    ]
    assert chain_starts, "no n=1 iterate links recorded; the test would be vacuous"

    for start_id in chain_starts:
        assert_area_preserved_along_chain(registry, start_id, rtol=1e-3)


def test_area_preserved_helper_rejects_a_broken_chain(henon_tangle_with_bridges):
    """The helper actually fires: corrupting one image's cdists trips it."""
    workbench, _fp = henon_tangle_with_bridges
    registry = workbench.intersection_registry

    # Skip the anchors: a periodic point's crossing sits at cdist (0, 0) and its
    # image is the next orbit point's anchor, also at (0, 0), so its area product
    # is 0 -> 0 and scaling a zero cdist corrupts nothing. Anchors are registered
    # first (Phase 6.3) and so are the lowest ids, which is why the plain "first
    # id with an n=1 link" used to land on a real crossing and now would not.
    start_id = next(
        iid
        for iid in registry.all_ids()
        if registry.iterate_table[iid, 1] is not None
        and registry[iid].unstable_cdist * registry[iid].stable_cdist != 0.0
    )
    image = registry[registry.iterate_table[start_id, 1]]
    original = image.unstable_cdist
    image.unstable_cdist = original * 2.0
    try:
        with pytest.raises(AssertionError):
            assert_area_preserved_along_chain(registry, start_id, rtol=1e-3)
    finally:
        image.unstable_cdist = original
