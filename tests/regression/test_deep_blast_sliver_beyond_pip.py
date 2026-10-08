"""
Regression (2026-10-07, stress campaign E11 error 1): a sliver past the pip.

Nine blasts of the k=2.8 zone (analytic seed, ``area_cutoff = 1e-7``,
``min_separation = 1e-5``) register a real crossing 1.1e-6 past the zone's pip,
so the sliver's midpoint lies within the registry's absolute ``cdist_tol`` of
the pip. The iterated partition used to find each element's homotopy parent by
TOLERANT midpoint ownership, so the sliver ``(pip, outer]`` was handed to the
closed pip singleton and the child-inside-parent assertion fired. Ownership of registry cdists is exact now.
No count of the build is pinned.
"""

from __future__ import annotations

from cases import BUILDERS
from helpers import laws


def test_deep_blast_sliver_beyond_pip():
    case = BUILDERS["k28_deep_blast_sliver"]()
    session = case.session
    trellis = session.trellis(case.fixed_points[0])
    tol = float(trellis.registry.cdist_tol)
    pip = trellis.intersection(case.pips[0])
    branch = trellis.branch(pip.manifold_b_key)
    pip_cdist = float(pip.stable_cdist)
    beyond = [
        float(trellis.intersection(iid).stable_cdist) - pip_cdist
        for iid in branch.intersection_ids
        if float(trellis.intersection(iid).stable_cdist) > pip_cdist
    ]
    # The build must still exercise the defect: a sliver past the pip whose
    # midpoint lies within the tolerance of it.
    assert any(gap < 2 * tol for gap in beyond), "the builder no longer makes the sliver"

    session.symbolic_dynamics()
    assert laws.check_iterated_child_inside_parent(case) > 0
    assert laws.check_iterated_unique_owner(case) > 0
    assert laws.check_partition_unique_owner(case) > 0
