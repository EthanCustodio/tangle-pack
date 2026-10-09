"""
Regression (2026-10-07, stress campaign E11 error 1): a sliver past the pip.

Nine blasts of the k=2.8 zone (analytic seed, ``area_cutoff = 1e-7``,
``min_separation = 1e-5``) register a real crossing 1.1e-6 past the zone's pip,
so the sliver's midpoint lies within the registry's absolute ``cdist_tol`` of
the pip. The iterated partition used to find each element's homotopy parent by
TOLERANT midpoint ownership, so the sliver ``(pip, outer]`` was handed to the
closed pip singleton and the child-inside-parent assertion fired. Ownership of
registry cdists is exact now, and since D10 (WP8b) the trellis is cut at the
pip, so the sliver crossing is no crossing of the trellis at all. No count of
the build is pinned.

The same build had two ambiguous classes and a walk contradicted by its
registered image until the walk left from the edges its landing SPAN meets
rather than from every edge of the iterated element (D13, WP9b).
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
    pip_cdist = float(pip.stable_cdist)
    beyond = {
        crossing_id: float(ix.stable_cdist) - pip_cdist
        for crossing_id, ix in trellis.registry
        if ix.manifold_b_key == pip.manifold_b_key and float(ix.stable_cdist) > pip_cdist
    }
    # The build must still exercise the defect: a computed crossing a sliver past
    # the pip whose midpoint lies within the tolerance of it ...
    slivers = [crossing_id for crossing_id, gap in beyond.items() if gap < 2 * tol]
    assert slivers, "the builder no longer makes the sliver"
    # ... which the cut at the pip orbit leaves out of the trellis (D10).
    assert not [crossing_id for crossing_id in slivers if crossing_id in trellis]

    dynamics = session.symbolic_dynamics()
    assert laws.check_iterated_child_inside_parent(case) > 0
    assert laws.check_iterated_unique_owner(case) > 0
    assert laws.check_partition_unique_owner(case) > 0

    assert laws.check_walk_starts_on_the_landing(case) > 0
    assert not [cd.letter for cd in dynamics.classes.values() if cd.ambiguous]
    assert laws.check_walks_agree_with_registered_images(case) > 0
