"""
Regression (WP9, NEST-1 / NEST-2, D7): a heteroclinic pair is a pseudoneighbor pair
of the tangle whose stable branch it sits on.

The ``heteroclinic`` case (period 1 + period 3, period 3 grown 16 unstable steps)
has a thin lobe of the inner orbit's unstable manifold A poking across the outer
saddle's stable manifold B inside B's fundamental window. Its two crossings are
consecutive on B's stable branch and on A's unstable branch, so they are a
reference pair of B: tested against A's crossings with A's branch cycle (the arc
is A's), holed in B's partition, and its backward chain runs along A's unstable
cycle into A's bridges. Such a pair exists only when both tangles are selected.
"""

from __future__ import annotations

import logging

from cases import BUILDERS


def test_heteroclinic_pairs_hole_the_stable_owners_partition(caplog):
    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.Pseudoneighbor"):
        case = BUILDERS["heteroclinic"]()
        trellis = case.session.topology().trellis
    assert not [r for r in caplog.records if "missing from the branch cycle" in r.message]

    def unstable_owner(pair):
        return trellis.intersection(pair.as_tuple()[0]).manifold_a_key[0]

    heteroclinic = [
        pair for pair in trellis.reference_pseudoneighbors
        if unstable_owner(pair) is not pair.branch_key[0]
    ]
    assert heteroclinic, "no heteroclinic reference pair"
    for pair in heteroclinic:
        owner, stable_owner = unstable_owner(pair), pair.branch_key[0]
        assert pair.hole is not None
        assert {opening[3] for opening in pair.hole.openings} == {pair.branch_key}
        backward = [h for h in trellis.holes if h.origin == pair.as_tuple() and h.iterate < 0]
        assert backward, "the heteroclinic reference propagates no hole"
        for hole in backward:
            assert trellis.bridge_between(*hole.bounding_ids).fixed_point is owner

    alone = case.session.topology([stable_owner]).trellis
    assert all(
        trellis.intersection(p.as_tuple()[0]).manifold_a_key[0] is stable_owner
        for p in alone.pseudoneighbors
    )
