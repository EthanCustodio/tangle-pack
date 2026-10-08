"""
Regression (INV-3, C1): a stable branch with no crossing past its anchor.

The k=10 inversion saddle grown only 4 or 5 unstable steps registers no crossing
past the anchor on some stable branch. Its partition had a single boundary and
no element at all, and ``_element_of_intersection`` asserted "owned by 0
partition elements". Under C1 such a branch is ONE element ``[anchor, end]``,
closed at the anchor and running to the end of the computed branch, and the
whole pipeline runs (the symbolic stage may stay provisional: grow further).
"""

from __future__ import annotations

import pytest

from cases import build_k10_inversion


@pytest.mark.parametrize("unstable_steps", [4, 5])
def test_an_anchor_only_branch_is_one_element(unstable_steps):
    case = build_k10_inversion(unstable_steps=unstable_steps)
    session = case.session
    session.symbolic_dynamics()

    trellis = session.trellis(case.fixed_point)
    bare = [
        branch
        for branch in trellis.stable_branches
        if all(trellis.intersection(i).stable_cdist == 0.0 for i in branch.intersection_ids)
    ]
    assert bare, "the shallow build must leave a stable branch with only its anchor"
    results = {(result.branch_key, result.side): result for result in trellis.stable_partitions}
    for branch in bare:
        for side in ("left", "right"):
            (element,) = results[branch.key, side].intervals
            assert element.lo_id == branch.intersection_ids[0]
            assert element.hi_id is None and element.closed_lo
            assert element.lo_cdist == 0.0 < element.hi_cdist
