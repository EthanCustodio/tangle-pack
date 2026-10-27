"""The refiner's cdist rule: a refined point lies between its neighbours.

Only the FIRM part of the rule is tested: the new point between two adjacent
nodes gets a canonical distance strictly between theirs, so refinement never
breaks the ordering along the curve.

Dev Notes:
    The current rule assigns the arithmetic MEAN of the neighbours' cdists.
    That is provisional (the geometric midpoint of the two pre-iterates need
    not have the mean canonical distance, and it may change when the wavy-lobe
    bug is addressed), so it is deliberately not pinned; the old
    ``test_refined_cdist_is_mean_of_neighbours`` was replaced by this test in
    the 2026-10-05 refactor.
"""

from __future__ import annotations

import numpy as np

from tanglepack import ManifoldView


def _adjacent_pair_with_preiterates(manifold):
    """Find adjacent rows (p0, p1) that both carry a pre-iterate and a cdist."""
    store, rows = manifold.store, manifold.ids.tolist()
    for p0, p1 in zip(rows, rows[1:]):
        if np.isnan(store.cdist[p0]) or np.isnan(store.cdist[p1]):
            continue
        if store.preimage[p0] < 0 or store.preimage[p1] < 0:
            continue
        return p0, p1
    return None


def test_refined_cdist_lies_strictly_between_its_neighbours(grown_unstable):
    """The batched refinement kernel (``_refine_layer``, no public single-pair
    route; a pair is its first point's chain position) places the new point
    between the pair on the chain, its cdist strictly inside theirs."""
    workbench, fp, manifold = grown_unstable
    machine = workbench._man_machine
    viewer = ManifoldView(manifold, machine.system)

    pair = _adjacent_pair_with_preiterates(manifold)
    assert pair is not None, "no suitable adjacent pair found"
    p0, p1 = pair
    store = manifold.store
    lo, hi = sorted((float(store.cdist[p0]), float(store.cdist[p1])))
    assert lo < hi

    chain = store.chain_of(p0)
    position = store.position_of(p0)
    assert chain.ids[position + 1] == p1

    machine.area_cutoff = 0.0  # refine this pair whatever its curvature area
    next_layer = machine._refine_layer(chain, np.array([position]), manifold, viewer)
    np.testing.assert_array_equal(next_layer, [position, position + 1])
    new_point = chain.ids[position + 1]
    assert chain.ids[position] == p0 and chain.ids[position + 2] == p1

    assert lo < float(store.cdist[new_point]) < hi
