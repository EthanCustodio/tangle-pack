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
    """Find an adjacent (p0, p1) where both carry a real pre-iterate and cdist."""
    nodes = manifold.get_point_array(return_nodes=True)
    for p0, p1 in zip(nodes, nodes[1:]):
        if getattr(p0, "cdist", None) is None or getattr(p1, "cdist", None) is None:
            continue
        if p0.prev_iterate is None or p1.prev_iterate is None:
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
    lo, hi = sorted((float(p0.cdist), float(p1.cdist)))
    assert lo < hi

    store = p0.store
    chain = store.chain_of(p0.row)
    position = store.position_of(p0.row)
    assert chain.node_at(position + 1) is p1

    machine.area_cutoff = 0.0  # refine this pair whatever its curvature area
    next_layer = machine._refine_layer(chain, np.array([position]), manifold, viewer, None)
    np.testing.assert_array_equal(next_layer, [position, position + 1])
    new_point = chain.node_at(position + 1)
    assert chain.node_at(position) is p0 and chain.node_at(position + 2) is p1

    assert lo < float(new_point.cdist) < hi
