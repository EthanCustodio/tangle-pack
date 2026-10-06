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
    """``_get_refined_point`` (the refinement kernel, no public single-point
    route) places the new point's cdist strictly inside its neighbours'."""
    workbench, fp, manifold = grown_unstable
    machine = workbench._man_machine
    viewer = ManifoldView(manifold, machine.system)

    pair = _adjacent_pair_with_preiterates(manifold)
    assert pair is not None, "no suitable adjacent pair found"
    p0, p1 = pair
    lo, hi = sorted((float(p0.cdist), float(p1.cdist)))
    assert lo < hi

    new_point = machine._get_refined_point(p0, p1, viewer, "unstable")

    assert lo < float(new_point.cdist) < hi
