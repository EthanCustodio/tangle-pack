"""Characterization of the refiner's cdist rule.

This is a *characterization* test, not a correctness gate: it documents that
``_get_refined_point`` currently assigns the new point the arithmetic mean of its
neighbours' cdists. The mean rule is itself suspect (the geometric midpoint of the
two pre-iterates need not have the mean canonical distance) and may change when
the wavy-lobe bug is addressed -- at which point this test should be updated to
match the corrected rule.
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


def test_refined_cdist_is_mean_of_neighbours(grown_unstable):
    workbench, fp, manifold = grown_unstable
    machine = workbench._man_machine
    viewer = ManifoldView(manifold, machine.system)

    pair = _adjacent_pair_with_preiterates(manifold)
    assert pair is not None, "no suitable adjacent pair found"
    p0, p1 = pair

    new_point = machine._get_refined_point(p0, p1, viewer, "unstable")
    expected = 0.5 * (float(p0.cdist) + float(p1.cdist))
    assert np.isclose(float(new_point.cdist), expected, rtol=1e-12), (
        f"refined cdist {new_point.cdist} != mean of neighbours {expected}"
    )
