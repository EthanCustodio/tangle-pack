"""A same-stability segment pair is never formed (CLAUDE.md fundamental invariant).

Two unstable manifolds (or two stable ones) can never cross, so the Tangle keeps
one rtree per stability and queries a segment only against the OTHER one
(WP10): a u x u or s x s pair is never formed, however the two polylines
straddle, while an unstable x stable pair over the same picture still resolves.
The geometric law that no registered crossing is same-stability runs on every
case in ``tests/invariants/test_law_crossings.py``.
"""

from __future__ import annotations

import pytest

from helpers.fakes import bare_fixed_point
from tanglepack import BaseManifold, Point
from tanglepack.numerics.Tangle import Tangle


def _polyline(points: list[tuple[float, float]], key: tuple) -> BaseManifold:
    """
    A manifold through ``points``, cdist = the index, keyed by ``key``.

    The walk leaves the fixed point along ``forward`` links on an unstable
    manifold and along ``backward`` links on a stable one, so the nodes are
    linked in the direction the key's stability walks.

    Args:
        points: The polyline's vertices in order.
        key: The manifold key ``(fixed_point, stability, orbit, branch)``.

    Returns:
        The manifold, root at the first point and tail at the last.
    """
    fixed_point, stability, _orbit, branch = key
    nodes = [Point(x, y, float(i)) for i, (x, y) in enumerate(points)]
    for node, following in zip(nodes, nodes[1:]):
        if stability == "unstable":
            node.insert_point_forward(following)
        else:
            node.insert_point_backward(following)
    return BaseManifold(
        nodes[0],
        stability,
        1.0,
        fixed_point=fixed_point,
        tail=nodes[-1],
        branch_index=branch,
        manifold_key=key,
    )


@pytest.mark.parametrize("bulk", [True, False], ids=["bulk", "incremental"])
@pytest.mark.parametrize("stability", ["unstable", "stable"])
def test_a_same_stability_pair_is_never_formed(stability: str, bulk: bool) -> None:
    """Two near-tangent same-stability curves straddling each other form no
    segment pair, in a bulk load or added one at a time; the curve of the other
    stability through the same picture still crosses both."""
    fp = bare_fixed_point(1)
    other = "stable" if stability == "unstable" else "unstable"
    # Two near-parallel polylines that straddle once (a polygonal artifact of a
    # near-tangency), and a transversal curve of the other stability.
    first = _polyline([(0.0, 0.0), (1.0, 0.001), (2.0, 0.0)], (fp, stability, 0, 0))
    second = _polyline([(0.0, 0.0005), (1.0, 0.0), (2.0, 0.0005)], (fp, stability, 0, 1))
    transversal = _polyline([(0.5, -1.0), (0.5, 1.0)], (fp, other, 0, 0))

    tangle = Tangle()
    if bulk:
        tangle.add_manifolds([first, second, transversal])
    else:
        for manifold in (first, second, transversal):
            tangle.add_manifold(manifold)

    stability_of = {sid: seg.manifold.stability for sid, seg in tangle._seg_lookup.items()}
    assert all(len({stability_of[sid] for sid in pair}) == 2 for pair in tangle._intersecting_segments)
    crossings = tangle.resolve_crossings()
    assert len(crossings) == 2
    for crossing in crossings:
        assert {crossing.manifold_a_key[1], crossing.manifold_b_key[1]} == {"unstable", "stable"}
