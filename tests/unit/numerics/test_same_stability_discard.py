"""A same-stability segment pair is never a crossing (CLAUDE.md fundamental invariant).

Two unstable manifolds (or two stable ones) can never cross, so a u x u or s x s
pair the spatial index turns up is a near-tangency of two polygonal
approximations: ``Tangle.resolve_crossings`` drops it and logs it at DEBUG on
the ``tanglepack.numerics.Tangle`` logger, while an unstable x stable pair over
the same picture still resolves. The geometric law that no registered crossing
is same-stability runs on every case in ``tests/invariants/test_law_crossings.py``.
"""

from __future__ import annotations

import logging

import pytest

from helpers.fakes import bare_fixed_point
from helpers.logs import assert_logged
from tanglepack import BaseManifold, Point
from tanglepack.numerics.Tangle import Tangle

LOGGER = "tanglepack.numerics.Tangle"


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


@pytest.mark.parametrize("stability", ["unstable", "stable"])
def test_a_same_stability_crossing_is_dropped_and_logged(
    stability: str, caplog: pytest.LogCaptureFixture
) -> None:
    """Two near-tangent same-stability curves straddling each other give no
    crossing and one DEBUG record; the curve of the other stability through the
    same picture still crosses both."""
    fp = bare_fixed_point(1)
    other = "stable" if stability == "unstable" else "unstable"
    # Two near-parallel polylines that straddle once (a polygonal artifact of a
    # near-tangency), and a transversal curve of the other stability.
    first = _polyline([(0.0, 0.0), (1.0, 0.001), (2.0, 0.0)], (fp, stability, 0, 0))
    second = _polyline([(0.0, 0.0005), (1.0, 0.0), (2.0, 0.0005)], (fp, stability, 0, 1))
    transversal = _polyline([(0.5, -1.0), (0.5, 1.0)], (fp, other, 0, 0))

    tangle = Tangle()
    tangle.add_manifolds([first, second, transversal])
    with caplog.at_level(logging.DEBUG, logger=LOGGER):
        crossings = tangle.resolve_crossings()

    assert_logged(caplog, logging.DEBUG, LOGGER)
    assert len(crossings) == 2
    for crossing in crossings:
        stabilities = {crossing.manifold_a_key[1], crossing.manifold_b_key[1]}
        assert stabilities == {"unstable", "stable"}
