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
from helpers.polylines import polyline_manifold
from tanglepack.numerics.Tangle import Tangle


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
    first = polyline_manifold([(0.0, 0.0), (1.0, 0.001), (2.0, 0.0)], (fp, stability, 0, 0))
    second = polyline_manifold([(0.0, 0.0005), (1.0, 0.0), (2.0, 0.0005)], (fp, stability, 0, 1))
    transversal = polyline_manifold([(0.5, -1.0), (0.5, 1.0)], (fp, other, 0, 0))

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
