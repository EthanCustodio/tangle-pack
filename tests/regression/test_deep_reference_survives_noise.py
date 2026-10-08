"""
Regression (BCG-5, WP8): cdist noise never removes a reference pair.

At ten outer blasts the nested case lost its only outer reference pair and
the class count fell (13 -> 9 on the base): a crossing that lies OUTSIDE the
image of the pair's unstable arc on its own polyline was scaled back eight
steps, and the scaling noise there (2e-6 relative) exceeded the gap, so it
landed a hair inside the pair. A landing is now read against the pair's
registered images. Rule-based: every reference lost between the two
recomputes has a witness the iterate table names, consistent in every frame
it is registered in (an order flip made by cdist noise -- 725 = M(439)
between M(37) and M(38) while 439 lies outside (37, 38) -- is no witness),
and the class count does not decrease (D1); no id or count is pinned.
"""

from __future__ import annotations

from cases import build_nested
from helpers.laws import reference_witness

#: Outer blasts before the recompute pair compared here: the loss came at 10.
_OUTER_BLASTS = 9


def _outer_blast(case, zone) -> None:
    """One outer blast, the henon_cases recipe (its own fixed point, 1e-4 apart),
    without the same-stability report (it reports, it never changes the state)."""
    case.session.blast_zone(
        zone,
        num_iterations=1,
        fixed_point=[zone.fixed_point],
        min_separation=1e-4,
        diagnose=False,
    )


def _references(case) -> dict:
    """The reference pairs of each fixed point, after a recompute."""
    session = case.session
    session.topology()
    return {
        fp: {pair.as_tuple() for pair in session.trellis(fp).reference_pseudoneighbors}
        for fp in case.fixed_points
    }


def test_a_deep_blast_loses_no_reference_to_noise():
    case = build_nested(outer_blasts=0, through="zones")
    outer = max(case.zones, key=lambda zone: zone.area)
    for _ in range(_OUTER_BLASTS):
        _outer_blast(case, outer)
    before, classes_before = _references(case), len(case.session.bridge_classes())

    _outer_blast(case, outer)
    after, classes_after = _references(case), len(case.session.bridge_classes())

    for fp in case.fixed_points:
        trellis = case.session.trellis(fp)
        for pair in before[fp] - after[fp]:
            assert reference_witness(trellis, pair) is not None, (fp, pair)
    assert classes_before <= classes_after
