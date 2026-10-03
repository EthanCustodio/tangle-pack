"""Period-3 regressions for the hole/partition layer (plan row 1.1).

A bridge's unstable branch is carried by ``Bridge.manifold_key``. Inferring it
from the endpoint intersections instead is unreliable on a period > 1 orbit:
every branch's anchor artifact sits at unstable cdist 0, so the after-the-fact
endpoint assignment can hand one branch's anchor to another branch's first
bridge, and the backward propagation then picks a container on the wrong
branch. The symptom is an orbit whose holes flip ``bridge_side`` partway down
the backward chain.

These tests run the real nested period-3 tangle (``henon_p3_session``), which is
the only fixture in the suite with a period > 1 orbit, and pin:

I1  all holes of one origin share their side of their bridge;
I2  a bridge whose two defining crossings lie on one stable branch approaches
    both from the same side of it;
III every propagated hole lands in a bridge on the unstable branch the
    branch cycle predicts from its iterate.

The deep runs (6 blasts, 15 unstable steps) pin I1 where a propagated hole's
image sub-arc and the nearest vertex of its whole containing bridge sit on
different folds (fixed 2026-10-02).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless: the session fixture touches the plotting stack
import numpy as np
import pytest

from tanglepack.examples.henon_cases import build_period3
from tanglepack.topology import (
    check_bridge_rows_consistent,
    check_holes_share_bridge_side,
)
from tanglepack.topology.Pseudoneighbor import forward_unstable_branch_cycle
from tanglepack.topology.StablePartition import _bridge_side_of, bridge_for_pair


@pytest.fixture
def p3_punched(henon_p3_session):
    """``(trellis, fp3)`` for the period-3 orbit with every hole punched."""
    session, fp3, _fp1, _zone = henon_p3_session
    trellis = session.trellis(fp3)
    trellis.classify_strong_pips()
    references = trellis.compute_pseudoneighbors()
    assert references, "the period-3 orbit must yield reference pairs"
    trellis.punch_holes()
    assert trellis.holes
    return trellis, fp3


def test_p3_holes_share_bridge_side(p3_punched):
    """I1: an orbit's backward chain never switches sides of its bridge.

    The period-3 Hénon map has b = 1 (det J = +1), so the side is carried
    unchanged to every backward image; a flip means the chain stepped into a
    container on the wrong unstable branch.
    """
    trellis, _fp3 = p3_punched
    check_holes_share_bridge_side(
        trellis.holes, orientation_preserving=trellis.orientation_preserving
    )


def test_p3_bridge_rows_consistent(p3_punched):
    """I2: every period-3 bridge meets its stable branch from one side."""
    trellis, _fp3 = p3_punched
    for bridge in trellis.bridges:
        check_bridge_rows_consistent(trellis, bridge)


def test_p3_propagated_holes_land_on_the_predicted_branch(p3_punched):
    """III: a hole's containing bridge sits on the cycle-predicted branch.

    One backward step moves the image to the PREVIOUS unstable branch of the
    forward cycle, so a hole at iterate ``n`` (negative) of the orbit whose
    reference bridge sits at cycle position ``pos_ref`` must land in a bridge
    keyed ``cycle[(pos_ref + n) % k]``.
    """
    trellis, fp3 = p3_punched
    cycle = forward_unstable_branch_cycle(fp3)
    k = len(cycle)
    assert k == 3

    reference_position: dict[tuple[int, int], int] = {}
    for pair in trellis.pseudoneighbors:
        if not pair.is_reference:
            continue
        bridge = bridge_for_pair(trellis, pair)
        assert bridge is not None and bridge.manifold_key in cycle
        reference_position[pair.as_tuple()] = cycle.index(bridge.manifold_key)
    assert reference_position

    bridges_by_ends = {
        frozenset((b.first_intersection, b.second_intersection)): b
        for b in trellis.bridges
    }

    checked = 0
    for hole in trellis.holes:
        if hole.pair is not None:
            continue  # a directly punched hole, not a propagated one
        assert hole.iterate is not None and hole.origin is not None
        pos_ref = reference_position[hole.origin]
        expected = cycle[(pos_ref + hole.iterate) % k]
        bridge = bridges_by_ends[frozenset(hole.bounding_ids)]
        assert bridge.manifold_key == expected, (
            f"hole of origin {hole.origin} at iterate {hole.iterate} landed in "
            f"bridge {hole.bounding_ids} on branch "
            f"{bridge.manifold_key[1:] if bridge.manifold_key else None}, "
            f"expected {expected[1:]}"
        )
        checked += 1
    assert checked, "backward propagation should have punched holes"


@pytest.mark.slow
@pytest.mark.regression
@pytest.mark.parametrize(
    "kwargs", [{"blasts": 6}, {"unstable_steps": 15}], ids=["6_blasts", "15_steps"]
)
def test_deep_p3_holes_share_bridge_side_through_symbolic_dynamics(kwargs):
    """I1 holds deep into the tangle, and the symbolic dynamics runs.

    Only that it runs: these tangles are not reliable yet (an unreachable
    class, a virtual ``new1``).
    """
    build = build_period3(**kwargs)
    (fp3,) = build.fixed_points
    trellis = build.session.trellis(fp3)
    check_holes_share_bridge_side(
        trellis.holes, orientation_preserving=trellis.orientation_preserving
    )
    assert build.session.symbolic_dynamics().classes


def _hole_sides(trellis) -> list:
    """Every hole as ``(iterate, bridge_side, sorted (which, row) openings)``, ids dropped."""
    return sorted(
        (hole.iterate, hole.bridge_side, sorted((which, row) for _id, which, row in hole.openings))
        for hole in trellis.holes
    )


_OPENS_LEFT = [("anchorward", "left"), ("outward", "left")]
_OPENS_RIGHT = [("anchorward", "right"), ("outward", "right")]


@pytest.mark.slow
def test_p3_hole_sides_are_pinned():
    """Sides and openings of every period-3 hole, as at 2315204 (2026-10-02).

    Direct holes (iterates 0..2) open their inward pair; each orbit's
    propagated holes keep its side. Ids are not reproducible, so none appear.
    """
    build = build_period3()
    (fp3,) = build.fixed_points
    expected = []
    for iterate in range(-3, 3):
        if iterate < 0:
            expected += [(iterate, "left", [("outward", "right")]), (iterate, "right", _OPENS_RIGHT)]
        else:
            expected += [(iterate, "left", _OPENS_RIGHT), (iterate, "right", _OPENS_LEFT)]
    assert _hole_sides(build.session.trellis(fp3)) == expected


@pytest.mark.slow
def test_k28_two_blast_hole_sides_are_pinned(k28_two_blasts_partitioned):
    """Sides and openings of the k=2.8 two-blast holes, as at 2315204."""
    session, fp = k28_two_blasts_partitioned
    assert _hole_sides(session.trellis(fp)) == [
        (-3, "left", [("outward", "right")]),
        (-2, "left", _OPENS_RIGHT),
        (-1, "left", _OPENS_RIGHT),
        (0, "left", _OPENS_LEFT),
    ]


@pytest.mark.slow
@pytest.mark.parametrize(
    "fixture",
    ["k10_partitioned", "k28_partitioned", "k28_two_blasts_partitioned", "p3_partitioned"],
)
def test_direct_hole_side_is_the_side_of_its_coordinates(fixture, request):
    """The crossing-sign rule agrees with measuring the hole against its bridge.

    On these cases the old nearest-vertex measurement was right, so the
    combinatorial rule (2026-10-02) must reproduce it hole for hole.
    """
    session, *fixed_points = request.getfixturevalue(fixture)
    checked = 0
    for fp in fixed_points:
        trellis = session.trellis(fp)
        for hole in trellis.holes:
            if hole.pair is None:
                continue
            bridge = bridge_for_pair(trellis, hole.pair)
            assert hole.bridge_side == _bridge_side_of(trellis, bridge, np.asarray(hole.coords))
            checked += 1
    assert checked
