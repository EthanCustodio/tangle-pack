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
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless: the session fixture touches the plotting stack
import pytest

from tanglepack.topology import (
    check_bridge_rows_consistent,
    check_holes_share_bridge_side,
)
from tanglepack.topology.Pseudoneighbor import forward_unstable_branch_cycle
from tanglepack.topology.StablePartition import bridge_for_pair


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
