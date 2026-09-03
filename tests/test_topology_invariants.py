"""Topological invariants that every punched trellis must satisfy.

Two facts are physically exact rather than incidental:

I1. All holes descending from one reference pseudoneighbor pair (same
    ``Hole.origin``) sit on the same side of their bridge. A propagated hole is
    the backward image of the reference hole, and the map carries the bridge's
    dynamical orientation to the image bridge's; an orientation-PRESERVING map
    (det J > 0, e.g. Hénon with b = 1) therefore keeps left/right fixed along
    the whole backward chain. Under an orientation-REVERSING map the side
    alternates with the parity of ``Hole.iterate``.
I2. A bridge whose two defining crossings lie on the SAME stable branch
    approaches both of them from the same side of that branch — the lobe lies
    entirely on one side of the stable manifold, because the bridge's endpoints
    are consecutive crossings so the arc between them never crosses it again.
    Bridges whose two crossings sit on different stable branches
    (heteroclinic / period > 1) are exempt.

The synthetic tests pin the checks themselves (they must FIRE on hand-built
violations); the Hénon tests run them on the real k = 10 tangle.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless: the fixtures touch the plotting stack
import numpy as np
import pytest

from tanglepack.numerics.FixedPoint import FixedPoint
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry
from tanglepack.topology.StablePartition import (
    bridge_row_violation,
    bridge_side_violations,
    check_bridge_rows_consistent,
    check_holes_share_bridge_side,
)
from tanglepack.topology.TopologyResults import Hole
from tanglepack.topology.Trellis import Trellis
from tanglepack.topology.TrellisBranch import TrellisBranch


# --------------------------------------------------------------------------- #
# I1 — holes of one origin share their side of the bridge
# --------------------------------------------------------------------------- #
def _hole(origin, iterate, side) -> Hole:
    """A bare hole carrying only the orbit identity the check reads."""
    return Hole(
        coords=(0.0, 0.0),
        near_intersection_id=0,
        bridge_side=side,
        bounding_ids=(0, 1),
        openings=[],
        iterate=iterate,
        origin=origin,
    )


def test_consistent_origin_sides_pass():
    """A backward chain that keeps one side is exactly the expected shape."""
    holes = [_hole((3, 4), i, "left") for i in (0, -1, -2, -3)]
    check_holes_share_bridge_side(holes)
    assert bridge_side_violations(holes) == []


def test_split_origin_sides_raise():
    """One stray side in an orbit is a violation, and it names the culprits."""
    holes = [
        _hole((3, 4), 0, "right"),
        _hole((3, 4), -1, "right"),
        _hole((3, 4), -2, "left"),
    ]
    with pytest.raises(AssertionError) as excinfo:
        check_holes_share_bridge_side(holes)
    message = str(excinfo.value)
    assert "(3, 4)" in message
    assert "-2" in message
    assert "left" in message and "right" in message


def test_distinct_origins_may_differ():
    """Two orbits landing on the same bridge legitimately differ in side."""
    holes = [_hole((3, 4), 0, "left"), _hole((5, 6), 0, "right")]
    check_holes_share_bridge_side(holes)


def test_holes_without_side_or_origin_are_skipped():
    """Degenerate geometry (no side) and unknown provenance are not violations."""
    holes = [
        _hole((3, 4), 0, "left"),
        _hole((3, 4), -1, None),
        _hole(None, -2, "right"),
    ]
    check_holes_share_bridge_side(holes)


def test_orientation_reversing_requires_alternating_sides():
    """Under det J < 0 the side flips every step, so a constant side violates."""
    constant = [_hole((3, 4), i, "left") for i in (0, -1, -2)]
    check_holes_share_bridge_side(constant, orientation_preserving=True)
    with pytest.raises(AssertionError) as excinfo:
        check_holes_share_bridge_side(constant, orientation_preserving=False)
    assert "parity" in str(excinfo.value)

    alternating = [
        _hole((3, 4), 0, "left"),
        _hole((3, 4), -1, "right"),
        _hole((3, 4), -2, "left"),
    ]
    check_holes_share_bridge_side(alternating, orientation_preserving=False)
    with pytest.raises(AssertionError):
        check_holes_share_bridge_side(alternating, orientation_preserving=True)


def test_orientation_reversing_skips_holes_without_iterate():
    """Parity is undecidable without an iterate, so such holes are skipped."""
    holes = [_hole((3, 4), 0, "left"), _hole((3, 4), None, "left")]
    check_holes_share_bridge_side(holes, orientation_preserving=False)


# --------------------------------------------------------------------------- #
# I2 — a bridge's two ends agree on their row
# --------------------------------------------------------------------------- #
class _StubNode:
    """A manifold node: a phase-space point at a canonical distance."""

    def __init__(self, point, cdist: float) -> None:
        self._point = np.asarray(point, dtype=np.float64)
        self.cdist = float(cdist)

    def get_point(self):
        return self._point


class _StubManifold:
    """A straight stable branch: nodes along y = 0 at cdist = x."""

    def __init__(self, nodes) -> None:
        self._nodes = list(nodes)

    def get_point_array(self, return_nodes: bool = False):
        if return_nodes:
            return self._nodes
        return np.array([n.get_point() for n in self._nodes])


class _StubBridge:
    """A minimal bridge stand-in: endpoint registry ids and a fixed polyline."""

    def __init__(self, first: int, second: int, points) -> None:
        self.first_intersection = first
        self.second_intersection = second
        self._points = np.asarray(points, dtype=np.float64)

    def get_point_array(self, return_nodes: bool = False):
        return self._points


def _row_trellis(*, split_branches: bool = False):
    """``(trellis, a_id, b_id)``: a straight stable branch with two crossings.

    The branch runs along ``y = 0`` with canonical distance equal to ``x``; the
    two crossings sit at ``x = 2`` and ``x = 4``. With ``split_branches`` the
    outward crossing is attributed to a second stable branch, which is what
    makes a bridge exempt from the row invariant.
    """
    fp = FixedPoint(1, 1)
    fp.unstable_eigenvalues = [4.0]
    fp.set_k_value()
    reg = IntersectionRegistry()
    stable = (fp, "stable", 0, 0)
    other = (fp, "stable", 0, 1)
    unstable = (fp, "unstable", 0, 0)
    a_id = reg.add_synthetic(
        (2.0, 0.0), unstable_cdist=1.0, stable_cdist=2.0,
        manifold_a_key=unstable, manifold_b_key=stable,
    )
    b_id = reg.add_synthetic(
        (4.0, 0.0), unstable_cdist=2.0, stable_cdist=4.0,
        manifold_a_key=unstable, manifold_b_key=other if split_branches else stable,
    )
    nodes = [_StubNode((float(x), 0.0), float(x)) for x in range(0, 8)]
    branches = {
        stable: TrellisBranch(
            key=stable, fixed_point=fp, stability="stable",
            orbit_index=0, branch_index=0, intersection_ids=[a_id],
        ),
        other: TrellisBranch(
            key=other, fixed_point=fp, stability="stable",
            orbit_index=0, branch_index=1, intersection_ids=[b_id],
        ),
    }
    if not split_branches:
        branches[stable].intersection_ids.append(b_id)
        branches.pop(other)
    trellis = Trellis(
        fixed_points=[fp], registry=reg, branches=branches, bridges=[],
        manifolds={stable: _StubManifold(nodes), other: _StubManifold(nodes)},
    )
    return trellis, a_id, b_id


# A lobe hanging above the stable line: both ends are approached from the same
# side, which is what a real bridge between consecutive crossings looks like.
_LOBE = [(2.0, 0.0), (2.5, 1.0), (3.5, 1.0), (4.0, 0.0)]
# The same arc with the far end pulled below the line: the arc now crosses the
# stable manifold between its own defining crossings, which cannot happen.
_CROSSING_LOBE = [(2.0, 0.0), (2.5, 1.0), (3.5, -1.0), (4.0, 0.0)]


def test_bridge_rows_consistent_on_a_one_sided_lobe():
    """A lobe entirely on one side of the stable branch passes the check."""
    trellis, a_id, b_id = _row_trellis()
    bridge = _StubBridge(a_id, b_id, _LOBE)
    check_bridge_rows_consistent(trellis, bridge)
    assert bridge_row_violation(trellis, bridge) is None


def test_bridge_rows_inconsistent_raise():
    """An arc that crosses its own stable branch fires the row assertion."""
    trellis, a_id, b_id = _row_trellis()
    bridge = _StubBridge(a_id, b_id, _CROSSING_LOBE)
    with pytest.raises(AssertionError) as excinfo:
        check_bridge_rows_consistent(trellis, bridge)
    message = str(excinfo.value)
    assert f"({a_id}, {b_id})" in message
    assert "left" in message and "right" in message


def test_bridge_rows_on_different_stable_branches_are_exempt():
    """Heteroclinic / period>1 bridges span two branches and are not checked."""
    trellis, a_id, b_id = _row_trellis(split_branches=True)
    bridge = _StubBridge(a_id, b_id, _CROSSING_LOBE)
    check_bridge_rows_consistent(trellis, bridge)
    assert bridge_row_violation(trellis, bridge) is None


# --------------------------------------------------------------------------- #
# Both invariants on the real k = 10 tangle
# --------------------------------------------------------------------------- #
@pytest.fixture
def henon_punched(henon_tangle_with_bridges):
    """A k = 10 trellis with pseudoneighbors computed and every hole punched."""
    workbench, fp = henon_tangle_with_bridges
    trellis = Trellis.from_workbench(workbench, fp)
    trellis.classify_strong_pips()
    trellis.compute_pseudoneighbors()
    trellis.punch_holes()
    return trellis


def test_henon_holes_share_bridge_side(henon_punched):
    """The k = 10 backward chains never switch sides (Hénon b = 1 preserves
    orientation, so no parity flip is expected)."""
    assert henon_punched.holes
    check_holes_share_bridge_side(henon_punched.holes)


def test_henon_bridge_rows_consistent(henon_punched):
    """Every k = 10 bridge lies on one side of its stable branch at both ends."""
    assert henon_punched.bridges
    for bridge in henon_punched.bridges:
        check_bridge_rows_consistent(henon_punched, bridge)


def test_henon_partition_still_built_after_wiring(henon_punched):
    """The production wiring is observational: partitioning is unaffected."""
    fp = henon_punched.fixed_points[0]
    results = henon_punched.partition_stable_manifold((fp, "stable", 0, 0))
    assert [r.side for r in results] == ["left", "right"]
    assert all(r.intervals for r in results)
