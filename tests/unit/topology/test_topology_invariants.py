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
    (heteroclinic / period > 1) are exempt. The rows are read off the crossing
    signs (D6), so the check fires exactly when the two signs are EQUAL; a
    polyline whose geometric rows disagree while the signs alternate is
    under-resolved and only flagged (``geometry_flags``), as is a propagated
    hole whose measured side disagrees with its carried one.

The synthetic tests pin the checks themselves (they must FIRE on hand-built
violations); both invariants run on every law case in
``tests/invariants/test_law_partition.py`` (I1) and ``test_law_classes.py`` (I2),
and ``tests/unit/topology/test_production_checks.py`` pins that ``punch_holes``
runs both checks itself.
"""

from __future__ import annotations

import logging

import matplotlib

matplotlib.use("Agg")  # headless: the fixtures touch the plotting stack
import numpy as np
import pytest

from helpers.fakes import bare_fixed_point
from helpers.logs import assert_logged
from tanglepack.numerics.FixedPoint import FixedPoint
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry
from tanglepack.topology.StablePartition import (
    bridge_row_diagnosis,
    bridge_row_violation,
    bridge_side_violations,
    check_bridge_rows_consistent,
    check_holes_share_bridge_side,
    geometry_flags,
    side_diagnoses,
)
from tanglepack.topology.TopologyResults import Hole
from tanglepack.topology.Trellis import Trellis
from tanglepack.topology.TrellisBranch import TrellisBranch


# --------------------------------------------------------------------------- #
# I1 — holes of one origin share their side of the bridge
# --------------------------------------------------------------------------- #
def _hole(origin, iterate, side, *, measured=None, source=None) -> Hole:
    """A bare hole carrying only the orbit identity and the sides the check reads."""
    return Hole(
        coords=(0.0, 0.0),
        near_intersection_id=0,
        bridge_side=side,
        bounding_ids=(0, 1),
        openings=[],
        iterate=iterate,
        origin=origin,
        measured_side=measured,
        side_source=source,
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
    with pytest.raises(AssertionError):
        check_holes_share_bridge_side(holes)
    assert len(bridge_side_violations(holes)) == 1


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
    with pytest.raises(AssertionError):
        check_holes_share_bridge_side(constant, orientation_preserving=False)
    assert bridge_side_violations(constant, orientation_preserving=False)

    alternating = [
        _hole((3, 4), 0, "left"),
        _hole((3, 4), -1, "right"),
        _hole((3, 4), -2, "left"),
    ]
    check_holes_share_bridge_side(alternating, orientation_preserving=False)
    with pytest.raises(AssertionError):
        check_holes_share_bridge_side(alternating, orientation_preserving=True)


@pytest.mark.parametrize(
    "periods, defaulted",
    [((1,), False), ((2,), True), ((), True)],
    ids=["odd_period_eigenvalues_decide", "even_period_only", "no_fixed_point"],
)
def test_orientation_without_a_jacobian(caplog, periods, defaulted):
    """With no jacobian an odd-period eigenvalue pair decides the orientation
    (det > 0 for a bare saddle); with only even periods, or no fixed point, the
    trellis assumes the map preserves orientation and says so at DEBUG."""
    from tanglepack.topology.Trellis import Trellis

    trellis = Trellis(
        fixed_points=[bare_fixed_point(period) for period in periods],
        registry=IntersectionRegistry(),
        branches={},
        bridges=[],
    )
    with caplog.at_level(logging.DEBUG, logger="tanglepack.topology.Trellis"):
        assert trellis.orientation_preserving is True
    assert_logged(caplog, logging.DEBUG, "tanglepack.topology.Trellis", count=int(defaulted))


def test_measured_side_disagreement_is_a_flag_not_an_error():
    """A propagated hole carries its orbit's side, so I1 holds for it by
    construction; its measured side is the geometric cross-check, and a
    disagreement is a resolution flag, never an error (D6)."""
    agreeing = [
        _hole((3, 4), 0, "left", source="signs"),
        _hole((3, 4), -1, "left", measured="left", source="carried"),
    ]
    contradicted = agreeing + [_hole((3, 4), -2, "left", measured="right", source="carried")]
    check_holes_share_bridge_side(contradicted)
    assert bridge_side_violations(contradicted) == []
    (diagnosis,) = side_diagnoses(contradicted)
    assert diagnosis.measured_disagrees and not diagnosis.orbit_disagrees
    trellis, _a_id, _b_id = _row_trellis()
    rows, sides = geometry_flags(trellis, [], contradicted)
    assert rows == [] and sides == [diagnosis]


def test_side_diagnoses_name_each_disagreeing_hole():
    """One diagnosis per disagreeing hole, saying which reading disagrees."""
    holes = [
        _hole((3, 4), 0, "left", source="signs"),
        _hole((3, 4), -1, "left", measured="right", source="carried"),
        _hole((3, 4), -2, "right", measured="right", source="carried"),
        _hole((5, 6), 0, "right", source="signs"),
    ]
    diagnoses = side_diagnoses(holes)
    assert [d.iterate for d in diagnoses] == [-1, -2]
    measured, orbit = diagnoses
    assert measured.measured_disagrees and not measured.orbit_disagrees
    assert orbit.orbit_disagrees and not orbit.measured_disagrees
    assert (orbit.expected_side, orbit.reference_iterate) == ("left", 0)
    assert len(bridge_side_violations(holes)) == 1  # the orbit check only


def test_side_diagnoses_expect_the_parity_flip():
    """Under det J < 0 the expected side of an odd iterate is the opposite one."""
    holes = [_hole((3, 4), 0, "left"), _hole((3, 4), -1, "left")]
    (diagnosis,) = side_diagnoses(holes, orientation_preserving=False)
    assert diagnosis.expected_side == "right"


def test_unmeasured_propagated_holes_are_logged(caplog):
    """A propagated hole with no measured side skips that check, at INFO."""
    holes = [_hole((3, 4), 0, "left", source="signs"), _hole((3, 4), -1, "left", source="carried")]
    with caplog.at_level(logging.INFO, logger="tanglepack.topology.StablePartition"):
        assert side_diagnoses(holes) == []
    assert_logged(caplog, logging.INFO, "tanglepack.topology.StablePartition", count=1)


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
    """A minimal bridge stand-in: endpoint registry ids, its BridgeId and a polyline.

    :func:`_row_trellis` registers its crossings in unstable-cdist order, so the
    BridgeId is the two ids sorted.
    """

    def __init__(self, first: int, second: int, points) -> None:
        self.first_intersection = first
        self.second_intersection = second
        self.id = tuple(sorted((first, second)))
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
    fp = FixedPoint(1)
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
# Its true crossing signs alternate: u+ leaves the near crossing up-right
# (cross(u+, s+) < 0) and enters the far one down-right (> 0).
_LOBE = [(2.0, 0.0), (2.5, 1.0), (3.5, 1.0), (4.0, 0.0)]
_LOBE_SIGNS = (-1, 1)
# The same arc with the far end pulled below the line: drawn on a bridge whose
# signs alternate it is an under-resolved polyline; a real arc like this would
# cross the stable manifold and carry EQUAL signs.
_CROSSING_LOBE = [(2.0, 0.0), (2.5, 1.0), (3.5, -1.0), (4.0, 0.0)]
_EQUAL_SIGNS = (-1, -1)


def _signed(trellis, a_id, b_id, signs):
    """Give the two crossings of :func:`_row_trellis` their signs."""
    trellis.intersection(a_id).crossing_sign, trellis.intersection(b_id).crossing_sign = signs


def test_bridge_rows_consistent_on_a_one_sided_lobe():
    """Alternating signs on one stable branch: one row at both ends, no flag."""
    trellis, a_id, b_id = _row_trellis()
    _signed(trellis, a_id, b_id, _LOBE_SIGNS)
    bridge = _StubBridge(a_id, b_id, _LOBE)
    check_bridge_rows_consistent(trellis, bridge)
    assert bridge_row_violation(trellis, bridge) is None
    assert geometry_flags(trellis, [bridge], []) == ([], [])


def test_bridge_rows_with_equal_signs_raise():
    """Equal signs on one stable branch put the two ends on opposite rows: the arc
    would cross its own stable branch, which the row assertion refuses (D6)."""
    trellis, a_id, b_id = _row_trellis()
    _signed(trellis, a_id, b_id, _EQUAL_SIGNS)
    bridge = _StubBridge(a_id, b_id, _CROSSING_LOBE)
    with pytest.raises(AssertionError):
        check_bridge_rows_consistent(trellis, bridge)
    assert bridge_row_violation(trellis, bridge) is not None


def test_contradicting_geometry_with_alternating_signs_is_a_flag():
    """NUM-1 in miniature: the signs alternate, so the rows agree and the check
    passes, while the polyline reads opposite rows; ``geometry_flags`` reports it."""
    trellis, a_id, b_id = _row_trellis()
    _signed(trellis, a_id, b_id, _LOBE_SIGNS)
    bridge = _StubBridge(a_id, b_id, _CROSSING_LOBE)
    check_bridge_rows_consistent(trellis, bridge)
    assert bridge_row_violation(trellis, bridge) is None
    (flag,), sides = geometry_flags(trellis, [bridge], [])
    assert sides == []
    assert flag.bridge_id == (a_id, b_id)
    assert flag.signs_alternate and flag.geometry_disagrees


def test_bridge_rows_on_different_stable_branches_are_exempt():
    """Heteroclinic / period>1 bridges span two branches and are not checked."""
    trellis, a_id, b_id = _row_trellis(split_branches=True)
    _signed(trellis, a_id, b_id, _EQUAL_SIGNS)
    bridge = _StubBridge(a_id, b_id, _CROSSING_LOBE)
    check_bridge_rows_consistent(trellis, bridge)
    assert bridge_row_violation(trellis, bridge) is None


@pytest.mark.parametrize("lobe", [_LOBE, _CROSSING_LOBE], ids=["one_sided", "crossing"])
def test_row_diagnosis_reads_every_end(lobe):
    """The diagnosis carries both crossings in BridgeId order, every row reading,
    the approach geometry, and the polyline the geometric rows were read on."""
    trellis, a_id, b_id = _row_trellis()
    _signed(trellis, a_id, b_id, _LOBE_SIGNS)
    bridge = _StubBridge(b_id, a_id, lobe[::-1])
    diagnosis = bridge_row_diagnosis(trellis, bridge)

    assert diagnosis.bridge_id == (a_id, b_id)  # ordered by unstable cdist
    assert diagnosis.branch_key == trellis.intersection(a_id).manifold_b_key
    assert diagnosis.n_nodes == len(lobe)
    assert diagnosis.chord == pytest.approx(2.0)
    assert diagnosis.signs_alternate
    first, second = diagnosis.ends
    assert (first.intersection_id, second.intersection_id) == (a_id, b_id)
    assert first.combinatorial_row == second.combinatorial_row  # alternating signs
    for end in diagnosis.ends:
        assert 0.0 < end.approach_angle <= 90.0
        assert end.displacement > 0.0
    same_row = first.geometric_row == second.geometric_row
    assert same_row == (lobe is _LOBE)
    assert diagnosis.geometry_disagrees == (lobe is _CROSSING_LOBE)
    assert bridge_row_violation(trellis, bridge) is None  # the signs decide


def test_row_diagnosis_without_a_sign_has_no_combinatorial_row():
    """A zero crossing sign leaves the combinatorial row undetermined, not an error."""
    trellis, a_id, b_id = _row_trellis()
    diagnosis = bridge_row_diagnosis(trellis, _StubBridge(a_id, b_id, _LOBE))
    assert [end.combinatorial_row for end in diagnosis.ends] == [None, None]
    assert not diagnosis.signs_alternate
