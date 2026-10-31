"""Topology readers on crossings whose plain cdists tie: the residuals order them.

Deep in a blast two crossings can share a plain unstable cdist and differ only
in their residuals (:mod:`tanglepack.numerics.precise`). No law case holds such
a tie, so these hand-built registries give the readers of
``Intersection.cdist_key`` their teeth (the Trellis branch sort, the container
bisect, near and far, the strong-pip box, the pseudoneighbor interval test and
a reference hole's side): each tie is built AGAINST its precise order, so a
reader on plain floats (a stable sort keeping insertion order, a bisect or a
comparison on equal floats) answers wrongly.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from helpers.fakes import bare_fixed_point
from tanglepack.numerics.Intersection import Intersection
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry
from tanglepack.topology.Pseudoneighbor import _interval_is_clear
from tanglepack.topology.StablePartition import _containing_bridge, _near_far, punch_holes
from tanglepack.topology.StrongPip import is_strong_pip
from tanglepack.topology.Trellis import Trellis
from tanglepack.topology.TopologyResults import PseudoneighborPair
from tanglepack.topology.TrellisBranch import TrellisBranch

DEEP = 1e13
STEP = 1e-4


def _crossing(fp, residual: float, stable_cdist: float, *, unstable_cdist: float = DEEP) -> Intersection:
    """A synthetic crossing of ``fp``'s home branches at plain unstable cdist ``unstable_cdist``."""
    return Intersection(
        coords=(stable_cdist, residual), unstable_cdist=unstable_cdist, stable_cdist=stable_cdist,
        unstable_cdist_lo=residual, manifold_a_key=(fp, "unstable", 0, 0), manifold_b_key=(fp, "stable", 0, 0),
    )


def _tied(fp, residuals) -> tuple[IntersectionRegistry, list[int]]:
    """Crossings at plain ``DEEP`` with the given residuals (in units of ``STEP``),
    registered from the last to the first; their ids in precise order."""
    registry = IntersectionRegistry()
    ids = [registry.add(_crossing(fp, r * STEP, float(i + 1))) for i, r in reversed(list(enumerate(residuals)))]
    return registry, ids[::-1]


class _StubBridge:
    """A bridge stand-in: its two endpoint ids on one unstable branch."""

    def __init__(self, first: int, second: int, manifold_key) -> None:
        self.first_intersection, self.second_intersection = first, second
        self.id = (first, second)
        self.manifold_key = manifold_key
        self.partial = False

    def get_point_array(self):
        return np.zeros((2, 2))


def test_a_branch_orders_tied_crossings_by_their_residuals():
    """``Trellis.from_workbench`` sorts each branch by the cdist key: the plain
    tie keeps no insertion order."""
    fp = bare_fixed_point(label="A")
    registry, (left, right) = _tied(fp, [1, 2])
    assert left > right, "registered against the precise order"
    keys = [(fp, "unstable", 0, 0), (fp, "stable", 0, 0)]
    workbench = SimpleNamespace(
        fixed_points=[fp], intersection_registry=registry, resolution_ratio_threshold=None,
        manifolds=dict.fromkeys(keys), bridges=[], dynamical_system=None, generation=0,
    )
    trellis = Trellis.from_workbench(workbench)
    assert trellis.branch(keys[0]).intersection_ids == [left, right]


def test_the_container_of_a_span_is_bisected_on_the_residuals():
    """Three crossings tie as plain floats: the span between the first two
    residuals lies in the first bridge, one across the second in none."""
    fp = bare_fixed_point(label="A")
    registry, ids = _tied(fp, [1, 2, 3])
    key = (fp, "unstable", 0, 0)
    bridges = [_StubBridge(a, b, key) for a, b in zip(ids, ids[1:])]
    branch = TrellisBranch(key=key, fixed_point=fp, stability="unstable", orbit_index=0, branch_index=0,
                           intersection_ids=ids)
    trellis = Trellis(fixed_points=[fp], registry=registry, branches={key: branch}, bridges=bridges)
    orders: dict = {}
    assert _containing_bridge(trellis, ((DEEP, 1.2 * STEP), (DEEP, 1.8 * STEP)), key, orders) is bridges[0]
    assert _containing_bridge(trellis, ((DEEP, 1.5 * STEP), (DEEP, 2.5 * STEP)), key, orders) is None


def test_near_and_far_of_a_tied_pair_follow_the_residuals():
    """Two crossings on one stable branch, tied as plain stable cdists: the
    nearer to the anchor is the one with the smaller residual."""
    fp = bare_fixed_point(label="A")
    registry = IntersectionRegistry()
    far = registry.add(Intersection(
        coords=(0.0, 1.0), unstable_cdist=2.0, stable_cdist=DEEP, stable_cdist_lo=2 * STEP,
        manifold_a_key=(fp, "unstable", 0, 0), manifold_b_key=(fp, "stable", 0, 0),
    ))
    near = registry.add(Intersection(
        coords=(0.0, 2.0), unstable_cdist=1.0, stable_cdist=DEEP, stable_cdist_lo=STEP,
        manifold_a_key=(fp, "unstable", 0, 0), manifold_b_key=(fp, "stable", 0, 0),
    ))
    trellis = Trellis(fixed_points=[fp], registry=registry, branches={}, bridges=[])
    assert _near_far(trellis, far, near) == (near, far)


def test_a_crossing_tied_on_the_unstable_side_falls_in_the_strong_pip_box():
    """``inside`` ties q0's plain unstable cdist with a smaller residual and has
    the smaller stable cdist: it lies in q0's open box and blocks it."""
    fp = bare_fixed_point(label="A")
    registry = IntersectionRegistry()
    q0 = registry.add(_crossing(fp, 1e-17, 1.0, unstable_cdist=1.0))
    inside = registry.add(_crossing(fp, -1e-17, 0.5, unstable_cdist=1.0))
    result = is_strong_pip(Trellis(fixed_points=[fp], registry=registry, branches={}, bridges=[]), q0)
    assert not result.is_strong_pip and result.blocking_intersection_id == inside


def _stable_crossing(fp, unstable_cdist: float, stable_cdist: float, residual: float) -> Intersection:
    """A crossing of ``fp``'s home branches with stable cdist ``stable_cdist`` (residual ``residual``)."""
    return Intersection(
        coords=(unstable_cdist, residual), unstable_cdist=unstable_cdist, stable_cdist=stable_cdist,
        stable_cdist_lo=residual, manifold_a_key=(fp, "unstable", 0, 0), manifold_b_key=(fp, "stable", 0, 0),
    )


def test_a_stable_interval_between_tied_ends_is_not_degenerate():
    """Two crossings tied as plain stable cdists with residuals that order them
    bound a real (if tiny) interval: with nothing in it, it is clear."""
    fp = bare_fixed_point(label="A")
    registry = IntersectionRegistry()
    x0 = registry.add(_stable_crossing(fp, 1.0, DEEP, -STEP))
    x1 = registry.add(_stable_crossing(fp, 2.0, DEEP, STEP))
    trellis = Trellis(fixed_points=[fp], registry=registry, branches={}, bridges=[])
    assert _interval_is_clear(trellis, x0, x1, "stable", [], 1e-2)


@pytest.mark.parametrize(("residual", "clear"), [(0.0, False), (5 * STEP, True)], ids=["inside", "outside"])
def test_a_table_landing_tied_with_the_ends_punctures_by_its_residual(residual, clear):
    """``r`` one step out maps (by the table) to ``z``, tied with the pair's ends
    as a plain stable cdist: ``z`` punctures the interval iff its residual lies
    between theirs."""
    fp = bare_fixed_point(label="A")  # beta 0.25: r's scaled landing is DEEP exactly
    registry = IntersectionRegistry()
    x0 = registry.add(_stable_crossing(fp, 1.0, DEEP, -STEP))
    x1 = registry.add(_stable_crossing(fp, 2.0, DEEP, STEP))
    z = registry.add(_stable_crossing(fp, 3.0, DEEP, residual))
    r = registry.add(_stable_crossing(fp, 0.75, 4 * DEEP, 0.0))
    registry.register_iterate(r, 1, z)
    trellis = Trellis(fixed_points=[fp], registry=registry, branches={}, bridges=[])
    assert _interval_is_clear(trellis, x0, x1, "stable", [r], 1e-2) is clear


def test_a_reference_hole_reads_its_side_off_the_residuals():
    """A reference pair tied as plain stable cdists: the second crossing (in
    unstable order) lies outward by its residual, so the segment from the first
    heads outward and the hole sits on the left of the bridge (sign +1)."""
    fp = bare_fixed_point(label="A")
    registry = IntersectionRegistry()
    first, second = (
        registry.add(Intersection(
            coords=(u, 0.0), unstable_cdist=u, stable_cdist=DEEP, stable_cdist_lo=residual,
            crossing_sign=sign, manifold_a_key=(fp, "unstable", 0, 0), manifold_b_key=(fp, "stable", 0, 0),
        ))
        for u, residual, sign in ((1.0, -STEP, 1), (2.0, STEP, -1))
    )
    branches = {
        (fp, stability, 0, 0): TrellisBranch(
            key=(fp, stability, 0, 0), fixed_point=fp, stability=stability, orbit_index=0, branch_index=0,
            intersection_ids=[first, second],
        )
        for stability in ("unstable", "stable")
    }
    bridge = _StubBridge(first, second, (fp, "unstable", 0, 0))
    trellis = Trellis(fixed_points=[fp], registry=registry, branches=branches, bridges=[bridge])
    pair = PseudoneighborPair(first, second, branch_key=(fp, "stable", 0, 0), iterate=0, origin=(first, second))
    (hole,) = punch_holes(trellis, [pair])
    assert hole.near_intersection_id == first and hole.bridge_side == "left"
