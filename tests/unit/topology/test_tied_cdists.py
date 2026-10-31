"""Topology readers on crossings whose plain cdists tie: the residuals order them.

Deep in a blast two crossings can share a plain unstable cdist and differ only
in their residuals (:mod:`tanglepack.numerics.precise`). No law case holds such
a tie, so these hand-built registries give the readers of
``Intersection.cdist_key`` their teeth: each crossing is registered AGAINST its
precise order, so a reader on plain floats (a stable sort keeping insertion
order, or a bisect on equal floats) answers wrongly.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from helpers.fakes import bare_fixed_point
from tanglepack.numerics.Intersection import Intersection
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry
from tanglepack.topology.StablePartition import _containing_bridge, _near_far
from tanglepack.topology.StrongPip import is_strong_pip
from tanglepack.topology.Trellis import Trellis
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
