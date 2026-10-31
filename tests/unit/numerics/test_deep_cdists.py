"""Deep canonical distances: neighbours closer than a float's ULP keep their order.

A hand-built unstable curve at cdist 1e13 whose nodes are 1e-4 apart
(``dc / c = 1e-17``, below float eps): every plain cdist ties, the residuals
do not (:mod:`tanglepack.numerics.precise`). Two crossings on it tie as plain
floats too; the registry, the brackets and the bridge cut order them by their
residuals (``Intersection.cdist_key``). The topology's readers of the key are
tested on tied crossings in ``tests/unit/topology/test_tied_cdists.py``.
"""

from __future__ import annotations

import numpy as np
import pytest

from helpers.fakes import bare_fixed_point
from helpers.invariants import assert_plain_ties_ordered
from helpers.polylines import TEST_STORE, polyline_manifold
from tanglepack import BaseManifold, IntersectionRegistry
from tanglepack.numerics.Tangle import Tangle

DEEP = 1e13
STEP = 1e-4


def _deep_tangle():
    """An unstable line at plain cdist ``DEEP`` (residual ``i * STEP``) crossed at
    x = 1.5 and x = 3.5 by one stable curve; its crossings, left to right."""
    fp = bare_fixed_point(label="A")
    points = np.column_stack([np.arange(6.0), np.zeros(6)])
    rows = TEST_STORE.append(points, np.full(6, DEEP), STEP * np.arange(6))
    TEST_STORE.new_chain(rows, 1.0)
    unstable = BaseManifold(
        int(rows[0]), "unstable", 1.0, fixed_point=fp, tail=int(rows[-1]), branch_index=0,
        manifold_key=(fp, "unstable", 0, 0), store=TEST_STORE,
    )
    stable = polyline_manifold([(1.5, -1.0), (1.5, 1.0), (3.5, 1.0), (3.5, -1.0)], (fp, "stable", 0, 0))
    tangle = Tangle()
    tangle.add_manifolds([unstable, stable])
    crossings = sorted(tangle.resolve_crossings(), key=lambda ix: ix.coords[0])
    assert [ix.unstable_cdist for ix in crossings] == [DEEP, DEEP], "the plain cdists tie"
    return tangle, unstable, crossings


def test_tied_crossings_keep_their_precise_order_in_the_registry():
    """Registered left to right, the right crossing would go first on a plain tie."""
    _tangle, _unstable, (left, right) = _deep_tangle()
    registry = IntersectionRegistry()
    left_id, right_id = registry.add(left), registry.add(right)
    assert left_id != right_id
    assert left.cdist_key("unstable") < right.cdist_key("unstable")
    assert [i for i in registry.by_unstable_cdist if i in (left_id, right_id)] == [left_id, right_id]


def test_brackets_carry_the_residuals():
    """Each crossing's unstable bracket is its two nodes, told apart by their residuals."""
    _tangle, _unstable, (left, right) = _deep_tangle()
    for crossing, below in ((left, 1), (right, 3)):
        hi, lo = np.array(crossing.precise_brackets)
        assert hi[:2].tolist() == [DEEP, DEEP]
        np.testing.assert_allclose(lo[:2], [below * STEP, (below + 1) * STEP], rtol=1e-12)
        _, lo_u = crossing.cdist_key("unstable")
        assert lo[0] < lo_u < lo[1]


def test_a_cut_between_tied_crossings_runs_in_precise_order():
    """The bridge between the two crossings runs from the node below the left one
    to the node above the right one, its id in precise unstable order."""
    tangle, unstable, (left, right) = _deep_tangle()
    registry = IntersectionRegistry()
    for crossing in (right, left):  # registered right first: ids against the curve order
        registry.add(crossing)
    (bridge,) = tangle.create_bridges([right, left], for_manifold=unstable)
    rows = unstable.ids
    assert bridge.id == (left.id, right.id)
    assert (bridge.root.id, bridge.tail.id) == (rows[1], rows[4])


def test_a_synthetic_crossing_has_no_brackets_and_a_zero_residual():
    """An anchor carries no segment: neither bracket reading has one, and its key is (0, 0)."""
    from tanglepack import Intersection

    anchor = Intersection.synthetic(coords=(0.0, 0.0), unstable_cdist=0.0, stable_cdist=0.0)
    assert anchor.brackets is None and anchor.precise_brackets is None
    assert anchor.cdist_key("unstable") == anchor.cdist_key("stable") == (0.0, 0.0)


def test_the_manifold_law_reads_the_residuals_at_a_plain_tie():
    """The deep curve passes the law's tie clause; the same curve with its
    residuals equal (or reversed) at a tie fails it."""
    _tangle, unstable, _crossings = _deep_tangle()
    assert_plain_ties_ordered(unstable)
    rows = unstable.ids
    TEST_STORE._cdist_lo[rows[3]] = TEST_STORE._cdist_lo[rows[2]]
    with pytest.raises(AssertionError):
        assert_plain_ties_ordered(unstable)
