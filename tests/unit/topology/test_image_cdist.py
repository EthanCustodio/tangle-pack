"""Phase C — images of crossings without the iterate table.

:meth:`Trellis.image_cdist` answers "where does this crossing land after ``n``
map steps, on this side" from the registry's iterate table where the table has
the entry and from the scaling law (:meth:`FixedPoint.advance_key` for the
branch, ``per_step_beta ** n`` for the distance) where it does not. The scaling
branch is checked against the table wherever both are available: the branch keys
must match exactly and the distances to the accuracy the scaling law has.

That accuracy is NOT the registry's ``cdist_tol``. A canonical distance is an arc
length measured on a polyline and the branches of one orbit are refined
independently, so the two answers agree to about 1e-3 relative, not 1e-6 — which
is why :data:`~tanglepack.topology.Trellis.SCALING_RTOL` (the agreement bound a
test may demand) exists. The one element image the pipeline uses,
:func:`~tanglepack.topology.DualWalk.land_element`, is tested in
``test_dual_walk.py`` (WP11 deleted the unused ``Trellis.image_of_element``).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless: the session fixture touches the plotting stack

import numpy as np
import pytest

from tanglepack.numerics.Intersection import Intersection
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry
from tanglepack.topology.Trellis import SCALING_RTOL, Trellis


#: Which iterate depths to sweep. The p3 fixture only registers +/-1, so the
#: tests assert that *some* n other than +1 was exercised rather than pinning a
#: depth the inference happens not to reach.
_STEPS = (1, -1, 2, -2, 3, -3)


# --------------------------------------------------------------------------- #
# fixtures
# --------------------------------------------------------------------------- #
def _bounded_elements(trellis):
    """Yield ``(result, interval)`` for every element with two real endpoints."""
    for result in trellis.stable_partitions:
        for interval in result.intervals:
            if interval.lo_id is None or interval.hi_id is None:
                continue
            yield result, interval


# --------------------------------------------------------------------------- #
# C.1 -- image_cdist against the iterate table
# --------------------------------------------------------------------------- #
def test_scaling_law_reproduces_the_iterate_table(p3_partitioned):
    """Where the table has the image, ``advance_key``/``per_step_beta`` find it too.

    The scaled answer is recomputed here from the fixed point alone — the branch
    from :meth:`FixedPoint.advance_key`, the distance from
    ``cdist * per_step_beta(stability) ** n`` — so this checks the law against
    the recorded iterate.
    The branch key must agree exactly; the distance agrees to
    :data:`SCALING_RTOL`, the polyline accuracy of a canonical distance.
    """
    session, _fp3, _fp1 = p3_partitioned
    trellis = session.trellis()

    checked = 0
    steps_seen: set[int] = set()
    for branch in trellis.branches.values():
        for iid in branch.intersection_ids:
            for n in _STEPS:
                if trellis.iterate(iid, n) is None:
                    continue
                for stability in ("unstable", "stable"):
                    crossing = trellis.intersection(iid)
                    key = (
                        crossing.manifold_a_key
                        if stability == "unstable"
                        else crossing.manifold_b_key
                    )
                    if key is None:
                        continue
                    table_key, table_cdist, from_table = trellis.image_cdist(
                        iid, n, stability
                    )
                    assert from_table

                    fixed_point = key[0]
                    cdist = (
                        crossing.unstable_cdist
                        if stability == "unstable"
                        else crossing.stable_cdist
                    )
                    scaled_key = fixed_point.advance_key(key, n)
                    scaled_cdist = (
                        float(cdist) * fixed_point.per_step_beta(stability) ** n
                    )

                    assert scaled_key == table_key, (
                        f"crossing {iid} after {n} steps: the scaling law puts it "
                        f"on {scaled_key[1:]}, the table on {table_key[1:]}"
                    )
                    # An anchor sits at cdist 0 on both sides and stays there,
                    # so the denominator is floored at the registry tolerance.
                    error = abs(scaled_cdist - table_cdist) / max(
                        abs(table_cdist), trellis.registry.cdist_tol
                    )
                    assert error < SCALING_RTOL, (
                        f"crossing {iid}, {stability}, n={n}: scaled "
                        f"{scaled_cdist!r} vs table {table_cdist!r}"
                    )
                    checked += 1
                    steps_seen.add(n)

    assert checked, "the p3 fixture must register some iterates"
    assert steps_seen - {1}, "at least one n other than +1 must be exercised"


def test_from_table_reports_which_branch_answered(p3_partitioned):
    """``from_table`` is True exactly when the iterate table holds the entry."""
    session, _fp3, _fp1 = p3_partitioned
    trellis = session.trellis()

    tabled = derived = 0
    for iid in trellis.registry.all_ids():
        crossing = trellis.intersection(iid)
        for stability in ("unstable", "stable"):
            key = (
                crossing.manifold_a_key
                if stability == "unstable"
                else crossing.manifold_b_key
            )
            if key is None:
                continue
            for n in _STEPS:
                _key, _cdist, from_table = trellis.image_cdist(iid, n, stability)
                if trellis.iterate(iid, n) is None:
                    assert not from_table
                    derived += 1
                else:
                    assert from_table
                    tabled += 1
    assert tabled, "some crossings have registered iterates"
    assert derived, "some crossings do not, which is what the fallback is for"


def test_image_cdist_at_zero_steps_is_the_identity(p3_partitioned):
    """``n = 0`` answers with the crossing's own key and distance, from the table."""
    session, _fp3, _fp1 = p3_partitioned
    trellis = session.trellis()

    checked = 0
    for iid in trellis.registry.all_ids():
        crossing = trellis.intersection(iid)
        for stability, key, cdist in (
            ("unstable", crossing.manifold_a_key, crossing.unstable_cdist),
            ("stable", crossing.manifold_b_key, crossing.stable_cdist),
        ):
            if key is None:
                continue
            assert trellis.image_cdist(iid, 0, stability) == (key, float(cdist), True)
            checked += 1
    assert checked


# --------------------------------------------------------------------------- #
# C.1 -- error cases
# --------------------------------------------------------------------------- #
def _keyless_trellis():
    """A one-crossing trellis whose crossing carries neither branch key."""
    registry = IntersectionRegistry()
    iid = registry.add(
        Intersection.synthetic(
            coords=(0.0, 0.0),
            unstable_cdist=1.0,
            stable_cdist=1.0,
        )
    )
    return Trellis(fixed_points=[], registry=registry, branches={}, bridges=[]), iid


def test_image_cdist_rejects_a_crossing_without_the_requested_key():
    """A crossing with no branch on that side cannot be mapped along it."""
    trellis, iid = _keyless_trellis()
    for stability in ("unstable", "stable"):
        with pytest.raises(ValueError):
            trellis.image_cdist(iid, 1, stability)
        with pytest.raises(ValueError):
            trellis.image_cdist(iid, 0, stability)


def test_image_cdist_rejects_an_unknown_stability():
    """The stability parameter is the two literals, not a free string."""
    trellis, iid = _keyless_trellis()
    with pytest.raises(ValueError):
        trellis.image_cdist(iid, 1, "sideways")


# --------------------------------------------------------------------------- #
# C.1 -- element ends move by the stable factor
# --------------------------------------------------------------------------- #
def test_element_image_endpoints_scale_by_the_stable_factor(p3_partitioned):
    """The image arc of an element is its own arc contracted by ``per_step_beta``.

    A direct check that the stable factor — not the reciprocal of the unstable
    one — is what moves an element's endpoints one step forward.
    """
    session, fp3, _fp1 = p3_partitioned
    trellis = session.trellis(fp3)
    beta = fp3.per_step_beta("stable")
    assert beta < 1.0

    checked = 0
    for result, interval in _bounded_elements(trellis):
        for end_id, cdist in (
            (interval.lo_id, interval.lo_cdist),
            (interval.hi_id, interval.hi_cdist),
        ):
            _key, image, _from_table = trellis.image_cdist(end_id, 1, "stable")
            assert image == pytest.approx(cdist * beta, rel=SCALING_RTOL)
            checked += 1
    assert checked


# --------------------------------------------------------------------------- #
# C.1 on the k=10 single saddle (cheap, and a period-1 orbit exercises the
# degenerate advance_key where every image key is the source key)
# --------------------------------------------------------------------------- #
def test_period_one_images_stay_on_their_own_branch(henon_tangle_with_bridges):
    """On a period-1 saddle without inversion every image key is the source key."""
    workbench, fp = henon_tangle_with_bridges
    trellis = Trellis.from_workbench(workbench, fp)

    checked = 0
    for iid in trellis.registry.all_ids():
        crossing = trellis.intersection(iid)
        for stability, key in (
            ("unstable", crossing.manifold_a_key),
            ("stable", crossing.manifold_b_key),
        ):
            if key is None:
                continue
            for n in (1, -1, 4):
                image_key, cdist, _from_table = trellis.image_cdist(iid, n, stability)
                assert image_key == key
                assert np.isfinite(cdist)
                checked += 1
    assert checked
