"""Phase C — images of crossings and elements without the iterate table.

Two lookups are pinned here.

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
test may demand) exists, and why :meth:`Trellis.image_of_element` snaps a scaled
endpoint that lands within the tighter
:data:`~tanglepack.topology.Trellis.SNAP_RTOL` of a partition boundary onto it
before computing the cover.

What that buys, precisely: forcing the scaling branch (``use_table=False``) on an
element the table CAN answer reproduces the table's elements exactly, which pins
the scaling law. It says nothing about the genuine fallback — an unregistered
iterate has no ground truth to compare against. The real fallback is pinned here
by coverage only: every bounded element gets a non-empty cover on the advanced
branch, spanning both of its endpoint images. Whether the snap ever fires on the
default path is deliberately NOT pinned (a fixture tripwire, removed 2026-10-05).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless: the session fixture touches the plotting stack

import numpy as np
import pytest

from tanglepack.numerics.Intersection import Intersection
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry
from tanglepack.topology.Trellis import (
    SCALING_RTOL,
    SNAP_RTOL,
    Trellis,
)


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
@pytest.mark.slow
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


@pytest.mark.slow
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


@pytest.mark.slow
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
# C.2 -- the image_of_element fallback
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_element_images_agree_whether_or_not_the_table_is_used(p3_partitioned):
    """The scaling fallback names the same elements the iterate table does.

    This is what licenses C.2: on every element whose two endpoint iterates ARE
    registered, forcing the scaling branch (``use_table=False``) must reproduce
    the table's answer exactly. It does because a scaled endpoint is snapped onto
    the partition boundary it lands within :data:`SNAP_RTOL` of.

    Read it as a test of the scaling LAW, not of the fallback in production: the
    elements swept here are exactly the ones the table already answers.
    """
    session, fp3, fp1 = p3_partitioned

    checked = 0
    for fixed_point in (fp3, fp1):
        trellis = session.trellis(fixed_point)
        for result, interval in _bounded_elements(trellis):
            if (
                trellis.iterate(interval.lo_id, 1) is None
                or trellis.iterate(interval.hi_id, 1) is None
            ):
                continue
            from_table = trellis.image_of_element(result, interval.element_id, 1)
            derived = trellis.image_of_element(
                result, interval.element_id, 1, use_table=False
            )
            assert derived == from_table, (
                f"element {interval.element_id} of {result.branch_key[1:]} "
                f"{result.side}: table says {from_table}, scaling says {derived}"
            )
            checked += 1
    assert checked, "the p3 partitions must have elements with registered iterates"


@pytest.mark.slow
def test_image_of_element_accepts_partitions_from_another_trellis(p3_partitioned):
    """The all-fixed-points trellis can map elements it holds no partition for.

    It carries every branch and every crossing but no partitions of its own, so
    the caller (the dual graph) hands it the per-fixed-point results it was built
    from. The answer must be the one that trellis gives itself.
    """
    session, fp3, fp1 = p3_partitioned
    combined = session.trellis()

    checked = 0
    for fixed_point in (fp3, fp1):
        trellis = session.trellis(fixed_point)
        for result, interval in _bounded_elements(trellis):
            assert combined.image_of_element(
                result, interval.element_id, 1, partitions=trellis.stable_partitions
            ) == trellis.image_of_element(result, interval.element_id, 1)
            checked += 1
    assert checked

    # Without the partitions it has none to find, so it says so.
    result, interval = next(_bounded_elements(session.trellis(fp3)))
    assert not combined.stable_partitions
    with pytest.raises(ValueError):
        combined.image_of_element(result, interval.element_id, 1)


@pytest.mark.slow
def test_image_of_element_covers_both_endpoint_images(p3_partitioned):
    """The covering elements span both endpoint images, table or fallback.

    The same invariant ``test_image_of_element_lands_on_the_advanced_branch``
    pins for tabled elements, extended to the ones the fallback now answers for.
    """
    session, fp3, _fp1 = p3_partitioned
    trellis = session.trellis(fp3)
    tol = trellis.registry.cdist_tol

    checked = 0
    for result, interval in _bounded_elements(trellis):
        image_key = fp3.advance_key(result.branch_key, 1)
        element_ids = trellis.image_of_element(result, interval.element_id, 1)
        assert element_ids, "a bounded element's image arc is always covered"

        image_result = next(
            r
            for r in trellis.stable_partitions
            if r.branch_key == image_key and r.side == result.side
        )
        covered = [image_result.element(e) for e in element_ids]
        span_lo = min(iv.lo_cdist for iv in covered)
        span_hi = max(iv.hi_cdist for iv in covered)
        for end_id in (interval.lo_id, interval.hi_id):
            _key, cdist, _from_table = trellis.image_cdist(end_id, 1, "stable")
            slack = tol + abs(cdist) * SCALING_RTOL
            assert span_lo - slack <= cdist <= span_hi + slack
        checked += 1
    assert checked


@pytest.mark.slow
def test_element_image_endpoints_scale_by_the_stable_factor(p3_partitioned):
    """The image arc of an element is its own arc contracted by ``per_step_beta``.

    A direct check that the stable factor — not the unstable one
    :meth:`Trellis.scale_cdist` would use for both — is what moves an element's
    endpoints one step forward.
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


def test_the_snap_window_is_tighter_than_the_agreement_bound():
    """SNAP_RTOL decides identity, SCALING_RTOL only bounds a comparison.

    The snap picks an element, so it takes the smallest window that still
    swallows the scaling error; the agreement bound is what a test may demand of
    two independently computed distances and is deliberately looser. Keeping the
    two apart is the point of having two constants.
    """
    assert 0 < SNAP_RTOL < SCALING_RTOL


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
