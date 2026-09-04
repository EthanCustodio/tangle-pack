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
iterate has no ground truth to compare against, and on these fixtures the snap
never fires on the default path at all. The real fallback is pinned here by
coverage only: every bounded element gets a non-empty cover on the advanced
branch, spanning both of its endpoint images.
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
    _snap_to_partition_boundary,
)


#: Which iterate depths to sweep. The p3 fixture only registers +/-1, so the
#: tests assert that *some* n other than +1 was exercised rather than pinning a
#: depth the inference happens not to reach.
_STEPS = (1, -1, 2, -2, 3, -3)


# --------------------------------------------------------------------------- #
# fixtures
# --------------------------------------------------------------------------- #
@pytest.fixture
def p3_partitioned_c(henon_p3_session):
    """``(session, fp3, fp1)`` with every trellis classified, punched, partitioned.

    A local twin of ``tests/test_partition_elements.py``'s ``p3_partitioned``:
    the same session fan-outs, kept here so this file owns its own fixture.
    """
    session, fp3, fp1, _zone = henon_p3_session
    session.classify_strong_pips()
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()
    assert session.trellis(fp3).stable_partitions
    return session, fp3, fp1


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
def test_scaling_law_reproduces_the_iterate_table(p3_partitioned_c):
    """Where the table has the image, ``advance_key``/``per_step_beta`` find it too.

    The scaled answer is recomputed here from the fixed point alone — the branch
    from :meth:`FixedPoint.advance_key`, the distance from
    ``cdist * per_step_beta(stability) ** n`` — so this checks the law and
    :meth:`Trellis._scaled_image_cdist` together against the recorded iterate.
    The branch key must agree exactly; the distance agrees to
    :data:`SCALING_RTOL`, the polyline accuracy of a canonical distance.
    """
    session, _fp3, _fp1 = p3_partitioned_c
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
                    assert trellis._scaled_image_cdist(iid, n, stability) == (
                        scaled_key,
                        scaled_cdist,
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
def test_from_table_reports_which_branch_answered(p3_partitioned_c):
    """``from_table`` is True exactly when the iterate table holds the entry."""
    session, _fp3, _fp1 = p3_partitioned_c
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
def test_image_cdist_at_zero_steps_is_the_identity(p3_partitioned_c):
    """``n = 0`` answers with the crossing's own key and distance, from the table."""
    session, _fp3, _fp1 = p3_partitioned_c
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


@pytest.mark.slow
def test_image_cdist_is_additive_in_n_on_the_scaling_branch(p3_partitioned_c):
    """Two derived steps land where one derived double step does.

    ``advance_key`` is additive in ``n`` and ``per_step_beta ** n`` is
    multiplicative, so the scaling branch composes exactly — the property that
    lets a caller ask for any depth without walking the chain.
    """
    session, _fp3, _fp1 = p3_partitioned_c
    trellis = session.trellis()

    checked = 0
    for iid in trellis.registry.all_ids():
        crossing = trellis.intersection(iid)
        if crossing.manifold_b_key is None:
            continue
        fixed_point = crossing.manifold_b_key[0]
        beta = fixed_point.per_step_beta("stable")
        key_1, cdist_1 = trellis._scaled_image_cdist(iid, 1, "stable")
        key_2, cdist_2 = trellis._scaled_image_cdist(iid, 2, "stable")
        assert fixed_point.advance_key(key_1, 1) == key_2
        assert cdist_2 == pytest.approx(cdist_1 * beta, rel=1e-12)
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
def test_element_images_agree_whether_or_not_the_table_is_used(p3_partitioned_c):
    """The scaling fallback names the same elements the iterate table does.

    This is what licenses C.2: on every element whose two endpoint iterates ARE
    registered, forcing the scaling branch (``use_table=False``) must reproduce
    the table's answer exactly. It does because a scaled endpoint is snapped onto
    the partition boundary it lands within :data:`SNAP_RTOL` of.

    Read it as a test of the scaling LAW, not of the fallback in production: the
    elements swept here are exactly the ones the table already answers, and the
    snap that makes them agree does not fire on the default path at all (see
    ``test_the_snap_does_not_fire_on_the_default_path``).
    """
    session, fp3, fp1 = p3_partitioned_c

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
def test_every_bounded_element_now_has_an_image(p3_partitioned_c):
    """No bounded element answers None any more, registered iterate or not.

    Coverage is all that can be asserted for the elements whose iterates are
    unregistered: nothing in the registry knows where they went, so there is no
    ground truth to compare the scaled answer against.
    """
    session, fp3, fp1 = p3_partitioned_c

    fell_back = 0
    for fixed_point in (fp3, fp1):
        trellis = session.trellis(fixed_point)
        for result, interval in _bounded_elements(trellis):
            images = trellis.image_of_element(result, interval.element_id, 1)
            assert images is not None
            if (
                trellis.iterate(interval.lo_id, 1) is None
                or trellis.iterate(interval.hi_id, 1) is None
            ):
                fell_back += 1
    assert fell_back, "the outermost elements have no registered forward iterate"


@pytest.mark.slow
def test_image_of_element_accepts_partitions_from_another_trellis(p3_partitioned_c):
    """The all-fixed-points trellis can map elements it holds no partition for.

    It carries every branch and every crossing but no partitions of its own, so
    the caller (the dual graph) hands it the per-fixed-point results it was built
    from. The answer must be the one that trellis gives itself.
    """
    session, fp3, fp1 = p3_partitioned_c
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
def test_image_of_element_covers_both_endpoint_images(p3_partitioned_c):
    """The covering elements span both endpoint images, table or fallback.

    The same invariant ``test_image_of_element_lands_on_the_advanced_branch``
    pins for tabled elements, extended to the ones the fallback now answers for.
    """
    session, fp3, _fp1 = p3_partitioned_c
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
def test_scaled_element_image_stays_on_the_advanced_branch(p3_partitioned_c):
    """``use_table=False`` still lands on ``advance_key(branch, n)``, not elsewhere."""
    session, fp3, _fp1 = p3_partitioned_c
    trellis = session.trellis(fp3)

    checked = 0
    for result, interval in _bounded_elements(trellis):
        image_key = fp3.advance_key(result.branch_key, 1)
        element_ids = trellis.image_of_element(
            result, interval.element_id, 1, use_table=False
        )
        assert element_ids
        image_result = next(
            r
            for r in trellis.stable_partitions
            if r.branch_key == image_key and r.side == result.side
        )
        assert all(
            image_result.element(e).branch_key == image_key for e in element_ids
        )
        checked += 1
    assert checked


@pytest.mark.slow
def test_element_image_endpoints_scale_by_the_stable_factor(p3_partitioned_c):
    """The image arc of an element is its own arc contracted by ``per_step_beta``.

    A direct check that the stable factor — not the unstable one
    :meth:`Trellis.scale_cdist` would use for both — is what moves an element's
    endpoints one step forward.
    """
    session, fp3, _fp1 = p3_partitioned_c
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


@pytest.mark.slow
def test_the_snap_does_not_fire_on_the_default_path(p3_partitioned_c):
    """No endpoint the table cannot answer lands within SNAP_RTOL of a boundary.

    The snap exists so that FORCED scaling (``use_table=False``) reproduces the
    table; on these fixtures it is never reached by the production path, so the
    genuine fallback returns the raw scaled span. Pinning that keeps the two
    claims apart: if a future fixture does start snapping real fallbacks, this
    test fails and the reasoning above has to be revisited rather than quietly
    inherited.
    """
    session, fp3, fp1 = p3_partitioned_c

    untabled = 0
    for fixed_point in (fp3, fp1):
        trellis = session.trellis(fixed_point)
        for result, interval in _bounded_elements(trellis):
            image_key = fixed_point.advance_key(result.branch_key, 1)
            image_result = next(
                (
                    r
                    for r in trellis.stable_partitions
                    if r.branch_key == image_key and r.side == result.side
                ),
                None,
            )
            if image_result is None:
                continue
            for end_id in (interval.lo_id, interval.hi_id):
                if trellis.iterate(end_id, 1) is not None:
                    continue
                _key, cdist, from_table = trellis.image_cdist(end_id, 1, "stable")
                assert not from_table
                snapped = _snap_to_partition_boundary(
                    cdist, image_result, trellis.registry.cdist_tol
                )
                assert snapped == cdist, (
                    f"crossing {end_id} has no registered iterate yet its scaled "
                    f"image {cdist!r} snapped to {snapped!r}"
                )
                untabled += 1
    assert untabled, "the p3 partitions must have endpoints the table cannot map"


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
