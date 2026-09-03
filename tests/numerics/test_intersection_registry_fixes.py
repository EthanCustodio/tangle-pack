"""Regression tests for the Phase 1 Intersection / IntersectionRegistry fixes.

Pins:
    1.9  ``Intersection.fixed_points`` must not raise IndexError when only
         ``manifold_b_key`` is set (the normal shape for crossings born on an
         iterated bridge).
    1.10 ``Intersection.synthetic`` must build with keyword arguments so the
         label does not land in the ``id`` slot, and must forward manifold keys.
    1.19 ``IntersectionRegistry._get_lambda_u`` must honour the stability
         argument and return a Python float.
"""

from __future__ import annotations

import numpy as np

from tanglepack.numerics.FixedPoint import FixedPoint
from tanglepack.numerics.Intersection import Intersection
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry


def _stub_fixed_point(lambda_u: float) -> FixedPoint:
    """A minimal period-1 FixedPoint carrying one unstable eigenvalue."""
    fp = FixedPoint(1, 1)
    fp.unstable_eigenvalues = [np.array([lambda_u])]
    fp.stable_eigenvalues = [np.array([1.0 / lambda_u])]
    return fp


# ── 1.9 Intersection.fixed_points ──────────────────────────────────────────


def test_fixed_points_with_only_b_key():
    fp = _stub_fixed_point(3.0)
    ix = Intersection(
        coords=(0.0, 0.0),
        unstable_cdist=1.0,
        stable_cdist=1.0,
        manifold_b_key=(fp, "stable", 0, 0),
    )
    assert ix.fixed_points == (fp,)


def test_fixed_points_with_only_a_key():
    fp = _stub_fixed_point(3.0)
    ix = Intersection(
        coords=(0.0, 0.0),
        unstable_cdist=1.0,
        stable_cdist=1.0,
        manifold_a_key=(fp, "unstable", 0, 0),
    )
    assert ix.fixed_points == (fp,)


def test_fixed_points_same_fixed_point_is_deduped():
    fp = _stub_fixed_point(3.0)
    ix = Intersection(
        coords=(0.0, 0.0),
        unstable_cdist=1.0,
        stable_cdist=1.0,
        manifold_a_key=(fp, "unstable", 0, 0),
        manifold_b_key=(fp, "stable", 0, 0),
    )
    assert ix.fixed_points == (fp,)


def test_fixed_points_distinct_fixed_points_in_ab_order():
    fp_a = _stub_fixed_point(3.0)
    fp_b = _stub_fixed_point(5.0)
    ix = Intersection(
        coords=(0.0, 0.0),
        unstable_cdist=1.0,
        stable_cdist=1.0,
        manifold_a_key=(fp_a, "unstable", 0, 0),
        manifold_b_key=(fp_b, "stable", 0, 0),
    )
    assert ix.fixed_points == (fp_a, fp_b)


def test_fixed_points_with_no_keys_is_empty():
    ix = Intersection(coords=(0.0, 0.0), unstable_cdist=1.0, stable_cdist=1.0)
    assert ix.fixed_points == ()


# ── 1.10 Intersection.synthetic ────────────────────────────────────────────


def test_synthetic_does_not_put_the_label_in_the_id_slot():
    ix = Intersection.synthetic((1.0, 2.0), 0.5, 0.25, label="anchor")

    assert ix.id is None
    assert ix.label == "anchor"
    assert ix.is_synthetic
    assert ix.coords == (1.0, 2.0)
    assert ix.unstable_cdist == 0.5
    assert ix.stable_cdist == 0.25


def test_synthetic_forwards_manifold_keys():
    fp = _stub_fixed_point(3.0)
    a_key = (fp, "unstable", 0, 0)
    b_key = (fp, "stable", 0, 1)
    ix = Intersection.synthetic(
        (0.0, 0.0),
        0.0,
        0.0,
        label="anchor",
        manifold_a_key=a_key,
        manifold_b_key=b_key,
    )

    assert ix.manifold_a_key == a_key
    assert ix.manifold_b_key == b_key
    assert ix.is_synthetic


def test_registry_add_synthetic_keeps_label_and_registry_id():
    fp = _stub_fixed_point(3.0)
    registry = IntersectionRegistry()
    iid = registry.add_synthetic(
        (0.0, 0.0),
        0.0,
        0.0,
        label="anchor",
        manifold_a_key=(fp, "unstable", 0, 0),
        manifold_b_key=(fp, "stable", 0, 0),
    )

    stored = registry[iid]
    assert stored.id == iid
    assert stored.label == "anchor"
    assert stored.is_synthetic


# ── 1.19 IntersectionRegistry._get_lambda_u ────────────────────────────────


def test_get_lambda_u_returns_a_python_float():
    fp = _stub_fixed_point(3.0)
    registry = IntersectionRegistry()
    ix = Intersection(
        coords=(0.0, 0.0),
        unstable_cdist=1.0,
        stable_cdist=1.0,
        manifold_a_key=(fp, "unstable", 0, 0),
    )

    lambda_u = registry._get_lambda_u(ix)
    assert type(lambda_u) is float
    assert lambda_u == 3.0


def test_get_lambda_u_honours_the_stability_argument():
    fp_a = _stub_fixed_point(3.0)
    fp_b = _stub_fixed_point(7.0)
    registry = IntersectionRegistry()
    ix = Intersection(
        coords=(0.0, 0.0),
        unstable_cdist=1.0,
        stable_cdist=7.0,
        manifold_a_key=(fp_a, "unstable", 0, 0),
        manifold_b_key=(fp_b, "stable", 0, 0),
    )

    assert registry._get_lambda_u(ix, "unstable") == 3.0
    assert registry._get_lambda_u(ix, "stable") == 7.0


def test_get_lambda_u_falls_back_to_the_other_key():
    fp_b = _stub_fixed_point(7.0)
    registry = IntersectionRegistry()
    ix = Intersection(
        coords=(0.0, 0.0),
        unstable_cdist=1.0,
        stable_cdist=1.0,
        manifold_b_key=(fp_b, "stable", 0, 0),
    )

    assert registry._get_lambda_u(ix, "unstable") == 7.0


def test_get_lambda_u_returns_none_without_keys():
    registry = IntersectionRegistry()
    ix = Intersection(coords=(0.0, 0.0), unstable_cdist=1.0, stable_cdist=1.0)
    assert registry._get_lambda_u(ix) is None


def test_on_interval_stable_uses_the_stable_side_eigenvalue():
    """f(p).stable_cdist = p.stable_cdist / lambda_u of the STABLE manifold's fp."""
    fp_a = _stub_fixed_point(3.0)
    fp_b = _stub_fixed_point(7.0)
    registry = IntersectionRegistry()
    ix = Intersection(
        coords=(0.0, 0.0),
        unstable_cdist=1.0,
        stable_cdist=7.0,
        manifold_a_key=(fp_a, "unstable", 0, 0),
        manifold_b_key=(fp_b, "stable", 0, 0),
    )
    registry.add(ix)

    # 7.0 / 7.0 = 1.0 with the stable-side eigenvalue; 7.0 / 3.0 = 2.33 with the
    # unstable-side one, which would fall outside the window.
    assert registry.on_interval(0.9, 1.1, stability="stable") == [ix]
    assert registry.on_interval(2.2, 2.5, stability="stable") == []
