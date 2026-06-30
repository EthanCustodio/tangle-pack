"""Strong-pip classification is restricted to one periodic point's own tangle.

A strong pip is a property of a single periodic point's homoclinic tangle, so
only crossings between two manifolds of that fixed point may disqualify a
candidate. A nested session also detects heteroclinic crossings (e.g. a period-1
outer tangle meeting a period-3 inner tangle); those are real but must not be
used to classify a strong pip. These tests pin that behaviour with a hand-built
registry so the heteroclinic disqualifier is present and unambiguous (the heavy
Henon fixture happens to grow no such crossing).
"""

from __future__ import annotations

from tanglepack.numerics.FixedPoint import FixedPoint
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry
from tanglepack.topology.StrongPip import is_strong_pip
from tanglepack.topology.Trellis import Trellis


def _fixed_point(period: int, lambda_u: float) -> FixedPoint:
    """A minimal period-`period` (no-inversion) fixed point with a positive
    unstable eigenvalue, enough for the strong-pip math (k_value, lambda_u,
    branch cycle)."""
    fp = FixedPoint(period, 1)
    fp.unstable_eigenvalues = [lambda_u] * period
    fp.set_k_value()
    return fp


def _trellis(registry: IntersectionRegistry, *fixed_points: FixedPoint) -> Trellis:
    return Trellis(
        fixed_points=list(fixed_points),
        registry=registry,
        branches={},
        bridges=[],
    )


def test_heteroclinic_crossing_does_not_disqualify_strong_pip():
    """A crossing whose unstable side is a *different* fixed point is ignored.

    q0 is a homoclinic crossing of the period-3 inner tangle. The only point
    inside its open box is a heteroclinic crossing (period-1 unstable × period-3
    stable). Since that point is not part of the period-3 tangle, q0 must still
    classify as a strong pip.
    """
    fp3 = _fixed_point(3, 4.0)
    fp1 = _fixed_point(1, 6.0)
    reg = IntersectionRegistry()

    stable_branch = (fp3, "stable", 0, 0)
    q0 = reg.add_synthetic(
        (0.0, 0.0), unstable_cdist=1.0, stable_cdist=1.0,
        manifold_a_key=(fp3, "unstable", 0, 0), manifold_b_key=stable_branch,
    )
    # Inside q0's open box (0, 1) x (0, 1) on the SAME stable branch, but its
    # unstable side belongs to fp1 — a heteroclinic crossing, so it is skipped.
    reg.add_synthetic(
        (0.1, 0.1), unstable_cdist=0.5, stable_cdist=0.5,
        manifold_a_key=(fp1, "unstable", 0, 0), manifold_b_key=stable_branch,
    )

    result = is_strong_pip(_trellis(reg, fp3, fp1), q0)

    assert result.is_strong_pip
    assert result.blocking_intersection_id is None


def test_homoclinic_crossing_still_disqualifies_strong_pip():
    """Control: the identical geometry, but with the disqualifier's unstable side
    on fp3 (a homoclinic crossing of the same tangle), correctly blocks q0."""
    fp3 = _fixed_point(3, 4.0)
    fp1 = _fixed_point(1, 6.0)
    reg = IntersectionRegistry()

    stable_branch = (fp3, "stable", 0, 0)
    q0 = reg.add_synthetic(
        (0.0, 0.0), unstable_cdist=1.0, stable_cdist=1.0,
        manifold_a_key=(fp3, "unstable", 0, 0), manifold_b_key=stable_branch,
    )
    blocker = reg.add_synthetic(
        (0.1, 0.1), unstable_cdist=0.5, stable_cdist=0.5,
        manifold_a_key=(fp3, "unstable", 0, 0), manifold_b_key=stable_branch,
    )

    result = is_strong_pip(_trellis(reg, fp3, fp1), q0)

    assert not result.is_strong_pip
    assert result.blocking_intersection_id == blocker


def test_missing_unstable_key_disqualifier_is_kept():
    """An iterated-bridge crossing has no unstable key (manifold_a_key is None);
    its stable side is the reliable discriminator, so it is treated as belonging
    to this tangle and may still disqualify q0."""
    fp3 = _fixed_point(3, 4.0)
    reg = IntersectionRegistry()

    stable_branch = (fp3, "stable", 0, 0)
    q0 = reg.add_synthetic(
        (0.0, 0.0), unstable_cdist=1.0, stable_cdist=1.0,
        manifold_a_key=(fp3, "unstable", 0, 0), manifold_b_key=stable_branch,
    )
    blocker = reg.add_synthetic(
        (0.1, 0.1), unstable_cdist=0.5, stable_cdist=0.5,
        manifold_a_key=None, manifold_b_key=stable_branch,
    )

    result = is_strong_pip(_trellis(reg, fp3), q0)

    assert not result.is_strong_pip
    assert result.blocking_intersection_id == blocker
