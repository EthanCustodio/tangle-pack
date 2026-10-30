"""
Canonical distances to double-double precision: the float a cdist is today, plus its residual.

A cdist is a PRODUCT, never an arc-length sum: a seed distance, scaled by the
stretch one map step at a time, and averaged or interpolated between two
neighbours. Deep in a blast it reaches 1e13 while two neighbouring nodes differ
by a few float ULPs, and there plain floats tie. The store therefore keeps each
cdist as two floats: ``hi``, computed by the SAME float expression as before
(so every plain cdist is bit-identical), and ``lo``, what that expression's
rounding dropped. ``hi + lo`` carries about 32 digits.

Each kernel below returns ``(hi, lo)`` with ``hi`` today's float. The pair is
NOT normalised (``|lo|`` may exceed half an ULP of ``hi``); :func:`normalise`
makes it the unique pair of its value, which is what an order or a comparison
reads. The error-free transforms are Knuth's ``two_sum`` and Dekker's
``two_prod`` (no fused multiply-add); they hold for finite cdists far from
overflow, which every cdist is.

Author's choice (2026-10-30, the deeper-blasts campaign's run 3): the residual
column, so plain cdists, merges and saves of the plain state need no re-basing.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

#: Dekker's splitter for float64, ``2**27 + 1``.
_SPLITTER = 134217729.0


def two_sum(a: ArrayLike, b: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """``s = fl(a + b)`` and the exact error ``e``: ``a + b == s + e``."""
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    s = a + b
    bb = s - a
    return s, (a - (s - bb)) + (b - bb)


def _split(a: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Dekker's split of ``a`` into two halves of 26 bits each."""
    c = _SPLITTER * a
    hi = c - (c - a)
    return hi, a - hi


def two_prod(a: ArrayLike, b: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """``p = fl(a * b)`` and the exact error ``e``: ``a * b == p + e``."""
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    p = a * b
    (ah, al), (bh, bl) = _split(a), _split(b)
    return p, ((ah * bh - p) + ah * bl + al * bh) + al * bl


def normalise(hi: ArrayLike, lo: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """The unique pair of the value ``hi + lo``: ``hi`` rounded, ``|lo| <= ulp(hi) / 2``."""
    return two_sum(hi, lo)


def dd_mul(hi: ArrayLike, lo: ArrayLike, s: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    ``(hi + lo) * s``: growth's ``stretch * cdist``, the graded seed's ``scale * cdist``.

    Returns:
        ``hi`` is ``fl(hi * s)``; ``lo`` its rounding error plus ``lo * s``.
    """
    p, e = two_prod(hi, s)
    return p, e + np.asarray(lo, dtype=np.float64) * s


def dd_div(hi: ArrayLike, lo: ArrayLike, s: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    ``(hi + lo) / s``: a stable bracket carried one map step forward.

    Returns:
        ``hi`` is ``fl(hi / s)``; ``lo`` the remainder over ``s``.
    """
    hi, s = np.asarray(hi, dtype=np.float64), np.asarray(s, dtype=np.float64)
    q = hi / s
    p, e = two_prod(q, s)
    return q, (((hi - p) - e) + lo) / s


def dd_mean(
    hi_0: ArrayLike, lo_0: ArrayLike, hi_1: ArrayLike, lo_1: ArrayLike
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    The mean of two cdists: refinement's ``0.5 * (c0 + c1)``.

    Returns:
        ``hi`` is ``fl(0.5 * fl(c0 + c1))``; ``lo`` half the sum's error and both residuals.
    """
    s, e = two_sum(hi_0, hi_1)
    return 0.5 * s, 0.5 * (e + lo_0 + lo_1)


def dd_lerp(
    t: ArrayLike, hi_a: ArrayLike, lo_a: ArrayLike, hi_b: ArrayLike, lo_b: ArrayLike
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    The cdist at parameter ``t`` between two ends: the Tangle's ``(1 - t) * ca + t * cb``.

    ``t`` stays the float it is (the projection of a point on a segment), so the
    plain value is unchanged; the residual is what interpolating the two precise
    ends at that ``t`` adds, ``1 - t`` taken exactly.

    Returns:
        ``hi`` is ``fl(fl((1 - t) * ca) + fl(t * cb))``; ``lo`` every rounding error of it.
    """
    t = np.asarray(t, dtype=np.float64)
    u, u_err = two_sum(1.0, -t)
    p_a, e_a = two_prod(u, hi_a)
    p_b, e_b = two_prod(t, hi_b)
    s, e = two_sum(p_a, p_b)
    return s, e + e_a + e_b + u_err * hi_a + u * lo_a + t * lo_b


def dd_sub(
    hi_a: ArrayLike, lo_a: ArrayLike, hi_b: ArrayLike, lo_b: ArrayLike
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """``a - b``, normalised: ``hi`` is the nearest float to the precise difference."""
    s, e = two_sum(hi_a, np.negative(hi_b))
    return normalise(s, e + (np.asarray(lo_a, dtype=np.float64) - lo_b))


def compare(
    hi_a: ArrayLike, lo_a: ArrayLike, hi_b: ArrayLike, lo_b: ArrayLike
) -> NDArray[np.int64]:
    """The sign of ``a - b`` (-1, 0 or +1), elementwise."""
    s, e = dd_sub(hi_a, lo_a, hi_b, lo_b)
    return np.where(s != 0.0, np.sign(s), np.sign(e)).astype(np.int64)


def rank(hi: ArrayLike, lo: ArrayLike) -> NDArray[np.int64]:
    """
    The dense rank of each value: equal values share one, the order is the precise one.

    A merge or a bisect that compares ranks compares the precise values, ties included.
    """
    hi, lo = normalise(hi, lo)
    order = np.lexsort((lo, hi))
    new = np.ones(len(order), dtype=bool)
    new[1:] = (hi[order][1:] != hi[order][:-1]) | (lo[order][1:] != lo[order][:-1])
    ranks = np.empty(len(order), dtype=np.int64)
    ranks[order] = np.cumsum(new) - 1
    return ranks


def key(hi: float, lo: float) -> tuple[float, float]:
    """One value's normalised pair as Python floats: a tuple compares as the precise value."""
    s, e = normalise(hi, lo)
    return float(s), float(e)
