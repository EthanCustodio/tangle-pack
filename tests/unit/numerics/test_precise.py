"""The double-double kernels of ``numerics/precise.py``, against exact rationals.

Each kernel's ``hi`` is the float expression the library used before the
residual column, bit for bit; ``hi + lo`` is the exact value to about 32
digits (relative 1e-30 here, far below a float's 1e-16).
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest

from tanglepack.numerics import precise

RNG = np.random.default_rng(20261030)
#: Values like deep cdists: 1e-3 .. 1e13, and stretches of a few.
VALUES = 10.0 ** RNG.uniform(-3, 13, 200)
STRETCHES = RNG.uniform(1.1, 6.0, 200)
#: Residuals as a pair carries them: below an ULP of their value.
LOWS = VALUES * RNG.uniform(-1e-16, 1e-16, 200)


def _exact(hi: float, lo: float) -> Fraction:
    return Fraction(float(hi)) + Fraction(float(lo))


def _close(pair: tuple, exact: Fraction, rtol: float = 1e-30) -> bool:
    return abs(_exact(*pair) - exact) <= rtol * abs(exact)


def test_two_sum_and_two_prod_are_exact() -> None:
    a, b = VALUES, STRETCHES * VALUES[::-1]
    for (s, e), op in ((precise.two_sum(a, b), Fraction.__add__), (precise.two_prod(a, b), Fraction.__mul__)):
        assert all(_exact(si, ei) == op(Fraction(ai), Fraction(bi))
                   for si, ei, ai, bi in zip(s, e, a, b))
    assert np.array_equal(precise.two_sum(a, b)[0], a + b)
    assert np.array_equal(precise.two_prod(a, b)[0], a * b)


def test_dd_mul_and_dd_div() -> None:
    hi, lo = precise.dd_mul(VALUES, LOWS, STRETCHES)
    assert np.array_equal(hi, STRETCHES * VALUES)
    assert all(_close(pair, _exact(v, l) * Fraction(s))
               for pair, v, l, s in zip(zip(hi, lo), VALUES, LOWS, STRETCHES))
    hi, lo = precise.dd_div(VALUES, LOWS, STRETCHES)
    assert np.array_equal(hi, VALUES / STRETCHES)
    assert all(_close(pair, _exact(v, l) / Fraction(s))
               for pair, v, l, s in zip(zip(hi, lo), VALUES, LOWS, STRETCHES))


def test_dd_mean_and_dd_lerp() -> None:
    a, b = VALUES, VALUES * (1 + RNG.uniform(0, 1e-12, 200))
    la, lb = LOWS, b * RNG.uniform(-1e-16, 1e-16, 200)
    hi, lo = precise.dd_mean(a, la, b, lb)
    assert np.array_equal(hi, 0.5 * (a + b))
    assert all(_close(pair, (_exact(x, lx) + _exact(y, ly)) / 2)
               for pair, x, lx, y, ly in zip(zip(hi, lo), a, la, b, lb))
    t = RNG.uniform(0, 1, 200)
    hi, lo = precise.dd_lerp(t, a, la, b, lb)
    assert np.array_equal(hi, (1 - t) * a + t * b)
    assert all(_close(pair, (1 - Fraction(ti)) * _exact(x, lx) + Fraction(ti) * _exact(y, ly))
               for pair, ti, x, lx, y, ly in zip(zip(hi, lo), t, a, la, b, lb))


def test_dd_sub_compare_and_rank_order_the_precise_values() -> None:
    """Values that tie as floats are told apart by their residuals."""
    hi = np.array([5.0, 5.0, 5.0, 3.0, 5.0])
    lo = np.array([1e-16, -1e-16, 0.0, 0.0, 1e-16])
    assert precise.rank(hi, lo).tolist() == [3, 1, 2, 0, 3]
    assert precise.compare(hi[:-1], lo[:-1], hi[1:], lo[1:]).tolist() == [1, -1, 1, -1]
    s, e = precise.dd_sub(5.0, 1e-16, 5.0, -1e-16)
    assert (float(s), float(e)) == (2e-16, 0.0)
    assert precise.key(5.0, 1e-16) < precise.key(5.0, 2e-16) < precise.key(np.nextafter(5.0, 6.0), 0.0)


@pytest.mark.parametrize("ulps", [1, 3])
def test_normalise_gives_the_unique_pair(ulps: int) -> None:
    """A residual of a few ULPs moves into ``hi``; the value is unchanged."""
    hi, lo = 7.0, ulps * np.spacing(7.0)
    s, e = precise.normalise(hi, lo)
    assert _exact(s, e) == _exact(hi, lo) and abs(e) <= np.spacing(s) / 2
