"""The closed-form curvature path must match the reference linear-algebra form.

:meth:`ManifoldMachine._parabolic_fit` (closed-form divided differences) must
reproduce the old Vandermonde-inverse coefficients.
"""

from __future__ import annotations

import numpy as np
import pytest

from tanglepack import ManifoldMachine


def _vandermonde_fit(points):
    x = points[:, 0]
    y = points[:, 1]
    A = np.array([[x[i] ** 2, x[i], 1] for i in range(3)], dtype=float)
    return np.linalg.solve(A, y)


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_parabolic_fit_matches_vandermonde(seed):
    rng = np.random.default_rng(seed)
    # Distinct, well-separated x-values keep the Vandermonde system conditioned.
    xs = np.sort(rng.uniform(-5, 5, size=3))
    xs = xs + np.array([-0.5, 0.0, 0.5])  # guarantee separation
    ys = rng.uniform(-5, 5, size=3)
    pts = np.column_stack([xs, ys])

    a, b, c = ManifoldMachine._parabolic_fit(pts)
    a_ref, b_ref, c_ref = _vandermonde_fit(pts)

    assert np.allclose([a, b, c], [a_ref, b_ref, c_ref], rtol=1e-9, atol=1e-12)


