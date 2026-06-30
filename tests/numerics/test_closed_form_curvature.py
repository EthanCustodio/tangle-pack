"""The closed-form curvature path must match the reference linear-algebra forms.

Two equivalences are pinned here:

* :meth:`ManifoldMachine._parabolic_fit` (closed-form divided differences) must
  reproduce the old Vandermonde-inverse coefficients.
* :meth:`ManifoldMachine._curvature_area_batch` (vectorized) must reproduce the
  scalar :meth:`ManifoldMachine._curvature_area` for every consecutive pair of a
  real grown manifold.
"""

from __future__ import annotations

import numpy as np
import pytest

from tanglepack import ManifoldView, ManifoldMachine


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


def test_curvature_area_batch_matches_scalar(grown_unstable):
    workbench, fp, manifold = grown_unstable
    machine = workbench._man_machine
    viewer = ManifoldView(manifold, machine.system)

    nan_row = np.array([np.nan, np.nan])
    p0_xy, p1_xy, left_xy, right_xy, scalar = [], [], [], [], []

    prev = manifold.root
    cur = manifold.walk_fwd(None, prev)
    while cur is not None:
        left = manifold.walk_back(cur, prev)
        right = manifold.walk_fwd(prev, cur)
        p0_xy.append(prev.get_point())
        p1_xy.append(cur.get_point())
        left_xy.append(left.get_point() if left is not None else nan_row)
        right_xy.append(right.get_point() if right is not None else nan_row)
        scalar.append(machine._curvature_area((prev, cur), viewer))
        if cur is manifold.tail:
            break
        nxt = manifold.walk_fwd(prev, cur)
        prev, cur = cur, nxt

    batch = ManifoldMachine._curvature_area_batch(
        np.array(p0_xy), np.array(p1_xy), np.array(left_xy), np.array(right_xy)
    )

    assert len(batch) == len(scalar)
    assert np.allclose(batch, np.array(scalar), rtol=1e-9, atol=1e-15)
