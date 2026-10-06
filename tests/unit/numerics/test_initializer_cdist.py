"""Pin the iterate/cdist law at its source: the manifold initializer.

The initial fundamental segment seeds the cdist of the first point as its
Euclidean distance from the fixed point, derives ``alpha`` (the per-step stretch)
from the ratio of consecutive iterate distances, and fills the cdist of pre-iterates
by dividing by ``alpha`` each backward step. With ``k_value = 1`` (the binary
horseshoe saddle) there are no fictitious pre-iterates, so the segment is just the
root fixed point plus two real points.
"""

from __future__ import annotations

import numpy as np
import pytest

from helpers.invariants import assert_iterate_relation
from tanglepack import FixedPoint, Point, TangleSession, TangleWorkbench
from tanglepack.examples import (
    HENON_P3,
    henon_jacobian,
    henon_map,
    henon_map_inverse,
    saddle_guesses,
)
from tanglepack.topology.Trellis import SCALING_RTOL


def _real_points(manifold):
    """Manifold's actual ``Point`` nodes (skips the root/intersection BranchPoints)."""
    return [
        n
        for n in manifold.get_point_array(return_nodes=True)
        if isinstance(n, Point) and n.cdist is not None
    ]


@pytest.mark.parametrize("stability", ["unstable", "stable"])
def test_first_point_cdist_is_euclidean_distance(initialized, stability):
    workbench, fp = initialized
    manifold = workbench.manifolds[(fp, stability, 0, 0)]
    first = _real_points(manifold)[0]

    fp_coord = np.asarray(fp.coordinates[0], dtype=float).ravel()[:2]
    euclid = float(np.linalg.norm(first.get_point() - fp_coord))
    assert np.isclose(first.cdist, euclid, rtol=1e-9), (
        f"first point cdist {first.cdist!r} != euclidean distance {euclid!r}"
    )


@pytest.mark.parametrize("stability", ["unstable", "stable"])
def test_iterate_law_holds_on_initial_segment(initialized, stability):
    """c_iterate = stretch_param * c on the freshly-seeded segment."""
    workbench, fp = initialized
    manifold = workbench.manifolds[(fp, stability, 0, 0)]
    assert_iterate_relation(manifold, rtol=1e-9)


def _p3_initialized() -> tuple[TangleWorkbench, FixedPoint]:
    """``(workbench, fp)``: the inner period-3 orbit with its fundamental segments.

    Returns:
        The workbench and the period-3 fixed point (``k_value = 3``).
    """
    session = TangleSession(
        henon_map(*HENON_P3), henon_map_inverse(*HENON_P3), henon_jacobian(*HENON_P3)
    )
    fp = session.construct_fixed_point(saddle_guesses(*HENON_P3)["period_3"])
    session.orient_eigenvectors(
        fp, {"unstable": np.array([0, -1]), "stable": np.array([-1, -1])}
    )
    session.initialize_both_manifolds(fp)
    return session.workbench, fp


@pytest.mark.parametrize("stability", ["unstable", "stable"])
@pytest.mark.parametrize("case", ["k10", "inversion", "p3"])
def test_alpha_is_the_per_step_factor(request, case, stability):
    """``alpha`` is the PER-MAP-STEP stretch, for k_value = 1, 2 and 3.

    The seed's two real points are ``k_value`` map steps apart (the fictitious
    pre-iterates fill the chain between them), so their cdist ratio is
    ``alpha ** k_value``; and ``alpha`` agrees with the eigenvalue's per-step
    factor (``per_step_beta``; its reciprocal on the stable side, which the
    initializer grows under the inverse map) to ``SCALING_RTOL``.
    """
    if case == "p3":
        workbench, fp = _p3_initialized()
    else:
        fixture = "initialized" if case == "k10" else "henon_inversion_initialized"
        workbench, fp = request.getfixturevalue(fixture)
    manifold = workbench.manifolds[(fp, stability, 0, 0)]
    outer = manifold.tail

    alpha = outer.stretch_param
    assert alpha is not None and alpha > 1.0

    # step k_value iterate links toward the fixed point to the other real point
    back_attr = "prev_iterate" if stability == "unstable" else "next_iterate"
    inner = outer
    for _ in range(fp.k_value):
        inner = getattr(inner, back_attr)
        assert isinstance(inner, Point) and inner.cdist is not None
        assert inner.stretch_param == pytest.approx(alpha)
    assert outer.cdist / inner.cdist == pytest.approx(alpha ** fp.k_value, rel=1e-9)

    eigen_factor = fp.per_step_beta(stability)
    if stability == "stable":
        eigen_factor = 1.0 / eigen_factor
    assert alpha == pytest.approx(eigen_factor, rel=SCALING_RTOL)
