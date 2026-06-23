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

from tanglepack import Point
from invariants import assert_iterate_relation


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


@pytest.mark.parametrize("stability", ["unstable", "stable"])
def test_alpha_matches_distance_ratio(initialized, stability):
    """alpha == (cdist ratio of adjacent iterates) ** (1 / k_value)."""
    workbench, fp = initialized
    manifold = workbench.manifolds[(fp, stability, 0, 0)]
    points = _real_points(manifold)

    # The growth iterate runs away from the fixed point; step toward it to find
    # the outer point's pre-image and compare cdists.
    back_attr = "prev_iterate" if stability == "unstable" else "next_iterate"
    outer = points[-1]
    inner = getattr(outer, back_attr, None)
    assert isinstance(inner, Point) and inner.cdist is not None, (
        "outer point has no real toward-fixed-point iterate"
    )

    alpha = outer.stretch_param
    assert alpha is not None and alpha > 0
    # one iterate step divides cdist by alpha
    assert np.isclose(inner.cdist, outer.cdist / alpha, rtol=1e-9)
    # and alpha is the k_value-th root of the full ratio across that step
    ratio = outer.cdist / inner.cdist
    assert np.isclose(alpha, ratio ** (1 / fp.k_value), rtol=1e-9)
