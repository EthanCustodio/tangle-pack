"""Regression: bridge boundary (straddle) points get the wrong cdist.

``Tangle._boundary_point`` places the boundary point at a 10% offset from the
crossing toward a segment endpoint (``_linear_interpolation(intersection, seg, 0.1)``)
but assigns its cdist as the value interpolated *at the crossing*, not at the
boundary point's actual offset location. The cdist a ``Point`` carries must reflect
where the point physically sits, so when the bridge is later iterated (its cdist
scaled by ``stretch_param``) this offset error is baked in and compounds.

Expected fix: assign ``_cdist_between(p0, p1, stability, boundary.get_point())``
(the cdist at the boundary's real location). Until then this xfails.
"""

from __future__ import annotations

import numpy as np
import pytest

from tanglepack import Point, Tangle


@pytest.mark.regression
@pytest.mark.parametrize("side", ["root", "tail"])
def test_boundary_cdist_matches_its_location(side):
    tangle = Tangle()

    # A simple horizontal segment with a linear cdist profile.
    p0 = Point(0.0, 0.0, cdist=1.0)
    p1 = Point(1.0, 0.0, cdist=2.0)
    intersection = (0.5, 0.0)  # cdist 1.5 here

    boundary, _t = tangle._boundary_point(p0, p1, "unstable", intersection, side)

    expected = Tangle._cdist_between(
        p0, p1, "unstable", np.asarray(boundary.get_point())
    )
    assert np.isclose(boundary.cdist, expected, rtol=1e-12), (
        f"boundary cdist {boundary.cdist} should equal cdist at its location "
        f"{expected}, not the cdist at the crossing"
    )
