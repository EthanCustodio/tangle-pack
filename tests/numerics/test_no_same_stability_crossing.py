"""Two unstable curves can never transversally cross (fundamental invariant).

A self-crossing of one unstable manifold, or a crossing between two unstable
manifolds, is physically impossible and signals a numerics defect (an early kink
sending the curve off-track). This checks the grown manifolds are geometrically
clean inside the tangle region, ignoring the huge under-resolved excursion of the
binary horseshoe.
"""

from __future__ import annotations

import numpy as np
import pytest


def _seg_cross(p1, p2, p3, p4):
    """Return True if segment p1-p2 transversally crosses p3-p4."""
    d1 = p2 - p1
    d2 = p4 - p3
    den = d1[0] * d2[1] - d1[1] * d2[0]
    if abs(den) < 1e-15:
        return False
    t = ((p3[0] - p1[0]) * d2[1] - (p3[1] - p1[1]) * d2[0]) / den
    u = ((p3[0] - p1[0]) * d1[1] - (p3[1] - p1[1]) * d1[0]) / den
    # strict interior crossing (exclude shared endpoints / grazing tangencies)
    return 1e-9 < t < 1 - 1e-9 and 1e-9 < u < 1 - 1e-9


def _region_segments(pts, lim=6.0):
    idx = np.where((np.abs(pts[:, 0]) < lim) & (np.abs(pts[:, 1]) < lim))[0]
    return [i for i in idx if i + 1 < len(pts)]


def _count_self_crossings(pts):
    segs = _region_segments(pts)
    count = 0
    for a in range(len(segs)):
        i = segs[a]
        for b in range(a + 1, len(segs)):
            j = segs[b]
            if abs(i - j) < 2:  # skip adjacent segments (shared vertex)
                continue
            if _seg_cross(pts[i], pts[i + 1], pts[j], pts[j + 1]):
                count += 1
    return count


@pytest.mark.slow
def test_period3_unstable_manifolds_do_not_self_cross(henon_p3_session):
    session, _fp3, _fp1, _zone = henon_p3_session
    manifolds = [
        np.asarray(M.get_point_array())
        for (kfp, stab, oi, bi), M in session.workbench.manifolds.items()
        if stab == "unstable"
    ]
    for pts in manifolds:
        if len(pts) > 1:
            assert _count_self_crossings(pts) == 0, (
                "unstable manifold transversally self-crosses in the tangle region"
            )
