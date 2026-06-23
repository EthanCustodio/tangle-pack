"""Intersection cdist interpolation and the fundamental same-stability invariant.

Every legitimate crossing is exactly one unstable + one stable segment (CLAUDE.md).
Each crossing's cdist along a manifold is interpolated between the two endpoint
cdists of the segment that produced it, so it must lie between them.
"""

from __future__ import annotations

import numpy as np


def test_intersections_have_both_cdists(small_tangle):
    workbench, fp = small_tangle
    tangle = workbench.Tangle
    assert len(tangle._intersections) > 0, "fixture produced no intersections"

    for ix in tangle._intersections:
        # both cdists must be populated; the fixed point, if detected as a
        # crossing, legitimately sits at cdist 0, so only require non-negative.
        assert ix.unstable_cdist is not None and ix.unstable_cdist >= 0
        assert ix.stable_cdist is not None and ix.stable_cdist >= 0


def test_each_crossing_is_one_unstable_one_stable(small_tangle):
    """The fundamental invariant: no same-stability pair survives as a crossing."""
    workbench, _ = small_tangle
    tangle = workbench.Tangle

    for pair in tangle._intersecting_segments:
        stabilities = sorted(
            tangle._seg_lookup[seg_id].manifold.stability for seg_id in pair
        )
        assert stabilities == ["stable", "unstable"], (
            f"same-stability crossing leaked through: {stabilities}"
        )


def test_crossing_cdist_lies_between_segment_endpoints(small_tangle):
    """Interpolated cdist is bracketed by the producing segment's endpoint cdists."""
    workbench, _ = small_tangle
    tangle = workbench.Tangle

    for ix in tangle._intersections:
        if ix.seg_ids is None:
            continue
        for seg_id in ix.seg_ids:
            seg = tangle._seg_lookup[seg_id]
            stability = seg.manifold.stability
            c0 = seg.p0.get_cdist(stability)
            c1 = seg.p0_seg1.get_cdist(stability)
            if c0 is None or c1 is None:
                continue
            lo, hi = sorted((float(c0), float(c1)))
            value = ix.unstable_cdist if stability == "unstable" else ix.stable_cdist
            # small tolerance for the offset/clamp arithmetic
            assert lo - 1e-9 <= value <= hi + 1e-9, (
                f"{stability} crossing cdist {value} not in [{lo}, {hi}]"
            )
