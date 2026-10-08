"""
Regression (INT-7, WP8): two builds of one case number their crossings alike.

``Tangle._ids`` is a process-wide counter, so the segment ids of a second
build carry an offset, and ``resolve_crossings`` used to hand its crossings
over in the iteration order of a set of segment-id pairs. Two builds of one
case in ONE process therefore numbered the same crossings differently (a
first build in a fresh process always matched, which is why it looked
reproducible between processes). Registration now runs along each unstable
branch in canonical-distance order. No id is pinned: the two builds are
compared with each other.
"""

from __future__ import annotations

from cases import build_k28


def _identities(case) -> dict:
    """``id -> (unstable key, stable key, unstable cdist, stable cdist)``, keys without the fixed point."""
    return {
        crossing_id: (
            crossing.manifold_a_key[1:],
            crossing.manifold_b_key[1:],
            round(crossing.unstable_cdist, 9),
            round(crossing.stable_cdist, 9),
        )
        for crossing_id, crossing in case.registry
    }


def test_two_builds_in_one_process_number_their_crossings_alike():
    first = _identities(build_k28(blasts=2))
    second = _identities(build_k28(blasts=2))
    assert first == second
