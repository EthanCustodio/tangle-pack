"""
Regression (INT-7, WP8): two builds of one case number their crossings alike.

Segment ids used to come from one process-wide counter, so the segment ids
of a second build carried an offset, and ``resolve_crossings`` used to hand
its crossings over in the iteration order of a set of segment-id pairs. Two
builds of one case in ONE process therefore numbered the same crossings
differently (a first build in a fresh process always matched, which is why it
looked reproducible between processes). Registration now runs along each
unstable branch in canonical-distance order. No id is pinned: the two builds
are compared with each other.

Since N10a each segment index draws its own ids from 0 (``Tangle._next_sid``).
A pair's operand order follows its two segment ids, so the offset also moved
the last bits of a second build's crossings; now two builds in one process are
byte-equal.
"""

from __future__ import annotations

import numpy as np

from cases import build_k28, build_period3


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


def _bits(case) -> tuple[dict, list]:
    """Every registry float as ``(id -> values)`` and every written store column."""
    floats = {
        crossing_id: (
            crossing.unstable_cdist, crossing.stable_cdist, *map(float, crossing.coords),
            crossing.resolution_ratio, crossing.chord_deviation,
        )
        for crossing_id, crossing in case.registry
    }
    store = case.workbench.point_store
    n = len(store)
    columns = [store._xy[:n], store._cdist[:n], store._image[:n], store._preimage[:n],
               store._chain[:n], store._position[:n]]
    return floats, columns


def test_two_builds_in_one_process_are_byte_equal():
    first_floats, first_store = _bits(build_period3())
    second_floats, second_store = _bits(build_period3())
    assert first_floats == second_floats
    assert all(a.dtype == b.dtype and a.tobytes() == b.tobytes() for a, b in zip(first_store, second_store))
