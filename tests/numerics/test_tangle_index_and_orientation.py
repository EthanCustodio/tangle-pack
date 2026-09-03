"""Segment-index bookkeeping and scale-invariant geometry predicates in ``Tangle``.

Two Phase 1 fixes are pinned here.

1.3 -- A bridge is cut from an already-indexed unstable manifold and therefore
reuses that manifold's ``Point`` objects. The edge-dedup table used to drop
those segments outright (``_insert_segment`` returned ``None``), so
``Tangle._manifold_segs[bridge]`` came back empty and every ``for_manifold=``
filter, trim and ``_segments_of``-based lookup saw the bridge as segment-less.
The dedup must instead map the edge to the segment id that already exists and
register that id under the bridge as well, and removing one owner must leave the
segment indexed while another owner still references it.

1.4 -- ``_orientation`` and ``_find_true_intersection`` compared area-like
quantities (a cross product, a 2x2 determinant) against ABSOLUTE epsilons, so on
fine segments (lengths ~1e-7) a shallow but perfectly real crossing fell under
the threshold and was either missed or aborted the whole intersection pass with
``ValueError``. Both thresholds must scale with the product of the two lengths,
which makes the predicates scale invariant, and a genuinely near-parallel pair
must be logged and discarded rather than raised on.
"""

from __future__ import annotations

import numpy as np
import pytest

from tanglepack.numerics.Point import Point
from tanglepack.numerics.Tangle import Tangle, _Segment


# --------------------------------------------------------------------------- #
# 1.3 -- bridge segments survive the edge dedup
# --------------------------------------------------------------------------- #
def _bridge_node_pairs(bridge) -> list[frozenset[int]]:
    """Walk a bridge root -> tail and return its consecutive node-id pairs."""
    pairs: list[frozenset[int]] = []
    prev_point, curr_point = None, bridge.root
    while curr_point is not None:
        if curr_point is bridge.tail:
            break
        next_point = bridge.walk_fwd(prev_point, curr_point)
        if next_point is None:
            break
        pairs.append(frozenset((id(curr_point), id(next_point))))
        prev_point, curr_point = curr_point, next_point
    return pairs


def test_bridges_are_registered_with_their_own_segments(henon_tangle_with_bridges):
    """Every cut bridge owns exactly the segments of its own node walk."""
    workbench, _fp = henon_tangle_with_bridges
    tangle = workbench.Tangle

    bridges = workbench.bridges
    assert bridges, "fixture must produce bridges"

    for bridge in bridges:
        seg_ids = tangle._manifold_segs.get(bridge, set())
        assert seg_ids, "bridge was dropped by the edge dedup"

        indexed = {
            frozenset(
                (id(tangle._seg_lookup[sid].p0), id(tangle._seg_lookup[sid].p0_seg1))
            )
            for sid in seg_ids
        }
        expected = _bridge_node_pairs(bridge)
        assert len(seg_ids) == len(expected)
        assert indexed == set(expected)


def test_readding_a_bridge_is_idempotent(henon_tangle_with_bridges):
    """Re-registering an existing bridge changes no segment, and the parent keeps its."""
    workbench, fp = henon_tangle_with_bridges
    tangle = workbench.Tangle
    bridge = workbench.bridges[0]

    parent = workbench.manifolds[(fp, "unstable", 0, bridge.branch_index)]
    parent_segs_before = set(tangle._manifold_segs[parent])
    bridge_segs_before = set(tangle._manifold_segs[bridge])
    lookup_size_before = len(tangle._seg_lookup)

    tangle.add_manifold(bridge, index_segments=False)

    assert set(tangle._manifold_segs[bridge]) == bridge_segs_before
    assert len(tangle._seg_lookup) == lookup_size_before
    # the bridge's segments belong to the parent too, and re-registering the
    # bridge must not evict them from the parent or from the lookup
    assert set(tangle._manifold_segs[parent]) == parent_segs_before
    assert bridge_segs_before <= set(tangle._seg_lookup)


# --------------------------------------------------------------------------- #
# 1.4 -- scale-invariant orientation / true-intersection
# --------------------------------------------------------------------------- #
class _StubManifold:
    """Minimal stand-in carrying only what ``Tangle`` reads off a manifold."""

    def __init__(self, stability: str) -> None:
        self.stability = stability
        # Every indexed unstable curve carries its branch key (plan 2.3/2.6); the
        # stub uses a placeholder fixed point since nothing here reads it.
        self.manifold_key = (None, stability, 0, 0)


def _shallow_pair(length: float, angle: float):
    """Two segments of ``length`` crossing at ``angle`` radians at (length/2, 0)."""
    half = 0.5 * length
    seg_a = ((0.0, 0.0), (length, 0.0))
    seg_b = (
        (half - half * np.cos(angle), -half * np.sin(angle)),
        (half + half * np.cos(angle), half * np.sin(angle)),
    )
    return seg_a, seg_b, np.array([half, 0.0])


def _as_segment(seg, stability: str, sid: int) -> _Segment:
    (x0, y0), (x1, y1) = seg
    return _Segment(sid, _StubManifold(stability), Point(x0, y0, 0.0), Point(x1, y1, 1.0))


def test_shallow_crossing_on_fine_segments_is_detected():
    """A 1e-3 rad crossing of two 1e-7-long segments is a real crossing."""
    seg_a, seg_b, expected = _shallow_pair(1e-7, 1e-3)

    assert Tangle._do_segments_intersect(seg_a, seg_b)

    tangle = Tangle()
    point = tangle._find_true_intersection(
        _as_segment(seg_a, "unstable", 0), _as_segment(seg_b, "stable", 1)
    )
    assert point is not None
    assert np.linalg.norm(np.asarray(point) - expected) <= 1e-12 * 1e-7


def test_crossing_predicates_are_scale_invariant():
    """Scaling the same configuration by 1e6 changes no boolean."""
    fine_a, fine_b, _ = _shallow_pair(1e-7, 1e-3)
    coarse_a, coarse_b, _ = _shallow_pair(1e-7 * 1e6, 1e-3)

    assert Tangle._do_segments_intersect(
        fine_a, fine_b
    ) == Tangle._do_segments_intersect(coarse_a, coarse_b)

    for fine, coarse in (
        ((fine_a[0], fine_a[1], fine_b[0]), (coarse_a[0], coarse_a[1], coarse_b[0])),
        ((fine_a[0], fine_a[1], fine_b[1]), (coarse_a[0], coarse_a[1], coarse_b[1])),
        ((fine_b[0], fine_b[1], fine_a[0]), (coarse_b[0], coarse_b[1], coarse_a[0])),
        ((fine_b[0], fine_b[1], fine_a[1]), (coarse_b[0], coarse_b[1], coarse_a[1])),
    ):
        assert Tangle._orientation(*fine) == Tangle._orientation(*coarse)


def test_parallel_offset_pair_yields_no_crossing_and_no_exception():
    """A truly parallel offset pair is discarded, never raised on."""
    length = 1e-7
    seg_a = ((0.0, 0.0), (length, 0.0))
    seg_b = ((0.0, 0.1 * length), (length, 0.1 * length))

    assert not Tangle._do_segments_intersect(seg_a, seg_b)

    tangle = Tangle()
    seg_1 = _as_segment(seg_a, "unstable", 0)
    seg_2 = _as_segment(seg_b, "stable", 1)
    assert tangle._find_true_intersection(seg_1, seg_2) is None

    # and the caller drops the pair instead of aborting the intersection pass
    tangle._seg_lookup[0] = seg_1
    tangle._seg_lookup[1] = seg_2
    pair = frozenset((0, 1))
    tangle._intersecting_segments.add(pair)
    assert tangle._resolve_crossing_pair(pair) is None
    assert pair not in tangle._intersecting_segments


def test_small_tangle_crossing_count_is_unchanged(small_tangle):
    """The relative epsilons must not change the k=10 fixture's crossing count."""
    workbench, _fp = small_tangle
    coords = {ix.coords for _iid, ix in workbench.intersection_registry}
    assert len(coords) == 2
    assert len(workbench.intersection_registry) == 2
