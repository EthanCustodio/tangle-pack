"""Collapse of noise-duplicated crossings on an iterated bridge.

A forward-mapped polyline's coordinate noise can make it zig-zag across the
stable manifold at a single transversal crossing, registering it several times
with unstable cdists agreeing to ~1e-12 relative (observed on the k=10
horseshoe at blast depth 10: one crossing detected three times, spawning
degenerate clone bridges and phantom pseudoneighbor pairs). Noise flips add
detections in pairs while a genuine crossing adds one, so odd runs keep their
median and even runs drop entirely; genuinely separate crossings (relative
gaps >= ~1e-3) are untouched.
"""

from __future__ import annotations

from tanglepack.numerics.Intersection import Intersection
from tanglepack.numerics.Tangle import Tangle


def _crossing(tangle: Tangle, seg_pair: tuple[int, int], u: float, s: float) -> Intersection:
    """Register a fabricated resolved crossing in the tangle's structures."""
    ix = Intersection(
        coords=(u, s),
        unstable_cdist=u,
        stable_cdist=s,
        seg_ids=frozenset(seg_pair),
    )
    pair = frozenset(seg_pair)
    tangle._intersections.append(ix)
    for seg_id in pair:
        tangle._intersection_by_seg[seg_id].append(ix)
    tangle._intersecting_segments.add(pair)
    tangle._intersecting_coords[pair] = ix.coords
    tangle._intersecting_points[pair] = None
    return ix


def test_odd_noise_run_keeps_median_and_purges_rest():
    tangle = Tangle()
    base = 8356210.1974
    trio = [
        _crossing(tangle, (1, 10), base + 0.0, 6.2207),
        _crossing(tangle, (2, 11), base + 7e-6, 6.2205),
        _crossing(tangle, (3, 12), base + 7.4e-6, 6.2207),
    ]
    genuine = _crossing(tangle, (4, 13), base * 1.001, 6.1)

    kept = tangle._collapse_noise_crossings(trio + [genuine])

    assert kept == [trio[1], genuine]  # the median of the run, then the real one
    assert trio[0] not in tangle._intersections
    assert trio[2] not in tangle._intersections
    assert frozenset((1, 10)) not in tangle._intersecting_segments
    assert frozenset((3, 12)) not in tangle._intersecting_coords
    # the survivor's bookkeeping is intact
    assert trio[1] in tangle._intersections
    assert frozenset((2, 11)) in tangle._intersecting_segments


def test_even_noise_run_is_a_grazing_artifact_and_drops():
    tangle = Tangle()
    grazing = [
        _crossing(tangle, (1, 10), 47831087.6, 1.08676),
        _crossing(tangle, (2, 11), 47831087.6 + 1e-5, 1.08678),
    ]

    kept = tangle._collapse_noise_crossings(list(grazing))

    assert kept == []
    assert tangle._intersections == []


def test_genuinely_separate_crossings_are_untouched():
    tangle = Tangle()
    crossings = [
        _crossing(tangle, (1, 10), 8343361.0, 6.25),
        _crossing(tangle, (2, 11), 8356210.2, 6.22),
        _crossing(tangle, (3, 12), 8361470.7, 6.20),
    ]

    kept = tangle._collapse_noise_crossings(list(crossings))

    assert sorted(kept, key=lambda ix: ix.unstable_cdist) == crossings
    assert len(tangle._intersections) == 3
