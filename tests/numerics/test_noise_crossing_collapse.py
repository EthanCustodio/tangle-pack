"""Collapse of noise-duplicated crossings on an iterated bridge.

A forward-mapped polyline's coordinate noise can make it zig-zag across the
stable manifold at a single transversal crossing, registering it several times
with unstable cdists agreeing to ~1e-12 relative (observed on the k=10
horseshoe at blast depth 10: one crossing detected three times, spawning
degenerate clone bridges and phantom pseudoneighbor pairs). Noise flips add
detections in pairs while a genuine crossing adds one, so odd runs keep their
median and even runs drop entirely; genuinely separate crossings (relative
gaps >= ~1e-3) are untouched.

The collapse runs on freshly resolved crossings BEFORE they reach the registry
(plan row 2.1), so the surviving crossings are exactly what it returns; the only
tangle state it has to unwind is the segment pair the discarded detection came
from, which must stop being reported by the index.
"""

from __future__ import annotations

from tanglepack.numerics.Intersection import Intersection
from tanglepack.numerics.Tangle import Tangle


def _crossing(tangle: Tangle, seg_pair: tuple[int, int], u: float, s: float) -> Intersection:
    """Fabricate a resolved crossing and mark its pair as detected."""
    ix = Intersection(
        coords=(u, s),
        unstable_cdist=u,
        stable_cdist=s,
        seg_ids=frozenset(seg_pair),
    )
    tangle._intersecting_segments.add(frozenset(seg_pair))
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
    # the discarded detections' pairs are no longer reported by the index
    assert frozenset((1, 10)) not in tangle._intersecting_segments
    assert frozenset((3, 12)) not in tangle._intersecting_segments
    # the survivor's bookkeeping is intact
    assert frozenset((2, 11)) in tangle._intersecting_segments
    assert frozenset((4, 13)) in tangle._intersecting_segments


def test_even_noise_run_is_a_grazing_artifact_and_drops():
    tangle = Tangle()
    grazing = [
        _crossing(tangle, (1, 10), 47831087.6, 1.08676),
        _crossing(tangle, (2, 11), 47831087.6 + 1e-5, 1.08678),
    ]

    kept = tangle._collapse_noise_crossings(list(grazing))

    assert kept == []
    assert tangle._intersecting_segments == set()


def test_genuinely_separate_crossings_are_untouched():
    tangle = Tangle()
    crossings = [
        _crossing(tangle, (1, 10), 8343361.0, 6.25),
        _crossing(tangle, (2, 11), 8356210.2, 6.22),
        _crossing(tangle, (3, 12), 8361470.7, 6.20),
    ]

    kept = tangle._collapse_noise_crossings(list(crossings))

    assert sorted(kept, key=lambda ix: ix.unstable_cdist) == crossings
    assert len(tangle._intersecting_segments) == 3
