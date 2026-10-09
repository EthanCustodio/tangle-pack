"""PK-4 / C19: the chords of a lobe image on a FOREIGN stable branch pair in stable order.

A hole's image chain can leave its base branch for an excursion across another
stable branch (another orbit branch or fixed point). Before WP10b those
crossings were left out of the chord pairing, so a crossing such as 44 on
``p3_6_blasts`` abutted no empty stretch and two cuts were skipped. A lobe image
holds no anchor, so walking the foreign branch from its anchor enters and leaves
it in turn (Jordan alternation): the stretches inside are the crossings paired
in stable order. The geometric check below is that statement, measured.
"""

from __future__ import annotations

import numpy as np
import pytest

import cases
from tanglepack.numerics.geometry import arc_polyline, point_in_polygon
from tanglepack.topology.BridgeClass import _image_chain
from tanglepack.topology.PartitionFamily import IteratedHomotopyPartition

CASES = ["p3_6_blasts", "p3_15", "p3_16", "heteroclinic"]


@pytest.fixture(scope="module", params=CASES)
def case(request: pytest.FixtureRequest) -> cases.Case:
    return cases.BUILDERS[request.param]()


def test_every_image_crossing_abuts_an_empty_stretch(case: cases.Case) -> None:
    cuts = case.session.iterated_partition().cuts
    assert cuts
    assert not [cut for cut in cuts if cut.reason == "no empty stretch"]


def _lobe_image(trellis, workbench, chain: list, base_key) -> np.ndarray:
    """The polygon of ``f(H)``: the image chain's arcs, closed along the base."""
    def arc(a: int, b: int) -> np.ndarray:
        bridge = trellis.bridge_between(a, b)
        assert bridge is not None, (a, b)
        ends = [np.ravel(trellis.intersection(c).coords).astype(float) for c in (a, b)]
        lo, hi = sorted(trellis.intersection(c).unstable_cdist for c in (a, b))
        return np.vstack([ends[0], arc_polyline(bridge, lo, hi), ends[1]])

    first, last = chain[0][0], chain[-1][1]
    lo, hi = sorted(trellis.intersection(c).stable_cdist for c in (first, last))
    base = arc_polyline(workbench.manifolds[base_key], lo, hi, stability="stable")
    if trellis.intersection(last).stable_cdist > trellis.intersection(first).stable_cdist:
        base = base[::-1]
    return np.vstack([*(arc(a, b) for a, b in chain), base])


def test_foreign_chords_are_inside_the_lobe_image(case: cases.Case) -> None:
    """Along every foreign branch the image chain crosses, from its anchor, the
    stretches between consecutive crossings lie outside, inside, outside, ... the
    lobe image, and every inside stretch is an empty stretch of the iterated cut."""
    session, workbench = case.session, case.session.workbench
    minimal = session.minimal_trellis()
    trellis = minimal.trellis
    empty = IteratedHomotopyPartition._empty_stretches(minimal, session.homotopy_partition())
    checked = 0
    for hole_bridge in minimal.hole_bridge_ids:
        chain = _image_chain(trellis, hole_bridge)
        if not chain:
            continue
        crossings = [chain[0][0]] + [pair[1] for pair in chain]
        (base_key,) = {
            trellis.intersection(c).manifold_b_key
            for c in (crossings[0], crossings[-1])
            if not trellis.intersection(c).is_anchor
        }
        foreign = {trellis.intersection(c).manifold_b_key for c in crossings[1:-1]} - {base_key}
        if not foreign:
            continue
        ring = _lobe_image(trellis, workbench, chain, base_key)
        for key in foreign:
            on = sorted(
                (c for c in crossings[1:-1] if trellis.intersection(c).manifold_b_key == key),
                key=lambda c: trellis.intersection(c).stable_cdist,
            )
            manifold = workbench.manifolds[key]
            bounds = [0.0, *(trellis.intersection(c).stable_cdist for c in on),
                      float(manifold.tail.get_cdist("stable"))]
            for i, (lo, hi) in enumerate(zip(bounds, bounds[1:])):
                stretch = arc_polyline(manifold, lo, hi, stability="stable")
                if not len(stretch):
                    continue
                inside = point_in_polygon(stretch[len(stretch) // 2], ring, tol=0.0)
                assert inside == bool(i % 2), (hole_bridge, key[1:], i, on)
                if inside:
                    assert any(other == on[i] for other, _ in empty.get(on[i - 1], [])), (
                        hole_bridge, on[i - 1], on[i]
                    )
            checked += 1
    assert checked, "the case must carry an image chain across a foreign stable branch"
