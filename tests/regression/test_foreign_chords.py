"""PK-4 / C19 / B7b: the chords of a lobe image pair in stable order, on every branch.

One statement, two branch kinds. A lobe image holds no anchor, so walking any
stable branch from its anchor enters and leaves it in turn (Jordan
alternation): the stretches inside are the crossings paired in stable order.

* A FOREIGN branch: a hole's image chain can leave its base branch for an
  excursion across another stable branch (another orbit branch or fixed point).
  Before WP10b those crossings were left out of the chord pairing, so a crossing
  such as 44 on ``p3_6_blasts`` abutted no empty stretch and two cuts were
  skipped.
* The BASE branch, the base counting as inside
  (``helpers.laws.base_chord_misreadings``): until B7b its chords paired in the
  order the chain visits them, which differs only where the image spirals back
  across its own base, at a blast step of ``k32_series`` (checked per step in
  ``test_class_count_under_blasting.py``); here the builds, the anchor-cornered
  inversion lobe among them.

The geometric checks below are that statement, measured.
"""

from __future__ import annotations

from functools import cache

import pytest

import cases
from helpers.laws import _lobe_image, base_chord_misreadings
from tanglepack.numerics.geometry import arc_polyline, point_in_polygon
from tanglepack.topology.BridgeClass import _image_chain
from tanglepack.topology.PartitionFamily import IteratedHomotopyPartition

CASES = ["p3_6_blasts", "p3_15", "p3_16", "heteroclinic"]
BASE_CASES = ["k10", "inversion", "inversion_period2", "heteroclinic", "p3_6_blasts"]
#: Builds whose lobes cross their base branch at most once: a visible skip.
NO_BASE_LOBE = {"p3_6_blasts"}

_build = cache(lambda name: cases.BUILDERS[name]())


@pytest.fixture(params=CASES)
def case(request: pytest.FixtureRequest) -> cases.Case:
    return _build(request.param)


def test_every_image_crossing_abuts_an_empty_stretch(case: cases.Case) -> None:
    cuts = case.session.iterated_partition().cuts
    assert cuts
    assert not [cut for cut in cuts if cut.reason == "no empty stretch"]


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
                      float(manifold.tail.cdist)]
            for i, (lo, hi) in enumerate(zip(bounds, bounds[1:])):
                stretch = arc_polyline(manifold, lo, hi)
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


@pytest.mark.parametrize("name", BASE_CASES)
def test_base_chords_are_inside_the_lobe_image(name: str) -> None:
    """Along a lobe image's own base branch, from its anchor, the gaps between
    its crossings lie outside, inside, ... the image (the base inside), and the
    empty stretches there are its inside gaps (``base_chord_misreadings``)."""
    checked, misreadings = base_chord_misreadings(_build(name).session)
    if name in NO_BASE_LOBE:
        assert not checked, f"{name} now has a base lobe: test it"
        pytest.skip("no stretched hole crosses its base branch twice")
    assert checked and misreadings == (), "\n".join(misreadings)
