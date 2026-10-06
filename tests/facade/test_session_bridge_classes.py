"""A.4 -- TangleSession.bridge_classes: describe, zones, warnings, the alphabet.

:meth:`TangleSession.bridge_classes` gathers the per-fixed-point partitions,
letters the active classes from its persistent :class:`BridgeAlphabet` and
records the resonance zone each class lies in. This module checks the
describe report, the zone annotation after trimming at a pip, the warnings
(a class straddling zones, an unpartitioned branch) and the alphabet on its
own. The cache contract and letter stability across rebuilds are
``tests/facade/test_session_caches.py``; session = direct build is
``tests/facade/test_session_equivalence.py``.
"""

from __future__ import annotations

import logging

import matplotlib

matplotlib.use("Agg")  # headless: the session fixtures touch the plotting stack
import pytest

from cases import build_k10
from helpers.logs import assert_logged
from tanglepack.loom.BridgeAlphabet import BridgeAlphabet, letter


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _repartition_at(session, fp, pip):
    """Trim at ``pip``, then redo the topological steps the trim invalidates."""
    zone = session.resonance_zone(pip)
    session.classify_strong_pips()
    session.set_strong_pip(fp, pip)
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()
    return zone


# --------------------------------------------------------------------------- #
# (a) describe and zones
# --------------------------------------------------------------------------- #
def test_describe_reports_letters_and_inertness(k10_partitioned):
    session, _fp = k10_partitioned

    report = session.describe_bridge_classes()

    assert report
    for entry in session.bridge_classes():
        assert entry.name in report


def test_active_class_lies_in_the_zone_after_trimming_at_the_image_pip():
    """Trimming at f(q0) makes the zone-side lobes interior and the exterior
    lobes exterior; the classes record that.

    Built at ten unstable steps (the ``scripts/henon_k10_dual_graph.py``
    recipe): at the fixture's nine, trimming at f(q0) leaves only four
    crossings and the partition collapses to singletons, so the two-class
    picture does not exist there.
    """
    case = build_k10(unstable_steps=10, through="pips")
    session, fp = case.session, case.fixed_point
    trellis = session.trellis(fp)
    pip = trellis.iterate(trellis.strong_pip, 1)
    assert pip is not None and pip in trellis.strong_pip_candidates

    zone = _repartition_at(session, fp, pip)
    table = session.bridge_classes()

    assert len(table) == 2
    active, inert = table.active[0], table.inert[0]
    assert active.letter == "a" and inert.letter is None
    assert active.zone_key == zone.key
    assert inert.zone_key is None
    assert active.describe()
    assert len(active.members) == 4 and len(inert.loops) == 1


def test_a_class_straddling_zones_warns_and_gets_no_zone(k10_partitioned, caplog):
    session, _fp = k10_partitioned
    table = session.bridge_classes()
    active = table.active[0]
    target = active.bridge_ids[0]

    class _FakeZone:
        key = ("fake", 0)

    def classify(bridge):
        return _FakeZone() if bridge.id == target else None

    session.classify_bridge = classify
    with caplog.at_level(logging.WARNING, logger="tanglepack.loom.TangleSession"):
        rebuilt = session.bridge_classes(rebuild=True)

    assert rebuilt.active[0].zone_key is None
    assert_logged(caplog, logging.WARNING, "tanglepack.loom.TangleSession")


# --------------------------------------------------------------------------- #
# (b) no partitions: ValueError from A.3, and the warning
# --------------------------------------------------------------------------- #
def test_no_partitions_raises_and_warns(k10_session, caplog):
    session, _fp = k10_session

    with caplog.at_level(logging.WARNING, logger="tanglepack.loom.TangleSession"):
        with pytest.raises(ValueError):
            session.bridge_classes()

    assert_logged(caplog, logging.WARNING, "tanglepack.loom.TangleSession")


# --------------------------------------------------------------------------- #
# (c) the alphabet on its own
# --------------------------------------------------------------------------- #
def test_letter_sequence():
    assert [letter(i) for i in range(4)] == ["a", "b", "c", "d"]
    assert letter(25) == "z"
    assert letter(26) == "aa"
    assert letter(27) == "ab"
    assert letter(26 + 26) == "ba"
    with pytest.raises(ValueError):
        letter(-1)


def test_alphabet_reuses_assigns_and_resets():
    alphabet = BridgeAlphabet()
    first, second = object(), object()

    assert alphabet.letter_for(first) == "a"
    assert alphabet.letter_for(second) == "b"
    assert alphabet.letter_for(first) == "a"
    assert len(alphabet) == 2 and first in alphabet
    assert alphabet.class_of("b") is second
    with pytest.raises(KeyError):
        alphabet.class_of("z")

    alphabet.reset()
    assert len(alphabet) == 0
    assert alphabet.letter_for(second) == "a"


