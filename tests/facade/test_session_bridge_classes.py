"""A.4 -- TangleSession.bridge_classes: describe, zones, warnings, the alphabet.

:meth:`TangleSession.bridge_classes` is a stage of the session's topology: it
letters the active classes from its persistent :class:`BridgeAlphabet` and
records the resonance zone each class lies in. This module checks the
describe report, the zone annotation after trimming at a pip, the warnings
(a class straddling zones, a pip trimmed away), that the products need no
manual stage call, and the alphabet on its own. The cache contract and letter stability across rebuilds are
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
from tanglepack.loom import define_resonance_zone
from tanglepack.loom.BridgeAlphabet import BridgeAlphabet


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

    assert table.active, "trimming at f(q0) keeps the zone's own lobes"
    for entry in table:
        assert entry.zone_key in (zone.key, None)
        assert (entry.letter is not None) == entry.active
        if entry.active:
            assert entry.zone_key == zone.key, "an active class lies in the zone"
        if entry.loops:
            assert entry.inert, "a class holding a loop is inert"
        assert entry.describe()
    assert any(entry.zone_key is None for entry in table), "the exterior lobes are recorded"


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
    with caplog.at_level(logging.WARNING, logger="tanglepack.loom.Topology"):
        rebuilt = session.bridge_classes(rebuild=True)

    assert rebuilt.active[0].zone_key is None
    assert_logged(caplog, logging.WARNING, "tanglepack.loom.Topology")


# --------------------------------------------------------------------------- #
# (b) D1: the products need only the numerics and a pip
# --------------------------------------------------------------------------- #
def test_products_need_only_the_numerics(k10_session):
    """No manual stage call: the default pip is taken, remembered, and the
    holes and partition are recomputed from it (D1, API-1)."""
    session, fp = k10_session

    assert session.bridge_classes()
    pip = session.strong_pip(fp)
    assert pip is not None
    assert session.topology().strong_pips == (pip,)


def test_a_pip_trimmed_away_warns_and_raises(k10_session, caplog):
    """A remembered pip that is no longer a candidate is never replaced silently."""
    session, fp = k10_session
    session.classify_strong_pips()
    q0 = session.strong_pip(fp)
    # Trimmed at f(q0) behind the session's back: q0 lies beyond the cut.
    define_resonance_zone(session.workbench, session.trellis(fp).iterate(q0, 1))

    with caplog.at_level(logging.WARNING, logger="tanglepack.loom.TangleSession"):
        with pytest.raises(ValueError, match="no longer a strong-pip candidate"):
            session.bridge_classes()

    assert_logged(caplog, logging.WARNING, "tanglepack.loom.TangleSession")


def test_a_pip_trimmed_away_is_not_replaced_by_classifying(k10_session, caplog):
    """Classifying again after the trim warns once and leaves the pip unset:
    the default never takes the chosen pip's place."""
    session, fp = k10_session
    session.classify_strong_pips()
    q0 = session.strong_pip(fp)
    define_resonance_zone(session.workbench, session.trellis(fp).iterate(q0, 1))

    with caplog.at_level(logging.WARNING, logger="tanglepack.loom.TangleSession"):
        session.classify_strong_pips()

    assert_logged(caplog, logging.WARNING, "tanglepack.loom.TangleSession", count=1)
    assert session.strong_pip(fp) is None
    with pytest.raises(ValueError, match="no longer a strong-pip candidate"):
        session.bridge_classes()


@pytest.mark.parametrize("plotter", ["plot_strong_pip", "plot_strong_pip_candidates"])
def test_a_pip_trimmed_away_is_not_replaced_by_plotting(k10_session, plotter):
    """The plot fan-outs classify through the session: drawing after the trim
    leaves the pip unset, so the products still refuse to run."""
    import matplotlib.pyplot as plt

    session, fp = k10_session
    session.classify_strong_pips()
    q0 = session.strong_pip(fp)
    define_resonance_zone(session.workbench, session.trellis(fp).iterate(q0, 1))

    fig, ax = plt.subplots()
    getattr(session, plotter)(ax=ax)
    plt.close(fig)

    assert session.strong_pip(fp) is None
    with pytest.raises(ValueError, match="no longer a strong-pip candidate"):
        session.bridge_classes()


# --------------------------------------------------------------------------- #
# (c) the alphabet on its own (the letter series is ``test_symbolic_dynamics``)
# --------------------------------------------------------------------------- #
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


