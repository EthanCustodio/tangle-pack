"""
Regression (2026-10-02, hole-side-own-blast): each nested resonance zone
blasts only its own tangle.

The outer zone's blasts must leave the inner period-3 tangle exactly as it is
alone, and the order of the two zones' blasts must not matter. Classes are
spelled by their homotopy elements, never by letters (letters differ between
sessions). The period-3 and nested cases come from the tests-only
:mod:`cases` builders (parameters frozen there), built fresh per test. The
names on these cases are ``tests/unit/topology/test_element_naming.py``; the
cross-branch chord rules of the iterated cut are
``tests/unit/topology/test_iterated_cut.py``.
"""

from __future__ import annotations

import logging


from helpers.logs import assert_logged
from cases import Case, build_nested, build_period3


# --------------------------------------------------------------------------- #
# Nested blasts: each zone blasts only its own tangle
# --------------------------------------------------------------------------- #
def _tangle_words(build: Case, fixed_point: object) -> dict:
    """Every class of one tangle mapped to its word, letters spelled out.

    Letters differ between sessions, so a class is named by its two homotopy
    elements (orbit code and short name, no fixed-point letter) and each token
    by its class's name and direction.
    """
    dynamics = build.session.symbolic_dynamics()

    def spelled(bridge_class):
        ends = (bridge_class.source, bridge_class.target)
        names = (dynamics.naming.homotopy_name(ref) for ref in ends)
        return tuple((name.orbit_code, name.short_text) for name in names)

    return {
        spelled(cd.bridge_class): [(spelled(s.bridge_class), s.direction) for s in cd.symbols]
        for cd in dynamics.classes.values()
        if cd.bridge_class.source.branch_key[0] is fixed_point
    }


def test_nested_inner_words_are_the_period3_words():
    """The outer blasts leave the inner tangle exactly as it is alone: the same
    words and the same bridges (2026-10-02: the outer zone's blasts no longer
    iterate the period-3 bridges)."""
    alone = build_period3(blasts=4)
    nested = build_nested(outer_blasts=2, inner_blasts=4)
    _outer, inner = nested.fixed_points
    (fp3,) = alone.fixed_points
    words = _tangle_words(alone, fp3)
    assert words and _tangle_words(nested, inner) == words
    assert len(nested.session.trellis(inner).bridges) == len(
        alone.session.trellis(fp3).bridges
    )


def test_nested_blast_order_does_not_matter(caplog):
    with caplog.at_level(logging.WARNING, logger="tanglepack.loom.Blast"):
        inner_first = build_nested(outer_blasts=4, inner_blasts=4)
    assert_logged(caplog, logging.WARNING, "tanglepack.loom.Blast", count=0)

    outer_first = build_nested(outer_blasts=4)
    session = outer_first.session
    inner_zone = min(session.resonance_zones.values(), key=lambda zone: zone.area)
    for _ in range(4):
        session.blast_zone(
            inner_zone, 1, fixed_point=[inner_zone.fixed_point], min_separation=1e-4
        )
        session.classify_strong_pips()
        for fp, pip in zip(outer_first.fixed_points, outer_first.pips):
            session.set_strong_pip(fp, pip)
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()

    for one, other in zip(inner_first.fixed_points, outer_first.fixed_points):
        words = _tangle_words(inner_first, one)
        assert words and _tangle_words(outer_first, other) == words
        assert len(inner_first.session.trellis(one).bridges) == len(
            session.trellis(other).bridges
        )
