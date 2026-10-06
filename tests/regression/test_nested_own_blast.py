"""
Higher-period and nested tangles: the shared case builders and the
cross-branch fixes they need. The circular-zone dual-graph cartoon is
``tests/plotting/``.

The period-3 and nested cases come from the tests-only :mod:`cases` builders
(parameters frozen there), built fresh per test. The cross-branch chord rules
of the iterated cut are ``tests/unit/topology/test_iterated_cut.py``.
"""

from __future__ import annotations

import logging

import pytest

from helpers.logs import assert_logged
from cases import Case, build_nested, build_period3
from tanglepack.topology.ElementNaming import ElementNaming


# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #
@pytest.fixture
def p3_built() -> Case:
    """The period-3 orbit alone."""
    return build_period3()


@pytest.fixture
def nested_built() -> Case:
    """The nested period-1 + period-3 tangle, outer zone blasted twice."""
    return build_nested()


# --------------------------------------------------------------------------- #
# Period 3
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_p3_names_carry_orbit_codes_without_a_letter(p3_built):
    naming = p3_built.session.symbolic_dynamics().naming
    assert not naming.letters
    codes = {name.branch_code for name in naming.names}
    assert codes == {"0.0", "1.0", "2.0"}


# --------------------------------------------------------------------------- #
# Nested
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_nested_names_carry_fixed_point_letters(nested_built):
    outer, inner = nested_built.fixed_points
    naming = nested_built.session.symbolic_dynamics().naming
    assert naming.letters == {id(outer): outer.label, id(inner): inner.label}
    assert {outer.label, inner.label} == {"A", "B"}
    for ref, name in naming.items():
        assert name.fixed_point_letter == ref.fixed_point.label
        assert name.mathtext.startswith(f"${{}}^{{{ref.fixed_point.label}}}")
        assert ref.label.startswith(f"{ref.fixed_point.label}:")
    assert isinstance(naming, ElementNaming)


# --------------------------------------------------------------------------- #
# Nested blasts: each zone blasts only its own tangle
# --------------------------------------------------------------------------- #
def _tangle_words(build, fixed_point) -> dict:
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


@pytest.mark.slow
@pytest.mark.regression
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


@pytest.mark.slow
@pytest.mark.regression
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
