"""Law tier: every law case really is the case it claims to be.

Each case builds through the full pipeline (every product non-empty), its
``k_value`` is ``2 * period`` exactly on an inversion point (the k=10 saddle
at ``(-2.3166, 2.3166)``: both eigenvalues negative, det J = +1, NOT
orientation reversing), every branch of every cycle is built and grown, and
every unstable and stable branch crosses something. The ``b = -1`` orientation-reversing
placeholder is a known issue: ``construct_fixed_point`` raises ``ValueError``
until det J < 0 is supported.

This module also guards the tier's wiring: every check of ``helpers.laws``
runs in some layer, every layer has its ``test_law_<layer>.py`` test, every
known issue of a layer has its own test, and every ``KNOWN_ISSUES`` /
``NOT_APPLICABLE`` key names a real law.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from cases import BUILDERS, KNOWN_ISSUES, NOT_APPLICABLE, issue_marks
from helpers import laws
from helpers.law_tier import assert_non_vacuous, build_products, known_issue_pairs, layer_test

test_case_sanity_laws = layer_test("case_sanity")

#: Known-issue ids that are tests of their own, not ``helpers.laws`` checks.
_STANDALONE_ISSUES = {"orientation_reversing_case_builds", "open_p3_deep_runs"}


@pytest.mark.parametrize(
    "name",
    [
        pytest.param(
            "orientation_reversing",
            marks=issue_marks("orientation_reversing_case_builds", "orientation_reversing"),
        )
    ],
)
def test_orientation_reversing_case_builds(name: str) -> None:
    """The ``b = -1`` map (det J = -1, one negative eigenvalue) runs the full pipeline.

    A placeholder: ``FixedPoint.set_k_value`` rejects the saddle today.
    """
    case = BUILDERS[name]()
    build_products(case)
    assert_non_vacuous(laws.check_case_builds, case)


def test_every_law_runs_in_the_tier() -> None:
    """Every check is in a layer, every layer (and known issue) is wired into a module,
    and the issue registries name real laws.

    A law written but never run over the cases would check nothing.
    """
    layers = {**laws.LAYERS, **laws.MUTATING_LAYERS}
    layered = {check for checks in layers.values() for check in checks}
    defined = {
        obj
        for name, obj in inspect.getmembers(laws, inspect.isfunction)
        if name.startswith("check_") and obj.__module__ == laws.__name__
    }
    assert defined == layered, sorted(c.__name__ for c in defined ^ layered)

    here = Path(__file__).parent
    sources = {path.stem: path.read_text() for path in here.glob("test_law_*.py")}
    everything = "".join(sources.values())
    missing = [layer for layer in layers if f'layer_test("{layer}")' not in everything]
    missing += [
        f"known_issue_test({layer!r})"
        for layer in layers
        if known_issue_pairs(layer) and f'known_issue_test("{layer}")' not in everything
    ]
    assert not missing, missing

    ids = {laws.law_id(check) for check in layered}
    unknown = [
        pair
        for pair in list(KNOWN_ISSUES) + list(NOT_APPLICABLE)
        if pair[0] not in ids and pair[0] not in _STANDALONE_ISSUES
    ]
    assert not unknown, unknown
