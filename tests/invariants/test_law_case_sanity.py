"""Law tier: every law case really is the case it claims to be.

Each case builds through the full pipeline (every product non-empty) at the
``area_cutoff`` its builder froze (never the library default), its
``k_value`` is ``2 * period`` exactly on an inversion point (the k=4 saddle at
``(1 - sqrt5, sqrt5 - 1)``: both eigenvalues negative, det J = +1, NOT
orientation reversing), the map is area preserving at every orbit point
(``|det J| = 1``, so ``construct_fixed_point`` never warns), every branch of
every cycle is built and grown, and every unstable and stable branch crosses
something. The ``b = -1`` orientation-reversing placeholder is a known issue:
det J < 0 is out of scope (C3) and ``construct_fixed_point`` raises
``ValueError``.

This module also guards the tier's wiring: every check of ``helpers.laws``
runs in some layer, every layer is wired into a ``test_law_<layer>.py``
module through ``layer_tests`` (which also gives every set-aside pair its own
test), every ``KNOWN_ISSUES`` / ``NOT_APPLICABLE`` key names a real law and a
real case, no pair is both, and every known issue names the WP that fixes it.
"""

from __future__ import annotations

import inspect
import re
from pathlib import Path

import pytest

from cases import BUILDERS, KNOWN_ISSUES, LAW_CASES, NOT_APPLICABLE, SERIES, issue_marks
from helpers import laws
from helpers.law_tier import assert_non_vacuous, build_products, layer_tests

test_case_sanity_laws, test_case_sanity_set_aside = layer_tests("case_sanity")

#: Known-issue ids that are tests of their own, not ``helpers.laws`` checks.
_STANDALONE_ISSUES = {
    "orientation_reversing_case_builds",
    "open_p3_deep_runs",
    "series_class_count_never_decreases",
    "series_class_count_constant_while_holes_constant",
    "series_alphabet_inside_minimal_trellis",
}


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
    """Every check is in a layer, every layer is wired into a module, and the issue
    registries name real laws and cases, apart, each known issue with its WP (or the
    decision that keeps it).

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
    missing = [layer for layer in layers if f'layer_tests("{layer}")' not in everything]
    assert not missing, missing

    ids = {laws.law_id(check) for check in layered}
    standalone_cases = {*BUILDERS, *SERIES, "hand_built"}
    unknown = [
        (law, case)
        for law, case in (*KNOWN_ISSUES, *NOT_APPLICABLE)
        if not (law in ids and case in LAW_CASES)
        and not (law in _STANDALONE_ISSUES and case in standalone_cases)
    ]
    assert not unknown, unknown
    assert not set(KNOWN_ISSUES) & set(NOT_APPLICABLE)
    untagged = [
        pair
        for pair, issue in KNOWN_ISSUES.items()
        if not re.fullmatch(r"WP\d+|C\d|DEC-\d+", issue.wp)
    ]
    assert not untagged, untagged
