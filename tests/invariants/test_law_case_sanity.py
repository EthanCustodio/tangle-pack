"""Law tier: every law case really is the case it claims to be.

Each case builds through the full pipeline (every product non-empty), its
``k_value`` is ``2 * period`` exactly on an inversion point (the k=10 saddle
at ``(-2.3166, 2.3166)``: both eigenvalues negative, det J = +1, NOT
orientation reversing), every branch of every cycle is built and grown, and
every unstable and stable branch crosses something. The ``b = -1`` orientation-reversing
placeholder is a known issue: ``construct_fixed_point`` raises ``ValueError``
until det J < 0 is supported.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from cases import BUILDERS, issue_marks
from helpers import laws
from helpers.law_tier import assert_non_vacuous, build_products, law_test

test_case_builds = law_test(laws.check_case_builds)
test_k_value_matches_inversion = law_test(laws.check_k_value_matches_inversion)
test_every_branch_is_built = law_test(laws.check_every_branch_is_built)
test_every_branch_crosses = law_test(laws.check_every_branch_crosses)


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
    """Every check of ``helpers.laws`` is wired into a ``test_law_*`` module.

    A law written but never parametrized over the cases would check nothing.
    """
    sources = "".join(
        path.read_text() for path in Path(__file__).parent.glob("test_law_*.py")
    )
    checks = [check for layer in laws.LAYERS.values() for check in layer] + list(laws.MUTATING)
    missing = [
        check.__name__
        for check in checks
        if not re.search(rf"\blaw_test\(\s*laws\.{check.__name__}\b", sources)
    ]
    assert not missing, missing
