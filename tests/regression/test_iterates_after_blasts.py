"""
Regression (BCG-3, WP2): blasting leaves the iterate table exact.

Ten period-3 blasts used to leave fold-born image crossings unlinked: the
bridge-endpoint path gated its matches by an ABSOLUTE ``1e-3`` phase-space
radius, smaller than these crossings' own polyline segments, and no
registry-wide inference ran after a blast. A later ``infer_iterates`` then
found links the topology had been built without. No count of the build is
pinned.
"""

from __future__ import annotations

import pytest

from cases import build_period3
from helpers import laws


@pytest.fixture(scope="module")
def p3_ten_blasts():
    """The inner period-3 orbit blasted ten times (~1 s)."""
    return build_period3(blasts=10)


def test_iterate_table_complete_after_blasts(p3_ten_blasts):
    assert laws.check_iterate_table_complete(p3_ten_blasts) > 0
    assert laws.check_iterate_links_pass_the_map_gate(p3_ten_blasts) > 0
    assert laws.check_iterate_table_partial_bijection(p3_ten_blasts) > 0


def test_blasts_leave_nothing_for_inference(p3_ten_blasts):
    assert p3_ten_blasts.workbench.infer_iterates() == 0
