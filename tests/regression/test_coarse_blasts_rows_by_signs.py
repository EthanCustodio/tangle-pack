"""
Regression (stress campaign E11, NUM-1): rows and openings come from the signs.

Six blasts of the period-3 zone refined at ``area_cutoff = 1e-5`` leave hole
bridges of 2-4 nodes. Read on such a polyline, the "approach" of a bridge end is
the vector to the far node, so the geometric row of an end can disagree with
its crossing sign even though the two signs alternate (the combinatorial rows
agree). The partition used to decide openings and the row assertion on that
geometry and raised "approaches its two defining crossings from opposite
sides". Under D6 the signs decide and the geometry only flags the bridge
(one WARNING per ``punch_holes``). Those bridge ends are below the resolution
threshold, so the build turns clause 1 of N off (``build_period3_coarse_blasts``):
under the default threshold they are no crossings (Q13 a, B1). No count of the
build is pinned.
"""

from __future__ import annotations

import logging

from cases import BUILDERS
from helpers import laws
from tanglepack.topology import bridge_row_diagnosis

_LOGGER = "tanglepack.topology.diagnostics"


def _warnings(caplog) -> int:
    """WARNING records of the diagnostics logger captured so far."""
    return sum(
        1 for r in caplog.records if r.name == _LOGGER and r.levelno == logging.WARNING
    )


def test_coarse_blasts_rows_by_signs(caplog):
    case = BUILDERS["p3_coarse_blasts"]()
    (fp,) = case.fixed_points
    trellis = case.session.trellis(fp)

    # The build must still exercise the defect: a hole bridge whose signs
    # alternate while its polyline reads opposite rows.
    diagnoses = [
        bridge_row_diagnosis(trellis, trellis.bridge_between(*hole.bounding_ids))
        for hole in trellis.holes
    ]
    assert any(
        d is not None and d.signs_alternate and d.geometry_disagrees for d in diagnoses
    ), "the builder no longer makes an under-resolved hole bridge"

    assert case.session.symbolic_dynamics().classes
    assert laws.check_partition_unique_owner(case) > 0
    assert laws.check_openings_follow_signs(case) > 0
    assert laws.check_bridge_rows_consistent(case) > 0

    # The geometric reading is a flag: punching with it logs more than without.
    caplog.set_level(logging.WARNING, logger=_LOGGER)
    before = _warnings(caplog)
    trellis.punch_holes(diagnose=False)
    quiet = _warnings(caplog) - before
    trellis.punch_holes()
    assert _warnings(caplog) - before > 2 * quiet
