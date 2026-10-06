"""The deep period-3 runs: invariant I1 holds, and the symbolic dynamics runs.

A bridge's unstable branch is carried by ``Bridge.manifold_key``; inferring it
from the endpoint intersections instead was unreliable on a period > 1 orbit
(every branch's anchor sits at unstable cdist 0), and the backward
propagation then picked a container on the wrong branch, so an orbit's holes
flipped ``bridge_side`` partway down the chain. The deep runs (6 blasts, 15
unstable steps) pin I1 -- all holes of one origin share their side of their
bridge -- where a propagated hole's image sub-arc and the nearest vertex of
its whole containing bridge sit on different folds (fixed 2026-10-02).

The same invariants on the law cases (I1, I2, the predicted-branch landing,
the direct-hole side and the openings law) live in
``tests/invariants/test_law_partition.py`` and ``test_law_classes.py``.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless: the session fixture touches the plotting stack
import pytest

from tanglepack.examples.henon_cases import build_period3
from tanglepack.topology import check_holes_share_bridge_side


@pytest.mark.slow
@pytest.mark.regression
@pytest.mark.parametrize(
    "kwargs", [{"blasts": 6}, {"unstable_steps": 15}], ids=["6_blasts", "15_steps"]
)
def test_deep_p3_holes_share_bridge_side_through_symbolic_dynamics(kwargs):
    """I1 holds deep into the tangle, and the symbolic dynamics runs.

    Only that it runs: these tangles are not reliable yet (an unreachable
    class, a virtual ``new1``).
    """
    build = build_period3(**kwargs)
    (fp3,) = build.fixed_points
    trellis = build.session.trellis(fp3)
    check_holes_share_bridge_side(
        trellis.holes, orientation_preserving=trellis.orientation_preserving
    )
    assert build.session.symbolic_dynamics().classes
