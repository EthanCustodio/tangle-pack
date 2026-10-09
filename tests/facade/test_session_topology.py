"""TangleSession.topology(): one recompute per (generation, selection, pips).

The topology is a pure function of the numerical state and the pip choice
(D1). The session holds ONE result per selection that every product is a
stage of (SIMP-8), carries the pip choice across generations (API-1) and
reads partitions, holes and pips from the selected fixed points only
(API-2). Freshness against a direct build after every event is
``test_session_caches.py``.
"""

from __future__ import annotations

import pytest

from cases import _blast, build_k10, build_k28, build_nested
from helpers.direct import direct_view, view_of

#: session accessor -> its stage on the Topology.
STAGES: dict[str, str] = {
    "homotopy_partition": "homotopy",
    "bridge_classes": "classes",
    "minimal_trellis": "minimal",
    "iterated_partition": "iterated",
    "dual_graph": "dual_graph",
    "symbolic_dynamics": "dynamics",
}


def test_every_product_is_a_stage_of_one_topology():
    session = build_k10().session
    topology = session.topology()

    assert session.topology() is topology
    for accessor, stage in STAGES.items():
        assert getattr(session, accessor)() is getattr(topology, stage), accessor


def test_a_blast_or_a_pip_change_recomputes():
    case = build_k10()
    session, fp = case.session, case.fixed_point
    first = session.topology()
    # The anchor-nearest other candidate (test_session_caches' rule): an inner
    # one's orbit cut can leave a bridge around the end of the cut stable branch,
    # which the row invariant rejects exactly as on a numerical trim at that pip.
    alternate = min(
        (c for c in session.strong_pip_candidates(fp) if c != session.strong_pip(fp)),
        key=lambda c: float(case.registry[c].stable_cdist),
    )
    session.set_strong_pip(fp, alternate)
    second = session.topology()
    assert second is not first and second.strong_pips == (alternate,)

    case = build_k28(blasts=1)
    session = case.session
    first = session.topology()
    _blast(session, case.zones[0], case.fixed_points)
    assert session.topology() is not first
    assert session.topology().generation == session.workbench.generation


def test_pip_choice_survives_a_blast():
    """API-1: no repin after a blast; the zone's pip is still the chosen one."""
    case = build_k28(blasts=0, through="zones")
    session, fp = case.session, case.fixed_point
    pip = session.strong_pip(fp)
    assert pip is not None

    _blast(session, case.zones[0], case.fixed_points)

    assert session.strong_pip(fp) == pip
    assert session.bridge_classes()


def test_one_snapshot_per_selection():
    session, fp = (case := build_k10()).session, case.fixed_point

    assert session.trellis(fp) is session.trellis([fp])
    assert session.trellis() is session.trellis(fp)
    assert session.topology(fp) is session.topology()


@pytest.mark.parametrize("first", ["all", "inner"])
def test_gathering_follows_the_selection(first: str):
    """API-2: the inner tangle's products hold its own partitions only, whatever
    was asked before."""
    case = build_nested(outer_blasts=0, through="zones")
    session = case.session
    _fp1, fp3 = case.fixed_points
    if first == "all":
        session.topology()
    graph = session.dual_graph(fp3)
    if first == "inner":
        session.topology()

    assert all(key[0] is fp3 for key in session.homotopy_partition(fp3).branch_keys)
    assert view_of("dual_graph", graph) == direct_view(case, "dual_graph", fp3)
