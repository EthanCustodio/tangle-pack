"""TangleSession.topology(): one recompute per (generation, selection, pips).

The topology is a pure function of the numerical state and the pip choice
(D1). The session holds ONE result per selection that every product is a
stage of (SIMP-8), carries the pip choice across generations (API-1) and
reads partitions, holes and pips from the selected fixed points only
(API-2). Freshness against a direct build after every event is
``test_session_caches.py``.
"""

from __future__ import annotations

import logging

import pytest

from cases import _blast, build_inversion, build_k10, build_k28, build_nested
from helpers.direct import direct_view, trellis_view, view_of
from helpers.logs import assert_logged
from tanglepack.loom import define_resonance_zone

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


def test_the_selection_trellis_is_cut_before_any_product():
    """D1, D10: the trellis of several fixed points is cut at their pips when it
    is read, so what it shows never depends on whether a product ran first."""
    case = build_nested(outer_blasts=0, through="zones")
    session = case.session
    zone = max(case.zones, key=lambda zone: zone.area)
    fp = zone.fixed_point
    _blast(session, zone, [fp])
    session.trellis()
    # The anchor-nearest other candidate cuts crossings the zone's pip keeps.
    session.set_strong_pip(fp, min(
        (c for c in session.strong_pip_candidates(fp) if c != session.strong_pip(fp)),
        key=lambda c: float(session.workbench.intersection_registry[c].stable_cdist),
    ))
    before = trellis_view(session.trellis())

    session.topology()

    assert trellis_view(session.trellis()) == before
    assert session.topology().trellis is session.trellis()


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


@pytest.mark.parametrize("builder", [lambda: build_k28(blasts=1), lambda: build_nested()], ids=["k28", "nested"])
def test_a_rebuild_tests_only_the_remembered_pip(builder, monkeypatch):
    """F7 (WP10): after a blast the recompute tests each remembered pip alone, not
    every crossing; the candidates, read later, are those of a fresh classification."""
    from tanglepack import Trellis
    from tanglepack.topology import StrongPip

    case = builder()
    session = case.session
    session.topology()
    _blast(session, case.zones[0], case.fixed_points)

    tested = []
    is_strong_pip = StrongPip.is_strong_pip
    monkeypatch.setattr(
        StrongPip, "is_strong_pip", lambda trellis, iid, **kw: tested.append(iid) or is_strong_pip(trellis, iid, **kw)
    )
    session.topology()
    remembered = [fp for fp in case.fixed_points if session.chosen_strong_pip(fp) is not None]
    assert remembered and len(tested) == len(remembered)
    monkeypatch.undo()

    for fp in case.fixed_points:
        fresh = Trellis.from_workbench(session.workbench, [fp])
        assert session.strong_pip_candidates(fp) == fresh.classify_strong_pips(choose_default=False)


def test_a_zone_at_a_non_candidate_keeps_the_chosen_pip(caplog):
    """B8: a zone whose pip is no strong-pip candidate is stored, warns once and
    leaves the chosen pip alone, so the topology still runs. The non-candidate
    lies on the other stable branch, so the zone does not trim the chosen pip away."""
    case = build_inversion(blasts=2)
    session, fp = case.session, case.fixed_point
    chosen = session.chosen_strong_pip(fp)
    candidates = set(session.strong_pip_candidates(fp))
    registry = case.registry
    outermost = max(
        (
            iid
            for iid in session.trellis(fp).own_intersection_ids
            if not registry[iid].is_synthetic
            and iid not in candidates
            and registry[iid].manifold_b_key != chosen.manifold_b_key
        ),
        key=lambda iid: float(registry[iid].stable_cdist),
    )

    with caplog.at_level(logging.WARNING, logger="tanglepack.loom.TangleSession"):
        zone = session.resonance_zone(outermost)

    assert_logged(caplog, logging.WARNING, "tanglepack.loom.TangleSession", count=1)
    assert session.resonance_zones[zone.key] is zone
    assert session.chosen_strong_pip(fp) is chosen
    session.topology()


def test_a_forgotten_pip_gives_the_default_back():
    """B8: forgetting a chosen pip uncuts its trellis; the next topology takes the default."""
    case = build_k10()
    session, fp = case.session, case.fixed_point
    default = session.strong_pip(fp)
    session.set_strong_pip(fp, min(c for c in session.strong_pip_candidates(fp) if c != default))

    session.forget_strong_pip(fp)
    session.topology()

    assert session.trellis(fp).strong_pip == default
    assert session.chosen_strong_pip(fp) is case.registry[default]


def test_a_stale_pip_error_names_the_remedy(k10_session):
    """B8: the topology's refusal of a pip trimmed away names ``forget_strong_pip``."""
    session, fp = k10_session
    session.classify_strong_pips()
    q0 = session.strong_pip(fp)
    define_resonance_zone(session.workbench, session.trellis(fp).iterate(q0, 1))

    with pytest.raises(ValueError, match="no longer a strong-pip candidate") as raised:
        session.topology()
    assert "forget_strong_pip" in str(raised.value)
