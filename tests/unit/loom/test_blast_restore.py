"""Loom-layer regressions: blast error handling and resonance-zone restore.

Two behaviours are pinned here, both on the cheap ``k=10`` binary-horseshoe flow
(single saddle at ``[4, -4]``, one resonance zone):

* :func:`tanglepack.loom.Blast.blast_zone` must let every iterate failure escape.
  Assertions are the correctness guard of this codebase (CLAUDE.md); swallowing
  the cdist-monotonicity assertion of ``ManifoldMachine.iterate_manifold`` turned a
  broken iterate into a silently skipped bridge, and the ``AttributeError`` it
  still skipped until WP10 was NUM-9.
* :meth:`tanglepack.loom.ResonanceZone.ResonanceZone.restore` must rebuild the
  bridges. Restoring only the tails left ``workbench._bridges`` cut against the
  trimmed (shorter) stable manifold, so the tangle disagreed with itself.

A third, on ``k=2.8`` two blasts: after stable growth (which re-finds every
physical end, un-trimming the zone) ``restore`` never shortens a branch (C22).
"""

from __future__ import annotations

import pytest

import cases
from helpers.zones import define_inner_zone
from tanglepack.numerics.geometry import polyline_midpoint


# --------------------------------------------------------------------------- #
# Helpers (the session is the shared ``k10_session`` fixture: function scoped,
# since both blasting and zone definition mutate it)
# --------------------------------------------------------------------------- #
def _interior_frontier(session, zone, fp):
    """The un-iterated bridges the blast would feed into its first step."""
    return [
        b
        for b in session.workbench.uniiterated_bridges
        if b.fixed_point is fp
        and (point := polyline_midpoint(b.get_point_array())) is not None
        and zone.contains_point(point)
    ]


def _bridge_endpoint_cdists(workbench):
    """Each bridge as its ``(first, second)`` unstable canonical distances.

    Registry ids are renumbered by a recompute, but the canonical distances of the
    bounding crossings are not, so this is the identity that must survive a
    trim/restore round trip.
    """
    registry = workbench.intersection_registry
    return sorted(
        (
            round(registry[b.first_intersection].unstable_cdist, 7),
            round(registry[b.second_intersection].unstable_cdist, 7),
        )
        for b in workbench.bridges
        if b.first_intersection is not None and b.second_intersection is not None
    )


# --------------------------------------------------------------------------- #
# 1.5 — blast_zone must not swallow AssertionError
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("error", [AssertionError, ValueError, AttributeError])
def test_every_iterate_failure_escapes(k10_session, monkeypatch, error: type):
    """Any error raised by a bridge's forward map escapes ``blast_zone`` (NUM-10):
    a failure to iterate is a bug, never a skipped bridge."""
    session, fp = k10_session
    zone = define_inner_zone(session, fp)
    assert _interior_frontier(session, zone, fp), "need a non-empty blast frontier"

    def boom(bridge):
        raise error("the forward map failed")

    monkeypatch.setattr(session.workbench, "iterate_bridge", boom)
    with pytest.raises(error):
        session.blast_zone(zone, num_iterations=1, fixed_point=fp)


# --------------------------------------------------------------------------- #
# 1.17 — ResonanceZone.restore must rebuild the bridges
# --------------------------------------------------------------------------- #
def test_restore_rebuilds_bridges(k10_session):
    """A trim/restore round trip returns the exact pre-trim bridge set."""
    session, fp = k10_session
    workbench = session.workbench

    before = _bridge_endpoint_cdists(workbench)
    assert before, "the fixture must produce bridges"

    zone = define_inner_zone(session, fp)
    trimmed = _bridge_endpoint_cdists(workbench)
    assert trimmed != before, "the trim must actually change the bridge set"

    zone.restore(workbench)
    after = _bridge_endpoint_cdists(workbench)

    assert len(after) == len(before)
    for got, want in zip(after, before):
        assert got[0] == pytest.approx(want[0], abs=1e-6)
        assert got[1] == pytest.approx(want[1], abs=1e-6)


def test_restore_preserves_intersection_ids(k10_session):
    """``restore`` must not renumber the registry underneath the ids callers hold.

    ``BridgeId`` is a pair of registry ids and ``rebuild_bridges`` carries per-bridge
    metadata across by that pair, so a restore that renumbered from zero would hand
    every held id to an unrelated crossing.
    """
    session, fp = k10_session
    workbench = session.workbench

    zone = define_inner_zone(session, fp)
    before = {
        iid: (
            round(ix.unstable_cdist, 9),
            round(ix.stable_cdist, 9),
        )
        for iid, ix in workbench.intersection_registry
    }
    assert before

    zone.restore(workbench)

    registry = workbench.intersection_registry
    for iid, cdists in before.items():
        assert iid in registry, f"id {iid} was renumbered away by restore"
        ix = registry[iid]
        assert (round(ix.unstable_cdist, 9), round(ix.stable_cdist, 9)) == cdists, (
            f"id {iid} now names a different crossing"
        )


def test_restore_keeps_bridges_within_restored_manifolds(k10_session):
    """Every restored bridge lies at or below its manifold's restored tail."""
    session, fp = k10_session
    workbench = session.workbench

    zone = define_inner_zone(session, fp)
    zone.restore(workbench)

    for bridge in workbench.bridges:
        manifold = workbench.manifolds[
            (bridge.fixed_point, "unstable", 0, bridge.branch_index or 0)
        ]
        tail_cdist = manifold.tail.cdist
        assert bridge.tail.cdist <= tail_cdist + 1e-9


def test_restore_without_recompute_leaves_bridges_stale(k10_session):
    """recompute=False is the documented batched path: bridges are left alone."""
    session, fp = k10_session
    workbench = session.workbench

    zone = define_inner_zone(session, fp)
    trimmed = [id(b) for b in workbench.bridges]
    zone.restore(workbench, recompute=False)

    assert [id(b) for b in workbench.bridges] == trimmed


# --------------------------------------------------------------------------- #
# C22 — growth un-trims a zone; restore never shortens a branch
# --------------------------------------------------------------------------- #
def test_restore_after_stable_growth_never_shortens():
    """Stable growth after a zone trim re-finds every branch's physical end (the
    trim is a view: the topology is cut at the pip orbit, D10), so a later
    ``restore`` must not cut the grown branch back to its pre-zone tail (C22)."""
    case = cases.build_k28(blasts=2, through="zones")
    workbench, fp, zone = case.session.workbench, case.fixed_points[0], case.zones[0]
    workbench.grow_n_times(fp, "stable", num_iterations=1)
    grown = {key: workbench.manifolds[key].tail.cdist for key in zone.previous_tails}
    assert any(
        grown[key] > tail.cdist for key, tail in zone.previous_tails.items()
    ), "precondition: the growth must run past a pre-zone tail"

    zone.restore(workbench)

    for key, before in grown.items():
        manifold = workbench.manifolds[key]
        _head, _store, rows = manifold._rows(stop_at_final=False)
        assert manifold.tail.cdist >= before
        assert manifold.ids[-1] == rows[-1]
