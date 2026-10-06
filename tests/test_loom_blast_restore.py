"""Loom-layer regressions: blast error handling and resonance-zone restore.

Two behaviours are pinned here, both on the cheap ``k=10`` binary-horseshoe flow
(single saddle at ``[4, -4]``, one resonance zone):

* :func:`tanglepack.loom.Blast.blast_zone` must let an ``AssertionError`` escape.
  Assertions are the correctness guard of this codebase (CLAUDE.md); swallowing
  the cdist-monotonicity assertion of ``ManifoldMachine.iterate_manifold`` turned a
  broken iterate into a silently skipped bridge. A ``ValueError`` is still a
  tolerable per-bridge failure under ``strict=False``.
* :meth:`tanglepack.loom.ResonanceZone.ResonanceZone.restore` must rebuild the
  bridges. Restoring only the tails left ``workbench._bridges`` cut against the
  trimmed (shorter) stable manifold, so the tangle disagreed with itself.
"""

from __future__ import annotations

import logging

import pytest

from helpers.logs import assert_logged
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
@pytest.mark.parametrize(
    "error, strict, escapes",
    [
        (AssertionError, False, True),
        (ValueError, False, False),
        (ValueError, True, True),
    ],
    ids=["assertion_escapes_lenient", "value_error_skipped_lenient", "value_error_raised_strict"],
)
def test_blast_zone_error_handling(
    k10_session, monkeypatch, caplog, error: type, strict: bool, escapes: bool
):
    """An AssertionError from a bridge's forward map always escapes; a ValueError
    is a per-bridge skip (with a warning) at ``strict=False`` and raised at
    ``strict=True``."""
    session, fp = k10_session
    zone = define_inner_zone(session, fp)
    frontier = _interior_frontier(session, zone, fp)
    assert frontier, "need a non-empty blast frontier"

    def boom(bridge):
        raise error("the forward map failed")

    monkeypatch.setattr(session.workbench, "iterate_bridge", boom)

    if escapes:
        with pytest.raises(error):
            session.blast_zone(zone, num_iterations=1, fixed_point=fp, strict=strict)
        return

    with caplog.at_level(logging.WARNING, logger="tanglepack.loom.Blast"):
        result = session.blast_zone(zone, num_iterations=1, fixed_point=fp, strict=strict)
    assert result.skipped == len(frontier)
    assert result.completed_iterations == 1
    assert_logged(caplog, logging.WARNING, "tanglepack.loom.Blast")


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
        tail_cdist = manifold.tail.get_cdist("unstable")
        assert bridge.tail.get_cdist("unstable") <= tail_cdist + 1e-9


def test_restore_without_recompute_leaves_bridges_stale(k10_session):
    """recompute=False is the documented batched path: bridges are left alone."""
    session, fp = k10_session
    workbench = session.workbench

    zone = define_inner_zone(session, fp)
    trimmed = [id(b) for b in workbench.bridges]
    zone.restore(workbench, recompute=False)

    assert [id(b) for b in workbench.bridges] == trimmed
