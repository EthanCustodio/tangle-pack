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

import numpy as np
import pytest

from tanglepack import TangleSession


# --------------------------------------------------------------------------- #
# Fixture: the cheapest flow that yields a resonance zone with interior bridges
# --------------------------------------------------------------------------- #
def _k10_map(point):
    k, b = 10, 1
    x, y = point
    return np.stack([y - k + x**2, -b * x], axis=0)


def _k10_map_inverse(point):
    k, b = 10, 1
    x, y = point
    return np.stack([-y / b, x + k - (y**2) / (b**2)], axis=0)


def _k10_jacobian(point):
    k, b = 10, 1
    x, y = point
    return np.array([[2 * x, 1], [-b, 0]])


@pytest.fixture
def k10_session():
    """``(session, fp)`` — the k=10 saddle with intersections and bridges built.

    Function-scoped: both blasting and zone definition mutate the session.
    """
    session = TangleSession(_k10_map, _k10_map_inverse, _k10_jacobian)
    fp = session.construct_fixed_point([4, -4])
    session.orient_eigenvectors(
        fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
    )
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=9)
    session.grow_until_turnaround(fp, "stable")
    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    return session, fp


def _define_zone(session, fp):
    """Trim at a strong pip that genuinely shortens the stable branch.

    The trellis' default strong pip on this flow happens to sit at the largest stable
    canonical distance, so trimming there is a geometric no-op. Pick instead the
    outermost candidate strictly inside the current stable extent, so the trim really
    removes crossings and the restore really has something to rebuild.
    """
    registry = session.workbench.intersection_registry
    trellis = session.trellis(fp)
    trellis.classify_strong_pips()
    outermost = max(registry[i].stable_cdist for i in registry.all_ids())
    inner = [
        c for c in trellis.strong_pip_candidates
        if registry[c].stable_cdist < outermost
    ]
    assert inner, "expected a strong-pip candidate inside the stable extent"
    pip = max(inner, key=lambda c: registry[c].stable_cdist)
    trellis.set_strong_pip(pip)
    session.add_resonance_zones([pip])
    return session.resonance_zones[(fp, 0)]


def _interior_frontier(session, zone, fp):
    """The un-iterated bridges the blast would feed into its first step."""
    return [
        b
        for b in session.workbench.uniiterated_bridges
        if b.fixed_point is fp
        and (point := session._bridge_test_point(b)) is not None
        and zone.contains_point(point)
    ]


def _bridge_signatures(workbench):
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
def test_blast_zone_propagates_assertion_error(k10_session, monkeypatch):
    """An AssertionError from a bridge's forward map escapes even at strict=False."""
    session, fp = k10_session
    zone = _define_zone(session, fp)
    assert _interior_frontier(session, zone, fp), "need a non-empty blast frontier"

    def boom(bridge):
        raise AssertionError("cdist monotonicity violated")

    monkeypatch.setattr(session.workbench, "iterate_bridge", boom)

    with pytest.raises(AssertionError, match="cdist monotonicity violated"):
        session.blast_zone(zone, num_iterations=1, fixed_point=fp, strict=False)


def test_blast_zone_skips_value_error_with_warning(k10_session, monkeypatch, caplog):
    """A ValueError is still a per-bridge skip (with a warning) at strict=False."""
    session, fp = k10_session
    zone = _define_zone(session, fp)
    frontier = _interior_frontier(session, zone, fp)
    assert frontier, "need a non-empty blast frontier"

    def boom(bridge):
        raise ValueError("under-resolved iterate")

    monkeypatch.setattr(session.workbench, "iterate_bridge", boom)

    with caplog.at_level(logging.WARNING, logger="tanglepack.loom.Blast"):
        result = session.blast_zone(
            zone, num_iterations=1, fixed_point=fp, strict=False
        )

    assert result.skipped == len(frontier)
    assert result.completed_iterations == 1
    assert any("skipping bridge" in record.getMessage() for record in caplog.records)


def test_blast_zone_reraises_value_error_when_strict(k10_session, monkeypatch):
    """strict=True still surfaces the tolerated failures."""
    session, fp = k10_session
    zone = _define_zone(session, fp)
    assert _interior_frontier(session, zone, fp), "need a non-empty blast frontier"

    def boom(bridge):
        raise ValueError("under-resolved iterate")

    monkeypatch.setattr(session.workbench, "iterate_bridge", boom)

    with pytest.raises(ValueError, match="under-resolved iterate"):
        session.blast_zone(zone, num_iterations=1, fixed_point=fp, strict=True)


# --------------------------------------------------------------------------- #
# 1.17 — ResonanceZone.restore must rebuild the bridges
# --------------------------------------------------------------------------- #
def test_restore_rebuilds_bridges(k10_session):
    """A trim/restore round trip returns the exact pre-trim bridge set."""
    session, fp = k10_session
    workbench = session.workbench

    before = _bridge_signatures(workbench)
    assert before, "the fixture must produce bridges"

    zone = _define_zone(session, fp)
    trimmed = _bridge_signatures(workbench)
    assert trimmed != before, "the trim must actually change the bridge set"

    zone.restore(workbench)
    after = _bridge_signatures(workbench)

    assert len(after) == len(before)
    for got, want in zip(after, before):
        assert got[0] == pytest.approx(want[0], abs=1e-6)
        assert got[1] == pytest.approx(want[1], abs=1e-6)


def test_restore_keeps_bridges_within_restored_manifolds(k10_session):
    """Every restored bridge lies at or below its manifold's restored tail."""
    session, fp = k10_session
    workbench = session.workbench

    zone = _define_zone(session, fp)
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

    zone = _define_zone(session, fp)
    trimmed = [id(b) for b in workbench.bridges]
    zone.restore(workbench, recompute=False)

    assert [id(b) for b in workbench.bridges] == trimmed
