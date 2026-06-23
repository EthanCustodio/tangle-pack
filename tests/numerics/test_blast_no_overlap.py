"""Single-copy bridge invariant and the loom's handling of already-known bridges.

A bridge is uniquely defined by the two intersections it connects, and the
workbench keeps exactly one copy of each computed bridge. Iterating a bridge always
returns its children, but a child that re-traces curve an existing bridge already
holds resolves to that existing persistent object instead of a duplicate -- most
visibly the bridge anchored at the fixed point, whose forward image re-covers
itself (the fixed point is invariant, cdist 0 -> 0). The blast then recognizes those
already-known bridges and does not re-iterate or double-count them.
"""

from __future__ import annotations

import pytest


def _fixed_point_bridge(workbench, fp):
    """The fp's bridge anchored at the fixed point (smallest unstable cdist root)."""
    bridges = [b for b in workbench.bridges if b.fixed_point is fp]
    return min(bridges, key=lambda b: b.root.get_cdist("unstable"))


@pytest.mark.slow
def test_iterating_fixed_point_bridge_returns_existing_copies(henon_p3_session):
    session, _fp3, fp1, _zone = henon_p3_session
    workbench = session.workbench

    fp_bridge = _fixed_point_bridge(workbench, fp1)
    before = list(workbench.bridges)
    children = workbench.iterate_bridge(fp_bridge)

    # Its forward image re-traces grown curve, so every child resolves to a bridge
    # already in the trellis -- nothing new is added (single-copy invariant).
    assert children, "iterate_bridge should always return the children"
    assert all(child in before for child in children)
    assert len(workbench.bridges) == len(before)


@pytest.mark.slow
def test_workbench_keeps_single_copy_per_bridge(henon_p3_session):
    session, _fp3, fp1, zone = henon_p3_session
    workbench = session.workbench
    session.blast_zone(zone, num_iterations=3, fixed_point=fp1)

    # No two registered bridges on the same manifold share a signature (the two
    # bounding-intersection cdists) -- i.e. there is exactly one copy of each
    # computed bridge. (cdist is per-manifold, so compare within a fixed point.)
    from collections import defaultdict

    by_manifold = defaultdict(list)
    for bridge in workbench.bridges:
        sig = workbench._bridge_signature(bridge)
        if sig is None:
            continue
        key = workbench._manifold_identity(bridge)
        for other in by_manifold[key]:
            assert not workbench._signatures_match(sig, other), (
                f"two registered bridges on one manifold share signature {sig}"
            )
        by_manifold[key].append(sig)


@pytest.mark.slow
def test_blast_recognizes_already_known_bridges(henon_p3_session):
    session, _fp3, fp1, zone = henon_p3_session
    result = session.blast_zone(zone, num_iterations=2, fixed_point=fp1)

    # The fixed-point bridge (and other re-tracing images) yield already-known
    # bridges that the loom recognizes rather than re-iterating.
    assert result.already_known > 0
    total = sum(len(step.already_known) for step in result.steps)
    assert total == result.already_known
