"""Single-copy bridge invariant and the loom's handling of already-known bridges.

A bridge is uniquely defined by the two intersections it connects, and the
workbench keeps exactly one copy of each computed bridge. Iterating a bridge always
returns its children, but a child that re-traces curve an existing bridge already
holds resolves to that existing persistent object instead of a duplicate -- most
visibly the bridge anchored at the fixed point, whose forward image re-covers
itself (the fixed point is invariant, cdist 0 -> 0). The blast then recognizes those
already-known bridges and does not re-iterate or double-count them. (The
single-copy and no-overlap laws on every law case, blasted ones included, run in
``tests/invariants/test_law_bridges.py``.)
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
def test_blast_recognizes_already_known_bridges(henon_p3_session):
    """Already-known children are counted once and never re-iterated.

    The fixed-point bridge (and other re-tracing images) yield bridges the blast
    has already enqueued; the loom must recognise them rather than feed them to
    a later step, so no bridge is a parent twice and no already-known child is
    a later parent.
    """
    session, _fp3, fp1, zone = henon_p3_session
    result = session.blast_zone(zone, num_iterations=2, fixed_point=fp1)

    assert result.already_known > 0
    total = sum(len(step.already_known) for step in result.steps)
    assert total == result.already_known

    parents = [b for step in result.steps for b in step.interior_parents]
    assert len({id(b) for b in parents}) == len(parents), "a bridge was iterated twice"
    for index, step in enumerate(result.steps):
        later = {id(b) for s in result.steps[index + 1:] for b in s.interior_parents}
        for child in step.already_known:
            assert id(child) not in later, f"already-known {child.id} was re-iterated"
            assert id(child) not in {id(b) for b in step.kept_interior}
