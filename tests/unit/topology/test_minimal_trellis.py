"""
The minimal trellis: hole bridges, images of those without a zeroth hole, nothing else.

These tests pin the WARNING paths (a pair no bridge spans, an unclassed hole
bridge, a hole no bridge spans); the structural invariants (hole bridges kept,
those without a zeroth hole mapped along their image chains, kept nodes, the sparse
arrangement) run on every law case in ``tests/invariants/test_law_iterated.py``.
"""

from __future__ import annotations

import logging
from dataclasses import replace

from helpers.logs import assert_logged
from tanglepack.topology.MinimalTrellis import image_chain, minimal_trellis


def _holes(session) -> list:
    """Every hole of every per-fixed-point trellis the session holds."""
    return [hole for fp in session.fixed_points for hole in session.trellis(fp).holes]


def test_a_pair_with_no_bridge_object_is_skipped_with_a_warning(
    k10_partitioned, monkeypatch, caplog
):
    session, _fp = k10_partitioned
    full = session.trellis()
    victim = session.minimal_trellis().image_bridge_ids[0]
    table = session.bridge_classes()
    original = full.bridge_between

    def hide(a, b):
        if frozenset((a, b)) == frozenset(victim):
            return None
        return original(a, b)

    monkeypatch.setattr(full, "bridge_between", hide)
    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.MinimalTrellis"):
        minimal = minimal_trellis(full, _holes(session), table)
    assert any(frozenset(p) == frozenset(victim) for _parent, p in minimal.skipped_pairs)
    assert not minimal.is_kept(victim)
    assert_logged(caplog, logging.WARNING, "tanglepack.topology.MinimalTrellis")
    assert minimal.describe(), "the report lists the skipped pair"


def test_an_unclassed_hole_bridge_is_kept_but_not_mapped(
    k10_partitioned, monkeypatch, caplog
):
    session, _fp = k10_partitioned
    table = session.bridge_classes()
    victim = next(iter(session.minimal_trellis().image_chains))
    original = table.entry_of

    def forget(bid):
        if frozenset(bid) == frozenset(victim):
            raise KeyError(bid)
        return original(bid)

    monkeypatch.setattr(table, "entry_of", forget)
    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.MinimalTrellis"):
        minimal = minimal_trellis(session.trellis(), _holes(session), table)
    assert victim in minimal.unclassed_hole_bridge_ids
    assert minimal.is_hole_bridge(victim) and victim not in minimal.image_chains
    assert_logged(caplog, logging.WARNING, "tanglepack.topology.MinimalTrellis")


def test_a_hole_with_no_spanning_bridge_contributes_nothing(k10_partitioned, caplog):
    session, _fp = k10_partitioned
    from tanglepack.topology.TopologyResults import Hole

    ghost = Hole(near_intersection_id=0, bounding_ids=(10**6, 10**6 + 1))
    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.MinimalTrellis"):
        minimal = minimal_trellis(
            session.trellis(), _holes(session) + [ghost], session.bridge_classes()
        )
    assert [frozenset(b) for b in minimal.hole_bridge_ids] == [
        frozenset(b) for b in session.minimal_trellis().hole_bridge_ids
    ]
    assert_logged(caplog, logging.WARNING, "tanglepack.topology.MinimalTrellis")


def test_a_hole_bridge_without_a_registered_image_is_unmapped(p3_partitioned):
    """A hole bridge without a zeroth hole whose image chain is not registered is
    listed as unmapped (kept, never mapped) instead of contributing image bridges;
    one carrying a zeroth hole is never mapped at all, so it is never unmapped.

    The unblasted nested tangle has hole bridges with no registered image, and
    they carry zeroth holes; the same holes read as propagated (no pair) make
    them unmapped. The law cases have none.
    """
    session, _fp3, _fp1 = p3_partitioned
    minimal = session.minimal_trellis()
    trellis = minimal.trellis
    zeroth = {
        frozenset(trellis.bridge_between(*hole.bounding_ids).id)
        for hole in trellis.holes
        if hole.pair is not None
    }
    for bid in minimal.unmapped:
        assert frozenset(bid) not in zeroth
    unchained = [bid for bid in minimal.hole_bridge_ids if image_chain(trellis, bid) is None]
    assert unchained, "the unblasted nested tangle has a hole bridge with no registered image"

    propagated = [replace(hole, pair=None) for hole in minimal.holes]
    as_propagated = minimal_trellis(trellis, propagated, minimal.classes)
    assert as_propagated.inert_hole_bridge_ids == []
    assert sorted(as_propagated.unmapped) == sorted(unchained)
    for bid in as_propagated.unmapped:
        assert as_propagated.is_hole_bridge(bid) and bid not in as_propagated.image_chains
    assert as_propagated.describe()
