"""
The minimal trellis: hole bridges, images of the active ones, nothing else.

These tests pin the WARNING paths (a pair no bridge spans, an unclassed hole
bridge, a hole no bridge spans); the structural invariants (hole bridges kept,
active ones mapped along their image chains, kept nodes, the sparse
arrangement) run on every law case in ``tests/invariants/test_law_iterated.py``.
"""

from __future__ import annotations

import logging

from helpers.logs import assert_logged
from minimal_helpers import build_pieces
from tanglepack.topology.MinimalTrellis import image_chain, minimal_trellis


def test_a_pair_with_no_bridge_object_is_skipped_with_a_warning(
    k10_partitioned, monkeypatch, caplog
):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    full = pieces.full
    victim = pieces.minimal.image_bridge_ids[0]
    original = full.bridge_between

    def hide(a, b):
        if frozenset((a, b)) == frozenset(victim):
            return None
        return original(a, b)

    monkeypatch.setattr(full, "bridge_between", hide)
    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.MinimalTrellis"):
        minimal = minimal_trellis(full, pieces.holes, pieces.table)
    assert any(frozenset(p) == frozenset(victim) for _parent, p in minimal.skipped_pairs)
    assert not minimal.is_kept(victim)
    assert_logged(caplog, logging.WARNING, "tanglepack.topology.MinimalTrellis")
    assert minimal.describe(), "the report lists the skipped pair"


def test_an_unclassed_hole_bridge_is_kept_but_not_mapped(
    k10_partitioned, monkeypatch, caplog
):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    victim = next(iter(pieces.minimal.image_chains))
    original = pieces.table.entry_of

    def forget(bid):
        if frozenset(bid) == frozenset(victim):
            raise KeyError(bid)
        return original(bid)

    monkeypatch.setattr(pieces.table, "entry_of", forget)
    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.MinimalTrellis"):
        minimal = minimal_trellis(pieces.full, pieces.holes, pieces.table)
    assert victim in minimal.unclassed_hole_bridge_ids
    assert minimal.is_hole_bridge(victim) and victim not in minimal.image_chains
    assert_logged(caplog, logging.WARNING, "tanglepack.topology.MinimalTrellis")


def test_a_hole_with_no_spanning_bridge_contributes_nothing(k10_partitioned, caplog):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    from tanglepack.topology.TopologyResults import Hole

    ghost = Hole(coords=(0.0, 0.0), near_intersection_id=0, bounding_ids=(10**6, 10**6 + 1))
    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.MinimalTrellis"):
        minimal = minimal_trellis(pieces.full, pieces.holes + [ghost], pieces.table)
    assert [frozenset(b) for b in minimal.hole_bridge_ids] == [
        frozenset(b) for b in pieces.minimal.hole_bridge_ids
    ]
    assert_logged(caplog, logging.WARNING, "tanglepack.topology.MinimalTrellis")


def test_an_active_hole_bridge_without_a_registered_image_is_unmapped(p3_partitioned):
    """An active hole bridge whose image chain is not registered is listed as
    unmapped (kept, never mapped) instead of contributing image bridges.

    The unblasted nested tangle has such bridges; the law cases do not.
    """
    session, _fp3, _fp1 = p3_partitioned
    minimal = session.minimal_trellis()
    table = session.bridge_classes()
    assert minimal.unmapped, "the unblasted nested tangle has an unmapped hole bridge"
    for bid in minimal.unmapped:
        assert minimal.is_hole_bridge(bid) and bid not in minimal.image_chains
        assert table.entry_of(bid).active
        assert image_chain(minimal.trellis, bid) is None
    assert minimal.describe()
