"""
The iterated homotopy partition: the homotopy partition cut by image bridges.

The empty-stretch rule (author, 2026-09-23): the holes' stretches and the
footprints of the holes' image lobes (the base ``(x_0, x_m)`` of a hole
bridge's image chain and its chords ``(x_1, x_2), ...``) are EMPTY on both
sides of the manifold; an image bridge's crossing stays with the element AWAY
from the empty stretch it abuts (the flank toward the empty stretch is
opened). The anchor's one-sided hole marks nothing. Existing hole boundaries
stand; the other side is untouched at an image endpoint.
"""

from __future__ import annotations

import dataclasses
import importlib
import logging

from helpers.logs import assert_logged
from minimal_helpers import build_pieces
from tanglepack.topology.PartitionFamily import IteratedHomotopyPartition

#: The module (the package attribute of the same name is the class).
family_module = importlib.import_module("tanglepack.topology.PartitionFamily")


def test_no_image_bridges_reproduces_the_homotopy_family(k10_partitioned):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    bare = dataclasses.replace(pieces.minimal, image_bridge_ids=[], image_chains={})
    iterated = IteratedHomotopyPartition.from_minimal(bare, pieces.homotopy)
    assert iterated.signature() == pieces.homotopy.signature()
    assert not iterated.cuts
    for result in iterated:
        for iv in result.intervals:
            assert iv.cut_by is None and iv.parent_element_id == iv.element_id


def test_a_row_invariant_violation_is_warned_and_skipped(
    k10_partitioned, monkeypatch, caplog
):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    victim = pieces.minimal.image_bridge_ids[0]
    original = family_module.row_of_end

    def twisted(trellis, bridge_id, endpoint):
        row = original(trellis, bridge_id, endpoint)
        if bridge_id == victim and endpoint == "second":
            return "right" if row == "left" else "left"
        return row

    monkeypatch.setattr(family_module, "row_of_end", twisted)
    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.PartitionFamily"):
        iterated = IteratedHomotopyPartition.from_minimal(pieces.minimal, pieces.homotopy)
    skipped = [cut for cut in iterated.cuts if cut.bridge_id == victim]
    assert skipped and all(cut.reason == "row invariant" for cut in skipped)
    assert_logged(caplog, logging.WARNING, "tanglepack.topology.PartitionFamily")
    assert sum(len(r.intervals) for r in iterated) < sum(
        len(r.intervals) for r in pieces.iterated
    )


def test_cuts_record_the_far_end_of_the_empty_stretch(k10_partitioned):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    # The empty-stretch kernel (an allowed private kernel) names, per crossing,
    # the far ends of the stretches it abuts; a cut's partner is one of them.
    empty = IteratedHomotopyPartition._empty_stretches(pieces.minimal, pieces.homotopy)
    applied = [cut for cut in pieces.iterated.cuts if cut.opened is not None]
    assert applied
    for cut in applied:
        assert cut.partner in {other for other, _kind in empty[cut.intersection_id]}
        assert repr(cut)
    for cut in pieces.iterated.cuts:
        if cut.opened is None:
            assert cut.reason in {
                "existing boundary", "row invariant", "no partition",
                "no empty stretch", "already cut",
            }


def test_a_crossing_abutting_no_empty_stretch_is_not_cut(k10_partitioned, monkeypatch, caplog):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    original = family_module.IteratedHomotopyPartition._empty_stretches

    def without_bases(minimal, homotopy):
        empty = original(minimal, homotopy)
        return {
            iid: [(other, kind) for other, kind in abutting if kind != "base"]
            for iid, abutting in empty.items()
        }

    monkeypatch.setattr(
        family_module.IteratedHomotopyPartition, "_empty_stretches", staticmethod(without_bases)
    )
    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.PartitionFamily"):
        iterated = IteratedHomotopyPartition.from_minimal(pieces.minimal, pieces.homotopy)
    skipped = [cut for cut in iterated.cuts if cut.reason == "no empty stretch"]
    assert skipped, "the base ends of the mapped hole's lobe abut nothing else"
    assert_logged(caplog, logging.WARNING, "tanglepack.topology.PartitionFamily")
    assert sum(len(r.intervals) for r in iterated) < sum(
        len(r.intervals) for r in pieces.iterated
    )


def test_describe_mentions_parents_and_cuts(k10_partitioned):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    assert pieces.iterated.describe()
    assert repr(pieces.iterated.cuts[0])
