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

import cases
from helpers.logs import assert_logged
from tanglepack.topology.PartitionFamily import IteratedHomotopyPartition

#: The module (the package attribute of the same name is the class).
family_module = importlib.import_module("tanglepack.topology.PartitionFamily")


def test_no_image_bridges_reproduces_the_homotopy_family(k10_partitioned):
    session, _fp = k10_partitioned
    homotopy = session.homotopy_partition()
    bare = dataclasses.replace(session.minimal_trellis(), image_bridge_ids=[], image_chains={})
    iterated = IteratedHomotopyPartition.from_minimal(bare, homotopy)
    assert iterated.signature() == homotopy.signature()
    assert not iterated.cuts
    for result in iterated:
        for iv in result.intervals:
            assert iv.cut_by is None and iv.parent_element_id == iv.element_id


def test_a_row_invariant_violation_is_warned_and_skipped(
    k10_partitioned, monkeypatch, caplog
):
    session, _fp = k10_partitioned
    minimal, homotopy = session.minimal_trellis(), session.homotopy_partition()
    unpatched = session.iterated_partition()
    victim = minimal.image_bridge_ids[0]
    original = family_module.row_of_end

    def twisted(trellis, bridge_id, endpoint):
        row = original(trellis, bridge_id, endpoint)
        if bridge_id == victim and endpoint == "second":
            return "right" if row == "left" else "left"
        return row

    monkeypatch.setattr(family_module, "row_of_end", twisted)
    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.PartitionFamily"):
        iterated = IteratedHomotopyPartition.from_minimal(minimal, homotopy)
    skipped = [cut for cut in iterated.cuts if cut.bridge_id == victim]
    assert skipped and all(cut.reason == "row invariant" for cut in skipped)
    assert_logged(caplog, logging.WARNING, "tanglepack.topology.PartitionFamily")
    assert sum(len(r.intervals) for r in iterated) < sum(
        len(r.intervals) for r in unpatched
    )


def test_cuts_record_the_far_end_of_the_empty_stretch(k10_partitioned):
    session, _fp = k10_partitioned
    iterated = session.iterated_partition()
    # The empty-stretch kernel (an allowed private kernel) names, per crossing,
    # the far ends of the stretches it abuts; a cut's partner is one of them.
    empty = IteratedHomotopyPartition._empty_stretches(
        session.minimal_trellis(), session.homotopy_partition()
    )
    applied = [cut for cut in iterated.cuts if cut.opened is not None]
    assert applied
    for cut in applied:
        assert cut.partner in {other for other, _kind in empty[cut.intersection_id]}
        assert repr(cut)
    for cut in iterated.cuts:
        if cut.opened is None:
            assert cut.reason in {
                "existing boundary", "row invariant", "no partition",
                "unmarked lobe", "no empty stretch", "already cut",
            }


def test_a_crossing_abutting_no_empty_stretch_is_not_cut(k10_partitioned, monkeypatch, caplog):
    session, _fp = k10_partitioned
    minimal, homotopy = session.minimal_trellis(), session.homotopy_partition()
    unpatched = session.iterated_partition()
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
        iterated = IteratedHomotopyPartition.from_minimal(minimal, homotopy)
    skipped = [cut for cut in iterated.cuts if cut.reason == "no empty stretch"]
    assert skipped, "the base ends of the mapped hole's lobe abut nothing else"
    assert_logged(caplog, logging.WARNING, "tanglepack.topology.PartitionFamily")
    assert sum(len(r.intervals) for r in iterated) < sum(
        len(r.intervals) for r in unpatched
    )


def test_describe_mentions_parents_and_cuts(k10_partitioned):
    session, _fp = k10_partitioned
    iterated = session.iterated_partition()
    assert iterated.describe()
    assert repr(iterated.cuts[0])


def _image_crossings(trellis, bridge_id):
    """Every crossing of the bridge's forward image arc: its ends' +1 iterates
    and the crossings of the advanced unstable branch between them."""
    images = [trellis.iterate(end, 1) for end in bridge_id]
    if any(image not in trellis for image in images):
        return set()
    key = trellis.unstable_key_of(bridge_id)
    branch = trellis.branch(key[0].advance_key(key, 1))
    lo, hi = sorted(trellis.intersection(i).unstable_cdist for i in images)
    return set(images) | {
        i for i in branch.intersection_ids if lo <= trellis.intersection(i).unstable_cdist <= hi
    }


def test_an_unmarked_lobe_is_skipped_by_design(caplog):
    """A crossing on the image lobe of a hole with no stretch of its own (the
    anchor's one-sided hole, or a lobe holding other holes) marks nothing by the
    author's 2026-09-23 rule: it is skipped as ``"unmarked lobe"``, at INFO."""
    session = cases.build_inversion(blasts=2).session
    minimal, homotopy = session.minimal_trellis(), session.homotopy_partition()
    trellis = minimal.trellis
    stretches = {
        frozenset((iv.lo_id, iv.hi_id))
        for result in homotopy
        for iv in result.intervals
        if iv.lo_id != iv.hi_id and not iv.closed_lo and not iv.closed_hi
    }
    on_unmarked = set().union(
        *(
            _image_crossings(trellis, hole_bridge)
            for hole_bridge in minimal.hole_bridge_ids
            if frozenset(hole_bridge) not in stretches
        )
    )
    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.PartitionFamily"):
        iterated = IteratedHomotopyPartition.from_minimal(minimal, homotopy)
    unmarked = [cut for cut in iterated.cuts if cut.reason == "unmarked lobe"]
    assert unmarked, "precondition: the anchor hole bridge's image lobe is cut by an image bridge"
    assert all(cut.intersection_id in on_unmarked for cut in unmarked)
    assert not [cut for cut in iterated.cuts if cut.reason == "no empty stretch"]
    assert not [
        r for r in caplog.records
        if r.name == "tanglepack.topology.PartitionFamily" and r.levelno >= logging.WARNING
    ]
