"""A.4 -- TangleSession.bridge_classes: gathering, caching, letters, zones.

:func:`tanglepack.topology.BridgeClass.bridge_classes` (A.3) needs a trellis to
read bridges from and a list of stable partitions to resolve elements against.
:meth:`TangleSession.bridge_classes` is the session-level wrapper that supplies
both: the trellis is the selection's own (default: the all-fixed-points one,
mirroring :meth:`TangleSession.arrangement`), and the partitions are gathered
from every cached, non-stale per-fixed-point trellis via the private
:meth:`TangleSession._gathered_partitions` helper. The session then letters the
active classes from its persistent :class:`BridgeAlphabet` and records the
resonance zone each class lies in. This module pins that composition, the
``workbench.generation``-keyed cache, letter stability across rebuilds, the
zone annotation after trimming at a pip, the warning logged for an
unpartitioned stable branch, and (on the nested fixture) that no class mixes
the two fixed points.
"""

from __future__ import annotations

import logging

import matplotlib

matplotlib.use("Agg")  # headless: the session fixtures touch the plotting stack
import numpy as np
import pytest

from tanglepack import TangleSession
from tanglepack.loom.BridgeAlphabet import BridgeAlphabet, letter
from tanglepack.topology.BridgeClass import bridge_classes as compute_bridge_classes


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _direct_table(session, fixed_point_list):
    """The same computation session.bridge_classes() should produce, done by hand."""
    partitions = []
    for fp in fixed_point_list:
        partitions.extend(session.trellis(fp).stable_partitions)
    return compute_bridge_classes(session.trellis(), partitions)


def _strip(table):
    """The table's structural content, without the session's annotations."""
    return [
        (entry.bridge_class, tuple(entry.members)) for entry in table
    ]


def _repartition_at(session, fp, pip):
    """Trim at ``pip``, then redo the topological steps the trim invalidates."""
    zone = session.resonance_zone(pip)
    session.classify_strong_pips()
    session.set_strong_pip(fp, pip)
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()
    return zone


# --------------------------------------------------------------------------- #
# (a) the session result matches a direct call over the gathered partitions
# --------------------------------------------------------------------------- #
def test_k10_session_result_matches_a_direct_call(k10_partitioned):
    session, fp = k10_partitioned

    expected = _direct_table(session, [fp])
    actual = session.bridge_classes()

    assert _strip(actual) == _strip(expected)
    assert len(actual) == 2


@pytest.mark.slow
def test_p3_session_result_matches_a_direct_call(p3_partitioned):
    session, fp3, fp1 = p3_partitioned

    expected = _direct_table(session, [fp3, fp1])
    actual = session.bridge_classes()

    assert _strip(actual) == _strip(expected)
    assert len(actual) > 0


def test_session_result_matches_gathered_partitions_helper(k10_partitioned):
    """The private gathering helper session.bridge_classes() actually uses."""
    session, fp = k10_partitioned
    session.bridge_classes()  # populate/observe the same cache the helper reads

    gathered = session._gathered_partitions()
    expected = compute_bridge_classes(session.trellis(), gathered)

    assert _strip(session.bridge_classes()) == _strip(expected)


# --------------------------------------------------------------------------- #
# (b) letters and symbols
# --------------------------------------------------------------------------- #
def test_active_class_is_lettered_a_and_the_inert_one_is_not(k10_partitioned):
    session, _fp = k10_partitioned

    table = session.bridge_classes()
    active, inert = table.active[0], table.inert[0]

    assert active.letter == "a"
    assert inert.letter is None
    assert session.bridge_alphabet.assigned == {active.bridge_class: "a"}
    assert session.bridge_alphabet.class_of("a") == active.bridge_class

    symbols = sorted(table.symbol(bid) for bid in active.bridge_ids)
    assert symbols == ["a", "a", "a^-1", "a^-1"]
    with pytest.raises(ValueError, match="inert"):
        inert.symbol(inert.bridge_ids[0])


def test_letters_are_stable_across_a_rebuild(k10_partitioned):
    session, _fp = k10_partitioned

    first = session.bridge_classes()
    second = session.bridge_classes(rebuild=True)

    assert second is not first
    assert second == first
    assert [e.letter for e in second] == [e.letter for e in first]
    assert len(session.bridge_alphabet) == 1


def test_describe_reports_letters_and_inertness(k10_partitioned):
    session, _fp = k10_partitioned

    report = session.describe_bridge_classes()

    assert "1 active, 1 inert" in report
    assert "a:" in report and "a^-1" in report
    assert "[inert, no zone]" in report and "[active, no zone]" in report


# --------------------------------------------------------------------------- #
# (c) zones
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_active_class_lies_in_the_zone_after_trimming_at_the_image_pip(
    henon_map, henon_map_inverse
):
    """Trimming at f(q0) makes the zone-side lobes interior and the exterior
    lobes exterior; the classes record that.

    Built at ten unstable steps (the ``scripts/henon_k10_dual_graph.py``
    recipe): at the fixture's nine, trimming at f(q0) leaves only four
    crossings and the partition collapses to singletons, so the two-class
    picture does not exist there.
    """
    session = TangleSession(henon_map, henon_map_inverse)
    fp = session.construct_fixed_point([4, -4])
    session.orient_eigenvectors(
        fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
    )
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=10)
    session.grow_until_turnaround(fp, "stable")
    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    session.infer_iterate_table()
    session.classify_strong_pips()
    trellis = session.trellis(fp)
    pip = trellis.iterate(trellis.strong_pip, 1)
    assert pip is not None and pip in trellis.strong_pip_candidates

    zone = _repartition_at(session, fp, pip)
    table = session.bridge_classes()

    assert len(table) == 2
    active, inert = table.active[0], table.inert[0]
    assert active.letter == "a" and inert.letter is None
    assert active.zone_key == zone.key
    assert inert.zone_key is None
    assert f"zone p{fp.period} branch {zone.branch_index}" in active.describe()
    assert len(active.members) == 4 and len(inert.loops) == 1


def test_a_class_straddling_zones_warns_and_gets_no_zone(k10_partitioned, caplog):
    session, _fp = k10_partitioned
    table = session.bridge_classes()
    active = table.active[0]
    target = active.bridge_ids[0]

    class _FakeZone:
        key = ("fake", 0)

    def classify(bridge):
        return _FakeZone() if bridge.id == target else None

    session.classify_bridge = classify
    with caplog.at_level(logging.WARNING, logger="tanglepack.loom.TangleSession"):
        rebuilt = session.bridge_classes(rebuild=True)

    assert rebuilt.active[0].zone_key is None
    assert "straddles resonance zones" in caplog.text


# --------------------------------------------------------------------------- #
# (d) caching
# --------------------------------------------------------------------------- #
def test_cache_hit_returns_the_same_object(k10_partitioned):
    session, _fp = k10_partitioned

    first = session.bridge_classes()

    assert session.bridge_classes() is first


def test_workbench_mutation_invalidates_the_cache(k10_partitioned):
    session, fp = k10_partitioned
    first = session.bridge_classes()

    # Bump the workbench generation and restore a fully partitioned state, so
    # the post-mutation call recomputes successfully rather than hitting the
    # "no partition" ValueError tested separately below.
    session.grow_n_times(fp, "unstable", num_iterations=1)
    session.compute_intersections([fp], preserve_ids=True)
    session.create_bridges(fp)
    session.classify_strong_pips()
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()

    assert session.bridge_classes() is not first


def test_repartitioning_at_the_same_generation_invalidates_the_cache(k10_partitioned):
    """Punching holes and partitioning are Trellis-level mutations: they never
    touch the workbench, so the generation alone cannot see a re-partition.
    Force one by re-running the four partition steps against a different
    strong pip (same registry, same generation, different element boundaries)
    and check the cache does not serve the stale table. A class whose element
    pair survived keeps its letter; a new pair takes a new one."""
    session, fp = k10_partitioned
    first = session.bridge_classes()
    generation_before = session.workbench.generation

    trellis = session.trellis(fp)
    trellis.clear_results()
    trellis.classify_strong_pips()
    alternates = [c for c in trellis.strong_pip_candidates if c != trellis.strong_pip]
    assert alternates, "the k10 fixture must offer more than one strong-pip candidate"
    trellis.set_strong_pip(alternates[0])
    session.compute_pseudoneighbors(fp)
    session.punch_holes(fp)
    session.partition_stable_manifold(fp)

    # The mutation happened entirely on the Trellis: the workbench itself
    # never moved.
    assert session.workbench.generation == generation_before

    second = session.bridge_classes()
    assert second is not first
    letters = {e.bridge_class: e.letter for e in second if e.letter is not None}
    for bridge_class, assigned in letters.items():
        assert session.bridge_alphabet.assigned[bridge_class] == assigned
    assert len(session.bridge_alphabet) >= 1


def test_reads_do_not_invalidate_the_cache(k10_partitioned):
    session, fp = k10_partitioned
    first = session.bridge_classes()

    session.workbench.bridges
    session.trellis(fp).stable_partitions

    assert session.bridge_classes() is first


def test_cache_is_kept_per_fixed_point_selection(k10_partitioned):
    session, fp = k10_partitioned

    whole = session.bridge_classes()
    single = session.bridge_classes(fp)

    # Different cache keys (None vs. the fp itself): both are valid, cheap to
    # ask for independently, and a hit on one must not disturb the other. A
    # cache that collapsed every selector to one key would still pass the
    # "is whole" / "is single" checks below (both would be the same cached
    # object), so pin the keys apart directly.
    assert single is not whole
    assert len(session._bridge_classes) == 2
    assert session.bridge_classes() is whole
    assert session.bridge_classes(fp) is single
    # Same classes either way, so the alphabet handed out one letter.
    assert len(session.bridge_alphabet) == 1


# --------------------------------------------------------------------------- #
# (e) no partitions: ValueError from A.3, and the warning
# --------------------------------------------------------------------------- #
def test_no_partitions_raises_and_warns(k10_session, caplog):
    session, _fp = k10_session

    with caplog.at_level(logging.WARNING, logger="tanglepack.loom.TangleSession"):
        with pytest.raises(ValueError):
            session.bridge_classes()

    assert "missing a partition" in caplog.text.lower()
    assert "stable branch" in caplog.text.lower()
    assert "left" in caplog.text.lower() and "right" in caplog.text.lower()


# --------------------------------------------------------------------------- #
# (f) on p3, no class mixes fixed points
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_p3_no_class_mixes_fixed_points(p3_partitioned):
    session, fp3, fp1 = p3_partitioned

    table = session.bridge_classes()

    seen = set()
    for cls in table.classes:
        owners = {id(cls.source.fixed_point), id(cls.target.fixed_point)}
        assert len(owners) == 1, f"class {cls} mixes two fixed points"
        seen |= owners
    assert seen == {id(fp3), id(fp1)}, "both tangles must contribute classes"


# --------------------------------------------------------------------------- #
# (g) the alphabet on its own
# --------------------------------------------------------------------------- #
def test_letter_sequence():
    assert [letter(i) for i in range(4)] == ["a", "b", "c", "d"]
    assert letter(25) == "z"
    assert letter(26) == "aa"
    assert letter(27) == "ab"
    assert letter(26 + 26) == "ba"
    with pytest.raises(ValueError):
        letter(-1)


def test_alphabet_reuses_assigns_and_resets():
    alphabet = BridgeAlphabet()
    first, second = object(), object()

    assert alphabet.letter_for(first) == "a"
    assert alphabet.letter_for(second) == "b"
    assert alphabet.letter_for(first) == "a"
    assert len(alphabet) == 2 and first in alphabet
    assert alphabet.class_of("b") is second
    with pytest.raises(KeyError):
        alphabet.class_of("z")

    alphabet.reset()
    assert len(alphabet) == 0
    assert alphabet.letter_for(second) == "a"
