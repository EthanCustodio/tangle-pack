"""A.4 -- TangleSession.bridge_classes: gathering, caching, warnings.

:func:`tanglepack.topology.BridgeClass.bridge_classes` (A.3) needs a trellis to
read bridges from and a list of stable partitions to resolve elements against.
:meth:`TangleSession.bridge_classes` is the session-level wrapper that supplies
both: the trellis is the selection's own (default: the all-fixed-points one,
mirroring :meth:`TangleSession.arrangement`), and the partitions are gathered
from every cached, non-stale per-fixed-point trellis via the private
:meth:`TangleSession._gathered_partitions` helper. This module pins that
composition, the ``workbench.generation``-keyed cache, the warning logged for
an unpartitioned stable branch, and (on the nested fixture) that no class mixes
the two fixed points.
"""

from __future__ import annotations

import logging

import matplotlib

matplotlib.use("Agg")  # headless: the session fixtures touch the plotting stack
import pytest

from tanglepack.topology.BridgeClass import bridge_classes as compute_bridge_classes


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _direct_classes(session, fixed_point_list):
    """The same computation session.bridge_classes() should produce, done by hand."""
    partitions = []
    for fp in fixed_point_list:
        partitions.extend(session.trellis(fp).stable_partitions)
    return compute_bridge_classes(session.trellis(), partitions)


# --------------------------------------------------------------------------- #
# (a) the session result matches a direct call over the gathered partitions
# --------------------------------------------------------------------------- #
def test_k10_session_result_matches_a_direct_call(k10_partitioned):
    session, fp = k10_partitioned

    expected = _direct_classes(session, [fp])
    actual = session.bridge_classes()

    assert actual == expected
    assert actual, "the fixture must produce at least one class"


@pytest.mark.slow
def test_p3_session_result_matches_a_direct_call(p3_partitioned):
    session, fp3, fp1 = p3_partitioned

    expected = _direct_classes(session, [fp3, fp1])
    actual = session.bridge_classes()

    assert actual == expected
    assert actual, "the fixture must produce at least one class"


def test_session_result_matches_gathered_partitions_helper(k10_partitioned):
    """The private gathering helper session.bridge_classes() actually uses."""
    session, fp = k10_partitioned
    session.bridge_classes()  # populate/observe the same cache the helper reads

    gathered = session._gathered_partitions()
    expected = compute_bridge_classes(session.trellis(), gathered)

    assert session.bridge_classes() == expected


# --------------------------------------------------------------------------- #
# (b) caching
# --------------------------------------------------------------------------- #
def test_cache_hit_returns_the_same_object(k10_partitioned):
    session, _fp = k10_partitioned

    first = session.bridge_classes()

    assert session.bridge_classes() is first


def test_rebuild_flag_forces_a_fresh_equal_dict(k10_partitioned):
    session, _fp = k10_partitioned

    first = session.bridge_classes()
    second = session.bridge_classes(rebuild=True)

    assert second is not first
    assert second == first


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
    and check the cache does not serve the stale grouping."""
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


# --------------------------------------------------------------------------- #
# (c) no partitions: ValueError from A.3, and the warning
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
# (d) on p3, no class mixes fixed points
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_p3_no_class_mixes_fixed_points(p3_partitioned):
    session, fp3, fp1 = p3_partitioned

    classes = session.bridge_classes()

    seen = set()
    for cls in classes:
        owners = {id(cls.first.fixed_point), id(cls.second.fixed_point)}
        assert len(owners) == 1, f"class {cls} mixes two fixed points"
        seen |= owners
    assert seen == {id(fp3), id(fp1)}, "both tangles must contribute classes"
