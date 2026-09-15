"""TangleSession.dual_graph wrapper.

:class:`~tanglepack.topology.DualGraph.DualGraph` needs a trellis's
arrangement, the stable partitions gathered across every cached
per-fixed-point trellis, and those trellises' strong pips.
:meth:`TangleSession.dual_graph` is the session-level wrapper that supplies
them, mirroring the composition
:meth:`~tanglepack.loom.TangleSession.TangleSession.bridge_classes` already
pins in ``tests/test_session_bridge_classes.py``: this module pins the same
shape one level up, plus the ``plot_*`` delegates that read from it.
"""

from __future__ import annotations

import logging

import matplotlib

matplotlib.use("Agg")  # headless: the session fixtures touch the plotting stack
import matplotlib.pyplot as plt
import pytest

from tanglepack.topology.DualGraph import DualGraph


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _direct_partitions(session, fixed_point_list):
    """The same gathered-partition list session.dual_graph() should use."""
    partitions = []
    for fp in fixed_point_list:
        partitions.extend(session.trellis(fp).stable_partitions)
    return partitions


def _direct_dual_graph(session, fixed_point_list):
    """The same computation session.dual_graph() should produce, done by hand."""
    partitions = _direct_partitions(session, fixed_point_list)
    pips = [session.trellis(fp).strong_pip for fp in fixed_point_list]
    return DualGraph(session.arrangement(), partitions, strong_pips=pips)


# --------------------------------------------------------------------------- #
# (a) dual_graph: the session result matches a direct build
# --------------------------------------------------------------------------- #
def test_k10_dual_graph_matches_a_direct_build(k10_partitioned):
    session, fp = k10_partitioned

    expected = _direct_dual_graph(session, [fp])
    actual = session.dual_graph()

    assert set(actual.arc_nodes) == set(expected.arc_nodes)
    assert len(actual.face_nodes) == len(expected.face_nodes)
    assert actual.fill_segments == expected.fill_segments
    assert actual.arc_nodes, "the fixture must produce at least one arc node"


@pytest.mark.slow
def test_p3_dual_graph_matches_a_direct_build(p3_partitioned):
    session, fp3, fp1 = p3_partitioned

    expected = _direct_dual_graph(session, [fp3, fp1])
    actual = session.dual_graph()

    assert set(actual.arc_nodes) == set(expected.arc_nodes)
    assert len(actual.face_nodes) == len(expected.face_nodes)
    assert actual.fill_segments == expected.fill_segments
    assert actual.arc_nodes, "the fixture must produce at least one arc node"


# --------------------------------------------------------------------------- #
# (b) dual_graph: caching
# --------------------------------------------------------------------------- #
def test_dual_graph_cache_hit_returns_the_same_object(k10_partitioned):
    session, _fp = k10_partitioned

    first = session.dual_graph()

    assert session.dual_graph() is first


def test_dual_graph_rebuild_flag_forces_a_fresh_equivalent_graph(k10_partitioned):
    session, _fp = k10_partitioned

    first = session.dual_graph()
    second = session.dual_graph(rebuild=True)

    assert second is not first
    assert set(second.arc_nodes) == set(first.arc_nodes)
    assert second.fill_segments == first.fill_segments


def test_dual_graph_workbench_mutation_invalidates_the_cache(k10_partitioned):
    session, fp = k10_partitioned
    first = session.dual_graph()
    generation_before = session.workbench.generation

    # Bump the workbench generation and restore a fully partitioned state, so
    # the post-mutation call recomputes successfully rather than raising.
    session.grow_n_times(fp, "unstable", num_iterations=1)
    session.compute_intersections([fp], preserve_ids=True)
    session.create_bridges(fp)
    session.classify_strong_pips()
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()

    assert session.workbench.generation != generation_before
    assert session.dual_graph() is not first


def test_dual_graph_repartitioning_at_the_same_generation_invalidates_the_cache(
    k10_partitioned,
):
    """Re-partitioning is a Trellis-level mutation: the workbench generation
    alone cannot see it, so the (generation, partition signature) cache key
    must."""
    session, fp = k10_partitioned
    first = session.dual_graph()
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

    assert session.workbench.generation == generation_before

    second = session.dual_graph()
    assert second is not first


# --------------------------------------------------------------------------- #
# (d) plotting delegates run headless and return Axes
# --------------------------------------------------------------------------- #
def test_session_plot_dual_graph_returns_axes(k10_partitioned):
    session, _fp = k10_partitioned

    plt.figure()
    try:
        ax = plt.gca()
        result = session.plot_dual_graph(ax=ax)
        assert result is ax
    finally:
        plt.close()


# --------------------------------------------------------------------------- #
# (e) the gathered strong pips drive the fill
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_p3_strong_pips_drive_the_fill_on_both_pip_branches(p3_partitioned, caplog):
    session, fp3, fp1 = p3_partitioned

    pip_branches = set()
    for fp in (fp3, fp1):
        strong_pip = session.trellis(fp).strong_pip
        assert strong_pip is not None
        pip_branches.add(session.trellis(fp).intersection(strong_pip).manifold_b_key)
    assert len(pip_branches) == 2

    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.DualGraph"):
        dg = session.dual_graph(rebuild=True)

    assert not caplog.records, "both pips must reach the constructor: no warning"
    assert set(dg.fill_segments) == pip_branches
    for low, high in dg.fill_segments.values():
        assert high > low
    filled_branches = {node.branch_key for node in dg.arc_nodes.values() if node.filled}
    assert filled_branches == pip_branches


def test_dual_graph_pip_change_at_the_same_generation_and_partition_invalidates_the_cache(
    k10_partitioned,
):
    """A new strong pip changes the fill, so it must be part of the cache key
    even when neither the workbench generation nor the partition moved."""
    session, fp = k10_partitioned
    first = session.dual_graph()
    generation_before = session.workbench.generation
    signature_before = session._partition_signature(session._gathered_partitions())

    trellis = session.trellis(fp)
    alt = trellis.iterate(trellis.strong_pip, 1)
    assert alt in trellis.strong_pip_candidates
    trellis.set_strong_pip(alt)

    assert session.workbench.generation == generation_before
    assert session._partition_signature(session._gathered_partitions()) == signature_before

    second = session.dual_graph()
    assert second is not first
    assert second.fill_segments != first.fill_segments
    assert session.dual_graph() is second
