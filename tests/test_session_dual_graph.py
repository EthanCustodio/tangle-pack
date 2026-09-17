"""
TangleSession.dual_graph / minimal_trellis / iterated_partition: results,
caching and invalidation.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from minimal_helpers import build_pieces
from tanglepack.topology import plotting
from tanglepack.topology.DualGraph import DualGraph
from tanglepack.topology.PartitionFamily import HomotopyPartition


def _same_graph(a: DualGraph, b: DualGraph) -> bool:
    return (
        set(a.stable_nodes) == set(b.stable_nodes)
        and [f.corners for f in a.face_nodes] == [f.corners for f in b.face_nodes]
        and a.fundamental_segments == b.fundamental_segments
        and {k: n.elements for k, n in a.stable_nodes.items()}
        == {k: n.elements for k, n in b.stable_nodes.items()}
    )


def test_k10_dual_graph_matches_a_direct_build(k10_partitioned):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    direct = DualGraph(pieces.minimal, pieces.iterated, strong_pips=pieces.strong_pips)
    assert _same_graph(session.dual_graph(), direct)


@pytest.mark.slow
def test_p3_dual_graph_matches_a_direct_build(p3_partitioned):
    session, fp3, fp1 = p3_partitioned
    pieces = build_pieces(session, [fp3, fp1])
    direct = DualGraph(pieces.minimal, pieces.iterated, strong_pips=pieces.strong_pips)
    assert _same_graph(session.dual_graph(), direct)


def test_dual_graph_is_built_over_the_cached_pieces(k10_partitioned):
    session, fp = k10_partitioned
    dual = session.dual_graph()
    assert dual.minimal is session.minimal_trellis()
    assert dual.partition is session.iterated_partition()
    assert session.iterated_partition().minimal is session.minimal_trellis()
    homotopy = session.homotopy_partition()
    assert isinstance(homotopy, HomotopyPartition)
    assert homotopy.signature() == session._partition_signature(
        session._gathered_partitions()
    )
    assert session.iterated_partition().homotopy.signature() == homotopy.signature()


def test_cache_hits_return_the_same_objects(k10_partitioned):
    session, fp = k10_partitioned
    assert session.dual_graph() is session.dual_graph()
    assert session.minimal_trellis() is session.minimal_trellis()
    assert session.iterated_partition() is session.iterated_partition()
    assert session.dual_graph(fp) is session.dual_graph(fp)


def test_rebuild_flag_forces_fresh_equivalent_objects(k10_partitioned):
    session, fp = k10_partitioned
    before = session.dual_graph()
    minimal_before = session.minimal_trellis()
    after = session.dual_graph(rebuild=True)
    assert after is not before and _same_graph(before, after)
    assert session.minimal_trellis() is not minimal_before
    assert after.minimal is session.minimal_trellis()


def test_workbench_mutation_invalidates_the_caches(k10_partitioned):
    session, fp = k10_partitioned
    before = session.dual_graph()
    minimal_before = session.minimal_trellis()
    iterated_before = session.iterated_partition()
    # Bump the workbench generation and restore a fully partitioned state, so
    # the post-mutation call recomputes successfully rather than raising.
    session.grow_n_times(fp, "unstable", num_iterations=1)
    session.compute_intersections([fp], preserve_ids=True)
    session.create_bridges(fp)
    session.classify_strong_pips()
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()
    assert session.minimal_trellis() is not minimal_before
    assert session.iterated_partition() is not iterated_before
    assert session.dual_graph() is not before


def test_repartitioning_at_the_same_generation_invalidates_the_caches(k10_partitioned):
    session, fp = k10_partitioned
    before = session.dual_graph()
    generation = session.workbench.generation
    trellis = session.trellis(fp)
    trellis.clear_results()
    trellis.classify_strong_pips()
    alternates = [c for c in trellis.strong_pip_candidates if c != trellis.strong_pip]
    assert alternates, "the k10 fixture must offer more than one strong-pip candidate"
    trellis.set_strong_pip(alternates[0])
    session.compute_pseudoneighbors(fp)
    session.punch_holes(fp)
    session.partition_stable_manifold(fp)
    assert session.workbench.generation == generation
    after = session.dual_graph()
    assert after is not before
    assert session.iterated_partition().homotopy.signature() == session._partition_signature(
        session._gathered_partitions()
    )


def test_pip_change_at_the_same_generation_and_partition_invalidates_the_graph(
    k10_partitioned,
):
    session, fp = k10_partitioned
    trellis = session.trellis(fp)
    alternatives = [c for c in trellis.strong_pip_candidates if c != trellis.strong_pip]
    if not alternatives:
        pytest.skip("the fixture has a single strong-pip candidate")
    before = session.dual_graph()
    minimal_before = session.minimal_trellis()
    trellis.set_strong_pip(alternatives[0])
    after = session.dual_graph()
    assert after is not before
    assert {n.edge_key for n in after.unified_nodes} != {
        n.edge_key for n in before.unified_nodes
    }
    # The minimal trellis and the partition do not depend on the pip choice.
    assert session.minimal_trellis() is minimal_before


def test_describe_iterated_partition(k10_partitioned):
    session, fp = k10_partitioned
    assert session.describe_iterated_partition() == session.iterated_partition().describe()


def test_session_plot_dual_graph_returns_axes(k10_partitioned):
    session, fp = k10_partitioned
    fig, ax = plt.subplots()
    try:
        assert session.plot_dual_graph(ax=ax) is ax
        assert session.plot_minimal_trellis(ax=ax) is ax
    finally:
        plt.close(fig)


def test_session_plot_delegates_to_plotting(k10_partitioned, monkeypatch):
    session, fp = k10_partitioned
    calls = []
    monkeypatch.setattr(
        plotting, "plot_dual_graph",
        lambda graph, ax=None, **kw: calls.append(("dual", graph, ax, kw)) or "dual",
    )
    monkeypatch.setattr(
        plotting, "plot_minimal_trellis",
        lambda minimal, ax=None, **kw: calls.append(("minimal", minimal, ax, kw)) or "minimal",
    )
    assert session.plot_dual_graph(show_labels=True) == "dual"
    assert session.plot_minimal_trellis(show_dropped=False) == "minimal"
    assert calls[0][1] is session.dual_graph() and calls[0][3] == {"show_labels": True}
    assert calls[1][1] is session.minimal_trellis() and calls[1][3] == {"show_dropped": False}
