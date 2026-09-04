"""D.6 -- TangleSession.dual_graph / symbolic_dynamics wrappers.

:class:`~tanglepack.topology.DualGraph.DualGraph` (Phase B) and
:func:`~tanglepack.topology.SymbolicDynamics.symbolic_dynamics` (Phase D) both
need a trellis's arrangement/bridge classes and the stable partitions gathered
across every cached per-fixed-point trellis. :meth:`TangleSession.dual_graph`
and :meth:`TangleSession.symbolic_dynamics` are the session-level wrappers that
supply those, mirroring the composition
:meth:`~tanglepack.loom.TangleSession.TangleSession.bridge_classes` already
pins in ``tests/test_session_bridge_classes.py``: this module pins the same
shape one (and two) levels up, plus the two ``plot_*`` delegates that read from
them.
"""

from __future__ import annotations

import logging

import matplotlib

matplotlib.use("Agg")  # headless: the session fixtures touch the plotting stack
import matplotlib.pyplot as plt
import pytest

from tanglepack.topology.DualGraph import DualGraph
from tanglepack.topology.SymbolicDynamics import (
    symbolic_dynamics as compute_symbolic_dynamics,
)


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
    return DualGraph(session.arrangement(), partitions)


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
# (c) symbolic_dynamics: the session result matches a direct build
# --------------------------------------------------------------------------- #
def test_k10_symbolic_dynamics_matches_a_direct_build(k10_partitioned):
    session, fp = k10_partitioned

    dg = _direct_dual_graph(session, [fp])
    classes = session.bridge_classes()
    expected = compute_symbolic_dynamics(dg, classes)

    actual = session.symbolic_dynamics()

    assert actual.describe() == expected.describe()


@pytest.mark.slow
def test_p3_symbolic_dynamics_matches_a_direct_build(p3_partitioned):
    session, fp3, fp1 = p3_partitioned

    dg = _direct_dual_graph(session, [fp3, fp1])
    classes = session.bridge_classes()
    expected = compute_symbolic_dynamics(dg, classes)

    actual = session.symbolic_dynamics()

    assert actual.describe() == expected.describe()


# --------------------------------------------------------------------------- #
# (b) symbolic_dynamics: caching
# --------------------------------------------------------------------------- #
def test_symbolic_dynamics_cache_hit_returns_the_same_object(k10_partitioned):
    session, _fp = k10_partitioned

    first = session.symbolic_dynamics()

    assert session.symbolic_dynamics() is first


def test_symbolic_dynamics_rebuild_flag_forces_a_fresh_equivalent_result(
    k10_partitioned,
):
    session, _fp = k10_partitioned

    first = session.symbolic_dynamics()
    second = session.symbolic_dynamics(rebuild=True)

    assert second is not first
    assert second.describe() == first.describe()


def test_symbolic_dynamics_workbench_mutation_invalidates_the_cache(k10_partitioned):
    session, fp = k10_partitioned
    first = session.symbolic_dynamics()
    generation_before = session.workbench.generation

    # A no-op recompute at the SAME growth extent still bumps the workbench
    # generation (a fresh registry swap-in): enough to invalidate the cache
    # without disturbing the bridge structure a walk needs to succeed (see
    # test_dual_graph_workbench_mutation_invalidates_the_cache for a mutation
    # that actually grows the manifold; symbolic_dynamics additionally walks
    # every bridge, which is pickier about mid-flight structural changes).
    session.compute_intersections([fp], preserve_ids=True)
    session.create_bridges(fp)
    session.infer_iterate_table()
    session.classify_strong_pips()
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()

    assert session.workbench.generation != generation_before
    assert session.symbolic_dynamics() is not first


def test_symbolic_dynamics_repartitioning_at_the_same_generation_invalidates_the_cache(
    k10_partitioned,
):
    session, fp = k10_partitioned
    first = session.symbolic_dynamics()
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

    second = session.symbolic_dynamics()
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


def test_session_plot_transition_graph_returns_axes(k10_partitioned):
    session, _fp = k10_partitioned

    plt.figure()
    try:
        ax = plt.gca()
        result = session.plot_transition_graph(ax=ax)
        assert result is ax
    finally:
        plt.close()


def test_session_plot_transition_graph_accepts_an_explicit_dynamics(k10_partitioned):
    session, _fp = k10_partitioned
    sd = session.symbolic_dynamics()

    plt.figure()
    try:
        ax = plt.gca()
        result = session.plot_transition_graph(sd, ax=ax)
        assert result is ax
    finally:
        plt.close()


# --------------------------------------------------------------------------- #
# (e) on p3, the strong pips are passed through to the dual graph
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_p3_strong_pips_are_passed_through_to_the_dual_graph(p3_partitioned, caplog):
    session, fp3, fp1 = p3_partitioned

    assert session.trellis(fp3).strong_pip is not None
    assert session.trellis(fp1).strong_pip is not None

    with caplog.at_level(
        logging.INFO, logger="tanglepack.topology.DualGraph"
    ):
        dg = session.dual_graph(rebuild=True)

    assert caplog.text.count("as expected") == 2, (
        "both fixed points' strong pips must reach the DualGraph constructor's "
        "fill check, one 'as expected' INFO line per pip branch"
    )

    for fp in (fp3, fp1):
        strong_pip = session.trellis(fp).strong_pip
        branch_key = session.trellis(fp).intersection(strong_pip).manifold_b_key
        lo, hi = dg.fill_segments[branch_key]
        assert hi > lo, f"fixed point {fp!r}'s strong-pip branch has no fill span"
