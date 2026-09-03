"""Session-level pseudoneighbor facade.

Covers the one-call API that fans compute/plot across a session's fixed
points, mirroring the strong-pip facade tests.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless: exercise the plot helpers without a display
import matplotlib.pyplot as plt
import numpy as np
import pytest

from tanglepack import TangleSession


@pytest.fixture
def henon_session(henon_map, henon_map_inverse):
    """A single-saddle k=10 session with bridges cut — the session analogue of
    the ``henon_tangle_with_bridges`` workbench fixture."""
    session = TangleSession(henon_map, henon_map_inverse)
    fp = session.construct_fixed_point([4, -4])
    session.orient_eigenvectors(
        fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
    )
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=9)
    session.grow_until_turnaround(fp, "stable")
    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    return session, fp


def test_compute_all_fixed_points_returns_dict(henon_session):
    """No-arg compute returns a per-fixed-point dict of reference pairs."""
    session, fp = henon_session

    result = session.compute_pseudoneighbors()

    assert set(result) == {fp}
    assert result[fp], "the tangle should contain reference pairs"
    assert result[fp] == session.trellis(fp).reference_pseudoneighbors


def test_compute_single_fixed_point_returns_list(henon_session):
    """Passing one fixed point returns just its reference list (not a dict)."""
    session, fp = henon_session

    references = session.compute_pseudoneighbors(fp)

    assert isinstance(references, list)
    assert all(p.is_reference for p in references)


def test_plot_helpers_draw_pairs_and_holes(henon_session):
    """plot_pseudoneighbors computes on demand; plot_holes draws punched holes."""
    session, fp = henon_session

    plt.figure()
    try:
        pair_handles = session.plot_pseudoneighbors()
        assert len(pair_handles) == 1

        session.trellis(fp).punch_holes()
        hole_handles = session.plot_holes()
        assert hole_handles
    finally:
        plt.close()


# --------------------------------------------------------------------------- #
# hole / partition fan-outs (plan Phase 5)
# --------------------------------------------------------------------------- #
def test_punch_and_partition_fan_out_to_every_fixed_point(henon_session):
    """Session-level punch/partition fill the per-fixed-point trellises."""
    session, fp = henon_session
    session.compute_pseudoneighbors()

    holes = session.punch_holes()
    assert set(holes) == {fp}
    assert holes[fp] == session.trellis(fp).holes

    results = session.partition_stable_manifold()
    assert set(results) == {fp}
    assert results[fp] == session.trellis(fp).stable_partitions
    assert results[fp], "a partitioned branch yields one result per side"


def test_punch_and_partition_single_fixed_point_return_lists(henon_session):
    """Passing one fixed point returns its own list, not a dict."""
    session, fp = henon_session
    session.compute_pseudoneighbors(fp)

    holes = session.punch_holes(fp)
    assert isinstance(holes, list) and holes

    results = session.partition_stable_manifold(fp)
    assert isinstance(results, list)
    assert [r.side for r in results] == ["left", "right"]


def test_describe_fan_outs_mention_every_fixed_point(henon_session):
    """The describe helpers concatenate the per-trellis reports under a header."""
    session, fp = henon_session
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()

    holes_report = session.describe_holes()
    assert repr(fp) in holes_report
    assert session.trellis(fp).describe_holes() in holes_report

    partition_report = session.describe_stable_partitions()
    assert repr(fp) in partition_report
    assert session.trellis(fp).describe_stable_partitions() in partition_report


def test_plot_stable_partition_fans_out(henon_session):
    """The session plot helper draws one row set per partitioned trellis."""
    session, fp = henon_session
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()

    plt.figure()
    try:
        handles = session.plot_stable_partition()
        assert len(handles) == 1
    finally:
        plt.close()
