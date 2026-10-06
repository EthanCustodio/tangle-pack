"""Session-level pseudoneighbor, hole and partition plot helpers.

The compute / punch / partition / describe fan-out shapes are
``tests/facade/test_session_fanouts.py``.
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
