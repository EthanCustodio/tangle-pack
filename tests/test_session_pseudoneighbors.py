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
