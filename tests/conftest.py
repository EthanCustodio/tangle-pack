"""Shared fixtures for the tanglepack test suite.

The default map is the ``k=10, b=1`` binary-horseshoe Hénon map, which has a
single clean saddle at ``[4, -4]`` and is used by all the fast invariant tests.
The heavy ``henon_p3_session`` fixture rebuilds the nested period-3 tangle from
``scripts/henon_blast_period_3.py`` and is used only by the blast regression test.
"""

from __future__ import annotations

import numpy as np
import pytest

from tanglepack import TangleWorkbench, TangleSession


# --------------------------------------------------------------------------- #
# Maps
# --------------------------------------------------------------------------- #
def _henon_map(point):
    k, b = 10, 1
    x, y = point
    return np.array([y - k + x**2, -b * x])


def _henon_map_inverse(point):
    k, b = 10, 1
    x, y = point
    return np.array([-y / b, x + k - (y**2) / (b**2)])


def _henon_jacobian(point):
    k, b = 10, 1
    x, y = point
    return np.array([[2 * x, 1], [-b, 0]])


@pytest.fixture
def henon_map():
    return _henon_map


@pytest.fixture
def henon_map_inverse():
    return _henon_map_inverse


@pytest.fixture
def henon_jacobian():
    return _henon_jacobian


# --------------------------------------------------------------------------- #
# Workbench build-up (k=10 single saddle at [4, -4])
# --------------------------------------------------------------------------- #
@pytest.fixture
def workbench(henon_map, henon_map_inverse):
    """A fresh workbench with no fixed points yet."""
    return TangleWorkbench(henon_map, henon_map_inverse)


@pytest.fixture
def fixed_point(workbench):
    """Return ``(workbench, fp)`` with the saddle constructed and oriented."""
    fp = workbench.construct_fixed_point([4, -4])
    workbench.orient_eigenvectors(
        fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
    )
    return workbench, fp


@pytest.fixture
def initialized(fixed_point):
    """``(workbench, fp)`` with both fundamental segments initialized."""
    workbench, fp = fixed_point
    workbench.initialize_both_manifolds(fp)
    return workbench, fp


def _unstable_manifold(workbench, fp, branch_index: int = 0):
    return workbench.manifolds[(fp, "unstable", 0, branch_index)]


def _stable_manifold(workbench, fp, branch_index: int = 0):
    return workbench.manifolds[(fp, "stable", 0, branch_index)]


@pytest.fixture
def grown_unstable(initialized):
    """``(workbench, fp, manifold)`` with the unstable manifold grown 5 times.

    Capped at 5: with the k=10 binary-horseshoe map, refinement currently
    explodes the point count beyond ~6 iterations (a known bug pinned by the
    regression suite), so the fast fixtures stay at 5 to remain snappy.
    """
    workbench, fp = initialized
    workbench.grow_n_times(fp, "unstable", num_iterations=5)
    return workbench, fp, _unstable_manifold(workbench, fp)


@pytest.fixture
def grown_both(initialized):
    """``(workbench, fp)`` with unstable grown 5x and stable grown to turnaround."""
    workbench, fp = initialized
    workbench.grow_n_times(fp, "unstable", num_iterations=5)
    workbench.grow_until_turnaround(fp, "stable")
    return workbench, fp


@pytest.fixture
def small_tangle(grown_both):
    """``(workbench, fp)`` after intersections are computed."""
    workbench, fp = grown_both
    workbench.compute_intersections([fp])
    return workbench, fp


# --------------------------------------------------------------------------- #
# Heavy nested period-3 session (for the blast regression test)
# --------------------------------------------------------------------------- #
def _p3_map(point):
    k, b = 2, 1
    x, y = point
    return np.array([y - k + x**2, -b * x])


def _p3_map_inverse(point):
    k, b = 2, 1
    x, y = point
    return np.array([-y / b, x + k - (y**2) / (b**2)])


def _p3_jacobian(point):
    k, b = 2, 1
    x, y = point
    return np.array([[2 * x, 1], [-b, 0]])


@pytest.fixture
def henon_p3_session():
    """Rebuild the nested period-3 tangle and its resonance zones.

    Mirrors ``scripts/henon_blast_period_3.py``. Returns ``(session, fp3, fp1,
    inner_zone)``. Function-scoped: blasting mutates the session (it iterates
    bridges), so each blast test needs its own clean build.
    """
    session = TangleSession(_p3_map, _p3_map_inverse, _p3_jacobian)

    session.workbench._man_machine.area_cutoff = 1e-7
    fp3 = session.construct_fixed_point([[0, 1], [-1, 0], [-1, 1]])
    session.orient_eigenvectors(
        fp3, {"unstable": np.array([0, -1]), "stable": np.array([-1, -1])}
    )
    session.initialize_both_manifolds(fp3)
    session.grow_n_times(fp3, "unstable", num_iterations=10)
    session.grow_n_times(fp3, "stable", num_iterations=6)

    session.workbench._man_machine.area_cutoff = 1e-4
    fp1 = session.construct_fixed_point([4, -4])
    session.orient_eigenvectors(
        fp1, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
    )
    session.initialize_both_manifolds(fp1)
    session.grow_n_times(fp1, "unstable", num_iterations=7)
    session.grow_until_turnaround(fp1, "stable")

    session.compute_intersections([fp3, fp1])
    session.trim_stable_manifolds(fp3)
    session.trim_stable_manifolds(fp1)
    session.create_bridges(fp3)
    session.create_bridges(fp1)
    session.infer_iterate_table()

    T1 = session.trellis(fp1)
    T1.classify_strong_pips()
    T3 = session.trellis(fp3)
    T3.classify_strong_pips()
    # Use the default strong pip each trellis chooses (smallest unstable cdist).
    # The intersection ids are not stable across refinement changes, so do not
    # hard-code a specific id here.
    session.add_resonance_zones([T1.strong_pip, T3.strong_pip])

    inner_zone = max(session.resonance_zones.values(), key=lambda z: z.area)
    return session, fp3, fp1, inner_zone
