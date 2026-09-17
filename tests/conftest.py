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
from tanglepack.examples import (
    HENON_K10,
    HENON_P3,
    henon_jacobian as _henon_jacobian_factory,
    henon_map as _henon_map_factory,
    henon_map_inverse as _henon_map_inverse_factory,
    saddle_guesses,
)


# --------------------------------------------------------------------------- #
# Maps
# --------------------------------------------------------------------------- #
# Batch-capable (coordinate on axis 0): a single (2,) point or a (2, N) batch.
_henon_map = _henon_map_factory(*HENON_K10)
_henon_map_inverse = _henon_map_inverse_factory(*HENON_K10)
_henon_jacobian = _henon_jacobian_factory(*HENON_K10)


@pytest.fixture(name="henon_map")
def _k10_map_fixture():
    return _henon_map


@pytest.fixture(name="henon_map_inverse")
def _k10_map_inverse_fixture():
    return _henon_map_inverse


@pytest.fixture(name="henon_jacobian")
def _k10_jacobian_fixture():
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
    fp = workbench.construct_fixed_point(saddle_guesses(*HENON_K10)["saddle"])
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
    """``(workbench, fp, manifold)`` with the unstable manifold grown 7 times.

    7 iterations of the k=10 binary-horseshoe map develops a small tangle that
    actually crosses the stable manifold (a handful of intersections/bridges)
    while staying snappy. The old "explodes beyond ~6 iterations" cap was an
    artifact of a mislocated orbit (the multipoint solver used to root-find
    element-wise); with the orbit located properly the count grows gradually.
    """
    workbench, fp = initialized
    workbench.grow_n_times(fp, "unstable", num_iterations=7)
    return workbench, fp, _unstable_manifold(workbench, fp)


@pytest.fixture
def grown_both(initialized):
    """``(workbench, fp)`` with unstable grown 7x and stable grown to turnaround."""
    workbench, fp = initialized
    workbench.grow_n_times(fp, "unstable", num_iterations=7)
    workbench.grow_until_turnaround(fp, "stable")
    return workbench, fp


@pytest.fixture
def small_tangle(grown_both):
    """``(workbench, fp)`` after intersections are computed."""
    workbench, fp = grown_both
    workbench.compute_intersections([fp])
    return workbench, fp


@pytest.fixture
def henon_tangle_with_bridges(grown_both):
    """``(workbench, fp)`` grown two extra unstable iterations (the 7 of
    ``grown_both`` leave too few crossings for a full pseudoneighbor reference
    window), with intersections, trimmed stable manifold, and bridges."""
    workbench, fp = grown_both
    workbench.grow_n_times(fp, "unstable", num_iterations=2)
    workbench.compute_intersections([fp])
    workbench.trim_stable_manifolds(fp)
    workbench.create_bridges(fp)
    return workbench, fp


# --------------------------------------------------------------------------- #
# Inversion saddle (k_value == 2 * period) -- same k=10, b=1 map
# --------------------------------------------------------------------------- #
# The map's OTHER fixed point, at (-2.3166, 2.3166), has BOTH eigenvalues
# negative (-4.406 and -0.227, product 1): an orientation-preserving saddle with
# inversion, so k_value = 2 * period and the two manifold branches are one chain.
#
# It is deliberately not a b < 0 Hénon map: there det J = b < 0, so exactly ONE
# eigenvalue is negative at every fixed point and the unstable and stable
# manifolds have different inversion status -- something a single k_value cannot
# represent (``FixedPoint.set_k_value`` rejects it).
#
# Eigenvector orientation is not applicable here: ``orient_eigenvectors`` is a
# documented no-op on an inversion point, because the second branch is the image
# of the first and so its direction is fixed by the map, not by the caller.
_INVERSION_GUESS = saddle_guesses(*HENON_K10)["inversion"]


def _build_inversion_workbench() -> tuple[TangleWorkbench, object]:
    """A workbench with the inversion saddle constructed and both manifolds seeded."""
    workbench = TangleWorkbench(_henon_map, _henon_map_inverse, _henon_jacobian)
    fp = workbench.construct_fixed_point(_INVERSION_GUESS)
    assert fp.check_inversion(), "the inversion fixture must have inversion"
    assert fp.k_value == 2 * fp.period
    workbench.initialize_both_manifolds(fp)
    return workbench, fp


@pytest.fixture
def henon_inversion_initialized():
    """``(workbench, fp)`` with the inversion saddle's fundamental segments only."""
    return _build_inversion_workbench()


@pytest.fixture(scope="module")
def henon_inversion():
    """``(workbench, fp)`` with the inversion saddle grown on both branches and
    its crossings computed.

    Module-scoped and read-only: growing the second branch of this saddle is the
    expensive part (its pieces refine hard), and the tests that use it only
    inspect the result. Anything that mutates the workbench must build its own.
    Computing the crossings leaves the manifold geometry untouched, so the growth
    invariants are still the growth invariants.
    """
    workbench, fp = _build_inversion_workbench()
    workbench.grow_n_times(fp, "unstable", num_iterations=6)
    workbench.grow_n_times(fp, "stable", num_iterations=5)
    workbench.compute_intersections([fp])
    return workbench, fp


# --------------------------------------------------------------------------- #
# Heavy nested period-3 session (for the blast regression test)
# --------------------------------------------------------------------------- #
_p3_map = _henon_map_factory(*HENON_P3)
_p3_map_inverse = _henon_map_inverse_factory(*HENON_P3)
_p3_jacobian = _henon_jacobian_factory(*HENON_P3)


@pytest.fixture
def henon_p3_session():
    """Rebuild the nested period-3 tangle and its resonance zones.

    Mirrors ``scripts/henon_blast_period_3.py``. Returns ``(session, fp3, fp1,
    inner_zone)``. Function-scoped: blasting mutates the session (it iterates
    bridges), so each blast test needs its own clean build.
    """
    session = TangleSession(_p3_map, _p3_map_inverse, _p3_jacobian)

    session.workbench._man_machine.area_cutoff = 1e-7
    fp3 = session.construct_fixed_point(saddle_guesses(*HENON_P3)["period_3"])
    session.orient_eigenvectors(
        fp3, {"unstable": np.array([0, -1]), "stable": np.array([-1, -1])}
    )
    session.initialize_both_manifolds(fp3)
    session.grow_n_times(fp3, "unstable", num_iterations=13)
    session.grow_n_times(fp3, "stable", num_iterations=9)

    session.workbench._man_machine.area_cutoff = 1e-4
    fp1 = session.construct_fixed_point(saddle_guesses(*HENON_P3)["saddle"])
    session.orient_eigenvectors(
        fp1, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
    )
    session.initialize_both_manifolds(fp1)
    session.grow_n_times(fp1, "unstable", num_iterations=11)
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


# --------------------------------------------------------------------------- #
# Session-level fixtures shared by the region and bridge-class tests
# --------------------------------------------------------------------------- #
@pytest.fixture
def k10_session(henon_map, henon_map_inverse):
    """A k=10 session grown far enough to close several faces."""
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
    session.infer_iterate_table()
    return session, fp


@pytest.fixture
def k10_partitioned(k10_session):
    """``(session, fp)`` with the k=10 trellis classified, punched and partitioned.

    Runs the whole flow through the session fan-outs, so the single saddle's
    stable branch carries holes and both of its sides' partitions.
    """
    session, fp = k10_session
    session.classify_strong_pips()
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()
    assert session.trellis(fp).stable_partitions
    return session, fp


@pytest.fixture
def p3_partitioned(henon_p3_session):
    """``(session, fp3, fp1)`` with every trellis classified, punched, partitioned.

    Runs the whole flow through the session fan-outs, so both the period-3 and the
    period-1 tangle of the nested fixture carry holes and partitions.
    """
    session, fp3, fp1, _zone = henon_p3_session
    session.classify_strong_pips()
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()
    assert session.trellis(fp3).stable_partitions
    return session, fp3, fp1


# --------------------------------------------------------------------------- #
# Blasted k=2.8 session (for the backward-only holes / inertness tests)
# --------------------------------------------------------------------------- #
HENON_K28 = (2.8, 1)
_k28_map = _henon_map_factory(*HENON_K28)
_k28_map_inverse = _henon_map_inverse_factory(*HENON_K28)
_k28_jacobian = _henon_jacobian_factory(*HENON_K28)


@pytest.fixture(scope="session")
def k28_partitioned():
    """``(session, fp)``: the k=2.8 case of ``scripts/henon_bridge_classes.py``.

    Ten unstable steps at ``area_cutoff = 1e-7``, the stable manifold trimmed at
    the pip f(q0) (a resonance zone), ONE single-iteration blast of that zone,
    then classify / pseudoneighbors / holes / partition. The blast registers the
    forward images of one reference pair as a blast-child bridge, which is what
    makes this the fixture for "holes propagate backward only": before that
    rule a direct hole was punched in the blast child. Session-scoped (about a
    minute to build); the tests that use it must not mutate the session.
    Registry ids are not reproducible between builds — locate crossings by
    stable-cdist order, never by id.
    """
    session = TangleSession(_k28_map, _k28_map_inverse, _k28_jacobian)
    session.workbench._man_machine.area_cutoff = 1e-7
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
    assert pip is not None, "the default strong pip must have a registered image"

    zone = session.resonance_zone(pip)
    session.classify_strong_pips()
    session.set_strong_pip(fp, pip)
    session.blast_zone(zone, num_iterations=1, fixed_point=[fp], min_separation=1e-5)
    session.classify_strong_pips()
    session.set_strong_pip(fp, pip)

    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()
    assert session.trellis(fp).stable_partitions
    return session, fp
