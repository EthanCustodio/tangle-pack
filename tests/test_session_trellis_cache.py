"""Phase 4: the Trellis cache is keyed by the workbench generation.

A Trellis is a snapshot: it buckets the registry's crossings per branch and
holds the bridge list of the moment. :meth:`TangleSession.trellis` may therefore
hand back a cached object only while nothing has changed underneath it. Since
Phase 4 that is a single comparison — the generation the snapshot was built at
against the workbench's current one — so this module pins the hit/miss table
rather than the individual staleness symptoms it replaced.
"""

from __future__ import annotations

import numpy as np
import pytest

from tanglepack import TangleSession


# --------------------------------------------------------------------------- #
# fixture: the cheapest session with crossings, bridges and a zone
# --------------------------------------------------------------------------- #
def _k10_map(point):
    k, b = 10, 1
    x, y = point
    return np.stack([y - k + x**2, -b * x], axis=0)


def _k10_map_inverse(point):
    k, b = 10, 1
    x, y = point
    return np.stack([-y / b, x + k - (y**2) / (b**2)], axis=0)


def _k10_jacobian(point):
    k, b = 10, 1
    x, y = point
    return np.array([[2 * x, 1], [-b, 0]])


@pytest.fixture
def k10_session():
    """``(session, fp)`` — the k=10 saddle with intersections and bridges built."""
    session = TangleSession(_k10_map, _k10_map_inverse, _k10_jacobian)
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


def _define_zone(session, fp):
    """Trim at a strong pip strictly inside the stable extent (see
    tests/test_loom_blast_restore.py, same reasoning)."""
    registry = session.workbench.intersection_registry
    trellis = session.trellis(fp)
    trellis.classify_strong_pips()
    outermost = max(registry[i].stable_cdist for i in registry.all_ids())
    inner = [
        c
        for c in trellis.strong_pip_candidates
        if registry[c].stable_cdist < outermost
    ]
    assert inner, "expected a strong-pip candidate inside the stable extent"
    pip = max(inner, key=lambda c: registry[c].stable_cdist)
    trellis.set_strong_pip(pip)
    session.add_resonance_zones([pip])
    return session.resonance_zones[(fp, 0)]


# --------------------------------------------------------------------------- #
# cache hits
# --------------------------------------------------------------------------- #
def test_same_generation_is_a_cache_hit(k10_session):
    session, fp = k10_session

    first = session.trellis(fp)

    assert session.trellis(fp) is first
    assert first._built_generation == session.workbench.generation


def test_reads_do_not_invalidate_the_cache(k10_session):
    session, fp = k10_session
    first = session.trellis(fp)

    session.workbench.bridges
    session.workbench.intersection_registry.by_unstable_cdist
    first.classify_strong_pips()

    assert session.trellis(fp) is first


def test_rebuild_flag_forces_a_miss(k10_session):
    session, fp = k10_session
    first = session.trellis(fp)

    assert session.trellis(fp, rebuild=True) is not first


# --------------------------------------------------------------------------- #
# cache misses, one per mutation path
# --------------------------------------------------------------------------- #
def test_growth_invalidates_the_cache(k10_session):
    session, fp = k10_session
    first = session.trellis(fp)

    session.grow_n_times(fp, "unstable", num_iterations=1)

    assert session.trellis(fp) is not first


def test_compute_intersections_invalidates_the_cache(k10_session):
    session, fp = k10_session
    first = session.trellis(fp)

    session.compute_intersections([fp])

    fresh = session.trellis(fp)
    assert fresh is not first
    assert fresh.registry is session.workbench.intersection_registry


def test_iterate_bridge_invalidates_the_cache(k10_session):
    session, fp = k10_session
    first = session.trellis(fp)

    session.workbench.iterate_bridge(session.workbench.uniiterated_bridges[0])

    assert session.trellis(fp) is not first


def test_rebuild_bridges_invalidates_the_cache(k10_session):
    session, fp = k10_session
    first = session.trellis(fp)

    session.workbench.rebuild_bridges(fp)

    assert session.trellis(fp) is not first


def test_add_resonance_zones_invalidates_the_cache(k10_session):
    session, fp = k10_session
    _define_zone(session, fp)
    first = session.trellis(fp)

    session.add_resonance_zones([first.registry.all_ids()[0]])

    assert session.trellis(fp) is not first


def test_restore_invalidates_the_cache(k10_session):
    """The gap the manual ``invalidate_trellises`` calls used to leave open."""
    session, fp = k10_session
    zone = _define_zone(session, fp)
    cached = session.trellis(fp)

    zone.restore(session.workbench)

    assert session.trellis(fp) is not cached


def test_invalidate_trellises_is_a_deprecated_no_op(k10_session):
    """Kept as an alias for one release; it warns and drops nothing."""
    session, fp = k10_session
    first = session.trellis(fp)

    with pytest.warns(DeprecationWarning, match="no-op"):
        session.invalidate_trellises()

    assert session.trellis(fp) is first
