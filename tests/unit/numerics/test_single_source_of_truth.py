"""Phase 2 Task A — single sources of truth in the numerics core.

These tests pin the registry-only invariants: resolved crossings live only in
the ``IntersectionRegistry``; a bridge's endpoints are the registry ids of the
two crossings it was cut at; every manifold and every bridge carries its
``manifold_key``; and the forward iterate of a bridge's two endpoint crossings
is registered explicitly rather than inferred.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

from helpers.logs import assert_logged

from tanglepack.numerics.BaseManifold import BaseManifold
from tanglepack.numerics.Bridge import Bridge
from helpers.polylines import TEST_STORE, store_rows
from tanglepack.numerics.ManifoldMachine import ManifoldMachine


# --------------------------------------------------------------------------- #
# 2.1 -- resolved crossings live in the registry only
# --------------------------------------------------------------------------- #
def test_compute_intersections_returns_the_registry_it_installs(grown_both, caplog):
    """The return value is the workbench's NEW registry; the old one is stale.

    A reference taken before the call keeps the old object (that is the trap the
    return value and the DEBUG record exist for), so the test holds one across
    the recompute and checks it is not what the workbench now serves.
    """
    workbench, fp = grown_both
    caplog.set_level(logging.DEBUG, logger="tanglepack.numerics.TangleWorkbench")
    first = workbench.compute_intersections([fp], infer_iterates=False)
    assert first is workbench.intersection_registry
    assert len(first) > 0

    second = workbench.compute_intersections([fp], infer_iterates=False)

    assert second is workbench.intersection_registry
    assert second is not first
    # one record: replacing the empty registry of the first call is not logged
    assert_logged(caplog, "DEBUG", "tanglepack.numerics.TangleWorkbench", count=1)


def test_trim_stable_manifolds_reads_the_registry(grown_both):
    """The new tail is the first node at or past the outermost crossing.

    The stable manifold grown to turnaround runs well past its outermost
    crossing, so the trim must actually remove points.
    """
    workbench, fp = grown_both
    workbench.compute_intersections([fp], infer_iterates=False)
    manifold = workbench.manifolds[(fp, "stable", 0, 0)]
    points_before = len(manifold.get_point_array())
    tail_cdist_before = manifold.tail.cdist

    outermost = max(ix.stable_cdist for _iid, ix in workbench.intersection_registry)
    workbench.trim_stable_manifolds(fp)

    tail = manifold.tail
    assert tail.cdist >= outermost
    assert tail.cdist < tail_cdist_before
    assert len(manifold.get_point_array()) < points_before
    # the FIRST such node: the node before it (root to tail) is short of it
    rows = manifold.ids
    assert manifold.store.cdist[rows[-1]] == tail.cdist
    assert manifold.store.cdist[rows[-2]] < outermost


# --------------------------------------------------------------------------- #
# 2.2 -- bridge endpoints set at cut time
# --------------------------------------------------------------------------- #
def test_partial_bridge_reports_its_missing_endpoint(henon_tangle_with_bridges):
    """A piece of an image bounded by fewer than two crossings is ``partial``."""
    workbench, fp = henon_tangle_with_bridges

    partials = []
    for parent in list(workbench.uniiterated_bridges):
        if parent.iterated:
            continue
        children = workbench.iterate_bridge(parent)
        partials.extend(c for c in children if c.partial)
        if partials:
            break

    assert partials, "expected at least one partial bridge on the k=10 blast frontier"
    for bridge in partials:
        assert bridge.first_intersection is None or bridge.second_intersection is None
        assert bridge.manifold_key is not None


# --------------------------------------------------------------------------- #
# 2.3 / 2.6 -- manifold keys are required and propagated
# --------------------------------------------------------------------------- #
def test_base_manifold_requires_a_manifold_key():
    point_a, point_b = (int(row) for row in store_rows([(0.0, 0.0), (1.0, 0.0)]))
    ManifoldMachine._manifold_from_rows(TEST_STORE, [point_a, point_b], "unstable", 1.0, None)

    with pytest.raises(TypeError):
        BaseManifold(point_a, "unstable", 1.0, fixed_point=None, tail=point_b, store=TEST_STORE)

    with pytest.raises(TypeError):
        Bridge(
            root=point_a,
            stability="unstable",
            stretch_param=1.0,
            fixed_point=None,
            tail=point_b,
            store=TEST_STORE,
        )


def test_register_manifold_rejects_a_mismatched_key(initialized):
    workbench, fp = initialized
    manifold = workbench.manifolds[(fp, "unstable", 0, 0)]

    with pytest.raises(AssertionError):
        workbench.register_manifold((fp, "stable", 0, 0), manifold)


def test_iterated_bridge_crossings_carry_the_advanced_key(henon_p3_session):
    """Every crossing born on an image bridge names the image's own branch."""
    session, fp3, _fp1, _zone = henon_p3_session
    workbench = session.workbench
    registry = workbench.intersection_registry

    candidates = [
        b
        for b in workbench.uniiterated_bridges
        if b.fixed_point is fp3 and not b.partial
    ]
    candidates.sort(key=lambda b: registry[b.second_intersection].unstable_cdist)

    for parent in candidates:
        expected = fp3.advance_key(parent.manifold_key, 1)
        known = set(registry.all_ids())
        children = workbench.iterate_bridge(parent)
        new_ids = set(registry.all_ids()) - known
        if not new_ids:
            continue  # image re-traced grown curve; nothing new was born

        for iid in new_ids:
            assert registry[iid].manifold_a_key == expected
        for child in children:
            assert child.manifold_key == expected
        return

    pytest.fail("no period-3 bridge produced new crossings when iterated")


# --------------------------------------------------------------------------- #
# 2.9 -- explicit iterate registration for bridge endpoints
# --------------------------------------------------------------------------- #
def test_bridge_endpoint_iterates_are_registered_explicitly(
    henon_tangle_with_bridges,
):
    """``iterate_bridge`` maps its parent's two endpoints forward and records them."""
    workbench, _fp = henon_tangle_with_bridges
    registry = workbench.intersection_registry
    table = registry.iterate_table
    dynamical_map = workbench.dynamical_system.map

    candidates = sorted(
        workbench.uniiterated_bridges,
        key=lambda b: -registry[b.second_intersection].unstable_cdist,
    )

    for parent in candidates:
        first_id, second_id = parent.first_intersection, parent.second_intersection
        if (first_id, 1) in table or (second_id, 1) in table:
            continue
        known = set(registry.all_ids())
        workbench.iterate_bridge(parent)
        new_ids = set(registry.all_ids()) - known
        if len(new_ids) < 2:
            continue

        images = []
        for src_id in (first_id, second_id):
            target = table[src_id, 1]
            assert target is not None, (
                f"endpoint {src_id} of the iterated bridge got no n=1 iterate"
            )
            assert target in new_ids
            images.append(target)

            predicted = np.asarray(
                dynamical_map(np.asarray(registry[src_id].coords, dtype=float)),
                dtype=float,
            ).ravel()
            distances = {
                iid: float(
                    np.linalg.norm(
                        np.asarray(registry[iid].coords, dtype=float) - predicted
                    )
                )
                for iid in new_ids
            }
            # The registered image is the crossing nearest f(source) -- and by a wide
            # margin: the runner-up is four orders of magnitude further away.
            assert min(distances, key=distances.get) == target
            scale = float(np.linalg.norm(predicted))
            assert distances[target] <= 1e-3 * scale, (
                f"registered image of {src_id} is {distances[target]:.3g} from "
                f"map(source), which is more than the image polyline's resolution"
            )
            assert min(d for i, d in distances.items() if i != target) > 0.1 * scale

            # the image sits on the branches one map step forward, not merely at
            # the right canonical distances
            source = registry[src_id]
            stable_fp = source.manifold_b_key[0]
            assert registry[target].manifold_b_key == stable_fp.advance_key(
                source.manifold_b_key, 1
            )

        # the map is injective, so one bridge's two endpoints cannot share an image
        assert images[0] != images[1]
        return

    pytest.fail("no bridge on the k=10 fixture produced two new crossings")
