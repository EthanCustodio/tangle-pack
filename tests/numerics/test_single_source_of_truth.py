"""Phase 2 Task A — single sources of truth in the numerics core.

Two kinds of test live here.

**Pin tests** record today's bridge-cutting result as expected data: for every
bridge, its unstable ``(orbit_index, branch_index)`` and the two canonical
distances of each of its two bounding crossings, rounded to six significant
digits, plus the registry size. They were captured against the pre-refactor
code (crossings resolved into ``Tangle._intersecting_coords``, bridge endpoints
assigned afterwards by nearest-cdist lookup) and must be reproduced exactly by
the registry-only cutting path.

**Acceptance tests** pin the new invariants: resolved crossings live only in the
``IntersectionRegistry``; a bridge's endpoints are the registry ids of the two
crossings it was cut at; every manifold and every bridge carries its
``manifold_key``; and the forward iterate of a bridge's two endpoint crossings is
registered explicitly rather than inferred.
"""

from __future__ import annotations

from math import floor, log10

import numpy as np
import pytest

from tanglepack.numerics.BaseManifold import BaseManifold
from tanglepack.numerics.Bridge import Bridge
from tanglepack.numerics.Point import Point


# --------------------------------------------------------------------------- #
# fingerprint helpers
# --------------------------------------------------------------------------- #
def _sig(value: float | None, digits: int = 6) -> float | None:
    """Round to ``digits`` significant digits (None passes through)."""
    if value is None:
        return None
    value = float(value)
    if value == 0.0:
        return 0.0
    return round(value, -int(floor(log10(abs(value)))) + (digits - 1))


def bridge_fingerprint(workbench) -> list[tuple]:
    """Sorted ``(orbit/branch, first cdists, second cdists)`` of every bridge."""
    registry = workbench.intersection_registry
    rows = []
    for bridge in workbench.bridges:
        first_id, second_id = bridge.first_intersection, bridge.second_intersection
        first = registry[first_id] if first_id in registry else None
        second = registry[second_id] if second_id in registry else None
        rows.append(
            (
                bridge.manifold_key[2:] if bridge.manifold_key is not None else None,
                _sig(first.unstable_cdist) if first else None,
                _sig(first.stable_cdist) if first else None,
                _sig(second.unstable_cdist) if second else None,
                _sig(second.stable_cdist) if second else None,
            )
        )
    return sorted(rows, key=repr)


# --------------------------------------------------------------------------- #
# expected data, captured against the pre-refactor code (2026-09-02)
# --------------------------------------------------------------------------- #
K10_REGISTRY_SIZE = 8
K10_BRIDGE_ENDPOINTS = [
    ((0, 0), 0.0, 0.0, 9.90728, 84.3686),
    ((0, 0), 631.412, 8.70701, 651.766, 76.5359),
    ((0, 0), 651.766, 76.5359, 697.225, 81.8741),
    ((0, 0), 697.225, 81.8741, 718.466, 1.16361),
    ((0, 0), 74.146, 74.146, 84.3686, 9.90728),
    ((0, 0), 84.3686, 9.90728, 631.412, 8.70701),
    ((0, 0), 9.90728, 84.3686, 74.146, 74.146),
]

P3_REGISTRY_SIZE = 34
P3_BRIDGE_ENDPOINTS = [
    ((0, 0), 0.0, 0.0, 0.615371, 1.48055),
    ((0, 0), 0.0, 0.0, 37.1964, 7.05211),
    ((0, 0), 0.615371, 1.48055, 1.19459, 0.770302),
    ((0, 0), 1.19459, 0.770302, 2.29603, 0.3968),
    ((0, 0), 155.287, 5.58106, 196.194, 1.33696),
    ((0, 0), 16.6409, 0.0552914, 31.9846, 0.0284734),
    ((0, 0), 2.29603, 0.3968, 4.45774, 0.206421),
    ((0, 0), 31.9846, 0.0284734, 62.1091, 0.0148128),
    ((0, 0), 37.1964, 7.05211, 155.287, 5.58106),
    ((0, 0), 4.45774, 0.206421, 8.57387, 0.106253),
    ((0, 0), 62.1091, 0.0148128, 119.366, 0.00762394),
    ((0, 0), 8.57387, 0.106253, 16.6409, 0.0552914),
    ((1, 0), 0.0, 0.0, 0.954405, 0.954604),
    ((1, 0), 0.954405, 0.954604, 1.8529, 0.496622),
    ((1, 0), 1.8529, 0.496622, 3.56148, 0.255809),
    ((1, 0), 13.288, 0.0685557, 25.8148, 0.0356415),
    ((1, 0), 25.8148, 0.0356415, 49.6083, 0.0183569),
    ((1, 0), 3.56148, 0.255809, 6.918, 0.13301),
    ((1, 0), 49.6083, 0.0183569, 96.3395, 0.00954955),
    ((1, 0), 6.918, 0.13301, 13.288, 0.0685557),
    ((1, 0), 96.3395, 0.00954955, 185.15, 0.00491442),
    ((2, 0), 0.0, 0.0, 1.48041, 0.61542),
    ((2, 0), 1.48041, 0.61542, 2.87409, 0.320164),
    ((2, 0), 10.7305, 0.0857516, 20.6125, 0.0441966),
    ((2, 0), 149.435, 0.00615608, 287.189, 0.003164),
    ((2, 0), 2.87409, 0.320164, 5.5245, 0.164912),
    ((2, 0), 20.6125, 0.0441966, 40.042, 0.0229778),
    ((2, 0), 40.042, 0.0229778, 76.9492, 0.0118349),
    ((2, 0), 5.5245, 0.164912, 10.7305, 0.0857516),
    ((2, 0), 76.9492, 0.0118349, 149.435, 0.00615608),
]


# --------------------------------------------------------------------------- #
# pin tests
# --------------------------------------------------------------------------- #
def test_bridge_cutting_pin_k10(henon_tangle_with_bridges):
    """Registry-only cutting agrees with the pre-refactor result (k=10)."""
    workbench, _fp = henon_tangle_with_bridges

    assert len(workbench.intersection_registry) == K10_REGISTRY_SIZE
    assert bridge_fingerprint(workbench) == K10_BRIDGE_ENDPOINTS


@pytest.mark.slow
def test_bridge_cutting_pin_period_3(henon_p3_session):
    """Registry-only cutting agrees with the pre-refactor result (nested p3)."""
    session, _fp3, _fp1, _zone = henon_p3_session
    workbench = session.workbench

    assert len(workbench.intersection_registry) == P3_REGISTRY_SIZE
    assert bridge_fingerprint(workbench) == P3_BRIDGE_ENDPOINTS


# --------------------------------------------------------------------------- #
# 2.1 -- resolved crossings live in the registry only
# --------------------------------------------------------------------------- #
def test_tangle_keeps_only_index_state(small_tangle):
    workbench, _fp = small_tangle
    tangle = workbench.Tangle

    for gone in (
        "_intersections",
        "_intersecting_coords",
        "_intersecting_points",
        "_intersection_by_seg",
        "iter_intersection_coords",
    ):
        assert not hasattr(tangle, gone), f"Tangle still owns {gone}"

    for kept in (
        "_intersecting_segments",
        "_seg_lookup",
        "_manifold_segs",
        "_seg_manifolds",
        "_edge_to_sid",
        "_processed_pairs",
        "_rtree",
    ):
        assert hasattr(tangle, kept), f"Tangle lost its index state {kept}"


def test_compute_intersections_returns_one_coord_per_registered_crossing(grown_both):
    """The returned coords are the registry's, read back independently.

    Independently: the expected list is rebuilt from the registry BY ID rather than
    by re-evaluating the comprehension the method returns, so the test would catch
    a return value that dropped, duplicated, or reordered a crossing.
    """
    workbench, fp = grown_both
    coords = workbench.compute_intersections([fp], infer_iterates=False)
    registry = workbench.intersection_registry

    assert len(coords) == len(registry) > 0
    for coord, iid in zip(coords, registry.all_ids()):
        assert coord == registry[iid].coords


def test_every_registered_crossing_carries_its_unstable_segment(small_tangle):
    """The bracketing points ``create_bridges`` cuts between are on the crossing."""
    workbench, _fp = small_tangle

    for _iid, ix in workbench.intersection_registry:
        assert ix.unstable_manifold is not None
        assert ix.unstable_segment is not None
        p0, p1 = ix.unstable_segment
        lo, hi = sorted((p0.get_cdist("unstable"), p1.get_cdist("unstable")))
        assert lo - 1e-9 <= ix.unstable_cdist <= hi + 1e-9


def test_trim_stable_manifolds_reads_the_registry(grown_both):
    """The new tail is the first node at or past the outermost crossing."""
    workbench, fp = grown_both
    workbench.compute_intersections([fp], infer_iterates=False)

    outermost = max(ix.stable_cdist for _iid, ix in workbench.intersection_registry)
    workbench.trim_stable_manifolds(fp)

    manifold = workbench.manifolds[(fp, "stable", 0, 0)]
    assert manifold.tail.get_cdist("stable") >= outermost


# --------------------------------------------------------------------------- #
# 2.2 -- bridge endpoints set at cut time
# --------------------------------------------------------------------------- #
def test_bridge_endpoints_are_the_ids_of_the_crossings_it_was_cut_at(
    henon_tangle_with_bridges,
):
    workbench, _fp = henon_tangle_with_bridges
    registry = workbench.intersection_registry

    assert not hasattr(workbench, "_assign_bridge_intersections")
    assert not hasattr(workbench, "_endpoint_candidates")

    for bridge in workbench.bridges:
        assert not bridge.partial
        assert isinstance(bridge.first_intersection, int)
        assert isinstance(bridge.second_intersection, int)
        first = registry[bridge.first_intersection]
        second = registry[bridge.second_intersection]
        assert first.unstable_cdist < second.unstable_cdist
        root_u = bridge.root.get_cdist("unstable")
        tail_u = bridge.tail.get_cdist("unstable")
        assert root_u <= first.unstable_cdist
        assert second.unstable_cdist <= tail_u


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
    point_a = Point(0.0, 0.0, 0.0)
    point_b = Point(1.0, 0.0, 1.0)
    point_a.insert_point_forward(point_b)

    with pytest.raises(TypeError):
        BaseManifold(point_a, "unstable", 1.0, fixed_point=None, tail=point_b)

    with pytest.raises(TypeError):
        Bridge(
            root=point_a,
            stability="unstable",
            stretch_param=1.0,
            fixed_point=None,
            tail=point_b,
        )


def test_register_manifold_rejects_a_mismatched_key(initialized):
    workbench, fp = initialized
    manifold = workbench.manifolds[(fp, "unstable", 0, 0)]

    with pytest.raises(AssertionError):
        workbench.register_manifold((fp, "stable", 0, 0), manifold)


@pytest.mark.slow
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
