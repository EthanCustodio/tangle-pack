"""Phase 6.7 — the resonance zone on the region layer.

``ResonanceZone`` describes its boundary as a list of
:class:`~tanglepack.topology.TopologyResults.Arc` and delegates area and
containment to :mod:`tanglepack.numerics.geometry`, which is the same code the
regions use.
"""

from __future__ import annotations


import numpy as np
import pytest

from cases import build_k10, build_period3
from helpers.fakes import bare_fixed_point
from tanglepack.loom.ResonanceZone import ring_closed, ring_cycles, ring_order
from tanglepack.numerics.Intersection import Intersection
from tanglepack.numerics.geometry import polyline_midpoint
from tanglepack.topology.TopologyResults import Arc


@pytest.fixture
def k10_zone_session() -> tuple:
    """``(session, fp)``: the k=10 tangle with its single resonance zone,
    trimmed at the default strong pip."""
    case = build_k10(through="pips")
    case.session.add_resonance_zones(case.pips)
    return case.session, case.fixed_point


# --------------------------------------------------------------------------- #
# area and containment
# --------------------------------------------------------------------------- #
def test_zone_area_is_positive_and_winding_independent(k10_zone_session):
    """The zone reports a magnitude, whichever way its ring happens to wind."""
    session, _fp = k10_zone_session
    (zone,) = session.resonance_zones.values()
    original = zone.boundary_vertices

    area = zone.area
    assert area > 0.0
    zone.boundary_vertices = original[::-1]
    assert zone.area == pytest.approx(area)
    zone.boundary_vertices = original


def test_zone_containment_includes_the_boundary(k10_zone_session):
    session, _fp = k10_zone_session
    (zone,) = session.resonance_zones.values()
    for vertex in zone.boundary_vertices[::7]:
        assert zone.contains_point(vertex)


def test_zone_of_an_uncaptured_boundary_is_empty(k10_zone_session):
    session, _fp = k10_zone_session
    (zone,) = session.resonance_zones.values()
    zone.boundary_vertices = None
    assert zone.area == 0.0
    assert not zone.contains_point((0.0, 0.0))


# --------------------------------------------------------------------------- #
# the boundary is now described as Arcs
# --------------------------------------------------------------------------- #
def test_boundary_arcs_are_topology_arcs(k10_zone_session):
    session, _fp = k10_zone_session
    (zone,) = session.resonance_zones.values()

    assert zone.boundary_arcs, "capture_boundary must record the arcs"
    for arc in zone.boundary_arcs:
        assert isinstance(arc, Arc)
        assert arc.kind in ("stable", "unstable")
        assert arc.branch_key in session.workbench.manifolds
    kinds = [arc.kind for arc in zone.boundary_arcs]
    assert "stable" in kinds and "unstable" in kinds


def test_boundary_arcs_run_from_the_anchor_to_the_pip(k10_zone_session):
    """Each arc spans its branch from the anchor (cdist 0) out to the cut point."""
    session, _fp = k10_zone_session
    (zone,) = session.resonance_zones.values()
    registry = session.workbench.intersection_registry
    cut = zone.boundary_intersections[0]

    for arc in zone.boundary_arcs:
        assert registry[arc.lo_id].label == "anchor"
        assert arc.lo_id != arc.hi_id
        which = "unstable_cdist" if arc.kind == "unstable" else "stable_cdist"
        assert getattr(registry[arc.hi_id], which) == pytest.approx(
            getattr(cut, which)
        )


def test_boundary_intersection_id_resolves_on_demand(k10_zone_session):
    """The pip is resolved against the live registry, on demand."""
    session, _fp = k10_zone_session
    (zone,) = session.resonance_zones.values()

    resolved = zone.resolve_boundary_intersection_id(session.workbench)
    assert resolved is not None
    registry = session.workbench.intersection_registry
    assert registry[resolved].coords == pytest.approx(
        np.asarray(zone.boundary_intersection.coords)
    )


# --------------------------------------------------------------------------- #
# the blast/classification path still uses the same test point
# --------------------------------------------------------------------------- #
def test_bridge_classification_uses_the_memoised_midpoint(k10_zone_session):
    """A bridge is classified by the midpoint of its points, read afresh on
    every call (the name dates from the walk memo, gone since N6): a bridge
    classified into the zone has its midpoint inside it."""
    session, _fp = k10_zone_session
    (zone,) = session.resonance_zones.values()

    classified = {b: session.classify_bridge(b) for b in session.workbench.bridges}
    assert classified, "the fixture must have bridges to classify"
    inside = [b for b, z in classified.items() if z is zone]
    assert inside, "some bridge must lie in the zone"
    for bridge in inside:
        point = polyline_midpoint(bridge.get_point_array())
        assert zone.contains_point(point)
        assert session.classify_bridge(bridge) is zone


# --------------------------------------------------------------------------- #
# 6.3 -- anchor ids survive a trim + preserve_ids recompute
# --------------------------------------------------------------------------- #
def _anchor_ids(workbench) -> dict:
    """Every anchor's registry id, keyed by the branch pair it joins."""
    return {
        (ix.manifold_a_key, ix.manifold_b_key): iid
        for iid, ix in workbench.intersection_registry
        if ix.label == "anchor"
    }


def test_k10_anchor_ids_survive_a_trim_and_recompute(k10_zone_session):
    """A resonance zone trims and recomputes with ``preserve_ids=True``.

    Anchors are registered first, so they hold the lowest ids; the point of
    registering them rather than detecting them is that those ids stay put. The
    zone fixture has already been through one trim + recompute, so the assertion
    below is on the far side of it.
    """
    session, fp = k10_zone_session
    workbench = session.workbench

    before = _anchor_ids(workbench)
    assert before, "the periodic point must be anchored"
    assert min(before.values()) == 0, "anchors take the lowest ids"

    workbench.compute_intersections([fp], preserve_ids=True)
    assert _anchor_ids(workbench) == before

    # ...and again after a further trim, which changes the crossing set.
    workbench.trim_stable_manifolds(fp)
    workbench.compute_intersections([fp], preserve_ids=True)
    assert _anchor_ids(workbench) == before


def test_p3_anchor_ids_survive_a_trim_and_recompute(henon_p3_session):
    """Same, on the nested period-3 fixture: four anchors across two orbits."""
    session, fp3, fp1, _zone = henon_p3_session
    workbench = session.workbench

    before = _anchor_ids(workbench)
    assert len(before) == fp3.period + fp1.period == 4
    assert sorted(before.values()) == list(range(len(before))), (
        "anchors are registered first, so they occupy ids 0..n-1"
    )

    workbench.compute_intersections([fp3, fp1], preserve_ids=True)
    assert _anchor_ids(workbench) == before

    workbench.trim_stable_manifolds(fp3)
    workbench.trim_stable_manifolds(fp1)
    workbench.compute_intersections([fp3, fp1], preserve_ids=True)
    assert _anchor_ids(workbench) == before


# --------------------------------------------------------------------------- #
# The ring, by full branch key (INV-5)
# --------------------------------------------------------------------------- #
def _cut(fp, unstable: tuple[int, int], stable: tuple[int, int]) -> Intersection:
    """A hand-made cut on ``u(orbit, branch) x s(orbit, branch)``."""
    return Intersection.synthetic(
        (0.0, 0.0), 1.0, 1.0,
        manifold_a_key=(fp, "unstable", *unstable),
        manifold_b_key=(fp, "stable", *stable),
    )


def test_a_period_three_ring_is_one_cycle():
    fp = bare_fixed_point(3)
    cuts = [_cut(fp, (2, 0), (0, 0)), _cut(fp, (0, 0), (1, 0)), _cut(fp, (1, 0), (2, 0))]
    assert ring_cycles(cuts) == [[0, 1, 2]]


def test_an_inversion_point_closes_one_loop_per_branch_pair():
    """Both cuts sit at orbit index 0; the branch index tells the loops apart."""
    fp = bare_fixed_point(inversion=True)
    cuts = [_cut(fp, (0, 0), (0, 0)), _cut(fp, (0, 1), (0, 1))]
    assert ring_cycles(cuts) == [[0], [1]]
    crossed = [_cut(fp, (0, 0), (0, 1)), _cut(fp, (0, 1), (0, 0))]
    assert ring_cycles(crossed) == [[0, 1]]


def test_a_missing_iterate_leaves_the_ring_open():
    fp = bare_fixed_point(3)
    cuts = [_cut(fp, (0, 0), (1, 0)), _cut(fp, (1, 0), (2, 0))]
    assert ring_cycles(cuts) == [[0, 1]]
    assert ring_order(cuts) == []


def test_a_heteroclinic_ring_closes_only_through_both_saddles():
    """The diamond: ``q0' in W^u(B) x W^s(A)`` and ``q0 in W^u(A) x W^s(B)`` close one
    ring by FULL key (fixed point included); either alone stays open."""
    A, B = bare_fixed_point(), bare_fixed_point()
    on_a = Intersection.synthetic((0.0, 0.0), 1.0, 1.0, manifold_a_key=(B, "unstable", 0, 0),
                                  manifold_b_key=(A, "stable", 0, 0))
    on_b = Intersection.synthetic((0.0, 0.0), 1.0, 1.0, manifold_a_key=(A, "unstable", 0, 0),
                                  manifold_b_key=(B, "stable", 0, 0))
    assert ring_closed([on_a, on_b]) and ring_cycles([on_a, on_b]) == [[0, 1]]
    assert not ring_closed([on_a]) and not ring_closed([on_b])
    assert ring_order([on_a]) == []


def test_p3_ring_order_visits_each_stable_branch_once_in_ring_order():
    """The ring runs from a cut's stable branch to the cut on the unstable branch of
    the same orbit point and branch index."""
    session = build_period3().session
    (zone,) = session.resonance_zones.values()
    order = ring_order(zone.boundary_intersections)
    assert sorted(order, key=lambda key: key[2:]) == zone.fixed_point.branch_cycle("stable")
    by_stable = {cut.manifold_b_key: cut for cut in zone.boundary_intersections}
    for key, following in zip(order, order[1:] + order[:1]):
        assert by_stable[following].manifold_a_key[2:] == key[2:]
