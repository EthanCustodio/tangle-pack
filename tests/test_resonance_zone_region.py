"""Phase 6.7 — the resonance zone on the region layer.

``ResonanceZone`` used to inline its own shoelace sum and its own ray cast, and
built its boundary as a bare vertex array. It now describes its boundary as a list
of :class:`~tanglepack.topology.TopologyResults.Arc` and delegates area and
containment to :mod:`tanglepack.numerics.geometry`, which is the same code the
regions use.

The point of this file is that NOTHING a caller can observe changed. The expected
values below were captured against the pre-6.7 implementation (2026-09-03), before
the anchors became deliberate and before the geometry moved, and are asserted
against the refactored one:

* ``zone.area`` to twelve significant digits;
* the in/out answer on an 11x11 grid over each zone's own bounding box, expanded
  by 10% — a signature dense enough that any change to the ray cast, the boundary
  stitching or the vertex list shows up as a different bit string.

The bounding box is derived from the zone's own captured boundary, so the grid
follows the geometry rather than hard-coding coordinates.
"""

from __future__ import annotations

import numpy as np
import pytest

from tanglepack import TangleSession
from tanglepack.loom.ResonanceZone import BoundaryArc
from tanglepack.topology.TopologyResults import Arc


# --------------------------------------------------------------------------- #
# expected data, captured 2026-09-03 against the pre-6.7 implementation
# --------------------------------------------------------------------------- #
K10_ZONE_AREA = 97.9399955713
K10_ZONE_MASK = (
    "0000000000000000000000001111000000011111100000111111100001111111100011111110"
    "001111111000011110000000100000000000000000000"
)

P3_INNER_ZONE_AREA = 0.538178460497
P3_INNER_ZONE_MASK = (
    "0000000000000011100000001111111000001111110000000111100000000111100000000111"
    "000000001110000000001000000000000000000000000"
)

P3_OUTER_ZONE_AREA = 21.0625044173
P3_OUTER_ZONE_MASK = (
    "0000000000000000100000000011100000001111100000011111000000111111000011111110"
    "000111111100001111110000110000000000000000000"
)


def bbox_grid(zone, n: int = 11, pad: float = 0.1) -> list[tuple[float, float]]:
    """An ``n`` x ``n`` grid over the zone's own bounding box, expanded by ``pad``."""
    verts = zone.boundary_vertices
    lo, hi = verts.min(axis=0), verts.max(axis=0)
    span = hi - lo
    lo, hi = lo - pad * span, hi + pad * span
    return [
        (float(x), float(y))
        for x in np.linspace(lo[0], hi[0], n)
        for y in np.linspace(lo[1], hi[1], n)
    ]


def containment_mask(zone) -> str:
    """The zone's in/out answer over :func:`bbox_grid`, as a bit string."""
    return "".join("1" if zone.contains_point(p) else "0" for p in bbox_grid(zone))


@pytest.fixture
def k10_zone_session(henon_map, henon_map_inverse):
    """A k=10 session with its single resonance zone defined."""
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
    trellis = session.trellis(fp)
    trellis.classify_strong_pips()
    session.add_resonance_zones([trellis.strong_pip])
    return session, fp


# --------------------------------------------------------------------------- #
# area and containment are unchanged
# --------------------------------------------------------------------------- #
def test_k10_zone_area_and_containment_are_unchanged(k10_zone_session):
    session, _fp = k10_zone_session
    (zone,) = session.resonance_zones.values()

    assert zone.area == pytest.approx(K10_ZONE_AREA, rel=1e-11)
    assert containment_mask(zone) == K10_ZONE_MASK


@pytest.mark.slow
def test_p3_zone_areas_and_containment_are_unchanged(henon_p3_session):
    session, _fp3, _fp1, _zone = henon_p3_session
    inner, outer = sorted(session.resonance_zones.values(), key=lambda z: z.area)

    assert inner.area == pytest.approx(P3_INNER_ZONE_AREA, rel=1e-11)
    assert containment_mask(inner) == P3_INNER_ZONE_MASK
    assert outer.area == pytest.approx(P3_OUTER_ZONE_AREA, rel=1e-11)
    assert containment_mask(outer) == P3_OUTER_ZONE_MASK


def test_zone_area_is_positive_and_winding_independent(k10_zone_session):
    """The zone reports a magnitude, whichever way its ring happens to wind."""
    session, _fp = k10_zone_session
    (zone,) = session.resonance_zones.values()
    original = zone.boundary_vertices

    assert zone.area > 0.0
    zone.boundary_vertices = original[::-1]
    assert zone.area == pytest.approx(K10_ZONE_AREA, rel=1e-11)
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
    """The stored id is gone; the pip is resolved against the live registry."""
    session, _fp = k10_zone_session
    (zone,) = session.resonance_zones.values()

    assert not hasattr(zone, "boundary_intersection_id")
    resolved = zone.resolve_boundary_intersection_id(session.workbench)
    assert resolved is not None
    registry = session.workbench.intersection_registry
    assert registry[resolved].coords == pytest.approx(
        np.asarray(zone.boundary_intersection.coords)
    )


def test_boundary_arc_dataclass_is_gone():
    """The zone-local arc record is replaced by the shared topology Arc."""
    assert BoundaryArc is Arc


# --------------------------------------------------------------------------- #
# the blast/classification path still uses the same test point
# --------------------------------------------------------------------------- #
def test_bridge_classification_uses_the_memoised_midpoint(k10_zone_session):
    session, _fp = k10_zone_session
    (zone,) = session.resonance_zones.values()

    classified = session.classify_bridges()
    assert classified, "the fixture must have bridges to classify"
    inside = [b for b, z in classified.items() if z is zone]
    assert inside, "some bridge must lie in the zone"
    for bridge in inside:
        point = session._bridge_test_point(bridge)
        assert zone.contains_point(point)


def test_shapely_is_not_a_dependency():
    """Phase 6.7 drops shapely: nothing imports it and it is not declared."""
    import pathlib
    import re

    root = pathlib.Path(__file__).resolve().parents[1]
    pyproject = (root / "pyproject.toml").read_text()
    assert "shapely" not in pyproject.lower()

    offenders = [
        path
        for path in (root / "src").rglob("*.py")
        if re.search(r"\bshapely\b", path.read_text())
    ]
    assert offenders == []


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


@pytest.mark.slow
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
