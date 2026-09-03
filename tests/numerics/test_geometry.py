"""Phase 6.1 — the geometry primitives, against analytic shapes.

``numerics/geometry.py`` consolidates the polyline/polygon helpers that used to
live in three places (the resonance-zone boundary closure, the stable-partition
midpoint helpers, ``ResonanceZone.area`` / ``contains_point``). They are pure
functions on plain arrays, so they are pinned here against shapes whose area and
containment are known exactly rather than against a tangle.
"""

from __future__ import annotations

import numpy as np
import pytest

from tanglepack.numerics.geometry import (
    arc_polyline,
    oriented_bridge_polyline,
    point_in_polygon,
    polygon_interior_point,
    polyline_midpoint,
    signed_polygon_area,
)


SQUARE = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
TRIANGLE = np.array([[0.0, 0.0], [4.0, 0.0], [0.0, 3.0]])
# An L: the unit 2x2 square with the top-right 1x1 corner bitten out (area 3).
CONCAVE = np.array(
    [[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [1.0, 1.0], [1.0, 2.0], [0.0, 2.0]]
)


# --------------------------------------------------------------------------- #
# signed_polygon_area
# --------------------------------------------------------------------------- #
def test_signed_area_of_a_ccw_square_is_positive_one():
    assert signed_polygon_area(SQUARE) == pytest.approx(1.0)


def test_signed_area_flips_with_winding():
    assert signed_polygon_area(SQUARE[::-1]) == pytest.approx(-1.0)


def test_signed_area_ignores_an_explicit_closing_vertex():
    closed = np.vstack([SQUARE, SQUARE[:1]])
    assert signed_polygon_area(closed) == pytest.approx(signed_polygon_area(SQUARE))


def test_signed_area_of_a_triangle():
    assert signed_polygon_area(TRIANGLE) == pytest.approx(6.0)


def test_signed_area_of_a_concave_polygon():
    assert signed_polygon_area(CONCAVE) == pytest.approx(3.0)


def test_signed_area_of_a_degenerate_polygon_is_zero():
    assert signed_polygon_area(np.array([[0.0, 0.0], [1.0, 1.0]])) == 0.0
    assert signed_polygon_area(np.empty((0, 2))) == 0.0


# --------------------------------------------------------------------------- #
# point_in_polygon
# --------------------------------------------------------------------------- #
def test_point_in_polygon_interior_and_exterior():
    assert point_in_polygon((0.5, 0.5), SQUARE)
    assert not point_in_polygon((1.5, 0.5), SQUARE)
    assert not point_in_polygon((0.5, -0.5), SQUARE)


def test_point_in_polygon_on_edge_and_on_vertex_count_as_inside():
    # Documented policy: the boundary is INCLUSIVE (see the module docstring) --
    # a bridge lying on its own zone's boundary arc must classify as inside.
    assert point_in_polygon((0.5, 0.0), SQUARE)
    assert point_in_polygon((0.0, 0.0), SQUARE)
    assert point_in_polygon((1.0, 1.0), SQUARE)


def test_point_in_polygon_respects_the_concavity():
    assert point_in_polygon((0.5, 1.5), CONCAVE)
    assert not point_in_polygon((1.5, 1.5), CONCAVE)


def test_point_in_polygon_is_winding_independent():
    for poly in (CONCAVE, CONCAVE[::-1]):
        assert point_in_polygon((0.5, 1.5), poly)
        assert not point_in_polygon((1.5, 1.5), poly)


def test_point_in_polygon_accepts_an_explicitly_closed_ring():
    closed = np.vstack([CONCAVE, CONCAVE[:1]])
    assert point_in_polygon((0.5, 1.5), closed)
    assert not point_in_polygon((1.5, 1.5), closed)


def test_point_in_polygon_of_a_degenerate_polygon_is_false():
    assert not point_in_polygon((0.0, 0.0), np.array([[0.0, 0.0], [1.0, 1.0]]))


# --------------------------------------------------------------------------- #
# polyline_midpoint
# --------------------------------------------------------------------------- #
def test_polyline_midpoint_is_the_middle_node():
    pts = np.array([[0.0, 0.0], [1.0, 0.0], [5.0, 0.0]])
    assert np.allclose(polyline_midpoint(pts), [1.0, 0.0])


def test_polyline_midpoint_of_an_empty_polyline_is_none():
    assert polyline_midpoint(np.empty((0, 2))) is None


def test_polyline_midpoint_returns_a_copy_not_a_view():
    """The input is routinely a memoised, shared, read-only point array."""
    pts = np.array([[0.0, 0.0], [1.0, 0.0], [5.0, 0.0]])
    pts.flags.writeable = False
    mid = polyline_midpoint(pts)
    mid[0] = 99.0  # must not raise, and must not touch the source
    assert pts[1][0] == 1.0


# --------------------------------------------------------------------------- #
# polygon_interior_point
# --------------------------------------------------------------------------- #
# An almost-closed annulus sector -- a "C". Its vertex mean sits at the origin, in
# the HOLE, so neither the mean of its vertices nor its centroid is inside it. This
# is the shape a tangle face takes when one short stable arc closes a long
# meandering bridge, and it is why the scanline construction exists at all.
_THETA = np.linspace(-0.9 * np.pi, 0.9 * np.pi, 60)
CRESCENT = np.vstack(
    [
        np.column_stack([np.cos(_THETA), np.sin(_THETA)]),
        0.6 * np.column_stack([np.cos(_THETA[::-1]), np.sin(_THETA[::-1])]),
    ]
)


@pytest.mark.parametrize(
    "polygon", [SQUARE, TRIANGLE, CONCAVE, CRESCENT], ids="square triangle L crescent".split()
)
def test_interior_point_is_inside(polygon):
    point = polygon_interior_point(polygon)
    assert point is not None
    assert point_in_polygon(point, polygon)


def test_interior_point_is_winding_independent():
    for polygon in (CONCAVE, CONCAVE[::-1], CRESCENT, CRESCENT[::-1]):
        assert point_in_polygon(polygon_interior_point(polygon), polygon)


def test_interior_point_beats_the_vertex_mean_on_a_crescent():
    """The cheap heuristic fails here; the scanline is why the fallback exists."""
    assert not point_in_polygon(CRESCENT.mean(axis=0), CRESCENT)
    assert point_in_polygon(polygon_interior_point(CRESCENT), CRESCENT)


def test_interior_point_of_a_degenerate_polygon_is_none():
    assert polygon_interior_point(np.array([[0.0, 0.0], [1.0, 1.0]])) is None
    # A horizontal sliver has only one distinct ordinate: no scanline exists.
    assert (
        polygon_interior_point(
            np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
        )
        is None
    )


# --------------------------------------------------------------------------- #
# arc_polyline / oriented_bridge_polyline on real manifolds
# --------------------------------------------------------------------------- #
def test_arc_polyline_clips_to_the_cdist_window(grown_unstable):
    _workbench, _fp, manifold = grown_unstable
    nodes = manifold.get_point_array(return_nodes=True)
    cdists = np.array([n.get_cdist("unstable") for n in nodes])
    lo, hi = float(cdists[2]), float(cdists[8])

    arc = arc_polyline(manifold, lo, hi)

    assert arc.shape[1] == 2
    assert len(arc) == 7
    assert np.allclose(arc[0], nodes[2].get_point())
    assert np.allclose(arc[-1], nodes[8].get_point())


def test_arc_polyline_reverse_flips_the_order(grown_unstable):
    _workbench, _fp, manifold = grown_unstable
    nodes = manifold.get_point_array(return_nodes=True)
    cdists = [n.get_cdist("unstable") for n in nodes]
    lo, hi = float(cdists[2]), float(cdists[8])

    forward = arc_polyline(manifold, lo, hi)
    backward = arc_polyline(manifold, lo, hi, reverse=True)

    assert np.allclose(forward, backward[::-1])


def test_arc_polyline_on_an_empty_window_is_empty(grown_unstable):
    _workbench, _fp, manifold = grown_unstable
    assert len(arc_polyline(manifold, -5.0, -4.0)) == 0


def test_oriented_bridge_polyline_runs_in_the_dynamical_direction(
    henon_tangle_with_bridges,
):
    workbench, _fp = henon_tangle_with_bridges
    registry = workbench.intersection_registry

    for bridge in workbench.bridges:
        if bridge.id is None:
            continue
        first, second = bridge.id
        poly = oriented_bridge_polyline(
            bridge,
            registry[first].unstable_cdist,
            registry[second].unstable_cdist,
            tol=registry.cdist_tol,
        )
        if poly is None:
            continue
        raw = bridge.get_point_array()
        # Either the storage order or its reverse, and always low-cdist first.
        assert np.allclose(poly, raw) or np.allclose(poly, raw[::-1])
        nodes = bridge.get_point_array(return_nodes=True)
        ordered = [n.get_cdist("unstable") for n in nodes]
        if ordered[0] > ordered[-1]:
            ordered = ordered[::-1]
        assert ordered == sorted(ordered)


def test_oriented_bridge_polyline_refuses_equal_cdist_endpoints(
    henon_tangle_with_bridges,
):
    workbench, _fp = henon_tangle_with_bridges
    bridge = next(b for b in workbench.bridges if b.id is not None)
    assert oriented_bridge_polyline(bridge, 1.0, 1.0, tol=1e-6) is None
