"""Pure geometric primitives shared by the numerics, topology and loom layers.

Dev Notes:

Before this module the same four operations existed in three places: the
resonance-zone boundary closure walked a manifold between two canonical
distances, :func:`topology.StablePartition._build_oriented_bridge_polyline`
normalized a bridge's storage order against its dynamical direction,
``ResonanceZone.area`` inlined a shoelace sum, and ``ResonanceZone.contains_point``
inlined an even-odd ray cast. Each copy had its own tolerance and its own
end-inclusion convention. They live here now as pure functions on plain arrays,
so a Region, a resonance zone and a partition midpoint all measure the same way.

Two conventions are fixed here and relied on everywhere downstream:

* **Winding.** :func:`signed_polygon_area` is positive for a counter-clockwise
  ring. The :class:`~tanglepack.topology.Arrangement.Arrangement` traverses its
  bounded faces clockwise, so a closed region's signed area is negative; callers
  that want a magnitude take ``abs``.
* **Boundary inclusion.** :func:`point_in_polygon` counts a point ON the boundary
  (within ``tol``) as inside. That is deliberate and load-bearing: a bridge that
  forms a resonance zone's own unstable boundary arc must be classified into the
  zone it bounds, and a region's representative point may land on a shared arc.

Everything here works on the memoised point arrays
(:meth:`~tanglepack.numerics.BaseManifold.BaseManifold.get_point_array`), so
repeated calls on an unchanged bridge cost one dict lookup.
"""

from __future__ import annotations

import logging
from typing import Optional, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from .Intersection import Stability

if TYPE_CHECKING:
    from .BaseManifold import BaseManifold
    from .Bridge import Bridge

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def arc_polyline(
    manifold: "BaseManifold",
    lo_cdist: float,
    hi_cdist: float,
    *,
    stability: Optional[Stability] = None,
    reverse: bool = False,
    tol: float = 0.0,
) -> NDArray[np.float64]:
    """
    The manifold nodes lying between two canonical distances, in curve order.

    The nodes of a manifold run root to tail in increasing canonical distance,
    so an arc between two crossings is a contiguous slice of them. Only real
    nodes are returned -- the crossings themselves are not manifold points, so a
    caller that wants an arc closed at its two crossings prepends/appends their
    coordinates (:meth:`~tanglepack.topology.TopologyResults.Arc.polyline` does).

    Args:
        manifold: The curve to walk.
        lo_cdist: Lower canonical distance of the window.
        hi_cdist: Upper canonical distance of the window.
        stability: Which canonical distance to read. Defaults to the manifold's
            own stability, which is what every caller but a diagnostic wants.
        reverse: If True the polyline is returned in decreasing cdist order.
        tol: Slack added to both ends of the window, so a node sitting exactly on
            a boundary is not lost to rounding. Zero by default.

    Returns:
        An ``(N, 2)`` array of coordinates in increasing canonical distance
        (decreasing when ``reverse``); empty when no node falls in the window.
        The ascending order is guaranteed regardless of the manifold's storage
        direction, so a bridge stored tail-first comes back the same way round as
        the parent it was cut from.
    """
    which = stability if stability is not None else manifold.stability
    nodes = manifold.get_point_array(return_nodes=True)
    lo, hi = (lo_cdist, hi_cdist) if lo_cdist <= hi_cdist else (hi_cdist, lo_cdist)
    inside = [
        (node.get_cdist(which), node.get_point())
        for node in nodes
        if lo - tol <= node.get_cdist(which) <= hi + tol
    ]
    if not inside:
        return np.empty((0, 2), dtype=np.float64)
    # A walk runs root to tail, which is monotone in cdist but not necessarily
    # ASCENDING: a bridge cut off a manifold can be stored either way round. One
    # end comparison is enough to normalize it, and it keeps the guarantee in the
    # Returns section true for every manifold-like input.
    if inside[0][0] > inside[-1][0]:
        inside.reverse()
    poly = np.asarray([point for _cdist, point in inside], dtype=np.float64)
    return poly[::-1] if reverse else poly


def oriented_bridge_polyline(
    bridge: "Bridge",
    first_cdist: float,
    second_cdist: float,
    *,
    tol: float = 0.0,
) -> Optional[NDArray[np.float64]]:
    """
    A bridge's polyline oriented by the unstable dynamical direction.

    A bridge is unstable manifold, so its dynamical direction -- the direction of
    forward flow -- points away from the fixed point, i.e. along increasing
    unstable canonical distance. Storage order (root to tail) need not match it;
    this is the single place that mismatch is normalized.

    The two endpoint canonical distances are passed in rather than looked up, so
    the function stays pure: the caller (a Trellis, an Arrangement) is the one
    that owns the registry the endpoint ids resolve against.

    Args:
        bridge: The bridge whose polyline to orient.
        first_cdist: Unstable canonical distance of the bridge's ``first``
            endpoint crossing.
        second_cdist: Unstable canonical distance of the ``second`` endpoint.
        tol: Canonical-distance tolerance below which the two endpoints count as
            equal and the orientation is undecidable.

    Returns:
        The ``(N, 2)`` polyline running from the lower-cdist end to the higher,
        or ``None`` when the bridge has fewer than two points or its endpoints
        have equal unstable canonical distance (logged as a warning -- a bridge
        of zero unstable span is degenerate, see the anchor Dev Notes in
        :mod:`tanglepack.numerics.TangleWorkbench`).
    """
    points = bridge.get_point_array()
    if points is None or len(points) < 2:
        return None
    if abs(first_cdist - second_cdist) <= tol:
        logger.warning(
            "Bridge (%s, %s) has endpoints of equal unstable cdist; dynamical "
            "orientation is undecidable",
            bridge.first_intersection,
            bridge.second_intersection,
        )
        return None
    poly = np.asarray(points, dtype=np.float64)
    return poly if first_cdist < second_cdist else poly[::-1]


def signed_polygon_area(points: NDArray[np.float64]) -> float:
    """
    The shoelace area of a polygon, signed by its winding.

    Args:
        points: An ``(N, 2)`` ring of vertices. A repeated closing vertex is
            harmless -- the degenerate final edge contributes nothing.

    Returns:
        Positive for a counter-clockwise ring, negative for a clockwise one, and
        ``0.0`` for fewer than three vertices.
    """
    verts = np.asarray(points, dtype=np.float64)
    if verts.ndim != 2 or len(verts) < 3:
        return 0.0
    x, y = verts[:, 0], verts[:, 1]
    return float(0.5 * (np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def point_in_polygon(
    point: "NDArray[np.float64] | tuple[float, float]",
    polygon: NDArray[np.float64],
    *,
    tol: float = 1e-9,
) -> bool:
    """
    Whether a point lies inside a polygon, boundary INCLUDED.

    Inside-ness is decided by an even-odd ray cast to the right of the point,
    which is winding-independent. A point on the boundary is what a ray cast
    cannot decide reliably, so it is detected separately (distance to the nearest
    edge within ``tol``) and always reported as contained -- see the module Dev
    Notes for why that inclusion is load-bearing.

    Args:
        point: The ``(x, y)`` to test.
        polygon: An ``(N, 2)`` ring. It need not be explicitly closed; the edge
            from the last vertex back to the first is always considered.
        tol: Distance below which the point counts as lying on the boundary.

    Returns:
        True if the point is strictly inside or on the boundary; False otherwise,
        including for a degenerate polygon of fewer than three vertices.
    """
    verts = np.asarray(polygon, dtype=np.float64)
    if verts.ndim != 2 or len(verts) < 3:
        return False
    px, py = float(point[0]), float(point[1])

    x1, y1 = verts[:, 0], verts[:, 1]
    x2, y2 = np.roll(x1, -1), np.roll(y1, -1)

    # On any edge (within tol)? The boundary is inclusive.
    dx, dy = x2 - x1, y2 - y1
    seg_len2 = dx * dx + dy * dy
    with np.errstate(invalid="ignore", divide="ignore"):
        t = np.where(seg_len2 > 0, ((px - x1) * dx + (py - y1) * dy) / seg_len2, 0.0)
    t = np.clip(t, 0.0, 1.0)
    cx, cy = x1 + t * dx, y1 + t * dy
    if np.any((cx - px) ** 2 + (cy - py) ** 2 <= tol * tol):
        return True

    # Even-odd ray cast to the right of the point.
    crosses = ((y1 > py) != (y2 > py)) & (
        px < (x2 - x1) * (py - y1) / np.where(y2 != y1, y2 - y1, 1.0) + x1
    )
    return bool(np.count_nonzero(crosses) % 2 == 1)


def polyline_midpoint(
    points: NDArray[np.float64],
) -> Optional[NDArray[np.float64]]:
    """
    The middle NODE of a polyline.

    Deliberately the index-middle vertex rather than the arclength midpoint: the
    result is guaranteed to lie exactly ON the curve, which is what makes a
    bridge's representative point land on a resonance zone's boundary arc when
    the bridge IS that arc (the classification relies on it). An interpolated
    arclength midpoint would sit a rounding error off the chord.

    Args:
        points: An ``(N, 2)`` polyline.

    Returns:
        A fresh ``(2,)`` array holding the middle vertex, or ``None`` for an empty
        polyline. It is a COPY: the input is routinely a memoised (and read-only)
        point array shared with every other caller, and a view into it would either
        refuse to be written or, for a writable input, let one caller corrupt the
        memo for the rest.
    """
    verts = np.asarray(points, dtype=np.float64)
    if verts.ndim != 2 or len(verts) == 0:
        return None
    return np.array(verts[len(verts) // 2], dtype=np.float64)


def polygon_interior_point(
    polygon: NDArray[np.float64],
) -> Optional[NDArray[np.float64]]:
    """
    A point guaranteed to lie strictly inside a simple polygon.

    Unlike a centroid or a mean of boundary points, this cannot fall outside a
    concave or crescent-shaped face -- and the faces of a tangle are routinely
    crescents. A horizontal scanline is placed strictly BETWEEN two consecutive
    distinct vertex ordinates (so no vertex sits on it and the crossing count is
    even), its crossings with the boundary are sorted, and the midpoint of the
    widest interior span is returned. Scanlines are tried from the middle of the
    polygon outward, which keeps the point away from the thin ends.

    The polygon must be simple (non-self-intersecting). A face of a planar
    arrangement always is: its boundary is a cycle of arcs that, by the
    fundamental invariant, never cross one another.

    Args:
        polygon: An ``(N, 2)`` ring. It need not be explicitly closed.

    Returns:
        A ``(2,)`` interior point, or ``None`` for a degenerate polygon (fewer
        than three vertices, or no scanline with an interior span).
    """
    verts = np.asarray(polygon, dtype=np.float64)
    if verts.ndim != 2 or len(verts) < 3:
        return None

    ys = np.unique(verts[:, 1])
    if len(ys) < 2:
        return None
    scanlines = 0.5 * (ys[:-1] + ys[1:])
    order = np.argsort(np.abs(scanlines - float(np.median(verts[:, 1]))))

    x1, y1 = verts[:, 0], verts[:, 1]
    x2, y2 = np.roll(x1, -1), np.roll(y1, -1)
    for y in scanlines[order]:
        straddles = (y1 > y) != (y2 > y)
        if not np.any(straddles):
            continue
        a_x, a_y = x1[straddles], y1[straddles]
        b_x, b_y = x2[straddles], y2[straddles]
        xs = np.sort(a_x + (y - a_y) * (b_x - a_x) / (b_y - a_y))
        if len(xs) < 2:
            continue
        # Interior spans of an even-odd ray cast are the pairs (xs[0], xs[1]),
        # (xs[2], xs[3]), ...; take the widest so the point sits well inside.
        widths = xs[1::2] - xs[0::2]
        if not len(widths):
            continue
        best = int(np.argmax(widths))
        if widths[best] <= 0.0:
            continue
        return np.array(
            [0.5 * (xs[2 * best] + xs[2 * best + 1]), float(y)], dtype=np.float64
        )
    return None
