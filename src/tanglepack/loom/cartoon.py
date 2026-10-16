"""
What the drawings need from the resonance zones.

:func:`cartoon_zones` reads the stored zones for the circular dual-graph
cartoon (which side of each stable branch faces its zone, which zone nests in
which, and how each zone's arcs go around it), and
:func:`plot_resonance_zones` shades them. Both take loom
:class:`~.ResonanceZone.ResonanceZone` objects (and ``ring_order``), which is
why they live here and not in :mod:`tanglepack.topology.plotting`: topology
never imports loom. The cartoon itself, and the enclosures it routes nested
arcs by (:func:`~tanglepack.topology.plotting.cartoon_enclosures`), are
topology's.

A circle cartoon passes all of it explicitly::

    interior, parents, order = cartoon_zones(session.resonance_zones.values(), session.workbench)
    extra = {"enclosures": plotting.cartoon_enclosures(
        dual.trellis, dual.minimal.kept_bridge_ids, parents, session.fixed_points)} if parents else {}
    plotting.plot_dual_graph_cartoon(dual, dynamics, shape="circle", interior_side=interior,
                                     parents=parents, arc_order=order, **extra)
"""

from __future__ import annotations

import logging
from typing import Iterable, Optional, TYPE_CHECKING

import numpy as np

from .ResonanceZone import ResonanceZone, ring_order

if TYPE_CHECKING:
    from ..numerics.TangleWorkbench import TangleWorkbench

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def cartoon_zones(
    zones: Iterable[ResonanceZone], workbench: "TangleWorkbench"
) -> tuple[dict, dict, dict]:
    """
    What the circular dual-graph cartoon needs from the resonance zones.

    For every stable branch bounding a zone, which side of it faces the
    zone's interior: a few points along the branch's boundary arc (anchor
    to cut point) are pushed a short distance to the LEFT of the stable
    dynamical direction (toward the anchor) and tested against the zone's
    captured boundary polygon; the majority decides. And which zone nests
    in which: a fixed point's circle goes inside the smallest OTHER zone
    containing its first orbit point. And how each zone's stable arcs go
    around it: the ring order of its boundary (pip_j, stable arc to z_j,
    unstable side from z_j to pip_{j+1}, as ``ResonanceZone`` stitches
    it), its rotational direction (the sign of the boundary polygon's
    area) and the angle of every orbit point about the orbit's centroid.

    Args:
        zones: The resonance zones (a session's ``resonance_zones.values()``).
        workbench: The workbench whose manifolds bound them.

    Returns:
        ``(interior_side, parents, arc_order)``: ``{branch_key: "left" |
        "right"}`` for every branch bounding a zone, ``{id(fixed_point):
        id(parent fixed point)}`` for the nested fixed points, and
        ``{id(fixed_point): ZoneArcOrder}``. All empty without zones (the
        cartoon then draws each branch's right side inside its circle, in
        ``branch_cycle`` order, and puts the circles side by side).
    """
    from ..numerics.geometry import arc_polyline

    tol = workbench.intersection_registry.cdist_tol
    interior: dict = {}
    zones = sorted(zones, key=lambda zone: zone.area)
    for zone in zones:
        vertices = zone.boundary_vertices
        if vertices is None or len(vertices) < 3:
            continue
        scale = float(np.linalg.norm(vertices.max(axis=0) - vertices.min(axis=0)))
        eps = 1e-3 * max(scale, 1e-9)
        for cut in zone.boundary_intersections:
            key = cut.manifold_b_key
            if key is None or key in interior or key not in workbench.manifolds:
                continue
            points = np.asarray(
                arc_polyline(
                    workbench.manifolds[key], 0.0, cut.stable_cdist,
                    stability="stable", tol=tol,
                ),
                dtype=float,
            )
            votes = {"left": 0, "right": 0}
            count = len(points)
            for i in sorted({int(count * f) for f in (0.3, 0.45, 0.6, 0.75)}):
                if not 1 <= i < count:
                    continue
                toward_anchor = points[i - 1] - points[i]
                norm = float(np.linalg.norm(toward_anchor))
                if norm == 0.0:
                    continue
                left = np.array([-toward_anchor[1], toward_anchor[0]]) / norm
                in_left = zone.contains_point(points[i] + eps * left)
                in_right = zone.contains_point(points[i] - eps * left)
                if in_left != in_right:
                    votes["left" if in_left else "right"] += 1
            if votes["left"] == votes["right"]:
                logger.info(
                    "cartoon zones: no side of stable branch %s is clearly inside "
                    "its zone; the cartoon default is used",
                    key[1:],
                )
                continue
            interior[key] = "left" if votes["left"] > votes["right"] else "right"
    from ..numerics.geometry import signed_polygon_area
    from ..topology.plotting import ZoneArcOrder

    arc_order: dict = {}
    for zone in zones:
        order = ring_order(zone.boundary_intersections)
        vertices = zone.boundary_vertices
        if not order or vertices is None or len(vertices) < 3:
            continue
        fixed_point = zone.fixed_point
        points = [np.asarray(c, dtype=float) for c in fixed_point.coordinates]
        centroid = np.mean(points, axis=0)
        angles = {}
        if len(order) > 1:
            for key in order:
                offset = points[key[2]] - centroid
                angles[key] = float(np.arctan2(offset[1], offset[0]))
        arc_order.setdefault(
            id(fixed_point),
            ZoneArcOrder(
                order=order,
                clockwise=signed_polygon_area(vertices) < 0,
                anchor_angles=angles,
            ),
        )
    parents: dict = {}
    for zone in zones:
        fixed_point = zone.fixed_point
        point = np.asarray(fixed_point.coordinates[0], dtype=float)
        containers = [
            other
            for other in zones
            if other.fixed_point is not fixed_point
            and other.area > zone.area
            and other.contains_point(point)
        ]
        if containers:
            parent = min(containers, key=lambda other: other.area).fixed_point
            parents[id(fixed_point)] = id(parent)
    return interior, parents, arc_order


def plot_resonance_zones(
    zones: Iterable[ResonanceZone],
    ax=None,
    colors: Optional[list] = None,
    alpha: float = 0.3,
    **fill_kwargs,
) -> list:
    """
    Shade resonance zones as filled regions, each a distinct color.

    Fills are translucent (low ``alpha``) so the red stable / blue unstable
    manifold lines and any nested inner zone remain visible on top. Zones are
    drawn in the order given (a session's ``resonance_zones.values()`` is
    insertion order).

    Args:
        zones: The zones to shade; each carries the boundary polygon captured
            when it was defined (``boundary_vertices``).
        ax: Optional matplotlib Axes. Defaults to the current axes (plt).
        colors: Optional list of fill colors cycled across zones. Defaults to a
            distinct qualitative palette.
        alpha: Fill transparency (default 0.3).
        **fill_kwargs: Forwarded to ``fill`` (e.g. ``zorder``, ``label``).

    Returns:
        List of the matplotlib polygon handles drawn (one per zone).
    """
    import matplotlib.pyplot as plt

    zones = list(zones)
    if not zones:
        logger.info("No resonance zones to plot; call resonance_zone() first.")
        return []

    palette = colors or [
        "#8c6bb1",  # purple
        "#41ab5d",  # green
        "#fd8d3c",  # orange
        "#6baed6",  # light blue
        "#df65b0",  # magenta
        "#a6761d",  # ochre
    ]
    target = ax if ax is not None else plt
    handles = []
    for i, rz in enumerate(zones):
        poly = rz.boundary_vertices
        if poly.size == 0:
            continue
        color = palette[i % len(palette)]
        handles.append(
            target.fill(
                poly[:, 0], poly[:, 1], color=color, alpha=alpha, **fill_kwargs
            )
        )
    return handles
