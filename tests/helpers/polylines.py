"""
Synthetic polylines for the numerics guards: hand-placed manifolds on a real
workbench.

Nothing here calls a map. :func:`polyline_manifold` builds a real
``BaseManifold`` through real ``Point`` nodes; :func:`polyline_workbench`
registers such curves on a real ``TangleWorkbench`` (identity maps), so
``compute_intersections``, the resolution ratio and the resolution report run
unchanged on geometry chosen by the test.
"""

from __future__ import annotations

from typing import Mapping, Optional, Sequence

from tanglepack import BaseManifold, BranchPoint, FixedPoint, Point, TangleWorkbench
from tanglepack.numerics import ManifoldMachine


def polyline_manifold(
    points: Sequence[Sequence[float]],
    key: tuple,
    *,
    root: Optional[BranchPoint] = None,
) -> BaseManifold:
    """
    A manifold through ``points``, cdist = the index, keyed by ``key``.

    The nodes are put on one curve in the order given, through the machine's
    own constructor (``ManifoldMachine._manifold_from_points``), and hung on
    ``root`` when there is one, as the initializer hangs a seed segment.

    Args:
        points: The polyline's vertices in order (after ``root``, if given).
        key: The manifold key ``(fixed_point, stability, orbit, branch)``.
        root: A shared periodic point the curve leaves from, on branch
            ``key[3]``; several curves may share one.

    Returns:
        The manifold, root at ``root`` (or the first point) and tail at the
        last point.
    """
    fixed_point, stability, _orbit, branch = key
    offset = 0 if root is None else 1
    nodes = [Point(x, y, float(i + offset)) for i, (x, y) in enumerate(points)]
    ManifoldMachine._manifold_from_points(nodes, stability, 1.0, fixed_point)
    if root is not None:
        hang = root.insert_point_forward if stability == "unstable" else root.insert_point_backward
        hang(nodes[0], branch)
    return BaseManifold(
        nodes[0] if root is None else root,
        stability,
        1.0,
        fixed_point=fixed_point,
        tail=nodes[-1],
        branch_index=branch,
        manifold_key=key,
    )


def polyline_workbench(
    fixed_point: FixedPoint,
    curves: Mapping[tuple, Sequence[Sequence[float]]],
    *,
    resolution_ratio_threshold: Optional[float] = 3.0,
) -> TangleWorkbench:
    """
    A workbench holding hand-placed curves of one fixed point.

    Args:
        fixed_point: The curves' fixed point (see ``helpers.fakes.bare_fixed_point``).
        curves: Manifold key -> polyline vertices.
        resolution_ratio_threshold: Passed to the workbench.

    Returns:
        The workbench, nothing computed yet.
    """
    identity = lambda point: point  # noqa: E731 - the map is never called
    workbench = TangleWorkbench(
        identity, identity, resolution_ratio_threshold=resolution_ratio_threshold
    )
    workbench.fixed_points.append(fixed_point)
    for key, points in curves.items():
        workbench.register_manifold(key, polyline_manifold(points, key))
    return workbench
