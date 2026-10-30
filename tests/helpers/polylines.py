"""
Synthetic polylines for the numerics guards: hand-placed manifolds on a real
workbench.

Nothing here calls a map. :func:`polyline_manifold` builds a real
``BaseManifold`` over rows of the one test store (:data:`TEST_STORE`); :func:`polyline_workbench`
registers such curves on a real ``TangleWorkbench`` (identity maps), so
``compute_intersections``, the resolution ratio and the resolution report run
unchanged on geometry chosen by the test.
"""

from __future__ import annotations

from typing import Mapping, Optional, Sequence

import numpy as np

from tanglepack import BaseManifold, BranchPoint, FixedPoint, TangleWorkbench
from tanglepack.numerics.PointStore import PointStore

#: The one store of every hand-built test curve (a Tangle asserts one store).
TEST_STORE = PointStore()


def store_rows(
    points: Sequence[Sequence[float]],
    cdists: Optional[Sequence[float]] = None,
    *,
    store: Optional[PointStore] = None,
):
    """Append hand-placed points (cdist = the index unless given) to the test
    store, on no chain yet: their rows, in order."""
    store = TEST_STORE if store is None else store
    cdists = range(len(points)) if cdists is None else cdists
    return store.append(points, [float(c) for c in cdists])


def polyline_manifold(
    points: Sequence[Sequence[float]],
    key: tuple,
    *,
    root: Optional[BranchPoint] = None,
    store: Optional[PointStore] = None,
) -> BaseManifold:
    """
    A manifold through ``points``, cdist = the index, keyed by ``key``.

    The rows are put on one chain in the order given and hung on
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
    store = TEST_STORE if store is None else store
    rows = store_rows(points, [i + offset for i in range(len(points))], store=store)
    chain = store.new_chain(rows, 1.0)
    if root is not None:
        root.attach(chain, stability, branch)
    return BaseManifold(
        int(rows[0]) if root is None else root,
        stability,
        1.0,
        fixed_point=fixed_point,
        tail=int(rows[-1]),
        branch_index=branch,
        manifold_key=key,
        store=store,
    )


def polyline_workbench(
    fixed_point: FixedPoint,
    curves: Mapping[tuple, Sequence[Sequence[float]]],
    *,
    resolution_ratio_threshold: Optional[float] = 3.0,
    anchor_resolution_factor: Optional[float] = None,
) -> TangleWorkbench:
    """
    A workbench holding hand-placed curves of one fixed point.

    Args:
        fixed_point: The curves' fixed point (see ``helpers.fakes.bare_fixed_point``).
        curves: Manifold key -> polyline vertices.
        resolution_ratio_threshold: Passed to the workbench.
        anchor_resolution_factor: Passed to the workbench. None by default:
            hand-placed curves put their fixed point anywhere (the origin) and
            their coarse polylines near it, so clause 3 of N is off unless a
            test of it asks for it.

    Returns:
        The workbench, nothing computed yet.
    """
    identity = lambda point: point  # noqa: E731 - the map is never called
    workbench = TangleWorkbench(
        identity,
        identity,
        resolution_ratio_threshold=resolution_ratio_threshold,
        anchor_resolution_factor=anchor_resolution_factor,
    )
    workbench.fixed_points.append(fixed_point)
    for key, points in curves.items():
        workbench.register_manifold(key, polyline_manifold(points, key))
    return workbench


def periodic_point(x: float = 0.0, y: float = 0.0) -> BranchPoint:
    """The orbit point of a hand-built period-1 fixed point at ``(x, y)``: its
    own image, cdist 0.0, every slot empty."""
    fixed_point = FixedPoint(1)
    fixed_point.coordinates = [np.array([x, y], dtype=np.float64)]
    return fixed_point.branch_points[0]
