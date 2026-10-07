"""Map and session helpers shared by the stress experiments.

Reuses :mod:`tanglepack.examples.henon` for the map itself; adds the analytic
period-1 Hénon fixed points for ANY ``(k, b)`` (the examples module only
records seeds for two parameter pairs), standard session recipes, and cheap
size snapshots of a workbench.
"""

from __future__ import annotations

import math
from typing import Any, Optional

import numpy as np

import tanglepack
from tanglepack.examples.henon import henon_jacobian, henon_map, henon_map_inverse

#: Orientation hints used by every period-1 recipe in tests/cases.py.
PERIOD1_ORIENTATION = {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}


def henon_fixed_points(k: float, b: float = 1.0) -> dict[str, Optional[list[float]]]:
    """Analytic period-1 fixed points of ``(x, y) -> (y - k + x^2, -b x)``.

    ``x = -b x - k + x^2`` gives ``x^2 - (1 + b) x - k = 0``; ``y = -b x``.

    Returns:
        ``{"saddle": [x+, y+], "other": [x-, y-]}`` (``x+`` the larger root,
        the regular saddle for ``b = 1``; the other root is the inversion
        saddle for large ``k`` and elliptic for small ``k``), or Nones when
        the discriminant is negative (no fixed points).
    """
    disc = (1 + b) ** 2 + 4 * k
    if disc < 0:
        return {"saddle": None, "other": None}
    xp = ((1 + b) + math.sqrt(disc)) / 2
    xm = ((1 + b) - math.sqrt(disc)) / 2
    return {"saddle": [xp, -b * xp], "other": [xm, -b * xm]}


def henon_session(
    k: float,
    b: float = 1.0,
    *,
    jacobian: bool = True,
    area_cutoff: Optional[float] = None,
    gpu: bool = False,
    min_batch_points: int = 64,
) -> tanglepack.TangleSession:
    """A TangleSession on the Hénon map, optionally GPU-wrapped."""
    session = tanglepack.TangleSession(
        henon_map(k, b),
        henon_map_inverse(k, b),
        henon_jacobian(k, b) if jacobian else None,
    )
    if area_cutoff is not None:
        session.workbench._man_machine.area_cutoff = area_cutoff
    if gpu:
        tanglepack.enable_gpu(session, min_batch_points=min_batch_points)
    return session


def saddle_session(
    k: float,
    b: float = 1.0,
    *,
    which: str = "saddle",
    area_cutoff: Optional[float] = None,
    gpu: bool = False,
    orient: bool = True,
) -> tuple[tanglepack.TangleSession, Any]:
    """Session + constructed, oriented, seeded period-1 saddle."""
    session = henon_session(k, b, area_cutoff=area_cutoff, gpu=gpu)
    seed = henon_fixed_points(k, b)[which]
    if seed is None:
        raise ValueError(f"no real period-1 fixed point for k={k}, b={b}")
    fp = session.construct_fixed_point(seed)
    if orient:
        session.orient_eigenvectors(fp, PERIOD1_ORIENTATION)
    session.initialize_both_manifolds(fp)
    return session, fp


def manifold_sizes(session_or_workbench: Any) -> dict[str, int]:
    """Point count per registered manifold, keyed ``"<stability>/<orbit>.<branch>"``."""
    wb = getattr(session_or_workbench, "workbench", session_or_workbench)
    sizes: dict[str, int] = {}
    for (fp, stability, orbit, branch), manifold in wb.manifolds.items():
        label = getattr(fp, "label", None) or "fp"
        key = f"{label}:{stability}/{orbit}.{branch}"
        sizes[key] = int(len(manifold.get_point_array()))
    return sizes


def total_points(session_or_workbench: Any, stability: Optional[str] = None) -> int:
    """Total manifold points, optionally restricted to one stability."""
    return sum(
        n
        for key, n in manifold_sizes(session_or_workbench).items()
        if stability is None or f":{stability}/" in key
    )


def max_abs_coord(session_or_workbench: Any) -> float:
    """Largest |coordinate| over every manifold (escape-to-infinity probe)."""
    wb = getattr(session_or_workbench, "workbench", session_or_workbench)
    out = 0.0
    for manifold in wb.manifolds.values():
        pts = manifold.get_point_array()
        if len(pts):
            out = max(out, float(np.nanmax(np.abs(pts))))
    return out


# --------------------------------------------------------------------------- #
# Geometry probes shared by e04 / e05 (additive)
# --------------------------------------------------------------------------- #
def distance_to_polyline(
    queries: np.ndarray, polyline: np.ndarray, k: int = 8
) -> tuple[np.ndarray, np.ndarray]:
    """Distance from each query to a polyline, via a KD-tree on its vertices.

    The ``k`` nearest vertices are candidates; the two segments meeting at
    each are checked exactly. Exact unless the nearest segment is very long
    compared with the vertex spacing elsewhere (the refinement keeps segments
    short, so ``k = 8`` is ample).

    Returns:
        ``(distance, beyond)``: ``beyond`` flags a query whose nearest point is
        the polyline's last vertex (a point mapped past the tail).
    """
    from scipy.spatial import cKDTree

    queries = np.asarray(queries, dtype=float)
    poly = np.asarray(polyline, dtype=float)
    n = len(poly)
    k = min(k, n)
    _, idx = cKDTree(poly).query(queries, k=k)
    idx = idx.reshape(len(queries), k)
    seg = np.concatenate([idx - 1, idx], axis=1).clip(0, n - 2)  # segment i = (i, i+1)
    a, b = poly[seg], poly[seg + 1]
    ab = b - a
    q = queries[:, None, :]
    t = np.einsum("ijk,ijk->ij", q - a, ab) / np.maximum(np.einsum("ijk,ijk->ij", ab, ab), 1e-300)
    t = t.clip(0.0, 1.0)
    d = np.linalg.norm(q - (a + t[..., None] * ab), axis=2)
    best = d.argmin(axis=1)
    rows = np.arange(len(queries))
    beyond = (seg[rows, best] == n - 2) & (t[rows, best] >= 1.0)
    return d[rows, best], beyond


def invariance_error(
    manifold: Any,
    step: Any,
    stability: str,
    samples: int = 2000,
    seed: int = 0,
    target: Any = None,
) -> dict[str, Any]:
    """How far one map step moves manifold points off the computed curve.

    The manifold is built by mapping its own vertices, so a vertex usually
    maps EXACTLY onto a vertex and tells nothing. The probe therefore samples
    a point at a random fraction ``t`` of a segment ``(v_i, v_{i+1})`` -- the
    polyline's interpolation, i.e. what every crossing is computed on -- and
    measures the distance from ``step(p)`` to the polyline. A segment is
    eligible when both its vertices map onto vertices (to 1e-9 relative),
    i.e. its image is still inside the computed curve; that also excludes the
    newest growth step and the root.

    Args:
        manifold: A ``BaseManifold``.
        step: The map for the unstable side, the inverse map for the stable
            side; takes and returns ``(2, N)``.
        stability: ``"unstable"`` or ``"stable"``.
        samples: Eligible segments sampled (all when fewer).
        seed: Sampling seed.
        target: The manifold ``step`` carries this one onto (the next branch
            of a period-k orbit, ``FixedPoint.advance_key``); defaults to
            ``manifold`` itself.

    Returns:
        ``{"err", "cdist", "beyond", "vertex_err"}`` arrays over the sample and
        ``"eligible"`` (eligible segment count, after sampling),
        ``"vertices"`` and ``"on_curve"`` (vertices whose image is a vertex).
        ``cdist`` interpolates the
        node cdists; ``vertex_err`` is the distance from ``step(v_i)`` to the
        nearest vertex (0 when the vertex image is itself a vertex).
    """
    from scipy.spatial import cKDTree

    nodes = manifold.get_point_array(return_nodes=True)
    poly = np.asarray([n.get_point() for n in nodes], dtype=float).reshape(-1, 2)
    cdist = np.array([n.get_cdist(stability) for n in nodes], dtype=float)
    image_poly = poly if target is None or target is manifold else (
        np.asarray(target.get_point_array(), dtype=float).reshape(-1, 2)
    )
    vertex_images = np.asarray(step(poly.T), dtype=float).T
    vertex_err, _ = cKDTree(image_poly).query(vertex_images)
    scale = max(1.0, float(np.abs(image_poly).max()))
    on_curve = vertex_err <= 1e-9 * scale
    eligible = np.flatnonzero(on_curve[:-1] & on_curve[1:])
    eligible = eligible[eligible > 0]
    rng = np.random.default_rng(seed)
    if len(eligible) > samples:
        eligible = np.sort(rng.choice(eligible, samples, replace=False))
    out: dict[str, Any] = {"eligible": int(len(eligible)), "vertices": len(poly),
                           "on_curve": int(on_curve.sum())}
    if not len(eligible):
        empty = np.zeros(0)
        return {**out, "err": empty, "cdist": empty, "beyond": empty.astype(bool),
                "vertex_err": empty}
    t = rng.uniform(0.1, 0.9, len(eligible))[:, None]
    sample = (1 - t) * poly[eligible] + t * poly[eligible + 1]
    images = np.asarray(step(sample.T), dtype=float).T
    err, beyond = distance_to_polyline(images, image_poly)
    sample_cdist = (1 - t[:, 0]) * cdist[eligible] + t[:, 0] * cdist[eligible + 1]
    return {**out, "err": err, "cdist": sample_cdist, "beyond": beyond,
            "vertex_err": vertex_err[eligible]}


def error_summary(values: np.ndarray, prefix: str = "") -> dict[str, float]:
    """Count, max, mean and the 50/90/99 percentiles of a non-negative sample."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return {f"{prefix}n": 0}
    p50, p90, p99 = np.percentile(values, [50, 90, 99])
    return {f"{prefix}n": int(len(values)), f"{prefix}max": float(values.max()),
            f"{prefix}mean": float(values.mean()), f"{prefix}p50": float(p50),
            f"{prefix}p90": float(p90), f"{prefix}p99": float(p99)}


def curve_arrays(session_or_workbench: Any) -> dict[Any, tuple[np.ndarray, np.ndarray]]:
    """``{manifold key: (cdist, points)}`` per manifold, ascending in cdist.

    Read once so arcs can be sliced with ``searchsorted`` instead of walking
    the linked list per arc (``geometry.arc_polyline`` is O(N) per call).
    """
    wb = getattr(session_or_workbench, "workbench", session_or_workbench)
    curves = {}
    for key, manifold in wb.manifolds.items():
        nodes = manifold.get_point_array(return_nodes=True)
        cdist = np.array([n.get_cdist(key[1]) for n in nodes], dtype=float)
        pts = np.asarray([n.get_point() for n in nodes], dtype=float).reshape(-1, 2)
        order = np.argsort(cdist, kind="stable")
        curves[key] = (cdist[order], pts[order])
    return curves


def curve_arc(curves: dict, key: Any, c_from: float, c_to: float) -> Optional[np.ndarray]:
    """The nodes of ``curves[key]`` strictly between two cdists, from ``c_from`` to ``c_to``.

    Returns None when the curve does not reach ``max(c_from, c_to)``.
    """
    if key not in curves:
        return None
    cdist, pts = curves[key]
    lo, hi = sorted((c_from, c_to))
    if not len(cdist) or cdist[-1] < hi:
        return None
    arc = pts[np.searchsorted(cdist, lo, "right"):np.searchsorted(cdist, hi, "left")]
    return arc if c_from <= c_to else arc[::-1]


def unstable_arc(session_or_workbench: Any, curves: dict, a: int, b: int) -> Optional[np.ndarray]:
    """The unstable arc from crossing ``a`` to crossing ``b`` (``a`` nearer the fixed point).

    Walks the chain of registered bridges ``(a, x1), (x1, x2), ... (xm, b)``
    first (it carries blast geometry beyond the grown manifold), and falls
    back to slicing the unstable manifold.
    """
    from tanglepack.numerics.geometry import oriented_bridge_polyline

    wb = getattr(session_or_workbench, "workbench", session_or_workbench)
    reg = wb.intersection_registry
    target = reg[b].unstable_cdist
    pieces, cur = [], a
    for _ in range(100_000):
        if cur == b:
            return np.vstack(pieces) if pieces else np.empty((0, 2))
        ahead = [bid for bid in wb.bridges_at(cur) if bid[0] == cur]
        if not ahead or reg[ahead[0][1]].unstable_cdist > target + 1e-12:
            break
        bid = ahead[0]
        poly = oriented_bridge_polyline(
            wb.bridge(bid), reg[bid[0]].unstable_cdist, reg[bid[1]].unstable_cdist
        )
        if poly is None:
            break
        pieces.append(poly)
        cur = bid[1]
    return curve_arc(curves, reg[a].manifold_a_key, reg[a].unstable_cdist, target)


def lobe_polygon(session_or_workbench: Any, curves: dict, a: int, b: int) -> Optional[np.ndarray]:
    """The closed lobe bounded by the unstable arc ``a -> b`` and the stable arc ``b -> a``.

    Returns None when the two crossings are not on one unstable branch and one
    stable branch, or when either arc is not computed.
    """
    wb = getattr(session_or_workbench, "workbench", session_or_workbench)
    ia, ib = wb.intersection_registry[a], wb.intersection_registry[b]
    if ia.manifold_a_key != ib.manifold_a_key or ia.manifold_b_key != ib.manifold_b_key:
        return None
    if ia.unstable_cdist > ib.unstable_cdist:
        a, b, ia, ib = b, a, ib, ia
    u = unstable_arc(wb, curves, a, b)
    s = curve_arc(curves, ia.manifold_b_key, ib.stable_cdist, ia.stable_cdist)
    if u is None or s is None:
        return None
    return np.vstack([ia.get_point(), u, ib.get_point(), s])
