"""Crossing detection against an independent brute force.

After ``compute_intersections`` (the rtree + Tangle path), every proper
crossing of every unstable polyline with every stable polyline is recomputed
here with exact floating-point orientation tests (no epsilon, no rtree), in
chunked vectorised O(N*M), and matched against the registry:

* ``matched`` -- a registry crossing on the same (unstable, stable) branch pair
  within ``1e-9 * max(1, |coord|)``;
* ``missed`` -- brute-force crossings the registry does not hold, each
  classified as ``near_parallel`` (``|sin| <= 1e-12``, the library's discard
  threshold), ``zigzag`` (another brute crossing of the same pair within
  ``1e-9`` relative unstable cdist: the noise runs ``_collapse_noise_crossings``
  folds), ``collision`` (a held crossing within the registry's absolute
  ``cdist_tol = 1e-6`` on both cdists: ``IntersectionRegistry.add`` merged it)
  or ``unexplained``;
* ``extra`` -- registry crossings with no brute partner: ``anchor`` (the
  synthetic ``(0, 0)`` crossing at the periodic point), ``touch`` (the brute
  force sees an exact zero orientation there, a shared or touching vertex the
  library's epsilon admits) or ``unexplained``.

Brute-force cdists are interpolated linearly along the segment from the
nodes' own cdists, the rule ``Tangle._cdist_between`` uses. As a side check the
unstable polylines are also tested against each other (non-adjacent segments
only): a proper u x u crossing is a polygon that breaks the fundamental
invariant at the discretisation level (``uu_crossings``; skipped above
``UU_MAX_SEGMENTS``).
"""

from __future__ import annotations

import sys
import time
from itertools import combinations

import numpy as np

import harness
from maps import saddle_session

sys.path.insert(0, str(harness.REPO / "tests"))

EXPERIMENT = "e09_intersections"
TIMEOUT = 600
WORKERS = 4
RSS_CAP_GB = 16

CHUNK_PAIRS = 4_000_000
UU_MAX_SEGMENTS = 40_000
NEAR_PARALLEL = 1e-12
MATCH_RTOL = 1e-9
CDIST_TOL = 1e-6


def configs(quick: bool) -> list[dict]:
    if quick:
        return [
            {"case": "k10", "unstable_steps": 7},
            {"case": "k28", "unstable_steps": 6, "area_cutoff": 1e-5},
            {"case": "p3", "unstable_steps": 9},
        ]
    return (
        [{"case": "k10", "unstable_steps": n} for n in (7, 8, 9, 10)]
        + [
            {"case": "k28", "unstable_steps": n, "area_cutoff": cutoff}
            for cutoff in (1e-5, 1e-6, 1e-7)
            for n in range(6, 11)
        ]
        + [{"case": "k28", "unstable_steps": n, "area_cutoff": 1e-7} for n in (11, 12)]
        + [{"case": "p3", "unstable_steps": n} for n in range(9, 14)]
    )


def build(config: dict):
    """Session and fixed point grown as the matching tests/cases.py recipe."""
    from cases import build_k10
    from tanglepack.examples import henon_cases

    steps = config["unstable_steps"]
    if config["case"] == "k10":
        case = build_k10(unstable_steps=steps, through="grown")
        return case.session, case.fixed_point
    if config["case"] == "k28":
        session, fp = saddle_session(2.8, 1, area_cutoff=config["area_cutoff"])
        session.grow_n_times(fp, "unstable", num_iterations=steps)
        session.grow_until_turnaround(fp, "stable")
        return session, fp
    session = henon_cases._session()
    return session, henon_cases._grow_period3(session, steps, 9, 1e-7)


# --------------------------------------------------------------------------- #
# Brute force
# --------------------------------------------------------------------------- #
def cross(ax, ay, bx, by):
    return ax * by - ay * bx


def segment_crossings(P: np.ndarray, Q: np.ndarray, *, skip_adjacent: bool = False):
    """Every proper crossing and every exact-zero touch of polylines P and Q.

    Returns ``(proper, touch)``: ``proper`` rows are ``(i, j, t, s, sin)`` for
    segment ``i`` of P at fraction ``t`` crossing segment ``j`` of Q at ``s``;
    ``touch`` rows are ``(i, j)`` pairs whose orientation test hit an exact zero
    while straddling otherwise (endpoint contact).
    """
    p0, d = P[:-1], np.diff(P, axis=0)
    q0, e = Q[:-1], np.diff(Q, axis=0)
    q_lo, q_hi = np.minimum(Q[:-1], Q[1:]), np.maximum(Q[:-1], Q[1:])
    chunk = max(1, CHUNK_PAIRS // max(len(q0), 1))
    proper, touch = [], []
    for start in range(0, len(p0), chunk):
        a, da = p0[start:start + chunk, None, :], d[start:start + chunk, None, :]
        lo = np.minimum(a, a + da)
        hi = np.maximum(a, a + da)
        boxes = np.all((lo <= q_hi[None]) & (hi >= q_lo[None]), axis=2)
        i, j = np.nonzero(boxes)
        if skip_adjacent:
            keep = np.abs(i + start - j) > 1
            i, j = i[keep], j[keep]
        if not len(i):
            continue
        A, D, C, E = p0[i + start], d[i + start], q0[j], e[j]
        o1 = cross(D[:, 0], D[:, 1], C[:, 0] - A[:, 0], C[:, 1] - A[:, 1])
        o2 = cross(D[:, 0], D[:, 1], C[:, 0] + E[:, 0] - A[:, 0], C[:, 1] + E[:, 1] - A[:, 1])
        o3 = cross(E[:, 0], E[:, 1], A[:, 0] - C[:, 0], A[:, 1] - C[:, 1])
        o4 = cross(E[:, 0], E[:, 1], A[:, 0] + D[:, 0] - C[:, 0], A[:, 1] + D[:, 1] - C[:, 1])
        is_proper = (o1 * o2 < 0) & (o3 * o4 < 0)
        is_touch = ~is_proper & (o1 * o2 <= 0) & (o3 * o4 <= 0) & ((o1 * o2 == 0) | (o3 * o4 == 0))
        denom = cross(D[:, 0], D[:, 1], E[:, 0], E[:, 1])
        norm = np.linalg.norm(D, axis=1) * np.linalg.norm(E, axis=1)
        with np.errstate(all="ignore"):
            t = o3 / (o3 - o4)
            s = o1 / (o1 - o2)
            sin = np.abs(denom) / norm
        rows = np.column_stack([i + start, j, t, s, sin])[is_proper]
        proper.append(rows)
        touch.append(np.column_stack([i + start, j])[is_touch])
    proper = np.vstack(proper) if proper else np.empty((0, 5))
    touch = np.vstack(touch) if touch else np.empty((0, 2))
    return proper, touch


def at(values: np.ndarray, index: np.ndarray, fraction: np.ndarray) -> np.ndarray:
    """Linear interpolation of per-node values at segment ``index`` + ``fraction``."""
    index = index.astype(int)
    return values[index] + fraction[:, None] * (values[index + 1] - values[index])


def brute_force(session) -> dict:
    """Brute-force crossings per (unstable key, stable key) pair."""
    wb = session.workbench
    polylines = {
        key: (m.get_point_array(), np.asarray(m.get_cdist_array()).reshape(-1, 1))
        for key, m in wb.manifolds.items()
    }
    unstable = [k for k in polylines if k[1] == "unstable"]
    stable = [k for k in polylines if k[1] == "stable"]
    found = {}
    for u_key in unstable:
        for s_key in stable:
            (U, u_cd), (S, s_cd) = polylines[u_key], polylines[s_key]
            proper, touch = segment_crossings(U, S)
            ui, si, t, s, sin = proper.T
            touch_points = U[touch[:, 0].astype(int)] if len(touch) else np.empty((0, 2))
            found[(u_key, s_key)] = {
                "coords": at(U, ui, t),
                "u_cdist": at(u_cd, ui, t)[:, 0],
                "s_cdist": at(s_cd, si, s)[:, 0],
                "sin": sin,
                "touch": touch_points,
            }
    return found


def uu_crossings(session) -> int | None:
    """Proper crossings between unstable polylines (non-adjacent segments only)."""
    lines = [m.get_point_array() for k, m in session.workbench.manifolds.items() if k[1] == "unstable"]
    if sum(len(line) for line in lines) > UU_MAX_SEGMENTS:
        return None
    count = sum(len(segment_crossings(line, line, skip_adjacent=True)[0]) // 2 for line in lines)
    count += sum(len(segment_crossings(a, b)[0]) for a, b in combinations(lines, 2))
    return count


# --------------------------------------------------------------------------- #
# Matching
# --------------------------------------------------------------------------- #
def close(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Boolean matrix: rows of a within MATCH_RTOL relative of rows of b."""
    tol = MATCH_RTOL * np.maximum(1.0, np.abs(a).max(axis=1))[:, None]
    return np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2) <= tol


def classify(session, found: dict, rec: harness.Recorder) -> None:
    registry = session.workbench.intersection_registry
    held: dict = {}
    anchors = 0
    for _id, ix in registry:
        if ix.label == "anchor":
            anchors += 1
            continue
        held.setdefault((ix.manifold_a_key, ix.manifold_b_key), []).append(ix)

    counts = {"matched": 0, "near_parallel": 0, "zigzag": 0, "collision": 0,
              "missed_unexplained": 0, "extra_touch": 0, "extra_unexplained": 0}
    extra_distance, missed_sin = [], []
    for pair in set(found) | set(held):
        brute = found.get(pair)
        lib = held.get(pair, [])
        n_brute = 0 if brute is None else len(brute["coords"])
        lib_xy = np.array([ix.coords for ix in lib]).reshape(-1, 2)
        lib_cd = np.array([(ix.unstable_cdist, ix.stable_cdist) for ix in lib]).reshape(-1, 2)
        hits = close(brute["coords"], lib_xy) if n_brute and len(lib) else np.zeros((n_brute, len(lib)), bool)
        matched_brute = hits.any(axis=1)
        counts["matched"] += int(matched_brute.sum())

        for k in np.nonzero(~matched_brute)[0]:
            u, s = brute["u_cdist"][k], brute["s_cdist"][k]
            missed_sin.append(float(brute["sin"][k]))
            others = np.delete(brute["u_cdist"], k)
            if brute["sin"][k] <= NEAR_PARALLEL:
                counts["near_parallel"] += 1
            elif np.any(np.abs(others - u) <= MATCH_RTOL * max(1.0, abs(u))):
                counts["zigzag"] += 1
            elif len(lib) and np.any((np.abs(lib_cd[:, 0] - u) < CDIST_TOL) & (np.abs(lib_cd[:, 1] - s) < CDIST_TOL)):
                counts["collision"] += 1
            else:
                counts["missed_unexplained"] += 1

        for k in np.nonzero(~hits.any(axis=0))[0]:
            point = lib_xy[k]
            touch = brute["touch"] if brute is not None else np.empty((0, 2))
            pool = np.vstack([brute["coords"], touch]) if brute is not None else touch
            distance = float(np.min(np.linalg.norm(pool - point, axis=1))) if len(pool) else None
            extra_distance.append(distance)
            near_touch = len(touch) and np.min(np.linalg.norm(touch - point, axis=1)) <= 1e-6
            counts["extra_touch" if near_touch else "extra_unexplained"] += 1

    rec.metric("library_crossings", len(registry))
    rec.metric("anchors", anchors)
    rec.metric("brute_crossings", sum(len(v["coords"]) for v in found.values()))
    rec.metric("brute_touches", sum(len(v["touch"]) for v in found.values()))
    for name, value in counts.items():
        rec.metric(name, value)
    rec.metric("missed_sin", missed_sin)
    rec.metric("extra_nearest_brute_distance", extra_distance)
    sines = np.concatenate([v["sin"] for v in found.values()]) if found else np.empty(0)
    rec.metric("brute_min_sin", float(sines.min()) if len(sines) else None)
    rec.metric("brute_n_sin_below_1e-9", int((sines < 1e-9).sum()))


def run(config: dict, rec: harness.Recorder) -> None:
    with rec.stage("build"):
        session, fp = build(config)
    manifolds = session.workbench.manifolds
    rec.metric("unstable_segments", sum(len(m.get_point_array()) - 1 for k, m in manifolds.items() if k[1] == "unstable"))
    rec.metric("stable_segments", sum(len(m.get_point_array()) - 1 for k, m in manifolds.items() if k[1] == "stable"))
    rec.checkpoint()

    t0 = time.perf_counter()
    session.compute_intersections([fp], infer_iterates=False)
    rec.metric("library_s", time.perf_counter() - t0)
    t0 = time.perf_counter()
    found = brute_force(session)
    rec.metric("brute_s", time.perf_counter() - t0)
    rec.checkpoint()

    classify(session, found, rec)
    rec.checkpoint()
    t0 = time.perf_counter()
    rec.metric("uu_crossings", uu_crossings(session))
    rec.metric("uu_s", time.perf_counter() - t0)


if __name__ == "__main__":
    harness.main(globals())
