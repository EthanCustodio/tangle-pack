"""The resolved-crossing oracle (B1 step 0): which detections the polylines resolve.

Usage (from a worktree root, against THAT worktree's library)::

    PYTHONPATH=src:tests MPLBACKEND=Agg env/bin/python tests/_tools/resolution_probe.py \\
        [--laws] [--builders [NAME ...]] [--series [NAME ...]] [--k28x8]

A detection is out of N (A.2, Q13 (a)) by

* clause 1: its resolution ratio is below 3;
* clause 2: a detection of the same ``(u_key, s_key)``, its neighbour in unstable
  cdist order, sits on the same or a node-sharing unstable segment closer along
  the polyline than the LARGER chord deviation of the two. Both are out, and so a
  run is out whole, median included.

The chord deviations are recomputed here, the way the library computes them
(``chord_deviation_batch`` over ``Tangle._outer_neighbours``), by wrapping
``Tangle._resolution_ratio`` before anything is built. The probe never reads a
deviation the library stored, so it stays an oracle. Deviations are recorded
under the crossing's four bracket ends (store rows, or the periodic point's
``BranchPoint``), the key a registered crossing carries since N5a. A cross-check
finds the same pairs by arc length along the whole curve object, with no
adjacency restriction.

Dev Notes:
    Not collected by pytest (no ``test_`` prefix). The printed numbers are
    reported, never pinned (test policy, 2026-10-05).
"""

from __future__ import annotations

import argparse
import sys
import time

import numpy as np

from tanglepack.numerics.BranchPoint import BranchPoint
from tanglepack.numerics.Tangle import Tangle
from tanglepack.numerics.geometry import chord_deviation_batch

THRESHOLD = 3.0

#: ``(bracket ends of the unstable and the stable segment) -> max(h_u, h_s)``.
DEVIATION: dict[tuple, float] = {}


def end(node):
    """A node as a bracket end: its store row, or the periodic point's BranchPoint."""
    return node if isinstance(node, BranchPoint) else node.row

_library_ratio = Tangle._resolution_ratio


def _recording_ratio(u_seg, s_seg, sin_theta):
    """The library's ratio, with the probe's own ``max(h_u, h_s)`` recorded beside it."""
    segments = (u_seg, s_seg)
    p0 = np.array([seg.p0.get_point() for seg in segments], dtype=float)
    p1 = np.array([seg.p0_seg1.get_point() for seg in segments], dtype=float)
    left, right = zip(*(Tangle._outer_neighbours(seg) for seg in segments))
    h = chord_deviation_batch(p0, p1, np.array(left), np.array(right))
    DEVIATION[end(u_seg.p0), end(u_seg.p0_seg1), end(s_seg.p0), end(s_seg.p0_seg1)] = float(h.max())
    return _library_ratio(u_seg, s_seg, sin_theta)


Tangle._resolution_ratio = staticmethod(_recording_ratio)


def deviation(ix) -> float | None:
    """The probe's ``max(h_u, h_s)`` of one detection."""
    return DEVIATION.get((*ix.unstable_segment, *ix.stable_segment))


def separation(ix, jx) -> float | None:
    """Distance along the unstable polyline, through a shared node; None when not adjacent."""
    p, q = ix.unstable_segment, jx.unstable_segment
    a, b = np.asarray(ix.coords, float), np.asarray(jx.coords, float)
    if set(p) == set(q):
        return float(np.hypot(*(a - b)))
    shared = [n for n in p if n in q]
    if len(shared) != 1:
        return None
    node = np.asarray(ix.segment_ends("unstable")[0][p.index(shared[0])], float)
    return float(np.hypot(*(a - node)) + np.hypot(*(node - b)))


def probe(registry) -> dict:
    """Both clauses over one registry, by the probe's own deviations."""
    detections = [(i, ix) for i, ix in registry if not ix.is_synthetic]
    clause1 = {i: ix.resolution_ratio for i, ix in detections
               if ix.resolution_ratio is not None and ix.resolution_ratio < THRESHOLD}
    groups: dict = {}
    for i, ix in detections:
        groups.setdefault((ix.manifold_a_key, ix.manifold_b_key), []).append((i, ix))
    clause2, missing = [], 0
    for members in groups.values():
        members.sort(key=lambda m: m[1].unstable_cdist)
        for (i, ix), (j, jx) in zip(members, members[1:]):
            hi, hj = deviation(ix), deviation(jx)
            if hi is None or hj is None:
                missing += 1
                continue
            sep = separation(ix, jx)
            if sep is not None and sep <= max(hi, hj):
                clause2.append((i, j, sep, hi, hj))
    out = set(clause1) | {k for i, j, *_ in clause2 for k in (i, j)}
    ratios = [ix.resolution_ratio for _, ix in detections if ix.resolution_ratio is not None]
    library = list(registry.unresolved_crossings(THRESHOLD))
    return {
        "registry": len(registry),
        "detections": len(detections),
        "N": len(registry) - len(out),
        "unresolved": sorted(out),
        "clause1": clause1,
        "clause2": clause2,
        "min_ratio": min(ratios, default=None),
        "missing_h": missing,
        "library_unresolved": library,
        "agrees": sorted(library) == sorted(out),
        "anycurve_extra": _anycurve_extra(detections, out),
    }


def _anycurve_extra(detections, out: set) -> list:
    """Same-curve pairs closer by arc length than their deviation that clause 2 missed."""
    by_curve: dict = {}
    for i, ix in detections:
        by_curve.setdefault((id(ix.unstable_manifold), ix.manifold_b_key), []).append((i, ix))
    extra = []
    for members in by_curve.values():
        if len(members) < 2:
            continue
        curve = members[0][1].unstable_manifold
        nodes = list(curve._walk_nodes(stop_at_final=False))
        points = np.array([n.get_point() for n in nodes], float)
        arc = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(points, axis=0).T))])
        index = {end(n): k for k, n in enumerate(nodes)}
        placed = []
        for i, ix in members:
            ks = [index.get(e) for e in ix.unstable_segment]
            if None in ks:
                continue
            k = min(ks)
            placed.append((arc[k] + float(np.hypot(*(np.asarray(ix.coords) - points[k]))), i, ix))
        placed.sort(key=lambda t: t[0])
        for (s, i, ix), (t, j, jx) in zip(placed, placed[1:]):
            h = max(deviation(ix) or 0.0, deviation(jx) or 0.0)
            if t - s <= h and not {i, j} <= out:
                extra.append((i, j, t - s, h))
    return extra


def line(name: str, result: dict) -> str:
    """One printable line per build."""
    clause2 = [f"({i},{j} sep {s:.2g} h {hi:.2g}/{hj:.2g})" for i, j, s, hi, hj in result["clause2"]]
    clause1 = {i: f"{r:.3g}" for i, r in result["clause1"].items()}
    ratio = result["min_ratio"]
    return (
        f"{name}: registry={result['registry']} detections={result['detections']} "
        f"|N|={result['N']} min_ratio={'-' if ratio is None else f'{ratio:.4g}'} "
        f"clause1={clause1} clause2={clause2} unresolved={result['unresolved']} "
        f"library_unresolved(3)={result['library_unresolved']} agrees={result['agrees']} "
        f"anycurve_extra={result['anycurve_extra']} missing_h={result['missing_h']}"
    )


def crossing_state(workbench):
    """``(N, bridge ids)``: the registry minus its unresolved detections, and the cut."""
    registry = workbench.intersection_registry
    unresolved = set(registry.unresolved_crossings(THRESHOLD))
    return (
        frozenset(i for i in registry.all_ids() if i not in unresolved),
        frozenset(b.id for b in workbench.bridges if b.id is not None),
    )


def recompute(workbench) -> None:
    """The plain recompute every path must agree with."""
    workbench.compute_intersections(workbench.fixed_points, preserve_ids=True, diagnose=False)
    workbench.rebuild_bridges()


def run_builds(names) -> None:
    """One line per named build."""
    from cases import BUILDERS

    for name in names:
        start = time.time()
        try:
            case = BUILDERS[name]()
        except Exception as exc:  # a set-aside case (orientation_reversing) raises
            print(f"{name}: build raised {type(exc).__name__}: {exc}", flush=True)
            continue
        print(line(name, probe(case.registry)) + f" ({time.time() - start:.0f}s)", flush=True)


def _series_session(name: str):
    """The series' build, its blasted zone and fixed point."""
    from cases import SERIES

    case = SERIES[name].build()
    zone = max(case.zones, key=lambda z: z.area)
    return case.session, zone, zone.fixed_point


def run_series(names) -> None:
    """One line per series step; a recompute at the last step."""
    from cases import SERIES, _blast

    for name in names:
        series = SERIES[name]
        session, zone, fp = _series_session(name)
        for blast in range(series.blasts + 1):
            if blast:
                _blast(session, zone, [fp], diagnose=False, min_separation=series.min_separation)
            session.topology()
            print(line(f"{name} blast {blast}", probe(session.workbench.intersection_registry)),
                  flush=True)
        before = crossing_state(session.workbench)
        recompute(session.workbench)
        after = crossing_state(session.workbench)
        print(f"  recompute: |N| {len(before[0])} -> {len(after[0])}, bridges "
              f"{len(before[1])} -> {len(after[1])}, N equal {before[0] == after[0]}, "
              f"bridges equal {before[1] == after[1]}", flush=True)


def run_k28x8() -> None:
    """k=2.8 at 8 blasts: the topology, the unresolved set, the holes that read it."""
    from cases import SERIES, _blast
    from tanglepack.topology.BridgeClass import _image_chain

    series = SERIES["k28_series"]
    session, zone, fp = _series_session("k28_series")
    for blast in range(9):
        if blast:
            _blast(session, zone, [fp], diagnose=False, min_separation=series.min_separation)
        session.topology()
    topo = session.topology()
    dynamics = session.symbolic_dynamics()
    print(f"k28 x8: holes={len(topo.trellis.holes)} seed={len(session.bridge_classes())} "
          f"closed={len(dynamics.classes)} reliable={dynamics.is_reliable}", flush=True)
    result = probe(session.workbench.intersection_registry)
    print(line("k28 x8", result), flush=True)

    before = crossing_state(session.workbench)
    recompute(session.workbench)
    after = crossing_state(session.workbench)
    recomputed = probe(session.workbench.intersection_registry)
    print(line("k28 x8 recomputed", recomputed), flush=True)
    print(f"  recompute: registry {result['registry']} -> {recomputed['registry']}, |N| "
          f"{len(before[0])} -> {len(after[0])}, N equal {before[0] == after[0]}, bridges "
          f"{len(before[1])} -> {len(after[1])} equal {before[1] == after[1]}", flush=True)

    registry = session.workbench.intersection_registry
    reasons = {i: f"clause 1 (ratio {r:.4g})" for i, r in recomputed["clause1"].items()}
    for i, j, sep, hi, hj in recomputed["clause2"]:
        for k, partner in ((i, j), (j, i)):
            reasons.setdefault(k, f"clause 2 (partner {partner}, sep {sep:.3g} <= max h {max(hi, hj):.3g})")
    for k in recomputed["unresolved"]:
        ix = registry[k]
        print(f"    id {k}: u={ix.unstable_cdist:.12g} ratio={ix.resolution_ratio:.4g} "
              f"coords=({ix.coords[0]:.12g}, {ix.coords[1]:.12g}) {reasons[k]}", flush=True)

    excluded = set(recomputed["unresolved"])
    for hole in topo.trellis.holes:
        bridge = tuple(hole.bounding_ids or ())
        chain = (_image_chain(topo.trellis, bridge) if bridge else None) or []
        by_id = excluded & set(bridge)
        by_chain = excluded & {k for pair in chain for k in pair}
        if by_id or by_chain:
            print(f"  hole {bridge} (iterate {hole.iterate}) reads excluded ids: by bounding id "
                  f"{sorted(by_id)}, by image chain {sorted(by_chain)}", flush=True)
    print("  (holes reading no excluded id are not listed)", flush=True)


def main(argv: list[str]) -> None:
    """Parse the flags and run the probes asked for."""
    from cases import BUILDERS, LAW_CASES, SERIES

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--laws", action="store_true", help="every law case")
    parser.add_argument("--builders", nargs="*", help="builders (default: every non-law one)")
    parser.add_argument("--series", nargs="*", help="blast series (default: all of them)")
    parser.add_argument("--k28x8", action="store_true", help="k=2.8 at 8 blasts in detail")
    args = parser.parse_args(argv)
    if args.laws:
        run_builds(LAW_CASES)
    if args.builders is not None:
        run_builds(args.builders or [name for name in BUILDERS if name not in LAW_CASES])
    if args.series is not None:
        run_series(args.series or list(SERIES))
    if args.k28x8:
        run_k28x8()


if __name__ == "__main__":
    main(sys.argv[1:])
