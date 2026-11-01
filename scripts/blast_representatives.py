"""
The representative blast against the full one: equivalence, cost, and the prune question.

Run 4 of the deeper-blasts campaign (``.refactor/campaigns/2026-10-29-deeper-blasts.md``).
``blast_zone(representatives=...)`` iterates one bridge per refined class; this
script measures what that changes. Subcommands (one process each;
``PYTHONPATH=src env/bin/python scripts/blast_representatives.py <cmd>``):

* ``equivalence`` -- every ``cases.SERIES`` and the law cases that blast
  (``nested_inner_blasts``: the inner zone twice, then the outer twice;
  ``heteroclinic``: the outer zone twice), one iteration at a time to the full
  series' depth, full against each policy. Per step it compares the alphabet
  identities, the closed alphabet, the words per identity, the refined
  transition matrix up to relabelling, ``pip_depth``, the
  pseudoneighbor pairs (references apart) and the holes, and records Theorem T
  (``helpers.laws.p1_itinerary_disagreements``) on both. The two builds number
  their crossings differently, so a crossing of one is matched to the other's
  by its branch keys and both cdists to ``DEDUP_RTOL`` (``registry.find``'s
  rule, the fixed points by their position). ``equivalence.csv``.
* ``cost --cutoff 1e-7 --blasts 12`` -- the k=2.5 scale series
  (``stress/profile/profile_numerics.scale_build``) blasted full and with
  representatives, at ``min_separation`` 1e-5 and None: bridges, points and
  seconds per blast. ``cost.csv``.
* ``lean --cutoff 1e-7`` -- lean plus representatives, unguarded, one
  iteration per call until the time budget, the point budget or an error
  (the cdist ceiling shows up as one). ``lean.csv``.
* ``prune --cutoff 1e-9 --separation none`` -- one prune run: representatives,
  one iteration per call under a time and point budget, a topology snapshot per
  blast (``prune/<cutoff>_<separation>.json``). ``prune-report`` compares the
  runs of one cutoff at their common pip depth (``prune.csv``).
* ``table`` -- ``table.txt`` from every CSV on disk.

Output: ``figures/blast_representatives/``. Pin a timing run (``taskset -c 28-31``).
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Optional

import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "tests"), str(ROOT / "stress" / "profile")]

from cases import SERIES, build_nested  # noqa: E402
from helpers.laws import p1_itinerary_disagreements  # noqa: E402
from tanglepack.numerics.IntersectionRegistry import DEDUP_RTOL  # noqa: E402

OUT = ROOT / "figures" / "blast_representatives"
POLICIES = ("fewest_points", "smallest_id", "largest_id")


# --------------------------------------------------------------------------- #
# Snapshots: what one blast left, in crossing descriptors two builds can match
# --------------------------------------------------------------------------- #
def _key(fixed_points, key) -> Optional[list]:
    """A branch key with its fixed point as a position."""
    return None if key is None else [fixed_points.index(key[0]), key[1], key[2], key[3]]


def crossings(session) -> dict[str, list]:
    """Every registered crossing as ``[unstable key, stable key, unstable cdist, stable cdist]``."""
    fps = list(session.workbench.fixed_points)
    return {
        str(i): [_key(fps, ix.manifold_a_key), _key(fps, ix.manifold_b_key),
                 float(ix.unstable_cdist), float(ix.stable_cdist)]
        for i, ix in session.workbench.intersection_registry
    }


def snapshot(session, result, fixed_point, *, with_crossings: bool = True) -> dict[str, Any]:
    """The topology of the whole session after one blast, ids kept (with every crossing, unless
    ``with_crossings`` is False: a run then stores its registry once, ids never change in a blast)."""
    workbench = session.workbench
    fps = list(workbench.fixed_points)
    topology = session.topology()
    dynamics = topology.dynamics

    def element(identity) -> list:
        branch_key, side, lo, hi = identity
        return [*_key(fps, branch_key), side, lo, hi]

    def identity(bridge_class) -> list:
        return sorted((element(topology.homotopy.element(ref).identity)
                       for ref in (bridge_class.source, bridge_class.target)), key=json.dumps)

    words = {
        json.dumps(identity(cls)): [[identity(s.bridge_class), s.direction] for s in cd.symbols]
        for cls, cd in dynamics.classes.items() if cd.itinerary is not None
    }
    names, matrix = dynamics.transition_matrix()
    trellis = topology.trellis
    tested, _skipped, disagreements = p1_itinerary_disagreements(session)
    snap = {
        "alphabet": [identity(entry.bridge_class) for entry in topology.classes],
        "closed": [identity(cls) for cls in dynamics.classes],
        "words": words,
        "matrix": matrix.tolist(),
        "pairs": sorted(list(p.as_tuple()) for p in trellis.pseudoneighbors),
        "references": sorted(list(p.as_tuple()) for p in trellis.reference_pseudoneighbors),
        "holes": sorted(([list(h.bounding_ids), h.iterate, h.bridge_side] for h in trellis.holes),
                        key=json.dumps),
        "pip_depth": result.pip_depth.get(fixed_point) if result else None,
        "theorem_t": [tested, len(disagreements)],
        "bridges": len(workbench.bridges),
        "points": len(workbench.point_store),
        "held": result.held if result else 0,
    }
    if with_crossings:
        snap["crossings"] = crossings(session)
    return snap


def _same(a: list, b: list) -> bool:
    """Two crossing descriptors are one crossing: equal keys, both cdists to ``DEDUP_RTOL``."""
    if a[0] != b[0] or a[1] != b[1]:
        return False
    return all(abs(x - y) <= DEDUP_RTOL * max(abs(x), abs(y)) for x, y in zip(a[2:], b[2:]))


def translate(snap: dict, into: dict) -> dict:
    """``snap`` with its crossing ids replaced by ``into``'s (``-id - 1`` where it has none)."""
    by_keys: dict[str, list[tuple[int, list]]] = {}
    for i, d in into["crossings"].items():
        by_keys.setdefault(json.dumps(d[:2]), []).append((int(i), d))
    ids = {}
    for i, d in snap["crossings"].items():
        match = next((j for j, e in by_keys.get(json.dumps(d[:2]), []) if _same(d, e)), None)
        ids[int(i)] = match if match is not None else -int(i) - 1

    def walk(value):
        if isinstance(value, list) and len(value) == 7 and isinstance(value[4], str):
            return [*value[:5], *(ids.get(x, x) if isinstance(x, int) else x for x in value[5:])]
        if isinstance(value, list) and len(value) == 2 and all(isinstance(x, int) for x in value):
            return sorted(ids[x] for x in value)
        if isinstance(value, list) and len(value) == 3 and isinstance(value[2], str):
            return [[ids[x] for x in value[0]], value[1], value[2]]
        return [walk(v) for v in value] if isinstance(value, list) else value

    out = dict(snap)
    for part in ("alphabet", "closed", "pairs", "references", "holes"):
        out[part] = walk(snap[part])
    out["words"] = {
        json.dumps(sorted(walk(json.loads(k)), key=json.dumps)): [
            [sorted(walk(i), key=json.dumps), d] for i, d in w
        ]
        for k, w in snap["words"].items()
    }
    return out


def _canon(items: list) -> set[str]:
    """A part of a snapshot as a set of texts (each identity's elements sorted)."""
    return {json.dumps(sorted(x, key=json.dumps) if x and isinstance(x[0], list) and len(x) == 2 and
                       all(isinstance(e, list) and len(e) == 7 for e in x) else x)
            for x in items}


def _isomorphic(a: list, b: list) -> bool:
    """Two transition matrices are one up to a relabelling of their symbols."""
    if len(a) != len(b):
        return False
    graphs = []
    for m in (np.array(a, dtype=int).reshape(len(a), -1), np.array(b, dtype=int).reshape(len(b), -1)):
        g = nx.DiGraph()
        g.add_nodes_from(range(len(m)))
        g.add_weighted_edges_from((i, j, int(m[i, j])) for i, j in zip(*np.nonzero(m)))
        graphs.append(g)
    return nx.is_isomorphic(*graphs, edge_match=lambda x, y: x["weight"] == y["weight"])


def compare(full: dict, other: dict) -> dict[str, Any]:
    """Every difference of ``other`` against ``full`` (ids of ``other`` translated into ``full``'s)."""
    other = translate(other, full)
    row: dict[str, Any] = {}
    for part in ("alphabet", "closed", "pairs", "references", "holes"):
        a, b = _canon(full[part]), _canon(other[part])
        row[f"{part}_only_full"], row[f"{part}_only_other"] = len(a - b), len(b - a)
    common = set(full["words"]) & set(other["words"])
    row["words_differ"] = sum(full["words"][k] != other["words"][k] for k in common)
    row["matrix_iso"] = _isomorphic(full["matrix"], other["matrix"])
    row["pip_depth_full"], row["pip_depth_other"] = full["pip_depth"], other["pip_depth"]
    row["T_tested_full"], row["T_fail_full"] = full["theorem_t"]
    row["T_tested_other"], row["T_fail_other"] = other["theorem_t"]
    row["differs"] = bool(
        any(row[f"{p}_only_{s}"] for p in ("alphabet", "closed", "pairs", "references", "holes")
            for s in ("full", "other"))
        or row["words_differ"] or not row["matrix_iso"]
        or full["pip_depth"] != other["pip_depth"]
    )
    return row


# --------------------------------------------------------------------------- #
# Equivalence
# --------------------------------------------------------------------------- #
def _runs() -> dict[str, tuple]:
    """``name -> (build, [(which zone, separation)] per blast)``: the series and the blasting law cases."""
    runs = {name: (s.build, [("largest", s.min_separation)] * s.blasts) for name, s in SERIES.items()}
    runs["nested_inner_blasts"] = (
        lambda: build_nested(outer_blasts=0, through="zones"),
        [("smallest", 1e-4)] * 2 + [("largest", 1e-4)] * 2,
    )
    runs["heteroclinic"] = (
        lambda: build_nested(outer_blasts=0, p3_unstable_steps=16, p3_stable_steps=10, through="zones"),
        [("largest", 1e-4)] * 2,
    )
    return runs


def _series(build, plan, policy) -> list[dict]:
    """One snapshot per blast of a run (the unblasted state first)."""
    case = build()
    session = case.session
    zones = sorted(case.zones, key=lambda zone: zone.area)
    snaps = [snapshot(session, None, None)]
    for which, separation in plan:
        zone = zones[-1] if which == "largest" else zones[0]
        result = session.blast_zone(zone, 1, fixed_point=[zone.fixed_point], min_separation=separation,
                                    diagnose=False, representatives=policy)
        snaps.append(snapshot(session, result, zone.fixed_point))
    return snaps


def equivalence(args) -> None:
    """Full against every policy on every run; ``equivalence.csv``."""
    rows = []
    for name, (build, plan) in _runs().items():
        if args.cases and name not in args.cases:
            continue
        t = time.perf_counter()
        full = _series(build, plan, None)
        for policy in POLICIES:
            other = _series(build, plan, policy)
            for blast, (a, b) in enumerate(zip(full, other)):
                rows.append({"case": name, "policy": policy, "blast": blast,
                             "bridges_full": a["bridges"], "bridges_other": b["bridges"],
                             "points_full": a["points"], "points_other": b["points"],
                             "held": b["held"], **compare(a, b)})
        print(f"{name}: {time.perf_counter() - t:.0f} s", flush=True)
    _write(OUT / "equivalence.csv", rows)


# --------------------------------------------------------------------------- #
# Cost and lean depth on the scale series
# --------------------------------------------------------------------------- #
def _scale(cutoff: float):
    from profile_numerics import scale_build
    return scale_build(cutoff)


def cost(args) -> None:
    """Bridges, points and seconds per blast, full against representatives; ``cost.csv``."""
    rows = []
    for separation in (1e-5, None):
        for policy in (None, "fewest_points"):
            session, fp, zone = _scale(args.cutoff)
            start = time.perf_counter()
            for blast in range(1, args.blasts + 1):
                t = time.perf_counter()
                result = session.blast_zone(zone, 1, fixed_point=[fp], min_separation=separation,
                                            diagnose=False, max_points=args.max_points,
                                            representatives=policy)
                rows.append({"cutoff": args.cutoff, "separation": separation, "policy": policy or "full",
                             "blast": blast, "bridges": len(session.workbench.bridges),
                             "points": len(session.workbench.point_store),
                             "iterated": sum(len(s.interior_parents) for s in result.steps),
                             "held": result.held, "pip_depth": result.pip_depth.get(fp),
                             "seconds": round(time.perf_counter() - t, 2), "stop": result.stop_reason})
                print(rows[-1], flush=True)
                if result.stop_reason == "point_budget" or time.perf_counter() - start > args.max_seconds:
                    rows[-1]["stop"] = result.stop_reason if result.stop_reason == "point_budget" else "time"
                    break
            _write(OUT / "cost.csv", rows)


def lean(args) -> None:
    """Lean plus representatives, unguarded, until a budget or an error; ``lean.csv``."""
    session, fp, zone = _scale(args.cutoff)
    rows, start = [], time.perf_counter()
    for blast in range(1, args.blasts + 1):
        t = time.perf_counter()
        try:
            result = session.blast_zone(zone, 1, fixed_point=[fp], min_separation=None, lean=True,
                                        max_points=args.max_points, representatives="fewest_points")
        except Exception as error:  # noqa: BLE001 - the stop is the measurement
            rows.append({"blast": blast, "stop": f"{type(error).__name__}: {str(error)[:160]}"})
            break
        store = session.workbench.point_store
        rows.append({"cutoff": args.cutoff, "blast": blast, "bridges": len(session.workbench.bridges),
                     "points": len(store), "released": result.released_rows,
                     "held": result.held, "pip_depth": result.pip_depth.get(fp),
                     "max_cdist": float(np.nanmax(store.cdist)),
                     "seconds": round(time.perf_counter() - t, 2), "stop": result.stop_reason})
        print(rows[-1], flush=True)
        _write(OUT / "lean.csv", rows)
        if result.stop_reason != "max_iterations":
            break
        if time.perf_counter() - start > args.max_seconds:
            rows[-1]["stop"] = "time"
            break
    _write(OUT / "lean.csv", rows)


# --------------------------------------------------------------------------- #
# The prune question
# --------------------------------------------------------------------------- #
def _separation(text: str):
    return None if text == "none" else text if text == "resolution" else float(text)


def prune(args) -> None:
    """One prune run: a snapshot per representative blast under the budgets."""
    session, fp, zone = _scale(args.cutoff)
    separation = _separation(args.separation)
    record: dict[str, Any] = {"cutoff": args.cutoff, "separation": args.separation, "blasts": [], "stop": None}
    path = OUT / "prune" / f"{args.cutoff:g}_{args.separation}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    for blast in range(1, args.blasts + 1):
        t = time.perf_counter()
        try:
            result = session.blast_zone(zone, 1, fixed_point=[fp], min_separation=separation,
                                        diagnose=False, max_points=args.max_points,
                                        representatives="fewest_points")
            snap = snapshot(session, result, fp, with_crossings=False)
        except Exception as error:  # noqa: BLE001 - the stop is the measurement
            record["stop"] = f"{type(error).__name__}: {str(error)[:200]}"
            break
        snap["blast"], snap["seconds"] = blast, round(time.perf_counter() - t, 2)
        record["blasts"].append(snap)
        print(f"{args.cutoff:g} {args.separation} blast {blast}: depth {snap['pip_depth']}, "
              f"{len(snap['alphabet'])} classes, {snap['points']} points, "
              f"{snap['seconds']} s", flush=True)
        record["crossings"] = crossings(session)
        path.write_text(json.dumps(record))
        if result.stop_reason != "max_iterations":
            record["stop"] = result.stop_reason
            break
        if time.perf_counter() - start > args.max_seconds:
            record["stop"] = "time"
            break
        if snap["points"] > args.max_total_points:
            record["stop"] = "total_points"
            break
    record["stop"] = record["stop"] or "max_iterations"
    path.write_text(json.dumps(record))


def prune_report(args) -> None:
    """Each cutoff's runs compared with the unguarded one, first at their common pip depth, then
    at the last blast all of them ran; ``prune.csv``."""
    rows = []
    runs: dict[str, dict[str, dict]] = {}
    for path in sorted((OUT / "prune").glob("*.json")):
        record = json.loads(path.read_text())
        for snap in record["blasts"]:
            snap["crossings"] = record.get("crossings", {})
        runs.setdefault(f"{record['cutoff']:g}", {})[record["separation"]] = record
    for cutoff, by_separation in runs.items():
        reached = {s: max((b["pip_depth"] or 0 for b in r["blasts"]), default=0)
                   for s, r in by_separation.items()}
        common = min(reached.values())
        at = {s: next((b for b in r["blasts"] if (b["pip_depth"] or 0) >= common), None)
              for s, r in by_separation.items()}
        last = min(len(r["blasts"]) for r in by_separation.values())
        for where, chosen in (("first at the common pip depth", at),
                              ("last common blast", {s: r["blasts"][last - 1] for s, r in by_separation.items()})):
            rows += _prune_rows(cutoff, where, chosen, by_separation, common, reached)
    _write(OUT / "prune.csv", rows)


def _prune_rows(cutoff, where, at, by_separation, common, reached) -> list[dict]:
    """One row per separation of one cutoff: its snapshot ``at`` compared with the unguarded one."""
    rows = []
    reference = at.get("none")
    for separation, snap in at.items():
        record = by_separation[separation]
        row = {"cutoff": cutoff, "compared": where, "separation": separation, "common_pip_depth": common,
               "deepest_pip": reached[separation],
               "blasts_run": len(record["blasts"]), "stop": record["stop"],
               "at_blast": snap["blast"] if snap else None,
               "classes": len(snap["alphabet"]) if snap else None}
        row["pip_depth"] = snap["pip_depth"] if snap else None
        if snap and reference and separation != "none":
            row.update({f"vs_none_{k}": v for k, v in compare(reference, snap).items()})
        rows.append(row)
    return rows


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #
def _write(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def table(args) -> None:
    """``table.txt``: every CSV, the equivalence differences listed row by row."""
    lines = []
    for name in ("equivalence", "cost", "lean", "prune"):
        path = OUT / f"{name}.csv"
        if not path.exists():
            continue
        rows = list(csv.DictReader(path.open()))
        if name == "equivalence":
            rows = [r for r in rows if r["differs"] == "True"]
            lines.append(f"== equivalence: {len(rows)} differing (case, policy, blast) rows ==")
        else:
            lines.append(f"== {name} ==")
        if rows:
            keys = [k for k in rows[0] if any(r[k] not in ("", "0", "False") for r in rows)
                    or k in ("case", "policy", "blast")]
            lines.append(" | ".join(keys))
            lines += [" | ".join(r[k] for k in keys) for r in rows]
        lines.append("")
    (OUT / "table.txt").write_text("\n".join(lines))
    print("\n".join(lines))


def main() -> None:
    logging.basicConfig(level=logging.ERROR)
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("equivalence").add_argument("--cases", nargs="*")
    for name in ("cost", "lean", "prune"):
        sub = commands.add_parser(name)
        sub.add_argument("--cutoff", type=float, default=1e-7)
        sub.add_argument("--blasts", type=int, default=12 if name == "cost" else 60)
        sub.add_argument("--max-points", type=int, default=6_000_000)
        sub.add_argument("--max-seconds", type=float, default=1800)
        if name == "prune":
            sub.add_argument("--separation", default="none")
            sub.add_argument("--max-total-points", type=int, default=15_000_000)
    commands.add_parser("prune-report")
    commands.add_parser("table")
    args = parser.parse_args()
    {"equivalence": equivalence, "cost": cost, "lean": lean, "prune": prune,
     "prune-report": prune_report, "table": table}[args.command](args)


if __name__ == "__main__":
    main()
