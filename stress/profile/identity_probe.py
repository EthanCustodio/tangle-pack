"""
Base-vs-branch identity probe of the numerics-arrays campaign (D3).

Every slice of the campaign keeps the points, the registry and the topology
byte-identical where it can, and lists every difference it cannot. This probe
is how the test engineer checks it: dump one checkout, dump the other, compare.

Usage (from the repository root)::

    TANGLEPACK_REPO=../tangle-pack-<slug>-base env/bin/python stress/profile/identity_probe.py dump --out DIR/base
    TANGLEPACK_REPO=../tangle-pack-<slug>      env/bin/python stress/profile/identity_probe.py dump --out DIR/branch
    env/bin/python stress/profile/identity_probe.py compare DIR/base DIR/branch   # exit 1 on any difference

``dump`` builds every ``cases.LAW_CASES`` case, every step of every
``cases.SERIES`` and both stages of every :data:`SCENARIOS` entry, each case,
series or scenario in a FRESH interpreter that imports ``tanglepack`` and
``cases`` from ``$TANGLEPACK_REPO`` (default: this checkout). A series runs
``cases.blast_series``'s loop: one blast per step at the series'
``min_separation``, the topology, then a plain recompute of the crossings and
bridges. A scenario blasts the nested zones unfiltered in its order and is
dumped after the blasts (``.blasted``) and after one recompute
(``.recomputed``). Every blast of the probe runs with ``diagnose=True``, which
changes no other item (N2, 2026-10-24: all 62 dumps of the base, compared). Per
build or step it writes ``<name>.npz`` (float arrays) and ``<name>.json``
(texts):

* ``curve:<key>`` -- ``(x, y, cdist)`` of every registered manifold (its whole
  list; the trimmed extent is a text), ``image:<i>`` of every blast image
  curve, ``bridge:<id>`` of every bridge (``partial:<i>`` for a partial piece);
* ``registry`` -- one row per crossing: id, both cdists, coordinates, sign,
  resolution ratio, chord deviation and the cdists of its four bracket ends;
* texts -- extents, the crossings' branch keys, N, the iterate table,
  ``unresolved_links``, the holes, ``pip_depth``, the class table, the
  dynamics and the refined transition matrix;
* ``same_stability_counts`` -- one line per blast since the previous dump (a
  law case: every blast of its build; a series step: its own blast; a
  scenario: its blasts, in ``.blasted``): the ``SameStabilityCount`` its
  resolution report carries, partial bridges named without their address.

``compare`` compares bytes. Each differing array gets one D3 line (rows
differing, length change, max |dxy|, max |dcdist|), each differing text a
unified diff; the last line reads ``identical: ...`` or ``differences in K items``.

Dev Notes:
    Every read of a curve's nodes, a crossing's brackets or the image curves
    goes through one reader below (``curve_arrays``, ``curve_extent``,
    ``bracket_cdists``, ``image_curves``, ``key_text``); registry fields,
    ``bridge.id`` and topology ``describe()`` are read directly; a blast's
    count is caught by wrapping ``TangleWorkbench.report_resolution``
    (``record_counts``), since ``blast_zone`` returns no report. A slice that
    changes an accessor edits only these, in its own diff. Nothing from
    ``src/`` is imported at module level.

    The scenarios (N1, 2026-10-24) are the builds where the edge dedup of an
    image curve against EARLIER image curves decides results; none of the law
    cases or series reaches it. ``nested_unfiltered`` is
    ``test_an_unfiltered_nested_blast_does_not_depend_on_order``'s outer-first
    order, ``shared_straddle`` the setup of
    ``test_blast_through_a_shared_straddle``.
"""

from __future__ import annotations

import argparse
import difflib
import json
import logging
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

REPO = Path(os.environ.get("TANGLEPACK_REPO", Path(__file__).resolve().parents[2])).resolve()

# Unfiltered one-iteration blasts of the nested build's zones, in order.
SCENARIOS = {
    "nested_unfiltered": ("outer",) * 4 + ("inner",) * 2,
    "shared_straddle": ("inner", "inner", "outer"),
}


# --------------------------------------------------------------------------- #
# Readers: the only places that touch a library object's shape
# --------------------------------------------------------------------------- #
def curve_arrays(curve: Any, whole: bool) -> np.ndarray:
    """``(x, y, cdist)`` per node, root first: the whole list, or root to tail."""
    nodes = (curve.get_cdist_array(return_nodes=True) if whole
             else curve.get_point_array(return_nodes=True))
    rows = [(*node.get_point(), node.get_cdist(curve.stability)) for node in nodes]
    return np.array(rows, dtype=np.float64).reshape(-1, 3)


def curve_extent(curve: Any) -> int:
    """How many nodes lie from the root to the (possibly trimmed) tail."""
    return len(curve.get_point_array(return_nodes=True))


def bracket_cdists(crossing: Any) -> list[float]:
    """The cdists of the crossing's two unstable then two stable bracket ends (NaN if none).

    Since N5a a bracket is store rows, read through ``segment_ends``; a tree
    before it holds nodes, read as before, so one probe dumps both.
    """
    if hasattr(crossing, "segment_ends"):
        if crossing.unstable_segment is None:
            return [np.nan] * 4
        return [float(c) for stability in ("unstable", "stable")
                for c in crossing.segment_ends(stability)[1]]
    ends = []
    for segment, stability in ((crossing.unstable_segment, "unstable"),
                               (crossing.stable_segment, "stable")):
        ends += ([node.get_cdist(stability) for node in segment] if segment is not None
                 else [np.nan, np.nan])
    return ends


def image_curves(workbench: Any) -> list:
    """The blast image curves the workbench keeps, in its order."""
    return list(workbench._image_curves)


def key_text(workbench: Any, key: Any) -> str:
    """A ``ManifoldKey`` as text, its fixed point by workbench position."""
    if key is None:
        return "None"
    fp, stability, orbit, branch = key
    return f"fp{workbench.fixed_points.index(fp)}.{stability}.{orbit}.{branch}"


def count_text(count: Any) -> str:
    """A ``SameStabilityCount`` as one line; a partial bridge's name loses its address."""
    pairs = [re.sub(r"partial image bridge 0x[0-9a-f]+", "partial image bridge", pair)
             for pair in count.curve_pairs]
    return (f"in_region={count.in_region} outside={count.outside} "
            f"cutoff={count.suggested_area_cutoff!r} pairs={pairs}")


def record_counts(workbench_class: Any) -> list:
    """Wrap ``report_resolution`` to collect every blast's ``same_stability_count``."""
    counts = []
    report_resolution = workbench_class.report_resolution

    def recording(self, *args, **kwargs):
        report = report_resolution(self, *args, **kwargs)
        if report.same_stability_count is not None:
            counts.append(report.same_stability_count)
        return report

    workbench_class.report_resolution = recording
    return counts


def _float(value: Any) -> float:
    return np.nan if value is None else float(value)


# --------------------------------------------------------------------------- #
# Dump (inside the per-case interpreter)
# --------------------------------------------------------------------------- #
def numerics_state(workbench: Any) -> tuple[dict, dict]:
    """The arrays and texts of the manifolds, bridges and registry."""
    arrays, extents = {}, []
    for key, manifold in workbench.manifolds.items():
        name = f"curve:{key_text(workbench, key)}"
        arrays[name] = curve_arrays(manifold, whole=True)
        extents.append(f"{name} {curve_extent(manifold)}")
    for i, curve in enumerate(image_curves(workbench)):
        arrays[f"image:{i}"] = curve_arrays(curve, whole=False)
    partial = 0
    for bridge in workbench.bridges:
        if bridge.id is None:
            name, partial = f"partial:{partial}", partial + 1
        else:
            name = f"bridge:{bridge.id[0]},{bridge.id[1]}"
        arrays[name] = curve_arrays(bridge, whole=False)

    registry = workbench.intersection_registry
    ids = sorted(registry.all_ids())
    rows, keys = [], []
    for i in ids:
        ix = registry[i]
        x, y = (np.nan, np.nan) if ix.coords is None else ix.coords
        rows.append([i, _float(ix.unstable_cdist), _float(ix.stable_cdist), x, y,
                     ix.crossing_sign, _float(ix.resolution_ratio), _float(ix.chord_deviation),
                     *bracket_cdists(ix)])
        keys.append(f"{i} {key_text(workbench, ix.manifold_a_key)} "
                    f"{key_text(workbench, ix.manifold_b_key)}")
    arrays["registry"] = np.array(rows, dtype=np.float64).reshape(-1, 12)
    texts = {
        "extents": "\n".join(extents),
        "crossing_keys": "\n".join(keys),
        "resolved_ids": " ".join(map(str, sorted(
            registry.resolved_ids(workbench.resolution_ratio_threshold)))),
        "iterate_table": "\n".join(f"{s} {n} {t}" for s, n, t in sorted(registry.iterate_table.items())),
        "unresolved_links": "\n".join(f"{s} {sorted(c)}" for s, c in sorted(registry.unresolved_links.items())),
    }
    return arrays, texts


def topology_texts(session: Any) -> dict:
    """The topology outputs as text; a build whose topology raises records the error."""
    try:
        topology = session.topology()
        labels, matrix = topology.dynamics.transition_matrix(refined=True)
        depths = []
        for fp in session.workbench.fixed_points:
            chosen = session.chosen_strong_pip(fp)
            if chosen is None:
                continue
            trellis = session.trellis(fp)
            pip, m = trellis.registry.find(chosen), 0
            while trellis.iterate(pip, m + 1) is not None:
                m += 1
            depths.append(f"fp{session.workbench.fixed_points.index(fp)} pip {pip} depth {m}")
        holes = sorted(
            f"{h.bounding_ids} side={h.bridge_side} iterate={h.iterate} "
            f"near={h.near_intersection_id} origin={h.origin} openings="
            f"{[(i, which, row, key_text(session.workbench, key)) for i, which, row, key in h.openings or ()]}"
            for h in topology.holes)
        return {
            "strong_pips": str(topology.strong_pips),
            "pip_depth": "\n".join(depths),
            "holes": "\n".join(holes),
            "classes": topology.classes.describe(),
            "dynamics": topology.dynamics.describe(),
            "transition_matrix": f"{labels}\n{matrix.tolist()}",
        }
    except Exception as exc:  # noqa: BLE001 -- a raising topology is an output too
        return {"topology_error": f"{type(exc).__name__}: {exc}"}


def save(out: Path, name: str, session: Any, counts: list) -> None:
    """Write ``<name>.npz`` and ``<name>.json`` for the session's current state,
    taking the blast counts recorded since the previous dump out of ``counts``."""
    arrays, texts = numerics_state(session.workbench)
    texts.update(topology_texts(session))
    texts["same_stability_counts"] = "\n".join(map(count_text, counts))
    counts.clear()
    np.savez(out / f"{name}.npz", **arrays)
    (out / f"{name}.json").write_text(json.dumps(texts, indent=1, sort_keys=True))
    print(f"  {name}: {len(arrays)} arrays, {len(texts)} texts", flush=True)


def recompute(workbench: Any) -> None:
    """The plain recompute of the crossings and bridges a series runs after each step."""
    workbench.compute_intersections(workbench.fixed_points, preserve_ids=True, diagnose=False)
    workbench.rebuild_bridges()


def run_worker(name: str, out: Path) -> None:
    """Dump one law case, every step of one series, or one scenario, in this interpreter."""
    sys.dont_write_bytecode = True
    sys.path[:0] = [str(REPO / "src"), str(REPO / "tests")]
    logging.basicConfig(level=logging.ERROR)
    import cases
    import tanglepack

    assert Path(tanglepack.__file__).resolve().is_relative_to(REPO / "src"), tanglepack.__file__
    counts = record_counts(tanglepack.TangleWorkbench)
    if name in cases.LAW_CASES:
        save(out, name, cases.BUILDERS[name]().session, counts)
        return
    if name in SCENARIOS:
        case = cases.build_nested(outer_blasts=0, through="zones")
        zones = sorted(case.zones, key=lambda zone: zone.area)
        zone = {"inner": zones[0], "outer": zones[-1]}
        for which in SCENARIOS[name]:
            case.session.blast_zone(zone[which], num_iterations=1, fixed_point=None,
                                    min_separation=1e-4, diagnose=True)
        save(out, f"{name}.blasted", case.session, counts)
        recompute(case.session.workbench)
        save(out, f"{name}.recomputed", case.session, counts)
        return
    series = cases.SERIES[name]
    case = series.build()
    session = case.session
    zone = max(case.zones, key=lambda zone: zone.area)
    for blast in range(series.blasts + 1):
        if blast:
            session.blast_zone(zone, num_iterations=1, fixed_point=[zone.fixed_point],
                               min_separation=series.min_separation, diagnose=True)
        save(out, f"{name}.b{blast:02d}", session, counts)
        recompute(session.workbench)


def dump(out: Path) -> int:
    """Run every law case, series and scenario in its own interpreter against ``REPO``."""
    out.mkdir(parents=True, exist_ok=True)
    names = subprocess.check_output(
        [sys.executable, "-B", "-c",
         "import sys; sys.path[:0] = [sys.argv[1]]; import cases; "
         "print(' '.join(cases.LAW_CASES + tuple(cases.SERIES)))", str(REPO / "tests")],
        env={**os.environ, "PYTHONPATH": str(REPO / "src")}, text=True).split() + list(SCENARIOS)
    env = {**os.environ, "TANGLEPACK_REPO": str(REPO), "MPLBACKEND": "Agg",
           "PYTHONDONTWRITEBYTECODE": "1"}
    start = time.perf_counter()
    for name in names:
        t = time.perf_counter()
        subprocess.run([sys.executable, str(Path(__file__).resolve()), "worker", name,
                        "--out", str(out)], env=env, check=True)
        print(f"{name}: {time.perf_counter() - t:.0f} s", flush=True)
    print(f"dumped {REPO} into {out} in {time.perf_counter() - start:.0f} s")
    return 0


# --------------------------------------------------------------------------- #
# Compare
# --------------------------------------------------------------------------- #
def _rows_differing(a: np.ndarray, b: np.ndarray) -> int:
    differ = (a != b) & ~(np.isnan(a) & np.isnan(b))
    return int(differ.reshape(len(a), -1).any(axis=1).sum())


def array_difference(where: str, a: np.ndarray, b: np.ndarray) -> str:
    """The D3 line of two arrays whose bytes differ."""
    n = min(len(a), len(b))
    head_a, head_b = a[:n], b[:n]
    line = (f"{where}: {_rows_differing(head_a, head_b)} of {n} common rows differ, "
            f"length {len(a)} -> {len(b)}")
    if a.ndim == 2 and a.shape[1:] == b.shape[1:] == (3,) and n:
        dxy = float(np.nanmax(np.abs(head_a[:, :2] - head_b[:, :2]), initial=0.0))
        dcdist = float(np.nanmax(np.abs(head_a[:, 2] - head_b[:, 2]), initial=0.0))
        line += f", max |dxy| {dxy!r}, max |dcdist| {dcdist!r}"
    elif a.shape[1:] == b.shape[1:] and n:
        line += f", max |d| {float(np.nanmax(np.abs(head_a - head_b), initial=0.0))!r}"
    return line


def compare(base: Path, branch: Path) -> int:
    """Compare two dumps item by item; print the differences and a last summary line."""
    differences, n_arrays, n_texts = [], 0, 0
    names = sorted({p.stem for p in base.glob("*.json")} | {p.stem for p in branch.glob("*.json")})
    for name in names:
        if not ((base / f"{name}.json").exists() and (branch / f"{name}.json").exists()):
            differences.append(f"{name}: dumped only in {base if (base / f'{name}.json').exists() else branch}")
            continue
        with np.load(base / f"{name}.npz") as old, np.load(branch / f"{name}.npz") as new:
            for item in sorted(set(old.files) | set(new.files)):
                n_arrays += 1
                if item not in old.files or item not in new.files:
                    differences.append(f"{name} {item}: only in {'branch' if item in new.files else 'base'}")
                elif old[item].dtype != new[item].dtype or old[item].tobytes() != new[item].tobytes():
                    differences.append(array_difference(f"{name} {item}", old[item], new[item]))
        old_texts = json.loads((base / f"{name}.json").read_text())
        new_texts = json.loads((branch / f"{name}.json").read_text())
        for item in sorted(set(old_texts) | set(new_texts)):
            n_texts += 1
            a, b = old_texts.get(item, ""), new_texts.get(item, "")
            if a != b:
                diff = difflib.unified_diff(a.splitlines(), b.splitlines(), f"base {name} {item}",
                                            f"branch {name} {item}", lineterm="")
                differences.append("\n".join(diff))
    for difference in differences:
        print(difference)
    if differences:
        print(f"differences in {len(differences)} items")
        return 1
    print(f"identical: {len(names)} builds, {n_arrays} arrays, {n_texts} texts")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("dump").add_argument("--out", type=Path, required=True)
    pair = commands.add_parser("compare")
    pair.add_argument("base", type=Path)
    pair.add_argument("branch", type=Path)
    worker = commands.add_parser("worker")
    worker.add_argument("name")
    worker.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "dump":
        return dump(args.out)
    if args.command == "compare":
        return compare(args.base, args.branch)
    run_worker(args.name, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
