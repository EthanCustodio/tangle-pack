"""
Length growth of the unstable curve inside each k=2.5 resonance zone, unpruned.

A geometric cross-check of the topological entropy: for a smooth map the
growth rate of the length of an iterated curve is the topological entropy
(Yomdin, Newhouse), so ``ln L_n`` of the unstable curve the zone keeps should
grow with slope ``h``. The entropy sweep's own ``length_growth`` figure reads
the blast frontier, which ``min_separation`` prunes from about blast 8 on;
this script measures the same quantity without pruning, two ways:

* ``blast`` -- the library: ``blast_zone(num_iterations=1, min_separation=None)``
  one step at a time, ``L_n`` the summed arc length of the step's
  ``kept_interior`` bridges, which is the step's NEW layer (earlier steps'
  bridges are never re-iterated), so the figure draws its running total
  ``L_0 + sum_{j<=n} L_j``, comparable with the curve's f^n(W); without the blast's same-stability diagnostic
  (``--diagnose`` restores it: the peak then grows ~2.3x per step instead of
  ~1.66x for p1). It stops on ``--blast-gb``; ``--blast-only`` reruns it and
  keeps the curve already in the file.
* ``curve`` -- independent of the library past the build: the zone's initial
  interior curve (the first blast step's ``interior_parents``, the same curve
  ``blast`` starts from) is mapped as a NumPy polyline. Each segment is
  refined until the midpoint of its image lies within ``--tol`` of the image
  chord: the Hénon image of a straight segment is a parabola, whose largest
  chord deviation is the one at its midpoint, so the test is exact for the
  polyline. The image is cut where it crosses the zone's stable boundary arcs
  (an unstable curve can leave a zone only across the stable manifold, the
  CLAUDE.md invariant) and a piece is kept iff a point of it lies inside the
  zone polygon and outside every smaller zone (the blast's own rule,
  ``Blast.in_zone``). About 20 bytes per point, so it runs far deeper
  (``--curve-gb``).

Cases are those of ``henon_k25_entropy.py``: ``p1`` (the outer zone), ``p3``
(the inner zone), ``nested`` (the outer zone of the nested build: the annulus,
the inner zone excluded as the nested blast excludes it).

Each (case, cutoff) run writes ``figures/henon_k25_entropy/length_growth/data/
<case>_a<cutoff>.json`` after every step; ``--plot-only`` redraws
``length_growth`` (ln L_n with the fitted slope and the topological h, and the
per-step rate) from the files.

Run with ``PYTHONPATH=src env/bin/python scripts/henon_k25_length_growth.py``
(``--cases p1 p3 nested``, ``--cutoff 1e-7``, ``--tol 1e-5``, ``--blast-gb 24``,
``--curve-gb 20``, ``--max-steps 60``, ``--max-blast-steps 30``, ``--fit-from 4``,
``--plot-only``). At ``--tol`` 1e-5 the p1 lengths agree with 1e-6 and 1e-7 to
five digits through step 14, piece counts identical, on a third of the points.
For parallel runs launch one process per case, then ``--plot-only``.
"""

from __future__ import annotations

import argparse
import gc
import json
import logging
import sys
import time
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.path import Path as PolygonPath
from numpy.typing import NDArray

sys.path.insert(0, str(Path(__file__).resolve().parent))

import henon_k25_entropy as entropy  # noqa: E402
from symbolic_figures import configure_logging, save_panel  # noqa: E402

from tanglepack.examples.henon import henon_map  # noqa: E402
from tanglepack.numerics.geometry import arc_polyline  # noqa: E402

logger = logging.getLogger(__name__)

OUT_DIR = entropy.OUT_DIR / "length_growth"
DATA_DIR = OUT_DIR / "data"
#: Points per chunk of the curve map (bounds the temporaries).
CHUNK = 20_000_000
#: Bisection rounds of one segment's refinement (2^60 pieces: never reached by a sane curve).
MAX_ROUNDS = 60
#: Bytes held per curve point (two float64 coordinates and an int64 piece id).
POINT_BYTES = 24
#: A probe within this distance of a zone's unstable boundary arc is ON it (the
#: polygon's chords cut corners by up to ~3e-6 at a 1e-7 cutoff).
BOUNDARY_TOL = 2e-5
#: Fraction of a stable arc's length it is extended by past each end.
ARC_EXTENSION = 0.05


# ── the zone's geometry ─────────────────────────────────────────────────────


def stable_arcs(session, zone) -> list[NDArray[np.float64]]:
    """The zone's stable boundary arcs (anchor -> pip), as the zone itself builds them."""
    workbench = session.workbench
    registry = workbench.intersection_registry
    arcs = []
    for ix in zone.cut_intersections or [zone.boundary_intersection]:
        manifold = workbench.manifolds[ix.manifold_b_key]
        points = list(arc_polyline(manifold, 0.0, ix.stable_cdist, tol=registry.cdist_tol))
        points.append(np.asarray(ix.coords, dtype=float))
        arcs.append(_extended(np.vstack(points)))
    return arcs


def unstable_arcs(session, zone) -> list[NDArray[np.float64]]:
    """The zone's unstable boundary arcs (anchor -> pip): the only boundary the curve can lie ON."""
    workbench = session.workbench
    registry = workbench.intersection_registry
    arcs = []
    for ix in zone.cut_intersections or [zone.boundary_intersection]:
        manifold = workbench.manifolds[ix.manifold_a_key]
        points = list(arc_polyline(manifold, 0.0, ix.unstable_cdist, tol=registry.cdist_tol))
        points.append(np.asarray(ix.coords, dtype=float))
        arcs.append(np.vstack(points))
    return arcs


def _extended(arc: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    The arc lengthened by :data:`ARC_EXTENSION` of its length past both ends, along the end tangents.

    The image can cross the stable manifold just past the pip, outside the arc;
    a cut there is needed, and a superfluous cut is harmless (every piece is
    classified by its own probe point).
    """
    length = np.linalg.norm(np.diff(arc, axis=0), axis=1).sum()
    ends = []
    for end, inner in ((arc[0], arc[1]), (arc[-1], arc[-2])):
        tangent = (end - inner) / np.linalg.norm(end - inner)
        ends.append(end + ARC_EXTENSION * length * tangent)
    return np.vstack([ends[0], arc, ends[1]])


class Graph:
    """A monotone piece of a stable arc, as ``v = s(u)`` in a rotated frame."""

    def __init__(self, points: NDArray[np.float64]) -> None:
        """Rotate onto the piece's principal axis; ``u`` must be strictly monotone."""
        direction = points[-1] - points[0]
        self.e_u = direction / np.linalg.norm(direction)
        self.e_v = np.array([-self.e_u[1], self.e_u[0]])
        u, v = points @ self.e_u, points @ self.e_v
        assert np.all(np.diff(u) > 0), "a graph piece must be monotone"
        self.u, self.v = u, v

    def side(self, points: NDArray[np.float64]) -> NDArray[np.float64]:
        """``v - s(u)``, NaN outside the piece's ``u`` range."""
        u, v = points @ self.e_u, points @ self.e_v
        out = v - np.interp(u, self.u, self.v)
        out[(u < self.u[0]) | (u > self.u[-1])] = np.nan
        return out


def monotone_pieces(arc: NDArray[np.float64]) -> list[Graph]:
    """Split an arc where its direction turns by more than 80 degrees, so each piece is a graph."""
    pieces, start = [], 0
    for i in range(1, len(arc) - 1):
        chord = arc[i] - arc[start]
        step = arc[i + 1] - arc[i]
        if np.linalg.norm(chord) > 0 and np.dot(chord, step) <= np.cos(np.radians(80)) * np.linalg.norm(chord) \
                * np.linalg.norm(step):
            pieces.append(arc[start:i + 1])
            start = i
    pieces.append(arc[start:])
    graphs = []
    for piece in pieces:
        if len(piece) < 2:
            continue
        try:
            graphs.append(Graph(piece))
        except AssertionError:  # not monotone along its chord: halve it until it is
            mid = len(piece) // 2
            graphs += monotone_pieces(piece[:mid + 1]) + monotone_pieces(piece[mid:])
    return graphs


# ── the curve ───────────────────────────────────────────────────────────────


class Curve:
    """
    A set of polylines as one point array and a piece id per point.

    A segment ``(i, i+1)`` exists iff ``ids[i] == ids[i+1]``.
    """

    def __init__(self, points: NDArray[np.float64], ids: NDArray[np.int64]) -> None:
        """Hold the arrays (``points`` is ``(N, 2)``)."""
        self.points, self.ids = points, ids

    @classmethod
    def from_polylines(cls, polylines: list[NDArray[np.float64]]) -> "Curve":
        """Concatenate polylines, one id each."""
        points = np.vstack(polylines)
        ids = np.concatenate([np.full(len(p), i, dtype=np.int64) for i, p in enumerate(polylines)])
        return cls(points, ids)

    def length(self) -> float:
        """Summed length of every segment."""
        same = self.ids[1:] == self.ids[:-1]
        return float(np.linalg.norm(np.diff(self.points, axis=0)[same], axis=1).sum())

    @property
    def pieces(self) -> int:
        """Number of polylines."""
        return int(np.count_nonzero(self.ids[1:] != self.ids[:-1]) + 1) if len(self.ids) else 0


def _map(fmap, points: NDArray[np.float64]) -> NDArray[np.float64]:
    """The map on an ``(N, 2)`` array (the library maps take ``(2, N)``)."""
    return fmap(points.T).T


def refine_and_map(fmap, points: NDArray[np.float64], ids: NDArray[np.int64],
                   tol: float) -> tuple[NDArray[np.float64], NDArray[np.int64]]:
    """
    Refine one chunk until every image segment's midpoint is within ``tol`` of its chord, then map it.

    Returns:
        The image points and their ids.
    """
    for _ in range(MAX_ROUNDS):
        same = ids[1:] == ids[:-1]
        a, b = points[:-1], points[1:]
        fa, fb = _map(fmap, a), _map(fmap, b)
        fm = _map(fmap, 0.5 * (a + b))
        deviation = np.linalg.norm(fm - 0.5 * (fa + fb), axis=1)
        split = np.nonzero(same & (deviation > tol))[0]
        del fa, fb, fm, deviation, a, b
        if not len(split):
            return _map(fmap, points), ids
        mids = 0.5 * (points[split] + points[split + 1])
        points = np.insert(points, split + 1, mids, axis=0)
        ids = np.insert(ids, split + 1, ids[split])
    raise RuntimeError(f"refinement did not converge in {MAX_ROUNDS} rounds")


def _distance_to_polyline(line: NDArray[np.float64], probe: NDArray[np.float64]) -> NDArray[np.float64]:
    """Distance from each probe point to a polyline (chunked over the probes)."""
    a, b = line[:-1], line[1:]
    ab = b - a
    denom = np.maximum((ab * ab).sum(axis=1), 1e-300)
    out = np.empty(len(probe))
    for lo in range(0, len(probe), 512):
        p = probe[lo:lo + 512, None, :]
        t = np.clip(((p - a) * ab).sum(axis=2) / denom, 0.0, 1.0)
        out[lo:lo + 512] = np.linalg.norm(p - (a + t[..., None] * ab), axis=2).min(axis=1)
    return out


def _contains(path: PolygonPath, probe: NDArray[np.float64],
              on_boundary: list[NDArray[np.float64]]) -> NDArray[np.bool_]:
    """
    Point in the zone polygon, a probe ON one of its unstable arcs counting inside.

    A piece lying along the zone's own unstable boundary arc (the image of the
    arc from the anchor) has its probe on the polygon edge, where the plain test
    answers by rounding; the blast counts the boundary inside
    (``ResonanceZone.contains_point``). Only the unstable arcs are tested: the
    curve is unstable manifold and can lie on no stable arc.
    """
    inside = path.contains_points(probe)
    outside = np.nonzero(~inside)[0]
    if len(outside):
        distance = np.min([_distance_to_polyline(arc, probe[outside]) for arc in on_boundary], axis=0)
        inside[outside[distance < BOUNDARY_TOL]] = True
    return inside


def cut_and_keep(curve: Curve, graphs: list[Graph], keep: PolygonPath, exclude: list[PolygonPath],
                 on_boundary: list[NDArray[np.float64]]) -> tuple[Curve, int]:
    """
    Cut the curve where it crosses a stable boundary arc and keep the pieces inside the region.

    Returns:
        The kept curve (pieces renumbered) and the number of boundary crossings.
    """
    points, ids = curve.points, curve.ids
    same = ids[1:] == ids[:-1]
    crossing = np.zeros(len(points) - 1, dtype=bool)
    where = np.zeros(len(points) - 1)
    for graph in graphs:
        side = graph.side(points)
        hit = same & (np.sign(side[:-1]) * np.sign(side[1:]) < 0)
        idx = np.nonzero(hit & ~crossing)[0]
        where[idx] = side[idx] / (side[idx] - side[idx + 1])
        crossing |= hit
    cuts = np.nonzero(crossing)[0]
    change = np.r_[0, (ids[1:] != ids[:-1]).astype(np.int64)]
    if len(cuts):
        x = points[cuts] + where[cuts, None] * (points[cuts + 1] - points[cuts])
        # Insert the crossing point twice: the end of one piece and the start of the next.
        points = np.insert(points, np.repeat(cuts + 1, 2), np.repeat(x, 2, axis=0), axis=0)
        change = np.insert(change, np.repeat(cuts + 1, 2), 0)
        change[cuts + 2 * np.arange(len(cuts)) + 2] = 1  # the second copy starts a new piece
    ids = np.cumsum(change)
    starts = np.r_[0, np.nonzero(ids[1:] != ids[:-1])[0] + 1]
    ends = np.r_[starts[1:], len(ids)]
    probe = np.where((ends - starts > 2)[:, None], points[(starts + ends) // 2],
                     0.5 * (points[starts] + points[ends - 1]))
    inside = _contains(keep, probe, on_boundary)
    for path in exclude:
        inside &= ~path.contains_points(probe)
    kept_point = np.repeat(inside, ends - starts)
    return Curve(points[kept_point], ids[kept_point]), int(len(cuts))


def iterate_curve(fmap, curve: Curve, tol: float, graphs, keep, exclude, on_boundary) -> tuple[Curve, int]:
    """One step: refine and map in chunks, then cut and keep."""
    out_points, out_ids = [], []
    n = len(curve.points)
    start = 0
    while start < n - 1:
        stop = min(start + CHUNK, n)
        # A chunk ends on a piece boundary or overlaps the next by one point.
        p, i = refine_and_map(fmap, curve.points[start:stop], curve.ids[start:stop], tol)
        if stop < n:
            p, i = p[:-1], i[:-1]  # the next chunk starts at stop-1 and owns that point
        out_points.append(p)
        out_ids.append(i)
        start = stop - 1 if stop < n else n
    image = Curve(np.vstack(out_points), np.concatenate(out_ids))
    del out_points, out_ids
    return cut_and_keep(image, graphs, keep, exclude, on_boundary)


# ── the runs ────────────────────────────────────────────────────────────────


def data_path(case: str, cutoff: float) -> Path:
    """Where one run's JSON goes."""
    return DATA_DIR / f"{case}_a{cutoff:g}.json"


def run(case: str, cutoff: float, tol: float, blast_gb: float, curve_gb: float, max_steps: int,
        max_blast_steps: int, blast_only: bool = False, diagnose: bool = False) -> None:
    """Both measurements of one case (only the blast with ``blast_only``, the curve kept from the file)."""
    path = data_path(case, cutoff)
    path.parent.mkdir(parents=True, exist_ok=True)
    record: dict = {"case": case, "k": entropy.K, "area_cutoff": cutoff, "tol": tol,
                    "blast": [], "curve": [], "blast_stopped": None, "curve_stopped": None}
    if blast_only and path.exists():
        record = json.loads(path.read_text())
        record.update({"blast": [], "blast_stopped": None})
    record["blast_diagnose"] = diagnose

    def save() -> None:
        path.write_text(json.dumps(record, indent=1))

    session, _, zones = entropy.build_case(case, cutoff)
    zone = max(zones, key=lambda z: z.area) if case == "nested" else zones[0]
    inner = [z for z in session.resonance_zones.values() if z.key is not zone.key and z.area < zone.area]
    keep = PolygonPath(zone.boundary_vertices)
    exclude = [PolygonPath(z.boundary_vertices) for z in inner]
    graphs = [g for z in [zone, *inner] for arc in stable_arcs(session, z) for g in monotone_pieces(arc)]
    on_boundary = unstable_arcs(session, zone)
    record["zone"] = {"area": zone.area, "excluded_areas": [z.area for z in inner], "stable_pieces": len(graphs)}
    fmap = henon_map(entropy.K, entropy.B)

    # The library, unpruned.
    start, previous_peak, initial = time.perf_counter(), entropy.resident_gb("VmHWM"), None
    for step in range(1, max_blast_steps + 1):
        peak = entropy.resident_gb("VmHWM")
        predicted = peak * peak / max(previous_peak, 1e-9)
        if step > 2 and predicted > blast_gb:
            record["blast_stopped"] = f"memory: step {step} predicted to peak at {predicted:.1f} GB"
            break
        previous_peak = peak
        result = session.blast_zone(zone, num_iterations=1, fixed_point=[zone.fixed_point], min_separation=None,
                                    diagnose=diagnose)
        if not result.steps:
            record["blast_stopped"] = f"frontier empty at step {step}"
            break
        kept = result.steps[-1].kept_interior
        if initial is None:
            initial = [b.get_point_array().copy() for b in result.steps[-1].interior_parents]
            record["initial_length"] = float(sum(entropy.arc_length(b) for b in result.steps[-1].interior_parents))
        record["blast"].append({
            "step": step, "length": float(sum(entropy.arc_length(b) for b in kept)), "bridges": len(kept),
            "points": int(result.points_created), "seconds": time.perf_counter() - start,
            "peak_gb": entropy.resident_gb("VmHWM"),
        })
        save()
        logger.info("%s blast %d: L=%.6g, %d bridges", case, step, record["blast"][-1]["length"], len(kept))
    del session, zones, zone, inner
    gc.collect()
    if blast_only:
        save()
        return

    # The independent curve.
    curve = Curve.from_polylines(initial)
    record["curve_initial_length"] = curve.length()
    start = time.perf_counter()
    for step in range(1, max_steps + 1):
        t0 = time.perf_counter()
        curve, crossings = iterate_curve(fmap, curve, tol, graphs, keep, exclude, on_boundary)
        n = len(curve.points)
        if not n:
            record["curve_stopped"] = f"curve left the zone at step {step}"
            break
        length = curve.length()
        record["curve"].append({
            "step": step, "length": length, "pieces": curve.pieces, "points": n, "crossings": crossings,
            "seconds": time.perf_counter() - start, "step_seconds": time.perf_counter() - t0,
            "rss_gb": entropy.resident_gb(),
        })
        save()
        logger.info("%s curve %d: L=%.6g, %d pieces, %d points, %.1f s", case, step, length, curve.pieces, n,
                    time.perf_counter() - t0)
        if len(record["curve"]) >= 2:
            growth = n / record["curve"][-2]["points"]
            # A step holds the curve, its image and one insertion copy, plus one chunk's temporaries.
            predicted = (n * growth * 3 * POINT_BYTES + CHUNK * 200) / 1e9
            if predicted > curve_gb:
                record["curve_stopped"] = f"memory: step {step + 1} predicted at {predicted:.1f} GB"
                break
    save()


def pruned_path(case: str, cutoff: float, separation: float) -> Path:
    """Where one pruned blast run's JSON goes (a subfolder: the figure globs the top folder for records)."""
    return DATA_DIR / "pruned" / f"{case}_a{cutoff:g}_s{separation:g}.json"


def run_pruned(case: str, cutoff: float, separation: float, blast_gb: float, max_blast_steps: int) -> None:
    """The library blast pruned at ``separation``, one step at a time; the step lengths, written after every step."""
    path = pruned_path(case, cutoff, separation)
    path.parent.mkdir(parents=True, exist_ok=True)
    record: dict = {"case": case, "area_cutoff": cutoff, "min_separation": separation, "blast": [], "stopped": None}
    session, _, zones = entropy.build_case(case, cutoff)
    zone = max(zones, key=lambda z: z.area) if case == "nested" else zones[0]
    start, previous_peak = time.perf_counter(), entropy.resident_gb("VmHWM")
    for step in range(1, max_blast_steps + 1):
        peak = entropy.resident_gb("VmHWM")
        predicted = peak * peak / max(previous_peak, 1e-9)
        if step > 2 and predicted > blast_gb:
            record["stopped"] = f"memory: step {step} predicted to peak at {predicted:.1f} GB"
            break
        previous_peak = peak
        result = session.blast_zone(zone, num_iterations=1, fixed_point=[zone.fixed_point],
                                    min_separation=separation, diagnose=False)
        if not result.steps:
            record["stopped"] = f"frontier empty at step {step}"
            break
        kept = result.steps[-1].kept_interior
        if "initial_length" not in record:
            record["initial_length"] = float(sum(entropy.arc_length(b) for b in result.steps[-1].interior_parents))
        record["blast"].append({
            "step": step, "length": float(sum(entropy.arc_length(b) for b in kept)), "bridges": len(kept),
            "points": int(result.points_created), "seconds": time.perf_counter() - start,
            "peak_gb": entropy.resident_gb("VmHWM"),
        })
        path.write_text(json.dumps(record, indent=1))
        logger.info("%s pruned %g step %d: L=%.6g, %d bridges", case, separation, step,
                    record["blast"][-1]["length"], len(kept))
    path.write_text(json.dumps(record, indent=1))


# ── the figure ──────────────────────────────────────────────────────────────


def fit_slope(steps: NDArray, lengths: NDArray, first: int) -> Optional[tuple[float, float, float]]:
    """``(slope, intercept, standard error)`` of ``ln L`` against the step, from ``first`` on."""
    mask = steps >= first
    if mask.sum() < 3:
        return None
    x, y = steps[mask], np.log(lengths[mask])
    (slope, intercept), cov = np.polyfit(x, y, 1, cov=True)
    return float(slope), float(intercept), float(np.sqrt(cov[0, 0]))


def topological_h(case: str) -> Optional[float]:
    """The sweep's last entropy at its finest cutoff."""
    best = entropy.finest(entropy.load_records()) if entropy.DATA_DIR.exists() else {}
    return best[case]["rows"][-1]["entropy"] if case in best else None


def pruned_series(case: str, cutoff: float) -> Optional[tuple[NDArray, NDArray]]:
    """The sweep's pruned frontier lengths at this cutoff, for contrast."""
    path = entropy.data_path(case, cutoff)
    if not path.exists():
        return None
    rows = [r for r in json.loads(path.read_text())["rows"] if r.get("frontier")]
    tangle = "p3" if case == "p3" else "p1"
    rows = [r for r in rows if tangle in r["frontier"] and r["frontier"][tangle]["length"] > 0]
    return (np.array([r["blasts"] for r in rows]), np.array([r["frontier"][tangle]["length"] for r in rows]))


def figure(fit_from: int) -> Path:
    """``ln L_n`` and its per-step rate per case: unpruned blast, independent curve, pruned sweep."""
    records = [json.loads(p.read_text()) for p in sorted(DATA_DIR.glob("*.json"))]
    records.sort(key=lambda r: (entropy.CASES.index(r["case"]), r["area_cutoff"]))
    fig, axes = plt.subplots(2, len(records), figsize=(6.2 * len(records), 9.5), squeeze=False)
    for col, record in enumerate(records):
        case, color = record["case"], entropy.CASE_COLORS[record["case"]]
        top, bottom = axes[0, col], axes[1, col]
        h = topological_h(case)
        series = []
        curve = record["curve"]
        if curve:
            x = np.array([r["step"] for r in curve])
            y = np.array([r["length"] for r in curve])
            series.append(("simple curve iteration: refine, map, cut at the zone's stable boundary", x, y, "-", "o", 1.0))
        # A blast step's kept bridges are the NEW layer of f^n(W) (earlier steps' bridges are never
        # re-iterated); the curve holds all of f^n(W), so the blasts are drawn as running totals.
        initial = record.get("initial_length", 0.0)
        blast = record["blast"]
        x_blast = [r["step"] for r in blast]
        if blast:
            x = np.array(x_blast)
            y = initial + np.cumsum([r["length"] for r in blast])
            series.append(("library blast, unpruned, min_separation None, running total", x, y, "none", "s", 1.0))
        pruned = pruned_series(case, record["area_cutoff"])
        if pruned is not None:
            series.append((r"library blast, pruned, min_separation $10^{-5}$, running total", pruned[0],
                           initial + np.cumsum(pruned[1]), ":", "x", 0.5))
        for extra in sorted(pruned_path(case, record["area_cutoff"], 1).parent.glob(
                f"{case}_a{record['area_cutoff']:g}_s*.json")):
            run_record = json.loads(extra.read_text())
            if run_record["blast"]:
                sep = run_record["min_separation"]
                exponent = int(round(np.log10(sep)))
                label = (rf"library blast, pruned, min_separation $10^{{{exponent}}}$, running total"
                         if np.isclose(sep, 10.0 ** exponent) else f"library blast, pruned, min_separation {sep:g}, running total")
                series.append((label, np.array([r["step"] for r in run_record["blast"]]),
                               run_record["initial_length"] + np.cumsum([r["length"] for r in run_record["blast"]]),
                               "--", "^", 0.7))
        for label, x, y, linestyle, marker, alpha in series:
            filled = marker != "s"
            top.plot(x, np.log(y), linestyle=linestyle, marker=marker, color=color if filled else "k",
                     markerfacecolor=None if filled else "none", markersize=5 if filled else 8, alpha=alpha,
                     label=label)
            rate = np.diff(np.log(y)) / np.diff(x)
            bottom.plot(x[1:], rate, linestyle=linestyle, marker=marker, color=color if filled else "k",
                        markerfacecolor=None if filled else "none", markersize=5 if filled else 8, alpha=alpha,
                        label=label)
        if curve:
            x = np.array([r["step"] for r in curve], dtype=float)
            y = np.array([r["length"] for r in curve])
            fit = fit_slope(x, y, fit_from)
            if fit is not None:
                slope, intercept, error = fit
                xs = np.array([fit_from, x[-1]])
                top.plot(xs, intercept + slope * xs, color="tab:orange", linewidth=2, alpha=0.8,
                         label=rf"Fit to the simple curve iteration: slope {slope:.4f} $\pm$ {error:.4f}")
                bottom.axhline(slope, color="tab:orange", linewidth=1.5, alpha=0.8)
        if h is not None:
            if curve:
                anchor_x, anchor_y = x[-1], np.log(y[-1])
                xs = np.array([1.0, anchor_x])
                top.plot(xs, anchor_y + h * (xs - anchor_x), color="tab:green", linestyle="--",
                         label=f"slope h = {h:.4f} (topology)")
            bottom.axhline(h, color="tab:green", linestyle="--", label=f"h = {h:.4f} (topology)")
        bottom.axhline(np.log(2), color="0.6", linestyle=":", linewidth=0.8, label=r"$\ln 2$")
        title = entropy.CASE_TITLES[case] + (" (the annulus: inner zone excluded)" if case == "nested" else "")
        top.set_title(f"{title}\ncutoff {record['area_cutoff']:g}, curve tol {record['tol']:g}", fontsize=10)
        top.set_ylabel(r"$\ln L_n$ (length kept in the zone)")
        bottom.set_ylabel(r"$\ln(L_n / L_{n-1})$")
        bottom.set_xlabel("step n")
        for ax in (top, bottom):
            ax.grid(alpha=0.3)
            ax.legend(fontsize=7, loc="lower right" if ax is top else "upper right")
        bottom.set_ylim(-0.2, 1.0)
        stops = [s for s in (record.get("blast_stopped"), record.get("curve_stopped")) if s]
        common = sorted(set(x_blast) & {r["step"] for r in curve})
        if common:
            last = common[-1]
            total = initial + sum(r["length"] for r in blast if r["step"] <= last)
            stops.insert(0, f"curve / blast total at step {last}: {curve[[r['step'] for r in curve].index(last)]['length'] / total:.4f}")
        if stops:
            top.text(0.02, 0.98, "\n".join(stops), transform=top.transAxes, ha="left", va="top", fontsize=7,
                     color="0.3")
    fig.suptitle(f"Hénon k={entropy.K:g}: length growth of the unstable curve inside each zone, unpruned and pruned")
    fig.tight_layout()
    return save_panel(fig, OUT_DIR, "length_growth")


def main() -> None:
    parser = argparse.ArgumentParser(description="Unpruned length growth inside the k=2.5 zones.")
    parser.add_argument("--cases", nargs="+", choices=entropy.CASES, default=list(entropy.CASES))
    parser.add_argument("--cutoff", type=float, default=1e-7)
    parser.add_argument("--tol", type=float, default=1e-5, help="chord tolerance of the independent curve")
    parser.add_argument("--blast-gb", type=float, default=24.0, help="peak budget of the library blast")
    parser.add_argument("--curve-gb", type=float, default=20.0, help="peak budget of the independent curve")
    parser.add_argument("--max-steps", type=int, default=60, help="last step of the independent curve")
    parser.add_argument("--max-blast-steps", type=int, default=30, help="last step of the library blast")
    parser.add_argument("--fit-from", type=int, default=4, help="first step of the slope fit")
    parser.add_argument("--blast-only", action="store_true", help="rerun the blast, keep the file's curve")
    parser.add_argument("--diagnose", action="store_true", help="the blast's same-stability count (~2.3x memory/step)")
    parser.add_argument("--pruned", type=float, nargs="+", default=None,
                        help="instead: the library blast pruned at these min_separations (data/pruned/)")
    parser.add_argument("--plot-only", action="store_true")
    args = parser.parse_args()
    configure_logging()
    plt.switch_backend("Agg")
    if args.pruned and not args.plot_only:
        for case in args.cases:
            for separation in args.pruned:
                run_pruned(case, args.cutoff, separation, args.blast_gb, args.max_blast_steps)
    elif not args.plot_only:
        for case in args.cases:
            run(case, args.cutoff, args.tol, args.blast_gb, args.curve_gb, args.max_steps, args.max_blast_steps,
                args.blast_only, args.diagnose)
    figure(args.fit_from)


if __name__ == "__main__":
    main()
