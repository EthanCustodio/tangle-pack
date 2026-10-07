"""Break the numerics on purpose: escape, deep blasts, deep growth, tangency.

Four sub-studies, one ``"study"`` key per config:

* ``escape`` -- the documented escape to infinity (ManifoldMachine Dev Notes):
  the k=2.1 period-3 orbit at ``area_cutoff = 1e-7``, its unstable arm oriented
  ``[0, 1]`` (stable ``[1, 1]``), grown ONE step at a time with a fixed count.
  The note records 478 -> ~9.9k -> ~1.5e6 points at iterations 8..10 and
  ``max|coord| ~ 2e26``; the RSS cap kills the run somewhere past that and the
  per-step checkpoint keeps the series. Mode ``arclength`` / ``predicate`` grow
  the same arm with the documented mitigation (a stop condition) instead.
* ``blasts`` -- k=2.8, period 3 and nested (outer zone) blasted one iteration
  at a time up to a cap, stopping at the first failure; then the topology
  stack (pseudoneighbors, holes, partition, symbolic dynamics) on the final
  state. The topology result is secondary and never fails the config.
* ``depth`` -- period 3 at 11..18 unstable steps and the k=10 inversion saddle
  at 4..8 (7 is documented as "does not finish"), through bridges and the
  iterate table, then the topology stack (secondary).
* ``tangency`` -- the period-1 ``b = 1`` saddle near the saddle-centre
  bifurcation ``k -> -1``, where the true separatrix splitting angle is
  exponentially small (``~ exp(-pi^2 / ln lambda)``, a heuristic reference
  only): crossings, near-parallel discards (counted by the harness warning
  templates) and the smallest crossing angle at the registered crossings.
  Probe (2026-10-07): every ``k`` from -0.9999 to 2.8 already has crossings
  when the arms are grown to ``20 * |x+ - x-|`` of canonical distance; there
  is no onset inside the existence range, so the scan is log-spaced in
  ``k + 1`` and asks how many of the near-``-1`` crossings are polygon
  artifacts.

Expected failures are data: an ``oom``/``timeout`` escape config, an
``exception`` or ``invariant_assert`` at a deep blast or depth.
"""

from __future__ import annotations

import math
import sys
import time
import traceback

import numpy as np
import psutil

import harness
from maps import henon_session, max_abs_coord, saddle_session, total_points

sys.path.insert(0, str(harness.REPO / "tests"))

EXPERIMENT = "e07_breaking"
TIMEOUT = 600
WORKERS = 4
RSS_CAP_GB = 20

P3_SEED = [[0, 1], [-1, 0], [-1, 1]]
#: The arm the ManifoldMachine Dev Note traced (it escapes); the cases grow the other one.
ESCAPE_ORIENTATION = {"unstable": np.array([0, 1]), "stable": np.array([1, 1])}


def configs(quick: bool) -> list[dict]:
    if quick:
        return [
            {"study": "escape", "mode": "fixed", "max_steps": 8},
            {"study": "escape", "mode": "arclength", "length": 0.5},
            {"study": "escape", "mode": "predicate", "radius": 2.5},
            {"study": "blasts", "case": "k28", "max_blasts": 2},
            {"study": "blasts", "case": "p3", "max_blasts": 1},
            {"study": "blasts", "case": "nested", "max_blasts": 1},
            {"study": "depth", "case": "p3", "unstable_steps": 11},
            {"study": "depth", "case": "inversion", "unstable_steps": 4},
            {"study": "tangency", "k": -0.99, "area_cutoff": 1e-7},
            {"study": "tangency", "k": 0.0, "area_cutoff": 1e-7},
        ]
    escape = [
        {"study": "escape", "mode": "fixed", "max_steps": 14, "_timeout": 600},
        {"study": "escape", "mode": "arclength", "length": 0.5},
        {"study": "escape", "mode": "arclength", "length": 3.0},
        {"study": "escape", "mode": "predicate", "radius": 2.5},
    ]
    blasts = [
        {"study": "blasts", "case": "k28", "max_blasts": 10, "_timeout": 900},
        {"study": "blasts", "case": "p3", "max_blasts": 10, "_timeout": 900},
        {"study": "blasts", "case": "nested", "max_blasts": 8, "_timeout": 900},
    ]
    depth = [
        {"study": "depth", "case": "p3", "unstable_steps": n, "_timeout": 480}
        for n in range(11, 19)
    ] + [
        {"study": "depth", "case": "inversion", "unstable_steps": n, "_timeout": 480}
        for n in range(4, 9)
    ]
    k_values = sorted({round(-1 + 10**e, 6) for e in np.linspace(-4, 0, 21)} | {0.5, 1.0, 1.5})
    tangency = [
        {"study": "tangency", "k": k, "area_cutoff": cutoff, "_timeout": 300}
        for cutoff in (1e-7, 1e-10, 1e-12)
        for k in k_values
    ]
    return escape + blasts + depth + tangency


def run(config: dict, rec: harness.Recorder) -> None:
    {"escape": run_escape, "blasts": run_blasts, "depth": run_depth, "tangency": run_tangency}[
        config["study"]
    ](config, rec)


# --------------------------------------------------------------------------- #
# Shared measurements
# --------------------------------------------------------------------------- #
def rss_mb() -> float:
    return psutil.Process().memory_info().rss / 2**20


def numerics_size(session) -> dict:
    """Crossings, bridges and points of a workbench, for one per-step snapshot."""
    wb = session.workbench
    bridges = wb.bridges
    return {
        "crossings": len(wb.intersection_registry),
        "bridges": len(bridges),
        "bridge_points": sum(len(b.get_point_array()) for b in bridges),
        "unstable_points": total_points(session, "unstable"),
        "stable_points": total_points(session, "stable"),
    }


def append_size(rec: harness.Recorder, session) -> None:
    for name, value in numerics_size(session).items():
        rec.append(name, value)


def run_topology(session, rec: harness.Recorder) -> None:
    """The topology stack on the final state: secondary, recorded, never raised."""
    t0 = time.perf_counter()
    try:
        session.compute_pseudoneighbors()
        session.punch_holes()
        session.partition_stable_manifold()
        dynamics = session.symbolic_dynamics()
        rec.metric("topology_outcome", "ok")
        rec.metric("is_reliable", bool(dynamics.is_reliable))
        rec.metric("n_classes", len(dynamics.classes))
        rec.metric("n_unresolved", len(dynamics.unresolved))
        rec.metric("n_virtual", len(dynamics.virtual_classes))
    except Exception as exc:  # noqa: BLE001 - the failure is the datum
        kind = "invariant_assert" if isinstance(exc, AssertionError) else "exception"
        rec.metric("topology_outcome", kind)
        rec.metric("topology_error", f"{type(exc).__name__}: {exc}"[:500])
        rec.metric("topology_traceback", traceback.format_exc()[-2000:])
    rec.metric("topology_s", time.perf_counter() - t0)


# --------------------------------------------------------------------------- #
# escape
# --------------------------------------------------------------------------- #
def escape_session():
    session = henon_session(2.1, 1, area_cutoff=1e-7)
    try:
        fp = session.construct_fixed_point(P3_SEED)
        seed_route = "direct"
    except ValueError:
        fp, seed_route = continue_period3(session, 2.1), "continuation"
    session.orient_eigenvectors(fp, ESCAPE_ORIENTATION)
    session.initialize_both_manifolds(fp)
    return session, fp, seed_route


def continue_period3(session, k_target: float, dk: float = 0.01):
    """Track the exact k=2 period-3 orbit to ``k_target`` in small k steps."""
    orbit = np.array(P3_SEED, dtype=float)
    for k in np.arange(2.0 + dk, k_target + dk / 2, dk):
        orbit = henon_session(float(k), 1).workbench._fp_solver.compute_fixed_point(orbit)
    return session.construct_fixed_point(orbit)


def arm_snapshot(session, fp) -> dict:
    manifold = session.workbench.manifolds[(fp, "unstable", 0, 0)]
    points = manifold.get_point_array()
    segments = np.linalg.norm(np.diff(points, axis=0), axis=1)
    return {
        "points": total_points(session, "unstable"),
        "max_abs": max_abs_coord(session),
        "tail_cdist": float(manifold.tail.cdist),
        "min_segment": float(segments.min()) if len(segments) else None,
        "rss_mb": rss_mb(),
    }


def run_escape(config: dict, rec: harness.Recorder) -> None:
    with rec.stage("setup"):
        session, fp, seed_route = escape_session()
    rec.metric("seed_route", seed_route)
    rec.metric("orbit", fp.coordinates)
    wb = session.workbench

    if config["mode"] == "fixed":
        for step in range(1, config["max_steps"] + 1):
            t0 = time.perf_counter()
            session.grow_n_times(fp, "unstable", num_iterations=1)
            rec.append("step_s", time.perf_counter() - t0)
            rec.append("step", step)
            for name, value in arm_snapshot(session, fp).items():
                rec.append(name, value)
            rec.metric("steps_done", step)
            rec.checkpoint()
        return

    if config["mode"] == "arclength":
        with rec.stage("grow"):
            wb.grow_until_arclength(fp, "unstable", length=config["length"], max_iterations=14)
    else:
        radius = config["radius"]

        def left_box(workbench) -> bool:
            # The predicate runs once per round, after the crossings are recomputed.
            snapshot = arm_snapshot(session, fp)
            for name, value in snapshot.items():
                rec.append(name, value)
            rec.checkpoint()
            return snapshot["max_abs"] > radius

        with rec.stage("grow"):
            rounds = wb.grow_until(fp, left_box, grow=("unstable",), max_iterations=14)
        rec.metric("rounds", rounds)
    for name, value in arm_snapshot(session, fp).items():
        rec.metric(f"final_{name}", value)


# --------------------------------------------------------------------------- #
# blasts
# --------------------------------------------------------------------------- #
def blast_setup(case: str):
    """Session, the fixed points to repin, their pips, and the zone to blast."""
    from cases import build_k28
    from tanglepack.examples import henon_cases as hc

    if case == "k28":
        built = build_k28(blasts=0, through="zones")
        return built.session, built.fixed_points, built.pips, built.zones[0]

    session = hc._session()
    fp3 = hc._grow_period3(session, 13, 9, 1e-7)
    if case == "p3":
        hc._intersect_and_bridge(session, [fp3])
        session.classify_strong_pips()
        pip = session.trellis(fp3).strong_pip
        zone = session.resonance_zone(pip)
        hc._repin(session, [fp3], [pip])
        return session, [fp3], [pip], zone

    fp1 = hc._grow_period1(session, 11, 1e-7)
    fixed_points = [fp1, fp3]
    hc._intersect_and_bridge(session, [fp3, fp1])
    session.classify_strong_pips()
    pips = [session.trellis(fp).strong_pip for fp in fixed_points]
    session.add_resonance_zones(pips)
    outer = max(session.resonance_zones.values(), key=lambda zone: zone.area)
    hc._repin(session, fixed_points, pips)
    return session, fixed_points, pips, outer


def run_blasts(config: dict, rec: harness.Recorder) -> None:
    case = config["case"]
    with rec.stage("setup"):
        session, fixed_points, pips, zone = blast_setup(case)
    # k28 blasts with min_separation 1e-5 (build_k28), the nested outer zone with 1e-4.
    min_separation = 1e-4 if case == "nested" else 1e-5
    append_size(rec, session)
    rec.append("blast_s", 0.0)
    rec.append("repin_s", 0.0)
    rec.append("interior_bridges", 0)
    rec.checkpoint()

    for blast in range(1, config["max_blasts"] + 1):
        rec.metric("attempted_blast", blast)
        t0 = time.perf_counter()
        result = session.blast_zone(
            zone, num_iterations=1, fixed_point=[zone.fixed_point], min_separation=min_separation
        )
        t1 = time.perf_counter()
        session.classify_strong_pips()
        for fp, pip in zip(fixed_points, pips):
            session.set_strong_pip(fp, pip)
        rec.append("blast_s", t1 - t0)
        rec.append("repin_s", time.perf_counter() - t1)
        rec.append("interior_bridges", len(result.all_interior_bridges()))
        append_size(rec, session)
        rec.metric("blasts_done", blast)
        rec.checkpoint()

    with rec.stage("topology"):
        run_topology(session, rec)


# --------------------------------------------------------------------------- #
# depth
# --------------------------------------------------------------------------- #
def run_depth(config: dict, rec: harness.Recorder) -> None:
    from cases import K10_INVERSION_SEED, PERIOD3_ORIENTATION

    steps = config["unstable_steps"]
    if config["case"] == "p3":
        session = henon_session(2, 1, area_cutoff=1e-7)
        fp = session.construct_fixed_point(P3_SEED)
        session.orient_eigenvectors(fp, {s: np.array(v) for s, v in PERIOD3_ORIENTATION.items()})
        stable_steps = 9
    else:
        session = henon_session(10, 1)
        fp = session.construct_fixed_point(list(K10_INVERSION_SEED))
        stable_steps = 5
    session.initialize_both_manifolds(fp)

    for step in range(1, steps + 1):
        t0 = time.perf_counter()
        session.grow_n_times(fp, "unstable", num_iterations=1)
        rec.append("step_s", time.perf_counter() - t0)
        rec.append("unstable_points", total_points(session, "unstable"))
        rec.append("max_abs", max_abs_coord(session))
        rec.append("rss_mb", rss_mb())
        rec.metric("steps_done", step)
        rec.checkpoint()
    with rec.stage("grow_stable"):
        session.grow_n_times(fp, "stable", num_iterations=stable_steps)
    rec.checkpoint()
    with rec.stage("intersections"):
        session.compute_intersections([fp])
    rec.checkpoint()
    with rec.stage("bridges"):
        session.trim_stable_manifolds(fp)
        session.create_bridges(fp)
        session.infer_iterate_table()
    for name, value in numerics_size(session).items():
        rec.metric(name, value)
    rec.metric("numerics_done", True)
    rec.checkpoint()
    with rec.stage("topology"):
        session.classify_strong_pips()
        if config["case"] == "p3":
            # build_period3: trim at the default pip's zone, then restore the pip.
            pip = session.strong_pip(fp)
            session.resonance_zone(pip)
            session.classify_strong_pips()
            session.set_strong_pip(fp, pip)
        run_topology(session, rec)


# --------------------------------------------------------------------------- #
# tangency
# --------------------------------------------------------------------------- #
def crossing_sines(session) -> list[float]:
    """|sin| of the angle between the two curves at every non-anchor crossing.

    The unstable direction is the crossing's own ``unstable_segment``; the stable
    one is the nearest segment of the stable polyline it was registered on.
    """
    wb = session.workbench
    stable_polylines = {}
    sines = []
    for _id, ix in wb.intersection_registry:
        if ix.label == "anchor" or ix.unstable_segment is None:
            continue
        key = ix.manifold_b_key
        if key not in stable_polylines:
            stable_polylines[key] = wb.manifolds[key].get_point_array()
        stable = stable_polylines[key]
        point = np.asarray(ix.coords)
        a, b = stable[:-1], stable[1:]
        ab = b - a
        t = np.clip(np.einsum("ij,ij->i", point - a, ab) / np.einsum("ij,ij->i", ab, ab), 0, 1)
        nearest = int(np.argmin(np.linalg.norm(a + t[:, None] * ab - point, axis=1)))
        u = ix.unstable_segment[1].get_point() - ix.unstable_segment[0].get_point()
        s = ab[nearest]
        sines.append(abs(u[0] * s[1] - u[1] * s[0]) / (np.linalg.norm(u) * np.linalg.norm(s)))
    return sines


def run_tangency(config: dict, rec: harness.Recorder) -> None:
    k = config["k"]
    with rec.stage("setup"):
        session, fp = saddle_session(k, 1, area_cutoff=config["area_cutoff"])
    lam = float(np.ravel(fp.unstable_eigenvalues[0])[0])
    separation = math.sqrt(4 + 4 * k)
    rec.metric("lambda", lam)
    rec.metric("separation", separation)
    rec.metric("splitting_reference", math.exp(-math.pi**2 / math.log(lam)))

    # Grow each arm to 20 fixed-point separations of canonical distance (the probe's
    # length), a stop condition rather than a step count: lambda -> 1 as k -> -1.
    wb = session.workbench
    for stability in ("unstable", "stable"):
        steps = 0
        with rec.stage(f"grow_{stability}"):
            while wb.manifolds[(fp, stability, 0, 0)].tail.cdist < 20 * separation and steps < 300:
                session.grow_n_times(fp, stability, num_iterations=1)
                steps += 1
        rec.metric(f"{stability}_steps", steps)
    rec.metric("unstable_points", total_points(session, "unstable"))
    rec.metric("stable_points", total_points(session, "stable"))
    rec.metric("max_abs", max_abs_coord(session))
    rec.checkpoint()

    with rec.stage("intersections"):
        session.compute_intersections([fp])
    sines = crossing_sines(session)
    rec.metric("crossings", len(wb.intersection_registry))
    rec.metric("non_anchor_crossings", len(sines))
    if sines:
        rec.metric("min_sin", min(sines))
        rec.metric("median_sin", float(np.median(sines)))
        rec.metric("n_sin_below_1e-6", sum(s < 1e-6 for s in sines))
        rec.metric("sines", sorted(sines)[:2000])


if __name__ == "__main__":
    harness.main(globals())
