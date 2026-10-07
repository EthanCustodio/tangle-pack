"""Accuracy versus cost of the refinement knob ``area_cutoff``.

Sweeps ``ManifoldMachine.area_cutoff`` over 1e-3 .. 1e-9 on two cases, each
at its canonical depth and at a deeper one (calibrated 2026-10-07: the
canonical depths finish in about a second even at 1e-9, so the deep variant
is what makes the cost curve visible; its 1e-9 run takes 1.5-3 min):

* ``k28``: the k=2.8 saddle, 10 (canonical) or 12 unstable steps, the stable
  manifold to turnaround (the arm escapes at 12 steps: 1.9M points at 1e-9);
* ``p3``: the period-3 orbit of the nested map (2, 1), 13/9 (canonical) or
  17/11 unstable/stable steps (2.7M points at 1e-9).

Each run records cost (stage times, per-step points and seconds, RSS) and
everything the figure-maker needs to measure convergence against the finest
run of the same case: crossing coordinates and cdists, bridge ids and their
lobe areas. Inside the run it also records a reference-free error proxy, the
INVARIANCE error: a point sampled inside a segment of the polyline is mapped
one step (the inverse map on the stable side) and its distance to the
polyline is measured (``maps.invariance_error``).

Record schema
-------------
stages: ``setup``, ``grow_unstable``, ``grow_stable``, ``intersections``,
``trim``, ``bridges``, ``iterate_table``, ``invariance``, ``lobes``.
metrics: ``case``, ``area_cutoff``, ``unstable_steps``, ``stable_steps``
(None = turnaround), ``points_unstable``, ``points_stable``,
``manifold_sizes``, ``max_abs_coord``, ``n_crossings`` (non-synthetic),
``n_anchors``, ``crossings`` (rows ``[id, x, y, unstable_cdist,
stable_cdist, unstable_branch, stable_branch]``, branch = ``"o.b"``),
``points_unstable_after_trim``, ``points_stable_after_trim``,
``n_bridges``, ``n_partial_bridges``, ``bridges`` (rows ``[first, second,
n_points, lobe_area]``, lobe area signed, None when not closable),
``n_iterate_links``; per manifold, tag ``inv_<stability>_<o.b>``:
``<tag>_{n,max,mean,p50,p90,p99}`` (invariance error), ``<tag>_eligible``
(sampled segments), ``<tag>_beyond`` (images past the tail, should be 0),
``<tag>_on_curve`` (``[vertices mapping onto a vertex, vertices]``).
series: ``unstable_points`` / ``unstable_step_s`` per unstable step,
``stable_points`` / ``stable_step_s`` per stable step,
``<tag>_sample`` (rows ``[cdist, err]``, up to 300 per manifold).
"""
from __future__ import annotations

import time

import numpy as np

import harness
from maps import (
    PERIOD1_ORIENTATION,
    curve_arrays,
    error_summary,
    henon_session,
    henon_fixed_points,
    invariance_error,
    lobe_polygon,
    manifold_sizes,
    max_abs_coord,
    total_points,
)
from tanglepack.examples import henon_cases
from tanglepack.examples.henon import HENON_P3, saddle_guesses
from tanglepack.numerics.geometry import signed_polygon_area

EXPERIMENT = "e04_area_cutoff"
TIMEOUT = 900
WORKERS = 1  # a timing experiment: one config at a time
RSS_CAP_GB = 16

CUTOFFS = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9]
#: (case, unstable steps, stable steps or None for turnaround)
DEPTHS = [("k28", 10, None), ("k28", 12, None), ("p3", 13, 9), ("p3", 17, 11)]


def configs(quick: bool) -> list[dict]:
    if quick:
        return [
            {"case": "k28", "unstable_steps": 10, "stable_steps": None, "area_cutoff": c}
            for c in (1e-4, 1e-6)
        ] + [{"case": "p3", "unstable_steps": 13, "stable_steps": 9, "area_cutoff": 1e-4}]
    return [
        {"case": case, "unstable_steps": nu, "stable_steps": ns, "area_cutoff": c}
        for case, nu, ns in DEPTHS
        for c in CUTOFFS
    ]


def _branch(key) -> str:
    return f"{key[2]}.{key[3]}"


def _grow(session, fp, stability: str, steps, rec) -> None:
    """Grow one step at a time so the cost per step is a series."""
    if steps is None:
        with rec.stage(f"grow_{stability}"):
            session.grow_until_turnaround(fp, stability)
        rec.append(f"{stability}_points", total_points(session, stability))
        return
    for _ in range(steps):
        t0 = time.perf_counter()
        with rec.stage(f"grow_{stability}"):
            session.grow_n_times(fp, stability, num_iterations=1)
        rec.append(f"{stability}_step_s", time.perf_counter() - t0)
        rec.append(f"{stability}_points", total_points(session, stability))


def run(config: dict, rec) -> None:
    case, cutoff = config["case"], config["area_cutoff"]
    with rec.stage("setup"):
        if case == "k28":
            session = henon_session(2.8, area_cutoff=cutoff)
            fp = session.construct_fixed_point(henon_fixed_points(2.8)["saddle"])
            session.orient_eigenvectors(fp, PERIOD1_ORIENTATION)
        else:
            session = henon_session(*HENON_P3, area_cutoff=cutoff)
            fp = session.construct_fixed_point(saddle_guesses(*HENON_P3)["period_3"])
            session.orient_eigenvectors(fp, henon_cases.PERIOD3_ORIENTATION)
        session.initialize_both_manifolds(fp)
    for name in ("case", "area_cutoff", "unstable_steps", "stable_steps"):
        rec.metric(name, config[name])

    _grow(session, fp, "unstable", config["unstable_steps"], rec)
    rec.checkpoint()
    _grow(session, fp, "stable", config["stable_steps"], rec)
    rec.metric("points_unstable", total_points(session, "unstable"))
    rec.metric("points_stable", total_points(session, "stable"))
    rec.metric("manifold_sizes", manifold_sizes(session))
    rec.metric("max_abs_coord", max_abs_coord(session))
    rec.checkpoint()

    # Invariance before the trim: the stable side then still holds its whole grown length.
    system = session.workbench.dynamical_system
    with rec.stage("invariance"):
        for key, manifold in session.workbench.manifolds.items():
            stability = key[1]
            forward = stability == "unstable"
            image_key = fp.advance_key(key, 1 if forward else -1)
            probe = invariance_error(
                manifold,
                system.map if forward else system.map_inv,
                stability,
                target=session.workbench.manifolds.get(image_key),
            )
            tag = f"inv_{stability}_{_branch(key)}"
            rec.metrics.update(error_summary(probe["err"], prefix=f"{tag}_"))
            rec.metric(f"{tag}_eligible", probe["eligible"])
            rec.metric(f"{tag}_beyond", int(probe["beyond"].sum()))
            rec.metric(f"{tag}_on_curve", [probe["on_curve"], probe["vertices"]])
            keep = np.linspace(0, len(probe["err"]) - 1, min(300, len(probe["err"]))).astype(int)
            rec.series[f"{tag}_sample"] = np.column_stack(
                [probe["cdist"][keep], probe["err"][keep]]
            ).tolist() if len(keep) else []
    rec.checkpoint()

    with rec.stage("intersections"):
        session.compute_intersections([fp])
    registry = session.workbench.intersection_registry
    crossings = [
        [iid, *ix.coords, ix.unstable_cdist, ix.stable_cdist,
         _branch(ix.manifold_a_key), _branch(ix.manifold_b_key)]
        for iid, ix in registry
        if not ix.is_synthetic
    ]
    rec.metric("n_crossings", len(crossings))
    rec.metric("n_anchors", len(registry) - len(crossings))
    rec.metric("crossings", crossings)
    rec.checkpoint()

    with rec.stage("trim"):
        session.trim_stable_manifolds(fp)
    rec.metric("points_unstable_after_trim", total_points(session, "unstable"))
    rec.metric("points_stable_after_trim", total_points(session, "stable"))
    with rec.stage("bridges"):
        session.create_bridges(fp)
    with rec.stage("iterate_table"):
        session.infer_iterate_table()
    table = session.workbench.intersection_registry.iterate_table
    rec.metric("n_iterate_links", sum(1 for _ in table.items()))
    rec.checkpoint()

    with rec.stage("lobes"):
        curves = curve_arrays(session)
        rows = []
        for bridge in session.workbench.bridges:
            if bridge.id is None:
                continue
            lobe = lobe_polygon(session, curves, *bridge.id)
            area = None if lobe is None else signed_polygon_area(lobe)
            rows.append([*bridge.id, len(bridge.get_point_array()), area])
    rec.metric("n_bridges", len(rows))
    rec.metric("n_partial_bridges", len(session.workbench.bridges) - len(rows))
    rec.metric("bridges", rows)


if __name__ == "__main__":
    harness.main(globals())
