"""Robustness of the numerics pipeline over the Hénon parameter plane.

One fixed recipe (``area_cutoff = 1e-6``; the regular saddle: 8 unstable
steps, the stable manifold to turnaround; the inversion saddle: 6 unstable /
5 stable steps, as ``cases.build_inversion``) run as far as it goes on

* ``b = 1``, the regular saddle, ``k`` geometric over 0.3 .. 15 (dense at
  small ``k``);
* ``b = 1``, the OTHER period-1 point: elliptic for ``k < 3`` (construction
  is expected to raise), the inversion saddle above;
* ``b = -1`` (orientation reversing: ``set_k_value`` is expected to raise);
* dissipative ``b in {0.9, 0.5, 0.3}`` (out of scope: whatever happens is
  recorded).

Each stage is checkpointed, so a failure, a timeout or an OOM still leaves
``stage_reached`` and every metric of the stages before it.

Record schema
-------------
stages: ``fixed_point``, ``eigen``, ``seeded``, ``grow_unstable``,
``grow_stable``, ``intersections``, ``bridges``, ``trellis_pips``.
metrics: ``k``, ``b``, ``which``, ``stage_reached`` (last stage completed:
``none``, ``fixed_point``, ``eigen``, ``seeded``, ``grown_unstable``,
``grown``, ``intersected``, ``bridges``, ``trellis_pips``),
``fp_x``, ``fp_residual``, ``lambda_u``, ``lambda_s``, ``k_value``,
``inversion``, ``points_unstable_grown``, ``points_unstable``,
``points_stable``, ``max_abs_coord``, ``n_crossings``, ``n_anchors``,
``n_bridges``, ``n_partial_bridges``, ``n_iterate_links``,
``strong_pip`` (id or None), ``n_pip_candidates``,
``trellis_bridges``.
series: ``unstable_points`` per unstable step.
"""
from __future__ import annotations

import numpy as np

import harness
from maps import (
    PERIOD1_ORIENTATION,
    henon_fixed_points,
    henon_session,
    max_abs_coord,
    total_points,
)

EXPERIMENT = "e06_parameter_sweep"
TIMEOUT = 240
WORKERS = 8
RSS_CAP_GB = 4

AREA_CUTOFF = 1e-6


def configs(quick: bool) -> list[dict]:
    if quick:
        return [
            {"k": 0.3, "b": 1.0, "which": "saddle"},
            {"k": 2.8, "b": 1.0, "which": "saddle"},
            {"k": 10.0, "b": 1.0, "which": "other"},
            {"k": 4.0, "b": -1.0, "which": "saddle"},
        ]
    saddle = np.round(np.geomspace(0.3, 15, 60), 4)
    other = np.round(np.geomspace(1, 15, 15), 4)
    reversing = [0.5, 1, 2, 4, 6, 8, 10, 15]
    dissipative = [1, 2.8, 5, 10, 15]
    return (
        [{"k": float(k), "b": 1.0, "which": "saddle"} for k in saddle]
        + [{"k": float(k), "b": 1.0, "which": "other"} for k in other]
        + [{"k": float(k), "b": -1.0, "which": "saddle"} for k in reversing]
        + [{"k": float(k), "b": b, "which": "saddle"} for b in (0.9, 0.5, 0.3) for k in dissipative]
    )


def run(config: dict, rec) -> None:
    k, b, which = config["k"], config["b"], config["which"]
    for name in ("k", "b", "which"):
        rec.metric(name, config[name])
    rec.metric("stage_reached", "none")

    def reached(stage: str) -> None:
        rec.metric("stage_reached", stage)
        rec.checkpoint()

    session = henon_session(k, b, area_cutoff=AREA_CUTOFF)
    seed = henon_fixed_points(k, b)[which]
    if seed is None:
        raise ValueError(f"no real period-1 fixed point for k={k}, b={b}")
    with rec.stage("fixed_point"):
        fp = session.construct_fixed_point(seed)
    x = np.ravel(fp.coordinates[0]).astype(float)
    image = np.asarray(session.workbench.dynamical_system.map(x.reshape(2, 1))).ravel()
    rec.metric("fp_x", x.tolist())
    rec.metric("fp_residual", float(np.linalg.norm(image - x)))
    rec.metric("lambda_u", float(np.ravel(fp.unstable_eigenvalues[0])[0]))
    rec.metric("lambda_s", float(np.ravel(fp.stable_eigenvalues[0])[0]))
    rec.metric("k_value", fp.k_value)
    inversion = bool(fp.check_inversion())
    rec.metric("inversion", inversion)
    reached("fixed_point")

    # The inversion saddle's second branch is the image of the first, so it is
    # not oriented (cases.build_inversion).
    with rec.stage("eigen"):
        if not inversion:
            session.orient_eigenvectors(fp, PERIOD1_ORIENTATION)
    reached("eigen")
    with rec.stage("seeded"):
        session.initialize_both_manifolds(fp)
    reached("seeded")

    for _ in range(6 if inversion else 8):
        with rec.stage("grow_unstable"):
            session.grow_n_times(fp, "unstable", num_iterations=1)
        rec.append("unstable_points", total_points(session, "unstable"))
    rec.metric("points_unstable_grown", total_points(session, "unstable"))
    reached("grown_unstable")
    with rec.stage("grow_stable"):
        if inversion:
            session.grow_n_times(fp, "stable", num_iterations=5)
        else:
            session.grow_until_turnaround(fp, "stable")
    rec.metric("points_unstable", total_points(session, "unstable"))
    rec.metric("points_stable", total_points(session, "stable"))
    rec.metric("max_abs_coord", max_abs_coord(session))
    reached("grown")

    with rec.stage("intersections"):
        session.compute_intersections([fp])
    registry = session.workbench.intersection_registry  # compute_intersections replaces it
    crossings = sum(1 for _, ix in registry if not ix.is_synthetic)
    rec.metric("n_crossings", crossings)
    rec.metric("n_anchors", len(registry) - crossings)
    reached("intersected")

    with rec.stage("bridges"):
        session.trim_stable_manifolds(fp)
        session.create_bridges(fp)
        session.infer_iterate_table()
    bridges = session.workbench.bridges
    rec.metric("n_bridges", sum(1 for br in bridges if br.id is not None))
    rec.metric("n_partial_bridges", sum(1 for br in bridges if br.id is None))
    table = session.workbench.intersection_registry.iterate_table
    rec.metric("n_iterate_links", sum(1 for _ in table.items()))
    reached("bridges")

    with rec.stage("trellis_pips"):
        session.classify_strong_pips()
        trellis = session.trellis(fp)
    rec.metric("strong_pip", trellis.strong_pip)
    rec.metric("n_pip_candidates", len(trellis.strong_pip_candidates))
    rec.metric("trellis_bridges", len(trellis.bridges))
    reached("trellis_pips")


if __name__ == "__main__":
    harness.main(globals())
