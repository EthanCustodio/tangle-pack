"""Physical and numerical invariants measured on the standard cases.

One config per frozen case of ``tests/cases.py`` (k10, k28 with 0/1/2 blasts,
p3, nested, inversion; each built in full, through the partition) plus one
CPU-versus-GPU comparison of the k28 growth. Per case:

(a) the fixed point: ``|f^p(x*) - x*|``, ``|lambda_u lambda_s| - 1`` (0 for an
    area-preserving map), ``|J v - lambda v|`` for both eigenvectors, and the
    stored full-cycle Jacobian against the analytic product along the orbit;
(b) manifold invariance (``maps.invariance_error``: a point inside a polyline
    segment mapped one step -- the inverse map on the stable side -- onto the
    image branch, distance to that polyline), summary, raw samples, error vs
    cdist;
(c) iterate consistency: ``|f^n(x_id) - x_target|`` for every registered
    iterate link (non-synthetic ends), by ``n`` and against the source's
    backward depth in the table;
(d) area preservation: the signed lobe area of every bridge (bridge + stable
    arc between its ends) against its image lobe (unstable arc between the
    ``+1`` iterates of its ends + stable arc), and the arrangement's own
    ``image_of`` (area-verified to 2 %) per closed region;
(e) cdist monotonicity: non-increasing and tied steps along every manifold
    and every bridge;
(f) WARNING counts (recorded by the harness).

Record schema
-------------
stages: ``build``, ``fixed_point``, ``invariance``, ``iterates``, ``lobes``,
``arrangement``, ``cdist`` (gpu config: ``cpu_grow``, ``gpu_grow``,
``compare``).
metrics (case configs):
``case``, ``points_unstable``, ``points_stable``, ``n_crossings``,
``n_bridges``;
(a) ``fp`` = list per fixed point of ``{label, period, k_value, orbit:
[{x, fp_residual, lambda_u, lambda_s, det_minus_1, evec_res_u, evec_res_s,
jac_vs_analytic}]}``;
(b) ``inv_<label>:<stability>/<o.b>_{n,max,mean,p50,p90,p99}``,
``..._eligible``, ``..._beyond``, ``..._on_curve`` ([vertices mapping onto a
vertex of the image branch, vertices]), ``..._missing_image`` (image branch
not registered);
(c) ``iter_n<n>_{n,max,...}``, ``iter_links`` (count), ``iter_skipped``;
(d) ``lobe_pairs`` (count compared), ``lobe_ratio_{n,max,...}`` (of
``|A_image / A_source - 1|``), ``lobe_sign_flips``, ``lobe_unclosable``;
``arr_regions``, ``arr_image_found``, ``arr_image_missing``,
``arr_ratio_{...}`` (of ``|A_image / A - 1|`` for the found ones);
(e) ``cdist_manifolds`` / ``cdist_bridges`` = ``{decreasing, ties, steps,
worst_drop}`` totals and ``cdist_per_manifold`` (``{key: [decreasing,
ties, steps]}``).
series: ``inv_<label>:<stability>/<o.b>`` (rows ``[cdist, err]``, <= 300),
``iter_links`` (rows ``[n, backward_depth, unstable_cdist, err]``, <= 2000),
``lobes`` (rows ``[first, second, area, image_area]``), ``arr_regions``
(rows ``[area, image_area or None]``, <= 500).
metrics (gpu config ``case = "gpu_k28"``): ``gpu_available``,
``cpu_points``/``gpu_points`` per manifold, ``same_lengths``,
``max_abs_diff`` (equal lengths) or ``hausdorff`` (polyline distance both
ways), ``cpu_crossings``, ``gpu_crossings``, ``crossing_max_dev``,
``crossing_unmatched``, ``cpu_bridges``, ``gpu_bridges``.
"""
from __future__ import annotations

import sys

import numpy as np

import harness
from maps import (
    curve_arrays,
    distance_to_polyline,
    error_summary,
    henon_session,
    henon_fixed_points,
    invariance_error,
    lobe_polygon,
    PERIOD1_ORIENTATION,
    total_points,
)
from tanglepack.examples.henon import henon_jacobian
from tanglepack.numerics.geometry import signed_polygon_area

sys.path.insert(0, str(harness.REPO / "tests"))
import cases  # noqa: E402

EXPERIMENT = "e05_invariants"
TIMEOUT = 600
WORKERS = 4
RSS_CAP_GB = 12

#: case -> (builder, Hénon parameters)
CASES = {
    "k10": (cases.build_k10, (10, 1)),
    "k28_unblasted": (lambda: cases.build_k28(blasts=0), (2.8, 1)),
    "k28_one_blast": (lambda: cases.build_k28(blasts=1), (2.8, 1)),
    "k28_two_blasts": (lambda: cases.build_k28(blasts=2), (2.8, 1)),
    "p3": (cases.build_period3, (2, 1)),
    "nested": (cases.build_nested, (2, 1)),
    "inversion": (cases.build_inversion, (10, 1)),
}


def configs(quick: bool) -> list[dict]:
    if quick:
        return [{"case": "k28_unblasted"}, {"case": "inversion"}, {"case": "gpu_k28", "steps": 6}]
    return [{"case": name} for name in CASES] + [{"case": "gpu_k28", "steps": 10}]


def _key_tag(key) -> str:
    fp, stability, orbit, branch = key
    return f"{fp.label or 'fp'}:{stability}/{orbit}.{branch}"


def _apply(system, x: np.ndarray, n: int) -> np.ndarray:
    """``f^n`` (``n < 0``: the inverse) of a ``(2, N)`` batch."""
    step = system.map if n > 0 else system.map_inv
    for _ in range(abs(n)):
        x = np.asarray(step(x), dtype=float)
    return x


def fixed_point_checks(system, fp, hk: float, hb: float) -> dict:
    """Residual, eigenvalue product and eigenvector residuals of one orbit."""
    jac = henon_jacobian(hk, hb)
    orbit = []
    for i in range(fp.period):
        x = np.ravel(fp.coordinates[i]).astype(float)
        image = _apply(system, x.reshape(2, 1), fp.period).ravel()
        analytic, y = np.eye(2), x.copy()
        for _ in range(fp.period):
            analytic = jac(y) @ analytic
            y = _apply(system, y.reshape(2, 1), 1).ravel()
        J = np.asarray(fp.jacobians[i], dtype=float)
        lu = float(np.ravel(fp.unstable_eigenvalues[i])[0])
        ls = float(np.ravel(fp.stable_eigenvalues[i])[0])
        vu = np.ravel(fp.unstable_eigenvectors[i]).astype(float)
        vs = np.ravel(fp.stable_eigenvectors[i]).astype(float)
        orbit.append({
            "x": x.tolist(),
            "fp_residual": float(np.linalg.norm(image - x)),
            "lambda_u": lu,
            "lambda_s": ls,
            "det_minus_1": abs(lu * ls) - 1.0,
            "evec_res_u": float(np.linalg.norm(J @ vu - lu * vu) / np.linalg.norm(vu)),
            "evec_res_s": float(np.linalg.norm(J @ vs - ls * vs) / np.linalg.norm(vs)),
            "jac_vs_analytic": float(np.abs(J - analytic).max()),
        })
    return {"label": fp.label, "period": fp.period, "k_value": fp.k_value, "orbit": orbit}


def invariance_checks(session, rec) -> None:
    """(b) one-step invariance of every manifold onto its image branch."""
    wb, system = session.workbench, session.workbench.dynamical_system
    for key, manifold in wb.manifolds.items():
        fp, stability = key[0], key[1]
        forward = stability == "unstable"
        image_key = fp.advance_key(key, 1 if forward else -1)
        tag = f"inv_{_key_tag(key)}"
        if image_key not in wb.manifolds:
            rec.metric(f"{tag}_missing_image", _key_tag(image_key))
            continue
        probe = invariance_error(
            manifold, system.map if forward else system.map_inv, stability,
            target=wb.manifolds[image_key],
        )
        rec.metrics.update(error_summary(probe["err"], prefix=f"{tag}_"))
        rec.metric(f"{tag}_eligible", probe["eligible"])
        rec.metric(f"{tag}_beyond", int(probe["beyond"].sum()))
        rec.metric(f"{tag}_on_curve", [probe["on_curve"], probe["vertices"]])
        keep = np.linspace(0, len(probe["err"]) - 1, min(300, len(probe["err"]))).astype(int)
        rec.series[tag] = (
            np.column_stack([probe["cdist"][keep], probe["err"][keep]]).tolist()
            if len(keep) else []
        )


def iterate_checks(session, rec) -> None:
    """(c) every registered iterate link against the real map."""
    wb = session.workbench
    registry, system = wb.intersection_registry, wb.dynamical_system
    table = registry.iterate_table
    by_n: dict[int, list[float]] = {}
    rows, skipped = [], 0
    for source, n, target in table.items():
        a, b = registry[source], registry[target]
        if a.is_synthetic or b.is_synthetic:
            skipped += 1
            continue
        err = float(np.linalg.norm(
            _apply(system, a.get_point().reshape(2, 1), n).ravel() - b.get_point()
        ))
        by_n.setdefault(n, []).append(err)
        rows.append([n, table.backward_depth(source), a.unstable_cdist, err])
    rec.metric("iter_links", len(rows))
    rec.metric("iter_skipped", skipped)
    for n, errs in sorted(by_n.items()):
        rec.metrics.update(error_summary(errs, prefix=f"iter_n{n}_"))
    step = max(1, len(rows) // 2000)
    rec.series["iter_links"] = rows[::step]


def lobe_checks(session, rec) -> None:
    """(d) signed lobe area of each bridge against its image lobe's."""
    wb = session.workbench
    table = wb.intersection_registry.iterate_table
    curves = curve_arrays(session)
    rows, ratios, flips, unclosable = [], [], 0, 0
    for bridge in wb.bridges:
        if bridge.id is None:
            continue
        a, b = bridge.id
        ia, ib = table[a, 1], table[b, 1]
        if ia is None or ib is None:
            continue
        lobe, image = lobe_polygon(session, curves, a, b), lobe_polygon(session, curves, ia, ib)
        if lobe is None or image is None:
            unclosable += 1
            continue
        area, image_area = signed_polygon_area(lobe), signed_polygon_area(image)
        rows.append([a, b, area, image_area])
        if area != 0:
            ratios.append(abs(image_area / area - 1))
            flips += int(np.sign(area) != np.sign(image_area))
    rec.metric("lobe_pairs", len(rows))
    rec.metric("lobe_sign_flips", flips)
    rec.metric("lobe_unclosable", unclosable)
    rec.metrics.update(error_summary(ratios, prefix="lobe_ratio_"))
    rec.series["lobes"] = rows


def arrangement_checks(session, rec) -> None:
    """(d') what the arrangement's area-verified ``image_of`` reports."""
    arrangement = session.arrangement()
    regions = list(arrangement.regions)
    rows, ratios, missing = [], [], 0
    for region in regions[:500]:
        image = arrangement.image_of(region)
        if image is None:
            missing += 1
            rows.append([region.area, None])
            continue
        rows.append([region.area, image.area])
        ratios.append(abs(image.area / region.area - 1))
    rec.metric("arr_regions", len(regions))
    rec.metric("arr_image_found", len(ratios))
    rec.metric("arr_image_missing", missing)
    rec.metrics.update(error_summary(ratios, prefix="arr_ratio_"))
    rec.series["arr_regions"] = rows


def _monotonicity(cdist: np.ndarray) -> list:
    steps = np.diff(np.asarray(cdist, dtype=float).ravel())
    if not len(steps):
        return [0, 0, 0, 0.0]
    return [int((steps < 0).sum()), int((steps == 0).sum()), int(len(steps)),
            float(-steps.min()) if steps.min() < 0 else 0.0]


def cdist_checks(session, rec) -> None:
    """(e) non-increasing and tied cdist steps along manifolds and bridges."""
    wb = session.workbench
    per = {_key_tag(key): _monotonicity(m.get_cdist_array()) for key, m in wb.manifolds.items()}
    bridges = [_monotonicity(b.get_cdist_array()) for b in wb.bridges]

    def total(rows):
        rows = rows or [[0, 0, 0, 0.0]]
        return {"decreasing": sum(r[0] for r in rows), "ties": sum(r[1] for r in rows),
                "steps": sum(r[2] for r in rows), "worst_drop": max(r[3] for r in rows)}

    rec.metric("cdist_per_manifold", {k: v[:3] for k, v in per.items()})
    rec.metric("cdist_manifolds", total(list(per.values())))
    rec.metric("cdist_bridges", total(bridges))


def run_case(name: str, rec) -> None:
    builder, (hk, hb) = CASES[name]
    with rec.stage("build"):
        case = builder()
    session = case.session
    wb = session.workbench
    rec.metric("case", name)
    rec.metric("points_unstable", total_points(session, "unstable"))
    rec.metric("points_stable", total_points(session, "stable"))
    rec.metric("n_crossings", sum(1 for _, ix in wb.intersection_registry if not ix.is_synthetic))
    rec.metric("n_bridges", len(wb.bridges))
    rec.checkpoint()
    with rec.stage("fixed_point"):
        rec.metric("fp", [fixed_point_checks(wb.dynamical_system, fp, hk, hb)
                          for fp in case.fixed_points])
    with rec.stage("invariance"):
        invariance_checks(session, rec)
    rec.checkpoint()
    with rec.stage("iterates"):
        iterate_checks(session, rec)
    with rec.stage("cdist"):
        cdist_checks(session, rec)
    rec.checkpoint()
    with rec.stage("lobes"):
        lobe_checks(session, rec)
    rec.checkpoint()
    with rec.stage("arrangement"):
        arrangement_checks(session, rec)


def _grow_k28(gpu: bool, steps: int):
    """The k28 recipe of ``cases.build_k28`` through the bridges, on CPU or GPU."""
    session = henon_session(2.8, area_cutoff=1e-7, gpu=gpu)
    fp = session.construct_fixed_point(henon_fixed_points(2.8)["saddle"])
    session.orient_eigenvectors(fp, PERIOD1_ORIENTATION)
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=steps)
    session.grow_until_turnaround(fp, "stable")
    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    return session


def run_gpu(config: dict, rec) -> None:
    """(g) CPU and GPU growth of one recipe, compared point by point."""
    rec.metric("case", "gpu_k28")
    rec.metric("steps", config["steps"])
    try:
        import cupy  # noqa: F401
        rec.metric("gpu_available", True)
    except ImportError:
        rec.metric("gpu_available", False)
        return
    with rec.stage("cpu_grow"):
        cpu = _grow_k28(False, config["steps"])
    with rec.stage("gpu_grow"):
        gpu = _grow_k28(True, config["steps"])
    rec.checkpoint()
    with rec.stage("compare"):
        cpu_m = {_key_tag(k): m.get_point_array() for k, m in cpu.workbench.manifolds.items()}
        gpu_m = {_key_tag(k): m.get_point_array() for k, m in gpu.workbench.manifolds.items()}
        rec.metric("cpu_points", {k: len(v) for k, v in cpu_m.items()})
        rec.metric("gpu_points", {k: len(v) for k, v in gpu_m.items()})
        same, diff, hausdorff = {}, {}, {}
        for tag in cpu_m.keys() & gpu_m.keys():
            p, q = np.asarray(cpu_m[tag]), np.asarray(gpu_m[tag])
            same[tag] = len(p) == len(q)
            if same[tag]:
                diff[tag] = float(np.abs(p - q).max())
            hausdorff[tag] = max(float(distance_to_polyline(p, q)[0].max()),
                                 float(distance_to_polyline(q, p)[0].max()))
        rec.metric("same_lengths", same)
        rec.metric("max_abs_diff", diff)
        rec.metric("hausdorff", hausdorff)

        def coords(session):
            return np.array([ix.coords for _, ix in session.workbench.intersection_registry
                             if not ix.is_synthetic], dtype=float).reshape(-1, 2)

        from scipy.spatial import cKDTree

        pc, pg = coords(cpu), coords(gpu)
        rec.metric("cpu_crossings", len(pc))
        rec.metric("gpu_crossings", len(pg))
        if len(pc) and len(pg):
            dev, _ = cKDTree(pg).query(pc)
            rec.metric("crossing_max_dev", float(dev.max()))
            rec.metric("crossing_unmatched", int((dev > 1e-6).sum()))
        rec.metric("cpu_bridges", len(cpu.workbench.bridges))
        rec.metric("gpu_bridges", len(gpu.workbench.bridges))


def run(config: dict, rec) -> None:
    if config["case"] == "gpu_k28":
        run_gpu(config, rec)
    else:
        run_case(config["case"], rec)


if __name__ == "__main__":
    harness.main(globals())
