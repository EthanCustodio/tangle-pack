"""Run-to-run determinism and timing spread, CPU and GPU.

One config per repeat, each in a fresh process: the k=2.8 saddle at
``area_cutoff = 1e-7``, 10 unstable steps, stable grown to turnaround,
crossings, trim, bridges and iterate table (the ``build_k28`` recipe through
``bridges``), on the CPU or through ``enable_gpu``. Records the stage times and
sha256 hashes of the final manifold point arrays, the bridge point arrays and
the sorted crossing coordinates, both of the exact bytes and rounded to
``1e-12``. The figure-maker reads the timing CV and the bit identity across
repeats (and CPU vs GPU).
"""

from __future__ import annotations

import hashlib

import numpy as np

import harness
from maps import saddle_session

EXPERIMENT = "e10_determinism"
TIMEOUT = 300
WORKERS = 1


def configs(quick: bool) -> list[dict]:
    repeats = 1 if quick else 5
    return [{"device": device, "repeat": r} for device in ("cpu", "gpu") for r in range(repeats)]


def hashes(arrays: list[np.ndarray]) -> dict[str, str]:
    exact, rounded = hashlib.sha256(), hashlib.sha256()
    for array in arrays:
        array = np.ascontiguousarray(array, dtype=np.float64)
        exact.update(array.tobytes())
        rounded.update((np.round(array, 12) + 0.0).tobytes())
    return {"exact": exact.hexdigest(), "rounded": rounded.hexdigest()}


def run(config: dict, rec: harness.Recorder) -> None:
    with rec.stage("setup"):
        session, fp = saddle_session(2.8, 1, area_cutoff=1e-7, gpu=config["device"] == "gpu")
    with rec.stage("grow_unstable"):
        session.grow_n_times(fp, "unstable", num_iterations=10)
    with rec.stage("grow_stable"):
        session.grow_until_turnaround(fp, "stable")
    with rec.stage("intersections"):
        session.compute_intersections([fp])
    with rec.stage("bridges"):
        session.trim_stable_manifolds(fp)
        session.create_bridges(fp)
        session.infer_iterate_table()

    wb = session.workbench
    for stability in ("unstable", "stable"):
        arrays = [m.get_point_array() for k, m in sorted(
            wb.manifolds.items(), key=lambda item: item[0][2:]) if k[1] == stability]
        rec.metric(f"{stability}_points", sum(len(a) for a in arrays))
        rec.metric(f"{stability}_hash", hashes(arrays))
    crossings = np.array(sorted(ix.coords for _id, ix in wb.intersection_registry))
    cdists = np.array(sorted((ix.unstable_cdist, ix.stable_cdist) for _id, ix in wb.intersection_registry))
    bridges = sorted(wb.bridges, key=lambda b: b.get_cdist_array()[0].item())
    rec.metric("crossings", len(crossings))
    rec.metric("bridges", len(bridges))
    rec.metric("crossing_hash", hashes([crossings]))
    rec.metric("cdist_hash", hashes([cdists]))
    rec.metric("bridge_hash", hashes([b.get_point_array() for b in bridges]))
    rec.metric("crossing_coords", crossings)


if __name__ == "__main__":
    harness.main(globals())
