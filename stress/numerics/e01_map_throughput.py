"""E1: raw Hénon map throughput, NumPy vs CuPy, and the real ``map_batch`` path.

One config per backend sweeps the batch size N (series indexed by ``n``):

* ``numpy``          -- ``henon_map`` on a host ``(2, N)`` array;
* ``cupy_e2e``       -- ``cp.asnumpy(f(cp.asarray(arr)))``, exactly what
  ``gpu.py``'s wrapper does (host -> device -> host included);
* ``cupy_compute``   -- ``f`` on data already on the device (synchronised);
* ``cupy_transfer``  -- the round trip alone, ``cp.asnumpy(cp.asarray(arr))``;
* ``map_batch_cpu`` / ``map_batch_gpu`` -- ``DynamicalSystem.map_batch`` on an
  ``(N, 2)`` array through a session, GPU disabled / enabled
  (``min_batch_points=1``): the code path growth actually takes, transposes
  included.

Two ``crossover`` configs bisect, in log N, the batch size where the GPU path
starts beating the CPU one (``raw``: numpy vs cupy_e2e; ``map_batch``: the two
map_batch backends).

Timing: one discarded warm-up call per N, then samples of ``inner`` calls each
(``inner`` sized so one sample lasts >= ~5 ms), at least 7 samples and >= 50 ms
in total; median and IQR of the per-call time.
"""

from __future__ import annotations

import gc
import math
import time

import numpy as np

import harness
from maps import henon_session
from tanglepack.examples.henon import henon_map

EXPERIMENT = "e01_map_throughput"
TIMEOUT = 300
WORKERS = 1  # timings: never share the machine with another worker

K = 2.8
BACKENDS = (
    "numpy",
    "cupy_e2e",
    "cupy_compute",
    "cupy_transfer",
    "map_batch_cpu",
    "map_batch_gpu",
)
SAMPLE_S = 0.005
TOTAL_S = 0.05
MIN_SAMPLES = 7


def configs(quick: bool) -> list[dict]:
    sizes = [100, 10_000, 1_000_000] if quick else [
        int(round(n)) for n in np.logspace(2, 7, 11)
    ]
    out = [{"kind": "sweep", "backend": b, "sizes": sizes} for b in BACKENDS]
    hi = 1_000_000 if quick else 10_000_000
    out += [{"kind": "crossover", "pair": p, "lo": 100, "hi": hi,
             "iterations": 4 if quick else 10} for p in ("raw", "map_batch")]
    return out


def time_per_call(call, sync=lambda: None) -> dict:
    """Median/IQR of one call's wall time (adaptive sample and repeat counts)."""
    call()  # warm-up: kernel compile, allocator, first-touch pages
    sync()
    gc.collect()
    t0 = time.perf_counter()
    call()
    sync()
    estimate = max(time.perf_counter() - t0, 1e-7)
    inner = max(1, math.ceil(SAMPLE_S / estimate))
    samples = max(MIN_SAMPLES, math.ceil(TOTAL_S / (inner * estimate)))
    per_call = []
    for _ in range(samples):
        t0 = time.perf_counter()
        for _ in range(inner):
            call()
        sync()
        per_call.append((time.perf_counter() - t0) / inner)
    q1, med, q3 = np.percentile(per_call, [25, 50, 75])
    return {"median_s": med, "q1_s": q1, "q3_s": q3, "samples": samples, "inner": inner}


def make_call(backend: str, n: int):
    """The zero-argument call timed for one backend at batch size n, and its sync."""
    rng = np.random.default_rng(0)
    points = rng.uniform(-3, 3, size=(2, n))  # (2, N): coordinate on axis 0
    f = henon_map(K)
    if backend == "numpy":
        return (lambda: f(points)), (lambda: None)
    if backend.startswith("map_batch"):
        session = henon_session(K, gpu=backend.endswith("gpu"), min_batch_points=1)
        system = session.workbench.dynamical_system
        coords = np.ascontiguousarray(points.T)  # (N, 2), as the growth code holds it
        return (lambda: system.map_batch(coords)), harness.gpu_sync
    import cupy as cp

    if backend == "cupy_e2e":
        return (lambda: cp.asnumpy(f(cp.asarray(points)))), (lambda: None)
    if backend == "cupy_transfer":
        return (lambda: cp.asnumpy(cp.asarray(points))), (lambda: None)
    if backend == "cupy_compute":
        on_device = cp.asarray(points)
        return (lambda: f(on_device)), harness.gpu_sync
    raise ValueError(backend)


def timed(backend: str, n: int) -> dict:
    call, sync = make_call(backend, n)
    return time_per_call(call, sync)


def run(config: dict, rec) -> None:
    if config["kind"] == "sweep":
        sweep(config, rec)
    else:
        crossover(config, rec)


def sweep(config: dict, rec) -> None:
    backend = config["backend"]
    rec.metric("backend", backend)
    rec.metric("k", K)
    for n in config["sizes"]:
        t = timed(backend, n)
        rec.append("n", n)
        for key, value in t.items():
            rec.append(key, value)
        rec.append("points_per_s", n / t["median_s"])
        rec.checkpoint()


def crossover(config: dict, rec) -> None:
    """Bisect log N for the size where the GPU path's median beats the CPU one."""
    cpu, gpu = {"raw": ("numpy", "cupy_e2e"),
                "map_batch": ("map_batch_cpu", "map_batch_gpu")}[config["pair"]]
    rec.metric("cpu_backend", cpu)
    rec.metric("gpu_backend", gpu)

    def gpu_ratio(n: int) -> float:
        ratio = timed(gpu, n)["median_s"] / timed(cpu, n)["median_s"]
        rec.append("n", n)
        rec.append("gpu_over_cpu", ratio)
        rec.checkpoint()
        return ratio

    lo, hi = config["lo"], config["hi"]
    if gpu_ratio(lo) < 1:
        rec.metric("crossover_n", lo)
        rec.metric("crossover_note", "GPU already faster at the smallest size")
        return
    if gpu_ratio(hi) >= 1:
        rec.metric("crossover_n", None)
        rec.metric("crossover_note", "GPU never faster up to the largest size")
        return
    for _ in range(config["iterations"]):
        mid = int(round(math.sqrt(lo * hi)))
        if gpu_ratio(mid) < 1:
            hi = mid
        else:
            lo = mid
    rec.metric("crossover_n", int(round(math.sqrt(lo * hi))))
    rec.metric("crossover_bracket", [lo, hi])


if __name__ == "__main__":
    harness.main(globals())
