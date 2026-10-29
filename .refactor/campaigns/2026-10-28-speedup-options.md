# Campaign: GPU and multicore CPU speedup options

- **Goal:** offer both GPU and multicore CPU speedup options (author, 2026-10-28), for one long series
  (threads, device) and for many independent jobs (a process pool), with expensive ODE maps in view.
- **Date:** 2026-10-28. **Base:** main `5ed9940`. **Cap:** 4 runs (this is 3; each one delivers a usable
  option, so the campaign may stop early).
- **Author decisions at planning:** `workers` defaults to 1 (opt-in); a sweep runner is in scope, in
  `scripts/` (D6 "no scheduler in the library" holds); maps will often be ODE integrations, so the map is a
  target; a lean-blast run (map only new bridges, drop older point data) was planned, then
  dropped by the author the same day.
- **Scout (Hénon, figures/numerics_profile/6860f99*):** 1e-7 k=2.5 blast series 32.5 s, of which the
  same-stability sweep (`diagnostics._hits` / `_candidate_pairs`) is 22.4 s, separation
  (`Blast._min_interior_distance`, a cKDTree per child) 4.8 s, and the map about 0%. At 1e-9 the sweep is
  30 of 61 s and separation 8.4 s. Hardware: 32 cores, 125 GB, RTX 5060 Ti 16 GB, CuPy 14.1.1.
- **Measuring rule:** time only on an idle machine or on `taskset`-pinned free cores. The author runs about
  16 sweeps at once, and a contended probe read 116 s for a 32 s series.
- **Every run:** `env/bin/python -m pytest -q -rxX` (only KNOWN_ISSUES xfail), `coverage_guard` against the
  base, `stress/profile/identity_probe.py` unchanged on CPU, `profile_numerics.py --cases scale` before and
  after.

## Outcome (2026-10-29, complete in 3 of 4 runs)

Both kinds of speedup now exist, all opt-in, nothing computed changed on the CPU path (identity probe identical
after runs 1 and 2). CPU: `workers` / `map_executor` split the map (solve_ivp Duffing 4.4x at 16 processes)
and the sweep's pair blocks (Hénon 1e-7 blast 1.55x). GPU: `enable_gpu(route_map=False)` runs the
same-stability sweep on the device with a byte-equal state (1e-7 scale blast 50.1 -> 17.4 s), and a map written
as one fused `--fmad=false` kernel is bit-identical and fast (Duffing 2000-step growth 1.93 -> 0.56 s). Many
jobs: `scripts/sweep_pool.py`, adopted by the period-3 zone sweep and the k25 entropy sweep. Open: chunked
zone sweeps differ slightly from one continuous run (warm-up path; as `--chunk` did); the 1e-5 seed crash in
`Tangle._query_curve` (Q7); a GPU flag for pool jobs; the hard time factor 2 is unmeasured.

## Run 1: cpu-workers (/feature)
- [x] done (5dec003, workers/map_executor opt-in: solve_ivp Duffing 4.4x at 16, Hénon 1e-7 blast 1.55x / 1e-9 1.26x capped by the serial whole-field sweep; separation unchanged; identity probe identical)

Goal: `workers=N` (default 1) on `TangleWorkbench` / `TangleSession` runs three things on N cores, with
nothing computed changing:
1. **The map.** `DynamicalSystem.map_batch` splits a large `(2, N)` batch into chunks, one per worker. The
   run's plan chooses threads for maps that release the GIL (vectorised integrators) and processes (fork,
   or picklable maps) for per-point `solve_ivp` maps. There is a `min_batch_points` floor.
2. **The same-stability sweep.** `_candidate_pairs` strips / pair blocks and their pair test in `_hits` run
   on a thread pool, merged in strip order so the listing is unchanged.
3. **Separation.** One cKDTree per generation plus an incremental sibling tree, instead of one per child;
   `query(..., workers=N)`.

Also new: `examples/pendulum.py`, a periodically forced pendulum's stroboscopic map. Leapfrog makes it
exactly symplectic and reversible, so it preserves area and has an exact inverse. It comes in two forms:
vectorised array-API, and per-point `solve_ivp`. It is the expensive-map benchmark for runs 1 and 2.

Files: `numerics/DynamicalSystem.py`, `numerics/TangleWorkbench.py`, `loom/TangleSession.py`,
`numerics/diagnostics.py`, `loom/Blast.py`, `examples/pendulum.py`, `numerics/persist.py` (`workers` not
saved), tests (worker-count invariance, the pendulum's area and inverse), CLAUDE.md rows.

Success:
- the identity probe and SameStabilityCount are identical at workers 1 and 8;
- the suite passes and `coverage_guard` holds;
- blast time of the 1e-7 and 1e-9 scale series at workers 1/4/8/16, with diagnose and separation split out;
- pendulum growth time at workers 1/8/16 for both map forms.

## Run 2: gpu-device (/perf)
- [x] done (3a35d74, device sweep + fused Duffing kernel + enable_gpu(route_map=False): 1e-7 scale blast 50.1 -> 17.4 s, diagnose 38.7 -> 5.9 s, state byte-equal; Duffing 2000-step growth 1.93 s CPU / 19.6 s old GPU -> 0.56 s)

Goal: `enable_gpu` speeds up both expensive maps and the same-stability sweep.
- **The map:** the vectorised pendulum on the device against CPU workers. Tune `min_batch_points` against
  refinement's layer batch sizes, and optionally keep a batch on the device across the integrator steps.
- **The sweep:** `_hits` / `_candidate_pairs` through `numerics/namespace.py`'s `array_namespace`, uploaded
  once per call when the GPU is enabled.

Candidates: (a) only the pair-test blocks on the device; (b) the whole sweep on the device.

Files: `numerics/gpu.py` (Measured result and What runs where rewritten), `numerics/diagnostics.py`,
`numerics/namespace.py`.

Success:
- SameStabilityCount is identical to the CPU on the scale series, or every difference is measured and
  explained (FMA at the margin);
- pendulum growth and blast-series times for GPU, CPU workers 1 and CPU workers 8;
- the map's break-even batch size is recorded.

## Run 3: sweep-runner (/feature)
- [x] done (dcbe49e, scripts/sweep_pool.py: processes x workers <= cores, per-job RSS/time caps killing the tree, resume per job / per k, one merged CSV; zone sweep --processes replaces --chunk/--merge; k25 `pool` subcommand; 29 tests)

Goal: a scripts-level process pool, `scripts/sweep_pool.py`, for independent jobs (k × cutoff × case) over
N processes, with:
- per-job RSS and time caps;
- resume from the CSV and npz on disk (`TangleWorkbench.save` / `load`);
- one merged CSV;
- `workers` per job, with `processes × workers <= cores`.

`henon_period3_zone_sweep.py` (replacing `--chunk` / `--merge`) and the k25 entropy sweep adopt it. It is
not in `src/`.

Files: `scripts/sweep_pool.py`, the two sweep scripts, CLAUDE.md Commands.

Success: one command reproduces a chunked zone-sweep range and the multi-cutoff entropy launch with
identical CSV rows; wall time is recorded against serial.
