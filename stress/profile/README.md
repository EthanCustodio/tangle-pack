# Numerics profile (campaign N0..Nk: linked lists to arrays)

`profile_numerics.py` measures where memory and time go in the numerics pipeline;
`profile_figures.py` draws its json. The baseline was taken on `a8bb8a0`; the last
slice of the campaign (N10b) reran the same commands on `dc1bd0c` and drew both runs on
one axes.

## Commands (from the repository root)

```bash
# the baseline (and, unchanged, the after-profile of the last slice)
env/bin/python stress/profile/profile_numerics.py
# -> figures/numerics_profile/<git-sha>/  (json per case and mode, png figures)

# the largest scaling point (k = 2.5 period 1 at area_cutoff 1e-9)
env/bin/python stress/profile/profile_numerics.py --cases scale --modes time \
    --scale-cutoff 1e-9 --out figures/numerics_profile/<git-sha>_scale1e-9

# the comparison figures: the branch run beside the baseline
env/bin/python stress/profile/profile_numerics.py --plot-only \
    --out figures/numerics_profile/<branch-sha> --compare figures/numerics_profile/a8bb8a0
```

Options: `--cases` (k10 p3 nested inversion heteroclinic scale), `--modes` (time memory),
`--scale-cutoff` (1e-7), `--scale-stop-gb` (8: stop after the blast whose peak RSS crosses
it), `--scale-max-blasts` (40; the series also stops on a plateau, two blasts without a new
point), `--scale-trace-blasts` (8: blasts run under tracemalloc), `--cap-gb` (12: the
parent kills a worker past it).

## What it measures

Each (case, mode) runs in a fresh interpreter. The law cases are the frozen
`tests/cases.py` builders; `scale` is the k = 2.5 period-1 saddle of
`scripts/henon_k25_entropy.py` (10 unstable steps, stable to turnaround, zone at f(q0)),
blasted one single-iteration step at a time with `min_separation = 1e-5` and
`diagnose=True`, the topology (and its refined transition matrix) recomputed after every
blast.

- **time mode**: library functions are wrapped (lists `TOP_STAGES` / `INNER_STAGES`); each
  call path gets its exclusive self time, so `blast > iterate_bridge > iterate > refine >
  map` splits a blast. RSS is sampled every 20 ms; its peak is kept per top-level call. A
  structural census (`census()`) sizes what the session holds by component, at the end and
  after every blast of the scaling case; the rtree's C memory is measured once at the end
  by dropping the two indices and reading RSS (glibc `malloc_trim`).
- **memory mode**: the same build under `tracemalloc` (imports done before the baseline):
  the traced growth is the reference the census must sum to, and the json keeps the top
  allocation sites.

Reruns reproduce memory to about 1% (census, traced) and times within noise (sub-second
law cases jitter by 20-40%).

## Renamed functions

A wrapped function that no longer exists fails the run (`AttributeError`). When a slice
renames one, add the new name under the SAME label so the figures stay comparable.

## Baseline

The baseline of `a8bb8a0` is committed, unchanged from the analyst's run:

- `figures/numerics_profile/a8bb8a0/` -- the five law cases (time and memory) and the 1e-7
  scaling series (time; memory = 8 traced blasts);
- `figures/numerics_profile/a8bb8a0_scale1e-9/` -- the largest scaling point (1e-9);
- `figures/numerics_profile/a8bb8a0_scale_traced/` -- the 1e-7 series run to its plateau under
  tracemalloc (json only).

## After-profile (N10b)

The last slice (N10b, 2026-10-28) reran the first two commands above, unchanged, on `dc1bd0c`
(whose `src` and `tests` are `08b1d12`'s, the N10a merge), reran the traced
series, and drew both cutoffs beside the baseline:

- `figures/numerics_profile/dc1bd0c/` -- the five law cases and the 1e-7 scaling series
  (12 json, 6 png, 1.3 MB);
- `figures/numerics_profile/dc1bd0c_scale1e-9/` -- the largest scaling point (1 json, 5 png,
  508 KB; a folder with no law case has no `stage_peak_rss.png`);
- `figures/numerics_profile/dc1bd0c_scale_traced/` -- the 1e-7 series to its plateau under
  tracemalloc (1 json, 264 KB).

```bash
env/bin/python stress/profile/profile_numerics.py --cases scale --modes memory \
    --scale-trace-blasts 40 --out figures/numerics_profile/dc1bd0c_scale_traced
env/bin/python stress/profile/profile_numerics.py --plot-only \
    --out figures/numerics_profile/dc1bd0c --compare figures/numerics_profile/a8bb8a0
env/bin/python stress/profile/profile_numerics.py --plot-only \
    --out figures/numerics_profile/dc1bd0c_scale1e-9 --compare figures/numerics_profile/a8bb8a0_scale1e-9
```

The result, `a8bb8a0` -> `dc1bd0c`: census at the plateau 1914.0 -> 231.9 MiB and blast
458.4 -> 52.2 s at 1e-7, 4683.2 -> 502.3 MiB and 541.4 -> 67.7 s at 1e-9; the traced growth
1880.6 -> 232.3 MiB. Every blast's `pip_depth`, `points_created` and stage calls are those of
the baseline (32 blasts at 1e-7, 29 at 1e-9). The peak RSS fell least (3034.7 -> 2015.0 MiB
at 1e-7, 6632.6 -> 1930.8 at 1e-9): it is a transient inside a blast and bimodal from run to
run (docs/numerics_arrays_plan_2026-10-23.md, Outcome row 10b).

## Diagnose on new segments (D1)

D1 (2026-10-28) runs the same-stability count's sweep inside strips across the sweep axis
(`numerics/diagnostics.py`), and reran the first two commands above, unchanged, on `6860f99`,
drawn with `--compare` beside `dc1bd0c` (the same `src` but for save / load, which the profiler
never runs): `figures/numerics_profile/6860f99/` (12 json, 6 png, 1.3 MB) and
`figures/numerics_profile/6860f99_scale1e-9/` (1 json, 5 png, 516 KB). The result, `dc1bd0c` ->
`6860f99`: `diagnose` 42.0 -> 22.4 s and blast 52.2 -> 32.5 s at 1e-7, `diagnose` 52.6 -> 30.0 s
and blast 67.7 -> 45.2 s at 1e-9. The peak RSS over the baseline fell 2015.0 -> 1054.6 MiB at
1e-7; at 1e-9 it moved from blast 12 to blast 19 and reads 1930.8 -> 2008.1 MiB (one run each).
Nothing computed moved: every blast's `same_stability_count`, `pip_depth` and `points_created`
and the census (231.9 MiB at 1e-7, 502.3 MiB at 1e-9) are `dc1bd0c`'s (32 blasts at 1e-7, 29 at
1e-9).

## Identity probe

`identity_probe.py` checks D3 for every slice: the points, the registry and the topology of
every `cases.LAW_CASES` build and of every step of every `cases.SERIES` (62 builds at
`08b1d12`), base against branch. Each case runs in a fresh interpreter importing the checkout
named by `TANGLEPACK_REPO`; nothing is written outside `--out`.

```bash
TANGLEPACK_REPO=../tangle-pack-<slug>-base env/bin/python stress/profile/identity_probe.py dump --out <scratch>/base
TANGLEPACK_REPO=../tangle-pack-<slug>      env/bin/python stress/profile/identity_probe.py dump --out <scratch>/branch
env/bin/python stress/profile/identity_probe.py compare <scratch>/base <scratch>/branch   # exit 1 on any difference
```

One dump takes about 27 s (129 MB) at `08b1d12`, about 60 s (126 MB; k32_series 26 s,
nested_outer_series 17 s, k28_series 8 s, every other case about 1 s) on `a8bb8a0`, and a
compare about 1 s. `compare` prints one line per differing array (rows differing, length
change, max |dxy|, max |dcdist|), a unified diff per differing text, and last
`identical: 58 builds, 5776 arrays, 638 texts` (at `a8bb8a0`; `identical: 62 builds, 6024
arrays, 744 texts` at `08b1d12`) or `differences in K items`.

Every blast of the probe runs with `diagnose=True`, and each dump carries one more text,
`same_stability_counts` (N2, 2026-10-24): one line per blast since the previous dump -- the
`SameStabilityCount` of the blast's resolution report (`in_region`, `outside`, the `repr` of
`suggested_area_cutoff`, the curve pairs with a partial bridge's address dropped), caught by
wrapping `TangleWorkbench.report_resolution` because `blast_zone` returns no report. So a slice
that touches the same-stability count is checked blast by blast, base against branch, with no
change to `compare`; `diagnose=True` leaves every other item as it was.

A slice that changes an accessor the probe reads edits only its readers (`curve_arrays`,
`curve_extent`, `bracket_cdists`, `image_curves`, `key_text`), in its own diff.
