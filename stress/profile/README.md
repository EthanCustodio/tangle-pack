# Numerics profile (campaign N0..Nk: linked lists to arrays)

`profile_numerics.py` measures where memory and time go in the numerics pipeline;
`profile_figures.py` draws its json. The baseline was taken on `a8bb8a0`; the LAST
slice of the campaign reruns the same command on its branch and draws both runs on one
axes.

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

The last slice (N10) reruns the first command above on its branch and draws it beside the
baseline:

```bash
env/bin/python stress/profile/profile_numerics.py --plot-only \
    --out figures/numerics_profile/<branch-sha> --compare figures/numerics_profile/a8bb8a0
```

## Identity probe

`identity_probe.py` checks D3 for every slice: the points, the registry and the topology of
every `cases.LAW_CASES` build and of every step of every `cases.SERIES` (58 builds), base
against branch. Each case runs in a fresh interpreter importing the checkout named by
`TANGLEPACK_REPO`; nothing is written outside `--out`.

```bash
TANGLEPACK_REPO=../tangle-pack-<slug>-base env/bin/python stress/profile/identity_probe.py dump --out <scratch>/base
TANGLEPACK_REPO=../tangle-pack-<slug>      env/bin/python stress/profile/identity_probe.py dump --out <scratch>/branch
env/bin/python stress/profile/identity_probe.py compare <scratch>/base <scratch>/branch   # exit 1 on any difference
```

One dump takes about 60 s (126 MB; k32_series 26 s, nested_outer_series 17 s, k28_series 8 s,
every other case about 1 s) and a compare about 1 s, measured on `a8bb8a0`. `compare` prints
one line per differing array (rows differing, length change, max |dxy|, max |dcdist|), a
unified diff per differing text, and last `identical: 58 builds, 5776 arrays, 638 texts` (at
`a8bb8a0`) or `differences in K items`.

Every blast of the probe runs with `diagnose=True`, and each dump carries one more text,
`same_stability_counts` (N2, 2026-10-24): one line per blast since the previous dump -- the
`SameStabilityCount` of the blast's resolution report (`in_region`, `outside`, the `repr` of
`suggested_area_cutoff`, the curve pairs with a partial bridge's address dropped), caught by
wrapping `TangleWorkbench.report_resolution` because `blast_zone` returns no report. So a slice
that touches the same-stability count is checked blast by blast, base against branch, with no
change to `compare`; `diagnose=True` leaves every other item as it was.

A slice that changes an accessor the probe reads edits only its readers (`curve_arrays`,
`curve_extent`, `bracket_cdists`, `image_curves`, `key_text`), in its own diff.
