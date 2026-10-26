# Numerics on arrays — the campaign plan (decided 2026-10-23)

Written in run N0 (`numerics-arrays-skeleton`) against `a8bb8a0`, from the author's request (D1-D6
below), the analyst's measured profile and inventory, and the planner's split. The author approved
the plan on 2026-10-23 with every default: "Approved as written; all Q1-Q9 defaults accepted (Q8:
same-stability count joins the campaign as N2; Q6 order N1..N10)."

The campaign moves `tanglepack.numerics` from the doubly linked `Point` / `BasePoint` lists to one
columnar point store with named views over it, slice by slice, one `tangle_refactor` run per slice
(section C). This doc carries the measurements and the inventory so that no later planner re-derives
them. **Every file:line below is at `a8bb8a0`**; a slice re-reads the lines it touches.

The tools this run committed:

- `stress/profile/profile_numerics.py`, `profile_figures.py` — the profiler (D2); its README states
  the commands.
- `stress/profile/identity_probe.py` — the base-vs-branch identity probe (D3), run by every slice.
- `figures/numerics_profile/a8bb8a0*/` — the baseline json and figures (Q1c).

---

## A. Goal and decisions

The author's goal, in their words: "change from the double linked list to the arrays. It is time. I
still want the intuitive class based structure and it to be structured how the topology is. That was a
nice feature of the double linked list. But we want to be able to have things take up far far less
data. And also be able to vectorize and GPU accelerate and scale to cluster computing more easily."

"Structured how the topology is": the topology layer is classes that are views over one store (the
`Trellis` over the `IntersectionRegistry`, with ids, not object graphs). Numerics gets the same shape:
one columnar point store, the single store of points as the registry is of crossings, with `Manifold`
/ `Bridge` / `BranchPoint` as named views over it.

### A.1 The author's decisions (2026-10-23, before the run)

- **D1. Skeleton, then split.** N0 profiles, scopes and splits; each later slice is one
  `tangle_refactor` run. No multi-slice driver.
- **D2. Profile before and after.** The profiler is tracked; the LAST slice reruns it unchanged on its
  branch and draws the same figures beside the baseline.
- **D3. Byte-identical where it can be, the laws always.**
  - The gate: every invariant, every law, and `KNOWN_ISSUES` unchanged (C3 alone, no XPASS).
  - The target: byte-identical points, registry and topology on the law cases and series. A
    difference is allowed when it is measured and explained: case, number of points, max coordinate
    and cdist difference, cause (a tie, a summation order, a vectorised refinement order).
  - Topology outputs (classes, words, matrices, `pip_depth`) that move are listed case by case and
    judged as discovery or bug, never silently accepted.
- **D4. Index-based public API, human-readable.** `Point` / `BasePoint` go; curves expose arrays and
  integer point ids, with no `PointRef` shim. Named accessors stay: `bridge.head` / `bridge.tail`,
  `manifold.root` / `manifold.tail`, `len(curve)`, `curve.points` / `curve.cdists` (section D).
- **D5. Be very careful with `BranchPoint`.** It stays a class with its readable API
  (`num_branches`, `get_cdist(stability)`, which branches attach in which slot). Its storage may
  change, in its own slice, last or near last, never folded into another (section G).
- **D6. Scale target: ready, not wired.** Arrays throughout; hot paths written against an array
  namespace (numpy now, CuPy swappable); workbench state serialisable. No distributed scheduler; GPU
  stays opt-in (`enable_gpu`).

### A.2 The questions, as answered (2026-10-23: every default)

- **Q1. The profiler.**
  - (a) It lives in `stress/profile/`, as the analyst built it, and does not reuse
    `stress/numerics/harness.py` (its stage timer is inclusive and double counts nested stages).
  - (b) Cases: the five law cases `k10`, `p3`, `nested`, `inversion`, `heteroclinic` plus the k=2.5
    series. The series stops at its plateau (the request's "~8 GB" is not reachable on `HEAD`), with
    `--scale-stop-gb 8` kept as a guard. 1e-7 is the default; 1e-9 is the largest point.
  - (c) The baseline json and png are committed under `figures/numerics_profile/` (1.3 MB), so N10
    draws beside them without rerunning the base.
- **Q2. Storage.** One shared append-only `PointStore` per workbench, int64 ids; curves are id
  arrays. Row ids never move, so brackets stay valid. Phantom preimages become ROWS on no curve (40 B
  each), not two columns on every row; the same "row on no curve" kind covers merge-dropped
  duplicates and the initializer's fictitious chain. int32 is a later measurement (est. -25 %), not
  a default.
- **Q3. Ties.** Keep the geometric splice semantics: merge by `searchsorted(..., side="right")` plus
  "existing wins" (all 19 measured ties are this case); refine by positional insert, never
  `argsort(cdist)`; unit tests pin the two tie branches no case ever reached.
- **Q4. The segment index.** Keep the two rtrees, one per stability, keyed by ints. No vectorised
  broad phase: the rtrees are 0.02 % of memory and 2-5 % of time.
- **Q5. Serialisation.** `save(path)` writes one npz plus a JSON sidecar (the store, curves, bridges
  with `iterated` / `dropped_at`, the registry with brackets as ids, the iterate table, fixed points
  and the machine's settings); `load(path, map, map_inv, jacobian)` recomputes the session and the
  topology and rebuilds the Tangle index. No HDF5.
- **Q6. Order.** N1 slim index, N2 same-stability count, N3 store, N4 order arrays, N5
  segments/brackets, N6 readers, N7 delete `Point`, N8 `BranchPoint`, N9 namespace, N10 save/load and
  after-profile. The index goes first although it lands before the store: it is the largest measured
  component and needs only an int identity per point (a stopgap if `id()` will not do).
- **Q7. Test churn.** Tests that pin list mechanics are deleted or rewritten as array-invariant tests
  (cdist non-decreasing along `ids`; image/preimage a partial bijection; no id twice on one curve).
  Law TEST functions stay unchanged; the law HELPERS (`tests/helpers/laws.py`, `invariants.py`)
  change in the slice that changes the shape they read (brackets in N5, nodes in N7), and the
  reviewer diffs them line by line. Each slice lists its collected-id diff.
- **Q8. The same-stability count** joins the campaign as **N2**: candidate pairs only from pairs that
  hold a new segment, counts unchanged. It is 82-87 % of blast time at scale and not a linked-list
  cost; without it the after-profile shows blast time almost flat.
- **Q9. `BranchPoint.num_branches`.** Both names and meanings stay: `BranchPoint.num_branches` is
  always 2 (eigendirection slots), `FixedPoint.num_branches` is 1 or 2 (branches one chain visits).
  N8 documents the difference in its docstrings; no rename. N8 retires roles 4 and 5 only with the
  probe's evidence.

---

## B. The baseline (`a8bb8a0`)

Committed under `figures/numerics_profile/`:

| folder | what | figures to open first |
|---|---|---|
| `a8bb8a0/` | the five law cases (time + memory) and the 1e-7 k=2.5 series (time; memory = 8 traced blasts) | `scaling_points.png`, `stage_time.png`, `memory_components.png`, `stage_peak_rss.png` |
| `a8bb8a0_scale1e-9/` | the largest scaling point (1e-9, time mode, census per blast) | `scaling_bytes_per_point.png` (the segment index is the amber band), `scaling_time_rss.png`, `memory_components.png`, `stage_time.png` |
| `a8bb8a0_scale_traced/` | the 1e-7 series to its plateau under tracemalloc (json only) | `scale_memory.json` |

GB are GiB (2^30 bytes), MB MiB; all memory is over the interpreter baseline after imports. The law
builds take 0.1-0.8 s and 2-21 MB: their times are noise (20-40 % run to run) and serve as a
regression check only. The signal is the scaling case.

**The answer in one paragraph.** The 45 GB peak that started the campaign no longer happens on `HEAD`
(it came from the old blast re-iterating dropped bridges, fixed by blast-steps-compose); the k=2.5
series plateaus at 2.9 GB (1e-7) and 6.5 GB (1e-9). At the 1e-9 plateau the segment index is 2.60 GB,
57 % of the census, almost all of it per-segment Python bookkeeping for 3.26 M query-only blast
segments; `Point` objects are 1.90 GB, 41 %, of which 0.65 GB are phantom preimages; the rtrees use
1 MB and the registry and topology less than 0.1 %. Arrays throughout would cut the held memory about
20x (est.). Arrays barely change time: at scale 82-87 % of blast time is the same-stability count
(N2), the linked list 10-15 %, the map 0.01 %.

The tables below are quoted from the analyst's findings (sections 3.2, 3.3, 3.5); "section N" inside
a quote refers to the findings, which this doc maps as 3 -> B, 4.1-4.2 -> E, 4.3 -> D, 5 -> F, 6 -> G.
The json named in a quote is under `figures/numerics_profile/`.

### B.1 Scaling, default cutoff (k = 2.5 p1, 1e-7): plateau at 2.18 M Point objects, 2.9 GB

`a8bb8a0/scale_time.json`; stopped "plateau: no new points since blast 30"; wall 562 s.

| blast | Point objects (image / phantom) | segments | RSS GB | step peak GB | B/pt RSS | B/pt census | Points (no phantoms) % | phantom % | segment index % | step s (blast / topology) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 8,848 (0 / 3,703) | 4,226 | 0.01 | 0.01 | 1581 | 817 | 29 | 20 | 48 | 0.0 |
| 4 | 32,154 (14,816 / 12,193) | 19,034 | 0.04 | 0.04 | 1281 | 889 | 28 | 16 | 53 | 0.2 (0.1 / 0.0) |
| 8 | 233,695 (142,041 / 86,509) | 146,195 | 0.25 | 0.42 | 1133 | 919 | 28 | 15 | 54 | 1.8 (1.8 / 0.1) |
| 10 | 603,430 (377,288 / 220,997) | 381,321 | 0.72 | 1.24 | 1272 | 956 | 27 | 15 | 57 | 6.1 (6.0 / 0.1) |
| 12 | 1,113,882 (704,747 / 403,990) | 708,610 | 1.15 | 2.44 | 1104 | 964 | 27 | 14 | 57 | 12.3 (12.1 / 0.1) |
| 14 | 1,543,122 (979,720 / 558,257) | 983,441 | 1.55 | 2.81 | 1076 | 933 | 28 | 15 | 56 | 17.4 (17.3 / 0.1) |
| 20 | 2,062,344 (1,316,749 / 740,450) | 1,320,305 | 2.37 | 2.85 | 1234 | 914 | 29 | 15 | 55 | 24.1 (24.0 / 0.1) |
| 30 | 2,177,786 (1,394,339 / 778,302) | 1,397,859 | 2.80 | 2.92 | 1383 | 912 | 29 | 15 | 55 | 25.1 (25.0 / 0.1) |

The request's "45 GB at blast 14, >10 KB per point" is the PRE-`a8bb8a0` sweep
(`figures/henon_k25_entropy/data/pre_blast_steps_compose/p1_a1e-07.json`, "killed by hand: blast 14
peaked past 45 GB", 1.98 M points at blast 13). At `HEAD` the same case holds 1.54 M Point objects at
blast 14 and peaks at 2.8 GB; the current `p1_a1e-07.json` (the 2026-10-21 rerun) agrees (peak 3.25 GB
for the whole 80-blast run). I did not find what made the old run's peak ~13x its resting memory
(believed: the pre-blast-steps-compose blast re-iterating dropped bridges; not reproduced).

### B.2 Scaling, largest point (k = 2.5 p1, 1e-9): plateau at 5.10 M Point objects, 6.5 GB

`a8bb8a0_scale1e-9/scale_time.json`; stopped "plateau: no new points since blast 27"; wall 780 s;
peak RSS 6.48 GB over baseline, 6.28 GB at rest.

| blast | Point objects (image / phantom) | segments | RSS GB | step peak GB | B/pt RSS | B/pt census | Points (no phantoms) % | phantom % | segment index % | step s (blast / topology) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 38,206 (0 / 15,952) | 18,256 | 0.04 | 0.05 | 1153 | 800 | 30 | 20 | 49 | 0.0 |
| 8 | 1,071,488 (654,000 / 395,234) | 672,184 | 1.07 | 1.41 | 1076 | 907 | 28 | 16 | 54 | 7.8 (7.5 / 0.3) |
| 12 | 3,626,633 (2,290,894 / 1,313,485) | 2,308,893 | 3.57 | 4.63 | 1057 | 924 | 28 | 15 | 55 | 25.3 (24.5 / 0.8) |
| 16 | 4,573,124 (2,906,841 / 1,644,029) | 2,924,770 | 4.71 | 5.66 | 1107 | 966 | 27 | 14 | 58 | 29.4 (28.6 / 0.9) |
| 20 | 4,969,813 (3,168,276 / 1,779,283) | 3,186,178 | 5.87 | 6.09 | 1269 | 958 | 27 | 14 | 57 | 33.3 (32.5 / 0.9) |
| 29 | 5,098,154 (3,258,777 / 1,817,123) | 3,276,670 | 6.28 | 6.28 | 1322 | 963 | 27 | 14 | 57 | 0.1 |

Components at the plateau (census 4.68 GB; RSS at rest 6.28 GB, peak 6.48 GB):

| component | GB | % of census | % of peak RSS |
|---|---|---|---|
| segment index, Python (`_edge_to_sid` 331 B/seg, owner sets 368 B/seg, `_Segment` 100 B/seg, `_seg_lookup` 51 B/seg) | 2.60 | 56.8 | 40.1 |
| Point objects on blast image curves (3.26 M) | 1.24 | 27.1 | 19.1 |
| phantom preimages (1.82 M) | 0.65 | 14.2 | 10.0 |
| walk memo (bridge point arrays and node lists) | 0.075 | 1.6 | 1.2 |
| Point objects on registered manifolds (22 k) + BranchPoints | 0.009 | 0.2 | 0.1 |
| rtree C memory (18,256 entries; measured 1.0 MB) | 0.001 | 0.0 | 0.0 |
| registry, topology caches, curve objects, rest | 0.004 | 0.1 | 0.1 |
| not in any Python object (allocator slack, freed-not-returned arenas, transients) | 1.80 (peak) | - | 27.8 |

Grouped as the request asks: **Point objects 41.5 % of the census / 29.3 % of peak RSS; segment
index 56.8 % / 40.1 %; bridges and image curves 28.7 % / 20.3 %** (the image-only points 27.1 % + walk
memo 1.6 % + curve objects < 0.01 %; the image Points are counted in both groups). Per Point object the
census is ~920-960 B and RSS ~1.05-1.32 KB at rest; per REAL point (phantoms excluded) the segment
index alone is ~800 B.

Why the segment index is so big: a blast indexes each image curve QUERY-ONLY
(`Tangle.add_manifold(iterated, index_segments=False)`, `BridgeIterator.py:359`), so its segments
never enter an rtree, yet each keeps a `_Segment`, a frozenset edge key of two `id()` ints in
`_edge_to_sid`, an owner set in `_seg_manifolds`, and a membership in the curve's `_manifold_segs`
set; the bridges cut from the curve claim the same segments again (`Tangle.create_bridges` ->
`add_manifold(bridge, index_segments=False, detect_crossings=False)`, `Tangle.py:966-968`). After the blast
step nothing queries them geometrically; they serve ownership / dedup (`release_manifold`, the
edge dedup of `add_manifold`, `compute_intersections`' "curve adds no segment" test).

### B.3 Time at scale

| stage (1e-7 series, total 562 s) | s | inner split (self time) |
|---|---|---|
| blast (30 calls) | 458.4 | diagnose 398.6 (same-stability count), walks 18.7, segment index 11.2, infer_iterates 9.6, refine 8.9, separation (KD-trees) 5.1, resolved_ids 1.8, map 0.07 |
| topology (31 calls) | 2.9 | walks 1.3, own 1.2, trellis build 0.08 |
| grow + intersections + zone | 0.14 | refine 0.04, index 0.07 |

1e-9 series: blast 541 s (diagnose 444 = 82 %, walks 31, index 27, refine 22, separation 8.8, map
0.06), topology 16 s. cProfile of blast 12 at 1e-7 (scratch `cprof_blast.py`): `_candidate_pairs`
8.2 s of 18.2 s (its own Python/numpy body), `_hits` 9.7 s cumulative; node-to-array conversion in
`_segment_rows` ~1 s of it. **The linked list costs ~10-15 % of blast time** (walks + index + refine's
node handling); the map ~0.01 %; the same-stability count ~80-87 %. GPU stays irrelevant here.

---

## C. Slices

Ordered by measured effect, provided each can land green on its own (Q6). The figures are from the
1e-9 plateau: census 4.68 GB, blast 541 s. **Est.** marks an estimate from the measured counts, not a
measurement.

| # | slug | what | expected effect | byte identity (D3) |
|---|---|---|---|---|
| N1 | `slim-segment-index` | Query-only curves (blast image curves, the bridges cut from them) hold no per-segment record. Ownership and dedup move to the curve level, keyed by int point identities (a stopgap per-`Point` int if `id()` will not do). The two rtrees are unchanged (WP10 rule). | Removes up to **2.58 GB of 4.68** (measured: 2.60 GB held by the index, of which the 18 k indexed segments are about 15 MB). | Expected exact. Risk: which crossings dedup, and so the registry ids. |
| N2 | `same-stability-count-new-only` | `count_same_stability_crossings` builds candidate pairs only for pairs that contain a new segment, instead of sweeping every unstable curve and filtering afterwards. The counts are unchanged. | Up to **-444 s of 541 s** blast time at 1e-9 (est.: the cost is measured, the new cost is not). | Exact: the counts per blast are equal. |
| N3 | `point-store` | One append-only `PointStore` per workbench (x, y, cdist, image, preimage). `Point` becomes a slotted object over its row and keeps its `forward`/`backward` links as an adapter. Phantom preimages become store rows on no curve (Q2). | Est. **-1.3 to -1.5 GB**: `Point` drops from 400 B to about 120 B, and a phantom from 376 B to 40 B. Walks may get slower until N4 (property reads); the profiler shows it. | Expected exact. Low risk; the trap is the tie rule. |
| N4 | `curves-on-order-arrays` | Each curve holds an int `ids` array in curve order. Grow, merge and refine work on positions (`searchsorted`, `side="right"`, "existing wins", per Q3; refine by positional insert, never a sort). `_image_chains` uses a per-curve position index. May split at its own plan. | Time: refine plus walks are 53 s of 541 at 1e-9. Removes the forward/backward links. | Target exact. The 19 measured ties (k28, heteroclinic) must reproduce. |
| N5 | `segments-and-brackets-on-ids` | Segments become `(ids[:-1], ids[1:])`. `Intersection.unstable_segment` / `stable_segment` become store id pairs, and every bracket reader follows. Split at its plan (2026-10-25) into N5a `brackets-on-rows` and N5b `segments-and-edge-keys-on-rows`. | Small: the bookkeeping left after N1. | Exact. The same floats are read. |
| N6 | `readers-on-views` | Every reader outside the engine moves to the D4 accessors. The walk memo goes. | -0.075 GB (walk memo), plus the node-to-array conversion in `_segment_rows`. | Exact. |
| N7 | `delete-point` | `Point` and `BasePoint` are deleted and leave `__all__`; every test, helper, stress and script touchpoint moves (Q7). Will likely split by test tier at its plan. | The remaining slotted objects go: census about 0.2-0.25 GB (est.). | Exact. Collected-id diff listed. |
| N8 | `branch-point-storage` | Alone, per D5. `BranchPoint` keeps its class and API and gets a store row plus a slot table. Retires roles 4 and 5 with evidence; deletes its dead fields and the stale docstring. | About 0 bytes. This slice is about correctness and clarity. | Exact. |
| N9 | `array-namespace` | Hot paths (map, refine, merge, segment build, crossing resolution) take an array namespace: numpy by default, CuPy when `enable_gpu` is on (D6). No GPU wiring beyond today's opt-in. | None on CPU; this is preparation. | Exact on numpy. |
| N10 | `save-load-and-after-profile` | `TangleWorkbench.save(path)` / `load(path, map, map_inv, jacobian)` (Q5) with a round-trip test, then the after-profile drawn beside the baseline (D2). | Cluster-ready state; the final before/after figures. | The round trip is exact. |

**Every slice's rubric includes:**

- the suite green, with `KNOWN_ISSUES` = C3 alone and no XPASS;
- the identity probe (`stress/profile/identity_probe.py`), base vs branch, either identical or every
  difference listed per D3 (case, points, max coordinate and cdist difference, cause);
- topology outputs that move listed case by case and judged;
- `profile_numerics.py --cases scale --modes time` (1e-7, about 10 min) rerun, with the census at the
  plateau stated against the baseline (1913.8 MB traced at the 1e-7 plateau; the census column of
  `figures/numerics_profile/a8bb8a0/scale_time.json`). This is how the "expected effect" column is
  checked, not trusted;
- `INNER_STAGES` in the profiler updated for any wrapped function the slice renames or deletes (the
  same label, so the figures stay comparable);
- the probe's readers (`curve_arrays`, `curve_extent`, `bracket_cdists`, `image_curves`,
  `key_text`) edited, in the slice's own diff, when the slice changes an accessor they call.

### N1 — `slim-segment-index`

- **Scope.** A blast indexes each image curve query-only (`Tangle.add_manifold(iterated,
  index_segments=False)`, `BridgeIterator.py:359`), yet each segment keeps a `_Segment`, a frozenset
  edge key of two `id()` ints in `_edge_to_sid`, an owner set in `_seg_manifolds` and a membership in
  the curve's `_manifold_segs`; the bridges cut from the curve claim the same segments again
  (`Tangle.py:966-968`). After the blast nothing queries them geometrically: they serve ownership and
  dedup (`release_manifold`, the edge dedup of `add_manifold`, `compute_intersections`' "curve adds
  no segment" test). Move that to the curve level, keyed by int point identities; keep both rtrees.
- **Files.** `numerics/Tangle.py`, `numerics/BridgeIterator.py`, maybe `numerics/TangleWorkbench.py`,
  `tests/unit/numerics/test_segment_index_bookkeeping.py`.
- **Rubric sketch.** Probe identical on all 58 builds (registry ids included); the 1e-7 plateau census
  drops by the segment-index band; one rtree per stability still asserted (u×s pairs only).
- **D3.** Expected exact. The risk is which crossings dedup, hence the registry ids.
- **Effect.** Up to 2.58 of 4.68 GB at the 1e-9 plateau (the index is 2.60 GB: `_edge_to_sid` 331
  B/seg, owner sets 368 B/seg, `_Segment` 100 B/seg, `_seg_lookup` 51 B/seg).

### N2 — `same-stability-count-new-only` (Q8)

- **Scope.** `TangleWorkbench.report_resolution(new_bridges=)` ->
  `numerics/diagnostics.count_same_stability_crossings` (:200) -> `_hits` (:314) ->
  `_candidate_pairs` (:401) builds rows for every unstable curve of the fixed point and sweeps all of
  them each blast; the `new` filter is applied only after the pairs are generated. Generate pairs
  only for pairs holding a new segment.
- **Files.** `numerics/diagnostics.py`, its unit test (`test_same_stability_diagnostic.py`).
- **Rubric sketch.** The counts per blast equal on the series (probe identical, plus the
  `ResolutionReport.same_stability_count` of every `blast_zone(diagnose=True)`, base vs branch); blast time at
  the 1e-7 plateau measured.
- **D3.** Exact.
- **Effect.** Up to -444 s of 541 s of blast time at 1e-9 (diagnose 444 s = 82 %; 398.6 of 458.4 s on
  the 1e-7 series). Est.: the cost is measured, the new cost is not.

### N3 — `point-store`

- **Scope.** A new `numerics/PointStore.py`: append-only rows `x, y, cdist, image, preimage` (int64,
  -1 none), a row id never reused. `Point` becomes a slotted object over its row and keeps its
  `forward` / `backward` links as an adapter for N4. Phantom preimages become rows on no curve
  (hard place 1); the tie rule's "has an iterate" reads `image >= 0 or preimage >= 0` and must count
  a phantom row as today.
- **Files.** new `numerics/PointStore.py`, `numerics/BasePoint.py`, `numerics/Point.py`,
  `numerics/ManifoldMachine.py` (`_cache_preiterate`, `_get_preiterate`, the tie rule),
  `numerics/TangleWorkbench.py`.
- **Rubric sketch.** Probe identical; the 15 of k28's 17 ties whose existing point holds a phantom
  keep their keeper; census per `Point` measured; walk time compared with the baseline.
- **D3.** Expected exact; low risk; the trap is the tie rule.
- **Effect.** Est. -1.3 to -1.5 GB (`Point` 400 B -> about 120 B; a phantom 376 B -> 40 B; 1.82 M
  phantoms and 3.28 M real points at the 1e-9 plateau).

### N4 — `curves-on-order-arrays`

- **Scope.** Each curve holds an int `ids` array in curve order. Grow, merge and refine work on
  positions: merge = `np.searchsorted(cdist_existing, cdist_new, side="right")` plus "existing wins"
  (hard place 2), refine = positional insert (`np.insert(ids, pair_index + 1, new_ids)`), never a
  sort. `_image_chains` reads a per-curve position index (hard place 3). Unit tests pin the two tie
  branches that never fired (`ManifoldMachine.py:570-589`). The stale merge comment
  (`ManifoldMachine.py:552-555`) and the per-point `stretch_param` (`Point.py:86-89`) are fixed here.
- **Files.** `numerics/ManifoldMachine.py`, `numerics/BaseManifold.py`,
  `numerics/ManifoldInitializer.py`, `numerics/ManifoldView.py`, the refinement tests. May split at
  its own plan.
- **Rubric sketch.** Probe identical; the 19 measured ties (k28 x17, heteroclinic x2) reproduce;
  `tests/regression/test_blast_through_a_shared_straddle.py` kept (the only test of several images on
  one curve) plus the array invariant; refine and walk time compared.
- **D3.** Target exact.
- **Effect.** Time: refine plus walks are 53 s of 541 at 1e-9 (walks 31, refine 22). Removes the
  forward/backward links.

### N5 — `segments-and-brackets-on-ids`

- **Scope.** Segments become `(ids[:-1], ids[1:])`; `Intersection.unstable_segment` /
  `stable_segment` become store id pairs (hard place 4). Readers: `Tangle.create_bridges` (:883, the
  separator logic :918-927), `IterateInference._brackets` (:67-74) and `_map_gated` (:281-290),
  `IntersectionRegistry._polyline_gap` (:48-62, a shared node by identity -> id equality), and the
  bracket readers in `tests/helpers/laws.py` (:578-580, :633-661), which the reviewer diffs line by
  line. `Tangle.populate_intersections_for_manifold` (:545-572), which scans the whole
  `_intersecting_segments` set per bridge, is looked at here.
- **Files.** `numerics/Tangle.py`, `numerics/Intersection.py`, `numerics/IntersectionRegistry.py`,
  `numerics/IterateInference.py`, `tests/helpers/laws.py`, the probe's `bracket_cdists`.
- **Rubric sketch.** Probe identical (the bracket cdists are in every registry row); the law helpers'
  diff read line by line.
- **D3.** Exact: the same floats are read.
- **Effect.** Small: the bookkeeping left after N1.
- **Split (2026-10-25, the N5 plan's Q7).** N5a `brackets-on-rows`: `Intersection` holds row
  pairs (the periodic point's `BranchPoint` in place of a row) and its `store`, every bracket reader
  moves to `segment_ends` / `brackets`, and the iterate search reads each bracket once and masks one
  bisected slice; no segment, edge key, sid or detection changes. N5b
  `segments-and-edge-keys-on-rows` (next): `_Segment` holds two rows (or the root `BranchPoint`)
  plus `manifold`, `id`, `in_rtree`; edge keys become `((lo + 1) << 32) | (hi + 1)` ints (a
  `BranchPoint` end as -1, never an index) for `_edge_to_sid` and `_query_edges`; the `Tangle`
  records its store on first add and asserts it after; the query-only `add_manifold` runs on arrays
  with one `rtree.intersection_v` per curve and keeps a `_Segment` only for an edge with a pair; the
  cut's edge registration is vectorised. N5b must not change the sid counter or its order,
  `frozenset((new_sid, cand_id))` in that element order, or the `_resolve_crossing_pair` operand
  order (`tuple(frozenset(pair))` decides float bits). Its files: `Tangle.py`,
  `tests/_tools/resolution_probe.py` (the `_resolution_ratio` wrapper),
  `test_tangle_index_and_orientation.py`, `test_segment_index_bookkeeping.py`,
  `test_same_stability_discard.py`, and `profile_numerics._tangle_bytes` if field names change.

### N6 — `readers-on-views`

- **Scope.** Every reader outside the engine moves to the D4 accessors (`curve.ids` / `points` /
  `cdists`, `root` / `tail`, `bridge.head` / `bridge.tail`, `manifold.end`): `Bridge` (and
  `joined`), `BridgeIterator`, `TangleWorkbench` (trim, anchors; `grow_until_turnaround` reading
  `_coords`, `TangleWorkbench.py:768-792`), `numerics/diagnostics.py`, `numerics/geometry.py`,
  `loom/ResonanceZone.py`, and the six topology lines (`StablePartition.py:835`,
  `topology/diagnostics.py:91,93-99,208`, `plotting.py:431,483,496-505`). The walk memo goes.
- **Files.** Those named.
- **Rubric sketch.** Probe identical; topology's diff limited to the six lines; no topology module
  calls the map (unchanged).
- **D3.** Exact.
- **Effect.** -0.075 GB (the walk memo) plus the node-to-array conversion in `_segment_rows`
  (~1 s of an 18 s blast at 1e-7).

### N7 — `delete-point`

- **Scope.** `Point` and `BasePoint` are deleted and leave `__all__` (`numerics/__init__.py:19,48`,
  `tanglepack/__init__.py:24,73`). The tests that read nodes (75 functions in 31 files, section E.2),
  `tests/helpers/{invariants,laws,polylines}.py`, `tests/cases.py:1052`, the stress and script
  touchpoints move (Q7).
- **Files.** `numerics/BasePoint.py`, `numerics/Point.py`, both `__init__.py`, the tests of E.2,
  `stress/numerics/{e07_breaking,maps}.py`, `scripts/refactor_proof/hole-side-and-own-blast/hole_sides.py`,
  the probe's `curve_arrays`. Will likely split by test tier at its plan.
- **Rubric sketch.** Probe identical; collected-id diff listed (deleted, rewritten, added); law helper
  diffs read line by line.
- **D3.** Exact.
- **Effect.** The remaining slotted objects go: census about 0.2-0.25 GB at the 1e-9 plateau (est.).

### N8 — `branch-point-storage` (D5: alone)

- **Scope.** `BranchPoint` keeps its class and API and gets a store row plus a slot table
  (section G). Roles 4 (the anchor's iterate cycle) and 5 (the walk-through toggle, a latent trap)
  retire with the probe's evidence; the dead fields (`edists`, `forward_stretch_params`,
  `backward_stretch_params`, `BranchPoint.py:91-97`) and the stale docstring (:30-36, crossings as
  BranchPoints) go. Q9: both `num_branches` keep their names and meanings, documented.
- **Files.** `numerics/BranchPoint.py`, `numerics/FixedPoint.py`, `numerics/FixedPointSolver.py`,
  `numerics/ManifoldInitializer.py`, `numerics/BaseManifold.py`, their tests.
- **Rubric sketch.** Probe identical on the inversion cases above all (roles 1, 3, 8, 9); every role
  of section G checked off as kept, moved or retired with evidence.
- **D3.** Exact.
- **Effect.** About 0 bytes. Correctness and clarity.

### N9 — `array-namespace`

- **Scope.** Hot paths (map, refine, merge, segment build, crossing resolution) take an array
  namespace: numpy by default, CuPy when `enable_gpu` is on (D6). No GPU wiring beyond today's opt-in.
- **Files.** `numerics/gpu.py` and the hot-path modules of N3-N5.
- **Rubric sketch.** Probe identical on numpy; `enable_gpu` tests unchanged.
- **D3.** Exact on numpy.
- **Effect.** None on CPU; preparation.

### N10 — `save-load-and-after-profile`

- **Scope.** `TangleWorkbench.save(path)` / `load(path, map, map_inv, jacobian)` (Q5): what a save
  holds is listed in the findings' Q5 paragraph (the store; per curve key, ids, `end`,
  `stretch_param`, version; bridges with `iterated`, `dropped_at`, partial pieces, `_image_curves`,
  `_half_iterated_bridges`; the registry with brackets as ids and `unstable_manifold` as a curve index,
  its orders, the iterate table, `unresolved_links`, generation counters, `_reported_flags` /
  `_reported_artifacts`; fixed points with the BranchPoint slot table; the machine's settings). Then
  the after-profile: the README's command, unchanged, drawn beside the baseline with `--compare
  figures/numerics_profile/a8bb8a0`.
- **Files.** `numerics/TangleWorkbench.py` (or a `numerics/persist.py`), a round-trip test, the
  after-profile folders under `figures/numerics_profile/`, and `stress/profile/profile_figures.py`:
  `figure_stage_peaks` takes no baseline, so `--compare` does not yet draw the stage peaks beside it.
- **Rubric sketch.** Round trip gives equal `compute_intersections` and `topology()`; probe identical;
  the after-profile figures show baseline and branch on one axes.
- **D3.** The round trip is exact.
- **Effect.** Cluster-ready state; the final before/after figures.


---

## D. The D4 vocabulary (findings 4.3; the slices settle the details)

One `PointStore` per workbench (append-only rows; a row id never changes or is reused), and views:

- `store.xy[id]` (float64 `(N, 2)`), `store.cdist[id]`, `store.image[id]` / `store.preimage[id]` (int,
  -1 none: the map's forward and backward iterate, exactly today's `next_iterate` / `prev_iterate`).
- `manifold.ids` (int array, curve order, root first), `manifold.points` = `store.xy[ids]`,
  `manifold.cdists`, `len(manifold)`, `manifold.root` / `manifold.tail` and `bridge.head` /
  `bridge.tail`: a small immutable record `(id, x, y, cdist)` (a VALUE, not a node: no `forward`, no
  links -- D4's "no PointRef shim"), `manifold.end` (the extent a trim or a zone sets; `tail` is the
  record at it), `manifold.iterated` (bool mask: the row has its image in the curve's own dynamical
  direction), `manifold.segments` (`(ids[:-1], ids[1:])`), `bridge.span = (manifold, lo, hi)` or
  `bridge.ids`.
- `intersection.unstable_segment` / `stable_segment`: `(id_lo, id_hi)`.
- `fp.branch_points[o]` stays a `BranchPoint` (D5): `num_branches`, `get_cdist(stability)`,
  `x`/`y`/`get_point()`, `branch(stability, b)` -> the curve in that slot (replacing
  `forward_branches[b]` / `backward_branches[b]`), `store_id` (its row).

---

## E. Touchpoint inventory (findings 4.1-4.2)

Grep (rubric 3's pattern) over the TRACKED files of `src/ tests/ scripts/ stress/` at `a8bb8a0`:
`\.forward\b|\.backward\b|next_iterate|prev_iterate|walk_fwd|walk_back|return_nodes|\.root\b|\.tail\b|BranchPoint|forward_branches|backward_branches`
gives **60 files** (src 23 files / 301 lines, tests 32 / 220, stress 4 / 16, scripts 1 / 1). Every hit is below, by
file and line, with its role and what it becomes under D4. Roles: **walk** (read a curve's nodes in
order), **splice** (insert / relink nodes), **iter** (iterate links), **bracket** (a crossing's two
nodes), **anchor** (the `BranchPoint` and slot logic), **tail** (setting / reading the extent),
**class** (the classes themselves, imports, exports, type hints), **doc** (comment / docstring only),
**unrelated**. The planner's D4 vocabulary used in the right column is proposed in section D.

A second grep for the node API the rubric pattern misses (`get_point()`, `get_cdist(`,
`insert_point_`, `unstable_segment`, `stable_segment`, `first_node(`, `_walk_nodes`, `_find_tail`,
`._coords`, `get_(non_)iterated_*`, `exists_*_iterate`, `get_*_iterate`, `insert_*_iterate`) adds the
sites marked (+).

### E.1 `src/` (23 files)

| file | lines | role | becomes (D4) |
|---|---|---|---|
| `numerics/BasePoint.py` | 6,37,38,64,65,87,118,137,140,142,145,146,148,167,170,172,175,176,178,196,200,218,222,230,235,243 | class + iter | deleted; `store.image[id]` / `store.preimage[id]` columns (int64, -1 none) |
| `numerics/Point.py` | 27,38,39,70,71,83,91,93,94,96,97,109,117,119,120,122,123 | class + splice | deleted; `curve.insert(position, ids)` on the order array |
| `numerics/BranchPoint.py` | 4,26,34,35,38,48,49,55,56,71,93,96,105,125,131,135,137,138,140,141,152,158,162,164,165,167,168 | class + anchor + splice | STAYS a class (D5), its own slice (N8); see section G |
| `numerics/BaseManifold.py` | 42,54,55,79,84,93,99,116,117,124,147,152,157,162,174,177,180,198,201,209,215,219,220,224,227,230,231,238,242,247,250,251,253,258,259,264,265,274,286,290,314,321,324,332,347,352,356,361,364,382,390,398,408,436,445,458,461,466,472,476,481,485,487,491,496,500,505,508,541,544,548,549,555,556,563,564,566,567,572,579,582,589,590,596,597,604,605,607,608,613,620 | walk (`_walk_nodes` :182, `_collect` :267, getters :323-459, `walk_fwd`/`walk_back` :461-508, `_branch_forward`/`_branch_backward` :540-620), tail (:147-180 setters, `_find_tail`), anchor (`first_node` :222) | `Manifold` view: `ids`, `points`, `cdists`, `len()`, `iterated_mask`; `root`/`tail` = end records; `end` (an extent, see "tail" below); the five getters become one `select(mask)` |
| `numerics/ManifoldMachine.py` | 97,226,287,291,300,335,344,346,366,368,409,413,462,503,504,510,514,518,519,524,525,531,532,540,544,545,546,562,563,566,567,571,572,578,581,582,588,597,602,646,647,650,652,691,729,730,796,804,809,811,815,824,833,835,839,840,859,865,867,873,881,885 (+ 360, 649) | walk (:344-346 node lists, :462 `_image_chains`, :646-652 layer 0), splice (merge :489-604, `_insert_point_geometrically` :838-887, `_manifold_from_points` :324-336), iter (:366-368 link images, `_get_preiterate` :796-811, `_cache_preiterate` :814-835), anchor (:226 `_branch_view`, :865-875 slot insert), tail (:335, :597, :602) | the engine on arrays: `map(store.xy[ids])`, `store.append(xy, cdist, preimage=ids)`, merge = positional `searchsorted` + tie rule, refine = positional insert; hard places 1-3 |
| `numerics/ManifoldInitializer.py` | 179,180,182,188,189,191,255,256,260,262,269,395,406,410 | iter (the "fictitious" k-step iterate chain :172-191, :255-269), anchor (:273-280 slot insert, :395, :406, :410 roots) | rows appended before they sit on a curve; slot assignment through the `BranchPoint` API |
| `numerics/ManifoldView.py` | 34,35,39,40 | walk aliases | `view.ids` / `view.points`; `walk_fwd` aliases go |
| `numerics/Tangle.py` | 62,99,100,380,386,387,391,392,713,722,724,731,896,897,905,1006 (+ `_Segment.bounds` :106-115, `intersects` :117-124, `_find_true_intersection` :418-420, `_resolve_crossing_pair` :293-294, `_cdist_between` :1037-1047) | walk (`_segments_of` :710-738), bracket (:293-294), anchor (`_outer_neighbours` :375-392, a BranchPoint has no neighbour, :896-897 `manifold.root is fp.branch_points[o]`, :1006 separator after a BranchPoint), splice (`_insert_crossing_separator` :972-1010) | segments = `(curve, position)` or `(id_lo, id_hi)` arrays; brackets = id pairs; `curve.starts_at_anchor`; separator = `store.append` + `curve.insert` |
| `numerics/Intersection.py` | (+) 49,54,55,99,102,119,120,137,140,153,156 | bracket | `unstable_segment: tuple[int, int]` / `stable_segment: tuple[int, int]` (store ids) |
| `numerics/IntersectionRegistry.py` | (+) 55,60,466 | bracket (`_polyline_gap` :48-62 shares a node by identity) | id equality |
| `numerics/IterateInference.py` | 37 (doc) (+) 69-71, 281, 290 | bracket (cdists and coords of the bracket nodes) | `store.cdist[ids]`, `store.xy[ids]` |
| `numerics/TangleWorkbench.py` | 784,786,1628,1629,1759,2091,2103 (+ 391, 710, 769, 771-772 and 788-789 `._coords`, 1632, 1666-1667, 2090) | walk (:1628-1632 anchor nodes, :2090 trim scan), tail (:710 `_find_tail` in `grow_n_times`' finally, :2103 trim, :784-786 tail and the node before it), anchor (:768-769 `branch_points[0]`, `first_node`), bracket-read (:1759 bridge root/tail cdists) | `curve.points[:2]`, `curve.points[-2:]`, `np.searchsorted(curve.cdists, c)` for the trim, `curve.end = ...` |
| `numerics/Bridge.py` | 82,151,155,191 | doc + tail (`joined` :127-161 builds from `first.root` / `last.tail`) | a bridge = `(manifold, lo, hi)` or its own ids; `joined` = the span `first.lo .. last.hi` |
| `numerics/BridgeIterator.py` | 344,348 | tail (wrap the machine's image as a Bridge from its root / tail) | `Bridge(ids=image.ids, ...)` |
| `numerics/diagnostics.py` | 263,278 (+ 385) | walk (`_curves` node lists, `_segment_rows` :381-398 converts nodes to arrays) | `curve.points` directly (the conversion disappears) |
| `numerics/geometry.py` | 80 (+ 83, 85) | walk (`arc_polyline` filters nodes by cdist) | `curve.points[mask]`, `mask` from `curve.cdists` |
| `numerics/FixedPoint.py` | 5,32,37,71,119,127 | class / anchor (:127 the `BranchPoint` per orbit point) | unchanged until the BranchPoint slice |
| `numerics/FixedPointSolver.py` | 179,180 (+ 165-168) | anchor + iter (coords, `cdist = 0.0`, the orbit cycle of `next_iterate` / `prev_iterate`) | BranchPoint slice |
| `numerics/__init__.py` | 19,48 | class export | `Point` / `BasePoint` leave `__all__` with N7; `BranchPoint` stays |
| `tanglepack/__init__.py` | 24,73 | class export | as above |
| `loom/ResonanceZone.py` | 30,67,114,364,365,516,521 | tail (zone trim captures `manifold.tail` :516-521, restore compares stable cdists and resets it :364-365), class (type hints :67, :114) | `previous_ends: dict[key, int]` (a point id or a position); `curve.end_cdist` |
| `topology/StablePartition.py` | 835 | tail-read (`tail.get_cdist("stable")`, the end of a branch with no crossing) | `curve.tail.cdist` (end record) |
| `topology/diagnostics.py` | 91,208 (+ 93-99) | walk (nodes, `node.cdist`, `node.get_point()` near a crossing) | `curve.points`, `curve.cdists` |
| `topology/plotting.py` | 431,483 (+ 496-505) | walk (stable-arc midpoint, bridge midpoint and tangent) | `curve.points`, `curve.cdists` |
| `topology/Arrangement.py` | 504,536,581 | **unrelated** (half-edge `twin.tail` / `edge.tail` / `edge.head`, ints) | none |

`topology/*.py` uses of `trellis.intersection(i).get_point()` (`Arrangement.py:616`, `DualGraph.py:891`,
`TopologyResults.py:557,560`, `topology/diagnostics.py:353,377`, `plotting.py:345,423-424,446-450`) are
`Intersection.get_point()`, **unrelated** to the lists. So topology's whole contact with the linked list
is 6 lines: `StablePartition.py:835`, `topology/diagnostics.py:91,93-99,208`, `plotting.py:431,483,496-505`;
all six read a curve's coordinates / cdists in order, i.e. `curve.points` / `curve.cdists` /
`curve.tail.cdist`. Loom's is `ResonanceZone.py` (tails) plus `Blast.py` through `get_point_array()`
(arrays already).

### E.2 `tests/`, `scripts/`, `stress/` (tracked)

| file | lines | role |
|---|---|---|
| `tests/unit/numerics/test_linked_list.py` | 1,3,6,16,24,36,37,42,43,47,48,49,50,62,63,64,65,74,75,76,77,87,88,98,99,102,106,112-115,120,126-130,135,141-144,149,155-158,164,165,168,170-173,177,178,181,183-186,190,193,196,200,203,206 | class + splice + iter: the whole file pins the lists (14 test functions, 18 ids) |
| `tests/unit/numerics/test_cleanup_walkers_and_examples.py` | 11,32,42,45,47,55,57,63,65,82,88,91,106,115,148,151,166,171,181-183,187,188,191,197,203,204,208-210,217,234 | walk + anchor (`first_node`, slots :170-210); 7 test functions |
| `tests/unit/numerics/test_manifold_initializer.py` | 87,101,102,113,115,119,132,134,167,169,173,202 | anchor (`root is branch_points[o]`) + iter |
| `tests/unit/numerics/test_inversion_fixture.py` | 87,88,92,112,114,116 | anchor + walk |
| `tests/unit/numerics/test_bridge_identity.py` | 453,458,459,463,485,512 | tail / walk (bridge root/tail nodes) |
| `tests/unit/numerics/test_stable_trim.py` | 15,33,46,53,58,78,79 | tail (trim sets `tail`) |
| `tests/unit/numerics/test_generation_and_caches.py` | 197,244,248,252 | walk memo (`return_nodes`) |
| `tests/unit/numerics/test_initializer_cdist.py` | 29,32,95,101 | iter + walk |
| `tests/unit/numerics/test_fixed_point.py` | 163,170,171 (+ 14 `branch_points[0].num_branches == 2`) | anchor (slot counts) |
| `tests/unit/numerics/test_geometry.py` | 160,174,217 | walk |
| `tests/unit/numerics/test_segment_index_bookkeeping.py` | 59,65,71 | walk + `_Segment` internals |
| `tests/unit/numerics/test_single_source_of_truth.py` | 58,63,68 | bracket / walk |
| `tests/unit/numerics/test_tangle_index_and_orientation.py` | 37,39,41 | walk |
| `tests/unit/numerics/test_workbench_bugfixes.py` | 38,42,66 | walk / tail |
| `tests/unit/numerics/test_invariant_helpers.py` | 47,57 | walk |
| `tests/unit/numerics/test_refinement.py` | 23,27 | walk |
| `tests/unit/numerics/test_same_stability_diagnostic.py` | 26,52 | class (builds a `BranchPoint(2, (0,0), 5.0, 5.0)` root by hand) |
| `tests/unit/numerics/test_resolved_crossings.py` | 157 | tail-read |
| `tests/unit/loom/test_blast_restore.py` | 147,148,173,183,184 | tail (zone restore) |
| `tests/unit/loom/test_blast_no_overlap.py` | 21 | walk |
| `tests/unit/topology/test_arrangement.py` | 82,83,427 | splice (hand-built two-node lists `root.forward, mid.backward = ...`) |
| `tests/unit/topology/test_arrangement_sparse.py` | 30,31 | splice (same) |
| `tests/unit/topology/test_stable_partition.py` | 120,521,522,701,702 | fake `get_point_array(return_nodes=)` stubs |
| `tests/unit/topology/test_topology_invariants.py` | 239,240,258 | fake `get_point_array(return_nodes=)` stubs |
| `tests/unit/topology/test_enclosure.py` | 110 | **unrelated** (a test name containing "next_iterate") |
| `tests/regression/test_boundary_straddle_cdist.py` | 54,59,63,67,68,73,74,102,103 | walk + iter |
| `tests/regression/test_blast_through_a_shared_straddle.py` | 28,36 | walk |
| `tests/regression/test_foreign_chords.py` | 83 | walk |
| `tests/helpers/invariants.py` | 32,34,46,54,60,63,268,269,277,304,305,309,310 (+ 162,197,233) | walk + iter + anchor (`node_cdist` special-cases `BranchPoint.cdists`) |
| `tests/helpers/laws.py` | 186,196,519,524,981,982,1625,2964,2966 (+ 187-201, 578-580, 633-635, 660-661, 1626-1627) | walk + bracket + tail: **the law tier's own helpers read nodes** |
| `tests/helpers/polylines.py` | 16,23 | class (builds a manifold from a polyline with a `BranchPoint` root) |
| `tests/cases.py` | 1052 | tail-read (`build_near_tangency` grows until `tail.cdist >= reach`) |
| `tests/_tools/resolution_probe.py` | (+) 50,51,63,69,76,128,133 | bracket + `_Segment` (a tool, not collected) |
| `stress/numerics/e07_breaking.py` | 189,376 | tail-read |
| `stress/numerics/maps.py` | 194,245 | walk |
| `stress/numerics/results/e02_gpu_end_to_end{,_quick}.jsonl` | 19-24 / 7-12 | **unrelated** (recorded traceback text) |
| `scripts/refactor_proof/hole-side-and-own-blast/hole_sides.py` | 74 | walk |

Untracked (staying so, the author's): `scripts/henon_k25_period3_class_focus.py:99`,
`scripts/henon_k25_period3_class_iterates.py:138,140,184,501`. No tracked notebook hits the rubric
pattern; `scripts/manifold_cutting_test.ipynb` and five tracked demo scripts
(`henon_period_3_manifolds.py`, `henon_basic_stable_and_unstable.py`, `henon_initial_intersections.py`,
`henon_period3_zone_sweep.py`, `henon_k10_dual_graph.py`, `tangle_workbench_test.py`) call
`get_point()` / `get_cdist()` on nodes (second grep).

Test churn, counted by AST (a test function whose body, or a same-file helper it calls, matches either
grep): **75 test functions in 31 test files** (of 1260 collected ids at `a8bb8a0`), plus the helpers in
`tests/helpers/{invariants,laws,polylines}.py` that many more tests reach. The laws in `laws.py`
(`check_bridge_endpoints_on_own_branch` :519-524, `check_crossing_cdist_bracketed` :578-580, the
bracket/gate readers :633-661, :981-982, :1625-1627, :2964-2966) are law-tier helpers: D3's "law
tests unchanged" is true of the test functions, NOT of these helpers, which must change in the slice
that removes `Point` (their logic must stay a faithful statement of the law; the reviewer should diff
them line by line).

---

## F. The five non-mechanical places (findings 5)
Counts are on the 11 LAW_CASES builders plus `k28_deep_blast_sliver`, `p3_6_blasts` and the
scaling case (12 blasts), instrumented in a scratch probe (`tie_probe.py`), verified.

**1. The refined point's phantom preimage** (`ManifoldMachine._cache_preiterate` :814-835).
Today: every refined point carries a `Point(x, y)` with `cdist` None in `prev_iterate` (unstable) /
`next_iterate` (stable), on no list. It costs 376 B (measured, tracemalloc, 1e5 objects) and there are
MANY: 1457 phantoms beside 1733 curve points on k10; 6471 beside 7858 + 2814 on nested; 1.82 M beside
3.28 M real points on the 1e-9 scaling case (36 % of all `Point` objects; section B). Readers: the next
refinement layer (`_get_preiterate` :796-811 for `pre_mid`), and -- easy to miss -- the merge tie rule
(:561-569), which counts a phantom as an iterate link (15 of k28's 17 ties have the existing point
holding a phantom).
- Array form (a): two columns `pre_x, pre_y` (NaN when the row has a real preimage or none): +16 B
  on every row, no extra row, the preimage of row `i` is `store.xy[preimage[i]]` if
  `preimage[i] >= 0` else `(pre_x[i], pre_y[i])`.
- Array form (b): phantom ROWS (cdist NaN, on no curve, `preimage[i]` points at them): 40 B (or 24 B
  with int32 links) per phantom only; the "row on no curve" kind also covers the two other off-curve
  rows that exist today (the initializer's k-step "fictitious" chain :172-191, adopted onto curves
  later; merge-dropped duplicates, kept alive by their preimage's `next_iterate`).
- Byte-identity risk: LOW for both. The midpoint `0.5 * (pre0 + pre1)` (:738) and the map of it
  are elementwise float ops, identical whether the operands come from objects or columns. The one
  trap is the tie rule: under (a) "has an iterate" must read `image >= 0 or preimage >= 0 or
  not isnan(pre_x)`, or a tie flips its keeper.

**2. The splice at cdist ties** (`merge_manifolds` :489-604, `_refine_layer` :671-786).
Measured: merge ties between DISTINCT nodes occur on k28_one/two_blasts (17 per build) and heteroclinic
(2), nowhere else; every one is in GROWTH (no blast calls `merge_manifolds` on any case measured:
`merge_calls` does not move during blasts), and in every one BOTH nodes carry iterate links
(`iter_11`), so the rule always keeps `head_1`, the node of the EXISTING image chain
(`merge_manifolds(chain, mapped)`, :383-384). The comment at :552-555 ("The newly generated point has no iterate links yet") is wrong in practice: the new image already holds `prev_iterate` (:366) before the
merge. The dropped node stays alive as its preimage's `next_iterate` (on no curve). Inputs are never
unsorted; no input list has an internal tie. Refinement midpoints that round onto an endpoint
(`mid_cd == c0 or c1`): **0** of ~0.6 M refined points (law cases + scaling + k28 sliver).
- Array form: merge = `np.searchsorted(cdist_existing, cdist_new, side="right")` -- new points go
  AFTER equal existing ones -- plus the tie rule: a new point whose cdist equals an existing point's is
  left OFF the curve (an off-curve row whose preimage's `image` still names it, exactly today's
  dropped duplicate), and the same-physical-node case (a row id in both
  inputs) deduplicated by id. Refinement = positional insert (`np.insert(ids, pair_index + 1,
  new_ids)`), never a sort: identity with today holds whatever the cdists do.
- Byte-identity risk: LOW on the measured corpus (19 ties, one keeper rule, all reproduced by
  "existing wins"), MEDIUM in general: the both-unlinked / only-new-linked branches (:570-589) never fired,
  so their array translation is untested by data; a unit test must pin them. A stable sort by cdist
  (Q3 alternative) is the same order as the merge whenever no tie exists, and differs exactly at
  those 19 ties (keeps both nodes -> the "spike" the memory `cdist-strict-monotonicity-fix` warns of).

**3. Several images on one curve** (`_image_chains` :421-466).
Measured: 0 calls returning more than one chain and 0 "short" chains (an image the walk skipped) over
every case above, including 382 calls on the scaling case. The two-curve situation the docstring
describes (WP10, a bridge sharing straddle nodes with an already-iterated neighbour) is covered by
`tests/regression/test_blast_through_a_shared_straddle.py`, not by any law case.
- Array form: each image id's `(curve, position)` from a per-curve inverse index (`pos_of[id]`,
  rebuilt lazily per curve version); a chain is a maximal run of images on one curve with
  consecutive positions in walk order bounded by the last image's cdist; ids on no curve (merge-
  dropped duplicates) are skipped as today.
- Byte-identity risk: LOW (it decides WHICH merges happen, not float values), but only the
  regression test exercises it; the slice must keep that test and add the array invariant.

**4. Brackets that must stay valid** (`Intersection.unstable_segment` / `stable_segment`,
`Tangle.py:293-294`). Readers: `create_bridges` (:883 reads the pair, the separator logic
:918-927 compares bracket nodes by identity), `IterateInference._brackets` :67-74 (the bracket cdists)
and `_map_gated` :281-290 (their coordinates), `IntersectionRegistry._polyline_gap` :48-62 (a SHARED
node, by identity), `tests/helpers/laws.py:578-580, 633-661`. Brackets never hold a `BranchPoint` (0
on every law case: the root-segment detections dedup into the anchor, which has no bracket).
- Array form: `(id_lo, id_hi)` store ids. Identity, cdists and coordinates read exactly as today,
  and they stay valid across growth, refinement and separators because a row id never moves; the
  CUT needs positions, found at cut time by `pos_of`.
- With per-curve arrays (Q2 alternative) a bracket is `(curve, position)`, and every insertion
  before it shifts it: each refinement and separator would have to remap every bracket on that
  curve, and a bracket would also have to name the curve object that `create_bridges` groups by. The
  shared store needs none of that. Byte-identity risk: NONE for the store form (the same floats are
  read); the per-curve form risks a stale position silently cutting one node off.

**5. `BranchPoint` is shared**: section G.

---

## G. Every role `BranchPoint` plays at `a8bb8a0` (findings 6; D5)
Verified by a scratch probe on k10, p3, inversion, inversion_period2, nested, heteroclinic
(`bp_roles.py`, `bp_probe.py`) and by reading; "believed" where marked.

1. **The one root of up to four branches.** One per orbit point (`FixedPoint.py:126-128`). Unstable
   branches hang in `forward_branches[b]`, stable ones in `backward_branches[b]`
   (`ManifoldInitializer.py:273-280`), `b` the EIGENDIRECTION (slot) index. Every registered manifold
   is rooted at one (8/8 on nested, 4/4 on inversion); an inversion point fills all four slots
   (`(2, 2)` on inversion and inversion_period2), any other point `(1, 1)`.
2. **Two slot counts that differ.** `BranchPoint.num_branches` is ALWAYS 2 (`_EIGENDIRECTIONS`,
   `FixedPoint.py:37-40,127`; pinned by `test_fixed_point.py:14`), while `FixedPoint.num_branches`
   is 1 or 2 (branches one chain visits). The author's `point.num_branches` must say which it means.
3. **The periodic point's coordinates and cdists.** `x`, `y`, `_coords` set by the solver
   (`FixedPointSolver.py:165-166`), a SCALAR `cdist = 0.0` (:168; read by every `node.cdist` loop
   over a curve, e.g. `topology/diagnostics.py:92`, `plotting.py:432`) AND a tuple `cdists = (0.0,
   0.0)` (`FixedPoint.py:127`) read by `get_cdist(stability)` (`Tangle._cdist_between` :1047-1048 on a
   root segment, `TangleWorkbench.py:1759`, `geometry.py:83`, `tests/helpers/invariants.py:63-66`).
   `edists`, `forward_stretch_params`, `backward_stretch_params` (`BranchPoint.py:91-97`) are never
   read anywhere (grep: 0 hits outside the class) -- dead.
4. **An iterate cycle.** `next_iterate` / `prev_iterate` link orbit point `i` to `i +- 1`
   (`FixedPointSolver.py:179-180`), so the anchor counts as "iterated" in `_collect`'s filters and
   would give refinement a preimage at the anchor. Believed unexercised: no growth view starts at a
   BranchPoint (`_branch_view` steps off it, `ManifoldMachine.py:226`), no anchor-rooted bridge is ever
   mapped on the law cases (its image already exists, `BridgeIterator._existing_image_bridges`), and
   the probe saw 0 refinement pairs starting at a BranchPoint.
5. **A walk-through toggle (legacy).** `_branch_forward` / `_branch_backward`
   (`BaseManifold.py:540-620`) switch slots when a walk passes THROUGH a BranchPoint, from when crossings
   were BranchPoints (the class docstring `BranchPoint.py:30-36` still says so; only `FixedPoint`
   constructs one now, plus two tests). Measured: 0 walks pass through a BranchPoint on the six cases;
   every call enters it as a ROOT (`prev is None`). A through-walk on an unstable curve would return
   the STABLE branch's first node in the same slot (`branches_in` is the other stability's list): a
   latent trap that arrays should not reproduce.
6. **`first_node(b)`** (`BaseManifold.py:222-253`) steps off it onto slot `b`: growth views
   (`ManifoldMachine.py:217-227`), the initializer (`:395`), `grow_until_turnaround`
   (`TangleWorkbench.py:768-772`).
7. **Slot insertion.** `insert_point_forward/backward(node, b)` (`BranchPoint.py:116-168`) through
   `ManifoldMachine._insert_point_geometrically` :865-875 attach each chain piece's first point
   (`ManifoldInitializer.py:405-407`) and would splice after the anchor in merge / refine.
8. **The anchor of the synthetic crossing.** The registry anchor (`_register_anchors`,
   `TangleWorkbench.py:1585-1643`) reads the BranchPoint's coordinates as `u_nodes[0].get_point()`
   (:1632) and its handedness from `u_nodes[1] - u_nodes[0]` and `s_nodes[1] - s_nodes[0]`
   (:1666-1667: the fundamental segments leave the BranchPoint along the oriented eigenvectors).
   `Tangle.create_bridges` cuts every unstable branch whose `root is fp.branch_points[orbit]`
   (:894-898) at `(nodes[0], nodes[1])`: the anchor bridge is ROOTED AT THE BRANCHPOINT (one per
   unstable branch: 1 on k10, 3 on p3, 2 on inversion, 4 on inversion_period2), which `Bridge`
   allows (`_check_input_types` checks only the tail, `Bridge.py:188-206`).
   `Tangle._outer_neighbours` (:375-392) treats a BranchPoint end as having no neighbour; the separator
   after one goes into its slot (:1006-1007).
9. **Inversion: left of 0.0 = right of 0.1.** In numerics this is ONE object holding slot 0 and slot
   1 of each stability (one chain covers both, `ManifoldInitializer.construct_kevin_way` :375-411) and
   ONE anchor id for all four branches (`IntersectionRegistry._find_collision` at exactly (0, 0)).
   The reading "an arc of `u(o, b)` meets `s(o, b)`" is topology's (`Trellis.stable_key_at`,
   `Trellis.py:737`), keyed by the anchor id, not by the object.
10. **Public class.** Exported by `tanglepack` and `tanglepack.numerics` (`__init__.py:24,73`;
    `numerics/__init__.py:19,48`); built by hand in `tests/helpers/polylines.py:23` and
    `test_same_stability_diagnostic.py:52`; `tests/helpers/invariants.py:56-67` special-cases it.

Consequence for the split: roles 1, 3, 8 and 9 survive any storage (a store row for the periodic
point that sits in up to four curves' id arrays, plus the slot table); 4 and 5 can be retired with
evidence; 2 is an API naming question for the author; 7 becomes "set the curve in slot b".

---

## H. Oddities assigned to slices

Noticed by the analyst (findings 8), each handled by the slice whose file it lives in:

- the stale `merge_manifolds` comment (`ManifoldMachine.py:552-555`, "the newly generated point has
  no iterate links yet": the new image already holds `prev_iterate`): **N4**;
- `Tangle.populate_intersections_for_manifold` (:545-572) scanning the whole
  `_intersecting_segments` set per bridge, `_processed_pairs` only growing (believed quadratic over a
  long series; 0.2-0.6 s here): **N5**; measured in N5: 0.23 s over 32 blasts, crossing pairs only,
  not quadratic; left;
- `BranchPoint`'s dead fields, stale docstring and through-walk (`BaseManifold._branch_forward` /
  `_branch_backward` returning the other stability's node): **N8**;
- `grow_until_turnaround` reading `_coords` and rebinding `first_point` to a manifold
  (`TangleWorkbench.py:768-792`): **N6**;
- `stretch_param` per point although it is per chain, and `Point.insert_point_forward` giving a
  stretch-less node the host's (`Point.py:86-89`): **N4**.

Not in this campaign: the near-anchor depth stall (the zone filter, a separate request); why the
pre-`a8bb8a0` run peaked at 45 GB (not reproduced on `HEAD`); int32 ids (an N10 measurement if
wanted, Q2).

---

## Outcome

Filled slice by slice: the merge commit, the probe's verdict (identical, or the D3 differences
listed), and the census and blast time at the 1e-7 plateau against the baseline (census of
`figures/numerics_profile/a8bb8a0/scale_time.json`; 1913.8 MB traced at the plateau; blast 458.4 s).

| # | Slice | Slug | Merge | Probe (D3) | Census at the 1e-7 plateau | Blast s (1e-7) | What changed |
|---|---|---|---|---|---|---|---|
| 1 | N1 | `slim-segment-index` | | identical: the 58 builds and the two new dedup scenarios (`nested_unfiltered`, `shared_straddle`, each after the blasts and after a recompute), 62 dumps, registry ids included; under option (d) the scenarios differ in 15 items | 990.6 MB (segment index 1045.4 -> 122.0: 1.39 M query edge keys 118.4, 6,584 records); peak RSS 3034.7 -> 2344.8 MiB over the baseline (3152 -> 2462 absolute; base and branch measured the same day, the test engineer's paired run) | 399.6 against 399.0 on the base re-measured the same day (flat; the 458.4 baseline differs by run drift, `diagnose` included); `index` + `cut` 10.7 -> 9.4 | A query-only curve (a blast image curve, its bridges, a partial piece) keeps one int edge key per queried edge (`Tangle._query_edges`) and a `_Segment` only when its query paired; it owns no segment. Same dedup keys; a separator spliced into a query-only bridge takes no sid (later sids shift; sid values decide nothing). `_query_edges` never forgets on `reset=False`: a re-added image curve is not re-queried (Q1's accepted risk; the registry collapsed those re-detections before). An indexed add neither indexes nor queries an edge a query-only curve queried (as before N1; it no longer claims it, which had no effect), so the five `reset=False` scenarios register what the base did, up to the last bits after a `reset=False` pass (and in a `preserve_ids` recompute that follows one): resolution ratios <= 1.6e-14 relative, cdists <= 5.1e-16, and coordinates on nested. `populate_intersections_for_manifold(curve)` -> `populate_intersections(new_sids)`; `add_manifold` returns the new sids. Branch-cut bridges keep their claims. e02/e03 report `indexed_segments` (the stress report still reads e03's old `segments`: `report/analysis.py:386`, feeding `e03_stage_fits`, `make_figures.py` `fig_e03_cost` and `make_numbers.py:232-233`; repoint them before e03 is rerun). Open, older than N1: the accumulating `reset=False` pass misses real crossings a fresh pass finds (k28 two blasts grown one step: 16 against 30; nested: 48 against 54), and an image curve over newly grown stable manifold is not re-queried (Tangle Dev Notes). |
| 2 | N2 | `same-stability-count-new-only` | | identical: 62 dumps, 6024 arrays, 744 texts, the new `same_stability_counts` item included (65 blasts counted, 6 with hits: k28_series b08, k32_series b08-b10, nested_outer_series b08-b09); the item has teeth (dropping new-against-new pairs moves k32_series b08-b10); on the base, `diagnose=True` on the series and scenarios changes no other item | 990.6 MB, unchanged; peak RSS 2396.2 -> 2374.0 MiB over the baseline (2513 -> 2491 absolute, both at blast 13; base and branch run one after the other the same day; a concurrent pair read 2248.6 -> 2427.2, which does not reproduce). The count's transient at blast 13 is about 1.55 GB on both (traced), its 12.65 M candidate pairs being the same set; the narrow phase over them is unblocked, so it sets the peak | 109.2 against 415.8 on the base the same day (sequential); `diagnose` 360.3 -> 53.7 (15 %); the count fields of all 32 blasts equal, 23 with hits | `_candidate_pairs` with `new` pairs an old segment only with the new segments inside its forward range (a pool of positions with a `start`/`count` per sorted position) instead of building every pair and filtering after; the axis, `counts`, the orientation, the `_PAIR_BLOCK` chunking and `new=None` are the base's, so every hit is bitwise the same. No library output changed. The probe runs every blast with `diagnose=True` and dumps `same_stability_counts`; each profiler scale row keeps `same_stability_count`. Left: the rows, the edge dedup and the two sorts over every segment (about 1.1-1.5 s a blast at the plateau). |
| 3 | N3 | `point-store` | | identical: 62 dumps, 6024 arrays, 744 texts, `same_stability_counts` included. The merge ties are record-equal base against branch (k28_one_blast 17, heteroclinic 2, k32_series 2; the keeper is `head_1` in all 21, which holds a phantom in 15, 2 and 2 of them and a real link in all 21). Forcing `head_1_has_iter = False` moves 90 items (heteroclinic 2, k28_one_blast 5, k28_two_blasts 5, k28_series 44, k32_series 33); not counting phantom rows moves none, so the phantom term is pinned by a unit tie (`test_merge_tie_keeps_a_point_whose_only_link_is_a_phantom`, which passes on the base too and fails under that mutant) | 398.5 MB (990.6 on the base): a `Point` object 408.0 -> 100.0 B (blast image points 542.5 -> 133.0 MB, 1,394,339 of them); the 778,302 phantoms 285.0 MB of objects -> rows with no object; the store 104.0 MB for 2,177,788 rows (capacity 2,273,629 at 40 B, 86.7 MB, plus the row -> `Point` list, 17.2 MB: 50.1 B a row). Peak RSS 2338.8 -> 1883.6 MiB over the baseline (2456 -> 2001 absolute, both at blast 13, where the RSS held after the blast is 767.8 -> 356.8 MiB over the baseline; base then branch, one after the other the same day) | 109.7 against 106.7 on the base the same day (sequential; +2.8 %): `walk` 17.9 -> 18.9, `infer_iterates` 8.2 -> 10.3, `diagnose` 52.1 -> 53.1 (property reads per node), `refine` 8.7 -> 8.0 and `iterate` 0.86 -> 0.37 (batched rows); `same_stability_count`, `pip_depth` and `points_created` equal on all 32 blasts, 23 with hits | One append-only `PointStore` per workbench (`TangleWorkbench.point_store`, handed to the initializer): columns xy, cdist (NaN none), image and preimage rows (-1 none), growth x1.5, a row never removed or reused. `Point` is a slotted object over its row (`store, row, forward, backward, stretch_param`); `x`/`y` are read-only properties, `cdist` a property with a setter, the iterate links properties over the columns (a link to another store or to a non-`Point` raises `TypeError`), `get_point()` a read-only view; `edist` (no reader) and `_coords` (its `grow_until_turnaround` reads now use `get_point()`) are gone from `Point`. `BasePoint.__slots__ = ()` and `has_iterate_link()`, which the tie rule reads on both heads (a phantom row counts). The machine appends a step's images in one call (`ValueError` if an iterate exists, as before) and a refine layer its points then their phantom rows, so `_cache_preiterate` is gone and `_get_preiterate` returns coordinates; the machine appends to the store of the curve it grows, a separator to `p1.store`. A point built without a store goes into the process-wide `loose_store()`. A `Point` now lives as long as its store: the initializer's detached fictitious tail (one per seeded chain, a `Point` with no cdist) is kept, so `n_points` reads 2,177,788 against 2,177,786. No library output changed. The census sizes a slotted point as object plus row int and reports `point_store` with `store_rows`, `store_capacity`, `phantom_rows` and `off_curve_rows` (0 at the plateau). Left: `profile_figures.MEMORY_GROUPS` has no group for `point_store` (the component stack omits it); the `id()` sites stay on `id()` (N5/N6). |
| 4 | N4 | `curves-on-order-arrays` | | identical: 62 dumps, 6024 arrays, 744 texts, `same_stability_counts` included. Every live curve's rows through `_walk_nodes(stop_at_final=False)`, and its extent, are equal base against branch at all 1,006 operation boundaries of the 18 probe builds (the 11 law cases, the 5 series, the 2 scenarios); the chain invariants (no row twice or on two chains, the columns agree with every chain, cdist non-decreasing along it, every curve's root and tail on one chain in order) hold at those boundaries and at the 11,186 of the suite. The merge ties are record-equal (k28_one_blast 17, heteroclinic 2, k32_series 2; keeper `head_1` in all 21). The tie branches no case reaches, (a)-(e), are pinned by `test_merge_branches.py`, which passes on the base and the branch; a keeper that is always `head_1` fails (b), a root tie that drops the new root fails (c), the vectorised rule on ties fails (e) | 394.6 MB against 398.5 on the base (same-day pair): `points_image` 133.0 -> 101.1 MB (a `Point` 100.0 -> 76.0 B, slots `store, row`), `point_store` 104.0 -> 121.3 MB (the two int32 columns at capacity 2,273,629), `curve_ids` 10.8 MB (the chains' `ids`, 708 live chains); 298.6 -> 295.7 B a non-phantom row all in (1,399,486 of them); `phantom_rows` 778,302 on both and the same two phantom `Point` objects (the detached fictitious tails), so no phantom `Point` is re-created. Peak RSS 1890.5 -> 1860.9 MiB over the baseline (2008 -> 1979 absolute, both at blast 13; base then branch at `1c7cf16`, one after the other the same day, the test engineer's pair; the coder's pair, its branch run on an uncommitted tree, read 1882.8 -> 1812.4) | 85.0 against 112.2 on the base the same day (sequential, the test engineer's pair on committed trees; -24.2 %): `walk` 19.2 -> 1.3, `refine` 8.2 -> 1.0, `iterate` 0.38 -> 0.28, `index` 9.3 -> 8.9, `infer_iterates` 10.3 -> 10.3 (it reads brackets, not walks: N5), `diagnose` 54.5 -> 53.8 (the coder's pair, its branch run on an uncommitted tree: 113.1 -> 85.0); `same_stability_count`, `pip_depth` and `points_created` equal on all 32 blasts, 23 with hits. A first pair, with `Point` reading the public column views (a fresh view per read, ~330 ns against ~50), read 103.6 against 111.7 with `infer_iterates` 20.4 and `diagnose` 61.1; `Point` now reads the private columns | The order along a curve is a `Chain` of the store (`PointStore.chains`): the rows of one physical curve in curve order, with the row -> chain and row -> position lookups two int32 store columns (`_chain`, `_position`; -1 off every chain, and a `_position` of -2 for a row a merge left off), kept exact by every chain write (`Chain.insert`, one `np.insert`; `Chain.absorb`, which retires the merged chain and carries its anchor). A `BranchPoint` is never a row: the chain hanging off one names it as its `anchor`, the left neighbour of position 0. Every curve (branch, bridge, partial piece, blast image curve) is a view: `BaseManifold` keeps `root` / `tail` and finds its chain from the root's row (the root slot's point for a `BranchPoint` root) on every read (`_span`, `_rows`; `_collect` gathers by masks; `get_cdist_array` and `_walk_nodes(stop_at_final=False)` still read to the chain's end; `_find_tail` is the chain's last point, so growth still adopts blast images merged onto a branch). `merge_manifolds` is one `searchsorted` when the runs lie on different chains and share no cdist, and the old zipper on positions otherwise (branches (a)-(e) as before, the tie log and its `head_1` / `head_2` locals kept); the states one chain cannot hold raise `AssertionError` naming both cdist ranges: a run placed second that starts after its chain's beginning or hangs off a periodic point (g), a run that runs out with its chain going on past its tail (f), both runs ending on a tie (h). Refinement takes chain positions and inserts each layer with one `chain.insert` (the anchor is the left neighbour of position 0, as the walk read it; NaN left of the anchor, where the walk read the other stability's first point through the retired toggle, which can only lower that pair's area, through the library (no caller passes `refine_manifold` a `branch_index`) the base raised `TypeError` wherever such a pair was flagged (`test_order_arrays::test_refining_the_pair_at_the_periodic_point_lays_the_new_first_point` pins the branch and fails so on the base); with an explicit `branch_index` the base completed from the other stability's point and the two differ (on that construction with a stable branch added: 9 against 7 points at a cutoff between the two areas, 25 against 21 below both), a difference no build reaches); `refine_manifold` returns None. `_image_chains` reads chains and positions (an off-chain image starting a run is a chain of one, one inside a run is skipped, a periodic point first starts at its slot). Gone: `Point.forward` / `backward` / `insert_point_*` and the `stretch_param` slot (now the chain's, read-only, None off a chain; `Point(...)` takes no stretch), `BaseManifold.walk_fwd` / `walk_back` and the through-walk toggle `_branch_forward` / `_branch_backward` (role 5, retired here on the 0-through-walk evidence, Q6, not in N8 as campaign Q9 had it), the `ManifoldView` aliases, `_insert_point_geometrically`, the per-step same-store guard and the false merge comment. `BranchPoint.insert_point_*` set the slot and the chain's anchor, raise `ValueError` on another chain's point and lose `only_forward`. Public read-only `PointStore.xy` / `cdist` / `image` / `preimage` (a fresh view per access) and the writers `link`, `link_phantoms`, `set_cdist`, `new_chain`; `Point` reads the private columns and `get_point` a private read-only alias. No library output changed. The census walks chains, reports `curve_ids` and the two columns in `point_store`; `MEMORY_GROUPS` shows the store ("point store + phantoms") and the ids ("curve objects + ids + memo"). Left: `BranchPoint`'s module and class docstrings still describe the linked list (D5 keeps the class to N8); readers outside the engine still take nodes from `store.nodes` and the bridge memo stays (N6); brackets and `infer_iterates` (N5) |
| 5 | N5 | `segments-and-brackets-on-ids` | | | | | |
| 6 | N6 | `readers-on-views` | | | | | |
| 7 | N7 | `delete-point` | | | | | |
| 8 | N8 | `branch-point-storage` | | | | | |
| 9 | N9 | `array-namespace` | | | | | |
| 10 | N10 | `save-load-and-after-profile` | | | | | |
