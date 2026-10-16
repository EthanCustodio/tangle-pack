# Change: B1 `resolved-crossings` (coder, job `change`, iteration 1)

*Coder, 2026-10-12. Plan `artifacts/plan_planner_human_2.md`, approved as written (`plan_ok_planner_human.json`: the three questions take their defaults).*

## Headline: steps 0-7 are done, but the run has to stop at rubric 1

**Rubric 1 and 2 fail.** Four tests in files outside the plan's list pass on the base and fail on the branch. The plan's section 6 says "If one fails, the run stops and brings it to you", so I have not edited them. All four fail because of what the rule decides, not because of an implementation fault I could fix inside the 14 paths (details in section 5). They need an author decision.

**k28 x8 is the other thing the author must see.** The branch gives **4 holes / 4 seed / 4 closed** classes. The base gives 24 / 15 / 23, and the uncollapsed build gives 40 / 23 / 31. It is a dependence nobody predicted:
- The strong pip's deepest registered orbit point, f^8(q0), is a clause-1 detection: base id 249, ratio 2.991.
- Once it leaves N, the unstable pip-orbit cut (D10) falls back from u = 2.787e8 to f^7(q0) at u = 4.869e7.
- Base's blast-8 reference pairs (185, 186) and (186, 187), at u ≈ 2.36e8, then lie past the cut. The reference never switches, and the hole set stays at the blast-2..7 value of 4.

The analyst's "no base hole depends on an excluded id, by bounding id or image chain" is true but incomplete: the dependence runs through the pip-orbit cut. The probe confirms it: no base hole reads an excluded id by bounding id or by image chain.

**What was not affected:**
- **Law cases:** unchanged. The probe finds both clauses empty on all 11, |N| equals the base's registry count on each, and heteroclinic's minimum ratio is 3.038.
- **Recompute:** it now keeps N and the bridges on all five series.

## 1. Where the work is

| | |
|---|---|
| Base commit | `1029fcc` (author's `HEAD` on `main`; their dirty `docs/` and `.refactor/` files untouched) |
| Branch | `refactor/resolved-crossings` |
| Worktree | `/home/dezhu/Development/TanglePack/tangle-pack-resolved-crossings` |
| Base worktree | `/home/dezhu/Development/TanglePack/tangle-pack-resolved-crossings-base` (detached at `1029fcc`, clean) |

**Commits** (oldest first):
```
3557e1e Probe: the resolved-crossing oracle, both clauses by its own chord deviations (B1 step 0)
8115f51 Numerics: every detection is registered with its chord deviation; the registry decides N, and the blast vote is gone (B1 steps 1, 2, 4)
fe8d3db Readers: every path that takes the registry as the crossing set reads N (B1 step 3)
ebb7eaa Report: one WARNING names every unresolved detection with its clause and the remedy 'refine'; the threshold is part of the definition of N (B1 steps 5, 6)
41f19b5 Tests: the zig-zag, the symmetric pair and a clean crossing; every series step recomputes and the law compares N and the bridges (B1 step 7)
d1396c9 Tests: the zig-zag asserts the cut and the links before N, so the base fails on its cut
e201ddc Tangle: a blast's detections are handed over in unstable cdist order, so ids stay reproducible without the vote's sort (INT-7)
```

**`git diff --stat 1029fcc`** lists exactly the 14 plan paths:
```
 src/tanglepack/numerics/BridgeIterator.py          |  22 +-
 src/tanglepack/numerics/Intersection.py            |  10 +
 src/tanglepack/numerics/IntersectionRegistry.py    | 132 +++++++++-
 src/tanglepack/numerics/IterateInference.py        |  16 +-
 src/tanglepack/numerics/Tangle.py                  | 131 +++-------
 src/tanglepack/numerics/TangleWorkbench.py         | 139 +++++-----
 src/tanglepack/numerics/diagnostics.py             |  19 +-
 src/tanglepack/topology/StrongPip.py               |   4 +-
 src/tanglepack/topology/Trellis.py                 |  15 +-
 tests/_tools/resolution_probe.py                   | 290 +++++++++++++++++++++
 tests/cases.py                                     |  31 ++-
 .../regression/test_class_count_under_blasting.py  |  18 +-
 .../unit/numerics/test_noise_crossing_collapse.py  |  81 ------
 tests/unit/numerics/test_resolved_crossings.py     | 107 ++++++++
 14 files changed, 739 insertions(+), 276 deletions(-)
```

## 2. What changed, step by step

0. **Probe.** `tests/_tools/resolution_probe.py` takes the flags `--laws`, `--builders [..]`, `--series [..]` and `--k28x8`.
   - It wraps `Tangle._resolution_ratio` and records its own `max(h_u, h_s)` (`chord_deviation_batch` over `_outer_neighbours`). It never reads `Intersection.chord_deviation`, and it runs unchanged on base and branch.
   - It applies both clauses by itself. Clause 2 is a symmetric pair rule: consecutive in unstable-cdist order within `(u_key, s_key)`, on the same or a node-sharing original segment, separation ≤ the larger deviation.
   - It cross-checks against `registry.unresolved_crossings(3)` (`agrees=`) and against an arc-length "any distance along the curve" pass (`anycurve_extra=`).
   - It contains no mention of the deleted vote, so rubric 5's grep stays clean. The cost: it no longer counts collapse firings, so the uncollapsed 40 / 23 / 31 is quoted from the analyst, not re-run.

1. **Chord deviation.**
   - `Tangle._resolution_ratio` returns `(ratio, max(h_u, h_s))`. The deviation is unfloored; only the ratio's denominator keeps the floor.
   - `_resolve_crossing_pair` stores `Intersection.chord_deviation`, which is a new optional field on `Intersection` and `from_segments`.

2. **The rule, in the registry.** `IntersectionRegistry` now has:
   - `run_separation(id)`: as in the plan.
   - `run_partner(id)`: the nearest neighbour the detection forms a run with, None when it is in no run. It is not in the plan; I added it so the WARNING can name the partner without restating the rule.
   - `resolved_ids(threshold)`.
   - `unresolved_crossings(threshold)`: widened to the complement, in unstable order. `threshold=None` turns clause 1 off and keeps clause 2 (Q2 default).
   - The runs are recomputed over all registered detections, not per call. They are cached on `_orders_version` (inserts and renumberings), not on `generation`: an iterate link does not move N, and `infer_iterates` would otherwise rebuild the cache thousands of times.

3. **Readers of N.**
   - `TangleWorkbench`: `create_bridges` (the cut and the partial-piece count), `_frozen_ids("all")` and `trim_stable_manifolds`.
   - `BridgeIterator.iterate_bridge`: one `resolved_ids` call gives both the cut and the inference sources.
   - `infer_iterates`: one call gives both the sources and the candidate images (passed into `_bracketed_image`).
   - `Trellis.from_workbench`: the branch orderings, with the set recorded as `Trellis.resolved_ids`. A hand-built `Trellis` (the tests, the `MinimalTrellis` sparse snapshot) defaults to every registered id.
   - `StrongPip.is_strong_pip`: skips `r_id not in trellis.resolved_ids`.

4. **Vote deleted.** `_collapse_noise_crossings` and `_purge_crossing` are gone. `populate_intersections_for_manifold` returns every detection, **sorted by `(unstable_cdist, stable_cdist)`**. The vote's sort was what kept the blast path's registration order reproducible. Without it, ids followed the iteration order of a set of segment ids, and `test_registry_ids_are_reproducible` and `test_the_example_inversion_recipe_is_the_law_cases[1]` failed. Commit `e201ddc` fixes this, and both pass now.

5. **`report_resolution`.**
   - **The WARNING:** one per call, only when a detection not reported under the current numbering is unresolved. It names every unresolved detection: id, coordinates, ratio, clause, and for clause 2 its partner and nearest separation. It ends with "refine: lower `area_cutoff` or grow the stretch; NUM-3/4".
   - **`threshold=None`:** the special cases that silenced the report are gone, so with None the report covers clause 2 only.
   - **`diagnostics.py`:** two docstrings, plus the module paragraph "Both are FLAGS. Nothing here discards".

6. **Threshold docstring** (constructor `Args`, with a pointer from the class docstring). It covers:
   - "Part of the definition of N";
   - "Derived bound 1";
   - "a margin 3";
   - "the corpus minimum, over the law cases, is 3.038 (heteroclinic, 2026-10-12)";
   - the None behaviour;
   - the Q3 note: re-run `compute_intersections` after changing the threshold.

7. **Tests.**
   - `test_noise_crossing_collapse.py` is deleted.
   - `test_resolved_crossings.py` has three tests: the zig-zag, the symmetric pair and a clean crossing.
   - `cases.SeriesStep` gains `crossings` / `recomputed`, each `(N, bridge ids)`. Every `blast_series` step ends with `compute_intersections(workbench.fixed_points, preserve_ids=True, diagnose=False)` plus `rebuild_bridges()`.
   - N is read as "registry minus `unresolved_crossings(threshold)`", which also exists on the base, so rubric 10 can run there.
   - `test_recompute_keeps_the_crossing_set[name]` runs over `SERIES`.

## 3. Suite

`MPLBACKEND=Agg env/bin/python -m pytest -q -rxX` on the branch at `e201ddc`:
```
4 failed, 1118 passed, 38 skipped, 5 deselected, 1 xfailed in 190.04s (0:03:10)
XFAIL tests/invariants/test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing] - C3: INV-11 ...   (the one KNOWN_ISSUES entry)
FAILED tests/regression/test_coarse_blasts_rows_by_signs.py::test_coarse_blasts_rows_by_signs
FAILED tests/regression/test_evidence_diagnosis.py::test_an_order_flip_gives_no_verdict
FAILED tests/regression/test_open_p3_deep_runs.py::test_every_law_holds_on_the_deep_runs[inversion_period2_10]
FAILED tests/regression/test_pip_orbit_cut.py::test_p2inv10_every_class_walks
```
The first full run, before `e201ddc`, had two more failures: the id-reproducibility pair named in step 4, both fixed. The four above pass at `1029fcc`.

## 4. Behaviour changes

**Law cases (11): none.**
- Probe on the branch: clause 1 and clause 2 are empty on every case.
- |N| equals the base's registry count (k10 8, k28_one_blast 6, k28_two_blasts 10, p3 30, nested 40, nested_inner_blasts 44, inversion 13, inversion_one_blast 19, inversion_two_blasts 25, inversion_period2 19, heteroclinic 116).
- Heteroclinic's minimum ratio is 3.038, as on the base.
- The law tier's fingerprint guard passes; it is part of the 1118 passing tests.

**k28 series.**

| step | base | branch |
|---|---|---|
| blasts 0-5 | — | identical |
| blast 6 | vote keeps the run's median; registry 86, recompute → 88 | all three detections registered (registry 88); the run of three (u = 8356210.197) is out of N whole, including the would-be median (ratio 30.2) |
| blast 7 | a run of three at u = 47831087.57 | no run: the blast-6 run is no longer cut, so its image is not split at an old median (the probe at blast 7 shows clause 2 holding only the blast-6 run) |
| blast 8 | — | clause 1 holds out f^8(q0) (ratio 2.991) and the detection at ratio 2.572 (base ids 249 and 148) |

- At blast 6, nothing is cut at the run and no link touches it.
- Recompute keeps N and the bridges at every step (the law passes).

**k28 x8, re-derived and not pinned.**

| | holes / seed / closed | reliable |
|---|---|---|
| base | 24 / 15 / 23 | True |
| uncollapsed (analyst) | 40 / 23 / 31 | True |
| **branch** | **4 / 4 / 4** | True |

- **Unresolved set on the branch** (ids from this run; probe and library agree):
  - 75, 76, 77: clause 2, the blast-6 run. 75–76 are 2.05e-8 apart, 76–77 are 8.51e-10 apart, and every max h is 4.67e-6.
  - 148: clause 1, ratio 2.572.
  - 245: clause 1, ratio 2.991. It is f^8(q0).
- **Recompute:** registry 244 → 244, |N| 239 → 239 equal, bridges 168 → 168 equal.
- **On the base's recompute numbering**, the same rule gives 8: {1, 2, 75, 250, 251, 125} by clause 2 and {148, 249} by clause 1. This matches the analyst's correction.
- **Why 4 / 4 / 4:** the pip-orbit cut (headline; probe of the trellis cuts):
  - base unstable cut `278699284.69` = base id 249;
  - branch unstable cut `48689459.58` = id 147, i.e. f^7(q0);
  - base's new blast-8 references (185, 186) and (186, 187), at u ≈ 2.359e8 .. 2.375e8, sit between those two cuts.
- **The B8 collision** (two sources claiming image 213) cannot recur in that form, because its blast-7 run no longer exists. I did not check the IterateInference logger for it.

**Clause-1 changes elsewhere** (the probe on the branch):

| build | detections out of N | min ratio | effect |
|---|---|---|---|
| k32 series, blasts 9-10 | 1 | 1.48 | recompute equal |
| nested_outer series, blasts 9-11 | 1 | 1.78 | recompute equal |
| `k28_deep_blast_sliver` | 1 | 1.56 | — |
| `p3_coarse_blasts` | 4 | 0.648 | its under-resolved hole bridge is gone |
| `inversion_period2_10` | 2 (38, 39) | 1.03 / 2.33 | **the base's reference pair IS (38, 39)**; on the branch it switches to (30, 31), (31, 32); holes 13 → 10, seed 14 → 6, closed 17 → 6 |
| `near_tangency` | 152 of 261 | — | clause 2 empty |

**Every path:** an unresolved detection no longer disqualifies a strong pip, sets the stable trim, or becomes a target of `grow_until_iterates_closed`.

**WARNING text:** the unresolved-crossing WARNING has new content (it names every detection and gives the remedy). With `resolution_ratio_threshold=None` it can now fire, for clause 2 only.

## 5. The four failing tests (outside the list, not edited)

1. **`test_coarse_blasts_rows_by_signs`.** The precondition fails: "the builder no longer makes an under-resolved hole bridge". That hole bridge was cut at clause-1 detections (ratios 0.65-1.05), which are out of N now. The fixture no longer exercises the defect it was written for.
2. **`test_an_order_flip_gives_no_verdict`** (`build_k28(blasts=9)`). The precondition fails: no order flip remains. Same cause as k28 x8: f^8(q0) is unresolved, so T(q0) stops at f^7(q0) and the deep lambda-lemma copy that made the flip is past the cut.
3. **`test_every_law_holds_on_the_deep_runs[inversion_period2_10]`.** `helpers.laws.check_trellis_is_cut_at_the_pip_orbit` reads `trellis.registry` as the crossing set, so it expects 38 (ratio 1.03, inside both cuts) on a branch. The oracle needs to read N, i.e. skip `registry.unresolved_crossings(threshold)`. That is a one-line change in `tests/helpers/laws.py`, which is outside the list.
4. **`test_p2inv10_every_class_walks`.** On the branch, `check_every_active_class_walks` counts 0. The reference moved off the unresolved pair (38, 39), and the topology fell from 14 / 17 classes to 6 / 6. Five of the six are trellis-path classes and the one walked class does not count. `is_reliable` stays True.

The author has to choose for each: rewrite the precondition or the build (1, 2, 4), widen the oracle to N (3), or reconsider the rule's reach. In particular, should an unresolved pip-orbit point shorten T(q0)? That is what drives k28 x8, (2) and, indirectly, the k28 numbers.

## 6. Rubric

1. **FAIL.** Four failures outside the list (section 5). The only xfail is `cases.KNOWN_ISSUES`'s. None of the tests named in section 6 was edited, and they all pass: `test_crossing_resolution.py` including its `inf`-threshold blast, `test_resolution_warning_follows_oracle.py`, `test_deep_blast_sliver_beyond_pip.py`. Two that are not named in section 6 fail: `test_coarse_blasts_rows_by_signs.py`, and `test_open_p3_deep_runs.py` on `inversion_period2_10` only.
2. **FAIL**, for the same four tests: they pass on the base and fail on the branch. Every other test outside the list passes on both.
3. **PASS.** `git diff --stat 1029fcc..refactor/resolved-crossings` lists exactly the 14 paths.
4. **PASS for the law cases.** The fingerprint guard passes. The branch probe prints empty clause 1 and clause 2 on all 11, with |N| equal to the base's registry count. Heteroclinic's minimum is 3.038; nothing is below 3 on a law case. Outside the law cases, the k28 x8 and `inversion_period2_10` changes go beyond section 5's "can differ" in size; they are in section 4.
5. **PASS.**
   - `grep -rn "_collapse_noise_crossings\|_purge_crossing\|Collapsed" src tests` returns nothing.
   - `grep -n "% 2" src/tanglepack/numerics/Tangle.py` returns nothing.
   - `Tangle.py` went from 1115 to 1052 lines (63 shorter, ≥ 55).
6. **PASS in substance; the strict grep has more hits than listed.**
   - **Readers:** `TangleWorkbench` ×3 (`_frozen_ids` :1102, `create_bridges` :1803, `trim_stable_manifolds` :2178), `BridgeIterator` :378, `IterateInference` :109, `Trellis.from_workbench` :422, `StrongPip` :169.
   - **Trellis:** it also carries the recorded set as an attribute. That adds a docstring, a constructor parameter, the attribute, and `resolved_ids=resolved` at :488: one reader, several hits.
   - **IntersectionRegistry:** besides the definition, it has its module docstring, the interface line and `unresolved_crossings`'s own call.
   - Every one of the seven sites the findings name iterates N only.
7. **PASS.** `chord_deviation is None` exactly when `is_synthetic`, over every law case's registry (`scratchpad/b1run/rubric7.py`: 0 violations on all 11).
8. **PASS.**
   - `test_a_pair_is_judged_by_the_larger_of_its_two_deviations` passes.
   - Measured: h_short 4.96e-6 < separation 4.47e-3 < h_long 0.745, and the short detection's ratio is 403 ≥ 3.
   - Mutation check (`scratchpad/b1run/mutant.py`): rewriting `_runs` so each detection is compared only with its own deviation makes the test fail with an AssertionError.
9. **PASS.**
   - On the base worktree the zig-zag test fails on the cut: `assert [(1, 2), (2, 3), (3, 4), (4, 5)] == [(1, 5)]`, i.e. all three detections are cut at.
   - On the branch it passes: 3 registered, none in N, one bridge `(1, 5)` across them, no iterate entry touching them, the WARNING logged once and not again after a `preserve_ids` recompute.
10. **PASS.**
    - On the base (the branch's `cases.py` and test file copied in temporarily, then restored), `test_recompute_keeps_the_crossing_set` fails only for `k28_series`, first at blast 6: `[(6, '86 -> 88', '61 -> 63'), (7, '148 -> 150', '105 -> 106')]`. The other four series pass (1 failed, 4 passed).
    - On the branch all five pass.
11. **PASS.**
    - The `--k28x8` probe on the branch lists 75, 76, 77 (clause 2), 148 and 245 (clause 1), each with its reason. Each one satisfies its clause by the probe's own h.
    - `anycurve_extra=[]`: no detection outside the set satisfies either clause.
    - `agrees=True`: the probe's set equals `registry.unresolved_crossings(3)`.
    - Recompute: N 239 → 239 and bridges 168 → 168, both equal.
12. **PASS** (`scratchpad/b1run/rubric12.py`, a fresh process, k28 series blasts with `diagnose=True`).
    - Blasts 1-5 log no WARNING.
    - Blast 6 logs exactly 1 WARNING. It names 75, 76, 77 and contains "refine".
    - The following `compute_intersections(preserve_ids=True)` logs 0 WARNINGs.
    - Blast 8, for reference: one resolution WARNING (5 members). Its following recompute logs one same-stability WARNING (82 hits, the full listing a blast only counts), the documented C5 behaviour and the same on the base.
13. **PASS.** The `resolution_ratio_threshold` docstring states "definition of N", "Derived bound 1", "margin 3" and "the corpus minimum, over the law cases, is 3.038".

## 7. Probe output, base (`1029fcc`)

Run as `PYTHONPATH=src:tests MPLBACKEND=Agg env/bin/python <probe> --laws --builders --series --k28x8` from the base worktree root. Lines over 900 characters (`near_tangency`'s 152-entry clause-1 dict) are truncated.
```
k10: registry=8 detections=7 |N|=8 min_ratio=522.4 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
k28_one_blast: registry=6 detections=5 |N|=6 min_ratio=5783 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
k28_two_blasts: registry=10 detections=9 |N|=10 min_ratio=2455 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
p3: registry=30 detections=27 |N|=30 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
nested: registry=40 detections=36 |N|=40 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
nested_inner_blasts: registry=44 detections=40 |N|=44 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
inversion: registry=13 detections=12 |N|=13 min_ratio=124 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
inversion_one_blast: registry=19 detections=18 |N|=19 min_ratio=19.55 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
inversion_two_blasts: registry=25 detections=24 |N|=25 min_ratio=12.18 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
inversion_period2: registry=19 detections=17 |N|=19 min_ratio=19.8 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
heteroclinic: registry=116 detections=112 |N|=116 min_ratio=3.038 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (1s)
heteroclinic_17u: registry=180 detections=176 |N|=180 min_ratio=3.038 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (29s)
nested_unblasted: registry=34 detections=30 |N|=34 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
orientation_reversing: build raised ValueError: orbit point 0 has eigenvalues of disagreeing sign (unstable 4.23606797749979, stable -0.2360679774997897): the map is orientation reversing there, so only one of the two manifolds inverts and a single k_value cannot model it.
p3_15: registry=46 detections=43 |N|=46 min_ratio=3.038 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
p3_16: registry=68 detections=65 |N|=68 min_ratio=3.038 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
p3_6_blasts: registry=46 detections=43 |N|=46 min_ratio=3.038 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
p3_16_12: registry=206 detections=203 |N|=206 min_ratio=3.038 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
inversion_period2_10: registry=47 detections=45 |N|=45 min_ratio=1.031 clause1={38: '1.03', 39: '2.33'} clause2=[] unresolved=[38, 39] library_unresolved(3)=[38, 39] agrees=True anycurve_extra=[] missing_h=0 (0s)
p3_coarse_blasts: registry=48 detections=45 |N|=44 min_ratio=0.6477 clause1={40: '0.981', 41: '0.648', 42: '1.05', 46: '0.793'} clause2=[] unresolved=[40, 41, 42, 46] library_unresolved(3)=[40, 41, 42, 46] agrees=True anycurve_extra=[] missing_h=0 (0s)
near_tangency: registry=262 detections=261 |N|=110 min_ratio=0.1901 clause1={1: '1.74', 2: '0.963', 4: '2.67', 6: '2.47', 8: '1.65', 10: '1.47', 12: '2.3', 14: '1.46', 16: '1.61', 18: '2.58', 20: '2.73', 22: '0.881', 23: '0.42', 24: '0.775', 26: '0.19', 28: '1.04', 29: '2.77', 30: '1.01', 32: '0.201', 34: '0.788', 35: '0.425', 36: '0.822', 38: '2.69', 40: '2.56', 42: '1.37', 44: '1.22', 46: '2.03', 48: '1.22', 50: '1.36', 52: '2.69', 54: '2.78', 56: '0.818', 57: '0.225', 58: '0.738', 60: '0.217', 62: '1.04', 63: '2.93', 64: '0.991', 66: '0.231', 68: '0.774', 69: '0.236', 70: '0.85', 72: '2.67', 74: '2.56', 76: '1.25', 78: '1.11', 80: '1.91', 82: '1.12', 84: '1.26', 86: '2.77', 88: '2.84', 90: '0.787', 91: '0.221', 92: '0.714', 94: '0.219', 96: '1.07', 97: '2.89', 98: '0.964', 100: '0.243', 102: '0.788', 103: '0.242', 104: '0.862', 106: '2.61', 108: '2.49', 110: '1.19', 112: '1.07', 114:  ...[truncated: the full line is 3530 chars]
k28_deep_blast_sliver: registry=422 detections=421 |N|=421 min_ratio=1.556 clause1={423: '1.56'} clause2=[] unresolved=[423] library_unresolved(3)=[423] agrees=True anycurve_extra=[] missing_h=0 (10s)
k28_series blast 0: registry=4 detections=3 |N|=4 min_ratio=7929 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k28_series blast 1: registry=6 detections=5 |N|=6 min_ratio=5783 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k28_series blast 2: registry=10 detections=9 |N|=10 min_ratio=2455 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k28_series blast 3: registry=18 detections=17 |N|=18 min_ratio=2455 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k28_series blast 4: registry=30 detections=29 |N|=30 min_ratio=616.4 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k28_series blast 5: registry=50 detections=49 |N|=50 min_ratio=211.2 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k28_series blast 6: registry=86 detections=85 |N|=86 min_ratio=18 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k28_series blast 7: registry=146 detections=145 |N|=146 min_ratio=6.133 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k28_series blast 8: registry=248 detections=247 |N|=246 min_ratio=2.572 clause1={148: '2.57', 249: '2.99'} clause2=[] unresolved=[148, 249] library_unresolved(3)=[148, 249] agrees=True anycurve_extra=[] missing_h=0
  recompute: |N| 246 -> 250, bridges 175 -> 179, N equal False, bridges equal False
p3_series blast 0: registry=30 detections=27 |N|=30 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
p3_series blast 1: registry=32 detections=29 |N|=32 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
p3_series blast 2: registry=34 detections=31 |N|=34 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
p3_series blast 3: registry=36 detections=33 |N|=36 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
p3_series blast 4: registry=38 detections=35 |N|=38 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
p3_series blast 5: registry=40 detections=37 |N|=40 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
p3_series blast 6: registry=46 detections=43 |N|=46 min_ratio=3.038 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
p3_series blast 7: registry=52 detections=49 |N|=52 min_ratio=3.038 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
p3_series blast 8: registry=58 detections=55 |N|=58 min_ratio=3.038 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
p3_series blast 9: registry=67 detections=64 |N|=67 min_ratio=3.038 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
  recompute: |N| 67 -> 67, bridges 49 -> 49, N equal True, bridges equal True
inversion_series blast 0: registry=13 detections=12 |N|=13 min_ratio=124 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
inversion_series blast 1: registry=19 detections=18 |N|=19 min_ratio=19.55 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
inversion_series blast 2: registry=25 detections=24 |N|=25 min_ratio=12.18 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
inversion_series blast 3: registry=35 detections=34 |N|=35 min_ratio=12.18 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
inversion_series blast 4: registry=47 detections=46 |N|=47 min_ratio=12.18 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
  recompute: |N| 47 -> 47, bridges 33 -> 33, N equal True, bridges equal True
k32_series blast 0: registry=4 detections=3 |N|=4 min_ratio=8092 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k32_series blast 1: registry=6 detections=5 |N|=6 min_ratio=6537 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k32_series blast 2: registry=12 detections=11 |N|=12 min_ratio=761.4 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k32_series blast 3: registry=22 detections=21 |N|=22 min_ratio=761.4 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k32_series blast 4: registry=36 detections=35 |N|=36 min_ratio=761.4 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k32_series blast 5: registry=58 detections=57 |N|=58 min_ratio=93.31 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k32_series blast 6: registry=100 detections=99 |N|=100 min_ratio=27.51 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k32_series blast 7: registry=174 detections=173 |N|=174 min_ratio=12.99 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k32_series blast 8: registry=300 detections=299 |N|=300 min_ratio=3.24 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k32_series blast 9: registry=522 detections=521 |N|=521 min_ratio=1.477 clause1={523: '1.48'} clause2=[] unresolved=[523] library_unresolved(3)=[523] agrees=True anycurve_extra=[] missing_h=0
k32_series blast 10: registry=898 detections=897 |N|=897 min_ratio=1.477 clause1={523: '1.48'} clause2=[] unresolved=[523] library_unresolved(3)=[523] agrees=True anycurve_extra=[] missing_h=0
  recompute: |N| 897 -> 897, bridges 639 -> 639, N equal True, bridges equal True
nested_outer_series blast 0: registry=34 detections=30 |N|=34 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 1: registry=36 detections=32 |N|=36 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 2: registry=40 detections=36 |N|=40 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 3: registry=48 detections=44 |N|=48 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 4: registry=60 detections=56 |N|=60 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 5: registry=80 detections=76 |N|=80 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 6: registry=116 detections=112 |N|=116 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 7: registry=176 detections=172 |N|=176 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 8: registry=276 detections=272 |N|=276 min_ratio=3.094 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 9: registry=448 detections=444 |N|=447 min_ratio=1.775 clause1={447: '1.78'} clause2=[] unresolved=[447] library_unresolved(3)=[447] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 10: registry=736 detections=732 |N|=735 min_ratio=1.775 clause1={447: '1.78'} clause2=[] unresolved=[447] library_unresolved(3)=[447] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 11: registry=1218 detections=1214 |N|=1217 min_ratio=1.775 clause1={447: '1.78'} clause2=[] unresolved=[447] library_unresolved(3)=[447] agrees=True anycurve_extra=[] missing_h=0
  recompute: |N| 1217 -> 1217, bridges 864 -> 864, N equal True, bridges equal True
k28 x8: holes=24 seed=15 closed=23 reliable=True
k28 x8: registry=248 detections=247 |N|=246 min_ratio=2.572 clause1={148: '2.57', 249: '2.99'} clause2=[] unresolved=[148, 249] library_unresolved(3)=[148, 249] agrees=True anycurve_extra=[] missing_h=0
k28 x8 recomputed: registry=252 detections=251 |N|=244 min_ratio=2.572 clause1={148: '2.57', 249: '2.99'} clause2=['(1,75 sep 2e-08 h 4.7e-06/4.7e-06)', '(75,2 sep 8.5e-10 h 4.7e-06/4.7e-06)', '(250,125 sep 6.1e-08 h 8.4e-06/2.9e-06)', '(125,251 sep 1.8e-07 h 2.9e-06/8.3e-06)'] unresolved=[1, 2, 75, 125, 148, 249, 250, 251] library_unresolved(3)=[148, 249] agrees=False anycurve_extra=[] missing_h=0
  recompute: registry 248 -> 252, |N| 246 -> 250, N equal False, bridges 175 -> 179 equal False
    id 1: u=8356210.1974 ratio=2297 coords=(1.92161081593, 1.87116945892) clause 2 (partner 75, sep 2.05e-08 <= max h 4.67e-06)
    id 2: u=8356210.19741 ratio=4642 coords=(1.92161083642, 1.87116946492) clause 2 (partner 75, sep 8.51e-10 <= max h 4.67e-06)
    id 75: u=8356210.19741 ratio=30.21 coords=(1.9216108356, 1.87116946468) clause 2 (partner 1, sep 2.05e-08 <= max h 4.67e-06)
    id 125: u=47831087.5662 ratio=12.03 coords=(2.76375752944, -1.92161006737) clause 2 (partner 250, sep 6.13e-08 <= max h 8.43e-06)
    id 148: u=226301907.053 ratio=2.572 coords=(2.94935770635, -2.94935221417) clause 1 (ratio 2.572)
    id 249: u=278699284.692 ratio=2.991 coords=(2.94935855712, -2.94935708394) clause 1 (ratio 2.991)
    id 250: u=47831087.5662 ratio=623.8 coords=(2.76375746982, -1.92161005297) clause 2 (partner 125, sep 6.13e-08 <= max h 8.43e-06)
    id 251: u=47831087.5663 ratio=632.6 coords=(2.76375770119, -1.92161010883) clause 2 (partner 125, sep 1.77e-07 <= max h 8.33e-06)
  (holes reading no excluded id are not listed)
done
```

## 8. Probe output, branch (`e201ddc`)

Same command from the branch worktree root.
```
k10: registry=8 detections=7 |N|=8 min_ratio=522.4 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
k28_one_blast: registry=6 detections=5 |N|=6 min_ratio=5783 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
k28_two_blasts: registry=10 detections=9 |N|=10 min_ratio=2455 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
p3: registry=30 detections=27 |N|=30 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
nested: registry=40 detections=36 |N|=40 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
nested_inner_blasts: registry=44 detections=40 |N|=44 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
inversion: registry=13 detections=12 |N|=13 min_ratio=124 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
inversion_one_blast: registry=19 detections=18 |N|=19 min_ratio=19.55 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
inversion_two_blasts: registry=25 detections=24 |N|=25 min_ratio=12.18 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
inversion_period2: registry=19 detections=17 |N|=19 min_ratio=19.8 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
heteroclinic: registry=116 detections=112 |N|=116 min_ratio=3.038 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (1s)
heteroclinic_17u: registry=180 detections=176 |N|=180 min_ratio=3.038 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (29s)
nested_unblasted: registry=34 detections=30 |N|=34 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
orientation_reversing: build raised ValueError: orbit point 0 has eigenvalues of disagreeing sign (unstable 4.23606797749979, stable -0.2360679774997897): the map is orientation reversing there, so only one of the two manifolds inverts and a single k_value cannot model it.
p3_15: registry=46 detections=43 |N|=46 min_ratio=3.038 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
p3_16: registry=68 detections=65 |N|=68 min_ratio=3.038 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
p3_6_blasts: registry=46 detections=43 |N|=46 min_ratio=3.038 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
p3_16_12: registry=206 detections=203 |N|=206 min_ratio=3.038 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0 (0s)
inversion_period2_10: registry=47 detections=45 |N|=45 min_ratio=1.031 clause1={38: '1.03', 39: '2.33'} clause2=[] unresolved=[38, 39] library_unresolved(3)=[38, 39] agrees=True anycurve_extra=[] missing_h=0 (0s)
p3_coarse_blasts: registry=48 detections=45 |N|=44 min_ratio=0.6477 clause1={40: '0.981', 41: '0.648', 42: '1.05', 46: '0.793'} clause2=[] unresolved=[40, 41, 42, 46] library_unresolved(3)=[40, 41, 42, 46] agrees=True anycurve_extra=[] missing_h=0 (0s)
near_tangency: registry=262 detections=261 |N|=110 min_ratio=0.1901 clause1={1: '1.74', 2: '0.963', 4: '2.67', 6: '2.47', 8: '1.65', 10: '1.47', 12: '2.3', 14: '1.46', 16: '1.61', 18: '2.58', 20: '2.73', 22: '0.881', 23: '0.42', 24: '0.775', 26: '0.19', 28: '1.04', 29: '2.77', 30: '1.01', 32: '0.201', 34: '0.788', 35: '0.425', 36: '0.822', 38: '2.69', 40: '2.56', 42: '1.37', 44: '1.22', 46: '2.03', 48: '1.22', 50: '1.36', 52: '2.69', 54: '2.78', 56: '0.818', 57: '0.225', 58: '0.738', 60: '0.217', 62: '1.04', 63: '2.93', 64: '0.991', 66: '0.231', 68: '0.774', 69: '0.236', 70: '0.85', 72: '2.67', 74: '2.56', 76: '1.25', 78: '1.11', 80: '1.91', 82: '1.12', 84: '1.26', 86: '2.77', 88: '2.84', 90: '0.787', 91: '0.221', 92: '0.714', 94: '0.219', 96: '1.07', 97: '2.89', 98: '0.964', 100: '0.243', 102: '0.788', 103: '0.242', 104: '0.862', 106: '2.61', 108: '2.49', 110: '1.19', 112: '1.07', 114:  ...[truncated: the full line is 3530 chars]
k28_deep_blast_sliver: registry=422 detections=421 |N|=421 min_ratio=1.556 clause1={423: '1.56'} clause2=[] unresolved=[423] library_unresolved(3)=[423] agrees=True anycurve_extra=[] missing_h=0 (10s)
k28_series blast 0: registry=4 detections=3 |N|=4 min_ratio=7929 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k28_series blast 1: registry=6 detections=5 |N|=6 min_ratio=5783 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k28_series blast 2: registry=10 detections=9 |N|=10 min_ratio=2455 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k28_series blast 3: registry=18 detections=17 |N|=18 min_ratio=2455 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k28_series blast 4: registry=30 detections=29 |N|=30 min_ratio=616.4 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k28_series blast 5: registry=50 detections=49 |N|=50 min_ratio=211.2 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k28_series blast 6: registry=88 detections=87 |N|=85 min_ratio=18 clause1={} clause2=['(75,76 sep 2e-08 h 4.7e-06/4.7e-06)', '(76,77 sep 8.5e-10 h 4.7e-06/4.7e-06)'] unresolved=[75, 76, 77] library_unresolved(3)=[75, 76, 77] agrees=True anycurve_extra=[] missing_h=0
k28_series blast 7: registry=146 detections=145 |N|=143 min_ratio=6.133 clause1={} clause2=['(75,76 sep 2e-08 h 4.7e-06/4.7e-06)', '(76,77 sep 8.5e-10 h 4.7e-06/4.7e-06)'] unresolved=[75, 76, 77] library_unresolved(3)=[75, 76, 77] agrees=True anycurve_extra=[] missing_h=0
k28_series blast 8: registry=244 detections=243 |N|=239 min_ratio=2.572 clause1={148: '2.57', 245: '2.99'} clause2=['(75,76 sep 2e-08 h 4.7e-06/4.7e-06)', '(76,77 sep 8.5e-10 h 4.7e-06/4.7e-06)'] unresolved=[75, 76, 77, 148, 245] library_unresolved(3)=[75, 76, 77, 148, 245] agrees=True anycurve_extra=[] missing_h=0
  recompute: |N| 239 -> 239, bridges 168 -> 168, N equal True, bridges equal True
p3_series blast 0: registry=30 detections=27 |N|=30 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
p3_series blast 1: registry=32 detections=29 |N|=32 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
p3_series blast 2: registry=34 detections=31 |N|=34 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
p3_series blast 3: registry=36 detections=33 |N|=36 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
p3_series blast 4: registry=38 detections=35 |N|=38 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
p3_series blast 5: registry=40 detections=37 |N|=40 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
p3_series blast 6: registry=46 detections=43 |N|=46 min_ratio=3.038 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
p3_series blast 7: registry=52 detections=49 |N|=52 min_ratio=3.038 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
p3_series blast 8: registry=58 detections=55 |N|=58 min_ratio=3.038 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
p3_series blast 9: registry=67 detections=64 |N|=67 min_ratio=3.038 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
  recompute: |N| 67 -> 67, bridges 49 -> 49, N equal True, bridges equal True
inversion_series blast 0: registry=13 detections=12 |N|=13 min_ratio=124 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
inversion_series blast 1: registry=19 detections=18 |N|=19 min_ratio=19.55 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
inversion_series blast 2: registry=25 detections=24 |N|=25 min_ratio=12.18 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
inversion_series blast 3: registry=35 detections=34 |N|=35 min_ratio=12.18 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
inversion_series blast 4: registry=47 detections=46 |N|=47 min_ratio=12.18 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
  recompute: |N| 47 -> 47, bridges 33 -> 33, N equal True, bridges equal True
k32_series blast 0: registry=4 detections=3 |N|=4 min_ratio=8092 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k32_series blast 1: registry=6 detections=5 |N|=6 min_ratio=6537 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k32_series blast 2: registry=12 detections=11 |N|=12 min_ratio=761.4 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k32_series blast 3: registry=22 detections=21 |N|=22 min_ratio=761.4 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k32_series blast 4: registry=36 detections=35 |N|=36 min_ratio=761.4 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k32_series blast 5: registry=58 detections=57 |N|=58 min_ratio=93.31 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k32_series blast 6: registry=100 detections=99 |N|=100 min_ratio=27.51 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k32_series blast 7: registry=174 detections=173 |N|=174 min_ratio=12.99 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k32_series blast 8: registry=300 detections=299 |N|=300 min_ratio=3.24 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
k32_series blast 9: registry=522 detections=521 |N|=521 min_ratio=1.477 clause1={523: '1.48'} clause2=[] unresolved=[523] library_unresolved(3)=[523] agrees=True anycurve_extra=[] missing_h=0
k32_series blast 10: registry=898 detections=897 |N|=897 min_ratio=1.477 clause1={523: '1.48'} clause2=[] unresolved=[523] library_unresolved(3)=[523] agrees=True anycurve_extra=[] missing_h=0
  recompute: |N| 897 -> 897, bridges 638 -> 638, N equal True, bridges equal True
nested_outer_series blast 0: registry=34 detections=30 |N|=34 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 1: registry=36 detections=32 |N|=36 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 2: registry=40 detections=36 |N|=40 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 3: registry=48 detections=44 |N|=48 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 4: registry=60 detections=56 |N|=60 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 5: registry=80 detections=76 |N|=80 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 6: registry=116 detections=112 |N|=116 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 7: registry=176 detections=172 |N|=176 min_ratio=8.704 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 8: registry=276 detections=272 |N|=276 min_ratio=3.094 clause1={} clause2=[] unresolved=[] library_unresolved(3)=[] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 9: registry=448 detections=444 |N|=447 min_ratio=1.775 clause1={447: '1.78'} clause2=[] unresolved=[447] library_unresolved(3)=[447] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 10: registry=736 detections=732 |N|=735 min_ratio=1.775 clause1={447: '1.78'} clause2=[] unresolved=[447] library_unresolved(3)=[447] agrees=True anycurve_extra=[] missing_h=0
nested_outer_series blast 11: registry=1218 detections=1214 |N|=1217 min_ratio=1.775 clause1={447: '1.78'} clause2=[] unresolved=[447] library_unresolved(3)=[447] agrees=True anycurve_extra=[] missing_h=0
  recompute: |N| 1217 -> 1217, bridges 863 -> 863, N equal True, bridges equal True
k28 x8: holes=4 seed=4 closed=4 reliable=True
k28 x8: registry=244 detections=243 |N|=239 min_ratio=2.572 clause1={148: '2.57', 245: '2.99'} clause2=['(75,76 sep 2e-08 h 4.7e-06/4.7e-06)', '(76,77 sep 8.5e-10 h 4.7e-06/4.7e-06)'] unresolved=[75, 76, 77, 148, 245] library_unresolved(3)=[75, 76, 77, 148, 245] agrees=True anycurve_extra=[] missing_h=0
k28 x8 recomputed: registry=244 detections=243 |N|=239 min_ratio=2.572 clause1={148: '2.57', 245: '2.99'} clause2=['(75,76 sep 2e-08 h 4.7e-06/4.7e-06)', '(76,77 sep 8.5e-10 h 4.7e-06/4.7e-06)'] unresolved=[75, 76, 77, 148, 245] library_unresolved(3)=[75, 76, 77, 148, 245] agrees=True anycurve_extra=[] missing_h=0
  recompute: registry 244 -> 244, |N| 239 -> 239, N equal True, bridges 168 -> 168 equal True
    id 75: u=8356210.1974 ratio=2297 coords=(1.92161081593, 1.87116945892) clause 2 (partner 76, sep 2.05e-08 <= max h 4.67e-06)
    id 76: u=8356210.19741 ratio=30.21 coords=(1.9216108356, 1.87116946468) clause 2 (partner 75, sep 2.05e-08 <= max h 4.67e-06)
    id 77: u=8356210.19741 ratio=4642 coords=(1.92161083642, 1.87116946492) clause 2 (partner 76, sep 8.51e-10 <= max h 4.67e-06)
    id 148: u=226301907.053 ratio=2.572 coords=(2.94935770635, -2.94935221417) clause 1 (ratio 2.572)
    id 245: u=278699284.692 ratio=2.991 coords=(2.94935855712, -2.94935708394) clause 1 (ratio 2.991)
  (holes reading no excluded id are not listed)
done
```

## 9. Scratch scripts

All are in `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b1run/`:

| script | checks |
|---|---|
| `rubric7.py` | rubric 7 |
| `mutant.py` | rubric 8 |
| `rubric12.py` | rubric 12 |
| `witness.py` | k28 references per blast, base and branch |
| `refs8.py` | the k28 x8 pip-orbit cuts |
| `p2inv.py`, `p2inv_b.py` | `inversion_period2_10` before and after |
| `base_rubric10.out` | the rubric-10 run on the base |
