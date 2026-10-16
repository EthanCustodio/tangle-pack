# Findings — B2 `pseudoneighbor-both-intervals` (analyst, iteration 1)

*2026-10-14. (The analyst's Write tool refused the artifact; the orchestrator saved this text verbatim.)*

**The working tree is dirty at `HEAD` `ba23a1f` (main).** The uncommitted files are documentation and run archives only: `.refactor/decisions.md`, `.refactor/knowledge/{loom,numerics,testing,topology}.md` (this run updated `topology.md` and `testing.md`), `docs/alignment_plan_2026-10-09.md`, `docs/code_fidelity_audit_2026-10-09.md`, the moved `artifacts/request_origin_analyst.md`, and the untracked `.refactor/runs/2026-10-1{2,3,4}-*/` and `artifacts/`. No file under `src/`, `tests/` or `stress/` is modified. The B2 branch starts from `HEAD` and should not touch any of these.

**What kind of request.** This is a **definition-completion refactor** of one function with one caller, in the topology layer: `src/tanglepack/topology/Pseudoneighbor.py` (`_unstable_interval_is_clear`, `compute_pseudoneighbors` and the Dev Notes). It adds a unit test file and a read-only stress script. Nothing is planned in the numerics or loom layers. I did not re-plan the slice. The sections below verify what it assumes and list where it differs from the code or the measurements.

## 1. Short answer

**No. The code tests only the open unstable interval against X.** It has no stable-side check: the "stable check" is just that the walk steps through *consecutive* crossings of the cut trellis on the pip's branch. The plan's S-a function would close that gap. **On every case and series step I measured, it changes no verdict.** The B0 law passes, my independent mirrored stable reading finds **0 punctures**, and this holds on 13 builds and 73 series steps. I also simulated a missing window representative, and even then the new clause never decides a pair that the unstable clause has not already decided. So B2 is purely definitional on this corpus. e13 should report empty differences, which means B2 gives DEC-25 no new input.

There are five plan details to fix before coding (section 5): the m-range when `log_lap < 0`, which β to use, the candidate set, the third unit test's "fails on base", and e13's depths and output file.

## 2. The code at `HEAD` (verified by reading `ba23a1f`)

| What | Where (`HEAD`) | Notes |
|---|---|---|
| Unstable reading | `Pseudoneighbor.py:329-471` (the plan says 329-470: the same function, plus one line) | Gate `:363-365` requires the same unstable key, otherwise False. `owner = u_key[0]`, `branch_position_map("unstable")`, `per_step_beta("unstable")` `:367-370`. The `lo <= 0` rejection is at `:373-374`. Closed-form `m` range `:380-401`. The three readings in order: **table landing by id** `:410-429`, **own frame** against `iterate(x, -steps)` `:434-458`, **scaled** with the 2-D guard `:460-470`. The guard applies only when `same_stable` (`:394`, D7). |
| Stable "check" | **none in code**. Consecutiveness comes from `trellis.next_intersection(..., "stable")` at `:209` and `:224` | The branch lists hold only the cut N: `Trellis.from_workbench` buckets `resolved_ids` only (`Trellis.py:433-446`, B1), and `cut_at_pip_orbit` (`:264`) removes crossings past the cut. The open stable interval of a consecutive pair therefore holds no crossing of the cut N, and nothing more is checked. |
| Caller | `compute_pseudoneighbors` `:149-227`. The single call is at `:217`, with candidates `trellis.own_intersection_ids` (`:186`) | The docstring at `:159-160` says "keep each pair whose open unstable interval contains no…". It must change in step 3 as well. |
| Dev Notes to edit (step 3) | "empty by construction; all the work is in the unstable-interval check" `:64-66` | **There is no "NUM-18" paragraph in `Pseudoneighbor.py`** (`grep NUM-18` finds nothing). NUM-18 is in `.refactor/knowledge/topology.md:215` and the audit (`docs/code_fidelity_audit_2026-10-09.md:155-169`). The nearest code paragraphs are the orphan-pair note `:50-58` (still true; it can stay), the unstable-only closed form `:83-89` and the unstable-only check-set note `:117-126`. |
| Sibling shape | `StrongPip.py:149-211` | Uses `branch_position_map("stable")`, table first (`iterate(r, -steps)`), and otherwise `image_cdist` on both sides with `same_point`. |
| 2-D guard | `Trellis.same_point` `Trellis.py:109-126` | |
| β | `FixedPoint.per_step_beta` `FixedPoint.py:252` | `"stable"` returns `1/per_step_beta("unstable")` when `λ_u·λ_s = 1` within 1e-9 (`:280-290`). The `Trellis.image_cdist` docstring adds that a stable distance is scaled by `per_step_beta("stable")`, **never** by the reciprocal of the unstable factor. |
| Other users of the name | none outside `Pseudoneighbor.py` | Tests do not import `_unstable_interval_is_clear` (`grep -rn`), so rubric 4's grep only constrains `src`. Docs, knowledge files and `docs/open_issues.md:84,473` still mention the old name. |

## 3. Measurements (all verified by running code at `HEAD`)

Probe scripts (session scratchpad, read-only): `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/{probe,missing2,unres}.py`.

`probe.py` implements the plan's stable reading, mirrored line for line from `:360-471`: the stable key asserted equal, `branch_position_map("stable")`, β = `per_step_beta("stable")`, candidates = `own_intersection_ids` gated by `pos_map.get(r.manifold_b_key)`, the 2-D guard only when the unstable fixed point is the same, and the three readings in order. It runs on **every consecutive window pair**, not just on the references. I checked that it does fire: on a pair that skips a window crossing it fires by table (k10, k28×2, p3, inversion).

**3a. The B0 law.** `env/bin/python -m pytest tests/invariants/test_law_partition.py`: 11 passed, 4 skipped (the skips are NOT_APPLICABLE). I also called `laws.check_reference_pairs_clear_of_X(case)` directly; it passes on all 11 law cases and on `p3_6_blasts` and `k28` x8 (count = 2 × the number of references).

**3b. Stable punctures from the mirrored reading = 0 on every build:**

| Build | window pairs | references | stable punctures |
|---|---|---|---|
| k10, k28_one_blast, k28_two_blasts, p3, nested, nested_inner_blasts | 4, 2, 4, 2, 6, 6 | 1, 2, 1, 2, 3, 3 | 0 |
| inversion, inversion_one_blast, inversion_two_blasts, inversion_period2, heteroclinic | 2, 4, 6, 2, 20 | 2, 1, 2, 2, 4 | 0 |
| p3_6_blasts, k28 x8 (`build_k28(blasts=8)`) | 6, 59 | 2, 1 | 0 |

Series: each `cases.SERIES` build blasted with `cases._blast` at its own `min_separation` and the recipe cutoff (1e-7; inversion 1e-5), **to the E12 depths**:

| Series | depth | stable punctures | references per step |
|---|---|---|---|
| `k28_series` | 0..13 | 0 at every step | 2,2,1,1,… |
| `p3_series` | 0..20 | 0 | 2 throughout |
| `inversion_series` | 0..10 | 0 | 2,1,2,2,2,2,**1**,1,1,1,1 (the blast-6 drop is DEC-25; the stable clause plays no part in it) |
| `k32_series` | 0..13 | 0 | 2…2, 6 from blast 8 |
| `nested_outer_series` | 0..12 | 0 | 4,4,3,… |

All 73 steps record `"stable_fires": []`. The trellis crossing count stops growing (k28 N = 143 from blast 7, k32 300 from 8, nested_outer 276 from 8, p3 119 from 14, inversion 167 from 9). That is the pip-orbit stall that B3 addresses, so past it e13 measures the same N again.

**3c. Where a stable rejection *could* fire, and whether it does.**
- **Unresolved detections (B1).** On `k28_series` from blast 6 on, a run of unresolved detections (75, 76, 77; later also 606-608) sits inside the open stable interval of the window pair (64, 14) or (64, 102). That pair is not a reference. These detections have unstable cdist 8.4e6, far past the unstable cut, and no table link (`unres.py`). Since B1, iterate inference links only crossings of N (`numerics/IterateInference.py:91,117,122`), so an unresolved detection is never a table landing. **Result: 0 fires.**
- **A missing window representative (simulated, `missing2.py`).** I deleted each interior window crossing `w` from the candidates **and** from the table (no link into or out of it, which is what B1 does to an unresolved detection), then judged the merged pair `(prev, next)`. The stable reading catches `w` through the **own-frame** reading in about half the cases (for example k28 x8: 34 of 58). **In every one of those cases the unstable reading already rejects the merged pair.** On these 13 builds the only merged pair the unstable reading keeps is one on `p3_6_blasts`, and the stable reading does not catch that one either: `w`'s orbit has no other point in T. So A.5's claim that the clause adds robustness to a missing representative holds in principle but **changes no verdict on the corpus**.
- **Dedup.** I found no case. Inside a window, a dedup would show up as a missing representative, which is covered by the previous point.

**3d. A trap for the candidate set (verified).** Many window pairs, **references included**, have resolved registered crossings strictly inside their open stable interval that lie **past the unstable pip-orbit cut** (not in T):
- the law case `inversion`: reference (2, 3) holds 23 and 24 (u = 28.3 against an unstable cut of 25.9; no table orbit from −6 to 6);
- `inversion_series`@4-5: reference (27, 26);
- `k32_series`@13: all 6 references;
- `nested_outer_series`@11-12: reference (38, 37);
- 40 to 981 such points per step on the deep k28, k32 and nested series.

These points are not in `X = ∪ M^n(N_T)`. Each one is its own orbit's single window visit and lies outside T, so no orbit point is in T. The code and the B0 law both use `own_intersection_ids` (the cut N), so they agree. **If the stable reading took its candidates from `registry` or `resolved_ids` instead, the law case `inversion` would lose reference (2, 3), and fingerprints would change** (rubric 2).

**3e. KNOWN_ISSUES.** `tests/cases.py:1321-1330` holds exactly one entry, `("orientation_reversing_case_builds", "orientation_reversing")`. **No entry exists for `reference_pairs_clear_of_X`**, so step 1's `tests/cases.py` edit is a no-op. NOT_APPLICABLE (`:1341-1390`) has none for it either.

## 4. `e12_reverify.py`, the model for e13 (read)

- **Shape.** A harness module (`stress/numerics/harness.py`): `EXPERIMENT` (`:56`), `configs(quick)` (`:70`), `run(config, rec)` (`:98`), `harness.main(globals())` (`:198`). The CLI is `--quick --workers --timeout --rss-cap-gb --out --only --list`. There is one subprocess per config. Output is one JSON line per config in `stress/numerics/results/<EXPERIMENT>[_quick].jsonl`, and an existing file is renamed to `.bak<time>`. Values are recorded per blast with `rec.append`, with a `rec.checkpoint()` after each blast. `TIMEOUT=900`, `WORKERS=8`, `RSS_CAP_GB=12`, `MAX_POINTS=4e6`, and a wall/memory budget stop (`stop_reason`).
- **Cases and depths.** `DEPTH` (`:65`) is k28 13, k32 13, p3 20, nested_outer 12, nested_inner 12, inversion 10, heteroclinic 8. `SEPARATION` (`:63`) matches each `cases.SERIES.min_separation`. Every series runs at **two cutoffs** (`session.area_cutoff` = 1e-7 and 1e-5, `:83-86`), with one iteration per call, plus `TWO_ITERATIONS` cells and a `depth` study. The builds come from `e11_blast_depth.build(case)` (`e11:224-254`): the same `cases` recipes as `SERIES`, plus `nested_inner` and `heteroclinic`, which are **not** in `SERIES`.
- **Report layout.** The `e11.record_topology` schema (`n_classes` seed, `n_closed`, `n_holes`, `n_trellis_crossings`, …). E12's report is built in `stress/numerics/report/{analysis,make_numbers,make_figures}.py` + `report.tex` (`make_numbers.py:579 e12()`), so wiring e13 into the PDF would touch more files than the slice allows.
- **The results file is committed.** `stress/numerics/results/e12_reverify.jsonl` was added in `a97fcb6`.

## 5. Where the plan section differs from the code and the measurements

1. **The m-range when `log_lap < 0`** (step 1, "the same code with `lap = β^{-k}` and the inequalities unchanged"). With β < 1, `log(lo/c)/log_lap > log(hi/c)/log_lap`, so `range(floor(a), ceil(b)+1)` runs from the larger bound to the smaller. It is correct only while `hi/lo < β^{-k}`, that is, while the interval is shorter than one lap. That holds for every window pair (the window is one lap), so the corpus would not show the bug. A wider hand-built interval (for example a=5.5, b=2.5) would leave the range empty and silently skip every landing. **Use `min/max` of the two bounds** (one line; my probe does this).
2. **β.** Step 1 says to scale by `1/β` via `per_step_beta("unstable")`. The documented rule (`Trellis.image_cdist` Note) and the B0 law (`_scaled`) use **`per_step_beta("stable")`**. Both give the same value on these maps (`FixedPoint.py:280-290`). Using `per_step_beta(stability)` makes the closed-form `value` equal `image_cdist(r, steps, stability)` exactly as on the unstable side, and the code stays one body.
3. **The candidate set** is not named in the plan. It must stay `trellis.own_intersection_ids` (the cut N), gated by `pos_map.get(r.manifold_b_key)`. The mirrored 2-D guard applies only when `r.manifold_a_key[0] is ix0.manifold_a_key[0]` (the plan says "the 2-D collision guard" without mirroring D7's `same_stable`). See 3d for why widening the set breaks a law case.
4. **Unit tests, "each fails on base" (step 4, rubric 3).** The third test (the pair's own orbit, a 2-D collision, is still a pair) **passes on base**, because base keeps any pair whose unstable interval is clear. It can only fail on base if it calls `_interval_is_clear(..., "stable")` directly, in which case it fails at import. The second test (a table link, rejected by id) needs the linked landing to be **off the trellis's branches**. Otherwise `x0` and `x1` are not consecutive, because `helpers.fakes.synthetic_trellis` (`tests/helpers/fakes.py:77-113`) puts every registered crossing on its branches. So the test must build the `Trellis` with `resolved_ids` that exclude the landing, or with a cut. On real builds this reading cannot fire: since B1 the table links only N, and A.5 shows that a linked window visit is in T.
5. **e13 (step 5).** Four points:
   - "every series of `cases.SERIES` and the E12 depths": E12 also runs `nested_inner` and `heteroclinic` (not in SERIES) and two cutoffs. The plan should say which, if any, to include.
   - "Before/after" can only be computed on the branch, from the one function: before = the unstable call True; after = both calls True. The old function no longer exists there (rubric 4).
   - `reference_witness` (`tests/helpers/laws.py:1721`) explains a pair lost *between recomputes* using the table alone, so it returns None for a pair killed only by the scaled reading. e13 should also record the stable reading that fired and its witness id.
   - Writing to `results/e13_both_intervals.jsonl` adds a **fifth file** if committed (E12 committed its own). That conflicts with rubric 1's "only the four files", so either the rubric admits the results file or it stays uncommitted.
6. **Behaviour "none on the corpus"** is confirmed (3a-3c), and it extends to all E12-depth series steps at the recipe cutoff. I did **not** run the 1e-5 cutoff cells, `nested_inner`, `heteroclinic` blasting or the two-iteration cells (not measured).

## 6. Recommendation (the smallest change)

Implement the plan's step 1 as one body with a `stability` parameter, with three amendments: the m-range by `min/max`, `per_step_beta(stability)`, and candidates `own_intersection_ids` with the mirrored D7 guard. Edit the docstring at `:159-160` and the Dev Notes `:64-66`, `:83-89` and `:117-126`; `:50-58` can stay. Leave `tests/cases.py` alone (there is no KNOWN_ISSUES entry). Expect no change in any fingerprint. Make the third unit test call `_interval_is_clear` directly so that rubric 3 holds. Keep e13 to one file that prints its summary, and say whether its jsonl is committed.

## 7. Out of scope (noticed, not recommended here)

- The trellis crossing count stops growing on every deep series (k28 from blast 7, k32 from 8, p3 from 14, inversion from 9, nested_outer from 8). Past that point e13 measures nothing new; that is B3's job.
- `CLAUDE.md`'s line "a crossing past a cut is a point of X" reads literally as "every registered crossing past a cut", which would kill reference (2, 3) on `inversion` (3d). A.5's `X = ∪ M^n(N)` is the reading the code and the B0 law use; the CLAUDE.md wording could say "iterates of N outside T".
- The heteroclinic build and the deep series have a table landing outside the interval while scaling put it inside (`table_outside`: 4 on heteroclinic, up to 30 per series step). This is the usual scaling noise, and the table-first ordering absorbs it.

---
Knowledge files updated:
- `.refactor/knowledge/topology.md`: new section "B2 analyst pass"; the stale NUM-18 "deferred" note at :215 corrected.
- `.refactor/knowledge/testing.md`: new section "B2 analyst notes".
