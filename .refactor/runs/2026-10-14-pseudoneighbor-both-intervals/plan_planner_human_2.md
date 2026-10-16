# Plan — B2 `pseudoneighbor-both-intervals`

*2026-10-14. Base: `main` at `ba23a1f`. This is slice B2 of `docs/alignment_plan_2026-10-09.md`, with the analyst's five corrections folded in.*

## 1. The answer

The definition (`Pseudoneighbor_Algorithm.pdf`) says a consecutive pair `(x0, x1)` on the window is a pseudoneighbor pair only when **neither** of its open intervals holds a point of X: not the unstable arc `W^u(x0, x1)` and not the stable arc `W^s(x0, x1)`. The code checks one of them and assumes the other. `_unstable_interval_is_clear` (`Pseudoneighbor.py:329-471`) checks the unstable arc against every iterate of every crossing of N. For the stable arc the only check is that the walk steps through *consecutive* crossings of the cut trellis (`:209`, `:224`), so no crossing of N sits between the two, but no iterate is ever looked at. The Dev Notes state the assumption outright (`:64-66`, "empty by construction; all the work is in the unstable-interval check").

The analyst ran an independent copy of the stable check at `HEAD`. It fired on no law case, not on k28 x8 or p3 x6, and not on any of 73 series steps out to the E12 depths. Simulating a missing window representative did not change that: the stable check then catches the missing point about half the time, but the unstable check has already rejected every one of those pairs. So on our cases this change only completes the definition and changes no result.

## 2. Recommendation

Change it, as you decided (Q1, S-a): one function tests either interval, and a pair is kept only when both are clear. Expect no new input for DEC-25: e13 should report no differences.

## 3. Scope

Slug **`pseudoneighbor-both-intervals`**. The work happens on branch `refactor/pseudoneighbor-both-intervals` in worktree `../tangle-pack-pseudoneighbor-both-intervals`, made from `ba23a1f`. The base worktree for comparison is `../tangle-pack-pseudoneighbor-both-intervals-base`.

Files (four):

1. `src/tanglepack/topology/Pseudoneighbor.py`
2. `tests/unit/topology/test_pseudoneighbor.py`
3. `stress/numerics/e13_both_intervals.py` (new, report only)
4. `stress/numerics/results/e13_both_intervals.jsonl` (new; only if you answer "yes" to the question below)

`tests/cases.py` is off the list. It has no `KNOWN_ISSUES` entry for `reference_pairs_clear_of_X`, so the planned edit there does nothing. Estimated size: about 60 changed lines in `Pseudoneighbor.py` (net shorter), about 90 lines of tests and about 150 lines of e13. Your uncommitted `docs/` and `.refactor/` files are not touched.

## 4. Steps

**Step 1. One function, `_interval_is_clear(trellis, x0, x1, stability, candidates, collision_rtol)`.** It replaces `_unstable_interval_is_clear`. The body stays one body and reads the arc's side through the `stability` argument:

```python
ix0, ix1 = trellis.intersection(x0), trellis.intersection(x1)
arc_key = attrgetter("manifold_a_key" if stability == "unstable" else "manifold_b_key")
cdist = attrgetter(f"{stability}_cdist")
key = arc_key(ix0)
if stability == "stable":
    assert key == arc_key(ix1)       # the walk steps along one stable branch
elif key is None or key != arc_key(ix1):
    return False                     # no single unstable arc connects the pair (as today)
owner = key[0]
pos_map, k = owner.branch_position_map(stability)
beta = owner.per_step_beta(stability)    # < 1 for "stable": = image_cdist's factor
lo, hi = sorted((cdist(ix0), cdist(ix1)))
...
lap, log_lap = beta**k, k * math.log(beta)
a, b = math.log(lo / c) / log_lap, math.log(hi / c) / log_lap
for m in range(math.floor(min(a, b)), math.ceil(max(a, b)) + 1):
    ...
```

The analyst's three corrections are part of the step:

- **The m-range takes `min/max` of its two bounds.** With β < 1 the bounds swap. On a window pair this does no harm, because a window is one lap, but a wider interval would silently skip every landing.
- **β is `per_step_beta(stability)`.** For `"stable"` that equals `1/per_step_beta("unstable")`, so the closed-form value is exactly what `image_cdist(r, steps, "stable")` returns. That keeps one body for both sides.
- **The candidates stay `trellis.own_intersection_ids`, which is the cut N.** Each candidate is gated by `pos_map.get(arc_key(r))`, so a foreign branch never lands on the pair's branch. The 2-D collision guard applies only when the *other* side's fixed point matches. For the unstable reading that is today's `same_stable`. For the stable reading it is `r.manifold_a_key[0] is ix0.manifold_a_key[0]` (D7 mirrored). Widening the candidates to `registry` or `resolved_ids` is not allowed: the law case `inversion` would lose its reference (2, 3) to crossings past the unstable cut that are not in X.

The three readings stay in today's order:

1. **Table landing, by id.** The warning about the wrong branch reads `arc_key(z)`, and the interval test reads `cdist(z)`.
2. **`r` in its own frame,** against the pair's registered images.
3. **Scaling,** with the 2-D guard.

The anchor skip becomes `cdist(r) <= 0`, and the `lo <= 0` rejection covers both sides. `same_point` always takes its arguments as `(unstable, stable)`, so the stable reading passes them swapped back into that order: in the scaled reading, `(image_cdist(r, steps, "unstable"), value)`. Each `logger.debug("Pair ... punctured by ...")` line names the stability and the reading. e13 reads those lines (step 4).

**Step 2. The caller.** In `compute_pseudoneighbors`, keep the pair when

```python
all(_interval_is_clear(trellis, x0, x1, s, candidates, collision_rtol)
    for s in ("unstable", "stable"))
```

holds. Its docstring (`:159-160`) changes from "open unstable interval" to "neither open interval holds a point of X".

**Step 3. Dev Notes.** Three changes:

- Delete `:64-66` ("empty by construction").
- Rewrite the closed-form note (`:83-89`) and the check-set note (`:117-126`) once for both intervals: β^k per lap, which shrinks the stable cdist, hence min/max; and the D7 guard mirrored.
- State the definition once: both open intervals, against X.

The orphan-pair note (`:50-58`) stays. No paragraph in the file mentions NUM-18, so there is nothing to remove for it.

**Step 4. Unit tests.** All three use hand-built registries and the existing `_reference_window` helper.

1. **`test_scaled_iterate_punctures_stable_interval`.** Add a registered crossing `r` below the window on the pair's stable branch, with no table links. Its backward-scaled stable cdist lands strictly inside `(s(x1), s(x0))`. Choose `r`'s unstable cdist so the unstable reading stays clear. Expect `compute_pseudoneighbors` to return no pair `(x0, x1)`.
2. **`test_table_landing_punctures_stable_interval`.** Add a crossing `w` strictly inside the pair's stable interval and a candidate `r` with `M^n(r) = w` in the table. Give `w` an unstable cdist outside the pair's unstable interval. Then remove `w` from the stable branch's `intersection_ids` (a cut), so that `x0` and `x1` stay consecutive and `w` is not a candidate. Expect no pair, decided by id.
3. **`test_own_orbit_does_not_puncture_stable_interval`.** Add an unlinked crossing that 2-D-collides with an iterate of `x0`. The test imports `_interval_is_clear` **inside the test body**, asserts that `_interval_is_clear(trellis, x0, x1, "stable", ...)` is True, and asserts that the pair is still returned.

Tests 1 and 2 go through the public path only, so on base they fail on their assertion. Test 3 fails on base only because the function does not exist there; the import is local so that tests 1 and 2 still fail for the right reason.

**Step 5. `stress/numerics/e13_both_intervals.py`.** A harness module shaped like E12.

- **Cases:** all seven E12 cases at E12's depths and separations: `k28` 13, `k32` 13, `p3` 20, `nested_outer` 12, `nested_inner` 12, `inversion` 10, `heteroclinic` 8. That covers `cases.SERIES` plus `nested_inner` and `heteroclinic`. Each runs at both cutoffs, 1e-7 and 1e-5, with one iteration per call: 14 cells. Builds come from `e11.build(case)`. E12's two-iteration cells and its depth study are left out.
- **What each blast measures:** e13 walks every consecutive window pair on the branch and calls the one function twice. *Before* is the set of pairs where the unstable call is True; *after* is the set where both calls are True. Both come from the same code, since the old function no longer exists. For each pair in *before* but not in *after*, e13 records:
  - the pair;
  - the stable reading that fired (table, own frame or scaled), the witness id `r` and its `steps`, all read from the debug line of step 1 through a log handler;
  - `reference_witness(trellis, pair)` from `tests/helpers/laws.py`, which returns None for a pair killed by scaling alone.
- **Per step it also records** the seed, closed and hole counts and `n_trellis_crossings`, so the report shows where N stops growing. Blasts beyond that point measure the same N again; that is B3's concern.
- **Output:** one JSON line per cell, as in E12, plus a printed summary of the cells with any difference. The script is not wired into the report PDF.

## 5. Behaviour

**None on our cases.** Every law case's reference pairs, holes, partition, classes and words stay as they are. The new check can only reject pairs, never add them, so any change e13 finds in a cell the analyst did not measure (the 1e-5 cutoff, `nested_inner` and `heteroclinic` blasting) is a lost reference. That is a finding to report, not a failure. The 2-D guard still uses `collision_rtol` as it is; no new tolerance.

## 6. Not doing

- Deciding DEC-25 (the inversion series' drop at blast 6). e13 is your input for it.
- Any change to the reference window or to `extend_pseudoneighbor_trajectories`.
- Wiring e13 into `stress/numerics/report/` (`make_numbers.py`, `report.tex`).
- E12's two-iteration cells and depth study.
- The N stall on the deep series (B3).
- Correcting the CLAUDE.md wording "a crossing past a cut is a point of X" (should say "iterates of N outside T").
- Updating the old function name in `docs/`, `.refactor/knowledge/` and `docs/open_issues.md`.

## 7. Rubric

1. `MPLBACKEND=Agg env/bin/python -m pytest -q` passes on the branch. `-rxX` lists only the existing `orientation_reversing_case_builds` xfail.
2. Every test in `tests/unit/topology/test_pseudoneighbor.py` at base passes on both the base and the branch. So does `tests/invariants/test_law_partition.py`, which runs `check_reference_pairs_clear_of_X` on every law case with no `KNOWN_ISSUES` entry for it.
3. `git diff --stat ba23a1f..refactor/pseudoneighbor-both-intervals` lists only the files in section 3.
4. No behaviour change beyond section 5. For every case in `cases.LAW_CASES`, `law_tier.fingerprint(case)` and the set of reference pairs from `compute_pseudoneighbors` are equal on the base and the branch, compared by a script run in both worktrees.
5. `grep -rn "_unstable_interval_is_clear" src tests stress` finds nothing. `_interval_is_clear` exists once and has a `stability` parameter. `compute_pseudoneighbors` calls it for both `"unstable"` and `"stable"`.
6. The m-range is computed with `min`/`max` of its two bounds. A direct call with a stable interval wider than one lap, holding a scaled landing, returns False. The reviewer runs this check by hand; it is not committed.
7. The stable reading scales with `per_step_beta("stable")`, not `1/per_step_beta("unstable")`, and its candidates are `own_intersection_ids`. Reading the diff shows no `registry` or `resolved_ids` used as the candidate set.
8. The three new tests pass on the branch. Run against base `src` (`PYTHONPATH` pointing at the base worktree), tests 1 and 2 fail on their assertion, not at collection, and test 3 fails because `_interval_is_clear` is missing.
9. `grep -n "empty by construction" src/tanglepack/topology/Pseudoneighbor.py` finds nothing. The `compute_pseudoneighbors` docstring names both intervals.
10. `env/bin/python stress/numerics/e13_both_intervals.py --quick` runs, and every record holds per-step before and after reference sets, the lost pairs with reading, witness id and `reference_witness`, and the seed, closed and hole counts. The full run covers 14 cells (7 cases × 2 cutoffs) at the depths in step 5.

## 8. Question

**Should e13's `results/e13_both_intervals.jsonl` be committed?** It would be the fourth file. **Default: yes.** E12 committed its results (`a97fcb6`), and this file is what you read before deciding DEC-25. If you answer no, the file stays untracked, the file list drops to three, and the reviewer reruns e13 to see the report.
