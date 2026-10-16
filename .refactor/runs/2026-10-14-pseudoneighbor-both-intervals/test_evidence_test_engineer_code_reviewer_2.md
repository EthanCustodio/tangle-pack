# Test evidence: B2 `pseudoneighbor-both-intervals` (test engineer, round 1)

*2026-10-14. Plan: `artifacts/plan_planner_human_2.md` (rubric §7, behaviour "none on our cases" §5). Coder report: `artifacts/change_coder_test_engineer.t1.md`.*

| | |
|---|---|
| Base | `ba23a1f`, detached worktree `/home/dezhu/Development/TanglePack/tangle-pack-pseudoneighbor-both-intervals-base` (clean) |
| Branch | `refactor/pseudoneighbor-both-intervals`, coder head `dfdfb2d`, **now `9cd29d5`** (one test commit of mine, below) |
| Worktree | `/home/dezhu/Development/TanglePack/tangle-pack-pseudoneighbor-both-intervals` (clean) |
| Python | `/home/dezhu/Development/TanglePack/tangle-pack/env/bin/python`, `MPLBACKEND=Agg`. Pytest reads each worktree's own `src` (`pythonpath = ["src", "tests"]`), and my scripts put `<worktree>/src` first on `sys.path` and assert `Pseudoneighbor.__file__` is under it. |
| Scratch | `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/te_b2/` (scripts `r4.py`, `r6.py`, `mutate.py`, and every output quoted here) |

## Verdict in one paragraph

Every rubric item I re-ran holds: 4, 8, 6, 5, 2, 3, 7, 9, 10 (quick run), and the full suite. The coder's deviation is sound. In each of the seven original fixtures, a crossing of N has a stable orbit that visits the window strictly inside a pair's stable interval. It does not 2D-collide with either end, and the table does not link it to either end. So no implementation of the definition can keep that pair. In all seven adjusted fixtures the failing pair is now decided by the unstable reading the test is named for. The adjusted tests catch exactly the same breaks as the originals do on base. Two of the seven had no teeth on base already, and still have none; I report those as pre-existing gaps. I added one test, which the plan did not ask for: the mirrored D7 guard had no test, and the break that drops it survived the coder's suite. One break survives the committed suite: unswapping the m-range bounds. Only the hand check of rubric 6 catches it, and the plan keeps that check uncommitted.

## 1. Rubric 4: no behaviour change on the law cases (re-run)

My script `te_b2/r4.py` covers every case in `cases.LAW_CASES` and runs once in each worktree root. For each case it writes:

- `law_tier.fingerprint(case)`. In `products`, each per-process `id()` is replaced by the product's type name; the other aspects are compared as they are, with memory addresses stripped.
- the sorted `compute_pseudoneighbors(topology.trellis)`;
- every pseudoneighbor pair with its iterate;
- every class's `(label, letter, word)`.

```
(cd <base>   && python te_b2/r4.py > te_b2/r4_base.jsonl)
(cd <branch> && python te_b2/r4.py > te_b2/r4_branch.jsonl)
diff <(tail -n +2 r4_base.jsonl) <(tail -n +2 r4_branch.jsonl) && echo IDENTICAL
-> IDENTICAL
```

Line 1 of each file records which module was loaded:

- base: `.../tangle-pack-pseudoneighbor-both-intervals-base/src/.../Pseudoneighbor.py`, `has_interval_is_clear: false`;
- branch: `.../tangle-pack-pseudoneighbor-both-intervals/src/...`, `true`.

The reference pairs are identical on both, and match the coder's table:

| case | references | pn pairs | classes |
|---|---|---|---|
| k10 | (5,6) | 1 | 2 |
| k28_one_blast | (3,4) (4,5) | 4 | 3 |
| k28_two_blasts | (9,10) | 1 | 4 |
| p3 | (3,4) (4,5) | 24 | 5 |
| nested | (7,8) (8,9) (37,38) | 25 | 9 |
| nested_inner_blasts | (7,8) (8,9) (41,42) | 29 | 9 |
| inversion | (1,2) (2,3) | 8 | 4 |
| inversion_one_blast | (26,27) | 1 | 7 |
| inversion_two_blasts | (23,24) (26,27) | 3 | 8 |
| inversion_period2 | (13,14) (14,15) | 10 | 6 |
| heteroclinic | (77,78) (78,79) (92,93) (113,114) | 16 | 75 |

**Pass.**

## 2. Rubric 8: the three new tests (re-run)

- **Branch:** `pytest -q tests/unit/topology/test_pseudoneighbor.py` gives `34 passed` at the coder's head, and `36 passed` with my test added.
- **Base `src`:** I copied the branch's test file into the base worktree, ran it, then ran `git checkout` on the file. The base worktree is clean afterwards. Output in `te_b2/r8_base.txt`:

```
E       assert [(0, 1), (1, 2)] == [(0, 1)]          # test_scaled_iterate_punctures_stable_interval
E       assert [(0, 1), (1, 2)] == [(1, 2)]          # test_table_landing_punctures_stable_interval
E       ImportError: cannot import name '_interval_is_clear' from 'tanglepack.topology.Pseudoneighbor'
        (/home/dezhu/.../tangle-pack-pseudoneighbor-both-intervals-base/src/...)   # test_own_orbit_does_not_puncture_stable_interval
3 failed, 31 passed
```

Tests 1 and 2 fail on their assertion and test 3 fails on the import, as the plan requires. The file collects on base, and all 31 other tests pass. With my test added, the run gives `4 failed, 32 passed`: its `other_unstable` case also fails on its assertion (`assert ((1, 2) in [(0, 1), (1, 2)]) is False`), as a fix should. **Pass.**

Test 3's second assertion (the pair is kept) also passes on base, because base ignores the stable side. Only its direct call is evidence. Break `m6` below shows that the call has teeth.

## 3. Rubric 2: the test files on both sides, and the seven adjusted fixtures

- The **branch's** `test_pseudoneighbor.py` on **base src** gives 31 of 31 pre-existing tests passing, as shown in §2. The only failures are the new tests.
- `tests/invariants/test_law_partition.py` gives `11 passed, 4 skipped` on the base and `11 passed, 4 skipped` on the branch.
- The **base's** literal `test_pseudoneighbor.py` on **branch src** gives `7 failed, 24 passed`. I ran it with `--log-level=DEBUG` (`te_b2/r2_orig_on_branch.txt`), then restored with `git checkout`. The seven failures are the ones the coder lists.

### Each original fixture holds a real stable-interval puncture

The window is the same in all seven. With λ = 4, a stable cdist scales by 1/4 per step:

- r_n is at s = 4, x1 at s = 3 and r_end at s = 1.
- The stable intervals are (3, 4) for (r_n, x1) and (1, 3) for (x1, r_end).

A stable orbit `s·4^m` meets [1, 4) exactly once per lap. So a crossing on the window's stable branch avoids the open intervals only if s = 4^j or 3·4^j, that is, only if it sits at the ends.

The deciding debug line in each failing test, logged by the branch code:

| test (original data) | witness: stable cdist on `s(0.0)` | orbit visit | pair killed | debug line (branch) | 2D collision with an end? | table link to an end? |
|---|---|---|---|---|---|---|
| `direct_intersection_punctures_interval` | s = 0.5 (u 1.5) | 0.5·4 = **2** ∈ (1,3) | (x1, r_end) | `Pair (1, 2): stable interval punctured by 3 at M^-1, scaled reading (cdist 2 in (1, 3))` | no: (0.375, 2) vs (2,3), (4,1) | none |
| `iterate_punctures_interval` | s = 0.5 (u 0.4) | **2** ∈ (1,3) | (x1, r_end) | same line | no: (0.1, 2) | none |
| `walk_stops_at_reference_window_end` | `below`, s = 0.5 (u 16) | **2** ∈ (1,3) | (x1, r_end) | same line | no: (4, 2) vs r_end (4, 1) is 100% off on s | none |
| `heteroclinic_candidate_ignored` | s = 0.5; its unstable side is the other saddle's | **2** ∈ (1,3) | (x1, r_end) | same line | not applicable: the other side is foreign, so D7 gives no collision skip | none |
| `table_linked_deep_iterate_does_not_disqualify` (unlinked half) | n = 6 chain point, s = 3·4⁻⁶·1.02 | 3·1.02 = **3.06** ∈ (3,4) | (r_n, x1) | `Pair (0, 1): stable interval punctured by 8 at M^-6, scaled reading (cdist 3.06 in (3, 4))` | no: 2% off x1 on s, twice `collision_rtol`, by the test's own design | none (unlinked half) |
| `table_linked_landing_uses_registered_cdist` (linked half) | z, z1, r at s = 0.5, 0.125, 0.03125 | all land at **2** ∈ (1,3) | (x1, r_end) | `Pair (1, 2): stable interval punctured by 3 at M^-1, scaled reading (cdist 2 in (1, 3))` | no: (1.125, 2) | z → z1 → r is linked, but to no end |
| `unlinked_landing_is_read_against_the_pairs_images[below_c]` | s = 0.01 (u 128.2) | 0.01·4⁴ = **2.56** ∈ (1,3) | (x1, r_end) | `Pair (1, 2): stable interval punctured by 9 at M^-4, scaled reading (cdist 2.56 in (1, 3))` | no: (0.50, 2.56) | the images of x1 and r_end stop at +3, so none at +4 |

In every row:

- the witness is a crossing of N (it is in the trellis);
- its orbit's visit to the window is strictly inside the interval;
- it does not 2D-collide with x0 or x1 within the library's `collision_rtol`;
- the iterate table does not identify it with either end.

By the definition (`X = ∪ M^n(N)`, both open intervals), that pair is not a pseudoneighbor pair, so **no correct implementation could keep it**. The original assertions require it to be kept, so the fixtures, not the code, were wrong.

The deep-iterate row is the only one that rests on the collision tolerance. Its docstring defines the n = 6 point as "drifted 2·collision_rtol", which is not x1's orbit when unlinked, so the verdict follows the test's own premise.

### The adjusted fixtures exercise what their names say

Run on the branch with DEBUG on, the deciding line of each adjusted test is the unstable reading it is named for. No stable puncture fires in any of them:

```
direct          Pair (0, 1): unstable interval punctured by 3 at M^0, table reading (cdist 1.5 in (1, 2))
iterate         Pair (0, 1): unstable interval punctured by 3 at M^1, scaled reading (cdist 1.6 in (1, 2))
walk_stops      (no puncture)
heteroclinic    (no puncture)
deep_iterate    Pair (1, 2): unstable interval punctured by 8 at M^-6, scaled reading (cdist 2.002 in (2, 4))   [unlinked half]
table_linked    Pair (1, 2): unstable interval punctured by 5 at M^-2, scaled reading (cdist 3 in (2, 4))      [unlinked half]
unlinked        Pair (1, 2): unstable interval punctured by 9 at M^-3, own frame reading (cdist 130 in (128.5, 256))  [between_c_f]
```

For each fixture, I broke the reading the test is named for, on the base (original tests) and on the branch (adjusted tests). `te_b2/mutate.py` applies one replacement, runs the unit file and restores with `git checkout`. The tables below show which of the seven fail:

| break | base, original data | branch, adjusted data |
|---|---|---|
| `a5`: skip landings at `steps == 0` (direct) | `direct` fails | `direct` fails |
| `a6`: skip landings at `steps > 0` (forward iterate) | `iterate`, `branch_residue` fail | `iterate`, `branch_residue` fail |
| `a2`: `landing = None` (no table lookup) | `table_linked_landing` fails | `table_linked_landing` fails |
| `a3`: own-frame reading off | `unlinked[below_c]` fails | `unlinked[below_c]` fails |
| `a4`: `while x1 is not None` (walk does not stop at r_end) | `endpoint_iterates`, `deep_iterate` fail; **`walk_stops` passes** | the same two fail; **`walk_stops` passes** |
| foreign-unstable gate off, unstable side only (`a1b` base / `a1c` branch) | **all pass** | **all pass** |

The adjusted data keeps the originals' teeth exactly, and the coder's change loses nothing. The table also exposes two **pre-existing gaps**, identical on base. Neither was introduced here, and I did not change either test, since they are outside this round:

- **`test_walk_stops_at_reference_window_end` has no teeth.** If the walk ran past r_end, the extra pair (r_end, below) would still be rejected, because x1's orbit at u = 8 lies inside its unstable interval (4, 16). This is true with s = 0.5 and with s = 0.25 alike.
- **`test_heteroclinic_candidate_ignored` does not isolate its gate.** Its only landing is at M⁰, so `trellis.iterate(r, 0)` is the candidate itself. The "wrong branch" table check then discards it, and the foreign-manifold gate is never the deciding line.
- **`deep_iterate` and break `a2`.** Under `a2`, `deep_iterate` still passes on both versions. The own-frame `r_id in images` check identifies the linked n = 6 point by id as well, so the table lookup is redundant there. This is not a gap in what the test asserts, but the lookup is not the only path that test exercises.

The adjusted `deep_iterate` drift (1 − 2·rtol) puts the n = 6 point at s = 2.94. That is inside (1, 3) of (x1, r_end), the pair the unstable reading already kills. The linked half still keeps both pairs by id on both readings. `walk_stops`' new `below` sits at M(r_end)'s cdists, so its orbit meets the window at r_end itself, an end. `OTHER_STABLE = 1` is a key outside the non-inversion saddle's stable cycle, as in the existing `test_pair_on_different_unstable_branches_rejected`.

**Judgement: the deviation is correct and the adjusted tests remain valid.** Rubric 2's literal wording ("every test at base passes on both") cannot hold for these seven. I recommend accepting the change: the coder changed data only, not names or assertions. The author may rule on it if they want.

## 4. Rubric 6: a stable interval wider than one lap (hand check, not committed)

`te_b2/r6.py` builds a hand-built λ = 4 saddle:

- x0 at s = 40 and x1 at s = 0.5, an interval about three laps wide;
- a candidate at (u 7, s 0.01), whose scaled landings 0.64, 2.56 and 10.24 lie inside;
- a call to `_interval_is_clear(..., "stable", own_intersection_ids, 1e-2)`.

```
branch:  with candidate: False | without: True   RUBRIC6 OK
```

**Pass.**

## 5. Breaks of the new code (all restored with `git checkout`; the worktree is clean after each)

| # | break | what fails |
|---|---|---|
| m1 | unswapped m-range: `range(floor(bounds[0]), ceil(bounds[1]) + 1)` | **rubric 6 check: `with candidate: True`, AssertionError.** The unit file and `test_law_partition.py` still pass (`47 passed, 4 skipped`): every window pair spans at most one lap, where the unswapped range is still non-empty. **The committed suite does not catch this break.** The plan keeps the check by hand (rubric 6 says "not committed"), so I left it uncommitted; the reviewer may want it committed. |
| m2 | `beta = owner.per_step_beta("unstable")` on both sides | 5 fail: `test_scaled_iterate_punctures_stable_interval`, `..._another_unstable_manifold..._stable_interval[other_unstable]`, and the pre-existing `endpoint_iterates_do_not_disqualify`, `single_reference_window_per_fixed_point`, `table_linked_deep_iterate_does_not_disqualify` |
| m3 | candidates = `sorted(trellis.resolved_ids)` | law case **`inversion`: references `[(1,2),(2,3)]` → `[]`**, so (2, 3) is lost as findings §3d predicts, and (1, 2) with it. Its holes and words change too (`r4.py inversion` under the break). `test_law_partition[inversion]` and `[inversion_period2]` fail (14 laws each), as do the unit tests `reference_window_starts_at_chosen_strong_pip` and `missing_pip_iterates_are_cut...[period_3_placed]` |
| m4 | mirrored D7 guard dropped: `same_other = stability == "stable" or ...` | **at the coder's head: nothing fails** (the three planned tests are all single-saddle). With my new test, `..._another_unstable_manifold..._stable_interval[other_unstable]` fails |
| m5 | `for stability in ("unstable",)` (stable call dropped) | `test_scaled_iterate_punctures_stable_interval`, `test_table_landing_punctures_stable_interval`, `[other_unstable]` fail |
| m6 | stable side never skips a 2D collision (`if stability == "unstable" and same_other and any(...)`) | `test_own_orbit_does_not_puncture_stable_interval`, `[own_orbit]` fail |

Outputs are in `te_b2/mutations_a.txt`, `te_b2/mutations_m3.txt`, `te_b2/mutations_adjusted.txt`, `te_b2/mutations_base.txt` and `te_b2/mutations_a6.txt`.

## 6. The test I added (commit `9cd29d5`, inside the plan's file 2)

`test_collision_on_another_unstable_manifold_still_punctures_stable_interval[own_orbit|other_unstable]` in `tests/unit/topology/test_pseudoneighbor.py`. It is the mirror of the existing `test_collision_on_another_stable_manifold_still_punctures`.

- **Setup:** a heteroclinic window of W^u(A) crossings on B's stable branch, with A at λ 4 and B at λ 6. The stable intervals are (3, 6) and (1, 3).
- **Witness:** a crossing on that branch at s = 0.4985. One step back it lands at 2.991, inside (1, 3). Its unstable cdist, scaled one step back by its own saddle, is 2.002, which matches x1's (2, 3) within `collision_rtol` by the numbers.
- **`own_orbit`:** the witness is on W^u(A). It is x1's orbit up to noise, so (x1, r_end) survives.
- **`other_unstable`:** the witness is B × B, on W^u(B). Its orbit never meets x1's, so the pair must go. The debug line is `Pair (1, 2): stable interval punctured by 3 at M^-1, scaled reading (cdist 2.991 in (1, 3))`.
- **What it guards:** the stable reading's `same_other` (rubric 7's "D7 mirrored"), against m4, and the stable collision skip, against m6.
- **On base:** `own_orbit` passes, because base ignores the stable side. `other_unstable` fails on its assertion, as a fix should.

Commit message: "Tests: the stable reading's 2D collision skips only an endpoint's orbit on the pair's own unstable manifold (B2, D7 mirrored)".

## 7. The other rubric items (re-run on the branch at `9cd29d5`)

- **1. Full suite.** `MPLBACKEND=Agg env/bin/python -m pytest -q -rxX` gives `1183 passed, 38 skipped, 5 deselected, 1 xfailed in 200.66s` (`te_b2/suite_branch.txt`). The 1183 are the coder's 1181 plus my two parameter cases. The only xfail is `test_orientation_reversing_case_builds[orientation_reversing]`. **Pass.**
- **3. Diff stat.** `git diff --stat ba23a1f..HEAD` lists four files and nothing else: `Pseudoneighbor.py` 238, `e13_both_intervals.py` 198, `results/e13_both_intervals.jsonl` 14 and `test_pseudoneighbor.py` 144. **Pass.**
- **5. Old name gone.**
  - `grep -rn "_unstable_interval_is_clear" src tests stress` finds nothing (exit 1).
  - `def _interval_is_clear(` appears once, at `Pseudoneighbor.py:343`, with a `stability` parameter.
  - It is called once, at `:229`, inside `all(... for stability in ("unstable", "stable"))`.
  - **Pass.**
- **7. β and candidates.**
  - In the diff's added `src` lines, `per_step_beta` appears only as `beta = owner.per_step_beta(stability)` (`:234`) and in the Dev Notes. The Dev Note mentioning `1/per_step_beta("unstable")` says this form is never used.
  - The candidates are still `trellis.own_intersection_ids` (`:208`).
  - `registry` appears in the diff only in the Dev Note "never the whole registry"; `resolved_ids` does not appear.
  - m2 and m3 show that both choices are load-bearing.
  - **Pass.**
- **9. Wording.** `grep -n "empty by construction"` finds nothing (exit 1). The docstring of `compute_pseudoneighbors` names "the unstable one W^U(x0, x1) and the stable one W^S(x0, x1)". **Pass.**
- **10. e13.**
  - `env/bin/python stress/numerics/e13_both_intervals.py --quick` exits 0 with 2 cells (p3 × 2, k28 × 3 at 1e-7) and prints `lost 0 stable fires 0 library mismatches 0 topology errors 0`. It writes an untracked `results/e13_both_intervals_quick.jsonl`, which I moved out to `te_b2/` so the branch stays clean.
  - **The committed full results** have 14 cells, 7 cases × {1e-5, 1e-7}, all `ok`. The step counts are depth + 1: p3 21, k28 and k32 14, nested_outer and nested_inner 13, inversion 11, heteroclinic 9.
  - Each `series` holds per-step `before`, `after`, `lost`, `n_stable_fires`, `matches_library`, `n_classes`, `n_closed`, `n_holes`, `n_trellis_crossings` and `topology_error`.
  - Summed over all cells: lost 0, stable fires 0, library mismatches 0, and `before == after` in every cell.
  - I did not re-run the full 14-cell sweep, which is about 100 s and 5 GB for k32 at 1e-7 alone and was not on my re-run list. I checked the committed file's shape and its totals instead.
  - **Pass.**

## 8. What I could not test, or chose not to

- **m1 (unswapped bounds)** is caught only by the rubric 6 hand check. Every interval the suite sees is at most one lap, so the committed suite cannot catch it. Committing the check is a one-function change in file 2 if the reviewer or author wants it; the plan says not to.
- **No real build exercises the stable clause.** e13 records zero stable fires on 182 steps, so the new path is evidenced by the hand-built tests alone. That is the expected result (§5 of the plan).
- **The two pre-existing gaps** of §3 (`walk_stops` and `heteroclinic_candidate_ignored`) were left alone. They are identical on base, and fixing them is outside this round.
- I did not re-run the full e13 sweep (see rubric 10).
