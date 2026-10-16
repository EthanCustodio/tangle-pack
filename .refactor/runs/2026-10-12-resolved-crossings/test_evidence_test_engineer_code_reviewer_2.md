# Test evidence: B1 `resolved-crossings` (test_engineer, round 2)

*Test engineer, 2026-10-12.*

**Sources:**
- Plan: `artifacts/plan_planner_human_2.md`.
- The author's round-1 decision at the end of `.refactor/decisions.md`: the rule stands, and the scope widens by four test files.
- The coder's round-2 report: `artifacts/change_coder_test_engineer.t2.md`.
- My round-1 evidence. This file replaces it, and §6 carries its unchanged findings forward.

**Ground rules this round:**
- I edited no library code and none of the coder's four rewritten tests.
- My new tests live in two files that are already on the list.

## 0. Summary

- **Full suite on the branch is green**, at both the coder's head and mine.

  | head | result |
  |---|---|
  | coder `b1a7e9c` | `1154 passed, 38 skipped, 5 deselected, 1 xfailed` |
  | mine `76d0710` | `1157 passed, 38 skipped, 5 deselected, 1 xfailed` |

  The only xfail is the `KNOWN_ISSUES` C3/INV-11 entry. Rubric 1 and 2 now hold.
- **The four rewritten tests and the deep-run law pass on the base too.** I copied them, with the branch's `cases.py`, into the base worktree. All 8 items pass (§2). So the reworked builds still describe behaviour the old code had.
- **Each rewritten test still guards its defect.** For each one I put back a form of the original defect on the branch, and the test failed (§4).
  - Coarse blasts: putting back the old geometric "approaches … from opposite sides" assertion fails the test. Flipping the sign rule for openings also fails it.
  - Order flip: if "order flip" is never diagnosed, the cause becomes `'walk'` and the test fails.
  - p2inv10 walk test: with the cut made a no-op, the test fails with `('n', 'unresolved', 'dual-graph walk unreachable')`.
- **`threshold = None` keeps clause 2 on.** Two new unit tests cover this, and each fails under the matching mutant. On a real build, the order-flip build at `None` still holds the clause-2 run {75, 76, 77} out of N. Those ids are in no trellis and end no bridge (§3.3).
- **Coverage gap closed: the oracle's new N clause was exercised nowhere in the suite.** No law case has an unresolved detection, and `inversion_period2_10` is now built at 1e-6, so it has none either. Dropping the clause left the whole deep-run law green.
  - New test `test_p2inv10_coarse_cut_holds_what_is_resolved` keeps the 1e-5 build, whose 38 and 39 are out of N and lie inside both pip-orbit cuts.
  - It fails when the oracle reads the registry, when the trellis holds unresolved detections, and when the cut is dropped.
- **Two weaker spots, reported:**
  - `test_p2inv10_every_class_walks` alone does not catch dropping the STABLE half of the cut. Four other tests in the same run do (§4, M3).
  - Tests (a) and (b) are not reachable under the default threshold. This is the coder's open question for the author, and §5 says the same.

## 1. Where

| | |
|---|---|
| Branch | `refactor/resolved-crossings` in `/home/dezhu/Development/TanglePack/tangle-pack-resolved-crossings` |
| Base | `/home/dezhu/Development/TanglePack/tangle-pack-resolved-crossings-base`, detached at `1029fcc`. It was clean before and after every copy (`git status --short` empty). |
| Coder head | `b1a7e9c` (round 2: `cb5cf92`, `e044f3d`, `e8c8105`, `f17a20f`, `b1a7e9c`) |
| My commit | `76d0710` "Tests: threshold None keeps clause 2 and admits a detection below the threshold; the pip-orbit cut at 1e-5 holds exactly T(q0) in N (B1 test evidence, round 2)" |
| `git diff --name-only 1029fcc \| wc -l` | **18**: the 14 planned paths plus the 4 the author added (`18 files changed, 961 insertions(+), 294 deletions(-)`) |
| Python | `/home/dezhu/Development/TanglePack/tangle-pack/env/bin/python` (written `$PY` below) |
| Scratch | `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/te2/` (written `$S/` below) |

## 2. Base run of the rewritten tests (pinning)

I copied the branch's `tests/cases.py` and the three rewritten regression files into the base. The base keeps its own `tests/helpers/laws.py`, because the branch oracle reads `Trellis.resolved_ids`, which the base does not have.

```
cd <base> && MPLBACKEND=Agg $PY -m pytest -q -p no:cacheprovider -rA \
  tests/regression/test_coarse_blasts_rows_by_signs.py tests/regression/test_evidence_diagnosis.py \
  tests/regression/test_pip_orbit_cut.py tests/regression/test_open_p3_deep_runs.py \
  -k "not p3_16 and not p3_6 and not p3_15"
```

Result: **8 passed** (`$S/base_rewritten.out`):
- `test_coarse_blasts_rows_by_signs`;
- `test_disagreeing_evidence_carries_its_cause`;
- `test_an_order_flip_gives_no_verdict`, the new hand-blasted build with `threshold=None`;
- the four `test_pip_orbit_cut` tests, with p2inv10 at 1e-6 and its new precondition;
- `test_every_law_holds_on_the_deep_runs[inversion_period2_10]` at 1e-6.

So:
- the reworked builds (clause 1 off, or refined to 1e-6) exercise the same defects on the old code as on the new;
- the 1e-6 deep run passes the base's own oracle.

I restored the base with `git checkout -- tests`; `git status --short` is empty.

## 3. New tests (commit `76d0710`)

Plain pytest functions on real objects (the hand-placed polyline workbench, `build_inversion_period2`, `laws`), no mocks, no pinned count, id or word.

### 3.1 `tests/unit/numerics/test_resolved_crossings.py`

| test | guards | base `1029fcc` | branch |
|---|---|---|---|
| `test_without_a_threshold_a_run_is_still_out_of_n` | **`None` turns clause 1 off only.** With `workbench.resolution_ratio_threshold = None`, the zig-zag's three detections are still `unresolved_crossings(None)`. The bridges are cut only at the two clean crossings, and `Trellis.from_workbench(...).resolved_ids == registry.resolved_ids(None)`. | FAIL: `TypeError: '<' not supported between instances of 'float' and 'NoneType'` (the base's `unresolved_crossings` cannot take `None`) | PASS |
| `test_without_a_threshold_a_detection_below_it_is_in_n` | **`None` does turn clause 1 off.** The grazing detection, which is out of N at 3, is in N at `None`, and `unresolved_crossings(None) == []`. | FAIL: the same `TypeError` | PASS |

### 3.2 `tests/regression/test_pip_orbit_cut.py`

| test | guards | base | branch |
|---|---|---|---|
| `test_p2inv10_coarse_cut_holds_what_is_resolved` | **The round's oracle change.** It runs on `build_inversion_period2(10, 10)` at the law case's 1e-5, partitioned. Preconditions: some detection is out of N, and a pip is chosen. It then checks that no unresolved detection is in `trellis(fp).own_intersection_ids`, and that `laws.check_trellis_is_cut_at_the_pip_orbit(case) > 0`. | FAIL: `assert not ({38, 39} & {0, 1, 2, …})` (the base trellis holds the unresolved 38 and 39) | PASS |

Why this test is needed:
- The probe `$S/probe.py` on the branch shows unresolved `[38, 39]` with ratios `[1.031, 2.328]`. Both lie inside both pip-orbit cuts:
  - stable `s(1,0)`: 1.645 and 1.663 against the cut at 1.749;
  - unstable `u(1,1)`: 3.172 and 3.211 against the cut at 3.377.

  So the oracle's skip is load-bearing on this build. The law checks 86 crossings there.
- No other test in the suite runs the oracle on a build with an unresolved detection. M1 below confirms this.

### 3.3 `threshold = None` on the two real builds (probe, not committed)

`$S/probe_none.py`, run on the branch:

```
p3_coarse_blasts  threshold None | registry 48  | out of N under None (clause 2) []           | below 3 only 4 | clause-2 ids in a trellis [] | bridge ends in runs []
order_flip_build  threshold None | registry 410 | out of N under None (clause 2) [75, 76, 77] | below 3 only 4 | clause-2 ids in a trellis [] | bridge ends in runs []
```

- In the order-flip build, clause 2 still removes the near-coincident run {75, 76, 77} under `None`. This is the same run the default-threshold k28 series flags at blast 6. It is in no trellis and ends no bridge.
- The coarse-blasts build has no run, so `None` changes only clause 1 there. Its 4 detections below 3 are the ones the test needs.

## 4. Break-it checks (branch)

Script: `$S/mutate.py <name>`. It makes one exact-text replacement (asserting the text occurs exactly once), runs the named tests, then restores with `git checkout -- <file>`. `git status --short src tests/helpers` was clean after every run. Output: `$S/mutants.out`.

| # | mutation | guards the defect of | tests run | result |
|---|---|---|---|---|
| M1 | `laws.check_trellis_is_cut_at_the_pip_orbit`: drop `or iid not in trellis.resolved_ids` (the oracle reads the registry again) | rewritten oracle (c) | new coarse-cut test; deep-run law `[inversion_period2_10]` | **new test FAILS** `(38, False, 3.1717…, 1.6453…)`, which is the reviewer's round-1 failure. The deep-run law **passes**: at 1e-6 nothing is out of N, so before my test the clause was unguarded. |
| M2 | `Trellis.cut_at_pip_orbit`: `kept = {… if inside(i)}` → keep every crossing (no cut) | (c), (d) | all of `test_pip_orbit_cut.py`; deep-run law | **6 FAIL**: k10 cut, p2inv10 cut=trim, the new coarse-cut test, the deep-run law `(7, True, 1.488…, 0.787…)`, mid-arc join, and **`test_p2inv10_every_class_walks` with `('n', 'unresolved', 'dual-graph walk unreachable')`**. This is a walk failure, stronger than the "no classes" the coder saw. |
| M3 | `cut_at_pip_orbit.inside`: drop the STABLE cut, keep the unstable one (the stable manifold running past f(S), the case the p2inv10 test was written for) | (d) | same | **4 FAIL**: k10 cut, the new coarse-cut test, mid-arc join, and the deep-run law `(11, True, 0.798…, 1.467…)`. **`test_p2inv10_every_class_walks` and `test_p2inv10_cut_equals_a_trim` PASS**, so the walk test alone does not catch a half-cut. The file and the deep-run law do. |
| M4 | `Trellis.from_workbench`: drop `if ix_id not in resolved: continue` | the oracle's premise | new coarse-cut test | **FAIL** `assert not ({38, 39} & {…})` |
| M5 | `IntersectionRegistry.resolved_ids`: `iid not in runs` → `(threshold is None or iid not in runs)` (`None` also turns clause 2 off) | "None keeps clause 2" | both new `None` tests | **FAIL**: `test_without_a_threshold_a_run_is_still_out_of_n` with `assert [] == [2, 3, 4]`. The grazing test passes, as it should, because it has no run. |
| M6 | `resolved_ids`: `None` read as 3.0 (`None` turns nothing off) | "None turns clause 1 off"; the premise of (a) and (b) | both `None` tests, coarse blasts, order flip | **3 FAIL**: the grazing `None` test `assert [1] == []`; coarse blasts "the builder no longer makes an under-resolved hole bridge"; order flip "precondition: a registered image reverses a member end's stable order". These are exactly the reviewer's round-1 failures (a) and (b). |
| M7 | `StablePartition` (`geometry_flags`): put back the pre-D6 hard assertion `assert not any(d.geometry_disagrees …), 'approaches its two defining crossings from opposite sides'` | rewritten (a), NUM-1 | coarse blasts | **FAIL** `AssertionError: approaches its two defining crossings from opposite sides`, the original defect's message |
| M8 | `StablePartition._hole_openings`: `(sign > 0) == (side == "left")` → `!=` (openings no longer follow the signs) | rewritten (a), D6 | coarse blasts | **FAIL** `((42, 41), 41, 'outward', 'right')`, on the hole bridge (42, 41), which the coder's rework targets |
| M9 | `SymbolicDynamics`: `if flip:` → `if False:` (an order flip is never diagnosed) | rewritten (b), C21 | order flip | **FAIL** `assert 'walk' == 'order flip'` |

Round-1 mutants 1-15 touch only `src/` files this round did not change (§6), so I re-ran them through the full suite only, not one by one.

## 5. Full suite

| run | command (from the worktree root) | result |
|---|---|---|
| branch `b1a7e9c` (coder) | `MPLBACKEND=Agg $PY -m pytest -q -rxX -p no:cacheprovider` | `1154 passed, 38 skipped, 5 deselected, 1 xfailed in 196.60s` (`$S/suite_b1a7e9c.out`) |
| branch `76d0710` (mine) | same | `1157 passed, 38 skipped, 5 deselected, 1 xfailed in 199.37s` (`$S/suite_mine.out`) |

On both runs the xfail is `test_orientation_reversing_case_builds[orientation_reversing]` (C3: INV-11), the one `KNOWN_ISSUES` entry. Round 1's base result stands: `1117 passed, 38 skipped, 5 deselected, 1 xfailed`.

No existing test was weakened or deleted by me.

The coder's four edits change builds and preconditions, not assertions. I checked each one in `git diff dc513f7..b1a7e9c`:
- the oracle gains one skip condition;
- p2inv10 gains a precondition and refines to 1e-6;
- coarse blasts and the order flip each set `resolution_ratio_threshold = None` in the build;
- every assertion line is unchanged.

**Whether tests (a) and (b) count as "weakened" is for the reviewer and the author.**
- **What they now prove:** NUM-1 and C21 still hold *on detections the default N excludes*.
- **Why that is not nothing:** M6 shows that under the default threshold neither defect is reachable. The coder's scan of `area_cutoff` and blasts found no build that reaches either one.
- **Why the guard is still real:** the guarded code (the partition's sign rule, the evidence verdict) runs downstream of N whatever the threshold. M7, M8 and M9 show both tests still fail when that code regresses.
- **Alternative:** the coder offered to mark them NOT_APPLICABLE or to re-guard them at unit level. That is an author choice, and I did not make it.

## 6. Carried forward from round 1 (unchanged `src/`)

Nothing under `src/` moved in round 2 (`git diff dc513f7..b1a7e9c -- src` is empty). So the following still hold:
- the base/branch byte-identity of all 11 law cases;
- the series comparison (k32 x9-10: holes 58 → 56, closed alphabet 39 → 45, `is_reliable` False → True; nested_outer x9-11: one fewer bridge; k28 x8: 24/15/23 → 4/4/4), which the author has now accepted as the rule's consequence;
- round-1 mutants 1-15 and the committed tests of `dc513f7`.

The law-case pin `test_every_law_case_crossing_is_in_n` passes in the green suite. The `cases.py` change keeps `build_inversion_period2`'s default at 1e-5, so the law case and its fingerprint guard are untouched (the guard passes inside the suite).

## 7. Not tested, and why

- **Defects (a) and (b) under the default threshold:** not reachable by any build the coder or I found (M6). See §5.
- **`test_p2inv10_every_class_walks` against a stable-only regression (M3):** it passes under that mutant. I did not strengthen it, because the rule says never edit the coder's rewritten assertions to change what they claim. The half-cut is caught by `test_k10_cut_at_every_candidate_equals_a_trim`, the new coarse-cut test, the mid-arc join test and the deep-run law.
- **The pre-D10 code itself** (the original p2inv10 defect) was not reconstructed. M2 and M3 stand in for it.
- **`Trellis.py:204` default for the sparse snapshot** (`resolved_ids=None` → the whole registry): not tested. As the coder notes, no reader of N runs on the sparse snapshot, and the fix belongs in `MinimalTrellis.py`, which is outside the list.
- **WARNING wording:** still not pinned, by policy.
