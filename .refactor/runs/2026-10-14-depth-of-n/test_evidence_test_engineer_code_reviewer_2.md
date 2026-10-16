# Test evidence: B3 `depth-of-n` (test engineer, round 1)

*2026-10-14. Plan: `artifacts/plan_planner_human_2.md` (option (a), Q2 "non-decreasing"). Coder report: `artifacts/change_coder_test_engineer.t1.md`.*

## Where

| | |
|---|---|
| Base | `9cd29d5`, worktree `/home/dezhu/Development/TanglePack/tangle-pack-depth-of-n-base` (detached, clean) |
| Branch | `refactor/depth-of-n`, worktree `/home/dezhu/Development/TanglePack/tangle-pack-depth-of-n` (clean after my commit) |
| My scripts and outputs | `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/te_b3/` (`measure.py`, `run.sh`, `brk.sh`, `break_law.py`, `{base,branch}_*.{txt,err}`, `brk_*`, `b6/`) |

I wrote my own measurement script (`te_b3/measure.py`) and did not reuse the coder's. It differs from the coder's in two ways:
- **Intersection objects** are canonicalised as `(id, unstable_cdist, stable_cdist)`, not by `repr` with the address stripped.
- **The law fingerprint is stronger.** It adds every bridge's node count and coordinate sum.

Every run used `PYTHONPATH=<worktree>/src:<branch>/tests`. Each `.err` file opens with the `tanglepack/__init__.py` it loaded, which confirms the base runs read the base library.

**My commit on the branch:**
```
7d3a294 Tests: the pip-orbit exemption skips only the separation test; every kept child is the blasted fixed point's, in the zone and never re-fed (B3, rubric 8)
```
`git diff --stat 9cd29d5..` after that commit:
```
 src/tanglepack/loom/Blast.py                       | 132 +++++++++++++++------
 tests/cases.py                                     |  15 +++
 .../regression/test_class_count_under_blasting.py  |  13 +-
 tests/unit/loom/test_blast_pip_orbit.py            |  84 +++++++++++++
 tests/unit/loom/test_blast_proximity_guard.py      |  43 +++++--
 5 files changed, 237 insertions(+), 50 deletions(-)
```
- Still the plan's five files.
- 287 changed lines. The cap is 300; my test used 23 lines of that margin.

## Rubric 1: full suite on the branch

`MPLBACKEND=Agg env/bin/python -m pytest -q -rxX` (from the branch worktree, with my test included):
```
XFAIL tests/invariants/test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing] - C3: INV-11: ...
1193 passed, 38 skipped, 5 deselected, 1 xfailed in 207.42s
```
- The single xfail is `cases.KNOWN_ISSUES[("orientation_reversing_case_builds", "orientation_reversing")]` (`tests/cases.py:1338`).
- 1193 = the coder's 1192 plus my one new test.

## Rubric 2: the branch's tests on the base source

**Setup.**
- `git archive` of the base `HEAD` into `te_b3/basesrc_branchtests/`.
- Over it, `git archive` of the branch `tests/` (including my test).
- A guard test, `test_zz_which_src.py`, asserts that `tanglepack.__file__` is under that directory. It passed, so the base library was the one under test.

**Command:** `pytest -q tests/unit/loom/test_zz_which_src.py tests/unit/loom/test_blast_pip_orbit.py tests/unit/loom/test_blast_proximity_guard.py "tests/regression/test_class_count_under_blasting.py::test_pip_orbit_never_shallower"`, run without my test (A) and then with it (B):
```
(A) FAILED test_blast_pip_orbit.py::test_pip_bearing_child_is_kept
    FAILED test_blast_pip_orbit.py::test_pip_depth_is_the_registered_orbit   (AttributeError: no attribute 'pip_depth')
    2 failed, 10 passed in 73.58s
(B) test_blast_pip_orbit.py only: 2 failed (the same two), 3 passed
```

| test | base | branch |
|---|---|---|
| 6(a) `test_pip_bearing_child_is_kept` | **FAIL** (the child ending at `f^4(q0)`, orbit `[4, 6, 35, 39, 47]` on base, is in `discarded_too_close`) | pass |
| 6(b) `test_pip_bearing_child_is_dropped_without_a_pip` | pass | pass |
| 6(c) `test_pip_depth_is_the_registered_orbit` | **FAIL** (AttributeError) | pass |
| `test_min_separation_drops_close_bridges[pip]` / `[no_pip]` | pass / pass | pass / pass |
| `test_min_separation_none_is_unchanged_default` | pass | pass |
| `test_pip_orbit_never_shallower` × 5 series | pass ×5 | pass ×5 |
| **new** `test_exemption_skips_only_the_separation_test` | pass (it pins old behaviour, so it is valid evidence) | pass |

This matches the rubric:
- 6(a) and 6(c) fail on base;
- 6(b), the restated guard test and the series law pass on both.

On 6(c): on base it fails only through the missing attribute. Its real teeth are shown by break 2 below.

## Rubric 3

As above: the five files, 287 changed lines.

## Rubric 4: law-case fingerprints (re-run)

**Method.** `te_b3/measure.py law` on base and on branch.
- It takes `helpers.law_tier.fingerprint` of every `cases.LAW_CASES` build.
- It drops `products`, which holds `id()`s of cached objects and differs per process by construction.
- It adds `bridge points`: id, node count and coordinate sum of every bridge.

**Result:**
- `cmp` of the summaries: identical.
- `cmp` of the full canonical JSON dumps (minus the header line naming the source path): identical.
```
k10 88a02638b9533872            k28_one_blast 6fe1477d564dadee   k28_two_blasts b25cd5c6fd692efb
p3 1606967469ac9718             nested b351fc5e4a66a204          nested_inner_blasts 282910d8d7da2edc
inversion 5c093e9c593c2da8      inversion_one_blast 13ae23074fb0bd09  inversion_two_blasts 0d9feb8cd2beaf6c
inversion_period2 ae5242f566187010  heteroclinic 340a485e001fa198
```
These hashes differ from the coder's because my canonicalisation and fields differ. Base and branch are equal under both.

## Rubric 5: every `SeriesStep` field, all five series

**Method.** `measure.py series <name>`, for each of `k28_series`, `p3_series`, `inversion_series`, `k32_series` and `nested_outer_series`.
- It dumps every dataclass field of every step as canonical JSON: `blasts, holes, classes, closed, ambiguous, active, inert, outside_minimal, references, unwitnessed, letters, pip_depth, history_free, crossings, recomputed`.
- The base side is measured with the branch's `cases.py`. Its only change is the test-side `pip_depth` composition, so on base it reads the base library.

**Result.** The per-step hash summaries and the full dumps are byte-identical on base and branch for all five series.

## Rubric 6: `pip_depth` per step (measured, base = branch)

| series | blasts | `pip_depth` |
|---|---|---|
| k28 | 0..8 | `[1, 2, 3, 4, 5, 6, 7, 8, 8]` |
| k32 | 0..10 | `[1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9]` |
| nested_outer | 0..11 | `[1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9, 9]` |
| p3 | 0..9 | `[14, 15, 16, 17, 18, 19, 20, 21, 22, 23]` |
| inversion | 0..4 | `[5, 6, 7, 8, 9]` |

These equal the plan's table, including the k28, k32 and nested_outer plateaus.

## Rubric 7: one 11-iteration blast on nested_outer

**Method.** `measure.py multi`:
- builds `SERIES["nested_outer_series"]`, takes the largest zone, and calls `blast_zone(zone, num_iterations=11, fixed_point=[fp], min_separation=1e-4)`;
- reads the depth before and after with `cases._pip_depth`, which composes through `trellis.iterate`;
- captures every INFO and WARNING record under the `tanglepack` logger.
```
base:   composed depth before 1 after 7 result.pip_depth <no attribute> completed 11 too_close 76 registry 486   exempt INFO lines 0
branch: composed depth before 1 after 9 result.pip_depth {'A': 9}       completed 11 too_close 78 registry 498   exempt INFO lines 8
```
- **Rubric 7 holds:** depth 7 on base, and `result.pip_depth[fp] == 9` on the branch, equal to the test-side composition.
- **The 8 exemption lines** name children `(34,35) ... (258,259)`, ending at `f^2(q0) ... f^9(q0)`. That is exactly one line per orbit point, with no duplicates.

## Rubric 8: reading `Blast.py` (branch, lines 420-470)

- **Order in the `for child in produced:` loop.** The checks run in this order:
  1. the `seen` check;
  2. the exterior check `if not (not child.iterated and fp_ok(child) and in_zone(child))`;
  3. `pip_end = next(...)` with its INFO branch;
  4. `elif min_separation is not None:`, the separation test.
- **Composition, not `forward_depth`.**
  - The orbit is composed by `_pip_orbits` through `iterate_table[x, 1]` link by link.
  - `grep forward_depth|forward_chain src/tanglepack/loom/Blast.py` finds only line 116, the docstring saying why `forward_depth` is not used.
- **Logging.** One `logger.info` per exempted child, and only when `min_separation` is set. This is confirmed by the 8 lines for 8 orbit points above.

**The loop split alone changes nothing.** With the exemption removed (break 1 below), the branch's restructured loop reproduces base exactly on the rubric-7 call:
- every log line is the same (`diff` of the full captured logs is empty);
- `too_close 76`, `registry 486`, depth 7.

So the move of the child loop after all parents are mapped has no effect of its own there. Rubrics 4 and 5 show the same on every law case and series.

## Rubric 9

Read in the diff:
- The `min_separation` docstring says the guard "never decides the depth of ``T(q0)``" and names the pip-orbit exemption.
- `BlastResult` documents `too_close`, `already_known` and `pip_depth`.

## Rubric 10: log lines (a finding the coder's report omits)

On every law case and series nothing changes (rubrics 4 and 5). Inside the rubric-7 call, however, the branch logs **WARNING lines the base does not log at all**. Base logs no TangleWorkbench record in that call. The branch logs:
```
WARNING 361 same-stability polyline crossing(s) of the new image bridges in the tangle region (A unstable 0.0 x image bridge (258, 259); ... (346, 347); (347, 348) ...): ... polygon artifacts ...
INFO    898 same-stability polyline crossing(s) of the new image bridges outside the tangle region
WARNING 1 of 494 detections are not crossings of N: ... 349 at (2.73205007717, -2.73204695504): ratio 1.775, ratio below 3 (clause 1)
```

**Why they appear.**
- The exempted `f^9(q0)` child `(258, 259)` and its images are among the same-stability hits. That child was exempted precisely because it runs within `min_separation` of existing curve, which is the re-tracing the guard exists to stop.
- The unresolved detection 349, ratio 1.775, is the (ii) stall point `f^{10}(q0)`. It has the same ratio as the series' nested_outer stall (plan: id 447, 1.775).

**My reading.** These are direct consequences of the intended section-5 change (the call now reaches depth 9), not a separate behaviour change. But rubric 10 says "the only new log lines are the exemption INFO lines", and within a multi-iteration call that is not literally true. The reviewer should decide whether section 5 covers it. The author may want to know that the exemption reintroduces near-coincident curve the guard used to stop: 361 polygon artifacts in the tangle region in this call.

## Break-it checks

Each break was applied in the branch worktree, run, and restored with `git checkout src/tanglepack/loom/Blast.py`. `git status` was clean after each one. The tests were `test_blast_pip_orbit.py` and `test_blast_proximity_guard.py`; the rubric-7 call is `brk.sh`.

| # | Break | Tests | Rubric-7 call |
|---|---|---|---|
| 1 | Remove the exemption: `pip_end = None` | **6(a) FAILS** (1 failed, 5 passed) | depth 7, `pip_depth {'A': 7}`, 0 INFO lines: **fails rubric 7** (needs 9). All log lines equal base. |
| 2 | Compose `pip_depth` wrongly: `iterate_table.forward_depth(orbit[0])` | **6(c) FAILS** (1 failed, 5 passed) | `pip_depth {'A': 1}` while the composed depth is 9: **fails rubric 7** |
| 3 | Exemption reads the direct table entry: orbit = `{table[q0, n]}`, which is None for n ≥ 2 (the table holds `+1` links only) | **6(a) FAILS** (1 failed, 5 passed) | depth 7, 0 INFO lines: children at `f^2(q0)` and deeper are no longer exempt, and the `f^7`/`f^8` children are dropped again: **fails rubric 7** |
| 4 | Series law fed a decreasing list (`break_law.py`): the real k28 series with its last `pip_depth` lowered by 1 | the law **FAILS**: `AssertionError: [1, 2, 3, 4, 5, 6, 7, 8, 7]` (the real series `[1..8, 8]` passes) | n/a. A law weakened to "accept a decrease" (`depths[-1] >= depths[0]`) would accept that list (`True`), so the `== sorted` form is what catches it. |
| 5 | Disable the guard: `if False and (too_near_existing or too_near_sibling)` | **6(b) FAILS** and **`test_min_separation_drops_close_bridges[no_pip]` FAILS** (2 failed, 4 passed) | n/a. The guard still demonstrably fires without a pip. |
| 6 | Exemption moved BEFORE the exterior check: a pip-bearing child skips the zone, fixed-point and already-iterated filter | **No existing test failed** (6 passed). Rubrics 4 and 5 are blind too: every law-case and series comparison was still equal (`te_b3/b6/`). The only trace was a 9th INFO line, `child (4, 5) ends at f^0(q0)`, in the rubric-7 call. With my new test: **`test_exemption_skips_only_the_separation_test` FAILS** (`zone.contains_point(array([15.20, 4.15]))` is False: an exterior pip-bearing child was kept). | 9 INFO lines instead of 8 |

**What break 6 shows.** The exemption's *position* (rubric 8, the first sub-item) was guarded only by reading the code. I added the test below to guard it.

## New test (committed, `7d3a294`, in `tests/unit/loom/test_blast_pip_orbit.py`)

`test_exemption_skips_only_the_separation_test(henon_p3_session)`:
- **Blast.** The same three-iteration, `fp1`, separation-0.02 blast as 6(a)-(c).
- **Assertions.** Every child in every step's `kept_interior`:
  - belongs to `fp1`;
  - has its polyline midpoint inside the zone, and outside every smaller stored zone (C4);
  - was never fed before, neither in the initial frontier nor as a kept child of an earlier step.
- **What it guards.** The exemption bypasses the separation test and nothing else (rubric 8, ordering).
- **Style.** No ids are pinned. It uses public API only: `polyline_midpoint` from `numerics.geometry` and `ResonanceZone.contains_point`.
- **Results:** passes on base (rubric 2 table), passes on the branch, fails under break 6.

## What I could not test, and why

- **The INFO line count per exempted child** is checked by measurement (8 lines for orbit points 2..9 in the rubric-7 call), not by a test. The 2026-10-05 test policy forbids pinning log wording.
- **The structural move of the child loop** was checked only on the law cases, the five series and the 11-iteration call (and break 1's split-only run equals base there).
  - I did not replay E12's two-iteration cells. The coder reports k28, k32 and p3 identical and nested_outer differing from call 3 on, which is the intended change.
  - **Not re-measured by me.**
- **`scripts/henon_blast_period_3.py`** raises before any blast on base and branch alike, according to the coder (stale pinned pip id 12). I did not re-run it, since it is outside the file list.
