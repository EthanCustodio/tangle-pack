# B1 `resolved-crossings` — change, review round 2 (coder)

Job: fix rubric 1 and 2 only, inside the 4 files the author added on 2026-10-12 ("Rule stands; widen
scope to fix the four tests"). I did not change the rule or any `src/` file.

- **Worktree / branch:** `/home/dezhu/Development/TanglePack/tangle-pack-resolved-crossings`, `refactor/resolved-crossings`
- **Base:** `1029fcc` (unchanged). Round 2 starts at the test engineer's `dc513f7`. HEAD is now `b1a7e9c`.

## Commits this round

| commit | what |
|---|---|
| `cb5cf92` | Laws: the pip-orbit-cut oracle reads N. An unresolved detection is in no trellis. |
| `e044f3d` | Tests: the p2inv10 walk test refines to 1e-6, so its reference window is in N. `build_inversion_period2` gains an `area_cutoff` keyword (default 1e-5, the law case). |
| `e8c8105` | Tests: the coarse-blasts build turns clause 1 of N off, so its 2-4 node hole bridges still reach the partition. |
| `f17a20f` | Tests: the order-flip build blasts with clause 1 of N off, so T(q0) keeps the pip orbit's ratio-2.99 point. |
| `b1a7e9c` | Cases: the `inversion_period2_10` deep run refines to 1e-6, so its reference window is in N and no law runs vacuous. |

`git diff --stat 1029fcc` (18 paths: the 14 planned plus the 4 the author added):

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
 tests/cases.py                                     |  56 +++-
 tests/helpers/laws.py                              |  14 +-
 .../regression/test_class_count_under_blasting.py  |  34 ++-
 .../regression/test_coarse_blasts_rows_by_signs.py |   5 +-
 tests/regression/test_evidence_diagnosis.py        |  18 +-
 tests/regression/test_pip_orbit_cut.py             |  15 +-
 .../unit/numerics/test_noise_crossing_collapse.py  |  81 ------
 tests/unit/numerics/test_resolved_crossings.py     | 213 +++++++++++++++
 18 files changed, 920 insertions(+), 294 deletions(-)
```

This round touched `tests/cases.py` too. It was already in the list, and the deep-run build and two of the reworked builds live there.

## Suite

`MPLBACKEND=Agg env/bin/python -m pytest -q -rxX` on `b1a7e9c` (the worktree's `pyproject` puts its own `src` first):

```
XFAIL tests/invariants/test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing] - C3: INV-11: ...
1154 passed, 38 skipped, 5 deselected, 1 xfailed in 197.89s (0:03:17)
```

The only xfail is the `KNOWN_ISSUES` C3/INV-11 entry. Before this round: 4 failed, 1150 passed.

Fixing the oracle exposed a fifth failure that the oracle's assertion had been hiding: `test_open_p3_deep_runs[inversion_period2_10]`, where the law `iterated_cut_provenance` "checked nothing". It has the same cause as (d): at 1e-5, crossings 38 and 39 are out of N and the reference window moves. The build lives in `tests/cases.py` (`BUILDERS["inversion_period2_10"]`), so I fixed it there with the same remedy (`b1a7e9c`). `test_open_p3_deep_runs.py` itself is not edited.

## The four failures: what each test guards now and how I checked it still has teeth

**(c) `laws.check_trellis_is_cut_at_the_pip_orbit`** (`cb5cf92`, one condition plus the docstring).
- The loop over `trellis.registry` now skips `iid not in trellis.resolved_ids`. The law reads "the trellis holds exactly the crossings of T(q0) in N".
- Teeth:
  - On the 1e-5 `inversion_period2` 10/10 build the law checks 86 crossings and passes.
  - Mutant: wrap `Trellis.cut_at_pip_orbit` so that, after cutting, it puts one resolved crossing from past the cut back on its branches. The law then fails with `(7, True, 1.488…, 0.787…)`.
  - Without the new condition the law fails on crossing 38, which is the reviewer's failure.

**(d) `test_pip_orbit_cut.py`** (`e044f3d`). Remedy: refine, as the author suggested.
- `_p2inv10` now builds at `area_cutoff=1e-6`.
- Scan of `area_cutoff` at 10/10 steps:

  | area_cutoff | unresolved | reference pair | classes | walks |
  |---|---|---|---|---|
  | 1e-5 | 2 | (30,31),(31,32) | 6 | 0 |
  | 5e-6 | 2 | (30,31),(31,32) | 6 | 0 |
  | 3e-6 | 0 (min ratio 3.49) | (38,39) | 17 | 12 |
  | 1e-6 | 0 (min ratio 5.43) | (38,39) | 17 | 12 |

  At 1e-6 the reference pair is the base's (38,39) again, and the words are reliable. I chose 1e-6 for the margin above 3.
- `test_p2inv10_every_class_walks` gains a precondition with a message: the build has no unresolved detection. Its two assertions are unchanged.
- `test_p2inv10_cut_equals_a_trim` uses the same helper and still passes.
- Teeth: with `Trellis.cut_at_pip_orbit` made a no-op, `test_p2inv10_every_class_walks` fails at `check_every_active_class_walks(case) > 0`. This is weaker evidence than reproducing the original pre-D10 defect: without the cut the dynamics ends up with no classes at all. I did not reconstruct the pre-D10 code.

**(a) `test_coarse_blasts_rows_by_signs.py`** (`e8c8105`). Remedy: turn clause 1 of N off for the build. Its assertions are unchanged.
- Refining cannot help here. This defect (NUM-1: rows and openings must come from the signs, not from a 2-4 node polyline) only appears on under-resolved bridges.
- I scanned `area_cutoff` ∈ {3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6} × 4-8 blasts under the default threshold. Not one hole bridge had alternating signs together with a geometric disagreement. On the base the ends of such bridges sit at ratios 0.65-1.05.
- So `build_period3_coarse_blasts`, whose only user is this test, sets `session.workbench.resolution_ratio_threshold = None` before its blasts. Clause 2 still applies.
- Result: 84 holes, and the hits are the hole bridges (42,41) and (43,42). The test is back to exercising the partition's sign rule on exactly those bridges.
- The builder and test docstrings now say why, citing Q13 (a) and B1.
- Teeth: the precondition `any(signs_alternate and geometry_disagrees)` fails under the default threshold, which is the reviewer's failure (a). The rest of the test is untouched.

**(b) `test_evidence_diagnosis.py::test_an_order_flip_gives_no_verdict`** (`f17a20f`). Remedy: turn clause 1 off for the blasts. Its assertions are unchanged.
- The flip itself is not an unresolved detection. I measured it on the base: the end, the bound and their images (247, 4, 419, 6) have ratios 7476, 7929, 14336 and 5783.
- It is lost only because f^9 of the default pip (ratio 2.991) is out of N, so T(q0) is shortened past the copy.
- I tried refining. None of these gave an order flip:
  - all blasts at 7e-8, 5e-8, 3e-8, 2e-8 or 1e-8;
  - growth and blasts at 3e-8;
  - only the last blast refined;
  - 9-12 blasts from either seed.

  Refining at 5e-8 or 3e-8 brings back the 23-class trellis, but it also moves the noise copy off the hole-bridge end.
- With `threshold = 2.5` there is no flip. With `threshold = None` (clause 2 still on) there is.
- So the test builds `build_k28(blasts=0, through="zones")`, sets the threshold to None, and does the 9 single-iteration blasts itself with `ZONE_MIN_SEPARATION`.
- The flip comes back on (243, 4, 407, 6), with ratios 7476, 7929, 14336 and 5783. Its precondition and assertions are unchanged.
- The trellis of that build holds 4 detections with ratio below 3 (148, 245, 246, 411). I wrote this in the module docstring's reasoning and report it here.
- Teeth: the precondition `flipped` fails under the default threshold, which is the reviewer's failure (b).

**For the reviewer and the author: tests (a) and (b) now opt out of clause 1.** I judged this to be reworking the build, not weakening an assertion, for two reasons:
- `resolution_ratio_threshold = None` is the documented public knob ("None turns clause 1 off, and clause 2 still applies").
- Both guarded defects live downstream of N: the partition's sign rule, and the evidence verdict.

Still, under the default threshold neither defect is reachable by any build I found. Say so if you would rather see these two tests marked NOT_APPLICABLE-style, or the defects re-guarded at unit level.

## Optional `Trellis.py:204` default (not done)
It cannot be done as a one-liner inside `Trellis.py`:
- `Trellis.__init__` has no threshold to compute N from.
- The safe fix is for `MinimalTrellis.py:411` to pass `resolved_ids=parent.resolved_ids` to the sparse snapshot, and that file is outside the list.
- It stays harmless today, because no reader of N runs on the sparse snapshot (the reviewer's reading).

## Behaviour changes this round
- None in `src/`.
- Tests:
  - `build_inversion_period2` gains `area_cutoff` (default 1e-5, so the law case is unchanged).
  - `BUILDERS["inversion_period2_10"]` and the p2inv10 regression build at 1e-6.
  - `build_period3_coarse_blasts` and the order-flip build blast with clause 1 off.
- No law case changes, so the fingerprint guard passes.
- The deep-run numbers the author listed (k28 x8 24/15/23 -> 4/4/4; k32 x9-10 holes 58->56, closed 39->45, reliable False->True) are still produced and are not pinned anywhere.

## Rubric, my reading
1. **Holds as amended by the author.** The full suite is green. The only tests edited outside the 14 are the 4 the author named. `tests/cases.py` was already in the list.
2. **Holds.** 1154 passed, 38 skipped, 5 deselected, 1 xfailed (KNOWN_ISSUES C3/INV-11 only).
3. **Holds.** `git diff --stat 1029fcc` shows exactly the 18 paths.
4-13. **Unchanged from the reviewer's pass.** Nothing under `src/` moved this round, and the law-case builds and fingerprints are untouched (`build_inversion_period2`'s default is still 1e-5).
