# Change — B2 `pseudoneighbor-both-intervals` (coder, iteration 1)

*2026-10-14. Job: `change`. Plan `artifacts/plan_planner_human_2.md`, approved as written (`plan_ok_planner_human.json`; e13's jsonl committed as the fourth file).*

## Where

| | |
|---|---|
| Base commit | `ba23a1f1c03e6dd6189fb109ad7a75fcd3070f94` (`main`, the author's HEAD) |
| Branch | `refactor/pseudoneighbor-both-intervals` |
| Worktree | `/home/dezhu/Development/TanglePack/tangle-pack-pseudoneighbor-both-intervals` |
| Base worktree (detached, `ba23a1f`) | `/home/dezhu/Development/TanglePack/tangle-pack-pseudoneighbor-both-intervals-base` |

The author's main tree was not touched: its uncommitted `docs/`, `.refactor/` and `artifacts/` are as they were, and HEAD is still `ba23a1f`. Note for anyone running scripts: the env's editable install points at the main repo's `src`. Pytest uses the worktree's own `src` (`pythonpath = ["src", "tests"]`). Scripts need `sys.path`/`PYTHONPATH` set to the worktree's `src`; the e13 harness does this itself.

## Commits (`git log --oneline ba23a1f..refactor/pseudoneighbor-both-intervals`)

```
dfdfb2d Stress: e13 results, 14 cells (7 cases x 1e-7 / 1e-5) at the E12 depths: the stable clause loses no reference on any step (B2)
530c1ba Stress: e13 re-measures the E12 series with both intervals tested, before / after per blast and every lost pair with its reading and witness (B2 step 5)
9a3e2fe Tests: the stable interval is punctured by a scaled iterate and by a table landing, never by the pair's own orbit (B2 step 4)
9a83036 Pseudoneighbor: Dev Notes state the definition once, both intervals, and drop 'empty by construction' (B2 step 3)
06fef3f Pseudoneighbor: one _interval_is_clear tests both open intervals of a window pair (B2 steps 1-2)
```

## `git diff --stat ba23a1f`

```
 src/tanglepack/topology/Pseudoneighbor.py        | 238 ++++++++++++-----------
 stress/numerics/e13_both_intervals.py            | 198 +++++++++++++++++++
 stress/numerics/results/e13_both_intervals.jsonl |  14 ++
 tests/unit/topology/test_pseudoneighbor.py       | 108 +++++++++-
 4 files changed, 438 insertions(+), 120 deletions(-)
```

These are exactly the plan's four files. `Pseudoneighbor.py` is +128/−110, so it is **not** net shorter as the plan estimated. The function body is about the same length as before; the Dev Notes and the two docstrings grew because they now describe both intervals and the mirrored D7 guard.

## Suite

`MPLBACKEND=Agg env/bin/python -m pytest -q -rxX` on the branch:

```
XFAIL tests/invariants/test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing] - C3: INV-11: det J < 0 is out of scope (C3); FixedPoint.set_k_value rejects a saddle with exactly one negative eigenvalue
1181 passed, 38 skipped, 5 deselected, 1 xfailed in 209.21s (0:03:29)
```

Other targeted runs:
- `tests/unit/topology/test_pseudoneighbor.py` on the branch: 34 passed.
- `tests/invariants/test_law_partition.py`: 11 passed, 4 skipped on base and on the branch.

## What changed

**Step 1: `_interval_is_clear(trellis, x0, x1, stability, candidates, collision_rtol)`** replaces `_unstable_interval_is_clear`. It is one body. It reads the arc side through `attrgetter`: `arc_key` / `other_key` = `manifold_a_key` / `manifold_b_key` (swapped for `"stable"`), and `cdist` = `{stability}_cdist`.
- **Stable side:** it asserts that both ends share the stable key. **Unstable side:** it keeps today's `return False` when the unstable keys differ.
- `beta = owner.per_step_beta(stability)`, so on the stable side `beta < 1`.
- The m-range is `range(floor(min(bounds)), ceil(max(bounds)) + 1)`.
- The anchor skip is `cdist(r) <= 0`, and `lo <= 0` rejects on both sides.
- The candidates are unchanged: the caller's `trellis.own_intersection_ids`, gated by `pos_map.get(arc_key(r))`.
- The 2-D guard is `same_other = other_key(r)[0] is other_key(ix0)[0]`. For `"unstable"` this is today's `same_stable`; for `"stable"` it is D7 mirrored.
- `same_point` always receives `(unstable, stable)`. A local `point(value, other_value)` puts them back in that order on the stable side.
- The three readings run in today's order: table by id, own frame, scaled.
- One debug line names everything a reader or e13 needs: `"Pair (%d, %d): %s interval punctured by %d at M^%d, %s reading (cdist ... in (..., ...))"`, with args `(x0, x1, stability, r_id, steps, reading, value, lo, hi)`. This replaces the three old debug wordings. No test pinned them.

**Step 2:** `compute_pseudoneighbors` keeps a pair iff `all(_interval_is_clear(..., s, ...) for s in ("unstable", "stable"))`. Its docstring now names both intervals.

**Step 3, Dev Notes:**
- The head paragraph states the definition once ("NEITHER open interval ... holds a point of X"; both are tested by one function, B2/Q1).
- The "empty by construction" bullet is deleted.
- The β note now says `per_step_beta(stability)` and never `1/per_step_beta("unstable")`.
- The closed-form note now says β^k per lap, with min/max because the bounds swap when β < 1.
- The check-set note now covers both intervals, gives the mirrored D7 guard, and says why the candidates are never the whole registry (law case `inversion`, reference (2, 3)).
- The heteroclinic bullet now names B's cycle for the stable arc. The lookup bullet now says "its cdist" and "r's own branch".
- The orphan-pair note stays. No NUM-18 paragraph existed.

**Step 4, three tests**, each confirmed to be decided by the reading it names (debug log):
- `test_scaled_iterate_punctures_stable_interval`: `(u=16, s=0.5)` below the window, no links. Its scaled stable landing is 2, inside `(1, 3)` of `(x1, r_end)`, a pair with no registered images. The unstable orbit only touches the interval ends. Expected `[(r_n, x1)]`. *My first draft targeted `(r_n, x1)`; that was decided by the own-frame reading, because `M(r_n) = r_end` is linked, so I retargeted it.*
- `test_table_landing_punctures_stable_interval`: `w` at `(u=8, s=3.5)` and `r = M(w)` at `(32, 0.875)`, linked. The trellis is rebuilt with `w` off every branch (`w not in trellis`, `r in trellis`, as past a cut). The table reading kills `(r_n, x1)` by id. Expected `[(x1, r_end)]`.
- `test_own_orbit_does_not_puncture_stable_interval`: an unlinked `(8.008, 0.7495)` scales back to `(2.002, 2.998)`, inside `(1, 3)` and 2-D-colliding with x1. It imports `_interval_is_clear` inside the body and asserts the stable call is True and the pair is kept. The guard is load-bearing: at `collision_rtol=0` the same call returns False.

**Step 5, e13** (`stress/numerics/e13_both_intervals.py`). It is a harness module like E12 and imports E12's `DEPTH`, `SEPARATION` and `MAX_POINTS`, so the depths cannot drift. There are 14 cells: 7 cases × {1e-7, 1e-5}, one iteration per call, each built by `e11.build`.

Per blast, one `session.topology()` runs, and then every consecutive pair of the reference window is judged by the one function, called twice. The record holds:
- `before` (unstable call clear) and `after` (both calls clear);
- `lost`: each pair with `reading`, `witness`, `steps`, `value`, `interval` and `reference_witness(trellis, pair)`;
- `n_stable_fires`: on every window pair, not only the references;
- `matches_library`: `after == trellis.reference_pseudoneighbors`;
- `n_classes` (seed), `n_closed`, `n_holes`, `n_trellis_crossings`, and `topology_error`.

The reading and witness come from the debug line's `record.args`, caught by a handler on the Pseudoneighbor logger. A printed summary follows the sweep. I checked the loss path on the scaled-puncture fixture: it reports `{'pair': (1, 2), 'reading': 'scaled', 'witness': 3, 'steps': -1, ..., 'reference_witness': None}`. `--quick` (k28 × 3 and p3 × 2 at 1e-7) runs in about 2 s. Its `_quick.jsonl` is **not** committed, because the plan admits only the full results file.

## Behaviour changes

**1. Seven existing unit-test fixtures. This is a deviation the reviewer must weigh against rubric 2.** With the stable clause in place, 7 tests of the base `test_pseudoneighbor.py` failed on the branch:
- `test_direct_intersection_punctures_interval`
- `test_iterate_punctures_interval`
- `test_walk_stops_at_reference_window_end`
- `test_heteroclinic_candidate_ignored`
- `test_table_linked_deep_iterate_does_not_disqualify`
- `test_table_linked_landing_uses_registered_cdist`
- `test_unlinked_landing_is_read_against_the_pairs_images[below_c]`

Each failure is a real puncture by the definition. The fixture adds a crossing on the window's stable branch whose stable orbit visits the window *strictly inside* a pair's stable interval. For example, `(u=1.5, s=0.5)` one step back sits at `s=2`, inside `(1, 3)`. In a period-1 window that is exactly one lap long, every point of X visits the window once. An extra crossing that is not an iterate of a window crossing therefore always punctures some stable interval. These fixtures were written when only the unstable side was read. The analyst measured real builds only, so the plan did not foresee this.

I fixed the **data only**. No assertion and no test name changed:
- Crossings that probe the unstable interval alone moved to the saddle's other stable branch, `OTHER_STABLE = 1`, so their stable orbit never visits the window. The base reading compares only the fixed point and the cdists, so base behaviour is unchanged. This covers direct, iterate, heteroclinic, table-linked-landing (its `z → z1 → r` chain) and unlinked-landing.
- `walk_stops`: `below` moved to `M(r_end)`'s cdists, `s` 0.5 → 0.25.
- `deep_iterate`: the drift of the sixth image's stable cdist flipped sign, from `1 + 2·rtol` to `1 − 2·rtol`. It still misses the collision, and it now lands inside `(x1, r_end)`, which that landing already punctures.
- The module docstring explains the rule.

**The adjusted file passes unchanged on base src (31 of 31) and on the branch.** So rubric 2 holds for the test file as it now stands on the branch. It does **not** hold for the literal base version of those 7 tests, which cannot pass on any correct implementation of the definition.

**2. Real builds: none.**
- **Rubric 4:** fingerprints, reference pairs, all pseudoneighbor pairs and class words are identical on base and branch for all 11 law cases (details below).
- **e13:** no reference lost on any of the 14 cells at any step.

**3. Debug wording.** The three puncture debug lines are now one format (above). This is DEBUG only, and nothing pinned the old wording.

## Rubric 4 evidence

Script: `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b2_rubric4.py`. It is run from each worktree root and imports that worktree's `src`; I checked this by printing `Pseudoneighbor.__file__` and `hasattr(_interval_is_clear)`. Per case in `cases.LAW_CASES` it writes:
- `law_tier.fingerprint(case)`, with `products` dropped and memory addresses stripped. `products` holds per-process `id()`s and cannot compare across processes; every other aspect is kept.
- sorted `compute_pseudoneighbors(topology.trellis)`;
- every pseudoneighbor pair with its iterate;
- every class's `(label, letter, word)`.

Outputs are `b2_r4_base.jsonl` and `b2_r4_branch.jsonl` (11 lines each). **`diff` is empty** (`b2_rubric4.diff`, 0 bytes, exit 0). Reference pairs per case (identical on both):

| case | refs |
|---|---|
| k10 | (5,6) |
| k28_one_blast | (3,4) (4,5) |
| k28_two_blasts | (9,10) |
| p3 | (3,4) (4,5) |
| nested | (7,8) (8,9) (37,38) |
| nested_inner_blasts | (7,8) (8,9) (41,42) |
| inversion | (1,2) (2,3) |
| inversion_one_blast | (26,27) |
| inversion_two_blasts | (23,24) (26,27) |
| inversion_period2 | (13,14) (14,15) |
| heteroclinic | (77,78) (78,79) (92,93) (113,114) |

## e13: what it found

`env/bin/python stress/numerics/e13_both_intervals.py --workers 8`: 14/14 `ok`, every cell to its full depth (`stop_reason = max_blasts`), with no topology errors. The slowest cell was k32 at 1e-7, 100 s and 5.0 GB.

**The stable clause rejected no window pair on any step of any cell.** `n_stable_fires = 0` holds on every window pair, references or not. So `before == after` everywhere, nothing is lost, and `after` equals the library's references on every step. This extends the analyst's result to the cells the analyst did not measure: the 1e-5 cutoff, `nested_inner` and `heteroclinic` blasting. **For DEC-25 this gives no new input.** The inversion series drops from 2 to 1 references at blast 6 at both cutoffs, and the stable clause plays no part in that drop.

References per blast (`after`), and N:

| cell | references per blast | N stops growing at |
|---|---|---|
| k28 1e-7 | 2,2,1,1,… (1 from blast 2) | 143 from blast 7 |
| k28 1e-5 | 2,2,1,… | 86 from blast 6 |
| k32 1e-7 | 2×8, then 6 from blast 8 | 300 from blast 8 |
| k32 1e-5 | 2×7, then 3 from blast 7 | 173 from blast 7 |
| p3 1e-7 | 2 throughout | 119 from blast 14 |
| p3 1e-5 | 2 throughout | 47 from blast 7 |
| nested_outer 1e-7 | 4,4,3,… | 276 from blast 8 |
| nested_outer 1e-5 | 4,4,3,… | 116 from blast 6 |
| nested_inner 1e-7 | 3 throughout | still growing at 12 (105) |
| nested_inner 1e-5 | 3 throughout | 57 from blast 7 |
| inversion 1e-7 | 2,1,2,2,2,2,**1**,1,1,1,1 | still growing at 10 (212) |
| inversion 1e-5 | 2,1,2,2,2,2,**1**,1,1,1,1 | 167 from blast 9 |
| heteroclinic 1e-7 | 4 throughout | 197 from blast 5 |
| heteroclinic 1e-5 | 4 throughout | 142 from blast 2 |

Seed / closed / hole counts are in the jsonl per step. For example, p3 at 1e-7 goes 5/5/8 → 23/43/42 at blast 6, and k32 at 1e-7 goes 3/3/4 → 34/45/56 at blast 8. Past each N stall the later blasts measure the same N again (B3). E12's committed results show larger N on some cells (k28 1e-7: 422). They predate B1 (N = resolved crossings), so they are not comparable; e13 is the current baseline.

## Rubric, my reading

1. **Pass.** Suite: 1181 passed, 38 skipped, 5 deselected, 1 xfailed. `-rxX` lists only `orientation_reversing_case_builds`.
2. **Pass for the branch's test file, with the deviation above.** The adjusted `test_pseudoneighbor.py` (minus the three new tests) passes on base and branch. `test_law_partition.py` gives 11 passed, 4 skipped on both. The literal base versions of 7 fixture tests fail on the branch, because their data puts points of X inside window stable intervals. The fixtures were corrected in data only. The reviewer must judge whether this is acceptable or whether the author should rule.
3. **Pass.** `git diff --stat ba23a1f..refactor/pseudoneighbor-both-intervals` lists exactly the four planned files.
4. **Pass.** The base and branch comparison is identical on all 11 law cases (fingerprint minus `products`, reference sets, all pairs, words). `products` cannot be compared across processes; it is the one aspect left out.
5. **Pass.** `grep -rn "_unstable_interval_is_clear" src tests stress` finds nothing. `_interval_is_clear` is defined once (`Pseudoneighbor.py:343`) with a `stability` parameter, and the caller calls it for `"unstable"` and `"stable"` (`:229`).
6. **Pass (hand check, not committed).** A stable interval `(0.5, 40)` (80×, three laps of 4) with a candidate at `s=0.01`, whose scaled landings include 0.64, returns **False**. With no candidate it returns True.
7. **Pass.** `beta = owner.per_step_beta(stability)`. The candidates remain the caller's `trellis.own_intersection_ids`, and the diff's only mention of `registry` is the Dev Note saying the candidates are never the whole registry. `resolved_ids` does not appear.
8. **Pass.** On the branch, all three pass. Against base src: test 1 fails `assert [(0, 1), (1, 2)] == [(0, 1)]`, test 2 fails `assert [(0, 1), (1, 2)] == [(1, 2)]`, and test 3 fails with `ImportError: cannot import name '_interval_is_clear'`. The test file collects on base: 31 others pass.
9. **Pass.** `grep -n "empty by construction"` finds nothing. The `compute_pseudoneighbors` docstring names W^U(x0, x1) and W^S(x0, x1).
10. **Pass.** `--quick` runs (2 cells, about 2 s). Every record holds, per step: `before`, `after`, `lost` (reading, witness, steps, `reference_witness`), `n_classes`, `n_closed`, `n_holes` and `n_trellis_crossings`. The full run covers the 14 cells at E12's depths. Its file is committed as `stress/numerics/results/e13_both_intervals.jsonl`.

## Files

- `/home/dezhu/Development/TanglePack/tangle-pack-pseudoneighbor-both-intervals/src/tanglepack/topology/Pseudoneighbor.py`
- `/home/dezhu/Development/TanglePack/tangle-pack-pseudoneighbor-both-intervals/tests/unit/topology/test_pseudoneighbor.py`
- `/home/dezhu/Development/TanglePack/tangle-pack-pseudoneighbor-both-intervals/stress/numerics/e13_both_intervals.py`
- `/home/dezhu/Development/TanglePack/tangle-pack-pseudoneighbor-both-intervals/stress/numerics/results/e13_both_intervals.jsonl`
- Scratch evidence (not committed), under `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/`: `b2_rubric4.py`, `b2_r4_base.jsonl`, `b2_r4_branch.jsonl`, `b2_rubric4.diff` (empty), `b2_suite_branch.txt`, `b2_e13_full.txt`, `b2_e13_digest.txt`.
