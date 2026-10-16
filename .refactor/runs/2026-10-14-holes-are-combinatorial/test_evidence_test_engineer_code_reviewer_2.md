# Test evidence: B4b `holes-are-combinatorial` (test engineer, round 1)

- **Plan:** `artifacts/plan_planner_human_2.md` (rubric in §7, behaviour in §5).
- **Coder report:** `artifacts/change_coder_test_engineer.t1.md`.
- **Base:** `f166573`, in the detached worktree `/home/dezhu/Development/TanglePack/tangle-pack-holes-are-combinatorial-base`.
- **Branch:** `refactor/holes-are-combinatorial`, in `/home/dezhu/Development/TanglePack/tangle-pack-holes-are-combinatorial`. It holds the coder's five commits `391f952..145455c` plus my test commit **`ba23a1f`**.
- **Python:** `MPLBACKEND=Agg PYTHONPATH=src:tests /home/dezhu/Development/TanglePack/tangle-pack/env/bin/python`, run from inside each worktree. pytest is `MPLBACKEND=Agg ../tangle-pack/env/bin/python -m pytest`, run from the worktree.
- **Probe scripts** are mine, written for this check and not copied from the coder. They live in the session scratchpad `.../scratchpad/te_b4b/`:
  - `probe.py` and `compare.py`: rubrics 4, 7, 10 and 11.
  - `nomap.py` and `nomap_compare.py`: the no-map fallback.
  - `rubric8.py`: anchor openings. It names the branch itself from `unstable_key_of` and does not call `stable_key_at`, so it does not depend on the code under test.
  - `lawwarn.py`: rubric 7's silence check.
  - `points.py` and `points_compare.py`: the break-it check on the `map_inv` carry.
- After every break the source was restored with `git checkout <file>`, and `git status` was clean apart from my uncommitted test.

**Verdict:** every rubric item I re-ran holds. Every requested break was caught by a test, by the base-vs-branch comparison, or by both. I added one test, because the no-map fallback had no test guarding it.

## 1. Equality: base vs branch (rubrics 4, 7, 10, 11)

I ran `probe.py` in each worktree. The base read each point from `Hole.coords`; the branch read it from `plotting.hole_point`. Per case the probe dumped:

- the law `fingerprint(case)`, minus the process-local `products` ids, after `build_products`;
- the full `topology().dynamics.describe()` text (every class, word, refinement and matrix; 1.2k to 57k characters per case);
- each hole's five fields `(origin, iterate, bounding_ids, bridge_side, openings)`;
- each hole's point;
- `diagnostics.measured_sides`;
- every WARNING on the `tanglepack.topology.StablePartition` logger during the build.

`compare.py` maps FixedPoint reprs to `FP0, FP1, …` per case. A point passes only if `max|d| <= 1e-12`.

```
k10                  holes=  3 (direct 1, propagated 2) fp=True words=True five=True keys=True pts=True measured=True max|d| direct=0 propagated=0 measured_n=2 SPwarn base=0 branch=0 new_kind=0
k28_one_blast        holes=  4 (direct 2, propagated 2) fp=True words=True five=True keys=True pts=True measured=True max|d| direct=0 propagated=0 measured_n=2 SPwarn base=0 branch=0 new_kind=0
k28_two_blasts       holes=  4 (direct 1, propagated 3) fp=True words=True five=True keys=True pts=True measured=True max|d| direct=0 propagated=0 measured_n=3 SPwarn base=0 branch=0 new_kind=0
p3                   holes=  8 (direct 2, propagated 6) fp=True words=True five=True keys=True pts=True measured=True max|d| direct=0 propagated=0 measured_n=6 SPwarn base=0 branch=0 new_kind=0
nested               holes= 12 (direct 3, propagated 9) fp=True words=True five=True keys=True pts=True measured=True max|d| direct=0 propagated=0 measured_n=9 SPwarn base=0 branch=0 new_kind=0
nested_inner_blasts  holes= 12 (direct 3, propagated 9) fp=True words=True five=True keys=True pts=True measured=True max|d| direct=0 propagated=0 measured_n=9 SPwarn base=0 branch=0 new_kind=0
inversion            holes=  6 (direct 2, propagated 4) fp=True words=True five=True keys=True pts=True measured=True max|d| direct=0 propagated=0 measured_n=4 SPwarn base=0 branch=0 new_kind=0
inversion_one_blast  holes=  7 (direct 1, propagated 6) fp=True words=True five=True keys=True pts=True measured=True max|d| direct=0 propagated=0 measured_n=6 SPwarn base=0 branch=0 new_kind=0
inversion_two_blasts holes=  8 (direct 2, propagated 6) fp=True words=True five=True keys=True pts=True measured=True max|d| direct=0 propagated=0 measured_n=6 SPwarn base=0 branch=0 new_kind=0
inversion_period2    holes= 10 (direct 2, propagated 8) fp=True words=True five=True keys=True pts=True measured=True max|d| direct=0 propagated=0 measured_n=8 SPwarn base=0 branch=0 new_kind=0
heteroclinic         holes= 73 (direct 4, propagated 69) fp=True words=True five=True keys=True pts=True measured=True max|d| direct=0 propagated=0 measured_n=69 SPwarn base=0 branch=0 new_kind=0
p3_coarse_blasts     holes= 84 (direct 4, propagated 80) fp=True words=True five=True keys=True pts=True measured=True max|d| direct=0 propagated=0 measured_n=80 SPwarn base=0 branch=0 new_kind=0
ALL EQUAL
```

- **Rubric 4: holds.** The fingerprints are equal on all 11 law cases. That covers crossings, bridge ids, holes, pips, partition signature, alphabet and generations. The five fields and the words are also equal.
- **Rubric 10: holds.** On every hole of every law case, `k28_two_blasts` included, `hole_point` on the branch equals the base's `Hole.coords` exactly (deviation 0.0). This holds for direct and propagated holes alike: 22 direct and 135 propagated over the 11 cases. `p3_coarse_blasts` also deviates by 0.0.
- **Rubric 11: holds.** `measured_sides` is equal on the 11 law cases and on `p3_coarse_blasts`.
- **Rubric 7 (silence): holds.** Zero `StablePartition` WARNINGs of any kind on the 11 law cases, on base and branch. `lawwarn.py` captured the logger again on the restored branch: `law cases logging the NUM-5 WARNING: 0 []`.

**No-map placement** (`nomap.py`). This is not a rubric item, but it is the path the fourth break targets. On each law case I set `trellis.dynamical_system = None` and re-punched (`punch_holes(diagnose=False)`):

- **Base and branch agree exactly:** base `Hole.coords` against branch `hole_point` gives `max|nomap base-branch| = 0` on all 11 cases (`ALL EQUAL`).
- **The fallback really runs:** without the map, holes move by 0.0034 to 0.284 against their with-map points.

## 2. Static checks (rubrics 3, 5, 8, 9)

- **Rubric 3: holds.** `git diff --stat f166573..refactor/holes-are-combinatorial` touches exactly the nine plan files plus `TangleSession.py`.
  - At the coder's head the total was 10 files, +423/−403.
  - With my test it is +442/−403, and the files are the same.
  - The `TangleSession.py` diff is one docstring line: `` ``epsilon``, `` is dropped from the `punch_holes` kwargs list.
- **Rubric 5: holds.**
  - `wc -l` gives 1091 on the branch (1354 on base), which is within ≤ 1,100.
  - The rubric grep `grep -n "epsilon\|coords\|_step_into\|np\.\|map_inv\|dynamical_system\|_map_backward\|_region_key\|punched_regions\|fallbacks" src/tanglepack/topology/StablePartition.py` prints nothing (exit 1).
- **Rubric 8: holds.**
  - **Static:** `StablePartition.py:697` is `key = trellis.stable_key_at(iid, trellis.unstable_key_of(bridge.id))`. The base's inline rebuild `key = (fp, "stable", orbit_index, branch_index)` is gone; the identification at `:701-702` only flips `branch_index` of the key it got.
  - **Counts:** `rubric8.py` gives identical output on base and branch, including the full per-opening listing (`diff` empty):
    ```
    inversion: anchor ends 4, opened on s(o,1-b) 2, on s(o,b) 2
    inversion_one_blast: anchor ends 2, opened on s(o,1-b) 2, on s(o,b) 0
    inversion_two_blasts: anchor ends 2, opened on s(o,1-b) 2, on s(o,b) 0
    inversion_period2: anchor ends 8, opened on s(o,1-b) 4, on s(o,b) 4
    ```
    These are the plan's 2/4, 2/2, 2/2 and 4/8.
- **Rubric 9: holds.**
  - `inspect.signature(Trellis.punch_holes)` is `(self, pairs=None, *, propagate=True, diagnose=True, verbose=False)`.
  - The `Hole` fields are `[near_intersection_id, pair, bridge_side, bounding_ids, openings, iterate, origin]`; there is no `coords`.

## 3. Break-it checks (each restored with `git checkout`)

| # | Break, on the branch | What caught it |
|---|---|---|
| 1 | Deleted the `assert all(p.iterate >= 0 for p in pairs)` line in `punch_holes`. | `test_a_backward_pair_is_never_punched` **FAILED**. Without the assertion, the `iterate=-1` pair just fell through to the "Skipped 1 pair(s) with no spanning bridge" WARNING. |
| 2a | The NUM-5 WARNING now fires on table ends too (`… in ends if not from_table` became `… in ends if True`). | `test_scaled_span_end_beside_a_crossing_warns[table]` **FAILED**. `lawwarn.py` over the 11 law cases (build plus products plus topology): **6 law cases log it**, namely `k28_one_blast` (1), `p3` (3), `nested` (3), `nested_inner_blasts` (3), `inversion` (2) and `inversion_period2` (4). This is the plan's "6 of 11", so the probe's logger capture can see the WARNING. Restored: 0 law cases. |
| 2b | The WARNING never fires (`if near is not None and False:`). | `test_scaled_span_end_beside_a_crossing_warns[scaled]` **FAILED**. |
| 3 | Dropped one `map_inv` step in `hole_point` (`range(steps)` became `range(steps - 1)`). | **The 1e-12 comparison FAILED** on every law case. Holes beyond 1e-12: k10 2/3 (max 0.355), k28_one_blast 1/4, k28_two_blasts 2/4, p3 6/8, nested 8/12, nested_inner_blasts 8/12, inversion 3/6, inversion_one_blast 4/7, inversion_two_blasts 4/8, inversion_period2 2/10, heteroclinic 26/73. `test_henon_propagated_hole_point_is_on_its_side` also **FAILED**. Variant `range(0)` (no carry at all): the same test **FAILED**. |
| 4 | No-map fallback: `point = _step_into(midpoint, stable_mid, epsilon, arc_chord)` became `point = midpoint`. | **The no-map comparison FAILED** on every law case (`max|nomap base-branch|` from 0.0027 to 0.179, `DIFFERENCES`). The existing suite (`test_stable_partition.py`, `tests/plotting`: 117 passed) did **not** catch it, so I added `test_a_hole_is_drawn_off_its_bridge` (§4). With the same break it **FAILED** (`[no_map]`). |
| 5a | Anchor block without the identification (`key = (fp, stability, orbit_index, 1 - branch_index)` became `… branch_index)`). | **The inversion counts changed** to 0/4, 0/2, 0/2 and 0/8 on `s(o,1-b)`; everything opened on `s(o,b)` instead. `tests/invariants tests/unit/topology` gave **15 failed**: `test_partition_laws[*]` on all 11 law cases and four `test_partition_set_aside` items (k10 and k28_two_blasts × `anchor_opening_faces_its_hole` and `partition_singletons`). |
| 5b | Anchor block always reads `s(o,1-b)` (`if not outward:` became `if True:`). | **The counts changed** to 4/4 on `inversion` and 8/8 on `inversion_period2`. The two blast cases stay 2/2, because all their anchor ends were already anchorward. **7 failed**: `test_partition_laws[k28_one_blast, p3, nested, nested_inner_blasts, inversion, inversion_period2, heteroclinic]`. |

## 4. New test (commit `ba23a1f`, `tests/unit/topology/test_stable_partition.py`)

`test_a_hole_is_drawn_off_its_bridge[map | no_map]`:

- **Setup:** it uses the real `henon_with_holes` fixture (the conftest Hénon tangle with bridges, pips, pseudoneighbors and holes). The `no_map` parametrization sets the documented `trellis.dynamical_system = None`; nothing is mocked.
- **What it asserts:** for every hole, `hole_point` is finite and lies strictly off its own bridge's polyline. The drawing must commit to a side. It also asserts that propagated holes exist.
- **What it guards:** the no-map fallback (break 4 puts the point on a bridge node, distance 0) and, in the `map` case, the normal nudge.
- **Measured margins:** 0.016 to 1.39 on the fixture.
- **Pinned on base:** I put a scratch copy into the base worktree that reads `Hole.coords` and re-punches with `dynamical_system = None` for `no_map`. Both cases pass (`2 passed`). The copy was then deleted.
- **On the branch:** it passes (`2 passed`). Its teeth are break 4 above, and it also fails under break 3's `range(0)` variant through the sibling test.

The coder's two new tests:

- `test_a_backward_pair_is_never_punched` (rubric 6) and `test_scaled_span_end_beside_a_crossing_warns[scaled|table]` (rubric 7) pass on the branch, and breaks 1, 2a and 2b each fail one of them.
- Neither can run on base: the assertion and the WARNING do not exist there.

## 5. Full suite (rubrics 1, 2)

- **Branch** at `ba23a1f`, `MPLBACKEND=Agg env/bin/python -m pytest -q -rxX`: **1178 passed, 38 skipped, 5 deselected, 1 xfailed** in 211 s. That is the coder's 1176 plus my two parametrizations.
- **Base** at `f166573`, same command: **1173 passed, 38 skipped, 5 deselected, 1 xfailed** in 211 s. This is the full tier set, so every hole-touching test (law tier, `test_stable_partition`, `test_topology_invariants`, `test_minimal_trellis`, `test_production_checks`, `tests/plotting`) passes on base too.
- **The only xfail on both** is `test_orientation_reversing_case_builds[orientation_reversing]` (C3: INV-11), a `cases.KNOWN_ISSUES` entry.
- **No existing test was weakened or deleted.** The coder renamed `test_henon_propagated_hole_side_matches_coords` to `..._point_is_on_its_side`. Its assertion is unchanged; it reads `hole_point` where it read `coords`.

## 6. Not tested, and remarks

- **Rubric 7's law-case silence is probe evidence, not a committed test.** The law tier builds each case once per session behind `law_case`, so no test captures build-time logs. Break 2a shows that the probe's capture would see the WARNING.
- **The NUM-5 WARNING on `k28_deep_blast_sliver`** and the disappearance of `near_tangency`'s 14 "cannot propagate" WARNINGs (§5 logs): I did not re-run these; they rest on the coder's `b4b_extra` output. Neither is a rubric item.
- **No test covers non-default `epsilon`** in `hole_point`. The coder notes that `epsilon` now also scales the propagated nudge, which used a literal 0.05 before. At the default the two agree, and nothing passes another value.
- **Two remarks on the coder's test code:**
  - `test_scaled_span_end_beside_a_crossing_warns` replaces `trellis.image_cdist` with a lambda on the instance. This is a stub, and I consider it justified: no cheap real build places a scaled end within `0.5 * SCALING_RTOL` of a crossing.
  - The block before `test_a_backward_pair_is_never_punched` opens with three blank lines. That is cosmetic, and I left it as it is.
