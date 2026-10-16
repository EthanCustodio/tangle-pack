# Outcome: B1 `resolved-crossings`

*Auditor (analyst, job `outcome`), 2026-10-12. The branch `refactor/resolved-crossings` is at `76d0710`, in the worktree `../tangle-pack-resolved-crossings`. Its base is `1029fcc`. Review passed in round 2, and there is no debrief verdict yet.*

**Tree state.**
- The branch worktree is clean.
- The main tree at `1029fcc` is dirty with uncommitted `docs/`, `.refactor/` and `artifacts/` files, as it was before the run. The branch does not touch them.
- Nothing is merged. Merging is your call.

## 1. Which detections are crossings of the trellis, now

A detection is in `N` iff `IntersectionRegistry.resolved_ids(threshold)` holds it (`IntersectionRegistry.py:441`). That is true when:
- it is synthetic; or
- **clause 1:** its `resolution_ratio` is at least `resolution_ratio_threshold` (3), **and**
- **clause 2:** it is in no run. A run is two neighbouring detections of the same `(u_key, s_key)`, on the same or a node-sharing unstable segment, that sit no farther apart along the polyline than the LARGER of their two chord deviations. A run is out whole, median included.

Every detection is still registered. The registry decides `N` over all its detections at once, so a run that straddles two `iterate_bridge` calls is still seen.

Seven readers now take `N` instead of the registry:
- `create_bridges` (`TangleWorkbench.py:1803`)
- the blast cut and its inference sources (`BridgeIterator.py:378`)
- `infer_iterates`, both sources and candidates (`IterateInference.py:109`, `:226`)
- `Trellis.from_workbench` (`Trellis.py:422`)
- `is_strong_pip` (`StrongPip.py:169`)
- `trim_stable_manifolds` (`TangleWorkbench.py:2178`)
- `_frozen_ids("all")` (`TangleWorkbench.py:1102`)

So an unresolved detection is never cut at, never linked and never trimmed at. It lies inside a bridge, and one WARNING names it with the remedy "refine". The blast-path vote is gone, and a plain recompute now leaves `N` and the bridges unchanged on all five series. That last point is *verified*: `test_recompute_keeps_the_crossing_set` passes on my run, and the `--k28x8` probe recompute gives 244 → 244, |N| 239 → 239, bridges 168 → 168.

**The one consequence you should weigh before merging:** every deep-run change below has a single cause. **The pip orbit's deepest detected point is itself unresolved.** It is never linked, so the pip-orbit cut stops one iterate earlier and `T(q0)` stays there under any further blasting until you refine (§7, point 1).

## 2. The most important hunk

Before, in `Tangle.populate_intersections_for_manifold`. The vote ran on the blast path only:
```python
        return self._collapse_noise_crossings(new_intersections)
    ...
            if len(run) % 2 == 1:
                keep = run[len(run) // 2]      # odd run -> keep the median
                kept.append(keep)
            else:
                keep = None                     # even run -> drop all
            ...
            for ix in run:
                if ix is not keep:
                    self._purge_crossing(ix)
```
After, the same method returns every detection, sorted by `(unstable_cdist, stable_cdist)` so ids stay reproducible (INT-7). The rule lives in the registry and every path reads it:
```python
    def resolved_ids(self, threshold: Optional[float]) -> frozenset[int]:
        ...
            runs = self._runs()[1]
            self._resolved_cache[threshold] = frozenset(
                iid
                for iid, ix in self._store.items()
                if iid not in runs
                and (threshold is None or ix.resolution_ratio is None or ix.resolution_ratio >= threshold)
            )
```
`_runs()` (`IntersectionRegistry.py:489`) compares each consecutive pair in a key group by `_polyline_gap <= max(a.chord_deviation, b.chord_deviation)`. That is the symmetric pair rule. The chord deviation is a new field, `Intersection.chord_deviation`, stored by `Tangle._resolution_ratio` (`Tangle.py:343`).

## 3. What got shorter or clearer

- **`Tangle.py`:** 1115 → 1052 lines (−63). `_collapse_noise_crossings` and `_purge_crossing` are deleted, and `populate_intersections_for_manifold` is down to 9 lines. The grep for `_collapse_noise_crossings|_purge_crossing|Collapsed` over `src tests` is empty.
- **`IntersectionRegistry.py`:** 860 → 968 (+108). This is where the rule now lives: `run_separation`, `run_partner`, `resolved_ids`, `_runs`, `_polyline_gap`, plus `unresolved_crossings` as its complement.
- **`TangleWorkbench.py`:** 2239 → 2260 (+21). Most of it is the threshold docstring ("Derived bound 1", "margin 3", "definition of N", corpus minimum 3.038) and the WARNING, which now names each detection.
- **`src` overall:** +295 / −193 across 9 files.
- **Clearer:**
  - One definition of the crossing set replaces one path's vote. The question "is this a crossing?" is now answered in one place.
  - `threshold=None` no longer has three special cases that silenced the report.
- **Tests:**
  - `test_noise_crossing_collapse.py` (−81) is replaced by `test_resolved_crossings.py` (+236). It covers the zig-zag, the symmetric pair, a clean crossing, the two `None` tests, and the N = registry pin on the law cases.
  - New probe: `tests/_tools/resolution_probe.py` (+290, not collected).

## 4. Every behaviour change

All numbers in this section come from my own runs. Script `scratchpad/audit/series.py` runs each series' D1 loop without a recompute between blasts, in both worktrees.

**Law cases (11): none.**
- The branch probe (`--laws`) finds clause 1 and clause 2 empty on every case, and |N| equals the registry on each.
- Minimum ratios: k10 522, k28 x1 5783, k28 x2 2455, p3 / nested / nested_inner 8.70, inversion 124 / 19.5 / 12.2, period2 19.8, heteroclinic **3.038**.
- The fingerprint guard passes inside the suite.

**k28 series.** Format: registry / crossings in `T(q0)` / holes / seed / closed.

| blast | base | branch |
|---|---|---|
| 0-5 | identical | identical |
| 6 | 86 / 86 / 4 / 4 / 4 | 88 / 85 / 4 / 4 / 4 (all three detections of the run registered; run {75, 76, 77} out of N whole, median ratio 30.2; 60 bridges vs 61) |
| 7 | 146 / 146 / 4 / 4 / 4 | 146 / 143 / 4 / 4 / 4 (no blast-7 run: the blast-6 run is no longer split at a median) |
| 8 | 248 / 248 / **24 / 15 / 23** | 244 / **143** / **4 / 4 / 4**, reliable on both |

- At blast 8 the unresolved set is {75, 76, 77} by clause 2 plus 148 (ratio 2.572) and 245 (ratio 2.991) by clause 1.
- 245 is the pip orbit's next point. The pip is id 3 and `iterate(3, 8) = 147` at u = 4.869e7. The scaled `f^9(q0)` is u = 2.78699e8, and 245 sits at u = 2.78699e8 (2.6e-6 relative).
- So `T(q0)` stops at `f^8(q0) = 147`. Branch blast 8 has exactly blast 7's 143 crossings, and the base's blast-8 reference switch, (185, 186) and (186, 187) at u ≈ 2.36e8, lies past the cut.
- **Label slip:** the 2026-10-12 decision entry and the coder call 245 "f^8(q0)". Counted from the pip it is `f^9(q0)`, and the test engineer's evidence agrees with my count. The mechanism is the same either way.

**k32 series, blasts 9-10.**
- 523 (ratio 1.477) is `f^10(q0)`: scaled u 2.409275e9, detected 2.409280e9.
- The branch trellis stays at blast 8's: 300 crossings, holes 56, seed 34, closed 45, `is_reliable` True.
- The base had 521 crossings, holes 58, seed 34, closed 39, `is_reliable` **False**.
- Blasts 0-8 are identical. So "reliable False → True" means **x9-10 now show the x8 topology**. No deeper topology was resolved.

**nested_outer series, blasts 9-11.**
- 447 (ratio 1.775) is `f^10(q0)`: scaled u 6.169431e8, detected 6.169446e8.
- The branch trellis stays at blast 8's 276 crossings; the base had 448.
- Holes / seed / closed are 12 / 9 / 9 on both, and both are reliable. There is one fewer workbench bridge per step (321 / 524 / 863 vs 322 / 525 / 864), because nothing is cut at 447.

**The two regression builds that now run with `resolution_ratio_threshold = None`** (clause 1 off, clause 2 on):
- **`cases.build_period3_coarse_blasts`** (`cases.py:594`), for `test_coarse_blasts_rows_by_signs`:
  - its 4 detections with ratio below 3 (minimum 0.65) stay in N;
  - so the 2-4 node hole bridges (42, 41) and (43, 42) still reach the partition;
  - under the default threshold, no build the coder scanned (6 cutoffs × 4-8 blasts) reaches NUM-1.
- **The order-flip build in `test_evidence_diagnosis.py:110`**:
  - it is now `build_k28(blasts=0, through="zones")` plus 9 hand blasts at `None`;
  - the flip is back on (243, 4, 407, 6), ratios above 5000;
  - clause 2 still holds out {75, 76, 77}.

**Refined test builds (1e-6 instead of 1e-5):**
- the `inversion_period2_10` deep run (`cases.py` BUILDERS) and `_p2inv10` in `test_pip_orbit_cut.py`;
- at 1e-5, crossings 38 and 39 (ratios 1.03, 2.33) are out of N, the reference moves to (30, 31), (31, 32), and the topology drops from 13 / 14 / 17 to 10 / 6 / 6 (coder's figures, *believed*);
- the law case `inversion_period2` keeps 1e-5 (new keyword, default unchanged);
- the 1e-5 deep build stays under test in `test_p2inv10_coarse_cut_holds_what_is_resolved`.

**Elsewhere** (coder's probe, *believed*): `k28_deep_blast_sliver` has 1 detection out of N; `near_tangency` has 152 of 261.

**Every path:**
- An unresolved detection no longer disqualifies a strong pip, sets the stable trim, or becomes a target of `grow_until_iterates_closed`.
- With `threshold=None` the report now fires for clause 2. Before, it was silent.
- The WARNING text is new: one line per detection, then "refine: lower `area_cutoff` or grow the stretch; NUM-3/4".

**Suite cost:** 162 s → 200 s (+38 s), mostly the per-step recompute in `blast_series`. Within the plan's estimate.

## 5. What was verified, and how

- **Full suite, run by me:**

  | where | passed | skipped | deselected | xfailed |
  |---|---|---|---|---|
  | branch `76d0710` | 1157 | 38 | 5 | 1 |
  | base `1029fcc` | 1117 | 38 | 5 | 1 |

  The only xfail on both is the KNOWN_ISSUES C3/INV-11 entry.
- **Diff:** I read the whole diff. `git diff --name-only 1029fcc..HEAD` lists exactly the 18 paths: the plan's 14 plus the four you added on 2026-10-12. There are no other paths, and no `src/` change after `dc513f7`.
- **The rule:**
  - The `--k28x8` probe, on both trees, recomputes each `h` independently. On the branch, the probe's set equals `registry.unresolved_crossings(3)` (`agrees=True`, `anycurve_extra=[]`).
  - On the base, the same rule applied to the recompute's numbering gives 8 ids, {1, 2, 75, 125, 250, 251} by clause 2 and {148, 249} by clause 1. The base's own recompute drift is 248 → 252.
- **The deep-run numbers and the pip-orbit identification of 245 / 523 / 447:** my scripts, as above.
- **Read, not re-run:**
  - the test engineer's mutants M1-M9 (each rewritten test fails when its defect is put back);
  - the base run of the four rewritten tests (8 passed);
  - the reviewer's break-its (i)-(iv);
  - rubric items 7, 8, 9, 12 and 13 (the coder's and reviewer's checks).

## 6. What did not go to plan

- **Round 1 stopped at rubric 1-2.** Four tests outside the list failed: coarse blasts, order flip, the p2inv10 walk test, and the oracle `check_trellis_is_cut_at_the_pip_orbit`. You widened the scope to 18 paths. Fixing the oracle exposed a fifth, `test_open_p3_deep_runs[inversion_period2_10]` checking nothing, which was fixed in `cases.py` (already on the list).
- **The pip-orbit dependence was not predicted.** The plan and the analyst expected k28 x8 to be decided by excluded ids at holes or image chains; none is. The plan expected k32 and nested to be untouched topologically ("none bounds a hole …"). They change too, through the cut.
- **Registration order.** Removing the vote also removed the sort that kept blast ids reproducible, and two id tests failed. `e201ddc` sorts in `populate_intersections_for_manifold`. That was not in the plan's steps.
- **Hand-built trellises.** `Trellis.__init__` defaults `resolved_ids` to the whole registry (`Trellis.py:204`). This applies to hand-built trellises and to the `MinimalTrellis` sparse snapshot. The fix belongs at `MinimalTrellis.py:411`, outside the list. It is harmless today, because no reader of N runs on the sparse snapshot (reviewer's reading, *believed*).
- **The B8 collision** (two sources claiming image 213) was expected to vanish. Nobody checked it.

## 7. Where I am least sure this is what you meant

1. **An unresolved pip-orbit point freezes discovery, not just shortens it.** Your decision says it "shortens `T(q0)` by one iterate".
   - In practice, the unresolved point is never linked, so the forward chain breaks there (`Trellis.py:322-323` walks `iterate`). Every later blast leaves `T(q0)` exactly where it was.
   - Measured: k32 blasts 9 and 10, and nested_outer blasts 9, 10 and 11, are each identical to blast 8.
   - The registry keeps growing (898 at k32 x10, 1218 at nested x11), but the topology no longer sees it until you rebuild at a finer `area_cutoff`.
   - That matches "refine" as the remedy. But it means a blast series past such a point is wasted work, and nothing says so beyond the one WARNING.
   - If you would rather place an unresolved orbit point by `image_cdist`, the way PK-6 places a missing one, so the cut keeps its depth, that is a different rule.
2. **Two regression tests now guard their defects only with clause 1 off.**
   - `test_coarse_blasts_rows_by_signs` guards NUM-1 / D6. `test_an_order_flip_gives_no_verdict` guards C21.
   - Under the default threshold, neither defect is reachable by any build anyone found.
   - The mutants show both tests still bite, so this is a build rework and not a weakening. Your options are to keep this form, mark them NOT_APPLICABLE-style, or re-guard them at unit level.
3. **`CLAUDE.md` is stale on this point until B12.**
   - It still says "flags only … nothing is discarded" (Resolution guards).
   - It still states the k=2.8 x8 fixture facts (15 → 23 classes) and "reliable at 8 blasts". On this branch those facts are 4 / 4 / 4.
   - The plan deferred both to B12. Anyone reading `CLAUDE.md` after the merge will be told the opposite of what the code does.

*Scratch evidence:*
- `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/`: `audit_head.txt`, `audit_base.txt`, `audit_probe_k28x8_{branch,base}.txt`
- the same directory's `audit/`: `series.py`, `piporbit.py`, `*.k28_series.txt`, `*.k32_series.txt`, `*.nested_outer_series.txt`, `pip.*.txt`
