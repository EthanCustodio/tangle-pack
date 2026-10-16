# Plan: B1 `resolved-crossings`

*Planner, 2026-10-12. This is slice B1 of `docs/alignment_plan_2026-10-09.md`, with the analyst's six corrections folded in (`artifacts/findings_analyst_planner.t1.md`). The base is `HEAD` = `1029fcc`. There is no earlier verdict on this run, so nothing new goes into `.refactor/decisions.md`.*

## 1. The answer

The blast path votes on near-coincident detections and no other path does: an odd run keeps its median, an even run is dropped. So a plain recompute changes `N`. On the k28 series this starts at **blast 6** (86 → 88 crossings), then 146 → 150 at blast 7 and 248 → 252 at blast 8.

Deleting the vote on its own turns k28 x8 from 24 holes / 15 classes / 23 closed into 40 / 23 / 31, and both versions report "reliable".

Under Q13 (a), a detection is in `N` only if the polylines resolve it. It must pass two clauses:
- **Clause 1:** its ratio is ≥ 3.
- **Clause 2:** no detection on the same two curves sits closer along the unstable polyline than the larger chord deviation of the pair. A run of such detections is out as a whole, median included.

The probe, run at `HEAD`, finds no unresolved detection on any of the 11 law cases. The lowest ratio is heteroclinic's 3.038.

## 2. Recommendation

Change it, as decided in Q13 (a). Three of the analyst's corrections shape the steps:
- **Clause 2 is a symmetric pair rule.** With a per-detection test, one member of a run can pass on a short segment.
- **Clause 2 is evaluated by the registry.** On the blast path a run can straddle two `iterate_bridge` calls, so a per-call check would miss it.
- **`N` is read at every place that uses the registry as the crossing set.** There are seven such places, not four.

## 3. Scope

**Slug** `resolved-crossings`.
- **Branch:** `refactor/resolved-crossings`, in the worktree `../tangle-pack-resolved-crossings`.
- **Base worktree:** `../tangle-pack-resolved-crossings-base`, detached at `1029fcc`.

**File list (14 paths). Nothing else is edited.**

| # | path | why |
|---|---|---|
| 1 | `src/tanglepack/numerics/Tangle.py` | vote deleted; chord deviation stored |
| 2 | `src/tanglepack/numerics/Intersection.py` | `chord_deviation` |
| 3 | `src/tanglepack/numerics/IntersectionRegistry.py` | `resolved_ids`, `run_separation`, `unresolved_crossings` widened |
| 4 | `src/tanglepack/numerics/TangleWorkbench.py` | `create_bridges`, `trim_stable_manifolds`, `_frozen_ids`, `report_resolution`, threshold docstring |
| 5 | `src/tanglepack/numerics/BridgeIterator.py` | the blast cut (:375) and its inference sources (:394) |
| 6 | `src/tanglepack/numerics/IterateInference.py` | sources (:106) and candidate images (:215) |
| 7 | `src/tanglepack/numerics/diagnostics.py` | two docstrings that define "unresolved" by the ratio alone (no code) |
| 8 | `src/tanglepack/topology/Trellis.py` | branch orderings from `N`; records `resolved_ids` |
| 9 | `src/tanglepack/topology/StrongPip.py` | disqualifiers from `N` only (:168) |
| 10 | `tests/unit/numerics/test_noise_crossing_collapse.py` | deleted |
| 11 | `tests/unit/numerics/test_resolved_crossings.py` | new |
| 12 | `tests/cases.py` | `blast_series` records the recompute |
| 13 | `tests/regression/test_class_count_under_blasting.py` | the recompute law |
| 14 | `tests/_tools/resolution_probe.py` | new: the step-0 probe (not collected) |

**Size.** This is above the usual targeted-run size: about 14 files and about 400 changed lines, of which the probe is about 150. I am not splitting it. The bug is precisely a definition of `N` that some readers apply and others do not, so a slice that reached only some of the readers would leave the same bug in place.

Compared with your eight files, it adds:
- the registry, where the rule lives;
- the blast cut and StrongPip, two readers the section missed;
- diagnostics, for two docstrings;
- the probe, which the section already named.

## 4. Steps

**0. Probe.** Move the analyst's scratchpad probe (`probe.py` and its drivers) into one script, `tests/_tools/resolution_probe.py`, with the flags `--laws`, `--series`, `--builders` and `--k28x8`.
- It recomputes each detection's `(h_u, h_s)` independently, the same way the library does, through a monkeypatch of `Tangle._resolution_ratio`. It never reads the new stored field, so it remains an oracle.
- It prints:
  - per build: the detection count, `|N|`, the clause-1 ids, the clause-2 pairs (each with separation and both deviations), and the minimum ratio;
  - for k28 x8: holes / seed / closed classes, and which of the base's holes depend on an excluded id, by bounding id or by image chain.
- Run it on the base worktree first, then on the branch. Attach both outputs to the report.

**1. Store the chord deviation.**
- `Tangle._resolution_ratio` returns `(ratio, max(h_u, h_s))`.
- `_resolve_crossing_pair` sets `Intersection.chord_deviation` beside `resolution_ratio`; it is None on a synthetic crossing.
- There is no `run_separation` field: a separation stored per call would be wrong when a run straddles two calls.

**2. The rule, in the registry.**
```python
def run_separation(self, id) -> Optional[float]:
    """Distance along the unstable polyline to the nearest other detection of the same
    (u_key, s_key) on the same or a node-sharing unstable segment; None when alone."""

def resolved_ids(self, threshold: Optional[float]) -> frozenset[int]:
    """N: synthetic, or (ratio >= threshold, when there is a threshold) and in no run.
    Two detections are a run when their separation <= max of their two chord deviations;
    a run is out whole (Q13 (a), A.2)."""

def unresolved_crossings(self, threshold) -> list[int]:   # existing name, now the complement, unstable order
```
- **Separation along the polyline:**
  - Same segment: `|coords_a - coords_b|`.
  - Adjacent segments: the distance from each detection to the shared node, summed.
  - "Segment" means the stored original segment (`unstable_segment`, captured before any separator is spliced in).
  - Adjacency means node identity.
- Neighbours are found among detections that are consecutive in unstable-cdist order within one `(u_key, s_key)` group.
- The result is cached per registry `generation` and threshold.

**3. Read `N` everywhere the registry is the crossing set.** Each place uses `resolved_ids(workbench.resolution_ratio_threshold)`:
- `TangleWorkbench.create_bridges`: the cut list and the partial-piece count.
- `BridgeIterator.iterate_bridge`: the cut at :375 and the inference sources at :394. An unresolved detection lies inside one bridge.
- `infer_iterates`: the sources and the candidate images. No link starts or ends at an unresolved detection.
- `Trellis.from_workbench`: the branch orderings, with the set recorded as `Trellis.resolved_ids`.
- `StrongPip.is_strong_pip`: skips any `r_id` that is not in `trellis.resolved_ids`.
- `TangleWorkbench.trim_stable_manifolds`: trims at the outermost crossing of `N`.
- `TangleWorkbench._frozen_ids("all")`: targets only crossings of `N`. Otherwise the driver would run to its cap on a detection it can never link.

**4. Delete the vote.** Delete `_collapse_noise_crossings` and `_purge_crossing`. `populate_intersections_for_manifold` returns every detection.

**5. `report_resolution`.**
- The WARNING names each unresolved detection with its id, coordinates, ratio, clause, and (for clause 2) its partner and separation.
- It ends with the remedy: "refine: lower `area_cutoff` or grow the stretch; NUM-3/4".
- The existing de-duplication stays: one WARNING per call, and only when the call adds a detection not reported under the current numbering. This settles the conflict the analyst raised: step 4 of B1 as written would repeat the warning on every recompute.
- `ResolutionReport.unresolved_crossings` now holds both clauses. The two docstrings in `diagnostics.py` say so.

**6. Threshold docstring.** The derived bound is 1, because separation equal to the polygon error means the crossing cannot be told apart. 3 is the margin. The threshold is part of the definition of `N` (A.2), and the corpus minimum is 3.038 (heteroclinic).

**7. Tests.**
- **Delete** `test_noise_crossing_collapse.py`.
- **New: the zig-zag.** `test_resolved_crossings.py` builds a hand-made unstable polyline that zig-zags across a stable segment, giving three detections within the chord deviation. It asserts:
  - all three are registered;
  - none is in `resolved_ids`;
  - one bridge runs across them;
  - no `iterate_table` entry touches them;
  - the WARNING is logged once.
- **New: the symmetric pair.** Two detections, one on a short segment (small deviation) and one on a long segment (large deviation). The separation is between the two deviations, and both must be out.
- **New: a clean crossing.** A clean transversal crossing with ratio ≥ 3 is resolved.
- **Series recompute.** In `cases.blast_series`, every step records `(N, bridge ids)` before and after `compute_intersections(…, preserve_ids=True)` plus `rebuild_bridges()`. `test_recompute_keeps_the_crossing_set(name)`, run over `cases.SERIES`, asserts both pairs are equal at every step.

## 5. Behaviour

**Law cases (11):** none. The probe at `HEAD` finds both clauses empty and the collapse fired 0 times.

**k28 series:** changes from blast 6. Every member of the runs at blasts 6 and 7 is out of `N`, so nothing is cut at them and no link touches them. Because of that, the branch's blasts 7 and 8 can differ from the base's (the blast-7 image is no longer split at the old median).

**k28 x8:** the unresolved set is **every detection with ratio < 3, plus every member of every same-curve run, the run's median included.**
- On the base's recompute numbering this is 8 detections: 2 by clause 1 and 6 by clause 2.
- On the branch the ids and possibly the set itself differ. It is re-derived and printed, never pinned.
- Holes / seed / closed classes are printed beside the base's 24 / 15 / 23 and the uncollapsed 40 / 23 / 31. The base's holes that depended on an excluded id are named. By bounding id, none did (analyst).
- The B8 collision (two sources claiming image 213) should disappear. This is reported, not asserted.

**Clause-1 changes elsewhere.** These detections leave `N` today. They are reported, never pinned:

| where | detections leaving `N` | ratio |
|---|---|---|
| k32 series, blasts 9-10 | one | 1.48 |
| nested_outer series, blasts 9-11 | one | 1.78 |
| `k28_deep_blast_sliver` | one | 1.56 |
| `p3_coarse_blasts` | four | min 0.65 |
| `inversion_period2_10` | two | min 1.03 |
| `near_tangency` | 152 of 261 | – |

On the series, none of them bounds a hole, sits in a reference pair, or ends a kept minimal-trellis bridge (analyst, verified). Any effect would come through image chains and pseudoneighbor trajectories, and is reported.

**Every path:** an unresolved detection no longer disqualifies a strong pip, sets the stable trim, or becomes a target of `grow_until_iterates_closed`.

## 6. Not doing

- **No local refinement.** Option (b), NUM-3/4, the remedy the WARNING points to, is not built.
- **No change** to `DEDUP_RTOL`, to the same-stability flags, or to the margin 3.
- **No new threshold setter.** See Q3.
- **`CLAUDE.md`'s "flags never discard" text and the k=2.8 x8 fixture fact are left alone.** They belong to B12, which rewrites `CLAUDE.md`. Until then they are stale on this one point.
- **Tests outside the list stay as they are.** These include `test_crossing_resolution.py` (its `inf` threshold blast), `test_resolution_warning_follows_oracle.py`, `test_deep_blast_sliver_beyond_pip.py`, `test_coarse_blasts_rows_by_signs.py` and `test_open_p3_deep_runs.py`. If one fails, the run stops and brings it to you (rubric 1).
- **The test-time cost of the per-step recompute** (estimated +30-60 s on the regression tier) is accepted and measured.

## 7. Rubric

1. On the branch, `MPLBACKEND=Agg env/bin/python -m pytest -q -rxX` passes, and the only xfails are `cases.KNOWN_ISSUES`. None of the tests named in section 6 is edited.
2. Every test in a file outside the list passes on the base commit and on the branch.
3. `git diff --stat 1029fcc..refactor/resolved-crossings` lists exactly the 14 paths of section 3.
4. No behaviour change beyond section 5:
   - the law tier's fingerprint guard passes unchanged;
   - the probe on the branch prints an empty unresolved set under both clauses for all 11 law cases, and `|N|` equals the base's detection count on each;
   - heteroclinic's minimum ratio is printed (3.038 on base), and any value below 3 is owned in the report.
5. `grep -rn "_collapse_noise_crossings\|_purge_crossing\|Collapsed" src tests` returns nothing. `grep -n "% 2" src/tanglepack/numerics/Tangle.py` returns nothing. `Tangle.py` is at least 55 lines shorter (the two deleted methods are 70).
6. `grep -rn "resolved_ids" src` hits:
   - the definition in `IntersectionRegistry.py`;
   - the seven readers of step 3, once each (`TangleWorkbench` ×3, `BridgeIterator`, `IterateInference`, `Trellis`, `StrongPip`);
   - nothing else.

   No reader iterates `registry` as the crossing set without it. The reviewer checks this by reading the seven sites named in the findings.
7. `Intersection` has `chord_deviation`. It is None exactly when `is_synthetic`, which is checked over every law case's registry.
8. The symmetric-pair test passes on the branch. Removing the `max` of the two deviations (comparing each detection to its own deviation only) makes it fail.
9. The zig-zag test fails on the base worktree, where all three detections are cut at, and passes on the branch.
10. `test_recompute_keeps_the_crossing_set[k28_series]` fails on the base at blast 6, the first unequal step, and passes on the branch. The other four series pass on both.
11. On the branch, the `--k28x8` probe lists the unresolved set together with each member's reason. Each member satisfies its clause by the probe's independent `h`, and no detection outside the set satisfies either clause. The probe's set equals `registry.unresolved_crossings(3)`. Recompute leaves `N` and the bridge ids unchanged.
12. In a fresh process, `blast_zone(diagnose=True)` logs one WARNING at the blast that first registers the k28 run (blast 6). It names every member and contains the text "refine". A following `compute_intersections(preserve_ids=True)` logs no new resolution WARNING.
13. The `resolution_ratio_threshold` docstring states "derived bound 1", "margin 3", "definition of N" and the 3.038 corpus minimum.

## 8. Questions

1. **A crossing demoted after use.** A later `iterate_bridge` in the same blast can register a detection that puts an already-cut, already-linked crossing into a run. *Default:* no repair mechanism in B1. Rubric 10 (crossings and bridges equal after a recompute) is the check. If any series step fails it on the branch, the run stops and brings that case back to you.
2. **`resolution_ratio_threshold = None`.** *Default:* clause 1 is off and clause 2 still applies. A run is unresolved whatever the margin, because it is judged against the chord deviation itself.
3. **Changing the threshold mid-session.** *Default:* no setter and no invalidation. The docstring says the threshold is read by every detection pass and every trellis built afterwards, and that changing it after crossings exist needs a `compute_intersections` to re-cut.
