# Findings: B1 `resolved-crossings`, the step-0 probe and the plan section checked against `HEAD`

*Analyst, 2026-10-12. Job `findings`. (The analyst's Write tool refused the artifact; the orchestrator saved this text verbatim.)*

**The working tree is dirty.** `HEAD` = `1029fcc` on `main`. Modified: `.refactor/decisions.md`, `.refactor/knowledge/{testing,topology}.md`, `docs/alignment_plan_2026-10-09.md`, `docs/code_fidelity_audit_2026-10-09.md`, and `.refactor/knowledge/numerics.md` (my edit). There is a staged rename of `artifacts/request_origin_analyst.md` into `.refactor/runs/2026-10-12-oracle-laws/` (the working copy is untracked under `artifacts/`). B0's run files are untracked. Nothing under `src/` or `tests/` is dirty. The change branches from `HEAD`.

**What kind of request.** A refactor that changes a definition (Q13 (a)). It touches numerics (the crossing set N) and the trellis that reads it. I was asked to verify the plan section, not to re-plan it, so everything below is a difference from the plan section or a fix to it.

## Short answer

The plan's direction holds up. On the 11 law cases nothing is held out under either clause. The 248 → 252 drift is real, and the vote fires only on the blast path. Six details in the plan are wrong or incomplete in ways that would make the rubric fail or hide behaviour changes:

1. **The rule as written holds out 8 detections on k28 x8, not 6.** Clause 2 also catches the runs' medians 75 and 125.
2. **Clause 2 is asymmetric when each detection is compared against its own `max(h)`.** One member of a run can read as resolved.
3. **On the blast path a run can straddle two `iterate_bridge` calls.** So `run_separation` cannot be set per call.
4. **There are more than four readers of the registry as the crossing set.** The blast cut (`BridgeIterator.py:375`, a file not in the list) is not among the four named. StrongPip, the stable trim and the growth-driver targets also read it.
5. **The k28 series drifts from blast 6, not 8.** k28 x8 is also not the only series or builder whose N changes.
6. **Two things are missing from the file list.** `IntersectionRegistry.py` (step 2's methods live there) is not listed. And `Intersection` stores no chord deviation for clause 2 to compare against.

## The step-0 probe (verified: run at `HEAD`)

**Where the script is.** It is not in the repo. It lives in `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b1/`: `probe.py` is the library, and `laws.py`, `series.py`, `k28_variants.py`, `nc_detail.py`, `incut.py`, `nt.py` are the drivers.

**How it works.** It monkeypatches `Tangle._resolution_ratio` to record `(h_u, h_s)` per detection, using the same `chord_deviation_batch` and `_outer_neighbours` the library uses. It also wraps `_collapse_noise_crossings` to log every firing.

- **Clause 1:** `resolution_ratio < 3`.
- **Clause 2:** the plan's text. Another detection with the same `(u_key, s_key)` sits on the same unstable segment or an adjacent one (sharing a node). Its distance along the polyline (through the shared node) is at most this detection's `max(h_u, h_s)`.
- **Cross-check:** the same pairs found by arc length along the whole curve object, with no adjacency restriction. Across the whole corpus it found nothing the adjacency rule missed.

**Law cases (all 11).** Clause 1 is empty, clause 2 is empty, and the collapse fired 0 times.

| case | detections | min ratio |
|---|---|---|
| k10 | 7 | 522.4 |
| k28_one_blast | 5 | 5783 |
| k28_two_blasts | 9 | 2455 |
| p3 | 27 | 8.704 |
| nested | 36 | 8.704 |
| nested_inner_blasts | 40 | 8.704 |
| inversion | 12 | 124 |
| inversion_one_blast | 18 | 19.55 |
| inversion_two_blasts | 24 | 12.18 |
| inversion_period2 | 17 | 19.8 |
| heteroclinic | 112 | **3.038** |

The plan's figure for period2 is 33.7; at `HEAD` it is 19.8, because B0 grew it to 10 stable steps. Heteroclinic's 3.04 is confirmed.

**Series.** Steps 0-5 of k28 are clean. "Registry" counts include anchors.

| series step | registry | min ratio | clause 1 | clause 2 | collapse fired | recompute |
|---|---|---|---|---|---|---|
| k28 b6 | 86 | 18 | – | – | yes, run of 3 at u=8356210.197 (drops ratios 2297, 4642) | **86 → 88** |
| k28 b7 | 146 | 6.13 | – | – | yes, run of 3 at u=47831087.57 (drops 624, 633) | **146 → 150** |
| k28 b8 | 248 | 2.57 | 148 (2.57), 249 (2.99) | – (only the medians are registered) | – | **248 → 252**, new ids 1, 2, 250, 251 |
| p3 b0-9 | 30 → 67 | 8.70; 3.038 from b6 | – | – | – | 67 → 67 |
| inversion b0-4 | 13 → 47 | 124 → 12.18 | – | – | – | 47 → 47 |
| k32 b0-8 | 4 → 300 | 3.24 at b8 | – | – | – | – |
| k32 b9, b10 | 522, 898 | 1.477 | **523** | – | – | 898 → 898 |
| nested_outer b0-8 | 34 → 276 | 3.094 at b8 | – | – | – | – |
| nested_outer b9-11 | 448 → 1218 | 1.775 | **447** | – | – | 1218 → 1218 |

**k28 x8 after the base recompute (252 crossings).** Clause 1 = {148, 249}. Clause 2 = {1, 2, **75**, 250, 251, **125**}.

| id | ratio | own max(h) | nearest same-key detection, along the polyline |
|---|---|---|---|
| 1 | 2297 | 4.7e-6 | 75 at 2.0e-8 |
| 75 (median) | **30.2** | 4.7e-6 | 2 at 8.5e-10 |
| 2 | 4642 | 4.7e-6 | 75 at 8.5e-10 |
| 250 | 624 | 8.4e-6 | 125 at 6.1e-8 |
| 125 (median) | **12.0** | 2.9e-6 | 250 at 6.1e-8 |
| 251 | 633 | 8.3e-6 | 125 at 1.8e-7 |

Measured along the polyline the members sit 8.5e-10 .. 1.8e-7 apart. The plan's "2.9e-7 .. 8.5e-5" reads like unstable-cdist gaps. The margin between separation and `h` is 1 to 4 orders of magnitude.

**Other builders** (not law cases, but tests use them). Clause 2 is empty on every one.

| builder | clause-1 detections |
|---|---|
| `k28_deep_blast_sliver` | {423}, ratio 1.56 |
| `p3_coarse_blasts` | 4, min ratio 0.65 |
| `inversion_period2_10` | 2, min ratio 1.03 |
| `near_tangency` | 152 of 261 |
| `heteroclinic_17u`, `p3_15`, `p3_16`, `p3_16_12`, `p3_6_blasts`, `nested_unblasted` | none |

`orientation_reversing` raises `ValueError` at build; it is a set-aside case.

## The plan section's claims, one by one

1. **"`_collapse_noise_crossings` is called only from populate ← `BridgeIterator.py:363`; `_purge_crossing`'s only caller is `Tangle.py:626`."** *Verified* by grep over `src`, `tests`, `scripts` and `stress`. The definitions are at `Tangle.py:569` and `:629`, and the call is at `:567`.
   - Tests call it at `test_noise_crossing_collapse.py:46,64,78`.
   - `stress/numerics/e09_intersections.py:13` names it in a docstring. That is outside the eight files and harmless, since rubric 2's grep covers only `src` and `tests`.
   - The vote groups runs across all stable keys. Clause 2 groups by `(u_key, s_key)`.

2. **"248 → 252 on k28 x8."** *Verified* (`len(registry)`, including the one anchor). But "the collapse fires ONCE, on k28 x8" is **wrong**.
   - It fires at blasts 6 and 7 of the k28 series, once per run of three.
   - The drift starts at blast 6: 86 → 88 at blast 6, 146 → 150 at blast 7 (verified, each in a fresh process).
   - So `test_recompute_keeps_the_crossing_set[k28_series]` fails on base at step 6, not at step 8.
   - Step 6 also recomputes inside the series. On base, that recompute registers the outer run members, the next blast maps bridges cut at them, and base's later numbers drift from today's. The rubric is still met, but the base run will not print "248 ≠ 252".

3. **"Uncollapsed 40 / 23 / 31; base 24 / 15 / 23; both reliable."** *Verified.* Holes / seed / closed are 24 / 15 / 23 on base and 40 / 23 / 31 with the vote disabled. `is_reliable` is True for both.

4. **"`unresolved_ids(3)` = {1, 2, 250, 251} by clause 2 and {148, 249} by clause 1; six exclusions; the medians 75 and 125 are kept."** **Wrong, in three ways.**
   - **(a) The medians are held out too.** By step 2's own definition, the medians 75 and 125 are clause-2 detections: their nearest neighbours sit 8.5e-10 and 6.1e-8 away, against `h` of 4.7e-6 and 2.9e-6. That gives **eight** exclusions on base's recompute numbering, which matches A.2's "a run is unresolved as a whole". 75 and 125 are in base's N, so B1 removes two crossings the base used. This agrees with B8: the 125/251 → 213 collision goes away, and 213 loses its preimage.
   - **(b) The ids in rubric 5 will not exist on the branch.** 1, 2, 250 and 251 are numbers the base recompute assigns via `reindex_from`. With no vote, the blast path registers the run members as it detects them. In the uncollapsed build the runs are {74, 75, 79} at blast 6 and {126, 128, 129} at blast 7, and clause 1 is {155, 251}. Rubric 5 has to state the set by rule (every member of each run, plus the detections with ratio < 3), never by id.
   - **(c) Clause 2 is asymmetric per detection.** In the uncollapsed build, 129 sits on a short junction segment between two image curves: its own `h` is 3.6e-10 and its neighbour is 5.6e-8 away. So `run_separation > max(h_u, h_s)` reads 129 as resolved, while 126 and 128 (whose neighbour is 129) read unresolved. That is "keep one member of the run" again, chosen by which segment happens to be short.

5. **"Four call sites, one name: `create_bridges`, `infer_iterates`, `Trellis.from_workbench` (branch orderings; `own_intersection_ids`)."** **Incomplete.** At `HEAD` the registry is read as N in these places:

   | place | what it reads | named in B1? | file in the list? |
   |---|---|---|---|
   | `TangleWorkbench.create_bridges` :1786, plus its partial-piece count :1804 (same list) | crossings to cut at | yes | yes |
   | `BridgeIterator.iterate_bridge` :375 | **the blast cut**: calls `Tangle.create_bridges(new_intersections, …)` directly | no, unless the filter sits inside `Tangle.create_bridges` | **no** |
   | `BridgeIterator.iterate_bridge` :394 | inference sources = the call's new detections | covered if `infer_iterates` filters | — |
   | `infer_iterates` :106 (sources) and `_bracketed_image` :215 (candidates, via `candidates_near_stable_cdist`) | crossings to link | yes (two places in one file) | yes |
   | `Trellis.from_workbench` :436 | branch orderings; `own_intersection_ids` :638 is derived from the branches, not read again | yes | yes |
   | `StrongPip.is_strong_pip` :168 | disqualifiers, looping over the whole `trellis.registry`; an unresolved detection could disqualify a pip | **no** | **no** |
   | `TangleWorkbench.trim_stable_manifolds` :2158 | the outermost registered crossing, which sets the trim | no | yes |
   | `TangleWorkbench._frozen_ids("all")` :1091 | targets of `grow_until_iterates_closed`; an excluded detection is never linked, so the driver runs to its cap and raises | no | yes |

   These only report and need no change: `ResonanceZone.py:561` and `TangleSession.py:1546` (zone `intersection_ids`), `diagnostics._tangle_box` :242, the `DualGraph` markers at :1134 (anchors only, which are synthetic), `plot_intersections`, and `Trellis.intersection_ids` / `by_*_cdist` (no caller in `src`).

6. **"Series p3, inversion, k32, nested_outer: recompute already exact; k28 x8 is the one case the rule decides."** The recompute is exact (*verified* at each series' last step: same count, same ids). But **N changes elsewhere too, under clause 1**:
   - k32 at blasts 9-10 (523), nested_outer at 9-11 (447), and the builders `k28_deep_blast_sliver`, `p3_coarse_blasts` and `inversion_period2_10`.
   - 523, 447, 148, 249, 75 and 125 are all in the cut trellis as bridge ends, for example `(521,523)`, `(446,447)`, `(148,149)`, `(248,249)`, `(74,75)`, `(75,76)` (*verified*).
   - None of them bounds a hole, sits in a reference pair, or ends a kept minimal-trellis bridge. So "the holes that depended on 148 / 249" comes out as **none by bounding id**. Any dependence is indirect: through image chains (148 ← 88, 249 ← 147), the evidence table and pseudoneighbor trajectories.
   - These tests run on the affected series and builders: `test_class_count_under_blasting.py` (the count never decreases; it stays constant while the holes do), `test_coarse_blasts_rows_by_signs.py`, `test_deep_blast_sliver_beyond_pip.py`, and `test_open_p3_deep_runs.py` (`inversion_period2_10`).
   - *Believed:* the effect is small, since none of these crossings is in the minimal trellis. *Not verified:* that those tests still pass.

7. **k28 x8 cannot be predicted from base (*believed*).** Under B1 nothing is cut at the blast-6 run, because none of it is resolved. One bridge then spans the run, and one image curve carries its image at blast 7. On base, the cut at 75 split that image over two curves, which is exactly why the blast-7 run straddles two curves (item 4c). So the branch's blasts 7-8 can differ from base. The plan already says "re-derived", but rubric 5's exact set cannot be fixed in advance.

## Two design holes the plan does not cover

Both were seen on the uncollapsed build (*verified*).

- **Runs straddle calls.** The blast-7 run is 126 on one image curve and 128, 129 on another, on adjacent segments that share a node. Each image curve comes from its own `iterate_bridge` → `populate` call. So a `run_separation` "set by `resolve_crossings` and `populate` alike" within one call's list misses the pair, and the blast path and the recompute disagree.
- **A crossing can be demoted after it is used.** If clause 2 is evaluated against the registry instead, a later call in the same blast can make an already-registered, already-cut, already-linked detection unresolved. A law that compares only crossing COUNTS would still pass, because N is a registry function, while the bridges differ.

## Recommendation: the smallest change to the plan, not a re-plan

1. **Make clause 2 a property of a pair, read symmetrically.** Two same-key detections on the same or adjacent unstable segments are both unresolved when they sit closer along the polyline than the larger `max(h_u, h_s)` of the two. A run of two or more is out whole. This closes 4a and 4c.
2. **Store the deviation.** `Intersection` gains `chord_deviation = max(h_u, h_s)`, set where `resolution_ratio` is set. Without it `resolved_ids` cannot be evaluated.
3. **Evaluate clause 2 in the registry, not per call.** `resolved_ids` becomes a registry function that recomputes the runs from the stored segment nodes over every registered detection. Both paths then agree on N by construction.
4. **Extend step 6's law to bridges.** Also compare the bridge-id set after `compute_intersections(preserve_ids=True)` plus `rebuild_bridges()`. This catches the order hazard of point 3 at no extra cost, since the recompute is already there.
5. **Fix the file list.**
   - Add `numerics/IntersectionRegistry.py`.
   - Add `numerics/BridgeIterator.py` to filter at :375 and :394, or state that the filter lives inside `Tangle.create_bridges`.
   - Add `topology/StrongPip.py`, so that :168 iterates the trellis's computed ids instead of the registry.
   - Cover `trim_stable_manifolds` :2158 and `_frozen_ids` :1091 in `TangleWorkbench`, which is already listed.
   - That is 10-11 files, not eight.
6. **Fix rubric 5 and the behaviour text.**
   - State the k28 x8 set by rule (run members plus ratio < 3), not as the ids {1, 2, 250, 251, 148, 249}.
   - Say the drift starts at blast 6.
   - Name k32 b9-10, nested_outer b9-11, `k28_deep_blast_sliver`, `p3_coarse_blasts` and `inversion_period2_10` as clause-1 behaviour changes, reported and not pinned.

## Out of scope (noticed, not for B1)

- **A test that sets the threshold to `inf`.** `test_crossing_resolution.py:195` (not in the file list) sets `resolution_ratio_threshold = inf` to flag every crossing. Under B1 that makes N empty, so its blast probably yields only partial pieces, and the test's meaning changes. Separately, `resolved_ids(None)` needs a meaning, since a None threshold disables the flag today.
- **The B0 warning test.** `test_resolution_warning_follows_oracle.py` checks "warn iff clause 1 or a same-stability hit". Clause 2 is empty on its three builds (verified), so it stays consistent.
- **How often the WARNING fires.** Step 4's "ONE WARNING per call" conflicts with the current "report each finding once" dedup (`_reported_flags`, :1489), and rubric 5's "fires once per call naming the six" relies on step 4. One of the two has to go.
- **A stale fixture fact.** CLAUDE.md's "k=2.8 at 8 blasts 15 -> 23 (3 flagged)" will go stale; that belongs to B12.
- **Test time (*believed*).** A recompute per series step costs about 1-4 s at the deep steps (NUM-12 measured 3.8 s at 9 k28 blasts), so the regression tier might grow by 30-60 s.

## Knowledge files

`.refactor/knowledge/numerics.md`: I re-dated the header and added the section "Resolved crossings — what B1 found at `1029fcc`". It records the readers with line numbers, the probe's numbers, the run geometry and the clause-2 asymmetry. I did not study the other knowledge files and left them alone.
