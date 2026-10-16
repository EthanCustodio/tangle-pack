# Findings: B4b `holes-are-combinatorial` (analyst, HEAD `f166573`, 2026-10-13)

*(The analyst's Write tool refused the artifact; the orchestrator saved this text verbatim.)*

**The tree is dirty.** `main` @ `f166573` has uncommitted changes in:
- `.refactor/decisions.md` and `.refactor/knowledge/{loom,numerics,testing,topology}.md`;
- `docs/alignment_plan_2026-10-09.md` (this run checks the B4b section there) and `docs/code_fidelity_audit_2026-10-09.md`;
- a staged rename of `artifacts/request_origin_analyst.md`, plus untracked `.refactor/runs/*` and `artifacts/`.

No tracked file under `src/` or `tests/` is dirty, so a branch from `HEAD` starts from exactly the code read here. `HEAD` holds B4a: `topology/diagnostics.py` exists and `Hole` has no `measured_side` or `side_source`.

**What kind of request.** A refactor (a move plus deletion of dead paths). Targets:
- in `topology/StablePartition.py`: `punch_holes`, `propagate_reference_holes`, `_punch_in_bridge`, `_hole_openings` and the placement helpers;
- `Hole.coords` in `TopologyResults.py`;
- `Trellis.punch_holes(epsilon=)`;
- a new `plotting.hole_point`.

As asked, this checks the plan against the code and does not re-plan the slice.

**Method.** I wrote a probe script (`b4b_probe.py`, in the session scratchpad). It monkeypatches `StablePartition.punch_holes` and `propagate_reference_holes` with instrumented wrappers. The wrapper runs a replica of today's propagation with counters next to the original and checks that the holes are identical (every field, `coords` included). Mismatches: **0 on every call**.

It ran over:
- all 22 `BUILDERS` cases (`orientation_reversing` raises its known `ValueError` before any hole);
- every step of the five `SERIES`: `k28_series` 9 steps, `p3_series` 10, `inversion_series` 5, `k32_series` 11, `nested_outer_series` 12.

On each built `session.topology().trellis` it then measured both placements and the anchor openings. Unless marked *believed*, every claim below comes from this probe.

## Short answer

Yes: a hole can be its five combinatorial fields. The two candidate placements measure as follows:

- **(a) Carried placement** (`hole_point` carries the reference placement back with `map_inv`): reproduces `Hole.coords` **bit for bit**. The max deviation is **0.0** on every hole of every builder, so rubric 6 (1e-12) holds with room to spare.
- **(b) No-map fallback** (`_step_into(midpoint, stable_mid, 0.05, arc_chord)`):
  - moves propagated holes by up to **0.328** (k28_one_blast);
  - puts propagated holes on the other side of their image sub-arc on **all 11** law cases, not 10 of them.

Five points differ from the plan section, and the first two would make the slice fail its own rubric if done as written:

- **A. The step-2 WARNING as specified fires on 6 of the 11 law cases** (rubric 5 fails). The fix is to warn on scaled span ends only. That keeps the `from_table` flag, so part (iii) becomes "use the flag", not "drop it".
- **B. Rubric 2's "≤ 700 lines" cannot be met.** About 1,060-1,090 lines is what B4b reaches.
- **C. The step-2 comment "never in the table, by construction" is false for one span end.**
- **D. `_near_far` is not placement-only.** It orders `bounding_ids`, a fingerprint field, so it stays.
- **E. The file list misses four files:**
  - `tests/helpers/laws.py` (three laws read `hole.coords`);
  - `topology/diagnostics.py`;
  - `test_topology_invariants.py`;
  - `test_minimal_trellis.py`.

## 1. The four non-PDF parts at `HEAD` (verified: Read and grep)

### (i) Placement in the plane

| piece | lines at `f166573` | note |
|---|---|---|
| `punch_holes(..., epsilon=0.05)` | :273; docstring :277-283, :297-298 | |
| direct placement | :337-353 | `_near_far` :340, `_stable_arc_midpoint` :346, `_nearest_arc_point` :347, `chord_len` with `np.linalg` :348, reference branch `_step_into` :350, **else-branch** (non-reference direct pair) :351-352 |
| skip when `arc_point is None` | :348-350 area | never fires (0 everywhere) |
| `Hole(coords=...)` | :367, :1018 | |
| carried point | :499 (`start.hole.coords`), :513 (`_map_backward`) | |
| `_step_into` | :794-812 | |
| `_stable_arc_midpoint` | :815-870 | |
| `_near_far` | :873-894 | **stays** (see D) |
| `_nearest_arc_point` | :897-913 | |
| `_map_backward` | :951-959 | the only `map_inv` / `dynamical_system` use in the module |
| `_punch_in_bridge` placement | :969-970 (`carried`, `span` params), :989-1015 | |
| `_image_arc_midpoint` | :1028-1061 | |
| debug log `"Hole at %s has no openings"` | :1090, in `_hole_marks_on_branch` | **also reads `hole.coords`**, not named in the plan; it must log `bounding_ids` instead, or rubric 2's grep is not empty |

### (ii) The direct punch of backward pairs, the two dedups and the fallback

- **Else-branch:** :351-352.
- **`seen`:** :436-440 (initial set), :526 and :539.
- **`punched_regions`:** :444-448 and :536-538.
- **`_region_key`:** :609-620. Its only users are the dedup lines, so it goes too; the plan does not name it.
- **Fallback:** :455-495, which includes the `fallbacks` sort, the "cannot propagate" WARNING (:480-485) and the "has no punched bridge; propagating backward" INFO (:486-492).

### (iii) The `from_table` flag

No name `from_table` has ever appeared in StablePartition (`git log -S`). The flag is the third element of `Trellis.image_cdist`'s return value. It is dropped by the `[1]` index at **:510**; the plan cites `:513`, which is now the `_map_backward` line.

Plan step 2's expression, `tuple(sorted(image_cdist(...)[1] for end in ids))`, **is today's code**.

### (iv) `_hole_openings`' inline copy of `stable_key_at`

- `_hole_openings` spans :746-791; the anchor block is :768-782.
- It rebuilds `(fp, "stable", orbit_index, branch_index)` from `trellis.unstable_key_of(bridge.id)`. That is `Trellis.stable_key_at` (`Trellis.py:740-764`) inlined. It then flips the branch index when `(sign > 0) != (side == "left")`.

## 2. No recorded pair has `iterate < 0` (verified)

`neg_iterate_pairs = 0` everywhere I looked:
- on all 22 builders;
- on every step of the 5 series: 33 + 193 + 24 + 62 + 343 pairs seen by `punch_holes`.

The iterates recorded are 0 and positive only. For example, `p3_series` runs 0..16 and `nested_outer_series` runs 0..11.

Over the same runs:
- `seen` and `punched_regions` rejected **0** holes;
- `_punch_in_bridge` returned None **0** times;
- **0** pairs were skipped for a missing `arc_point`;
- every origin has exactly one reference (`ref_per_origin_max = 1`).

**Why** (believed, consistent with A.0): `extend_pseudoneighbor_trajectories` does walk `direction = -1` (`Pseudoneighbor.py:256-283`), but a reference lies in the window `S \ f(S)`. One of its two preimages therefore leaves T, and the walk breaks at the first step. The step-1 assertion makes this falsifiable, as the plan intends.

**The fallback is not quite never entered.** On `near_tangency` (a non-law builder), 14 references have no spanning bridge. The fallback runs, finds no punched pair, and logs **14 "Orbit … cannot propagate" WARNINGs**; it creates no holes. Deleting it removes those 14 WARNINGs. `punch_holes`' own "Skipped 14 pair(s)" WARNING still names the pairs. That is a log change on one non-law builder; no test reads it (`test_resolution_warning_follows_oracle` listens only to the workbench logger).

## 3. The placement decision (verified, measured on the 11 law cases and more)

What was measured:
- **Direct holes:** the reference placement recomputed from `hole.bounding_ids` (already in near/far order) against `Hole.coords`.
- **Propagated holes, (a):** the recomputed reference point is carried back `-iterate` times by `trellis.dynamical_system.map_inv`, using `np.asarray(..., float64).ravel()` per step exactly as `_map_backward` does. The span is the reference bridge's ids moved by `image_cdist`, and then `_punch_in_bridge`'s nudge is applied.
- **Propagated holes, (b):** the no-map fallback.
- **Side:** read by `diagnostics._bridge_side_of(trellis, container, point, span)`.

| case | holes (direct / propagated) | max \|Δ\| (a) | max \|Δ\| (b) | sub-arc side changes (a) | sub-arc side changes (b) | wrong side vs `bridge_side`: today / (b) |
|---|---|---|---|---|---|---|
| k10 | 1 / 2 | 0.0 | 0.032 | 0 | 1 | 0 / 1 |
| k28_one_blast | 2 / 2 | 0.0 | **0.328** | 0 | 1 | 0 / 1 |
| k28_two_blasts | 1 / 3 | 0.0 | 0.178 | 0 | 1 | 0 / 1 |
| p3 | 2 / 6 | 0.0 | 0.010 | 0 | 3 | 0 / 3 |
| nested | 3 / 9 | 0.0 | 0.099 | 0 | 4 | 0 / 4 |
| nested_inner_blasts | 3 / 9 | 0.0 | 0.099 | 0 | 4 | 0 / 4 |
| inversion | 2 / 4 | 0.0 | 0.057 | 0 | 2 | 0 / 2 |
| inversion_one_blast | 1 / 6 | 0.0 | 0.038 | 0 | 1 | 0 / 1 |
| inversion_two_blasts | 2 / 6 | 0.0 | 0.038 | 0 | 1 | 0 / 1 |
| inversion_period2 | 2 / 8 | 0.0 | 0.004 | 0 | 1 | 4 / 5 |
| heteroclinic | 4 / 69 | 0.0 | 0.099 | 0 | **30** | 13 / 29 |
| p3_coarse_blasts | 4 / 80 | 0.0 | 0.002 | 0 | 24 | 43 / 59 |

- **Direct holes:** the reference placement is also exact (deviation 0.0 on every builder).
- **(a) everywhere:** it reproduces today exactly on all 22 builders, including heteroclinic_17u, the p3 deep builds and k28_deep_blast_sliver.
- **No map available:** `hole_point` must keep the fallback branch for `trellis.dynamical_system is None`. That is today's behaviour exactly (`_map_backward` returns None, which leads to `_step_into`), so it is not a choice to make.
- **"Figures: identical" holds only under (a).** Under (b), propagated holes move on every law case.
- **The B4a numbers are reproduced:** heteroclinic 29 wrong with (b); this table adds the 30 side changes relative to today. One correction: the request's "10 of the 11 law cases" is **11 of 11** by side change.

## 4. What `diagnostics.measured_sides` needs once `coords` goes (verified)

`measured_sides` (`diagnostics.py:429-482`) reads `reference.coords` at :469. That is the reference placement only; it never needs the nudged propagated point. The measurement itself (carry back, then `_bridge_side_of` on the sub-arc) does not change.

Recomputing the reference point gives the same values exactly (0.0 deviation above), so every measured side stays identical.

**Smallest fix.** Line :469 becomes the reference placement from `plotting.hole_point(trellis, reference)`, imported inside the function. `plotting` imports only `TopologyResults` at module level, and `topology/__init__` already imports both modules, so there is no cycle and no new import cost.

The docstring note at :455-457 (`StablePartition._map_backward, until B4b`) also goes. So **`diagnostics.py` joins the file list** (about 3 lines).

**Its tests do not join.** `test_production_checks.py` drives `measured_sides` on a real trellis, and `test_topology_invariants.py:174` monkeypatches it. Both are unaffected. `test_topology_invariants.py` joins only for its `Hole(coords=)` (see §6).

## 5. Inversion-anchor identification counts (verified)

| case | anchor ends among hole bridges | open on `s(o,1-b)` | open on `s(o,b)` | open **nothing** |
|---|---|---|---|---|
| inversion | 4 | **2** | 2 | 0 |
| inversion_one_blast | 2 | **2** | 0 | 0 |
| inversion_two_blasts | 2 | **2** | 0 | 0 |
| inversion_period2 | 8 | **4** | 4 | 0 |
| k10 | 1 | 0 | 0 | **1** |
| k28_one_blast | 2 | 0 | 1 | **1** |
| k28_two_blasts | 1 | 0 | 0 | **1** |
| p3 | 6 | 0 | 3 | **3** |
| nested / nested_inner_blasts | 7 / 7 | 0 | 3 / 3 | **4 / 4** |
| heteroclinic | 10 | 0 | 6 | **4** |

The plan's inversion counts (2/4, 2/2, 2/2, 4/8) are **confirmed**.

The plan's "non-inversion cases take `s(o,b)` at every anchor end" is **imprecise**. On every non-inversion case, roughly half the anchor ends get the D6 verdict "anchorward". The identification then names `s(o,1-b)`, which is not a branch (a point seeded on one eigendirection). So today's `key in trellis.branches` guard (:780) **records no opening** there.

So step 3's rewrite must keep "opens nothing when the identified branch is absent". It is today's code reworded: name the branch with `stable_key_at`; when the opening would be anchorward, take `s(o,1-b)` and open its outward half; skip if that branch is absent. With that it is behaviour-identical.

## 6. What differs from the plan section

1. **Step 2's WARNING fires on law cases (rubric 5 fails).**
   - Span ends that come from the iterate table are **exactly** a crossing of the container's branch (distance 0). Each is the container's own end crossing, which is in N.
   - This happens on **6 of the 11 law cases** (k28_one_blast 1, p3 3, nested 3, nested_inner_blasts 3, inversion 2, inversion_period2 4) and in every series.
   - Examples: k28_one_blast `iterate(5, -1) = 3`, container `(0, 3)`; inversion_period2 reference `(14, 15)`, where 15 maps to 5, 9, 2 and 13, each the end of its container.
   - Restricted to **scaled** ends (`from_table` False), the WARNING fires on **no** law case. The nearest scaled end sits at a relative 0.0199 (inversion blasts) and 0.022 (k28_two_blasts), against `SCALING_RTOL = 1e-2`: a factor-2 margin. It does fire where NUM-5 lives: `k28_deep_blast_sliver` (18 ends) and `k32_series` (242 ends, 5 per-branch firings, minimum relative distance 1.2e-4). No test reads topology WARNINGs on those.
   - So step 2 must keep the flag: `key, c, from_table = trellis.image_cdist(...)`. Part (iii) becomes "**use** the flag for the WARNING", not "drop it".
   - Rubric 5's unit test must build its near-crossing end as a **scaled** end.
2. **Step 2's comment is false as worded.** "Never in the table, by construction" fails for one end, exactly as INT-2 says ("at most one of its two ends a registered crossing"). A correct one-liner: *"the span's interior is never a crossing of N (it lies past the pip cut); one end may be a registered crossing, which is then the container's own end."*
3. **`_near_far` stays.**
   - It orders `bounding_ids` and `near_intersection_id` at :340 and :985. Both are hole fields: `bounding_ids` is in the law fingerprint, and `near_intersection_id` is read by `MinimalTrellis.py:303`, `DualGraph.py:1150` and `laws.check_holes_are_classified`.
   - Step 1's "delete `_near_far` where only placement used it" therefore deletes nothing.
   - The `_near_far` tests (`test_stable_partition.py:529-549`) stay put.
   - `hole_point` does not need `_near_far`, because `bounding_ids` is already in near/far order.
4. **The file list misses four files:**
   - **`tests/helpers/laws.py`.** Three laws read `hole.coords`:
     - `check_holes_are_classified` :1213 (asserts finite coords);
     - `check_backward_holes_keep_bridge_side` :1254 (carries `origin.coords` back by `map_inv`);
     - `check_direct_hole_side_is_coordinate_side` :1280.

     The last two are geometric oracles on the reference placement. Each becomes `plotting.hole_point(trellis, hole)` (values identical, §3). The finite-coords clause either uses `hole_point` or is dropped.
   - **`src/tanglepack/topology/diagnostics.py`** (§4).
   - **`tests/unit/topology/test_topology_invariants.py`** (:68, `_hole` builds `Hole(coords=...)`).
   - **`tests/unit/topology/test_minimal_trellis.py`** (:71, `ghost = Hole(coords=...)`).

   Also inside the listed files:
   - `test_stable_partition.py:82`, `:392` (`Hole(coords=)`) and `:355-373` (reads propagated `hole.coords`, which becomes `hole_point`);
   - the `TangleSession` docstring line is now **:1321**, not :1312.

   That makes **9 files** plus the one docstring line.
5. **Rubric 2's ≤ 700 lines is unreachable** (estimate, not built).
   - B4b deletes about:
     - 145 lines of helpers (`_step_into`, `_stable_arc_midpoint`, `_nearest_arc_point`, `_map_backward`, `_image_arc_midpoint`);
     - about 43 of `_punch_in_bridge`;
     - about 25 of `punch_holes` placement;
     - about 50 of the propagation dedups and fallback (net of the new WARNING);
     - 12 for `_region_key`;
     - about 20 of Dev Notes;
     - the numpy imports.
   - That takes the file from 1,354 to **about 1,060-1,090** lines. The B4a analyst estimated about 1,090.
   - Reaching 700 would mean moving the 216-line Dev Notes or the partition code, which is out of B4b's scope.
   - Rubric 2's grep is achievable once :1090's log is changed. No other `np.` use survives, so `numpy` and `NDArray` drop out of the module.
6. **Smaller points:**
   - The plan's step 4 describes propagated placement as "nudged toward the chord". That is the fallback; today's placement nudges toward the carried point (§3).
   - Step 1 says "punch the REFERENCE pairs only (`iterate == 0`)". Today every punched pair already has iterate 0, and `is_reference` holds exactly when `iterate == 0` (`Pseudoneighbor.py:219, :279`), so this is the same set.

## 7. Recommendation (the smallest change that makes the slice hold)

1. **Placement: default to (a), the carried placement in `plotting.hole_point`.**
   - Reference holes: `_stable_arc_midpoint` + `_nearest_arc_point` + `_step_into` on `hole.bounding_ids`.
   - Propagated holes: find the origin's reference hole in `trellis.holes`, compute its point, carry it back `-iterate` times by `trellis.dynamical_system.map_inv`, take the span from the reference bridge's ids by `image_cdist`, and apply `_punch_in_bridge`'s nudge.
   - Fall back to `_step_into(midpoint, stable_mid, 0.05, arc_chord)` when there is no map, as today.
   - This meets rubric 6 exactly (0.0) and "Figures: identical". The algorithm module then calls no map.
   - The planner puts (b) to the author as the alternative: a stated figure change of up to 0.33, with propagated holes on the other side of their sub-arc on all 11 law cases.
2. **Step 2 WARNING on scaled ends only:** `if not from_table and any(abs(c - x) <= SCALING_RTOL * abs(x) for x in branch_cdists)`, once per branch. Reword the comment as in §6.2.
3. **Keep `_near_far`** in StablePartition and drop it from the deletions. Delete `_region_key` with the dedups. Change the :1090 debug log to name `hole.bounding_ids`.
4. **Step 3:** state the identification as in §5, keeping the "branch absent → no opening" guard.
5. **File list:** add `tests/helpers/laws.py`, `topology/diagnostics.py`, `test_topology_invariants.py` and `test_minimal_trellis.py`. `measured_sides` and the three laws read the reference point via `plotting.hole_point`.
6. **Rubric 2:** change the line target to "≤ 1,100" (or "−260 or more"); keep the grep.

## Out of scope (noticed, not for the recommendation)

- The `direction = -1` loop in `Pseudoneighbor.extend_pseudoneighbor_trajectories` (:256-283) produces no pair on any build. It is the upstream of the dead direct-backward path, and B4b's assertion sits downstream of it.
- Today's placement (and therefore (a)) puts propagated holes on the wrong side of their sub-arc on heteroclinic (13/69), inversion_period2 (4/8), p3_coarse_blasts (43/80), p3_16_12 (14/40) and inversion_period2_10 (9/12). This affects plots only, never the topology.
- The I1/I2 cores still live in `diagnostics` (B4a's open point).
- `CLAUDE.md` still describes `coords`, `measured_side` and `side_source` (B12).

## Verified vs believed

- **Verified (probe, Read, grep, `git log -S`):** everything in §§1-6 except the two items below.
- **Verified baseline:** 140 tests passed at `HEAD` in the touched files (`test_stable_partition`, `test_topology_invariants`, `test_minimal_trellis`, `test_production_checks`, `test_plot_topology`, `test_plot_smoke`). I did not run the full suite.
- **Believed:**
  - the line estimate in §6.5 (arithmetic over spans, not a built file);
  - the reason no backward pair exists (§2, consistent with A.0 and with the probe).

## Files

- Knowledge files updated:
  - `.refactor/knowledge/topology.md` (header; a new B4b section)
  - `.refactor/knowledge/testing.md` (a new B4b section)
- Probe script and outputs (session scratchpad, not in the repo):
  - `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b4b_probe.py`
  - `.../scratchpad/builders.out`, `.../scratchpad/series.out`
