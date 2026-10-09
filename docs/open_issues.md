# Open issues — tanglepack (catalog of 2026-10-07)

This catalog is the output of a read-only audit workflow (2026-10-07). It harvested the docs, the `.refactor` notes, the stress report (second edition), the source Dev Notes, the tests and the memory. Three agents pushed inversion, period-k and nested cases through the whole pipeline and reproduced the bridge-class growth, and a chalkboard review compared the code with the mathematics (see [`pipeline_chalkboard.md`](pipeline_chalkboard.md)). Line numbers are as of main `684d70b`.

Every issue carries the work package (WP) that owns it; the plan is in the "Work packages" section at the end. Author decisions taken on 2026-10-07 are recorded in the "Author decisions" section and in `.refactor/decisions.md`. When an issue is closed, change its status to `closed (<commit>)` and leave the entry in place.

## Summary

| Theme | High | Medium | Low | Total |
|---|---|---|---|---|
| Inversion end-to-end (INV) | 6 | 4 | 1 | 11 |
| Period-k orbits (PK) | 2 | 2 | 5 | 9 |
| Nested / heteroclinic (NEST) | 2 | 3 | 1 | 6 |
| Bridge-class growth (BCG) | 4 | 3 | 0 | 7 |
| Intent vs code (INT) | 1 | 8 | 2 | 11 |
| Unnecessary complexity (SIMP) | 0 | 1 | 12 | 13 |
| Resolution & numerics (NUM) | 0 | 11 | 10 | 21 |
| Performance (PERF) | 1 | 1 | 3 | 5 |
| Tests & coverage (TEST) | 0 | 2 | 5 | 7 |
| API traps (API) | 0 | 2 | 5 | 7 |
| Docs drift (DOC) | 0 | 0 | 5 | 5 |
| **Total** | **16** | **37** | **49** | **102** |

The audit raised 18 decision questions (DEC-1..18); all are answered or deferred below. Five issues are marked "possibly stale — verify": INV-10, PK-7, NUM-19, DOC-3 and DOC-4.

The most serious situations:
1. **Class count grows by bookkeeping, not dynamics (BCG-1/2/3).** Every bridge in the trellis is classed, including blast children with no image. Loops are folded through preimage links that blasts never register. Result: +1 active and +1 inert class per blast, while holes, elements and the minimal trellis stay constant.
2. **The hole set depends on depth (BCG-4, PK-1).** Deep tip or blast-child crossings land in the fundamental window and replace the primary pseudoneighbor pairs. The new reference punches about log(c_u)/log β holes, so p3 goes 12→42 holes at blast 6 or at 15 steps. That makes p3 unreliable: an unreachable class and a virtual `new1`. **WP8:** every switch is a definitional discovery with a witness the iterate table names; the losses driven by cdist-scaling noise are fixed; the decreases left (witnesses inside the placement error bar, noise not ruled out, but one) are DEC-22. **WP8b:** the decreases DEC-22 left are gone: every witness lay deeper than the deepest registered point of the pip's own orbit, and the trellis is now cut there (D12, the trellis of the pip); the class count is non-decreasing on every series and E11 replay measured.
3. **Inversion is wrong by construction, and no fixture can show it (INV-1/2/4/5/6).** Four anchors sit at one point. The pseudoneighbor lap factor is λ_u instead of λ_u². The zone ring is keyed by orbit index. The side and partition model is per branch. The only fixture has no holes.
4. **The heteroclinic path has never run (NEST-1/2).** No fixture has an A(u)×B(s) crossing. A hand-built case is unreliable, and heteroclinic pseudoneighbors are rejected by construction. **WP9:** the law case `heteroclinic` runs end to end with heteroclinic pairs holing both partitions (D7); it is unreliable only through DEC-23 (equal shortest walks around the lobe). **WP9b:** the non-crossing rule (D13) picks the side of every class with a member lobe; the one tie left is the member-less connecting class `bn`, which no closed curve with a known image goes through (NUM-21, DEC-24). **WP9c:** `bn` is no construction error: it is `f^2` of the computed connecting bridges `(91, 92)` / `(93, 94)` (class `bi`), member-less only by D12's cut, and the rule that decides it is D13's own wording, walls (C17): the decided image walks of arcs on other unstable branches, which `f(b)` cannot cross. A heteroclinic bridge is a bridge like any other (D14); `heteroclinic` is reliable, and no tie is left on any law case, deep run or series step.
5. **The dual-graph walk fails where it must succeed (INT-1, PK-1).** Walks come back "unreachable" or "ambiguous" and virtual `new*` classes appear. Under D3 the walk IS the word and must always exist, so each of these is a construction bug (WP7). **WP7:** the virtual classes were the alphabet not being closed under the word map. It is now closed over the classes registered bridges span, and a walk pair NO bridge spans is a construction error: its class is unresolved with a WARNING and the dynamics is not reliable, never a lettered class. The deep runs (p3 15/16 steps, 6 blasts, nested with 6 inner blasts, k=2.8 at 8 blasts) still each name one or two such pairs; whether they must get their bridges is DEC-21. The one unreachable walk left (inversion x period 2 at 10 steps, not a law case) needs the stable manifold cut at the pip orbit (DEC-20). **WP8b:** with S cut at the pip orbit (D10) inversion x period 2 at 10 steps walks every class, and a walk pair no computed bridge spans is a class flagged "no computed bridge yet" (D11): the deep p3 and nested runs are reliable; k=2.8 at 8 blasts is not (two ambiguous walks, one evidence disagreement: WP10). **WP9b:** those ties were walks starting on edges the landing span does not meet; walking between the spans, k=2.8 at 8 blasts is reliable.

---

## Inversion end-to-end

### INV-1: Four coincident (0,0) anchors at an inversion point instead of one degree-4 node
**Status:** closed (WP3) | **Severity:** high | **Related:** INV-3, INV-6, DEC-4, API-5, SIMP-12 | **WP:** WP3

**WP3:** a periodic point is ONE anchor, registered under its home pair `(u(o,0), s(o,0))` with that pair's sign and put first on every branch of its point by the trellis and the registry's graph view; at exactly `(0, 0)` the registry matches any two keys of the same point, so the root segments of the other branch pairs resolve to it. `create_bridges` cuts every unstable branch rooted at the point at its root (`Tangle._both_anchors` is gone), and the arrangement draws one degree-4 node (branch 1 in the `u-`/`s-` slots). Bridge ends, image chains and the iterated cut ask `Trellis.stable_key_at` which branch an arc meets the anchor on (C7). Inversion: 1 anchor (was 4), no virtual class and `is_reliable` at 0/1/2 blasts, classes 4 / 10 / 10 (was 6 / 11 / 11), holes 1 / 9 / 9. Two of its old strict xfails had other roots and are re-tagged: `anchor_bridge_class_leads_its_tangle` x `inversion` (PK-3, WP5) and `arrangement_euler` x `inversion_two_blasts` (INT-9, WP7).

`_register_anchors` registers one anchor per (unstable branch, stable branch) pair. The only filter is on orbit index, so an inversion point gets two anchors per unstable branch, four in all. The author's rule is one.

This cascades:
- 5 strict xfails: `one_anchor_per_unstable_branch`, `anchor_bridge_class_leads_its_tangle`, `arrangement_euler` (V−E+F=0), `itinerary_pairs_same_side`, `landings_contained`.
- `_image_chain` collects both anchors, which creates the phantom pair (0,1). On the k=4 probe this gives `new1 = {L_(0.0;1), L_(0.1;1)}`, and anchor bridges `b`/`h` are split by which anchor id they got.
- The singleton `[0,0]` owns two crossings, and the dual graph builds a zero-length stable edge.

**Evidence:** `src/tanglepack/numerics/TangleWorkbench.py:1537-1616` (:1574-1599, filter :1598, docstring :1575-1580); `src/tanglepack/numerics/Tangle.py:9-16`; `src/tanglepack/numerics/Tangle.py:1050-1072` (`_both_anchors`); `src/tanglepack/topology/Arrangement.py:80-88`; `src/tanglepack/topology/BridgeClass.py:748-758`; `tests/cases.py:796-815`; `tests/conftest.py:99-104`; `docs/regions_bridge_identity_and_cleanup_plan.md:334-337`; `stress/numerics/report/report.tex:2012-2014`; probe `inv_probe.py` (zone+1 blast: `new1 b u^-1 c^-1` vs walk `h u^-1 c^-1`)

### INV-2: Pseudoneighbor uses bare λ_u as the branch-return factor (should be β^k_value = λ_u²)
**Status:** closed (WP2) | **Severity:** high | **Related:** SIMP-3, INV-4 | **WP:** WP2

- `_resolve_reference_end` scales by `/lambda_u` and `*lambda_u`.
- `_unstable_interval_is_clear` uses `c*lambda_u**m` per k-step lap, with `log(lambda_u)` for the m-range.

A lap is k_value = 2·period steps, which is β^k = λ_u². Patching this in memory changes the pairs:

| Variant | Code as is | With β^k |
|---|---|---|
| No zone | `[(23,4)]`, 7 holes | `[(11,9),(19,4),(25,17)]`, 16 holes |
| Zone | `(8,5)` | `(8,7)` |

The current fixture never reaches this code because it has no pairs. `TrellisBranch.eigenvalue` and `IntersectionRegistry.on_interval` have the same defect.

**WP2:** the all-iterates check laps by `per_step_beta("unstable") ** k_value`, and r_{n+p} is the table's `M^k(r_n)` (the λ_u cdist fallback is gone). The inversion law cases have 1 / 9 / 9 holes at 0 / 1 / 2 blasts.

**Evidence:** `src/tanglepack/topology/Pseudoneighbor.py:51-60` (OPEN Dev Note), `:491-492`, `:551`, `:571-576`; `src/tanglepack/topology/TrellisBranch.py:79`; `src/tanglepack/numerics/IntersectionRegistry.py:656`; `docs/regions_bridge_identity_and_cleanup_plan.md:330`; probe `inv_probe.py`

### INV-3: An anchor-only stable branch at 4-5 steps raises AssertionError "owned by 0 partition elements" (WP4-B)
**Status:** closed (WP3) | **Severity:** high | **Related:** INV-1, DEC-5 | **WP:** WP3

**WP3:** C1: a stable branch with no crossing past its anchor is one element `[anchor, end]` (closed at the anchor, `hi_id None`, to the end of the computed branch), logged at INFO on the `StablePartition` logger. `build_k10_inversion(unstable_steps=4|5)` runs through `symbolic_dynamics` (regression `test_inversion_before_the_tangle_forms`).

`_build_intervals` gets a single boundary, so `_element_of_intersection` asserts. This reproduces at k=10 with 4-5 steps, and at k=3.5 and k=4 with 4-8 steps. There are two options: a closed singleton at the anchor, or one element `[anchor, end]`. Either one has to come after INV-1 and dedupe by point. Once decided, add `test_inversion_before_the_tangle_forms`.

**Evidence:** `src/tanglepack/topology/StablePartition.py:214-235` (DEFERRED Dev Note), `:662`→`:2147`; `stress/numerics/report/report.tex:941, 956-967, 2009-2014`; `stress/numerics/report/numbers_extra.tex:297`; `stress/numerics/e07_breaking.py:312-318`; probe `inv.py` (steps 4, 5 ERR)

### INV-4: No usable inversion fixture: the k=10 saddle escapes, has no pseudoneighbors and no holes, so most laws are N/A
**Status:** closed (WP1) | **Severity:** high | **Related:** INV-2, INV-11, TEST-1, DEC-4 | **WP:** WP1

**WP1:** the k=4 other-root saddle is the law case `inversion` at 0, 1 and 2 blasts (`tests/cases.py` `build_inversion`); the k=10 saddle survives as the numerics fixture `build_k10_inversion`. Holes were 1 / 7 / 0 (the 0 was INV-2); they are 1 / 9 / 9 since WP2, and the vacuous-law xfails are gone.

`build_inversion` grows k=10 to 6/5 steps:
- 7 steps raises PointBudgetExceeded, and points jump from 196 to 49,417 between steps 5 and 6. The window between "too shallow" and "escaped" is one step wide.
- At 6 steps: 0 holes, 0 kept minimal-trellis bridges, 4 classes, `new1`, `is_reliable=False`.
- 15 or more laws are NOT_APPLICABLE on inversion, so StrongPip, Pseudoneighbor, StablePartition, ResonanceZone and Arrangement all say "inversion unvalidated".
- No E11 case and no E6 sweep point gives a nontrivial inversion tangle. In the sweep, k≤5.70 has 0 crossings, 6.92 has 1, 8.40 has 7, and k≥10.19 hit OOM.

**Candidate fixture (Case B):** the k=4 other root `(1-√5, √5-1)`, λ -1.96/-0.51, `area_cutoff=1e-5`, 11/11 steps or `grow_until_bounded`, a zone at f(q0), 0/1/2 blasts. Laws pass 70/75, 71/75 and 70/75; step 12 escapes.

**Evidence:** `tests/cases.py:568-636`, `:843-875` (`_NO_HOLES`); `src/tanglepack/topology/StrongPip.py:63-72`; `src/tanglepack/topology/Pseudoneighbor.py:37-40`; `src/tanglepack/topology/StablePartition.py:236-238`; `src/tanglepack/loom/ResonanceZone.py:37-39`; `stress/numerics/report/report.tex:210, 830-832, 940-991, 1150-1156`; `stress/numerics/report/numbers_extra.tex:219-221, 287, 301`; `stress/numerics/results/e06_parameter_sweep.jsonl` (which='other'); probes `inv.py`, `inv_probe.py`

### INV-5: Resonance-zone ring is matched by orbit index only, so the inversion zone polygon is wrong; the logic is duplicated in TangleSession
**Status:** closed (WP3) | **Severity:** high | **Related:** INV-4, DEC-4, NEST-4 | **WP:** WP3

**WP3:** `ResonanceZone.ring_cycles` follows a cut's stable branch `s(o,b)` in and `u(o,b)` out (the FULL key); `_build_boundary` and `TangleSession._zone_ring_order` both use it. The inversion zone is the union of its two loops through the anchor (area 1.412087 = 0.706 + 0.706, was 1.232); p3, k28 and nested zone areas are unchanged to 12 digits. One zone per fixed point (C8): `ResonanceZone.key` is the fixed point, `branch_index` is gone, `zone_label` reads "zone p<period>". The zone changes which bridges the inversion blasts iterate (classes 4 / 10 / 10).

`_build_boundary` builds `by_u_orbit = {manifold_a_key[2]: i}`. Both inversion cuts have orbit 0, so u0.1 is used twice and u0.0[p,q0] is never used.

- The probe zone area is 1.232. The correct single-pip loops are 0.705 and 0.704 (sum 1.41).
- The zone's own boundary bridge (1,28) is classified "no zone", so the choice of blast interior bridges is suspect.
- Zones are keyed `(fixed_point, branch_index)`, which gives "two zones" for one chain.
- `TangleSession._zone_ring_order` re-implements the same hop.

**Evidence:** `src/tanglepack/loom/ResonanceZone.py:223-252`, `:121-133`, `:37-39`; `src/tanglepack/loom/TangleSession.py:2052-2076`, `:155`; probe `inv_probe.py`

### INV-6: The same-side law and the partition are per branch, while L(0.0) and R(0.1) are one global half-plane through p
**Status:** closed (WP3) | **Severity:** high | **Related:** INV-1, DEC-4, INT-8 | **WP:** WP3

**WP3:** `cross_side` is set only for a pair on one branch (D4). The elements next to p are NOT glued across branches: every branch keeps its own two sides (C7), and the zone boundary through p separates them. `landings_contained` passes on every inversion case.

Both inversion stable branches are oriented anchorward, so they point in opposite directions along one curve. Several pieces have no inversion-aware side map:
- `cross_side = first.side != second.side`
- the `itinerary_pairs_same_side` law
- `oriented_class`
- the cross-branch row checks

As a result, legitimate `(R(0.1), L(0.0))` pairs are flagged as cross-side, and the WARNING "class a: word b^-1 disagrees with the registered image" appears. The elements next to p on one global side (L(0.0;1) and R(0.1;1)) are not glued. Landing spans also straddle cuts (`contained=False`).

**Evidence:** `src/tanglepack/topology/SymbolicDynamics.py:474`; `tests/helpers/laws.py:1885`; `tests/cases.py:808-815`; `CLAUDE.md` "a bridge only connects same-side elements"; probes `inv.py`, `inv_probe.py`

### INV-7: `branch_index` means two things (inversion half-chain vs an independent second chain); the initializer duplicates branch ordering
**Status:** closed (WP3) | **Severity:** medium | **Related:** INV-9, SIMP-9 | **WP:** WP3

**WP3:** `construct_kevin_way` is one loop: a chain of `k_value` pieces per requested eigendirection, walked by `advance_key`, skipping a direction an earlier chain covers (branch 1 of an inversion point). `FixedPoint.get_iterable_array` and `get_branch_array` are deleted (`test_num_branches_agrees_with_get_branch_array` with them). A branch index names an eigendirection. Behaviour-preserving: `test_kevin_way_seeds_one_chain_per_eigendirection` passes on main too.

- For inversion, branch 1 is the second half of one k_value chain.
- Otherwise it is a separate chain seeded with `num_branches=2`.

Both use ManifoldKey slot 3, so every consumer special-cases `check_inversion`, and `num_branches` means different counts in FixedPoint and in the initializer. `construct_kevin_way` still uses the old `get_iterable_array` ordering ("0, p−1, …, 1") beside `advance_key`, with two near-identical loops.

**Evidence:** `src/tanglepack/numerics/FixedPoint.py:151-183, 453-456, 484-492, 497`; `src/tanglepack/numerics/ManifoldMachine.py:286-298`; `src/tanglepack/numerics/ManifoldInitializer.py:25-30, 153, 432, 483-598, 527, 530-598`; `src/tanglepack/numerics/TangleWorkbench.py:352, 521-548, 810, 1000`

### INV-8: `orient_eigenvectors` is a no-op for inversion; per-orbit-point eigenvector propagation is skipped for period>1 inversion orbits
**Status:** closed (WP3) | **Severity:** medium | **Related:** INV-7, TEST-1 | **WP:** WP3

**WP3:** the early return is gone: orienting an inversion point chooses which half of each manifold is branch 0, propagated along the orbit as without inversion. No law case orients an inversion point, so no case output moves.

Period-1 inversion is fine. For a period>1 inversion orbit, the orientation along the orbit is never set.

**Evidence:** `docs/regions_bridge_identity_and_cleanup_plan.md:331-332`; `src/tanglepack/numerics/TangleWorkbench.py:~529-536`; `src/tanglepack/numerics/FixedPoint.py:10-21` (inversion eigenvector storage Dev Note)

### INV-9: `grow_n_times(branch_index=0)` default logs "branch_index is ignored" on every inversion growth step
**Status:** closed (WP3) | **Severity:** low | **Related:** INV-7 | **WP:** WP3

**WP3:** the machine walks the one inversion chain from branch 0 whatever is asked, without a warning; the workbench's `growth_branch` special cases are gone. `build_inversion(blasts=1)` logs no WARNING from `ManifoldMachine` (was 22).

This is log noise; the probe saw it three times. Fix by defaulting to `branch_index=None` or skipping the warning for inversion.

**Evidence:** `src/tanglepack/numerics/TangleWorkbench.py:590`; `src/tanglepack/numerics/ManifoldMachine.py:286-289`

### INV-10: `iterate_bridge` on an anchor-rooted inversion bridge raised TypeError in `_insert_after` (branch_index None)
**Status:** closed (WP3) | **Severity:** medium | **Related:** INV-1 | **WP:** WP3

**WP3:** verified: no TypeError on main, but `iterate_bridge` on an anchor-rooted inversion bridge returned ONE PARTIAL child (its image starts at an anchor the image's crossing set did not name). With one anchor every child is identified and the first starts at the anchor (k10: `(0,1) -> [(0,7)]`, `(0,7) -> [(0,1), (1,3), (3,2)]`); regression `test_iterating_an_anchor_rooted_inversion_bridge_cuts_identified_children`.

`create_bridges` now refuses the zero-span cut between two anchors (inversion bridges went 9→7). Nobody re-checked whether a non-degenerate anchor-rooted bridge still fails.

**Evidence:** `docs/regions_bridge_identity_and_cleanup_plan.md:196, 332-334`

### INV-11: Orientation-reversing maps (det J < 0) are rejected, so every side-flip path is unexercised
**Status:** closed (WP11) | **Severity:** medium | **Related:** SIMP-1, DEC-10, TEST-1 | **WP:** WP11

`set_k_value` rejects disagreeing eigenvalue signs, so all 8 b=-1 configs are rejected and there is a strict xfail with `raises=ValueError`. The parity flip on carried hole sides (WP7-D2), `check_holes_share_bridge_side(orientation_preserving=False)`, `land_element` and `image_of_element` flips have never run on real data.

Note that inversion has det J = +1. Earlier notes conflate "inversion" with "orientation-reversing", so the plan must say which one is meant. Supporting det<0 would need a per-stability k_value.

**WP11:** closed by C3: det J < 0 stays out of scope and `orientation_reversing_case_builds` stays a strict xfail. The flip rule is ONE function, `topology.side_after(side, steps, preserving)`, called by `land_element`, the dual-graph markers, backward hole propagation, `_orbit_side` and `side_diagnoses`, and pinned by `test_side_after.py`. The orientation is read off the FixedPoint (C16): `preserves_orientation` is the sign of the solver's single-step Jacobian (any period), and `Trellis.orientation_preserving` no longer calls the map's Jacobian. A period-2 orbit of b = -1 passes `set_k_value` (det DM^2 > 0) and could exercise the rule end to end if the author ever admits det J < 0.

**Evidence:** `src/tanglepack/numerics/FixedPoint.py:185-225`; `src/tanglepack/numerics/FixedPointSolver.py:16-22`; `tests/cases.py:36-39, 826-830`; `tests/conftest.py:108-110`; `src/tanglepack/topology/Trellis.py:32-40`; `stress/numerics/report/report.tex:833-834, 1906`; `docs/dual_graph_symbolic_dynamics_plan.md:338`

---

## Period-k orbits

### PK-1: p3 at ≥15 unstable steps or ≥6 blasts is unreliable: an unreachable class, a virtual `new1`, and "crossing abuts no hole"
**Status:** closed (WP8b: every class walks, virtual classes deleted (WP7); a walk pair no computed bridge spans is a flagged class, D11) | **Severity:** high | **Related:** BCG-4, PK-3, PK-5, NEST-3, DEC-1, DEC-21 | **WP:** WP7, WP8b

At 11-14 steps and 0-5 blasts: 9 classes, 3 active, reliable. At 15 steps there are four tip crossings with neither a −1 nor a +1 link inside the ray-1.0 window:
- They break the adjacency of `(18,34)`.
- They disqualify `(4,20)` through `M^-17` of tip crossing 24.
- The only references left are the tip pairs `(31,9),(9,15)`.

The anchor image chain `0→pip→4→20` then finds no hole at 20/4, so `PartitionFamily.py:800` warns "abuts no hole" and does not cut. That produces wrong elements, then the unreachable walk `0.0/R#2→R#3`, then `new1`.

E11 has exactly 1 unresolved and 1 virtual class at every blast from 6 to 20 (p3 and nested_inner); blasting never resolves them. After the bursts, every unresolved class is "dual-graph walk unreachable" (127 entries, 0 singleton). k28 ends with 12 unresolved, 1 virtual and 3 ambiguous; k32 with 9/0/2. The report suggests stable growth instead of blasting (untested).

**WP4:** the unreachable class `t` was the geometric hole openings: with the signs deciding them, the deep runs gain 3 elements (81 -> 84, 89 -> 92, 81 -> 84) and every class walks (new test `test_deep_p3_every_class_walks`); the virtual `new1` remains (strict xfail `open_p3_deep_runs`, WP7). Replayed E11 p3 / nested_inner at 1e-5: 0 unresolved / 1 virtual at blasts 6-8, then 3 unresolved / 2 virtual from blast 9 (WP7). Inversion x period 2 at 10/10 steps no longer fails `openings_*` or the library row check; `backward_holes_keep_bridge_side` (its oracle re-measures a carried point that escapes the 3-node reference lobe), `bridge_endpoints_on_own_branch` (WP7-D1) and the symbolic laws still fail there.

**WP6 (for WP7): the alphabet is not closed under the word map.** In deep runs an image-bridge class that is not a hole-bridge class has a word naming the class of the image's image, which the minimal trellis does not keep: p3 from blast 7 (`u = {(40,41), (44,45)}` maps to `w = {(46,47), (50,51)}`), k28 at blast 8 (`e` maps to `g`). Those were lettered blast children before WP6 and are virtual classes now (p3 deep 1 `new1`, k28 blast 8 2, was 3). Rooting it out (keep the images of active image-bridge classes, or show the walk never needs them) changes the minimal trellis's definition: D3, WP7.

**WP7: root cause and fix.** Every walk exists on the deep runs since WP4; the only failing symbolic law was "no virtual class". The walk does not read the wrong elements: the virtual pair is the class of the IMAGE of an image-bridge class, a piece of `f(B)` that the minimal trellis (hole bridges + the FIRST images of active hole bridges) does not keep and that, on most runs, no computed bridge spans. On `p3_6_blasts` class `t` (a hole bridge at iterate −1 of the window references) reads `aa u^-1 ab u aa^-1`; `aa = {R(0.0)#1, R(1.0)#27}` is the class of the image bridges `(40,41)`, `(44,45)` of hole bridge `(38,39)`, and its trivial walk names `{R(1.0)#1, R(2.0)#25}`, which was `new1`. On `p3_16` the first such classes are in the evidence table: their members are exactly the registered images of the image bridges. Fix (C9): `symbolic_dynamics` closes the alphabet under the word map over REGISTERED bridges only. A class a word names that some bridge spans copies its evidence entry, is walked in turn, lettered by the session (`letter_for`) and records `named_by`. A pair no bridge spans is a construction error (D3): the class whose itinerary names it is unresolved ("grow or blast"), one WARNING names the pairs and the members whose image predicts them, nothing is lettered, and `is_reliable` is False. Virtual-class support is deleted. A first version of the closure also lettered bare pairs (no members, active by the no-evidence convention, INFO only, reliable): on the deep runs 82 of its 90 closure classes had no bridge behind them (p3_6_blasts 20/20, p3_15 20/20, p3_16 19/22, nested 6 inner blasts 20/20, k=2.8 at 8 blasts 3/8), a transport chain predicted from the walk alone; the code review rejected it against D2/D3. Now: `p3_15`, `p3_6_blasts` and nested (6 inner blasts) add no class and leave one class unresolved (`aa`, the image-bridge class); `p3_16` closes 27 -> 30 (two closure classes verify against their registered images, the third names a pair no bridge spans); k=2.8 at 8 blasts 15 -> 20 with two unresolved. Whether those pairs must get their bridges is DEC-21. "Crossing abuts no hole" on `p3_6_blasts` is crossing 44, an excursion of the image chain of `(38,39)` across s(1.0); it blocks no walk (PK-4, slice 2). Laws still failing on the deep runs, also on base and not WP7's: `partition_covers_branch` (BCG-4 hole chains, WP8) and `dual_face_side_is_geometry` (representative points of long anchor faces, INT-10 / WP10).

**WP3 evidence (inversion x period 2):** the `inversion_period2` law case (k=4.1, 9/9 steps) passes every law. One step deeper (10/10 steps, 47 crossings) eleven laws fail on the WP3 branch, every one also failing on main (which fails sixteen): `bridge_endpoints_on_own_branch`, `bridge_rows_consistent`, `partition_covers_branch`, `backward_holes_keep_bridge_side`, `openings_missing_only_at_anchor_or_tail`, `alphabet_inside_minimal_trellis`, `iterated_cut_provenance`, `every_active_class_walks` (an unreachable walk), `no_virtual_class`, `landings_contained`, `refined_children_inherit_word`. Deep period-k evidence for WP5 / WP7 / WP10; the law case stays at 9 steps. **WP5:** at 10/10 steps it is not a WP5 gap. Seven laws fail on main `3e12d53` (`bridge_endpoints_on_own_branch`, `partition_covers_branch`, `backward_holes_keep_bridge_side`, `alphabet_inside_minimal_trellis`, `every_active_class_walks`, `no_virtual_class`, `landings_contained`). `partition_covers_branch` fails on an element of s(0.0) left opened at both ends by two holes of ONE origin a lap apart (iterates -3 and -7), the -3 one on a cross-branch container through the inversion point: the INT-2 whole-bridge convention meeting the law's one-hole-per-open-element reading, not the missing `+3` hole. The landing failures are outermost elements whose image runs past the next branch's partition (INT-8). Owners WP7 / WP8; not a law case.

**WP8b (D11):** a walk pair no computed bridge spans is a class all the same: `symbolic_dynamics` adds a member-less `BridgeClassEntry` (active by the no-evidence convention, `tangle` from its elements), letters it by the session, walks it from its elements' landings and flags it "no computed bridge yet" (`BridgeClassEntry.spanned` False, `describe()`); `_unspanned` is gone and nothing is unresolved by such a pair. Measured (`wp08b_impl/deep.py`, `lawsig.py`): `p3_15` and `p3_6_blasts` 23 -> 43 classes (20 flagged), `p3_16` 27 -> 49, nested with 6 / 7 inner blasts 27 -> 47 / 28 -> 49 (20 flagged each); all `is_reliable` True, no unresolved, ambiguous or straddling class (regression `test_open_p3_deep_runs`). On `p3_series` the closed alphabet is 43, 45, 47, 47 at blasts 6..9 against 23, 25, 27, 28 on base, where it grew one class per blast at a constant hole set (law `series_closed_alphabet_constant_while_holes_constant`). k=2.8 at 8 blasts: 15 -> 23 classes (3 flagged), nothing unresolved, but still NOT reliable: ambiguous `q`, `aa` and an evidence disagreement on `r` (WP7's u x u artifacts, WP10). Inversion x period 2 at 10 / 10 steps (D10): 14 seed / 17 closed classes, every class walks, reliable (`test_pip_orbit_cut`). "Crossing abuts no hole" on `p3_6_blasts` is PK-4 (untouched).

**Evidence:** `tests/cases.py:831-839`; `tests/regression/test_open_p3_deep_runs.py:1-18, 58-67`; `src/tanglepack/examples/henon_cases.py:30-39`; `src/tanglepack/topology/PartitionFamily.py:800`; `src/tanglepack/topology/DualWalk.py:757-760`; `stress/numerics/report/report.tex:1477-1527, 1734-1735, 2025-2033`; `stress/numerics/report/numbers_extra.tex` (NumDeepUnresolvedUnreachable=127); `stress/numerics/report/numbers.tex:468-472, 576-580, 678-680`; `stress/numerics/results/e11_blast_depth.jsonl` (p3 1e-5 `invariant_assert` blasts 6-19); `.refactor/knowledge/topology.md` (OPEN p3 ≥6 blasts); `.refactor/decisions.md` 2026-10-02

### PK-2: The provisional rule that punches forward iterates +1..+(k_value−1) (`_is_forward_beyond_fundamental`)
**Status:** closed (WP5) | **Severity:** medium | **Related:** PK-3, BCG-4, DEC-6 | **WP:** WP5

This is how the other k−1 branches get their window holes. Only the "≥k_value not punched" half is tested. The chalkboard alternative: "a reference pair is any consecutive pair inside a window F_j on any ray", which needs no table and no exemption.

A side effect on the baseline is the inert chain `u(1.0)→v(2.0)→x(0.0)→loop`. It is not Z3-equivariant (x does not map to u) because the image of the hole at iterate k−1 is not punched. That is by design.

**WP5:** decided by D5 and kept: `_is_forward_beyond_fundamental(pair)` is `iterate >= k_value` alone, and the positive half is now tested by the law `window_holes_on_every_branch` (a direct hole at every `+j`, `j < k_value`, whose image pair is registered and spanned, on `advance_key(branch, j)`; counts k10 1, k28_one 2, k28_two 1, p3 6, nested 7, nested_inner 7, inversion 1, inversion_one 3, inversion_two 3, inversion_period2 2). The window-walking alternative is rejected: see PK-3.

**Evidence:** `src/tanglepack/topology/StablePartition.py:173-190, 1712-1723`; `tests/conftest.py:72-77`; `CLAUDE.md` "Holes propagate BACKWARD only"; memory `holes-backward-only-2026-09`

### PK-3: When `extend_pseudoneighbor_trajectories` finds no table links, the other k−1 rays silently get no window holes
**Status:** closed (WP5 + DEC-19: asymmetry accepted, 2026-10-08) | **Severity:** high | **Related:** PK-2, BCG-3, PK-1, NUM-18, DEC-19 | **WP:** WP5

Probe p3 at 6 blasts: `pairs=2` and zero forward pairs, so rays 0.0 and 2.0 got no direct window holes (at 5 blasts they had +1/+2 holes). The partition becomes asymmetric across the orbit, and the provisional rule is never even applied.

**WP3 evidence:** the `inversion` law case (0 blasts) is a second instance. Its one reference pair `(10, 22)` on `s(0,0)` has no registered `+1` iterate, so its `+1` appearance on `s(0,1)` is never punched; stable branch 0.1 has no hole, the anchor bridge's image `(0, 7)` on `u(0,1)` is a loop in `L_(0.1;1)` and the anchor class is inert. The strict xfail `anchor_bridge_class_leads_its_tangle` x `inversion` is re-tagged here (WP5); at 1 blast the `+1` pair exists and the law passes.

**WP5 (slice 1):** root cause found. With WP2's exact table every missing `+j` link is an image crossing that is NOT COMPUTED, never a missed match: its `image_cdist` landing lies past the last crossing of the unstable branch it would land on. A period-k chain is grown one map step at a time, so its branches are staggered (each reaches β further than the one before it in the cycle), and a reference on the outermost lap of the longest branch has no forward image. Measured landings: inversion `+1` at u 55.5 on u(0.0), whose last crossing is at 13.2; p3 at 6 blasts `+1` at 4097 against 1662; inversion_period2 `+1..+3` at 1.13 / 1.55 / 2.11 against 0.38 / 0.52 / 0.71. Gaps by builder: inversion (s(0.1)), inversion_one_blast (s(0.1)), inversion_period2 (s(0.1), s(1.1), s(0.0)), p3_15 / p3_6_blasts / p3_coarse_blasts (s(2.0), s(0.0)); k10, k28 x3, p3, nested, nested_inner_blasts and p3_16 have none; inversion_two_blasts has a registered `+1` pair no Bridge spans (INT-9, WP7). `extend_pseudoneighbor_trajectories` now logs ONE WARNING per stable branch left without a reference's appearance, naming the references and each missing member's predicted landing against the branch's last crossing; nothing is invented (regression `test_period_k_window_gaps`), and `punch_holes`' no-spanning-bridge WARNING names each skipped pair, its iterate, origin and branch. Walking the gap branches' windows does find pairs there (p3_6_blasts s(2.0) `(15,23)`, s(0.0) `(18,20)`; inversion s(0.1) `(17,3)`), but they are ORPHANS: each one's open stable interval holds the predicted images of the tip crossings, so by the definition (X = every iterate) they are no pairs and pass only because the stable-side check reads registered crossings alone (NUM-18). Symmetric per-branch hole counts are NOT a WP5 property: p3 is 4/4/4 at 0-5 blasts and 16/14/12 at 6, and p3_16, which has no gap, is 14/16/16 (the backward chain of 20 steps per origin is not a multiple of k = 3: BCG-4, WP8). Filling the gap needs the author (DEC-19).

**Evidence:** chalkboard review §4 probe (p3 6 blasts); `src/tanglepack/topology/StablePartition.py:1712`

### PK-4: Cross-branch image-bridge cuts are skipped with a WARNING and unvalidated
**Status:** closed (WP10b) | **Severity:** medium | **Related:** PK-1, INV-6 | **WP:** WP7

These cases only warn and do nothing:
- an empty stretch whose far end is on another stable branch ("no flank")
- a chain that ends on another branch (marks nothing)
- an odd number of interior crossings
- "no empty stretch"

`BridgeClass` cannot place image endpoints that sit on different unstable branches. The only checks are ownership invariants. This matters for period-k and inversion, where image chains switch branches.

**WP7 evidence (slice 2, needs the author's rule):** on `p3_6_blasts` "crossing abuts no hole" is crossing 44. The image chain of hole bridge `(38,39)` (`40..45`) leaves its base branch s(0.0) for an excursion across s(1.0) (41-44, crossings with `S \ f(S)`, no `−1` link), and `_empty_stretches` pairs chords on the base branch only. Sorted along s(1.0) the excursion reads 43 (0.973) < 42 (1.025) < 41 (1.158) < 44 (1.285). By Jordan alternation `f(lobe) ∩ s(1.0) = [43,42] ∪ [41,44]`, where `[43,42]` is the stretch of the window hole (42,43), so the missing empty stretch is the chord `(41,44)` and 44 belongs to the element outward of it. That pairs a foreign branch by STABLE order, while the author's base-branch rule pairs by unstable order (`(x1,x2),(x3,x4)`). It changes refinement labels only (`aa_1` / `aa_2` stay distinct, and the trellis word of `t` verifies either way) and blocks no walk.

**WP10b (C19):** implemented. The rule is the clean one: a lobe H (its bridge and the stable segment between the bridge's ends) holds no point of any stable branch of the trellis, so f(H) holds no anchor, and walking any stable branch other than the base from its anchor enters and leaves f(H) at successive crossings with the image chain. `_empty_stretches` therefore pairs the chain's crossings per stable branch: chain order on the base branch (the author's rule), STABLE order on any other. Measured geometrically (point-in-polygon on the lobe image polygon, `test_foreign_chords.py`): on `p3_6_blasts`, `p3_15`, `p3_16` and the law case `heteroclinic` every foreign branch reads out, IN, out, IN, out from its anchor, the IN stretches being the window hole `(43,42)` and the new chord `(41,44)` (ids of `p3_6_blasts`). The other bullets: "a chain that ends on another branch" was the anchor (an anchor lies on every branch, D4, but its C7 branch was read as the base: three WARNINGs on `inversion_period2`); the base is now read at a non-anchor end and asserted to lie on one branch. "No flank" can only be reached by an anchor (every empty stretch lies on one branch) and is asserted for a non-anchor. An odd count on one branch still warns; it fires on no law case, deep run or series. "No empty stretch" now fires only at the anchor lobe's fold (`inversion_two_blasts`, crossing 2 of `(1,2)`), which by the author's 2026-09-23 rule marks nothing: by design. Effect: only the iterated partition of `heteroclinic`, `p3_15`, `p3_16` and `p3_6_blasts` changes (+3 iterated elements each, so iterated names, iterated itineraries and refined-pair labels); class tables, holes, words, refined rules, matrices and reliability are identical. PartitionFamily WARNINGs 2 -> 0 on those four, 3 -> 0 on `inversion_period2`. The `BridgeClass` cross-branch skips (`_image_chain` off-branch image end, `_collect_image_evidence`'s uncovered branch) fire on nothing measured and are left for the D8 API pass.

**Evidence:** `src/tanglepack/topology/PartitionFamily.py:714-717` (Dev Note), `:727, 774, 799, 813, 941, 955, 976`; `src/tanglepack/topology/BridgeClass.py:735, 798`; `CLAUDE.md` "Partition families"; memory `minimal-dual-graph-2026-09` ("cross-branch cut is the open point")

### PK-5: The fundamental segment is unified only on the pip's branch; the other k−1 branches are walls everywhere
**Status:** closed (WP8b, D10) | **Severity:** low | **Related:** PK-1, DEC-6, DEC-20 | **WP:** WP7, WP8b

A walk that needs to cross another branch's window becomes "unreachable". This is one plausible contributor to the p3 unreachable class.

**WP7 evidence (DEC-20):** the only unreachable walk left anywhere is inversion x period 2 at 10/10 steps (class `n = {L(0.1)#1, R(0.0)#5}`: its start face has no unified exit). There the stable branches run PAST the strong pip's orbit (s(1.0) ends at 4.41, q0 = 1.75), so the pip window `(f^4 q0, q0] = (0.504, 1.749]` is not where `f(W^s)` stops. The image of a bridge crosses the trellis's stable manifold only at the images of its ends and in `S \ f(S)`. Per branch that is `(c_s(f(end of the previous branch)), end]`: `(0.577, 4.413]` on s(1.0), empty elsewhere. Unifying `S \ f(S)` instead makes `n` walk. On every law case and the three deep runs `S \ f(S)` equals the pip window exactly, because every builder trims the stable manifold at the pip orbit. Either fix (unify `S \ f(S)`, or cut S at the pip orbit in the partition and dual graph as the chalkboard's stage 3 does) reverses the 2026-09-14 "the pip moves the unified set without trimming" design, so WP7 only states the precondition: law `dual_unified_is_where_f_of_s_stops`, and one DualGraph WARNING per branch where the two differ (regression `test_unified_segment_precondition`). The choice is DEC-20.

**WP8b (D10, with D12):** root cause: only `Pseudoneighbor.compute_pseudoneighbors` cut N at the pip orbit (`_strong_pip_cuts`, `_within_cut`); holes, the partition, the arrangement and the dual graph read the UNCUT branches, so on inversion x period 2 at 10 steps (pip 30 on s(1.0); s(1.0) ends at 4.413 against its cut 1.749, s(0.1) 3.233 / 1.280, s(1.1) 1.075 / 0.938, s(0.0) 0.788 / 0.688; 15 of 47 crossings past a cut, 9 of them mid-arc) the partition kept elements beyond the pip, the backward chain landed holes on bridges ending beyond it, and the dual graph had walls where `S \ f(S)` needs a doorway (class `n` unreachable). Fix: choosing the pip cuts the trellis to `T(q0)` in ONE place, `Trellis.cut_at_pip_orbit` (the `strong_pip` setter; None restores), so every stage reads the same N; the pseudoneighbor stage's private cut is deleted; the arc across a dropped crossing joins its two bridges (`Bridge.joined`, never a workbench bridge: `Topology._zone_key` and `cartoon_enclosures` resolve bridges through the trellis); `DualGraph._check_where_f_of_s_stops` is gone, its WARNING vacuous (`S \ f(S)` is the pip window by construction; the other k-1 branches walls as a consequence, not DEC-6's per-branch unification). The cut equals a numerical trim at the pip: on inversion x period 2 at 10 steps and on every k10 candidate the holes (by cdists, iterate, side), class count, words in element names and reliability agree with `resonance_zone(pip)` on a fresh build, or both raise the same assertion (k10's second candidate: the anchor bridge wraps the end of the cut branch, I2) (`test_pip_orbit_cut`). Every other builder trims at the pip orbit already, so no law case, series or deep run loses a crossing to D10. New law `trellis_is_cut_at_the_pip_orbit`.

**Evidence:** `src/tanglepack/topology/DualGraph.py:583-611`; `src/tanglepack/topology/DualWalk.py:757-760`

### PK-6: A strong pip with fewer cut points than k_value only warns
**Status:** closed (WP5) | **Severity:** low | **Related:** PK-3 | **WP:** WP5

If the pip's iterates are not registered on all branches, the windows on the remaining branches stay uncut.

**WP5:** the cut is a stable-cdist threshold, so `_strong_pip_cuts` places each of the k cuts by `Trellis.image_cdist(pip, j, "stable")` (table first, then scaling, lookup only) and flags a placed cut in one WARNING; no branch is left uncut. The unused cut ids and `strong_pip_cut_points` are gone from the module. Cuts are complete on every builder, so no case output changes.

**Evidence:** `src/tanglepack/topology/Pseudoneighbor.py:444-450`

### PK-7: Crossings with `manifold_a_key is None` are tested at every cycle residue (can reject a true pair)
**Status:** closed (WP2) | **Severity:** low | **Related:** SIMP-12, DOC-3 | **WP:** WP2

Verified 2026-10-07: blast children carry their key (`BridgeIterator.py:334-335`, `:388`), and `Tangle._key_of` asserts a key on every unstable curve. So the `None` branches look dead, not conservative. The Dev Notes still say "legitimately None".

**WP2:** confirmed dead on every builder; `Tangle._key_of` now asserts a key on every curve and the None-key branches (Pseudoneighbor residues, StrongPip, BridgeClass, IterateInference) and their Dev Notes are deleted.

**Evidence:** `src/tanglepack/topology/Pseudoneighbor.py:84-93, 561-567, 591`; `src/tanglepack/topology/StrongPip.py:84-87, 232`; `src/tanglepack/numerics/IterateInference.py:89-94, 254-361`; `src/tanglepack/numerics/Tangle.py:220-228`

### PK-8: Noisy WARNING "inert class u has the non-empty word v" on period-3 runs
**Status:** closed (WP5) | **Severity:** low | **Related:** PK-2 | **WP:** WP5

For period k it is expected that an inert class maps to the next branch's inert class, so this warning is noise.

**WP5:** `_edges` warns only when an inert class's word names a non-inert (active or virtual) symbol, which contradicts the table's inertness; a word of inert symbols is the definition holding and is silent. p3, nested and nested_inner_blasts lose their transition-graph WARNINGs; no law case has an inert class naming an active one.

**Evidence:** `src/tanglepack/topology/SymbolicDynamics.py:898`; `.refactor/decisions.md` 2026-10-02

### PK-9: The hole-chain truncation at `max_steps=50` is logged only at DEBUG
**Status:** closed (WP5) | **Severity:** low | **Related:** BCG-4 | **WP:** WP5

A deeper reference chain would be cut silently.

**WP5:** the guard is deleted (`max_steps` had no caller). Every backward step either stops (no container, or a visited one) or adds a new bridge id to `visited`, a subset of the trellis's finite bridge set, so the chain ends within `len(trellis.bridges) + 1` steps and no truncation can happen. The longest chain seen is 21 steps (p3_6_blasts, p3_16, p3_coarse_blasts).

**Evidence:** `src/tanglepack/topology/StablePartition.py:437-438`, `~583`

---

## Nested / heteroclinic tangles

### NEST-1: No fixture has a heteroclinic crossing; connecting classes, the merged outer node and cross-tangle walks have never run, and a hand-built case is unreliable
**Status:** closed (WP9c) | **Severity:** high | **Related:** NEST-2, INT-10, TEST-1, DEC-7, DEC-23, DEC-24 | **WP:** WP9

**WP9:** the heteroclinic path runs end to end. `henon_cases.build_heteroclinic()` (p3 16u / 10s, p1 11u, two outer blasts, 1e-7, ~0.7 s) is the law case `heteroclinic`: four `W^u(p3) x W^s(p1)` crossings, two inside both pip-orbit cuts (a thin lobe across p1's window); the other two are past p3's unstable cut (points of X). With D7 (NEST-2) the lobe pair is a reference of p1, its chain holes p3's partition, and connecting classes (`entry.tangle is None`), the merged outer face (the sparse arrangement is ONE component), cross-tangle walks and `HETEROCLINIC_COLORS` all run; every law passes but two strict xfails: `every_active_class_walks` (DEC-23: two shortest walks of equal length pass either side of the lobe's hole on p1's window, so three classes are ambiguous and the dynamics unreliable; slice 2, needs the author) and `partition_covers_branch` (inherited from p3 at 16 unstable steps alone, WP10). Open: DEC-23; escape-lobe images (none of the four crossings has an iterate link: the preimage lies past p1's stable trim, the image past p3's unstable cut; WP10); the `W^u(p1) x W^s(p3)` recipe (p3 13u / 15s) hits a NUM-4 near-duplicate crossing (WP10). **WP10:** the escape-lobe images are NUM-21 (accepted, D9's category); the "near-duplicate" was two crossings on one segment (NUM-19), and the recipe builds now; its heteroclinic iterate links fail the beta and area laws (NUM-20, deferred).

**WP9b:** with the non-crossing rule (D13) the classes `a` and `c` keep one walk each (`a` the same word, `c` `a d^-1 a^-1`, its registered word), and `bh`, whose unique shortest walk was trivial and contradicted by its registered image, reads `bi bj bi^-1` (its registered word, through a two-crossing detour around the lobe; C13). Every registered image with a verdict agrees with its walk on the law case. `bn` = {A:R_1, B:L_(2.0;1)} stays ambiguous: it is member-less (D11, named by `bl`), and the marker law needs a closed curve with a known image through it -- for a connecting class, an arc closed through a heteroclinic crossing h of N with `f(h)` or `f^-1(h)` in N. The trellis has none: its two heteroclinic crossings in N (and the two past p3's unstable cut) carry no iterate link at all; the image of each lies past p3's unstable cut and the preimage past p1's stable trim (NUM-21). So the side of the B peninsula `bn`'s chord passes is not determined by T(q0); the strict xfail stays, tagged DEC-24. p3 at 17 unstable steps holds eight linked heteroclinic crossings in N but takes 27 s per build and still shows one connecting tie.

**WP9c (D14):** traced, `bn` is no construction error. `named_by` gives `bn <- bl <- bi`, and `bi = {B:L_(0.0;1), A:R_(0.0;5)}` (A = p1, B = p3) has two COMPUTED members, `(91, 92)` and `(93, 94)`: arcs of `W^u(B 2.0)` from B's stable branch to the heteroclinic crossings 92 and 93 on A's window. `bl` is `f(bi)` and `bn` is `f^2(bi)`; they are member-less only because `f(92)` and `f(93)` lie past B's unstable cut (D12, NUM-21), D11's situation exactly. Every landing is table-exact and contained; face 0 is one open face whose ring holds B's left-side edges and A's right-side edges, and the B arms run through it, so a piece of `f(b)` crossing it from a B element to an A element is legitimate (the pair (A:R_1, B:L_(2.0;1)) is the start/goal of `bl`'s trivial walk); A's partition is sound. The defect was that D13 was implemented only as the marker law on a member LOBE, which needs a member with both ends on one stable branch, so a class joining two stable branches had no rule (a connecting class, with or without members: at 17u the connecting class `bm` has two registered members and is still ambiguous on base). D13's own words decide it: "f(bridge) cannot cross the images of the unstable arcs". Images of disjoint arcs are disjoint, so a walk may not cross, inside a face, the decided walk of an earlier class whose arcs lie on other unstable branches (C17: a class's arcs lie on its members' branches, or one map step past its namer's). On `heteroclinic` this rejects exactly `bn`'s second walk (its chord in face 0 interleaves with those of `bi`'s and `bl`'s walks; the kept one reads the word base printed first); at 17u it resolves `bm` to its registered word. Left elsewhere: NUM-20 and NUM-21 (no longer needed for any word).

`build_nested()` has 0 A×B crossings (only A×A 30 and B×B 10) and 0 connecting classes. Never run on real data:
- the merged outer node, which is "over-generous once two tangles are linked heteroclinically"
- `HETEROCLINIC_COLORS`
- connecting-class ordering
- cross-tangle walks

A real case (p3 16u/12s, p1 11u, 1e-7, about 1 s, 4 A(u)×B(s) crossings) runs end to end but gives:
- 42 classes, 2 virtual
- a cross-side pair
- a standalone "connecting" loop
- 3 WARNINGs "Unstable branch … missing from the branch cycle"
- `ah` disagreeing with its registered image
- `ai = {A:p3@0.0/L#0, B:p1@0.0/R#2}` with no member image

At 1e-5 the same case raises an AssertionError (NUM-1). The p1 arm hits the budget at 14 steps, so the fixture recipes sit where no heteroclinic crossings exist.

**Evidence:** probe `nested.py`; period-k/nested investigation §2; `docs/dual_graph_symbolic_dynamics_plan.md:331`; `docs/regions_bridge_identity_and_cleanup_plan.md:172`; `src/tanglepack/topology/Arrangement.py:60-75`

### NEST-2: Heteroclinic pseudoneighbor pairs are rejected by construction, with a WARNING per candidate
**Status:** closed (WP9) | **Severity:** high | **Related:** NEST-1, DEC-7 | **WP:** WP9

**WP9:** two root causes. (1) The pair test built the cycle, `k` and β from the WINDOW's fixed point and dropped an unstable key outside it with the WARNING; it now reads them from the ARC's fixed point (`_unstable_interval_is_clear`), the scaled 2D endpoint collision applies only to a candidate on the pair's own stable manifold, and the WARNING is gone. (2) Stages 4-6 ran on each fixed point's own trellis, which holds neither the other tangle's arcs nor its bridges; they run on the selection's trellis (`loom.Topology`), and a crossing is in a trellis only when both its branches are selected (`Trellis.from_workbench`; a tangle's own trellis joins its arcs across heteroclinic crossings). Regression `test_heteroclinic_pairs_hole_the_stable_owners_partition`.

Two consecutive heteroclinic crossings on B's stable branch carry A's unstable key, which is not in B's branch cycle, so they are rejected. No hole can ever sit in a heteroclinic lobe.

**Evidence:** `src/tanglepack/topology/Pseudoneighbor.py:535-540`

### NEST-3: Nested is reliable only for outer ≥2 and inner ≤4 blasts; nested_outer loses classes and holes as it is blasted
**Status:** partly closed (WP7) | **Severity:** medium | **Related:** PK-1, BCG-5 | **WP:** WP7

**WP7:** inner ≥6 blasts is exactly PK-1 and shares its state: no virtual class, every class walks, but one class names a pair no bridge spans and is unresolved with a WARNING, so the run is not reliable (6 inner blasts: 27 classes, no closure; 7: 28 -> 29); DEC-21. The unresolved outer class at 0 outer blasts is a singleton landing whose members have no registered image: D3's f(bridge) path needs the registered chain, so "grow or blast" stays, by design. The nested_outer loss under blasting is BCG-5 (WP8).

- **0 outer blasts:** the outer class `e` is a singleton with no registered chain.
- **Inner ≥6 blasts, or p3 at 15 steps:** the run inherits PK-1 exactly (27 classes, 42 p3 holes); E11 nested_inner goes 13→55.
- **nested_outer:** classes 12→13→14→11, holes 16→12→6, active 6→3 and `minimal_image` 3→0 at blast 10. Unexplained.

**Evidence:** `src/tanglepack/examples/henon_cases.py:40-47`; `.refactor/knowledge/loom.md` ("Reliable iff outer ≥2 and inner ≤4"); `CLAUDE.md` ("Nested needs two outer blasts"); `stress/numerics/results/e11_blast_depth.jsonl` (nested_outer); `stress/numerics/report/report.tex:1286-1288`

### NEST-4: Blast containment is by zone polygon only; the library still lets an outer blast iterate the inner tangle
**Status:** closed (WP9, C4) | **Severity:** medium | **Related:** NUM-9, DEC-11 | **WP:** WP9

**WP9:** `blast_zone`'s `in_zone` is `session.classify_bridge(bridge) is zone` (the innermost containing zone). 4 outer + 2 inner unfiltered blasts: p3 holes 4 (outer first) vs 42 (inner first) on base, now equal both ways and equal to the own-fixed-point recipe. An unfiltered outer blast still iterates the inner tangle's bridges in the annulus between the zones (they are the outer zone's, C4 allows it), and one of them can still hit NUM-9 (WP10); the recipes keep `fixed_point=[zone.fixed_point]`. Regressions `test_an_unfiltered_outer_blast_iterates_only_its_innermost_bridges`, `test_an_unfiltered_nested_blast_does_not_depend_on_order`.

`build_nested` passes `fixed_point=[zone.fixed_point]`, but `blast_zone` uses `zone.contains_point` plus whatever filter the caller passes. A caller can reproduce three failures:
- order dependence
- the `merge_manifolds` crash
- silent p3 corruption at ≥6 outer blasts

**Evidence:** `src/tanglepack/loom/Blast.py:278-292`; `.refactor/knowledge/loom.md` (Smells); memory `hole-side-own-blast-2026-10-02`

### NEST-5: `rebuild_bridges(fp)` clears every fixed point's bridges
**Status:** closed (WP9) | **Severity:** medium | **Related:** NEST-1 | **WP:** WP9

**WP9:** `clear_bridges(fixed_point=None)` drops only that fixed point's bridges and partial pieces and rebuilds `_bridges_at` from what is left; `rebuild_bridges(fp)` clears and recuts only `fp` while the ids are current, and every fixed point after a renumbering recompute (the survivors' ids would name other crossings). Partial pieces do not accumulate across rebuilds. Unit `test_rebuilding_one_fixed_points_bridges_keeps_the_others`.

The docstring says "The CLEAR is global even when fixed_point is given". `clear_bridges()` is at :1849, so a per-fixed-point rebuild in a nested session drops the other tangle's bridges. Whether partial bridges accumulate per re-blast was not re-checked.

**Evidence:** `src/tanglepack/numerics/TangleWorkbench.py:1796-1851` (:1831-1833, :1849); `docs/regions_bridge_identity_and_cleanup_plan.md:333`

### NEST-6: `entry.tangle` follows the combined trellis order `[A(p3), B(p1)]`; `TangleBuild.fixed_points` is "outermost first" `[B, A]`
**Status:** closed (WP9) | **Severity:** low | **Related:** — | **WP:** WP9

**WP9:** the recipes construct the outermost fixed point first (`henon_cases._orbit` / `_grow`; `tests/cases.py::_nested_through` the same), so p1 is `A` and tangle 0 and `entry.tangle` indexes `TangleBuild.fixed_points`; the orbits are still grown p3 first, each at its own cutoff. Holes, references, classes and reliability of nested / nested_inner_blasts are unchanged; labels and registry ids move (p1's classes take the first letters). Law `tangles_in_case_order`.

The two orders are opposite, which invites mix-ups. Align them or document the difference.

**Evidence:** period-k/nested investigation §3.5; `src/tanglepack/examples/henon_cases.py`

---

## Bridge-class growth

### BCG-1: Every trellis bridge is classed (blast children included), and members with no image count as active, giving +1 active class per blast
**Status:** closed (WP6) | **Severity:** high | **Related:** BCG-2, BCG-3, BCG-6, DEC-2, TEST-2 | **WP:** WP6

`bridge_classes` defaults to all non-partial bridges, and the session passes no `bridge_ids`. The whole table reaches `symbolic_dynamics` and the transition matrix.

From blast 8, p3 has these counts constant: holes 44, homotopy elements 85, iterated 90, `minimal_kept` 28, `minimal_image` 4. Yet classes rise 25→51, +2 per blast, which is +1 active (this issue) and +1 inert (BCG-2).
- The newest class is the image of the previous blast's newest: (89,90)←(79,80)←(71,72)←(63,64).
- These march one element pair inward: `1.0/R#21↔R#25 → 2.0/R#19↔R#23 → 0.0/R#16↔R#20`.
- The canonical rules form a shift chain `A3→A4→…→A37→new1`.

Classes with a minimal-trellis member stay at 25 (p3 blasts 8-16, while the total goes 27→43); k28 stays at 14-20 between switches. The fix is to build letters from `minimal.kept_bridge_ids` and keep blast children as evidence only. E7 shows the same thing with depth: p3 classes 9 at 11-14 steps, then 23/30/36 at 15/16/17.

**WP6:** `bridge_classes(trellis, partitions)` is the EVIDENCE table (every bridge); the session's minimal trellis is built from it (activity is a fixed point over ALL bridges' images, so the order is evidence -> minimal trellis -> alphabet) and keeps it as `MinimalTrellis.classes`; `MinimalTrellis.alphabet` (`BridgeClassTable.restricted_to(kept_bridge_ids)`, same entries, members whole) is what `session.bridge_classes()` letters and what the words and the matrix use (D2). A second leak is closed with it: `registered_evidence` translated member images with the SHARED letters, so a blast child's image naming a class outside the table became a virtual class and a matrix node; it now translates on a copy. Series (`cases.SERIES`, total classes): `k28_series` 3,3,4,4,4,4,6,7,23 -> 3,3,4,4,4,4,4,4,15 (matrix 4 = the active refined symbols at blasts 2-7, 0 virtual); `p3_series` 9 x6,23,25,27,28 -> 9 x6,23,24,25,25; `inversion_series` 4,10,10 -> 1,9,9. E11 small on the branch: p3 alphabet 25 at every blast 8-16 (evidence 27 -> 35), matrix 23 and 1 virtual over 9-16; k28 alphabet 4 (blasts 2-7), 15 (8-10), 19 (11-12). Law cases unchanged except inversion (4 -> 1 class: the inert hole-bridge class alone, no active class until DEC-19), inversion_one_blast and inversion_two_blasts (10 -> 9: an active class with no kept member leaves). Members with no image still count as active (CLAUDE.md's evidence rule); with a minimal alphabet that no longer grows anything.

**Evidence:** `src/tanglepack/topology/BridgeClass.py:616-622, 789-792, 833`; `src/tanglepack/loom/TangleSession.py:391`; `src/tanglepack/topology/SymbolicDynamics.py:1251`; `stress/numerics/results/e11_blast_depth.jsonl` (p3/nested_inner 1e-7); `stress/numerics/report/report.tex:1238-1246, 1284-1291`; `stress/numerics/report/numbers.tex:646-666, 854-866`; `stress/numerics/report/tables/e11_summary.tex`; `e07_breaking.jsonl`; probes `probe_p3.py`, `minimal_classes.py`

### BCG-2: Loops are folded by following preimage links, so a loop born from a fold stands alone as an inert BridgeClass(x,x), giving +1 inert class per blast
**Status:** closed (WP6) | **Severity:** high | **Related:** BCG-3, INT-1, DEC-2 | **WP:** WP6

`_non_loop_ancestor` walks `iterate(x,−1)`. Crossings created by a fold have no −1 link, so the chain breaks. Example: loop (61,62) with preimages `[None, 55]` goes to its own class, and its images (69,70), (77,78) each stand alone in a different element (R#21, R#19, R#16…).

The parent's `_image_chain` does contain (61,62), and the author's rule ("the image of a bridge that did connect two elements") is the forward relation, so the intent did not make it into the code. Fix: fold via the bridge whose image chain contains the loop, falling back to `BridgeIterator.preimage_bridges`. Put true orphans in one bucket per orbit.

**WP6:** WP2's exact table had already removed every standalone loop but inversion's four at 0 blasts (`(9,6)`, `(17,25)`, `(7,12)`, `(8,24)`, under the anchor bridge whose `-1` link breaks the backward walk). A loop now folds FORWARD: `bridge_classes` computes every bridge's image chain once, indexes each image pair by the one bridge whose chain holds it (asserted: the images of disjoint arcs are disjoint; 0 double parents measured) and climbs those parents through loops to the first bridge that connects two elements (`_non_loop_parent`). It agrees with the old walk wherever the old walk found an ancestor and leaves 0 orphans on every law case, the deep p3 runs and every series step; inversion's four fold into one class (evidence 4 -> 2 classes). The `preimage_bridges` fallback is redundant with a partial-bijection table and is not used. A true orphan is logged once and listed on `BridgeClassTable.orphans`, never classed or lettered (D2), so `_non_loop_ancestor`, the standalone `BridgeClass(x, x)` entry and its inert seed in `_resolve_inertness` are gone. New law `loops_fold_forward`.

**Evidence:** `src/tanglepack/topology/BridgeClass.py:636-647, 706, 842-875` (:863-866); E11 p3 `n_inert` 1,1,1,2,3,…,13 over blasts 6..20; probe `probe_p3.py`; `CLAUDE.md` "a loop with no registered ancestor stands alone"

### BCG-3: Blasts never re-run registry-wide iterate inference, so interior image crossings stay unlinked
**Status:** closed (WP2) | **Severity:** high | **Related:** BCG-2, INT-6, PK-3 | **WP:** WP2

- `iterate_bridge` only calls `register_endpoint_iterates`.
- `infer_iterates` runs only inside `compute_intersections`.
- The recipe calls it once before the blasts, and `Blast.py` never calls it.

**WP2:** the root cause was twofold: the endpoint path's ABSOLUTE map gate `1e-3 * max(1, |f|)` rejected the six true links (`|f(x) - y|` = 1.2e-3 .. 3.5e-3 on ~1e-2 segments), and no registry-wide pass ran after a blast. The exact bracket rule with a segment-scale gate now runs in every `iterate_bridge` and at the end of every blast; after `build_period3(blasts=10)` a further `infer_iterates()` records 0 (regression `test_iterates_after_blasts`), and three more near-anchor links the 5% window missed are found.

Example: p3 crossing 61's true preimage 54 matches to 0.3-0.5% but `iterate(54,1)=None`. Running `session.workbench.infer_iterates()` after `build_period3(blasts=10)` adds 6 links: 31→29 classes and 2→0 standalone loops (another run: 29/27 active → 28/24 active).

**Evidence:** `src/tanglepack/numerics/BridgeIterator.py:393`; `src/tanglepack/numerics/TangleWorkbench.py:1380`; `src/tanglepack/numerics/IterateInference.py:119, 377`; `src/tanglepack/loom/Blast.py:195`; `src/tanglepack/examples/henon_cases.py:146`

### BCG-4: The reference pair switches when deep crossings fold onto the fundamental window, and the new backward hole chain has length ~ log(c_u)/log β
**Status:** closed (WP8: verified legitimate under D1) | **Severity:** high | **Related:** PK-1, PK-2, PK-9, INT-2, BCG-5, DEC-1 | **WP:** WP8

**p3, blast 6:** blast-child crossings 41-43 (c_u ≈ 2.6e3-3.1e3, ray 2.0) land inside the window (0.397, 1.481] of ray 1.0.
- The primary pairs `(18,29),(17,18)` stop being consecutive.
- The references become `(41,42),(42,43)`, each propagated back 20-21 bridges with opposite sides on the same bridges.
- Holes go 12→42 (origins 23 + 21), elements 24→81, with about 25 right-row elements per ray and a staircase `d→e→…→t` of 17 one-letter classes.
- `(42,41)` is a sub-arc of f¹⁷ of the old reference bridge `(17,18)`. So the switch also undoes the backward-only rule: old forward iterates (≥k) now get holes as negative iterates.

**k28, blast 8:** new references `(184/185,185/186)` (unstable cdist ≈2.36e8, resolution ratio ≈320, so not polygon noise). Their preimages are not registered because they are trimmed away, so the membership test runs on scaled cdists only. Holes 4→24 (3 direct + 21 carried), classes 4→20. Another switch at blast 11 `(892,893),(899,900)`.

**k32, blast 8:** classes 14→37, holes 16→56.

Classes at fixed depth also depend on the cutoff (NUM-5).

**WP8:** every switch follows the pseudoneighbor definition on the recomputed trellis, and every lost reference has a witness the iterate table names without scaling: a crossing strictly between the pair on one of its branches, or strictly between its registered images. k28 blast 2: `(3,4)` lost to 8, 9 between it on the stable branch, `(2,3)` punctured in the frame of its +2 images; new reference `(8,9)` (base numbering). k28 blast 8 adds `(184,185)`, `(185,186)` (unstable cdist 2.4e8, ratio 320); blast 11 loses `(8,9)` to 1203 and `(184,185)` to 1026, 893, 892, 1027 between, `(185,186)` in its image frame; new `(892,893)`, `(899,900)`. p3 blast 6: `(18,29)` lost to 41-44 between, `(17,18)` punctured by 41 in the frame of its +17 images; new `(41,42)`, `(42,43)` (41-44 at ratios 3.0-8.2, near the 3.0 threshold: NUM-3/4). The burst of holes after a switch is the backward-only chain of a deep reference (stage 5), by design. Law `series_lost_references_are_witnessed` (`helpers.laws.reference_witness`) checks every loss along `cases.SERIES`.

**Evidence:** `src/tanglepack/topology/StablePartition.py:433, 547-600, 1712-1723, 1771`; `src/tanglepack/topology/Pseudoneighbor.py:131, 505-620` (all-iterates test from :530); `stress/numerics/report/report.tex:1222-1230, 1277-1295`; `stress/numerics/report/numbers.tex:438-456, 540-564`; E11 k28 `n_holes [4x8,24,24,20,26,…]`; probes `p3probe2.py`, `probe_k28.py`, `p3_holes.py`, `k28_classes.py`; memory `hole-side-own-blast-2026-10-02:15`

### BCG-5: Holes and classes decrease under further blasting (the partition is non-monotone in depth)
**Status:** closed (WP8b, D12: the decreases WP8 left; WP8 closed the scaling-noise punctures); p3 at 1e-5 -> NUM-5 (WP12); p3 blast 23 not reproduced | **Severity:** medium | **Related:** BCG-4, NEST-3, DEC-22 | **WP:** WP8, WP8b

- nested_outer: classes 13→11, holes 16→12, `minimal_image` 3→0 at blast 10.
- k32: holes 58,50,42,34,13 over blasts 9-13 (classes 50→34, iterated elements 116→42), then 0, then 36.
- p3: holes 44→23 at blast 23.

This is consistent with the confirmed full-orbit disqualification, but the report says the burst and drop positions are "not checked".

**WP8: two kinds of decrease.** (1) A BUG, fixed: a false puncture from cdist noise. `_unstable_interval_is_clear` compared an UNLINKED crossing scaled back `beta**steps` with the pair's REGISTERED cdists, two frames. nested_outer blast 10: crossing 592 (unlinked, u 2907393829.12, ratio 5357) sits 6178 BELOW `M^8(37) = 593` on its own polyline, outside the image of `W^u(37,38)`, but scaled back eight steps it lands 1.3e-5 inside `(37,38)`, where 593's own scaled value misses 37 by 0.0103 (2e-6 relative); the outer tangle lost its only reference (holes 4 -> 0, classes 13 -> 9). k32 blast 12 the same (`(79,80)` punctured by 2094 at m=-6; holes 42 -> 34, classes 29 -> 24). An unlinked landing is now read against the pair's registered images on its own polyline (lookup first on the pair's side too), and scaling is the fallback only when neither image is registered. After the fix: nested at 10 outer blasts keeps `(37,38)` (holes `[4, 12]`, 13 classes); k32 at 12 blasts 42 holes, 29 classes; every law case, deep run and series unchanged (id-free signatures identical). (2) NOT root-caused, noise NOT ruled out (WP8 review): the remaining decreases are witnessed by a crossing within a hair of a window crossing on the stable branch, and the witness is unlinked in both directions, so no other frame can confirm its order. The resolution ratio (200-9000 here) measures the crossing ANGLE against polygon noise, not the cdist placement: the placement error bar along W^s is `min(l_u, l_s) / ratio` (= `max(h) / |sin|`), and the closest witnesses sit inside it. Measured on this branch (numbering of this branch; probes `wp08fix/k32_witness.py`, `wp08fix/errbar.py`, review `wp08rev/blast11.py`, `wp08rev/witness_probe2.py`): nested_outer blast 11: 1215 and 740 lie 2.8e-6 and 1.0e-6 inside 38 on s(0.0), on 38's own stable segment (l_s 1.93e-2, ratio 5684), error bar 3.4e-6 on all three; the cross-frame order inconsistency one blast earlier on the same pair is 4e-6 (725 = M(439) sits between M(37), M(38) in frame 1 while 439 lies outside (37,38) in frame 0); the real map cannot settle it either: M(38), M(1215), M(740) all land within 6e-6 of crossing 1209, not of 44 = M(38) by the table, and M(37) misses 43 by 1.1e-5 (the outer tangle is left with no reference, 13 -> 9 classes). k32 blast 11 (30 -> 29 classes): `(93,94)` lost to 1461 at 1.55e-5 from an end, error bars 0.8e-5..2.1e-5. k32 blast 13 (29 -> 22): `(79,80)` lost to 4334 at 1.1e-6, `(449,450)` to 4031 at 2.2e-6, error bars 2.3e-6..4.5e-6. Only k32 blast 10 (34 -> 30) is clearly definitional: `(217,218)` lost to 647 at 1.25e-2 on the stable branch and `(218,219)` to 647 at 2.4e7 on the unstable one, against error bars of 2e-6. Under D1 a decrease is a bug unless root-caused: these stay open as DEC-22, which asks first whether a witness inside the placement error bar may kill a reference at all. p3 at `area_cutoff = 1e-5`, blast 7 (84 -> 64 holes): `f(43)` lies 0.0018 from crossing 50 and 0.0037 from 51 on a 0.017 segment (ratios ~5), an under-resolved image at the blasting cutoff (NUM-5 / DEC-14, WP12). p3 blast 23 (44 -> 23 holes in the report): not reproduced; main at 1e-7 runs 24 blasts at 44 holes / 25 classes. Regression `test_deep_reference_survives_noise`.

**WP8b (D12): root cause = an unbalanced N.** Every decrease WP8 left happens exactly when N holds crossings deeper than the deepest registered point of the pip's own orbit. k32 (E11 build): the deepest crossing IS `f^m(q0)` at every blast 0..9; at blast 9 the bridge at `f^10(q0)` (stable cdist 4.3e-7, under the blast's `min_separation = 1e-5`) is dropped as too close and never fed again, so the pip orbit stalls at `f^10(q0)` (c_u 2.41e9) while the lobes away from the anchor keep being iterated (deepest crossing 1.43e10 at blast 10); every k32 witness lies beyond (blast 10 at c_u 1.2e10, blast 11 at 8.5e10, blast 13 at ~3e12). nested_outer: balanced through blast 9, stalls at blast 10, the blast-11 witnesses (1.7e10, 1.36e10) lie beyond. New witness: the inversion series decreased at blast 4 on base (9 -> 6 classes, holes 9 -> 5). The inversion law cases were unbalanced from the start (11 grown steps leave a partial lap on u(0.1) past `f^5(q0)` carrying their only reference). So DEC-22's "witness inside the placement error bar" question is moot: under D12 no such witness is a crossing of N. Rule (history-free: registry, iterate table, pip): the trellis of the pip `T(q0) = (W^u[P, f^m(q0)], W^s[P, q0])`, `m` the last registered forward iterate of `q0`; the unstable branch of `f^i(q0)` is cut at `c_u(f^i q0)` for `i = m-k+1..m` (`Trellis.cut_at_pip_orbit`, with D10's stable cut). Measured (`wp08b_impl/replay.py`): k32 seed classes 34 at every blast 8..13 (base 34, 34, 30, 29, 29, 22); nested_outer 12, 12, then 13 at blasts 2..11 (base 13 -> 9 at 11); inversion series 6, 7, 8, 9, 9, 9 at blasts 0..5 (base 1, 9, 9, 9, 6, 6). Law-case changes: only the three inversion cases (`inversion` 1 -> 6 classes, holes 1 -> 8, references `(23,24)` -> two pairs, the D9 gap on s(0.1) gone, two classes active; `inversion_one_blast` 9 -> 7; `inversion_two_blasts` 9 -> 8; every law passes on all three). Costs (for the author): discovery stops at the depth the pip's lap is registered to, which the blast's `min_separation` sets: k28 E11's blast-11 switch (15 -> 19 classes on base) is no longer reached (15 from blast 8), k32 and nested_outer freeze from blast 9 / 10; a smaller separation reaches deeper. Monotonicity of the class count in `m` is measured, not proven; the strict series law `series_class_count_never_decreases` (now on `k32_series` 10 blasts and `nested_outer_series` 11, `inversion_series` 4) is the guard. The CLOSED alphabet can still drop by one across a witnessed reference switch inside the balanced depth (k32 blast 8 -> 9: 44 -> 43, seed 34 constant); D1's laws count the seed.

**Evidence:** `stress/numerics/report/report.tex:1247-1252, 1285-1293`; `stress/numerics/report/numbers.tex:560-564, 770-774`; E11 k32/nested_outer series; memory `pseudoneighbor-blast-semantics`

### BCG-6: The class table is not quotiented by the map or the orbit shift; one class per hole-ladder rung, and period-3 classes come in shifted triples
**Status:** closed (WP6, decided D2: no quotient by the map or the orbit shift, no transport-chain collapse; what remains is one class per hole-chain bridge, WP8's BCG-4) | **Severity:** medium | **Related:** BCG-1, BCG-4, DEC-2 | **WP:** WP6

A ladder `E_n→E_{n+1}` produces a class `a_n` for each rung. Transport chains (`a1→a2→…→a20`, single-token rules) could be shown as one symbol plus a delay. Period-3 words run `A0→A1→A2→A0`. The growth is bounded by the hole-chain depth (about 20 blasts to close), but it is linear by construction.

**WP6:** with the minimal alphabet the per-blast growth is gone (BCG-1); p3 keeps 23-25 alphabet classes for 26-28 kept bridges after its reference switch, one per hole-chain bridge class, bounded by the hole chain (BCG-4, WP8). No code change.

**Evidence:** class-growth investigation §4.5; chalkboard §4; E11 p3 blast-20 `canonical_rules`

### BCG-7: Class identity is not stable across blasts: positional `element_id` re-indexes, so every name and letter churns
**Status:** closed (WP8: letters by element identity; names positional by design) | **Severity:** medium | **Related:** BCG-4, DEC-9, TEST-2 | **WP:** WP8

A new hole renumbers every element. None of the 9 class names from blast 0 survive p3 at 6 blasts, the words change letter every blast, and k28 letters reach `cz`. The report did not check whether the unresolved class is the same class across blasts. The chalkboard identity of an element is its bounding holes `(origin, iterate, side)`.

**WP8:** an element is identified across recomputes by its BOUNDING CROSSINGS, `PartitionInterval.identity = (branch_key, side, lo_id, hi_id)`: each branch is cut at every hole bridge's endpoints, and registry ids persist across the recomputes of one session (C10, refining C2). The `(origin, iterate, side)` tuple of C2 does not survive a switch: across p3 blast 5 -> 6, 3 of 9 classes keep their span identity and 0 their hole-tuple identity (k28 blast 10 -> 11: 10/15 vs 0/15); between switches the span identity keeps every class (p3 6 -> 7: 23/23, positional ids 22/23). The session alphabet letters a class by the pair of its elements' identities (`Topology.identity`), so on `p3_series` the 3 classes the base re-lettered at blast 6 keep their letters (law `series_letters_follow_element_identity`). `ElementRef` stays the within-recompute handle (putting the identity inside its equality broke 58 unit tests built on hand-made refs, for no gain), and element NAMES (`R_(o.b;i)^j`) stay positional.

**Evidence:** `.refactor/knowledge/topology.md` (OPEN); `src/tanglepack/loom/TangleSession.py:375-382`; E11 k28 `letters_final`; `stress/numerics/report/report.tex:1478-1479`

---

## Intent vs code (chalkboard review)

Cross-references, catalogued elsewhere: loop folding by preimage (BCG-2), element identity (BCG-7), window reference rule (PK-2), inversion anchor (INV-1), λ_u lap factor (INV-2).

### INT-1: The word is the shortest dual-graph walk between midpoint landings, not the homotopy class of the registered image chain
**Status:** closed (WP7) | **Severity:** high | **Related:** INT-8, PK-1, DEC-3 | **WP:** WP7

**WP7:** under D3 a disagreement between the walk and a registered image is a construction or evidence fault, and it now says which. The WARNING is one per disagreeing member, and `ImageEvidence.diagnosis` is an `EvidenceDiagnosis`: the first differing token, the walked and registered occurrences there, and a cause. The cause is `"stable end"` when two consecutive crossings of the registered chain on one stable branch share a crossing sign, else `"walk"`. It still counts against `is_reliable`. The measured case is k28 at 6 / 7 blasts: class `a` member `(44,45)` has image chain `[(74,75),(75,76),(76,77)]` with three consecutive crossings on s(0.0) of sign −1. The image arc passes around the end of the trimmed stable manifold, its crossings with the trimmed-away part are not registered, and the pairing shifts by one. Resolution ratios are 30-8400, so this is not under-resolution: the walk is right and the evidence is incomplete. Regression `test_evidence_diagnosis`.

The homotopy class of f(bridge) relative to the holes is determined by its crossing sequence. That is `trellis_itinerary`, which is used only for singletons and as "evidence". The answer comes instead from `land_element` midpoint probes plus `shortest_walks`. Shortest equals homotopic only where the dual graph is tree-like; elsewhere the result is `ambiguous`, `unreachable` or virtual `new*`, and the code admits "the shortest walk can be topologically wrong where nothing is registered".

**Resolved by D3 (2026-10-07): rejected.** The shortest walk stays the definition and must always exist; the registered chain is the cross-check and the word for singletons. The audit's original proposal is kept for the record: the trellis itinerary of `_image_chain` becomes THE word, tightened by cancelling a a⁻¹. Classes with no registered image become unresolved, and a loom `resolve_classes()` iterates one representative per class. That would delete the walk search, the ambiguity handling and `registered_evidence/verified`, and shrink `PartitionFamily`/`MinimalTrellis`.

**Evidence:** `src/tanglepack/topology/DualWalk.py:273, 350, 672, 827`; `src/tanglepack/topology/SymbolicDynamics.py:583, 1094, 1169`; `src/tanglepack/topology/PartitionFamily.py:684, 900`; `src/tanglepack/topology/MinimalTrellis.py:243-416`

### INT-2: A hole is punched on the whole containing bridge, not on the backward-image sub-arc
**Status:** closed (WP4, by design: whole containing bridge, evidence in StablePartition Dev Notes) | **Severity:** medium | **Related:** BCG-4, PERF-4, DEC-1 | **WP:** WP4

**WP4:** kept, by design. Instrumenting every containing lookup on the ten law cases, the three deep p3 runs and four deep E11 states (600+ lookups) finds every span a STRICT sub-arc of its container, none needing slack, and at most one of its two ends linked by the iterate table: the other end is a crossing of the unstable branch with `f^-1(W^s)`, beyond the trimmed stable manifold, so the trellis has no node to cut the sub-arc at and the region it bounds is, in the trellis, the region beside the containing bridge. Chalkboard stage 5 says the same ("that bridge gets the hole"). The container is now an exact bisect (PERF-4).

Every ancestor bridge becomes a hole bridge and every ancestor endpoint a partition boundary (42 holes → 85 elements). `_containing_bridge` scans with `rtol=1e-3` slack and "tightest span wins"; bridges tile a ray, so a bisect would be exact.

**Evidence:** `src/tanglepack/topology/StablePartition.py:1771-1800`

### INT-3: Hole openings and partition rows use the GEOMETRIC row while classes, dual graph and walks use the COMBINATORIAL row
**Status:** closed (WP4) | **Severity:** medium | **Related:** INT-4, INT-5, NUM-1, NUM-2, DEC-8 | **WP:** WP4

`_hole_openings` uses `_row_at`, `_stable_frame`, `_bridge_end_geometry`, `_arc_side_of`, nearest-node tests and special trimmed-end/anchor/`inward=True` branches, about 450 lines. Sign arithmetic decides the same thing: `row = row_of_end`, `outward ⇔ (sign>0) == (side==left)`. `measured_side`, `SideDiagnosis` and `RowDiagnosis` exist only to diagnose that geometry.

The two readings disagree on coarse polylines. E11 nested_inner with 2 iterations, blasts 16-20, fails "approaches its two defining crossings from opposite sides" while the combinatorial rows agree (2-node polylines at about 0.2°).

**Evidence:** `src/tanglepack/topology/StablePartition.py:714, 755, 790, 916, 982, 1086-1475, 1236-1285, 1476-1582`; `stress/numerics/e11_blast_depth.py:76-95`

**WP4:** openings are sign arithmetic (`outward iff (crossing_sign > 0) == (bridge_side == "left")`, row `row_of_end`), and so is the I2 assertion (equal signs on one stable branch). The polyline readings are flags (`geometry_flags`, one WARNING per kind per `punch_holes`). The geometry was not only missing at depth but wrong: on p3 at 6 blasts and 15-16 steps 4 of 42-44 holes opened differently (a hole opening nothing, its twin opening `anchorward` at both ends), and at 1e-5 18 of 84-88. Law cases are unchanged; the deep p3 partitions gain 3 elements and lose the unreachable class (PK-1).

### INT-4: `_hole_openings` uses the nearest vertex of the whole polyline (the same defect already fixed for propagated hole sides)
**Status:** closed (WP4; not latent: it misread the deep p3 runs and inversion x period 2 at 10 steps, and went with the geometry) | **Severity:** low | **Related:** INT-3 | **WP:** WP4

The same kind of misread broke p3 hole sides at 15 steps / 6 blasts.

**Evidence:** `src/tanglepack/topology/StablePartition.py:~1480-1500`; `.refactor/knowledge/topology.md` ("Latent (believed)"); `.refactor/decisions.md` 2026-10-02

### INT-5: `_arc_side_of` returns None for a 2-point span slice, giving hole `bridge_side` None with a misleading DEBUG
**Status:** closed (WP4: observed off the law cases, 19 of 26 propagated holes of inversion x period 2 at 10/10 steps; the measured side is read at the nearest segment, which answers on two points, and it is a flag now) | **Severity:** low | **Related:** INT-3 | **WP:** WP4

Fix: widen the slice to at least 3 points.

**Evidence:** `src/tanglepack/topology/StablePartition.py:794-795`

### INT-6: The iterate table is filled by a 5% relative-cdist heuristic, not registered exactly when a crossing is created
**Status:** closed (WP2) | **Severity:** medium | **Related:** BCG-3, NUM-8, SIMP-3 | **WP:** WP2

`infer_iterates` uses `cdist_rtol=0.05` and has no map gate (a 1e-3 phase-space gate exists only on the endpoint path). Deep blasts put real crossings much closer than 5%, so two crossings can be mis-linked, and pips, pairs, holes and chains all trust the link as exact. The topology re-implements the scaling fallback three more times.

The machine already links each mapped point to its image, but refined points cache only a preiterate. So exact registration is a numerics change, not a tweak.

**WP2:** each crossing records both brackets (`unstable_segment`, `stable_segment`); y = M(x) iff y is on the advanced branches and both its brackets meet x's carried one step by the manifolds' own `stretch_param`, uniquely, confirmed by the real map within the image of x's unstable segment plus y's own (one batched call per pass). No 5% window, no absolute gate. No mis-link was ever observed; the old rule's harm was missed links. Topology's own matchers (`_iterate_intersection`, `_resolve_reference_end`) are gone.

**Evidence:** `src/tanglepack/numerics/IterateInference.py:119, 151, 436`; `src/tanglepack/numerics/ManifoldMachine.py:372-374, 878-899`; `src/tanglepack/topology/Pseudoneighbor.py` (`_iterate_intersection` `match_rtol 0.05`, `:491`); `StablePartition._backward_endpoint`

### INT-7: Crossing identity is re-derived by tolerance matching on every `compute_intersections`; ids are not reproducible
**Status:** closed (WP8) | **Severity:** medium | **Related:** NUM-8, API-1 | **WP:** WP2, WP8

The registry is replaced and `reindex_from` re-matches by `cdist_tol` (absolute 1e-6). `Tangle` keys segment edges on `id()`, so ids differ between builds in one process. `ResonanceZone` re-finds ids by tolerance (`_live_id`, `_anchor_id`), and `reindex_from`'s remap is never consumed. Crossing identity is really combinatorial: segment pair plus ordinal along each ray.

**WP2:** exact iterate registration does NOT make `reindex_from` unnecessary: `compute_intersections` replaces the registry, and only the rematch gives a re-detected crossing its old id. Every rematch on every builder is a bit-identical re-detection, so it stays correct under the relative dedup. Ids are not reproducible because `Tangle.resolve_crossings` registers in the iteration order of a `set` of segment-id pairs (`Tangle.py:265`); sorting by `(keys, unstable_cdist)` first is ~3 lines but renumbers every crossing and the id-based fixture facts, so it belongs with WP8's identity work.

**WP8:** ids WERE reproducible between processes (the first build of a fresh process) but not between two builds in one process: `Tangle._ids` is a class-level counter, so the segment ids of a second build carry an offset, and the set of segment-id pairs iterates in another order (`build_period3(unstable_steps=16)` twice: digests `649f51b74e`, `452e48c6b3`, the same two in a second process). So the WP7 coder (two builds in one pytest process) and reviewer (one build) were both right. `compute_intersections` now registers along each unstable branch in unstable-cdist order (then by stable branch and cdist): two builds in one process give one digest, and nothing but the ids moves (every law case, deep run and series has the same id-free signature). Regression `test_registry_ids_are_reproducible`. Left for WP11: `reindex_from`'s returned remap is unused; `ResonanceZone._live_id`.

**Evidence:** `src/tanglepack/numerics/IntersectionRegistry.py:241-330, 863`; `src/tanglepack/topology/Trellis.py:24-30`; `docs/dual_graph_symbolic_dynamics_plan.md:324`; `.refactor/knowledge/testing.md` ("Registry ids are NOT reproducible")

### INT-8: A landing span that straddles a cut still feeds a word (midpoint owner, `contained=False`) instead of being unresolved
**Status:** closed (WP7) | **Severity:** medium | **Related:** INT-1, INV-6 | **WP:** WP7

**WP7:** a straddling landing still resolves to its midpoint owner, so the walk exists (D3), but it counts against `is_reliable` (a construction inconsistency). No law case straddles. The only straddles measured are off the law cases, on the under-resolved `p3_coarse_blasts`: `R(2.0)#1 -> (0, c(47))` covers four s(0.0) elements because crossing 40, an iterated cut, lies inside the image.

Under a correct partition, a homotopy element's image lies inside one iterated element, so a straddle signals a partition or cdist-estimate problem. The current output is a wrong word with a WARNING. Strict xfail on inversion.

**Evidence:** `src/tanglepack/topology/DualWalk.py:270-380`; `tests/cases.py:812-815`; `CLAUDE.md` ("Landing")

### INT-9: Pairs with no spanning Bridge object drop holes and image bridges; the trellis has no notion of an "unstable arc" independent of Bridge objects
**Status:** closed (WP7; filling escape-lobe images deferred) | **Severity:** medium | **Related:** NUM-16, BCG-1 | **WP:** WP7

**WP7: root cause.** Over the law cases every consecutive unstable pair with no Bridge has its two crossings on DIFFERENT computed curves (`Intersection.unstable_manifold`): k28_one 1, k28_two 2, nested_inner_blasts 4, inversion_one 2, inversion_two 5, and 0 on k10 / p3. Each is the uncomputed image of a bridge the blast did not iterate (outside the zone; e.g. inversion_two_blasts `(32,33) = f((26,27))`, where `(26,27)` was never iterated), or the gap between the grown tip and a blast image. `Trellis.is_computed_arc(lo, hi)` says whether a Bridge spans the pair or both crossings lie on one computed curve (an anchor lies on every curve rooted at it). The dense arrangement draws an unstable arc only there; a gap's two ends dangle. All seven arrangement strict xfails pass (`arrangement_euler` x `inversion_two_blasts` and the six NUM-16 entries). Law-case words, matrices and reliability are unchanged. Dense region counts drop on the blasted cases: k28_one 5 -> 4, k28_two 9 -> 7, nested 33 -> 31, nested_inner 37 -> 33, inversion_one 16 -> 15, inversion_two 20 -> 17. The `punch_holes` "no spanning bridge" WARNING on inversion_two_blasts is the `+1` appearance `(32,33)` of reference `(26,27)`, whose own bridge escapes the zone. That is D9's category (an image not computed), not a construction error; the WARNING stays and now names the cause, and so does the minimal trellis's. Computing those images means iterating escaping lobes (WP9 / WP10). The virtual classes this issue listed were PK-1's closure gap.

After a blast there is no Bridge between a parent crossing and a blast-child crossing:
- `punch_holes` skips the pair, and propagation restarts from the most-backward punched pair ("open gap in the punching layer").
- `MinimalTrellis` skips the image pair.
- Virtual classes appear.
- Zero-area regions form (NUM-16).

`BridgeIterator` also creates partial or duplicate Bridge objects over one arc.

**WP3 evidence:** on `inversion_two_blasts` the DENSE arrangement draws unstable arcs between consecutive blast-image crossings no bridge spans (the gaps between images of different parents are not computed manifold), and `V - E + F = 0`; dropping exactly those arcs gives 2 (probe `dbg6.py`, 1 and 2 blasts). The strict xfail `arrangement_euler` x `inversion_two_blasts` is re-tagged here (WP7). Dropping uncut arcs to stubs would also move the k28 blasted arrangements.

**Evidence:** `src/tanglepack/topology/StablePartition.py:190-196, 239-244, 428`; `src/tanglepack/topology/MinimalTrellis.py:280, 323`; `src/tanglepack/topology/SymbolicDynamics.py:494`; `src/tanglepack/numerics/BridgeIterator.py:141-156`; `docs/regions_bridge_identity_and_cleanup_plan.md:196`

### INT-10: Face nesting is decided by geometry (point-in-polygon and signed area), so OuterFaceError fires on slivers
**Status:** closed (WP10) | **Severity:** medium | **Related:** NEST-1, NUM-1 | **WP:** WP10

`_merge_faces`/`_innermost_container` test every node with `point_in_polygon` and signed areas. E11 k28 at 1e-6, blast 16: one component has three positively wound faces, which a planar subdivision cannot have; a self-intersecting sliver is suspected. The rotation system already determines the faces, and the chalkboard needs one containment test per fixed point. There is also a "component with no closed face" WARNING.

**Evidence:** `src/tanglepack/topology/DualGraph.py:597-604, 771, 808-810, 842-875`; `stress/numerics/report/report.tex:1667-1673, 1912`; `stress/numerics/report/numbers.tex:412-416`

**WP4 evidence:** with the openings and row check on the signs, the cells that used to stop in the partition now stop here: the E11 2-iteration cells (p3 at 1e-7 from blast 15, nested_inner from blast 17) and k28 at 1e-5 / 1e-6 blasts 16-17.

**WP10:** a component's outer face is read off the rotation system (`DualGraph._outer_face`: the face whose sector at the component's easternmost polyline vertex holds due east), never an area; nesting stays one point-in-polygon test per component. `OuterFaceError` is gone; `test_outer_faces_do_not_read_areas` fakes every traversal area (all positive, all flipped) and the face nodes of p3 and nested do not change (base raises). Same outer faces as the area rule on every law case and p3 16u. k28 at 1e-5, 17 single blasts, no longer reproduces the error on main (D10 / D12 cut the trellis at the pip orbit); the E11 2-iteration cells were not replayed (WP12). The pip restore cost that made those deep recomputes slow is fixed too: a rebuild tests only the remembered pip (k28 at 1e-5, 12 blasts: recompute 4.6 -> 0.8 s).

### INT-11: At depth a registered image ends outside its own element's landing span, or crosses the window one element away from the walk
**Status:** open | **Severity:** medium | **Related:** INT-1, INT-8, DEC-23 | **WP:** WP12 (found by WP9b)

Evidence disagreements left on non-law deep builds once the walk runs between the landing spans under the non-crossing rule (D13); none is a tie. Letters below are those of a session that computes the dynamics once, at the given blast (a series that recomputes at every blast letters differently). (1) k28 at 9 blasts (class `l`, member `(246, 247)`) and k32 at 9 blasts (`o`, `(298, 299)`): the member's registered image ends OUTSIDE its own element's landing span at the target end, so the landing and the partition disagree at depth (a landing or partition inconsistency, not a walk one). (2) k32 at 9 blasts, `k` and `s`: the registered image crosses the window one element away from the walk, both ends inside their spans, and both curves pass the marker law. These keep k28 at 9-13 and k32 at 9-10 blasts unreliable. `inversion_period2` is unreliable through `e`, a singleton landing with no registered chain (pre-existing, unrelated).

**Evidence:** WP9b plan F7 (probe `evcheck.py`); WP9b probe `f7.py` (`cases.SERIES["k28_series"]` / `["k32_series"]` blasted to 9, the chain ends against the landing spans)

---

## Unnecessary complexity / simplification

### SIMP-1: The orientation side-flip rule is coded at six sites and cannot be reached while det<0 is rejected
**Status:** closed (WP11) | **Severity:** low | **Related:** INV-11, DEC-10 | **WP:** WP11

One helper (`image_side(side, n)`) would do, or delete the flips until det<0 is admitted.

**WP11:** one helper, `side_after` (see INV-11); `image_of_element`'s flip went with SIMP-2, and no `OPPOSITE_SIDE[...]` stands next to an orientation test any more.

**Evidence:** `src/tanglepack/topology/StablePartition.py:563, 596-598, 1287-1300, 1409-1412`; `src/tanglepack/topology/Trellis.py:292, 1423-1424`; `src/tanglepack/topology/DualGraph.py:934-940`; `src/tanglepack/topology/DualWalk.py:317`

### SIMP-2: Two element-image implementations with different tolerance rules
**Status:** closed (WP11) | **Severity:** medium | **Related:** INT-8, SIMP-3 | **WP:** WP11

- `Trellis.image_of_element`: cover plus `_snap_to_partition_boundary` with `SNAP_RTOL`/`SCALING_RTOL`. It feeds `StableNode.images`, which only plots read.
- `land_element`: midpoint owner.

Same rule, possibly two answers.

**WP11:** `Trellis.image_of_element` (with `_snap_to_partition_boundary` and `SNAP_RTOL`) fed only `StableNode.images`, which nothing read (not the plots, not the walk); both are deleted, with `DualGraph._link_element_images`. `land_element` is the one element image.

**Evidence:** `src/tanglepack/topology/Trellis.py:42-105, 128, 143, 1315-1440`; `src/tanglepack/topology/DualWalk.py:270-380`; `src/tanglepack/topology/DualGraph.py:932-956`

### SIMP-3: Cdist scaling is implemented several ways; `scale_cdist` is unused
**Status:** closed (WP2) | **Severity:** low | **Related:** INV-2, INT-6 | **WP:** WP2

- `scale_cdist` divides by the unstable β for stable cdists. Verified 2026-10-07: its only reference is its own docstring.
- `_scaled_image_cdist` uses `per_step_beta(stability)`.
- `_backward_endpoint` divides by β.
- `Pseudoneighbor` uses `lambda_u**m`.
- `land_element` scales the open end by `per_step_beta('stable')`.

`Trellis.lambda_u` and `TrellisBranch.eigenvalue` are extra accessors; the latter has 0 uses and a wrong docstring. `image_cdist` should be the one path.

**WP2:** a crossing is moved only by `Trellis.image_cdist`, a bare distance with no crossing only by `FixedPoint.per_step_beta` (`_backward_endpoint`, `land_element`'s unbounded end); nothing in topology reads an eigenvalue (`scale_cdist`, `lambda_u`, `TrellisBranch.eigenvalue`, `IntersectionRegistry.on_interval` / `_get_lambda_u` deleted; the orientation fallback reads `FixedPoint.preserves_orientation`).

**Evidence:** `src/tanglepack/topology/Trellis.py:679-763` (:690, :725); `src/tanglepack/topology/StablePartition.py:1744-1768`; `src/tanglepack/topology/Pseudoneighbor.py:491, 551-576`; `src/tanglepack/topology/TrellisBranch.py:72-84`; `src/tanglepack/topology/DualWalk.py:330-334`

### SIMP-4: The 2-D same-point collision test is duplicated in StrongPip and Pseudoneighbor; three branch-cycle wrappers sit over one memo
**Status:** closed (WP11) | **Severity:** low | **Related:** INT-6, SIMP-12 | **WP:** WP11

The collision test exists only for unlinked crossings and would vanish with an exact table. The wrappers are `forward_stable_branch_cycle`, `forward_unstable_branch_cycle` and `_branch_position_map`.

**WP11:** one 2-D collision test, `Trellis.same_point(u, s, u0, s0, rtol)` (relative only, C15), used by StrongPip and twice by Pseudoneighbor; the three branch-cycle wrappers and their exports are gone (callers read `FixedPoint.branch_position_map`).

**Evidence:** `src/tanglepack/topology/StrongPip.py:92, 260-280`; `src/tanglepack/topology/Pseudoneighbor.py:114-128, 609-610`

### SIMP-5: PartitionFamily keeps a local copy of `_build_intervals` ("concurrent edit" reason is stale)
**Status:** closed (WP11) | **Severity:** low | **Related:** DOC-3 | **WP:** WP11

**WP11:** `StablePartition.intervals_from_marks` is the shared tail; `_build_intervals` adds the branch's ends (outermost crossing, anchor, the C1 branch) and calls it, and `PartitionFamily` imports it instead of its copy. The marks round trip test stays the guard.

**Evidence:** `src/tanglepack/topology/PartitionFamily.py:54-67, 518`; `src/tanglepack/topology/StablePartition.py:1968`

### SIMP-6: Two letter generators (active base-26 in loom, inert u..z,uu.. in topology), with fragile collision-skipping
**Status:** closed (WP6) | **Severity:** low | **Related:** BCG-1 | **WP:** WP6

Active letters can reach `u` after 20 classes, and E11 reaches `bw`/`cz`.

**WP6:** one generator, `topology.letter(index, inert=False)`: a bijective numeral over `a..t` (active) or `u..z` (inert). The two series are disjoint by construction, so `inert_letters(table)` (no `used=`) and the active fallback skip nothing, and `symbolic_dynamics` asserts one letter per class. `BridgeAlphabet` uses it; `loom.letter` is gone. Active letters change only past the 20th active class of one session (`u` -> `aa`); no law case gets there.

**Evidence:** `src/tanglepack/loom/BridgeAlphabet.py:10-25, 74`; `src/tanglepack/topology/SymbolicDynamics.py:330-345`

### SIMP-7: `Blast.py` tidy-up
**Status:** closed (WP11) | **Severity:** low | **Related:** NEST-4 | **WP:** WP11

- It hand-rolls the fixed-point normalisation that `TangleSession._resolve_fixed_points` already does.
- `cKDTree` is imported inside the loop, with a dead scipy-missing fallback.
- `BlastResult.fixed_point` is typed `Optional[FixedPoint]` but stores lists.

**WP11:** `blast_zone` normalises its filter with `session._resolve_fixed_points` (`BlastResult.fixed_points`, a list; the mistyped `fixed_point` had no reader), `cKDTree` and `ResonanceZone` are imported at module top (the scipy-missing fallback and `_resolve_zone` are gone), `max_depth_reached` is deleted (SIMP-12).

**Evidence:** `src/tanglepack/loom/Blast.py:95-99, 160, 226-236, 281-289, 373`; `src/tanglepack/loom/TangleSession.py:538`

### SIMP-8: Four or five hand-rolled session caches keyed on `(generation, signature[, pips])` for stages that each cost under 0.2 s
**Status:** closed (WP8) | **Severity:** low | **Related:** API-1, API-2 | **WP:** WP8

`trellis(fp)` and `trellis([fp])` are separate cache entries. One topology object per generation (classes → holes → partition → words) would replace them.

**WP8:** `TangleSession.topology(fixed_points)` holds one `loom.Topology` per `(generation, selection, pips)`: it recomputes pseudoneighbors, holes and the partition of every selected fixed point from its chosen pip (D1), and the classes, minimal trellis, iterated partition, dual graph and dynamics are its lazy stages; the six accessors are one-line delegates. The five hand-rolled caches, `_partition_signature`, the three `_gathered_*` scans and `_annotate_bridge_classes` are gone, and so is the unpartitioned-branch WARNING (every selected branch is partitioned). A selection is a set: `trellis(fp) is trellis([fp])`, fixed points in workbench order.

**Evidence:** `src/tanglepack/loom/TangleSession.py:291-400, 379-388, 442-560, 603-611, 643-651, 677-860, 729-740, 828-839, 1313`

### SIMP-9: ManifoldInitializer open notes
**Status:** partly closed (WP11) | **Severity:** low | **Related:** INV-7 | **WP:** WP11

- Ethan vs Kevin construction.
- A duplicate `_insert_point_geometrically` (it also calls the machine's).
- `get_initial_fundamental_segment` is too long.
- Requests: orient at any iterate, or one direction alone.

**WP11:** one construction is kept and the Dev Note says so (the per-orbit-point seeding never existed in code); `get_initial_fundamental_segment` is one path; `construct_manifold_from_point_list` and its duplicate `_insert_point_geometrically` moved to `ManifoldMachine._manifold_from_points` (the machine no longer builds a throw-away initializer through a lazy import). Open: the two feature requests (orient at any iterate, one direction alone) -> the API plan.

**Evidence:** `src/tanglepack/numerics/ManifoldInitializer.py:9-45, 391, 557, 588`; `src/tanglepack/numerics/ManifoldMachine.py:902`

### SIMP-10: Two stable trims
**Status:** closed (WP11) | **Severity:** low | **Related:** INT-7 | **WP:** WP11

`TangleWorkbench.trim_stable_manifolds` and `ResonanceZone._trim_stable_at`.

**WP11:** one numerics primitive, `TangleWorkbench.trim_stable_at(key, cdist)` (the first node past the root at or beyond `cdist`); `trim_stable_manifolds` and the zone trims call it, and `ResonanceZone._trim_stable_at` (which read the rtree's private segment tables with a `cdist_tol` slack) is gone. The two rules were measured identical on every zone trim of 15 builders before the change.

**Evidence:** `src/tanglepack/numerics/TangleWorkbench.py:2045`; `src/tanglepack/loom/ResonanceZone.py:429`

### SIMP-11: The cartoon bridge-drawing block is duplicated between the line and circle layouts; `plotting.py` is 3275 lines
**Status:** partly closed (WP11) | **Severity:** low | **Related:** — | **WP:** WP11

**WP11:** the node scatter, unified S-curves, bridge arcs and walks of both layouts are one `_draw_cartoon_graph`, parametrised by the shape's `arc` / `cross` / `track`; every drawn artist is byte-identical on four fixtures in both shapes and colourings. Open: splitting `plotting.py` (import churn for TangleSession, scripts and tests) -> a later slice.

**Evidence:** `src/tanglepack/topology/plotting.py:2140, 2360-2375, 2922, 3060-3075`

### SIMP-12: Dead or redundant public API, kept and whitelisted by the coverage guard
**Status:** closed (WP11) | **Severity:** low | **Related:** SIMP-3, SIMP-4, PK-7, API-4 | **WP:** WP11

- `Trellis.scale_cdist`, `_get_lambda_u`, `IntersectionRegistry.on_interval`
- `TangleSession.invalidate_trellises`
- scalar `ManifoldMachine._curvature_area`
- the branch-cycle wrappers
- `IterateTable.forward_depth/backward_chain/all_registered_ids`
- `IntersectionRegistry.unstable_order_array/stable_order_array`
- `Symbol.is_inverse`, `SymbolicDynamics.dynamics_of`
- `TrellisBranch.anchor_coords/eigenvalue`
- `BridgeClassEntry/Table.member_of`
- `BlastResult.max_depth_reached`
- the `reindex_from` remap
- `get_iterable_array`
- the `manifold_a_key is None` branches

WP2 removed `Trellis.scale_cdist`, `_get_lambda_u`, `IntersectionRegistry.on_interval`, `TrellisBranch.eigenvalue` and the `manifold_a_key is None` branches.

**WP11:** deleted `invalidate_trellises`, `BlastResult.max_depth_reached`, `IntersectionRegistry.unstable_order_array` / `stable_order_array`, `IterateTable.backward_chain` / `all_registered_ids`, `Symbol.is_inverse`, `SymbolicDynamics.dynamics_of`, `TrellisBranch.anchor_coords`, the scalar `ManifoldMachine._curvature_area` and its three helpers, `BridgeClassTable.member_of`, the branch-cycle wrappers, three unused `Intersection.synthetic` keywords, `bridge_classes(bridge_ids=)`, the `infer_iterate_table` forwarder (scripts and the stress harness call `infer_iterates`) and the stateless `IterateInference` class (a module function `infer_iterates`); `reindex_from` returns None. The coverage guard's whitelist went with them (TEST-6).

**Evidence:** `tests/conftest.py:85-93`; `src/tanglepack/numerics/IterateTable.py:81, 98, 107`; `src/tanglepack/numerics/IntersectionRegistry.py:792, 800`; `src/tanglepack/topology/SymbolicDynamics.py:167, 787`; `src/tanglepack/topology/TrellisBranch.py:72, 96`; `src/tanglepack/topology/BridgeClass.py:410`; `src/tanglepack/loom/Blast.py:178`; `src/tanglepack/topology/Trellis.py:24-30, 690`; `docs/dual_graph_symbolic_dynamics_plan.md:326`

### SIMP-13: Tolerance and warning density
**Status:** closed (WP11) | **Severity:** low | **Related:** INT-6, INT-7 | **WP:** WP11

- 86 tolerance sites in `topology/`.
- 73 WARNINGs against 38 asserts across the library.
- About 130 lines of `tangle_of`/`min_cdist`/`anchor_outward_key` ordering just to make letters deterministic.

**WP11:** counted by an AST counter over `topology/` (identifier sites matching `tol` / `*_tol` / `rtol`, `cdist_tol` sites, `logger.warning` calls; the catalog's 86 was an unrecorded grep): before 115 / 19 / 56, after 59 / 2 / 51. The decisive measurement: forcing the `cdist_tol` topology and loom read to 0 (or 1e-12) moved no output on any law case, deep run or series, and the suite passed. So (C15) anchors are told by identity, ownership / containment / the pip cut / the fundamental-segment overlap / `interval_at` are exact, and only RELATIVE slack is left (`SCALING_RTOL`, `collision_rtol`); `cdist_tol` remains for the degenerate bridge cut and polyline clipping. Two WARNINGs of `land_element` (an image on a branch other than `advance_key`'s, which the iterate table rules out) are asserts; the cartoon's duplicated warnings went with SIMP-11. The walk's WARNINGs and the class-ordering code (the author's rules) stay.

**Evidence:** chalkboard §8(b); `src/tanglepack/topology/BridgeClass.py:465-560, 659-693`

---

## Resolution & numerics

### NUM-1: The bridge-row AssertionError fires on under-resolved bridges even when its own diagnosis says the combinatorial rows agree
**Status:** closed (WP4) | **Severity:** medium | **Related:** INT-3, NEST-1, DEC-18 | **WP:** WP4

Seen in 39 snapshots across 4 configs: p3 and nested_inner at 1e-5 from blast 6, and the 2-iteration cells. Example: p3 bridge (41,42) has signs +1/−1 and the combinatorial rows agree (both left), but the geometric row is read on a 2-node polyline at 0.16°. The heteroclinic case at 1e-5 asserts with "signs alternate … likely below resolution" (8 nodes, 0.3-0.6°). Asserting here contradicts the 2026-10-07 rule that resolution problems are flags. Relaxing only the geometric reading exposes the hole-side check next.

**Evidence:** `src/tanglepack/topology/StablePartition.py:1269`; `stress/numerics/report/report.tex:235-237, 1674-1690, 1728-1729, 1910-1912, 1995-2002`; `stress/numerics/report/numbers.tex:424-428`

**WP4:** root cause: on a 2-node polyline the "approach displacement" is the vector to the FAR end node. The assertion reads the signs now (it fires on equal signs only) and the geometric row is a flag. Replayed E11: p3 and nested_inner at 1e-5 (one iteration) run all 20 blasts (were AssertionError from blast 6); regression `test_coarse_blasts_rows_by_signs` (builder `p3_coarse_blasts`). The 2-iteration cells now pass the partition and stop at `OuterFaceError` (p3 blast 15, nested_inner blast 17: INT-10, WP10).

### NUM-2: Hole-side AssertionError: carried "left" vs measured "right"
**Status:** closed (WP4; k28 then stops at INT-10, WP10) | **Severity:** medium | **Related:** INT-3, NUM-1 | **WP:** WP4

3 snapshots: k28 at 1e-5/1e-6, blasts 16-17, several iterates of one orbit. Probably measurement on under-resolved sub-arcs; the recommendation is to re-measure on refined sub-arcs.

**Evidence:** `stress/numerics/report/report.tex:240, 1654-1666, 1730-1732, 1910`; `stress/numerics/report/numbers.tex:418-422`

**WP4:** the measured side is the geometric cross-check of the carried one (D6): a disagreement is one WARNING per `punch_holes` (`geometry_flags`), never an assertion; `check_holes_share_bridge_side` asserts the orbit invariant only. The measured side is unreliable for a second reason: the carried point is the reference hole mapped back by the real inverse map, and when the reference lobe is under-resolved the point is not in the true lobe and escapes (inversion x period 2 at 10/10 steps: origin `(3,31)` reaches 1e15 by iterate -11, its reference bridge has 3 nodes). Replayed on the branch, the k28 1e-5 / 1e-6 cells pass blast 15 (the measured-side disagreements are flagged, not asserted; the analyst's prototype counted 14-17 per snapshot) and raise `OuterFaceError` at blasts 16-17 (INT-10, WP10), as k28 1e-6 blast 17 already did; so "k28 runs through" is reached past the partition, not end to end.

### NUM-3: Registered crossings below the resolution ratio appear from blast 9 (period 1) or 15 (period 3)
**Status:** open (flags only; auto-refinement of flagged stretches is open) | **Severity:** medium | **Related:** NUM-4, BCG-4 | **WP:** deferred

Affects 52/53 E11 configs, up to 96 crossings in one blast. Last blast: k28 15, k32 15, nested_outer 21, p3 3. These are the depths where classes explode and topology errors appear.

**Evidence:** `stress/numerics/report/report.tex:1539-1547, 1601-1603, 2003-2005`; `stress/numerics/report/numbers_extra.tex` (NumDeepUnresolvedCrossings*)

### NUM-4: Near-tangency artefact crossings are flagged but still registered and fed to the topology
**Status:** decided 2026-10-07 | **Severity:** medium | **Related:** DEC-15 | **WP:** deferred

At k+1 = 1e-4: 25/153/499 registered crossings at cutoffs 1e-7/1e-10/1e-12. Flagged in 46/48 configs below k+1 = 0.16. The threshold of 3 was calibrated on one map family (p3 margin 7.85). Near-parallel pairs are discarded with a WARNING. A validated decision would need interval or multiprecision arithmetic.

**Evidence:** `src/tanglepack/numerics/Tangle.py:347-356, 470, 1074`; `stress/numerics/report/report.tex:1036-1127, 1899, 2018-2024`; `stress/numerics/report/numbers_extra.tex:303-345`

### NUM-5: Symbolic words (and class counts) depend on the blasting cutoff before crossing counts do
**Status:** decided 2026-10-07 | **Severity:** medium | **Related:** BCG-4, DEC-14 | **WP:** WP12

- Words leave the reference first in 14/15 cutoff cells (p3 @1e-5: words from blast 3, crossings from 6).
- Final p3 classes are 48/41/51/49 at cutoffs 1e-5..1e-8.
- k28 at blast 17: 63 @1e-5, 56 @1e-6.
- The 1e-8 runs stopped earlier.

**Evidence:** `stress/numerics/report/report.tex:1295-1300, 1387-1402, 1456-1472, 2029-2036`; `stress/numerics/report/numbers.tex:508, 518, 614, 624, 702-730`

### NUM-6: Two iterations per blast ≠ two single blasts in period-1 zones
**Status:** open | **Severity:** medium | **Related:** NUM-5, API-1 | **WP:** WP12

Period 3 agrees 10/10. k28 is 4/7 (differs at n = 5, 6, 7), k32 4/7, nested_outer 6/7. Crossing counts agree only up to n = 3-5. Suspected cause: different refine and repin moments. The differing pairs have not been examined.

**Evidence:** `stress/numerics/report/report.tex:1403-1416, 1462-1472`; `stress/numerics/report/numbers.tex:522-526, 628-632, 732-734, 838-842`

### NUM-7: Lobe-area excess |A'/A−1| up to 2.29 with 4 sign flips in the two-iteration cells
**Status:** open | **Severity:** medium | **Related:** NUM-1 | **WP:** WP12

Reference cells peak at 0.0514 (median 1.2e-5). The excess falls on the first blasts where the row error fires; it could be a self-intersecting probe polygon or a real area defect.

**Evidence:** `stress/numerics/report/report.tex:1625-1633, 1716-1719, 2037-2040`; `stress/numerics/report/numbers.tex:384-394`

### NUM-8: The registry `cdist_tol` is an absolute 1e-6, so dedup can merge distinct deep crossings
**Status:** closed (WP2) | **Severity:** low | **Related:** INT-6, INT-7 | **WP:** WP2

68% of consecutive k28 gaps fall below it by blast 14, and it caused E11 error 1. The fix is a relative tolerance. The near-ULP vs exact merge dedup question is also undecided.

**WP2:** dedup is relative, `DEDUP_RTOL = 1e-9` (every registry collision on every builder, k28 at 11 blasts included, is a bit-identical re-detection; the closest distinct same-key pair is 4e-3 apart; an anchor matches only an exact 0.0). `cdist_tol` stays as the topology's absolute slack (its ~86 consumers are SIMP-13, WP11). **WP11:** removed from topology: it reads no absolute cdist slack (C15, SIMP-13).

**Evidence:** `src/tanglepack/numerics/IntersectionRegistry.py:85-95`; `src/tanglepack/topology/StablePartition.py:200-206`; `stress/numerics/report/report.tex:1645-1649, 2003-2005`

### NUM-9: `ManifoldMachine.merge_manifolds` crashes on non-contiguous `old_iterated_points` (None head when a bridge reaches the grown tip)
**Status:** closed (WP10) | **Severity:** medium | **Related:** NUM-10, NEST-4 | **WP:** WP10

`build_nested` no longer reaches it, but the bug remains.

**Evidence:** `src/tanglepack/numerics/ManifoldMachine.py:391-404, 516`; `.refactor/knowledge/topology.md`; `.refactor/decisions.md` 2026-10-02; memory `loom-and-subpackages:76-77`

**WP10:** root cause: a bridge sharing its two straddle nodes with a neighbour that was iterated first has its existing images on TWO curves (`IIIIIII.....II`: the grown manifold's tip, then the start of the neighbour's image); `_iterate_manifold` merged them as ONE chain and walked off the first curve's tail into `None` (repro: `build_nested(outer_blasts=0, through="zones")`, two inner blasts, one unfiltered outer blast, bridge `(31,32)`). `ManifoldMachine._image_chains` splits the images where the walk from the first one breaks (images the walk passes over without reaching are the duplicates a merge dropped, skipped as before; "one chain per run of iterated nodes" changed k28 at 6 blasts and was rejected). Regression `test_blast_through_a_shared_straddle` (fails on base). Law cases, deep runs and series unchanged.

### NUM-10: Blast silently counts AttributeError/ValueError iterate failures as "skipped"
**Status:** closed (WP10) | **Severity:** medium | **Related:** NUM-9 | **WP:** WP10

About 5 failures per blast in Phase 6 ("'NoneType' object has no attribute 'cdist'"). The count varies by run and was never examined, so it can hide NUM-9.

**Evidence:** `docs/regions_bridge_identity_and_cleanup_plan.md:339`; `.refactor/knowledge/loom.md`

**WP10:** with NUM-9 fixed no law case, series or deep run skips a bridge; the tolerance, `strict`, `BlastStep.skipped` and `BlastResult.skipped` are gone and every iterate failure propagates (`test_every_iterate_failure_escapes`, three error types).

### NUM-11: Escape to infinity: `grow_n_times` has no notion of escape, and PointBudgetExceeded leaves a session that cannot be resumed
**Status:** closed (WP10: no recipe escapes into its budget; measured) | **Severity:** medium | **Related:** INV-4, TEST-1 | **WP:** WP10

One step multiplies points by up to 40× once |x| > ~100. p3 at 18 steps and inversion at 7-8 steps hit the budget. Recipes and tests still use fixed step counts instead of `grow_until_bounded`.

**Evidence:** `stress/numerics/report/report.tex:404-452, 899-934, 2003-2004`; `stress/numerics/report/numbers_extra.tex:63-131, 259-265`

**WP10:** no recipe or test grows into its point budget (`escape.py`): the escaping arms reach |x| 2574 (k10 u, 1503 points), 154 (k28, 4398), 235 (nested / heteroclinic p1, 6071), 478 (inversion, 1182); p3 and inversion_period2 stay inside |x| <= 10. The escape stop (`grow_until_bounded`), `escape_radius` and the point budget landed 2026-10-07; switching the frozen builders would change law cases for no measured gain. The stress recipes that hit the budget (E3/E7) are rerun by WP12.

### NUM-12: Growth after a blast recuts the bridge set and clears the `iterated` flags (blast genealogy lost)
**Status:** closed (WP10b) | **Severity:** low | **Related:** BCG-3 | **WP:** WP10

**Evidence:** `src/tanglepack/numerics/TangleWorkbench.py:1071-1080`

**WP10 slice 1 (plan):** `compute_intersections` re-indexes MANIFOLDS only, so every crossing born on a blast image curve vanishes from the new registry, and `rebuild_bridges` cuts only manifolds. The fix keeps the blast's image curves, re-indexes them query-only against the stable index inside `compute_intersections` and lets `create_bridges` cut them; flags carry by id under `preserve_ids`. Medium size, touches the recompute every recipe runs: slice 2. The WARNING in `_recut_bridges_after_growth` stays.

**WP10b (C20):** worse than catalogued: ANY recompute after a blast, with no growth at all, lost every crossing on a blast image curve, its bridges and the blast frontier, and class counts dropped (plain `compute_intersections(preserve_ids=True)` + `rebuild_bridges()`: k28_two_blasts 10 -> 4 crossings, 7 -> 3 bridges, classes 4 -> 3, reliable -> unreliable; nested 13 -> 12, inversion_two_blasts 8 -> 6, heteroclinic 57 -> 53; deep: p3_6_blasts 23 -> 9, k28_deep_blast_sliver 422 -> 4 crossings). Zone (re)definition and `ResonanceZone.restore` after a blast hit it too; no law case recomputed after a blast. Fixed as WP10 proposed: the workbench keeps every blast image curve (`_image_curves`), every `compute_intersections` re-indexes it query-only (a curve that adds no segment of its own lies on a grown branch and is dropped; `Tangle.add_manifold` returns the count), and `create_bridges` cuts it like any curve (one crossed fewer than twice is a partial piece, the curve itself, so its `iterated` flag survives). `_recut_bridges_after_growth`, its flag clearing and its WARNING are gone (`grow_until` calls `rebuild_bridges`). Measured: 0 crossings lost by position after a recompute, an unstable step or a stable step (k28: 10 / 10 / 30 crossings); blast / recompute / blast equals blasting straight (k28 and p3, id-free `describe()`); on the 11 law cases, `k28_deep_blast_sliver` and `p3_6_blasts` a recompute + rebuild leaves the bridges (ids and flags), holes, class table and exact `describe()` unchanged. Law `recompute_keeps_the_topology` (new) and `recompute_preserves_ids` (every id back) fail on base on k28_one_blast, k28_two_blasts, nested, nested_inner_blasts, inversion_one_blast, inversion_two_blasts and heteroclinic. Cost: the recompute now carries the blast (k28 after 9 blasts 1.4 -> 3.8 s; base had thrown it away). Found, not fixed: growing the STABLE manifold after a zone trim un-trims the zone (`grow_n_times` resets every tail of the grown stability; k28 two blasts 10 -> 30 crossings, classes 4 -> 4: topology unaffected by D10, the zone object goes stale): a zone-semantics question for WP12 / the D8 API plan.

### NUM-13: `grow_until_arclength` stops on canonical distance, not arc length
**Status:** closed (WP10) | **Severity:** low | **Related:** — | **WP:** WP10

**Evidence:** `src/tanglepack/numerics/TangleWorkbench.py:21-26, 823`

**WP10:** `grow_until_arclength` stops on the polyline arc length of the grown branch (`test_arclength_stops_on_arc_length`: on base the tail's cdist passes the length a step before the polyline does).

### NUM-14: `iterate_bridge` does not re-assert Bridge invariants
**Status:** closed (WP10) | **Severity:** low | **Related:** INT-9 | **WP:** WP10

**Evidence:** `src/tanglepack/numerics/ManifoldMachine.py:80-90`

**WP10:** every bridge enters through `TangleWorkbench._register_bridge`, which asserts a non-anchor end lies on the bridge's own branch inside its node span (to `DEDUP_RTOL`) and an anchor end comes first and is the bridge's own periodic point. Holds on every law case, deep run and series; it fired on base for the separator builds of NUM-19 (`test_a_bridge_off_its_branch_is_refused`).

### NUM-15: Small-k recipes: MaxIterationsReached for k ≤ 0.545, and 26/50 b=1 saddles give trivial tangles (fixed 8-step recipe)
**Status:** closed (WP10) | **Severity:** low | **Related:** NUM-11 | **WP:** WP10

Recipe depth should scale with 1/ln λ_u. These are first-campaign numbers.

**Evidence:** `stress/numerics/report/report.tex:823-826, 864-871, 883-889`; `stress/numerics/report/numbers_extra.tex:199-213`

**WP10:** `grow_until_turnaround` / `_bounded` / `_arclength` / `_intersection` default to `max_iterations=None` = `max(10, ceil(log_stretch / ln stretch) + 1)` with the driver's own needed stretch, never below the old 10, so no recipe changes (`test_default_cap_scales_with_the_stretch`: the k = 0.3 saddle turns around). The fixed 8-step recipe is the stress harness's E6 (WP12).

### NUM-16: Arrangement: a zero-area region after a blast reports a representative point inside its neighbour
**Status:** closed (WP7, root cause INT-9) | **Severity:** medium | **Related:** INT-9 | **WP:** WP10

It is bounded by an unstable arc that no Bridge spans; this is a sparse/full arrangement mismatch at a blast child.

**Evidence:** `tests/cases.py:816-825`; `docs/regions_bridge_identity_and_cleanup_plan.md:196`

**WP10:** verified: `arrangement_regions_disjoint` passes on all 11 law cases, with no KNOWN_ISSUES or NOT_APPLICABLE entry.

### NUM-17: Arrangement `image_of` misses 39 of 86 regions; union-valued images are undefined
**Status:** decided 2026-10-07 | **Severity:** low | **Related:** DEC-17 | **WP:** deferred

None found for k10, unblasted k28 or inversion (first-campaign numbers). `arrangement_image_of`/`preimage` are N/A on k10 and inversion. The hole set is not forward invariant at the computed extent.

**Evidence:** `tests/cases.py:874-880`; `stress/numerics/report/report.tex:205, 722-728, 2023-2027`; `stress/numerics/report/numbers.tex:190-192`; `docs/regions_bridge_identity_and_cleanup_plan.md:338, 342`

### NUM-18: The pseudoneighbor stable-side check walks only registered crossings (asymmetric with the unstable all-iterates check)
**Status:** deferred (author chose no change) | **Severity:** low | **Related:** BCG-5 | **WP:** deferred

**Evidence:** memory `pseudoneighbor-blast-semantics:8`

### NUM-19: Small deferred plan findings
**Status:** closed (WP10b) | **Severity:** low | **Related:** — | **WP:** WP0 → owners (verify in the WP touching each file)

- `_insert_crossing_separator` stale edge bounds. **Closed (WP10):** two or more crossings on one segment are bracketed on the current polyline before cutting (a separator after each previous crossing's lower bracket), and bridge `(a, b)` runs from below `a` to above `b`; before, the separator became the TAIL of `(a, b)` (crossing `b` outside the bridge: WP7-D1, `bridge_endpoints_on_own_branch` on inversion x period 2 at 10 / 10 steps, bridge `(18,19)`) and a third crossing put its separator before the second (cdist decreasing: p3 13u / 16s). The p1(u) x p3(s) recipe (p3 13u / 15s), whose "near-duplicate crossing" was two crossings on one segment, builds now (`test_crossings_sharing_a_segment_are_bracketed`).
- `_validate_saddle` may reject near-parabolic period-k saddles. **WP10b: not a bug.** Exact seeds (k = 10, 2.8, 2, with and without an analytic Jacobian) pass, so `ier != 1` does not fire at a zero residual; period-1 saddles at k = -1 + d pass down to d = 1e-12; on a rotated `diag(1+e, 1/(1+e))` with additive noise n the validator first rejects at `|lambda - 1| ~ n`, i.e. only where the orbit is numerically parabolic and has no growable tangle (NUM-15). No change.
- `bridge_classes` silently accepts a reversed BridgeId. **WP10b: stale since WP11** (`bridge_classes(bridge_ids=)` is deleted: classes read the trellis's own bridges, cut in unstable order; law `bridge_id_in_unstable_order`); closed at the source with one assertion in `_register_bridge` (`test_a_reversed_bridge_id_is_refused`, fails on base).
- The walk and the class decide ownership by different rules. **WP10b: stale since C15 (WP11).** The class reads a crossing's owner (`element_of_intersection`: id at a boundary, exact cdist elsewhere); the walk reads the same table for singletons (`owner_of_intersection`) and the exact `interval_at` for a span midpoint; the only slack anywhere is the relative error at a scaled landing end (C15). No change.

**Evidence:** `docs/regions_bridge_identity_and_cleanup_plan.md:319-320, 335`; `docs/dual_graph_symbolic_dynamics_plan.md:328, 332`


### NUM-20: Heteroclinic iterate links on the p1(u) x p3(s) recipes fail the beta and area laws
**Status:** closed (WP10b) | **Severity:** low | **Related:** NEST-1, INT-6 | **WP:** deferred (found by WP10)

`build_nested(p3_unstable_steps=13, p3_stable_steps=15)` and `(13, 16)` build since WP10's separator fix, but `area_along_iterate_links` and `iterate_link_scales_by_beta` fail on links `7 -> 215` (area 3.1e5 -> 1.06e6) and `5 -> 9`: a `W^u(p1) x W^s(p3)` crossing is linked to a wrong image by the exact rule. Not a law case. Needs its own root cause in `IterateInference` for a `W^u(A) x W^s(B)` crossing.

**WP10b: the root cause was the laws, not the inference.** On 13u/15s (4 heteroclinic links) and 13u/16s (22) every link scales the unstable cdist by exactly beta_A = 5.27451 (p1) and the stable cdist down by beta_B = 1.55114 (p3's per-step factor): each side moves on its own saddle, which `IterateInference._bracketed_image` already does (U by the unstable curve's stretch, S by the stable manifold of `manifold_b_key`), as do `Trellis._scaled_image_cdist`, `StrongPip` and `Pseudoneighbor`. `check_iterate_link_scales_by_beta` predicted BOTH sides with A's factor, and `check_area_along_iterate_links` asserted the product preserved, whereas on a heteroclinic link it scales by beta_A / beta_B (link 7 -> 215: ratio 3.4005 = 5.27451 / 1.55114). Both laws now read each side's own saddle (homoclinic factor exactly 1, so no law case moves; none has a heteroclinic link, NUM-21); `IterateInference.py` is untouched. CLAUDE.md's "Area Preservation" states the homoclinic case as such. Regression `test_heteroclinic_iterate_links` on the unblasted 13u / 16s crossings (12 heteroclinic links; fails on base on 4 -> 8; the blasted recipes fail on 7 -> 215 and 5 -> 9 and pass now with 190 and 274 links).

**Evidence:** WP10 plan F13; `alllaws.py "build_nested(p3_unstable_steps=13, p3_stable_steps=15)"`

### NUM-21: Escape-lobe images are never computed
**Status:** accepted (D9 category) | **Severity:** low | **Related:** NEST-1, D9 | **WP:** — (documented by WP10)

A blast never iterates a bridge outside its zone (by design: exterior bridges grow exponentially), and the heteroclinic crossings of the law case have no iterate links (the preimage lies past p1's stable trim, the image past p3's unstable cut). This is D9's category, an image not computed: the topology never invents it, and growing or blasting further registers it. No code.

**WP9b:** this is exactly what the non-crossing rule needs for a member-less connecting class (`bn`, DEC-24): a heteroclinic crossing h of N with `f^{±1}(h)` in N closes a connecting arc into a curve with a known image, `S(B)[y, z_B] ∪ W^u(B)[z_B, h] ∪ S(A)[h, x]`.

**WP9c:** no longer needed to decide connecting ties: the walls (C17) decide `bn` without a closed curve through it.

**Evidence:** WP7 / WP9 notes (NEST-1); WP10 plan F14
---

## Performance

### PERF-1: The blast-end same-stability diagnostic does not scale past ~10 blasts
**Status:** closed (WP10, C5) | **Severity:** high | **Related:** DEC-12 | **WP:** WP10

`blast_zone(diagnose=True)` is the default and materialises every u×u hit:
- k28: 104,760 hits (1.58× the blast time)
- k32: 170,192 (1.76×)
- nested_outer: 589,693 (1.96×)

It cut reachable depth by 3-4 blasts. The hits are near-parallel image arcs at about 1e-4. Options: count instead of storing, check only new image bridges, or default `diagnose=False` when blasting.

**Evidence:** `src/tanglepack/numerics/TangleWorkbench.py:1495-1525`; `src/tanglepack/numerics/diagnostics.py`; `stress/numerics/report/report.tex:232-234, 1532-1564, 1594-1603, 1986-1997`; `.refactor/knowledge/testing.md` (OPEN)

**WP10:** the blast-end pass only COUNTS the same-stability hits of the bridges the blast added (`count_same_stability_crossings`, `report_resolution(new_bridges=)`, `ResolutionReport.same_stability_count`): one vectorised core shared with the full listing, testing only the candidate pairs holding a new segment and building no per-hit object, remembering nothing (every hit involves a curve new to the blast). Measured before (`perf1.py`, k28 zone at 1e-7): free to blast 7, then 1 405 hits at blast 8, 27 850 at 9, 189 341 at 10 (3.9 s against a 5.6 s blast), 826 877 at 11 (11.8 s against 9.9 s); every hit involved an image bridge and the old x old ones were the previous blast's, re-found. The full listing stays for on-demand use (`check_resolution`). Also (F6) the Tangle keeps one rtree per stability, so the 269 k of 270.5 k exact tests on `k32_series` that were same-stability pairs are never formed (`k32_series` 18.0 -> 14.5 s, `nested_outer_series` 19.3 -> 15.5 s, outputs identical).

### PERF-2: Blast cost grows geometrically, about 1.6-1.7× per blast in period-1 zones
**Status:** open | **Severity:** medium | **Related:** BCG-1, TEST-7 | **WP:** WP12

k28: ×1.61 crossings, ×1.65 points, ×1.71 time; the last blast takes 82 s and 9.0 GB, and memory stops it at 15 blasts. Coarser cutoffs hit the 780 s budget at 16-17. Blasting iterates every interior bridge, including those outside the fundamental domain (the same bridges that inflate the class table).

**Evidence:** `stress/numerics/report/report.tex:207, 1307-1364`; `stress/numerics/report/numbers.tex:458-466, 566-574, 776-784`

### PERF-3: Growth cost is dominated by the Python object graph (one Point per vertex); GPU gives no gain
**Status:** deferred | **Severity:** low | **Related:** — | **WP:** deferred

The map is at most 0.84% of growth time; refinement plus the linked list take 75-83%, at about 1 KB per point. GPU runs 0.96-1.00×.

**Evidence:** `stress/numerics/report/report.tex:103-108, 241-243, 523-565`

### PERF-4: The hole-propagation `_containing_bridge` is a linear scan over all bridges per backward step; termination uses `id(bridge)`
**Status:** closed (WP4: a bisect over the branch's crossing order, the span moved by `image_cdist` from the start bridge's crossing ids, termination on `BridgeId`; same containers as the scan everywhere) | **Severity:** low | **Related:** INT-2 | **WP:** WP4

**Evidence:** `src/tanglepack/topology/StablePartition.py:571-588, 1771-1800`

### PERF-5: Harness minors
**Status:** open | **Severity:** low | **Related:** PERF-1 | **WP:** WP12

- The E7 harness runs the same-stability pass twice.
- E11's cheap-path "unresolved" count is registry-wide.

**Evidence:** memory `stress-recs-plan-2026-10-07` (harness minors)

---

## Tests & coverage gaps

### TEST-1: No law case combines inversion with period-k, no nested case blasts its inner zone, and there is no heteroclinic or det<0 case
**Status:** partly closed (WP1) | **Severity:** medium | **Related:** INV-4, INV-11, NEST-1, NEST-3 | **WP:** WP1

**WP1:** the k=4 inversion cases and `nested_inner_blasts` (2 outer + 2 inner blasts) are law cases. Inversion × period-k waits for WP3: the period-2 orbits with both eigenvalues of DM² negative (k ≥ 4.1) cross only at their 8 anchors at 3-5 steps, so the partition raises INV-3, and they escape at about 6 steps for k=4.5. Heteroclinic is WP9 (NEST-1); det<0 is out of scope by C3.

**WP3:** inversion x period-k is the law case `inversion_period2`: the k=4.1 period-2 orbit `(-2.0488, -0.0488) <-> (0.0488, 2.0488)` with both DM^2 eigenvalues negative (-1.863 / -0.537), `k_value = 4`, `area_cutoff = 1e-5`, 9/9 steps, no zone (8 crossings at 3-5 steps grew into a real tangle once INV-1/INV-3 were fixed). Every law passes; not applicable: `iterated_cut_provenance`, `refined_children_inherit_word`, `member_matching_consistent`. On main it fails seven laws.

`LAW_CASES` are k10, k28_one_blast, k28_two_blasts, p3 (13 steps, unblasted), nested (2 outer, 0 inner) and inversion (no holes). Inner blasting goes untested even though E11 nested_inner grows 13→55. There is no inversion or det<0 case in E11 or in the deep sweep.

**Evidence:** `tests/cases.py:385-388, 493-501, 745-752`

### TEST-2: No law asserts structural stability under blasting
**Status:** closed (WP1) | **Severity:** medium | **Related:** BCG-1, BCG-4, BCG-7 | **WP:** WP1

**WP1:** `tests/regression/test_class_count_under_blasting.py` runs the D1 loop on `cases.SERIES` (k28 8 blasts, p3 9 blasts, inversion 2 blasts): the TOTAL class count never decreases, is constant while the hole set is constant, and every class has a member in the minimal trellis. Their strict xfails name WP2 (inversion) and WP6 (BCG-1/2). The count is total, not per kind: per kind it moves legitimately (k28 blast 0→1 turns two active classes inert once their image loops are registered; reference switches merge inert classes). **WP6:** all three laws pass on every series (the count is the alphabet's: k28 3,3,4,4,4,4,4,4,15; p3 9 x6,23,24,25,25; inversion 1,9,9).

Nothing checks that "more blasts refine, never contradict": a stable hole set between switches, class count independent of blast count between switches, or letters kept. The p3 burst and the linear growth pass a green suite (833 pass). `laws.py` reads the class count but never bounds it. Under the no-golden policy and the "topology change is expected" rule, this has to be a rule-based law.

**Evidence:** `tests/helpers/laws.py:1974`; `tests/regression/test_blast_monotonicity.py`; `tests/conftest.py:61-66`; memory `topology-change-is-expected`

### TEST-3: NOT_APPLICABLE pairs are dropped inside the layer tests, so pytest reports 0 skips
**Status:** closed (WP1) | **Severity:** low | **Related:** INV-4 | **WP:** WP1

**WP1:** every not-applicable pair has its own set-aside test, which runs the law, insists it still checks nothing and reports a visible skip (19 today); known issues share the same test as strict xfails.

The 15+ laws switched off on inversion, and the k10/p3/k28 exclusions, never appear in the summary. Last run: 833 passed, 3 deselected, 11 xfailed, 0 skipped, 100.88 s.

**Evidence:** `tests/helpers/law_tier.py:186-202`; `tests/cases.py:849-889`

### TEST-4: Provisional, untested rules
**Status:** provisional (WP1: no decision; anchor faces answered by WP3, see API-5) | **Severity:** low | **Related:** API-5 | **WP:** WP1

- The refined-point cdist is the mean of its neighbours.
- Anchor faces never close (the `grow_until_faces_closed` cap).

**Evidence:** `tests/conftest.py:68-83`; `docs/regions_bridge_identity_and_cleanup_plan.md:250`

### TEST-5: Test-policy open items
**Status:** closed (WP1) | **Severity:** low | **Related:** NUM-8 | **WP:** WP1

**WP9:** (c) closed: a third `TANGLE_COLOR_FAMILIES` family (violet and plum, colormap `Purples`), disjoint from the other two and from `HETEROCLINIC_COLORS` and at least CIE76 15 from every colour of theirs (the review rejected a first pick whose pinks sat 7-8 from `#c2185b` and `#cab2d6`; `test_families_do_not_look_alike` now holds every pair of families to 13, the old closest pair being 13.4); the strict xfail is gone.

**WP1:** (a) the `generation_bump` cache cell runs on `build_k28(blasts=1)` plus one blast and asserts the direct view changed; (b) the seed-seam stand-in is accepted (the real k=2.1 run never reaches near-ULP) and the private `_man_machine` poke is gone; (c) is a strict xfail `three_tangles_take_three_families` for WP9; (d) was already done in `e963893`.

- `generation_bump` cache cells are blind on k10.
- The p3 near-ULP regression uses a seed-seam stand-in.
- A third tangle has no `TANGLE_COLOR_FAMILIES` base.
- The colour-family-per-tangle rule test is still to do.

**Evidence:** memory `test-suite-policy-2026-10-05:52`

### TEST-6: `coverage_guard.py` trips when lines shift
**Status:** closed (WP11) | **Severity:** low | **Related:** — | **WP:** WP11

**WP11:** `coverage_guard.py` keys every line by (enclosing qualname, stripped text, occurrence inside the qualname) and reads the base source from git (`--base-ref`, default HEAD), so a shifted line is neither a false miss nor a mask for a real one (`tests/_tools/test_coverage_guard.py`). The dead-API whitelist is gone: dead API is deleted instead.

**Evidence:** `tests/_tools/coverage_guard.py`; memory `stress-recs-plan-2026-10-07`

### TEST-7: `blast_zone(max_points=...)` (WP5) is never exercised by the harness
**Status:** closed (WP10) | **Severity:** low | **Related:** DEC-13, PERF-2 | **WP:** WP10

**Evidence:** `stress/numerics/report/report.tex:1970`; `.refactor/knowledge/testing.md`

**WP10:** `test_the_budget_is_per_call` pins C6: a budget an earlier call spent in full does not stop the next call's iteration.

---

## API traps

### API-1: The pip choice is lost on every generation bump, so `_repin` after each blast is dead work
**Status:** closed (WP8) | **Severity:** medium | **Related:** SIMP-8, NUM-6 | **WP:** WP8

A rebuilt trellis starts with empty strong-pip slots. The session should carry the choice, as an Intersection, across generations.

**WP8:** the session remembers each fixed point's pip as an `Intersection` (`set_strong_pip`, `resonance_zone` — a zone's pip becomes the chosen one —, and the first default choice) and restores it on every rebuilt single-fixed-point trellis through `registry.find`; a pip that is no longer a candidate is left unset with a WARNING and `topology()` raises. `henon_cases._repin` / `_partition` and the test builders' `_repin` are gone; the stress harness's `repin()` still runs and now returns the restored pips (WP12 rerun note).

**Evidence:** `src/tanglepack/examples/henon_cases.py:142-153, 192, 199, 260`; `.refactor/knowledge/loom.md`; `.refactor/decisions.md` 2026-10-02; memory `strong-pip-choice-gotchas`

### API-2: The session gathers partitions from whichever per-fixed-point trellises happen to be cached (order-dependent)
**Status:** closed (WP8) | **Severity:** medium | **Related:** SIMP-8, NEST-1 | **WP:** WP8

`_gathered_partitions` scans every non-stale cached trellis without filtering by the requested fixed points ("first built wins"). The holes and pips gatherers filter only single-fixed-point trellises. Bridges come from the all-fixed-point trellis, so `trellis()` and `trellis(fp)` give different results.

**WP8:** a `Topology` reads partitions, holes and pips from the selected fixed points only: `homotopy_partition(fp3)` on the nested session holds fp3's branches alone, whatever was asked before. A heteroclinic bridge's class now needs both its fixed points selected (`bridge_classes` raises on an end no selected partition covers); WP9 builds the first such fixture.

**Evidence:** `src/tanglepack/loom/TangleSession.py:479-505, 528-545, 755-779`

### API-3: `PointBudgetExceeded.blast_context` is set dynamically and not declared; `MaxIterationsReached` and `TURNAROUND_LOG_STRETCH` are not exported
**Status:** closed (WP11) | **Severity:** low | **Related:** — | **WP:** WP11

**WP11:** `PointBudgetExceeded.blast_context` is declared (None outside a blast); `MaxIterationsReached` and `TURNAROUND_LOG_STRETCH` are exported from `tanglepack.numerics` and `tanglepack`.

**Evidence:** `src/tanglepack/loom/Blast.py:262, 388`; `src/tanglepack/numerics/budget.py` (no `blast_context`); `src/tanglepack/__init__.py:30, 78`; `src/tanglepack/numerics/__init__.py:27, 56`

### API-4: The deprecated no-op `invalidate_trellises` is still present
**Status:** closed (WP11) | **Severity:** low | **Related:** SIMP-12 | **WP:** WP11

**WP11:** deleted (no caller anywhere).

**Evidence:** `src/tanglepack/loom/TangleSession.py:1290-1306`

### API-5: `grow_until_faces_closed` never converges (it returns 0 immediately or hits the cap); anchors never close under the four-anchor model
**Status:** closed (WP3) | **Severity:** low | **Related:** INV-1, TEST-4 | **WP:** WP3

**WP3:** the four-anchor cause is gone: an inversion anchor is one degree-4 node with no virtual stub. What remains is a fact, not a bug: a crossing on the outer envelope never closes (the inversion anchor at 1 blast lies on the open outer face), and a non-inversion anchor seeded on one eigendirection keeps its `u-`/`s-` stubs. The Note is rewritten; the check-with-a-cap semantics are public API (WP11).

**Evidence:** `src/tanglepack/numerics/TangleWorkbench.py:1204-1265`

### API-6: Coordinate shape is (2,1) in FixedPoint and (2,) in BasePoint, so callers ravel everywhere
**Status:** closed (WP11) | **Severity:** low | **Related:** — | **WP:** WP11

Also: `partial_jacobians` is unexplained, and BranchPoint has a `_description_` placeholder.

**WP11:** the solver stores `(2,)` orbit points and eigenvectors and float eigenvalues, the FixedPoint placeholders match, and FixedPoint, ManifoldInitializer and TangleSession stop ravelling. `partial_jacobians` is documented (it propagates the eigenvector orientation and gives the map's orientation, C16); `BranchPoint`'s placeholder is filled in.

**Evidence:** `src/tanglepack/numerics/FixedPoint.py:10-21`; `src/tanglepack/numerics/BasePoint.py:9-11`; `src/tanglepack/numerics/BranchPoint.py:8-12, 41`; `src/tanglepack/numerics/Point.py:8-10`

### API-7: `build_nested` always returns the slug `nested_p1_p3`, so figures overwrite each other whatever the blast counts
**Status:** closed (WP11) | **Severity:** low | **Related:** DOC-4 | **WP:** WP11

**WP11:** `build_nested`'s slug adds `_outer<n>` / `_inner<m>` for non-default blast counts (the defaults keep `nested_p1_p3` / `heteroclinic_p1_p3`, so committed folders keep their names); `build_heteroclinic` derives its slug from it.

**Evidence:** `src/tanglepack/examples/henon_cases.py:282`

---

## Docs drift / stale notes

### DOC-1: Knowledge files and memory hold outdated suite counts and line numbers
**Status:** open | **Severity:** low | **Related:** — | **WP:** WP12

- `testing.md` says 782/795; the current count is 833.
- `loom.md`/`topology.md` line numbers come from 3b33baa.
- The `pytest-preexisting-failures` memory says 738.
- The `numerics-test-suite` memory cites an xfail `test_cdist_collision_growth` that no longer exists.
- The audit-polish `henon_basic_manifolds` decision is resolved (the script was deleted).

**Evidence:** `.refactor/knowledge/testing.md:3-4`; `.refactor/knowledge/loom.md:3`; memory `numerics-test-suite:127`, `pytest-preexisting-failures`, `audit-polish-2026-07:34`

### DOC-2: The report's framing contradicts the author's rule
**Status:** open | **Severity:** low | **Related:** BCG-1 | **WP:** WP12

- "Unreliable" legends and macros remain.
- "Fixed number of new classes per blast" is left unexplained.
- The burst cause is marked "not checked".

**Evidence:** `stress/numerics/make_numbers.py:428`; `stress/numerics/make_figures.py:816, 888`; `stress/numerics/report/report.tex:1284-1293`; `.refactor/knowledge/testing.md` ("Report framing smells")

### DOC-3: Stale Dev Notes
**Status:** possibly stale — verify | **Severity:** low | **Related:** SIMP-5, PK-7 | **WP:** WP12

- MinimalTrellis says the inertness definition is "changing" (settled 2026-09-16). **Fixed by WP6:** the note now says activity is read from the evidence table, never the alphabet.
- TrellisBranch's `return_period` question is resolved.
- ResonanceZone says "no inversion example" (`build_inversion` exists).
- PartitionFamily cites a "concurrent edit".
- TopologyResults asks whether a Hole should be a coordinate.
- ManifoldMachine has an old TODO.
- StrongPip/Pseudoneighbor say the key is "legitimately None".
- The topology-layer and loom memories say inversion is "not validated (no example)".

**WP3:** the ResonanceZone "no inversion example", the StrongPip / Pseudoneighbor / StablePartition / Arrangement / Trellis / Tangle inversion and anchor notes and the ManifoldInitializer `branch_index` question are rewritten. The rest stays with WP12.

**Evidence:** `src/tanglepack/topology/MinimalTrellis.py:36-40`; `src/tanglepack/topology/TrellisBranch.py:21-25`; `src/tanglepack/loom/ResonanceZone.py:37-39`; `src/tanglepack/topology/PartitionFamily.py:60-66`; `src/tanglepack/topology/TopologyResults.py:22-24`; memory `topology-layer:66`, `loom-and-subpackages:64`

### DOC-4: Stale artifacts
**Status:** possibly stale — verify | **Severity:** low | **Related:** API-7 | **WP:** WP12

- The committed `figures/henon_symbolic_itineraries/nested_p1_p3/partitions.png` has been stale since before 8148d48.
- The notebooks carry their own Hénon copies and now refine at 1e-7.

**Evidence:** `.refactor/knowledge/testing.md` (Deferred/open)

### DOC-5: CLAUDE.md says p3 at 13 steps is CLOSED and "blasting changes nothing", while p3 bursts at blast 6
**Status:** open | **Severity:** low | **Related:** PK-1, BCG-4 | **WP:** WP12

The same paragraph also notes that runs at 6+ blasts are unreliable, so the "blasting changes nothing" claim holds only up to 5 blasts.

**Evidence:** `CLAUDE.md` ("Period 3 at 13 unstable steps is CLOSED"); E11 p3 series

---

## Decision questions raised by the audit

Every question of the audit is now answered or deferred; see "Author decisions" below. DEC-19 was raised later by WP5 and answered on 2026-10-08 (D9). DEC-20 and DEC-21 were raised by WP7 and DEC-22 by WP8; the author answered all three on 2026-10-08 (D10-D12), implemented in WP8b. DEC-23 was raised by WP9 and answered on 2026-10-09 (D13), implemented in WP9b. DEC-24 was raised by WP9b and answered on 2026-10-09 (D14), implemented in WP9c.

| ID | Question | Blocks | Resolved by |
|---|---|---|---|
| DEC-1 | **Hole policy under growth and blasting:** (a) holes and partition fixed on the pre-blast fundamental trellis, with blasts only revealing words; (b) the current recompute, accepting refinement with depth; (c) holes on the true backward-image sub-arc only. Should old pairs be refined rather than replaced? May a crossing with neither a −1 nor a +1 link (a tip or fold crossing) change the topology? Should the ≥k exemption follow the physical orbit? | BCG-4, PK-1, INT-2, BCG-5 | D1 |
| DEC-2 | **Alphabet:** only the minimal trellis's classes, with blast children as evidence only? Quotient classes by the map or orbit shift, and collapse transport chains? Fold orphan loops into the bridge whose image chain contains them, or one loop bucket per tangle? | BCG-1, BCG-2, BCG-6 | D2 |
| DEC-3 | **Word definition:** the homotopy class of the actual registered image (drop the shortest-walk search; resolve by iterating one representative), or keep the dual-graph walk as the definition? | INT-1, INT-8, PK-1 | D3 |
| DEC-4 | **Inversion model:** one anchor per point (a degree-4 node) or one per unstable branch? Glue the elements on the same global side of p? Global "same side" (or demote it to INFO)? Zone = the union of the two loops at p, or each loop its own zone? Is the k=4 other-root saddle acceptable as the fixture? | INV-1, INV-4, INV-5, INV-6 | D4 |
| DEC-5 | **Anchor-only branch:** a closed singleton at the anchor, or one element `[anchor, end]` (WP4-B)? | INV-3 | C1 |
| DEC-6 | **Period-k references:** "a direct hole for every pair in any window F_j" instead of the provisional +1..+(k−1) table rule? Unify the window on all k branches? | PK-2, PK-3, PK-5 | D5 |
| DEC-7 | **Heteroclinic lobes:** pseudoneighbors and holes in B's trellis, A's, or neither? May a walk pass through both fixed points' fundamental segments? | NEST-1, NEST-2 | D7 |
| DEC-8 | **Hole openings and rows by the sign rule**, with geometry (`measured_side`) kept only as an optional check? | INT-3, NUM-1, NUM-2 | D6 |
| DEC-9 | **Element identity by bounding holes `(origin, iterate, side)`** instead of index? | BCG-7 | C2 |
| DEC-10 | **det J < 0 (b=−1) in scope?** If not, may the orientation-flip code go? | INV-11, SIMP-1 | C3 |
| DEC-11 | **Blast `in_zone` by innermost-zone classification**, so an outer blast never iterates inner bridges? | NEST-4 | C4 |
| DEC-12 | **Blast-end same-stability pass:** count only, new image bridges only, or `diagnose=False` by default? | PERF-1 | C5 |
| DEC-13 | **Confirm `blast_zone(max_points)` is a per-call budget.** | TEST-7 | C6 |
| DEC-14 | **Which blasting cutoff defines the symbolic result?** | NUM-5 | deferred |
| DEC-15 | **Crossings below the resolution ratio:** accept as undecidable, or invest in validated methods (interval/multiprecision, parameterisation-method manifolds)? | NUM-4 | deferred |
| DEC-16 | **Dissipative maps (b ≠ ±1) in scope?** Today they run unvalidated with one \|det J\| WARNING, and the area laws silently do not hold. | NUM-15 | deferred |
| DEC-17 | **`image_of` misses:** expected (unregistered image crossings), or should union-valued images be supported? | NUM-17 | deferred |
| DEC-18 | **Demote the bridge-row assertion to a resolution flag when the signs alternate?** Keep the assertion when the signs are equal. | NUM-1 | D6 |
| DEC-20 | **Where the walk may cross S** (raised by WP7): the image of a bridge crosses the stable manifold only at the images of its ends and in `S \ f(S)`, but the dual graph unifies the strong pip's window. They agree when S is trimmed at the pip orbit (every law case); inversion x period 2 at 10 steps runs past it and a walk is unreachable. (a) unify `S \ f(S)` per branch (pip-free); (b) cut S at the pip orbit in the partition and dual graph (the chalkboard's stage 3); (c) keep the 2026-09-14 design and warn (WP7 slice 1). (a) and (b) are no-ops on every law case. | PK-5 | D10 (implemented WP8b) |
| DEC-21 | **A walk pair no bridge spans** (raised by WP7): past the reference switch the word of an image-bridge class names the class of ITS image, a piece of `f(B)` no computed bridge spans (deep p3 runs, nested with 6 inner blasts, k=2.8 at 8 blasts). Today that class is unresolved with a WARNING and the dynamics is not reliable (D3: a construction error). (a) The minimal trellis also keeps the images of active image-bridge classes, iterating them when no bridge exists (map calls in loom; the minimal trellis grows by layers); (b) keep it a warned construction error and grow or blast until a bridge spans the pair; (c) letter the pair from the walk alone (the rejected first WP7 closure, against D2). | PK-1 | D11 (implemented WP8b) |
| DEC-22 | **Class-count decreases at depth: noise or definition?** (raised by WP8, reworded after its review): the remaining decreases are reference losses to a witness within a hair of a window crossing on the stable branch, unlinked both ways. Only k32 blast 10 is clearly definitional (a lambda-lemma copy 1.25e-2 inside the pair, error bar 2e-6). The others sit INSIDE the placement error bar `min(l_u, l_s) / ratio` (the ratio measures the crossing angle, not cdist placement): nested_outer blast 11 (1.0e-6 and 2.8e-6 from 38, bar 3.4e-6; the outer tangle is left with no reference; 4e-6 cross-frame order inconsistency measured on the same pair at blast 10), k32 blast 11 (1.55e-5, bars up to 2.1e-5), k32 blast 13 (1.1e-6 and 2.2e-6, bars 2.3e-6..4.5e-6), so noise is NOT ruled out, and under D1 a decrease is a bug unless root-caused. (a) A witness closer to a pair end than the end's placement error bar does not kill the reference (a tolerance, read off the registry; keeps the noise cases); (b) accept: a decrease is legitimate exactly when every lost reference is witnessed by an order-consistent witness (the law `series_lost_references_are_witnessed` checks this); (c) keep monotonicity by persistence: a reference stays one while its only witnesses are deeper (in unstable cdist) than the trellis it was found in -- phrased as "pairs of the truncation `{c_u <= U}`" to stay history-free (needs a mathematical check); (d) balance the trellis: cut N at an unstable depth matched to the stable cut at the pip, so blasting stops re-shaping the window at depth. A real-map refinement of the witness (Newton on the two curves) would settle each case but calls the map, so it belongs in numerics. | BCG-5 | D12 (implemented WP8b) |
| DEC-23 | **Equal shortest walks around a heteroclinic lobe** (raised by WP9): with D7 the lobe's hole splits B's window element, and two shortest walks of EQUAL length pass the lobe on either side (crossing B's window at the edge above vs below it); they read different words, so three classes are ambiguous and `is_reliable` is False on `heteroclinic`. The faces on each side of the window are the same face (the lobe and the A bubble across it are bays of the stable manifold); the walks differ by a loop around the lobe's hole. The registered member images do not agree on a side (one class's image matches the first walk, another's the second; the closure class has none). D3 does not say how to break the tie: (a) by the registered member image, unreliable where none is registered; (b) keep it ambiguous until growth or blasting registers the images; (c) enrich the walk with the non-crossing constraint of the images of the A arcs bounding the lobe (f(b) cannot cross f(W^u_A)). WP9 slice 1 keeps (b) as a strict xfail (`every_active_class_walks` x `heteroclinic`). WP10: the two ambiguous classes of k28 at 8 blasts are this kind (two walks of equal length crossing the window at different edges; none of the six faces they visit holds a u x u hit point, so they are not resolution artifacts). | NEST-1 | D13 (implemented WP9b; residue DEC-24, closed WP9c) |
| DEC-24 | **A member-less connecting class around another tangle** (raised by WP9b): on `heteroclinic` the class `bn` = {A:R_1, B:L_(2.0;1)} (no member, D11; a connecting class) has two shortest walks that differ by a loop around B's peninsula on A's window. The non-crossing rule needs a closed curve with a known image through the class, and a connecting arc closes only through a heteroclinic crossing of N whose image or preimage is in N; none of the trellis's heteroclinic crossings has an iterate link (NUM-21). So the side is not determined by T(q0). (a) Undecidable by data, like D9: the class stays flagged and `every_active_class_walks` is re-scoped to exclude a member-less connecting class with no linked heteroclinic crossing (the xfail is removed); (b) keep the strict xfail until slice 2 (closing through a linked heteroclinic crossing) and a fixture with linked heteroclinic crossings exist (p3 at 17 unstable steps: eight linked crossings in N, 27 s per build, still one connecting tie); (c) slice 2 alone, with the fixture moved to p3 17u under `perf`. | NEST-1 | D14 (implemented WP9c) |
| DEC-19 | **Period-k stagger** (raised by WP5): a reference trajectory on the outermost lap of the longest unstable branch has no `+j` image in the trellis. (a) accept and warn (WP5 slice 1); (b) trim N to crossings whose `+1..+(k-1)` images are registered (an unstable-side cut, making every reference Z_k-complete; changes the deep p3 / inversion references, BCG-4); (c) virtual partition boundaries at the predicted landings (invents holes, rejected). Orphan window pairs are not pairs (NUM-18). | PK-3 | D9 (accept) |

## Author decisions (2026-10-07)

| # | Decision | Resolves |
|---|---|---|
| D1 | **Full recompute.** Topology runs only after the numerics; ANY growth or blast means the topology is recomputed from scratch. A reference-pair switch on the recomputed trellis is legitimate. Planned loop: blast → topology → blast → compare. The class count is non-decreasing: it plateaus and may jump again after much more growth. A DECREASE, or linear growth over many iterates, is a bug. | DEC-1 |
| D2 | **Alphabet = the minimal trellis's bridges only.** Blast children are image evidence only. A loop folds FORWARD into the bridge whose image chain contains it, never into a standalone `BridgeClass(x,x)`. No collapsing of transport chains. | DEC-2 |
| D3 | **The shortest dual-graph walk IS the word.** A walk must always exist: an unreachable walk or a virtual `new*` class is a construction bug to root out, not a state to model. f(bridge) (the trellis itinerary) makes the word only for singletons that cannot walk; elsewhere it is the cross-check. | DEC-3 |
| D4 | **Inversion.** An anchor point has two stable and two unstable branches, each with its own left and right side. L of stable branch 0.0 is the same region as R of 0.1. Zones partition the plane and must be built correctly at an inversion point. A cross-side flag applies only to a pair on the SAME branch, never to a bridge between two branches. | DEC-4 |
| D5 | **Period-k.** Keep the +1..+(k−1) iterate rule for the other branches' holes, but close the gap: a missing table link never silently leaves a branch without holes. | DEC-6 |
| D6 | **Signs decide** hole openings and bridge rows. Geometry stays as an optional cross-check that flags under-resolution; it never asserts when the signs alternate. | DEC-8, DEC-18 |
| D7 | **Heteroclinic.** A heteroclinic crossing pair on B's stable branch can be a pseudoneighbor pair and gets a hole in B's partition. | DEC-7 |
| D9 | **Period-k stagger accepted** (2026-10-08, DEC-19). A branch whose reference images were never computed simply has no window hole yet; the WARNING naming it stays, and the topology never invents a hole or trims the trellis to complete orbits. | DEC-19 |
| D10 | **Cut S at the pip orbit** (2026-10-08, DEC-20). The partition and the dual graph use the stable manifold cut at the strong pip's orbit, as the pseudoneighbor stage already does (chalkboard stage 3); the pip still decides the window. | DEC-20 |
| D11 | **Words need no computed bridge** (2026-10-08, DEC-21). Words come from the dual-graph walk (f(bridge) only for singletons), so a pair a word names that no computed bridge spans is still a legitimate class: it gets a normal letter, its own word is walked from its elements' landings, and it does NOT make the dynamics unreliable. Such classes are marked "no computed bridge yet" (informational only). Supersedes the second half of C9. | DEC-21 |
| D12 | **Balanced trellis** (2026-10-08, DEC-22). The crossing set is cut at an unstable depth matched to the stable cut at the pip, so blasting cannot reshape the fundamental segment at depth through lambda-lemma copies; the class count then never drops under blasting. | DEC-22 |
| D13 | **Non-crossing constraint breaks walk ties** (2026-10-09, DEC-23). When several shortest dual-graph walks of equal length read different words, the walk is enriched with the fact that f(bridge) cannot cross the images of the unstable arcs bounding a lobe (f(W^u) is a wall); that picks the side combinatorially. It stays a walk; no registered image is needed. | DEC-23 |
| D14 | **Every bridge is the same** (2026-10-09, DEC-24). A bridge connects two crossings; each lies in a stable partition element, and side matters for those elements. A bridge between elements of two different tangles is no different from one within a tangle: no special case for connecting classes. A class exists only where bridges exist, so "a connecting class with no computed bridge" means the construction that produced it is wrong. | DEC-24 |
| D8 | **Simplify freely** toward the chalkboard design with a strong architecture. Aligning the public API with how the author uses it is a later, separate plan. | — |

Taken under D8 (Claude's call, revisit if wrong):

| # | Decision | Resolves |
|---|---|---|
| C1 | An anchor-only stable branch is ONE element `[anchor, end]`. | DEC-5 |
| C2 | An element is identified by its bounding holes `(origin, iterate, side)` rather than its positional index, so recomputes can be compared. | DEC-9 |
| C3 | det J < 0 stays out of scope; the six orientation-flip sites collapse into one helper (kept, not deleted). | DEC-10 |
| C4 | A blast's `in_zone` classifies by innermost zone, so an outer blast never iterates an inner tangle's bridges. | DEC-11 |
| C5 | The blast-end same-stability pass checks only the new image bridges and counts hits instead of storing them. | DEC-12 |
| C6 | `blast_zone(max_points)` is a per-call budget. | DEC-13 |
| C7 | At an anchor a bridge of `u(o, b)` meets `s(o, b)` (same orbit point and branch index: equivariant under the map, and the pair carries the anchor's own sign). Elements are NOT glued across the periodic point: every branch keeps its own two sides (WP3). | DEC-4 |
| C8 | One resonance zone per fixed point: `ResonanceZone.key` is the fixed point; an inversion zone is the union of its two loops (WP3). | DEC-4 |
| C10 | An element is identified across recomputes by its BOUNDING CROSSINGS `(branch_key, side, lo_id, hi_id)` (the hole-bridge endpoints that cut it; registry ids persist in a session), which refines C2: the `(origin, iterate, side)` of its bounding holes changes at every reference switch even where the hole bridges stay. The alphabet letters by it; `ElementRef` and element names stay positional within one recompute (WP8). | DEC-9 |
| C17 | (under D13 / D14, WP9c) Walls: a class's arcs lie on its members' unstable branches, or one map step past its namer's (D11); the decided walk (`unique` / `trivial`) of a class walked EARLIER whose branches are disjoint from this class's is a wall, and a tie loses the walks whose chords interleave with a wall's on a face ring. A tie filter only; when every walk crosses a wall the tie stands. | DEC-24 |
| C18 | (under D14, WP9c) D14 is tested as "a class connects two tangles only where a computed bridge in the trellis connects them", not "every connecting class has a computed member", which contradicts D11 (`bl`, `bn` are images of computed connecting bridges, member-less by D12's cut, computed at p3 17u). | DEC-24 |
| C19 | (under D8, PK-4, WP10b) Chords of a lobe image pair per stable branch: chain order on the base branch (the author's rule), STABLE order on any other (f(H) holds no anchor, so a branch walked from its anchor enters and leaves it in turn); the base is read at a non-anchor end. Only refinement labels change. | PK-4 |
| C20 | (under D1, NUM-12, WP10b) A blast's image curves are computed manifold: kept by the workbench, re-indexed query-only by every recompute until growth merges them into a branch, cut like any curve (crossed fewer than twice: one partial piece). Grow, blast and recompute commute. | NUM-12 |
| C9 | The alphabet (the minimal trellis's classes, D2) is CLOSED under the word map over REGISTERED bridges: a class a word names that some bridge spans (the evidence table) is added, walked, lettered by the session and records `named_by`. A pair no bridge spans is never a class: its namer is unresolved with a WARNING and `is_reliable` is False (D3); whether it must get its bridge is DEC-21 (WP7). (Second half superseded by D11, WP8b: a pair no computed bridge spans is a flagged class.) | D2, D3 |

Deferred: DEC-14 (which blasting cutoff defines the result; re-measured in WP12), DEC-15 (validated arithmetic), DEC-16 (dissipative maps), DEC-17 (union-valued `image_of`).

## Work packages

Each WP runs on its own branch `wp/<n>-<slug>` through analyst → planner → test-engineer → coder → code-reviewer and merges into main in this order. Every WP starts with rule-based laws or strict xfails (no golden pins) and ends with `env/bin/python -m pytest -q -rxX` green and `KNOWN_ISSUES` no longer than before.

| WP | Scope | Issues |
|---|---|---|
| WP0 | This catalog, `pipeline_chalkboard.md`, decisions recorded | — |
| WP1 | Laws encoding the decisions (walk existence, same-branch cross-side, class monotonicity over a blast series, alphabet ⊆ minimal trellis), the k=4 inversion fixture, NOT_APPLICABLE as skips | INV-4, TEST-1..5 |
| WP2 | Exact iterate bookkeeping: inference after every blast, map-gated, relative `cdist_tol`, `image_cdist` the one scaling path (fixes the λ_u lap factor), dead `None`-key branches | INV-2, PK-7, BCG-3, INT-6, INT-7, NUM-8, SIMP-3 |
| WP3 | Inversion end to end: one anchor point with four branches, anchor-only branch, zone ring by full key, same-branch cross-side, one branch-ordering rule | INV-1, 3, 5..10, API-5 |
| WP4 | Sign rule for hole openings and rows, geometry as a flag, bisect `_containing_bridge` | INT-2..5, NUM-1, NUM-2, PERF-4 |
| WP5 | Period-k hole gap closed, explicit flags | PK-2, 3, 6, 8, 9 |
| WP6 | Alphabet from the minimal trellis, forward loop folding, one letter generator | BCG-1, 2, 6, SIMP-6 |
| WP7 | The walk always exists: root-cause unreachable / `new*` / "abuts no hole" / cross-branch skips; delete virtual classes | PK-1, 4, 5, INT-1, 8, 9, NEST-3 |
| WP8 | Recompute semantics: verify reference switches, root-cause every decrease, element identity by holes, one `Topology` result per generation, pip choice carried by the session | BCG-4, 5, 7, API-1, 2, SIMP-8 |
| WP8b | The author's D10-D12: the trellis of the pip (S cut at the pip orbit, N cut at the deepest registered pip-orbit point), words need no computed bridge | PK-1, PK-5, BCG-5 |
| WP9 | Nested and heteroclinic: A×B fixture, heteroclinic holes in B, innermost-zone blasting, per-fixed-point `rebuild_bridges` | NEST-1, 2, 4..6 |
| WP10 | Numerics and performance leftovers | PERF-1, NUM-9..16, INT-10, TEST-7 |
| WP10b | Leftovers | NUM-12, NUM-19, NUM-20, PK-4, WP11 review nits |
| WP11 | Simplification sweep and API traps | INV-11, SIMP-*, API-3, 4, 6, 7, TEST-6 |
| WP12 | Rerun E7/E11 and the figures, confirm class counts plateau, update CLAUDE.md / knowledge / memory, close issues | NUM-5..7, PERF-2, PERF-5, DOC-1..5 |
| deferred | — | NUM-3, NUM-4, NUM-17, NUM-18, PERF-3 |

### Laws added by WP1

Each strict xfail below is a `cases.KNOWN_ISSUES` entry naming the WP that turns it green; an XPASS fails the run, so that WP removes the entry.

| Law | Cases that xfail | WP |
|---|---|---|
| `cross_side_only_within_a_branch` (D4; replaces `itinerary_pairs_same_side`) | none (WP3) | — |
| `every_active_class_walks` (D3) | none (not applicable on p3 and k28_one_blast: singleton words) | — |
| `no_virtual_class` (D3; replaced by `alphabet_closed_under_words`, WP7) | none (WP3) | — |
| `alphabet_inside_minimal_trellis` (D2; an iff over the evidence table since WP6) | none (WP6) | — |
| `series_class_count_never_decreases` (D1) | none (WP2) | — |
| `series_class_count_constant_while_holes_constant` (D1, D2) | none (WP6) | — |
| `series_alphabet_inside_minimal_trellis` (D2) | none (WP6) | — |
| `three_tangles_take_three_families` (TEST-5(c)) | none (WP9) | — |
| `trellis_is_cut_at_the_pip_orbit` (D10, D12; WP8b) | none (fails on base on the three inversion cases) | — |
| `series_class_count_never_decreases` on `k32_series`, `nested_outer_series`, `inversion_series` to 4 blasts (D12; WP8b) | none (fails on base on all three) | — |
| `series_closed_alphabet_constant_while_holes_constant` (D11; WP8b) | none (fails on base on `p3_series` and `k32_series`) | — |
| `trellis_crossings_belong_to_its_selection` (WP9) | none (fails on base on `heteroclinic`) | — |
| `tangles_in_case_order` (NEST-6; WP9) | none (fails on base on `nested`, `nested_inner_blasts`, `heteroclinic`) | — |
| `case_heteroclinic` (WP9) | none | — |
| `every_active_class_walks` on the new case `heteroclinic` | heteroclinic | WP9 (DEC-23) |
| `partition_covers_branch` on the new case `heteroclinic` | none (WP10: the law's last clause checks each element end against the halves the holes open, derived from their bounds, bridge sides and crossing signs) | — |
| `walks_agree_with_registered_images` (D3 + D13; WP9b) | none (fails on base on `heteroclinic`: `c`, `bh`) | — |
| `member_lobe_classes_have_one_walk` (D13; WP9b: a walked class with a member lobe has a unique or trivial walk) | none (fails on base on `heteroclinic`: `a`, `c`) | — |
| `walk_starts_on_the_landing` (D3, C14; WP9b) | none (a no-op on the law cases; fails on base on `k28_deep_blast_sliver`) | — |
| `lobe_enclosure_is_geometry` (D13; WP9b: the flood against point-in-polygon) | none | — |
| `registered_images_keep_markers` (D13's physics; WP9b) | none | — |
| `series_no_walk_is_ambiguous` (D13; WP9b) | none (fails on base on `k28_series` and `k32_series`) | — |
| `every_active_class_walks` on `heteroclinic` | none (WP9c: XPASS, removed; the walls decide `bn`) | — |
| `walks_avoid_disjoint_images` (D13, D14; WP9c: decided walks of classes on disjoint unstable branches never interleave in a face) | none (NOT_APPLICABLE on the seven law cases with no such chord pair: k10, k28 x2, inversion x3, inversion_period2) | — |
| `alphabet_closed_under_words`, D14 clause (WP9c: a class connects two tangles only where a computed bridge does, C18) | none | — |

The new law cases also put existing laws in KNOWN_ISSUES: on the inversion cases `one_anchor_per_unstable_branch`, `arrangement_euler`, `anchor_bridge_class_leads_its_tangle`, `partition_singletons` and `stable_nodes_never_own_a_singleton` (WP3, INV-1), `arrangement_regions_disjoint` (WP10, NUM-16) and, on `inversion_two_blasts`, the 15 hole-dependent laws that are vacuous with 0 holes (WP2, INV-2); on `nested_inner_blasts` `arrangement_regions_sound` and `arrangement_regions_disjoint` (WP10, NUM-16). WP2 must add `stable_nodes_never_own_a_singleton` on `inversion_two_blasts` (WP3) once holes appear there.

**WP2** (KNOWN_ISSUES 52 -> 34): with the β^k lap the 15 vacuous hole laws and the two `inversion_series` entries XPASSed and were removed, as were `anchor_bridge_class_leads_its_tangle` and `cross_side_only_within_a_branch` on `inversion_two_blasts` (WP3 tags; both still xfail on `inversion`). It added `stable_nodes_never_own_a_singleton` x `inversion_two_blasts` (INV-1, WP3); `partition_singletons` x `inversion_two_blasts` is no longer vacuous and shares that entry. WP2 also added the crossing-layer laws `iterate_table_partial_bijection`, `iterate_links_pass_the_map_gate` and `iterate_table_complete` (pass on every law case).

**WP3** (KNOWN_ISSUES 34 -> 21): the inversion anchor is one point with four branches. The 15 WP3 entries were resolved: 13 pairs pass and were removed (`one_anchor_per_unstable_branch` x 3, replaced by the stronger `one_anchor_per_periodic_point`; `arrangement_euler` x `inversion`, `inversion_one_blast`; `no_virtual_class` x 3; `partition_singletons` and `stable_nodes_never_own_a_singleton` x `inversion_one_blast`, `inversion_two_blasts`; `cross_side_only_within_a_branch` x `inversion`), and two with other roots were re-tagged with evidence: `anchor_bridge_class_leads_its_tangle` x `inversion` -> WP5 (PK-3) and `arrangement_euler` x `inversion_two_blasts` -> WP7 (INT-9). New laws (pass on every law case): `one_anchor_per_periodic_point`, `anchor_rays_are_its_branches`, `anchor_opening_faces_its_hole` (not applicable where no hole touches an anchor: k10, k28_two_blasts, inversion). New law case `inversion_period2`. Hole opening records are 4-tuples `(id, which, row, branch_key)`.

**WP4** (KNOWN_ISSUES 21 -> 21, none was tagged WP4): signs decide hole openings and the row assertion (D6), the geometry is a flag (`geometry_flags`, `Trellis.punch_holes(diagnose=True)`), the backward chain is `image_cdist` from the start bridge's crossing ids with a bisect container, and the measured side reads the nearest segment. INT-2 (by design), INT-3, INT-4, INT-5, NUM-1, NUM-2 and PERF-4 are closed. The law `openings_follow_signs` replaces `openings_missing_only_at_anchor_or_tail`; new regression `test_coarse_blasts_rows_by_signs` (builder `p3_coarse_blasts`) and `test_deep_p3_every_class_walks`. Law cases are unchanged (holes, partitions, classes, words); the deep p3 runs lose their unreachable class; E11 p3 / nested_inner at 1e-5 run 20 blasts. The next blocker of the deep cells is `OuterFaceError` (INT-10, WP10).

**WP5** (KNOWN_ISSUES 21 -> 21, one reason re-tagged): D5's period-k rule is kept and stated (`_is_forward_beyond_fundamental(pair)`), and no branch is left without its window holes silently. A missing `+j` appearance is a crossing past the computed unstable branch (the stagger of a period-k chain); `extend_pseudoneighbor_trajectories` names each such stable branch in one WARNING with the `image_cdist` prediction, and nothing is invented. The strong pip's k cuts are placed by `image_cdist` (PK-6), the hole-chain step cap is gone (PK-9), and an inert class mapping onto inert classes is silent (PK-8). PK-2, PK-6, PK-8 and PK-9 are closed; PK-3 is slice 1 done, with DEC-19 open for slice 2. New law `window_holes_on_every_branch` (passes on every law case) and regression `test_period_k_window_gaps`. Measured gap WARNINGs per build: inversion 1, inversion_one_blast 1, inversion_period2 3, p3_15 / p3_6_blasts / p3_coarse_blasts 2 each; k10, k28 x2, p3, nested, nested_inner_blasts and inversion_two_blasts 0. Law cases are unchanged (holes with openings, partitions, classes, words). `("anchor_bridge_class_leads_its_tangle", "inversion")` stays a strict xfail tagged WP5: its hole cannot be made without inventing one (DEC-19).

**WP6** (KNOWN_ISSUES 21 -> 13): the alphabet is the minimal trellis's classes (D2). `bridge_classes` is the evidence table, the minimal trellis reads activity from it and keeps it (`MinimalTrellis.classes`), and `MinimalTrellis.alphabet` is what the session letters; registered evidence never adds a symbol. Loops fold FORWARD into their image-chain parent's class (one parent per pair, asserted); orphans are logged and unclassed. One letter generator, `topology.letter`, over disjoint digit sets. BCG-1, BCG-2, BCG-6 (decided D2) and SIMP-6 are closed. All eight WP6 entries XPASSed and were removed (`alphabet_inside_minimal_trellis` x 3 inversion cases, `series_class_count_constant_while_holes_constant` x k28_series / p3_series, `series_alphabet_inside_minimal_trellis` x 3 series). New law `loops_fold_forward`; the classing laws read the evidence table, `alphabet_inside_minimal_trellis` is an iff; two NOT_APPLICABLE entries on `inversion` (`every_active_class_walks`, `matrix_is_token_counts`: its alphabet is one inert class until DEC-19). Regression `test_blast_children_are_evidence_only`. Law cases change only on the inversion cases (4 -> 1, 10 -> 9, 10 -> 9 classes). Open for WP7: the alphabet is not closed under the word map in deep runs (PK-1).

**WP7** (KNOWN_ISSUES 13 -> 2): the walk exists wherever the stable manifold stops at the pip orbit. Two root causes explain every walk-side failure. (1) The alphabet was not closed under the word map: the virtual `new1` of the deep p3 runs is the class of the image of an image bridge, a piece of `f(B)` the minimal trellis does not keep (PK-1, NEST-3). `symbolic_dynamics` now closes it over the classes registered bridges span (C9), virtual classes are deleted, and a pair no bridge spans unresolves its namer with a WARNING (`is_reliable` False): the deep runs keep one or two such classes, DEC-21. (2) The dense arrangement drew arcs over manifold that was never computed (INT-9 -> NUM-16); `Trellis.is_computed_arc` decides, and a gap dangles. Evidence disagreements carry an `EvidenceDiagnosis` (INT-1), and a straddling landing is unreliable (INT-8). The dual graph warns where its unified span is not `S \ f(S)` (PK-5, DEC-20). Removed from KNOWN_ISSUES: `arrangement_euler` x `inversion_two_blasts`; `arrangement_regions_disjoint` x k28_one_blast / k28_two_blasts / inversion_one_blast / inversion_two_blasts / nested_inner_blasts and `arrangement_regions_sound` x nested_inner_blasts (NUM-16, root cause INT-9); `open_p3_deep_runs` x 3; `anchor_bridge_class_leads_its_tangle` x `inversion`, which became NOT_APPLICABLE. That law is re-scoped by D9 to tangles whose every stable branch carries a hole: on `inversion` s(0.1) has none, the anchor bridge of u(0.1) is a loop, and the anchor bridge of u(0.0) maps onto it. New laws: `alphabet_closed_under_words` (replaces `no_virtual_class`: every symbol a table class, every closure class spanned by a registered bridge, an unresolved class makes the dynamics unreliable) and `dual_unified_is_where_f_of_s_stops`; both pass on all 10 law cases. Law-case words, refined words, matrices, letters and reliability are identical to base. Left for later: PK-4 (slice 2, foreign-branch chords), DEC-20, escape-lobe images (WP9 / WP10), k28 at 8+ blasts (2 ambiguous closure walks over 173 u x u artifacts, WP10), and the deep-run residual laws (WP8 / WP10).

**WP8** (KNOWN_ISSUES 2 -> 2, no entry flips or is added): recompute semantics (D1) and stable identity. (1) INT-7: crossings are registered along each unstable branch in unstable-cdist order, so registry ids are reproducible within and between processes; every id-bearing output is renumbered (describe reports, figure labels, logs, the CLAUDE.md k=2.8 id facts, regenerated), nothing else moves. (2) BCG-5: an unlinked landing is read against the pair's registered images on its own polyline; the noise punctures of nested at 10 outer blasts (holes `[0, 12]` -> `[4, 12]`, classes 9 -> 13) and k32 at 12 blasts (34 -> 42 holes, 24 -> 29 classes) are gone; the decreases left (witnesses inside the placement error bar, noise not ruled out; only k32 blast 10 is clearly a lambda-lemma copy) are DEC-22. (3) BCG-4 verified: every switch is a definitional discovery with a table-only witness. (4) BCG-7: letters by element identity, the bounding crossings (C10). (5) SIMP-8 / API-1 / API-2: `TangleSession.topology()` is one `loom.Topology` per `(generation, selection, pips)` that recomputes holes and partitions from the pips; the session carries the pip as an `Intersection`; gathering follows the selection. Behaviour changes: products need no manual stage call (the default pip is taken and remembered) and a manual `punch_holes(...)` / `partition_stable_manifold(...)` on a per-fixed-point trellis never feeds them; `session.strong_pip(fp)` after a blast is the chosen pip (base: None); `set_strong_pip` to another candidate recomputes the holes; `homotopy_partition(fp3)` holds fp3's partitions only; `trellis(fp) is trellis([fp])`; the unpartitioned-branch WARNING is gone (a heteroclinic bridge needs both its fixed points selected). Law cases, `p3_15`, `p3_16`, `p3_6_blasts`: hole, reference and class counts, the id-free `describe()` of the dynamics (letters included) and `is_reliable` identical to base; series counts unchanged (k28 3,3,4,4,4,4,4,4,15; p3 9x6,23,24,25,25; inversion 1,9,9). New laws `series_lost_references_are_witnessed` and `series_letters_follow_element_identity` (the latter fails on base: 3 classes re-lettered on `p3_series` from blast 6); regressions `test_registry_ids_are_reproducible`, `test_deep_reference_survives_noise` (~11 s), facade `test_session_topology`, unit `test_unlinked_landing_is_read_against_the_pairs_images`. Ids quoted in this WP's notes (BCG-4, BCG-5) are in the base `083f910` numbering. Left: DEC-22; element NAMES by identity (positional by design); `partition_element_for` still scans every snapshot, `reindex_from`'s unused remap and `ResonanceZone._live_id` (WP11); heteroclinic selection (WP9); the E11 harness `repin()` now returns the restored pips (WP12); `scripts/refactor_proof/hole-side-and-own-blast/*` import the deleted `hc._repin` / `hc._partition` (historical, not maintained); the deep-run law `partition_covers_branch` still fails on `p3_6_blasts` / `p3_15` / `p3_16` (open elements that are no hole's stretch, 6 / 6 / 5 spans, where hole bridges of one chain nest on a branch -- the law assumes disjoint hole stretches; not a law case, not examined here). Review fixes (WP8): an end of the pair placed by `image_cdist` in r's frame gets the 2D endpoint-collision guard (an endpoint's own unlinked orbit is no puncture); `reference_witness` rejects a witness whose registered chain falls outside the pair in another frame (order is preserved along a branch); `classify_strong_pips()` never replaces a lost chosen pip by the default (one WARNING; `chosen_strong_pip(fp)` reads the carried pip); stages 3-6 and the zone key live in `Topology`; DEC-22 reworded (noise not ruled out, measured error bars).

**WP8b** (KNOWN_ISSUES 2 -> 2, no entry flips or is added; NOT_APPLICABLE 28 -> 21, the seven removed pairs all on `inversion`): the author's D10-D12. (1) D10 + D12, the trellis of the pip: choosing a strong pip cuts the trellis to `T(q0) = (W^u[P, f^m(q0)], W^s[P, q0])` in ONE place, `Trellis.cut_at_pip_orbit` (the `strong_pip` setter; None restores; classification un-chooses first): stable branch `j` at `c_s(f^j q0)`, the unstable branch of `f^i(q0)` at `c_u(f^i q0)` for the `k` deepest registered `i`. A crossing past a cut is a point of X: it leaves the branches, the arc across it joins its bridges (`Bridge.joined`), `iterate` / `image_cdist` stay raw and membership is `crossing in trellis`. Pseudoneighbors lost their private cut (`_strong_pip_cuts`, `_within_cut`); trajectories, image chains, singleton landings and region images stop at an image outside T; `DualGraph._check_where_f_of_s_stops` is gone (vacuous); `TangleSession.trellis` cuts a selection of several fixed points at their pips on every read (review fix: the cut was in `Topology`, so `session.trellis()` showed the computed N until a product ran: nested_outer_series after 10 outer blasts, 736 crossings / 525 bridges before `topology()`, 448 / 322 after), `Topology` resolves members through it, and `cartoon_enclosures` takes the plotted dual graph's trellis (a missing bridge asserts). PK-5 and BCG-5 closed. (2) D11: a walk pair no computed bridge spans is a member-less, lettered, walked class flagged "no computed bridge yet" (`BridgeClassEntry.spanned`); `_unspanned` deleted; PK-1 closed. Behaviour changes: law cases other than the inversion three are identical (id-free `describe()` signatures, `PYTHONHASHSEED=0`); `inversion` 1 -> 6 classes (holes 1 -> 8), `inversion_one_blast` 9 -> 7, `inversion_two_blasts` 9 -> 8, every law passing on them (seven NOT_APPLICABLE pairs apply now; `every_active_class_walks` x `inversion` is reworded: its active classes land on singletons); the deep runs close to 43 / 49 / 43 classes, reliable; `session.trellis(fp)` with a chosen pip holds `T(q0)` only (its `bridges` may hold joined arcs the workbench does not); the inversion series is 6, 7, 8, 9, 9 (base 1, 9, 9, 9, 6). New laws `trellis_is_cut_at_the_pip_orbit`, `series_closed_alphabet_constant_while_holes_constant` and `series_session_trellis_is_history_free` (plus the facade test `test_the_selection_trellis_is_cut_before_any_product`); new series `k32_series` (10 blasts, ~18 s) and `nested_outer_series` (11, ~20 s), `inversion_series` to 4 blasts; new regression `test_pip_orbit_cut`, unit `test_trellis_cut`; `test_unified_segment_precondition` deleted (its WARNING is gone). Left: k=2.8 at 8 blasts unreliable (ambiguous walks, evidence disagreement: WP10); inversion x period 2 at 10 steps is no law case (`bridge_endpoints_on_own_branch` on workbench bridge `(18,19)`, WP7-D1); E11 k32 to 13 blasts is a script, not a test (2 min, ~5 GB); re-classifying every candidate on every rebuild is O(N^2) (PERF, WP10); PK-4 untouched; for the author: D12 stops discovery at the registered depth of the pip's orbit (a smaller blast `min_separation` reaches deeper), and the closed alphabet can drop by one across a witnessed switch (k32 blast 8 -> 9, seed constant).

**WP9** (KNOWN_ISSUES 2 -> 3; NOT_APPLICABLE unchanged): nested and heteroclinic tangles. (1) NEST-1 / NEST-2 / D7: a crossing is in a trellis only when both its branches are selected (a tangle alone holds no heteroclinic crossing; its arcs are joined across them), the pair test reads the cycle, `k` and β of the ARC's fixed point (its two 2D endpoint collisions, scaled and own-frame, skip only a candidate on the pair's own stable manifold), and stages 4-6 run on the selection's trellis (`loom.Topology`; stage 3, the pip, stays per fixed point). The selection cut at both pips is `T(q0_A) ∪ T(q0_B)`; D10 / D12 needed no change. (2) NEST-4 / C4: a blast iterates the bridges whose innermost zone it is, read against the stored zones of the OTHER fixed points, so a zone the session does not store (`define_resonance_zone`, a replaced or copied handle) still blasts its own bridges (review: comparing `classify_bridge(bridge) is zone` blasted nothing there). (3) NEST-5: `clear_bridges` / `rebuild_bridges` per fixed point. (4) NEST-6: the recipes construct the outermost fixed point first. (5) TEST-5(c): a third colour family. New law case `heteroclinic` (`henon_cases.build_heteroclinic`: p3 16u / 10s, p1 11u, two outer blasts, 1e-7, 0.7 s): on base refs 3, holes 50, classes 31, closed 53, 0 connecting, reliable; on the branch refs 4 (one heteroclinic reference on p1's window), holes 77, classes 57, closed 81, 2 connecting, NOT reliable (DEC-23). New laws `trellis_crossings_belong_to_its_selection` (partition layer; fails on base on `heteroclinic`), `tangles_in_case_order` (classes; fails on base on `nested`, `nested_inner_blasts`, `heteroclinic`), `case_heteroclinic` (case sanity); `trellis_is_cut_at_the_pip_orbit` now skips a crossing unless BOTH its keys are branches of the trellis (the selection rule). Removed from KNOWN_ISSUES: `three_tangles_take_three_families` x `hand_built`; added, on the new case: `every_active_class_walks` (WP9 slice 2, DEC-23) and `partition_covers_branch` (WP10: fails on p3 at 16 unstable steps alone). Regressions `test_heteroclinic_pairs_hole_the_stable_owners_partition`, `test_an_unfiltered_outer_blast_iterates_only_its_innermost_bridges`, `test_an_unfiltered_nested_blast_does_not_depend_on_order` (all fail on base), `test_an_unstored_zone_blasts_the_bridges_inside_it`, `test_a_copy_of_a_stored_zone_blasts_like_the_stored_zone` (review), unit `test_rebuilding_one_fixed_points_bridges_keeps_the_others` (fails on base), `test_collision_on_another_stable_manifold_still_punctures` (review; each parameter fails without its gate). Behaviour changes: every other law case and `p3_15` / `p3_16` / `p3_6_blasts` are identical to base (id-free `describe()` signatures, `PYTHONHASHSEED=0`); `nested` / `nested_inner_blasts` keep refs 3, holes 16, classes 13, closed 13, reliable, but p1 is now `A` and tangle 0 (its classes first, registry ids renumbered); after `topology()` of several fixed points the per-fixed-point trellises carry no stage output (`topology(fp)` fills one); the "missing from the branch cycle" WARNING is gone. Suite 1017 passed / 21 skipped / 3 xfailed, 177.6 s (base 1000 / 21 / 2, 162.1 s, same machine); after the review round 1022 / 21 / 3 (the five new tests take ~0.6 s; 211 s measured under a load average of 2.6). Left: DEC-23 (slice 2); escape-lobe images, `partition_covers_branch` and `dual_face_side_is_geometry` at p3 16u, NUM-9 through an unfiltered outer blast's annulus bridges, the p1(u) x p3(s) recipe (WP10).

Verification of the "possibly stale" items (WP0): DOC-3 is confirmed open (the MinimalTrellis "definition is changing", TrellisBranch `return_period`, ResonanceZone "no inversion example", PartitionFamily "concurrent edit" and StrongPip "legitimately None" notes are all still present). PK-7 is confirmed dead code (WP2 deletes it). INV-10, NUM-19 and DOC-4 need a run; they are checked by the WP that owns them.

**WP10 slice 1** (KNOWN_ISSUES 3 -> 2; NOT_APPLICABLE unchanged): numerics and performance leftovers. (1) F6: the Tangle keeps one rtree per stability and queries a segment only against the other, so no same-stability pair is ever formed (`_discard_same_stability` and its INFO line are gone). (2) F3 / NUM-19: crossings sharing one segment are bracketed before cutting, so a bridge contains both its crossings and cdist stays monotone (inversion x period 2 at 10 / 10, p3 13u / 16s, and the p1(u) x p3(s) recipe, which builds now). (3) NUM-9: a curve's existing images are merged chain by chain. (4) NUM-10: iterate failures propagate; `strict` and `skipped` are gone. (5) NUM-14: `_register_bridge` asserts the bridge invariants. (6) PERF-1 / C5: a blast's resolution report counts its new image bridges' same-stability hits only. (7) F7: the strong-pip candidates are classified on first read and a rebuild tests only the remembered pip (k28 1e-5 at 12 blasts: recompute 4.6 -> 0.8 s). (8) INT-10: outer faces from the rotation system; `OuterFaceError` is gone. (9) Laws: `partition_covers_branch`'s last clause is now an element end is open iff a hole opens that half there, the halves derived in the test from each hole's bounds, `bridge_side` and the registry crossing signs (D6), never read from `hole.openings` (the old clause over-claimed on p3 16u, where two adjacent reference pairs hole both sides of the same bridges; a first restatement read `hole.openings` with the partition's own filter and could not fail: the review's break-it, every non-anchor opening flipped in `_hole_openings`, passed it on k10, p3 and k28_two_blasts, and fails it now), `dual_face_side_is_geometry` a local probe (the half-plane test failed on slivers, non-convex faces and a stable hairpin); `(partition_covers_branch, heteroclinic)` XPASSed and was removed. (10) C11 (proposed, under D3): stable-end evidence is incomplete and gives no verdict (k28 at 6-7 blasts reliable; at 8+ the DEC-23 ties remain). (11) NUM-13 / NUM-15: `grow_until_arclength` stops on arc length; the drivers' caps scale with the stretch. (12) TEST-7: C6 pinned. (13) Suite: the law tier builds each case once per run (C12, proposed); `test_deep_reference_survives_noise` read the per-fixed-point trellises, which carry no stage output since WP9 (an empty witness loop), and is rewritten on the selection's trellis (the WP8 frame-fix break-it fails only `test_unlinked_landing_is_read_against_the_pairs_images`, no series test, so it is kept). New tests: `test_blast_through_a_shared_straddle`, `test_crossings_sharing_a_segment_are_bracketed`, `test_every_law_holds_on_the_deep_runs` (`p3_15`, `p3_16`, `p3_6_blasts` and the new builders `p3_16_12`, `inversion_period2_10`), `test_the_count_is_the_listing_of_the_new_curves`, `test_a_blast_counts_and_remembers_nothing` (plus the ninth-blast k28 check under `perf`), the rewritten `test_a_same_stability_pair_is_never_formed`, `test_a_bridge_off_its_branch_is_refused`, `test_every_iterate_failure_escapes`, `test_a_rebuild_tests_only_the_remembered_pip`, `test_outer_faces_do_not_read_areas`, `test_arclength_stops_on_arc_length`, `test_default_cap_scales_with_the_stretch`, `test_the_budget_is_per_call`. Behaviour: law cases, deep runs and the five series are identical to base (id-free signatures, `PYTHONHASHSEED=0`) except `k28_deep_blast_sliver`, where one class's `verified` reads True instead of False (its only disagreement was stable-end evidence; still unreliable). Suite 1042 passed / 21 skipped / 4 deselected / 2 xfailed in 137.3 s (base 1022 / 21 / 3 / 3 in 184.3 s, back to back on the same machine: 0.75x). Deferred: NUM-12 (slice 2), NUM-20, NUM-21 (accepted), `row_of_end_is_geometry` at 16 stable steps (resolution), the E11 2-iteration cells (WP12), vectorised rtree queries (PERF-2 / PERF-3), DEC-23.

**WP9b** (KNOWN_ISSUES 2 -> 2, nothing flips or is added; NOT_APPLICABLE unchanged): the non-crossing rule (D13, DEC-23), on `DualWalk` / `DualGraph`, lookup only. (1) C14: the walk leaves from and arrives at the stable edges the landing SPAN meets (the hole bridges bounding X map onto the kept image bridges at the ends of f(X), which f(b) cannot cross); every tie of k28 at 8-13 blasts and k32 at 8-10 was two walks starting on edges the span does not meet. (2) The marker law: a kept non-loop member `b` with both ends on one stable branch and the stable segment between them bound a disk whose image holds exactly the images of its markers (periodic points at their anchors, holes beside their own bridges on their bridge sides); the word is the shortest walk whose curve, closed along S, holds them (`DualGraph.markers` / `image_marker` / `enclosed`, `DualWalk.lobe_image_markers` / `holds_markers`; detours up to `MAX_DETOUR = 4` crossings; status `no_admissible_walk` when none, `truncated` when the `max_walks` budget, which counts `holds_markers` judgements only, runs out first). C13: it applies to every class with a member lobe, not only ties. `MinimalTrellis.holes`, `Walk.start_node` / `goal_node`, `WalkSearch.constrained` / `rejected`. Behaviour: law cases and the deep runs (`p3_15`, `p3_16`, `p3_6_blasts`, `p3_16_12`, `inversion_period2_10`) have id-free `describe()` signatures identical to base (`PYTHONHASHSEED=0`) except `heteroclinic` (`a` unambiguous with the same word; `c` `a u^-1 a^-1` -> `a d^-1 a^-1` and `bh` empty -> `bi bj bi^-1`, both their registered words; no class `verified` False; ambiguous `[a, c, bn]` -> `[bn]`; still unreliable, DEC-24) and `k28_deep_blast_sliver` (ambiguous `[o, s]` and `s` contradicted -> reliable). Series: k28 at 8 blasts reliable (was ambiguous `[ab, af]`); k32 at 6-7 reliable (`m`'s unique walk started on the wrong edge and its image contradicted it), 8 reliable (closed 44 -> 45), 9-10 no ambiguous class (closed 43 -> 39), still unreliable through INT-11; p3, inversion and nested_outer series unchanged. The two non-tie changes (`bh`, k32's `m`) break the request's "no word changes except where it was ambiguous"; both new words are the registered ones (C13, C14). New laws (pass on all 11, non-vacuous): `walks_agree_with_registered_images`, `member_lobe_classes_have_one_walk`, `walk_starts_on_the_landing` (symbolic), `lobe_enclosure_is_geometry`, `registered_images_keep_markers` (dual graph); series test `test_no_walk_is_ambiguous` (`SeriesStep.ambiguous`); the deep-sliver regression checks its landings, ties and evidence; unit tests of the span, the constrained search and its budget, the `sigma'` rule, and the enclosure (`test_enclosure.py`: interleaved chords cross, the orbit-marker image follows `advance_key`). Fails on base: `heteroclinic`'s symbolic layer, `test_no_walk_is_ambiguous[k28_series, k32_series]`, the deep sliver (the dual-graph laws need the new API). `(every_active_class_walks, heteroclinic)` stays strict, re-tagged DEC-24: the member-less connecting class `bn` has no closed curve with a known image (no heteroclinic crossing has an iterate link, NUM-21). New catalog entry INT-11 (deep evidence disagreements, WP12). `symbolic_dynamics`: heteroclinic 0.04 -> 0.11 s. Suite 1067 passed / 21 skipped / 4 deselected / 2 xfailed in 139.3 s after the review fix (1055 in 139.4 s before it; base 1042 / 21 / 4 / 2 in 138.6 s, back to back on the same machine). Left: DEC-24 (slice 2: close a connecting arc through a linked heteroclinic crossing); classes whose ends lie on two branches (an inversion anchor bridge: sigma would run through the anchor); a walk crossing one edge twice (beyond the midpoint model, kept as a candidate, not observed); INT-11.

**WP11 slice 1** (KNOWN_ISSUES 2 -> 2, nothing flips or is added; NOT_APPLICABLE unchanged): the simplification sweep, a refactor whose proof is that no output moved. Before every commit a harness re-ran the id-free and exact `describe()` signatures of the dynamics, partitions, class table and holes, `is_reliable` and per-logger WARNING counts on the 11 law cases, `p3_15`, `p3_16`, `p3_6_blasts` and five deep builders, every `SeriesStep` of the five series, and every artist of the dual-graph cartoon (both shapes, both colourings, four fixtures): all identical to base. (1) The decisive measurement (SIMP-13, C15): forcing the `cdist_tol` that topology and loom read to 0 or 1e-12 moved nothing anywhere, so the topology now reads no absolute cdist slack (anchors by identity, exact ownership / containment / cuts, relative error only for scaled cdists); tolerance sites in `topology/` 115 -> 59, `cdist_tol` sites 19 -> 2, `logger.warning` calls 56 -> 51. (2) One rule each: `side_after` (INV-11, SIMP-1; orientation read off the FixedPoint, C16), `land_element` as the element image (SIMP-2), `Trellis.same_point` (SIMP-4), `StablePartition.intervals_from_marks` (SIMP-5), `TangleWorkbench.trim_stable_at` (SIMP-10), `_draw_cartoon_graph` (SIMP-11, first half). (3) Dead API deleted (SIMP-12, API-4), `IterateInference` a function, Blast tidy-up (SIMP-7), API traps (API-3, API-6, API-7), the coverage guard keyed by content (TEST-6), and the handed-over nits (`owner_of_intersection` needs its branch, `from_branch` holds the anchor on every branch, `_anchor_branch` inlined, `bridge_classes` without `bridge_ids`, `groupby` in `_diagnose`, a tautological law assertion removed). API changes (D8): removed `invalidate_trellises`, `Trellis.image_of_element`, `SNAP_RTOL`, `StableNode.images`, `infer_iterate_table`, the `IterateInference` class, the branch-cycle wrappers, the dead accessors of SIMP-12, the `tol=` keywords of the strong-pip / pseudoneighbor entry points and `interval_at`, `PartitionFamily.tol`, `bridge_classes(bridge_ids=)`, three `Intersection.synthetic` keywords; changed `reindex_from` (returns None), `owner_of_intersection` (branch required), `from_branch` (anchor on every branch), `BlastResult.fixed_point` -> `fixed_points`, `PointBudgetExceeded.blast_context` (always set), nested slugs (non-default blast counts), FixedPoint shapes, `preserves_orientation` (any period), `land_element` (a non-advanced image branch asserts), `interval_at` (raises instead of a tolerant fallback); added `trim_stable_at`, `ring_order`, `side_after`, `same_point`, `intervals_from_marks`, exports `MaxIterationsReached` / `TURNAROUND_LOG_STRETCH` / `infer_iterates`. Tests rewritten to the new rules: the four orientation tests set the single-step Jacobian, `interval_at` is exact for an estimated cdist, a boundary near 0 gets a synthetic anchor below it, an image on a foreign branch asserts; deleted with their code: the `image_of_element` / `SNAP_RTOL` tests and the `bridge_ids=` test. New tests (each fails on a copy of the base tree): `test_coverage_guard.py`, `test_side_after.py`, the orientation and shape tests of `test_fixed_point.py` (a period-2 orbit included), `test_from_branch_holds_the_anchor_on_every_branch`, the two API-3 tests, `test_stable_trim.py`, `test_henon_slugs.py`. Suite 1088 passed / 21 skipped / 4 deselected / 2 xfailed in 132.0 s (base 1067 / 21 / 4 / 2 in 139.3 s). The coverage guard (base 4f7e815 vs branch, both `--cov`) reports no newly missed line; the old raw-line guard reported 87 on the same pair. Deferred: the `plotting.py` split (SIMP-11), the SIMP-9 feature requests, `partition_element_for` reading the topology and `grow_until_faces_closed`'s cap semantics (D8 API plan), the geometric polyline-clipping tolerances, WARNING -> assert in the walk (frozen until DEC-24), the class-ordering code (the author's rules).

**WP9c** (KNOWN_ISSUES 2 -> 1; NOT_APPLICABLE +5): every bridge is the same (D14). Traced, the member-less connecting class `bn` of `heteroclinic` is no construction error: `bn <- bl <- bi`, and `bi` has two computed members, the arcs of `W^u(p3)` to the heteroclinic crossings on p1's window; `bl` and `bn` are its first two images, member-less by D12's cut (D11), and computed at p3 17u (NEST-1). The defect was that D13 had a rule only for a class with a member lobe. C17 implements D13's own words for every class: a class's arcs lie on its members' unstable branches, or one map step past its namer's; the decided walk of an earlier class on disjoint branches is a wall, and a tie loses the walks whose chords interleave with a wall's on a face ring (`walk_chords`, `DualGraph.chords_cross`, `shortest_walks(walls=)`; `WalkSearch.rejected` counts them); only a walk of multiplicity 1 (one face path) is a wall or is judged by one, in the library and in the law. A member-less class's `tangle` reads its namer's fixed points (the table's rule; a no-op). C18: D14 is tested as "a class connects two tangles only where a computed bridge connects them" (a clause of `alphabet_closed_under_words`), not "every connecting class has a computed member", which contradicts D11. Behaviour: id-free and exact `describe()` signatures identical to base on every law case and deep run except `heteroclinic` (`bn` ambiguous -> resolved, same word; reliable; one DualWalk WARNING fewer; the `describe()` differs on exactly those two lines); the five series identical; p3 17u: the connecting class `bm` resolved to its registered word. `(every_active_class_walks, heteroclinic)` XPASSed and was removed. New: law `walks_avoid_disjoint_images` (symbolic layer; counts class pairs, NOT_APPLICABLE on k10, k28 x2, inversion, inversion_period2), the D14 clause, regression `test_a_tie_no_lobe_decides_is_broken_by_the_images_of_disjoint_arcs`, perf regression `test_registered_connecting_images_agree_with_the_wall_rule` (builder `heteroclinic_17u`), unit tests of the walls, the chord test and the arc branches (all fail on a copy of base). Suite 1102 passed / 26 skipped / 5 deselected / 1 xfailed in 132.0 s after the review fix round (base 1088 / 21 / 4 / 2 in 132 s). Deferred: walls from LATER classes (never needed), same-branch disjointness (measured identical), NUM-20 / NUM-21 (no longer needed for any word).

**WP10b** (KNOWN_ISSUES 1 -> 1, nothing flips or is added; NOT_APPLICABLE unchanged): the code leftovers. (1) NUM-12 / C20: a blast's image curves are computed manifold: the workbench keeps them, every `compute_intersections` re-indexes them query-only (a curve growth has merged into a branch adds no segment and is dropped; `Tangle.add_manifold` returns the count), `create_bridges` cuts them like any curve and holds one crossed fewer than twice as a partial piece (the curve itself, so its `iterated` flag survives); `clear_bridges` releases only id-carrying bridges; `_recut_bridges_after_growth`, its flag clearing and its WARNING are gone. Before, ANY recompute after a blast (plain, `grow_until`, zone define / restore) lost every crossing on an image curve, its bridges and the frontier, and classes dropped (k28 two blasts 4 -> 3, p3_6_blasts 23 -> 9); now a recompute + rebuild changes nothing on any law case, `k28_deep_blast_sliver` or `p3_6_blasts`, and blast / recompute / blast equals blasting straight. (2) NUM-20: the link laws read each side's own saddle; `IterateInference` untouched. (3) PK-4 / C19: foreign chords pair in stable order, the base is read at a non-anchor end, the two "cannot happen" branches are assertions. (4) WP11 review nits: `ManifoldMachine._parabolic_fit` and its pinning test deleted; `CartoonLayout.tol` / `ZoneLayout.tol` / `_cluster` deleted (every cartoon artist identical); `span_contains` / `owns_cdist` docstrings; `test_stable_trim` exact-node and past-the-end cases. (5) NUM-19: one assertion (`_register_bridge`: a `BridgeId` runs in unstable order); the other bullets closed with evidence. Behaviour: id-free and exact `describe()`, partitions, class tables, holes, reliability and WARNING counts identical to base on every law case (`sig.py`, `PYTHONHASHSEED=0`) except `heteroclinic` (iterated partition only: +3 iterated elements, so iterated names, iterated itineraries and refined-pair labels; class table, holes, words, refined rules, matrix and reliability identical; PartitionFamily WARNINGs 2 -> 0) and `inversion_period2` (WARNINGs 7 -> 4, nothing else); `p3_15`, `p3_16`, `p3_6_blasts` the same kind of iterated-name change, 2 WARNINGs fewer each. After a blast, every recompute keeps the blast (above). API (D8): removed `CartoonLayout.tol`, `ZoneLayout.tol`, `ManifoldMachine._parabolic_fit`; `Tangle.add_manifold` returns the count of new segments; `_register_bridge` asserts the id order; a partial piece is the indexed image curve. New tests (each fails on a `git archive` copy of base with the new test files copied over, `.../tmp/wp10b/basecheck`): mutating law `recompute_keeps_the_topology` and the strengthened `recompute_preserves_ids` (every id back; fail on k28_one_blast, k28_two_blasts, nested, nested_inner_blasts, inversion_one_blast, inversion_two_blasts, heteroclinic), `test_growth_after_a_blast.py` (3), `test_a_growth_round_keeps_the_blast` (replaces the flag-clearing test), `test_heteroclinic_iterate_links.py` (fails with the base laws), `test_foreign_chords.py` (8), the rewritten `test_iterated_cut.py` (4), `test_a_reversed_bridge_id_is_refused`; `test_stable_trim.py` +2 (pass on base: they strengthen). Deleted with its code: `test_closed_form_curvature.py`. Suite 1113 passed / 26 skipped / 5 deselected / 1 xfailed in 143.8 s (base 1102 / 26 / 5 / 1 in 137.9 s, back to back on the same machine); `-m perf` 5 passed. Deferred: stable growth un-trims a resonance zone (`grow_n_times` resets every tail of the grown stability; topology unaffected by D10, the zone object goes stale; WP12 / D8 API plan); the "no empty stretch" wording at the anchor lobe (by the author's rule; WP12 docs); the never-firing `BridgeClass` cross-branch skips (D8 API pass); `inversion_period2`'s unreliable class (D9); `row_of_end_is_geometry` at nested 13u / 15s (resolution).
