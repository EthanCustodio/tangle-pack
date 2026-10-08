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
| Intent vs code (INT) | 1 | 7 | 2 | 10 |
| Unnecessary complexity (SIMP) | 0 | 1 | 12 | 13 |
| Resolution & numerics (NUM) | 0 | 11 | 8 | 19 |
| Performance (PERF) | 1 | 1 | 3 | 5 |
| Tests & coverage (TEST) | 0 | 2 | 5 | 7 |
| API traps (API) | 0 | 2 | 5 | 7 |
| Docs drift (DOC) | 0 | 0 | 5 | 5 |
| **Total** | **16** | **36** | **47** | **99** |

The audit raised 18 decision questions (DEC-1..18); all are answered or deferred below. Five issues are marked "possibly stale — verify": INV-10, PK-7, NUM-19, DOC-3 and DOC-4.

The most serious situations:
1. **Class count grows by bookkeeping, not dynamics (BCG-1/2/3).** Every bridge in the trellis is classed, including blast children with no image. Loops are folded through preimage links that blasts never register. Result: +1 active and +1 inert class per blast, while holes, elements and the minimal trellis stay constant.
2. **The hole set depends on depth (BCG-4, PK-1).** Deep tip or blast-child crossings land in the fundamental window and replace the primary pseudoneighbor pairs. The new reference punches about log(c_u)/log β holes, so p3 goes 12→42 holes at blast 6 or at 15 steps. That makes p3 unreliable: an unreachable class and a virtual `new1`.
3. **Inversion is wrong by construction, and no fixture can show it (INV-1/2/4/5/6).** Four anchors sit at one point. The pseudoneighbor lap factor is λ_u instead of λ_u². The zone ring is keyed by orbit index. The side and partition model is per branch. The only fixture has no holes.
4. **The heteroclinic path has never run (NEST-1/2).** No fixture has an A(u)×B(s) crossing. A hand-built case is unreliable, and heteroclinic pseudoneighbors are rejected by construction.
5. **The dual-graph walk fails where it must succeed (INT-1, PK-1).** Walks come back "unreachable" or "ambiguous" and virtual `new*` classes appear. Under D3 the walk IS the word and must always exist, so each of these is a construction bug (WP7).

---

## Inversion end-to-end

### INV-1: Four coincident (0,0) anchors at an inversion point instead of one degree-4 node
**Status:** open | **Severity:** high | **Related:** INV-3, INV-6, DEC-4, API-5, SIMP-12 | **WP:** WP3

`_register_anchors` registers one anchor per (unstable branch, stable branch) pair. The only filter is on orbit index, so an inversion point gets two anchors per unstable branch, four in all. The author's rule is one.

This cascades:
- 5 strict xfails: `one_anchor_per_unstable_branch`, `anchor_bridge_class_leads_its_tangle`, `arrangement_euler` (V−E+F=0), `itinerary_pairs_same_side`, `landings_contained`.
- `_image_chain` collects both anchors, which creates the phantom pair (0,1). On the k=4 probe this gives `new1 = {L_(0.0;1), L_(0.1;1)}`, and anchor bridges `b`/`h` are split by which anchor id they got.
- The singleton `[0,0]` owns two crossings, and the dual graph builds a zero-length stable edge.

**Evidence:** `src/tanglepack/numerics/TangleWorkbench.py:1537-1616` (:1574-1599, filter :1598, docstring :1575-1580); `src/tanglepack/numerics/Tangle.py:9-16`; `src/tanglepack/numerics/Tangle.py:1050-1072` (`_both_anchors`); `src/tanglepack/topology/Arrangement.py:80-88`; `src/tanglepack/topology/BridgeClass.py:748-758`; `tests/cases.py:796-815`; `tests/conftest.py:99-104`; `docs/regions_bridge_identity_and_cleanup_plan.md:334-337`; `stress/numerics/report/report.tex:2012-2014`; probe `inv_probe.py` (zone+1 blast: `new1 b u^-1 c^-1` vs walk `h u^-1 c^-1`)

### INV-2: Pseudoneighbor uses bare λ_u as the branch-return factor (should be β^k_value = λ_u²)
**Status:** open, confirmed by probe | **Severity:** high | **Related:** SIMP-3, INV-4 | **WP:** WP2

- `_resolve_reference_end` scales by `/lambda_u` and `*lambda_u`.
- `_unstable_interval_is_clear` uses `c*lambda_u**m` per k-step lap, with `log(lambda_u)` for the m-range.

A lap is k_value = 2·period steps, which is β^k = λ_u². Patching this in memory changes the pairs:

| Variant | Code as is | With β^k |
|---|---|---|
| No zone | `[(23,4)]`, 7 holes | `[(11,9),(19,4),(25,17)]`, 16 holes |
| Zone | `(8,5)` | `(8,7)` |

The current fixture never reaches this code because it has no pairs. `TrellisBranch.eigenvalue` and `IntersectionRegistry.on_interval` have the same defect.

**Evidence:** `src/tanglepack/topology/Pseudoneighbor.py:51-60` (OPEN Dev Note), `:491-492`, `:551`, `:571-576`; `src/tanglepack/topology/TrellisBranch.py:79`; `src/tanglepack/numerics/IntersectionRegistry.py:656`; `docs/regions_bridge_identity_and_cleanup_plan.md:330`; probe `inv_probe.py`

### INV-3: An anchor-only stable branch at 4-5 steps raises AssertionError "owned by 0 partition elements" (WP4-B)
**Status:** decided 2026-10-07 / deferred | **Severity:** high | **Related:** INV-1, DEC-5 | **WP:** WP3

`_build_intervals` gets a single boundary, so `_element_of_intersection` asserts. This reproduces at k=10 with 4-5 steps, and at k=3.5 and k=4 with 4-8 steps. There are two options: a closed singleton at the anchor, or one element `[anchor, end]`. Either one has to come after INV-1 and dedupe by point. Once decided, add `test_inversion_before_the_tangle_forms`.

**Evidence:** `src/tanglepack/topology/StablePartition.py:214-235` (DEFERRED Dev Note), `:662`→`:2147`; `stress/numerics/report/report.tex:941, 956-967, 2009-2014`; `stress/numerics/report/numbers_extra.tex:297`; `stress/numerics/e07_breaking.py:312-318`; probe `inv.py` (steps 4, 5 ERR)

### INV-4: No usable inversion fixture: the k=10 saddle escapes, has no pseudoneighbors and no holes, so most laws are N/A
**Status:** open | **Severity:** high | **Related:** INV-2, INV-11, TEST-1, DEC-4 | **WP:** WP1

`build_inversion` grows k=10 to 6/5 steps:
- 7 steps raises PointBudgetExceeded, and points jump from 196 to 49,417 between steps 5 and 6. The window between "too shallow" and "escaped" is one step wide.
- At 6 steps: 0 holes, 0 kept minimal-trellis bridges, 4 classes, `new1`, `is_reliable=False`.
- 15 or more laws are NOT_APPLICABLE on inversion, so StrongPip, Pseudoneighbor, StablePartition, ResonanceZone and Arrangement all say "inversion unvalidated".
- No E11 case and no E6 sweep point gives a nontrivial inversion tangle. In the sweep, k≤5.70 has 0 crossings, 6.92 has 1, 8.40 has 7, and k≥10.19 hit OOM.

**Candidate fixture (Case B):** the k=4 other root `(1-√5, √5-1)`, λ -1.96/-0.51, `area_cutoff=1e-5`, 11/11 steps or `grow_until_bounded`, a zone at f(q0), 0/1/2 blasts. Laws pass 70/75, 71/75 and 70/75; step 12 escapes.

**Evidence:** `tests/cases.py:568-636`, `:843-875` (`_NO_HOLES`); `src/tanglepack/topology/StrongPip.py:63-72`; `src/tanglepack/topology/Pseudoneighbor.py:37-40`; `src/tanglepack/topology/StablePartition.py:236-238`; `src/tanglepack/loom/ResonanceZone.py:37-39`; `stress/numerics/report/report.tex:210, 830-832, 940-991, 1150-1156`; `stress/numerics/report/numbers_extra.tex:219-221, 287, 301`; `stress/numerics/results/e06_parameter_sweep.jsonl` (which='other'); probes `inv.py`, `inv_probe.py`

### INV-5: Resonance-zone ring is matched by orbit index only, so the inversion zone polygon is wrong; the logic is duplicated in TangleSession
**Status:** open, confirmed by probe | **Severity:** high | **Related:** INV-4, DEC-4, NEST-4 | **WP:** WP3

`_build_boundary` builds `by_u_orbit = {manifold_a_key[2]: i}`. Both inversion cuts have orbit 0, so u0.1 is used twice and u0.0[p,q0] is never used.

- The probe zone area is 1.232. The correct single-pip loops are 0.705 and 0.704 (sum 1.41).
- The zone's own boundary bridge (1,28) is classified "no zone", so the choice of blast interior bridges is suspect.
- Zones are keyed `(fixed_point, branch_index)`, which gives "two zones" for one chain.
- `TangleSession._zone_ring_order` re-implements the same hop.

**Evidence:** `src/tanglepack/loom/ResonanceZone.py:223-252`, `:121-133`, `:37-39`; `src/tanglepack/loom/TangleSession.py:2052-2076`, `:155`; probe `inv_probe.py`

### INV-6: The same-side law and the partition are per branch, while L(0.0) and R(0.1) are one global half-plane through p
**Status:** open | **Severity:** high | **Related:** INV-1, DEC-4, INT-8 | **WP:** WP3

Both inversion stable branches are oriented anchorward, so they point in opposite directions along one curve. Several pieces have no inversion-aware side map:
- `cross_side = first.side != second.side`
- the `itinerary_pairs_same_side` law
- `oriented_class`
- the cross-branch row checks

As a result, legitimate `(R(0.1), L(0.0))` pairs are flagged as cross-side, and the WARNING "class a: word b^-1 disagrees with the registered image" appears. The elements next to p on one global side (L(0.0;1) and R(0.1;1)) are not glued. Landing spans also straddle cuts (`contained=False`).

**Evidence:** `src/tanglepack/topology/SymbolicDynamics.py:474`; `tests/helpers/laws.py:1885`; `tests/cases.py:808-815`; `CLAUDE.md` "a bridge only connects same-side elements"; probes `inv.py`, `inv_probe.py`

### INV-7: `branch_index` means two things (inversion half-chain vs an independent second chain); the initializer duplicates branch ordering
**Status:** open | **Severity:** medium | **Related:** INV-9, SIMP-9 | **WP:** WP3

- For inversion, branch 1 is the second half of one k_value chain.
- Otherwise it is a separate chain seeded with `num_branches=2`.

Both use ManifoldKey slot 3, so every consumer special-cases `check_inversion`, and `num_branches` means different counts in FixedPoint and in the initializer. `construct_kevin_way` still uses the old `get_iterable_array` ordering ("0, p−1, …, 1") beside `advance_key`, with two near-identical loops.

**Evidence:** `src/tanglepack/numerics/FixedPoint.py:151-183, 453-456, 484-492, 497`; `src/tanglepack/numerics/ManifoldMachine.py:286-298`; `src/tanglepack/numerics/ManifoldInitializer.py:25-30, 153, 432, 483-598, 527, 530-598`; `src/tanglepack/numerics/TangleWorkbench.py:352, 521-548, 810, 1000`

### INV-8: `orient_eigenvectors` is a no-op for inversion; per-orbit-point eigenvector propagation is skipped for period>1 inversion orbits
**Status:** open | **Severity:** medium | **Related:** INV-7, TEST-1 | **WP:** WP3

Period-1 inversion is fine. For a period>1 inversion orbit, the orientation along the orbit is never set.

**Evidence:** `docs/regions_bridge_identity_and_cleanup_plan.md:331-332`; `src/tanglepack/numerics/TangleWorkbench.py:~529-536`; `src/tanglepack/numerics/FixedPoint.py:10-21` (inversion eigenvector storage Dev Note)

### INV-9: `grow_n_times(branch_index=0)` default logs "branch_index is ignored" on every inversion growth step
**Status:** open | **Severity:** low | **Related:** INV-7 | **WP:** WP3

This is log noise; the probe saw it three times. Fix by defaulting to `branch_index=None` or skipping the warning for inversion.

**Evidence:** `src/tanglepack/numerics/TangleWorkbench.py:590`; `src/tanglepack/numerics/ManifoldMachine.py:286-289`

### INV-10: `iterate_bridge` on an anchor-rooted inversion bridge raised TypeError in `_insert_after` (branch_index None)
**Status:** possibly stale — verify | **Severity:** medium | **Related:** INV-1 | **WP:** WP3

`create_bridges` now refuses the zero-span cut between two anchors (inversion bridges went 9→7). Nobody re-checked whether a non-degenerate anchor-rooted bridge still fails.

**Evidence:** `docs/regions_bridge_identity_and_cleanup_plan.md:196, 332-334`

### INV-11: Orientation-reversing maps (det J < 0) are rejected, so every side-flip path is unexercised
**Status:** decided 2026-10-07 | **Severity:** medium | **Related:** SIMP-1, DEC-10, TEST-1 | **WP:** WP11

`set_k_value` rejects disagreeing eigenvalue signs, so all 8 b=-1 configs are rejected and there is a strict xfail with `raises=ValueError`. The parity flip on carried hole sides (WP7-D2), `check_holes_share_bridge_side(orientation_preserving=False)`, `land_element` and `image_of_element` flips have never run on real data.

Note that inversion has det J = +1. Earlier notes conflate "inversion" with "orientation-reversing", so the plan must say which one is meant. Supporting det<0 would need a per-stability k_value.

**Evidence:** `src/tanglepack/numerics/FixedPoint.py:185-225`; `src/tanglepack/numerics/FixedPointSolver.py:16-22`; `tests/cases.py:36-39, 826-830`; `tests/conftest.py:108-110`; `src/tanglepack/topology/Trellis.py:32-40`; `stress/numerics/report/report.tex:833-834, 1906`; `docs/dual_graph_symbolic_dynamics_plan.md:338`

---

## Period-k orbits

### PK-1: p3 at ≥15 unstable steps or ≥6 blasts is unreliable: an unreachable class, a virtual `new1`, and "crossing abuts no hole"
**Status:** open, strict xfail `open_p3_deep_runs` (`p3_15`, `p3_16`, `p3_6_blasts`) | **Severity:** high | **Related:** BCG-4, PK-3, PK-5, NEST-3, DEC-1 | **WP:** WP7

At 11-14 steps and 0-5 blasts: 9 classes, 3 active, reliable. At 15 steps there are four tip crossings with neither a −1 nor a +1 link inside the ray-1.0 window:
- They break the adjacency of `(18,34)`.
- They disqualify `(4,20)` through `M^-17` of tip crossing 24.
- The only references left are the tip pairs `(31,9),(9,15)`.

The anchor image chain `0→pip→4→20` then finds no hole at 20/4, so `PartitionFamily.py:800` warns "abuts no hole" and does not cut. That produces wrong elements, then the unreachable walk `0.0/R#2→R#3`, then `new1`.

E11 has exactly 1 unresolved and 1 virtual class at every blast from 6 to 20 (p3 and nested_inner); blasting never resolves them. After the bursts, every unresolved class is "dual-graph walk unreachable" (127 entries, 0 singleton). k28 ends with 12 unresolved, 1 virtual and 3 ambiguous; k32 with 9/0/2. The report suggests stable growth instead of blasting (untested).

**Evidence:** `tests/cases.py:831-839`; `tests/regression/test_open_p3_deep_runs.py:1-18, 58-67`; `src/tanglepack/examples/henon_cases.py:30-39`; `src/tanglepack/topology/PartitionFamily.py:800`; `src/tanglepack/topology/DualWalk.py:757-760`; `stress/numerics/report/report.tex:1477-1527, 1734-1735, 2025-2033`; `stress/numerics/report/numbers_extra.tex` (NumDeepUnresolvedUnreachable=127); `stress/numerics/report/numbers.tex:468-472, 576-580, 678-680`; `stress/numerics/results/e11_blast_depth.jsonl` (p3 1e-5 `invariant_assert` blasts 6-19); `.refactor/knowledge/topology.md` (OPEN p3 ≥6 blasts); `.refactor/decisions.md` 2026-10-02

### PK-2: The provisional rule that punches forward iterates +1..+(k_value−1) (`_is_forward_beyond_fundamental`)
**Status:** provisional / decided 2026-10-07 | **Severity:** medium | **Related:** PK-3, BCG-4, DEC-6 | **WP:** WP5

This is how the other k−1 branches get their window holes. Only the "≥k_value not punched" half is tested. The chalkboard alternative: "a reference pair is any consecutive pair inside a window F_j on any ray", which needs no table and no exemption.

A side effect on the baseline is the inert chain `u(1.0)→v(2.0)→x(0.0)→loop`. It is not Z3-equivariant (x does not map to u) because the image of the hole at iterate k−1 is not punched. That is by design.

**Evidence:** `src/tanglepack/topology/StablePartition.py:173-190, 1712-1723`; `tests/conftest.py:72-77`; `CLAUDE.md` "Holes propagate BACKWARD only"; memory `holes-backward-only-2026-09`

### PK-3: When `extend_pseudoneighbor_trajectories` finds no table links, the other k−1 rays silently get no window holes
**Status:** open (suspected bug) | **Severity:** high | **Related:** PK-2, BCG-3, PK-1 | **WP:** WP5

Probe p3 at 6 blasts: `pairs=2` and zero forward pairs, so rays 0.0 and 2.0 got no direct window holes (at 5 blasts they had +1/+2 holes). The partition becomes asymmetric across the orbit, and the provisional rule is never even applied.

**Evidence:** chalkboard review §4 probe (p3 6 blasts); `src/tanglepack/topology/StablePartition.py:1712`

### PK-4: Cross-branch image-bridge cuts are skipped with a WARNING and unvalidated
**Status:** open | **Severity:** medium | **Related:** PK-1, INV-6 | **WP:** WP7

These cases only warn and do nothing:
- an empty stretch whose far end is on another stable branch ("no flank")
- a chain that ends on another branch (marks nothing)
- an odd number of interior crossings
- "no empty stretch"

`BridgeClass` cannot place image endpoints that sit on different unstable branches. The only checks are ownership invariants. This matters for period-k and inversion, where image chains switch branches.

**Evidence:** `src/tanglepack/topology/PartitionFamily.py:714-717` (Dev Note), `:727, 774, 799, 813, 941, 955, 976`; `src/tanglepack/topology/BridgeClass.py:735, 798`; `CLAUDE.md` "Partition families"; memory `minimal-dual-graph-2026-09` ("cross-branch cut is the open point")

### PK-5: The fundamental segment is unified only on the pip's branch; the other k−1 branches are walls everywhere
**Status:** decided 2026-10-07 (kept 2026-09-30) | **Severity:** low | **Related:** PK-1, DEC-6 | **WP:** WP7

A walk that needs to cross another branch's window becomes "unreachable". This is one plausible contributor to the p3 unreachable class.

**Evidence:** `src/tanglepack/topology/DualGraph.py:583-611`; `src/tanglepack/topology/DualWalk.py:757-760`

### PK-6: A strong pip with fewer cut points than k_value only warns
**Status:** open | **Severity:** low | **Related:** PK-3 | **WP:** WP5

If the pip's iterates are not registered on all branches, the windows on the remaining branches stay uncut.

**Evidence:** `src/tanglepack/topology/Pseudoneighbor.py:444-450`

### PK-7: Crossings with `manifold_a_key is None` are tested at every cycle residue (can reject a true pair)
**Status:** possibly stale — verify | **Severity:** low | **Related:** SIMP-12, DOC-3 | **WP:** WP2

Verified 2026-10-07: blast children carry their key (`BridgeIterator.py:334-335`, `:388`), and `Tangle._key_of` asserts a key on every unstable curve. So the `None` branches look dead, not conservative. The Dev Notes still say "legitimately None".

**Evidence:** `src/tanglepack/topology/Pseudoneighbor.py:84-93, 561-567, 591`; `src/tanglepack/topology/StrongPip.py:84-87, 232`; `src/tanglepack/numerics/IterateInference.py:89-94, 254-361`; `src/tanglepack/numerics/Tangle.py:220-228`

### PK-8: Noisy WARNING "inert class u has the non-empty word v" on period-3 runs
**Status:** open | **Severity:** low | **Related:** PK-2 | **WP:** WP5

For period k it is expected that an inert class maps to the next branch's inert class, so this warning is noise.

**Evidence:** `src/tanglepack/topology/SymbolicDynamics.py:898`; `.refactor/decisions.md` 2026-10-02

### PK-9: The hole-chain truncation at `max_steps=50` is logged only at DEBUG
**Status:** open | **Severity:** low | **Related:** BCG-4 | **WP:** WP5

A deeper reference chain would be cut silently.

**Evidence:** `src/tanglepack/topology/StablePartition.py:437-438`, `~583`

---

## Nested / heteroclinic tangles

### NEST-1: No fixture has a heteroclinic crossing; connecting classes, the merged outer node and cross-tangle walks have never run, and a hand-built case is unreliable
**Status:** open | **Severity:** high | **Related:** NEST-2, INT-10, TEST-1, DEC-7 | **WP:** WP9

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
**Status:** open / decided 2026-10-07 | **Severity:** high | **Related:** NEST-1, DEC-7 | **WP:** WP9

Two consecutive heteroclinic crossings on B's stable branch carry A's unstable key, which is not in B's branch cycle, so they are rejected. No hole can ever sit in a heteroclinic lobe.

**Evidence:** `src/tanglepack/topology/Pseudoneighbor.py:535-540`

### NEST-3: Nested is reliable only for outer ≥2 and inner ≤4 blasts; nested_outer loses classes and holes as it is blasted
**Status:** open | **Severity:** medium | **Related:** PK-1, BCG-5 | **WP:** WP7

- **0 outer blasts:** the outer class `e` is a singleton with no registered chain.
- **Inner ≥6 blasts, or p3 at 15 steps:** the run inherits PK-1 exactly (27 classes, 42 p3 holes); E11 nested_inner goes 13→55.
- **nested_outer:** classes 12→13→14→11, holes 16→12→6, active 6→3 and `minimal_image` 3→0 at blast 10. Unexplained.

**Evidence:** `src/tanglepack/examples/henon_cases.py:40-47`; `.refactor/knowledge/loom.md` ("Reliable iff outer ≥2 and inner ≤4"); `CLAUDE.md` ("Nested needs two outer blasts"); `stress/numerics/results/e11_blast_depth.jsonl` (nested_outer); `stress/numerics/report/report.tex:1286-1288`

### NEST-4: Blast containment is by zone polygon only; the library still lets an outer blast iterate the inner tangle
**Status:** open (the 2026-10-02 fix was recipe-only) | **Severity:** medium | **Related:** NUM-9, DEC-11 | **WP:** WP9

`build_nested` passes `fixed_point=[zone.fixed_point]`, but `blast_zone` uses `zone.contains_point` plus whatever filter the caller passes. A caller can reproduce three failures:
- order dependence
- the `merge_manifolds` crash
- silent p3 corruption at ≥6 outer blasts

**Evidence:** `src/tanglepack/loom/Blast.py:278-292`; `.refactor/knowledge/loom.md` (Smells); memory `hole-side-own-blast-2026-10-02`

### NEST-5: `rebuild_bridges(fp)` clears every fixed point's bridges
**Status:** open (verified 2026-10-07: documented, not fixed) | **Severity:** medium | **Related:** NEST-1 | **WP:** WP9

The docstring says "The CLEAR is global even when fixed_point is given". `clear_bridges()` is at :1849, so a per-fixed-point rebuild in a nested session drops the other tangle's bridges. Whether partial bridges accumulate per re-blast was not re-checked.

**Evidence:** `src/tanglepack/numerics/TangleWorkbench.py:1796-1851` (:1831-1833, :1849); `docs/regions_bridge_identity_and_cleanup_plan.md:333`

### NEST-6: `entry.tangle` follows the combined trellis order `[A(p3), B(p1)]`; `TangleBuild.fixed_points` is "outermost first" `[B, A]`
**Status:** open | **Severity:** low | **Related:** — | **WP:** WP9

The two orders are opposite, which invites mix-ups. Align them or document the difference.

**Evidence:** period-k/nested investigation §3.5; `src/tanglepack/examples/henon_cases.py`

---

## Bridge-class growth

### BCG-1: Every trellis bridge is classed (blast children included), and members with no image count as active, giving +1 active class per blast
**Status:** open | **Severity:** high | **Related:** BCG-2, BCG-3, BCG-6, DEC-2, TEST-2 | **WP:** WP6

`bridge_classes` defaults to all non-partial bridges, and the session passes no `bridge_ids`. The whole table reaches `symbolic_dynamics` and the transition matrix.

From blast 8, p3 has these counts constant: holes 44, homotopy elements 85, iterated 90, `minimal_kept` 28, `minimal_image` 4. Yet classes rise 25→51, +2 per blast, which is +1 active (this issue) and +1 inert (BCG-2).
- The newest class is the image of the previous blast's newest: (89,90)←(79,80)←(71,72)←(63,64).
- These march one element pair inward: `1.0/R#21↔R#25 → 2.0/R#19↔R#23 → 0.0/R#16↔R#20`.
- The canonical rules form a shift chain `A3→A4→…→A37→new1`.

Classes with a minimal-trellis member stay at 25 (p3 blasts 8-16, while the total goes 27→43); k28 stays at 14-20 between switches. The fix is to build letters from `minimal.kept_bridge_ids` and keep blast children as evidence only. E7 shows the same thing with depth: p3 classes 9 at 11-14 steps, then 23/30/36 at 15/16/17.

**Evidence:** `src/tanglepack/topology/BridgeClass.py:616-622, 789-792, 833`; `src/tanglepack/loom/TangleSession.py:391`; `src/tanglepack/topology/SymbolicDynamics.py:1251`; `stress/numerics/results/e11_blast_depth.jsonl` (p3/nested_inner 1e-7); `stress/numerics/report/report.tex:1238-1246, 1284-1291`; `stress/numerics/report/numbers.tex:646-666, 854-866`; `stress/numerics/report/tables/e11_summary.tex`; `e07_breaking.jsonl`; probes `probe_p3.py`, `minimal_classes.py`

### BCG-2: Loops are folded by following preimage links, so a loop born from a fold stands alone as an inert BridgeClass(x,x), giving +1 inert class per blast
**Status:** open | **Severity:** high | **Related:** BCG-3, INT-1, DEC-2 | **WP:** WP6

`_non_loop_ancestor` walks `iterate(x,−1)`. Crossings created by a fold have no −1 link, so the chain breaks. Example: loop (61,62) with preimages `[None, 55]` goes to its own class, and its images (69,70), (77,78) each stand alone in a different element (R#21, R#19, R#16…).

The parent's `_image_chain` does contain (61,62), and the author's rule ("the image of a bridge that did connect two elements") is the forward relation, so the intent did not make it into the code. Fix: fold via the bridge whose image chain contains the loop, falling back to `BridgeIterator.preimage_bridges`. Put true orphans in one bucket per orbit.

**Evidence:** `src/tanglepack/topology/BridgeClass.py:636-647, 706, 842-875` (:863-866); E11 p3 `n_inert` 1,1,1,2,3,…,13 over blasts 6..20; probe `probe_p3.py`; `CLAUDE.md` "a loop with no registered ancestor stands alone"

### BCG-3: Blasts never re-run registry-wide iterate inference, so interior image crossings stay unlinked
**Status:** open, confirmed by probe | **Severity:** high | **Related:** BCG-2, INT-6, PK-3 | **WP:** WP2

- `iterate_bridge` only calls `register_endpoint_iterates`.
- `infer_iterates` runs only inside `compute_intersections`.
- The recipe calls it once before the blasts, and `Blast.py` never calls it.

Example: p3 crossing 61's true preimage 54 matches to 0.3-0.5% but `iterate(54,1)=None`. Running `session.workbench.infer_iterates()` after `build_period3(blasts=10)` adds 6 links: 31→29 classes and 2→0 standalone loops (another run: 29/27 active → 28/24 active).

**Evidence:** `src/tanglepack/numerics/BridgeIterator.py:393`; `src/tanglepack/numerics/TangleWorkbench.py:1380`; `src/tanglepack/numerics/IterateInference.py:119, 377`; `src/tanglepack/loom/Blast.py:195`; `src/tanglepack/examples/henon_cases.py:146`

### BCG-4: The reference pair switches when deep crossings fold onto the fundamental window, and the new backward hole chain has length ~ log(c_u)/log β
**Status:** open / decided 2026-10-07 | **Severity:** high | **Related:** PK-1, PK-2, PK-9, INT-2, BCG-5, DEC-1 | **WP:** WP8

**p3, blast 6:** blast-child crossings 41-43 (c_u ≈ 2.6e3-3.1e3, ray 2.0) land inside the window (0.397, 1.481] of ray 1.0.
- The primary pairs `(18,29),(17,18)` stop being consecutive.
- The references become `(41,42),(42,43)`, each propagated back 20-21 bridges with opposite sides on the same bridges.
- Holes go 12→42 (origins 23 + 21), elements 24→81, with about 25 right-row elements per ray and a staircase `d→e→…→t` of 17 one-letter classes.
- `(42,41)` is a sub-arc of f¹⁷ of the old reference bridge `(17,18)`. So the switch also undoes the backward-only rule: old forward iterates (≥k) now get holes as negative iterates.

**k28, blast 8:** new references `(184/185,185/186)` (unstable cdist ≈2.36e8, resolution ratio ≈320, so not polygon noise). Their preimages are not registered because they are trimmed away, so the membership test runs on scaled cdists only. Holes 4→24 (3 direct + 21 carried), classes 4→20. Another switch at blast 11 `(892,893),(899,900)`.

**k32, blast 8:** classes 14→37, holes 16→56.

Classes at fixed depth also depend on the cutoff (NUM-5).

**Evidence:** `src/tanglepack/topology/StablePartition.py:433, 547-600, 1712-1723, 1771`; `src/tanglepack/topology/Pseudoneighbor.py:131, 505-620` (all-iterates test from :530); `stress/numerics/report/report.tex:1222-1230, 1277-1295`; `stress/numerics/report/numbers.tex:438-456, 540-564`; E11 k28 `n_holes [4x8,24,24,20,26,…]`; probes `p3probe2.py`, `probe_k28.py`, `p3_holes.py`, `k28_classes.py`; memory `hole-side-own-blast-2026-10-02:15`

### BCG-5: Holes and classes decrease under further blasting (the partition is non-monotone in depth)
**Status:** open | **Severity:** medium | **Related:** BCG-4, NEST-3 | **WP:** WP8

- nested_outer: classes 13→11, holes 16→12, `minimal_image` 3→0 at blast 10.
- k32: holes 58,50,42,34,13 over blasts 9-13 (classes 50→34, iterated elements 116→42), then 0, then 36.
- p3: holes 44→23 at blast 23.

This is consistent with the confirmed full-orbit disqualification, but the report says the burst and drop positions are "not checked".

**Evidence:** `stress/numerics/report/report.tex:1247-1252, 1285-1293`; `stress/numerics/report/numbers.tex:560-564, 770-774`; E11 k32/nested_outer series; memory `pseudoneighbor-blast-semantics`

### BCG-6: The class table is not quotiented by the map or the orbit shift; one class per hole-ladder rung, and period-3 classes come in shifted triples
**Status:** open / decided 2026-10-07 | **Severity:** medium | **Related:** BCG-1, BCG-4, DEC-2 | **WP:** WP6

A ladder `E_n→E_{n+1}` produces a class `a_n` for each rung. Transport chains (`a1→a2→…→a20`, single-token rules) could be shown as one symbol plus a delay. Period-3 words run `A0→A1→A2→A0`. The growth is bounded by the hole-chain depth (about 20 blasts to close), but it is linear by construction.

**Evidence:** class-growth investigation §4.5; chalkboard §4; E11 p3 blast-20 `canonical_rules`

### BCG-7: Class identity is not stable across blasts: positional `element_id` re-indexes, so every name and letter churns
**Status:** open | **Severity:** medium | **Related:** BCG-4, DEC-9, TEST-2 | **WP:** WP8

A new hole renumbers every element. None of the 9 class names from blast 0 survive p3 at 6 blasts, the words change letter every blast, and k28 letters reach `cz`. The report did not check whether the unresolved class is the same class across blasts. The chalkboard identity of an element is its bounding holes `(origin, iterate, side)`.

**Evidence:** `.refactor/knowledge/topology.md` (OPEN); `src/tanglepack/loom/TangleSession.py:375-382`; E11 k28 `letters_final`; `stress/numerics/report/report.tex:1478-1479`

---

## Intent vs code (chalkboard review)

Cross-references, catalogued elsewhere: loop folding by preimage (BCG-2), element identity (BCG-7), window reference rule (PK-2), inversion anchor (INV-1), λ_u lap factor (INV-2).

### INT-1: The word is the shortest dual-graph walk between midpoint landings, not the homotopy class of the registered image chain
**Status:** open / decided 2026-10-07 | **Severity:** high | **Related:** INT-8, PK-1, DEC-3 | **WP:** WP7

The homotopy class of f(bridge) relative to the holes is determined by its crossing sequence. That is `trellis_itinerary`, which is used only for singletons and as "evidence". The answer comes instead from `land_element` midpoint probes plus `shortest_walks`. Shortest equals homotopic only where the dual graph is tree-like; elsewhere the result is `ambiguous`, `unreachable` or virtual `new*`, and the code admits "the shortest walk can be topologically wrong where nothing is registered".

**Resolved by D3 (2026-10-07): rejected.** The shortest walk stays the definition and must always exist; the registered chain is the cross-check and the word for singletons. The audit's original proposal is kept for the record: the trellis itinerary of `_image_chain` becomes THE word, tightened by cancelling a a⁻¹. Classes with no registered image become unresolved, and a loom `resolve_classes()` iterates one representative per class. That would delete the walk search, the ambiguity handling and `registered_evidence/verified`, and shrink `PartitionFamily`/`MinimalTrellis`.

**Evidence:** `src/tanglepack/topology/DualWalk.py:273, 350, 672, 827`; `src/tanglepack/topology/SymbolicDynamics.py:583, 1094, 1169`; `src/tanglepack/topology/PartitionFamily.py:684, 900`; `src/tanglepack/topology/MinimalTrellis.py:243-416`

### INT-2: A hole is punched on the whole containing bridge, not on the backward-image sub-arc
**Status:** open / decided 2026-10-07 | **Severity:** medium | **Related:** BCG-4, PERF-4, DEC-1 | **WP:** WP4

Every ancestor bridge becomes a hole bridge and every ancestor endpoint a partition boundary (42 holes → 85 elements). `_containing_bridge` scans with `rtol=1e-3` slack and "tightest span wins"; bridges tile a ray, so a bisect would be exact.

**Evidence:** `src/tanglepack/topology/StablePartition.py:1771-1800`

### INT-3: Hole openings and partition rows use the GEOMETRIC row while classes, dual graph and walks use the COMBINATORIAL row
**Status:** open / decided 2026-10-07 | **Severity:** medium | **Related:** INT-4, INT-5, NUM-1, NUM-2, DEC-8 | **WP:** WP4

`_hole_openings` uses `_row_at`, `_stable_frame`, `_bridge_end_geometry`, `_arc_side_of`, nearest-node tests and special trimmed-end/anchor/`inward=True` branches, about 450 lines. Sign arithmetic decides the same thing: `row = row_of_end`, `outward ⇔ (sign>0) == (side==left)`. `measured_side`, `SideDiagnosis` and `RowDiagnosis` exist only to diagnose that geometry.

The two readings disagree on coarse polylines. E11 nested_inner with 2 iterations, blasts 16-20, fails "approaches its two defining crossings from opposite sides" while the combinatorial rows agree (2-node polylines at about 0.2°).

**Evidence:** `src/tanglepack/topology/StablePartition.py:714, 755, 790, 916, 982, 1086-1475, 1236-1285, 1476-1582`; `stress/numerics/e11_blast_depth.py:76-95`

### INT-4: `_hole_openings` uses the nearest vertex of the whole polyline (the same defect already fixed for propagated hole sides)
**Status:** open (latent) | **Severity:** low | **Related:** INT-3 | **WP:** WP4

The same kind of misread broke p3 hole sides at 15 steps / 6 blasts.

**Evidence:** `src/tanglepack/topology/StablePartition.py:~1480-1500`; `.refactor/knowledge/topology.md` ("Latent (believed)"); `.refactor/decisions.md` 2026-10-02

### INT-5: `_arc_side_of` returns None for a 2-point span slice, giving hole `bridge_side` None with a misleading DEBUG
**Status:** open (latent; 175/175 observed reads had ≥3 points) | **Severity:** low | **Related:** INT-3 | **WP:** WP4

Fix: widen the slice to at least 3 points.

**Evidence:** `src/tanglepack/topology/StablePartition.py:794-795`

### INT-6: The iterate table is filled by a 5% relative-cdist heuristic, not registered exactly when a crossing is created
**Status:** open | **Severity:** medium | **Related:** BCG-3, NUM-8, SIMP-3 | **WP:** WP2

`infer_iterates` uses `cdist_rtol=0.05` and has no map gate (a 1e-3 phase-space gate exists only on the endpoint path). Deep blasts put real crossings much closer than 5%, so two crossings can be mis-linked, and pips, pairs, holes and chains all trust the link as exact. The topology re-implements the scaling fallback three more times.

The machine already links each mapped point to its image, but refined points cache only a preiterate. So exact registration is a numerics change, not a tweak.

**Evidence:** `src/tanglepack/numerics/IterateInference.py:119, 151, 436`; `src/tanglepack/numerics/ManifoldMachine.py:372-374, 878-899`; `src/tanglepack/topology/Pseudoneighbor.py` (`_iterate_intersection` `match_rtol 0.05`, `:491`); `StablePartition._backward_endpoint`

### INT-7: Crossing identity is re-derived by tolerance matching on every `compute_intersections`; ids are not reproducible
**Status:** open | **Severity:** medium | **Related:** NUM-8, API-1 | **WP:** WP2

The registry is replaced and `reindex_from` re-matches by `cdist_tol` (absolute 1e-6). `Tangle` keys segment edges on `id()`, so ids differ between builds in one process. `ResonanceZone` re-finds ids by tolerance (`_live_id`, `_anchor_id`), and `reindex_from`'s remap is never consumed. Crossing identity is really combinatorial: segment pair plus ordinal along each ray.

**Evidence:** `src/tanglepack/numerics/IntersectionRegistry.py:241-330, 863`; `src/tanglepack/topology/Trellis.py:24-30`; `docs/dual_graph_symbolic_dynamics_plan.md:324`; `.refactor/knowledge/testing.md` ("Registry ids are NOT reproducible")

### INT-8: A landing span that straddles a cut still feeds a word (midpoint owner, `contained=False`) instead of being unresolved
**Status:** open | **Severity:** medium | **Related:** INT-1, INV-6 | **WP:** WP7

Under a correct partition, a homotopy element's image lies inside one iterated element, so a straddle signals a partition or cdist-estimate problem. The current output is a wrong word with a WARNING. Strict xfail on inversion.

**Evidence:** `src/tanglepack/topology/DualWalk.py:270-380`; `tests/cases.py:812-815`; `CLAUDE.md` ("Landing")

### INT-9: Pairs with no spanning Bridge object drop holes and image bridges; the trellis has no notion of an "unstable arc" independent of Bridge objects
**Status:** open | **Severity:** medium | **Related:** NUM-16, BCG-1 | **WP:** WP7

After a blast there is no Bridge between a parent crossing and a blast-child crossing:
- `punch_holes` skips the pair, and propagation restarts from the most-backward punched pair ("open gap in the punching layer").
- `MinimalTrellis` skips the image pair.
- Virtual classes appear.
- Zero-area regions form (NUM-16).

`BridgeIterator` also creates partial or duplicate Bridge objects over one arc.

**Evidence:** `src/tanglepack/topology/StablePartition.py:190-196, 239-244, 428`; `src/tanglepack/topology/MinimalTrellis.py:280, 323`; `src/tanglepack/topology/SymbolicDynamics.py:494`; `src/tanglepack/numerics/BridgeIterator.py:141-156`; `docs/regions_bridge_identity_and_cleanup_plan.md:196`

### INT-10: Face nesting is decided by geometry (point-in-polygon and signed area), so OuterFaceError fires on slivers
**Status:** open | **Severity:** medium | **Related:** NEST-1, NUM-1 | **WP:** WP10

`_merge_faces`/`_innermost_container` test every node with `point_in_polygon` and signed areas. E11 k28 at 1e-6, blast 16: one component has three positively wound faces, which a planar subdivision cannot have; a self-intersecting sliver is suspected. The rotation system already determines the faces, and the chalkboard needs one containment test per fixed point. There is also a "component with no closed face" WARNING.

**Evidence:** `src/tanglepack/topology/DualGraph.py:597-604, 771, 808-810, 842-875`; `stress/numerics/report/report.tex:1667-1673, 1912`; `stress/numerics/report/numbers.tex:412-416`

---

## Unnecessary complexity / simplification

### SIMP-1: The orientation side-flip rule is coded at six sites and cannot be reached while det<0 is rejected
**Status:** open | **Severity:** low | **Related:** INV-11, DEC-10 | **WP:** WP11

One helper (`image_side(side, n)`) would do, or delete the flips until det<0 is admitted.

**Evidence:** `src/tanglepack/topology/StablePartition.py:563, 596-598, 1287-1300, 1409-1412`; `src/tanglepack/topology/Trellis.py:292, 1423-1424`; `src/tanglepack/topology/DualGraph.py:934-940`; `src/tanglepack/topology/DualWalk.py:317`

### SIMP-2: Two element-image implementations with different tolerance rules
**Status:** open | **Severity:** medium | **Related:** INT-8, SIMP-3 | **WP:** WP11

- `Trellis.image_of_element`: cover plus `_snap_to_partition_boundary` with `SNAP_RTOL`/`SCALING_RTOL`. It feeds `StableNode.images`, which only plots read.
- `land_element`: midpoint owner.

Same rule, possibly two answers.

**Evidence:** `src/tanglepack/topology/Trellis.py:42-105, 128, 143, 1315-1440`; `src/tanglepack/topology/DualWalk.py:270-380`; `src/tanglepack/topology/DualGraph.py:932-956`

### SIMP-3: Cdist scaling is implemented several ways; `scale_cdist` is unused
**Status:** open | **Severity:** low | **Related:** INV-2, INT-6 | **WP:** WP2

- `scale_cdist` divides by the unstable β for stable cdists. Verified 2026-10-07: its only reference is its own docstring.
- `_scaled_image_cdist` uses `per_step_beta(stability)`.
- `_backward_endpoint` divides by β.
- `Pseudoneighbor` uses `lambda_u**m`.
- `land_element` scales the open end by `per_step_beta('stable')`.

`Trellis.lambda_u` and `TrellisBranch.eigenvalue` are extra accessors; the latter has 0 uses and a wrong docstring. `image_cdist` should be the one path.

**Evidence:** `src/tanglepack/topology/Trellis.py:679-763` (:690, :725); `src/tanglepack/topology/StablePartition.py:1744-1768`; `src/tanglepack/topology/Pseudoneighbor.py:491, 551-576`; `src/tanglepack/topology/TrellisBranch.py:72-84`; `src/tanglepack/topology/DualWalk.py:330-334`

### SIMP-4: The 2-D same-point collision test is duplicated in StrongPip and Pseudoneighbor; three branch-cycle wrappers sit over one memo
**Status:** open | **Severity:** low | **Related:** INT-6, SIMP-12 | **WP:** WP11

The collision test exists only for unlinked crossings and would vanish with an exact table. The wrappers are `forward_stable_branch_cycle`, `forward_unstable_branch_cycle` and `_branch_position_map`.

**Evidence:** `src/tanglepack/topology/StrongPip.py:92, 260-280`; `src/tanglepack/topology/Pseudoneighbor.py:114-128, 609-610`

### SIMP-5: PartitionFamily keeps a local copy of `_build_intervals` ("concurrent edit" reason is stale)
**Status:** open | **Severity:** low | **Related:** DOC-3 | **WP:** WP11

**Evidence:** `src/tanglepack/topology/PartitionFamily.py:54-67, 518`; `src/tanglepack/topology/StablePartition.py:1968`

### SIMP-6: Two letter generators (active base-26 in loom, inert u..z,uu.. in topology), with fragile collision-skipping
**Status:** open | **Severity:** low | **Related:** BCG-1 | **WP:** WP6

Active letters can reach `u` after 20 classes, and E11 reaches `bw`/`cz`.

**Evidence:** `src/tanglepack/loom/BridgeAlphabet.py:10-25, 74`; `src/tanglepack/topology/SymbolicDynamics.py:330-345`

### SIMP-7: `Blast.py` tidy-up
**Status:** open | **Severity:** low | **Related:** NEST-4 | **WP:** WP11

- It hand-rolls the fixed-point normalisation that `TangleSession._resolve_fixed_points` already does.
- `cKDTree` is imported inside the loop, with a dead scipy-missing fallback.
- `BlastResult.fixed_point` is typed `Optional[FixedPoint]` but stores lists.

**Evidence:** `src/tanglepack/loom/Blast.py:95-99, 160, 226-236, 281-289, 373`; `src/tanglepack/loom/TangleSession.py:538`

### SIMP-8: Four or five hand-rolled session caches keyed on `(generation, signature[, pips])` for stages that each cost under 0.2 s
**Status:** open | **Severity:** low | **Related:** API-1, API-2 | **WP:** WP8

`trellis(fp)` and `trellis([fp])` are separate cache entries. One topology object per generation (classes → holes → partition → words) would replace them.

**Evidence:** `src/tanglepack/loom/TangleSession.py:291-400, 379-388, 442-560, 603-611, 643-651, 677-860, 729-740, 828-839, 1313`

### SIMP-9: ManifoldInitializer open notes
**Status:** open | **Severity:** low | **Related:** INV-7 | **WP:** WP11

- Ethan vs Kevin construction.
- A duplicate `_insert_point_geometrically` (it also calls the machine's).
- `get_initial_fundamental_segment` is too long.
- Requests: orient at any iterate, or one direction alone.

**Evidence:** `src/tanglepack/numerics/ManifoldInitializer.py:9-45, 391, 557, 588`; `src/tanglepack/numerics/ManifoldMachine.py:902`

### SIMP-10: Two stable trims
**Status:** open | **Severity:** low | **Related:** INT-7 | **WP:** WP11

`TangleWorkbench.trim_stable_manifolds` and `ResonanceZone._trim_stable_at`.

**Evidence:** `src/tanglepack/numerics/TangleWorkbench.py:2045`; `src/tanglepack/loom/ResonanceZone.py:429`

### SIMP-11: The cartoon bridge-drawing block is duplicated between the line and circle layouts; `plotting.py` is 3275 lines
**Status:** open | **Severity:** low | **Related:** — | **WP:** WP11

**Evidence:** `src/tanglepack/topology/plotting.py:2140, 2360-2375, 2922, 3060-3075`

### SIMP-12: Dead or redundant public API, kept and whitelisted by the coverage guard
**Status:** open | **Severity:** low | **Related:** SIMP-3, SIMP-4, PK-7, API-4 | **WP:** WP11

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

**Evidence:** `tests/conftest.py:85-93`; `src/tanglepack/numerics/IterateTable.py:81, 98, 107`; `src/tanglepack/numerics/IntersectionRegistry.py:792, 800`; `src/tanglepack/topology/SymbolicDynamics.py:167, 787`; `src/tanglepack/topology/TrellisBranch.py:72, 96`; `src/tanglepack/topology/BridgeClass.py:410`; `src/tanglepack/loom/Blast.py:178`; `src/tanglepack/topology/Trellis.py:24-30, 690`; `docs/dual_graph_symbolic_dynamics_plan.md:326`

### SIMP-13: Tolerance and warning density
**Status:** open | **Severity:** low | **Related:** INT-6, INT-7 | **WP:** WP11

- 86 tolerance sites in `topology/`.
- 73 WARNINGs against 38 asserts across the library.
- About 130 lines of `tangle_of`/`min_cdist`/`anchor_outward_key` ordering just to make letters deterministic.

**Evidence:** chalkboard §8(b); `src/tanglepack/topology/BridgeClass.py:465-560, 659-693`

---

## Resolution & numerics

### NUM-1: The bridge-row AssertionError fires on under-resolved bridges even when its own diagnosis says the combinatorial rows agree
**Status:** deferred (WP7-D1: ≥3 nodes per bridge) | **Severity:** medium | **Related:** INT-3, NEST-1, DEC-18 | **WP:** WP4

Seen in 39 snapshots across 4 configs: p3 and nested_inner at 1e-5 from blast 6, and the 2-iteration cells. Example: p3 bridge (41,42) has signs +1/−1 and the combinatorial rows agree (both left), but the geometric row is read on a 2-node polyline at 0.16°. The heteroclinic case at 1e-5 asserts with "signs alternate … likely below resolution" (8 nodes, 0.3-0.6°). Asserting here contradicts the 2026-10-07 rule that resolution problems are flags. Relaxing only the geometric reading exposes the hole-side check next.

**Evidence:** `src/tanglepack/topology/StablePartition.py:1269`; `stress/numerics/report/report.tex:235-237, 1674-1690, 1728-1729, 1910-1912, 1995-2002`; `stress/numerics/report/numbers.tex:424-428`

### NUM-2: Hole-side AssertionError: carried "left" vs measured "right"
**Status:** open | **Severity:** medium | **Related:** INT-3, NUM-1 | **WP:** WP4

3 snapshots: k28 at 1e-5/1e-6, blasts 16-17, several iterates of one orbit. Probably measurement on under-resolved sub-arcs; the recommendation is to re-measure on refined sub-arcs.

**Evidence:** `stress/numerics/report/report.tex:240, 1654-1666, 1730-1732, 1910`; `stress/numerics/report/numbers.tex:418-422`

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
**Status:** open | **Severity:** low | **Related:** INT-6, INT-7 | **WP:** WP2

68% of consecutive k28 gaps fall below it by blast 14, and it caused E11 error 1. The fix is a relative tolerance. The near-ULP vs exact merge dedup question is also undecided.

**Evidence:** `src/tanglepack/numerics/IntersectionRegistry.py:85-95`; `src/tanglepack/topology/StablePartition.py:200-206`; `stress/numerics/report/report.tex:1645-1649, 2003-2005`

### NUM-9: `ManifoldMachine.merge_manifolds` crashes on non-contiguous `old_iterated_points` (None head when a bridge reaches the grown tip)
**Status:** open | **Severity:** medium | **Related:** NUM-10, NEST-4 | **WP:** WP10

`build_nested` no longer reaches it, but the bug remains.

**Evidence:** `src/tanglepack/numerics/ManifoldMachine.py:391-404, 516`; `.refactor/knowledge/topology.md`; `.refactor/decisions.md` 2026-10-02; memory `loom-and-subpackages:76-77`

### NUM-10: Blast silently counts AttributeError/ValueError iterate failures as "skipped"
**Status:** open | **Severity:** medium | **Related:** NUM-9 | **WP:** WP10

About 5 failures per blast in Phase 6 ("'NoneType' object has no attribute 'cdist'"). The count varies by run and was never examined, so it can hide NUM-9.

**Evidence:** `docs/regions_bridge_identity_and_cleanup_plan.md:339`; `.refactor/knowledge/loom.md`

### NUM-11: Escape to infinity: `grow_n_times` has no notion of escape, and PointBudgetExceeded leaves a session that cannot be resumed
**Status:** open | **Severity:** medium | **Related:** INV-4, TEST-1 | **WP:** WP10

One step multiplies points by up to 40× once |x| > ~100. p3 at 18 steps and inversion at 7-8 steps hit the budget. Recipes and tests still use fixed step counts instead of `grow_until_bounded`.

**Evidence:** `stress/numerics/report/report.tex:404-452, 899-934, 2003-2004`; `stress/numerics/report/numbers_extra.tex:63-131, 259-265`

### NUM-12: Growth after a blast recuts the bridge set and clears the `iterated` flags (blast genealogy lost)
**Status:** open | **Severity:** low | **Related:** BCG-3 | **WP:** WP10

**Evidence:** `src/tanglepack/numerics/TangleWorkbench.py:1071-1080`

### NUM-13: `grow_until_arclength` stops on canonical distance, not arc length
**Status:** open | **Severity:** low | **Related:** — | **WP:** WP10

**Evidence:** `src/tanglepack/numerics/TangleWorkbench.py:21-26, 823`

### NUM-14: `iterate_bridge` does not re-assert Bridge invariants
**Status:** open | **Severity:** low | **Related:** INT-9 | **WP:** WP10

**Evidence:** `src/tanglepack/numerics/ManifoldMachine.py:80-90`

### NUM-15: Small-k recipes: MaxIterationsReached for k ≤ 0.545, and 26/50 b=1 saddles give trivial tangles (fixed 8-step recipe)
**Status:** open | **Severity:** low | **Related:** NUM-11 | **WP:** WP10

Recipe depth should scale with 1/ln λ_u. These are first-campaign numbers.

**Evidence:** `stress/numerics/report/report.tex:823-826, 864-871, 883-889`; `stress/numerics/report/numbers_extra.tex:199-213`

### NUM-16: Arrangement: a zero-area region after a blast reports a representative point inside its neighbour
**Status:** open (strict xfail on k28_one_blast and k28_two_blasts) | **Severity:** medium | **Related:** INT-9 | **WP:** WP10

It is bounded by an unstable arc that no Bridge spans; this is a sparse/full arrangement mismatch at a blast child.

**Evidence:** `tests/cases.py:816-825`; `docs/regions_bridge_identity_and_cleanup_plan.md:196`

### NUM-17: Arrangement `image_of` misses 39 of 86 regions; union-valued images are undefined
**Status:** decided 2026-10-07 | **Severity:** low | **Related:** DEC-17 | **WP:** deferred

None found for k10, unblasted k28 or inversion (first-campaign numbers). `arrangement_image_of`/`preimage` are N/A on k10 and inversion. The hole set is not forward invariant at the computed extent.

**Evidence:** `tests/cases.py:874-880`; `stress/numerics/report/report.tex:205, 722-728, 2023-2027`; `stress/numerics/report/numbers.tex:190-192`; `docs/regions_bridge_identity_and_cleanup_plan.md:338, 342`

### NUM-18: The pseudoneighbor stable-side check walks only registered crossings (asymmetric with the unstable all-iterates check)
**Status:** deferred (author chose no change) | **Severity:** low | **Related:** BCG-5 | **WP:** deferred

**Evidence:** memory `pseudoneighbor-blast-semantics:8`

### NUM-19: Small deferred plan findings
**Status:** possibly stale — verify | **Severity:** low | **Related:** — | **WP:** WP0 → owners (verify in the WP touching each file)

- `_insert_crossing_separator` stale edge bounds.
- `_validate_saddle` may reject near-parabolic period-k saddles.
- `bridge_classes` silently accepts a reversed BridgeId.
- The walk and the class decide ownership by different rules.

**Evidence:** `docs/regions_bridge_identity_and_cleanup_plan.md:319-320, 335`; `docs/dual_graph_symbolic_dynamics_plan.md:328, 332`

---

## Performance

### PERF-1: The blast-end same-stability diagnostic does not scale past ~10 blasts
**Status:** decided 2026-10-07 | **Severity:** high | **Related:** DEC-12 | **WP:** WP10

`blast_zone(diagnose=True)` is the default and materialises every u×u hit:
- k28: 104,760 hits (1.58× the blast time)
- k32: 170,192 (1.76×)
- nested_outer: 589,693 (1.96×)

It cut reachable depth by 3-4 blasts. The hits are near-parallel image arcs at about 1e-4. Options: count instead of storing, check only new image bridges, or default `diagnose=False` when blasting.

**Evidence:** `src/tanglepack/numerics/TangleWorkbench.py:1495-1525`; `src/tanglepack/numerics/diagnostics.py`; `stress/numerics/report/report.tex:232-234, 1532-1564, 1594-1603, 1986-1997`; `.refactor/knowledge/testing.md` (OPEN)

### PERF-2: Blast cost grows geometrically, about 1.6-1.7× per blast in period-1 zones
**Status:** open | **Severity:** medium | **Related:** BCG-1, TEST-7 | **WP:** WP12

k28: ×1.61 crossings, ×1.65 points, ×1.71 time; the last blast takes 82 s and 9.0 GB, and memory stops it at 15 blasts. Coarser cutoffs hit the 780 s budget at 16-17. Blasting iterates every interior bridge, including those outside the fundamental domain (the same bridges that inflate the class table).

**Evidence:** `stress/numerics/report/report.tex:207, 1307-1364`; `stress/numerics/report/numbers.tex:458-466, 566-574, 776-784`

### PERF-3: Growth cost is dominated by the Python object graph (one Point per vertex); GPU gives no gain
**Status:** deferred | **Severity:** low | **Related:** — | **WP:** deferred

The map is at most 0.84% of growth time; refinement plus the linked list take 75-83%, at about 1 KB per point. GPU runs 0.96-1.00×.

**Evidence:** `stress/numerics/report/report.tex:103-108, 241-243, 523-565`

### PERF-4: The hole-propagation `_containing_bridge` is a linear scan over all bridges per backward step; termination uses `id(bridge)`
**Status:** open | **Severity:** low | **Related:** INT-2 | **WP:** WP4

**Evidence:** `src/tanglepack/topology/StablePartition.py:571-588, 1771-1800`

### PERF-5: Harness minors
**Status:** open | **Severity:** low | **Related:** PERF-1 | **WP:** WP12

- The E7 harness runs the same-stability pass twice.
- E11's cheap-path "unresolved" count is registry-wide.

**Evidence:** memory `stress-recs-plan-2026-10-07` (harness minors)

---

## Tests & coverage gaps

### TEST-1: No law case combines inversion with period-k, no nested case blasts its inner zone, and there is no heteroclinic or det<0 case
**Status:** open | **Severity:** medium | **Related:** INV-4, INV-11, NEST-1, NEST-3 | **WP:** WP1

`LAW_CASES` are k10, k28_one_blast, k28_two_blasts, p3 (13 steps, unblasted), nested (2 outer, 0 inner) and inversion (no holes). Inner blasting goes untested even though E11 nested_inner grows 13→55. There is no inversion or det<0 case in E11 or in the deep sweep.

**Evidence:** `tests/cases.py:385-388, 493-501, 745-752`

### TEST-2: No law asserts structural stability under blasting
**Status:** open | **Severity:** medium | **Related:** BCG-1, BCG-4, BCG-7 | **WP:** WP1

Nothing checks that "more blasts refine, never contradict": a stable hole set between switches, class count independent of blast count between switches, or letters kept. The p3 burst and the linear growth pass a green suite (833 pass). `laws.py` reads the class count but never bounds it. Under the no-golden policy and the "topology change is expected" rule, this has to be a rule-based law.

**Evidence:** `tests/helpers/laws.py:1974`; `tests/regression/test_blast_monotonicity.py`; `tests/conftest.py:61-66`; memory `topology-change-is-expected`

### TEST-3: NOT_APPLICABLE pairs are dropped inside the layer tests, so pytest reports 0 skips
**Status:** open | **Severity:** low | **Related:** INV-4 | **WP:** WP1

The 15+ laws switched off on inversion, and the k10/p3/k28 exclusions, never appear in the summary. Last run: 833 passed, 3 deselected, 11 xfailed, 0 skipped, 100.88 s.

**Evidence:** `tests/helpers/law_tier.py:186-202`; `tests/cases.py:849-889`

### TEST-4: Provisional, untested rules
**Status:** provisional | **Severity:** low | **Related:** API-5 | **WP:** WP1

- The refined-point cdist is the mean of its neighbours.
- Anchor faces never close (the `grow_until_faces_closed` cap).

**Evidence:** `tests/conftest.py:68-83`; `docs/regions_bridge_identity_and_cleanup_plan.md:250`

### TEST-5: Test-policy open items
**Status:** open | **Severity:** low | **Related:** NUM-8 | **WP:** WP1

- `generation_bump` cache cells are blind on k10.
- The p3 near-ULP regression uses a seed-seam stand-in.
- A third tangle has no `TANGLE_COLOR_FAMILIES` base.
- The colour-family-per-tangle rule test is still to do.

**Evidence:** memory `test-suite-policy-2026-10-05:52`

### TEST-6: `coverage_guard.py` trips when lines shift
**Status:** open | **Severity:** low | **Related:** — | **WP:** WP11

**Evidence:** `tests/_tools/coverage_guard.py`; memory `stress-recs-plan-2026-10-07`

### TEST-7: `blast_zone(max_points=...)` (WP5) is never exercised by the harness
**Status:** open | **Severity:** low | **Related:** DEC-13, PERF-2 | **WP:** WP10

**Evidence:** `stress/numerics/report/report.tex:1970`; `.refactor/knowledge/testing.md`

---

## API traps

### API-1: The pip choice is lost on every generation bump, so `_repin` after each blast is dead work
**Status:** open | **Severity:** medium | **Related:** SIMP-8, NUM-6 | **WP:** WP8

A rebuilt trellis starts with empty strong-pip slots. The session should carry the choice, as an Intersection, across generations.

**Evidence:** `src/tanglepack/examples/henon_cases.py:142-153, 192, 199, 260`; `.refactor/knowledge/loom.md`; `.refactor/decisions.md` 2026-10-02; memory `strong-pip-choice-gotchas`

### API-2: The session gathers partitions from whichever per-fixed-point trellises happen to be cached (order-dependent)
**Status:** open | **Severity:** medium | **Related:** SIMP-8, NEST-1 | **WP:** WP8

`_gathered_partitions` scans every non-stale cached trellis without filtering by the requested fixed points ("first built wins"). The holes and pips gatherers filter only single-fixed-point trellises. Bridges come from the all-fixed-point trellis, so `trellis()` and `trellis(fp)` give different results.

**Evidence:** `src/tanglepack/loom/TangleSession.py:479-505, 528-545, 755-779`

### API-3: `PointBudgetExceeded.blast_context` is set dynamically and not declared; `MaxIterationsReached` and `TURNAROUND_LOG_STRETCH` are not exported
**Status:** open (verified 2026-10-07) | **Severity:** low | **Related:** — | **WP:** WP11

**Evidence:** `src/tanglepack/loom/Blast.py:262, 388`; `src/tanglepack/numerics/budget.py` (no `blast_context`); `src/tanglepack/__init__.py:30, 78`; `src/tanglepack/numerics/__init__.py:27, 56`

### API-4: The deprecated no-op `invalidate_trellises` is still present
**Status:** open | **Severity:** low | **Related:** SIMP-12 | **WP:** WP11

**Evidence:** `src/tanglepack/loom/TangleSession.py:1290-1306`

### API-5: `grow_until_faces_closed` never converges (it returns 0 immediately or hits the cap); anchors never close under the four-anchor model
**Status:** open | **Severity:** low | **Related:** INV-1, TEST-4 | **WP:** WP3

**Evidence:** `src/tanglepack/numerics/TangleWorkbench.py:1204-1265`

### API-6: Coordinate shape is (2,1) in FixedPoint and (2,) in BasePoint, so callers ravel everywhere
**Status:** open | **Severity:** low | **Related:** — | **WP:** WP11

Also: `partial_jacobians` is unexplained, and BranchPoint has a `_description_` placeholder.

**Evidence:** `src/tanglepack/numerics/FixedPoint.py:10-21`; `src/tanglepack/numerics/BasePoint.py:9-11`; `src/tanglepack/numerics/BranchPoint.py:8-12, 41`; `src/tanglepack/numerics/Point.py:8-10`

### API-7: `build_nested` always returns the slug `nested_p1_p3`, so figures overwrite each other whatever the blast counts
**Status:** open | **Severity:** low | **Related:** DOC-4 | **WP:** WP11

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

- MinimalTrellis says the inertness definition is "changing" (settled 2026-09-16).
- TrellisBranch's `return_period` question is resolved.
- ResonanceZone says "no inversion example" (`build_inversion` exists).
- PartitionFamily cites a "concurrent edit".
- TopologyResults asks whether a Hole should be a coordinate.
- ManifoldMachine has an old TODO.
- StrongPip/Pseudoneighbor say the key is "legitimately None".
- The topology-layer and loom memories say inversion is "not validated (no example)".

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

Every question is now answered or deferred; see "Author decisions" below.

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
| WP9 | Nested and heteroclinic: A×B fixture, heteroclinic holes in B, innermost-zone blasting, per-fixed-point `rebuild_bridges` | NEST-1, 2, 4..6 |
| WP10 | Numerics and performance leftovers | PERF-1, NUM-9..16, INT-10, TEST-7 |
| WP11 | Simplification sweep and API traps | INV-11, SIMP-*, API-3, 4, 6, 7, TEST-6 |
| WP12 | Rerun E7/E11 and the figures, confirm class counts plateau, update CLAUDE.md / knowledge / memory, close issues | NUM-5..7, PERF-2, PERF-5, DOC-1..5 |
| deferred | — | NUM-3, NUM-4, NUM-17, NUM-18, PERF-3 |

Verification of the "possibly stale" items (WP0): DOC-3 is confirmed open (the MinimalTrellis "definition is changing", TrellisBranch `return_period`, ResonanceZone "no inversion example", PartitionFamily "concurrent edit" and StrongPip "legitimately None" notes are all still present). PK-7 is confirmed dead code (WP2 deletes it). INV-10, NUM-19 and DOC-4 need a run; they are checked by the WP that owns them.
