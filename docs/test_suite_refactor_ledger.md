# Test-suite refactor ledger (2026-10-05)

The running record of the approved test-suite refactor
(`.refactor/runs/2026-10-05-test-suite/plan.md`, binding "Author decisions";
detail in `planner_detail.md`). Every phase appends its section: what it
deleted (node id → reason → regression guard checked), what moved where, and
any deviation from the plan. The baseline node-id list is
`.refactor/runs/2026-10-05-test-suite/nodeids_base.txt`; the final
verification diffs every surviving and deleted node id against it.

---

## §0 Phase 0: baseline and analyst probes (no test edits)

### Baseline

| Item | Value |
|---|---|
| Branch / HEAD at baseline | `test-suite-refactor` @ `02b11ef` |
| Collected | **796** (795 passed, 1 skipped) |
| Skip | `tests/numerics/test_gpu.py::test_enable_gpu_without_cupy_raises_clearly` ("CuPy is installed; the no-CuPy error path is not exercised") — Phase 9 makes it run via `sys.modules` monkeypatch |
| Wall | 84 s with `--cov` (45 s without) |
| Total coverage (`src/tanglepack`) | 90 % |
| Slowest | 3.5 s `test_session_bridge_classes::test_active_class_lies_in_the_zone_after_trimming_at_the_image_pip`; 3.1 s `test_higher_period_cartoon::test_nested_blast_order_does_not_matter`; 2.5 s `test_grow_until::test_faces_closed_raises_at_the_cap_for_the_anchor` |

Run-folder artifacts (`.refactor/runs/2026-10-05-test-suite/`):
`nodeids_base.txt`, `cov_base.json`, `junit_base.xml`, `durations_base.txt`,
`baseline_run.txt`, `cov_p0.json` (Phase 0 verification run, identical
coverage), the audit reports (`synthesis.md`, `critique.md`, `infra.md`,
`durations.txt`, `catalog_journal.jsonl`), `plan.md`, `planner_detail.md`, and
`probes/` (one script and one `*_output.txt` per probe).

### Coverage guard

`tests/_tools/coverage_guard.py BASE.json CANDIDATE.json` compares the
coverage.py JSON reports module by module and exits 1 when a line the baseline
executed is missed in the candidate, unless the line lies inside a whitelisted
dead-API qualname (spans resolved with `ast` on the current source, so they
survive line shifts; a stale entry is reported, not failed). Whitelist (plan
decision 4): `IntersectionRegistry.on_interval`,
`IntersectionRegistry._get_lambda_u`, `ManifoldMachine._curvature_area`
(scalar), `TangleSession.invalidate_trellises`,
`Pseudoneighbor.forward_unstable_branch_cycle`,
`StrongPip.forward_stable_branch_cycle` (`BoundaryArc` no longer exists in
`src/`, so it needs no entry). Self-tested: a non-whitelisted line removed from
a candidate report fails the guard; lines inside `on_interval` pass as
whitelisted. `cov_base.json` vs `cov_p0.json`: OK.

### P1: the hole quirk (decision 11). Verdict: rows really flip, but it is NOT a bug; the openings law replaces the snapshots, no xfail

Probe: `probes/p1_hole_quirk.py` (output `p1_output.txt`), on
`build_period3()` and the k=2.8 two-blast build. It records, per hole, the
`bounding_ids`, near/far bound, whether the near bound is an anchor, the
openings, the row of the hole's bridge at each bound, the hole's own geometric
side of the nearest stable segment, and for each bound without an opening
whether it is an anchor or a trimmed tail. It also checks which bounds are
registered iterates of the origin's bounds.

Findings:

- **Rows really flip.** In the k=2.8 two-blast build the origin `(9,10)` (iterate
  0, side `left`) opens the LEFT pair, but its backward holes at -1, -2, -3
  open the RIGHT row. In p3 the `right`-side orbit's direct holes (iterates
  0..2) open the LEFT pair, but its backward holes (-1..-3) open the RIGHT
  pair. The `left`-side p3 orbit keeps its row (RIGHT) and loses only the
  anchorward opening at the anchor.
- **The flip is geometrically correct.** Every flipped hole sits in a
  containing bridge whose bounds are NOT registered iterates of the origin's
  bounds. The origin's preimages lie beyond the trimmed stable manifold:
  `T.iterate(9, -1)` and `T.iterate(10, -1)` are None for k=2.8, and likewise
  for p3's right orbit. The backward image lobe is therefore bounded by a
  stable segment that was never computed. Its hole lands in a larger bridge
  that meets the computed stable manifold from the other side. In every
  flipped case the recorded row equals that bridge's `_row_at` row and the
  hole's own geometric side of the nearest computed stable segment. For
  example, the k=2.8 iterate -1 hole is on the right, 5.1 away from the
  nearest segment. An orientation-preserving map preserves the row of a TRUE
  image lobe. Where such a lobe exists, the row is preserved.
- **Linked-bound law** (`--law` sweep, all five non-inversion cases): when a
  propagated hole's bound IS the registered iterate of its origin's direct-hole
  bound, the opening there has the origin's `(which, row)`. Results: k10
  checked 0 (4 unlinked bounds); k28 one blast 1/0 violations; k28 two blasts
  0 checked (6 unlinked); p3 3/0; nested 3/0. No violations anywhere.
- **Every other difference** is an opening dropped at an anchor bound, which
  is the documented `_hole_openings` rule ("at the anchor artifact a hole
  facing the far half … opens nothing"). The baseline shows no tail drops.

Applying the decision rule:
- Literally, (a) fails: not every difference is a drop at an anchor or tail.
- (b)'s condition "rows really flip" holds.
- (b)'s premise "then it is a bug" does not hold: the flips are the correct
  geometry of a hole whose image lobe's stable side was never computed.

Recording `test_propagated_holes_open_their_origins_pair` as
`xfail(strict=True)` would encode a false law as a known bug, so **no xfail is
added**. The closest action consistent with the decisions:

- **Phase 4 openings law** (in `test_law_partition`), three parts:
  1. a direct hole opens its inward pair on its own bridge's row (already a
     planned law);
  2. linked-bound law: at a bound that is the registered iterate of the
     origin's bound, a propagated hole keeps the origin's `(which, row)`;
  3. an opening is missing only at an anchor bound (stable cdist 0) or a
     trimmed tail.
- `NOT_APPLICABLE` for the linked-bound part where it checks 0 items: k10 and
  k28 two blasts (and inversion, which has no holes).
- The two snapshots `test_p3_hole_sides_are_pinned` and
  `test_k28_two_blast_hole_sides_are_pinned` are deleted in the Phase 4 commit
  that adds the law.

**AUTHOR ITEM (E7):** confirm that a backward hole's row may differ from its
origin's when the image lobe's stable side lies beyond the trimmed manifold.
If the author instead wants the row of the uncomputed true lobe, that is a
source change (out of scope), and the linked-bound law still holds.

### P2: the k=10 inversion saddle through the full pipeline (decision 9)

Probe: `probes/p2_inversion.py` (output `p2_output.txt`) and
`probes/p2b_inversion_depth.py` (`p2b_output.txt`). The saddle is at
`(-2.3166, 2.3166)`: period 1, `k_value` 2, λu = -4.406, λs = -0.227
(det J = +1, orientation-preserving; NOT orientation-reversing). Recipe: 6
unstable and 5 stable steps (the `henon_inversion` fixture), then the full
session pipeline.

**Stages:** none raises. initialize, grow, compute_intersections, trim,
create_bridges, infer_iterate_table, classify_strong_pips,
compute_pseudoneighbors, punch_holes, partition_stable_manifold, arrangement,
bridge_classes, minimal_trellis, iterated_partition, dual_graph and
symbolic_dynamics all complete. The whole run takes about 1.3 s, almost all of
it growth and intersection.

**Counts:** 11 crossings, 7 bridges, **0 pseudoneighbors, 0 holes**, 4
partitions, 4 bounded regions, 4 bridge classes (2 active, 2 inert), and 1
virtual symbol `new1`. `is_reliable` is False.

**Anchors:** 4 in total, one per (unstable branch, stable branch) pair, which
gives **2 per unstable branch** (`(0,0)`: 2, `(0,1)`: 2). The author's rule is
exactly 1, so this is a KNOWN ISSUE.

**Law results on inversion:**

| Law | Result |
|---|---|
| one anchor per unstable branch | **FAIL** (2 per branch) → `xfail(strict=True)` |
| manifold cdist monotone (non-strict), and strictly | PASS (4 manifolds) |
| no geometric spikes / iterate law / one-to-one | PASS (4) |
| every crossing is unstable × stable | PASS (11) |
| area product preserved along iterate chains | PASS (3 chains) |
| BridgeId unstable order; partial ⇔ id None | PASS (7) |
| bridge endpoints on the bridge's own branch | PASS (14 ends) |
| `bridges_at` exact | PASS (11) |
| single copy per BridgeId | PASS (7) |
| I2 bridge rows consistent | PASS (7) |
| partition one owner per crossing per side | PASS (22) |
| I1 holes share bridge side; no direct hole at iterate ≥ k_value | VACUOUS (0 holes) → `NOT_APPLICABLE` |
| arrangement `image_of` preserves area | VACUOUS (`image_of` returned None for all 4 regions) → `NOT_APPLICABLE` |
| every bridge in exactly one class | PASS (7) |
| itineraries even | PASS (4) |
| no virtual symbols / `is_reliable` | **FAIL** (`new1`) → not a law; inversion is not in golden |

**Depth sweep:** 7 unstable steps does not finish within 240 s (refinement
blow-up), so 6 is the practical depth and the hole/partition laws stay vacuous
on inversion.

**Draft `KNOWN_ISSUES`:**
- `("one_anchor_per_unstable_branch", "inversion")`: 2 per branch.
- `("orientation_reversing_case_builds", "orientation_reversing")`:
  ValueError (P4).
- `("open_p3_deep_runs", "p3_15|p3_16|p3_6_blasts")`: P3.

**Draft `NOT_APPLICABLE`:**
- `(hole/partition-hole laws, "inversion")`: no pseudoneighbors at the only
  feasible depth.
- `("arrangement_image_of", "inversion")`: no closed images.
- `("openings_linked_bound", "k10"|"k28_two_blasts"|"inversion")`: see P1.

### P3: the deep p3 runs (decision 10)

Probe: `probes/p3_deep_p3.py` (`p3_output.txt`; full rules in
`p3_output_full.txt`).

| Run | Holes | I1 | Classes (inert) | Virtual | Unresolved | `is_reliable` |
|---|---|---|---|---|---|---|
| default (13 steps) | 12 | holds | 9 (6) | none | none | True |
| 4 blasts / 5 blasts | 12 | holds | 9 (6) | none | none | True |
| **15 steps** | 42 | holds | 23 (1) | `new1` | `A:p3@2.0/R#2 <-> A:p3@2.0/R#3`: "dual-graph walk unreachable" | False |
| **16 steps** | 46 | holds | 30 (6) | `new1` | the same class, unreachable | False |
| **6 blasts** | 42 | holds | 23 (1) | `new1` | the same class, unreachable | False |

Exact assertions for `regression/test_open_p3_deep_runs.py`:
- I1 (`check_holes_share_bridge_side`) holds. This is a plain passing test.
- `xfail(strict=True)`: `not dyn.virtual_classes`, and no
  `cd.unresolved_reason`, for each of the 15-step, 16-step and 6-blast runs.

Each build takes ≤ 0.3 s.

### P4: the b = -1 placeholder (decision 9)

Probe: `probes/p4_orientation_reversing.py` (`p4_output.txt`). For b = -1 the
fixed points are x = y = ±√k, and the Jacobian `[[2x,1],[1,0]]` always gives a
real saddle with det = -1 (exactly one negative eigenvalue).

**Chosen: k = 4, b = -1, seed `[2, 2]`** (eigenvalues 4.236 and -0.236).
`construct_fixed_point` raises **ValueError from `FixedPoint.set_k_value`**,
as do `[-2, -2]` and k = 10 (`[±3.1623, ±3.1623]`). `saddle_guesses(4, -1)`
raises KeyError, so `build_orientation_reversing` passes the seed explicitly.
This is the placeholder `xfail(strict=True, raises=ValueError)`.

### P5: regression-guard audit (decision 4 / rules)

Probe: `probes/p5_guard_audit.py` (`p5_output.txt`). It greps the nine
guard notes for every Phase 2 and Phase 7 deletion candidate's name and
subject keywords. Every hit was read. **Verdict: the planner's Guard column
("none" for the Phase 2 deletions) stands.** Some candidates touch a real past
bug, but each of those bugs keeps a guard elsewhere:

| Candidate | Bug in memory | Kept guard |
|---|---|---|
| `test_per_step_beta_matches_the_expression_it_replaces` | k_value-root `per_step_beta` (regions-refactor) | `test_per_step_beta_is_the_k_th_root_only_without_inversion` and inversion `test_one_map_step_scales_cdist_by_per_step_beta` (KEEP both) |
| `test_fixed_point_takes_no_branch_count` | `has_inversion()` crash at num_branches = 2 (codebase-audit) | `test_kevin_way_builds_both_branches_in_opposite_directions` / two-branch init (KEEP → inversion case sanity) |
| Phase 7: mixed-sign rejection ×2 | `set_k_value` raises on mixed signs (regions-refactor) | the b = -1 placeholder (P4), which must land in the SAME commit as the deletion |
| Phase 7: `_region_key` test | singleton bridge, holes on either side (pseudoneighbor-partition) | rewrite as "both holes of a singleton bridge survive propagation" (planned) |
| Phase 7: `test_sparse_kept_endpoints_are_degree_three…` | inversion arrangement has four coincident anchors (regions "least confident") | none needed; the anchor xfail (P2) covers the inversion anchor issue |
| `_stable_arc_midpoint` ×2, `…hug_the_stable_manifold` | hole display placement (pseudoneighbor-partition) | display only, not a bug; coordinates are still covered by `test_direct_hole_side_is_the_side_of_its_coordinates` and `test_henon_propagated_hole_side_matches_coords` (KEEP) |
| `test_forward_pairs_beyond…[1..k-1]` | provisional exemption (holes-backward-only) | the firm iterate ≥ k_value parameters stay |
| `visualizations/`, `test_workbench_split.py`, the tombstones, `_get_lambda_u`/`on_interval`, `invalidate_trellises`, `_curvature_area` | none (layout, dead API or deprecation) | – |

### P6: builder parity (decision 7, risk R2)

Probe: `probes/p6_builder_parity.py` (`p6_output.txt`). It rewrites
`tests/conftest.py` temporarily and restores it in a `finally`; `git diff`
confirms the conftest is unchanged.

| Variant | Result | Consequence for Phase 1 |
|---|---|---|
| **A**: `henon_tangle_with_bridges` = 9 unstable steps + turnaround + infer, instead of 7 + turnaround + 2 without infer | 795 passed, **0 failing** | `build_k10(unstable_steps=9)` with infer can back every k10 fixture; no 7+2 parameter is needed |
| **B4**: `henon_p3_session` = `build_nested(outer_blasts=0, p1_area_cutoff=1e-4)` (adds `_repin` and partition) | 1 failing: `test_partition_elements::test_image_of_element_requires_the_image_branchs_partition` | That test needs an UNPARTITIONED trellis. `build_nested` needs a `through` stage before `"partitioned"` (or the test builds its own). `_repin` itself changes nothing |
| **B7**: as B4 with `p1_area_cutoff=1e-7` | 4 failing: the B4 test; `test_bridge_cutting_pin_period_3` and `test_p3_zone_areas_and_containment_are_unchanged` (both Phase 2 deletions); `test_dual_graph::test_k10_a_different_pip_moves_the_unified_set` (order-dependent, see below) | The p1 cutoff changes the outer tangle, so keep `p1_area_cutoff` explicit (E5: 1e-4 for the unblasted fixture) |
| **C**: k=2.8 fixtures function-scoped | 795 passed, 0 failing; suite +2 s | Function scope is safe for the k=2.8 cases |

**Finding (affects the isolation spot-checks):**
`test_k10_a_different_pip_moves_the_unified_set` fails when run alone, as
already known. The root cause is that **registry ids are not deterministic
across k10 builds in one process.** Three consecutive `build_k10()` give
strong-pip candidates `[1,2,3,4,7]`, `[2,3,4,6,7]` and `[1,3,5,6,7]`, with the
same crossings under permuted ids. The test picks `alternatives[0]` (the
lowest id), so whether it gets the candidate whose k-th iterate is
unregistered (stable cdist 8.707, no `+1` link) depends on id assignment.

- The regions memory note says "p3 only"; it is k10 too.
- Fix in Phase 1/7: choose the alternative by stable cdist among the
  candidates whose `T.iterate(c, k_value)` exists, never by id.
- Until then, this node id is expected to fail in the isolation checks.

### P7: k=2.8 two-blast cut facts → id-free relations

Probe: `probes/p7_k28_cut.py` (`p7_output.txt`), three builds in one process.
The CLAUDE.md numbers are registry ids from one build. They do NOT reproduce:
today `R_1^1 = [0, 11]` and `R_1^2 = (11, 8)`, and the ids permute between
builds (build 3's `L_1^1.hi` is 3, not 2). The structure is identical
throughout.

**Today's mapping** (CLAUDE.md id → this build's id): 10 → 11, 7 → 8, 6 → 7,
3 → 2, 8 → 9, 9 → 10, 4 → 1. The hole `(6,5)` is today's iterate -1 hole
`(7,6)`. The hole `(1,3)` is the iterate -2 hole.

All 18 relations below hold on all three builds (54/54). These are the golden
k28 relations:
- R_1^1 begins at the anchor (lo cdist 0, closed). R_1^1.hi = R_1^2.lo and
  R_1^2.hi = R_1^3.lo.
- Brackets: R_1^1 `[ ]`, R_1^2 `( )`, R_1^3 `[ ]`; R_5^1 `[ ]`, R_5^2 `( )`,
  R_5^3 `[ ]`; L_1^1 `[ )`, L_1^2 `[ ]`. R_5^1.hi = R_5^2.lo and
  R_5^2.hi = R_5^3.lo.
- L_1^2 has exactly R_5^1's two ends. The L_1^1/L_1^2 boundary = R_5^1.lo,
  and it is a bound of the iterate -2 hole (the anchor lobe's fold abuts that
  hole).
- Holes sit at iterates exactly {0, -1, -2, -3}. R_2's ends are the iterate -1
  hole's bounds.
- **Lobe base:** R_1^2's ends are the `+1` images of the iterate -1 hole's
  bounds.
- **Chord:** R_5^2's ends (= L_2's ends) are the iterate 0 (reference) hole's
  bounds.

### P8: runtime model (E1, R1)

Probe: `probes/p8_runtime.py` (`p8_output.txt`). Median of 3 full builds,
including every downstream product:

| Case | Build | Downstream | Total |
|---|---|---|---|
| k10 | 0.05 s | 0.00 s | 0.05 s |
| k28 one blast | 0.23 s | 0.00 s | 0.23 s |
| k28 two blasts | 0.28 s | 0.01 s | 0.29 s |
| p3 | 0.07 s | 0.01 s | 0.07 s |
| nested | 0.41 s | 0.04 s | 0.45 s |
| inversion | 1.25 s | 0.00 s | 1.25 s |
| **all six** | | | **2.34 s** |

With one module-scoped read-only build per law module (decision 7), the ~11
law modules cost about 26 s of builds, inside the ≤ 40 s budget. Per-check
fresh builds (about 330 × 0.4 s average) would not fit, which confirms the
decision 7 exception. Inversion dominates (54 %), almost entirely growth and
intersection.

### Phase 0 verification

- `tests/` is unchanged except for the new, uncollected `tests/_tools/coverage_guard.py`.
- The collected node-id list is identical to `nodeids_base.txt` (796).
- Suite: 795 passed, 1 skipped in 84 s; coverage guard OK against `cov_p0.json`.
- `-rxX` shows no xfail or xpass.
- Isolation, 5 node ids:
  - `test_dual_graph.py::test_k10_a_different_pip_moves_the_unified_set`: FAILS alone (known; root cause in P6).
  - `test_dual_graph.py::test_no_pips_warns_and_unifies_nothing[pips2]`: passes.
  - `test_symbolic_dynamics.py::test_transition_graph_leaves_inert_classes_out`: passes.
  - `numerics/test_workbench_split.py::test_new_modules_import_cleanly`: passes.
  - `test_branch_point.py::test_insert_next_iterate_error`: passes.

### Deletion ledger, Phase 0

None. Phase 0 makes no test edits.

### Deviations, Phase 0

1. **P1:** literally neither (a) nor (b). The rows flip, so (b)'s condition
   holds, but the flips are geometrically correct, so (b)'s premise ("bug")
   fails. Applied instead: a three-part openings law and no xfail. Author
   sign-off needed (E7).
2. **P2:** the inversion saddle completes the whole pipeline. The only law
   failure is the anchor count. Several laws are vacuous (no holes at 6 steps;
   7 steps is infeasible) and go into `NOT_APPLICABLE` rather than
   `KNOWN_ISSUES`.
3. **P6:** registry-id nondeterminism also affects k10, and it is the root
   cause of the known isolation failure. Recorded for Phase 1/7.
4. The coverage-guard whitelist also includes the two `forward_*_branch_cycle`
   wrappers (dead API under decision 4) beyond the four planner entries.

---

## §1 Phase 1: infrastructure (additions only)

### What landed

- **`tests/cases.py`** (tests only, parameters frozen there): `Case`
  (session, fixed points outermost first, pips, zones, stage, `CaseExpect`),
  `Stage` (`empty` … `partitioned`), builders `build_k10`, `build_k28(blasts=)`,
  `build_period3`, `build_nested(through=)`, `build_inversion`,
  `build_orientation_reversing`; `LAW_CASES`, frozen `BUILDERS` (the six law
  cases plus `nested_unblasted`, `orientation_reversing`, `p3_15`, `p3_16`,
  `p3_6_blasts`), `KNOWN_ISSUES` (`KnownIssue(reason, raises)`),
  `NOT_APPLICABLE`, `law_params(law_id)` (xfail strict / skip marks).
  k10, k28 and inversion are reimplemented; p3 and the partitioned nested build
  delegate to `henon_cases` with every keyword explicit (E4 default).
- **`tests/helpers/`**: `invariants.py` (moved from `tests/invariants.py`;
  `assert_cdist_monotonic` now defaults to `strict=False`; the spike message
  prints the real median segment and ratio), `fakes.py` (moved from
  `tests/walk_helpers.py`, plus `bare_fixed_point()` returning a REAL
  `FixedPoint`), `laws.py` (manifold, crossing, anchor and bridge checks, each
  returning its item count; `LAYERS`, `run_layer`), `names.py` (letter-free
  spellers: `class_pair`, `homotopy_pair`, `word_in_names`,
  `refined_children`, `itinerary_pairs`, `brackets`, `matrix_by_classes`),
  `logs.py` (`assert_logged(caplog, level, logger, count=)`).
  `tests/invariants.py` and `tests/walk_helpers.py` are one-line re-export
  shims until Phase 10.
- **`tests/conftest.py`**: `matplotlib.use("Agg")`; every fixture is a thin,
  function-scoped builder call (`k28_partitioned`,
  `k28_two_blasts_partitioned` lose session scope, `henon_inversion` loses
  module scope); a new indirect `law_case` fixture. Fixture names and return
  shapes are unchanged.
- **`pyproject.toml`**: `golden` and `perf` markers registered;
  `addopts = "-ra --strict-markers -m 'not perf'"`. `slow` and `regression`
  stay registered until Phase 10.

### Builder parameters kept (from P6)

| Fixture | Builder call | Why |
|---|---|---|
| `grown_unstable` / `grown_both` / `small_tangle` | `build_k10(unstable_steps=7, through=…)` | the 7-step numerics stages are what those tests were written for |
| `henon_tangle_with_bridges`, `k10_session` | `build_k10(through="bridges")` (9 steps, infer on) | P6 variant A: interchangeable with the old 7 + 2 recipe |
| `henon_p3_session` | `build_nested(outer_blasts=0, p1_area_cutoff=1e-4, through="zones")` | P6 B4/B7: the unpartitioned stage and the 1e-4 outer cutoff both matter |
| `p3_partitioned` | `build_nested(outer_blasts=0, p1_area_cutoff=1e-4)` | P6 B4 (delegated, repin included; equivalent) |
| `k28_*` | `build_k28(blasts=1|2)` function-scoped | P6 variant C |

The k10 builder keeps the historical recipe without an analytic Jacobian
(checked: seed `[4,-4]` vs `saddle_guesses` give identical crossings without
the Jacobian; WITH the Jacobian the crossing cdists differ at the bit level).

### Verification

- Collected **796**, node ids identical to `nodeids_base.txt`.
- `795 passed, 1 skipped` in 111 s with coverage (baseline 84 s). The +27 s
  is the inversion fixture going function-scoped (about 2.3 s per user under
  coverage, 10 users), above the planner's +15 s estimate; Phase 4 (law-tier
  module build) and Phase 7 consolidation remove most of those users.
- `-rxX`: no xfail, no xpass.
- Coverage guard vs `cov_base.json`: OK (`cov_p1.json`). A first run missed
  `TangleSession.partition_element_for`'s stale-trellis `continue`: it was
  covered only through the ORDER of the session's trellis cache, which the old
  `henon_p3_session` filled outermost first (`T1` then `T3`). The local nested
  recipe classifies per trellis in that order (not via the fan-out, which
  goes in construction order, fp3 first), restoring the coverage.
- Isolation (5 node ids): `test_dual_graph.py::test_k10_a_different_pip_moves_the_unified_set`
  FAILS alone (known, P6 root cause; fixed in Phase 7);
  `test_session_bridge_classes.py::test_k28_letters_only_the_anchor_class`,
  `test_partition_elements.py::test_element_for_resolves_a_bridge_missing_from_the_snapshot`,
  `numerics/test_geometry.py::test_arc_polyline_clips_to_the_cdist_window`,
  `numerics/test_map_step_and_graph.py::test_branch_cycle_is_the_advance_key_chain[stable-fp3]`
  pass.
- Scratch check of `helpers.laws` on the six law cases: every layer passes
  except `check_one_anchor_per_unstable_branch` on inversion (2 anchors per
  branch, the KNOWN_ISSUE); `build_orientation_reversing()` raises
  `ValueError`. `helpers.names` on k10 reproduces the CLAUDE.md rows, words,
  refinement and `[[1,1],[1,1]]` matrix in element names.

### Deletion ledger, Phase 1

None. Phase 1 deletes no test. (`tests/invariants.py` and
`tests/walk_helpers.py` were moved into `tests/helpers/` and replaced by
re-export shims; no node id changed.)

### Deviations, Phase 1

1. **`build_nested(through=…)` is not a pure delegation.** `henon_cases` has
   no stage parameter and `src/` is out of scope, so a build stopped before
   the partition (`through` ≤ `"zones"`, no blasts) runs a local copy of the
   recipe; `"partitioned"` delegates. The local copy classifies per trellis,
   outermost first (see the coverage note).
2. **Stages beyond the plan's list:** `empty`, `fixed_point`,
   `grown_unstable` were added so the k10 numerics fixtures (`workbench`,
   `fixed_point`, `grown_unstable`) are builder calls too.
3. **`KNOWN_ISSUES` values are `KnownIssue(reason, raises)`**, not bare
   strings, so the b = -1 placeholder carries `raises=ValueError`.
4. **Area law is checked link by link** in `helpers.laws` with
   `AREA_LINK_RTOL = 1e-2` (the `collision_rtol` default of `StrongPip` /
   `Pseudoneighbor`). Comparing every chain member against the head with
   `rtol = 1e-3` (the `assert_area_preserved_along_chain` default) FAILS on
   p3 and nested: one link reaches 1.4e-3 (k10 1.8e-4, k28 two blasts
   2.3e-4). Phase 4 decides the final tolerance.
5. **`helpers.fakes` is moved, not yet rewired** (planner §A Phase 1:
   "Phase 7 rewires it"): `FakeFixedPoint`, the hand-written
   `FakeTrellis.image_cdist` and the derived-owner default of `make_result`
   stay until Phase 7; `bare_fixed_point()` is added now.
6. **`laws.py` covers the numerics layers only**; the topology layers land
   with the law-tier modules in Phase 4 (the task allowed a scaffold).
