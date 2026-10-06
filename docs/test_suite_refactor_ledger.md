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

---

## §2 Phase 2: pure deletions

### Deletion ledger, Phase 2

Every candidate was checked against the regression-guard notes (P5 verdict,
§0): none guards a past bug that lacks another kept guard. 74 collected node
ids removed, 2 added (796 → 724); `tests/visualizations/` (2 uncollected
scripts) deleted too.

| Node id(s) | Reason | Guard checked |
|---|---|---|
| `tests/visualizations/viz_nested_bridges.py`, `viz_unstable_artifacts.py` (uncollected) | dead scripts: not collected, crash with KeyError None, match a dead log string | none |
| `numerics/test_workbench_split.py` (whole file, 20 ids: `collaborators_are_wired…`, `delegating_signatures_unchanged[×7]`, `moved_bodies_live_in_their_new_module[×9]`, `new_modules_import_cleanly`, `top_level_import_still_works`, `visualize_signature_unchanged`, `iterate_bridge_still_patchable…`) | code-layout and signature pins (decision 4); the patchability guard is redundant (the blast error tests break loudly if patching stops working) | none |
| `test_manifold_machine.py` (whole file: `machine_initialization`, `new_grow_manifold_matches_old_period_one/three`) | frozen pre-1.8 grow loop (migration snapshot) plus a tautology | none |
| `numerics/test_bridge_identity.py::test_deleted_genealogy_attributes_are_gone` | tombstone | none |
| `numerics/test_workbench_bugfixes.py::test_grown_until_intersection_is_gone` | tombstone | none |
| `numerics/test_workbench_bugfixes.py::test_advance_key_returns_the_next_orbit_point` | tombstone half (`_advance_key_forward` gone) + duplicate of the `test_fixed_point` advance_key tests | advance_key kept in `test_fixed_point` |
| `numerics/test_map_step_and_graph.py::test_per_step_beta_matches_the_expression_it_replaces` | migration pin | k_value-root bug: `test_per_step_beta_is_the_k_th_root_only_without_inversion` + inversion `test_one_map_step_scales_cdist_by_per_step_beta` kept |
| `numerics/test_map_step_and_graph.py::test_k10_intersection_graph_is_unchanged` | migration snapshot (counts 8/18/7/7/4) | none |
| `numerics/test_map_step_and_graph.py::test_branch_cycle_reproduces_the_topology_branch_orders` | dead API (`forward_*_branch_cycle` wrappers, whitelisted) | none |
| `numerics/test_single_source_of_truth.py::test_bridge_cutting_pin_k10`, `::test_bridge_cutting_pin_period_3` | migration snapshot (6-digit bridge-cutting cdists); fingerprint helpers and pinned constants removed with them | none (refactor finished) |
| `numerics/test_single_source_of_truth.py::test_tangle_keeps_only_index_state` | tombstone | none |
| `numerics/test_cleanup_walkers_and_examples.py::test_string_dispatched_iter_method_is_gone`, `::test_stability_alias_is_defined_once_in_numerics`, `::test_no_test_or_script_defines_its_own_henon` | tombstone / code-layout | none; `test_collect_is_the_only_traversal` KEPT |
| `numerics/test_tangle_index_and_orientation.py::test_small_tangle_crossing_count_is_unchanged` | fixture-count snapshot | none |
| `test_fixed_point.py::test_fixed_point_takes_no_branch_count` | signature tombstone; `num_branches` is derived and tested | has_inversion bug: two-branch init / kevin-way tests kept |
| `numerics/test_intersection_registry_fixes.py::test_get_lambda_u_*` (×4), `::test_on_interval_stable_uses_the_stable_side_eigenvalue` | dead API (whitelisted) | none |
| `numerics/test_closed_form_curvature.py::test_curvature_area_batch_matches_scalar` | dead API (scalar `_curvature_area`, whitelisted) | `test_parabolic_fit_matches_vandermonde` kept |
| `test_session_trellis_cache.py::test_invalidate_trellises_is_a_deprecated_no_op` | dead API (whitelisted) | restore-invalidates guard kept |
| `test_resonance_zone_region.py::test_boundary_arc_dataclass_is_gone`, `::test_shapely_is_not_a_dependency` | tombstone / layout | none |
| `test_resonance_zone_region.py::test_k10_zone_area_and_containment_are_unchanged`, `::test_p3_zone_areas_and_containment_are_unchanged` | migration snapshot (areas to 1e-11, containment bit masks); masks, `bbox_grid`, `containment_mask` removed | area code stays covered by the winding test (rewritten, see below) |
| `test_manifold_initializer.py::test_initialization_unstable`, `::test_initialization_stable` | migration snapshot ("exactly 3 init points") | none |
| `numerics/test_machine_iterate_invariants.py::test_cdists_are_positive_and_increasing` | implied by the monotonicity tests | none |
| `test_minimal_trellis.py::test_image_chain_is_the_bridge_class_function`, `test_partition_family.py::test_kinds_are_distinct`, `test_element_naming.py::test_names_are_deterministic_across_builds` | tautologies | none |
| `test_topology_plotting.py::test_every_plotter_is_a_module_function_in_plotting`, `::test_stable_partition_module_holds_no_drawing_code`, `::test_hole_style_conventions_live_once`, `::test_session_exposes_one_fanout_helper_per_shape` | code-layout | none |
| `test_topology_invariants.py::test_henon_partition_still_built_after_wiring` | historical smoke | none |
| `test_image_cdist.py::test_image_cdist_is_additive_in_n_on_the_scaling_branch`, `::test_scaled_element_image_stays_on_the_advanced_branch`, `::test_every_bounded_element_now_has_an_image` | private / tautology / duplicate | KEPT `test_partition_elements::test_image_of_element_falls_back_when_an_iterate_is_missing` |
| `test_partition_elements.py::test_image_of_element_lands_on_the_advanced_branch` | subsumed | KEPT `test_image_cdist::test_image_of_element_covers_both_endpoint_images` |
| `test_session_bridge_classes.py::test_session_result_matches_gathered_partitions_helper` | private tautology (`_gathered_partitions`) | session-vs-direct tests kept |
| `test_session_bridge_classes.py::test_p3_no_class_mixes_fixed_points` | verbatim twin | KEPT `test_bridge_class::test_p3_classes_cover_both_tangles_and_mix_neither` |
| `test_session_bridge_classes.py::test_describe_reports_the_image_evidence` | string pin; evidence asserted structurally elsewhere | covered by the new report smoke test |
| `test_stable_partition.py::test_henon_reference_holes_hug_the_stable_manifold`, `::test_stable_arc_midpoint_walks_the_branch_between_same_branch_bounds`, `::test_stable_arc_midpoint_falls_back_to_the_chord_across_branches` | private display placement | coordinates still covered by KEPT `test_direct_hole_side_is_the_side_of_its_coordinates` and `test_henon_propagated_hole_side_matches_coords` |
| `test_stable_partition.py::test_forward_pairs_beyond_the_fundamental_segment_get_no_hole[3-1-False]`, `[3-2-False]` | the provisional `+1..+(k-1)` exemption (decision 3) | holes-backward-only: the firm parameters (references/backward punched; iterate ≥ k_value never) stay; docstring now names the provisional half as untested |
| `test_dual_graph.py::test_k10_summary_reports_counts`, `test_dual_walk.py::test_walk_search_repr_and_walk_repr`, `test_element_naming.py::test_k10_describe_lists_every_parent_with_its_children` | string pins (decision 5) | replaced in this commit by `test_report_smoke.py` |

### Added / changed in the same commit

- **`tests/test_report_smoke.py::test_reports_are_non_empty_and_mention_what_they_report`**
  (new): on k10, `DualGraph.summary()`/repr, a `StableNode` and a `FaceNode`
  repr, `ElementNaming.describe()` (mentions every homotopy element's name),
  `describe_bridge_classes()` (mentions every lettered class), every
  `WalkSearch`/`Walk` repr, `SymbolicDynamics.describe()`/repr: all non-empty.
- **`tests/numerics/test_bridge_identity.py::test_standalone_collaborators_agree_with_the_workbench`**
  (new): the deleted `test_workbench_split` was the only coverage of the public
  `BridgeIterator.workbench` / `IterateInference.workbench` properties
  (coverage guard flagged `BridgeIterator.py:80`, `IterateInference.py:72`).
  Rewritten behaviourally: standalone collaborators over a workbench give the
  same image/preimage genealogy as its delegates, and a second inference pass
  adds nothing.
- **`test_resonance_zone_region.py::test_zone_area_is_positive_and_winding_independent`**:
  compared the reversed-ring area with the deleted `K10_ZONE_AREA` pin at
  rel=1e-11; now compares it with the forward-ring area (`pytest.approx`).
- Orphaned imports, module constants, fingerprint helpers and empty section
  headers removed; module docstrings of `test_single_source_of_truth`,
  `test_resonance_zone_region`, `test_cleanup_walkers_and_examples` and
  `test_intersection_registry_fixes` no longer describe the deleted pins.

### Verification

- Collected **724** (796 − 74 + 2), node-id diff in
  `.refactor/runs/2026-10-05-test-suite/p2_deleted.txt` (`nodeids_p2.txt`).
- `723 passed, 1 skipped` in 108 s with coverage; `-rxX`: no xfail, no xpass.
- Coverage guard vs `cov_base.json`: OK (`cov_p2.json`); the only newly missed
  lines are whitelisted dead API (`on_interval`, `_get_lambda_u`,
  `invalidate_trellises`).
- Isolation: `test_dual_graph.py::test_k10_a_different_pip_moves_the_unified_set`
  now PASSES alone (see deviation 1), `test_pseudoneighbor.py::test_table_linked_deep_iterate_does_not_disqualify`,
  `numerics/test_generation_and_caches.py::test_point_arrays_reflect_growth`,
  `numerics/test_grow_until.py::test_grow_until_preserves_crossing_ids_across_rounds`,
  `test_report_smoke.py::test_reports_are_non_empty_and_mention_what_they_report`: all pass.

### Deviations, Phase 2

1. **`test_k10_a_different_pip_moves_the_unified_set` fixed now, not in
   Phase 7.** After the deletions it failed in the FULL run too (fewer k10
   builds before it changed the registry-id permutation, P6 root cause). Its
   alternative pip is now chosen by stable cdist among the candidates whose
   `k_value`-th iterate is registered, never by id (the P6 fix). Test-only
   change; it passes alone and in the suite.
2. **Two tests added beyond the planned smoke test** (the collaborator test
   and the rewritten winding assertion), so the coverage guard holds and the
   winding test no longer depends on a deleted migration pin.
3. The plan's summary list and planner §A differ slightly; the union was
   applied (§A adds `test_advance_key_returns_the_next_orbit_point`,
   `test_cdists_are_positive_and_increasing`,
   `test_describe_reports_the_image_evidence`). 74 ids, matching the ~75 estimate.
4. Not touched (later phases): `test_p3_forward_holes_stop_at_the_branch_return`
   still asserts the provisional `+1, +2` forward holes on p3 (decision 3);
   flagged for Phase 3b/4.

---

## §3a Phase 3a: loosen in place (`tests/numerics/`, `tests/regression/`)

### Deletion ledger, Phase 3a

None. No test function is deleted; only assertions change. 2 node ids added
(`test_walk_back_inverts_walk_fwd_through_the_root[unstable|stable]`, see
below), 724 → 726.

### What changed

- **Exceptions, type only (decision 5).** Every `match=` dropped in
  `test_cleanup_walkers_and_examples` (`collect_rejects_a_broken_iterate_link`,
  `first_node_needs_a_branch_index`, `saddle_guesses_refuses_unknown_parameters`),
  `test_grow_until` (×9), `test_workbench_bugfixes` (×3), `test_map_step_and_graph`
  (×2), `test_minimal_solver` (×2), `test_inversion_fixture` (mixed-sign) and
  `test_gpu` (no-CuPy). `test_collect_is_the_only_traversal` now raises a local
  sentinel exception class from the sabotaged `_collect` instead of matching its
  text. Where the message was what told a validation error from the cap error
  (`grow_until_intersection_honours_branch_index`, the empty-grow-set /
  unknown-id / unknown-direction / empty-id-set rejections), the test now also
  asserts that no manifold grew, so the type-only check keeps its meaning.
- **Logs, level and logger only.** `test_rebuild_drops_metadata_after_an_unpreserved_recompute`
  (DEBUG, `tanglepack.numerics.TangleWorkbench`) and
  `test_grow_until_clears_the_iterated_flag_it_invalidates` (WARNING, same
  logger) go through `helpers.logs.assert_logged`.
- **Tolerances (decision 5).** `per_step_beta` `rel=0, abs=0` → default
  `pytest.approx`; the area-preserving reciprocal and the dissipative
  non-reciprocal checks compare with `approx` instead of exact `==`/`!=`;
  `test_one_map_step_scales_cdist_by_per_step_beta` states its tolerance as
  `20 * _MIN_SEED_STEP` (= the old 1e-4); `test_refinement_essentially_converges_to_cutoff`
  measures with `_curvature_area_batch` against the machine's `area_cutoff`.
- **Strictness.** `test_machine_iterate_invariants` (both tests) and
  `test_growth_integration` use the default `strict=False` and now also run
  `assert_no_geometric_spikes` (decision 5: the no-spike check runs
  everywhere); the only strict cdist test left is
  `test_low_stretch_growth_keeps_cdist_injective`.
- **Private → public.**
  - `_intersection_registry` → `intersection_registry` everywhere in
    `numerics/` and `regression/`.
  - `test_batched_map` (×3): `_map_batchable` replaced by a recording wrapper
    around the map — a batch-capable map is called once per batch on `(2, N)`,
    a scalar-only map once per `(2,)` point, and the two growths really took
    those two paths.
  - `test_partial_pieces_are_held_outside_the_id_registry`: `_bridges` /
    `_partial_bridges` → `workbench.bridges`, `workbench.bridge(id)`,
    `id is None`.
  - `test_clear_bridges_releases_segment_ownership`: the Tangle's segment
    ownership tables → weak references (no discarded bridge survives
    `gc.collect()` after `clear_bridges`) plus a fresh `create_bridges`
    reproducing the same id set. Checked to bite: with
    `Tangle.release_manifold` monkeypatched to a no-op the test fails.
  - `test_grow_until_clears_the_iterated_flag_it_invalidates`: `_bridges` →
    `workbench.bridges`.
  - `test_iterates_closed_all_is_scoped_to_the_grown_fixed_point`:
    `_frozen_ids` → observed through the driver (see deviation 1).
  - `test_grow_until_intersection_stops_once_a_crossing_exists`:
    `Tangle._intersecting_segments` → the registry holds a non-anchor crossing.
  - `regression/test_near_vertical_refinement`: scalar `_curvature_area` (dead
    API) → the `_curvature_area_batch` kernel; it now also asserts the area is
    positive (a degenerate row returns 0, which would pass vacuously).
- `test_closed_form_curvature` module docstring no longer describes the
  scalar-vs-batch test deleted in Phase 2.

### Verification

- Collected **726** (`nodeids_p3a.txt`).
- `725 passed, 1 skipped` in 108 s with coverage; `-rxX`: no xfail, no xpass.
- Coverage guard vs `cov_base.json`: OK (`cov_p3a.json`), after deviations 2
  and 3.
- Isolation: `test_dual_graph.py::test_k10_a_different_pip_moves_the_unified_set`,
  `numerics/test_invariant_helpers.py::test_fundamental_segments_have_injective_cdist`,
  `numerics/test_geometry.py::test_interior_point_is_inside[square]`,
  `numerics/test_closed_form_curvature.py::test_parabolic_fit_matches_vandermonde[4]`,
  `numerics/test_geometry.py::test_polyline_midpoint_of_an_empty_polyline_is_none`
  pass alone, as does every rewritten test above.

### Deviations, Phase 3a

1. **`test_iterates_closed_all_is_scoped_to_the_grown_fixed_point` is not
   observed by a full driver run.** The two-saddle workbench cannot close the
   k=10 saddle's own crossings within any affordable budget (one of them never
   gets a forward image in 3 rounds; the stable growth makes each round cost
   seconds). The test instead closes them by hand through the public
   `IntersectionRegistry.register_iterate` (stand-in `+1` links; the predicate
   reads only the table) and asserts that `grow_until_iterates_closed(simple)`
   with `ids="all"` returns 0 without growing, while the other saddle's
   crossing is still open. Checked to bite: passing the whole registry's ids
   instead runs to the cap and raises.
2. **Coverage-guard whitelist extended** with `ManifoldMachine._chord_frame`,
   `_linear_fit` and `_compute_single_area`: they are called only by the scalar
   `_curvature_area` (dead API, decision 4), which the near-vertical
   regression test no longer reaches.
3. **New test `test_walk_back_inverts_walk_fwd_through_the_root`** (numerics,
   parametrized unstable/stable): the scalar curvature path was the only
   caller covering `BaseManifold.walk_back` through a root `BranchPoint`
   (`_branch_backward`), which is live public API. The test asserts the
   structural inverse `walk_back(walk_fwd(prev, node), node) is prev`, at the
   root slots and along the ordinary nodes; it pins no slot semantics.
4. **Out of 3a scope, left for later phases:** private access not named by
   the plan's Phase 3 list (`Tangle._seg_lookup` / `_manifold_segs` in the trim
   test and `test_tangle_index_and_orientation`, `_iter_manifolds`,
   `_man_machine.area_cutoff`) stays for Phase 4 (twins deleted) / Phase 7
   (trim oracle, strengthening); `test_workbench_has_no_private_key_advance`
   is upgraded to the advance_key spy in Phase 7 per planner §A. The
   numerics-layer files at the top level (`test_fixed_point.py`, …) belong to
   Phase 3b by directory.

---

## §3b Phase 3b: loosen in place (top-level topology, loom, session files)

### Deletion ledger, Phase 3b

None. No test function is deleted and no node id changes (726 → 726,
`nodeids_p3b.txt` identical to `nodeids_p3a.txt`); only assertions change.
Every rewritten assertion was checked against the regression-guard notes: the
guards keep their bite (collision-rtol unlinked halves unchanged in value;
blast `AssertionError` still propagates; restore, repunching, own-key and
backward-endpoint guards untouched).

### What changed

- **Exceptions, type only (decision 5).** All 42 `match=` arguments dropped in
  `test_arrangement` (2), `test_bridge_class` (6), `test_dual_graph` (6),
  `test_dual_walk` (4), `test_element_naming` (4), `test_fixed_point` (3),
  `test_higher_period_cartoon` (2), `test_loom_blast_restore` (2),
  `test_partition_family` (5), `test_session_bridge_classes` (1),
  `test_topology_plotting` (8). `excinfo` message checks replaced by
  structure: `test_split_origin_sides_raise` → one entry in
  `bridge_side_violations`; `test_orientation_reversing_requires_alternating_sides`
  → `bridge_side_violations(..., orientation_preserving=False)` non-empty;
  `test_bridge_rows_inconsistent_raise` → `bridge_row_violation(...) is not None`;
  `test_a_missing_partition_names_the_branch` → the same call succeeds once
  the missing partition is back.
- **Logs, level and logger only**, through `helpers.logs.assert_logged`:
  `test_dual_walk` (11 checks; "no warning" checks become `count=0`, the
  start/goal pair `count=2`), `test_minimal_trellis` (3),
  `test_dual_graph::test_no_pips_warns_and_unifies_nothing`,
  `test_iterated_partition` (row invariant, no empty stretch),
  `test_pseudoneighbor` (uncut-branch pip warning, E6; silent case `count=0`),
  `test_higher_period_cartoon` (foreign crossing INFO, other-branch lobe
  WARNING, blast-order `count=0`), `test_session_bridge_classes` (zone
  straddle, missing partition), `test_symbolic_dynamics` (virtual class, warns
  no more → `count=0`, cross side, inert-with-word, no evidence),
  `test_stable_partition::test_partition_warns_without_pseudoneighbors` (exact
  count dropped: each of the two calls warns), `test_bridge_class` (unresolved
  loop), `test_topology_plotting` verbose tests (INFO from
  `tanglepack.topology.Trellis`; the print fallback asserts non-empty stdout,
  E6), `test_loom_blast_restore` (skipped bridge).
- **Reason strings and reprs.** `ElementLanding.reason` is asserted present,
  not worded; landing / `Cut` / `ElementName` / partition reprs non-empty.
- **`describe()` / `summary()` (decision 5): non-empty and mentions.**
  `test_describe_reports` (stable partition), `test_describe_mentions_parents_and_cuts`,
  `test_describe_has_one_line_per_element`, `test_table_lookups_symbols_and_report`
  and `test_describe_reports_letters_and_inertness` (every class's `name`),
  `test_homotopy_only_naming_is_its_own_parent` and
  `test_branch_codes_always_and_letters_only_with_two_fixed_points` (every
  element name), `test_is_reliable_and_describe` and
  `test_describe_symbolic_dynamics` (every class letter),
  `test_active_class_lies_in_the_zone…` (`active.describe()` non-empty),
  `test_k28_letters_only_the_anchor_class` (counts from the table, not the
  report), the minimal-trellis summary/describe, the sparse-arrangement
  summary (`"sparse"` dropped; the `+2` / `+4 dangling ends` pins become
  virtual-node counts). `ElementRef` labels: only "equal refs share a label,
  different refs differ" (`test_element_ref_hashes_compares_and_labels`,
  `test_fixed_point_label_wins_…`). Delegation equalities
  (`session.describe_* == product.describe()`) are kept: they pin no wording.
- **Evidence (decision 2).** `source ==` removed everywhere
  (`test_k10_active_class_word`, `test_k28_one_blast_singleton_path` — its
  `search is None` twin too — and `test_k28_two_blasts_structure`).
  `is_reliable` / `verified` removed from the synthetic tests:
  `test_is_reliable_and_describe` now checks the unresolved listing and the
  report, `test_ambiguous_class_is_unreliable` checks
  `ClassDynamics.status == "ambiguous"`. The fact tests (k10 word, k28 one
  and two blasts, p3 equivariance, k10 session word) keep them until Phase 5.
- **Cache identity in fact tests (decision 6).** Removed
  `dyn is session.symbolic_dynamics([fp])` (`test_k28_two_blasts_structure`),
  the `dyn.table is …` / `dyn.naming.iterated is …` half of
  `test_k10_symbolic_dynamics_is_built_over_the_cached_pieces`, and
  `len(session._bridge_classes) == 2` (`test_cache_is_kept_per_fixed_point_selection`).
- **Tolerances.** `min_unstable_cdist == 0.0` → `approx(0, abs=cdist_tol)`
  (`test_classes_and_members_come_out_in_the_documented_order`,
  `test_k28_two_blasts_structure`). The table-linked tests state their drift as
  `2 * collision_rtol`, read from the `compute_pseudoneighbors` /
  `is_strong_pip` signature defaults (values unchanged: the unlinked halves are
  NOT loosened, critic note 7).
- **Private → public.** `_built_generation` (cache-hit test: the hit does not
  move the workbench generation; minimal-trellis check dropped, registry and
  manifold sharing kept); `TangleSession._bridge_test_point` →
  `geometry.polyline_midpoint(bridge.get_point_array())` (zone midpoint test,
  blast frontier helper).
- **Provisional rule (decision 3).** `test_p3_forward_holes_stop_at_the_branch_return`
  no longer asserts that the +1 and +2 forward holes exist; it asserts the
  firm half (no hole at iterate ≥ `k_value`, while such pairs are recorded)
  and names the exemption as untested (flagged in Phase 2 deviation 4).

### Verification

- Collected **726** (`nodeids_p3b.txt`, identical to 3a).
- `725 passed, 1 skipped` in 108 s with coverage (`p3b_run.txt`); `-rxX`: no
  xfail, no xpass.
- Coverage guard vs `cov_base.json`: OK (`cov_p3b.json`), after deviation 2.
- Isolation: `test_dual_graph.py::test_k10_a_different_pip_moves_the_unified_set`,
  `test_bridge_class.py::test_p3_classes_cover_both_tangles_and_mix_neither`,
  `test_session_trellis_cache.py::test_same_generation_is_a_cache_hit`,
  `test_stable_partition.py::test_p3_forward_holes_stop_at_the_branch_return`,
  `test_dual_walk.py::test_no_start_node_and_no_goal_node`,
  `test_pseudoneighbor.py::test_table_linked_deep_iterate_does_not_disqualify`,
  `test_stable_partition.py::test_partition_warns_without_pseudoneighbors`:
  all pass alone.

### Deviations, Phase 3b

1. **Private access not on the Phase 3 list is left for later phases**, where
   the tests holding it are deleted, moved or rewritten: `build_pieces` /
   `_gathered_partitions` / `_partition_signature` (Phase 6/7), the
   arrangement `_nodes` / `_half_edges` tables (including the new
   `_dangling_ends` helper that replaces the summary pin) and
   `DualGraph._nodes_of_edge` (Phase 4/7), `_row_at`, `_row_polyline`,
   `_image_chain`, `_scaled_image_cdist`, `_snap_to_partition_boundary`,
   `_strong_pip_cuts` (E6), `_inert_letter` / `_active_letter`,
   `_bridge_side_of`, `_region_key`, `_is_forward_beyond_fundamental`,
   `_intervals_from_marks` / `_marks_of` / `_build_intervals`, plotting
   `_cartoon_name` / `_anchorward_look` (Phases 7/8). The allowed kernels
   (`_side_of`, `_resolve_inertness`, `_backward_endpoint`,
   `_bridge_unstable_span`, `_containing_bridge`, `_near_far`, the chord
   kernel `_empty_stretches`) stay.
2. **One repr kept for coverage.** Dropping `"toward {partner}" in repr(cut)`
   left `Cut.__repr__`'s partner branch (`PartitionFamily.py:574`) uncovered;
   `test_cuts_record_the_far_end_of_the_empty_stretch` now asserts each
   partnered cut's repr is non-empty.
3. **Kept as specifications, not wording:** element-name texts
   (`R_(1.0;1)^2`, mathtext), symbol texts (`a^-1`, `a_2`), the inert letter
   series and cut `reason` codes are author-defined formats in CLAUDE.md; the
   `StableNode.label` join (`"L | R"`) is left for Phase 8 (label-separator
   tests). The session-vs-direct `_same_dynamics` helper still compares
   `is_reliable` between the two builds (an equivalence, not a reliability
   pin); Phase 6 rewrites it.
4. Test names describing the old wording checks (`test_describe_mentions_parents_and_cuts`,
   `test_is_reliable_and_describe`, `test_a_missing_partition_names_the_branch`,
   …) are kept so no node id moves; Phases 4–10 rename or delete them.
