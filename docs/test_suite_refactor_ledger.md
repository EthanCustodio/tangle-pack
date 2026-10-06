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

---

## §4 Phase 4: the physical-law tier (`tests/invariants/`)

### What landed

- **`tests/invariants/`** (11 modules, one test per (law × case)):
  `test_law_case_sanity` (4 laws + the b = -1 placeholder + a wiring meta
  test), `test_law_manifolds` (4), `test_law_crossings` (9, two of them on
  fresh builds), `test_law_anchors` (2), `test_law_bridges` (7),
  `test_law_partition` (14, the P1 openings law among them),
  `test_law_arrangement` (5), `test_law_classes` (9), `test_law_iterated` (7),
  `test_law_dual_graph` (5), `test_law_symbolic` (8). 74 laws, 446 cases.
- **`tests/invariants/conftest.py`**: module-scoped indirect `law_case` (one
  read-only build per module per case, every product built inside the
  fixture), function-scoped `fresh_law_case` for the two mutating laws
  (`recompute_preserves_ids`, `recompute_is_idempotent`), and the autouse
  FINGERPRINT GUARD (workbench and registry generation, crossings, bridge
  ids, holes with openings, pips, partition signature, session alphabet, the
  identity of every cached product). Checked to bite: a scratch test growing
  the shared k10 build errors in teardown.
- **`tests/helpers/laws.py`**: every check is `check_<law>(case) -> int`; the
  law id is the name without `check_`. Topology layers added (partition,
  arrangement, classes, iterated, dual graph, symbolic). Tolerances read from
  library defaults: `AREA_LINK_RTOL` = `compute_pseudoneighbors`'
  `collision_rtol` (1e-2; settles Phase 1 deviation 4), `LINK_SCALING_RTOL` =
  `infer_iterate_table`'s `cdist_rtol`, `REGION_AREA_RTOL` =
  `Arrangement.image_of`'s `rtol`, cdists against `cdist_tol`. The geometric
  oracles (crossing sign, row of a bridge end, a direct hole's side, a face's
  side) are measured on the curves with public geometry and the `_side_of`
  kernel. `run_layer` (unused) removed.
- **`tests/helpers/law_tier.py`**: `build_products`, `fingerprint`,
  `assert_non_vacuous` (count > 0, else the pair belongs in `NOT_APPLICABLE`),
  `law_test(check)` (the parametrized test factory).
- **`tests/cases.py`**: `KNOWN_ISSUES` / `NOT_APPLICABLE` keyed by real law ids
  (the generic `("holes", "inversion")` key replaced per law); `issue_marks`
  factored out of `law_params`.
- **`tests/conftest.py`**: the function-scoped `law_case` removed (superseded
  by the module-scoped one in `tests/invariants/conftest.py`); docstring
  updated.

### P1 outcome: the openings law (no xfail)

Four laws in `test_law_partition`, all passing where applicable:
`direct_hole_opens_inward_pair` (outward of the near bound, anchorward of the
far, on its own bridge's row), `openings_on_own_bridge_row` (every opening of
every hole sits on the row of the hole's bridge at that bound: this is what
makes the P1 "flips" legitimate), `openings_linked_bound` (at a registered
iterate of the origin's bound the opening is the origin's; not applicable on
k10, k28 two blasts, inversion), `openings_missing_only_at_anchor_or_tail`.
The two snapshots `test_p3_hole_sides_are_pinned` and
`test_k28_two_blast_hole_sides_are_pinned` are deleted in this commit.
AUTHOR ITEM E7 stands (confirm the row of a backward hole may differ from its
origin's when its image lobe's stable side was never computed).

### KNOWN_ISSUES (all `xfail(strict=True)`)

| (law, case) | Reason | New? |
|---|---|---|
| `one_anchor_per_unstable_branch`, inversion | two (0,0) anchors per unstable branch (P2) | P2 |
| `orientation_reversing_case_builds`, orientation_reversing | `set_k_value` raises `ValueError` (P4), `raises=ValueError` | P4 |
| `arrangement_euler`, inversion | V - E + F = 0 on the one component (3 open faces) | **Phase 4** |
| `anchor_bridge_class_leads_its_tangle`, inversion | the bridge leaving the second anchor sits in an inert class | **Phase 4** |
| `itinerary_pairs_same_side`, inversion | a walked pair `L` of branch 0.1 with `R` of branch 0.0 (even length holds) | **Phase 4** |
| `landings_contained`, inversion | a landing straddles a cut (`contained=False`) | **Phase 4** |
| `arrangement_regions_disjoint`, k28_one_blast and k28_two_blasts | a zero-area region bounded by an unstable arc no `Bridge` spans reports a representative point inside its neighbour | **Phase 4** |

The four inversion findings are probably downstream of the anchor issue
(unverified; source fixes are out of scope). The k28 one is NOT on inversion:
the arrangement of a blasted trellis classifies a zero-area, bridgeless
2-corner face as a region (AUTHOR ITEM). The P3 deep-p3 issues stay listed for
Phase 9.

### NOT_APPLICABLE (skip with reason; 31 cases)

Inversion: the 16 hole/refinement/image-bridge laws (no pseudoneighbors at the
only feasible depth), `dual_face_side_is_geometry` (the minimal arrangement
closes no region), `arrangement_image_of` / `_preimage_inverts_image`. k10:
`arrangement_image_of` / `_preimage_inverts_image` (its one carrier is a
sub-face), `partition_singletons`, `openings_linked_bound`. k28 two blasts:
`partition_singletons`, `openings_linked_bound`. k28 one blast and p3:
`iterated_cut_provenance` (no image bridges), `refined_children_inherit_word`,
`member_matching_consistent` (no split class).

### Deletion ledger, Phase 4

95 collected node ids removed, 449 added (446 law-tier cases and 3 outside
it): 726 → 1080. Node-id diff in
`.refactor/runs/2026-10-05-test-suite/p4_deleted.txt` (`nodeids_p4.txt`).
Every candidate was checked against the regression-guard notes: each guarded
bug keeps a guard as a law that runs on the case where it broke (column 3).

| Node id(s) | New home (law, cases) | Guard checked |
|---|---|---|
| `test_arrangement::test_every_detected_crossing_has_a_definite_sign`, `::test_crossing_sign_is_the_cross_product_of_the_two_directions` | `crossing_sign_is_cross_product` (all six; blast crossings read on the bridge that ends there) | none |
| `test_arrangement::test_anchor_sign_matches_the_oriented_eigendirections` | `anchor_sign_matches_eigendirections` | none |
| `test_arrangement::test_signs_alternate_along_a_stable_branch` | deleted: its own docstring says it is not an invariant | none |
| `test_arrangement::test_k10_/test_p3_arrangement_is_*_component*_and_euler_holds`, `::test_k10_/test_p3_regions_are_pairwise_disjoint`, `::test_k10_/test_p3_arrangement_regions_are_geometrically_sound`, `::test_k10_/test_p3_region_images_agree_with_the_dynamics`, `::test_preimage_inverts_image` (9) | `arrangement_euler`, `_regions_disjoint`, `_regions_sound`, `_image_of`, `_preimage_inverts_image` | regions memory (nested swallow bug): disjoint law on nested |
| `numerics/test_machine_iterate_invariants.py` (whole file, 7) | `manifold_cdist_monotone` / `_no_spikes` / `_iterate_law` / `_one_to_one` | cdist-strict: non-strict + no-spike on all cases; the strict low-stretch test is kept |
| `numerics/test_growth_integration::test_bridges_individually_satisfy_invariants` | the manifold laws now include every bridge | none |
| `numerics/test_tangle_intersection_cdist.py` (whole file, 3) | `crossings_cdists_defined`, `crossings_unstable_by_stable` + `no_same_stability_crossing`, `crossing_cdist_bracketed` (unstable side, through the public `unstable_segment`) | straddle cdist: `regression/test_boundary_straddle_cdist` kept |
| `numerics/test_no_same_stability_crossing.py` (1) | `no_same_stability_crossing` (geometric, manifolds + blast image bridges, all branches and fixed points, u×u and s×s) | numerics-test-suite (self-crossing): law on all six |
| `numerics/test_bridge_identity::test_bridge_id_is_endpoints_in_unstable_order`, `::test_bridge_id_ordering_on_period_three`, `::test_bridges_at_indexes_exactly_the_two_endpoints`, `::test_bridges_at_consistency_on_period_three`, `::test_image_and_preimage_round_trip` | `bridge_id_in_unstable_order`, `bridges_at_exact`, `bridge_image_round_trip` | none |
| `numerics/test_bridge_identity::test_bridge_identity_on_the_inversion_saddle` | id/order → inversion law params + anchor xfail; its iterate half kept as the NEW `test_iterating_an_inversion_bridge_lands_on_the_other_branch` | inversion advance_key path kept |
| `numerics/test_workbench_bugfixes::test_bridge_endpoints_stay_on_the_bridges_own_branch` | `bridge_endpoints_on_own_branch` (p3, nested) | codebase-audit root cause: law on p3 and nested |
| `numerics/test_single_source_of_truth::test_every_registered_crossing_carries_its_unstable_segment`, `::test_bridge_endpoints_are_the_ids_of_the_crossings_it_was_cut_at` | `crossing_cdist_bracketed`, `bridge_id_in_unstable_order` + `bridge_endpoints_on_own_branch` | none |
| `numerics/test_blast_no_overlap::test_workbench_keeps_single_copy_per_bridge` | `bridge_single_copy`, `bridges_do_not_overlap` (blasted k28 and nested) | numerics-test-suite single copy: `test_iterating_fixed_point_bridge_returns_existing_copies` kept |
| `numerics/test_inversion_fixture::test_the_fixture_really_is_an_inversion_point`, `::test_growth_preserves_the_invariants_on_both_branches` (×4), `::test_forward_iterates_are_registered_on_the_inversion_tangle`, `::test_intersections_carry_keys_on_both_branches` | `k_value_matches_inversion`, the manifold laws, `iterate_link_scales_by_beta`, `every_branch_is_built` / `every_branch_crosses` (inversion) | has_inversion: kevin-way two-branch test kept; k_value-root: `test_one_map_step_scales_cdist_by_per_step_beta` kept |
| `numerics/test_invariant_helpers::test_area_preserved_along_every_recorded_iterate_chain` | `area_along_iterate_links` | none |
| `test_bridge_class::test_row_of_end_agrees_with_the_geometry_on_k10/_p3` | `row_of_end_is_geometry` | none |
| `test_bridge_class::test_same_branch_bridges_never_mismatch_on_k10/_p3` | `bridge_rows_consistent` (I2) | none |
| `test_bridge_class::test_every_non_partial_bridge_lands_in_exactly_one_class` | `every_bridge_in_one_class` | none |
| `test_bridge_class::test_k10_anchor_bridge_runs_source_to_target` | `anchor_bridge_class_leads_its_tangle` | none |
| `test_bridge_class::test_direction_matches_the_element_order_at_the_ends`, `::test_oriented_class_orders_anchor_outward` | `class_orientation_anchor_outward` | none |
| `test_bridge_class::test_classes_and_members_come_out_in_the_documented_order` | `class_table_order` (tangle grouping, then min cdist, then `class_sort_key`) | none |
| `test_bridge_class::test_p3_classes_cover_both_tangles_and_mix_neither` | `classes_do_not_mix_tangles` | none |
| `test_bridge_class::test_p3_period_three_classes_use_every_stable_branch` | `classes_use_every_orbit_branch` | planner guard: law on p3 and nested |
| `test_dual_graph::test_k10_/test_p3_node_structure`, `::test_k10_/test_p3_face_side_agrees_with_the_geometry`, `::test_k10_/test_p3_payload`, `::test_k10_unified_nodes_…`, `::test_p3_unifies_…`, `::test_k10_graph_is_bipartite` | `dual_node_structure` (now with the wall raises, face corners/bridge ids), `dual_face_side_is_geometry`, `dual_face_nodes` (with `image_face`), `dual_unified_is_pip_segment`, `dual_bipartite_degree`; payload deleted as tautological | none |
| `test_partition_elements::test_element_ids_are_positional_…`, `::test_element_of_intersection_covers_…`, `::test_simple_tangle_elements_own_…`, `::test_singleton_elements_own_…` | `partition_covers_branch`, `partition_unique_owner`, `partition_singletons` | none |
| `test_stable_partition::test_henon_holes_are_classified`, `::test_henon_partition_covers_branch`, `::test_henon_direct_hole_side_reproduces_inward_pair` | `holes_are_classified`, `partition_covers_branch`, `direct_hole_opens_inward_pair` | pseudoneighbor-partition: inward pair law on all hole cases |
| `test_stable_partition::test_p3_propagation_terminates_at_periodicity` | `propagation_terminates` | pseudoneighbor-partition (07-06 termination): law on p3 and nested |
| `test_stable_partition::test_k28_blast_child_gets_no_forward_hole`, `::test_p3_forward_holes_stop_at_the_branch_return` | `no_direct_hole_beyond_fundamental` (firm half only) | holes-backward-only: firm half on all hole cases |
| `test_stable_partition_period3::test_p3_holes_share_bridge_side`, `::test_p3_bridge_rows_consistent`, `test_topology_invariants::test_henon_holes_share_bridge_side`, `::test_henon_bridge_rows_consistent` | `holes_share_bridge_side` (I1), `bridge_rows_consistent` (I2) | hole-side 2026-10-02: I1 law on p3, nested; deep-run I1 test kept |
| `test_stable_partition_period3::test_p3_propagated_holes_land_on_the_predicted_branch` | `propagated_holes_land_on_predicted_branch` | none |
| `test_stable_partition_period3::test_direct_hole_side_is_the_side_of_its_coordinates[×4]` | `direct_hole_side_is_coordinate_side` (critic note 1; nested added) | hole-side 2026-10-02 |
| `test_stable_partition_period3::test_p3_hole_sides_are_pinned`, `::test_k28_two_blast_hole_sides_are_pinned` | the openings law (P1) | none |
| `test_pseudoneighbor::test_henon_reference_pairs_are_structurally_valid` | `reference_pairs_valid` (pip window) | none |
| `test_iterated_partition::test_k10_/p3_/k28_iterated_partition_invariants` | `iterated_child_inside_parent`, `_keeps_homotopy_boundaries`, `_unique_owner`, `_cut_provenance`; the test-side `_empty_stretches` oracle deleted | none |
| `test_minimal_trellis::test_k10_/p3_/k28_minimal_trellis_invariants`, `::test_k10_minimal_trellis_drops_something` | `minimal_trellis_bridges`, `minimal_trellis_nodes_and_arrangement` (public API; the private half-edge check dropped); "drops something" is not a law | none |
| `test_symbolic_dynamics::test_k10_every_landing_contained` | `landings_contained` | none |
| `test_symbolic_dynamics::test_k10_itineraries_even_and_matrix` | law half → `itineraries_even`, `itinerary_pairs_same_side`, `inert_classes_outside_transitions`; fact half renamed `test_k10_transition_matrix` (Phase 5 golden) | none |
| `test_element_naming::test_k10_names_agree_with_the_partition_structure` | `names_agree_with_structure` | none |

### Added / changed in the same commit (outside the tier)

- `numerics/test_bridge_identity::test_iterating_an_inversion_bridge_lands_on_the_other_branch`
  (the iterate half of the deleted inversion test; the coverage guard flagged
  `BridgeIterator` 239/246/357, `IterateInference` 415, `ManifoldMachine`
  317–352 and `Tangle` 797, reached only through it).
- `test_minimal_trellis::test_an_active_hole_bridge_without_a_registered_image_is_unmapped`
  (unblasted nested: no law case has an unmapped hole bridge; covers
  `MinimalTrellis.describe`'s unmapped line).
- `test_symbolic_dynamics::test_k10_transition_matrix` (rename, see above).
- `test_iterated_partition::test_cuts_record_the_far_end_of_the_empty_stretch`
  checks the cut partner against the library kernel
  `IteratedHomotopyPartition._empty_stretches` (allowed kernel) instead of the
  deleted test-side reimplementation.
- `test_minimal_trellis::test_a_pair_with_no_bridge_object_…` also asserts the
  report is non-empty; `test_report_smoke` adds the minimal trellis's
  summary/repr/describe.
- Module docstrings of every file that lost tests now point to the law that
  replaced them; empty section headers removed; orphaned helpers and imports
  removed.

### Verification

- Collected **1080** (726 − 95 + 449; `nodeids_p4.txt`); the law tier alone
  is 446 cases.
- `1040 passed, 32 skipped, 8 xfailed` in 145 s with coverage (`p4_run.txt`),
  about 78 s without. `-rxX` lists exactly the 8 `KNOWN_ISSUES` xfails above;
  no XPASS. Skips: 31 `NOT_APPLICABLE` + the GPU no-CuPy skip (Phase 9).
- **Law tier wall time: 37 s** without coverage (`pytest tests/invariants`,
  407 passed, 31 skipped, 8 xfailed), within the ~40 s budget; one
  module-scoped build per case per module (~2.3 s per module, inversion 54 %)
  plus 12 fresh builds for the two recompute laws.
- Coverage guard vs `cov_base.json`: OK (`cov_p4.json`), after the two added
  tests and the dual-graph law additions above.
- Isolation: `test_dual_graph.py::test_k10_a_different_pip_moves_the_unified_set`,
  `invariants/test_law_partition.py::test_openings_linked_bound[p3]`,
  `invariants/test_law_crossings.py::test_recompute_is_idempotent[nested]`,
  `numerics/test_bridge_identity.py::test_iterating_an_inversion_bridge_lands_on_the_other_branch`,
  `test_element_naming.py::test_element_name_is_hashable_and_compares_by_value`,
  `test_loom_blast_restore.py::test_restore_rebuilds_bridges`: all pass alone;
  `pytest tests/invariants` alone passes.

### Deviations, Phase 4

1. **Per-check form** (plan.md: one test per (law × case)), not the
   planner's per-layer fallback: the module-scoped shared build (decision 7)
   keeps it at 37 s. It brings the suite to 1080 cases for now (446 in the
   tier); Phases 5–7 delete the remaining fact pins and duplicates.
2. **Kept, not deleted:** `numerics/test_growth_integration::test_unstable_invariants_hold_at_every_step`
   (the per-step check is not subsumed by a law on finished builds; planner §B
   moves it to `unit/numerics/test_growth.py`), and
   `test_single_source_of_truth::test_partial_bridge_reports_its_missing_endpoint`
   (no law case carries a partial bridge, so `bridge_partial_iff_no_id` cannot
   cover partials).
3. **Deleted beyond the planner's §A list** because their assertions moved
   into a law in this commit: the crossing-sign and anchor-sign tests,
   `test_preimage_inverts_image`, `test_image_and_preimage_round_trip`,
   `test_p3_period_three_classes_use_every_stable_branch`, and the four I1/I2
   re-check tests (planned for Phase 9; the production-wiring spy still lands
   there).
4. **Kept from the "twins" list:** `test_dual_graph::test_k10_face_nodes_…`
   and `::test_p3_merges_…` (not in the planner's list; the face law is more
   general but the nested merge fact stays), and the iterate half of the
   inversion bridge test (re-added, see above).
5. **New known issues** (5 laws, 6 cases) beyond P2's draft; see the table.
   The k28 `arrangement_regions_disjoint` finding is not inversion-related.
6. **Crossing cdist bracket** is checked on the unstable side only (through
   the public `unstable_segment`); the old stable-side half read the private
   `Tangle._seg_lookup`.
7. **Crossing sign** is rebuilt from the curve the crossing lies on: a blast
   registers crossings on image bridges, which the manifold nodes do not pass
   through (the old test ran on the unblasted 7-step k10 only).
8. **Same-stability law** includes the blast's image bridges (curves sharing
   no node with a manifold) and runs inside the crossings' box padded by 10 %
   (the old p3 test used a fixed |x|, |y| < 6 box).
9. **The anchor law** now requires exactly one anchor on EVERY unstable
   manifold (the Phase 1 check only looked at branches that had one).
10. **The b = -1 placeholder** is a standalone test with `issue_marks`
    (through `law_params` it would also run the six law cases).
11. `tests/conftest.py`'s function-scoped `law_case` is removed (superseded).

---

## §5 Phase 5: the golden tier (`tests/golden/`, `@pytest.mark.golden`)

### What landed

- **`tests/golden/`** (no `__init__.py`; module-level `pytestmark =
  pytest.mark.golden`; every module docstring says AUTHOR-GATED, change only
  with the author's sign-off, no re-record mechanism). Each test makes a fresh,
  function-scoped build from `tests/cases.py` (decision 7). 21 cases, ~4 s.
  - `test_golden_k10.py` (7, `build_k10()`): iterated rows with closedness
    (and the homotopy names, branch code `0.0`); exactly two classes, X =
    (R_1 → R_3) active and mapping over itself, U = (L_1 → L_3) inert; U rests
    on its folded loop (E8); itinerary `(R_1^1,R_3^3) (L_3,L_1^2)
    (R_3^1,R_1^3)`, not ambiguous, word [X, U⁻¹, X⁻¹], U a trivial loop;
    refinement X₁ = (R_1^1,R_3^3), X₂ = (R_1^3,R_3^1), both words
    [X₁, U⁻¹, X₂⁻¹], no unmatched member; refined matrix all ones over
    (X₁, X₂), U absent; `X.verified`, `is_reliable`.
  - `test_golden_k28.py` (9). One blast: classes {X = (R_1→R_5) active,
    U = (L_1→L_3), V = (R_2→R_4) inert}; X's two landings are singletons; word
    [X, U⁻¹, V⁻¹]; U and V trivial loops; X verified; reliable; the exterior
    class U is inert through a VIRTUAL loop (its one image loop's pair has no
    `Bridge`); 0 image bridges (E8). Two blasts: active classes in table order
    A = (R_1→R_5), B = (R_3→R_5), C = (R_1→R_3), inert U = (L_1→L_3), A first
    with min cdist ≈ 0 (`cdist_tol`); words A → [A, U⁻¹, B⁻¹], B → [C],
    C → [A, U⁻¹, A⁻¹], U → []; only A refines, A₁ = (R_1^1,R_5^3),
    A₂ = (R_1^3,R_5^1), equal words, C's refined word [A₁, U⁻¹, A₂⁻¹], every
    member matched; unrefined matrix over (A,B,C) `[[1,1,0],[0,0,1],[2,0,0]]`;
    the cut's closedness plus every P7 id-free relation (anchor start, lobe
    base = images of the iterate −1 hole's bounds, chord = reference hole =
    L_2, L_1^2 = R_5^1's ends, L_1^1/L_1^2 boundary = R_5^1.lo in the iterate
    −2 hole, holes at {0,−1,−2,−3}, R_2 = the iterate −1 hole); every class
    resolves and is verified; reliable.
  - `test_golden_period3.py` (5): p3 words (below); closed (0 image bridges);
    the same words after 4 blasts; every class resolves and is verified, no
    virtual symbol, reliable; nested (2 outer blasts): every class resolves,
    at least one `unique` walk and no `unreachable` one, image bridges exist,
    reliable.
- **`tests/helpers/names.py`**: `refined_symbol_pair`,
  `refined_words_in_names` (refined rules spelled letter-free: a split
  class's token by its child's iterated pair), `class_by_pair`.
- `source == "walk"/"trellis"` is asserted NOWHERE (decision 2, overriding the
  planner's E2 default); the k28 one-blast "resolved via the trellis path"
  fact is represented by its two singleton landings. `is_reliable` /
  `verified` now appear only in `tests/golden/` (the one remaining
  `a.is_reliable == b.is_reliable` is the session-vs-direct equivalence of
  `test_session_symbolic_dynamics._same_dynamics`, Phase 6).

### p3 element names pinned (AUTHOR ITEM E8, CLOSED: golden tier removed by the author)

Full names, recorded from today's build (CLAUDE.md states the p3 words in
letters only: `a -> b -> c -> a u^-1 w^-1`):

| CLAUDE.md | Pinned oriented pair |
|---|---|
| a (X₀) | (R_(0.0;1) → R_(1.0;5)) |
| b (X₁) | (R_(1.0;1) → R_(2.0;5)) |
| c (X₂) | (R_(2.0;1) → R_(0.0;5)) |
| u (U) | (L_(1.0;1) → L_(1.0;3)), inert |
| w (W) | (R_(1.0;2) → R_(1.0;4)), inert |

Words: X₀ → [X₁], X₁ → [X₂], X₂ → [X₀, (L_(1.0;3), L_(1.0;1)),
(R_(1.0;4), R_(1.0;2))]. The same pairs and words hold after 4 blasts and
for the inner tangle of the nested build (letter `A:` prefixed).

### Flagged for author sign-off (E8) -- CLOSED: removed by the author

**CLOSED (2026-10-05, author follow-up):** the author removed the golden tier
entirely, so none of the items below is pinned any more; see "Author
follow-up (2026-10-05): golden tier removed" at the end of this ledger.


1. **p3 names** above (`test_golden_period3.py`, constants `X0..X2`, `U`, `W`).
2. **k10 folded loop**: `test_golden_k10.py::test_k10_inert_class_rests_on_its_folded_loop`
   (U inert, exactly one member with no registered image, no image class,
   image loops = its folded loops).
3. **k28 one blast, 0 image bridges**:
   `test_golden_k28.py::test_k28_one_blast_has_no_image_bridges` (from the
   minimal-dual-graph memory note, not CLAUDE.md).
4. (Added) **which k28 one-blast class is "exterior"**: pinned as U =
   (L_1 → L_3), the inert class whose only evidence is the virtual loop (taken
   from the deleted `test_k28_has_one_active_and_two_inert_classes`); the
   other inert class (R_2 → R_4) holds the folded blast-child loop, which is
   no longer pinned (not a CLAUDE.md fact).

### Deletion ledger, Phase 5

19 collected node ids removed, 21 added (the golden tier): 1080 → 1082
(`nodeids_p5.txt`, `p5_deleted.txt`). Every deleted name and subject was
grepped in the regression-guard notes: no hit; the nested own-blast guards
(`test_nested_inner_words_are_the_period3_words`,
`test_nested_blast_order_does_not_matter`,
`test_nested_outer_blasts_leave_the_inner_bridges_alone`) are kept.

| Node id | New home | Guard checked |
|---|---|---|
| `test_bridge_class::test_k10_has_one_active_and_one_inert_class` | `golden/test_golden_k10::test_k10_classes`, `::test_k10_inert_class_rests_on_its_folded_loop` (member direction counts dropped: not a CLAUDE.md fact; fold mechanics stay in `test_loop_folds_into_its_preimage_class`) | none |
| `test_bridge_class::test_k10_inert_class_rests_on_its_folded_loop_despite_an_unresolved_member` | `golden/test_golden_k10::test_k10_inert_class_rests_on_its_folded_loop` (E8), `::test_k10_classes` (X maps over itself) | none |
| `test_bridge_class::test_k28_has_one_active_and_two_inert_classes` | `golden/test_golden_k28::test_k28_one_blast_classes_and_word`, `::test_k28_one_blast_exterior_class_is_inert_through_a_virtual_loop` | none |
| `test_session_bridge_classes::test_active_class_is_lettered_a_and_the_inert_one_is_not` | law `only_active_classes_lettered` (Phase 4); letters are presentation | none |
| `test_session_bridge_classes::test_k28_letters_only_the_anchor_class` | law `only_active_classes_lettered`; golden k28 one-blast classes | none |
| `test_element_naming::test_k10_homotopy_names`, `::test_k10_iterated_names` | `golden/test_golden_k10::test_k10_iterated_rows` (children/parent round trips stay in the synthetic naming tests) | none |
| `test_iterated_partition::test_k10_empty_stretch_cut_reads_as_expected` | `golden/test_golden_k10::test_k10_iterated_rows`; "every crossing owned closedly" is the law `iterated_unique_owner` | none |
| `test_symbolic_dynamics::test_k10_active_class_word` | `golden/test_golden_k10::test_k10_itinerary_and_word`, `::test_k10_evidence` | none |
| `test_symbolic_dynamics::test_k10_refinement_matches_every_member` | `golden/test_golden_k10::test_k10_refinement` | none |
| `test_symbolic_dynamics::test_k10_transition_matrix` | `golden/test_golden_k10::test_k10_transition_matrix` | none |
| `test_symbolic_dynamics::test_k28_one_blast_singleton_path` | `golden/test_golden_k28::test_k28_one_blast_classes_and_word` | none |
| `test_session_symbolic_dynamics::test_k28_two_blasts_structure` | `golden/test_golden_k28::test_k28_two_blasts_*` (6 tests) | none |
| `test_session_symbolic_dynamics::test_k10_symbolic_dynamics_is_built_over_the_cached_pieces` | after Phase 3b only the k10 word substring and `is_reliable` were left: golden k10; session = direct build stays in `test_k10_symbolic_dynamics_matches_a_direct_build` | none |
| `test_dual_walk::test_k10_landings_and_the_active_class_walk` | golden k10 itinerary / `verified`; law `landings_contained` | none |
| `test_higher_period_cartoon::test_p3_words_are_equivariant_under_the_orbit_shift` | `golden/test_golden_period3::test_p3_words_form_the_orbit_shift_chain`, `::test_p3_resolves_and_is_reliable` | none |
| `test_higher_period_cartoon::test_p3_words_survive_four_blasts` | `golden/test_golden_period3::test_p3_words_survive_four_blasts` (letter-free) | none |
| `test_higher_period_cartoon::test_nested_resolves_everything_with_real_walks` | `golden/test_golden_period3::test_nested_resolves_everything_after_two_outer_blasts` | none |
| `test_higher_period_cartoon::test_nested_default_words_are_pinned` | deleted, no new home: a commit snapshot (2315204), not a CLAUDE.md fact (decisions 2, 3) | hole-side-own-blast: the three own-blast guards are kept |

**Assertions removed inside surviving tests:**
- `test_topology_plotting::test_cartoon_brackets_match_closedness`: the k=10
  right-row closedness pin (golden k10 rows); the glyph = closedness property
  stays.
- `test_topology_plotting::test_plot_transition_graph_draws_every_symbol`:
  the letter `u` pin (now: no inert class's letter is drawn) and the k=10
  self-loop fact.
- `test_topology_plotting::test_plot_itinerary_table_lists_every_class`: the
  `plain["u"]` k=10 inert-loop row pin (not in the planner's list; same kind
  of fact half).
- Orphaned helpers removed: `k10_dynamics` / `_dual_and_table` / `_names`
  (`test_symbolic_dynamics`), `k10_naming` / `_homotopy`
  (`test_element_naming`), `_homotopy_names` / `_class_named`
  (`test_session_symbolic_dynamics`), `_P3_WORDS` / `_active_words`
  (`test_higher_period_cartoon`), three unused imports in `test_dual_walk`;
  module docstrings now point to `tests/golden/`.

### Verification

- Collected **1082** (`nodeids_p5.txt`).
- `1042 passed, 32 skipped, 8 xfailed` in 149 s with coverage (`p5_run.txt`);
  `-rxX` lists exactly the 8 `KNOWN_ISSUES` xfails; no XPASS.
- Coverage guard vs `cov_base.json`: OK (`cov_p5.json`).
- `pytest -m golden` selects exactly the 21 `tests/golden/` cases
  (1061 deselected); `pytest tests/golden` passes repeatedly (fresh builds,
  permuted registry ids).
- Isolation: `test_dual_graph.py::test_k10_a_different_pip_moves_the_unified_set`
  (passed alone this time; the P6 id-nondeterminism makes it flaky alone, fix
  still planned for Phase 7), `golden/test_golden_k28.py::test_k28_two_blasts_cut`,
  `golden/test_golden_period3.py::test_p3_words_survive_four_blasts`,
  `numerics/test_intersection_registry_fixes.py::test_synthetic_forwards_manifold_keys`,
  `test_stable_partition.py::test_plot_stable_partition_smoke`,
  `numerics/test_generation_and_caches.py::test_workbench_generation_bumps_on_every_mutation[trim_stable_manifolds]`:
  all pass alone.

### Deviations, Phase 5

1. **E2 resolved by decision 2, not the planner default:** the k28 one-blast
   class's "trellis path" is not pinned (`source` and `search is None` are
   asserted nowhere); its singleton landings are.
2. **Not re-pinned (not CLAUDE.md facts, snapshot detail):** the k10 active
   class's member directions `[-1,-1,1,1]`, the k28 one-blast folded
   blast-child loop and `zone_key` of the right-side classes, the k28 active
   class's single anchor member, and the k10 landing details of the deleted
   dual-walk test (anchor element lands on element 0, image on the same
   branch).
3. **Kept as is:** `test_session_symbolic_dynamics::test_describe_symbolic_dynamics`
   (Phase 3b already reduced it to "non-empty and mentions every class";
   no word substring remained).
4. **Kept the nested "at least one `unique` walk" assertion** (planner §D) —
   it reads `search.status`, not `source`.
5. **Extra derived pins:** k10's unrefined matrix `{X: {X: 2}}` (follows from
   the word) and "X maps over itself" (CLAUDE.md: "the anchor bridge maps over
   itself").
6. One fact half beyond the planner's list was removed from the plotting
   itinerary-table test (see above).

## §6 Phase 6: the facade tier (`tests/facade/`)

### What landed

- **`tests/facade/test_session_caches.py`** (50 cases, ~3.5 s): the cache
  contract as one literal `EXPECTED` table, products (`trellis`,
  `arrangement`, `bridge_classes`, `minimal_trellis`, `iterated_partition`,
  `dual_graph`, `symbolic_dynamics`) × events (`hit` = reads only, `rebuild`
  = `rebuild=True`, `generation_bump`, `partition_signature`, `pip_change`),
  each cell a fresh `build_k10()`:

  | product | hit | rebuild | gen bump | signature | pip |
  |---|---|---|---|---|---|
  | trellis, arrangement | HIT | MISS | MISS | HIT | HIT |
  | bridge_classes, minimal_trellis, iterated_partition | HIT | MISS | MISS | MISS | HIT |
  | dual_graph, symbolic_dynamics | HIT | MISS | MISS | MISS | MISS |

  - The generation bump is one `iterate_bridge` + the partition redone
    through the fan-outs; the event asserts the partition signature and the
    strong pip come back unchanged, so the generation alone is the cause.
  - The signature event re-partitions at the same generation against another
    pip and then restores the original pip (asserted), so only the signature
    moves. This fills the gap: `symbolic_dynamics` invalidated by a
    signature change (it was only tested through the pip before).
  - `rebuild` asserts only the product's own rebuild (no cascade).
  - Extra trellis rows (all MISS): growth, `iterate_bridge`,
    `rebuild_bridges`, `add_resonance_zones`, a zone's `restore` (**guard**:
    the stale-trellis gap `invalidate_trellises` used to leave), and
    `compute_intersections` on the two-fixed-point nested session (the cached
    trellis follows the swapped registry).
  - Per-selection keys for every product (`product()` vs `product(fp)`), and
    the alphabet across a rebuild (equal table, same letters, same alphabet
    size) and a re-partition (every lettered class carries the alphabet's
    letter).
- **`tests/facade/test_session_equivalence.py`** (16 cases, ~4 s): session
  result = direct build for `trellis`, `arrangement`, `homotopy_partition`,
  `bridge_classes`, `minimal_trellis`, `iterated_partition`, `dual_graph`,
  `symbolic_dynamics`, on `k10` and `nested` (2 outer blasts). The direct
  build uses only public constructors (`Trellis.from_workbench`,
  `bridge_classes`, `minimal_trellis`, `HomotopyPartition.from_results`,
  `IteratedHomotopyPartition.from_minimal`, `DualGraph`,
  `symbolic_dynamics`) over the partitions, holes and pips read off the
  per-fixed-point trellises; no `_gathered_partitions`, no `build_pieces`.
  Dynamics compared letter-free (`helpers.names`). Partition families are
  compared order-free: the session gathers in trellis-cache order, which on
  nested differs from the case's outermost-first order.
- **`tests/facade/test_session_fanouts.py`** (4 cases): no-argument fan-outs
  (`classify_strong_pips`, `compute_pseudoneighbors`, `punch_holes`,
  `partition_stable_manifold`) return a dict over both nested fixed points
  equal to each trellis's own result; single-fixed-point calls return that
  trellis's list; the describe fan-outs hold every per-trellis report under
  its fixed point; workbench attributes and drivers are reachable on the
  session (delegation smoke).

### Deletion ledger, Phase 6

40 collected node ids removed, 70 added: 1082 → 1112 (`nodeids_p6.txt`,
`p6_deleted.txt`). Every deleted name and subject was grepped in the
regression-guard notes: no hit; the one guard in this area (a restore must
drop the cached trellis) is the `restore` row.

| Node id | New home | Guard checked |
|---|---|---|
| `test_session_trellis_cache::test_same_generation_is_a_cache_hit`, `::test_reads_do_not_invalidate_the_cache` | `facade/test_session_caches::test_session_cache_table[trellis-hit]` (reads in the `hit` event) | none |
| `test_session_trellis_cache::test_rebuild_flag_forces_a_miss` | `…cache_table[trellis-rebuild]` | none |
| `test_session_trellis_cache::test_growth_invalidates_the_cache`, `::test_iterate_bridge_…`, `::test_rebuild_bridges_…`, `::test_add_resonance_zones_…`, `::test_restore_invalidates_the_cache` | `…test_trellis_misses_after_every_mutation_path[growth / iterate_bridge / rebuild_bridges / add_resonance_zones / restore]` | restore = the 2026-07 stale-trellis gap: kept as a row |
| `test_session_trellis_cache::test_compute_intersections_invalidates_the_cache`, `test_session_strong_pips::test_trellis_auto_rebuilds_after_registry_change` | `…test_trellis_misses_after_a_recompute_on_the_nested_session` (registry identity kept) | none |
| `test_arrangement::test_arrangement_is_cached_by_generation` | `…cache_table[arrangement-hit / arrangement-generation_bump]` | none |
| `test_session_bridge_classes::test_cache_hit_returns_the_same_object`, `::test_reads_do_not_invalidate_the_cache`, `::test_workbench_mutation_invalidates_the_cache`, `::test_repartitioning_at_the_same_generation_invalidates_the_cache` | `…cache_table[bridge_classes-*]`; the letter half → `…test_bridge_class_letters_survive_rebuilds_and_repartitions` | none |
| `test_session_bridge_classes::test_cache_is_kept_per_fixed_point_selection` | `…test_cache_is_kept_per_fixed_point_selection[*]` (all products) | none |
| `test_session_bridge_classes::test_letters_are_stable_across_a_rebuild` | `…test_bridge_class_letters_survive_rebuilds_and_repartitions` (+ `…cache_table[bridge_classes-rebuild]`) | none |
| `test_session_bridge_classes::test_k10_session_result_matches_a_direct_call`, `::test_p3_session_result_matches_a_direct_call` | `facade/test_session_equivalence::…[bridge_classes-k10 / -nested]` (nested now the 2-outer-blast case, not the unblasted 1e-4 one) | none |
| `test_session_dual_graph::test_cache_hits_return_the_same_objects`, `::test_rebuild_flag_forces_fresh_equivalent_objects` (cascade dropped), `::test_workbench_mutation_invalidates_the_caches`, `::test_repartitioning_at_the_same_generation_invalidates_the_caches`, `::test_pip_change_at_the_same_generation_and_partition_invalidates_the_graph` | `…cache_table[minimal_trellis-* / iterated_partition-* / dual_graph-*]`; "a different pip moves the unified set" stays in `test_dual_graph::test_k10_a_different_pip_moves_the_unified_set` | none |
| `test_session_dual_graph::test_dual_graph_is_built_over_the_cached_pieces` | deleted (private `_partition_signature` / `_gathered_partitions`, cache identity); content equality → equivalence `[homotopy_partition-*]`, `[iterated_partition-*]`, `[dual_graph-*]` | none |
| `test_session_dual_graph::test_k10_dual_graph_matches_a_direct_build`, `::test_p3_dual_graph_matches_a_direct_build` | equivalence `[dual_graph-k10 / -nested]` | none |
| `test_session_symbolic_dynamics::test_cache_hits_return_the_same_object`, `::test_rebuild_flag_forces_a_fresh_equivalent_object` (cascade dropped), `::test_workbench_mutation_invalidates_the_cache`, `::test_pip_change_at_the_same_generation_invalidates_the_cache` | `…cache_table[symbolic_dynamics-*]` (+ the new signature cell); the "itineraries stay even" half is the law `itineraries_even` | none |
| `test_session_symbolic_dynamics::test_k10_symbolic_dynamics_matches_a_direct_build` | equivalence `[symbolic_dynamics-k10 / -nested]` (letter-free; `is_reliable` no longer compared outside golden) | none |
| `test_partition_family::test_from_results_signature_matches_the_session_signature` | equivalence `[homotopy_partition-*]` (public `signature()`, no `_partition_signature` / `_gathered_partitions`) | none |
| `test_session_strong_pips::test_classify_strong_pips_all_fixed_points`, `::test_classify_single_fixed_point_returns_list` | `facade/test_session_fanouts::test_no_argument_fanouts_…`, `::test_single_fixed_point_fanouts_…` | none |
| `test_session_pseudoneighbors::test_compute_all_fixed_points_returns_dict`, `::test_compute_single_fixed_point_returns_list`, `::test_punch_and_partition_fan_out_to_every_fixed_point`, `::test_punch_and_partition_single_fixed_point_return_lists`, `::test_describe_fan_outs_mention_every_fixed_point` | `facade/test_session_fanouts` (dict shapes now over the two nested fixed points instead of the single k10 one) | none |

**Assertions removed inside surviving tests:**
- `numerics/test_grow_until::test_session_exposes_the_drivers`: the cache half
  (`trellis(fp) is not stale` after growth → the `growth` row) and the
  `callable(getattr(session, …))` loop (→ the fan-out delegation smoke).
- Orphaned helpers and imports removed (`_direct_table`, `_strip`,
  `_same_graph`, `_rule_texts`, `_refined_rule_texts`, `_same_dynamics`,
  `build_pieces` in the two session files, unused `DualGraph`,
  `HomotopyPartition`, `symbolic_dynamics`, `numpy`, `pytest`,
  `TangleSession` imports); module docstrings point to `tests/facade/`.

### Verification

- Collected **1112** (`nodeids_p6.txt`).
- `1072 passed, 32 skipped, 8 xfailed` in 150 s with coverage (`p6_run.txt`); `-rxX`
  lists exactly the 8 `KNOWN_ISSUES` xfails; no XPASS.
- Coverage guard vs `cov_base.json`: OK (`cov_p6.json`). The first run
  flagged `BridgeClassTable.__eq__` (lines 445/447), covered only by the
  deleted `test_letters_are_stable_across_a_rebuild`; the letters test now
  asserts `rebuilt == first`.
- Isolation (each alone): `test_dual_graph.py::test_k10_a_different_pip_moves_the_unified_set`,
  `facade/test_session_caches.py::test_session_cache_table[symbolic_dynamics-partition_signature]`,
  `facade/test_session_caches.py::test_trellis_misses_after_every_mutation_path[restore]`,
  `facade/test_session_equivalence.py::test_session_product_equals_a_direct_build[dual_graph-nested]`,
  `facade/test_session_fanouts.py::test_no_argument_fanouts_return_one_entry_per_fixed_point`,
  `test_session_bridge_classes.py::test_no_partitions_raises_and_warns`,
  `numerics/test_grow_until.py::test_session_exposes_the_drivers`: all pass.

### Deviations, Phase 6

1. **Wall time not reduced** (~151 s vs ~149 s with coverage). The deleted
   tests were cheap (k10 builds take ~0.05 s); the table adds 50 cells and
   the equivalence suite 16, each a fresh function-scoped build (decision 7).
   The facade tier costs ~8 s in all.
2. **The generation-bump event is cleaner than planned:** `iterate_bridge` on
   k10 registers no new crossing, so after re-partitioning the signature and
   pip are asserted unchanged, isolating the generation as the cause.
3. **The signature event restores the original pip** so the `pip` component
   of the dual-graph/dynamics key stays put; the old repartition tests moved
   both at once.
4. **Fan-out shape tests moved in this phase** (planner §B puts
   `test_session_pseudoneighbors` / `test_session_strong_pips` fan-outs in
   `F/fanouts`); their plot halves stay for Phase 8.
5. `homotopy_partition` added as an eighth equivalence product (it replaces
   the deleted `from_results` signature test); the partition families are
   compared order-free (see above).
6. Kept `test_session_symbolic_dynamics::test_p3_symbolic_dynamics_smoke` and
   the plot-delegate mocks: not cache tests; they are Phase 7/8 items.

## §7a Phase 7a: unit consolidation, numerics

### What changed

- **Linked-list points** (`tests/test_point.py`, `tests/test_branch_point.py`
  merged into it and deleted): the iterate tests are parametrized over the two
  point kinds (`point`, `branch_point`): 8 → 4 functions. The `Point(num_branches)`
  misuse (it set `x = 2`) is gone; the tautological `test_get_point`
  (`a.all() == b.all()`) is folded into `test_point_creation` /
  `test_branch_point_creation` as `assert_array_equal(get_point(), [x, y])`.
- **Geometry** (`test_geometry.py`): the three area tests and the winding test
  → `test_signed_area_of_a_known_shape_flips_with_winding[square|triangle|L]`;
  concavity + winding + closed-ring containment →
  `test_point_in_polygon_respects_the_concavity[ccw|cw|closed_ring]`;
  interior-point winding merged into `test_interior_point_is_inside` (`L_cw`,
  `crescent_cw` added). **Strengthened:**
  `test_oriented_bridge_polyline_runs_in_the_dynamical_direction` now asserts
  `poly[0]` is the low-cdist end (nearer the lower-cdist registered crossing,
  equal to the lower-cdist end node), non-vacuously.
- **Registry generation** (`test_generation_and_caches.py`): four bump tests →
  `test_registry_generation_bumps_on_every_mutation[add|add_synthetic|register_iterate|reindex_from]`
  (the two multi-step mutators assert their own step's bump). Memo `is`
  asserts dropped (decision: no cache-identity pins outside the facade):
  `test_bridge_point_array_is_memoised` → `…_agrees_with_a_walk`,
  `test_branch_position_map_is_memoised` → `…_is_stable_across_calls`.
- **Growth drivers** (`test_grow_until.py`, `test_workbench_bugfixes.py`):
  - caps → `test_driver_raises_at_the_cap_after_growing[grow_until|grow_until_intersection|grow_until_arclength]`
    (each asserts the cap is reached by growing);
  - rejections → `test_driver_rejects_a_bad_request_before_growing[...]` (6
    rows: empty grow set, uninitialized manifold, unknown id, unknown
    direction, empty faces id set, uninitialized branch = the
    `honours_branch_index` **guard**, audit 2026-07), each asserting no point
    moved;
  - the anchor faces-closed limitation moved to the module's Dev Notes
    (decision 3, provisional).
- **`Intersection.fixed_points`** (`test_intersection_registry_fixes.py`): five
  tests → `test_fixed_points_table[only_b|only_a|same_deduped|distinct_ab_order|no_keys]`
  (**guard**: the IndexError, codebase-audit-2026-09); the two
  `Intersection.synthetic` tests → one, keeping the label-slot assertion
  (**guard**); the local `_stub_fixed_point` replaced by
  `helpers.fakes.bare_fixed_point` (a real `FixedPoint`).
- **One-rule guard upgraded** (`test_map_step_and_graph.py`, open item E9):
  `test_workbench_has_no_private_key_advance` (a `hasattr` tombstone) →
  `test_every_map_step_reaches_fixed_point_advance_key`: a counting spy on
  `FixedPoint.advance_key`, asserting iterate inference, `iterate_bridge` and
  `Trellis.image_cdist` (40 steps, off the table, `from_table=False`) each
  reach it. `test_per_step_beta_is_the_k_th_root_only_without_inversion` is
  merged as the explicit inversion branch of
  `test_per_step_beta_unstable_is_the_period_th_root_of_the_eigenvalue`
  (**guard**: the k_value-root bug).
- **Strengthened:**
  - `test_kevin_way_period_three_orbit_chain[unstable|stable]` asserts the
    chain order read off the iterate links (unstable 0, 1, 2; stable 0, 2, 1)
    at strictly increasing cdist (**guard**: kevin-way ordering,
    codebase-audit-2026-09).
  - `test_period_three_orbit_at_k_two_is_accepted`: each orbit point maps to
    the next, `f^3(p) = p`, the orbit is not a fixed point, and
    `lambda_u * lambda_s = 1` per orbit point (**guard**: solver root cause).
  - Initializer alpha: `test_alpha_matches_distance_ratio` (tautological for
    k = 1) → `test_alpha_is_the_per_step_factor[k10|inversion|p3 × stability]`:
    the two real seed points are `k_value` iterate links apart with cdist ratio
    `alpha ** k_value`, every fictitious link carries `alpha`, and
    `alpha ≈ per_step_beta` (reciprocal on the stable side) to `SCALING_RTOL`.
  - `test_trim_stable_manifolds_reads_the_registry` asserts the trim removed
    points and the tail is the FIRST node at or past the outermost crossing.
  - `test_min_separation_drops_close_bridges` asserts the distances: every
    kept child's interior is ≥ `min_separation` from the curve accumulated
    before it, every dropped one < (**guard**: zig-zag).
  - `test_blast_recognizes_already_known_bridges`: no bridge is a parent
    twice, no already-known child is a later parent (gap 6).
  - `test_low_stretch_growth_keeps_cdist_injective` absorbs the
    fundamental-segment case (still the only strict cdist test).
- **Replaced:** `test_refined_cdist_is_mean_of_neighbours` →
  `test_refined_cdist_lies_strictly_between_its_neighbours` (decision 3; the
  mean rule is in the module's Dev Notes). It still calls the refinement kernel
  `_get_refined_point` (no public single-point route).

### Deletion ledger, Phase 7a

54 collected node ids removed, 50 added: 1112 → 1108 (`nodeids_p7a.txt`,
`p7a_deleted.txt`, `p7a_added.txt`). Every deleted subject was grepped in the
regression-guard notes; the hits (mixed signs, `_advance_key_forward`,
kevin-way, `Intersection.fixed_points`, label slot) are all kept as noted.

| Node id | New home / reason | Guard checked |
|---|---|---|
| `test_point::test_get_point`, `::test_insert_next_iterate`, `::test_insert_prev_iterate`, `::test_insert_next_iterate_error`, `::test_insert_prev_iterate_error`; `test_branch_point::*` (9) | `test_point::test_point_creation` / `test_branch_point_creation` / `test_branch_point_insert_point_*` / `test_insert_{next,prev}_iterate[*]` / `test_insert_{next,prev}_iterate_twice_raises[*]` | none |
| `test_geometry::test_signed_area_of_a_ccw_square_is_positive_one`, `::test_signed_area_flips_with_winding`, `::test_signed_area_of_a_triangle`, `::test_signed_area_of_a_concave_polygon` | `::test_signed_area_of_a_known_shape_flips_with_winding[*]` | none |
| `test_geometry::test_point_in_polygon_respects_the_concavity` (unparametrized), `::test_point_in_polygon_is_winding_independent`, `::test_point_in_polygon_accepts_an_explicitly_closed_ring` | `::test_point_in_polygon_respects_the_concavity[ccw|cw|closed_ring]` | none |
| `test_geometry::test_interior_point_is_winding_independent` | `::test_interior_point_is_inside[L_cw|crescent_cw]` | none |
| `test_generation_and_caches::test_registry_generation_bumps_on_{add,add_synthetic,register_iterate,reindex_from}` | `::test_registry_generation_bumps_on_every_mutation[*]` | none |
| `test_generation_and_caches::test_bridge_point_array_is_memoised`, `::test_branch_position_map_is_memoised` | renamed `…_agrees_with_a_walk` / `…_is_stable_across_calls`, memo `is` dropped | none |
| `test_grow_until::test_grow_until_raises_at_the_cap`, `test_workbench_bugfixes::test_grow_until_intersection_raises_at_the_cap`, `::test_grow_until_arclength_raises_at_the_cap` | `test_grow_until::test_driver_raises_at_the_cap_after_growing[*]` | none |
| `test_grow_until::test_grow_until_rejects_an_empty_grow_set`, `::test_grow_until_requires_an_initialized_manifold`, `::test_iterates_closed_rejects_an_unknown_id`, `::test_iterates_closed_rejects_an_unknown_direction`, `::test_faces_closed_rejects_an_empty_id_set`, `test_workbench_bugfixes::test_grow_until_intersection_honours_branch_index` | `test_grow_until::test_driver_rejects_a_bad_request_before_growing[*]` | `honours_branch_index` (audit-polish-2026-07): kept as the `uninitialized_branch` row |
| `test_grow_until::test_faces_closed_raises_at_the_cap_for_the_anchor` | deleted (decision 3/4: provisional anchor limitation → module Dev Notes; the cap is the parametrized cap test) | regions-refactor mentions `grow_until_faces_closed` as a feature only |
| `test_workbench_bugfixes::test_trim_stable_manifolds_cuts_just_past_the_outermost_crossing` | deleted (private R-tree oracle: `Tangle._intersecting_segments`, `_manifold_segs`, `_seg_lookup`); the property → strengthened `test_single_source_of_truth::test_trim_stable_manifolds_reads_the_registry` | none |
| `test_intersection_registry_fixes::test_fixed_points_*` (5) | `::test_fixed_points_table[*]` | IndexError (codebase-audit-2026-09): every row kept |
| `test_intersection_registry_fixes::test_synthetic_does_not_put_the_label_in_the_id_slot`, `::test_synthetic_forwards_manifold_keys` | `::test_synthetic_keeps_the_label_out_of_the_id_slot_and_forwards_keys` | label slot (codebase-audit-2026-09): kept |
| `test_invariant_helpers::test_fundamental_segments_have_injective_cdist` | `::test_low_stretch_growth_keeps_cdist_injective` | cdist-strict-monotonicity-fix: strict stays on low stretch only |
| `test_inversion_fixture::test_mixed_eigenvalue_signs_are_rejected`, `test_fixed_point::test_set_k_value_rejects_disagreeing_eigenvalue_signs` | deleted (decision 9/14): the b = -1 placeholder `test_law_case_sanity::test_orientation_reversing_case_builds[orientation_reversing]`, `xfail(strict=True, raises=ValueError)`, encodes it and still executes the raise (coverage guard OK) | regions-refactor ("`set_k_value` raises on mixed signs"): covered by the placeholder |
| `test_map_step_and_graph::test_per_step_beta_is_the_k_th_root_only_without_inversion` | merged into `::test_per_step_beta_unstable_is_the_period_th_root_of_the_eigenvalue` | k_value-root bug: kept |
| `test_map_step_and_graph::test_workbench_has_no_private_key_advance` | `::test_every_map_step_reaches_fixed_point_advance_key` (spy) | `_advance_key_forward` (codebase-audit, numerics-test-suite): behavioural guard now |
| `test_initializer_cdist::test_alpha_matches_distance_ratio[*]` | `::test_alpha_is_the_per_step_factor[*]` | none |
| `test_refinement::test_refined_cdist_is_mean_of_neighbours` | `::test_refined_cdist_lies_strictly_between_its_neighbours` | none |
| `test_manifold_initializer::test_kevin_way_period_three_orbit_chain[unstable|stable]` | same name, new ids `[unstable|stable]` from the `(stability, chain)` parametrization | kevin-way (codebase-audit): strengthened |

### Verification

- Collected **1108** (`nodeids_p7a.txt`).
- `1068 passed, 32 skipped, 8 xfailed` in 161 s with coverage (`p7a_run.txt`);
  `-rxX` lists exactly the 8 `KNOWN_ISSUES` xfails; no XPASS.
- Coverage guard vs `cov_base.json`: OK (`cov_p7a.json`).
- Isolation (each alone): `test_dual_graph.py::test_k10_a_different_pip_moves_the_unified_set`,
  `numerics/test_map_step_and_graph.py::test_every_map_step_reaches_fixed_point_advance_key`,
  `numerics/test_initializer_cdist.py::test_alpha_is_the_per_step_factor[p3-stable]`,
  `numerics/test_grow_until.py::test_driver_rejects_a_bad_request_before_growing[grow_until_intersection-uninitialized_branch]`,
  `numerics/test_blast_proximity_guard.py::test_min_separation_drops_close_bridges`,
  `test_manifold_initializer.py::test_kevin_way_period_three_orbit_chain[stable]`,
  `numerics/test_generation_and_caches.py::test_workbench_generation_bumps_on_every_mutation[compute_intersections]`:
  all pass.

### Deviations, Phase 7a

1. **Spy reaches inference through `infer_iterates`, not
   `infer_iterate_table`.** `infer_iterate_table` only visits endpoints of
   iterated bridges, and `iterate_bridge` already links those, so on any real
   build it predicts nothing (0 spy calls). The spy test builds k10 with
   `compute_intersections(infer_iterates=False)` and calls `infer_iterates()`,
   which runs the same `_image_prediction` step.
2. **Scope split.** Items in the Phase 7 list that live in topology/loom
   files (fakes rewiring and `FakeFixedPoint` retirement, `build_pieces`,
   dual-graph assembly copies, `k10_session` shadows, `_define_zone`, the
   dual-walk / strong-pip / blast-error merges, the pseudoneighbor
   strengthenings, `test_forward_unstable_branch_cycle_matches_orbit_order`,
   the snap tripwire, the half-edge and `_region_key` tests) are Phase 7b.
   Two blast tests in `tests/numerics/` (proximity guard, already-known) were
   done here because they live in the numerics directory.
3. **Files not renamed.** The merged Point/BranchPoint tests stay in
   `tests/test_point.py` (planned destination `unit/numerics/test_linked_list.py`
   is a Phase 10 move).
4. **Finding, not fixed (regression tier, Phase 9):**
   `regression/test_high_stretch_period3_growth.py` is VACUOUS on the current
   code: four unstable growth steps of the k = 2.1 period-3 orbit leave 7
   nodes per branch with adjacent cdist gaps of ~3e15 ULP (no refinement
   happens), so the run never reaches the near-ULP regime it guards. Planner
   §B's "add a non-vacuity check" would fail today; the parameters need
   re-deriving (author item). `regression/test_cdist_collision_growth.py`
   merge and the blast-monotonicity strengthening were also left for that
   phase.
5. **Dissipative `per_step_beta` test kept** (E6 lists it as a deletion
   candidate pending the author's answer; until then current behaviour stays).

## §7b Phase 7b: unit consolidation, topology and loom

### What changed

- **Fakes delegate to the real rules** (`tests/helpers/fakes.py`):
  `FakeFixedPoint` is gone; every synthetic fixed point is
  `bare_fixed_point(...)`, a real `FixedPoint` (new `coordinates=` keyword for
  hand-built `TrellisBranch` anchors). `FakeTrellis` now borrows
  `Trellis.image_cdist`, `_scaled_image_cdist` and `_keyed_cdist` and supplies
  data only (a `_FakeRegistry` for `registry[id]` / `cdist_tol`; crossings carry
  no unstable key, so an unstable lookup raises as the real trellis does).
  `make_result` requires explicit `owners` (the derived-ownership default is
  removed); `test_dual_walk.py` spells every family's owners next to its
  intervals (`Family = (specs, owners)`).
- **Local fake fixed points retired** in `test_bridge_class`, `test_arrangement`,
  `test_arrangement_sparse`, `test_symbolic_dynamics`, `test_element_naming`
  (`_FakeFixedPoint` → `bare_fixed_point`), and the three real-but-duplicated
  `_fixed_point(period, lambda_u)` helpers in `test_pseudoneighbor`,
  `test_strong_pip_periodic_point`, `test_stable_partition` now delegate to it.
- **`minimal_helpers.build_pieces` retired** (all call sites in
  `test_dual_graph`, `test_iterated_partition`, `test_minimal_trellis`,
  `test_element_naming`, `test_higher_period_cartoon`, `test_topology_plotting`)
  in favour of `session.minimal_trellis()`, `homotopy_partition()`,
  `iterated_partition()`, `bridge_classes()`, `dual_graph()`,
  `symbolic_dynamics()`; `session._gathered_partitions()` in
  `test_partition_family` → `session.homotopy_partition()`.
  `tests/minimal_helpers.py` is deleted now (it had no users left, and the
  phase's grep for `build_pieces|_gathered_partitions` must be empty).
  `grep -rn "build_pieces\|_gathered_partitions" tests/` → empty.
- **Dual-graph assembly copies retired**: `_k10_graph` / `_p3_graph`
  (`test_dual_graph`), `_k10_dynamics` (`test_higher_period_cartoon`),
  `_k10_dual_graph` / `_k10_dynamics` (`test_topology_plotting`) all read the
  session's own products. Direct `DualGraph(...)` / `from_minimal(...)` /
  `minimal_trellis(...)` construction stays only where the test varies an
  input (other pips, a homotopy-only family, a patched table).
  `test_dual_graph` no longer reads `dual._nodes_of_edge` (edges are grouped
  from the public `stable_nodes`).
- **Session shadows retired**: the local `k10_session` of
  `test_loom_blast_restore` (→ the conftest fixture), `henon_session` in
  `test_session_pseudoneighbors` (→ `k10_session`) and `test_topology_plotting`
  (→ `build_k10(through="bridges")`), `k10_zone_session` in
  `test_resonance_zone_region` (→ `build_k10(through="pips")` + the zone), the
  hand-grown 10-step session in `test_session_bridge_classes`
  (→ `build_k10(unstable_steps=10, through="pips")`), `p3_partitioned_c` in
  `test_image_cdist` (→ conftest `p3_partitioned`), and the session-scoped
  `p3_built` / `nested_built` fixtures plus the direct `henon_cases` imports of
  `test_higher_period_cartoon` and `test_stable_partition_period3`
  (→ function-scoped `cases.build_period3` / `cases.build_nested`).
- **`_define_zone` ×2 retired** (`test_loom_blast_restore`,
  `facade/test_session_caches`) → `tests/helpers/zones.py::define_inner_zone`.
- **Merges / parametrizations**: dual-walk trivial cases (×2 → 1 ×2), the
  entry-then-exit crossing and its mirror (×2 → 1 ×2), the three table landings
  (×3 → 1 ×3); strong-pip disqualifier ownership (heteroclinic / homoclinic /
  no unstable key, ×3 → 1 ×3; the table-linked collision guard untouched);
  blast error handling (×3 → 1 ×3, **guard** "assertions never swallowed"
  kept as the `assertion_escapes_lenient` row);
  `test_nested_outer_blasts_leave_the_inner_bridges_alone` merged into
  `test_nested_inner_words_are_the_period3_words` (same words AND same inner
  bridge count; **guard** hole-side-own-blast-2026-10-02).
- **Strengthened**: `test_pair_on_different_unstable_branches_rejected` uses a
  clean fixture (the standard window, middle crossing on the saddle's other
  unstable branch, plus a same-branch control); trajectory extension split
  into `…forward_stops_at_the_end_of_the_chain` (asserts iterates 1, 2) and
  `test_trajectory_extension_deduplicates[against_a_reference|against_another_trajectory]`
  (REAL duplicates: two references on one trajectory).
- **Private access rewritten through the public API**: the strong-pip
  uncut-branch warning through `compute_pseudoneighbors` (was
  `_strong_pip_cuts`; E6 Q13: level/logger/count only); the inert letter
  series through `inert_letters` (was `SD._inert_letter` / `_active_letter`);
  `test_image_cdist` drops its `trellis._scaled_image_cdist` equality (the test
  recomputes the law itself); the `not hasattr(zone, "boundary_intersection_id")`
  tombstone line is dropped; the zone midpoint classification additionally
  asserts the public `session.classify_bridge(bridge) is zone`.
- **Coverage-gap tests added in the same commit** (the deleted tests were the
  only route to these lines): `test_unlettered_active_classes_get_fallback_letters`
  (`_active_letter` via `symbolic_dynamics` over an unlettered table),
  `test_unresolved_class_carries_a_reason_and_warns` (the singleton landing
  with no registered member chain, on the unblasted nested tangle; replaces the
  p3 smoke test), `test_orientation_without_a_jacobian[*]` (the E6 orientation
  default, DEBUG count only; it was covered by the old duck-typed fake fixed
  points that lacked stable eigenvalues).

### Deletion ledger, Phase 7b

24 collected node ids removed, 27 added: 1108 → 1111 (`nodeids_p7b.txt`,
`p7b_deleted.txt`, `p7b_added.txt`). Every deleted subject was grepped in the
regression-guard notes; hits are noted.

| Node id | New home / reason | Guard checked |
|---|---|---|
| `test_dual_walk::test_start_equals_goal_is_trivial`, `::test_start_face_equal_to_goal_face_is_trivial` | `::test_trivial_walks_stay_on_the_start_face[start_is_goal|start_face_is_goal_face]` | none |
| `test_dual_walk::test_single_crossing_records_entry_then_exit`, `::test_crossing_from_the_right_side_records_right_then_left` | `::test_single_crossing_records_entry_then_exit[left|right]` (the mirror now asserts every step field) | none |
| `test_dual_walk::test_contained_landing_from_the_table`, `::test_unbounded_low_end_maps_to_zero`, `::test_outermost_element_lands_by_the_table` | `::test_landing_from_the_table_is_contained[interior|unbounded_low_end_maps_to_zero|outermost]` | none |
| `test_strong_pip_periodic_point::test_heteroclinic_crossing_does_not_disqualify_strong_pip`, `::test_homoclinic_crossing_still_disqualifies_strong_pip`, `::test_missing_unstable_key_disqualifier_is_kept` | `::test_only_the_own_tangle_disqualifies_a_strong_pip[*]` | regions-refactor (heteroclinic): kept as a row |
| `test_loom_blast_restore::test_blast_zone_propagates_assertion_error`, `::test_blast_zone_skips_value_error_with_warning`, `::test_blast_zone_reraises_value_error_when_strict` | `::test_blast_zone_error_handling[*]` | pseudoneighbor-partition "swallowed" hit is about cdist slack, not this; the assertion-escape guard is the first row |
| `test_higher_period_cartoon::test_nested_outer_blasts_leave_the_inner_bridges_alone` | merged into `::test_nested_inner_words_are_the_period3_words` | hole-side-own-blast-2026-10-02: kept (bridge-count assertion) |
| `test_pseudoneighbor::test_pair_on_different_unstable_branches_rejected` | `::…rejected[same_branch|other_branch]` (clean fixture) | none |
| `test_pseudoneighbor::test_trajectory_extension_forward_with_dedup` | `::test_trajectory_extension_forward_stops_at_the_end_of_the_chain` + `::test_trajectory_extension_deduplicates[*]` | dedup hits (codebase-audit `_edge_seen`, pseudoneighbor-partition region dedup) are other code; nothing lost |
| `test_pseudoneighbor::test_strong_pip_cuts_warn_when_branches_are_left_uncut`, `::test_strong_pip_cuts_are_silent_when_every_branch_is_cut` | `::test_strong_pip_cuts_warn_when_branches_are_left_uncut[period_1_cut|period_3_uncut]` (public `compute_pseudoneighbors`) | none |
| `test_pseudoneighbor::test_forward_unstable_branch_cycle_matches_orbit_order` | deleted (decision 13: dead-API wrapper; `branch_cycle` is tested in numerics) | pseudoneighbor-partition / regions-refactor mention `branch_cycle` as a feature only |
| `test_stable_partition::test_region_key_separates_the_two_sides_of_one_bridge` | `::test_both_holes_of_a_singleton_bridge_survive_propagation[k28_one_blast|p3]` (real data: a bridge with a propagated hole on either side keeps both) | pseudoneighbor-partition "SINGLETON BRIDGE" (user-confirmed): kept as the property |
| `test_arrangement_sparse::test_sparse_kept_endpoints_are_degree_three_without_unstable_stubs` | deleted (private half-edges; no public degree query); "no unstable stub" → `_dangling_ends == 2` in `::test_sparse_faces_close_and_euler_holds` | regions-refactor mentions half-edges as design only |
| `test_image_cdist::test_the_snap_does_not_fire_on_the_default_path` | deleted (E6: private fixture tripwire on `_snap_to_partition_boundary`); module docstring says it is deliberately unpinned | none ("snap" hits are "snapshot") |
| `test_session_symbolic_dynamics::test_p3_symbolic_dynamics_smoke` | deleted (planner §B); its laws run on nested in `invariants/test_law_symbolic.py`; the unresolved path → `test_symbolic_dynamics::test_unresolved_class_carries_a_reason_and_warns` | none |
| `test_dual_graph::test_k10_stable_node_names_are_the_expected_iterated_names` | deleted (planner §B: scattered fact pin of the k10 names; pinned once in `golden/test_golden_k10.py::test_k10_iterated_rows`) | none |

Unchanged ids whose body changed: `test_symbolic_dynamics::test_inert_letter_series`
(public `inert_letters`), `test_image_cdist::*` (fixture `p3_partitioned`),
every rewritten `build_pieces` caller.

### Verification

- Collected **1111** (`nodeids_p7b.txt`).
- `1071 passed, 32 skipped, 8 xfailed` in 169 s with coverage (`p7b_run.txt`);
  `-rxX` lists exactly the 8 `KNOWN_ISSUES` xfails; no XPASS.
- Coverage guard vs `cov_base.json`: OK (`cov_p7b.json`). The first run
  flagged `SymbolicDynamics` 337 / 1122 / 1144 / 1148–1149 and `Trellis`
  316 / 320; the three gap tests above close them.
- `grep -rn "build_pieces\|_gathered_partitions" tests/` → empty.
- Isolation (each alone): `test_dual_graph.py::test_k10_a_different_pip_moves_the_unified_set`,
  `test_dual_walk.py::test_single_crossing_records_entry_then_exit[right]`,
  `test_loom_blast_restore.py::test_blast_zone_error_handling[value_error_skipped_lenient]`,
  `test_pseudoneighbor.py::test_trajectory_extension_deduplicates[against_another_trajectory]`,
  `test_stable_partition.py::test_both_holes_of_a_singleton_bridge_survive_propagation[p3]`,
  `test_symbolic_dynamics.py::test_unresolved_class_carries_a_reason_and_warns`,
  `test_higher_period_cartoon.py::test_nested_inner_words_are_the_period3_words`,
  `test_topology_plotting.py::test_plot_minimal_trellis_draws_every_kept_bridge`:
  all pass.

### Deviations, Phase 7b

1. **`tests/minimal_helpers.py` deleted in this phase** instead of Phase 10:
   it had no users left, and the phase's grep gate covers the file itself.
   Phase 10 still deletes the two remaining shims.
2. **`_region_key` rewrite is on real data**, not synthetic. On k28 one blast
   and p3 the both-sided bridges carry holes of two DIFFERENT origins (one
   origin on both sides would break I1), so the test pins the public outcome
   (both holes survive propagation) rather than the private key.
3. **More shadows retired than the plan named** (`henon_session` ×2,
   `k10_zone_session`, the 10-step session, `p3_partitioned_c`, the
   session-scoped `p3_built`/`nested_built`, two direct `henon_cases`
   imports): all were the same recipes as `cases.py` builders (decision 7).
   `test_session_pseudoneighbors` and `test_topology_plotting` now run on the
   recipe with the iterate table inferred (P6: interchangeable).
4. **Gap tests added** (fallback letters, unresolved class, orientation
   default) because the deleted / rewritten tests were the only route to those
   lines; each asserts behaviour, levels and loggers only.
5. **Private access left for later phases**: `_bridge_side_of`
   (`test_stable_partition`, no public single-bridge route; candidate for a
   kernel exception or a Phase 9 rewrite), `_dangling_ends` reading
   `arrangement._nodes` (no public dangling-end count; Phase 3b kept it),
   `_is_forward_beyond_fundamental` (firm-only since Phase 2), and the plotting
   privates (`_side_of`, `_anchorward_look`, `_cartoon_name`: Phase 8). The
   `_empty_stretches` and `_resolve_inertness` uses are allowed kernels.
6. `test_k10_stable_node_label_joins_the_names` and the plot-delegate mocks
   stay for Phase 8 (they are on its delete list).

## §8 Phase 8: the plotting tier (`tests/plotting/`)

### What landed

- **`tests/plotting/test_plot_smoke.py`**: every plotter and session delegate
  runs on Agg and returns its type: the five `Trellis.plot_*` on real k=10
  data; the session `plot_*` fan-outs on the two-fixed-point nested session
  (one handle per fixed point, one for a single fixed point); the
  compute-on-demand `plot_pseudoneighbors`; `plot_stable_partition` with row
  labels and one element label per interval; `plot_dual_graph` under every
  label/clip option (canvas drawn, so the mathtext must parse);
  `plot_minimal_trellis`, `plot_walk` (style overrides), `plot_transition_graph`
  and `plot_bridges_by_class` (refined / unrefined); `plot_itinerary_table`
  (one row per class; plain, unrefined, a column subset on a wide axes, the
  font-shrink path on a 2.5 in axes); `itinerary_table_rows` column subset in
  the caller's order; the cartoon under every option (anchor left, recoloured,
  outside labels, uniform width, class-coloured, circle); every legend helper
  returns a non-empty list of artists; the session delegates (forwarded
  keywords checked by their EFFECT, e.g. `refined=False` gives the unrefined
  graph's nodes, not by a monkeypatched spy); the cartoon toggles; and one
  `ValueError` table (type only) for reserved style keywords and unknown
  option values (12 rows).
- **`tests/plotting/test_plot_topology.py`**: face point inside its region;
  face-point copy (regression); unbounded node outside the midpoint bbox; side
  nodes on their side through public `stable_node_point` + the arc polyline +
  the allowed `_side_of` kernel (the private `_anchorward_look` is gone);
  every stable node drawn at its point; a REAL walk threads face -> node ->
  face and ends in one arrowhead; trivial / empty walks (public `Walk`
  dataclass instead of the `_FakeWalk` shim); every kept minimal-trellis bridge
  drawn; transition-graph nodes/edges = the dynamics' own, all active;
  `plot_bridges_by_class` draws every classed bridge (and drops exactly the
  inert ones with `show_inert=False`); an unmatched member is still drawn and
  keeps a colour of its own; `bridge_ids` restricts to the minimal trellis
  (nested); `class_colors` never repeats (k10 and nested, refined and not);
  cartoon one node per element per side, ordinal coordinate strictly
  increasing, brackets = closedness, kept bridges drawn and walks element to
  element; circle interior side inside the circle (k10 default and flipped);
  p3 circle keeps the zone's ring order and direction (cyclic order of the
  arc anchors, no angle pins); nested inner circle inside the outer one and
  its keep-out disc, with the enclosing lobes reported by
  `cartoon_enclosures`.
- **Moved**: `test_name_mathtext` (11 -> 8, one canonical case per notation
  form: sub+superscript, inverse symbol, plain word, tagged label, branch code,
  fixed-point letter, already-mathtext, empty) and
  `test_name_mathtext_agrees_with_element_name_mathtext` ->
  `tests/test_element_naming.py`; the `test_verbose_*` tests (5) ->
  `tests/unit/numerics/test_logging_policy.py` (level, logger, non-empty only).

### Deletion ledger, Phase 8

78 collected node ids removed, 80 added: 1111 -> 1113 (`nodeids_p8.txt`,
`p8_deleted.txt`, `p8_added.txt`). Every deleted subject was grepped in the
regression-guard notes (codebase-audit-2026-09, pseudoneighbor-collision-fragility,
pseudoneighbor-partition, numerics-test-suite, regions-refactor-2026-09,
cdist-strict-monotonicity-fix, hole-side-own-blast-2026-10-02,
audit-polish-2026-07): the only plotting hits are design notes (hole marker
per orbit, plot placement), no fixed bug; the face-point copy regression is
kept.

| Node id | New home / reason | Guard checked |
|---|---|---|
| `test_topology_plotting::test_trellis_plotters_delegate_to_plotting` | monkeypatch delegate; behaviour in `P/test_plot_smoke::test_trellis_plotters_draw_on_a_real_tangle` | none |
| `test_topology_plotting::test_plotters_draw_on_a_real_tangle` | `P/test_plot_smoke::test_trellis_plotters_draw_on_a_real_tangle` | none |
| `test_topology_plotting::test_verbose_logs_and_never_prints[*]` ×3, `::test_verbose_prints_when_logging_is_unconfigured`, `::test_verbose_survives_an_application_level_above_info` | `unit/numerics/test_logging_policy.py` (same names; the level-restored check now runs before the test's own restore, so it is no longer tautological) | pseudoneighbor-partition (verbose flags: design only) |
| `test_topology_plotting::test_session_plot_fanouts_go_through_fanout_plot[*]` ×5 | monkeypatch of `_fanout_plot`; behaviour: `P/test_plot_smoke::test_session_plot_fanouts_draw_every_fixed_point` | none |
| `test_topology_plotting::test_session_call_fanouts_go_through_fanout_call[*]` ×4 | monkeypatch of `_fanout_call`; shapes already in `facade/test_session_fanouts.py` | none |
| `test_topology_plotting::test_plot_dual_graph_scatters_exactly_the_stable_nodes` | collection-order + facecolor pin; property `P/test_plot_topology::test_plot_dual_graph_draws_every_stable_node_at_its_point` | none |
| `test_topology_plotting::test_side_nodes_sit_on_their_side_and_solid_nodes_on_the_edge` | `P/test_plot_topology::test_side_nodes_sit_on_their_side_and_unified_nodes_on_the_edge` (no private `_anchorward_look`) | none |
| `test_topology_plotting::test_face_point_of_unbounded_node_is_outside_the_edge_bbox`, `::test_face_point_of_a_region_is_inside_it`, `::test_face_point_of_a_region_returns_a_copy_not_the_cached_array` | `P/test_plot_topology` (same names) | face-point copy: kept as regression |
| `test_topology_plotting::test_scatter_kwargs_reject_facecolors_and_c` | `P/test_plot_smoke::test_plotters_reject_bad_arguments[dual_graph_facecolors|dual_graph_c]` + `test_plot_dual_graph_draws[unclipped_restyled]` | none |
| `test_topology_plotting::test_clip_to_arcs_true_keeps_the_axes_within_the_padded_edge_bbox`, `::test_clip_to_arcs_false_leaves_the_axes_to_autoscale` | deleted: reimplemented the padding (planner §A Phase 8); both options smoke-run in `test_plot_dual_graph_draws` | none |
| `test_topology_plotting::test_show_labels_annotates_every_stable_and_face_node`, `::test_show_labels_ref_style_keeps_the_element_labels`, `::test_stable_node_label_falls_back_to_the_element_label_without_a_name` | deleted: label-separator / wording pins; both label styles smoke-run with the canvas drawn (`test_plot_dual_graph_draws[names|refs]`), `label_style="bogus"` in the ValueError table | none |
| `test_topology_plotting::test_name_mathtext[*]` ×11, `::test_name_mathtext_agrees_with_element_name_mathtext` | `test_element_naming.py::test_name_mathtext[*]` ×8 (one per form; `R_3`, `a_2^-1`, `u` were repeats of a form), `::test_name_mathtext_agrees_with_element_name_mathtext` | none |
| `test_topology_plotting::test_dual_graph_legend_handles_match_the_plotter`, `::test_walk_legend_handles_match_the_plotter`, `::test_cartoon_legend_handles_match_the_plotter` | deleted: legend wording / colour pins; `P/test_plot_smoke::test_legend_handles_are_artists` | none |
| `test_topology_plotting::test_walk_zorder_sits_between_dual_graph_edges_and_nodes` | deleted: z-order pin | none |
| `test_topology_plotting::test_plot_minimal_trellis_draws_every_kept_bridge` | `P/test_plot_topology` (same name; legend labels dropped) | none |
| `test_topology_plotting::test_plot_stable_partition_accepts_row_labels`, `::test_plot_stable_partition_draws_element_labels_at_midpoints` | `P/test_plot_smoke::test_plot_stable_partition_rows_and_element_labels` (one label per interval; the tick-text and midpoint-placement pins dropped); mismatched labels -> ValueError table | none |
| `test_topology_plotting::test_plot_walk_threads_face_node_face`, `::test_plot_walk_honours_overrides_and_trivial_walks` | `P/test_plot_topology::test_plot_walk_threads_face_node_face` (a real walk), `::test_plot_walk_of_a_trivial_or_empty_walk`; overrides smoke `P/test_plot_smoke::test_plot_walk_accepts_style_overrides` (colour/zorder pins dropped) | none |
| `test_topology_plotting::test_class_colors_are_assigned_in_fixed_order` | palette order deleted; no-repeat half -> `P/test_plot_topology::test_class_colors_never_repeat[*]` | none |
| `test_topology_plotting::test_plot_bridges_by_class_draws_every_classed_bridge` | `P/test_plot_topology` (same name; linestyle/legend-text pins dropped; reserved kwargs -> ValueError table) | none |
| `test_topology_plotting::test_unmatched_member_gets_its_parent_letter_slot` | palette-slot pin deleted; property `P/test_plot_topology::test_an_unmatched_member_is_still_drawn_under_its_own_colour` | none |
| `test_topology_plotting::test_plot_transition_graph_draws_every_symbol[*]` | `P/test_plot_topology::test_transition_graph_nodes_are_the_active_symbols[*]` (title/axis/linestyle/label-text pins dropped) | none |
| `test_topology_plotting::test_plot_itinerary_table_lists_every_class`, `::test_plot_itinerary_table_columns_font_and_rows` | width/font/row-height/cell-text pins deleted; `P/test_plot_smoke::test_plot_itinerary_table_draws_one_row_per_class[*]`, `::test_itinerary_table_rows_follow_the_requested_columns`; bad columns -> ValueError table | none |
| `test_topology_plotting::test_cartoon_layout_has_one_node_per_element_per_side`, `::test_cartoon_ordinal_coordinate_is_strictly_increasing`, `::test_cartoon_brackets_match_closedness`, `::test_cartoon_draws_every_kept_bridge_and_walks_element_to_element` | `P/test_plot_topology` (same names; collection-count, colour and gid-format pins dropped; `anchor="top"` -> ValueError table) | none |
| `test_topology_plotting::test_cartoon_toggles_and_session_delegate` | `P/test_plot_smoke::test_session_cartoon_delegate_honours_its_toggles` | none |
| `test_topology_plotting::test_cartoon_bridge_colour_labels_and_walkless_legend` | colour/alpha/label-alignment/font/legend pins deleted; options smoke-run in `P/test_plot_smoke::test_plot_dual_graph_cartoon_draws[*]`; `label_position="above"` -> ValueError table | none |
| `test_higher_period_cartoon::test_circle_layout_puts_every_node_on_its_sides_circle`, `::test_circle_layout_honours_interior_side_and_rejects_bad_shapes` | angle/sweep/radius/legend pins deleted; `P/test_plot_topology::test_circle_interior_side_lies_inside_the_circle`; bad shape / label position -> ValueError table | none |
| `test_higher_period_cartoon::test_line_cartoon_labels_every_row_with_its_branch_code` | deleted: label wording | none |
| `test_higher_period_cartoon::test_p3_circle_follows_the_zone_boundary` | `P/test_plot_topology` (same name; the `cartoon_zones` facts kept, the sweep/gap/anchor-miss angle pins replaced by the cyclic order of the arc anchors) | none |
| `test_higher_period_cartoon::test_nested_inner_circle_sits_inside_the_outer_one`, `::test_nested_outer_arcs_go_around_the_inner_circle` | `P/test_plot_topology::test_nested_inner_circle_sits_inside_the_outer_one` (containment; the enclosure + keep-out facts merged in, gid parsing and polar-route sweep deleted: `ZoneLayout` exposes no routed paths; label texts dropped) | hole-side-own-blast (nested blasts): unrelated to drawing; the blast guards stay in `test_higher_period_cartoon` |
| `test_higher_period_cartoon::test_bridges_by_class_can_be_restricted_to_the_minimal_trellis` | `P/test_plot_topology` (same name) | none |
| `test_higher_period_cartoon::test_nested_classes_are_lettered_and_coloured_tangle_by_tangle` | grouping/cdist order = law `check_class_table_order`; active letters in table order = `test_symbolic_dynamics` (letters a, b, … test); no-repeat -> `test_class_colors_never_repeat[*-nested]`; `TANGLE_COLOR_FAMILIES` membership deleted (colour pin) | none |
| `test_session_dual_graph::test_session_plot_delegates_to_plotting`, `test_session_symbolic_dynamics::test_session_plot_delegates_forward_to_plotting` | monkeypatch delegates; `P/test_plot_smoke::test_session_plot_delegates_return_their_types` (forwarding checked by effect) | none |
| `test_session_dual_graph::test_session_plot_dual_graph_returns_axes`, `test_session_symbolic_dynamics::test_session_plot_delegates_draw_on_the_given_axes` | `P/test_plot_smoke::test_session_plot_delegates_return_their_types` | none |
| `test_session_pseudoneighbors::test_plot_helpers_draw_pairs_and_holes` | `P/test_plot_smoke::test_session_plot_pseudoneighbors_computes_on_demand` (file deleted, empty) | none |
| `test_session_pseudoneighbors::test_plot_stable_partition_fans_out`, `test_session_strong_pips::test_plot_helpers_cover_every_tangle` | `P/test_plot_smoke::test_session_plot_fanouts_draw_every_fixed_point` (both files deleted, empty) | none |
| `test_stable_partition::test_plot_stable_partition_smoke` | tick-count pin dropped; smoke on real data in `P/test_plot_smoke::test_plot_stable_partition_rows_and_element_labels`; its non-plot half (the trellis stores a one-branch partition) -> `test_stable_partition::test_trellis_partition_of_one_branch_is_stored` | none |
| `test_stable_partition::test_henon_plot_helpers_smoke` | same recipe (`build_k10(through="bridges")`) as `P/test_plot_smoke::test_trellis_plotters_draw_on_a_real_tangle` | none |

### Verification

- Collected **1113** (`nodeids_p8.txt`).
- `1073 passed, 32 skipped, 8 xfailed` in 166 s with coverage (`p8_run.txt`);
  `-rxX` lists exactly the 8 `KNOWN_ISSUES` xfails; no XPASS.
- Coverage guard vs `cov_base.json`: OK (`cov_p8.json`); `plotting.py`
  91.39% -> 91.86%.
- Isolation (each alone): `test_dual_graph.py::test_k10_a_different_pip_moves_the_unified_set`,
  `plotting/test_plot_topology.py::test_nested_inner_circle_sits_inside_the_outer_one`,
  `plotting/test_plot_smoke.py::test_plotters_reject_bad_arguments[circle_label_position]`,
  `unit/numerics/test_logging_policy.py::test_verbose_prints_when_logging_is_unconfigured`,
  `test_element_naming.py::test_name_mathtext[fixed_point_letter]`,
  `test_stable_partition.py::test_trellis_partition_of_one_branch_is_stored`,
  `test_image_cdist.py::test_the_snap_window_is_tighter_than_the_agreement_bound`,
  `invariants/test_law_iterated.py::test_iterated_unique_owner[p3]`: all pass.
- `ruff check --select F` on every touched test file: clean.

### Deviations, Phase 8

1. **Case count above the planner's estimate**: the plotting tier collects 65
   cases (41 smoke, 24 property) against ≈ 28, because option sweeps and the
   `ValueError` table are parametrized (one id per option, each build ≈ 0.05 s)
   rather than folded into loops. Same assertions either way; the run is still
   green in ≈ 8 s for the whole directory.
2. **New files created in their final directories** (`tests/plotting/`,
   `tests/unit/numerics/test_logging_policy.py`), as Phases 4–6 did; Phase 10
   moves only the survivors of the old flat layout. The logging-policy file
   sits under `unit/numerics/` as planned although the verbose reports it
   checks are the topology `Trellis`'s (the policy is cross-cutting).
3. **`test_nested_outer_arcs_go_around_the_inner_circle` not rewritten as a
   keep-out property**: `ZoneLayout` exposes no routed paths, so per the
   planner it is deleted; its public half (`cartoon_enclosures` reports the
   enclosing lobes, and `keepout` encloses the inner circle) is merged into
   the containment test. `_polar_route` stays covered by the nested circle
   draw.
4. **The p3 ring-order fact** ("0.0 -> 1.0 -> 2.0 clockwise, interior on the
   right", CLAUDE.md) stays in the plotting tier as a cyclic-order property;
   it is a fixture fact not in `tests/golden/`. Author item: move it to
   `golden/test_golden_period3.py` (needs sign-off) or leave it here.
   CLOSED (author follow-up, golden removed): the ring-order / clockwise /
   right-interior pins were dropped; the test now checks the circle against
   whatever order and direction `cartoon_zones()` reports.
5. **Non-plot half of `test_plot_stable_partition_smoke` kept** as
   `test_stable_partition::test_trellis_partition_of_one_branch_is_stored`
   (the one-branch `Trellis.partition_stable_manifold` store had no other
   test).
6. `test_session_strong_pips.py` and `test_session_pseudoneighbors.py` are
   deleted outright: every test in them was a plot helper.

---

## §9 Phase 9: production-check wiring, gap tests, open deep-p3 runs, GPU, perf

### What landed

- **`tests/unit/topology/test_production_checks.py`** (decision 12): spies on
  `StablePartition.check_bridge_rows_consistent` / `check_holes_share_bridge_side`
  (module-attribute patches; both are imported inside `Trellis.punch_holes`).
  One `punch_holes` call on p3 runs I2 exactly once per hole-bearing bridge
  (9 bridges for 12 holes) on that trellis and I1 once on all holes with
  `orientation_preserving=trellis.orientation_preserving`. Negative test: a
  propagated hole's side flipped inside `propagate_reference_holes` makes
  `punch_holes` raise `AssertionError` (type only). `DualGraph` construction
  runs its self-check (`_check`) exactly once. The Arrangement's
  consecutive-bridge guard keeps its existing test
  (`test_arrangement::test_every_bridge_must_be_consecutive_on_its_branch`).
- **`tests/regression/test_open_p3_deep_runs.py`** (decision 10, P3), runs
  `p3_15`, `p3_16`, `p3_6_blasts` (fresh builds, ~0.2 s each):
  `test_deep_p3_holes_share_bridge_side` (I1 holds, holes non-empty, the
  symbolic dynamics runs; passes), `test_deep_p3_every_class_resolves_without_virtual_symbols`
  (`xfail(strict=True)` through `issue_marks("open_p3_deep_runs", ...)`, the
  `KNOWN_ISSUES` entries Phase 1 registered), and
  `test_virtual_symbols_are_transition_sinks` (real data on `p3_15`: a virtual
  node has out-degree 0, in-degree > 0 and a zero matrix row; skips if the
  run ever stops producing one).
- **Gap tests:**
  - `unit/numerics/test_same_stability_discard.py` (u×u and s×s): two
    near-tangent same-stability polylines straddling each other plus a
    transversal curve of the other stability; `Tangle.resolve_crossings`
    returns only the two unstable × stable crossings and logs DEBUG on
    `tanglepack.numerics.Tangle` (level and logger only).
  - `test_bridge_class::test_a_cross_branch_class_runs_from_the_anchor_nearer_element`
    (hand-built `ElementRef`s on two orbit branches: smaller element id is the
    source whichever branch, a tie falls back to the run-stable order) and
    `::test_classes_group_by_tangle_with_connecting_classes_last` (a duck
    trellis with two fixed points: tangle 0, tangle 1, then the connecting
    class with `tangle is None`, even though it has the smallest cdist). The
    nested law half already runs (`class_table_order`,
    `classes_do_not_mix_tangles` on nested); nested has no connecting class,
    hence the hand-built case.
  - `unit/topology/test_iterated_cut.py`: the two synthetic 09-30 rules
    (chords only on the base branch; a chain ending on another branch marks
    nothing) MOVED from `test_higher_period_cartoon.py`, plus the new
    `test_a_far_end_on_another_branch_is_no_flank` (nested, public
    `from_minimal`; one foreign far end injected per image crossing: more
    WARNINGs than the unpatched run, identical partition signature, every
    applied cut's partner on its own branch). Covers `PartitionFamily` 746.
  - `test_arrangement::test_hand_built_sub_face_is_rejected_by_area` (gap 12):
    a deterministic lobe-in-a-lobe (`_nested_lobe_trellis`): the one face
    carrying the mapped cycle falls short of the source area by more than
    `REGION_AREA_RTOL`, the sub-face plus the little face recover it, and
    `image_of` returns None. The k=10 fixture test stays.
- **GPU** (`unit/numerics/test_gpu.py`, replaces `numerics/test_gpu.py`):
  parity runs by default (`importorskip("cupy")`, skips with no CUDA device;
  here CuPy 14.1.1 with one device, so it RUNS), now through
  `enable_gpu(session)` and checks `disable_gpu` restores the CPU map;
  `test_enable_gpu_without_cupy_raises_clearly` hides CuPy with
  `monkeypatch.setitem(sys.modules, "cupy", None)` and RUNS (ImportError, map
  unchanged), retiring the suite's one non-NOT_APPLICABLE skip; new
  `test_enable_gpu_rejects_a_target_without_a_system` (TypeError). `gpu.py`
  70.9 % -> 92.7 %.
- **Perf:** `test_registry_insert_is_near_linear` moved to
  `unit/numerics/test_registry_perf.py` under `@pytest.mark.perf`;
  `pytest -m perf` selects exactly it, the default run deselects it.
- **Regression tier (left over from Phase 7a):**
  `regression/test_cdist_collision_growth.py` and
  `test_high_stretch_period3_growth.py` merged into
  `regression/test_high_stretch_growth.py` (two functions, same assertions;
  the vacuity of the period-3 run is recorded in its Dev Notes);
  `test_blast_monotonicity` strengthened: after the blast every manifold and
  every bridge has non-decreasing cdist and no geometric spike.

### Deletion ledger, Phase 9

9 node ids removed, 24 added (1113 -> 1128 collected, 1 of them perf,
deselected by default). Lists: `p9_deleted.txt`, `p9_added.txt`
(`nodeids_p9.txt`).

| Node id(s) | New home | Guard checked |
|---|---|---|
| `test_stable_partition_period3.py::test_deep_p3_holes_share_bridge_side_through_symbolic_dynamics[15_steps]`, `[6_blasts]` (file now empty, deleted) | `regression/test_open_p3_deep_runs.py::test_deep_p3_holes_share_bridge_side[p3_15\|p3_16\|p3_6_blasts]` (16 steps added) | hole-side-own-blast-2026-10-02: I1 on the deep runs kept |
| `test_higher_period_cartoon.py::test_chords_pair_only_the_base_branchs_crossings`, `::test_a_lobe_ending_on_another_branch_marks_nothing` | `unit/topology/test_iterated_cut.py` (same assertions) | higher-period 09-30 rules kept |
| `numerics/test_gpu.py::test_gpu_growth_matches_cpu`, `::test_enable_gpu_without_cupy_raises_clearly` | `unit/numerics/test_gpu.py` (the no-CuPy test now runs) | numerics-speedups: GPU parity kept |
| `numerics/test_generation_and_caches.py::test_registry_insert_is_near_linear` | `unit/numerics/test_registry_perf.py` (`@pytest.mark.perf`, body unchanged incl. the `gc.collect()` flake fix) | none |
| `regression/test_cdist_collision_growth.py::test_growth_keeps_geometry_smooth`, `regression/test_high_stretch_period3_growth.py::test_period3_high_stretch_growth_is_not_scrambled` | `regression/test_high_stretch_growth.py::test_k10_growth_keeps_geometry_smooth`, `::test_period3_high_stretch_growth_is_not_scrambled` | cdist-strict-monotonicity-fix: same geometric no-spike + non-strict checks |

The planned deletions of the I1/I2 re-check tests
(`test_topology_invariants::test_henon_holes_share_bridge_side`,
`test_stable_partition_period3::test_p3_holes_share_bridge_side`,
`::test_p3_bridge_rows_consistent`) were already done in Phase 4 (deviation
4.3); nothing left to delete there.

### Verification

- Collected **1128** (1127 by default + 1 perf).
- `1085 passed, 31 skipped, 1 deselected, 11 xfailed` in 172 s with coverage
  (`p9_run.txt`). `-rxX` lists exactly the `KNOWN_ISSUES`: the 8 law-tier
  xfails plus the 3 `open_p3_deep_runs`; no XPASS. Skips: the 31
  `NOT_APPLICABLE` cases only.
- Coverage guard vs `cov_base.json`: OK (`cov_p9.json`).
- `pytest -m perf`: 1 selected, passes (0.4 s).
- Isolation (each alone): `test_dual_graph.py::test_k10_a_different_pip_moves_the_unified_set`,
  `unit/topology/test_iterated_cut.py::test_a_far_end_on_another_branch_is_no_flank`,
  `regression/test_open_p3_deep_runs.py::test_virtual_symbols_are_transition_sinks`,
  `unit/numerics/test_gpu.py::test_enable_gpu_without_cupy_raises_clearly`,
  `test_bridge_class.py::test_classes_group_by_tangle_with_connecting_classes_last`,
  `regression/test_blast_monotonicity.py::test_blast_completes_without_monotonicity_failure`,
  `test_fixed_point.py::test_advance_key_rejects_a_foreign_fixed_point`: all pass.
- `ruff check --select F` on every new/touched file: clean (one pre-existing
  F841 in `test_bridge_class.py::test_table_lookups_symbols_and_report`, not
  touched).

### Deviations, Phase 9

1. **Gaps already covered, no new test:** `preserve_ids` / idempotent
   recompute and area preservation (laws since Phase 4); `min_separation`
   (`numerics/test_blast_proximity_guard::test_min_separation_drops_close_bridges`);
   refined cdist strictly between its neighbours (`numerics/test_refinement.py`);
   ambiguous walks / unmatched members / evidence (synthetic tests in
   `test_dual_walk.py` and `test_symbolic_dynamics.py`; P3 found none on real
   data).
2. **Tangle grouping with connecting classes is hand-built only:** no law case
   (nested included) has a class with `tangle is None`, so the nested law
   cannot exercise the "connecting last" half; the duck trellis
   (`_ClassTrellis`, crossings only, no iterates) is local to the test.
3. **Cross-branch flank on real data, not synthetic:** the flank rule lives in
   `from_minimal`, so the test runs the public path on nested with a wrapped
   `_empty_stretches` (allowed kernel) injecting foreign far ends.
4. **The DualGraph self-check spy wraps the private `DualGraph._check`**: it
   is the production-wiring assertion the plan asks for (that construction
   runs it), not a test of the method's internals.
5. **New files created in their final directories** (`unit/`, `regression/`),
   as Phases 4-8 did; the regression merge and blast strengthening (Phase 7a
   deviation 4) are done here. The period-3 high-stretch run is still vacuous
   (AUTHOR ITEM from Phase 7a: re-derive its parameters); recorded in the
   merged file's Dev Notes instead of a failing non-vacuity check.
6. No new `slow` / `regression` markers on new tests (Phase 10 drops both).

## §10 Phase 10: move once, finalize

Two commits: `117dfb3` (move only) and the finalize commit that carries this
section.

### Commit 1: move only (`git mv`, 48 renames, zero content change)

Every survivor outside the tiered layout moved into it; git records each as a
100 % rename. Basenames stay unique (pytest's default import mode needs that,
since no directory gets an `__init__.py`). Four files were renamed on the way
to their planned names: `test_point` → `test_linked_list`,
`test_strong_pip_periodic_point` → `test_strong_pip`, `test_loom_blast_restore`
→ `test_blast_restore`, `test_resonance_zone_region` → `test_resonance_zone`;
`test_higher_period_cartoon` → `regression/test_nested_own_blast`. Map
(`p10a_moves.tsv`):

| Old path | New path |
|---|---|
| `tests/test_report_smoke.py` | `tests/facade/test_report_smoke.py` |
| `tests/test_session_bridge_classes.py` | `tests/facade/test_session_bridge_classes.py` |
| `tests/test_session_dual_graph.py` | `tests/facade/test_session_dual_graph.py` |
| `tests/test_session_symbolic_dynamics.py` | `tests/facade/test_session_symbolic_dynamics.py` |
| `tests/test_higher_period_cartoon.py` | `tests/regression/test_nested_own_blast.py` |
| `tests/numerics/test_blast_no_overlap.py` | `tests/unit/loom/test_blast_no_overlap.py` |
| `tests/numerics/test_blast_proximity_guard.py` | `tests/unit/loom/test_blast_proximity_guard.py` |
| `tests/test_loom_blast_restore.py` | `tests/unit/loom/test_blast_restore.py` |
| `tests/test_resonance_zone_region.py` | `tests/unit/loom/test_resonance_zone.py` |
| `tests/numerics/test_batched_map.py` | `tests/unit/numerics/test_batched_map.py` |
| `tests/numerics/test_bridge_identity.py` | `tests/unit/numerics/test_bridge_identity.py` |
| `tests/numerics/test_cleanup_walkers_and_examples.py` | `tests/unit/numerics/test_cleanup_walkers_and_examples.py` |
| `tests/numerics/test_closed_form_curvature.py` | `tests/unit/numerics/test_closed_form_curvature.py` |
| `tests/test_fixed_point.py` | `tests/unit/numerics/test_fixed_point.py` |
| `tests/numerics/test_generation_and_caches.py` | `tests/unit/numerics/test_generation_and_caches.py` |
| `tests/numerics/test_geometry.py` | `tests/unit/numerics/test_geometry.py` |
| `tests/numerics/test_grow_until.py` | `tests/unit/numerics/test_grow_until.py` |
| `tests/numerics/test_growth_integration.py` | `tests/unit/numerics/test_growth_integration.py` |
| `tests/numerics/test_initializer_cdist.py` | `tests/unit/numerics/test_initializer_cdist.py` |
| `tests/numerics/test_intersection_registry_fixes.py` | `tests/unit/numerics/test_intersection_registry_fixes.py` |
| `tests/numerics/test_invariant_helpers.py` | `tests/unit/numerics/test_invariant_helpers.py` |
| `tests/numerics/test_inversion_fixture.py` | `tests/unit/numerics/test_inversion_fixture.py` |
| `tests/test_point.py` | `tests/unit/numerics/test_linked_list.py` |
| `tests/test_manifold_initializer.py` | `tests/unit/numerics/test_manifold_initializer.py` |
| `tests/numerics/test_map_step_and_graph.py` | `tests/unit/numerics/test_map_step_and_graph.py` |
| `tests/numerics/test_minimal_solver.py` | `tests/unit/numerics/test_minimal_solver.py` |
| `tests/numerics/test_noise_crossing_collapse.py` | `tests/unit/numerics/test_noise_crossing_collapse.py` |
| `tests/numerics/test_refinement.py` | `tests/unit/numerics/test_refinement.py` |
| `tests/numerics/test_segment_index_bookkeeping.py` | `tests/unit/numerics/test_segment_index_bookkeeping.py` |
| `tests/numerics/test_single_source_of_truth.py` | `tests/unit/numerics/test_single_source_of_truth.py` |
| `tests/numerics/test_tangle_index_and_orientation.py` | `tests/unit/numerics/test_tangle_index_and_orientation.py` |
| `tests/numerics/test_workbench_bugfixes.py` | `tests/unit/numerics/test_workbench_bugfixes.py` |
| `tests/test_arrangement.py` | `tests/unit/topology/test_arrangement.py` |
| `tests/test_arrangement_sparse.py` | `tests/unit/topology/test_arrangement_sparse.py` |
| `tests/test_bridge_class.py` | `tests/unit/topology/test_bridge_class.py` |
| `tests/test_dual_graph.py` | `tests/unit/topology/test_dual_graph.py` |
| `tests/test_dual_walk.py` | `tests/unit/topology/test_dual_walk.py` |
| `tests/test_element_naming.py` | `tests/unit/topology/test_element_naming.py` |
| `tests/test_image_cdist.py` | `tests/unit/topology/test_image_cdist.py` |
| `tests/test_iterated_partition.py` | `tests/unit/topology/test_iterated_partition.py` |
| `tests/test_minimal_trellis.py` | `tests/unit/topology/test_minimal_trellis.py` |
| `tests/test_partition_elements.py` | `tests/unit/topology/test_partition_elements.py` |
| `tests/test_partition_family.py` | `tests/unit/topology/test_partition_family.py` |
| `tests/test_pseudoneighbor.py` | `tests/unit/topology/test_pseudoneighbor.py` |
| `tests/test_stable_partition.py` | `tests/unit/topology/test_stable_partition.py` |
| `tests/test_strong_pip_periodic_point.py` | `tests/unit/topology/test_strong_pip.py` |
| `tests/test_symbolic_dynamics.py` | `tests/unit/topology/test_symbolic_dynamics.py` |
| `tests/test_topology_invariants.py` | `tests/unit/topology/test_topology_invariants.py` |

1128 collected before and after; `nodeids_p10a.txt` is `nodeids_p9.txt` with
the paths renamed (497 ids moved, none added or removed).

### Commit 2: finalize

- **Markers:** `@pytest.mark.slow` (24 uses) and `@pytest.mark.regression`
  (5 uses) removed from 16 files; their registrations removed from
  `pyproject.toml`. Only `golden` and `perf` remain (`--strict-markers`).
  The `import pytest` lines left unused were removed (ruff F401, which also
  dropped three pre-existing unused imports: `numpy` in
  `test_arrangement_sparse`, `Trellis` in `test_partition_elements`, `pytest`
  in `test_workbench_bugfixes`).
- **Shims deleted:** `tests/invariants.py` and `tests/walk_helpers.py`; the
  three remaining `from invariants import` lines now import
  `helpers.invariants`. `tests/minimal_helpers.py` was already gone (Phase 7b,
  deviation 1).
- **Split finished:** the two naming tests of
  `regression/test_nested_own_blast.py` moved to
  `unit/topology/test_element_naming.py` (planner §B: names →
  `test_element_naming`, words/blast order → `R/test_nested_own_blast`); same
  assertions, built inline from `build_period3()` / `build_nested()`. The
  regression module docstring now states the 2026-10-02 bug it guards.
- **Docstrings:** `tests/conftest.py` gained the tier layout, the markers and
  the **Dev Notes** (author decision 14 / E3): the untested provisional rules
  (the `+1..+(k-1)` exemption, the mean-of-neighbours refined cdist, the
  anchor faces-closed limitation, each with where its firm part is tested),
  the dead API kept without tests (whitelisted in `coverage_guard.py`), and
  the `KNOWN_ISSUES` index. The stale `p3_partitioned` fixture docstring
  (`henon_cases.build_nested` → `cases.build_nested`), `helpers/fakes.py`
  ("shim stays until Phase 10"), and five stale test-path references
  (`tests/invariants.py`, `test_high_stretch_period3_growth.py`,
  `tests/test_image_cdist.py`, a deleted test named in `test_image_cdist`)
  fixed.
- **CLAUDE.md** "Run tests": `env/bin/python -m pytest` commands for the
  default run, `-rxX`, one tier, one test, `-m golden`, `-m perf`, the
  coverage guard, and one paragraph on the tier layout. Nothing else in
  CLAUDE.md changed.
- **`tests/_tools/nodeid_ledger.py`** (new, not collected): walks
  `nodeids_base.txt → p2 → … → p9 → p10a → p10` in the run folder and
  requires every removed id to be in that phase's `<phase>_deleted.txt` or
  carried by its `<phase>_moves.tsv`, and every listed deletion to name its
  test function verbatim in this ledger. The appendix below lists verbatim
  the 75 ids earlier phases recorded only in grouped form.

### Deletion ledger, Phase 10

2 node ids removed, 2 added (1128 → 1128; `p10_deleted.txt`,
`p10_added.txt`, `nodeids_p10.txt`). Commit 1 removed none (moves only).

| Node id | New home | Guard checked |
|---|---|---|
| `tests/regression/test_nested_own_blast.py::test_p3_names_carry_orbit_codes_without_a_letter` | `tests/unit/topology/test_element_naming.py::test_p3_names_carry_orbit_codes_without_a_letter` (same assertions) | higher-period-figures-2026-09-30: branch codes on every name, kept |
| `tests/regression/test_nested_own_blast.py::test_nested_names_carry_fixed_point_letters` | `tests/unit/topology/test_element_naming.py::test_nested_names_carry_fixed_point_letters` (same assertions) | higher-period-figures-2026-09-30: fixed-point letters, kept |

### Node-id accounting (end-to-end check)

`env/bin/python tests/_tools/nodeid_ledger.py .refactor/runs/2026-10-05-test-suite docs/test_suite_refactor_ledger.md`:

| Transition | Removed | Deleted (listed) | Moved | Added | Collected |
|---|---|---|---|---|---|
| base → p2 | 74 | 74 | 0 | 2 | 724 |
| p2 → p3a | 0 | 0 | 0 | 2 | 726 |
| p3a → p3b | 0 | 0 | 0 | 0 | 726 |
| p3b → p4 | 95 | 95 | 0 | 449 | 1080 |
| p4 → p5 | 19 | 19 | 0 | 21 | 1082 |
| p5 → p6 | 40 | 40 | 0 | 70 | 1112 |
| p6 → p7a | 54 | 54 | 0 | 50 | 1108 |
| p7a → p7b | 24 | 24 | 0 | 27 | 1111 |
| p7b → p8 | 78 | 78 | 0 | 80 | 1113 |
| p8 → p9 | 9 | 9 | 0 | 24 | 1128 |
| p9 → p10a | 497 | 0 | 497 | 497 | 1128 |
| p10a → p10 | 2 | 2 | 0 | 2 | 1128 |

Of the 796 baseline ids, 3 survive at their original path and 793 are gone:
each is either deleted (listed in a phase's deletion list and in this ledger)
or carried, renamed, by the Phase 10 move. 0 unaccounted, 0 unexplained.
(Phase 1 kept every id and has no snapshot; Phases 3a/3b removed none.)

#### Appendix: ids earlier phases recorded in grouped form

Each id below is covered by a grouped row of its phase's deletion table
(whole-file rows, `×N`, `*`, `k10_/p3_` twins, `{a,b}` sets); listed here
verbatim so the accounting is mechanical.

**§2** (26 ids; `p2_deleted.txt`):

- `tests/numerics/test_intersection_registry_fixes.py::test_get_lambda_u_falls_back_to_the_other_key`
- `tests/numerics/test_intersection_registry_fixes.py::test_get_lambda_u_honours_the_stability_argument`
- `tests/numerics/test_intersection_registry_fixes.py::test_get_lambda_u_returns_a_python_float`
- `tests/numerics/test_intersection_registry_fixes.py::test_get_lambda_u_returns_none_without_keys`
- `tests/numerics/test_workbench_split.py::test_collaborators_are_wired_to_their_workbench`
- `tests/numerics/test_workbench_split.py::test_delegating_signatures_unchanged[build_intersection_graph-params6]`
- `tests/numerics/test_workbench_split.py::test_delegating_signatures_unchanged[image_bridges-params2]`
- `tests/numerics/test_workbench_split.py::test_delegating_signatures_unchanged[infer_iterate_table-params5]`
- `tests/numerics/test_workbench_split.py::test_delegating_signatures_unchanged[infer_iterates-params4]`
- `tests/numerics/test_workbench_split.py::test_delegating_signatures_unchanged[iterate_all_bridges-params1]`
- `tests/numerics/test_workbench_split.py::test_delegating_signatures_unchanged[iterate_bridge-params0]`
- `tests/numerics/test_workbench_split.py::test_delegating_signatures_unchanged[preimage_bridges-params3]`
- `tests/numerics/test_workbench_split.py::test_iterate_bridge_still_patchable_on_the_workbench`
- `tests/numerics/test_workbench_split.py::test_moved_bodies_live_in_their_new_module[build_intersection_graph-graphviz]`
- `tests/numerics/test_workbench_split.py::test_moved_bodies_live_in_their_new_module[image_bridges-BridgeIterator]`
- `tests/numerics/test_workbench_split.py::test_moved_bodies_live_in_their_new_module[infer_iterate_table-IterateInference]`
- `tests/numerics/test_workbench_split.py::test_moved_bodies_live_in_their_new_module[infer_iterates-IterateInference]`
- `tests/numerics/test_workbench_split.py::test_moved_bodies_live_in_their_new_module[iterate_all_bridges-BridgeIterator]`
- `tests/numerics/test_workbench_split.py::test_moved_bodies_live_in_their_new_module[iterate_bridge-BridgeIterator]`
- `tests/numerics/test_workbench_split.py::test_moved_bodies_live_in_their_new_module[preimage_bridges-BridgeIterator]`
- `tests/numerics/test_workbench_split.py::test_moved_bodies_live_in_their_new_module[visualize_intersection_graph-graphviz]`
- `tests/numerics/test_workbench_split.py::test_top_level_import_still_works`
- `tests/numerics/test_workbench_split.py::test_visualize_signature_unchanged`
- `tests/test_manifold_machine.py::test_machine_initialization`
- `tests/test_manifold_machine.py::test_new_grow_manifold_matches_old_period_one`
- `tests/test_manifold_machine.py::test_new_grow_manifold_matches_old_period_three`

**§4** (33 ids; `p4_deleted.txt`):

- `tests/numerics/test_machine_iterate_invariants.py::test_stable_growth_preserves_invariants[1]`
- `tests/numerics/test_machine_iterate_invariants.py::test_stable_growth_preserves_invariants[2]`
- `tests/numerics/test_machine_iterate_invariants.py::test_stable_growth_preserves_invariants[4]`
- `tests/numerics/test_machine_iterate_invariants.py::test_unstable_growth_preserves_invariants[1]`
- `tests/numerics/test_machine_iterate_invariants.py::test_unstable_growth_preserves_invariants[2]`
- `tests/numerics/test_machine_iterate_invariants.py::test_unstable_growth_preserves_invariants[4]`
- `tests/numerics/test_machine_iterate_invariants.py::test_unstable_growth_preserves_invariants[5]`
- `tests/numerics/test_no_same_stability_crossing.py::test_period3_unstable_manifolds_do_not_self_cross`
- `tests/numerics/test_tangle_intersection_cdist.py::test_crossing_cdist_lies_between_segment_endpoints`
- `tests/numerics/test_tangle_intersection_cdist.py::test_each_crossing_is_one_unstable_one_stable`
- `tests/numerics/test_tangle_intersection_cdist.py::test_intersections_have_both_cdists`
- `tests/test_arrangement.py::test_k10_arrangement_is_one_component_and_euler_holds`
- `tests/test_arrangement.py::test_k10_arrangement_regions_are_geometrically_sound`
- `tests/test_arrangement.py::test_k10_region_images_agree_with_the_dynamics`
- `tests/test_arrangement.py::test_k10_regions_are_pairwise_disjoint`
- `tests/test_arrangement.py::test_p3_arrangement_is_two_components_and_euler_holds`
- `tests/test_bridge_class.py::test_row_of_end_agrees_with_the_geometry_on_p3`
- `tests/test_bridge_class.py::test_same_branch_bridges_never_mismatch_on_p3`
- `tests/test_dual_graph.py::test_k10_face_side_agrees_with_the_geometry`
- `tests/test_dual_graph.py::test_k10_node_structure`
- `tests/test_dual_graph.py::test_k10_payload`
- `tests/test_dual_graph.py::test_k10_unified_nodes_are_the_edges_between_the_pip_and_its_image`
- `tests/test_dual_graph.py::test_p3_unifies_the_pip_segment_on_each_pips_own_branch`
- `tests/test_iterated_partition.py::test_k10_iterated_partition_invariants`
- `tests/test_iterated_partition.py::test_k28_iterated_partition_invariants`
- `tests/test_iterated_partition.py::test_p3_iterated_partition_invariants`
- `tests/test_minimal_trellis.py::test_k10_minimal_trellis_invariants`
- `tests/test_minimal_trellis.py::test_k28_minimal_trellis_invariants`
- `tests/test_minimal_trellis.py::test_p3_minimal_trellis_invariants`
- `tests/test_partition_elements.py::test_element_ids_are_positional_and_carry_their_branch`
- `tests/test_partition_elements.py::test_element_of_intersection_covers_every_crossing_on_the_branch`
- `tests/test_partition_elements.py::test_simple_tangle_elements_own_every_crossing_once`
- `tests/test_partition_elements.py::test_singleton_elements_own_exactly_their_own_point`

**§6** (3 ids; `p6_deleted.txt`):

- `tests/test_session_trellis_cache.py::test_add_resonance_zones_invalidates_the_cache`
- `tests/test_session_trellis_cache.py::test_iterate_bridge_invalidates_the_cache`
- `tests/test_session_trellis_cache.py::test_rebuild_bridges_invalidates_the_cache`

**§7a** (13 ids; `p7a_deleted.txt`):

- `tests/numerics/test_generation_and_caches.py::test_registry_generation_bumps_on_add`
- `tests/numerics/test_generation_and_caches.py::test_registry_generation_bumps_on_add_synthetic`
- `tests/numerics/test_generation_and_caches.py::test_registry_generation_bumps_on_register_iterate`
- `tests/numerics/test_generation_and_caches.py::test_registry_generation_bumps_on_reindex_from`
- `tests/numerics/test_intersection_registry_fixes.py::test_fixed_points_distinct_fixed_points_in_ab_order`
- `tests/numerics/test_intersection_registry_fixes.py::test_fixed_points_same_fixed_point_is_deduped`
- `tests/numerics/test_intersection_registry_fixes.py::test_fixed_points_with_no_keys_is_empty`
- `tests/numerics/test_intersection_registry_fixes.py::test_fixed_points_with_only_a_key`
- `tests/numerics/test_intersection_registry_fixes.py::test_fixed_points_with_only_b_key`
- `tests/test_branch_point.py::test_insert_point_backward`
- `tests/test_branch_point.py::test_insert_point_backward_connected`
- `tests/test_branch_point.py::test_insert_point_forward`
- `tests/test_branch_point.py::test_insert_point_forward_connected`

### Verification

- Collected **1128** (1127 by default + 1 `perf`), unchanged by Phase 10.
- `env/bin/python -m pytest -q -rxX --durations=15 --cov=tanglepack`:
  `1085 passed, 31 skipped, 1 deselected, 11 xfailed` in 159 s with coverage
  (93 s without; `p10_run.txt`). `-rxX` lists exactly the 11 `KNOWN_ISSUES`
  (8 law tier + 3 `open_p3_deep_runs`); no XPASS. Skips: the 31
  `NOT_APPLICABLE` cases only.
- Coverage guard vs `cov_base.json`: OK (`cov_p10.json`).
- `-m golden` selects exactly the three golden files (7 + 9 + 5 = 21 cases);
  `-m perf` selects exactly `unit/numerics/test_registry_perf.py::test_registry_insert_is_near_linear`
  (passes); the default run deselects it.
- Each tier alone: invariants 407 passed / 31 skipped / 8 xfailed (36.5 s);
  unit/numerics 233; unit/topology 249; unit/loom 20; facade 79; golden 21;
  regression 11 + 3 xfailed; plotting 65 — together the 1127 default cases.
- Isolation (each alone): `unit/topology/test_dual_graph.py::test_k10_a_different_pip_moves_the_unified_set`,
  `unit/topology/test_element_naming.py::test_nested_names_carry_fixed_point_letters`,
  `plotting/test_plot_smoke.py::test_plotters_reject_bad_arguments[bridges_by_class_color]`,
  `facade/test_session_caches.py::test_trellis_misses_after_every_mutation_path[add_resonance_zones]`,
  `unit/numerics/test_generation_and_caches.py::test_workbench_generation_bumps_on_every_mutation[create_bridges]`,
  `unit/topology/test_production_checks.py::test_dual_graph_construction_runs_its_self_check`;
  `invariants/test_law_bridges.py` alone: all pass.
- `ruff check --select F tests`: only the four pre-existing F841 unused
  locals remain (not touched).

### Deviations, Phase 10

1. **No file merges in the move.** Planner §C names consolidated targets
   (`unit/numerics/test_growth.py`, `test_registry.py`, `test_bridges.py`,
   `unit/loom/test_blast.py`, …) that would each merge several survivors.
   Merging is a content change, which a move-only commit cannot carry and
   which would break git's rename detection, so each survivor moved whole
   under its own basename (four renamed to their planned 1:1 names). The
   tier of every file matches §B; the per-test splits across tiers were
   already done in Phases 4–9.
2. **Session report/delegate tests in `facade/`:** `test_report_smoke`,
   `test_session_bridge_classes`, `test_session_dual_graph`,
   `test_session_symbolic_dynamics` (§B sends their surviving halves to
   F/ and U/t; what is left after Phases 5–8 is session-level).
3. **Final size and wall time over the plan's target.** 1128 cases (plan:
   ~550–650) and 93 s without coverage (plan: ≤ ~90 s): the law tier is one
   test per (law × case), 446 cases in ~37 s, the "per-check" option
   planner finding 1 / E1 priced at ~860 total; Phase 4 chose it so each
   `KNOWN_ISSUES` xfail stays precise. Not reduced here (out of a move
   phase's scope); flagged for the author.
4. **`tests/_tools/nodeid_ledger.py` added** as the diff script the plan's
   end-to-end verification asks for; its ledger-name check needed the
   appendix above because earlier phases recorded 75 ids in grouped form
   (all checked by hand against their grouped rows: `_get_lambda_u ×4`,
   `test_workbench_split.py` whole file, `test_manifold_machine.py` whole
   file, `test_machine_iterate_invariants.py` whole file, the k10/p3 twins of
   Phase 4, `test_session_trellis_cache` mutation rows, the registry-bump and
   `test_fixed_points_*` merges, `test_branch_point::*`).
5. **Mutation sanity check (plan "Verification")** not run in this phase: it
   needs a scratch branch with a `src/` edit, which this phase's rules
   forbid on the working branch. Left for the author or a follow-up.

## Author follow-up (2026-10-05): law tier per layer x case, the hole rule

### What landed

- **Law tier collapsed from one test per (law x case) to one test per
  (LAYER x case).** `tests/invariants/test_law_<layer>.py` each hold
  `test_<layer>_laws = layer_test("<layer>")`, parametrized over
  `cases.LAW_CASES`. `helpers/law_tier.run_layer` runs every check of the
  layer, insists each applicable check counted at least one item, and raises
  ONE `AssertionError` listing every failing or vacuous check with its message
  and failing source line. `helpers.laws` / `helpers.invariants` are now
  assert-rewritten (`pytest.register_assert_rewrite` at the top of
  `tests/conftest.py`), so the messages carry the compared values.
- **Known issues stay precise.** A `(law, case)` pair in `KNOWN_ISSUES` is left
  out of its layer run and gets its own `test_<layer>_known_issue[<law>-<case>]`
  test, `xfail(strict=True)` (`known_issue_test`), so a fix flips to XPASS.
  `NOT_APPLICABLE` pairs are skipped inside the layer and not counted; a layer
  with nothing applicable on a case would be a skip (none today).
- **Builds:** module-scoped read-only `law_case` + the fingerprint guard kept.
  The two mutating laws became the `recompute` layer
  (`laws.MUTATING_LAYERS`), run in order on ONE fresh build per case
  (`test_recompute_laws`, in `test_law_crossings.py`) instead of one fresh
  build per law.
- **The hole rule (author).** A hole maps backward onto the same side of the
  BRIDGE it is punched in (I1), not the same side of the stable manifold.
  `check_holes_share_bridge_side` (which delegated to the library's
  `bridge_side_violations`) is replaced by
  `check_backward_holes_keep_bridge_side`: every hole of an origin carries the
  reference hole's `bridge_side` (parity-flipped under an orientation-reversing
  map), and the reference hole's coordinates carried back `|iterate|` steps by
  the real inverse map lie on that side of the propagated hole's own bridge
  (an independent geometric oracle). Runs on every case with holes (k10, k28
  one/two blasts, p3, nested; inversion has none). `check_openings_linked_bound`
  (a propagated hole keeps its origin's `(which, row)` at a linked bound) is
  DELETED: it asserted row preservation, which the author's rule says is not a
  law (closes E7). No law compares stable-manifold rows along an orbit;
  `openings_on_own_bridge_row` only checks each hole against ITS OWN bridge.
- **Wiring guard** `test_every_law_runs_in_the_tier` rewritten: every
  `check_*` of `helpers.laws` is in exactly the layers, every layer has its
  `layer_test`, every layer with known issues has its `known_issue_test`, and
  every `KNOWN_ISSUES` / `NOT_APPLICABLE` key names a real law.
- `cases.law_params` removed (unused); `NOT_APPLICABLE` keys renamed /
  dropped with the two laws above.

### Deletion / change ledger

- `tests/invariants/test_law_*.py::test_<law>[<case>]` (446 ids: 407 pass,
  31 skip, 8 xfail) -> `test_<layer>_laws[<case>]` (11 layers x 6 cases) +
  `test_recompute_laws[<case>]` (6) + `test_<layer>_known_issue[<law>-<case>]`
  (7) + `test_orientation_reversing_case_builds` + `test_every_law_runs_in_the_tier`
  = 81 ids (73 pass, 8 xfail). Every check still runs on every case it ran on.
- Deleted law `openings_linked_bound` (row preservation along the backward
  orbit; author's hole rule). Its `NOT_APPLICABLE` entries (k10,
  k28_two_blasts, inversion) removed.
- Replaced law `holes_share_bridge_side` -> `backward_holes_keep_bridge_side`
  (explicit I1 against the reference hole + geometric carried-point oracle).
- `tests/facade/test_session_caches.py::_alternate_pip` picks the
  anchor-nearest alternate candidate by stable cdist instead of `alternates[0]`
  (registry-id order is not reproducible between builds, so the old choice
  flipped `DualGraph._pip_segment`'s no-registered-iterate branch in and out of
  coverage and failed the coverage guard on one run).

### Verification

- `env/bin/python -m pytest -q -rxX`: 751 passed, 1 deselected, 11 xfailed in
  90.5 s (`followup_run.txt`); the 11 xfails are exactly `KNOWN_ISSUES`
  (8 law tier + 3 `open_p3_deep_runs`), no XPASS.
- Coverage guard vs `cov_base.json`: OK (`cov_followup.json`).
- Law-tier wall time: 36.0-38.6 s over three runs, against 36.3-40.9 s for
  the per-law tier on the same machine (not grown).
- A probe layer with one failing and one vacuous check raised one error naming
  both (`bad: AssertionError: assert 1 == 2 [...]`, `empty: checked nothing`).

### Deviations

1. CLAUDE.md line "`invariants/` (physical laws, one test per law x case ...)"
   still says per law: CLAUDE.md is the author's to edit; it should read "one
   test per layer x case; a known issue has its own xfail".
2. The `_alternate_pip` fix in `facade/` is outside the law tier; it was needed
   to make the coverage guard deterministic.

## Author follow-up (2026-10-05): golden tier removed

The author removed the golden tier entirely ("incredibly arbitrary and
irrelevant"). The CLAUDE.md fixture facts stay in CLAUDE.md as documentation;
no test pins a case-specific word, element name, transition matrix, letter or
count of a real build. `is_reliable`, `verified` and `source` are asserted
nowhere (the only remaining `.source` assertion is `ElementLanding.source`, the
landing's homotopy element, which is a different attribute).

### What changed

- **`tests/golden/` deleted** (3 files, 21 cases).
- **`pyproject.toml`**: the `golden` marker is no longer registered; `perf` is
  the only marker (`--strict-markers` would now reject a stray `golden`).
- **`tests/conftest.py`** docstring: the `golden/` tier line and the `golden`
  marker removed; a paragraph states the no-case-fact rule (a specific word or
  name appears only in a synthetic test whose input was built by hand to
  produce it).
- **`tests/helpers/names.py`**: the golden-only helpers `brackets`,
  `refined_children` and `class_by_pair` deleted, and the golden-only
  `short=` parameter dropped from every speller (always the full
  `ElementName.text`). What is left serves `facade/test_session_equivalence.py`
  (session product = direct build), which compares two builds of one case and
  pins nothing. Module docstring rewritten to say so.
- **`tests/cases.py`**: `CaseExpect.has_image_bridges` deleted (a per-case fact
  that no test read once golden was gone; k10 True, p3 False, nested True).
  The remaining `CaseExpect` fields are construction parameters (periods,
  inversion, number of tangles, blasts).
- **Scattered pins turned into rules** (found by grepping for letters, names,
  matrices and counts on real builds):
  - `facade/test_session_bridge_classes.py::test_active_class_lies_in_the_zone_after_trimming_at_the_image_pip`
    pinned the k=10 (10-step) table: 2 classes, letter `"a"`, 4 members,
    1 loop. Now: there is an active class; every class's `zone_key` is the
    zone's or None; a class is lettered iff active; every active class lies
    in the zone; a class holding a loop is inert; some class is recorded
    exterior (`zone_key is None`).
  - `plotting/test_plot_topology.py::test_p3_circle_follows_the_zone_boundary`
    pinned the p3 ring as clockwise with every interior on the right and the
    orbit order 0 -> 1 -> 2. Now: the ring holds each stable branch of the
    cycle once, interiors are `left`/`right`, and the circle keeps whatever
    cyclic order and direction `cartoon_zones()` reports, interiors inside.
  - `unit/topology/test_element_naming.py::test_p3_names_carry_orbit_codes_without_a_letter`
    pinned the codes `{"0.0", "1.0", "2.0"}`. Now: the codes are exactly the
    `branch_code` of each key in `fp.branch_cycle("stable")`, `k_value` of
    them.
  - `unit/topology/test_bridge_class.py::test_unresolved_loop_stands_alone_and_warns`
    (k10 build) pinned the k=10 class count `len(table) == 3` (arrived with
    the Phase 10a move, missed by the first grep; found by the verifier).
    The line is deleted: the rule is already asserted (the loop now stands
    alone as an inert `BridgeClass(x, x)` with no `folded_from`, a WARNING is
    logged, and the class it used to fold into is now active). Its k10
    precondition (an inert class holding a folded loop) and that of
    `test_loop_folds_into_its_preimage_class` are now looked up among all
    inert classes (not `inert[0]`) and the tests skip when the build has none.
  - `unit/topology/test_symbolic_dynamics.py::test_is_reliable_and_describe`
    renamed `test_describe_lists_every_class` (it never asserted
    `is_reliable`; the name was a leftover).
- **Docstrings pointing at `tests/golden/`** rewritten in
  `facade/test_session_symbolic_dynamics.py`, `unit/topology/test_dual_walk.py`,
  `unit/topology/test_bridge_class.py`, `unit/topology/test_element_naming.py`,
  `unit/topology/test_symbolic_dynamics.py`.
- **CLAUDE.md** "Run tests": `-m golden` line and `golden` in the tier list
  removed; the tier paragraph drops `golden/`, says the law tier is one test
  per layer x case (closing deviation 1 of the previous follow-up), and states
  the no-case-fact rule. The fixture facts elsewhere in CLAUDE.md are
  untouched.
- **E8 CLOSED** (removed by the author): the p3 names, the k10 folded loop,
  the k28 one-blast "0 image bridges" and the "exterior class" pins no longer
  exist. The Phase 8 deviation 4 item (p3 ring-order fact) is closed the same
  way.

Kept on purpose: synthetic tests that spell words, names or matrices
(`test_symbolic_dynamics.py`'s hand-made layout, `test_element_naming.py`'s
hand-built families, `test_bridge_class.py`'s hand-built references): those
outputs are the rule applied to an input the test constructed, not a fact of
a real build. `test_nested_names_carry_fixed_point_letters` keeps
`{outer.label, inner.label} == {"A", "B"}` (the workbench stamps `A, B, ...`
in construction order, a rule).

### Deletion / change ledger

| Deleted id | Where its rule is checked now | Coverage lost |
|---|---|---|
| `golden/test_golden_k10.py::test_k10_iterated_rows` | law `iterated_unique_owner`, `names_agree_with_structure`; naming rules in `unit/topology/test_element_naming.py` | none |
| `golden/test_golden_k10.py::test_k10_classes` | laws `anchor_bridge_class_leads_its_tangle`, `class_orientation_anchor_outward`, `only_active_classes_lettered` | none |
| `golden/test_golden_k10.py::test_k10_inert_class_rests_on_its_folded_loop` | `unit/topology/test_bridge_class.py` loop folding / inertness on hand-built tables | none |
| `golden/test_golden_k10.py::test_k10_itinerary_and_word` | symbolic laws `itineraries_even`, `itinerary_pairs_same_side`, `itinerary_pairs_are_classes` | none |
| `golden/test_golden_k10.py::test_k10_refinement` | laws `refined_children_inherit_word`, `member_matching_consistent` + refinement on the synthetic layout in `test_symbolic_dynamics.py` | none |
| `golden/test_golden_k10.py::test_k10_transition_matrix` | law `matrix_is_token_counts` | none |
| `golden/test_golden_k10.py::test_k10_evidence` | (asserted `verified`/`is_reliable`; dropped by policy) | none |
| `golden/test_golden_k28.py::test_k28_one_blast_classes_and_word` | class + symbolic laws on `k28_one_blast` | none |
| `golden/test_golden_k28.py::test_k28_one_blast_exterior_class_is_inert_through_a_virtual_loop` | inertness rules in `unit/topology/test_bridge_class.py` | none |
| `golden/test_golden_k28.py::test_k28_one_blast_has_no_image_bridges` | none (a case count) | none |
| `golden/test_golden_k28.py::test_k28_two_blasts_classes_in_table_order` | laws `class_table_order`, `anchor_bridge_class_leads_its_tangle` | none |
| `golden/test_golden_k28.py::test_k28_two_blasts_words` | symbolic laws on `k28_two_blasts` | none |
| `golden/test_golden_k28.py::test_k28_two_blasts_refinement` | laws `refined_children_inherit_word`, `member_matching_consistent` | none |
| `golden/test_golden_k28.py::test_k28_two_blasts_transition_matrix` | law `matrix_is_token_counts` | none |
| `golden/test_golden_k28.py::test_k28_two_blasts_cut` | laws `iterated_child_inside_parent`, `iterated_keeps_homotopy_boundaries`, `iterated_unique_owner`, `iterated_cut_provenance` + `unit/topology/test_iterated_partition.py` | none |
| `golden/test_golden_k28.py::test_k28_two_blasts_evidence` | (asserted `verified`/`is_reliable`; dropped by policy) | none |
| `golden/test_golden_period3.py::test_p3_words_form_the_orbit_shift_chain` | symbolic laws on `p3` | none |
| `golden/test_golden_period3.py::test_p3_is_closed` | none (a case count) | none |
| `golden/test_golden_period3.py::test_p3_words_survive_four_blasts` | `regression/test_nested_own_blast.py` (inner words = p3 words, relational) | none |
| `golden/test_golden_period3.py::test_p3_resolves_and_is_reliable` | (asserted `is_reliable`; dropped by policy); resolution checked by symbolic laws | none |
| `golden/test_golden_period3.py::test_nested_resolves_everything_after_two_outer_blasts` | symbolic laws on `nested` | none |

Changed (not deleted): the five tests listed under "Scattered pins turned
into rules" above; `helpers/names.py` (3 helpers and the `short=` parameter
removed); `cases.CaseExpect.has_image_bridges` removed.

### Verification

- `env/bin/python -m pytest -q -rxX`: 730 passed, 1 deselected, 11 xfailed in
  86.8 s (751 -> 730: the 21 golden cases). The 11 xfails are exactly
  `KNOWN_ISSUES` (8 law tier + 3 `open_p3_deep_runs`), no XPASS.
- Coverage guard vs `cov_base.json`: OK (`cov_followup.json`). Against the
  previous follow-up's `cov_followup.json` (HEAD): total 90.24 % -> 90.24 %,
  and no line of any module changed state (zero newly missed, zero newly
  hit): everything the golden tier exercised is exercised elsewhere.
- Re-verified after the count-pin fix: 730 passed, 1 deselected, 11 xfailed
  (exactly `KNOWN_ISSUES`, no XPASS) in 86.1 s; coverage guard OK. A second
  grep for `len(...) == <int>`, nested-list matrices, letters and element
  names finds only synthetic inputs, construction parameters (nested has two
  fixed points) and rule counts.
- `grep -rn golden tests/` finds only the conftest sentence recording the
  removal.

### Deviations

1. CLAUDE.md's tier paragraph also had the stale "one test per law x case";
   since the paragraph was being edited for golden anyway, it now reads "one
   test per layer x case; a known issue has its own xfail" (deviation 1 of
   the previous follow-up).
2. The colour-family-per-tangle test (`TANGLE_COLOR_FAMILIES` /
   `HETEROCLINIC_COLORS` as a rule) is not part of this change.

## Author follow-up (2026-10-05): cache freshness, not hit/miss policy

`tests/facade/test_session_caches.py` no longer pins whether an event is
served from the cache or rebuilt. For every product (`trellis`,
`arrangement`, `bridge_classes`, `minimal_trellis`, `iterated_partition`,
`dual_graph`, `symbolic_dynamics`) x every event (`no_change`, `rebuild`,
`generation_bump`, `partition_signature`, `pip_change`) it asserts FRESHNESS:
after the event the session's product equals, through a structural
letter-free view, a direct cache-free build from the current state. Object
identity is asserted only for `no_change` (same object) and `rebuild=True`
(the product itself is a new object). Every product is asked for once before
the event, so every cache holds an entry the event could make stale.

### What changed

- New `tests/helpers/direct.py`: `DirectBuild` (moved from
  `test_session_equivalence._Direct`, now with a fixed-point selection and
  its classes lettered from a FRESH `BridgeAlphabet` instead of the session's
  table, so the direct dynamics reads no session cache), the structural views
  (moved unchanged), `view_of`, `direct_view`, `PRODUCT_NAMES`.
- `tests/facade/test_session_equivalence.py` uses the helper; same 16 cases,
  same assertions.
- `tests/conftest.py` layout docstring updated (facade = cache freshness;
  helpers include the direct build).

### Deletion / change ledger

| Old | New | Note |
|---|---|---|
| `test_session_caches::test_expected_table_is_complete` | deleted | the `EXPECTED` hit/miss table is gone |
| `::test_session_cache_table[<product>-<event>]` (35) | `::test_session_product_is_fresh_after_every_event[<product>-<event>]` (35) | event `hit` renamed `no_change`; HIT/MISS pins replaced by freshness; reuse only on `no_change`, new object only on `rebuild`. The `generation_bump` / `pip_change` events no longer assert the partition signature stays put (they assert only that the generation moved / the pip changed) |
| `::test_trellis_misses_after_every_mutation_path[*]` (5) | `::test_trellis_is_fresh_after_every_mutation_path[*]` (5) | `fresh is not first` replaced by trellis view = `Trellis.from_workbench` view |
| `::test_trellis_misses_after_a_recompute_on_the_nested_session` | `::test_trellis_is_fresh_after_a_recompute_on_the_nested_session` | registry identity kept, plus view = direct build |
| `::test_cache_is_kept_per_fixed_point_selection[*]` (7) | `::test_every_fixed_point_selection_is_fresh[*]` (7) | identity pins replaced by: interleaved selections each equal a direct build of that selection |
| `::test_bridge_class_letters_survive_rebuilds_and_repartitions` | unchanged | the alphabet rule, not cache policy |

Observability (measured on `k10`, recorded here, not asserted):
`partition_signature` changes every view from `bridge_classes` on and
`pip_change` changes `dual_graph` and `symbolic_dynamics`, so a stale cache
there fails the test (checked by monkeypatching the session to ignore the
partition signature and the pips: 6 cells fail). `generation_bump`
(`iterate_bridge`) and the growth / `iterate_bridge` / `rebuild_bridges`
mutation paths register no crossing on `k10`, so their cells are consistency
checks that a stale cache would also pass; `add_resonance_zones`, `restore`
(the 2026-07 regression) and the nested recompute do change the trellis.

### Verification

- `env/bin/python -m pytest -q -rxX`: 729 passed, 1 deselected, 11 xfailed in
  87.5 s (730 -> 729: the deleted table-completeness test). The 11 xfails are
  exactly `KNOWN_ISSUES`, no XPASS.
- `test_session_caches.py` alone: 49 passed in 3.95 s (was 50 in 3.72 s).
- Coverage guard vs `cov_base.json`: OK (`cov_followup.json`); total 90.24 %
  -> 90.24 %.

### Deviations

1. A view-changing cheap generation bump was not found: one more unstable
   growth step + recompute broke the arrangement's bridge/branch agreement
   and a trim at an inner pip + re-partition tripped
   `check_bridge_rows_consistent` (both test-recipe problems, not pursued);
   the `iterate_bridge` event was kept and documented as unobservable.

## Author follow-up (2026-10-06): the period-3 high-stretch test made non-vacuous

`regression/test_high_stretch_growth.py::test_period3_high_stretch_growth_is_not_scrambled`
was vacuous (Phase 7a finding): four unstable steps from the default seed left
7 nodes per branch with adjacent relative cdist gaps of ~3e15 eps, nowhere near
the one-ULP collisions it guards (cdist-strict-monotonicity-fix memory).

### Re-deriving the parameters

Measured on k=2.1, `P3_ORBIT_SEED`, `area_cutoff = 1e-7` (minimum over branches
of `gap / cdist` between adjacent nodes, in units of `np.finfo(float).eps`):

- Default seed (`max(fixed_point.accuracy, 5e-6)`; the fsolve orbit makes
  `accuracy` tiny): unstable 9 steps -> 2.5e11 eps (8451 nodes); step 10 is
  the escape (~1.5e6 nodes); stable 11 steps -> 3.9e13 eps. Unreachable:
  growth scales a gap and its cdist by the same factor, so only refinement
  halves the relative gap, and the arm escapes long before ~52 halvings.
- Coarser seed via `fp3.accuracy` (the documented dynamic seed term):
  `2e-3`/`3e-3` reach it only on some branches; `5e-3` and `7e-3` reach
  0.5-0.8 eps (adjacent cdists one ULP apart) on every stable branch in
  5 steps, 0.01-0.03 s, with zero hairpins on the current code; `1e-2`
  reaches it too but the current code already shows 4-15 hairpins there
  (near-duplicate seam points), so it was not used.
- Mechanism: with a 5e-3 seed the fundamental segment is curved enough that
  the images of its two ends land ~1e-5 apart in space while their cdists
  agree to one ULP: the same "two distinct points the cdist cannot order"
  collision, produced at the seed seam rather than at a high-stretch fold.

### What changed

- `test_period3_high_stretch_growth_is_not_scrambled` ->
  `test_period3_near_ulp_growth_is_not_scrambled`: k=2.1 period 3,
  `fp3.accuracy = 5e-3`, `area_cutoff = 1e-7`, 5 STABLE steps (was 4
  unstable steps from the default seed). New non-vacuity assertion per branch:
  `min(gap / cdist) < 4 * eps` (`NEAR_ULP_EPS`). Kept: no spike, cdist
  non-decreasing; added: no hairpin reversal. Runtime 0.02 s.
- `helpers/invariants.py`: new `find_reversals` / `assert_no_reversals`
  (a joint whose consecutive segments have turning cosine below -0.5, i.e.
  two nearby points stored out of order). The spike check alone has no
  teeth here: a seam scramble swaps points ~1e-5 apart, far below the
  30x-median segment threshold.
- `test_k10_growth_keeps_geometry_smooth` also asserts `assert_no_reversals`
  (passes; zero hairpins on k=10).

### Teeth (scratch worktree at HEAD, src edited there only)

Re-introduced the old strict tie-breaking in `ManifoldMachine.py` (the
strictify code was never committed, so it was re-written from the memory note):
(1) nudge tied image cdists up by `nextafter` in `iterate_manifold`;
(2) in `merge_manifolds`, a collision of two distinct points nudges the second
up by `nextafter` and keeps both instead of deduplicating; (3) the
`representable` guard in `_refine_layer` (refine only where `c0 < mid < c1`).

| build | old test (4 unstable steps) | new test |
|---|---|---|
| all three edits | passes (vacuous) | FAILS: "22 hairpin reversal(s) on the 'stable' manifold ... turning cosine -1.000" |
| (1) image strictify only | - | passes |
| (2) merge nudge only | - | FAILS (same 22 reversals) |
| (3) representable guard only | - | passes |
| current code | passes | passes |

So the test guards the merge's exact-tie dedup (the part of the fix that
matters in this regime); edits (1) and (3) produce no observable scramble in
5 steps from this seed (the scan above showed 20-30 hairpins per branch for
the full broken build at `accuracy` 3e-3..7e-3, all from the merge nudge).

### Deletion / change ledger

| old | new | why |
|---|---|---|
| `regression/test_high_stretch_growth.py::test_period3_high_stretch_growth_is_not_scrambled` | `::test_period3_near_ulp_growth_is_not_scrambled` | re-derived parameters reach one-ULP adjacent cdists; non-vacuity + no-hairpin assertions; fails on the strictify build |
| `::test_k10_growth_keeps_geometry_smooth` | unchanged name | + `assert_no_reversals` |

### Verification

- `env/bin/python -m pytest -q -rxX`: 729 passed, 1 deselected, 11 xfailed in
  85.7 s; the 11 xfails are exactly `KNOWN_ISSUES`, no XPASS.
- Coverage guard vs `cov_base.json`: OK (`cov_followup.json`); total
  90.24 % (unchanged; src untouched).
