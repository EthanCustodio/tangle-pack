# Test-suite refactor plan (tanglepack)

## Context

The suite has 796 collected cases (~640 functions, 64 files, 18.3k lines) and runs in ~45 s. Speed is not the problem. The author worries the suite overtests and enforces things it shouldn't.

A read-only audit read every test: 14 catalog agents, then an infrastructure agent, a synthesis agent and a critic. Reports (scratchpad; copy into the run folder in Phase 0):
- `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/d02d75de-a9cd-403d-9360-8fb8fc84071e/scratchpad/{synthesis,critique,infra}.md`
- the per-test catalog: `~/.claude/projects/-home-dezhu-Development-TanglePack-tangle-pack/d02d75de-a9cd-403d-9360-8fb8fc84071e/subagents/workflows/wf_3e805c13-402/journal.jsonl`

**What the tests are:**

| Kind | Share |
|---|---|
| Algorithm rules | 52% |
| Physical laws | 9% |
| Implementation details | 9% |
| Session layer / caches | 8% |
| Fixed-bug regressions | 7% |
| Plotting | 6% |
| Exact Hénon-case results | 5% |

**What is wrong:**
- About 40 cases only check that old code is gone or where code lives.
- About 8 cases freeze the output of a finished migration.
- About 75 cases pin wording: `describe()`, log and exception text.
- About 45 cases reach into private internals.
- About 35 cases pin exact case results, often the same fact in 3 to 6 files.
- About 14 cases enforce rules the author marked provisional.
- About 20 cases freeze plot geometry or style.
- There are 26 clusters where several files check the same fact; the k=10 build recipe is copied 9 times.
- `slow` and `regression` markers do not match reality and select nothing.
- `tests/visualizations/` is dead.

**Outcome:** about 600 cases in a tiered layout.
- Each physical law is checked on every case, about six cases in all.
- Each of the author's fixture facts is pinned once and changes only with the author's sign-off.
- No pins on wording, style, private internals or provisional rules.
- Known bugs are recorded as `xfail(strict=True)`.

## Author decisions (binding)

**1. Scale.** Goal is "right things", ~600 cases is fine.

**2. Golden tier: what is pinned and how.**
- CLAUDE.md fixture facts are pinned ONCE in `tests/golden/`, marked `@pytest.mark.golden`, runs by default.
- Author-gated: no re-record mechanism, and every change needs the author's sign-off.
- Words are spelled in element names / oriented element pairs, NEVER letters.
- `is_reliable` / `verified` are asserted only in golden. `source == "walk"/"trellis"` is asserted nowhere, not even k28 one-blast.

**3. Provisional rules.** Test only the firm part. The provisional part goes into the Dev Notes. Provisional rules are:
- the `+1..+(k-1)` forward-hole exemption (the firm part is that iterate ≥ k_value is never punched);
- the mean-of-neighbours refined cdist (replace with "strictly between neighbours");
- the anchor faces-closed limitation.

**4. Deletion policy.**
- **Tombstones and code-layout tests are deleted.** KEEP the behavioural one-rule guards: `test_collect_is_the_only_traversal`, and `test_workbench_has_no_private_key_advance`, upgraded to a spy showing that `infer_iterate_table`, `iterate_bridge` and `Trellis.image_cdist` all reach `FixedPoint.advance_key`.
- **Private internals:** only pure kernels with no public route; rewrite the rest through the public API.
  - Allowed kernels: `_side_of`, `_resolve_inertness`, `_backward_endpoint`, `_bridge_unstable_span`, `_containing_bridge`, `_near_far`, `_collapse_noise_crossings`, `_do_segments_intersect`/orientation predicates, `_parabolic_fit`, `_curvature_area_batch`, the iterated-cut chord kernel.
- **Migration snapshots are deleted:**
  - the frozen pre-1.8 grow loop;
  - the 6-digit bridge-cutting cdists;
  - zone areas pinned to 1e-11;
  - the k10 graph counts 8/18/7/7/4;
  - "exactly 3 init points".
- **Dead API (code stays, listed in the Dev Notes): drop the tests only.** This covers `on_interval`/`_get_lambda_u`, `invalidate_trellises`, `BoundaryArc`, scalar `_curvature_area` and the `forward_*_branch_cycle` wrappers.

**5. Assertion style.**
- **Plotting:** smoke plus topological properties only. No angles, fonts, colours, z-order, palette order or legend wording. Keep "`class_colors` never repeats".
- **Text:** exceptions are checked by type only, logs by level + logger only, and `describe()`/`summary()`/repr only for being non-empty and mentioning each class or fixed point.
- **cdist:** the helper defaults to `strict=False`. `strict=True` only in `test_low_stretch_growth_keeps_cdist_injective`; the no-spike check runs everywhere.
- **Tolerances:** keep them, stated against library constants (`collision_rtol`, `SCALING_RTOL`, `area_cutoff`, `_MIN_SEED_STEP`) and using `pytest.approx` instead of exact `==`.
  - Do NOT loosen the `collision_rtol` unlinked halves: they guard the 2026-09-15 collision bug.

**6. Caches.** One table-driven suite: product × {hit, rebuild, gen bump, partition-signature change, pip change}, plus one check per product that the session result equals a direct build. Remove cache-identity asserts from the fact tests.

**7. Fixtures.**
- Builders live in a tests-only `tests/cases.py`. Their parameters are frozen there and never read from `henon_cases` defaults.
- Fixtures are function-scoped. **Exception:** the physical-law tier uses one read-only build per module per case, protected by a fingerprint guard (below).

**8. Layout.**
- Tiered directories: `invariants/`, `unit/{numerics,topology,loom}/`, `facade/`, `golden/`, `regression/`, `plotting/`.
- Consolidate and delete first, then move files once.
- Markers: drop `slow` and `regression`; add `golden` (runs by default) and `perf` (opt-in, `addopts -m 'not perf'`).

**9. Cases in the law tier.**
- The cases are k10, k28 one blast, k28 two blasts, p3, nested, and the **k10 inversion saddle (−2.3166, 2.3166)** through the FULL pipeline. Whatever fails on inversion today becomes a `xfail(strict=True)` known issue.
- Add a placeholder **b = −1 orientation-reversing** case. It is `xfail(strict=True, raises=ValueError)` until det J<0 is supported. Note that the inversion saddle has det J = +1 and is not orientation-reversing.

**10. Known bugs, recorded as `xfail(strict=True)`; fixing them is out of scope.**
- **One anchor per unstable branch.** An inversion point has two unstable branches, each attached to its one anchor. Today the code gives two (0,0) anchors per branch, so this is a bug; it gets an xfail on inversion. The source fix is a separate follow-up.
- **The deep p3 runs.** At 15–16 steps or 6+ blasts there is an unreachable class and a virtual `new1`.

**11. Hole quirk.** The p3 and k28 snapshots show a 'left' hole at iterate < 0 opening only ('outward','right'). INVESTIGATE before deleting those snapshots (Phase 0, P1).

**12. Extras.**
- **GPU parity** runs by default and skips without CuPy or a device. Make the no-CuPy error test actually run by monkeypatching `sys.modules`.
- **Production-check wiring:** assert that `punch_holes` runs its I1/I2 checks, instead of repeating those checks in tests.
- **Coverage gap tests.**
- NOT included: a smoke test of the figure scripts.

**13. Critic corrections stand.** These are KEPT:
- `test_direct_hole_side_is_the_side_of_its_coordinates` (becomes a law, with nested added);
- `test_henon_propagated_hole_side_matches_coords`;
- the `Intersection.fixed_points` IndexError tests;
- the synthetic label-slot test;
- `test_bridge_endpoints_stay_on_the_bridges_own_branch`, which must run on p3/nested;
- `test_iterating_fixed_point_bridge_returns_existing_copies`.

Strengthened: `test_kevin_way_period_three_orbit_chain` asserts the chain order (unstable 0,1,2; stable 0,2,1).

**14. Ledger and notes.** `docs/test_suite_refactor_ledger.md` records every deleted test and why. The `tests/conftest.py` docstring holds the standing Dev Notes: the untested provisional rules, the dead API and the `KNOWN_ISSUES` index.

## Rules for every phase

- **One phase = one coder commit.** The commit body carries that phase's deletion ledger, which is also appended to the ledger file.
- **Every deleted test is checked first** against the memory notes for the bug it guards (`codebase-audit-2026-09`, `pseudoneighbor-collision-fragility`, `pseudoneighbor-partition`, `numerics-test-suite`, `regions-refactor-2026-09`, `cdist-strict-monotonicity-fix`, `hole-side-own-blast-2026-10-02`, `audit-polish-2026-07`).
- **A test whose assertions move** is deleted in the same commit that adds their new home.
- **Verify each phase:**
  - `env/bin/python -m pytest -q -rxX --durations=15 --cov=tanglepack --cov-report=json:$RUN/cov_pN.json`
  - `env/bin/python tests/_tools/coverage_guard.py $RUN/cov_base.json $RUN/cov_pN.json`: fails on any newly missed `src/` line outside the whitelisted dead-API qualnames.
  - Record the collected count.
  - Run 5 random node ids in isolation, always including `test_dual_graph.py::test_k10_a_different_pip_moves_the_unified_set` (known to fail when run alone).
  - The `-rxX` listing must show only `KNOWN_ISSUES`.

## Phases

### Phase 0: baseline and analyst probes (no test edits)

**Baseline.** In `.refactor/runs/2026-10-xx-test-suite/`, save:
- the node-id list, `cov_base.json` and durations;
- copies of the audit reports.

**Add `tests/_tools/coverage_guard.py`.** It compares each module's `missing_lines`, with whitelist spans taken from `ast`.

**Probes.** Each is a scratch script; findings go into ledger §0.

| Probe | What it does | Output |
|---|---|---|
| P1 hole quirk | For each hole in p3 and k28 two blasts: `bounding_ids`, whether the near bound is an anchor, the origin's openings against its backward images' openings, and what each missing half faces | (a) every difference is an opening dropped at an anchor or tail → an openings law. (b) rows really flip → `xfail` `test_propagated_holes_open_their_origins_pair`. Note: the k28 iterate-0 left hole opens the LEFT pair while iterates −1 and −2 open RIGHT, which is suspicious |
| P2 inversion | Run the full pipeline stage by stage | The first failing stage, a law-by-law result and the anchor count; these become the `KNOWN_ISSUES` draft |
| P3 deep p3 | Run 15 and 16 steps and 6 blasts | The exact xfail assertions |
| P4 b = −1 | Choose (k, −1) and a seed | Confirm that `construct_fixed_point` raises `ValueError` |
| P5 regression-guard audit | Check every deletion candidate against the memory notes | The Guard column |
| P6 builder parity | Repoint the fixtures at the planned builders | Which parameters matter. The k10 7+2 vs 9-step recipes and `infer_iterate_table` on/off differ; nested has p1 cutoff 1e-4 vs 1e-7 and repin vs none |
| P7 k28 cut | Translate the CLAUDE.md cut facts (`R_1^1=[0,10]`…) | Id-free relations |
| P8 runtime | Time each law-case build | Timings |

### Phase 1: infrastructure (additions only; same node ids, green)

**`tests/cases.py`:**
- `Case` dataclass: session, fixed points (outermost first), pips, expectations.
- Builders: `build_k10`, `build_k28(blasts=1|2)`, `build_period3`, `build_nested`, `build_inversion`, `build_orientation_reversing`.
- `LAW_CASES`, a frozen `BUILDERS` mapping, `KNOWN_ISSUES`, `NOT_APPLICABLE`, and `law_params(law_id)`, which applies the xfail/skip marks.
- Delegate p3 and nested to `henon_cases` with every keyword explicit; reimplement k10, k28 and inversion.

**`tests/helpers/`:**
- `invariants.py`: moved from `tests/invariants.py`, default `strict=False`, spike message fixed.
- `laws.py`: `check_<law>(case) -> int` returns the number of items checked, grouped into layers.
- `fakes.py`: from `walk_helpers.py`.
  - `bare_fixed_point()` returns a REAL `FixedPoint`.
  - `FakeTrellis` borrows `Trellis.image_cdist`.
  - `make_result` requires explicit owners.
- `names.py`: letter-free spellers `word_in_names`, `brackets`, `matrix_by_classes`.
- `logs.py`: `assert_logged(caplog, level, logger)`.
- Shims at the old paths until Phase 10.

**`conftest.py` and config:**
- `conftest.py`: Agg backend; every fixture becomes a thin call to a builder, and all are function-scoped.
- `pyproject.toml`: register `golden` and `perf`.

### Phase 2: pure deletions (~75 cases)

Delete:
- `tests/visualizations/`.
- Whole files: `numerics/test_workbench_split.py`, `test_manifold_machine.py`.
- Tombstones:
  - `test_deleted_genealogy_attributes_are_gone`, `test_grown_until_intersection_is_gone`, `test_tangle_keeps_only_index_state`;
  - `test_string_dispatched_iter_method_is_gone`, `test_stability_alias_is_defined_once_in_numerics`, `test_no_test_or_script_defines_its_own_henon`;
  - `test_fixed_point_takes_no_branch_count`, `test_boundary_arc_dataclass_is_gone`, `test_shapely_is_not_a_dependency`;
  - the 4 plotting layout tests (`test_every_plotter_is_a_module_function_in_plotting`, `test_stable_partition_module_holds_no_drawing_code`, `test_hole_style_conventions_live_once`, `test_session_exposes_one_fanout_helper_per_shape`).
- Migration snapshots:
  - `test_bridge_cutting_pin_k10` and `test_bridge_cutting_pin_period_3`;
  - `test_per_step_beta_matches_the_expression_it_replaces`;
  - `test_k10_intersection_graph_is_unchanged`, `test_small_tangle_crossing_count_is_unchanged`;
  - the 2 zone-area tests (k10 and p3);
  - `test_initialization_unstable`/`stable`.
- Dead-API tests: `_get_lambda_u` ×4, `on_interval`, `test_curvature_area_batch_matches_scalar`, `test_invalidate_trellises_is_a_deprecated_no_op`, `test_branch_cycle_reproduces_the_topology_branch_orders`.
- Tautologies and duplicates:
  - `test_image_chain_is_the_bridge_class_function`, `test_kinds_are_distinct`, `test_names_are_deterministic_across_builds`;
  - `test_session_result_matches_gathered_partitions_helper`, `test_session_bridge_classes::test_p3_no_class_mixes_fixed_points` (keep its twin in test_bridge_class);
  - 3 tests in `test_image_cdist`, but KEEP `test_partition_elements::…falls_back_when_an_iterate_is_missing` and `test_image_cdist::…covers_both_endpoint_images`;
  - `test_partition_elements::test_image_of_element_lands_on_the_advanced_branch`;
  - `test_henon_reference_holes_hug_the_stable_manifold` and the `_stable_arc_midpoint` ×2 tests;
  - `test_henon_partition_still_built_after_wiring`.
- From `test_forward_pairs_beyond_the_fundamental_segment_get_no_hole`, the 1..k−1 parameters (provisional).

In the same commit, replace the summary/repr/describe string tests with one non-empty smoke test.

### Phase 3: loosen in place (3a numerics + regression, 3b topology, loom, session)

- Remove every `match=` on exceptions (the full list is in the planner's §A, Phase 3).
- Every caplog check goes through `assert_logged`.
- Bring `describe()` assertions down to "non-empty plus mentions".
- Remove `is_reliable`/`verified`/`source` from all non-fact tests.
- Remove the cache-identity asserts from the fact tests (`dyn is session.symbolic_dynamics(...)`, `_bridge_classes` sizes).
- Tolerances:
  - use `approx` against the constants;
  - `min_unstable_cdist == 0.0` → `approx(0, abs=cdist_tol)`;
  - switch `strict=True` to `strict=False` in `test_machine_iterate_invariants` and `test_growth_integration`.
- Rewrite private access through the public API: `_map_batchable`, `_built_generation`, `_partial_bridges`/`_bridges`, `_intersection_registry`, segment ownership.
- `regression/test_near_vertical_refinement` uses `_curvature_area_batch`.

### Phase 4: physical-law tier (`tests/invariants/`)

**Structure.**
- One test per (law × case). Each module has a module-scoped `law_case` fixture that builds the case once with every topology step done inside the fixture.
- An autouse **fingerprint guard** compares before and after each test: generation, holes, partition signature, pips, alphabet, bridge ids. A test that mutates the shared build fails.
- Laws that mutate (`preserve_ids`/idempotent recompute) use fresh builders instead.

**Marks and non-vacuity.**
- Every law asserts that it checked more than 0 items, unless the `(law, case)` pair is listed in `NOT_APPLICABLE` (a skip with a reason).
- Known issues get `xfail(strict=True)` from `KNOWN_ISSUES`, always including:
  - one anchor per unstable branch on inversion;
  - the b = −1 placeholder.

**Files and their laws:**
- `test_law_case_sanity`
- `test_law_manifolds`:
  - non-strict monotone cdist;
  - no spikes;
  - the iterate law;
  - one-to-one.
- `test_law_crossings`:
  - every crossing is one unstable a-key and one stable b-key;
  - no u×u and no s×s, across branches and fixed points;
  - cdist bracketed;
  - sign = cross product;
  - area preserved along iterate chains;
  - `preserve_ids` / idempotent recompute (fresh builds).
- `test_law_anchors`
- `test_law_bridges`:
  - BridgeId order;
  - consecutive crossings on the bridge's own branch (the p3/nested guard);
  - `bridges_at` exact;
  - a single copy;
  - image/preimage round trip;
  - partial ⇔ id None.
- `test_law_partition`:
  - covers the branch;
  - unique owner;
  - singletons;
  - direct-hole side = coordinate side;
  - inward pair;
  - propagation lands on the predicted branch and terminates;
  - no direct hole at iterate ≥ k;
  - reference pairs valid;
  - the openings law from P1.
- `test_law_arrangement`:
  - Euler per component;
  - regions disjoint and sound;
  - `image_of` agrees with the dynamics and area (non-vacuous);
  - preimage∘image.
- `test_law_classes`:
  - every bridge in one class;
  - `anchor_outward_key` orientation;
  - row = geometry;
  - I2 on all bridges;
  - tangle grouping, then min-cdist order, connecting classes last;
  - no mixing of fixed points;
  - only active classes are lettered.
- `test_law_iterated`:
  - a child lies inside its parent;
  - homotopy boundaries are kept;
  - unique owner;
  - names agree with the structure;
  - the public minimal-trellis invariants.
- `test_law_dual_graph`:
  - wall/unified structure;
  - face side = geometry;
  - unified nodes = the pip segment on the pip's branch only;
  - bipartite degree.
- `test_law_symbolic`:
  - even, same-side pairs;
  - each pair is a class or its inverse;
  - landings contained;
  - refined children inherit the word;
  - inert classes are outside the graph and matrix;
  - matrix = token counts.

**Twins deleted in the same commit (~70 cases):**
- the k10/p3 twin pairs in `test_arrangement`, plus `test_signs_alternate_along_a_stable_branch`;
- `test_machine_iterate_invariants`, `test_growth_integration`, `test_tangle_intersection_cdist`, `test_no_same_stability_crossing`;
- the bridge-order, `bridges_at` and single-copy tests;
- the inversion growth/keys tests;
- the bridge-class row/same-branch/direction/order tests;
- the dual-graph structure/face/payload tests;
- the partition coverage/hole tests;
- the iterated and minimal-trellis invariant tests (including the test-side `_empty_stretches` reimplementation);
- the symbolic landing/itinerary-law tests.

The full list is in the planner's §A, Phase 4, step 5.

**Budget:** the law tier takes at most ~40 s.

### Phase 5: golden tier (`tests/golden/`, `@pytest.mark.golden`)

Module docstring: *"author-gated: change only with author sign-off"*.

Facts per case, letter-free (classes as oriented element pairs, anchor outward):

**k10** (`test_golden_k10.py`):
- **Rows:**
  - right: `R_1^1[ ] R_1^2( ) R_1^3[ ] R_2( ) R_3^1[ ] R_3^2( ) R_3^3[ ]`;
  - left: `L_1^1[ ) L_1^2[ ] L_2( ) L_3[ ]`.
- **Classes:** two in all, X = (R_1→R_3) active and U = (L_1→L_3) inert.
- **Itinerary:** `R_1^1 R_3^3 | L_3 L_1^2 | R_3^1 R_1^3`, a unique walk.
- **Word:** X → X U⁻¹ X⁻¹.
- **Refinement:**
  - X₁ = (R_1^1, R_3^3);
  - X₂ = (R_1^3, R_3^1);
  - both children's words are X₁ U⁻¹ X₂⁻¹;
  - every member matches.
- **Matrix:** `[[1,1],[1,1]]`, with U absent.
- **Evidence:** `verified` and `is_reliable` hold.

**k28** (`test_golden_k28.py`):
- **One blast:**
  - one active class, singleton → singleton;
  - word X U⁻¹ V⁻¹ (U and V inert trivial loops);
  - the exterior class is inert through a virtual loop.
- **Two blasts:**
  - A = (R_1→R_5), B = (R_3→R_5), C = (R_1→R_3), U = (L_1→L_3);
  - A is first;
  - words: A → A U⁻¹ B⁻¹, B → C, C → A U⁻¹ A⁻¹;
  - A refines into A₁ = (R_1^1, R_5^3) and A₂ = (R_1^3, R_5^1);
  - matrix over (A, B, C) = `[[1,1,0],[0,0,1],[2,0,0]]`;
  - cut closedness and the id-free relations from P7;
  - `is_reliable`.

**p3** (`test_golden_period3.py`):
- 3 active classes in an orbit-shift chain, X₀ → X₁ → X₂ → X₀ U⁻¹ W⁻¹;
- no image bridges;
- the same words after 4 blasts;
- `is_reliable`.

**nested** (`test_golden_period3.py`): with 2 outer blasts, everything resolves.

**Author sign-off needed (E8):** the element names for p3, the k10 "folded loop" fact, and the k28 one-blast "0 image bridges" fact.

**Deleted:** the scattered fact pins in `test_bridge_class`, `test_session_bridge_classes`, `test_element_naming`, `test_iterated_partition`, `test_symbolic_dynamics`, `test_session_symbolic_dynamics`, `test_dual_walk` and `test_higher_period_cartoon` (including `test_nested_default_words_are_pinned`, which pins a commit snapshot, not a CLAUDE.md fact), and the fact halves of the plotting tests.

### Phase 6: facade tier (`tests/facade/`)

- **`test_session_caches.py`:** a literal hit/miss table over products (trellis, arrangement, bridge_classes, minimal_trellis, iterated_partition, dual_graph, symbolic_dynamics) × events.
  - The generation bump uses a cheap `iterate_bridge`, not a regrow.
  - Extra trellis rows: growth, compute on the two-fixed-point nested session, `rebuild_bridges`, `add_resonance_zones`, restore.
  - Covers the gap: symbolic dynamics is invalidated by a signature change.
  - `rebuild=True` asserts that the product itself is rebuilt; no cascade.
- **`test_session_equivalence.py`:** session result = direct build, per product, on k10 and nested.
- **`test_session_fanouts.py`:** fan-out shapes, plus one delegation smoke test.
- **Deleted:** `test_session_trellis_cache.py`; the cache tests in the three `test_session_*` files; the cache halves of `test_session_strong_pips`, `test_arrangement` and `test_grow_until`; `test_from_results_signature_matches_the_session_signature`.

### Phase 7: unit consolidation (7a numerics, 7b topology + loom)

**Infrastructure clean-up.**
- Fakes delegate to the real rules.
- Retire `minimal_helpers.build_pieces` (37 call sites) in favour of the session API.
- Retire the 4 dual-graph assembly copies, the local `k10_session` shadows and the 2 copies of `_define_zone`.

**Merges and parametrizations:**
- the Point and BranchPoint tests, with the tautological `test_get_point` rewritten;
- geometry areas;
- registry bumps;
- the driver caps and rejections;
- the `test_fixed_points_*` tests, into one table;
- the dual_walk trivial and mirror cases;
- the blast error tests.

**Strengthened:**
- kevin-way chain order;
- period-3 solver: f³(p) = p and the eigenvalue product;
- initializer alpha for k > 1;
- the oriented polyline test checks orientation;
- the trim test asserts something is trimmed;
- the min_separation distances;
- known bridges are not re-iterated;
- the dedup test builds a real duplicate.

**Deleted:**
- `test_faces_closed_raises_at_the_cap_for_the_anchor`;
- `test_refined_cdist_is_mean_of_neighbours` (replaced by the strictly-between test);
- the mixed-sign rejection tests ×2 (replaced by the b = −1 placeholder);
- `test_forward_unstable_branch_cycle_matches_orbit_order`;
- the private R-tree trim oracle;
- the snap tripwire;
- the private half-edge test.

**Verify:** grepping `build_pieces|_gathered_partitions` in tests/ returns nothing.

### Phase 8: plotting tier (`tests/plotting/`)

`test_plot_smoke.py` checks that every plotter and session delegate runs and returns its type.

`test_plot_topology.py` keeps these properties:
- face point inside its region;
- unbounded node outside the bbox;
- side nodes on their side;
- walks thread face → node → face;
- one cartoon node per element;
- strictly increasing ordinal coordinate;
- brackets = closedness;
- every kept bridge is drawn;
- transition nodes = refined active symbols;
- `class_colors` never repeats;
- the nested circle contains its inner circle;
- the interior side lies inside the circle.

Delete:
- the monkeypatch fan-out/delegate tests;
- the legend wording tests;
- the z-order, palette-order, font and column-width tests;
- the clip-padding reimplementation;
- the label-separator tests;
- the angle and radius pins;
- the gid parsing.

Move `test_name_mathtext` (one canonical case per form) to the naming tests, and `test_verbose_*` to `unit/numerics/test_logging_policy.py`.

### Phase 9: gaps, production wiring, GPU, perf

**`unit/topology/test_production_checks.py`:**
- Spy on `check_holes_share_bridge_side` and `check_bridge_rows_consistent` as called from `Trellis.punch_holes`.
- Add one negative test (`AssertionError` type only).
- Delete the redundant I1/I2 re-check tests.

**Coverage gaps:**
- a synthetic u×u near-tangency is discarded and logged;
- cross-branch anchor-outward orientation, and tangle grouping with connecting classes (hand-built refs, plus the nested law);
- the 2026-09-30 iterated-cut rules: chords only on the base branch, a chain ending on another branch marks nothing, a cross-branch flank is skipped at WARNING;
- a hand-built sub-face rejection in `Arrangement.image_of`;
- real-data virtual symbols are sinks.

**Known issues and the rest:**
- `regression/test_open_p3_deep_runs.py`: I1 still holds, plus a `xfail(strict=True)` test that every class resolves without virtual symbols.
- GPU parity test.
- `test_registry_insert_is_near_linear` → `@pytest.mark.perf`.

### Phase 10: move once, finalize

- **Move-only `git mv` commit** of the survivors into the tiered layout. Basenames stay unique and no `__init__.py` files are added.
- **Second commit:**
  - drop the `slow` and `regression` markers;
  - delete the shims and `minimal_helpers.py`;
  - write the conftest Dev Notes and finish the ledger;
  - update the CLAUDE.md "Run tests" section (`-m golden`, `-m perf`) and the memory notes.

## Verification (end to end)

- `env/bin/python -m pytest -q -rxX`: green. Only `KNOWN_ISSUES` show as xfail, with no XPASS. About 550–650 cases, wall time ≤ ~90 s.
- `pytest -m golden` selects exactly the golden files; `pytest -m perf` selects exactly 1 test; the default run excludes perf.
- `coverage_guard.py` against the Phase 0 baseline: no module loses a line outside the dead-API whitelist.
- The ledger lists every deleted node id against the baseline node-id list (a script diffs them, so none are unaccounted for).
- Isolation spot-checks pass, including running the `invariants/` modules alone.
- A deliberate mutation sanity check: temporarily break one law (e.g. flip the sign in `crossing_sign`) on a scratch branch and confirm that the law tier fails on every case.

## Open items to settle during execution

- P1 outcome: the openings law or an xfail.
- P2 decides the inversion `KNOWN_ISSUES`.
- E8: golden name sign-off.
- E6: questions left over from the audit; until answered, keep current behaviour with assertions checking only the exception type, log level or logger:
  - the uncut-branch pip warning;
  - the dissipative `per_step_beta` test (delete);
  - defensive fallbacks (orientation default, missing parent id, `_near_far`);
  - the verbose print fallback.
