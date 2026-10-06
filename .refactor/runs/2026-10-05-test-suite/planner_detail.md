# Executable plan: refactoring the tanglepack pytest suite

I could not save a plan file because this session has no Write tool and Bash writes are forbidden here. The full plan is below. I created and changed nothing.

## Three findings that shape the plan

1. **Runtime conflicts with decision 10 (every case built per test).**
   - Measured build costs: k10 about 0.05 s, k28 with one or two blasts about 0.35 s, p3 about 0.1 s, nested about 0.5 s. The inversion case takes 1.3 s for growth alone, before the full pipeline.
   - With one test per (law check × case) and a fresh build per test, about 55 law checks × about 3 s of builds comes to about 160 s for the law tier alone.
   - My default is one test per (layer × case) that runs every named check of that layer and reports each check's name and how many items it covered. Checks with a known issue get their own test, so `xfail(strict=True)` stays precise. That gives about 70 law cases in about 35 s.
   - Phase 0 measures this first. This is open item E1.
2. **The hole quirk (decision 21) is probably two different things.** In the snapshots, rows flip between a hole and its backward images, which an orientation-preserving map should not do:
   - k=2.8 two blasts: the iterate-0 'left' hole opens the LEFT pair, but its iterates −1 and −2 open the RIGHT pair.
   - p3: the direct 'right' holes open the LEFT pair, but the backward 'right' holes open the RIGHT pair.
   - The memory note hole-side-own-blast-2026-10-02 lists "latent `_hole_openings` nearest-vertex" as still open.
   - Separately, the single missing anchorward opening matches the documented `_hole_openings` rule: "at the anchor artifact a hole facing the far half … opens nothing".
   - Phase 0 has a dedicated analyst probe with a decision rule.
3. **The catalog recommends deleting both tests of two near-identical pairs.** It marks `test_image_cdist::test_every_bounded_element_now_has_an_image` and `test_partition_elements::test_image_of_element_falls_back_when_an_iterate_is_missing` both "delete". Likewise `test_image_of_element_lands_on_the_advanced_branch` (delete) and `test_image_of_element_covers_both_endpoint_images` (merge). The plan keeps one test of each pair.

---

## A. Phased execution sequence

**Conventions for every phase:**
- `RUN=.refactor/runs/<date>-test-suite`
- Run the suite with `env/bin/python -m pytest -q --durations=15 --cov=tanglepack --cov-report=json:$RUN/cov_pN.json -rxX`.
- Run the coverage check with `env/bin/python tests/_tools/coverage_guard.py $RUN/cov_base.json $RUN/cov_pN.json`.
- Record the count from `pytest --collect-only -q | tail -1`.
- Spot-check isolation by running five random node ids one at a time. Always include `test_dual_graph.py::test_k10_a_different_pip_moves_the_unified_set`, which the memory note footprint-cut-cdist-letters-2026-09-21 records as failing when run alone.
- Each phase is one coder commit. The commit body carries that phase's deletion ledger, which is also appended to `docs/test_suite_refactor_ledger.md` (open item E3).
- **Rule:** a test whose assertions move somewhere else is deleted in the same commit that adds its replacement, never earlier.

### Phase 0: baseline and analyst probes (no test edits)

1. **Baseline.** Write `$RUN/nodeids_base.txt` (from `--collect-only -q`), plus `cov_base.json`, a junit file and durations.
2. **Coverage tool.** Add `tests/_tools/coverage_guard.py`. Pytest does not collect it.
   - For each `src/tanglepack/**` module it compares `missing_lines`.
   - A newly missed line passes only if it falls inside a whitelisted qualname. The spans come from `ast`.
   - Whitelist: `IntersectionRegistry.on_interval`, `IntersectionRegistry._get_lambda_u`, `ManifoldMachine._curvature_area` (the scalar version), `TangleSession.invalidate_trellises`.
   - It prints the per-module percentage change and exits 1 on any other new miss.
3. **Analyst probes.** These are scratch scripts only, and their findings go into ledger §0.
   - **P1, hole quirk (decision 21).** For every hole in `build_period3()` and in the k=2.8 two-blast build, record:
     - `bounding_ids`;
     - whether the near bound is an anchor (stable cdist 0, label "anchor");
     - the origin's direct-hole openings compared with each backward image's openings;
     - for each missing opening, whether its half faces another stable branch or a removed tail.

     Decision rule:
     - **(a) Legitimate:** every difference from "same `(which, row)` pair as the origin's direct hole" is an opening dropped at an anchor or tail. Then the snapshots are replaced by an openings law (Phase 4).
     - **(b) Rows really flip** between the origin and its backward images. Then it is a bug: add `test_propagated_holes_open_their_origins_pair` as `xfail(strict=True)` on the failing cases. The source fix is out of scope.
   - **P2, inversion through the full pipeline.** Run the inversion case stage by stage: intersect, trim, bridges, infer, pips, pseudoneighbors, holes, partition, classes, minimal trellis, iterated partition, dual graph, symbolic dynamics. Record:
     - the first stage that raises, and the exception type;
     - the result of every law check on the stages that do complete;
     - the number of anchors per unstable branch (2 is expected).

     The output is the draft `KNOWN_ISSUES` table.
   - **P3, known open problems.** Run `build_period3(unstable_steps=15)`, `(16)` and `(blasts=6)`. Record which classes are unreachable, which `new*` virtual symbols appear, and whether I1 holds. The output is the exact assertions for the `xfail(strict=True)` tests.
   - **P4, b = −1 placeholder.** Choose (k, b = −1) with a real saddle and a guess (`saddle_guesses` refuses unknown parameters). Confirm that `construct_fixed_point` raises `ValueError` from `set_k_value`.
   - **P5, regression-guard audit (decision 23).** Check every deletion candidate in §B against `memory/codebase-audit-2026-09.md`, `pseudoneighbor-collision-fragility.md`, `pseudoneighbor-partition.md`, `numerics-test-suite.md`, `regions-refactor-2026-09.md`, `topology-layer.md`, `cdist-strict-monotonicity-fix.md`, `hole-side-own-blast-2026-10-02.md` and `audit-polish-2026-07.md`. The output is the Guard column, already filled in §B.
   - **P6, builder parity.** Locally, point the conftest fixtures at the planned `cases.py` builders:
     - k10 at 9 steps with infer, in place of the `henon_tangle_with_bridges` 7+2 recipe;
     - nested with `_repin`, in place of the no-repin `henon_p3_session`;
     - the k=2.8 cases function-scoped.

     Record which node ids fail. This tells Phase 1 which builder parameters must remain.
   - **P7, k=2.8 cut numbers.** Confirm that the numbers in `R_1^1=[0,10]` and the other CLAUDE.md cut facts are registry ids from the build at writing time. Then translate them into the id-free relations given in §D.
   - **P8, runtime model.** Time each law case's full build. This decides E1.
4. **Verify.** No changes to `tests/`, and the node-id list is identical.

### Phase 1: infrastructure, additions only

- **New `tests/cases.py`** (§C). Its builders reproduce today's recipes exactly. Their default parameters are frozen in `cases.py`, and they never read `henon_cases` defaults.
- **New `tests/helpers/` package:**
  - `invariants.py` is the current file moved, with the default changed to `strict=False` (decision 18) and the spike message fixed (print the median).
  - `laws.py` holds the checks. Each returns how many items it checked.
  - `fakes.py` is the current `walk_helpers.py` moved. Phase 7 rewires it.
  - `names.py` holds the letter-free spellers.
- `tests/invariants.py` and `tests/walk_helpers.py` become one-line re-export shims until Phase 10.
- **`pyproject.toml`:** register the `golden` and `perf` markers. Set `addopts = "-ra --strict-markers -m 'not perf'"`. Leave `slow` and `regression` registered until Phase 10.
- **`conftest.py`:** every existing fixture becomes a thin call to a `cases.py` builder, and everything is function-scoped:
  - `k28_partitioned`, `k28_two_blasts_partitioned` and `henon_inversion` lose their session or module scope (decision 10);
  - `henon_p3_session` calls `build_nested(outer_blasts=0, p1_area_cutoff=1e-4, through="zones")`.
- Only builder parameters that P6 showed to matter are kept.
- **Verify:**
  - same node ids;
  - suite green;
  - coverage check passes with no whitelist hits;
  - durations recorded. Expect about +15 s, mostly the 11 inversion users at about 1.3 s each.
- **Deletions:** none.

### Phase 2: pure deletions

Each deletion below has a justification. The Guard column shows the result of the P5 check.

| Test(s) | Count | Reason | Guard |
|---|---|---|---|
| `tests/visualizations/viz_nested_bridges.py`, `viz_unstable_artifacts.py` | – | Decision 24: not collected; crash with KeyError None; match a dead log string | none |
| `numerics/test_workbench_split.py` (whole file) | 20 | Decision 5: layout and signature pins. The patchability guard is redundant because the blast error tests break loudly if patching stops working | none |
| `test_manifold_machine.py` (whole file) | 3 | Decision 17: frozen pre-1.8 grow loop; plus a tautology | none |
| `test_bridge_identity::test_deleted_genealogy_attributes_are_gone` | 1 | Tombstone | none |
| `test_workbench_bugfixes::test_grown_until_intersection_is_gone`, `::test_advance_key_returns_the_next_orbit_point` | 2 | Tombstone; duplicate of the `test_fixed_point` advance_key tests | none |
| `test_map_step_and_graph::test_per_step_beta_matches_the_expression_it_replaces` | 1 | Migration pin | none |
| `test_map_step_and_graph::test_k10_intersection_graph_is_unchanged` | 1 | Decision 17 (counts 8/18/7/7/4) | none |
| `test_map_step_and_graph::test_branch_cycle_reproduces_the_topology_branch_orders` | 1 | Decision 13: forward wrappers | none |
| `test_single_source_of_truth::test_bridge_cutting_pin_k10`, `::…_period_3` | 2 | Decision 17 | none (the refactor is finished) |
| `test_single_source_of_truth::test_tangle_keeps_only_index_state` | 1 | Tombstone | none |
| `test_cleanup_walkers_and_examples::test_string_dispatched_iter_method_is_gone`, `::test_stability_alias_is_defined_once_in_numerics`, `::test_no_test_or_script_defines_its_own_henon` | 3 | Decision 5 | none |
| `test_tangle_index_and_orientation::test_small_tangle_crossing_count_is_unchanged` | 1 | Fixture count | none |
| `test_fixed_point::test_fixed_point_takes_no_branch_count` | 1 | Signature tombstone; `num_branches` is derived and tested | none |
| `test_intersection_registry_fixes::test_get_lambda_u_*` and `::test_on_interval_stable_uses_the_stable_side_eigenvalue` | 5 | Decision 13 (whitelisted) | none |
| `test_closed_form_curvature::test_curvature_area_batch_matches_scalar` | 1 | Decision 13 | none |
| `test_session_trellis_cache::test_invalidate_trellises_is_a_deprecated_no_op` | 1 | Decision 13 | none |
| `test_resonance_zone_region::test_boundary_arc_dataclass_is_gone`, `::test_shapely_is_not_a_dependency` | 2 | Decision 5 | none |
| `test_resonance_zone_region::test_k10_zone_area_and_containment_are_unchanged`, `::test_p3_zone_areas_…` | 2 | Decision 17 (1e-11 pins). Area code stays covered by the winding test | none |
| `test_manifold_initializer::test_initialization_unstable`, `::…_stable` | 2 | Decision 17 ("exactly 3 init points") | none |
| `test_machine_iterate_invariants::test_cdists_are_positive_and_increasing` | 1 | Implied by monotonicity | none |
| `test_minimal_trellis::test_image_chain_is_the_bridge_class_function`, `test_partition_family::test_kinds_are_distinct`, `test_element_naming::test_names_are_deterministic_across_builds` | 3 | Tautologies | none |
| `test_topology_plotting::test_every_plotter_is_a_module_function_in_plotting`, `::test_stable_partition_module_holds_no_drawing_code`, `::test_hole_style_conventions_live_once`, `::test_session_exposes_one_fanout_helper_per_shape` | 4 | Decision 5 | none |
| `test_topology_invariants::test_henon_partition_still_built_after_wiring` | 1 | Historical smoke | none |
| `test_image_cdist::test_image_cdist_is_additive_in_n_on_the_scaling_branch`, `::test_scaled_element_image_stays_on_the_advanced_branch`, `::test_every_bounded_element_now_has_an_image` | 3 | Private / tautology / duplicate. **Keep** `test_partition_elements::test_image_of_element_falls_back_when_an_iterate_is_missing` | none |
| `test_partition_elements::test_image_of_element_lands_on_the_advanced_branch` | 1 | Subsumed. **Keep** `test_image_cdist::test_image_of_element_covers_both_endpoint_images` | none |
| `test_session_bridge_classes::test_session_result_matches_gathered_partitions_helper`, `::test_p3_no_class_mixes_fixed_points` | 2 | Private tautology; verbatim twin. **Keep** `test_bridge_class::test_p3_classes_cover_both_tangles_and_mix_neither` | none |
| `test_session_bridge_classes::test_describe_reports_the_image_evidence` | 1 | String pin; the evidence is asserted structurally | none |
| `test_stable_partition::test_henon_reference_holes_hug_the_stable_manifold`, `::test_stable_arc_midpoint_walks_…`, `::test_stable_arc_midpoint_falls_back_…` | 3 | Private display placement. The coordinates are still exercised by the kept cross-check `test_direct_hole_side_is_the_side_of_its_coordinates` | none |
| `test_stable_partition::test_forward_pairs_beyond_the_fundamental_segment_get_no_hole` | 10 params → firm only | Decision 4: drop the iterate 1..k−1 parameters; keep those with iterate ≥ k_value | holes-backward-only: the firm rule stays |
| `test_dual_graph::test_k10_summary_reports_counts`, `test_dual_walk::test_walk_search_repr_and_walk_repr`, `test_element_naming::test_k10_describe_lists_every_parent_with_its_children` | 3 | String pins. **Replace in the same commit** with one smoke test asserting `summary()`, `repr()` and `describe()` are non-empty, so coverage holds | none |

About 75 cases are removed. **Verify:** suite green, and the coverage check passes with only whitelisted misses.

### Phase 3: loosen in place (two commits)

**Phase 3a** covers `numerics/` and `regression/`. **Phase 3b** covers topology, loom and session. No test functions are deleted; only individual assertions change.

- **Exceptions (decision 8): type only.** Drop every `match=`, including:
  - `test_collect_rejects_a_broken_iterate_link` (no `'NoneType'` match), `test_split_origin_sides_raise`, `test_bridge_rows_inconsistent_raise`;
  - `test_every_bridge_must_be_consecutive_on_its_branch`, `test_trellis_itinerary_errors`, `test_bad_parent_raises_value_error`;
  - `test_saddle_guesses_refuses_unknown_parameters`, `test_grow_until_raises_at_the_cap`, `test_a_missing_partition_names_the_branch`, `test_unstamped_interval_refuses_to_name_itself`;
  - `test_advance_key_rejects_an_out_of_range_branch`, `test_non_convergence_raises_with_the_minpack_message`.
- **Logs (decision 8): level and logger only.** Add a helper `helpers.logs.assert_logged(caplog, level, logger, count=None)`. It applies to:
  - about 13 `test_dual_walk` caplog checks;
  - `test_rebuild_drops_metadata_after_an_unpreserved_recompute`, `test_grow_until_clears_the_iterated_flag_it_invalidates`;
  - `test_partition_warns_without_pseudoneighbors` (drop the exact count), `test_unresolved_loop_stands_alone_and_warns`;
  - `test_no_partitions_raises_and_warns`, `test_a_class_straddling_zones_warns_and_gets_no_zone`;
  - `test_virtual_class_named_new1_with_warning`, `test_cross_side_pair_flagged_with_warning`, `test_inert_class_with_word_warns_and_stays_out`;
  - `test_no_pips_warns_and_unifies_nothing`, `test_a_row_invariant_violation_is_warned_and_skipped`;
  - `test_blast_zone_skips_value_error_with_warning`, `test_a_pair_with_no_bridge_object_is_skipped_with_a_warning`.
- **`describe()`, `summary()` and repr (decision 8): non-empty, and mentioning each class or fixed point.** Applies to:
  - `test_describe_reports` (stable partition), `test_describe_mentions_parents_and_cuts`, `test_describe_has_one_line_per_element`;
  - `test_table_lookups_symbols_and_report`, `test_describe_reports_letters_and_inertness`, `test_homotopy_only_naming_is_its_own_parent`;
  - `test_branch_codes_always_and_letters_only_with_two_fixed_points` (the describe phrase only), `test_is_reliable_and_describe`, `test_sparse_and_dense_agree_when_every_bridge_is_kept`;
  - `test_element_name_is_hashable_and_compares_by_value` (repr), `test_element_ref_hashes_compares_and_labels` (assert only that labels differ), `test_fixed_point_label_wins_…`.
- **Evidence (decision 2):** remove `is_reliable`, `verified` and `source==` from every non-fact test. That includes the synthetic `test_ambiguous_class_is_unreliable`, which keeps `status == "ambiguous"`. The fact tests keep these assertions until Phase 5 moves them to golden.
- **Cache identity inside fact tests (decision 20):** remove `dyn is session.symbolic_dynamics(...)` (`test_k28_two_blasts_structure`), the identity part of `test_k10_symbolic_dynamics_is_built_over_the_cached_pieces`, and `session._bridge_classes` sizes.
- **Tolerances (decision 19).** Use `pytest.approx` and state each tolerance against its library constant:
  - per_step_beta `rel=0` → `approx` (`test_per_step_beta_unstable_is_…`, `…stable_is_the_reciprocal…`);
  - `test_scaled_end_tolerance_is_relative` → against `SCALING_RTOL`;
  - the unlinked halves of `test_table_linked_*` → against `collision_rtol`, **not loosened** (critic note 7);
  - `test_refinement_essentially_converges_to_cutoff` → against `area_cutoff`, using `_curvature_area_batch`;
  - `test_one_map_step_scales_cdist_by_per_step_beta` → comment tying it to `_MIN_SEED_STEP`;
  - `min_unstable_cdist == 0.0` → `approx(0, abs=registry.cdist_tol)`.
- **Strictness:** switch `strict=True` to `strict=False` in `test_machine_iterate_invariants` and `test_growth_integration`. The only strict test left is `test_low_stretch_growth_keeps_cdist_injective` (decision 18).
- **Private access → public:**
  - `_map_batchable` (test_batched_map ×3); `_intersection_registry` → `intersection_registry`;
  - `_built_generation` (`test_same_generation_is_a_cache_hit`);
  - `_partial_bridges`/`_bridges` → `.bridges`, `bridge(id)`, `id is None`;
  - `test_clear_bridges_releases_segment_ownership` → a fresh `create_bridges` reproduces the same set;
  - `test_iterates_closed_all_is_scoped_to_the_grown_fixed_point` → observe ids through the driver;
  - `test_grow_until_intersection_stops_once_a_crossing_exists` → the registry.
- `test_near_vertical_refinement` and `test_refinement_essentially_converges_to_cutoff` switch from the scalar `_curvature_area` to `_curvature_area_batch` (a kernel, so allowed).
- **Verify:** green, coverage check, no XPASS.

### Phase 4: the physical-law tier (`tests/invariants/`), then delete the twins it subsumes

1. Add `tests/invariants/test_law_*.py` (§C), parametrized over `LAW_CASES` = k10, k28_one_blast, k28_two_blasts, p3, nested, inversion.
2. **Non-vacuity:** every law test asserts that each check covered `> 0` items, unless `(law, case)` appears in `NOT_APPLICABLE`, which gives a `pytest.skip` with a reason. Example: image-bridge laws on k28_one_blast and p3, which have 0 image bridges.
3. **Known issues:** every `(law, case)` that P2 or P3 found failing becomes `pytest.param(..., marks=xfail(strict=True, reason=…))` from `cases.KNOWN_ISSUES`. Always included:
   - `one_anchor_per_unstable_branch` on inversion (decision 15);
   - the b = −1 placeholder `test_orientation_reversing_case_runs_the_full_pipeline`, with `xfail(strict=True, raises=ValueError)` (decision 14).
4. Add the openings law from P1, or the xfail if P1 found a bug.
5. **Twins deleted in the same commit** (about 70 cases):
   - **`test_arrangement`:** the k10/p3 Euler, disjoint, sound and region-image tests (8); `test_every_detected_crossing_has_a_definite_sign` (merged into the cross-product law); `test_signs_alternate_along_a_stable_branch` (its own docstring says it is not an invariant).
   - **numerics:** `test_machine_iterate_invariants` (the remaining 7); `test_growth_integration` (2); `test_tangle_intersection_cdist` (3); `test_no_same_stability_crossing` (1).
   - **Bridges and crossings:**
     - `test_bridge_identity::test_bridge_id_is_endpoints_in_unstable_order`, `…ordering_on_period_three`, `test_bridges_at_*` (×2), `test_bridge_identity_on_the_inversion_saddle` (now inversion law parameters plus the anchor law);
     - `test_workbench_bugfixes::test_bridge_endpoints_stay_on_the_bridges_own_branch`. **Guard:** the codebase-audit root cause. The law runs on p3 and nested, so the guard is preserved.
     - `test_single_source_of_truth::test_every_registered_crossing_carries_its_unstable_segment`, `…bridge_endpoints_are_the_ids…`, `…partial_bridge_reports…`;
     - `test_blast_no_overlap::test_workbench_keeps_single_copy_per_bridge`.
   - **Inversion and growth:** `test_inversion_fixture` growth ×4, forward iterates, keys on both branches, `test_the_fixture_really_is_an_inversion_point` (moves to case sanity); `test_invariant_helpers::test_area_preserved_along_every_recorded_iterate_chain`.
   - **Classes:** `test_bridge_class` row_of_end k10/p3, same_branch k10/p3, every bridge in one class, `test_k10_anchor_bridge_runs_source_to_target`, `test_direction_matches_the_element_order_at_the_ends` (the law uses `anchor_outward_key`; critic note 3), `test_oriented_class_orders_anchor_outward`, `test_classes_and_members_come_out_in_the_documented_order` (rewritten as the tangle-grouping law), `test_p3_classes_cover_both_tangles_and_mix_neither`.
   - **Dual graph:** node structure ×2, face side ×2, payload ×2 (tautological, deleted), unified nodes, bipartite.
   - **Partition:**
     - `test_partition_elements`: covers every crossing ×2, singletons, positional element ids;
     - `test_stable_partition`: `test_henon_holes_are_classified`, `…partition_covers_branch`, `…direct_hole_side_reproduces_inward_pair`, `test_p3_propagation_terminates_at_periodicity` (**guard:** higher-period termination, 07-06; kept as a law), `test_k28_blast_child_gets_no_forward_hole` and `test_p3_forward_holes_stop_at_the_branch_return` (firm part only);
     - `test_stable_partition_period3`: `test_p3_propagated_holes_land_on_the_predicted_branch` (using `fp.branch_cycle`), `test_direct_hole_side_is_the_side_of_its_coordinates` (**kept as a law** on all cases, nested added; critic note 1), `test_p3_hole_sides_are_pinned` and `test_k28_two_blast_hole_sides_are_pinned` (per P1);
     - `test_topology_invariants::test_henon_bridge_rows_consistent` (I2 law on all bridges);
     - `test_pseudoneighbor::test_henon_reference_pairs_are_structurally_valid`.
   - **Iterated partition and minimal trellis:** `test_iterated_partition` k10/p3/k28 invariants (the test-side `_empty_stretches` reimplementation dies); `test_minimal_trellis` k10/p3/k28 invariants plus `drops_something`.
   - **Symbolic and naming:** `test_symbolic_dynamics::test_k10_every_landing_contained` and the law half of `test_k10_itineraries_even_and_matrix`; `test_element_naming::test_k10_names_agree_with_the_partition_structure`.
6. **Verify:**
   - green, coverage check, and the `-rxX` listing shows only `KNOWN_ISSUES` entries;
   - the law-tier wall time is at most 40 s;
   - **gate:** if the per-check form is chosen and the law tier exceeds 60 s, fall back to the per-layer form (E1).

### Phase 5: golden tier (`tests/golden/`, `@pytest.mark.golden`), then delete the fact pins

1. Add the three golden files with the §D facts.
2. Every golden module docstring states: "author-gated: change only with author sign-off; no re-record mechanism".
3. `is_reliable` and `verified` assertions now live only here.
4. **Deleted:**
   - `test_bridge_class::test_k10_has_one_active_and_one_inert_class`, `::test_k10_inert_class_rests_on_its_folded_loop_despite_an_unresolved_member`, `::test_k28_has_one_active_and_two_inert_classes`;
   - `test_session_bridge_classes::test_active_class_is_lettered_a_and_the_inert_one_is_not` (the active-only-lettering contract becomes a law), `::test_k28_letters_only_the_anchor_class`;
   - `test_element_naming::test_k10_homotopy_names`, `::test_k10_iterated_names`;
   - `test_iterated_partition::test_k10_empty_stretch_cut_reads_as_expected`;
   - `test_symbolic_dynamics::test_k10_active_class_word`, `::test_k10_refinement_matches_every_member`, `::test_k28_one_blast_singleton_path`;
   - `test_session_symbolic_dynamics::test_k28_two_blasts_structure` and the word substrings in `…built_over_the_cached_pieces` and `test_describe_symbolic_dynamics`;
   - `test_dual_walk::test_k10_landings_and_the_active_class_walk`;
   - `test_higher_period_cartoon::test_p3_words_survive_four_blasts`, `::test_p3_words_are_equivariant_under_the_orbit_shift`, `::test_nested_resolves_everything_with_real_walks`, `::test_nested_default_words_are_pinned` (a commit snapshot, not a CLAUDE.md fact; decisions 2 and 3);
   - the second half of `test_topology_plotting::test_cartoon_brackets_match_closedness`;
   - the k10 facts in `test_plot_transition_graph_draws_every_symbol`.
5. **Verify:** `pytest -m golden` selects exactly the golden files; the full suite is green; coverage check.

### Phase 6: facade and cache tier (`tests/facade/`)

1. Add `test_session_caches.py`. It is table-driven over products (trellis, arrangement, bridge_classes, minimal_trellis, iterated_partition, dual_graph, symbolic_dynamics) × events (hit, rebuild, generation bump, partition-signature change, pip change). The expected hit or miss per cell is a literal table.
   - The generation bump uses a cheap mutation (one `iterate_bridge`), not the 1.4 s regrow.
   - Extra rows: trellis × (growth, compute_intersections on the two-fixed-point nested session, iterate_bridge, rebuild_bridges, add_resonance_zones, restore). **Guard:** blast staleness, 07-06.
   - Fills gap 10: symbolic_dynamics invalidated by a partition-signature change.
2. Add `test_session_equivalence.py`: session result equals the direct build, one case per product on k10 and nested.
3. Add `test_session_fanouts.py`: dict and list fan-out shapes; describe mentions every fixed point; one delegation smoke check.
4. **Deleted:**
   - `test_session_trellis_cache.py` (all remaining tests);
   - the cache, mutation, rebuild and reads tests in `test_session_bridge_classes`, `test_session_dual_graph` and `test_session_symbolic_dynamics`;
   - `test_session_strong_pips::test_trellis_auto_rebuilds_after_registry_change` (now the nested row of the table);
   - `test_arrangement::test_arrangement_is_cached_by_generation`;
   - `test_grow_until::test_session_exposes_the_drivers` (cache half);
   - `test_partition_family::test_from_results_signature_matches_the_session_signature` (private; replaced by the equivalence test);
   - `test_session_dual_graph::test_dual_graph_is_built_over_the_cached_pieces` (private);
   - `test_rebuild_flag_forces_fresh_equivalent_objects` cascade assertions (only the product's own rebuild is asserted).
5. **Verify:** green, coverage check, wall time reduced by about 5 s.

### Phase 7: unit consolidation (7a numerics, 7b topology and loom)

- **Fakes (§C):** fakes delegate to real `FixedPoint` and `Trellis` rules. `FakeFixedPoint` is replaced by `bare_fixed_point(...)`.
- Retire `minimal_helpers.build_pieces` (37 call sites) in favour of the session API (`session.minimal_trellis()`, `iterated_partition()`, `dual_graph()`, `bridge_classes()`). Keep direct assembly only inside `facade/test_session_equivalence.py`.
- Retire the 4 duplicated dual-graph assembly helpers, the local `k10_session` shadows in `test_session_trellis_cache` and `test_loom_blast_restore`, and the 2 copies of `_define_zone`.
- **Merges and parameterizations:**
  - Point and BranchPoint iterate tests (×8 → 4 parametrized); fix the `Point(num_branches)` misuse; rewrite the `test_get_point` tautology;
  - geometry areas (×3 → 1), concave and winding merges;
  - registry generation bumps (×4 → 1); driver caps and rejections across `grow_until`, `grow_until_intersection`, `grow_until_arclength` and `faces_closed`;
  - `test_fixed_points_*` (×5 → 1 table; **guard:** `Intersection.fixed_points` IndexError);
  - synthetic constructor (×2 → 1, keeping the label-slot assertion; **guard**);
  - dual_walk trivial, mirror and landing cases; strong-pip disqualifier cases;
  - blast error handling (×3).
- **Upgrade `test_workbench_has_no_private_key_advance` to a behavioural one-rule guard:** wrap `FixedPoint.advance_key` with a counting spy and assert that `infer_iterate_table`, `iterate_bridge` and `Trellis.image_cdist` all reach it. Keep `test_collect_is_the_only_traversal` (decision 5). Open item E9.
- **Strengthen:**
  - `test_kevin_way_period_three_orbit_chain`: assert the chain order (unstable 0,1,2; stable 0,2,1). **Guard:** the Kevin-way bug.
  - `test_period_three_orbit_at_k_two_is_accepted`: assert f³(p) = p and the eigenvalue product (gap 11).
  - Initializer alpha/k for k > 1 (gap 11).
  - `test_oriented_bridge_polyline_runs_in_the_dynamical_direction`: assert that `poly[0]` is the low-cdist end.
  - `test_trim_stable_manifolds_reads_the_registry`: assert that something was trimmed.
  - `test_pair_on_different_unstable_branches_rejected`: use a clean fixture.
  - `test_trajectory_extension_forward_with_dedup`: build a real duplicate.
  - `test_blast_recognizes_already_known_bridges`: assert known bridges are not re-iterated (gap 6).
- **Deleted:**
  - `test_workbench_bugfixes::test_trim_stable_manifolds_cuts_just_past_the_outermost_crossing` (private R-tree oracle);
  - `test_pseudoneighbor::test_forward_unstable_branch_cycle_matches_orbit_order` (decision 13 wrapper);
  - `test_inversion_fixture::test_mixed_eigenvalue_signs_are_rejected` and `test_fixed_point::test_set_k_value_rejects_disagreeing_eigenvalue_signs` (the b = −1 xfail placeholder now encodes this; decision 14);
  - `test_grow_until::test_faces_closed_raises_at_the_cap_for_the_anchor` (decision 4; Dev Note; the cap is covered by the parametrized driver cap);
  - `test_refinement::test_refined_cdist_is_mean_of_neighbours` (decision 4; Dev Note), replaced in the same commit by `test_refined_cdist_lies_strictly_between_its_neighbours`;
  - `test_symbolic_dynamics::test_k10_itineraries_even_and_matrix` remainder;
  - the `_region_key` test (rewritten: both holes of a singleton bridge survive propagation);
  - `test_sparse_kept_endpoints_are_degree_three_without_unstable_stubs` (private half-edges; rewrite only if a public degree query exists, otherwise delete);
  - `test_image_cdist::test_the_snap_does_not_fire_on_the_default_path` (private fixture tripwire; E6).
- **Kernel exceptions kept on private APIs (decision 6):**
  - `_side_of`, `_resolve_inertness`, `_backward_endpoint`, `_bridge_unstable_span` and `_containing_bridge` (**guard:** root cause), `_near_far`;
  - `_collapse_noise_crossings` (**guard**), `_do_segments_intersect` and the orientation predicates (**guard:** `_orientation` eps);
  - `_parabolic_fit`, `_curvature_area_batch`, `PartitionFamily` marks round-trip;
  - the iterated-cut chord kernel (`test_chords_pair_only_the_base_branchs_crossings`, `test_a_lobe_ending_on_another_branch_marks_nothing`), `_inert_letter` → public `inert_letters`.
- **Verify:** green, coverage check, `grep -r "build_pieces\|_gathered_partitions" tests/` is empty.

### Phase 8: plotting tier (`tests/plotting/`)

- **Keep:**
  - smoke: `test_plotters_draw_on_a_real_tangle`, the `test_session_pseudoneighbors` and `test_session_strong_pips` plot helpers, `test_plot_stable_partition_smoke` (no tick count), `test_henon_plot_helpers_smoke`, the session plot delegates (return types only);
  - topology: face point inside its region; unbounded node outside the bbox; face-point copy (regression); side nodes on their side (through public `stable_node_point`); walk threads face → node → face; cartoon has one node per element per side; ordinal coordinate strictly increasing; brackets match closedness (first half); every kept bridge drawn and walks go element to element; transition-graph node set equals the active refined symbols; `plot_bridges_by_class` draws every classed bridge; `plot_minimal_trellis` draws every kept bridge; `class_colors` never repeats; scatter kwargs reject facecolors and `c` (type only); nested inner circle inside outer (containment only); circle interior side inside the circle.
- **Deleted:**
  - `test_trellis_plotters_delegate_to_plotting`, `test_session_plot_fanouts_go_through_fanout_plot` ×5 and `test_session_call_fanouts_go_through_fanout_call` ×4. These are replaced in the same commit by a behavioural `session.plot_*` smoke test per fixed point, so session coverage holds;
  - `test_session_dual_graph::test_session_plot_delegates_to_plotting`, `test_session_symbolic_dynamics::test_session_plot_delegates_forward_to_plotting`;
  - legend wording (`test_dual_graph_legend_handles_match_the_plotter`, `test_cartoon_legend_handles…`, `test_walk_legend…`), `test_walk_zorder_sits_between…`, palette order (`test_class_colors_are_assigned_in_fixed_order` keeps only the no-repeat half);
  - `test_plot_itinerary_table_columns_font_and_rows`, the width and font parts of `test_plot_itinerary_table_lists_every_class` (keeps "one row per class");
  - `test_clip_to_arcs_*` (padding reimplementation), `test_show_labels_*` separators, `test_stable_node_label_falls_back…`, `test_k10_stable_node_label_joins_the_names`, `test_unmatched_member_gets_its_parent_letter_slot`;
  - `test_circle_layout_puts_every_node…` angle, sweep and radius pins; `test_line_cartoon_labels…`; `test_nested_outer_arcs_go_around_the_inner_circle` (gid parsing; rewrite as a keepout property only if `ZoneLayout` exposes routed paths publicly); `test_cartoon_bridge_colour_labels_and_walkless_legend` alignment pins.
- `test_name_mathtext` ×11 moves to `unit/topology/test_element_naming.py`, reduced to one canonical case per notation form. The `test_verbose_*` tests ×5 move to `unit/numerics/test_logging_policy.py` (level and logger only; E6).
- **Verify:** green under the Agg backend, coverage check.

### Phase 9: gap tests, production-check wiring, GPU, perf

- **Production-check wiring (decision 22)** in `unit/topology/test_production_checks.py`:
  - Spy on `tanglepack.topology.StablePartition.check_holes_share_bridge_side` and `check_bridge_rows_consistent`. Both are imported locally inside `Trellis.punch_holes`, so a module-attribute patch works. Assert that one `punch_holes` call checks every hole once with the trellis orientation, and every hole-bearing bridge once.
  - One negative test: a corrupted side raises `AssertionError` (type only).
  - Also assert that `DualGraph` construction runs its self-check, and that the Arrangement's consecutive-bridge guard fires (existing test).
  - **Delete** `test_topology_invariants::test_henon_holes_share_bridge_side` and `test_stable_partition_period3::test_p3_holes_share_bridge_side` and `::test_p3_bridge_rows_consistent`. The I2 law on all bridges already covers what remains.
- **Gaps from synthesis §5:**
  - a synthetic u×u near-tangency is discarded and logged at DEBUG by the right logger (unit);
  - `preserve_ids` keeps ids and cdists and `compute_intersections` is idempotent (law, mutating);
  - area preserved on all cases (law, already in Phase 4);
  - cross-branch anchor-outward orientation with hand-built `ElementRef`s on two branches, and tangle grouping with connecting classes (`tangle=None` last), as unit tests plus a law on nested;
  - the 09-30 iterated-cut rules: chords only on the base branch, a chain ending on another branch marks nothing, a cross-branch flank is skipped at WARNING (`unit/topology/test_iterated_cut.py`);
  - `min_separation`: every kept child is at least `min_separation` from the accumulated interiors (memory numerics-test-suite);
  - refined cdist strictly between its neighbours;
  - a hand-built sub-face rejection in `Arrangement.image_of` (gap 12);
  - real-data symbolic paths where they exist: virtual symbols are sinks on `build_period3(unstable_steps=15)` (passes even though that run is unreliable); ambiguous / unmatched / evidence-length mismatch only where P3 finds them, otherwise synthetic.
- **Known open problems (decision 16)** in `regression/test_open_p3_deep_runs.py`. Parameters are 15 steps, 16 steps and 6 blasts:
  - `test_deep_p3_holes_share_bridge_side` passes; it moves here from `test_deep_p3_holes_share_bridge_side_through_symbolic_dynamics`;
  - `test_deep_p3_every_class_resolves_without_virtual_symbols` is `xfail(strict=True)` with P3's reason.
- **GPU:** `unit/numerics/test_gpu.py`.
  - `test_gpu_growth_matches_cpu` uses `pytest.importorskip("cupy")` and skips when no device is available.
  - `test_enable_gpu_without_cupy_raises_clearly` monkeypatches `sys.modules["cupy"] = None`, so it actually runs.
- **Perf:** `test_registry_insert_is_near_linear` moves to `unit/numerics/test_registry_perf.py` with `@pytest.mark.perf`.
- **Verify:** `pytest -m perf` runs 1 test; the default run deselects it; green; coverage check.

### Phase 10: move files once and finalize

- `git mv` only. Survivors move into the layout in a **move-only commit** so git detects the renames. A separate commit then:
  - removes the `slow` and `regression` markers (decision 12) and their registrations;
  - deletes the shims `tests/invariants.py`, `walk_helpers.py` and `minimal_helpers.py`;
  - fixes the stale conftest and k28 docstrings;
  - adds a "Dev Notes" section to the `tests/conftest.py` docstring (E3) covering: untested provisional rules (the +1..+(k−1) exemption, the mean-of-neighbours refined cdist, the anchor faces-closed limitation); dead API kept without tests (`on_interval`/`_get_lambda_u`, `invalidate_trellises`, `BoundaryArc`, scalar `_curvature_area`, `forward_*_branch_cycle` wrappers); and the index of `KNOWN_ISSUES`.
- Every test basename must be unique, or each directory gets an `__init__.py` (E10).
- **Verify:** final node-id count (expect about 550–650); coverage check against the Phase 0 baseline; total wall time at most 90 s; `-m golden` and `-m perf` selection works.

---

## B. Disposition of all 64 current files

Approximate current case counts are in brackets. Abbreviations: **L/** = `invariants/test_law_*`, **U/n**, **U/t**, **U/l** = `unit/numerics`, `unit/topology`, `unit/loom`, **F/** = `facade/`, **G/** = `golden/`, **R/** = `regression/`, **P/** = `plotting/`.

### numerics/ (27 files)

| File | Destination | Action | Notable tests |
|---|---|---|---|
| test_batched_map [4] | U/n/test_growth.py | keep+loosen | Drop the `_map_batchable` asserts; refinement test stated against `area_cutoff` and uses the batch curvature |
| test_blast_no_overlap [3] | L/bridges + U/l/test_blast.py | split | single_copy → law; `test_iterating_fixed_point_bridge_returns_existing_copies` **guard** (single-copy, numerics-test-suite); already_known rewritten (gap 6) |
| test_blast_proximity_guard [2] | U/l/test_blast.py | rewrite+merge | `test_min_separation_drops_close_bridges` asserts distances (**guard:** zig-zag) |
| test_bridge_identity [17] | L/bridges, U/n/test_bridges.py, F/fanouts | split | Delete `test_deleted_genealogy_attributes_are_gone`; id-order and bridges_at → law; inversion-saddle test → law params plus anchor xfail |
| test_cleanup_walkers_and_examples [18] | U/n/test_linked_list.py, U/n/test_solver_and_initializer.py | split | **Keep** `test_collect_is_the_only_traversal` (one rule); delete 3 lints; broken-iterate test is type-only |
| test_closed_form_curvature [6] | U/n/test_growth.py | keep 5 / delete 1 | Delete `test_curvature_area_batch_matches_scalar` (dead API) |
| test_generation_and_caches [24] | U/n/test_generation.py, U/n/test_registry.py, U/n/test_registry_perf.py | split+parametrize | Registry bumps → 1; near-linear → perf; memo `is` asserts dropped |
| test_geometry [27] | U/n/test_geometry.py | parametrize+rewrite | Oriented polyline asserts orientation |
| test_gpu [2] | U/n/test_gpu.py | keep+rewrite | Parity with importorskip; no-CuPy via monkeypatch |
| test_grow_until [19] | U/n/test_growth_drivers.py, F/fanouts | parametrize | Delete `test_faces_closed_raises_at_the_cap_for_the_anchor` (decision 4); **guard:** `test_grow_until_recuts_existing_bridges` |
| test_growth_integration [2] | L/manifolds + U/n/test_growth.py | move | Per-step non-strict check in unit |
| test_initializer_cdist [6] | U/n/test_solver_and_initializer.py | keep+merge | Alpha tautology dropped; k > 1 added |
| test_intersection_registry_fixes [13] | U/n/test_registry.py | parametrize / delete 5 | **Guards:** `test_fixed_points_*` (IndexError), `test_synthetic_does_not_put_the_label_in_the_id_slot`; delete `_get_lambda_u` ×4 and `on_interval` |
| test_invariant_helpers [4] | U/n/test_law_helpers.py + L/crossings | split | **Keep** `test_low_stretch_growth_keeps_cdist_injective` (the only strict test) with the fundamental-segment case merged in; `…rejects_a_broken_chain` is a meta-test |
| test_inversion_fixture [14] | L/case_sanity, L/manifolds, L/crossings, U/n/test_fixed_point.py | split | Mixed-sign rejection becomes the b = −1 xfail placeholder |
| test_machine_iterate_invariants [8] | L/manifolds | merge/delete | Delete `test_cdists_are_positive_and_increasing` |
| test_map_step_and_graph [26] | U/n/test_fixed_point.py, U/n/test_registry.py, U/t/test_trellis.py | split | **Keep** `test_workbench_has_no_private_key_advance` (upgraded); `k_th_root_only_without_inversion` merged as an explicit inversion case (**guard:** k_value-root bug); delete 3 (Phase 2) |
| test_minimal_solver [6] | U/n/test_solver_and_initializer.py | keep | Strengthen period-3 (**guard:** solver root cause) |
| test_no_same_stability_crossing [1] | L/crossings + U/n/test_crossing_kernel.py | rewrite | u×u and s×s across branches and fixed points |
| test_noise_crossing_collapse [3] | U/n/test_crossing_kernel.py | keep | **Guard:** noise collapse, 07-07 |
| test_refinement [1] | U/n/test_growth.py | replace | Strictly-between-neighbours test; mean rule → Dev Note |
| test_segment_index_bookkeeping [4] | U/n/test_linked_list.py, U/n/test_crossing_kernel.py | split | **Guards:** stretch_param (audit-polish), `_edge_seen` |
| test_single_source_of_truth [12] | L/crossings, L/bridges, U/n/test_bridges.py, U/n/test_registry.py | split / delete 3 | Delete `bridge_cutting_pin` ×2 and `tangle_keeps_only_index_state` |
| test_tangle_index_and_orientation [6] | U/n/test_crossing_kernel.py | keep / delete 1 | **Guards:** scale-invariant predicates, own segments |
| test_tangle_intersection_cdist [3] | L/crossings | move | – |
| test_workbench_bugfixes [10] | L/bridges, U/n/test_growth_drivers.py | split / delete 3 | **Guard** (p3 law): `test_bridge_endpoints_stay_on_the_bridges_own_branch`; **guard:** `…honours_branch_index` |
| test_workbench_split [20] | – | delete | – |

### regression/ (5 files)

| File | Destination | Action | Notable tests |
|---|---|---|---|
| test_blast_monotonicity | R/test_blast_monotonicity.py | keep+strengthen | Assert manifolds monotone after the blast |
| test_boundary_straddle_cdist | R/ | keep | **Guard:** straddle cdist |
| test_cdist_collision_growth | R/test_high_stretch_growth.py | merge | **Guard:** cdist-strict memory |
| test_high_stretch_period3_growth | R/test_high_stretch_growth.py | keep | Add a non-vacuity check that the run really reaches near-ULP spacing |
| test_near_vertical_refinement | R/ | keep | Switch to `_curvature_area_batch` |

### Top-level (32 files)

| File | Destination | Action | Notable tests |
|---|---|---|---|
| test_arrangement [28] | L/arrangement, L/crossings, U/t/test_arrangement.py, F/caches | split | Delete `test_signs_alternate_along_a_stable_branch`; sub-face test rewritten hand-built; `…survives_a_blast` kept (**guard**) |
| test_arrangement_sparse [6] | U/t/test_arrangement.py | loosen | Summary strings → counts |
| test_branch_point [9] | U/n/test_linked_list.py | parametrize | – |
| test_bridge_class [30] | L/classes, U/t/test_bridge_class.py, G/ | split | **Guard:** `test_p3_period_three_classes_use_every_stable_branch` (law); `_resolve_inertness` is a kernel |
| test_dual_graph [27] | L/dual_graph, U/t/test_dual_graph.py, U/t/test_trellis.py | split | Delete summary and expected-names tests; keep pip-branch-only unification |
| test_dual_walk [34] | U/t/test_dual_walk.py | loosen logs | Delete the repr test; k10 walk → G/ |
| test_element_naming [19] | U/t/test_element_naming.py, G/, L/iterated | split | Delete deterministic and describe tests; `missing_parent_id` → E6 |
| test_fixed_point [16] | U/n/test_fixed_point.py | keep | Delete `takes_no_branch_count` and `set_k_value_rejects_…` |
| test_higher_period_cartoon [19] | U/t/test_iterated_cut.py, U/t/test_element_naming.py, P/, R/test_nested_own_blast.py, G/ | split | **Guards:** `test_nested_inner_words_are_the_period3_words` (letter-free) plus `…blast_order_does_not_matter`, with `…outer_blasts_leave_the_inner_bridges_alone` merged in; delete `test_nested_default_words_are_pinned` |
| test_image_cdist [15] | U/t/test_trellis.py, U/t/test_partition_family.py | split / delete 4 | Keep `covers_both_endpoint_images` |
| test_iterated_partition [9] | L/iterated, U/t/test_iterated_cut.py, G/ | rewrite | The test-side `_empty_stretches` oracle dies |
| test_loom_blast_restore [7] | U/l/test_blast.py | keep | **Guard:** assertions never swallowed |
| test_manifold_initializer [7] | U/n/test_solver_and_initializer.py, L/case_sanity | split | Delete 3-point ×2; kevin-way p3 strengthened (**guard**); two-branch init → inversion sanity (**guard:** has_inversion) |
| test_manifold_machine [3] | – | delete | – |
| test_minimal_trellis [8] | L/iterated, U/t/test_minimal_trellis.py | split | Public API in place of `_check_minimal` internals |
| test_partition_elements [15] | L/partition, U/t/test_partition_family.py, F/fanouts | split | Keep the fallback test (twin of the deleted image_cdist test) |
| test_partition_family [10] | U/t/test_partition_family.py, F/equivalence | keep / delete 1 | – |
| test_point [10] | U/n/test_linked_list.py | parametrize+rewrite | – |
| test_pseudoneighbor [21] | U/t/test_pseudoneighbor.py, L/partition | keep | **Guards:** table-linked ×2 with the `collision_rtol` half, backward-iterate puncture |
| test_resonance_zone_region [13] | U/l/test_resonance_zone.py | delete 4 / loosen | Midpoint classification via the public `classify_bridge` |
| test_session_bridge_classes [19] | F/caches, F/equivalence, U/t/test_bridge_class.py, G/ | split / delete 3 | – |
| test_session_dual_graph [11] | F/, P/ | split | – |
| test_session_pseudoneighbors [7] | F/fanouts, P/ | merge | – |
| test_session_strong_pips [4] | F/fanouts, F/caches, P/ | merge | – |
| test_session_symbolic_dynamics [11] | F/, G/test_golden_k28.py, P/ | split | Delete `p3_symbolic_dynamics_smoke` and the forward_to_plotting mock |
| test_session_trellis_cache [10] | F/caches | merge / delete 1 | **Guard:** restore invalidates |
| test_stable_partition [45] | U/t/test_stable_partition.py, U/t/test_trellis.py, L/partition, R/test_hole_side_cross_checks.py, P/ | split | **Keep** `test_henon_propagated_hole_side_matches_coords` (critic note 2); **guards:** `bridge_span_own_key`, `containing_bridge`, `repunching`, `partition_ignores_strong_pip_cut`, `backward_endpoint`; provisional parameters deleted |
| test_stable_partition_period3 [11] | L/partition, R/test_open_p3_deep_runs.py, wiring | split | **Keep** `test_direct_hole_side_is_the_side_of_its_coordinates` as a law; hole-side snapshots per P1 |
| test_strong_pip_periodic_point [4] | U/t/test_strong_pip.py | parametrize | **Guards:** self-disqualification, collision |
| test_symbolic_dynamics [22] | U/t/test_symbolic_dynamics.py, L/symbolic, G/ | split | `is_reliable` assertions removed (decision 2) |
| test_topology_invariants [12] | U/t/test_topology_invariants.py, L/partition | split / delete 2 | Keep the orientation-reversing parity (synthetic) |
| test_topology_plotting [65] | P/test_plot_smoke.py, P/test_plot_topology.py, U/t/test_element_naming.py, U/n/test_logging_policy.py | split / delete about 30 | – |

### Non-test files

| File | Disposition |
|---|---|
| conftest.py | Rewritten (§C) |
| invariants.py | → `helpers/invariants.py` |
| walk_helpers.py | → `helpers/fakes.py` |
| minimal_helpers.py | Deleted (Phase 10) |
| visualizations/ | Deleted |

---

## C. New layout, builders, helpers, conftest

| File | Content | Estimated cases |
|---|---|---|
| invariants/test_law_case_sanity.py | Each law case builds through the full pipeline; inversion is real (k_value = 2·period); seeded inversion builds both branches; b = −1 placeholder xfail | 8 |
| invariants/test_law_manifolds.py | Non-strict cdist monotonicity, no spikes, iterate law, one-to-one on every manifold and bridge | 6 (layer×case) or 24 (per check) |
| invariants/test_law_crossings.py | Every crossing is one unstable a-key and one stable b-key; no u×u or s×s; cdists bracketed; crossing sign equals the cross product; area preserved along chains; preserve_ids and idempotent recompute | 6 + 6 |
| invariants/test_law_anchors.py | One anchor per unstable branch (xfail on inversion) | 6 |
| invariants/test_law_bridges.py | BridgeId in unstable order; consecutive crossings on its own branch; bridges_at exact; single copy, no overlap; image/preimage round trip; partial ⇔ id None | 6 / 36 |
| invariants/test_law_partition.py | Covers each branch per side; unique owner; singletons; direct-hole side equals coordinate side; direct holes open the inward pair; propagated holes land on the predicted branch and stop at periodicity; no direct hole at iterate ≥ k_value; reference pairs valid; openings law (P1) | 6 / 54 |
| invariants/test_law_arrangement.py | Euler per component; disjoint; sound; `image_of` agrees with dynamics and area (non-vacuous); preimage∘image | 6 / 30 |
| invariants/test_law_classes.py | Every bridge in one class; anchor-outward orientation; row-of-end equals geometry; I2 on all bridges; tangle grouping, then min-cdist order, connecting classes last; no fixed-point mixing unless connecting; only active classes lettered | 6 / 42 |
| invariants/test_law_iterated.py | Child inside parent; homotopy boundaries kept; unique owner; every crossing in a closed element; names agree with structure; minimal trellis public invariants | 6 / 30 |
| invariants/test_law_dual_graph.py | Wall and unified structure; one face per region with a single outer; face side equals geometry; unified nodes equal the pip segment on the pip's branch; bipartite degree rule | 6 / 30 |
| invariants/test_law_symbolic.py | Even itineraries of same-side pairs; each pair is a class or its inverse; landings contained; refined children inherit the word; member matching consistent; inert classes outside the graph and matrix; matrix equals the word token counts | 6 / 42 |
| unit/numerics/ | test_geometry (22), test_linked_list (20), test_fixed_point (28), test_solver_and_initializer (22), test_growth (16), test_growth_drivers (22), test_crossing_kernel (14), test_registry (22), test_registry_perf (1, perf), test_bridges (18), test_generation (14), test_law_helpers (2), test_gpu (2), test_logging_policy (4) | ≈ 207 |
| unit/topology/ | test_pseudoneighbor (18), test_strong_pip (4), test_stable_partition (30), test_topology_invariants (9), test_production_checks (4), test_arrangement (16), test_bridge_class (22), test_partition_family (22), test_iterated_cut (8), test_minimal_trellis (4), test_element_naming (18), test_dual_graph (10), test_dual_walk (26), test_symbolic_dynamics (16), test_trellis (12) | ≈ 220 |
| unit/loom/ | test_blast.py (10), test_resonance_zone.py (8) | 18 |
| facade/ | test_session_caches.py (≈ 32), test_session_equivalence.py (≈ 10), test_session_fanouts.py (≈ 10) | ≈ 52 |
| golden/ | test_golden_k10.py (≈ 6), test_golden_k28.py (≈ 7), test_golden_period3.py (≈ 4) | ≈ 17 |
| regression/ | test_blast_monotonicity, test_boundary_straddle_cdist, test_high_stretch_growth (2), test_near_vertical_refinement, test_nested_own_blast (2), test_hole_side_cross_checks (1), test_open_p3_deep_runs (6) | ≈ 14 |
| plotting/ | test_plot_smoke.py (≈ 12), test_plot_topology.py (≈ 16) | ≈ 28 |

**Totals:** about 620 cases with the per-layer law tier, or about 860 with the per-check law tier.

### `tests/cases.py` (tests only; parameters frozen here, never read from `henon_cases` defaults)

```python
Stage = Literal["seeded", "grown", "intersected", "bridges", "pips", "zones", "partitioned"]

@dataclass
class Case:
    name: str
    session: TangleSession
    fixed_points: list[FixedPoint]          # outermost first
    pips: list[Optional[int]]
    expect: CaseExpect                       # period(s), inversion, n_tangles, has_image_bridges, blasts
    @property
    def workbench(self) -> TangleWorkbench: ...
    @property
    def registry(self) -> IntersectionRegistry: ...

def build_k10(*, unstable_steps: int = 9, through: Stage = "partitioned") -> Case
def build_k28(*, blasts: Literal[1, 2], through: Stage = "partitioned") -> Case
    # 10 unstable steps, area_cutoff 1e-7, pip = f(q0), zone trim, min_separation 1e-5
def build_period3(*, unstable_steps=13, stable_steps=9, area_cutoff=1e-7, blasts=0,
                  min_separation=1e-5) -> Case
def build_nested(*, outer_blasts=2, inner_blasts=0, p1_area_cutoff=1e-7, p3_unstable_steps=13,
                 p3_stable_steps=9, p1_unstable_steps=11, min_separation=1e-4,
                 through: Stage = "partitioned") -> Case
def build_inversion(*, unstable_steps=6, stable_steps=5, through: Stage = "partitioned") -> Case
def build_orientation_reversing(*, k: float = <P4>, through: Stage = "partitioned") -> Case
    # raises ValueError today

LAW_CASES = ("k10", "k28_one_blast", "k28_two_blasts", "p3", "nested", "inversion")
BUILDERS: Mapping[str, Callable[[], Case]]          # MappingProxyType; the frozen defaults
KNOWN_ISSUES: Mapping[tuple[str, str], str]         # (law_id, case) -> reason, from P2/P3
NOT_APPLICABLE: Mapping[tuple[str, str], str]       # (law_id, case) -> reason (skip)

def law_params(law_id: str) -> list[ParameterSet]   # applies xfail(strict=True) and skip marks
```

Whether `build_period3` and `build_nested` delegate to `henon_cases` with every keyword passed explicitly, or reimplement the recipe, is open item E4.

### Helpers

- **`helpers/invariants.py`:** the existing manifold checks, with `strict=False` as the default and the message fixed. `assert_no_cdist_collision` is documented as "sanctioned low-stretch only".
- **`helpers/laws.py`:**
  - each check is `def check_<law>(case: Case) -> int` and returns the number of items it checked;
  - `LAYERS: dict[str, list[Callable]]`;
  - `run_layer(case, layer) -> dict[str, int]` runs every check, collects the failures and raises one `AssertionError` that names each failing check.
- **`helpers/fakes.py`:**
  - `bare_fixed_point(period, *, inversion=False, beta=0.25, label=None) -> FixedPoint` returns a **real** `FixedPoint`: it sets the eigenvalues (λu = (1/β)^period, sign negative when `inversion=True`), calls `set_k_value()`, and stamps the label. This replaces `FakeFixedPoint`, so `advance_key`, `per_step_beta` and `num_branches` are the production rules.
  - `FakeTrellis` borrows `image_cdist = Trellis.image_cdist` and `_scaled_image_cdist = Trellis._scaled_image_cdist`. It supplies only data (`iterate`, `_keyed_cdist`, `intersection`, `registry`), so the "table first, else `advance_key` and `per_step_beta`" rule is no longer reimplemented.
  - `make_result` requires explicit `owners`; the derived-ownership default is removed.
  - `FakeDual` stays as is (risk R7).
- **`helpers/names.py`:**
  - `homotopy_pair(dyn, cd) -> tuple[str, str]` gives the oriented source and target in short homotopy names;
  - `word_in_names(dyn, cd) -> list[tuple[str, str]]` gives each symbol as its oriented element pair, with direction applied;
  - `brackets(result) -> list[tuple[str, str]]` gives (name, `"[ ]" | "( )" | "[ )" | "( ]"`);
  - `matrix_by_classes(dyn, refined) -> dict[frozenset, dict[frozenset, int]]`.
- **`helpers/logs.py`:** `assert_logged(caplog, level, logger, *, count=None)`.

### `tests/conftest.py`

- Set `matplotlib.use("Agg")` once.
- Every fixture is function-scoped and calls a builder: `k10`, `k10_partitioned`, `k28_one_blast`, `k28_two_blasts`, `p3`, `nested`, `nested_unblasted`, `inversion`, `inversion_seeded`, and the numerics stage fixtures (`k10_fixed_point`, `k10_initialized`, `k10_grown`, `k10_bridges`, all `build_k10(through=...)`).
- `law_case` is an indirect fixture: `request.param` is a case name, and it returns `BUILDERS[name]()`.
- The Dev Notes docstring (Phase 10).
- No session or module scope anywhere (decision 10).
- `pyproject.toml` registers the markers `golden` (run by default) and `perf` (excluded by default through `addopts -m 'not perf'`).

---

## D. Golden facts per case, spelled without letters

Each fact is pinned once. Classes are written as oriented element pairs (source → target, anchor outward), and words as sequences of oriented pairs (a reversed pair means the inverse). Letters never appear.

### `test_golden_k10.py` (`build_k10()` defaults)

1. **Iterated rows, anchor outward, with closedness.** The name-to-bracket pairing is to be confirmed against today's output in Phase 5 and signed off by the author.
   - Right row: R_1^1 `[ ]`, R_1^2 `( )`, R_1^3 `[ ]`, R_2 `( )`, R_3^1 `[ ]`, R_3^2 `( )`, R_3^3 `[ ]`.
   - Left row: L_1^1 `[ )`, L_1^2 `[ ]`, L_2 `( )`, L_3 `[ ]`.
   - There are 11 iterated elements. R_1 and R_3 each have 3 children, L_1 has 2, and the others are unsplit.
2. **Classes.** Exactly two:
   - active X = (R_1 → R_3);
   - inert U = (L_1 → L_3).
3. **Walk itinerary of X.** R_1^1 R_3^3 | L_3 L_1^2 | R_3^1 R_1^3. A unique walk, not ambiguous.
4. **Image of X:** [(R_1 → R_3), (L_3 → L_1), (R_3 → R_1)], that is X, U⁻¹, X⁻¹.
5. **Refinement.** X splits into exactly two children:
   - X₁ = (R_1^1, R_3^3);
   - X₂ = (R_1^3, R_3^1).

   Both children's words are [X₁, U⁻¹, X₂⁻¹]. Every non-loop member matches a child, and the unmatched set is empty.
6. **Matrix.** The refined transition matrix over (X₁, X₂) is [[1,1],[1,1]]. U is in neither the graph nor the matrix.
7. **Inert evidence.** U is inert, has exactly one member with no registered image, and rests on its folded loop (the folded-loop part is to be confirmed by the author; see E8).
8. **Evidence.** `X.verified is True`, and `dyn.is_reliable` is true.

### `test_golden_k28.py`

**One blast:**
1. Exactly one active class X, singleton to singleton (both landings are singletons), resolved through the trellis path (E2).
2. Word of X = [X, U⁻¹, V⁻¹], where U and V are distinct inert classes and both are trivial loops (empty words).
3. The exterior class is inert through a virtual loop: an image pair that no bridge spans.
4. The minimal trellis has 0 image bridges. This comes from the memory note minimal-dual-graph, not CLAUDE.md; to be confirmed by the author.
5. `X.verified`.

**Two blasts:**
1. **Classes.** Active A = (R_1 → R_5), B = (R_3 → R_5), C = (R_1 → R_3). Inert U = (L_1 → L_3). Exactly 3 active and 1 inert.
2. **Order.** A is first in the table (it holds the anchor bridge; minimum member cdist ≈ 0 within `cdist_tol`). The table is sorted by minimum cdist.
3. **Words:**
   - A → [(R_1→R_5), (L_3→L_1), (R_5→R_3)];
   - B → [(R_1→R_3)];
   - C → [(R_1→R_5), (L_3→L_1), (R_5→R_1)].
4. **Refinement.** Only A refines:
   - A₁ = (R_1^1, R_5^3);
   - A₂ = (R_1^3, R_5^1).

   A₁ and A₂ have equal words. C's refined word is [A₁, U⁻¹, A₂⁻¹]. Every member matches a child.
5. **Matrix.** The unrefined matrix over (A, B, C) is [[1,1,0],[0,0,1],[2,0,0]]. U is absent from it.
6. **Cut closedness:**
   - R_1^1 `[ ]`, R_1^2 `( )`, R_1^3 `[ ]`;
   - R_5^1 `[ ]`, R_5^2 `( )`, R_5^3 `[ ]`;
   - L_1^1 `[ )`, L_1^2 `[ ]`.

   Id-free relations from P7:
   - R_1^1 begins at the anchor;
   - R_1^1.hi = R_1^2.lo and R_1^2.hi = R_1^3.lo, so R_1^2 is the open base;
   - R_5^1.hi = R_5^2.lo and R_5^2.hi = R_5^3.lo, so R_5^2 is the open chord;
   - L_1^2 has the same two ends as R_5^1;
   - the L_1^1 / L_1^2 boundary is R_5^1.lo.
7. Every class resolves with an even itinerary; `verified` holds for every class; `is_reliable`.

### `test_golden_period3.py`

**p3, `build_period3()` defaults:**
1. Exactly 3 active classes X₀, X₁, X₂, each sourced at an anchor element. There are two distinct inert classes U and W.
2. Words form the orbit-shift chain X₀ → [X₁], X₁ → [X₂], X₂ → [X₀, U⁻¹, W⁻¹]. X₂'s word starts with X₀ in the forward direction.
3. The minimal trellis has no image bridges (the case is closed).
4. `build_period3(blasts=4)` gives the same words spelled in element names.
5. Everything resolves, nothing is unreachable, and `is_reliable` holds.

The element names of X₀, X₁ and X₂ are recorded in Phase 5 for the author's sign-off (E8).

**Nested, `build_nested()` defaults (2 outer blasts):**
1. Every class resolves, no walk is unreachable, at least one walk is unique, and the minimal trellis has image bridges. This is CLAUDE.md's "nested needs two outer blasts".

---

## E. Risks and open items for the author

### Risks

- **R1, runtime.** Decision 10 plus per-check parameterization costs about 160 s for the law tier. Taking the k=2.8 and inversion fixtures off shared scope adds about 15 s even before Phase 4. This is gated by P8 and the Phase 4 gate.
- **R2, unifying the builders can shift facts.**
  - The k10 7+2 and 9-step recipes may differ in stable growth.
  - `infer_iterate_table` is on in some copies and off in others.
  - The nested p1 cutoff is 1e-4 in one recipe and 1e-7 in the other, and the two differ on the repin step.

  P6 measures all of this before Phase 1. A difference is either investigated or kept as an explicit builder parameter, never hidden.
- **R3, inversion signal.** If the inversion case fails early in the pipeline, most of its law parameters become xfail and carry little signal. Strict xfail still turns any fix into a loud XPASS.
- **R4, p3 ids are not deterministic across builds in one process** (regions memory). Never compare ids across builds; golden uses names and relations only.
- **R5, XPASS on the deep p3 runs.** A numerics change could make the known-open deep-p3 xfails pass. Strict mode then fails the suite on purpose: the coder removes the mark and updates the Dev Note.
- **R6, coverage.** Deleting string-pin tests and fan-out mocks drops coverage unless their smoke replacements land in the same commit, which the plan requires.
- **R7, remaining fakes and kernel tests depend on internals.** `FakeDual` still hard-codes `StableNode` internals. The kernel tests that monkeypatch `_image_chain` and `_empty_stretches` depend on the call graph. These are accepted exceptions under decision 6.
- **R8, GPU.** CuPy is installed but a GPU device may not be available. The test skips on device absence, not only on a missing import.
- **R9, blame.** Merged files lose their history regardless. The move-only commit keeps renames detectable for the rest.
- **R10, decision 2 side effect.** Taking `is_reliable` out of the synthetic `test_ambiguous_class_is_unreliable` leaves "ambiguous ⇒ not reliable" untested outside golden.

### Open items

| # | Question | Default if not answered |
|---|---|---|
| **E1** | Law-tier granularity under decision 10 | Per-layer aggregated tests: about 70 cases in about 35 s. The alternative is per-check tests (about 330 cases, about 160 s), or relaxing decision 10 to module scope for the read-only law tier only |
| **E2** | CLAUDE.md says the k=2.8 one-blast class is "resolved via the trellis path", which means `source == "trellis"`. Pin it in golden, despite decision 2's "source not asserted"? | Pin it in golden only |
| **E3** | Where the ledger and Dev Notes live | `docs/test_suite_refactor_ledger.md`, plus a Dev Notes section in the `tests/conftest.py` docstring. The alternative is mirroring into `src` Dev Notes, though decision 13 says leave the code alone |
| **E4** | `cases.py` delegates to `henon_cases` with all keyword arguments explicit, or reimplements the recipes | Delegate for p3 and nested; reimplement k10, k28 and inversion |
| **E5** | Keep the unblasted-nested p1 `area_cutoff=1e-4` variant, or unify on 1e-7 | Keep it as a builder parameter |
| **E6** | Undecided synthesis questions that §B touches: Q13 (uncut-branch pip warning), Q17 (dissipative `per_step_beta`), Q20 (defensive fallbacks: orientation default, missing parent id, `_near_far` cross-branch), the snap tripwire, the verbose print fallback, a designed error for a broken iterate link | Keep the current behaviour with type-, level- or logger-only assertions; delete the snap tripwire and the dissipative test |
| **E7** | The hole-quirk outcome from P1, which decides between the openings law and the xfail | Decided after P1 |
| **E8** | Sign-off on the golden names for p3's X₀..X₂ and on the k10 "folded loop" and k=2.8 one-blast "0 image bridges" facts | – |
| **E9** | Keep the one-rule `advance_key` guard as today's hasattr check, or upgrade it to the behavioural spy | Upgrade |
| **E10** | `__init__.py` packages in the test subdirectories, or unique basenames | Unique basenames, keeping `pythonpath = ["src", "tests"]` |

### Critical files for implementation

- /home/dezhu/Development/TanglePack/tangle-pack/tests/conftest.py
- /home/dezhu/Development/TanglePack/tangle-pack/tests/invariants.py
- /home/dezhu/Development/TanglePack/tangle-pack/tests/walk_helpers.py
- /home/dezhu/Development/TanglePack/tangle-pack/src/tanglepack/examples/henon_cases.py
- /home/dezhu/Development/TanglePack/tangle-pack/tests/test_stable_partition_period3.py
- /home/dezhu/Development/TanglePack/tangle-pack/src/tanglepack/topology/StablePartition.py
- /home/dezhu/Development/TanglePack/tangle-pack/pyproject.toml