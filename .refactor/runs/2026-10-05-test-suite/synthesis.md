# Audit of the tanglepack test suite: what it tests and where it over-tests

Scope: 796 collected test cases in 64 files (about 640 test functions, 18.3k lines). Summed test time is about 43.5 s; the wall clock is about 90 s if every file runs on its own. No single test takes more than 1.5 s. This audit was read-only.

Checks I made against the source myself:
- `IntersectionRegistry.on_interval` is never called anywhere in `src/`.
- The scalar `ManifoldMachine._curvature_area` is defined, and the only remaining callers are tests.
- `BoundaryArc = Arc` is still exported from `tanglepack.loom`.
- `_is_forward_beyond_fundamental` carries the PROVISIONAL note in `StablePartition.py`.
- `uniiterated_bridges` (with the typo) is public API, and `Blast.py` uses it.
- `invalidate_trellises` still exists.

---

## 1. Headline numbers

### 1a. Category (catalog entries, about 640 functions; counts rounded)

| Category | ~Count | % |
|---|---|---|
| algorithm-contract | 330 | 52% |
| implementation-detail | 60 | 9% |
| physical-invariant | 55 | 9% |
| api-facade-caching | 50 | 8% |
| regression-pin | 45 | 7% |
| plotting-visual | 37 | 6% |
| fixture-fact-pin | 30 | 5% |
| smoke | 15 | 2% |
| numerical-accuracy | 8 | 1% |
| test-infra (meta / fixture guards / lints) | 8 | 1% |

### 1b. Recommendation (per function; collected cases in brackets)

| Recommendation | ~Functions | ~Cases |
|---|---|---|
| keep | 285 | 360 |
| keep-but-loosen | 125 | 150 |
| merge | 60 | 65 |
| parametrize | 60 | 70 |
| delete | 48 | 70 |
| ask-author | 30 | 40 |
| rewrite | 25 | 30 |
| move-to-invariant-suite | 22 | 25 |

About 45% of tests can stay exactly as they are. About 35% stay with changes. About 10% go, and the remaining roughly 5% wait on your answers.

### 1c. Brittleness / cost

| Brittleness | ~% | | Cost | ~% |
|---|---|---|---|---|
| low | 50% | | fast (<0.5 s) | 96% |
| medium | 35% | | medium (0.5–1.5 s) | 3% |
| high | 15% | | slow (>1.5 s) | <1% (2 tests) |

**Cost is not the problem.** The problems are brittleness and duplication. The `slow` marker does not track cost: 20 of the 81 tests marked slow take under 0.1 s, and the 6 most expensive tests are not marked.

### 1d. Per-file table (collected cases; "% flagged" = any recommendation other than keep)

| File | n | Dominant category | % flagged |
|---|---|---|---|
| test_topology_plotting | 65 | plotting-visual / impl-detail | 75% |
| test_stable_partition | 48 | algorithm-contract | 55% |
| test_dual_walk | 34 | algorithm-contract (synthetic) | 50% (mostly log-string loosening) |
| test_bridge_class | 30 | algorithm-contract | 60% |
| test_arrangement | 28 | invariant + contract | 50% |
| test_dual_graph | 27 | algorithm-contract | 55% |
| numerics/test_geometry | 27 | algorithm-contract | 25% |
| numerics/test_map_step_and_graph | 26 | algorithm-contract | 45% |
| numerics/test_generation_and_caches | 24 | algorithm-contract | 55% |
| test_symbolic_dynamics | 22 | algorithm-contract | 50% |
| test_pseudoneighbor | 21 | algorithm-contract (synthetic) | 30% |
| numerics/test_workbench_split | 20 | impl-detail | **100% (delete)** |
| test_session_bridge_classes | 19 | facade/caching | 75% |
| test_higher_period_cartoon | 19 | plotting / invariants | 70% |
| test_element_naming | 19 | algorithm-contract (strings) | 60% |
| numerics/test_grow_until | 19 | algorithm-contract | 40% |
| numerics/test_cleanup_walkers_and_examples | 18 | algorithm-contract | 50% |
| numerics/test_bridge_identity | 17 | algorithm-contract | 60% |
| test_fixed_point | 16 | algorithm-contract | 25% |
| test_partition_elements | 15 | algorithm-contract | 55% |
| test_image_cdist | 15 | algorithm-contract | 65% |
| numerics/test_inversion_fixture | 14 | physics / contract | 20% |
| test_resonance_zone_region | 13 | algorithm-contract + golden | 70% |
| numerics/test_intersection_registry_fixes | 13 | algorithm-contract | 85% |
| test_topology_invariants | 12 | algorithm-contract (synthetic) | 40% |
| numerics/test_single_source_of_truth | 12 | algorithm-contract | 70% |
| test_stable_partition_period3 | 11 | invariants + golden | 80% |
| test_session_symbolic_dynamics | 11 | facade/caching | 75% |
| test_session_dual_graph | 11 | facade/caching | 60% |
| test_session_trellis_cache | 10 | facade/caching | 70% |
| test_point / test_branch_point | 10 / 9 | algorithm-contract (linked list) | 70% |
| test_partition_family | 10 | algorithm-contract | 60% |
| numerics/test_workbench_bugfixes | 10 | algorithm-contract | 70% |
| test_iterated_partition | 9 | algorithm-contract (rule reimplemented) | 90% |
| test_minimal_trellis | 8 | physical invariant battery | 60% |
| numerics/test_machine_iterate_invariants | 8 | physical-invariant | 100% (merge/delete) |
| test_manifold_initializer | 7 | fixture-fact | 85% |
| test_session_pseudoneighbors | 7 | facade | 40% |
| test_loom_blast_restore | 7 | algorithm-contract | 40% |
| test_arrangement_sparse | 6 | algorithm-contract (synthetic) | 85% |
| numerics/test_tangle_index_and_orientation | 6 | numerical-accuracy | 50% |
| numerics/test_minimal_solver | 6 | physics | 0% |
| numerics/test_initializer_cdist | 6 | physics | 35% |
| numerics/test_closed_form_curvature | 6 | impl-detail | 100% |
| ~15 small files (1–4 cases each; regression/, gpu, refinement, noise, …) | 38 | mixed | ~55% |

**Worst offenders by share of tests flagged:** test_workbench_split (all of it), test_iterated_partition, test_intersection_registry_fixes, test_stable_partition_period3, test_topology_plotting, test_session_* (all five), test_resonance_zone_region.

---

## 2. What the suite actually tests, layer by layer

| Layer / feature | Primary files | Secondary (overlapping) files |
|---|---|---|
| Geometry kernel (area, point-in-polygon, interior point, arc clipping) | numerics/test_geometry | test_resonance_zone_region, test_arrangement |
| Linked-list data model (Point, BranchPoint, BaseManifold walkers) | test_point, test_branch_point, test_cleanup_walkers | test_segment_index_bookkeeping (insert params) |
| FixedPoint: advance_key, num_branches, per_step_beta, branch_cycle | test_fixed_point, test_map_step_and_graph | test_inversion_fixture, test_workbench_bugfixes, test_image_cdist, test_pseudoneighbor |
| Solver, initializer, growth, refinement | test_minimal_solver, test_initializer_cdist, test_manifold_initializer, test_machine_iterate_invariants, test_growth_integration, test_batched_map, test_closed_form_curvature, test_refinement, regression/* | test_inversion_fixture, test_invariant_helpers, test_gpu |
| Crossings, registry, R-tree index | test_tangle_index_and_orientation, test_tangle_intersection_cdist, test_noise_crossing_collapse, test_intersection_registry_fixes, test_single_source_of_truth, test_segment_index_bookkeeping | test_no_same_stability_crossing, test_generation_and_caches |
| Bridges: identity, genealogy, blasting | test_bridge_identity, test_blast_no_overlap, test_blast_proximity_guard, test_loom_blast_restore | test_single_source_of_truth, test_workbench_bugfixes, regression/test_blast_monotonicity |
| Growth drivers | test_grow_until, test_workbench_bugfixes | test_session_trellis_cache |
| Generation counter and caches | test_generation_and_caches, test_session_trellis_cache | test_session_{bridge_classes,dual_graph,symbolic_dynamics}, test_arrangement, test_grow_until, test_session_strong_pips |
| Strong pips, pseudoneighbors | test_pseudoneighbor, test_strong_pip_periodic_point | test_session_strong_pips, test_session_pseudoneighbors |
| Holes and stable partition | test_stable_partition, test_stable_partition_period3, test_topology_invariants | test_partition_elements |
| Element images (image_cdist / image_of_element) | test_image_cdist | test_partition_elements (near-verbatim twin) |
| Arrangement / Regions | test_arrangement, test_arrangement_sparse | test_minimal_trellis |
| Bridge classes, alphabet | test_bridge_class, test_session_bridge_classes | test_higher_period_cartoon (tangle ordering) |
| Minimal trellis, partition families, iterated cut | test_minimal_trellis, test_partition_family, test_iterated_partition | test_higher_period_cartoon (`_empty_stretches`), test_topology_plotting (brackets) |
| Dual graph | test_dual_graph, test_session_dual_graph | test_topology_plotting |
| Walks, landing, naming | test_dual_walk, test_element_naming | test_dual_graph, test_symbolic_dynamics, test_higher_period_cartoon |
| Symbolic dynamics | test_symbolic_dynamics, test_session_symbolic_dynamics, test_higher_period_cartoon | test_topology_plotting |
| Plotting | test_topology_plotting, test_higher_period_cartoon | test_session_* plot delegates, test_stable_partition |
| Repo hygiene (lints written as tests) | test_workbench_split, test_cleanup_walkers, test_resonance_zone_region (shapely), test_topology_plotting (Phase 8) | — |

### Cross-file duplication clusters (the same fact asserted in several places)

1. **advance_key / branch_cycle:** test_fixed_point (×6), test_map_step_and_graph::test_branch_cycle_is_the_advance_key_chain, test_inversion_fixture::test_advance_key_flips…, test_workbench_bugfixes::test_advance_key_returns_the_next_orbit_point, test_image_cdist::test_period_one_images_stay_on_their_own_branch, test_image_cdist::test_image_cdist_is_additive…, test_pseudoneighbor::test_forward_unstable_branch_cycle_matches_orbit_order.
2. **Growth invariants (monotone cdist, iterate law, one-to-one):** test_machine_iterate_invariants (×7 params plus one extra test), test_growth_integration::test_unstable_invariants_hold_at_every_step, regression/test_cdist_collision_growth, test_invariant_helpers::test_low_stretch_growth_keeps_cdist_injective, test_initializer_cdist::test_iterate_law…, test_inversion_fixture::test_growth_preserves…. These files **disagree on strict versus non-strict monotonicity.**
3. **BridgeId ordering, bridge endpoints and bridges on their own branch:** test_bridge_identity::test_bridge_id_is_endpoints_in_unstable_order, …_ordering_on_period_three, test_single_source_of_truth::test_bridge_endpoints_are_the_ids…, test_workbench_bugfixes::test_bridge_endpoints_stay_on_the_bridges_own_branch, regression/test_boundary_straddle_cdist, test_growth_integration::test_bridges_individually_satisfy_invariants, test_geometry::test_oriented_bridge_polyline_runs….
4. **Partial bridges:** test_bridge_identity::test_partial_pieces_are_held_outside…, test_single_source_of_truth::test_partial_bridge_reports_its_missing_endpoint.
5. **Crossing is one unstable plus one stable, and its cdist is bracketed:** test_tangle_intersection_cdist (×3), test_single_source_of_truth::test_every_registered_crossing_carries_its_unstable_segment, test_no_same_stability_crossing.
6. **"Grow cap raises" and argument validation:** test_grow_until (×6), test_workbench_bugfixes (×3).
7. **Generation bumps / cache invalidation:** test_generation_and_caches (MUTATIONS ×8 + 4 registry-bump tests); test_session_trellis_cache (6 mutation tests); test_session_strong_pips::test_trellis_auto_rebuilds… (an exact duplicate); test_grow_until::test_session_exposes_the_drivers; test_arrangement::test_arrangement_is_cached_by_generation; plus the three 1.3–1.6 s "grow and repartition" tests in session_bridge_classes, session_dual_graph and session_symbolic_dynamics.
8. **Cache hit, rebuild, pip change, repartition at the same generation:** the same 4–5-test pattern in test_session_bridge_classes, test_session_dual_graph, test_session_symbolic_dynamics and test_session_trellis_cache.
9. **Session result equals direct build:** session_bridge_classes (×3, including the private `_gathered_partitions`), session_dual_graph (×2), session_symbolic_dynamics (×1), partition_family::test_from_results_signature….
10. **Session plot delegate wiring (monkeypatch):** test_session_dual_graph::test_session_plot_delegates_to_plotting, test_session_symbolic_dynamics::…forward_to_plotting, test_topology_plotting::test_session_plot_fanouts_go_through_fanout_plot (×5), …call_fanouts… (×4), test_trellis_plotters_delegate_to_plotting.
11. **Orbit hole-side invariant (I1):** test_topology_invariants::test_henon_holes_share_bridge_side (already run inside `punch_holes`), test_stable_partition_period3::test_p3_holes_share_bridge_side, …deep_p3…, test_stable_partition::test_henon_propagated_hole_side_matches_coords.
12. **Same-branch bridge rows (I2):** test_topology_invariants::test_henon_bridge_rows_consistent, test_stable_partition_period3::test_p3_bridge_rows_consistent, test_bridge_class::test_same_branch_bridges_never_mismatch_on_{k10,p3}.
13. **Direct-hole side:** test_stable_partition::test_henon_direct_hole_side_reproduces_inward_pair, test_stable_partition_period3::test_direct_hole_side_is_the_side_of_its_coordinates (×4), the p3/k28 hole snapshots.
14. **Provisional forward-hole exemption:** test_stable_partition::test_forward_pairs_beyond_the_fundamental_segment_get_no_hole (×10), test_p3_forward_holes_stop_at_the_branch_return, test_stable_partition_period3::test_p3_hole_sides_are_pinned.
15. **Each crossing has exactly one owner, and the partition covers the branch:** test_partition_elements (×2 fixtures), test_stable_partition::test_henon_partition_covers_branch, test_partition_family::test_owner_of_intersection…, the `_check_iterated` helper.
16. **image_of_element (bounded elements have an image, image covers both endpoint images, fallback):** test_image_cdist (×4) against test_partition_elements (×3), described as a "local twin".
17. **Region soundness, disjointness, Euler, image agreement, k10 against p3:** test_arrangement (4 twin pairs).
18. **Dual graph twins (k10 against p3):** structure, face side, payload, unification, marks round trip (partition_family), iterated invariants (×3 fixtures).
19. **The k=10 class picture (2 classes, 1 active, 4 members, folded loop, letter 'a'):** test_bridge_class (×3), test_session_bridge_classes (×5, including the 2.4 s zone test).
20. **The k=10 iterated element set (11 names) and bracket rows:** test_element_naming::test_k10_iterated_names, …_k10_homotopy_names, …_k10_describe_lists…, test_dual_graph::test_k10_stable_node_names…, test_iterated_partition::test_k10_empty_stretch_cut_reads_as_expected, test_topology_plotting::test_cartoon_brackets_match_closedness (second half).
21. **The k=10 word `a → a u^-1 a^-1` and the [[1,1],[1,1]] matrix:** test_symbolic_dynamics (synthetic ×3, fixture ×3), test_session_symbolic_dynamics (×3, including describe), test_dual_walk::test_k10_landings_and_the_active_class_walk, test_topology_plotting::test_plot_transition_graph… / …itinerary_table….
22. **Letters / naming notation strings:** test_element_naming (×6), test_higher_period_cartoon::test_nested_names_carry_fixed_point_letters, test_p3_names…, test_bridge_class::test_element_ref_hashes… (ElementRef.label), test_topology_plotting::test_name_mathtext (×11).
23. **Period-3 / nested words:** test_higher_period_cartoon (equivariance, four blasts, nested default, nested inner equals period 3, blast order), test_session_symbolic_dynamics::test_p3_symbolic_dynamics_smoke (on a *different* nested recipe).
24. **Cross-tangle "no class mixes fixed points":** test_bridge_class::test_p3_classes_cover_both_tangles_and_mix_neither and test_session_bridge_classes::test_p3_no_class_mixes_fixed_points (verbatim).
25. **Refactor tombstones ("X is gone"):** test_bridge_identity::test_deleted_genealogy_attributes_are_gone, test_workbench_bugfixes::test_grown_until_intersection_is_gone, the hasattr check in …advance_key_returns…, test_map_step_and_graph::test_workbench_has_no_private_key_advance, test_single_source_of_truth::test_tangle_keeps_only_index_state and the hasattr checks in …bridge_endpoints…, test_cleanup_walkers::test_string_dispatched_iter_method_is_gone, test_resonance_zone_region::test_boundary_arc_dataclass_is_gone / boundary_intersection_id, test_fixed_point::test_fixed_point_takes_no_branch_count, test_topology_plotting::test_hole_style_conventions_live_once / test_stable_partition_module_holds_no_drawing_code.
26. **Fixture recipes:** the k=10 "9u + turnaround + intersect + trim + bridges" recipe exists **9 times**. Two of those are `k10_session` fixtures that shadow conftest's, and the recipes differ silently in `infer_iterate_table`. There are two different "nested p3" builds (conftest p1 cutoff 1e-4 with 0 blasts, against henon_cases 1e-7 with 2 blasts), and facts are pinned on both. The dual-graph assembly helper is copied 4 times. `_define_zone` is copied twice.

---

## 3. Tests that should not be enforced (grouped by theme)

| Theme | ~Count (cases) | Representative tests |
|---|---|---|
| **A. Refactor tombstones and code-organisation pins** (hasattr-absent, source greps, line counts, signatures, module layout) | ~40 | all of test_workbench_split (20); test_deleted_genealogy_attributes_are_gone; test_grown_until_intersection_is_gone; test_workbench_has_no_private_key_advance; test_tangle_keeps_only_index_state; test_collect_is_the_only_traversal; test_string_dispatched_iter_method_is_gone; test_stability_alias_is_defined_once_in_numerics; test_fixed_point_takes_no_branch_count; test_boundary_arc_dataclass_is_gone; test_every_plotter_is_a_module_function_in_plotting; test_stable_partition_module_holds_no_drawing_code; test_hole_style_conventions_live_once; test_image_chain_is_the_bridge_class_function; test_kinds_are_distinct; test_branch_cycle_reproduces_the_topology_branch_orders |
| **B. Migration golden masters whose job is finished** (bit-exact or 6-digit floats, area to 1e-11, bit masks) | ~8 | test_bridge_cutting_pin_k10, …_period_3; test_new_grow_manifold_matches_old_period_{one,three} (60 lines of frozen legacy code); test_k10/p3_zone_area_and_containment_are_unchanged; test_per_step_beta_matches_the_expression_it_replaces; test_per_step_beta_is_the_k_th_root_only_without_inversion |
| **C. Superseded rules or superseded oracles** | ~10 | test_direction_matches_the_element_order_at_the_ends (uses `element_sort_key` and not `anchor_outward_key` from 2026-09-30); test_classes_and_members_come_out_in_the_documented_order (no tangle grouping); test_henon_propagated_hole_side_matches_coords (whole-bridge reading, replaced by the 10-02 image sub-arc); test_direct_hole_side_is_the_side_of_its_coordinates ×4 (old geometric oracle against the new crossing-sign rule); test_iterated_partition `_empty_stretches` test reimplementation (missing the 09-30 branch restriction on chords and cross-branch flanks), used by test_{k10,p3,k28}_iterated_partition_invariants and test_cuts_record_the_far_end…; test_low_stretch_growth_keeps_cdist_injective plus the `strict=True` default in `invariants.py` (the memory notes say ties are legitimate); test_bridge_identity_on_the_inversion_saddle (calls the per-branch-pair anchors a "defect") |
| **D. Provisional rules** | ~14 | test_forward_pairs_beyond_the_fundamental_segment_get_no_hole (the 1..k-1 cases); test_p3_forward_holes_stop_at_the_branch_return (`== [1, 2]`); test_p3_hole_sides_are_pinned; test_refined_cdist_is_mean_of_neighbours (its own docstring calls the rule suspect); test_faces_closed_raises_at_the_cap_for_the_anchor (pins the deferred anchor limitation); test_mixed_eigenvalue_signs_are_rejected / test_set_k_value_rejects_disagreeing_eigenvalue_signs (a limitation or a design rule?) |
| **E. Exact fixture facts (counts, sets, words, letters, snapshots)** | ~35 | test_small_tangle_crossing_count_is_unchanged; test_k10_intersection_graph_is_unchanged (8/18/7/7/4); test_initialization_{unstable,stable} (exactly 3 points); test_k28_two_blast_hole_sides_are_pinned; test_k28_blast_child_gets_no_forward_hole (exactly 4 holes, 2 refs); test_k10_has_one_active_and_one_inert_class; test_k10_iterated_names; test_k10_empty_stretch_cut_reads_as_expected; test_k10_active_class_word; test_k28_two_blasts_structure (letters a/b/c/u, against its own docstring); test_p3_words_survive_four_blasts ('w'); test_nested_default_words_are_pinned (d/e/f/'uu' at commit 2315204); test_signs_alternate_along_a_stable_branch; test_p3_arrangement_is_two_components…; test_k10_rejects_an_image_that_is_only_a_sub_face; test_the_snap_does_not_fire_on_the_default_path; test_active_class_lies_in_the_zone… (counts) |
| **F. Exact strings: describe(), repr, summary, legend, title, labels, ElementRef.label** | ~40 | test_describe_reports (stable partition); test_table_lookups_symbols_and_report; test_describe_reports_letters_and_inertness; test_describe_reports_the_image_evidence; test_k10_summary_reports_counts; test_walk_search_repr_and_walk_repr; test_describe_has_one_line_per_element; test_describe_mentions_parents_and_cuts; test_k10_describe_lists_every_parent…; test_homotopy_only_naming_is_its_own_parent (double-spaced describe line); test_is_reliable_and_describe (including `'spectral' not in text`); test_sparse_and_dense_agree… ('+2 dangling ends'); test_dual_graph_legend_handles…; test_cartoon_legend_handles…; test_line_cartoon_labels… ('anchor  0.0'); test_plot_transition_graph… (title); test_k10_stable_node_label_joins_the_names (' \| ') |
| **G. Log-message and exception-message text** | ~35 | test_rebuild_drops_metadata… ('no longer name the same crossings'); test_grow_until_clears_the_iterated_flag… ('recut the bridge set'); test_partition_warns_without_pseudoneighbors (exact count 2); test_strong_pip_cuts_warn… / …_are_silent…; test_multiple_start_faces_are_all_searched ('faces 2 faces'); test_singleton_owned_by_a_non_singleton… ('NON-singleton', 'backward only'); about 10 more test_dual_walk caplog checks; test_split_origin_sides_raise / test_bridge_rows_inconsistent_raise (message substrings); test_every_bridge_must_be_consecutive… ('not consecutive crossings'); test_collect_rejects_a_broken_iterate_link (matches an accidental 'NoneType'); test_no_partitions_raises_and_warns; test_a_class_straddling_zones_warns…; test_unreadable_image_chain…; test_virtual_class_named_new1… |
| **H. Private internals** (private attributes, monkeypatched private helpers, private index dicts) | ~45 | `_map_batchable` ×3; `_get_lambda_u` ×4; `Tangle._manifold_segs/_seg_lookup/_intersecting_segments` (×6 tests, including the trim oracle); `_frozen_ids`; `_built_generation`; `session._bridge_classes` dict size; `_gathered_partitions` / `_partition_signature`; `_fanout_plot` / `_fanout_call` (×10); Arrangement `_nodes/_half_edges/slots` (test_sparse_kept_endpoints…, `_check_minimal`); `_region_key`; `_stable_arc_midpoint` ×2; `_snap_to_partition_boundary`; `_scaled_image_cdist` ×2; monkeypatched `_empty_stretches` / `_image_chain` / `_inert_letter` |
| **I. Memoisation-strategy and cache-policy pins** | ~8 | `is`-identity in test_bridge_point_array_is_memoised / test_branch_position_map_is_memoised; test_bridge_point_arrays_survive_an_iterate (every version must bump); test_rebuild_flag_forces_fresh… (rebuild cascades upstream); test_restore_without_recompute_leaves_bridges_stale (id identity) |
| **J. Plot geometry, palette and font metrics** | ~20 | test_circle_layout_puts_every_node… (alpha=0, sweep=π, offset radii, gids); test_p3_circle_follows_the_zone_boundary (π/3); test_nested_outer_arcs_go_around… (gid parsing, len(vertices)>4); test_clip_to_arcs_true… (reimplemented padding formula); test_plot_itinerary_table_columns_font_and_rows (60 pt ±5%, floor 4); test_plot_itinerary_table_lists_every_class (width heuristic); test_walk_zorder_sits_between…; test_class_colors_are_assigned_in_fixed_order (slot order); test_plot_dual_graph_scatters… (artist index order); test_nested_classes_are_lettered_and_coloured… (colour families) |
| **K. Exact float tolerances and wall-clock checks** | ~8 | test_registry_insert_is_near_linear (wall clock; has flaked); test_refinement_essentially_converges_to_cutoff (2% / 5×); per_step_beta exact `==` ×2; test_one_map_step_scales… (rel 1e-4 tied to seed step); test_scaled_end_tolerance_is_relative (assumes SCALING_RTOL ≥ 1e-3); the unlinked halves of the table-linked tests (collision_rtol < 2%) |
| **L. Questionable or dead API pinned** | ~8 | `on_interval` / `_get_lambda_u` ×5 (no caller in src; uses full λ, not per_step_beta; includes a dubious cross-saddle fallback); test_per_step_beta_stable_uses_the_stable_eigenvalue… (dissipative maps); test_invalidate_trellises_is_a_deprecated_no_op; test_curvature_area_batch_matches_scalar (the scalar path is used only by tests); test_orientation_defaults_to_preserving_without_evidence; test_missing_parent_id_reads_as_own_parent; test_stable_node_label_falls_back… |
| **M. Repo-hygiene lints written as tests** | 4 | test_no_test_or_script_defines_its_own_henon; test_shapely_is_not_a_dependency; test_stability_alias_is_defined_once…; test_iterate_bridge_still_patchable_on_the_workbench (guards another test's technique) |

---

## 4. Weak, tautological and dead tests

**Tautological:**
- test_point::test_get_point compares `.all() == .all()`, which is `True == True`.
- test_scaled_element_image_stays_on_the_advanced_branch selects by the same key it then asserts.
- test_kinds_are_distinct.
- test_image_chain_is_the_bridge_class_function.
- test_names_are_deterministic_across_builds (builds twice from the same objects).
- test_session_result_matches_gathered_partitions_helper.
- test_machine_initialization.
- test_alpha_matches_distance_ratio, second assertion (k_value = 1).
- The class_sort_key concatenation check in test_classes_and_members….
- The neighbour block in test_regions_at_and_bounded_by….
- test_image_cdist_is_additive_in_n….
- test_henon_reference_holes_hug_the_stable_manifold (re-runs the production helpers).
- test_k10_payload (re-derives each field with the same functions).
- test_bridge_classification_uses_the_memoised_midpoint (classifies against its own helper).

**Weaker than the name claims:**
- test_oriented_bridge_polyline_runs_in_the_dynamical_direction: never checks the orientation, and skips on None.
- test_min_separation_drops_close_bridges: no distance check.
- test_blast_recognizes_already_known_bridges: only checks counters.
- test_trim_stable_manifolds_reads_the_registry: a no-op trim would pass.
- test_kevin_way_period_three_orbit_chain: the chain order is never asserted.
- test_trajectory_extension_forward_with_dedup: no duplicate input is ever built.
- test_pair_on_different_unstable_branches_rejected: one negative assertion on a muddled fixture.
- test_no_same_stability_crossing: only self-crossings, only on p3.
- test_k10_region_images_agree_with_the_dynamics: can pass vacuously (no `with_image > 0` guard).
- test_p3_merges_the_inner_outer_face…: does not check which face it merged into.
- test_row_of_end_answers_at_an_anchor: only checks membership in {left, right}.
- test_k28_iterated_partition_invariants: the one-blast fixture has 0 image bridges, so most checks are vacuous.
- test_clip_to_arcs_false…: checks inequality only.
- test_compute_single_fixed_point_returns_list: `all()` over a possibly empty list.
- test_k10_stable_node_names…: a subset check.
- test_blast_completes_without_monotonicity_failure: asserts nothing about monotonicity.
- test_high_stretch_period3_growth: nothing proves it still reaches the near-ULP regime.

**Dead or never run:**
- test_enable_gpu_without_cupy_raises_clearly always skips here because CuPy is installed.
- `tests/visualizations/*.py` is never collected, crashes with `KeyError: None`, and matches a log string that no longer exists.
- The skip guards in test_k10_a_different_pip_moves… and test_iterate_returns_none… never fire.
- The `skipif` scaffolding in test_session_plot_delegates_draw… is unreachable.
- The `slow` and `regression` markers are never selected by anything.

**Redundant because the library already checks it:** test_henon_holes_share_bridge_side. `punch_holes` runs the same check internally.

---

## 5. Coverage gaps worth closing (important ones only)

1. **The fundamental invariant is checked weakly.** Nothing asserts that every registered crossing is (unstable `manifold_a_key`, stable `manifold_b_key`) on the p3, nested and blasted fixtures. Nothing checks unstable-vs-unstable crossings across branches or fixed points, or stable-vs-stable crossings. Nothing tests that same-stability pairs are discarded and logged, for example with a synthetic u×u near-tangency.
2. **Area preservation is applied only to k=10 iterate chains.** It is never checked on p3, inversion or blasted fixtures, and never against the area of a ResonanceZone or Region.
3. **The cross-branch anchor-outward orientation (2026-09-30)** and the **tangle grouping / connecting classes (tangle=None)** have no test. A hand-built pair of ElementRefs on two branches would cover them cheaply.
4. **The iterated cut's 2026-09-30 rules** have no test: chords are paired only on the base branch, a chain ending on another branch marks nothing, and a cross-branch flank is skipped with a WARNING. The **k=2.8 two-blast cut facts** (R_1^1=[0,10]…) are documented but nothing tests them.
5. **Orientation-reversing maps end to end.** This is covered only by synthetic I1/landing tests. No real det J<0 run flips hole sides, rows or words. Whether this is wanted depends on Q16.
6. **min_separation guarantee:** no test checks that kept siblings are at least min_separation apart. **already_known:** no test checks that known bridges are not re-iterated.
7. **preserve_ids:** no direct check that ids and cdists are identical before and after `compute_intersections(preserve_ids=True)`. No check that `compute_intersections` is idempotent.
8. **Refined-point cdist lies strictly between its neighbours.** This is the real invariant behind the suspect mean-of-neighbours rule.
9. **Real-data symbolic paths:** an ambiguous walk, a virtual class as a sink in the graph or matrix, a non-empty unmatched_members, an evidence-length-mismatch WARNING, and refinement with three or more occurrences.
10. **The symbolic_dynamics cache** is not invalidated by a partition-signature change at the same generation.
11. **Period-k orbits:** `f^k(p)=p` and the eigenvalue product are not checked for period-k solver output, and the initializer's alpha/k relation is not checked for k>1.
12. **A deterministic hand-built sub-face case** for the area rejection in `Arrangement.image_of`. Today only a k=10 fixture accident exercises it.
13. **A shared invariant checker** that runs on every fixture: even itineraries, same-side pairs, a unique owner per crossing, partition coverage, I1/I2, and the BridgeId/bridges_at consistency checks.

---

## 6. Questions for you

Each question lists candidate answers and roughly how many test cases it decides.

### A. How to treat fixture facts and golden values

**Q1. Golden masters in general.** For exact numeric or structural snapshots captured at a commit (bridge-endpoint cdists, zone areas and masks, k10 graph counts, hole-side snapshots): (a) delete them all and rely on invariants; (b) keep a small "golden" tier marked `@pytest.mark.golden` and re-record it when numerics change (with a `--update-golden` helper); (c) keep only the ones you have verified scientifically. *Decides about 15 cases:* bridge_cutting_pin_k10/p3, zone_area ×2, k10_intersection_graph, small_tangle_crossing_count, p3/k28 hole-side snapshots, initialization 3-point ×3, test_k10_rejects_an_image….

**Q2. Your CLAUDE.md fixture facts** (k=10 brackets `[ ] ( ) …`, the 11 k=10 iterated names, `a → a u^-1 a^-1`, the k=2.8 two-blast words and matrix, the p3 `a→b→c→a u^-1 w^-1`): (a) these are acceptance facts, so pin each one **once** in a dedicated `tests/golden/` file; (b) pin them, but as a separately marked suite you expect to update when a rule changes; (c) keep them only in CLAUDE.md and drop them from tests. *Decides about 25 cases spread over 8 files* (duplication clusters 19–23).

**Q3. Letters in pinned words.** test_k28_two_blasts_structure pins a/b/c/u although its own docstring says "never pin a letter". Are the cdist-ordered letters (and the inert series u, v, w, …, uu) a contract? (a) Yes, pin letters; (b) no, compare words spelled in element names (letter-free); (c) pin letters only for the anchor class 'a'. *Decides about 6 cases:* k28_two_blasts, p3_words_survive_four_blasts, nested_default_words_are_pinned, k10_active_class_word, k28_letters_only…, inert_letter_series.

**Q4. The nested words at commit 2315204** (`d → d uu^-1 e^-1, e → f, f → d uu^-1 d^-1`): (a) a verified scientific result, so keep it; (b) only a snapshot, so replace it with letter-free structural checks (inner words equal the p3 words, which already exist); (c) delete it. *Decides 2 tests.*

**Q5. Pinned hole openings.** The p3 and k28 hole snapshots freeze a single `('outward','right')` opening on the iterate<0 'left' holes. (a) Understood and correct; (b) unexplained, so investigate before pinning anything; (c) drop the snapshots and keep the rule-derived checks (inward pair for direct holes, one side per orbit). *Decides 2 tests, and may expose a bug.*

**Q6. k=10 class picture.** (2 classes, 4+3 members, folded loop, unresolved inert member.) (a) One consolidated golden test; (b) incidental to the 9-step build, so drop it; (c) keep only the folded-loop and inertness evidence. *Decides about 8 cases.*

### B. Provisional and in-flux rules

**Q7. The period-k +1..+(k-1) forward-hole exemption.** (a) Pin only the firm rule (iterate ≥ k_value never punched); (b) pin the 1..k-1 cases as `xfail(strict=False)` with a "provisional" reason; (c) pin all of it as now and update when you decide. *Decides about 13 cases* (10 parametrized cases, p3_forward ==[1,2], the p3 snapshot, test_single_reference_window_per_fixed_point).

**Q8. Empty-stretch cut oracle.** The test helper `_empty_stretches` reimplements the rule, and it already lacks the 09-30 refinement. (a) Generic tests check only rule-independent invariants (refinement inside the parent, homotopy boundaries kept, a unique owner), and the rule itself is pinned through fixture facts plus a few hand-built unit cases; (b) keep the reimplementation and update it to 09-30; (c) test the production `_empty_stretches` directly with synthetic lobes only. *Decides about 8 cases* (the iterated-invariants tests ×3, cuts_record_the_far_end, abutting_no_empty_stretch, two cartoon-file chord tests, plus moving them).

**Q9. Refined-point cdist.** The mean-of-neighbours rule is called "suspect": (a) delete the test; (b) replace it with "refined cdist lies strictly between its neighbours"; (c) the rule is settled, so keep it. *Decides 1 test plus 1 new one.*

**Q10. cdist strictness.** Your memory notes say ties are legitimate. (a) Every growth invariant uses `strict=False`, and the helper default changes to False; (b) strict only for k=10 low-stretch unstable growth; (c) also drop the injectivity check (`assert_no_cdist_collision`). *Decides about 12 cases* (machine_iterate_invariants ×7, growth_integration, low_stretch_injective, fundamental_segments_injective, inversion growth ×4 that use the strict default).

**Q11. Inversion-saddle anchors.** Each unstable branch carries two (0,0) anchors, so a zero-length anchor-to-anchor bridge gets cut. (a) Acceptable by design, so assert it and drop the "defect" wording; (b) a bug, so the test should fail on it (xfail until fixed); (c) remove those bridges in `create_bridges`. *Decides 1 test, possibly a code change.*

**Q12. Anchor faces-closed limitation.** (a) Keep a test pinning that an anchor never closes its faces; (b) test the cap with a cheap impossible target and document the limitation in the Dev Note only. *Decides 1 test (2 s).*

**Q13. A period-k pip that leaves k-1 branches uncut.** Is that still worth a WARNING after the pip-branch-only rule? (a) Yes, and test it behaviourally; (b) no, demote it to DEBUG and delete the 2 tests; (c) keep the warning but don't test its text. *Decides 2 tests, possibly a code change.*

### C. Plotting tests

**Q14. Plot test depth.** (a) Smoke only (runs, returns the right object, draws N artists for N inputs); (b) smoke plus topological properties (node on the correct side, face point inside its region, bracket glyph matches closedness, one node per element, walk goes element to element); (c) freeze geometry and styles as well. My recommendation is (b). *Decides about 50 cases* across test_topology_plotting, test_higher_period_cartoon and the session plot delegates.

**Q15. Plot strings and style constants.** Legend labels, titles, label separators, z-order ladder, palette slot order, colour families per tangle, itinerary column widths and font shrink: (a) drop all such pins; (b) keep only "class_colors never repeats" plus "legend matches the plotted artists (count and colour)"; (c) keep as is. *Decides about 20 cases.*

### D. Scope of the library (dead or unsupported API)

**Q16. Orientation-reversing maps** (mixed-sign eigenvalues are rejected; det J<0 code paths exist): (a) a permanent decision, so keep the rejection tests and drop the reversing-path tests; (b) planned support, so mark the rejection as an xfail limitation and add an end-to-end b<0 fixture; (c) keep the synthetic reversing-path tests as they are. *Decides about 8 cases.*

**Q17. Dissipative maps.** Is `per_step_beta` on a |det J|≠1 pair a contract? (a) No, area-preserving only, so delete; (b) yes. *Decides 1 test.*

**Q18. `IntersectionRegistry.on_interval` / `_get_lambda_u`** (no caller in src; uses the full λ, not per_step_beta; falls back to the other saddle's eigenvalue): (a) remove the API and its 5 tests; (b) keep it, fix it to per_step_beta, and raise or skip instead of the cross-saddle fallback; (c) keep as is. *Decides 5 tests.*

**Q19. Deprecated and compatibility aliases.** (`invalidate_trellises`, `BoundaryArc`, the `forward_*_branch_cycle` wrappers, the scalar `_curvature_area` used only by tests, `use_table=False` plus the snap tolerance as a supported mode.) (a) Remove them all and their tests; (b) keep them, test only that they work; (c) decide one by one. *Decides about 10 cases.*

**Q20. Defensive fallbacks.** Should these be silent defaults or errors? Orientation when there is no evidence, an iterated element with no parent_element_id, a StableNode side with no name, the cross-branch `_near_far` fallback. (a) Errors (and the tests check the raise); (b) silent defaults (keep the tests); (c) warn. *Decides about 5 cases.*

### E. Regression tests and the history of a change

**Q21. Refactor tombstones and code-organisation tests** (cluster 25 plus test_workbench_split): (a) delete them all; (b) keep the few that protect a real bug (say which). *Decides about 40 cases.*

**Q22. Repo-hygiene rules** (one Henon definition, no shapely, one Stability alias, no pyplot in topology logic modules): (a) drop them; (b) move them to a lint or CI step (ruff, import-linter, a grep script); (c) keep them as a single `tests/test_repo_hygiene.py`. *Decides 4–5 tests.*

**Q23. Migration-equivalence tests** (frozen pre-1.8 grow loop, pre-Phase-2 cdist fingerprints, "matches the expression it replaces"): may I delete them? *Decides about 6 cases.*

**Q24. Bug-fix regressions.** For regression tests that go through private helpers (`_bridge_unstable_span`, `_containing_bridge`, `_backward_endpoint`, `_region_key`, noise collapse): (a) keep the private unit test as long as the bug was subtle; (b) replace each with a public end-to-end check where one exists (for example `test_p3_propagated_holes_land_on_the_predicted_branch`), and keep the private test only when no public route exists. *Decides about 8 cases.*

**Q25. The `regression` marker and the `tests/regression/` directory.** (a) Keep the directory, drop the marker; (b) keep the marker, dissolve the directory into the per-module files; (c) drop both and name regression tests `test_regression_<bug>_…` inside the module files. *Decides 9 marked tests plus about 25 unmarked bug-fix tests.*

### F. Numerical accuracy

**Q26. Tolerance-dependent tests.** (Refinement leaves pairs below area_cutoff within 2%/5×, per_step_beta measured to 1e-4, SCALING_RTOL ≥ 1e-3, collision_rtol < 2%, the 1e-3/0.1 margins for iterate images.) (a) These are guarantees, so keep them; (b) they are diagnostics, so keep them with tolerances loosened by about 10× and a comment saying where each number comes from; (c) delete them. *Decides about 8 cases.*

**Q27. Exact float equality** (`per_step_beta == |λ|^(1/p)` with rel=0): may I switch these to `approx`? *Decides 3 cases.*

**Q28. Wall-clock performance guard** (registry near-linear insert): (a) move it to an opt-in `@pytest.mark.perf`; (b) delete it; (c) keep it in the default run. *Decides 1 test.*

### G. Cache and facade tests

**Q29. Cache contract.** (a) One parametrized cache-contract suite over every cached session product (trellis, arrangement, bridge_classes, minimal_trellis, iterated_partition, dual_graph, symbolic_dynamics) × {hit, rebuild, generation bump through a cheap mutation, partition-signature change, pip change}, with an expected hit/miss table; (b) keep per-file tests. *Decides about 40 cases, and saves about 5 s.*

**Q30. rebuild=True cascade.** Must `dual_graph(rebuild=True)` and `symbolic_dynamics(rebuild=True)` also rebuild upstream products? (a) Yes, it is a contract; (b) only the product itself, so don't test the cascade. *Decides 2 tests.*

**Q31. Memo invalidation.** (a) Only "memos never go stale" (equality with a fresh walk); (b) also the `is`-identity and "every version bumps" pins. *Decides about 4 cases.*

**Q32. Facade delegation tests** (session method forwards kwargs to plotting or the workbench via monkeypatch): (a) delete them and keep only behavioural session tests; (b) keep one generic delegation smoke test. *Decides about 15 cases.*

**Q33. Generation bump on a deduplicated add, and the list of read accessors by name.** Is "no spurious bump" a requirement, or is "every mutation bumps" enough? *Decides 3 cases.*

### H. Report, log and error text

**Q34. `describe()` / `summary()` / repr output.** (a) Smoke only (non-empty, mentions each fixed point or class); (b) a stable format, so pin it with snapshot files; (c) pin key phrases only. *Decides about 20 cases.*

**Q35. Log assertions.** (a) Assert only the level plus the logger (and a structured field when there is one), never the wording; (b) assert one keyword per message; (c) as now. *Decides about 35 cases.* A related point: should the `verbose` print fallback when logging is unconfigured be guaranteed? *Decides 3 more.*

**Q36. Exception messages.** (a) Assert type only; (b) type plus one keyword; (c) as now. Also, should a broken iterate link raise a *designed* error instead of the accidental 'NoneType' one? *Decides about 25 cases plus 1 code change.*

**Q37. Notation stability.** Is `R_(o.b;i)^j` with mathtext `$R_{o.b;i}^{j}$`, the `A:` prefix and `ElementRef.label` (`A:p3@1.0/L#2`) now stable? (a) Yes, pin one canonical example per form; (b) pin only the structural parts (side, index, superscript, branch, letter); (c) `ElementRef.label` is only a debug string, so don't pin it. *Decides about 20 cases.*

### I. Slow tests and markers

**Q38.** The whole suite takes about 45 s. (a) Drop `slow` entirely; (b) define slow as more than 0.5 s, applied mechanically (about 10 tests), and add a CI or `make test-fast`; (c) keep it as is. *Decides 81 marked tests.*

**Q39. Fixture scope.** May read-only built fixtures (k10_partitioned, p3_partitioned, k28 ×2, inversion, nested) become session-scoped and *frozen* (tests that mutate must request a `fresh_*` factory)? This saves about 15–20 s. Doing it safely needs a guard that fails if a test mutates a shared build (for example by checking `workbench.generation` at teardown). Is that acceptable?

### J. Test layout and fixtures

**Q40. One case factory.** May I add `build_k10(...)`, `build_k28(blasts=n)` and `build_inversion()` to `examples/henon_cases.py`, and retire the 9 copied k10 recipes and the conftest `henon_p3_session` nested build in favour of `build_nested`? Or should test fixtures stay *independent* of the figure-script builders, so that changing a figure default never moves a test? Options: (a) shared builders with defaults frozen for tests; (b) a tests-only `tests/cases.py`; (c) as now. *Decides the fixture strategy for roughly 580 fixture-using tests.*

**Q41.** `henon_p3_session.inner_zone` is actually the period-1 *outer* zone. Rename it, or drop it with the fixture?

**Q42. Fakes** (`walk_helpers.FakeFixedPoint.advance_key` and `FakeTrellis.image_cdist` reimplement the one-rule functions): (a) the fakes delegate to the real `FixedPoint` and a real minimal trellis; (b) keep the fakes but add a contract test that fake equals real; (c) as now.

**Q43. `minimal_helpers.build_pieces`** (hand assembly of the session pipeline; it returns an un-lettered table and the tests work around that): (a) replace it with the session API everywhere; (b) keep it for direct-assembly tests and fix the lettering at its source.

**Q44. Directory layout.** Do you want the test tree to mirror `src/` (`tests/numerics/`, `tests/topology/`, `tests/loom/`, `tests/plotting/`) with separate tiers (`invariants/`, `golden/`, `regression/`)? Or a flat layout with markers?

**Q45. `tests/visualizations/`.** It is broken and never collected. (a) Delete it; (b) move it to `scripts/archive/`.

**Q46. Private-API tests in general.** (a) Allowed only for pure algorithmic kernels (for example `_resolve_inertness`, `_do_segments_intersect`, `_side_of`, `_backward_endpoint`); (b) never; (c) as now. *Decides about 45 cases.*

---

## 7. Proposed target layout for the rewritten suite

```
tests/
  conftest.py              # markers, frozen session-scoped cases, fresh_* factories, mutation guard
  cases.py                 # (or henon_cases) build_k10, build_k28(blasts), build_period3, build_nested, build_inversion, synthetic builders
  helpers/
    invariants.py          # strict=False default; every check_* below lives here
    fakes.py               # fakes delegating to real FixedPoint rules
  invariants/              # TIER 1: property suite, parametrized over ALL_CASES
    test_manifold.py       # monotone cdist (non-strict), no spikes, iterate law, one-to-one, area along chains
    test_crossings.py      # every crossing unstable x stable, cdist bracketed, signs = cross product, no u×u/s×s
    test_bridges.py        # BridgeId order, consecutive on own branch, single copy, non-overlap, bridges_at exact, image/preimage round trip
    test_partition.py      # covers branch, unique owner per side, singletons, I1 hole sides, I2 rows, propagation lands on predicted branch
    test_arrangement.py    # Euler per component, regions disjoint/sound, image_of agrees with map + area, preimage∘image
    test_classes.py        # every bridge in one class, anchor-outward orientation, row=geometry, no tangle mixing (when no heteroclinic)
    test_dual_graph.py     # wall/unified structure, bipartite degree rule, face side = geometry, unified = pip segment
    test_symbolic.py       # even itineraries, same-side pairs, refined keep identity, children inherit word, inert out of matrix
  unit/                    # TIER 2: contract tests per module, synthetic inputs, no Hénon build where avoidable
    numerics/  (geometry, point, fixed_point, solver, initializer, refinement, tangle_kernel, noise_collapse, registry, iterate_table, generation)
    topology/  (pseudoneighbor, strong_pip, stable_partition_intervals, holes, arrangement_hand_built, arrangement_sparse,
                bridge_class_inertness, partition_family, iterated_cut_synthetic, element_naming, dual_walk, landing, trellis_itinerary, symbolic_translate_refine)
    loom/      (blast_errors, restore, resonance_zone, alphabet)
  facade/                  # TIER 3: session = direct build (one param per product), cache-contract table, fan-out shapes
    test_session_equivalence.py
    test_session_caches.py # products × {hit, rebuild, gen bump, partition sig, pip change}
    test_growth_drivers.py
  golden/                  # TIER 4: author-approved fixture facts, each pinned ONCE, @pytest.mark.golden
    test_k10.py            # brackets, iterated names, word, matrix
    test_k28.py            # one-blast singleton path; two-blast cut + words + matrix
    test_period3_nested.py # orbit-shift equivariance, inner = p3, blast-order commutes (+ letter pins per Q3/Q4)
  regression/              # TIER 5: one file per bug, public route preferred, named test_<bug>_…
  plotting/                # TIER 6: smoke + topological properties, Agg backend
    test_smoke.py          # every plot_* / session delegate runs, returns expected type
    test_cartoon.py        # one node per element, brackets = closedness, ordinal ranks, walks element-to-element, nested keepout
    test_dual_graph_plot.py
```

**Fixture strategy.**
- One builder per case. Every `ALL_CASES` member (k10, k28_1, k28_2, p3, nested, inversion) is session-scoped and read-only. A teardown guard asserts that the generation has not changed.
- Tests that mutate request `fresh_case("k10")`.
- Synthetic cases (stable_line, hand-built lobes, FakeDual) live in `cases.py`.
- `build_pieces` is removed in favour of the session API.

**Markers.**
- `golden`: run by default, and `--update-golden` re-records snapshot files.
- `perf`: opt-in.
- `slow` is dropped, or defined as more than 0.5 s, depending on Q38.
- Plain `regression/` directory, no marker.
- Hygiene rules move to ruff / import-linter (Q22).

**Estimated test count after the refactor.**

| Tier | Functions | Collected cases |
|---|---|---|
| Invariants (~35 checks × ~6 cases) | 35 | ~200 |
| Unit contracts | ~230 | ~300 |
| Facade / cache (table-driven) | ~15 | ~50 |
| Golden | ~12 | ~15 |
| Regression | ~20 | ~25 |
| Plotting | ~25 | ~35 |
| **Total** | **~335** | **~600–630** |

That is about 25% fewer collected cases. Invariant coverage across fixtures goes **up** (each invariant moves from one or two fixtures to six). There are about 45% fewer test functions and an estimated 25–30 s of wall time instead of about 45 s, and much less brittleness: under the default answers (deleting categories A, B, F, G and H, and loosening J and K) about 180 cases that break on harmless changes go away.

**Plan sequencing, once the questions are answered:**
1. Infrastructure: builders, frozen fixtures, helper fixes, delete `visualizations/`.
2. Delete tombstones and migration pins (Q21, Q23).
3. Build the invariant tier and delete the per-fixture twins it subsumes.
4. Collapse the facade and cache tests into tables.
5. Consolidate golden facts (Q1–Q6).
6. Loosen strings and logs (Q34–Q37).
7. Plotting tier (Q14–Q15).
8. Write the gap tests from §5.

Run the suite after each step, and check that coverage (`pytest --cov`) does not drop for any module.