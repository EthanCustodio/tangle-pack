# testing — knowledge file

Last updated 2026-10-02 (auditor second pass, branch `refactor/hole-side-and-own-blast` at `f8a6d73`).

- Suite at `2315204`: 782 passed, 1 skipped (CuPy); on the branch 795 passed, 1 skipped; ~50 s, `MPLBACKEND=Agg env/bin/python -m
  pytest -q`. No pytest-xdist in env.
- Higher-period fixtures: `tests/test_higher_period_cartoon.py` (`p3_built` = `build_period3()`,
  `nested_built` = `build_nested()` = outer 2 / inner 0), all `@pytest.mark.slow` but fast (~2 s).
  Since the branch: p3 at 4 blasts (words), p3 at 6 blasts and 15 steps (runs, I1 holds, not
  reliable), nested inner words == p3 alone (2 outer + 4 inner), nested 4+4 order independence
  and no "skipping bridge" warning, pinned hole sides (p3, k28 two blasts), direct-hole side ==
  old measured side (4 fixtures). Compare across sessions through element short names + orbit
  code + direction (`_tangle_words`), never letters; select a tangle's classes by
  `bridge_class.source.branch_key[0] is fp`, not `entry.tangle` (session order differs from
  `TangleBuild.fixed_points`).
- Registry ids were NOT reproducible between builds in one process until WP8 (15-step failing
  origin seen as (12,33), (21,35), (3,17)); they are now, but tests still pin iterates and
  sides, never ids.
- Worktrees have no env: run `../tangle-pack/env/bin/python -m pytest` inside the worktree;
  `pyproject.toml` `pythonpath = ["src", "tests"]` makes it import that worktree's src.
- Probe technique that works without touching the repo: put a module on PYTHONPATH that
  monkeypatches (e.g. `topology.StablePartition._punch_in_bridge`, `TangleSession.blast_zone`)
  and load it with `pytest -p <module>`; `Trellis.punch_holes` imports `punch_holes` from the
  module at call time, so patching the module attribute takes effect.
- Every henon_cases build is < 1 s without blasts; p1 at 8 blasts ~4 s.
- Proof scripts (branch only, commit `f8a6d73`): `scripts/refactor_proof/hole-side-and-own-blast/`
  `words_grid.py` (nested outer x inner grid vs each tangle alone, per tree, `--tree/--label`,
  `--blast-both` for base's both-saddles recipe, `--plot`), `hole_sides.py` (wraps
  `_punch_in_bridge`, reads each propagated hole's side old vs new), `panels.py` (six
  `draw_panels` per case). Each inserts its tree's `src` and asserts the import path: the env's
  editable install points at the MAIN tree. NOTE `words_grid.py` runs p3-alone at
  `min_separation=1e-4` (nested's default), `build_period3` defaults to 1e-5; the 0..5 stable /
  6+ jump to 42 holes holds at both.
- Grid fact (base `2315204`, recipe both saddles): 12 of 20 nested cells raise (inner>=2 & outer>=4,
  (2,4),(2,6),(0,6)) plus outer-first (4,4); (2,2) order-dependent; (6,0),(8,0) silently corrupt
  period 3 (a,b unresolved, `c -> b^-1`). Branch: all 20 match each tangle alone; order checks at
  (2,2),(4,4),(8,4) identical (itinerary tables identical apart from the title).

## Numerics stress campaign (branch `numerics-stress`, `stress/numerics/`, audited 2026-10-07)

- Map: `harness.py` (subprocess per config, RSS watchdog, outcome ok/exception/invariant_assert/
  oom/timeout/crash), `e01..e10_*.py`, `results/*.jsonl`, `report/` (`report.tex`,
  `make_figures.py`, `analysis.py`, `make_numbers*.py` -> `numbers*.tex` macros, `tables/`).
- E7 `blasts` (`e07_breaking.py:261`): ONE config per case (k28, p3 max 10; nested 8), one cutoff,
  fixed `min_separation`; topology probed ONLY on the final state (`run_topology`, :131) ->
  no per-blast series of classes/words. Final: k28 32 classes / 7 unresolved / 2 virtual;
  p3 31/1/1; nested 13/0/0. E7 `depth`: p3 classes 9 (11-14 steps) -> 23, 30, 36 (15-17),
  is_reliable False from 15; inversion 4,5 = genuine topology AssertionError ("owned by 0
  partition elements"), 6 = ok but is_reliable False (1 virtual), 7,8 OOM.
- SECOND EDITION of the report (WP10, merge 70222ee): E3, E7, E8, E9 and E11 rerun in full on
  the fixed library (harness commits efccdd4 / 61ce364); E1, E2, E4, E5, E6, E10 keep their
  first-campaign numbers, labelled as such. Harness changes: `_worker` records
  `PointBudgetExceeded` details (outcome stays `exception`), `maps.resolution_metrics`
  summarises `workbench.report_resolution`, E9 times `check_resolution` separately
  (`compute_intersections(diagnose=False)`), E7 depth cells pin `area_cutoff=1e-4`,
  `session.area_cutoff` setter replaces the private poke. Rerun: 0 OOM kills, 9 budget stops;
  E8 0 lower-period accepted; E9 154/154 matched; E11 error 1 gone, 43 error snapshots in 6
  configs (39 bridge rows = D1 resolution effect, 3 hole side, 1 outer face).
- OPEN (author decision): the blast-end same-stability pass `blast_zone` runs by default keeps
  1e5-1e6 hits per blast from blast 7-8, costs 1.6-2x the blast from blast 10-11 and cut E11's
  reachable depth by 3-4 blasts, so E11's main grid blasts with `diagnose=False` and the pass is
  measured only in five `resolution_pass` cells (report sec. 10.6). Bound it (count hits, check
  only new image bridges) or default `diagnose=False` for blasting.
- AUTHOR RULE 2026-10-07 for any write-up: topology changing under blasting / deeper growth
  (new crossings, classes, holes, splits, words, letters) is EXPECTED, never a failure.
  Three tiers: (a) topology change = new information; (b) library provisional flags
  (`is_reliable` False, ambiguous, unresolved, virtual `newN`) = symbolic stage incomplete,
  open question; (c) genuine errors = exception, AssertionError, OOM/timeout, physical-law
  violation. `is_reliable` (`SymbolicDynamics.py:826`) = every class has an itinerary, none
  ambiguous, none contradicted by registered evidence.
- Report framing smells (as of 2026-10-07): scorecard Blasting row "topology unreliable",
  taxonomy row "Unreliable topology", `make_numbers.py:428` "(unreliable)" macro suffix,
  `make_figures.py:816,888` legend strings "unreliable"/"topology unreliable", section
  "Failure modes and breaking points" holding blasts + depth.

## Stress-report recommendations integrated (2026-10-07, WP1-WP10, main 12d5600)

- Ten packages from `~/.claude/plans/great-okay-so-there-purrfect-turtle.md` ran as coder +
  reviewer pairs on their own branches (`wp1-…` … `wp10-…`, worktrees under
  `.claude/worktrees/wf_20fd0a4d-e33-*`, all KEPT) and were merged `--no-ff` onto main in the
  plan's order, the full suite green after each: WP1 faf3940, WP2 cd22646, WP3 3ad5195,
  WP4-A 1b3e0f4, WP7 958e552, WP8 a6b782d, WP9 0c9e702, WP5 ded56db, WP6 41e5c7d,
  WP10 70222ee. Nothing pushed. Suite 738 -> 833 passed, 3 perf deselected (pass with
  `-m perf`), 11 xfail = exactly `cases.KNOWN_ISSUES`, no XPASS, ~100 s (was ~85 s; the WP4
  deep-blast regression costs ~7.5 s).
- What the tests gained: unit `test_crossing_resolution`, `test_same_stability_diagnostic`
  (+1 perf), `test_fixed_point_minimal_period`, `test_point_budget` (+1 perf),
  `test_exact_ownership`, `test_api_traps`, `test_blast_budget`, production-check tests for
  `RowDiagnosis`/`SideDiagnosis`/`OuterFaceError`; regression
  `test_resolution_warning_follows_oracle` (WARNING iff the independent oracle finds a hit),
  `test_deep_blast_sliver_beyond_pip` (k28, 9 blasts, ANALYTIC seed; asserts its own
  precondition loudly); `case_sanity` laws `check_fixed_point_minimal_period`,
  `check_map_is_area_preserving`, `check_case_area_cutoff`; `dual_graph` law
  `check_stable_nodes_never_own_a_singleton`. Builders: `near_tangency` (k=-0.999, 1e-10) and
  `k28_deep_blast_sliver` in BUILDERS, not LAW_CASES.
- Rule extension (WP6): builders never read LIBRARY defaults either (`CaseExpect.area_cutoff`
  is required; `build_k10`/`build_inversion`/`build_orientation_reversing` freeze 1e-4, k28 1e-7,
  near_tangency 1e-10). `DEFAULT_AREA_CUTOFF` is 1e-7 now.
- Review process lessons: pyproject sets `pythonpath=["src","tests"]`, so break-it mutants must
  be run on a COPIED tree (src+tests+pyproject), never via PYTHONPATH; `coverage_guard.py`
  compares raw line numbers and trips on every edited covered file (a content-mapped diff found
  0 newly missed lines) -> rebase BASE.json after each merge or make it content-mapped.
- Deferred / open after the round: WP4-B anchor-only stable branch (inversion 4-5 steps, see
  topology.md); WP7-D1 minimum nodes per bridge (39 remaining E11 bridge-row snapshots);
  the parity flip on carried hole sides is untested end to end (no det J < 0 case builds);
  the blast-end same-stability pass does not scale past ~10 blasts (author decision, see the
  campaign section above); `MaxIterationsReached`/`TURNAROUND_LOG_STRETCH` not exported from
  the package `__init__`s; `PointBudgetExceeded.blast_context` set dynamically, undeclared on
  the class; `max_points` on `blast_zone` is a PER-CALL budget (author to confirm);
  `scripts/*.ipynb` now refine at 1e-7; committed `figures/.../nested_p1_p3/partitions.png`
  was already stale at 8148d48 (crossing-id labels only).

## Open-issues round WP1 (2026-10-07, branch `oi/wp01-decision-laws`, base `d76ef7c`)

- Tests only; `src/` byte-identical. Law cases now `k10, k28_one_blast, k28_two_blasts, p3,
  nested, nested_inner_blasts, inversion, inversion_one_blast, inversion_two_blasts`.
  `inversion*` = the k=4 OTHER root `(1-sqrt5, sqrt5-1)` (eigenvalues -1.96/-0.51),
  `area_cutoff=1e-5`, 11/11 steps (12 escapes), zone at `f(q0)` via the shared k28 recipe
  `cases._zone_at_image_pip`, 0/1/2 blasts, each build ~0.1-0.2 s. Holes 1 / 7 / 0; the 0 at two
  blasts is INV-2 (patching `Trellis.lambda_u` to beta^k gives 1 / 9 / 9). The old k=10 saddle is
  `build_k10_inversion`, only behind the numerics fixtures `henon_inversion*` (moving those to
  k=4 is WP3). `nested_inner_blasts` = `build_nested(inner_blasts=2)`; other blast counts are
  named `nested_<outer>_<inner>`.
- New laws: symbolic `cross_side_only_within_a_branch` (D4; replaced `itinerary_pairs_same_side`),
  `every_active_class_walks` (D3; walk status unique/trivial unless a landing is a singleton;
  not applicable on p3 and k28_one_blast, whose active classes all land on singletons),
  `no_virtual_class` (D3); classes `alphabet_inside_minimal_trellis` (D2).
- Blast series (`cases.SERIES`, `blast_series`, `tests/regression/test_class_count_under_blasting.py`):
  the D1 loop on one session from the zone stage (k28 8 blasts, p3 9, inversion 2;
  `min_separation=1e-5`, `diagnose=False`), ~5 s total, cached per run in immutable
  `SeriesStep`s. Counts today: k28 classes 3,3,4,4,4,4,6,7,23 (holes switch at 2 and 8); p3
  9 x6, 23, 25, 27, 29 (switch at 6); inversion 6, 11, 6 (holes 1, 7, 0). The law compares the
  TOTAL class count, not per kind: per kind moves legitimately (k28 0->1 active 3->1 as image
  loops register, inert drops at reference switches). Author may overrule (D1 said "the class
  count"); the per-kind counts are recorded, so it is a one-line change.
- Set-aside tests: `helpers.law_tier.layer_tests(layer) -> (test_<layer>_laws,
  test_<layer>_set_aside)`; the set-aside test is a strict xfail for a known issue and, for a
  not-applicable pair, runs the law, asserts count 0 and skips with
  "not applicable: <law> on <case>: <reason>" (19 skips). `KnownIssue.wp` names the fixing WP
  (or `C3`), reasons lead with the catalog id; the meta test checks `WP\d+|C\d`, KI/NA disjoint,
  law-tier pairs on LAW_CASES.
- KNOWN_ISSUES 11 -> 52 (6 kept, re-tagged; 5 k=10 inversion entries removed; 46 new naming
  WP2/WP3/WP6/WP9/WP10, none WP1). NOT_APPLICABLE 28 -> 19. WP2 must drop the 15 vacuous
  `inversion_two_blasts` entries when holes appear and ADD `stable_nodes_never_own_a_singleton`
  there (WP3).
- TEST-5: `generation_bump` cache cells run on `build_k28(blasts=1)` + one blast and assert the
  direct view changed; the p3 near-ULP stand-in is accepted (constructor `area_cutoff`, no
  `_man_machine` poke); three tangles -> strict xfail `three_tangles_take_three_families` (WP9).
- Break-it (copied tree): `cross_side = True` reddens symbolic on p3/nested/inversion_one_blast;
  DualWalk status "unreachable" for one walk reddens `every_active_class_walks` on k10; a class
  decrease after blast 1 (NOT on the last step: k28 jumps 7 -> 23 there, so -1 is invisible)
  reddens `never_decreases`; a false NA pair fails its set-aside test with "now checks".
- Inversion x period-k: no fixture builds today (period-2 orbits with negative DM^2 eigenvalues
  at k >= 4.1 cross only at their 8 anchors at 3-5 steps -> INV-3; k=4.5 escapes ~6 steps); WP3.
- Suite after WP1: 871 passed, 19 skipped, 3 deselected, 52 xfailed, ~113 s (base on the same
  machine 833 passed / 11 xfailed, ~106 s). The k=10 inversion law builds (1.3 s per module)
  are gone; `nested_inner_blasts` (~0.6 s per module), the three k=4 cases, the series (~5 s)
  and the k28 generation_bump cells more than offset that.

## Open-issues round WP2 (2026-10-08, branch `oi/wp02-exact-iterates`, base `f39f617`)

- Suite: 877 passed, 19 skipped, 3 deselected, 34 xfailed, ~106 s (base on the same machine 871
  passed / 52 xfailed, 108 s). KNOWN_ISSUES 52 -> 34: the 15 `_HOLE_LAWS` x `inversion_two_blasts`
  and the two `inversion_series` entries (all WP2) XPASSed with the beta**k lap and went, as did
  `anchor_bridge_class_leads_its_tangle` and `cross_side_only_within_a_branch` on
  `inversion_two_blasts` (WP3 tags; still xfail on `inversion`). Added:
  `stable_nodes_never_own_a_singleton` x `inversion_two_blasts` (INV-1, WP3), folded with
  `partition_singletons` (no longer vacuous) into the `inversion_one_blast` entry. No WP2 entry
  remains; no remaining entry is vacuous for want of holes.
- New crossing-layer laws (all nine LAW_CASES, no set-aside): `iterate_table_partial_bijection`,
  `iterate_links_pass_the_map_gate` (the real map as the oracle), `iterate_table_complete` (the
  law recomputes both brackets from the nodes). `iterate_link_scales_by_beta` and
  `area_along_iterate_links` now take `SCALING_RTOL` / `AREA_LINK_RTOL` plus the linked
  crossings' own bracket widths as slack (`LINK_SCALING_RTOL`, read off the deleted
  `cdist_rtol`, is gone).
- New tests: unit `test_iterate_inference.py` (one `map` call per pass from an empty table,
  rebuilding the same links; the gate refuses a displaced image with one WARNING), registry
  `test_dedup_is_relative` / `test_iterate_table_is_a_partial_bijection`, pseudoneighbor
  `test_inversion_branch_return_scales_by_beta_to_the_k` / `test_reference_end_is_a_table_iterate`,
  regression `test_iterates_after_blasts.py` (p3 at ten blasts, ~1 s). All fail on the base.
- Adapted, not weakened: `test_grow_until` stand-in links now point at keyless stand-in crossings
  (a self-link gave a crossing two preimages, now asserted); `test_point_budget`'s image-cost spy
  counts the machine's evaluations only (the gate's one map call creates no points); the
  pseudoneighbor residue / single-window / cut-warning tests register the table links they used
  to get from the cdist fallback. Deleted: `test_unknown_unstable_branch_is_conservative`,
  `test_reference_end_fallback_by_cdist_match`, StrongPip `missing_unstable_key_blocks` (PK-7 and
  the removed fallback), `test_scale_cdist_is_in_map_steps` (API gone).
- Series after WP2: inversion holes [1, 9, 9], classes [6, 11, 11]; p3 classes
  9 x6, 23, 25, 27, 28 (last step 29 -> 28); k28 unchanged. p3 at ten blasts: 62 detected
  forward links (+3 anchor links), 29 classes (31 on base).

## Open-issues round WP3 (2026-10-08, branch `oi/wp03-inversion`, base `c16e75d`)

- Suite: 905 passed, 25 skipped, 3 deselected, 21 xfailed, ~103 s (base on the same machine
  877 passed / 19 skipped / 34 xfailed, 104 s). KNOWN_ISSUES 34 -> 21, no WP3 entry left: 13
  pairs XPASSed and were removed; `anchor_bridge_class_leads_its_tangle` x `inversion` re-tagged
  WP5 (PK-3) and `arrangement_euler` x `inversion_two_blasts` WP7 (INT-9).
- New law case `inversion_period2` (`build_inversion_period2`: k=4.1, `HENON_K41`, frozen seed
  `K41_INVERSION_PERIOD2_SEED`, 1e-5, 9/9 steps, no zone, default strong pip); in LAW_CASES and
  BUILDERS; NOT_APPLICABLE `iterated_cut_provenance`, `refined_children_inherit_word`,
  `member_matching_consistent`. Main fails 7 laws on it. At 10 steps it fails 11 laws on the
  branch (all failing on main too): deep period-k, PK-1 evidence.
- Laws: `one_anchor_per_unstable_branch` REPLACED by `one_anchor_per_periodic_point` (one anchor
  per point, first on every trellis branch); new `anchor_rays_are_its_branches` (anchors layer)
  and `anchor_opening_faces_its_hole` (partition layer; NA on k10, k28_two_blasts, inversion).
  Adapted, not weakened: `bridge_id_in_unstable_order` (no two-anchor exemption),
  `bridge_endpoints_on_own_branch` (an anchor end is the bridge's own point's anchor, still in
  span), `_element_ends` / `bridge_rows_consistent` / iterated laws read an anchor end through
  the independent `_stable_branch_at` (C7), opening laws check the record's branch
  (`_opening_row`, `_anchor_row` from the home sign), `partition_covers_branch` matches row AND
  branch.
- New tests (all fail on base except the INV-7 characterization, which is behaviour-preserving):
  regression `test_inversion_before_the_tangle_forms` (C1), `test_inversion_zone` (INV-5);
  unit `test_a_periodic_point_is_one_anchor_of_all_its_branches`,
  `test_iterating_an_anchor_rooted_inversion_bridge_cuts_identified_children` (INV-10),
  `test_growing_an_inversion_point_logs_no_branch_index_warning` (INV-9),
  `test_orient_eigenvectors_acts_on_an_inversion_point` (INV-8),
  `test_kevin_way_seeds_one_chain_per_eigendirection` (INV-7, passes on base),
  `ring_cycles` tests in `test_resonance_zone.py` (import fails on base),
  `test_cross_side_is_flagged_only_within_one_branch` (INV-6),
  `test_the_inversion_anchor_is_one_node_with_four_real_rays`.
- Fakes borrow the new rules (`FakeTrellis` / `_ClassTrellis` / iterated-cut `_FakeTrellis`:
  `unstable_key_of`, `stable_key_at`, `is_anchor=False`; `_OwnerFamily.owner_at_end`); hand-built
  openings are 4-tuples naming the crossing's branch; `test_generation_and_caches` fake keys are
  4-tuples. Deleted: `test_num_branches_agrees_with_get_branch_array` (API gone).
- WP1 deferred (b): the `henon_inversion*` numerics fixtures STAY on k=10 (no tangle needed, or
  bridges at 6 steps, which k=4 only reaches at 10-11; the k=4 saddle is already the law case).

## Open-issues round WP4 (2026-10-08, branch `oi/wp04-sign-rule`, base `dad72d2`)

- Suite: 910 passed, 25 skipped, 3 deselected, 21 xfailed, ~103 s (base 905 passed / 25
  skipped / 21 xfailed). KNOWN_ISSUES stays 21 (no entry was tagged WP4); the
  `open_p3_deep_runs` reason now names only the virtual `new1`.
- Law `openings_follow_signs` (partition layer) REPLACES `openings_missing_only_at_anchor_or_tail`:
  exactly one opening per non-anchor bound, at most one at an anchor bound; a non-anchor
  opening is on the crossing's own branch, on the row the sign gives (left iff
  `(sign > 0) == (end == first)`, written in the law, not via `row_of_end`) and `outward` iff
  `(sign > 0) == (bridge_side == "left")`. Fails on the base on `inversion_period2` (tail
  openings the geometry dropped), passes on all 10 law cases. No NOT_APPLICABLE.
- Builder `p3_coarse_blasts` (`cases.build_period3_coarse_blasts`, BUILDERS only, not a law
  case): the p3 zone at 1e-7, then 6 single blasts at `area_cutoff = 1e-5`,
  `min_separation = 1e-5`, `diagnose=False`, repinned after each, partitioned (the E11 cell).
  Regression `test_coarse_blasts_rows_by_signs`: precondition (a hole bridge whose signs
  alternate while `geometry_disagrees`), the symbolic dynamics runs, laws
  `partition_unique_owner` / `openings_follow_signs` / `bridge_rows_consistent`, and
  `punch_holes()` logs more StablePartition WARNINGs than `punch_holes(diagnose=False)`
  (level + logger only). Base: AssertionError "approaches its two defining crossings" in the
  builder.
- `test_open_p3_deep_runs`: `test_deep_p3_every_class_walks` (no `unresolved_reason`, D3, no
  xfail) split off the strict xfail, which keeps only `not dyn.virtual_classes`. Fails on the
  base (the unreachable class) on all three deep runs.
- Unit, rewritten to the new rule (not weakened): the I2 stubs carry a BridgeId and the TRUE
  crossing signs of their lobe (`_LOBE_SIGNS = (-1, 1)`); equal signs raise, alternating signs
  with `_CROSSING_LOBE` pass and come back from `geometry_flags`
  (`test_contradicting_geometry_with_alternating_signs_is_a_flag`, NUM-1 in miniature);
  `test_measured_side_disagreement_is_a_flag_not_an_error` replaces `..._raises`;
  `bridge_side_violations` counts the orbit check only. Production checks:
  `test_punch_holes_flags_a_contradicting_measured_side` (one WARNING, no raise) replaces
  `..._raises_...`, new `test_punch_holes_raises_on_equal_signs`. Stable partition:
  `test_containing_bridge_brackets_by_bisect` replaces the key-filter test,
  `test_side_at_nearest_segment_reads_two_nodes` (INT-5) is new; deleted with their API:
  `test_bridge_span_reads_the_branch_from_the_bridges_own_key` (covered by law
  `propagated_holes_land_on_predicted_branch`), `test_backward_endpoint_prefers_the_table_then_scales`
  (`image_cdist`'s table-first rule is pinned in `test_image_cdist.py`).

## Open-issues round WP5 (2026-10-08, branch `oi/wp05-period-k-holes`, base `3e12d53`)

- Suite: 921 passed, 25 skipped, 3 deselected, 21 xfailed, ~104 s (base 910 passed / 25
  skipped / 21 xfailed). KNOWN_ISSUES stays 21: `("anchor_bridge_class_leads_its_tangle",
  "inversion")` keeps `wp="WP5"` with its reason re-tagged (PK-3 / DEC-19: the `+1` image lies
  past the computed u(0.0), a stagger; its hole cannot be made without inventing one). No entry
  added or removed.
- Law `window_holes_on_every_branch` (partition layer, after `no_direct_hole_beyond_fundamental`;
  D5's positive half): for every reference and `j < k_value` whose `+j` image pair is registered
  and spanned, a direct (`side_source == "signs"`) hole with that origin, iterate and bounding
  pair exists, on `advance_key(branch, j)`. Coverage that was missing, not a fix: passes on base
  and branch on all 10 law cases (counts 1/2/1/6/7/7/1/3/3/2); `j = 0` keeps it non-vacuous, no
  NOT_APPLICABLE. `test_law_partition.py` and the `no_direct_hole_beyond_fundamental` docstring
  no longer call the `+1..+(k-1)` rule provisional; `conftest.py` Dev Notes drop it from
  "Provisional rules".
- Regression `test_period_k_window_gaps.py` (PK-3): gaps by rule from the table (first `j` with
  an unlinked member -> branches `advance_key(branch, j..k-1)`); on `inversion`,
  `inversion_period2`, `p3_6_blasts` a rerun of `trellis.compute_pseudoneighbors()` logs exactly
  one Pseudoneighbor WARNING per gap branch and no direct hole of a gapped reference sits on its
  gap branch; on `p3`, `k28_one_blast` (no gap, asserted) it logs none. Base: 0 records, the
  gap test fails.
- Unit, rewritten (FAIL on base, pass on branch): `test_missing_pip_iterates_are_cut_at_their_predicted_cdist`
  (PK-6; replaces `test_strong_pip_cuts_warn_when_branches_are_left_uncut`; a crossing on
  `s(1,0)` past the predicted cut 2.0 would puncture the first window pair if uncut),
  `test_missing_forward_appearance_warns_once_per_branch` (PK-3; 1 / 2 / 0 WARNINGs for a
  period-3 `+2` gap, a period-3 unlinked reference, a period-1 reference),
  `test_inert_class_mapping_onto_inert_classes_is_silent` (PK-8; replaces
  `test_inert_class_with_word_warns_and_stays_out`). New, passes on both:
  `test_inert_class_naming_an_active_class_warns_and_stays_out`. Stable partition: the
  forward-pair kernel takes `(pair)` and gains the firm rows `(3, 1, False)`, `(3, 2, False)`.
- `p3_6_blasts` / `inversion_period2` builds take under a second each here; the regression adds
  ~2 s.

## Open-issues round WP6 (2026-10-08, branch `oi/wp06-minimal-alphabet`, base `b90b718`)

- Suite: 938 passed, 27 skipped, 3 deselected, 13 xfailed, ~108 s (base 921 passed / 25
  skipped / 21 xfailed). KNOWN_ISSUES 21 -> 13: all eight WP6 entries XPASSed and were removed
  (`alphabet_inside_minimal_trellis` x inversion / inversion_one_blast / inversion_two_blasts,
  `series_class_count_constant_while_holes_constant` x k28_series / p3_series,
  `series_alphabet_inside_minimal_trellis` x the three series), with `_BCG1` and the now unused
  `_INVERSION` tuple. NOT_APPLICABLE +2: `every_active_class_walks` and `matrix_is_token_counts`
  on `inversion` (its alphabet is one inert class until DEC-19); they report skips and flip
  loudly when s(0.1) gets its hole.
- Laws (classes layer). `_evidence(case)` = `minimal_trellis().classes`. The classing laws read
  it (`every_bridge_in_one_class`, `class_orientation_anchor_outward`,
  `anchor_bridge_class_leads_its_tangle`, `class_table_order`, `classes_do_not_mix_tangles`,
  `classes_use_every_orbit_branch`; `minimal_trellis_bridges` reads `minimal.classes`): on base
  they read this same full table, so no weakening. `every_bridge_in_one_class` counts
  `table.orphans` (each a loop by independent `_element_ends`); `class_orientation_anchor_outward`
  asserts no class is a loop (was "or inert"). NEW `loops_fold_forward` (D2): each loop's
  `folded_from` is a non-loop member of its class, reached from the loop up independently indexed
  `image_chain` parents through loops only; returns entries + loops, non-vacuous everywhere.
  `alphabet_inside_minimal_trellis` is an iff over the evidence table and an evidence class
  outside the alphabet is unlettered. The lettering law `only_active_classes_lettered` reads the
  alphabet.
- Unit (`test_bridge_class.py`): `test_loop_folds_into_the_class_whose_image_holds_it` (was the
  preimage test), `test_orphan_loop_is_logged_and_unclassed` (was "stands alone and warns"),
  `test_restricted_to_keeps_the_entries_with_a_member_in_table_order`, the hand-built BCG-2
  kernel `test_a_fold_born_loop_folds_forward` and `test_an_image_pair_has_one_parent` (the
  bijection assertion); `_ClassTrellis` gains `bridges`, `iterates` and `branch()`. Deleted
  `test_an_unresolved_loop_class_is_inert_outright` (the rule is gone).
  `test_symbolic_dynamics.py`: `test_letter_series_are_disjoint` replaces
  `test_inert_letters_skip_active_letters` (the skipping is gone; `test_inert_letter_series`
  kept); new `test_registered_evidence_never_extends_the_letters`. Facade:
  `test_letter_sequence` deleted (covered by the disjoint-series test); `helpers/direct.py` has
  `evidence()`, `minimal()` over it and `classes()` = `minimal().alphabet`; equivalence adds
  `inversion_one_blast` (alphabet strictly inside the evidence table).
- Regression `test_blast_children_are_evidence_only.py` (k28 zone, 6 single blasts, ~3 s): the
  alphabet is exactly the evidence classes with a kept member (and strictly smaller there), no
  virtual class, every word symbol in the alphabet, matrix nodes = the active alphabet's refined
  symbols. Fails on base.
- `test_class_count_under_blasting.py` unchanged in code: 9 passed, 0 xfailed.

## Open-issues round WP7 (2026-10-08, branch `oi/wp07-walk-exists`, base `e883720`)

- Suite: 952 passed, 28 skipped, 3 deselected, 2 xfailed, ~104 s after the review fix (base 938 passed / 27 skipped
  / 13 xfailed). KNOWN_ISSUES 13 -> 2: removed `arrangement_euler` x `inversion_two_blasts`
  (WP7), `_NUM16` with its six entries (`arrangement_regions_disjoint` x k28_one_blast /
  k28_two_blasts / inversion_one_blast / inversion_two_blasts / nested_inner_blasts,
  `arrangement_regions_sound` x nested_inner_blasts: root cause INT-9), `open_p3_deep_runs` x 3,
  `anchor_bridge_class_leads_its_tangle` x `inversion`. Left: `orientation_reversing_case_builds`
  (C3), `three_tangles_take_three_families` (WP9). NOT_APPLICABLE +1:
  `anchor_bridge_class_leads_its_tangle` on `inversion` (D9: s(0.1) has no hole).
- Laws: `alphabet_closed_under_words` replaces `no_virtual_class` (every rule symbol in
  `dyn.table`, `dyn.classes` in table order, seed iff `named_by is None`, a closure class is
  in `minimal_trellis().classes` and named in its `named_by`'s word, and any class without an
  itinerary makes `is_reliable` False; returns symbols + classes).
  `matrix_is_token_counts` lost its virtual-sink branch (every active token is a node). NEW
  `dual_unified_is_where_f_of_s_stops` (dual-graph layer; helper
  `helpers.laws.misplaced_unified_crossings` reads `S \ f(S)` off the iterate table: a crossing
  of a pip-owning branch is unified iff it has no registered preimage -- independent of
  `DualGraph._check_where_f_of_s_stops`' cdist arithmetic, review fix). `anchor_bridge_class_leads_its_tangle`
  checks only tangles whose every stable branch carries a hole (`Hole.openings` branch keys,
  D9 / DEC-19), the rule the probe found (scoping by the anchor's own branch does not work:
  s(0.0) has the window hole and u(0.0)'s anchor class is still inert). All pass on the 10 law
  cases; the anchor law is non-vacuous on the 9 others.
- Regression: `test_open_p3_deep_runs.py` has no xfail (`test_deep_p3_alphabet_is_closed`
  x3: every rule symbol in the table, every closure class spanned by a registered bridge, and
  every class whose kept walk names a pair outside table + evidence (computed independently in
  the test) is unresolved, warned, and the dynamics unreliable; `test_deep_p3_every_class_walks`:
  every class has an itinerary or a kept walk; the sink test is deleted);
  NEW `test_closure_predicts_registered_images.py` (p3_16: a closure class with a registered
  member chain exists and every such class is verified), `test_evidence_diagnosis.py` (k28 at 6
  blasts, ~3 s: each disagreeing `ImageEvidence` has a diagnosis whose cause matches an
  independent equal-sign scan of its chain; agreeing ones have none; unreliable),
  `test_unified_segment_precondition.py` (`build_inversion_period2(unstable_steps=10,
  stable_steps=10)` under `np.errstate(over="ignore")`: one DualGraph WARNING per independently
  differing branch, naming it; p3 / k10 silent). `test_blast_children_are_evidence_only.py`
  reads `dyn.table` and checks seed iff `named_by is None`. `helpers/names.py` no longer looks
  up virtual names.
- Unit: `test_symbolic_dynamics.py` — `test_a_class_outside_the_table_is_named_by_its_label`
  (was the virtual `new1` test), `test_a_straddling_landing_is_unreliable`, and a monkeypatched
  driver harness (`_drive`: `ElementNaming`, `land_element`, `_resolve_itinerary`,
  `_member_pair`, `_image_chain`, `trellis_itinerary`, a `_Crossings` trellis stand-in) for
  `test_a_word_naming_a_new_class_closes_the_alphabet` (active via `letter_for`; inert from an
  evidence entry) and `test_evidence_disagreement_names_a_stable_end` (both causes).
  `test_arrangement.py` — `test_is_computed_arc_truth_table`,
  `test_a_gap_between_two_curves_is_two_dangling_ends`, and `test_arrangement_survives_a_blast`
  rewritten (a gap exists after the k10 blast, it is no arc, every drawn arc is computed, Euler
  per component).
- `build_inversion_period2(*, unstable_steps=9, stable_steps=9, through=...)` (law case
  unchanged).
- Base check: the branch tests copied onto a detached base worktree
  (`../tp-oi-wp07-walk-exists-base`) fail 31 there (every new test above, the arrangement laws
  on the five blasted cases, the symbolic laws on all ten).

## Open-issues round WP8 (2026-10-08, branch `oi/wp08-recompute-identity`, base `083f910`)

- Suite: 974 passed, 28 skipped, 3 deselected, 2 xfailed, ~120 s after the second review fixes (972
  after the first, 969 before them; base 952 passed / 28 skipped / 2 xfailed). KNOWN_ISSUES 2 -> 2 (C3, WP9); NOT_APPLICABLE unchanged.
- Registry ids are reproducible (INT-7): `test_registry_ids_are_reproducible` (two
  `build_k28(blasts=2)` in one process, same `id -> (keys, cdists)` map). Tests still never pin ids.
- Builders: `cases._repin` / `_partition` gone; `_blast(session, zone, fixed_points, *, diagnose)`
  takes no pips (the session carries them); the `partitioned` stage is `session.topology()`.
  `SeriesStep` gains `references`, `unwitnessed` (via `helpers.laws.reference_witness`: a crossing
  strictly between the pair, or between its registered images, on either branch; table only) and
  `letters` (`(identity, letter)`, the identity built IN `cases._class_identity` from
  `homotopy_partition().element(ref)` lo/hi ids, never `Topology.identity`).
- New laws (series): `series_lost_references_are_witnessed` (passes on base too: a
  characterization of BCG-4), `series_letters_follow_element_identity` (fails on base: 3 classes
  re-lettered on `p3_series` from blast 6).
- New: unit `test_unlinked_landing_is_read_against_the_pairs_images` x2 (frame fix; `below_c`
  fails on base), regression `test_deep_reference_survives_noise` (nested zone stage, 9 outer
  blasts, recompute, 1 more, recompute; `diagnose=False` blasts, ~11 s, kept in the default run;
  base: `(37,38)` lost with no witness, 13 -> 9 classes), facade `test_session_topology.py`
  (every product is a stage of one `topology()`; a blast / `set_strong_pip` recomputes; the pip
  survives a blast; `trellis(fp) is trellis([fp])`; gathering follows the selection, asked before
  or after the full topology).
- Adapted, not weakened: `helpers.direct.DirectBuild` builds its OWN per-fixed-point trellises
  for the SELECTED fixed points (classify, the session's pip, pseudoneighbors, holes, partition)
  and letters by element bounds on a fresh alphabet. `test_session_caches.py`:
  `generation_bump` blasts with no repin or manual partition; `partition_signature` ->
  `manual_repartition` (the per-fp trellis's own partition changes, no product follows);
  `pip_change` through `session.set_strong_pip`; the letters test reads
  `bridge_alphabet.assigned[identity]` after a pip change. `test_session_bridge_classes.py`:
  `test_no_partitions_raises_and_warns` -> `test_products_need_only_the_numerics` (D1: no manual
  stage call) + `test_a_pip_trimmed_away_warns_and_raises`. `test_unified_segment_precondition`
  wraps its `dual_graph(rebuild=True)` in `np.errstate(over="ignore")` (a rebuild recomputes the
  holes, whose geometric check maps the escaping arms).
- Base check: the new / changed tests copied onto `../tp-oi-wp08-recompute-identity-base`
  (detached `083f910`) fail 24 there (the frame unit test for the right reason; the rest mostly
  on the missing `topology()` / pip memory). Right-reason probes on base: `wp08/ids.py 16` two
  digests, `nested10.py 10` holes `[0, 12]` / 9 classes, `letters.py p3_series` 3 re-lettered.
- Review fixes: `reference_witness` rejects a frame-n witness whose own registered chain lies
  OUTSIDE the pair's images in another frame (the map preserves order along a branch; nested
  blast 10's `725 = M(439)` between `M(37)`, `M(38)` was a 4e-6 noise flip), so with the frame fix
  reverted `test_deep_reference_survives_noise` fails on the WITNESS assertion, not only on the
  class count. New unit `test_endpoint_own_unlinked_orbit_is_not_a_puncture` x2 (an endpoint's
  unlinked orbit a 1e-6 hair inside a placed frame end; both signs fail before the fix). New facade
  `test_a_pip_trimmed_away_is_not_replaced_by_classifying` (one WARNING, pip stays None,
  `bridge_classes()` raises; failed before: two warnings, the default took over). Second review: `test_a_pip_trimmed_away_is_not_replaced_by_plotting` x2 (`plot_strong_pip`,
  `plot_strong_pip_candidates`; the plot fan-outs classify through the session; the `plot_strong_pip`
  case failed before: the trellis took the default pip and `topology()` adopted it). The straddling-zone
  WARNING moved to the `tanglepack.loom.Topology` logger (its test follows). Law-case and deep-run
  signatures (`wp08/lawsig.py`) identical before and after the review fixes.
- Probe trap: `stress/numerics/harness.py` puts MAIN's `src` first on `sys.path`; import
  `tanglepack` before the harness or a probe silently runs main.

## Open-issues round WP8b (2026-10-08, branch `oi/wp08b-d10-d12`, base `014ea4b`)

- Suite: 1000 passed, 21 skipped, 3 deselected, 2 xfailed, ~161 s (994 before the review fixes) (base 974 passed / 28 skipped /
  2 xfailed). KNOWN_ISSUES 2 -> 2. NOT_APPLICABLE 28 -> 21: removed on `inversion`
  `anchor_bridge_class_leads_its_tangle`, `matrix_is_token_counts`, `partition_singletons`,
  `anchor_opening_faces_its_hole`, `backward_holes_keep_bridge_side`,
  `propagated_holes_land_on_predicted_branch`, `propagation_terminates` (all pass now);
  `every_active_class_walks` x `inversion` reworded (its active classes land on singletons).
- Series: `cases.SERIES` values are `Series(build, blasts, min_separation)` (NamedTuple);
  `blast_series` blasts the build's LARGEST zone, of its own fixed point only. New
  `k32_series` (`build_k32_zone`: the k28 recipe at k = 3.2 from the analytic seed, 10 blasts,
  ~18 s) and `nested_outer_series` (`build_nested(outer_blasts=0, through="zones")`, 11 blasts at
  `min_separation = 1e-4`, ~20 s); `inversion_series` to 4 blasts. `SeriesStep.closed` =
  `len(session.symbolic_dynamics().classes)`. `build_k28` shares `_saddle_zone` with k32;
  `_blast(..., min_separation=)`.
- Laws: `trellis_is_cut_at_the_pip_orbit` (partition layer; helper `helpers.laws.pip_orbit_cuts`
  reads the cuts off `+1` links and the registry independently; checks the per-fixed-point AND the
  selection trellis; fails on base on the three inversion cases); series
  `series_closed_alphabet_constant_while_holes_constant` (fails on base on `p3_series` and
  `k32_series`); `series_class_count_never_decreases` now fails on base on `inversion_series`,
  `k32_series`, `nested_outer_series`. Changed: `alphabet_closed_under_words` (a closure class is
  in the evidence table iff `entry.spanned`), `misplaced_unified_crossings` (no preimage IN the
  trellis).
- New: regression `test_pip_orbit_cut` (k10 every candidate: cut vs `resonance_zone` trim on a
  fresh build, same holes by cdists / class count / words in names / reliability or same raised
  type; p2inv10 the same and every class walks; a joined arc), unit `test_trellis_cut`
  (`helpers.fakes.synthetic_trellis`, now also `test_pseudoneighbor`'s `_trellis`).
- Rewritten to the new rules: `test_partition_ignores_strong_pip_cut` ->
  `test_partition_of_a_cut_trellis_ends_at_the_pip`; `test_backward_iterate_punctures_interval`
  registers `M^2(q0)` so its candidate is a crossing of T (D12); the placed-cut WARNING test reads
  the Trellis logger; `test_a_pair_no_bridge_spans_is_a_flagged_class` (was `..._is_warned_and_
  unreliable`); `test_deep_p3_alphabet_is_closed` (flagged class exists, reliable, no "no bridge
  spans" WARNING); `test_period_k_window_gaps` (`inversion` moved to the complete list);
  `test_deep_blast_sliver_beyond_pip` (the sliver crossing is not in the trellis);
  `test_trellis_is_fresh_after_every_mutation_path` (the direct build is cut at the session's
  pip); `test_a_blast_or_a_pip_change_recomputes` (the anchor-nearest alternate pip, the caches
  test's rule); `test_a_pair_with_no_bridge_object_is_skipped_with_a_warning` (classes built
  before the monkeypatch). Fakes `FakeTrellis` and `_ClassTrellis` gained `__contains__`.
  `helpers.direct.DirectBuild` cuts its selection trellis at the pips. Deleted:
  `test_unified_segment_precondition`.
- Base check: the new / changed tests copied onto `../tp-oi-wp08b-d10-d12-base` (detached
  `014ea4b`) fail there: 3 series monotone, 2 closed-alphabet, 6 deep-run, the D11 unit test, the
  placed-cut unit test, the cut law on the 3 inversion cases, `test_pip_orbit_cut` x4 (k10
  candidates 2 and 3 differ, p2inv10 differs and `n` unreachable, no joined arc).
- Review fixes: facade `test_the_selection_trellis_is_cut_before_any_product` (nested zones, one
  outer blast, the anchor-nearest alternate pip: `trellis_view(session.trellis())` before
  `topology()` equals the view after it) and series law `series_session_trellis_is_history_free`
  (`SeriesStep.history_free`: `session.trellis()` read before the recompute vs the topology's
  trellis after it; fails on the previous head on `nested_outer_series` at blasts 10, 11; the
  facade test fails there too). `test_nested_inner_circle_sits_inside_the_outer_one` passes the
  minimal trellis's own trellis to `cartoon_enclosures` and checks a sub-selection too.
- Probes: `/home/dezhu/.claude/jobs/c5c98371/tmp/wp08b_impl/` (`replay.py` rubric 10,
  `lawsig.py` id-free signatures, `deep.py` deep runs).


## Open-issues round WP9 (2026-10-08, branch `oi/wp09-heteroclinic`, base `b460363`)

- Suite: 1017 passed, 21 skipped, 3 deselected, 3 xfailed, 177.6 s (base on the same machine:
  1000 / 21 / 3 / 2 xfailed, 162.1 s). KNOWN_ISSUES 2 -> 3: removed
  `three_tangles_take_three_families` x `hand_built`; added on the new law case `heteroclinic`
  `every_active_class_walks` (WP9 slice 2, DEC-23) and `partition_covers_branch` (WP10).
- New law case `heteroclinic` (`cases.build_heteroclinic` = frozen `build_nested(16u / 10s,
  every keyword)` renamed, `CaseExpect.heteroclinic=True`); ~0.8 s per law module.
- New laws: `trellis_crossings_belong_to_its_selection` (partition: every crossing of every
  per-fixed-point trellis and of the selection has both branches in that trellis, every bridge
  ends at its crossings), `tangles_in_case_order` (classes: `session.trellis().fixed_points` IS
  `case.fixed_points`; a class of tangle t has elements and members on `case.fixed_points[t]`),
  `case_heteroclinic` (case sanity). `trellis_is_cut_at_the_pip_orbit` skips a crossing unless
  BOTH keys are branches (`<=`, was `&`).
- Stages 4-6 run on the selection trellis, so every hole / partition / pair read moved to
  `case.session.topology().trellis` (`helpers.laws._staged`; `_holes`, `_partitions`,
  `check_case_builds` "every fixed point is partitioned" = the partitions' fixed points,
  `backward_holes_keep_bridge_side`, `propagated_holes_land_on_predicted_branch` (the cycle of
  the reference BRIDGE's fixed point), `propagation_terminates`, `no_direct_hole_beyond_
  fundamental` (`k` of the origin's / pair's stable fixed point), `window_holes_on_every_branch`
  (`fp = ref.branch_key[0]`), `reference_pairs_valid` (`trellis.pips[fp]`)); `law_tier.
  fingerprint` hashes the topology's holes; `cases.blast_series` reads holes and references off
  `session.topology().trellis`; `helpers.direct.DirectBuild` picks pips per fixed point, then runs
  stages 4-6 on a fresh selection trellis; `conftest.p3_partitioned` calls `session.topology(fp)`
  per fixed point; `test_image_of_element_accepts_partitions_from_another_trellis` takes a fresh
  `Trellis.from_workbench` for its "no partitions" half. `_nested_through` constructs p1 first.
- New tests (all fail on a detached `b460363` worktree with the new tests copied in):
  regression `test_heteroclinic` (D7), `test_nested_blast_containment` (C4: one unfiltered outer
  blast's parents all classify to the outer zone and no inner-zone bridge is iterated; 4 outer +
  2 inner unfiltered blasts give the same holes per period / class count / spelled words in both
  orders, base p3 holes 4 vs 42, ~2.5 s), unit `test_rebuilding_one_fixed_points_bridges_keeps_
  the_others` (NEST-5, base KeyError). The three-family plotting test is a plain test now.
- Review round: `test_nested_blast_containment` also blasts an UNSTORED zone (k10,
  `define_resonance_zone`: parents = every un-iterated bridge whose midpoint it holds) and a
  copy of the stored outer nested zone (parents = the bridges classifying to the stored one);
  both fail on the first cut. Unit `test_collision_on_another_stable_manifold_still_punctures`
  (scaled / own frame): a hand-built A x A candidate whose numbers collide with an endpoint on
  W^s(B) still punctures; each parameter fails with its gate removed. Plotting
  `test_families_do_not_look_alike`: every pair of families (and `HETEROCLINIC_COLORS`) at CIE76
  >= 13 (the old closest pair is 13.4; the first third family had 7.1 and 8.1).
- Probes: `/home/dezhu/.claude/jobs/c5c98371/tmp/wp09_impl/` (`het_sig.py` base vs branch,
  `sig_base.txt` / `sig_branch.txt` from the planner's `wp09/sig.py`).

## Open-issues round WP10 slice 1 (2026-10-08, branch `oi/wp10-numerics-leftovers`, base `3f27456`)

- KNOWN_ISSUES 3 -> 2: `(partition_covers_branch, heteroclinic)` removed (XPASS after the law
  restatement). Left: C3 `orientation_reversing`, WP9 slice 2 `every_active_class_walks` x
  `heteroclinic` (DEC-23).
- `law_case` is session scoped (C12): law tier 35.7 -> 11.2 s. Fingerprint guard unchanged.
- `tests/regression/test_deep_reference_survives_noise.py` was vacuous (per-fixed-point trellises
  carry no stage output since WP9, empty witness loop); rewritten on
  `session.topology().trellis`. Break-it (WP8 frame fix reverted, `if False and any(image ...)` in
  `Pseudoneighbor`): the full suite fails only
  `test_unlinked_landing_is_read_against_the_pairs_images[below_c]`; no `nested_outer_series`
  test fails, so the regression is kept (plan rubric 16). It duplicates the series at blasts 9-10.
- New builders `p3_16_12`, `inversion_period2_10` (`cases._named` keeps them off the law-case
  names their builders would give).
- New tests, each failing on a detached `3f27456` worktree with the new tests copied in (base
  laws, branch builders where needed): `test_a_same_stability_pair_is_never_formed` (all four
  params), `test_crossings_sharing_a_segment_are_bracketed[2,3]`,
  `test_a_bridge_off_its_branch_is_refused`, `test_blast_through_a_shared_straddle`,
  `test_every_iterate_failure_escapes[ValueError, AttributeError]` (the AssertionError param
  passes on base by design), `test_a_rebuild_tests_only_the_remembered_pip[k28, nested]`,
  `test_outer_faces_do_not_read_areas` (4 params), `test_every_law_holds_on_the_deep_runs` (5
  params), `test_arclength_stops_on_arc_length`, `test_default_cap_scales_with_the_stretch`,
  `test_the_count_is_the_listing_of_the_new_curves` and `test_a_blast_counts_and_remembers_nothing`
  (the count function / report field do not exist on base). `test_the_budget_is_per_call` passes
  on base (it pins C6). The ninth-blast k28 count check is `perf` (~9 s).
- Reworded to new rules: `test_same_stability_discard.py` (never formed), the evidence tests
  (C11), `test_blast_restore` / `test_blast_monotonicity` / `test_blast_budget` (no `strict`).
