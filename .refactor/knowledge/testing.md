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
- Registry ids are NOT reproducible between builds (15-step failing origin seen as (12,33),
  (21,35), (3,17)); pin iterates and sides, never ids.
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
