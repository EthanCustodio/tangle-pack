# testing — knowledge file

Last updated 2026-10-16 (B12b OUTCOME, branch `refactor/laws-and-loose-ends` @ `e163f7a`: 1208 passed / 43 skipped / 5 deselected / 1 xfailed (C3) in 195 s, 1257 collected, 97 laws; see "B12b outcome" at the end. Before that B12b analyst at `a550a5d` = B12a merged: 1207 passed / 47 skipped / 5 deselected / 1 xfailed (C3) in 200 s, 1260 collected; see "B12b analyst" at the end. Before that B12 analyst at `94fca6b`: law classification, id facts, sphinx -- see "B12 analyst" at the end. Before that, B9b OUTCOME, branch `refactor/one-element-at-end` @ `94fca6b`: 1207 passed / 47 skipped / 5 deselected / 1 xfailed (C3) in 209 s, see "B9b outcome" at the end. Before that, B9b analyst at `HEAD` `5f0db5f` = B9a merged: 1208 passed / 47 skipped / 5 deselected / 1 xfailed (C3) in 199 s). Before that (B9a OUTCOME, branch `refactor/session-three-concerns` @ `5f0db5f`: 1208 passed / 47 skipped / 5 deselected / 1 xfailed (C3); see "B9a outcome" at the end). Before that (B9 analyst at `HEAD` `517aa01` (B10 merged): 1215 passed / 47 skipped / 5 deselected / 1 xfailed (C3) in 204 s; see the last section). Before that 2026-10-15 (B6b `word-finished` OUTCOME on branch `65edec4`, base `63477e3`: 1223 passed / 47 skipped / 5 deselected / 1 xfailed (C3), see the last section). Before that: B6b analyst at HEAD `63477e3` (B7b merged): 1211 passed / 47 skipped / 5 deselected / 1 xfailed (C3) in 214 s; see the last section). Before that: 2026-10-14 (B7 `refined-level-defined` outcome: 1200 passed / 46 skipped / 5 deselected / 2 xfailed on branch `a6a4f68`, base `3fb310a`, see the last section; before that B5 `inertness-is-transience` outcome (MERGED, main `3fb310a`): 1193 passed / 38 / 5 / 1 xfailed on branch `3fb310a`, see the last section; before that B3 `depth-of-n` outcome: 1193 passed / 38 / 5 / 1 xfailed on branch `7d3a294`, see the B3 bullets; before that 2026-10-12 B1 `resolved-crossings` outcome audit, branch @ `76d0710`: 1157 passed / 38 skipped / 5 deselected / 1 xfailed in ~200 s, base `1029fcc` 1117 / 38 / 5 / 1 in ~162 s; see the last section. Before that: B0 `oracle-laws`, branch `refactor/oracle-laws` HEAD `1029fcc`, base `a3bd8f0`; NOT merged at writing): suite 1117 passed, 38 skipped, 5 deselected, 1 xfailed (C3). Base `a3bd8f0` (after D15 `74a19da`): 1117 / 27 / 5 / 1 -- the 1122 / 26 of WP12 predates D15.
Run: `MPLBACKEND=Agg env/bin/python -m pytest -q -rxX -p no:cacheprovider` (~150 s; `-m perf` adds the 5 timing tests).

- Suite history: 782 at `2315204` (2026-10-02); 833 / 11 xfail after the stress integration
  (2026-10-07); 1114 / 26 skipped / 5 deselected / 1 xfailed at main `a2e0e65` (WP10b); 1122 after
  WP12 (+8 tests: restore, unmarked lobe, image end asserts, order flip, slugs, inversion recipe).
  KNOWN_ISSUES = 1 (C3, det J < 0 out of scope). No pytest-xdist in env.
- Higher-period fixtures (2026-10-02 notes; `tests/test_higher_period_cartoon.py` became `tests/regression/test_nested_own_blast.py` in the tiered-layout move `117dfb3`, B12 analyst 2026-10-15; the fixture names below are historical): (`p3_built` = `build_period3()`,
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
- Scratch-copy prototyping (B10 analyst, 2026-10-15): `git archive HEAD | tar -x -C <scratch>/proto`,
  symlink `env`, `git init` + commit, then edit and run the suite there with `PYTHONPATH=src`; the
  repo stays untouched and `git diff --stat` in the copy gives honest line counts. Running the copy's
  suite concurrently with the repo's roughly doubles both (~150 s alone).
- B10 (`dead-api`, 2026-10-15) removed 8 test ids: the five registry-graph tests of
  `test_map_step_and_graph.py` (the file keeps "graph" in its name; renaming would move ids) and
  `test_from_branch_holds_the_anchor_on_every_branch[None|unstable|stable]`; renamed
  `test_remove_iterate_forgets_the_link_and_its_twin`. Rewrites keep strength: iterate links read
  `registry.iterate_table[x, 1] == y`; `test_refinement.py` calls
  `machine._refine_layer([(p0, p1)], manifold, viewer, None, set())` at `area_cutoff = 0.0` (strict
  `lo < cdist < hi`; the function-scoped `grown_unstable` fixture absorbs the mutation);
  `test_blast_budget.py` asserts `stop_reason` only. Suite after B10: 1215 passed / 47 skipped /
  5 deselected / 1 xfailed (C3). The D4 anchor-on-every-branch rule is still covered at trellis level by the
  law at `tests/helpers/laws.py:905` (run by `tests/invariants/test_law_anchors.py`) and at registry
  level by the dedup test `at_point` in `test_intersection_registry_fixes.py`; no REGISTRY-QUERY
  test of it remains (there is no registry query by branch any more).
- Probe harness for "no number changed" (B10 test engineer): `<scratchpad>/b10/probe/{run.sh,probe.py}`
  takes `REPO` and an out dir, dumps per item (16 builds + 47 series points = 63) classes, seed, rules,
  refined rules, both transition matrices, per-class status/word/reason and `is_reliable` as sorted
  JSON; diff base vs branch. The B6b probe no longer runs on post-B6b trees (`WalkSearch.truncated`).

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
- E7 `blasts` (`e07_breaking.py` `run_blasts`): ONE config per case (k28, p3 max 10; nested 8), one cutoff,
  fixed `min_separation`; topology probed ONLY on the final state (`run_topology`) ->
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
  violation. `is_reliable` (`SymbolicDynamics.is_reliable`) = every class has an itinerary, none
  ambiguous, none contradicted by registered evidence (a `"walk"` diagnosis), no straddling
  landing; D11 flags and C11 / C21 evidence never count.
- Report framing smells (as of 2026-10-07; CLOSED, DOC-2): the "unreliable" legends and macros
  were renamed "(provisional)" in the second edition; the last smell, the "why the bursts fall
  there, not checked" / "fixed number of new classes per blast" paragraph, is answered by the
  third edition (WP12: reference switches at a pip-orbit lap; the per-blast step was BCG-1
  bookkeeping, gone since WP6).

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
- B3 analyst (2026-10-14): every series blasts with `num_iterations=1`, so `min_separation` never
  changes series state. Series stalls of the pip depth `m` within `SERIES` depths: k28 at 8, k32 at
  9, nested_outer at 9 (p3 stalls at 15, inversion at 10: beyond their 9 / 4). The only too-close
  child of `test_blast_proximity_guard`'s fixture (nested, 3 iterations, sep 0.02) is (46,47),
  ending at `f^4(q0)`: a "keep pip-bearing children" rule makes that test fail.
- B3 outcome (branch `refactor/depth-of-n` @ `7d3a294`, base `9cd29d5`, not merged): suite 1193
  passed / 38 skipped / 5 deselected / 1 xfailed (~245 s; base 1183). New: `SeriesStep.pip_depth`
  (`cases._pip_depth`, composed through `trellis.iterate`), law `test_pip_orbit_never_shallower`
  (`depths == sorted(depths)`, five series, passes on base too: the stalls are flat, never
  shallower), `tests/unit/loom/test_blast_pip_orbit.py` (kept with pip / dropped without /
  `pip_depth` = composed orbit / exemption skips ONLY the separation test -- the last one the only
  guard of the exemption's position after the exterior check), and the proximity-guard test
  parametrized `pip` / `no_pip` (`too_close > 0` asserted in `no_pip` only). Smell: the test-side
  orbit composition is written three times (`cases._pip_depth`, `_pip_orbit` in both unit files).
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

- KNOWN_ISSUES 3 -> 2: `(partition_covers_branch, heteroclinic)` removed (XPASS after the law's
  last clause was rewritten). Review round: the first rewrite read `hole.openings` with the
  partition builder's own filter and restated `StablePartition` (a break-it flipping every
  non-anchor opening in `_hole_openings` passed it on k10, p3, k28_two_blasts); the clause now
  derives each hole's opened halves from its bounds, `bridge_side` and crossing signs
  (`_sign_openings`, `_anchor_row` at an anchor) and fails that break-it on all three. The
  derived halves equal `hole.openings` on every law case (anchor openings exercised on 9 of 11).
  Left: C3 `orientation_reversing`, WP9 slice 2 `every_active_class_walks` x
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
- Suite (back to back, same machine, load ~1.7): branch 1042 passed / 21 skipped / 4 deselected /
  2 xfailed in 137.3 s; base `3f27456` 1022 / 21 / 3 / 3 in 184.3 s (0.75x). Savings: law tier
  setups (C12, ~25 s), one rtree per stability (~7 s on the series), lazy pip candidates.
  Slowest left: the `k32_series` and `nested_outer_series` builds (~15 s each, shared by the
  series tests), the deep-reference regression (~9 s, a duplicate of the series at 9-10 blasts,
  kept per rubric 16) and `test_deep_blast_sliver_beyond_pip` (~10 s).

## Open-issues round WP9b (2026-10-09, branch `oi/wp09b-noncrossing`, base `40175fb`)

- KNOWN_ISSUES 2 -> 2: `(every_active_class_walks, heteroclinic)` stays strict (only `bn` is ambiguous now:
  member-less, connecting, NUM-21), its reason rewritten and its tag `WP9` -> `DEC-24`. The registry test
  accepts an open decision as a tag (`WP\d+|C\d|DEC-\d+`; `KnownIssue.wp` docstring says so).
- New laws (all 11 law cases pass, every one non-vacuous): symbolic `walks_agree_with_registered_images`
  (every evidence verdict agrees, no class `verified` False), `member_lobe_classes_have_one_walk` (a walked
  class with a member lobe -- kept non-loop member, both ends on one branch by `_stable_branch_at` -- has a
  `unique` or `trivial` walk: strengthened during the run from "not ambiguous" after the global-side-flip
  break-it left classes `no_admissible_walk` unnoticed), `walk_starts_on_the_landing` (the first / last face is
  beside a node the span meets, the span error derived in the test); dual graph `lobe_enclosure_is_geometry`
  (flood vs `point_in_polygon` on the lobe polylines, periodic points within 1e-6 skipped),
  `registered_images_keep_markers` (a kept bridge's lobe vs its fully kept image chain closed along S).
  `_closed_curve` closes at an anchor along the other end's branch, so inversion_period2 counts (36 checks;
  0 with the C7 branch).
- Series: `SeriesStep.ambiguous` and `test_no_walk_is_ambiguous` (fails on base on `k28_series` [blast 8: 2]
  and `k32_series` [8: 2, 9-10: 3]). Regression `test_deep_blast_sliver_beyond_pip` also checks the landing,
  ties and evidence (fails on base at `walk_starts_on_the_landing`). Unit (`test_dual_walk.py`, FakeDual,
  `holds_markers` monkeypatched on the module): span picks the start node (3 params: inside, touching the
  neighbour, scaled end past it), span meeting nothing falls back with one WARNING, the markers break a tie,
  a None verdict next to an admissible walk stays a tie, a detour through a face twice, no admissible walk.
- Fails on base (detached `40175fb` with the new tests copied in): `test_symbolic_laws[heteroclinic]`
  (walks_agree: `c`; ambiguous_only: `a`), `test_no_walk_is_ambiguous[k28_series, k32_series]`,
  `test_deep_blast_sliver_beyond_pip`; every `test_dual_graph_laws[*]` with AttributeError (`enclosed`).
- Break-its (`.../tmp/wp09b_impl/breakit.py`, `breakit2.py`): holes at the trellis-lobe corner fail
  `registered_images_keep_markers` on 9 of 11 (not p3, inversion_period2); dropping the span restriction fails
  the k28 series test (blast 8: 2 ambiguous) and the deep-sliver regression; flipping every hole side fails
  `member_lobe_classes_have_one_walk` on k10, p3, nested, heteroclinic; orbit points mapped by +2 pass
  every law (not observable on these fixtures; the unit test below guards it).
- Review fix: the law first named `ambiguous_only_without_a_member_lobe` is `member_lobe_classes_have_one_walk`
  (it asks for `unique` / `trivial`, so it also catches `no_admissible_walk` and `truncated`). New unit tests on
  rules the laws could not see (the reviewer's mutants all survived the law tiers): `test_enclosure.py`
  (real `Arc` / `Region` / `FaceNode` objects around one bridge lobe on a `DualGraph` made with
  `object.__new__`: the lobe holds the hole beside its bridge and not the anchor, nested chords give one
  routing, interleaved chords `[]`, an orbit marker maps to the orbit index `FixedPoint.advance_key` names on
  periods 1-3 with and without inversion, a hole marker to its next iterate only when placed), and in
  `test_dual_walk.py` the `sigma'` rule (a `FakeDual` whose enclosure holds nothing, so only `sigma'` can
  reject: the walk across `[1, 2]` is False and loses to the one round the goal), a budget spent with nothing
  admitted is `truncated`, and an admitted itinerary's other walks are counted, not judged. `FakeDual` nodes
  now carry real `Arc`s (no geometry) and the dual has `stable_nodes_on`. Mutants killed: `sigma'` -> True,
  `enclosed` crossing -> `[{}]`, `_cut_ring` crossing check removed, orbit `+0` / `+2`. (`return []` ->
  `continue` in `enclosed` is equivalent: `_cut_ring` does not depend on the routing.)

## Open-issues round WP11 slice 1 (2026-10-09, branch `oi/wp11-simplify`, base `4f7e815`)

- KNOWN_ISSUES 2 -> 2 (C3 det J < 0; DEC-24). NOT_APPLICABLE unchanged.
- Proof that nothing moved: `/home/dezhu/.claude/jobs/c5c98371/tmp/wp11/check_all.sh <root> <tag>` runs
  `sig.py` (14 cases + 5 deep builders, exact and id-free `describe()` hashes, holes, partitions,
  classes, reliability, WARNING counts per logger; `PYTHONHASHSEED=0`), `cartoon_sig.py` (every
  PathPatch / scatter / line / text of both cartoon shapes and colourings on four fixtures),
  `series_sig.py` (every `SeriesStep` of the five series) and `tol_count.py`; all identical to base
  after every step. ~55 s.
- The tolerance experiment behind C15: `sig_tol.py`, `series_sig_tol.py` and the pytest plugin
  `plug/tolzero.py` force the `cdist_tol` read from topology / loom frames to 0 (numerics keeps 1e-6):
  every fingerprint and the full suite were unchanged on base.
- `coverage_guard.py` keys lines by (qualname, stripped text, occurrence) and reads the base source
  with `git show <base-ref>:path`; its own tests are `tests/_tools/test_coverage_guard.py` (collected).
  The dead-API whitelist and the conftest "Dead API" paragraph are gone.
- Tests rewritten to new rules (say so in any later review): `test_stable_partition.py` orientation
  tests set `fp.partial_jacobians[0]` (and one checks a contradicting dynamical-system jacobian is
  ignored); `test_topology_invariants.py::test_orientation_from_the_single_step_jacobian` (the DEBUG
  default fires for every unsolved fixed point); `test_exact_ownership.py`: `interval_at` exact for an
  estimated cdist, a boundary at TOL/2 gets a synthetic anchor below it, no INFO;
  `test_dual_walk.py::test_image_on_an_unexpected_branch_asserts`.
- New-test base check: copy `src`, `tests`, `pyproject.toml` of the detached base worktree
  (`../tp-oi-wp11-simplify-base`) to `.../tmp/wp11/basecopy`, overwrite the new / extended test files,
  run there (never via PYTHONPATH).

## Open-issues round WP9c (2026-10-09, branch `oi/wp09c-every-bridge-the-same`, base `5c7070e`)

- KNOWN_ISSUES 2 -> 1: `(every_active_class_walks, heteroclinic)` XPASSed (the walls decide `bn`) and was
  removed; left `(orientation_reversing_case_builds, orientation_reversing)` (C3). NOT_APPLICABLE +5.
- New symbolic law `walks_avoid_disjoint_images` (D13, D14): two classes with decided walks whose arcs lie on
  disjoint unstable-branch sets (decided = `unique` / `trivial` of `multiplicity == 1`, the library's
  rule; re-derived in the test, `laws._arc_branches`: members' branches, or one
  `advance_key` step past the namer's) have no chords interleaving on a face ring (`walk_chords`, checked
  with `DualGraph.chords_cross` AND an independent ring-position test, `_chords_interleave`). It counts the
  CLASS pairs judged, not the chord pairs: `test_every_law_holds_on_the_deep_runs` runs every law on
  `inversion_period2_10`, which has 54 such class pairs and none sharing a face, and NOT_APPLICABLE may only
  name law cases (`test_law_case_sanity`). NOT_APPLICABLE on k10, k28_one_blast, k28_two_blasts, inversion,
  inversion_period2 (no such class pair). Non-vacuous on p3 (1 pair), nested, nested_inner_blasts,
  heteroclinic, inversion_one_blast / _two_blasts (pairs with no shared face).
- D14 clause in `alphabet_closed_under_words` (C18): every pair of fixed points a class connects is connected
  by a computed bridge of the selection trellis (`_staged(case).bridges`, ends read by `_stable_branch_at`);
  exercised on heteroclinic only (3 connecting classes, one fixed-point pair).
- Regressions (`test_heteroclinic.py`): `test_a_tie_no_lobe_decides_is_broken_by_the_images_of_disjoint_arcs`
  (every class across two fixed points walks one way, some class without a lobe has `rejected > 0`, no class
  ambiguous; `is_reliable` is not pinned, conftest policy) and, under `perf`,
  `test_registered_connecting_images_agree_with_the_wall_rule` on the new builder `heteroclinic_17u`
  (`build_heteroclinic(p3_unstable_steps=17)`, 27 s; reproduces `wp09c_17u_base.txt`: `bm` ambiguous on base).
- Unit: `test_dual_walk.py` (a wall breaks the `tie`, a wall crossing every walk keeps it, a walk of
  multiplicity 2 is never rejected (fails with the multiplicity guard removed), no walls / walls in
  an unvisited face change nothing, `walk_chords` end to end; `chords_cross` stubbed on the `FakeDual`
  instance), `test_enclosure.py` (`chords_cross` interleaved / nested / shared end edge / apart, both orders,
  and the same verdict as `enclosed(...) == []`), `test_symbolic_dynamics.py` (`_arc_branches`: a member-less
  class one step past its namer on period 3 and inversion periods 1 / 2, namer = one branch so every
  parameter fails without the step; a member class on its members').
  Scaffolding updated: the synthetic `_Crossings` trellis answers `unstable_key_of`, and the `_resolve_itinerary`
  stub takes `walls`.
- Base check: the new and extended test files plus `cases.py` copied over a `git archive` of base
  (`.../tmp/wp09c_base_tree`): all 14 new tests fail there (12 on the missing API, the two regressions on
  `bn` / `bm` ambiguous); the existing tests of those files pass.
- Suite 1102 passed / 26 skipped / 5 deselected / 1 xfailed in 132.0 s after the review fix round (WP11 base: 1088 / 21 / 4 / 2 in 132 s).

## Open-issues round WP10b (2026-10-09, branch `oi/wp10b-leftovers`, base `f9d1be0`)

- KNOWN_ISSUES 1 -> 1, NOT_APPLICABLE unchanged.
- Mutating law tier (`MUTATING_LAYERS["recompute"]`): NEW first law `recompute_keeps_the_topology`
  (D1, NUM-12: snapshot of sorted `(bridge id, iterated)`, the number of partial pieces, sorted
  `(hole.bounding_ids, iterate, bridge_side)`, `bridge_classes().describe()` and
  `symbolic_dynamics().describe()`, around `compute_intersections(preserve_ids=True)` +
  `rebuild_bridges()`; counts bridges); `recompute_preserves_ids` STRENGTHENED to every crossing
  back under its own id (`set(after) == set(before)` and `_same_crossing` per id; the "may drop out"
  caveat is gone). On base both fail on k28_one_blast, k28_two_blasts, nested, nested_inner_blasts,
  inversion_one_blast, inversion_two_blasts, heteroclinic (the layer reports the first; the ids law
  alone was checked with `.../tmp/wp10b/ids_law_base.py`).
- Laws `area_along_iterate_links` / `iterate_link_scales_by_beta` (NUM-20) read each side's own
  saddle: product factor `beta_u(A) * beta_s(B)` (exactly 1 homoclinic), stable prediction
  `beta_s(B) * s`. No law case has a heteroclinic link (NUM-21).
- Regressions: `test_growth_after_a_blast.py` (k28 grown in both stabilities, nested in its stable
  manifolds only: p1's unstable arm escapes and costs 2 s; every pre-growth crossing re-found within
  its own unstable segment length, class / closed counts non-decreasing; one blast + recompute +
  blast == two blasts, id-free `describe()`), `test_heteroclinic_iterate_links.py` (unblasted
  `build_nested(outer_blasts=0, 13u, 16s, through="intersected")`: 12 heteroclinic links, base laws
  fail on 4 -> 8), `test_foreign_chords.py` (`p3_6_blasts`, `p3_15`, `p3_16`, `heteroclinic`: no
  "no empty stretch" cut; Jordan check by point-in-polygon on the lobe image, each IN stretch an
  empty stretch; non-vacuous on all four).
- Unit: `test_iterated_cut.py` rewritten to C19 (`_FakeTrellis` crossings carry `stable_cdist`,
  `is_anchor`, `manifold_a_key`): stable-order foreign chords with the p3_6_blasts cdists, the base
  read at a non-anchor end (inversion_period2 pattern), a base across two branches asserts, the
  foreign-far-end injection on nested now asserts. `test_grow_until.py`:
  `test_a_growth_round_keeps_the_blast` replaces `test_grow_until_clears_the_iterated_flag_it_invalidates`
  (the old test pinned the WARNING that disappears). `test_bridge_identity.py`:
  `test_a_reversed_bridge_id_is_refused`. `test_stable_trim.py`: exact-node and past-the-end cases
  (pass on base). `test_plot_topology.py`: ranks strictly increasing (no `layout.tol`). Deleted:
  `test_closed_form_curvature.py` (with `_parabolic_fit`).
- Base check: every new / edited test file copied over a `git archive` of `f9d1be0`
  (`.../tmp/wp10b/basecheck`): 24 failures (7 recompute-layer cases, growth x2 + order, foreign
  chords x8, grow_until, iterated cut x4, reversed id) and the heteroclinic-link regression with
  the base `laws.py`; the stable-trim and plotting tests pass there.
- Suite 1113 passed / 26 skipped / 5 deselected / 1 xfailed in 143.8 s (base 1102 / 26 / 5 / 1 in
  137.9 s back to back; a first branch run at 152.7 s led to slimming the new regressions);
  `-m perf` 5 passed in 44.9 s.


## Open-issues round WP12 (2026-10-09, branch `oi/wp12-reverify`, base `a2e0e65`)

- The stress harnesses were BROKEN on main: `e11_blast_depth.build` called the deleted
  `hc._grow_period3` / `hc._grow_period1`, `topology_snapshot` / `e07.run_topology` read the
  deleted `dynamics.virtual_classes`, `e07.blast_setup` used the deleted growers and
  `run_depth` called `session.infer_iterates()`. Fixed: every case is a frozen `tests/cases.py`
  builder (`e11.build` also builds `inversion` and `heteroclinic`), `repin` deleted (the session
  carries the pips), snapshots read ONE `session.topology()` and record `n_closed`,
  `n_trellis_crossings` (N), `n_flagged` (D11), `n_disagree`, `n_straddle`; PERF-5: the cheap
  `unresolved` counts the blasted fixed point's crossings only, `e07 run_tangency` runs the
  resolution pass once (`compute_intersections(diagnose=False)` then `check_resolution`). The
  `--quick` grids of e07 / e11 / e12 end all `ok` (run to scratch `--out` files so the second
  edition's `results/e07_*.jsonl`, `e11_*.jsonl` stay untouched).
- E12 (`stress/numerics/e12_reverify.py`, results `results/e12_reverify.jsonl`, 31 cells, all
  ok, ~4 min wall at 8 workers): series to depth at 1e-7 and 1e-5, NUM-6 / NUM-7 two-iteration
  cells, growth depth. Report: `analysis.e12_*`, `make_figures.fig_e12_classes` /
  `fig_e12_cost`, `make_numbers.e12` (`tables/e12_summary.tex`, `e12_iterations.tex`,
  `e12_depth.tex`, `\NumWpTwelve*` macros), the third-edition section of `report.tex`; built with
  `~/.local/bin/tectonic report.tex` (no latexmk / pdflatex on this machine).
- New tests (each fails on a `git archive` copy of base with the new tests copied in):
  `test_restore_after_stable_growth_never_shortens`, `test_an_unmarked_lobe_is_skipped_by_design`,
  `test_an_image_end_off_its_branch_asserts`, `test_an_order_flip_gives_no_verdict` (~8 s,
  `build_k28(blasts=9)`), `test_every_period3_build_gets_its_own_folder`,
  `test_the_example_inversion_recipe_is_the_law_cases[0, 1]`,
  `test_every_inversion_build_gets_its_own_folder`; extended (passes on base):
  `test_disagreeing_evidence_carries_its_cause` (independent order-flip scan) and the allowed
  reasons of `test_cuts_record_the_far_end_of_the_empty_stretch`. The synthetic `_drive` harness
  of `test_symbolic_dynamics.py` gives its homotopy elements no bounds (no order flip to read).
- No law changed, no series extended: an `inversion_series` to 6 blasts would need a new
  KNOWN_ISSUES entry naming DEC-25 (BCG-8). KNOWN_ISSUES stays 1.

## B0 `oracle-laws` (2026-10-12, branch `refactor/oracle-laws`, base `a3bd8f0`, HEAD `1029fcc`)

Tests only (`git diff a3bd8f0 -- src` empty). Four ORACLE laws that re-derive a definition without
calling the code under test (fidelity audit section 5: a third of the laws could only fail if the
code reported its own failure). Line numbers are the branch's.

- **Map** (`tests/helpers/laws.py`):
  - `check_chosen_pip_is_strong` (`:1566`, layer `partition`) -- A.4's open-box test from the
    registry + `FixedPoint` alone: every homoclinic non-anchor crossing of the pip's fixed point on
    its stable cycle, carried back `(pos(r) - pos(q0)) mod k` steps by `per_step_beta` (stable cdist
    grows, unstable shrinks), must lie outside `(0, c_s(q0)) x (0, c_u(q0))`; q0's own orbit skipped
    by a 2-D collision within `SCALING_RTOL`. No iterate table. Guards `StrongPip.is_strong_pip`.
  - `check_reference_pairs_clear_of_X` (`:1632`, `partition`) -- both open intervals of every
    reference pair (Q1: the stable one too), each with its OWN key's cycle and beta (D7), scaled
    landings of every non-anchor crossing of the cut N (`own_intersection_ids`) via `_laps`; a landing
    is allowed only if `_table_image` links it to an end or `_scaled` collides with an end on both
    cdists. Guards `Pseudoneighbor._interval_is_clear` (both sides since B2).
  - `check_propagated_hole_inside_container` (`:1672`, `partition`) -- every hole with
    `iterate < 0`: the reference bridge's ends (`hole.origin`) carried back by the REAL
    `map_inv_batch`, projected on the container's (`bridge_between(*bounding_ids)`) node polyline,
    strictly inside its crossing cdists; CLOSED at an end the table links to a container crossing
    (`:1690`; the real map overshoots an exact crossing by 1e-7..7e-4 relative: k28_one_blast,
    inversion, inversion_period2 at 9/9). Guards the backward chain / `_containing_bridge`.
    Bounds NO distance to the polyline (open point, see Smells).
  - `check_k10_is_the_full_shift` (`:3024`, `symbolic`) -- Q11, the ONE textbook pin: returns 0
    unless `case.name == "k10"`; there one active class, two self-tokens of opposite direction, two
    occurrences, class matrix `[[2]]`, refined `[[1,1],[1,1]]`. The only law guarding the CLASS
    matrix (`check_matrix_is_token_counts` reads the refined one, all 1s on law cases).
  - helpers `_key` `:1591`, `_cdist` `:1596`, `_scaled` `:1601`, `_laps` `:1612`, `_table_image`
    `:1619` (direct entry else composed +-1 links).
  - `LAW_CASES_ONLY` (`:3172`) = `("reference_pairs_clear_of_X", "k10_is_the_full_shift")`:
    `tests/regression/test_open_p3_deep_runs.py:94-98` loops `assert_non_vacuous` over every
    `LAYERS` check except these (scaling alone is wrong at depth, A.5).
- **Cases** (`tests/cases.py`): `build_inversion_period2` default `stable_steps=10` (`:825`, Q12):
  6 classes, 10 holes, `is_reliable` True (9/9: same counts, False). NOT_APPLICABLE +11:
  `k10_is_the_full_shift` on the 10 other law cases (`:1339`), `every_active_class_walks` on
  `inversion_period2` (`:1318`, every active class lands on a singleton). KNOWN_ISSUES unchanged (1).
- **Policy:** `tests/conftest.py:77-79` records the Q11 exception to the no-pin rule.
- **Teeth (scratch, not committed):** check 1 red only via a driver forcing a rejected pip
  (k10 pips 5, 6; the default `min(unstable_cdist)` still picks a strong one); check 2 red on 7/11
  with `_unstable_interval_is_clear` (B2: `_interval_is_clear`) forced True (k28_one_blast, p3, inversion, inversion_period2
  have no punctured pair, so unobservable there); check 3 red on 9/11 (vacuous, hence red, on the
  other 2) with the container shifted; check 4 red on k10 for an unsplit refine and for
  `transition_matrix` weight -> 1.
- **Smells:** `laws.py:3021-3023` three blank lines, `:3042` one (PEP 8: two); `_scaled` used once
  (inlinable), the nested `any(all(...))` at check 2 wants a named local; the deep-run loop
  (`test_open_p3_deep_runs.py:94-98`) bypasses `run_layer`, so it stops at the first failing law
  and ignores NOT_APPLICABLE / KNOWN_ISSUES (none exist for deep builds today); check 3's missing
  distance bound -- measured max nearest distance of a non-closed carried end / polyline length:
  k10 1.1e-3, heteroclinic 6.1e-3, so a bound needs a constant ~1e-2 of length; the tolerance-free
  alternative "container = nearest bridge to the carried end" is the author's call.
- A.0/A.6's "the table never holds these points" is FALSE on k28_one_blast, inversion and
  inversion_period2 9/9 (`iterate(5,-1) = 3`, `iterate(3,-2) = 1`); noted for B4b.

## B1 `resolved-crossings` (2026-10-12, branch `refactor/resolved-crossings` @ `76d0710`, base `1029fcc`, NOT merged at writing)

- **Deleted:** `tests/unit/numerics/test_noise_crossing_collapse.py`. **New:**
  `tests/unit/numerics/test_resolved_crossings.py` (hand-placed polyline workbench: zig-zag of three
  detections registered / out of N / one bridge across / no link / one WARNING; symmetric pair
  judged by the larger deviation; clean crossing resolved; `threshold=None` keeps clause 2 and
  admits a detection below 3; law cases: N = registry).
- **`cases.py`:** `SeriesStep.crossings` / `recomputed` = `(N, bridge ids)` before / after a
  `compute_intersections(preserve_ids=True)` + `rebuild_bridges()` at EVERY step (the series now
  mutate state between blasts; exact on all 5 series). `build_inversion_period2(area_cutoff=1e-5)`
  keyword (law case unchanged); `BUILDERS["inversion_period2_10"]` at 1e-6;
  `build_period3_coarse_blasts` sets `resolution_ratio_threshold = None` (clause 1 off).
- **Regression:** `test_class_count_under_blasting.py` + `test_recompute_keeps_the_crossing_set`,
  `test_every_bridge_ends_at_crossings_of_n` over `SERIES`. `test_evidence_diagnosis.py`'s order
  flip hand-blasts 9 times from `build_k28(blasts=0, through="zones")` at threshold None.
  `test_pip_orbit_cut.py`: `_p2inv10` at 1e-6 with a "resolves every detection" precondition;
  new `test_p2inv10_coarse_cut_holds_what_is_resolved` keeps the 1e-5 build (38, 39 out of N inside
  both cuts) -- the ONLY test exercising the oracle's N skip.
- **Tool:** `tests/_tools/resolution_probe.py` (`--laws`, `--builders`, `--series`, `--k28x8`;
  monkeypatches `Tangle._resolution_ratio` to recompute `h`, runs on base and branch alike;
  `PYTHONPATH=src:tests`).
- **Weak spots (test engineer, read):** the two threshold-None builds guard NUM-1 / C21 only with
  clause 1 off (unreachable at the default); `test_p2inv10_every_class_walks` alone misses a
  stable-only cut regression (four sibling tests catch it).

## B6a `walk-bounds-honest` (2026-10-12, audited on branch @ `8a1cdc3`)

- `test_a_span_meeting_no_edge_falls_back_to_every_edge` is REPLACED by
  `test_a_span_meeting_no_edge_leaves_the_walk_undecided` (`span_meets_no_node`, one WARNING); new
  `test_a_span_narrower_than_its_error_starts_from_the_edge_it_lies_in` (second tier) in `test_dual_walk.py`.
- `test_symbolic_dynamics.py`: `test_members_holding_different_markers_leave_the_class_unwalked` (stub dual,
  spy on `shortest_walks`, 0 calls), `test_a_truncated_search_decides_nothing_and_is_no_wall` (FakeDual tie
  at `max_walks=1`, real driver), `test_a_budget_too_small_to_decide_a_walk_leaves_its_class_unresolved[1|2]`
  (real heteroclinic, ~1 s, asserts its own precondition). `_drive` gained `dual_graph=` / `max_walks=`;
  with `walks=None` it walks for real and MUTATES the dual it is given (sets trellis / partition / minimal).
- Real-data guard of the second tier: `test_closed_alphabet_constant_while_holes_constant[k32_series]`
  fails 42 == 45 without it. Suite on the branch: 1162 passed, 38 skipped, 5 deselected, 1 xfailed.
- Base-vs-branch probe scripts: `.../scratchpad/te/pin.py` (TE), `.../scratchpad/audit/r6.py`, `k32.py` (auditor).

## B8 `unresolved-sets` analyst (2026-10-12, HEAD 8a1cdc3)

- Iterate-table laws live in `tests/helpers/laws.py`: `check_iterate_table_partial_bijection` :742,
  `check_iterate_links_pass_the_map_gate` :757, `check_iterate_table_complete` :770 (an unlinked
  detected crossing with a bracketed, gated image fails: it would fail on any multi-match or
  collision source, so it must exempt `unresolved_links` once those exist); registered under
  layer `"crossings"` at :3066-3076 (`tests/invariants/test_law_crossings.py`), and called directly
  by `tests/regression/test_iterates_after_blasts.py:26-29` (p3 ten blasts).
- The inference unit module is `tests/unit/numerics/test_iterate_inference.py` (2 tests, fixture
  `henon_tangle_with_bridges`, `assert_logged` from `helpers.logs`); also
  `test_intersection_registry_fixes.py:125` (`test_iterate_table_is_a_partial_bijection`).
- Facade pip tests: `tests/facade/test_session_topology.py` (`test_a_blast_or_a_pip_change_recomputes`
  :38, `test_pip_choice_survives_a_blast` :61, `test_a_rebuild_tests_only_the_remembered_pip` :122).
- `near_tangency` is the only corpus case where the multi-match fires (25 N sources); only
  `tests/regression/test_resolution_warning_follows_oracle.py:43` builds it.

## B8 `unresolved-sets` outcome (2026-10-13, branch `refactor/unresolved-sets` @ `ae18859`, base `8a1cdc3`)

- Suite on the branch (auditor rerun): 1171 passed, 38 skipped, 5 deselected, 1 xfailed (C3 only); base 1162.
- `check_iterate_table_complete` (`tests/helpers/laws.py` ~:770) now first checks every
  `registry.unresolved_links` entry (>= 1 candidate, no `+1` link, a `(j,)` entry with no `(j, -1)`
  link) and skips unresolved sources in the completeness scan. The set is empty on every law case and
  series step, so these clauses are guarded ONLY by the review probe
  (`scratchpad/te/law.py`, `lawbreaks.py`), not by a committed test (plan's choice). Still scans the
  whole registry, not N (fails on `near_tangency` crossing 1, base too; not a law case).
- New unit tests (`tests/unit/numerics/test_iterate_inference.py`, now 8): collision across passes
  (hand-linked `z -> y`, then `infer_iterates([x])`), three claimants in one pass (monkeypatched
  `_bracketed_image` / `_map_gated`), multi-match recorded + silent second pass (monkeypatched
  `_brackets`), `remove_iterate` forgets twin and edge, a later pass does not relink a contested image,
  a pass re-decides its own sources. One known equivalent mutant: an in-order loop that keeps the live
  `earlier` term gives the same outcome as grouping.
- New facade tests (`tests/facade/test_session_topology.py`): `test_a_zone_at_a_non_candidate_keeps_the_chosen_pip`
  (inversion x2, crossing chosen by rule: outermost non-candidate off the pip's branch),
  `test_a_forgotten_pip_gives_the_default_back` (k10), `test_a_stale_pip_error_names_the_remedy`.
- Fingerprint comparison (rubric 4) must drop `products` (tuple of `id()`s) and read
  `workbench generation` as an observation counter: it moves +1 per extra trellis read
  (B8 zone path: 2 -> 3 / 2 -> 4) with every state aspect equal.

## B4a `holes-algorithm-vs-diagnostics` analyst (2026-10-13, HEAD `ae18859`)

Tests that read the names B4a moves, the two `Hole` fields, or the StablePartition LOGGER:
- imports by module path: `tests/helpers/laws.py:70` (`_side_of`, `bridge_row_violation`),
  `tests/unit/topology/test_stable_partition.py:24` (`_bridge_side_of`, `_side_of`),
  `tests/plotting/test_plot_topology.py:25` (`_side_of`), `tests/unit/topology/test_topology_invariants.py:44`
  (seven names, all public, could come from `tanglepack.topology`).
- `Hole(measured_side=, side_source=)`: `test_topology_invariants.py:61-73` (`_hole`, used by
  `test_measured_side_disagreement_is_a_flag_not_an_error` :155, `test_side_diagnoses_name_each_disagreeing_hole`,
  `test_unmeasured_propagated_holes_are_logged` :197).
- field reads: `test_production_checks.py:118-133` (`test_propagated_holes_carry_their_orbit_side`),
  `:136-154` (`test_punch_holes_flags_a_contradicting_measured_side`, mutates `hole.measured_side`),
  `tests/helpers/laws.py:1209` (`check_holes_are_classified` asserts `side_source`).
- logger name `"tanglepack.topology.StablePartition"` for the moved WARNING/INFO:
  `test_production_checks.py:29` (WARNING count 1), `test_topology_invariants.py:200-202` (INFO from
  `side_diagnoses`), `tests/regression/test_coarse_blasts_rows_by_signs.py:25` (counts geometry
  WARNINGs, asserts `> 2 * quiet`) -- all break when the warnings move to `topology.diagnostics`'s
  logger. `test_exact_ownership.py:150` captures an INFO that stays.
- spies `stable_partition.check_*` (`test_production_checks.py:61-92`): survive iff the checks stay
  in StablePartition and `Trellis.punch_holes` keeps importing them from there at call time.
- outside tests: `scripts/refactor_proof/hole-side-and-own-blast/hole_sides.py` uses
  `sp._bridge_side_of`, `sp._oriented_bridge_polyline` (archived proof script);
  `stress/numerics/report/analysis.py:681-689` classifies failures by the I1/I2 MESSAGE texts.

## B4a outcome audit (2026-10-13, branch @ `f166573`)

- Logger for geometry WARNINGs and the skip INFO is now `tanglepack.topology.diagnostics`
  (`test_coarse_blasts_rows_by_signs._LOGGER`, `test_production_checks._LOGGER`,
  `test_unmeasured_propagated_holes_are_logged`).
- `test_topology_invariants._hole(origin, iterate, side, *, direct=False)`: direct = has a
  `PseudoneighborPair`; measured sides are passed as `side_diagnoses(..., measured={(origin,
  iterate): side})` or by monkeypatching `diagnostics.measured_sides`.
- New (test engineer): `test_side_source_names_the_step_that_made_the_hole` (spies which step made
  each hole, not a restatement of `side_source`) and `test_measured_sides_stay_off_the_assertion_path`
  (`diagnose=False` -> 0 calls of `measured_sides`, default -> 1).
- `laws.check_holes_are_classified` no longer asserts `side_source` (it restated the definition).
- The `window_holes_on_every_branch` note above says "`side_source == "signs"`": read it as
  `diagnostics.side_source(hole) == "signs"`, i.e. `hole.pair is not None`.
- Suite at `f166573`: 1173 passed / 38 skipped / 5 deselected / 1 xfailed.

## B4b `holes-are-combinatorial` analyst (2026-10-13, HEAD `f166573`)

- Removing `Hole.coords` (a required, first positional field) breaks the following:
  - `Hole(coords=...)` in `test_topology_invariants.py:68` (`_hole`), `test_stable_partition.py:82`
    and `:392`, and `test_minimal_trellis.py:71` (`ghost`);
  - in `tests/helpers/laws.py`, `check_holes_are_classified` (:1213, finite coords),
    `check_backward_holes_keep_bridge_side` (:1254, carries `origin.coords` back by `map_inv`) and
    `check_direct_hole_side_is_coordinate_side` (:1280); the last two are geometric ORACLES on the
    reference placement;
  - `test_stable_partition.py:355-373`, which reads propagated `hole.coords`.
- `test_production_checks.py` drives `measured_sides` on a real trellis, so it is unaffected once
  `measured_sides` computes the reference point itself.
- Nothing tests the propagation fallback, `seen`, `punched_regions` or `_region_key`; deleting
  them removes only uncovered lines. The `_near_far` tests (`test_stable_partition.py:529-549`)
  stay with the helper, which stays.
- No series or regression test inspects topology WARNINGs. The proposed scaled-end WARNING would
  fire on `k32_series` and `k28_deep_blast_sliver`, with nothing in the suite to see it.
- Baseline at HEAD: 140 passed in 9 s (test_stable_partition, test_topology_invariants,
  test_minimal_trellis, test_production_checks, test_plot_topology, test_plot_smoke).

## B4b outcome audit (2026-10-14, branch @ `ba23a1f`, base `f166573`)

- Branch suite: 1178 passed / 38 skipped / 5 deselected / 1 xfailed (base 1173; +3 coder, +2 TE).
- New: `test_a_backward_pair_is_never_punched` (hand-built `iterate=-1` pair -> AssertionError in
  `punch_holes`); `test_scaled_span_end_beside_a_crossing_warns[scaled|table]` (stubs
  `trellis.image_cdist` on the instance, synthetic 4-crossing branch; scaled end at
  `2.0*(1-0.5*SCALING_RTOL)` warns once, table end does not); TE's
  `test_a_hole_is_drawn_off_its_bridge[map|no_map]` (every `hole_point` strictly off its own bridge
  polyline; guards the no-map fallback, which nothing else did).
- Renamed: `test_henon_propagated_hole_side_matches_coords` -> `..._point_is_on_its_side` (reads
  `hole_point`). `Hole(coords=)` gone from `_hole` helpers in three test files.
- The three geometric laws now read `plotting.hole_point`: the law oracle and the figure share code,
  so a placement bug would move both together (the base-vs-branch 0.0 comparison is the guard,
  probe-only, not committed).
- Not committed as tests: law-case silence of the NUM-5 WARNING (law tier builds once per session,
  no build-time log capture) and the 1e-12 placement equality (needs the base worktree).

## B2 analyst notes (2026-10-14, HEAD ba23a1f)

- `KNOWN_ISSUES` (`tests/cases.py:1321-1330`) holds ONE entry (`orientation_reversing_case_builds`);
  none for `reference_pairs_clear_of_X` (B2's cases.py edit is a no-op). `SERIES` `:1134-1144`
  (k28 8, p3 9, inversion 4, k32 10, nested_outer 11; recipe cutoff 1e-7, inversion 1e-5);
  `blast_series` `:1147`. `check_reference_pairs_clear_of_X` is `tests/helpers/laws.py:1647` (helpers
  `:1606-1645`), `reference_witness` `:1721` (table only, no scaling), `LAW_CASES_ONLY` `:3187`.
- `tests/invariants/test_law_partition.py`: 11 passed / 4 skipped (law ids are not -k selectable; the
  tier runs one test per layer x case).
- `helpers.fakes.synthetic_trellis` (`fakes.py:77-113`) puts EVERY registered crossing on its branches and
  passes no `resolved_ids`: a test needing a table-linked landing strictly inside a consecutive pair's
  interval must build the `Trellis` with `resolved_ids` excluding it (or cut).
- `stress/numerics/e12_reverify.py`: harness module (`EXPERIMENT :56`, `configs :70`, `run :98`,
  `harness.main :198`), `DEPTH :65` (k28 13, k32 13, p3 20, nested_outer 12, nested_inner 12, inversion
  10, heteroclinic 8), `SEPARATION :63`, two cutoffs 1e-7 / 1e-5; builds via `e11_blast_depth.build`
  (`e11:224`); results `stress/numerics/results/e12_reverify.jsonl` are COMMITTED (a97fcb6); its PDF
  numbers live in `report/make_numbers.py:579`.

## B2 outcome audit (2026-10-14, branch `refactor/pseudoneighbor-both-intervals` @ `9cd29d5`)

- `tests/unit/topology/test_pseudoneighbor.py` 36 tests. Module docstring: under both intervals every
  extra crossing on the window's stable branch punctures some stable interval unless it is an iterate
  of a window crossing, so unstable-only probes sit on `OTHER_STABLE = 1` (the no-inversion saddle's
  second stable branch). Seven base fixtures changed in DATA only (direct, iterate, walk_stops,
  heteroclinic_candidate_ignored, table_linked_deep_iterate (drift 1+2rtol -> 1-2rtol),
  table_linked_landing, unlinked_landing[below_c]); the literal base file on branch src: 7 failed /
  24 passed (audited). New: `test_scaled_iterate_punctures_stable_interval`,
  `test_table_landing_punctures_stable_interval` (builds `Trellis(...)` with `w` removed from every
  branch, as past a cut), `test_own_orbit_does_not_puncture_stable_interval` (local import of
  `_interval_is_clear`), `test_collision_on_another_unstable_manifold_still_punctures_stable_interval`
  `[own_orbit|other_unstable]` (TE, D7 mirrored; kills break m4).
- Gaps: unswapped m-range bounds survive the committed suite (no interval wider than one lap);
  `walk_stops_at_reference_window_end` and `heteroclinic_candidate_ignored` have no teeth for their
  named rule (identical on base).
- Suite on the branch: 1183 passed, 38 skipped, 5 deselected, 1 xfailed (only
  `orientation_reversing_case_builds`).

## B5 outcome audit (2026-10-14, branch `refactor/inertness-is-transience` @ `3fb310a`, base `7d3a294`)

- Suite on the branch (auditor run): 1193 passed, 38 skipped, 5 deselected, 1 xfailed (only C3).
  Test-count bookkeeping vs base 1193: -1 `test_resolve_inertness_fixed_point`, -2 PK-8 tests
  (`..._naming_an_active_class_warns_and_stays_out`, `..._mapping_onto_inert_classes_is_silent`:
  a self-cycle is active under transience), +2 `inert_classes` cases (`a->b, b->(), c->c`;
  `a->b, b unresolved`), +1 TE `test_the_session_alphabet_carries_the_dynamics_verdicts`.
- New laws (layer `symbolic`, `tests/helpers/laws.py`): `check_inert_is_transient` :3025 (own DFS,
  never `inert_classes`; also `entry.inert is (kind == "inert")`), `check_zeroth_hole_classes_have_empty_words`
  :3052 (itinerary, no symbol, inert -- in that order). `check_minimal_trellis_bridges` :2296 now
  asserts `inert_hole_bridge_ids == zeroth-hole bridges` and `hole_set - zeroth == unclassed |
  mapped | unmapped`. `anchor_bridge_class_leads_its_tangle` reads `dyn.classes[cls].kind`.
- On the base with the branch's `laws.py`: the empty-word half passes everywhere; transience and the
  kind half fail on exactly inversion_two_blasts (`g`) and heteroclinic (`d`, `bj`).
- GAP the TE found: the law tier is ORDER-DEPENDENT for "verdict before letters" --
  `anchor_bridge_class_leads_its_tangle` runs the dynamics on the shared case before
  `only_active_classes_lettered` reads `session.bridge_classes()`, so dropping `self.dynamics` from
  `Topology.classes` (break 5) passes the law tier; only the facade test and the TE's snapshot test
  catch it.
- `test_minimal_trellis.py::test_a_hole_bridge_without_a_registered_image_is_unmapped` (renamed,
  `p3_partitioned`): the fixture's `unmapped` is now EMPTY (both bridges carry zeroth holes), so the
  test strips `pair` from the holes (`dataclasses.replace`) and re-runs `minimal_trellis` to keep
  the unmapped path on real data.
- `tests/helpers/direct.py` `classes()` returns the unlettered alphabet (no `BridgeAlphabet`).
- Rubric-10 mutation `word is not None -> True` fails by `TypeError` (`SymbolicDynamics.py:516`),
  not on the assertion; the semantic mutation `word or ()` fails on the assertion.

## B7 analyst (2026-10-14, HEAD 3fb310a)

- Symbolic law tier + `test_symbolic_dynamics.py` + `test_partition_family.py`: 55 passed / 28 skipped (9.6 s).
- `check_k10_is_the_full_shift` docstring already cites A.8d / in-splitting (laws.py:3093-3096).
- A T law (`members_share_the_P1_itinerary`) would apply on k10, inversion_two_blasts, heteroclinic only
  (8 law cases NOT_APPLICABLE); the set-asides and any per-step series field live in `tests/cases.py`
  (`NOT_APPLICABLE` :1356, `blast_series` :1151 keeps no session, `SeriesStep` :1085).
- Mutating one child's `refined_rules` to a sub-word changes `transition_matrix()` live (the k10 `_Layout`
  goes `[[1,1],[1,1]] -> [[1,1],[1,0]]`), so an in-splitting law has teeth.

## B7 outcome audit (2026-10-14, branch `refactor/refined-level-defined` @ `a6a4f68`, base `3fb310a`, NOT merged)

- Suite (auditor's run): 1200 passed, 46 skipped, 5 deselected, 2 xfailed in 203 s. KNOWN_ISSUES = 2:
  C3 (orientation_reversing) and `("series_members_share_the_P1_itinerary", "k32_series")`, wp `B7b`
  (cases.py:1354), a strict xfail of `test_members_share_the_P1_itinerary_under_blasting`
  (regression/test_class_count_under_blasting.py:189). With `--runxfail` it fails at steps 6 and 7 on
  class `k A:p1@0.0/R#7 <-> A:p1@0.0/R#9` (= `{R_(0.0;8), R_(0.0;10)}`), member `(6, 7)`.
- New laws (tests/helpers/laws.py): helper `p1_itinerary_disagreements(session)` :2998 (non-loop members
  with a registered `image_chain`, stable-end chains skipped by `(manifold_b_key, crossing_sign)` equality
  of consecutive crossings, `ValueError` skip never fires; returns `(tested, skipped, disagreements)`),
  `check_members_share_the_P1_itinerary` :3083 (in `LAW_CASES_ONLY` :3380 because the deep-run test's
  `assert_non_vacuous` would fail on p3_15 / p3_6_blasts / inversion_period2_10, which test 0 classes),
  `check_refined_matrix_is_an_in_splitting` :3091 (all 11 law cases + deep runs), both in
  `LAYERS["symbolic"]`. `check_refined_children_inherit_word` :2980 kept (docstring only).
- `tests/cases.py`: `SeriesStep.p1_tested` / `p1_disagreements` (:1127-1128, filled in `blast_series`
  :1188 -- `cases.py` imports the helper from `helpers.laws`); `_NO_TWO_CHAINED_MEMBERS` set-asides
  (:1370, :1413) on k28_one_blast, k28_two_blasts, p3, nested, nested_inner_blasts, inversion,
  inversion_one_blast, inversion_period2.
- **Series KNOWN_ISSUES must use a `series_` id** in `_STANDALONE_ISSUES`
  (invariants/test_law_case_sanity.py:39-43); a layered law id on a series case makes
  `law_tier.set_aside_pairs` build a set-aside asking `law_case` for a non-law case (a sham xfail by
  KeyError). The wp regex now accepts `WP\d+|C\d|DEC-\d+|B\d+[a-z]?` (:100).
- Tested-class counts (T): k10 1, inversion_two_blasts 2, heteroclinic 7; series p1_tested k28
  [0,0,0,2,3,4,4,4,4], p3 [0,0,0,0,0,0,0,2,3,3], inversion [0,0,2,3,4], k32 [0,0,0,2,3,4,9,13,13,13,13],
  nested_outer [0,0,0,2,3,4,4,4,4,4,4,4] (141 total).
- Unit tests: `test_partition_family.py:150` (stub chain, three folds on the base branch ->
  `incomplete_footprints == [((1,2), branch, 3)]`, one INFO, no WARNING); `test_symbolic_dynamics.py`
  `_k28_dynamics` :332 (A.8d's k=2.8 shape by hand), `:370` (M and M_1), `:387` (the law fails on a
  sub-word). Teeth checked by the test engineer and reviewer (stable-end skip, sub-word, dropped column,
  dropped node, INFO->WARNING, B7b patch flips k32 to XPASS(strict)).
- When B7b lands, remove the KNOWN_ISSUE (strict: the series test XPASSes and fails until removed).
  [DONE on `refactor/base-chords-in-stable-order` @ `63477e3`: KNOWN_ISSUES = C3 only.]

## B7b analyst notes (2026-10-15, HEAD `a6a4f68`)

- `test_iterated_cut.py:92-110` (`test_foreign_crossings_pair_in_stable_order`) and `:112-123` give
  every OWN crossing the default `stable_cdist=1.0`, so their base-branch assertions cannot tell chain
  from stable order (both pass under a stable-order patch, verified with the scratch pytest plugin
  `scratchpad/b7b/b7b_patch_plugin.py`, `-p b7b_patch_plugin`). The module docstring :5-7 also states the
  chain-order rule. `test_partition_family.py:150` uses `stable_cdist=float(iid)` (stable order =
  chain order), also indifferent.
- With the stable-order patch the only failure in the affected files is the strict xfail
  `test_members_share_the_P1_itinerary_under_blasting[k32_series]` turning XPASS(strict); p1_tested per
  step `[0,0,0,2,3,4,9,13,13,13,13]`.
- `_STANDALONE_ISSUES` (`test_law_case_sanity.py:36-43`) is an ALLOW-list, not tied to a live entry:
  four of its six ids already have no KNOWN_ISSUES entry. Dropping `series_members_share_the_P1_itinerary`
  is safe either way; the guard (`:91-95`) still rejects any unknown id; the wp regex is at `:100`.
- `test_foreign_chords.py` reads the stretch midpoint as `stretch[len // 2]` of `arc_polyline` (real
  nodes only, the crossings excluded) and skips a stretch with no node; on k32 the base-branch gaps
  `(11,10)` has 2 nodes and `(79,80)` / `(81,78)` have 0, so a base-branch sibling must close each gap
  at its crossing coordinates and take an arc-length midpoint (as `b7b/probe.py` does).
- [B7b outcome: the toothless chain-order half of `test_foreign_crossings_pair_in_stable_order` (the
  tie `(12,"chord") in empty[11]`) is deleted; the k32 hand-built shape is the new unit pin
  `test_base_crossings_pair_in_stable_order`. The foreign test now uses a function-scoped fixture over a
  module `functools.cache` of `cases.BUILDERS`; `_lobe_image` moved to `helpers/laws.py`. Its own
  node-midpoint reading (skips node-less gaps) is unchanged.]

## B7b outcome audit (2026-10-15, branch `refactor/base-chords-in-stable-order` @ `63477e3`, base `a6a4f68`, NOT merged)

- Suite (auditor's run on the branch, 207 s; test engineer and reviewer agree): 1211 passed, 47
  skipped, 5 deselected, 1 xfailed (C3 only). Base: 1200 / 46 / 5 / 2 xfailed.
  Collected ids: +11, none removed.
- New reading `helpers.laws.base_chord_misreadings(session) -> (lobes checked, texts)` (laws.py:3111,
  beside `_lobe_image` :3085 and `_arclength_midpoint` :3102): every stretched hole whose CONSECUTIVE
  image chain crosses its base branch at 2+ interior crossings (a non-consecutive chain is skipped
  silently); gaps anchor -> chain crossings in stable order -> branch tail, each closed at the crossings'
  coordinates, read at its arc-length midpoint by `point_in_polygon(tol=0.0)`. Clauses: (a) the base is
  one gap, (b) alternation with the base inside, (c) every inside gap is an empty stretch of the full
  kernel, (d) every base/chord the kernel declares for THIS hole bridge alone (`SimpleNamespace` minimal)
  is one inside gap. Imports the private `_hole_stretches`. Reads polylines only, never the map.
- Where it runs: `SeriesStep.chords_checked` / `chord_misreadings` (cases.py:1110-1133, filled in
  `blast_series` :1195); `test_base_chords_are_inside_the_lobe_image_under_blasting`
  (regression/test_class_count_under_blasting.py:204, id `series_base_chords_inside_the_lobe_image`,
  prints per-step counts: k28 [0,0,1,...], p3 all 0, inversion [1]*5, k32 [0,0,1,1,1,1,2,2,2,2,2],
  nested_outer [0,0,1,...]); `test_base_chords_are_inside_the_lobe_image` (regression/test_foreign_chords.py:99)
  over k10, inversion, inversion_period2, heteroclinic (checked 1 each) and p3_6_blasts (`NO_BASE_LOBE`
  :37, asserts `not checked` then a visible skip). Total 36 series + 4 build readings = 40 by the tests;
  the helper reads all 47 analyst lobes clean (11 builds), 7 of them on builds the test does not list.
- Teeth: on base the series test fails on k32 steps 2-10 (clause (c) at 2-7, (d) at 2-10), the unit pin
  fails (`(10,'chord') in [(12,'chord')]`, auditor re-ran it), T fails at 6-7. Outside the k32 lobe every
  base lobe has one inside gap that is ALSO another hole's stretch, so dropping chords (breaks A, B) goes
  unnoticed there; a wrongly added chord (break C) is caught everywhere by (d). The reviewer's reversed
  pairing break turned 8 tests red.


## B6b `word-finished` analyst notes (2026-10-15, HEAD `63477e3`)

- HEAD suite: 1211 passed, 47 skipped, 5 deselected, 1 xfailed (C3) in 214 s.
- Scratch plugin `scratchpad/b6b/b6b_reduce_plugin.py` (`-p b6b_reduce_plugin`, `PYTHONPATH` = the
  b6b dir; `B6B_RULE=occ` for the same-occurrence rule) patches a free reduction into
  `SD.translate_itinerary` before test modules import it. Class-level: 3 failures, all hand-built unit
  fixtures (`test_loops_dropped_from_word_but_recorded` pins `["a","a^-1"]`;
  `test_evidence_disagreement_names_a_stable_end[stable end|walk]` -- registered `a a^-1`, position 1
  -> 0), 3 hits total; every law, series and T test passes. Same-occurrence: 0 hits, suite green.
- WalkSearch fine-status pins: `test_dual_walk.py` 29 assertion lines (unique 12, ambiguous 9,
  unreachable 2, one each trivial / no_start_node / no_goal_node / no_admissible_walk / truncated /
  span_meets_no_node) plus `truncated` flag pins :254-260 (status `ambiguous` WITH the flag) and :450;
  `test_heteroclinic.py` :76, :94; `laws.py` :2775, :2817, :2932; `test_symbolic_dynamics.py` :754
  (`first.search.truncated and first.search.walk is not None`). `cd.status` pins (:443, :445) are the
  ClassDynamics four values, untouched by a WalkSearch collapse.
- Logger: evidence WARNING/INFO currently on `tanglepack.topology.SymbolicDynamics`;
  `test_symbolic_dynamics.py:673` counts WARNINGs on that logger name (`assert_logged`), so a move to
  the diagnostics logger changes it; `test_evidence_diagnosis.py` filters by message (believed to pass
  unchanged).

## B6b `word-finished` outcome (2026-10-15, branch `refactor/word-finished` @ `65edec4`, base `63477e3`; NOT merged)

- **Suite.** Branch: 1223 passed / 47 skipped / 5 deselected / 1 xfailed (C3), no XPASS (TE and
  reviewer, re-run by the auditor). Base: 1211. The +12 are new tests.
- **New tests.**
  - In `test_dual_walk.py`:
    - `test_every_failure_is_unresolved_with_its_reason`: 6 ids on the `tie` graph, all fail on base;
    - `test_a_cap_cutting_a_tie_short_keeps_its_candidate_unresolved`: fails on base, where the
      capped tie reads `unique` with `truncated=True`, NOT `ambiguous` as the plan's rubric said.
  - In `test_symbolic_dynamics.py`:
    - `test_free_reduction_cancels_adjacent_inverses_only`;
    - `test_a_walk_reading_c_cinv_is_the_empty_word`;
    - `test_a_class_named_only_inside_a_cancelled_pair_is_not_added`;
    - the TE's `test_the_registered_word_is_reduced_before_it_is_compared` (fails on base) and
      `test_an_ambiguous_search_marks_its_class_and_is_no_wall` (a pin, passes on base).
- **Pins rewritten, never weakened.**
  - `status in ("unique","trivial")` -> `== "resolved"`. This is stricter: a capped search used to
    pass and now fails.
  - `not search.truncated` -> `reason is None`.
  - The fine failures -> `status == "unresolved" and reason == <value>`.
  - Two tests were renamed: `test_zero_crossing_walks_stay_on_the_start_face` and
    `..._dedupe_to_one_walk`.
  - The loops test still pins the unreduced `translate_itinerary` (`["a","a^-1"]`) and adds
    `free_reduce(...) == []`.
  - The evidence fixture registered word is `a a` (it was `a a^-1`). Position 1 and the
    occurrences are kept.
- **Patching traps.**
  - `_drive` now patches `diagnostics._image_chain` / `diagnostics.trellis_itinerary`. A patch on
    `SD.*` no longer reaches the evidence path: the TE broke it on purpose and got 7 FAIL.
  - `_drive`'s fake naming needs `homotopy` / `iterated` attributes.
  - The evidence logger is `tanglepack.topology.diagnostics`. `test_evidence_diagnosis.py` filters
    WARNINGs by that name.
- **Break-its (TE, `scratchpad/te_b6b/breaks*.txt`).**
  - B1-B11 each fail at least one test after the TE's commit.
  - B12 (`plotting._class_status` returns `-`) fails nothing: the label is untested by design (it
    is wording).
- **Corpus coverage gap.** The 63 items hold only `resolved` searches (553 = 402 zero-crossing +
  151 unique) and no search at all for 349 classes. The `ambiguous`, failure and truncated rows of
  the status mapping and any free-reduction cancellation are reached by hand-built unit fixtures
  only.
- **Probe scripts.** `scratchpad/coder/probe.py` / `run.sh` / `compare.py` (the TE's copies, with a
  naive-reduction oracle, are in `scratchpad/te_b6b/`). The auditor's spot check is
  `scratchpad/auditor/spot.py REPO case...`.
- To run branch tests against base code, use a `git archive` copy with the tests laid over it
  (pyproject's `pythonpath` puts the rootdir `src` first; a `PYTHONPATH` override does not work).

## B9 analyst (2026-10-15, `HEAD` `517aa01`)

- Suite: 1215 passed / 47 skipped / 5 deselected / 1 xfailed (C3), no XPASS, 204 s.
- Figure determinism: the symbolic-itinerary figure scripts regenerate byte-identical PNGs run to run
  (42 PNGs: k10, k28_blast1/2, p3, nested_p1_p3, heteroclinic_p1_p3, inversion; exported HEAD copy,
  Agg). k10 / k28 equal the COMMITTED PNGs; the committed p3 / nested / heteroclinic / inversion PNGs
  are STALE (HEAD draws different ones), so a base-vs-branch figure check must compare two
  regenerations, never the committed files. Script times 2-14 s.
- Tests reading session plot delegates / cartoon geometry / fan-outs: `test_plot_smoke.py:59-71,
  74-80, 277-308, 311-320`; `test_plot_topology.py:381-405, 408-431`;
  `test_session_bridge_classes.py:36-38 (_repartition_at), 151-168`; `test_session_fanouts.py`
  (all but :86-95); `test_session_caches.py:135-137`; `test_nested_own_blast.py:78-80`;
  `test_arrangement.py:316-318`; `test_partition_elements.py:162,186`. The law tier
  (`laws.py`, `cases.py`, `invariants/`, `conftest.py`) reads none of them (it reads
  `resonance_zones`, `add_resonance_zones`, `row_of_end`, `element_of_intersection`).

## B9a outcome (2026-10-15, branch `refactor/session-three-concerns` @ `5f0db5f`, base `517aa01`)

- 1215 -> 1208 passed (47 skipped / 5 deselected / 1 xfailed C3 unchanged). Deleted ids:
  `plotting/test_plot_smoke.py::test_session_plot_fanouts_draw_every_fixed_point`,
  `::test_session_plot_pseudoneighbors_computes_on_demand`, `::test_session_plot_delegates_return_their_types`,
  `facade/test_session_bridge_classes.py::test_a_pip_trimmed_away_is_not_replaced_by_plotting[plot_strong_pip|plot_strong_pip_candidates]`,
  `facade/test_session_fanouts.py::test_describe_fanouts_mention_every_fixed_point`,
  `unit/topology/test_partition_elements.py::test_session_closes_the_cross_trellis_gap`. Renamed:
  `test_session_fanouts.py` -> `test_classify_strong_pips_covers_every_fixed_point` /
  `test_classify_strong_pips_of_one_fixed_point_is_its_list`; `test_plot_smoke.py` ->
  `test_cartoon_honours_its_toggles`.
- `test_session_caches._repartition_by_hand` now asserts the repartitioned trellis HOLDS a partition
  (before, `clear_results()` alone satisfied `!= partition`; TE commit `5f0db5f`).
- Gap: `tanglepack.loom.plot_resonance_zones` has no test (it had none as a session method either).
- `test_plot_topology.py::_circle_cartoon` is the explicit circle-cartoon recipe the tests use.

## B9b analyst (2026-10-15, `HEAD` `5f0db5f`)

- Suite on an exported HEAD copy: 1208 passed / 47 skipped / 5 deselected / 1 xfailed (C3), 199 s.
- `test_partition_elements.py` (168 lines, 5 tests): :42 checks every WORKBENCH bridge end whose
  `manifold_b_key` is the result's branch resolves on that result's side (`element_for` ==
  `element_of_intersection`) plus `set(elements_at_bridge) <= expected` (:74); :78 is the field's own
  test; :96 both ends of a branch-wrapping bridge, both sides; :119 a bridge dropped from
  `trellis.bridges` before `partition_stable_manifold` still resolves both ends, both sides (and is
  absent from the field, :155); :161 `element_for(..., "third", ...)` -> ValueError.
- A full-strength rewrite (scratchpad `b9b/cand_b_tests.diff`, 149 lines, 4 tests) reads
  `HomotopyPartition.from_results(trellis.stable_partitions, trellis=trellis).owner_of_intersection(
  id, side, trellis.stable_key_at(id, trellis.unstable_key_of(bid)))` and compares whole
  `ElementRef`s (branch, side, id) -- stronger than the old `element_id` compare. Teeth verified:
  `owner_of_intersection` returning the neighbour element fails the three both-side tests;
  `endpoint_index` accepting "third" fails the ValueError test.
- `tanglepack.loom.plot_resonance_zones` smoke: `k10_partitioned` has NO resonance zone;
  `k28_partitioned` has one (draft passes in 0.23 s, scratchpad `b9b/smoke_draft.py`).
- Probes: scratchpad `b9b/probe/` (element-at-end agreement, runs on base and branch, per-item
  readings sha256) and `b9b/probe63/` (the B10 63-item outcome probe copy).

## B9b outcome (2026-10-15, branch `refactor/one-element-at-end` @ `94fca6b`)

- Suite (auditor run): 1207 passed / 47 skipped / 5 deselected / 1 xfailed (C3), 209 s; collected
  1256 -> 1255.
- `test_partition_elements.py` is 149 lines, 4 tests, all on `p3_partitioned` /
  `henon_p3_session`: `test_every_workbench_bridge_touching_the_branch_resolves`,
  `test_both_ends_of_a_branch_wrapping_bridge_resolve`,
  `test_a_bridge_missing_from_the_snapshot_resolves` (both sides, both ends, via the helper `_owner`
  = `owner_of_intersection` on the `stable_key_at` branch, against `_expected` read off the end's
  own `manifold_b_key` result; whole `ElementRef`s), and
  `test_owner_at_end_rejects_an_unknown_endpoint_name`. Deleted:
  `test_elements_at_bridge_entries_agree_with_the_ownership_table` (with the field; recover from
  `5f0db5f`). Renamed from `test_element_for_*`.
- What this file does NOT catch (test engineer's planted breaks): `stable_key_at` ignoring C7 and
  `owner_at_end` reading the opposite row both pass here; the law tier (`test_classes_laws`,
  `test_symbolic_laws`, `test_partition_laws`, ...) catches both. It never tested the row: the
  old `element_for` took the side as an argument.
- The class rule's error paths: `test_bridge_class.py:311` (missing partition -> ValueError naming
  the branch, now raised by `owner_of_intersection`) and `:332` (duplicate partitions, still
  `_index_partitions`) unchanged and passing. Nothing pins the unowned-crossing path (now ValueError).
- Still a gap: `tanglepack.loom.plot_resonance_zones` has no test (B12; draft `b9b/smoke_draft.py`).

## B12 analyst (2026-10-15, `HEAD` `94fca6b`, read-only)

- Suite at HEAD re-run: 1207 passed / 47 skipped / 5 deselected / 1 xfailed (C3), 203 s.
- `tests/helpers/laws.py`: 106 laws (+ the helper `law_id`, which the request counted as the
  107th): LAYERS 103 + MUTATING 3. **Deleting a law changes NO law-tier test id**: the tier runs
  ONE test per (layer, case) (`helpers.law_tier.run_layer`); only the law's set-aside ids
  (`test_<layer>_set_aside[<law>-<case>]`, one per `KNOWN_ISSUES` / `NOT_APPLICABLE` pair) go.
- Classification (conservative; full table in the B12 findings artifact): 10 TAUTOLOGY candidates
  that re-implement or restate the producing code line for line --
  `bridge_partial_iff_no_id` (both properties derived from the same two end ids, `Bridge.py:161-185`),
  `manifold_iterate_law` (`c_iter = stretch*c` is the assignment at `ManifoldMachine.py:355`),
  `openings_follow_signs` (= `StablePartition._hole_openings:684-704`; `partition_covers_branch`'s
  `_sign_openings` clause re-derives the same halves), `openings_on_own_bridge_row` (non-anchor rows
  via `row_of_end`, the function `_hole_openings` calls), `class_table_order` (= `table_key`,
  `BridgeClass.py:642-648`), `classes_do_not_mix_tangles` (= `tangle_of`, `BridgeClass.py:631-640`;
  one residual "every tangle has a class" clause), `dual_unified_is_pip_segment` (= `_pip_segment`
  + `_meets_fundamental_segment`; the oracle is `dual_unified_is_where_f_of_s_stops`),
  `dual_bipartite_degree` (restates `DualGraph.graph()`; degrees implied by `dual_node_structure`),
  `itinerary_pairs_are_classes` (= `translate_itinerary`, and compares to the FREE-REDUCED
  `cd.symbols`: it would fail falsely the first time a cancellation fires), and
  `refined_children_inherit_word` (= `refined_rules[child] = list(refined_word)`,
  `SymbolicDynamics.py:722`; T and the in-splitting law are the oracles since B7). PROBE-ONLY:
  `iterate_table_complete`'s three `unresolved_links` clauses (B8). Everything else ORACLE or
  MIXED (keep). Only `refined_children_inherit_word` has set-aside ids (4 NOT_APPLICABLE:
  k28_one_blast, p3, inversion, inversion_period2; `cases.py:1395-1398`). Two direct callers outside
  the tier: `tests/unit/loom/test_blast_budget.py:39` (`check_bridge_partial_iff_no_id`) and
  `tests/regression/test_coarse_blasts_rows_by_signs.py:52` (`check_openings_follow_signs`).
- `tests/test_higher_period_cartoon.py` no longer exists (the 2026-10-02 bullet above is stale; its
  tests were absorbed into `tests/plotting/test_plot_topology.py`); CLAUDE.md:90 still names it.
- `plot_resonance_zones` smoke draft (`scratchpad/b9b/smoke_draft.py`) re-run 2026-10-15 from a
  scratch conftest that loads `tests/conftest.py`'s `k28_partitioned`: 1 passed, 0.23 s.
- Sphinx: `env/bin/sphinx-build` 8.1.3 is installed, but `conf.py` needs `myst_parser` and `furo`
  (neither installed) and `_static/` does not exist; with those overridden
  (`-D extensions=sphinx.ext.autodoc,sphinx.ext.napoleon -D html_theme=alabaster`) the build
  finishes with 1,206 WARNINGs (duplicate cross-reference targets from the package re-exports,
  docutils indentation). A `-W` build is not a feasible rubric; an automodule-list check is.

## B12b analyst (2026-10-16, `HEAD` `a550a5d`, read-only; scratch copies under scratchpad `b12b/`)

- Baseline at HEAD: 1260 collected (5 deselected), 1207 passed / 47 skipped / 1 xfailed (C3), 200 s.
- The nine TAUTOLOGY laws (B12a Q2) deleted in a scratch copy (`b12b/delete_laws.py`): 106 -> 97
  `check_*` functions, laws.py -209 lines (9 functions + trailing blanks, 9 `LAYERS` lines, and the
  two imports only they used: `assert_iterate_relation` (helpers.invariants, still used by two unit
  tests) and `class_sort_key` (BridgeClass)). The `NOT_APPLICABLE` "four entries" are ONE dict
  comprehension (`cases.py:1395-1399`): the edit is `for law in ("member_matching_consistent",)`.
  Collected ids: exactly the 4 `test_symbolic_set_aside[refined_children_inherit_word-{inversion,
  inversion_period2,k28_one_blast,p3}]` go (they were skips). With the smoke test: 1257 collected,
  1208 passed / 43 skipped / 1 xfailed (C3). Rubric script `b12b/laws_rubric.py` (deleted / LAYERS /
  imports / no surviving body changed by AST dump).
- `DualGraph.graph()` has exactly one caller in src/tests/scripts/stress/notebooks:
  `laws.py:2716` (`check_dual_bipartite_degree`).
- `test_coarse_blasts_rows_by_signs.py:52`: replacing `check_openings_follow_signs` by
  `check_partition_covers_branch` passes (count 6 on `p3_coarse_blasts`); `_sign_openings`
  (`laws.py:1086`) uses the same D6 formula but checks the PARTITION, not `hole.openings`.
- Smoke test for `loom.plot_resonance_zones` on `k28_partitioned`: 1 passed, 0.22 s setup.
- PNG reproducibility: matplotlib SVGs carry `<dc:date>` and random ids (`svg.hashsalt` unset);
  pickles differ byte-wise on every run but re-render identically. With `SOURCE_DATE_EPOCH` and an
  rc `svg.hashsalt` the SVGs reproduce; PNGs reproduce always (same interpreter).
- 63-item probe (`b10/probe/run.sh REPO OUTNAME`, writes `b10/probe/OUTNAME/`): 34 s for two trees.

## B12b outcome (2026-10-16, branch `refactor/laws-and-loose-ends` @ `e163f7a`, base `a550a5d`; auditor)

- Suite at head (auditor's own run): 1208 passed / 43 skipped / 5 deselected / 1 xfailed in 195 s;
  `-rxX` lists only C3, no XPASS. Collected 1257 (1252 + 5 deselected); diff vs base = exactly the
  four `test_symbolic_set_aside[refined_children_inherit_word-*]` out and
  `tests/plotting/test_plot_smoke.py::test_plot_resonance_zones_shades_every_zone` in.
- `tests/helpers/laws.py`: 97 laws. Gone (TAUTOLOGY, each restated an implementation):
  `manifold_iterate_law`, `bridge_partial_iff_no_id`, `openings_on_own_bridge_row`,
  `openings_follow_signs`, `class_table_order`, `classes_do_not_mix_tangles`,
  `dual_unified_is_pip_segment`, `itinerary_pairs_are_classes` (= `translate_itinerary`'s loop,
  `SymbolicDynamics.py:511-518`), `refined_children_inherit_word`. `check_dual_bipartite_degree`
  stays and is still the ONE caller of `DualGraph.graph()`. `NOT_APPLICABLE` (`tests/cases.py`)
  is now `member_matching_consistent` only. `test_blast_budget.py` lost its call to the
  partial-iff law; `test_coarse_blasts_rows_by_signs.py:52` calls `check_partition_covers_branch`.
- Smoke test `test_plot_resonance_zones_shades_every_zone` (`k28_partitioned`): pins
  `len(plot_resonance_zones(zones)) == len(zones)` and `[] -> []` only. Weak: a stub returning
  `[None] * len(zones)` stays green (TE break); breaking the `poly.size == 0` branch turns it red.
- Rubric tools that worked: `scratchpad/b12b/laws_rubric.py` (AST dump per function, LAYERS order,
  imports), `scratchpad/b12/stripped_ast.py` (docstring-stripped AST per src file), the 63-item probe
  `b10/probe/run.sh`. One that did not: `b12b/png/svgnorm.py` misses the salted SVG marker ids
  (6/18 equal; 18/18 with `te/svgnorm2.py`).

## `blast-steps-compose` analyst (2026-10-21, `HEAD` `c1ca6d0`)

- `too_close` per single-iteration call (base, ZONE_MIN_SEPARATION 1e-5 unless noted): every
  LAW_CASE 0 (no law case drops a child); `k28_deep_blast_sliver` 0x8 then 24 (last blast);
  k28_series 0x7, 2; p3_series / inversion_series 0; k32_series 0x7, 5, 36, 97; nested_outer_series
  (1e-4) 0x6, 5, 22, 58, 119, 225. Drops before the last blast leak (re-iterated next call), so only
  k32 steps 9-10 and nested_outer steps 8-11 carry leaked bridges at HEAD.
- Running the suite against a scratch copy of `src`: `pytest -o pythonpath="<copy> tests"`
  (pyproject's `pythonpath = ["src", "tests"]` overrides `PYTHONPATH`). 1208 passed, 43 skipped,
  5 deselected, 1 xfailed (C3) at `c1ca6d0`, ~3.5 min.
- Composition-test candidates (n x `blast_zone(...,1)` vs `blast_zone(...,n)`): `henon_p3_session`
  fp1 sep 0.02 `forget_strong_pip` n=4 (0.4 s), inversion zone sep 0.02 n=4 (0.3 s): fail at HEAD,
  pass with drops remembered. k28 zone sep 1e-5 n=11 (~22 s) also fails with memory alone (the
  reference-cloud leak) and passes with either cloud rule.
