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
