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
