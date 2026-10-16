# Test evidence — B10 `dead-api` (test engineer, round 1, 2026-10-15)

## In short

- **Base** `65edec4` (detached worktree `/home/dezhu/Development/TanglePack/tangle-pack-dead-api-base`). **Branch** `refactor/dead-api` at head `517aa01` (worktree `/home/dezhu/Development/TanglePack/tangle-pack-dead-api`).
- **The two trees compute the same thing.** The 63-item probe output is **byte-identical** on base and branch: md5 `bca69acaa8b8bc3be4438463ff174eef` for both. It covers every class's word, both transition matrices, the class counts and `is_reliable`.
- **The tests that disappeared are exactly the intended ones.** The collected-id diff is the five graph tests, the three `from_branch` parametrisations and one rename. Nothing else changed.
- **The suite passes at every one of the eight commits.** I ran all eight myself, in scratch worktrees for C1–C7. Base gives 1223/47/5/1, and C4 onwards gives 1215/47/5/1. The only xfail is C3 (INV-11), and there is no XPASS anywhere.
- **I added no test and made no commit.** Every rewritten test fails under a break of the code it guards (8 breaks, all caught), so none needed strengthening.
- **No rubric item failed.** There is one note: `ruff F401` flags the 56 names that now sit in `topology/__init__.py`'s import block but not in `__all__`. This is the expected effect of Q8, which keeps every import line and leaves the import block to B12. No other edited file has an F401/F841/F811/F821 finding.

Commands use `SUITE` = `MPLBACKEND=Agg PYTHONPATH=src ../tangle-pack/env/bin/python -m pytest -q -rxX -p no:cacheprovider`, run from inside the named tree. Scratch-worktree runs used the absolute path of `env/bin/python`. Scratch root: `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b10/`. Below, `te/` and `probe/` are relative to it.

## Rubric, item by item

| # | Item | Ran by me | Result |
|---|---|---|---|
| 1 | Suite on the branch | yes | **1215 passed, 47 skipped, 5 deselected, 1 xfailed in 210.23s**. Only xfail: `test_orientation_reversing_case_builds[orientation_reversing]` (C3: INV-11). No XPASS. Exit 0 (`te/suite_head_517aa01.txt`). Base: **1223 passed, 47 skipped, 5 deselected, 1 xfailed**, same xfail, exit 0 (`te/suite_base_65edec4.txt`). |
| 2 | Exactly the intended tests gone | yes | See "Collected ids" below. 1271 → 1263 collected (5 deselected on both). Gone: the 5 graph tests and the 3 `from_branch` parametrisations. Renamed: 1. Nothing new. |
| 3 | No law-tier pin moved | yes | `git diff --stat 65edec4 517aa01 -- tests/helpers/laws.py tests/cases.py tests/invariants/` is empty. `pytest tests/invariants` gives **133 passed, 46 skipped, 1 xfailed** on base (18.75s) and **133 passed, 46 skipped, 1 xfailed** on the branch (18.72s), with the same C3 xfail (`te/inv_base.txt`, `te/inv_branch.txt`). |
| 4 | Diff inside the file list | yes | `git diff --stat 65edec4 517aa01`: **27 files, +63/−1097**. The 27 names (`--name-only`) are exactly rows 1–27 of the plan's table. Nothing under `CLAUDE.md`, `docs/`, `.refactor/`, `artifacts/`, `notebooks/` or `figures/` changed. The main tree's tracked `git status` is the same as at session start (only this artifact is new). The branch and base worktrees show 0 changes after all my runs. |
| 5 | Eight commits, each alone green | yes, all eight | Subjects are `Dead API: <name> (B10)`, one per row C1–C8. Each commit's file set matches its row, and all 8 messages carry `Co-Authored-By`. Suite results in the table below. |
| 6 | Net lines over src+tests | yes | `git diff --shortstat 65edec4 22feb61 -- src tests` (C7): **18 files, +35/−992 = −957** (target ≤ −950). `git diff --shortstat 65edec4 517aa01 -- src tests` (head): **19 files, +35/−1048 = −1013** (target ≤ −1000). |
| 7 | No deleted name in `*.py` | yes | The plan's grep over `src tests scripts stress` on the branch returns nothing (exit 1). |
| 8 | No graph in numerics | yes | `grep -rn "networkx\|_graph" src/tanglepack/numerics` returns nothing (exit 1). `src/tanglepack/numerics/graphviz.py` does not exist. |
| 9 | Imports from outside the repo | yes | From `/tmp`, `PYTHONPATH=<worktree>/src`, a fresh interpreter per statement, each asserting that `tanglepack.__file__` lies under that worktree. All five statements succeed on both trees. `print('networkx' in sys.modules)` after `import tanglepack` prints **True** on base and **False** on the branch. |
| 10 | `topology.__all__` = the 48 | yes | Branch: `len(__all__) == 48`, 48 unique, every name resolves (`hasattr`). The plan's C8 list typed out as 48 names equals `set(__all__)` with an empty symmetric difference, so all 15 (a) names are present. Base has 104. The `topology/__init__.py` diff is `0 +, 56 −`, and every removed line is an `__all__` string. No import line changed. |
| 11 | Refinement test keeps its teeth | yes | Read: the test sets `area_cutoff = 0.0`, calls `_refine_layer([(p0, p1)], manifold, viewer, None, set())`, asserts `after is new_point`, and asserts `lo < float(new_point.cdist) < hi`, which is strict with no tolerance. Mutations M1–M3 below all fail it. |
| 12 | Classify test keeps its teeth | yes | Read: the `test_resonance_zone.py` diff is one line, `classified = {b: session.classify_bridge(b) for b in session.workbench.bridges}`. The midpoint-in-zone and `classify_bridge` assertions are untouched. Break B3 is caught. |
| 13 | Iterate links still asserted | yes | Read: the old `:85` became `assert registry.iterate_table[z, 1] == y`. The old `:147` (now `:146`) became `assert registry.iterate_table[x, 1] == y`. The dropped `not graph().has_edge` lines are covered by the twin asserts that stay: `table[x, 1] is None and table[z, 1] is None and table[y, -1] is None`, and `iterate_table[x, 1] is None and iterate_table[y, -1] is None`. Breaks B1 and B2 are caught. |
| 14 | Blast tests keep a stop reason | yes | Read: `assert second.stop_reason == "point_budget"` stays. The no-budget test keeps `stop_reason in {"max_iterations", "frontier_empty"}` and `== "frontier_empty"` for the empty filter. Only the `terminated_early` lines (a restatement of the line above, and a tautology) went. Break B4 is caught. |
| 15 | 63 probe items agree | yes | **0 differences** on all 63 items (16 builds + 5 series = 47 series points). See "Probe" below. |
| 16 | Rewritten scripts print the same | yes, independently | The 7 scripts were run on both trees and the normalised outputs diffed. The only differences are the three allowed ones; see "Scripts" below. All exit 0. |
| 17 | Dead script dead the same way | yes | `henon_blast_period_3.py` raises on both trees at `line 95`: `ValueError: Intersection 12 is not a strong-pip candidate. Candidates: [31, 32, 33]. …`, exit 1. `git diff --numstat` for it is `1 1`. **The script is dead on base and branch alike.** |
| 18 | Long scripts compile | yes | `py_compile` (doraise) passes on the branch for `henon_symbolic_itineraries{,_inversion,_nested,_period3}.py`, `show_figure.py` and `stress/numerics/e12_reverify.py`. It also passes for the two untracked `henon_period3_{k,zone}_sweep.py`, compiled from the main tree because they are not on the branch. None was run, per the plan. |

## Suite at every commit (rubric 1, 5)

Each commit was checked out as a detached scratch worktree (`git worktree add --detach te/wt_<hash> <hash>`) and `SUITE` was run in it. The scratch worktrees were removed afterwards (`git worktree remove`, then `prune`), and none was left dirty.

| Commit | Row | Result (all exit 0, only xfail C3, no XPASS) | Log |
|---|---|---|---|
| `65edec4` base | — | 1223 passed, 47 skipped, 5 deselected, 1 xfailed | `te/suite_base_65edec4.txt` |
| `1867f41` | C1 | 1223 passed, 47 skipped, 5 deselected, 1 xfailed | `te/suite_C1_1867f41.txt` |
| `88d5a75` | C2 | 1223 passed, 47 skipped, 5 deselected, 1 xfailed | `te/suite_C2_88d5a75.txt` |
| `217920a` | C3 | 1223 passed, 47 skipped, 5 deselected, 1 xfailed | `te/suite_217920a.txt` |
| `7eaf64d` | C4 | 1218 passed, 47 skipped, 5 deselected, 1 xfailed (−5 graph tests) | `te/suite_7eaf64d.txt` |
| `3a9db04` | C5 | 1215 passed, 47 skipped, 5 deselected, 1 xfailed (−3 `from_branch`) | `te/suite_3a9db04.txt` |
| `1348e8c` | C6 | 1215 passed, 47 skipped, 5 deselected, 1 xfailed | `te/suite_C6_1348e8c.txt` |
| `22feb61` | C7 | 1215 passed, 47 skipped, 5 deselected, 1 xfailed | `te/suite_22feb61.txt` |
| `517aa01` | C8 head | 1215 passed, 47 skipped, 5 deselected, 1 xfailed | `te/suite_head_517aa01.txt` |

## Collected ids (rubric 2)

Command, in each worktree: `MPLBACKEND=Agg PYTHONPATH=src ../tangle-pack/env/bin/python -m pytest --collect-only -q -p no:cacheprovider`. The `::` lines were sorted and diffed (`te/ids_base.txt`, `te/ids_branch.txt`).

- base: `1271/1276 tests collected (5 deselected)`
- branch: `1263/1268 tests collected (5 deselected)`

`diff ids_base.txt ids_branch.txt`, verbatim:

```
< tests/unit/numerics/test_intersection_registry_fixes.py::test_from_branch_holds_the_anchor_on_every_branch[None]
< tests/unit/numerics/test_intersection_registry_fixes.py::test_from_branch_holds_the_anchor_on_every_branch[stable]
< tests/unit/numerics/test_intersection_registry_fixes.py::test_from_branch_holds_the_anchor_on_every_branch[unstable]
< tests/unit/numerics/test_iterate_inference.py::test_remove_iterate_forgets_the_link_its_twin_and_its_edge
> tests/unit/numerics/test_iterate_inference.py::test_remove_iterate_forgets_the_link_and_its_twin
< tests/unit/numerics/test_map_step_and_graph.py::test_build_intersection_graph_does_not_mutate_the_registry_graph
< tests/unit/numerics/test_map_step_and_graph.py::test_registry_graph_falls_back_to_per_branch_unstable_order
< tests/unit/numerics/test_map_step_and_graph.py::test_registry_graph_takes_the_bridges_it_is_given
< tests/unit/numerics/test_map_step_and_graph.py::test_stable_edges_join_consecutive_crossings_on_one_branch
< tests/unit/numerics/test_map_step_and_graph.py::test_unstable_edges_are_exactly_the_bridges
```

This is exactly the expected set: 8 gone, 1 renamed, nothing new.

## Probe (rubric 15)

**Harness.** I copied `b6b/{run.sh,probe.py,summ.py}` to `probe/` and left the b6b copies untouched.

- `run.sh` takes the worktree and an output dir as arguments (`REPO=$1`, `O=$D/$2`, `cd $REPO`) and exports `REPO`. `probe.py:4` reads `REPO = os.environ["REPO"]`.
- **Adaptation (deliberate).** The B6b probe measured free-reduction counts. It also reads `WalkSearch.truncated`, which B6b itself removed, so at `65edec4` it would raise an `AttributeError`. I kept its build/series driver unchanged: the same 16 builds, the same 5 series, and the same `cases._blast(..., min_separation=series.min_separation)` loop. I replaced its `probe()` with a dump of every outcome rubric 15 names, written as one sorted-key JSON line per item:
  - `classes` (closed table) and `seed` (`len(session.bridge_classes())`) counts, plus `kinds`
  - `rules` and `refined_rules` (every word, unrefined and refined)
  - `transition_matrix(refined=False)` and `(refined=True)`, names and entries
  - per class: letter, label, kind, status, word, `unresolved_reason`, `ambiguous`, `verified`, walk-search status and reason
  - `len(unmatched_members)` and `is_reliable`
- I confirmed with `tanglepack.__file__` / `cases.__file__` that each run loads its own worktree's `src/` and `tests/cases.py` (the branch shows `len(__all__)` 48, base 104).
- `summ.py` was copied but not used; it summarises the old reduction fields.

**Runs.** `bash probe/run.sh /home/dezhu/Development/TanglePack/tangle-pack-dead-api-base out_base` and `bash probe/run.sh /home/dezhu/Development/TanglePack/tangle-pack-dead-api out_branch`. Both printed `ALLDONE`, wrote 63 lines each with no `error` rows, and produced no stderr beyond timing lines.

**Diff.** `cat out_base/*.jsonl > all_base.jsonl`, the same for the branch, then `diff all_base.jsonl all_branch.jsonl` gives **empty**. md5 is `bca69acaa8b8bc3be4438463ff174eef` for both.

**Coverage of the 63 items.** The 16 builds are k10, k28_one_blast, k28_two_blasts, p3, nested, nested_inner_blasts, inversion, inversion_one_blast, inversion_two_blasts, inversion_period2, heteroclinic, p3_16, p3_6_blasts, k28_8, nested_inner6 and p3_15. The series points are k28_series@0–8, p3_series@0–9, inversion_series@0–4, k32_series@0–10 and nested_outer_series@0–11. Over the 63 items (identical on both trees) there are 902 classes, 899 rules and 988 refined rules, and 60 items are reliable.

## Break-it checks (rubric 11–14)

Each break was made on the branch worktree, the named test run, and the file restored with `git checkout -- <file>`. `git status --short` was empty afterwards. The unbroken baseline was `test_refinement.py` 1 passed, and the three other files together 27 passed.

| # | Break | Guarded by | Result |
|---|---|---|---|
| M1 | Rubric 11's mutation, applied to the test file and then restored: `assert lo < float(new_point.cdist) < hi` → `float(p0.cdist)` | `test_refined_cdist_lies_strictly_between_its_neighbours` | 1 failed |
| M2 | `ManifoldMachine._refine_layer`: `mid_cd = 0.5 * (c0_arr + c1_arr)` → `mid_cd = c0_arr` (new point at its left neighbour's cdist) | same | 1 failed |
| M3 | `_refine_layer`: `flagged = np.flatnonzero(areas >= self.area_cutoff)` → `areas > np.inf` (refine nothing) | same | 1 failed |
| B1 | `IntersectionRegistry.remove_iterate`: `self.iterate_table.remove(source_id, n)` → `pass` | `test_remove_iterate_forgets_the_link_and_its_twin` (renamed) | 3 failed: the renamed test, `test_an_image_claimed_across_passes_links_neither_claimant`, `test_a_later_pass_does_not_relink_a_contested_image` |
| B2 | `IntersectionRegistry.register_iterate`: `self.iterate_table.register_iterate(...)` → `pass` | `test_an_image_claimed_across_passes_links_neither_claimant` (the new `iterate_table[z, 1] == y` at the old `:85`) | 1 failed |
| B3 | `TangleSession.classify_bridge`: final `return min(containing, ...)` → `return None` | `test_bridge_classification_uses_the_memoised_midpoint` (rewritten comprehension) | 1 failed |
| B4 | `Blast.py:381`: `result.stop_reason = "point_budget"` → `"max_iterations"` | `test_blast_budget.py` | 3 failed: `test_budget_stops_between_iterations_and_resumes`, `test_stored_factor_predicts_the_next_call`, `test_without_prediction_the_budget_stops_once_reached` |

No rewritten test needed strengthening, so **I made no edit and no commit on the branch**.

## Scripts (rubric 16, 17)

**Run.** On each tree I ran `cd <worktree> && MPLBACKEND=Agg PYTHONPATH=src timeout 900 env/bin/python scripts/<s>.py` for:

- `henon_basic_stable_and_unstable`
- `henon_initial_intersections`
- `henon_period_3_manifolds`
- `henon_new_intersection_structures`
- `henon_period_3_new_intersection_structures`
- `henon_bridge_classification`
- `tangle_workbench_test`
- `henon_blast_period_3`

Logs are in `te/scripts/{base,branch}/`.

**Normalisation:** worktree paths → `REPO`, `HH:MM:SS(,ms)`, `0x…` addresses, and the line number in `.py:N: UserWarning`.

The whole diff:

```
=== henon_new_intersection_structures.log
< REPO/src/tanglepack/numerics/graphviz.py:N: UserWarning: FigureCanvasAgg is non-interactive, and thus cannot be shown
<   plt.show()
< Graph: 8 nodes, 17 edges
=== henon_period_3_new_intersection_structures.log
< Graph: 207 nodes, 582 edges
```

The other five scripts are identical. All seven exit 0. `henon_blast_period_3.py` exits 1 on both trees with the identical `ValueError` at line 95 (rubric 17).

## Could not check / notes

- **Probe harness.** It is an adaptation, not a verbatim copy; the reason is given above. The B6b probe cannot run at `65edec4` (`WalkSearch.truncated` is gone), and it did not record words or matrices.
- **`ruff --select F401,F841,F811,F821` on the 26 edited Python files.** Base: "All checks passed". Branch: 56 F401 findings, all in `src/tanglepack/topology/__init__.py`. They are the names kept importable but dropped from `__all__`, which is Q8's stated design; the import block is B12's. There are no other findings. This is not a rubric item; I report it so the reviewer is not surprised by it.
- **Fixture left without a user.** `test_map_step_and_graph.py` lost the only users of `henon_p3_session` in that file. I did not check whether the fixture still has other users in `conftest.py`. This is not a rubric item, and fixtures were never named for deletion.
- Nothing failed.
