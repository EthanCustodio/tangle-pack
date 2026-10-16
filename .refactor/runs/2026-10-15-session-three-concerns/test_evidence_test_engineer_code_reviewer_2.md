# Test evidence: B9a `session-three-concerns` (test engineer, traversal 1, 2026-10-15)

## Summary

All 13 rubric items pass. I ran every check myself. Nothing computed or drawn differs between base and branch.

I added one commit on the branch:
- **Commit:** `5f0db5f` "Tests: the manual repartition event must leave a partition, not a cleared trellis".
- **What it changes:** one assertion in `tests/facade/test_session_caches.py::_repartition_by_hand` is made stronger (see R4).
- **Effect on the scope:** none. The test ids are unchanged, and the diff still touches the same 24 files.

## Setup

| Item | Value |
|---|---|
| Base | `517aa01`, worktree `/home/dezhu/Development/TanglePack/tangle-pack-session-three-concerns-base` |
| Branch | `refactor/session-three-concerns`, worktree `/home/dezhu/Development/TanglePack/tangle-pack-session-three-concerns` |
| Branch commits | C1 `d38dd3f`, C2 `86b453d`, C3 `7944f2e`, plus my test commit `5f0db5f` |
| Suite command | `MPLBACKEND=Agg PYTHONPATH=src ../tangle-pack/env/bin/python -m pytest -q -rxX -p no:cacheprovider` (`PY` below is `/home/dezhu/Development/TanglePack/tangle-pack/env/bin/python`) |
| Scratch directory | `$S = /tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b9/te/`. Logs: `suite_*.txt`, `ids_*.txt`, `breaks.py`, `breaks.log`, `canvas/`, `probe/`, `exp/` |

- I did not edit the main tree.
- I did not touch `tests/helpers/laws.py`, `tests/cases.py` or `tests/invariants/`.
- The branch worktree's `git status` is clean.

## R1: the suite passes on the branch

| Commit | Result | Time |
|---|---|---|
| `7944f2e` (C3) | **1208 passed, 47 skipped, 5 deselected, 1 xfailed** | 222 s |
| `5f0db5f` (after my commit) | **1208 passed, 47 skipped, 5 deselected, 1 xfailed** | 201 s |

- The only xfail is `test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing]` (C3).
- Grepping the logs for `XPASS|FAILED` finds nothing.

## R2: exactly the planned test ids are gone

Command, run in each worktree:

```
pytest --collect-only -q -p no:cacheprovider | grep '::' | sort
```

followed by `comm -23` / `comm -13` on the two sorted lists.

- **Base:** `1263/1268 tests collected (5 deselected)`. That is 1215 passed + 47 skipped + 1 xfailed.
- **Branch:** `1256/1261 tests collected (5 deselected)`. That is 1208 + 47 + 1.
- The ids are unchanged after my commit (`cmp` is equal).

**Only on base (10 ids):**
- `facade/test_session_bridge_classes.py::test_a_pip_trimmed_away_is_not_replaced_by_plotting[plot_strong_pip]` (deleted)
- `…[plot_strong_pip_candidates]` (deleted)
- `facade/test_session_fanouts.py::test_describe_fanouts_mention_every_fixed_point` (deleted)
- `facade/test_session_fanouts.py::test_no_argument_fanouts_return_one_entry_per_fixed_point` (renamed)
- `facade/test_session_fanouts.py::test_single_fixed_point_fanouts_return_that_trellis_list` (renamed)
- `plotting/test_plot_smoke.py::test_session_cartoon_delegate_honours_its_toggles` (renamed)
- `plotting/test_plot_smoke.py::test_session_plot_delegates_return_their_types` (deleted)
- `plotting/test_plot_smoke.py::test_session_plot_fanouts_draw_every_fixed_point` (deleted)
- `plotting/test_plot_smoke.py::test_session_plot_pseudoneighbors_computes_on_demand` (deleted)
- `unit/topology/test_partition_elements.py::test_session_closes_the_cross_trellis_gap` (deleted)

**Only on branch (3 ids):**
- `test_session_fanouts.py::test_classify_strong_pips_covers_every_fixed_point`
- `test_session_fanouts.py::test_classify_strong_pips_of_one_fixed_point_is_its_list`
- `test_plot_smoke.py::test_cartoon_honours_its_toggles`

That is exactly the plan's 7 deletions and 3 renames, with nothing else different.

## R13: each commit is green, full suite at every commit

I made detached scratch worktrees at `$S/c1` (`d38dd3f`) and `$S/c2` (`86b453d`) and ran the full suite at each, alongside the C3 run.

| Commit | Result | Time |
|---|---|---|
| C1 `d38dd3f` | 1215 passed, 47 skipped, 5 deselected, 1 xfailed (C3) | 222 s |
| C2 `86b453d` | 1210 passed, 47 skipped, 5 deselected, 1 xfailed (C3) | 222 s |
| C3 `7944f2e` | 1208 passed, 47 skipped, 5 deselected, 1 xfailed (C3) | 222 s |

- C2's count is 1215 − 3 smoke tests − 2 parametrized ids.
- C3's count is a further −2 (the describe test and the cross-trellis-gap test).
- There is no XPASS and no FAILED at any commit.
- I removed both scratch worktrees afterwards (`git worktree remove --force`, then `prune`).
- The commits are in the order C1 → C2 → C3, and each message names its concern.

## R3: the rewritten tests pass on base code

- **The base copy:** `git archive 517aa01` was extracted to `$S/basecopy`, a copy, not the base worktree.
- **What I copied in:** the branch's versions of `test_session_fanouts.py`, `test_session_caches.py`, `test_arrangement.py`, `test_nested_own_blast.py`, `test_session_bridge_classes.py`, `test_partition_elements.py` and `test_plot_smoke.py`.
- **The command:**

  ```
  pytest -q tests/facade/test_session_fanouts.py tests/facade/test_session_caches.py \
    tests/unit/topology/test_arrangement.py::test_stable_arcs_report_their_partition_elements \
    tests/regression/test_nested_own_blast.py tests/facade/test_session_bridge_classes.py \
    tests/unit/topology/test_partition_elements.py tests/plotting/test_plot_smoke.py::test_cartoon_honours_its_toggles
  ```

- **Result: 68 passed** on base code.
- **Base code really was imported:** `tanglepack.__file__` resolves to `$S/basecopy/src`, and `tanglepack.loom` there has no `cartoon_zones`.
- **Every rewritten test is in the 68:**
  - both `classify_strong_pips` tests;
  - all 49 of `test_session_caches`;
  - `test_stable_arcs_report_their_partition_elements`;
  - both tests of `test_nested_own_blast`;
  - all 7 remaining tests of `test_session_bridge_classes`, `_repartition_at`'s users included;
  - all 5 remaining tests of `test_partition_elements`;
  - `test_cartoon_honours_its_toggles`.
- **My strengthened `test_session_caches.py` (commit `5f0db5f`) also passes on the base copy:** 49 passed.
- **Exempt, as the plan allows:** the two `test_plot_topology` circle tests, which need `loom.cartoon`. R9 pins their figures instead.

## R4: the rewritten tests still have teeth

**Reading each rewrite against the base diff:**
- **`test_cartoon_honours_its_toggles`** keeps all three assertions: `isinstance(layout, plotting.CartoonLayout)`, `bridges_drawn == 0 and walks_drawn == 0`, and `not ax.patches`. Only the call changed.
- **`test_p3_circle_follows_the_zone_boundary`** keeps all of its asserts: `parents == {}`, the interior keys, the interior values, the arc-order permutation (two asserts), the cyclic order, and `_assert_interior_inside`.
- **`test_nested_inner_circle_sits_inside_the_outer_one`** keeps all of its asserts: `parents`, enclosures non-empty, every enclosure equal to the inner fixed point, the `alone` sub-selection non-empty, `layout.parents`, the two circle inequalities, and the node-inside loop. Only the call lines changed.
- **The `classify_strong_pips` halves** have verbatim bodies. Only the names and the first docstring changed.
- **`test_arrangement`'s element loop is not vacuous.**
  - Measured with a throwaway test placed in a base copy and a branch copy, then deleted: the loop checks **36 element entries** on the branch.
  - The `(branch, side, element_id)` list is identical in three runs: base with the old three stage calls, base with `session.topology(fp)`, and the branch.
  - Removing `session.topology(fp)` on the branch makes the test fail (`assert []`, the "partition must have produced a result" line).
- **`test_session_caches` had no teeth, and this is not new.**
  - Dropping the three stage calls in `_repartition_by_hand` still gave 49 passed. `clear_results()` alone satisfies `_own_partition(trellis) != partition`, because an empty tuple differs from the original partition.
  - Base has the same weakness: dropping `session.compute_pseudoneighbors(fp)` and the other two calls in the base copy also gave 49 passed.
  - **My fix (commit `5f0db5f`):** the test now asserts that the repartitioned trellis holds a partition (`repartitioned = _own_partition(trellis); assert repartitioned, …; assert repartitioned != partition`).
  - With the fix, the strengthened test passes on the branch (49) and on the base copy (49).
  - With the fix, dropping the three calls makes it fail: **7 failed** (every `manual_repartition` id).
  - The normal run has 2 partition results before and 2 after (`PARTS 2 2`, measured).

**Breaks planted on the branch** (`$S/breaks.py`, log at `$S/breaks.log`). Each break was applied alone, the 9 rewritten files or tests were run, and the break was restored with `git checkout -- src`. `git status` was empty after every restore.

| # | Break | Rewritten tests that fail |
|---|---|---|
| B1 | `plotting.cartoon_enclosures` returns `{}` | `test_nested_inner_circle_sits_inside_the_outer_one` |
| B2 | `loom.cartoon.cartoon_zones` returns `{}, {}, arc_order` | `test_nested_inner_circle…`, `test_p3_circle_follows_the_zone_boundary` |
| B3 | `Trellis.compute_pseudoneighbors` returns `[]` at once | 3 `test_session_bridge_classes` tests (`test_a_class_straddling_zones…`, `test_active_class_lies_in_the_zone…`, `test_products_need_only_the_numerics`); `test_session_caches::test_session_product_is_fresh_after_every_event`; `test_nested_inner_circle…`; both `test_nested_own_blast` tests |
| B4 | `Trellis.partition_stable_manifold` returns `[]` at once | 18 tests, among them `test_stable_arcs_report_their_partition_elements`, all 5 `test_partition_elements`, both circle tests, the toggles test, all 3 `test_session_caches` functions, both `test_nested_own_blast` tests and 4 `test_session_bridge_classes` tests |
| B5 | the cartoon's `if show_bridges:` becomes `if True:` | `test_cartoon_honours_its_toggles` |
| B6 | `session.classify_strong_pips` always returns the dict | `test_classify_strong_pips_of_one_fixed_point_is_its_list` |
| B7 | `cartoon_enclosures` keeps only the first of `fixed_points` | `test_nested_inner_circle…` (the new `fixed_points` argument is exercised) |

## R5: no plot, cartoon, fan-out or describe method is left on the session

The grep `def plot_\|def cartoon_\|_fanout\|def partition_element_for\|def describe_holes\|def describe_stable_partitions\|def _describe` on `src/tanglepack/loom/TangleSession.py` **prints nothing.**

## R6: the session is under its cap

`wc -l src/tanglepack/loom/TangleSession.py` gives **843**, under the 870 cap (base 1744).

## R7: no caller of a deleted method is left

Run literally, the grep is not empty. I classified every hit (31 lines in all), and the coder's reading is confirmed.

- **Workbench methods reached through `__getattr__`** (`TangleWorkbench.py`: `def plot_tangle` :718, `def plot_intersections` :1672, `def plot_all_bridges` :2044):
  - `scripts/henon_minimal_dual_graph.py:110`
  - `henon_strong_pip_nested.py:74-76`
  - `henon_pseudoneighbors_nested.py:94-96`
  - `henon_bridge_classes.py:248,251,265,266`
  - `henon_dual_graph.py:173-174`
  - `henon_k10_dual_graph.py:212-213`
  - `tangle_workbench_test.ipynb` json lines 163, 208, 253, 307, 308, 365, 366, 367, 426, 427, 428, 501, 502

  None of these is a deleted method.
- **Data, not code:** `stress/numerics/results/e07_breaking.jsonl:13,15` and `e07_breaking_quick.jsonl:7`, `session.partition_stable_manifold` inside the `traceback_tail` text of the committed result `480a2d1`.
  - The script that produced it no longer calls it: a grep of `stress/**/*.py` for every deleted name finds only `Pseudoneighbor.compute_pseudoneighbors` in `e13_both_intervals.py`, which is a module function.
- **Refined grep:** the same grep, without `plot_tangle|plot_intersections|plot_all_bridges` and without `.jsonl`, prints nothing.
- **Wider grep:** every deleted method name on any receiver in `scripts src tests notebooks`. The only hits are `Trellis` objects (`T`, `T1`, `T3`) and docstrings.
- **Two stale docstring cross-references** remain outside the scope, as the coder reported: `Trellis.py:1413` and `TopologyResults.py:427` both cite `TangleSession.partition_element_for`. These are for B9b or B12.

## R8: the scope holds

- **`git diff --stat 517aa01..`** at `7944f2e`: 24 files, +410/−1161. At `5f0db5f`: the same 24 files, +413/−1162.
- **The file list** is the plan's §3, exactly: 4 src (`cartoon.py` new), 8 tests, 11 scripts and the notebook.
- **Net lines per group, at `7944f2e`:**

  | Group | Added | Deleted | Net |
  |---|---|---|---|
  | src | +296 | −909 | −613 |
  | tests | +49 | −194 | −145 |
  | scripts | +65 | −58 | +7 |
  | **Total** | | | **−751** |

  The plan said about −800. My commit moves tests by +3/−1.
- **The notebook diff** is one source line: `session.plot_resonance_zones()` becomes `tanglepack.loom.plot_resonance_zones(session.resonance_zones.values())`. The notebook imports `tanglepack`.
- **`plotting.py`** goes from 3404 to 3479 lines, **+75**, under the +85 cap.
- **No topology module imports loom:** `grep -n loom src/tanglepack/topology/*.py | grep -i import` is empty. The only `loom` hits in `plotting.py` are three docstring citations, at :2519, :2895 and :3017.
- **Imports:**
  - `from tanglepack.loom import cartoon_zones, plot_resonance_zones` gives the same objects as in `tanglepack.loom.cartoon`.
  - `plotting.cartoon_enclosures` exists.
  - `topology.__all__` still has 48 names, without `cartoon_enclosures`.
  - `cartoon.py` is 209 lines.

## R9: byte-identical canvases

My script is `$S/canvas/canvas_te.py`. It extends the coder's `scratchpad/b9a/canvas.py`.
- **Base** draws through the session methods; **the branch** draws through the moved functions or `Trellis.plot_*`.
- For each figure it compares the sha256 of `fig.canvas.buffer_rgba()` and the sha256 of an artist dump (path vertices, offsets, texts).
- It runs with `REPO=<worktree> MPLBACKEND=Agg PY canvas_te.py {a,b,c}`.
- All stderr is empty.

**32 of 32 pairs are equal, canvas and artists both.**
- **(a) The plan's list:**
  - the circle cartoon on p3, nested and k10;
  - the line cartoon on k10 with every toggle off;
  - `plot_resonance_zones` on p3 and on nested;
  - `plot_dual_graph` on k10;
  - per fixed point of nested, after `topology(fp)`: `plot_strong_pip`, `plot_strong_pip_candidates`, `plot_holes`, `plot_pseudoneighbors` and `plot_stable_partition`.

  That is 17 figures.
- **(b) Beyond the plan:** the circle cartoon, the zones and the toggles-off line cartoon on `heteroclinic`, `inversion_period2` (the case Q2 (b) would have changed), `inversion` and `nested_inner_blasts`. That is 12 figures.
- **(c) The natural fan-out path, with no `topology(fp)` first:**
  - base: `session.plot_strong_pip_candidates()`, `plot_strong_pip()` and `plot_pseudoneighbors()` on nested;
  - branch: the per-fixed-point `session.trellis(fp).plot_*` loop.

  That is 3 figures.

## R10: the 42 symbolic-itinerary PNGs are identical

- **Exports:** `git archive 517aa01` and `git archive 7944f2e`, extracted under `$S/exp/{base,branch}/<run>/tree`. One fresh export per script run, outside the repo.
- **Command:** `MPLBACKEND=Agg PYTHONPATH=src $PY …`.
- **The PNGs:** written after a marker file, then sha256-compared. I never compared them with the committed PNGs.

| Run | PNGs | Exit codes | stdout / stderr | Serial time, base / branch |
|---|---|---|---|---|
| `henon_symbolic_itineraries.py` | 18 | 0 / 0 | identical | 6.3 / 6.2 s |
| `_period3.py` | 6 | 0 / 0 | identical | 2.1 / 2.0 s |
| `_nested.py` | 6 | 0 / 0 | identical | 3.7 / 3.6 s |
| `_nested.py --heteroclinic` | 6 | 0 / 0 | identical | **14.5 / 14.5 s** |
| `_inversion.py` | 6 | 0 / 0 | identical | 1.8 / 1.8 s |

- **All 42 sha256 are equal**, in the first run and again in a serial rerun. The digest of the sorted list is `b590e7d7…` on both sides.
- **The heteroclinic run takes 14.46 s on both sides,** 0.46 s over the plan's "≤ 14 s" figure. It is identical on base, so it is not a regression; I am stating the measurement as it is.
- **Normalisation:** stdout and stderr are compared after replacing `0x…` addresses and the export path.

## R11: the 11 rewritten scripts behave the same

These were run in the same per-script exports.
- A `sitecustomize.py` placed on `PYTHONPATH` (`$S/hook`) replaces `plt.show` with a dump of the sha256 of every open figure's canvas. This way the show-only scripts are compared too.
- `symbolic_figures.py` is a module with no `__main__`. It is run through the five R10 scripts.

| Script | Exit codes | stdout (normalised) | PNGs written | Canvases shown |
|---|---|---|---|---|
| `henon_pseudoneighbors` | 0 / 0 | identical (2 lines) | — | 1, equal |
| `henon_pseudoneighbors_nested` | 0 / 0 | identical (5) | — | 2, equal |
| `henon_strong_pip_nested` | 0 / 0 | identical (4) | — | 1, equal |
| `henon_minimal_dual_graph` | 0 / 0 | identical (72) | 1, equal | — |
| `henon_dual_graph` | 0 / 0 | identical (9) | 2, equal | — |
| `henon_k10_dual_graph` | 0 / 0 | identical (6) | 1, equal | — |
| `henon_bridge_classes` | 0 / 0 | identical (21) | 1, equal | — |
| `henon_resonance_zone` | 0 / 0 | identical (4) | — | 1, equal |
| `henon_bridge_classification` | 0 / 0 | identical (13) | — | 1, equal |
| `henon_blast_period_3` | **1 / 1** | identical (0) | — | — |
| `symbolic_figures` (via R10) | 0 / 0 | identical | 42, equal | — |

- **stdout normalisation:** raw stdout differs only in `henon_strong_pip_nested`, by `FixedPoint object at 0x…` addresses. It is identical after `0x…` → `0xADDR`.
- **`henon_blast_period_3.py` is still dead, identically on both sides.**
  - Traceback at `scripts/henon_blast_period_3.py", line 95, in <module>`.
  - Last line: `ValueError: Intersection 12 is not a strong-pip candidate. Candidates: [31, 32, 33]. Run classify_strong_pips() first, or choose one of these.`
- **stderr** is identical in every script except `henon_k10_dual_graph.py`, which has 32 lines on base and 31 on the branch. This is the change the coder declared. The exact diff:

  ```
  < INFO:__main__:<tanglepack.numerics.FixedPoint.FixedPoint object at 0xADDR>:
  < left partition (p1, orbit 0):
  ---
  > INFO:__main__:left partition (p1, orbit 0):
  ```

  The `_describe` fan-out's fixed-point header line is gone. The report text is otherwise the same.
- **The logger-name change** (`tanglepack.loom.TangleSession` → `tanglepack.loom.cartoon` for "no cartoon side" and "No resonance zones to plot") appears in **no** script's stderr in these runs. It cannot be observed here, and no test asserts on it.

## R12: no number changes on the corpus

- **The harness:** copied to `$S/probe/` (`probe.py` unchanged, `run.sh` with only its directory re-pointed). The b10 copies are untouched.
- **The command:** `./run.sh <base worktree> out_base` and `./run.sh <branch worktree> out_branch`. Both runs ended `ALLDONE`.
- **Each tree imports its own code:** with the probe's `sys.path` insert, `tanglepack.__file__` resolves to each worktree's own `src`, and `hasattr(TangleSession, 'plot_holes')` is True on base and False on the branch.
- **The result:**
  - 63 lines each, 0 `error` rows;
  - `cmp all_base.jsonl all_branch.jsonl` is equal: **0 differences**;
  - both files have sha256 `e2f4d14d…3902`;
  - they are byte-equal to B10's `all_branch.jsonl` too.

## Not checked, or worth knowing

- **`tanglepack.loom.plot_resonance_zones` has no test**, on base or on branch: a grep of `tests/` finds none in either tree. It is pinned here only by R9 (zones on 6 cases) and by the R11 show-dumps.
  - A smoke test would change the collected ids, which rubric 2 forbids. So I did not add one. It is a candidate for a later slice.
- **The heteroclinic itinerary run** takes 14.46 s on both sides (see R10).
- **The scratch exports** were deleted after the comparison; the hashes, logs and stdout/stderr files remain in `$S/exp/`.
