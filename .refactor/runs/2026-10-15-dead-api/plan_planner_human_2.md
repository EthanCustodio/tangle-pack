# Plan — B10 `dead-api` (planner, 2026-10-15, base `65edec4`)

## In short

- **The answer.** At `HEAD`, every name in the plan's B10 table has no user on any computed path except `release_manifold`. That one is called three times from `TangleWorkbench` (`:1360`, `:1761`, `:1861`), and `test_bridge_identity.py:257` covers it.
- **The evidence.** The analyst deleted everything else in a scratch copy of `HEAD` and ran the full suite: 1223 → 1215 passed. The 8 missing tests are exactly the ones that existed only for the deleted names (5 registry-graph tests, plus `from_branch` × 3). The only xfail is still C3, with no XPASS. The 7 scripts that call a deleted name and still run print the same output after a one-to-three-line rewrite, apart from the removed `Graph: N nodes` line.
- **Recommendation: change it.** Delete, one commit per row.
- **What goes:**
  - `grow_x_times`
  - `_get_iterate` and `_get_refined_point`
  - `trim_stable_at_intersection`
  - the registry's networkx graph, with `build_intersection_graph` / `visualize_intersection_graph` and the whole of `graphviz.py`
  - five registry query methods, plus the two `IterateTable` array methods they orphan
  - `classify_bridges`
  - `BlastResult.terminated_early`
  - then trim `topology.__all__` from 104 to 48 names
- **What stays:**
  - `release_manifold`
  - every law-tier oracle: `euler_characteristic`, `Arrangement.preimage_of`, `StableNode.other_face`, `preimage_bridges`, and `DualGraph.graph()`, which `laws.py:2716` reads
  - `refine_manifold`
  - the registry's `filter`, `from_fixed_point` and `iterate_orbit`
  - every import line in `topology/__init__.py`
- **Size.** Measured net over src + tests: **−950 lines** (src −822, tests −128), plus −56 from the `__all__` trim. Scripts change by a further −21.
- **Behaviour.** Nothing on a computed path changes. Visible effects:
  - `import tanglepack` no longer loads networkx.
  - Two demo scripts lose their graph sections.
- **What you decide.** Q1–Q8 at the end. Each has a default, so "approved" accepts all eight. Q3, the dead script `henon_blast_period_3.py`, is the one most worth reading.

## Scope

- **Slug:** `dead-api`. Branch `refactor/dead-api` in the worktree `../tangle-pack-dead-api`, made from `65edec4`. The base probe is a detached worktree `../tangle-pack-dead-api-base` at `65edec4`.
- **What the branch never touches:** `CLAUDE.md`, `docs/`, `.refactor/`, `artifacts/`, notebooks, `figures/`, and the two untracked sweep scripts.
- **Suite command** (run inside a worktree): `MPLBACKEND=Agg PYTHONPATH=src ../tangle-pack/env/bin/python -m pytest -q -rxX -p no:cacheprovider`
- **Size:** 27 paths, about 24 of them edits. That is more than the five-file guideline, but every edit is a deletion or a one-to-three-line caller rewrite and nothing changes logic. The request already flagged this as larger, and splitting it would just mean eight runs of the same thing.

**File list** (nothing outside it is edited; `[Qn]` = only if that default stands):

| # | Path | Commits |
|---|---|---|
| 1 | `src/tanglepack/numerics/ManifoldMachine.py` | C1, C2 |
| 2 | `src/tanglepack/numerics/IntersectionRegistry.py` | C4, C5 |
| 3 | `src/tanglepack/numerics/IterateTable.py` | C5 [Q7] |
| 4 | `src/tanglepack/numerics/TangleWorkbench.py` | C4 |
| 5 | `src/tanglepack/numerics/graphviz.py` (deleted) | C4 |
| 6 | `src/tanglepack/numerics/__init__.py` | C4 |
| 7 | `src/tanglepack/loom/ResonanceZone.py` | C3 |
| 8 | `src/tanglepack/loom/__init__.py` | C3 |
| 9 | `src/tanglepack/__init__.py` | C3 |
| 10 | `src/tanglepack/loom/TangleSession.py` | C6 [Q2] |
| 11 | `src/tanglepack/loom/Blast.py` | C7 |
| 12 | `src/tanglepack/topology/__init__.py` | C8 [Q8] |
| 13 | `tests/unit/numerics/test_refinement.py` | C2 |
| 14 | `tests/unit/numerics/test_map_step_and_graph.py` | C4 |
| 15 | `tests/unit/numerics/test_iterate_inference.py` | C4 |
| 16 | `tests/unit/numerics/test_generation_and_caches.py` | C4, C5 |
| 17 | `tests/unit/numerics/test_intersection_registry_fixes.py` | C5 [Q4] |
| 18 | `tests/unit/loom/test_resonance_zone.py` | C6 |
| 19 | `tests/unit/loom/test_blast_budget.py` | C7 |
| 20 | `scripts/henon_basic_stable_and_unstable.py` | C1 |
| 21 | `scripts/henon_initial_intersections.py` | C1 |
| 22 | `scripts/henon_period_3_manifolds.py` | C1 |
| 23 | `scripts/henon_new_intersection_structures.py` | C4, C5 [Q1] |
| 24 | `scripts/henon_period_3_new_intersection_structures.py` | C4, C5 [Q1] |
| 25 | `scripts/henon_bridge_classification.py` | C6 [Q2] |
| 26 | `scripts/henon_blast_period_3.py` | C7, line `:123` only [Q3] |
| 27 | `stress/numerics/e02_gpu_end_to_end.py` | C2, line `:172` only [Q5] |

## Commits, in order

Every commit leaves the suite green. Line numbers are at `65edec4`, from the findings. Each commit message names its row: `Dead API: <name> (B10)`.

### C1 — `grow_x_times`
- **Delete:** `ManifoldMachine.py:490-508`.
- **Docstring:** the class docstring `:128` points at `new_grow_manifold` instead.
- **Rewrite the three scripts:** `henon_period_3_manifolds.py:62,63,69,70`, `henon_initial_intersections.py:50,54` and `henon_basic_stable_and_unstable.py:51,54`. Each call becomes:
  ```python
  for _ in range(num_iterations):
      man_machine.new_grow_manifold(fp, "unstable")
  ```
  These scripts drive a bare `ManifoldMachine` with no workbench, so the plan's `grow_n_times` does not apply.
- **Notebooks:** left alone (Q6).

### C2 — `_get_iterate` + `_get_refined_point`
- **Delete:** `ManifoldMachine.py:815-851`. `_get_preiterate` and `_cache_preiterate` stay, because `_refine_layer` uses them.
- **Rewrite the test** `test_refinement.py:33-49`:
  ```python
  machine.area_cutoff = 0.0
  (p0_, new), (_, p1_) = machine._refine_layer([(p0, p1)], manifold, viewer, None, set())
  assert lo < new.cdist < hi          # unchanged strength
  ```
  The exact call is in `<scratchpad>/proto_refine.py`.
- **Stress file:** drop `"_get_iterate"` from the profiler category set in `e02_gpu_end_to_end.py:172` (Q5).

### C3 — `trim_stable_at_intersection`
- **Delete:** `ResonanceZone.py:473-495`.
- **Exports:** remove `loom/__init__.py:26,39` and `tanglepack/__init__.py:62,108`.

### C4 — the registry graph and `graphviz.py`
- **Delete the file** `numerics/graphviz.py` (358 lines), with its import and `__all__` entry in `numerics/__init__.py:35,67`.
- **`TangleWorkbench.py`:**
  - delete the graph-view section `:2047-2094`
  - delete `import networkx` (`:41`) and `from . import graphviz` (`:43`)
  - remove the module docstring bullet at `:15`
  - reword the `:173` docstring sentence so it names only `infer_iterates`
- **`IntersectionRegistry.py`:**
  - delete `graph` / `_rebuild_adjacency_edges` / `_consecutive_on_branch` (`:530-656`)
  - delete the `_graph*` fields (`:178-182`)
  - delete the node adds in `add` (`:236-246`) and in `reindex_from` (`:350`, `:363-373`)
  - delete the edges in `register_iterate` (`:674-682`) and `remove_iterate` (`:693-696`)
  - delete `import networkx` (`:23`)
  - reword the docstrings at `:91`, `:120-122`, `:184-189`, `:201`, about `:239`, `reindex_from`, `register_iterate` and `remove_iterate`, so no sentence mentions a graph, edge or node of the registry
- **Tests:**
  - `test_map_step_and_graph.py`: delete the five graph tests (`:201-298`) and the docstring paragraph at `:8-11`.
  - `test_iterate_inference.py`:
    - `:85` becomes `registry.iterate_table[z, 1] == y`
    - `:147` becomes `registry.iterate_table[x, 1] == y`
    - drop `:95` and `:153` (the table asserts at `:90` and `:152` already say the same)
    - rename `test_remove_iterate_forgets_the_link_its_twin_and_its_edge` to `…_the_link_and_its_twin`
  - `test_generation_and_caches.py:116`: drop the `registry.graph()` read.
- **Scripts** `henon_new_intersection_structures.py:76-95` and `henon_period_3_new_intersection_structures.py:102-…`: drop the "Live graph" and "Intersection graph" sections (Q1).

### C5 — registry query methods (+ `IterateTable` arrays, Q7)
- **Delete in `IntersectionRegistry.py`:**
  - `on_cdist_range` (`:763-784`; and the `:95` and `:264` mentions)
  - `from_branch` (`:816-849`; and `:97`)
  - `nearest_by_unstable_cdist` (`:851-873`)
  - `as_forward_array` / `as_backward_array` (`:875-886`)
  - the `Iterable` / `NDArray` imports this leaves unused
- **Delete in `IterateTable.py`:** `as_forward_array` / `as_backward_array` (`:98-125`), plus its `numpy` / `NDArray` imports if they become unused.
- **Tests:**
  - `test_generation_and_caches.py:114,115` (the read-only probe) and `:142`: drop the reads of the deleted methods.
  - `test_intersection_registry_fixes.py:168-186`: delete the `from_branch` test (Q4).
- **Scripts:**
  - `on_cdist_range(a, b, stability="stable")` becomes `registry.filter(lambda ix: a <= ix.stable_cdist <= b)`.
  - `from_branch(0, stability="stable")` becomes `registry.filter(lambda ix: ix.is_anchor or (ix.manifold_b_key is not None and ix.manifold_b_key[3] == 0))`.
  - Drop the array calls (`:98-99`, `:128-129`).

### C6 — `classify_bridges` (Q2)
- **Delete:** `TangleSession.py:1618-1634`.
- **Rewrite** `henon_bridge_classification.py:83` and `test_resonance_zone.py:115` to `{b: session.classify_bridge(b) for b in session.workbench.bridges}`.

### C7 — `terminated_early`
- **Delete:** `Blast.py:201-204`.
- **Rewrite `test_blast_budget.py`:**
  - drop `:54` (it follows from `:53`'s `stop_reason == "point_budget"`)
  - drop `:95` (a tautology)
  - reword the docstrings at `:6-7` and `:90-91` to speak of `stop_reason`
- **Dead script:** in `henon_blast_period_3.py:123`, print `result.stop_reason == "frontier_empty"` instead (Q3). The script stays dead at `:95` and is not otherwise touched.

### C8 — `topology.__all__` to 48 names (Q8)
- **Keep these 48:**
  - **(a) names imported from outside the subpackage (15):**
    - the ten re-exported by `tanglepack/__init__.py:45`
    - `side_after`, `letter`, `bridge_row_diagnosis`, `check_holes_share_bridge_side`, `plot_stable_partition`
  - **(b) stage objects (15):** `MinimalTrellis`, `minimal_trellis`, `HomotopyPartition`, `IteratedHomotopyPartition`, `BridgeClass`, `BridgeClassEntry`, `BridgeClassTable`, `DualGraph`, `Walk`, `ClassDynamics`, `SymbolicDynamics`, `inert_letters`, `symbolic_dynamics`, `bridge_classes`, `zone_label`
  - **plotting (18):**
    - drawing functions and layouts: `plot_dual_graph`, `plot_minimal_trellis`, `plot_walk`, `plot_bridges_by_class`, `plot_transition_graph`, `plot_itinerary_table`, `plot_dual_graph_cartoon`, `CartoonLayout`, `ZoneArcOrder`
    - helpers: `dual_graph_legend_handles`, `face_point`, `minimal_trellis_legend_handles`, `stable_node_point`, `walk_legend_handles`, `dual_graph_cartoon_layout`, `dual_graph_cartoon_legend_handles`, `ZoneLayout`, `dual_graph_zone_layout`
- **Leave alone:** every `from .X import …` line. The 56 names that leave `__all__` stay importable by name. The import block is B12's.
- **Cost:** −56 lines.

## Behaviour

**None** on any computed path. Words, matrices, class counts and `is_reliable` are identical on all 63 probe items. Outputs that do change:
- `import tanglepack` no longer imports `networkx`. It stays a dependency in `pyproject.toml`, and topology imports it lazily.
- `henon_new_intersection_structures.py` and `henon_period_3_new_intersection_structures.py` no longer print `Graph: N nodes, M edges` or open the graph figure.
- `henon_blast_period_3.py:123` would print a bool from `stop_reason` instead of `terminated_early`, but the script never gets there (it still raises at `:95`).

## Not doing

- `release_manifold`: it is live.
- The oracles listed above.
- `refine_manifold`.
- `partition_element_for` and the session fan-outs and plot delegates: these are B9's.
- The plotting split (SIMP-11).
- The `topology/__init__` import lines, and the `CLAUDE.md` rows (`:136` graphviz, `:154` `terminated_early`): these are B12's.
- Notebooks (plan Deferred list, Q6).
- Reviving `henon_blast_period_3.py` (Q3).
- The analyst's out-of-scope notes:
  - the period-3 demo's anchor label
  - the `[2, 8]` / `[2, 6]` print mismatch
  - `IterateTable.forward_depth` / `forward_chain` reading direct entries only

## Rubric

The reviewer checks each item by running or reading, never from a report. The figures assume all defaults; the variants follow in brackets.

1. **The suite passes on the branch.** `-rxX` shows exactly **1215 passed, 47 skipped, 5 deselected, 1 xfailed (C3)** and no XPASS. Base is 1223/47/5/1. [Q4 = keep: 1218 passed.]
2. **Exactly the intended tests are gone.** Diff the set of collected node ids (`--collect-only -q`) between base and branch:
   - The only ids that disappear are the five graph tests in `test_map_step_and_graph.py` and the three `from_branch` parametrisations in `test_intersection_registry_fixes.py`.
   - The only id that changes is the rename in `test_iterate_inference.py`.
   - No id appears that base lacks, apart from that rename.
3. **No law-tier pin moved, on either side.**
   - `git diff --stat 65edec4..refactor/dead-api -- tests/helpers/laws.py tests/cases.py tests/invariants/` is empty.
   - `pytest tests/invariants` passes on both worktrees with the same pass/skip/xfail counts.
4. **The diff stays inside the file list.** `git diff --stat 65edec4..refactor/dead-api` lists only paths from the 27-row table, minus any a declined Q-default removes. Nothing under `CLAUDE.md`, `docs/`, `.refactor/`, `artifacts/`, `notebooks/` or `figures/` is touched, and `git status` in the main tree is unchanged by the run.
5. **Eight commits, one per row (C1–C8)**, and each commit alone passes the suite. At minimum, the reviewer runs the suite at the C4 and C5 commits.
6. **Net lines over `src tests`** (from `git diff --shortstat 65edec4 -- src tests`):
   - at the C7 commit: **≤ −950** (measured −950 in the prototype; the reworded docstrings are net 0, and the dropped unused imports only lower it)
   - for the whole branch: **≤ −1000** (C8 adds −56)
   - [Q4 = keep: ≤ −895 / ≤ −950; Q7 = keep: +29 to each]
7. **No deleted name survives in Python files.** This returns nothing:
   ```
   grep -rnw --include='*.py' -e grow_x_times -e _get_iterate -e _get_refined_point -e trim_stable_at_intersection -e build_intersection_graph -e visualize_intersection_graph -e _rebuild_adjacency_edges -e _consecutive_on_branch -e graphviz -e on_cdist_range -e from_branch -e nearest_by_unstable_cdist -e as_forward_array -e as_backward_array -e classify_bridges -e terminated_early src tests scripts stress
   ```
   - The grep covers `*.py` files only. That deliberately excludes the notebooks (Q6) and the untouched `docs/`.
   - Exemptions: none under the defaults. [Q3 = (b): `scripts/henon_blast_period_3.py:123` is exempt. Q5 = keep: `stress/numerics/e02_gpu_end_to_end.py:172` is exempt. Q4 = keep: `from_branch` is exempt.]
8. **No graph left in numerics.** `grep -rn "networkx\|_graph" src/tanglepack/numerics` returns nothing, and `src/tanglepack/numerics/graphviz.py` does not exist.
9. **The imports still work from outside the repo.** From `/tmp`, with `PYTHONPATH=<worktree>/src`, each of these succeeds in a fresh interpreter:
   - `import tanglepack`
   - `from tanglepack.numerics import *`
   - `from tanglepack.topology import *`
   - `from tanglepack.loom import *`
   - `from tanglepack.topology import StableNode, free_reduce, MembersDisagree, plotting` (names dropped from `__all__` stay importable)

   Also, `python -c "import sys, tanglepack; print('networkx' in sys.modules)"` prints `False` on the branch and `True` on base.
10. **`topology.__all__` is exactly the 48 names listed under C8.**
    - `len(tanglepack.topology.__all__) == 48`.
    - All 15 (a) names are in it.
    - Every name in it resolves: `all(hasattr(tanglepack.topology, n) for n in __all__)`.
11. **The refinement test keeps its teeth.** `test_refinement.py` calls `_refine_layer` with `area_cutoff = 0.0` and still asserts `lo < cdist < hi` on the new point. Read it: it must not loosen to `<=` or use a tolerance. On the branch, swapping the new point for `p0` must make it fail.
12. **The classify test keeps its teeth.** `test_resonance_zone.py` still asserts, for every in-zone bridge, that its midpoint lies in the zone and that `classify_bridge` agrees. Only the comprehension changes; the diff on that file is one line.
13. **The iterate links are still asserted.** `test_iterate_inference.py` still asserts the forward link through `iterate_table[·, 1]` at the old `:85` and `:147`, and its twin assertions stay.
14. **Blast tests keep a stop reason.** `test_blast_budget.py` still asserts `stop_reason == "point_budget"` at the budget stop, and still asserts a `stop_reason` in the no-budget test.
15. **Base and branch agree on the 63 probe items.**
    - Copy `<scratchpad>/b6b/{run.sh,probe.py,summ.py}` to `<scratchpad>/b10/`.
    - Make `REPO` and the `out/` directory arguments: `REPO` in `probe.py:4` and the `cd` and `D` in `run.sh`.
    - Run it once against `../tangle-pack-dead-api-base` and once against `../tangle-pack-dead-api`.
    - Diffing the two outputs gives **0 differences on all 63 items** (16 builds + 5 series): words, matrices, class counts, `is_reliable`.
16. **Rewritten scripts print the same.** With `MPLBACKEND=Agg`, run each of these on the base worktree and on the branch worktree, capturing stdout and stderr:
    - the six rewritten scripts (rows 20–25)
    - `scripts/tangle_workbench_test.py`

    Then normalise paths, `HH:MM:SS` timestamps, `0x…` addresses and the line numbers in `UserWarning` lines. The only remaining differences allowed:
    - the `Graph: … nodes, … edges` line in the two `*_new_intersection_structures.py` scripts
    - the graphviz `plt.show` warning in `henon_new_intersection_structures.py`

    Any other differing line fails.
17. **The dead script is still dead in the same way.** `scripts/henon_blast_period_3.py` raises the same `ValueError: Intersection 12 is not a strong-pip candidate` at `:95` on base and on branch, and its diff is one line [Q3 = (a)]. The run report says it is dead.
18. **The long CLAUDE.md scripts compile but are not run.** These are `henon_symbolic_itineraries*.py`, the two `henon_period3_*_sweep.py`, `show_figure.py` and `stress/numerics/e12_reverify.py`. They take minutes and write into `figures/`, and none imports a deleted name (rubric 7 covers that), so they are **not run**. Their computed path is covered by rubric 15. The check is that `python -m py_compile` passes on each on the branch.

## Questions

Each has a default. "Approved" accepts all eight.

- **Q1.** The two demo scripts `henon_new_intersection_structures.py` and `henon_period_3_new_intersection_structures.py` exist to show off the deleted API. Rewrite them or delete them?
  - **Default: rewrite.** Use `filter` for the two queries and drop the graph and array sections. Verified to run.
- **Q2.** `classify_bridges` has a live script caller (`henon_bridge_classification.py`). Delete the method and rewrite both the script and its test as a comprehension over `classify_bridge`, or keep the method?
  - **Default: delete and rewrite.**
- **Q3.** `scripts/henon_blast_period_3.py` has been dead since B3 (pip 12 is no longer a candidate, so it raises at `:95`), and its `:123` reads `terminated_early`. The standing rule is that a dead script is reported, never silently fixed or deleted.
  - **Default (a):** change only `:123` to `stop_reason == "frontier_empty"` so no deleted name survives. The script stays dead, and the run report says so.
  - **(b):** leave the file untouched and exempt it from rubric 7.
  - Reviving or deleting the script is not on offer in this run.
- **Q4.** `from_branch`'s three-way test is the only test of the D4 anchor rule written as a registry query. The rule itself lives on in `Trellis` and `Arrangement`; the analyst did not trace their tests.
  - **Default: delete it with the method.**
  - Saying "keep" keeps both `from_branch` and its test, and adjusts rubrics 1, 6 and 7 as bracketed.
- **Q5.** `e02_gpu_end_to_end.py:172` names `"_get_iterate"` as a profiler category string. It is not a call.
  - **Default: drop the string.**
- **Q6.** Two notebooks call `grow_x_times`, and `notebooks/archive/period_three.ipynb` is already broken (old signature).
  - **Default: leave them** (plan Deferred list). Rubric 7 greps `*.py` only.
- **Q7.** `IterateTable.as_forward_array` / `as_backward_array` are not in the plan's table, and row 6 orphans them.
  - **Default: delete them in C5.**
- **Q8.** `topology.__all__`: trim only the list, or the import block too?
  - **Default: trim `__all__` to the 48 names in C8 and keep every import line,** so no name stops being importable. The import block is B12's.
