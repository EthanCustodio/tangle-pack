# Findings — B10 `dead-api` (analyst, 2026-10-15, `HEAD` = `65edec4`)

**The working tree is dirty.** It carries uncommitted edits to:
- `.refactor/decisions.md`
- `.refactor/knowledge/*.md` (this run adds to four of them)
- `CLAUDE.md`
- `docs/alignment_plan_2026-10-09.md` and `docs/code_fidelity_audit_2026-10-09.md`

It also has an `artifacts/request_origin_analyst.md` rename, untracked run dirs, figures and two sweep scripts. The slice branches from `HEAD` and must touch none of these.

## What kind of request this is

This is a deletion refactor of unused API: one commit per row of the plan's B10 table (`docs/alignment_plan_2026-10-09.md:810-845`). I read it as: check every row again at `HEAD`, give each a verdict, measure what deleting it removes, and sort the 104 `topology.__all__` names by who uses them. I have not re-planned the slice.

## Short answer

- **Every row the table calls dead is dead, except `release_manifold`.** It is called from three places on the computed path and a regression test covers it. Keep it.
- **Deleting all the other rows leaves the suite passing.** I deleted them in a scratch copy of `HEAD`:
  - At `HEAD`: 1223 passed, 47 skipped, 5 deselected, 1 xfailed (C3).
  - In the copy: 1215 / 47 / 5 / 1 xfailed (C3).
  - The 8 missing tests are exactly the deleted ones: 5 registry-graph tests, plus `from_branch` × 3 parametrisations.
- **Every script that called a deleted name runs after a one-to-three-line rewrite.** Output is identical except the removed `Graph: N nodes` line.
- **Net is −950 lines (src + tests), measured, before any `__all__` trim.** That is src −822 and tests −128; scripts are −21 on top. So −900 holds.
- **Two items fall outside the plan's table:**
  - `IterateTable.as_forward_array` / `as_backward_array` lose their only callers.
  - `stress/numerics/e02_gpu_end_to_end.py:172` names `_get_iterate` as a string.
- **Q1–Q3 below need the author;** Q4–Q8 have defaults the planner can propose.

**How I checked (repo untouched):**
- Extracted `git archive HEAD` into `<scratchpad>/proto`.
- Applied the deletions with `<scratchpad>/delete.py`, `delete_tests.py` and `delete_scripts.py`.
- Ran the full suite there (`suite_proto.log`) and in the repo (`suite_base.log`).
- Ran the 8 affected scripts on `HEAD` (`run_*.log`) and on the copy (`p_*.log`).
- Ran fresh-interpreter import checks from `/tmp`.

`<scratchpad>` is `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad`.

## Row by row

### 1. `grow_x_times` — delete, rewrite 3 scripts

- **What it is.** `src/tanglepack/numerics/ManifoldMachine.py:490-508`, 20 lines. It just loops `self.new_grow_manifold(...)` `num_times`. The class docstring mentions it at `:128`; point that line at `new_grow_manifold` instead.
- **Callers.** 8 lines in 3 scripts:
  - `scripts/henon_period_3_manifolds.py:62,63,69,70`
  - `scripts/henon_initial_intersections.py:50,54`
  - `scripts/henon_basic_stable_and_unstable.py:51,54`
- **The plan's rewrite to `grow_n_times` does not fit.** All three drive a bare `ManifoldMachine` with no workbench.
  - The minimal rewrite is `for _ in range(num_iterations): man_machine.new_grow_manifold(fp, "<stab>")`.
  - Verified: the rewritten copies run and print exactly what `HEAD` prints (`rw_*.py`, `rw_*.log`).
- **Notebooks also call it.**
  - `scripts/manifold_cutting_test.ipynb` uses the current signature.
  - `notebooks/archive/period_three.ipynb` uses an old signature `(stable_manifold, fixed_point, n)`, so it is already broken.
  - The plan defers notebooks. Rubric 2's `grep -rn … scripts` would match the `.ipynb` unless limited to `*.py` (Q6).
- **Size:** src −20 (+1/−1 docstring); scripts +16/−8.

### 2. `release_manifold` — KEEP

- **Location.** `src/tanglepack/numerics/Tangle.py:648` (the plan said `:711`).
- **Callers, all live:**
  - `TangleWorkbench.py:1360`: `compute_intersections` dropping an image curve that adds no segment.
  - `TangleWorkbench.py:1761`: `_register_bridge`, when a duplicate bridge object arrives.
  - `TangleWorkbench.py:1861`: `clear_bridges`.
- **Probe (`probe_release.py`, on `build_k28()` through `topology()`):**
  - On a full build it is called 5 times from `clear_bridges`, all no-ops because `clear_all` ran first.
  - On the direct `clear_bridges` path it does real work. `tests/unit/numerics/test_bridge_identity.py:257` covers that path.
- **Making it private does not help.** The call crosses objects (workbench → Tangle), so a leading underscore would just be private access from outside. No change.

### 3. `trim_stable_at_intersection` — delete with its exports

- **What it is.** `src/tanglepack/loom/ResonanceZone.py:473-495`, 25 lines; a one-line delegate to `workbench.trim_stable_at`.
- **Exports:** `src/tanglepack/loom/__init__.py:26,39` and `src/tanglepack/__init__.py:62,108`.
- **Callers:** none in src, tests, scripts or stress.
- **Size:** src −29.

### 4. `_get_iterate` + `_get_refined_point` — delete, rewrite one test

- **What they are.** `ManifoldMachine.py:815-845` and `:847-851`, 38 lines. `_get_iterate`'s only caller is `:837`, inside `_get_refined_point`.
  - `_get_preiterate` and `_cache_preiterate` stay: `_refine_layer` uses them at `:739,740,805`.
- **The test.** `tests/unit/numerics/test_refinement.py:33-49` calls `_get_refined_point` at `:46`.
- **Rewrite, verified in `proto_refine.py`:**
  - Set `machine.area_cutoff = 0.0`, then call `next_layer = machine._refine_layer([(p0, p1)], manifold, viewer, None, set())`. It returns `[(p0, new), (new, p1)]`.
  - On the k10 7-step fixture the new point is bit-identical to the old one: same cdist, max coordinate difference 0.0.
  - The assertion `lo < cdist < hi` keeps its full strength. Test diff +5/−2.
- **Stress file.** `stress/numerics/e02_gpu_end_to_end.py:172` lists `"_get_iterate"` in a set of profiler category names. Harmless, but rubric 2's grep (which covers `stress`) would hit it (Q5).
- **Size:** src −38.

### 5. Registry networkx graph, `build_intersection_graph` / `visualize_intersection_graph`, `graphviz.py` — delete all

Nothing computed reads any of this, and no oracle does. The law tier reads only `dual.graph()` at `tests/helpers/laws.py:2716`, which is the DualGraph's own graph and stays.

**What goes:**
- `graphviz.py`, the whole file, 358 lines. It holds `build_intersection_graph` (`:38`) and `visualize_intersection_graph` (`:62`). Also its export in `numerics/__init__.py:35,67`.
- In `TangleWorkbench.py`, 51 lines:
  - the section header and both methods (`:2047-2094`; `:2049`, `:2060`)
  - `import networkx` (`:41`) and `from . import graphviz` (`:43`)
  - docstring mentions at `:15` and `:173`
- In `IntersectionRegistry.py`, about 176 lines:
  - `graph` (`:530-569`), `_rebuild_adjacency_edges` (`:571-626`), `_consecutive_on_branch` (`:628-656`)
  - the `_graph*` fields (`:178-182`) and docstring lines `:91,120-122`
  - the node adds in `add` (`:236-246`) and `reindex_from` (`:350,363-373`)
  - the edge in `register_iterate` (`:674-682`) and `remove_iterate` (`:693-696`)
  - `import networkx` (`:23`)

**Readers and what happens to them:**
- **`tests/unit/numerics/test_map_step_and_graph.py:201-298`:** five tests that exist only for the graph. Delete them, plus the module-docstring paragraph at `:8-11`. Diff −105/+2. The tests:
  - `test_stable_edges_join_consecutive_crossings_on_one_branch`
  - `test_unstable_edges_are_exactly_the_bridges`
  - `test_registry_graph_falls_back_to_per_branch_unstable_order`
  - `test_registry_graph_takes_the_bridges_it_is_given`
  - `test_build_intersection_graph_does_not_mutate_the_registry_graph`
- **`tests/unit/numerics/test_iterate_inference.py:85,95,147,153`:** these lines read graph edges. B8 added them (commits `a001961`, `ae18859`).
  - `:85` becomes `registry.iterate_table[z, 1] == y`.
  - `:147` becomes `registry.iterate_table[x, 1] == y`.
  - `:95` and `:153` are dropped; the table assertions on the line before (`:90`, `:152`) already say the same.
  - Rename `test_remove_iterate_forgets_the_link_its_twin_and_its_edge` to `…_the_link_and_its_twin`.
  - Nothing about the iterate relation is weakened.
- **`tests/unit/numerics/test_generation_and_caches.py:116`:** goes with `:113-115` (see row 6).
- **Scripts:** `scripts/henon_new_intersection_structures.py:76,95` and `scripts/henon_period_3_new_intersection_structures.py:102`. Their graph sections are dropped (Q1).

**Side effect (verified):** after the deletion, `import tanglepack` no longer loads `networkx`. The topology layer imports it lazily, and `pyproject.toml:19` still lists it as a dependency. Only import time changes.

### 6. Registry query methods — delete

| Name | Definition (`IntersectionRegistry.py` unless noted) | Library callers | Tests | Scripts |
|---|---|---|---|---|
| `on_cdist_range` | `:763-784`, plus `:95` and the `add_synthetic` mention `:264` | none | `test_generation_and_caches.py:115` | `henon_new_intersection_structures.py:63`, `henon_period_3_new_intersection_structures.py:89` |
| `from_branch` | `:816-849`, plus `:97` | none | `test_intersection_registry_fixes.py:168-186`, a 3-param test that exists only for it (D4 fix from WP11) | `henon_new_intersection_structures.py:70`, `henon_period_3_new_intersection_structures.py:96` |
| `nearest_by_unstable_cdist` | `:851-873` | none | `test_generation_and_caches.py:114,142` | none |
| `as_forward_array` / `as_backward_array` | `:875-886` | none | none | `henon_new_intersection_structures.py:98-99`, `henon_period_3_new_intersection_structures.py:128-129` |
| Not in the plan: `IterateTable.as_forward_array` / `as_backward_array` | `IterateTable.py:98-125` | only the registry wrappers above | none | none |

- **Size:** registry about −97; IterateTable −29, including its `numpy` / `NDArray` imports, which become unused. The registry's own `Iterable` and `NDArray` imports (`:21`, `:24`) also become unused.
- **Script rewrites (verified, same counts printed):**
  - `on_cdist_range(a, b, stability="stable")` becomes `registry.filter(lambda ix: a <= ix.stable_cdist <= b)`.
  - `from_branch(0, stability="stable")` becomes `registry.filter(lambda ix: ix.is_anchor or (ix.manifold_b_key is not None and ix.manifold_b_key[3] == 0))`.
  - The array calls are dropped.
- **These stay, because they have live callers:**
  - `filter`
  - `from_fixed_point`: used at `tests/unit/numerics/test_grow_until.py:425` and by both scripts
  - `iterate_orbit`: used at `Trellis.py:1457`, `ResonanceZone.py:539` and `test_dual_graph.py:115`

### 7. `classify_bridges` — delete, rewrite its script and its test (Q2)

- **What it is.** `src/tanglepack/loom/TangleSession.py:1618-1634`, 18 lines: a dict comprehension over `classify_bridge` (`:1592`). `classify_bridge` stays; `Topology.py:165` and the blast code use it.
- **Its script is a live recipe.** `scripts/henon_bridge_classification.py:83` runs at `HEAD` and prints a zone × fixed-point table.
  - Rewrite: `{b: session.classify_bridge(b) for b in session.workbench.bridges}`.
  - Verified: the output is identical.
- **Its test should be rewritten, not deleted** (the plan said delete it).
  - `tests/unit/loom/test_resonance_zone.py:111-122` really tests that each in-zone bridge's midpoint lies in the zone and that `classify_bridge` agrees.
  - Change only line `:115` to the same comprehension.

### 8. `terminated_early` — delete

- **What it is.** `src/tanglepack/loom/Blast.py:201-204`, 5 lines.
- **Tests.** `tests/unit/loom/test_blast_budget.py:54,95`, plus docstring mentions at `:6-7` and `:90-91`. Diff +3/−6.
  - `:54` follows from `:53`, which asserts `stop_reason == "point_budget"`.
  - `:95` is a tautology.
- **Script.** `scripts/henon_blast_period_3.py:123` reads it, but that script is dead at `HEAD` (verified). It raises at `:95` on `T1.set_strong_pip(12)`: `ValueError: Intersection 12 is not a strong-pip candidate. Candidates: [31, 32, 33]` (Q3).
- **CLAUDE.md:154** mentions it; that is B12's to fix.

### 9. `partition_element_for`

Belongs to B9; not touched here.

### 10. Oracles stay

`euler_characteristic`, `Arrangement.preimage_of`, `StableNode.other_face` and `preimage_bridges` stay. The prototype diff touches no `tests/helpers/laws.py`, `tests/cases.py` or `tests/invariants/`, and every law test passes in the copy.

### 11. `topology.__all__` (104 names, `topology/__init__.py:150-251`)

I scanned imports and `topology.X` reads in src, tests, scripts and stress (with an AST scan), plus notebooks (with a regex). **Only 15 names are ever imported through the `tanglepack.topology` namespace from outside the subpackage.**

Everything else is reached through its submodule:
- loom uses `from ..topology.X import …`;
- tests use `importlib` module aliases (`SD`, `DG`, `family_module`);
- `plot_*` is reached via `from tanglepack.topology import plotting`.

Nothing reads `__all__`, and nothing does `import *`.

- **(a) Must stay importable from the namespace (15):**
  - The ten re-exported at `src/tanglepack/__init__.py:45`: Trellis, TrellisBranch, Arrangement, Arc, Region, Hole, PartitionInterval, PseudoneighborPair, StablePartitionResult, StrongPipResult.
  - `side_after` (`tests/unit/topology/test_side_after.py:9`)
  - `letter` (`tests/unit/topology/test_symbolic_dynamics.py:869`)
  - `bridge_row_diagnosis` (`tests/regression/test_coarse_blasts_rows_by_signs.py:23`)
  - `check_holes_share_bridge_side` (`tests/regression/test_open_p3_deep_runs.py:33`)
  - `plot_stable_partition` (`scripts/henon_pseudoneighbors_nested.py:29`)
- **(b) Used by loom through submodules — the natural "stage" exports (24):**
  - Stage objects: MinimalTrellis, minimal_trellis, HomotopyPartition, IteratedHomotopyPartition, BridgeClass, BridgeClassEntry, BridgeClassTable, DualGraph, Walk, ClassDynamics, SymbolicDynamics, inert_letters, symbolic_dynamics, bridge_classes, zone_label.
  - Plotting: plot_dual_graph, plot_minimal_trellis, plot_walk, plot_bridges_by_class, plot_transition_graph, plot_itinerary_table, plot_dual_graph_cartoon, CartoonLayout, ZoneArcOrder.
- **(c) Only tests, scripts or stress use them, through submodules (47):**
  - Dual graph, naming and walks: StableNode, image_chain, PartitionFamily, BridgeMember, FaceNode, ElementName, ElementNaming, ElementLanding, WalkSearch, WalkStep, holds_markers, land_element, shortest_walks, trellis_itinerary, walk_chords.
  - Symbolic dynamics and results: EvidenceDiagnosis, Symbol, translate_itinerary, ElementRef, canonical_corners.
  - Pips and pseudoneighbors: is_strong_pip, compute_pseudoneighbors, extend_pseudoneighbor_trajectories.
  - Rows and bridges: bridge_for_pair, bridge_row_violation, bridge_side_violations, side_diagnoses, check_bridge_rows_consistent, geometry_flags, class_sort_key, anchor_outward_key, oriented_class, owns_cdist.
  - Plotting helpers: dual_graph_legend_handles, face_point, minimal_trellis_legend_handles, stable_node_point, walk_legend_handles, dual_graph_cartoon_layout, dual_graph_cartoon_legend_handles, ZoneLayout, dual_graph_zone_layout.
  - Partition: partition_stable_manifold, propagate_reference_holes, punch_holes, row_of_end, rows_of_bridge.
- **(d) No name import from outside topology (18):**
  - Used only inside topology: Cut, Marker, MembersDisagree, WalkReason, lobe_image_markers, ImageEvidence, RefinedClass, free_reduce, refine, classify_strong_pips, RowDiagnosis, SideDiagnosis, element_sort_key.
  - Not used by name anywhere: SIDE_LETTER, WalkStatus, registered_evidence, RowEnd, span_contains.
  - Some of (d) are read through module aliases (`SD.free_reduce`, `DG.registered_evidence`). Those reads survive any `__all__` trim.

**What this means for the ≤ 70 bound:**
- (a)+(b) is 39 names, so ≤ 70 is easy and needs no re-derived bound.
- The only hard constraint is that the (a) names stay in the import block.
- Size: 104 → 70 is at least −34 lines; 104 → 39 is at least −65. Dropping the matching import lines too roughly doubles that.

## Line counts (measured in the prototype)

| Row | src | tests | scripts |
|---|---|---|---|
| grow_x_times | −20 (+1 doc) | 0 | +16/−8 |
| _get_refined_point + _get_iterate | −38 | +5/−2 | 0 |
| trim_stable_at_intersection | −29 | 0 | 0 |
| graph: registry −176, workbench −51, graphviz.py −358, numerics init −2 | −587 | −105+2; +4/−6 (iterate_inference); −1 | in the two demo scripts |
| registry queries −97, IterateTable −29 | −126 | −19; −3 | two demo scripts +8/−37 (with the graph row) |
| classify_bridges | −18 | ±1 | ±1 |
| terminated_early | −5 | +3/−6 | (dead script, Q3) |
| **Total** | **−822 (+1)** | **−128 (+15/−143)** | **−21 (+25/−46)** |

- **src + tests = −950**, not counting the `__all__` trim (at least −34 more at 70 names) or 2–4 newly unused import lines.
- **Docstrings still to reword** (net 0 lines):
  - `IntersectionRegistry.py`: `:184-189`, `:201`, `:264`, `reindex_from`, `register_iterate`, `remove_iterate`
  - `TangleWorkbench.py:173`
  - `ManifoldMachine.py:128`
- **Honest targets:** ≤ −950 with no `__all__` trim, about ≤ −980 with the trim to 70. The plan's ≤ −900 holds either way.

## What I verified, and how

| Claim | How |
|---|---|
| Every caller list above | `grep -rnw` over src, tests, scripts and stress (`*.py`) plus notebooks; AST scan for `__all__` (`<scratchpad>/importers.py`, `all_users.py`) |
| Nothing computed reads a deleted name | Full suite in the prototype: 1215/47/5/1 xfailed vs `HEAD` 1223/47/5/1; −8 = the deleted tests; only xfail is C3, no XPASS |
| `release_manifold` is live | Call probe on `build_k28()`, plus `test_bridge_identity.py:257` |
| The refinement rewrite gives the exact same point | `proto_refine.py` |
| Scripts | At `HEAD`, 7 of 8 run (about 1 s each) and `henon_blast_period_3.py` fails at `:95`. In the prototype, the 6 rewritten scripts and `tangle_workbench_test.py` run with identical output, except the removed `Graph:` line |
| Imports | In the prototype, from `/tmp`: `import tanglepack` and `from tanglepack.{topology,numerics,loom} import *` succeed |
| Rubric 2 grep in the prototype | Clean except the e02 string (Q5), `henon_blast_period_3.py:123` (Q3), and two docstrings the prototype left unreworded (`IntersectionRegistry.py` ~`:239`, `TangleWorkbench.py` ~`:170`) |

Not verified by me:
- **The 63-item probe.** Harness is `<scratchpad>/b6b/run.sh` + `probe.py`. It hardcodes `REPO` to the main tree, so it needs re-pointing at a worktree before use.
- **The minutes-long CLAUDE.md scripts.** None imports a deleted name, but they write into the repo's `figures/`, so I did not run them on the main tree.

## Where the plan's table is wrong or out of date

1. `release_manifold` sits at `:648`, not `:711`, and is live: KEEP.
2. `grow_x_times` has 3 scripts and 8 lines, not 10. The scripts cannot use `grow_n_times`; use the `new_grow_manifold` loop. Two notebooks also call it.
3. The `classify_bridges` test should be rewritten, not deleted, and its script caller is a live recipe.
4. Five `registry.graph()` reads outside the graph test file (added in B8) move to `iterate_table`.
5. New candidates not in the table: the `IterateTable` array methods and the e02 string.
6. `__all__` has 104 names, not 101.
7. The file list grows to about 24 files (the request said about 15). It must add:
   - `IterateTable.py`
   - `test_intersection_registry_fixes.py`
   - the six rewritten scripts
   - and, depending on Q3 and Q5, `scripts/henon_blast_period_3.py` and `stress/numerics/e02_gpu_end_to_end.py`

## Questions for the author (my default in each)

- **Q1.** `scripts/henon_new_intersection_structures.py` and `scripts/henon_period_3_new_intersection_structures.py` exist to show off the API being deleted. Rewrite them or delete them? Default: **rewrite** — use `filter` for the two queries and drop the graph and array sections (verified to run).
- **Q2.** `classify_bridges` has a live script caller. Delete the method and rewrite its two callers, or keep the method? Default: **delete and rewrite**.
- **Q3.** `scripts/henon_blast_period_3.py` has been dead since B3, and its `:123` reads `terminated_early`. Default: **(a)** rewrite only `:123` to `stop_reason == "frontier_empty"`, and report that the script still fails at `:95`. The alternative, (b), is to leave it and exempt it from rubric 2.
- **Q4.** `from_branch`'s test is the only test of the D4 anchor rule written as a registry query. It goes with the method; the rule lives on in Trellis and Arrangement, which I did not trace. Default: **delete**.
- **Q5.** The e02 string is a profiler category name, not a call. Default: **drop the string**.
- **Q6.** The two notebooks that call `grow_x_times`. Default: **leave them** (plan Deferred), and limit rubric 2's grep to `*.py`.
- **Q7.** `IterateTable`'s array methods are orphaned by row 6. Default: **delete them in the same commit**.
- **Q8.** Trim `__all__` only, or the import block too? Default: trim `__all__` to (a)+(b) plus the plotting helpers (about 48 names) and **keep every import**, so nothing can break; dropping imports fits B12.

## Recommendation

- Delete as tabled, one commit per row.
- KEEP `release_manifold` and every oracle.
- Rewrite, don't delete, the refinement and classify tests.
- Rewrite the six scripts as prototyped.
- Extend the file list as above.
- State rubric 1 as ≤ −950 (src + tests, before the `__all__` trim) and rubric 2's grep as `*.py` only.

## Out of scope

- `scripts/henon_period_3_new_intersection_structures.py` prints `Forward chain from 3: [3, 3]`. Id 3 is the period-1 anchor, so that is correct, but the script labels it as period 3.
- `scripts/henon_new_intersection_structures.py:64` prints "[2, 8]" for a `[2, 6]` query.
- `IterateTable.forward_depth` / `forward_chain` read direct entries only (noted since B3).
- CLAUDE.md `:136` (`graphviz`, networkx views) and `:154` (`terminated_early`) are B12's.

Evidence is in `<scratchpad>`: `suite_base.log`, `suite_proto.log`, `proto/` (the prototype), `delete*.py`, `proto_refine.py`, `probe_release.py`, `importers.json`, `all_users.json`, `run_*.log`, `p_*.log`, `rw_*`. Knowledge files updated: `.refactor/knowledge/numerics.md`, `loom.md`, `topology.md`, `testing.md`.
