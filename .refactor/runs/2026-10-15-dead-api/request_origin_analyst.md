# Request: B10 — `dead-api`

*From the author, 2026-10-15 ("continue"). Input for one `tangle_refactor` run. The thirteenth slice of `docs/alignment_plan_2026-10-09.md` (§B, "B10 — `dead-api`", ~line 810); B0, B1, B6a, B8, B4a, B4b, B2, B3, B5, B7, B7b and B6b are merged (main `65edec4`; runs archived under `.refactor/runs/2026-10-1{2,3,4,5}-*/`).*

## The question

**Which API has no user, and does deleting it change anything?** The plan's B10 table (2026-10-09, callers verified by grep THEN) lists the zero-user and test-only names with a disposition each: delete the zero-user ones with their exports, tests and doc mentions; keep test-only names that serve as ORACLES for the law tier (`euler_characteristic`, `Arrangement.preimage_of`, `StableNode.other_face`, `preimage_bridges`; plan ~`:59-62`) or move them into `tests/helpers`; keep `refine_manifold` (27 script users); leave `partition_element_for` to B9; trim `topology.__all__` (104 names at `HEAD`) toward the stage exports plus plotting (B12 finishes). Behaviour: none, by construction — nothing deleted is on any computed path.

The table is six days and twelve slices old. A pre-check at `HEAD` shows callers have MOVED since it was written, so every row must be re-derived before anything is deleted:

| Name (plan row) | Plan said | At `HEAD` (grep, counts of lines: src / tests / scripts / stress / docs) |
|---|---|---|
| `grow_x_times` | scripts only (10) | 2 / 0 / 8 / 0 / 2 |
| `release_manifold` | workbench-internal | 4 / 0 / 0 / 0 / 4 |
| `trim_stable_at_intersection` | none, exported twice | 5 / 0 / 0 / 0 / 2 |
| `_get_iterate` | none | 2 / 0 / 0 / **1** / 3 |
| `_get_refined_point` | `test_refinement.py` | 1 / 2 / 0 / 0 / 11 |
| `visualize_intersection_graph` | none | 6 / 0 / **1** / 0 / 4 |
| `build_intersection_graph` + `IntersectionRegistry.graph` / `_rebuild_adjacency_edges` / `_consecutive_on_branch` | `test_map_step_and_graph.py` | 7 / 3 / **2** / 0 / 10; **`registry.graph()` is ALSO read by `tests/unit/numerics/test_generation_and_caches.py:116` and `test_iterate_inference.py:85,95,147,153`** (B8 or earlier added them) |
| `on_cdist_range`, `from_branch`, `nearest_by_unstable_cdist`, `as_forward_array` / `as_backward_array` | tests only / none | 3/1/**2**, 2/2/**2**, 1/2/0, 3/0/**2**, 3/0/**2** |
| `classify_bridges` | `test_resonance_zone.py:115` | 1 / 1 / **1** / 0 / 2 (the test is `tests/unit/loom/test_resonance_zone.py`) |
| `terminated_early` | `test_blast_budget.py` | 1 / 4 / **1** / 0 / 2 + one CLAUDE.md mention (`BlastResult` row: "`terminated_early` IS `frontier_empty`") |
| `partition_element_for` | B9, not here | 3 / 3 / 0 / 0 / 10 (`tests/unit/topology/test_partition_elements.py`) |
| `graphviz.py` | likely delete (358 lines) | 358 lines; `src` 10 mentions, one CLAUDE.md table row |

`docs/` counts are the uncommitted plan/audit/knowledge files and are NOT touched by the slice (B12). The `scripts/` and `stress/` hits are what the analyst must resolve: a script that calls a deleted name is updated in the slice (the plan's `grow_x_times` → `grow_n_times` rule generalised: every script caller is rewritten to the surviving API, or the script is listed as dead with the author asked), never left broken.

## What "done" means

Exactly slice B10 of the alignment plan, re-derived at `HEAD`:

1. **One commit per table row**, each deleting the name, its exports (`numerics`/`loom`/`topology` `__init__`), its tests (a test that exists only to exercise the deleted name goes with it; a test that reads the name as a convenience is rewritten to the surviving API — e.g. the five `registry.graph()` readers in `test_generation_and_caches.py` and `test_iterate_inference.py` read the iterate table or `iterate_table[id, n]` instead), its script callers (rewritten), and its mention in the module Dev Notes / docstrings. `CLAUDE.md` table rows are B12's (the tree's CLAUDE.md carries uncommitted edits; the branch never touches it).
2. **The registry's networkx graph** (`IntersectionRegistry.graph`, `_rebuild_adjacency_edges`, `_consecutive_on_branch`, ~230 lines at `:416-545` then) and `graphviz.py` (358 lines) go, unless the analyst finds a computed path or an oracle reading them (then the planner says what stays and why). `tests/helpers/laws.py:2716` reads `dual.graph()` — that is the DUAL GRAPH's networkx view, not the registry's; it stays.
3. **Oracles stay.** `euler_characteristic`, `Arrangement.preimage_of`, `StableNode.other_face`, `preimage_bridges` and anything else the law tier reads (`laws.py` is NOT in the diff, rubric 4); if the planner moves one into `tests/helpers`, the law tier's read must not change.
4. **`topology.__all__`** trimmed from 104 toward the stage exports + plotting: the analyst lists which of the 104 names have a user outside `tanglepack.topology` (tests, scripts, loom, stress) and the planner keeps exactly those plus the public stage classes; rubric: `len(topology.__all__) ≤ 70` (re-derived if the count of names with users is higher — say so). `from tanglepack.topology import *` and `import tanglepack` succeed in a fresh interpreter.
5. **`release_manifold`:** delete if the workbench path that calls it is itself unused (the analyst verifies the four src mentions), else keep private. **`_get_refined_point`:** the one kernel test in `test_refinement.py` reads the vectorised kernel's output for the same pair instead. **`terminated_early`:** `test_blast_budget.py`'s four reads become `stop_reason == "frontier_empty"`. **`classify_bridges`:** deleted with its one test unless the script caller is a live recipe (then rewritten).

Files (the plan's list, re-derived): `src/tanglepack/numerics/ManifoldMachine.py`, `IntersectionRegistry.py`, `TangleWorkbench.py`, `graphviz.py` (deleted), `src/tanglepack/loom/ResonanceZone.py`, `Blast.py`, `TangleSession.py`, the three `__init__.py`, the tests named above (`test_refinement.py`, `test_map_step_and_graph.py`, `test_generation_and_caches.py`, `test_iterate_inference.py`, `test_resonance_zone.py`, `test_blast_budget.py`), and every script/stress caller the analyst finds. Flagged larger: deletions across ~15 files, NO logic change.

**Behaviour.** None: every law case, deep run and series step gives identical words, matrices, counts, `is_reliable` (rubric: the B6b/B7b probe scripts in the session scratchpad, rerun, show 0 differences on all 63 items); every script listed in CLAUDE.md "Run a scripted example" still runs (the ones that take minutes may be smoke-run with reduced parameters; say which).

**Not doing.** Session API shrink (B9: `partition_element_for`, the fan-outs, plot delegates); the plotting split (SIMP-11); anything the law tier reads; CLAUDE.md rows (B12); `refine_manifold`; the notebooks' private copies of the Hénon map (plan Deferred list).

**Rubric** (the plan's four, re-derived by the planner):
1. Suite passes; `pytest -q -rxX` lists exactly one xfail (C3), no XPASS; `git diff --stat` touches only the listed files; **net ≤ −900 lines** (src + tests; the planner re-derives the figure from the analyst's line counts and states it).
2. For every deleted name, `grep -rn "\b<name>\b" src tests scripts stress` → nothing (docs are B12's).
3. `import tanglepack`, `from tanglepack.topology import *`, `from tanglepack.numerics import *`, `from tanglepack.loom import *` succeed in a fresh interpreter from `/tmp`; `len(tanglepack.topology.__all__) ≤ 70` (or the re-derived bound).
4. The law tier is unchanged: no `tests/helpers/laws.py`, `tests/cases.py` or `tests/invariants/` edit in the diff; base-vs-branch probe equal on all 63 items.
5. Every script in CLAUDE.md "Run a scripted example" and every `scripts/*.py` that imported a deleted name runs (smoke) on the branch.

## Where to start

`docs/alignment_plan_2026-10-09.md` (B10 `:810-845`; the oracle note `:59-62`; §C `:1117`, `:1161`), `tests/conftest.py` (the test policy docstring; the "dead API" list if still there), `src/tanglepack/numerics/graphviz.py`, `src/tanglepack/numerics/IntersectionRegistry.py` (`graph` `:530`, the networkx block), `src/tanglepack/numerics/TangleWorkbench.py` (`build_intersection_graph`, `visualize_intersection_graph`, `release_manifold` callers), `src/tanglepack/numerics/ManifoldMachine.py` (`grow_x_times`, `_get_iterate`, `_get_refined_point`), `src/tanglepack/loom/ResonanceZone.py` (`trim_stable_at_intersection`), `src/tanglepack/loom/Blast.py` (`terminated_early`), `src/tanglepack/loom/TangleSession.py` (`classify_bridges`), the three `__init__.py`, `scripts/` (every caller), `stress/numerics/` (the `_get_iterate` mention), the tests named above.

## Binding context

- Test policy (`tests/conftest.py`, 2026-10-05): no case-specific pins; a test that exists only to exercise a deleted name is deleted with it, a test that USES the name incidentally is rewritten to the surviving API with the same assertion strength (never weakened); `KNOWN_ISSUES` is C3 alone and may only shrink; the law tier is untouched.
- `CLAUDE.md` invariants are physics and none is touched: the registry stays THE single store of crossings; `iterate_table[id, n]` and `Trellis.iterate` are the lookup paths (the networkx graph was a VIEW over them); `numerics`/`topology` never import `loom`.
- Author's rules in force: a script that is dead on main (B3 found `scripts/henon_blast_period_3.py` dead: pinned pip id no longer a candidate) is reported, not silently fixed or deleted — the planner lists dead scripts as a question for the author with a default of "leave, note in the report".
- The analyst stage is NOT skipped: the analyst re-derives EVERY table row at `HEAD` (= `65edec4`) with exact file:line callers in src, tests, scripts, stress (docs excluded), decides per row "delete / rewrite callers / keep (oracle) / keep (user)", measures the lines each deletion removes (so the planner can state the net figure honestly), lists the 104 `topology.__all__` names with their external users, checks `release_manifold`'s path, runs the scripts that import any name to be deleted to see which are live, and reports what differs from the plan's table. It does not re-plan the slice.
- The working tree may carry uncommitted documentation changes (`docs/`, `.refactor/`, `CLAUDE.md`); the change branches from `HEAD` and never touches them.
