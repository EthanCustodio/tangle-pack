# Change — B10 `dead-api` (coder, iteration 1, 2026-10-15)

## Where

- **Base commit:** `65edec4` (the author's `HEAD`, `main`).
- **Branch:** `refactor/dead-api`, worktree `/home/dezhu/Development/TanglePack/tangle-pack-dead-api`, head `517aa01`.
- **Base probe:** detached worktree `/home/dezhu/Development/TanglePack/tangle-pack-dead-api-base` at `65edec4`.
- The main tree was not edited. `git status` there is the same as at the start, plus this artifact.
- Scratch evidence: `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b10/`
  - `suite_217920a.txt`, `suite_22feb61.txt`: the C3 and C7 suite runs
  - `scripts_base/`, `scripts_branch/`: script outputs for rubric 16 and 17

## Commits

Eight commits, C1–C8, in the plan's order. Every message is `Dead API: <name> (B10)`, has a short body, and ends with the two attribution lines.

| # | Hash | Subject | Files (+/−) |
|---|---|---|---|
| C1 | `1867f41` | grow_x_times | ManifoldMachine.py +1/−21; henon_basic_stable_and_unstable.py +4/−2; henon_initial_intersections.py +4/−2; henon_period_3_manifolds.py +8/−4 |
| C2 | `88d5a75` | _get_iterate and _get_refined_point | ManifoldMachine.py −37; test_refinement.py +5/−2; e02_gpu_end_to_end.py +1/−1 |
| C3 | `217920a` | trim_stable_at_intersection | ResonanceZone.py −26; loom/__init__.py −2; tanglepack/__init__.py −2 |
| C4 | `7eaf64d` | the registry graph and graphviz.py | graphviz.py −358 (deleted); IntersectionRegistry.py +6/−187; TangleWorkbench.py +9/−60; numerics/__init__.py −2; test_map_step_and_graph.py +1/−106; test_iterate_inference.py +5/−5; test_generation_and_caches.py −1; the two `*_new_intersection_structures.py` scripts −14 / −11 |
| C5 | `3a9db04` | registry query methods and IterateTable arrays | IntersectionRegistry.py +1/−98; IterateTable.py −31; test_generation_and_caches.py −3; test_intersection_registry_fixes.py −21; the two demo scripts +4/−6 each |
| C6 | `1348e8c` | classify_bridges | TangleSession.py −18; test_resonance_zone.py +1/−1; henon_bridge_classification.py +2/−2 |
| C7 | `22feb61` | terminated_early | Blast.py −5; test_blast_budget.py +5/−7; henon_blast_period_3.py +1/−1 (line 123 only) |
| C8 | `517aa01` | topology.__all__ to 48 names | topology/__init__.py −56 |

`git diff --stat 65edec4` gives **27 files changed, 63 insertions(+), 1097 deletions(-)**. These are exactly the 27 paths in the plan's table.

## Suite

Command, run in the worktree: `MPLBACKEND=Agg PYTHONPATH=src ../tangle-pack/env/bin/python -m pytest -q -rxX -p no:cacheprovider`

| Commit | Result |
|---|---|
| C3 `217920a` (scratch worktree) | 1223 passed, 47 skipped, 5 deselected, 1 xfailed |
| C4 `7eaf64d` | 1218 passed, 47 skipped, 5 deselected, 1 xfailed (−5: the graph tests) |
| C5 `3a9db04` | 1215 passed, 47 skipped, 5 deselected, 1 xfailed (−3: `from_branch` × 3) |
| C7 `22feb61` (scratch worktree) | 1215 passed, 47 skipped, 5 deselected, 1 xfailed |
| C8 `517aa01` (branch head) | **1215 passed, 47 skipped, 5 deselected, 1 xfailed in 199.80s** |

- In every run the only xfail is `test_orientation_reversing_case_builds[orientation_reversing]` (C3: INV-11), and there is no XPASS.
- C1, C2 and C6 were not run as whole suites. They were covered by targeted runs:
  - `test_refinement.py` at C2: 1 passed
  - `test_resonance_zone.py` at C6: 13 passed
  - `test_blast_budget.py` at C7: 6 passed
- The suites at C3 and C7 bracket C1–C2 and C6, since each of those commits only deletes code with no caller left.

## Behaviour changes

Nothing on a computed path changed.

- **`import tanglepack` no longer loads `networkx`.** Checked: `'networkx' in sys.modules` is `False` on the branch and `True` on base. `networkx` stays in `pyproject.toml`, and topology imports it lazily.
- **Two demo scripts no longer print the graph line.** `henon_new_intersection_structures.py` stops printing `Graph: 8 nodes, 17 edges` and stops calling the graphviz `plt.show()`. `henon_period_3_new_intersection_structures.py` stops printing `Graph: 207 nodes, 582 edges`.
- **`henon_blast_period_3.py` is still dead (Q3 (a)).** Line 123 now reads `result.stop_reason == 'frontier_empty'`. The script still raises before it gets there, on base and branch alike:
  - where: `File ".../scripts/henon_blast_period_3.py", line 95`
  - error: `ValueError: Intersection 12 is not a strong-pip candidate. Candidates: [31, 32, 33]. Run classify_strong_pips() first, or choose one of these.`

  It was not revived; that is outside this run.
- **`topology.__all__` drops 56 names.** They are still importable by name; only `from tanglepack.topology import *` now brings in 48 names instead of 104.

## Deviations from the plan

All are small and inside the file list.

1. **C3 also drops a now-unused `TYPE_CHECKING` import.** `from ..numerics.BaseManifold import BaseManifold` in `ResonanceZone.py` was used only by the deleted function's signature.
2. **The unused `Iterable` import moved from C5 to C4.** The plan listed it under C5, but the graph code deleted in C4 was its last user, so I removed it there. `NDArray` went in C5 as planned.
3. **`remove_iterate` no longer reads `target_id`.** That local existed only to name the graph edge. `IterateTable.remove` still raises `KeyError` on a missing link, as the docstring says.
4. **`henon_bridge_classification.py` changes in 2 lines, not 1.** The trailing comment `# {bridge: ResonanceZone | None}` moved onto the section header line, so the comprehension line stays at 84 columns. Rubric 12 is about the test file, and that diff is one line.
5. **Docstrings that the plan or findings did not list were also reworded,** so that no sentence mentions a registry graph:
   - `IntersectionRegistry.generation`: "the graph's adjacency edges, the rank maps" becomes "the rank maps, the candidate index"
   - `TangleWorkbench.generation` (about `:406`): "the registry graph's adjacency edges, the rank maps" becomes "the registry's rank maps"
   - The `TangleWorkbench` module docstring now says "Two collaborators", since the graphviz bullet went.
   - `add_synthetic` now points at `filter()` instead of `on_cdist_range()`.
6. **`test_blast_budget.py` module docstring.** The plan said to reword it to speak of `stop_reason`, but the prototype only deleted the sentence. I wrote the sentence: "``BlastResult.stop_reason`` says which stop ended a blast."
7. **`test_intersection_registry_fixes.py`:** the trailing blank lines left after deleting the last test were stripped.

## Rubric items I ran

- **R1 (suite):** 1215 / 47 / 5 deselected / 1 xfailed C3, no XPASS (table above).
- **R3 (law-tier pins):** `git diff --stat 65edec4 -- tests/helpers/laws.py tests/cases.py tests/invariants/` is empty. I did not run `pytest tests/invariants` on both worktrees separately; the full-suite runs cover it on the branch.
- **R4 (diff stays in the file list):** exactly the 27 listed paths. Nothing under `CLAUDE.md`, `docs/`, `.refactor/`, `artifacts/`, `notebooks/` or `figures/` changed, and the main-tree status is unchanged.
- **R5 (one commit per row):** eight commits, one per row. Suite runs at C3, C4, C5, C7 and C8 all pass.
- **R6 (net lines over `src tests`):**

  | Point | `--shortstat` | Net | Target |
  |---|---|---|---|
  | at C7 | 18 files, +35 / −992 | **−957** | ≤ −950 |
  | whole branch | 19 files, +35 / −1048 | **−1013** | ≤ −1000 |

- **R7 (no deleted name in `*.py`):** the grep returns nothing (exit 1).
- **R8 (no graph in numerics):** `grep -rn "networkx\|_graph" src/tanglepack/numerics` returns nothing, and `graphviz.py` does not exist.
- **R9 (imports from outside the repo):** run from `/tmp` with `PYTHONPATH=<worktree>/src`, in a fresh interpreter each time.
  - All five statements succeed on both worktrees: `import tanglepack`, `from tanglepack.numerics import *`, `from tanglepack.topology import *`, `from tanglepack.loom import *`, and `from tanglepack.topology import StableNode, free_reduce, MembersDisagree, plotting`.
  - The `networkx` check prints `False` on the branch and `True` on base.
- **R10 (`__all__`):** `len(__all__) == 48`, every name resolves, and all 15 (a) names are present. Base has 104.
- **R11 (refinement test keeps its teeth):**
  - The test sets `machine.area_cutoff = 0.0` and calls `_refine_layer([(p0, p1)], manifold, viewer, None, set())`.
  - It asserts `lo < float(new_point.cdist) < hi`, strict, with no tolerance. It also asserts that both child pairs share the new point.
  - Mutation check: swapping `new_point` for `p0` in a temporary copy makes the test fail (1 failed). The copy was then deleted.
- **R12 (classify test):** the `test_resonance_zone.py` diff is the one line of the comprehension. The midpoint-in-zone and `classify_bridge` assertions are untouched.
- **R13 (iterate links):** at the old `:85` and `:147`, `test_iterate_inference.py` now asserts `registry.iterate_table[z, 1] == y` and `registry.iterate_table[x, 1] == y`. The twin and `unresolved_links` assertions stay. The test is renamed `test_remove_iterate_forgets_the_link_and_its_twin`.
- **R14 (blast stop reasons):** `test_blast_budget.py` still asserts `second.stop_reason == "point_budget"`. The no-budget test still asserts `stop_reason in {"max_iterations", "frontier_empty"}` and that the empty filter gives `"frontier_empty"`.
- **R16 (rewritten scripts print the same):** the seven scripts were run with `MPLBACKEND=Agg` on both worktrees. Normalised for paths, `HH:MM:SS`, `0x…` and `UserWarning` line numbers. The only differences:

  | Script | Difference |
  |---|---|
  | `henon_new_intersection_structures.py` | the graphviz `plt.show` UserWarning (2 lines) and `Graph: 8 nodes, 17 edges` |
  | `henon_period_3_new_intersection_structures.py` | `Graph: 207 nodes, 582 edges` |
  | the other five | identical |

- **R17 (dead script dead the same way):** `henon_blast_period_3.py` raises the same `ValueError` at line 95 on base and branch, and its diff is 1 line.
- **R18 (long scripts compile):** `py_compile` passes for:
  - `henon_symbolic_itineraries*.py` (4 files)
  - `show_figure.py`
  - `stress/numerics/e12_reverify.py`
  - the two untracked `henon_period3_*_sweep.py`, compiled from the main tree (they are not on the branch)

## Left to the test engineer

- **Rubric 15:** the 63-item probe, base against branch. The harness is `<scratchpad>/b6b/{run.sh,probe.py,summ.py}`; re-point `REPO` and the output dir to `<scratchpad>/b10/`.
- **Rubric 2:** the collected-id diff (`--collect-only -q`) between base and branch. Expected:
  - gone: the five graph tests in `test_map_step_and_graph.py` and `test_from_branch_holds_the_anchor_on_every_branch[None|unstable|stable]`
  - renamed: `test_remove_iterate_forgets_the_link_its_twin_and_its_edge` → `…_the_link_and_its_twin`
- **Rubric 3, second half:** `pytest tests/invariants` on each worktree, comparing counts.
- **Rubric 5:** suite runs at C1, C2 and C6 if a run at every commit is wanted.
