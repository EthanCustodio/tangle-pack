# Outcome — B10 `dead-api` (auditor, 2026-10-15)

## In plain words

- **Your question was:** which API has no user, and does deleting it change anything?
- **Answer, now that the change is done:** every name in the B10 table had no user on a computed path except `release_manifold`, which is live and stays. Everything else is gone, along with its exports, its tests and its script callers. Deleting it changed nothing that the library computes. On all 63 corpus items, the words, matrices, class counts and `is_reliable` match base byte for byte.
- **What you would notice:**
  - `import tanglepack` no longer loads networkx.
  - `from tanglepack.topology import *` brings in 48 names instead of 104. The other 56 can still be imported by name.
  - Two demo scripts no longer print `Graph: N nodes, M edges`.
- **Not done, on purpose:**
  - CLAUDE.md rows and `docs/` (B12)
  - the `topology/__init__` import block (B12)
  - the notebooks (Q6)
  - reviving `henon_blast_period_3.py` (Q3)
- **Before merging:** nothing blocks the merge. The one item worth a look is the stale Sphinx stub, `docs/source/tanglepack.numerics.rst:183-189`, which no earlier report named (see "Open points").

**What to merge:** branch `refactor/dead-api`, head `517aa01`, 8 commits (C1–C8) off `65edec4`, fast-forward. I re-ran the suite at head: **1215 passed, 47 skipped, 5 deselected, 1 xfailed** (only C3 INV-11, no XPASS). Base gives 1223/47/5/1. Net **−1013 lines over src + tests** (+35/−1048), and −1034 over all 27 paths (+63/−1097).

The main tree is unchanged apart from this artifact, `.refactor/decisions.md` and the four knowledge files. The branch worktree is clean.

## Did the branch meet the request's "done"?

| Request item | Done? | How I know |
|---|---|---|
| One commit per row | yes | `git log 65edec4..517aa01`: 8 commits titled `Dead API: <name> (B10)`, in the plan's order. The test engineer ran the suite at all 8 and each passed. |
| Deleted names gone, with exports, tests and script callers | yes | The plan's rubric-7 grep over `*.py` in `src tests scripts stress` finds nothing (test engineer's and reviewer's runs; I re-checked the wider tree, see the notebooks below). The `loom`, `numerics` and top-level `__init__` exports are dropped (diff read). |
| Oracles kept | yes | No diff touches `tests/helpers/laws.py`, `tests/cases.py` or `tests/invariants/`. `euler_characteristic`, `preimage_of`, `other_face`, `preimage_bridges` and `DualGraph.graph()` are untouched. Invariants: 133/46/1 on both trees (test engineer). |
| `release_manifold` kept, and why | yes | It is called from `TangleWorkbench` `compute_intersections`, `_register_bridge` and `clear_bridges`. It does the real work on the direct `clear_bridges` path, which `test_bridge_identity.py:257` pins. |
| `topology.__all__` ≤ 70 | yes, 48 | I ran from `/tmp` against the branch: `len == 48`, every name resolves, `from tanglepack.topology import *` works, and `StableNode, free_reduce, MembersDisagree, plotting` still import by name. |
| No behaviour change on the 63 items | yes | The probe outputs `probe/out_{base,branch}` have the same md5, `bca69aca…`. I recomputed it: 63 lines each. |
| Dead script reported, not fixed | yes | `henon_blast_period_3.py` changes on one line (`:123`). It still raises `ValueError: Intersection 12 is not a strong-pip candidate. Candidates: [31, 32, 33]` at `:95` on both trees (test engineer's logs). |
| Rewritten scripts still run | yes | 7 scripts were run on both trees. The only output differences are the `Graph:` lines and the graphviz `plt.show` warning, as the plan allowed. |

**Asked for in the request but NOT done, and why:**

- **CLAUDE.md rows.** Three rows now describe deleted code:
  - `:133`, the registry row: "per-branch `graph(bridges=)`"
  - `:136`, the `BridgeIterator, infer_iterates, graphviz` row
  - `:154`, `terminated_early` IS `frontier_empty`

  These are B12's: the branch never touches CLAUDE.md. Row `:133` was not on the analyst's or the planner's list. I found it.
- **The `topology/__init__.py` import block** is left whole (Q8). Because of that, `ruff F401` reports 56 findings there, one for each name dropped from `__all__` but still imported. B12 should either trim those imports or add `# noqa`.
- **Notebooks (Q6).** `scripts/manifold_cutting_test.ipynb` (2 calls) and `notebooks/archive/period_three.ipynb` (3 calls, already broken) still call `grow_x_times`. The rubric greps `*.py` only, by design.
- **`docs/`** is B12's.

## What changed

Line counts, base → branch:

| File | Before | After |
|---|---|---|
| `IntersectionRegistry.py` | 988 | 712 |
| `TangleWorkbench.py` | 2260 | 2209 |
| `ManifoldMachine.py` | 956 | 899 |
| `topology/__init__.py` | 251 | 195 |
| `graphviz.py` | 358 | deleted |
| `IterateTable.py` | — | −31 |

The most important hunk is that the registry stops keeping a second copy of the iterate relation. `IntersectionRegistry.remove_iterate` is a good example.

Before:

```python
target_id = self.iterate_table[source_id, n]
self.iterate_table.remove(source_id, n)
self._generation += 1
edge = (source_id, target_id, f"iter_{source_id}_{n}")
if self._graph.has_edge(*edge):
    self._graph.remove_edge(*edge)
```

After:

```python
self.iterate_table.remove(source_id, n)
self._generation += 1
```

`register_iterate` and `add` lose their matching edge and node bookkeeping. With that, `iterate_table[id, n]` is the only place an iterate relation is stored, as CLAUDE.md already says.

**Behaviour changes, all of them:**

1. `import tanglepack` no longer imports networkx. The reviewer and the test engineer both ran the check: `False` on the branch, `True` on base. networkx stays in `pyproject.toml`, and topology imports it lazily.
2. `henon_new_intersection_structures.py` no longer prints `Graph: 8 nodes, 17 edges` or opens the graphviz figure. `henon_period_3_new_intersection_structures.py` no longer prints `Graph: 207 nodes, 582 edges`.
3. `from tanglepack.topology import *` exports 48 names instead of 104.
4. If `henon_blast_period_3.py` ever reached `:123`, it would print the bool from `stop_reason`. It never gets there.

Nothing on a computed path changed. I read every src hunk: they are deletions, docstring rewordings, and the dropped `target_id` local above.

## What was verified, and how

- **Ran myself:**
  - the full suite at `517aa01` in the branch worktree: 1215/47/5/1, C3 only
  - the probe md5 on both output sets: equal
  - the `__all__` and import checks from `/tmp`
  - the git status of both trees
- **Read myself:** the whole diff, all 27 paths. Two test rewrites keep their strength:
  - `test_refinement.py` now calls `_refine_layer([(p0, p1)], …)` with `area_cutoff = 0.0` and still asserts the strict `lo < cdist < hi`. The fixture `grown_unstable` is function-scoped (`tests/conftest.py:196`), so mutating the machine does not leak into other tests.
  - In `test_iterate_inference.py`, each graph-edge read became `iterate_table[·, 1] == y`.
- **Taken from the test engineer's evidence, not repeated by me:**
  - the suite at each of the 8 commits
  - the collected-id diff: 8 ids gone (5 graph tests, `from_branch` × 3) and 1 renamed
  - the 7 break-it mutations, all caught
  - the script output diffs
  - `py_compile` of the long scripts

## What did not go to plan

The coder listed seven deviations. All of them stay inside the file list, and the reviewer accepted them:

- An unused `TYPE_CHECKING` import was dropped in C3.
- The `Iterable` import went in C4 rather than C5, because its last user was the graph code.
- `remove_iterate` no longer reads `target_id`.
- Extra docstrings were reworded so that none mentions a registry graph.
- `henon_bridge_classification.py` changes on 2 lines, not 1: a comment moved.
- The `test_blast_budget` docstring sentence was written out in full.
- Trailing blank lines were stripped.

Two more things:

- The test engineer rebuilt the probe harness. The B6b probe can no longer run on these trees because `WalkSearch.truncated` is gone. The replacement dumps more than the old one (every word, both matrices, statuses), so it is a stronger check, not a weaker one.
- The commit trailers say `Claude Fable 5.1`, as in B7b and B6b. This is cosmetic.

## What the deletion removed that you might want back

| Removed | Why you might miss it | Where it lives |
|---|---|---|
| The registry's networkx graph and `graphviz.visualize_intersection_graph` (layouts, label and colour modes, iterate edges) | A quick picture of crossings and iterate links when debugging a tangle | `git show 65edec4:src/tanglepack/numerics/graphviz.py`; deleted in `7eaf64d` |
| `IntersectionRegistry.graph()` and `_consecutive_on_branch`, plus 5 tests asserting that adjacency is per branch (a global cdist sort joins crossings that are not neighbours on a period > 1 orbit) | The per-branch rule now lives only in `Trellis`'s branch orderings | `git show 65edec4:src/tanglepack/numerics/IntersectionRegistry.py` and `:tests/unit/numerics/test_map_step_and_graph.py` |
| `from_branch`, whose 3-way test was the only registry-query test of D4 (the anchor lies on every branch) | D4 itself is still enforced: by the trellis-level law at `tests/helpers/laws.py:905` (run by `test_law_anchors.py`) and by the registry dedup test `at_point` (`test_intersection_registry_fixes.py`). What is gone is only a branch query on the registry. | deleted in `3a9db04`; the scripts' replacement is `filter(lambda ix: ix.is_anchor or ix.manifold_b_key[3] == b)` |
| `_get_refined_point`, the scalar reference version of the refinement kernel | A readable one-pair version of `_refine_layer`. The analyst found the two bit-identical on k10. | deleted in `88d5a75` |
| `as_forward_array` / `as_backward_array` | Dense iterate matrices for export | deleted in `3a9db04` |

## Open points for the author (seen, out of scope)

1. **Stale Sphinx stub.** No earlier report mentioned this. `docs/source/tanglepack.numerics.rst:183-189` still has `.. automodule:: tanglepack.numerics.graphviz`, so a Sphinx build will warn or fail on the missing module. It is a tracked `docs/` file, so it falls under B12.
2. **CLAUDE.md `:133`, `:136`, `:154`** are stale (B12).
3. **`henon_blast_period_3.py` is still dead** at `:95` (pinned pip 12 is no longer a candidate). Revive it or delete it? That is your call.
4. **The analyst's notes:**
   - The period-3 demo's anchor label is wrong.
   - `henon_new_intersection_structures.py` filters s-cdist in `[2, 6]` but prints "in [2, 8]". This predates the branch, and the rewrite kept it.
   - `IterateTable.forward_depth` / `forward_chain` read direct entries only.
5. **`test_map_step_and_graph.py` keeps "graph" in its name** although its graph tests are gone. Renaming it would move test ids, so it belongs to a later slice.
6. **56 F401 findings in `topology/__init__.py`** (Q8): B12 should trim the import block or add `__all__`-aware noqa.

## The three points I am least sure were done the way you meant

1. **Q8's "keep every import line".** `__all__` and the import block now disagree by 56 names, and ruff flags every one of them. That is exactly the plan, but if you expected the trim to make the namespace smaller, it did not: only `import *` shrank.
2. **The two demo scripts were rewritten, not retired (Q1).** With the graph and array sections gone, `henon_new_intersection_structures.py` mostly shows `filter`, which plain list comprehensions already do. If those scripts existed to show the registry API, they now show less of it.
3. **Losing the per-branch adjacency tests.** They tested a view nobody reads, so deleting them is consistent with the request. But they were the only tests that said in words "a global cdist sort is meaningless for period > 1". I believe, without having traced it, that the `Trellis` branch-order tests cover the same rule.

Knowledge files updated to describe the post-B10 code: `.refactor/knowledge/numerics.md`, `loom.md`, `testing.md`, `topology.md`. The B10 entry (Q1–Q8 with defaults, consequences, open items) is appended to `.refactor/decisions.md`.
