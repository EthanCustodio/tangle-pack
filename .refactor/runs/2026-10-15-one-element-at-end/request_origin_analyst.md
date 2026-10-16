# Request: B9b — `one-element-at-end`

*From the author, 2026-10-15 ("continue"). Input for one `tangle_refactor` run. The fifteenth slice of `docs/alignment_plan_2026-10-09.md`: the second half of skeleton row B9 (~line 1102; the element-at-end note ~`:246-249`; B7's "Not doing" ~`:1038`), split off at the B9a plan assay (Q0, 2026-10-15; the B9a plan's "B9b skeleton" section and the B9a debrief's "B9b next" paragraph in `.refactor/runs/2026-10-15-session-three-concerns/` are the scoping this request is written from). B0 … B10 and B9a are merged (main `5f0db5f`; runs archived under `.refactor/runs/2026-10-1{2,3,4,5}-*/`).*

## The question

**Four functions answer "which partition element owns this bridge end on this row"; can three of them fold into `PartitionFamily.owner_at_end` without changing a single class, word or figure?** At `HEAD`:

| Implementation | Where | Reads | Side | Callers at `HEAD` |
|---|---|---|---|---|
| `PartitionFamily.owner_at_end` | `PartitionFamily.py:343` (calls `owner_of_intersection` `:369`) | its own family (homotopy or iterated), the trellis passed | the ROW (`row_of_end`, D6) on the branch `stable_key_at` names (C7) | `SymbolicDynamics.py:642-643`, `DualWalk.py:1366`, `plotting.py:2303`; `tests/unit/topology/test_dual_walk.py:877` binds it onto a fake |
| `Trellis.element_for` | `Trellis.py:1381` | `self.stable_partitions` (homotopy) | ANY side, passed in | tests only since B9a: `test_partition_elements.py:68,105,106,158,168`; docstring cites at `StablePartition.py:31`, `TopologyResults.py:424` |
| `StablePartition._elements_at_bridge` | `StablePartition.py:1059` | the result being built, `trellis.bridges` | each side | written into the never-read-in-src field `StablePartitionResult.elements_at_bridge` (`TopologyResults.py:415-424`, `:442`) at `StablePartition.py:536` and `PartitionFamily.py:898` (import `:140`); read by `test_partition_elements.py:74,84,155` |
| `BridgeClass._element_at_end` | `BridgeClass.py:763` | the homotopy `partitions` list | the row | `BridgeClass.py:597-598` (the class rule's `ends`; `_non_loop_parent` `:740` via `ends`) |

The B9a analyst's agreement probe (scratchpad `b9/probe/`, `all_base.jsonl`) found **0 disagreements** on the 63 corpus items: 7,090 of 7,090 bridge ends on their row (all four + the law oracle `laws._element_ends` + the since-deleted `partition_element_for`), 14,180 of 14,180 (end, side) pairs (`owner_of_intersection` = `element_for` = `elements_at_bridge`), 1,776 of 1,776 iterated ends. So the fold is the identity on the corpus; the slice is a deletion with one shared core.

**The import cycle.** `PartitionFamily.py:136` imports `_image_chain, _index_partitions` from `BridgeClass`, so `BridgeClass` cannot import `PartitionFamily` at module level. The shared core therefore sits beside `row_of_end` in `StablePartition.py` (`:545`; `BridgeClass.py:97` already imports from there) or behind a local import — the planner decides and says why. The error paths differ: `_element_at_end` raises `AssertionError` on an unowned crossing where `owner_of_intersection` raises `ValueError`; no corpus item reaches either path; the planner states which one survives.

## What "done" means

Exactly the B9b skeleton of the B9a plan:

1. **One implementation on every computed path.** `Trellis.element_for` deleted (its docstring `:1413` cites the deleted session method); `_elements_at_bridge` deleted with the `elements_at_bridge` field, its two writers and the `PartitionFamily.py:140` import (the field is written and never read in src; `TopologyResults.py:427` cites the deleted session method too); `BridgeClass._element_at_end` becomes one call into the shared core (or into `owner_at_end` through a homotopy family, if the planner finds that cleaner and cycle-free). `owner_at_end` keeps calling `self.owner_of_intersection` (the dual-walk fake at `test_dual_walk.py:877` depends on it).
2. **Tests keep their strength.** `tests/unit/topology/test_partition_elements.py` (168 lines, 5 tests at `HEAD`): the both-side checks (`test_every_workbench_bridge_touching_the_branch_resolves` `:42`, `test_element_for_resolves_both_ends_of_a_branch_wrapping_bridge` `:96`, `test_element_for_resolves_a_bridge_missing_from_the_snapshot` `:119`) are rewritten with `family.owner_of_intersection(id, side, trellis.stable_key_at(id, trellis.unstable_key_of(bridge_id)))` — never through `owner_at_end`, which would weaken them to the row only; `test_elements_at_bridge_entries_agree_with_the_ownership_table` `:78` and the `elements_at_bridge` clauses at `:74`, `:155` go with the field; `test_element_for_rejects_an_unknown_endpoint_name` `:161` becomes `owner_at_end(…, "third")` raising `ValueError` through `endpoint_index`. The module docstring `:10` stops naming deleted API.
3. **Stale docstrings in the touched files** (`StablePartition.py:27,31,503`, `TopologyResults.py:424-427`, `Trellis.py:1413`) reworded; `ResonanceZone.py:295,353` (cite `TangleSession` for what is now `loom.cartoon.cartoon_zones`) is a planner question — default: fix the two lines here, since B12 is docs-only (CLAUDE.md), or leave for B12; say which.
4. **No behaviour change.** Every class, word, matrix, count and `is_reliable` identical on the 63 items; the agreement probe re-pointed at the survivor gives the same 7,090 / 14,180 / 1,776 on the branch.

Files (re-derived): `src/tanglepack/topology/Trellis.py`, `StablePartition.py`, `PartitionFamily.py`, `BridgeClass.py`, `TopologyResults.py`, `tests/unit/topology/test_partition_elements.py`; `src/tanglepack/loom/ResonanceZone.py` only under the docstring default. Small: ~−130 src lines, ~200 test lines rewritten, no new file.

**Behaviour.** None (rubric: the rebuilt 63-item probe at scratchpad `b10/probe/`, rerun base vs branch, 0 differences; the agreement probe at `b9/probe/probe.py`, re-pointed from `partition_element_for` and `element_for` to the survivors, 0 disagreements on the branch).

**Not doing.** Any change to `row_of_end`, `stable_key_at`, `unstable_key_of`, `owner_of_intersection`, `interval_at`, the class rule's orientation (`anchor_outward_key`), loop folding or inertness; a new "any side" API on `owner_at_end` (Q6 of B9a: no); the `plot_resonance_zones` smoke test the B9a debrief asked for (B12, or a question for the author here — default: add it in this slice ONLY if the planner says the rubric can carry one new collected id, else B12); CLAUDE.md rows (B12); anything the law tier reads (`laws.py:75,1177-1197,1323,1868,1874,1979,2200,2226,2252` read `row_of_end`, `element_of_intersection` and its own `_element_ends` oracle — none of the four — and stay untouched).

**Rubric** (the planner re-derives):
1. Suite passes; `pytest -q -rxX` lists exactly one xfail (C3), no XPASS; `git diff --stat` touches only the listed files; collected ids differ from base only by the deleted/renamed `test_partition_elements` tests the plan names (and the one new smoke test, if chosen).
2. `grep -rn "element_for\b\|_elements_at_bridge\|elements_at_bridge\|_element_at_end" src tests` → nothing (or exactly the one-line forward the plan keeps); `grep -rn "partition_element_for" src tests` → nothing.
3. The agreement probe on the branch: `owner_at_end` (homotopy and iterated) equals the law oracle `_element_ends` on every bridge end × row of all 63 items, and `owner_of_intersection` on both sides equals base's `element_for`/`elements_at_bridge` readings recorded in `all_base.jsonl`.
4. Base-vs-branch 63-item probe equal; the law tier untouched (`laws.py`, `cases.py`, `tests/invariants/` not in the diff).
5. No import cycle: `python -c "import tanglepack.topology.BridgeClass, tanglepack.topology.PartitionFamily"` from `/tmp` in a fresh interpreter succeeds in both orders.
6. Teeth: each rewritten both-side test fails under a planted break of `owner_of_intersection` (e.g. return the neighbour element) and the `ValueError` test fails if `endpoint_index` accepts "third".

## Where to start

`.refactor/runs/2026-10-15-session-three-concerns/plan_planner_human_2.md` ("B9b skeleton") and `outcome_auditor_human_2.md` ("B9b next"); `.refactor/runs/2026-10-15-session-three-concerns/findings_analyst_planner.t1.md` §4 (the four implementations, the probe, the five things to know before folding); scratchpad `b9/probe/{probe.py,run.sh,all_base.jsonl}`; `src/tanglepack/topology/PartitionFamily.py:343-400`, `Trellis.py:1381-1420`, `StablePartition.py:545` (`row_of_end`), `:1059` (`_elements_at_bridge`), `:536`, `BridgeClass.py:590-600`, `:740`, `:763-791`, `TopologyResults.py:410-445`; `tests/unit/topology/test_partition_elements.py`; `tests/unit/topology/test_dual_walk.py:870-885`; `tests/helpers/laws.py:1868` (`_element_ends`).

## Binding context

- Test policy (`tests/conftest.py`, 2026-10-05): no case-specific pins; a test that exists only for a deleted path goes with it; a test that uses it incidentally is rewritten with the SAME assertion strength (the both-side checks must stay both-side); `KNOWN_ISSUES` is C3 alone; the law tier is untouched.
- `CLAUDE.md` invariants are physics and none is touched: rows are the signs' (`row_of_end`, D6); an anchor end reads `stable_key_at` (C7); exactly one element owns each crossing per side (by id at a boundary, by exact cdist elsewhere); `ElementRef` carries no family tag and the two families are never mixed; `numerics`/`topology` never import `loom`.
- Author's rules in force: one slice per run; topology code never calls the map; a dead script is reported, not fixed (`scripts/henon_blast_period_3.py` still dead at `:95`, untouched here).
- The analyst stage is NOT skipped: the analyst re-derives the four implementations and their callers at `HEAD` (= `5f0db5f`, after B9a), re-runs the agreement probe on `HEAD` (it last ran on `517aa01`), traces the import cycle and proposes where the shared core lives, reads what each of the five tests actually asserts, and lists every docstring that cites a function to be deleted. It does not re-plan the slice.
- The working tree may carry uncommitted documentation changes (`docs/`, `.refactor/`, `CLAUDE.md`); the change branches from `HEAD` and never touches them.
