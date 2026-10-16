# Debrief: B9b `one-element-at-end` (auditor, 2026-10-15)

> **Tree state.** The main tree at `5f0db5f` has uncommitted documentation changes (`CLAUDE.md`, `docs/`, `.refactor/`, `artifacts/`) and untracked sweep scripts and figures. The branch was cut from `HEAD` and touches none of them. Its `src` and `tests` are clean against `HEAD`. Both worktrees are clean.

## In short

**The answer to your question is now yes.** Before this slice, four functions answered "which partition element owns this bridge end on this row". Now one does: `PartitionFamily.owner_at_end`.
- `Trellis.element_for` is deleted.
- `StablePartition._elements_at_bridge` is deleted, together with the `elements_at_bridge` field that nothing read.
- `BridgeClass._element_at_end` is deleted. The class rule now calls `owner_at_end` itself.

No class, word, matrix, count or `is_reliable` value changed on any of the 63 corpus items. I re-ran both probes on the branch myself. The only behaviour change is on error paths that nothing reaches.

**What to merge:**
- Branch `refactor/one-element-at-end`, head `94fca6b`, 3 commits on `5f0db5f`, a fast-forward.
- Suite: **1207 passed / 47 skipped / 5 deselected / 1 xfailed (C3)**, no XPASS (my run, 209 s).
- Size: **net −168 lines**. Source is +33/−182 (−149); tests are +43/−62 (−19).
- Nothing blocks the merge.

**What remains of the plan:** B12 `one-design-note` is the last slice. It is docs plus tidying: it moves the decision ids out of `src/` into one design note, rewrites CLAUDE.md, and replaces or deletes the tautological laws. "Open points" below lists what it now inherits. B11 `manifold-arrays` stays deferred (Q7: later, its own plan).

## Did the branch do what "done" asked?

| "Done" item | Status | How I know |
|---|---|---|
| One implementation on every computed path (`owner_at_end`) | **yes** | `ends` calls `family.owner_at_end` at `BridgeClass.py:598-599`. The other callers (SymbolicDynamics, DualWalk, plotting) already used it. The greps for `element_for\b`, `_elements_at_bridge`, `elements_at_bridge`, `_element_at_end`, `partition_element_for` and `def test_element_for` over `src` and `tests` print nothing (I ran them). |
| The three others gone, plus the never-read field | **yes** | Read in the diff: the field and both writers (`StablePartition.py:536`, `PartitionFamily.py:898`) and the import at `PartitionFamily.py:140` are gone. `StablePartitionResult` now has four fields. |
| The class rule reads the survivor through a function-level import | **yes** | `BridgeClass.py:582`, with a one-line comment giving the reason: `PartitionFamily.py:136` imports `_image_chain, _index_partitions` from `BridgeClass`, so a module-level import would be circular. I checked from `/tmp` with `python -I`: importing BridgeClass then PartitionFamily works, the reverse order works, and `tanglepack.loom` imports. `BridgeClass` no longer imports anything from `StablePartition`. |
| Tests rewritten through `owner_of_intersection`, both sides, no loss of strength | **yes** | Read below. |
| Stale docstrings in the touched files reworded | **yes** | The `StablePartition` module notes and the `partition_stable_manifold` docstring were reworded. In `TopologyResults`, "two lookup tables" became one table and the cite of `partition_element_for` was removed. Trellis's cite went with the method. In `bridge_classes`, `manifold_b_key` became `stable_key_at` (C7), and the Raises and Note sections were updated. A grep for the five stale phrases in the six files finds nothing. |
| No behaviour change | **yes** | **Agreement probe** (my rerun, `b9b/aud/probe/`): 7,090/7,090 row ends, 14,180/14,180 (end, side) pairs and 1,776/1,776 iterated ends agree, with 0 errors and 0 bad readings. `diff -r` of all 63 readings files against the test engineer's base run is empty. **63-item outcome probe** (my rerun, `b9b/aud/probe63/`): `cmp` against base exits 0, with 0 errors. |

Every file the diff touches is on the plan's list: `Trellis.py`, `StablePartition.py`, `PartitionFamily.py`, `BridgeClass.py`, `TopologyResults.py` and `tests/unit/topology/test_partition_elements.py`. The diff touches nothing under `laws.py`, `cases.py` or `tests/invariants/`. `ruff --select F` on the six files gives the same single pre-existing finding on base and on the branch: the unused `Intersection` import at `StablePartition.py:238`.

### The most important hunk: the class rule's `ends`

Before (`5f0db5f`, `BridgeClass.py:583, 595-599` and a private reader at `:763-791`):

```python
by_branch_side = _index_partitions(partitions)
...
def ends(bridge_id):
    return (
        _element_at_end(trellis, by_branch_side, bridge_id, "first"),
        _element_at_end(trellis, by_branch_side, bridge_id, "second"),
    )
# _element_at_end: row_of_end -> stable_key_at(unstable_key_of) -> by_branch_side[(branch, row)]
#                  -> element_of_intersection, with its own ValueErrors and an AssertionError
```

After (`94fca6b`, `BridgeClass.py:581-600`; the private reader is gone):

```python
# PartitionFamily imports this module (``_image_chain``), so not at module level.
from .PartitionFamily import HomotopyPartition

family = HomotopyPartition.from_results(partitions, trellis=trellis)
...
def ends(bridge_id):
    return (
        family.owner_at_end(trellis, bridge_id, "first"),
        family.owner_at_end(trellis, bridge_id, "second"),
    )
```

`owner_at_end` performs the same steps in the same order: `row_of_end`, then `stable_key_at`, then the `(branch, side)` result, then `element_of_intersection`. Its body is unchanged. The diff of `PartitionFamily.py` has only two hunks, the deleted import and the deleted writer.

`from_results` is `cls(results, trellis=)`, and its constructor calls the same `_index_partitions`, so duplicate partitions are still a `ValueError`. `test_bridge_class.py:332` still passes.

### What got shorter

| File | Base | Branch |
|---|---|---|
| `Trellis.py` | 1626 | 1570 |
| `StablePartition.py` | 1091 | 1048 |
| `TopologyResults.py` | 922 | 906 |
| `BridgeClass.py` | 791 | 761 |
| `PartitionFamily.py` | 1077 | 1073 |
| `test_partition_elements.py` | 168 lines, 5 tests | 149 lines, 4 tests |

The answer to "who owns this end" used to need three readers kept in agreement plus a stored table. Now there is one function, and the class rule calls it by name.

### Tests: is the strength kept?

I read all four tests in the new file:
- `test_every_workbench_bridge_touching_the_branch_resolves` still enumerates every workbench bridge, checks both sides and asserts `checked > 0`.
- The other two "resolves" tests check both ends and both sides. The branch-wrapping test now also asserts that the two ends lie on different branches.

All three compare whole `ElementRef`s (branch, side, id) where the old tests compared bare element ids, so they are **stronger** on that point. They go through `owner_of_intersection` on the `stable_key_at` branch, never through `owner_at_end`, as the request required.

Two old checks went with the deleted field, by design:
- that the table never invents a key;
- that the dropped bridge is absent from it.

The "third" test now goes through `owner_at_end` → `endpoint_index`.

Teeth:
- Test engineer: planting "owner_of_intersection returns the neighbour element" fails 3 of the 4 tests, and planting "endpoint_index accepts third" fails the fourth.
- Reviewer: re-planted the first break independently and got the same result.

## Not done, and why

- **`plot_resonance_zones` smoke test.** Not added (Q4). It goes to B12. The draft is at scratchpad `b9b/smoke_draft.py` on `k28_partitioned`; `k10_partitioned` has no zone.
- **`ResonanceZone.py:295,353`.** Untouched (Q5). The analyst found both still accurate: `TangleSession.classify_bridge` still orders zones by area, and `add_resonance_zones` exists. I did not re-check this beyond reading the lines.
- **`CLAUDE.md`'s `Trellis` row still lists `element_for`.** This is `:137` at `5f0db5f` and `:139` in your uncommitted working copy. It goes to B12.
- **The unused `Intersection` import at `StablePartition.py:238`.** It is pre-existing, ruff reports it on base too, and the plan put it out of scope.

## Every behaviour change

All four are on error paths. No corpus item, test or law reaches any of them.

1. **Unowned crossing in the class rule** (Q2, approved). It used to raise `AssertionError`; it now raises `ValueError` from `owner_of_intersection`, which is what every other caller already got. Nothing pinned the old type.
2. **Missing partition.** Still a `ValueError`, and the message still names the branch (`branch_key[1:]`), so the Raises docstring holds. The text is now `owner_of_intersection`'s ("crossing N is on no 'left'-partitioned branch (o, b) of this homotopy partition family"). It drops the old hint "partition that branch (or pass its trellis's stable_partitions)". `test_bridge_class.py:311` pins only the type and the recovery, and it passes.
3. **The `branch_key is None` `ValueError` is gone.** That path was dead: every registered crossing is unstable × stable.
4. **The build-time assertion in `_elements_at_bridge` is gone.** It checked that an on-branch bridge end has an owner. Three things still cover it:
   - the laws `check_partition_unique_owner` and `check_iterated_unique_owner`;
   - the partition's own build assertions;
   - `owner_of_intersection`'s `ValueError` at read time.

## What was removed that someone might want back

| What | Where to find it | Why you might want it |
|---|---|---|
| `Trellis.element_for(bridge_id, endpoint, side) -> Optional[int]` | `git show 5f0db5f:src/tanglepack/topology/Trellis.py`, lines 1379-1430 | It was the only reader that took **any side** and returned `None` instead of raising when the branch was not partitioned. The replacement is `HomotopyPartition.from_results(trellis.stable_partitions, trellis=trellis).owner_of_intersection(id, side, trellis.stable_key_at(id, trellis.unstable_key_of(bid)))`. It raises where the old method returned `None`. |
| `StablePartition._elements_at_bridge` and the field `StablePartitionResult.elements_at_bridge` | `git show 5f0db5f:src/tanglepack/topology/StablePartition.py`, lines 1059-1091; `TopologyResults.py`, line 442 | A per-result index from bridge to its `(first, second)` elements, with `None` for an end on another branch. It was written but never read in `src`. |
| `BridgeClass._element_at_end` | `git show 5f0db5f:src/tanglepack/topology/BridgeClass.py`, lines 763-791 | — |
| Test `test_elements_at_bridge_entries_agree_with_the_ownership_table` | Deleted in commit `cfd89e6`; recover with `git show 5f0db5f:tests/unit/topology/test_partition_elements.py` | — |
| The three renamed tests | Same file, under their old names | They were `test_element_for_*`. |

## What did not go to plan

Nothing material went off plan. Small deviations:
- **Commit trailers** read `Co-Authored-By: Claude Opus 5.5`. The spawn asked for `Claude Fable 5.1`. It is cosmetic: B9a's commits say Fable 5.1, these say Opus 5.5.
- **Docstring wording** goes a little beyond the analyst's candidate:
  - the `bridge_classes` Note;
  - a singular "it" in the `StablePartitionResult` Note;
  - the `:540-544` sentence, which the plan asked for.

  As a result, source came out −149 rather than −152. That is inside the rubric's window.
- **Commit 1** briefly checked the field across every result before commit 2 deleted it. This is harmless.

## Open points for you

1. **The rewritten test file does not catch two plausible breaks.** These are the test engineer's extra breaks:
   - (c) `stable_key_at` ignoring the anchor rule C7;
   - (e) `owner_at_end` reading the opposite row.

   Both pass `test_partition_elements.py`. The suite still catches both: (c) gives 16 failures and 19 errors, and (e) gives 57 failures, mostly in the law tier through the oracle `laws._element_ends`. The old file never tested the row either, so this is no regression. But now that `owner_at_end` is the only reader, its row has no dedicated unit test. Adding one is a B12-or-later question, not a defect.
2. **Commit trailers** say Claude Opus 5.5 (see above).
3. **Four stale committed higher-period PNGs** (`p3`, `nested`, `heteroclinic`, `inversion`) predate this branch. Nothing compared against them.
4. **`scripts/henon_blast_period_3.py:95`** is still dead (`T1.set_strong_pip(12)`). It is reported, not fixed, per your rule.
5. **What B12 now has to cover:**
   - CLAUDE.md rows for B10 (`graphviz`, `graph(bridges=)`, `terminated_early`), for B9a (the `TangleSession` row, a `cartoon.py` row) and for B9b (the `Trellis` row's `element_for`, and the `StablePartitionResult` wording);
   - the `plot_resonance_zones` smoke test;
   - the 56 F401 findings in `topology/__init__.py` (ruff, re-counted on the branch);
   - the Sphinx stub `docs/source/tanglepack.numerics.rst:183-189`, which still runs `automodule` on the deleted `tanglepack.numerics.graphviz`.

## How I verified this

**Ran myself:**
- the full suite on `94fca6b`;
- the agreement probe and the 63-item outcome probe on the branch, compared against the test engineer's base outputs (which equal the analyst's);
- the greps;
- the import-cycle checks;
- `ruff --select F` on base and on the branch;
- `git diff` / `--stat` / `--numstat` against `5f0db5f`, all of which I read.

**Taken from the test engineer's logs:**
- the per-commit suites (1208 / 1207 / 1207);
- the collected-id diff;
- the six planted breaks.

**Believed, not re-checked:** the `ResonanceZone.py` citations are accurate (the analyst checked them).

My scratch outputs are in `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b9b/aud/`.

I updated the knowledge files `.refactor/knowledge/topology.md` ("B9b outcome"), `loom.md` and `testing.md`, and appended this run's decisions Q1-Q5 to `.refactor/decisions.md`.
