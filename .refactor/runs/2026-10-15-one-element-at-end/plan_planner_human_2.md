# Plan: B9b `one-element-at-end`

## In short

Four functions answer the same question: which partition element owns a bridge end, on that end's row. This slice keeps one of them, `PartitionFamily.owner_at_end`, and deletes the other three:

- `Trellis.element_for` has no caller in `src`. It is deleted.
- `StablePartition._elements_at_bridge` fills a field that nothing in `src` reads. It is deleted, along with the `StablePartitionResult.elements_at_bridge` field, the two places that write it, and the import of it.
- `BridgeClass._element_at_end` is deleted. The class rule gets its owner from `owner_at_end` itself: inside `bridge_classes`, a function-level import builds a `HomotopyPartition` from the partitions it was given.
  - The import has to sit inside the function because `PartitionFamily` already imports from `BridgeClass`, so a module-level import would be circular.
  - This is the analyst's option B.

Nothing changes in the output:

- On `HEAD`, the four functions agree on every bridge end of all 63 corpus items: 7,090 ends on their row, 14,180 (end, side) pairs, and 1,776 iterated ends, with 0 disagreements.
- The analyst built the whole slice on a scratch copy (`b9b/cand_b`):
  - Each item's readings hash equals base's (63/63).
  - The 63-item outcome probe is byte-identical.
  - The suite gives 1207 passed / 47 skipped / 5 deselected / 1 xfail (C3).
- Source shrinks by about 152 lines net (+23 / −175).

**What you decide.** Five questions are at the end, each with a default. "Approved" accepts all five defaults:

1. Option B as described above.
2. An unowned crossing in the class rule raises `ValueError` instead of `AssertionError`.
3. Rename the three tests whose names cite `element_for`.
4. No `plot_resonance_zones` smoke test in this slice.
5. Leave `ResonanceZone.py` alone: its two citations are still accurate.

## Recommendation

Change it. This is a pure deletion, and one implementation remains on every computed path.

## Scope

- **Slug:** `one-element-at-end`
- **Branch:** `refactor/one-element-at-end`
- **Worktree:** `../tangle-pack-one-element-at-end`, branched off `HEAD` = `5f0db5f`
- **Base worktree:** `../tangle-pack-one-element-at-end-base`, detached at `5f0db5f`

**Files the change may touch** (under the defaults):

1. `src/tanglepack/topology/Trellis.py`
2. `src/tanglepack/topology/StablePartition.py`
3. `src/tanglepack/topology/PartitionFamily.py`
4. `src/tanglepack/topology/BridgeClass.py`
5. `src/tanglepack/topology/TopologyResults.py`
6. `tests/unit/topology/test_partition_elements.py`

**Under a non-default answer:**

- **Q1 = A:** the same six files. `StablePartition.py` gains `end_on_row` and `owner_in` next to `row_of_end` (src −133 instead of −152).
- **Q4 = yes:** `tests/plotting/test_plot_smoke.py` joins the list, with one new test.
- **Q5 = fix:** `src/tanglepack/loom/ResonanceZone.py` joins the list.

**Never touched:**

- `tests/helpers/laws.py`, `tests/cases.py`, `tests/invariants/`
- `CLAUDE.md`, `docs/`, `.refactor/`, `artifacts/`, `figures/`, and the two untracked sweep scripts

The size is about 175 lines removed and 23 added in source, and +43 / −62 in the test file. That is targeted.

## Steps: three commits, each green

There are three commits, one per deletion. Each one stands on its own and can be reverted alone. The order is chosen so that every intermediate tree passes the suite:

1. The test rewrite lands with the first deletion that needs it.
2. The field test goes with the field.
3. The class rule moves last.

The analyst's `b9b/cand_b_src.diff` and `b9b/cand_b_tests.diff` are the reference for every hunk. The coder applies them per commit as split below, plus the one docstring clause that `cand_b` missed (commit 3).

All line numbers refer to `HEAD`.

### Commit 1: delete `Trellis.element_for`

`src/tanglepack/topology/Trellis.py`:
- Delete `:1379-1430`. This is the `# ── partition elements` header and the whole `element_for` method, including its docstring at `:1413`, which cites the deleted session method.
- Drop the now-unused imports `Endpoint` (`:77`), `Side` (`:80`) and `endpoint_index` (`:83`) from the `TopologyResults` import block. Keep `Optional` if anything else still uses it; ruff decides.

`tests/unit/topology/test_partition_elements.py`: rewrite the four tests that call `element_for`, using `cand_b_tests.diff`.
- Add two module helpers:
  ```python
  def _owner(trellis, family, bid, end, side):
      return family.owner_of_intersection(
          end, side, trellis.stable_key_at(end, trellis.unstable_key_of(bid)))

  def _expected(trellis, end, side):   # the independent derivation, unchanged in spirit
      result = <the result of the end's own manifold_b_key branch and side>
      return result.ref(result.element_of_intersection[end])
  ```
- Build `family = HomotopyPartition.from_results(trellis.stable_partitions, trellis=trellis)`, which reads exactly the results `element_for` read.
- Compare whole `ElementRef`s, not bare ids.
- Use `owner_of_intersection` on both sides. Never route through `owner_at_end`, which would weaken a both-side test to the row only.
- `:42 test_every_workbench_bridge_touching_the_branch_resolves`:
  - Keeps its name, the workbench enumeration, both sides and `checked > 0`.
  - The `elements_at_bridge` clause at `:74` **stays in this commit** and goes in commit 2.
- `:96`: becomes `test_both_ends_of_a_branch_wrapping_bridge_resolve`.
  - Asserts both sides and both ends.
  - Asserts `first.branch_key != second.branch_key`.
- `:119`: becomes `test_a_bridge_missing_from_the_snapshot_resolves`.
  - Keeps the drop-and-repartition mechanism.
  - Rebuilds `family` after `partition_stable_manifold()`.
  - The `:155` clause stays until commit 2.
- `:161`: becomes `test_owner_at_end_rejects_an_unknown_endpoint_name`.
  - Asserts `pytest.raises(ValueError)` on `family.owner_at_end(trellis, bridge.id, "third")`.
- Module docstring `:9-11`: stop naming `element_for`. The `elements_at_bridge` mention goes in commit 2.

### Commit 2: delete `_elements_at_bridge` and the `elements_at_bridge` field

`src/tanglepack/topology/StablePartition.py`:
- Delete the writer at `:536-538`.
- Delete the helper `_elements_at_bridge` at `:1059-1091` (end of file).
- Reword the module notes `:23-32`: one lookup table (`element_of_intersection`). Stop citing `elements_at_bridge` and `Trellis.element_for`.
- Reword the `partition_stable_manifold` docstring at `:499-504` the same way.

`src/tanglepack/topology/PartitionFamily.py`:
- Delete the import at `:140`.
- Delete the iterated writer at `:898-900`.

`src/tanglepack/topology/TopologyResults.py`:
- Delete the field `:442-444` and its doc `:415-427`, including the `partition_element_for` cite.
- Reword `:400-403` ("the two lookup tables …") and `:429` ("Both tables") to the one table.
- The field is the fifth and last field, so `test_exact_ownership.py:72`, which builds the result positionally with 4 args, is unaffected.

`tests/unit/topology/test_partition_elements.py`:
- Delete `test_elements_at_bridge_entries_agree_with_the_ownership_table` (`:78`). It exists only for the field.
- Delete the clauses at `:74` and `:155`.
- Remove the last `elements_at_bridge` mention from the module docstring and `:47-51`.

### Commit 3: the class rule reads `owner_at_end`; delete `_element_at_end`

`src/tanglepack/topology/BridgeClass.py`:
- At `:583`, replace `by_branch_side = _index_partitions(partitions)` with:
  ```python
  # PartitionFamily imports this module (``_image_chain``), so not at module level.
  from .PartitionFamily import HomotopyPartition

  family = HomotopyPartition.from_results(partitions, trellis=trellis)
  ```
- In `ends` (`:595-599`), call `family.owner_at_end(trellis, bridge_id, "first")` and `family.owner_at_end(trellis, bridge_id, "second")`.
- Delete `_element_at_end` (`:763-791`).
- Keep `_index_partitions`. `from_results` and `PartitionFamily` use it, so the duplicate-partition `ValueError` at `test_bridge_class.py:332-337` is unchanged.
- Drop the unused imports `row_of_end` (`:97`), `Endpoint` and `endpoint_index` (`:100-101`).

The `bridge_classes` docstring:
- `:540-544`: "the endpoint's `manifold_b_key` picks the stable branch" has been stale since C7. Change it to "`Trellis.stable_key_at` picks the stable branch (C7)", and say that each end is read by `HomotopyPartition.owner_at_end`.
- Raises `:569-574`: drop "if a crossing carries no stable manifold key". Add "if a covering result does not own a crossing" to the `ValueError`. **`cand_b` missed this clause.**
- Note `:575-581`: use `cand_b`'s wording. The unowned crossing becomes `owner_at_end`'s `ValueError`, and the image pair held by two chains stays an `AssertionError`.

`PartitionFamily.owner_at_end` itself is **not edited**. It keeps calling `self.owner_of_intersection`, which the `_OwnerFamily` fake at `test_dual_walk.py:877` relies on.

## Behaviour

None on any computed path.

The one change is an error path that no corpus item, test or law reaches. In `bridge_classes`, a crossing on a partitioned branch that the covering result does not own now raises `ValueError` (from `owner_of_intersection`) instead of `AssertionError` (Q2). The invariant is still enforced in three places:
- the partition's own assertions when it is built;
- the laws `check_partition_unique_owner` and `check_iterated_unique_owner`;
- `owner_of_intersection` at read time.

The dead `branch_key is None` path disappears. It cannot be reached, because every registered crossing is unstable × stable.

## Not doing

- Any change to `row_of_end`, `stable_key_at`, `unstable_key_of`, `owner_of_intersection`, `interval_at`, `anchor_outward_key`, loop folding or inertness.
- An "any side" variant of `owner_at_end`. B9a's Q6 said no.
- Changing `bridge_classes`'s signature to take a `HomotopyPartition` directly. That would touch `Topology.py:140`, `helpers/direct.py:102` and many unit tests.
- Moving `_image_chain` out of `BridgeClass` to break the cycle (option C). Four test files and four modules import it from there.
- The `plot_resonance_zones` smoke test (Q4). It goes to B12; the draft is at `b9b/smoke_draft.py`, on `k28_partitioned`.
- `CLAUDE.md:139`, which still lists `element_for` (B12).
- The pre-existing unused `Intersection` import at `StablePartition.py:238`.
- `scripts/henon_blast_period_3.py:95`, which is dead. It is reported, not fixed.

## Rubric

The reviewer checks each item by running or reading. Paths under `b9b/` are in the scratchpad `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b9b/`.

1. **Suite on the branch.** Run `MPLBACKEND=Agg PYTHONPATH=src ../tangle-pack/env/bin/python -m pytest -q -rxX -p no:cacheprovider` in the worktree.
   - Expected: **1207 passed / 47 skipped / 5 deselected / 1 xfailed**.
   - The one xfail is C3. There is no XPASS.
   - Each of the three commits also passes on its own: check out each commit and run the suite; all three are green.
2. **Every pinning test passes on base and on the branch.** These tests pin behaviour the slice must keep:
   - `test_bridge_class.py:311-337` (the missing-partition and duplicate-partition `ValueError`s)
   - `test_partition_family.py:109-110` (the `owner_of_intersection` `ValueError`)
   - `test_dual_walk.py` (the `_OwnerFamily` fake)

   They pass unmodified in both worktrees.
3. **Collected ids.** Run `pytest --collect-only -q` in base and in branch. The sorted lists differ by exactly the following, 1256 → 1255:
   - deleted: `test_elements_at_bridge_entries_agree_with_the_ownership_table`
   - renamed: `test_element_for_resolves_both_ends_of_a_branch_wrapping_bridge` → `test_both_ends_of_a_branch_wrapping_bridge_resolve`
   - renamed: `test_element_for_resolves_a_bridge_missing_from_the_snapshot` → `test_a_bridge_missing_from_the_snapshot_resolves`
   - renamed: `test_element_for_rejects_an_unknown_endpoint_name` → `test_owner_at_end_rejects_an_unknown_endpoint_name`
4. **The diff stays inside the file list.** `git diff --stat 5f0db5f..refactor/one-element-at-end` lists exactly the six files above. In particular, `laws.py`, `cases.py`, `tests/invariants/`, `CLAUDE.md`, `docs/` and `.refactor/` are absent. The worktree started clean, so the author's uncommitted changes in the main tree are not involved.
5. **The deleted names are gone.**
   - `grep -rn "element_for\b\|_elements_at_bridge\|elements_at_bridge\|_element_at_end" src tests` prints nothing.
   - `grep -rn "partition_element_for" src tests` prints nothing.
   - `grep -rn "def test_element_for" tests` prints nothing. This check is needed because `\b` does not stop at `_`, so the first grep alone would miss the old test names.
6. **`Trellis.element_for` no longer exists.** `hasattr(Trellis, "element_for")` is False. `"elements_at_bridge" not in {f.name for f in dataclasses.fields(StablePartitionResult)}`. `BridgeClass.py` defines no `_element_at_end`.
7. **The class rule calls the survivor.** In `bridge_classes`, `ends` calls `family.owner_at_end` twice. `BridgeClass.py` imports nothing from `StablePartition`. `PartitionFamily.owner_at_end`'s body is identical to base's (`git diff` shows no hunk inside `:343-367`).
8. **Agreement probe on the branch.** Run `b9b/probe/run.sh <worktree> out_branch`.
   - Counts: **7,090 / 7,090** row ends (`owner_at_end`, homotopy, = the law oracle `laws._element_ends`), **14,180 / 14,180** (end, side) pairs (`owner_of_intersection`), and **1,776 / 1,776** iterated ends.
   - 0 exceptions.
   - All 63 per-item `readings_sha256` equal those in `b9b/probe/all_out_base.jsonl`.
9. **63-item outcome probe.** Run `b9b/probe63/run.sh` on the base worktree and on the branch. The sorted outputs are byte-identical (`cmp` exits 0), with 0 errors. This covers every class, word, matrix, count, walk status and `is_reliable` value.
10. **No import cycle.** Run from `/tmp` with a fresh interpreter:
    - `-I -c "import sys; sys.path.insert(0, '<worktree>/src'); import tanglepack.topology.BridgeClass, tanglepack.topology.PartitionFamily"`
    - the same with the two imports in the opposite order
    - `import tanglepack.loom`

    All three succeed.
11. **The rewritten tests have teeth.** Plant each break in a throwaway copy, never in a commit:
    - Make `owner_of_intersection` return the neighbouring element (`result.ref(element_id + 1)` where it exists). Then `test_every_workbench_bridge_touching_the_branch_resolves`, `test_both_ends_of_a_branch_wrapping_bridge_resolve` and `test_a_bridge_missing_from_the_snapshot_resolves` all fail. Expected: 3 failed, 1 passed in the file.
    - Make `endpoint_index` accept `"third"` (return 1). Then `test_owner_at_end_rejects_an_unknown_endpoint_name` fails. Expected: 1 failed, 3 passed.
12. **Line counts.**
    - `git diff --numstat` over `src/` gives about +23 / −175. The net must be between −140 and −165.
    - `test_partition_elements.py` ends at about 149 lines, with 4 tests (5 at base).
13. **No stale citations.** No docstring or comment in the six files names `element_for`, `elements_at_bridge`, `_element_at_end`, `partition_element_for`, or "two lookup tables" (grep).
    - The `bridge_classes` Raises clause no longer mentions "carries no stable manifold key".
    - Its Note no longer calls an unowned crossing an `AssertionError`.
    - `:540-544` names `stable_key_at`, not `manifold_b_key`.
14. **No behaviour change beyond the Behaviour section.** Items 8 and 9 cover this, plus a read of the three commits: the only semantic change is the `bridge_classes` error type on an unreachable path.

## Questions

Each question has a default. "Approved" accepts all five.

- **Q1. Where does the class rule get its owner?**
  - (B) A function-level import of `HomotopyPartition` in `bridge_classes`, which calls `owner_at_end`. This deletes `_element_at_end`, adds no function, and removes 152 src lines net. Function-level imports already exist in `StablePartition.py:635,660`.
  - (A) Two shared helpers, `end_on_row` and `owner_in`, beside `row_of_end`. This needs no local import but adds two new functions, removes 133 lines net, and the class rule still does not call `owner_at_end`.
  - *Default: B.*
- **Q2. Unowned crossing in the class rule: `AssertionError` becomes `ValueError`?** Nothing pins the `AssertionError`. The `ValueError` is what every other caller already gets, and it is pinned. *Default: yes.*
- **Q3. Rename the three surviving tests whose names cite `element_for`?** *Default: rename* (rubric 3 lists the names).
- **Q4. Add the `plot_resonance_zones` smoke test here?** *Default: no, leave it for B12.* It is unrelated to this slice and would widen the file list.
- **Q5. `ResonanceZone.py:295,353`?** *Default: no change.* Both are accurate at `HEAD`: `TangleSession.classify_bridge` still orders by area, and `add_resonance_zones` exists at `TangleSession.py:696`.
