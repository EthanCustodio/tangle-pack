# Test evidence: B9b `one-element-at-end` (test engineer, round 1, 2026-10-15)

## Setup

- **Base:** `5f0db5f`, detached worktree `/home/dezhu/Development/TanglePack/tangle-pack-one-element-at-end-base`.
- **Branch:** `refactor/one-element-at-end` @ `94fca6b` (commits `ee9fe8a`, `cfd89e6`, `94fca6b`), worktree `/home/dezhu/Development/TanglePack/tangle-pack-one-element-at-end`.
- **Scratch:** `TE=/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b9b/te/`. Every log named below is there.
- **Suite command,** run in each tree: `MPLBACKEND=Agg PYTHONPATH=src <main>/env/bin/python -m pytest -q -rxX -p no:cacheprovider`.
- **No test commit.** No rewritten test needed strengthening (see rubric 11), so I made no commit on the branch. Branch HEAD is still `94fca6b`.
- **Clean trees.** Both worktrees are clean at the end. `git status --short -- src tests` in the main tree is empty. The analyst's `b9b/probe*` copies were not edited.
- **Who ran what.** I ran every result in this report myself, unless it says otherwise.

## Rubric, item by item

| # | Item | Ran | Result |
|---|---|---|---|
| 1 | Suite at head and at each commit | yes | **pass**: head 1207 / 47 / 5 / 1 xfailed (C3); `ee9fe8a` 1208 (the field test still present); `cfd89e6` and `94fca6b` 1207; never an XPASS |
| 2 | Pinning tests, unmodified, on base and on branch | yes | **pass**: 59 passed in both worktrees; 84 passed for the three whole files in both |
| 3 | Collected ids, base vs branch | yes | **pass**: 1256 → 1255, one deletion and the three renames, nothing else |
| 4 | Diff is the six files | yes | **pass** |
| 5 | The three greps print nothing | yes | **pass**: all three exit 1 |
| 6 | Deleted names are gone (hasattr / fields) | yes | **pass** |
| 7 | `ends` calls `owner_at_end`; `owner_at_end` itself untouched | yes | **pass** |
| 8 | Agreement probe | yes, base and branch | **pass**: 7,090 / 14,180 / 1,776 agreements, 0 exceptions, 63/63 digests equal |
| 9 | 63-item outcome probe | yes, base and branch | **pass**: `cmp` exits 0, 0 errors |
| 10 | No import cycle | yes | **pass**: all three imports succeed |
| 11 | Teeth | yes, 2 planned breaks and 4 extra | **pass**: (a) 3F/1P, (b) 1F/3P; all four extra breaks are caught by the suite |
| 12 | Line counts | yes | **pass**: src +33/−182, net −149; the test file is 149 lines with 4 tests (5 at base) |
| 13 | No stale citations | yes | **pass** |
| 14 | No behaviour beyond the plan's Behaviour section | yes, read all three commits | **pass**, with one cosmetic message difference noted |

## Details

### 1. Full suite

I ran four suites in parallel:
- head, in the branch worktree;
- each of the three commits, in detached worktrees `TE/wt_<hash>` made by `git worktree add --detach`. I removed these afterwards and pruned them; `git worktree list` no longer shows them.

| Tree | Result | Log |
|---|---|---|
| head `94fca6b` (branch worktree) | 1207 passed, 47 skipped, 5 deselected, 1 xfailed in 221 s | `suite_head.log` |
| `ee9fe8a` | 1208 passed, 47 skipped, 5 deselected, 1 xfailed in 227 s | `suite_ee9fe8a.log` |
| `cfd89e6` | 1207 passed, 47 skipped, 5 deselected, 1 xfailed in 224 s | `suite_cfd89e6.log` |
| `94fca6b` (scratch worktree) | 1207 passed, 47 skipped, 5 deselected, 1 xfailed in 222 s | `suite_94fca6b.log` |

- In every run, the one `-rxX` line is `XFAIL tests/invariants/test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing] - C3: INV-11 ...`.
- No run has an XPASS, a FAILED or an ERROR.
- `ee9fe8a` has one extra test, `test_elements_at_bridge_entries_agree_with_the_ownership_table`, which commit 2 deletes along with the field. This matches the coder's table.

### 3. Collected ids

Command, in each worktree:

```
MPLBACKEND=Agg PYTHONPATH=src ../tangle-pack/env/bin/python -m pytest --collect-only -q -p no:cacheprovider | grep '::' | sort
```

This gives 1256 ids on base and 1255 on the branch. `diff base branch` prints exactly:

```
< test_partition_elements.py::test_element_for_rejects_an_unknown_endpoint_name
< test_partition_elements.py::test_element_for_resolves_a_bridge_missing_from_the_snapshot
< test_partition_elements.py::test_element_for_resolves_both_ends_of_a_branch_wrapping_bridge
< test_partition_elements.py::test_elements_at_bridge_entries_agree_with_the_ownership_table
> test_partition_elements.py::test_a_bridge_missing_from_the_snapshot_resolves
> test_partition_elements.py::test_both_ends_of_a_branch_wrapping_bridge_resolve
> test_partition_elements.py::test_owner_at_end_rejects_an_unknown_endpoint_name
```

That is one deletion and three renames, as the plan lists. The two id lists are `TE/ids_*.txt`.

### 2. Pinning tests on base and branch

`git diff --stat 5f0db5f..94fca6b -- test_bridge_class.py test_partition_family.py test_dual_walk.py` is empty, so the three files are unmodified.

In each worktree I ran:
- `test_bridge_class.py::test_a_missing_partition_names_the_branch` (`:311`)
- `test_bridge_class.py::test_duplicate_partitions_are_rejected` (`:332`)
- `test_partition_family.py::test_owner_of_intersection_matches_the_result_table` (`:109-110`)
- the whole of `test_dual_walk.py`, which includes the `_OwnerFamily` fake at `:873`

| Run | Base | Branch |
|---|---|---|
| The four selections above | 59 passed | 59 passed |
| The three whole files | 84 passed | 84 passed |

The missing-partition `ValueError` comes from a different place on each side:
- on base, from `_element_at_end`;
- on the branch, from `PartitionFamily.owner_of_intersection`, which runs when no result covers `(branch, row)`.

The test pins the type and the recovery, and both hold.

### 8. Agreement probe

I copied the harness to `TE/probe/`. The only edit is the output directory `D=` in `run.sh`; `probe.py` is byte-identical to the analyst's (`diff` is empty). The copy was needed so that the base run could not overwrite `b9b/probe/all_out_base.jsonl`.

Commands:
- `bash TE/probe/run.sh <branch> out_branch` (17 s)
- `bash TE/probe/run.sh <base> out_base_te` (17 s)
- `summ.py` on each result

| Run | Items | Errors | Bad | `(has_element_for, has_element_at_end)` | agree / ends | sided_agree / sided | it_agree / it_ends | it_all_ok |
|---|---|---|---|---|---|---|---|---|
| branch | 63 | 0 | 0 | (False, False) | 7090 / 7090 | 14180 / 14180 | 1776 / 1776 | 7090 |
| base (my run) | 63 | 0 | 0 | (True, True) | 7090 / 7090 | 14180 / 14180 | 1776 / 1776 | 7090 |
| analyst's `all_out_base.jsonl` | 63 | 0 | 0 | (True, True) | 7090 / 7090 | 14180 / 14180 | 1776 / 1776 | 7090 |

- There are no `A_exc`, `oracle_exc` or `*_disagree` counters in any run, so there were 0 exceptions.
- **Digests.** Per-item `readings_sha256`: analyst base == my base == branch, 63/63, no mismatch.
- **Raw readings.** `diff -r out_branch/readings out_base_te/readings` is empty, so all 63 readings files are identical.
- **What the branch run checks.** On the branch, `owner_at_end` (row) and `owner_of_intersection` (both sides) are checked against the law oracle `laws._element_ends`. On base, they are also checked against `element_for`, `_element_at_end` and `elements_at_bridge`.

### 9. 63-item outcome probe

I copied the harness to `TE/probe63/`; again only `D=` differs. I ran `run.sh <base> out_base` (19 s) and `run.sh <branch> out_branch` (19 s), then sorted each output with `cat */*.jsonl | sort`.

| Comparison | Result |
|---|---|
| `wc -l` | 63 lines each |
| `cmp out_base.all out_branch.all` | exit 0 |
| `cmp out_base.all <analyst>/probe63/out_base.all` | exit 0 |
| `grep -c '"error"'` | 0 on both |
| stderr | timing lines only |

### 11. Teeth

I made six throwaway copies by running `git archive 94fca6b | tar -x` into `TE/teeth_<n>`, planted one break in each, and deleted all six afterwards. Nothing was committed.

| Break | Planted change | `test_partition_elements.py` | Full suite |
|---|---|---|---|
| (a) neighbouring element | `owner_of_intersection` returns `result.ref(element_id + 1)` where that element exists | **3 failed, 1 passed** (as planned): `test_every_workbench_bridge_touching_the_branch_resolves`, `test_both_ends_of_a_branch_wrapping_bridge_resolve`, `test_a_bridge_missing_from_the_snapshot_resolves` | 51 failed |
| (b) `"third"` accepted | `endpoint_index` returns 1 for `"third"` | **1 failed, 3 passed** (as planned): `test_owner_at_end_rejects_an_unknown_endpoint_name` | 1 failed (the same test) |
| (c) extra: anchor rule dropped | `stable_key_at` ignores C7 (`if True: return ix.manifold_b_key`) | 4 passed, not caught here | 16 failed + 19 errors |
| (d) extra: a partition dropped | `HomotopyPartition.from_results` drops the first result | 3 failed, 1 passed | 314 failed + 167 errors |
| (e) extra: wrong row | `owner_at_end` reads the opposite row | 4 passed, not caught here | 57 failed |
| (f) extra: one end read twice | `bridge_classes`'s `ends` reads `"first"` twice | 4 passed, not caught here | 114 failed |

Which tests catch the three breaks the file misses:
- **(c):** `test_every_law_case_crossing_is_in_n`, `test_partition_laws`, `test_classes_laws`, `test_symbolic_laws`, `test_base_chords_are_inside_the_lobe_image` (x2), the `p2inv10` tests, `test_members_share_the_P1_itinerary_under_blasting`, the deep-run laws, plus setup errors across the law tier for one case.
- **(e):** `test_classes_laws` (x11), `test_symbolic_laws` (x11), `test_iterated_laws`, the deep-run laws, `test_classes_are_keyed_by_element_pairs`, the `test_dual_walk` itinerary tests, among others.
- **(f):** `test_classes_laws`, `test_iterated_laws`, `test_symbolic_laws` and `test_case_sanity_laws` (x11 each), `test_every_image_crossing_abuts_an_empty_stretch`, `test_table_lookups_symbols_and_report`, among others.

**Assessment.** The rewritten file pins the any-side lookup (`owner_of_intersection` on the `stable_key_at` branch) and the endpoint-name check, and it bites on both.
- The **row** reading of `owner_at_end` (e) and the **anchor branch** rule (c) are not this file's job. At base, `element_for` took the side as an argument, so the old file never tested the row either.
- Both (c) and (e) are caught across the law tier (the classes law reads the oracle `laws._element_ends`) and the unit tests.
- I added no test: anything I added would duplicate those guards.

Logs: `TE/suite_teeth_{a..f}.log`. After cleanup, both worktrees show an empty `git status --short`, the branch is at `94fca6b` and base is at `5f0db5f`.

### Rubrics 4, 5, 6, 7, 10, 12, 13

**4. File list.** `git diff --stat 5f0db5f..94fca6b`:

```
BridgeClass.py 66 | PartitionFamily.py 4 | StablePartition.py 57 | TopologyResults.py 32 | Trellis.py 56 | tests/unit/topology/test_partition_elements.py 105
6 files changed, 76 insertions(+), 244 deletions(-)
```

These are exactly the plan's six files. `laws.py`, `cases.py`, `tests/invariants/`, `CLAUDE.md`, `docs/` and `.refactor/` are absent.

**5. Greps, run in the branch worktree.** All three print nothing (exit 1):
- `grep -rn "element_for\b\|_elements_at_bridge\|elements_at_bridge\|_element_at_end" src tests`
- `grep -rn "partition_element_for" src tests`
- `grep -rn "def test_element_for" tests`

**6. Deleted names.** Run with `python -I` from `/tmp`, with the branch's `src` first on `sys.path` (confirmed via `BC.__file__`):
- `hasattr(Trellis, "element_for")` → False
- `"elements_at_bridge" in fields(StablePartitionResult)` → False
- `hasattr(BridgeClass module, "_element_at_end")` → False
- `_index_partitions` is still present → True

**7. The class rule calls the survivor.**
- `BridgeClass.py:598-599` reads `family.owner_at_end(trellis, bridge_id, "first")` and `family.owner_at_end(trellis, bridge_id, "second")`.
- `:582` holds the function-level `from .PartitionFamily import HomotopyPartition`, and `:584` holds `HomotopyPartition.from_results(partitions, trellis=trellis)`.
- `BridgeClass.py` has no `from .StablePartition` import. The only `StablePartition` mentions left are docstring cross-references (`:6`, `:541`) and `StablePartitionResult` from `TopologyResults`.
- The `PartitionFamily.py` diff has two hunks, at `@@ -137` and `@@ -895`.
- `owner_at_end` is unchanged: its source segment, found by `ast`, has the same sha256 on base (`:343-367`) and on the branch (`:342-366`), `fe690ebff0313db8…`. Only the line numbers moved, by the one deleted import line.

**10. No import cycle.** Each runs with `python -I` from `/tmp` and `sys.path.insert(0, '<branch>/src')`, and all three succeed:
- `import tanglepack.topology.BridgeClass, tanglepack.topology.PartitionFamily` → ok
- the reverse order → ok
- `import tanglepack.loom` → ok, loaded from the branch's `src`

**12. Line counts.**
- `git diff --numstat 5f0db5f..94fca6b -- src` gives:
  - BridgeClass.py +18/−48
  - PartitionFamily.py +0/−4
  - StablePartition.py +7/−50
  - TopologyResults.py +8/−24
  - Trellis.py +0/−56
- The total is **+33 / −182, net −149**, inside [−140, −165].
- `wc -l` on the test file gives 149 lines. `grep -c '^def test_'` gives 4 on the branch and 5 on base.

**13. Stale citations.** In the six files, a case-insensitive grep for each of the following finds nothing in prose:
- `element_for`
- `elements_at_bridge`
- `_element_at_end`
- `partition_element_for`
- `two lookup tables`
- `both tables`
- `carries no stable manifold key`

The remaining `manifold_b_key` hits are code reads, none of them in `BridgeClass.py`. Reading the `bridge_classes` docstring:
- the body names `Trellis.stable_key_at` (C7) and `owner_at_end`;
- the Raises clause says "if a covering result does not own a crossing";
- the Note calls an unowned crossing `owner_of_intersection`'s ValueError.

### 14. Read of the three commits

- **`ee9fe8a`:** deletes `Trellis.element_for` (52 lines with its header) and three imports, and rewrites the test file as the plan says. No `src` path besides the deleted method changes.
- **`cfd89e6`:** deletes the `elements_at_bridge` field, its helper, both writers and the import, and rewords the docstrings to describe one table. The field is the last one, so positional construction is unaffected.
- **`94fca6b`:** `bridge_classes` builds a `HomotopyPartition` and reads both ends through `owner_at_end`. The base's `_element_at_end` did the same steps in the same order: `row_of_end`, then `stable_key_at` on `unstable_key_of`, then the `(branch, row)` result, then `element_of_intersection`. `owner_at_end` performs exactly those steps.
  - **Duplicate partitions.** These still raise `ValueError`, because the `PartitionFamily` constructor calls the same `_index_partitions`.
  - **Generator input.** `partitions` is consumed once in both versions, so a generator argument behaves alike.

Semantic differences, all on error paths no corpus item reaches:

1. **Unowned crossing (plan Q2, approved):** `AssertionError` becomes `ValueError`.
2. **Dead `branch_key is None` branch removed:** it was unreachable.
3. **Not listed in the plan, cosmetic:** the missing-partition `ValueError` text is now `PartitionFamily.owner_of_intersection`'s, "crossing N is on no 'row'-partitioned branch (o, b) of this homotopy partition family". It still names the branch, as the docstring promises. It drops `_element_at_end`'s hint "partition that branch (or pass its trellis's stable_partitions)". The type is unchanged, and the pinning test checks only the type.

Rubrics 8 and 9 confirm that no computed output changed: 63/63 identical readings and byte-identical outcomes.

## Not checked / caveats

- **Coverage guard.** I did not run `tests/_tools/coverage_guard.py`. It is not in the rubric.
- **Commit trailers.** The coder's commits carry `Co-Authored-By: Claude Opus 5.5` rather than the `Fable 5.1` the spawn message named. The coder reports this as deviation 1. I made no commit, so the trailer question did not arise for me.
