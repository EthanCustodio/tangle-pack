# Findings: B9b `one-element-at-end` (analyst, 2026-10-15, HEAD = `5f0db5f`)

**Tree state.** The tree is dirty, but only outside the code. Modified: `.refactor/`, `CLAUDE.md`, `docs/`. Untracked: run dirs, `artifacts/`, figures, two sweep scripts. `git diff HEAD -- src tests scripts` is empty. The change branches from `HEAD` and must leave all of that alone. *(verified)*

**Kind of request.** A refactor, specifically a deletion. Three of the four readers of "which element owns this bridge end" collapse into `PartitionFamily.owner_at_end`.

**Answer.** Yes, the three can fold into `owner_at_end` without changing anything.
- I re-ran the agreement probe on `HEAD`: 0 disagreements.
- I built the whole slice on a scratch copy. The suite passes, the probe readings are identical item by item, and the 63-item outcome probe is byte-identical.
- The real choice is where the class rule gets `owner_at_end` from. **I recommend a function-level import in `bridge_classes` (option B, §4).** It deletes `_element_at_end` outright, adds no new function and removes about 152 src lines.

Scratch artifacts are under `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b9b/` (called `b9b/` below):
- `head/`: an exported copy of `5f0db5f`.
- `cand_b/` and `cand_a/`: the two options, fully applied.
- `cand_b_src.diff`, `cand_b_tests.diff`, `cand_a_src.diff`.
- `probe/`: the agreement probe.
- `probe63/`: a copy of the B10 outcome probe.
- `smoke_draft.py`.

---

## 1. The four implementations at `HEAD` *(verified: grep + read)*

| Implementation | Lines | Reads | Side | Callers at `HEAD` |
|---|---|---|---|---|
| `PartitionFamily.owner_at_end` | `PartitionFamily.py:343-367` → `self.owner_of_intersection` `:369-398` | its own family (homotopy or iterated); the trellis passed | the row: `row_of_end` (D6) on the branch `stable_key_at(id, unstable_key_of(bid))` names (C7) | `SymbolicDynamics.py:642-643`, `DualWalk.py:1366`, `plotting.py:2303`; `test_dual_walk.py:877` binds it onto `_OwnerFamily` |
| `Trellis.element_for` | `Trellis.py:1379-1430` (with its `# ── partition elements` header) | `self.stable_partitions` (homotopy); first `(branch, side)` match | ANY side, passed in; returns `Optional[int]` | **no src caller**; tests `test_partition_elements.py:68,105,106,158,168` |
| `StablePartition._elements_at_bridge` | `StablePartition.py:1059-1091` (end of file) | the result being built; `trellis.bridges` (iterated: `minimal.sparse`) | each side; `None` for an end on another branch | writers `StablePartition.py:536-538` and `PartitionFamily.py:898-900`; import `PartitionFamily.py:140`; field `TopologyResults.py:415-427` (doc) and `:442-444`; **no src reader**; tests `:74,84,155` |
| `BridgeClass._element_at_end` | `BridgeClass.py:763-791` | the homotopy `partitions`, indexed by `_index_partitions` (`:583`) | the row | `ends` `:595-599`, used at `:604` and by `_non_loop_parent` `:740` |

No script, notebook, stress file or example references any of the three that are to be deleted. *(verified: grep)* `CLAUDE.md:139` still lists `element_for`; that is for B12.

**What `_element_at_end` does that `owner_at_end` does not, and the reverse.** Both compute the same triple: the endpoint id, `row_of_end`, and `stable_key_at(id, unstable_key_of(bid))`. Both then look the crossing up in `element_of_intersection` of the `(branch, row)` result. They differ only on the error paths (§5):
1. `_element_at_end` raises `ValueError("carries no stable manifold key")` when `branch_key is None`.
   - This path is dead on computed crossings: `stable_key_at` returns `manifold_b_key`, and `compute_intersections` registers every crossing as unstable × stable.
   - `owner_at_end` would hit `TypeError` there instead, because its message formats `branch_key[1:]` on `None`.
   - No test reaches this path *(verified: suite green on `cand_b`, where it is gone)*.
2. A missing `(branch, row)` result raises `ValueError` in both, with different wording.
3. An unowned crossing raises `AssertionError` in `_element_at_end` and `ValueError` in `owner_of_intersection`.

`owner_at_end` additionally serves the iterated family and image pairs (`DualWalk.py:1366`). `element_for` is the only one that answers on a side other than the row. `owner_of_intersection(id, side, stable_key_at(...))` already gives that, and B9a's Q6 ruled out new API for it.

**The one guard that disappears.** `_elements_at_bridge` asserted `element_id is not None` for every end of every trellis bridge on the branch. It did this on every partition build, iterated builds included. Deleting it loses no coverage:
- The law tier checks a stronger statement, on every law case: `check_partition_unique_owner` (`laws.py:1168-1183`, `element_of_intersection` keys == the branch's `intersection_ids`) and `check_iterated_unique_owner` (`:2214-2226`).
- At read time, `owner_of_intersection` raises.
- *(verified: read)*

## 2. The agreement probe, re-run on `HEAD` *(verified)*

- **Harness.** I copied the b9 harness to `b9b/probe/{probe.py,run.sh,summ.py}`; the b9 copies are untouched.
  - The deleted `session.partition_element_for` reading is dropped.
  - Every reading of a function B9b deletes is guarded with `hasattr`, so **the same script runs on base and branch**.
  - Each item also emits a sha256 of its readings and writes them out in full: row owners, both-side owners and iterated owners.
  - One fix was needed: `import tanglepack.topology.BridgeClass as BC` binds the *class*, because `topology/__init__` shadows the submodule name. The probe now uses `importlib.import_module`.
- **HEAD** (exported `5f0db5f`, 63 items: 16 builds + 47 series steps):
  - **7,090 / 7,090** ends agree on the row: `owner_at_end` = `_element_at_end` = `element_for` = `elements_at_bridge` = the law oracle `laws._element_ends`, each on the `stable_key_at` branch and the `row_of_end` side.
  - **14,180 / 14,180** (end, side) pairs: `owner_of_intersection` = `element_for` = `elements_at_bridge`.
  - **1,776 / 1,776** iterated ends: `owner_at_end` = the iterated `elements_at_bridge`. The full and sparse trellises read the same.
  - Iterated `owner_at_end` succeeds on all 7,090 trellis-bridge ends.
  - 0 exceptions, 0 "bad" rows. The counts equal the b9 run on `517aa01`.
- **Branch candidates** (`cand_b` and `cand_a`): the same 7,090 / 14,180 / 1,776 with the deleted functions absent, and **all 63 per-item reading digests equal base's**.
- **Rubric 3 needs adjusting.** The rubric says to compare against "base's `element_for`/`elements_at_bridge` readings recorded in `all_base.jsonl`". The b9 `all_base.jsonl` holds counts only, not readings. Point rubric 3 at `b9b/probe/` instead: run `run.sh <repo> out_base` and `run.sh <repo> out_branch`, then compare `readings_sha256` per item. Base digests are in `b9b/probe/all_out_base.jsonl`.

## 3. Outcome on the 63 items *(verified)*

I copied the B10 outcome probe to `b9b/probe63/`. It records words, matrices, class counts, `is_reliable`, walk statuses, kinds and unmatched members. I ran it on `head/` and on `cand_b/`. The sorted outputs are **byte-identical** (`cmp`), with 0 errors.

## 4. The import cycle, and where the shared reading lives

**The cycle.**
- `PartitionFamily.py:136` does `from .BridgeClass import _image_chain, _index_partitions`.
- `BridgeClass.py:97` does `from .StablePartition import row_of_end`.
- `topology/__init__.py` imports `BridgeClass` (`:59`) before `PartitionFamily` (`:82`).
- So any module-level `from .PartitionFamily import …` in `BridgeClass` fails. I added one to a scratch copy and both import orders raised `ImportError: cannot import name '_image_chain' from partially initialized module`. *(verified)*

Python always runs the package `__init__` first, so rubric 5's "both orders" go through the same `__init__` order. It is still the right gate, because it catches this failure.

**Option B (recommended): a function-level import in the class rule.** This is `cand_b`.
- In `bridge_classes`, replace `by_branch_side = _index_partitions(partitions)` with `from .PartitionFamily import HomotopyPartition` plus `family = HomotopyPartition.from_results(partitions, trellis=trellis)`.
- `ends` becomes `family.owner_at_end(trellis, bridge_id, "first"/"second")`.
- `_element_at_end` is deleted.
- `BridgeClass` then imports nothing from `StablePartition` (`row_of_end`, `Endpoint` and `endpoint_index` become unused, ruff F401).
- `from_results` runs `_index_partitions`, so the duplicate-partition `ValueError` (`test_bridge_class.py:332-337`) is unchanged.
- **Pros:**
  - It is literally one implementation: the class rule calls `owner_at_end`.
  - No new function.
  - `owner_at_end` is untouched, so the `_OwnerFamily` fake works unchanged.
  - Biggest deletion: src −152 net.
- **Cons:**
  - A function-level import. The codebase already has these, e.g. `StablePartition.py:635,660` (`from .diagnostics import …`).
  - The class rule now depends on the family module at call time. That matches CLAUDE.md, where the class rule is defined over the homotopy family.
  - `HomotopyPartition` is right because `bridge_classes` is only ever given homotopy results: `Topology.py:140`, `helpers/direct.py:102`, `test_symbolic_dynamics.py:873`. Its error message then says "of this homotopy partition family".

**Option A: a shared core beside `row_of_end` in `StablePartition.py`.** This is `cand_a`.
- Add `end_on_row(trellis, bid, endpoint) -> (id, row, branch_key)` and `owner_in(results, id, side, branch_key, what)` to `StablePartition.py`.
- `owner_at_end` becomes `self.owner_of_intersection(*end_on_row(...))`, so the fake still works.
- `owner_of_intersection` becomes `owner_in(self.results, …)`.
- The class rule calls `owner_in(by_branch_side, *end_on_row(...))`.
- **Pros:** no function-level import, and the module graph is unchanged.
- **Cons:**
  - Two new module functions in a 1,048-line file.
  - The class rule still does not call `owner_at_end`, so "one implementation" holds only at the level of the two helpers.
  - Smaller deletion: src −133 net.

**Option C, rejected: move `_image_chain` out of `BridgeClass` to break the cycle.** `_image_chain` is imported by `MinimalTrellis`, `diagnostics`, `SymbolicDynamics` and `PartitionFamily`, and by four test files from `BridgeClass` directly (`test_foreign_chords.py:31`, `test_evidence_diagnosis.py:35`, `test_closure_predicts_registered_images.py:13`, `_tools/resolution_probe.py:223`). Moving it would widen the slice well beyond its file list.

**Both A and B were verified on scratch copies.**
- Fresh interpreter, run with `-I` from an empty directory, `sys.path` set to the copy: `import tanglepack.topology.BridgeClass, tanglepack.topology.PartitionFamily` succeeds in both orders, as do `tanglepack` and `tanglepack.loom`.
- `cand_b` full suite: **1207 passed / 47 skipped / 5 deselected / 1 xfailed (C3)**, no XPASS. HEAD gives 1208, and the difference is the deleted field test.
- `cand_a`: 88 passed on `test_partition_elements`, `test_bridge_class`, `test_dual_walk` and `test_partition_family`; probe digests equal.

## 5. Error paths *(verified: read + grep of every `pytest.raises`)*

| Situation | `_element_at_end` (HEAD) | `owner_at_end` → `owner_of_intersection` | Pinned by |
|---|---|---|---|
| no `(branch, row)` result | `ValueError` ("no 'left' partition covers stable branch …") | `ValueError` ("crossing N is on no 'left'-partitioned branch … of this homotopy partition family") | `test_bridge_class.py:311-329` (bare `ValueError`, no message match): passes on both |
| crossing absent from the covering result | **`AssertionError`** (documented in the `bridge_classes` Note `:575-581`) | **`ValueError`** | `ValueError` pinned for `owner_of_intersection` at `test_partition_family.py:109-110` (`10**6`). The `AssertionError` is pinned **nowhere**, in tests or the law tier |
| `branch_key is None` | `ValueError` | `TypeError` (formatting `None[1:]`) | nothing; unreachable on computed crossings |

**What survives:** `ValueError`, through `owner_of_intersection`. It is already pinned, and it is what the other three src callers already see.
- The `bridge_classes` docstring must change to match: Raises `:569-574` drops "if a crossing carries no stable manifold key", and the Note at `:575-581` drops the unowned-crossing `AssertionError` while keeping the image-pair `AssertionError`.
- `cand_b` rewrote the Note but not the Raises clause. The planner should do both.
- `:542-544` still says "the endpoint's `manifold_b_key` picks the stable branch". That has been stale since C7 (`stable_key_at`). Reword it in the same edit.

## 6. `tests/unit/topology/test_partition_elements.py`: what each test asserts, and the rewrite

What each test asserts at `HEAD` *(read)*:
- **`:42` `test_every_workbench_bridge_touching_the_branch_resolves`.**
  - Setup: for every WORKBENCH bridge end whose `manifold_b_key` is a result's branch.
  - Asserts `element_for(bid, name, result.side) == result.element_of_intersection[end]`, on both sides because every result is visited.
  - Asserts `checked > 0`.
  - Plus `:74`: `set(result.elements_at_bridge) <= expected_keys`.
- **`:78` `test_elements_at_bridge_entries_agree_with_the_ownership_table`.** The field's own test: each entry is the owner, or `None` off-branch. **It goes with the field.**
- **`:96` `test_element_for_resolves_both_ends_of_a_branch_wrapping_bridge`.** For a bridge whose ends' `manifold_b_key` differ, it asserts on both sides that both ends are non-`None` and equal their own branch result's table entry.
- **`:119` `test_element_for_resolves_a_bridge_missing_from_the_snapshot`.**
  - Mechanism: run stages 3-5 by hand on `trellis(fp3)`, drop one bridge from `trellis.bridges`, then re-partition.
  - On both sides and both ends it asserts `dropped.id not in result.elements_at_bridge` (`:155`, which goes) and that `element_for` equals the table entry.
- **`:161` `test_element_for_rejects_an_unknown_endpoint_name`.** Asserts `element_for(bid, "third", "left")` raises `ValueError`.

**The rewrite** (`b9b/cand_b_tests.diff`; 149 lines, 4 tests):
- Module docstring `:9-11` now names `owner_of_intersection` / `stable_key_at` instead of `elements_at_bridge` / `element_for`.
- Two helpers:
  - `_owner(trellis, family, bid, end, side)` = `family.owner_of_intersection(end, side, trellis.stable_key_at(end, trellis.unstable_key_of(bid)))`.
  - `_expected(trellis, end, side)` = `result.ref(result.element_of_intersection[end])` for the result of the end's own `manifold_b_key` branch. It stays the independent derivation.
- `family = HomotopyPartition.from_results(trellis.stable_partitions, trellis=trellis)`: exactly the results `element_for` read.
- The rewritten tests compare whole `ElementRef`s (branch, side, id), which is stronger than the old `element_id`. Per test:
  - `:42` keeps its name, workbench enumeration, both sides and `checked`. Line `:74` goes.
  - `:96` asserts both sides, both ends, and `first.branch_key != second.branch_key`.
  - `:119` keeps the drop-and-repartition mechanism and rebuilds the family after `partition_stable_manifold()`. `:155` goes.
  - `:161` becomes `family.owner_at_end(trellis, bridge.id, "third")` raising `ValueError`. It raises through `endpoint_index` on the first line of `owner_at_end`.
- **Teeth** *(verified on two planted copies)*:
  - `owner_of_intersection` returning the neighbour element fails the three both-side tests (3 failed / 1 passed).
  - `endpoint_index` returning 1 for "third" fails the `ValueError` test (1 failed / 3 passed).
- **Collected ids in `cand_b` vs HEAD:**
  - Deleted: `test_elements_at_bridge_entries_agree_with_the_ownership_table`.
  - Renamed: `test_element_for_resolves_both_ends_of_a_branch_wrapping_bridge` → `test_both_ends_of_a_branch_wrapping_bridge_resolve`; `test_element_for_resolves_a_bridge_missing_from_the_snapshot` → `test_a_bridge_missing_from_the_snapshot_resolves`; `test_element_for_rejects_an_unknown_endpoint_name` → `test_owner_at_end_rejects_an_unknown_endpoint_name`.
  - Nothing else changes (1256 → 1255 ids).
  - Renaming is a choice (Q3). Rubric 2's `grep "element_for\b"` does **not** match `test_element_for_resolves…`, because `_` is a word character, so the old names would pass the grep.
- **`test_dual_walk.py:877`.** `_OwnerFamily` binds `PartitionFamily.owner_at_end` and supplies only `owner_of_intersection(crossing_id, side, branch_key)`, which asserts `branch_key == self.key`. Any option in which `owner_at_end` calls `self.owner_of_intersection(id, row, branch_key)` keeps it green. A and B both do; verified in both copies.

## 7. Docstrings citing deleted API *(verified: read)*

To fix in this slice (`cand_b` shows wording for all but the `bridge_classes` Raises clause):
- `StablePartition.py:23-32`: module notes ("two lookup tables … `elements_at_bridge` … `Trellis.element_for`").
- `StablePartition.py:499-504`: the `partition_stable_manifold` docstring.
- `TopologyResults.py:400-403` ("the two lookup tables … which elements sit at this bridge's two ends"), `:415-427` (the field doc, including the `partition_element_for` cite) and `:429` ("Both tables").
- `Trellis.py:1413`: goes with `element_for`.
- `BridgeClass.py:542-544` and `:569-581` (§5).
- `test_partition_elements.py:9-11, :47-51`.

**`ResonanceZone.py:295` and `:353` are NOT stale**, contrary to the request:
- `:295` says "the nested-zone classifier in `TangleSession` orders zones by [area]". `TangleSession.classify_bridge` (`TangleSession.py:749-773`) still does exactly that: `min(containing, key=area)`.
- `:353` cites `TangleSession.add_resonance_zones`, which exists at `TangleSession.py:696`.
- Recommendation: **leave both alone**, so `ResonanceZone.py` stays out of the file list.

## 8. What the law tier reads *(verified: grep of `tests/helpers/*.py`, `cases.py`, `invariants/`, `conftest.py`)*

- **None of the four, and nothing B or A touches.**
- The law tier reads `row_of_end` (`laws.py:75,1323,1874,1979,2200,2252-2253`), `element_of_intersection` (`:1177,1182,1197,1876,2226`), its own `_element_ends` (`:1868`) and `bridge_classes()` through the session.
- `helpers/direct.py:31,102` calls `bridge_classes(trellis, partitions)` directly. Its signature is unchanged in both options.
- `helpers/fakes.py:164` builds `StablePartitionResult` by keyword without `elements_at_bridge`.
- `test_exact_ownership.py:72` builds it positionally with 4 args. Since `elements_at_bridge` is the 5th field, this still works.

## 9. Line counts *(verified: `wc -l`, `cand_b` vs `head`)*

| File | HEAD | Option B | Option A |
|---|---|---|---|
| `Trellis.py` | 1626 | 1570 (−56: `element_for` + 3 unused imports) | −56 |
| `StablePartition.py` | 1091 | 1048 (−43: helper −33, writer −3, docstrings) | 1077 (−14) |
| `PartitionFamily.py` | 1077 | 1073 (−4) | 1066 (−11) |
| `BridgeClass.py` | 791 | 758 (−33) | 756 (−35) |
| `TopologyResults.py` | 922 | 906 (−16) | −16 |
| **src total** | | **+23 / −175 = −152** | **−133** |
| `test_partition_elements.py` | 168 | 149 (+43 / −62) | same |

Two pre-existing issues stay out of scope. `StablePartition.py:238` imports `Intersection` unused (ruff F401 on `HEAD`). Option A leaves `endpoint_index` unused in `PartitionFamily`.

## 10. The `plot_resonance_zones` smoke test

The rubric can carry one new id: it is deterministic and quick. I verified a draft in `cand_b` (`b9b/smoke_draft.py`):

```python
def test_plot_resonance_zones_shades_every_zone(k28_partitioned, fig_ax) -> None:
    session, _fp = k28_partitioned
    zones = list(session.resonance_zones.values())
    assert zones
    _fig, ax = fig_ax
    assert len(plot_resonance_zones(zones, ax=ax)) == len(zones)
    assert plot_resonance_zones([], ax=ax) == []
```

`k10_partitioned` has **no** resonance zone, and the draft fails on it. `k28_partitioned` has one, and the test passes in 0.23 s. Adding it would put `tests/plotting/test_plot_smoke.py` in the file list and has nothing to do with this slice's concern. Recommendation: **leave it for B12** (Q4).

## 11. For the planner

1. **Shape:** option B (§4). Delete `Trellis.element_for`, `_elements_at_bridge` with its field, writers and import, and `_element_at_end`. Add the function-level `HomotopyPartition` import in `bridge_classes`. Drop the unused imports in `Trellis.py:77,80,83` and `BridgeClass.py:97,100-101`. `cand_b_src.diff` is a complete reference, short of the `bridge_classes` Raises clause.
2. **Error path:** `ValueError` survives; update the `bridge_classes` docstring as in §5.
3. **Tests:** `cand_b_tests.diff` as in §6.
4. **Proof:**
   - The `b9b/probe/` digests base vs branch (63/63 equal) plus counts 7,090 / 14,180 / 1,776.
   - `b9b/probe63/` byte-identical.
   - Suite 1207 / 47 / 5 / 1 xfail (C3); 1208 if Q4 is answered yes.
   - Import check as in rubric 5.
   - Teeth as in §6.

**Questions only the author can answer (each with my default):**
- **Q1. Where does the class rule get its owner?**
  - (B) function-level import, `family.owner_at_end`: src −152, one implementation literally.
  - (A) `end_on_row` / `owner_in` beside `row_of_end`: no local import, −133, two new helpers.
  - *Default: B.*
- **Q2. Is the unowned-crossing `AssertionError` becoming `ValueError` acceptable?** Nothing pins the `AssertionError`; the `ValueError` is already pinned for `owner_of_intersection`. The invariant itself stays enforced: by `_element_of_intersection`'s own assertions when the partition is built, and by the law tier. *Default: yes.*
- **Q3. Rename the three surviving tests that cite `element_for` in their names?** This changes 3 collected ids. *Default: rename,* so no test name cites deleted API.
- **Q4. Add the `plot_resonance_zones` smoke test here?** *Default: no,* leave it for B12. It is unrelated to this slice and would widen the file list. If yes, use the `k28_partitioned` draft in §10.
- **Q5. `ResonanceZone.py:295,353`?** *Default: no change.* Both citations are still accurate (§7).

## Out of scope (noticed, not recommended here)
- `CLAUDE.md:139` lists `Trellis.element_for` (B12).
- `StablePartition.py:238` imports `Intersection` unused, pre-existing.
- `henon_blast_period_3.py:95` is still dead (author's rule: report, don't fix).
- `bridge_classes`'s `partitions: Iterable[StablePartitionResult]` could take a `HomotopyPartition` directly. That would change its signature for `Topology.py:140`, `helpers/direct.py:102` and many unit tests; not doing.

*Knowledge files updated:* `.refactor/knowledge/topology.md` (header, the users list under C7, and a new "B9b analyst" section) and `.refactor/knowledge/testing.md` (header, and a new "B9b analyst" section).
