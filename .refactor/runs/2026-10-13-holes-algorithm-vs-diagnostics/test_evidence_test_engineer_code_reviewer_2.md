# Test evidence: B4a `holes-algorithm-vs-diagnostics` (test engineer, round 1)

## Where

| | |
|---|---|
| Base | `ae18859` (worktree `../tangle-pack-holes-algorithm-vs-diagnostics-base`, detached, clean) |
| Branch | `refactor/holes-algorithm-vs-diagnostics` in `../tangle-pack-holes-algorithm-vs-diagnostics` |
| Coder's commits | `95ea780` (move), `728dcbb` (measured side on demand), `fc9f9a4` (tests) |
| My commit | **`f166573`**: 2 new tests in `tests/unit/topology/test_production_checks.py` |
| Python | `/home/dezhu/Development/TanglePack/tangle-pack/env/bin/python` |
| Scratch files | `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/` (`b4a_compare.py`, `cmp_coarse.py`, `rubric5.py`, `rubric9.py`, `breaks.py`, outputs) |

How code is picked up:
- pytest run inside a worktree imports that worktree's `src` (`pythonpath = ["src", "tests"]`).
- Scripts run with `PYTHONPATH=<wt>/src:<wt>/tests`. I checked `tanglepack.topology.__file__` for each worktree; each resolves to its own tree.

Both worktrees are clean at the end. The base has nothing left over, and the branch has nothing uncommitted.

## Verdict

**The move is an exact equality.** I re-ran every rubric check myself, and all 13 pass.

Two gaps in the coder's tests, both now closed by `f166573`:

1. **`side_source` had no teeth.** Breaking it (always `"carried"`, or inverted) failed no test. The coder's tests use it only as a filter, and both mutations still pass that filter.
2. **Rubric 9 had no test.** It existed only as a coder script.

Full suite on the branch tip `f166573`:

```
MPLBACKEND=Agg env/bin/python -m pytest -q -rxX
XFAIL tests/invariants/test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing] - C3: INV-11 ...
1173 passed, 38 skipped, 5 deselected, 1 xfailed in 196.94s
```

- At `fc9f9a4`, before my commit, the same run gave 1171 passed, 38 skipped, 5 deselected, 1 xfailed (198.75 s).
- My commit adds exactly 2 tests.
- The only xfail is the `cases.KNOWN_ISSUES` entry `orientation_reversing`.
- I weakened or deleted no existing test.

## Rubric, re-checked

### 1. Full suite, and `-rxX` shows only KNOWN_ISSUES

**Pass.** Results are above.

### 2. Pinning files pass on base and branch

**Pass.**

On base, with the base versions of the files:

```
cd ../tangle-pack-holes-algorithm-vs-diagnostics-base && MPLBACKEND=Agg ../tangle-pack/env/bin/python -m pytest -q \
  tests/unit/topology/test_topology_invariants.py tests/unit/topology/test_production_checks.py \
  tests/unit/topology/test_stable_partition.py tests/plotting/test_plot_topology.py \
  tests/regression/test_coarse_blasts_rows_by_signs.py tests/invariants
227 passed, 38 skipped, 1 xfailed in 21.11s
```

- The `tests/invariants` tier exercises `helpers/laws.py`.
- On the branch, the same set is inside the green full suite.
- The branch run used for the break checks gives 229 passed, 38 skipped, 1 xfailed. That is the same set plus my 2 tests.
- Files 6-11 differ between the two sides only in imports, the `_hole` helper and logger names.

### 3. `git diff --name-only ae18859..` lists exactly the eleven files

**Pass.** This includes my commit:
- `StablePartition.py`, `TopologyResults.py`, `Trellis.py`, `__init__.py`, `diagnostics.py`
- `laws.py`, `test_plot_topology.py`, `test_coarse_blasts_rows_by_signs.py`, `test_production_checks.py`, `test_stable_partition.py`, `test_topology_invariants.py`

`git diff --stat ae18859..`: 11 files changed, 853 insertions, 704 deletions.

### 4. No behaviour change: equality

**Pass.** I re-ran the coder's script, extracted verbatim from its report into `scratchpad/b4a_compare.py`.

What it dumps, over the 11 `LAW_CASES` plus `p3_coarse_blasts`, on `topology().trellis` and on every per-fixed-point trellis:
- every hole `(origin, iterate, bounding_ids, near_intersection_id, bridge_side, openings, coords)`;
- the measured side of every propagated hole: base reads `h.measured_side`, branch calls `measured_sides(trellis, trellis.holes)`;
- `bridge_row_diagnosis` of every hole bridge, field by field;
- the `geometry_flags` and `side_diagnoses` `describe()` strings;
- the law-tier fingerprint, minus `products` (object ids).

```
PYTHONPATH=<base>/src:<base>/tests     MPLBACKEND=Agg python b4a_compare.py base.txt
PYTHONPATH=<branch>/src:<branch>/tests MPLBACKEND=Agg python b4a_compare.py branch.txt
sed -E 's/ at 0x[0-9a-f]+>/>/g' {base,branch}.txt > *.norm.txt ; diff base.norm.txt branch.norm.txt
diff exit 0
74385b0439e6dec2619241cef8c9ca2c  base.norm.txt
74385b0439e6dec2619241cef8c9ca2c  branch.norm.txt   (1258 lines; 435 hole lines, 382 measured, 240 row diagnoses)
```

- The md5 equals the coder's.
- The raw, un-normalised diff is 24 lines: the 12 fingerprint lines on each side. Each differs only by a `FixedPoint ... at 0x...` address.

Per-case `(row flags, side flags, measured)` counts, printed identically on both sides:

| Case | Counts |
|---|---|
| k10 | (0,0,2)×2 |
| k28_one_blast | (0,0,2)×2 |
| k28_two_blasts | (0,0,3)×2 |
| p3 | (0,0,6)×2 |
| nested, nested_inner_blasts | (0,0,9), (0,0,3), (0,0,6) |
| inversion | (0,0,4)×2 |
| inversion_one_blast | (0,0,6)×2 |
| inversion_two_blasts | (0,0,6)×2 |
| inversion_period2 | (0,0,8)×2 |
| heteroclinic | (0,0,69), (0,0,3), (0,0,40) |
| **p3_coarse_blasts** | **(2,20,80)×2** |

**The comparison can fail.** It is non-empty where it matters: 2 row flags and 20 measured flags on `p3_coarse_blasts`. I also broke the branch and re-ran the comparison on `p3_coarse_blasts` alone (`cmp_coarse.py`), diffing against the base run alone:

| Break (branch `diagnostics.py`) | Lines differing from base | What shows |
|---|---|---|
| none (control) | 0 (`diff -q` equal) | |
| `side_source` inverted | 160 | the side diagnoses read `(signs)` for propagated holes |
| `measured_sides`: `for _ in range(steps)` → `range(0)` (no backward map) | 408 | measured entries change; flags 20 → **40** sides |

Each break was restored with `git checkout -- src/tanglepack/topology/diagnostics.py`, then checked with `git status --short -- src` (empty).

**Caveat on comparing runs.** A single-case run and the all-cases run of `p3_coarse_blasts` differ in the last ulp of a few crossing cdists. That moves the fingerprint and 4 `rows` lines per side. The cause is process history (which cases ran earlier in the process), not the code: base alone equals branch alone, and base all-cases equals branch all-cases. So a comparison must always pair runs of the same shape.

### 5. The 16 names are byte-identical at commit 1

**Pass.** `scratchpad/rubric5.py` runs `ast.get_source_segment` over:
- `git show ae18859:.../StablePartition.py`
- `git show 95ea780:.../diagnostics.py`
- `git show 95ea780:.../StablePartition.py`

Decorated classes have their decorator lines included.

```
identical: 16 /16; differing: []
still defined in StablePartition@95ea780: []
same order: True
```

### 6. `StablePartition.py` is at most 1,370 lines

**Pass.** `wc -l` gives 1354 (base 1947).

### 7. `diagnostics.py` names and imports

**Pass.**

Top-level names:
- the 16 names, with `_warn_geometry_flags` renamed `warn_geometry_flags` as plan step 2 says;
- plus `side_source`, `measured_sides`, `logger`, and the `SideSource` alias (moved from `TopologyResults`).

Imports:
- The only import from `StablePartition` is the top-level `from .StablePartition import bridge_for_pair, row_of_end, rows_of_bridge` (line 38).
- `StablePartition.py` imports `.diagnostics` only at `:713` (in `check_bridge_rows_consistent`) and `:738` (in `check_holes_share_bridge_side`). The other grep hits are docstring text.

Fresh interpreter: `python -I -c "import tanglepack.topology.diagnostics"` from `/tmp` with the branch `src` works.

### 8. `measured_side|side_source` grep

**Pass.**
- `grep -rn "measured_side\|side_source" src` hits only `topology/diagnostics.py`, 17 lines.
- The same grep on `TopologyResults.py` is empty.

### 9. The map stays off the assertion path

**Pass.** `scratchpad/rubric9.py` sets `diagnostics.measured_sides` to raise, then runs on `p3_coarse_blasts`:

```
diagnose=False ok, 84 holes
diagnose=True raised: the map was called
```

This is now also a test (new test 2 below).

### 10. `map_inv|dynamical_system` in `StablePartition.py`

**Pass.** The grep hits only `:955` and `:958`, both inside `_map_backward` (`:951`), which still exists.

### 11. `__all__` is unchanged

**Pass.** Sorted `__all__` has 101 names on both sides, md5 `94e1d08e27819322e67881712068f511` on both.

### 12. The test's logger line has teeth

**Pass.** I reverted `_LOGGER` in `test_coarse_blasts_rows_by_signs.py` to `"tanglepack.topology.StablePartition"` on the branch:
- With the revert, the test fails: `test_coarse_blasts_rows_by_signs.py:61: assert (0 - 0) > (2 * 0)`.
- After `git checkout --`, it passes again.

### 13. The I1/I2 phrases appear in `diagnostics.py`

**Pass.** "its orbit expects", "on its image sub-arc" and "approaches its two defining crossings" each appear exactly once.

## New tests (`f166573`, `tests/unit/topology/test_production_checks.py`, in the file list)

### 1. `test_side_source_names_the_step_that_made_the_hole`

- **What it does.** On p3, it wraps `StablePartition.punch_holes` and `propagate_reference_holes` (both called at call time from `Trellis.punch_holes`) to record the holes each returns. It then asserts both lists are non-empty and that `side_source` is `{"signs"}` for the first and `{"carried"}` for the second.
- **What it guards.** D2's meaning: direct holes read their side off the signs, and propagated holes carry it. It checks against the step that made each hole, not by re-stating `side_source`'s one-line definition. That avoids the tautology the plan dropped from `laws.py`.
- **Base.** It passes on `ae18859`. I copied it into the base worktree as a temporary `tests/unit/topology/test_b4a_base_pin.py`, reading the old field `hole.side_source` in place of `diagnostics.side_source(hole)`: 1 passed. The file was removed afterwards, and the base is clean.
- **Branch.** It passes.

### 2. `test_measured_sides_stay_off_the_assertion_path`

- **What it does.** It spies `diagnostics.measured_sides` with the existing `_spy` helper:
  - `punch_holes(diagnose=False)` makes 0 calls;
  - the default `punch_holes()` makes exactly 1 call.
- **What it guards.** Rubric 9: the inverse map runs only on the flag path, once.
- **Base.** Not runnable there. This property is NEW in commit 2: on the base the measured side was computed inside `_punch_in_bridge` on every punch, and `diagnostics` does not exist.
- **Branch.** It passes.

## Break-it checks (branch, each restored by `git checkout -- <file>`)

**Harness.** `scratchpad/breaks.py` applies one edit and runs the six targeted test sets (the five files plus `tests/invariants`) with `-q --tb=no -rf`. It then runs `git checkout -- <file>` and asserts `git status --short -- src` is empty.

**Clean run.** Before any break: 229 passed, 38 skipped, 1 xfailed.

| # | Break | Result: failing tests |
|---|---|---|
| B1 | `side_source` returns `"carried"` always | **before `f166573`: 0 failed (gap)**; after: `test_side_source_names_the_step_that_made_the_hole` |
| B2 | `side_source` inverted | **before: 0 failed (gap)**; after: `test_side_source_names_the_step_that_made_the_hole` |
| B3 | `side_diagnoses` ignores `measured` (`measured_side=None`) | `test_measured_side_disagreement_is_a_flag_not_an_error`, `test_side_diagnoses_name_each_disagreeing_hole`, `test_punch_holes_flags_a_contradicting_measured_side` |
| B4 | `side_diagnoses` never counts unmeasured (`0 if True else`) | `test_unmeasured_propagated_holes_are_logged` |
| B5 | counts unmeasured on the I1 path too (None read as `{}`) | `test_unmeasured_propagated_holes_are_logged` |
| B6 | `geometry_flags` passes `measured=None` | `test_measured_side_disagreement_is_a_flag_not_an_error`, `test_punch_holes_flags_a_contradicting_measured_side`, `test_measured_sides_stay_off_the_assertion_path` |
| B7 | `measured_sides` maps no step back | `test_propagated_holes_carry_their_orbit_side` (and +20 flags in the rubric-4 comparison, above) |
| B8 | lazy import hoisted: `from .diagnostics import bridge_row_violation, bridge_side_violations` at `StablePartition` module top | collection error: `ImportError: cannot import name 'bridge_for_pair' from partially initialized module 'tanglepack.topology.StablePartition' (most likely due to a circular import)` |
| B9 | `Trellis.punch_holes`'s call-time `StablePartition` imports moved to module level | `test_punch_holes_checks_every_hole_bridge_once_and_every_orbit`, `test_punch_holes_raises_on_a_corrupted_hole_side` (the spies miss) |
| B10 | `diagnostics.logger` named `"tanglepack.topology.StablePartition"` | `test_unmeasured_propagated_holes_are_logged`, `test_punch_holes_flags_a_contradicting_measured_side`, `test_coarse_blasts_rows_by_signs` |
| B11 | `Trellis.punch_holes`: `if diagnose:` → `if True:` before `warn_geometry_flags` | `test_measured_sides_stay_off_the_assertion_path`, `test_coarse_blasts_rows_by_signs` |
| B12 | `bridge_row_violation` always returns None (the I2 core behind the lazy import) | `test_bridge_rows_with_equal_signs_raise`, `test_punch_holes_raises_on_equal_signs` |
| R12 | test-side `_LOGGER` reverted (rubric 12) | `test_coarse_blasts_rows_by_signs` |

**A note on B4 and B5.** My first versions of B4 and B5 (`0 * sum(...)`, and a bare `sum(...)`) iterated `None` and crashed: 22 failed, 199 errors. That proves only that the `None` guard is needed, so I replaced them with the quiet mutations in the table, each caught by one test.

## Not tested, and why

- **The fallback-orbit difference** (coder's behaviour note 4). Propagated holes of an orbit started from the propagation fallback would be "unmeasured" on the branch where the base measured them. The analyst found that the fallback never runs on any probed build, so no case can exercise it; B4b removes it.
- **Base pin for new test 2.** As noted above, it cannot be pinned on the base: the property was introduced by commit 2.
- **Commit 2 alone is red** (12 tests, per the coder). That follows the plan's split (src in 2, tests in 3). I checked commit 1 only through rubric 5's source equality, not by re-running the suite at `95ea780`. The coder reports 1171/38/5/1 there.
- **Out-of-scope breakage, left as the plan says:**
  - `scripts/refactor_proof/hole-side-and-own-blast/hole_sides.py` imports moved names.
  - The `numerics/geometry.py:7` docstring still names `StablePartition._build_oriented_bridge_polyline`.
