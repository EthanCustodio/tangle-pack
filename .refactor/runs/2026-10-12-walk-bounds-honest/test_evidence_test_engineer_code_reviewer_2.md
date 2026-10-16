# Test evidence: B6a `walk-bounds-honest` (test engineer, round 1)

*2026-10-12. Plan: `artifacts/plan_planner_human_2.md`. Coder report: `artifacts/change_coder_test_engineer.t1.md`.*

- **Base:** `76d0710`, detached worktree `/home/dezhu/Development/TanglePack/tangle-pack-walk-bounds-honest-base`. It is clean; every copy I made into it was restored with `git checkout`.
- **Branch:** `refactor/walk-bounds-honest`, worktree `/home/dezhu/Development/TanglePack/tangle-pack-walk-bounds-honest`. The head is now `8a1cdc3`, which is my one test commit on top of the coder's `3465c60`.
- **Python:** `/home/dezhu/Development/TanglePack/tangle-pack/env/bin/python`, with `MPLBACKEND=Agg`. The worktrees have no `env/` of their own.
- **Scripts:** in `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/te/`:
  - `pin.py`, the dump;
  - `run.sh`, which runs every target;
  - `breaks.py`, the named source mutations;
  - `base/` and `branch/`, the outputs.

## Summary

| Rubric item | Result |
|---|---|
| 1. Suite | PASS: `1162 passed, 38 skipped, 5 deselected, 1 xfailed`. The only xfail is KNOWN_ISSUES `orientation_reversing` (C3). |
| 2. New tests fail on base and pass on branch | PASS: 6 of 6 fail on base; the 78 others pass on base. |
| 3. Diff stat is exactly the four files | PASS: still the four files after my commit. |
| 4. Behaviour pin | PASS: the base-vs-branch dump diff is empty on every target. The only diffs are (a) `heteroclinic` at `max_walks=1` and `2`, which is rubric 6, and (b) which nodes the k32 blast-6 span keeps, which is the second tier at work, with the word unchanged. |
| 5. k32 regression pair | PASS on the branch (`2 passed`). It FAILS (`42 == 45`) when the second tier is removed. |
| 6. heteroclinic at `max_walks=2` | PASS: True on base, False on the branch. The reason is "search budget reached before uniqueness was established". |
| 7–10. Unit tests | PASS. Each has a break-it check below. |
| 11. Greps | PASS. |

## 1. Pinning the old behaviour (rubric 4)

### How the dump works

`pin.py ROOT TARGET` puts `ROOT/src` and `ROOT/tests` first on `sys.path`. It asserts that `cases`, `DualWalk` and `SymbolicDynamics` are imported from ROOT. For each target it prints:
- `is_reliable`;
- `len(sd.table)`, the closed-table size;
- a sha1 over the sorted, letter-free per-class rules, then every rule in full.

Each rule is the tuple below. Everything in it is a positional `ElementRef` / `BridgeClass` label, never a session letter.

```
(class label, itinerary labels | None, word as ((class label, direction), ...),
 unresolved_reason, ambiguous, (search.status, search.truncated) | None)
```

**Wrapping `_nodes_meeting`.** `pin.py` also wraps `DualWalk._nodes_meeting` in both worktrees. For every call where NO node overlaps the span by more than its tolerance (the path B6a changes), it records:
- the span and the tolerance;
- every overlap;
- the nodes the function returned.

On base that is the old "every node" fallback; on the branch it is the second tier.

**Targets** (`run.sh`):
- the 11 `cases.LAW_CASES`;
- `p3_6_blasts`, `p3_16`, `k28_8` (`build_k28(blasts=8)`);
- `heteroclinic_mw2`, which runs `symbolic_dynamics(..., max_walks=m)` for m in 1, 2, 8 and 64;
- every step of all five `cases.SERIES`.

The series loop is the `cases.blast_series` loop, including the `_recomputed(workbench)` after each step. It therefore runs exactly the topology the regression tier tests, and it differs from the coder's script, which skipped the recompute.

### Commands

```
./run.sh /home/dezhu/Development/TanglePack/tangle-pack-walk-bounds-honest-base base
./run.sh /home/dezhu/Development/TanglePack/tangle-pack-walk-bounds-honest      branch
diff -r base branch -x '*.err'
```

Neither side's stderr had a traceback.

### Result: the full diff

There is nothing for any law case, `p3_6_blasts`, `p3_16`, `k28_8`, or any step of `k28_series`, `p3_series`, `inversion_series` or `nested_outer_series`. On `k32_series`, every step 0–10 has an identical header and rules hash, and identical rules. The diff contains exactly the following.

**`heteroclinic_mw2`, which is rubric 6 and the intended change:**

```
< heteroclinic max_walks=1: reliable=False closed=75 rules=a311352571fb
> heteroclinic max_walks=1: reliable=False closed=75 rules=ba9bbffb433c
<     ('A:p1@0.0/R#0 <-> A:p1@0.0/R#4', (...6 elements...), (...word...), None, False, ('unique', True))
<     ('A:p1@0.0/R#0 <-> B:p3@2.0/L#0', (...6 elements...), (...word...), None, False, ('unique', True))
>     ('A:p1@0.0/R#0 <-> A:p1@0.0/R#4', None, (), 'search budget reached before uniqueness was established', False, ('unique', True))
>     ('A:p1@0.0/R#0 <-> B:p3@2.0/L#0', None, (), 'search budget reached before uniqueness was established', False, ('unique', True))
< heteroclinic max_walks=2: reliable=True closed=75 rules=11ede46a52fc
> heteroclinic max_walks=2: reliable=False closed=75 rules=678f5f8a146c
<     ('B:p3@2.0/L#0 <-> B:p3@2.0/L#2', (...6 elements...), (...word...), None, False, ('unique', True))
>     ('B:p3@2.0/L#0 <-> B:p3@2.0/L#2', None, (), 'search budget reached before uniqueness was established', False, ('unique', True))
```

- At `max_walks=8` and `64` the hashes are identical (`1c7fdf65dabb`, which equals the default `heteroclinic` run).
- At `max_walks=1`, the two classes whose status is `"truncated"` (no walk) were already unresolved on base and stay so.
- The closed size stays 75 at every budget.

**`series:k32_series` blast 6: only the `SPAN` probe lines differ.** These are the nodes `_nodes_meeting` returned for the twice-landed span `A:p1@0.0/L#0`, `[1.270716632, 1.273790549]`, tolerance `0.012707166`. The element has 22 stable nodes and the span overlaps exactly one of them by more than zero (`0.003073918`, node `[1.270241706, 1.273790549]`):

```
< SPAN (... ((0.0, 0.001030459), (0.001030459, 0.005311691), ... 22 nodes ..., (1.273790549, 6.190159312)))
> SPAN (... ((1.270241706, 1.273790549),))
```

Base walked from all 22 nodes and the branch walks from the one the span lies in. The resulting rule for that class (itinerary, word, status) is byte-identical.

### Where the changed path fires

| | Base | Branch |
|---|---|---|
| k32 blast 6 | 2 | 2 |
| k32 blast 8 | 8 | 8 |
| Any other target | 0 | 0 |

- The 8 firings at blast 8 are on `R#48` and `R#38`, elements with ONE node each, so base's "every node" and the branch's second tier return the same node. That is why the blast-8 lines do not differ.
- Every one of the 10 spans has exactly one node with positive overlap. The empty-list (`span_meets_no_node`) branch never fires on real data.

**Headers** (branch, identical on base), with law-case and series closed sizes for the record:
- k10 2;
- k28_one_blast 3; k28_two_blasts 4; k28_8 4;
- p3 5; p3_6_blasts 43; p3_16 43;
- nested 9; nested_inner_blasts 9;
- inversion 4 / 7 / 8; inversion_period2 6;
- heteroclinic 75;
- k32_series: 3, 3, 9, 9, 9, 9, 18, 18, 45, 45, 45.

All are reliable except where base is unreliable too: `k28_series` and `k32_series` blast 0, and `nested_outer_series` blast 0.

## 2. New and rewritten tests, and what each guards

**The coder's four** (kept unchanged):
- `tests/unit/topology/test_dual_walk.py::test_a_span_narrower_than_its_error_starts_from_the_edge_it_lies_in` (rubric 8). Guards the second tier of `_nodes_meeting`: a scaled span `(1.2, 1.205)` with tol 0.012 starts from `[f3]` only.
- `tests/unit/topology/test_dual_walk.py::test_a_span_meeting_no_edge_leaves_the_walk_undecided` (rubric 7). It replaces `test_a_span_meeting_no_edge_falls_back_to_every_edge`, which is gone. Guards the `span_meets_no_node` return: no walks, one WARNING.
- `tests/unit/topology/test_symbolic_dynamics.py::test_members_holding_different_markers_leave_the_class_unwalked` (rubric 9). Guards the `MembersDisagree` raise and the catch before walking. The spy on `shortest_walks` counts zero calls.
- `tests/unit/topology/test_symbolic_dynamics.py::test_a_truncated_search_decides_nothing_and_is_no_wall` (rubric 10). Guards both the truncation check and the wall gate.

**Mine**, added in commit `8a1cdc3`:
- `tests/unit/topology/test_symbolic_dynamics.py::test_a_budget_too_small_to_decide_a_walk_leaves_its_class_unresolved[1|2]`. This is the truncation rule on REAL data: `build_heteroclinic()`, then `symbolic_dynamics(..., max_walks=1|2)`.
  - Every class whose search is `truncated` has `itinerary is None` and no symbols.
  - A truncated class that has a candidate walk carries the "budget" reason.
  - The table is not reliable.
  - It asserts its own precondition ("the heteroclinic tie needs more than two walks"), so it cannot pass vacuously.
  - It pins no letter, word or count. Heteroclinic is the only cheap case with a truncation at a budget of 1–3: the probe over k10, k28 one and two blasts, p3, inversion and heteroclinic found none elsewhere.
  - Runtime is about 1 s.

### On base

Command:

```
cp <branch>/tests/unit/topology/test_{symbolic_dynamics,dual_walk}.py <base>/tests/unit/topology/
pytest -q tests/unit/topology/test_symbolic_dynamics.py tests/unit/topology/test_dual_walk.py
git checkout tests/unit/topology
```

Result: `6 failed, 78 passed`.

```
FAILED test_symbolic_dynamics.py::test_members_holding_different_markers_leave_the_class_unwalked   (AttributeError: the spy returned None, the class WAS walked)
FAILED test_symbolic_dynamics.py::test_a_truncated_search_decides_nothing_and_is_no_wall             (itinerary is not None)
FAILED test_symbolic_dynamics.py::test_a_budget_too_small_to_decide_a_walk_leaves_its_class_unresolved[1]  (truncated class has an itinerary)
FAILED test_symbolic_dynamics.py::test_a_budget_too_small_to_decide_a_walk_leaves_its_class_unresolved[2]  (same)
FAILED test_dual_walk.py::test_a_span_narrower_than_its_error_starts_from_the_edge_it_lies_in       ([f0, f3] != [f3])
FAILED test_dual_walk.py::test_a_span_meeting_no_edge_leaves_the_walk_undecided                     ('unique' == 'span_meets_no_node')
```

**On the branch:** the two modules pass (`84 passed`).

## 3. Break-it checks (on the branch, each restored with `git checkout -q src`, `git diff --quiet src` confirmed)

`breaks.py NAME` makes one asserted-unique text replacement. Unless noted, each break is followed by `pytest -q tests/unit/topology/test_symbolic_dynamics.py tests/unit/topology/test_dual_walk.py`.

| Break | Mutation | Result |
|---|---|---|
| `no_second_tier` | `meeting = [... o > tol]` (second tier dropped) | 1 failed: `test_a_span_narrower_than_its_error_starts_from_the_edge_it_lies_in`. **Real data:** the k32 pair below fails too. |
| `no_first_tier` | `meeting = [... o > 0]` (tolerance gate dropped) | 1 failed: `test_the_walk_leaves_from_the_edges_the_landing_span_meets[scaled_end_past_the_neighbour]` (existing test) |
| `old_fallback` | `return meeting or nodes` (the base behaviour) | 1 failed: `test_a_span_meeting_no_edge_leaves_the_walk_undecided` |
| `no_empty_return` | delete `if not start_nodes or not goal_nodes: return WalkSearch(..., "span_meets_no_node")` | 1 failed: `test_a_span_meeting_no_edge_leaves_the_walk_undecided` |
| `no_raise` | `lobe_image_markers` no longer raises `MembersDisagree` (reads the first lobe, as on base) | 1 failed: `test_members_holding_different_markers_leave_the_class_unwalked` |
| `walk_on_disagree` | catch `MembersDisagree`, set `wanted = None` and walk anyway (the trap the plan names) | 1 failed: `test_members_holding_different_markers_leave_the_class_unwalked` |
| `no_truncated_check` | `elif search.truncated:` becomes `elif False:` | 3 failed: `test_a_truncated_search_decides_nothing_and_is_no_wall`, `test_a_budget_too_small_..._unresolved[1]`, `[2]` |
| `no_wall_gate` | drop `and cd.itinerary is not None` from the wall gate | 1 failed: `test_a_truncated_search_decides_nothing_and_is_no_wall` (`walls_b == []`) |

### `no_second_tier` on the real regression pair (rubric 5's guard)

Command:

```
pytest -q "tests/regression/test_class_count_under_blasting.py::test_closed_alphabet_constant_while_holes_constant[k32_series]" \
          "tests/regression/test_class_count_under_blasting.py::test_class_count_constant_while_holes_constant[k32_series]"
```

- Unbroken branch: `2 passed in 27.89s`.
- With `no_second_tier`: `1 failed, 1 passed`.

```
FAILED ...::test_closed_alphabet_constant_while_holes_constant[k32_series]
E   ... 8: 47 holes, 34 (30/4), 42 closed ...; 9: 47 holes, 34 (30/4), 45 closed ...
E   assert 42 == 45
```

This is exactly the plan's prediction that `return []` alone would take blast 8 from 45 to 42. The regression tier catches a second tier that goes missing on real data.

## 4. Full suite (branch, at `8a1cdc3`)

```
cd /home/dezhu/Development/TanglePack/tangle-pack-walk-bounds-honest
MPLBACKEND=Agg /home/dezhu/Development/TanglePack/tangle-pack/env/bin/python -m pytest -q -rxX -p no:cacheprovider
XFAIL tests/invariants/test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing] - C3: INV-11: ...
1162 passed, 38 skipped, 5 deselected, 1 xfailed in 197.95s
```

The coder had 1160; my parametrised test adds 2. No existing test was weakened or deleted, apart from the plan-mandated replacement of `test_a_span_meeting_no_edge_falls_back_to_every_edge` (rubric 7).

## 5. Diff stat and greps (rubric 3, 11)

```
git diff --stat 76d0710
 src/tanglepack/topology/DualWalk.py           |  89 +++++++++++-------
 src/tanglepack/topology/SymbolicDynamics.py   |  51 ++++++++--
 tests/unit/topology/test_dual_walk.py         |  21 ++++-
 tests/unit/topology/test_symbolic_dynamics.py | 129 +++++++++++++++++++++++---
 4 files changed, 231 insertions(+), 59 deletions(-)
```

- `grep -rn "may leave from any of them" src`: nothing.
- `grep -rn "the non-crossing rule reads" src`: nothing.
- `grep -n span_meets_no_node DualWalk.py`: the Literal at `:566`, the return at `:1123`, and docs at `:68`, `:590` and `:1098`.

## 6. What I could not test, and why

- **`MembersDisagree` on real data.** Site (iii) never fires on any law case, deep run or series step; the plan and analyst say so, and my dumps found no "different markers" reason anywhere. Its only evidence is the hand-built stub test (rubric 9) and its two break checks.
- **`span_meets_no_node` on real data.** It never fires: every one of the 10 real firings of the empty-first-tier path has exactly one node with positive overlap. Only the hand-built graph test exercises it.
- **A truncated `trivial` search.** It is covered by the same `elif search.truncated` line as `unique`, but no fixture produces a truncated `trivial`. `no_truncated_check` is caught through the `unique` path only.
- **The wall gate on real data.** No real case has a truncated class followed by a disjoint-branch class whose tie a wall would break, so only the unit test (rubric 10) and its break check guard it.
- **`MembersDisagree` is not exported** from `tanglepack.topology` (`__init__.py` is out of scope). The tests reach it through `SD.MembersDisagree` / `tanglepack.topology.DualWalk`. This is a note for the reviewer, not a test gap.
