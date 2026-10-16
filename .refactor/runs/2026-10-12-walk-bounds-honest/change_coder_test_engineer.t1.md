# Change: B6a `walk-bounds-honest` (coder, iteration 1)

*2026-10-12. Job `change`. Plan: `artifacts/plan_planner_human_2.md`, approved as written (`plan_ok_planner_human.json`; question 1 takes its default: `lobe_image_markers` raises `MembersDisagree`, caught in `_resolve_itinerary` before walking).*

- **Base commit:** `76d0710dd07b9fb9d0a6f2459a01561d7884de20` (author's `HEAD`, `main`). The author's uncommitted `docs/`, `.refactor/` and `artifacts/` edits were not touched.
- **Branch:** `refactor/walk-bounds-honest`, worktree `/home/dezhu/Development/TanglePack/tangle-pack-walk-bounds-honest`.
- **Base worktree:** `/home/dezhu/Development/TanglePack/tangle-pack-walk-bounds-honest-base` (detached at `76d0710`, clean).

## Commits (one per plan step)

```
3465c60 Docs: is_reliable, Dev Notes and WalkSearch say an unestablished walk leaves its class unresolved (B6a step 4)
438fa10 SymbolicDynamics: a truncated search decides nothing and lends no wall; tests for B6a steps 2-3 (B6a step 3)
acc31de Walks: member lobes holding different markers raise MembersDisagree; the class is unresolved before any walk (B6a step 2)
131f39f DualWalk: a landing span keeps the edges it overlaps, or is span_meets_no_node; never a walk from every edge (B6a step 1)
```

Each commit passes the two unit modules (80, 80, 82 tests).

## `git diff --stat 76d0710`

```
 src/tanglepack/topology/DualWalk.py           |  89 ++++++++++++---------
 src/tanglepack/topology/SymbolicDynamics.py   |  51 +++++++++---
 tests/unit/topology/test_dual_walk.py         |  21 ++++-
 tests/unit/topology/test_symbolic_dynamics.py | 110 +++++++++++++++++++++++---
 4 files changed, 212 insertions(+), 59 deletions(-)
```

Source: DualWalk +53/-36 lines and SymbolicDynamics +42/-9, 140 changed lines in all, mostly docstrings, so under the plan's "well under 150" only narrowly. The tests account for the other 131 changed lines (rubric items 7 to 10).

## Suite

`MPLBACKEND=Agg env/bin/python -m pytest -q -rxX` on the branch:

```
XFAIL tests/invariants/test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing] - C3: INV-11: det J < 0 is out of scope (C3); FixedPoint.set_k_value rejects a saddle with exactly one negative eigenvalue
1160 passed, 38 skipped, 5 deselected, 1 xfailed in 201.53s (0:03:21)
```

The only xfail is `cases.KNOWN_ISSUES`' `orientation_reversing`. The base count would be 1157 passed: there are three net new tests (7 is replaced, 8, 9 and 10 are new).

## What changed

1. **Step 1, `DualWalk._nodes_meeting`: two tiers.**
   - First tier: the nodes the span overlaps by more than its end tolerance.
   - Second tier, if the first is empty: the nodes it overlaps at all (`> 0`).
   - If both are empty: `[]` and one WARNING ("(resolution flag): the class is unresolved").
   - `shortest_walks` then returns `WalkSearch(start, goal, "span_meets_no_node")`.
   - The new status is in the `WalkStatus` Literal, the `WalkSearch.status` docstring, the `shortest_walks` Returns, and the Dev Notes bullet (1).
   - The overlap is kept as a list of `(node, overlap)` pairs, not a dict, so test doubles need not be hashable.
2. **Step 2, `MembersDisagree(ValueError)` in `DualWalk.py`.**
   - `lobe_image_markers` raises it instead of the WARNING-and-read-first behaviour. Its message names both bridges, the class and the markers in conflict. Its Returns and Raises sections are updated.
   - `_resolve_itinerary` computes `wanted` before `shortest_walks`. It catches the exception and sets the reason "members hold different markers: the class is not one homotopy class (resolution; <detail>)". It logs the existing "class %s unresolved: %s" WARNING and returns without walking.
   - `MembersDisagree` is NOT exported from `tanglepack.topology.__init__`, because that file is outside the plan's list. It is importable from `tanglepack.topology.DualWalk` and re-imported by `SymbolicDynamics`. If it should be public, a follow-up must add it to `__init__`.
3. **Step 3, truncation in `_resolve_itinerary`.**
   - After the status check: `elif search.truncated: cd.unresolved_reason = "search budget reached before uniqueness was established"`. Both branches share the one WARNING, now worded `class %s unresolved: walk from %s to %s: <reason>`, where it used to end `is <status>`. The reason for a failed status is still `"dual-graph walk <status>"`.
   - The wall gate in `symbolic_dynamics` also requires `cd.itinerary is not None`, so an unresolved (truncated) class lends no wall.
   - `shortest_walks` statuses and the default `max_walks=64` are unchanged.
4. **Step 4, docs.**
   - `is_reliable` keeps its formula. Its docstring now says that `itinerary is None` covers a truncated search, a span that meets no node, and members that disagree.
   - Also updated: the `ClassDynamics.unresolved_reason` docstring, the `DualWalk` Dev Notes truncation bullet, `WalkSearch.truncated`, and the `SymbolicDynamics` Dev Notes non-crossing bullet.

### Tests

These are in the plan's two test files. I wrote them so that the suite stays green and the rubric is checkable; the test engineer may extend them.

- `test_dual_walk.py`
  - `test_a_span_meeting_no_edge_falls_back_to_every_edge` is **removed**. It is replaced on the same graph by `test_a_span_meeting_no_edge_leaves_the_walk_undecided`: status `span_meets_no_node`, no walks, one WARNING.
  - New: `test_a_span_narrower_than_its_error_starts_from_the_edge_it_lies_in`. It uses the `two_starts` graph with the scaled span `(1.2, 1.205)`; tol 0.012 is greater than its width. The walk starts from `[f3]` only and is `unique`.
- `test_symbolic_dynamics.py`
  - New: `test_members_holding_different_markers_leave_the_class_unwalked`. A SimpleNamespace dual makes the REAL `lobe_image_markers` read two lobes with different verdicts: `enclosed` returns `{m: True}` and then `{m: False}`. It asserts that `_resolve_itinerary` never calls `shortest_walks` (spy count 0), that the itinerary and search are None, that the reason contains "different markers", and that there is one WARNING.
  - New: `test_a_truncated_search_decides_nothing_and_is_no_wall`. It runs the real driver over a `FakeDual` tie at `max_walks=1`, with two classes on different unstable branches. Class 1 has `search.truncated` and a walk, but `itinerary is None`, a reason containing "budget", and `is_reliable` False. Class 2 is walked with `walls == []`, and with those walls it sees both walks of its tie. The `chords_cross` stub makes any wall reject the `through_1` walk.
  - Mutation check: I dropped the new `cd.itinerary is not None` gate and the test failed on `walls_b == []`. The gate is restored.
  - `_drive` gained `dual_graph=` and `max_walks=`. With `walks=None` it does not stub `_resolve_itinerary`. Existing callers are unchanged.
  - A test pins a reason by substring only ("different markers", "budget"); no log wording is pinned, per test policy.
- **On base**, with both test files copied into the base worktree and then restored with `git checkout`: the 4 new or rewritten tests FAIL and the other 78 pass. No import breaks on base, because `MembersDisagree` is reached as `SD.MembersDisagree`, never imported at module level. Failures:
  - span narrower: start faces `[f0, f3]`, not `[f3]`;
  - no edge: `'unique' == 'span_meets_no_node'`;
  - members: `AttributeError` (the spy returned None, so the class WAS walked);
  - truncated: the itinerary is not None.

## Behaviour changes

On every case and series step checked: **none in rules, `is_reliable` or closed-table size.** That is all 11 `LAW_CASES`, `p3_6_blasts`, `p3_16`, `build_k28(blasts=8)`, and every step of all five `SERIES`: k32 0 to 10, k28 0 to 8, p3 0 to 9, inversion 0 to 4, nested_outer 0 to 11.

There are two observable differences, both expected:

- **k32_series blasts 6 and 8: the "meets none of its" WARNING no longer fires.** Base fires it 2 and 8 times. Every one of those spans is now kept by the second tier, and the per-class itineraries, words and reasons are identical. This is a log change only.
- **`heteroclinic` at `max_walks=2`:**
  - base: `is_reliable=True`, with `bh` `('unique', truncated=True)` and no reason;
  - branch: `is_reliable=False`, with `bh` unresolved and the reason "search budget reached before uniqueness was established";
  - at 64, both are True.

On hand-built graphs, as the plan states: a span that overlaps nothing now gives `span_meets_no_node` (it used to give `unique` from every node), and a WARNING wording changed (see step 3).

### Rubric 4 / 6 comparison

The scripts are mine, adapted from the analyst's `probe.py`, `sim.py`, `k32.py` and `cmp.py` so that they run against a ROOT worktree; each script asserts its imports come from ROOT. They are in `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b6a_check/`: `rubric4.py`, `rubric6.py`, `run.sh`, and per-worktree output in `base/` and `branch/`.
- `rubric4.py` prints, per target and step: `is_reliable`, `len(sd.table)` (closed), the class count, and a hash plus the full list of letter-free per-class rules. Each rule is (class label, itinerary labels, symbols as (class label, direction), unresolved reason).
- It also prints the count of "meets none" WARNINGs.

`diff -r base branch`, verbatim:

```
diff -r base/rubric6.txt branch/rubric6.txt
1,2c1,2
< max_walks=2: is_reliable=True
<     bh ('unique', True) | None
---
> max_walks=2: is_reliable=False
>     bh ('unique', True) | search budget reached before uniqueness was established
diff -r base/series:k32_series.txt branch/series:k32_series.txt
49c49
< series:k32_series blast 6: reliable=True closed=18 classes=18 rules=4e98cda7f7e4 span_flags=2
---
> series:k32_series blast 6: reliable=True closed=18 classes=18 rules=4e98cda7f7e4 span_flags=0
87c87
< series:k32_series blast 8: reliable=True closed=45 classes=45 rules=4b63537bf15b span_flags=8
---
> series:k32_series blast 8: reliable=True closed=45 classes=45 rules=4b63537bf15b span_flags=0
```

## Rubric

1. **Suite passes; `-rxX` lists only KNOWN_ISSUES.** PASS: 1160 passed, 38 skipped, 1 xfailed (`orientation_reversing`, C3).
2. **New and rewritten tests fail on base and pass on the branch.** PASS: all 4 fail on base (failures listed above) and pass on the branch.
3. **The diff stat lists exactly the four files.** PASS (see above).
4. **Rules, `is_reliable` and closed-table size are equal on base and branch.** PASS for every law case, `p3_6_blasts`, `p3_16`, `k28_8`, and every k32 step (6 and 8 included), plus the other four series. The only diff lines are the span-flag counter and rubric 6.
5. **The two `k32_series` constant-while-holes-constant tests pass.** PASS: run by name, 2 passed; they are also in the full suite.
6. **`heteroclinic` at `max_walks=2`.** PASS: True on base, False on the branch; the one unresolved class is `bh`, with the budget reason.
7. **Old fallback test gone; replacement asserts `span_meets_no_node`, no walks, one WARNING.** PASS.
8. **Narrow scaled span starts from one node.** PASS: `start_faces == [f3]` on the branch, `[f0, f3]` on base.
9. **Members disagree.** PASS: `itinerary is None`, reason contains "different markers", spy counts zero `shortest_walks` calls. The graph is a hand-built stub; the real `lobe_image_markers` raises.
10. **Two-walk tie at `max_walks=1`.** PASS: unresolved with the budget reason; the later class on a disjoint branch gets no wall and sees both walks.
11. **Greps.** PASS:
    - `grep -n "may leave from any of them" src` prints nothing;
    - `grep -n "the non-crossing rule reads" src` prints nothing;
    - `grep -n "span_meets_no_node" DualWalk.py` finds the Literal (`:566`) and the return (`:1123`), plus docs at `:68`, `:590` and `:1098`.

## Notes for the test engineer and reviewer

- `MembersDisagree` is not in `tanglepack.topology.__all__`, because `__init__.py` is outside scope (see step 2).
- Nothing outside the four files was needed.
