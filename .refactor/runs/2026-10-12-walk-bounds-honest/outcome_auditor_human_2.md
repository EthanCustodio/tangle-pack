# Outcome: B6a `walk-bounds-honest`

*2026-10-12, auditor. Branch `refactor/walk-bounds-honest` @ `8a1cdc3` (5 commits), base `76d0710`. Review passed in round 1. Not merged; merging is your call. Your main checkout still has uncommitted changes in `docs/`, `.refactor/` and `artifacts/`; the branch does not touch them.*

## The answer now

**Question:** where does the walk return a word it has not established?

**Answer:** nowhere now. There were three such places, and each one now leaves the class **unresolved** with a stated reason. Unresolved means `itinerary is None`, so `is_reliable` turns False through its unchanged formula.

| Situation | What it means mathematically | Before | Now |
|---|---|---|---|
| (i) The `max_walks` budget ran out after a walk was found | Uniqueness of the shortest admissible walk was never checked | `unique` + `truncated`, counted as reliable | Unresolved: "search budget reached before uniqueness was established" |
| (ii) The landing span `f(X)` overlaps no stable edge of `X'` | The walk has no edge it is known to leave from | Walked from *every* edge, with a WARNING | Status `span_meets_no_node`, unresolved |
| (iii) Two member lobes disagree on which marker images `f(b)+f(σ_b)` must hold | The members are not one homotopy class | WARNING, then read the first lobe | `MembersDisagree` is raised and caught before any walk; unresolved |

One addition goes with these: an unresolved class no longer acts as a wall for later classes. Before, a truncated walk could still constrain the walks that came after it.

## The most important hunk (`DualWalk._nodes_meeting`)

Before:
```python
meeting = [node for node in nodes if min(hi, node.hi_cdist) - max(lo, node.lo_cdist) > tol]
if meeting:
    return meeting
logger.warning("... meets none of its %d stable node(s); the walk may leave from any of them", ...)
return nodes
```
After:
```python
overlaps = [(node, min(hi, node.hi_cdist) - max(lo, node.lo_cdist)) for node in nodes]
meeting = [n for n, o in overlaps if o > tol] or [n for n, o in overlaps if o > 0]
if not meeting:
    logger.warning("... meets none of its %d stable node(s) (resolution flag): the class is unresolved", ...)
return meeting      # shortest_walks: empty -> WalkSearch(start, goal, "span_meets_no_node")
```

The tolerance `SCALING_RTOL·|c|` is meant to drop a *neighbouring* edge that a scaled end's error just grazes. On k32 at blasts 6 and 8 the span is narrower than that tolerance: about 0.001 to 0.006 wide, against a tolerance of about 0.012. So it "met" nothing, even though it lies inside one edge.

The second tier ("overlaps at all") keeps that one edge. A plain `return []`, as the alignment plan wrote it, would have made 10 k32 classes unresolved and dropped the closed alphabet at blast 8 from 45 to 42.

The other two sites are small additions to `SymbolicDynamics._resolve_itinerary` (`SymbolicDynamics.py:1183-1216`): a `try/except MembersDisagree` before the walk, and `elif search.truncated:` after it.

## What got shorter or clearer

The source diff is +95 / −45 lines, and most of the additions are docstrings and Dev Notes. Nothing got shorter. What got clearer:
- **The fallback is gone.** It was a silent "walk from anywhere" path. `grep "may leave from any of them" src` and `grep "the non-crossing rule reads" src` both return nothing.
- **Every unresolved walk has a reason.** All of them now share one WARNING line, which states the reason instead of a raw status.
- **The docs match the code.** `is_reliable`, `WalkSearch.truncated`, `ClassDynamics.unresolved_reason` and both modules' Dev Notes now say what the code does.

The tests add +136 / −14 lines.

## Behaviour changes, all of them

1. **None on any law case, deep run (`p3_6_blasts`, `p3_16`, k28 x8) or series step.** The words, `is_reliable` and the closed-table size are identical on base and branch.
2. **`heteroclinic` at `max_walks=2`:** reliable goes from True to False. The one unresolved class is `B:p3@2.0/L#0 <-> B:p3@2.0/L#2`, with the budget reason.
3. **`heteroclinic` at `max_walks=1`:** two more classes become unresolved (`R#0<->R#4` and `R#0<->B:L#0`). The table was already unreliable at this budget on base. Budgets of 8 and 64 are unchanged.
4. **k32 series:** the "meets none of its" WARNING fired 10 times on base (blasts 6 and 8) and fires 0 times on the branch. The words are identical, so this is a log change only.
5. **Log wording:** the walk-path WARNING now reads `class X unresolved: walk from A to B: <reason>`. It used to read `... is <status>`.
6. **Hand-built graphs:** a span that overlaps nothing now gives `span_meets_no_node`. On base it gave `unique` from every edge.

## What was verified, and how

- **Full suite, which I ran on the branch:** `1162 passed, 38 skipped, 5 deselected, 1 xfailed`. The only xfail is KNOWN_ISSUES `orientation_reversing` (C3).
- **Scope, which I checked myself:** `git diff --stat 76d0710` lists exactly the four planned files, and the worktree is clean. I read the whole diff.
- **`heteroclinic` budgets, which I ran myself** (`scratchpad/audit/r6.py`, in both worktrees):
  - base: 1 → False, 2 → True, 8 → True, 64 → True;
  - branch: 1 → False, 2 → False, 8 → True, 64 → True;
  - the closed size is 75 at every budget.
- **k32 series, which I ran myself** (`scratchpad/audit/k32.py`): the class and closed counts are identical at every blast from 0 to 10 (…18, 18, 45, 45, 45), with 0 ambiguous. The span WARNINGs were 10 on base and 0 on the branch.
- **Checked by the test engineer and the reviewer; I did not re-run these:**
  - the full base-vs-branch dump of every rule over all law cases, the deep runs and all five series (`scratchpad/te/pin.py`);
  - 6 new or rewritten tests that fail on base and pass on the branch;
  - 8 mutation checks, each caught. The one that matters most: removing the second tier fails the real k32 regression test with 42 == 45.

## What did not go to plan

**The plan was followed as approved.** It already carried the analyst's four corrections to the alignment-plan section:
- the second tier;
- raising before the walk instead of setting a flag after it;
- the wall gate;
- the `max_walks=2` probe.

**The test engineer added one test the plan did not list:** `test_a_budget_too_small_..._unresolved[1|2]`. It runs on the real heteroclinic case and takes about 1 s.

**Sites (ii)'s empty path and (iii) never fire on real data.** They are guarded only by hand-built unit tests.

## Where I am least sure this is what you meant

1. **`MembersDisagree` is not exported from `tanglepack.topology`**, even though `lobe_image_markers`, which raises it, is. The reason is that `__init__.py` was outside the four files. A caller has to import it from `tanglepack.topology.DualWalk`. It needs a one-line follow-up if you want it public.
2. **The second tier is a judgement about the tolerance**, not a rule from your decisions. A span that overlaps an edge by less than its own error is taken to lie in that edge.
   - The alternative reading is "the span is too uncertain to know, so the class is unresolved". That reading would cost k32 its 45-class alphabet at blast 8.
   - Whether `SCALING_RTOL·|c|` is simply too coarse for such spans is left as a numerics question.
3. **`CLAUDE.md` is now stale on this point.** It still says "a span meeting no node falls back to all of them with a WARNING", and its `WalkSearch.status` list lacks `span_meets_no_node`. The changed WARNING wording is pinned by no test, but anyone grepping logs for `is truncated` will no longer find it.

*Knowledge files updated:* `.refactor/knowledge/topology.md` and `.refactor/knowledge/testing.md` (B6a sections). The outcome entry is appended to `.refactor/decisions.md`.
