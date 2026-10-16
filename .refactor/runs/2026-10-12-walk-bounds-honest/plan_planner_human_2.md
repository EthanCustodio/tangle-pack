# Plan: B6a `walk-bounds-honest`

*2026-10-12, planner, first plan. Base: `main` at `76d0710`. Source: section B6a of `docs/alignment_plan_2026-10-09.md`, corrected by the analyst's findings (`artifacts/findings_analyst_planner.t1.md`).*

## 1. The answer

There are three places where the walk claims a word it has not established:

- **(i) A truncated search.** If the `max_walks` budget runs out after a walk is found, the status stays `unique` (or `trivial`), and `is_reliable` never reads `truncated` (`SymbolicDynamics.py:839-848`, `:1180`).
- **(ii) A landing span that meets no stable node.** The search then "may leave from any of them" (`DualWalk.py:729-757`).
- **(iii) Two member lobes that hold different marker images.** This logs a WARNING and reads the first lobe (`DualWalk.py:815-828`).

The line numbers are the same at `HEAD` as in the plan section. (i) and (iii) never fire on any case or series step, at a budget of 64 or 8.

The section got two facts wrong, and the analyst ran both to check:

1. **Site (ii) does fire.** It fires on the `k32_series` at blasts 6 and 8: 10 searches. Each time it is a tolerance artifact. The landing span is 0.001 to 0.006 wide and lies wholly inside one node's edge. But the test asks for an overlap greater than `SCALING_RTOL * |c|`, which is about 0.012, so the span "meets" nothing. The fallback then happens to pick the right node.
   - If step 1 is done as written (return `[]`), 10 classes become unresolved and the closed alphabet at blast 8 drops from 45 to 42.
   - That would almost certainly break `test_closed_alphabet_constant_while_holes_constant[k32_series]`. The analyst computed the compared quantities but did not run the test under the patch.
2. **Rubric 2's probe is already False on base.** At `max_walks=1`, `heteroclinic` has two classes that run out of budget with no walk. The probe that is True on base and False after the change is `max_walks=2`.

## 2. Recommendation

Make the change. It is section B6a with the analyst's four corrections folded in:

- a two-tier `_nodes_meeting`;
- the marker disagreement is detected before walking;
- a truncated class is no wall;
- the rubric 2 probe uses `max_walks=2`.

Each of these keeps the section's purpose (stop claiming what is unproven), and none widens it.

## 3. Scope

**Slug:** `walk-bounds-honest`. The branch is `refactor/walk-bounds-honest` in `../tangle-pack-walk-bounds-honest`, made from `76d0710`. The base worktree is `../tangle-pack-walk-bounds-honest-base`.

**Files (four, unchanged from the section):**
- `src/tanglepack/topology/DualWalk.py`
- `src/tanglepack/topology/SymbolicDynamics.py`
- `tests/unit/topology/test_dual_walk.py`
- `tests/unit/topology/test_symbolic_dynamics.py`

None of the corrections needs another file. The `k32_series` regression test already guards step 1, and it is read, not edited. The change should come to well under 150 changed lines. Your uncommitted `docs/` and `.refactor/` edits are not touched.

## 4. Steps

### 1. A span meets the nodes it really overlaps; if it overlaps none, the class is unresolved (site ii)

`_nodes_meeting` works in two tiers. The tolerance is there to drop a neighbour that a scaled end's error barely touches. It should not reject the one edge a span lies in when the span is narrower than its own error.

```python
overlap = {node: min(hi, node.hi_cdist) - max(lo, node.lo_cdist) for node in nodes}
meeting = [n for n in nodes if overlap[n] > tol] or [n for n in nodes if overlap[n] > 0]
if not meeting:
    logger.warning("the landing span [%.6g, %.6g] of %s meets none of its %d stable node(s) "
                   "(resolution flag): the class is unresolved", lo, hi, ref.label, len(nodes))
return meeting
```

- In `shortest_walks`, if either list comes back empty, return `WalkSearch(start, goal, "span_meets_no_node")`.
- Add `"span_meets_no_node"` to the `WalkStatus` Literal and to the `WalkSearch.status` docstring. This is a new status, not a rename.
- `_resolve_itinerary` already treats any other status as unresolved (reason `"dual-graph walk span_meets_no_node"`), so `SymbolicDynamics.py` needs no change for this step.
- Rewrite the `_nodes_meeting` docstring to match.

### 2. Members that disagree on a marker leave the class unresolved, and the class is not walked (site iii)

- `lobe_image_markers` raises `MembersDisagree` (a small `ValueError` subclass in `DualWalk.py` whose message names the two bridges and the markers in conflict) instead of logging a WARNING and returning the first lobe.
- `_resolve_itinerary` computes `wanted` before calling `shortest_walks`, catches the exception, and returns unresolved with the reason "members hold different markers: the class is not one homotopy class (resolution)". It logs through the existing "class %s unresolved: %s" WARNING.
- The disagreement must be caught before walking. Otherwise a `None` would run the unconstrained search and pick a word under fewer constraints.
- This replaces the section's `members_disagree` flag on `WalkSearch`. That flag would be set only after a search that should never run.
- Update the `Returns:` of `lobe_image_markers` and add a `Raises:`.

### 3. A truncated search decides nothing (site i)

In `_resolve_itinerary`, after the status check:

```python
if walk is None or search.status not in ("unique", "trivial", "ambiguous"):
    cd.unresolved_reason = f"dual-graph walk {search.status}"
elif search.truncated:
    cd.unresolved_reason = "search budget reached before uniqueness was established"
```

- Both branches log the one "class unresolved" WARNING and return.
- This covers a truncated `trivial` too, so it reads "never decided" rather than only "never unique".
- `shortest_walks` keeps its statuses, and the default `max_walks=64` is kept.
- The wall gate (`SymbolicDynamics.py:1410`) also requires `cd.itinerary is not None`. Otherwise a truncated class that is now unresolved would still lend its walk to later classes as a wall.

### 4. Docs

- `is_reliable`: same formula, and the docstring now says that `itinerary is None` covers a truncated search, a span that meets no node, and members that disagree.
- Fix the text that describes the old behaviour:
  - the `DualWalk.py` module Dev Notes (`:91-99`);
  - `WalkSearch.truncated` (`:590-592`);
  - the Dev Notes in `SymbolicDynamics.py:94-99`.

## 5. Behaviour

**None on any law case, deep run or series step.**
- (i) and (iii) never fire there, at a budget of 64 or 8.
- Every firing of (ii) (k32 at blasts 6 and 8) overlaps exactly one node, and the second tier keeps that node. The analyst simulated this and got identical per-class itineraries at every k32 step from 0 to 10.

**On hand-built graphs and at small budgets, the statuses change as stated:**
- `heteroclinic` at `max_walks=2` stops being reliable.
- A span that overlaps nothing gives `span_meets_no_node`. Before, it gave `unique` from every node.

## 6. Not doing

- Statuses are not renamed. B6b collapses the list.
- The evidence comparison is not moved (B6b).
- An unplaced marker image (`image_marker` None) stays informational.
- No change to `MAX_DETOUR` or to `SCALING_RTOL`. Whether the landing span's error bound is too coarse next to spans a tenth of its size is a numerics question for later.
- `shortest_walks` does not change status on truncation. The decision stays in `_resolve_itinerary`.

## 7. Rubric

1. `MPLBACKEND=Agg env/bin/python -m pytest -q` passes on the branch, and `-rxX` lists only `cases.KNOWN_ISSUES`.
2. Every new or rewritten test in rubric items 7 to 10 fails on the base commit and passes on the branch. The reviewer checks this by copying the two test files into the base worktree.
3. `git diff --stat 76d0710..refactor/walk-bounds-honest` lists exactly the four files.
4. Letter-free per-class `rules`, `is_reliable` and the closed-table size are equal on base and branch for:
   - every `cases.LAW_CASES` case;
   - `p3_6_blasts`, `p3_16`, and `build_k28(blasts=8)`;
   - **every `k32_series` step, 6 and 8 included.**

   This is checked by one script run in both worktrees, and the outputs are diffed. That is "no behaviour change beyond section 5".
5. `test_closed_alphabet_constant_while_holes_constant[k32_series]` and `test_class_count_constant_while_holes_constant[k32_series]` pass on the branch.
6. `symbolic_dynamics(topo.dual_graph, topo.classes, max_walks=2, letter_for=topo.letter_for).is_reliable` on `heteroclinic` is True on base and False on the branch. The one unresolved class has the reason "search budget reached before uniqueness was established".
7. `test_a_span_meeting_no_edge_falls_back_to_every_edge` (`test_dual_walk.py:333`) no longer exists. In its place, a test on the same graph asserts status `"span_meets_no_node"`, no walks, and one WARNING.
8. A new `test_dual_walk.py` test: a scaled-end span narrower than its tolerance and lying inside one of two nodes starts from that node only (`start_faces` has one face). On base it starts from both.
9. A new `test_symbolic_dynamics.py` test: two member lobes holding different marker images (hand-built graph, or `enclosed` monkeypatched) give `itinerary is None` with the members-disagree reason, and `shortest_walks` is never called (a spy counts zero calls).
10. A new `test_symbolic_dynamics.py` test: the two-walk tie at `max_walks=1` gives a class with `itinerary is None` and the budget reason. That class also adds no wall: a later class on disjoint branches still sees both of its walks.
11. `grep -n "may leave from any of them" src` prints nothing. `grep -n "the non-crossing rule reads" src` prints nothing. `grep -n "span_meets_no_node" src/tanglepack/topology/DualWalk.py` finds both the Literal and the return.

## 8. Questions

1. **Raise or return for site (iii)?** By default, `lobe_image_markers` raises `MembersDisagree` and `_resolve_itinerary` catches it. The alternative is a sentinel return value. I chose the exception because `None` already means "no member lobe". "Yes" keeps the exception.
