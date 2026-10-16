# Request: B6a — `walk-bounds-honest`

*From the author, 2026-10-12 ("go ahead and archive and run the request, then I'll do the command again"). Input for one `tangle_refactor` run. The third slice of `docs/alignment_plan_2026-10-09.md`; B0 (`oracle-laws`) and B1 (`resolved-crossings`) are merged (main `76d0710`; runs archived in `.refactor/runs/2026-10-12-oracle-laws/` and `.refactor/runs/2026-10-12-resolved-crossings/`).*

## The question

**Where does the dual-graph walk return a word it has not established?** Three places (plan B6a, audit of A.8c): (i) a `max_walks` budget that runs out AFTER a walk was found leaves status `unique` with `truncated=True`, and `is_reliable` reads only `ambiguous` / `unresolved` (`SymbolicDynamics.py:839-848`; `max_walks=1` makes `heteroclinic` "reliable"); (ii) a landing span meeting no stable node "may leave from any of them" (`DualWalk.py:729-757`); (iii) two member lobes holding different marker images are a WARNING and the first is read (`DualWalk.py:815-828`). Under A.8c each is `unresolved` with its reason. This is W-a's decision-free half (Q4): it changes no definition, it stops claiming what is unproven.

## What "done" means

Exactly slice B6a of the alignment plan (`docs/alignment_plan_2026-10-09.md`, section B, "B6a — `walk-bounds-honest`"): its scope, four steps, behaviour ("none on any law case, deep run or series step"), not-doing list and five rubric items are the plan. In short:

1. `_nodes_meeting` returns `[]` when nothing meets; `shortest_walks` returns `WalkSearch(start, goal, "span_meets_no_node")` (one WARNING kept, worded as a resolution flag).
2. `lobe_image_markers` returns `None` and sets a flag the search records (`members_disagree`) when two member lobes hold different marker images; `_resolve_itinerary` makes the class unresolved with reason "members hold different markers: the class is not one homotopy class (resolution)".
3. A `WalkSearch` with `truncated=True` is never `unique`: `_resolve_itinerary` treats `search.truncated` as unresolved ("search budget reached before uniqueness was established"). The default `max_walks=64` is kept (nothing truncates on any case or series; the audit's `max_walks=8` changes nothing).
4. `is_reliable` docstring: unchanged formula (it reads `itinerary is None`), now covering the three.

Files (four): `src/tanglepack/topology/DualWalk.py`, `src/tanglepack/topology/SymbolicDynamics.py`, `tests/unit/topology/test_dual_walk.py`, `tests/unit/topology/test_symbolic_dynamics.py`.

**Behaviour.** None on any law case, deep run or series step: none of the three fires there (the WARNING strings "meets none of its", "hold different markers" and the truncated path appear in no recorded log; rubric 3). On a hand-built graph the statuses change as stated.

**Not doing.** Statuses are not renamed (B6b collapses the enumeration); the evidence comparison is not moved (B6b); an UNPLACED marker image (`image_marker` None) stays informational; no change to `MAX_DETOUR`.

**Rubric** (the plan's five): suite passes on the branch and the diff touches only the four files; `symbolic_dynamics(..., max_walks=1)` on `heteroclinic` gives `is_reliable` False on the branch and True on base; letter-free `rules` of every law case and of `p3_6_blasts`, `p3_16`, `k28` at 8 blasts equal on base and branch; three new unit tests on hand-built dual graphs (span meeting no node, two members disagreeing on a marker, a budget of 1 with two candidate walks), each `unresolved` on the branch and failing on base; `grep -n "may leave from any of them" src` is empty.

## Where to start

`docs/alignment_plan_2026-10-09.md` (B6a; A.8c for the word of a class; E Q4), `docs/code_fidelity_audit_2026-10-09.md` (the walk findings), `src/tanglepack/topology/DualWalk.py` (`_nodes_meeting`, `shortest_walks`, `lobe_image_markers`, `WalkSearch`), `src/tanglepack/topology/SymbolicDynamics.py` (`_resolve_itinerary`, `is_reliable`), `tests/unit/topology/test_dual_walk.py`, `tests/unit/topology/test_symbolic_dynamics.py`.

## Binding context

- The author's decisions on the alignment plan are in `.refactor/decisions.md` and plan section E; Q4 is decided W-a (the walk is the word, the registered chain is an oracle), and this slice is its decision-free half.
- `CLAUDE.md` invariants are physics. Topology never calls the map (lookup only); every itinerary has even length and splits into disjoint pairs; `is_reliable` is the one reliability flag and provisional / resolution flags never feed it except through `itinerary is None`.
- Test policy (`tests/conftest.py` docstring, 2026-10-05): rule-based laws, no case-specific pins; rubric 3 is an equality between base and branch, not a pinned word.
- Since B1 (`76d0710`), `N` is the registry's `resolved_ids`; k28 at 8 blasts is 4 / 4 / 4 (holes / seed / closed) on the current main, not the 24 / 15 / 23 the plan quotes. Rubric 3's equality is between base and branch at the SAME commit lineage, so this changes nothing for B6a; the analyst notes it when it records the line numbers above, which were taken before B1 and may have moved.
- The analyst stage is NOT skipped: the analyst verifies the three sites against the code at `HEAD` (current line numbers, that none of the three fires on any law case, deep run or series step, and that `max_walks=8` changes nothing), and reports what differs from the plan section. It does not re-plan the slice.
- The working tree may carry uncommitted documentation changes (`docs/`, `.refactor/`); the change branches from `HEAD` and never touches them.
