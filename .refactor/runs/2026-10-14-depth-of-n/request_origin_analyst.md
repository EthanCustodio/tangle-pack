# Request: B3 — `depth-of-N`

*From the author, 2026-10-14 ("continue"). Input for one `tangle_refactor` run. The eighth slice of `docs/alignment_plan_2026-10-09.md`; B0, B1, B6a, B8, B4a, B4b and B2 are merged (main `9cd29d5`; runs archived under `.refactor/runs/2026-10-1{2,3,4}-*/`).*

## The question

**Who fixes the depth `m` of the trellis of the pip, `T_m(q0)`?** Today `m` is "the deepest registered forward iterate of the pip", and two mechanisms can freeze it under blasting while N keeps deepening elsewhere: (i) the blast's `min_separation` (passed as 1e-5 by every recipe) can drop the child that carries `f^{m+1}(q0)` and never feed it again (the plan's B3 diagnosis of BCG-5: k32 blast 9, `f^10(q0)` dropped, the orbit stalls); (ii) since B1 (`resolved-crossings`), the deepest detected orbit point can itself be an unresolved detection (ratio under 3) and so leave N and never be linked (B1's outcome: k28 x8 id 245 = `f^9(q0)` at ratio 2.991, k32 x9 id 523 = `f^10(q0)` at ratio 1.48, nested_outer x9 id 447 at 1.78), after which every later blast leaves `T(q0)` where it was. Author decision Q2 (N-a, 2026-10-09): a child holding a pip-orbit point is never dropped; `m` = the deepest registered iterate. No explicit `depth=` parameter (N-b not chosen).

## What "done" means

Exactly slice B3 of the alignment plan (`docs/alignment_plan_2026-10-09.md`, section B, "B3 — `depth-of-N`"): its scope, four steps, behaviour, not-doing list and four rubric items are the plan. In short:

1. In the blast's child filter (`Blast.py:397-412` at the plan's numbering), before the `min_separation` test: a child is `pip_bearing` iff one of its end crossings is a forward iterate of a selected fixed point's chosen pip (`session.chosen_strong_pip(fp)` → the registry's forward iterate chain; the point `f^{m+1}(q0)` is born as an END of a child of the bridge ending at `f^m(q0)`); a pip-bearing child is kept regardless of separation (one INFO naming it).
2. `BlastResult.pip_depth: dict[FixedPoint, int]` = `forward_depth(pip)` after the blast (data).
3. `cases.blast_series` records `pip_depth` per step; regression law `test_pip_orbit_deepens_under_blasting(name)`: `pip_depth` is non-decreasing and increases by one on every step whose frontier contained the bridge ending at `f^m(q0)` (read off the result).
4. Docstring of `blast_zone`: `min_separation` is a cost guard on children that re-trace existing curve; it never decides the depth of `T(q0)` (A.4).

Files (four): `src/tanglepack/loom/Blast.py`, `tests/unit/loom/test_blast_budget.py` (or a new `test_blast_pip_orbit.py`), `tests/regression/test_class_count_under_blasting.py`, `tests/cases.py`.

**Behaviour.** Law cases: none (no law case blasts past a stall; k28 x1/x2 and inversion x1/x2 keep their frontiers). k32 series: the plan says `forward_depth(pip)` on base reads `[1, 2, …, 9, 10, 10]` over blasts 0..10 with the stall at blast 10, and on the branch the orbit keeps deepening (11 at blast 10), so N and the hole/class series at blast 10 may change (a discovery; the numbers go in the report beside E12's); nested_outer at blast 10 likewise. THOSE NUMBERS PREDATE B1. Since B1 the k32 stall is at blast 9 and is caused by mechanism (ii), an unresolved `f^10(q0)`, which N-a's "never drop the child" does not fix (the child may be kept and the point still unresolved). The analyst must re-measure `forward_depth(pip)` per step on every series at `HEAD`, separate the two mechanisms, and say what step 1 alone achieves; the planner then states the behaviour honestly and puts to the author whether mechanism (ii) is B3's (the remedy B1 named is "refine"; the alternative B1's auditor raised is placing an unresolved pip-orbit point by `image_cdist`, as PK-6 places a missing one, so the cut keeps its depth) or a later slice's.

**Not doing.** N-b (`depth=` parameter): not chosen; the sharp D1 law stays a measured series law.

**Rubric** (the plan's four, to be re-derived against `HEAD`): suite passes and the diff touches only the four files; on `k32_series`, `pip_depth` over blasts 0..10 deepens on the branch where it plateaus on base (the plan's `[1, …, 10, 11]` vs `[1, …, 10, 10]` must be re-measured; the law fails on base at the plateau); fingerprints of every law case equal on base and branch; a unit test: a child whose end is `f^{m+1}(q0)` with interior points under `min_separation` of an existing curve is KEPT, the same child without the pip end is dropped, failing on base.

## Where to start

`docs/alignment_plan_2026-10-09.md` (B3; A.4 the trellis of the pip; E Q2, Q8), `.refactor/runs/2026-10-12-resolved-crossings/outcome_auditor_human_2.md` (section 7 point 1: the pip-orbit freeze after B1) and that run's `test_evidence_test_engineer_code_reviewer_2.md` (the three instances), `src/tanglepack/loom/Blast.py` (`blast_zone`, the child filter, `BlastResult`), `src/tanglepack/numerics/BridgeIterator.py` (`iterate_bridge`, where a child's ends are registered), `src/tanglepack/topology/Trellis.py` (`cut_at_pip_orbit`, `image_cdist`, the PK-6 placement of a missing orbit point), `src/tanglepack/loom/TangleSession.py` (`chosen_strong_pip`), `tests/cases.py` (`blast_series`, `SERIES`, `_blast`), `tests/regression/test_class_count_under_blasting.py`, `tests/unit/loom/test_blast_budget.py`, `stress/numerics/results/e13_both_intervals.jsonl` (per-step `n_trellis_crossings` on every series: where N stops growing).

## Binding context

- The author's decisions on the alignment plan are in `.refactor/decisions.md` and plan section E; Q2 is decided N-a; Q8 restates D1 (at fixed depth `m` the hole set and class count are non-decreasing under blasting; a deeper `m` is a discovery). B1's revise-round decision (2026-10-12): the rule stands, an unresolved pip-orbit point leaves N, the remedy is "refine".
- `CLAUDE.md` invariants are physics: the trellis of the pip is cut at the deepest REGISTERED (and, since B1, resolved) orbit point; `image_cdist` places a missing orbit point by scaling and flags it (PK-6); blasting iterates the bridges whose innermost zone it is; the point budget is a resource guard, not an invariant. Topology never calls the map.
- Test policy (`tests/conftest.py` docstring, 2026-10-05): rule-based laws, no case-specific pins; the series laws are measured laws over `cases.SERIES`.
- The analyst stage is NOT skipped: the analyst verifies the child filter at `HEAD` (current line numbers), measures `forward_depth(pip)` per step on all five series (and `nested_inner`, `heteroclinic` if cheap) at the recipe cutoff, attributes every plateau to mechanism (i) or (ii) with the dropped child or the unresolved detection named, says what step 1 alone changes (which plateaus it removes), and reports what differs from the plan section. It does not re-plan the slice.
- The working tree may carry uncommitted documentation changes (`docs/`, `.refactor/`); the change branches from `HEAD` and never touches them.
