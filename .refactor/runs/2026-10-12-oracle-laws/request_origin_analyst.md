# Request: B0 — `oracle-laws` (tests only)

*From the author, 2026-10-12 ("Go ahead and run B0 through the planner and coder agents"). Input for one `tangle_refactor` run. The first slice of `docs/alignment_plan_2026-10-09.md`.*

## The question

**Can the suite check the three PDF definitions (strong pip, pseudoneighbor pair, hole containment) and the one textbook result (k=10 is the full 2-shift) on real builds by computations independent of the code under test?** Today no law does; about a third of the 98 laws can only fail if the code reports its own failure (`docs/code_fidelity_audit_2026-10-09.md`, section 5).

## What "done" means

Exactly slice B0 of the alignment plan (`docs/alignment_plan_2026-10-09.md`, section B, "B0 — `oracle-laws`"): its scope, five steps, behaviour ("none in the library"), not-doing list and six rubric items are the plan. In short:

1. `check_chosen_pip_is_strong` — the open-box test re-derived from the registry and `FixedPoint` alone.
2. `check_reference_pairs_clear_of_X` — both open intervals of every reference pair, scaling only, law cases only; failing cases on base become strict `KNOWN_ISSUES` entries citing NUM-18 with `wp="B2"` (the orchestrator's probe of 2026-10-09 found zero punctures on every law case, so the expected set is empty; the slice verifies on base).
3. `check_propagated_hole_inside_container` — the reference bridge's ends carried back with the real `map_inv`, projected onto the container's polyline.
4. `check_k10_is_the_full_shift` — letter-free: one active class, word two occurrences of itself with opposite directions, class matrix `[[2]]`, refined matrix `[[1,1],[1,1]]`; `NOT_APPLICABLE` on every other law case (decision Q11: the single exception to the no-pin policy).
5. `build_inversion_period2(stable_steps=10)` (decision Q12).

Files: `tests/helpers/laws.py`, `tests/cases.py`, `tests/invariants/test_law_partition.py`, `tests/invariants/test_law_symbolic.py`. No `src/` file.

## Where to start

`docs/alignment_plan_2026-10-09.md` (B0, and A.4-A.6 for the definitions), `docs/code_fidelity_audit_2026-10-09.md` (section 5, the test-suite findings), `tests/helpers/laws.py`, `tests/cases.py`, `tests/conftest.py` (the test policy).

## Binding context

- The author's decisions of 2026-10-09 on the alignment plan are in `.refactor/decisions.md` (Q1-Q13; Q11 and Q12 gate this slice and are "yes").
- `CLAUDE.md` invariants are physics. Topology never calls the map; a TEST may (step 3 uses `map_inv` on the test side).
- Test policy (`tests/conftest.py` docstring, 2026-10-05): rule-based laws, `KNOWN_ISSUES` strict xfail citing an issue id and the slice that fixes it, `NOT_APPLICABLE` visible skips; the k=10 pin is the one authorised exception.
- The analyst stage is skipped for this run: the alignment plan's B0 section and the audit are the findings (probes of 2026-10-09 recorded in both).
- The working tree may carry uncommitted documentation changes (`docs/`, `.refactor/decisions.md`); the change branches from `HEAD` and never touches them.
