# Request: B2 — `pseudoneighbor-both-intervals`

*From the author, 2026-10-14 ("continue"). Input for one `tangle_refactor` run. The seventh slice of `docs/alignment_plan_2026-10-09.md`; B0, B1, B6a, B8, B4a and B4b are merged (main `ba23a1f`; runs archived under `.refactor/runs/2026-10-12-*`, `2026-10-13-*` and `2026-10-14-holes-are-combinatorial/`).*

## The question

**Does the code test both open intervals of a candidate pseudoneighbor pair against X, as the definition does?** The PDF definition tests BOTH open intervals (the unstable one and the stable one) against `X`; the code tests the unstable one against all iterates and the stable one against registered crossings only (`Pseudoneighbor.py:329-470` at the plan's numbering; Dev Notes "all the work is in the unstable-interval check"). Author decision Q1 (S-a, 2026-10-09): test the open STABLE interval against X exactly as the unstable one is tested. One function for both intervals finishes the definition; what it kills at depth is a discovery, not a failure.

## What "done" means

Exactly slice B2 of the alignment plan (`docs/alignment_plan_2026-10-09.md`, section B, "B2 — `pseudoneighbor-both-intervals`"): its scope, five steps, behaviour ("none on the corpus"), not-doing list and five rubric items are the plan. In short:

1. `_interval_is_clear(trellis, x0, x1, stability, candidates, collision_rtol)` generalises `_unstable_interval_is_clear`: for `"unstable"` exactly today's body; for `"stable"` the arc is the stable branch `ix0.manifold_b_key` (asserted equal for both), the owner is ITS fixed point, the cycle `branch_position_map("stable")`, a forward step scales the stable cdist by `1/β` (`per_step_beta("unstable")`), the candidates are the crossings whose STABLE side is that fixed point (D7 mirrored: for a heteroclinic pair on B's branch, A×B and B×B crossings), and the same three readings in the same order: table landing by id, `r` read against the pair's registered images in its own frame, β-scaling with the 2-D collision guard. The closed-form `m`-range uses `β^{-k}` per lap (the stable cdist shrinks along forward iterates): the same code with `lap = β^{-k}` and the inequalities unchanged.
2. `compute_pseudoneighbors` keeps a consecutive pair iff BOTH calls return True; the `lo <= 0` rejection (anchor) applies to both.
3. Dev Notes: delete the paragraphs "the open stable interval between a candidate pair is empty by construction" and the NUM-18 asymmetry; state the definition once.
4. Unit tests (the only place the clause bites today): a hand-built trellis whose window lacks the representative of a deeper registered crossing, whose scaled backward image lands inside the stable interval of a consecutive pair → no pair; the same with a table link → no pair, by id; the pair's own orbit (2-D collision) → still a pair. Each fails on base.
5. `stress/numerics/e13_both_intervals.py`: for every series of `cases.SERIES` and the E12 depths, the reference set per blast before/after, each lost pair with `reference_witness`, and the seed/closed/holes series; its output is the report the author reads before DEC-25 (BCG-8) is decided.

Files (four): `src/tanglepack/topology/Pseudoneighbor.py`, `tests/unit/topology/test_pseudoneighbor.py`, `tests/cases.py` (B0's `KNOWN_ISSUES` entries for `reference_pairs_clear_of_X` removed, if any exist: B0's outcome found the stable clause vacuous on every law case, so the expected set is empty and this may be a no-op), and the read-only re-measurement script `stress/numerics/e13_both_intervals.py` (series to depth, report only).

**Behaviour.** None on the corpus: B0's law `check_reference_pairs_clear_of_X` found the stable clause vacuous on every law case, k28 x8 and p3 x6 (A.5: on the pip branch consecutiveness implies it when `N` is complete inside the cut), so every fingerprint is unchanged. The clause bites only where a window representative is missing (an unresolved detection under Q13 since B1, a dedup): the hand-built unit tests, and possibly the deep series, which e13 re-measures and reports. No tolerance added (`collision_rtol` reused as is).

**Not doing.** DEC-25 itself (the author decides on the e13 report); any change to the window or to trajectory extension.

**Rubric** (the plan's five): suite passes and the diff touches only the four files, and the law `reference_pairs_clear_of_X` passes on every law case with no `KNOWN_ISSUES` entry; fingerprints of every law case equal on base and branch; the three unit tests fail on base; `grep -n "_unstable_interval_is_clear" src` is empty, one function with one `stability` parameter; the e13 report exists and lists, per series step, the references lost or gained with witnesses.

## Where to start

`docs/alignment_plan_2026-10-09.md` (B2; A.5 for the pair definition; E Q1), `docs/code_fidelity_audit_2026-10-09.md` (the NUM-18 finding), `src/tanglepack/topology/Pseudoneighbor.py` (`_unstable_interval_is_clear`, `compute_pseudoneighbors`, `extend_pseudoneighbor_trajectories`, the Dev Notes), `src/tanglepack/topology/StrongPip.py` (the sibling interval test, for the shape of the stable reading), `tests/helpers/laws.py` (`check_reference_pairs_clear_of_X`, B0), `tests/unit/topology/test_pseudoneighbor.py`, `tests/cases.py` (`KNOWN_ISSUES`, `SERIES`), `stress/numerics/e12_reverify.py` (the model for e13: depths, report layout).

## Binding context

- The author's decisions on the alignment plan are in `.refactor/decisions.md` and plan section E; Q1 is decided S-a. D7 (a heteroclinic pair on B's branch is tested on A's cycle, k and β for the unstable side) stands and is mirrored for the stable side. DEC-25 (the inversion series' drop at blast 6) is NOT decided here; e13's report is its input.
- `CLAUDE.md` invariants are physics: topology never calls the map (table first, `per_step_beta` scaling only when no chain; never the map); the iterate table is exact; a 2-D collision (both cdists) decides "same point", never the product. Since B1, `N` is the registry's `resolved_ids`, so "registered crossings" means crossings of N.
- Test policy (`tests/conftest.py` docstring, 2026-10-05): rule-based laws, no case-specific pins; rubric 2 is a base-vs-branch fingerprint equality.
- The analyst stage is NOT skipped: the analyst verifies `_unstable_interval_is_clear` and the stable-side check at `HEAD` (current line numbers), confirms by running the B0 law and a probe that the stable clause is vacuous on every law case and on the series steps the plan names, checks whether any `KNOWN_ISSUES` entry for `reference_pairs_clear_of_X` exists in `tests/cases.py`, reads `e12_reverify.py` for the depths and layout e13 should copy, and reports what differs from the plan section. It does not re-plan the slice.
- The working tree may carry uncommitted documentation changes (`docs/`, `.refactor/`); the change branches from `HEAD` and never touches them.
