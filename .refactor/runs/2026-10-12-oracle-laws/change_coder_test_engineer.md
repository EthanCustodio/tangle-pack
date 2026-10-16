# Change: B0 `oracle-laws` (tests only), coder report

*Coder, 2026-10-12. Job `change`, first round plus the author's addendum (round 2). Plan: `artifacts/plan_planner_human_2.md`, approved in
`artifacts/plan_ok_planner_human.json` (Q1 answered yes: the conftest sentence is in, so five files; the addendum adds `tests/regression/test_open_p3_deep_runs.py`, six in all).*

## Round 2: the deep-run law test (author addendum, 2026-10-12)

**What happened.** The first pass stopped with 5 failures. `test_open_p3_deep_runs.py::test_every_law_holds_on_the_deep_runs`
runs every `LAYERS` law on five builds that are not law cases: `p3_15`, `p3_16`, `p3_6_blasts`,
`p3_16_12` and `inversion_period2_10`. On those builds the k10 law returns 0 by design, so the test
flagged it as vacuous.

**The fix.** The author widened the file list by exactly that file (addendum in
`plan_ok_planner_human.json`), and commit `1029fcc` implements the rule:

- `tests/helpers/laws.py` declares `LAW_CASES_ONLY = ("reference_pairs_clear_of_X", "k10_is_the_full_shift")`
  beside `LAYERS`. Its comment gives the reason: the scaling-only oracle and the textbook pin run on law
  cases only (A.5, plan section 6).
- The deep-run test's loop now runs every check of every layer except those two, through
  `law_tier.assert_non_vacuous`. Its module docstring bullet and its own docstring say so.
- There are no `NOT_APPLICABLE` entries for the deep builds.

**One side effect.** The loop no longer goes through `run_layer`, which had no hook for leaving a law
out. As a result:

- A failure is reported at the first failing check, not collected per layer.
- `run_layer` used to skip `NOT_APPLICABLE` / `KNOWN_ISSUES` pairs, and the loop no longer looks at
  them. That has no effect: none exist for these five builds.
- Every check still has to be non-vacuous.

**Result.** The deep-run file passes: 14 passed. The full suite passes with 0 failures and 38 skips
(numbers in "Suite result" below).

## Workspace

- **Base commit:** `a3bd8f06056e0f9357bf388500edf83c82393b3e` (`a3bd8f0`, the author's `HEAD` on `main`)
- **Branch:** `refactor/oracle-laws`
- **Worktree:** `/home/dezhu/Development/TanglePack/tangle-pack-oracle-laws`
- **The author's tree:** untouched, apart from this artifact. Its uncommitted `docs/`,
  `.refactor/decisions.md` and `artifacts/` files are as they were.
- **Base worktree:** I did not create `../tangle-pack-oracle-laws-base`. No `src/` file changed, so the
  branch's `src/` is the base's. That worktree is for the reviewer's rubric-2 run.

**Commits:**

```
1029fcc Laws: LAW_CASES_ONLY keeps the scaling-only pair oracle and the k=10 pin off the deep runs
70bf44b Cases: the module docstring names inversion_period2's 10 stable steps
ff14ef8 Docs: the oracle laws in the layer docstrings; the Q11 exception in the test policy (B0 step 6)
49a9561 Cases: inversion_period2 grows 10 stable steps (B0 step 5, Q12)
b9d562b Laws: k=10 is the full 2-shift, the one textbook pin (B0 step 4, Q11)
4657284 Laws: the pip box, both pair intervals and hole containment as oracles (B0 steps 1-3)
```

**`git diff --stat a3bd8f0`:** six files, as the addendum allows.

```
 tests/cases.py                             |  24 ++--
 tests/conftest.py                          |   3 +
 tests/helpers/laws.py                      | 169 +++++++++++++++++++++++++++++
 tests/invariants/test_law_partition.py     |   3 +
 tests/invariants/test_law_symbolic.py      |   3 +-
 tests/regression/test_open_p3_deep_runs.py |  15 ++-
 6 files changed, 201 insertions(+), 16 deletions(-)
```

**How the tests ran:** from the worktree, as
`MPLBACKEND=Agg PYTHONPATH=src /home/dezhu/Development/TanglePack/tangle-pack/env/bin/python -m pytest -q -rxXs`.
I confirmed that `tanglepack.__file__` resolves to the worktree's `src/`.

## What was done

**`tests/helpers/laws.py`.** Four checks, plus five small private readers. They are registered in
`LAYERS`: three appended to `"partition"` and one to `"symbolic"`.

- **`check_chosen_pip_is_strong`.** Written as the plan's code.
  - It reads the whole registry. It skips anything that is an anchor, not homoclinic to the pip's fixed
    point, or not on its stable cycle.
  - Each crossing is carried onto the pip's stable branch by `steps = (pos(r) - pos(q0)) mod k` backward
    steps: stable cdist divided by `per_step_beta("stable")**steps`, unstable cdist divided by
    `per_step_beta("unstable")**steps`.
  - q0's own orbit is skipped by a collision on both cdists within `SCALING_RTOL`.
  - It asserts that nothing falls in the open box. It never uses the iterate table.
- **`check_reference_pairs_clear_of_X`.** For every reference pair and both stabilities:
  - It takes the interval's own key, `branch_cycle` and `per_step_beta` (D7).
  - It tests every non-anchor crossing of `trellis.own_intersection_ids` whose key of that stability
    lies on the cycle.
  - `_laps` lists every `n ≡ n0 (mod k)` with `lo < c·β^n < hi`, solved between two logarithms.
  - A landing passes if `_table_image` (from `registry.iterate_table`: the direct entry, else composed
    ±1 links) is an end of the pair. It also passes if `_scaled(r, n)` collides with `x0` or `x1`: the
    keys come from `advance_key`, each cdist is scaled by its own side's β, and both must match within
    `SCALING_RTOL`.
  - Any other landing raises.
  - The count is 1 per pair and stability.
- **`check_propagated_hole_inside_container`.** For every hole with `iterate < 0`:
  - The container is `trellis.bridge_between(*hole.bounding_ids)`.
  - The ends are `hole.origin`. Each is carried back `|iterate|` times by
    `dynamical_system.map_inv_batch`, the REAL inverse map.
  - Each carried point is projected onto the container's node polyline (nearest segment, cdist
    interpolated). It must lie strictly between the container's two crossing cdists.
  - The [changed] closure: an end that `_table_image` links to one of the container's two crossings
    counts as that shared end and is not projected.
  - **Small deviation from the plan text, stricter in spirit (kept as implemented, per the orchestrator).** The plan says "if the table links the
    end, the linked crossing MUST be a container crossing". I wrote it as "linked to a container
    crossing → closed; otherwise the strict projection test applies". So an end linked to some other
    crossing is still tested geometrically rather than failing outright. The probe shows no such case.
- **`check_k10_is_the_full_shift`.** Returns 0 unless `case.name == "k10"`. On k10 it makes the plan's
  six letter-free assertions and returns the number of active tokens (2).

**`tests/cases.py`:**

- `NOT_APPLICABLE` gains `("k10_is_the_full_shift", name)` for the 10 law cases other than k10, and adds
  `"inversion_period2"` to the `every_active_class_walks` tuple.
- `build_inversion_period2` now takes `stable_steps: int = 10`. Its docstring says "9 unstable and 10
  stable steps … (author Q12)", and "past the pip orbit" is gone.
- The module docstring's `inversion_period2` bullet now reads "9 unstable / 10 stable … (at 10 / 10 …)"
  instead of the stale 9/9. The parenthetical about the 10/10 regression case is kept.
- `KNOWN_ISSUES` is unchanged: check 2 fails on no case on base (see rubric 11).

**Docstrings:** the plan's sentences in `test_law_partition.py` and `test_law_symbolic.py`, and one
sentence on the Q11 exception in `tests/conftest.py`'s no-pin paragraph.

## Suite result

`MPLBACKEND=Agg PYTHONPATH=src env/bin/python -m pytest -q -rxXs` on `1029fcc`:

`1117 passed, 38 skipped, 5 deselected, 1 xfailed in 156.10s`

- **Failures:** none.
- **XFAIL:** only `test_orientation_reversing_case_builds[orientation_reversing]`, which is the base's.
- **Skips:** 38, which is 11 more than base. Ten are `not applicable: k10_is_the_full_shift on <case>`
  and one is `not applicable: every_active_class_walks on inversion_period2`.
- **Round 1, for the record:** `5 failed, 1112 passed, ...`, all five in the deep-run test, fixed by
  `1029fcc`.

## Behaviour changes

- **Library:** none. `git diff a3bd8f0 -- src` is empty.
- **The deep-run regression test:** no longer runs `reference_pairs_clear_of_X` or
  `k10_is_the_full_shift` on its five builds. Every other law still runs there, each still required to be
  non-vacuous.
- **The suite:**
  - Four new laws inside the existing layer tests.
  - 11 new visible set-aside skips.
  - No new xfail.
  - The law case `inversion_period2` is now the 9-unstable / 10-stable build: alphabet 6 classes,
    10 holes, `is_reliable` True, 0 unresolved classes. This was measured on the branch; no test pins it.

## Break-it checks

These were scratch scripts in my session scratchpad, not committed (plan §6).

- **Rubric 7 (check 1).**
  - Patched `StrongPip.is_strong_pip` to return the original result with `is_strong_pip=True`.
  - On k10 the unpatched function rejects crossings `[5, 6]`.
  - For each, I built k10, ran `set_strong_pip(fp, pip)` and then `topology()`. The trellis is cut at
    that pip, and `check_chosen_pip_is_strong` **RAISES** on both: pip 5 on crossing 2 (s 74.15 < 76.54,
    u 74.15 < 651.8), and pip 6 on crossing 2 as well.
- **Rubric 8 (check 2).**
  - Patched `Pseudoneighbor._unstable_interval_is_clear` to return True and built `k28_two_blasts`.
  - `check_reference_pairs_clear_of_X` **RAISES**: `pair (3, 10): M^1(3) lies inside its open unstable
    interval`.
- **Rubric 9 (check 3).**
  - Patched `StablePartition._containing_bridge` to return the bridge after the real container
    (bisect, then `ids[i+1], ids[i+2]`), and `_map_backward` to return None. Then built `k10`.
  - `check_propagated_hole_inside_container` **RAISES** on its first hole: origin (5, 6), end 5,
    iterate -1, container (4, 3). The carried cdist is 83.94, not in (84.37, 631.4).

## Rubric, my reading

1. **Holds** (round 2): 0 failures. The only xfail is the base's
   `orientation_reversing_case_builds`, and no check-2 entry was needed.
2. **Partly checked.** `git diff a3bd8f0 -- src` is empty, so the library under the new laws is the base
   library. I ran all four laws on all 11 law cases on base `src`: they pass. Under the 9/9
   `inversion_period2` build, check 3 passes with the closure and check 4 is not applicable. I did not
   run the reviewer's "base worktree with the branch's `tests/`" run.
3. **Holds.** The diff stat lists exactly the six files: the plan's five plus the deep-run test the addendum allows.
4. **Holds.** The only changed build is `inversion_period2`: 6 classes, 10 holes, `is_reliable` True.
5. **Holds.** I grepped the bodies of the four new checks and the five helpers. `StrongPip`,
   `is_strong_pip`, `Pseudoneighbor`, `_unstable_interval_is_clear`, `StablePartition`,
   `_containing_bridge`, `bridge_for_pair` and `image_cdist` each have 0 hits.
   - Check 1 has no `iterate`.
   - Checks 2 and 3 contain "iterate" only as `hole.iterate` and in prose about the "iterate table".
   - The capitalised names are absent. Lowercase `pseudoneighbor` appears in a docstring and in
     `reference_pseudoneighbors`.
6. **Holds.** The three checks are listed under `"partition"` and the one under `"symbolic"`. The test
   ids are still per case, so they are unchanged.
7. **Holds** (RAISES; see above).
8. **Holds** (RAISES; see above).
9. **Holds** (RAISES; see above).
10. **Holds.** With the closure clause turned off (scratch toggle), check 3 fails on base for
    `k28_one_blast`: hole (0, 3), end 5, iterate -1, cdist 42.250102 against hi 42.250098. It also
    fails on `inversion` (end 3, iterate -2, 0.88968 against 0.88903), and on `inversion_period2` at
    9/9. It passes on the other cases. With the closure it passes on all 11. At 10 stable steps
    `inversion_period2` passes even without the closure.
11. **Holds.** Check 2 passes on base on all 11 law cases, so there are no `KNOWN_ISSUES` entries.
    Landings tested per case:

    | Case | Landings |
    |---|---|
    | `k10` | 0 |
    | `k28_one_blast` | 8 |
    | `k28_two_blasts` | 0 |
    | `p3` | 34 |
    | `nested` | 34 |
    | `nested_inner_blasts` | 40 |
    | `inversion` | 15 |
    | `inversion_one_blast` | 0 |
    | `inversion_two_blasts` | 2 |
    | `inversion_period2` | 11 (28 at 10 steps) |
    | `heteroclinic` | 22 |

    Every landing found was its pair's own end: the table linked it, or it collided within
    `SCALING_RTOL`. The stable half found no puncture anywhere, which matches A.5's argument.
12. **Holds.** `check_k10_is_the_full_shift` passes on k10, and `-rs` lists exactly the 10
    `not applicable: k10_is_the_full_shift on <case>` skips plus
    `every_active_class_walks on inversion_period2`.
13. **Holds.** The signature reads `stable_steps: int = 10`, and the function's docstring no longer
    contains "past the pip orbit".

## Non-vacuity counts on base

Each item is the count of things checked.

| Case | Pip box | Pair intervals | Hole containment |
|---|---|---|---|
| `k10` | 7 | 2 | 4 |
| `k28_one_blast` | 5 | 4 | 4 |
| `k28_two_blasts` | 9 | 2 | 6 |
| `p3` | 27 | 4 | 12 |
| `nested` | 36 | 6 | 18 |
| `nested_inner_blasts` | 40 | 6 | 18 |
| `inversion` | 12 | 4 | 8 |
| `inversion_one_blast` | 18 | 2 | 12 |
| `inversion_two_blasts` | 24 | 4 | 12 |
| `inversion_period2` (9/9) | 6 | 4 | 16 |
| `heteroclinic` | 108 | 8 | 138 |
