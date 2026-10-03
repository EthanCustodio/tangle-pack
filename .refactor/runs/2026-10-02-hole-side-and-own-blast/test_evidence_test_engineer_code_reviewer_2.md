# Test evidence: hole-side-and-own-blast

*Test engineer, 2026-10-02. Plan: `artifacts/plan_planner_human_2.md`. Coder report: `artifacts/change_coder_test_engineer.t1.md`.*

- **Base:** `2315204`, detached worktree `/home/dezhu/Development/TanglePack/tangle-pack-hole-side-and-own-blast-base` (left clean).
- **Branch:** `refactor/hole-side-and-own-blast`, worktree `/home/dezhu/Development/TanglePack/tangle-pack-hole-side-and-own-blast`, now at `3b33baa`. That is the coder's 5 commits plus my one test commit.
- **Interpreter:** `../tangle-pack/env/bin/python`. The worktrees have no `env/`. `pyproject.toml` sets `pythonpath = ["src", "tests"]`, so pytest in each worktree imports that worktree's own `src`. For my probe scripts I checked this with `sys.path.insert(0, "src")`: `tanglepack.__file__` resolved to the base worktree's `src`.

## Summary

| Run | Command (in the worktree) | Result |
|---|---|---|
| Base, full suite | `MPLBACKEND=Agg ../tangle-pack/env/bin/python -m pytest -q -p no:cacheprovider` | **782 passed, 1 skipped** |
| Base, rubric 2 pinning files | same, with `tests/test_stable_partition_period3.py tests/test_higher_period_cartoon.py tests/test_loom_blast_restore.py` | **24 passed** |
| Branch, full suite, coder's commits only | same | **786 passed, 1 skipped** |
| Branch, full suite, with my tests (`3b33baa`) | same | **795 passed, 1 skipped** (782 + 4 from the coder + 9 from me) |
| Branch test files run against base code | see "At base" below | **25 passed, 5 failed**. The 5 failures are exactly the bug-fix tests. |

`git diff --stat 2315204 --name-only` lists only the five planned files: `CLAUDE.md`, `henon_cases.py`, `StablePartition.py`, `test_higher_period_cartoon.py`, `test_stable_partition_period3.py`. My tests went into the two existing test modules, so rubric 5 still holds.

## New tests I added (commit `3b33baa`)

The code changes are A (propagated-hole side read on its image sub-arc), A' (direct-hole side from crossing signs) and B (each nested zone blasts only its own fixed point). In the regression cases the plan expects A and A' to change nothing, and B to leave the default nested words unchanged.

**Probe first.** I took an id-free signature of every hole: `(iterate, bridge_side, sorted (which, row) of openings)`, as a Counter. It was identical at base and on the branch for:

- k28 one blast
- k28 two blasts
- `build_period3()`
- `build_period3(blasts=4)`
- both tangles of `build_nested()`

The only probe differences were in nested bridge counts and blast sizes: period-3 bridges go from 31 to 27, and `blast_sizes` from `[32, 7]` to `[3, 3]`.

**Words at base and branch.** Printed per letter, these are also identical for p3, p3 with 4 blasts, and nested:

- active: `a->b, b->c, c->a u^-1 w^-1`, plus nested `d->d uu^-1 e^-1, e->f, f->d uu^-1 d^-1`
- inert words also identical
- `is_reliable` is True in all three runs

New tests:

`tests/test_stable_partition_period3.py`:

- `test_p3_hole_sides_are_pinned` (pin, step A and A'). Checks every hole of `build_period3()`: iterates -3..2, one left and one right each, with the exact opening rows. Ids are dropped because they do not reproduce between builds.
- `test_k28_two_blast_hole_sides_are_pinned` (pin, step A and A'). The four k=2.8 two-blast holes: iterates -3..0, all `left`, with their openings.
- `test_direct_hole_side_is_the_side_of_its_coordinates[k10_partitioned | k28_partitioned | k28_two_blasts_partitioned | p3_partitioned]` (pin, step A'). Every directly punched hole (`hole.pair is not None`) has the `bridge_side` that the old measurement `_bridge_side_of(trellis, bridge_for_pair(...), hole.coords)` gives. This means the crossing-sign rule reproduces the old side hole for hole wherever the old rule was right. `p3_partitioned` covers both fixed points of the conftest nested session.

`tests/test_higher_period_cartoon.py`:

- `test_p3_words_survive_four_blasts` (pin, rubric 9). `build_period3(blasts=4)` gives active words `{a: b, b: c, c: a u^-1 w^-1}`.
- `test_nested_default_words_are_pinned` (pin, step B and rubric 9). `build_nested()` gives the p3 words plus `{d: d uu^-1 e^-1, e: f, f: d uu^-1 d^-1}`.
- `test_nested_outer_blasts_leave_the_inner_bridges_alone` (bug fix, step B). The inner tangle of `build_nested()` has as many bridges as `build_period3()` alone.

Together the new tests run in about 1.5 s.

## At base: the branch's test files run against `2315204`

I copied both branch test files into the base worktree as `tests/test_zz_branch_period3.py` and `tests/test_zz_branch_cartoon.py` and ran them. I removed both copies afterwards; `git status` in the base worktree is clean.

```
MPLBACKEND=Agg ../tangle-pack/env/bin/python -m pytest -q -p no:cacheprovider tests/test_zz_branch_period3.py tests/test_zz_branch_cartoon.py
5 failed, 25 passed in 4.56s
```

All 25 passes include every pinning test above: the 8 new pins and all the pre-existing tests in both files. The 5 failures are exactly the bug-fix tests, each failing for the reason the plan predicts:

| Test | Failure at base |
|---|---|
| `test_deep_p3_holes_share_bridge_side_through_symbolic_dynamics[6_blasts]` | `AssertionError: Holes of origin (41, 42) disagree on bridge_side: iterate 0 is right, iterate -3 is left; ... -5 ...; ... -6` (rubric 3: origin matches) |
| `...[15_steps]` | Same assertion at iterates -3/-5/-6, but origin **`(3, 17)`**. The plan says `(12, 33)` and the coder saw `(21, 35)`. Registry ids are not reproducible between builds (conftest says so), so only the assertion and the iterates are stable. The plan's id is not. Rubric 4 holds in substance. |
| `test_nested_inner_words_are_the_period3_words` | `check_holes_share_bridge_side` on origin `(51, 52)`, iterates -3/-5/-6, with "blast: skipping bridge that failed to iterate" WARNINGs |
| `test_nested_blast_order_does_not_matter` | same origin `(51, 52)` assertion |
| `test_nested_outer_blasts_leave_the_inner_bridges_alone` | `assert 31 == 27` |

All five pass on the branch.

## Break-it checks (on the branch, each restored with `git checkout -- <file>`)

| # | Break | Suite run | Result |
|---|---|---|---|
| B1 | A' sign flipped: `outward = -1 if second.stable_cdist > first.stable_cdist else 1` (`StablePartition.py:355`) | full | **20 failed, 219 errors**, 556 passed. Failures include all four `test_direct_hole_side_is_the_side_of_its_coordinates[...]`, `test_p3_hole_sides_are_pinned`, `test_p3_words_survive_four_blasts`, both deep tests, and both coder nested tests. `test_k28_two_blast_hole_sides_are_pinned` and `test_nested_default_words_are_pinned` ERROR in their fixtures on `check_holes_share_bridge_side` (direct and propagated holes disagree). |
| B2 | A undone: `_punch_in_bridge` stops passing `span` (`:1547`) | full | **2 failed**, 793 passed: `[6_blasts]` and `[15_steps]`. As the plan expects, on today's passing runs A changes nothing, so only the deep tests guard it. `test_nested_inner_words_are_the_period3_words` (4 inner blasts) still passes. |
| B3 | Span slice half-implemented: `poly = poly[lo:]`, ignoring the span's far end (`:823`) | the three hole/cartoon files | **2 failed** (`[6_blasts]`, `[15_steps]`), 76 passed |
| B4 | B reverted: `fixed_point=fixed_points` in `build_nested` (`henon_cases.py:264`) | full | **3 failed**, 792 passed: `test_nested_inner_words_are_the_period3_words`, `test_nested_blast_order_does_not_matter` (the "skipping bridge" WARNING is caught), `test_nested_outer_blasts_leave_the_inner_bridges_alone` |
| B5 | Nested default recipe changed: `outer_blasts: int = 1` (`henon_cases.py:220`) | `test_higher_period_cartoon.py` | **2 failed**: `test_nested_default_words_are_pinned`, `test_nested_resolves_everything_with_real_walks` |

After the checks, `git status --short` on the branch is clean at `3b33baa`.

## Rubric items I checked directly

1. Pass. 795 passed, 1 skipped (at least 782 + 1 skipped + new tests).
2. Pass. At base, the three pinning files give 24 passed. On the branch they pass in the full suite.
3. Pass. `[6_blasts]` fails at base on origin `(41, 42)` and passes on the branch.
4. Pass in substance. `[15_steps]` fails at base on the same assertion and iterates, but the origin id is `(3, 17)` in my run, not `(12, 33)`, because ids are build-dependent. It passes on the branch.
5. Pass. Five files.
9. Pass. The k=10/k=2.8 fixture word tests pass; p3 with 4 blasts and the nested default words are now pinned and identical at base and on the branch.
   - Unstated behaviour change: `build_nested().blast_sizes` goes from `[32, 7]` to `[3, 3]`, as the coder reported. Same cause as 31 -> 27. I did not pin it.
10. Pass. The coder's `test_nested_blast_order_does_not_matter` asserts that no "skipping bridge" warning is logged. B4 shows the assertion bites.

## Not tested, and why

- **The 16-step period-3 case.** It is not pinned (the plan names only 15 steps). The coder reports it now runs. I did not add it, because it duplicates the 15-step mechanism.
- **Reliability of the deep runs.** `is_reliable` is False there by design (Not doing 1). The deep tests check only that the run completes.
- **A' on the deep runs where the old measurement was wrong** (the +1/+2 lobes). There is no ground truth to pin beyond the invariant. That case is covered indirectly: `[6_blasts]`/`[15_steps]` pass `check_holes_share_bridge_side` on the branch.
- **`_hole_openings` whole-polyline nearest vertex** (Not doing 3). It is latent and has no failing case, so it is not tested.
