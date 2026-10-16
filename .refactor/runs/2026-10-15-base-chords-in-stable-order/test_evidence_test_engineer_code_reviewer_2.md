# Test evidence: B7b `base-chords-in-stable-order` (test engineer, traversal 1)

**To the code reviewer:** this is the evidence you grade. Read it against the plan's rubric (`artifacts/plan_planner_human_2.md` §7) and grade **"pass"** or **"revise"**. Everything below I re-measured myself. I did not reuse the coder's outputs: the probe was rerun into a fresh directory, and both suites were rerun.

## Where and how

| | |
|---|---|
| Base | `a6a4f68` (main), worktree `/home/dezhu/Development/TanglePack/tangle-pack-base-chords-in-stable-order-base`, detached. Clean before and after; files were copied in only for the run under rubric 7-9 and then restored with `git checkout -- .` |
| Branch | `refactor/base-chords-in-stable-order` @ `63477e3`, worktree `/home/dezhu/Development/TanglePack/tangle-pack-base-chords-in-stable-order`. Clean; every break was restored with `git checkout` |
| Command | From inside each worktree: `MPLBACKEND=Agg PYTHONPATH=src ../tangle-pack/env/bin/python -m pytest -q -rxX -p no:cacheprovider` |
| Scratch | `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/te/`: `suite_*.txt`, `ids_*.txt`, `out/` (probe), `check.py`, `base_copied.txt`, `branch_new.txt`, `bc_*.txt`, `k28.py` |

**No test was added or changed by me, and nothing was committed.** Every rubric item already has a test from the coder. One teeth limitation is reported under "Break-it checks".

## Rubric 1: branch suite

```
XFAIL tests/invariants/test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing] - C3: INV-11: det J < 0 is out of scope (C3); ...
1211 passed, 47 skipped, 5 deselected, 1 xfailed in 218.11s (0:03:38)
```

The only xfail is C3, and there is no XPASS. **Met.**

## Rubric 2: base suite

```
XFAIL tests/invariants/test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing] - C3 ...
XFAIL tests/regression/test_class_count_under_blasting.py::test_members_share_the_P1_itinerary_under_blasting[k32_series] - B7b: B7: at 6 and 7 blasts ...
1200 passed, 46 skipped, 5 deselected, 2 xfailed in 218.10s (0:03:38)
```

**Collected node ids, base vs branch** (`pytest --co`): 1248 on base and 1259 on the branch. The `diff` shows exactly 11 added ids and none removed:
- `test_base_chords_are_inside_the_lobe_image_under_blasting[k28_series|p3_series|inversion_series|k32_series|nested_outer_series]`
- `test_base_chords_are_inside_the_lobe_image[k10|inversion|inversion_period2|heteroclinic|p3_6_blasts]`
- `test_iterated_cut.py::test_base_crossings_pair_in_stable_order`

The accounting closes:
- 1200 passed on base, + 10 new passing tests, + 1 former xfail (k32 T) now passing = 1211.
- 46 skipped on base, + 1 (`p3_6_blasts`) = 47.

Every test the branch leaves unchanged passes on base. The base's second xfail is the k32 T case whose known issue the branch deletes.

One existing test lost an assertion: `test_foreign_crossings_pair_in_stable_order` dropped `(12,"chord") in empty[11]`. The plan directs this (§4 step 7). I checked that the assertion was toothless: every `OWN` crossing there has the default `stable_cdist=1.0` (`_crossing`, `test_iterated_cut.py:40`). Python's sort is stable, so on a tie it keeps chain order, and the assertion would hold under either rule. It was not removed to make the suite pass. **Met.**

## Rubric 3: diff stat

```
 src/tanglepack/topology/PartitionFamily.py         |  89 +++++++-------
 src/tanglepack/topology/SymbolicDynamics.py        |   7 +-
 tests/cases.py                                     |  19 +--
 tests/helpers/laws.py                              | 133 ++++++++++++++++++++-
 .../regression/test_class_count_under_blasting.py  |  18 ++-
 tests/regression/test_foreign_chords.py            |  70 ++++++-----
 tests/unit/topology/test_iterated_cut.py           |  24 ++--
 7 files changed, 261 insertions(+), 99 deletions(-)
```

These are exactly the seven files of §3, with nothing under `docs/`, `.refactor/`, `CLAUDE.md` or `artifacts/`. **Met.** The plan estimated about 200 changed lines; the diff is 360. Most of the excess is the helper (133 lines) and the docstrings.

## Rubric 4: base-vs-branch probe over 63 items

**How I ran it.** I reran the coder's `probe2.py` (copied unchanged to `te/probe2.py`) through `te/run.sh`. That is `run2.sh` with its output redirected to `te/out/`. It made 42 processes: 16 builds and 5 series in each worktree, with no errors. Every row's `cases.__file__` resolved to its own worktree:

```
worktrees: {'tangle-pack-base-chords-in-stable-order-base'} {'tangle-pack-base-chords-in-stable-order'}
```

**How I compared.** I wrote my own comparison, `te/check.py`, which ASSERTS rather than prints:
- 63 items, the same keys on both sides, no errors;
- identity on the plan's list;
- the set of items with any difference in elements, `unmatched_members`, `refined_rules` or T equals exactly `{k32_series@2..10}`.

It passed:

```
identical on all 63: ['rules', 'M', 'M1', 'classes', 'active', 'inert', 'split_active', 'split_inert', 'reliable', 'incomplete', 'insplit']
incomplete_footprints non-empty anywhere: []
items with any difference: ['k32_series@2', ..., 'k32_series@10']
k32_series@2: elements 23->24 singletons -[] +['R_(0.0;12)^2[9,9]']; unmatched 0->0; T 0->0 tested 0/0; refined_rules 11->11 []
k32_series@3: elements 23->24 singletons -[] +['R_(0.0;12)^2[9,9]']; unmatched 0->0; T 0->0 tested 2/2; refined_rules 11->11 []
k32_series@4: elements 23->24 (same); tested 3/3
k32_series@5: elements 23->24 (same); tested 4/4
k32_series@6: elements 40->39 singletons -['L_(0.0;1)^4[10,10]'] +[]; unmatched 1->0; T 1->0 tested 9/9; refined_rules 20->20 []
k32_series@7: elements 40->39 singletons -['L_(0.0;1)^4[10,10]'] +[]; unmatched 1->0; T 1->0 tested 13/13; refined_rules 20->20 []
k32_series@8: elements 113->111 singletons -['L_(0.0;3)^5[10,10]', 'L_(0.0;5)^2[9,9]'] +[]; unmatched 0->0; T 0->0; refined_rules 55->54 [L#2<->L#4 child "4" gone; word of R#23<->R#27 changed]
k32_series@9: same as 8
k32_series@10: same as 8
```

Labelled `M` (class names and entries) and `M_1` are inside `M`/`M1` above and are identical.

**The one word change at steps 8-10**, token by token:

```
base   ... ['A:p1@0.0/L#2 <-> A:p1@0.0/L#4', 4, 1] ...
branch ... ['A:p1@0.0/L#2 <-> A:p1@0.0/L#4', 3, 1] ...
changed tokens: [(3, [.. L#2<->L#4, 4, 1], [.. L#2<->L#4, 3, 1])]
```

That is the inert `v_4` → `v_3`, and the only change in that word.

**Agreement with plan §5:**
- steps 2-5: 23 → 24 elements, crossing 9 pinched;
- steps 6-7: 40 → 39, `[10,10]` gone; `unmatched` 1 → 0; T 1 → 0;
- steps 8-10: 113 → 111, left `[10,10]` and `[9,9]` gone; `refined_rules` 55 → 54; `v_1..v_4` → `v_1..v_3`; one token `v_4` → `v_3`.

Nothing differs on any other build or step, and `incomplete_footprints` is empty everywhere. **Met.**

## Rubric 5: the grep and the line count

- `grep -n "branch_key != base_branch" src/tanglepack/topology/PartitionFamily.py` is empty (rc=1).
- The code diff of `_empty_stretches`, comments and docstring excluded, from `git diff -U0`:

  ```
  -            (base_branch,) = ends
  -            for branch_key, on_branch in folds.items():
  -                if branch_key != base_branch:
  -                    on_branch.sort(key=lambda c: trellis.intersection(c).stable_cdist)
  +            for on_branch in folds.values():
  +                on_branch.sort(key=lambda c: trellis.intersection(c).stable_cdist)
  ```

  That is 4 removed lines and 2 added. Counted as old lines touched it is 4, which meets "at most 4 lines". Counted as `+` and `-` lines together it is 6, which exceeds it.

  The plan's "at most 4" did not foresee deleting `(base_branch,) = ends`. That line became dead once the test went: `base_branch` had no other reader, and `assert len(ends) == 1` stays. The plan allowed renaming it "only if it becomes unused", and it did become unused. **Deviation from the letter of the plan, for the reviewer to judge.** I read it as within its intent, since the result is shorter than the plan's own sketch.

## Rubric 6: Theorem T on k32

**Branch:**
```
k32_series p1_tested per step: [0, 0, 0, 2, 3, 4, 9, 13, 13, 13, 13]
```
The test passes with 0 disagreements on steps 0-10.

**Base, clean:**
```
k32_series p1_tested per step: [0, 0, 0, 2, 3, 4, 9, 13, 13, 13, 13]
XFAIL ...test_members_share_the_P1_itinerary_under_blasting[k32_series] - B7b: ...
1 xfailed in 27.73s
```

The `p1_tested` lists are equal. On base with the branch's `cases.py` copied in (no known issue), it fails at steps 6 and 7 only, on class `k = R#7 <-> R#9`. Member `[(6, 7)]` lands in `v(L_(0.0;1)^4, L_(0.0;3)^2) l^-1(R_(0.0;5)^2, R_(0.0;19)^4)`, while its siblings land in `v(L_(0.0;1)^5, L_(0.0;3)^1) l^-1(R_(0.0;5)^2, R_(0.0;19)^3)`.

The other series' T runs pass on the branch with these counts: k28 `[0,0,0,2,3,4,4,4,4]`, p3 `[0,0,0,0,0,0,0,2,3,3]`, inversion `[0,0,2,3,4]`, nested_outer `[0,0,0,2,3,4,4,4,4,4,4,4]`. **Met.**

## Rubric 7: the series test across the two worktrees

**Branch** (`-s`; all five pass):

```
k28_series chords_checked per step: [0, 0, 1, 1, 1, 1, 1, 1, 1]
p3_series chords_checked per step: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
inversion_series chords_checked per step: [1, 1, 1, 1, 1]
k32_series chords_checked per step: [0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2]
nested_outer_series chords_checked per step: [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
27 passed, 1 skipped   (this run also held the T series, the build test and test_iterated_cut.py)
```

On `k32_series`, `chords_checked` is at least 1 on every step 2-10. The series sum is 7 + 0 + 5 + 14 + 10 = **36**.

**Totals against the analyst's 47.** The build test reads 4 lobes (rubric 8), so the tests make 36 + 4 = **40 readings**. To reconcile, I ran `base_chord_misreadings` on the branch over all 16 builds (`te/one_build.py`, `bc_*.txt`). It found 11 lobes and 0 misreadings:
- 1 lobe each: k10, k28_two_blasts, k28_8, nested, nested_inner_blasts, nested_inner6, inversion, inversion_one_blast, inversion_two_blasts, inversion_period2, heteroclinic;
- 0 lobes: k28_one_blast, p3, p3_15, p3_16, p3_6_blasts.

11 + 36 = **47, exactly the analyst's count**. The 7 lobes the tests do not read are on builds that `test_base_chords_are_inside_the_lobe_image` does not parametrise, and the helper reads all 7 without a misreading.

**Base** (the branch's `tests/helpers/laws.py`, `tests/cases.py`, `test_class_count_under_blasting.py`, `test_foreign_chords.py` and `test_iterated_cut.py` copied in, library untouched, restored afterwards). The `chords_checked` lists are identical to the branch's. The test FAILS on `k32_series` only:

```
step 2: hole bridge (6, 7) on ('stable', 0, 0): (c) the inside gap (9,12) is not empty; (anchor,13)out (13,8)B (8,11)out (11,10)in (10,9)out (9,12)in (12,tail)out
step 2: hole bridge (6, 7) on ('stable', 0, 0): (d) its chord (11,12) is not one inside gap; <same readings>
steps 3, 4, 5: the same two lines
step 6: (c) the inside gap (9,12) is not empty; (d) its chord (11,12) ...; (d) its chord (9,10) is not one inside gap
step 7: the same three lines
step 8: (d) its chord (11,12) is not one inside gap; (d) its chord (9,10) is not one inside gap
steps 9, 10: the same two lines
```

The readings are identical at every step. The plan predicted clause (c) at 2-7 and clause (d) at 8-10.

**Deviation from the prediction, confirmed:** clause (d) fires at **2-10**, not only at 8-10.
- Chain chord `(11,12)` spans three gaps, so it is never one inside gap, from step 2 on.
- Chain chord `(9,10)` is flagged at 6-10. At 2-5 it is not flagged: as the analyst's findings note, `(10,9)` there is another hole's stretch, which (d) does not test. (d) tests only this lobe's own base and chords, read through `_empty_stretches` over the hole bridge alone.

This makes the test stricter than the plan predicted, not weaker. The second k32 lobe at steps 6-10 (`chords_checked` 2) is read without a misreading on base too, because its orders agree.

## Rubric 8: the build test

**Branch** (`base_chord_misreadings` on each build):
```
k10 (1, ())
inversion (1, ())
inversion_period2 (1, ())
heteroclinic (1, ())
p3_6_blasts (0, ())
SKIPPED [1] tests/regression/test_foreign_chords.py:106: no stretched hole crosses its base branch twice
```

The four builds pass with `checked > 0`, and the SKIP is visible. Before skipping, the test asserts `not checked`, so the skip fails loudly if `p3_6_blasts` ever gains a base lobe.

On base (copied in), the build test also passes, 4 passed and 1 skipped. No build has a lobe where the two orders differ. The plan expects this (§8 Q2), and the series test is what carries teeth on base. **Met.**

## Rubric 9: the unit pin

On base (copied in) it FAILS:
```
test_base_crossings_pair_in_stable_order: AssertionError: assert ((10, 'chord') in [(12, 'chord')])
```

The full `empty` map of that chain on each side, from a throwaway test file that was deleted afterwards:
```
base   {1: [(2,'hole')], 2: [(1,'hole')], 8: [(13,'base')], 9: [(10,'chord')], 10: [(9,'chord')], 11: [(12,'chord')], 12: [(11,'chord')], 13: [(8,'base')]}
branch {1: [(2,'hole')], 2: [(1,'hole')], 8: [(13,'base')], 9: [(12,'chord')], 10: [(11,'chord')], 11: [(10,'chord')], 12: [(9,'chord')], 13: [(8,'base')]}
```

Base pairs the chain's `(9,10),(11,12)`, and the branch pairs the stable order's `(11,10),(9,12)`. It passes on the branch. **Met.**

## Rubric 10: the six call sites

Call sites:
- `test_iterated_cut.py:83` (helper), `:166-177` (wrapper);
- `test_iterated_partition.py:70`, `:90-100` (wrapper);
- `test_foreign_chords.py:60`;
- `test_partition_family.py:169-170` (with `_incomplete_footprints`).

`test_iterated_partition.py` and `test_partition_family.py` are byte-identical to base (`git diff --quiet a6a4f68..`). Running the four files on the branch gives `32 passed, 1 skipped`. `test_an_odd_fold_count_is_an_incomplete_footprint` is among them.

**`_incomplete_footprints` has no code change.** Its only diff is one docstring word, "last" → "outermost" (`PartitionFamily.py:1018`). **That edit is outside the plan's list**, which named the module docstring and the `_empty_stretches` Note for this word, not this docstring. It is accurate: with stable order on every branch, the crossing left unpaired on an odd count is the outermost one. As the coder says, its log text and `describe()` text still say "last"; they are code and were left untouched. **Met, with that docstring deviation reported.**

## Rubric 11: the greps (branch)

- `grep -n -i "chain order" src/tanglepack/topology/PartitionFamily.py` returns `91:  footprint (2026-09-21) -> empty stretch, base chords in chain order`. That is the History line only. The coder reworded a second mention in the same bullet to "pairing along the chain".
- `grep -n -i "chain order" tests/unit/topology/test_iterated_cut.py` is empty.
- `grep -n "B7b\|known issue" src/tanglepack/topology/SymbolicDynamics.py` returns `30:order since B7b (``PartitionFamily``).`, which makes no known-issue claim.
- `grep -rn "R_7 <-> R_9" tests/` is empty.

**Met.**

## Rubric 12: the k=2.8 example

`session.iterated_partition()` on `k28_two_blasts` at the branch, run from the branch worktree (`tanglepack.__file__` = the branch's `src`):

```
left  L_1^1=[anchor,4) L_1^2=[4,9] L_2=(9,10) L_3=[10,3]
right R_1^1=[anchor,11] R_1^2=(11,8) R_1^3=[8,7] R_2=(7,6) R_3=[6,5] R_4=(5,4) R_5^1=[4,9] R_5^2=(9,10) R_5^3=[10,3]
```

This matches the new module docstring (`PartitionFamily.py:53-58`):
- hole `(6,7)`, base `(8,11)`, chord `(9,10)`;
- `R_1^1=[0,11] R_1^2=(11,8) R_1^3=[8,7]` and `R_5^1=[4,9] R_5^2=(9,10) R_5^3=[10,3]`;
- `L_1^1=[0,4) L_1^2=[4,9]`.

`describe()` gives the same intervals with their cuts, e.g. `#1 [5.993, 6.517] (id 4 -> id 9) (parent #0, cut by (3, 4))` on the left. **Met.**

## Break-it checks (branch, each restored with `git checkout -- src/tanglepack/topology/PartitionFamily.py`)

Each run covered the five series cases, the five build cases and the unit pin.

| Break in `_empty_stretches` | Result |
|---|---|
| **Base code** (chain order on the base branch), i.e. the base commit | Series `[k32_series]` FAILS (c)+(d) at 2-10; unit pin FAILS; T `[k32_series]` FAILS at 6-7 (rubric 6, 7, 9 above) |
| **A**: no chords at all (`stretches.setdefault(..., "chord")` → `pass`) | 2 failed, 8 passed, 1 skipped: series `[k32_series]` and the unit pin fail. The other four series and the four builds PASS |
| **B**: pairing shifted by one (`zip(on_branch[1::2], on_branch[2::2])`) | 2 failed, 8 passed, 1 skipped: the same two fail, the rest pass |
| **C**: a wrong chord, `frozenset((crossings[0], y))` instead of `(x, y)` | 10 failed, 1 skipped: all five series, all four builds, the unit pin. k10 fails on `(d) its chord (4,6) is not one inside gap; (anchor,7)out (7,4)B (4,5)out (5,6)in (6,tail)out`. `p3_series` (0 lobes read) fails collaterally on a library assertion (`PartitionFamily.py:841`, `((40, 41), 40, 42, 'chord')`), not on a misreading |

**Teeth limitation (observation, not a rubric gap).** Outside the k32 lobe, the base and build lobes have exactly two interior crossings on the base branch. Their one inside gap is also another hole's stretch (k10's lobe folds through the reference hole `(4,5)`: `(5,6)in`). So clause (c) is satisfied by the `hole` stretch even when no chord is declared, and breaks A and B (which only REMOVE chords there) do not change those cases' partitions. The tests correctly pass on them. A wrongly ADDED chord (break C) is caught everywhere through clause (d).

The plan already expected that only the k32 series lobe distinguishes chain from stable order (§8 Q2). I did not add a test, since no rubric item lacks one, and a hand-built lobe whose inside gap is not a hole's stretch is already covered by the unit pin.

## Deviations from the plan, in one place

1. **Rubric 5 line count.** The `_empty_stretches` code diff is 4 removed and 2 added lines, 6 counting both. The extra removal is the dead `(base_branch,) = ends`, which the plan allowed removing once unused.
2. **`_incomplete_footprints` docstring** "last" → "outermost". The plan did not list it. It is docstring only and accurate, and its log/`describe()` text still says "last".
3. **Rubric 7 base failure clause (d)** fires at steps 2-10, not only 8-10 as predicted: `(11,12)` at 2-10, `(9,10)` at 6-10. This is stricter than predicted.
4. **Diff size** is 360 changed lines against the plan's estimate of about 200, in the same seven files.
5. **A plan-directed assertion removal** in `test_foreign_crossings_pair_in_stable_order` (the tie `(12,"chord") in empty[11]`). I verified it was toothless: all cdists are equal at 1.0.

Nothing could not be tested.
