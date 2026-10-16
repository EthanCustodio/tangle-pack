# Outcome: B7b `base-chords-in-stable-order` (auditor)

Branch `refactor/base-chords-in-stable-order` @ `63477e3`, base `a6a4f68`. Five commits. **Not merged**: merging is your call. Both worktrees are clean. On main, the working tree carries only the uncommitted records (`.refactor/`, `docs/`, `CLAUDE.md`, `artifacts/`), as before.

## Your question, answered now

**Which stretches of a lobe image's own (base) stable branch are empty?** The interior crossings paired in **stable order**, the same rule as on every foreign branch.

- **The code.** `IteratedHomotopyPartition._empty_stretches` now has one chord rule for every branch (`PartitionFamily.py:988`). Your 2026-09-23 chain-order clause for the base branch is reversed.
- **Why.** The lobe image `f(H)` holds no anchor. So walking any branch from its anchor enters and leaves `f(H)` in turn, and the base branch is no exception (Jordan, as in C19).
- **The measurement, from the analyst.** There are 47 base-branch lobe readings in the corpus (16 builds and 47 series steps). On all 47:
  - the gaps alternate out/in, counting the base as inside;
  - the covered gaps are exactly the stable-order pairs.
- **Where the two orders actually differ.** Only on one lobe: `k32_series` hole `(6,7)`, at blast steps 2-10 (not only 6-7).
  - Chain order there declared `(9,10)` empty, but that gap lies outside `f(H)`.
- **The stop rule did not fire.**
- **Theorem T** now holds on every series step. The strict xfail on k32 is deleted, and `KNOWN_ISSUES` is down to C3 alone.

## What changed

Seven files, +261 / −99, exactly the plan's list. Nothing under `docs/`, `.refactor/`, `CLAUDE.md` or `artifacts/` changed on the branch.

The most important hunk is the whole library change (`PartitionFamily.py`, `_empty_stretches`):

```python
# before (a6a4f68)
            (base_branch,) = ends
            ...
            for branch_key, on_branch in folds.items():
                if branch_key != base_branch:
                    on_branch.sort(key=lambda c: trellis.intersection(c).stable_cdist)
                for x, y in zip(on_branch[::2], on_branch[1::2]):
                    stretches.setdefault(frozenset((x, y)), "chord")
# after (63477e3)
            for on_branch in folds.values():
                on_branch.sort(key=lambda c: trellis.intersection(c).stable_cdist)
                for x, y in zip(on_branch[::2], on_branch[1::2]):
                    stretches.setdefault(frozenset((x, y)), "chord")
```

**Shorter and clearer:**
- `PartitionFamily.py` goes from 1,082 to 1,077 lines.
- The kernel loses one binding, one branch test and one level of indentation.
- The module docstring states one rule, where it used to state two. It no longer says T has a known issue.
- The k=2.8 example in the docstring used pre-WP8 ids and was wrong at `HEAD`. It now uses current ids: hole `(6,7)`, base `(8,11)`, chord `(9,10)`. The test engineer checked it against `iterated_partition()`.
- The History note is one line: closed cut → footprint → chain order → C19 → B7b.
- `SymbolicDynamics.py:27-30` (docstring only) now says T and the in-splitting are laws "on the law cases and every blast series". The old text said "every build", which overstated it.

**Tests (+10 passing, +1 former xfail now passing, +1 visible skip):**
- **New helper** `helpers.laws.base_chord_misreadings` (`laws.py:3111`). It reads base-branch gaps by point-in-polygon, from existing polylines only; it never calls the map. It checks four clauses:
  - (a) the base is one gap;
  - (b) the readings alternate;
  - (c) every inside gap is empty;
  - (d) every base or chord of this lobe is one inside gap.
- **Where the helper runs:**
  - in every `blast_series` step (`SeriesStep.chords_checked` / `chord_misreadings`), through the new test `test_base_chords_are_inside_the_lobe_image_under_blasting` over the five series;
  - in `test_foreign_chords.py::test_base_chords_are_inside_the_lobe_image` on `k10`, `inversion`, `inversion_period2` and `heteroclinic`, with a visible skip for `p3_6_blasts`.
- **New unit pin** `test_iterated_cut.py::test_base_crossings_pair_in_stable_order`, using the k32 shape built by hand.
- **Deleted:** the toothless chain-order half of the old foreign pin, and the k32 `KNOWN_ISSUES` entry.

## Every behaviour change

All of these are on `k32_series` only. I re-diffed the test engineer's base and branch probe output (`scratchpad/te/out/*.jsonl`, 63 items) with my own script. It gives exactly this list and nothing on any other build or step:

| k32 step | iterated elements | other |
|---|---|---|
| 2-5 | 23 → 24 (crossing 9 is pinched to a singleton `R_(0.0;12)^2 = [9,9]`) | none |
| 6-7 | 40 → 39 (the singleton `L_(0.0;1)^4 = [10,10]` is gone) | `unmatched_members` 1 → 0; T disagreements 1 → 0 |
| 8-10 | 113 → 111 (left singletons `[10,10]` and `[9,9]` are gone) | `refined_rules` 55 → 54 keys: inert `v = {L_(0.0;3), L_(0.0;5)}` has 3 children instead of 4, and one active word's token `v_4` becomes `v_3` |

Element names and intervals around those changes shift accordingly.

These are **identical on all 63** (verified in the same diff):
- letter-free `rules`, `M` and `M_1`;
- class, active and inert counts;
- the number of split active and split inert classes;
- `is_reliable`;
- `incomplete_footprints` (empty everywhere).

**Correction to your request.** The request expected the number of split inert classes to drop by one on 15 builds. That does not happen. The analyst traced that figure to B7's scratch patch, which also dropped the `empty.setdefault(b, ...)` line, and reproduced the whole list from that buggy patch. On the landed code, the split-inert count is equal on both sides at every item (for example 1 and 1 at k32 steps 2-7, 3 and 3 at steps 8-10).

## What was verified, and how

What I ran myself:
- **Full suite on the branch:** `1211 passed, 47 skipped, 5 deselected, 1 xfailed in 207 s`. The only xfail is C3, and there is no XPASS. This matches the test engineer's and the reviewer's runs. The base gives 1200 / 46 / 5 / 2 xfailed (test engineer). The collected ids differ by +11, with none removed.
- **The new unit pin on base:** I copied the branch's `test_iterated_cut.py` into the base worktree and restored it afterwards. It fails with `assert ((10, 'chord') in [(12, 'chord')])`; the other 4 tests pass.
- **The probe diff** described in the previous section.
- **A full read of the diff.**
  - The library change is the 4 → 2 line hunk above plus docstrings.
  - `base_chord_misreadings` reads only `arc_polyline`, root and tail nodes, and crossing coordinates.
  - The six `_empty_stretches` call sites are untouched (`test_iterated_partition.py` and `test_partition_family.py` are byte-identical to base).

From the test engineer and reviewer, read but not rerun by me:
- On base, the series test fails on k32 at steps 2-10: clause (c) at 2-7 and clause (d) at 2-10.
- Break-it runs:
  - the reviewer's reversed pairing turns 8 tests red;
  - a wrong added chord (break C) is caught on every case;
  - dropped or shifted chords (breaks A and B) are caught only on k32 (see the third doubt below).
- The tests read 40 lobes (36 series + 4 builds). The helper reads all 47 analyst lobes clean, including 7 on builds that the build test does not list.

## What did not go to plan

1. **Two of the request's predictions were wrong,** and the planner corrected them before coding:
   - the split-inert "−1 on 15" figure (see the correction above);
   - `refined_rules` is not identical, because one inert token changes at k32 steps 8-10.
2. **Clause (d) fires on base at steps 2-10, not only 8-10 as the plan predicted.** The chain chord `(11,12)` spans three gaps from step 2 on. This makes the test stricter than predicted, not weaker.
3. **The `_empty_stretches` code diff is 4 removed and 2 added lines**, against the plan's "at most 4". The extra removal is the binding `(base_branch,) = ends`, which became dead.
4. **One docstring word outside the plan's list:** in `_incomplete_footprints`, "last" became "outermost". It is accurate, since under stable order the unpaired crossing is the outermost one.
5. **The diff is 360 lines against the plan's estimate of about 200.** Most of the difference is the 133-line helper.

## Where I am least sure this is what you meant

1. **Is stable order now the rule, or just the measured fact?** On base, B7b's teeth come from one lobe: k32 `(6,7)`. Every other base lobe in the corpus has two interior crossings, so the two orders coincide there. On those lobes the single inside gap is also another hole's stretch, so a dropped chord would go unnoticed outside k32.
   - The Jordan argument covers the general case.
   - The reviewer checked that pairing in stable order never pairs across the base: each side of the base has an even number of crossings.
   - Still, if you wanted a proof stated in the docstring rather than "measured on 47 lobes", it is not there.
2. **The base test's clause (c) reads "inside ⇒ empty for any reason".** An inside gap that happens to be another hole's stretch passes, even when this lobe declares no chord there. Clause (d) is the half that pins this lobe's own chords.
3. **The build test reads 4 of the 11 build lobes.** It skips `k28_two_blasts`, `nested`, `nested_inner_blasts` and the blasted inversions. Those are read clean by hand, but no test pins them. This follows the plan's parametrisation; say if you wanted all of them.

## Left open (not in this slice)

- `CLAUDE.md` "Partition families" still says base-branch chords pair in chain order. This is for B12; I did not touch `CLAUDE.md`.
- `_incomplete_footprints` still says "last crossing" in three places, while its docstring now says "outermost":
  - its INFO log, `PartitionFamily.py:1040`;
  - the `describe()` line, `:690`;
  - the class attribute doc, `:653`.
- `_incomplete_footprints` still repeats `_empty_stretches`' chain walk. You declined the fold (Q1 default).
- **Reviewer's nits:**
  - the `misread` closure is redefined inside the loop;
  - `gaps = {}` has no type annotation;
  - gap keys mix `'anchor'`/`'tail'` strings with integer ids;
  - `_lobe_image` keeps its untyped parameters (moved verbatim);
  - the helper imports the private `_hole_stretches`;
  - the `functools.cache` in `test_foreign_chords.py` now keeps four big builds alive for the whole run.
- `base_chord_misreadings` silently skips a non-consecutive image chain. None occurs on the corpus.
- **Commit trailers:** all five commits end with `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`.

**Records updated (uncommitted, on main):**
- `.refactor/decisions.md`: the B7b outcome entry.
- `.refactor/knowledge/topology.md`:
  - the C19 note is marked superseded;
  - a B7b outcome section, with the inert-split correction confirmed on the landed code.
- `.refactor/knowledge/testing.md`:
  - the known issue is marked done;
  - a B7b outcome-audit section (the helper, where it runs, its teeth).
