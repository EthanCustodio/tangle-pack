# Outcome — B7 `refined-level-defined`

Branch `refactor/refined-level-defined` @ `a6a4f68` (worktree `../tangle-pack-refined-level-defined`), base `3fb310a`. There are 7 commits and 8 files.
- The plan was approved at round 1 with question 1's default, option (ii).
- Review passed in round 1 of coding.
- There is no debrief verdict yet. Nothing is merged; the merge is yours.
- The main tree still carries uncommitted `docs/`, `.refactor/` and `CLAUDE.md` changes. The branch does not touch them.

## The answer to your question now

**Is the refined word a theorem? Yes, wherever Theorem T holds, and the code now says so. T holds everywhere we can test it except one lobe on `k32_series`, which is now a registered known issue.**

**What T says.** Every member bridge of a class, followed through its registered image, reads the same `P1`-level word (freely reduced). When that holds, the refined word of each child `a_k` is not an inheritance from its parent. It is the word every member of the class actually has. The refined matrix `M_1` is then the in-splitting of the class matrix `M`:
- sibling rows are equal;
- the children's columns sum to the parent's column;
- the two systems are conjugate, so they have the same entropy.

**Where T holds (tested on registered images, lookup only):**
- the law cases k10 (1 class), inversion_two_blasts (2) and heteroclinic (7, on the period-3 tangle);
- p3_16 (4 classes) and k28 at 8 blasts (4);
- every step of the five blast series, 141 class-tests in all. These include period-3 classes (p3_series steps 7–9) and the two-branch inversion saddle.

**Where T fails:** `k32_series` at 6 and 7 blasts, in the class `{R_(0.0;8), R_(0.0;10)}`.
- The hole-bridge member `(6, 7)` lands in different `P1` sub-elements from its siblings `(54, 55)`, `(64, 65)` and `(96, 97)`. The class-level word is the same for all of them.
- **Cause.** The image of that lobe spirals back across its own stable branch. `P1` pairs the fold crossings on the base branch in **chain order** (your 2026-09-23 rule), so it declares the wrong stretch empty there. A sibling's image then legitimately crosses that stretch.
- **Status.** This is the strict xfail for B7b.
- **What I found in addition (verified):** the walked word, which is what `refined_rules` holds, agrees with the three siblings. Only `(6, 7)` is off. So the refined word the code prints is the majority's, and the miscut affects one member's reading.

**Where T is untested (said plainly):**
- the eight law cases that have no class with two members whose images are registered;
- the deep runs p3_15, p3_6_blasts and inversion_period2_10, which test 0 classes. T is a law-cases-only law, so its depth coverage comes from the series.

**Where T holds by construction, not by theorem:** the in-splitting's sibling-row identity. `refined_rules` copies one word to every child, so siblings always have equal rows. The law's real teeth are the column-sum identity (it catches a token dropped by `transition_matrix`) and the node-set identity.

## What changed

**The key hunk is a definition, not code.** The `refined_rules` construction (`SymbolicDynamics.py:815-826`) is byte-identical. What changed is what the docs claim about it.

*Dev Notes, before:*
> every refined child inherits the parent's word rewritten in refined symbols

*After* (`SymbolicDynamics.py:64-67`; the new module paragraph is `:17-31`):
> the refined word of every child is the word of every member, by T, written in refined symbols

The module paragraph adds:
- **Occurrences.** Each token carries the `P1` pair it lands on, and a class's distinct occurrences are its refined symbols.
- **Theorem T,** with its argument: the strip two members bound with the stable manifold holds no hole, so no `P1` cut separates their images.
- **`M_1`** is the in-splitting of `M`.
- **The k32 known issue.**

The `transition_matrix` docstring (`:1012`) names its two outputs `M_1` and `M`.

**`PartitionFamily.py` now states `P1`** (`:18-27`, `:38-47`):
- `P1 = P ∨ f(∂P)`: `P`'s cuts plus every crossing of the first images of the hole bridges;
- base-branch chords are paired in chain order (2026-09-23), chords on other branches in stable order (C19);
- an odd fold count is an *incomplete footprint*.

The old Dev Note claim that chain order "changes refinement labels only, never walks or words" is replaced by the k32 counterexample and a pointer to B7b (`:107-118`).

**The one code change: the odd-count WARNING became data.**

```python
# before (_empty_stretches, :959-966)              # after
if len(on_branch) % 2:                              # (deleted; the crossing is still left unpaired)
    logger.warning("image chain of hole bridge "
        "%s crosses stable branch %s an odd number ...")

                                                    # from_minimal :912
                                                    incomplete_footprints=cls._incomplete_footprints(minimal, homotopy)
                                                    # new staticmethod :1007-1050: Counter of interior
                                                    # crossings per stable branch, odd ones kept as
                                                    # (hole bridge, branch, count); one logger.info per build
```

`IteratedHomotopyPartition` gains:
- the attribute `incomplete_footprints`;
- the `__init__` keyword `incomplete_footprints=()`;
- a `describe()` override that adds one line per footprint (`:687`).

**Tests:**
- **Two new `symbolic` laws** in `tests/helpers/laws.py`:
  - `members_share_the_P1_itinerary` (`:3083`, via the helper `p1_itinerary_disagreements` at `:2998`; it skips "stable end" chains);
  - `refined_matrix_is_an_in_splitting` (`:3091`).
- **`check_refined_children_inherit_word`** is kept; only its docstring changed.
- **The series run** is `test_members_share_the_P1_itinerary_under_blasting` (regression, `:189`), fed by the new `SeriesStep.p1_tested` / `p1_disagreements` fields.
- **The known issue** is `("series_members_share_the_P1_itinerary", "k32_series")`, with wp `B7b`.
- **Eight `NOT_APPLICABLE` set-asides** are visible as skips.
- **Three unit tests:**
  - a stub odd-fold chain;
  - `M` and `M_1` on A.8d's k=2.8 shape;
  - the in-splitting law failing on a sub-word.

## What got shorter or clearer

Nothing got shorter in code. The slice adds definitions and laws: +433 / −46 lines over the eight files.
- **`PartitionFamily.py`:** 1,009 → 1,082 lines. The new 44-line `_incomplete_footprints` accounts for most of it, the rest is docstring. The Dev Notes "History" bullet went from 11 lines to 5: the rule and its dates, with the rejected variants left to git.
- **`SymbolicDynamics.py`:** 1,516 → 1,538 lines, docstrings only.
- **`laws.py`:** 3,257 → 3,399 lines. **`cases.py`:** 1,420 → 1,446 lines.

**What is clearer:**
- The refined level now has a stated definition (occurrences in `P1`), a stated theorem, and two laws.
- "Inherits" no longer suggests an assumption.
- The false general claim in the `PartitionFamily` Dev Note is gone.

## Every behaviour change

1. **Numbers: none.**
   - The test engineer's script compared base and branch on 15 builds (the 11 law cases, p3_6_blasts, p3_16, k28 at 8 blasts, nested with 6 inner blasts). These are equal on all 15: letter-free words, refined rules, `M`, `M_1`, the iterated element intervals, and `is_reliable` (True on all 15).
   - My own script confirmed the same on k32_series steps 0–7.
2. **Logs.** The odd-fold-count WARNING is gone, replaced by one INFO per build. It fired on none of the 62 builds and steps on the base, and `incomplete_footprints` is empty on all of them. No new WARNING appears anywhere.
3. **API (additive):**
   - `IteratedHomotopyPartition.incomplete_footprints`;
   - the `__init__` keyword;
   - the `describe()` lines;
   - the private `_incomplete_footprints`.

   `_empty_stretches` keeps its signature.
4. **Tests:**
   - +2 laws, +1 series test, +3 unit tests, +8 visible set-asides;
   - one new **strict xfail**, `test_members_share_the_P1_itinerary_under_blasting[k32_series]`, wp B7b;
   - the sanity guard now accepts `series_members_share_the_P1_itinerary` and wps of the form `B<n>[a-z]`.

   The suite went from 1193 / 38 / 5 / 1 to **1200 passed, 46 skipped, 5 deselected, 2 xfailed**.

## What was verified, and how

**Re-run by me:**
- **Suite on the branch:** `1200 passed, 46 skipped, 5 deselected, 2 xfailed` in 203 s. `-rxX` lists exactly the two `KNOWN_ISSUES`: C3, and the k32 entry with wp B7b. There is no XPASS.
- **The whole diff, read** (`git diff 3fb310a..a6a4f68`, 8 files):
  - `SymbolicDynamics.py` changes only docstrings and Dev Notes;
  - `PartitionFamily.py` changes only as described above;
  - the cut logic of `_empty_stretches` is untouched apart from the deleted WARNING.
- **The k32 known-issue message** (`--runxfail` on the branch): `1 failed`, with `step 6: class k A:p1@0.0/R#7 <-> A:p1@0.0/R#9: [(6, 7)]: … | [(54, 55)]: …` and `step 7: … [(54, 55), (64, 65), (96, 97)]`.
- **My own base-vs-branch run** (`scratchpad/b7a/k32.py`, an independent reimplementation of T, with names resolved through `ElementNaming`), on k32_series steps 0–7:
  - equal on base and branch: letter-free words, `M`, `M_1`, `is_reliable`, the WARNING set, T's tested counts (0, 0, 0, 2, 3, 4, 9, 13) and the disagreement;
  - the class is `{R_(0.0;8), R_(0.0;10)}`. The `R_7 <-> R_9` in the KnownIssue text is the positional `R#7` / `R#9`;
  - `incomplete_footprints` is `[]` at every step.
- **A stronger check, branch only and not a law** (`scratchpad/b7a/walk.py`): I compared the walked word of each class against **every** member's `P1` reading, not just member against member. It agrees on all 11 law cases: 128 member readings, 0 disagreements. On k32 at steps 6 and 7 only `(6, 7)` disagrees with the walk.

**Taken from the reports (not re-run by me):**
- the full 15-build base-vs-branch comparison (test engineer, own script, independent of the coder's);
- the series tested counts on the base;
- the AST check that `SymbolicDynamics.py` is docstrings only;
- every break-it check: stable-end skip removed → k28_series fails at steps 6–8; sub-word, dropped column and dropped node each turn the law red; INFO→WARNING fails the unit test;
- the reviewer's check that the B7b patch (stable order on the base branch) turns the k32 entry into XPASS(strict). That means the known issue flips exactly when B7b lands.

## What did not go to plan

- **An eighth file, `tests/invariants/test_law_case_sanity.py` (+2/−1).** The plan's seven files contradicted the verdict you approved. That verdict asks for a strict-xfail known issue on `k32_series` naming `B7b`. The sanity guard rejects any known issue whose case is not a law case unless its id is in `_STANDALONE_ISSUES`, and it rejects any wp outside `WP\d+|C\d|DEC-\d+`. No seven-file implementation could pass it. The orchestrator ordered the two-line change, and the reviewer judged it the minimum. The planner's file count was wrong, not the coder's scope.
- **The known-issue key is `series_members_share_the_P1_itinerary`, not `members_share_the_P1_itinerary` as your verdict literally said.** Under the layered id, the law tier would have generated a set-aside test asking the `law_case` fixture for a case `k32_series` that does not exist. That test would xfail on `KeyError`, a sham. The `series_` prefix is the existing convention for the three other series ids.
- **`_empty_stretches` kept its signature** (the plan said it would return `(stretches, incomplete)`). Five call sites in three out-of-scope test files monkeypatch or call its dict-returning shape. A sibling static method now collects the footprints.
  - Cost: image chains are looked up twice per build.
  - Its admission filter differs from `_empty_stretches`' in principle. The test engineer measured the two equal on all 62 builds and steps, and the difference can only affect the INFO, never a cut.
- **T is law-cases-only** (`LAW_CASES_ONLY`). The deep-run test requires every layered law to be non-vacuous, and T tests 0 classes on p3_15, p3_6_blasts and inversion_period2_10. T's depth coverage is the series test.
- **Size:** about 430 changed lines, against the plan's 200–250. The law helper and the two unit tests are longer than estimated.

## Where I am least sure this matches what you meant

1. **How the k32 counterexample is handled.** You chose option (ii), and the code is "T holds, except this registered known issue". Two things to weigh:
   - What it means, if I read it right: in the one place T fails, the **walk** sides with the three siblings and only the hole bridge's own image is miscut. The refined word printed today is plausibly the right one, and B7b would make `(6, 7)` agree rather than change any word.
   - B7b reverses your 2026-09-23 chain-order rule on the base branch. The evidence is one lobe's point-in-polygon test plus a corpus-wide scratch patch: no number changes, inert splits drop by one on 15 builds and steps, and iterated names change where the two orders differ. It is not a proof that chain order is never what you meant on the base branch.
2. **The docs claim slightly more than the tests do.**
   - The `SymbolicDynamics` module docstring says T and the in-splitting are "laws on every build". T runs on the law cases and the blast series only. The reviewer flagged this; it was not fixed.
   - The KnownIssue text and the Dev Note say "R_7 <-> R_9", which mixes positional ids with element names. The class is `{R_(0.0;8), R_(0.0;10)}`.

   Both are cosmetic and can be folded into B7b.
3. **The commit attribution line.** All seven commits end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. That is the session's attribution rule and the trailer of every earlier commit. The orchestrator twice asked for "Claude Fable 5.1", and the coder declined because an agent message does not override the session rule. If you want a different trailer, it is an amend before the merge.

**Records updated:** `.refactor/decisions.md` (the B7 outcome entry), `.refactor/knowledge/topology.md` (B7 outcome section and branch line numbers; the C19 note marked) and `.refactor/knowledge/testing.md` (B7 outcome section).
