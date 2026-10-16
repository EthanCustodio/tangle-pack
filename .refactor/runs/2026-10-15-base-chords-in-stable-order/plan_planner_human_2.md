# Plan: B7b `base-chords-in-stable-order`

Base: `main` at `a6a4f68`. Branch: `refactor/base-chords-in-stable-order`, worktree `../tangle-pack-base-chords-in-stable-order`. Base worktree for the evidence: `../tangle-pack-base-chords-in-stable-order-base`, detached at `a6a4f68`. Your uncommitted `docs/`, `.refactor/`, `CLAUDE.md` and `artifacts/` are not touched.

## 1. The answer

On a lobe image's own (base) branch, the stretches that are empty are the interior crossings paired in **stable order**: the same rule C19 already uses on foreign branches.

The analyst read all 47 lobes in the corpus that have 2 or more interior crossings on their base branch: 16 builds and 47 series steps. They used point-in-polygon on each gap, closing the gap at its two crossings. In every lobe the gaps alternate out/in, walking from the anchor, with the base counted as inside. In every lobe the covered gaps are exactly the stable-order pairs.

Only one lobe in the corpus makes the two orders differ: `k32_series` hole `(6,7)`, at blast steps 2 to 10. Its chain runs `8→9→10→11→12→13` and its interior stable order is `11,10,9,12`.
- Covered (inside the lobe): `(11,10)` and `(9,12)`, the stable-order pairs.
- Your chain-order pairs: `(9,10)` lies outside the lobe, and `(11,12)` spans three gaps.

The stop rule did not fire.

## 2. Recommendation

**Change it.** Delete the one test `if branch_key != base_branch:` in `IteratedHomotopyPartition._empty_stretches`, so the chords on every branch are sorted by stable cdist. This reverses your 2026-09-23 chain-order rule on the base branch.

Two of the request's predictions were wrong. Both came from a bug in B7's scratch patch, which dropped one `setdefault` line. The correct patch:
- changes no split counts anywhere (the "−1 split inert on 15 builds" does not happen);
- changes `refined_rules` in one inert detail on k32 steps 8-10;
- changes iterated element names on k32 steps 2-10, not just 6-7.

Section 5 has the full list.

## 3. Scope

Slug: **`base-chords-in-stable-order`**. Seven files, about 200 changed lines:

1. `src/tanglepack/topology/PartitionFamily.py`: the one-line change, plus its docstrings.
2. `src/tanglepack/topology/SymbolicDynamics.py`: module docstring only.
3. `tests/cases.py`:
   - delete the k32 `KNOWN_ISSUES` entry;
   - add two fields to `SeriesStep`, filled in by `blast_series`.
4. `tests/helpers/laws.py`: one geometric reading, `base_chord_misreadings(session)`, placed beside `p1_itinerary_disagreements`.
5. `tests/regression/test_class_count_under_blasting.py`: one series test.
6. `tests/regression/test_foreign_chords.py`: one base-branch test over builds, plus the module docstring.
7. `tests/unit/topology/test_iterated_cut.py`: re-pin with the k32 shape, plus the module docstring.

Files 4 and 5 are two more than the request listed. They are needed because the lobe that matters only exists at a series step, and series steps carry records, not sessions. The reading therefore has to run inside `blast_series`, which is how T is already checked there. See Question 2.

Not touched:
- `test_law_case_sanity.py`: the standalone id stays in `_STANDALONE_ISSUES`. It is an allow-list, and four of its six ids have no entry today.
- `test_partition_family.py` and `test_iterated_partition.py`: there is no fold (Question 1).

## 4. Steps

1. **The rule** (`PartitionFamily.py:985-996`). Delete the test so that the sort runs on every branch:
   ```python
   for on_branch in folds.values():
       on_branch.sort(key=lambda c: trellis.intersection(c).stable_cdist)
       for x, y in zip(on_branch[::2], on_branch[1::2]):
           stretches.setdefault(frozenset((x, y)), "chord")
   ```
   - Rewrite the comment above the loop as one sentence: the lobe image holds no anchor, so walking any branch from its anchor enters and leaves it in turn, and the base branch is no exception (Jordan; B7b, 2026-10-15).
   - `base_branch` is still needed for the assertion. Rename it only if it becomes unused.
2. **Docstrings in `PartitionFamily.py`.**
   - Module bullets `:39-41` state one rule: stable order on every branch.
   - The odd-count clause `:44-45` says "outermost crossing", not "last".
   - The `_empty_stretches` Note gets the same two edits.
   - The k=2.8 example `:49-58` is restated in current ids, which I checked at `HEAD`:
     - hole `(6,7)`, base `(8,11)`, chord `(9,10)`;
     - right row `R_1^1=[0,11] R_1^2=(11,8) R_1^3=[8,7]` and `R_5^1=[4,9] R_5^2=(9,10) R_5^3=[10,3]`;
     - anchor lobe fold `(3,4)` abuts hole `(4,5)` at `4`, giving `L_1^1=[0,4) L_1^2=[4,9]`.
   - Dev Notes "History" gets one condensed line: closed cut `2a2e35a` → footprint 2026-09-21 → empty-stretch with base chords in chain order 2026-09-23 → foreign chords in stable order C19 2026-10-09 → stable order everywhere B7b 2026-10-15. That last change rests on the k32 `(6,7)` lobe, where chain order declared `(9,10)` empty although it lies outside `f(H)`.
   - Dev Notes "Foreign chords" becomes "Chords". The known-issue and "slice B7b" sentences `:109-116` go, and the class is named `{R_(0.0;8), R_(0.0;10)}` if it is mentioned at all.
3. **`SymbolicDynamics.py:27-31`.**
   - "laws on every build" becomes "laws on the law cases and every blast series".
   - The k32/B7b sentence becomes one line: "the base-branch chords pair in stable order since B7b (`PartitionFamily`)".
4. **The geometric reading** (`helpers/laws.py`). The signature is `base_chord_misreadings(session) -> tuple[int, tuple[str, ...]]`, returning `(lobes checked, one text per misreading)`. It reads only polylines that numerics already computed, and it calls the map never.

   **Which lobes it reads.** For every hole bridge with a stretch whose consecutive image chain has 2 or more interior crossings on its base branch:
   - Order the chain crossings on that branch by stable cdist. The gaps are anchor→first crossing, then between each consecutive pair, then last crossing→branch tail.
   - Read each gap's midpoint by arc length on the stable polyline, closed at the two crossings' own coordinates. Use the anchor's root point and the tail node at the two ends. Do not take a raw node: on k32 some gaps have 0 to 2 nodes.
   - Test the midpoint with `point_in_polygon(tol=0.0)` against the lobe image. The lobe image is the chain's arcs closed along the base, as `_lobe_image` builds it today; that function moves into the helper.

   **What it checks.**
   - (a) The base is a single gap.
   - (b) The readings alternate, with the base counted as inside. When the anchor is a chain crossing, the first gap is the base (the inversion lobes).
   - (c) Every inside gap that is not the base is an empty stretch of the kernel, `_empty_stretches(minimal, homotopy)`.
   - (d) Every chord the kernel declares for THIS lobe is one single inside gap. The lobe's own chords come from the kernel called with only this hole bridge: `SimpleNamespace(trellis=..., hole_bridge_ids=[hb])`, the same harness the unit tests use.

   Clause (d) is the analyst's corrected form, "no outside gap is this lobe's base or chord", stated so that it can fail on either side. The request's form, "no outside gap is empty", is false on the branch too (inversion `(3,2)` and k32 2-5 `(10,9)` are other holes' stretches) and is not used.
5. **Series** (`cases.py`, `test_class_count_under_blasting.py`).
   - `SeriesStep` gains `chords_checked: int` and `chord_misreadings: tuple`, filled in `blast_series` next to `p1_*`.
   - New test `test_base_chords_are_inside_the_lobe_image_under_blasting`, with id `series_base_chords_inside_the_lobe_image`, using `_params` like its neighbours. It prints `chords_checked` per step and asserts that there are no misreadings.
   - Delete `KNOWN_ISSUES[("series_members_share_the_P1_itinerary", "k32_series")]`.
6. **Builds** (`test_foreign_chords.py`). Add `test_base_chords_are_inside_the_lobe_image`, parametrised over `k10`, `inversion` (the anchor-cornered lobe), `inversion_period2`, `heteroclinic` and `p3_6_blasts`.
   - It asserts that `misreadings == ()` and `checked > 0`.
   - `p3_6_blasts` has no base lobe, so it takes a visible `pytest.skip("no stretched hole crosses its base branch twice")`.
   - Builds are shared with the foreign test through a module `functools.cache` over `cases.BUILDERS`, so `heteroclinic` and `p3_6_blasts` are built once.
   - The module docstring says "one statement, two branch kinds".
   - The foreign test itself is unchanged, apart from importing `_lobe_image` back from the helper.
7. **Unit pin** (`test_iterated_cut.py`).
   - New test `test_base_crossings_pair_in_stable_order`: hole `(1,2)`, chain `8..13` on `OWN`, stable cdists `8:0.176, 9:6.786, 10:6.665, 11:6.571, 12:7.418, 13:0.036`.
     - It asserts `(10,"chord") in empty[11]`, `(12,"chord") in empty[9]`, `(10,"chord") not in empty[9]` and `(12,"chord") not in empty[11]`.
     - The analyst's `handbuilt.py` gives `(9,10),(11,12)` at HEAD, so this test fails on base.
   - Strip the existing `test_foreign_crossings_pair_in_stable_order` of its base half. That half is the toothless "keep the chain order" sentence and the tie assertion `(12,"chord") in empty[11]`, where all cdists are 1.0.
   - The module docstring's first bullet says stable order on every branch.

## 5. Behaviour

The 16 builds (the 15 law/deep builds plus `p3_15`) are identical. Changes on series steps only, all on `k32_series`:

- **Iterated partition, steps 2-10.** The `(6,7)` lobe's chords move from `(9,10),(11,12)` to `(10,11),(9,12)`:
  - steps 2-5: 23 → 24 elements (crossing 9 is pinched to a singleton);
  - steps 6-7: 40 → 39 (singleton `[10,10]` gone);
  - steps 8-10: 113 → 111 (left singletons `[10,10]` and `[9,9]` gone).
  - Element names and intervals around them shift accordingly.
- **`unmatched_members`, steps 6-7:** 1 → 0. Member `(6,7)` now matches its child.
- **Theorem T, steps 6-7:** 1 disagreement → 0. The strict xfail is deleted.
- **`refined_rules`, steps 8-10:** the inert class `v = {L_(0.0;3), L_(0.0;5)}` has children `v_1..v_4` → `v_1..v_3`.
  - Member `(9,10)` joins `(93,94)` in `v_3`.
  - One active word's token `v_4` becomes `v_3`.
  - The key count goes 55 → 54.
  - `M` and `M_1` are unchanged, because inert classes are not in them.

Unchanged everywhere (63 builds and steps): letter-free `rules`, `M`, `M_1`, class/active/inert counts, split active and split inert counts, `is_reliable`, and `incomplete_footprints` (the odd count occurs nowhere).

## 6. Not doing

- **Folding `_incomplete_footprints` into `_empty_stretches`** (Question 1).
- **Which holes mark a base and chords.** The anchor's one-sided hole and "unmarked lobe" (WP12) stay as they are.
- **How a cut is applied.** `open_anchorward`/`open_outward`, "existing boundary" and the row invariant stay as they are.
- **Out-of-slice rules:** B9's "one element at a bridge end", `P2` (Q16), the minimal trellis, the walk, `refine`, `match_members`, and the `_image_chain` lookup (K-3).
- **The foreign test's node-midpoint reading**, which skips node-less gaps. It could use the new closed-gap reading, but that would put a new strictness on the foreign branches outside this slice's evidence. It is a candidate for later.
- **A law-tier version of the base-chord reading** (layer × case with NOT_APPLICABLE set-asides). The regression placement mirrors `test_foreign_chords.py` and T's series run.
- **The `CLAUDE.md` "Partition families" bullet**, which still says chain order. That is B12.

## 7. Rubric

1. `MPLBACKEND=Agg env/bin/python -m pytest -q -rxX` passes on the branch and lists exactly one xfail, `orientation_reversing_case_builds`/C3, with no XPASS.
2. Every test the branch leaves unchanged passes on the base commit and on the branch. Run the base worktree's suite and show its count.
3. `git diff --stat a6a4f68..` names exactly the seven files of section 3, and nothing under `docs/`, `.refactor/`, `CLAUDE.md` or `artifacts/`.
4. **Base-vs-branch probe over 63 items** (the analyst's `scratchpad/b7b/probe.py`/`compare.py`, rerun, or an equivalent the test engineer writes): the 16 builds and the 47 series steps.
   - Identical on all 63: letter-free `rules`, labelled `M`, `M_1` (names and entries), class/active/inert counts, split-active and split-inert counts, `is_reliable`, `incomplete_footprints`.
   - Any difference in iterated elements, `unmatched_members` or `refined_rules` is exactly the section 5 list, and occurs on no other build or step.
5. `grep -n "branch_key != base_branch" src/tanglepack/topology/PartitionFamily.py` is empty. The code diff of `_empty_stretches` (not counting its docstring and comments) is at most 4 lines.
6. `test_members_share_the_P1_itinerary_under_blasting[k32_series]` passes with 0 disagreements on steps 0-10. Its printed `p1_tested` equals the base's `[0,0,0,2,3,4,9,13,13,13,13]`. Every other series passes as before.
7. **`test_base_chords_are_inside_the_lobe_image_under_blasting` across the two worktrees:**
   - On the branch it passes for all five series, and `chords_checked` is printed per step. It is at least 1 on every `k32_series` step 2-10, so the lobe where the orders differ is actually read.
   - The test engineer adds these per-step counts to the build counts from item 8 and reports the total next to the analyst's 47 readings.
   - On the base worktree (a copy of the new helper and test only) it FAILS on `k32_series`, with misreadings at steps 2-10: clause (c) at 2-7 and clause (d) at 8-10.
   - The test engineer shows both outputs.
8. **`test_base_chords_are_inside_the_lobe_image`** passes on `k10`, `inversion`, `inversion_period2` and `heteroclinic`, each with `checked > 0`. It shows a visible SKIP for `p3_6_blasts`.
9. **`test_base_crossings_pair_in_stable_order`** fails on the base commit, showing `(9,10),(11,12)`, and passes on the branch.
10. **Six `_empty_stretches` call sites.** The call sites in `test_iterated_cut.py:85,158-170`, `test_iterated_partition.py:70,90-100`, `test_foreign_chords.py` and `test_partition_family.py:169-170` pass unchanged. `_incomplete_footprints` has no code change, and its unit test passes.
11. **Docstring greps.**
    - `grep -n -i "chain order" src/tanglepack/topology/PartitionFamily.py` finds only the History line.
    - `grep -n -i "chain order" tests/unit/topology/test_iterated_cut.py` is empty.
    - `grep -n "B7b\|known issue" src/tanglepack/topology/SymbolicDynamics.py` finds no claim that T has a known issue.
    - `grep -n "R_7 <-> R_9" tests/` is empty.
12. **The k=2.8 example matches.** The new example in the module docstring matches `session.iterated_partition()` on `k28_two_blasts` at the branch: right row `[·,11] (11,8) [8,7] … [4,9] (9,10) [10,3]` and left row `[·,4) [4,9] (9,10) [10,3]`. Print it.

## 8. Questions

1. **Fold `_incomplete_footprints` into `_empty_stretches`?** Default: **no**.
   - Against folding:
     - The footprint count does not depend on order (it is a per-branch `Counter`) and is zero on the whole corpus.
     - The only clean fold changes `_empty_stretches`' return shape, which means editing all six call sites in four test files. Two of those sites are monkeypatch wrappers that `from_minimal` must keep reading.
     - The other fold, two wrappers over one private walk, adds exactly the indirection house style rejects.
   - Cost of keeping both: about eight duplicated lines of the chain lookup.
   - "Yes" would add `test_iterated_partition.py` and `test_partition_family.py` to the file list (nine files).
2. **Test placement: the two extra files (`helpers/laws.py` and the series test) and the regression tier.** Default: **yes**, as planned.
   - Alternative: drop the series test and check only the builds. But no build carries a lobe where the orders differ, so the test would have no teeth on base.

No earlier `plan_ok` verdict exists for this run, so `.refactor/decisions.md` has no notes to record.
