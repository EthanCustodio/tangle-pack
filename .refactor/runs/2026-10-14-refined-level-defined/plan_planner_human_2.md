# Plan — B7 `refined-level-defined` (base `3fb310a`)

## 1. The answer

**Is the refined word a theorem? Almost everywhere: Theorem T holds on all 15 builds and on 45 of the 47 series steps.** Theorem T says every member of a class has the same `P1`-level itinerary.
- **Where it was tested.** On the law cases T applies to k10 (1 class), inversion_two_blasts (2) and heteroclinic (7), and it holds with 0 disagreements. It also holds on p3_16 (4 classes) and on 141 class-tests across the five blast series. Period-3 classes (heteroclinic, p3_16, p3_series 7–9) and the two-branch inversion saddle are included.
- **The in-splitting identities hold on all 62 builds and series steps.** Sibling rows of `M_1` are equal, and the children's columns sum to the parent's.
- **One counterexample: `k32_series` blasts 6 and 7.** It is described in the next paragraph and is the one decision this plan needs from you.
- **Not a counterexample: k28 at 8 blasts.** One member's image there runs past the computed stable branch, the existing "stable end" case (C11). This is incomplete evidence, and the law skips such members.

**The k32 disagreement, in plain words.** On `k32_series` at 6 and 7 blasts:
- One class (`R_7 ↔ R_9`) has four member bridges. Their images all spell the same class-level word, but one member, the hole bridge `(6, 7)`, lands in differently named `P1` sub-elements from its three siblings.
- The cause is how `P1` cuts the stable manifold under a lobe's image. That image folds back across its own stable branch twice, in a spiral. Today those fold crossings are paired off in the order the curve visits them (your 2026-09-23 rule for the base branch). The analyst's point-in-polygon test on that lobe shows the stretches the image really covers are paired in the order they sit along the stable manifold. That is the same rule `P1` already uses on every other branch (C19).
- So on this one lobe `P1` marks the wrong stretch empty, and a sibling's image legitimately crosses it.
- The analyst tried stable-order pairing on the base branch as well, as a scratch patch only. T then holds on all 62 builds and steps, and `M`, `M_1`, class counts and `is_reliable` are unchanged. What changes: inert splits (−1 on 15 builds and steps, all outside both matrices) and iterated element names where the two orders differ.

## 2. Recommendation

**Make B7 as specified. Register the k32 counterexample as a known issue that the next slice fixes.**
- B7 gives the refined level its definition and its two laws without changing a number.
- The pairing change goes into its own slice, **B7b `base-chords-in-stable-order`**, run next. That change reverses one of your rules and renames iterated elements.
- The strict xfail on `(members_share_the_P1_itinerary, k32_series)` becomes B7b's failing test, and it flips to a pass when B7b lands.

This follows Q15's spirit: B7 stops short of deciding the premise, and the decision is yours (question 1).

## 3. Scope

Slug **`refined-level-defined`**. The branch is `refactor/refined-level-defined` in `../tangle-pack-refined-level-defined`, made from `3fb310a`. The base worktree is `../tangle-pack-refined-level-defined-base`. The uncommitted `docs/`, `.refactor/` and `CLAUDE.md` changes in your tree are not the branch's.

Files (seven):

1. `src/tanglepack/topology/PartitionFamily.py`
2. `src/tanglepack/topology/SymbolicDynamics.py`
3. `tests/helpers/laws.py`
4. `tests/cases.py`: the `NOT_APPLICABLE` set-asides, a `SeriesStep` field, and the one `KNOWN_ISSUES` entry. This is the seventh file; the request listed six.
5. `tests/regression/test_class_count_under_blasting.py`
6. `tests/unit/topology/test_symbolic_dynamics.py`
7. `tests/unit/topology/test_partition_family.py`

Size: about 200–250 changed lines, mostly docstrings and the two laws.

## 4. Steps

1. **`PartitionFamily.py`, module docstring.**
   - State `P1 = P ∨ f(∂P)` and the empty-stretch rule as in A.8d: holes, lobe bases and chords; base chords paired in chain order (2026-09-23), foreign chords in stable order (C19).
   - Condense the Dev Notes' three-rewrite history (`:94-106`) to the rule and its dates.
   - Drop the claim that chain order "changes refinement labels only, never walks or words" as a general fact. Point to B7b and the k32 known issue instead.
   - The cut itself is not touched.

2. **`PartitionFamily.py`, `incomplete_footprints`.**
   - `_empty_stretches` (`:892`) returns `(stretches, incomplete)`, where `incomplete` is a `list[tuple[BridgeId, ManifoldKey, int]]` (hole bridge, branch, odd count).
   - `from_minimal` (`:748`) passes it to `IteratedHomotopyPartition.__init__` (`:645`). The `__init__` gains the keyword `incomplete_footprints=()` and the docstring line.
   - The WARNING at `:959-966` goes. After the loop, one `logger.info` per build names every entry when there is any:

     ```python
     if incomplete:
         logger.info("an odd number of fold crossings, last one unpaired (incomplete footprint): %s", incomplete)
     ```

   - `IteratedHomotopyPartition.describe()` overrides the base: `super().describe()` plus one line per footprint.
   - This is data, never a reliability failure (Q14).

3. **`SymbolicDynamics.py`, documentation only.**
   - The module docstring states Theorem T and the in-splitting (A.8d): the refined symbols are occurrences (the `P1` pair a token lands on); under T, `M_1` is the in-splitting of `M`, a conjugacy with the same entropy.
   - The Dev Notes bullet (`:47-58`) and the `refined_rules` docstring (`:758-760`) change "every refined child inherits the parent's word" to "the word of every member, by T".
   - The `transition_matrix(refined=True)` docstring names it `M_1`, the in-splitting of `M`.
   - No arithmetic changes: `refine`, `_refined_symbol`, `match_members` and the `refined_rules` construction (`:796-809`) are untouched.

4. **`laws.py`: Theorem T.**
   - A helper `p1_itinerary_disagreements(session) -> (tested, skipped, disagreements)`. For every class of `session.symbolic_dynamics().classes`, it collects each non-loop member whose `_image_chain` on `topology().trellis` is non-empty, using `trellis_itinerary(..., member.direction)` then `translate_itinerary`, as `(bridge_class, direction, occurrence)` tuples.
   - It skips a member whose chain has two consecutive crossings on one stable branch with one `crossing_sign` (the stable-end criterion). The law recomputes this criterion from the chain and never reads `evidence.diagnosis`.
   - It freely reduces each sequence. A class with two or more remaining members is tested, and a class whose sequences differ is a disagreement, recorded with its members and their sequences.
   - The law `members_share_the_P1_itinerary` (layer `symbolic`) asserts there is no disagreement, and its message quotes the class and both itineraries. It returns `tested`.
   - The eight law cases with nothing to test get `NOT_APPLICABLE` entries in `cases.py`: "no class has two members with registered image chains".

5. **The series run.**
   - `SeriesStep` gains `p1_tested: int` and `p1_disagreements: tuple`, filled in `blast_series` by the same helper. It records and never asserts, so one bad step does not kill the other series tests.
   - A new `test_members_share_the_P1_itinerary_under_blasting(name)` asserts that no step has disagreements and prints the per-step `p1_tested`.
   - `KNOWN_ISSUES[("members_share_the_P1_itinerary", "k32_series")]` is a `KnownIssue` citing the base-branch chain-order chord pairing, with wp `"B7b"`. It runs as a strict xfail.

6. **`laws.py`: in-splitting.**
   - The law `refined_matrix_is_an_in_splitting` (layer `symbolic`) checks three identities:
     - sibling rows of `M_1` are equal;
     - for every row and every parent column `d`, `Σ_k M_1[r][d_k] = M[parent(r)][d]`. This is the identity that catches a token whose occurrence matched no child and was silently dropped (`transition_matrix` `:1012`);
     - `M_1`'s node set is the children of the split active classes plus the unsplit active letters.
   - It is applicable wherever `M` is non-empty. A case with no active class gets a set-aside.
   - `check_refined_children_inherit_word` keeps its name and assertions; its docstring says "the word of every member, by T". It stays because it must keep failing if the inheritance breaks.
   - Step 6 of the request is already done: `check_k10_is_the_full_shift` cites A.8d.

7. **Unit tests.**
   - `test_symbolic_dynamics.py`: an in-splitting test on a hand-built table of A.8d's k=2.8 shape. It checks `M = [[1,1,0],[0,0,1],[2,0,0]]`, with `a` split into `a_1` and `a_2`, so that `M_1`'s sibling rows are equal and its columns sum to `M`.
   - `test_partition_family.py`: `_empty_stretches` on a stub trellis whose hole-bridge chain has three interior crossings on one branch. It checks that `incomplete_footprints == [(bridge, branch, 3)]`, that the INFO is logged and that no WARNING is.

## 5. Behaviour

**None in any number.** Words, refined children, `refined_rules`, both matrices, iterated elements and `is_reliable` are the same on base and branch.
- **Logs:** the odd-count WARNING becomes one INFO per build. It fires on none of the 62 builds and steps, so `incomplete_footprints` is empty on the whole corpus.
- **API:** `IteratedHomotopyPartition.incomplete_footprints` is new, and `describe()` lists it.
- **Tests:** one new strict xfail, `(members_share_the_P1_itinerary, k32_series)`, citing B7b.

## 6. Not doing

- **The base-branch pairing change.** That is B7b, or this slice if question 1 takes the alternative.
- **A second level `P2` (Q16).**
- **Any other change to `from_minimal`'s cuts.**
- **"One element at a bridge end"** (B9).
- **The refinement arithmetic**: `refine`, `_refined_symbol`, `match_members`.
- **Inert-only splitting.** `refine` also counts occurrences that only inert words name (`:625`), so a child can have no incoming edge in `M_1`. No instance was measured.
- **Docs.** `CLAUDE.md` and the alignment plan are left to B12.

## 7. Rubric

1. **The suite.** `MPLBACKEND=Agg env/bin/python -m pytest -q` passes on the branch. `-rxX` lists only `cases.KNOWN_ISSUES`, which now includes the k32 entry.
2. **Pinning tests on both commits.** Copy the branch's `laws.py` and `cases.py` into the base worktree at `3fb310a`.
   - On the base, `refined_matrix_is_an_in_splitting` passes on every applicable law case.
   - On the base, `members_share_the_P1_itinerary` passes on k10, inversion_two_blasts and heteroclinic.
   - On the base, the series test xfails on k32_series. Its message names class `R_7 ↔ R_9` and member `(6, 7)` at step 6.
   - All of these hold on the branch too.
3. **The diff.** `git diff --stat 3fb310a..` touches exactly the seven files of section 3.
4. **No behaviour change.** One script runs base and branch and the report quotes its output. It covers the 11 law cases, `p3_6_blasts`, `p3_16`, k28 at 8 blasts and nested with 6 inner blasts, and it shows equal on both:
   - letter-free `rules` and `refined_rules`;
   - `transition_matrix(refined=False)` and `transition_matrix(refined=True)`;
   - the iterated partition's element intervals;
   - `is_reliable`.

   No WARNING record appears on the branch that is absent on the base.
5. **Tested-class counts.** The report prints the per-build tested-class counts. Expected: k10 1, inversion_two_blasts 2, heteroclinic 7. It also prints the per-step `p1_tested` of every series, and confirms whether period-3 classes were among them (expected: yes, p3_series 7–9).
6. **The stable-end skip has teeth.** With the skip removed, `members_share_the_P1_itinerary_under_blasting` fails on `k28_series`, at steps 6–8. With it in place, it passes.
7. **The in-splitting unit test has teeth.** Replacing `refined_rules["a_2"]` by its first token makes the test fail.
8. **The odd-count INFO.**
   - `grep -n "an odd" src/tanglepack/topology/PartitionFamily.py` hits a `logger.info`, and no `logger.warning` call in that file mentions an odd count.
   - The three-crossing unit test passes and fails if the INFO is turned back into a WARNING.
9. **No arithmetic touched.** `git diff 3fb310a.. -- src/tanglepack/topology/SymbolicDynamics.py` changes only docstrings and comments, with no line inside a function body outside a docstring.

## 8. Questions

1. **How should the k32 counterexample be handled?**

   **Default (yes = this):** option (ii). Keep chain order in B7. Register `(members_share_the_P1_itinerary, k32_series)` as a strict-xfail known issue naming B7b. Write B7b's request (switch base-branch chords to stable order, with a base-branch geometric test like `test_foreign_chords.py`) as the next slice.

   **Alternative, option (i):** widen B7 to switch the base branch to stable order now. Section 4 step 1 drops its `branch_key != base_branch` test, and there is no known issue. Section 5 then has a behaviour change: inert splits −1 on 15 builds and steps, and iterated element names renamed wherever the orders differ. That reverses your 2026-09-23 chain-order rule, and the evidence for it is one lobe plus the corpus-wide patch run, not a proof.

## Where the coder starts

`PartitionFamily.py:892-975` (`_empty_stretches`) and `:633-760` (`IteratedHomotopyPartition`); `laws.py:2977-3100` and `LAYERS["symbolic"]` `:3218`; `cases.py:1085` (`SeriesStep`), `:1151` (`blast_series`), `:1336` (`KNOWN_ISSUES`) and `:1356` (`NOT_APPLICABLE`). The analyst's probe is `scratchpad/b7/probe.py` under `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/`; it is the reference for the law's itinerary and stable-end criterion.
