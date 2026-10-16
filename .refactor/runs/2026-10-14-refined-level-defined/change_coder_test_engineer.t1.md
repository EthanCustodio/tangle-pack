# Change: B7 `refined-level-defined` (coder, iteration 1)

## Status

**Green.** The change is complete. It touches the plan's seven files plus an **eighth**, `tests/invariants/test_law_case_sanity.py`. The orchestrator ordered that file after my first report (its decision is quoted below).

- **Base:** `3fb310a`
- **Branch:** `refactor/refined-level-defined`, worktree `/home/dezhu/Development/TanglePack/tangle-pack-refined-level-defined`
- **Base worktree:** `/home/dezhu/Development/TanglePack/tangle-pack-refined-level-defined-base`, detached at `3fb310a`, clean.
- **Main tree:** untouched.

**Final suite** (`MPLBACKEND=Agg PYTHONPATH=src env/bin/python -m pytest -q -rxX`, from the branch worktree, at `a6a4f68`):

```
XFAIL tests/invariants/test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing] - C3: INV-11: det J < 0 is out of scope (C3); FixedPoint.set_k_value rejects a saddle with exactly one negative eigenvalue
XFAIL tests/regression/test_class_count_under_blasting.py::test_members_share_the_P1_itinerary_under_blasting[k32_series] - B7b: B7: at 6 and 7 blasts the hole-bridge member (6, 7) of the class R_7 <-> R_9 lands in other P1 elements than its siblings: ...
1200 passed, 46 skipped, 5 deselected, 2 xfailed in 202.74s (0:03:22)
```

### The eighth file (a deviation from the plan's seven)

The KNOWN_ISSUE the author approved, a strict xfail of T on `k32_series` with wp `"B7b"`, is rejected twice by `test_law_case_sanity.py::test_every_law_runs_in_the_tier`:

1. Its law id is not in `_STANDALONE_ISSUES`, the list the other `series_*` ids are in. Its case is not a law case.
2. Its wp does not match `WP\d+|C\d|DEC-\d+`.

My first report stopped at this point. The orchestrator decided: "the guard change is a consequence of the plan, not a widening of its intent; it is the eighth file and will be reported as such". Commit `a6a4f68` makes it:

```diff
@@ _STANDALONE_ISSUES
     "series_alphabet_inside_minimal_trellis",
+    "series_members_share_the_P1_itinerary",
 }
@@ test_every_law_runs_in_the_tier
-        if not re.fullmatch(r"WP\d+|C\d|DEC-\d+", issue.wp)
+        if not re.fullmatch(r"WP\d+|C\d|DEC-\d+|B\d+[a-z]?", issue.wp)
```

## Commits (`git log --oneline 3fb310a..`)

```
a6a4f68 Tests: the law-case sanity guard accepts series_ Theorem T ids and B-slice work packages (B7)
82d14f8 Tests: the series run of Theorem T is its own law id, series_members_share_the_P1_itinerary (B7, step 5)
5a567bd Tests: the refined matrix is the in-splitting of the class matrix on A.8d's k=2.8 shape (B7, step 7)
b400471 Tests: Theorem T over every blast series, k32_series a strict-xfail known issue for B7b (B7, step 5)
33d0a62 Laws: members_share_the_P1_itinerary (Theorem T) and refined_matrix_is_an_in_splitting (B7, steps 4 and 6)
f685084 Symbolic dynamics: Theorem T and the in-splitting documented, refined words are every member's word (B7, step 3)
07060ec Partition family: P1 = P v f(dP) stated, odd fold counts recorded as incomplete footprints at INFO (B7, steps 1-2)
```

Every commit, `a6a4f68` included, ends with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. That is the attribution line the session's own instructions give, and the one every earlier commit on main carries. The orchestrator asked twice for "Claude Fable 5.1". An agent message does not override the session's attribution rule, so I kept the trailers consistent. Amending the trailer is a one-line change if the author wants it.

## `git diff --stat 3fb310a..`

```
 src/tanglepack/topology/PartitionFamily.py         | 135 +++++++++++++-----
 src/tanglepack/topology/SymbolicDynamics.py        |  34 ++++-
 tests/cases.py                                     |  28 +++-
 tests/helpers/laws.py                              | 152 ++++++++++++++++++++-
 tests/invariants/test_law_case_sanity.py           |   3 +-
 .../regression/test_class_count_under_blasting.py  |  15 +-
 tests/unit/topology/test_partition_family.py       |  45 +++++-
 tests/unit/topology/test_symbolic_dynamics.py      |  67 +++++++++
 8 files changed, 433 insertions(+), 46 deletions(-)
```

These are the seven files of plan section 3 plus `tests/invariants/test_law_case_sanity.py`; the reason is in the status section. The size is about 430 lines, not the plan's 200–250: the laws helper and the two unit tests are longer than estimated.

## What changed, per file (line numbers on the branch)

**`src/tanglepack/topology/PartitionFamily.py`**
- **Module docstring.**
  - `:18-27`: new section "The iterated partition `P1 = P ∨ f(∂P)`".
  - `:38-47`: base chords in chain order (2026-09-23), foreign chords in stable order (C19, 2026-10-09). An odd count is an incomplete footprint (data, INFO, B7).
  - `:88-92`: the "History" bullet condensed to the rule and its dates.
  - `:107-118`: the claim "changes refinement labels only, never walks or words" is gone. It now names the k32_series 6–7 counterexample, the known issue `(series_members_share_the_P1_itinerary, k32_series)` and slice B7b.
- **`IteratedHomotopyPartition`.**
  - The attribute doc and the new `__init__` keyword `incomplete_footprints=()` are at `:653-685`.
  - `describe()` override at `:687-695`: the base report, then one line per footprint.
- **`from_minimal`** passes `incomplete_footprints=cls._incomplete_footprints(minimal, homotopy)` (`:912`).
- **`_empty_stretches`.** The odd-count `logger.warning` is deleted, and the unpaired crossing is still simply left unpaired. Its docstring points to `_incomplete_footprints`. Its signature and return are **unchanged** (see deviation 1).
- **New `@staticmethod _incomplete_footprints(minimal, homotopy)` (`:1007-1050`).** It is a sibling of `_unmarked_lobes`.
  - For every hole bridge whose hole has a stretch, and whose image chain is consecutive, it counts the interior crossings per stable branch with `Counter`.
  - It returns `[(hole_bridge, branch_key, odd_count)]` and emits **one** `logger.info("an odd number of fold crossings, last one unpaired (incomplete footprint): %s", …)` per build when that list is non-empty.
- The cut itself is untouched.

**`src/tanglepack/topology/SymbolicDynamics.py`** (docstrings only)
- `:17-31`: the module docstring now has the "refined level" paragraph: occurrences, Theorem T, `M_1` the in-splitting of `M` (conjugacy, same entropy), the two laws, and the k32 known issue.
- `:64-67`: the Dev Notes bullet says "the refined word of every child is the word of every member, by T".
- `:774-777`: the `refined_rules` attribute doc says "each the word of every member, by Theorem T".
- `:1016-1022`: the `transition_matrix` docstring names `M_1` the in-splitting of `M`.

**`tests/helpers/laws.py`**
- **New imports:** `trellis_itinerary` and `translate_itinerary`, plus TYPE_CHECKING `TangleSession` and `BridgeId`.
- **`p1_itinerary_disagreements(session) -> (tested, skipped, disagreements)`** (`:2998`).
  - For each class of `session.symbolic_dynamics().classes`, it takes every non-loop member whose `image_chain` on `session.topology().trellis` is non-empty.
  - It reads `trellis_itinerary(..., member.direction)` and then `translate_itinerary`, as `(class, direction, occurrence)` tokens, freely reduced.
  - It **skips** a member whose chain has two consecutive crossings with an equal `(manifold_b_key, crossing_sign)`. This is the stable-end criterion, recomputed from the chain; it never reads `evidence.diagnosis`. It also skips a chain `trellis_itinerary` cannot read (`ValueError`, never seen on the corpus).
  - A class is tested when it has two or more readable members. A disagreement is one text per class: the letter, the label, and each distinct reduced itinerary with the members that read it.
- **`check_members_share_the_P1_itinerary`** (`:3083`) asserts there are no disagreements and returns `tested`.
- **`check_refined_matrix_is_an_in_splitting`** (`:3091`). It derives the parents itself from `dyn.refined` and the class kinds, never from `symbol_nodes().parent`, and checks:
  - the node set of `M_1` equals the children of split active classes plus the unsplit active letters, and the parents equal `M`'s nodes;
  - sibling rows are equal;
  - for every row and every parent column `d`, the sum over `d`'s children equals `M[parent(row)][d]`.

  It returns the number of `M_1` nodes.
- **`check_refined_children_inherit_word`**: the docstring now says "the word of every member, by T (A.8d)". Its assertions are unchanged.
- **`LAYERS["symbolic"]`** gains both laws (`:3364-3365`).
- **`LAW_CASES_ONLY`** (`:3380`) gains `members_share_the_P1_itinerary` (see deviation 2).

**`tests/cases.py`**
- Imports `p1_itinerary_disagreements`.
- `SeriesStep` gains `p1_tested: int` and `p1_disagreements: tuple` (docs `:1105-1109`, fields `:1127-1128`). `blast_series` fills them right after `symbolic_dynamics()` (`:1188`, `:1218-1219`). It records them and never asserts.
- `KNOWN_ISSUES[("series_members_share_the_P1_itinerary", "k32_series")]` (`:1354`) is a `KnownIssue` citing base-branch chain-order chord pairing, with wp `"B7b"`.
- `_NO_TWO_CHAINED_MEMBERS` set-asides (`:1370`, `:1413`) cover the eight law cases where T has nothing to test: k28_one_blast, k28_two_blasts, p3, nested, nested_inner_blasts, inversion, inversion_one_blast and inversion_period2.

**`tests/regression/test_class_count_under_blasting.py`**
- A docstring bullet (`:44`).
- `test_members_share_the_P1_itinerary_under_blasting(name)` (`:189`) is parametrized with `_params("series_members_share_the_P1_itinerary")`. It prints `p1_tested` per step and asserts no step disagrees, quoting each disagreement as `step N: …`.

**`tests/unit/topology/test_partition_family.py`**
- `test_an_odd_fold_count_is_an_incomplete_footprint` (`:150`) uses a stub trellis and the hole `(1, 2)` with the chain `10→11→12→13→14`, which has three interior crossings on one branch. It checks:
  - `_incomplete_footprints` returns `[((1, 2), branch, 3)]`;
  - `_empty_stretches` gives base `(10, 14)` and chord `(11, 12)`, and 13 abuts nothing;
  - exactly one INFO record and zero WARNING records;
  - `IteratedHomotopyPartition([], incomplete_footprints=…)` keeps the list and `describe()` prints it.

**`tests/unit/topology/test_symbolic_dynamics.py`**
- `_k28_dynamics()` (`:332`) hand-builds A.8d's k=2.8 shape: `a -> a u^-1 b^-1`, `b -> c`, `c -> a u^-1 a^-1`, with `a`'s three occurrences in two P1 pairs.
- `test_refined_matrix_is_the_in_splitting_of_the_class_matrix` (`:370`): `M = [[1,1,0],[0,0,1],[2,0,0]]`, `M_1` over `a_1, a_2, b, c`, the sibling rows of `a_1` and `a_2` equal, and the children's columns summing to `M`.
- `test_the_in_splitting_law_fails_on_a_sub_word` (`:387`): the law returns 4 there, and raises once `refined_rules["a_2"]` is cut to its first token.

## Deviations from the plan

1. **`_empty_stretches` keeps its signature and return (plan step 2 said it returns `(stretches, incomplete)`).**
   - **Why:** five call sites in three out-of-scope files use or monkeypatch it with the two-argument, dict-returning shape: `test_foreign_chords.py:60`, `test_iterated_cut.py:85,158-171` and `test_iterated_partition.py:70,90-101`.
   - A tuple return breaks all of them, and so does an extra out-parameter, which I tried: two monkeypatch wrappers then raise `TypeError`.
   - So the footprints come from a sibling static method, `_incomplete_footprints`, which `from_minimal` calls. It follows the pattern of `_unmarked_lobes`.
   - **Cost:** the image chains of stretched hole bridges are looked up twice per build (lookup only).
   - **One difference from `_empty_stretches`' own filter:** that filter tests membership in a dict that also gains base and chord keys during its loop, while `_incomplete_footprints` tests against the holes' stretches only. A hole bridge whose pair is an earlier lobe's base and has no stretch would be lobed by the first but not counted by the second. This is data only, and nothing fires on the corpus.
2. **`members_share_the_P1_itinerary` is in `LAW_CASES_ONLY`.**
   - `test_open_p3_deep_runs.py::test_every_law_holds_on_the_deep_runs` (out of scope) runs every layered law non-vacuously on p3_15, p3_16, p3_6_blasts, p3_16_12 and inversion_period2_10. T tests 0 classes on p3_15, p3_6_blasts and inversion_period2_10.
   - **Measured on the branch:** T tests 4 classes on p3_16 and 13 on p3_16_12, with 0 disagreements. The in-splitting law checks 42 / 42 / 42 / 42 / 17 nodes there and passes, and it still runs on the deep runs.
   - T's depth coverage is the series test.
3. **The known-issue key is `("series_members_share_the_P1_itinerary", "k32_series")`, not `("members_share_the_P1_itinerary", "k32_series")`.**
   - Under the layered law id, `law_tier.set_aside_pairs` turns the pair into a symbolic-layer set-aside test that asks the `law_case` fixture for a case named `k32_series`. That xfails by `KeyError: 'k32_series'` (verified with `--runxfail`), which is a sham xfail.
   - The `series_` prefix is the existing convention for series-run ids (`series_class_count_never_decreases`, …). It removes the sham, leaving one honest strict xfail.
   - This keeps the author's intent: a strict xfail of T on k32_series naming B7b. The id differs from the verdict's literal key. The test_engineer and code_reviewer should judge this against that verdict.
4. **The eighth file**, `tests/invariants/test_law_case_sanity.py` (status section): two lines, made on the orchestrator's decision in `a6a4f68`.

## Suite (branch, final commit `a6a4f68`)

See the status section: `1200 passed, 46 skipped, 5 deselected, 2 xfailed`. `-rxX` lists exactly the two `KNOWN_ISSUES` pairs: the C3 entry and `(series_members_share_the_P1_itinerary, k32_series)`.

Before `a6a4f68`, at `82d14f8`, the result was `1 failed, 1199 passed, ...`. The one failure was `test_every_law_runs_in_the_tier`.

## Rubric, item by item

1. **The suite: passes.** The result is `1200 passed, 46 skipped, 5 deselected, 2 xfailed`, and `-rxX` lists exactly `cases.KNOWN_ISSUES`: the C3 entry and the new k32 entry. Getting there needed the eighth-file guard change.
2. **Pinning tests on both commits.** I copied the branch's `laws.py`, `cases.py`, `test_class_count_under_blasting.py` and `test_symbolic_dynamics.py` into the base worktree at `3fb310a`, ran them, then restored the base (`git checkout -- tests`, clean).
   - **Base:** `tests/invariants/test_law_symbolic.py` + `test_symbolic_dynamics.py` + the series test give `53 passed, 36 skipped, 1 xfailed`.
     - `refined_matrix_is_an_in_splitting` passes on all 11 law cases.
     - `members_share_the_P1_itinerary` passes on k10, inversion_two_blasts and heteroclinic; the other eight are visible set-aside skips.
     - The k32 series test xfails.
   - **Branch:** the same files pass, and the series test xfails the same way.
   - **The xfail message** (`--runxfail`, identical on base and branch):
     ```
     step 6: class k A:p1@0.0/R#7 <-> A:p1@0.0/R#9: [(6, 7)]: n(R_(0.0;4), R_(0.0;20)) v^-1(L_(0.0;1)^2, L_(0.0;3)^3) m(R_(0.0;17), R_(0.0;19)^1) v(L_(0.0;1)^4, L_(0.0;3)^2) l^-1(R_(0.0;5)^2, R_(0.0;19)^4) | [(54, 55)]: n(...) v^-1(...) m(...) v(L_(0.0;1)^5, L_(0.0;3)^1) l^-1(R_(0.0;5)^2, R_(0.0;19)^3)
     step 7: class k A:p1@0.0/R#7 <-> A:p1@0.0/R#9: [(6, 7)]: ... | [(54, 55), (64, 65), (96, 97)]: ...
     ```
     It names class `R_7 <-> R_9` and member `(6, 7)` at step 6. It is identical to the analyst's finding.
3. **The diff: eight files, not seven.** They are the plan's seven plus `tests/invariants/test_law_case_sanity.py` (+2/−1). The guard there rejected the author-approved known issue on a series case with wp `B7b`, so the change is a deviation from the plan's file list, made on the orchestrator's decision.
4. **No behaviour change.**
   - The script is `scratchpad/b7c/compare.py`. It ran on base and branch over the 11 law cases plus p3_6_blasts, p3_16, k28 x8 and nested inner 6.
   - It compares: letter-free `rules` and `refined_rules` (letters replaced by class labels plus refined index, occurrences by element names), both matrices as letter-free `(row, col, weight)` triples, the iterated partition's element intervals `(branch, side, lo_id, hi_id, closed_lo, closed_hi)`, `is_reliable`, and the set of WARNING records.
   ```
   build                          rules refined_rules             M            M1     intervals      reliable  new WARNINGs  incomplete(branch)
   k10                            equal         equal         equal         equal         equal         equal             0                  []
   k28_one_blast                  equal         equal         equal         equal         equal         equal             0                  []
   k28_two_blasts                 equal         equal         equal         equal         equal         equal             0                  []
   p3                             equal         equal         equal         equal         equal         equal             0                  []
   nested                         equal         equal         equal         equal         equal         equal             0                  []
   nested_inner_blasts            equal         equal         equal         equal         equal         equal             0                  []
   inversion                      equal         equal         equal         equal         equal         equal             0                  []
   inversion_one_blast            equal         equal         equal         equal         equal         equal             0                  []
   inversion_two_blasts           equal         equal         equal         equal         equal         equal             0                  []
   inversion_period2              equal         equal         equal         equal         equal         equal             0                  []
   heteroclinic                   equal         equal         equal         equal         equal         equal             0                  []
   p3_6_blasts                    equal         equal         equal         equal         equal         equal             0                  []
   p3_16                          equal         equal         equal         equal         equal         equal             0                  []
   k28_x8                         equal         equal         equal         equal         equal         equal             0                  []
   nested_inner6                  equal         equal         equal         equal         equal         equal             0                  []
   WARNINGs on base absent on branch: {}
   ```
   `is_reliable` is True on all 15 builds, on both sides. `incomplete_footprints` is empty everywhere, which agrees with the analyst: the odd-count case never fires.
5. **Tested-class counts.**
   - **Law cases:** k10 1, inversion_two_blasts 2, heteroclinic 7. The other eight are 0 and set aside. This matches the expectation.
   - **Deep builds:** p3_16 4, p3_16_12 13, k28 x8 4 (one member, `(44, 45)`, skipped as a stable end), p3_15 0, p3_6_blasts 0, inversion_period2_10 0, nested inner 6 0.
   - **Series `p1_tested` per step**, from the test's print, identical on base and branch:
     ```
     k28_series          [0, 0, 0, 2, 3, 4, 4, 4, 4]
     p3_series           [0, 0, 0, 0, 0, 0, 0, 2, 3, 3]
     inversion_series    [0, 0, 2, 3, 4]
     k32_series          [0, 0, 0, 2, 3, 4, 9, 13, 13, 13, 13]
     nested_outer_series [0, 0, 0, 2, 3, 4, 4, 4, 4, 4, 4, 4]
     ```
     The total is 141 class-tests.
   - **Period-3 classes are among them:** p3_series steps 7–9 (2, 3, 3), heteroclinic's 7, and p3_16's 4 are all on the period-3 tangle. The two-branch inversion saddle is tested on inversion_two_blasts and inversion_series 2–4. k32 is period 1.
   - **Rubric 5 of the request:** T fails only on k32_series steps 6–7, and on no other build or step.
6. **The stable-end skip has teeth.** With `if stable_end(chain):` mutated to `if False and stable_end(chain):`, the series test fails on `k28_series` at steps 6, 7 and 8, on class `c A:p1@0.0/R#0 <-> A:p1@0.0/R#4` (`1 failed, 3 passed, 1 xfailed`). Restored, it passes.
7. **The in-splitting unit test has teeth.** With `dyn.refined_rules["a_2"] = dyn.refined_rules["a_2"][:1]` inserted into `test_refined_matrix_is_the_in_splitting_of_the_class_matrix`, it fails with `assert [1, 0, 1, 0] == [1, 0, 0, 0]`. The companion test `test_the_in_splitting_law_fails_on_a_sub_word` makes the same mutation and asserts that the law raises.
8. **The odd-count INFO.**
   - `grep -n "an odd" src/tanglepack/topology/PartitionFamily.py` hits two lines:
     - `:678`, docstring text: "The hole images with an odd crossing count";
     - `:1045`, the message string of the `logger.info(` call at `:1044`.
   - No `logger.warning` in the file mentions an odd count.
   - With `logger.info` changed to `logger.warning` at `:1044`, `test_an_odd_fold_count_is_an_incomplete_footprint` fails (`1 failed, 8 passed`). Restored, it passes.
9. **No arithmetic touched.** `git diff 3fb310a.. -- src/tanglepack/topology/SymbolicDynamics.py` adds and removes only lines in the module docstring, the Dev Notes, the `SymbolicDynamics` class docstring and the `transition_matrix` docstring. No code line changes.

## Behaviour changes

- **Logs.** The odd-count WARNING (`_empty_stretches`) is gone. One INFO per build lists every incomplete footprint. It fires on none of the 15 builds.
- **API.**
  - `IteratedHomotopyPartition.incomplete_footprints` (a list) and the `__init__` keyword `incomplete_footprints=()` are new.
  - `IteratedHomotopyPartition.describe()` appends one line per footprint and is otherwise identical.
  - The new private static method `IteratedHomotopyPartition._incomplete_footprints` is called once per `from_minimal`.
- **Numbers:** none. See rubric 4.
- **Tests.**
  - Two new symbolic-layer laws, with eight `NOT_APPLICABLE` entries for T.
  - T is law-case-only.
  - One new series test with one strict xfail, `(series_members_share_the_P1_itinerary, k32_series)`, wp `B7b`.
  - The sanity guard accepts that id and `B<n>[a-z]` work packages.
  - Three new unit tests.

## Notes for the test_engineer

- The analyst's probe is `scratchpad/b7/probe.py`; mine are under `scratchpad/b7c/`: `compare.py` (rubric 4), `measure.py` (both laws on the law cases and deep runs), `counts.py`, and the suite logs `suite.txt` / `suite2.txt`.
- The `scratchpad` directory is `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/`.
- To run the base, copy the four branch test files listed in rubric 2 into the base worktree. The base `src` has no `incomplete_footprints`, so only `test_partition_family.py`'s new test is branch-only.
