# Test evidence — B7 `refined-level-defined` (test_engineer, traversal 1)

## Verdict

**Every rubric item (1–9) holds, and every break-it check I ran made the guarding test fail.** The coder's four deviations are sound; the corpus measurements are below.

- **New tests:** none added. I found no gap that a test could close without pinning a case-specific count, which `tests/conftest.py` policy forbids (see "Gaps").
- **Branch:** `refactor/refined-level-defined` at `a6a4f68`, worktree `/home/dezhu/Development/TanglePack/tangle-pack-refined-level-defined`.
- **Base:** detached worktree `/home/dezhu/Development/TanglePack/tangle-pack-refined-level-defined-base` at `3fb310a`.
- **Final state:** both worktrees are clean at those heads, and the main tree was not touched.
- **Python:** `/home/dezhu/Development/TanglePack/tangle-pack/env/bin/python`, run with `MPLBACKEND=Agg`.
- **Import check:** every run asserted `tanglepack.__file__` (pytest with `PYTHONPATH=src`; scripts with `-I` and an explicit `sys.path`). Base runs printed `.../tangle-pack-refined-level-defined-base/src/tanglepack/__init__.py` and branch runs printed `.../tangle-pack-refined-level-defined/src/tanglepack/__init__.py`.
- **Scripts:** all in `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b7t/`: `fp.py`, `cmp.py`, `laws_run.py`, `series_probe.py`, `skips.py`, `astcmp.py`, with their outputs. I wrote them myself; none is the coder's `compare.py`.

## Rubric 1 — the suite on the branch

```
cd tangle-pack-refined-level-defined && PYTHONPATH=src MPLBACKEND=Agg ../tangle-pack/env/bin/python -m pytest -q -rxX -p no:cacheprovider
XFAIL tests/invariants/test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing] - C3: INV-11: det J < 0 is out of scope (C3); ...
XFAIL tests/regression/test_class_count_under_blasting.py::test_members_share_the_P1_itinerary_under_blasting[k32_series] - B7b: B7: at 6 and 7 blasts the hole-bridge member (6, 7) of the class R_7 <-> R_9 lands in other P1 elements than its siblings: ...
1200 passed, 46 skipped, 5 deselected, 2 xfailed in 208.47s (0:03:28)
```

- No XPASS and no failure.
- `-rxX` lists exactly the two `cases.KNOWN_ISSUES` entries: the existing C3 entry and the new `("series_members_share_the_P1_itinerary", "k32_series")` with wp `B7b`. The id differs from the verdict's literal key; see deviation 3.

## Rubric 2 — the pinning tests on base and branch

**Setup.** I copied five branch files into the base worktree:
- `tests/helpers/laws.py`
- `tests/cases.py`
- `tests/regression/test_class_count_under_blasting.py`
- `tests/unit/topology/test_symbolic_dynamics.py`
- `tests/invariants/test_law_case_sanity.py`

After the runs I restored the base with `git checkout -- tests`; `git status --short` printed nothing.

**pytest on the base:**
```
PYTHONPATH=src ... pytest -q -rxXs -s tests/invariants/test_law_symbolic.py tests/unit/topology/test_symbolic_dynamics.py \
  tests/regression/test_class_count_under_blasting.py::test_members_share_the_P1_itinerary_under_blasting tests/invariants/test_law_case_sanity.py
k28_series p1_tested per step: [0, 0, 0, 2, 3, 4, 4, 4, 4]
p3_series p1_tested per step: [0, 0, 0, 0, 0, 0, 0, 2, 3, 3]
inversion_series p1_tested per step: [0, 0, 2, 3, 4]
k32_series p1_tested per step: [0, 0, 0, 2, 3, 4, 9, 13, 13, 13, 13]
nested_outer_series p1_tested per step: [0, 0, 0, 2, 3, 4, 4, 4, 4, 4, 4, 4]
XFAIL ...test_members_share_the_P1_itinerary_under_blasting[k32_series] - B7b: ...
XFAIL ...test_orientation_reversing_case_builds[orientation_reversing] - C3: ...
SKIPPED x8: not applicable: members_share_the_P1_itinerary on {inversion, inversion_one_blast, inversion_period2, k28_one_blast, k28_two_blasts, nested, nested_inner_blasts, p3}: no class has two members with registered image chains
65 passed, 36 skipped, 2 xfailed in 81.43s
```
The eight not-applicable entries are visible set-aside skips. Their set-aside tests also run the law and assert it still checks 0 items.

**Each law per law case, called directly** (`laws_run.py`; the return value is the count of items checked). The output was identical line for line on base and branch:
```
k10                    members_share_the_P1_itinerary=1  refined_matrix_is_an_in_splitting=2
k28_one_blast          members_share_the_P1_itinerary=0  refined_matrix_is_an_in_splitting=1
k28_two_blasts         members_share_the_P1_itinerary=0  refined_matrix_is_an_in_splitting=4
p3                     members_share_the_P1_itinerary=0  refined_matrix_is_an_in_splitting=3
nested                 members_share_the_P1_itinerary=0  refined_matrix_is_an_in_splitting=7
nested_inner_blasts    members_share_the_P1_itinerary=0  refined_matrix_is_an_in_splitting=7
inversion              members_share_the_P1_itinerary=0  refined_matrix_is_an_in_splitting=2
inversion_one_blast    members_share_the_P1_itinerary=0  refined_matrix_is_an_in_splitting=7
inversion_two_blasts   members_share_the_P1_itinerary=2  refined_matrix_is_an_in_splitting=7
inversion_period2      members_share_the_P1_itinerary=0  refined_matrix_is_an_in_splitting=4
heteroclinic           members_share_the_P1_itinerary=7  refined_matrix_is_an_in_splitting=74
```

- **In-splitting:** passes on all 11 law cases on both commits.
- **T:** passes, non-vacuously, on k10, inversion_two_blasts and heteroclinic. It checks 0 classes on the other eight, which are exactly the eight `NOT_APPLICABLE` entries.

**The k32 message.** I ran it on each commit with `--runxfail`:
```
pytest -q --runxfail "tests/regression/test_class_count_under_blasting.py::test_members_share_the_P1_itinerary_under_blasting[k32_series]"
```
Both commits print `1 failed`, and the output is byte-identical apart from the timing line:
```
E   AssertionError: step 6: class k A:p1@0.0/R#7 <-> A:p1@0.0/R#9: [(6, 7)]: n(R_(0.0;4), R_(0.0;20)) v^-1(L_(0.0;1)^2, L_(0.0;3)^3) m(R_(0.0;17), R_(0.0;19)^1) v(L_(0.0;1)^4, L_(0.0;3)^2) l^-1(R_(0.0;5)^2, R_(0.0;19)^4) | [(54, 55)]: n(...) v^-1(...) m(...) v(L_(0.0;1)^5, L_(0.0;3)^1) l^-1(R_(0.0;5)^2, R_(0.0;19)^3)
E     step 7: class k A:p1@0.0/R#7 <-> A:p1@0.0/R#9: [(6, 7)]: ... | [(54, 55), (64, 65), (96, 97)]: ...
```

- It names class `R#7 <-> R#9` and member `(6, 7)` at step 6, as the rubric requires. Only steps 6 and 7 disagree; steps 8–10 test 13 classes each with no disagreement.
- **Cosmetic, not a rubric failure.** `R#7` / `R#9` is the positional `ElementRef.label` (`element_id`). In the `ElementNaming` convention (subscript = `element_id + 1`) the elements are `R_(0.0;8)` and `R_(0.0;10)`. The KnownIssue text, the plan and the PartitionFamily Dev Note write "R_7 <-> R_9", which mixes the two conventions. I inferred this from the convention in CLAUDE.md and did not resolve the names through `ElementNaming`.

## Rubric 3 — the diff

```
git diff --stat 3fb310a..
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

- There are **eight files, not seven**. The eighth is `tests/invariants/test_law_case_sanity.py` (+2/−1, commit `a6a4f68`).
- That commit makes two changes:
  - it adds `"series_members_share_the_P1_itinerary"` to `_STANDALONE_ISSUES`;
  - it widens the work-package regex from `WP\d+|C\d|DEC-\d+` to also accept `B\d+[a-z]?`.
- Both are required for the approved `B7b` known issue to pass the sanity guard.
- The size is about 430 lines, against the plan's estimate of 200–250.

## Rubric 4 — no behaviour change (my own script)

**Method.** `fp.py ROOT BRANCH_ROOT` runs once per tree and writes JSON lines; `cmp.py` compares the two runs. On both trees the Theorem T helper is loaded from the branch's `laws.py` file; it uses library API only.

**Identities, letter-free.**
- A class is the sorted pair of its two homotopy elements' `(branch, side, lo_id, hi_id, closed_lo, closed_hi)`.
- An occurrence is the pair of iterated-element identities.
- A refined child is its parent class plus its occurrence.
- Matrices are sets of `(row identity, col identity, weight)` triples.
- WARNING records are a multiset of full messages, captured by a handler on the `tanglepack` logger and taken per build.

```
build                          rules refined_rules             M            M1     intervals      reliable  warn b/br newWARN tested b/br  skip dis filt= odd incompl #cls #ref
k10                            equal         equal         equal         equal         equal         equal    0/0           0     1/1         0   0  True   0       0    2    3
k28_one_blast                  equal         equal         equal         equal         equal         equal    0/0           0     0/0         0   0  True   0       0    3    3
k28_two_blasts                 equal         equal         equal         equal         equal         equal    0/0           0     0/0         0   0  True   0       0    4    5
p3                             equal         equal         equal         equal         equal         equal    0/0           0     0/0         0   0  True   0       0    5    5
nested                         equal         equal         equal         equal         equal         equal    0/0           0     0/0         0   0  True   0       0    9   10
nested_inner_blasts            equal         equal         equal         equal         equal         equal    0/0           0     0/0         0   0  True   0       0    9   10
inversion                      equal         equal         equal         equal         equal         equal    0/0           0     0/0         0   0  True   0       0    4    4
inversion_one_blast            equal         equal         equal         equal         equal         equal    0/0           0     0/0         0   0  True   0       0    7    8
inversion_two_blasts           equal         equal         equal         equal         equal         equal    0/0           0     2/2         0   0  True   0       0    8    9
inversion_period2              equal         equal         equal         equal         equal         equal    0/0           0     0/0         2   0  True   0       0    6    6
heteroclinic                   equal         equal         equal         equal         equal         equal    0/0           0     7/7         0   0  True   0       0   75   80
p3_6_blasts                    equal         equal         equal         equal         equal         equal    0/0           0     0/0         0   0  True   0       0   43   45
p3_16                          equal         equal         equal         equal         equal         equal    0/0           0     4/4         0   0  True   0       0   43   45
k28_x8                         equal         equal         equal         equal         equal         equal    4/4           0     4/4         1   0  True   0       0    4    5
nested_inner6                  equal         equal         equal         equal         equal         equal    0/0           0     0/0         0   0  True   0       0   47   50
reliable: all 15 True on both
base-only WARNINGs: none on any build
```

- All six compared quantities are equal on all 15 builds: letter-free rules, refined rules, `M`, `M_1`, iterated element intervals and `is_reliable`.
- No WARNING appears on the branch that is absent on the base.
- The only WARNINGs anywhere are the same 4 on k28 x8, identical on both commits and all pre-existing:
  - two "detections are not crossings of N" (one with 3 of 87 detections, one with 5 of 243);
  - one "82 same-stability polyline crossing(s)";
  - one cross-side itinerary pair.
- **On base, the odd-count WARNING never fired**: `odd` is 0 everywhere. So deleting it removes nothing that was observable.
- **On branch, `incomplete_footprints` is empty on all 15 builds.**

## Rubric 5 — tested-class counts

**Per build** (`tested`, the same on base and branch):
- Law cases: k10 **1**, inversion_two_blasts **2**, heteroclinic **7**. This matches the expectation; the other eight law cases are 0.
- Deep builds: p3_16 4, k28 x8 4, p3_6_blasts 0, nested inner 6 0.

**Every skip is a stable end.** I wrapped `laws.trellis_itinerary` to count `ValueError`s (`skips.py`, `series_probe.py`):
- inversion_period2 skips `(0, 5)` and `(1, 9)`;
- k28 x8 skips `(44, 45)`;
- k28_series skips 1 member at each of steps 6–8.

None of these is a `ValueError`: the law's `except ValueError` skip never fires on the corpus.

**Series `p1_tested` per step** (the test's print on base, re-measured on branch by `series_probe.py`):
```
k28_series          [0, 0, 0, 2, 3, 4, 4, 4, 4]           skipped [0,0,0,0,0,0,1,1,1]
p3_series           [0, 0, 0, 0, 0, 0, 0, 2, 3, 3]
inversion_series    [0, 0, 2, 3, 4]
k32_series          [0, 0, 0, 2, 3, 4, 9, 13, 13, 13, 13]  disagreements [0,0,0,0,0,0,1,1,0,0,0]
nested_outer_series [0, 0, 0, 2, 3, 4, 4, 4, 4, 4, 4, 4]
```

- The steps total 141 class-tests.
- **Period-3 classes are covered:** p3_series steps 7–9 (2, 3, 3), heteroclinic (7) and p3_16 (4).
- **The two-branch inversion saddle is covered:** inversion_two_blasts and inversion_series steps 2–4.

## Rubrics 6–8 — break-it checks (each run on the branch, then restored with `git checkout`)

| # | Break | Command | Result | Restored |
|---|---|---|---|---|
| 6 | `laws.py:3060` `if stable_end(chain):` → `if False and stable_end(chain):` | `pytest -q -s ...::test_members_share_the_P1_itinerary_under_blasting` | **`1 failed, 3 passed, 1 xfailed`**: `FAILED [k28_series]`, at `step 6`, `step 7`, `step 8`, class `c A:p1@0.0/R#0 <-> A:p1@0.0/R#4`. k32 still xfails, and the other three series pass | `git checkout -- tests/helpers/laws.py`, clean |
| 7a (the rubric's) | Insert `dyn.refined_rules["a_2"] = dyn.refined_rules["a_2"][:1]` after `dyn = _k28_dynamics()` in `test_refined_matrix_is_the_in_splitting_of_the_class_matrix` | `pytest -q tests/unit/topology/test_symbolic_dynamics.py` | **`1 failed, 37 passed`**: `assert [1, 0, 1, 0] == [1, 0, 0, 0]` | `git checkout -- tests`, clean |
| 7b (library) | `SymbolicDynamics.py:824`: every child after the first gets `refined_word[:1]` | `pytest -q test_symbolic_dynamics.py test_law_symbolic.py` | **`10 failed`**. The in-splitting law fails on 5 law cases (k10, k28_two_blasts, nested, nested_inner_blasts, heteroclinic): "sibling rows of a differ: a_1 [1, 0, 1, 0] vs a_2 [1, 0, 0, 0]". The old `refined_children_inherit_word` also fails, and so do 5 unit tests | `git checkout -- src`, clean |
| 7c (column-sum identity, the dropped-token case) | `SymbolicDynamics.py:1033` drops every edge into a `*_2` column, so sibling rows stay equal | `pytest -q test_law_symbolic.py` | **`7 failed, 4 passed`**. `refined_matrix_is_an_in_splitting` fails on k10, k28_two_blasts, nested, nested_inner_blasts, inversion_one_blast, inversion_two_blasts and heteroclinic, e.g. "row c: the columns of a's children sum to 1, M[c][a] = 2" | `git checkout -- src`, clean |
| 7d (node-set identity) | `transition_matrix(refined=True)` drops its last node | `pytest -q test_law_symbolic.py` | **`10 failed, 1 passed`**. The in-splitting law fails on every law case except k28_one_blast, whose `M_1` has one node so the mutation is a no-op there, e.g. `(['a_1'], ['a_1', 'a_2'])` | `git checkout -- src`, clean |
| 8 | `PartitionFamily.py:1044` `logger.info(` → `logger.warning(` | `pytest -q tests/unit/topology/test_partition_family.py` | **`1 failed, 8 passed`**: "0 INFO record(s) from 'tanglepack.topology.PartitionFamily', expected 1" | `git checkout -- src`, clean |

**Rubric 8 grep:**
```
grep -n "an odd" src/tanglepack/topology/PartitionFamily.py
678:            incomplete_footprints: The hole images with an odd crossing count      (docstring)
1045:                "an odd number of fold crossings, last one unpaired (incomplete "   (message of logger.info at :1044)
```
`grep -n "logger.warning" -A3 ... | grep -i odd` returns nothing. The only other `logger.warning` left in `_empty_stretches` is the pre-existing "not consecutive" one.

**Beyond the rubric.** The strict xfail on k32 is itself a teeth check on T: if T's helper became vacuous, k32 would XPASS, and a strict xfail turns an XPASS into a failure.

## Rubric 9 — `SymbolicDynamics.py` is docstrings only

- `git diff 3fb310a.. -- src/tanglepack/topology/SymbolicDynamics.py` has four hunks, at `@@ -14`, `@@ -46`, `@@ -756` and `@@ -996`, with 32 changed lines.
- **The AST check.** `astcmp.py` parses the base and branch files, replaces every leading docstring with `pass`, and compares `ast.dump`. It printed **`AST equal after stripping docstrings: True`**.
- **The check is sensitive.** Applied to a copy mutated `index += 1` → `index += 2` (line 476), the same script prints `False`.
- Comments are not in the AST. All four hunks were read by eye: they are module docstring, Dev Notes and class or method docstring text.

## The coder's four deviations

1. **`_empty_stretches` keeps its signature; a sibling `_incomplete_footprints` collects the footprints.**
   - **The difference.** `_empty_stretches` admits a hole bridge when `frozenset(hb)` is in a dict that grows (holes, then each lobe's base and chords as they are found). `_incomplete_footprints` tests against the holes' stretches only.
   - **Measurement.** `fp.py` and `series_probe.py` re-implement the growing filter and compare the admitted hole-bridge lists. They are **identical on all 15 builds and all 47 series steps** (`filter_equal` True everywhere). An odd fold count occurs **nowhere** under either filter: `odd` 0 and `incomplete` 0 on all 62.
   - **Conclusion.** The two cannot diverge on the corpus. They could diverge in principle, but only in what the INFO reports, never in a cut.
   - **The reason given is true.** `test_iterated_cut.py` and `test_iterated_partition.py` monkeypatch `_empty_stretches` with the dict-returning shape.
   - **Remaining cost:** the image chains are looked up twice per build.
2. **T is in `LAW_CASES_ONLY`, and this is forced.**
   - `test_open_p3_deep_runs.py::test_every_law_holds_on_the_deep_runs` calls `assert_non_vacuous` on every layered law not in `LAW_CASES_ONLY`, with no per-run set-aside mechanism.
   - `assert_non_vacuous` fails on a count of 0, and I measured T = 0 on p3_6_blasts myself.
   - T's depth is covered by the series test.
3. **The id `series_members_share_the_P1_itinerary` is the convention.**
   - The three existing series ids are `series_class_count_never_decreases`, `series_class_count_constant_while_holes_constant` and `series_alphabet_inside_minimal_trellis`. All are in `_STANDALONE_ISSUES` and are read by `_params(...)` in the same test file.
   - Under the layered id, `law_tier.set_aside_pairs` would pick up `("members_share_the_P1_itinerary", "k32_series")`, because the pair's law id is in the symbolic layer. The resulting set-aside test would ask the `law_case` fixture for a case that does not exist: a sham xfail.
   - The author's intent is kept: one honest strict xfail of T on k32_series, naming B7b.
4. **The eighth file** is two lines, and both are needed (see rubric 3).

## Gaps and what I could not verify

- **No new tests.**
  - `from_minimal`'s passing of `incomplete_footprints` is exercised on every real build, but its value is always `[]`. Asserting that on a case would pin a case-specific count, which the test policy forbids.
  - The odd case is covered by the stub unit test. A stub is justified because no real case has an odd fold count.
- **The law's `except ValueError` skip never fires** on the 62 builds and steps. It is a documented incomplete-evidence path (`trellis_itinerary` raises when a crossing has no owner). The reviewer may judge it defensive.
- **The `R_7 <-> R_9` wording** is cosmetic (rubric 2). I did not resolve the elements through `ElementNaming`; the inference is from the naming convention.
- **Not re-run by me:** p3_15, p3_16_12 and inversion_period2_10 as standalone T measurements. The deep-run test in the full suite passed on the branch with T excluded, and the in-splitting law included.
- **Commit attribution:** the coder's commits carry `Co-Authored-By: Claude Opus 5.5`. I made no commit, so I added no trailer.
