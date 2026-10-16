# Test evidence: B0 `oracle-laws` (tests only)

*Test engineer, 2026-10-12. Plan: `artifacts/plan_planner_human_2.md`, approved in
`artifacts/plan_ok_planner_human.json` (with the addendum for the deep-run test). Coder report:
`artifacts/change_coder_test_engineer.md`.*

- **Base:** `a3bd8f0`, in the detached worktree `../tangle-pack-oracle-laws-base`, which I created.
- **Branch:** `refactor/oracle-laws` at `1029fcc`, in the worktree `../tangle-pack-oracle-laws`.
- **New tests:** I added none and committed nothing. This slice's tests are the deliverable, and I
  found no gap in the rubric that needed a test.
- **Final state:** both worktrees are clean (`git status --short` is empty). The author's tree was
  not touched apart from this file.

Every run used this command, from the worktree named:

`MPLBACKEND=Agg PYTHONPATH=src /home/dezhu/Development/TanglePack/tangle-pack/env/bin/python -m pytest -q -rxXs -p no:cacheprovider …`

I checked `tanglepack.__file__` first. It resolves to `<worktree>/src/tanglepack/__init__.py` in both
worktrees, so each run imports that worktree's own library.

## 1. Full suite, base against branch

| Tree | Result |
|---|---|
| base `a3bd8f0` | **1117 passed, 27 skipped, 5 deselected, 1 xfailed** in 157 s |
| branch `1029fcc` | **1117 passed, 38 skipped, 5 deselected, 1 xfailed** in 157 s |

**The base numbers differ from the job's expected 1122 / 26.** Those figures date from `69371fd`.
Commit `74a19da` (D15, the base's parent) says so itself: its message reads "Suite 1117 pass, 1 xfail".
It also added `p3` to the `walks_avoid_disjoint_images` NOT_APPLICABLE tuple, which is the 27th skip.
So 1117 / 27 is the correct base.

**Every difference is accounted for.**
- **Collection.** I diffed `--collect-only` ids: the base collects 1145 tests and the branch 1156.
- **The 11 extra ids** are all `test_symbolic_set_aside` parameters:
  - `[every_active_class_walks-inversion_period2]`;
  - `[k10_is_the_full_shift-<case>]` for the 10 law cases other than k10.
- **Skips.** All 11 are skips, so skips go 27 → 38 and passes stay at 1117.
- **Unchanged.** No other id appears or disappears, and the `test_partition_laws[*]` /
  `test_symbolic_laws[*]` ids are identical (rubric 6).
- **Skip lists.** The diff of the `-rs` skip lists is exactly those 11 lines.
- **Xfail.** On both trees the only xfail is
  `test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing]` (C3). No
  `reference_pairs_clear_of_X` KNOWN_ISSUES entry exists or is needed (rubrics 1 and 11).

**Other repository facts.**
- `git diff a3bd8f0 -- src` is empty: 0 lines.
- `git diff --stat a3bd8f0` lists exactly six files, as the addendum allows:

  | File | Lines changed |
  |---|---|
  | `tests/cases.py` | 24 |
  | `tests/conftest.py` | +3 |
  | `tests/helpers/laws.py` | +169 |
  | `tests/invariants/test_law_partition.py` | +3 |
  | `tests/invariants/test_law_symbolic.py` | 3 |
  | `tests/regression/test_open_p3_deep_runs.py` | 15 |
  | **Total** | 201 insertions, 16 deletions |

## 2. The new laws, run on the OLD library

I put the branch's tests on the base worktree with `git -C ../tangle-pack-oracle-laws-base checkout
1029fcc -- tests/`. That means `tests/cases.py` there also builds `inversion_period2` at 9/10 steps,
which is intended. I then ran:

`pytest tests/invariants/test_law_partition.py tests/invariants/test_law_symbolic.py tests/regression/test_open_p3_deep_runs.py`

**Result: 36 passed, 32 skipped.**
- **Law tier.** `test_partition_laws` and `test_symbolic_laws` pass on all 11 law cases.
- **Deep runs.** All 14 deep-run tests pass.
- **Skips.** The 32 skips are the set-aside tests, including the 11 new ones.
- **Verdict.** All four new laws hold on the old code for every law case. They pin the base behaviour.

**Rubric 10: the shared-end closure, toggled on base.** I made a scratch edit in the base worktree,
`laws.py:1690` `if _table_image(...)` → `if False and _table_image(...)`, and ran
`test_law_partition.py`. Result: **2 failed, 9 passed**.

| Case | Hole (origin) | End | Iterate | Container | Carried cdist | Container's `hi` |
|---|---|---|---|---|---|---|
| `test_partition_laws[k28_one_blast]` | (4, 5) | 5 | -1 | (0, 3) | 42.250102 | 42.250098 |
| `test_partition_laws[inversion]` | (2, 3) | 3 | -2 | (0, 1) | 0.88968 | 0.88903 |

I also ran a scratch driver, `check_propagated_hole_inside_container(build_inversion_period2(stable_steps=N))`, on base:

| Build | Closure off | Closure on |
|---|---|---|
| 9/9 | RAISES: hole (6, 7), end 7, iterate -1, container (0, 3), cdist 0.712531 > hi 0.712482 | passes, count 16 |
| 9/10 | passes, count 16 | passes, count 16 |

This matches the plan's [changed] note and the coder's rubric-10 claim.

**Restoring the base tests.** I used `git checkout HEAD -- tests/`; the base worktree is clean again.
The job's suggested `git checkout -- tests/` would NOT restore it. The earlier `checkout 1029fcc --
tests/` staged the branch's files, so a plain `checkout -- tests/` restores from that index, which is
the branch version.

**Rubric 4 (no test pins these numbers).** On the branch, `build_inversion_period2` gives:

| Build | Classes | Holes | `is_reliable` |
|---|---|---|---|
| 9/10 (new default) | 6 | 10 | True |
| 9/9 | 6 | 10 | False |

**Rubric 13.** `tests/cases.py:826` reads `stable_steps: int = 10`. The function's docstring no
longer says "past the pip orbit". The phrase does remain in the module docstring at
`tests/cases.py:52`, where it describes the separate 10/10 build, and that is still accurate.

## 3. The new laws and what each guards

All four are in `tests/helpers/laws.py` and are registered in `LAYERS` (rubric 6): three under
`"partition"` (`:3105-3107`) and one under `"symbolic"` (`:3165`).

- **`check_chosen_pip_is_strong`** (`:1566`).
  - **What it computes.** It re-derives A.4's box test from the registry and `FixedPoint`
    (`branch_cycle`, `per_step_beta`) alone. Every homoclinic crossing of the pip's fixed point is
    carried back onto the pip's stable branch by scaling.
  - **What it asserts.** None lies in the open box `(0, c_s(q0)) x (0, c_u(q0))` (`:1587`).
  - **What it guards.** `StrongPip.is_strong_pip` and the pip choice.
- **`check_reference_pairs_clear_of_X`** (`:1632`).
  - **What it computes.** It tests both open intervals of every reference pair, using the cycle and β
    of the interval's own fixed point (D7). Every non-anchor crossing of the cut N is scaled to every
    lap landing inside the interval.
  - **What it allows.** A landing passes only if the iterate table links it to an end of the pair, or
    if it collides with an end on both cdists within `SCALING_RTOL`. Anything else raises (`:1665`).
  - **What it guards.** `Pseudoneighbor._unstable_interval_is_clear` and the pair definition.
  - **Where it runs.** On the law cases only: it is listed in `LAW_CASES_ONLY` (`:3172`).
- **`check_propagated_hole_inside_container`** (`:1672`).
  - **What it computes.** It takes every hole with `iterate < 0` and carries the reference bridge's
    ends back with the REAL `map_inv_batch`. Each carried end is projected onto the container's node
    polyline.
  - **What it asserts.** The projection lies strictly inside the container's crossing cdists
    (`:1701`). It is closed at an end that the table links to a container crossing (`:1690`).
  - **What it guards.** The backward hole chain, `StablePartition._containing_bridge` and the bisect.
- **`check_k10_is_the_full_shift`** (`:3024`).
  - **What it computes.** It returns 0 unless the case is k10 (the set-aside contract). On k10 it is
    letter-free.
  - **What it asserts.**
    - There is one active class.
    - Its active tokens are that class twice, with directions `[-1, 1]`.
    - The two tokens have two distinct occurrences.
    - The class matrix is `[[2]]`.
    - The refined matrix is `[[1, 1], [1, 1]]`.
  - **What it guards.** The words, `refine`, and both transition matrices: the Q11 textbook pin.
    The pin is recorded in `tests/conftest.py:78` and `test_law_symbolic.py:15`.

**The deep-run test.** `tests/regression/test_open_p3_deep_runs.py:92-98` now loops
`assert_non_vacuous` over every check except `LAW_CASES_ONLY`, where it used to call `run_layer`.
- **What changed.** The loop stops at the first failure instead of collecting them.
- **What did not weaken.** The deep builds have no NOT_APPLICABLE or KNOWN_ISSUES pairs, so dropping
  `run_layer`'s skip loses nothing. Every other law is still required to be non-vacuous.
- **Result.** All 14 deep-run tests pass on the branch and on base with the branch's tests.

## 4. Break-it checks (on the branch worktree's `src/`)

Each break was made with `sed` and then reverted with `git -C ../tangle-pack-oracle-laws checkout --
src/`. After each one, `git status --short` was empty (0 lines).

**Check 1: `chosen_pip_is_strong`. Red when a rejected pip is chosen.**

- **Edit (a).** `StrongPip.py:214`: `if best_id is not None:` → `if False and best_id is not None:`.
  With it, `is_strong_pip` accepts every crossing.
- **Driver run.** I ran the scratch script `pip_driver.py` on k10. It builds k10, sets the pip and
  calls `topology()`, then calls the check.

  | Run | Candidates / pip | Check result |
  |---|---|---|
  | unbroken, `list` | candidates `[1, 2, 3, 4, 7]` | count 7 |
  | edit (a), `list` | candidates `[1, …, 7]` | — |
  | edit (a), pip 5 | 5 | **RAISES** at `laws.py:1587`: `(2, 0, s=74.146 < s0=76.536, u=74.146 < u0=651.77)` |
  | edit (a), pip 6 | 6 | **RAISES** on crossing 2: `s0=81.874`, `u0=697.23` |

- **Edit (a) alone is NOT observable in the law tier.** `test_law_partition.py` gives 11 passed. The
  default policy is `min(unstable_cdist)` (`Trellis.py:1106`), and it still picks a genuinely strong
  candidate. The break only shows once a rejected crossing is chosen.
- **Edit (a) plus (b), source only.** Edit (b) is `Trellis.py:1106` `chosen = min(` → `chosen = max(`.
  - In the law tier, `chosen_pip_is_strong` goes red on `test_partition_laws[inversion_period2]`:
    `AssertionError: (4, 3, 0.938 < 1.075, 0.150 < 1.090)`.
  - The run is noisy: 6 failed and 7 setup errors. An outermost pip empties the reference window, so
    many other laws become vacuous. On k10 the outermost pip (7) is genuinely strong.
  - The driver run above is the clean evidence.

**Check 2: `reference_pairs_clear_of_X`. Red on 7 of 11 law cases.**

- **Edit.** `Pseudoneighbor.py:217`: `if _unstable_interval_is_clear(…):` → `if True or
  _unstable_interval_is_clear(…):`. Every consecutive pair in the window is then accepted.
- **Command.** `pytest tests/invariants/test_law_partition.py`. Result: **7 failed, 4 passed**.
- **Where `reference_pairs_clear_of_X` fails** (`laws.py:1665`):

  | Case | First failing pair and landing |
  |---|---|
  | k10 | pair (1, 6), `M^1(1)` |
  | k28_two_blasts | pair (3, 10), `M^1(3)` |
  | nested | pair (4, 38), `M^1(4)` |
  | nested_inner_blasts | pair (4, 42), `M^1(4)` |
  | inversion_one_blast | pair (1, 27), `M^2(1)` |
  | inversion_two_blasts | pair (1, 27), `M^2(1)` |
  | heteroclinic | pair (4, 114), `M^1(4)` |

  Each of these lies inside its open unstable interval.
- **Other laws.** On the same cases `reference_pairs_valid` and
  `propagated_holes_land_on_predicted_branch` also go red. So this break is not caught by the new law
  alone.
- **Not observable on `k28_one_blast`, `p3`, `inversion`, `inversion_period2`.** There every
  consecutive pair in the window is already clear, so the break changes no pair. For example, the
  driver gives p3's pairs as `[(3, 4), (4, 5)]` with count 4 either way. That is a fact about those
  cases, not a weakness of the law.

**Check 3: `propagated_hole_inside_container`. Red on 9 of 11 law cases.**

- **Edits.**
  - `StablePartition.py:1535`: `return trellis.bridge_between(ids[i], ids[i + 1])` → `return
    trellis.bridge_between(ids[i + 1], ids[i + 2]) if i + 2 < len(ids) else None`.
  - `:1542`: `if coords is None or …` → `if True or coords is None or …`. This makes `_map_backward`
    return None, as in plan rubric 9.
- **Command.** `pytest tests/invariants/test_law_partition.py`. Result: **11 failed**. The law fails
  (`laws.py:1701`) on 9 cases:

  | Case | Hole (origin) | End | Container | Carried cdist | Container's `lo` |
  |---|---|---|---|---|---|
  | k10 | (5, 6) | 5 | (4, 3) | 83.938 | 84.369 |
  | k28_one_blast | (4, 5) | 4 | (4, 3) | 42.1996 | 42.2501 |
  | p3 | (3, 4) | — | (22, 21) | — | — |
  | nested | (7, 8) | — | (26, 25) | — | — |
  | nested_inner_blasts | (7, 8) | — | (26, 25) | — | — |
  | inversion | (1, 2) | — | (14, 13) | 5.909 | 1.745 |
  | inversion_two_blasts | (26, 27) | — | (21, 23) | — | — |
  | inversion_period2 | (13, 14) | — | (6, 5) | — | — |
  | heteroclinic | (92, 93) | — | (55, 56) | — | — |

- **The new law is the only real red on k10, k28_one_blast, inversion and inversion_period2.**
  - Every other entry for those cases is "checked nothing" for `anchor_opening_faces_its_hole`.
  - `backward_holes_keep_bridge_side` catches the break only on p3, nested, nested_inner_blasts,
    inversion_two_blasts and heteroclinic.
  - So check 3 adds coverage that did not exist.
- **The other two cases.** On k28_two_blasts and inversion_one_blast the shifted bisect leaves no
  propagated hole. The law is then vacuous there, and `assert_non_vacuous` still turns it red.
- **Control.** The `_map_backward` edit alone gives 11 passed, so the red comes from the container
  shift.

**Check 4: `k10_is_the_full_shift`. Red on k10, two breaks.**

- **Break (a).** `SymbolicDynamics.py:590` (`refine`): `if len(seen) < 2:` → `if len(seen) < 3:`. The
  k=10 class is then not split.
  - Command: `pytest tests/invariants/test_law_symbolic.py`. Result: 6 failed, 5 passed.
  - `test_symbolic_laws[k10]`: `k10_is_the_full_shift: AssertionError: assert [[2]] == [[1, 1], [1, 1]]`.
  - The other five failures only make `refined_children_inherit_word` and `member_matching_consistent`
    vacuous.
- **Break (b).** `SymbolicDynamics.py:973` (`transition_matrix`): `= weight` → `= 1`.
  - Result: **1 failed, 10 passed**.
  - `test_symbolic_laws[k10]`: `k10_is_the_full_shift: AssertionError: assert [[1]] == [[2]]`.
  - **It is the only law to catch this break.** `check_matrix_is_token_counts` compares the REFINED
    matrix only, and every refined entry on the law cases is 1. So the class matrix (`refined=False`)
    was unguarded before this slice.

## 5. Independence (rubric 5)

I extracted the bodies of the four checks and their five helpers (`laws.py:1566-1702`, `:3024-3041`)
and grepped them. Each of the following has **0 hits**: `StrongPip`, `is_strong_pip`,
`Pseudoneighbor`, `_unstable_interval_is_clear`, `StablePartition`, `_containing_bridge`,
`bridge_for_pair`, `image_cdist`, `import`.

- **Check 1** (`:1566-1588`) has 0 hits for `iterate`.
- **Checks 2 and 3** use the word only as `registry.iterate_table` (numerics data, read in
  `_table_image`), `hole.iterate` and docstring prose.

**Module-level imports of `laws.py`.** These predate the slice and the new code does not use them:
- `:46-51`: `inspect`, `defaultdict`, `typing`, `numpy`, `pytest`;
- `:53` `helpers.invariants`;
- `:59-78`: `geometry`, `Arrangement`, `BridgeClass`, `DualGraph.face_side`, `DualWalk.walk_chords`,
  `MinimalTrellis.image_chain`, `Pseudoneighbor.compute_pseudoneighbors`, `StablePartition`,
  `Trellis.SCALING_RTOL`.

**What the new checks read from the library.** Only the products under test (`trellis.pips`,
`reference_pseudoneighbors`, `holes`, `bridge_between`, `own_intersection_ids`), the registry, and
`FixedPoint`'s `branch_cycle`, `per_step_beta` and `advance_key`. Check 3 also calls the real
`map_inv_batch`.

## 6. What I could not test, and nits

- **Check 2 on four cases.** Its break cannot be observed on `k28_one_blast`, `p3`, `inversion` and
  `inversion_period2`: no consecutive pair there is punctured, so forcing acceptance changes nothing
  (§4).
- **Check 1's break needs a driver.** The `is_strong_pip` break shows only through a scratch driver
  that chooses a rejected pip, or noisily with the default policy also flipped (§4). No law case's
  default choice is affected by `is_strong_pip` accepting everything.
- **No break tests are committed.** This follows plan §6: they patch private library code. They were
  scratch scripts in my session scratchpad (`pip_driver.py`, `pairs_driver.py`, `ip2_9.py`).
- **Style nit** in `tests/helpers/laws.py`. There are three blank lines before
  `def check_k10_is_the_full_shift` (`:3021-3023`) and one before the `#: The read-only law checks`
  comment (`:3042`). PEP 8 wants two in both places. It is cosmetic and I left it unchanged.
