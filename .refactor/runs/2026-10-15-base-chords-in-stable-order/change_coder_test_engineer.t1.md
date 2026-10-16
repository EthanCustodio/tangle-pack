# Change: B7b `base-chords-in-stable-order` (coder, iteration 1)

**To the test engineer:** you read this report and write the test evidence. `code_reviewer` then grades your work "pass" or "revise" against the plan's rubric (`artifacts/plan_planner_human_2.md` §7). Below I give my own reading of each rubric item. Where my numbers or behaviour differ from the plan, I say so (rubric 5 and rubric 7). Re-measure everything; do not take my word for it.

## Where

| | |
|---|---|
| Base commit | `a6a4f683a3d214d712243110a31e8e2b1a787a49` (main) |
| Branch | `refactor/base-chords-in-stable-order` @ `63477e3` |
| Worktree | `/home/dezhu/Development/TanglePack/tangle-pack-base-chords-in-stable-order` |
| Base worktree | `/home/dezhu/Development/TanglePack/tangle-pack-base-chords-in-stable-order-base`, detached at `a6a4f68` and clean (I copied files in temporarily and restored them with `git checkout -- .`) |
| How I ran things | from inside each worktree: `MPLBACKEND=Agg PYTHONPATH=src ../tangle-pack/env/bin/python -m pytest -q -rxX -p no:cacheprovider` (the main repo's env, with each worktree's `src` first on the path) |

I did not touch the main checkout's uncommitted `docs/`, `.refactor/`, `CLAUDE.md` or `artifacts/`. The one exception is this report, written to `artifacts/` as instructed.

## Commits (`git log --oneline a6a4f68..`)

```
63477e3 Tests: base-branch chords inside the lobe image on the builds, one statement for two branch kinds (B7b, step 6)
c2b570b Tests: the base crossings of the k32 lobe pair in stable order, the toothless chain-order half goes (B7b, step 7)
11ead2f Tests: base-branch chords read by point-in-polygon after every blast; the k32 Theorem T known issue goes (B7b, steps 4 and 5)
69ee644 SymbolicDynamics: T runs on the law cases and every blast series, no known issue left (B7b)
9d9f368 Iterated cut: base-branch chords pair in stable order (B7b, author 2026-10-15)
```

## `git diff --stat a6a4f68..`

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

These are exactly the seven files of plan §3.

## What changed, step by step

1. **The rule** (`PartitionFamily.py`, `_empty_stretches`). The `if branch_key != base_branch:` test is gone, so the loop is now `for on_branch in folds.values(): on_branch.sort(key=stable_cdist)`, then pairing. That made `(base_branch,) = ends` dead, so I deleted it; `assert len(ends) == 1` stays. The comment above the loop is now one sentence: Jordan, the base branch no exception, C19/B7b, 2026-10-15.

2. **PartitionFamily docstrings.**
   - The module bullets state one rule: stable order on every branch.
   - The odd-count clause says "outermost", in the module docstring and in the `_empty_stretches` Note.
   - The k=2.8 example is restated in current ids. It matches the branch (rubric 12).
   - History is one condensed chain: `2a2e35a` → footprint → chain-order base chords → C19 → B7b. It cites the k32 `(6,7)` lobe and `(9,10)` lying outside `f(H)`.
   - "Foreign chords" is renamed "Chords". The known-issue and "slice B7b" sentences are gone, and no class name is mentioned.
   - **One edit beyond the plan's list, in the same file:** in the `_incomplete_footprints` docstring I changed "last" to "outermost", for consistency. That is docstring only. Its log text and `describe()` text still say "last" (code, left untouched per rubric 10).

3. **SymbolicDynamics.py** (docstring only). It now says "laws on the law cases and every blast series", followed by one line: "the base-branch chords pair in stable order since B7b (`PartitionFamily`)".

4. **`tests/helpers/laws.py`.**
   - `_lobe_image` moved here verbatim from `test_foreign_chords.py`.
   - New `_arclength_midpoint`.
   - New `base_chord_misreadings(session) -> (checked, texts)`, placed right after `p1_itinerary_disagreements`. What it does:
     - It takes every hole bridge with a stretch (`_hole_stretches`) whose consecutive `image_chain` has two or more interior crossings on its base branch.
     - The gaps run anchor root → the chain's crossings on that branch (anchors included), in stable order → the branch tail. Each gap is closed at the crossings' own coordinates, read at its arc-length midpoint, and tested with `point_in_polygon(tol=0.0)`. A gap of zero length before an anchor that is itself a chain crossing is skipped.
     - Clauses (a) to (d) are as in plan step 4. Clause (b): consecutive readings differ, and the first gap is inside iff it is the base. Clause (d) calls the kernel over `SimpleNamespace(trellis=..., hole_bridge_ids=[hb])` and checks every `"base"`/`"chord"` pair whose non-anchor ends lie on the base branch.
     - Each misreading text carries the hole bridge, the branch, the clause and the full gap readings.
   - New imports: `SimpleNamespace`, `arc_polyline`, `IteratedHomotopyPartition`, `_hole_stretches`.

5. **`tests/cases.py`.**
   - `SeriesStep` gains `chords_checked` and `chord_misreadings`, documented and filled in `blast_series` next to `p1_*`.
   - The `KNOWN_ISSUES[("series_members_share_the_P1_itinerary", "k32_series")]` entry is deleted.
   - `test_law_case_sanity.py` is untouched.

6. **`test_class_count_under_blasting.py`.** New `test_base_chords_are_inside_the_lobe_image_under_blasting`, parametrised by `_params("series_base_chords_inside_the_lobe_image")`. It prints `chords_checked` per step. The module docstring gains a bullet for it.

7. **`test_foreign_chords.py`.**
   - The module docstring now reads "One statement, two branch kinds".
   - `_build = cache(lambda name: cases.BUILDERS[name]())` is shared by the `case` fixture (now function-scoped, via `_build`) and the new test, so each build is made once per run.
   - New `test_base_chords_are_inside_the_lobe_image` over `k10`, `inversion`, `inversion_period2`, `heteroclinic` and `p3_6_blasts`. It asserts `checked and misreadings == ()`.
   - `p3_6_blasts` is in `NO_BASE_LOBE`. It asserts `not checked` and then takes the visible skip "no stretched hole crosses its base branch twice", so the skip fails loudly if that build ever gains a base lobe.
   - The foreign test itself is unchanged; it imports `_lobe_image` from the helper. `numpy` is no longer imported.

8. **`test_iterated_cut.py`.**
   - New `test_base_crossings_pair_in_stable_order`, with the k32 shape on `OWN` (cdists 8:0.176, 9:6.786, 10:6.665, 11:6.571, 12:7.418, 13:0.036). It asserts `(13,"base") in empty[8]`, `(10,"chord") in empty[11]`, `(12,"chord") in empty[9]` and the two `not in` assertions from the plan.
   - The foreign test lost its base half: the sentence about keeping the chain order, and `(12,"chord") in empty[11]`.
   - The module docstring's first bullet says stable order on every branch, base included.

## Behaviour changes (measured, base vs branch, 63 items)

I adapted the analyst's probe as `scratchpad/b7b/probe2.py` + `run2.sh` + `compare2.py`. The scratchpad is `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b7b/`; raw output is in `out2/`.
- **How it differs from `probe.py`:** the repo comes from `$REPO`, it runs unpatched, it runs once in each worktree, and it records `incomplete_footprints`.
- **Coverage:** 16 builds (the 11 law cases, `p3_16`, `p3_6_blasts`, `k28_8`, `nested_inner6`, `p3_15`) and 47 series steps. There were no errors, and every import resolved to its own worktree (checked via `cases.__file__`).

**Identical on all 63:**
- letter-free `rules`, `M` and `M_1`;
- class, active and inert counts;
- split-active and split-inert counts;
- `is_reliable`;
- `incomplete_footprints` (empty everywhere);
- the in-splitting verdict.

**Different, only on `k32_series`:**

| step | iterated elements | other |
|---|---|---|
| 2-5 | 23 → 24 | none |
| 6-7 | 40 → 39 | `unmatched_members` 1 → 0; Theorem T 1 → 0 disagreements |
| 8-10 | 113 → 111 | `refined_rules` 55 → 54 keys: inert `L#2<->L#4` loses child `4`, and in the word of `R#23<->R#27` that class's token goes from refined index 4 to 3 (`v_4` → `v_3`) |

Nothing differs on any other build or step. This is exactly plan §5.

Other behaviour changes: none in library code beyond the iterated cut. The k32 strict xfail is gone and that test now passes.

## Suite

- **Branch:** `1211 passed, 47 skipped, 5 deselected, 1 xfailed in 220.66s`. The only xfail is `test_orientation_reversing_case_builds[orientation_reversing]` (C3). No XPASS.
- **Base** (`a6a4f68`, unchanged): `1200 passed, 46 skipped, 5 deselected, 2 xfailed in 219.97s`. The xfails are C3 plus the k32 T case (B7b).
- **The difference accounts exactly:**
  - +10 new tests passed: 5 series, 4 builds and 1 unit;
  - +1 former xfail now passed (k32 T);
  - +1 skip (`p3_6_blasts`);
  - 1200 + 11 = 1211.

## Rubric, my own reading

1. **Met.** The branch suite passes, and `-rxX` lists only the C3 xfail. No XPASS.
2. **Met.** The base suite is 1200 passed with the two expected xfails, so every test the branch leaves unchanged passes on base. You should re-run it.
3. **Met.** The diff names the seven files and nothing under `docs/`, `.refactor/`, `CLAUDE.md` or `artifacts/`.
4. **Met.** See the behaviour section: the identical set holds on all 63, and the differences are exactly the §5 list.
5. **Grep met; the line count needs your judgement.** The grep is empty. The code diff of `_empty_stretches` (excluding comments) is 4 removed lines and 2 added:
   - removed: `(base_branch,) = ends`, `for branch_key, on_branch in folds.items():`, the `if`, the old indented sort;
   - added: the new `for`, the dedented sort.

   That is 4 changed lines if you count removals, or 6 if you count `+` and `-` together. The extra line is the dead `(base_branch,) = ends` I deleted (the plan: "rename only if it becomes unused"; it did become unused).
6. **Met.** `test_members_share_the_P1_itinerary_under_blasting[k32_series]` passes on the branch with 0 disagreements. `p1_tested` = `[0, 0, 0, 2, 3, 4, 9, 13, 13, 13, 13]` on both branch and base. Every other series passes.
7. **Met, with one deviation from the plan's prediction.**
   - Branch `chords_checked` per step (all five pass):

     | series | per step | sum |
     |---|---|---|
     | `k28_series` | `[0,0,1,1,1,1,1,1,1]` | 7 |
     | `p3_series` | all 0 | 0 |
     | `inversion_series` | `[1,1,1,1,1]` | 5 |
     | `k32_series` | `[0,0,1,1,1,1,2,2,2,2,2]` | 14; at least 1 on every step 2-10 |
     | `nested_outer_series` | `[0,0,1,1,1,1,1,1,1,1,1,1]` | 10 |

     The series total is 36. Adding the 4 build lobes from item 8 gives **40 readings by the tests**, against the analyst's 47. The analyst also read builds the test does not parametrise, such as `k28_two_blasts`, `nested` and `nested_inner_blasts`.
   - **Base** (helper, `cases.py` and the series test copied in, then reverted): FAILS on `k32_series`, with misreadings at every step 2-10.
     - Clause (c) "the inside gap (9,12) is not empty" fires at steps 2-7.
     - Clause (d) fires at **2-10**, not only 8-10 as the plan predicted. Chain chord `(11,12)` spans three gaps, so it is never one inside gap, from step 2 on. At 6-10 the chord `(9,10)` is also flagged.
     - The readings are identical at every step: `(anchor,13)out (13,8)B (8,11)out (11,10)in (10,9)out (9,12)in (12,tail)out`.
     - This is stricter than predicted, not weaker: (d) is the plan's clause, and it catches both wrong chain chords.
   - Please reproduce both outputs.
8. **Met.** `checked` = 1 on `k10`, `inversion`, `inversion_period2` and `heteroclinic`, with no misreadings. `p3_6_blasts` shows a visible SKIP with checked = 0.
9. **Met.** On base, `test_base_crossings_pair_in_stable_order` fails with `assert ((10, 'chord') in [(12, 'chord')])`: base gives `empty[11] == [(12,'chord')]`, i.e. the chain pairs `(11,12)` and `(9,10)`. It passes on the branch.
10. **Met.** All six `_empty_stretches` call sites pass unchanged, in `test_iterated_cut.py`, `test_iterated_partition.py`, `test_foreign_chords.py` and `test_partition_family.py`. `_incomplete_footprints` has no code change; its docstring changed one word ("last" → "outermost"). Its unit test passes.
11. **Met.**
    - `chain order` in `PartitionFamily.py` matches only line 91, inside the History bullet. I reworded a second mention in that bullet to "pairing along the chain" so the grep finds the one line.
    - It is empty in `test_iterated_cut.py`.
    - In `SymbolicDynamics.py`, `B7b` appears only in "the base-branch chords pair in stable order since B7b". No known-issue claim remains.
    - `R_7 <-> R_9` is empty in `tests/`.
12. **Met.** `session.iterated_partition().describe()` on `k28_two_blasts` at the branch:

    ```
    left : #0 [0,5.993) anchor->4 (p#0) | #1 [4->9] (p#0) | #2 (9->10) (p#1) | #3 [10->3] (p#2)
    right: #0 [anchor->11] (p#0) | #1 (11->8) (p#0) | #2 [8->7] (p#0) | #3 (7->6) (p#1) | #4 [6->5] (p#2)
           | #5 (5->4) (p#3) | #6 [4->9] (p#4) | #7 (9->10) (p#4) | #8 [10->3] (p#4)
    ```

    So the right row reads `R_1^1=[0,11] R_1^2=(11,8) R_1^3=[8,7]` and `R_5^1=[4,9] R_5^2=(9,10) R_5^3=[10,3]`. The left row reads `[·,4) [4,9] (9,10) [10,3]`, i.e. `L_1^1=[0,4) L_1^2=[4,9]`. This matches the new docstring. The print script is `scratchpad/b7b/k28_example.py`.

## Notes for the reviewer

- `base_chord_misreadings` imports the private `_hole_stretches` from `PartitionFamily`; the helpers already import private kernels (`_side_of`). The helper reads only polylines numerics already computed, and calls the map never.
- The foreign test's `case` fixture moved from module scope to function scope backed by `functools.cache`. That keeps the builds alive for the run, as the plan asks.
- Nothing needed a file outside the seven.
