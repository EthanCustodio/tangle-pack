# Change: B3 `depth-of-n` (coder, iteration 1)

*2026-10-14. Plan `artifacts/plan_planner_human_2.md`, approved as written: Q1 = (a), Q2 = "non-decreasing".*

## Where

| | |
|---|---|
| Base commit | `9cd29d5fb146f10cced2e95f134cb141b9c8e04a` (the author's `HEAD`, `main`) |
| Branch | `refactor/depth-of-n` |
| Worktree | `/home/dezhu/Development/TanglePack/tangle-pack-depth-of-n` |
| Base worktree | `/home/dezhu/Development/TanglePack/tangle-pack-depth-of-n-base` (detached at `9cd29d5`) |

I did not touch the author's working tree (its dirty `docs/`, `.refactor/` and `artifacts/` files). The only file I wrote there is this artifact.

**Commits (`git log --oneline 9cd29d5..`):**
```
33f6e60 Tests: SeriesStep.pip_depth composed test-side and the series law 'the pip orbit never gets shallower' (B3, Q2)
aa70900 Tests: the proximity guard exempts a child ending on the pip orbit; pip_depth is the composed orbit (B3)
c2e3bd2 Blast: a child ending on the chosen pip's orbit is kept regardless of min_separation; BlastResult.pip_depth (B3, N-a)
```

**`git diff --stat 9cd29d5`:**
```
 src/tanglepack/loom/Blast.py                       | 132 +++++++++++++++------
 tests/cases.py                                     |  15 +++
 .../regression/test_class_count_under_blasting.py  |  13 +-
 tests/unit/loom/test_blast_pip_orbit.py            |  61 ++++++++++
 tests/unit/loom/test_blast_proximity_guard.py      |  43 +++++--
 5 files changed, 214 insertions(+), 50 deletions(-)
```
- These are exactly the plan's five files, with 264 changed lines (the cap is 300).
- Ignoring whitespace (`git diff -w --shortstat`), the diff is 181 insertions and 17 deletions. Most of `Blast.py`'s count is the child loop moving out one indentation level (see step 1).

## What was done, step by step

1. **Pip-orbit exemption (`Blast.py`).**
   - A new module helper, `_pip_orbits(session, fixed_points)`, returns `[q0, f(q0), ..., f^m(q0)]` per selected fixed point with a chosen pip. It composes the chain link by link with `registry.iterate_table[x, 1]`. It never calls `forward_depth`.
   - The pip's current id comes from `registry.find(session.chosen_strong_pip(fp))`, the same lookup `TangleSession._restore_strong_pip` uses. This differs from the plan's snippet, which read `.id` off the remembered `Intersection`; that id can be stale after a recompute without `preserve_ids`.
   - I left out the plan snippet's `x not in orbit` cycle guard: a pip is not periodic, so its chain cannot cycle.
   - The per-iteration loop is split in two. The first loop maps every frontier bridge (`iterate_bridge`, which links `f^m -> f^{m+1}`). Then the orbit is read once, as `{id: (fp, n)}`. Then the second loop sorts the children (`for child in produced:`), in the same order and with the same tests as before.
   - Order in the child loop: seen check, exterior check, then the exemption (`pip_end = next((x for x in child.id or () if x in orbit), None)`), then the separation test (`elif min_separation is not None:`).
   - One INFO line is logged per exempted child, only when `min_separation` is set: `blast: child %s ends at f^%d(q0) of %s (crossing %d), kept regardless of separation`.
   - Partial children have `id is None`, so they are never pip-bearing.
2. **`BlastResult.pip_depth: dict[FixedPoint, int]`.** It is filled after the final `infer_iterates()`, as `len(orbit) - 1` per fixed point with a chosen pip. The docstring now documents `too_close`, `already_known` and `pip_depth`.
3. **`min_separation` docstring.** It now says:
   - it is a cost guard that tests curve already in the tangle as well as kept siblings (fixing the old "sibling" wording);
   - a child ending on a selected fixed point's chosen-pip orbit is exempt (N-a, Q2) and logged at INFO;
   - the guard never decides the depth of `T(q0)`; that is decided by whether `f^{m+1}(q0)` is detected and resolved (A.4).
4. **Series data and law.**
   - `cases.SeriesStep.pip_depth: int`, filled by `cases._pip_depth(session, fp)`. That helper calls `registry.find(session.chosen_strong_pip(fp))`, then steps `trellis.iterate(pip, m + 1)` until it returns None. It never reads the library's count.
   - New law `test_pip_orbit_never_shallower(name)` in `tests/regression/test_class_count_under_blasting.py`: `depths == sorted(depths)`. It pins no number.
   - The module docstring gets one bullet for the law.
5. **Proximity-guard test (`test_min_separation_drops_close_bridges`).**
   - It is now parametrized over `pip_chosen` (`pip` / `no_pip`); in `no_pip` the pip is forgotten (`session.forget_strong_pip(fp1)`) before the blast.
   - Both variants check that every kept child keeps its distance unless it ends on the chosen pip's orbit (composed test-side after the blast), and that every dropped child was too close.
   - `too_close > 0` is asserted in the `no_pip` variant only. With the pip chosen, the fixture's one too-close child, which ends at `f^4(q0)`, is exempt, so the guard drops nothing.
6. **New `tests/unit/loom/test_blast_pip_orbit.py`.** All three tests use `henon_p3_session`, `fp1`, 3 iterations and separation 0.02.
   - (a) `test_pip_bearing_child_is_kept`: at least one kept child ends on the composed orbit, and no too-close child does.
   - (b) `test_pip_bearing_child_is_dropped_without_a_pip`: the pip is read, then forgotten, then the blast runs. Some too-close child ends on that pip's composed orbit.
   - (c) `test_pip_depth_is_the_registered_orbit`: `result.pip_depth == {fp1: len(orbit) - 1}`.
   - Deviation in (b): the plan said to match "the same child (by bridge-id end pair)" as in (a). (a) and (b) run on separate fixture builds, so (b) finds the child by the same rule instead, as an end on the pip's composed orbit. No id is pinned or carried between builds.

## Suite

**Branch:** `MPLBACKEND=Agg env/bin/python -m pytest -q -rxX`, run in the branch worktree:
```
XFAIL tests/invariants/test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing] - C3: INV-11: det J < 0 is out of scope (C3); FixedPoint.set_k_value rejects a saddle with exactly one negative eigenvalue
1192 passed, 38 skipped, 5 deselected, 1 xfailed in 203.52s (0:03:23)
```
- The one xfail is `cases.KNOWN_ISSUES[("orientation_reversing_case_builds", "orientation_reversing")]`.

**Base, same command:** `1183 passed, 38 skipped, 5 deselected, 1 xfailed in 217.48s`.
- The +9 tests are: 5 series laws, 3 unit tests in the new module, and the second parameter of the proximity test.

**Branch tests against the base source (rubric 2).** I made a scratch copy of the base `src/` with the branch's `tests/` and `pyproject.toml`; pytest's `pythonpath=["src", "tests"]` picks up the base library. Selection: `test_blast_pip_orbit.py`, `test_blast_proximity_guard.py`, `test_pip_orbit_never_shallower`.
```
FAILED tests/unit/loom/test_blast_pip_orbit.py::test_pip_bearing_child_is_kept
FAILED tests/unit/loom/test_blast_pip_orbit.py::test_pip_depth_is_the_registered_orbit   (AttributeError: 'BlastResult' object has no attribute 'pip_depth')
2 failed, 9 passed in 73.47s
```
- 6(a) and 6(c) fail on base.
- These pass on base: 6(b), both `test_min_separation_drops_close_bridges` variants, `test_min_separation_none_is_unchanged_default`, and `test_pip_orbit_never_shallower` on all 5 series.

## Measurements (rubrics 4 to 7), base and branch

**How they were run.**
- Scripts are in `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b3measure/`: `measure.py`, `run.sh`, `two_iter.py`, `p3script.py`, `run2.sh`. Outputs are `{base,branch}_*.txt`, with the full JSON dumps in `*.err`.
- Each was run as `PYTHONPATH=<worktree>/src:<branch>/tests`.
- The branch's `cases.py` is used for both worktrees: its `pip_depth` is test-side composition, so on base it measures the base library. The `.err` headers confirm which `tanglepack/__init__.py` was loaded.

### Rubric 4: law-case fingerprints
- Method: `helpers.law_tier.fingerprint` of each `cases.LAW_CASES` build.
- Dropped from the comparison:
  - the `products` aspect, which holds `id()`s of objects and differs per process;
  - `0x…` addresses in object reprs, which were stripped. My first run did not strip them, and every hash differed for that reason alone.
- Result: the hashes are identical on base and branch (`diff` empty). The full canonical JSON in `*_fingerprints.err` is identical too.
```
k10 c4a754225a4b6a9f 1444
k28_one_blast 6157ebf2fc37ebcd 1585
k28_two_blasts 33a1aef77f71ebe9 1817
p3 ac323fe04d94ea76 4212
nested 40f53c21eaa6b7f8 5901
nested_inner_blasts 99173edb9a6a6c51 6134
inversion 62626ed48e88b716 2853
inversion_one_blast 8b5dedb978178036 3487
inversion_two_blasts 783724971bdbfd65 4086
inversion_period2 7b50b0585b38482a 4627
heteroclinic de606d5056d3539c 27063
```

### Rubrics 5 and 6: every `SeriesStep` field, all five series
- Every field of every step was compared: `blasts, holes, classes, closed, ambiguous, active, inert, outside_minimal, references, unwitnessed, letters, pip_depth, history_free, crossings, recomputed`.
- They were dumped as canonical JSON, one line per step. Base and branch are **byte-identical** on all five series (`diff` empty on the summary lines and on the full dumps).
- Summary per step (`sha` = hash of the full step):
```
== k28_series (branch; base identical)
    {"blasts": 0, "pip_depth": 1, "holes": 4, "classes": 3, "closed": 3, "references": 2, "letters": 3, "N": 4, "bridges": 3, "sha": "b78bd2d7883fc75c"}
    {"blasts": 1, "pip_depth": 2, "holes": 4, "classes": 3, "closed": 3, "references": 2, "letters": 1, "N": 6, "bridges": 4, "sha": "f6f18ca7c8a5797c"}
    {"blasts": 2, "pip_depth": 3, "holes": 4, "classes": 4, "closed": 4, "references": 1, "letters": 3, "N": 10, "bridges": 7, "sha": "8a8e93dc3c3a7ee0"}
    {"blasts": 3, "pip_depth": 4, "holes": 4, "classes": 4, "closed": 4, "references": 1, "letters": 3, "N": 18, "bridges": 13, "sha": "d682aea55b8213a3"}
    {"blasts": 4, "pip_depth": 5, "holes": 4, "classes": 4, "closed": 4, "references": 1, "letters": 3, "N": 30, "bridges": 21, "sha": "acad3f766b833ebc"}
    {"blasts": 5, "pip_depth": 6, "holes": 4, "classes": 4, "closed": 4, "references": 1, "letters": 3, "N": 50, "bridges": 35, "sha": "8e7c5081e0f40e86"}
    {"blasts": 6, "pip_depth": 7, "holes": 4, "classes": 4, "closed": 4, "references": 1, "letters": 3, "N": 85, "bridges": 60, "sha": "065f39e38fb117a4"}
    {"blasts": 7, "pip_depth": 8, "holes": 4, "classes": 4, "closed": 4, "references": 1, "letters": 3, "N": 143, "bridges": 101, "sha": "35ef129dcfa9d66a"}
    {"blasts": 8, "pip_depth": 8, "holes": 4, "classes": 4, "closed": 4, "references": 1, "letters": 3, "N": 239, "bridges": 168, "sha": "94dadc78d70be2a4"}
== p3_series (branch; base identical)
    {"blasts": 0, "pip_depth": 14, "holes": 8, "classes": 5, "closed": 5, "references": 2, "letters": 3, "N": 30, "bridges": 27, "sha": "92356dde35d0ba1d"}
    {"blasts": 1, "pip_depth": 15, "holes": 8, "classes": 5, "closed": 5, "references": 2, "letters": 3, "N": 32, "bridges": 28, "sha": "4c9415b5f150eb5b"}
    {"blasts": 2, "pip_depth": 16, "holes": 8, "classes": 5, "closed": 5, "references": 2, "letters": 3, "N": 34, "bridges": 29, "sha": "c011efa2d4e69066"}
    {"blasts": 3, "pip_depth": 17, "holes": 8, "classes": 5, "closed": 5, "references": 2, "letters": 3, "N": 36, "bridges": 30, "sha": "da8eae66c28038f5"}
    {"blasts": 4, "pip_depth": 18, "holes": 8, "classes": 5, "closed": 5, "references": 2, "letters": 3, "N": 38, "bridges": 31, "sha": "51eb2cb4c84841b1"}
    {"blasts": 5, "pip_depth": 19, "holes": 8, "classes": 5, "closed": 5, "references": 2, "letters": 3, "N": 40, "bridges": 32, "sha": "4aa24f9c7436b9b5"}
    {"blasts": 6, "pip_depth": 20, "holes": 42, "classes": 23, "closed": 43, "references": 2, "letters": 22, "N": 46, "bridges": 37, "sha": "c5826d9c5471de67"}
    {"blasts": 7, "pip_depth": 21, "holes": 42, "classes": 23, "closed": 43, "references": 2, "letters": 21, "N": 52, "bridges": 40, "sha": "0b47770367fc33a3"}
    {"blasts": 8, "pip_depth": 22, "holes": 42, "classes": 23, "closed": 43, "references": 2, "letters": 21, "N": 58, "bridges": 43, "sha": "308627d023dffe23"}
    {"blasts": 9, "pip_depth": 23, "holes": 42, "classes": 23, "closed": 43, "references": 2, "letters": 21, "N": 67, "bridges": 49, "sha": "ce7b4cb6a2e66e97"}
== inversion_series (branch; base identical)
    {"blasts": 0, "pip_depth": 5, "holes": 6, "classes": 4, "closed": 4, "references": 2, "letters": 2, "N": 13, "bridges": 12, "sha": "c56ee44c08eb257a"}
    {"blasts": 1, "pip_depth": 6, "holes": 7, "classes": 7, "closed": 7, "references": 1, "letters": 6, "N": 19, "bridges": 16, "sha": "30a4ac0846e5bcd6"}
    {"blasts": 2, "pip_depth": 7, "holes": 8, "classes": 8, "closed": 8, "references": 2, "letters": 7, "N": 25, "bridges": 19, "sha": "23ae55d497cedc72"}
    {"blasts": 3, "pip_depth": 8, "holes": 8, "classes": 8, "closed": 8, "references": 2, "letters": 6, "N": 35, "bridges": 26, "sha": "4f7ff702b695a422"}
    {"blasts": 4, "pip_depth": 9, "holes": 8, "classes": 8, "closed": 8, "references": 2, "letters": 6, "N": 47, "bridges": 33, "sha": "b875b083e9c1065d"}
== k32_series (branch; base identical)
    {"blasts": 0, "pip_depth": 1, "holes": 4, "classes": 3, "closed": 3, "references": 2, "letters": 3, "N": 4, "bridges": 3, "sha": "b78bd2d7883fc75c"}
    {"blasts": 1, "pip_depth": 2, "holes": 4, "classes": 3, "closed": 3, "references": 2, "letters": 1, "N": 6, "bridges": 4, "sha": "f6f18ca7c8a5797c"}
    {"blasts": 2, "pip_depth": 3, "holes": 8, "classes": 6, "closed": 9, "references": 2, "letters": 5, "N": 12, "bridges": 9, "sha": "46ca7720a5d3eb47"}
    {"blasts": 3, "pip_depth": 4, "holes": 8, "classes": 6, "closed": 9, "references": 2, "letters": 4, "N": 22, "bridges": 16, "sha": "b4b731ff3711b895"}
    {"blasts": 4, "pip_depth": 5, "holes": 8, "classes": 6, "closed": 9, "references": 2, "letters": 4, "N": 36, "bridges": 25, "sha": "f0862ba4f2c1e73d"}
    {"blasts": 5, "pip_depth": 6, "holes": 8, "classes": 6, "closed": 9, "references": 2, "letters": 4, "N": 58, "bridges": 40, "sha": "e7117badcd88d7b3"}
    {"blasts": 6, "pip_depth": 7, "holes": 13, "classes": 13, "closed": 18, "references": 2, "letters": 11, "N": 100, "bridges": 71, "sha": "784c50706f5ff3ca"}
    {"blasts": 7, "pip_depth": 8, "holes": 13, "classes": 13, "closed": 18, "references": 2, "letters": 10, "N": 174, "bridges": 124, "sha": "52a7d68878c29390"}
    {"blasts": 8, "pip_depth": 9, "holes": 47, "classes": 34, "closed": 45, "references": 6, "letters": 30, "N": 300, "bridges": 213, "sha": "f9ee18688961c7a5"}
    {"blasts": 9, "pip_depth": 9, "holes": 47, "classes": 34, "closed": 45, "references": 6, "letters": 30, "N": 521, "bridges": 371, "sha": "54668c7f6ea4aca3"}
    {"blasts": 10, "pip_depth": 9, "holes": 47, "classes": 34, "closed": 45, "references": 6, "letters": 30, "N": 897, "bridges": 638, "sha": "1c87a9493e4fbcee"}
== nested_outer_series (branch; base identical)
    {"blasts": 0, "pip_depth": 1, "holes": 12, "classes": 8, "closed": 8, "references": 4, "letters": 6, "N": 34, "bridges": 30, "sha": "9f56c2d62a3c0305"}
    {"blasts": 1, "pip_depth": 2, "holes": 12, "classes": 8, "closed": 8, "references": 4, "letters": 4, "N": 36, "bridges": 31, "sha": "6fb158ba8e93e0fc"}
    {"blasts": 2, "pip_depth": 3, "holes": 12, "classes": 9, "closed": 9, "references": 3, "letters": 6, "N": 40, "bridges": 34, "sha": "d11a6e3c9bd3b415"}
    {"blasts": 3, "pip_depth": 4, "holes": 12, "classes": 9, "closed": 9, "references": 3, "letters": 6, "N": 48, "bridges": 40, "sha": "7b9855e9f680c31c"}
    {"blasts": 4, "pip_depth": 5, "holes": 12, "classes": 9, "closed": 9, "references": 3, "letters": 6, "N": 60, "bridges": 48, "sha": "d0feaa8d3b51e2b8"}
    {"blasts": 5, "pip_depth": 6, "holes": 12, "classes": 9, "closed": 9, "references": 3, "letters": 6, "N": 80, "bridges": 62, "sha": "40797d6d690f632e"}
    {"blasts": 6, "pip_depth": 7, "holes": 12, "classes": 9, "closed": 9, "references": 3, "letters": 6, "N": 116, "bridges": 88, "sha": "6617db4f3743ac6b"}
    {"blasts": 7, "pip_depth": 8, "holes": 12, "classes": 9, "closed": 9, "references": 3, "letters": 6, "N": 176, "bridges": 130, "sha": "0a561405d0b88ab8"}
    {"blasts": 8, "pip_depth": 9, "holes": 12, "classes": 9, "closed": 9, "references": 3, "letters": 6, "N": 276, "bridges": 200, "sha": "824cf9ced5223b08"}
    {"blasts": 9, "pip_depth": 9, "holes": 12, "classes": 9, "closed": 9, "references": 3, "letters": 6, "N": 447, "bridges": 321, "sha": "87bfcc8be7388986"}
    {"blasts": 10, "pip_depth": 9, "holes": 12, "classes": 9, "closed": 9, "references": 3, "letters": 6, "N": 735, "bridges": 524, "sha": "b138d45a74782e7a"}
    {"blasts": 11, "pip_depth": 9, "holes": 12, "classes": 9, "closed": 9, "references": 3, "letters": 6, "N": 1217, "bridges": 863, "sha": "7295e5b190604b6f"}
```
`pip_depth` per step is equal on base and branch, and equal to the plan's table:

| series | blasts | `pip_depth` (base = branch) |
|---|---|---|
| k28 | 0..8 | `[1, 2, 3, 4, 5, 6, 7, 8, 8]` |
| k32 | 0..10 | `[1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9]` |
| nested_outer | 0..11 | `[1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9, 9]` |
| p3 | 0..9 | `[14, 15, 16, 17, 18, 19, 20, 21, 22, 23]` |
| inversion | 0..4 | `[5, 6, 7, 8, 9]` |

### Rubric 7: one 11-iteration `blast_zone` on nested_outer
- Call: `blast_zone(zone, num_iterations=11, fixed_point=[fp], min_separation=1e-4)` on `SERIES["nested_outer_series"].build()`, largest zone.
- Depth composed by `cases._pip_depth`.
```
base:   composed depth before 1 after 7 result.pip_depth <no attribute> completed 11 too_close 76 registry 486
branch: composed depth before 1 after 9 result.pip_depth {<FixedPoint A>: 9} completed 11 too_close 78 registry 498
```
- The branch's exemption INFO lines for that call, one per pip-bearing child (captured on `tanglepack.loom.Blast`):
```
  INFO blast: child (34, 35) ends at f^2(q0) of A (crossing 35), kept regardless of separation
  INFO blast: child (38, 39) ends at f^3(q0) of A (crossing 39), kept regardless of separation
  INFO blast: child (46, 47) ends at f^4(q0) of A (crossing 47), kept regardless of separation
  INFO blast: child (58, 59) ends at f^5(q0) of A (crossing 59), kept regardless of separation
  INFO blast: child (78, 79) ends at f^6(q0) of A (crossing 79), kept regardless of separation
  INFO blast: child (114, 115) ends at f^7(q0) of A (crossing 115), kept regardless of separation
  INFO blast: child (174, 175) ends at f^8(q0) of A (crossing 175), kept regardless of separation
  INFO blast: child (258, 259) ends at f^9(q0) of A (crossing 259), kept regardless of separation
```

### Other multi-iteration callers (plan section 5)
**E12's two-iteration cells** (`TWO_ITERATIONS`: k28 6 calls, k32 6, nested_outer 6, p3 10).
- Replayed by `two_iter.py`: `e11.build`, `area_cutoff = 1e-7`, `num_iterations=2`, E12's `SEPARATION` and `MAX_POINTS`, then a topology recompute after each call.
- k28, k32 and p3 are identical on base and branch on every call (pip depth, registry, N, bridges, holes, classes, `too_close`).
- nested_outer differs, from call 3 on:
```
           base                                                   branch
call 3: pip_depth=7 reg=116 N=116 bridges=88  too_close=1   |  pip_depth=7 reg=116 N=116 bridges=88  too_close=0
call 4: pip_depth=8 reg=256 N=256 bridges=186 too_close=20  |  pip_depth=9 reg=260 N=260 bridges=189 too_close=20
call 5: pip_depth=9 reg=506 N=506 bridges=362 too_close=77  |  pip_depth=9 reg=512 N=511 bridges=365 too_close=77
call 6: pip_depth=9 reg=902 N=901 bridges=641 too_close=166 |  pip_depth=9 reg=908 N=907 bridges=645 too_close=167
```
- Holes (12) and classes (9) are equal on every call.
- The branch reaches the (ii) stall depth 9 one call earlier and carries a few more crossings and bridges from there on.

**`scripts/henon_blast_period_3.py`** raises before any blast, on base and branch alike: `ValueError: Intersection 12 is not a strong-pip candidate. Candidates: [31, 32, 33]` at line 95, `T1.set_strong_pip(12)`.
- This is the stale pinned id its own comment warns about. It predates this change, and the script is outside the file list, so it was not touched and has no numbers.

## Behaviour changes

1. **Inside one multi-iteration `blast_zone` call:** a child ending on the chosen pip's orbit is kept even when `min_separation` would drop it. This is the intended N-a change. Measured effects:
   - nested_outer with 11 iterations reaches depth 9 instead of 7;
   - the `henon_p3_session` 3-iteration blast keeps the child ending at `f^4(q0)` (test 6a);
   - E12's nested_outer two-iteration cell differs from call 3 on (table above).
2. **New INFO lines:** one per pip-bearing child that the separation test was skipped for (`min_separation` set). Single-iteration blasts log them too: there is typically one such child per blast. No other log line was added or changed.
3. **New data:** `BlastResult.pip_depth`.
4. **None** on any law case (rubric 4) and on any step of the five `SERIES` (rubrics 5 and 6).
5. **Structural, no observable effect:** the child loop runs after all frontier bridges of an iteration are mapped, instead of interleaved per parent. The filter reads nothing `iterate_bridge` changes for other parents (`seen` already holds every frontier bridge, and `in_zone` / `fp_ok` are pure). Rubrics 4 and 5 bear this out.

## My reading of the rubric

1. **Pass.** The branch suite passes; the only xfail is the `KNOWN_ISSUES` `orientation_reversing` entry.
2. **Pass.** On base, 6(a) and 6(c) fail; 6(b), the restated proximity test (both variants) and `test_pip_orbit_never_shallower` on all five series pass. All of these pass on the branch.
3. **Pass.** Exactly the five files, 264 changed lines (≤ 300).
4. **Pass.** All 11 law-case fingerprints are identical (with the per-process `products` ids and repr addresses dropped).
5. **Pass.** Every `SeriesStep` field of all five series is byte-identical.
6. **Pass.** `pip_depth` per step is identical on base and branch and matches the plan's table, including the k28, k32 and nested_outer plateaus.
7. **Pass.** Composed depth 7 on base; `result.pip_depth[fp] == 9` on the branch.
8. **Pass.**
   - The exemption sits after the exterior check and before the separation test.
   - The new code composes `iterate_table[x, 1]` and never calls `forward_depth` (only the docstring names it, to say why it is not used).
   - One INFO line per exempted child.
9. **Pass.**
   - The `blast_zone` docstring says `min_separation` "never decides the depth of `T(q0)`" and names the pip-orbit exemption.
   - The `BlastResult` docstring documents `pip_depth`, `too_close` and `already_known`.
10. **Pass, with notes.**
    - The only new log lines are the exemption INFO lines.
    - The behaviour changes are the multi-iteration ones of section 5, measured above (nested_outer 11 iterations, and E12 nested_outer two-iteration from call 3).
    - The period-3 script cannot be measured: it fails identically on base.

## Notes for the test engineer
- `_pip_orbits` resolves the pip by `registry.find`, not by the remembered `Intersection.id` (see step 1). The test side does the same (`cases._pip_depth`, and the helpers in both test files).
- The restated proximity test is parametrized, so its node ids are `test_min_separation_drops_close_bridges[pip]` and `[no_pip]`.
