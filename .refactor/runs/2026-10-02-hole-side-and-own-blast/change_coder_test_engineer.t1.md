# Change: hole-side-and-own-blast

*Coder, 2026-10-02. Plan: `artifacts/plan_planner_human_2.md`. Verdict: `plan_ok_planner_human.json` (approved, no notes; the plan's defaults are used: a combinatorial direct-hole rule, and the one CLAUDE.md sentence updated).*

- **Base commit:** `231520415a42badc257b63abe69b51fcebb137c4` (`2315204`, `main`)
- **Branch:** `refactor/hole-side-and-own-blast`
- **Worktree:** `/home/dezhu/Development/TanglePack/tangle-pack-hole-side-and-own-blast`
- **Base worktree** (detached at `2315204`, clean, used for the at-base checks): `/home/dezhu/Development/TanglePack/tangle-pack-hole-side-and-own-blast-base`

I did not touch the author's working tree or their uncommitted changes.

## Commits

```
5d51413 Hole sides: propagated holes read on their image sub-arc, direct holes from crossing signs
cf75bc0 build_nested: each zone blasts only its own fixed point
f68c706 StablePartition: reflow the side note, split a long line
9bbfbda Tests: deep period-3 hole sides; nested inner words and blast order
490c6e7 CLAUDE.md: the period-3 side trip is fixed, deep runs still unreliable
```

## `git diff --stat 2315204`

```
 CLAUDE.md                                  |  2 +-
 src/tanglepack/examples/henon_cases.py     | 22 ++++++----
 src/tanglepack/topology/StablePartition.py | 38 +++++++++++++-----
 tests/test_higher_period_cartoon.py        | 64 ++++++++++++++++++++++++++++++
 tests/test_stable_partition_period3.py     | 25 ++++++++++++
 5 files changed, 134 insertions(+), 17 deletions(-)
```

Only the five files in the plan were changed.

## What changed

**A. Propagated hole side (`StablePartition.py`).** `_bridge_side_of(trellis, bridge, point, span=None)` takes an optional unstable-cdist `span`. When a span is given, the oriented polyline is cut down to the image sub-arc (with one node of margin each side, found by `bisect` over the sorted node cdists) before `_arc_side_of` runs. `_punch_in_bridge` is the only caller that passes it: `_bridge_side_of(trellis, bridge, carried, span)`. The openings still come from the side through `_hole_openings`, which is unchanged (Not doing 3).

**A'. Direct hole side (`punch_holes`).** The side now comes from the crossing signs instead of `_bridge_side_of`:
```python
first, second = sorted((near, far), key=lambda x: x.unstable_cdist)
outward = 1 if second.stable_cdist > first.stable_cdist else -1
bridge_side = "left" if first.crossing_sign * outward > 0 else "right"
```
Why it holds: the bridge leaves `first` along u+. The pair's stable segment runs from `first` along s+ when `second` is outward. `crossing_sign = sign(cross(u+, s+))`, so s+ (and with it the lobe) is on the left of u+ exactly when `crossing_sign > 0`. The comment in the code says this. The openings are still `inward=True` by construction. The module Dev Note on "Side" now describes both rules.

**B. `build_nested` (`henon_cases.py`).** `fixed_point=[zone.fixed_point]` is used in the blast loop. The `inner_blasts` docstring now says the order does not matter because each zone's blast touches only its own tangle. The Dev Notes line about 16 steps now says the trip at 15/16 steps and at 6+ blasts is fixed but those runs are unreliable. I also added a Dev Note bullet explaining why each zone blasts only its own fixed point.

**CLAUDE.md:** one sentence. "16 steps trips `check_holes_share_bridge_side` (open)." became: the trip at 15–16 steps and 6+ blasts is fixed (2026-10-02, with both rules named), but those runs are unreliable (an unreachable class, a virtual `new1`; open).

**Tests:**
- `tests/test_stable_partition_period3.py`: `test_deep_p3_holes_share_bridge_side_through_symbolic_dynamics[6_blasts | 15_steps]`. It builds `build_period3(**kwargs)`, runs `check_holes_share_bridge_side` on the trellis holes, and asserts that `symbolic_dynamics().classes` is non-empty. It only checks that the run completes, not that the result is reliable. Marked `slow` and `regression`, about 0.2 s each.
- `tests/test_higher_period_cartoon.py`:
  - The helper `_tangle_words(build, fixed_point)` maps each class of one tangle to its word. A class is named by its two homotopy elements as `(orbit_code, short_text)`, with no fixed-point letter. A token is `(class name, direction)`. Classes are selected by `bridge_class.source.branch_key[0] is fixed_point`. I did not use `entry.tangle`: it indexes the session's fixed-point order, which is `[fp3, fp1]` here, not `TangleBuild.fixed_points`.
  - `test_nested_inner_words_are_the_period3_words`: the inner words of `build_nested(outer_blasts=2, inner_blasts=4)` equal those of `build_period3(blasts=4)`, all 9 classes (inert ones included). About 0.6 s.
  - `test_nested_blast_order_does_not_matter`: `build_nested(4, 4)`, run under caplog, logs no `skipping bridge`. A second run does `build_nested(outer_blasts=4)`, then 4 inner blasts through `session.blast_zone(inner_zone, 1, fixed_point=[inner_zone.fixed_point], min_separation=1e-4)` with a repin after each, then re-partitions. For both fixed points, the words and the bridge counts are equal. About 1.8 s.

## Suite

- Branch, full suite (`MPLBACKEND=Agg env/bin/python -m pytest -q`, run in the worktree, importing the worktree's `src`): **786 passed, 1 skipped in 51.62s**. That is the 782 from base plus 4 new.
- After commits A and B only, before any new tests: 782 passed, 1 skipped. So no existing test changed result.

## Behaviour changes

1. **Direct-hole `bridge_side`:** no change in any run I compared (see the side comparison under rubric item 9). The combinatorial rule did not disagree with the suite anywhere, so I did not need to stop over Question 1.
2. **Propagated-hole side:** the period-3, period-3 at 4 blasts and nested default runs give identical `(bridge_side, openings)` on every hole at base and on the branch. Period 3 at 6 blasts, 15 steps and 16 steps now pass `check_holes_share_bridge_side` and reach `symbolic_dynamics()` with `is_reliable == False`, as the plan expected.
3. **`build_nested()` defaults:** the words are identical to base. `blast_sizes` goes from `[32, 7]` to `[3, 3]`, because the outer blasts no longer count blasted period-3 bridges. The period-3 bridge count goes from 31 to 27 and the period-1 count stays at 7. The plan named the 31 -> 27 change but **not the `blast_sizes` change**, which follows from the same cause. Registry ids in the nested run shift (the period-1 holes `(37,36),(39,40)` became `(35,34),(37,38)`), with the same sides and openings.
4. **Nested runs with inner blasts:** no "blast: skipping bridge" warning, and the result does not depend on blast order. At base, the new nested tests fail on `check_holes_share_bridge_side` (origin `(51, 52)`) and log that warning.

## Rubric, my reading

1. **Pass.** 786 passed, 1 skipped (at least 782 plus 4 new).
2. **Pass.** At `2315204`, `tests/test_stable_partition_period3.py tests/test_higher_period_cartoon.py tests/test_loom_blast_restore.py` gives 24 passed. On the branch they pass inside the full suite.
3. **Pass.** I ran the branch's test file against the base worktree's `src`, by copying it into a temporary `tests_branch/` in the base worktree (now removed; the base worktree is clean). `[6_blasts]` fails with `AssertionError: Holes of origin (41, 42) disagree on bridge_side: iterate 0 is right, iterate -3 is left; ... -5 ...; ... -6 ...`. It passes on the branch.
4. **Pass, with one deviation.** `[15_steps]` fails at base with the same assertion, but the origin is **`(21, 35)`, not the plan's `(12, 33)`**. Iterates -3/-5/-6 match the plan. My guess is that registry ids depend on the exact build path (the analyst may have built 15 steps differently); this is unverified. It passes on the branch.
5. **Pass.** The stat above lists exactly the five planned files.
6. **Pass.** `grep -n "_bridge_side_of(" StablePartition.py` shows only the definition (`:789`) and `_punch_in_bridge` (`:1547`, which passes `span`). No other caller remains: the direct-hole call is gone.
7. **Pass.** The direct-hole `bridge_side` in `punch_holes` uses only `crossing_sign` and the stable cdists.
8. **Pass.** No `fixed_point=fixed_points` is left in `henon_cases.py`, and the "so the outer blasts do not consume its bridges" wording is gone.
9. **Pass, with the extra item in Behaviour 3.**
   - The k=10 and k=2.8 fixture tests pass.
   - `build_period3(blasts=4)` gives `a -> b`, `b -> c`, `c -> a u^-1 w^-1`, with blast sizes `[16, 2, 2, 2]` at base and on the branch.
   - `build_nested()` gives the same words as base (`a->b, b->c, c->a u^-1 w^-1, d->d uu^-1 e^-1, e->f, f->d uu^-1 d^-1`, and the inert words).
   - Period-3 bridges go from 31 to 27 as stated, but `blast_sizes` goes from `[32, 7]` to `[3, 3]`, which section 5 did not state.
10. **Pass.** `test_nested_blast_order_does_not_matter` asserts that `build_nested(4, 4)` logs no `skipping bridge` warning.

## Notes for the test engineer

- The new tests are marked `slow`, as the existing nested and period-3 tests are, so they run in the default suite. Together they take about 3 s.
- `_tangle_words` compares inert classes too. If a stricter or looser comparison is wanted (for example, active classes only), the helper is the one place to change it.
- The 16-step case also runs on the branch now (unreliable). It is not pinned by a test, because the plan names only 15 steps.
