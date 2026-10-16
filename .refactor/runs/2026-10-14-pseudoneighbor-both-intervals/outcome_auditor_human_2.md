# Outcome — B2 `pseudoneighbor-both-intervals`

*2026-10-14. Auditor: analyst, `outcome` job. Branch `refactor/pseudoneighbor-both-intervals` @ `9cd29d5`, base `ba23a1f`, worktree `../tangle-pack-pseudoneighbor-both-intervals` (clean). Not merged; the merge is yours. Your main tree still has its uncommitted `docs/`, `.refactor/` and `artifacts/` changes. The branch never touched them.*

## The answer now

**Yes, the code now tests both open intervals against X, as the definition does.** The PDF keeps a consecutive window pair `(x0, x1)` only when neither `W^U(x0, x1)` nor `W^S(x0, x1)` holds a point of `X = ∪ M^n(N)`. Before this branch, the code checked the unstable interval against every iterate. For the stable interval it relied on consecutiveness, so only crossings of N were excluded and no iterate was ever looked at.

Now one function, `_interval_is_clear(trellis, x0, x1, stability, candidates, collision_rtol)` (`Pseudoneighbor.py:343`), tests either interval. It uses the same three readings in the same order: the table landing by id, `r` read in its own frame against the pair's registered images, and β-scaling with the 2-D collision guard. `compute_pseudoneighbors` keeps a pair only when both calls return True (`:229`).

On the stable side the function uses:
- the stable branch's own fixed point, cycle and `per_step_beta("stable")` (β < 1);
- the same candidates (`own_intersection_ids`, the cut N);
- a mirrored D7 guard: the 2-D collision skip applies only when the candidate's *unstable* side belongs to the pair's unstable fixed point.

On every case we have, the result is the same as before. The stable clause never fires on a real build.

## What changed

The four planned files, and nothing else (`git diff --stat ba23a1f..9cd29d5`, verified):

| file | ± |
|---|---|
| `src/tanglepack/topology/Pseudoneighbor.py` | +128 / −110 |
| `tests/unit/topology/test_pseudoneighbor.py` | +136 / −8 (incl. the test engineer's commit `9cd29d5`) |
| `stress/numerics/e13_both_intervals.py` | new, 198 |
| `stress/numerics/results/e13_both_intervals.jsonl` | new, 14 cells |

The most important hunk is the head of the function and its m-range.

Before (`ba23a1f`):
```python
u_key = ix0.manifold_a_key
if u_key is None or u_key != ix1.manifold_a_key:
    return False
owner = u_key[0]
pos_map, k = owner.branch_position_map("unstable")
beta = owner.per_step_beta("unstable")
lo, hi = sorted((ix0.unstable_cdist, ix1.unstable_cdist))
...
same_stable = r.manifold_b_key[0] is ix0.manifold_b_key[0]
m_lo = math.floor(math.log(lo / c) / log_lap)
m_hi = math.ceil(math.log(hi / c) / log_lap)
for m in range(m_lo, m_hi + 1):
```
After (`9cd29d5`):
```python
arc_key = attrgetter("manifold_a_key" if stability == "unstable" else "manifold_b_key")
other_key = attrgetter("manifold_b_key" if stability == "unstable" else "manifold_a_key")
cdist = attrgetter(f"{stability}_cdist")
key = arc_key(ix0)
if stability == "stable":
    assert key == arc_key(ix1), "the walk steps along one stable branch"
elif key is None or key != arc_key(ix1):
    return False
owner = key[0]
pos_map, k = owner.branch_position_map(stability)
beta = owner.per_step_beta(stability)  # < 1 for the stable arc
lo, hi = sorted((cdist(ix0), cdist(ix1)))
...
same_other = other_key(r)[0] is other_key(ix0)[0]          # D7, and mirrored
bounds = (math.log(lo / c) / log_lap, math.log(hi / c) / log_lap)
for m in range(math.floor(min(bounds)), math.ceil(max(bounds)) + 1):
```
The rest of the body only swaps `.unstable_cdist` for `cdist(...)`. A local `point()` restores `(unstable, stable)` order for `same_point`, and the three debug wordings become one.

## What got shorter or clearer

- **Not shorter.** `Pseudoneighbor.py` went from 471 to 489 lines, and the function from 143 to 147. The plan estimated "net shorter"; it is not. The growth is in the Dev Notes and docstrings, which now describe both intervals and the mirrored D7 guard.
- **Clearer.**
  - The Dev Notes state the definition once: neither open interval may hold a point of X.
  - "Empty by construction" is gone (verified by grep).
  - The β note says `per_step_beta(stability)` and never uses `1/per_step_beta("unstable")`.
  - The check-set note explains why the candidates are never the whole registry. With the registry as candidates, the law case `inversion` would lose its references; the test engineer's break m3 shows this.
  - One debug line names the stability, the witness, the step count and the reading. e13 parses it.

## Behaviour changes

1. **Results: none on the corpus.**
   - Law cases: the fingerprint, reference pairs, all pseudoneighbor pairs and every class word are identical on base and branch for all 11. This was checked by two independent base-vs-branch scripts (coder and test engineer). I did not run a third.
   - e13: all 14 cells (7 E12 cases × 1e-7 / 1e-5, at E12's depths) are `ok`, with 0 stable fires over 182 blast steps, `before == after`, and `after ==` the library's references at every step. Verified: I re-read the committed jsonl.
2. **Seven unit-test fixtures changed in data only.** Names and assertions are unchanged. Their data put a crossing of N whose stable orbit visits the window strictly inside a pair's stable interval. That is a real puncture under the definition you decided on, so no correct implementation could keep those pairs. The seven tests:
   - `direct_intersection_punctures_interval`
   - `iterate_punctures_interval`
   - `walk_stops_at_reference_window_end`
   - `heteroclinic_candidate_ignored`
   - `table_linked_deep_iterate_does_not_disqualify`
   - `table_linked_landing_uses_registered_cdist`
   - `unlinked_landing_is_read_against_the_pairs_images[below_c]`

   The fixes:
   - Five probes moved to the saddle's other stable branch (`OTHER_STABLE = 1`).
   - `walk_stops` now has `below` at M(r_end)'s cdists.
   - `deep_iterate` flips its drift from `1 + 2·rtol` to `1 − 2·rtol`.

   Verified: the literal base test file run against branch `src` gives **7 failed, 24 passed**, exactly these seven. The adjusted file passes on base `src` as well (test engineer).
3. **Debug message format.** The three puncture lines are now one: `Pair (x0, x1): <stability> interval punctured by r at M^n, <table|own frame|scaled> reading (cdist v in (lo, hi))`. This is DEBUG level only, and no test pinned the old wording.

## What was verified, and how

All checks below were run by me on `9cd29d5` unless marked otherwise.

- **Full suite:** `MPLBACKEND=Agg env/bin/python -m pytest -q -rxX` gives **1183 passed, 38 skipped, 5 deselected, 1 xfailed**. The only xfail is `orientation_reversing_case_builds`.
- **Scope:** the diff stat shows the four files; the worktree is clean.
- **Diff:** read in full.
  - β is `per_step_beta(stability)` (`:394`).
  - The candidates are still `own_intersection_ids` (`:197`).
  - The guard is `same_other` (`:424`).
  - The m-range takes `min/max` (`:430`).
  - No `registry` or `resolved_ids` is used as the candidate set.
- **Rubric 6 (a stable interval wider than one lap, by hand):** I ran the test engineer's `r6.py`. With the candidate the call returns False; without it, True.
- **The original fixtures under the new rule:** the base test file on branch `src` gives 7 failed, as listed above.
- **e13 jsonl:** 14 cells, all `ok`, 0 lost, 0 stable fires, 0 library mismatches. It was generated at `530c1ba`, whose `src` is the final one; `9cd29d5` only adds a test.
- **One extra probe:** I called the new function on every consecutive pair of every stable branch of `p3_6_blasts`.
  - On the two non-pip branches, 6 pairs each have a clear unstable interval and a punctured stable one.
  - These are the "orphan window pairs" of the DEC-19 gap note. The new clause now rejects them, as the definition says.
  - `compute_pseudoneighbors` walks only the pip branch's window, so this changes nothing today. It does matter if gap-filling (DEC-19) ever walks other branches.
- **Believed, not re-run by me:** the law-case fingerprint equality (two independent agents' scripts) and the test engineer's mutation table (m1–m6, a1–a6).

## What did not go to plan

- **Rubric 2's literal wording failed.** It required every base test to pass on both base and branch, and the seven fixtures above did not. The analyst had measured only real builds, so nobody foresaw that hand-built fixtures would contradict the definition. The reviewer accepted the data-only rewrite.
- **The function is not net shorter:** +128 / −110 in the file.
- **The test engineer added a fourth test**, `test_collision_on_another_unstable_manifold_still_punctures_stable_interval[own_orbit|other_unstable]`, which the plan did not ask for. Without it, dropping the mirrored D7 guard went undetected (break m4).
- **`tests/cases.py` was correctly left off the list:** it has no `KNOWN_ISSUES` entry for `reference_pairs_clear_of_X`.

## e13 and DEC-25

**No new input.** The inversion series still drops from 2 references to 1 at blast 6 at both cutoffs: `(23,24),(26,27)` becomes `(56,57)`, holes go 8 → 5 and seed classes 8 → 6. `before == after` at every step, so the stable clause plays no part. The drop is decided by the unstable clause and the window, as before.

One more fact from e13: N stops growing on most cells (k28 at 143 from blast 7, k32 at 300 from blast 8, and so on), so the later blasts measure the same N again. That is B3's concern.

## Where I am least sure it was done as you meant

1. **The fixture rewrites.** They are correct by the definition, and the test engineer showed they keep the originals' teeth (breaks a2–a6). But they changed the meaning of seven of your tests' data without your ruling: five probes now live on a different stable branch. If you would rather those tests probe the window's own branch, each one needs a crossing that is an iterate of a window end, not a move.
2. **The m-range `min/max` is guarded only by an uncommitted hand check.** Every interval the window walk asks about is at most one lap, so no committed test or real build fails if the bounds are un-swapped (break m1). It is a one-function test if you want it committed; the plan deliberately left it out.
3. **Two older tests have no teeth, and CLAUDE.md now disagrees with the code's own note.**
   - `walk_stops_at_reference_window_end` and `heteroclinic_candidate_ignored` pass even when the rule they are named for is broken. This was identical on base and was left alone.
   - CLAUDE.md says "a crossing past a cut is a point of X, not of N". The new Dev Note says the opposite, and correctly: "a registered crossing past the unstable cut is no point of X unless its orbit meets T". The plan deferred this wording fix (B12), but until it lands a reader of CLAUDE.md would widen the candidates and lose inversion's references.

## Files

- `/home/dezhu/Development/TanglePack/tangle-pack-pseudoneighbor-both-intervals/src/tanglepack/topology/Pseudoneighbor.py`
- `/home/dezhu/Development/TanglePack/tangle-pack-pseudoneighbor-both-intervals/tests/unit/topology/test_pseudoneighbor.py`
- `/home/dezhu/Development/TanglePack/tangle-pack-pseudoneighbor-both-intervals/stress/numerics/e13_both_intervals.py`
- `/home/dezhu/Development/TanglePack/tangle-pack-pseudoneighbor-both-intervals/stress/numerics/results/e13_both_intervals.jsonl`
- Updated: `.refactor/decisions.md` (B2 outcome entry), `.refactor/knowledge/topology.md` (B2 outcome section; stale `_unstable_interval_is_clear` references and the orphan-pair note corrected), `.refactor/knowledge/testing.md` (B2 outcome audit).
