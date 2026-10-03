# Outcome (second pass, after figures): hole-side-and-own-blast

*Auditor, 2026-10-02. Branch `refactor/hole-side-and-own-blast` in worktree `../tangle-pack-hole-side-and-own-blast`. The branch head is now `f8a6d73`: the figure maker added three proof scripts on top of the reviewed head `3b33baa`. Base `2315204`. The review passed on its first round. This pass replaces my first outcome.*

Your `main` tree is still dirty: `.claude/` edits, `artifacts/`, `.refactor/` and two figure folders. Nothing on the branch touches them. The branch is not merged. In the branch worktree, `figures/refactor_proof/` is untracked, as intended.

## What the figures confirm and what they change

I opened every summary PNG myself, plus the panels that carry the before/after comparisons. Overall, **the figures confirm the first outcome and add three things to it.**

**Confirmed:**
- **Nested equals each tangle alone, and blast order does not matter** (`nested_grid.png`, right-hand grid). On the branch, all 20 cells (outer {0,2,4,6,8} × inner {0,2,4,6}) read the same as each tangle alone. The three order checks at (2,2), (4,4) and (8,4) all say "same".
  - I compared `panels/branch/nested_o4_i4/itineraries.png` with `nested_o4_i4_outer_first/itineraries.png`. Classes, iterated itineraries, words and refined words are identical row by row. Only the title differs.
  - The period-3 rows match `p3_b4`: `a->b, b->c, c->a u^-1 w^-1`. The period-1 rows read `d -> d uu^-1 e^-1, e -> f, f -> d uu^-1 d^-1`, with `d` refined into `d_1, d_2`. These are the k=2.8 two-blast words, relettered.
- **Fix (A) fixes the actual bug** (`hole_sides.png`). Exactly 3 of the 40 propagated holes read differently under the old and new rules, at iterates −3, −5 and −6 of origin (41, 42). That is the same origin and the same iterates as the base `AssertionError`, quoted verbatim in `panels/base/p3_b6/ERROR.txt`.
  - In every panel, the old reference vertex lies outside the hole's own sub-arc, on a neighbouring fold: cdist 887 for a sub-arc spanning 707.7–745.0, 369 for 294.1–309.6, and 238 for 189.6–199.6.
  - The new vertex lies inside the sub-arc and gives `right`, which agrees with iterate 0.
  - Panels −5 and −6 show the two folds clearly apart. In panel −3 they nearly overlap, so the picture alone does not show it; the cdist label is what settles it.
- **Period 3 alone is stable for 0 to 5 blasts and breaks at 6** (`p3_blast_sweep.png`). For 0 to 5 blasts, base and branch draw the same line: 12 holes, 9 classes, 3 active. At 6 blasts the branch jumps to 42 / 23 / 22, and the word for `c` changes at every later blast. Base raises at 6, 7 and 8 blasts.
  - So the open failure of your "never contradict" criterion is confirmed. It is still open.

**New, from the figures:**
1. **Before the change, outer blasts alone silently corrupted period 3.** At base, (6,0) and (8,0) run without error, but period 3 has collapsed. Compare `panels/base/nested_o6_i0/itineraries.png` with the branch version:
   - classes `a` and `b` are unresolved;
   - `c -> b^-1`;
   - every other period-3 class is a trivial loop;
   - there are 8 period-3 classes instead of 9.

   My first outcome only said the old recipe "depends on order". This is worse: wrong answers with no error raised. Fix (B) is what removes it.
2. **Inner = 6 inside nested still matches period 3 alone.** Both are unreliable, but they agree. So even the broken deep regime is not made worse by nesting: the open problem belongs to period 3, not to nesting.
3. **The p3 sweep was run at `min_separation=1e-4`** (`words_grid.py:8,155`, chosen to match the nested recipe). `build_period3` defaults to `1e-5`, which is what the tests and my first probe used. Both settings give the same jump at 6 blasts (12 → 42 holes), so nothing changes, but the figure is not at the recipe default.

**Two small corrections to the figure report:**
- Its summary says "at base, 13 of 20 cells raise". `grid_base.json` has **12** of the 20 grid cells raising. The 13th failure is the extra outer-first (4,4) run, which is not a grid cell. The report's own detailed list ("all 11 cells with inner ≥ 2 and outer ≥ 2, plus (0,6)") says 12.
- Fix (A′), the direct-hole rule, has **no figure of its own**. The only figure evidence is the outer-first (4,4) crash at base, plus base `p3_7`/`p3_8` in the JSON. My point 2 below therefore still stands.

## Your question, answered now

**Yes, with one exception. Each tangle alone, and the nested case in any order, now run through to `symbolic_dynamics()` with every asserted invariant holding. The nested words equal each tangle's words alone.**

- The exception is **period 3 past 5 blasts** (or at 15 or 16 unstable steps). It runs, but it is unreliable and its words change with every further blast. Nesting inherits exactly that and adds nothing.
- Reliable regime: outer ≥ 2 and inner ≤ 5. With outer = 0, the outer class `e` stays unresolved, as documented.
- The figure maker's 20 cells add inner = 6 and outer = 6 to my earlier 16 (which had inner and outer in {0,2,4,8}). There are no new failures.

## What changed

The library change is unchanged since the review. There are 6 commits up to `3b33baa`, touching the five planned files:

| File | Change |
|---|---|
| `StablePartition.py` | +29 / −9 |
| `henon_cases.py` | +15 / −7 |
| `CLAUDE.md` | one sentence |
| two test files | +188 / −1 |

Commit `f8a6d73` adds only `scripts/refactor_proof/hole-side-and-own-blast/{words_grid,hole_sides,panels}.py` (+537 lines). It touches no `src/` or `tests/` file. I checked this with `git diff --stat 3b33baa f8a6d73`.

The most important hunk is (A), `StablePartition.py:789-823` and `:1547`:

```python
# before
bridge_side = (
    _bridge_side_of(trellis, bridge, carried) if carried is not None else None
)
# after: the side is read only on the hole's image sub-arc
if span is not None:
    nodes = bridge.get_point_array(return_nodes=True)
    cdists = sorted(node.cdist for node in nodes)
    lo = max(bisect_left(cdists, span[0]) - 1, 0)
    poly = poly[lo : bisect_right(cdists, span[1]) + 1]
...
_bridge_side_of(trellis, bridge, carried, span)
```

(A′), at `StablePartition.py:349-356`:

```python
bridge_side = "left" if first.crossing_sign * outward > 0 else "right"
```

(B), at `henon_cases.py:264`:

```python
fixed_point=[zone.fixed_point]
```

## Shorter or clearer

- **Nothing got shorter.** The source grew by 28 net lines, mostly docstrings and Dev Notes.
- **(A′) reads like the maths in `row_of_end`.** It replaces a geometric measurement that its own comment called "estimated".
- **The `build_nested` docstring gives the real reason** blast order does not matter.

## Behaviour changes, all of them

1. **Propagated-hole sides:**
   - No change in any run that passed before.
   - In the runs that failed, the holes at iterates −3, −5 and −6 flip to agree with iterate 0. That is exactly 3 of 40 at p3 b6 (figure 3).
2. **Direct-hole sides:** no change on any fixture (pinned hole by hole on k10, k28, k28 with two blasts, and p3).
3. **`build_nested()` defaults:**
   - same words;
   - period-3 bridges go from 31 to 27;
   - `blast_sizes` goes from `[32, 7]` to `[3, 3]` (not listed in the plan, but disclosed, with the same cause);
   - registry ids shift.
4. **Runs that newly complete:** period 3 at 6 or more blasts, and at 15 or 16 steps. Also 12 nested cells that raised at base. Period 3 is unreliable in all of these.
5. **Runs that newly give correct answers:** nested with 6 or 8 outer blasts. These ran at base but gave wrong period-3 words (new finding 1). Nested at (2,2) is no longer order-dependent and no longer logs "skipping bridge".

## How it was verified

- **Full suite at `f8a6d73`, run by me this pass:** 795 passed, 1 skipped (50 s). Base gives 782 passed, 1 skipped (my first pass).
- **The new tests fail on the base code.** The test engineer ran the branch's tests against base: 5 failed and 25 passed, and the 5 failures are exactly the bug-fix tests. The reviewer reproduced this.
- **Each new test catches the change it covers.** Undoing (A) turns 2 tests red. Flipping (A′) gives 20 failures and 219 errors. Reverting (B) gives 3 failures.
- **Figures and data:**
  - I read all three summary PNGs and five panels: both order-variant tables, the base and branch `nested_o6_i0` tables, and the `p3_b6` cartoon.
  - I checked `grid_base.json` and `grid_branch.json`: the 12 base exceptions, and branch p3 reliable for 0 to 5 blasts and unreliable for 6 to 8.
  - I read the base `ERROR.txt`.
- **I did not re-run the figure scripts.** The PNGs and JSON are the figure maker's output. I checked them against each other and against my own first-pass probe, and they agree everywhere they overlap.

## What did not go to plan

- **Rubric 5:** `git diff --stat 2315204` now lists 8 files, because of the proof scripts that were added on instruction. Against `3b33baa`, the reviewed head, it is still exactly the five planned files. If you merge, decide whether `scripts/refactor_proof/` should come along. The scripts hard-code `../tangle-pack/env/bin/python` and the worktree layout.
- **Rubric 4** named the 15-step origin as `(12, 33)`. Registry ids are not reproducible between builds; the failing origin has been seen as `(12,33)`, `(21,35)` and `(3,17)`.
- **The test diff is about 230 lines,** against the planned ~80, because of the 9 pinning tests.
- **The 16-step case is not pinned by any test.**

## Where I am least sure this is what you meant

1. **Fix (B) lives in the recipe, not the library.** `Blast.blast_zone` still blasts whatever fixed points the caller passes. The figures make the cost of passing both visible: at base, outer blasts silently corrupt period 3. Any script other than `build_nested` that passes both saddles will hit it again. If "a zone blasts only its own tangle" is meant to be a rule, it belongs in `Blast.py` or `TangleSession.blast_zone`.
2. **(A′) is a rule, and nothing independent confirms it.** It matches the old measurement wherever that measurement was right. No figure or test shows it is right where the two disagree, on the tiny +1/+2 lobes in deep runs. The only evidence there is that the side invariant now holds.
3. **`CLAUDE.md` calls the deep period-3 runs "unreliable", which understates the problem.** The sweep figure shows them breaking your stability criterion outright:
   - 30 new holes, all on the R rows (see the `p3_b6` cartoon);
   - the 0-blast classes are renumbered;
   - a different word for `c` at 6, 7 and 8 blasts.

   That is the next thing to investigate, not a reliability footnote.

## Records updated

- `.refactor/decisions.md`: appended your "figures" verdict and what it led to.
- `.refactor/knowledge/testing.md`: the proof scripts, the `min_separation` caveat and the base/branch grid facts.
- `.refactor/knowledge/topology.md`: what the `p3_b6` cartoon shows about the open deep period-3 issue.
