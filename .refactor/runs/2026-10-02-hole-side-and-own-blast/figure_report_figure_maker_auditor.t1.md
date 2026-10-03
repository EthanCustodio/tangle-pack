# Figure report: hole-side-and-own-blast

*Figure maker, 2026-10-02. Branch `refactor/hole-side-and-own-blast`, worktree `../tangle-pack-hole-side-and-own-blast`. Base `2315204`, detached worktree `../tangle-pack-hole-side-and-own-blast-base`, which already existed. Your main tree was not touched.*

## In short

- **The nested words now match each tangle on its own.** On the branch, every nested cell from outer {0,2,4,6,8} × inner {0,2,4,6} runs. In every cell, the period-1 words equal period 1 alone and the period-3 words equal period 3 alone. Blasting in either order gives the same result at (2,2), (4,4) and (8,4).
- **The old code failed in most of the grid.** At base, 13 of the 20 cells raise. 3 more give period-3 words that differ from period 3 alone, or depend on blast order.
- **Period 3 is stable up to 5 blasts, then changes.** On the branch it gives identical, reliable words for 0 to 5 blasts. At 6 blasts it now runs, but the hole count jumps from 12 to 42, the classes are renumbered and the words change at every further blast. This is the open issue; the base code raises there.

## Commit on the branch, and the effect on rubric 5

There is one new commit, `f8a6d73`, which adds three scripts under `scripts/refactor_proof/hole-side-and-own-blast/`. My instructions put them on the branch.

This means `git diff --stat 2315204` now lists these three files on top of the five planned ones. When you re-audit rubric 5, compare against `3b33baa` (the reviewed head), or exclude `scripts/refactor_proof/`.

The figures are **not** committed. They sit untracked under `figures/refactor_proof/hole-side-and-own-blast/` in the branch worktree, because each panel also writes `.svg` and `.fig.pickle` files.

All three scripts use `../tangle-pack/env/bin/python` with `MPLBACKEND=Agg`. Each one puts the tree it runs into `sys.path` and asserts that `tanglepack` was imported from that tree. The env's editable install points at your main tree, so without this the wrong code would run.

## Figures

All paths are under `/home/dezhu/Development/TanglePack/tangle-pack-hole-side-and-own-blast/figures/refactor_proof/hole-side-and-own-blast/`.

### 1. `nested_grid.png`: nested vs each tangle alone, before and after

**What it proves.**
- After the change, at every (outer, inner) cell, each tangle of the nested case reads exactly as it does alone, and the blast order does not matter.
- Before the change, most cells raised or disagreed.

**What to look at.** There are two grids: base on the left, branch on the right. Rows are outer (period-1) blasts {0,2,4,6,8}; columns are inner (period-3) blasts {0,2,4,6}. Each cell says, for P1 and P3, whether that tangle's words equal the same tangle alone at the same blast count. Words are compared spelled out in element names (orbit code plus short name, no fixed-point letter) with directions, because letters differ between sessions.

Colours:
- green: both tangles match alone, and both are reliable;
- pale yellow: both match, but not reliable;
- red: a tangle reads differently from alone, or the blast order matters;
- grey: the run raised.

Three cells also show "order", the result of rerunning them outer first.

**Base (left):**
- All 11 cells with inner ≥ 2 and outer ≥ 2, plus (0,6), are grey. They raise in `check_holes_share_bridge_side`:
  - with inner ≥ 4 or outer ≥ 4, the origin's iterates −3, −5 and −6 read `left` against `right` at iterate 0;
  - at (2,6) the same happens at iterates −3, −5 and −6;
  - the outer-first variant of (4,4) raises on a direct-hole pair (iterate 0 against iterates 1 and 2), which is fix A′'s case.
- (2,2) is red: the order matters, and one "skipping bridge" is logged.
- (6,0) and (8,0) are red: P3 differs from alone. Outer blasts that pass both fixed points ate the period-3 tangle.

**Branch (right):**
- Every cell matches alone, and the three order checks say "same".
- Pale yellow cells are expected:
  - outer = 0: class `e` stays unresolved, as documented;
  - inner = 6: the period-3 issue described under figure 2.
- No cell logs "skipping bridge".

Data: `grid_base.json` and `grid_branch.json` hold, per case, the words, the letters, reliability, hole/class/bridge counts, blast sizes, warnings and any exception.

### 2. `p3_blast_sweep.png`: period 3 alone, 0 to 8 blasts (the open issue)

**What it proves.**
- From 0 to 5 blasts, period 3 gives the same reliable words on both trees.
- From 6 blasts, the branch runs where base raises, but the result is unreliable and does not stay fixed as blasts are added.

**What to look at.**

Left panel, counts per blast:
- Holes, classes and active classes are 12, 9 and 3 at every count from 0 to 5, and base (dotted) and branch (solid) lie on top of each other.
- At 6 blasts the branch jumps to 42 holes, 23 classes and 22 active, then rises to 44/27/26 by 8 blasts.
- At 6, 7 and 8 blasts, base has only a "base raises" mark.

Right table (branch only):
- 0 to 5 blasts: "same as 0 blasts", reliable, with words `a->b, b->c, c->a u^-1 w^-1`.
- 6, 7 and 8 blasts: "DIFFERENT words", not reliable. The word for `c` is `a w^-1 d^-1`, then `a y^-1 d^-1`, then `a uu^-1 d^-1`, followed by a long chain `d->e->f->…`.

This is the criterion failure the outcome described: past 5 blasts, more blasts do not just refine the answer.

### 3. `hole_sides.png`: why fix A is right (period 3, 6 blasts)

**What it proves.**
- The old side reading measured each of these holes against a vertex on a different fold of its containing bridge.
- The new reading uses the hole's own image sub-arc, which puts the hole on the same side as iterate 0.

**How it was made.** The script records every propagated punch: the bridge, the carried point and the sub-arc span. It does this by wrapping `StablePartition._punch_in_bridge` for one run. It then reads each hole's side both ways: with `_bridge_side_of(..., carried)` (old) and with `..., span)` (new).

**Result.** Exactly 3 of the 40 propagated holes read differently. All three belong to origin `(41, 42)`, at iterates −3, −5 and −6. That is the same origin, and the same iterates, that the base `AssertionError` names for this run.

**What to look at, in each panel:**
- the whole bridge in grey and its image sub-arc in black;
- the carried reference point as a purple star;
- the vertex each reading measures against, with its unstable cdist: red for old, green for new.

In every panel the old vertex lies outside the sub-arc's span, on the neighbouring fold:

| Iterate | Sub-arc span | Old vertex | New vertex |
|---|---|---|---|
| −3 | 707.7–745.0 | 887 | 720 |
| −5 | 294.1–309.6 | 369 | 299 |
| −6 | 189.6–199.6 | 238 | 193 |

In each case the old reading says `left` and the new one says `right`, which agrees with iterate 0. In panel −3 the two folds nearly coincide, so the red vertex looks as if it is on the arc; its cdist (887) shows it is not. The script also checks `check_holes_share_bridge_side` on the branch's holes.

There is no matching figure for fix A′ (direct holes). Its only evidence here is the outer-first (4,4) crash at base in figure 1, and the base `p3_7` / `p3_8` errors in `grid_base.json` (origin `(42,43)`, iterate 0 `left` against iterates 1 and 2 `right`).

### 4. Six symbolic-itinerary panels per case (`panels/<tree>/<case>/`)

Each folder holds `trellis`, `partitions`, `cartoon`, `itineraries`, `transitions` and `transitions_refined`, as `.png`, `.svg` and `.fig.pickle`, drawn by the existing `scripts/symbolic_figures.draw_panels`. Every run has its own folder, so `build_nested`'s fixed slug does not overwrite anything. Titles are prefixed `[base]` or `[branch]`.

**Branch:**
- `p3_b4`: reliable, words `a->b->c->a u^-1 w^-1`, blast sizes `[16,2,2,2]`.
- `p3_b6`: runs but is not reliable (`t` unresolved), blast sizes `[16,2,2,2,2,4]`. Its `cartoon.png` makes the open issue visible: about 25 elements per branch, compared with 5 at 4 blasts.
- `nested_o2_i0`: the defaults. Reliable, blast sizes `[3,3]` (base had `[32,7]`).
- `nested_o2_i4`: reliable. Base raises here.
- `nested_o4_i4` and `nested_o4_i4_outer_first`: blasted in both orders.
  - Both reliable, with the same tables.
  - Period 3 (A) reads `a->b, b->c, c->a u^-1 w^-1`, exactly the `p3_b4` words.
  - Period 1 (B) reads `d->d uu^-1 e^-1, e->f, f->d uu^-1 d^-1`.
  - Compare the two `itineraries.png`.
- `nested_o6_i0`: reliable, and period 3 still reads `a->b->c->a u^-1 w^-1`.

**Base, for the before/after:**
- `nested_o6_i0`: runs at base, but period 3 has collapsed.
  - `a` and `b` are unresolved, `c -> b^-1`, and every other period-3 class is an inert loop.
  - The table has only 8 period-3 classes, against 9 on the branch.
  - It is not reliable.
  - Compare `panels/base/nested_o6_i0/itineraries.png` with `panels/branch/nested_o6_i0/itineraries.png`. This is the old outer blast iterating the period-3 bridges.
- `nested_o2_i0`: the defaults at base. Same words as the branch; blast sizes `[32,7]`.
- `nested_o2_i4` and `p3_b6`: no panels, only `ERROR.txt`, which holds the base traceback: `check_holes_share_bridge_side`, origin `(51,52)` and `(41,42)` respectively, iterates −3, −5 and −6.

## Scripts

All three are on the branch, under `scripts/refactor_proof/hole-side-and-own-blast/`, and are run from the branch worktree.

- **`words_grid.py`** produces figures 1 and 2. It works in two passes:
  1. Collect the data in one tree, which writes `grid_<label>.json`:
     - base: `--tree ../tangle-pack-hole-side-and-own-blast-base --label base --blast-both`;
     - branch: `--label branch`.

     `--blast-both` makes the outer-first variant blast the inner zone the way the base recipe did, with every fixed point.
  2. Draw: `--plot`.

  The references are period 3 alone at `min_separation=1e-4` and period 1 alone, using the auditor's probe recipe: 11 steps, the zone at the default pip, blasting only fp1.
- **`hole_sides.py`** produces figure 3, on the branch only.
- **`panels.py`** produces figure 4. It takes `--tree`, `--label` and `--cases`.

## Notes for the re-audit

- The branch's grid confirms the auditor's matrix. It also adds outer = 6 and inner = 6, which the auditor had not run. There are no new failures.
- **A difference from the outcome:** the auditor found the nested case reliable for outer ≥ 2 and inner ≤ 4. Inner = 6 also matches period 3 alone, but is not reliable. This is consistent with the outcome, since period 3 alone is unreliable from 6 blasts.
- A finding the outcome did not state: **at base, outer blasts alone (6 or 8, inner 0) already corrupted the period-3 words without raising.** Fix B is what removes that.
