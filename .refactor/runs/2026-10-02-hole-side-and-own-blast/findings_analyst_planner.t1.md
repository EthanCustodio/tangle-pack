# Findings: blasting nested tangles through to symbolic dynamics

*Analyst, 2026-10-02, edge `analyst_planner` iteration 1. Base commit `2315204` (branch `main`).*

**Working tree is dirty.** There are modified and deleted `.claude/` files, plus untracked `artifacts/`, `.claude/agents/analyst.md`, `.claude/athanor/state/` and figure folders (`figures/henon_symbolic_itineraries/p3_blast4/`, `p3_blast5/`, a PDF). None of them is under `src/` or `tests/`. The change branches from `HEAD`. I edited nothing tracked. I wrote only `.refactor/knowledge/{loom,topology,testing}.md`, which are new.

## What kind of request this is
This is an integration question with a bug hunt in it, followed by a refactor lens. The target is the path `TangleSession.blast_zone` → `loom/Blast.py` → `Trellis` (pips, pseudoneighbors, `StablePartition.punch_holes` / `propagate_reference_holes`) → `MinimalTrellis` → `PartitionFamily` → `DualGraph` → `DualWalk` → `SymbolicDynamics`. It was driven through `examples/henon_cases.py` on Hénon `(2, 1)`.

## Short answer
1. **No, not once blasting gets deep.** One bug in `StablePartition` breaks the pipeline: a hole's side of its bridge is read from the nearest vertex of the *whole* bridge polyline.
   - Period 3 alone survives 0 to 5 blasts, with the same words at every count. At **6 blasts it raises `AssertionError` in `check_holes_share_bridge_side`**.
   - **This is the same bug as the open "16 unstable steps" failure.** At 15 steps and at 13 steps + 6 blasts the failing orbit has identical cdist spans, iterates −3/−5/−6 fail in both, and both have the same cause.
   - At 16 steps and at 8 blasts a second form of the same mistake also fires, on direct holes at the forward iterates +1/+2.
   - Period 1 alone survives 0 to 8 blasts and is stable from blast 2 on.
2. **Nested: the bug fires in every cell that reaches ≥ 6 period-3 blasts, counting blasts received indirectly.** `build_nested` blasts both zones with `fixed_point=[fp1, fp3]`, so an **outer blast also blasts every period-3 bridge inside the outer polygon**, including period 3's *exterior*. That is not harmless:
   - the result depends on blast order;
   - inner-then-outer hits a numerics merge failure, which the blast tolerates and logs as "skipped";
   - ≥ 6 outer blasts collapse the period-3 pseudoneighbors (34 → 2), and the period-3 classes become unreachable.

   With **each zone blasting only its own fixed point** plus the side fix, **every nested cell equals each tangle alone** (16 of 16 cells, both tangles) and is order-independent.
3. **Smallest change** (details below):
   - (A) measure a propagated hole's side against its *image sub-arc*;
   - (A′) give a direct hole the side its lobe has by construction;
   - (B) one line in `build_nested`: `fixed_point=[zone.fixed_point]`.

   With all three, period 3 at ≥ 6 blasts **runs** but is **not reliable** (one unreachable class, one virtual class). That is a separate, newly exposed issue. It is out of scope here.

## Evidence
Probe scripts are in the session scratchpad (`/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/71e21633-ea51-40b5-81a9-56015ac501ee/scratchpad/`):

| Script | What it does |
|---|---|
| `probe.py` | Builds a case, times each stage, captures every WARNING, dumps words/refinements/matrix as JSON |
| `fixpatch.py` | Monkeypatches the candidate fix A + A′ |
| `ownfp.py` | Pytest plugin for fix B |
| `sides.py`, `direct.py`, `combi.py` | Hole-side diagnosis |
| `compare.py` | Nested vs alone, compared on normalised class labels |

Outputs are in `out/`. The repo was never written to.

### Matrix (unpatched HEAD unless marked)
| Case | Result at HEAD | First failure / note |
|---|---|---|
| p3 alone b0, b1, b2, b4, b5 | pass, reliable | words `a→b, b→c, c→a u⁻¹ w⁻¹`, identical at every count; 0 image bridges. **Verified.** |
| p3 alone b6, b7, b8 | **AssertionError** | `Holes of origin (41, 42) disagree on bridge_side: iterate 0 is right, iterate -3 is left; … -5; … -6`; from b7 on, also `(42, 43)` at iterates 1/2 |
| p3 alone 14 steps | pass | |
| p3 alone 15 steps | **AssertionError** | `origin (12, 33)`, iterates −3/−5/−6. Same spans as b6 → **same bug, verified** |
| p3 alone 16 steps | **AssertionError** | (8, 44) at −3/−5/−6, and (8, 52) at +1/+2 |
| p1 alone b0 | runs, not reliable | class `b` singleton, no chain. Expected |
| p1 alone b1 | reliable | `a → a u⁻¹ v⁻¹` |
| p1 alone b2, b4, b8 | reliable, identical | `a → a u⁻¹ b⁻¹, b → c, c → a u⁻¹ a⁻¹`; `a = {a_1, a_2}`; matrix unchanged. **Verified.** `min_separation` drops 1/6/23 children at b6–b8 |
| nested o0 × i{0,2,4} | runs, not reliable | outer class `e` unresolved (documented, expected) |
| nested o0 i8, o2 i4, o2 i8, o4 i{2,4,8}, o8 i{2,4,8}, o4 i4 outer-first | **AssertionError** | same hole-side bug. o2 i4 fails although p3 alone passes at b4: the outer blasts added period-3 blasts |
| nested o2 i0, o2 i2, o4 i0, o8 i0 | pass | |

After the side fix, still recipe-as-is (both fixed points):

- o2 i0, o2 i2 (both orders) and o4 i0: reliable, words = alone.
- o5 i0 reliable. **o6, o7, o8 i0: p3 classes `a`, `b` unreachable.**
  - The p3 trellis has 2 holes and 2–4 pseudoneighbors where it had 12 and 34.
  - The partition reads `[1,1,1,1,3,3]`.
  - There are 12 "loop bridge … has no registered non-loop preimage" warnings.

  There are no heteroclinic crossings: I counted the registry per blast and found only (3,3) and (1,1) pairs. The collapse is the documented full-orbit pseudoneighbor disqualification, set off by period-3 bridges that the outer blast iterated *outside* the period-3 zone.
- Order dependence: o2 i2 vs reversed gives 32 vs 33 p3 bridges, and the reversed run has two extra inert loop classes. o4 i4 vs reversed gives 50 vs 47 bridges and different words. **Verified.**
- Every inner-first both-fixed-point run that has outer blasts (18 runs, 1–4 each) logs `blast: skipping bridge that failed to iterate (AttributeError: 'NoneType' object has no attribute 'cdist')`. With `strict=True` the traceback is `ManifoldMachine.iterate_manifold:295 → merge_manifolds:438`, on p3 bridge (15, 7), an exterior bridge that the outer blast re-visits.

After the side fix **and** own-fixed-point blasts:

- **All 16 cells: nested inner words = p3 alone at the same inner count, and nested outer words = p1 alone at the same outer count** (normalised labels + directed tokens). **Verified.**
- Reliable wherever outer ≥ 2 and inner ≤ 4.
- o4 i4 and o8 i4 run in either order with identical words, refinements and bridge counts. **Verified.**
- Zero Blast "skipped" warnings.
- With B alone (unpatched), o8 i0/i2/i4 also pass. So with B the nested case only fails where p3 alone fails.

### Root cause of the assertion (verified)
`_punch_in_bridge` sets `bridge_side = _bridge_side_of(trellis, bridge, carried)` (`StablePartition.py:1526-1527`). That calls `_arc_side_of` (`:719-727`): the nearest vertex of the **whole** containing-bridge polyline, with a tangent from its two neighbours. For origin (12, 33) at iterate −3:

- the image sub-arc spans unstable cdist 707.7–745.0 inside bridge (34, 15), which spans 551.6–1076;
- the carried point's nearest vertex is at cdist **887**, on another fold, 0.0039 away → "left";
- against the sub-arc it is "right", like iterate 0 and every other iterate.

Exactly 3 holes change (−3, −5, −6) and none changes in any run that passes today.

Direct holes (`StablePartition.py:349`) have a second form of the same mistake.
- The +1/+2 pairs of the period-k exemption sit on tiny 17-node lobes.
- The hole sits beside the bridge end. The nearest vertex is index 1, whose neighbour `poly[0]` lies past the crossing. The answer comes out "right", while the lobe polygon (bridge + stable arc) says "left", the same side as iterate 0.
- Direct-hole openings are inward by construction (`inward=True`), so only I1 is tripped; the partition is unaffected.

**The invariant was right and the measurement was wrong.** The +1..+(k−1) exemption only exposes it, so it stays as the constraints require.

The full suite with A + A′ monkeypatched gives 782 passed / 1 skipped, the same as HEAD. With A + A′ + B, `test_higher_period_cartoon.py`, `test_loom_blast_restore.py` and `test_stable_partition_period3.py` give 24/24 pass. **Verified.**

### Warnings, and which point to bugs
| Warning | Where | Verdict |
|---|---|---|
| `inert class u has the non-empty word 'v'` (4× per p3 run) | every p3 run | Expected for period k (inert → next branch's inert class). Noisy, not a bug |
| class `e` / `b` unresolved, singleton | nested o0, p1 b0 | Expected ("grow or blast") |
| `Skipped N pair(s) with no spanning bridge` | deep blasts | Documented blast behaviour |
| `blast: skipping bridge … NoneType … cdist` | inner-first + outer, both fixed points | **Bug** (numerics merge). Avoided by B |
| loop bridge with no preimage, DualGraph component with no closed face, walk unreachable (a, b) | outer ≥ 6, both fixed points | Consequence of blasting p3's exterior. Avoided by B |
| `image bridge … abuts no hole …`, class `t` unreachable, virtual `new1`, `is_reliable=False` | p3 ≥ 6 blasts / ≥ 15 steps after fix A | **Open.** A new reference orbit appears (20 backward holes → chain of ~18 active classes d→…→t). 12 blasts does not resolve it. Out of scope |
| straddling landing (`contained=False`), evidence disagreement | only o4 i4 outer-first, both fixed points | Consequence of the both-fixed-point recipe |
| same-stability crossings | none in any run | — |
| `unmatched_members` | 0 in every run | — |

`image_of` was run on all regions of six cells. It never asserts; it returns None for many regions, by design. This is evidence only, not a guard on this path.

### Cost
- 0.1–5.3 s per cell. Blasting is ≥ 80 % (p1 8 blasts: 4.1 s; nested 8+8: 5.1 s, 403 crossings).
- All topology stages together take < 0.2 s.
- cProfile of the blast: `Tangle.add_manifold`/`_insert_segment` + rtree queries ~50 %, `refine_manifold` ~30 %.
- Intermediate `_repin`s are dead work: identical classes and words without them (verified, nested 4+4).

## Recommendation (smallest change)
**A. `topology/StablePartition.py` `_punch_in_bridge` (`:1526`).** Classify `carried` against the polyline restricted to the image `span` (nodes with cdist in span, plus one node either side). Then derive openings from that side as now. The cleanest form is an optional `span` argument on `_bridge_side_of`, with `_punch_in_bridge` the one caller that passes it. This changes 3 holes in the failing cases and none elsewhere.

**A′. Same file, `:349`, direct hole side by construction.** Stop measuring and use the lobe's side. Pick one of two forms; both agreed on 22/22 direct holes (p3 b0/b4/b8, p1 b0/b2/b6, k=10):
- **Combinatorial**, in the spirit of `row_of_end`: with `first` the lower-unstable-cdist end, `side = left iff crossing_sign(first) · (+1 if second is outward of first on the stable branch else −1) > 0`.
- **Geometric**: the sign of the lobe polygon area.

Either way I1 keeps its meaning for propagated holes, which are still measured.

**B. `examples/henon_cases.py:253-258`.** Change the argument to `fixed_point=[zone.fixed_point]`. This is a behaviour change to the author's recipe.
- At the defaults (outer 2, inner 0) the words are unchanged, but p3 bridges go 31 → 27.
- Reword the "inner first so the outer blasts do not consume its bridges" docstring (`:229-230`).

**Tests:** fast regression tests, each < 2 s.
- `build_period3(blasts=6)` and `build_period3(unstable_steps=15)` reach `symbolic_dynamics()` without raising (assert it runs, not `is_reliable`).
- Nested 4+4 gives the same words in both blast orders.
- The inner classes of `build_nested(outer_blasts=2, inner_blasts=4)` equal `build_period3(blasts=4)` on normalised labels.

## Out of scope
1. **p3 ≥ 6 blasts after A is unreliable.** A new reference orbit's 20 backward holes give a chain of ~18 active classes, an unreachable `t`, a virtual `new1`, and "crossing abuts no hole" in `PartitionFamily`. It needs the author's view.
2. `numerics/ManifoldMachine.py:285-295, 438`: a bridge whose end nodes were iterated by both neighbours builds a non-contiguous `old_iterated_points` and crashes in `merge_manifolds`. B avoids it; the bug remains. Mechanism believed, traceback verified.
3. `_hole_openings` (`StablePartition.py:1211-1217`) uses the same whole-polyline nearest vertex. Latent and believed; no failure seen.
4. The session forgets the strong-pip choice on every generation bump, hence `_repin` after each blast. A pip kept as an `Intersection` and re-resolved via `registry.find` (as zones already do) would fix that.
5. `Blast.py:226-236` duplicates `_resolve_fixed_points`. `:289` re-imports cKDTree in the loop, and `_min_interior_distance` keeps an unreachable scipy fallback. `BlastResult.fixed_point` is typed `Optional[FixedPoint]` but holds lists.
6. `build_nested` slug is constant (`henon_cases.py:274`), so figure runs overwrite each other.
7. The "inert class has a non-empty word" WARNING fires on every p3 run; consider INFO for inert → inert words.
8. Blast speed: per-segment Python bookkeeping in `Tangle._insert_segment` and the rtree query loop dominate (numerics layer).

## Claims marked
- **Verified** (probe/pytest):
  - all matrix cells;
  - the identical 15-step and 6-blast failure signature;
  - nearest-vertex cdist 887 vs sub-arc 707–745;
  - direct-hole vertex-1 vs lobe side;
  - the combinatorial rule = lobe area on 22/22;
  - nested = alone 16/16 and order independence with B;
  - the both-fixed-point order dependence and collapse;
  - suite 782/1 with A+A′;
  - 24/24 higher-period tests with A+A′+B;
  - `_repin` irrelevance;
  - the timings and profile.
- **Believed:**
  - the `merge_manifolds` crash mechanism;
  - the latent `_hole_openings` misfire;
  - the combinatorial direct-hole rule beyond the 22 holes tested.
