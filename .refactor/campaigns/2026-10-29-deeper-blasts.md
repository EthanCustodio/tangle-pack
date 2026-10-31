# Campaign: deeper blasts

**Outcome (2026-10-31, complete in 4 runs, main 7015483).** Run 1 brought the lean blast and the opt-in width guard: lean alone frees about 1% because exact topology keeps the frontier's lineage, lean plus width ends 3-5x smaller, and the width guard stays off by default. Run 2 made the near-anchor guard (clause 3 of N, c=5) the default, so `pip_depth` equals `trusted_depth`, and added the opt-in graded mesh `near_orbit_radius`, which reaches depth 18 at 1e-7 and 1e-9 (the xy floor). Run 3 added the cdist residual column (double-double): the 773 plain neighbour ties of the deep unguarded series are all split, nothing computed changed, the merge stays plain, and it costs +17 MiB and about 15% deep. Run 4 added representative blasts: the topology matches the full series except a delayed k32 reference switch, but they are only 2-2.6x cheaper while the pip depth is stalled. The prune question is answered: at trusted depth 9-12 (1e-9 to 1e-12), 1e-5, "resolution" and None agree on everything. Open for a later effort: holding unclassed bridges past the cut, or representatives with the graded mesh; a precise merge; the heteroclinic `blast_zone(fixed_point=None)` TypeError; the forked-pool futex hang in the suite.

- **Goal:** blast deeper, to more iterates and more pip depth, with less memory (author, 2026-10-29).
- **Date:** 2026-10-29. **Base:** main `5dec003` (planned on `5ed9940`; the speedup campaign's run 1, cpu-workers,
  merged in between and touched `loom/Blast.py`'s separation, `diagnostics.py` and `DynamicalSystem.py`; its
  separation time is unchanged). Each run branches from main as it stands when the run starts.
- **Cap:** 4 runs, all used.
- **Author decisions at planning:**
  - Item 1 (lean r-tree) is dropped. The scout showed the index is already lean: the rtrees hold 0.3 MiB at
    the 1e-7 plateau. The 110 MiB "segment index" is the `Tangle._query_edges` key set, and it stays.
  - Item 3 (lean blast) is IN, and first: "much more important" (author, 2026-10-29). This reverses the
    2026-10-28 rejection of a lean blast. To fit the cap it is merged with item 2's width guard. Lean mode
    cannot run the geometric guard, which needs every interior point, but the width guard needs only one
    number per bridge.
  - Added 2026-10-29 at the author's request: run 3 cdist-precision, before class-representatives (run 4).
- **Scout** (the k=2.5 scale series, `profile_numerics.scale_build`, 4 pinned cores):
  - **Separation** (`Blast._min_interior_distance` plus the per-generation `cKDTree` over every interior point,
    about 2M points) takes 4.8 of 32.5 s per series at 1e-7 and 8.4 of 45 s at 1e-9.
    - Its share of the peak RSS (1172 MiB peak against a 232 MiB census) is unmeasured. The estimate is about
      100 MiB.
    - D1 found that diagnose's whole-field arrays dominate the transient.
  - **Depth stalls at m = 9 (1e-7) and m = 12 (1e-9).**
    - f^m(q0) walks into P by about β ≈ 5.4 per level, and f^10(q0) is 5.2e-7 from P.
    - The f^9 bridge (248, 249) is iterated, but its image's nearest node is 6.5e-3 from P. f^10(q0) is never
      detected.
    - `ManifoldMachine._refine_layer` uses an absolute `area >= area_cutoff` and an absolute `apart >= 1e-8`
      floor. So features below about √area_cutoff near P are invisible.
    - The densify run (2026-10-21) found that the zone filter (one midpoint) stops depth 1-4 levels later.
  - **The cdist precision ceiling.**
    - Unstable cdists reach 6e12 at blast 14, and the smallest relative step dc/c is 4e-12 at 1e-7 (1e-12 at
      1e-9).
    - Per blast (1e-7, 1e-5 guard): 5.4e-7 at blast 4, 5.3e-10 at 8, 3.3e-11 at 12, 4.3e-12 at 14, then flat
      because the guard stops new curve. Before the plateau it shrinks about β per blast, a decade every
      1-1.5 blasts.
    - About 4 decades remain to float64's 2.2e-16, roughly 5-7 blasts once the guard is gone.
    - `DEDUP_RTOL` (1e-9 relative) has been beaten since blast 8: crossings on one branch are told apart only
      by their stable cdist.
    - The coordinate limit behind it: xy near P resolves to about 6e-16, and f^m(q0) reaches it at
      m ≈ 21-22. That needs local coordinates around P, outside this campaign.
  - **Item 4's lever.**
    - Without separation the series grows about 1.68× per blast: 839 bridges and 1.1M points at blast 11.
    - The alphabet is 15 classes over 23 minimal-trellis bridges at blasts 8-12.
- **Error analysis (2026-10-29, after planning).** Source: `figures/henon_k25_entropy_post_array_refactor/errors/`
  (`errors.png` and its README), from `scripts/henon_k25_errors.py`. It reruns every pruned k=2.5 sweep, p1 / p3 /
  nested at 1e-7 to 1e-14 with `min_separation` 1e-5, until the frontier is exhausted, and records the numerical
  errors at every blast. Every sweep reproduces the original pip depth exactly.
  - **The deepest orbit points are unresolved.** The registered orbit contracts toward P by β per step (5.56 for
    p1; λ per period for p3) only until |f^j(q0) − P| reaches about the image curves' chord deviation (1e-9 to
    1e-10 at fine cutoffs).
    - Past that, the crossing positions and their stable order are noise. The resolution ratio (about 1e4) does
      not catch it, because it never compares a crossing's distance to the anchor, or to its stable neighbour,
      against h.
    - Trustworthy depth against reported depth: p1 7–13 against 9–22; p3 16–33 against 17–40.
    - Consequences: p1's 1e-13 rise to 0.5281272 is an artifact (trustworthy h is 0.5280489), and so are p3's
      values above 0.1877778. Nested 0.5281517 stands: it is reached at a trustworthy depth.
  - **Two stall mechanisms, by cutoff.**
    - **1e-7, and p3 at 1e-8:** the image curve reaches f^{m+1}(q0), but the crossing's ratio is 0.7–2.1. The
      unstable chord deviation (h ≈ 1.5·A/L on long, nearly straight segments) exceeds the short stable segment
      near P. This is the case the scale-relative cutoff addresses.
    - **1e-9 and finer:** the bridge ending at f^m(q0) is never iterated. `TangleSession.classify_bridge` tests
      its index-middle node (`polyline_midpoint`) against the frozen zone polygon. The bridge is a long arc
      (length 15.6 to 53.5) whose middle lies outside the zone, so its zone is None. Confirmed on p1 at 1e-10,
      1e-12 and 1e-13, and on p3 at 1e-12. At 1e-7 the same 15.6-long bridge happened to classify inside.
    - So in the fine regime, "a finer cutoff gives more depth" is mostly a change in which bridge ends at f^m,
      not resolution.
  - **Stable bracket near P.** Every deep orbit point lies on the FIRST stable segment from the anchor (about
    2e-6 long for p1). The stable side of `infer_iterates` therefore cannot tell those points apart, and links
    there rest on the unstable cdist alone.
  - **cdist precision.** The smallest relative step dc/c reaches float64 eps at p1 1e-12 and 1e-13, blasts 13–20.
  - **Same-stability crossings.** p1 at 1e-13 has 37k–420k per blast at blasts 14–21 (the deep, exempt
    pip-orbit children), from near-coincident unstable curves near P.
  - **The prune.** Unpruned (no `min_separation`) against pruned gives more classes and a higher h at p1 1e-11
    (depth 14) and p3 1e-12 (depth 33). Both lie beyond the trustworthy depth (12 and 29), so there is no evidence
    yet that the prune hides real topology. At 1e-9 the two are identical to blast 17–45.
- **Overlap:** run 1 of `2026-10-28-speedup-options.md` includes the separation cKDTree. When run 1 here
  merges, note there that the width guard leaves that piece for the geometric opt-in only. Lean mode also
  changes what that campaign's sweep runner can save and resume (no save/load in lean mode).
- **With the speedup campaign:** run its runs and these one at a time, never side by side. They share
  `diagnostics.py` and `Blast.py`, the ff-only merge would force a rebase, and both measure blast time, which
  needs an idle machine.
- **Every run:**
  - `env/bin/python -m pytest -q -rxX` (only KNOWN_ISSUES xfail) and `coverage_guard` against the base;
  - `stress/profile/identity_probe.py` unchanged on the default path;
  - `profile_numerics.py --cases scale` before and after, on an idle machine or pinned cores.

## Run 1: lean-blast (/feature)
- [x] done (876f241, lean keeps the frontier lineage, so alone it frees <=1.3%; the opt-in width guard reaches the 1e-5 guard's pip depth at 1/19 (1e-7) and 1/4.5 (1e-9) of the rows, lean+width a further 3.5-4.8x; width guard off by default, WIDTH_SEPARATION_FACTOR=1.0 uncalibrated, author's choice)

Goal: an opt-in lean blast (`lean=True`, default off) holds only what the topology reads, plus the frontier.
Its separation is the analytic width guard (`min_separation="width"`), which is also available to ordinary
blasts.

**Measure first.** A blast's peak RSS and census, split into point store, frontier, separation cloud and tree,
and diagnose's whole-field arrays, at 1e-7 and 1e-9, on the scale series and an unguarded series. This sets the
target and shows what lean mode removes.

**Lean mode.**
- **Kept after each iterate:**
  - crossings, cdists, brackets and iterate links;
  - bridge ends, as `End` values;
  - full polylines for the minimal-trellis bridges and the pip-orbit bridges;
  - the frontier's polylines.
- **Released:** every other bridge interior becomes ends only, its store rows freed. A bracket's two rows are
  kept, because the registry pins them.
- **Chunked frontier:** map a chunk, intersect it (query-only against the stable index), register, infer
  iterates, release interiors, then the next chunk. The peak is bounded by the chunk, not the generation.
- **Forced off, each with a clear error:**
  - the geometric `min_separation` (the width guard is allowed);
  - diagnose;
  - the dense arrangement;
  - `is_computed_arc` re-indexing of released curves;
  - `save` / `load`;
  - plotting of released curves.
- **Store row release** (`PointStore`): a compaction or a free list, chosen by the plan. Chains and positions
  stay O(1). A released bridge reads its ends and refuses `get_point_array`.
- **Minimal trellis between steps.** Needs a topology recompute per step, or a cheaper rule ("hole bridges and
  their image chains") the plan defines. The plan decides.

**Width guard.**
- A child's width is carried down its image chain: `w_child = w_parent / β_s` per step (`per_step_beta`). The
  seed is lobe area / unstable length (bridge plus stable segment ring).
- A child is dropped when `w < c · h`, with h the chord deviation at `area_cutoff`
  (`geometry.chord_deviation_batch`) and c chosen in the plan.
- The pip-orbit exemption stays.
- A drop is remembered through `Bridge.dropped_at`, so blasts compose.
- The geometric guard stays as an option for non-lean blasts.
- Near-misses between different lobes are left to the same-stability diagnose (off in lean mode; known).

Files: `loom/Blast.py`, `loom/TangleSession.py`, `numerics/PointStore.py` (row release),
`numerics/TangleWorkbench.py`, `numerics/BaseManifold.py` / `numerics/Bridge.py` (released state),
`numerics/persist.py` (refuse lean), `topology/MinimalTrellis.py` if read between steps, tests, CLAUDE.md rows.

Success:
- **Topology equal.** Lean against full blast, with the same guard, on the scale series and `cases.SERIES`, to
  the full series' depth: crossings of N, iterate table, alphabet identities, words, refined matrix and
  `pip_depth`.
- **Memory and time.** Census and peak RSS lean against full, at 1e-7 and 1e-9, with and without the width
  guard. Blast time.
- **The width guard** against `min_separation=1e-5`: depth, bridges, points, separation time ≈ 0.
- **Defaults unchanged.** The identity probe is identical on the default (non-lean) path.
- Suite and `coverage_guard`.

## Run 2: anchor-depth (/feature)
- [x] done (5990bb9, near-anchor guard on at c=5 (orbit-distance half only, floored at the xy resolution; the stable-pair rule stays out, author) cuts pip_depth to trusted_depth; pip-end arriving-bridge zone rule; trusted_depth as data; opt-in graded mesh r0=1e-3 gives depth 18 = trusted at 1e-7 and 1e-9 (base 9/7, 12/9), at the xy floor, with no more points; k25 errors sweep not rerun)

Goal: pip depth stops being set by the absolute cutoff near the periodic orbit, and every reported level is a
resolved one.

The order of the items is set by the error analysis above.
1. **Pip-orbit bridges classified by their pip end** in the zone filter, not by their midpoint. This is the main
   fix: from 1e-9 down, every stall is this filter, so it goes first.
   - Any bridge ending on a selected pip orbit is iterated whatever its midpoint says.
   - Check it against the errors sweep: the bridge ending at f^m(q0) (p1 1e-10/1e-12/1e-13, p3 1e-12) gets
     iterated, and the depth moves.
2. **A near-anchor resolution guard** (new; it must land WITH item 1, or item 1 deepens into noise).
   - A crossing whose distance to its anchor, or to its stable-order neighbour on the same stable segment, is below
     `c · max(h_u, h_s)` is not a crossing of N. It joins clause 2's family, is reported with the remedy
     "refine", and stops the pip orbit like any unresolved f^{m+1}(q0). The plan chooses c (3, like the ratio, is
     the default).
   - Validate against the errors sweep: on every record the guard should cut the orbit at the trustworthy depth
     (the β-contraction test), not later.
3. **A scale-relative refinement cutoff** near the orbit, `area_cutoff · min(1, (r / r0)²)`. The 1e-8 floor
   scales the same way. The plan decides r0, and whether this is opt-in or the default.
   - It covers the coarse-cutoff ratio stall (1e-7, p3 1e-8).
   - Its exponent must make the chord deviation fall FASTER than r, so the guard of item 2 does not stop the orbit
     right away. Check h/r along the orbit with the errors figure (panels c, f).
   - It covers the stable branch too: the first stable segment from the anchor (about 2e-6 for p1) holds every
     deep orbit point, so the stable bracket in `infer_iterates` cannot tell them apart.
4. **The trustworthy depth as data.** `BlastResult` (or a diagnostic) reports the depth up to which the registered
   orbit contracts by β per step (λ per period), inside the linear regime. This is the test in
   `scripts/henon_k25_errors_figures.py:trustworthy_depth`. It is a flag, never a decision (the topology never
   reads the map), but every series report quotes it beside `pip_depth`.

Files: `numerics/ManifoldMachine.py`, `numerics/TangleWorkbench.py`, `loom/Blast.py`, `loom/TangleSession.py`,
tests, design note §A.1 / §A.10.

Success:
- `pip_depth` past 9 (1e-7) and 12 (1e-9) on the scale series, and on the five B3 series, AND equal to the
  trustworthy depth: no orbit point inside the guard;
- the k=2.5 errors sweep rerun (`scripts/henon_k25_errors.py`): trustworthy depth per cutoff before and after,
  and the entropy at that depth. Before: p1 7–13 with h 0.5280489; p3 16–33 with h 0.1877778.
- points per extra level;
- the dc/c floor at the new depth;
- the same-stability count at the deep blasts (37k–420k per blast at p1 1e-13 before);
- the laws pass (topology change is discovery).

## Run 3: cdist-precision (/refactor)
- [x] done (1bff84d, residual column cdist_lo (double-double); probe identical on 62 builds; deep unguarded 1e-7 series: 773 plain neighbour ties at blast 15 all split by the residuals, same rows/links/depth 16 as base; merge, crossing sign, BRACKET_PAD 1e-12, ownership and dual walks stay plain (a precise merge changes k28/inversion: a /feature if wanted); +17 MiB at the 1e-7 plateau, deep blasts ~15% slower)

Goal: canonical distances keep their full relative precision however deep a curve is. Nothing computed changes
on today's cases.

Confirmed by the error analysis: dc/c already reaches float64 eps at p1 1e-12 and 1e-13 (blasts 13–20). Run 2
takes depth past that, so this run is needed there. It fixes cdists only. Positions near P (|f^m − P| ≈ 6e-16 at
p1 m ≈ 21–22) stay a float64 limit, but the trustworthy depth (13 for p1) is far from it, so local coordinates
around P stay deferred.

**The fix in plain terms.** Today every point stores its cdist as one big number, like a street address written
as "4,000,000,000,017 metres from the start of the road". Two neighbouring houses then differ only in the last
digit, and float64 runs out of digits. The fix is to write it like a real address: each chain (one stretch of
curve) stores where it starts, once and precisely, and each point stores only how far it sits from that start,
"17 metres past the corner". Two points on the same chain compare by those short offsets, so they keep every
digit. Two points on different chains compare by their starts first.

- **Storage.** Per chain, a precise base cdist: a double-double, or an exact integer exponent of the stretch
  factor plus a mantissa; the plan chooses. Per row, a float64 offset from its chain's base.
- **Writes.** `PointStore.append`, the refinement midpoint, merge (`_image_chains`, `merge_manifolds`) and
  image linking, where the base is multiplied by the chain's `stretch_param`.
- **Reads.** Through ONE accessor for "compare two cdists" and "difference of two cdists", in place of raw
  `store.cdist` reads:
  - the crossing order and dedup (`DEDUP_RTOL` applied to the offset within a chain);
  - the brackets;
  - `infer_iterates`;
  - `image_cdist` / `SCALING_RTOL`;
  - the pseudoneighbor intervals.
- A plain float `cdist` view stays for plotting and for any reader that needs no precision.
- `persist` format bumped (format 2), and format 1 still loads.

Files: `numerics/PointStore.py`, `numerics/ManifoldMachine.py`, `numerics/BaseManifold.py`,
`numerics/IntersectionRegistry.py`, `numerics/Tangle.py`, the iterate-inference module,
`numerics/persist.py`, the cdist readers in `topology/` (`Trellis`, `Pseudoneighbor`), tests, design note §A.0
(data model).

Success:
- the identity probe identical on today's cases;
- a deep unguarded scale series (from run 2's depth, `min_separation=None`, as far as memory allows) shows no
  tie between neighbouring nodes and no crossing-order flip, against the float64 build where ties appear;
- census rows' cost reported (offset column plus a per-chain base);
- suite and `coverage_guard`.

## Run 4: class-representatives (/feature)
- [x] done (7015483, representatives= opt-in; probe identical on 62 builds; full-series topology everywhere but a delayed k32 reference switch at blasts 6-9/10 (a held witness) and 1-4 fewer trajectory pairs; 1.9-2.6x cheaper at k=2.5 blast 12, not 10x: past the stalled pip depth the frontier is unclassed; lean+rep unguarded to blast 18 at 15.3M rows; the prune hides nothing at trusted depth 9-12)

Goal: an opt-in blast that iterates one bridge per (class, landing element) plus every pip-orbit child, so cost
scales with the alphabet and not with e^{hn}.
- Classes come from the current `Topology` between steps.
- Non-representatives stay held and un-iterated, with a recorded reason, so a later full blast can take them.
- **Validation first:** against the full unpruned series on `cases.LAW_CASES` and `cases.SERIES`, to the full
  series' depth, compare alphabet identities, words, the refined matrix and `pip_depth`.
  - Theorem T is checked, not proved.
  - Every difference is reported per case.
  - The effect of a thinner X on pseudoneighbors is measured.
- **Does the prune hide topology?** (The author's question, 2026-10-29.) Compare `min_separation` 1e-5,
  `"resolution"` and None at the same TRUSTWORTHY depth (run 2's guard on), at 1e-9 to 1e-12.
  - Today's unpruned differences, at p1 1e-11 and p3 1e-12, lie beyond the trustworthy depth and so decide
    nothing.
  - Run 4's representatives are cheap enough to make the unpruned side affordable.

Files: `loom/Blast.py`, `loom/TangleSession.py`, `numerics/Bridge.py`, a regression test over the series, a
`scripts/` deep-series experiment, CLAUDE.md.

Success:
- the equivalence table is clean, or every difference is explained;
- bridges and points at least 10× below the full series at blast 12;
- run in lean mode as well (it composes with run 1): the deepest blast reached before time or the cdist
  ceiling.
