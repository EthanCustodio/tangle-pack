# tanglepack design note

*Verified against `HEAD` `ad1d96a` on 2026-10-16.*

The library's design, one section per pipeline stage, as the code at `HEAD` implements it. It is
section A of `docs/alignment_plan_2026-10-09.md` rewritten after the slices B0-B9b closed: the
planning text is gone, every rule is the one the code runs, and each section names the functions that
run it. The `src/` docstrings point here as "design note §A.x"; `CLAUDE.md` gives one paragraph per
stage.

How to read it:

- **N** is the set of RESOLVED transversal crossings (A.2). **X** is N together with every iterate
  `M^n` of its points; a point of X may lie outside the trellis (past a cut). **T(q0)** is the
  trellis of the strong pip q0 (A.4). `M` and `f` both name one step of the map.
- **plan Qn** is question n of `docs/alignment_plan_2026-10-09.md` section E.
- An anchor in parentheses, "(D13)", is a decision or issue id. Its full statement and history live
  in `.refactor/decisions.md` (D1-D15, C1-C22 and every slice outcome) and `docs/open_issues.md`
  (the issue catalog, the D/C tables and the WP outcomes); the statement is folded in here so the id
  is never needed to read the rule. Section I lists every id with its section.
- "The author's rule" marks a rule the author stated; everything else is a call made under D8
  ("simplify freely") and revisable.
- Names in backticks exist at `HEAD`. Deleted names appear only in section H.

## A.0 Global rules

**Topology reads the table, never the map.** Every iterate the topology needs is a lookup in the
iterate table (`Trellis.iterate`, which composes the `±1` links) or, when the table has no link, a
canonical-distance scaling by the branch's per-step factor (`Trellis.image_cdist`, the one scaling
path; `FixedPoint.advance_key` and `FixedPoint.per_step_beta`). Nothing in `topology/` reads an
eigenvalue directly. There are exactly three exceptions, none of them part of an algorithm:

- `Arrangement.image_of` carries a region's corners by the iterate table (a corner whose image is
  not in the trellis gives None) and VERIFIES the candidate face by its area and by mapping ONE
  interior point (`representative_point`) with the real map (an acceptance filter, `_verify_image`);
- `diagnostics.measured_sides` carries a reference hole back with the real inverse map to measure a
  propagated hole's side, as a resolution flag; it runs on the default
  `Trellis.punch_holes(diagnose=True)` path;
- `plotting.hole_point` carries the same point back to place a hole in a figure.

**Flags never decide.** A geometric reading (polyline rows, measured hole sides, same-stability
polyline crossings, resolution ratios outside the definition of N) is reported, never used to change
the topology. The one place a resolution criterion enters a definition is N itself (A.2).

**A class exists iff a computed bridge spans it or the word of an existing class names it** (the
closure from the spanned classes; plan Q9 reads D11 and D14 together). A pair no computed bridge spans is a
class all the same, flagged "no computed bridge yet" (A.8c).

**Completeness is not a theorem.** The walk constraints (landing spans, the marker law, walls) are
the facts the real image of a bridge satisfies; nobody has proved they always single out one walk. A
tie that survives them is honestly "ambiguous", and a search that cannot establish its walk is
"unresolved" (A.8c). Neither is ever broken by an arbitrary order.

**Discovery is not failure.** Blasting or growing further may change the topology (a new reference
pair, more holes, more classes); that is expected (the author's rule). A decrease is examined, never
assumed to be a bug (A.10, O).

## A.1 Stage 0 — manifolds

**Definition.** For each branch, a polyline approximation of the manifold with canonical distance `c`
non-decreasing along it (asserted after every merge). Growth maps the fundamental segment `[x, M^k x]`
one step at a time and refines so that every chord deviates from the image curve by less than
`area_cutoff` (the chord-deviation rule `h = 1.5 A / L`). One fundamental segment is seeded per chain
of `k_value` pieces and mapped around it; an inversion chain covers both branches.

`FixedPoint` is the one source of branch bookkeeping: `advance_key` (what n map steps do to a
manifold key), `per_step_beta`, `branch_cycle` and `preserves_orientation`. Orientation is a fact of
the orbit: the sign of the solver's single-step Jacobian determinant, read by the topology and never
computed there (C16). Orbit points and eigenvectors are `(2,)` arrays and eigenvalues floats, so no
reader ravels (API-6).

The image of a curve whose points already have images may lie on several curves: the existing images
are merged chain by chain, a chain breaking where the walk from the first image breaks (NUM-9).

The `grow_until_*` drivers cap their rounds at `max(10, ceil(log_stretch / ln stretch) + 1)` by
default, so a small stretch gets the depth it needs (NUM-15); `grow_until_arclength` stops on the
polyline's arc length, not on the tail's canonical distance (NUM-13). Every top-level growth call
counts the points it creates and raises `PointBudgetExceeded` before the map call that would exceed
`max_points_per_step` (a resource guard, never an invariant; S). The default `area_cutoff` is 1e-7
(S). The opt-in `workers` (default 1) splits a large map batch across a thread or forked process
pool and concatenates the chunks in order, so for an elementwise map it changes no result
(cpu-workers, 2026-10-29).

**The graded mesh near the orbit (anchor-depth, 2026-10-29, opt-in).** With `near_orbit_radius = r0`
set (before growth; None, the default, is the refinement above byte for byte), a pair within `r0` of
an orbit point of its curve's own fixed point refines when its curvature area is at least
`area_cutoff * (r / r0) ** 2` (so `h ~ r ** 2` falls faster than `r` and the near-anchor guard of
A.2 does not stop the pip orbit) or when its chord is longer than `r`, down to a spacing floor
`max(1e-8 * r / r0, xy_resolution(P))` (`xy_resolution = 64 eps max(1, |P|)`). Two choices are not
the plan's and were forced. `r` is how near the pair's CHORD passes the orbit point, not its
midpoint's distance: a deep blast image curve runs straight past `P` between nodes 0.04 apart with
its closest approach at 8e-6 (1e-7 scale series, blast 7), and the midpoint rule never refined it.
And growth refines only the images it makes, so the seed's segment at the orbit point would never
be graded, and a point refined into it would have no iterate (the merge of its image fails). So
`ManifoldMachine.grade_seed` grades the first fundamental domain of each seed chain once, at
initialisation, and carries its points inward one piece of the chain per map step, each linked to
its source as its iterate, placed on the ray from the target piece's orbit point through its
innermost row at the distance its cdist gives (the seed's own linear approximation), down to the xy
resolution. They are not mapped: the forward map along a stable seed multiplies the seed's distance
from the true manifold by `lambda_u / lambda_s` per step (measured: the fourth generation was off by
its own size). The pair at the orbit point keeps the plain rule. The point budget is unchanged.

The data model (one `PointStore`, chains of its rows, curves as views of a chain whose row ends read
as `End` records, and `BranchPoint` a view of its `FixedPoint`: the orbit point's coordinates, a cdist
of 0.0, the next orbit point as its image, and the fixed point's slot table
`(stability, orbit_index, branch) -> Chain`, on no row) is a representation choice, not part of the
definition. Rows are appended and never reused; only a lean blast removes any (`PointStore.compact`,
lean-blast 2026-10-29): it renumbers the survivors in order, and every holder of a row (curve ends,
the registry's brackets, the segment index's records and edge keys) is remapped through the
old-to-new map (`TangleWorkbench.release_curves`, A.10).

**A cdist is a product with a residual (cdist-precision, 2026-10-30).** A canonical distance is a
PRODUCT, never an arc-length sum: the seed's Euclidean distance, scaled by the chain's stretch one map
step at a time (`stretch * c`, the graded seed's `scale * c`), averaged by refinement
(`0.5 * (c0 + c1)`) and interpolated at a separator or a crossing (`(1 - t) * ca + t * cb`). Deep in a
blast it reaches 3e13 with neighbours `dc / c` below float eps, where plain floats tie. The store
keeps beside the plain float a residual column, `cdist_lo`: what each of those expressions' rounding
dropped (`numerics/precise.py`: error-free `two_sum` and fma-free `two_prod`, so `hi` is the plain
expression bit for bit and `hi + lo` carries about 32 digits; 0.0 for a seed row and the
initializer's fictitious tails, NaN with no cdist; 8 bytes a row). `Intersection` carries
`unstable_cdist_lo` / `stable_cdist_lo`, and `cdist_key(stability)` is a crossing's ORDER key, the
pair `(cdist, residual)` compared as a tuple (a scaled cdist is `(c, 0.0)`): the plain float decides,
the residual only breaks a plain tie. It is not the precise value, which the merge's chain can
contradict (below); the key never does, so a bridge cut in key order runs forward along its chain
(`Tangle.create_bridges` asserts it) and its id is in plain order too (the sparse `Arrangement`
asserts that). The key orders the registry, the bridge cut and the `_register_bridge` checks, each
branch's crossing order, the pip-orbit cuts, the hole bridge's ends, the containing bridge, the
pseudoneighbor intervals and the strong-pip box. The precise value (`hi + lo`) compares in dedup
and in the iterate inference's carried brackets (`dd_mul` / `dd_div`). `BRACKET_PAD` stays
1e-12: at 1e-24 the k=2.8 deep-blast sliver build links three crossings fewer, so the carry does
not meet every image to double-double precision. The PLAIN float still orders the merge, its run
cuts and the chain check (a precise merge keeps 17 plain ties the zipper drops on the k=2.8
one-blast build and reorders neighbours one ULP apart on the k=2.8 and inversion builds: a change
of what is computed, the author's to make), the crossing sign (the segment's chain order), the
stable partition's intervals and ownership, and the dual-graph nodes. Scaled cdists and plotting
stay float. Save format 2 holds the residuals; a format-1 save loads with 0.0. On the unguarded
k=2.5 scale series at 1e-7 (`near_orbit_radius = 1e-3`, `min_separation=None`) neighbours first tie
as plain floats at blast 10 (8 ties; 773 at blast 15, 9.1M rows); the residuals tell every one of
them apart, and none runs against the chain.

**Implemented by.** `ManifoldInitializer.construct_kevin_way`, `ManifoldMachine.iterate_manifold`,
`ManifoldMachine._image_chains`, `TangleWorkbench.grow_n_times`, `TangleWorkbench.grow_until_bounded`,
`TangleWorkbench.grow_until_arclength`, `FixedPoint.advance_key`, `FixedPoint.per_step_beta`,
`FixedPoint.preserves_orientation`, `budget.PointBudgetExceeded`, `ManifoldMachine.grade_seed`.

## A.2 Stage 1 — the crossing set N

**Definition.** Every crossing of an unstable polyline with a stable polyline is found by ONE
procedure on every path (first compute, recompute, blast image): each unstable segment is queried
against the index of the OTHER stability only, so no same-stability pair is ever formed; each hit is
resolved once, deduplicated by branch keys and by both canonical distances to the relative tolerance
`DEDUP_RTOL` (an absolute slack would merge distinct deep crossings, NUM-8), and bracketed on the
current polyline. A run of crossings on one segment gets a separator after each crossing's lower
bracket, and a `BridgeId` runs in chain order (NUM-19): the cut sorts a branch's crossings by their
lower bracket's chain position, then by cdist key, since a lerped cdist can step back by an ULP across
a fold of tied nodes (2026-10-30). Crossings are registered along each
unstable branch in unstable-cdist order, so registry ids are reproducible within and between
processes (INT-7).

Every detection is registered and reported, but N is the RESOLVED subset
(`IntersectionRegistry.resolved_ids`, plan Q13 (a), B1). A detection is in N iff

1. its `resolution_ratio` (separation from its brackets over the larger chord deviation) is at least
   `resolution_ratio_threshold` (default 3: ratio 1 is the derived bound where the displacement
   equals the polygon error, and 3 is a margin; the corpus minimum over the law cases is 3.04, on
   `heteroclinic`), and
2. it is in no RUN: no other detection of the same two curves lies closer to it along the unstable
   polyline than their chord deviation (`IntersectionRegistry.run_separation`). A near-coincident run
   is one unresolved stretch whatever its parity.
3. it clears the NEAR-ANCHOR GUARD (anchor-depth, 2026-10-29): its distance to the orbit point of
   either of its branches is at least `anchor_resolution_factor * max(h, xy_resolution(P))`
   (`IntersectionRegistry.orbit_clearance`). Near the orbit point the pip orbit contracts by the
   stable multiplier, and a polyline straying `h` puts an error of about `h / d` on each step's
   contraction; the trusted band allows 20%, so `c = 1 / (1 - 0.8) = 5`
   (`ANCHOR_RESOLUTION_FACTOR`, on by default, None off). The floor at the xy resolution matters only
   with the graded mesh (A.1), whose `h` falls below what float64 resolves near `P`. The registry
   holds the factor (`IntersectionRegistry.anchor_resolution_factor`; the workbench sets its own on
   every registry it holds and persists it as an optional setting), so every reader of N reads one
   set. The plan's second half, a detection closer to its stable-order neighbour than `c` times the
   pair's larger deviation (both out), is NOT part of N: measured, it takes out real, well-resolved
   crossings (ratio 4000-6000) of a tight fold on `heteroclinic`, three on one unstable curve, and on
   the k28 series a crossing of the original trellis whose stable neighbour is a deep blast crossing of
   another curve; an arc then loses an odd number of crossings, the arrangement's regions overlap,
   the heteroclinic case loses its heteroclinic crossings and a hole bridge's signs stop
   alternating. Whether a pair rule returns, and in what form, is the author's call.

An unresolved detection is not a point of N and not of X: no bridge ends at it, no iterate link
starts or ends at it, the arc across it is one arc, and `TangleWorkbench.report_resolution` names it
with the remedy "refine" (lower `area_cutoff` or grow the stretch). The registry decides over ALL its
detections, because on the blast path one run can straddle two image curves. This decides NUM-4 and
leaves NUM-3 as the resolution floor: from blast 9 (period 1) or 15 (period 3) sub-ratio detections
appear that only refinement can resolve.

**The anchor.** A periodic point is ONE crossing of all its branches (D4): one synthetic anchor at
canonical distance (0, 0) per periodic point, registered on every branch attached to it and stored
under its home pair `(u(o, 0), s(o, 0))` with that pair's crossing sign. An inversion point has two
stable and two unstable branches, each with its own left and right; the left of stable branch 0.0 is
the right of 0.1. A cross-side flag applies only to a pair on ONE branch. `Intersection.is_anchor`
tells an anchor by identity, never by a cdist slack.

**Bridges.** A bridge is the unstable arc between two consecutive crossings of N on one unstable
branch, cut at the exact registry crossings; every bridge's invariants are asserted where it enters
the workbench (NUM-14): a non-anchor end lies on the bridge's own branch inside its node span to
`DEDUP_RTOL`, an anchor end comes first and is the bridge's own periodic point.

**Resolution flags.** `TangleWorkbench.report_resolution` logs ONE WARNING for new unresolved
detections and one for new same-stability polyline crossings in the padded crossing box (INFO
outside it), each finding once. A blast only COUNTS the same-stability hits its new image bridges
take part in, because a blast adds no stable curve (C5). Same-stability polyline crossings (two
unstable or two stable curves, or a curve crossing itself) are polygon artifacts of under-resolution,
never crossings: the fundamental invariant forbids them.

**Implemented by.** `Tangle.resolve_crossings`, `Tangle.create_bridges`,
`TangleWorkbench.compute_intersections`, `TangleWorkbench._register_anchors`,
`TangleWorkbench._register_bridge`, `IntersectionRegistry.resolved_ids`,
`IntersectionRegistry.run_separation`, `IntersectionRegistry.orbit_clearance`,
`IntersectionRegistry.unresolved_crossings`, `geometry.xy_resolution`,
`TangleWorkbench.report_resolution`, `diagnostics.same_stability_crossings`,
`diagnostics.count_same_stability_crossings`.

## A.3 Stage 2 — the iterate table

**Definition.** A partial bijection `N -> N`: `y = M(x)` iff `y` lies on the advanced branches
(`advance_key`), both its brackets (the node cdists either side of it on each curve) meet `x`'s
carried one step by the manifolds' `stretch_param`, UNIQUELY, and the real map confirms it within the
image of `x`'s unstable segment plus `y`'s own (one batched call per pass). The two cdists are
compared individually, never their product (two different iterate chains can share a product).

Only crossings of N are linked. A source whose brackets meet several images, and every claimant of
an image named by two sources, is an UNRESOLVED LINK: recorded in `registry.unresolved_links`, never
linked, never resolved by registry order (B8); the set is disjoint from the table's `+1` sources. The
inference is registry-wide and runs in `compute_intersections`, in every `iterate_bridge` and at the
end of every `blast_zone`, so a crossing born early in a blast is linked once its image is registered
(BCG-3). Topology only looks links up (A.0).

On a heteroclinic link `W^u(A) x W^s(B)` each side scales by its own saddle: the unstable cdist by
A's factor and the stable by B's. The `p1(u) x p3(s)` recipes still fail the beta and area laws on
their heteroclinic links (NUM-20); no recipe in the corpus depends on those links.

**Implemented by.** `infer_iterates`, `TangleWorkbench.infer_iterates`,
`IntersectionRegistry.unresolved_links`, `IterateTable`, `Trellis.iterate`.

## A.4 Stage 3 — the strong pip and the trellis of the pip

**Definition (strong pip).** q0 is a strong pip iff no crossing of q0's fixed point, mapped onto
q0's stable branch by the minimal number of backward steps in the branch cycle (table first, else
scaling), lands strictly inside the open box `(0, c_s(q0)) x (0, c_u(q0))`; q0's own orbit is
recognised by id (table) or by a 2-D collision (`same_point`, `collision_rtol`). One
representative per crossing suffices because the box lies inside the computed trellis. An unresolved
detection is no crossing and disqualifies nothing (`is_strong_pip` skips it). The candidates are
classified on first read; restoring a remembered pip tests that one crossing only (F7). The default
pip is the candidate of smallest unstable cdist.

**Definition (the trellis of the pip).** Choosing q0 CUTS the trellis to
`T(q0) = (W^u computed, W^s[P, q0])` (D10; D12 reversed 2026-10-31, the author's rule: use every
computed manifold): stable branch `j` ends at c_s(f^j q0), `j < k_value`, and no unstable branch is
cut, so every computed unstable crossing whose stable end lies inside the stable cuts is a crossing
of N. `f` maps its stable half into itself and `f^-1` the computed unstable set (every blast child's
parent is held), so it is a genuine trellis. A deeper crossing may refute a reference pair and the
class count may then fall under blasting; the entropy may not (the series law
"entropy_never_decreases"). D12's unstable cut at the `k_value` deepest orbit points hid such
crossings (well resolved, ratio 128-9650, persistent at 10x and 100x finer `area_cutoff`) and froze
the topology whenever the pip orbit stalled (k28 blast 8: 4 classes cut, 15 uncut; k32 blasts
8-10: 13 cut, 34 / 34 / 30 uncut with h 0.540 -> 0.546 -> 0.548). The cut is the ONE
place N is decided for every later stage. A crossing past a cut is a point of X, not of N: it leaves
every branch, and the unstable arc across it joins its two bridges into one (`Bridge.joined`, a
bridge the workbench never holds; resolve bridges through `Trellis.bridge_between`).

**`m`, the depth.** `m` is the registered `+1` chain of the CHOSEN pip, as far as it consists of
crossings of N (N-a, plan Q2, B3); `BlastResult.pip_depth` reports it per fixed point. It is a
report, no cut: an unresolved orbit point stalls `m`, never the trellis. One consequence:

- a blast child ending on the chosen pip's composed `+1` orbit is exempt from the `min_separation`
  guard (one INFO), so the guard never decides the depth. The depth is decided by whether
  `f^{m+1}(q0)` is detected and resolved; the three ways it stalls are open (O). The k32 blast 9
  stall is an UNRESOLVED `f^10(q0)` (ratio 1.48), not a child dropped by `min_separation`.

A cut whose orbit point the table does not hold is PLACED by `image_cdist` and flagged in one
WARNING, never left uncut (PK-6). A crossing is past a cut by any amount: no absolute slack. Restoring
None uncuts. At an anchor, an arc of `u(o, b)` meets `s(o, b)`: the same orbit point and branch
index, equivariant under the map, and the pair carries the anchor's own sign (C7,
`Trellis.stable_key_at`). That every crossing of a cut trellis lies inside its cuts holds by
construction and is checked by the law "trellis_is_cut_at_the_pip_orbit".

**Several fixed points.** A trellis holds a crossing only when BOTH its branches are selected: a
tangle's own trellis holds no heteroclinic crossing (its arcs are joined across them), the
selection's does. Cut at every pip it is `T(q0_A) ∪ T(q0_B)`; `Trellis.pips` records the pip per
fixed point.

**Implemented by.** `StrongPip.is_strong_pip`, `Trellis.strong_pip_candidates`,
`Trellis.cut_at_pip_orbit`, `Trellis.image_cdist`, `Trellis.stable_key_at`,
`Trellis.unstable_key_of`, `Trellis.bridge_between`, `blast_zone`, `BlastResult.pip_depth`.

## A.5 Stage 4 — pseudoneighbor pairs

**Definition (the pseudoneighbor PDF).** Consecutive crossings `x, x'` on a stable branch of `T(q0)`
are a pair iff neither the open stable interval `W^s(x, x')` nor the open unstable interval
`W^u(x, x')` holds a point of X. Both intervals are tested by ONE function, `_interval_is_clear(...,
stability, ...)` (plan Q1, B2): for each crossing `r` and each `n` that lands `M^n(r)` on the interval's
branch (closed form over the branch cycle), the landing is the table's crossing when linked (exact,
by id), else `r` read against the pair's registered images in `r`'s own frame, else placed by
scaling with the 2-D collision guard for the pair's own orbit. A branch return is `k_value` map steps
and scales a cdist by `beta^k_value` (INV-2). A pair touching the anchor is no pair (the anchor is a
point of N on the arc).

The stable clause adds robustness, not new pairs: under the cut the reference window is one lap of
stable cdist on the pip's branch, so the orbit of every crossing of N visits it once, by backward
steps, as a registered crossing when resolved. On the whole corpus it fired 0 times in 182 series
steps. The "orphan" window pairs on the OTHER branches of a period-k orbit (NUM-18) were never
observed to die under it, because those windows are no longer walked.

**References and trajectories.** Reference pairs are those in the window `(c_s(f^k q0), c_s(q0)]`
on the pip's branch; trajectories are their table iterates while inside `T(q0)`
(`extend_pseudoneighbor_trajectories`). Forward pairs are recorded (they are real pairs) but never
holed (A.6). Deep crossings that fold onto the window replace the references, and the new reference
punches about `log(c_u) / log beta` holes (BCG-4): a discovery, not a bug.

A HETEROCLINIC pair (two consecutive crossings of `W^u(A)` on B's stable branch, inside B's window)
is a pair of B's tangle, tested with the cycle, `k` and beta of the ARC's fixed point (D7); it gets a
hole in B's partition and its backward chain runs along A's unstable cycle into A's bridges, so it
exists only on a trellis holding both tangles (the selection's).

**Superseded.** D5 (keep the `+1 .. +(k_value - 1)` iterate rule on the other branches of a
period-k orbit) is SUPERSEDED by D15: no forward iterate gets a hole on any orbit. With it the PK-3
WARNING (a branch whose `+j` appearance is missing) is superseded too: a missing `+j` costs nothing
and is not warned about.

**Implemented by.** `Pseudoneighbor._interval_is_clear`, `Trellis.compute_pseudoneighbors`,
`extend_pseudoneighbor_trajectories`, `same_point`.

## A.6 Stage 5 — holes

**Definition.** Holes are combinatorial (B4a, B4b). Each REFERENCE pair's bridge gets a hole beside
it, on its `bridge_side`, read off the crossing signs. A hole is the data
`Hole(near_intersection_id, pair, bridge_side, bounding_ids, openings, iterate, origin)`: the pair
(None for a propagated hole), its side of its OWN bridge in the bridge's dynamical orientation, the
bridge's two crossings, the stable halves it opens, its step along the reference orbit (0 for the
reference, negative backward) and the reference it descends from. Where it is drawn is plotting
(`plotting.hole_point`).

**Holes propagate BACKWARD only** (the author's rule, on every orbit, D15). Only references are
punched directly; `punch_holes` asserts no backward pair is ever recorded (a reference lies in the
window `S \ f(S)`, so a preimage leaves T at once). A forward iterate gets NO hole even when the
trellis has registered its image crossings and a bridge spans them: its region is the image of a
hole still attached to the stable manifold. The other branches of a period-k orbit get their holes
from the backward chain, which walks the whole unstable branch cycle (`_is_forward_iterate`, law
"no_forward_hole").

**Propagation.** `propagate_reference_holes` maps each reference bridge backward: both crossing ids
move by `Trellis.image_cdist` onto `advance_key(key, -1)`, and the bridge whose consecutive crossings
bracket that span (a bisect over the branch's crossing order, `_containing_bridge`, PERF-4) gets a
hole on the WHOLE bridge, by design: the image sub-arc's far end is a crossing with `f^-1(W^s)` past
the trimmed stable manifold, so the trellis has no node to cut it at (INT-2). The chain stops when a
container repeats (periodic) or none exists (the image left the trellis); it needs no step cap
because every step visits a new bridge (PK-9). The side is CARRIED from the reference, flipped every
step only under an orientation-reversing map (`side_after`, the one rule; det J < 0 stays out of
scope, C3). Containment is decided in the canonical coordinate; a SCALED span end within
`SCALING_RTOL` of a crossing of the container's branch is ONE WARNING per branch, because the
container could be its neighbour (NUM-5). A table end is exact and not watched. No dedup, no
fallback.

**Signs decide** (D6). At a crossing with `sign = crossing_sign(x)` the hole on `bridge_side` opens
the `outward` half iff `(sign > 0) == (bridge_side == "left")`; at the anchor the signs pick the
branch (`_hole_openings`, via `Trellis.stable_key_at`). No polyline decides an opening; the nearest
vertex of the whole polyline is never read (INT-4).

**Geometry is a flag.** The geometric readings live in `topology/diagnostics.py`: the polyline row of
a bridge end (`bridge_row_diagnosis`; on a 2-4 node bridge the approach is the vector to the far
node, NUM-1) and the measured side of a propagated hole (`diagnostics.measured_sides`, the reference
point carried back by the real inverse map, read on the image sub-arc, which escapes the lobe when
the reference lobe is under-resolved, NUM-2). A disagreement is one WARNING per kind per
`Trellis.punch_holes(diagnose=True)` (`geometry_flags`), never an assertion; `side_source(hole)`
says whether a side was read off the signs or carried. The unreachable class and virtual class of
the deep period-3 runs were geometric openings (PK-1); with the signs deciding, every class walks.

**Implemented by.** `punch_holes`, `propagate_reference_holes`, `_containing_bridge`,
`_hole_openings`, `_is_forward_iterate`, `side_after`, `diagnostics.geometry_flags`,
`diagnostics.measured_sides`, `diagnostics.side_source`, `plotting.hole_point`.

## A.7 Stage 6 — the stable partition

**Definition (the partition PDF + D6).** Each stable branch is cut at the endpoints of every hole
bridge, plus the anchor and the branch's outermost crossing. An interval is open at an end iff a hole
opens that half there; a point excluded by both neighbouring intervals is a closed singleton
`[x, x]` (two holes flanking one bridge pinch both its crossings; the fixed point and the strong pip
pinch by the same rule); each point belongs to exactly one piece. A stable branch with no crossing
past its anchor is ONE element `[anchor, end]` (C1). Rows: `row(first) = left iff sign > 0`,
`row(second) = left iff sign < 0` (`row_of_end`).

**Ownership is exact** (C15). Every partition boundary is a registry float, so a crossing is owned
by id at a boundary and by its exact stable cdist elsewhere. The topology reads no absolute cdist
slack anywhere: anchors by identity, ownership and containment exact, and the only slack is the
RELATIVE error a scaled cdist carries (`SCALING_RTOL` at a scaled landing end, `collision_rtol` in
`same_point`). `cdist_tol` survives only as geometry (degenerate bridge cuts, clipping a
polyline to its arc).

**One reader of a bridge end's element.** `PartitionFamily.owner_at_end` reads the owner of a bridge
end on its row, on the branch `Trellis.stable_key_at` names; `owner_of_intersection` reads any side
of a named branch. Every computed path (classes, landings, the iterated partition) goes through them.

The invariants asserted here: (I1) all holes of one origin share `bridge_side` (parity-flipped under
an orientation-reversing map), `check_holes_share_bridge_side`; (I2) a hole bridge whose crossings
share a stable branch has alternating crossing signs, `check_bridge_rows_consistent`; one owner per
crossing per side.

**Implemented by.** `partition_stable_manifold`, `_hole_marks_on_branch`, `intervals_from_marks`,
`row_of_end`, `rows_of_bridge`, `HomotopyPartition`, `PartitionFamily.owner_at_end`,
`PartitionFamily.owner_of_intersection`, `PartitionFamily.interval_at`.

## A.8 Stage 7 — bridge classes and inertness

**Definition (class).** The class of a bridge is the UNORDERED pair `{X, Y}` of partition elements
its two ends sit in on their rows (`owner_at_end`), oriented anchor outward (`anchor_outward_key`):
the homotopy class of arcs between `X` and `Y` relative to the stable manifold and the holes. A
member traversed in the unstable direction runs `source -> target` (symbol `a`) or back (`a^-1`).

**Loops fold forward** (the author's rule, D2). A bridge with both ends in one element is the image
of a bridge that did connect two elements: it joins the class of its IMAGE-CHAIN PARENT, the one
bridge whose `_image_chain` holds it (one parent per pair, asserted: images of disjoint arcs are
disjoint), climbing through loop parents. Walking preimages instead breaks at every fold-born
crossing (BCG-2). A loop no chain reaches is an ORPHAN, logged once and left unclassed
(`BridgeClassTable.orphans`).

**Evidence and alphabet** (D2). `bridge_classes` is the EVIDENCE table, every bridge's class, grouped
by TANGLE (tangle 0's classes, then tangle 1's, then the connecting ones; unstable cdists of
different saddles are never interleaved) and within a group by each class's smallest member
first-crossing unstable cdist, so the class holding a tangle's anchor bridge leads its group (the
author's rules). The ALPHABET is its restriction to the minimal trellis's kept bridges
(`MinimalTrellis.alphabet`, `restricted_to`); a blast child is image evidence only.

**Identity and letters.** Element ids are positional within one recompute; across recomputes an
element is its bounding crossings `(branch_key, side, lo_id, hi_id)` (`PartitionInterval.identity`;
C2, refined by C10), and a class is the pair of its elements' identities (`Topology.identity`). The
session's `BridgeAlphabet` letters ACTIVE classes by that identity, so a class keeps its letter while
its elements keep their bounds (BCG-7). One letter generator serves both series over disjoint digit
sets: active `a..t, aa, ...`, inert `u..z, uu, ...`, no skipping (`letter`, SIMP-6).

**Definition (inertness)** (plan Q3, I-a, B5). A class is INERT iff it is TRANSIENT in the graph of the
class-level words: no cycle is reachable from it, i.e. its word is empty or names only inert classes
(`inert_classes`). An unresolved class (no word) is active, and so is every class whose word names
one. `BridgeClassEntry.inert` is None until `symbolic_dynamics` sets it; only then are classes
lettered. The minimal trellis decides BEFORE the words which hole bridges to map, by the zeroth-hole
rule (A.8b). Inert classes keep their symbols in the words and are no nodes of the transition graph
(A.9).

**Implemented by.** `bridge_classes`, `oriented_class`, `anchor_outward_key`, `_image_chain`,
`BridgeClassTable.restricted_to`, `PartitionInterval.identity`, `Topology.identity`,
`BridgeAlphabet.letter_for`, `letter`, `inert_letters`, `inert_classes`.

## A.8b Stage 8 — the minimal trellis

**Definition.** The minimal trellis keeps every hole bridge, and maps every hole bridge EXCEPT one
carrying a ZEROTH (reference) hole (the author's rule, B5): the bridges tiling a mapped bridge's first
forward image are kept too. The zeroth-hole rule, read from `holes`, decides which bridges to map;
the class table only places them (`classes.entry_of`): a hole bridge in no class is kept with a
WARNING but not mapped. A zeroth-hole bridge's lobe holds no point of X and one marker (its own hole),
whose image is no marker, so its image is homotopic to a stable segment and its class word is empty;
the law "zeroth_hole_classes_have_empty_words" checks it. The rule is SUFFICIENT, not the definition
of inertness: the heteroclinic class `u` is inert with no hole at all, and transience over the words
stays the verdict (A.8).

The stable manifold is kept whole; the surviving nodes are the kept bridges' endpoints, the anchors
and each stable branch's outermost crossing. Its arrangement is SPARSE: unstable arcs come from the
bridge list, and an absent unstable ray is an absent wall (the two sectors merge), not a slit.

An image pair no `Bridge` spans is manifold that was never computed: two consecutive crossings on two
different computed curves are a gap, not an arc (`Trellis.is_computed_arc`, INT-9). It is skipped
with a WARNING, never invented; the dense arrangement draws its ends as dangling stubs, because a
closed arc there would make zero-area faces and break Euler's formula (NUM-16). A period-k branch
whose reference images were never computed simply has no window hole yet (D9, the stagger accepted).

**Implemented by.** `minimal_trellis`, `MinimalTrellis.kept_bridge_ids`, `MinimalTrellis.alphabet`,
`MinimalTrellis.inert_hole_bridge_ids`, `image_chain`, `Arrangement`.

## A.8c Stage 8 — the word of a class

**Definition (the object).** Let `b` be any member of class `C = {X, Y}`. `f(b)` is an unstable arc
from `f(X)` to `f(Y)`; it meets the stable manifold `S` only at its ends and inside the window
`S \ f(S)`. The ITINERARY of `f(b)` lists, at each crossing, the element on the entry row then the
element on the exit row; consecutive pairs are classes or their inverses; pairs in one element are
loops and drop, and an adjacent `c c^-1` cancels (free reduction: the two sub-arcs are homotopic with
sliding ends). The WORD `w(C)` is that reduced itinerary, a class invariant.

**Definition (the algorithm).** The shortest walk of the dual graph IS the word (D3, the author's
rule); a registered image chain is only its cross-check. A walk must exist; an unreachable walk is a
construction bug to root out, not a state to model. The walk:

- runs from the stable nodes the landing SPAN `f(X)` meets to those `f(Y)` meets, not every node of
  the iterated element (C14): the hole bridges bounding `X` map onto the kept image bridges at the
  ends of `f(X)`, which `f(b)` cannot cross. A landing span overlapping no stable node is unresolved
  ("span_meets_no_node", B6a);
- alternates face, unified node, face, crossing `S` only at unified nodes (the window, A.4) and never
  an unstable edge, and records BOTH sides of every crossing, so the itinerary has even length;
- holds the marker law (D13, the author's rule): `f` is a homeomorphism and `W^u` never crosses
  `W^u`, so a member bridge and its stable segment bound a disk whose image holds the images of
  exactly the markers (periodic points, holes beside their own bridges) the lobe holds. It binds every
  class with a member lobe, not only ties (C13). Two member lobes holding different markers raise
  `MembersDisagree`, and the class is unresolved (B6a). Detours up to `MAX_DETOUR` crossings; none
  admissible is unresolved;
- respects WALLS (C17, under D13 and D14): a class's arcs lie on its members' unstable branches, or
  one map step past those of the class that named it; images of disjoint arcs are disjoint, so the
  decided walk of a class walked earlier on disjoint branches is a wall, and a tie loses every walk
  whose chords cross one (`walk_chords`, `DualGraph.chords_cross`). Every bridge is the same (D14): a
  class joining two branches or two tangles has no special case.

A search ends "resolved" (one decided walk; zero crossings included), "ambiguous" (a tie no rule
decides; the first candidate in deterministic order feeds the rest of the table, and the dynamics is
unreliable) or "unresolved" with one of six `WalkReason`s: "unreachable", "no_start_node",
"no_goal_node", "no_admissible_walk", "truncated" (the `max_walks` budget ran out: the walks found are
kept but decide nothing) or "span_meets_no_node". A class whose WORD is empty is "trivial"; a capped
search is unresolved, never a decided walk. An unresolved class is no wall.

**Singletons.** A singleton element owns no dual-graph node and cannot walk: its word is read off the
registered image chain of its first non-loop member (`trellis_itinerary`); no registered chain is
unresolved ("grow or blast").

**Closure** (C9, D11, plan Q9). Every pair of a walk is the class of one piece of `f(B)`, so a word naming
a class outside the alphabet adds it, walked in turn, with `named_by`. A pair no computed bridge spans
is a class all the same (D11, the author's: "we do not need computed bridges to find words"): a
member-less entry walked from its elements' landings and flagged "no computed bridge yet"
(`BridgeClassEntry.spanned` False), informational only. The closure ends because element pairs are
finite, and it is a function of the dual graph alone.

**Free reduction** (B6b). `free_reduce` applies to the class word (`ClassDynamics.symbols`), to the
closure and inertness words and to the evidence word; `translate_itinerary` stays UNREDUCED (Theorem
T reads it). It cancels nothing on the current corpus.

**Evidence** (B6b). `diagnostics.verify_words` compares each member's registered image, read with
the same translation and reduction, at WORD level, once every word is final, and never changes one.
A disagreement carries an `EvidenceDiagnosis` (INT-1: the word is the walk, not the registered
chain): "stable end" (two consecutive chain crossings on one branch share a sign; the image passed
the end of the computed stable branch) is INCOMPLETE evidence (C11); "order flip" (a member end and a
bound of its own element have registered images in reverse stable order, which the map forbids) is
NOISE (C21); neither gives a verdict. "walk" is a WARNING and counts against `is_reliable`. A landing
that straddles a cut still walks from its midpoint owner but counts against `is_reliable` (INT-8).

**Dual graph.** One `FaceNode` per merged face (the outer face of a component is read off the
rotation system: the face whose sector at the component's easternmost vertex holds due east, never an
area, INT-10); two open `StableNode`s per stable edge, one solid unified node on the fundamental
segment `S \ f(S)`.

**Implemented by.** `DualGraph`, `DualGraph.markers`, `DualGraph.enclosed`, `land_element`,
`shortest_walks`, `lobe_image_markers`, `holds_markers`, `walk_chords`, `trellis_itinerary`,
`WalkSearch`, `WalkReason`, `MembersDisagree`, `translate_itinerary`, `free_reduce`,
`symbolic_dynamics`, `diagnostics.verify_words`, `EvidenceDiagnosis`.

## A.8d Stage 8 — the iterated partition and the refined level

**Definition (P1).** `P1 = P ∨ f(∂P)`: each stable branch is cut, beyond the homotopy partition's
cuts, at every crossing of the kept image bridges (the first images of the hole bridges) with `S`,
by the EMPTY-STRETCH rule (the author's: "all of the stretched and folded shape is inclusive, the
empty parts between the bends are open"). A hole's stretch is empty, and so is the footprint of its
image lobe `f(H)`: its BASE `(x_0, x_m)` and one CHORD per pair of interior crossings. **Chords pair
in STABLE order on every stable branch, the base branch included** (B7b): `f(H)` holds no anchor, so
walking any branch from its anchor enters and leaves it in turn (C19's Jordan argument). An image
crossing belongs, on the row the bridge cuts, to the element AWAY from the empty stretch it abuts. A
crossing on the image chain of a hole with no stretch of its own (the anchor's one-sided hole, or a
lobe holding other holes' crossings) marks nothing ("unmarked lobe", INFO). Cross-branch pieces are
cut like any other (PK-4). Iterated elements record `parent_element_id` (exact containment) and
`cut_by`.

An image chain crossing a branch an ODD number of times has left the computed stable manifold there
(a "stable end"): its footprint is INCOMPLETE, the outermost crossing stays unpaired, and the branch is
recorded in `IteratedHomotopyPartition.incomplete_footprints`, data with one INFO, never a
reliability failure (plan Q14, B7).

**Superseded.** C19's base-branch half (chords in CHAIN order on the base branch, the author's
2026-09-23 rule) is SUPERSEDED by B7b; it misread one lobe, k32_series hole `(6, 7)`, which spirals
back across its base branch.

**Definition (refined symbols).** The walk of every class is read at the `P1` level: each token
carries its class and the pair of `P1` elements it lands on, its OCCURRENCE. The occurrences of a
class across all words are its refined children `c_1, c_2, ...` anchor outward (`refine`); members
are matched by their own endpoint owners (`SymbolicDynamics.match_members`), and an unmatched member
is informational.

**Theorem T.** Every member of a class has the same `P1`-level freely reduced itinerary (the strip
two members bound with the stable manifold holds no hole, so no `P1` cut separates their images).
Then the refined word of c_k is the class word with every token written as the child its
occurrence names: the word of every member, not an inheritance. The refined matrix M_1 is the
IN-SPLITTING of `M`: siblings have equal rows and the children's columns sum to the parent's, a
conjugacy with the same entropy. Both are LAWS, not on every build:
T ("members_share_the_P1_itinerary", in `LAW_CASES_ONLY`) runs on the law cases and after every
blast of the series ("series_members_share_the_P1_itinerary"), not on the deep builds; the
in-splitting ("refined_matrix_is_an_in_splitting") runs on the law cases and the deep builds, with
no series test.
Under the zeroth-hole rule `P1`, the minimal trellis and the refined classes are what they were under
the evidence rule on every build.

**Implemented by.** `IteratedHomotopyPartition.from_minimal`, `IteratedHomotopyPartition._empty_stretches`,
`IteratedHomotopyPartition.incomplete_footprints`, `Cut`, `refine`, `RefinedClass`,
`SymbolicDynamics.match_members`, `SymbolicDynamics.refined_rules`.

## A.9 Stage 9 — the transition matrix

**Definition.** Over the active classes, `M[C][D]` counts the tokens `D` or `D^-1` in the FREELY
REDUCED word of `C`; over the refined symbols it is M_1, the in-splitting of `M` (A.8d), shown by
default. Inert classes (transient, A.8) are not nodes; a token naming one is skipped, and an inert
word names only inert classes by definition. No spectral radius is computed or shown (the author's
rule, kept by plan Q6); it would be equal for `M` and M_1.

**Implemented by.** `SymbolicDynamics.transition_graph`, `SymbolicDynamics.transition_matrix`.

## A.10 Stage 10 — recompute and blasting

**Definition (D1, the author's rule, restated by plan Q8).** The topology is a pure function of
(registry, bridges, pips): any growth or blast means a full recompute (`TangleSession.topology`
holds one `Topology` per generation, selection and pips; SIMP-8). Grow, blast and recompute commute:
there is no collapse path (N is the registry's verdict on every path, A.2), and a blast's image curves
are computed manifold, kept by the workbench, re-indexed query-only by every recompute until growth
merges them into a branch and cut like any curve (one crossed fewer than twice is a partial piece),
so a recompute keeps every crossing, bridge, partial piece and `iterated` flag (NUM-12, C20). A
blast only adds computed manifold, and what it finds is a discovery and may change anything, the
class count included (it may fall when a deeper crossing refutes a reference); the series laws
measure that the entropy and the depth `pip_depth` (a report, A.4) never decrease. D1 does not bind
the closed alphabet.

**Blasting.** `blast_zone` iterates the bridges whose INNERMOST zone it is (C4): a bridge belongs to
the smallest zone holding its midpoint, so an outer blast never iterates a bridge inside an inner
zone and nested blasting is order independent. One exception (anchor-depth, 2026-10-29): an
un-iterated bridge of a selected fixed point ARRIVING at `f^j(q0)`, `j >= 1` (its second crossing on
the registered orbit of the blasted zone's own corner pip `ResonanceZone.boundary_intersection`,
not the session's chosen pip, which may differ, read once per iteration), is the zone's whatever
its midpoint says, in the initial frontier and among the children: near `P` the frozen zone polygon
cannot place a thin lobe's midpoint. Only the arriving arc: `f^j` maps the zone's corner at `q0`
(between `W^u[P, q0]` arriving and `W^s[P, q0]`) onto the corner at `f^j(q0)` on the zone's side of
`W^s[P, q0]`, the side an orientation-preserving map keeps, so the arc arriving lies inside and the
one leaving lies outside (at `q0` it is the turnstile's exit lobe). The plan's "ends on a crossing of
the orbit" iterated that exit lobe and changed every blasted law case. Other fixed points' pip-orbit
bridges keep the innermost-zone rule; `classify_bridge` and the class zones are unchanged.
`BlastResult.trusted_depth` is a flag beside `pip_depth`: the deepest `j` to which the registered
orbit contracts like the real one (inside `TRUST_LINEAR_RADIUS = 1e-3` each period divides the
distance to the orbit point, read off `coords`, by `per_step_beta("stable") ** -period` within
`TRUST_BAND = (0.8, 1.25)`; the test of `scripts/henon_k25_errors_figures.py`), logged at INFO when
smaller; nothing reads it. Any failure to map a bridge forward propagates
(NUM-10). `max_points` is an opt-in budget checked between iterations, a clean resumable stop.
A `min_separation` drop is remembered on the bridge (`Bridge.dropped_at`, carried by
`rebuild_bridges` beside `iterated`), and the zone's later calls skip it unless they pass a smaller
separation or None; the guard's reference is every curve the workbench holds, dropped and exterior
children included. So `n` single-iteration blasts equal one `n`-iteration blast (blast-steps-compose,
2026-10-21). The separation follows the resolution (the author's rule, same run):
`min_separation="resolution"` is `RESOLUTION_SEPARATION = 1e-5` at the default cutoff, scaled as the
polylines' chord deviation, `area_cutoff ** (2/3)` (measured on k=2.5: the median deviation falls by
`10^(2/3)` per decade). A finer cutoff then keeps curves that run closer; the recipes and law cases
keep their explicit constants. `workers` changes no blast result: the map batches and the
resolution sweep's pair blocks run on that many cores, and every verdict and count is the serial
one (cpu-workers, 2026-10-29).

**Lean blasting (lean-blast, 2026-10-29).** `blast_zone(lean=True)` frees, at the END of each
iteration (the author's choice: nothing can be freed before that iteration's topology is known),
the curve neither the topology nor a later blast reads. It links the iterate table, computes
`session.topology()` of the selection (and of every fixed point when the selection is not all of
them) and calls `TangleWorkbench.release_curves(keep)`. Kept: every bridge on a grown branch, every
un-iterated undropped bridge some stored zone (or the blasted one) holds, the minimal trellis, the
pieces of every joined arc a trellis holds, every bridge ending on a pip-orbit crossing, and the
LINEAGE of the kept frontier: its ancestors (the bridges its rows' preimages lie inside) and every
bridge their images tile. The lineage is not in the plan; it is forced: a later reference pair is
cut from a frontier descendant and its holes propagate backward through those ancestors, whose
first images the minimal trellis maps (k32_series: without the ancestors a later topology reads a
freed hole bridge, without their images a freed image bridge). Rows a reader still needs are never
freed: a kept bridge's rows and one neighbour past each end, each crossing's four bracket rows and
one neighbour of each, the ends of every curve, and the images an un-iterated kept bridge already
has through a straddle neighbour (the run of their chain and its phantoms, which
`ManifoldMachine._image_chains` merges into and refinement reads). A released bridge keeps its ends,
id, key, flags and width, its midpoint node (`Bridge.test_point`, zone membership) and its point
count; a row read raises `CurveReleased`. A lean workbench (`TangleWorkbench.lean`) refuses whatever
would re-read freed curve: recomputes, cuts, growth, trims, zones, save, the same-stability readings,
the dense arrangement and a blast that is not lean. The geometric guard reads every interior, so a
lean blast takes `min_separation` None or `"width"` and never diagnoses. Measured on the blast series
(the lean-vs-full law, `tests/regression/test_lean_blast.py`): unguarded, every lineage survives and
lean frees 0.2-1.3% of the rows; with the width guard, whose drops end lineages, the store ends
2.6-4.8x smaller. On the k=2.5 scale series (`profile_numerics.py --cases scale`) lean unguarded
freed 0.5% (1e-7) and 0.01% (1e-9) and grew 1.7x per blast as before; lean width ended at 24k / 320k
rows against full width's 116k / 1.13M (peak RSS 141 / 290 MiB against 165 / 591). Lean mode buys a
constant factor, and only where a guard ends lineages; it is not a new exponent.

**Representative blasting (class-representatives, 2026-10-30).** `blast_zone(representatives=policy)`
iterates one bridge per refined class instead of every bridge: by Theorem T (A.8d) the members of a
class share their reduced itinerary, so iterating all of them mostly repeats one fact. Each step
links the iterate table, computes `session.topology()` of the selection and sorts the step's
candidates, its frontier plus every bridge held earlier in the call (`Blast._representatives`).
Always iterated: a bridge ending on a pip-orbit crossing, one the minimal trellis maps (kept and not
a zeroth-hole bridge: its image feeds the words), and one the class table does not know (a partial
piece, a bridge outside the cut trellis, an orphan). Every other bridge goes into the bucket of its
refined child (`member_refinement`), else its unmatched pair (`unmatched_members`), else its class;
a loop by its class and element. One bridge per bucket is iterated, by the policy: `fewest_points`
(ties by smallest id), `smallest_id` or `largest_id`; the others are HELD (`BlastStep.held`, held ->
representative; `BlastResult.held`, at the call's end; data only). Nothing is remembered on a held
bridge (the author's choice): each later step of the call re-offers it with the new frontier, and a
later call sees an ordinary un-iterated bridge, so `n` one-step calls equal one `n`-step call and a
full blast afterwards takes every held bridge. It composes with `lean` (a held bridge is an
un-iterated, undropped frontier bridge, kept by `_lean_keep`), `min_separation` (a held bridge was
never dropped), `max_points` (the prediction counts only the bridges to iterate) and `diagnose`;
`None` is the full blast byte for byte (identity probe). Measured
(`scripts/blast_representatives.py`, `figures/blast_representatives/`): against the full series on
every `cases.SERIES` and the blasting law cases, the alphabet, the closed alphabet, the words, the
refined matrix, the holes, the references and both depths agree at every step except k32_series
blasts 6-9 (`fewest_points`; 6 under `smallest_id`, 6-10 under `largest_id`), where the full blast
switches its reference pair one step after the witness crossing appears on the image of a bridge the
representative blast holds, and the representative blast switches when that bridge is re-offered
and iterated; the pseudoneighbor TRAJECTORY pairs are thinner (1-4 fewer) on every other series
from blast 4, never a reference. On the k=2.5 scale series at 1e-7, blast 12, the saving is 1.9x (1e-5 guard) to 2.6x (none) in bridges and points, not the
10x the plan hoped: once the pip depth stops (7, the near-anchor guard), the frontier lies past the
cut trellis, is unclassed and is always iterated, so it still grows like `e^{hn}`. At 1e-9 to 1e-12,
the separations 1e-5, `"resolution"` and None, all representative, agree on every compared item at
their common trusted depth and at their last common blast (A.4's depth equals the trusted depth
there); the prune hides nothing at trusted depth on that series.

**The width guard.** `min_separation="width"` needs no point cloud. A bridge's width is seeded, when
it is first iterated, as `|area| / L` of the ring it and its stable segment bound (`L` its unstable
arc length; a bridge whose crossings lie on two stable branches has no ring and is never dropped),
and its children carry `w * per_step_beta("stable")` (the lobe contracts across as it stretches
along: area is preserved; the plan's `w / beta_s` read the stable factor as a magnitude above 1). The
factor is the saddle's whose stable branch holds the seed ring (`Bridge.width_saddle`, carried with
the width): a lobe of `W^u(A)` based on `W^s(B)` contracts by B's, as NUM-20 reads each side of a
heteroclinic link by its own saddle. A
child is dropped when `w < WIDTH_SEPARATION_FACTOR * h`, `h = 1.5 * area_cutoff / L` with `L` its
median chord: a lobe thinner than its polygon error is not resolved. The pip-orbit exemption holds,
a width drop is remembered as `dropped_at = WIDTH_DROP` (-1.0, so only `None` re-admits it), widths
are saved (a column only when some bridge has one), and `n` single-iteration width blasts equal one
`n`-iteration blast. Measured on the scale series against `min_separation=1e-5`: the same pip
depth (9 at 1e-7, 12 at 1e-9), 116k / 1.13M rows against 2.18M / 5.10M, blast 0.5 / 4.8 s against
54 / 71 s; its drops start earlier (blast 6 against 8 at 1e-7). `WIDTH_SEPARATION_FACTOR` stays at
the plan's 1.0: 0.01 matches the 1e-5 guard's first drop (blast 8) but drops nothing in the coarse
k28 test builds; the calibration is the author's call.

**Zones.** A stable trim is a VIEW (C22): growth re-finds every branch's physical end and so
un-trims a zone, the topology never depends on the trim (the trellis is cut at the pip orbit), and
`ResonanceZone.restore` never shortens a branch growth has run past. Zones, trims and letters feed no
topology output. Simplify freely toward the chalkboard design; aligning the public API with the
author's use is a later plan (D8).

**Implemented by.** `Topology`, `TangleSession.topology`, `blast_zone`, `BlastResult`,
`TangleWorkbench.compute_intersections`, `ResonanceZone.restore`, `TangleWorkbench.release_curves`,
`PointStore.compact`, `Bridge.release`, `Blast._lean_keep`, `Blast._ring_width`, `Blast._too_thin`,
`Blast._trusted_depth`, `Blast._representatives`.

## A.11 Session and API

The session has three concerns: hold the workbench, cache one `Topology` per (generation, selection,
pips), and remember pips and letters. There are no stage fan-outs (each stage takes ONE selection and
computes it once; `classify_strong_pips(fixed_point=None)` still runs per fixed point, to keep each
pip's memory) and no plot delegates on the session (B9a): `cartoon_zones` and `plot_resonance_zones`
live in `loom/cartoon.py` (exported from tanglepack.loom), `cartoon_enclosures` in
`topology/plotting.py`; per-fixed-point plots go through `session.trellis(fp).plot_*`. The cache key
has no depth.

**Pips.** The chosen strong pip of each fixed point is remembered as a crossing (cdists and keys),
never as an id, and found again on every rebuilt trellis (API-1). `resonance_zone` adopts its pip
only when it is a candidate on the trimmed trellis (else one WARNING, the previous choice kept). A
remembered pip that stops being a candidate is never replaced by the default; `topology()` raises
naming both ways out, and `forget_strong_pip` is the remedy (B8). Partitions, holes and pips come
from the SELECTED fixed points only, and stages 4-6 run on the selection's trellis (API-2).

**Zones and tangles.** One resonance zone per fixed point, keyed by it; an inversion zone is the union
of its two loops (C8). Fixed points are constructed OUTERMOST first, so the period-1 saddle is `A`
and tangle 0 (NEST-6). The inversion series loses classes at blast 6 by a genuine reference switch
(BCG-8), still open as DEC-25 (O).

**Implemented by.** `TangleSession`, `TangleSession.trellis`, `TangleSession.topology`,
`TangleSession.set_strong_pip`, `TangleSession.forget_strong_pip`, `TangleSession.chosen_strong_pip`,
`TangleSession.resonance_zone`, `TangleSession.classify_bridge`, `BridgeAlphabet`,
`cartoon_zones`, `plot_resonance_zones`, `cartoon_enclosures`.

## H. History

**The alignment slices** (plan `docs/alignment_plan_2026-10-09.md`, one per `/tangle_refactor` run;
outcomes in `.refactor/decisions.md`; the run folders `.refactor/runs/` were removed with the
Athanor workflow on 2026-10-28 and live in git history up to `ab213db`), with merge commits:

- B0 `oracle-laws` (`1029fcc`): oracle laws on the law cases; the k=10 full 2-shift pinned (plan Q11).
- B1 `resolved-crossings` (`76d0710`): N is the resolved set (plan Q13 (a)); the blast-path parity vote
  `Tangle._collapse_noise_crossings` deleted; the trellis cut at the deepest resolved orbit point.
- B6a `walk-bounds-honest` (`8a1cdc3`): truncated searches, spans meeting no node and members
  disagreeing on markers are unresolved.
- B8 `unresolved-sets` (`ae18859`): contested and multi-match links go to `unresolved_links`;
  `forget_strong_pip`; `resonance_zone` adopts only candidates.
- B4a `holes-algorithm-vs-diagnostics` (`f166573`): the geometric readings moved to
  `topology/diagnostics.py`; `Hole.measured_side` and `Hole.side_source` removed.
- B4b `holes-are-combinatorial` (`ba23a1f`): references only; the direct backward punch, its two
  dedups, the fallback, `epsilon` and `_step_into` deleted; the backward span keeps `image_cdist`'s
  table flag (`from_table`) so the one NUM-5 WARNING watches scaled ends only.
- B2 `pseudoneighbor-both-intervals` (`9cd29d5`): both open intervals tested (plan Q1, S-a).
- B3 `depth-of-n` (`7d3a294`): N-a, the `min_separation` exemption, `BlastResult.pip_depth`. N-b (an
  explicit depth parameter) was offered and NOT chosen.
- B5 `inertness-is-transience` (`3fb310a`): transience plus the zeroth-hole rule; the evidence fixed
  point (`_resolve_inertness`, `image_classes`, `image_loops`) and PK-8's WARNING deleted.
- B7 `refined-level-defined` (`a6a4f68`): Theorem T and the in-splitting stated and made laws;
  `incomplete_footprints` replaced a WARNING.
- B7b `base-chords-in-stable-order` (`63477e3`): stable-order chords on the base branch too.
- B6b `word-finished` (`65edec4`): three statuses and six reasons; `free_reduce`; evidence moved to
  `diagnostics.verify_words`.
- B10 `dead-api` (`517aa01`): `grow_x_times`, the registry's networkx graph, `graphviz.py`,
  `BlastResult.terminated_early` and other dead API deleted.
- B9a `session-three-concerns` (`5f0db5f`): fan-outs, plot delegates and `partition_element_for`
  left the session; cartoon code moved.
- B9b `one-element-at-end` (`94fca6b`): `Trellis.element_for`, `_elements_at_bridge`,
  `StablePartitionResult.elements_at_bridge` and `BridgeClass._element_at_end` deleted; one reader,
  `owner_at_end`.
- B12a `one-design-note` (`a550a5d`): this note; campaign ids out of `src/`; `CLAUDE.md` rewritten to `HEAD`.
- B12b `laws-and-loose-ends` (`ad1d96a`): nine tautological laws deleted; `plot_resonance_zones` smoke
  test; `topology.__all__` = 51; rst lists and `conf.py`; chalkboard closed; recipe reasons; three
  figure folders; narration reworded.

**The work packages** of the open-issues round (2026-10-07 .. 10-09; table at
`docs/open_issues.md` "Work packages", outcomes after it):

- WP0: the issue catalog, the chalkboard, decisions D1-D8 recorded.
- WP1: laws encoding the decisions; the k=4 inversion fixture.
- WP2: exact iterate bookkeeping (inference after every blast, map-gated), `image_cdist` the one
  scaling path, relative dedup.
- WP3: inversion end to end: one anchor with four branches, anchor-only branch, zone ring by full key.
- WP4: signs decide openings and rows; geometry a flag; bisect container.
- WP5: the period-k hole gap closed under D5 (later superseded by D15).
- WP6: the alphabet from the minimal trellis; loops fold forward; one letter generator.
- WP7: the walk always exists; the alphabet closed under the word map.
- WP8: recompute semantics; reproducible ids; element identity by bounds; one `Topology` per state.
- WP8b: the trellis of the pip (D10, D12; D12 reversed 2026-10-31); words need no computed bridge (D11).
- WP9: nested and heteroclinic tangles; innermost-zone blasting.
- WP9b: the non-crossing rule (D13): landing spans, markers.
- WP9c: every bridge is the same (D14): walls.
- WP10: numerics and performance leftovers: one rtree per stability, bracketed cuts, bridge
  invariants asserted, outer face from the rotation system.
- WP10b: image curves are computed manifold (C20); foreign chords in stable order (C19).
- WP11: the simplification sweep; orientation a FixedPoint fact; no absolute cdist slack.
- WP12: re-verification (E12); order flips are noise (C21); restore never shortens (C22).

**Superseded rules.**

- D5 (period-k `+1 .. +(k_value - 1)` holes) -> superseded by D15 (backward only, every orbit).
- PK-3's WARNING (a missing `+j` appearance) -> superseded with D5: gone.
- C19's base-branch chain-order chords (2026-09-23) -> superseded by B7b (stable order everywhere).
- C9's second half (a pair no bridge spans is never a class) -> superseded by D11.
- "Flags never discard" -> narrowed by plan Q13 (a): an unresolved detection is not in N.
- Inertness as a fixed point over registered evidence (2026-09-16) -> transience (B5).
- Statuses "unique / trivial / truncated" of a walk search -> "resolved / ambiguous / unresolved".
- The "known asymmetry" of the stable interval (2026-09-15) -> both intervals (B2).
- The resonance-zone association of holes (July 2026) and the earlier partition rules (outermost
  beyond-cut hole, piece containment) -> holes know only their side of their own bridge.
- The initializer's `for orbit_index: for branch_index:` walk (before 2026-09), which visited
  `(0, 1)` before `(1, 0)` -> the physical chain of `FixedPoint.advance_key`.

**Audit divergences closed without a chalkboard row** (the 2026-10-09 fidelity audit,
`docs/code_fidelity_audit_2026-10-09.md`):

- the parity vote on near-coincident runs (B1);
- first-in-registry-order iterate links (B8);
- the stable interval left untested, NUM-18 (B2);
- the depth `m` set by `min_separation` and the zone filter (B3);
- inertness read off registered evidence (B5);
- the refined word assumed rather than derived (B7);
- base chords in chain order (B7b);
- walk bounds that claimed a word under a truncated budget (B6a);
- no free reduction of words (B6b).

## S. Stress-campaign findings that still bind

From `stress/numerics/report/report.pdf` (experiments E1-E12; the E12 re-verification is
`stress/numerics/e12_reverify.py`):

- **The resolution floor (NUM-3).** Deep blasts (from about 9 period-1 or 15 period-3 blasts) produce
  detections below the resolution ratio; they are held out of N and named with the remedy "refine".
  What C21 cannot prove (an image crossing 1e-7..1e-6 from a hole-bridge end with no registered
  preimage, k=3.2 from 9 blasts) is this floor too.
- **Image curves survive growth (NUM-12)**: blast / recompute / blast equals blasting straight (A.10).
- **The default `area_cutoff` is 1e-7** (E4: the cheapest cutoff converged on every case; at 1e-4 the
  period-3 crossing count still moves at fixed depth). It costs about 8-12x the points of 1e-4 per
  step; a recipe tuned at 1e-4 passes it explicitly. At 1e-5 the seed decreases on several series, so
  coarse blasting cutoffs are for exploration only (DEC-14's recommendation).
- **The point budget** (E3, E7, E11): past escape a step's point count can jump 5-8x in one refinement
  step, so the budget is checked inside the step; `grow_until_bounded` is the one stop condition
  shown to stop before take-off.
- **The GPU gains on the sweep and on a fused map, not on a cheap map** (E1, E2: routing Hénon
  gives 0.96-1.00x end to end; re-measured 2026-10-27: 0.98-1.00x, the map, refinement, merge, segment
  build and crossing resolution together 3.5% of the 1e-7 scale series, a cap of 1.037x). It stays
  opt-in. Since 2026-10-29 `enable_gpu` also builds and pre-tests the same-stability sweep's candidate
  pairs on the device, the host re-testing the survivors exactly (byte-identical); `route_map=False`
  keeps the map on the host (the 1e-7 scale series 50.1 -> 17.4 s, state byte-equal). An ODE map
  written as one `--fmad=false` kernel is bit-identical and fast (`duffing_map`). The two geometry kernels `curvature_area_batch` and `chord_deviation_batch`
  compute in their inputs' namespace (`numerics/namespace.py`), and the user's batched map calls run on the device
  under `enable_gpu`; merge, refinement's layer loop, the segment build and crossing resolution are
  host-only, each for the reason `gpu.py` gives.
- **Open for the author:** the blast-end same-stability count does not scale past about 10 blasts.

## O. Open items

- **DEC-25** (BCG-8): the inversion series loses classes at blast 6 by a genuine reference switch
  (2 -> 1 references); B2 gave it no new input. Options in `docs/open_issues.md`.
- **The pip-orbit stalls** (B3): the depth `m` stops growing when `f^{m+1}(q0)` is (ii) detected but
  unresolved (k28, k32, nested_outer, inversion), (iii) not detected (p3, heteroclinic), or (iv) on
  a bridge the zone filter never iterates. Options: refine, place by `image_cdist`, a finer image near
  the anchor.
- `p3_partitioned`'s zeroth-hole class (bridge `(4, 5)`) is UNRESOLVED and therefore active while its
  bridge is not mapped (outside the zeroth-hole theorem's hypothesis).
- `inert_classes` is not exported from `tanglepack.topology`.
- `scripts/henon_blast_period_3.py:95` is dead (its pinned pip is no candidate); reported, not fixed.
- `topology/diagnostics.py` is 988 lines, under its 990-line cap.
- The `p1(u) x p3(s)` heteroclinic links fail the beta and area laws (NUM-20, A.3).

## I. Index: every id -> section

Decisions (`.refactor/decisions.md`, `docs/open_issues.md` "Author decisions"):
D1 A.10 · D2 A.8 · D3 A.8c · D4 A.2 · D5 A.5 (superseded) · D6 A.6 · D7 A.5 · D8 A.10 · D9 A.8b ·
D10 A.4 · D11 A.8c · D12 (reversed 2026-10-31) A.4, A.10 · D13 A.8c · D14 A.8c · D15 A.6 · DEC-14 S · DEC-25 O, A.11.

Claude's calls under D8: C1 A.7 · C2 A.8 · C3 A.6 · C4 A.10 · C5 A.2 · C7 A.4 · C8 A.11 · C9 A.8c
· C10 A.8 · C11 A.8c · C13 A.8c · C14 A.8c · C15 A.7 · C16 A.1 · C17 A.8c · C19 A.8d (base half
superseded) · C20 A.10 · C21 A.8c · C22 A.10.

Issues (`docs/open_issues.md`): API-1 A.11 · API-2 A.11 · API-6 A.1 · BCG-2 A.8 · BCG-3 A.3 ·
BCG-4 A.5 · BCG-5 A.4 · BCG-7 A.8 · BCG-8 A.11, O · INT-1 A.8c · INT-2 A.6 · INT-4 A.6 · INT-7 A.2 ·
INT-8 A.8c · INT-9 A.8b · INT-10 A.8c · INV-2 A.5 · NEST-6 A.11 · NUM-1 A.6 · NUM-2 A.6 · NUM-3 A.2,
S · NUM-4 A.2 · NUM-5 A.6 · NUM-8 A.2 · NUM-9 A.1 · NUM-10 A.10 · NUM-12 A.10, S · NUM-13 A.1 ·
NUM-14 A.2 · NUM-15 A.1 · NUM-16 A.8b · NUM-18 A.5, H · NUM-19 A.2 · NUM-20 A.3, O · PERF-4 A.6 ·
PK-1 A.6 · PK-3 A.5 (superseded) · PK-4 A.8d · PK-6 A.4 · PK-9 A.6 · SIMP-6 A.8 · SIMP-8 A.10.

Plan questions (`docs/alignment_plan_2026-10-09.md` section E): Q1 A.5 · Q2 A.4 · Q3 A.8 · Q6 A.9 ·
Q8 A.10 · Q9 A.0, A.8c · Q11 H · Q13 A.2 · Q14 A.8d · N-a A.4 · N-b H.

Slices: B0-B12b H (B1 A.2, A.4 · B2 A.5 · B3 A.4 · B4a, B4b A.6 · B5 A.8, A.8b · B6a, B6b A.8c ·
B7, B7b A.8d · B8 A.3, A.11 · B9a A.11 · B9b A.7).

Work packages: WP0-WP12, WP8b, WP9b, WP9c, WP10b H. Stress experiments: E1, E2, E3, E4, E7, E11,
E12 S. F7 (a WP10 plan finding: candidates classified on first read) A.4.
