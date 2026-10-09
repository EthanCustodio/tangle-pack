# topology — knowledge file

Last updated 2026-10-02 (auditor, outcome of `hole-side-and-own-blast`; line numbers are on branch
`refactor/hole-side-and-own-blast` at `3b33baa`). Only the hole/partition path was studied in
depth; the rest of the layer is described in `CLAUDE.md`.

## Map (studied parts)

- `topology/StablePartition.py` — `punch_holes` :256 (direct holes, one per recorded pair except
  forward iterates >= k_value; side from crossing signs :349-356), `propagate_reference_holes`
  :378 (backward chain per reference, containing bridge by unstable-cdist span, carried hole point
  by the real inverse map), `_punch_in_bridge` :1488 (passes `span` to the side read :1546),
  `_image_arc_midpoint` :1562, side helpers `_arc_side_of` :726, `_bridge_side_of(trellis,
  bridge, point, span=None)` :789, `_hole_openings` :1148, invariants
  `check_bridge_rows_consistent` :1036 and `check_holes_share_bridge_side` :1124, combinatorial
  rows `row_of_end` :901.
- `Trellis.punch_holes` (Trellis.py ~1080-1150) wraps punch + propagate, then asserts I2 per hole
  bridge and I1 over all holes.

## How it works

- Hole side (`Hole.bridge_side`) = side of the hole's bridge, bridge oriented by increasing
  unstable cdist.
  - Direct holes (since 2026-10-02, author-approved): COMBINATORIAL, no geometry. `first` = lower
    unstable cdist end; `outward = +1 if second.stable_cdist > first.stable_cdist else -1`;
    `left iff first.crossing_sign * outward > 0`. Reason: the bridge leaves `first` along u+, the
    lobe's stable segment leaves `first` along s+ when it heads outward, and
    `crossing_sign = sign(cross(u+, s+))`. Openings are still `inward=True`, so this side only
    feeds I1 and plots. Reproduces the old measured side on every direct hole of k10, k28,
    k28 two blasts and p3 (test `test_direct_hole_side_is_the_side_of_its_coordinates`).
  - Propagated holes: measured from the backward-carried reference point against the IMAGE
    SUB-ARC only (`span` slice of the oriented polyline, one node of margin each side, by
    `bisect` over sorted node cdists — valid because node cdists are monotone along a bridge,
    probed by the reviewer). The side chooses the openings -> partition.
- `_arc_side_of(poly, p)`: nearest polyline VERTEX, tangent from its two neighbours, cross sign;
  returns None for fewer than 3 points.

## Invariants

- I1: holes of one origin share `bridge_side` (orientation-preserving map). I2: same row at both
  ends of a same-branch hole bridge. Both still asserted; the 2026-10-02 fix changed the
  measurement, not the invariant.

## Smells / open issues

- FIXED 2026-10-02: whole-bridge nearest vertex for propagated holes (p3 15 steps / 6 blasts,
  iterates -3/-5/-6 read "left" from a far fold) and the nearest-vertex misread of direct holes on
  tiny +1/+2 lobes. Guarded by `test_deep_p3_holes_share_bridge_side_through_symbolic_dynamics`
  [6_blasts, 15_steps]; removing `span` turns exactly those two red.
- CLOSED (WP4, INT-5): if a `span` lies strictly inside one polyline segment the slice is 2
  points and `_arc_side_of` returned None. Observed off the law cases (inversion x period 2 at
  10/10 steps: 19 of 26 propagated holes). The measured side is now read at the nearest
  SEGMENT, which answers on two points.
- CLOSED (WP4, INT-4): `_hole_openings` used the same whole-polyline nearest vertex for the
  stable nodes flanking each defining crossing; not latent, it opened nothing (or one half at
  both ends) on p3 at 6 blasts / 15-16 steps. The openings are sign arithmetic now.
- OPEN: p3 >= 6 blasts (or >= 15 steps) now runs but is unreliable. A new reference orbit
  appears (holes 12 -> 42 at 6 blasts), the partition is re-indexed so NONE of the 9 b0 class
  names survive, ~23-27 classes, class `t` unreachable (`R#2 -> R#3` on branch 0.0), a virtual
  `new1`, `PartitionFamily` warns "crossing abuts no hole". Words at b6, b7, b8 all differ.
  Verified by the 2026-10-02 audit probe. Needs the author's view.
  The figure maker's `p3_b6/cartoon.png` (branch) shows the shape: ~25 elements on every R row,
  1 (or 3 on 1.0) on the L rows, i.e. the new holes are all on the right side; the word for `c`
  changes letter each blast (`a w^-1 d^-1`, `a y^-1 d^-1`, `a uu^-1 d^-1`) with a long chain
  `d->e->f->...`. Same at `min_separation` 1e-4 and 1e-5.
- `numerics/ManifoldMachine.py:285-295,438`: non-contiguous `old_iterated_points` crash in
  `merge_manifolds` when a bridge's end nodes were iterated by neighbour bridges. No longer
  reached by `build_nested`; bug remains.
- CLOSED (WP5, PK-8): `SymbolicDynamics` warned "inert class u has the non-empty word 'v'" on
  period-3 runs (inert -> next branch's inert class); it now warns only for an inert word naming
  a non-inert symbol.

## Features

- p3 alone 0..5 blasts: words `a -> b -> c -> a u^-1 w^-1`, reliable, identical at every count
  (verified 2026-10-02 on the branch; pinned at 4 blasts by `test_p3_words_survive_four_blasts`,
  hole sides pinned by `test_p3_hole_sides_are_pinned`).
- p1 alone (probe recipe, `min_separation=1e-4`): b0 class unresolved (singleton), b1
  `a -> a u^-1 v^-1`, b2/b4/b8 identical (verified).

## Deferred / open

- CLOSED (WP3, 2026-10-08, C1): a stable branch with ONLY its anchor is one element
  `[anchor, end]` (`_build_intervals`, INFO on the `StablePartition` logger); regression
  `test_inversion_before_the_tangle_forms` (k10 inversion at 4 and 5 steps).
- FIXED 2026-10-07 (WP4 part A): tolerance-based ownership of registry cdists. Partition
  boundaries are registry floats, so crossings are owned by id / exact cdist, the iterated parent
  by exact containment (`PartitionFamily._parent_interval`) and `interval_at` is exact first. The
  E11 "zero-width iterated element" was a real sliver `(pip, outer]` of width 1.1e-6 beyond the
  pip whose midpoint fell within `cdist_tol` of the closed pip singleton (regression
  `tests/regression/test_deep_blast_sliver_beyond_pip.py`, builder `k28_deep_blast_sliver`:
  the analytic seed reproduces it, the `(4,-4)` seed does not). Still tolerant by design:
  `Trellis.image_of_element` and `DualWalk._end_tolerance` (estimated cdists). Open: a relative
  `cdist_tol` for deduplication (a numerics question).

## Open-issues round WP2 (2026-10-08, branch `oi/wp02-exact-iterates`)

- Scaling rule: a crossing is moved only by `Trellis.image_cdist` (table first, then
  `advance_key` / `per_step_beta`); a bare distance with no crossing behind it only by
  `FixedPoint.per_step_beta` (`StablePartition._backward_endpoint`, `land_element`'s unbounded
  end keep it by design). Nothing in topology reads an eigenvalue: `Trellis.lambda_u`,
  `Trellis.scale_cdist`, `TrellisBranch.eigenvalue` are gone, the orientation fallback reads
  `FixedPoint.preserves_orientation`.
- Pseudoneighbor (INV-2): a branch return is `k_value` map steps and laps a cdist by
  `beta ** k_value` (lambda_u squared under inversion, not lambda_u). r_{n+p} and every
  trajectory image are `Trellis.iterate` lookups; `_resolve_reference_end` /
  `_iterate_intersection` and every `match_rtol` are gone. A missing link means a missing image
  crossing; predicting it so no branch loses its holes is PK-3 / D5, WP5. An unlinked landing's
  stable cdist comes from `image_cdist`.
- StrongPip: the unlinked backward representative is two `image_cdist` calls (same numbers).
- PK-7 closed: every detected crossing carries both keys (`Tangle._key_of` asserts), so the
  "unknown unstable branch, test every residue" path and the StrongPip / BridgeClass None-key
  branches are deleted.
- Inversion law cases after WP2: holes 1 / 9 / 9 at 0 / 1 / 2 blasts (was 1 / 7 / 0).

## Open-issues round WP3 (2026-10-08, branch `oi/wp03-inversion`, base `c16e75d`)

- The anchor model (D4): a periodic point is ONE crossing, stored under its home pair
  `(u(o,0), s(o,0))` with that pair's sign, and the trellis puts it FIRST on every branch of
  its point (`Trellis.from_workbench`). Every bridge-end reader asks two helpers instead of
  `crossing.manifold_b_key`: `Trellis.unstable_key_of(pair)` (the arc's own unstable branch,
  read off the non-anchor end) and `Trellis.stable_key_at(id, unstable_key)` (C7: at an anchor
  an arc of `u(o,b)` meets `s(o,b)`, equivariant under the map, sign = the anchor's own since
  `cross(-u,-s) = cross(u,s)`). Users: `BridgeClass._element_at_end` / `_image_chain` (image
  branch = `advance_key(unstable_key_of(bid), 1)`), `PartitionFamily.owner_at_end` (used by
  `DualWalk.trellis_itinerary`, `SymbolicDynamics._member_pair`, `plotting._bridge_end_elements`),
  the iterated cut, `StablePartition._elements_at_bridge` / `_row_mismatch_message`,
  `Trellis.element_for`, `MinimalTrellis._bridge_sort_key`. `image_cdist(..., branch=)` maps an
  anchor to the next point's anchor on the advanced branch (`land_element`, `_land_singleton`,
  `image_of_element` pass it). `owner_of_intersection(..., branch_key=)` is needed at an anchor:
  without it the first result knowing the id answers.
- The arrangement: branch 1's arcs at an anchor fill `u-`/`s-` (`_add_edge`), so an inversion
  anchor is one degree-4 node (law `anchor_rays_are_its_branches`).
- Openings at the anchor by signs (D6): `Hole.openings` records are `(id, which, row,
  branch_key)`; at an anchor the hole opens the OUTWARD half of the stable branch on which the
  bridge's row is `bridge_side` (`_anchor_branch`: `s(o,b)` iff `(home_sign > 0) == (row ==
  "left")`, else `s(o,1-b)`; None when that branch is not computed). At the anchor the record's
  row is always `bridge_side`. An anchor end has no geometric row (no frame of its own). The
  geometry agrees on every case (law `anchor_opening_faces_its_hole`). Non-anchor ends: CLOSED
  by WP4, signs everywhere.
- Cross-side is flagged only for a pair on ONE branch (`SymbolicDynamics.translate_itinerary`).
  Elements are not glued across the point (C7; author may revisit, it would change `ElementRef`).
- Re-tagged, not fixed: `anchor_bridge_class_leads_its_tangle` x `inversion` (PK-3: the
  reference pair `(10,22)` has no `+1` link, so `s(0,1)` gets no hole and the anchor class is
  inert; passes at 1 blast) and `arrangement_euler` x `inversion_two_blasts` (INT-9: dense arcs
  between consecutive blast-image crossings no bridge spans).
- Inversion law cases after WP3: anchors 1, holes 1 / 9 / 9, classes 4 / 10 / 10 at
  0 / 1 / 2 blasts, no virtual class, `is_reliable` True throughout.

## Open-issues round WP4 (2026-10-08, branch `oi/wp04-sign-rule`, base `dad72d2`)

- Signs decide (D6, chalkboard stage 6). `_hole_openings(trellis, bridge, bridge_side)`: at a
  non-anchor bound `which = "outward" iff (crossing_sign > 0) == (bridge_side == "left")`, row
  `row_of_end`, branch the crossing's own; the anchor keeps `_anchor_branch`. Derivation: the
  bridge passes every bound along u+, and s+ is left of u+ iff `cross(u+, s+) > 0`. With
  alternating signs this IS the inward pair of a direct hole, so `inward=True` is gone; so are
  the `bridge_side is None` fallback and `SideSource "fallback"` (a propagated hole always
  carries a side). Deleted: `_stable_frame`, `_arc_side_of`, `_cross_sign`, `_row_polyline`,
  `_row_at`, `_approach`, `_row_mismatch_message` (folded into `bridge_row_violation`),
  `_bridge_end_geometry` (now `_approach_displacement`, displacement only),
  `_backward_endpoint`, `_bridge_unstable_span`. StablePartition.py 2223 -> 1963 lines.
- I2 is the sign form: `bridge_row_violation` fires exactly on EQUAL signs on one stable
  branch (C7 via `stable_key_at`); I1 (`check_holes_share_bridge_side`,
  `bridge_side_violations`) is the orbit check only. The geometry is a flag:
  `geometry_flags(trellis, bridges, holes) -> (row diagnoses with geometry_disagrees, side
  diagnoses with measured_disagrees)`, logged once per kind on the StablePartition logger by
  `Trellis.punch_holes(diagnose=True)` (silent on every law case and on p3_15/p3_16/p3_6_blasts;
  on `p3_coarse_blasts` 2 row flags (bridges (41,42) 2 nodes, (42,43) 4 nodes) and 20
  measured-side flags). The geometric row reader is `_anchorward_look` (bisect over node
  cdists) + `_approach_displacement`.
- Backward chain (PERF-4): `span = sorted(image_cdist(id, -step, "unstable", branch=start_key))`
  over the START bridge's two ids, container = bisect over
  `trellis.branch(advance_key(start_key, -step)).intersection_ids` (`orders` memo per
  propagation), termination `containing.id in visited`. Same containers as the old scan on every
  law case and deep run (probe `chain_ids.py`).
- INT-2 closed BY DESIGN: 600+ container lookups (law cases, deep p3, four deep E11 states) are
  all strict sub-arcs with at most one linked end; the other end is a crossing with
  `f^-1(W^s)`, beyond the trimmed stable manifold. The whole containing bridge gets the hole
  (chalkboard stage 5 says the same).
- Measured side at the nearest SEGMENT (foot point clipped), same rule as the law oracle
  `helpers.laws._point_side_of_bridge` (INT-5). Unchanged on every law case and deep run.
- C1 nit: the bare-branch element asserts no hole opens outward of its anchor.
- Behaviour: law cases unchanged (holes incl. measured sides, partitions, class counts, words:
  probes `holeset.py`, `partsig.py`). `inversion_period2`: 3 holes gain the tail opening the
  geometry dropped, partition unchanged. Deep p3 (15, 16 steps, 6 blasts): 81 -> 84 / 89 -> 92
  / 81 -> 84 elements, the unreachable class goes (unresolved 1 -> 0), the virtual `new1`
  stays (WP7). E11 p3 / nested_inner at 1e-5 run all 20 blasts (were AssertionError from
  blast 6); from blast 9 they hold 3 unresolved / 2 virtual (WP7) and the hole set dips at
  blast 7 (84 -> 64, BCG-5, WP8).
- Open after WP4: `OuterFaceError` on the 2-iteration E11 cells and k28 1e-5/1e-6 blasts 16-17
  (INT-10, WP10); inversion x period 2 at 10 steps still fails
  `backward_holes_keep_bridge_side` (the law re-measures a carried point that escapes a
  3-node reference lobe; a law-oracle question), `bridge_endpoints_on_own_branch` (WP7-D1) and
  the PK-1 symbolic laws (WP7). Hole PLACEMENT (`_stable_arc_midpoint`, `_nearest_arc_point`,
  `_step_into`, `_image_arc_midpoint`) stays geometric: it places plotted coordinates and the
  carried point, nothing decides on it.

## Open-issues round WP5 (2026-10-08, branch `oi/wp05-period-k-holes`, base `3e12d53`)

- Period-k stagger (PK-3 root cause). With WP2's exact table every missing `+j` link is an
  image crossing NOT COMPUTED: a period-k chain is grown one map step per branch, so each
  unstable branch reaches beta further than the one before it in the cycle, and a reference on
  the outermost lap of the longest branch has no forward image. Landings (`image_cdist`,
  unstable) vs the target branch's last crossing: inversion `+1` 55.5 vs 13.2 on u(0.0);
  inversion_period2 `+1..+3` 1.13 / 1.55 / 2.11 vs 0.38 / 0.52 / 0.71; p3_6_blasts `+1` 4097 vs
  1662, `+2` 6355 vs 2579. Every gap on every builder is "past its last crossing".
- `extend_pseudoneighbor_trajectories` walks `for steps in count(1)`; a FORWARD walk that stops at
  `steps < k_value` records the branches `advance_key(ref.branch_key, j)`, `j = steps .. k-1`
  (the `key in seen` break is not a gap). After the loop: ONE WARNING per stable branch on
  `tanglepack.topology.Pseudoneighbor`, naming each reference that misses it and, per unlinked
  member, its predicted landing against the unstable branch's last crossing (helpers
  `_missing_appearance`, `_branch_text` -> `A:s(2.0)`; StablePartition imports `_branch_text`).
  Counts: inversion 1, inversion_one_blast 1, inversion_period2 3, p3_15 / p3_6_blasts /
  p3_coarse_blasts 2; 0 on k10, k28, p3, nested, nested_inner_blasts, inversion_two_blasts.
- Orphan window pairs (why the gap is not filled by walking each branch's window): the gap
  branches do have consecutive "pairs" in their windows (p3_6_blasts s(2.0) `(15,23)`, s(0.0)
  `(18,20)`; inversion s(0.1) `(17,3)`), but each one's open STABLE interval holds the
  predicted images of the tip crossings (s(2.0) `(0.4966, 0.9546)` holds 0.6275, 0.661,
  0.7466), so they are no pairs by the definition; they pass only because the stable-side check
  reads registered crossings alone (NUM-18, deferred). Filling is DEC-19 (WP5 slice 2; option
  (b), an unstable-side cut of N, changes the deep references, BCG-4).
- `_strong_pip_cuts(trellis) -> {stable key: cdist}` from `image_cdist(pip, j, "stable")`,
  `j < k_value`: one cut per branch always; a cut not from the table is flagged in one WARNING
  (PK-6). Never fires on a builder (cuts complete everywhere). `_cut_ids` and both
  `_pair_fixed_point` copies are deleted.
- `_is_forward_beyond_fundamental(pair)`: `iterate is not None and iterate >= k_value` (D5, kept;
  every forward pair carries `branch_key`). `punch_holes` partitions `forward` / `pairs` by two
  comprehensions; the no-spanning-bridge WARNING names each skipped pair, iterate, origin and
  branch (INT-9 stays WP7's).
- PK-9: `propagate_reference_holes` has no `max_steps`; it loops `count(1)` and ends because
  every step stops or adds a new bridge id to `visited` (finite). Longest chain seen: 21 steps.
- PK-8: `_edges` warns only when an inert class's word holds a non-inert symbol.
- Symmetric per-branch hole counts are NOT a WP5 property: p3 4/4/4 at 0-5 blasts, 16/14/12 at
  6; p3_16 (no gap) 14/16/16 because the backward chain has 20 steps per origin, not a multiple
  of k = 3 (BCG-4, WP8).
- Inversion x period 2 at 10 steps is not a WP5 gap: 7 laws fail on main; `partition_covers_branch`
  is an element of s(0.0) opened at both ends by two holes of one origin a lap apart (iterates
  -3 and -7), the -3 one on a cross-branch container (INT-2 convention vs the law's reading);
  the landing failures are outermost elements whose image runs past the next branch (INT-8).
  Owners WP7 / WP8. Not a law case.
- Behaviour: law cases unchanged (holes with openings and sides, partitions, class counts,
  words; probes `wp04/holeset.py`, `wp04/partsig.py`, `wp05/impl/holeopen.py`).

## Open-issues round WP6 (2026-10-08, branch `oi/wp06-minimal-alphabet`, base `b90b718`)

- Evidence vs alphabet (D2). Activity is a fixed point over the images of ALL bridges
  (`_resolve_inertness`), and the minimal trellis maps only active hole bridges, so the alphabet
  cannot feed the minimal trellis. The order is: `bridge_classes(trellis, partitions)` = the
  EVIDENCE table (every bridge) -> `minimal_trellis(trellis, holes, evidence)` (kept as
  `MinimalTrellis.classes`) -> `MinimalTrellis.alphabet` = `classes.restricted_to(kept_bridge_ids)`
  (same entry objects, table order, members whole: a class's other members are its image
  evidence). Hole bridges are never loops (the hole's empty stretch separates their ends), and an
  image bridge's image-chain parent is its hole bridge, so a kept bridge is never an orphan.
- Forward folding (`_non_loop_parent`): every bridge's `_image_chain` is computed once per
  `bridge_classes` call (also the inertness evidence: `_collect_image_evidence(chains, ...)`), each
  image pair indexed by the one bridge whose chain holds it (asserted), and a loop climbs those
  parents through loops to the first non-loop bridge. Probe `wp06/fold_probe.py` over the 10 law
  cases, the deep p3 runs and every series step: 0 double parents, 0 disagreements with the old
  backward walk where it found an ancestor, 0 orphans. Inversion at 0 blasts had the only
  standalone loops left after WP2 (four, under the anchor bridge `(0,7)` whose `-1` link breaks the
  backward walk); they fold into one class (evidence 4 -> 2). `BridgeIterator.preimage_bridges`
  as a fallback is redundant with a partial-bijection table (it finds a parent only when both ends
  have `-1` links, and then the forward index already holds the pair). Orphans: one WARNING,
  `BridgeClassTable.orphans`, no entry, never lettered. `_resolve_inertness` starts from nothing
  (no loop-class seed).
- The evidence-letters leak: `registered_evidence` translated with the SHARED `letters`, so an
  image naming a class outside the table became `new<n>` in the driver's `virtual_classes` and a
  matrix node (k28 series blasts 6-7: 2 virtual and a 6-node matrix with the restriction alone).
  It translates on a copy now; `virtual_classes` comes from words only.
- Letters (SIMP-6): `SymbolicDynamics.letter(index, inert=False)`, one bijective numeral over
  `a..t` / `u..z`; `inert_letters(table)` and `_active_letters` are one comprehension each;
  `symbolic_dynamics` asserts one letter per class. Exported from `tanglepack.topology`.
- Measured (series, alphabet size): k28 3,3,4,4,4,4,4,4,15 (evidence 3,3,4,4,4,4,6,7,23); p3
  9 x6,23,24,25,25; inversion 1,9,9. E11 small (`wp06/e11_proto.py`): p3 alphabet 25 at blasts
  8-16 (evidence 27 -> 35), matrix 23 / 1 virtual over 9-16; k28 4 (2-7), 15 (8-10), 19 (11-12).
- OPEN for WP7 (PK-1 / D3): the alphabet is not closed under the word map. An image-bridge class
  that is not a hole-bridge class can name the class of its image's image, which the minimal
  trellis does not keep (p3 from blast 7: `u = {(40,41), (44,45)}` -> `w = {(46,47), (50,51)}`;
  k28 blast 8 `e` -> `g`); those symbols are virtual now (p3 deep 1, k28 blast 8 2, was 3).
- Inversion (0 blasts) has no active class: its alphabet is the inert class of the hole bridge
  `(22,10)`; the anchor bridges carry no hole because s(0.1) has no window hole (PK-3 / DEC-19).

## Open-issues round WP7 (2026-10-08, branch `oi/wp07-walk-exists`, base `e883720`)

- Closure (C9, PK-1). Every walk already existed on the deep runs after WP4; the virtual `new1`
  was the class of the IMAGE of an image-bridge class (`p3_6_blasts`: `t -> aa u^-1 ab u aa^-1`,
  `aa = {R(0.0)#1, R(1.0)#27}` = the images `(40,41)`, `(44,45)` of hole bridge `(38,39)`, whose
  trivial walk names `{R(1.0)#1, R(2.0)#25}`). A walk pair is the class of a real piece of
  `f(B)`, so `symbolic_dynamics` closes the alphabet under the word map OVER REGISTERED
  BRIDGES: the driver iterates `entries` while appending; a named class that
  `MinimalTrellis.classes` holds copies its entry (`replace(..., letter=None)`, members and
  inertness kept, so its members verify it); active ones are lettered by `letter_for` (the
  session passes `BridgeAlphabet.letter_for`), else the next unused letter of their series;
  `ClassDynamics.named_by`. A walk pair NO bridge spans is a construction error (D3, review of
  the first WP7 version, which lettered such bare pairs: 82 of 90 closure classes on the deep
  runs had no bridge): `_unspanned` unresolves the class whose itinerary names it (itinerary and
  source None, `search` kept, reason "grow or blast"), ONE WARNING naming the pairs and the
  non-loop members; nothing is added from that itinerary; `is_reliable` False (SUPERSEDED by D11, WP8b: such a pair is a flagged class). DEC-21 asked
  whether such a pair must get its bridge. Translation runs once the table is closed
  (`SymbolicDynamics.table` = seed then closure in discovery order); `session.bridge_classes()`
  stays the seed. Measured (`wp07fix/probe.py`): p3_15 / p3_6_blasts / nested 6 inner blasts:
  no closure, 1 unresolved (`aa`, the image-bridge class); p3_16 27 -> 30 (`ab`, `ac` verified,
  `ad` unresolved); nested 7 inner blasts 28 -> 29 (1 unresolved); k28 8 blasts 15 -> 20 (2
  unresolved). Law cases: no closure class; words / matrices / letters / reliability identical.
  Registry crossing ids of one build differ between processes (same element pairs).
- Virtual classes deleted: `SymbolKind` is `active | inert`; a class outside the closed table
  (only evidence can name one) is named by its label and counted active; `transition_graph`
  asserts every edge target is a node; `plot_transition_graph` has no dashed sinks.
- `EvidenceDiagnosis` (INT-1): `ImageEvidence.chain` / `.diagnosis`; cause `"stable end"` = the
  first run of >= 2 consecutive chain crossings with equal `(manifold_b_key, crossing_sign)`
  (`crossings` = that run), else `"walk"`. k28 at 6 blasts: `a` member `(44,45)`, run
  `(74,75,76)` on s(0.0), ratios 30-8400 (not under-resolution). One WARNING per disagreeing
  member (the old length-mismatch WARNING is folded in).
- `is_reliable` also needs every non-None landing `contained` (INT-8; no law case straddles;
  `p3_coarse_blasts` does).
- INT-9 / NUM-16: `Trellis.is_computed_arc(lo, hi)` = a Bridge spans the pair, or both
  crossings' `unstable_manifold` is one object (None = an anchor, on every curve). Dense
  `Arrangement._build_unstable_arcs` skips other pairs (DEBUG lists them); their slots get the
  usual stubs. Gaps per law case: k28_one 1, k28_two 2, nested_inner_blasts 4, inversion_one 2,
  inversion_two 5. Sparse mode unchanged. The two "no spanning bridge" WARNINGs
  (`punch_holes`, minimal trellis) name the cause (D9: the image of a bridge no blast iterated).
- PK-5 / DEC-20 (the WARNING below was DELETED in WP8b: the trellis is cut at the pip orbit, D10): `DualGraph._check_where_f_of_s_stops` WARNs per stable branch of a pip's fixed
  point whose unified span is not `S \ f(S) = (c_s(f(end of previous branch)), end]`
  (`image_cdist(last, 1, "stable", branch=previous)`, tol `cdist_tol` from the table, else
  `SCALING_RTOL * end`). Equal on every law case and deep run (all from the table, exact);
  inversion x period 2 at 10/10 steps warns on s(1.0) (`(0.577, 4.413]` vs `(0.504, 1.749]`),
  the one remaining unreachable walk (class `n`).
- PK-4 (slice 2, author): crossing 44 on `p3_6_blasts` is a foreign-branch chord of the lobe
  excursion 41-44 across s(1.0); by Jordan alternation the empty stretch is `(41,44)` (stable
  order), not the base-branch unstable-order pairing. Moves refinement labels only.
- Residual laws on the deep runs (also on base): `partition_covers_branch` (BCG-4, WP8),
  `dual_face_side_is_geometry` (INT-10, WP10). Inversion x period 2 at 10 steps also fails
  `partition_covers_branch`, `bridge_endpoints_on_own_branch`, `backward_holes_keep_bridge_side`
  (WP8 / WP10).
- Deferred: escape-lobe images (WP9 / WP10); DEC-21 (a walk pair no bridge spans; answered D11, WP8b).

## Open-issues round WP8 (2026-10-08, branch `oi/wp08-recompute-identity`)

- Pseudoneighbor endpoint frame (BCG-5): in `_unstable_interval_is_clear`, after the linked-landing
  branch, an unlinked candidate r is compared with the pair's images `iterate(x0|x1, -steps)` on
  r's own unstable branch when at least one is registered (the other by `image_cdist`); r punctures
  iff its registered cdist lies strictly between them (r itself being an image: skip). Scaling +
  2D collision only when neither image is registered. Root cause of nested 10 outer blasts
  (`(37,38)` lost to 592, 6178 below `M^8(37)` yet 1.3e-5 inside after scaling, noise 0.0103) and
  k32 12 blasts. Law cases / deep runs / series: hole sets identical.
- The decreases left were DEC-22 (answered D12; root cause and fix in WP8b below), NOT shown definitional: the witnesses (nested 11 outer
  blasts, k32 11 and 13) sit within 1e-6..1.6e-5 of a window crossing, unlinked both ways, INSIDE
  the placement error bar `min(l_u, l_s) / ratio` (the ratio is the crossing angle's resolution,
  not cdist placement). Only k32 blast 10 is a clear lambda-lemma copy (1.25e-2 inside).
- Element identity (C10): `PartitionInterval.identity = (branch_key, side, lo_id, hi_id)`; registry
  ids persist across one session's recomputes. `ElementRef` stays positional (an identity inside its
  equality broke 58 unit tests built on hand-made refs). Names stay positional.
- Open: `partition_covers_branch` fails on `p3_6_blasts` / `p3_15` / `p3_16` (open elements that are
  no hole's stretch: 6 / 6 / 5 spans, e.g. `(0,35)` on p3_6_blasts) -- hole bridges of one deep chain
  nest on a branch and the law assumes disjoint hole stretches; unexamined (not a law case).

## Open-issues round WP8b (2026-10-08, branch `oi/wp08b-d10-d12`, base `014ea4b`)

- The trellis of the pip (D10 + D12): `Trellis.cut_at_pip_orbit(pips)` cuts N to `T(q0) = (W^u[P,
  f^m(q0)], W^s[P, q0])`: stable branch `j` at `c_s(f^j q0)` (`image_cdist(pip, j, "stable")`, `j <
  k`), the unstable branch of `f^i q0` at `c_u(f^i q0)` for `i = m-k+1..m` (`m` = the largest `n`
  with `iterate(pip, n)` registered: `forward_depth` then `+1` composition). Exact comparison
  against a registered cut, `+ cdist_tol` against a PLACED one (one WARNING on the Trellis logger,
  PK-6 moved there). `Trellis.__init__` keeps `_computed` (branch orderings, bridges); the
  `strong_pip` PROPERTY setter cuts, None restores; `classify_strong_pips` sets None FIRST, so
  candidates are always classified on the computed trellis; `clear_results` restores. Same cut
  twice is a no-op; stage outputs are never cleared by the cut.
- Joined arcs: `_bridges_inside` keeps a bridge whose present ends are kept and, for a run of
  consecutive bridges of one curve through dropped crossings between two kept ones, ONE
  `Bridge.joined(run)` (root of the first, tail of the last; the parent curve's points, no map
  call), memoised in `_joined` by id. A run never crosses a gap (bridges exist only within one
  curve). Joined arcs are NOT workbench bridges: `Topology._zone_key` and
  `TangleSession.cartoon_enclosures` resolve through `trellis.bridge_between`. Without them
  `arrangement_regions_*` failed on inversion x period 2 at 10 steps (planner's prototype).
- Membership vs raw lookups: `Trellis.__contains__` = a crossing of the cut N. `iterate` and
  `image_cdist` stay raw (X). Consumers that need an iterate AS A CROSSING: pseudoneighbor
  trajectory extension and its gap note, `BridgeClass._image_chain` (None = no evidence),
  `DualWalk._land_singleton` / `_debug_endpoint_owners`, `Arrangement.image_of`. The
  pseudoneighbor all-iterates test keeps raw landings (`_unstable_interval_is_clear`).
- `Pseudoneighbor._strong_pip_cuts` / `_within_cut` deleted (check set = `own_intersection_ids`).
  `DualGraph._check_where_f_of_s_stops` deleted (vacuous: the unified window is `S \ f(S)`).
  `Topology.__init__`: stages first, then the selection trellis, which `TangleSession.trellis`
  cuts at every selected pip on each read (review fix: no longer cut by `Topology`).
- Why D12 is the BCG-5 fix: every decrease happened when N held crossings deeper than the
  deepest registered pip-orbit point (the pip lap stalls when its bridge falls under the blast's
  `min_separation`: k32 at blast 9 / `f^10(q0)`, nested_outer at blast 10). Measured
  (`wp08b_impl/replay.py`): k32 34 classes at blasts 8..13, nested_outer 13 at 2..11, inversion
  series 6, 7, 8, 9, 9, 9. Only the three inversion law cases change (1 -> 6, 9 -> 7, 9 -> 8
  classes); everything else bit-identical (id-free describe). N freezes once the pip orbit stops
  (k32 N = 521 from blast 9 while the registry grows to 4366 at 13).
- D11: `symbolic_dynamics.add` builds a member-less `BridgeClassEntry(named, tangle=...)` for a
  pair outside the evidence table; `BridgeClassEntry.spanned = bool(members)`; describe's class
  line ends "; no computed bridge yet", the closure line counts them; `_unspanned` deleted. A
  member-less class landing on a singleton has no member chain and stays unresolved (not seen on
  any run). Deep runs: p3_15 / p3_6_blasts 23 -> 43 (20 flagged), p3_16 27 -> 49, nested 6 / 7
  inner 27 -> 47 / 28 -> 49, all reliable; k28 at 8 blasts 15 -> 23, unreliable (ambiguous `q`,
  `aa`, evidence `r`: WP10); inversion x period 2 at 10 steps 14 seed / 17 closed, reliable.
- k10 candidates under the cut: [1, 2, 3, 4, 7] -> 1 fine; 2 raises the I2 row assertion (the
  joined anchor bridge `(0, 2)` wraps the end of the cut branch; a numerical trim at 2 raises the
  same); 3 has 4 holes / 3 classes, not reliable; 4 and 7 have no hole. The cut equals a trim at
  every candidate (`test_pip_orbit_cut`).
- Open / for later: k32 dynamics is unreliable from blast 6 of its series (not examined; WP10's
  evidence/ambiguity work); PK-4 untouched; re-classification is O(N^2) per rebuild (PERF).


## Open-issues round WP9 (2026-10-08, branch `oi/wp09-heteroclinic`, base `b460363`)

- Selection rule (`Trellis.from_workbench`): a non-anchor crossing is bucketed only when BOTH its
  keys are branches of the selection. A tangle's own trellis therefore holds no heteroclinic
  crossing; `Trellis.__init__` passes the given bridges through `_bridges_inside(bridges,
  crossings)` so an arc ending at a foreign crossing is JOINED across it (same rule and memo as a
  crossing past a pip cut; `_computed` stores the joined list). No-op on every trellis without
  heteroclinic crossings, the minimal trellis's synthetic one included. A run through foreign
  crossings that ends in a partial tip is dropped, like a run past a cut (only an uncut trellis
  could notice; none of the cases do).
- `Trellis.pips: {fixed point: pip}` is set by `cut_at_pip_orbit` (the `strong_pip` setter too).
  Cut at both pips the selection is `T(q0_A) ∪ T(q0_B)`: `cut_at_pip_orbit` already merged the
  cuts per branch key, so D10/D12 needed no change for two fixed points.
- D7 pair test (`Pseudoneighbor._unstable_interval_is_clear`): the cycle, `k` and β are the ARC's
  fixed point's (`u_key[0]`), so a heteroclinic pair on B's stable branch is tested against the
  crossings whose unstable side is A's with A's cycle; the "missing from the branch cycle"
  WARNING is gone (it cannot happen). The scaled-fallback 2D endpoint collision applies only when
  the candidate's stable fixed point is the pair's (stable cdists of two manifolds never compare);
  the own-frame branch is unchanged. `_reference_branch` reads `trellis.pips.get(fp)`.
  `extend_pseudoneighbor_trajectories`, `punch_holes` and `propagate_reference_holes` needed no
  change: the window's `k` is the pair's STABLE fixed point's (D5), the backward chain follows the
  bridge's own `advance_key` (A's cycle), `_containing_bridge` reads A's branch of the selection.
- `heteroclinic` (p3 16u/10s): refs 3 -> 4 (one heteroclinic reference on p1's window), holes
  50 -> 77 (its backward chain, ~23 steps along A's cycle, holes p3's anchor bridges), alphabet
  31 -> 57, closed 53 -> 81, connecting classes 0 -> 2, one sparse-arrangement component.
  DEC-23: with the lobe's hole, two equal-length shortest walks pass either side of it on p1's
  window (the faces on each side are one face: the lobe and the bubble across it are bays of the
  stable manifold); three classes are ambiguous and the dynamics unreliable (strict xfail
  `every_active_class_walks` x `heteroclinic`, WP9 slice 2). The member images disagree between
  classes, so no evidence rule breaks the tie for all of them.
- Open for later: `partition_covers_branch` on every p3 16u build (WP10, inherited);
  `dual_face_side_is_geometry` on p3 16u/9s and 16u/12s (half-plane test against a non-convex
  face's representative point; 16/10 avoids it, not fixed); escape-lobe images (the four
  heteroclinic crossings carry no iterate links, WP10).

## Open-issues round WP10 slice 1 (2026-10-08, branch `oi/wp10-numerics-leftovers`, base `3f27456`)

- `Trellis.strong_pip_candidates` is a property (`_candidates`, None until read): classified with
  the default slack over the COMPUTED crossings (`_computed[0]`, never the cut set; is_strong_pip
  reads the raw registry/iterate so the cut does not matter for the test itself).
  `classify_strong_pips()` keeps its eager behaviour; `clear_results` resets to None.
- `DualGraph._outer_face(members)`: easternmost vertex of the members' traversal rings
  (`_traversal_ring`: arcs stitched, consecutive duplicates and the closing vertex dropped), the
  face whose sector there (clockwise from `next - v` to `prev - v`; the whole circle at a stub tip)
  holds due east; asserted unique. `OuterFaceError` and the area-sign rule are deleted;
  `_traversal_area` survives only for `_innermost_container`'s magnitude.
- `SymbolicDynamics._compare_evidence`: `"stable end"` evidence -> `word_agrees` None + INFO;
  `"walk"` -> False + WARNING; `verified` over the evidence with a verdict (C11).
- Laws (tests/helpers/laws.py): `check_partition_covers_branch` last clause: an end is open iff a
  hole opens that half there, the halves derived by `_sign_openings` from the hole's bounds,
  `bridge_side` and crossing signs (never `hole.openings`, which the partition is built from;
  pinched singletons skipped);
  `check_dual_face_side_is_geometry` = a local probe off the stable edge's middle segment pushed by
  a quarter of its distance to every OTHER boundary segment of the region (its own arc included:
  `inversion_period2_10` has a stable hairpin region `(3, 4)` where the other leg is the nearest
  wall). `_anchorward_tangent` is gone.
