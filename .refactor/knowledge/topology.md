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
- Latent: if a `span` lies strictly inside one polyline segment the slice is 2 points,
  `_arc_side_of` returns None and the hole gets `bridge_side=None` with a misleading DEBUG.
  Never observed (175/175 span reads had >= 3 points). Widening to >= 3 points would close it.
- Latent (believed): `_hole_openings` :1148 uses the same whole-polyline nearest vertex for the
  stable nodes flanking each defining crossing.
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
- `SymbolicDynamics` warns "inert class u has the non-empty word 'v'" on period-3 runs (inert ->
  next branch's inert class); expected, noisy.

## Features

- p3 alone 0..5 blasts: words `a -> b -> c -> a u^-1 w^-1`, reliable, identical at every count
  (verified 2026-10-02 on the branch; pinned at 4 blasts by `test_p3_words_survive_four_blasts`,
  hole sides pinned by `test_p3_hole_sides_are_pinned`).
- p1 alone (probe recipe, `min_separation=1e-4`): b0 class unresolved (singleton), b1
  `a -> a u^-1 v^-1`, b2/b4/b8 identical (verified).

## Deferred / open

- DEFERRED (author, 2026-10-07; stress campaign WP4 part B): a stable branch with ONLY anchors.
  The k=10 inversion saddle grown 4 or 5 unstable steps (E7) registers no crossing past the
  anchors of a stable branch. `_build_intervals` then gets a single boundary (the anchor is also
  the outermost crossing) and yields zero intervals, and `_element_of_intersection` asserts
  "crossing 0 at stable cdist 0.0 ... is owned by 0 partition elements" (reproduced on `main` at
  `8148d48` via `cases.build_inversion(unstable_steps=4|5)`). Neither check is over-strict; what
  such a branch partitions into is undecided. Candidates: (i) one closed singleton at the anchor
  (consistent with the partition spanning anchor .. outermost crossing, but a singleton owns no
  stable edge and no dual-graph node, so the branch drops out of the dual graph); (ii) one closed
  element `[anchor, end]` to the branch's end, like an unbounded high end in `land_element`, which
  gives the branch dual-graph nodes with no crossing to own. Either must work with ONE anchor per
  branch: the inversion saddle's duplicate (0,0) anchors are the strict xfail
  `one_anchor_per_unstable_branch` (`tests/cases.py` `KNOWN_ISSUES`) and must stay failing until
  that is fixed. Evidence: the planning probe monkeypatched (i) and ran the full stack (4 steps:
  0 classes; 5 steps: 1 unresolved class, provisional). When decided, log at INFO on the
  `StablePartition` logger ("no crossing past its anchor(s); grow the unstable manifold"), never
  on `PartitionFamily` (whose WARNING count `test_iterated_cut.py` watches), and add the
  regression `test_inversion_before_the_tangle_forms` (4 and 5 steps run through
  `symbolic_dynamics`). Recorded in the `StablePartition` Dev Notes too.
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
