# topology — knowledge file

(2026-10-15 B12 analyst: see "B12 analyst" at the end -- section A vs HEAD, citations, `__init__` trim facts. Earlier 2026-10-15 addendum: B9b `one-element-at-end` OUTCOME, branch `refactor/one-element-at-end` @ `94fca6b` off `5f0db5f` (not merged when written): ONE reader of a bridge end's element remains, `PartitionFamily.owner_at_end`; `Trellis.element_for`, `StablePartition._elements_at_bridge` + the `StablePartitionResult.elements_at_bridge` field, and `BridgeClass._element_at_end` are deleted -- see "B9b outcome" at the very end, which supersedes the "B9 analyst" / "B9b analyst" sections' line numbers. Before that: B9b ANALYST at `5f0db5f`: the four element-at-end readers re-derived, the probe re-run, the import-cycle options measured. Earlier: B9a OUTCOME, branch `refactor/session-three-concerns` @ `5f0db5f`, since merged: `plotting.py` is 3,479 lines there and gains `cartoon_enclosures(trellis, bridge_ids, parents, fixed_points)` (moved from the session; reads only the trellis's registry / `manifolds` / `bridge_between` and the fixed points passed in; not in `topology.__all__`); its module docstring no longer says the session keeps delegates; `cartoon_zones` is now `tanglepack.loom.cartoon_zones`. The element-at-end functions are untouched (B9b). Earlier addendum: B9 `session-three-concerns` ANALYST at `HEAD` `517aa01` -- the "element at a bridge end" map and agreement probe in the section "B9 analyst: one element at a bridge end" at the very end.)

Last updated 2026-10-15 (B6b `word-finished` OUTCOME audit at the very end: branch `refactor/word-finished` @ `65edec4`, base `63477e3`, NOT merged -- see "B6b outcome"; it supersedes the WalkStatus / truncated / evidence-location statements of every earlier section (e.g. the INT-1 / C11 / B6a bullets that name `SymbolicDynamics._compare_evidence` or `search.truncated`). Before that: B6b ANALYST section, HEAD `63477e3` = B7b merged: the symbolic-dynamics walk/evidence map and the free-reduction probe). Before that: 2026-10-14 (B7 OUTCOME audit at the end: branch `refactor/refined-level-defined` @ `a6a4f68`, base `3fb310a`, NOT merged -- see "B7 outcome"; B7 analyst section before it; B5 is MERGED (main `3fb310a`): the "Inertness" map below is main; B5 analyst sections at the end; top sections still the B4b outcome audit: branch `refactor/holes-are-combinatorial` @ `ba23a1f`,
base `f166573`; NOT merged -- the Map below is the BRANCH, `main` still has `Hole.coords` and the
placement in StablePartition). The top sections describe the branch; the dated sections below are
each WP's notes, kept as the record. The whole layer is described in `CLAUDE.md` and
`docs/pipeline_chalkboard.md`.

## Map (the hole / partition path, at `ba23a1f`)

- `topology/StablePartition.py` (1,091 lines; 1,354 at B4a, 1,947 before B4a) -- the ALGORITHM, no
  numpy, no map, no plane: `bridge_for_pair` :244, `punch_holes` :268 (REFERENCE pairs only;
  `assert all(p.iterate >= 0 ...)` :304; forward pairs skipped, D15; side from crossing signs),
  `propagate_reference_holes` :376 (backward chain from each reference that has a hole; span =
  sorted `image_cdist(...)[1]` of the reference bridge's ids, the `from_table` flag kept; NUM-5
  WARNING once per branch when a SCALED end is within `SCALING_RTOL*|x|` of a crossing of the
  container's branch :449-466; side CARRIED by `side_after`), `partition_stable_manifold` :482,
  `row_of_end` :545, `rows_of_bridge` :599, `check_bridge_rows_consistent` :621 (I2) and
  `check_holes_share_bridge_side` :641 (I1) -- wrappers over cores still in diagnostics,
  `_hole_openings` :668 (D6 on every end; at an anchor the branch is `Trellis.stable_key_at`, and an
  "anchorward" verdict opens the OUTWARD half of `s(o,1-b)` -- one curve through the point; nothing
  when that branch is not computed), `_near_far` :715 (orders `bounding_ids`, a hole field: stays),
  `_containing_bridge` :750, `_punch_in_bridge` :774 (ids + side + openings -> `Hole`, never None),
  `intervals_from_marks` :886, `owns_cdist` :972. Gone in B4b: `_map_backward`, `_region_key`,
  `seen`, `punched_regions`, the fallback start, the 4 placement helpers, `epsilon`.
- `topology/plotting.py` (3,408 lines) -- `plot_holes` :238 reads `origin`, `iterate` and
  `hole_point`; `hole_point(trellis, hole, epsilon=0.05)` :312 is THE placement (a view): reference =
  `_step_into(stable-arc middle, nearest bridge node)`; propagated = image sub-arc middle nudged
  toward the reference point carried back `-iterate` steps by `dynamical_system.map_inv` (no map:
  nudged toward the stable-arc middle); `epsilon` now scales BOTH nudges (base used a literal 0.05
  for the propagated one). A propagated hole whose origin has no reference hole in `trellis.holes`
  raises a bare `StopIteration` (`next()` without default, :343). Helpers `_step_into` :375,
  `_stable_arc_midpoint` :396, `_nearest_arc_point` :454, `_image_arc_midpoint` :473 (verbatim moves).
- `topology/diagnostics.py` (713 lines) -- flags, never decisions; own logger. The 16 B4a names,
  `side_source(hole)`, `measured_sides(trellis, holes)` :429 (imports `plotting.hole_point` inside the
  function for the reference point, then carries it with `map_inv`; runs only on
  `punch_holes(diagnose=True)`, which is the DEFAULT, so `topology()` still calls the map there, on the
  flag path).
- `Hole` (TopologyResults) = `near_intersection_id, pair, bridge_side, bounding_ids, openings,
  iterate, origin`; no `coords` / `measured_side` / `side_source`.
- `Trellis.punch_holes(pairs=None, *, propagate=True, diagnose=True, verbose=False)` (no `epsilon`);
  imports the checks / `punch_holes` / `propagate_reference_holes` at CALL time (spies in
  `test_production_checks.py`); asserts I2 per hole bridge, I1 over all holes; `diagnose=True` calls
  `warn_geometry_flags`. Stale comment :1250-1251 (and StablePartition :306-308) still says
  propagation "seeds already-punched orbits and regions" from `pair.hole`.
- Laws `check_holes_are_classified`, `check_backward_holes_keep_bridge_side`,
  `check_direct_hole_side_is_coordinate_side` (tests/helpers/laws.py) read their point from
  `plotting.hole_point`: the geometric oracle and the figure share one code path.

## Map (inertness, at the B5 branch `3fb310a`; `main` = `7d3a294` still has the evidence fixed point)

- Inert = TRANSIENT in the class word graph (Q3, I-a): a class is inert iff its word is empty or
  names only inert classes; an unresolved class (no itinerary) is active, and so is every class
  whose word reaches one. `SymbolicDynamics.inert_classes(words)` :499 (six-line fixed point,
  public, NOT re-exported from `tanglepack.topology`), called ONCE at :1471 inside
  `symbolic_dynamics` :1361, after the closure loop (which records `words[cls]` :1465 = the
  class-level tokens, loops dropped, `None` when unresolved) and before `translate_itinerary`.
  Then: `entry.inert` set on every closed-table entry, letters in one pass (`_active_letters`
  :487 = an entry's own letter, else `letter_for`, else `letter(i)`; `inert_letters` :472),
  `entry.letter` on active entries only, `cd.letter` / `cd.kind` set. During the walk a
  `ClassDynamics` carries `letter = cls.label`, `kind = "active"`, so the "class %s unresolved"
  WARNINGs name the class by LABEL. The closure's `add` :1413 copies evidence entries with
  `letter=None, inert=None` and letters nothing. PK-8's WARNING in `_edges` is gone.
- `BridgeClass.py`: `BridgeClassEntry.inert: Optional[bool] = None` :206 (None until the dynamics;
  stays None on evidence entries outside the alphabet); `active` :213 = `not inert` (so EVERY
  entry is "active" before the dynamics -- `BridgeClassTable.active` returns all of them);
  `describe()` takes no argument, prints `undecided` for None. Deleted: `_resolve_inertness`,
  `_collect_image_evidence`, `has_image_evidence`, `image_classes`, `image_loops`, `unresolved`.
  `_image_chain` stays (loop folding, minimal trellis, trellis itinerary). Dev Note :51.
- `MinimalTrellis.minimal_trellis` :268 maps every hole bridge EXCEPT one carrying a ZEROTH hole
  (`hole.pair is not None`; author 2026-10-14), decided from `holes` alone; the class table is read
  only for the unclassed WARNING. `inert_hole_bridge_ids` = the zeroth-hole bridges; `unmapped` =
  the other hole bridges with no registered image. Dev Note :38 carries the argument: the lobe
  `b + sigma` of a zeroth hole holds one marker (its hole), whose image is no marker (backward-only
  holes), so `f(b)` is homotopic to `f(sigma)`: empty word. Proven under two conditions (no
  selected orbit inside the lobe, no other hole beside `b`); true on all 15 builds; SUFFICIENT,
  not the definition (heteroclinic `u`, hole-free, is inert by transience).
- `loom/Topology.py`: `_alphabet` :144 (cached; zones `minimal.alphabet`, no letters), `classes`
  :152 (plain property: touches `self.dynamics`, returns `_alphabet` -- `session.bridge_classes()`
  now RUNS THE DYNAMICS), `dynamics` :188 passes `_alphabet` and `letter_for`.
- Boundary case: the unblasted nested fixture (`conftest.p3_partitioned`) has a zeroth-hole class
  (bridge `(4,5)`) that is UNRESOLVED (singleton landing, no registered chain) and so ACTIVE while
  its bridge is not mapped; its sibling `(5,6)` goes active `c` -> inert `u`; `is_reliable` False
  on base and branch alike (verified by the auditor, 2026-10-14).

## How it works

- Hole side (`Hole.bridge_side`) = side of the hole's OWN bridge, in the dynamical orientation.
  A direct hole reads it off its crossing signs (2026-10-02, author-approved); a propagated hole
  CARRIES it from the hole it propagates from (2026-10-07; `diagnostics.side_source(hole)` =
  `"carried"`), flipped only under an orientation-reversing map (`side_after`, C16). The geometry
  (`diagnostics.measured_sides`, read at the nearest segment of the image sub-arc, computed on
  demand since B4a) is a resolution flag, never an assertion (NUM-2, WP4).
- Openings: at each non-anchor bound the hole opens the `outward` half iff `(crossing_sign > 0)
  == (bridge_side == "left")` (D6, WP4).

## Invariants

- I1: holes of one origin share `bridge_side` (parity-flipped under an orientation-reversing
  map). I2: a hole bridge with both ends on one stable branch has alternating signs (one row).

## Smells / open issues (current)

- The 2026-10-02 list is closed: whole-bridge nearest vertex (fixed 2026-10-02), INT-4 / INT-5
  (WP4), p3 >= 6 blasts unreliable (every class walks since WP7, flagged classes since WP8b,
  reliable to 20 blasts in E12), `merge_manifolds` non-contiguous images (NUM-9, WP10), the
  noisy inert-word WARNING (PK-8, WP5; the WARNING itself deleted in B5).
- Open (deferred / waiting): BCG-8 / DEC-25 (inversion loses classes at blast 6 by a genuine
  reference switch whose `+1` pair is never spanned, the reference arc lying outside the zone);
  INT-11 kind 2 = NUM-3 / NUM-4 (an image crossing the window 1e-7..1e-6 from a hole-bridge end
  at a crossing with no registered preimage: k32 from 9 blasts); the k32 closed alphabet drops
  45 -> 39 at blast 9 while the seed and holes grow (no law binds the closed alphabet across a
  switch; question for the author, with DEC-25).
- `topology.__all__` holds 48 names since B10 (commit `517aa01`, was 104): the 15 names imported
  THROUGH `tanglepack.topology` from outside the subpackage (the ten `tanglepack/__init__.py:45`
  re-exports Trellis, TrellisBranch, Arrangement, Arc, Region, Hole, PartitionInterval,
  PseudoneighborPair, StablePartitionResult, StrongPipResult; plus `side_after`, `letter`,
  `bridge_row_diagnosis`, `check_holes_share_bridge_side`, `plot_stable_partition`), 15 stage
  objects and 18 plotting names. Every `from .X import` line was KEPT, so the 56 dropped names stay
  importable by name (`from tanglepack.topology import StableNode` works) and `ruff F401` flags all
  56 in `topology/__init__.py`; trimming the import block is B12's. No test reads `__all__`.

## Features

- Re-measured by E12 (WP12, 1e-7): p3 9 classes to blast 5, 23 / 24 / 25 at blasts 6-8, then 25
  flat to 20, reliable throughout; k28 3, 3, 4 x 6, 15 from blast 8, reliable from blast 1;
  heteroclinic 57 seed / 81 closed flat to 8 blasts.

## Deferred / open

- CLOSED (WP3, C1): a stable branch with ONLY its anchor is one element `[anchor, end]`.
- FIXED 2026-10-07 (WP4 part A) and WP11 (C15): ownership, containment and `interval_at` are
  exact; the topology reads no absolute cdist slack.

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
  crossing; since D15 (2026-10-09) a missing `+j` appearance costs no hole. An unlinked landing's
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
  `cross(-u,-s) = cross(u,s)`). Users: `BridgeClass._image_chain` (image
  branch = `advance_key(unstable_key_of(bid), 1)`), `PartitionFamily.owner_at_end` (THE one
  bridge-end reader since B9b; used by the class rule `bridge_classes.ends`,
  `DualWalk.trellis_itinerary`, `SymbolicDynamics._member_pair`, `plotting._bridge_end_elements`),
  the iterated cut, `_row_mismatch_message`,
  `MinimalTrellis._bridge_sort_key`. `image_cdist(..., branch=)` maps an
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
  `_step_into`, `_image_arc_midpoint`) lives in `plotting.hole_point` since B4b: it places drawn
  points, the law oracles' points and `measured_sides`' reference point; nothing decides on it.

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
  reads registered crossings alone (NUM-18). Since B2 (branch `refactor/pseudoneighbor-both-intervals`
  @ 9cd29d5) `_interval_is_clear(..., "stable", ...)` rejects them: audited on p3_6_blasts, the two
  non-pip stable branches have 6 consecutive pairs each with a clear unstable and a punctured stable
  interval (0 kept by both); only the pip branch is walked, so nothing changes today. Filling is DEC-19 (WP5 slice 2; option
  (b), an unstable-side cut of N, changes the deep references, BCG-4).
- `_strong_pip_cuts(trellis) -> {stable key: cdist}` from `image_cdist(pip, j, "stable")`,
  `j < k_value`: one cut per branch always; a cut not from the table is flagged in one WARNING
  (PK-6). Never fires on a builder (cuts complete everywhere). `_cut_ids` and both
  `_pair_fixed_point` copies are deleted.
- `_is_forward_iterate(pair)`: `iterate is not None and iterate > 0` (D15, 2026-10-09, reversing
  D5's `+1 .. +(k_value-1)` exemption: no forward hole on any orbit). `punch_holes` partitions `forward` / `pairs` by two
  comprehensions; the no-spanning-bridge WARNING names each skipped pair, iterate, origin and
  branch (INT-9 stays WP7's).
- PK-9: `propagate_reference_holes` has no `max_steps`; it loops `count(1)` and ends because
  every step stops or adds a new bridge id to `visited` (finite). Longest chain seen: 21 steps.
- PK-8: `_edges` warned only when an inert class's word held a non-inert symbol (deleted in B5: under transience an inert word names only inert classes).
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

- Evidence vs alphabet (D2). [SUPERSEDED by B5: activity is transience over the words, the minimal trellis skips zeroth-hole bridges.] Activity is a fixed point over the images of ALL bridges
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
  **Resolved in WP10b (C19).**
- Residual laws on the deep runs (also on base): `partition_covers_branch` (BCG-4, WP8),
  `dual_face_side_is_geometry` (INT-10, WP10). Inversion x period 2 at 10 steps also fails
  `partition_covers_branch`, `bridge_endpoints_on_own_branch`, `backward_holes_keep_bridge_side`
  (WP8 / WP10).
- Deferred: escape-lobe images (WP9 / WP10); DEC-21 (a walk pair no bridge spans; answered D11, WP8b).

## Open-issues round WP8 (2026-10-08, branch `oi/wp08-recompute-identity`)

- Pseudoneighbor endpoint frame (BCG-5): in `_unstable_interval_is_clear` (since B2 `_interval_is_clear`, both sides), after the linked-landing
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
  pseudoneighbor all-iterates test keeps raw landings (`_interval_is_clear`, B2).
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
- D7 pair test (`Pseudoneighbor._unstable_interval_is_clear`, since B2 `_interval_is_clear`; the stable arc mirrors it with B's cycle): the cycle, `k` and β are the ARC's
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

## Open-issues round WP9b (2026-10-09, branch `oi/wp09b-noncrossing`, base `40175fb`)

The non-crossing rule (D13). Plan and probes: `/home/dezhu/.claude/jobs/c5c98371/tmp/wp09b-noncrossing_plan.md`,
`.../tmp/wp09b/` (planner: `nc.py` + `proto3.py` prototype, `truth2.py`, `geomcheck.py`, `amb.py`, `bh.py`,
`k28_struct.py`, `span_series.py`, `hetlinks.py`, `evcheck.py`), `.../tmp/wp09b_impl/` (coder: `cmp_lobes.py`
flood vs prototype, `sig.py`, `series.py`, `series_describe_all.py`, `breakit.py` / `breakit2.py`, `f7.py`,
`bn_missing.py`).

- **Span ends (C14).** `shortest_walks(..., landings=)` keeps only the start / goal nodes whose edge the
  landing span overlaps beyond its end error (`_nodes_meeting`: 0 at a table end, `SCALING_RTOL * |c|` at a
  scaled one; since B6a (branch `refactor/walk-bounds-honest` @8a1cdc3) TWO tiers: overlap > tol, else
  overlap > 0 (a span narrower than its own error keeps the one edge it lies in, k32 blasts 6/8), else
  `[]` + WARNING and `shortest_walks` returns status `"span_meets_no_node"` -- never all nodes). This alone removes every tie of k28 at
  8-13 blasts and k32 at 8-10 (both walks of each tie started on edges of the element the span does not meet)
  and fixes k32 6-7's `m` (a UNIQUE walk from the wrong edge, contradicted by its image). No-op on every law
  case and deep run.
- **Markers (`DualGraph.markers`).** Periodic point: `("orbit", position in trellis.fixed_points, orbit
  index)` at its anchor vertex (`("vertex", id)`). Hole: `("hole", origin, iterate)` beside its own bridge's
  unstable arc on `bridge_side` (`("beside", edge key, side)`; `face_side(arc)` gives the side a face lies on:
  right of each half-edge, dynamical direction hi->lo on S, lo->hi on U). (origin, iterate) is unique on every
  law case and deep run (asserted). A hole whose image is no hole (reference iterate 0; the last +j) gets that
  image placed beside the image chain's pair touching `iterate(near, 1)` (near is always a bound), same side
  (opposite under a reversing map), only when that pair is a kept arc (k10: not). The trellis-LOBE-corner
  placement (the side of the containing bridge facing its own trellis lobe) breaks the registered-image law on
  9 of 11 law cases (`breakit.py corner`): a propagated hole sits beside its container's image SUB-ARC (INT-2).
  A GLOBAL side flip passes the registered-image law (it is equivariant: f preserves sides) but makes `a` on
  k10 / nested / p3 / heteroclinic `no_admissible_walk`, which `member_lobe_classes_have_one_walk` catches.
  Orbit +2 passes every law (no lobe with a kept image chain separates the points of one orbit).
- **Enclosure (`DualGraph.enclosed`).** Pieces = ring arcs, a split stable edge (a chord ends at its middle)
  as two halves named by their end vertex. A chorded face: the ring holding both ends of every chord (else
  None), segments between chord ends glued across each chord, `len(regions) == m + 1` or the chords cross
  (`[]`: no simple routing). Every other boundary component of a chorded face (other rings of a merged face,
  `"infinity"` on the unbounded node) may sit in any of the `m + 1` regions: the product of placements, one
  dict per routing. Pieces are joined across every edge / half the curve does not run along. Outside = the root
  of `("component", unbounded.index, "infinity")`. Two chord ends on one side of one edge -> None (one crossing
  per edge, at its middle). Matches the prototype's flood on all 131 lobes of the 11 law cases
  (`cmp_lobes.py`), and point-in-polygon on every periodic point (law `lobe_enclosure_is_geometry`).
- **Lobe (`DualWalk._lobe`).** Both ends on ONE stable branch by `stable_key_at` (C7 at an anchor end), the
  same branch the class's elements and so the walk's ends lie on. The laws close a curve along the non-anchor
  end's branch instead (any stable path from the anchor is a Jordan curve), which adds inversion_period2's
  anchor bridges to the physics check; the library gains nothing from it (those classes land on two branches,
  so the walk curve cannot close: verdict None).
- **Verdicts (`holds_markers`).** True: some routing holds every decided wanted marker exactly. False: no
  routing does, the walk crosses an edge strictly inside sigma' (`f(b)` meets `f(sigma_b)` only at its ends),
  or chords cross (`[]`). None: ends on two branches, an edge crossed twice, a step through the start / goal
  edge (the midpoint model cannot draw it). A start node == goal node with no steps holds nothing. In
  `_admissible_walks` an itinerary is kept when one of its realizations (start node, goal node, face path) is
  True or None; None-only itineraries stay candidates (never preferred by order, so they can only make a tie).
  The prototype DROPPED None verdicts; with "crossing chords = False" (it returned None there) the outputs are
  identical (bh's second walk crosses its own chords in face 0).
- **Search.** Without `wanted`: the old BFS + backtrack (unchanged results; `start == goal` is now just the
  shared-face case, identical without landings). With `wanted`: L0 from the BFS (0 when a start face is a goal
  face), then `_walks(length)` for L0..L0+MAX_DETOUR (DFS, faces may repeat, a node crossed once, pruned by
  the BFS distance to the goal faces). The `max_walks` budget counts JUDGEMENTS (`holds_markers` calls): once an
  itinerary has a True realization its other realizations only add to its multiplicity (review fix: the first
  cut charged every realization and truncated the k32 series' trivial `u` at 6 blasts, 132 node pairs, and `ba`
  at 10, words unchanged; now never hit on the law cases, deep runs or series). Out of budget with nothing
  admitted = status `truncated`, no walk (was mislabelled `no_admissible_walk`). Status `trivial` iff the chosen
  length is 0.
- **Measured.** Law cases and deep runs: id-free signatures identical to base except `heteroclinic` (`a`
  unambiguous, `c` -> `a d^-1 a^-1`, `bh` -> `bi bj bi^-1` via a 2-crossing detour, ambiguous -> [bn], verified
  False -> none) and `k28_deep_blast_sliver` (reliable, hash 794bc6ac4f -> 0e8153ad32), both equal to the
  prototype. Series: k28 8 (unreliable [ab, af] -> reliable), k32 6-7 (`m` -> reliable), 8 (44 -> 45 closed,
  reliable), 9-10 (43 -> 39 closed, ambiguous [ck, br, ct] -> none, still unreliable through INT-11); p3,
  inversion, nested_outer unchanged; k28 9-13 and k32 7-10 describes identical to the prototype's.
  `symbolic_dynamics`: heteroclinic 0.04 -> 0.11 s, deep sliver 0.11 -> 0.15 s.
- **Residue (DEC-24) -- resolved in WP9c by the walls (below).** `bn` (member-less, connecting, named by `bl`):
  its two walks differ by a loop around the B peninsula; none of the heteroclinic crossings (92, 93 in N; 108,
  109 past p3's unstable cut) has any iterate link (`bn_missing.py`), so no closed curve with a known image goes
  through it (NUM-21). The walls need no such curve.

## Open-issues round WP11 slice 1 (2026-10-09, branch `oi/wp11-simplify`, base `4f7e815`)

- **No absolute cdist slack (C15).** Anchors by `is_anchor` (StrongPip x2, Pseudoneighbor,
  MinimalTrellis, `Arc._interior_polyline`); exact: `cut_at_pip_orbit` (placed cuts too),
  `_anchorward_look`, `DualGraph._meets_fundamental_segment`, `Arc.elements`, `interval_at` (no
  fallback, no `tol=`, `PartitionFamily.tol` gone); `_build_intervals`' near-anchor INFO branch gone;
  `land_element._end_tolerance` = 0 at a table end, `SCALING_RTOL * |c|` at a scaled one. Public
  `tol=` of `is_strong_pip` / `classify_strong_pips` / `compute_pseudoneighbors` (+ Trellis wrappers)
  removed. Kept (geometry): `Arc._interior_polyline` and `_oriented_bridge_polyline` clipping. The
  cartoon layouts identify boundaries exactly (`tol` 0; boundaries are registry floats). Counter
  (`tol_count.py`): 115 / 19 / 56 -> 59 / 2 / 51 (tolerance sites / `cdist_tol` sites / warnings).
- **One orientation rule (C3, C16).** `TopologyResults.side_after(side, steps, preserving)`, exported;
  `Trellis.orientation_preserving` reads `FixedPoint.preserves_orientation` (no jacobian call, no
  cache; DEBUG default True when no fixed point is solved).
- **One collision test.** `Trellis.same_point(u, s, u0, s0, rtol)`; the branch-cycle wrappers are
  gone (use `FixedPoint.branch_position_map`).
- **One interval builder.** `StablePartition.intervals_from_marks` (+ `Marks`), shared with
  `PartitionFamily`; `_build_intervals` only adds the branch's ends.
- **One element image.** `land_element`; `Trellis.image_of_element`, `SNAP_RTOL`,
  `StableNode.images`, `_link_element_images` deleted. In `land_element` (and its singleton path) an
  image end on a branch other than `advance_key`'s is an assert (was WARNING + unresolved).
- **One cartoon drawer.** `plotting._draw_cartoon_graph(target, layout, dual, dynamics, arc=, cross=,
  track=, ...)`; module split deferred.
- Smaller: `owner_of_intersection(id, side, branch_key)` (required, a direct lookup);
  `_anchor_branch` inlined in `_hole_openings`; `bridge_classes(trellis, partitions)` (no
  `bridge_ids`); `_diagnose` uses `groupby`; dead `Symbol.is_inverse`, `dynamics_of`,
  `TrellisBranch.anchor_coords`, `BridgeClassTable.member_of` deleted.
- Proof harness (`/home/dezhu/.claude/jobs/c5c98371/tmp/wp11/check_all.sh <root> <tag>`): `sig.py`
  (14 + 5 builders), `cartoon_sig.py`, `series_sig.py` (5 series), `tol_count.py`; every hash identical
  to base after every step.

## Open-issues round WP9c (2026-10-09, branch `oi/wp09c-every-bridge-the-same`, base `5c7070e`)

- **Where `bn` comes from (D14).** `named_by`: `bn <- bl <- bi`; `bi = {B:L_(0.0;1), A:R_(0.0;5)}` (A = p1,
  B = p3) has two COMPUTED members `(91, 92)` +1 / `(93, 94)` -1 on u_B(2.0), arcs from B's stable branch to the
  heteroclinic crossings 92 / 93 on A's window. `bl = f(bi)`, `bn = f^2(bi)`, member-less only because f(92),
  f(93) lie past B's unstable cut (D12, NUM-21): D11. Landings table-exact and contained (A:R_(0.0;5) = [5,4] ->
  A:R_3 = [110,6] -> inside A:R_1 = [0,111]; B:L_(0.0;1) -> B:L_(1.0;1) -> B:L_(2.0;1) -> B:L_(0.0;1)^1). Face 0
  is one open face, one ring of 59 arcs, holding B's left-side edges on all three B branches and A's right-side
  edges, with corners 91..94 (the B arms run through it): a cross-tangle pair there is legitimate. A's
  partition is sound. Probes `wp09c_bn.py`, `wp09c_trace.py`, `wp09c_part.py`, `wp09c_face0.py`.
- **The gap.** D13 was implemented only as the marker law on a member LOBE (`DualWalk._lobe` is None across two
  stable branches), so a class joining two stable branches or two tangles had no rule -- members or not: at
  17u the connecting class `bm` has two registered members and is ambiguous on base.
- **The rule (C17, walls).** `SymbolicDynamics._arc_branches` (`SymbolicDynamics.py:1090`): a class's arcs lie on
  its members' unstable branches (`Trellis.unstable_key_of`, loops included), or one `advance_key` step past
  its namer's. In `symbolic_dynamics` (`:1347-1352`) the decided walks (`unique` / `trivial`) of the classes
  walked EARLIER whose branches are disjoint from this class's are its walls (`walk_chords`,
  `DualWalk.py:831`); `shortest_walks(walls=)` (`DualWalk.py:1229`) drops, from a tie only, the walks with a
  chord crossing a wall in its face (`DualGraph.chords_cross` `DualGraph.py:1372`: four
  distinct ends on one ring, alternating; `_ring_index` shared with `_cut_ring`), INFO, counted in
  `WalkSearch.rejected`; when every walk crosses a wall the tie stands (never observed). The walls apply to
  member-lobe ties too (after the marker filter). Only a walk of `multiplicity == 1` is a wall or is judged
  by one (fix round, reviewer: a trivial walk's face is `shared[0]` of `len(shared)`, a unique itinerary
  counts its other paths, so the chord read is an arbitrary representative; heteroclinic had u, bj at x5 and
  v..y at x2 as walls, latent). The filter is a comprehension in `shortest_walks` (`_crosses_a_wall`
  removed); `DualGraph._ring_of` reads `_ring_index`. The member-less entry's `tangle` takes the namer's fixed
  points into the table's rule (a no-op everywhere).
- **Measured.** Id-free and exact `describe()` signatures (`wp11/sig.py`, `PYTHONHASHSEED=0`) of the 11 law
  cases, `p3_15`, `p3_16`, `p3_6_blasts`, `p3_16_12`, `inversion_period2_10`, `k28_deep_blast_sliver` identical
  to base except `heteroclinic` (exact 60bc15796d -> 972f82e7c1, the prototype's hash: `bn` ambiguous ->
  resolved with the same word `a u^-1 bi^-1`, reliable, DualWalk WARNINGs 1 -> 0); `series_sig.py` identical on
  the five series. The law prototype found no crossing pair of decided walks anywhere. p3 17u (`heteroclinic_17u`,
  27 s): `bm` resolved to `bh uu^-1 bi^-1`, its registered word, reliable. `symbolic_dynamics` on heteroclinic
  0.104 -> 0.105 s; k32 / nested_outer series 14.0 / 14.7 s unchanged.
- **Deferred.** Walls from LATER classes (a second pass; never needed); same-branch disjointness (distinct kept
  bridges of one branch; measured identical, not needed); NUM-20 / NUM-21 untouched, no longer needed for words.

## Open-issues round WP10b (2026-10-09, branch `oi/wp10b-leftovers`, base `f9d1be0`)

- PK-4 / C19: `IteratedHomotopyPartition._empty_stretches` pairs a lobe image's interior crossings
  per stable branch: chain order on the base branch (author, 2026-09-23), STABLE order on any other.
  [SUPERSEDED by B7b (branch `63477e3`, 2026-10-15): STABLE order on every branch, base included.]
  Proof: a lobe H (bridge + stable segment between its ends) holds no point of any stable branch of
  the trellis (a bridge crosses none in its interior), so f(H) holds no anchor (f(H) holds z_j iff
  H holds f^-1(z_j)); a branch walked from its anchor starts outside f(H) and enters / leaves at
  each crossing with f(b). At an inversion anchor that is a corner of f(H) the other half of the
  stable curve leaves on the opposite side of the sector, so it starts outside too. Measured by
  point-in-polygon (`test_foreign_chords.py`, `pk4_geom.py`): p3_6_blasts, p3_15, p3_16,
  heteroclinic all read out, IN, out, IN, out from the anchor.
- BASE branch is NOT checked by `test_foreign_chords.py` and chain order can contradict the plane
  (B7 analyst, 2026-10-14): k32_series blast 6, hole bridge `(6,7)` image `8..13` all on `s(0,0)`,
  stable order `13,8,11,10,9,12`; point-in-polygon says `(11,10)` and `(9,12)` are inside `f(L)`,
  chain order marks `(9,10)`/`(11,12)` empty; the sibling `(54,55)`'s image crosses at 93/94 inside
  `(10,9)` -> Theorem T fails at P1 (see the B7 section at the end).
  [B7 branch `a6a4f68`: recorded as the strict xfail `(series_members_share_the_P1_itinerary,
  k32_series)`, wp B7b; the fix (stable order on the base branch too) is slice B7b. B7b analyst
  (2026-10-15): the same lobe `8..13` exists from blast 2 to blast 10, and its chain and stable
  orders differ at every one of those steps; T fails only at steps 6 and 7. B7b LANDED on
  `refactor/base-chords-in-stable-order` @ `63477e3`: the base branch is now checked by
  `test_foreign_chords.py::test_base_chords_are_inside_the_lobe_image` and the series test; known issue gone.]
- The base is read at a NON-anchor end (an anchor's `stable_key_at` is the C7 row rule, which on
  `inversion_period2` named the other half: three "ends on another stable branch" WARNINGs) and
  asserted to lie on one branch; the "no flank" skip in `from_minimal` is an assertion for a
  non-anchor crossing (every empty stretch lies on one branch).
- Effect: only the iterated partitions of heteroclinic / p3_15 / p3_16 / p3_6_blasts (+3 elements
  each: a chord cut at its outward end and a pinched singleton at the inward end), hence iterated
  names, iterated itineraries and refined-pair labels; class tables, holes, words, refined rules,
  matrices and reliability identical (`sig.py`, `rr.py`). PartitionFamily WARNINGs 2 -> 0 on those
  four, 3 -> 0 on inversion_period2; inversion_two_blasts keeps its one "no empty stretch" (the
  anchor lobe's fold, by the author's rule).
- Left: the `BridgeClass` cross-branch skips (`_image_chain` off-branch image end,
  `_collect_image_evidence`'s uncovered branch) fire on nothing measured; assertion candidates for
  the D8 API pass. The "no empty stretch" text is misleading at the anchor lobe (WP12 docs pass).
- NUM-19 bullet 4 (walk vs class ownership) is stale since C15: both read
  `element_of_intersection` / `owner_of_intersection` or the exact `interval_at`; the only slack is a
  scaled landing end's relative error.
- Plotting: `CartoonLayout.tol` / `ZoneLayout.tol` and `_cluster` deleted; ranks are
  `sorted({float(c) ...})` of the partition boundaries, `_rank_of(representatives, cdist)` exact
  (every cartoon artist identical, `wp11/cartoon_sig.py`). `span_contains` / `owns_cdist` docstrings
  name their real slack.

## Open-issues round WP12 (2026-10-09, branch `oi/wp12-reverify`, base `a2e0e65`)

- C21 (INT-11 kind 1): `SymbolicDynamics._order_flip(trellis, homotopy, cls, evidence)`, called by
  `_diagnose` after the stable-end test: a member end and a bound of its OWN element (the class's
  homotopy element at that end, read through `homotopy.element`) whose registered `+1` images
  run in the reverse stable order on one branch. The map preserves order along a branch, so the
  registered cdists contradict the dynamics: noise, `cause = "order flip"`, `crossings = (end,
  bound)`, `word_agrees` None, INFO with the two gaps. Lookup only. Measured: k28 member
  (246, 247) at 9 blasts, end 8.5e-6 outward of bound 4, images 1.2e-6 the other way, both
  resolution ratios ~1e4. k28 at 9-13 blasts becomes reliable; k32 at 9+ stays unreliable through
  kind 2. Test `test_an_order_flip_gives_no_verdict` (`build_k28(blasts=9)`, fails on base) and
  the independent cause scan of `test_disagreeing_evidence_carries_its_cause`. The synthetic
  `_drive` harness gives its homotopy elements no bounds.
- Unmarked lobe (`PartitionFamily`): a crossing abutting no empty stretch that lies on the image
  chain of a hole bridge whose hole has no stretch (`_unmarked_lobes`, computed lazily on that
  path; `_hole_stretches` is the shared stretch set) is skipped as `"unmarked lobe"` at INFO, the
  author's 2026-09-23 rule; anything else stays the `"no empty stretch"` WARNING without the old
  "grow or blast" advice (the fold IS registered). Every occurrence on the law cases and deep
  builders is this kind (inversion_two_blasts 1, k28_deep_blast_sliver 13; 13 / 47 per recompute
  on k28 / k32 at depth). Test `test_an_unmarked_lobe_is_skipped_by_design` (fails on base).
- `BridgeClass._image_chain` asserts an image end on the advanced branch (was a WARNING + None,
  never fired: table links pass the map gate onto `advance_key`); `_collect_image_evidence` lets
  `ends(pair)` raise like a member's (the try/except never fired). Test
  `test_an_image_end_off_its_branch_asserts` (fails on base).
- DOC-3: the `TopologyResults` "is a Hole a coordinate" note is settled (a hole is a marker beside
  its own bridge, D13).

## Oracle laws over the three PDF definitions (B0, 2026-10-12, branch `refactor/oracle-laws`)

- Strong pip, pseudoneighbor pair and hole containment are now checked by laws that never call
  `StrongPip`, `Pseudoneighbor` or `StablePartition` (`tests/helpers/laws.py:1566`, `:1632`, `:1672`;
  details in `testing.md`). All pass on the 11 law cases on base `a3bd8f0`: no code change was needed.
- Physics learned: a reference bridge end's preimage can BE a crossing of N that is an end of the
  container (k28_one_blast `iterate(5,-1) = 3`; inversion `iterate(3,-2) = 1`), so containment is
  closed there; the real inverse map overshoots it by 1e-7..7e-4 relative. The alignment plan's A.0/A.6 record this
  (patched 2026-10-09: "table first" is real at a reference bridge's far end, scaling only along the span).
- Check 2 found no punctured reference interval on any law case, the stable interval included (Q1's
  stable-side test, slice B2, should change nothing on law cases).

## B1 `resolved-crossings` (2026-10-12, branch `refactor/resolved-crossings` @ `76d0710`, NOT merged at writing)

- `Trellis.resolved_ids` (`Trellis.py:204`): N as the registry decided at snapshot time;
  `from_workbench` :422/:447 puts only crossings of N on a branch (an unresolved detection is
  registered, in `trellis.registry`, on NO branch, so not in `own_intersection_ids` or `in trellis`).
  Hand-built trellises (and the `MinimalTrellis` sparse snapshot) default to the whole registry.
- `StrongPip.is_strong_pip` :169 skips any `r_id` not in `trellis.resolved_ids`.
- The pip-orbit cut reads `iterate`, and an unresolved detection is never linked: when the pip's
  next orbit point is unresolved, `m` stops one iterate early and STAYS there under further blasting
  (k28 x8 = x7, k32 x9-10 = x8, nested_outer x9-11 = x8). Remedy: rebuild at a finer `area_cutoff`.
  Open author question: freeze vs placing such a point by `image_cdist` (PK-6 style).
- Law `check_trellis_is_cut_at_the_pip_orbit` (`tests/helpers/laws.py:1506`) now reads N (skips
  `iid not in trellis.resolved_ids`).

## B6a `walk-bounds-honest` (2026-10-12; findings at `76d0710`, outcome audited on `refactor/walk-bounds-honest` @ `8a1cdc3`, not merged)

A walk is a word only when established. Three sites now leave the class unresolved (`itinerary is
None`, so `is_reliable` False through its unchanged formula):
- **Truncated search** (`SymbolicDynamics._resolve_itinerary`, the `elif search.truncated` after the
  status check): reason "search budget reached before uniqueness was established", whatever the status
  (unique / trivial / ambiguous). `shortest_walks` statuses and `max_walks=64` unchanged. The wall gate in
  `symbolic_dynamics` also needs `cd.itinerary is not None`, so an unresolved class lends no wall.
- **Span meeting no node** (`DualWalk._nodes_meeting`, two tiers, see "Span ends" above): new
  `WalkStatus` `"span_meets_no_node"`; reason "dual-graph walk span_meets_no_node".
- **Members disagreeing on a marker**: `lobe_image_markers` raises `DualWalk.MembersDisagree(ValueError)`
  (NOT exported from `tanglepack.topology`; `lobe_image_markers` is), caught in `_resolve_itinerary`
  BEFORE `shortest_walks` (never walk unconstrained instead); reason "members hold different markers:
  the class is not one homotopy class (resolution; <detail>)". None still means "no member lobe".
- The walk-path WARNING now reads `class X unresolved: walk from A to B: <reason>` (was `... is <status>`).

Measured (auditor + TE, base vs branch): no rule / reliability / closed-size change on any law case,
p3_6_blasts, p3_16, k28 x8 or any step of the five SERIES. heteroclinic at `max_walks=2`: reliable True
-> False (`B:p3@2.0/L#0 <-> L#2`); at 1 two more classes become unresolved (already False); 8 and 64
unchanged. The span WARNING fired 10 times on k32 blasts 6+8 at base, 0 on the branch (second tier keeps
the one node; words identical). Neither the empty-list path nor `MembersDisagree` fires on real data:
only hand-built unit tests guard them.
- Stale outside scope: `CLAUDE.md` still says "a span meeting no node falls back to all of them with a
  WARNING" and lists `WalkSearch.status` without `span_meets_no_node`.
- Open numerics question (plan "not doing"): `SCALING_RTOL * |c|` (~0.012) is coarse next to spans of
  0.001-0.006; the second tier papers over it.

## B4a `holes-algorithm-vs-diagnostics` analyst (2026-10-13, HEAD `ae18859`)

- The 18 names + `_map_backward` all exist at the plan's names (spans in the findings
  `artifacts/findings_analyst_planner.t1.md`). Moving them verbatim, minus the two checks the plan
  keeps, removes ~599 lines: StablePartition.py lands at ~1,348, not <= 1,100 (that figure needs B4b
  too).
- Coupling the plan does not show: the kept checks DELEGATE to moved names (I2 ->
  `bridge_row_violation` -> `rows_of_bridge` + `bridge_row_diagnosis` for the message; I1 ->
  `bridge_side_violations` -> `side_diagnoses` -> `_orbit_side`, `SideDiagnosis.describe`), and the
  new module needs `row_of_end` / `rows_of_bridge` from StablePartition: a module-level import both
  ways is circular. `owns_cdist` is not used by any moved name.
- `_map_backward`'s carried point also PLACES every propagated hole (`_punch_in_bridge` :1585-1603).
  Without the map the existing fallback (`_step_into(midpoint, stable_mid)`) puts 1-29 propagated
  holes per law case on the wrong side of their image sub-arc (today: 0, except inversion_period2 4/8
  and heteroclinic 13/69 by the same reading). Fingerprints do not include `coords`.
- Measured side is reproducible on demand (verified): for a propagated hole, carried =
  `map_inv^(-iterate)(reference.coords)`, span = sorted `image_cdist(id, iterate, "unstable",
  branch=ref_bridge.manifold_key)` over the reference bridge's ids, read by `_bridge_side_of` on
  `bridge_between(*bounding_ids)` -- equals the stored `measured_side` on all 128 propagated holes of
  the 11 law cases.
- `side_source`: `"signs" iff hole.pair is not None` (today's assignment) and `iff iterate == 0`
  agree on every hole of the 11 law cases (no direct hole has iterate != 0).
- Baseline: `geometry_flags(trellis, hole_bridges, holes)` = `([], [])` and `side_diagnoses` = `[]`
  on all 11 law cases (topology trellis), so the rubric-5 equality is vacuous on them; only
  `p3_coarse_blasts` exercises the flags (2 row, 20 measured, WP4 notes).
- `plotting.plot_holes` (:237-300) reads `Hole.origin`, `.coords`, `.iterate` -- not `bridge_side` or
  `openings` (the plan's wording); no src reader of `measured_side` / `side_source` outside
  StablePartition.py.

## B4a outcome audit (2026-10-13, branch @ `f166573`, base `ae18859`)

- Verified: commit 1 (`95ea780`) moves the 16 names byte for byte (ast segments, decorators
  included, same order, none left behind; the only other StablePartition edits there are the two
  lazy imports in the checks and one in `_punch_in_bridge`, plus import narrowing). Full suite
  1173 passed / 38 skipped / 5 deselected / 1 xfailed (orientation_reversing only). Base-vs-branch
  dump (holes incl. coords, measured sides, row diagnoses, flags, fingerprints; 11 law cases +
  p3_coarse_blasts) identical, md5 74385b04...
- Smells left: (1) I1/I2 assertions run code that lives in diagnostics (`bridge_side_violations` ->
  `side_diagnoses` -> `_orbit_side`; `bridge_row_violation` -> `bridge_row_diagnosis` for the
  message) -- "nothing in diagnostics feeds a decision" is true of flags, not of these cores.
  (2) With `diagnose=True` the inverse map now runs twice per propagation step (placement +
  `measured_sides` replay). (3) `measured_sides` counts from iterate 0 and keys references by
  origin: an orbit started from the propagation FALLBACK (reference with no spanning bridge,
  StablePartition.py ~:455) is unmeasured; never seen on any build, B4b removes the fallback.
  (4) Stale outside the file list: `numerics/geometry.py:7` docstring,
  `scripts/refactor_proof/hole-side-and-own-blast/hole_sides.py` (`sp._bridge_side_of`,
  `sp._oriented_bridge_polyline`, `sp.bisect_left`), `CLAUDE.md` :35 and the `Hole` row (:139).

## B4b `holes-are-combinatorial` analyst (2026-10-13, HEAD `f166573`)

Verified by `b4b_probe.py` (session scratchpad): an instrumented replica of
`propagate_reference_holes` checked against the original on every call (identical holes, coords
included), over every `BUILDERS` case and every step of the five `SERIES`.

- **Dead paths.** No recorded pair has `iterate < 0`, on any build or any series step. The
  `seen` and `punched_regions` dedups never reject anything. `_punch_in_bridge` never returns None,
  and no pair is skipped for a missing `arc_point`. Each origin has one reference. The fallback
  (`propagate_reference_holes` :455-495) is entered only on `near_tangency`: 14 references there
  have no bridge, no orbit has a punched pair, so it logs 14 "cannot propagate" WARNINGs and
  nothing else. Why: a reference lies in the window `S \ f(S)`, so the `direction = -1` walk in
  `Pseudoneighbor.extend_pseudoneighbor_trajectories` (:256-283) breaks at once, because one of
  the two preimages is not in T. That walk is itself dead (out of scope for B4b).
- **Table span ends are real.** The reference bridge's far end has a REGISTERED preimage. It is a
  crossing of N, and it is exactly the container's own end crossing, on 6 of the 11 law cases
  (k28_one_blast 1, p3 3, nested 3, nested_inner_blasts 3, inversion 2, inversion_period2 4) and
  on every series. Examples: `iterate(5, -1) = 3` with container `(0, 3)`; inversion_period2
  `(14, 15)` gives `15 -> 5, 9, 2, 13`. So "never in the table" is false for one end (INT-2's "at
  most one end registered"), and `from_table` is what separates an exact end from a scaled one.
- **Scaled-end proximity.** A scaled end lies within `SCALING_RTOL * |c|` of a crossing of the
  container's branch on no law case (closest: 0.0199 on the inversion blasts, 0.022 on
  k28_two_blasts). It does so on `k28_deep_blast_sliver` (18 ends) and `k32_series` (242 ends, 5
  per-branch firings, minimum relative distance 1.2e-4).
- **Placement.** `hole_point` with the reference placement carried back by `map_inv` (the same
  per-step `asarray(...).ravel()`) reproduces `Hole.coords` BIT FOR BIT, max deviation 0.0, on
  every hole of every builder, direct holes included. The no-map fallback `_step_into(midpoint,
  stable_mid, 0.05, arc_chord)` moves holes by up to 0.328 (k28_one_blast) and flips the side of
  the image sub-arc on every law case: k10 1/2, k28_one 1/2, k28_two 1/3, p3 3/6, nested 4/9,
  nested_inner 4/9, inversion 2/4, inv_one 1/6, inv_two 1/6, inv_p2 1/8, heteroclinic 30/69;
  p3_coarse_blasts 24/80.
- **Anchor openings.** On the inversion cases `_hole_openings` names `s(o,1-b)` at 2/4, 2/2, 2/2
  and 4/8 anchor ends, as the plan says. On the non-inversion cases the D6 formula says
  "anchorward" at about half the anchor ends (k10 1/1, k28_one 1/2, k28_two 1/1, p3 3/6, nested
  4/7, heteroclinic 4/10), and there NO opening is recorded, because `s(o,1)` is not a branch. The
  `key in trellis.branches` guard is load-bearing.
- `_near_far` is NOT placement-only: it orders `bounding_ids` / `near_intersection_id` (a
  fingerprint field, and `laws.check_holes_are_classified`). It stays.
- Readers of `Hole.coords` outside StablePartition / plotting:
  - `diagnostics.measured_sides` :469 (the reference point);
  - `StablePartition._hole_marks_on_branch` :1090 (a debug log);
  - three laws in `tests/helpers/laws.py` (:1213, :1254, :1280);
  - `test_stable_partition.py:370`;
  - `Hole(coords=...)` constructions in `test_topology_invariants.py:68`,
    `test_stable_partition.py:82, :392` and `test_minimal_trellis.py:71`.

## B4b outcome audit (2026-10-14, branch @ `ba23a1f`, base `f166573`)

- Verified: full suite on the branch 1178 passed / 38 skipped / 5 deselected / 1 xfailed
  (orientation_reversing only), 205 s. Diff = the nine plan files + one `TangleSession.py`
  docstring line (`git diff --stat`). StablePartition.py grep for `epsilon|coords|np.|map_inv|
  dynamical_system|_map_backward|...` empty; imports are logging/bisect/itertools/typing + local.
  Own probe (`scratchpad/aud_b4b/fig.py`): `plot_holes` scatter offsets and iterate labels identical
  (max|d| = 0.0) on k28_two_blasts (4), inversion_period2 (10), heteroclinic (73); five fields equal.
- Edge behaviour changes, unreachable on every build (B4b analyst probe): a bridge with no nodes
  used to skip its hole (`arc_point is None`, `_punch_in_bridge` -> None) and now gets one; a pair
  with `iterate=None` now raises TypeError in the assertion. `hole_point` on a node-less bridge would
  fail.
- Open: I1/I2 cores in diagnostics (B4a); plotting calls the map (view); CLAUDE.md :35 / :139 still
  describe `measured_side` / `side_source` as hole state (B12); `hole_point`'s bare StopIteration.

## B2 analyst pass (2026-10-14, `HEAD` ba23a1f; line numbers are BASE's, see "B2 outcome" below)

- **Map at HEAD.** `Pseudoneighbor.py` (471 lines): Dev Notes `:7-127` (stable "empty by construction"
  `:64-66`; orphan pairs `:50-58`; unstable-only closed form `:83-89`; check set `:117-126`; NO "NUM-18"
  text in the file, NUM-18 lives at `:215` above and in the audit); `compute_pseudoneighbors` `:149-227`
  (one call `:217`, candidates `own_intersection_ids` `:186`); `_unstable_interval_is_clear` `:329-471`
  (table by id `:410-429`, own frame `:434-458`, scaled + `same_stable` 2-D guard `:460-470`). There is
  NO stable-side code: consecutiveness via `next_intersection` (`:209`, `:224`) over branches holding the
  cut N only (`Trellis.from_workbench:433-446` resolved only; `cut_at_pip_orbit:264`).
- **Measured: the mirrored stable reading never punctures** (every consecutive window pair, 11 law cases +
  p3_6_blasts + k28 x8, and the five `cases.SERIES` to the E12 depths k28 13 / p3 20 / inversion 10 /
  k32 13 / nested_outer 12 at the recipe cutoff: 73 steps, 0 fires). Simulating a missing window
  representative (dropped from candidates AND table, as B1 does) the stable reading catches it via the
  OWN-FRAME reading about half the time, but the unstable reading already rejects every such merged pair;
  B2 is definitional on the corpus.
- **Candidate-set trap.** Resolved crossings PAST the unstable pip-orbit cut sit inside the open stable
  interval of real references (law case inversion: (2,3) holds 23, 24 at u=28.3 vs cut 25.9, no table
  orbit; k32_series@13 all 6 references; nested_outer@11-12 (38,37)). They are not in
  `X = U M^n(N_T)` (their orbit's only window visit is themselves). Candidates must stay
  `own_intersection_ids`; reading the registry or `resolved_ids` would kill law-case references.
- **m-range with beta < 1.** Mirroring `:399-401` with `lap = beta_s**k` makes `log_lap < 0`, so the
  bounds swap; correct only while `hi/lo < beta^-k` (true on the one-lap window). Use min/max.
- Since B1 the iterate table links only N (`numerics/IterateInference.py:91,117,122`): an unresolved
  detection is never a table landing; on k28_series from blast 6 the run 75-77 (u=8.4e6, unlinked) sits
  in the non-reference window pair (64,102) and fires nothing.
- N plateaus on every deep series (k28 143 from blast 7, k32 300 from 8, nested_outer 276 from 8, p3 119
  from 14, inversion 167 from 9): the pip-orbit stall (B3). Attributed at `9cd29d5` (B3 analyst):
  k28 / k32 / nested_outer / inversion = an unresolved `f^{m+1}` (ids 245 / 523 / 447 / 238);
  p3 / heteroclinic = `f^{m+1}` never detected (the image chord misses `s(0,0)` past its anchor);
  never `min_separation` at one iteration per call (see loom.md).

## B2 outcome (2026-10-14, branch `refactor/pseudoneighbor-both-intervals` @ `9cd29d5`, base `ba23a1f`; not merged)

- **Map.** `Pseudoneighbor.py` 489 lines (471 at base). `_interval_is_clear(trellis, x0, x1, stability,
  candidates, collision_rtol)` `:343` (one body, 147 lines vs 143): `attrgetter` picks `arc_key` /
  `other_key` / `cdist` from `stability`; stable side asserts both ends share `manifold_b_key`;
  `beta = owner.per_step_beta(stability)` `:394` (< 1 on the stable side); m-range
  `range(floor(min(bounds)), ceil(max(bounds)) + 1)`; 2-D guard `:424` `same_other = other_key(r)[0] is
  other_key(ix0)[0]` (D7 and mirrored); local `point()` re-orders to `same_point`'s (unstable, stable);
  one debug format `"Pair (%d, %d): %s interval punctured by %d at M^%d, %s reading (...)"` with reading
  `table` / `own frame` / `scaled` (e13 parses `record.args`). Caller `compute_pseudoneighbors` `:229` (candidates `:197`):
  `all(... for stability in ("unstable", "stable"))`. Candidates still `own_intersection_ids`.
- **Behaviour.** None on any law case (fingerprint, references, all pn pairs, words identical, two
  independent base/branch scripts). e13 (`stress/numerics/e13_both_intervals.py`, results jsonl committed,
  produced at `530c1ba`): 14 cells (7 E12 cases x 1e-7/1e-5), 0 stable fires on 182 steps, before ==
  after, after == library everywhere. Inversion still 2 -> 1 references at blast 6 ((23,24),(26,27) ->
  (56,57)), holes 8 -> 5, seed 8 -> 6, both cutoffs: DEC-25 is not the stable clause.
- **Invariant relied on.** Every interval the window walk asks about is <= one lap, so min/max of the
  m-range is untested by any real build; only the uncommitted hand check (TE `r6.py`) covers it.
- **Smells.** The stable branch's `assert key == arc_key(ix1)` would fire on a direct call with two
  branches (no caller does); the `other_key(r)[0] is other_key(ix0)[0]` guard is now the only place D7
  is encoded for both sides -- a third reader of the same rule is `tests/helpers/laws.py`
  `check_reference_pairs_clear_of_X` (independent by design).

## B5 `inertness-is-transience` analyst (2026-10-14, HEAD `7d3a294`)

Inertness map at `7d3a294` (what B5 deletes or moves):
- `BridgeClass.py`: module Dev Notes "INERT = MAPS TO NO ACTIVE CLASS" :51-71; slots `inert`,
  `image_classes`, `image_loops`, `unresolved` :230-233 (docstring :201-213); `active` :239-242;
  `has_image_evidence` :250-253; `describe()` prints the evidence :312-353; call sites in
  `bridge_classes` :712-714; `_collect_image_evidence` :783-823; `_resolve_inertness` :826-849.
  `_image_chain` :725 stays (loop folding, the minimal trellis, the trellis itinerary).
- `MinimalTrellis.py`: `inert_hole_bridge_ids` field :91/:117, `summary` :216, `describe` :229-233,
  docstring :266-276, the mapping loop :310-334 with the `entry.active` test :326 (an inert hole
  bridge goes to `inert` and is not mapped); constructor arg :425. DEBUG text "active hole bridge %s
  has no registered image" :331-333.
- `loom/Topology.py`: `classes` :143-152, lettering loop :148-150 (`if entry.active: entry.letter =
  self.letter_for(...)`), `zone_key` :151; `dynamics` :182-185 passes `letter_for`.
- `SymbolicDynamics.py`: `inert_letters` :466-478 (reads `table.inert`), `_active_letters` :481-483,
  `translate_itinerary` kind from `entry.inert` :551, `symbol_nodes` :871-899 (active only), `_edges`
  :904-933 with PK-8's WARNING :919-925 (fires only when an inert class's word names a NON-inert
  symbol; never fires on the corpus), driver `symbolic_dynamics` :1336: initial letters :1383-1386,
  closure `add` :1389-1415 (letters by the copied entry's evidence inertness), kind :1426, the
  class-level tokens `named` :1446-1449 (the same `oriented_class(parent_of x, parent_of y)` as
  `translate_itinerary`: the inertness fixed point can run on them BEFORE translation), translate
  loop :1462-1466.
- Readers of evidence-table / session-table inertness OUTSIDE B5's seven files: `tests/facade/
  test_session_bridge_classes.py:73-80` (`session.bridge_classes()` letter iff active, "holds a loop
  => inert", without running the dynamics), `tests/unit/topology/test_minimal_trellis.py:94`
  (unmapped hole bridges are active), `tests/helpers/direct.py:107-125` (letters `entry.active`
  before the dynamics), `tests/cases.py:1192-1193` (SeriesStep active/inert, read after the
  dynamics ran), `tests/regression/test_blast_children_are_evidence_only.py:53` (dynamics.table:
  fine), `scripts/henon_bridge_classes.py:275-318`, `stress/numerics/e11_blast_depth.py:359-393`.

Probe (scratch, monkeypatching `loom.Topology.bridge_classes` to mark every evidence entry active,
`topology(rebuild=True)`; transience = fixed point "word empty or every token inert", unresolved
active). 11 law cases + p3_6_blasts, p3_16, k28 x8, nested inner 6 (verified, run 2026-10-14):
- Class-level letter-free words IDENTICAL on all 15; closed-alphabet class sets identical;
  `is_reliable` True on all 15 both ways; no unresolved class anywhere (so "unresolved is active"
  is never exercised on the corpus); the transience computed from HEAD's own words equals I-a's.
- Inert verdict changes (all active -> inert, all "active for lack of evidence" at HEAD):
  inversion_two_blasts `g`; heteroclinic `d`, `bj`; p3_6_blasts `ab`; nested inner 6 `ae`. None
  inert -> active. Class matrices equal on the surviving nodes (the moved classes drop out).
- Image bridges / iterated elements: k28_one +1 / +3; p3 +2 / +5; nested and nested_inner_blasts
  +2 / +5; inversion +2 / +5; inversion_period2 +2 / +5; heteroclinic +2 / +6; p3_16 +2 / +6; the
  rest 0 / 0. Refined ACTIVE nodes change only by a new split on heteroclinic `bk` (1->2) and p3_16
  `ab` (1->2), plus the dropped inert classes.
- Side effects: newly mapped inert hole bridges with no registered image land in `unmapped`
  (k10, k28_two, nested, nested_inner, inversion_one, p3_6 +1; nested inner 6 +2) and a new
  "image pair ... has no Bridge object ... skipped" WARNING appears on k28_one_blast,
  inversion_two_blasts and k28 x8 (one each, D9).
- Every alphabet class holding a loop member is inert under I-a on all 15 builds.

### B5 round 2: the author's zeroth-hole rule (2026-10-14, HEAD `7d3a294`, verified on the 15 builds)

Author: "every bridge that has an initial hole punched in it should be inert" (a zeroth hole =
`Hole.pair is not None` = `iterate == 0`; its bridge = `trellis.bridge_between(*hole.bounding_ids)`).
- Every zeroth-hole class has an EMPTY word (walk `trivial`, one loop pair) on all 15 builds. The
  five classes I-a flips (`g`, `d`, `bj`, `ab`, `ae`) are exactly the zeroth-hole classes HEAD
  calls active for lack of registered images. The set of HEAD-inert hole classes is contained in
  the zeroth-hole classes, so (I-a inert) = (zeroth-hole classes) ∪ (HEAD inert). No propagated-only
  hole class is inert. Every loop-holding class is a zeroth-hole class. No bridge carries both a
  zeroth and a propagated hole. On the corpus, inert means "the word is empty".
- The converse is false: heteroclinic `u` (bridge `(4,5)`, the middle piece of the anchor bridge's
  image chain `(0,4) -> [(0,4),(4,5),(5,6)]`) holds no hole, maps to the loop `(6,110)`, and is
  inert. So the rule is sufficient, not the definition.
- Why it holds (sketch; a theorem given two conditions): the pseudoneighbor clarity of both open
  arcs leaves the lobe `b ∪ σ` with no arc of T and exactly one marker, its zeroth hole, provided
  no selected orbit's whole tangle lies inside it. The image `(origin, +1)` is no hole
  (`DualGraph.image_marker` returns None, `DualGraph.py:1172`), so `f(b)` is homotopic rel the
  markers to `f(σ)`. `f(σ)` has no T crossing inside it, because the preimage of such a crossing
  would be an X point on `σ`. So the class's image is a loop.
- Prototype (monkeypatch `loom.Topology.bridge_classes`): with inert = (zeroth-hole OR HEAD
  evidence), the minimal trellis, `P1`, the image bridges, the refined classes, the words and
  `is_reliable` are IDENTICAL to HEAD on all 15 builds. The verdicts equal I-a's, there are 0 new
  WARNINGs, and `unmapped` shrinks by 1 or 2 on 4 builds. With zeroth-hole alone, the verdicts
  equal I-a's except heteroclinic `u`, which becomes active.
- Recommended B5 shape: the minimal trellis maps every hole bridge except the zeroth-hole ones,
  read from `holes`, not from the class table. Transience in the dynamics stays the verdict for
  every class. A new law, `zeroth_hole_classes_have_empty_words`, checks the theorem.
  Scratch scripts: `.../scratchpad/b5r2/{holes,probe2,pre}.py`.


## B7 analyst (2026-10-14, `HEAD` 3fb310a): refined level, Theorem T, in-splitting

- Code: `refined_rules` built at `SymbolicDynamics.py:796-809` (one copy of the parent's refined word per
  child), `refine` :606 (occurrences from ALL words, inert too), `_refined_symbol` :656 (unmatched token
  returned bare), `match_members` :813, `_edges` :940, `transition_graph` :960, `transition_matrix` :995
  (silently drops a non-node target :1012), `translate_itinerary` :525, `DualWalk.trellis_itinerary`
  :1301. Odd-count WARNING `PartitionFamily.py:959-966` inside staticmethod `_empty_stretches` :892
  (called once at :748); base-branch chain order at :957; `IteratedHomotopyPartition` has no `describe`
  override (base :413).
- In-splitting (siblings' rows equal; children's columns sum to the parent's `M` column; node set =
  refined active symbols) holds on all 15 builds + 47 series steps (verified). Free reduction never fires
  on member itineraries.
- Theorem T (members' P1 occurrence sequences equal): tested classes k10 1, inversion_two_blasts 2,
  heteroclinic 7 (p3 tangle), k28 x8 4, p3_16 4; series total 141 class-tests (p3_series 7-9, inversion
  2-4, k28 3-8, k32 3-10, nested_outer 3-11). FAILS: (1) k28 x8 / k28_series 6-8 member `(44,45)` --
  a "stable end" chain (C11), excluded by the `_diagnose` run criterion; (2) k32_series 6 and 7, class
  `R#7<->R#9`, hole-bridge member `(6,7)` vs `(54,55)`/`(64,65)`/`(96,97)`: base-branch chain-order chords
  (above). Monkeypatching stable order on the base branch too removes (2) and leaves M, M_1, names,
  counts, `is_reliable` identical on all 62. [CORRECTED by the B7b analyst, 2026-10-15: the "inert
  splits drop on 15 builds/steps" figure came from a scratch patch (`b7/patched_empty.py`) that also
  lost the `empty.setdefault(b, ...)` line, so every empty stretch was indexed at its smaller id only.
  The correct patch changes no split count anywhere; see the B7b section at the end.] Odd-count WARNING
  fires on none of the 62. Scripts: `.../scratchpad/b7/{probe,probe_patched,dig,jordan}.py`.

## B7 outcome (2026-10-14, branch `refactor/refined-level-defined` @ `a6a4f68`, base `3fb310a`, NOT merged)

- Map at the branch (supersedes the line numbers of the B7 analyst section above for these names):
  `PartitionFamily.py` 1,009 -> 1,082 lines. Module docstring states `P1 = P ∨ f(∂P)` (:18-27) and the
  base-chain / foreign-stable chord rule + incomplete footprints (:38-47); Dev Notes History condensed
  (:88-92); the "changes refinement labels only, never walks or words" claim is REPLACED by the k32
  counterexample and B7b (:107-118). `IteratedHomotopyPartition` :645, `incomplete_footprints` kwarg /
  attribute (:669, :685), `describe()` override :687, `from_minimal` :698 passes
  `cls._incomplete_footprints(minimal, homotopy)` :912, `_empty_stretches` :927 (signature unchanged,
  dict return; base-branch chain order at :993; the odd-count WARNING deleted), NEW staticmethod
  `_incomplete_footprints` :1007 (Counter over `chain[:-1]` second crossings per stable branch; one
  `logger.info` :1044 "an odd number of fold crossings ...").
- `SymbolicDynamics.py` 1,516 -> 1,538, docstrings only (AST-equal with docstrings stripped): module
  "refined level" paragraph :17-31 (occurrences, Theorem T, `M_1` = in-splitting of `M`), Dev Notes :64-67,
  `refined_rules` attribute doc :774-777, `transition_matrix` :1012 docstring names `M_1`. Construction
  unchanged: `refined_rules` :815-826 (one copy of the refined class word per child), `refine` :622,
  `_refined_symbol` :672, `match_members` :830, `translate_itinerary` :541.
- Invariants now stated: (T) every member of a class has the same freely reduced `P1` itinerary;
  (in-split) sibling rows of `M_1` equal, children's columns sum to the parent's `M` column, nodes =
  children of split active classes + unsplit active letters. Both are TEST laws (tests/helpers/laws.py),
  not library assertions. The in-splitting's sibling-row identity holds by construction of
  `refined_rules`; the column-sum identity catches a token dropped by `transition_matrix`.
- Measured by the auditor (`.../scratchpad/b7a/{k32,walk}.py`): k32_series 0-7 identical on base and
  branch (letter-free words, `M`, `M_1`, reliable, T tests 0,0,0,2,3,4,9,13, one disagreement at 6 and
  7, `incomplete_footprints` empty). The disagreeing class is `{R_(0.0;8), R_(0.0;10)}` (ElementNaming;
  `R#7`/`R#9` positional). The WALKED word agrees with the siblings `(54,55)` / `(64,65)` / `(96,97)`;
  only the hole-bridge member `(6,7)` differs (its `L_1^4 L_3^2`, `R_19^4` vs the walk's `L_1^5 L_3^1`,
  `R_19^3`). Stronger check (walk vs EVERY readable member at P1 level) agrees on all 11 law cases (128
  member readings, 0 disagreements); not a law.
- Smells: `_incomplete_footprints` repeats `_empty_stretches`' chain walk (two lookups per build) and
  its admission filter differs in principle (holes' stretches only vs the growing dict; measured equal
  on 62 builds/steps) -- fold together in B7b, updating the monkeypatches in `test_iterated_cut.py` /
  `test_iterated_partition.py`; module docstring of SymbolicDynamics says T and the in-splitting are
  "laws on every build" (T is law-cases + series only).
- Features: `IteratedHomotopyPartition.incomplete_footprints` / `describe()` (unit test with a stub
  chain; never populated on the corpus); refined matrix documented as `M_1` (laws + unit test on A.8d's
  k=2.8 shape).

## B7b analyst (2026-10-15, HEAD `a6a4f68`, slice `base-chords-in-stable-order`)

- Map at HEAD: `PartitionFamily.py` module docstring bullets :31-47 (the chain/stable clause :39-41,
  the odd-count clause :43-47), stale k=2.8 example :49-58 (pre-WP8 ids: "hole (6, 5), base (10, 7),
  chord (8, 9)"; at HEAD the k28_two_blasts lobe is hole `(6,7)`, base `(8,11)`, chord `(9,10)`, as
  CLAUDE.md says), Dev Notes History :88-92, Foreign chords :101-118; `from_minimal` calls
  `cls._empty_stretches` :778 and `cls._incomplete_footprints` :912; `_empty_stretches` :927 (Note
  :946-956, chord comment :985-988, `if branch_key != base_branch:` :993); `_incomplete_footprints`
  :1007; `_hole_stretches` :1070. `SymbolicDynamics.py` "laws on every build" :27-28, the k32/B7b
  sentence :29-31.
- Callers of `_empty_stretches` outside the module (all tests): two monkeypatch WRAPPERS
  (`test_iterated_cut.py:158-170` `with_foreign_far_ends`, `test_iterated_partition.py:90-100`
  `without_bases`; both call the original and rebuild a dict from `.items()`), four direct calls
  (`test_iterated_cut.py:85` harness used by :106/:120/:129, `test_iterated_partition.py:70`,
  `test_foreign_chords.py:60`, `test_partition_family.py:169`, with `_incomplete_footprints` at :170).
- Base-branch lobes measured on all 16 builds (11 law cases + p3_15/p3_16/p3_6_blasts + k28 x8 +
  nested inner x6) and 47 series steps: 47 stretched-hole lobe readings carry 2+ interior crossings on their
  BASE branch; only ONE lobe carries 4+: hole `(6,7)` on `k32_series` steps 2-10 (9 of the 47) (chain `8..13`, stable
  order of the interior `11,10,9,12`; the base `(8,13)` is adjacent in stable order, as it must be:
  `f(b)` meets `f(sigma)` only at the ends). Point-in-polygon (`b7b/probe.py`, arc-length midpoints,
  `tol=0`, distance-to-ring >= 5% of the gap length everywhere) reads every lobe: the gaps alternate
  out/in from the anchor with the base counted as inside (the inversion lobes start AT the anchor, so
  their first gap is the base); the covered gaps are the stable-order pairs on every lobe, = the chain
  pairs on the 38 two-crossing readings; chain order is never the covered set where they differ.
- Outside gaps can still be empty for another reason (another hole's stretch): inversion `(3,2)`,
  inversion_period2 `(15,14)`, inversion_series@0, k32 steps 2-5 `(10,9)`. A law "no outside gap is
  empty" is false; "no outside gap is THIS lobe's base or chord" is the right form.
- Correct patch (drop the test only), base vs patch on all 63: letter-free `rules`, `M`, `M_1`, class /
  active / inert counts, split-class counts and `is_reliable` identical; T 0 disagreements everywhere
  (k32 6-7 from 1); in-splitting ok everywhere. Differences: iterated elements on k32 2-10 only (2-5:
  +1 element, the right-row crossing 9 pinched between hole `(10,9)` and chord `(9,12)`; 6-7: -1, the
  singleton `[10,10]` goes; 8-10: -2); `unmatched_members` 1 -> 0 on k32 6-7 (member `(6,7)`);
  `refined_rules` on k32 8-10: inert class `{L_(0.0;3), L_(0.0;5)}` (positional `L#2<->L#4`) has 4
  children at HEAD and 3 patched (member `(9,10)` joins `(93,94)`), so one active word's inert token
  index changes (`v_4 -> v_3`); `M_1` unaffected (inert not a node).
- Teeth: the HEAD unit pin `test_foreign_crossings_pair_in_stable_order` (`test_iterated_cut.py:92-110`)
  gives all OWN crossings `stable_cdist=1.0`, so a stable sort keeps chain order and the pin PASSES under
  the patch (toothless for B7b); the k32 hand-built shape (cdists 8:0.176, 9:6.786, 10:6.665, 11:6.571,
  12:7.418, 13:0.036) gives chords `(9,10),(11,12)` at HEAD and `(10,11),(9,12)` patched.
- Odd counts: every base-branch interior count on the corpus is even (0, 2 or 4); under stable order the
  unpaired crossing of an odd count is the OUTERMOST, not the chain's last.

## B7b outcome (2026-10-15, branch `refactor/base-chords-in-stable-order` @ `63477e3`, base `a6a4f68`, NOT merged)

- Map at the branch (supersedes the B7b analyst line numbers above): `PartitionFamily.py` 1,082 -> 1,077
  lines. Module docstring states ONE chord rule (:39-47: stable order on every branch, C19 off the base,
  B7b on it; odd count leaves the OUTERMOST crossing unpaired); the k=2.8 example now in current ids
  (:53-58: hole `(6,7)`, base `(8,11)`, chord `(9,10)`, verified against `iterated_partition()` on
  k28_two_blasts); Dev Notes History one condensed chain (:88-95, the only "chain order" mention); Dev
  Notes "Chords" (:104-116, both branch kinds, no known issue). `_empty_stretches` :925, the chord loop
  `for on_branch in folds.values(): on_branch.sort(stable_cdist)` :988 (no base-branch exemption; the
  `(base_branch,) = ends` binding is gone, `assert len(ends) == 1` stays); `_incomplete_footprints` :1002
  (code unchanged; docstring says "outermost"); `_hole_stretches` :1065. `SymbolicDynamics.py:27-30`:
  T and the in-splitting are laws "on the law cases and every blast series"; "the base-branch chords pair
  in stable order since B7b".
- Behaviour (auditor re-diffed the test engineer's probe JSONL, `scratchpad/te/out/`, 63 items): only
  `k32_series` steps 2-10 differ -- iterated elements 23->24 (2-5), 40->39 (6-7), 113->111 (8-10);
  `unmatched_members` 1->0 and T disagreements 1->0 at 6-7; `refined_rules` 55->54 keys at 8-10 (inert
  `{L_(0.0;3), L_(0.0;5)}` loses its 4th child; one active word's token `v_4 -> v_3`). Split-inert
  counts are EQUAL on every item (1 at k32 2-7, 3 at 8-10, both sides): the B7 "-1 split inert on 15
  builds" claim was a buggy-patch artefact, now confirmed on the landed code.
- Smells left: `_incomplete_footprints`' INFO log (:1040 "last one unpaired"), `describe()` line (:690
  "the last crossing unpaired") and the class attribute doc (:653 "its last crossing") still say "last";
  `_incomplete_footprints` still repeats `_empty_stretches`' chain walk (fold declined, Q1 default:
  it would change the return shape at six test call sites); CLAUDE.md "Partition families" bullet still
  says chain order on the base branch (B12).
- Features: the iterated cut `P1` now marks the k32 lobe's covered gaps `(11,10)`, `(9,12)` empty, so
  Theorem T holds on every series step (verified: series test passes, p1_tested
  `[0,0,0,2,3,4,9,13,13,13,13]`).


## B6b `word-finished` analyst (2026-10-15, HEAD `63477e3`; B7b merged)

- Map (verified at HEAD). `DualWalk.py` 1,370 lines: module docstring :1-105 (the truncation rule
  :94-104), `WalkStatus` :557-567 (nine values), `WalkSearch` :575-634 (no `reason` field; `truncated`
  flag separate from status), `MembersDisagree` :790 (raised by `lobe_image_markers` :795; NO
  WalkSearch is built for it, `cd.search` stays None), `shortest_walks` :1027, status returns :1109
  no_start_node, :1118 no_goal_node, :1123 span_meets_no_node, :1153 trivial, :1164 unreachable, :1224
  truncated (constrained, no walk), :1234 no_admissible_walk, :1280-1294 trivial/unique/ambiguous
  (+ `truncated=True` flag with a candidate kept), `trellis_itinerary` :1301.
  `SymbolicDynamics.py` 1,537 lines: module docstring :1-154 (evidence bullet :80-99), `Symbol` :211,
  `EvidenceDiagnosis` :267-296, `ImageEvidence` :299-333, `ClassDynamics` :336-413 (`status` :404-413,
  four values), `inert_classes` :514, `translate_itinerary` :540-610, `_signature` :613, `refine` :621,
  `registered_evidence` :687-736, `_member_pair` :739-745, `SymbolicDynamics` :751-1149 (`match_members`
  :829 CALLS `_member_pair` :852; `is_reliable` :899; `describe` :1070 prints `cd.status`, never
  `search.status`), `_resolve_itinerary` :1171-1265 (status mapping :1250-1253, `cd.ambiguous` :1265),
  `_diagnose` :1268, `_order_flip` :1299, `_compare_evidence` :1323-1379, `symbolic_dynamics` :1382-1537
  (wall gate `search.status in ("unique","trivial")` :1477; class-level `words` for inertness and
  closure built from RAW pairs :1481-1487, NOT from `translate_itinerary`; translate :1515;
  evidence calls :1518-1519).
- Evidence names: the six of the plan (`EvidenceDiagnosis`, `ImageEvidence`, `registered_evidence`,
  `_diagnose`, `_order_flip`, `_compare_evidence`) = 223 lines of defs; with separators and the
  section header removing them leaves SD at 1,297. `_member_pair` is ALGORITHM (refinement's
  `match_members`) and is monkeypatched by `test_symbolic_dynamics._drive` :567 -- it does not move.
  Consumers outside SD: `topology/__init__.py` :108-116 / :179-187 only; tests call
  `SD.registered_evidence` (:470, :496) and patch `SD._image_chain` / `SD.trellis_itinerary` to drive
  it (:466-467, :492-493, `_drive` :567-569).
- Import graph (verified on a scratch copy): `diagnostics` may import `BridgeClass._image_chain` and
  `DualWalk.trellis_itinerary` at module level and SD may import `diagnostics` at module level (no
  cycle from any entry module); `diagnostics` importing SD at module level is a cycle (ImportError
  from `import tanglepack`). `translate_itinerary` / `_signature` must reach diagnostics at call time
  (function-level import, the pattern `measured_sides` already uses for `plotting.hole_point`) or be
  passed in.
- Corpus facts (63 items, `scratchpad/b6b/probe.py`): free reduction fires 0 times under either rule
  (class-level adjacent `c c^-1`, same-occurrence) over 902 class words and 1,830 evidence readings;
  WalkSearch statuses seen are only `unique` (151) and `trivial` (402); 279 trivial (zero-crossing)
  walks give a ONE-TOKEN word, so WalkSearch "trivial" != ClassDynamics "trivial" (empty word, 151);
  `truncated` flag never set at `max_walks=64`; the 3 unresolved classes (k28/k32/nested_outer series
  step 0) are singleton landings with no chain (`cd.search` None).
- Smells: the class-level word graph (inertness, closure) and `cd.symbols` are read by two code paths
  (raw pairs vs `translate_itinerary`); a reduction added to one only would split them. The two
  hand-built fixtures `test_loops_dropped_from_word_but_recorded` and
  `test_evidence_disagreement_names_a_stable_end` are exactly `a a^-1` with DIFFERENT occurrences:
  class-level reduction cancels them, same-occurrence does not.

## B6b outcome (2026-10-15, branch `refactor/word-finished` @ `65edec4`, base `63477e3`; NOT merged)

What landed. The analyst section above describes `HEAD` before the change; this one describes the
branch.

- **Map.**
  - **`DualWalk.py` (1,383 lines).**
    - The module docstring opens with A.8c's two definitions (object, algorithm) and the list of
      unresolved causes.
    - `WalkStatus = Literal["resolved", "ambiguous", "unresolved"]` (:569).
    - `WalkReason` (:572) holds six values: `unreachable`, `no_start_node`, `no_goal_node`,
      `no_admissible_walk`, `truncated`, `span_meets_no_node`.
    - `WalkSearch` (:587) has `reason: Optional[WalkReason]`, set iff the status is unresolved. The
      `truncated` field is GONE; `_Budget.truncated` stays internal.
    - `MembersDisagree` (:799) is now exported.
    - `shortest_walks` (:1036) returns:
      - a zero crossing or one itinerary -> `resolved`;
      - several itineraries -> `ambiguous` (WARNING);
      - a budget cut with candidates -> `unresolved`/`truncated` with `walks` KEPT (WARNING "the
        search decides nothing", no ambiguous WARNING);
      - a constrained search out of budget with no walk -> `unresolved`/`truncated`;
      - the other five failures -> `unresolved`/<reason>.
  - **`SymbolicDynamics.py` (1,312 lines).**
    - The module docstring opens with A.8c's definitions.
    - `free_reduce(tokens) -> kept indices` (:545) is a one-pass stack reduction to the fixed point,
      DEBUG per cancelled pair. It is exported.
    - `_resolve_itinerary` (:1070) sets `unresolved_reason` from `status == "unresolved"`: the budget
      sentence for `truncated`, else `dual-graph walk <reason>`. `cd.ambiguous = status ==
      "ambiguous"`.
    - The wall gate is `search.status == "resolved" and multiplicity == 1`.
    - The closure/inertness `words` are freely reduced: a class named only inside a cancelled pair
      is NOT added.
    - `cd.symbols = [symbols[i] for i in free_reduce(_signature(symbols))]`.
    - `verify_words(...)` is called once after all symbols (:1301).
    - `_member_pair` stays at :638.
    - `translate_itinerary` is UNREDUCED (loops dropped only).
  - **`diagnostics.py` (990 lines, exactly at the plan's cap).**
    - `EvidenceDiagnosis` :733, `ImageEvidence` :765, `registered_evidence` :801 (its `symbols` are
      freely reduced; `itinerary` stays raw), `_diagnose` :857 (takes `position` from the caller),
      `_order_flip` :887.
    - `verify_words(records, trellis, naming, table, letters)` :911 is the folded `_compare_evidence`.
      It reads `naming.iterated` and `naming.homotopy`.
    - Imports: `_image_chain` / `trellis_itinerary` at module level; `_signature`, `free_reduce` and
      `translate_itinerary` inside the functions; SD types under `TYPE_CHECKING`. There is no
      module-level import of SD (that would be a cycle).
    - The logger is `tanglepack.topology.diagnostics`.
  - **`plotting._class_status`** (:1808) returns `unresolved: <reason>`, else `trellis`, else
    `cd.status`. The `-` fallback is gone. It is untested.
- **Invariants.**
  - Free reduction removes whole adjacent pairs, so the word parity and the even-length itinerary are
    untouched; it acts on tokens, never on the itinerary.
  - The corpus has 0 cancellations (63 items, 899 class words, 1,830 evidence readings; three oracles
    agree). Refine and inertness are therefore unchanged.
  - Same-occurrence vs class-level differ only on hand-built `a a^-1` fixtures.
- **Smells.**
  - `free_reduce` runs twice per class (closure tokens from `oriented_class` of homotopy parents,
    then `cd.symbols` from `translate_itinerary`): two derivations of one word, already present
    before B6b.
  - `diagnostics.py` has zero headroom under its cap.
  - `verify_words` uses one squashed multi-argument logger call and the name `e` in `verdicts` to
    meet the cap.
- **Features / verification.**
  - Statuses:
    - six-reason parametrized test `test_every_failure_is_unresolved_with_its_reason`;
    - `test_a_cap_cutting_a_tie_short_keeps_its_candidate_unresolved`;
    - `test_an_ambiguous_search_marks_its_class_and_is_no_wall`.
  - Reduction:
    - `test_free_reduction_cancels_adjacent_inverses_only`;
    - `test_a_walk_reading_c_cinv_is_the_empty_word`;
    - `test_a_class_named_only_inside_a_cancelled_pair_is_not_added`;
    - `test_the_registered_word_is_reduced_before_it_is_compared`.
  - Base vs branch on all 63 items: equal (TE probe). The auditor's own spot check on k10 and
    k28_two_blasts gave equal words, statuses, `verified` and `is_reliable`; only `search.status`
    differs (`unique`/`trivial` -> `resolved`). On k28 two blasts, class `b` (word `c`) had walk
    status `trivial` on base: the two-meanings case.
- **Stale elsewhere:**
  - `CLAUDE.md`: the symbolic-dynamics bullet and the `DualWalk` / `SymbolicDynamics` table rows
    still describe the nine statuses, `WalkSearch.truncated` and the evidence in SD (B12).
  - `scripts/symbolic_figures.py:268` prints the new walk statuses unedited.

## B9 analyst: one element at a bridge end (2026-10-15, `HEAD` `517aa01`)

Four implementations answer "which element owns this bridge end", plus the law tier's own oracle:
- `PartitionFamily.owner_at_end` :343 -- any family (homotopy or iterated), the row (`row_of_end`),
  branch `stable_key_at(id, unstable_key_of(bid))`, any BridgeId pair (image pairs too), ValueError.
  Callers `SymbolicDynamics.py:642-643`, `DualWalk.py:1366`, `plotting.py:2303`.
- `Trellis.element_for` :1381 -- homotopy (`self.stable_partitions`), ANY side passed in (not just
  the row), `Optional[int]`; only src caller `TangleSession.partition_element_for` :1409; tests
  `test_partition_elements.py` x8.
- `StablePartition._elements_at_bridge` :1059 -- stored as `StablePartitionResult.elements_at_bridge`
  (`TopologyResults.py:442`) by `StablePartition.py:536` (homotopy) and `PartitionFamily.py:898`
  (iterated, over `minimal.sparse.bridges`); WRITTEN ONLY: no src reader, tests
  `test_partition_elements.py:77,87,158`.
- `BridgeClass._element_at_end` :763 -- homotopy `partitions` list, the row; the class rule
  (`ends` :595-598, folding `_non_loop_parent` :740); ValueError / AssertionError.
- Test oracle `tests/helpers/laws.py:1868` `_element_ends` (own `_stable_branch_at` :126).
- Agreement probe (scratchpad `b9/probe/probe.py`, 63 items = 16 builds + 47 series steps, the B10
  set): 7,090 / 7,090 ends agree on the row (all four + the oracle + `partition_element_for`),
  14,180 / 14,180 (end, side) pairs (`owner_of_intersection` = `element_for` = `elements_at_bridge`),
  1,776 / 1,776 iterated ends (`owner_at_end` = iterated `elements_at_bridge`; full vs sparse trellis
  identical). No disagreement.
- Folding: `PartitionFamily.py:136` imports from `BridgeClass`, so `BridgeClass` cannot import
  `PartitionFamily` at module level (cycle); the shared core belongs beside `row_of_end` in
  `StablePartition.py`. `element_for`'s any-side reading = `owner_of_intersection(id, side,
  stable_key_at(...))`.

## B9b analyst: one element at a bridge end (2026-10-15, `HEAD` `5f0db5f`, B9a merged)

Supersedes the line numbers of "B9 analyst" above (same facts, re-derived after B9a):
- `PartitionFamily.owner_at_end` :343-367 -> `self.owner_of_intersection` :369-398 (ValueError on a
  missing (branch, side) OR an unowned crossing, one message; pinned by
  `test_partition_family.py:109-110`). Callers: `SymbolicDynamics.py:642-643`, `DualWalk.py:1366`,
  `plotting.py:2303`; `test_dual_walk.py:877` binds it onto `_OwnerFamily` (needs only
  `owner_of_intersection(id, side, branch_key)` on self).
- `Trellis.element_for` :1379-1430 (with its section header): NO src caller; tests
  `test_partition_elements.py:68,105,106,158,168`; stale cite of the deleted session method :1413.
- `StablePartition._elements_at_bridge` :1059-1091 (end of file): writers `StablePartition.py:536`
  (homotopy) and `PartitionFamily.py:898` (iterated, `minimal.sparse`), import `PartitionFamily.py:140`;
  field `TopologyResults.py:415-427` doc + `:442-444`; readers tests only (`:74,84,155`). Its
  `assert element_id is not None` is the only thing it guards; the law tier covers it
  (`check_partition_unique_owner` laws.py:1168-1183: `element_of_intersection` == the branch's ids;
  `check_iterated_unique_owner` :2214-2226).
- `BridgeClass._element_at_end` :763-791: callers `ends` :595-599 (used at :604 and by
  `_non_loop_parent` :740). Differences from `owner_at_end`: ValueError on `branch_key is None`
  (dead on computed crossings; `owner_of_intersection` would raise TypeError formatting `None[1:]`),
  a distinct missing-partition message, and AssertionError (not ValueError) on an unowned crossing.
  No test pins the AssertionError; `test_bridge_class.py:311-329` pins only `ValueError` on a
  missing partition (both raise it).
- Probe re-run at `5f0db5f` (scratchpad `b9b/probe/`, the b9 harness minus `partition_element_for`,
  guarded so it runs on base and branch, with a per-item sha256 of every reading):
  7,090 / 7,090, 14,180 / 14,180, 1,776 / 1,776, 0 exceptions, 0 bad.
- Import cycle: `topology/__init__.py` imports `BridgeClass` (:59) before `PartitionFamily` (:82);
  a module-level `from .PartitionFamily import ...` in `BridgeClass` fails in EITHER import order
  (ImportError on `_image_chain`, verified). Two cycle-free shapes, both verified on scratch copies
  (suite 1207 + 47 + 1 xfail with the rewritten `test_partition_elements.py`; probe sha-identical on
  63 items; 63-item outcome probe byte-identical): (B) a function-level
  `from .PartitionFamily import HomotopyPartition` in `bridge_classes`, `ends` = `family.owner_at_end`,
  `_element_at_end` deleted, `BridgeClass` then imports nothing from `StablePartition` (src -152);
  (A) `end_on_row` + `owner_in` beside `row_of_end` in `StablePartition.py`, called by
  `owner_at_end`/`owner_of_intersection` and by the class rule (src -133).
- Docstrings: `StablePartition.py:23-32`, `:499-504`; `TopologyResults.py:400-403` ("two lookup
  tables"), `:415-427`, `:429` ("Both tables"); `Trellis.py:1413`; `BridgeClass.py:542-544`
  (`manifold_b_key` picks the branch -- stale since C7) and `:569-581` (Raises / Note on the
  AssertionError). `ResonanceZone.py:295` and `:353` are NOT stale: `TangleSession.classify_bridge`
  (:749-773) still orders zones by area and `TangleSession.add_resonance_zones` (:696) exists.
- Smells left: `StablePartition.py:238` imports `Intersection` unused (pre-existing, ruff F401).

## B9b outcome: one element at a bridge end (2026-10-15, branch `refactor/one-element-at-end` @ `94fca6b`, base `5f0db5f`)

Current state (line numbers at `94fca6b`; the two analyst sections above are history):
- **One reader.** `PartitionFamily.owner_at_end` (`PartitionFamily.py:342-366`, body unchanged)
  = `row_of_end` for the side, `stable_key_at(id, unstable_key_of(bid))` for the branch (C7), then
  `self.owner_of_intersection` (`:368`), which raises ValueError on an uncovered (branch, side) OR
  an unowned crossing. The any-side reading is `owner_of_intersection` called directly on the
  `stable_key_at` branch (what `test_partition_elements._owner` does).
- **The class rule** (`BridgeClass.bridge_classes` :529) builds
  `family = HomotopyPartition.from_results(partitions, trellis=trellis)` behind a FUNCTION-LEVEL
  import (:582; `PartitionFamily.py:136` imports `_image_chain, _index_partitions` from
  `BridgeClass`, so a module-level import is a cycle in either order) and `ends` (:596-600) calls
  `family.owner_at_end` for "first" and "second". `BridgeClass` imports nothing from
  `StablePartition` any more. `_index_partitions` (:748) stays: `PartitionFamily.__init__` uses it,
  so duplicate partitions are still a ValueError.
- **Gone** (recover with `git show 5f0db5f:<path>`): `Trellis.element_for` (Trellis.py:1379-1430,
  the only ANY-side reader returning `Optional[int]`, None when uncovered),
  `StablePartition._elements_at_bridge` (:1059-1091) and its two writers (StablePartition.py:536,
  PartitionFamily.py:898), the field `StablePartitionResult.elements_at_bridge`
  (TopologyResults.py:442), `BridgeClass._element_at_end` (:763-791).
  `StablePartitionResult` has four fields; positional construction with four args is unaffected.
- **Error paths changed** (no corpus item, test or law reaches them): an unowned crossing in the
  class rule is now ValueError (was AssertionError, Q2); the missing-partition ValueError carries
  `owner_of_intersection`'s message ("crossing N is on no 'row'-partitioned branch (o, b) of this
  homotopy partition family"; it no longer suggests "pass its trellis's stable_partitions"); the
  dead `branch_key is None` ValueError is gone; the build-time assertion in `_elements_at_bridge`
  is gone (unique ownership is the laws `check_partition_unique_owner` /
  `check_iterated_unique_owner`, plus the partition's own build assertions).
- **Verified** identity on the 63-item corpus: agreement probe 7,090 / 14,180 / 1,776, readings
  byte-identical to base; 63-item outcome probe byte-identical (auditor rerun, scratchpad
  `b9b/aud/`).
- Smells left: `StablePartition.py:238` unused `Intersection` import (pre-existing, ruff F401);
  `topology/__init__.py` 56 F401s (B12); CLAUDE.md `Trellis` row still lists `element_for` (B12).

## B12 analyst (2026-10-15, `HEAD` `94fca6b` = B9b merged, read-only)

Where the code stands against the alignment plan's section A after fifteen slices (full list in the
B12 findings artifact, `artifacts/findings_analyst_planner.t1.md`):
- A.7's chord rule is now STABLE order on every branch (B7b); A.8b's "maps every hole bridge" is
  "every hole bridge except one carrying a zeroth (reference) hole" (B5, `MinimalTrellis.py:13-47`,
  `:274-321`); inertness is `inert_classes(words)` (`SymbolicDynamics.py:441`), `entry.inert` None
  until the dynamics runs; `WalkStatus = resolved | ambiguous | unresolved` + six `WalkReason`s
  (`DualWalk.py:569-572`); `free_reduce` (`SymbolicDynamics.py:545`) on `cd.symbols`, closure and
  evidence, `translate_itinerary` unreduced; evidence comparison in `diagnostics.verify_words`;
  `PartitionFamily.owner_at_end` (`:342`) is the one element-at-end reader (B9b).
- "Topology never calls the map" has three exceptions at HEAD: `Arrangement.py:839`
  (`image_of` verification), `diagnostics.py:486` (`measured_sides`, on the default
  `punch_holes(diagnose=True)` path) and `plotting.py:364` (`hole_point`).
- Decision-id citations in `src/`: 313 regex hits on 264 lines in 30 files (271 in docstrings, 35
  in comments, 7 in LOG STRINGS: `TangleWorkbench.py:1528` NUM-3, `BridgeClass.py:620` D2,
  `MinimalTrellis.py:355` D9, `StablePartition.py:365` D9 and `:457` NUM-5,
  `diagnostics.py:714/:722` D6 -- these are code, a stripped-AST rubric sees them). Outside the
  regex: 71 slice / question / plan-section / experiment ids (`B1`, `Q13`, `A.8c`, `E11`, `F7`) on
  ~56 more lines, and 67 dated lines; 361 lines in 34 files carry some campaign reference.
- `topology/__init__.py`: 105 imported names, `__all__` 48, 56 imported names neither in `__all__`
  nor used (+ `annotations`). The names read through the namespace are 15, all in `__all__`
  (10 via `tanglepack/__init__.py`, plus `bridge_row_diagnosis`, `check_holes_share_bridge_side`,
  `letter`, `plot_stable_partition`, `side_after`); submodule imports (`DualWalk`, `diagnostics`,
  `plotting`, ...) are unaffected by a trim. BUT `MembersDisagree` (author-approved export, B6b Q5),
  `WalkReason` and `free_reduce` (B6b exports) are importable today only through those unused
  import lines: a trim to `__all__` un-exports them.

## N6 `readers-on-views` outcome audit (2026-10-26, numerics campaign, branch @ `e25adf1`)

- Six topology reader lines moved off nodes (no map call added): `StablePartition.py:835` reads
  `manifold.cdists[-1]`; `diagnostics._anchorward_look` bisects `manifold.cdists` and indexes
  `manifold.points`; `_bridge_side_of` sorts `bridge.cdists`; `plotting._stable_arc_midpoint` masks
  `points` by `cdists` (`if len(inside):`), `_image_arc_midpoint` uses `np.argmin` (same first-minimum
  rule). `Trellis.py:8,16` prose now says "the point store and its chains".
- Smell: `Trellis.cached_bridge_polyline` (`Trellis.py:571-604`) now caches a WRITABLE shared array (the
  bridge memo used to hand read-only arrays); `diagnostics._oriented_bridge_polyline`'s Note "must not be
  written to" is unenforced. Its docstring :578-579 ("the same invalidation as the point array
  underneath it") is stale: there is no point-array memo.


## N7 `delete-point` analyst note (2026-10-26, main `d279815`)

- Topology has NO code touchpoint for deleting `Point` / `BasePoint` (grep of `src/tanglepack/topology`:
  every `get_point()` is `Intersection.get_point()`; the six curve readers read `points` / `cdists` since N6).
- `Trellis.cached_bridge_polyline` (`Trellis.py:595-605`) shares a WRITABLE array since N6; one line
  (`value.flags.writeable = False` when not None) before the cache store restores the base's read-only
  guard: verified in the N7 prototype, the whole suite's failure set unchanged (no caller writes; readers
  `topology/diagnostics.py:149,201,343` only slice). Stale prose: `Trellis.py:578-579` ("the same
  invalidation as the point array underneath it").

