# Regions, Bridge Identity, Single Sources of Truth, and Cleanup — Implementation Plan

Date: 2026-09-02. Baseline: commit `df1afa9`, suite 123 passed / 1 skipped.

## Goal

Set up every piece of machinery needed to build a dual graph over the regions of a
tangle, without building the dual graph itself. Concretely:

1. Bridges keyed by their topological definition (ordered endpoint intersection pair).
2. Stable partition elements as identifiable objects, with bridge-to-element and
   element-to-image lookups, cross-trellis at the session level.
3. Enumeration of every minimal face bounded purely by manifold arcs, with lazy
   geometry (representative point, area, containment) and combinatorial adjacency
   and image lookups. No polygons are built unless asked for.
4. A growth driver in the style of `grow_until_turnaround` that grows until the
   requested crossings have computed iterates / the requested faces are closed.
5. One source of truth for every fact that currently has several, and a single
   generation counter that invalidates every derived cache.
6. The confirmed bugs from the 2026-09-02 audit fixed, with regression tests, and
   the dead code, style drift, and stale docs cleaned up.

Explicitly out of scope: the dual graph construction rules, any GUI, README.

## Definitions used throughout

- **Crossing / intersection**: one unstable x stable crossing, registry id `int`.
- **Anchor**: the periodic point itself, registered as a synthetic intersection at
  cdist `(0, 0)` on each (unstable branch, stable branch) pair that meets there.
  Today this happens by accident; after this plan it is deliberate.
- **Stable arc**: two consecutive crossings along one stable branch, ordered by
  stable cdist. Identified by `(lo_id, hi_id)`. This is the geometric edge; the
  partition elements that live on it are labels carried by the arc.
- **Unstable arc / bridge**: the piece of unstable manifold between two
  consecutive crossings along one unstable branch. `BridgeId = (first_id, second_id)`
  ordered by unstable cdist, i.e. in the unstable dynamical direction.
- **Face / region**: a cycle of arcs, traversable strictly along manifolds, that
  contains no other manifold piece. Any face whose boundary passes a node where a
  manifold continues but has no further computed crossing (a dangling end) is
  **open** and is not a region.
- **Generation**: monotone counter on the workbench, bumped by any mutation of the
  registry, the bridge set, or manifold geometry.

## Phase 0 — Guardrails and dead-code removal (do first, shrinks the refactor surface)

**0.1 Regression-marker semantics.** ✅ DONE (2026-09-02; `assert_area_preserved_along_chain` re-targeted to the registry iterate table because no BranchPoint iterate chains exist). `pyproject.toml` marker text becomes
"pins a previously fixed bug"; strip the three "this xfails" docstrings in
`tests/regression/`. Wire the two unused invariant helpers
(`assert_no_cdist_collision`, `assert_area_preserved_along_chain`) into
`tests/invariants.py` consumers.

**0.2 Delete verified-dead code** ✅ DONE (2026-09-02; `intersections_for_segment` had one script caller, script migrated; `Bridge.parent/children/next_bridge/prev_bridge` deferred to Phase 3 by decision so the suite stays green; `_linear_interpolation`/`_fractional_position` also deleted as they became caller-less) (zero callers in src/tests/scripts):
`ManifoldMachine.cut_manifold`, `_check_bridge_readiness`, `iterate_x_times`,
`_shift_list`; `TangleWorkbench.create_resonance_zone`, `index_manifolds`,
`populate_registry`; `IntersectionRegistry.find_by_cdist`;
`Tangle._same_unstable_branch`, `intersections_for_segment`, `_boundary_point`,
`_cache_boundary_preiterate` (and the one regression test that only exists to
exercise `_boundary_point`, replaced by an equivalent test through `create_bridges`);
`BaseManifold.plot_colormap`, `_iter_attr`; `ManifoldInitializer.get_all_initial_segments`;
the 58-line commented-out initializer body; `Bridge.next_bridge/prev_bridge`
(never read, wired wrong); `Bridge.parent/children` (replaced by derived queries in
Phase 3, deleted here so Blast is rewritten once).

**0.3 Two new invariant assertions** ✅ DONE as `check_holes_share_bridge_side` / `check_bridge_rows_consistent` in StablePartition.py, asserted in `Trellis.punch_holes` / `_hole_openings`. DEVIATION: wired as warnings until Phase 1.1 because period-3 violates both at baseline (origin (20,33) flips at iterate -2; bridge (11,33) rows differ); I1 takes an `orientation_preserving` flag (side alternates with iterate parity under det J < 0).
- All holes with the same `origin` share `bridge_side`.
- A bridge whose two crossings share a stable branch yields the same `row` at both ends.

## Phase 1 — Confirmed bug fixes (each with a regression test)

✅ DONE 2026-09-02 (all 20 rows; suite 211 passed / 1 skipped). Deviations: 1.2 uses branch-filtered bracketing with a logged nearest-cdist fallback (a 1e-6 relative match is impossible because root/tail are real nodes; the fallback is deleted by 2.2); 1.13 raises `ValueError` (input error) instead of asserting, and the Hénon sign flip turned out never to have existed in the repo history (only its comment did) so the `orient` hook is a new documented extension point; 1.14 needed NO code change — the plan text had the two ends swapped, `_build_intervals` already consults the existing neighbour on each side and the missing side is vacuously excluded (two pin tests added); 1.16 resets `pair.hole` for all recorded pairs, not just those in scope, because `Trellis.punch_holes` clears `self.holes` wholesale; 1.7's `FixedPoint.advance_key` landed here (plan 2.5). The Phase 0 invariants are now hard assertions in `Trellis.punch_holes` with `orientation_preserving` derived from det J.

| # | Location | Fix |
|---|---|---|
| 1.1 | `StablePartition._containing_bridge`, `_bridge_unstable_span` | Filter on `bridge.manifold_key`, never on endpoint keys. Root cause of the period>1 partition instability. Superseded again in Phase 2 when endpoint ids come from the cut itself, but fixed first so the topology suite is trustworthy. |
| 1.2 | `TangleWorkbench._assign_bridge_intersections` | Deleted in Phase 2 (endpoints set at cut time). Interim: branch filter + relative tolerance. |
| 1.3 | `Tangle._insert_segment` `_edge_seen` | Map edge to existing segment id; register that id under the new manifold instead of dropping it. |
| 1.4 | `Tangle._orientation`, `_find_true_intersection` | Scale both epsilons by the product of segment lengths; on a near-parallel pair log and discard instead of raising. |
| 1.5 | `Blast.blast_zone` except clause | Stop catching `AssertionError`. |
| 1.6 | `ManifoldInitializer` `has_inversion` | Call `check_inversion`. |
| 1.7 | `construct_kevin_way` vs `_advance_key_forward` | One (orbit, branch) advance rule, on `FixedPoint.advance_key` (Phase 2.5); both callers use it. |
| 1.8 | `ManifoldMachine.new_grow_manifold` | Single loop over `get_branch_array()`; drop the unused `iterated_manifold`. |
| 1.9 | `Intersection.fixed_points` | No `IndexError` when only `manifold_b_key` is set. |
| 1.10 | `Intersection.synthetic` | Keyword args; `label` no longer lands in `id`. |
| 1.11 | `trim_stable_manifolds` | Skip branches with no crossings; pass `"stable"` to `get_cdist`. |
| 1.12 | `grow_until_arclength` | Add `max_iterations`; rename `grown_until_intersection` to `grow_until_intersection` and give it `branch_index`. |
| 1.13 | `FixedPointSolver` | `fsolve(full_output=True)`, raise on `ier != 1`; assert real eigenvalues with one modulus above 1 and one below; move the Hénon-specific sign flip behind a documented `orient` hook with the eigenvector-orientation step as the intended mechanism. |
| 1.14 | `StablePartition._build_intervals` singleton loop | Consult `open_outward` at the last boundary and `open_anchorward` at the first, so the documented end pinches fire. |
| 1.15 | `StablePartition` region dedup key | Include `bridge_side`. |
| 1.16 | `Trellis.punch_holes` | Reset `pair.hole = None` for pairs in scope. |
| 1.17 | `ResonanceZone.restore` | Call `rebuild_bridges()`. |
| 1.18 | `Pseudoneighbor._strong_pip_cuts` | Warn when fewer cut points than `k_value`. |
| 1.19 | `IntersectionRegistry._get_lambda_u` | Honour the stability argument; return a float. |
| 1.20 | `StablePartition._stable_arc_midpoint`, `_near_far` | Chord fallback when the two bounds are on different stable branches. |

## Phase 2 — Single sources of truth (for review: each row is one change)

✅ DONE 2026-09-02 (rows 2.1, 2.2, 2.3, 2.5, 2.6, 2.7, 2.8, 2.9; 2.4 is Phase 3). Deviations: `Intersection` carries `unstable_manifold` (object) and `unstable_segment` (bracketing Points) because `create_bridges` must group by manifold OBJECT (a bridge and its parent share a key); a bridge bounded by fewer than two crossings is `Bridge.partial` (endpoint id None) and never dedupes — the blast script registers 141 instead of 140 bridges (a partial that used to get fabricated endpoints); `manifold_key` is keyword-required but the transient forward image inside `ManifoldMachine.iterate_manifold` is None until the workbench stamps the advanced key; `per_step_beta` uses the PERIOD root of the full-cycle eigenvalue (one branch return is k_value steps and costs λ^num_branches), and `per_step_factor` keeps the initializer's measured k_value-step ratio; `branch_cycle("stable")` is forward-map order for both stabilities (what StrongPip/Pseudoneighbor consume); `registry.graph` became a method taking `bridges=`; `IterateTable.items()` has no production caller yet; BranchPoints always allocate two eigendirection slots; `set_k_value` raises when the unstable and stable eigenvalue signs disagree. INVERSION FIXTURE DEVIATION: not a negative-b Hénon — with b<0 det J<0 so exactly one eigenvalue is negative at every fixed point and the two manifolds have different inversion status, which one k_value cannot represent; the fixture is the k=10, b=1 map's other fixed point (-2.3166, 2.3166), eigenvalues -4.406/-0.227, a genuine inversion saddle (k_value = 2) validated on all four manifold pieces.

| Fact | Today | After |
|---|---|---|
| 2.1 Resolved crossings | `Tangle._intersections`, `_intersecting_coords`, `_intersecting_points`, `_intersection_by_seg` **and** `IntersectionRegistry` | Registry only. `_resolve_crossing_pair` returns the `Intersection`; the workbench adds it. `create_bridges` consumes the registry grouped by unstable key. `grow_until_intersection` counts `len(registry)`. `trim_stable_manifolds` reads the registry. Tangle keeps only index state (`_intersecting_segments`, `_seg_lookup`, `_manifold_segs`). |
| 2.2 Bridge endpoints | Assigned after the fact by nearest-cdist lookup | Set by `Tangle.create_bridges` from the exact crossings it cut at. `first_intersection`/`second_intersection` typed `int`. `_assign_bridge_intersections` deleted. |
| 2.3 Bridge unstable branch | `manifold_key`, `_manifold_identity` fallback, two endpoint-inference paths | `Bridge.manifold_key` required at construction and propagated into every crossing born on that bridge (`_key_of` never returns None for a bridge). `_manifold_identity` fallback and keyless paths in `_containing_bridge` deleted. `ManifoldMachine.iterate_bridge` stops writing the un-advanced key. |
| 2.4 Bridge identity | Object identity, no id | `BridgeId` tuple (Phase 3). |
| 2.5 One map step | `_advance_key_forward`, kevin-way ordering, beta = lambda_u**(1/k) in four places, `Trellis.scale_cdist` in branch-return units | `FixedPoint.advance_key(key, n)`, `FixedPoint.per_step_beta(stability)`, `FixedPoint.branch_cycle(stability)`. `Trellis.scale_cdist` takes map steps. All callers switched. |
| 2.6 Manifold key | `BaseManifold.manifold_key`, dict key in `workbench.manifolds`, `Tangle._key_of` | `BaseManifold.__init__` requires the key; `workbench.manifolds` setter asserts equality; `_key_of` is a direct attribute read. |
| 2.7 Intersection graph | Registry graph (global-cdist adjacency, wrong for period>1) and workbench graph | Registry graph rebuilt from per-branch adjacency (stable: `TrellisBranch` order; unstable: bridge ids). `build_intersection_graph` decorates a copy of it. `IterateTable.items()` exposed; no private `_forward` access. |
| 2.8 Branch count | Hardcoded 2, default 1, derived from k | `FixedPoint.num_branches` property = 2 if inversion else 1; constructor argument removed. |
| 2.9 Iterate step for bridge endpoints | Heuristic `infer_iterate_table` on bridge boundary cdists | `iterate_bridge` maps the parent's two endpoint crossings forward explicitly and registers `(a -> f(a))`, `(b -> f(b))` in the `IterateTable` by collision on both cdists with tolerance. The heuristic remains for non-endpoint crossings. |

## Phase 3 — Bridge identity and derived genealogy

✅ DONE 2026-09-02. Deviations: `_existing_image_bridges` requires the existing bridges to tile the image exactly (first starts at f(a), last ends at f(b), consecutive share endpoints) rather than mere presence; `image_manifold` is not stored (iterate_bridge never kept it), so `iterated` is the only metadata carried by `rebuild_bridges`; a period-1 anchor at cdist (0,0) registers itself as its own forward image (`_register_anchor_self_iterate`) because the heuristic's "image is not its source" guard can never find it — to be absorbed by the deliberate anchor registration in 6.3; the k=10 graph gains one self-loop iterate edge; Blast needed no genealogy rewrite (it never read parent/children), only a BridgeId-keyed cycle guard; `Tangle.release_manifold` added for `clear_bridges`; `ResonanceZone.restore` now recomputes with `preserve_ids=True` (it renumbered the registry under every held id).

- `BridgeId = tuple[int, int]` ordered by unstable cdist. `Bridge.id` property.
- `TangleWorkbench._bridges: dict[BridgeId, Bridge]`; `bridges` still returns a list;
  new `bridge(bid)`, `bridges_at(intersection_id) -> list[BridgeId]` reverse index
  maintained in `create_bridges` / `clear_bridges`.
- `rebuild_bridges` carries per-id metadata (`iterated`, `image_manifold`) across
  by id; ids survive because registry ids survive `reindex_from`.
- Genealogy derived, not stored: `workbench.image_bridges(bid, n=1) -> list[BridgeId]`
  = bridges whose endpoints lie between `f^n(a)` and `f^n(b)` in the image branch
  order; `preimage_bridges(bid, n=1)` symmetric. Returns `None` when an endpoint has
  no registered iterate. Blast's frontier uses these.
- `clear_bridges` removes the discarded bridges' segments from the Tangle index.
- `_find_existing_bridge` / `_existing_image_bridges` become dict lookups by id.

## Phase 4 — Generation counter and caches

✅ DONE 2026-09-02 (suite 322 passed / 1 skipped; three smoke scripts byte-identical;
blast script 3.61s -> 3.59s). Deviations: `TangleWorkbench.generation` is a monotone
token derived from a composite state (own mutation counter, registry generation, the
tuple of manifold versions) rather than a plain sum, so a registry SWAP -- whose fresh
registry restarts its own counter at zero -- can never cancel out; the point-array memo
is enabled on `Bridge` ONLY (a growable manifold gains points through `Point`, which the
manifold cannot observe, and `ManifoldMachine`/`Point` were outside this phase's file
list), with every registered bridge's version bumped by the workbench on growth, a
re-cut and an iterate; memoised arrays are handed out READ-ONLY; the registry keeps a
second counter (`_orders_version`, bumped only by an insert or a renumbering) that the
ordering-derived views key on, so one `infer_iterates` pass does not rebuild them per
crossing; `bridge_for_pair` resolves through a trellis-local endpoint index
(`Trellis.bridge_between`) rather than `workbench.bridge`, because a trellis's bridge
list includes iterated children the workbench dict does not key the same way;
`_registry_id_epoch`/`_bridge_id_epoch` are KEPT and documented as orthogonal (they
answer "do ids still name the same crossings", not "is this snapshot stale");
`invalidate_trellises()` warns `DeprecationWarning` on the way to being removed;
`StrongPip`'s per-call `_cache` parameter is gone, folded into the `FixedPoint` memo.
20k synthetic inserts: 17.78s -> 0.13s. On the available fixtures `iterate_bridge`
never refines an already-cut arc (probed on k=10 including a 3-generation blast), so
its walk invalidation is pinned on the contract (every registered bridge's version
moves) rather than on a geometric side effect; growth does refine one (359 -> 360
points), and that is pinned directly.

- `IntersectionRegistry.generation`, bumped in `add`, `reindex_from`, `register_iterate`.
- `TangleWorkbench.generation` = own counter (bumped by bridge mutation, manifold
  root/tail reassignment, `manifolds` insertion) combined with the registry's.
- `Trellis._built_generation`; `TangleSession.trellis` staleness is one comparison.
  `_built_registry_size` and the identity check go away. Manual
  `invalidate_trellises()` calls become unnecessary (kept as a no-op alias for one release).
- `BaseManifold._version` bumped on any point insertion / tail change; `get_point_array`
  and friends memoise on it (bridges are immutable once cut, so this is free for them).
- Registry: maintained sorted key arrays plus an id-to-rank map; `_cdist_index`
  becomes the collision prefilter (or is deleted if the sorted arrays suffice).
  `_find_collision`, `_insert_into_*_order`, `nearest_by_unstable_cdist`, `*_rank`
  all leave O(N).
- `_register_forward_iterate`: candidates bucketed by branch key, bisect on cdist.
- Per-fixed-point memo of `_branch_position_map` on `FixedPoint`, invalidated by `set_k_value`.
- Per-trellis memo of `bridge_for_pair` (now a dict lookup by id) and oriented bridge polylines.

## Phase 5 — Partition elements

- `PartitionInterval` gains `element_id` (index within its result), `branch_key`, `side`.
- `StablePartitionResult` gains `element_of_intersection: dict[int, int]` and
  `elements_at_bridge: dict[BridgeId, tuple[Optional[int], Optional[int]]]`, built
  in `partition_stable_manifold` right after `_build_intervals`. Endpoints on a
  foreign branch are `None`.
- `Trellis.element_for(bridge_id, endpoint, side)` scans its own results.
- `TangleSession.partition_element_for(bridge_id, endpoint, side)` scans every
  cached trellis, closing the heteroclinic gap. `TangleSession` also gains
  `punch_holes`, `partition_stable_manifold`, `plot_stable_partition`,
  `describe_*` fan-outs, mirroring the pseudoneighbor ones.
- `StablePartitionResult.image_of_element(element_id, n=1)`: maps `lo_id`, `hi_id`
  through the iterate table to the arc on the image branch, then to the element(s)
  on that arc. Returns `None` if either iterate is missing.
- Result maps are keyed by intersection ids, so they survive everything but a
  reindex; `reindex_from` returns the remap needed to migrate them, applied in
  `Trellis` on generation change.

## Phase 6 — Geometry module, arrangement, regions

**6.1 `numerics/geometry.py`** consolidates the four primitives that exist today
in three places: `arc_polyline(manifold, lo, hi, reverse)` (from the resonance-zone
closure and the partition midpoint helper), `oriented_bridge_polyline` (moved from
StablePartition), `signed_polygon_area` (from `ResonanceZone.area`),
`point_in_polygon` (from `ResonanceZone.contains_point`), plus `polyline_midpoint`.
All work on the memoised point arrays from Phase 4.

**6.2 Crossing sign.** `Intersection.crossing_sign: int` = sign of the cross
product of the unstable segment direction and the stable segment direction, both
taken in increasing cdist, computed in `_resolve_crossing_pair`. Anchors get it
from the oriented eigenvectors.

**6.3 Anchors.** `compute_intersections` registers one synthetic anchor crossing
per (unstable branch, stable branch) pair meeting at each periodic point, via
`Intersection.synthetic`, keyed like any crossing. The accidental cdist-0 crossing
is removed by the same collision test.

**6.4 `topology/Arrangement.py`.** Built from the all-fixed-points trellis
(`session.trellis(None)`), cached by generation.
- Nodes: registry ids including anchors.
- Half-edges per node, four slots: `s-` (stable toward anchor), `s+`, `u-`, `u+`.
  Stable arcs from consecutive ids in `TrellisBranch.ordered_ids`; unstable arcs are
  exactly the bridge ids (asserted against the per-branch order).
- A missing slot is filled with a **virtual half-edge** to a virtual tail node.
- Rotation at a node from `crossing_sign` alone: `(u+, s+, u-, s-)` for positive
  sign, `(u+, s-, u-, s+)` for negative. No angle sorting.
- Face traversal: next half-edge = the one after the reverse of the incoming
  half-edge in the rotation. Linear in the number of arcs, no networkx.
- Faces containing a virtual half-edge are `open`; the rest are regions.

**6.5 `Region`** (`topology/TopologyResults.py`): `corners: tuple[int, ...]`
canonicalised to start at the smallest id with fixed orientation, `arcs: list[Arc]`,
`is_closed`, `stable_arcs`, `bridge_ids`, `neighbor_across(arc)`, and lazy
`boundary_points`, `area`, `representative_point` (mean of arc midpoints;
`verify_representative_point()` runs the ray cast on demand), `contains(point)`.
`Arc` = `(kind, lo_id, hi_id, branch_key, bridge_id | None, reverse)`; a stable
`Arc` reports the partition elements on it via the Phase 5 maps.

**6.6 Combinatorial lookups on `Arrangement`**: `region(corners)`,
`regions_at(intersection_id)`, `regions_bounded_by(arc)`, `image_of(region, n)`
(corners through the iterate table, `None` if any is missing or the image face is
open), `preimage_of`. These are what the dual graph will be built from.

**6.7 `ResonanceZone` on `Region`.** `_build_boundary` returns `list[Arc]`;
`area`/`contains_point` delegate to geometry; `boundary_intersection_id` dropped in
favour of resolving on demand; trimming, `restore`, `previous_tails` stay on the zone.
`Blast.in_zone` and `classify_bridge(s)` use the region `contains` with the memoised
bridge test point. `shapely` removed from `pyproject.toml`; the stale "cached shapely
polygon" docstring deleted.

## Phase 7 — Growth driver

`TangleWorkbench.grow_until(fixed_point, predicate, *, grow=("unstable","stable"),
max_iterations=10, branch_index=0)`: loop grow one iteration on the listed
stabilities, `compute_intersections(preserve_ids=True)`, `infer_iterates`, stop when
`predicate(workbench)` is true; raise `ValueError` at the cap like the siblings.
Two named wrappers:
- `grow_until_iterates_closed(fixed_point, ids="all", direction="forward", ...)`:
  forward images need the unstable manifold grown, backward the stable; the wrapper
  picks the stability. Predicate: every id in the set has an `n=±1` entry.
- `grow_until_faces_closed(fixed_point, ids, ...)`: predicate: every face incident
  to the given crossings is closed in the arrangement.
Both exposed on `TangleSession`.

## Phase 8 — Cleanup

- Duplication: one `_collect` walker in `BaseManifold` (replaces five getters and
  the string-dispatched `_iter_method`); one branch walker; `_require_manifold`
  guard; `_iter_manifolds` used everywhere; one `_fanout_plot` in `TangleSession`;
  `tanglepack.examples.henon` used by scripts, tests, and visualisations;
  `Stability`/`Side` aliases in one place.
- Splits: `TangleWorkbench` into orchestrator + `graphviz.py` + `BridgeIterator.py`
  + `IterateInference.py`; `StablePartition` into partition logic +
  `topology/plotting.py` (shared with Trellis plotters); geometry already moved.
- Style: `from __future__ import annotations` everywhere; type hints on all public
  numerics signatures; class docstrings for `TangleWorkbench` and `ManifoldMachine`;
  Dev Notes moved above imports as module docstrings; inline TODOs into Dev Notes;
  `Trellis` `verbose` prints removed; `Tangle` gets a `NullHandler`; the stale
  "should be fixed now" WARNING block removed; the contradictory refinement notes reconciled.
- Repo: untrack `notebooks/nothing.txt`, `Untitled-1.ipynb`; ignore `profile.out`;
  delete the four dead notebooks or move them to `notebooks/archive/`; fix the one
  stale-API notebook; re-run the three scripts-folder notebooks to clear stale outputs.
- Docs: CLAUDE.md table rewritten for the three subpackages with the new modules
  and the `BridgeId` / `Region` / generation concepts; the two partial plans in
  `docs/` annotated with what was implemented and what was abandoned; Sphinx stubs
  regenerated.

## Tests added (per phase)

- P1: one regression test per row of the table, plus the two Phase 0 invariants
  run on the period-3 fixture.
- P2: registry-only bridge cutting agrees with today's on the k=10 and p3 fixtures;
  a keyed-crossing test on iterated bridges; `advance_key` round-trips for period 1,
  3, and an inversion fixture (new: a map with a negative-eigenvalue saddle, so the
  `k_value = 2*period` path is exercised for the first time).
- P3: `BridgeId` stable across `rebuild_bridges`; `image_bridges` matches a
  hand-checked blast on p3.
- P4: generation bumps on every mutation path; trellis cache hit/miss table; registry
  insert is linear (timing bound on 20k synthetic crossings).
- P5: element maps on p3 including the heteroclinic endpoints via the session lookup;
  `image_of_element` against the iterate table.
- P6: face count and closedness on a hand-drawn arrangement, on k=10 period 1, and on
  p3; every closed region's representative point passes `contains`; region image
  agrees with mapping the representative point forward and locating it by ray cast;
  ResonanceZone area and containment unchanged from today.
- P7: drivers terminate on the fixtures and raise at the cap.

## Order and checkpoints

0 -> 1 -> 2 -> 3 -> 4 -> 5 -> 6 -> 7 -> 8. Phases 0 through 2 are a coherent
correctness checkpoint (suite green, topology tests trustworthy). Phases 3 through
5 are the identity checkpoint. Phase 6 and 7 deliver the region machinery. Phase 8
is mechanical and can be interleaved once the interfaces have settled.

## Deferred (found during implementation, not in the plan)

- ~~Phase 1: `clear_bridges`/`rebuild_bridges` leave discarded bridges as segment owners~~ fixed in Phase 3.
- Phase 1: `_insert_crossing_separator` leaves the parent's pre-split edge registered with stale bounds (pre-existing).
- Phase 1: `_validate_saddle` may reject a near-parabolic period-k saddle (complex pair from rounding) and `ier != 1` rejects `ier=5` at machine-precision residuals; gate on the residual if it ever bites.
- Phase 1: style items for Phase 8 — `TangleWorkbench.py` lacks `from __future__ import annotations` (PEP 604 return annotation in `_endpoint_candidates`), `Tangle` logger has no `NullHandler`, unannotated params in `Tangle._orientation` / `ManifoldMachine._branch_view`, `Tangle._segments_touching` and `StablePartition._bridge_midpoint` are caller-less, `FixedPoint.get_iterable_array(shift=)` has no production caller.
- Phase 1: p3 registry ids are not deterministic across two builds in one process (id()-ordered iteration somewhere in `Tangle.create_bridges`/`_intersecting_segments`); no test hard-codes ids, but worth pinning down in Phase 2/4.
- Notebooks `scripts/tangle_workbench_test.ipynb` and `scripts/tangle_workbench_higher_period_test.ipynb` call the renamed `grown_until_intersection`; Phase 8 notebook pass.
- `TangleSession.blast_zone`'s `strict` docstring should say AssertionError always propagates; a `TangleSession.restore_resonance_zone` wrapper invalidating trellises would close a staleness gap (moot after Phase 4).

- Notebooks `scripts/tangle_workbench_test.ipynb` (`create_resonance_zone`), `scripts/manifold_cutting_test.ipynb` (`intersections_for_segment`) call names deleted in Phase 0; fix or archive with the Phase 8 notebook pass.
- `docs/nested_tangle_bridge_cutting_fix_plan.md` and `docs/higher_period_bridge_bug_plan.md` describe `_boundary_point` / `index_manifolds` as current; annotate in the Phase 8 docs pass.
- Pre-existing unused imports in TangleWorkbench.py (`Annotated`, `Callable`, `npt`) and ManifoldInitializer.py (`Optional`); Phase 8 style pass.
- `Tangle._cdist_at_point` / `_cdist_between` near-duplicates; merge with Phase 6.1 geometry consolidation.
- Phase 2: `Pseudoneighbor.py` still uses bare `lambda_u` as the one-full-return factor in two places (wrong for inversion: the factor is beta**k_value); the module says inversion is unvalidated there.
- Phase 2: `orient_eigenvectors` is a no-op on inversion points (pre-existing); fine for period 1 but the per-orbit-point eigenvector propagation is skipped for period>1 inversion points.
- Phase 2: `grow_n_times` defaults `branch_index=0` so every call on an inversion point logs the ignored-branch warning; growing the second eigendirection of a simple saddle is untested; `Trellis.scale_cdist` has no production caller; `registry.by_*_cdist`/`*_rank` are still global-order views (Phase 6 wants per-branch).
- Phase 2: `tests/visualizations/viz_nested_bridges.py` and `viz_unstable_artifacts.py` crash with `KeyError: None` (strong pip None on their shorter growth) — pre-existing, Phase 8.
- Phase 2: the `_insert_crossing_separator`-related reused-edge crossings are registered but cut only on their own curve (Dev Note in TangleWorkbench).
- Phase 3: `rebuild_bridges(fp)` calls the global `clear_bridges()`, dropping every other fixed point's bridges (pre-existing); partial bridges still accumulate one object per re-blast of the same arc (no identity for arcs with one crossing).
- Phase 3: INVERSION (found by the new `henon_inversion_initialized` bridge test, both pre-existing and both for 6.3/the inversion work, not Phase 3): the periodic point registers FOUR crossings at cdist (0, 0) — one per (unstable branch, stable branch) pair — so each unstable branch carries two distinct anchor crossings and `Tangle.create_bridges` cuts a degenerate zero-length bridge between them; and `iterate_bridge` on such an anchor-rooted bridge raises `TypeError: list indices must be integers or slices, not NoneType` in `ManifoldMachine._insert_after` (`branch_index=None` on a BranchPoint). Every non-degenerate inversion bridge iterates and derives its image correctly.
- Phase 3: id validity is tracked by `TangleWorkbench._registry_id_epoch` (bumped by a `compute_intersections` that renumbers) and `_bridge_id_epoch` (its value when the bridges were cut); `rebuild_bridges` refuses to carry metadata across a mismatch. Phase 4's generation counter should absorb both.
- Phase 3: on the inversion fixture the periodic point registers four cdist-(0,0) crossings (one per (unstable branch, stable branch) pair), so each unstable branch carries two anchors and `create_bridges` cuts a degenerate zero-length bridge between them; iterating such an anchor-rooted bridge raises in `ManifoldMachine._insert_after` (branch_index None on a BranchPoint). Both belong to 6.3's deliberate anchor registration.
