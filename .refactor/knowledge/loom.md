# loom — knowledge file

Last updated 2026-10-02 (auditor, outcome of `hole-side-and-own-blast`, branch at `3b33baa`).

## Map

- `loom/TangleSession.py` (~2000 lines) — the facade. Trellis cache per fixed-point selection
  keyed by `workbench.generation` (`trellis` :158); downstream caches (`bridge_classes`,
  `minimal_trellis`, `iterated_partition`, `dual_graph` :618, `symbolic_dynamics` :721) keyed on
  `(generation, partition signature[, gathered strong pips])`; zones (`resonance_zone` :1737,
  `add_resonance_zones` :1772, `classify_bridge` :1825 = innermost containing zone by area);
  `blast_zone` :1871 is a pure forward to `Blast.blast_zone`; `cartoon_zones` :1061.
- `loom/Blast.py` — `blast_zone(session, zone, n, fixed_point=, strict=, min_separation=)`,
  `BlastResult`/`BlastStep` genealogy, `_seen_key` (BridgeId or id()), proximity guard
  (`_interior_points`, cKDTree on every existing bridge interior), `_resolve_zone`.
- `loom/ResonanceZone.py` — `ResonanceZone` (cut intersections kept as `Intersection` objects,
  frozen `boundary_vertices` polygon, `contains_point`, `restore`), `define_resonance_zone`
  (trim each stable branch at pip + k-1 iterates, recompute with `preserve_ids`, rebuild bridges,
  capture boundary), `_anchor_id`/`_live_id` = `registry.find(...)` against the live registry.
- `examples/henon_cases.py` — `build_period3`, `build_nested`, `TangleBuild`; helpers `_repin`
  (classify + set_strong_pip per fp) and `_partition` (pseudoneighbors, holes, partition).
  Since 2026-10-02 `build_nested` blasts each zone with `fixed_point=[zone.fixed_point]` (:264),
  so blast order does not matter; default `blast_sizes` is `[3, 3]` (was `[32, 7]`), p3 bridges 27
  (was 31), words unchanged.

## How it works

- Blast step: frontier = un-iterated bridges passing `fp_ok` (fixed-point filter) and `in_zone`
  (bridge midpoint inside the zone's FROZEN polygon, boundary inclusive). Each parent ->
  `workbench.iterate_bridge`; children already `seen` are skipped; children out of zone are
  discarded; children within `min_separation` of any existing bridge interior are dropped.
  `AttributeError`/`ValueError` from the forward map are tolerated (WARNING, counted `skipped`).
- Containment uses `zone.contains_point` only — NOT `classify_bridge`'s innermost-zone rule. An
  outer-zone blast therefore includes every bridge inside the inner zone (and in the annulus).
- Any registry mutation moves the generation; the rebuilt trellis has empty strong-pip slots, hence
  `_repin` after blasts. Intermediate `_repin`s between blasts do not change the result (verified
  2026-10-02, nested 4+4 own-fp: identical classes and words); only the final one matters.

## Invariants relied on

- The zone polygon is frozen at definition; later trims/recomputes do not move it.
- Registry ids survive `compute_intersections(preserve_ids=True)`; zones never store ids.

## Smells

- FIXED 2026-10-02: `build_nested` blasted both zones with both fixed points, so an outer blast
  iterated p3's exterior (order dependence, a `merge_manifolds` crash logged as "skipped", p3
  pseudoneighbor collapse at >= 6 outer blasts). Guarded by `test_nested_blast_order_does_not_matter`,
  `test_nested_inner_words_are_the_period3_words`, `test_nested_outer_blasts_leave_the_inner_bridges_alone`.
  The underlying `Blast` behaviour (containment by polygon only, any fixed point the caller passes)
  is unchanged: a caller passing both fixed points to an outer blast still gets the old effect.
- `henon_cases.py:282` constant slug `nested_p1_p3` whatever the blast counts (figure overwrite).
- `henon_cases.py:142-146,192,260` `_repin` after every blast is dead work; the real issue is that
  the session forgets the pip choice on every generation bump (design: trellis slots start empty).
- `Blast.py:226-236` hand-rolled fixed-point normalisation duplicates
  `TangleSession._resolve_fixed_points`. `Blast.py:289` imports cKDTree inside the loop although
  `_min_interior_distance` already imports it (with a scipy-missing fallback that cannot happen).
- `BlastResult.fixed_point` stores the raw filter argument (FixedPoint, list or None) — typed
  `Optional[FixedPoint]`, wrong for lists.

## Features

- Blasting a single zone, repeated single-iteration blasts: p1 alone 0..8 and p3 alone 0..5
  reach reliable symbolic dynamics (verified by probe, 2026-10-02); p3 >= 6 runs but is
  unreliable, see topology.md.
- Nested (branch, verified by the 2026-10-02 audit probe): all 16 cells outer {0,2,4,8} x inner
  {0,2,4,8} run; outer words = p1 alone and inner words = p3 alone at the same counts in every
  cell; outer-first = inner-first for (2,2), (4,4), (8,4). Reliable iff outer >= 2 and inner <= 4.
- Resonance zone define/restore: `tests/test_loom_blast_restore.py`, `test_resonance_zone_region.py`.
- Cost: blasting dominates (p1 8 blasts ~4 s; nested 8+8 ~5 s); every topology stage after it is
  < 0.2 s even at 400 crossings. Inside `iterate_bridge`: `Tangle.add_manifold` segment
  bookkeeping + rtree queries ~50 %, `refine_manifold` ~30 % (cProfile, p1 8 blasts).
