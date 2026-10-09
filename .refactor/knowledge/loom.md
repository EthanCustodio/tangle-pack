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

## Open-issues round WP3 (2026-10-08, branch `oi/wp03-inversion`)

- `ResonanceZone.ring_cycles(cuts)`: from a cut's stable branch `s(o,b)` to the cut on
  `u(o,b)` (FULL key, not orbit index). A period-k orbit closes one cycle; an inversion point
  one loop per branch pair, both through its anchor; a run reaching a branch with no cut stops
  open. `_build_boundary` emits each cycle `u` arc out, `s` arc back; non-inversion polygons are
  the same cycle (areas equal to 12 digits; the ring starts at an anchor instead of a pip).
  `TangleSession._zone_ring_order` reuses it (empty when a cycle does not close).
- One zone per fixed point (C8): `ResonanceZone.key` is the fixed point,
  `session.resonance_zones[fp]`, `_blast_stats[fp]`, `zone_label` "zone p<period>";
  `ResonanceZone.branch_index` deleted. Inversion zone area 1.412087 (was 1.232), which changes
  the bridges the inversion blasts iterate (classes 4 / 10 / 10, was 6 / 11 / 11).

## Open-issues round WP6 (2026-10-08, branch `oi/wp06-minimal-alphabet`)

- `TangleSession.minimal_trellis` builds the EVIDENCE table itself
  (`_bridge_classes(trellis, gathered partitions)`, after `_warn_unpartitioned_branches`) and
  passes it to `minimal_trellis`; `bridge_classes()` is `minimal_trellis().alphabet`, lettered
  (`BridgeAlphabet`) and zoned, cached as `(minimal, table)` per selection key and rebuilt
  whenever the minimal trellis object changes (so it can never hold entries of a stale minimal
  trellis). `symbolic_dynamics` takes `bridge_classes(fixed_points)` after the dual graph (no
  second rebuild of the minimal trellis).
- `BridgeAlphabet.letter_for` uses `topology.SymbolicDynamics.letter` (active `a..t, aa, ...`);
  `loom.letter` is gone (no shim, D8). Session caches stay hand-rolled (SIMP-8, WP8).

## Open-issues round WP7 (2026-10-08, branch `oi/wp07-walk-exists`)

- `TangleSession.symbolic_dynamics` passes `letter_for=self.bridge_alphabet.letter_for`, so an
  active class the closure adds (a class a word names outside `bridge_classes()`, C9) keeps its
  letter across rebuilds; `bridge_classes()` stays the seed alphabet, `dyn.table` is closed.

## Open-issues round WP8 (2026-10-08, branch `oi/wp08-recompute-identity`, base `083f910`)

- `loom/Topology.py`: `Topology(session, fixed_points)` = one recompute (D1, SIMP-8). Its
  constructor runs `Topology._stages_3_to_6(fp)` per selected fixed point (classify when no pip,
  pseudoneighbors, holes, partition, on the per-fixed-point trellis, from scratch) and copies the
  holes / partitions into its own lists (a later manual stage call clears the trellis slots in
  place). Lazy `cached_property` stages: `homotopy`, `minimal`, `classes` (alphabet lettered by
  `identity`, zoned by `Topology._zone_key` through the public `session.classify_bridge`), `iterated`, `dual_graph` (pips with
  None filtered), `dynamics` (`letter_for=self.letter_for`). `identity(bridge_class)` =
  frozenset of the two elements' `PartitionInterval.identity`.
- `TangleSession`: `_trellises` / `_topologies` keyed by `frozenset(id(fp))` of the resolved
  selection (`_selection(fixed_points)` returns the key and the fixed points in workbench order; `trellis(fp) is trellis([fp])`, a one-fixed-point
  session's `trellis() is trellis(fp)`), the trellis built with the selection in WORKBENCH order.
  `topology()` rebuilds on `rebuild`, a new generation, or `topology.strong_pips !=
  tuple(self.strong_pip(fp) ...)`. The six product accessors are one-line delegates.
- Pip memory (API-1): `_strong_pips[fp] = Intersection`, written by `set_strong_pip`,
  `resonance_zone` (the zone's pip, read from the registry BEFORE the trim; `add_resonance_zones`
  goes through it), `classify_strong_pips` (fan-out: re-applies a remembered pip with `choose_default=False`, so a
  lost one stays unset after ONE WARNING and is never replaced by the default; else remembers the
  default) and `Topology._stages_3_to_6` (a pip set directly on the trellis is adopted through
  `set_strong_pip`; a lost chosen pip, read by the public `chosen_strong_pip(fp)`, raises). A rebuilt
  single-fixed-point trellis is classified (`choose_default=False`) and the remembered pip found
  by `registry.find`; not a candidate -> WARNING, left unset, `topology()` raises ValueError.
- Deleted: `_bridge_classes`, `_dual_graphs`, `_minimal_trellises`, `_iterated_partitions`,
  `_symbolic_dynamics`, `_partition_signature`, `_gathered_partitions`, `_gathered_holes`,
  `_gathered_strong_pips`, `_annotate_bridge_classes`, `_warn_unpartitioned_branches` (dead: every
  selected branch is partitioned). Kept: `partition_element_for` (scans snapshots; WP11).
- `henon_cases`: `_repin` / `_partition` gone; builders end at `session.topology()`. The E07 / E11
  harness lost its `hc._repin` calls; E11's `repin()` helper still runs and now returns the
  restored pips, not the defaults (WP12 rerun). `scripts/henon_bridge_classes.py` ends at
  `session.topology()`. `scripts/refactor_proof/hole-side-and-own-blast/*` still import the
  deleted `hc._repin` / `hc._partition` (historical proof scripts, not maintained).
- A `rebuild=True` on any product recomputes the holes too, whose geometric cross-check calls the
  map (escaping builds need `np.errstate(over="ignore")` around it).

## Open-issues round WP8b (2026-10-08, branch `oi/wp08b-d10-d12`, base `014ea4b`)

- `Topology.__init__` runs stages 3-6 per fixed point FIRST (each per-fixed-point trellis is cut
  at its pip by `set_strong_pip`; a default pip is chosen there), then takes the selection
  trellis. `TangleSession.trellis` cuts a selection of SEVERAL fixed points at their pips
  (`self.trellis(fp).strong_pip`) on EVERY read, a no-op while they are unchanged; a
  per-fixed-point trellis is cut as soon as its pip is set or restored (rebuild path). Review fix:
  the cut used to live in `Topology.__init__`, so `session.trellis()` of several fixed points
  showed the computed N until a product ran (nested_outer_series after 10 outer blasts: 736
  crossings / 525 bridges before `topology()`, 448 / 322 after; probe
  `/home/dezhu/.claude/jobs/c5c98371/tmp/wp08b_review/hist.py`).
- A class member can be a JOINED arc of the cut (not a workbench bridge): `Topology._zone_key`
  and `TangleSession.cartoon_enclosures` resolve bridges through `trellis.bridge_between`
  (`cartoon_enclosures(trellis, bridge_ids, parents)`: the plotted dual graph's own trellis,
  `graph.trellis`, so a missing bridge asserts; it used to read `self.trellis()` and skip)
  (`workbench.bridge` would raise KeyError).
- Choosing an inner pip candidate on k10 can now raise the I2 row assertion in `topology()`
  (candidate 2: the anchor bridge wraps the end of the cut stable branch), exactly as a
  numerical trim at that pip does. Tests choose the anchor-nearest alternate.


## Open-issues round WP9 (2026-10-08, branch `oi/wp09-heteroclinic`, base `b460363`)

- `Topology.__init__`: stage 3 per fixed point (`_choose_pip(fp)`, on its own trellis; same
  ValueError text when a chosen pip is lost), then `session.trellis(selection)` (cut at every pip)
  and stages 4-6 ON IT; `holes` / `partitions` are copies of that trellis's slots. For one fixed
  point the selection trellis IS its trellis (`trellis(fp) is trellis([fp])`), so nothing
  changes there. After `topology()` of several fixed points the per-fixed-point trellises carry
  no stage output; `topology(fp)` fills one (tests: `conftest.p3_partitioned`). `topology()` and
  `topology(fp)` differ exactly where heteroclinic pairs exist.
- C4 (`Blast.py`): `in_zone(b)` = the blasted zone contains b's midpoint and no SMALLER stored
  zone of ANOTHER fixed point does (innermost containing zone), for the frontier and the kept
  children. The first cut, `classify_bridge(b) is zone`, blasted nothing for a zone the session
  does not store (`define_resonance_zone`, a replaced or copied handle; review). Unfiltered nested blasting is order independent and
  equals the own-fixed-point recipe; an unfiltered outer blast still iterates the inner tangle's
  annulus bridges (its own by C4), one of which can hit NUM-9 (WP10).
- `henon_cases`: `_orbit(session, guess, orientation)` + `_grow(session, fp, u, s, cutoff)`
  (`s=None` = turnaround) replace `_grow_period3/_grow_period1`; `build_nested` constructs p1
  THEN p3 (NEST-6: p1 = `A` = tangle 0, the order of `TangleBuild.fixed_points`), grows p3 then
  p1 (unchanged cutoffs, session ends at p1's), `compute_intersections([fp1, fp3])`.
  `build_heteroclinic(**overrides)` = `build_nested(p3_unstable_steps=16, p3_stable_steps=10)`,
  slug `heteroclinic_p1_p3` (scan in its Dev Notes).
- `scripts/symbolic_figures.py::draw_trellis` draws the topology trellis's holes
  (`trellis.plot_holes`), not the per-fixed-point fan-out (which is empty for nested now).
- Later (D8 API plan): defaulting `blast_zone`'s filter to the zone's own fixed point;
  `partition_element_for` could read the selection trellis directly (WP11).

## Open-issues round WP10 slice 1 (2026-10-08, branch `oi/wp10-numerics-leftovers`, base `3f27456`)

- `Blast.blast_zone`: no `strict`, no `BlastStep.skipped` / `BlastResult.skipped`; any iterate
  failure propagates (NUM-10). The growth factor divides by `frontier_points` (nothing skipped).
  It collects the bridges it ADDED (the `held` bookkeeping) and ends with
  `report_resolution(allowed_fps, new_bridges=added)` (C5).
- `TangleSession.trellis` no longer classifies on a rebuild: `_restore_strong_pip` finds the
  remembered crossing and tests it alone (`Trellis._is_candidate`), then sets `trellis.strong_pip`
  directly (set_strong_pip would test it a second time). `Topology._choose_pip` ->
  `session.set_strong_pip` -> `Trellis.set_strong_pip` skips validation when the id is already the
  chosen pip. One `is_strong_pip` call per remembered pip per recompute
  (`test_a_rebuild_tests_only_the_remembered_pip`).
