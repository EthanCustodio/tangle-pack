# loom — knowledge file

Last updated 2026-10-09 (open-issues round WP12, branch `oi/wp12-reverify`, base `a2e0e65`; line
numbers on that branch). The dated sections below are each WP's notes, kept as the record.

## Map

- `loom/TangleSession.py` (1735 lines) — the facade. One trellis snapshot per fixed-point SET
  keyed by `workbench.generation` (`trellis` :200; a selection of several is cut at their pips on
  every read); `topology` :289 = one `loom.Topology` per `(generation, selection, pips)`, every
  product accessor a delegate; pip memory (`classify_strong_pips` :1048, `_restore_strong_pip`
  :1092, `set_strong_pip` :1109); zones (`resonance_zone` :1474, `add_resonance_zones` :1512,
  `classify_bridge` :1565 = innermost containing zone); `blast_zone` :1611 forwards to
  `Blast.blast_zone`; `cartoon_zones` :780.
- `loom/Topology.py` (205 lines) — one recompute (D1): stage 3 per fixed point, stages 4-6 on the
  selection's trellis, lazy `homotopy`, `minimal`, `classes`, `iterated`, `dual_graph`,
  `dynamics`; `identity` / `letter_for` (C10).
- `loom/Blast.py` — `blast_zone(session, zone, n, fixed_point=, min_separation=, diagnose=,
  max_points=, predict=)` :185, `BlastResult` :134 (`stop_reason`, `points_created`,
  `fixed_points`), `in_zone` :271 (innermost zone, C4), proximity guard (`_interior_points` :85,
  `_min_interior_distance` :100), `_seen_key` :72.
- `loom/ResonanceZone.py` — `ResonanceZone` (cut crossings kept as `Intersection` objects, frozen
  polygon, `contains_point`, `restore` :328), `ring_cycles` :383 / `ring_order` :417 (full keys,
  C7), `define_resonance_zone` :498 (trim each stable branch at the pip orbit with
  `workbench.trim_stable_at`, recompute with `preserve_ids`, rebuild bridges).
- `examples/henon_cases.py` — `TangleBuild`; `build_period3` :177 (slug `_period3_slug` :226),
  `build_nested` :250 (`_nested_slug` :239), `build_heteroclinic` :327, `build_inversion` :350
  (WP12: the law-case recipe; raises when f(q0) is not registered); helpers `_session`, `_orbit`,
  `_grow`, `_intersect_and_bridge`. Every builder ends at `session.topology()`; nothing repins.

## How it works

- Blast step: frontier = un-iterated bridges of the filter's fixed points whose innermost zone is
  the blasted one (C4). Each parent -> `workbench.iterate_bridge`; children already `seen` are
  skipped, children out of zone discarded, children within `min_separation` of an existing bridge
  interior dropped; any iterate failure propagates (NUM-10). The image curves stay computed
  manifold (C20): every later recompute keeps the blast.
- The session carries each fixed point's pip as an `Intersection` across generations (API-1);
  a rebuilt trellis finds it by `registry.find`, tests it alone, and is cut at its orbit (D10).
- A zone's trim is a VIEW (C22, WP12): `grow_n_times` re-finds every branch's physical end, so
  stable growth un-trims the zone; the topology does not care (the trellis is cut at the pip
  orbit), and `ResonanceZone.restore` never shortens a branch growth has run past (before WP12
  it cut k=2.8 two blasts back from 362.9 to the pre-zone tail 42.25). A zone that re-trims
  itself after growth is D8 API work.

## Invariants relied on

- The zone polygon is frozen at definition; later trims/recomputes do not move it.
- Registry ids survive `compute_intersections(preserve_ids=True)` and are reproducible (WP8);
  zones never store ids.

## Smells

- Fixed since the 2026-10-02 audit: the constant nested slug (API-7, WP11), `_repin` after every
  blast (API-1, WP8), the hand-rolled fixed-point normalisation and the in-loop cKDTree import in
  `Blast.py` (SIMP-7, WP11), `BlastResult.fixed_point` (now `fixed_points`, WP11), the
  polygon-only blast containment (C4, WP9), the period-3 slug collision (WP12).
- PERF-2 (WP12, accepted): blasting past the freeze of N costs most of the time. k28 at 1e-7:
  blasts 10-13 take 88 % of 13 blasts' blasting time (x1.67 per blast, 3.6 GB peak) while N stays
  at 422 from blast 9; k32 the same from blast 9 (4.9 GB at 13). A "blast until T(q0) stops
  growing" driver, or blasting only T(q0)'s bridges, belongs to the D8 API plan.
- Open for the D8 API plan: `blast_zone`'s filter could default to the zone's own fixed point;
  `partition_element_for` could read the topology instead of scanning snapshots.

## Features

- Blasting a single zone (re-measured by E12, WP12, at 1e-7): k28 to 13 blasts reliable from
  blast 1 (C21); k32 reliable to blast 8, then INT-11 kind 2 (resolution floor); p3 to 20 blasts
  reliable throughout, the alphabet plateaus at 25 from blast 8; inversion drops at blast 6
  (BCG-8). (The 2026-10-02 notes said p3 >= 6 blasts was unreliable: fixed by WP4-WP8b.)
- Nested: nested_outer 12 -> 13 classes, flat to 12 blasts, reliable from blast 1; nested_inner
  (2 outer blasts, then the inner zone) 13 -> 29 at blasts 6-8, flat to 12, reliable throughout.
  (The old "reliable iff outer >= 2 and inner <= 4" no longer holds: every class walks since
  WP7 / WP8b.)
- Resonance zone define/restore: `tests/unit/loom/test_blast_restore.py`, `tests/unit/loom/test_resonance_zone.py`.
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

## Open-issues round WP11 slice 1 (2026-10-09, branch `oi/wp11-simplify`, base `4f7e815`)

- `Blast.py`: `fixed_points = session._resolve_fixed_points(fixed_point)` feeds `fp_ok` and
  `report_resolution`; `BlastResult.fixed_points` (list) replaces the unread `fixed_point`;
  `cKDTree` / `ResonanceZone` imported at module top; `_resolve_zone` and `max_depth_reached` gone.
- `ResonanceZone`: `_trim_stable_at` (rtree internals + `cdist_tol`) replaced by
  `workbench.trim_stable_at`; `_live_id` inlined as `registry.find`; `ring_order(cuts)` (was the
  session's `_zone_ring_order`) sits beside `ring_cycles`.
- `TangleSession`: `invalidate_trellises` deleted; coordinate readers stop ravelling (API-6).
- `henon_cases._nested_slug(outer, inner)`: non-default blast counts in the slug (API-7).
- Deferred: `partition_element_for` reading the topology instead of scanning snapshots (changes a
  public answer, D8 API plan).

## Open-issues round WP12 (2026-10-09, branch `oi/wp12-reverify`, base `a2e0e65`)

- C22: `ResonanceZone.restore` restores a tail only when it lies past the current one (INFO
  otherwise); `grow_n_times` / `trim_stable_at` docstrings say growth undoes a trim. Test
  `test_restore_after_stable_growth_never_shortens` (fails on base).
- `henon_cases.build_inversion` (the k=4 saddle, `HENON_K4` / `saddle_guesses(4, 1)["inversion"]`
  in `examples/henon.py`) and `_period3_slug` (`p3[_u<n>][_s<n>][_blast<n>]`; the 16-step build
  used to land in the 13-step folder). `test_henon_slugs.py` checks the folders and that the
  example inversion recipe's `describe()` and holes equal the law case's.
- E12 (`stress/numerics/e12_reverify.py`): every case from the frozen `tests/cases.py` builders
  (`e11_blast_depth.build`), one `session.topology()` per blast; all 31 cells ok. Series at 1e-7:
  non-decreasing seed counts and plateaus on every case but inversion (one drop at blast 6, BCG-8 /
  DEC-25); nested_outer 13 flat to 12 blasts, nested_inner 29 from blast 8 to 12, both reliable.
