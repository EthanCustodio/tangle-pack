# loom — knowledge file

Last updated 2026-10-15 (B9a `session-three-concerns` OUTCOME, branch `5f0db5f` off `517aa01`: the
TangleSession / cartoon map below and the "B9" section at the end). Before that 2026-10-14 (B5
OUTCOME, B3 OUTCOME; `Blast.py` line numbers in the B3 bullets are that branch's). The dated sections
below are each WP's notes, kept as the record.

## Map

- `loom/TangleSession.py` — 843 lines on the B9a branch `refactor/session-three-concerns`
  (`5f0db5f`, not yet merged; 1744 at `517aa01`). After B9a it does ONE job, the facade: `trellis`
  (one snapshot per fixed-point SET keyed by `workbench.generation`; a selection of several is cut
  at their pips on every read), `arrangement`, `topology` (one `loom.Topology` per `(generation,
  selection, pips)`), the stage readers `bridge_classes` .. `symbolic_dynamics` and the three
  stage describes (`describe_bridge_classes`, `describe_iterated_partition`,
  `describe_symbolic_dynamics`), the pip memory (`classify_strong_pips` -- still dict-or-list --,
  `_restore_strong_pip`, `set_strong_pip`, `forget_strong_pip`, `strong_pip`,
  `chosen_strong_pip`, `strong_pip_candidates`), `resonance_zone`, `add_resonance_zones`,
  `classify_bridge` (innermost containing zone), `blast_zone` -> `Blast.blast_zone`,
  `_resolve_fixed_points` (also read by `Blast.py:309`), `__getattr__` delegation and the three
  setters. GONE in B9a (base bodies at `git show 517aa01:src/tanglepack/loom/TangleSession.py`):
  all 13 `plot_*` delegates, `cartoon_enclosures`, `cartoon_zones`, `_fanout_call`, `_fanout_plot`,
  `compute_pseudoneighbors`, `punch_holes`, `partition_stable_manifold`, `partition_element_for`,
  `describe_holes`, `describe_stable_partitions`, `_describe`. Callers: `session.trellis(fp).plot_*`,
  `plotting.plot_*(stage object)`, `session.topology(fp)`. A removed name now raises
  `AttributeError` through `__getattr__`.
- `loom/cartoon.py` (209 lines, NEW in B9a) — `cartoon_zones(zones, workbench)` -> `(interior_side,
  parents, arc_order)` (verbatim move; reads the STORED zones, never `Trellis.pips`, author Q2 (a))
  and `plot_resonance_zones(zones, ax=, colors=, alpha=, **fill)` (dropped the unreachable
  `boundary_polygon` fallback: a hand-built zone with `boundary_vertices=None` now raises
  AttributeError). In loom because both take `ResonanceZone` / `ring_order`. Both exported from
  `tanglepack.loom`. Its INFO messages log as `tanglepack.loom.cartoon`. A circle cartoon is three
  explicit calls: `cartoon_zones` -> `plotting.cartoon_enclosures(dual.trellis,
  dual.minimal.kept_bridge_ids, parents, session.fixed_points)` only when `parents` ->
  `plotting.plot_dual_graph_cartoon(..., interior_side=, parents=, arc_order=, enclosures=)`
  (`scripts/symbolic_figures.py::draw_cartoon`, `tests/plotting/test_plot_topology.py::_circle_cartoon`).
- `loom/Topology.py` (207 lines; 212 on the B5 branch) — one recompute (D1): stage 3 per fixed point, stages 4-6 on the
  selection's trellis, lazy `homotopy`, `minimal`, `iterated`, `dual_graph`,
  `dynamics`; `identity` / `letter_for` (C10). B5 branch (`3fb310a`): `_alphabet` :144 (cached,
  zoned, unlettered), `classes` :152 a plain property that runs `self.dynamics` first (the
  verdicts and letters come from `symbolic_dynamics`), so `session.bridge_classes()` costs the
  whole dynamics; no lettering loop here any more.
- `loom/Blast.py` (522 lines on the B3 branch) — `blast_zone(session, zone, n, fixed_point=,
  min_separation=, diagnose=, max_points=, predict=)` :215, `BlastResult` :156 (`stop_reason`,
  `points_created`, `fixed_points`, `too_close`, `already_known`, `pip_depth` (B3) — all
  documented since B3), `_pip_orbits(session, fixed_points)` :109 (`[q0, ..., f^m(q0)]` per fixed
  point with a chosen pip: `registry.find(chosen)` then `iterate_table[x, 1]` link by link),
  `in_zone` (innermost zone, C4), proximity guard (`_interior_points` :85,
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
  the blasted one (C4, rebuilt from `workbench.uniiterated_bridges` on EVERY call). Since B3 an
  iteration first maps EVERY parent (`workbench.iterate_bridge`, which links `f^m -> f^{m+1}`),
  then reads the pip orbits once (`_pip_orbits`, `Blast.py:420-426`), then sorts `produced`
  (`:427-471`): already `seen` -> `already_known`; out of zone / wrong fixed point / iterated ->
  discarded (BEFORE the exemption); an end on a selected pip's orbit -> kept regardless of
  separation, one INFO when `min_separation` is set (`:440-447`, N-a); else within
  `min_separation` of existing curve or a kept sibling -> `too_close`. The verdict feeds only
  `frontier = kept` and the counters, so at `num_iterations=1` (every recipe, every series,
  E12/E13) `min_separation` changes NOTHING in the workbench: a dropped child stays un-iterated
  and in zone and the next call re-feeds it (verified B3). `result.pip_depth` is read after the
  final `infer_iterates` (`:515-517`). Any
  iterate failure propagates (NUM-10). The image curves stay computed
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
- Open for the D8 API plan: `blast_zone`'s filter could default to the zone's own fixed point.
  (`partition_element_for` was deleted in B9a.)

- B10 `dead-api` (2026-10-15) DELETED `trim_stable_at_intersection` (and its two exports; commit
  `217920a`; `TangleWorkbench.trim_stable_at(key, cdist)` is the one stable trim),
  `TangleSession.classify_bridges` (`1348e8c`; callers write
  `{b: session.classify_bridge(b) for b in session.workbench.bridges}`) and
  `BlastResult.terminated_early` (`22feb61`; read `stop_reason == "frontier_empty"`).
  `scripts/henon_blast_period_3.py` is still dead at `:95` (pinned pip 12 is no candidate since B3,
  candidates `[31, 32, 33]`); B10 changed only its `:123` print.

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
- Pip-orbit depth (B3, branch `refactor/depth-of-n`): `BlastResult.pip_depth` and the exemption
  from `min_separation` — `tests/unit/loom/test_blast_pip_orbit.py` (4 tests, `henon_p3_session`,
  3 iterations, sep 0.02), `test_blast_proximity_guard.py` (`pip` / `no_pip`), the series law
  `test_pip_orbit_never_shallower` over `cases.SERIES` (`SeriesStep.pip_depth`, composed
  test-side by `cases._pip_depth`).
- Cost: blasting dominates (p1 8 blasts ~4 s; nested 8+8 ~5 s); every topology stage after it is
  < 0.2 s even at 400 crossings. Inside `iterate_bridge`: `Tangle.add_manifold` segment
  bookkeeping + rtree queries ~50 %, `refine_manifold` ~30 % (cProfile, p1 8 blasts).

## Open-issues round WP3 (2026-10-08, branch `oi/wp03-inversion`)

- `ResonanceZone.ring_cycles(cuts)`: from a cut's stable branch `s(o,b)` to the cut on
  `u(o,b)` (FULL key, not orbit index). A period-k orbit closes one cycle; an inversion point
  one loop per branch pair, both through its anchor; a run reaching a branch with no cut stops
  open. `_build_boundary` emits each cycle `u` arc out, `s` arc back; non-inversion polygons are
  the same cycle (areas equal to 12 digits; the ring starts at an anchor instead of a pip).
  `ring_order` (moved beside it, WP10) feeds `loom.cartoon.cartoon_zones` (empty when a cycle does not close).
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
  `resonance_zone` :1483-1537 (the zone's pip, read from the registry BEFORE the trim; since B8 it
  pops the previous choice, builds `self.trellis(rz.key)` WITHOUT it, and adopts the pip only when
  `registry.find(pip)` is a candidate there (`_is_candidate`), setting `trellis.strong_pip = found`;
  else one WARNING "is not a strong-pip candidate; the zone is stored, the chosen strong pip is
  unchanged", the previous choice put back and `_restore_strong_pip`-ed; with no previous choice the
  session holds none and the next `topology()` takes the default); `add_resonance_zones` goes through
  it with `recompute=False`, so the test reads the trimmed but NOT yet recomputed registry (author's
  default, B8 Q1; every recipe's zone pip passes there, nested and heteroclinic verified 2026-10-13);
  `forget_strong_pip(fp)` :1119 pops `_strong_pips[fp]` AND sets `trellis(fp).strong_pip = None`
  (uncuts), so the next `topology()` classifies and takes the default (B8), `classify_strong_pips` (fan-out: re-applies a remembered pip with `choose_default=False`, so a
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

- Pip facts (B8, measured 2026-10-12 at 8a1cdc3; outcome audited 2026-10-13 at `ae18859`):
  * `forget_strong_pip` exists; deleting `_strong_pips[fp]` alone is still NOT a forget (the cached
    trellis keeps `strong_pip` and `Topology._choose_pip` re-adopts it). The stale-pip WARNING
    (`_restore_strong_pip`) and `ValueError` (`Topology._choose_pip`) both keep "no longer a
    strong-pip candidate" and end "choose another with set_strong_pip() or drop it with
    forget_strong_pip()".
  * A zone at a non-candidate usually trims the chosen pip AWAY (inward on its branch), so
    `topology()` still raises (k10 5/6, k28 x2 9/10, nested 37/38, heteroclinic, inversion x1);
    `forget_strong_pip` is the remedy the message names. Off the chosen pip's branch
    (inversion_two_blasts 33 on (0, 1); pip 1 on (0, 0)) B8 keeps pip 1 and `topology()` runs.
  * Each zone call now builds one more trellis, so `workbench.generation` (an OBSERVATION counter)
    reads +1 per zone on the zone-defining law cases (2 -> 3, nested/heteroclinic 2 -> 4); no state differs.
  * OPEN: some candidates make `topology()` ASSERT once set (`check_bridge_rows_consistent`, equal
    signs on a bridge joined across the cut): k10 2, k28 x1/x2 4, inversion 2/14; reachable through
    the public `set_strong_pip`. A test choosing a non-default pip must avoid them.
  * The stale-pip tests: `test_session_bridge_classes.py:119-170` and (B8)
    `test_session_topology.py::test_a_stale_pip_error_names_the_remedy`.

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

## B3 `depth-of-N` analyst (2026-10-14, `HEAD` `9cd29d5`) and outcome (branch `7d3a294`)

- Where `m` (the depth of `T_m(q0)`) moves: `f^{m+1}(q0)` is born and linked while the PARENT
  bridge ending at `f^m` is iterated (`BridgeIterator.py:363-405`: detect, `registry.add` every
  detection `:371`, cut only at `resolved_ids` `:378-379`, `infer_iterates` of the parent's ends
  `:398`), i.e. before and independent of the blast's child filter.
- `m` per single-iteration blast at the recipe cutoff (composed `+1` chain from
  `chosen_strong_pip`): +1 per blast until the first stall, flat after. Stalls and causes:
  k28 b8 (ii) id 245 ratio 2.991; k32 b9 (ii) 523 / 1.477; nested_outer b9 (ii) 447 / 1.775;
  inversion b10 (ii) 238 / 1.571; p3 b15 and heteroclinic b6 (iii) `f^{m+1}` NOT DETECTED: the
  image's last chord (0.0138 long) passes ~1e-6 on the wrong side of the anchor (0,1) of `s(0,0)`
  while the true point is 8e-6 from it (the unresolved 119/197 there is another fold, u=108186 vs
  207917); nested_inner no stall to 12. (ii) = unresolved detection, never linked (B1).
- Mechanism (i) (`min_separation` drops the pip-bearing child) occurs ONLY inside a
  multi-iteration call: one `blast_zone(nested_outer, 11, min_separation=1e-4)` drops (114,115) at
  iteration 6 and ends at m=7 (series: 9). Pre-B1 (`1029fcc`) the k32 / nested_outer stall at blast
  10 was (iv), the ZONE FILTER: (522,523) / (446,447) discarded exterior at blast 9.
- A prototype of "never drop a pip-bearing child" (scratch copy) leaves every series bit-identical
  and breaks `test_blast_proximity_guard.py::test_min_separation_drops_close_bridges` (its only
  too-close child (46,47) ends at `f^4(q0)`).
- Outcome (option (a)): the analyst's prototype landed as N-a (map above). Series and law cases
  bit-identical on base and branch; `pip_depth` per step k28 `[1..8, 8]`, k32 `[1..9, 9, 9]`,
  nested_outer `[1..9, 9, 9, 9]`, p3 `[14..23]`, inversion `[5..9]`. One
  `blast_zone(nested_outer, 11, min_separation=1e-4)`: m 7 -> 9, and it now fires 361
  same-stability WARNINGs + 1 unresolved-detection WARNING (id 349, ratio 1.775) base does not:
  the exemption keeps near-coincident curve (verified by the auditor). The stalls (ii)/(iii)/(iv)
  are UNADDRESSED: options refine (B1), place by `image_cdist` (PK-6-like), finer image near the
  anchor for (iii). The smells the analyst listed (undocumented `too_close`/`already_known`,
  "sibling" wording) are fixed.
- Smell: `IterateTable.forward_depth` / `forward_chain` are misleading (the table holds `+1` links
  only, depth always 1); untouched, `_pip_orbits` and `Trellis._pip_orbit_cuts` compose instead.
- Smell: `scripts/henon_blast_period_3.py:95` pins `set_strong_pip(12)`, not a candidate (31-33):
  the script raises on `main` before any blast.

## B9 `session-three-concerns` analyst (2026-10-15, `HEAD` `517aa01`)

- Removing concerns 2-4 leaves `TangleSession.py` at ~898 lines (~856 with the describe fan-outs),
  before ~8 import lines; keeping the seven symbolic plot delegates with docstrings would add ~230.
- `cartoon_zones` reads ONLY the stored `resonance_zones` (area, `boundary_vertices`,
  `boundary_intersections`, `contains_point`, `fixed_point.coordinates`), `workbench.manifolds`,
  the registry's `cdist_tol`, loom `ring_order`, topology `ZoneArcOrder`: no pip memory.
  `cartoon_enclosures` reads `workbench.fixed_points` (children may lie outside the trellis's
  selection: `test_plot_topology.py:418-419` needs them), `workbench.manifolds` (= `trellis.manifolds`,
  same live dict), the trellis's registry and `bridge_between`, and `parents`: no pip memory.
- Zones rebuilt from `topology().trellis.pips` (`registry.iterate_orbit(pip, k_value)` + the zone
  boundary builder, live) equal the stored zones BIT FOR BIT (same pip, `np.array_equal` polygon, same
  `cartoon_zones` output) on every law case / figure case / series step that stores a zone (41 of 43
  checked). With no stored zone they ADD geometry: k10's circle cartoon is unchanged, inversion_period2's
  changes (ring order, anchor angles). Probe: scratchpad `b9/probe/zones.py`, `zones_fig.py`.
- `plot_strong_pip` / `plot_strong_pip_candidates` fan-outs classify via the session when a trellis has
  no pip; `Trellis.strong_pip_candidates` classifies WITHOUT choosing, so `Trellis.plot_*` never
  touches the pip memory.
- Every script calling a leaving method runs at `HEAD` (1-14 s each, exported copy, Agg) except the
  known-dead `henon_blast_period_3.py` (`:95`). Nested scripts that call the session stage fan-outs
  (`henon_dual_graph.py:134-136`) need `session.topology(fp)` per fixed point to keep the
  per-fixed-point trellises filled; `topology()` of the selection leaves them empty.

### B9a outcome (2026-10-15, branch `refactor/session-three-concerns` @ `5f0db5f`)

- Done as planned (C1 `d38dd3f` cartoon, C2 `86b453d` plot delegates, C3 `7944f2e` fan-outs, TE
  `5f0db5f`); suite 1208/47/5/1 xfail (C3); 63-item probe, 32 canvases, 42 itinerary PNGs, 11
  scripts byte-identical base vs branch. Net -749 lines.
- Smells left: `ResonanceZone.py:295` (`area`) and `:353` cite `TangleSession` -- B9b's analyst
  found both still accurate (`TangleSession.classify_bridge` orders zones by area;
  `add_resonance_zones` exists), left alone (B9b Q5); the `partition_element_for` cites at
  `Trellis.py:1413`, `TopologyResults.py:427` went with B9b; `tanglepack.loom.plot_resonance_zones`
  has no test (B12, B9b Q4; draft at scratchpad `b9b/smoke_draft.py` on `k28_partitioned`);
  `henon_k10_dual_graph.py` lost the `<FixedPoint ...>:` INFO header the `_describe` fan-out printed.
- Deleted tests (recover from `517aa01`): `test_plot_smoke.py::test_session_plot_fanouts_draw_every_fixed_point`,
  `::test_session_plot_pseudoneighbors_computes_on_demand`, `::test_session_plot_delegates_return_their_types`
  (C2); `test_session_bridge_classes.py::test_a_pip_trimmed_away_is_not_replaced_by_plotting` (C2);
  `test_session_fanouts.py::test_describe_fanouts_mention_every_fixed_point` and
  `test_partition_elements.py::test_session_closes_the_cross_trellis_gap` (C3).
- B9b (2026-10-15, branch `refactor/one-element-at-end` @ `94fca6b`): the fold is done -- topology
  only, no loom file touched; `PartitionFamily.owner_at_end` is the one bridge-end reader (see
  topology.md "B9b outcome"). Loom never read any of the three deleted readers.

## B12 analyst (2026-10-15, `HEAD` `94fca6b`, read-only)

- Q10 recipe constants, measured (scratchpad `b12/q10_probe.py`, `tests/cases._saddle_zone` with the
  zone at `q0` instead of `f(q0)`): on k28 and k32 the zone at `f(q0)` after n blasts has exactly
  the topology (holes, seed / closed classes, active count, reliability) of the zone at `q0` after
  n-1 blasts (n = 1, 2; at n = 0 the `f(q0)` zone has one unresolved class `b`, the `q0` zone none).
  No reason for `f(q0)` is recorded anywhere (decisions, knowledge, open_issues, scripts: only
  "one step in along the stable branch", `scripts/henon_k10_dual_graph.py:51-56`). The k28 / k32
  recipes live in `tests/cases.py:399-418`, not `henon_cases.py`; `henon_cases.py` already states the
  heteroclinic 16 / 10 reason (`:50-62`) and `outer_blasts=2` (`:41-45`); its inversion `f(q0)` has none.
- `tanglepack.loom.plot_resonance_zones` (`loom/cartoon.py:155`) still has no committed test; the
  draft passes on `k28_partitioned` (one zone; `k10_partitioned` has none).
- Committed symbolic-itinerary figures vs a HEAD regeneration (git-archive export, env matplotlib
  3.10.0): content-stale = `heteroclinic_p1_p3` (2.5 M px), `nested_p1_p3_inner6` (1.8 M px),
  `inversion_blast2` (96 k px), all in the itinerary / transition panels (B5's five inert verdicts);
  `p3`, `nested_p1_p3`, `inversion` and the other 7 folders are PIXEL-identical and differ only in
  the PNG `Software` chunk (committed with system `python3`'s matplotlib 3.10.9); k10 / k28 are
  byte-identical.
