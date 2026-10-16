# numerics — knowledge file

Created 2026-10-08 (open-issues round WP2, branch `oi/wp02-exact-iterates`, base `f39f617`); last
updated 2026-10-13 (B8 `unresolved-sets` OUTCOME audit, branch `refactor/unresolved-sets` @
`ae18859`, NOT merged at writing: the iterate-table map entries; the older sections' line numbers are at
their own branches).
Only the crossing / iterate-table path is described; growth, refinement and budgets are in
`CLAUDE.md`.

## Map (crossings and the iterate table)

- `numerics/Tangle.py` — one rtree of segments PER STABILITY (WP10). `resolve_crossings` :220
  resolves each u x s segment pair; `_resolve_crossing_pair` builds each `Intersection` with both
  keys (`_key_of` asserts a key on EVERY curve since WP2), `unstable_segment` and
  `stable_segment` (the two nodes either side of the crossing on each curve), `crossing_sign`,
  `resolution_ratio` and (B1) `chord_deviation` = unfloored `max(h_u, h_s)` (`_resolution_ratio`
  :343 returns both); registration is in unstable-cdist order per branch (reproducible ids, WP8);
  `populate_intersections_for_manifold` :548 (blast path) returns EVERY detection sorted by
  `(unstable_cdist, stable_cdist)` (B1: the noise vote is deleted; the sort keeps INT-7).
- `numerics/Intersection.py` — the crossing record; `stable_segment` added in WP2;
  `chord_deviation` (B1), None exactly when synthetic.
- `numerics/IntersectionRegistry.py` — the store. `DEDUP_RTOL = 1e-9` :35 and `_find_collision`
  :762 (keys equal AND both cdists equal to 1e-9 relative; an anchor at exactly `(0, 0)` matches
  only another exact `(0, 0)`); `register_iterate` :546 asserts the partial bijection;
  `reindex_from` :261 re-gives re-detected crossings their old ids (rematch by `_find_collision`).
  `cdist_tol` (absolute 1e-6) is only the slack of degenerate bridge cuts and polyline clipping
  (C15, WP11). B1: N lives here -- `resolved_ids(threshold)` :441 (synthetic, or ratio >= threshold
  [clause 1, off at None] and in no run [clause 2]), `unresolved_crossings` its complement in
  unstable order, `run_separation` :421 / `run_partner`, `_runs` :489 (per `(u_key, s_key)` group,
  consecutive in unstable order, `_polyline_gap` :50 on the same / node-sharing ORIGINAL segment,
  run iff gap <= max of the PAIR's two deviations), cached on `_orders_version` (not on iterate links).
- `numerics/IterateInference.py` — the EXACT rule (module docstring is the statement):
  `infer_iterates(workbench, sources=None)` :82 (a module function since WP11; reads N via
  `resolved_ids` since B1: pending = sources in N with no `+1` link, candidates in N only) → anchors by
  `_register_anchor_iterate` :142 (advanced keys at `(0, 0)`, `registry.find`), everything else by
  `_bracketed_image` (multi-match = INFO + recorded in `unresolved_links` since B8; carry `U(x)` by the unstable curve's `stretch_param`, `S(x)` by
  `1/stretch_s` of the stable manifold; candidates = a bisected window of the advanced stable
  branch widened by its widest stable bracket `_widest_stable_brackets` :185; unique match on both
  brackets with `BRACKET_PAD = 1e-12`), then `_map_gated` :244 (ONE `map_batch` call over
  `[x, p0, p1]` per source; keep iff `|f(x) - y| <= |f(p1) - f(p0)| + |q1 - q0|`).
- Undecidable links are DATA (B8, branch `refactor/unresolved-sets` @ `ae18859`):
  `IntersectionRegistry.unresolved_links: dict[int, tuple[int, ...]]` (source -> candidate images,
  none linked; reset in `__init__` and `reindex_from`). `_bracketed_image` records a multi-match
  (`unresolved_links[src] = tuple(images)`, INFO kept). `infer_iterates` (IterateInference.py
  ~:118-156) snapshots the set, pops its own sources' entries (each pass re-decides them), groups the
  gated matches BY IMAGE, and for each image `j` forms `contested` = this pass's claimants + the held
  `table[j, -1]` + every earlier entry equal to `(j,)`; one member -> linked, else ONE WARNING
  "image %d is claimed by crossings %s; none linked", `remove_iterate(held, 1)` and every claimant
  recorded as `(j,)`. Order-free (grouping happens before any link); the generation is bumped
  (`registry._generation += 1`, a private write from another module) only when the set changed, so a
  no-op re-pass is silent. The early `return` on "no detected source" is gone (the compare must run).
  Measured: fires NOWHERE in the corpus after B1 (every builder, all five SERIES, k28 x8); the set is
  empty on every law case and series step; `near_tangency` records 25 multi-match entries (e.g.
  25, 27 -> (61, 59)), which nothing reads.
- `IterateTable.forward_depth(x)` / `forward_chain(x)` (`IterateTable.py:81-96`) read DIRECT
  entries `[x, n]` only; only `n = 1` links are ever registered (`IterateInference.py:147,202`), so
  they return 1 and `[x, f(x)]` for any linked crossing. The real depth is the composed `+1` chain
  (`Trellis._pip_orbit_cuts`, `Trellis.py:322-324`). Misleading names (B3 analyst, 2026-10-14).
- `IterateTable.remove(source, n)` pops the forward link and its backward twin (KeyError if absent);
  `IntersectionRegistry.remove_iterate(source, n)` calls it, bumps the generation and drops the graph
  edge `iter_{src}_{n}` (added by `register_iterate`).
- `reindex_from` runs on a FRESH registry (`compute_intersections` replaces it), so any per-registry
  set (B8's `unresolved_links`) starts empty after a recompute; `reset=False` keeps it.
- Callers: `TangleWorkbench.compute_intersections` :1264, `BridgeIterator.iterate_bridge` :264
  (sources = the parent's endpoints + the crossings just born, searched against the whole
  registry), `iterate_all_bridges` :419 and `loom/Blast.blast_zone` (:463, registry-wide, after
  the loop, whatever the stop reason; not after a raised `PointBudgetExceeded`).
- Growth and trims (C22, WP12): `grow_n_times` re-finds every branch's physical end in its
  `finally` (`_find_tail`), so stable growth undoes a stable trim (`trim_stable_at` moves only the
  tail pointer); a resonance zone's `restore` never shortens a branch growth has run past.

## Why it works (and the probes that showed it)

- The machine builds every node by mapping another: `cdist(f(p)) = stretch_u * cdist(p)` on an
  unstable branch, `cdist(f^-1(q)) = stretch_s * cdist(q)` on a stable one; refined nodes take
  the MEAN of their neighbours
  (`ManifoldMachine.py` :360, :841), and the mean commutes with the scaling: the brackets are
  carried exactly. `stretch_param` is uniform over a chain (probed: p3, k28, inversion, nested).
- `per_step_beta` must NOT replace `stretch_param`: they differ 3e-7 (k28) .. 6e-6 (inversion)
  relative, while deep k28 brackets are ~3e-8 wide (links at 9 k28 blasts drop 248 -> 81 with it).
- The map is a gate, never the matcher: at 9 k28 blasts crossings at u = 1.1e3 and u = 1.6e9 sit
  3e-6 apart in the plane (lambda lemma).
- BCG-3 root cause: the old endpoint path's ABSOLUTE gate `1e-3 * max(1, |f|)` rejected six true
  links on p3 at 10 blasts (`|f(x) - y|` 1.2e-3 .. 3.5e-3 on ~1e-2 segments) and no pass ran
  after a blast. Three more near-anchor links (stable cdist 5-21 % off the beta prediction because
  the stable segment there spans `[7e-5, 2.6e-4]`) failed the old 5 % window on both paths.
  Measured in segment units every true link on every builder is <= 0.29 of the gate.
- Refined points cache only a preiterate, so point-level `next_iterate` links cannot register
  crossings; the cdist parameter can.

## Invariants

- The iterate table is a partial bijection (`register_iterate` asserts; law
  `iterate_table_partial_bijection`). Every detected `+1` link passes the map gate (law
  `iterate_links_pass_the_map_gate`) and no unlinked detected crossing has a bracket+gate image
  (law `iterate_table_complete`; regression `test_iterates_after_blasts` on p3 at ten blasts),
  except the sources in `registry.unresolved_links` (B8): each entry has >= 1 candidate and no `+1`
  link, and a single-candidate entry `(j,)` never coexists with a `(j, -1)` link. The table never
  holds a link chosen by the order of the passes.
- Every detected crossing carries both branch keys (PK-7 closed).

## Smells / deferred

- Absolute `cdist_tol` consumers (~86: anchor detection, StrongPip/Pseudoneighbor slack,
  PartitionFamily fallback, DualWalk `_end_tolerance`, `Tangle.create_bridges`,
  `BridgeIterator._image_bridge_ids`) — SIMP-13, WP11. WP2 only removed the dedup use.
- `reindex_from`'s remap return value is unused (SIMP-12, WP11). Exact iterate registration does
  NOT make `reindex_from` unnecessary: `compute_intersections` replaces the registry and only the
  rematch gives a re-detected crossing its old id. Every rematch on every builder is a
  bit-identical re-detection, so the relative dedup keeps it correct (INT-7, answered).
- (Fixed in WP8, see below.) Ids were not reproducible across builds: `resolve_crossings` registered in the iteration order of
  a `set` of segment-id pairs (`Tangle.py` :265). Sorting by `(keys, unstable_cdist)` before
  registering is ~3 lines but renumbers every crossing and the id-based fixture facts in
  `CLAUDE.md` — deferred to WP8 (identity across recomputes).
- `_widest_stable_brackets` is an O(N) scan per inference pass (one per `iterate_bridge`); the
  suite time did not move (106 s vs 108 s base), but a deep blast series could feel it.
- `IterateInference` is a stateless holder; module functions would read better (API pass, D8).

- B10 `dead-api` (branch `refactor/dead-api`, 8 commits `1867f41..517aa01` off `65edec4`, 2026-10-15)
  DELETED: the registry's networkx graph (`IntersectionRegistry.graph`, `_rebuild_adjacency_edges`,
  `_consecutive_on_branch`, the `_graph*` fields, node adds in `add`/`reindex_from`, edges in
  `register_iterate`/`remove_iterate`), `numerics/graphviz.py` (358 lines) with
  `TangleWorkbench.build_intersection_graph` / `visualize_intersection_graph` (commit `7eaf64d`);
  the registry queries `on_cdist_range`, `from_branch`, `nearest_by_unstable_cdist`,
  `as_forward_array` / `as_backward_array` and `IterateTable.as_forward_array` / `as_backward_array`
  (`3a9db04`); `ManifoldMachine.grow_x_times` (`1867f41`, callers loop `new_grow_manifold`);
  `_get_refined_point` / `_get_iterate` (`88d5a75`; the kernel test now calls `_refine_layer` at
  `area_cutoff = 0`). To recover any of them: `git show 65edec4:<path>`.
  - The registry now has NO graph: iterate relations live only in `iterate_table[id, n]`
    (`IterateTable`), per-branch order only in the sorted cdist arrays and `Trellis`. The query
    surface is `filter`, `from_fixed_point`, `iterate_orbit` (all live). A `from_branch` query is
    `filter(lambda ix: ix.is_anchor or (ix.manifold_b_key is not None and ix.manifold_b_key[3] == b))`
    (the two demo scripts use exactly that).
  - `import tanglepack` no longer loads `networkx` (topology imports it lazily for
    `DualGraph.graph()` and the transition graph; it stays in `pyproject.toml`).
  - `Tangle.release_manifold` stays (live): `compute_intersections`, `_register_bridge` and
    `clear_bridges` call it; a no-op on a full build, real work on the direct `clear_bridges` path,
    pinned by `test_bridge_identity.py:257`.
  - Stale outside src after B10 (B12's): `CLAUDE.md` registry row ("per-branch `graph(bridges=)`")
    and `BridgeIterator ... graphviz` row; `docs/source/tanglepack.numerics.rst:183-189` still
    automodules `tanglepack.numerics.graphviz` (a Sphinx build will warn);
    `scripts/manifold_cutting_test.ipynb` and `notebooks/archive/period_three.ipynb` still call
    `grow_x_times` (plan Deferred list).

## Open-issues round WP3 (2026-10-08, branch `oi/wp03-inversion`)

- `TangleWorkbench._register_anchors` registers ONE `Intersection.synthetic` per
  `(fixed point, orbit index)`, keyed `(u(o,0), s(o,0))`, sign from the two branch-0 root
  segments, NO `unstable_manifold` / segments (`Intersection.is_anchor`: synthetic at exactly
  (0, 0)). A point whose branch-0 manifolds are missing or have < 2 nodes gets none.
- `IntersectionRegistry._find_collision`: at EXACTLY (0, 0) two keys match when they name the
  same fixed point and orbit index, so the root segments of the other three branch pairs (which
  ARE detected at exactly (0, 0)) and `find` probes with advanced keys resolve to the one anchor.
  `IterateInference._register_anchor_iterate` needed no change. `_consecutive_on_branch` (graph
  view) puts each anchor first on every branch group of its point.
- `Tangle.create_bridges`: an anchor (no curve) cuts every unstable manifold ROOTED at its
  periodic point (`manifold.root is fp.branch_points[orbit]`) at `(nodes[0], nodes[1])`.
  `_both_anchors` is gone. In `iterate_bridge` the image's own (0, 0) detection carries the
  anchor id (registry collision), so an anchor-rooted image is cut like any other (INV-10).
- One chain rule: `ManifoldInitializer.construct_kevin_way` seeds a chain of `k_value` pieces per
  requested eigendirection, walked by `advance_key`, skipping one an earlier chain covers;
  `FixedPoint.get_iterable_array` / `get_branch_array` deleted. `orient_manifolds` no longer
  returns early under inversion. `ManifoldMachine.new_grow_manifold` walks an inversion chain from
  branch 0 without warning; the workbench's `growth_branch` special cases are gone
  (`_growth_key` still maps an inversion point's stats to one key).

## Open-issues round WP8 (2026-10-08)

- Registration order (INT-7): `compute_intersections` registers `sorted(resolve_crossings(), key=
  (unstable fp position, orbit, branch, unstable_cdist, stable fp position, orbit, branch,
  stable_cdist))`, positions from `self.fixed_points`. Before, the order was the iteration order of
  a set of frozensets of segment ids, and `Tangle._ids` is a process-wide counter: reproducible for
  the first build of a process, not for a second build in the same process. Anchors are still
  registered first; blast registration was ordered by the vote's sort until B1, now by
  `populate_intersections_for_manifold`'s own sort (`e201ddc`).

## Open-issues round WP9 (2026-10-08, branch `oi/wp09-heteroclinic`)

- NEST-5: `clear_bridges(fixed_point=None)` releases and drops only that fixed point's bridges
  and partial pieces (in place), rebuilds `_bridges_at` from the survivors, bumps the generation.
  `rebuild_bridges(fp)`: ids current (`_bridge_id_epoch == _registry_id_epoch`) -> clear and
  recut `fp` only, carrying `iterated` by id; ids stale and bridges held -> DEBUG line and every
  fixed point is rebuilt (the survivors' ids would name other crossings). `create_bridges`
  stamps the epoch only from an empty set, which the stale path guarantees. Recutting a
  fixed point twice reproduces its ids and adds no partial piece.

## Open-issues round WP10 slice 1 (2026-10-08, branch `oi/wp10-numerics-leftovers`, base `3f27456`)

- `Tangle` keeps one rtree per stability (`_rtrees`, `_segment_index`): `_insert_segment` inserts
  into its own stability's index and queries only the other; `add_manifolds` stream-loads both and
  queries each unstable segment against the stable index (each u x s pair once, the lower-id dedup
  is gone). `_resolve_crossing_pair` asserts u x s; `_discard_same_stability`, its counter and the
  INFO summary are deleted. `k32_series` 18.0 -> 14.5 s with identical outputs.
- `create_bridges` brackets every crossing on the CURRENT polyline first: a run of crossings on one
  ORIGINAL segment gets a separator spliced after the previous crossing's lower bracket
  (`_insert_crossing_separator(after, p0, p1, ...)`, cdist still interpolated on `[p0, p1]`), then
  bridge `(k, k+1)` = `(below(k), above(k+1))`. Consecutive bridges overlap on the two bracket
  nodes of their shared crossing, as for any crossing. No law case has two crossings on one
  segment (bridge node totals identical).
- `ManifoldMachine._image_chains(manifold, images)`: the existing images as contiguous chains,
  split where the walk from the first image breaks (to the curve's end or past the last image's
  cdist); unreached images in between are merge-dropped duplicates, skipped as the one-chain merge
  did. Both the partly mapped and the all-images-exist paths merge chain by chain. The rejected
  alternative (one chain per run of iterated nodes) changed k28 at 6 blasts.
- `TangleWorkbench._register_bridge` asserts the bridge invariants (NUM-14). `DEDUP_RTOL` slack.
- `diagnostics._hits` is the one array core; `same_stability_crossings` builds objects from it,
  `count_same_stability_crossings(workbench, bridges)` counts (only pairs with a new segment:
  `_candidate_pairs(a, b, new)` filters `new[i] | new[j]` inside each block). `_curves` returns
  `(name, nodes, bridge)`; a partial image bridge is named `partial image bridge 0x...` (unique).
  `report_resolution(new_bridges=)` reports the count (`_report_count`) and stores nothing.
- Growth drivers: `_stretch_per_round`, `_default_cap` (`max(10, ceil(log_stretch / ln s) + 1)`),
  `grow_until_arclength` on polyline arc length; module Dev Note rewritten.
- Probes: `/home/dezhu/.claude/jobs/c5c98371/tmp/wp10/` (planner) and `.../tmp/wp10c/` (coder:
  `sig.py`, `series_time.py`, `alllaws.py`, `prof_topo.py` with worktree-relative paths).

## Open-issues round WP11 slice 1 (2026-10-09, branch `oi/wp11-simplify`, base `4f7e815`)

- `IterateInference.py` is a module function `infer_iterates(workbench, sources=None)` plus private
  helpers (`_register_anchor_iterate(registry, ...)`, `_widest_stable_brackets`, `_bracketed_image`,
  `_map_gated`); the stateless class, the workbench's `_iterate_inference` attribute and the
  `infer_iterate_table` forwarder are gone (scripts / stress call `infer_iterates`; the stress stage
  label `"infer_iterate_table"` is kept so recorded results stay comparable). Logger name unchanged.
  `_widest_stable_brackets` scanning the registry per pass was measured negligible (0.10 s of 13 s of
  blasting on `k28_deep_blast_sliver`).
- `TangleWorkbench.trim_stable_at(key, cdist)` is the one stable trim (first node past the root at
  or beyond `cdist`; no node -> WARNING, untrimmed); `trim_stable_manifolds` and the loom zone trims
  call it.
- `FixedPoint`: `(2,)` orbit points and eigenvectors, float eigenvalues (solver and placeholders);
  `partial_jacobians` placeholders None; `preserves_orientation` = sign of `det(partial_jacobians[0])`
  (any period, None unsolved, C16).
- `ManifoldMachine._manifold_from_points(points, stability, stretch, fp)` replaces the initializer's
  `construct_manifold_from_point_list` (no lazy initializer in `iterate_manifold`);
  `get_initial_fundamental_segment` is one code path.
- Dead code deleted: `IterateTable.backward_chain` / `all_registered_ids`, the registry's
  `unstable_order_array` / `stable_order_array`, the scalar `_curvature_area` and its three helpers
  (`_parabolic_fit` stayed then; deleted in WP10b, it had no caller), three `Intersection.synthetic` keywords;
  `reindex_from` returns None. (`from_branch`, which held the anchor for every b (D4), was deleted in B10.)
- `PointBudgetExceeded.blast_context` declared; `MaxIterationsReached`, `TURNAROUND_LOG_STRETCH`
  exported.
- `cdist_tol` (registry) is documented as geometry-only: degenerate bridge cuts and polyline clipping.

## Open-issues round WP10b (2026-10-09, branch `oi/wp10b-leftovers`, base `f9d1be0`)

- NUM-12 / C20: a blast's image curves are computed manifold. `iterate_bridge` wraps the machine's
  image as a `Bridge` with the advanced key BEFORE indexing it and appends it to
  `TangleWorkbench._image_curves`; a partial piece is that curve object itself (so its `iterated`
  flag lives on the curve and survives a rebuild). `compute_intersections` re-indexes every image
  curve of a selected fixed point query-only right after the bulk load and drops a curve that adds
  no segment of its own (`Tangle.add_manifold` now returns the count of new segments: 0 means growth
  merged the curve into a branch, edges are deduped by point identity; the dropped curve's claims
  are released). That reading holds only on a FRESH index: under `reset=False` the bridges cut out
  of the curve still own its edges, so the count is 0 for a live curve and every curve is kept
  (review fix; before it, `reset=False` after a blast lost k28 two blasts 10 -> 8 crossings on
  the next recompute). `create_bridges` cuts image
  curves like any curve and holds a curve of an indexed fixed point (one with an anchor in the
  registry) crossed fewer than twice as one partial piece. `clear_bridges` releases only id-carrying
  bridges: a partial piece is computed manifold and stays indexed. `_recut_bridges_after_growth`
  (flag clearing + WARNING) is gone; `grow_until` calls `rebuild_bridges`.
- Before: ANY recompute after a blast dropped every crossing on an image curve (k28 two blasts
  10 -> 4 crossings, classes 4 -> 3; p3_6_blasts 23 -> 9 classes); now a recompute + rebuild leaves
  bridges, flags, holes, classes and the exact `describe()` unchanged on every law case,
  `k28_deep_blast_sliver` and `p3_6_blasts`, and blast / recompute / blast equals blasting straight.
  Cost: recompute after 9 k28 blasts 1.4 -> 3.8 s.
- Deferred (found here): growing the STABLE manifold after a zone trim un-trims the zone
  (`grow_n_times` ends with `_find_tail()` on every manifold of the grown stability); topology
  unaffected (D10), the zone object goes stale. WP12 / D8 API plan.
- NUM-19: `_register_bridge` asserts a `BridgeId` runs in unstable-cdist order (an anchor end
  first). `_validate_saddle` near-parabolic: not a bug (re-measured, `num19_solver.py`).
- NUM-20: `IterateInference` was right (each side advances on its own saddle); the laws
  `area_along_iterate_links` / `iterate_link_scales_by_beta` now read each side's own
  `per_step_beta` (product factor `beta_u(A) * beta_s(B)`, exactly 1 homoclinic).
- `ManifoldMachine._parabolic_fit` deleted with its pinning test (no caller).

## Resolved crossings — B1 (2026-10-12, branch `refactor/resolved-crossings` @ `76d0710`, base `1029fcc`)

- **Rule** (A.2, Q13 (a)): every detection is registered; N = `registry.resolved_ids(
  workbench.resolution_ratio_threshold)`. Seven readers: `TangleWorkbench._frozen_ids` :1102,
  `create_bridges` :1803, `trim_stable_manifolds` :2178, `BridgeIterator.iterate_bridge` :378 (cut
  and inference sources), `infer_iterates` :109 / `_bracketed_image` :226, `Trellis.from_workbench`
  :422 (records `Trellis.resolved_ids`), `StrongPip.is_strong_pip` :169. An unresolved detection
  lies INSIDE a bridge, is never linked, sets no trim, disqualifies no pip.
- **Report:** `_report_unresolved` -- ONE WARNING per call when a not-yet-reported id is
  unresolved, naming each (id, coords, ratio, clause, partner + separation) and ending "refine: lower
  `area_cutoff` or grow the stretch; NUM-3/4". `threshold=None` keeps clause 2, so it can fire.
- **Threshold docstring** (`TangleWorkbench.__init__`): derived bound 1, margin 3, part of the
  definition of N, corpus minimum 3.038 (heteroclinic); re-run `compute_intersections` after a change.
- **Measured (outcome audit):** law cases unchanged (min ratios k10 522, k28 x1 5783 / x2 2455, p3
  family 8.70, inversion 124 / 19.5 / 12.2, period2 19.8, heteroclinic 3.038). k28 series: run
  {75,76,77} at blast 6 (u = 8356210.197, sep 8.5e-10 / 2e-8 vs h 4.7e-6, median ratio 30.2) out
  whole; blast 8 adds 148 (2.572) and 245 = `f^9(q0)` (2.991). Recompute exact on all 5 series.
- **The pip-orbit freeze** (not predicted by plan or analyst): an unresolved pip-orbit point is
  never linked, so `Trellis.cut_at_pip_orbit`'s `m` (`Trellis.py:322-323`, `iterate` walk) stops
  before it and no later blast can extend `T(q0)` until a finer rebuild. k28 x8 = x7 (143 crossings,
  4/4/4 vs base 24/15/23); k32 x9-10 = x8 (300; 56/34/45 reliable vs base 58/34/39 unreliable;
  523 = `f^10(q0)`, 1.477); nested_outer x9-11 = x8 (276; 12/9/9 both; 447 = `f^10(q0)`, 1.775).
- Suite cost of the per-step recompute in `cases.blast_series`: +38 s (162 -> 200 s).

## Smells (crossing set)

- `Trellis.py:204`: a hand-built `Trellis` (tests, the `MinimalTrellis` sparse snapshot) defaults
  `resolved_ids` to the whole registry; harmless while no N reader runs on the sparse snapshot; the
  fix is `MinimalTrellis.py:411` passing `parent.resolved_ids`.
- `_runs` neighbours are consecutive in unstable order within a KEY group, but blast image curves
  share their branch's key; a run interleaved with a same-key detection on another curve would be
  missed (believed; the probe's any-curve cross-check found none on the law cases or k28 x8).
- A crossing already cut and linked can be demoted later in the same blast (plan Q1, no repair);
  the recompute law over `SERIES` is the guard.
- `CLAUDE.md` still says resolution flags "never discard" and states the old k=2.8 x8 fixture
  facts; deferred to B12.
