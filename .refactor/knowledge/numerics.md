# numerics — knowledge file

Latest: 2026-10-28, D1 `diagnose-on-new-segments` OUTCOME audit, branch @ `7d9d009`, base `88c62ea`, NOT merged at writing (the section "D1 outcome" at the end describes the branch's `diagnostics.py`; the analyst section before it describes the base).
Created 2026-10-08 (open-issues round WP2, branch `oi/wp02-exact-iterates`, base `f39f617`); last
updated 2026-10-27 (N7b `delete-point` OUTCOME audit, branch @ `647786b`, base `4c9c3cb`, NOT merged at
writing: the section "After N7b `delete-point`" at the end describes the branch; every `Point` /
`BasePoint` / `store.nodes` / `node_at` / `return_nodes` / `_walk_nodes` / `first_node` /
`_manifold_from_points` / `loose_store` mention in the older sections is history). Before that
2026-10-26 (N7a `engine-on-rows` OUTCOME audit, branch @ `88f1a47`, base `d279815`, NOT merged
at writing: the section "After N7a `engine-on-rows`" at the end describes the branch). Before that 2026-10-26 (N7 `delete-point` ANALYST at main `d279815` = N6 merged fast-forward: the section
"Deleting Point (N7 analyst)" at the end describes main and two scratch prototypes). Before that
2026-10-26 (N6 `readers-on-views` OUTCOME audit, branch @ `e25adf1`, base `9dc702b`, NOT merged
at writing: the section "After N6 `readers-on-views`" at the end describes the branch; the N6 analyst
section before it describes the base). Before that 2026-10-26 (N6 `readers-on-views` ANALYST at main `9dc702b`, N5b merged: the section "Readers
on views (N6 analyst)" at the end, measured with a scratch prototype; it describes main). Before that 2026-10-26 (N5b `segments-and-edge-keys-on-rows` OUTCOME audit, branch @ `c944cf6`, base `1c25afb`, NOT
merged at writing: the section "After N5b `segments-and-edge-keys-on-rows`" at the end; the `Tangle.py` map
entry below is corrected to it). Before that 2026-10-26 (N5b ANALYST at main `1c25afb`, N5a merged: the
section "Segments and edge keys on rows (N5b analyst)", measured with a scratch prototype; it describes the base). Before that 2026-10-25 (N5a `brackets-on-rows` OUTCOME audit, branch @ `9b72985`, base `01e991a`, NOT
merged at writing: the section "After N5a `brackets-on-rows`" at the end; the map entries for
`Intersection` / `IterateInference` below are corrected to it). Before that 2026-10-25 (N5
`segments-and-brackets-on-ids` ANALYST at main `01e991a`, N4 merged: the
section "Segments, brackets and edge keys (N5 analyst)", measured). Before that 2026-10-25 (N4 `curves-on-order-arrays` OUTCOME audit, branch @ `ab7bd20`, base `8fc3933`, NOT
merged at writing: the section "After N4 `curves-on-order-arrays`"; the N4 analyst section and the
N3 sections describe the base's linked lists). Before that 2026-10-25 (N4 ANALYST at `8fc3933`, N3
merged: the section "Curves, splices and shared nodes", measured). Before that 2026-10-24 (N3 `point-store` OUTCOME audit, branch `refactor/point-store` @ `5ce5685`, base
`138b34c`, NOT merged at writing: the section "After N3 `point-store`"; the N0 data-model section
describes the base). Before that 2026-10-24 (N3 ANALYST at `138b34c`: the section "Points, phantoms and
the store (N3 analyst)", measured with a scratch prototype). Before that 2026-10-24 (N1 `slim-segment-index` ANALYST at `41056f8`: the section "The segment index's
bookkeeping", measured). Before that 2026-10-23 (N0 `numerics-arrays-skeleton` OUTCOME audit, branch `refactor/numerics-arrays-skeleton`
@ `41056f8`, base `a8bb8a0`, NOT merged at writing: the profiler, probe and campaign plan doc are
committed on the branch; before that the N0 ANALYST wrote the "Data model and where memory / time go"
section, measured). Before that 2026-10-21 (densify-near-anchor ANALYST at
`a8bb8a0`: the "Refinement" section and the reindex / stall smells). Before that 2026-10-13 (B8 `unresolved-sets` OUTCOME audit, branch `refactor/unresolved-sets` @
`ae18859`, NOT merged at writing: the iterate-table map entries; the older sections' line numbers are at
their own branches).
The crossing / iterate-table path, plus (2026-10-21) how refinement creates points and how
accurate they are; growth drivers and budgets are in `CLAUDE.md`.

## Map (crossings and the iterate table)

- `numerics/Tangle.py` — one rtree of segments PER STABILITY (WP10). Since N5b a `_Segment` is
  `(id, manifold, store, row0, row1, in_rtree)` (rows; the root `BranchPoint` in `row0`), every edge
  key is one int of rows (`_edge`, `_edge_keys`), the query-only add is `_query_curve` (one
  `intersection_v` per curve) and the Tangle holds one store (`_store`); see "After N5b" at the end.
  (Line numbers in this entry are older.) `resolve_crossings` :220
  resolves each u x s segment pair; `_resolve_crossing_pair` builds each `Intersection` with both
  keys (`_key_of` asserts a key on EVERY curve since WP2), `unstable_segment` and
  `stable_segment` (the two nodes either side of the crossing on each curve), `crossing_sign`,
  `resolution_ratio` and (B1) `chord_deviation` = unfloored `max(h_u, h_s)` (`_resolution_ratio`
  :343 returns both); registration is in unstable-cdist order per branch (reproducible ids, WP8);
  `populate_intersections_for_manifold` :548 (blast path) returns EVERY detection sorted by
  `(unstable_cdist, stable_cdist)` (B1: the noise vote is deleted; the sort keeps INT-7).
- `numerics/Intersection.py` — the crossing record; `stable_segment` added in WP2;
  `chord_deviation` (B1), None exactly when synthetic. Since N5a the two brackets are STORE ROWS
  (the periodic point's `BranchPoint` in place of a row), with `store`, `segment_ends(stability)`
  (coords + cdists, one view of each column per call) and the `brackets` memo (sorted, unpadded).
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
- (Fixed in N5a.) `_widest_stable_brackets` was an O(N) scan per inference pass re-reading every
  crossing's Points (1.7 s at scale); it is now one `np.fmax.reduce` per stable bucket.
- `IterateInference` is a stateless holder; module functions would read better (API pass, D8).
- (D1 analyst, 2026-10-28, main `88c62ea`) LATENT CRASH: `Tangle._query_curve` (`Tangle.py:610`) calls rtree
  1.4.0's `intersection_v`, whose resize branch (`rtree/index.py:1106` `ids = ids.resize(...)`, which returns
  None) raises `TypeError` as soon as one call's results exceed 2x its query boxes (reproduced: one box over three
  entries). u x s queries stay under it so far: max results/queries 1.31 (p3_series), <= 0.4 elsewhere, 0.17 at
  scale (scratch `d1/iv_margin.py`). A short query-only curve over a dense stable bundle would crash a blast.
- (D1 analyst) `diagnostics._candidate_pairs`' docstring (`:440-443` at the base; `:478-481` on the D1 branch, now
  with both reasons) gives the right reason to avoid `intersection_v` (the same resize bug), but a u x u query would
  hit it on every call (~80 overlaps per box). The `_query_curve` crash above is STILL OPEN after D1 (author's Q7:
  not in the diagnose campaign; a small `query-curve-guard` run is the option).

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
  `rebuild_bridges` (`TangleWorkbench.py:1866`) makes NEW Bridge objects and carries ONLY
  `iterated` (`:1904-1923`); `create_bridges` on top of a held set returns the HELD object for a
  known id (`_register_bridge` :1753), so any other attribute survives `create_bridges` but not
  `rebuild_bridges` (verified 2026-10-21, blast-steps-compose probe_q1). On branch
  `refactor/blast-steps-compose` (`a6d2ffe`) `Bridge.dropped_at` (the loom blast's `min_separation`
  drop, `Bridge.py:122`) is carried beside `iterated` in one `(iterated, dropped_at)` dict under the
  same id-epoch guard (`:1904-1926`), so it is lost exactly when `iterated` is.

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
  `count_same_stability_crossings(workbench, bridges)` counts (only pairs with a new segment; since
  N2 `_candidate_pairs(a, b, new)` BUILDS only those: a partner table `pool`/`start`/`count` per sorted
  position, an old position taking only the new positions of its forward range). `_curves` returns
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

## Refinement: how points are made and how good they are (densify analyst, 2026-10-21, `a8bb8a0`)

- `ManifoldMachine.refine_manifold` (:606) -> `_refine_layer` (:671): breadth-first; a pair is
  skipped when shorter than 1e-8 ABSOLUTE (:715) or when either end has no one-step preiterate
  (:718-726: `prev_iterate` on an unstable curve); flagged iff curvature AREA >= `area_cutoff`
  (:765; `h = 1.5 A / L`, so a fixed area is a length-dependent h); `_charge` (:769) BEFORE the one
  batched `map_fwd_batch(pre_mid)` (:771); the new point is `f(mid(pre0, pre1))` with the cdist
  mean, its preimage cached by `_cache_preiterate` (:814-833) as a PHONEY `Point` (coords only, no
  cdist, on no list, `next_iterate` None) in `new_point.prev_iterate`; spliced by
  `_insert_point_geometrically` (:786). It NEVER inserts into or links the preimage curve, never sets
  `next_iterate` on an existing point: no point becomes iterated or un-iterated. A new point has no
  `next_iterate`, so a later growth or bridge iteration maps it (verified: densified p3 branches, one
  `grow_n_times` maps the 56 extra points, +99 created).
- Stored one-step preiterates equal `map_inv(p)` to the bit (verified on the k=2.5 stall segments).
  The prev-iterate chain is short: a refined point's phoney preimage ends it (chain depth 1-5 seen).
- ACCURACY LAW (area preservation): the new point lies on `f(preimage chord)`, whose lens with the
  true arc has the preimage lens area, so its error ~ `h_pre * l_pre / l` (preimage chord deviation
  times preimage length over image length). Measured predicted = actual (dense reference) to 10%.
  Near a period-1 anchor the preimage of a near-anchor pass is straight (`h_pre` ~1e-16): curvature
  is made in ONE step by the map, the rule is exact (k=2.5 p1 / p3 stall segments: error 2e-17 at
  depth 1, WORSE at depth 2-3: 1e-8 / 1e-6). Elsewhere the curvature is inherited: on k=2.5 p3,
  13 of 42 new points within 1e-3 of the orbit had true error 1.3e-7..1.6e-6 (0.05-0.65 of the
  pair's prior h), so a polyline refined by the rule converges to `f(pre chord)`, not to the curve,
  and its vertex-based h (hence `resolution_ratio`) then OVERSTATES the resolution. Growth's own
  refinement has the same law at its own cutoff.
- `Tangle._insert_crossing_separator` (:972-1010) puts a vertex at the midpoint of two crossings ON
  THE CHORD (off the curve by up to h) with no preiterate: the refiner skips its pairs.

## Data model and where memory / time go (N0 analyst, 2026-10-23, `a8bb8a0`, measured)

Map (the linked-list data model the arrays campaign replaces):
- `BasePoint.py` (x, y, `_coords` (2,) array, cdist, edist, `next_iterate` / `prev_iterate` = the map's
  forward / backward iterate), `Point.py` (+ `forward` / `backward` geometric links,
  `stretch_param`), `BranchPoint.py` (the periodic point: 2 slots per stability,
  `forward_branches` = unstable roots, `backward_branches` = stable roots, scalar `cdist = 0.0` AND
  tuple `cdists = (0, 0)`). No `__slots__`: one `Point` = 400 B (tracemalloc, 1e5 linked; object 48 +
  `__dict__` 144 + `_coords` 128 + x, y `np.float64` 32 each + cdist 24).
- A stable curve links AWAY from its root through `backward` (`BaseManifold.walk_fwd` :485 reads
  `node.backward` when stable); unstable through `forward`.
- `BaseManifold.tail` is an EXTENT, not the physical end: a trim (`trim_stable_at` :2068) or a zone
  moves it inward; `grow_n_times` resets it to the physical end (`_find_tail`, :708-711).
- Refined points hold a coordinate-only PHANTOM preimage (`_cache_preiterate` :814-835). They are
  15-45 % of all `Point` objects (k10 1457 of 3190; 1e-7 scaling case 778 k of 2.18 M).
- Off-curve points also exist as merge-dropped duplicates (kept alive by their preimage's
  `next_iterate`: exactly one per merge tie, 17 on k28, 2 on heteroclinic / k32, 0 on the scale case)
  and as the initializer's k-step "fictitious" chain (`ManifoldInitializer.py:172-191`, adopted onto
  curves later, except its detached LAST point, cdist None, which is garbage: 1 per seeded chain).
  Nothing else is ever freed (N3 analyst: created == alive + those tails on every case and at scale).

Measured (profiler committed on the N0 branch as `stress/profile/profile_numerics.py` +
`profile_figures.py`, baseline under `figures/numerics_profile/a8bb8a0*/`; campaign plan with the
full tables, inventory, five hard places and ten BranchPoint roles in
`docs/numerics_arrays_plan_2026-10-23.md`, sections B, E, F, G):
- The census matches tracemalloc to 1.01-1.02 on the scaling case (8 blasts and the 2.18 M-object
  plateau) and 0.96-0.99 on nested / heteroclinic. Census bytes at rest are ~810-960 B per `Point` object (phantoms counted as points), RSS 1.1-1.6 KB.
  Shares at the 1e-7 k = 2.5 p1 plateau (2.18 M Point objects, 1.40 M segments, 30 blasts):
  `Tangle`'s Python bookkeeping 55 % (`_edge_to_sid` frozenset keys 339 B/seg, `_seg_manifolds` /
  `_manifold_segs` owner sets 375 B/seg, `_Segment` 100 B/seg, `_seg_lookup` 59 B/seg), `Point` objects
  on blast images 28 %, phantom preimages 15 %, walk memo + curve objects 2 %, rtree C memory ~0
  (only 4.2 k of 1.40 M segments at 1e-7, 18 k of 3.28 M at 1e-9, are in an rtree: blast images are query-only, yet every one of their
  segments keeps a `_Segment` + edge key + owner set), registry and topology < 0.3 %.
- Law cases (tiny: 1-22 k points, 0.1-0.8 s): the same picture, segment index 42-54 %.
- Time at scale: `blast` is 99 % of the series, and inside it the same-stability COUNT
  (`report_resolution(new_bridges=)` -> `diagnostics._candidate_pairs`, a sweep over EVERY unstable
  segment of the fixed point each blast) is 87 % (404 of 458 s on the 1e-7 series); the map is 0.07 s.
  Growth, crossings and topology are < 1 % there. `_candidate_pairs` is numpy, not linked-list work:
  arrays do not remove it.
- `merge_manifolds` ties between distinct nodes: 17 on each k28 law build, 2 on heteroclinic, 0
  elsewhere; all in growth, all with both nodes iterate-linked, keeper always the existing chain's
  (`head_1`). 0 refinement midpoints round onto an endpoint (~0.6 M refined). `_image_chains` never
  returned more than one chain on any LAW case or series; it returns two once each in the probe's
  scenarios `shared_straddle` and `nested_unfiltered` (N4 analyst, 2026-10-25) and in
  `test_blast_through_a_shared_straddle.py`.
- No walk ever passes THROUGH a BranchPoint (only roots); each unstable branch has one anchor bridge
  rooted AT its BranchPoint. (Corrected by the N5 analyst: a registered UNRESOLVED bracket can hold one
  -- k32_series id 503, scale id 415; "no bracket holds one" is true of the law cases only.)

Smells (data model):
- `BaseManifold._branch_forward/_branch_backward` (:540-620): walking THROUGH a BranchPoint returns the
  OTHER stability's first node in the same slot; legacy from when crossings were BranchPoints; never
  exercised (0 through-walks on six law cases).
- `BranchPoint.edists`, `forward_stretch_params`, `backward_stretch_params` (:91-97): never read.
- The `merge_manifolds` comment (:552-555) says the newly generated point has no iterate links; it
  always has `prev_iterate` by then (:366), so the tie rule's branches other than "keep head_1" are
  dead on every measured case.
- (FIXED by N1 on `refactor/slim-segment-index` @ `8ded634`, not merged at writing) `Tangle` kept
  full per-segment bookkeeping for query-only image curves (55 % of memory at scale). Still open:
  `populate_intersections(new_sids)` scans ALL `_intersecting_segments` per image curve (N5).
- (FIXED by N2 on `refactor/same-stability-count-new-only` @ `6bd37cf`, not merged at writing)
  `count_same_stability_crossings` built every box-overlap pair and filtered after (86 % of blast
  time). Still open: the narrow phase of `_hits` (`diagnostics.py:344-366`) runs UNBLOCKED over all
  candidate pairs (12.65 M at scale blast 13 -> ~1.55 GB traced transient), which sets the series'
  peak RSS; the rows, dedup and two argsorts over every segment remain (~1.1-1.5 s a blast).

Features (N0 branch `41056f8`, all verified by the TE and re-checked by the auditor):
- `stress/profile/profile_numerics.py` (README `stress/profile/README.md`): rerun from a clean
  worktree reproduces every `Point` / segment count EXACTLY at every blast of the 1e-7 series
  (plateau "no new points since blast 30", 2,177,786 Point objects) and the 1e-9 point (5,098,154),
  census / traced bytes within 0.27 %, peak RSS within 0.7 %. Full run ~10 min, 1e-9 ~13 min.
- `--plot-only --compare figures/numerics_profile/a8bb8a0` draws 5 of 6 figures with both series
  ("base" / "this run"); `stage_peak_rss.png` and the right (by-component) panel of
  `scaling_bytes_per_point.png` draw this run only (`profile_figures.py:195` `figure_stage_peaks`
  takes no base). N10 owns that fix (add `profile_figures.py` to N10's Files).
- `stress/profile/identity_probe.py` (D3 probe, 309 lines): `dump` (58 builds = 11 law cases + 47
  series steps, ~60-70 s, 126 MB) per checkout via `TANGLEPACK_REPO`, `compare` byte-equal -> exit 0
  `identical: 58 builds, 5776 arrays, 638 texts` at `a8bb8a0`. Readers: `curve_arrays`,
  `curve_extent`, `bracket_cdists`, `image_curves` (private `_image_curves`), `key_text`; registry
  fields, `bridge.id` and topology `describe()` are read directly in `numerics_state` /
  `topology_texts` (its Dev Note over-claims "every read"; fix at merge or in N1). It dumps per-point
  (x, y, cdist) only: NOT the point-level iterate links (`next_iterate`) nor phantom preimages, so a
  slice that changes those (N3, N4) is checked only through their downstream effects (registry,
  bridges, topology) unless it adds a reader. Teeth shown: one ulp, one text char, one src constant
  (419 diffs).

## Smells (crossing set)
- `IntersectionRegistry.reindex_from` (:270): the docstring says unmatched crossings get "fresh ids
  above every preserved one", the code (:309-314) gives the LOWEST unused ids, so an id an old
  crossing held (one that moved and was not re-matched) can name a DIFFERENT crossing after the
  recompute. Seen live after a densify recompute (k=2.5 p3: old 87 and old 88 = f^17(q0) swapped
  ids; p1: f^9 / f^10 became ids 1 / 2, the gaps of the build's pip candidates). `rebuild_bridges`
  carries `iterated` / `dropped_at` by id under an unchanged epoch, so a mis-carry is possible
  (believed; none observed). Matching is exact (`DEDUP_RTOL`): ANY crossing on a refined segment
  moves and gets a new id (2 per near-anchor densify, 47 of 414 at radius 1e-2).
- `Tangle._resolution_ratio` near an anchor: the stable polyline's head is straight chords at a
  geometric spacing set by growth, not by `area_cutoff` (k=2.5 p1 `s(0,0)`: 0, 1.98e-6, 1.1e-5,
  6.1e-5), never refined, so a deep crossing there scores ~ l_s |sin| / h_u (h_s = 0). This is what
  leaves `f^{m+1}(q0)` unresolved on k=2.5 p1 / p3 at 1e-7 (ratios 0.88 / 2.10) and why a finer cutoff
  (smaller h_u) buys depth (blast-steps-compose analyst, 2026-10-21; possibly over-strict when h_s = 0).
  Densify analyst (2026-10-21): refining ONLY the unstable segments within 2|f^m(q0) - P| of the
  orbit to h <= l_s1 / 6 resolves every such stall measured (k=2.5 p1 415: 0.877 -> 14.0, 3 points;
  p3 107: 2.10 -> 542, 7-12 points; k28 245, k32 503, nested_outer 349, inversion 238), and the
  detection then sits 1.3e-7 / 4.5e-8 from `f(f^m(q0))` (was 1.25e-6 / 9.3e-6). Depth then gains
  +1 (period 1) or +4 (p3) before the ZONE FILTER stops it (see loom.md).

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

## The segment index's bookkeeping (N1 analyst, 2026-10-24, `41056f8`, measured)

Probes in `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/7a98b0ff-df74-46bf-9769-0728d7d299d5/scratchpad/n1/an/`
(`instrument.py` monkeypatches the Tangle; `n1plugin2.py` runs it under pytest; `sim_census.py`
simulates the slim end state and runs the profiler's `census()`; `optd_compare.py` + `optd_plugin.py`
show the curve-vs-curve dedup is load-bearing). Lines at `41056f8`.

Map (who touches the per-segment records; all in `numerics/Tangle.py` unless named):
- Fields `__init__` :189-195 (`_edge_to_sid`, `_seg_lookup`, `_manifold_segs`, `_seg_manifolds`,
  `_intersecting_segments`, `_processed_pairs`), cleared `clear_all` :208-218. `_Segment.in_rtree`
  :104 is read only by `_remove_segment` :698.
- Writers: `add_manifolds` bulk :506-538 (branches; its own inline edge dedup :515-525),
  `_insert_segment` :580-643 (edge dedup :611-618 -> `_claim_segment`; new record :627-636; query
  :638-641), `_claim_segment` :575-578, `add_manifold` purge :478-481 (re-add of a known manifold).
- Readers: `_resolve_crossing_pair` :252-264 (processed-pair early return, stale check, both
  `_Segment`s, the u x s assert :264), `populate_intersections_for_manifold` :563-567 (filter by
  `_manifold_segs[curve]`), `add_manifold` :484/:491 (the new-segment COUNT),
  `release_manifold` :667-669 / `_remove_segment` :685-708, `add_manifolds` :506 (`if
  self._seg_lookup` = "index non-empty").
- Callers outside: `BridgeIterator.py:359` (image curve, query-only), `:362`;
  `TangleWorkbench.py:1349-1356` (recompute re-indexes `_image_curves`, drops a count-0 curve and
  releases it), `:1756` (`_register_bridge` releases a duplicate), `:1856` (`clear_bridges`);
  `Tangle.create_bridges` :967-968 claims every cut bridge's segments (`detect_crossings=False`).
  Tests: `test_tangle_index_and_orientation.py:58-90` (branch-cut bridges own their segments),
  `test_segment_index_bookkeeping.py:18-32`, `test_same_stability_discard.py:72-75`,
  `test_crossing_resolution.py:121`, `test_bridge_identity.py:257` (weakref: no bridge leaks).
  Tools: `stress/profile/profile_numerics.py:333-354` `_tangle_bytes` (the census reads the Tangle
  ONLY through these five fields: anything new on the Tangle is invisible to it),
  `stress/numerics/e02_gpu_end_to_end.py:135`, `e03_scaling.py:148` (`len(_seg_lookup)` metric).
  `Intersection.seg_ids` (:113) is read only for None-ness (`is_anchor` :238, `is_synthetic` :243).

How it behaves (verified on the 11 law cases, the 5 series in the probe's loop, the k=2.5 scale
case to its 1e-7 plateau, and the whole suite under instrumentation):
- The query path's edge dedup (`_insert_segment` :611) NEVER fires on the 58 probe builds or the
  scale case (0 of 4.7 M query edges; every image edge has an endpoint the blast just created). It
  fires only in tests: 3,794 blast-time duplicate edges (shared-straddle, unfiltered nested blasts,
  inversion bridge iterate, partial pieces, growth round), 49 of them against EARLIER IMAGE CURVES,
  plus 3,550 at recompute (reset=False 3,534; `test_a_growth_round_keeps_the_blast` 16).
- Every duplicate edge on the blast path has BOTH endpoints pre-existing, and no edge with both
  endpoints pre-existing is new (0 counter-examples over 81,672 (suite) + 2.27 M (corpus + scale) new blast edges): "query an image
  edge iff an endpoint was created by this iterate" reproduces the dedup on all data (believed in
  general; not by construction).
- Dropping the dedup against query-only records (dedup vs indexed only) CHANGES outputs on the
  unfiltered nested blast (outer x4, inner x2): partial pieces, bridges, an inert class's members,
  9 registry rows at 1e-11. The probe would not see it (0 firings there). The suite still passes.
- Recompute: 0 of 1,323 curve re-adds on the corpus drops a curve; the one drop in the suite
  (`test_growth_keeps_every_crossing_and_class[k28_two_blasts]`) is a curve whose every edge is a
  branch edge; no curve partly overlaps a branch.
- `release_manifold` releases 0 claims on every corpus build (every call follows `clear_all`);
  `_remove_segment` pops a query-only record once in the whole suite (reset=False test).
- Separators (`_insert_crossing_separator`) make claim-only NEW edges (registered under bridges
  only, never queried): 0 on the corpus, 43 in the suite.
- `populate` can filter on the add's NEW segment ids: at every call every unprocessed pair holds one.
  No stale pair is ever dropped. No two detections tie on the registration sort key (0 on the
  corpus), so segment-id VALUES decide nothing (the registry order is set by cdists).
- At the 1e-7 plateau (32 blasts): 1,397,859 segments, 4,226 indexed; 1,391,935 have 2 owners
  (curve + bridge). 1e-9 (29 blasts): 3,276,670, 18,256 indexed.
- Memory, simulated end state (census()): no query-only record: 1e-7 1,902.7 -> 861.8 MiB, segment
  index 1,045.4 -> 4.5 MiB; 1e-9 4,657.3 -> 2,011.6 MiB. A persistent set of `id()`-packed 128-bit
  edge ints instead costs 81-88 B/edge (+117 MiB at 1e-7, +252 MiB at 1e-9). A frozenset edge dict
  entry costs 350 B. Adding one attribute to `Point` in `__init__` grows its `__dict__` 144 -> 232 B
  (Python 3.10).
- Time (32 blasts at 1e-7, no diagnose, 57.6 s): query-only `add_manifold` 11.9 s, bridge claims
  1.3 s, populate 0.6 s; the key/dict/owner work is ~36 % of the query (bench on 200 curves).


## After N1 `slim-segment-index` (auditor, 2026-10-24, branch @ `8ded634`, base `41056f8`; lines are the branch's)

Read the section above as the BASE; this is what the branch changed (not merged at writing).
- Two kinds of curve (`Tangle` class docstring :148-176). INDEXED (branches, and bridges cut from
  one): `_Segment` + `_edge_to_sid` frozenset key + owners in `_manifold_segs` / `_seg_manifolds` +
  rtree entry, as before. QUERY-ONLY (blast image curves, bridges cut from them, partial pieces):
  one int per edge in `Tangle._query_edges` (`_edge(p0, p1)` :93-97, `(min id << 64) | max id`,
  same identity as the frozenset key), and a `_Segment` in `_seg_lookup` ONLY when its query paired
  (kept forever after, ~2,400 at scale; N5's cleanup). No owner, so `release_manifold(curve)` is a
  no-op for it (the call in `compute_intersections` stays for indexed curves).
- `add_manifold(curve, index_segments=False)` (query path :513-536, early return) skips an edge in
  `_query_edges` or in `_edge_to_sid` (asserting same stability :522-523), else adds the key, takes
  `next(Tangle._ids)`, queries `_rtrees[other]`. Returns `set[int]` of new sids on both paths
  (truthiness = the old count's meaning; `compute_intersections` drops a curve whose set is empty).
- `populate_intersections(new_sids)` (renamed from `populate_intersections_for_manifold`, :602) filters
  pairs on `pair & new_sids`; caller `BridgeIterator.py:359-362`.
- `_insert_segment` :665-668: an edge already in `_query_edges` returns None BEFORE the `_edge_to_sid`
  claim (fix A of review round 1): an indexed add over a query-only edge neither indexes nor queries
  it. The base ALSO claimed it (an owner on a non-rtree record); the branch does not. No observable
  effect over two reset=False rounds (TE). Returns None for any already-registered edge.
- `create_bridges` step 3 (:1020-1035, loop :1028): a bridge whose parent is a key of `_manifold_segs` claims
  (Q5); otherwise its edges (separators included) go into `_query_edges`. Corner case (believed, not
  observed): an indexed branch all of whose edges are already query-only edges never becomes a key of
  `_manifold_segs`, so bridges cut from it are treated as query-only.
- `_query_edges` is cleared only by `clear_all`; it never forgets on `reset=False` (Q1's accepted
  risk). `id()` reuse would alias keys if a curve's Points were ever freed mid-index (no path does).
- Measured at the 1e-7 k=2.5 plateau (TE paired runs, auditor re-read the JSONs): census 1914.0 ->
  990.6 MiB; segment-index band 1045.4 -> 122.0 MiB (query keys ~117, 6,584 records vs 1,397,859);
  peak RSS over baseline 3034.7 -> 2344.8 MiB (absolute 3151.5 -> 2462.1); blast 399.05 -> 399.57 s
  same day (flat; the committed 458.4 baseline is run drift); `index`+`cut` 10.70 -> 9.40 s. Census
  fell 923 MiB but peak RSS only 690: the peak is not set by held state alone (believed: transients
  of the same-stability count, N2's code).
- Open, older than N1 (Tangle Dev Note :31-42): the accumulating `reset=False` pass misses crossings
  a fresh pass finds (k28 2 blasts grown u+s: 16 vs 30; nested grown s: 48 vs 54; two rounds 24 vs 50,
  84 vs 98), identical on base and branch.
- Features: probe now 62 dumps (58 + `nested_unfiltered` / `shared_straddle`, each `.blasted` and
  `.recomputed`; these are where the image-over-image dedup decides, 15 diffs under "dedup on
  `_edge_to_sid` only"); probe and profiler ran blasts with `diagnose=False` / recorded no
  `same_stability_count` until N2, which added that reader (see the N2 auditor section). e02/e03 report `indexed_segments`; the stress
  report (`report/analysis.py:386`, `make_figures.py` `fig_e03_cost`, `make_numbers.py:232-233`) still
  reads e03's old `segments` key.
- Tests: `tests/unit/numerics/test_segment_index_bookkeeping.py` (9 tests: the paired-records
  invariant fails on the base by design; the image-over-image, separator, merged-curve and
  branch-over-query-only tests pass on base and branch).


## N2 analyst (2026-10-24, main `0b1db28`): the same-stability count, measured

- Path: `Blast.py:554` `report_resolution(fixed_points, new_bridges=added)` (the returned
  `ResolutionReport` is DISCARDED; `BlastResult` has no field for it, so any reader must wrap
  `TangleWorkbench.report_resolution`) -> `TangleWorkbench.py:1454` -> `diagnostics.py:200`
  `count_same_stability_crossings` -> `_curves` (:257, manifolds of the new bridges' fixed points,
  then `workbench.bridges` = full bridges in insertion order then partial pieces, minus any bridge
  sharing a node with a manifold) -> `_hits` (:314: `_segment_rows` per curve, a Python
  `get_point()` per node; first-owner edge dedup `np.unique(..., return_index)` :337-338; pairs
  :343; shared-node drop; narrow phase; `in_region` from `point = a[i] + t (b[i]-a[i])` :368).
- "New" is a CURVE flag (:223-227) carried to segments by the first-owner dedup: a segment is new
  iff the first curve listing its edge is a new bridge. A new bridge sharing a node with a manifold
  is no curve (3 on `nested_unfiltered`, 1 on `shared_straddle`); an edge a new bridge shares with
  an earlier old curve is old (1 segment, `nested_unfiltered` blast 5). Self-crossings of a new
  curve and new-vs-new pairs count; old-vs-old never.
- Orientation matters for exactness: `curve_pairs` are ORDERED `(owner[i], owner[j])` with `i` the
  earlier segment in the stable argsort on the sweep axis, and the axis is chosen over ALL segments
  (fewer pairs). Orienting "new first" keeps the WARNING/INFO split on every corpus hit blast but
  changes `curve_pairs` and hit points' last bits (6 of 6 corpus hit blasts).
- Cost at the 1e-7 k=2.5 series (32 blasts, in-process): the count 344.1 s; at the plateau a blast's
  count is ~20 s of which the pair generation is 18-19.6 s, rows 0.26, dedup 0.40, narrow <0.2,
  `_curves` walks 0.07 (profiler label `walk`, not `diagnose`). Late blasts have 1-50 k new segments
  of 1.39 M, yet 0.1-12.6 M candidate pairs hold a new one (segments in folds overlap ~100 boxes).
- A new-only enumeration (scratch `n2lab.py`: keep the axis choice; (A) each new segment's forward
  range; (B) each old segment's forward range restricted to new positions by `searchsorted`; then
  the other-axis test) gives the SAME pair set and orientation, so bitwise-equal hits, on all 65
  corpus blasts (law cases, `cases.SERIES`, both scenarios) and all 32 scale blasts; stage sum
  344 -> 63 s (+1.4 s walks). Floor per late blast ~1.1-1.5 s (rows, dedup, two global argsorts).
- `diagnose=True` changes no probe dump (k28_series, nested_outer_series, nested_unfiltered: 23
  builds identical against `diagnose=False`).
- Hits are rare in the probe corpus: 6 of 65 blasts (k28_series b8, k32_series b8-10,
  nested_outer_series b8-9), 113 k hits, NO self-crossing of a new curve. The scale series has
  hits on 23 of 30 blasts, 2,191 self hits.



## After N2 `same-stability-count-new-only` (auditor, 2026-10-24, branch @ `6bd37cf`, base `0b1db28`; lines are the branch's)

- `diagnostics._candidate_pairs` (:403-457): the axis choice and `counts` are over ALL segments (two
  stable argsorts, one per axis, :424-428); then per sorted position `p` the partners are
  `pool[start[p] : start[p] + count[p]]`: `new=None` -> `pool = arange`, `start = p + 1` (the base's
  arithmetic); with `new`, a new `p` keeps its whole forward range, an old `p` gets the new positions
  in `[p+1, p+1+counts[p])` by two `searchsorted`s. Orientation (first = earlier sorted position) and
  the per-first ascending order are the full sweep's, so hits are bitwise the base's (verified: 3,750
  harness calls, `scratchpad/aud2/eq.py`; probe 62 dumps; 32 scale blasts). Blocks are cut on
  `cumsum(count)` by `_PAIR_BLOCK` (the same rule, different boundaries).
- Measured (sequential same-day pair, 1e-7 scale): `diagnose` 360.3 -> 53.7 s, blast 415.8 -> 109.2 s,
  peak RSS over baseline 2396.2 -> 2374.0 MiB. A concurrent pair read +178.6 MiB: never run profiler
  pairs concurrently. At blast 13 the count sees 866,911 segments (160,710 new) and 12.65 M candidate
  pairs, 5.5 s, traced peak ~1.55 GB on base and branch: the narrow phase, not held state, sets peak.
- Features: the identity probe runs EVERY blast with `diagnose=True` and dumps `same_stability_counts`
  (one line per blast since the previous dump, partial-bridge addresses stripped; `record_counts`
  wraps `TangleWorkbench.report_resolution` because `Blast.py:554` discards the report); 65 blast
  lines, only 6 with hits and NO self-crossing of a new curve. `profile_numerics.py` scale rows keep
  `same_stability_count` (23 of 32 blasts with hits, 2,191 self hits): compare both in later slices.
- Tests: `test_the_count_is_the_full_listing_of_the_new_curves_on_random_curves[1..5]`
  (`test_same_stability_diagnostic.py`): random walks, image chains with shared edges, a Z-kinked new
  bridge, odd seeds rounded for ties; the only test that catches dropping a new curve's self-crossings.
  It misses an upper-end off-by-one in branch B (`hi_off`; seed 10 would catch it); the older C5 test
  `test_crossing_resolution.py::test_the_count_is_the_listing_of_the_new_curves` does.
- Identity-keyed sites N3 must keep valid (Point identity): `Tangle._edge` (:93-97, `id()` pair),
  `diagnostics._curves` (:274-279) and `_segment_rows` (:386), `ManifoldMachine._image_chains`
  (`reached.add(id(node))` :461).

## Points, phantoms and the store (N3 analyst, 2026-10-24, main `138b34c`, measured)

Scratch: `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/7a98b0ff-df74-46bf-9769-0728d7d299d5/scratchpad/n3/` (`instrument.py` + `probe_cases.py` /
`probe_series.py` / `probe_scale.py`: creation site and fate of every `Point`; `tie_probe.py`: per-tie
record; `proto/`: a full prototype of the store on a copy of `138b34c`; `prof_{base,proto}/`: a
sequential scale profile pair).

- Where points are made (`Point(` in src): `ManifoldMachine._iterate_manifold` :359-362 (images),
  `_refine_layer` :780 (refined), `_cache_preiterate` :830 (phantom), `ManifoldInitializer`
  :178/:187 (fictitious chain), :225/:246 (seed pair), `Tangle._insert_crossing_separator` :1071
  (separator; 0 created on every law case, series, scenario and the scale case -- tests only).
- Who reads a phantom: only `_get_preiterate` (:808-811, for `pre_mid` :738) and the merge tie rule
  (:561-568). A phantom's link is ONE-WAY (`next_iterate` of the phantom is None). Every refined point
  has exactly one phantom, on its PREITERATE side (unstable: `prev_iterate` = its preimage; stable:
  `next_iterate` = its IMAGE; k10: 1270 + 187 of 1457); no link in the other direction is ever a phantom.
  (Calling it the "preimage side" on a stable curve is what made N3's prose wrong in three places.)
- The tie rule: 17 ties on k28 (one / two blasts, k28_series), 2 on heteroclinic, 2 on k32_series, 0
  elsewhere (scale case included); keeper always `head_1`. `head_1` holds a phantom in 15/17 (k28) and
  2/2, but it ALSO holds a real image-direction link in every tie: no tie depends on counting the
  phantom (a mutant that ignores phantom rows leaves the probe identical). Flipping the keeper moves
  7 probe items (teeth).
- No `Point` ever links to a `BranchPoint` (or the reverse) in library code (0 on 11 law cases, 7
  series / scenarios, scale); only `test_linked_list.py`'s four `[branch_point]` ids do.
- Costs (tracemalloc / census): today a curve `Point` 392-408 B, a phantom 368-384 B. A slotted
  `Point` over a row: 72 B object (+8 with a `__weakref__` slot) + 28 B row int + 8 B node-table entry;
  a row 40 B (xy 16, cdist 8, image 8, preimage 8). So ~156 B per real point, 48 B per phantom.
- Prototype (`proto/`, store owned by the workbench, handed to the initializer; the machine appends to
  the store of the curve's own points; separators to `p1.store`; phantoms rows with NaN cdist and no
  `Point`; a node table `store.nodes` materialising a row's `Point` on demand): probe `identical: 62
  builds, 6024 arrays, 744 texts`; suite 1222 pass / 43 skip / 1 xfail (C3) with exactly the four
  `test_linked_list[branch_point]` failures (+1, `test_collect_rejects_a_broken_iterate_link`, if
  `BasePoint.__slots__ = ()`: its `__class__` swap needs a `__weakref__` slot); ties equal base/proto.
  Scale 1e-7 pair (base then proto, same day): census 990.6 -> 409.2 MiB (`points_image` 542.5 ->
  143.6, `points_phantom` 285.0 -> 0, store 104.0 incl. 4 % capacity slack and the 17 MiB node table),
  peak RSS over baseline 2372.2 -> 1819.5 MiB (both at blast 13), blast 105.8 -> 110.4 s (walk 17.6 ->
  19.2, infer_iterates 8.2 -> 10.3, diagnose 51.5 -> 53.4, refine 8.8 -> 8.0), same-stability count and
  `pip_depth` equal on all 32 blasts.
- Per-node property reads cost 4-7x an attribute (get_point view 122 ns vs 33, cdist 85 vs 19,
  next_iterate 97 vs 17); a vectorised row gather is 25 ns/node against 338 for today's `vstack`.
- Smell: `refine_manifold` evaluates `manifold.get_point_array()` as a `logger.debug` ARGUMENT (:634;
  :636 after N3)
  on every call, a whole walk even with DEBUG off.
- (FIXED by N3's profiler commit `0acf066`) The profiler's census read `p.__dict__` and classed a
  Point as a phantom by `cdist is None`; it now sizes a slotted point and reports `point_store` +
  `phantom_rows` (see "After N3").

## After N3 `point-store` (auditor, 2026-10-24, branch @ `5ce5685`, base `138b34c`; lines are the branch's)

- Map: `numerics/PointStore.py` (new, 170 lines): `PointStore` with private columns `_xy (cap, 2)`,
  `_cdist` (NaN = none), `_image` / `_preimage` (int64, -1 = none), `_n`; `append(xy, cdist) -> rows`
  (vectorised, `cdist=None` -> NaN), `append_one`, `nodes` (row -> `Point` or None), `node(row)`
  (materialises a phantom's Point once; no library path calls it on any case), `xy` (read-only alias,
  rebuilt by `_reserve` on each x1.5 resize), `nbytes`; `loose_store()` = one process-wide store for
  hand-built points (tests, notebooks), never freed. NOT exported (D4 API is N7's).
- `Point` (`Point.py`): `__slots__ = ("store", "row", "forward", "backward", "stretch_param")`;
  `Point(x, y, cdist=None, stretch_param=None, *, store=None)` appends one row; `Point._over(store, row,
  sp)` wraps a batch row; both register in `store.nodes` (so one persistent Point per real row: `id()`
  sites `Tangle._edge`, `diagnostics._curves` / `_segment_rows`, `_image_chains` untouched and valid,
  and an `id()` is never reused). `x` / `y` read-only `np.float64`, `cdist` float with setter, all three
  read a stored NaN as None (a base `nan` read `nan`: no library reader depends on it).
  `next_iterate` / `prev_iterate` are properties over the columns; `_row_of` raises `TypeError` for a
  non-Point or another store's Point. `get_point()` = read-only row view (Q9; a row's xy is never
  rewritten, so a pre-resize view stays right). `edist` and `_coords` gone. `BasePoint.__slots__ = ()`
  (BranchPoint, unslotted, keeps its `__dict__`), `BasePoint.has_iterate_link()` (attribute path) /
  `Point.has_iterate_link()` (columns: a phantom row counts).
- `ManifoldMachine`: `_iterate_manifold` (:353-377) one `store.append` per step into
  `non_iterated_points[0].store`, a per-step same-store guard (`TypeError`, :360), the "already exists"
  `ValueError` AFTER the append (:362 vs :370-371: dead rows on that non-resumable error path), then two
  vectorised column writes. Tie rule (:569-570) `head_k.has_iterate_link()`, no isinstance.
  `_refine_layer` (:777-782) appends the layer's points, then their phantom rows (`cdist=None`), and a
  one-way column write on the preiterate side. `_get_preiterate` (:799-814) returns coordinates (a
  WRITABLE view of `store._xy` -- reviewer (c)) or None; BranchPoint keeps the attribute read.
  `_cache_preiterate` deleted; no `Point(` left in the file. The machine reads the store's private
  columns (accepted until public names exist).
- Ownership: `TangleWorkbench.point_store` (`TangleWorkbench.py:254`) handed to
  `ManifoldInitializer(..., store=)` (a fresh store when omitted: the four hand-wired scripts
  `henon_basic_stable_and_unstable`, `henon_initial_3_points`, `henon_initial_intersections`,
  `henon_period_3_manifolds` build their own initializer + machine and still run, auditor-verified);
  separators `Point(..., store=p1.store)` (`Tangle.py:1071`).
- Rows on no curve: phantoms (778,302 at the 1e-7 plateau; k10 1457), the initializer's detached
  fictitious tail (one per seeded chain, a Point with NaN cdist holding one link: k10 2, heteroclinic /
  nested 4; the only rows the base would have freed), merge-dropped duplicates (alive in the base too,
  named by their preimage's link; `off_curve_rows` = 0 at the plateau).
- Measured (1e-7 scale, sequential same-day pair, coder): census 990.6 -> 398.5 MiB; `points_image`
  408.0 -> 100.0 B/object; store 104.0 MiB = 86.7 columns at capacity 2,273,629 x 40 B + 17.2 node
  list = 50.1 B/row; all point bytes 829.6 -> 237.4 MiB (~151 B per real point all-in); peak RSS over
  baseline 2338.8 -> 1883.6 MiB (blast 13; held after blast 13 767.8 -> 356.8); blast 106.7 -> 109.7 s
  (`walk` 17.9 -> 18.9, `infer_iterates` 8.2 -> 10.3, `diagnose` 52.1 -> 53.1, `refine` 8.7 -> 8.0,
  `iterate` 0.86 -> 0.37). Per-blast count / `pip_depth` / `points_created` equal on all 32.
- Tests: `tests/unit/numerics/test_point_store.py` (8 tests; the phantom-only tie test runs on the
  base too, the rest `importorskip` the store): mutants A (phantoms ignored), C (resize drops rows),
  D (phantoms materialised), E (phantom link both ways) each fail a test. The rewritten
  `test_collect_rejects_a_broken_iterate_link` monkeypatches `BaseManifold._iterate_of`.
- Smells: the store's columns are private but read by `Point` and `ManifoldMachine` (name them in N4);
  `_get_preiterate`'s writable view; the error-path dead rows; `profile_figures.MEMORY_GROUPS` (eight
  fixed slots) has no `point_store` group, so `memory_components.png` omits 104 MiB on the branch;
  CLAUDE.md's data-model paragraph does not mention the store.

## Curves, splices and shared nodes (N4 analyst, 2026-10-25, main `8fc3933`, measured)

(Describes the BASE `8fc3933`: the links, walkers and `_insert_point_geometrically` named here are
gone on the N4 branch; see "After N4" below.)

Scratch: `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/7a98b0ff-df74-46bf-9769-0728d7d299d5/scratchpad/n4/` (`instr.py` instrumentation:
every splice, merge, refine layer, `_image_chains`, `_outer_neighbours`, and a per-operation walk of
every live curve; `structure_probe.py <case|series|scenario> [--light]`; `n4plugin.py` = the same as a
pytest plugin; `scale_structure.py`; `beyond_tail.py`; `skipped_probe.py`; `stretch_probe.py`; `micro.py`).

- How it works today: a curve is `root` + `tail` over ONE physical list; a walk follows `forward`
  (unstable) / `backward` (stable). Restricted to live-curve nodes the links form DISJOINT PATHS: 0
  nodes with two successors or predecessors (BranchPoint excluded), 0 `walk_back` disagreements, 0
  unreached tails, on 11 law cases + 5 series + the probe's 2 scenarios at every operation boundary,
  and at the 1e-7 scale plateau. So "one int array per path, every curve a view (path, root, tail)"
  is exact by construction (the whole suite under `n4plugin.py` too: 0 conflicts, 0 backlink
  disagreements, 22,711 merges). On the probe corpus bridges, partial pieces and image curves never
  change after they are made; in the suite, growth AFTER a blast refines inside cut bridges and an
  image curve (`test_growth_after_a_blast`, `test_session_caches[growth]`,
  `test_generation_and_caches`: 7 curves) and absorbs an image curve as a second chain
  (`[4396, 1133]`), and separators splice into branches (`p3_16_12`, `inversion_period2_10`,
  `test_pip_orbit_cut`, `test_crossings_sharing_a_segment_are_bracketed`: 24; N3's "tests only" was
  about the probe corpus). So a per-bridge COPY would go stale; a view does not.
- Lists run past curve ends and the code reads there: trimmed stable branches (462-974 nodes past
  the tail); in `shared_straddle` / `nested_unfiltered` unstable branches run 63 / 208 nodes past
  their tails into blast images (dumped by the probe, which reads the WHOLE list), image curves start
  on branch lists (2 / 6) and run past their tails (40 / 182). Refinement's left neighbour of a merged
  image root lies outside it in every growth refine (the `BranchPoint` itself on k > 1 chains, 44-57
  reads per build); `_outer_neighbours` reads past trimmed stable tails and image-curve ends. A
  per-curve COPY of ids loses these; a view over the path keeps them.
- `merge_manifolds`: equal to `searchsorted(c1, c2, side="right")` + "a new point tying an existing
  one is left off" on all 1,078 merges (1,059 probe builds + 19 scale; 22,711 in the suite with
  1,118 tie values), incl. 55 tie values and 146
  new points interleaved inside the existing run (28 merges; k28 26 interior inserts per build).
  Never fired anywhere: (a) tie, neither linked; (b) tie, only the new linked (the EXISTING point
  leaves the curve); (c) root tie (both kept); (d) same node in both lists; (e) internal ties; (f) the
  existing run exhausted while its list continues (`only_forward` would orphan it); (g) new root
  before an existing root that has a predecessor; (h) both lists ending on a tie (splices None). A
  pure "drop on equal" differs from today on (b), (c), (e).
- `_refine_layer`: pairs always adjacent; every left / right read consistent with the path (0 of
  3.77 M pairs: 18 probe builds + scale); `refine_manifold`'s returned `modified_points` is never read.
- `_image_chains`: one chain everywhere except once in each scenario (two chains); an image on no
  list starting a run is its own chain (the initializer's fictitious points, 2-6 per build); off-list
  images inside a run (merge-dropped duplicates, 9 on k28, stale links) are skipped.
- Blast images at scale are ALL new (0 of 706 image curves hold an existing image; 0 blast merges);
  blast merges, path concatenation and multi-chain images occur only in the scenarios and in
  `test_blast_through_a_shared_straddle.py`.
- Counts at the 1e-7 scale plateau (32 blasts): 2 branches (4,228 nodes root..tail), 1,654 bridges
  (1,395,483 node references), 706 image curves (1,394,339 nodes), 1,398,566 distinct curve nodes, of
  which 1,393,589 lie on 2-3 curves; store 2,177,788 rows.
- `stretch_param`: every curve's points carry exactly the curve's value (7 law cases); readers of the
  per-point value: `ManifoldMachine.py:291` and `tests/helpers/invariants.py:283` only.
- Readers of `forward` / `backward` in `src`: the engine only (`BaseManifold`, `ManifoldMachine`,
  `ManifoldInitializer`, `ManifoldView`, `Tangle` `_outer_neighbours` / `_segments_of` / separator,
  `TangleWorkbench` :716, :775-789, :2090, `BranchPoint` inserts writing the Point's links); every
  other module reads nodes through `get_point_array(return_nodes=True)`. Outside `src`: 20 test files
  (22 test functions + hand-built-list helpers, `tests/helpers/polylines.py` among them), the
  profiler's census `_walk_ids`, `tests/_tools/resolution_probe.py`, three scripts + one notebook
  calling `_find_tail()`.
- Time (N3 branch profile, 1e-7): blast `walk` 18.9 s = Blast reads 7.0 + separation 5.3 + diagnose
  3.8 + iterate 2.5 + topology 1.2; `infer_iterates` has NO walk under it (bracket reads: N5's). A
  walk is ~627 ns a node with its `vstack`, a gather ~5 ns, a node list from ids ~27 ns (`micro.py`).
- Smells: `merge_manifolds` comment `:559-563` false (the new image already holds its preimage);
  branch (h) above is a latent crash; `refine_manifold :636` debug argument walks every call;
  `modified_points` dead; the profiler census sizes curve `__dict__`s but not arrays held in them (a
  future `ids` array would be invisible).

## After N4 `curves-on-order-arrays` (auditor, 2026-10-25, branch @ `ab7bd20`, base `8fc3933`; lines are the branch's)

- Map: `PointStore.py` (406 lines) adds `Chain` (slots `store, index, ids, stretch_param, anchor,
  version`; `node_at(pos)` = anchor at -1, None off the ends; `insert(positions, rows)` one `np.insert`
  + renumber from the first position, asserts rows on no chain; `absorb(other, ids, start)` takes the
  merged order, retires `other` (empty ids, keeps its index), carries its anchor, marks dropped rows
  `_position = LEFT_OFF` (-2)), the store columns `_chain` / `_position` (int32, -1 off every chain),
  `chain_of(row)`, `position_of(row)`, `new_chain(rows, stretch)`, public read-only per-access views
  `xy` / `cdist` / `image` / `preimage` (`[:n]`), writers `link`, `link_phantoms`, `set_cdist`; a private
  `_xy_read` alias for `Point.get_point`. `Point.__slots__ = ("store", "row")`; `stretch_param` a
  read-only property of the chain (None off a chain); `Point(x, y, cdist, *, store=None)`; Point reads
  the PRIVATE columns (a public view per node read costs ~330 vs ~50 ns). `BranchPoint` is never a row:
  a chain hanging off one names it `anchor`; `insert_point_forward/backward(node, b)` set the slot and
  `chain.anchor`, `ValueError` if the slot holds another chain's point.
- How it works: a curve (`BaseManifold`, `Bridge`, partial piece, image curve) is a VIEW: `_span`
  (`BaseManifold.py:185`) resolves `(head, chain, lo, hi)` from the root's row (the slot's point for a
  BranchPoint root) on every read, nothing cached; `final_node` off the chain or before the root reads to
  the chain's end; `_rows` (:227) gives the row slice (a lone off-chain root = its own row);
  `_walk_nodes` (:240) = head + `store.nodes[rows]`; `_collect` (:307) gathers `store.xy[rows]` /
  `store.cdist[rows]` and masks iterates over the link column (`_iterate_rows`), the BranchPoint head
  read per node. `get_cdist_array` and `_walk_nodes(stop_at_final=False)` read to the chain's end;
  `_walk_nodes()` / `(None)` reads to the TAIL (measured 176 vs 232 on k10's trimmed stable branch; the
  docstring at :248 says otherwise). `_find_tail` = the chain's last row.
- `ManifoldMachine`: `_iterate_manifold` appends a step's images, `store.link`, a new chain via
  `_manifold_from_points`. `_image_chains` (:423) reads `store._chain` / `_position` (private) and
  cdists: a run from image i takes every later image on the same chain at a position >= its own with
  cdist <= the last image's; an off-chain start becomes a one-row chain, an off-chain image inside a run
  is skipped, a LEFT_OFF start asserts; a BranchPoint first image (anchor bridge) starts at its slot.
  `merge_manifolds` (:516): vectorised `np.insert(ids_1, searchsorted(c1, c2, "right"), ids_2)` when the
  runs are on DIFFERENT chains and no cdist ties (:557); else the base zipper on positions (branches
  (a)-(e), the `head_1`/`head_2` debug line), then `first.absorb(second, prefix + merged + suffix,
  start)`; (f) `done_hi == len(done)`, (g) `second_lo == 0 and second.anchor is None`, (h) both runs
  exhausted on a tie: `AssertionError` naming both cdist ranges. `refine_manifold` (:621) returns None;
  layer 0 = pair positions `arange(lo - has_head, hi - 1)`; `_refine_layer` (:680) gathers neighbours
  from the chain (anchor = left of position 0, NaN left of the anchor and past the end), one
  `chain.insert(refined + 1, rows)` per layer, re-points the anchor's slot when it inserts at 0; next
  layer `(q + k, q + k + 1)`.
- `Tangle`: `_segments_of` = `zip(nodes, nodes[1:])` over `_walk_nodes()`; `_outer_neighbours` =
  `chain.node_at(pos -+ 1)` (a BranchPoint end has none); separator = `chain.insert` after `after`, or at
  0 plus `insert_point_forward` after a BranchPoint. `grow_until_turnaround` = `points[-1] -
  points[-2]` of `get_point_array()`.
- Invariants (asserted): a row on at most one chain (`Chain.insert`); equal stretch and at most one
  anchor on absorb; merge (f)-(h); `_manifold_from_points` rows on no chain; cdist non-decreasing after
  every growth step (`_iterate_manifold`). Held by construction: the columns agree with every chain
  (shadow-checked 0 violations at 1,006 probe + 11,186 suite boundaries, coder + TE); a curve's root and
  tail on one chain, root not after tail. A view breaks only if an END leaves its chain: the merge's
  keeper branch (b) on an existing point (never reached).
- Measured (1e-7 scale, two sequential pairs): blast 113.1 / 112.2 -> 85.0 s; `walk` 19.2 -> 1.3,
  `refine` 8.2 -> 1.0, `infer_iterates` 10.3 unchanged (bracket reads: N5), `index` 9.3 -> 8.9;
  census 398.5 -> 394.6 MiB (`points_image` 133.0 -> 101.1, `point_store` 104.0 -> 121.3 with the two
  int32 columns at capacity 2,273,629, `curve_ids` 10.8, 708 live chains); peak RSS over baseline
  -70 MiB (coder, dirty tree) / -30 MiB (TE, clean). Chains: k10 19 total / 2 live, k28 two blasts 23 /
  4 live (retired chains stay as empty objects). Probe identical, also against N3's pre-N3 dump.
- Tests: `test_order_arrays.py` (27 ids: chain invariants on a grown / blasted / regrown k10, (f)-(h)
  raises for both runs, slot hanging, separator and refine at the anchor, one-node curves, the iterate
  list), `test_merge_branches.py` (6 ids: tie branches (a)-(e) + the anchor-bridge image, all passing on the
  base too, auditor-verified), `test_point_store.py` (+3), `tests/helpers/polylines.py` builds through
  `_manifold_from_points`.
- Smells: `_walk_nodes` / `_span` docstrings say None reads to the chain's end (it reads to the tail);
  merge docstring omits "different chains"; `test_every_getter_matches_an_independent_walk` (F3) has
  an oracle sharing `_span`; `_image_chains` reads private `_chain` / `_position` (no public view);
  `Chain.version` written, never read (N6 memo); `refine_manifold(final_node=)` unused;
  `BranchPoint.insert_point_*` on an unchained node -> `AttributeError`; `BranchPoint.py:5-12, 29-37`,
  `BasePoint`, `gpu.py`, `Trellis.py:8,16`, `design_note.md:82`, `tests/helpers/invariants.py:24,301`
  still say "linked list"; `trim_stable_at` still scans `_walk_nodes`. Segments and brackets still hold
  `Point` objects and `Tangle._edge` keys by `id()` (`Tangle.py:93-97`): N5, whose hard case is a root
  segment starting at a BranchPoint (no row until N8).

## Segments, brackets and edge keys (N5 analyst, 2026-10-25, main `01e991a`, measured)

N4 is merged: the "After N4" section above describes `HEAD` (its smell list's F3, BranchPoint docstrings
etc. still open). Scratch: `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/7a98b0ff-df74-46bf-9769-0728d7d299d5/scratchpad/n5/`
(`n5instr.py` corpus instrumentation + `n5plugin.py` suite plugin; `infer_breakdown.py`; `batched_infer.py`
and `batched_query.py` shadow prototypes; `bp_bracket.py`, `widest_bp.py`, `swap_order.py`).

- Map (lines at `01e991a`): bracket capture `Tangle._resolve_crossing_pair` :323-324 (`(seg.p0,
  seg.p0_seg1)` of both `_Segment`s); readers `create_bridges` :914/:928 (+ anchor pseudo-bracket
  `(nodes[0], nodes[1])` :941-943, separator test `previous[2] is p0 and previous[3] is p1` :967,
  `_insert_crossing_separator` :1024 with the `isinstance(after, BranchPoint)` slot path :1060),
  `IterateInference._brackets` :67, `_widest_stable_brackets` :206 (EVERY detection, unresolved ones
  included), `_bracketed_image` :220 (listcomp :243 calls `_brackets(y)` TWICE per candidate), `_map_gated`
  :281/:290 (unstable bracket coords), `IntersectionRegistry._polyline_gap` :48 (identity) via `_runs`
  :466. Outside src: `tests/helpers/laws.py` :578-579, :627-661 (`_brackets`, `_passes_map_gate`; used by
  `check_area_along_iterate_links`, `check_iterate_link_scales_by_beta`, `check_iterate_table_complete`
  which reads EVERY detected crossing's brackets), `test_iterate_inference.py:31-32`,
  `test_growth_after_a_blast.py:31`, `tests/_tools/resolution_probe.py:50-54, 63, 69, 133` (also wraps
  `Tangle._resolution_ratio` reading `seg.p0`), `identity_probe.bracket_cdists` :103-110, `e07_breaking.py:343-354`.
  `id()` edge keys: `_edge` :93 (int) and frozenset keys :519, :570, :663, :766.
- How it behaves (corpus = 11 law cases, 5 series, the probe's 2 scenarios; suite = 16.0 M edges):
  - id-keys -> row-keys (`-1` for a BranchPoint end) is a BIJECTION per Tangle: 0 conflicts either way,
    corpus and suite; every Tangle sees ONE store (polyline tests: the loose store, not
    `workbench.point_store`). One `Point` per row and rows never reused make it so by construction.
  - Root segments: 2-8 indexed per build; root x root detections (`_do_segments_intersect` counts the
    shared BranchPoint endpoint as a crossing) resolve to EXACTLY (0, 0) and dedup into the anchor
    (2,342 in the suite). BUT a deep blast image can cross the STABLE root segment: a registered,
    unresolved detection whose STABLE bracket starts at the BranchPoint exists on k32_series (id 503,
    blast 9 on, ratio 1.48) and on the 1e-7 scale series (id 415, blast 9 to the plateau, ratio 0.88);
    in the suite also `test_deep_blast_sliver_beyond_pip`, `test_an_order_flip_gives_no_verdict` (x2),
    plus a hand-built UNSTABLE bracket at the BranchPoint cut at in
    `test_a_separator_after_the_periodic_point_becomes_the_first_point`. N0's "0" holds on the LAW CASES
    only. Readers that see it: `_widest_stable_brackets`, the probe's `bracket_cdists`, the law
    `check_iterate_table_complete` (`BranchPoint.get_cdist` = 0.0). It is never in N, never linked, never
    the widest bracket (k32: widest 0.0511 is an ordinary one). A `-1` row read as `cdist[-1]` would
    silently take the LAST row.
  - `_polyline_gap`: identity == row equality on 431,785 calls (suite; 334 one-shared, 1,856
    two-shared), 0 differences. Brackets are adjacent rows of one chain at every corpus cut (0 separators
    on the corpus; 24 in the suite, 1 after a BranchPoint).
  - `_find_true_intersection(seg1, seg2)` takes its operand order from `tuple(frozenset(pair))`, i.e.
    from the sid VALUES (and insertion order when they collide mod 8). Forcing (u, s) or (s, u) moves 3
    / 8 of nested's 40 registry rows by <= 4e-16 relative (ratio 1.6e-14): sid values DO decide float
    bits (N1's "decide nothing" was about registration order only). Keep the sid counter and build
    `frozenset((new_sid, cand_id))` in that order.
  - `populate_intersections`' scan is 0.23 s over 706 calls at scale: `_intersecting_segments` holds only
    CROSSING pairs (2,362 at the plateau, = `_processed_pairs` at rest), not segments.
- Time (scale 1e-7, 32 blasts, diagnose off, under load; same-process ratios): `infer_iterates` 11.2 s =
  `_bracketed_image` 9.4 (22,875 searches, ~140 candidates each; 7.44 M `_brackets` calls) +
  `_widest_stable_brackets` 1.7 (738 calls over the whole registry) + `_map_gated` 0.06; 8.3 s of it is
  the 32 registry-wide passes re-searching every unlinked source. A batched version (padded brackets
  read once per crossing from its rows, per-branch numpy arrays per `_orders_version`, one vectorised
  mask per search) gave the SAME image lists in the same order on all 18 builds + scale (0 mismatches,
  22,883 searches) at 0.91 s against 11.2. The query-only `add_manifold` on rows (one public `store.xy`
  view per curve, row edge keys, ONE `rtree.intersection_v` per curve, `intersects` on candidates only)
  found the same new edges and pairs (0 mismatches; scale, k28, k32, nested_outer) at 0.66 s against 9.84
  (scale). Micro: a property read 111 ns, private column 63, public view per read 358, view once then
  scalar 66, gather 0.9 ns/row.
- Memory: the query-edge set with row keys `((lo+1) << 32) | (hi+1)` (32 B ints) instead of id-keys
  (40 B) is 106.5 against 117.1 MiB for 1.39 M edges (synthetic, tracemalloc).

## After N5a `brackets-on-rows` (auditor, 2026-10-25, branch @ `9b72985`, base `01e991a`; lines are the branch's)

Scratch: `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/7a98b0ff-df74-46bf-9769-0728d7d299d5/scratchpad/aud5/`
(`pairs.py` five-pair table, `probe_m4` / `probe_m8` mutant dumps, `mut/` trap-(a) mutant).

- Map:
  - `Intersection.py:32-38` `SegmentEnd = Union[int, BranchPoint]`, `Bracket`; `:178` `segment_ends`
    (an end is tested with `isinstance(end, int)`: importing BranchPoint there is circular); `:205`
    `brackets` (`cached_property`; None iff `unstable_segment is None`).
  - `Tangle.py:310-331` the capture (`node.row` or the BranchPoint; asserts ONE store for the four
    ends); `:285` `tuple(seg_id_pair)` and the `_find_true_intersection` operand order untouched;
    `:934-940` the cut turns rows into nodes once (`store.node`), the rest of the cut (separator test
    `previous[2] is p0 ...`, `_insert_crossing_separator`, `_cdist_between`) is the base's.
  - `IntersectionRegistry.py:54` `_polyline_gap` (row equality; returns None before any array when no
    end is shared; asserts one store); `:666` `stable_brackets(branch_key)` -> `(cdists list, int64 ids,
    (m, 4) brackets with NaN rows)`, rebuilt per `_orders_version` in `_stable_branch_buckets`.
  - `IterateInference.py:67` `_brackets` pads `crossing.brackets`; `:74` `_padded` (numpy, same IEEE
    ops); `:207` `_widest_stable_brackets` (reads the PRIVATE `registry._stable_branch_buckets()`,
    accepted); `:221` `_bracketed_image` (same bisect, one mask, survivors in slice order, then `j in
    resolved` and `manifold_a_key`); `:271` `_map_gated` (coords via `segment_ends`, one call per
    decided match).
- Invariants relied on: a row's cdist is never written after detection (writers: BasePoint init,
  `ManifoldInitializer:275`, `FixedPointSolver:168`) -- the `brackets` memo depends on it, NOTHING
  asserts it; one `Point` per row for the store's life (node identity == row equality in the cut).
- Measured (five sequential pairs, scale 1e-7, 32 blasts, all under a concurrent user game): blast
  85.9-91.4 -> 75.9-82.3 s; `infer_iterates` 10.5-11.7 -> 1.10-1.34; `resolved_ids` -0.62 every pair;
  `index`/`cut`/`resolve`/`separation` flat; census 394.6 -> 395.4 MiB (registry +0.83); `segment_index`
  122.0 both; peak RSS over baseline = noise (base alone 1805-1883 MiB at blast 13).
- Found by the audit: the probe is IDENTICAL with the two stable comparisons dropped (m4) and with the
  candidates reversed (m8): on the 18 probe builds the stable half of the bracket rule never changes a
  decision (the map gate and the unstable half decide); only
  `test_an_image_whose_stable_bracket_misses_the_carried_one_is_no_image` and the stable-order assert
  in `test_a_multi_match_is_recorded_and_a_second_pass_changes_nothing` pin them.
- Smells: `Tangle.py:274` docstring still says "bracketing points"; `IterateInference.py:241` window
  comment elliptical; `segment_ends` always builds both arrays (brackets reads only cdists,
  `_map_gated` only coords); the brackets memo's immutability assumption is unasserted.
- Left for N5b: `_Segment` holds Points (`:110`); `_edge` (`:93`) keys `_query_edges` by `id()`;
  `_edge_to_sid` keys by `frozenset((id, id))` (`:525, 576, 669, 772`); EVERY new query edge takes a sid
  (`:531`) paired or not, and pairs are `frozenset((seg.id, cand_id))` (`:534, 605, 705`) -- sid values
  decide float bits through `tuple(seg_id_pair)` (`:285`); the cut's `_query_edges.update` (`:1031`).
- Features: brackets on rows -- verified by the probe (62 dumps), the trace (TE/reviewer) and
  `test_a_bracket_end_at_the_periodic_point_reads_its_point_and_zero` (fails under a -1 sentinel,
  auditor-verified).

## Segments and edge keys on rows (N5b analyst, 2026-10-26, main `1c25afb`, measured)

N5a is merged: the "After N5a" section above describes `HEAD` (its line numbers are `9b72985`'s; at
`1c25afb` they shift by one from the query path on). Scratch:
`/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/7a98b0ff-df74-46bf-9769-0728d7d299d5/scratchpad/n5b/`
(`base/`, `new/` = `git archive HEAD` copies; `tools/patch_new.py` writes the prototype into `new/`;
`tools/trace.py` the call-by-call trace; `tools/cmp_trace.py`; `tools/dump_all.sh` / `run_all.sh`;
mutants `new_m1` (pair built `(cand, new)`), `new_m2` (sid only for a paired edge), `new_m3` (candidates
reversed), `new_m4` (cdists as np.float64); `prof/` two profiler pairs; `prototype.diff`).

- Map (lines at `1c25afb`, `numerics/Tangle.py`): `_edge` :93 (`(min id << 64) | max id`); `_Segment`
  :110 (`id, manifold, p0, p0_seg1, in_rtree`; `bounds` :126, `intersects` :137); `Tangle._ids =
  count(0)` :186; `_edge_to_sid` :217 (frozenset-of-`id()` keys built at :526, :577, :670, :773);
  query-only path :519-541 (`next(Tangle._ids)` :532 for EVERY new edge, `frozenset((seg.id, cand_id))`
  :535); bulk `add_manifolds` :556 (pairs :606); `_insert_segment` :640 (key :670-671, sid :688, pairs
  :706, the `else seg.id` branch at :688 is dead: `_segments_of` always yields `id=None`);
  `_remove_segment` :736 (key :773); `_segments_of` :775 (one `_Segment` per edge from `_walk_nodes`);
  the cut's `_query_edges.update` :1031-1033. `_resolve_crossing_pair` :266 (`tuple(seg_id_pair)` :286,
  the bracket capture :311-315 turns the four nodes back into rows). Callers of the query-only add:
  `BridgeIterator.py:359` (blast) and `TangleWorkbench.py:1348`, which reads the returned set's
  TRUTHINESS to drop an image curve on a fresh recompute.
- Readers of a `_Segment`: none needs a `Point`. Coordinates: `bounds`, `intersects`,
  `_find_true_intersection` :454, `_crossing_sign` :366 (+ both cdists), `_resolution_ratio` :403,
  `_cdist_at_point` :1080 (+ cdists). A row: `_outer_neighbours` :422 (chain and position of the row;
  `node_at` for the neighbour, the anchor at -1), the capture :311, the edge key. `manifold` only:
  `_key_of`, the stability asserts, `_remove_segment`'s rtree, tests. Outside: `resolution_probe.py:51-63`
  (coords + the `end()` mapping that becomes the identity), `profile_numerics._tangle_bytes` :367-375
  (ITERATES each `_edge_to_sid` key: an int key raises `TypeError`, so it changes with the keys, not
  with the field names), `test_tangle_index_and_orientation.py:34-37, 50-56, 108-110` (id-frozensets,
  `_Segment(sid, stub, Point, Point)`), `test_segment_index_bookkeeping.py:66-74` (`_RecordingIndex`
  implements `intersection` only: the batched add needs `intersection_v`).
- Measured with the prototype (`_Segment(id, manifold, store, row0, row1, in_rtree)`, one row key
  `((lo+1) << 32) | (hi+1)` for both maps, the query-only add as `_query_curve`: rows from `_rows()`,
  keys as one int64 array, the new-edge filter in a Python loop over the keys (same order, same assert),
  sids `list(islice(Tangle._ids, n))`, bounds from ONE `store.xy` view, ONE `intersection_v`, the scalar
  `_do_segments_intersect` per candidate, per edge the same set comprehension `frozenset((sid, cand))`
  and `|=`, a `_Segment` only for a paired edge; the cut's update from the bridge's row keys):
  - the trace (per query-only add: the ordered (row key, sid) list, the pairs as tuples in frozenset
    order, `new_sids`, `_intersecting_segments`' iteration order after the add; per
    `_resolve_crossing_pair` the unpacked tuple in call order; per cut the keys added; final
    `repr(Tangle._ids)`) is EQUAL on the 18 probe builds and the scale series: 2,536 adds, 4.76 M new
    edges, 8,396 pairs, 10,275 resolve calls. The cut adds NO new key anywhere on the corpus (its update
    re-keys edges the curve's own query already keyed; only a separator's edges would be new).
  - probe: identical (62 dumps, 6024 arrays, 744 texts), traced and untraced.
  - suite: 1244 pass / 6 fail / C3 xfail only; the 6 are the three orientation tests (Point-built
    `_Segment`, id frozensets) and the three bookkeeping tests on `_RecordingIndex`.
  - profiler, two sequential pairs (base `1c25afb` then prototype): blast 73.70 -> 64.88, 73.96 -> 65.14
    s; `index` 8.77 -> 0.57, 8.86 -> 0.56; `cut` 0.61 -> 0.07 (both); `resolve` 0.48 -> 0.40, 0.50 ->
    0.38 (a view per accessor call costs nothing measurable); census 395.4 -> 383.9 MiB, `segment_index`
    122.0 -> 110.4 (edge keys 118.43 -> 106.80; records 0.63 -> 0.68 for the `store` slot); peak RSS over
    baseline 1889 -> 1820 and 1823 -> 1870 MiB (noise); per-blast counts, `pip_depth`,
    `points_created` equal on 32 blasts.
  - micro: `islice` of 1.39 M sids 13 ms, a `next()` loop 44 ms: the sid allocation is not a cost.
  - at scale the 706 `intersection_v` calls take 0.22 s for 1.39 M boxes and return 4,004 candidates
    (2,982 edges have one): candidate work is negligible; a public view per candidate would cost ~3 ms.
  - m3 (each edge's candidates reversed) leaves the trace, the set iteration order and the probe EQUAL;
    m4 (cdists as np.float64, so the crossing's cdists are np.float64) leaves the probe IDENTICAL: the
    probe does not guard the Python-float type.
  - the suite catches NEITHER m1 nor m2 (only the six expected failures): trap (b) is guarded by the
    probe and the trace alone.
  - trap (b) has teeth: m1 (`frozenset((cand, sid))`) moves the trace on 12 of 19 builds and 30 probe
    items (registry rows, <= 1.5e-8 absolute in cdists of ~1e8); m2 (no sid for an unpaired edge) moves
    46 items, among them two image curves (and their partial pieces) DROPPED on recompute
    (`nested_unfiltered.recomputed`, `shared_straddle.recomputed`): a curve whose new edges pair nothing
    returns an empty set and `TangleWorkbench.py:1348` drops it.
- Invariants relied on: one store per Tangle (asserted in the prototype; every suite Tangle has one);
  rows < 2^31 for the int64 key (5.1 M at the 1e-9 plateau); a chain holds no row twice, so a curve's
  keys are distinct and the filter may test the pre-call sets; `intersection_v` returns each box's
  candidates in `intersection`'s order (equal set iteration order on every add).
- Smells: `_insert_segment`'s `else seg.id` (dead); `_Segment.bounds` docstring names `p0` /
  `p0_seg1`; the cut's re-keying adds no key on the corpus (keep it anyway; believed, by reading: an
  image-curve edge that is INDEXED is skipped by the query without entering `_query_edges`, and only
  the cut's update puts it there, which changes what `_insert_segment` :671 decides for it later).

## After N5b `segments-and-edge-keys-on-rows` (auditor, 2026-10-26, branch @ `c944cf6`, base `1c25afb`; lines are the branch's)

Scratch: `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/7a98b0ff-df74-46bf-9769-0728d7d299d5/scratchpad/aud5b/`
(`suite_branch.log`, `basecopy/` = base + the branch's bookkeeping test, `m7` / `m8` / `m9` mutants,
`inst/` instrumented copy, `dump_br` / `dump_m9` / `dump_basesub` probe dumps of 8 builds).

- Map (`numerics/Tangle.py`): `_edge(end0, end1)` :94 (rows; a `BranchPoint` end is -1 through the +1
  offset; NO 31-bit assert); `_edge_keys(r0, r1)` :103 (vector twin, asserts rows < 2**31);
  `_curve_ends(manifold)` :109 -> `(store, r0, r1, head)` from `manifold._rows()` (-1 = the root
  BranchPoint); `_Segment` :129 `(id, manifold, store, row0: int | BranchPoint, row1: int, in_rtree)`,
  `bounds` / `intersects` read `coords()`, `coords(xy=None)` :152 ((2, 2), the caller's `store.xy`
  view or one per call), `cdists(stability)` :159 (Python floats); `Tangle._ids = count(0)` :216;
  `clear_all` :269 resets `_store`; `seg1_id, seg2_id = tuple(seg_id_pair)` :319 (unchanged); capture
  :344-360 (asserts `u_seg.store is s_seg.store`, brackets = the records' rows); `_outer_neighbours`
  :442 (rows -> `store.chain_of` / `position_of`, neighbour still a `Point` via `chain.node_at`);
  `_adopt_store` :564 (first store wins, later curves asserted); `_query_curve` :572-622 (membership
  loop in curve order, `islice` sids, one `xy` view, BranchPoint end masked by
  `a[e0 < 0] = head.get_point()` on a `np.maximum(e0, 0)` gather, one `intersection_v`, scalar
  `_do_segments_intersect` per candidate on `coords(xy)`, `frozenset((sids[k], cand_id))` then `|=`, a
  record only for a pair, returns `set(sids)`); `_insert_segment` (one key; dead `else seg.id` :755);
  `_segments_of` :842 (row records from `_curve_ends`, adopts the store); the cut's re-keying
  :1100-1102 (`_edge_keys(...).tolist()`); `_interpolate` :1179 shared by `_cdist_at_point` and
  `_cdist_between` (the latter keeps node arguments for the separator).
- Invariants: one store per Tangle (asserted at every add; unpinned: m6 survives); rows < 2**31 for
  vector keys; trap (b) -- float bits of `_find_true_intersection` follow `tuple(frozenset(pair))`, so
  the sid values, one sid per NEW query edge (paired or not), and the new-sid-first pair order must
  not change: pinned by `test_a_crossing_pair_keeps_the_new_sid_first` (also pins Python-float cdists)
  and `test_a_curve_crossing_nothing_takes_one_sid_per_new_edge` (a curve whose add returns no sid is
  dropped by `TangleWorkbench.py:1345-1351`).
- Verified by the auditor: suite 1252 / 43 skipped / C3 alone; both pins pass on `1c25afb`; 14 own
  branch dumps identical to the coder's base dumps; profiler JSONs recomputed (blast 74.9 -> 65.2 s,
  `index` 9.0 -> 0.63, `cut` 0.63 -> 0.08, census 395.2 -> 383.9 MiB, only `segment_index` moves).
  Mutants: cut re-keying dropped (m7) fails the separator test; unsorted `_edge` (m8) fails
  `test_readding_a_bridge_is_idempotent`; BranchPoint mask in `_query_curve` dropped (m9) fails
  NOTHING and leaves the dumps identical -- the path is reached once in the suite
  (`test_order_arrays.py::test_a_separator_after_the_periodic_point_becomes_the_first_point`),
  deciding nothing.
- Smells: m9 (the BranchPoint root edge of a query-only curve is unpinned); scalar `_edge` lacks the
  31-bit assert; capture comment :355-357 says "nodes" / "never off the index"; the periodic point has
  four spellings now (object, -1 key offset, clip-and-mask in `_query_curve`, `Chain.anchor`; N8);
  `loom/ResonanceZone.py:27-30` still explains the trim with `segment.p0_seg1` (N6);
  `diagnose` is 82 % of the blast (53 of 65 s at scale).
- Features: the query-only add at scale costs 0.63 s per 32 blasts (was 9.0); verified by the
  profiler pairs and the shadow trace (`n5b/tools/trace.py`, equal on 19 targets).

## Readers on views (N6 analyst, 2026-10-26, main `9dc702b`, measured)

Scratch: `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/7a98b0ff-df74-46bf-9769-0728d7d299d5/scratchpad/n6/`
(`base/` = main archive, `proto/` = the prototype, `patch_proto.py` builds it, `prototype_src.diff`,
`memo_shadow.py`, `reader_shadow.py`, `rss_probe.py`, `prof/` profiler pairs, `v_*` single-change variants).

- Map at main: the bridge walk memo is `BaseManifold._walk_cache` (:125) behind `_memoise_walks`
  (:85 False, `Bridge.py:63` True), served in `get_point_array` (:403-419, key `(final_node,
  return_nodes)`, arrays read-only); `bump_version` (:144-154) clears it; `TangleWorkbench._bump_generation(
  invalidate_walks=True)` (:433-445) bumps EVERY bridge's version (callers `BridgeIterator.py:332,416`,
  `TangleWorkbench.py:720,1810`). The bridge `version` has a second reader, `Trellis.cached_bridge_polyline`
  (`Trellis.py:570-603`), so the bump stays when the memo goes. `Chain.version` (`PointStore.py:85,88,107,
  170,179`) is written, never read.
- Readers that still walk `Point` objects outside the engine (all at main): `numerics/diagnostics._curves`
  :259-284 (`id(node)` sets :276/:281) and `_segment_rows` :383-400 (`id(node)` keys :388, `get_point()`
  per node :387); `geometry.arc_polyline` :79-95; `TangleWorkbench._register_anchors` :1625-1629 and
  `_anchor_crossing_sign` :1643-1665; `trim_stable_at` :2083-2091 (`_walk_nodes`, the LAST src caller);
  `grow_until_turnaround` :773-777 (`first_node` + `get_point()`; the §H `_coords` / rebinding defect is
  ALREADY gone since N3/N4); topology `StablePartition.py:835`, `topology/diagnostics.py:91-99,208-209`,
  `plotting.py:428-432,443,483-505`. Ends (root/tail as `Point`) are read by `Bridge.joined`,
  `BridgeIterator.py:344,348`, `ResonanceZone.py:364-365,521`: those stay nodes until N7's end records.
  Engine-side node reads left: `Tangle.py:1023` (an anchor cut builds the whole node list for two nodes),
  `ManifoldMachine.py:343-345`.
- `points_image` is NOT created by readers: the engine makes every `Point` (`ManifoldMachine.py:370,790`,
  `Tangle.py:1140` separators, `PointStore.node` for a phantom link); readers only read `store.nodes`.
  So N6 cannot shrink `points_image` (101.06 MiB both trees, measured).
- `store.nodes[row]` is one object per row, so `id(node)` keys partition nodes exactly as rows do within
  one store; hand-built test curves live in several stores (`loose_store()` and the workbench's), so a
  row-keyed field must offset rows per store (a one-store assert in `_curves` fails 22 tests), and
  periodic points need DISTINCT keys (they have no row). N5b's `_edge` (-1 for every BranchPoint,
  store-relative) is not a node key.
- Measured (memo shadow, all 18 probe builds + the 1e-7 scale series at 32 blasts): ~120k bridge reads,
  0 where the memo differs from a fresh `_collect` (arrays byte-equal, node lists element-identical);
  scale: 88,958 reads, 29,050 hits, memo path 1.00 s vs fresh 1.12 s. `_curves` reads never hit (the bump
  between blasts clears them). The memo is worth <= 0.12 s and costs 33.48 MiB of census.
- Prototype (`proto/`, ~200 changed src lines in 8 files): memo deleted; `BaseManifold.points` /
  `cdists` (to the TAIL; `get_cdist_array` still reads to the chain's end); `_curves` on `_rows()`;
  `_segment_rows` on `store.xy[rows]` with keys row + per-store offset and -1, -2, ... per periodic point;
  `arc_polyline`, anchors, trim (`np.flatnonzero(store.cdist[rows] >= c)`), turnaround and the six
  topology readers on arrays. Probe identical (62 dumps); per-reader shadow equal at every call on 19
  targets (anchorward 54,833, image-arc midpoints 112,614, counts 97 incl. 23 with hits at scale, ...).
  Suite: only 7 tests fail, all through three hand-built stubs (`_StubManifold` x2, `_NodeBridge`) lacking
  `points` / `cdists`; adding them passes on both trees (64 passed each).
- Profiler pairs (base, proto, base, proto; sequential, scratch trees): blast 62.95 / 63.14 -> 52.46 /
  52.48 s; `diagnose` 51.4 / 51.6 -> 41.9 / 41.8; `walk` 1.21 -> 0.63; topology 1.73 -> 1.17; census
  383.86 -> 350.66 MiB (`walk_memo` 33.48 -> 0.14); per-blast counts / pip_depth / points_created equal.
  Peak RSS over baseline 1831 / 1834 -> 2186 / 2127 MiB, and RSS held after each blast +220..+360 MiB:
  caused by the diagnostics float gather alone (v_memo 1818, v_diag 2126, old node builder over new keys
  1857), and it is glibc's dynamic mmap threshold, not live memory: with `MALLOC_MMAP_THRESHOLD_=131072`
  base 1721.7 vs v_diag 1732.8; with `malloc_trim(0)` after every blast the branch holds LESS (blast 13:
  303 vs 313 MiB; end 460 vs 474). A one-gather field (`v_diag2`) does not avoid it (1940).
- Teeth: none of "every periodic point keyed -1", "no store offset", "no BranchPoint-head check in
  `_curves`" fails the suite tiers or the probe (62 identical each); only a per-reader shadow (key
  equivalence) sees the first two. N5b's m9 (the BranchPoint root edge of `_query_curve`) is pinned by a
  10-line prototype test (`n6/prototype_pin_test.py`, a decoy row 0): passes on main, fails under m9.
- Smells added: `cdists` (to tail) vs `get_cdist_array` (to chain end) differ on a trimmed stable
  branch; `_anchor_crossing_sign`'s docstring names node lists; design note :82-84 and `Trellis.py:8,16`
  still say linked list / walk caches; `ResonanceZone.py:26-30` explains the trim with `p0_seg1`.


## After N6 `readers-on-views` (outcome audit, 2026-10-26, branch @ `e25adf1`, base `9dc702b`)

Scratch: `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/7a98b0ff-df74-46bf-9769-0728d7d299d5/scratchpad/aud6/`.

- Map (branch): the bridge walk memo is GONE (`_memoise_walks`, `_walk_cache`, `Chain.version`).
  `BaseManifold.get_point_array` is one `_collect` line and returns a fresh WRITABLE array per call.
  New public `BaseManifold.points` / `cdists` (root to TAIL, the periodic point first; `cdists` reads
  the head's `.cdist`, NOT `get_cdist(stability)`), both through `_collect` (profiler label `walk`).
  A bridge's `version` is still bumped by `_bump_generation(invalidate_walks=True)`; its ONLY library
  reader is `Trellis.cached_bridge_polyline` (`TangleWorkbench.generation` sums MANIFOLD versions only).
- `numerics/diagnostics.py`: `_curves` returns `(name, (head, store, rows), bridge)`; a bridge is "cut
  from a manifold" iff a manifold of its store holds one of its rows (bool mask per store object) or its
  `BranchPoint` head is a manifold's head. `_hits` builds `offsets` (per store, by first appearance) and
  `periodic` (-1, -2, ... per `BranchPoint`); `_segment_rows(owner, read, offsets, periodic)` gathers
  `store.xy[rows]` once per curve. `diagnose` 52.3 -> 42.4 s self at scale.
- `geometry.arc_polyline` filters `manifold.points` by `manifold.cdists`; `stability=` is accepted and
  NOT read (N7 drops it). `TangleWorkbench`: anchors and `_anchor_crossing_sign(u_points, s_points)` read
  `.points`; the turnaround reads `manifold.points[:2]` (gathers the whole branch to read two rows);
  `trim_stable_at` = `_rows(stop_at_final=False)`, `rows[1:]` when the root is a Point, first
  `store.cdist >= c`, then `store.nodes[row]` as the new tail. Topology readers: `StablePartition.py:835`
  `cdists[-1]`; `topology/diagnostics._anchorward_look` bisects `manifold.cdists`; `_bridge_side_of`
  `np.sort(bridge.cdists)`; `plotting._stable_arc_midpoint` / `_image_arc_midpoint` (argmin).
- Measured (auditor, verified): suite 1255 / 43 / 5 / 1 xfailed (C3); own probe dumps of 8 builds
  (14 dumps) against the coder's base dumps identical (629 arrays, 168 texts); malloc_trim control
  313 -> 302 MiB (blast 13), 474 -> 440 (blast 32); pinned-mmap pair peak 1808.6 -> 1791.5; per-blast
  `same_stability_count` (curve-pair lists included), `pip_depth`, `points_created` equal on six profiler
  runs (23 blasts with hits). Default allocator (coder's pairs): peak over baseline 1838 -> 2056 and
  1871 -> 2027; RSS HELD after blast 13 337 -> 736 and after blast 32 596 -> 957/966 MiB (retained free
  heap, Q7 accepted). Census 383.86 -> 350.65 (`walk_memo` 33.48 -> 0.14, the profiler's empty-dict
  defaults). Walk paths other than `diagnose/walk` rose 0.36 -> 0.63 s (`separation/walk` 0.11 -> 0.33):
  the fresh-read cost of the memo is +0.27 s, not the analyst's +0.12.
- Smells (branch):
  - `arc_polyline` / `BaseManifold.cdists` on a curve rooted at a `BranchPoint` whose `.cdist` is None (a
    FixedPoint not built by `FixedPointSolver`, e.g. tests' `bare_fixed_point`): the base read
    `get_cdist(stability)` (the `(0.0, 0.0)` tuple `FixedPoint.__init__` sets) and worked; the branch
    raises `TypeError` ('>=' None vs float) and `cdists` is an object array (auditor, `aud6/bp_cdist.py`).
    No library path builds such a point. A BranchPoint has TWO cdist fields (`cdist`, `cdists`): N8.
  - `Trellis.cached_bridge_polyline` now caches a WRITABLE shared array (base: read-only, inherited from
    the memo); `topology/diagnostics.py:146` says "must not be written to" with nothing enforcing it.
  - Unpinned (auditor mutants survive the whole suite): `trim_stable_at` with `candidates = rows` (root
    may become tail; unreachable, stable roots are BranchPoints); `arc_polyline` without its reversal
    (the storage-reversed branch is reached 0 of 27,236 calls in the suite); `_curves` without the
    `head in heads` test (= mC / mH).
  - Stale prose outside N6's files: `CLAUDE.md` "bridges memoise their point arrays (read-only)",
    `TangleSession.py:743-745`, `Trellis.py:578-579`, `test_geometry.py:102`, the test names
    `test_bridge_classification_uses_the_memoised_midpoint` and `test_overrun_drops_every_bridge_walk`;
    `profile_numerics.py:64,487` still describe and measure `walk_memo`; design note :82-85 still says
    "(B11, a later optional slice)".
- Node readers left in src (rubric 7 grep): `ManifoldInitializer.py:400`, `ManifoldMachine.py:230,343,345`,
  `Tangle.py:1023`, `ResonanceZone.py:364`, `BaseManifold` own definitions; `store.nodes[...]` 24 sites;
  ends (`Bridge.joined`, `BridgeIterator` wraps, `ResonanceZone` tails) still `Point`s: N7.


## Deleting Point (N7 analyst, 2026-10-26, main `d279815`, measured with two prototypes)

Scratch: `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/7a98b0ff-df74-46bf-9769-0728d7d299d5/scratchpad/n7/`
(`base/` = `git archive d279815`; `proto/` = Point and BasePoint deleted, ends `End` records; `proto_a/` =
engine on rows with `Point` kept as the lazy read adapter; `probe.py` (the probe with an `ids`-aware
`curve_arrays` / `curve_extent`, dumps both trees), `storedump.py` (every store column, byte-compared),
`slotcheck.py` (slot invariant on the base), `touch.py` (AST touchpoint counter), `prof/` profiler runs,
`rss_*.json` malloc_trim control, `proto_src.diff`, `proto_a_src.diff`, `suite_*.txt`).

- Map at main (node readers left in src, file:line): `BaseManifold.py:177` (`_find_tail`, `store.nodes`),
  `:204,252-288` (`first_node`), `:235-250` (`_walk_nodes`), `:350-352` (`_collect(return_nodes=)`), the five
  getters' `return_nodes` (:375-505); `ManifoldMachine.py:230` (`view.root = view.first_node()`), `:295`
  (`view.root.stretch_param` = the chain's), `:343,345` (`return_nodes=True` on the iterate path), `:358-370`
  (rows from nodes, `Point._over` per image: the ONLY bulk Point creation), `:423-492` (`_image_chains` over
  nodes), `:561,584,591` (merge tie on `store.nodes` + `has_iterate_link`), `:788-796` (`Point._over` per
  refined point; slot re-point), `:804-818` (`_get_preiterate`, now only ever called with the anchor);
  `ManifoldInitializer.py:172-199` (the fictitious k-step chain via `insert_prev/next_iterate`), `:233-280`
  (`Point(...)`, links, cdist fill, slot hang), `:400,416,419`; `Tangle.py:454-463` (`_outer_neighbours` via
  `Chain.node_at`), `:1004-1010` (rows -> nodes via `store.node`), `:1022-1024` (anchor cut builds the whole
  node list for two nodes), `:1050` (`previous[2] is p0`), `:1139-1147` (separator `Point(...)`, slot re-point
  after a BranchPoint), `:1175-1176` (`_cdist_between` on nodes); `Bridge.py:149,153,186-204` (`joined`,
  `isinstance(tail, Point)`); `BridgeIterator.py:344,348`; `TangleWorkbench.py:1755` (`get_cdist("unstable")`
  on bridge ends), `:2086` (`store.nodes[row]` as the new tail); `ResonanceZone.py:364-365,516-521`;
  `PointStore.py:112-126` (`Chain.node_at`), `:204,262,286,348-363` (`nodes`, `node()`), `:394-402`
  (`loose_store`, used only by `Point()`). Topology: ZERO code touchpoints (N6 moved its six readers).
- `points_image` (101.06 MiB) is created by exactly two engine lines: `ManifoldMachine.py:370` (images)
  and `:790` (refined points); `store.nodes` is a further 17.2 MiB list.
- Slot invariant (verified, `slotcheck.py`, 74,114 checks on 9 base builds incl. both inversion cases,
  nested_outer_series, p3/k28/k32 series, shared_straddle): a non-empty BranchPoint slot ALWAYS names
  position 0 of a chain whose `anchor` is that BranchPoint. So a slot can hold the Chain itself and the two
  re-point sites (`ManifoldMachine.py:794-796`, `Tangle.py:1144`) become no-ops.
- Prototype `proto/` (src +265 / -626 in 14 files, 443 of the deletions are Point.py + BasePoint.py):
  `End(id, x, y, cdist)` NamedTuple in PointStore + `PointStore.end(row)`; `BaseManifold` stores `_root` /
  `_tail` as rows (or the BranchPoint) plus `_store` (ctor kw `store=`, derived from the slot chain for a
  BranchPoint root), `root` / `tail` build an `End` per read; `ids` (rows root..tail); `_select(kind)` (head,
  store, rows) behind `_collect`; `_slot_chain()` replaces `first_node`; slots hold chains; BranchPoint
  inlines BasePoint (x, y, cdist, edist, links, get_point, set_x/y, get/exists iterate); engine on rows
  (`_manifold_from_rows`, tie rule on the link columns, initializer chain by column writes re-fetched after
  each `append_one` because `_reserve` replaces the arrays); Tangle brackets compared by `==`.
  RESULT: probe identical on 62 dumps; whole PointStore (all six columns + chain ids) byte-identical on
  k10, k28_two_blasts, p3, nested_inner_blasts, inversion_two_blasts, inversion_period2, heteroclinic.
- Profiler (1e-7 scale, sequential, prototype trees; base = `d279815`): census 350.65 -> 232.00 MiB
  (`points_image` 101.06 -> 0, `points_curve` 0.37 -> 0, `point_store` 121.30 -> 104.08 = the nodes list);
  blast 52.02 / 53.39 / 52.16 -> 52.35 / 52.10 / 51.88 s (flat), `iterate` 0.30 -> 0.08, `refine` 0.92 ->
  0.74, `diagnose` unchanged (41.5-42.7); peak RSS over baseline (default allocator) 2053 / 2115 / 1980 ->
  2017 / 2023 / 2078 MiB (noise); `MALLOC_MMAP_THRESHOLD_=131072`: 1760.6 -> 1658.2; malloc_trim(0) after
  every blast: held at blast 13 302 -> 221, at blast 32 440 -> 315 MiB; per-blast count / pip_depth /
  points_created equal in every run. `proto_a` (Point kept as a lazy adapter): census 249.54 (nodes list
  stays, 0.32 MiB of lazily made end Points), same times.
- Which cdist a BranchPoint end is read through (trap a, unchanged by the prototype): `_collect` / `points`
  / `cdists` / `arc_polyline` / `_image_chains` read `.cdist` (solver sets 0.0, hand-built None);
  `_Segment.cdists`, `segment_ends`, `_cdist_between`, `_register_bridge` lo/hi read `get_cdist(s)` (the
  `cdists` tuple, `(0.0, 0.0)` from `FixedPoint.__init__`). `ResonanceZone.restore` read `get_cdist("stable")`
  on Point tails (= `.cdist`); a stable tail is never a BranchPoint on a library build.
- Smells added: `ManifoldView.root` / `.tail` aliases are read by nobody; `refine_manifold(final_node=)`
  unused; `_get_preiterate`'s Point branch is dead (only the anchor reaches it); `loose_store()` exists only
  for `Point()` without a store; `docs/source/tanglepack.numerics.rst` documents `BasePoint` and `Point`
  (and not `PointStore`); `stress/numerics/e02_gpu_end_to_end.py:171` filters profiles by `"Point.py"`.


## After N7a `engine-on-rows` (outcome audit, 2026-10-26, branch @ `88f1a47`, base `d279815`, NOT merged at writing)

- Map (branch): `BaseManifold` stores `_root` / `_tail` as a row of `_store`, the `BranchPoint`, or None
  (`__init__` sets every plain attribute FIRST, then the ends: R1); ctor kw `store=` (required for a row
  root; a `Point` end brings its store; a `BranchPoint` root adopts its slot chain's through the `store`
  property). `_stored(end)` (Point -> row, asserts one store), `_end(end)` (row -> `store.node(row)`, the
  adapter). New: `ids` (fresh int64 rows root..tail, periodic point excluded), `store`, `_slot_chain(b)`
  (the slot's Chain; raises the old ValueError with no branch), `_select(kind, n, final, stop)` -> (head,
  store, rows) behind `_collect`. `_span` / `_rows` / `_find_tail` on rows; `final_node` a row or Point.
  Adapters with no src caller (N7b deletes): `_walk_nodes`, `first_node`, `return_nodes=True`,
  `Chain.node_at` (through `store.node`), `ManifoldMachine._manifold_from_points` (over the new
  `_manifold_from_rows(store, rows, ...)`).
- `BranchPoint.forward_branches[b]` / `backward_branches[b]` hold the branch's `Chain` (or None);
  `insert_point_*(chain_or_point, b)` (a Point stands for its chain until N7b); refinement at position 0
  and a separator after the periodic point no longer touch the slot.
- `PointStore`: `Chain.xy_at(pos)` (anchor's `get_point()` at -1, a read-only `_xy_read` row view, None off
  the ends) is what `Tangle._outer_neighbours` reads. `nodes[row]` is None for every row nobody read via
  `node()`. No `unlink` (R2 deleted it: the initializer's k-step `store.link` is what detaches the last
  fictitious row, `ManifoldInitializer.py:246-248` comment).
- Engine: `_iterate_manifold` reads `_select("non_iterated")` / `_select("iterated")` (no `Point._over`);
  `_image_chains(manifold, head, store, rows)` (slot chain `.index` for a BranchPoint head, `head.cdist`);
  merge tie = `(store._image[pair] >= 0) | (store._preimage[pair] >= 0)`, keep_1 = linked[0] or not
  linked[1]; `_get_preiterate(BranchPoint)` only; `Tangle.create_bridges` keeps brackets as rows, anchor
  cut `(0.0, anchor, head, int(rows[0]))`, separator test with `==`; `_insert_crossing_separator` appends a
  row; `_cdist_between(store, p0, p1, stability, point)` (BranchPoint via `get_cdist`); `trim_stable_at`'s
  new tail an int. `Bridge(store=)`, `joined` / `BridgeIterator` wrap pass `_root` / `_tail` / `store`.
  `ManifoldView` has no `root` / `tail`. `Trellis.cached_bridge_polyline` sets `writeable = False`.
- Who still makes `Point`s (auditor, `aud7a/whomakes.py`): `TangleWorkbench._register_bridge` (597 of 598
  on k28_deep_blast_sliver; it reads `bridge.root` / `.tail` `get_cdist`), `ResonanceZone.define_resonance_zone`
  (1), and test builders reading ends (`cases.build_near_tangency` 102). Scale at 32 blasts: 3,312 of
  2,177,788 rows (0.15 %). These are N7b's end readers.
- Verified (auditor): whole state equal base vs branch on 10 cases beyond the team's 7 (k28_one_blast,
  nested, inversion, p3_6_blasts, heteroclinic_17u (1.16 M rows), k28_deep_blast_sliver, p3_coarse_blasts,
  inversion_period2_10, near_tangency, p3_16_12): six column hashes, chains (ids, stretch, anchor), every
  manifold and bridge end row, version, slot chain index and `workbench.generation`. Suite 1258 / 43 / 5 /
  1 xfailed (C3). malloc_trim control 302 -> 236 (blast 13) and 439 -> 332 MiB (blast 32).
- Measured (team): census 350.51 -> 249.32 MiB (`points_image` 101.06 -> 0.24, `point_store` 121.30 both:
  the `nodes` list stays to N7b, ~17.2 MiB); blast ~52.4 -> ~51.9 s; `iterate` 0.30 -> 0.07, `refine` 0.92
  -> 0.73; the growth's row reads moved from `walk` to `iterate` (`_select` is not an `INNER_STAGES` target).
- Smells (branch): src is +142 lines net (adapters kept + docstrings; BaseManifold 534 -> 641);
  `refine_manifold(branch_index=)` / `_refine_layer(..., branch_index)` unused (`test_refinement.py:55`
  passes it positionally); `Point.py`'s module docstring is false since N7a ("a phantom row gets one only
  when `node` reads it", "every reader ... gets them from `store.nodes`", "read once per point on the hot
  paths") -- the file goes in N7b; `PointStore._alias_xy` docstring names only `Point.get_point` as its
  reader; `_image_chains` with a BranchPoint head and no rows would fail on `store=None` (base failed too,
  unreachable); `profile_numerics.py` has two pre-existing F401s and its census reads `n_points` 0 on a
  store-less tree.


## N7b `delete-point` analyst (2026-10-27, main `4c9c3cb`; prototype `n7b/proto` in the scratchpad)

- Prototype src: +217/-633 in 15 files (443 = Point.py + BasePoint.py) before the docstring sweep (57
  lines naming Point / the node API in 7 src files, 32 of them in `BaseManifold.py`). Probe identical
  (62 dumps), whole store byte-equal (7 cases), end shadow equal (`_register_bridge` 5,819 calls,
  `previous_tails` 22 zones; no probe build calls `restore`).
- Design that works: `End(id, x, y, cdist)` NamedTuple in `PointStore.py`, `PointStore.end(row)`;
  `BaseManifold._stored` accepts End / int / np.integer / BranchPoint (raises TypeError otherwise,
  ValueError with no store after slot adoption); `_end` -> `store.end`; `Bridge.head` alias;
  `BranchPoint` inlines BasePoint verbatim (no base class); `_register_bridge` reads `get_cdist` on a
  BranchPoint else `End.cdist`; `ResonanceZone.restore` the same through a `_stable_cdist` helper
  (needs a runtime `BranchPoint` import in ResonanceZone).
- `arc_polyline(stability=)` has 10 tracked callers, not 6: also `TopologyResults.py:571,597`,
  `topology/plotting.py:3043`, `loom/cartoon.py:85` (three in topology).
- The probe needs no change for N7b (it reads `_rows` when a curve has `ids`, heads via `get_cdist`).
- Census: `point_store` loses exactly the `nodes` list (heteroclinic 1,441,528 -> 1,261,632 B); a zone's
  `End` tails are counted under `topology` (+688 B on nested); needs a `BasePoint` import fallback and
  `getattr(store, "nodes", None)`. `profile_numerics.py` F401s: `Optional` (:102, 0 uses), `Intersection`
  (:407).
- `scripts/refactor_proof/hole-side-and-own-blast/hole_sides.py` is broken on main since B4a
  (`sp._oriented_bridge_polyline`, `sp._bridge_side_of`, `sp.bisect_left` no longer exist).
- `refine_manifold(branch_index=, final_node=)` / `_refine_layer(..., branch_index)`: no caller passes
  them (src `ManifoldMachine.py:395,416`, test `test_order_arrays.py:332`; `test_refinement.py:55`
  passes `None` positionally to `_refine_layer`).


## After N7b `delete-point` (outcome audit, 2026-10-27, branch @ `647786b`, base `4c9c3cb`, NOT merged at writing)

- Map (branch): `Point.py` and `BasePoint.py` are GONE. A curve end reads as `End(id, x, y, cdist)`
  (`PointStore.py:61-71`, a `NamedTuple`, exported from `tanglepack.numerics` only, not `tanglepack`),
  built per read by `PointStore.end(row)` (`:354-357`), or as the periodic point's `BranchPoint`.
  `Bridge.head` (`Bridge.py:168-171`, property at :169) is a read-only alias of `root` (no caller, no test).
  `BaseManifold._stored` (`:215-230`): `End` -> its `id`; `int` / `np.integer` (not `bool`) -> row;
  `BranchPoint` / None kept; anything else `TypeError`; a row with `self.store is None` (after slot
  adoption) `ValueError`. `_end` -> `store.end(row)`. Getters return arrays only (`final_node: End | int`).
  `tanglepack.__all__` 36.
- `BranchPoint` (346 lines) has no base class (`__mro__ == (BranchPoint, object)`): `BasePoint`'s
  public members inlined (x, y, `cdist` scalar, dead `edist` / `edists` / `*_stretch_params`, links,
  `get_point`, `set_x` / `set_y`, `insert_next_iterate` / `insert_prev_iterate`, `get_` / `exists_*_iterate`,
  `has_iterate_link` (no reader), `_set_coords`). Both cdist fields: `.cdist` read by `_collect`,
  `points`, `cdists`, `arc_polyline`, `_image_chains`, `_refine_layer`; `get_cdist(s)` by
  `_Segment.cdists`, `segment_ends`, `_cdist_between`, `_register_bridge` (`TangleWorkbench.py:1756-1759`)
  and `ResonanceZone._stable_cdist` (`ResonanceZone.py:84-86`). `insert_point_*` take a `Chain` only.
- Gone: `PointStore.nodes` / `node`, `Chain.node_at`, `loose_store`, `_walk_nodes`, `first_node`,
  `return_nodes`, `_manifold_from_points`, `arc_polyline(stability=)`, `refine_manifold(branch_index=,
  final_node=)`, `_refine_layer`'s `branch_index`.
- Verified (auditor, `n7b_audit/` in the scratchpad): fresh probe dumps of BOTH trees, compare
  `identical: 62 builds, 6024 arrays, 744 texts`; whole store byte-equal on the 7 cases (own dumps);
  suite on the tip (see testing.md); `End` is a value (`==` True, `is` False); `tail=2.7` / `tail=True`
  raise `TypeError`.
- Measured (team): census 249.32 -> 231.85 MiB (`point_store` 121.30 -> 104.08 = the `nodes` list;
  `points_*` keys absent, not 0); blast flat (+1.2 % / -0.5 %); mmap control 1722.9 -> 1639.1, trim
  control 235.1 -> 224.1 (blast 13) / 332.3 -> 315.4 MiB (blast 32).
- Smells (branch):
  - `geometry.py:34` `from .Intersection import Stability` is now unused (F401, new: its only use was
    the dropped `stability=` annotation).
  - `_stored` checks neither the row's range nor the `End`'s store: an `End` of ANOTHER store, or an int
    past `len(store)`, is stored silently and `root` / `tail` then read uninitialised capacity
    (`PointStore.end` indexes `_xy[row]` up to capacity). The base asserted a `Point`'s store and raised
    `IndexError` on `nodes[row]`. No library path does it (hand-built curves only).
  - A NaN-cdist row end now sorts / compares silently (`End.cdist` NaN; base `Point.cdist` None raised).
    Unreached: bridge and zone ends always carry a cdist.
  - `BranchPoint.py` module + class docstrings still say it is "a fixed point or a crossing" and a
    "doubly linked list"; Dev Notes name `only_forward` (gone since N4) -- N8's prose.
  - `_image_chains`' `-1` head placeholder (the plan cites `ManifoldMachine.py:490` at base; never an index), N8's.
  - `stress/profile/profile_numerics.py:418,420` keeps a `getattr(numerics, "BasePoint", BranchPoint)`
    fallback so the branch census can size older trees.


## N8 `branch-point-storage` analyst (2026-10-27, main `d0d78c5`, read-only; scratch `n8/` in the scratchpad)

- Map at main (`BranchPoint.py`, 346 lines, no base class): built by `FixedPoint.py:126-128`
  (`BranchPoint(_EIGENDIRECTIONS, (0.0, 0.0))`); the solver writes `set_x` / `set_y` (`FixedPointSolver.py:165-166`),
  `.cdist = 0.0` (:168) and the iterate cycle (:179-180). Slots: `forward_branches[b]` / `backward_branches[b]`
  (a `Chain`), read by `BaseManifold._slot_chain` (:302-321) and `ManifoldMachine._image_chains` (:477-478),
  written by `insert_point_*` from `ManifoldInitializer.py:256-258,388-393` (and `polylines.py`, tests).
  `.cdist` read by `BaseManifold._collect` (:373), `_image_chains` (:482), `_refine_layer` (:732);
  `get_cdist` by `Tangle` (:161, :1172), `_register_bridge` (`TangleWorkbench.py:1757`), `ResonanceZone._stable_cdist`
  (:86), `laws.py:975`, the probe. Topology reads NO `BranchPoint` member (grep: 0).
- Who reads what (instrumented `BranchPoint`, `n8/instr`, counts by caller): on the 18 probe builds `.cdist` is read
  only by `_collect` (4,903), `get_cdist` by `Tangle` (946) / `_register_bridge` (205); `next_iterate` /
  `prev_iterate`, `_has_iterate` / `_iterate_of` with a head, `_get_preiterate`: 0 reads. In the SUITE role 4 IS
  read: `_select` with a periodic head and `kind != "all"` from `_iterate_manifold` (4 tests: an anchor bridge
  mapped in a blast, `test_merge_branches::test_an_anchor_bridge_image_is_its_branch_from_the_periodic_point`,
  `test_order_arrays::{test_a_blast_holds_the_invariants, test_growth_after_a_blast_holds_the_invariants,
  test_ids_are_the_rows_from_root_to_tail}`) and `_get_preiterate` from `_refine_layer` (6 calls, the same tests
  plus the periodic-pair refine test). So the anchor's iterate cycle is real semantics (the image of orbit point o
  is orbit point o+1), not dead: it can move to a view, not be deleted.
- Slot bijection (verified on k10, p3, inversion, inversion_period2, heteroclinic, nested): every chain with an
  `anchor` sits in exactly one slot of that anchor, and every such slot roots one registered manifold; inversion
  fills `(2, 2)` slots, others `(1, 1)` in slot 0.
- Store ROW experiment (`n8/rowshift`: `construct_fixed_point` appends one row per orbit point, nothing reads it):
  probe identical (62), but every other row shifts by P (all fixed points are constructed before growth on the
  law cases: rows 0..P-1), so the store is equal only after dropping those rows and renumbering
  (`n8/cmp_remap.py`: 7 cases equal after remap; P = 1,1,3,4,1,2,4). A row buys nothing while N5 Q1 keeps the end
  the object and trap (c) keeps it off every chain.
- View prototype (`n8/proto`, `proto_src.diff`): `BranchPoint(fixed_point, orbit_index)` stores nothing:
  `x`/`y`/`get_point()` read `fixed_point.coordinates[o]` (read-only view), `get_cdist` returns 0.0, slots in
  `FixedPoint._slots[(stability, orbit, branch)]` via `chain()` / `attach()`, `next_iterate` / `prev_iterate` /
  `iterate(n)` derived from the orbit. src +120/-345 in 6 files (BranchPoint 346 -> 133). Probe identical (62),
  whole store BYTE-equal on the 7 cases with no added rows, suite green after 8 test files move (5 ids removed),
  scale profile two pairs: blast flat, census equal once `workbench_other` stops at `Chain` (else +42 KB: the
  slot dict on `FixedPoint` lets the walk reach the branch chains' ids).
- Smells (main): `_image_chains`' `[-1]` head row (:479) is never an index (`rows[i:i+1]` only runs for an
  off-chain image); `Tangle.py:1025` comment names `head.next_iterate`; `Intersection.py:33`, `PointStore.py:43`,
  `BaseManifold.py:19-20` and `design_note.md:83-85` say the storage is N8's; `invariants.assert_one_to_one`'s
  periodic back-link check becomes a tautology once the cycle is derived. `diagnostics._curves`' `head in heads`
  treats a bridge on slot 1 as cut from a slot-0 manifold of the same periodic point (the mH state; no build
  reaches it: every slot holding a chain roots a registered manifold).


## After N8 `branch-point-storage` (outcome audit, 2026-10-27, branch @ `91677b8`, base `d0d78c5`, NOT merged at writing)

The N8 analyst section above describes `d0d78c5`; on the branch it reads as follows (lines are the branch's).
- `BranchPoint` (`BranchPoint.py`, 133 lines, no base class) is a VIEW of its `FixedPoint` and stores only
  `fixed_point` and `orbit_index`: `x` / `y` (`float`) and `get_point()` (a fresh read-only view per call) read
  `fixed_point.coordinates[orbit_index]`; `get_cdist(stability)` is 0.0 (`ValueError` on any other string);
  `chain(stability, branch)` / `attach(chain, stability, branch)` read / write `FixedPoint._slots[(stability,
  orbit_index, branch)]` through `FixedPoint.slot_chain` / `attach_chain` (`FixedPoint.py:526-560`; `IndexError`
  on a branch outside 0..1, `ValueError` on a second chain); `forward_branches` / `backward_branches` are
  read-only tuples; `iterate(n)` = `branch_points[(o + n) % period]`, `next_iterate` / `prev_iterate` its +-1.
  `num_branches = 2` is a class constant, and `FixedPoint._EIGENDIRECTIONS` reads it.
- Gone: the old constructor, `set_x` / `set_y`, `.cdist`, `cdists`, `edist(s)`, the stretch lists, the stored
  iterate links and their seven methods, `has_iterate_link`, `BaseManifold._has_iterate`, the solver's writes
  (`FixedPointSolver.py:161-172` now writes eigendata only), the "Iterate computed incorrectly" error.
- Readers: `BaseManifold._slot_chain` -> `root.chain(...)`; `_iterate_of` -> `node.iterate(+-n)` (stable: back);
  `_select` maps a periodic head for "iterated" and drops it for "non_iterated" (a periodic point always has its
  image); `_collect`, `_image_chains` (`ManifoldMachine.py:477-481`) and `_refine_layer` (`:727-730`) read
  `get_cdist(stability)`; `_get_preiterate` (`:791-801`) = `iterate(-1 if unstable else 1).get_point()`.
- `_stored` (`BaseManifold.py:216-234`) raises `IndexError` for a row outside `0 <= row < len(store)`; a
  cross-store `End` is still undetectable (D4).
- Census (`stress/profile/profile_numerics.py`): no `BasePoint` fallback, no `points_*` branch, no `_point_bytes`;
  `workbench_other`'s walk stops at `Chain` (the slot table would otherwise reach the chains' ids, +42 KB).
  Scale: `branch_points` 880 -> 152 B, `workbench_other` +415 B, total 231.85 MiB unchanged.
- Verified (auditor, scratch `n8_audit/`): fresh probe dumps of both trees `identical: 62 builds, 6024 arrays,
  744 texts`; whole store byte-equal on the 7 cases, no row added; suite 1255 / 43 / 5 / C3 alone.
- Still there (by plan): the six periodic-point spellings (`Tangle._edge` -1, clip-and-mask, `Chain.anchor`,
  `_segment_rows`' negative keys, position -1, `_image_chains`' `[-1]` head placeholder); `diagnostics._curves`'
  `heads` rule (pinned by the mH test).
- Smells (branch):
  - `BranchPoint.attach` / `FixedPoint.attach_chain` validate the branch but not the stability: a misspelt
    stability (`"Unstable"`) writes a phantom slot, anchors the chain and leaves `forward_branches` empty, silently
    (auditor, checked). `get_cdist` does raise on it. Hand-built only; library callers pass a key's stability.
  - `BranchPoint.attach`'s docstring lists only `ValueError` (the `IndexError` on a bad branch is unlisted).
  - `BaseManifold.py:20-21` (Dev Notes) says the head's image is "the next orbit point"; along a stable curve `_iterate_of`
    reads the previous one.
  - `FixedPoint` now holds mutable state (the slot table) though its module docstring calls it "pure data";
    N10 must serialise `_slots`.
  - Stale prose outside N8's files: `IterateInference.py:37` ("`next_iterate` links"), `gpu.py:27` ("Python
    `Point` objects, pointer walks", N9), `test_manifold_initializer.py:101-102`.

## N9 `array-namespace` analyst (2026-10-27, main `a2a2c4b`, read-only; scratch `n9/` in the scratchpad)

- Where time goes at scale (N8's four scale profiles): the five §C N9 hot paths (map 0.035, iterate 0.08, merge
  0.001, image_chains 0.001, refine 0.77, index 0.70, resolve 0.39 s) are 3.5 % of 53-56 s -> Amdahl cap 1.037x;
  `diagnose` 78.5 %, `separation` 9 %. GPU end to end on main (`n9/e2e.py`, 3 interleaved, medians): k28 11 steps
  0.990x, p3 14 steps 0.981x, k10 0.983x; k28 at 12 steps now hits `PointBudgetExceeded` (escape), so no batch
  reaches the 6e5-4e6 break-even. Map body 0.56 % of k28 growth (cProfile); refine (`_refine_layer`, its
  `coords` closure, the geometry kernels) is the rest.
- What is host by nature: merge (`ManifoldMachine.py:519-620`, chain/id bookkeeping), segment build
  (`Tangle._query_curve` :571-621 -> `rtree.intersection_v` :606, CPU libspatialindex), crossing resolution
  (`_resolve_crossing_pair` :298, `_find_true_intersection` :465-506, per-pair 2x2 `np.linalg.solve` :505),
  refine's store reads/writes (`_refine_layer` :704-731, :778-784). Pure array kernels: `geometry.curvature_area_batch`
  (:292-375) and `chord_deviation_batch` (:378-406); callers `ManifoldMachine.py:750/811`, `Tangle.py:435`,
  `diagnostics.py:384-385`, all host. The map seam is `DynamicalSystem._apply_batch` (:129-144) + `gpu._wrap_for_gpu`
  (:104-118); user maps reach CuPy by NEP-18.
- Today the geometry kernels FAIL on CuPy inputs only at `np.asarray` ("Implicit conversion ..."); every other
  call (`np.hypot`, `np.where`, `np.maximum`, `np.isfinite`) already dispatches. `cp.errstate` does not exist
  (`np.errstate` is a harmless no-op for device math).
- Device != numpy bitwise: CuPy `x**2`, `x**3`, `hypot` differ in the last bit on 25-29 % of elements; `+ - * /`
  equal. On real refine inputs (k28_two_blasts 17,132 pairs, p3 3,386) 66 % of areas differ, worst rel 2.8e-7 /
  2.9e-12, 0 flips across `area_cutoff`. Random ill-conditioned inputs reach 1.1e-4 rel. A device refine is a
  tolerance path only.
- Prototype (scratch, removed; `n9/proto_src.diff`, `n9/namespace.py`): `array_namespace(*arrays)` reading
  `sys.modules.get("cupy")` (returns `numpy` itself; raises `TypeError` on mixed host/device) + the two kernels on
  `xp`. Probe identical (62 / 6024 / 744) vs `n8_merge/dump_tip`; whole store byte-equal on the 7 cases; 60/60
  random numpy calls byte-equal to main's kernels; suite 1258 / 43 / 5 / C3; scale blast 54.38/54.53 base vs
  54.47/54.77 proto, census 231.85 MiB equal, per-blast triple equal; lookup overhead 0.5 us/call.
- Stale prose for N9: `gpu.py:8, 26-27, 32-34` (linked list / `Point` objects); `design_note.md:700` numbers.

## After N9 `array-namespace` (outcome audit, 2026-10-27, branch @ `4dd587f`, base `a2a2c4b`, NOT merged at writing)

The N9 analyst section above describes `a2a2c4b`; on the branch (lines are the branch's):
- `numerics/namespace.py` (39 lines, NOT exported, in the Sphinx rst): `array_namespace(*arrays)` returns the
  `numpy` module itself unless every input is a `cupy.ndarray` (then `cupy`); a mix (a list counts as host) raises
  `TypeError`; it reads `sys.modules.get("cupy")` and never imports CuPy (`import tanglepack` leaves `cupy` out of
  `sys.modules`, pinned by `test_tanglepack_never_imports_cupy`).
- `geometry.curvature_area_batch` (`:294-387`) and `chord_deviation_batch` (`:390-426`) are the ONLY namespaced
  functions: `xp = array_namespace(p0, p1, left, right)`, then `xp.asarray / hypot / where / abs / isfinite /
  maximum`; `np.errstate` stays (a no-op for CuPy). CuPy in -> CuPy out, never a transfer. Every library caller
  (`ManifoldMachine.py:750,811`, `Tangle.py:435`, `diagnostics.py:384-385`) passes host arrays, so the library
  always computes them in numpy, byte-equal to the base (auditor: 113 outputs incl. warnings, with and without CuPy
  imported). `chord_deviation_batch` resolves the namespace twice (its own and the inner call): negligible.
- `enable_gpu` / `disable_gpu` / `_wrap_for_gpu` unchanged (AST equal); the map stays the only host<->device
  transfer. `enable_gpu`'s log message still quotes "map <= 0.84% of growth" (left to keep the AST equal).
- `gpu.py`'s docstring "What runs where": device-capable = the user's batched map (NEP-18) and the two kernels;
  host-only with reasons = `merge_manifolds` (`ManifoldMachine.py:519`), `_refine_layer` (`:674`), the segment
  build named as `Tangle._query_curve` (`:571`; the bulk load is `add_manifolds` `:623`, unnamed), crossing
  resolution `_find_true_intersection` (`:465`, `np.linalg.solve` `:505`). Measured result re-stated 2026-10-27
  (0.98-1.00x end to end; map 0.56-0.60 % of k28 growth, 0.06 % of the scale series; hot paths 3.5 %, cap 1.037x;
  `diagnose` 78 %).
- Device parity (auditor, k28_two_blasts stored chains, 9,151 quadruples): 6,041 not bit-equal, max rel 7.9e-11,
  0 flips at 1e-7. Refine layers reach 2.8e-7 (coder `t4.py`).
- Smells: the "break-even out of reach" line rests on total points of a step that exploded at the escape, not on a
  measured batch size; `stress/numerics/harness.py:60` inserts ITS OWN tree's `src` first, so a copied script that
  imports main's harness silently loads main's `tanglepack` (bit the coder's first `e2e` rerun); `diagnose`'s narrow
  phase (`diagnostics._hits` `:353-387`, ~80 % of blast) is the one place a device could matter, owned by no slice.

## N10 `save-load-and-after-profile` analyst (2026-10-27, main `ebc9806`, read-only; scratch `n10/` in the scratchpad)

- **The sid counter is process-global and decides float bits.** `Tangle._ids = count(0)` is a CLASS attribute
  (`Tangle.py:215`; drawn at `:599`, `:652`, `:754`); a pair's operand order is `tuple(frozenset)` (`:318`), so the
  sid VALUES pick `seg_1` / `seg_2` in `_find_true_intersection`. Measured (`n10/twice2.py`): building a law case a
  SECOND time in one process moves registry floats on 6 of 7 cases (k10 1, p3 3, heteroclinic 29 of 116,
  nested_inner_blasts 10, inversion_two_blasts 2, inversion_period2 2 crossings; max rel 8.3e-15 in xy, 6.3e-15 in
  stable cdist), N / iterate table / dynamics equal. WP8 (`docs/open_issues.md:579`) fixed only the ids. A PER-INSTANCE
  counter (`self._ids = count(0)` in `Tangle.__init__`, `clear_all` leaves it; `n10/per_instance_counter.diff`): probe
  identical 62 / 6024 / 744 vs `n9_merge/dump_tip`, whole store byte-equal on the 7 cases, second in-process builds
  byte-equal, suite 1259 pass / 2 fail = the two tests that monkeypatch the CLASS attribute
  (`test_segment_index_bookkeeping.py:207-243`: `monkeypatch.setattr(Tangle, "_ids", ...)` at :210, :218, :240 and `next(Tangle._ids)` at :239, :243).
- **What a workbench holds (every field the prototype had to save for a byte-equal resave):** store columns + capacity
  (`PointStore.py:206-213`); chains `ids`, `stretch_param`, `anchor`, retired ones too (`:97-117`); fixed points'
  eigendata, `partial_jacobians` (None per orbit point), `k_value`, `label`, `_slots` in insertion order
  (`FixedPoint.py:132-161`); every curve (`BaseManifold.py:139-150`: `_root`/`_tail` row or BranchPoint, `_store`
  None or set, `_version`, name, key; `Bridge.py:125-129`: `iterated`, `dropped_at`, first/second ids, a partial may
  hold one); `_manifolds` / `_bridges` / `_partial_bridges` / `_image_curves` / `_half_iterated_bridges` /
  `_bridges_at` (`TangleWorkbench.py:268-299`); the registry `_store` (insertion order), `_next_id`, `cdist_tol`,
  iterate table both dicts, `unresolved_links`, `_generation`, `_orders_version`, the four order lists
  (`IntersectionRegistry.py:142-175`; the caches rebuild lazily); each `Intersection` incl. `seg_ids`, label, VALUE
  TYPES (detected coords are `np.float64`, anchors `float`; numpy 2.2 repr differs) and `unstable_manifold` as a curve
  (identity read by `Tangle.create_bridges` `:988` and `Trellis.is_computed_arc` `topology/Trellis.py:544-545`); the
  index: sid counter, `_seg_lookup` records, `_edge_to_sid`, `_query_edges`, `_intersecting_segments`,
  `_processed_pairs`, `_manifold_segs`, `_seg_manifolds`, `_store` (`Tangle.py:243-254`; rtrees rebuilt from the
  `in_rtree` records); `_growth_stats` (`TangleWorkbench.py:260`, read at `:672-702`; NOT in the plan doc's list);
  `_reported_flags` / `_reported_artifacts` (artifact keys hold a float point), the epochs, `_mutations`,
  `_generation_state`, `_generation`; machine `area_cutoff`, `max_points_per_step`, `escape_radius`,
  `resolution_ratio_threshold`. No dead curve (reachable only from a crossing or the index) on any of 20 saves.
- **Traps measured:** (1) the counter must be restored EXACTLY: restoring "past the largest saved sid" moves 11 of 59
  resumed series steps (query-only edges take sids and keep no record); (2) rebuilding the index fresh (recompute
  order, counter continuing) moves 51 of 59 (835 registry rows, max rel 1.5e-14), no topology text; (3) a
  `frozenset` of two colliding sids iterates in INSERTION order, so a loader must rebuild each pair in its saved
  order (`fs = frozenset((a, b))`, else `(b, a)`): p3, nested, inversion_one/two_blasts, heteroclinic had such pairs;
  (4) save must not consume the counter (law-tier fingerprint guard): keep an int `_next_sid`.
- **Prototype** (`n10/persist_proto.py`, 466 lines, no docstrings, private reaches): one npz (float64 / int64 / int32
  only, `allow_pickle=False` loads it) + one JSON (structure, names, ints). Resave byte-equal on all 11 law cases;
  probe state equal before and after a recompute; grow 1+1 steps + recompute + topology equal on 6 cases
  (inversion_two_blasts raises `PointBudgetExceeded` on both); series resumed from 9 save points (k28 0/3, p3 0/5,
  inversion 2, k32 0/5, nested_outer 0/6) in a fresh process: 59 steps byte-identical to the probe dumps, counts and
  `pip_depth` included; also on main's global counter (load re-seeds the class attribute: process-global side effect).
  Scale 1e-7: save 0.12 s / 76.6 MiB at blast 13, 0.22 s / 122.2 MiB at blast 32 (store 107 MiB of it, query edges
  10.6), load 0.19 / 0.33 s, blasts 14-16 after a load byte-equal. `savez_compressed` ~3x smaller, ~13x slower.

## After N10a `save-load` (outcome audit, 2026-10-28, branch `refactor/save-load` @ `90aba08`, base `ebc9806`, NOT merged at writing)

The N10 analyst section above describes `ebc9806`; on the branch (lines are the branch's):
- **Segment ids are per index.** `Tangle._ids` (class `count`) is GONE; `Tangle.__init__` sets `self._next_sid = 0`
  and `_new_sids(n)` returns a `range` and advances it (`Tangle.py:706-710`; draws at `:603`, `:656`, `:764`).
  One `Tangle` per workbench (`TangleWorkbench.py:266`), `clear_all` keeps the counter, so a workbench's sids are
  what they were in a fresh process; only a SECOND workbench in one process changes (now draws from 0, byte-equal
  to a fresh process; on `ebc9806` 6 of 7 cases moved last bits). Probe identical 62/6024/744.
- **`numerics/persist.py`** (730 lines, not exported, Sphinx stanza): `FORMAT_VERSION = 1`; `save_workbench`,
  `encode_workbench`, `write`, `read`, `restore_workbench`; `Encoder` (fixed points by position, curves and
  crossings numbered once by object identity) / `Decoder`; `typed_floats` / `floats_of` (every float with a type
  tag `f` float / `d` np.float64 / `i` int / `-` None, NaN in the npz). npz = plain `np.savez`, float64 / int64 /
  int32 only; JSON = structure + ints, NO float (auditor: 0 floats in 33 JSONs). Load: `np.load(allow_pickle=False)`,
  `"format" != 1` -> `ValueError` naming both; a format-1 file missing a key -> `KeyError` (no schema check).
  Sets are saved SORTED (pair sets, `_query_edges`, owner sets, `_half_iterated_bridges`, `_reported_*`); each
  segment pair keeps its operand order and is rebuilt by `_pair(a, b)` (frozenset of colliding ints iterates in
  insertion order). Fixed points: `ndarray_fields` records which of `coordinates` / eigenvectors / `jacobians` were
  ndarrays (the solver leaves `coordinates` `(period, 2)`), slots restored only through `attach_chain`. The registry's
  `_unstable_keys` / `_stable_keys` are rebuilt from the crossings; caches, `Intersection.brackets`, FixedPoint memos,
  `_map_batchable` and the rtrees (bulk-loaded from `in_rtree` records) are not saved.
- `TangleWorkbench.save(path)` / classmethod `load(path, map, map_inv, jacobian=None)` (`TangleWorkbench.py:2049-2097`)
  forward to persist; `load` builds `cls(map, ...)` then `restore_workbench`. The map is NEVER checked against the
  save: a wrong map loads silently and only the future is wrong.
- `FixedPoint.attach_chain` raises `ValueError` on a stability other than "unstable"/"stable" (`FixedPoint.py:556-557`),
  pinned with the `IndexError` in `test_persist.py`.
- Measured: fresh-process round trip byte-equal on 11 law cases (auditor rerun); 59 resumed series steps 0 differ
  (auditor rerun of 4 save points: 16 steps, 32 comparisons, 0 differ); scale 1e-7 save at blast 13 -> blasts 14-16
  byte-equal; save 0.14 / 0.23 s, load 0.22 / 0.36 s, npz 77 / 123 MiB (blast 13 / 32); census and blast unchanged.
- Smells: `FixedPoint.py:81` documents `coordinates (List[np.ndarray])` (the solver sets an ndarray); a save writes the
  npz then the JSON, non-atomically; `save("x.npz")` writes `x.npz.npz`; zones' shared `intersection_ids` list is
  un-aliased by a load (unobservable). Any change to what a save holds must bump `FORMAT_VERSION` after the merge.

## N10b `after-profile` analyst (2026-10-28, main `08b1d12`; measured in a temporary worktree, removed)

- **The profiler imports the tree it sits in**: `profile_numerics.py:112` `REPO = TANGLEPACK_REPO or parents[2]`, then
  `src` first on `sys.path` (`:113`). Verified: a worktree run's tracemalloc `top_files` read `src/...` relative to the
  worktree. `git_sha()` (`:701-708`) marks `-dirty` only for `src/`; `tests/cases.py` is imported too, so "clean" must
  be checked with `git status --porcelain` by hand. The parent's end-of-run `profile_figures.draw(out, args.compare)`
  (`:769`) draws without a baseline unless `--compare` is passed; `--plot-only --compare` then overwrites the pngs.
- **After-profile at `08b1d12`** (scratch `n10b/prof_1e-7`, `prof_1e-9`, `prof_traced`): 1e-7 census 1914.0 -> 231.9
  MiB (922 -> 112 B/pt), blast 458.4 -> 52.7 s, peak RSS 3034.7 -> 2012.1 MiB (ebc9806 read 1931.9: the peak is
  bimodal), wall 562 -> 62 s, 32 blasts, plateau at 30 on both; 1e-9 census 4683.2 -> 502.3 (963 -> 103 B/pt), blast
  541.4 -> 68.1, peak 6632.6 -> 1932.8, wall 780 -> 83, 29 blasts; traced 1e-7 series (40 blasts, `--modes memory
  --scale-trace-blasts 40`, 67 s) 1880.6 -> 232.3 MiB, census/traced 1.018 -> 0.998. Same-stability counts,
  `pip_depth`, `points_created` and the top-level stage sequence equal to `ebc9806`'s on every blast at both cutoffs.
  What is left at the 1e-7 plateau: segment index 110.4 MiB (48 %), point store 104.1 (45 %), curve ids 10.8; blast
  time is `diagnose` 42.6 of 52.7 s (81 %) and `separation` 4.8.
- **Rows vs objects**: every case reads +2 points per fixed point (k10 3190 -> 3192, nested 17143 -> 17147): the
  census's `n_phantom` reads `no_cdist_rows`, which is `phantom_rows` + 2 per fixed point (rows with no cdist but a
  link; believed the initializer's fictitious tails, which the base held as no object).
- **`profile_figures.py`** (as found at `08b1d12`; FIXED on `refactor/after-profile`, see the outcome section below):
  `figure_stage_peaks(runs, out)` (`:195-216`) took no base and on a scale-only folder drew
  an empty axes with a "No artists with labels" warning; the right panel of `scaling_bytes_per_point.png` (`:245-247`
  `if ls == "-"`) and the `stage_time.png` titles (`:119-120`, this run's wall only) never show the base. Every
  top-level stage sequence matches the baseline's call for call (k10 10, p3 14, nested 24, inversion 10,
  heteroclinic 24, scale 75 / 27 memory). Without `--compare`, the six pngs are byte-reproducible run to run.
- Probe dump of `08b1d12`: 27 s, 129 MB, two dumps identical (62 builds, 6024 arrays, 744 texts). Suite at
  `08b1d12`: 1380 passed / 43 skipped / 5 deselected / 1 xfailed (C3), 103 s.

## N10b `after-profile` outcome (auditor, 2026-10-28; branch `refactor/after-profile` @ `9aa741e`, base `08b1d12`, not merged at writing)

- **What the branch holds.** Three commits: `dc1bd0c` (`profile_figures.py` +26 -4), `361d61e` (the three folders),
  `9aa741e` (README + plan doc row 10b and a "Closed 2026-10-28" line). `src`, `tests`, `profile_numerics.py` and
  `identity_probe.py` have the same tree / blob hashes on all four commits; the `a8bb8a0*` folders have no diff.
- **`figure_stage_peaks(runs, base, out) -> Optional[Path]`** (`profile_figures.py:195-235`): returns None (no file)
  when the folder has no law case; draws each law case's baseline peak dashed (same colour, width 1.5) only when its
  top-level `stage` list equals the baseline's call for call, else prints one "differ" note and draws the run alone;
  the title gains "; baseline peak dashed" iff a dashed line was drawn. `draw` (`:336`) appends it only when not None.
  Without `--compare` every png is byte-equal to the old file's, except that a scale-only folder no longer gets an
  empty `stage_peak_rss.png` (verified by drawing both files on copies of the json). Still not beside the base, by
  decision (Q5): the right panel of `scaling_bytes_per_point.png` and `stage_time.png`'s wall in the title.
- **Committed after-profile** (`figures/numerics_profile/dc1bd0c/` 12 json + 6 png 1.3 MB, `dc1bd0c_scale1e-9/` 1 json
  + 5 png 508 KB, `dc1bd0c_scale_traced/` 1 json 264 KB, every json `git_sha` `dc1bd0c`). Read with the auditor's own
  reader: 1e-7 census 1914.05 -> 231.85 MiB, blast (sum of `blast_s`) 458.41 -> 52.16 s (by_stage total 52.15),
  `diagnose` 398.59 -> 42.01, `separation` 5.09 -> 4.83 (now the second-largest blast part, 9 %), peak RSS over the
  baseline 3034.7 -> 2015.0 (both at blast 13), wall 562.1 -> 60.6; 1e-9 4683.22 -> 502.31, 541.43 -> 67.69 s,
  `diagnose` 443.98 -> 52.56, `separation` 8.82 -> 8.29, peak 6632.6 -> 1930.8; traced growth 1880.6 -> 232.3
  (traced peak 2055.0 -> 410.9). Segment index 110.4 MiB (48 %) and point store 104.1 (45 %) of the 1e-7 census.
  `pip_depth`, `points_created` and the top-level stage list equal to `a8bb8a0`'s on all 14 json pairs. A fresh
  `--plot-only --compare` of copies of the committed json gives 11 pngs `cmp`-equal to the committed ones.
- **Features:** `profile_numerics.py --plot-only --out X --compare B` draws all six figures beside B, the stage peaks
  included (verified: fresh draw equal to the committed pngs; a popped baseline event drops that case's dashed line,
  test engineer and reviewer).

## D1 `diagnose-on-new-segments` analyst (2026-10-28, main `88c62ea`; scratch `d1/` in the scratchpad, worktree removed)

- **Where the 42 s goes** (1e-7 scale series, 30 diagnosed blasts, in-process split `d1/lab_split.py`, sum 41.06 s):
  pair generation of the new-holding pairs 23.5 (the 1-D sweep's forward ranges, ~1.3 G candidates, ~16x the
  2-D box overlaps: 198 M vs 12.65 M at blast 13), narrow phase 8.8 (92.5 M pairs after the shared-node drop,
  0.90 M hits), the two global stable argsorts 4.9, `np.unique(axis=0)` dedup 3.3, rows + concat 0.75, curves +
  box 0.14. The unresolved-detection half of `report_resolution` is 0.001 s (registry caches). So ~78 % is already
  the new segments' own pairs; the whole-field part is ~9 s. The sweep axis is always y (overlap ratio 0.69-0.72).
- **Offline replay** (`d1/lab_record.py` records each blast's curve order, rows, new flags, box and base count;
  `d1/offline.py` replays `_hits` on a fake store): reproduces the base count on every blast (1e-7: 30, 1e-9: 27).
- **The rtree is no option**: the workbench's unstable tree holds manifold segments only (blast image curves are
  query-only, `BridgeIterator.py:360`, `TangleWorkbench.py:1350`), so it misses every earlier image; inserting them
  would change the crossing index (trap b). A separate tree: bulk load 1.4-2.2 s per blast, or 14 us per insert
  (~20 s over the series); a fixed `intersection_v` wrapper answers 14.3 M results in 2.2 s at blast 13.
- **Exact faster enumerations measured** (count equal and hits bitwise equal as multisets on every blast):
  uniform 2-D grid in numpy 54.6 s (fewer candidates, dearer each); the base sweep on arrays permuted into sweep
  order with the narrow phase per block 31.2 s; **the sweep inside strips of the other axis (`d1/v2.py`) 21.6 s**
  (1e-9: 51.7 -> 28.0); v2 plus a 2-D cell prefilter 21.2 / 24.3 (not worth its code).
- **v2's rules**: int64 edge-key `np.unique` (same first occurrences as `axis=0`); the axis from sorted arrays
  (sum of overlaps = sum(searchsorted(sort(mins), sort(maxs), right)) - N(N+1)/2, exactly the base's `counts`
  sum); one stable argsort on that axis; each box in every strip of width 4 x median new other-axis extent it
  spans (at most 4096), only strips a new box spans; within a strip the base's forward ranges (new: all, old:
  the new ones); a pair kept only in the strip of max(other-axis minima) (monotone `floor`, so exactly once);
  orientation = sweep order, so `curve_pairs` and points are the base's; `new=None` lexsorts hits by (first,
  second) sweep position = the base's listing order. Narrow phase per `_PAIR_BLOCK` block.
- **Verified for v2**: probe 62 dumps identical to main's (`n10_merge/dump_tip`) incl. `same_stability_counts`; an
  orientation mutation moves 3 probe items (k32_series), so the probe sees it; `same_stability_crossings` equal in
  order and every float on 88 corpus listings (11 law cases x 2 regions + 4 series every step; 305,614 hits) and
  on the scale field at blasts 13 / 30 (415,572 / 900,464 hits; base 12.4 / 36.3 s, v2 4.6 / 11.1 s); 3,000 random
  polyline sets (ties, shared edges, a shared root, long segments, repeated points) 0 mismatches; suite 1396
  passed / 43 skipped / 1 xfail (C3), perf 5 passed. Profiler (default allocator, two sequential pairs): `diagnose`
  41.8 / 42.2 -> 22.5 / 22.4 s, blast 51.9 / 52.3 -> 32.6 / 32.5, step peak RSS over baseline 2000 / 1901 MiB
  (blast 13) -> 1207 / 1194 (blast 21), per-blast (count, `pip_depth`, `points_created`) equal on all 32 blasts.
  Traced peak of the count at blast 13: 1538 -> 638 MiB.
- **Tests that see what**: the oracle test (`test_agrees_with_the_test_tier_oracle_on_random_polylines`) catches lost
  pairs (wrong reference strip, dropped entries); NOTHING in the suite sees an orientation swap or a changed
  listing order (both pass all 34 diagnostics tests); the probe sees orientation, only a scratch harness the order.
- **What stays per blast in v2** (1e-7 / 1e-9 series sums): pairs 11.4 / 8.1 s, fixed part 10.2 / 20.3 s (rows,
  dedup, axis, and the strip selection, which keeps most of the field because a blast's images spread along the
  folds). Cost in proportion to the blast needs state kept across blasts (sorted keys, axis sums, strip arrays),
  with an invalidation rule (the scale series' curve list only grows: 0 curves gone or changed on 30 blasts).

## D1 outcome (auditor, 2026-10-28; branch `refactor/diagnose-on-new-segments` @ `7d9d009`, base `88c62ea`, not merged at writing)

- **What the branch holds.** Five commits: `1f688c6` (oracle test), `6860f99` (the strip sweep, `diagnostics.py`
  481 -> 537 lines), `1e22b63` (`figures/numerics_profile/6860f99{,_scale1e-9}/`), `c0ca13d` (README paragraph),
  `7d9d009` (the test engineer's tie test). `src` diff = `numerics/diagnostics.py` only; `Tangle.py`,
  `TangleWorkbench.py`, `loom/` untouched; `git rev-parse 6860f99:src` = tip's.
- **`_hits` (`:333-401`) now, in order:** rows (unchanged) -> edge dedup `np.unique(_edge_keys(node_a, node_b),
  return_index=True)` (`_edge_keys` `:436`: `((min+2**20) << 32) | (max+2**20)`, asserted in range; same first
  occurrences as the old `axis=0` unique) -> `_sweep_axis` (`:445`: overlap pairs per axis =
  `sum(searchsorted(sort(mins), sort(maxs), right)) - n(n+1)/2`, axis 0 on a tie) -> ONE stable argsort along it;
  arrays permuted into sweep positions, `order` maps a position back to its row -> narrow phase once per block from
  `_candidate_pairs` -> with `new=None` a `lexsort((j, i))` restores the listing order (by first, then second sweep
  position) -> curvature areas through `order`, owners remapped at the end.
- **`_candidate_pairs(mins, maxs, is_new, axis)` (`:461`), a generator of `(p, q)` blocks, `p < q`:** strips across
  the sweep axis, `_STRIP_WIDTH` (4) x the median across-extent of the NEW boxes wide, at most `_MAX_STRIPS` (4096);
  only strips some new box spans; one entry per (box, spanned strip it spans); inside a strip the base's forward
  ranges (new box: all later entries it reaches; old box: only the new ones); a pair kept only in strip
  `max(first_strip[p], first_strip[q])` = strip of the larger across-minimum (monotone floor), so exactly once.
  The strip width changes cost, never pairs (W1-W5 extremes pass, test engineer). Stateless: no module mutable,
  no workbench attribute.
- **Contracts now pinned by tests** (`test_same_stability_diagnostic.py`): the full listing's ORDER and pair
  ORIENTATION (earlier segment along the sweep axis first, not curve order: `same_stability_crossings` Returns was
  false before D1 and now says "sweep order"), `count.curve_pairs` orientation, the axis-0 tie
  (`test_on_a_tie_the_sweep_runs_along_x`). Any later index (persistent cache, rtree) must reproduce them or change
  these tests on purpose.
- **Measured (committed folders, auditor's reader):** 1e-7 `diagnose` 42.01 -> 22.43 s, blast 52.16 -> 32.55 s,
  peak RSS over baseline 2015.0 (blast 13) -> 1054.6 MiB (blast 21); 1e-9 `diagnose` 52.57 -> 30.00 s, blast
  67.69 -> 45.25 s, peak 1930.8 (blast 12) -> 2008.1 MiB (blast 19): a REAL +60-95 MiB rise at 1e-9 (test engineer,
  five runs, pinned mmap too, and tracemalloc of the call at blasts 13-20), from whole-field permuted copies and
  strip entries when a late blast adds few segments. Census and every per-blast count / `pip_depth` /
  `points_created` equal. Auditor's own pair at 1e-7 (88c62ea vs 7d9d009, load 3.0-3.2): `diagnose` 43.62 ->
  22.92 s (0.525x), blast 54.15 -> 33.43 s (0.617x), triple equal on 33 steps.
- **What is left in `diagnose`** (coder's in-process split, 1e-7, 30 blasts, tip): rows 0.91, dedup 0.52, axis +
  sweep sort + permutation 5.14, strip choice 3.60, pairs 4.46, narrow 7.88 s; i.e. ~10 s whole-field per series
  (~20 s at 1e-9) that only a persistent cache across blasts removes (rejected, Q2). `separation` (4.8 s at 1e-7)
  is the next blast cost after `diagnose`.
- **Smells left:** `count_same_stability_crossings` still uses `np.unique(..., axis=0)` for distinct owner pairs
  (`:241`; tiny input, outside D1's scope); M5 (the shared-node drop) is an equivalent mutant on every reachable
  input (test engineer's argument: a store row is on one chain), kept as base code; the `_candidate_pairs`
  docstring's "about a thousand candidates for some 75 overlaps" is blast 13 at 1e-7 (1,231 / 78.7 per new segment;
  the series runs 176 / 13.7 to 2,132 / 104).
- **Features:** `report_resolution(new_bridges=)` / `blast_zone(diagnose=True)` same-stability count (verified: probe
  identical 62 builds 6024 arrays 744 texts, re-run by the auditor; 88 + 4 listings and 3,000 random sets equal in
  order; warning text byte-identical, test engineer).
