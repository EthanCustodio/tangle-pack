# numerics — knowledge file

Created 2026-10-08 (open-issues round WP2, branch `oi/wp02-exact-iterates`, base `f39f617`). Only
the crossing / iterate-table path is described; growth, refinement and budgets are in `CLAUDE.md`.

## Map (crossings and the iterate table)

- `numerics/Tangle.py` — rtree of segments. `resolve_crossings` :247 iterates a `set` of segment-id
  pairs (:265), `_resolve_crossing_pair` builds each `Intersection` with both keys
  (`_key_of` asserts a key on EVERY curve since WP2), `unstable_segment` and `stable_segment`
  (the two nodes either side of the crossing on each curve), `crossing_sign`, `resolution_ratio`.
- `numerics/Intersection.py` — the crossing record; `stable_segment` added in WP2.
- `numerics/IntersectionRegistry.py` — the store. `DEDUP_RTOL = 1e-9` :35 and `_find_collision`
  :774 (keys equal AND both cdists equal to 1e-9 relative; an anchor at exactly `(0, 0)` matches
  only another exact `(0, 0)`); `register_iterate` :546 asserts the partial bijection;
  `reindex_from` :256 re-gives re-detected crossings their old ids (rematch by `_find_collision`).
  `cdist_tol` (absolute 1e-6) is no longer the dedup rule; it is the topology's slack.
- `numerics/IterateInference.py` — the EXACT rule (module docstring is the statement):
  `infer_iterates(sources=None)` :107 → anchors by `_register_anchor_iterate` :161 (advanced keys
  at `(0, 0)`, `registry.find`), everything else by `_bracketed_image` :220 (carry `U(x)` by the
  unstable curve's `stretch_param`, `S(x)` by `1/stretch_s` of the stable manifold; candidates =
  a bisected window of the advanced stable branch widened by its widest stable bracket
  `_widest_stable_brackets` :207; unique match on both brackets with `BRACKET_PAD = 1e-12`), then
  `_map_gated` :258 (ONE `map_batch` call over `[x, p0, p1]` per source; keep iff
  `|f(x) - y| <= |f(p1) - f(p0)| + |q1 - q0|`).
- Callers: `TangleWorkbench.compute_intersections` :1380, `BridgeIterator.iterate_bridge` :394
  (sources = the parent's endpoints + the crossings just born, searched against the whole
  registry), `iterate_all_bridges` :437 and `loom/Blast.blast_zone` :498 (registry-wide, after the
  loop, whatever the stop reason; not after a raised `PointBudgetExceeded`).
  `TangleWorkbench.infer_iterate_table()` is a forwarder kept for ~15 scripts.

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
  (law `iterate_table_complete`; regression `test_iterates_after_blasts` on p3 at ten blasts).
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
  registered first; blast registration was already ordered (`_collapse_noise_crossings`).
