# Request: B1 — `resolved-crossings`

*From the author, 2026-10-12 ("I want to just call tangle_refactor"). Input for one `tangle_refactor` run. The second slice of `docs/alignment_plan_2026-10-09.md`; B0 (`oracle-laws`) is merged (main `1029fcc`, run archived in `.refactor/runs/2026-10-12-oracle-laws/`).*

## The question

**Which detections are crossings of the trellis?** Today the blast path votes on near-coincident detections (`_collapse_noise_crossings`: an odd run collapses to its median, an even run is dropped) and no other path does, so a plain recompute after deep blasting changes `N` (k28 at 8 blasts: 248 → 252). Deleting the vote alone turns k28 x8 from 24 holes / 15 classes into 40 / 23 with diverging words, both "reliable". The pipeline has no rule for a detection the polylines cannot resolve once it is registered. Decision Q13 (a) states the rule (plan A.2): a detection is in `N` iff the polylines resolve it; an unresolved one is registered, reported with the remedy "refine", and excluded on every path, the way a crossing past the pip cut is.

## What "done" means

Exactly slice B1 of the alignment plan (`docs/alignment_plan_2026-10-09.md`, section B, "B1 — `resolved-crossings`"): its scope, eight steps (0-7), behaviour, not-doing list and seven rubric items are the plan. In short:

0. **Probe** (test-side script under `tests/_tools/`): for every law case and series step, the same-curve detection pairs closer along the unstable polyline than the chord deviation there, and the minimum per-crossing ratio; its output is attached to the report (rubric 4). The expected answer on the law cases is none; any hit is a law-case behaviour change the report names.
1. `Intersection.resolution_ratio` stays; add `Intersection.run_separation` (distance along the unstable polyline to the nearest other detection of the same `(u_key, s_key)` on the same or an adjacent unstable segment, None when alone), set by `Tangle.resolve_crossings` and `populate_intersections_for_manifold` alike through one helper.
2. `IntersectionRegistry.resolved_ids(threshold)` = ids that are synthetic or (`resolution_ratio ≥ threshold` and `run_separation` None or `> max(h_u, h_s)`); `unresolved_ids(threshold)` the complement. `N` is `resolved_ids` wherever the registry is read as the crossing set: `create_bridges`, `infer_iterates`, `Trellis.from_workbench` (branch orderings, `own_intersection_ids`). Four call sites, one name.
3. Delete `_collapse_noise_crossings` and `_purge_crossing`; `populate_intersections_for_manifold` returns every detection.
4. `report_resolution`: the unresolved set is ONE WARNING per call naming the detections (ids, coordinates, ratio, separation) and the remedy "refine: lower `area_cutoff` or grow the stretch; NUM-3/4". It stays a flag; the exclusion is the definition.
5. Delete `tests/unit/numerics/test_noise_crossing_collapse.py`; `test_resolved_crossings.py` builds a hand-made zig-zag across a stable segment (three detections within the chord deviation) and asserts all three registered, none resolved, one bridge across them, no iterate link to them, the WARNING once; a clean transversal crossing with ratio ≥ 3 is resolved.
6. `cases.blast_series` records `recompute_crossings` per step; `test_recompute_keeps_the_crossing_set(name)` over `cases.SERIES` asserts it equals the step's `len(N)` at every step.
7. `resolution_ratio_threshold` keeps 3 and gains a docstring: the derived bound is 1, 3 is the margin, part of the definition of `N` (A.2).

Files (eight): `src/tanglepack/numerics/Tangle.py`, `src/tanglepack/numerics/Intersection.py`, `src/tanglepack/numerics/TangleWorkbench.py`, `src/tanglepack/numerics/IterateInference.py`, `src/tanglepack/topology/Trellis.py`, `tests/unit/numerics/test_noise_crossing_collapse.py` (deleted) → `tests/unit/numerics/test_resolved_crossings.py`, `tests/cases.py`, `tests/regression/test_class_count_under_blasting.py`; plus the step-0 probe under `tests/_tools/`.

**Behaviour.** Law cases: unchanged (clause 1 holds nothing out; the collapse fired 0 times there; clause 2 is what the probe measures). Series p3, inversion, k32, nested_outer: recompute already exact. **k28 at 8 blasts is the one case the rule decides**: `unresolved_ids(3)` = `{1, 2, 250, 251}` by clause 2 and `{148, 249}` by clause 1. Its holes / seed / closed classes are RE-DERIVED and printed beside the base's 24 / 15 / 23 and the uncollapsed 40 / 23 / 31, with the holes that depended on 148 / 249 named. The claim is "the topology of the resolved crossings", whatever the numbers; nothing is pinned. The corpus minimum ratio (heteroclinic, 3.04) is recorded in the rubric.

**Not doing.** No local refinement (option (b), NUM-3/4); no change to `DEDUP_RTOL`; no change to the same-stability flags.

## Where to start

`docs/alignment_plan_2026-10-09.md` (B1; A.2 for the definition of `N`; E Q13), `docs/code_fidelity_audit_2026-10-09.md` (the noise-collapse finding), `src/tanglepack/numerics/Tangle.py` (`_collapse_noise_crossings`, `_purge_crossing`, `resolve_crossings`), `src/tanglepack/numerics/IntersectionRegistry.py` (`unresolved_crossings`), `tests/cases.py` (`blast_series`, `SERIES`), `tests/unit/numerics/test_noise_crossing_collapse.py`.

## Binding context

- The author's decisions on the alignment plan are in `.refactor/decisions.md` and plan section E; Q13 is decided: (a), margin 3, corpus minimum recorded.
- `CLAUDE.md` invariants are physics. Topology never calls the map; the iterate table is exact; flags never discard — this slice makes the resolution criterion part of the DEFINITION of `N`, which is the one authorised change to that rule.
- Test policy (`tests/conftest.py` docstring, 2026-10-05): rule-based laws, no case-specific pins; k28 x8's numbers are re-derived and reported, never asserted.
- The analyst stage is NOT skipped: the analyst runs step 0 (the probe) and verifies the B1 section's claims against the code at `HEAD` (the four readers of the registry, the collapse's call sites, the 248 → 252 recompute drift on k28 x8), reporting what differs. It does not re-plan the slice.
- The working tree may carry uncommitted documentation changes (`docs/`, `.refactor/`); the change branches from `HEAD` and never touches them. B0's branch `refactor/oracle-laws` is already in `HEAD`.
