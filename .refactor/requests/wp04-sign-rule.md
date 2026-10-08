# Request: WP4 — Signs decide hole openings and rows

*From the author, 2026-10-07. Input for one `tangle_refactor` run (copy to `artifacts/request_origin_analyst.md`).*

## The question

**Can hole openings and partition rows be decided by crossing signs alone (D6)?** About 450 lines of polyline geometry (`_hole_openings`, `_row_at`, `_stable_frame`, `_bridge_end_geometry`, `_arc_side_of`, nearest-vertex tests) decide what `row_of_end` and `outward ⇔ (sign>0) == (side==left)` decide. The two readings disagree on coarse polylines, and that disagreement is every stress-campaign assertion failure in NUM-1 and NUM-2.

## Issues in scope

INT-2, INT-3, INT-4, INT-5, NUM-1, NUM-2, PERF-4 — see `docs/open_issues.md` for the evidence behind each.

## What "done" means

- Openings and rows come from signs; the geometric reading becomes an optional check that flags under-resolution (it never asserts when the signs alternate, and still asserts when they are equal).
- INT-4 and INT-5 disappear with the geometry, or are fixed if any of it stays.
- `_containing_bridge` becomes an exact bisect over the ray's crossing order (PERF-4). Say whether the hole should sit on the whole containing bridge or the backward-image sub-arc (INT-2), with evidence; change it only if the chalkboard demands it.
- The E11 cells that assert today (p3 / nested_inner at 1e-5 from blast 6, k28 1e-5/1e-6 blasts 16-17) run through.

## Where to start

`topology/StablePartition.py`, `topology/TopologyResults.py`, `topology/PartitionFamily.py`, `stress/numerics/e11_blast_depth.py` (to reproduce)

## Binding context

- **Catalog:** `docs/open_issues.md` holds the evidence (path:line, probe numbers) for every issue ID named here. `docs/pipeline_chalkboard.md` is the one-board pipeline the code is judged against. Read both first.
- **Author decisions of 2026-10-07** (D1-D8, C1-C6 in `docs/open_issues.md` and `.refactor/decisions.md`) are binding. In short: topology is recomputed from scratch after any growth or blast, and the class count never decreases (D1); the alphabet comes from the minimal trellis only (D2); the shortest dual-graph walk IS the word and must always exist (D3); an inversion anchor is one point with two stable and two unstable branches, and cross-side is flagged only within one branch (D4); keep the period-k iterate rule but close its gap (D5); signs decide holes and rows (D6); heteroclinic pairs give holes in B's partition (D7); simplify freely (D8).
- `CLAUDE.md` invariants are physics. Topology never calls the map. `numerics` and `topology` never import `loom`.
- No golden pins in tests (policy 2026-10-05): new laws are rule-based; an expected failure is a strict xfail in `tests/cases.py` `KNOWN_ISSUES`, and `KNOWN_ISSUES` may only shrink.
- If the work is bigger than one targeted run, plan the first slice and name the rest as later runs under the same WP number.
