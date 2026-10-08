# Request: WP10 — Numerics and performance leftovers

*From the author, 2026-10-07. Input for one `tangle_refactor` run (copy to `artifacts/request_origin_analyst.md`).*

## The question

**What stops deep blasting and growth, now that topology is sound?** The blast-end same-stability pass costs up to 1.96× the blast and cut reachable depth by 3-4 blasts (PERF-1); `merge_manifolds` crashes on non-contiguous iterated points (NUM-9); blast failures are counted as "skipped" (NUM-10); face nesting by point-in-polygon raises `OuterFaceError` on slivers (INT-10).

## Issues in scope

PERF-1, NUM-9, NUM-10, NUM-11, NUM-12, NUM-13, NUM-14, NUM-15, NUM-16, INT-10, TEST-7 — see `docs/open_issues.md` for the evidence behind each.

## What "done" means

- Same-stability pass checks only new image bridges and counts hits (C5).
- NUM-9, NUM-10, NUM-12, NUM-13, NUM-14 fixed; NUM-16 verified gone after WP7.
- Face nesting from the rotation system plus one containment test per fixed point.
- Recipes and tests use `grow_until_bounded` where fixed step counts escape (NUM-11, NUM-15); `blast_zone(max_points)` exercised as a per-call budget (C6, TEST-7).

## Where to start

`numerics/TangleWorkbench.py`, `numerics/ManifoldMachine.py`, `numerics/diagnostics.py`, `loom/Blast.py`, `topology/DualGraph.py`, `examples/henon_cases.py`

## Binding context

- **Catalog:** `docs/open_issues.md` holds the evidence (path:line, probe numbers) for every issue ID named here. `docs/pipeline_chalkboard.md` is the one-board pipeline the code is judged against. Read both first.
- **Author decisions of 2026-10-07** (D1-D8, C1-C6 in `docs/open_issues.md` and `.refactor/decisions.md`) are binding. In short: topology is recomputed from scratch after any growth or blast, and the class count never decreases (D1); the alphabet comes from the minimal trellis only (D2); the shortest dual-graph walk IS the word and must always exist (D3); an inversion anchor is one point with two stable and two unstable branches, and cross-side is flagged only within one branch (D4); keep the period-k iterate rule but close its gap (D5); signs decide holes and rows (D6); heteroclinic pairs give holes in B's partition (D7); simplify freely (D8).
- `CLAUDE.md` invariants are physics. Topology never calls the map. `numerics` and `topology` never import `loom`.
- No golden pins in tests (policy 2026-10-05): new laws are rule-based; an expected failure is a strict xfail in `tests/cases.py` `KNOWN_ISSUES`, and `KNOWN_ISSUES` may only shrink.
- If the work is bigger than one targeted run, plan the first slice and name the rest as later runs under the same WP number.
