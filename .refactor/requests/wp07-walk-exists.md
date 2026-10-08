# Request: WP7 — The dual-graph walk always exists

*From the author, 2026-10-07. Input for one `tangle_refactor` run (copy to `artifacts/request_origin_analyst.md`).*

## The question

**Why do walks come back unreachable or ambiguous, and where do virtual `new*` classes come from?** D3: "There should always be a walk on the dual graph for the iterate. Otherwise something is wrong with how it is constructed. … the virtual 'new1' class thing … is a very clear bug. In general how the period 3 stuff is working seems buggy." Today `open_p3_deep_runs` (p3 at 15/16 steps, 6 blasts) is a strict xfail: "crossing abuts no hole" (`PartitionFamily.py:800`), cross-branch "no flank" skips, fundamental-segment unification on the pip branch only, straddling landings, and pairs no `Bridge` spans.

## Issues in scope

PK-1, PK-4, PK-5, INT-1 (resolved by D3: walk stays), INT-8, INT-9, NEST-3 — see `docs/open_issues.md` for the evidence behind each.

## What "done" means

- Root-cause each failure on p3 15/16 steps, p3 6+ blasts and the nested runs, and fix the construction so a walk exists.
- Delete virtual-class support: a pair no bridge spans gets its Bridge, or is a construction error.
- Singletons keep the f(bridge) path. Elsewhere, a walk disagreeing with the registered image chain at word level is a WARNING with a diagnosis object and counts against `is_reliable`.
- `open_p3_deep_runs` leaves `KNOWN_ISSUES`; WP1's walk-existence law turns green on every case.

## Where to start

`topology/PartitionFamily.py`, `topology/DualGraph.py`, `topology/DualWalk.py`, `topology/SymbolicDynamics.py`, `topology/MinimalTrellis.py`, `numerics/BridgeIterator.py`, `tests/regression/test_open_p3_deep_runs.py`

## Binding context

- **Catalog:** `docs/open_issues.md` holds the evidence (path:line, probe numbers) for every issue ID named here. `docs/pipeline_chalkboard.md` is the one-board pipeline the code is judged against. Read both first.
- **Author decisions of 2026-10-07** (D1-D8, C1-C6 in `docs/open_issues.md` and `.refactor/decisions.md`) are binding. In short: topology is recomputed from scratch after any growth or blast, and the class count never decreases (D1); the alphabet comes from the minimal trellis only (D2); the shortest dual-graph walk IS the word and must always exist (D3); an inversion anchor is one point with two stable and two unstable branches, and cross-side is flagged only within one branch (D4); keep the period-k iterate rule but close its gap (D5); signs decide holes and rows (D6); heteroclinic pairs give holes in B's partition (D7); simplify freely (D8).
- `CLAUDE.md` invariants are physics. Topology never calls the map. `numerics` and `topology` never import `loom`.
- No golden pins in tests (policy 2026-10-05): new laws are rule-based; an expected failure is a strict xfail in `tests/cases.py` `KNOWN_ISSUES`, and `KNOWN_ISSUES` may only shrink.
- If the work is bigger than one targeted run, plan the first slice and name the rest as later runs under the same WP number.
