# Request: WP8 — Recompute semantics and stable identity

*From the author, 2026-10-07. Input for one `tangle_refactor` run (copy to `artifacts/request_origin_analyst.md`).*

## The question

**Under full recompute (D1), is every change in the topology a legitimate discovery, and can two recomputes be compared?** The author: "any growth that happens requires entirely recomputing the topology", with an eventual loop of blast → topology → blast → compare. Reference-pair switches (p3 blast 6, k28 blasts 8 and 11) are allowed; a DECREASE is a bug (k32 blasts 9-13, nested_outer blast 10, p3 blast 23). Elements are named by position, so every name churns (BCG-7).

## Issues in scope

BCG-4, BCG-5, BCG-7, API-1, API-2, SIMP-8 — see `docs/open_issues.md` for the evidence behind each.

## What "done" means

- Confirm with evidence that the switches follow the pseudoneighbor definition on the recomputed trellis; root-cause and fix every decrease.
- Elements are identified by their bounding holes `(origin, iterate, side)` (C2), so letters and classes compare across recomputes.
- The session holds ONE topology result per (generation, fixed points, pips), built in one pass, replacing the four caches (SIMP-8); it carries the pip choice across generations so `_repin` goes (API-1); gathering filters by the requested fixed points (API-2).

## Where to start

`loom/TangleSession.py`, `loom/BridgeAlphabet.py`, `topology/TopologyResults.py`, `topology/ElementNaming.py`, `topology/StablePartition.py`, `topology/Pseudoneighbor.py`, `examples/henon_cases.py`

## Binding context

- **Catalog:** `docs/open_issues.md` holds the evidence (path:line, probe numbers) for every issue ID named here. `docs/pipeline_chalkboard.md` is the one-board pipeline the code is judged against. Read both first.
- **Author decisions of 2026-10-07** (D1-D8, C1-C6 in `docs/open_issues.md` and `.refactor/decisions.md`) are binding. In short: topology is recomputed from scratch after any growth or blast, and the class count never decreases (D1); the alphabet comes from the minimal trellis only (D2); the shortest dual-graph walk IS the word and must always exist (D3); an inversion anchor is one point with two stable and two unstable branches, and cross-side is flagged only within one branch (D4); keep the period-k iterate rule but close its gap (D5); signs decide holes and rows (D6); heteroclinic pairs give holes in B's partition (D7); simplify freely (D8).
- `CLAUDE.md` invariants are physics. Topology never calls the map. `numerics` and `topology` never import `loom`.
- No golden pins in tests (policy 2026-10-05): new laws are rule-based; an expected failure is a strict xfail in `tests/cases.py` `KNOWN_ISSUES`, and `KNOWN_ISSUES` may only shrink.
- If the work is bigger than one targeted run, plan the first slice and name the rest as later runs under the same WP number.
