# Request: WP2 — Exact iterate bookkeeping

*From the author, 2026-10-07. Input for one `tangle_refactor` run (copy to `artifacts/request_origin_analyst.md`).*

## The question

**Is the iterate table — the central combinatorial object — exact, and is cdist scaling done one way?** Blasts never re-run `infer_iterates`, so fold-born crossings stay unlinked (BCG-3; running it after `build_period3(blasts=10)` removed 2 classes and 2 orphan loops). Inference matches on a 5% relative cdist with no map gate (INT-6). Scaling is implemented five ways, and the pseudoneighbor lap factor is `lambda_u` instead of β^k (INV-2).

## Issues in scope

BCG-3, INT-6, INT-7, NUM-8, SIMP-3, PK-7, INV-2 — see `docs/open_issues.md` for the evidence behind each.

## What "done" means

- `infer_iterates` (or exact registration) runs after every blast, gated by one forward map call in numerics; the 5% heuristic is no longer the only test.
- `cdist_tol` is relative (NUM-8).
- `Trellis.image_cdist` is the ONE scaling path; `scale_cdist`, `TrellisBranch.eigenvalue`, `_get_lambda_u` and the topology-side re-implementations go; the lap factor is β^k (`Pseudoneighbor.py:491,551,576`, `IntersectionRegistry.on_interval`).
- The dead `manifold_a_key is None` branches go (PK-7).
- Crossing ids: say whether exact registration makes `reindex_from` matching unnecessary (INT-7); fix only if small.

## Where to start

`numerics/IterateInference.py`, `numerics/IntersectionRegistry.py`, `numerics/BridgeIterator.py`, `loom/Blast.py`, `topology/Trellis.py`, `topology/TrellisBranch.py`, `topology/Pseudoneighbor.py`, `topology/StrongPip.py`

## Binding context

- **Catalog:** `docs/open_issues.md` holds the evidence (path:line, probe numbers) for every issue ID named here. `docs/pipeline_chalkboard.md` is the one-board pipeline the code is judged against. Read both first.
- **Author decisions of 2026-10-07** (D1-D8, C1-C6 in `docs/open_issues.md` and `.refactor/decisions.md`) are binding. In short: topology is recomputed from scratch after any growth or blast, and the class count never decreases (D1); the alphabet comes from the minimal trellis only (D2); the shortest dual-graph walk IS the word and must always exist (D3); an inversion anchor is one point with two stable and two unstable branches, and cross-side is flagged only within one branch (D4); keep the period-k iterate rule but close its gap (D5); signs decide holes and rows (D6); heteroclinic pairs give holes in B's partition (D7); simplify freely (D8).
- `CLAUDE.md` invariants are physics. Topology never calls the map. `numerics` and `topology` never import `loom`.
- No golden pins in tests (policy 2026-10-05): new laws are rule-based; an expected failure is a strict xfail in `tests/cases.py` `KNOWN_ISSUES`, and `KNOWN_ISSUES` may only shrink.
- If the work is bigger than one targeted run, plan the first slice and name the rest as later runs under the same WP number.
