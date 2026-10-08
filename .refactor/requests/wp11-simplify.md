# Request: WP11 — Simplification sweep

*From the author, 2026-10-07. Input for one `tangle_refactor` run (copy to `artifacts/request_origin_analyst.md`).*

## The question

**Where is the code still more complex than the chalkboard?** With the correctness WPs merged, sweep what remains: six orientation-flip sites (SIMP-1, C3: one helper), two element-image implementations (SIMP-2), a copied `_build_intervals` (SIMP-5), two stable trims (SIMP-10), the `Blast.py` tidy-up (SIMP-7), a duplicated cartoon block (SIMP-11), dead API (SIMP-12), API traps (API-3, 4, 6, 7).

## Issues in scope

INV-11, SIMP-1, SIMP-2, SIMP-4, SIMP-5, SIMP-7, SIMP-9, SIMP-10, SIMP-11, SIMP-12, SIMP-13, API-3, API-4, API-6, API-7, TEST-6 — see `docs/open_issues.md` for the evidence behind each.

## What "done" means

- Each item is shorter or clearer after, and the suite passes.
- Report tolerance-site and WARNING counts in `topology/` before and after (SIMP-13).
- `coverage_guard.py` stops tripping when lines shift (TEST-6).

## Where to start

spread across `topology/`, `numerics/`, `loom/Blast.py`, `topology/plotting.py`, `tests/conftest.py`, `tests/_tools/coverage_guard.py`

## Binding context

- **Catalog:** `docs/open_issues.md` holds the evidence (path:line, probe numbers) for every issue ID named here. `docs/pipeline_chalkboard.md` is the one-board pipeline the code is judged against. Read both first.
- **Author decisions of 2026-10-07** (D1-D8, C1-C6 in `docs/open_issues.md` and `.refactor/decisions.md`) are binding. In short: topology is recomputed from scratch after any growth or blast, and the class count never decreases (D1); the alphabet comes from the minimal trellis only (D2); the shortest dual-graph walk IS the word and must always exist (D3); an inversion anchor is one point with two stable and two unstable branches, and cross-side is flagged only within one branch (D4); keep the period-k iterate rule but close its gap (D5); signs decide holes and rows (D6); heteroclinic pairs give holes in B's partition (D7); simplify freely (D8).
- `CLAUDE.md` invariants are physics. Topology never calls the map. `numerics` and `topology` never import `loom`.
- No golden pins in tests (policy 2026-10-05): new laws are rule-based; an expected failure is a strict xfail in `tests/cases.py` `KNOWN_ISSUES`, and `KNOWN_ISSUES` may only shrink.
- If the work is bigger than one targeted run, plan the first slice and name the rest as later runs under the same WP number.
