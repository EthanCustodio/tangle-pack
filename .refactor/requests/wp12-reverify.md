# Request: WP12 — Re-verify and document

*From the author, 2026-10-07. Input for one `tangle_refactor` run (copy to `artifacts/request_origin_analyst.md`).*

## The question

**After WP1-WP11, does every case run end to end, and does the class count plateau?** Rerun the stress subset and the figures, and bring the docs up to date.

## Issues in scope

NUM-5, NUM-6, NUM-7, PERF-2, PERF-5, DOC-1, DOC-2, DOC-3, DOC-4, DOC-5 — see `docs/open_issues.md` for the evidence behind each.

## What "done" means

- E7 and E11 subsets (k28, k32, p3, nested_inner, nested_outer, inversion k=4, heteroclinic): class counts non-decreasing and flat between switches; no +1/+2-per-blast staircase; `is_reliable` True, no `new*`, no unreachable walks. Re-measure NUM-5, NUM-6, NUM-7, PERF-2.
- Figures regenerated for k10, k28 (1/2 blasts), p3 (13 and 16 steps, 0 and 8 blasts), nested (inner blasts > 4), inversion k=4 (0/1/2 blasts) and the heteroclinic case.
- CLAUDE.md, `.refactor/knowledge/*.md`, memory and the report framing updated (DOC-1..5); each closed issue marked in `docs/open_issues.md`.

## Where to start

`stress/numerics/`, `scripts/`, `CLAUDE.md`, `.refactor/knowledge/`, `docs/open_issues.md`

## Binding context

- **Catalog:** `docs/open_issues.md` holds the evidence (path:line, probe numbers) for every issue ID named here. `docs/pipeline_chalkboard.md` is the one-board pipeline the code is judged against. Read both first.
- **Author decisions of 2026-10-07** (D1-D8, C1-C6 in `docs/open_issues.md` and `.refactor/decisions.md`) are binding. In short: topology is recomputed from scratch after any growth or blast, and the class count never decreases (D1); the alphabet comes from the minimal trellis only (D2); the shortest dual-graph walk IS the word and must always exist (D3); an inversion anchor is one point with two stable and two unstable branches, and cross-side is flagged only within one branch (D4); keep the period-k iterate rule but close its gap (D5); signs decide holes and rows (D6); heteroclinic pairs give holes in B's partition (D7); simplify freely (D8).
- `CLAUDE.md` invariants are physics. Topology never calls the map. `numerics` and `topology` never import `loom`.
- No golden pins in tests (policy 2026-10-05): new laws are rule-based; an expected failure is a strict xfail in `tests/cases.py` `KNOWN_ISSUES`, and `KNOWN_ISSUES` may only shrink.
- If the work is bigger than one targeted run, plan the first slice and name the rest as later runs under the same WP number.
