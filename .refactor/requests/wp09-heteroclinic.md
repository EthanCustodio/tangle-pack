# Request: WP9 — Nested and heteroclinic tangles

*From the author, 2026-10-07. Input for one `tangle_refactor` run (copy to `artifacts/request_origin_analyst.md`).*

## The question

**Does a nested case with heteroclinic crossings run end to end?** No fixture has an A(u)×B(s) crossing (NEST-1); heteroclinic pseudoneighbor pairs are rejected by construction (NEST-2), while D7 says they give holes in B's partition. The audit's hand-built case (p3 at 16u/12s plus p1 at 11u, area_cutoff 1e-7, about 1 s, 4 heteroclinic crossings) ran but gave 2 virtual classes, a cross-side pair and a standalone connecting loop; at 1e-5 it asserted.

## Issues in scope

NEST-1, NEST-2, NEST-4, NEST-5, NEST-6 — see `docs/open_issues.md` for the evidence behind each.

## What "done" means

- A recipe with heteroclinic crossings in `henon_cases.py`, and a law case for it.
- Heteroclinic pairs accepted into B's partition (`Pseudoneighbor.py:535-540`).
- Blasting classifies by innermost zone (C4); `rebuild_bridges(fp)` clears only that fixed point's bridges (NEST-5); `entry.tangle` order agrees with `TangleBuild.fixed_points` (NEST-6).
- Connecting classes, the merged outer face and cross-tangle walks run and satisfy the laws.

## Where to start

`examples/henon_cases.py`, `topology/Pseudoneighbor.py`, `topology/DualGraph.py`, `loom/Blast.py`, `numerics/TangleWorkbench.py`, `tests/cases.py`

## Binding context

- **Catalog:** `docs/open_issues.md` holds the evidence (path:line, probe numbers) for every issue ID named here. `docs/pipeline_chalkboard.md` is the one-board pipeline the code is judged against. Read both first.
- **Author decisions of 2026-10-07** (D1-D8, C1-C6 in `docs/open_issues.md` and `.refactor/decisions.md`) are binding. In short: topology is recomputed from scratch after any growth or blast, and the class count never decreases (D1); the alphabet comes from the minimal trellis only (D2); the shortest dual-graph walk IS the word and must always exist (D3); an inversion anchor is one point with two stable and two unstable branches, and cross-side is flagged only within one branch (D4); keep the period-k iterate rule but close its gap (D5); signs decide holes and rows (D6); heteroclinic pairs give holes in B's partition (D7); simplify freely (D8).
- `CLAUDE.md` invariants are physics. Topology never calls the map. `numerics` and `topology` never import `loom`.
- No golden pins in tests (policy 2026-10-05): new laws are rule-based; an expected failure is a strict xfail in `tests/cases.py` `KNOWN_ISSUES`, and `KNOWN_ISSUES` may only shrink.
- If the work is bigger than one targeted run, plan the first slice and name the rest as later runs under the same WP number.
