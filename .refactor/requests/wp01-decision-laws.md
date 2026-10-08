# Request: WP1 — Laws that encode the 2026-10-07 decisions

*From the author, 2026-10-07. Input for one `tangle_refactor` run (copy to `artifacts/request_origin_analyst.md`).*

## The question

**Can the test suite tell us, by rule, when the decisions of 2026-10-07 hold?** Today 833 tests pass while p3's class count grows by two per blast, inversion has four anchors and walks come back unreachable. This run is tests first: it adds the laws later WPs must turn green, and nothing in `src/` changes except what a law needs to read.

## Issues in scope

INV-4, TEST-1, TEST-2, TEST-3, TEST-4, TEST-5 — see `docs/open_issues.md` for the evidence behind each.

## What "done" means

- **Walk existence (D3):** every active, non-singleton class has a walk whose status is `unique` or `trivial`; no `unreachable`, no virtual `new*` class.
- **Cross-side within one branch only (D4):** replaces `itinerary_pairs_same_side`; a pair across two branches is never flagged.
- **Class monotonicity (D1, D2):** in `tests/regression/`, along a blast series (k28 and p3, enough blasts to cover the p3 switch at blast 6), the active and inert class counts never decrease, and they are constant while the hole set is constant.
- **Alphabet ⊆ minimal trellis (D2):** every lettered class has a member among `minimal_trellis().kept_bridge_ids`.
- **New inversion fixture (INV-4):** the k=4 other-root saddle `(1-sqrt5, sqrt5-1)`, `area_cutoff=1e-5`, a zone at f(q0), at 0/1/2 blasts, frozen in `tests/cases.py` and added to `LAW_CASES` in place of the hole-free k=10 case.
- **NOT_APPLICABLE as visible skips (TEST-3).**
Each law that fails today goes in `KNOWN_ISSUES` with the WP that will fix it.

## Where to start

`tests/cases.py`, `tests/helpers/laws.py`, `tests/helpers/law_tier.py`, `tests/conftest.py`, `tests/regression/test_blast_monotonicity.py`, `tests/invariants/`

## Binding context

- **Catalog:** `docs/open_issues.md` holds the evidence (path:line, probe numbers) for every issue ID named here. `docs/pipeline_chalkboard.md` is the one-board pipeline the code is judged against. Read both first.
- **Author decisions of 2026-10-07** (D1-D8, C1-C6 in `docs/open_issues.md` and `.refactor/decisions.md`) are binding. In short: topology is recomputed from scratch after any growth or blast, and the class count never decreases (D1); the alphabet comes from the minimal trellis only (D2); the shortest dual-graph walk IS the word and must always exist (D3); an inversion anchor is one point with two stable and two unstable branches, and cross-side is flagged only within one branch (D4); keep the period-k iterate rule but close its gap (D5); signs decide holes and rows (D6); heteroclinic pairs give holes in B's partition (D7); simplify freely (D8).
- `CLAUDE.md` invariants are physics. Topology never calls the map. `numerics` and `topology` never import `loom`.
- No golden pins in tests (policy 2026-10-05): new laws are rule-based; an expected failure is a strict xfail in `tests/cases.py` `KNOWN_ISSUES`, and `KNOWN_ISSUES` may only shrink.
- If the work is bigger than one targeted run, plan the first slice and name the rest as later runs under the same WP number.
