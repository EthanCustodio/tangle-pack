# Request: WP5 — Period-k: no branch silently without holes

*From the author, 2026-10-07. Input for one `tangle_refactor` run (copy to `artifacts/request_origin_analyst.md`).*

## The question

**Does every branch of a period-k orbit get its fundamental-segment holes?** D5 keeps the +1..+(k−1) iterate rule but closes the gap: at p3 blast 6, `extend_pseudoneighbor_trajectories` found no table links, so branches 0.0 and 2.0 got no direct holes and the partition became asymmetric across the orbit (PK-3).

## Issues in scope

PK-2, PK-3, PK-6, PK-8, PK-9 — see `docs/open_issues.md` for the evidence behind each.

## What "done" means

- A missing +j link is placed by `Trellis.image_cdist` (table first, then scaling; still lookup only), never skipped silently; if even that fails, a WARNING names the branch left without holes.
- A strong pip with fewer than k cut points is an explicit flag (PK-6); the `max_steps` hole-chain truncation logs a WARNING (PK-9).
- The expected "inert class maps to the next branch's inert class" warning is silenced (PK-8).
- Hole counts are symmetric across the orbit on p3 at 0-8 blasts.

## Where to start

`topology/StablePartition.py` (`extend_pseudoneighbor_trajectories`, `_is_forward_beyond_fundamental`), `topology/Pseudoneighbor.py`, `topology/SymbolicDynamics.py`

## Binding context

- **Catalog:** `docs/open_issues.md` holds the evidence (path:line, probe numbers) for every issue ID named here. `docs/pipeline_chalkboard.md` is the one-board pipeline the code is judged against. Read both first.
- **Author decisions of 2026-10-07** (D1-D8, C1-C6 in `docs/open_issues.md` and `.refactor/decisions.md`) are binding. In short: topology is recomputed from scratch after any growth or blast, and the class count never decreases (D1); the alphabet comes from the minimal trellis only (D2); the shortest dual-graph walk IS the word and must always exist (D3); an inversion anchor is one point with two stable and two unstable branches, and cross-side is flagged only within one branch (D4); keep the period-k iterate rule but close its gap (D5); signs decide holes and rows (D6); heteroclinic pairs give holes in B's partition (D7); simplify freely (D8).
- `CLAUDE.md` invariants are physics. Topology never calls the map. `numerics` and `topology` never import `loom`.
- No golden pins in tests (policy 2026-10-05): new laws are rule-based; an expected failure is a strict xfail in `tests/cases.py` `KNOWN_ISSUES`, and `KNOWN_ISSUES` may only shrink.
- If the work is bigger than one targeted run, plan the first slice and name the rest as later runs under the same WP number.
