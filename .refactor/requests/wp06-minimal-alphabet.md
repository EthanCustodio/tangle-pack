# Request: WP6 — The alphabet comes from the minimal trellis

*From the author, 2026-10-07. Input for one `tangle_refactor` run (copy to `artifacts/request_origin_analyst.md`).*

## The question

**Why does the class count grow by about two per blast, and does building the alphabet from the minimal trellis stop it (D2)?** `bridge_classes` classes every bridge, blast children included, so each blast adds one class with no image (BCG-1). Loops are folded through endpoint preimages that fold-born crossings do not have, so each blast adds a standalone inert `BridgeClass(x,x)` (BCG-2). Over p3 blasts 8-16 the holes, elements and minimal trellis are constant while classes go 27 → 43.

## Issues in scope

BCG-1, BCG-2, BCG-6, SIMP-6 — see `docs/open_issues.md` for the evidence behind each.

## What "done" means

- Letters, words and the transition matrix are built from classes of `minimal_trellis().kept_bridge_ids`; other bridges are image evidence only.
- A loop folds FORWARD into the bridge whose `_image_chain` contains it (fallback `BridgeIterator.preimage_bridges`); true orphans are logged, never lettered.
- One letter generator for the active and inert series (SIMP-6).
- Rerun the E11 p3 and k28 series (small): flat between reference switches. WP1's monotonicity and alphabet laws turn green for this mechanism.

## Where to start

`topology/BridgeClass.py`, `topology/MinimalTrellis.py`, `topology/SymbolicDynamics.py`, `loom/TangleSession.py`, `loom/BridgeAlphabet.py`

## Binding context

- **Catalog:** `docs/open_issues.md` holds the evidence (path:line, probe numbers) for every issue ID named here. `docs/pipeline_chalkboard.md` is the one-board pipeline the code is judged against. Read both first.
- **Author decisions of 2026-10-07** (D1-D8, C1-C6 in `docs/open_issues.md` and `.refactor/decisions.md`) are binding. In short: topology is recomputed from scratch after any growth or blast, and the class count never decreases (D1); the alphabet comes from the minimal trellis only (D2); the shortest dual-graph walk IS the word and must always exist (D3); an inversion anchor is one point with two stable and two unstable branches, and cross-side is flagged only within one branch (D4); keep the period-k iterate rule but close its gap (D5); signs decide holes and rows (D6); heteroclinic pairs give holes in B's partition (D7); simplify freely (D8).
- `CLAUDE.md` invariants are physics. Topology never calls the map. `numerics` and `topology` never import `loom`.
- No golden pins in tests (policy 2026-10-05): new laws are rule-based; an expected failure is a strict xfail in `tests/cases.py` `KNOWN_ISSUES`, and `KNOWN_ISSUES` may only shrink.
- If the work is bigger than one targeted run, plan the first slice and name the rest as later runs under the same WP number.
