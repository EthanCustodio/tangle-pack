# Request: WP3 — Inversion end to end

*From the author, 2026-10-07. Input for one `tangle_refactor` run (copy to `artifacts/request_origin_analyst.md`).*

## The question

**Does a saddle with inversion go through the whole pipeline correctly?** The author's model (D4): "Each anchor point should have two stable branches and two unstable branches. It is that simple. Each of those individual branches of course have a left side and right side." The left of one branch is the right of the other, and zones break the plane up. Today the code registers four coincident anchors (INV-1), an anchor-only branch asserts (INV-3), the zone ring is keyed by orbit index (INV-5), cross-side is flagged between branches (INV-6), and branch ordering has two rules (INV-7).

## Issues in scope

INV-1, INV-3, INV-5, INV-6, INV-7, INV-8, INV-9, INV-10, API-5 — see `docs/open_issues.md` for the evidence behind each.

## What "done" means

- One anchor per periodic point, with its four branches; `Tangle._both_anchors` goes; the arrangement, image chains and the `[0,0]` singleton see one node.
- An anchor-only stable branch is one element `[anchor, end]` (C1).
- The zone ring is keyed by the full `ManifoldKey`; the inversion zone is built correctly from both loops; `TangleSession._zone_ring_order` reuses `ResonanceZone`.
- Cross-side is flagged only for a pair on one branch.
- One branch-ordering rule (`advance_key`); `get_iterable_array` goes; `orient_eigenvectors` works for period>1 inversion; the `grow_n_times(branch_index)` log noise stops.
- Re-check INV-10 (`iterate_bridge` on an anchor-rooted inversion bridge).
- WP1's k=4 inversion laws turn green.

## Where to start

`numerics/TangleWorkbench.py` (`_register_anchors`), `numerics/Tangle.py`, `numerics/ManifoldInitializer.py`, `numerics/FixedPoint.py`, `topology/StablePartition.py`, `topology/Arrangement.py`, `topology/BridgeClass.py`, `topology/SymbolicDynamics.py`, `loom/ResonanceZone.py`, `loom/TangleSession.py`

## Binding context

- **Catalog:** `docs/open_issues.md` holds the evidence (path:line, probe numbers) for every issue ID named here. `docs/pipeline_chalkboard.md` is the one-board pipeline the code is judged against. Read both first.
- **Author decisions of 2026-10-07** (D1-D8, C1-C6 in `docs/open_issues.md` and `.refactor/decisions.md`) are binding. In short: topology is recomputed from scratch after any growth or blast, and the class count never decreases (D1); the alphabet comes from the minimal trellis only (D2); the shortest dual-graph walk IS the word and must always exist (D3); an inversion anchor is one point with two stable and two unstable branches, and cross-side is flagged only within one branch (D4); keep the period-k iterate rule but close its gap (D5); signs decide holes and rows (D6); heteroclinic pairs give holes in B's partition (D7); simplify freely (D8).
- `CLAUDE.md` invariants are physics. Topology never calls the map. `numerics` and `topology` never import `loom`.
- No golden pins in tests (policy 2026-10-05): new laws are rule-based; an expected failure is a strict xfail in `tests/cases.py` `KNOWN_ISSUES`, and `KNOWN_ISSUES` may only shrink.
- If the work is bigger than one targeted run, plan the first slice and name the rest as later runs under the same WP number.
