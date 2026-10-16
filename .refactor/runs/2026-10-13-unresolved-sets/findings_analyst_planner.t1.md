# Findings — B8 `unresolved-sets` (analyst, 2026-10-12, HEAD `8a1cdc3`)

*(The analyst's Write tool refused the artifact; the orchestrator saved this text verbatim.)*

**The tree is dirty.** Modified: `.refactor/decisions.md`, `.refactor/knowledge/{loom,numerics,testing,topology}.md`, `docs/alignment_plan_2026-10-09.md`, `docs/code_fidelity_audit_2026-10-09.md`. A rename is staged in the index (`artifacts/request_origin_analyst.md -> .refactor/runs/2026-10-12-oracle-laws/…`). Untracked: run archives and `artifacts/`. None of these touch the files B8 edits. The change branches from `HEAD` and must leave all of them alone, and must not commit the staged rename.

**What kind of request.** Two warn-and-continue sites become recorded data. Targets:
- `numerics/IterateInference.py`: `infer_iterates`, `_bracketed_image`.
- `numerics/IntersectionRegistry.py`, together with `numerics/IterateTable.py`.
- `loom/TangleSession.py`: `resonance_zone`, and the new `forget_strong_pip`.
- `loom/Topology.py`: `_choose_pip`.

My job is to check the plan against `HEAD`, not to re-plan.

## 1. The two sites at HEAD (verified, read)

| Plan says | HEAD | What it does |
|---|---|---|
| `IterateInference.py:125-132` collision | **`IterateInference.py:128-136`** | For each gated match `i -> j`: if `(j, -1)` is already in the table, it logs the WARNING "crossing %d: its image %d already has the preimage %d; not linked" and `continue`s. The earlier link stays. "Earlier" means either a previous pass's link, or a source of the same pass that comes first in `matches` dict order (detected-id order). |
| `:227` multi-match | **`IterateInference.py:232-240`** | Logs INFO "… meet its carried brackets (under-resolved); not linked", returns None, and records nothing. |
| `TangleSession.py:1509` | **`TangleSession.py:1509`** (in `resonance_zone` :1474-1510) | `self._strong_pips[rz.key] = pip` runs with no candidacy test. `pip` is the `Intersection` read before the trim. |
| `Topology._choose_pip` | **`Topology.py:113-127`**, raise at :121-125 | `ValueError("… is no longer a strong-pip candidate; choose another with set_strong_pip()")`. |

Context:
- `infer_iterates` reads N at :109 (`resolved_ids`, B1), so sources and candidates are both in N.
- Callers of `infer_iterates`:
  - `compute_intersections` (`TangleWorkbench.py:1389`), registry-wide.
  - `iterate_bridge`, **partial** (`BridgeIterator.py:398`: the bridge's two ends plus the crossings just born).
  - `iterate_all_bridges` (`BridgeIterator.py:441`), registry-wide.
  - the end of `blast_zone` (`Blast.py:463`), registry-wide.

## 2. Does a collision or a multi-match fire after B1? (verified, ran it)

**Method.** A probe attaches a counting handler to the `tanglepack.numerics.IterateInference` logger (DEBUG) and counts the collision WARNING, the multi-match INFO and the map-gate miss WARNING. It ran on:
- every `cases.BUILDERS` case: build, then a `compute_intersections(preserve_ids=True)` recompute;
- every step of all five `cases.SERIES`: k28 to 8 blasts, p3 to 9, inversion to 4, k32 to 10, nested_outer to 11, each with its per-step plain recompute;
- `build_k28(blasts=8)`: build and recompute.

**Results.**
- **Collisions: 0 everywhere.**
- **Multi-matches: 0 everywhere except `near_tangency`.** There, 25 distinct N sources fire it once per pass, 50 log lines over build plus recompute. Example: sources 25 and 27 both see `[61, 59]`.
  - `near_tangency` is the resolution-guard case at stage `intersected`. No law runs on it; only `tests/regression/test_resolution_warning_follows_oracle.py:43` builds it.
  - After B8 step 2 its registry will hold about 25 `unresolved_links` entries. Nothing reads it.
- **Map-gate misses: 0.**
- `orientation_reversing` raises at build. That is its intended `ValueError` and is irrelevant here.

**An order-free check as well.** For every N source I computed its pre-gate bracket images and counted claims per image:
- No image is claimed twice, and there is no multi-match, on k10, k28 x1/x2/x8 or near_tangency.
- This holds even with N replaced by the **whole registry**, unresolved detections included.
- So the collision is absent, not just hidden by the order the passes run in.

**The plan's k28 x8 collision (125 and 251 claiming 213) is stale.** Those ids are pre-B1 numbering. At HEAD, k28 x8 has 244 crossings, |N| = 239, unresolved `[75, 76, 77, 148, 245]` and 140 `+1` links; crossing 125 now links to 209 and 213 has no link. The plan's claim "the collision disappears after B1" is **true**, but it should be stated as "no collision anywhere (probed)", not tied to ids 125/251/213. The B1 outcome itself noted that nobody had checked this (`.refactor/runs/2026-10-12-resolved-crossings/outcome_auditor_human_2.md:161`).

**Positive control (verified).** On k28 x2:
1. Empty the table.
2. Hand-link a wrong `z -> 5`.
3. Run `infer_iterates([3])`; the true link is `3 -> 5`.

Result: the WARNING "crossing 3: its image 5 already has the preimage 8" fires, 3 stays unlinked, and 8 keeps 5. That is the order-decided state rubric 3 targets, built without stubs. This covers the cross-pass collision. The same-pass collision (two sources in one `matches` dict) needs a monkeypatched `_bracketed_image` (and `_map_gated`) returning one `j` for two sources.

## 3. Where the tests go (verified, read)

- **The iterate-table law is in `tests/helpers/laws.py`.**
  - `check_iterate_table_partial_bijection` :742, `check_iterate_links_pass_the_map_gate` :757, `check_iterate_table_complete` :770.
  - All are registered under layer `"crossings"` at :3066-3076, run by `tests/invariants/test_law_crossings.py:20`, and also called directly by `tests/regression/test_iterates_after_blasts.py:26-29`.
  - Step 5's law belongs in **`check_iterate_table_complete`**. Today it asserts that every unlinked detected crossing has no bracketed, gated image. It would therefore fail on any recorded multi-match or collision source, so B8 must exempt `unresolved_links` sources there in any case.
- **The inference test module is `tests/unit/numerics/test_iterate_inference.py`.** It exists, holds 2 tests, uses the fixture `henon_tangle_with_bridges` and `helpers.logs.assert_logged`. The "or" in the plan resolves to this file.
- **Facade:** `tests/facade/test_session_topology.py` exists (7 tests).
- Baseline: these modules plus `test_session_bridge_classes.py` and `test_iterates_after_blasts.py` give 44 passed at HEAD.

## 4. What differs from the plan section

1. **The file list is short by three.**
   - `tests/helpers/laws.py` is needed for step 5; the request already concedes this.
   - `src/tanglepack/numerics/IterateTable.py` is needed because the table has **no removal method**. `IterateTable`'s own docstring (:73-75) forbids touching `_forward`/`_backward` from outside. `IntersectionRegistry.register_iterate` (:654-674) also adds the graph edge `iter_{src}_{n}`, so `remove_iterate` must undo both.
   - `src/tanglepack/loom/Topology.py` is needed because step 4 changes `_choose_pip`'s `ValueError` text (:121-125) and its docstring (:97-99).
   - Rubric 1's "only the five files" is therefore unmeetable as written; it should be eight.
2. **Step 3 as worded is order-dependent for a third claimant (believed, by reading the plan).** Within one pass:
   - the 1st source links;
   - the 2nd hits `(j, -1)`, so both are recorded and the 1st is unlinked;
   - the 3rd then finds `(j, -1)` empty and **links**.

   Across passes, a partial pass (`iterate_bridge`, sources = bridge ends and new crossings) that clears its own sources' entries can relink an unlinked loser while the other claimant is outside the pass. Every library path ends with a registry-wide pass, which would catch that; a user's bare `iterate_bridge` would not. The step 3 text also says `remove_iterate(i, 1)` where `i` is the new source; the earlier source is `table[j, -1]`.
3. **Step 1's generation bump.** "Cleared at the start of each pass … bumping `_generation` when it changes" would bump on every pass that clears and re-records the same entries (every near_tangency pass). The bump should compare before and after. This does not affect law-case fingerprints: nothing is recorded on them. `compute_intersections` builds a fresh registry, so the set starts empty after a recompute and nothing stale survives renumbering.
4. **Step 4's `forget_strong_pip` must also reset the trellis (verified).** Deleting `_strong_pips[fp]` alone does not forget: the cached trellis keeps `strong_pip`, and `_choose_pip` (:118-120) re-adopts it through `set_strong_pip`. On k10, set 3, 4 or 7, then a dict-only "forget", then `topology()` gives the same pip back. With `self.trellis(fp).strong_pip = None` as well, it gives the default 1.
5. **Rubric 4 works only on a specific geometry (verified).** A zone at a non-candidate usually trims the chosen pip **away**, because the non-candidate lies inward of the pip on the same branch. `topology()` then raises with or without B8's adoption:
   - verified on k10 (5, 6), k28 x2 (9, 10), nested (37, 38), nested_inner (41, 42), heteroclinic, inversion x1 (26, 27) and inversion x2 (23-27);
   - on k10 at stage `bridges`, with no pip chosen, `topology()` even asserts.

   The one realisation found: **`inversion_two_blasts`, a non-candidate on the other stable branch from the chosen pip** (32 or 33 on `(stable, 0, 1)`; pip 1 is on `(0, 0)`).
   - On base, `session.resonance_zone(32)` adopts it and `topology()` raises `ValueError`.
   - With the adoption simply skipped, `topology()` runs and the pip stays 1.

   Choose it by rule (a non-candidate whose stable branch differs from the chosen pip's), not by id, per the test policy.
6. **Rubric 5 needs a candidate that `topology()` can run on (verified).** Some candidates assert in `check_bridge_rows_consistent` once set: equal signs on a bridge joined across the cut. These are k10 2, k28 x1/x2 4, and inversion 2 and 14. k10 3, 4 or 7 are safe.
7. **Every recipe's zone pip is a candidate both before and after its zone (verified).** I wrapped `resonance_zone` across all BUILDERS. So the candidacy test changes no recipe, and the plan's "no recipe adopts a non-candidate" holds.
   - The test that matches `_restore_strong_pip` and guarantees no later raise is the **post-zone** one: `t = self.trellis(rz.key); found = t.registry.find(pip); found is not None and t._is_candidate(found)`.
   - The existing stale-pip tests (`test_session_bridge_classes.py:119-170`) call `define_resonance_zone` directly and match "no longer a strong-pip candidate", so the new `ValueError` text must keep that phrase.

## 5. Answer and smallest change

**Answer.** Both sites are where the plan says, a few lines lower for the collision and the multi-match. After B1:
- the collision fires nowhere in the corpus, not even counting unresolved detections, so B8's collision rule is a safety net;
- the multi-match fires only on `near_tangency`, which no law reads;
- no recipe adopts a non-candidate pip.

The pip remedy is real: on base there is no way to forget a pip, and a zone at a non-candidate makes `topology()` raise.

**Recommendation.** Keep slice B8, with these adjustments:
- **(a) File list:** add `numerics/IterateTable.py` (`remove`/`__delitem__` for both directions), `loom/Topology.py` (the error text) and `tests/helpers/laws.py`. Rubric 1 becomes eight files.
- **(b) Step 3, order-free:** after gating, group the confirmed matches by image. An image named by two or more sources of the pass, or by one source while `(j, -1)` is already held, records every claimant (the pass's sources plus `table[j, -1]`) as `(j,)`, and unlinks the held one. No claimant links.
- **(c) Step 5 law:** add one clause, "no unresolved source with a single candidate `j` coexists with a `(j, -1)` link". Together with "disjoint from linked sources" and "≥ 1 candidate", it catches the third-claimant and partial-pass relinks. Put it in `check_iterate_table_complete`, which must also exempt unresolved sources.
- **(d) `forget_strong_pip`:** pop `_strong_pips[fp]` **and** set `self.trellis(fp).strong_pip = None`.
- **(e) `resonance_zone`:** test candidacy post-zone, as `_restore_strong_pip` does.
- **(f) Rubric 4 and 5 fixtures:** as in §4.5 and §4.6.
- **(g) Generation bump:** compare the dict before and after the pass.

## Out of scope (noticed, not for B8)
- Setting some legitimate candidates as the pip makes `topology()` **assert**: equal crossing signs on a bridge joined across the pip-orbit cut. Affected: k10 2, k28 x1/x2 4, inversion 2 and 14. That is an assertion, not a flag, reachable through the public `set_strong_pip`.
- A zone at a non-candidate inward of the chosen pip trims the pip away, and `topology()` still raises after B8. `forget_strong_pip` is then the remedy, which is what step 4's message should say.
- `check_iterate_table_complete` scans the whole registry, not N (`laws.py:779-786`). Harmless on the law cases, where N equals the registry.

## Probes
All in `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/`:
- `probe.py` (log counts): output in `builders.log`, `series1.log` … `series3.log`, `k28x8.log`;
- `static.py` (order-free census);
- `pip.py`, `pip2.py`, `pip3.py`, `scan_nc.py` (pip behaviour).

Knowledge files updated:
- `.refactor/knowledge/numerics.md`: collision and multi-match status, HEAD line numbers, IterateTable has no removal, fresh registry on recompute.
- `.refactor/knowledge/loom.md`: pip facts (forget must reset the trellis, non-candidate zone geometry, candidates that assert).
- `.refactor/knowledge/testing.md`: law and test locations for B8.

During the work, three probe logs were briefly written into the repo root by a cwd slip. They have been moved out, and no tracked file was touched.
