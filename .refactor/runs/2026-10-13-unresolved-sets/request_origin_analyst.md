# Request: B8 — `unresolved-sets`

*From the author, 2026-10-12 ("I say so"). Input for one `tangle_refactor` run. The fourth slice of `docs/alignment_plan_2026-10-09.md`; B0, B1 and B6a are merged (main `8a1cdc3`; runs archived under `.refactor/runs/2026-10-12-{oracle-laws,resolved-crossings,walk-bounds-honest}/`).*

## The question

**Where does the code decide something the definition leaves undecided, and continue?** Two warn-and-continue sites (plan B8): an iterate-table collision is resolved by registry order (`IterateInference.py:125-132` at the plan's numbering: two sources whose brackets name one image, and the first registered wins), and a zone's pip is adopted without a candidacy test (`TangleSession.py:1509`) while a pip that stops being a candidate makes `topology()` raise with no API to forget it (`Topology._choose_pip`). Both become DATA instead of a silent choice: `registry.unresolved_links` and `TangleSession.forget_strong_pip`.

## What "done" means

Exactly slice B8 of the alignment plan (`docs/alignment_plan_2026-10-09.md`, section B, "B8 — `unresolved-sets`"): its scope, five steps, behaviour, not-doing list and six rubric items are the plan. In short:

1. `IntersectionRegistry.unresolved_links: dict[int, tuple[int, ...]]` (source id → the candidate image ids), cleared at the start of each `infer_iterates` pass over its sources, bumping `_generation` when it changes.
2. `_bracketed_image`: a multi-match records `unresolved_links[src] = tuple(images)` (INFO kept).
3. The collision: when `(j, -1)` is already in the table, record BOTH sources under `j`'s candidates (`unresolved_links[i] = (j,)`, `unresolved_links[table[j, -1]] = (j,)`) and UNLINK the earlier one (`register_iterate` gains `remove_iterate(i, 1)`), so the table never holds a link decided by order. One WARNING naming the three ids.
4. `TangleSession.forget_strong_pip(fp)`: drops the remembered pip; the next `topology()` chooses the default again. `resonance_zone()` remembers the zone's pip only when `trellis._is_candidate(pip)` and WARNS otherwise (the zone is still defined). `Topology._choose_pip`'s `ValueError` names `forget_strong_pip`.
5. Law (layer `crossings`, in the existing iterate-table check): `unresolved_links` is disjoint from the linked sources, and every unresolved source has ≥ 1 candidate.

Files (five): `src/tanglepack/numerics/IterateInference.py`, `src/tanglepack/numerics/IntersectionRegistry.py`, `src/tanglepack/loom/TangleSession.py`, `tests/unit/numerics/test_iterate_inference.py` (or the existing inference test module), `tests/facade/test_session_topology.py`. The law of step 5 lives in `tests/helpers/laws.py` if the existing iterate-table check is there; the analyst says where, and the planner lists the file.

**Behaviour.** Law cases: none (no collision fires; no recipe adopts a non-candidate pip). Since B1 (`N` = the registry's `resolved_ids`), the k28 x8 collision the plan describes (sources 125 and 251 claiming image 213) disappears with the unresolved detections, and the plan expects NO collision anywhere in the corpus: B8's collision rule is then a safety net the definition keeps honest, not a fix; the analyst verifies that. `resonance_zone()` on a non-candidate: a WARNING and no pip adoption instead of a later `ValueError` from `topology()`.

**Not doing.** Odd foreign-chord counts and I2 on image bridges (`PartitionFamily.py:760, 960-967`) are B7's; marker disagreement was B6a.

**Rubric** (the plan's six): suite passes on the branch and the diff touches only the listed files; `grep -n "already has the preimage" src` finds the one WARNING followed by `remove_iterate` / recording, no `continue` that leaves the earlier link in place; unit test: two sources whose brackets name one image → neither linked, both in `unresolved_links`, failing on base; facade test: `resonance_zone(id)` with a non-candidate `id` → `chosen_strong_pip(fp)` unchanged, one WARNING, `topology()` does not raise, failing on base; facade test: `set_strong_pip`, then `forget_strong_pip`, then `topology()` → the default pip; fingerprints of every law case equal on base and branch.

## Where to start

`docs/alignment_plan_2026-10-09.md` (B8; A.3 for the iterate table; A.4 for the pip), `docs/code_fidelity_audit_2026-10-09.md` (the collision and pip findings), `src/tanglepack/numerics/IterateInference.py` (`infer_iterates`, `_bracketed_image`, the collision branch), `src/tanglepack/numerics/IntersectionRegistry.py` (`register_iterate`, `IterateTable`), `src/tanglepack/loom/TangleSession.py` (`resonance_zone`, `set_strong_pip`, `chosen_strong_pip`), `src/tanglepack/loom/Topology.py` (`_choose_pip`), `tests/helpers/laws.py` (the iterate-table law), `tests/facade/test_session_topology.py`.

## Binding context

- The author's decisions on the alignment plan are in `.refactor/decisions.md` and plan section E; none gates B8 (Q9 records D11 and D14 as the class definition, unrelated here).
- `CLAUDE.md` invariants are physics: the iterate table is a partial bijection (`register_iterate` asserts it); the map is a bijection, so an image pair lies in at most one bridge's image chain; topology never calls the map; the chosen pip is remembered as an `Intersection` and restored via `registry.find`, a pip no longer a candidate warns once and is never replaced by the default (`classify_strong_pips()` included), and `topology()` raises (B8 gives that raise a remedy, `forget_strong_pip`, and stops the zone path from adopting a non-candidate).
- Test policy (`tests/conftest.py` docstring, 2026-10-05): rule-based laws, no case-specific pins; rubric 6 is a base-vs-branch fingerprint equality.
- Since B1 the registry is read as `N` through `resolved_ids`; B8's `unresolved_links` is a second such set (links the brackets cannot decide) and must stay disjoint from the table. Since B6a the walk leaves a class unresolved rather than guess; B8 does the same for the iterate table and the pip.
- The analyst stage is NOT skipped: the analyst verifies the two sites at `HEAD` (current line numbers), confirms by running that no collision and no multi-match fires on any law case, deep run or series step after B1, finds the existing iterate-table law and the inference test module the plan leaves open ("or the existing inference test module"), and reports what differs from the plan section. It does not re-plan the slice.
- The working tree may carry uncommitted documentation changes (`docs/`, `.refactor/`); the change branches from `HEAD` and never touches them.
