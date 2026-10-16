# Findings: B6a `walk-bounds-honest` (analyst, iteration 1)

*2026-10-12. Job `findings`. (The analyst's Write tool refused the artifact; the orchestrator saved this text verbatim.)*

**Dirty tree at `HEAD` (`76d0710`, main).** Uncommitted: `.refactor/decisions.md`, `.refactor/knowledge/{numerics,testing,topology}.md` (testing and topology now also carry this run's notes), `docs/alignment_plan_2026-10-09.md`, `docs/code_fidelity_audit_2026-10-09.md`, the staged move of `artifacts/request_origin_analyst.md` into `.refactor/runs/2026-10-12-oracle-laws/`, and untracked `.refactor/runs/2026-10-12-oracle-laws/*`, `.refactor/runs/2026-10-12-resolved-crossings/` and `artifacts/`. None of it is in B6a's four files. The branch starts from `HEAD` and must not pick any of it up.

**Kind and target.** This is a bug fix: three places where the walk claims a word it has not established. Targets: `DualWalk._nodes_meeting`, `DualWalk.lobe_image_markers` and `DualWalk.shortest_walks`; also `SymbolicDynamics._resolve_itinerary` and `SymbolicDynamics.is_reliable`. I did not re-plan the slice. Below is what holds, what does not, and the smallest corrections.

## Short answer

- The line numbers are unchanged. B0 and B1 did not touch either source file (`git log a3bd8f0..HEAD` is empty for both).
- **Two of the plan's factual claims are wrong at `HEAD`. Both were already wrong at `a3bd8f0`, so B1 did not cause them:**
  1. **(ii) does fire.** The "meets none of its" WARNING appears on the **k32 series at blasts 6 and 8**. Step 1 as written ("return `[]`") changes behaviour on that series step. It makes 10 classes unresolved and drops the closed alphabet at blast 8 from 45 to 42. That breaks an existing regression test.
  2. **Rubric 2's probe is already False on base.** `max_walks=1` on `heteroclinic` gives `is_reliable` False today, because two classes run out of budget with no walk. The probe that is True on base and False on the branch is **`max_walks=2`**.
- (i) and (iii) are as the plan says:
  - No search truncates at 64 or at 8 anywhere.
  - `max_walks=8` changes no word anywhere.
  - Members disagreeing on a marker never fires.

## The three sites at `HEAD` (verified by reading)

| Site | Plan said | At `HEAD` |
|---|---|---|
| (i) `is_reliable` ignores `truncated` | `SymbolicDynamics.py:839-848` | same: `:839-848`. It reads `itinerary is None`, `ambiguous`, `verified`, `contained`. |
| (i) where the status becomes the itinerary | — | `_resolve_itinerary` `SymbolicDynamics.py:1168-1192`. Only `status not in (unique, trivial, ambiguous)` is unresolved (`:1180`), and `search.truncated` is never read. |
| (i) where truncation is set | — | `DualWalk.py:1218` (`truncated = budget.truncated`); status chosen at `:1262-1264`, independent of `truncated`. The constrained no-walk case returns status `"truncated"` at `:1196-1207`, which is already unresolved. |
| (ii) span meets no node | `DualWalk.py:729-757` | same: `:729-757`. WARNING at `:747-754`, `return nodes` at `:756`; called at `:1104-1106`. |
| (iii) members disagree | `DualWalk.py:815-828` | same: the WARNING loop is `:815-828` inside `lobe_image_markers` `:778-829`; it returns `first` at `:829`. |
| walls built from decided walks | (not in plan) | `SymbolicDynamics.py:1408-1411` reads `search.status in ("unique", "trivial")` and `multiplicity == 1`. |

## What fires where (verified by running)

**How I ran it.** One probe per build, in a scratchpad script, using the frozen builders in `tests/cases.py`:
- every `BUILDERS` case except `orientation_reversing` (a known issue that raises): the 11 law cases, `p3_15`, `p3_16`, `p3_6_blasts`, `p3_16_12`, `inversion_period2_10`, `heteroclinic_17u`, `k28_deep_blast_sliver`, `nested_unblasted`, `p3_coarse_blasts` and `near_tangency`;
- `build_k28(blasts=8)`, which gives 4 classes at `HEAD`, as the request says;
- every step of all five `SERIES`, replaying the `blast_series` loop.

On each `session.topology()` I re-ran `symbolic_dynamics(dual_graph, classes, max_walks=m, letter_for=topo.letter_for)` for `m` = 64, 8 and 1. A wrapper on `shortest_walks` recorded `truncated`, and a log handler on the `DualWalk` and `SymbolicDynamics` loggers caught the three WARNING strings.

| Check | Result |
|---|---|
| Any `truncated` search at `max_walks=64` | none, on every case and series step |
| Any `truncated` search at `max_walks=8` | none; per-class itineraries at 8 equal those at 64 everywhere (the audit's "nothing changes" holds at `HEAD`) |
| "hold different markers" WARNING | never, at any budget |
| "meets none of its" WARNING | **k32_series blast 6** (2 searches, landing of `L#0`) and **blast 8** (8 searches, landings into `R#38` and `R#48`); nowhere else. It fires identically on the pre-B1 tree `a3bd8f0`. |
| `max_walks=1`, `heteroclinic` | `is_reliable` **False on `HEAD` and on `a3bd8f0`**. Classes `c` and `bh` are status `"truncated"` (constrained, no walk), so already unresolved. Classes `a` and `bn` are `unique` with `truncated=True`. `heteroclinic_17u` behaves the same. |
| `max_walks=2`, `heteroclinic` | `is_reliable` **True** on `HEAD`, with `bh` `unique` and `truncated=True`. This is the case B6a step 3 flips. From 3 up, nothing truncates. |
| Truncated walks vs the 64-budget words | `a`, `bn` (at 1) and `bh` (at 2) give the same itinerary as at 64. The claim was unproven, not wrong. |

### Site (ii) is a tolerance artifact, and the walk it hands on is right

All ten firings have the same cause.
- `_nodes_meeting` asks for overlap `> tol`, where `tol = SCALING_RTOL * |c|` at a scaled end and `SCALING_RTOL = 1e-2` (`Trellis.py:106`).
- At c ≈ 1.1–1.3 that is about 0.011–0.013. The landing spans are narrower than that (0.001–0.006 wide).
- So a span lying wholly inside one node's edge "meets" nothing. Examples:
  - blast 8, `R#48`: span `[1.25083, 1.25361]` inside its only node `[1.25083, 1.27024]`, overlap 0.0028 < tol 0.0125;
  - blast 6, `L#0`: span `[1.27072, 1.27379]` inside node `[1.27024, 1.27379]` of 22 nodes, overlap 0.0031 < tol 0.0127.
- `R#38` and `R#48` have one node each, so "any of them" is the one correct node.
- Blasts 7 and 9 give the same words, and there the end comes from the table, so the WARNING does not fire. All classes involved are `is_reliable` True, and `u`, `v`, `aj` and `as` are `verified` True.

**What step 1 as written does (simulated by monkeypatching `_nodes_meeting` to return `[]`).**
- k32 blast 6 becomes unreliable: `u` and `aa` are unresolved.
- k32 blast 8 becomes unreliable: `as`, `bf`, `bh`, `bj`, `br` and more are unresolved.
- At blast 8 the closed alphabet goes from 45 to 42, because unresolved classes no longer name their closure classes.
- The hole set is equal from blast 8 to 9 (checked with `blast_series`), and blast 9 still closes 45. So `test_closed_alphabet_constant_while_holes_constant[k32_series]` (`tests/regression/test_class_count_under_blasting.py:95`, run on every series, no known issue) would fail. The test is *believed* to fail; I did not run it under the patch. The quantities it compares were computed.
- This contradicts the plan's Behaviour line ("none on … any series step"). Rubric 3 does not list k32, so the rubric would not catch it.

**A two-tier rule keeps every word (simulated).** First keep the nodes that overlap by more than `tol`. If there are none, keep the nodes that overlap by more than 0. Only if there are still none, return `[]`.
- With this rule, the letter-free per-class itinerary hash is identical to base at every k32 step (0–10).
- Every firing above resolves to exactly one node.
- On every other case and series step the first tier already succeeds, so nothing else can change.

## Other differences from the plan section

1. **Rubric 2 must use `max_walks=2`, not 1.** At 1 the base is already False. At 2 the base is True, and the branch is False (simulated: `bh` becomes the only unresolved class).
2. **The wall gate (`SymbolicDynamics.py:1410`).** Step 3 puts the truncation check in `_resolve_itinerary` and leaves `search.status` as `"unique"`. A truncated, now-unresolved class would then still contribute its walk as a wall to later classes. That is again a claim from an unproven walk. On heteroclinic at `max_walks` 1 and 2, gating it changes no other class (simulated). It is a consistency point, not an observed effect.
3. **Step 2's mechanism is underspecified.**
   - `lobe_image_markers` returning `None` already means "no member lobe". `shortest_walks(wanted=None)` then runs the unconstrained search (`DualWalk.py:1152`).
   - A "flag the search records" can only be set after the search has run. If `None` is returned for a disagreement, the class gets an unconstrained walk, a word chosen under fewer constraints, unless the caller checks first.
   - The disagreement must be visible to `_resolve_itinerary` before it calls `shortest_walks`. Two ways: return a distinct value, or have `_resolve_itinerary` test the members itself. In either case the class becomes unresolved without walking. A `members_disagree` field on `WalkSearch` is then optional.
4. **"Never `unique`" should read "never decided".** In the constrained path a truncated search can also come back `trivial` (length 0, `DualWalk.py:1262`). Step 3 as worded (`search.truncated` → unresolved in `_resolve_itinerary`) covers that already. Only the plan's wording is narrower.
5. **A new status, not a rename.** Step 1 adds `"span_meets_no_node"` to the `WalkStatus` Literal (`DualWalk.py:552-561`) and to the `WalkSearch.status` docstring (`:579-583`). `_resolve_itinerary` (`:1180`) already treats any other status as unresolved (`"dual-graph walk span_meets_no_node"`), and `plotting._class_status` (`plotting.py:1607-1617`) shows the unresolved reason first. So `SymbolicDynamics.py` needs no change for step 1.
6. **One existing test pins the old behaviour.** `tests/unit/topology/test_dual_walk.py:333` `test_a_span_meeting_no_edge_falls_back_to_every_edge` asserts status `unique` with both start faces. It must be rewritten, not only supplemented; it is a natural home for rubric 4's span test. Its span `(5, 6)` against nodes `[0, 1]` and `[1, 2]` overlaps by nothing, so it is `[]` under the two-tier rule too.
   - `test_max_walks_cap_sets_truncated` (`:244`) and the truncated-status test (`:428`) read `shortest_walks` statuses. Step 3 leaves those alone.
   - Both unit modules pass at `HEAD`: 79 tests in 0.2 s.
7. **Docs describing the old semantics.** These are in the four files and would contradict the change: `DualWalk.py:91-99` (module Dev Notes: "the status still says unique …"), `:590-592` (`WalkSearch.truncated`), the `_nodes_meeting` docstring `:732-738`, `:808-810` (`lobe_image_markers` Returns), and `SymbolicDynamics.py:94-99`.
8. **Rubric 5 will work.** The phrase "may leave from any of them" is split across `DualWalk.py:748-749`, and `grep -n "may leave from any of them" src` matches `:749` today.

## Recommendation (the smallest change; it amends steps 1 and 3 and the rubric, it does not re-plan)

- **Step 1.** Make `_nodes_meeting` two-tier: first overlap `> tol`, then overlap `> 0`. Only when both are empty, return `[]`, and let `shortest_walks` return `WalkSearch(start, goal, "span_meets_no_node")` with one resolution-flag WARNING.
  - The second tier is not a fallback to "any of them". It keeps only nodes the span really overlaps. The tolerance exists to drop a neighbour touched by a scaled end's error; it should not reject the one edge a span narrower than its own error lies in.
  - This keeps "none on any law case, deep run or series step" true. As written, step 1 changes k32 blasts 6 and 8 and breaks a series test.
- **Step 2.** Have `_resolve_itinerary` detect the disagreement before walking: a distinct return from `lobe_image_markers`, or a direct check. Then mark the class unresolved with the plan's reason, and never call the unconstrained search for it.
- **Step 3, as planned.** Also gate the wall at `SymbolicDynamics.py:1410` on `cd.itinerary is not None` (one condition).
- **Rubric 2.** Use `max_walks=2` (base True, branch False).
- **Rubric 3.** Add `k32_series` steps 6 and 8, the only places site (ii) is exercised, or state the Behaviour claim as "every series step".
- **Steps 4 and 5, the not-doing list, the four files and the default `max_walks=64`:** unchanged.

## Verified vs believed

- *Verified (ran):*
  - none of the three fires at 64 or 8 except site (ii) on k32 6 and 8;
  - words at 8 equal words at 64 everywhere;
  - heteroclinic `is_reliable` at budgets 1 to 10, 16, 32 and 64 on `HEAD`, and at 1, 8 and 64 on `a3bd8f0`;
  - the k32 firings and their tolerance cause;
  - simulated "return `[]`": k32 6 and 8 unreliable, closed alphabet 45 → 42 at blast 8;
  - simulated two-tier: identical words at every k32 step;
  - simulated truncated → unresolved with and without the wall gate: no other class changes on heteroclinic;
  - k32 hole sets equal from 8 to 9;
  - the 79 unit tests pass.
- *Believed:*
  - that `test_closed_alphabet_constant_while_holes_constant[k32_series]` fails under step 1 as written (not run under the patch; the compared quantities were computed);
  - that "recorded logs" in the plan's Behaviour line did not include a k32 series run.

## Out of scope

- `SCALING_RTOL = 1e-2` is coarse next to landing spans one-tenth its size at c ≈ 1.2. Whether the landing span's own error bound (`land_element`) should be tighter is a numerics question, not B6a.
- `nested_unblasted`, `p3_coarse_blasts`, k28 at blast 0 and k32 at blast 0 are already `is_reliable` False, for reasons unrelated to the walk (for example, a singleton landing with no image chain).
- The k32 series now closes 45 classes at blasts 8–10. The knowledge file recorded a drop from 45 to 39 at blast 9 before B1, so B1 apparently removed it. I have not investigated this.

## Files

- Knowledge files updated:
  - `.refactor/knowledge/topology.md`: corrected the WP9b "never observed" note on `_nodes_meeting`, and added a B6a findings section.
  - `.refactor/knowledge/testing.md`: the pinned fallback test, the k32 exposure, the `max_walks=2` probe.
- Probe scripts are in `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b6a/`: `probe.py`, `sweep.py`, `cmp.py`, `k32.py`, `k32b.py`, `sim.py`, `walls.py`; per-case results are in `out/`.
