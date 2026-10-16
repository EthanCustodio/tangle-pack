# Findings: B3 `depth-of-N` — who fixes the depth `m` of `T_m(q0)`, at `HEAD` `9cd29d5`

*2026-10-14. (The analyst's Write tool refused the artifact; the orchestrator saved this text verbatim.)*

**The working tree is dirty.**
- Modified: `.refactor/decisions.md`, `.refactor/knowledge/{loom,numerics,testing,topology}.md` (this run added notes to all four), `docs/alignment_plan_2026-10-09.md`, `docs/code_fidelity_audit_2026-10-09.md`.
- The old request was moved into `.refactor/runs/2026-10-12-oracle-laws/`. Several `.refactor/runs/*` directories and `artifacts/` are untracked.
- Nothing under `src/` or `tests/` is dirty. The change branches from `HEAD`.

**Kind of request.** An integration question that ends in a small loom change. Targets:
- `loom/Blast.py`: the `blast_zone` child filter and `BlastResult`.
- `numerics/BridgeIterator.py`: `iterate_bridge`, where a child's ends are born.
- `topology/Trellis.py`: `_pip_orbit_cuts`, which decides `m`.
- `tests/cases.py`: `blast_series` and `SERIES`.

**The reading I chose.** "`forward_depth(pip)`" means the depth the cut actually uses: the composed `+1` chain (`Trellis.py:322-324`). The literal `IterateTable.forward_depth` is always 1 (section 1).

**Scripts and raw output** are in `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b3/`:
- probes: `depth.py`, `dbg2.py`, `p3img2.py`, `multi.py`, `slim.py`, plus their `*.txt` output;
- `base1029/`: an archive of `1029fcc`;
- `step1/`: an archive of `HEAD` with step 1 prototyped.

---

## 0. The answer in five points

1. **`min_separation` causes no plateau of any blast series at `HEAD`.**
   - Every series and recipe blasts with `num_iterations=1`.
   - At one iteration per call, the filter's verdict only shapes the frontier of the next iteration in the same call, and with one iteration there is none.
   - A dropped child stays un-iterated and inside the zone. The next call's frontier is rebuilt from un-iterated bridges (`Blast.py:285-287`), so it takes the child up again.
   - *Verified* (section 2).
2. **Every plateau at the recipe cutoff is mechanism (ii) or a new mechanism (iii).**
   - (ii), `f^{m+1}(q0)` detected but unresolved:
     - k28 at blast 8: id 245, ratio 2.991.
     - k32 at blast 9: id 523, ratio 1.477.
     - nested_outer at blast 9: id 447, ratio 1.775.
     - inversion at blast 10: id 238, ratio 1.571.
   - (iii), `f^{m+1}(q0)` **never detected**: the image polyline's last chord passes on the wrong side of the stable anchor. This happens on p3 at blast 15 and on heteroclinic at blast 6, which is the same point.
   - nested_inner does not stall through 12 blasts.
3. **Plan step 1 alone removes no series plateau.**
   - I prototyped it in a scratch copy. The `m` series and the registry sizes on k28, k32 and nested_outer are identical to `HEAD`.
   - It only removes a plateau inside one multi-iteration call. One `blast_zone(num_iterations=11)` on nested_outer ends at `m = 7` at `HEAD` and at `m = 9` with step 1.
4. **The plan's k32 numbers predate B1, and their cause was misattributed.**
   - At `HEAD` the k32 series reads `[1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9]` over blasts 0..10.
   - The plan's `[1, …, 10, 10]` is the pre-B1 base; I reproduced it on `1029fcc`.
   - That plateau was **not** `min_separation` either. The child ending at `f^10(q0)`, (522, 523), was discarded as **exterior** by the zone filter, which N-a keeps.
5. **Step 1 breaks a test outside the plan's four files.**
   - The test is `tests/unit/loom/test_blast_proximity_guard.py::test_min_separation_drops_close_bridges`.
   - Its fixture's only too-close child, (46, 47), ends at `f^4(q0)`.
   - Under the prototype the test fails with `too_close == 0`.

---

## 1. The code at `HEAD` (verified by reading; current line numbers)

**`Blast.py`, `blast_zone` (line 185):**
- `:285-287`: the frontier is rebuilt on every call from `workbench.uniiterated_bridges`, filtered by `fp_ok` and `in_zone`.
- The child loop (`:385-420`):
  - `:388-390`: a child already seen is skipped.
  - `:392-395`: a child that is not (un-iterated, own fixed point, in zone) is **discarded as exterior**. This is the zone filter, and it runs *before* the separation test.
  - `:397-412`: the `min_separation` test. A child too near the existing curve (`ref_tree`) or a kept sibling goes to `too_close`, is added to `seen`, and the loop continues.
  - `:414-420`: kept, and added to the reference clouds.
- `:459`: `frontier = kept`. This is the only consumer of the keep/drop verdict besides the `BlastStep` / `BlastResult` counters.
- `:463`: `workbench.infer_iterates()` runs at the end of every call.
- The plan's `:397-412` numbering is still exact.

**`BlastResult` (`:133-178`):**
- Fields: `zone`, `fixed_points`, `num_iterations_requested`, `steps`, `interior_bridges_by_iteration`, `completed_iterations`, `stop_reason`, `points_created`, `too_close`, `already_known`.
- The last two are missing from the docstring.
- There is no `pip_depth`.

**`BridgeIterator.iterate_bridge` (`numerics/BridgeIterator.py:264`) — where a child's end crossings are born:**
- `:321-324`: if the image re-traces existing curve, the existing bridges are returned and nothing is registered.
- `:363`: detection, via `populate_intersections_for_manifold(iterated)`.
- `:369-371`: `ix.id = registry.add(ix)` registers **every** detection, resolved or not.
- `:378-379`: `cut_at` keeps only the detections in `registry.resolved_ids(threshold)` (B1).
- `:380-394`: `create_bridges(cut_at)`, or a single partial piece.
- `:398-405`: `infer_iterates(sources = the parent's two ends + cut_at)`. This is where `f^m → f^{m+1}` is linked. `IterateInference.py:116-122` links only crossings of N, so an unresolved `f^{m+1}` is never linked.
- `:413-414`: `_register_bridge(child)`, then `bridge.iterated = True`.

So `f^{m+1}(q0)` is born, and linked, when its **parent** (the bridge ending at `f^m`) is iterated. That happens before the child filter runs. The filter only decides whether the child ending at `f^{m+1}` is iterated **later in the same call**.

**`m`:**
- `Trellis._pip_orbit_cuts` (`Trellis.py:320-333`) does `m = iterate_table.forward_depth(pip)`, then `while self.iterate(pip, m + 1): m += 1`.
- `IterateTable.forward_depth` (`:81-83`) reads direct entries only, and only `n = 1` links are ever registered (`IterateInference.py:147,202`). So it returns 1 for any linked crossing. *Verified*: `fd=1` at every step of every series.
- `IterateTable.forward_chain` (`:89-96`) has the same limitation.

---

## 2. Measured: the pip depth `m` per blast at the recipe cutoff

**How it was measured** (*verified*; `depth.py`):
- Builds as E12 does (`e11.build`), keeping each recipe's own `area_cutoff`.
- Each step is one `blast_zone(num_iterations=1, fixed_point=[fp], min_separation=SEP, diagnose=False)`.
- It blasts the largest zone, or the smallest for nested_inner and heteroclinic.
- `m` is the composed chain from `session.chosen_strong_pip(fp)`.

**Two cross-checks:**
- Adding the per-step `compute_intersections(preserve_ids=True)` + `rebuild_bridges()` of `cases.blast_series` (`_recomputed`) gives identical `m` and registry sizes on k28, k32 and nested_outer.
- The stalls line up with the `n_trellis_crossings` plateaus in `stress/numerics/results/e13_both_intervals.jsonl`:

| series | `n_trellis_crossings` plateau | `m` stalls at |
|---|---|---|
| k28 | 143 from blast 7 | blast 8 |
| k32 | 300 from 8 | 9 |
| p3 | 119 from 14 | 15 |
| nested_outer | 276 from 8 | 9 |
| inversion | 167 at 10 | 10 |
| heteroclinic | 197 from 5 | 6 |

**`m` per blast:**

| series (cutoff, separation, blasts) | `m` over blasts 0..N | first stall | registry / unresolved at the stall |
|---|---|---|---|
| k28 (1e-7, 1e-5, 13) | 1,2,3,4,5,6,7,8, **8**,8,8,8,8,8 | blast 8 | 244 / 5 |
| k32 (1e-7, 1e-5, 13) | 1,2,3,4,5,6,7,8,9, **9**,9,9,9,9 | blast 9 | 522 / 1 |
| p3 (1e-7, 1e-5, 20) | 14,15,…,28, **28**,…,28 | blast 15 | 133 / 1 |
| inversion (1e-5, 1e-5, 10) | 5,6,…,14, **14** | blast 10 | 247 / 2 |
| nested_outer (1e-7, 1e-4, 12) | 1,…,9, **9**,9,9,9 | blast 9 | 448 / 1 |
| nested_inner (1e-4, 12) | 14,15,…,26 | none | 105 / 0 at blast 12 |
| heteroclinic (1e-4, 8) | 26,…,31, **31**,31,31 | blast 6 | 216 / 1 |

- `m` rises by exactly one per blast until the first stall, then stays flat.
- On every series, the stall is the first blast whose registry holds an unresolved detection on the pip orbit's branches. On k28, ids {75, 76, 77} from blast 6 are not orbit points.

### Attribution

At every first stall, the bridge ending at `f^m(q0)` **was in the frontier and was iterated**. That rules out mechanism (i).

| series | stall | mechanism | evidence |
|---|---|---|---|
| k28 | blast 8, `f^8 = 147` | **(ii)** | Id **245**, ratio **2.991**. Its unstable cdist is 278 699 285 against 147's cdist × β_u = 278 699 207 (2.8e-7 apart), on the same branches. The other unresolved ids (75/76/77 at u = 8.36e6; 148 at u = 2.26e8, ratio 2.572) are not orbit points. |
| k32 | blast 9, `f^9 = 301` | **(ii)** | Id **523**, ratio **1.477**. Unstable cdist 2.40928036e9 against 2.40927986e9 (2.1e-7). Its stable cdist is 4.25e-7, against 3.44e-7 from one β_s step and 1.40e-7 from `image_cdist`. |
| nested_outer | blast 9, `f^9 = 275` | **(ii)** | Id **447**, ratio **1.775**. Unstable cdist 616 944 647 against 616 944 486 (2.6e-7). |
| inversion | blast 10, `f^14 = 190` | **(ii)** | Id **238**, ratio **1.571**. Unstable cdist 21 948.97 against 21 948.85 (5.7e-6). Id 203 (ratio 2.168, u = 19 108) is not the orbit point. |
| p3 | blast 15, `f^28 = 118` | **(iii)** | `f^29` was never detected; see the notes below the table. |
| heteroclinic | blast 6, `f^31 = 194` | **(iii)** | Same point and same chord as p3 (its p3 orbit). Unresolved id 197 is the same off-orbit fold as p3's 119. |
| nested_inner | none through 12 | — | — |

**Why p3's `f^29` was never detected:**
- The true image under the real map, `f(f^28 q0)`, is (−6.97e-6, 1.000004). That is 8.05e-6 from the anchor (0, 1) of `s(0,0)`.
- The image of bridge (117, 118) ends with a 0.0138-long chord, (−0.00793, 0.99542) → (0.00582, 1.00336).
- That chord passes about 1e-6 *below* the anchor. `s(0,0)` leaves the anchor up and to the left, so the chord misses it.
- Hand-iterating the bridge registers only ids 119-121, all on `s(1,0)`.
- The step's single unresolved id, 119 (ratio 1.499, u = 108 186), is a different fold. `f^29` is expected at u ≈ 207 917 from one full return of `f^26`.
- Source: `p3img2.py`.

**Mechanism (i) never fires under single-iteration blasting.** In every case below, the pip-orbit child was too close, yet `m` still rose by one on the next blast:
- nested_inner, blasts 10-12: (83, 84), (91, 92), (103, 104).
- heteroclinic, blast 5: (193, 194).
- nested_outer on the pre-B1 base, blasts 6-8: (114, 115), (174, 175), (274, 275).

**Mechanism (i) is real inside one multi-iteration call** (`multi.py`):
- One `blast_zone(nested_outer, num_iterations=11, min_separation=1e-4)` drops (114, 115), which ends at `f^7`, at iteration 6. The call ends at `m = 7`, where the series reaches 9.
- Single calls on k32 (10 iterations), k28 (8) and p3 (20) end at 9, 8 and 28, the same stalls as their series. None has a pip-orbit drop.

### The pre-B1 base (`1029fcc`, run from an archive)

*Verified:*
- k32 reads `m` = 1,…,9,**10,10** over blasts 0..10, the plan's numbers. At blast 9 the child **(522, 523)** lands in `discarded_exterior`.
- nested_outer reads 1,…,9,**10,10,10**. **(446, 447)** is exterior at blast 9.

So the pre-B1 plateau was the **zone filter**. Call it mechanism **(iv)**; it sits at `Blast.py:392-395`, before the separation test, and step 1 does not reach it.

*Believed:* once (ii) is gone, k32 would stall next on (iv) at blast 10.

---

## 3. What plan step 1 alone achieves

*Verified with a prototype* in a scratch copy of `HEAD`. A child whose end lies on the composed forward chain of a selected fixed point's chosen pip skips the separation test and logs one INFO line.

- **Series: nothing changes.**
  - k28 to 8, k32 to 10 and nested_outer to 11 give identical `m` and registry sizes per blast (k32 522 / 898; nested_outer 448 / 736 / 1218; k28 244).
  - Law-case fingerprints: unchanged (*believed*; no builder uses more than one iteration).
- **It removes:** plateaus inside one multi-iteration call. nested_outer with a single 11-iteration call goes from `m = 7` to 9, then stops on (ii).
  - Callers with multi-iteration blasts: `scripts/henon_blast_period_3.py`, E12's two-iteration cells, and the tests `test_blast_proximity_guard`, `test_crossing_resolution`, `test_arrangement`, `test_blast_monotonicity` and `test_blast_no_overlap`.
- **It cannot remove:**
  - any series plateau: (ii) on k28, k32, nested_outer and inversion; (iii) on p3 and heteroclinic;
  - the pre-B1 plateau (iv).
- **Collateral.** `tests/unit/loom` plus the multi-iteration blast tests give 71 passed and 1 failed under the prototype. The failure is `test_min_separation_drops_close_bridges`: at `HEAD` its fixture's only too-close child is (46, 47), and 47 = `f^4(q0)`.

---

## 4. What differs from the plan section (B3)

1. **The diagnosis is wrong, both at `HEAD` and pre-B1.** The plan says `min_separation` drops the child carrying `f^{m+1}` and never feeds it again.
   - At one iteration per call, a dropped child is always re-fed by the next call.
   - The pre-B1 stall was the zone filter.
   - The post-B1 stalls are (ii) and (iii).
2. **The behaviour section** says the k32 series reads `[1, …, 9, 10, 10]` and that the branch reaches 11.
   - At `HEAD` k32 reads `[1, …, 9, 9, 9]`, and step 1 leaves it exactly so.
   - No N, hole or class change at blast 10; the same holds for nested_outer.
3. **Rubric 2 cannot be met with step 1.** It expects `[1, …, 10, 11]` on the branch. The branch equals the base on every series.
4. **Step 2's `forward_depth(pip)` is 1 for every pip.** `pip_depth` must compose the `+1` links as `_pip_orbit_cuts` does. `forward_chain` does not compose either.
5. **Step 3's law fails on base and branch alike.** It asks for "+1 on every step whose frontier contained the bridge ending at `f^m`".
   - It fails on `k28_series` at blast 8, `k32_series` at blast 9 and `nested_outer_series` at blast 9. In each case that bridge was in the frontier and was iterated, but `m` did not move.
   - `p3_series` (9 blasts) and `inversion_series` (4 blasts) end before their stalls.
   - As written, the law needs `KNOWN_ISSUES` entries for (ii) on three of the five series, or a condition that reads the outcome.
6. **Rubric 4's unit test can only be observed inside a multi-iteration call.** At `num_iterations=1`, keeping or dropping the child leaves the same state.
7. **The file list misses `tests/unit/loom/test_blast_proximity_guard.py`**, which step 1 breaks.
8. **"Law cases: none" holds**: no builder uses multi-iteration blasts.

---

## 5. Recommendation: the smallest change

N-a's step 1 is correct, but on today's recipes it is almost inert: it never decides `m` on any series. A B3 that states this honestly:

1. **Keep step 1 as decided (N-a).** Test whether a child is pip-bearing on the **composed** chain: repeat `registry.iterate_table[x, 1]` from `session.chosen_strong_pip(fp).id` for each fixed point in `fixed_points`. Place the test after the exterior check, as the plan says. About 10 lines.
2. **`BlastResult.pip_depth` = the composed depth.** Add `test_blast_proximity_guard.py` to the file list. Either give it a fixture or separation whose drops are not pip-bearing, or exclude pip-bearing kept children from its distance assertion.
3. **The rubric-4 unit test uses a multi-iteration call.** For example nested_outer with `num_iterations ≥ 7`: at `HEAD`, (114, 115) is dropped at iteration 6.
4. **Reword the step-3 law**, or carry it as a measured law with `KNOWN_ISSUES` entries for (ii).
5. **Restate the behaviour section and rubric 2:** B3 leaves the series unchanged.

**For the planner to put to the author, not decided here:**
- Whether (ii) belongs to B3. B1's remedy is "refine"; the alternative is to place the point by `image_cdist`, as PK-6 does for a missing one.
- (iii), which has no detection to refine. A PK-6-style placement of `f^{m+1}` would cover both (ii) and (iii).
- (iv), the zone filter. It is next on k32 once (ii) is gone, and N-a keeps it.

---

## 6. Out of scope (noticed, not recommended)

- `BlastResult.too_close` and `already_known` (`Blast.py:168-169`) are undocumented.
- The names `IterateTable.forward_depth` / `forward_chain` (`:81-96`) are misleading under a table that holds only `+1` links.
- At s ≈ 1e-7, the stable cdist of k32's (ii) detection disagrees with both scalings (4.25e-7, 3.44e-7, 1.40e-7). A placement at that depth would carry a large stable error.
- The `min_separation` docstring (`:211-216`) says "sibling", but the code tests all prior curve as well (`:399-409`). Plan step 4 rewrites this docstring anyway.

**What was verified and what is believed:**
- *Verified* by running code: every `m` series; each stall's id and ratio; the p3 / heteroclinic geometry under the real map; the pre-B1 base numbers; the step-1 prototype's series and test results; the fixture's pip-bearing too-close child.
- *Believed*:
  - law-case fingerprints are unchanged under step 1 (from reading the code and grepping the builders);
  - id 503 in the k32 10-iteration call is `f^10`;
  - k32 stalls next on (iv) once (ii) is gone.

**Knowledge files updated:** `.refactor/knowledge/loom.md` (the blast-filter line numbers and a new B3 analyst section), `numerics.md` (`forward_depth` / `forward_chain` read direct entries only), `topology.md` (the stall attribution), `testing.md` (series stalls; the proximity-guard fixture's pip-bearing child).
