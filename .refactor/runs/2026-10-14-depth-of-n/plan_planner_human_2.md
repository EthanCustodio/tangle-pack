# Plan: B3 `depth-of-n`. Who decides the depth `m` of `T_m(q0)`

*2026-10-14. Base `9cd29d5` (main). Findings: `artifacts/findings_analyst_planner.t1.md`.*

**What the plan expected and what we measured.** Section B3 assumed the blast's `min_separation` drops the child carrying `f^{m+1}(q0)` and never feeds it again, and that this freezes `m` while N keeps deepening. That is not what happens.

- **Single-iteration blasts.** Every recipe and every series blasts one iteration per call. A dropped child stays un-iterated inside the zone, and the next call rebuilds its frontier from the un-iterated bridges, so the child is iterated anyway.
- **The real plateaus.** The analyst measured every series at `HEAD`. The depth of `T(q0)` is decided by whether `f^{m+1}(q0)` is detected and resolved, not by the child filter. Each plateau is one of two things:
  - **(ii) Detected but unresolved:**
    - k28 at blast 8: id 245, ratio 2.991.
    - k32 at blast 9: id 523, ratio 1.477.
    - nested_outer at blast 9: id 447, ratio 1.775.
    - inversion at blast 10: id 238, ratio 1.571.
  - **(iii) Never detected:** p3 at blast 15, and heteroclinic at blast 6, which is the same point. The image polyline's last chord passes about 1e-6 on the wrong side of the stable anchor.
- **The pre-B1 k32 plateau** quoted in B3 (`[…, 10, 10]`) was a third mechanism, (iv): the zone filter (`Blast.py:392-395`) discarding child (522, 523) as exterior. N-a keeps that filter.
- **What step 1 does.** The prototype changes no series and no law case. It only matters inside one multi-iteration call: a single 11-iteration blast of nested_outer ends at `m = 7` today and at `m = 9` with step 1.
- **Corrections to the section's other details:**
  - It breaks `test_blast_proximity_guard.py::test_min_separation_drops_close_bridges`, whose only too-close child ends at `f^4(q0)`.
  - `IterateTable.forward_depth` is always 1, so the depth has to be composed link by link.
  - Its step-3 law fails on base and branch alike.
  - Its rubric 2 (`[1, …, 10, 11]`) cannot be met.

## 1. The answer

`m` is the length of the registered `+1` chain from the chosen pip, composed as `Trellis._pip_orbit_cuts` does. Blasting deepens it by exactly one per blast until `f^{m+1}(q0)` is born unresolved (ii) or not detected at all (iii). After that `m` stays flat on every later blast. `min_separation` never decides `m` on any series. It can only delay the pip orbit inside a multi-iteration call. N-a's step 1 is correct, but on today's recipes it is inert. The stalls belong to (ii) and (iii), which are numerical (detection and resolution), not part of the blast's bookkeeping.

## 2. Recommendation

Do the slice as decided but restated, option **(a)** in question 1: N-a's step 1, a composed `pip_depth`, and the series law weakened to "never shallower". Behaviour is unchanged on every series and law case. This implements Q2 as you decided it. It also makes the stalls visible as data (`BlastResult.pip_depth`, `SeriesStep.pip_depth`) without claiming they are fixed. Whether (ii), (iii) or (iv) get a slice is a separate decision, listed under "Not doing".

The rest of this plan is written for (a). If you pick (b) or (c), say so and the plan is redone or the run closes.

## 3. Scope

Slug `depth-of-n`: branch `refactor/depth-of-n`, worktree `../tangle-pack-depth-of-n` off `9cd29d5`.

The change touches five files, about 120 changed lines:

1. `src/tanglepack/loom/Blast.py`
2. `tests/cases.py`
3. `tests/regression/test_class_count_under_blasting.py`
4. `tests/unit/loom/test_blast_proximity_guard.py`
5. `tests/unit/loom/test_blast_pip_orbit.py` (new)

## 4. Steps

1. **Exempt pip-bearing children (`Blast.py`, child loop).**
   - After the exterior check (`:392-395`) and before the separation test (`:397`), a child whose end crossing lies on a selected fixed point's pip orbit is kept regardless of separation.
   - Each kept-by-exemption child gets one INFO line naming it and the orbit point.
   - The orbit is the composed `+1` chain, recomputed once per iteration after the children are produced. `iterate_bridge` links `f^m -> f^{m+1}` before the filter runs.
   ```python
   def pip_orbit() -> set[int]:
       """The chosen pips of the selected fixed points and every registered forward iterate."""
       table = workbench.intersection_registry.iterate_table
       orbit = set()
       for fp in fixed_points:
           x = getattr(session.chosen_strong_pip(fp), "id", None)
           while x is not None and x not in orbit:
               orbit.add(x)
               x = table[x, 1]
       return orbit
   ...
   pip_bearing = {child.id[0], child.id[1]} & orbit   # partial children have id None: never pip-bearing
   if min_separation is not None and not pip_bearing: ...   # the existing test
   ```
2. **`BlastResult.pip_depth: dict[FixedPoint, int]`.**
   - Read after the final `infer_iterates()` (`:463`): the length of the composed chain minus one, per selected fixed point that has a chosen pip.
   - Documented, together with the undocumented `too_close` and `already_known` fields of the same dataclass.
   - It never calls `IterateTable.forward_depth`.
3. **Docstring of `blast_zone` (`min_separation`, `:211-216`).**
   - `min_separation` is a cost guard against children that re-trace existing or sibling curve.
   - A child ending on the pip orbit is exempt.
   - It never decides the depth of `T(q0)`. That is decided by whether `f^{m+1}(q0)` is detected and resolved (A.4).
4. **Series data and law.**
   - `cases.SeriesStep.pip_depth: int` is composed on the test side from `session.chosen_strong_pip(fp)` through `trellis.iterate(pip, n)`, so the law does not trust the library's own count.
   - New regression law `test_pip_orbit_never_shallower(name)`: `pip_depth` is non-decreasing over the series. It pins no number.
5. **Proximity-guard test.** `test_min_separation_drops_close_bridges` is restated to the rule as it now reads:
   - every kept child keeps its distance unless it is pip-bearing;
   - every dropped child was too close;
   - `too_close > 0` is asserted with the pip forgotten (`session.forget_strong_pip(fp1)`), so the guard still demonstrably fires.
6. **New unit tests (`test_blast_pip_orbit.py`).** They use the `henon_p3_session` fixture (`fp1`, 3 iterations, separation 0.02):
   - (a) with the pip chosen, the child whose end lies on the pip orbit is in `kept_interior` and not in `discarded_too_close`;
   - (b) after `forget_strong_pip(fp1)`, the same child (by bridge-id end pair) is in `discarded_too_close`;
   - (c) `result.pip_depth[fp1]` equals the depth composed by `trellis.iterate` after the blast.

   The child is found by rule (an end on the composed orbit), never by pinned ids.

## 5. Behaviour

- **Unchanged:** every law case, and every step of all five `SERIES` (holes, classes, N, bridges, pip depth).
- **One intended change, inside a single multi-iteration call (N-a):** a pip-bearing child the guard used to drop is now iterated in the same call. Measured:
  - the `henon_p3_session` 3-iteration blast keeps the child ending at `f^4(q0)`;
  - one `blast_zone(nested_outer, num_iterations=11, min_separation=1e-4)` reaches `m = 9` instead of 7.
- **Other callers that may change:** `scripts/henon_blast_period_3.py` and E12's two-iteration cells can differ the same way. Their numbers go in the report.

## 6. Not doing

- **(ii) Unresolved `f^{m+1}(q0)`.** B1's remedy is "refine", and nothing here changes it. Placing the point by `image_cdist` is option (b) in question 1.
- **(iii) Undetected `f^{m+1}(q0)`** (p3 blast 15, heteroclinic blast 6). There is no detection to refine. Only a placement or a finer image near the anchor reaches it.
- **(iv) The zone filter** discarding a pip-bearing child as exterior (pre-B1 k32 blast 9). After (ii), it is next to stall k32 (believed, not verified).
- **N-b:** no `depth=` parameter, as decided.
- **Renaming** `IterateTable.forward_depth` / `forward_chain`, which are misleading under a table of `+1` links only (noted, untouched).
- **Recipes:** no change to their `num_iterations=1`.

## 7. Rubric

The reviewer checks each item by running or reading, never from the report.

1. `MPLBACKEND=Agg env/bin/python -m pytest -q -rxX` passes on the branch, and only `cases.KNOWN_ISSUES` appear as xfails.
2. Run the branch's tests against the base source (`9cd29d5`):
   - these pass: `test_pip_orbit_never_shallower` on all five series, and the restated `test_min_separation_drops_close_bridges`;
   - these fail: unit tests 6(a) and 6(c);
   - 6(b) passes on both.
3. `git diff --stat 9cd29d5..` lists exactly the five files of section 3 and at most 300 changed lines.
4. The fingerprint of every `cases.LAW_CASES` case is equal on base and branch.
5. The series are unchanged: every `SeriesStep` field (holes, classes, closed, references, letters, crossings, recomputed) of all five series is equal on base and branch.
6. The measured `pip_depth` per step is the same on base and branch. These lists are measured by a script, not pinned in a test; the k28, k32 and nested_outer lists end in a plateau:

   | series | blasts | `pip_depth` |
   |---|---|---|
   | k28 | 0..8 | `[1..8, 8]` |
   | k32 | 0..10 | `[1..9, 9, 9]` |
   | nested_outer | 0..11 | `[1..9, 9, 9, 9]` |
   | p3 | 0..9 | `[14..23]` |
   | inversion | 0..4 | `[5..9]` |
7. One `blast_zone(nested_outer zone, num_iterations=11, fixed_point=[fp], min_separation=1e-4)` ends with composed depth 7 on base and `result.pip_depth[fp] == 9` on the branch.
8. In `Blast.py`:
   - the pip-orbit exemption sits after the exterior check and before the separation test;
   - `forward_depth` is not called in the new code;
   - one INFO line is logged per exempted child.
9. The `blast_zone` docstring says `min_separation` never decides the depth of `T(q0)` and names the pip-orbit exemption. The `BlastResult` docstring documents `pip_depth`, `too_close` and `already_known`.
10. No behaviour change beyond section 5: the only new log lines are the exemption INFO lines.

## 8. Questions

**1. What is B3?** Default: **(a)**.

- **(a) N-a as decided, restated (this plan).**
  - Contents:
    - step 1 on the composed chain;
    - a composed `pip_depth`;
    - the series law "non-decreasing";
    - rubric 2 restated as "the series are unchanged";
    - the unit test on a multi-iteration call;
    - the proximity-guard test added to the file list and adjusted.
  - Five files; behaviour is none on the series.
  - What you get: Q2 is implemented, the stalls are visible as data, and the multi-iteration delay is gone.
  - What you do not get: any deeper `T(q0)` on any series.
- **(b) Widen B3: place an undetected or unresolved `f^{m+1}(q0)` by `image_cdist`, as PK-6 places a missing orbit point, so the cut keeps its depth.**
  - Covers (ii) and (iii).
  - Changes every deep series from its stall onward: k28 from 8, k32 and nested_outer from 9, inversion from 10, p3 from 15, heteroclinic from 6. All of them would be re-derived.
  - It reverses B1's "refine" decision.
  - It needs a rule for how far to place past the last registered point, which drifts toward N-b, which you did not choose.
  - At k32's depth the stable cdist disagrees with both scalings (4.25e-7 detected against 1.40e-7 placed), so the placed cut carries a large stable error.
  - A bigger slice; better as its own run after (a).
- **(c) `answered`.**
  - The findings answer the question and the code is unchanged.
  - N-a stays unimplemented, so the multi-iteration delay remains (nested_outer 11 iterations: `m = 7`).
  - The B3 section in `docs/` keeps its wrong diagnosis until B12.
  - Whether (ii), (iii) or (iv) get a slice is decided later.

**2. (Only under (a).) The series law's form.** Default: **"non-decreasing" only**, as in step 4. The alternative keeps B3's sharper law, "+1 on every step whose frontier held the bridge ending at `f^m`", carried as a measured law with three `KNOWN_ISSUES` entries for (ii): k28, k32 and nested_outer. It would fail strictly until some later slice fixes (ii). Choose it if you want the stalls on the xfail list rather than only in the data.
