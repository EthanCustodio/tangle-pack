# Outcome: B3 `depth-of-n`, option (a)

*2026-10-14. Branch `refactor/depth-of-n` @ `7d3a294` (worktree `../tangle-pack-depth-of-n`), base `9cd29d5`. The review passed in round 1. There is no debrief verdict and no figure report. Not merged: the merge is yours.*

> **Your working tree is dirty** (`main`): `CLAUDE.md`, `docs/`, `.refactor/` and the untracked run folders and figure sweeps. The branch was cut from `HEAD` and touches none of them. I have also edited `.refactor/decisions.md` and `.refactor/knowledge/loom.md` there.

**What happened.** You asked who decides the depth `m` of `T_m(q0)`. The plan section said `min_separation` drops the child that carries `f^{m+1}(q0)`, and that this freezes `m` while N keeps deepening. The measurements disagree.

- **The child filter is not what sets the depth.** Every recipe and every series blasts one iteration per call. A dropped child stays un-iterated inside the zone, so the next call iterates it anyway. The depth of `T(q0)` is set by whether `f^{m+1}(q0)` is detected and resolved.
- **Every series stall is an orbit point that is either unresolved or never detected:**
  - **(ii) Detected but unresolved:** k28 at blast 8, k32 at blast 9, nested_outer at blast 9, inversion at blast 10.
  - **(iii) Never detected:** p3 at blast 15, and heteroclinic at blast 6.
- **The pre-B1 k32 stall the plan quoted** (`[…, 10, 10]`) came from a third cause, **(iv) the zone filter**: it threw the pip-bearing child away as lying outside the zone.

So B3 did N-a as you decided it, restated:
1. A child ending on the chosen pip's orbit is exempt from `min_separation`.
2. `BlastResult.pip_depth` reports the depth as data.
3. A series law says the depth never gets shallower.
4. Nothing changes on any series.

## 1. The answer to your question now

- **What sets `m`:** the length of the registered `+1` chain from the chosen pip. It is composed link by link, exactly as `Trellis._pip_orbit_cuts` cuts `T(q0)`.
- **How it moves:** each blast deepens it by one, until `f^{m+1}(q0)` is born unresolved (ii) or not detected at all (iii). After that it stays flat.
- **`min_separation` no longer decides it**, even inside a call of several iterations (`Blast.py:237-247`, `:440-447`).
- **The depth is now visible** at two places:
  - `BlastResult.pip_depth` (`Blast.py:181-199`, `:515-517`);
  - `SeriesStep.pip_depth`, composed on the test side (`tests/cases.py:1215`).
- **The stalls themselves are untouched.** Section 7 covers them.

## 2. What changed (five files, the plan's; 237+/50−, 204+/17− ignoring whitespace)

**Most important hunk: `Blast.py`'s child loop.** Before, each parent's children were sorted as soon as that parent had been mapped:
```python
            produced.extend(children)
            for child in children:
                if _seen_key(child) in seen: ...
                if not (not child.iterated and fp_ok(child) and in_zone(child)): ... discarded
                if min_separation is not None:
                    ... too_near_existing / too_near_sibling -> too_close; continue
                seen.add(...); kept.append(child) ...
```
After: every parent of the iteration is mapped first. Then the orbit is read once, and only then are the children sorted (`Blast.py:420-471`):
```python
        orbit = {x: (fp, n) for fp, chain in _pip_orbits(session, fixed_points).items()
                 for n, x in enumerate(chain)}
        for child in produced:
            if _seen_key(child) in seen: ...                                  # unchanged
            if not (not child.iterated and fp_ok(child) and in_zone(child)): ...  # unchanged
            pip_end = next((x for x in child.id or () if x in orbit), None)
            if min_separation is not None and pip_end is not None:
                logger.info("blast: child %s ends at f^%d(q0) of %s (crossing %d), kept regardless of separation", ...)
            elif min_separation is not None:
                ...                                                       # the old test, verbatim
            seen.add(...); kept.append(child) ...
```

**`_pip_orbits` (`Blast.py:109-128`)** returns `[q0, f(q0), …, f^m(q0)]` for each selected fixed point.
- It finds the pip with `registry.find(chosen_strong_pip(fp))`.
- It then follows `iterate_table[x, 1]` link by link, and never calls `forward_depth`.

**Tests:**
- `tests/cases.py`: `SeriesStep.pip_depth` and `_pip_depth`.
- `test_class_count_under_blasting.py`: the law `test_pip_orbit_never_shallower`, over five series.
- `test_blast_proximity_guard.py`: the guard test restated and parametrized `pip` / `no_pip`.
- New `test_blast_pip_orbit.py`, four tests: kept with a pip, dropped without one, `pip_depth` equals the composed orbit, and the exemption skips only the separation test.

## 3. What got shorter or clearer

Nothing got shorter. `Blast.py` grew from 468 to 522 lines; most of the raw diff is the child loop moving out one indentation level. What got clearer:
- **The `min_separation` docstring** now says what the guard is: a cost guard against curve already in the tangle or a kept sibling. The old text said "sibling" only, which was wrong. It also says the guard never decides the depth of `T(q0)`.
- **`BlastResult`** now documents `too_close` and `already_known`, which had no docstring, as well as `pip_depth`.

## 4. Every behaviour change

1. **None on any law case or series.**
   - The test engineer measured all 11 `LAW_CASES` fingerprints equal on base and branch, and the coder did too, with an independent canonicalisation.
   - Every `SeriesStep` field of all five series is byte-identical on base and branch.
   - `pip_depth` reads:

   | series | `pip_depth` |
   |---|---|
   | k28 | `[1..8, 8]` |
   | k32 | `[1..9, 9, 9]` |
   | nested_outer | `[1..9, 9, 9, 9]` |
   | p3 | `[14..23]` |
   | inversion | `[5..9]` |

2. **Inside one call of several iterations, a pip-bearing child is kept and iterated in the same call.** One `blast_zone(nested_outer, num_iterations=11, min_separation=1e-4)` now ends at depth **9 instead of 7**: too_close 78 instead of 76, registry 498 instead of 486. E12's nested_outer two-iteration cell differs from call 3 onward (the coder's replay).
3. **That same call now fires two existing WARNINGs that base does not fire:**
   - `361 same-stability polyline crossing(s) of the new image bridges in the tangle region`. The exempted `f^9` child (258, 259) is among them.
   - `1 of 494 detections are not crossings of N`, ratio 1.775. That detection is the (ii) stall point.

   The exemption brings back near-coincident curve that the guard used to stop. The coder's report left this out; the test engineer caught it.
4. **One new INFO line per exempted child**, whenever `min_separation` is set.
   - The rubric-7 call logs 8 of them, `f^2 … f^9`.
   - Single-iteration blasts usually log one each, so a series run with the recipes' `1e-5` now logs one INFO per blast.
5. **The child filter now runs after every bridge of the iteration is mapped.** The reviewer checked by instrumentation that this changes nothing on its own: 0 children altered across 611 + 47 + 67 + 79 children, and identical kept and too-close lists with the pip forgotten. The test engineer's break 1 (exemption removed) reproduces base's logs exactly.

## 5. What was verified, and how

**Verified by me:**
- **Suite on the branch:** `MPLBACKEND=Agg env/bin/python -m pytest -q -rxX` gives `1193 passed, 38 skipped, 5 deselected, 1 xfailed` (245 s). The only xfail is the `KNOWN_ISSUES` `orientation_reversing` entry.
- **Diff scope:** `git diff --stat 9cd29d5..` lists exactly the plan's five files. I read the whole diff.
- **The 11-iteration nested_outer call**, my own script (`scratchpad/aud_b3/multi.py`):
  - branch: depth 1 → 9, `pip_depth {A: 9}`, 8 exemption INFO lines, and the two WARNINGs above;
  - base: depth 1 → 7, no WARNING.
- **`scripts/henon_blast_period_3.py`** raises `ValueError: Intersection 12 is not a strong-pip candidate. Candidates: [31, 32, 33]` at `:95` on both base and branch. The script was last touched in `5f072ee`.

**Verified by the reviewer and the test engineer** (their runs; I did not repeat them):
- The branch's tests run against the base source:
  - fail: 6(a) `test_pip_bearing_child_is_kept` and 6(c) `pip_depth` (AttributeError);
  - pass: everything else, including the series law on all five series.
- Law-case fingerprints and all five series are equal on base and branch.
- Six break-it checks, each red where it should be. Moving the exemption before the exterior check was caught only by the test the test engineer added.

## 6. What did not go to plan

- **The diagnosis.** The plan section's premise, that `min_separation` freezes `m` on the series, was wrong (summary above). As a result:
  - its rubric 2 (`[1, …, 10, 11]` on k32) could not be met and was restated as "series unchanged";
  - its sharper law ("+1 on every step whose frontier held the `f^m` bridge") was replaced by "non-decreasing", your default for question 2;
  - `IterateTable.forward_depth` turned out to be always 1, so the depth is composed link by link.
- **Five files instead of four.** `test_blast_proximity_guard.py` had to change, because its only too-close child ends at `f^4(q0)` and is now exempt.
- **Deviations the reviewer accepted:**
  - the pip is found by `registry.find(chosen)` rather than the stored `.id`, which can go stale;
  - there is no cycle guard in the orbit walk (a strong pip is not periodic, and the table is an asserted partial bijection);
  - 6(b) finds its child by rule on its own build rather than by carrying ids over from (a);
  - the child loop was moved after the parents are mapped.
- **The test engineer's added test,** `test_exemption_skips_only_the_separation_test` (`7d3a294`). Before it, the exemption's position, after the exterior check, was guarded only by reading the code.
- **Plan section 5 asked for the numbers of `scripts/henon_blast_period_3.py`.** There are none: it fails before blasting, on base too.
- **Minor, not graded:** the test-side orbit composition is written three times (`cases._pip_depth` and the two unit-test helpers).

## 7. Where I am least sure this is what you meant

1. **The stalls are now visible as data, but nothing fixes them, so no series gets a deeper `T(q0)`.** They are:
   - (ii) unresolved `f^{m+1}(q0)` on k28, k32, nested_outer and inversion;
   - (iii) undetected on p3 and heteroclinic;
   - (iv) the zone filter, believed to be the next stall on k32 once (ii) is fixed.

   Your options:
   - **Refine:** B1's decision, and the remedy for (ii). Re-grow the parent's image at a finer `area_cutoff` near the orbit point.
   - **Place:** put an unresolved or missing `f^{m+1}(q0)` by `image_cdist`, as PK-6 already places a missing orbit point, so the cut keeps its depth.
     - This would reverse B1's "refine" for the pip orbit.
     - It needs a rule for how far past the last registered point to place, which drifts toward N-b.
     - At k32's depth the placed stable cdist is off by about 3x (1.40e-7 placed against 4.25e-7 detected).
   - **For (iii), a finer image near the anchor:** there is no detection to refine. The image's last chord passes about 1e-6 on the wrong side of the stable anchor.

   Each is its own slice. The sharper law, as `KNOWN_ISSUES` xfails, would put (ii) on the xfail list if you want it tracked there.
2. **The exemption overrides a resolution guard.** Inside a call of several iterations it keeps a child that runs within `min_separation` of existing curve, and that call then raises 361 same-stability WARNINGs (item 3 in section 4). This is N-a as written: the guard never decides the depth. But you may have meant the guard to keep winning when the curve really is under-resolved. If so, the exemption should go to the WARNING path rather than to the kept path.
3. **Stale documentation and a broken script:**
   - The B3 section of `docs/alignment_plan_2026-10-09.md` still carries the wrong diagnosis: `:149-154` (decision N), `:903-943` (B3: `forward_depth`, `[…, 10, 11]`, the four files) and `:1166-1167` ("the pip orbit never stalls on k32"). That is B12's job, or yours now.
   - `scripts/henon_blast_period_3.py` pins pip id 12, which is no longer a candidate, so the script is dead on `main` regardless of B3.
