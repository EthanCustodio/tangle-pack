# Outcome: B4b `holes-are-combinatorial`

*Branch `refactor/holes-are-combinatorial` @ `ba23a1f` (worktree `../tangle-pack-holes-are-combinatorial`), base `f166573`. Review passed in round 1. Not merged; merging is your call. Your main checkout has uncommitted `docs/`, `.refactor/` and `artifacts/` changes. The branch does not touch them.*

**What a hole is now.** A hole is its combinatorics and nothing more: `(origin, iterate, bounding_ids, bridge_side, openings)`, plus `near_intersection_id` and `pair`. `punch_holes` punches the reference pairs. Each one's side is read off its crossing signs, the cross-product rule. `propagate_reference_holes` then walks each reference bridge's crossing ids backward through `image_cdist` and finds the bridge that contains them. It carries the side along with `side_after`, which flips it only under an orientation-reversing map. That is the whole algorithm. `StablePartition.py` no longer imports numpy, holds no coordinates and calls no map. Drawing is a separate concern handled by one function, `plotting.hole_point(trellis, hole)`. It computes a point on demand. For a propagated hole it carries the reference hole's point back with `map_inv`, which a view is allowed to do. The algorithm does not do that.

## The answer now

Yes, a hole is its five fields, and where to draw it is a plotting matter. With today's builds there is nothing fancier left in the algorithm. *Verified:* the rubric grep on `StablePartition.py` (for `epsilon`, `coords`, `np.`, `map_inv`, `dynamical_system`, `_map_backward`, `_region_key`, `punched_regions` and `fallbacks`) prints nothing. The file's imports are now logging, bisect, itertools, typing and local modules only.

## What changed

The most important hunk is the propagation loop in `propagate_reference_holes` (`StablePartition.py:419-472`).

Before (`f166573`), with the fallback start (about 40 lines), `seen` and `punched_regions` elided:
```python
start_iterate = start.iterate or 0
carried = np.asarray(start.hole.coords, dtype=np.float64)
for step in count(1):
    iterate = start_iterate - step
    span = tuple(sorted(trellis.image_cdist(end_id, -step, "unstable", branch=start_key)[1] for end_id in start_ids))
    carried = _map_backward(trellis, carried)          # the map, inside the algorithm
    containing = _containing_bridge(trellis, span, advance(start_key, -step), orders)
    ...
    key = (origin, iterate)
    if key not in seen:
        hole = _punch_in_bridge(..., carried=carried, span=span)   # placed the hole; could return None
        if hole is not None:
            region = _region_key(origin, hole)
            if region not in punched_regions: ...
```
After:
```python
for step in count(1):
    ends = [trellis.image_cdist(end_id, -step, "unstable", branch=start_key) for end_id in start_ids]
    span = tuple(sorted(c for _, c, _ in ends))
    branch_key = advance(start_key, -step)
    containing = _containing_bridge(trellis, span, branch_key, orders)
    ...  # NUM-5 WARNING: a SCALED end within SCALING_RTOL of a crossing, once per branch
    new_holes.append(_punch_in_bridge(
        trellis, containing, bridge_side=side_after(start_side, -step, preserving),
        iterate=-step, origin=origin))
```

The other changes:
- `punch_holes(trellis, pairs=None)` asserts `p.iterate >= 0` (`:304`). Its else-branch placement is gone.
- The anchor block of `_hole_openings` (`:695-705`) names the branch with `stable_key_at`. It then applies the one identification: the anchorward half of `s(o,b)` is the outward half of `s(o,1-b)`.
- `Hole.coords` is removed, and so is `Trellis.punch_holes(epsilon=)`.
- `plotting.hole_point` is new (`plotting.py:312`). The four placement helpers moved into `plotting.py` unchanged.
- `measured_sides` and three laws now read their point from `hole_point`.

## What got shorter

| | before B4a (`95ea780^`) | B4a (`f166573`) | B4b (`ba23a1f`) |
|---|---|---|---|
| `StablePartition.py` | 1,947 | 1,354 | **1,091** |
| `topology/diagnostics.py` | — | 712 | 713 |
| `topology/plotting.py` | 3,207 | 3,207 | 3,408 |

In B4b, `src` is net −64 lines (+321/−385) and tests are +103. Most of B4b's reduction is code that moved, not code that was deleted. About 200 lines of placement went to `plotting.py`. The real deletions are the dead backward-pair path, the two dedups, the fallback and `_map_backward`, about 110 lines. Across the three files together, the line count is flat (5,154 → 5,212). What changed is that the algorithm file now holds only the algorithm.

## Behaviour changes

- **Topology output: none.** *Verified* by the test engineer's base-vs-branch probe, which I read in full. All 11 law cases and `p3_coarse_blasts` match on the five fields, the law fingerprint, the full words report and `measured_sides`.
- **Figures: identical.** *Verified* by my own independent probe. I ran `plot_holes` on base and on the branch and compared the rendered scatter offsets and iterate labels on `k28_two_blasts` (4 holes), `inversion_period2` (10) and `heteroclinic` (73): max |d| = 0.0, labels equal. Each run imported its own worktree's `src`.
- **Logs:**
  - The 14 "Orbit … cannot propagate" WARNINGs on the non-law `near_tangency` builder are gone. Its "Skipped 14 pair(s)" WARNING remains.
  - There is one new WARNING kind (NUM-5, a scaled span end beside a crossing). It fires once on `k28_deep_blast_sliver` and on no law case. These rest on the coder's probe output; I did not re-run them.
- **`epsilon` now also scales the propagated nudge** in `hole_point`. Base hard-coded 0.05 there. The output is the same at the default, and no caller passes another value.
- **`hole_point` raises a bare `StopIteration`** for a propagated hole whose origin has no reference hole in `trellis.holes` (`plotting.py:343`). `Trellis.punch_holes` can never produce that case.
- **Edge cases no build reaches, found reading the diff:**
  - A bridge with no nodes used to have its hole skipped (`arc_point is None`, `_punch_in_bridge` → None). It now gets the hole, and `hole_point` would then fail on it.
  - A hand-built pair with `iterate=None` now raises `TypeError` in the assertion.
  - The B4b analyst's probe found neither condition on any builder or series step.
- **API:** `Hole(coords=)` and `punch_holes(epsilon=)` are gone. No script, stress file or notebook used either. *Verified* by grep: every remaining `.coords` is `Intersection.coords`.

## What was verified, and how

- **Full suite on the branch** (I ran it: `MPLBACKEND=Agg env/bin/python -m pytest -q -rxX` in the worktree): **1178 passed, 38 skipped, 5 deselected, 1 xfailed** in 205 s. The only xfail is `orientation_reversing` (C3: INV-11, a KNOWN_ISSUE). Base was 1173 with the same xfail. The difference is the 3 new coder tests plus 2 from the test engineer. The run used the branch's `src`: the new tests import `hole_point`, which base does not have.
- **Diff scope** (`git diff --stat f166573`): exactly the nine plan files plus `TangleSession.py`. The `TangleSession.py` change is one docstring line that drops `epsilon`.
- **The whole diff read.** The algorithm side matches the plan. The anchor rewrite is algebraically the base rule. The NUM-5 check reads `orders[branch_key]` after `_containing_bridge` has filled it. `diagnostics` importing `plotting` adds no new import, because `Trellis` already imports `plotting` at module level and pyplot is loaded by `import tanglepack`.
- **Break-it checks** (test engineer, all caught):
  - deleting the assertion;
  - the WARNING watching table ends (6 law cases would log it);
  - the WARNING never firing;
  - one `map_inv` step dropped (the 1e-12 comparison fails on every law case);
  - the no-map fallback broken (caught only by the test engineer's new test);
  - the anchor identification removed or forced (15 and 7 law failures).

## What did not go to plan

- **Step ordering.** The carried point, `_map_backward` and `epsilon` were deleted in step 5, not in steps 1 and 2. `Hole.coords` stays a required field until step 5, so every intermediate commit needed a placement to stay green. Between commits 4 and 5 the four helpers existed twice; the coder checked commit 4 against base and found it bit-identical. The end state is the plan's. There are five commits, in plan order.
- **One extra test from the test engineer**, `test_a_hole_is_drawn_off_its_bridge[map|no_map]`. Breaking the no-map fallback passed the whole existing suite, so this test now guards it.
- **The diff is larger than the plan's estimate:** +442/−403 against the estimated +120/−300, because of the verbatim moves and the new tests.

## Where I am least sure this is what you meant

1. **`plotting.py` now calls the map** (`hole_point` → `map_inv`), and `measured_sides` calls it through `hole_point`. "The algorithm calls no map" holds for `StablePartition.py`. However, `punch_holes(diagnose=True)` is the default, so `topology()` still runs `map_inv` on the flag path, as it did after B4a. If you meant "no map anywhere in `topology/` on the default path", that is not true yet.
2. **The law oracles and `measured_sides` read their point from `plotting.hole_point`.** The geometric laws (`check_backward_holes_keep_bridge_side`, `check_direct_hole_side_is_coordinate_side`) and the resolution flag now share code with the figure. A placement bug would therefore move the oracle and the drawing together. The only guard against that is the base-vs-branch 0.0 comparison, which is a probe, not a committed test.
3. **Loose ends, outside the slice:**
   - The I1/I2 deciding cores still live in `diagnostics` (B4a's open point).
   - `CLAUDE.md` is stale. Its holes bullet (`:35`) and `Hole` row (`:139`) still describe `measured_side` and `side_source` as hole state, and it says nothing of `hole_point` or the NUM-5 WARNING. It never names `coords` or `epsilon`, which I checked by grep. The "span falls back to all nodes" text is B6a's. All of this is B12's job.
   - Two comments still say propagation "seeds already-punched orbits and regions" from `pair.hole`: `Trellis.py:1250-1251` and `StablePartition.py:306-308`. The reviewer flagged only the first.

Updated: `.refactor/decisions.md` (B4b outcome entry), `.refactor/knowledge/topology.md` (Map rewritten for the branch, B4b audit section), `.refactor/knowledge/testing.md` (B4b tests).
