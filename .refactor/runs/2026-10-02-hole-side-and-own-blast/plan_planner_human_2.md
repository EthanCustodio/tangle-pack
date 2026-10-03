# Plan: blasting nested tangles through to symbolic dynamics

*Planner, 2026-10-02. Base commit `2315204` (`main`). Source: `artifacts/findings_analyst_planner.t1.md`.*

## 1. The answer

No, not once blasting goes deep. Period 3 alone keeps the same words from 0 to 5 blasts (`a -> b -> c -> a u^-1 w^-1`, no image bridges). At 6 blasts it raises in `check_holes_share_bridge_side`. The bug is in how holes are measured, not in the invariant. `_punch_in_bridge` decides a propagated hole's side from the nearest vertex of the whole containing bridge. For origin (12, 33) at iterate -3, that vertex sits at unstable cdist 887, on a different fold from the hole's own sub-arc (707.7–745.0). It is the same failure as the open "16 unstable steps" case: same spans, and iterates -3/-5/-6 fail in both. A direct hole has a second version of this mistake: it is measured one node past its crossing. Period 1 alone runs from 0 to 8 blasts and its words do not change from blast 2 on.

Nested fails for a second reason. `build_nested` passes both fixed points to every blast, so an outer blast also blasts the period-3 bridges, including period 3's exterior. Because of that, the result depends on blast order, a numerics merge crash is caught and logged as "skipped", and 6 or more outer blasts collapse the period-3 pseudoneighbors from 34 to 2. With the side fix, and each zone blasting only its own fixed point, all 16 nested cells give exactly the words of each tangle run alone, in either order.

## 2. Recommendation

Change it, with three small fixes: two in `StablePartition` (the hole side) and one line in the nested recipe.

## 3. Scope

Slug: **`hole-side-and-own-blast`**. Branch `refactor/hole-side-and-own-blast`, worktree `../tangle-pack-hole-side-and-own-blast`.

Files:

1. `src/tanglepack/topology/StablePartition.py`
2. `src/tanglepack/examples/henon_cases.py`
3. `tests/test_stable_partition_period3.py`
4. `tests/test_higher_period_cartoon.py`
5. `CLAUDE.md`: one sentence (the "16 steps trips `check_holes_share_bridge_side` (open)" claim). See Question 2.

About 80 changed lines.

## 4. Steps

**A. A propagated hole's side is read against its image sub-arc.** `_bridge_side_of` gains an optional `span`. `_punch_in_bridge` is the only caller that passes it.

```python
def _bridge_side_of(trellis, bridge, point, span=None):
    poly = _oriented_bridge_polyline(trellis, bridge)
    if poly is None:
        return None
    if span is not None:
        # the image sub-arc only: a far fold of the same bridge can be nearer (2026-10-02)
        cdists = sorted(node.cdist for node in bridge.get_point_array(return_nodes=True))
        lo = max(bisect_left(cdists, span[0]) - 1, 0)
        poly = poly[lo : bisect_right(cdists, span[1]) + 1]
    sign = _arc_side_of(poly, np.asarray(point, dtype=np.float64))
    ...
```

`_punch_in_bridge` (`:1526`) calls `_bridge_side_of(trellis, bridge, carried, span)`. The openings still come from the side, as they do now.

**A'. A direct hole takes its lobe's side by construction, not by measuring.** At `:349`, the side comes from the crossing signs, in the same spirit as `row_of_end`:

```python
first, second = sorted((near, far), key=lambda x: x.unstable_cdist)
outward = 1 if second.stable_cdist > first.stable_cdist else -1
bridge_side = "left" if first.crossing_sign * outward > 0 else "right"
```

The openings stay `inward=True` and are unchanged. Replace the "estimated side" comment with the reason for the rule.

**B. Each zone blasts only its own fixed point.** In `build_nested` (`henon_cases.py:253-258`), change the argument to `fixed_point=[zone.fixed_point]`. Reword the docstring (`:229-230`): the order no longer matters, and the reason is that a zone's blast touches only its own tangle. Update the Dev Notes line about 16 steps the same way.

**Tests** (each under ~2 s):

- `test_stable_partition_period3.py`:
  - `build_period3(blasts=6)` reaches `symbolic_dynamics()` and `check_holes_share_bridge_side` passes;
  - the same for `build_period3(unstable_steps=15)`.
  - Both check that it runs, not that `is_reliable` holds (see Not doing 1).
- `test_higher_period_cartoon.py`:
  - the inner-tangle words of `build_nested(outer_blasts=2, inner_blasts=4)` equal `build_period3(blasts=4)`'s;
  - nested 4 + 4 gives the same words and bridge count inner-first and outer-first. The outer-first run is driven through `session.blast_zone(zone, 1, fixed_point=[zone.fixed_point])`.
  - Words are compared through each class's element pair in short names, without the fixed-point letter, plus the token direction. Letters differ between sessions, so they are not compared directly.

## 5. Behaviour

- **A:** in today's passing runs, no hole changes side. In the failing runs, exactly the 3 propagated holes at iterates -3/-5/-6 flip to agree with iterate 0.
- **A':** `Hole.bridge_side` of a direct hole can change only where the nearest-vertex reading was wrong (the +1/+2 lobes at 16 steps and at 7 or more blasts). Openings and partitions do not change.
- **B:** at the defaults (outer 2, inner 0), the nested words are unchanged and the period-3 bridges go from 31 to 27. Nested runs with outer blasts no longer depend on blast order and no longer log "blast: skipping bridge".
- **Newly runs:** period 3 at 6 or more blasts, and at 15 steps, now runs instead of raising. It is unreliable (one unreachable class, one virtual `new1`).

## 6. Not doing

1. Period 3 at 6 or more blasts after the fix is unreliable: a new reference orbit with 20 backward holes, a chain of about 18 active classes, an unreachable `t`, a virtual `new1`, and "image bridge abuts no hole".
2. The numerics crash in `ManifoldMachine.iterate_manifold:295 -> merge_manifolds:438` (non-contiguous `old_iterated_points`). B avoids it but does not fix it.
3. `_hole_openings` (`:1211-1217`) uses the same whole-polyline nearest vertex. It is latent; no failure has been seen.
4. The session forgets the strong-pip choice on every generation bump (`_repin` after each blast). Keeping the pip as an `Intersection` would remove the repin calls.
5. `Blast.py` tidy-up: the duplicate `_resolve_fixed_points`, the cKDTree import inside the loop, the dead scipy fallback, and `BlastResult.fixed_point` typed `Optional[FixedPoint]` while it holds a list.
6. The `build_nested` slug is the same whatever the blast counts, so figure runs overwrite each other.
7. The "inert class has a non-empty word" WARNING on every period-3 run (could be INFO for inert -> inert).
8. Blast speed (`Tangle._insert_segment` and rtree queries are about 50%, refinement about 30%).

## 7. Rubric

1. `MPLBACKEND=Agg env/bin/python -m pytest -q` passes on the branch with at least 782 passed, 1 skipped, plus the new tests.
2. The pinning tests (`test_stable_partition_period3.py`, `test_higher_period_cartoon.py`, `test_loom_blast_restore.py`, as they are at base) pass at `2315204` and on the branch.
3. The new test for `build_period3(blasts=6)` **fails at base** with `AssertionError` from `check_holes_share_bridge_side` naming origin `(41, 42)`, and passes on the branch. Run the test file from the branch against the base worktree.
4. The new test for `build_period3(unstable_steps=15)` fails at base (origin `(12, 33)`) and passes on the branch.
5. `git diff --stat 2315204` lists only the five files in section 3.
6. `_bridge_side_of` has a `span` parameter and `_punch_in_bridge` passes it. `grep -n "_bridge_side_of(" src/tanglepack/topology/StablePartition.py` shows no other caller passing a span.
7. In `punch_holes`, the direct-hole `bridge_side` no longer calls `_bridge_side_of`. It is computed from `crossing_sign` and the stable cdists.
8. `henon_cases.py` has no `fixed_point=fixed_points` left in `build_nested`, and the docstring no longer says the inner blasts run first "so the outer blasts do not consume its bridges".
9. No behaviour change beyond section 5:
   - the k=10 and k=2.8 fixture words in `CLAUDE.md` are unchanged (the existing tests);
   - `build_period3(blasts=4)` words are still `a -> b -> c -> a u^-1 w^-1`;
   - `build_nested()` at the defaults gives the same words as at base, and `blast_sizes`/bridge counts move only as stated (period-3 bridges 31 -> 27).
10. `build_nested(outer_blasts=4, inner_blasts=4)` logs no `blast: skipping bridge` warning (check with caplog or by running it).

## 8. Questions

1. **Direct-hole rule.** Should it be combinatorial (crossing signs, no geometry, in the style of `row_of_end`) or geometric (the sign of the lobe polygon's area)? Both agreed on all 22 direct holes tested, but only the geometric form was run through the full suite. *Default: combinatorial. If the suite disagrees anywhere, the coder stops and reports rather than switching.*
2. **CLAUDE.md.** Change "16 steps trips `check_holes_share_bridge_side` (open)" to say the bug is fixed, and that 15 steps and 6 or more blasts run but are unreliable (open)? *Default: yes, that one sentence only.*
