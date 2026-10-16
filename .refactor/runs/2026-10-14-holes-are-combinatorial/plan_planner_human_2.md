# Plan — B4b `holes-are-combinatorial`

*Base `main` @ `f166573` (B4a merged). First plan for this slice; there is no earlier verdict to answer.*

After this slice a hole is only its combinatorics: `(origin, iterate, bounding_ids, bridge_side,
openings)`. `Hole.coords` no longer exists, and `StablePartition.py` no longer places anything in
the plane or calls the map. Drawing a hole becomes the job of one plotting function,
`plotting.hole_point(trellis, hole)`. It computes the point when a figure, a law oracle or the
measured-side diagnostic needs one, using the same arithmetic as today, so every figure stays the
same. The direct punch of backward pairs, which no build has ever produced, goes, together with its
two dedups and its fallback. A new assertion makes that deletion falsifiable.

## 1. The answer

Yes, a hole can be its five fields. The analyst checked this at `HEAD` with a probe that ran today's
propagation next to an instrumented copy: there were 0 mismatches on all 22 builders and every step
of the five series.

- **No backward pairs.** No recorded pseudoneighbor pair has `iterate < 0` on any of those runs. The
  dedups `seen` and `punched_regions` rejected 0 holes, and every origin has exactly one reference
  hole. That makes the else-branch, both dedups and the fallback dead code.
- **Placement in plotting is exact.** Moving placement into plotting and carrying the reference point
  back with `map_inv` reproduces today's `Hole.coords` **bit for bit**: a deviation of 0.0 on every
  hole of every builder.
- **The no-map alternative moves holes.** Placing propagated holes without the map moves them by up
  to 0.328, and puts them on the other side of their image sub-arc on all 11 law cases.

## 2. Recommendation

Change it. Slice B4b as decided, with the analyst's five corrections folded in:

- **(A) The new WARNING looks at scaled span ends only.** An end that comes from the iterate table
  sits exactly on the container's own end crossing. Watching every end would fire the WARNING on 6
  of the 11 law cases. So the `from_table` flag is **used**, not dropped.
- **(B) The line target is ≤ 1,100.** The original ≤ 700 cannot be reached inside the slice.
- **(C) The span comment is reworded.** It currently says "never in the table", which is false for
  one end.
- **(D) `_near_far` stays.** It orders `bounding_ids` and `near_intersection_id`, which are both hole
  fields. `_region_key` is deleted along with the dedups.
- **(E) Four more files join the list.** One debug log line that reads `coords` changes too.

## 3. Scope

**Slug** `holes-are-combinatorial`. The branch is `refactor/holes-are-combinatorial` in
`../tangle-pack-holes-are-combinatorial`, made from `f166573`. The base worktree is
`../tangle-pack-holes-are-combinatorial-base`, detached at `f166573`. Your uncommitted `docs/` and
`.refactor/` changes are not touched.

**Files (nine), plus one docstring line:**

1. `src/tanglepack/topology/StablePartition.py`: the punch, the propagation and the openings; the
   placement helpers leave this file.
2. `src/tanglepack/topology/plotting.py`: new `hole_point`, which takes in the placement helpers;
   `plot_holes` calls it.
3. `src/tanglepack/topology/TopologyResults.py`: `Hole.coords` is removed.
4. `src/tanglepack/topology/Trellis.py`: `punch_holes(epsilon=)` is removed.
5. `src/tanglepack/topology/diagnostics.py`: `measured_sides` reads the reference point from
   `hole_point` (about 3 lines, plus a stale docstring note about `_map_backward` "until B4b").
6. `tests/helpers/laws.py`: three laws read `hole.coords` (`:1213`, `:1254`, `:1280`); they read
   `hole_point` instead.
7. `tests/unit/topology/test_stable_partition.py`: drops `Hole(coords=)` at `:82` and `:392`;
   `:355-373` reads `hole_point`; adds the two new unit tests.
8. `tests/unit/topology/test_topology_invariants.py`: `Hole(coords=)` at `:68`.
9. `tests/unit/topology/test_minimal_trellis.py`: `Hole(coords=)` at `:71`.

The **docstring line** is `src/tanglepack/loom/TangleSession.py:1321`: `punch_holes` names
`epsilon`, and that mention goes. It is the only change in that file.

The diff is roughly −300 / +120 lines, which is within a targeted change. Most of it is helpers
moving from `StablePartition.py` to `plotting.py`.

## 4. Steps

Commit each step on its own, in this order.

**Step 1. Reference pairs only (`punch_holes`).**
- New signature: `punch_holes(trellis, pairs=None)`.
- It punches the pairs with `iterate == 0` and skips forward pairs, as today (D15). This is the same
  set as today, because `is_reference` holds exactly when `iterate == 0`.
- It asserts that no pair has `iterate < 0`:
  ```python
  assert all(p.iterate >= 0 for p in pairs), "holes propagate backward only; no backward pair is recorded"
  ```
- Deleted: the else-branch placement; `seen`, `punched_regions`, `_region_key`, the `fallbacks`
  path, and its "cannot propagate" WARNING and "has no punched bridge" INFO.
- `_near_far` stays, because it orders `bounding_ids`.

**Step 2. Propagation by ids alone (`propagate_reference_holes`, `_punch_in_bridge`).**
- The span keeps the table flag:
  ```python
  # The span's interior is never a crossing of N (it lies past the pip cut); one end
  # may be a registered crossing, which is then the container's own end.
  ends = [trellis.image_cdist(end_id, -step, "unstable", branch=start_key) for end_id in start_ids]
  span = tuple(sorted(c for _, c, _ in ends))
  ```
- After `_containing_bridge` returns, it logs one WARNING per branch when a **scaled** end
  (`from_table` False) lies within `SCALING_RTOL * |x|` of a crossing `x` on the container's branch.
  In that case the container could be its neighbour (NUM-5).
- `_map_backward` and the `carried` point are deleted. `_punch_in_bridge` loses its `carried` and
  `span` parameters and builds the `Hole` from ids, side and openings alone.
- The `:1090` debug log ("Hole at %s has no openings") names `hole.bounding_ids` instead of
  `coords`.

**Step 3. Anchor openings (`_hole_openings`).**
- The rule is the D6 formula plus one identification, stated in the comment:
  - name the branch with `trellis.stable_key_at(id, unstable_key)`;
  - when D6 says "anchorward", open the **outward** half of `s(o, 1-b)`, because *the two stable
    branches of an inversion point are one smooth curve through it, so the anchorward half of
    `s(o,b)` at the anchor is the outward half of `s(o,1-b)`*;
  - open nothing when that branch is not computed (a point seeded on one eigendirection).
- This is today's code reworded, so the counts stay as they are. The inversion cases open on
  `s(o,1-b)` at 2/4, 2/2, 2/2 and 4/8 anchor ends. On non-inversion cases, about half the anchor ends
  open nothing.

**Step 4. `plotting.hole_point(trellis, hole, epsilon=0.05)`.**
- Moved here from `StablePartition`: `_step_into`, `_stable_arc_midpoint`, `_nearest_arc_point` and
  `_image_arc_midpoint`, verbatim.
- **Reference hole:** take the chord midpoint of `bounding_ids` (already in near/far order) and
  nudge it toward the bridge.
- **Propagated hole:**
  1. compute the reference point of the origin's reference hole in `trellis.holes`;
  2. carry it back `-iterate` steps with `trellis.dynamical_system.map_inv`;
  3. take the span from the reference ids by `image_cdist`;
  4. nudge the sub-arc midpoint toward the carried point.
- **No map** (`dynamical_system is None`): fall back to the sub-arc midpoint nudged toward the chord,
  exactly as today.
- `plot_holes` calls `hole_point` and reads only `origin`, `iterate` and the point.
- `diagnostics.measured_sides` and the laws `check_holes_are_classified`,
  `check_backward_holes_keep_bridge_side` and `check_direct_hole_side_is_coordinate_side` call
  `hole_point` where they read `coords`. `diagnostics` imports `plotting` inside the function; the
  analyst checked that this creates no import cycle.

**Step 5. Remove the field and the parameter.**
- Delete `Hole.coords` and `Trellis.punch_holes(epsilon=)`, and the `epsilon` line in the
  `TangleSession` docstring.
- Remove `Hole(coords=)` from the three test files.
- Drop the `numpy` / `NDArray` imports in `StablePartition.py`, which nothing uses any more.

## 5. Behaviour

- **Topology output:** none. Every hole's five fields, every partition, class, word and law-case
  fingerprint stays the same.
- **Figures:** identical under placement (a), the default in question 1.
- **Logs:**
  - one new WARNING kind (a scaled span end near a neighbour crossing). It fires on no law case. It
    does fire on `k28_deep_blast_sliver` and `k32_series`, where NUM-5 lives; no test reads it there.
  - the 14 "Orbit … cannot propagate" WARNINGs on the non-law builder `near_tangency` disappear with
    the fallback. That builder created no holes from them, and its "Skipped 14 pair(s)" WARNING still
    names the same pairs.

## 6. Not doing

- "One element at a bridge end" (B7).
- The `_containing_bridge` bisect itself, and any change to which bridge gets a hole.
- Moving the I1/I2 deciding cores back from `diagnostics` (B4a's open point).
- The `direction = -1` loop in `Pseudoneighbor.extend_pseudoneighbor_trajectories`. It is the
  upstream of the dead path; the new assertion guards it from downstream.
- Today's placement puts some propagated holes on the wrong side of their sub-arc (heteroclinic
  13/69, `inversion_period2` 4/8, `p3_coarse_blasts` 43/80). This affects drawing only, and (a) keeps
  it as it is.
- `StablePartition.py` down to 700 lines. That would mean moving its Dev Notes or the partition code.
- `CLAUDE.md`, which still describes `coords`, `measured_side` and `side_source` (B12).

## 7. Rubric

1. The suite passes on the branch (`MPLBACKEND=Agg env/bin/python -m pytest -q -rxX`), and only
   `cases.KNOWN_ISSUES` appear as xfails.
2. Every test that touches holes (the law tier, `test_stable_partition`, `test_topology_invariants`,
   `test_minimal_trellis`, `test_production_checks`, `tests/plotting`) passes on the base worktree
   and on the branch. The exceptions are the two new unit tests (items 6 and 7) and the tests edited
   only to drop `Hole(coords=)`; each of those must pass on the branch.
3. `git diff --stat f166573..refactor/holes-are-combinatorial` touches only the nine files and
   `TangleSession.py`. The `TangleSession.py` diff is one docstring line.
4. No behaviour change beyond section 5. The law-case fingerprints, holes compared by their five
   fields, are equal on base and branch for all 11 law cases.
5. `StablePartition.py` is at most 1,100 lines (it is 1,354 at base), and this grep prints nothing:
   `grep -n "epsilon\|coords\|_step_into\|np\.\|map_inv\|dynamical_system\|_map_backward\|_region_key\|punched_regions\|fallbacks" src/tanglepack/topology/StablePartition.py`.
6. `punch_holes` contains the `iterate >= 0` assertion, and a unit test that passes a hand-built
   `iterate=-1` pair gets the `AssertionError`.
7. A unit test triggers the new WARNING. It builds a container whose neighbour crossing sits at
   `0.5 * SCALING_RTOL` (relative) from a **scaled** span end. The same setup with a table end
   (`from_table` True) does not trigger it. No law-case build logs the WARNING (the reviewer
   captures the `StablePartition` logger over the 11 law cases).
8. In `_hole_openings`, the anchor block calls `trellis.stable_key_at` and no longer rebuilds
   `(fp, "stable", orbit_index, branch_index)` inline. The anchor openings are equal to base on the
   four inversion cases (2/4, 2/2, 2/2, 4/8 on `s(o,1-b)`).
9. `Hole` has no `coords` field, and `Trellis.punch_holes` has no `epsilon` parameter (it is
   `inspect.signature(Trellis.punch_holes)` that shows this).
10. **Placement (a):** on every hole of `k28_two_blasts` and of all 11 law cases, `hole_point` on the
    branch equals the base's `Hole.coords` within 1e-12.
11. `diagnostics.measured_sides` returns the same sides on base and branch on the 11 law cases and
    on `p3_coarse_blasts`.

## 8. Questions

1. **Where do figures put a propagated hole?**
   - (a) Carry the reference point back with `trellis.dynamical_system.map_inv` inside
     `plotting.hole_point`. Plotting is a view, so the algorithm still calls no map. Figures stay
     identical (0.0 deviation), and the laws and `measured_sides` keep their oracle point.
   - (b) Use the no-map fallback everywhere. This is a stated figure change: propagated holes move by
     up to 0.33 and sit on the other side of their sub-arc on all 11 law cases. Rubric 10 becomes
     "the stated change", and the two geometric laws would need a new oracle point.

   **Default: (a).**
