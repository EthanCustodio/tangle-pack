# Change: B4b `holes-are-combinatorial` (coder, iteration 1)

- **Base commit:** `f166573e08d64d48ae540176495ea1b0bcd41ba8` (the author's `main` HEAD; their uncommitted `docs/`, `.refactor/` and `artifacts/` changes were not touched).
- **Branch:** `refactor/holes-are-combinatorial`, worktree `/home/dezhu/Development/TanglePack/tangle-pack-holes-are-combinatorial`.
- **Base worktree:** `/home/dezhu/Development/TanglePack/tangle-pack-holes-are-combinatorial-base`, detached at `f166573`.
- **Run the code in a worktree with** `MPLBACKEND=Agg PYTHONPATH=src:tests ../tangle-pack/env/bin/python`. The `env/` editable install points at the main checkout's `src`. pytest picks up `src` from the worktree through `pyproject` `pythonpath`.

## Commits (one per plan step, in plan order)

```
391f952 Holes: punch the reference pairs only and assert no backward pair is recorded (B4b step 1)
a29aa0e Holes: the backward span keeps image_cdist's table flag and a scaled end beside a crossing warns (B4b step 2)
4029f95 Holes: an anchor opening is D6 on stable_key_at's branch plus one identification (B4b step 3)
26ed58c Plotting: hole_point places a hole on demand; plot_holes, measured_sides and the laws read it (B4b step 4)
145455c Holes: Hole.coords and punch_holes(epsilon=) go; StablePartition places nothing and calls no map (B4b step 5)
```

### Where the commits differ from the plan's step text

I wanted every intermediate commit to stay green on the tiers it touches. I ran `tests/unit/topology`, `tests/plotting` and `tests/facade` after each step, and the full suite at the end. That forced two things to move:

- **Step 1 keeps the `epsilon=` keyword**, and **step 2 keeps the carried point, `_map_backward` and `_punch_in_bridge`'s `carried`/`span` parameters.** `Hole.coords` is a required field until step 5, so propagated holes still need a placement until then. All of these go in **step 5**, together with the field.
- **Step 4 copies the four placement helpers into `plotting.py`**, and step 5 deletes the originals in `StablePartition.py`. Between those two commits the helpers exist twice. At the step-4 commit I ran the probe and compared that commit's `hole_point` against the base's `Hole.coords`: 0.0 on every hole. The same result holds at the final commit (below). The end state is exactly the plan's.

I also made one small wording fix to a docstring line written in step 1. It was amended into the step-5 commit. It is docstring-only and came after the full-suite run.

## `git diff --stat f166573`

```
 src/tanglepack/loom/TangleSession.py            |   2 +-
 src/tanglepack/topology/StablePartition.py      | 471 ++++++------------------
 src/tanglepack/topology/TopologyResults.py      |   9 +-
 src/tanglepack/topology/Trellis.py              |   4 +-
 src/tanglepack/topology/diagnostics.py          |   7 +-
 src/tanglepack/topology/plotting.py             | 213 ++++++++++-
 tests/helpers/laws.py                           |  13 +-
 tests/unit/topology/test_minimal_trellis.py     |   2 +-
 tests/unit/topology/test_stable_partition.py    | 104 +++++-
 tests/unit/topology/test_topology_invariants.py |   1 -
 10 files changed, 423 insertions(+), 403 deletions(-)
```

Numstat (+/−):

| file | + | − |
|---|---|---|
| `StablePartition.py` | 104 | 367 |
| `plotting.py` | 207 | 6 |
| `test_stable_partition.py` | 94 | 10 |
| `laws.py` | 7 | 6 |
| `diagnostics.py` | 4 | 3 |
| `TopologyResults.py` | 4 | 5 |
| `Trellis.py` | 1 | 3 |
| `TangleSession.py` | 1 | 1 (the docstring line) |
| `test_minimal_trellis.py` | 1 | 1 |
| `test_topology_invariants.py` | 0 | 1 |

The plan estimated −300/+120. The real totals are larger because:
- about 150 lines are the four helpers moving verbatim into `plotting.py`;
- about 85 lines are the two new unit tests.

Source alone is net about −60 lines.

## Suite

- **Branch:** `MPLBACKEND=Agg env/bin/python -m pytest -q -rxX` gives **1176 passed, 38 skipped, 5 deselected, 1 xfailed in 204.31s**. The only xfail is `test_orientation_reversing_case_builds[orientation_reversing]` (C3: INV-11), which is a `cases.KNOWN_ISSUES` entry.
- **Base, same command:** **1173 passed, 38 skipped, 5 deselected, 1 xfailed in 204.67s**, with the same xfail.
- The +3 on the branch are the new tests:
  - `test_a_backward_pair_is_never_punched`
  - `test_scaled_span_end_beside_a_crossing_warns[scaled]`
  - `test_scaled_span_end_beside_a_crossing_warns[table]`
- I also ran the touched tiers after each step (`tests/unit/topology`, `tests/plotting`, `tests/facade`). Every run was green.

## What changed (final state)

- **`StablePartition.punch_holes(trellis, pairs=None)`**
  - Asserts `all(p.iterate >= 0 for p in pairs)`, with the plan's message.
  - Skips forward pairs as today, so it punches the references only.
  - Builds `Hole` from `near_intersection_id`, `pair`, `bridge_side`, `bounding_ids`, `openings`, `iterate` and `origin`.
  - Deleted: the else-branch placement and the `arc_point is None` skip (placement-only).
- **`propagate_reference_holes`**
  - A reference without a hole is skipped.
  - Deleted: `seen`, `punched_regions`, `_region_key`, the `fallbacks` path, its "cannot propagate" WARNING and its "has no punched bridge" INFO.
  - The span keeps `image_cdist`'s table flag (`ends = [...]`, `span = tuple(sorted(c for _, c, _ in ends))`) under the plan's corrected comment.
  - **New WARNING** (one per branch per call, logger `tanglepack.topology.StablePartition`). It fires when a **scaled** end (`from_table` False) lies within `SCALING_RTOL * |x|` of a crossing `x` of the container's branch. It is checked whenever a container is found.
  - `_punch_in_bridge(trellis, bridge, *, bridge_side, iterate, origin) -> Hole` builds the hole from ids, side and openings alone.
- **`_hole_openings` anchor block:** `key = trellis.stable_key_at(iid, trellis.unstable_key_of(bridge.id))`. When D6 says anchorward, the comment states the identification and the key becomes `s(o, 1-b)`. If that branch was never computed, nothing opens.
- **`_hole_marks_on_branch`:** the debug line names `hole.bounding_ids`.
- **Removed from `StablePartition.py`:**
  - `numpy` and `NDArray`;
  - `_step_into`, `_stable_arc_midpoint`, `_nearest_arc_point`, `_image_arc_midpoint` and `_map_backward`.
- **Kept in `StablePartition.py`:** `_near_far`, because it orders `bounding_ids`. The Dev Notes are updated: holes are combinatorial, the new WARNING, the propagation start, and the carried point as a view only.
- **`plotting.hole_point(trellis, hole, epsilon=0.05)`:**
  - **Reference hole:** steps off the stable-arc middle of `bounding_ids` toward the nearest bridge node.
  - **Propagated hole:** finds the origin's reference hole in `trellis.holes`, carries its point back `-iterate` steps with `trellis.dynamical_system.map_inv`, reads the span from the reference bridge's ids by `image_cdist`, and nudges the sub-arc midpoint toward the carried point.
  - **No map:** falls back to the sub-arc midpoint nudged toward the chord, as today.
  - The four helpers are verbatim copies.
  - `plot_holes` calls `hole_point` and reads only `origin`, `iterate` and the point.
- **`Hole.coords` removed** (its docstring points at `hole_point`). **`Trellis.punch_holes(epsilon=)` removed**, along with the `epsilon` mention in the `TangleSession.punch_holes` docstring.
- **`diagnostics.measured_sides`:** imports `hole_point` inside the function and reads the reference point from it. The stale "until B4b" note is gone.
- **`laws.py`:** `check_holes_are_classified`, `check_backward_holes_keep_bridge_side` and `check_direct_hole_side_is_coordinate_side` call `hole_point`.
- **Tests:**
  - The three `Hole(coords=)` go.
  - `test_henon_propagated_hole_side_matches_coords` is renamed `..._point_is_on_its_side` and reads `hole_point`.
  - The `_is_forward_iterate` parametrize comment for `iterate=-1` is reworded (the function is unchanged).
  - Two new tests are added (rubrics 6 and 7).

## Behaviour changes

- **Topology output: none.** On all 11 law cases the following are equal on base and branch: every hole's five fields, the full law fingerprint (crossings, bridges, holes, pips, partition signature, alphabet, generations), and `measured_sides`. `p3_coarse_blasts` also matches on holes, points and measured sides. Probe output is below.
- **Figures: identical.** `hole_point` differs from the base's `Hole.coords` by **0.0** on every hole of every probed case.
- **`epsilon` (API detail):** `hole_point`'s `epsilon` now scales the propagated-hole nudge too. Before, that nudge used a literal `0.05` whatever `punch_holes(epsilon=)` was. At the default the two are identical, and no caller passed a non-default value.
- **API:**
  - `Hole(coords=...)` and `punch_holes(..., epsilon=...)` are gone. No script, stress file or notebook used them (grep).
  - `StablePartition.punch_holes` raises `AssertionError` on a pair with `iterate < 0`.
- **Logs** (as in plan section 5):
  - `near_tangency` (non-law builder) logs the 14 "Orbit … has no punched pair … cannot propagate" WARNINGs on base and **none** on the branch. Its "Skipped 14 pair(s)" WARNING is still there.
  - The new NUM-5 WARNING fires once on `k28_deep_blast_sliver`, and on no law case.
  - The "Reference … has no punched bridge" INFO path is gone. It never fired.

## Probe outputs (attached)

The scripts are in the session scratchpad `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/`:
- `b4b_rubric.py` runs in a worktree and dumps the fingerprint, five fields, points, `measured_sides` and `StablePartition` WARNINGs for the 11 law cases plus `p3_coarse_blasts`. On base the point is `Hole.coords`; on the branch it is `plotting.hole_point`.
- `b4b_compare.py` replaces FixedPoint reprs with `FP0`, `FP1`, … per case and compares.
- `b4b_extra.py` covers rubrics 8 and 9 and the non-law WARNINGs.
- Dumps: `b4b_base.json`, `b4b_branch.json`.

`b4b_compare.py b4b_base.json b4b_branch.json`:

```
k10                  holes=  3 fingerprint_equal=True five_fields_equal=True same_holes=True max|dpoint|=0 measured_equal=True (n=2) SP_warnings base=0 branch=0 new_kind=0
k28_one_blast        holes=  4 fingerprint_equal=True five_fields_equal=True same_holes=True max|dpoint|=0 measured_equal=True (n=2) SP_warnings base=0 branch=0 new_kind=0
k28_two_blasts       holes=  4 fingerprint_equal=True five_fields_equal=True same_holes=True max|dpoint|=0 measured_equal=True (n=3) SP_warnings base=0 branch=0 new_kind=0
p3                   holes=  8 fingerprint_equal=True five_fields_equal=True same_holes=True max|dpoint|=0 measured_equal=True (n=6) SP_warnings base=0 branch=0 new_kind=0
nested               holes= 12 fingerprint_equal=True five_fields_equal=True same_holes=True max|dpoint|=0 measured_equal=True (n=9) SP_warnings base=0 branch=0 new_kind=0
nested_inner_blasts  holes= 12 fingerprint_equal=True five_fields_equal=True same_holes=True max|dpoint|=0 measured_equal=True (n=9) SP_warnings base=0 branch=0 new_kind=0
inversion            holes=  6 fingerprint_equal=True five_fields_equal=True same_holes=True max|dpoint|=0 measured_equal=True (n=4) SP_warnings base=0 branch=0 new_kind=0
inversion_one_blast  holes=  7 fingerprint_equal=True five_fields_equal=True same_holes=True max|dpoint|=0 measured_equal=True (n=6) SP_warnings base=0 branch=0 new_kind=0
inversion_two_blasts holes=  8 fingerprint_equal=True five_fields_equal=True same_holes=True max|dpoint|=0 measured_equal=True (n=6) SP_warnings base=0 branch=0 new_kind=0
inversion_period2    holes= 10 fingerprint_equal=True five_fields_equal=True same_holes=True max|dpoint|=0 measured_equal=True (n=8) SP_warnings base=0 branch=0 new_kind=0
heteroclinic         holes= 73 fingerprint_equal=True five_fields_equal=True same_holes=True max|dpoint|=0 measured_equal=True (n=69) SP_warnings base=0 branch=0 new_kind=0
p3_coarse_blasts     holes= 84 fingerprint_equal=True five_fields_equal=True same_holes=True max|dpoint|=0 measured_equal=True (n=80) SP_warnings base=0 branch=0 new_kind=0
```

The `SP_warnings` count is every WARNING on the `StablePartition` logger during the build (`build_products` plus `topology()`). `new_kind` counts the new NUM-5 WARNING.

`b4b_extra.py near_tangency k28_deep_blast_sliver`, **base**:

```
punch_holes signature: (self, pairs: 'Optional[Iterable[PseudoneighborPair]]' = None, *, epsilon: 'float' = 0.05, propagate: 'bool' = True, diagnose: 'bool' = True, verbose: 'bool' = False) -> 'list[Hole]'
Hole fields: ['coords', 'near_intersection_id', 'pair', 'bridge_side', 'bounding_ids', 'openings', 'iterate', 'origin']
inversion: anchor ends 4, opened on s(o,1-b) 2
inversion_one_blast: anchor ends 2, opened on s(o,1-b) 2
inversion_two_blasts: anchor ends 2, opened on s(o,1-b) 2
inversion_period2: anchor ends 8, opened on s(o,1-b) 4
near_tangency: 15 StablePartition WARNING(s): {'Skipped 14': 1, 'Orbit (3,': 1, 'Orbit (5,': 1, 'Orbit (7,': 1, 'Orbit (9,': 1, 'Orbit (11,': 1, 'Orbit (13,': 1, 'Orbit (15,': 1, 'Orbit (17,': 1, 'Orbit (19,': 1, 'Orbit (21,': 1, 'Orbit (25,': 1, 'Orbit (27,': 1, 'Orbit (31,': 1, 'Orbit (33,': 1}
k28_deep_blast_sliver: 0 StablePartition WARNING(s): {}
```

**Branch**:

```
punch_holes signature: (self, pairs: 'Optional[Iterable[PseudoneighborPair]]' = None, *, propagate: 'bool' = True, diagnose: 'bool' = True, verbose: 'bool' = False) -> 'list[Hole]'
Hole fields: ['near_intersection_id', 'pair', 'bridge_side', 'bounding_ids', 'openings', 'iterate', 'origin']
inversion: anchor ends 4, opened on s(o,1-b) 2
inversion_one_blast: anchor ends 2, opened on s(o,1-b) 2
inversion_two_blasts: anchor ends 2, opened on s(o,1-b) 2
inversion_period2: anchor ends 8, opened on s(o,1-b) 4
near_tangency: 1 StablePartition WARNING(s): {'Skipped 14': 1}
k28_deep_blast_sliver: 1 StablePartition WARNING(s): {'Backward image': 1}
    Backward image -1 of reference (185, 186) on A:u(0.0): its scaled end 41203463.1 lies within SCALING_RTOL of crossing 106, so its container (110, 111) could be the neighbour bridge (NUM-5); grow or blast to register the image
```

## Rubric, my reading

1. **Pass.** 1176 passed, 38 skipped, 1 xfailed. The xfail is a `KNOWN_ISSUES` entry (C3 INV-11), the same one as on base.
2. **Pass.**
   - The base suite passes in full (1173), and that includes every hole-touching tier: law tier, `test_stable_partition`, `test_topology_invariants`, `test_minimal_trellis`, `test_production_checks`, `tests/plotting`.
   - The branch suite passes in full.
   - The two new tests (three items) and the tests edited to drop `Hole(coords=)` pass on the branch. One base test was renamed and now reads `hole_point`: `test_henon_propagated_hole_point_is_on_its_side`.
3. **Pass.** The diff touches exactly the nine files plus `TangleSession.py`, and the `TangleSession.py` diff is the single docstring line (`-`/`+` shown by `git diff`).
4. **Pass.**
   - The law fingerprint, without the process-local `id()` of the products, is equal on base and branch for all 11 law cases.
   - The holes' five fields are equal as well.
   - The only behaviour changes are the section-5 logs.
5. **Pass.** `StablePartition.py` is **1091** lines (≤ 1,100). The rubric grep prints nothing (exit 1).
6. **Pass.** `punch_holes` has `assert all(p.iterate >= 0 for p in pairs), "holes propagate backward only; no backward pair is recorded"`. `test_a_backward_pair_is_never_punched` passes a hand-built `iterate=-1` pair and gets the `AssertionError`.
7. **Pass.**
   - `test_scaled_span_end_beside_a_crossing_warns`, scaled case: the container `(ids[0], ids[1])` has its far crossing at `2.0` and the span end sits at `2.0 * (1 - 0.5 * SCALING_RTOL)`. Exactly 1 WARNING.
   - Table case (`from_table` True): 0 WARNINGs.
   - None of the 11 law-case builds logs it (`new_kind=0` above, captured on the `StablePartition` logger).
8. **Pass.**
   - The anchor block calls `trellis.stable_key_at` and no longer rebuilds `(fp, "stable", orbit_index, branch_index)` from `unstable_key_of`.
   - The identification unpacks the key it got and flips `branch_index`, with the comment stating why.
   - Anchor openings on `s(o,1-b)` are 2/4, 2/2, 2/2 and 4/8, equal on base and branch.
   - All openings are equal, because the five fields are equal.
9. **Pass.**
   - The `Hole` fields are `[near_intersection_id, pair, bridge_side, bounding_ids, openings, iterate, origin]`.
   - `inspect.signature(Trellis.punch_holes)` is `(self, pairs=None, *, propagate=True, diagnose=True, verbose=False)`.
10. **Pass.** `max|dpoint| = 0` on every hole of `k28_two_blasts`, of all 11 law cases, and of `p3_coarse_blasts` (branch `hole_point` against base `Hole.coords`).
11. **Pass.** `measured_sides` is equal on base and branch for the 11 law cases and `p3_coarse_blasts`: 2, 2, 3, 6, 9, 9, 4, 6, 6, 8, 69 and 80 entries.

## Notes for the test engineer / reviewer

- `hole_point` of a propagated hole needs that orbit's reference hole in `trellis.holes`, which `punch_holes` always provides. On a hole whose reference is missing it raises `StopIteration`. No caller does that.
- `CLAUDE.md` still names `Hole.coords`. That is left for B12, as the plan says.
