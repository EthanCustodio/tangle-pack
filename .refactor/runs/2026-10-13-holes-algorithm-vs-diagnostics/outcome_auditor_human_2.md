# Outcome: B4a `holes-algorithm-vs-diagnostics`

**Main has uncommitted changes**, but the branch does not touch them. They are in `.refactor/`, `docs/` and `artifacts/`, plus untracked run folders.

- **Branch.** `refactor/holes-algorithm-vs-diagnostics` @ `f166573`, worktree `../tangle-pack-holes-algorithm-vs-diagnostics`, clean.
- **Base.** `ae18859`.
- **Review.** Passed in round 1. There is no debrief verdict yet.
- **Merge.** Not merged; that is your call.

## The answer now

Yes. The algorithm is now a module of its own, and the diagnostics are a second module.

**`StablePartition.py` is the hole and partition algorithm**, and nothing else:

- It punches direct holes and reads their side off the crossing signs.
- It carries that side backward along the orbit with `side_after`. This is your "cross product decided once at the reference, carried to every preiterate".
- It opens the holes by sign arithmetic (`_hole_openings`, D6) and builds the intervals.
- It asserts I1 and I2.
- It no longer measures anything on a polyline.

**`topology/diagnostics.py`** (new, the counterpart of `numerics/diagnostics.py`) holds every reading that re-measures on the polylines what the signs already decided:

- the geometric row of a bridge end;
- the hole side measured by carrying the reference hole back with the real inverse map.

They come out only as flags: one WARNING per kind, from `punch_holes(diagnose=True)`.

**One exception, which you approved:** `_map_backward` stays in the algorithm until B4b. The same carried point also places propagated holes on the plot.

## What changed

Four commits:

1. `95ea780`: a verbatim move of 16 names.
2. `728dcbb`: the measured side is computed on demand, and `Hole` loses two fields.
3. `fc9f9a4`: the coder's test updates.
4. `f166573`: two tests added by the test engineer.

The most important hunk is in `_punch_in_bridge`. The algorithm stops measuring, and `Hole` stops storing the measurement.

Before (`ae18859`, `StablePartition.py:1605`):
```python
    measured_side = (
        _bridge_side_of(trellis, bridge, carried, span) if carried is not None else None
    )
    return Hole(..., iterate=iterate, origin=origin,
                measured_side=measured_side, side_source="carried")
```

After (`f166573`, `StablePartition.py:1017`; the measurement moves to `diagnostics.measured_sides` at `:429`):
```python
    return Hole(..., openings=_hole_openings(trellis, bridge, bridge_side),
                iterate=iterate, origin=origin)
```

`measured_sides(trellis, holes)` replays the measurement. For each propagated hole it carries the orbit's reference hole back `-iterate` steps with `map_inv` and reads that point with `_bridge_side_of` on the hole's bridge. It returns `{(origin, iterate): side}`.

`side_source(hole)` is now a one-line function: `"signs"` if the hole has a pair, otherwise `"carried"`.

I1 and I2 stay as three-line assertions in `StablePartition` (`:699`, `:719`). Each lazily imports its core from `diagnostics`, because a module-level import would be circular.

## What got shorter

| file | lines (base → branch) | code lines, no docstrings or comments |
|---|---|---|
| `StablePartition.py` | 1,947 → **1,354** | 926 → 608 |
| `diagnostics.py` | new, 712 | 365 |
| `TopologyResults.py` | 938 → 923 | 309 → 306 |

The plan's target was at most 1,370 lines; your original target of 1,100 needs B4b as well. `StablePartition` went from 42 top-level definitions to 26. `git diff --stat ae18859`: 11 files, +853 / −704.

## Behaviour changes

No result changes. Every hole (including its `coords`), every measured side, every row diagnosis, every flag and every law-case fingerprint is identical. What does change:

1. **Logger name.** The geometry-flag WARNINGs and the "Hole-side check skips …" INFO now come from `tanglepack.topology.diagnostics`, with the same text. Anything filtering on `tanglepack.topology.StablePartition` will miss them.
2. **Hand-built holes.** A hand-built `Hole` with no pair now describes itself as `carried` in `SideDiagnosis.describe()`; before, it said `None`.
3. **`SideSource` moved** from `TopologyResults` to `diagnostics`. It was not in `__all__`, and nothing else imported it.
4. **API names.**
   - `Hole(measured_side=, side_source=)` no longer exists.
   - `_warn_geometry_flags` is now `warn_geometry_flags`.
   - `side_diagnoses` gains an optional `measured=` argument.
   - `tanglepack.topology.__all__` is unchanged (101 names).
5. **The I1 path calls no map.** `check_holes_share_bridge_side` no longer counts unmeasured holes and never calls the map. The map runs only on the flag path, so it appears at most once per `punch_holes(diagnose=True)`.
   - Side effect: with `diagnose=True` (the default), `map_inv` now runs **twice** per propagation step, once to place the hole and once to replay the measurement. It used to run once. This is a negligible cost for Hénon, and it ends when B4b removes the placement.

## What I verified, and how

- **Full suite**, run by me on `f166573`: 1173 passed, 38 skipped, 5 deselected, 1 xfailed. The only xfail is the `orientation_reversing` KNOWN_ISSUE.
- **Pinning files on base** (topology invariants, production checks, stable partition, plot topology, coarse-blast rows, `tests/invariants`): 227 passed, 38 skipped, 1 xfailed.
- **Commit 1 is a byte-identical move.** My own `ast` script, with decorators included, compared each name's source at `ae18859` with `95ea780`:
  - all 16 names are identical, in the same order, and none is left behind;
  - the only other `StablePartition` edits in that commit are the three planned lazy imports and narrowed module imports (`bisect_left`, `dataclass` and `oriented_bridge_polyline` are no longer needed there).
- **File list.** `git diff --stat` lists exactly the eleven planned files.
- **Base vs branch output.** I re-ran the coder's comparison script myself in both worktrees. It covers the 11 law cases plus `p3_coarse_blasts`, on the topology trellis and on every per-fixed-point trellis. It dumps holes with `coords`, measured sides, row diagnoses, `geometry_flags` / `side_diagnoses` and fingerprints.
  - After stripping object addresses, the two dumps are identical: 1,258 lines each, md5 `74385b04…`.
  - The comparison can fail: `p3_coarse_blasts` flags 2 rows and 20 sides on both sides.
- **Greps.**
  - `measured_side|side_source` in `src` hits only `diagnostics.py`.
  - `map_inv|dynamical_system` in `StablePartition.py` hits only `_map_backward` (`:955`, `:958`).
- **Read by me, not re-run:** the test engineer's break-it table (13 mutations, each caught) and the reviewer's checks of R9 (map off the assertion path) and R12 (the logger line has teeth).

## What did not go to plan

- **Commit 2 alone is red.** 12 tests fail there until commit 3. The plan split source and tests this way and only required commit 1 to be green; that commit is green (1171 passed, per the coder).
- **Four small things went beyond the plan's letter, all inside the file list:**
  - `SideSource` moved to `diagnostics`.
  - Two `Trellis.py` docstrings were repointed to `diagnostics`.
  - The test engineer added two tests to close real gaps: `side_source` had no test that could fail, and the "map off the assertion path" rule had no test.
  - The `side_source` clause in `laws.py` was dropped as a tautology, as planned.
- **Your own "18 names" became 16 moved, plus `side_source` and `measured_sides`.** The two checks stayed in the algorithm, and `_map_backward` stayed until B4b, as approved.

## Where I am least sure this is what you meant

1. **The map is still in the algorithm.** `_map_backward` (`StablePartition.py:951`, called at `:513`) still places every propagated hole. So "the topology never calls the map" is still broken in `StablePartition` until B4b. With `diagnose=True` the map is now called twice per step. B4b's plan has to say which placement it means: today's carried point, or the no-map midpoint fallback. The fallback misplaces holes on 10 of the 11 law cases.
2. **I1 and I2 run code that lives in `diagnostics`.** They stay assertions in `StablePartition`, but each is a three-line wrapper. The logic that decides the assertion (`bridge_side_violations` → `side_diagnoses` → `_orbit_side`, and `bridge_row_violation`) now lives in the module that is supposed to hold "flags, never decisions". You asked for about 45 lines of checks in the algorithm, which reads as the cores too. If so, a follow-up would move `bridge_row_violation`, `bridge_side_violations` and `_orbit_side` back. That needs care with the circular import.
3. **Loose ends outside the file list:**
   - **Stale references, left alone:**
     - `numerics/geometry.py:7` still names `StablePartition._build_oriented_bridge_polyline`.
     - `scripts/refactor_proof/hole-side-and-own-blast/hole_sides.py` calls `sp._bridge_side_of`, `sp._oriented_bridge_polyline` and `sp.bisect_left`, so it will break.
     - `CLAUDE.md` still describes `Hole.measured_side` and `side_source` in the holes bullet (`:35`) and in the `Hole` table row (`:139`).
   - **A case `measured_sides` cannot replay:** an orbit whose backward chain started from the propagation fallback. That is a reference with no spanning bridge, where propagation starts from another pair. `measured_sides` starts only from iterate 0 at the reference, so those holes would show as "unmeasured" (INFO) where base measured them. The fallback never runs on any build we have probed, and B4b removes it.

## Files

- Report: `/home/dezhu/Development/TanglePack/tangle-pack/artifacts/outcome_auditor_human_2.md`
- Decision appended: `/home/dezhu/Development/TanglePack/tangle-pack/.refactor/decisions.md`
- Knowledge updated: `/home/dezhu/Development/TanglePack/tangle-pack/.refactor/knowledge/topology.md` (Map rewritten for the branch, plus an audit section) and `.../testing.md`
