# Plan — B4a `holes-algorithm-vs-diagnostics`

*Base `main` @ `ae18859`. First plan for this slice; no earlier verdict to answer.*

## 1. The answer

Yes, with one adjustment. The rule that decides a hole's side is three lines: the crossing sign
(`StablePartition.py:365`), `side_after` (`:535`) and `_hole_openings` (`:1333`). Around it the
module carries about 600 lines that re-measure on the polylines what the signs already decided.
Those are 16 names, from `_side_of` to `_warn_geometry_flags`, and they move verbatim to a new
`topology/diagnostics.py`. `Hole` loses `measured_side` and `side_source`: the measured side
becomes a reading that diagnostics computes on demand, and `side_source` becomes a one-line
function.

The adjustment is that `_map_backward` stays in the algorithm for now. It does more than fill
`measured_side`: it also places every propagated hole, nudging it toward the carried point
(`:1585-1603`). Without the map, the fallback placement moves propagated holes on 10 of the 11 law
cases (for example, 29 of the 69 heteroclinic holes end up on the wrong side of their sub-arc). So
"no map in `StablePartition`" belongs to B4b, which owns placement. If it stays in B4a, the
behaviour changes.

## 2. Recommendation

Make the move, with the four corrections the analyst found folded in:

- the map stays where it places holes;
- the two invariant checks stay in the algorithm module and import their cores lazily;
- the line target becomes about 1,360 lines, not 1,100;
- three more test files go on the file list.

## 3. Scope

**Slug** `holes-algorithm-vs-diagnostics`. Branch `refactor/holes-algorithm-vs-diagnostics` in
`../tangle-pack-holes-algorithm-vs-diagnostics`, made from `ae18859`. The base worktree is
`../tangle-pack-holes-algorithm-vs-diagnostics-base`.

**Files (eleven):**

1. `src/tanglepack/topology/StablePartition.py`
2. `src/tanglepack/topology/diagnostics.py` (new)
3. `src/tanglepack/topology/Trellis.py`: the imports in `punch_holes`
4. `src/tanglepack/topology/TopologyResults.py`: `Hole` loses two fields
5. `src/tanglepack/topology/__init__.py`: the re-exports
6. `tests/helpers/laws.py`: imports, plus the `side_source` clause at `:1209`
7. `tests/unit/topology/test_stable_partition.py`: imports
8. `tests/plotting/test_plot_topology.py`: imports
9. `tests/unit/topology/test_topology_invariants.py`: imports, the `_hole` helper, the logger name
10. `tests/unit/topology/test_production_checks.py`: two measured-side tests, the logger name
11. `tests/regression/test_coarse_blasts_rows_by_signs.py`: the logger name

**Size.** About 600 lines move and about 60 change. The size comes from relocating code, not from
new logic, so I am not splitting the slice further.

## 4. Steps

Three commits, so the reviewer can read the move as a move.

**Commit 1, "move".** Cut these 16 names out of `StablePartition.py` byte for byte and paste them
into `diagnostics.py` in the same order:

> `_side_of`, `_anchorward_look`, `_approach_displacement`, `_oriented_bridge_polyline`,
> `_build_oriented_bridge_polyline`, `_bridge_side_of`, `RowEnd`, `RowDiagnosis`,
> `bridge_row_diagnosis`, `bridge_row_violation`, `_orbit_side`, `SideDiagnosis`,
> `side_diagnoses`, `bridge_side_violations`, `geometry_flags`, `_warn_geometry_flags`

`diagnostics.py` gets its own `logger = logging.getLogger(__name__)`. At module top it imports only
`row_of_end` and `rows_of_bridge` from `StablePartition`, plus `bridge_for_pair` in commit 2.

The two checks stay in `StablePartition.py` and import their cores inside the body, which is how
`Trellis.punch_holes` already does it. That avoids a circular import:

```python
def check_bridge_rows_consistent(trellis, bridge) -> None:
    from .diagnostics import bridge_row_violation
    message = bridge_row_violation(trellis, bridge)
    assert message is None, message

def check_holes_share_bridge_side(holes, *, orientation_preserving=True) -> None:
    from .diagnostics import bridge_side_violations
    messages = bridge_side_violations(holes, orientation_preserving=orientation_preserving)
    assert not messages, "; ".join(messages)
```

Two things do not change: the docstrings, and the I1/I2 message texts (`stress/numerics/report/analysis.py`
matches on them). `_bridge_side_of`, which the algorithm still calls at `:1606`, is reached through
a lazy import in this commit and is gone from the algorithm in commit 2.

`Trellis.punch_holes` keeps importing `check_*`, `punch_holes` and `propagate_reference_holes` from
`StablePartition` at call time; the spies in `test_production_checks.py` depend on that. It imports
the warning from `diagnostics`. `__init__.py` re-exports the moved names from `diagnostics` and the
checks from `StablePartition`, and `__all__` is unchanged. The test imports in files 6-11 follow,
and so do the logger names: `"tanglepack.topology.diagnostics"` for the moved WARNINGs and the skip
INFO. The suite passes at this commit.

**Commit 2, "measured side on demand".**

- In `_punch_in_bridge`, the `measured_side` lines go (`:1605-1607`). `carried` still places the
  hole, and `_map_backward` stays.
- The `side_source="signs"` / `"carried"` arguments go, `Hole` drops the two fields and their
  docstring paragraphs, and the Dev Notes at `:45` and `:61` are reworded to name the functions.
- `diagnostics.py` gains:

```python
def side_source(hole: Hole) -> SideSource:
    """A direct hole reads its side off its crossing signs; a propagated one carries it."""
    return "signs" if hole.pair is not None else "carried"

def measured_sides(trellis, holes) -> dict[tuple[tuple[int, int], int], Side]:
    """Each propagated hole's side MEASURED on the curves (a flag, D6): its orbit's
    reference hole carried back by the real inverse map, read on the image sub-arc.
    Empty without an inverse map. Keyed by (origin, iterate), unique per propagation."""
    # for each propagated hole h with reference ref (pair set, same origin):
    #   carried = ref.coords mapped back -h.iterate times by trellis.dynamical_system.map_inv
    #   rb = bridge_for_pair(trellis, ref.pair)
    #   span = sorted(trellis.image_cdist(i, h.iterate, "unstable", branch=rb.manifold_key)[1] for i in rb.id)
    #   side = _bridge_side_of(trellis, trellis.bridge_between(*h.bounding_ids), carried, span)
```

- `side_diagnoses(holes, *, orientation_preserving=True, measured=None)` keeps its signature and
  adds an optional `measured` mapping. `SideDiagnosis` keeps its `side_source` and `measured_side`
  fields, now filled from `side_source(hole)` and `measured.get((origin, iterate))`. It counts
  "unmeasured" only when `measured` is given, so the I1 path neither logs that INFO nor calls the
  map.
- `geometry_flags(trellis, bridges, holes)` is the one caller that builds `measured_sides(...)`.
  `_warn_geometry_flags` becomes the public `warn_geometry_flags`, with the same arguments, because
  `Trellis` calls it across modules.

The map therefore runs only when `diagnose=True`, and never on the assertion path.

**Commit 3, tests.**

- `laws.py:1209`: drop the `side_source` clause. It restated the field's assignment and would now
  test a function against its own definition.
- `test_topology_invariants._hole` builds a `Hole` without the two fields. Its three measured-side
  tests pass `measured={(origin, iterate): side}` to `side_diagnoses`.
- `test_production_checks`: "carry their orbit side" checks `side_source(hole)` and
  `measured_sides(trellis, trellis.holes)` instead of the fields. "Flags a contradicting measured
  side" monkeypatches `diagnostics.measured_sides` to flip one entry and still expects exactly one
  WARNING, now on the diagnostics logger.

## 5. Behaviour

None. Every `Hole` keeps its coordinates, ids, side, openings, iterate and origin. Every row reading
and every geometry flag is unchanged.

Two small things do change, and neither is a result:

- The flag WARNINGs and the skip INFO now come from the logger `tanglepack.topology.diagnostics`,
  with the same text.
- A hand-built `Hole` with no pair now describes itself as `carried` where it used to say `None`.

## 6. Not doing

- **Removing `_map_backward` and the carried placement.** That is B4b. B4b's step 4 says
  `hole_point` "reproduces today's placement", but what it describes is the no-map fallback, not
  today's placement. B4b's plan has to say which one it means.
- **Rest of B4b:** the direct-pair path and its fallback start (`:458-495`, which never runs on any
  build probed; `measured_sides` replays from the reference hole only), and B4b's one WARNING.
- **Unchanged here:** `row_of_end`, `_hole_openings` and `intervals_from_marks`.
- **Stale references outside the file list:**
  - The docstring at `numerics/geometry.py:7` names `StablePartition._build_oriented_bridge_polyline`.
  - The archived `scripts/refactor_proof/hole-side-and-own-blast/hole_sides.py` imports moved
    names from `StablePartition` and will break.
- **Placement already off its sub-arc today** (inversion_period2 4/8, heteroclinic 13/69). This
  affects plots only and goes to B4b.

## 7. Rubric

1. `MPLBACKEND=Agg env/bin/python -m pytest -q` passes on the branch. `-rxX` lists only `cases.KNOWN_ISSUES`.
2. Every pinning test, meaning files 6-11, passes on the base commit and on the branch.
3. `git diff --stat ae18859..` lists exactly the eleven files in section 3.
4. There is no behaviour change beyond section 5. A script run on base and branch over the 11
   `cases.LAW_CASES` plus `p3_coarse_blasts`, on each case's `session.topology().trellis` and on
   each per-fixed-point trellis, gives equal output on both sides for:
   - every hole as `(origin, iterate, bounding_ids, near_intersection_id, bridge_side, openings, coords)`;
   - the measured side of every propagated hole (base `hole.measured_side`, branch
     `measured_sides(trellis, trellis.holes)`);
   - `bridge_row_diagnosis` of every hole bridge, field by field;
   - `geometry_flags(...)` and `side_diagnoses(...)` as `describe()` strings.

   On `p3_coarse_blasts` that output is not empty (2 row flags, 20 measured flags), so the
   comparison can fail. Law-case fingerprints are equal on both sides.
5. In commit 1, each of the 16 names has source text in `diagnostics.py` (`ast.get_source_segment`)
   identical to its text in `StablePartition.py` at `ae18859`. The 16 names are no longer defined in
   `StablePartition.py`.
6. `wc -l src/tanglepack/topology/StablePartition.py` is at most 1,370 (base: 1,947).
7. `diagnostics.py` defines the 16 names plus `side_source`, `measured_sides` and
   `warn_geometry_flags`. Its only import from `StablePartition` is the top-level
   `row_of_end, rows_of_bridge, bridge_for_pair`. `StablePartition.py` imports from `diagnostics`
   only inside the bodies of the two `check_*` functions. `import tanglepack.topology.diagnostics`
   works in a fresh interpreter.
8. `grep -rn "measured_side\|side_source" src` hits only `topology/diagnostics.py`.
   `grep -n "measured_side\|side_source" src/tanglepack/topology/TopologyResults.py` is empty.
9. The map stays off the assertion path. With `diagnostics.measured_sides` monkeypatched to raise,
   `trellis.punch_holes(diagnose=False)` succeeds on `p3_coarse_blasts`, and
   `punch_holes(diagnose=True)` raises.
10. `grep -n "map_inv\|dynamical_system" src/tanglepack/topology/StablePartition.py` hits only
    `_map_backward`, which still exists. Under question 1's non-default answer this item becomes
    "empty" instead.
11. `set(tanglepack.topology.__all__)` is equal on base and branch.
12. `test_coarse_blasts_rows_by_signs.py` passes with `_LOGGER = "tanglepack.topology.diagnostics"`,
    and fails on the branch if that line is reverted.
13. The I1/I2 message phrases "its orbit expects", "on its image sub-arc" and "approaches its two
    defining crossings" still appear verbatim in `diagnostics.py`.

## 8. Questions

1. **Does the map's placement of propagated holes stay in `StablePartition` until B4b?**
   - **Default, yes.** `_map_backward` stays, "Behaviour: none" holds, and the "no map in
     `StablePartition`" rubric moves to B4b.
   - **If no,** B4a removes `_map_backward` and propagated holes take the fallback placement. That
     is the one behaviour change: propagated `coords` move on 10 of the 11 law cases (up to 29 of
     69 heteroclinic holes are drawn on the wrong side of their sub-arc; the largest move is 0.33 at
     k28 one blast). Rubric 4 then drops `coords`, and rubric 10 reads "empty".
