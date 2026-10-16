# Request: B4a — `holes-algorithm-vs-diagnostics`

*From the author, 2026-10-13 ("onward!"). Input for one `tangle_refactor` run. The fifth slice of `docs/alignment_plan_2026-10-09.md`; B0, B1, B6a and B8 are merged (main `ae18859`; runs archived under `.refactor/runs/2026-10-12-{oracle-laws,resolved-crossings,walk-bounds-honest}/` and `.refactor/runs/2026-10-13-unresolved-sets/`).*

## The question

**Can `StablePartition.py` read as the hole and partition algorithm (plan A.6, A.7) alone, with every geometric self-diagnosis moved out?** Today the module (about 1,950 lines, roughly 990 of code) holds about 250 lines of algorithm and about 340 of readings that re-measure, on the polylines, what the crossing signs already decided: `_side_of`, `_anchorward_look`, `_approach_displacement`, `_oriented_bridge_polyline`, `_build_oriented_bridge_polyline`, `_bridge_side_of`, `RowEnd`, `RowDiagnosis`, `bridge_row_diagnosis`, `bridge_row_violation`, `check_bridge_rows_consistent`, `_orbit_side`, `SideDiagnosis`, `side_diagnoses`, `bridge_side_violations`, `check_holes_share_bridge_side`, `geometry_flags`, `_warn_geometry_flags`, plus `_map_backward`, the topology calling the real inverse map to fill `Hole.measured_side`, a flag the law `backward_holes_keep_bridge_side` re-measures independently. None of it feeds a decision. The readings go to a new `topology/diagnostics.py`, mirroring `numerics/diagnostics.py` (flags, never discards).

The author's own words on 2026-10-12, which prompted this slice's placement: the side of a hole is a cross product of the manifold's cdist direction with the direction to the hole, decided once at the reference hole and carried to every preiterate; "it doesn't have to be any fancier than that". The decisive code already is that (the crossing sign at `StablePartition.py:365`, `side_after` at `:535`, `_hole_openings` at `:1333`); the complexity is the diagnostic scaffolding around it. B4a moves it; B4b (next) takes out placement and the dead path.

## What "done" means

Exactly slice B4a of the alignment plan (`docs/alignment_plan_2026-10-09.md`, section B, "B4a — `holes-algorithm-vs-diagnostics`"): its scope, four steps, behaviour ("none"), not-doing list and five rubric items are the plan. In short:

1. Create `topology/diagnostics.py` with the 18 names above moved verbatim: one commit "move", then one commit "tidy", so the reviewer can diff the move as a rename.
2. `Hole` loses `measured_side` and `side_source`: the measured side is returned by `diagnostics.side_diagnoses(trellis)` on demand, computed there by `_map_backward`; `side_source` is `iterate == 0 → "signs", else "carried"`, a function not a field.
3. `Trellis.punch_holes(diagnose=True)` calls `diagnostics.warn_geometry_flags(self)` after the algorithm; the two invariant checks `check_bridge_rows_consistent` and `check_holes_share_bridge_side` stay assertions in the algorithm module (they are I1/I2, physics, not geometry), about 45 lines.
4. `topology/__init__.py` re-exports the moved names from `diagnostics` (no caller outside tests imports them by module path except `Trellis.py`).

Files (eight): `src/tanglepack/topology/StablePartition.py`, `src/tanglepack/topology/diagnostics.py` (new), `src/tanglepack/topology/Trellis.py` (`punch_holes(diagnose=)` imports), `src/tanglepack/topology/TopologyResults.py` (`Hole`), `src/tanglepack/topology/__init__.py`, and the test imports `tests/helpers/laws.py` (`_side_of`, `bridge_row_violation`), `tests/unit/topology/test_stable_partition.py` (`_bridge_side_of`, `_side_of`), `tests/plotting/test_plot_topology.py` (`_side_of`). A MOVE with import updates, no logic change.

**Behaviour.** None. `Hole` loses two fields; no library reader outside the diagnostics uses them (`plotting.py` reads `coords`, `bridge_side`, `openings` only; the rubric verifies by grep). Since B1 nothing in the hole algorithm changed, so the plan's line numbers may have moved only by B1's `Trellis.py` edits.

**Not doing.** Placement, the direct-pair path and the WARNING (B4b); any change to `row_of_end`, `_hole_openings`, `intervals_from_marks`.

**Rubric** (the plan's five): suite passes on the branch and the diff touches only the eight files; `StablePartition.py` ≤ 1,100 lines, and `topology/diagnostics.py` holds the 18 names and imports nothing from `StablePartition` except `row_of_end`, `rows_of_bridge`, `owns_cdist`; `grep -rn "measured_side\|side_source" src` hits only `topology/diagnostics.py`; `grep -rn "map_inv\|dynamical_system" src/tanglepack/topology/StablePartition.py` is empty (the algorithm never calls the map); fingerprints of every law case equal on base and branch, and `geometry_flags` / `side_diagnoses` output identical on the 11 law cases.

## Where to start

`docs/alignment_plan_2026-10-09.md` (B4a; A.6 holes, A.7 partition; B4b for what is deferred), `docs/code_fidelity_audit_2026-10-09.md` (the stage 5-6 findings), `src/tanglepack/topology/StablePartition.py` (module Dev Notes "Side", "Signs decide", "Invariants"; the function list), `src/tanglepack/topology/TopologyResults.py` (`Hole`), `src/tanglepack/topology/Trellis.py` (`punch_holes`), `src/tanglepack/numerics/diagnostics.py` (the model for the new module), `tests/unit/topology/test_topology_invariants.py` and `tests/unit/topology/test_production_checks.py` (the tests that read `measured_side`, `side_diagnoses`, `geometry_flags`; the analyst says whether they need an import change the plan's file list missed).

## Binding context

- The author's decisions on the alignment plan are in `.refactor/decisions.md` and plan section E; none gates B4a. Earlier author rules on holes stand: a hole's side is the side of its OWN bridge in the bridge's dynamical orientation, decided from the crossing signs at the reference hole and CARRIED to every preiterate (flipped per step only under an orientation-reversing map); holes propagate backward only; the signs decide the openings (D6); geometry only flags (NUM-1, NUM-2).
- `CLAUDE.md` invariants are physics: the topology never calls the map (B4a moves the one violation, `_map_backward`, into the diagnostics module, which is the flags layer; a TEST may call the map); I1 (holes of one origin share `bridge_side`) and I2 (equal signs on one stable branch are impossible) stay assertions in the algorithm.
- Test policy (`tests/conftest.py` docstring, 2026-10-05): no wording or private pins; rubric 5 is a base-vs-branch equality of fingerprints and of the diagnostics' output.
- The analyst stage is NOT skipped: the analyst verifies the 18 names and their line spans at `HEAD`, lists every reader of `measured_side` and `side_source` in `src/` and `tests/` (the plan names three test files; `test_topology_invariants.py` and `test_production_checks.py` also read them, so the file list may need widening and the analyst says so), confirms `plotting.py` reads only `coords`, `bridge_side`, `openings`, and reports what differs from the plan section. It does not re-plan the slice.
- The working tree may carry uncommitted documentation changes (`docs/`, `.refactor/`); the change branches from `HEAD` and never touches them.
