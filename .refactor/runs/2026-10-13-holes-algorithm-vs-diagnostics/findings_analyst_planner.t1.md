# Findings — B4a `holes-algorithm-vs-diagnostics` (analyst, HEAD `ae18859`, 2026-10-13)

*(The analyst's Write tool refused the artifact; the orchestrator saved this text verbatim.)*

**Dirty tree.** `main` @ `ae18859` has uncommitted changes in `.refactor/decisions.md`, `.refactor/knowledge/{loom,numerics,testing,topology}.md`, `docs/alignment_plan_2026-10-09.md` and `docs/code_fidelity_audit_2026-10-09.md`. There is also a staged rename (`artifacts/request_origin_analyst.md -> .refactor/runs/2026-10-12-oracle-laws/`) and untracked `.refactor/runs/*` and `artifacts/`. No tracked file under `src/` or `tests/` is dirty, so a branch from `HEAD` starts from exactly the code read here. The B4a section I checked is in the uncommitted plan doc.

**What kind of request this is.** A refactor (a move). Targets: `topology/StablePartition.py`, a new `topology/diagnostics.py`, `Hole` in `TopologyResults.py`, `Trellis.punch_holes` and `topology/__init__.py`. As asked, this checks the plan against the code and does not re-plan the slice.

**Short answer.** The direction is sound and all 18 names exist under the plan's names. Four facts at `HEAD` stop the slice being done exactly as written:
- **(A)** `_map_backward` does more than fill `measured_side`: it also places every propagated hole.
- **(B)** The two checks the plan keeps in `StablePartition.py` call into the names that move. That makes a circular import, and it puts the map on the I1 path.
- **(C)** The ≤ 1,100-line target cannot be reached by B4a; the file lands at about 1,348 lines.
- **(D)** The file list misses three test files, and `laws.py` needs more than an import change.

On top of that, the rubric-5 baseline is empty on all 11 law cases, so that check cannot fail as written.

## 1. The 18 names and `_map_backward` at `HEAD` (verified: `grep -n` and Read)

All 19 exist as top-level definitions. Spans run from the `def`/`@dataclass` line to the end of the body.

| # | name | span | note |
|---|---|---|---|
| 1 | `_side_of` | 626-638 | |
| 2 | `_anchorward_look` | 641-665 | |
| 3 | `_approach_displacement` | 668-682 | |
| 4 | `_oriented_bridge_polyline` | 685-714 | memoised via `trellis.cached_bridge_polyline` |
| 5 | `_build_oriented_bridge_polyline` | 717-734 | named in a docstring at `numerics/geometry.py:7`, outside the file list |
| 6 | `_bridge_side_of` | 737-781 | the algorithm also calls it at :1606 (`measured_side`) |
| — | `row_of_end` / `rows_of_bridge` | 784-835 / 838-857 | these STAY, and sit between the two moved blocks |
| 7 | `RowEnd` | 860-900 | |
| 8 | `RowDiagnosis` | 903-961 | |
| 9 | `bridge_row_diagnosis` | 964-1023 | uses `row_of_end` |
| 10 | `bridge_row_violation` | 1026-1070 | the I2 core: `rows_of_bridge` + `stable_key_at`; its message ends with `bridge_row_diagnosis(...).describe()` |
| 11 | `check_bridge_rows_consistent` | 1073-1088 | asserts `bridge_row_violation(...) is None` |
| 12 | `_orbit_side` | 1091-1096 | |
| 13 | `SideDiagnosis` | 1099-1149 | has fields `side_source` and `measured_side` |
| 14 | `side_diagnoses` | 1152-1221 | signature `(holes, *, orientation_preserving=True)`; reads the two `Hole` fields at :1195 and :1213-1215 |
| 15 | `bridge_side_violations` | 1224-1245 | the I1 core: `side_diagnoses` filtered by `orbit_disagrees` |
| 16 | `check_holes_share_bridge_side` | 1248-1270 | asserts `not bridge_side_violations(...)` |
| 17 | `geometry_flags` | 1273-1308 | `(trellis, bridges, holes)` |
| 18 | `_warn_geometry_flags` | 1311-1330 | logs on the StablePartition logger |
| + | `_map_backward` | 1538-1546 | called once, at :516 |

- **Size.** The file is 1,947 lines; the blocks are 626-783 and 860-1332, plus 1538-1548.
- **What the moved block uses.** `row_of_end`, `rows_of_bridge`, `side_after`, `oriented_bridge_polyline`, `bisect_*`, `Hole`, `Side`, `SideSource` and `logger`. No moved name uses `owns_cdist`, so rubric 2's allowance for it is unused but harmless.
- **Map calls (verified).** `grep "map_inv|dynamical_system|\.map(|jacobian"` hits only :1542 and :1545, both inside `_map_backward`. It is the only map call in the module.

## 2. Who reads `Hole.measured_side` and `Hole.side_source` (verified: grep over src, tests, scripts, stress)

**In `src/`:** only `StablePartition.py`, at:
- Dev Notes :45 and :61;
- `side_source="signs"` at :377;
- the diagnosis block, :1108-1215;
- `_punch_in_bridge`, :1565-1617.

The `Hole` definition is at `TopologyResults.py:131-174`. `SideSource` (`TopologyResults.py:61`) is used by those two files only.

**In `tests/`, six files must change, not the plan's three:**

| file | what it uses | change needed |
|---|---|---|
| `tests/helpers/laws.py` | imports `_side_of`, `bridge_row_violation` (:70); `check_holes_are_classified` **asserts `hole.side_source` (:1209)** | the import, plus a logic edit at :1209. The law `check_backward_holes_keep_bridge_side` (:1219-1265) re-measures with its own `map_inv` and `_point_side_of_bridge` and reads neither field. |
| `tests/unit/topology/test_stable_partition.py` | `_bridge_side_of`, `_side_of` (:24-26; used at :104-142, :373, :696) | import only |
| `tests/plotting/test_plot_topology.py` | `_side_of` (:25, :110) | import only |
| **`tests/unit/topology/test_topology_invariants.py`** | imports 7 moved names from `StablePartition` (:44-52); the `_hole` helper builds `Hole(measured_side=, side_source=)` (:61-73), used by 3 tests (:155, :174, :197); INFO logger `"tanglepack.topology.StablePartition"` (:200-202) | imports, the helper, the 3 tests, the logger name |
| **`tests/unit/topology/test_production_checks.py`** | reads both fields (:118-133); flips `hole.measured_side` to trigger the flag (:136-154); `_LOGGER` set to StablePartition (:29, :152-154) | rewrite 2 tests and the logger name |
| **`tests/regression/test_coarse_blasts_rows_by_signs.py`** (the request did not name it) | `_LOGGER = "tanglepack.topology.StablePartition"` (:25); asserts `punch_holes()` logs more than twice as many WARNINGs (`> 2 * quiet`) | logger name. It is in the default suite and fails once the geometry WARNINGs move to `topology.diagnostics`. |

- **The spies.** `test_production_checks.py:61-92` spies on `stable_partition.check_*`. They survive only if the checks stay in `StablePartition` and `Trellis.punch_holes` keeps importing them from there at call time, as it does today (`Trellis.py:1240`).
- **Covered by the re-export.** `test_open_p3_deep_runs.py:33` and `test_coarse_blasts_rows_by_signs.py:23` import from `tanglepack.topology`, so step 4 covers them.
- **No change needed.** `test_exact_ownership.py:150` captures an INFO that stays in StablePartition.

**Outside `src/` and `tests/`:**
- `scripts/refactor_proof/hole-side-and-own-blast/hole_sides.py` calls `sp._bridge_side_of`, `sp._oriented_bridge_polyline` and `sp.check_holes_share_bridge_side`. It is archived and will break.
- `stress/numerics/report/analysis.py:681-689` sorts failures by the I1/I2 message texts ("its orbit expects", "on its image sub-arc", "approaches its two defining crossings"). The "tidy" commit must not reword them.

## 3. What `plotting.py` reads from a `Hole` (verified: grep and Read of :237-300)

`plot_holes` reads `Hole.origin` (:268, :280, :284), `Hole.coords` (:285, :298) and `Hole.iterate` (:294-297). It reads neither `bridge_side` nor `openings`, so the plan's wording is wrong. The conclusion still holds: no `src` reader of the two fields exists outside `StablePartition.py`, and the rubric-3 grep is achievable.

## 4. The two invariant checks (verified: Read)

Both are assertions the algorithm needs:
- `check_bridge_rows_consistent` is I2 (equal crossing signs on one stable branch);
- `check_holes_share_bridge_side` is I1 (one side per origin, flipped by parity under an orientation-reversing map).

`Trellis.punch_holes` runs them at `Trellis.py:1271-1275`. Each is a two-line assert over a moved name:
- **I2** imports `bridge_row_violation`. The decision uses signs only (`rows_of_bridge`, `stable_key_at`), but the failure message appends `bridge_row_diagnosis(...).describe()`, which is geometry.
- **I1** imports `bridge_side_violations`, which goes through `side_diagnoses`, `_orbit_side` and `SideDiagnosis.describe`. The decision uses the orbit only.

Two consequences:
- **(B1) Circular import.** If the checks stay in `StablePartition.py` with top-level imports from `diagnostics`, and `diagnostics` imports `row_of_end` / `rows_of_bridge` from `StablePartition` (as rubric 2 allows), the two modules import each other. The plan's "~45 lines" counts only the two thin wrappers.
- **(B2) The map on the I1 path.** Step 2 has `side_diagnoses(trellis)` compute the measured side via `_map_backward`. I1 calls `side_diagnoses`, so the map would run on every `punch_holes`, including `diagnose=False`, inside the assertion path. The rubric-4 grep would still pass, but the algorithm would still be calling the map through `diagnostics`.

## 5. `_map_backward` also places holes (verified)

`propagate_reference_holes` maps the reference hole's coordinates back one step at a time (:502, :516) and hands the result, `carried`, to `_punch_in_bridge`. There it is used twice:
- for the measured side (:1605-1607);
- **for placement** (:1585-1603): the hole is nudged from the image-arc midpoint toward `carried`, with the step capped by the distance to it.

Without the map, the existing fallback places the hole instead (:1601/:1603, `_step_into(midpoint, stable_mid, 0.05, chord)`). I counted propagated holes on the wrong side of their image sub-arc, read with `_bridge_side_of(..., span)`:

| case | propagated | wrong side today | wrong side with the no-map fallback |
|---|---|---|---|
| k10 | 2 | 0 | 1 |
| k28_one_blast | 2 | 0 | 1 |
| k28_two_blasts | 3 | 0 | 1 |
| p3 | 6 | 0 | 3 |
| nested / nested_inner_blasts | 9 / 9 | 0 / 0 | 4 / 4 |
| inversion / one_blast / two_blasts | 4 / 6 / 6 | 0 / 0 / 0 | 2 / 1 / 1 |
| inversion_period2 | 8 | 4 | 5 |
| heteroclinic | 69 | 13 | 29 |

The largest move is 0.33 (k28 one blast).

So keeping rubric 4 ("no map in StablePartition") changes `Hole.coords` for propagated holes, and "Behaviour: none" is false. The fingerprint (`tests/helpers/law_tier.py:58-87`) does not include `coords`, so rubric 5 would not catch it; the plots would.

This also affects B4b. Its step 4 says `hole_point` "reproduces today's placement (… propagated: image sub-arc midpoint nudged toward the chord)". That description is the no-map fallback, not today's placement.

## 6. `side_source` as a function (verified on the 11 law cases and `p3_coarse_blasts`)

Today the field is `"signs"` exactly when `hole.pair is not None` (:377 vs :1617). The plan's rule, `iterate == 0 → "signs"`, gives the same answer on every hole probed: no direct hole has `iterate != 0`. The two rules differ only for a direct backward-iterate pair, which B4b asserts never exists.

## 7. The measured side can be computed on demand (verified)

For a propagated hole `h`, let `ref` be its origin's reference hole (`pair is not None`, `iterate == 0`) and `rb = bridge_for_pair(trellis, ref.pair)`:

```
carried  = map_inv applied (-h.iterate) times to ref.coords
span     = sorted(trellis.image_cdist(i, h.iterate, "unstable", branch=rb.manifold_key)[1] for i in rb.id)
measured = _bridge_side_of(trellis, trellis.bridge_between(*h.bounding_ids), carried, span)
```

- **Result.** This equals the stored `measured_side` on all 128 propagated holes of the 11 law cases. It also matches all 80 of `p3_coarse_blasts`, where 20 disagree with `bridge_side`.
- **Key.** `(origin, iterate)` is unique among propagated holes, by the propagation's `seen` set (80/80 unique on coarse).
- **Caution 1.** This replay relies on the reference hole's `coords`, which B4b deletes. After B4b the diagnostics need the reference placement from somewhere: `plotting.hole_point` (topology importing plotting) or a copy of it.
- **Caution 2.** The replay ignores the fallback start at :458-495, which never runs on the builds probed.

## 8. Baseline for rubric 5 (verified)

Script in the session scratchpad (`b4a_baseline.py`). It ran on `session.topology().trellis`, with hole bridges selected the way `Trellis.punch_holes` selects them:

| case | holes | propagated | hole bridges | `geometry_flags` (rows, sides) | `side_diagnoses` |
|---|---|---|---|---|---|
| k10 | 3 | 2 | 3 | ([], []) | [] |
| k28_one_blast | 4 | 2 | 3 | ([], []) | [] |
| k28_two_blasts | 4 | 3 | 4 | ([], []) | [] |
| p3 | 8 | 6 | 5 | ([], []) | [] |
| nested | 12 | 9 | 9 | ([], []) | [] |
| nested_inner_blasts | 12 | 9 | 9 | ([], []) | [] |
| inversion | 6 | 4 | 4 | ([], []) | [] |
| inversion_one_blast | 7 | 6 | 7 | ([], []) | [] |
| inversion_two_blasts | 8 | 6 | 8 | ([], []) | [] |
| inversion_period2 | 10 | 8 | 6 | ([], []) | [] |
| heteroclinic | 73 | 69 | 50 | ([], []) | [] |

- **All 11 cases.** Every case is orientation-preserving, and no propagated hole has `measured_side is None`.
- **The only non-empty baseline.** `p3_coarse_blasts` (in `BUILDERS`) has 84 holes on both the per-fixed-point and the topology trellis. It shows row flags on bridges `(41, 42)` and `(42, 43)`, 20 measured-side flags and 20 `side_diagnoses`.
- **Why rubric 5 is vacuous.** On the 11 cases, "identical on base and branch" holds even if the branch computes no measured side at all.
- **Tests run.** The touched subset passes at `HEAD`: 100 passed in 5.2 s (`test_topology_invariants`, `test_production_checks`, `test_stable_partition`, `test_plot_topology`, `test_coarse_blasts_rows_by_signs`, `test_exact_ownership`). I did not run the full suite.

## 9. Where the plan section differs from the code

1. **Name count.** "18 names moved" includes the two checks step 3 keeps. In effect 16 names move, plus `_map_backward`, and rubric 2's "holds the 18 names" should say 16.
2. **Line target.** A verbatim move leaves `StablePartition.py` at about 1,348 lines (1,947 − 599), not ≤ 1,100. About 1,090 is the figure after B4b, which removes roughly another 170 + 90.
3. **Placement.** `_map_backward` places holes (§5), so rubric 4 and "Behaviour: none" cannot both hold in B4a.
4. **Imports.** There is a circular import, and the map lands on the I1 path (§4 B1, B2).
5. **Signature change.** `side_diagnoses(trellis)` replaces the signature `(holes, orientation_preserving=)` used today by the I1 chain, by `test_topology_invariants.py` (5 tests) and, through the check, by `test_open_p3_deep_runs.py`.
6. **File list.** Add `test_topology_invariants.py`, `test_production_checks.py` and `test_coarse_blasts_rows_by_signs.py`, for 11 files. `laws.py` changes logic at :1209, not only imports.
7. **`plotting.py`.** It reads `origin`, `coords` and `iterate`, not `bridge_side` or `openings`.
8. **Rubric 5.** It is vacuous on the 11 law cases (§8).

## 10. Recommendation (the smallest changes that make the slice hold)

1. **Keep `_map_backward` and the carried placement in `StablePartition.py` until B4b.**
   - In B4a, remove only the `measured_side` reading (:1605-1607) and the two `Hole` fields, and move rubric 4 to B4b, which owns placement.
   - `diagnostics` gets its own few-line inverse-map helper for the measured side.
   - B4b step 4 should then state what propagated placement becomes, since reproducing today's placement needs the map.
   - If the author wants rubric 4 now instead: accept the fallback placement and record it as the one behaviour change (propagated hole coordinates; up to 29 of 69 heteroclinic holes drawn on the wrong side of their sub-arc).
2. **Keep both checks in `StablePartition.py` and import their cores lazily.**
   - Put `from .diagnostics import bridge_row_violation` / `bridge_side_violations` inside the function body, as `Trellis.punch_holes` already does.
   - Keep `side_diagnoses(holes, *, orientation_preserving=True, measured=None)`. Here `measured` is an optional `{(origin, iterate): Side}`, built by a new `diagnostics.measured_sides(trellis)` (§7).
   - `geometry_flags(trellis, bridges, holes)` builds that mapping. Only the flag path then calls the map, and every existing call signature stays.
   - Alternative: keep an I2 core that reads signs only and an I1 core that reads the orbit only in StablePartition (about 60 lines). That rewords the messages `stress/.../analysis.py` matches on.
3. **Define `side_source(hole)` in `diagnostics` as `"signs" if hole.pair is not None else "carried"`.** This is exactly today's assignment; the plan's `iterate == 0` rule agrees on every hole probed. `laws.py:1209` then calls it or drops the clause.
4. **Widen the file list by the three tests in §2.** They need:
   - the logger name `"tanglepack.topology.diagnostics"` for the moved WARNING/INFO;
   - the `_hole` helper without the two fields;
   - the measured-side tests driven through `measured_sides`.
5. **Make rubric 5 able to fail.** Compare, base vs branch:
   - the measured side of every propagated hole (base: `hole.measured_side`; branch: `measured_sides(trellis)`);
   - the full `bridge_row_diagnosis` of every hole bridge.

   Do this on the 11 law cases plus `p3_coarse_blasts`, the one non-empty baseline (2 row flags, 20 measured flags).
6. **Correct rubric 2** to 16 names and a file of about 1,350 lines or fewer.

## Out of scope (noticed, not for B4a)

- **Today's placement is not always on the right side.** By the same nearest-segment reading it is off its sub-arc on inversion_period2 (4 of 8) and heteroclinic (13 of 69). That affects plots only.
- **A possible B4b placement (believed, not verified).** A no-map nudge toward `bridge_side`, with the normal taken along the bridge's cdist direction, falls in the same half-plane as today's placement on 9 of the 11 cases, but not on inversion_period2 (4/8) or heteroclinic (64/69). It may be worth a look in B4b.
- **Stale docstring.** `numerics/geometry.py:7` names `topology.StablePartition._build_oriented_bridge_polyline`, which will be stale after the move; numerics should not name topology internals anyway.
- **Archived script.** `scripts/refactor_proof/hole-side-and-own-blast/hole_sides.py` will break.

## Knowledge files updated

- `.refactor/knowledge/topology.md`: the Map's line numbers are now at `ae18859`, plus a new B4a section.
- `.refactor/knowledge/testing.md`: a new B4a section listing every test that reads the moved names or fields, or pins the logger name.
