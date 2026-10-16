# Findings — B9 `session-three-concerns` (analyst, scoping run, 2026-10-15, HEAD = `517aa01`)

**Tree state.**
- `HEAD` = `517aa01` (main, B10 merged).
- The working tree is dirty, but only in documents and run directories:
  - Modified: `CLAUDE.md`, `docs/alignment_plan_2026-10-09.md`, `docs/code_fidelity_audit_2026-10-09.md`, `.refactor/{decisions.md,knowledge/*}`.
  - Untracked: `.refactor/runs/*`, `artifacts/`, `figures/henon_period3_*`, `scripts/henon_period3_{k,zone}_sweep.py`.
- No tracked `src/` or `tests/` file is modified. An AST scan shows the two untracked sweep scripts call nothing this slice touches.
- The change branches from `HEAD` and should leave all of these alone.

**Kind of request.** This is a refactor (moving and deleting code), scoped before it is planned. The target is `loom/TangleSession.py` plus the four "element at a bridge end" implementations: `PartitionFamily.py:343`, `Trellis.py:1381`, `StablePartition.py:1059` and `BridgeClass.py:763`.

**Short answer.**
- About 860 of the 1,744 lines can leave the session:
  - 438 lines of plot delegates.
  - 180 lines of cartoon geometry.
  - 204 lines of fan-outs and `partition_element_for`.
  - 39 more lines of describe fan-outs, which the request did not list.
- **Nothing computed changes:**
  - The four element-at-end readings agree everywhere: 7,090 bridge ends on their row, 14,180 (end, side) pairs and 1,776 iterated ends, with 0 disagreements across 63 items.
  - The cartoon geometry needs no session pip memory.
  - Building the zones from `Trellis.pips` gives bit-identical figures on every case that stores a zone and on every figure script. It changes the circle cartoon only where a pip exists but no zone was stored (the law case `inversion_period2`). That is an author question.
- **Recommendation:** split the slice. B9a takes the session's three concerns; B9b folds the element-at-end implementations into one.

---

## 1. Every `TangleSession` method mapped to a concern

Line ranges come from the AST at `517aa01` and include decorators and docstrings. The last number is the method's line count.

**Concern 1 — facade proper (stays): 700 lines.**
- Construction and delegation: `__init__` 96-146 (51), `__getattr__` 150-159 (10), the three property/setter pairs 166-191 (21), `fixed_points` 195-198 (4).
- Snapshots and topology: `trellis` 200-254 (55), `arrangement` 256-287 (32), `topology` 289-337 (49).
- Stage readers: `bridge_classes` 339-368 (30), `describe_bridge_classes` 370-384 (15), `homotopy_partition` 386-400 (15), `minimal_trellis` 402-422 (21), `iterated_partition` 424-442 (19), `describe_iterated_partition` 444-456 (13), `dual_graph` 458-484 (27), `symbolic_dynamics` 486-511 (26), `describe_symbolic_dynamics` 513-528 (16).
- Private helpers: `_is_stale` 947-949 (3), `_selection` 951-955 (5), `_resolve_fixed_points` 957-965 (9; `Blast.py:309` also reads it).
- Pip memory:
  - `classify_strong_pips` 1049-1091 (43). It has the shape of a fan-out, but it is the pip memory.
  - `_restore_strong_pip` 1093-1109 (17), `set_strong_pip` 1111-1117 (7), `forget_strong_pip` 1119-1124 (6), `strong_pip` 1126-1128 (3), `chosen_strong_pip` 1130-1133 (4), `strong_pip_candidates` 1135-1137 (3).
- Zones and blasting: `resonance_zone` 1483-1537 (55), `add_resonance_zones` 1539-1576 (38), `_bridge_test_point` 1580-1590 (11), `classify_bridge` 1592-1616 (25), `blast_zone` 1620-1686 (67).

**Concern 2 — plot delegates (leave): 438 lines.**
- Forwards to `topology.plotting`:
  - `plot_transition_graph` 530-560 (31), `plot_walk` 562-594 (33), `plot_bridges_by_class` 596-628 (33), `plot_itinerary_table` 630-661 (32).
  - `plot_dual_graph_cartoon` 663-707 (45). With `shape="circle"` it also fills `interior_side`, `parents`, `arc_order` and `enclosures` from concern 3 (:697-706).
  - `plot_dual_graph` 891-917 (27), `plot_minimal_trellis` 919-945 (27).
- `_fanout_plot` forwards to `Trellis.plot_*`:
  - `plot_strong_pip_candidates` 1139-1177 (39) and `plot_strong_pip` 1179-1214 (36). Their `prepare` hook classifies through `self.classify_strong_pips`, i.e. through the pip memory.
  - `plot_pseudoneighbors` 1244-1275 (32); its `prepare` runs `compute_pseudoneighbors`.
  - `plot_holes` 1277-1300 (24), `plot_stable_partition` 1414-1437 (24).
- `plot_resonance_zones` 1688-1742 (55) has its own drawing body of about 30 lines. It reads `self.resonance_zones` and `self.workbench`; its `boundary_polygon(workbench)` fallback is dead for a stored zone.
- Twelve of the thirteen just forward to existing functions and are deleted. Only `plot_resonance_zones` has a body that needs a new home.

**Concern 3 — cartoon geometry (leave): 180 lines.** `cartoon_enclosures` 709-779 (71) and `cartoon_zones` 781-889 (109).

**Concern 4 — fan-outs and `partition_element_for` (deleted): 204 lines.**
- `_fanout_call` 967-997 (31), `_fanout_plot` 999-1039 (41).
- `compute_pseudoneighbors` 1218-1242 (25), `punch_holes` 1304-1333 (30), `partition_stable_manifold` 1335-1360 (26).
- `partition_element_for` 1362-1412 (51).
- **Not in the request, same shape:** `describe_holes` 1439-1452 (14), `describe_stable_partitions` 1454-1467 (14), `_describe` 1469-1479 (11).

**What is left.** The other 183 lines are the module header and imports (1-54), the class docstring (55-95), section comments and blank lines.

**Simulated size after removal (deleting the ranges and their trailing blank lines):**
- Without concerns 2-4: 898 lines.
- Also without the describe fan-outs: 856 lines.
- Pruning the now-unused imports (`plotting`, `Callable`, `endpoint_index`, `ring_order`, `CartoonLayout`, `Walk`, `Table`, `Line2D`, `networkx`) saves about 8 more.
- So a cap of ≤ 900 (or ≤ 870 with the describe fan-outs gone) is honest.
- Keeping the seven symbolic-dynamics delegates with their docstrings adds about 230 lines (about 1,080), which breaks the ≤ 1,000 cap unless each is cut to one line.

---

## 2. External callers of every leaving method, and what each becomes

*Verified* by:
- an AST attribute scan over `src tests scripts stress notebooks` (receivers resolved);
- a grep for `getattr` strings;
- a grep over every `.ipynb`.

`stress/` and `notebooks/` have no callers. `scripts/tangle_workbench_test.ipynb` has one (`plot_resonance_zones`).

### Plot delegates

| Caller | Becomes |
|---|---|
| `scripts/symbolic_figures.py:181` `session.plot_strong_pip(fp, ax=ax, s=120, …)` | `session.trellis(fp).plot_strong_pip(ax=ax, …)` (the pip is always set after `topology()`) |
| `scripts/symbolic_figures.py:225` `session.plot_dual_graph_cartoon(dynamics, …, shape=…)` | `plotting.plot_dual_graph_cartoon(session.dual_graph(), dynamics, …)`, with the circle inputs from the moved functions (§3) |
| `scripts/henon_pseudoneighbors.py:95-97` | `trellis.plot_holes` / `plot_pseudoneighbors` / `plot_strong_pip` (`trellis` is already a local at :89) |
| `scripts/henon_pseudoneighbors_nested.py:97-98` | a loop over `(fp3, fp1)`: `session.trellis(fp).plot_pseudoneighbors(s=40)` and `.plot_holes()` (both trellises are partitioned by hand at :83-89) |
| `scripts/henon_strong_pip_nested.py:81-82`, plus the docstring at :6 and the comment at :80 | a loop over the per-fp trellises (`classify_strong_pips()` already ran at :67) |
| `scripts/henon_minimal_dual_graph.py:118-119` | `session.trellis(fp).plot_holes` / `.plot_strong_pip` |
| `scripts/henon_minimal_dual_graph.py:149` | `plotting.plot_dual_graph(dual, ax=, clip_to_arcs=False)` |
| `scripts/henon_dual_graph.py:244` | `plotting.plot_dual_graph(…)` |
| `scripts/henon_k10_dual_graph.py:219` (dual graph), `:225`, `:226`, `:230`, `:239` | `plotting.plot_dual_graph`; `session.trellis(fp).plot_*` |
| `scripts/henon_bridge_classes.py:253`, `:267`, `:332` | `session.trellis(fp).plot_strong_pip` / `.plot_stable_partition` |
| `session.plot_resonance_zones(alpha=0.3)` in `henon_resonance_zone.py:92`, `henon_bridge_classification.py:145`, `henon_blast_period_3.py:136` (dead at :95) and `tangle_workbench_test.ipynb` | the new zones plotter (Q3), called with `session.resonance_zones.values()` |
| `tests/plotting/test_plot_smoke.py:59-71` | delete (fan-out only; the trellis plotters are covered at :48-56) |
| `tests/plotting/test_plot_smoke.py:74-80` | delete (tests the `prepare` hook only) |
| `tests/plotting/test_plot_smoke.py:277-308` | delete (each wrapped function has its own smoke test at :116-245) |
| `tests/plotting/test_plot_smoke.py:311-320` | rewrite on `plotting.plot_dual_graph_cartoon` with the same three assertions |
| `tests/plotting/test_plot_topology.py:400`, `:422` | `plotting.plot_dual_graph_cartoon` plus the moved functions |
| `tests/facade/test_session_bridge_classes.py:151-168` | exists only for the deleted `prepare → classify` path; delete (Q5) |

`Trellis.strong_pip_candidates` (`Trellis.py:622-637`) classifies WITHOUT choosing a pip, so `Trellis.plot_*` never touches the pip memory.

### Cartoon geometry

| Caller | Becomes |
|---|---|
| `TangleSession.py:699`, `:704` | goes with the delegate |
| `tests/plotting/test_plot_topology.py:391`, `:412`, `:415`, `:419` | the moved functions |

At `:419` the sub-selection `alone` needs the CHILD fixed point, which is not in `alone.trellis.fixed_points`. A moved `cartoon_enclosures` must therefore take the fixed points (or the parents as objects), not just a trellis.

The docstrings at `plotting.py:6`, `:2519` and `:2896` cite the session and need a wording change only.

### Fan-outs and `partition_element_for`

| Caller | Becomes |
|---|---|
| `tests/facade/test_session_fanouts.py:17-44`, `:47-65` | keep the `classify_strong_pips` halves (:23-27, :51-53); delete the rest |
| `tests/facade/test_session_fanouts.py:68-83` | delete if the describe fan-outs go (Q4) |
| `tests/facade/test_session_bridge_classes.py:36-38` (`_repartition_at`) | delete; they are no-ops for the assertion because the product recomputes (`Topology.py:106-108`) |
| `tests/facade/test_session_caches.py:135-137` | `trellis.compute_pseudoneighbors()` / `punch_holes()` / `partition_stable_manifold()` on the local `trellis` (:127); this is exactly what the one-fp fan-out ran |
| `tests/regression/test_nested_own_blast.py:78-80` | delete (`_tangle_words` reads `symbolic_dynamics()`, which recomputes) |
| `tests/unit/topology/test_arrangement.py:316-318` | `session.topology(fp)` |
| `scripts/henon_pseudoneighbors.py:83` | `trellis.compute_pseudoneighbors(verbose=True)` |
| `scripts/henon_dual_graph.py:83-85` (one fp) | `session.topology()` |
| `scripts/henon_dual_graph.py:134-136` (NESTED) | `for fp in fps: session.topology(fp)` |
| `scripts/henon_k10_dual_graph.py:186-188` | `session.topology(fp)` |
| `scripts/henon_k10_dual_graph.py:196` | `session.trellis(fp).describe_stable_partitions()` |
| `partition_element_for` at `test_partition_elements.py:162`, `:186` | see §4 |

For the nested script, `topology()` of the whole selection would leave the per-fp trellises empty (`Topology.py:103-110`). Calling `topology(fp)` per fixed point fills them exactly as the fan-out did.

**Totals for one branch.**
- src, 9-11 files: `TangleSession.py`; a new module for the cartoon geometry and the zones plotter; `ResonanceZone.py` only under option (b) of §3; `plotting.py` (docstrings); `Trellis.py`, `StablePartition.py`, `PartitionFamily.py`, `BridgeClass.py`, `TopologyResults.py`; and `loom/__init__.py` if the new functions are exported.
- tests, 8 files: `test_plot_smoke.py`, `test_plot_topology.py`, `test_session_fanouts.py`, `test_session_bridge_classes.py`, `test_session_caches.py`, `test_nested_own_blast.py`, `test_arrangement.py`, `test_partition_elements.py`.
- scripts: 11 `.py` files plus 1 notebook.

That is about 30 files and about 60 call sites.

---

## 3. Cartoon geometry: what it reads, and whether zones can come from `Trellis.pips`

*Verified* by reading :709-889 and by running prototypes.

**What `cartoon_zones` reads:**
- `self.resonance_zones.values()`, sorted by area. From each zone: `boundary_vertices`, `boundary_intersections` (`manifold_a_key`, `manifold_b_key`, `stable_cdist`), `contains_point`, `fixed_point.coordinates`.
- `workbench.manifolds`, through `arc_polyline` and `key in …`.
- The registry's `cdist_tol`.
- loom `ring_order` (`ResonanceZone.py:417`), `signed_polygon_area` and topology's `ZoneArcOrder`.

It reads no pip memory, no trellis and no topology.

**What `cartoon_enclosures` reads:**
- `workbench.fixed_points`, which can include fixed points outside the selection.
- `workbench.manifolds`. This is the same live dict as `trellis.manifolds` (`Trellis.py:203`, `:486`).
- The trellis it is given: `trellis.registry` and `trellis.bridge_between`.
- `parents`.

It also reads no pip memory.

**Where they can go under the dependency rule:**
- `cartoon_enclosures` can move to `topology.plotting` as a function of `(trellis, bridge_ids, parents, fixed_points)`.
- `cartoon_zones` as written takes loom `ResonanceZone` objects and loom `ring_order`, so it cannot go to topology. It goes verbatim to a loom module (e.g. `loom/cartoon.py`) as a function of `(zones, manifolds, cdist_tol)`, with zero figure change by construction.
- `plot_resonance_zones` has the same constraint.

**Zones from pips: prototype.** Script: scratchpad `b9/probe/zones.py`. Method:
- Build one zone per entry of `topology().trellis.pips`.
- Take its cuts from `registry.iterate_orbit(pip, k_value)`.
- Build the polygon LIVE with the same boundary builder.
- Rerun `cartoon_zones` on those zones and compare.

It ran on the 11 law cases, the 3 `henon_symbolic_itineraries` cases, the 4 `henon_cases` recipes, and the k28, nested_outer and inversion series at every step.

**Results:**
- **Every case that stores a zone (41 of 43):** same pip, `np.array_equal` polygons, identical `cartoon_zones` output.
- **The two cases with no stored zone gain geometry:**
  - `k10`: nodes and artists are identical, because the implied interior side is "right", which is the default.
  - `inversion_period2`: the circle cartoon **changes** (the ring order and the anchor angles). No figure script draws it, and no plotting test draws it as a circle.
- **By code reading, not in the corpus:** a zone stored at a pip other than the chosen one (the B8 WARNING path, or `set_strong_pip` after a zone) would be drawn at the chosen pip.

Doing this in topology means moving `ResonanceZone._build_boundary` (:164-253) and `ring_cycles` / `ring_order` (:383-431) into topology, which puts `ResonanceZone.py` into the slice. This is Q2.

---

## 4. The four "element at a bridge end", what each reads, and the probe

| Implementation | Family it reads | Trellis it reads | Side | Domain | Returns / on missing |
|---|---|---|---|---|---|
| `PartitionFamily.owner_at_end` :343 | its own: homotopy or iterated | the one passed | the row (`row_of_end`) | any `BridgeId` pair, including image pairs (`DualWalk.py:1366`) | `ElementRef`; `ValueError` |
| `Trellis.element_for` :1381 | `self.stable_partitions` (homotopy) | self | **any side, passed in** | any `BridgeId` | `Optional[int]`; `None` |
| `StablePartition._elements_at_bridge` :1059 | the result being built | `trellis.bridges` (iterated: `minimal.sparse`) | each side | that trellis's bridges only | dict; asserts |
| `BridgeClass._element_at_end` :763 | the homotopy `partitions` list | the one passed | the row | trellis bridge ids | `ElementRef`; `ValueError` / `AssertionError` |

- `_elements_at_bridge` is stored as `elements_at_bridge` at `:536` (homotopy) and `PartitionFamily.py:898` (iterated).
- `_element_at_end` is called at `BridgeClass.py:597-598` and, through `ends`, from `_non_loop_parent` at :740.
- The law tier has its own independent oracle, `laws._element_ends` (`tests/helpers/laws.py:1868`).

**Probe.** Script: scratchpad `b9/probe/probe.py` + `run.sh`; output in `all_base.jsonl`. It ran at `HEAD` on the 63 items of the B10 probe:
- 16 builds: the 11 law cases plus `p3_16`, `p3_6_blasts`, `k28_8`, `nested_inner6`, `p3_15`.
- 47 series steps: `k28` 0-8, `p3` 0-9, `inversion` 0-4, `k32` 0-10, `nested_outer` 0-11.

All 63 run on `session.topology()`.

**Three checks:**
1. On each end's row, all four implementations plus the law oracle plus `partition_element_for` agree, and the ref has the `stable_key_at` branch and the `row_of_end` side: **7,090 of 7,090**.
2. On both sides, `owner_of_intersection` = `element_for` = `elements_at_bridge`: **14,180 of 14,180**.
3. Iterated: `owner_at_end` = the iterated results' `elements_at_bridge`: **1,776 of 1,776**. Reading on the full and on the sparse trellis also agree.

There were 0 exceptions. **No disagreement, so there is no finding for the author.** *(verified)*

**Before folding, the planner should know five things:**
1. **`element_for` answers on ANY side**, not just the row.
   - Its only src caller dies with `partition_element_for`; its other readers are tests.
   - `test_partition_elements.py` checks both sides (:45-78, :99-119, :122-162).
   - Rewriting those via `owner_at_end` would weaken them. Rewriting them as `family.owner_of_intersection(id, side, trellis.stable_key_at(id, trellis.unstable_key_of(bid)))` keeps their strength.
2. **`elements_at_bridge` is written but never read in src.**
   - Its readers are tests only: `test_partition_elements.py:77`, `:87`, `:158`.
   - The field (`TopologyResults.py:415-424`, `:442`), both writers and the import at `PartitionFamily.py:140` can go together, with the test at :81-96.
3. **Import cycle.** `PartitionFamily.py:136` imports from `BridgeClass`, so `BridgeClass` cannot import `PartitionFamily` at module level. The shared core belongs beside `row_of_end` in `StablePartition.py`, or in a local import. This is a planner choice.
4. **Error-path exceptions differ.** `_element_at_end` raises `AssertionError` for an unowned crossing where `owner_of_intersection` raises `ValueError`. No probe item reaches either path.
5. **Law-tier reads.** The law tier reads none of the four. It reads `row_of_end` and `element_of_intersection` directly (`laws.py:75`, `:1177-1197`, `:1323`, `:1874`, `:1979`, `:2200`, `:2226`, `:2252`), and B9 touches neither.

**`test_partition_elements.py` (201 lines, 6 tests):**
- :45-78 — rewrite with `owner_of_intersection` and drop the clause at :77.
- :81-96 — delete.
- :99-119 — rewrite (both sides).
- :122-162 — keep the mechanism, rewrite, and drop :158 and :162.
- :165-191 — delete. Optionally replace it with "the selection trellis answers both tangles".
- :194-201 — `owner_at_end(…, "third")` raises `ValueError` through `endpoint_index`.

`test_dual_walk.py:877` binds `PartitionFamily.owner_at_end` onto a fake. It still works as long as `owner_at_end` calls `self.owner_of_intersection`.

---

## 5. What the law tier reads of anything that moves or folds

*Verified* by grep over `tests/helpers/*.py`, `tests/cases.py`, `tests/invariants/*.py` and `tests/conftest.py`.

- **Moving or folding pieces: none.** The law tier reads no plot delegate, `cartoon_*`, `_fanout_*`, session stage call or `describe_*`, nor `partition_element_for`, `element_for`, `elements_at_bridge`, `_element_at_end` or `owner_at_end`.
- **Reads that survive unchanged:**
  - `session.resonance_zones` and `add_resonance_zones`: `cases.py:536`, `:666-667`, `:744`; `helpers/zones.py:45-46`; `conftest.py:279`.
  - `Trellis` stage methods: `helpers/direct.py:90-92`.
  - `row_of_end` and `element_of_intersection`.

---

## 6. Plotting tests and a base-vs-branch figure check

**Tests touching the delegates or the cartoon:**
- `test_plot_smoke.py:59-71`, `:74-80`, `:277-308`, `:311-320`.
- `test_plot_topology.py:381-405`, `:408-431`.
- `test_class_color_families.py`: none.
- Outside `tests/plotting/`: `test_session_bridge_classes.py:151-168`.

**Figure check.** *(verified)*
- In a `git archive HEAD` export, the figure scripts regenerate 42 byte-identical PNGs run to run: `k10`, `k28_blast1/2`, `p3`, `nested_p1_p3`, `heteroclinic_p1_p3` and `inversion`, six panels each.
- `k10` and `k28` equal the committed PNGs. The committed `p3`, `nested`, `heteroclinic` and `inversion` PNGs are STALE: HEAD draws different ones.
- So compare base-regenerated against branch-regenerated, never against the committed files.
- Run times: `henon_symbolic_itineraries.py` 7 s, `_period3` 2 s, `_nested` 3 s, `_nested --heteroclinic` 14 s, `_inversion` 2 s. All of them can be run in full; none needs to fall back to `py_compile`.

**Suggested rubric:**
1. sha256 equality of the 42 PNGs, base vs branch, both regenerated outside the repo.
2. A scratch script drawing the circle cartoon, `plot_resonance_zones` and every rewritten fan-out site on the `tests/plotting` fixtures. Compare the sha256 of `fig.canvas.buffer_rgba()` plus an artist dump; scratchpad `b9/probe/zones_fig.py` is a working example.
3. Smoke-run all 11 touched scripts.

---

## 7. State at HEAD

- **Suite** (`MPLBACKEND=Agg PYTHONPATH=src env/bin/python -m pytest -q -rxX -p no:cacheprovider`): 1215 passed, 47 skipped, 5 deselected, 1 xfailed (C3), no XPASS, 204 s.
- **Scripts:** every script touched here exits 0 in an export, except `henon_blast_period_3.py`. That one is still dead at :95 (`ValueError: Intersection 12 is not a strong-pip candidate. Candidates: [31, 32, 33]`). Its :136 `plot_resonance_zones` call is never reached but must still be rewritten (or reported).
- **63-item probe:** the "no number changes" check applies unchanged. The only change on the words path is the element-at-end fold, which the agreement probe shows is the identity on the whole corpus.

---

## 8. For the planner

**Line counts:**

| Part | Lines |
|---|---|
| Concern 2 (plot delegates) | 438 (about 30 of body move) |
| Concern 3 (cartoon geometry) | 180, moved |
| Concern 4 (fan-outs, `partition_element_for`) | 204, plus 39 describe fan-outs |
| `TangleSession.py` after the slice | 898 or 856, minus about 8 for imports |
| Element-at-end fold | about 130 src lines deleted: `element_for` 55, `_elements_at_bridge` 37 plus about 12 for the field and its writers, `_element_at_end` 31 down to a one-line call |
| `test_partition_elements.py` | about 200 lines rewritten |

The plan's "−750 (≈250 moved)" holds: roughly −1,000 deleted and +230 moved.

**Recommendation: split.**
- **B9a `session-three-concerns`:**
  - Scope: `TangleSession.py`; the new loom module for the cartoon geometry and the zones plotter; `plotting.py` docstrings; 11 scripts plus the notebook; 7 test files.
  - It deletes `partition_element_for`, and in `test_partition_elements.py` only :162 and :165-191.
  - It ends with `TangleSession.py` ≤ 900.
  - Proof: the 42 PNG hashes, the 63-item probe and the suite.
- **B9b `one-element-at-end`:**
  - Scope: `Trellis.py`, `StablePartition.py`, `PartitionFamily.py`, `BridgeClass.py`, `TopologyResults.py`, and the rest of `test_partition_elements.py`.
  - It is small, but it touches the class rule (`bridge_classes`' `ends`) and an import cycle, which needs a different kind of review.
  - Proof: the agreement probe rerun on the branch, plus the 63-item probe.

The split falls at `partition_element_for`: B9a deletes it, leaving `Trellis.element_for` read only by tests, and B9b deletes `element_for`. As one branch it is about 30 files across loom, topology, tests and scripts.

**Questions only the author can answer** (suggested defaults in brackets):
1. **Do thin plot delegates stay on the session for scripts?** [none, per the plan; delegates with docstrings break the ≤ 1,000 cap]
2. **Zones from the stored zones, or from pips?**
   - (a) Move `cartoon_zones` verbatim to a loom module that reads the stored zones: zero figure change, `ResonanceZone.py` untouched.
   - (b) Compute it in topology from `trellis.pips`: identical wherever a zone is stored and on every figure script, but it changes `inversion_period2`'s circle cartoon, follows the chosen pip rather than a stored zone, and moves the boundary builder plus `ring_cycles` / `ring_order` into topology.
   - [(a)]
3. **Where does `plot_resonance_zones` go?** [the same new loom module, taking zones]
4. **Delete the describe fan-outs too?** These are `describe_holes`, `describe_stable_partitions` and `_describe`, called from `test_session_fanouts.py:68-83` and `henon_k10_dual_graph.py:196`. [yes]
5. **`test_a_pip_trimmed_away_is_not_replaced_by_plotting` (`test_session_bridge_classes.py:151-168`):** delete it, or keep it rewritten on `Trellis.plot_*` as a guard that drawing never chooses a pip? [delete]
6. **Keep "element on any side" as API** (e.g. a `side=` on `owner_at_end`), or are `owner_of_intersection` plus `stable_key_at` enough? [no new API]

## Out of scope
- The committed higher-period figures are stale: HEAD draws different PNGs.
- `Trellis.plot_*` (`Trellis.py:1465-1575`) are themselves thin delegates (SIMP-11 / B12).
- The CLAUDE.md `TangleSession` row still lists the delegates and `cartoon_zones` (B12).
- `henon_blast_period_3.py` is still dead at :95.
- `classify_strong_pips` still returns a dict or a list depending on its argument (the fan-out shape); it is pip memory and stays.

Probe and prototype files are in /tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b9/probe/: probe.py, run.sh, all_base.jsonl, zones.py, zones_lib.py, zones_fig.py, zout/. Knowledge files updated: `.refactor/knowledge/loom.md`, `topology.md`, `testing.md`.
