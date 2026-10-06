# Test infrastructure audit (conftest, helpers, visualizations, pytest config, henon_cases)

## 0. Headline numbers (measured)

- 796 tests collected across 64 files and 18.3k lines. 583 of them use at least one fixture.
- I ran every file on its own with `--durations=0`. The test time adds up to **43.5 s**: 13.2 s in tests marked `slow` (81 tests) and 30.4 s in unmarked ones. With pytest startup per file the wall time is about 90 s. No single test takes more than 1.5 s.
- **Fixture builds are all cheap now.** The bulk-loaded rtree changed the cost picture, so the `slow` marker and several docstrings no longer match reality (details in §1 and §4).

## 1. Fixture inventory

Build cost was measured by calling each fixture body directly in one process. Usage comes from `pytest --fixtures-per-test`.

### conftest.py

| Fixture | Scope | Builds | Cost | Tests / files |
|---|---|---|---|---|
| `henon_map` / `henon_map_inverse` | function | return module globals (k=10) | 0 | 309 / 41 each |
| `henon_jacobian` | function | module global | 0 | 1 / 1 |
| `workbench` | function | empty `TangleWorkbench` k=10 | 0 | 127 / 26 |
| `fixed_point` | function | saddle [4,-4], oriented | <5 ms | 126 / 26 |
| `initialized` | function | + fundamental segments | <5 ms | 121 / 26 |
| `grown_unstable` | function | + 7 unstable steps | <5 ms | 9 / 5 |
| `grown_both` | function | 7 unstable + stable to turnaround | 10 ms | 85 / 21 |
| `small_tangle` | function | + `compute_intersections` | 10 ms | 28 / 8 |
| `henon_tangle_with_bridges` | function | grown_both + 2 more steps (9 total), intersect, trim, bridges (workbench) | 40 ms | 52 / 15 |
| `henon_inversion_initialized` | function | inversion saddle, seeded | <5 ms | 3 / 2 |
| `henon_inversion` | **module** | inversion: 6 u + 5 s steps, intersections | **1.24 s** (most expensive fixture) | 11 / 1 |
| `henon_p3_session` | function | nested p3 + p1 (p1 `area_cutoff=1e-4`), zones, no blast | 90 ms (34 crossings) | **83 / 24** |
| `k10_session` | function | session: 9 u steps, turnaround, intersect, trim, bridges, `infer_iterate_table` | 40 ms (8 crossings) | 147 / 14 |
| `k10_partitioned` | function | + classify / pseudoneighbors / holes / partition | 40 ms | 136 / 13 |
| `p3_partitioned` | function | henon_p3_session + partition | 110 ms | 31 / 11 |
| `k28_partitioned` | **session** | k=2.8, 10 steps, trim at f(q0), 1 blast, partition | 0.34 s | 6 / 6 |
| `k28_two_blasts_partitioned` | **session** | same with 2 blasts | 0.31 s | 2 / 2 |

### Local fixtures in test files

| Fixture | File | Scope | Notes |
|---|---|---|---|
| `k10_session` | test_session_trellis_cache.py | function | **Shadows the conftest name.** Same recipe minus `infer_iterate_table`. |
| `k10_session` | test_loom_blast_restore.py | function | **Shadows the conftest name.** Byte-identical to the one above. |
| `k10_zone_session` | test_resonance_zone_region.py | function | Same k10 recipe + zone. |
| `henon_session` | test_session_pseudoneighbors.py | function | Same k10 recipe, no infer. |
| `henon_session` / `populated_trellis` | test_topology_plotting.py | function | Same k10 recipe again; then partition. |
| `p3_built` / `nested_built` | test_higher_period_cartoon.py | session | `build_period3()` 70 ms and `build_nested()` 0.40 s, from henon_cases. |
| `p3_partitioned_c` | test_image_cdist.py | function | Wraps henon_p3_session. |
| `p3_punched` | test_stable_partition_period3.py | function | Wraps henon_p3_session. |
| `henon_partitioned` | test_partition_elements.py | function | Wraps henon_tangle_with_bridges. |
| `henon_with_holes` | test_stable_partition.py | function | Same wrapping. |
| `henon_punched` | test_topology_invariants.py | function | Same wrapping. |
| `k10_naming`, `k10_dynamics` | element_naming, symbolic_dynamics | function | Wrap k10_partitioned. |
| `stable_line`, `stable_line_with_manifold`, `two_branch_bridges` | test_stable_partition.py | function | Synthetic. |
| `layout` | test_symbolic_dynamics.py | function | Synthetic. |
| `key`, `landing_setup`, `chain_setup` | test_dual_walk.py | function | Synthetic. |

### Where the time really goes

- **Function-scoped rebuilds of `henon_p3_session`:** 83 builds × about 0.1 s ≈ 8–9 s, about 20 % of all test time.
- **k10 fixtures:** ~280 uses × 0.04 s ≈ 10 s.
- **A repeated "mutate then repartition" block** (grow one more step, then re-intersect, re-bridge and re-partition) costs 1.3–1.5 s each time. It is copy-pasted into three unmarked tests:
  - `test_session_bridge_classes.py::test_workbench_mutation_invalidates_the_cache`
  - `test_session_symbolic_dynamics.py::test_workbench_mutation_invalidates_the_cache`
  - `test_session_dual_graph.py::test_workbench_mutation_invalidates_the_caches`

  These are the three most expensive unmarked tests. They all assert the same generation-invalidates-cache fact, on three caches that share one cache key.
- **The inversion tangle is built twice:** once by the `henon_inversion` module fixture (1.24 s) and again inline in `test_bridge_identity.py::test_bridge_identity_on_the_inversion_saddle` (1.35 s, unmarked, same 6 u / 5 s recipe from `henon_inversion_initialized`).

## 2. Duplicated fixture-building code

**The k=10 "9 unstable, turnaround, intersect, trim, bridges" recipe exists 9 times:**
- conftest `k10_session`
- conftest `henon_tangle_with_bridges` (the workbench flavour, as 7+2)
- the local fixtures in test_session_trellis_cache, test_loom_blast_restore, test_resonance_zone_region, test_session_pseudoneighbors and test_topology_plotting
- inline in `test_session_bridge_classes.py:154` (the 10-step variant)
- inline in `numerics/test_grow_until.py:290` (the 7-step variant)

They differ only in whether `infer_iterate_table` runs, which is a silent behavioural difference between fixtures that look identical.

**Other duplicated pieces:**
- **`_define_zone(session, fp)`** is duplicated between test_session_trellis_cache.py:51 and test_loom_blast_restore.py:61. One docstring says "same reasoning" as the other.
- **The k28 fixtures:** `k28_partitioned` and `k28_two_blasts_partitioned` are about 45 lines each and differ only in the blast count (1 vs `range(2)`). Their pip-repin loop is the same logic as `henon_cases._repin`. There is no `build_k28(blasts=n)` in henon_cases, although the script `henon_bridge_classes.py` has the recipe.
- **Two different "nested p3" fixtures:**
  - conftest `henon_p3_session`: p1 `area_cutoff=1e-4`, 0 blasts.
  - `henon_cases.build_nested`: p1 `1e-7`, 2 outer blasts by default.

  Tests pin facts on both, and henon_cases' Dev Notes even cite the conftest fixture as its parameter source. The two `viz_*.py` scripts carry a third copy (10 u / 6 s / 7 u steps).
- **`henon_p3_session`'s `inner_zone` is misnamed.** It is `max(area)`, and I checked that it is the period-1 OUTER zone (area 21.06 against 0.54 for the period-3 zone). henon_cases calls the smallest zone `inner`. The fixture name contradicts the library's convention.
- **The dual-graph assembly `DualGraph(pieces.minimal, pieces.iterated, strong_pips=pieces.strong_pips)`** is repeated as a private helper in four files: `_k10_graph`/`_p3_graph` (test_dual_graph), `_k10_dual_graph` and `_k10_dynamics` (test_topology_plotting), and `_k10_dynamics` (test_higher_period_cartoon).
  - The two `_k10_dynamics` copies are identical, including the workaround `if any(entry.letter is None ...): pieces.table = session.bridge_classes([fp])`.
  - That workaround means `build_pieces` hands back an un-lettered table, which is an API smell the tests paper over.
- **Raw `construct_fixed_point` / `TangleSession(` / `TangleWorkbench(` builds** appear outside conftest in about 18 files. Most are legitimate (solver, GPU or batched-map tests that need custom maps). The ones listed above are not.
- **Private access in fixtures:** `session.workbench._man_machine.area_cutoff` is set in 4 test files, conftest and henon_cases. `minimal_helpers.build_pieces` calls the private `session._gathered_partitions()`.

## 3. Helpers that encode rules

### tests/invariants.py

**Good, physics-level checks:** spikes, iterate law, one-to-one, area along a chain. Issues:

- **`assert_cdist_monotonic` defaults to `strict=True`.** The module docstring and the cdist memory note both say ties are legitimate and strict is opt-in for low-stretch growth only.
  - `test_inversion_fixture.py:106` calls it with the default, so strict.
  - The inversion branch growth is the hard-refining case, so this default could produce false failures.
  - The default should be `strict=False`.
- **`assert_no_cdist_collision` encodes cdist injectivity**, which the memory note "cdist strict-monotonicity fix" says is NOT the correct invariant (the correct one is geometric no-spike). It is used only by `test_invariant_helpers.py`, on low-stretch growth. Low value; candidate for removal.
- **`assert_area_preserved_along_chain` only runs inside `test_invariant_helpers.py`.** That file is a test of the helper plus one real application. Area preservation is not checked on any blasted or period-3 fixture. If it is a core invariant, it is under-applied.
- **Cosmetic bug in `assert_no_geometric_spikes`:** its message computes `worst[1] / length_ratio` and labels it "= 30x the median". The median is never printed, so the arithmetic in the message is meaningless.
- **Unused exports:** `find_geometric_spikes` is used only inside invariants.py. `walk_nodes` is used by one outside file, `node_cdist` by none.

### tests/walk_helpers.py

Fakes used by 7 files: FakeFixedPoint in bridge_class, element_naming, arrangement_sparse, symbolic_dynamics, arrangement and dual_walk; FakeTrellis in dual_walk, symbolic_dynamics and higher_period_cartoon; FakeDual in dual_walk only. Concerns:

- **`FakeFixedPoint.advance_key` re-implements the ONE rule CLAUDE.md says lives only in `FixedPoint.advance_key`** (the branch flips on an orbit wrap with inversion). If the rule changes, the fake silently diverges and the synthetic tests keep passing.
- **`per_step_beta` is reduced to `beta` / `1/beta`.**
- **`FakeTrellis.image_cdist` re-implements the "table first, else `advance_key` + `per_step_beta`" rule.** `iterate` re-implements the inverse-table composition.
- **`FakeTrellis.intersection` returns a `SimpleNamespace` with `stable_cdist`, `crossing_sign` and `manifold_b_key`.** That couples the fake to Intersection field names.
- **`make_result` derives default owners by "the interval closed at an end owns it"**, which is a partition-ownership rule encoded in the test helper.
- **`FakeDual` hard-codes `StableNode` internals.** It builds the `edge_key` tuple layout `("stable", lo, hi, branch_key)` and `key=(edge_key, "both"|side)`, and sets `traversable` / `fundamental_span` directly. Any refactor of the StableNode identity breaks every synthetic walk test, while the science stays the same.

### tests/minimal_helpers.py

- **`build_pieces` duplicates the assembly in `TangleSession.minimal_trellis` / `iterated_partition`** (no caches, private `_gathered_partitions`). It is used in 11 files (37 call sites). Tests built on it assert facts about a hand-assembled pipeline, not the session API. If the session's assembly changes (for example, how holes or pips are gathered), the two can diverge with no test noticing.

### src/tanglepack/examples/henon_cases.py

- Clean builders, session-scoped in test_higher_period_cartoon.
- The Dev Notes pin fixture facts (`a -> b -> c -> a u^-1 w^-1`, "the default stays at 13"). They are also asserted in tests (`_P3_WORDS`, `test_p3_words_survive_four_blasts`). Those are raw fixture facts likely to change if the growth recipe changes.
- The module is shared library code that tests depend on. A change for figure-making (for example, the default `outer_blasts`) silently changes the test fixtures.

## 4. Markers (pyproject: `--strict-markers`, `slow`, `regression`)

**Nothing ever selects them.** There is no CI config, Makefile or doc that uses `-m slow`, `-m "not slow"` or `-m regression`, and there is no `.github`. Both markers are decorative.

**`slow` is inconsistent and stale.** It is used on 81 tests in 22 files.
- **Marked slow but cheap:** 20 tests take under 0.1 s (setup included), for example all of `test_growth_integration` (0.00–0.01 s) and most of the `nested_built` tests (0.00 s).
- **Unmarked but among the slowest:**
  - The three mutate/repartition tests (1.3–1.5 s).
  - `test_bridge_identity_on_the_inversion_saddle` (1.35 s).
  - `test_faces_closed_raises_at_the_cap_for_the_anchor` (1.38 s).
  - `test_the_fixture_really_is_an_inversion_point` (1.23 s; it absorbs the module fixture build).
  - The `henon_inversion` users: 0 of 11 marked.
- **Same fixture, mixed marks:**
  - `k28_partitioned`: 3 of 6 users marked. Unmarked: `test_bridge_class::test_k28_has_one_active_and_two_inert_classes`, `test_session_bridge_classes::test_k28_letters_only_the_anchor_class`, `test_stable_partition::test_k28_blast_child_gets_no_forward_hole`.
  - `henon_p3_session`: 52 of 83 marked. `p3_partitioned`: 18 of 31 marked.
- **Stale docstrings:**
  - The conftest module docstring says `henon_p3_session` is "used only by the blast regression test"; it is used by 83 tests in 24 files.
  - The `k28_partitioned` docstring says "about a minute to build"; it takes 0.34 s.

**`regression`** is on 9 tests (the 5 in `tests/regression/`, 3 in test_higher_period_cartoon, 1 in test_stable_partition_period3). About 22 test files describe fixed bugs ("used to", "bug", "previously") without it, including `numerics/test_workbench_bugfixes.py` (10 tests) and `test_intersection_registry_fixes.py` (13 tests), both with 0 markers. The `tests/regression/` directory and the marker are two overlapping, both partial, ways of saying the same thing.

**Recommendation:** either drop both markers, or redefine `slow` as "> 0.5 s" (about 10 tests) and apply it mechanically. The suite at about 45 s total does not need a fast lane.

## 5. tests/visualizations/

- **Not collected.** `viz_nested_bridges.py` and `viz_unstable_artifacts.py` don't match `test_*.py` and contain no test functions. They are scripts that write PNGs to `/tmp` and print.
- **Dead or broken:**
  - The deferred list in `docs/regions_bridge_identity_and_cleanup_plan.md:333` records that both crash with `KeyError: None` (a None strong pip on their shorter growth). It was deferred to "Phase 8" and never fixed; the last commit touching them is Phase 8 cleanup `66f5688`.
  - `viz_unstable_artifacts.py` captures log records containing `"Discarding impossible"`. That string no longer exists: `Tangle._discard_same_stability` now logs `"Dropping same-stability (...) segment pair ... as a near-tangency"` at DEBUG. So even if the script ran, it would capture nothing. It also reads the private `tangle._seg_lookup`.
  - `viz_nested_bridges.py` investigates a "nested cdist span" bridge-cutting bug through "bridge-signature helpers" that no longer exist by that name.
  - Both carry a third, private copy of the nested-p3 recipe.
- **Stray bytecode:** `__pycache__` holds `-pytest-8.3.4.pyc` files, so pytest imported them at some point. Stale bytecode only.
- **Recommendation:** delete them, or move them to `scripts/` (archive) if the investigations still matter. They do not belong under `tests/`.

## 6. Other config notes

- `pythonpath = ["src", "tests"]` is what lets test files do `from invariants import ...`, `from walk_helpers import ...` and `from minimal_helpers import ...`. Fine, but making `tests/` a package with relative imports would be clearer.
- The `henon_map`, `henon_map_inverse` and `henon_jacobian` fixtures just return module globals. Plain imports from `tanglepack.examples` would do, and several files already bypass them.

## 7. Suggested direction for the refactor plan (infrastructure only)

1. **One case factory.** Add `build_k10(unstable_steps=9, partition=False, infer=True)` and `build_k28(blasts=n)` to `henon_cases` next to `build_period3` / `build_nested`. Every k10/k28/p3 fixture and local fixture becomes a one-line call. Delete the 7 duplicate k10 recipes, the shadowing `k10_session`s, the duplicated `_define_zone`, and the second nested recipe (or keep `henon_p3_session` explicitly as "nested, p1 cutoff 1e-4" and rename `inner_zone` to `outer_zone`).
2. **Scope.** Make read-only fixtures session- or module-scoped: k10_partitioned, p3_partitioned, the k28 cases and the inversion case. Give mutating tests an explicit fresh-build factory fixture. This is an estimated 15–20 s saving on a 45 s suite.
3. **Cache-invalidation tests.** Fold the three mutate/repartition tests into one parametrized test over (bridge_classes, minimal_trellis, iterated_partition, dual_graph, symbolic_dynamics) on a single mutated session.
4. **Fakes.** Have `FakeFixedPoint` delegate to the real `FixedPoint.advance_key` / `per_step_beta` logic, for example via a tiny real `FixedPoint` or a shared pure function, instead of re-implementing it. Have `FakeDual` construct nodes through a public `StableNode` factory if one exists.
5. **`build_pieces`.** Replace it with the session API, or keep it only where a test deliberately checks "direct assembly == session". Fix the un-lettered-table workaround at its source.
6. **Invariants.** Change `assert_cdist_monotonic` to default `strict=False`. Drop or quarantine `assert_no_cdist_collision`. Apply `assert_area_preserved_along_chain` to the blasted and period-3 fixtures as a real invariant test. Fix the spike error message.
7. **Markers.** Either remove `slow` and `regression`, or define `slow` by measured cost (> 0.5 s) and keep `regression` only as the directory. Fix the stale conftest and k28 docstrings.
8. **Visualizations.** Delete `tests/visualizations/`, or move it to `scripts/archive/`.

Files referenced:
- /home/dezhu/Development/TanglePack/tangle-pack/tests/conftest.py
- /home/dezhu/Development/TanglePack/tangle-pack/tests/invariants.py
- /home/dezhu/Development/TanglePack/tangle-pack/tests/walk_helpers.py
- /home/dezhu/Development/TanglePack/tangle-pack/tests/minimal_helpers.py
- /home/dezhu/Development/TanglePack/tangle-pack/tests/visualizations/viz_nested_bridges.py
- /home/dezhu/Development/TanglePack/tangle-pack/tests/visualizations/viz_unstable_artifacts.py
- /home/dezhu/Development/TanglePack/tangle-pack/src/tanglepack/examples/henon_cases.py
- /home/dezhu/Development/TanglePack/tangle-pack/pyproject.toml