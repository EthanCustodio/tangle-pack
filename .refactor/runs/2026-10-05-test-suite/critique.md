## Corrections

I checked 18 of the synthesis's claims against the code. Most hold: test_get_point is tautological, CuPy is installed so the GPU no-CuPy test always skips, test_k10_region_images is vacuous, the oriented_bridge_polyline test never checks orientation, the test `_empty_stretches` lacks the 09-30 base-branch chord rule, `punch_holes` already runs `check_holes_share_bridge_side` (Trellis.py:1117), on_interval and `_get_lambda_u` have no callers in src, `henon_p3_session.inner_zone` is the max-area outer zone, the refinement test's docstring calls the mean rule suspect, the inversion test uses the word "defect", the five test_fixed_points_* tests exist, and the two "no class mixes fixed points" tests are verbatim twins. These claims are wrong or overstated:

1. **test_direct_hole_side_is_the_side_of_its_coordinates (×4) is not a "superseded oracle" (theme C).** Its docstring says the crossing-sign rule of 2026-10-02 must reproduce the old geometric measurement hole for hole. Two independent methods agreeing is the strongest kind of test in the suite. Keep it, and consider adding nested to its parametrization.
2. **test_henon_propagated_hole_side_matches_coords** is also an independent geometric cross-check of I1 on k=10, not an obsolete rule. At most, extend it to p3 using the image sub-arc reading. Do not delete it.
3. **test_direction_matches_the_element_order_at_the_ends is not superseded.** `anchor_outward_key`'s own Note says that on one branch and side it equals `element_sort_key`, and k=10 has one branch. The test is correct, only narrow. Switching it to `anchor_outward_key` is a one-line change. The real gap is cross-branch orientation (already §5.3).
4. **test_low_stretch_growth_keeps_cdist_injective is not superseded.** `tests/invariants.py` (lines 12–17) says ties are legitimate only at high-stretch folds and to "pass strict=True only for low-stretch growth known to stay injective". This test is exactly that sanctioned case. The "files disagree on strictness" finding is mostly by design. The real problem is only that the helper defaults to `strict=True`. Option (c) of Q10 would remove a check the author added on purpose.
5. **test_k28_two_blasts_structure: no docstring says "never pin a letter".** A grep for never/not/don't pin finds nothing relevant. The section header says "pinned by element names", and the test pins names and then letters too. Q3 is still a fair question, but the framing that the test contradicts its own documentation is unsupported.
6. **numerics/test_intersection_registry_fixes "85% flagged" is too high.** The test_fixed_points_* tests (×5) guard the real `Intersection.fixed_points` IndexError when manifold_a_key is None after a blast. The test_synthetic_* tests (×3) guard the bug where a label landed in the id slot. Both bugs are listed in memory/codebase-audit-2026-09. Only the four `_get_lambda_u` tests and the one on_interval test are candidates.
7. **The "unlinked halves of the table-linked tests (collision_rtol < 2%)" in theme K guard the 2026-09-15 pseudoneighbor collision bug.** In that bug, pair (10,9) was dropped at 8 blasts (memory/pseudoneighbor-collision-fragility). The unlinked half shows that the 2-D collision path alone fails, so its tolerance dependence is the point of the test. Restate it against the `collision_rtol` constant. Do not loosen it by 10× (Q26b would break it).
8. **test_bridge_endpoints_stay_on_the_bridges_own_branch (cluster 3) guards the root cause of the period>1 partition bug** (`_assign_bridge_intersections` with no branch filter, all anchors at cdist 0). Any merge must keep it running on a period-3 fixture. A k=10 version cannot fail this way.
9. **test_kevin_way_period_three_orbit_chain guards the known Kevin-way ordering bug.** That bug's (orbit, branch) order disagreed with the old `_advance_key_forward` (codebase-audit memory). It should be strengthened to assert the chain order, not treated as weak filler.
10. **Two tombstones in cluster 25 enforce "one rule" design invariants stated in CLAUDE.md, not dead history.**
    - test_collect_is_the_only_traversal enforces "one `_collect` walker behind every getter". It patches `_collect` behaviourally rather than grepping source, and it relates to the `_segments_of` tail bug.
    - test_workbench_has_no_private_key_advance enforces "`advance_key` is the ONE rule". A second key-advance rule caused the Kevin-way bug.

    Deleting them under Q21(a) with no flag is wrong. List them under Q21(b).
11. **test_session_strong_pips::test_trellis_auto_rebuilds is a near duplicate, not an exact one.** test_compute_intersections_invalidates_the_cache does the same thing on k=10 with one fixed point. The strong-pips version runs on the nested two-fixed-point session. The duplication is real but costs little.
12. **The Q39 strawman's mutation guard is insufficient.** Checking `workbench.generation` at teardown misses most of the shared-state mutations tests actually perform. `trellis.classify_strong_pips()`, `compute_pseudoneighbors()`, `punch_holes()`, `partition_stable_manifold()`, `set_strong_pip`, and the `session._bridge_classes` / alphabet state all mutate cached Trellis and session objects without bumping the generation. A frozen session-scoped fixture needs a deeper fingerprint (holes, partition signature, pips, alphabet), or the topology steps must run inside the fixture and the tests must be forbidden from calling them.

## Missed overtesting

- **Cache-identity assertions inside fact tests.** For example, test_k28_two_blasts_structure line 2 asserts `dyn is session.symbolic_dynamics([fp])`. These are hidden cache-policy pins, and the cache tier will not catch them when consolidated.
- **Evidence treated as contract.**
  - Four tests assert `verified is True`.
  - Ten assert `is_reliable`.
  - Two assert `cd.source == "walk"`.

  CLAUDE.md says registered evidence is "evidence only, the walk is the answer", and a mismatch is INFO or WARNING. Pinning `verified` turns a diagnostic into a requirement. `source == "walk"` pins which resolution path ran, not the result.
- **Exact float pin:** `table.entries[0].min_unstable_cdist == 0.0`.
- **Seven `print(` calls in test files** (for example test_k10_region_images_agree).
- **The 1.6 s test_session_plot_delegates_draw / fan-out tests use monkeypatch.** Ten test files use monkeypatch. The synthesis counts the plotting ones but not the monkeypatched `_inert_letter`, `_image_chain` and `_empty_stretches` in algorithm tests. Those make the tests depend on the internal call graph.
- **`henon_p3_session` is function-scoped because blast tests mutate it, but non-mutating tests also use it.** The scope split is a cost issue, not only Q39's.
- **Old parallel snapshot builds.** test_new_grow_manifold_matches_old (60 lines of frozen legacy code) is cited, but `minimal_helpers.build_pieces` and the copied dual-graph assembly helpers run a second, slightly different pipeline. Tests then pin facts on both the session pipeline and the hand-assembled one.

## Missed questions

1. **Known-open failures:** p3 at 15–16 steps or 6+ blasts has an unreachable class and a virtual `new1`. Should these be `xfail(strict=True)` tests that document the open issue, or left untested?
2. **Production-side invariant checks:** `punch_holes` runs I1 and I2 internally. Should the invariant tier rely on production asserts (and test that they are wired), or should invariants live only in tests? This decides whether test_henon_holes_share_bridge_side and test_henon_partition_still_built_after_wiring are redundant.
3. **Fixtures that some invariants do not cover:** the inversion path was historically "unvalidated everywhere". Which invariants are expected to hold on the inversion fixture (partition, classes, symbolic dynamics)? The ALL_CASES strawman assumes all of them.
4. **Golden re-recording authority:** for a scientific code, may an `--update-golden` flag re-record values, or must every golden change be reviewed by the author? Automatic re-recording defeats the purpose.
5. **Running the GPU tests:** CuPy is installed on this machine. Is GPU parity a supported contract, and should test_gpu_growth_matches_cpu run in the default suite?
6. **Smoke tests for the figure scripts:** `scripts/henon_symbolic_itineraries*.py` and `symbolic_figures.draw_panels` are the main consumers. Should one cheap end-to-end figure run be a test?
7. **The "one rule" design tests:** should "single-source-of-truth" structural rules (one walker, one advance_key, one Hénon definition) be enforced as tests, as lint, or not at all? This is separate from deleting dead-attribute tombstones.
8. **CI:** does CI exist, and what wall time is acceptable? Most marker and perf decisions depend on this.
9. **Success measure:** the user's concern is the quantity of tests. Is 600–630 cases (only 25% fewer) acceptable, or is the target a much smaller suite? The invariant tier multiplies cases.
10. **Is `verified`/`is_reliable` meant to be guaranteed on the fixtures?** See "Missed overtesting".

## Verdict on the strawman architecture

The tiering is sound: invariants, then unit, then facade, then a once-only golden tier. It matches CLAUDE.md, and so do the fixed-point rule, the frozen-builder idea and the cache table. Before planning on it:

- **Q39's frozen fixtures are unsafe as specified.** See correction 12. They need a state fingerprint beyond the generation, or a hard rule that topology steps run only inside fixtures.
- **ALL_CASES needs per-case applicability.** The inversion fixture, the one-blast k28 with 0 image bridges, and nested will each need skips or xfails. Otherwise checks pass vacuously, which is the same "vacuous on k28" weakness the synthesis flags elsewhere. Every parametrized invariant needs a non-vacuity guard (count > 0).
- **"Unit tier, no Hénon build where avoidable" pushes toward more fakes.** That conflicts with the Q42 concern about fakes that reimplement the one-rule functions. Prefer small real builds or fakes that delegate to the real rules.
- **The suite shrinks less than the user may expect.** The function count roughly halves, but collected cases drop only about 25%. Lead with that. The 25–30 s estimate depends on Q39 holding, and building the nested and two-blast fixtures once per session is the floor.
- **Moving every file into the mirrored layout costs git-blame continuity.** Consider deleting and consolidating first, then moving files once.
- **Before anything is deleted, add a "Regression guard" column** naming the bug each candidate protects (corrections 6–10), using memory/codebase-audit-2026-09.md and pseudoneighbor-collision-fragility.md.