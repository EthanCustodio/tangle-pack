# Test evidence — B8 `unresolved-sets` (test engineer, traversal 1, 2026-10-13)

## Where

- **Base:** `8a1cdc3`, detached worktree `/home/dezhu/Development/TanglePack/tangle-pack-unresolved-sets-base`. It is clean: test files were copied in for each run and then restored with `git checkout -- tests`.
- **Branch:** `refactor/unresolved-sets`, worktree `/home/dezhu/Development/TanglePack/tangle-pack-unresolved-sets`. The coder's head was `2b21439`. My test commit on top is **`ae18859`**.
- **Author's tree:** untouched, apart from this artifact.
- **Scripts:** every probe and break script is in `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/te/`.
  - `pin.py` writes the dumps; `diff.py` compares them.
  - `zone.py` is rubric 11 on base and branch.
  - `law.py` is the rubric-14 probe; `lawbreaks.py` breaks the law's clauses.
  - `b8breaks.py` is every other break.
  - `gen.py` explains the counter.
  - `invwarn.py` checks for a spurious warning.
  - Outputs: `*.json`, `*.log`, `breaks*.out`, `lawbreaks*.out`.

## Verdict

- **Behaviour is preserved where the plan says it is.**
  - Fingerprints match on every law case. The one exception is the observation counter `workbench generation`, explained below.
  - The `+1` link sets are identical to base on every law case and on every step of every series.
  - `registry.unresolved_links` is empty on the branch everywhere.
- **Every guarded path has a test with teeth**, with two exceptions, both explained under the break-it checks:
  - one mutant is behaviour-equivalent;
  - the law's new clauses are guarded by the rubric-14 probe only, as the plan prescribes.
- **I added 3 tests and strengthened 1.** Without them, three breaks went undetected: dropping the `others`/`earlier` term, the per-pass pop, and the generation bump on recording alone.

## 1. Full suite

| tree | command | result |
|---|---|---|
| branch `ae18859` | `MPLBACKEND=Agg env/bin/python -m pytest -q -rxX` | **1171 passed, 38 skipped, 5 deselected, 1 xfailed** in 197.6 s |
| base `8a1cdc3` | same command | 1162 passed, 38 skipped, 5 deselected, 1 xfailed in 196.4 s |

- On both trees the single XFAIL is `test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing]` (C3), the only `cases.KNOWN_ISSUES` entry. **Rubric 1: PASS.**
- The difference of 9 tests is exactly the coder's 6 new tests plus my 3. No existing test was weakened or deleted.

## 2. New tests: pass on the branch, fail on base (rubric 2)

**Command (base):**
```
cp ../tangle-pack-unresolved-sets/tests/unit/numerics/test_iterate_inference.py tests/unit/numerics/
cp ../tangle-pack-unresolved-sets/tests/facade/test_session_topology.py tests/facade/
MPLBACKEND=Agg ../tangle-pack/env/bin/python -m pytest -q -p no:cacheprovider tests/unit/numerics/test_iterate_inference.py tests/facade/test_session_topology.py
git checkout -- tests
```

**Results:**
- **Base:** `9 failed, 11 passed`. The 11 that pass are the existing tests of the two modules.
- **Branch:** `20 passed`. That is 8 in the unit file and 12 in the facade file, counted from the break runs; the full suite is green too.

| test (file) | guards | base failure |
|---|---|---|
| `test_an_image_claimed_across_passes_links_neither_claimant` (unit, coder; rubric 8) | a held link `z->y` meets a new match `x->y`: both unlinked, both `(y,)`, one WARNING naming `y`, the graph edge gone | `assert (None is None and 3 is None)` — `z` keeps `y` |
| `test_an_image_claimed_three_times_in_one_pass_links_none` (unit, coder; rubric 9) | grouping by image: three claimants in one pass, none linked, all recorded | `assert 2 is None` — the first claimant is linked |
| `test_a_multi_match_is_recorded_and_a_second_pass_changes_nothing` (unit, coder; **strengthened by me**; rubric 10) | multi-match recorded with `>1` candidates and unlinked; the first pass bumps the generation, the second does not | `AttributeError: ... unresolved_links` |
| `test_remove_iterate_forgets_the_link_its_twin_and_its_edge` (unit, **new**; rubric 7) | `remove_iterate(x,1)` on a real inferred link: no `(x,1)`, no `(y,-1)`, no edge `iter_x_1`, generation bumped | `AttributeError: ... remove_iterate` |
| `test_a_later_pass_does_not_relink_a_contested_image` (unit, **new**; plan step 3, "the others term") | after the cross-pass collision, a second partial pass `infer_iterates(wb, [x])` leaves `x` unlinked, the set unchanged and the generation unchanged | `assert (None is None and 2 is None)` — on base `z` holds `y` |
| `test_a_pass_re_decides_its_own_sources` (unit, **new**; plan step 3, "drops the entries of its own sources") | entries recorded under wide brackets vanish on the next pass with real brackets, and the original links come back | `AttributeError: ... unresolved_links` |
| `test_a_zone_at_a_non_candidate_keeps_the_chosen_pip` (facade, coder; rubric 11) | `resonance_zone` at a non-candidate warns once, keeps the chosen pip, `topology()` runs | `0 WARNING record(s) ... expected 1` |
| `test_a_forgotten_pip_gives_the_default_back` (facade, coder; rubric 12) | `forget_strong_pip` uncuts, the default comes back on trellis and session | `AttributeError: ... forget_strong_pip` |
| `test_a_stale_pip_error_names_the_remedy` (facade, coder; rubric 13) | the `ValueError` matches "no longer a strong-pip candidate" and names `forget_strong_pip` | `assert 'forget_strong_pip' in '... choose another with set_strong_pip()'` |

**My change to the coder's multi-match test.** I strengthened it; nothing was weakened.
- It now passes `detected` as the sources of both passes. With `sources=None`, the first pass also re-links the anchors, and the generation moves anyway.
- I added `assert generation > start`, so recording alone must be visible to a cache.
- Its other assertions are unchanged.

**Rubric 11 on base, checked directly.** `zone.py` builds `build_inversion(blasts=2)` and picks the outermost non-candidate off the chosen pip's stable branch: crossing **33** on `(stable,0,1)`; the chosen pip is 1 on `(stable,0,0)`.

| tree | chosen pip unchanged after `resonance_zone(33)` | `topology()` | WARNING |
|---|---|---|---|
| base | no (adopted 33) | **raises `ValueError`** "... no longer a strong-pip candidate; choose another with set_strong_pip()" | one, from `_restore_strong_pip` |
| branch | **yes** | **runs** | one: "the zone pip 33 ... is not a strong-pip candidate; the zone is stored, the chosen strong pip is unchanged" |

## 3. Pinning the old behaviour (rubrics 4 and 5)

**Method.** `pin.py <worktree> <out.json> laws|series` runs in four separate processes, each with `-I` and its worktree's `src` and `tests` first on the path.
- **`laws`:** for each of the 11 `cases.LAW_CASES`, it dumps `law_tier.fingerprint(case)`, canonicalised: fixed points by label and sets sorted. `products` is dropped because it is a tuple of `id()`s. It also dumps the **full sorted `+1` link set** (not only its count) and `unresolved_links`.
- **`series`:** it wraps `cases._recomputed` and records the `+1` link set and `unresolved_links` at every step of every `cases.SERIES`, both **after the blast** and **after the step's plain recompute**. It also dumps the canonicalised `SeriesStep` tuple.
- `diff.py base.X.json branch.X.json` compares the dumps aspect by aspect.

**Law cases (`diff.py base.laws.json branch.laws.json`).**

| case | `+1` link set equal (count) | branch `unresolved_links` | aspects differing |
|---|---|---|---|
| k10 | yes (4) | {} | none |
| k28_one_blast | yes (4) | {} | workbench generation 2->3 |
| k28_two_blasts | yes (6) | {} | workbench generation 2->3 |
| p3 | yes (28) | {} | none |
| nested | yes (34) | {} | workbench generation 2->4 |
| nested_inner_blasts | yes (38) | {} | workbench generation 2->4 |
| inversion | yes (9) | {} | none |
| inversion_one_blast | yes (13) | {} | workbench generation 2->3 |
| inversion_two_blasts | yes (19) | {} | workbench generation 2->3 |
| inversion_period2 | yes (16) | {} | none |
| heteroclinic | yes (94) | {} | workbench generation 2->4 |

On every case these aspects are equal: `registry generation`, `crossings`, `bridge ids`, `holes`, `pips`, `partition signature` and `alphabet`.

**The `workbench generation` difference is a read artifact, not a state change.**
- `TangleWorkbench.generation` counts *observed* changes ("compare for equality only; the value counts observed changes, not mutations"). Step 4's candidacy test reads it once more per zone, by building the trimmed trellis.
- `gen.py` proves this on **base**. It wraps `resonance_zone` with nothing but an extra `self.trellis(rz.key)` read, and base moves the same way:

  | case (base) | as is | with one extra read |
  |---|---|---|
  | k28_one_blast | 2 | 3 |
  | nested | 2 | 4 |

- **Rubric 4:** PASS on every state-carrying aspect. On the raw counter it is a strict FAIL, and the cause is fully accounted for above. It is the reviewer's call.

**Series (`diff.py base.series.json branch.series.json`).** Branch `+1` counts per step are shown; the link sets equal base at every step, both after the blast and after the recompute.

| series | steps | `+1` link sets equal at every step, both points | branch `unresolved_links` empty at every step, both points | `SeriesStep`s equal | branch `+1` counts |
|---|---|---|---|---|---|
| k28_series | 9 | yes | yes | yes | 2,4,6,10,18,30,50,84,140 |
| p3_series | 10 | yes | yes | yes | 28,30,32,34,36,38,40,46,52,57 |
| inversion_series | 5 | yes | yes | yes | 9,13,19,25,35 |
| k32_series | 11 | yes | yes | yes | 2,4,6,12,22,36,58,100,174,299,517 |
| nested_outer_series | 12 | yes | yes | yes | 30,32,34,38,46,58,78,114,174,273,443,729 |

**Rubric 5: PASS.**

## 4. Break-it checks (branch; each break is restored with `git checkout -- <file>`)

Run with `b8breaks.py`, which applies one string edit, runs the guarding test file(s) and restores. The tree was clean after every run. Rows marked † were confirmed by a second run after my last strengthening of the multi-match test (`breaks2.out`).

| # | break | file | tests run | result |
|---|---|---|---|---|
| 1 | **grouping removed**: walk the matches in order and check only the held link (base's shape) | IterateInference | unit | **FAIL** — `three_times_in_one_pass`, `later_pass_does_not_relink` |
| 2 | grouping removed, but the in-order loop keeps the live `earlier` term | IterateInference | unit | passes — **an equivalent mutant**: c1 links, c2 unlinks c1 and records both, c3 then sees `earlier=[c1,c2]` and is not linked. The outcome is identical by construction, so no test can tell it apart. The plan's "would link the third" holds for break 1. |
| 3 | **`others`/`earlier` term dropped** (`earlier = []`) | IterateInference | unit | **FAIL** — `later_pass_does_not_relink` (new; the coder's tests all passed under this break) |
| 4 | collision does not call `remove_iterate(held, 1)` | IterateInference | unit | **FAIL** — `across_passes`, `later_pass` |
| 5 | **`IterateTable.remove` leaves the backward twin** | IterateTable | unit | **FAIL** — `across_passes`, `remove_iterate_forgets`, `later_pass` |
| 6 | **`remove_iterate` leaves the graph edge** | IntersectionRegistry | unit | **FAIL** — `across_passes`, `remove_iterate_forgets` |
| 7 | `remove_iterate` does not bump the generation | IntersectionRegistry | unit | **FAIL** — `remove_iterate_forgets` |
| 8† | multi-match not recorded | IterateInference | unit | **FAIL** — `multi_match`, `re_decides` |
| 9 | a pass does not pop its own sources' entries | IterateInference | unit | **FAIL** — `re_decides` (new; nothing caught it before) |
| 10† | generation bumped unconditionally | IterateInference | unit | **FAIL** — `multi_match`, `later_pass` |
| 11† | generation never bumped for `unresolved_links` | IterateInference | unit | **FAIL** — `multi_match` (only after my strengthening; with the coder's version it passed) |
| 12 | **`forget_strong_pip` without `trellis(fp).strong_pip = None`** | TangleSession | facade | **FAIL** — `forgotten_pip_gives_the_default_back` |
| 13 | **candidacy test in `resonance_zone` removed** (`if found is not None:`) | TangleSession | facade | **FAIL** — `zone_at_a_non_candidate` |
| 14 | `resonance_zone` back to base (adopt and return before any test) | TangleSession | facade | **FAIL** — `zone_at_a_non_candidate` |
| 15 | non-candidate path does not put the previous choice back | TangleSession | facade | **FAIL** — `zone_at_a_non_candidate` |
| 16 | the plan's literal step 4 (`get` instead of `pop` of the previous choice, the coder's deviation 1) | TangleSession | facade + `tests/invariants` + `tests/regression` | **FAIL** — `zone_at_a_non_candidate` (1 failed, 257 passed); see the note below |
| 17 | `Topology._choose_pip` message without `forget_strong_pip` | Topology | facade | **FAIL** — `stale_pip_error_names_the_remedy` |

**Note on break 16.** `invwarn.py` counts `TangleSession` WARNINGs while building `build_inversion(blasts=2)`:

| tree | WARNINGs |
|---|---|
| base | 0 |
| branch | 0 |
| branch with break 16 | **1** (the spurious "chosen strong pip ... is no longer a strong-pip candidate") |

This confirms the coder's deviation 1 is needed. No committed test asserts a silent inversion build, so only the zone test's warning count guards it.

**The law's new clauses (rubric 14; `lawbreaks.py`, probe `law.py` on k10 at the bridges stage).**

Intact branch, each hand-built registry against `check_iterate_table_complete`:

| hand-built registry | result |
|---|---|
| clean | passes, 7 checked |
| `{z:(y,)}` while `y` is linked from `x` | **rejected**: "image 3 is contested by 4 but linked to 1" |
| `{x:(y,)}` with `x` linked | **rejected**: "crossing 1 is unresolved but linked" |
| `{z:()}` | **rejected**: "crossing 4 is unresolved with no candidate image" |
| `{x:(y,z)}` with `x` linked | **rejected**: "crossing 1 is unresolved but linked" |
| `x` unlinked and not recorded | **rejected** (the old completeness clause): "crossing 1 is unlinked, but [3] is its image" |
| `x` unlinked and recorded as `(y,z)` | **accepted** (the new skip) |

Each clause broken in turn:

| clause removed | what the probe shows | committed law callers (`tests/invariants -k iterate` + `test_iterates_after_blasts.py`) |
|---|---|---|
| no-candidate | `{z:()}` -> **NOT rejected** | 13 passed, 4 skipped |
| unresolved-but-linked | `{x:(y,z)}` -> **NOT rejected**. `{x:(y,)}` is still caught by the contested clause, so the two overlap for single-candidate entries. | 13 passed |
| contested image | `{z:(y,)}` -> **NOT rejected** | 13 passed |
| skip of unresolved sources | an unlinked `x` recorded as `(y,z)` -> **rejected** ("unlinked, but [3] is its image") | 13 passed |

The committed callers pass under every law break, as expected: the set is empty on every law case and series step (section 3). So, as the plan prescribes, the law's clauses are guarded by the probe only, not by a committed test.

## 5. Static checks

- **Rubric 3:** `git diff --stat 8a1cdc3` on the branch gives `8 files changed, 301 insertions(+), 31 deletions(-)`, exactly the plan's eight files. My commit touches only `tests/unit/numerics/test_iterate_inference.py`. **PASS.**
- **Rubric 6:** **PASS.**
  - `grep -rn "already has the preimage" src` finds nothing (exit 1).
  - `grep -n "none linked" src/tanglepack/numerics/IterateInference.py` finds line 150 only, the WARNING. Its `else` branch calls `registry.remove_iterate(held, 1)` and `registry.unresolved_links.update(...)`.
  - There is no `continue` in lines 140-160.
- **Rubric 7:** **PASS.** Both methods exist. The new direct test and break checks 5-7 show that `(x,1)`, `(y,-1)` and the edge `iter_x_1` are all gone.

## 6. Not tested, and why

- **A real collision on a built case.** None fires on any law case, any series step or k28 x8, per the analyst's probe and section 3. The collision paths are therefore tested only on fixtures arranged by hand (a hand-linked table, a monkeypatched `_bracketed_image`/`_map_gated`).
- **`near_tangency`'s ~25 multi-match entries.** These are behaviour outside the law cases that nothing reads. I did not pin them. A count is not a law, and the test policy forbids case-specific pins.
- **`add_resonance_zones` with a non-candidate pip** (plan question 1). There, `recompute=False` runs the test before the batch's single recompute. No recipe has a non-candidate zone pip, and the plan's default is to leave the batch unchanged.
- **`IntersectionRegistry.reindex_from` resetting `unresolved_links`.** This is a no-op in practice, because it runs before inference. No test reaches a state where it matters.
