# Change — B8 `unresolved-sets` (coder, iteration 1, 2026-10-13)

## Where

- **Base commit:** `8a1cdc3c697ae0e631984cf9085386009efd0593` (the author's `HEAD`, `main`).
- **Branch:** `refactor/unresolved-sets`. Worktree: `/home/dezhu/Development/TanglePack/tangle-pack-unresolved-sets`. Head is `2b21439`.
- **Base worktree:** `/home/dezhu/Development/TanglePack/tangle-pack-unresolved-sets-base`, detached at `8a1cdc3`. I copied files into its `tests/` for the fail-on-base runs and then restored them with `git checkout -- tests`, so it is clean.
- **Author's tree:** untouched. The uncommitted `docs/` and `.refactor/` edits, the staged rename and the untracked files are as I found them. The only new file there is this artifact.

## Commits (oldest first)

1. `6848552` Registry: `unresolved_links` holds the links the brackets cannot decide; `IterateTable.remove` / `remove_iterate` forget a link (step 1).
2. `1a6cd3b` IterateInference: a multi-match, and every claimant of a contested image, go to `unresolved_links`; none is linked, and no order decides (steps 2-3).
3. `960674b` Session: `forget_strong_pip`; `resonance_zone` adopts its pip only if it is a candidate on the trimmed trellis; the stale-pip error names the remedy (step 4).
4. `c65c712` Laws: the iterate table is complete except for `unresolved_links`. Each entry holds a candidate, is unlinked, and leaves a contested image without a preimage (step 5).
5. `99c55aa` Session: `resonance_zone` tests its pip on a trellis built without the previous choice, so a trim that cuts the old pip away logs nothing (step 4, a fix to commit 3; see "Deviations").
6. `a001961` Tests: rubric 8-10.
7. `2b21439` Tests: rubric 11-13.

## `git diff --stat 8a1cdc3`

```
 src/tanglepack/loom/TangleSession.py            | 35 +++++++++++--
 src/tanglepack/loom/Topology.py                 |  6 ++-
 src/tanglepack/numerics/IntersectionRegistry.py | 20 ++++++++
 src/tanglepack/numerics/IterateInference.py     | 65 +++++++++++++++---------
 src/tanglepack/numerics/IterateTable.py         |  9 ++++
 tests/facade/test_session_topology.py           | 61 ++++++++++++++++++++++-
 tests/helpers/laws.py                           | 16 +++++-
 tests/unit/numerics/test_iterate_inference.py   | 66 +++++++++++++++++++++++++
 8 files changed, 247 insertions(+), 31 deletions(-)
```

These are exactly the plan's eight files. Source: +104/-31 lines (most of it docstrings). Tests: about 143 lines.

## Suite

- **Branch:** `MPLBACKEND=Agg env/bin/python -m pytest -q -rxX` gives `1168 passed, 38 skipped, 5 deselected, 1 xfailed in 202.63s`. The one XFAIL is `test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing]` (C3), which is the only entry in `cases.KNOWN_ISSUES`.
- **Base `8a1cdc3`, same command:** `1162 passed, 38 skipped, 5 deselected, 1 xfailed in 200.51s`. The difference is exactly the 6 new tests.

## What changed

- **`IterateTable.remove(source, n)`** pops `_forward[source][n]` and `_backward[target][n]`. It raises `KeyError` when there is no such link.
- **`IntersectionRegistry`**
  - New attribute `unresolved_links: dict[int, tuple[int, ...]]`, documented in Attributes. It starts empty in `__init__` and is also reset in `reindex_from`, next to the `IterateTable()` reset there. Ids are renumbered there and the table is emptied too; `reindex_from` runs before inference, so this is a no-op in practice.
  - New `remove_iterate(source, n)` removes the table link, bumps `_generation`, and drops the graph edge `iter_{source}_{n}` when present.
- **`infer_iterates`**
  - At the start of each pass it snapshots `before = dict(registry.unresolved_links)` and pops its own sources' entries, so each pass re-decides its own sources.
  - The early return on "no detected sources" is gone. `matches` is simply `{}` in that case, so the generation compare at the end covers every path.
  - After the map gate it groups matches by image. For each image `j`, the contested set is: this pass's claimants, plus `table[j, -1]` if one is held, plus every earlier entry equal to `(j,)`.
  - If that set has exactly one member, the link is registered. Otherwise there is one WARNING ("image %d is claimed by crossings %s; none linked"), the held link is removed with `remove_iterate`, and every claimant is recorded as `(j,)`.
  - At the end, `registry._generation += 1` runs only if `unresolved_links != before`.
  - The docstring's description and Note say this.
- **`_bracketed_image`:** a multi-match records `unresolved_links[src] = tuple(images)` and keeps its INFO line.
- **`TangleSession.forget_strong_pip(fp)`** pops `_strong_pips[fp]` and sets `self.trellis(fp).strong_pip = None`, which uncuts the trellis.
- **`TangleSession.resonance_zone`:** after the zone is stored, it pops the previous choice, builds `self.trellis(rz.key)`, and runs `registry.find(pip)` and `_is_candidate`, the same test `_restore_strong_pip` uses.
  - Candidate: it sets `trellis.strong_pip = found` and remembers the pip.
  - Otherwise: one WARNING. The zone stays stored, and the previous choice is put back and restored on that trellis.
  - The docstring says this.
- **Error texts:**
  - `_restore_strong_pip`'s WARNING and `Topology._choose_pip`'s `ValueError` keep "no longer a strong-pip candidate" and now end "choose another with set_strong_pip() or drop it with forget_strong_pip()".
  - The `Raises:` sections of `TangleSession.topology` and `Topology.__init__` name `forget_strong_pip`.
- **Law (`check_iterate_table_complete`):**
  - Sources in `unresolved_links` are skipped, the same way linked ones are.
  - Every entry must hold at least one candidate and must have no `+1` link.
  - An entry `(j,)` must not coexist with a `(j, -1)` link.

## Deviations from the plan text (all inside the eight files)

1. **`resonance_zone` pops the previous choice before building the trellis (commit 5).**
   - **Why:** taken literally, step 4's `t = self.trellis(rz.key)` builds the trellis with the OLD remembered pip restored. In recipes whose zone trims that pip away, `_restore_strong_pip` then logs a spurious "no longer a strong-pip candidate" WARNING. The inversion build is one: the default pip on `(0,1)` at cdist 13.13, zone at `f(q0)`. Base never logs that, because it never builds a trellis inside `resonance_zone`.
   - **Fix:** popping first means only the zone's pip is tested. In the non-candidate case the previous choice is put back and restored, which warns only if it really is stale.
   - **Also needed:** on the candidate path, `trellis.strong_pip = found` is required. Otherwise `Topology._choose_pip`'s "one set on the trellis is adopted" would re-adopt whatever pip the cached trellis holds.
2. **Step 3's variable names** are `claimed_by`, `claimants`, `held`, `earlier` and `contested`, and the single-claimant branch links `contested[0]`. The plan's `i := sources[0]` would shadow the `sources` parameter.

## Behaviour changes

**On every law case and every series step: none in state.** The rubric 4/5 probes are below. The one measurable difference is the workbench's observed-generation counter: `resonance_zone` now reads the workbench generation once more, because it builds a trellis. Details under rubric 4.

**Outside the law cases:**
- **`near_tangency`** now holds 25 `unresolved_links` entries (e.g. `25 -> (61, 59)`, `27 -> (61, 59)`, `31 -> (67, 65)`), all sources in N. Nothing reads them.
- **A collision**, if one ever fired, would now link neither claimant. On base the first one won.
- **`resonance_zone()` at a non-candidate** now logs one WARNING and keeps the previous pip. On base it adopted the pip and `topology()` later raised `ValueError`; verified on `inversion_two_blasts` at crossing 33.
- **New API:** `forget_strong_pip`.
- **Error and warning texts** now name `forget_strong_pip`.

## Rubric, my reading

1. **Suite passes, only KNOWN_ISSUES xfail — PASS.**
   - 1168 passed, 1 xfailed, and that xfail is the single `KNOWN_ISSUES` entry.

2. **New tests pass on branch and fail on base; existing tests pass on both — PASS.**
   - All 6 new tests pass on the branch. I copied the two test modules into the base worktree and ran them there: the 6 new tests fail and the 11 existing ones pass. The base failure reasons:
     - cross-pass: `assert (None is None and 3 is None)`, i.e. `z` keeps `y`;
     - three claimants: `assert 2 is None`, i.e. the first claimant is linked;
     - multi-match: `AttributeError ... unresolved_links`;
     - zone: `0 WARNING record(s) ... expected 1`;
     - forget: `AttributeError ... forget_strong_pip`;
     - remedy: `'forget_strong_pip' in '... choose another with set_strong_pip()'` is false.
   - For the zone test, base never reaches `topology()`; a separate probe confirmed that base's `topology()` raises `ValueError` after `resonance_zone(33)`.
   - Existing tests: the base suite passes in full (1162), and so does the branch.

3. **Diff touches exactly the eight files — PASS** (stat above).

4. **Fingerprints equal on base and branch — PASS on every state aspect, with one caveat: the counter `workbench generation`.**
   - **Method.** The probe dumps `law_tier.fingerprint(case)` for each of the 11 `LAW_CASES` in base and branch processes, canonicalised: FixedPoints by label, sets sorted. The `products` entry is dropped because it is a tuple of `id()`s, which differs between any two processes.
   - **Equal on all 11 cases:** `registry generation`, `crossings`, `bridge ids`, `holes`, `pips`, `partition signature`, `alphabet`.
   - **`workbench generation` differs** on the 7 cases that build a zone: k28 one/two blasts 2→3, inversion one/two blasts 2→3, nested, nested_inner_blasts and heteroclinic 2→4. It is equal on k10, p3, inversion (unblasted) and inversion_period2.
   - **Why the counter moves.** `TangleWorkbench.generation` "counts observed changes, not mutations; compare for equality only". It advances when a READ sees a changed state. Step 4's post-zone candidacy test has to build the trimmed trellis, so it reads the generation between the zone and the next blast. One more observation means a count one higher; nested and heteroclinic make two zone calls, hence 2→4.
   - **Why it is harmless.** Every state-carrying aspect is identical. Within one process the law tier only compares the counter for equality before and after a law, which is unaffected.
   - **What it would take to remove it.** Avoiding the extra read would mean not building a trellis in `resonance_zone`, which contradicts the approved step 4. I leave it to the reviewer whether rubric 4 should be read on state or on the raw counter.

5. **`unresolved_links` empty, `+1` link count equal to base — PASS.**
   - **Law cases:** `unresolved_links == {}` on all 11. The `+1` link counts are equal on base and branch: k10 4, k28 x1 4, k28 x2 6, p3 28, nested 34, nested_inner 38, inversion 9/13/19, inversion_period2 16, heteroclinic 94.
   - **Series:** on every step of all five `cases.SERIES` (k28 9 steps, p3 10, inversion 5, k32 11, nested_outer 12), `unresolved_links == {}` both after the blast and after the step's plain recompute. The `+1` counts are equal to base at both points; for example k28 runs 2,4,6,10,18,30,50,84,140 and nested_outer runs 30 … 729.
   - **Whole steps:** every `SeriesStep` (holes, classes, letters, references, crossings, recomputed, …) is also equal to base on every step, canonicalised.

6. **Old warning gone; one "none linked" WARNING whose branch removes and records; no `continue` — PASS.**
   - `grep -n "already has the preimage" src` finds nothing.
   - `grep -n "none linked" src/tanglepack/numerics/IterateInference.py` finds line 150 only. Its `else` branch calls `registry.remove_iterate(held, 1)` when `held` is set and writes `registry.unresolved_links.update(...)`.
   - The loop is an `if`/`else` with no `continue`.

7. **`remove` exists and undoes a link — PASS.**
   - `IterateTable.remove` and `IntersectionRegistry.remove_iterate` both exist.
   - The cross-pass test does `register_iterate(z, 1, y)` and asserts the graph has edge `iter_{z}_1`. After the collision's `remove_iterate(z, 1)` it asserts that `table[z, 1]` and `table[y, -1]` are both None and that the edge is gone.

8. **Collision across passes — PASS.** This is `test_an_image_claimed_across_passes_links_neither_claimant`.
   - Fixture choice, by rule: `x -> y` is the smallest detected link, and `z` is the smallest other detected crossing.
   - The test empties the table, hand-links `z -> y`, and runs `infer_iterates(workbench, [x])`.
   - Expected and asserted: none of the three is linked, `unresolved_links == {x: (y,), z: (y,)}`, and exactly one WARNING whose first argument is `y`.
   - On base, `z` keeps `y`.

9. **Collision within one pass — PASS.** This is `test_an_image_claimed_three_times_in_one_pass_links_none`.
   - It monkeypatches `_bracketed_image` so that three sources name one image, and makes `_map_gated` the identity.
   - None of the three is linked, the image has no preimage, and all three are recorded as `(image,)`.
   - On base the first claimant is linked.

10. **Multi-match — PASS.** This is `test_a_multi_match_is_recorded_and_a_second_pass_changes_nothing`.
    - It monkeypatches `_brackets` so every crossing on the advanced branches meets.
    - Every detected source is recorded with more than one candidate and none is linked.
    - A second `infer_iterates` leaves `registry.generation` unchanged.

11. **Non-candidate zone — PASS.** This is `test_a_zone_at_a_non_candidate_keeps_the_chosen_pip`, on `build_inversion(blasts=2)`.
    - Choice, by rule: the outermost (largest stable cdist) non-synthetic non-candidate of the trellis whose stable branch differs from the chosen pip's. That is crossing 33 on `(stable, 0, 1)`; the chosen pip 1 is on `(0, 0)`.
    - Exactly one WARNING from `tanglepack.loom.TangleSession` during `resonance_zone`.
    - The zone is stored, `chosen_strong_pip(fp) is chosen`, and `topology()` runs.
    - Expected degenerate output: `topology()` logs Pseudoneighbor/Trellis WARNINGs that the trellis is too small for a reference window. Those come from other loggers, outside the call the test counts.
    - On base, `topology()` raises `ValueError`.

12. **Forget — PASS.** This is `test_a_forgotten_pip_gives_the_default_back`, on k10.
    - It runs `set_strong_pip(fp, c)`, with `c` the smallest-id candidate other than the default, then `forget_strong_pip(fp)`, then `topology()`.
    - Afterwards `trellis(fp).strong_pip == default` and `chosen_strong_pip(fp) is registry[default]`.
    - No topology runs at `c`.

13. **The error names the remedy — PASS.** This is `test_a_stale_pip_error_names_the_remedy`.
    - It reuses the existing stale-pip setup from `test_session_bridge_classes.py`, copied into `test_session_topology.py`, since that module is outside the file list: `k10_session`, then `classify_strong_pips`, then `define_resonance_zone` at `f(q0)`.
    - `topology()` raises a `ValueError` that matches "no longer a strong-pip candidate" and contains "forget_strong_pip".

14. **The law rejects bad data — PASS (one-off probe, not committed).**
    - On k10 at the bridges stage, the law passes clean (7 checked).
    - **`{z: (y,)}`**, with `z` unlinked and `y` held by `x`, is rejected: "image 3 is contested by 4 but linked to 1".
    - **`{x: (y,)}`**, with `x` linked, is rejected: "crossing 1 is unresolved but linked".
    - **`{z: ()}`** is rejected: "crossing 4 is unresolved with no candidate image".

## Probes (scratchpad, not committed)

All in `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b8/`:
- `probe_b8.py <worktree> <out.json> laws|series`: fingerprints and per-step stats. Outputs `unresolved-sets{,-base}.{laws,series}.json` and `.log`.
- `law_probe.py <worktree>`: rubric 14.
- `inv.py <worktree> [index]`: the inversion non-candidate zone, on base and branch.
- `suite.txt`: the branch suite output.

## Noticed, not done

- **`check_iterate_table_complete` still scans the whole registry, not N.** On `near_tangency` it fails at crossing 1, which is outside N, with images `[35, 36]`. This is the same on base and is the plan's "Not doing" item; it is not a law case.
- **`registry._generation += 1` is written from `IterateInference`,** as the plan specifies. That is a write to a private attribute from another module. A public bump method on the registry would stay inside the eight files, but the plan did not ask for one.
