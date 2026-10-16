# Outcome: B0 `oracle-laws` (tests only)

*Analyst (auditor), 2026-10-12. Job `outcome`. The branch is `refactor/oracle-laws` at `1029fcc`, in the worktree
`/home/dezhu/Development/TanglePack/tangle-pack-oracle-laws`. Its base is `a3bd8f0`. Review passed in one round
(`cycle_ok_test_engineer_code_reviewer.json`). There is no debrief verdict and no figure report.*

**Your working tree is dirty.** `docs/alignment_plan_2026-10-09.md`, `docs/code_fidelity_audit_2026-10-09.md`,
`artifacts/request_origin_analyst.md` and `.refactor/decisions.md` are uncommitted, plus this run's artifacts.
The branch was made from `HEAD` and touches none of them. I appended the 2026-10-12 entry to `.refactor/decisions.md`.

## The answer, now

**Yes.** Four new laws check the definitions with computations that never call the code they guard.
Three of them run on all 11 law cases; the fourth is the k=10 pin. All four pass on the unchanged library.

| Law (`tests/helpers/laws.py`) | What it re-derives, and from what | Where it runs | What turned it red |
|---|---|---|---|
| `chosen_pip_is_strong` (`:1566`) | The open-box test (A.4), from the registry and `FixedPoint` only. It scales by `per_step_beta` and never reads the iterate table. | All 11 law cases. It checks 5 to 108 crossings per case. | `is_strong_pip` made to accept everything, plus a forced pip it would reject (k10, pips 5 and 6). |
| `reference_pairs_clear_of_X` (`:1632`) | Both open intervals of every reference pair (the pair definition), by scaling only. It uses each interval's own cycle and β (D7). | The 11 law cases only (`LAW_CASES_ONLY`). | `_unstable_interval_is_clear` forced to True: red on 7 of 11 cases. |
| `propagated_hole_inside_container` (`:1672`) | The reference bridge's ends, carried back by the **real** `map_inv`, must land inside the container. | All 11 law cases. | Container shifted to the next bridge: red on 9 of 11, and vacuous (so also red) on the other 2. |
| `k10_is_the_full_shift` (`:3024`) | Q11. Letter-free: one active class whose word is itself twice in opposite directions, class matrix `[[2]]`, refined matrix `[[1,1],[1,1]]`. | k10. The other 10 cases are NOT_APPLICABLE. | Two breaks: `refine` not splitting, and `transition_matrix` weight forced to 1. |

Check 4 is the first law that guards the **class** matrix: `matrix_is_token_counts` reads only the refined one.

## What changed

No `src/` file changed (`git diff a3bd8f0 -- src` is empty). `--stat` lists exactly these six files (*verified*):

| File | + | − |
|---|---|---|
| `tests/helpers/laws.py` | 169 | 0 |
| `tests/cases.py` | 15 | 9 |
| `tests/regression/test_open_p3_deep_runs.py` | 9 | 6 |
| `tests/conftest.py` | 3 | 0 |
| `tests/invariants/test_law_partition.py` | 3 | 0 |
| `tests/invariants/test_law_symbolic.py` | 2 | 1 |
| **Total** | **201** | **16** |

Nothing got shorter. This slice only adds. The most important hunk is the containment test with its shared-end
closure (`laws.py:1689-1701`).

**Before:** no law checked that a propagated hole's container actually holds the carried reference bridge. The
backward hole chain was guarded only by laws that read the chain's own output.

**After:**

```python
for end in hole.origin:
    if _table_image(registry, end, hole.iterate) in hole.bounding_ids:
        count += 1          # the preimage IS the container's end crossing: shared end
        continue
    point = registry[end].coords; for _ in range(-hole.iterate): point = map_inv(point)
    ... nearest segment of the container's node polyline, cdist interpolated ...
    assert lo < cdist < hi
```

## Behaviour changes

- **Library:** none.
- **`build_inversion_period2` now defaults to `stable_steps=10`** (Q12, `tests/cases.py:825`). It gives 6 classes,
  10 holes and `is_reliable` True. *Verified:* I built it on the branch. The 9/9 build had the same counts but
  `is_reliable` False.
- **The suite:**
  - Four laws were added inside the existing layer tests.
  - There are 11 new visible set-aside skips: `k10_is_the_full_shift` on the 10 other cases, and
    `every_active_class_walks` on `inversion_period2`.
  - `KNOWN_ISSUES` is unchanged, still only C3.
- **The deep-run test no longer runs the two `LAW_CASES_ONLY` laws** on its five builds.

## What was verified and how

- **Branch suite (*verified*, my own run).** `MPLBACKEND=Agg PYTHONPATH=src ... pytest -q -rxXs` gave **1117 passed,
  38 skipped, 5 deselected, 1 xfailed**. I first checked that `tanglepack.__file__` resolves to the worktree's `src`.
- **Base suite (reported by the test engineer).** 1117 / 27 / 5 / 1. The 11 extra collected ids are exactly the 11
  new set-aside parameters.
- **Base library with the branch's tests (*verified*).** The library is byte-identical, since the `src` diff is
  empty, so my branch run is that run. The test engineer also checked out the branch's `tests/` on the base
  worktree and got 36 passed, 32 skipped.
- **Independence (*verified*).** I grepped the bodies of the four checks and five helpers (`:1566-1702`,
  `:3024-3041`). None of `StrongPip`, `is_strong_pip`, `Pseudoneighbor`, `_unstable_interval_is_clear`,
  `StablePartition`, `_containing_bridge`, `bridge_for_pair`, `image_cdist` or `import` appears. Check 1 has no
  `iterate`.
- **Break-it table (reported).** Run by the test engineer with scratch edits, then reverted; the reviewer re-ran
  checks 1, 3 and 4 independently. The results are in the table above.
- **Shared-end closure switched off (reported, re-run by the reviewer by monkeypatch).** The law is then red on
  two cases:
  - `k28_one_blast`: 42.250102 against `hi` 42.250098.
  - `inversion`: 0.88968 against 0.88903.

  With the closure on, all 11 cases pass. `inversion_period2` at 9/9 also needs the closure; at 9/10 it does not.

## What did not go to plan

- **The base tally is 1117, not 1122.** The 1122 figure dates from `69371fd`, before your D15 commit `74a19da`.
- **The scope widened to a sixth file,** `test_open_p3_deep_runs.py`. Its law loop read `LAYERS` and failed the k10
  pin as vacuous on the deep builds. The fix was `LAW_CASES_ONLY`, which keeps "no deep-series oracle" true.
  - A side effect: the deep-run loop no longer goes through `run_layer`. It now stops at the first failing law
    and bypasses `NOT_APPLICABLE` / `KNOWN_ISSUES`.
  - None of those entries exist for the deep builds today, so a future one would be silently ignored there
    (`test_open_p3_deep_runs.py:94-98`).
- **Check 1's break shows only through the forced-pip driver.** The default pip policy `min(unstable_cdist)` still
  picks a truly strong crossing on every law case, so `is_strong_pip` accepting everything goes unnoticed in the
  law tier.
- **Check 2's break is invisible on 4 cases** (`k28_one_blast`, `p3`, `inversion`, `inversion_period2`). On those,
  every consecutive window pair is already clear.
- **Small deviation from the plan in check 3.** An end that the table links to a NON-container crossing is tested
  geometrically, not failed outright. The reviewer judged this no looser; no case exercises it.

## Where I am least sure it matches what you meant

1. **The shared-end closure in check 3.** It turns "strictly between" into "between, or exactly at a container
   crossing the table names". I believe that is physics: the preimage is a crossing of N. But it is a relaxation
   of A.6's text that you have not ruled on.
2. **Check 3 bounds no distance from the carried end to the container polyline.** A wrong container on another
   branch could still take an interior projection and pass. See the measurement below. Whether to close this
   gap is your call; nothing was changed.
3. **The deep-run loop bypassing `run_layer`** (above). It may be acceptable as a cost of keeping the file list
   short, or you may prefer the `LAW_CASES_ONLY` filter inside `run_layer` (`tests/helpers/law_tier.py:139`, a
   seventh file).

### The distance gap, measured

The orchestrator measured this on the branch on 2026-10-12, using the real inverse map. Each row gives the
maximum nearest-distance from a carried end to its container polyline, and the same divided by the polyline
length. Only ends not closed by the shared-end rule are counted.

| Case | Max distance | ÷ polyline length |
|---|---|---|
| k10 | 9.0e-3 | 1.1e-3 (2 steps) |
| k28_one_blast | 4.0e-6 | 2.4e-7 |
| k28_two_blasts | 4.7e-5 | 2.9e-6 |
| p3 / nested / nested_inner_blasts | 1.7e-5 | 4.2e-5 |
| inversion | 3.6e-4 | 4.5e-4 |
| inversion_one / two_blasts | 1.1e-3 | 1.4e-3 (6 steps) |
| inversion_period2 | 2.0e-4 | 1.5e-3 |
| heteroclinic | 2.6e-4 | 6.1e-3 (up to 26 steps) |

There are two ways to close the gap:

- **A distance bound.** It would need a tuned constant of order **1e-2 of the polyline length**.
- **A tolerance-free test:** "the container is the NEAREST bridge to the carried end". You may prefer this or
  reject it.

## Optional touch-up (reviewer's notes, not rubric)

- Blank lines: `tests/helpers/laws.py:3021-3023` has three before `check_k10_is_the_full_shift`, and `:3042` has
  one before the `#: The read-only law checks` comment. PEP 8 wants two in both places.
- `_scaled` (`:1601`) is called once and could be inlined. The nested `any(all(...))` in check 2 (`:1653-1660`)
  would read better with a named local, for example `is_own_end`.
- The open distance gap in check 3 (above).

## House loop

Merging `refactor/oracle-laws` into `main` is yours to do. Your debrief answer goes in
`artifacts/debrief_ok_auditor_human.json`, as one of:

- `"confirmed"`
- `"touchup"` (for example the optional list above)
- `"figures"`

## Files written

- `.refactor/decisions.md`: the 2026-10-12 B0 entry.
- `.refactor/knowledge/testing.md`: the header and a new B0 section covering the laws, `LAW_CASES_ONLY`, the Q11
  exception, `inversion_period2` at 10 steps, the teeth and the smells.
- `.refactor/knowledge/topology.md`: an oracle-laws note and the A.0/A.6 correction.
