# Change report: B6b `word-finished` (coder, iteration 1)

**Test engineer:** this report is your input. `code_reviewer` grades your work afterwards as
"pass" or "revise", against the plan's eleven-item rubric (`artifacts/plan_planner_human_2.md`
§7). I implemented the plan's §4 steps exactly as the author approved them
(`plan_ok_planner_human.json`, all six questions at their defaults), in the ten files of §3.
Nothing under `docs/`, `.refactor/`, `CLAUDE.md` or `artifacts/` was touched, apart from this
report.

## Where it is

- **Base commit:** `63477e3b587eb7c147fed19d68e11bac5e5d81ee` (main HEAD, B7b merged).
- **Branch:** `refactor/word-finished`, in the worktree
  `/home/dezhu/Development/TanglePack/tangle-pack-word-finished`. HEAD is `3486039`.
- **Base worktree:** `/home/dezhu/Development/TanglePack/tangle-pack-word-finished-base`,
  detached at `63477e3` and clean.
- The main checkout is untouched: still on `main` at `63477e3`, with the same uncommitted docs as
  at the start.

## Commits (`git log --oneline 63477e3..refactor/word-finished`)

All three are green on their own. Each intermediate tree was run through the full suite in a
scratch copy:

| Commit | Message | Suite |
|---|---|---|
| `cc603d4` | Dual walk: three statuses and a reason; the truncated flag goes (B6b) | 1218 passed |
| `01bc087` | Symbolic dynamics: free reduction, one function for both readers (B6b, A.8c) | 1221 passed |
| `3486039` | Diagnostics: the evidence comparison moves out of the algorithm, verify_words (B6b) | 1221 passed |

Each run also had 47 skipped, 5 deselected and 1 xfailed (C3). Every commit ends with
`Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`.

## `git diff --stat 63477e3..refactor/word-finished`

```
 src/tanglepack/topology/DualWalk.py           | 143 +++++-----
 src/tanglepack/topology/SymbolicDynamics.py   | 367 +++++---------------------
 src/tanglepack/topology/__init__.py           |  12 +-
 src/tanglepack/topology/diagnostics.py        | 285 +++++++++++++++++++-
 src/tanglepack/topology/plotting.py           |   8 +-
 tests/helpers/laws.py                         |  18 +-
 tests/regression/test_evidence_diagnosis.py   |  12 +-
 tests/regression/test_heteroclinic.py         |   4 +-
 tests/unit/topology/test_dual_walk.py         | 107 ++++++--
 tests/unit/topology/test_symbolic_dynamics.py |  77 ++++--
 10 files changed, 598 insertions(+), 435 deletions(-)
```

These are exactly the ten files of plan §3.

## Suite

The command, run from inside each worktree, was
`MPLBACKEND=Agg PYTHONPATH=src ../tangle-pack/env/bin/python -m pytest -q -rxX -p no:cacheprovider`.

- **Branch (HEAD `3486039`):** 1221 passed, 47 skipped, 5 deselected, 1 xfailed in 201.67 s.
  - The only `-rxX` line is
    `XFAIL tests/invariants/test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing] - C3: ...`.
  - There is no XPASS.
- **Base (`63477e3`):** 1211 passed, 47 skipped, 5 deselected, 1 xfailed in 206.38 s. The same
  single C3 xfail.
- The +10 are new tests:
  - `test_dual_walk.py`: 6 parametrized reason cases plus the cap/tie test;
  - `test_symbolic_dynamics.py`: `free_reduce`, `c c^-1` gives the empty word, and the cancelled
    pair adds no class.

## What changed, step by step

**Step 1: statuses (`DualWalk.py`).**
- `WalkStatus = Literal["resolved", "ambiguous", "unresolved"]`, and `WalkReason` holds the six
  reasons.
- `WalkSearch.reason: Optional[WalkReason]` is set iff the status is `unresolved`. The
  `truncated` field is gone.
- Return sites:
  - the six failures return `"unresolved"` with their reason;
  - a zero-crossing walk and a unique walk return `"resolved"`;
  - a tie returns `"ambiguous"`;
  - a search cut short by `max_walks` returns `"unresolved"`/`"truncated"` and KEEPS its `walks`;
  - the constrained out-of-budget case returns `"unresolved"`/`"truncated"` with no walk.
- Positional `truncated` arguments became keywords. `__repr__` prints the reason.
- The module docstring opens with A.8c's two definitions and lists what leaves a class
  unresolved. The `truncated` Dev Note is rewritten.

**Step 2: one reduction (`SymbolicDynamics.free_reduce`).** It is a stack pass over
`(class, direction)` tokens and returns the kept indices. It has three readers:
1. the closure and inertness `words` in the driver (a class named only inside a cancelled pair is
   no longer added);
2. `cd.symbols` (`[symbols[i] for i in free_reduce(_signature(symbols))]`);
3. the registered evidence word in `registered_evidence`, so a word is compared with a word
   (question 3). `ImageEvidence.itinerary` stays the raw record.

`translate_itinerary` stays unreduced, so the T law keeps its same-occurrence reading.

**Step 3: readers of the new form.**
- `_resolve_itinerary` sets `cd.unresolved_reason` from `search.status == "unresolved"`. It gives
  the budget sentence for reason `truncated`, else `f"dual-graph walk {search.reason}"`.
- `cd.ambiguous = search.status == "ambiguous"`.
- The wall gate is `search and search.status == "resolved" and search.walk.multiplicity == 1`.
- `MembersDisagree` is unchanged and lives on the class only (question 4).

**Step 4: the move.**
- `ImageEvidence`, `EvidenceDiagnosis`, `registered_evidence`, `_diagnose` and `_order_flip` move
  to `diagnostics.py`, together with the evidence bullet of the module docstring (now a short
  `**Words.**` bullet).
- `_compare_evidence` (called once) is folded into the new entry point. Its signature is trimmed
  to what it reads:

  ```python
  verify_words(records, trellis, naming, table, letters)
  ```

  It reads `naming.iterated` and `naming.homotopy`. The driver calls it once, after all symbols
  are final.
- Imports:
  - `diagnostics` imports `_image_chain` and `trellis_itinerary` at module level;
  - `translate_itinerary`, `_signature` and `free_reduce` are imported inside `registered_evidence`
    and `verify_words`;
  - type names sit under `TYPE_CHECKING`.
- `SymbolicDynamics` imports `verify_words` at module level, and `ImageEvidence` under
  `TYPE_CHECKING` for the `ClassDynamics.evidence` annotation.
- `_member_pair` stays in `SymbolicDynamics.py`.
- `_diagnose` now takes the first-difference `position` from its caller, so it needs no
  `_signature` import of its own.

**Step 5: docstrings and plotting.** `plotting._class_status` returns `unresolved: <reason>`, else
`trellis`, else `cd.status` (`resolved` / `trivial` / `ambiguous`).

**Step 6: tests.** Pins are rewritten, never weakened. The pin-level details are listed under
rubric 2.

**`__init__.py`.**
- The moved names are re-exported from `.diagnostics`.
- `MembersDisagree` is exported.
- I also exported `WalkReason` and `free_reduce`. They are the public siblings of the already
  exported `WalkStatus` and `translate_itinerary`; the plan did not list them, so this is
  flagged here as a small addition.

## Behaviour changes (all of them)

No number changes on any of the 63 items (rubric 4 below). What reads differently:

1. **`WalkSearch.status` values.**
   - `unique` and zero-crossing walks read `resolved`.
   - Failures read `unresolved` plus a `reason`.
   - A capped search reads `unresolved`/`truncated` and keeps its walks. On base it read
     `unique`/`ambiguous` with `truncated=True`.
   - `WalkSearch.truncated` no longer exists.
2. **The plot's class status label.** It shows the class status. A one-token zero-crossing class
   that showed "trivial" now shows "resolved"; "trivial" appears only for an empty word. The `-`
   fallback is gone: a hand-built class with no search and no reason now shows its `cd.status`.
3. **A constrained search that ran out with no walk** now gets the budget sentence as its
   `unresolved_reason`; on base it read "dual-graph walk truncated". It never occurs on the corpus.
4. **Free reduction on hand-built data.** `a a^-1` is now the empty word (class `trivial`), and a
   class named only inside a cancelled pair is not added to the closure. Corpus: 0 cancellations.
5. **Logger.** Evidence WARNING/INFO lines now come from `tanglepack.topology.diagnostics`.
6. **Two log-only changes.** I did not list these in the plan:
   - (a) The per-class DEBUG line `"class %s: %s, word %r, verified %s"` of the driver is gone. It
     would have needed a third loop, and `diagnostics.py` sits exactly at its cap.
   - (b) The `max_walks` WARNING now ends "..., the search decides nothing". A capped search with
     several walks no longer also logs the "ambiguous" WARNING, because it is unresolved, not
     ambiguous.
7. **Exports.** `MembersDisagree`, `WalkReason` and `free_reduce` are importable from
   `tanglepack.topology`.

## Rubric, item by item (my own reading)

**1. PASS.** The branch passes: 1221 passed, 47 skipped, 5 deselected. `-rxX` lists exactly one
xfail (C3) and no XPASS.

**2. PASS.** The base suite passes in full (1211 passed). The rewritten pins (step 6) pass on base
in their old form, because they are the base tests. The new form of each:
- `test_dual_walk.py`:
  - `status == "trivial"/"unique"` → `"resolved"` (plus `reason is None` where the old pin said
    `not .truncated`);
  - the failure pins → `status == "unresolved" and reason == <value>`;
  - `:254-260` → `full.reason is None`, and capped reads `unresolved`/`truncated` with 4 walks
    kept;
  - `:450` → `unresolved and reason == "truncated" and walk is None`.
  - Two test functions were renamed to drop the old vocabulary:
    `test_trivial_walks_stay_on_the_start_face` → `test_zero_crossing_walks_stay_on_the_start_face`,
    and `..._dedupe_to_unique` → `..._dedupe_to_one_walk`.
- `laws.py` `:2775, :2817, :2932` and `test_heteroclinic.py:76, :94`: `in ("unique", "trivial")`
  → `== "resolved"`.
- `test_symbolic_dynamics.py`:
  - `:754` → `status == "unresolved" and reason == "truncated"`, with `walk is not None`;
  - `:831` → filters on `reason == "truncated"`, asserts the `unresolved` status, and asserts
    "budget" in the reason UNCONDITIONALLY. This is stricter than base, which only checked when
    a walk was kept; it follows behaviour change 3.
- **Fixtures.**
  - The loops test (`:214`) still pins `translate_itinerary` → `["a", "a^-1"]`. That function
    stays unreduced by design (question 3). It now ALSO asserts that
    `free_reduce(...) == []`, the empty word.
    - **This deviates from the plan's literal "now expects the empty word".** That wording
      assumed the reduction would sit inside `translate_itinerary`.
    - The driver-level empty word is pinned by the new rubric-8 test.
  - The evidence fixture (`:656`) is `a a` (last two refs swapped). It keeps position 1 and the
    same `walked`/`registered` occurrences.
- `_drive` patches `diagnostics._image_chain` / `diagnostics.trellis_itinerary`. Its fake naming
  gets the `homotopy`/`iterated` attributes the real `ElementNaming` has.
- The two `registered_evidence` unit tests are retargeted to `diagnostics`.
- The logger name is `diagnostics` in `:673` and in `test_evidence_diagnosis.py`. The latter now
  also filters its WARNING count by `r.name == "tanglepack.topology.diagnostics"`.

**3. PASS.** `git diff --stat` touches only the ten files of §3 (above).

**4. PASS.** The probe ran on all 63 items: the 16 builds, plus every step of the 5 `cases.SERIES`
(9 + 10 + 5 + 11 + 12 = 47).
- Each series step is blasted as `cases.blast_series` does: `cases._blast(session, largest zone,
  [its fp], diagnose=False, min_separation=series.min_separation)`.
- It ran once against the base worktree and once against the committed branch. A cross-check
  asserts `tanglepack.__file__` lies under the given repo.
- **Equal on all 63, 0 mismatches:**
  - letter-free `rules` and `refined_rules` (tokens as `(class label, direction, refined index)`);
  - `transition_matrix()` at both levels (node names letter-free);
  - `is_reliable`, every `verified` and every `itinerary_verified`;
  - every `ClassDynamics.status` and every `unresolved_reason`;
  - the closed-table size, the inert set, and every evidence word.
- **Walk statuses**, base `(status, truncated)` → branch `(status, reason)`, per class over all
  63 items:

  ```
  ('trivial', False)  -> ('resolved', None)   402
  ('unique', False)   -> ('resolved', None)   151
  None (no search)    -> None                 349
  mapping violations: none
  ```

  Every branch value equals the plan's mapping of the base value. The number of walks per search
  is equal too. No `ambiguous`, no failure value and no `truncated=True` occurs on the corpus. As
  the analyst found, the failure rows of the mapping are exercised only by the unit tests.

**5. PASS, with one remark on the first grep.**

```
$ grep -n "EvidenceDiagnosis\|ImageEvidence\|_order_flip\|_diagnose\|_compare_evidence\|registered_evidence" src/tanglepack/topology/SymbolicDynamics.py
188:    from .diagnostics import ImageEvidence
311:    evidence: list[ImageEvidence] = field(default_factory=list)
```

- Line 188 is the `TYPE_CHECKING` import. Line 311 is the `ClassDynamics.evidence` annotation the
  plan itself names ("and `ImageEvidence` for its annotation"); it is not evidence code.
- The single `verify_words(...)` call is at the end of the driver's translate loop. The grep
  pattern does not match it.
- The other checks:
  - `grep -n "^from .SymbolicDynamics" src/tanglepack/topology/diagnostics.py` finds nothing;
  - `_member_pair` is still defined in `SymbolicDynamics.py` (`:638`);
  - `ImageEvidence.__module__ == "tanglepack.topology.diagnostics"`;
  - `from tanglepack.topology import ImageEvidence, EvidenceDiagnosis, registered_evidence, MembersDisagree`
    succeeds;
  - `import tanglepack` succeeds in a fresh interpreter (cwd `/tmp`), as does a fresh import of
    each of `topology.diagnostics`, `topology.SymbolicDynamics`, `topology.DualWalk` and
    `tanglepack.loom` alone.

**6. PASS.** All three counts are within their caps:

| File | `wc -l` | Cap |
|---|---|---|
| `SymbolicDynamics.py` | **1,312** | ≤ 1,330 |
| `DualWalk.py` | **1,383** | ≤ 1,390 |
| `diagnostics.py` | **990** | ≤ 990, exactly at the cap |

`diagnostics.py` reached 990 only after `_compare_evidence` was folded into `verify_words` and
two log calls were compacted. **Any later addition there must remove a line elsewhere.**

**7. PASS, with one remark on the second grep.**
- `grep -rn '"unique"\|"trivial"' src/tanglepack/topology/DualWalk.py` finds nothing.
- `grep -rn "search\.truncated\|\.truncated and\|not .*\.truncated" src tests` finds exactly one
  line: `DualWalk.py:1214  if not found_walks and budget.truncated:`. That is `_Budget.truncated`,
  which the rubric explicitly exempts; the line is also present on base.
- The new parametrized `test_every_failure_is_unresolved_with_its_reason` builds on the `tie`
  graph. It reaches all six reasons with `status == "unresolved"` (ids `no_start_node`,
  `no_goal_node`, `span_meets_no_node`, `unreachable`, `no_admissible_walk`, `truncated`).
- The new `test_a_cap_cutting_a_tie_short_keeps_its_candidate_unresolved` runs the tie at
  `max_walks=1`. It asserts `unresolved`/`truncated`/`walk is not None` and FAILS on base.
  - **On base the status there is `'unique'` (with `truncated=True`), not `ambiguous` as the
    rubric text says.** The budget refuses the second face path before it becomes a second
    itinerary.

**8. PASS.**
- `test_free_reduction_cancels_adjacent_inverses_only`:
  - `a, b, b^-1, a^-1` → `[]`;
  - `a, b, a^-1` → `[0, 1, 2]`;
  - `a^-1, a, b` → `[2]`.
- `test_a_walk_reading_c_cinv_is_the_empty_word` (driver): the walk is `(R_1^1,R_3^3,R_3^1,R_1^3)`,
  which reads `a a^-1`. The result is `cd.symbols == []`, `cd.status == "trivial"` and
  `cd.loops == []`.
  - Run against base with `_drive`'s patches pointed back at `SD`, it fails with
    `assert [Symbol(...), Symbol(...)] == []`.

**9. PASS.** `test_a_class_named_only_inside_a_cancelled_pair_is_not_added`:
- The alphabet is `{a}`. `a` walks `L_1^2 L_3 | L_3 L_1^2`, which reads `u u^-1`, and `u` is
  outside the alphabet. The closed table stays `[a]`, and `a` is trivial.
- `u` is given a walk of its own so that base fails on the assertion itself:
  `assert [a, u] == [a]` ("Left contains one more item"), not on a KeyError.

**10. PASS.** `test_evidence_disagreement_names_a_stable_end` still asserts:
- `position == 1`;
- `walked == (L_1^2, L_3)` and `registered == (R_1^3, R_3^1)`;
- one WARNING for "walk" and none for "stable end", counted on `DG.__name__`
  (`tanglepack.topology.diagnostics`).

On a scratch copy of the branch with the two `_drive` patch targets reverted to
`SD._image_chain` / `SD.trellis_itinerary`, both parametrizations FAIL with
`AttributeError: '_Crossings' object has no attribute 'iterate'`. The real lookups run, so the
patch demonstrably reaches the evidence path.

**11. PASS: 0 cancellations, the slice continues.** On the branch, for every item, each class
word and each evidence reading is re-translated UNREDUCED from its itinerary
(`translate_itinerary(cd.itinerary)` / `translate_itinerary(ev.itinerary)`). It is then reduced
by three counters:
- the library's `SD.free_reduce` (`lib`);
- an independent stack oracle at class level (`oracle`);
- the oracle at same occurrence (`occ`).

The base run applied the same two oracles.

```
reduction totals (branch): class_words 899, class_lib 0, class_oracle 0, class_occ 0,
                           evidence 1830, ev_lib 0, ev_oracle 0, ev_occ 0
reduction totals (base):   class_oracle 0, class_occ 0, ev_oracle 0, ev_occ 0
CANCELLATION lines: none
```

- 899 class words = the 902 classes minus the 3 unresolved ones (`k28_series@0`, `k32_series@0`,
  `nested_outer_series@0`, singletons with no chain). The 1,830 evidence readings match the
  analyst's count.
- Per item, every `lib/oracle/occ` column is 0 on all 63. The full table is in
  `compare_final.txt` (scripts below).

## Scripts (for the test engineer to re-run, not to trust)

All are under
`/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/coder/`.

- `probe.py REPO build|series NAME`: one JSON line per item, with all the rubric-4 fields and the
  rubric-11 counters.
  - It imports `SymbolicDynamics` via `importlib.import_module`.
  - **Note:** `from tanglepack.topology import SymbolicDynamics` gives the CLASS, not the module.
- `run.sh SIDE REPO`: all 63 items in parallel, into `out/SIDE/`. It takes about 2 minutes.
- `compare.py`: equality, the walk-status mapping table, and the reduction totals and per-item
  table. Its output is in `compare_final.txt`.
- `states.py` and `states/{A,B}/`: the intermediate file states the first two commits were built
  from. `tree_A/` and `tree_B/` with `suite_A.txt` / `suite_B.txt` hold their full-suite runs.
- Base and branch suite logs: `../base_suite.txt`, `suite_final.txt`.

**pytest's own `pythonpath = ["src", "tests"]`** (pyproject) puts the rootdir's `src` first. To
run branch tests against base code, use a scratch copy of the base tree with the test files laid
over it. A `PYTHONPATH` override does not work. `../basecopy/` is such a copy, with
`_drive`'s patches pointed at `SD`.

## Not done / open

- The stale `CLAUDE.md` symbolic-dynamics bullet (B12; plan §6). It still lists the old statuses,
  `WalkSearch.truncated` and "falls back to all of them".
- `scripts/symbolic_figures.py:268` prints the new walk values unedited, as planned.
- `stress/numerics/report/make_numbers.py` was not touched. It reads archived reason strings, and
  `"dual-graph walk unreachable"` keeps its text.
