# Test evidence — B6b `word-finished` (test engineer, round 1)

**Reviewer:** you read this evidence and grade the slice **"pass" or "revise"** against the
eleven-item rubric of `artifacts/plan_planner_human_2.md` §7. Every number below was re-measured
by me; nothing is copied from the coder's report without a rerun.

## Where

- Base: `63477e3` in `/home/dezhu/Development/TanglePack/tangle-pack-word-finished-base` (detached,
  left clean: `git status --short` empty).
- Branch: `refactor/word-finished` in `/home/dezhu/Development/TanglePack/tangle-pack-word-finished`.
  The coder's HEAD was `3486039`. I added ONE test commit, so HEAD is now **`65edec4`**.
- The main checkout is untouched: still `main` @ `63477e3`, with the same uncommitted
  `docs/`, `.refactor/`, `CLAUDE.md`.
- My scratch area is `.../scratchpad/te_b6b/`:
  - `basecopy/`: a `git archive 63477e3` tree with the branch tests laid over it;
  - `brcopy/`: a `git archive` of the branch, git-initialised, used for the break-it runs;
  - `out/`: fresh probe output;
  - `compare.txt`, `breaks.txt`, `breaks2.txt`, `suite_*.txt`.

## Rubric 1: full suite

Command, run from inside each worktree:
`MPLBACKEND=Agg PYTHONPATH=src ../tangle-pack/env/bin/python -m pytest -q -rxX -p no:cacheprovider`

| Tree | Result |
|---|---|
| base `63477e3` | 1211 passed, 47 skipped, 5 deselected, 1 xfailed (211.6 s) |
| branch `3486039` (the coder's) | 1221 passed, 47 skipped, 5 deselected, 1 xfailed (211.5 s) |
| branch `65edec4` (+ my 2 tests) | **1223 passed**, 47 skipped, 5 deselected, 1 xfailed (202.0 s) |

- On every run the only `-rxX` line is
  `XFAIL tests/invariants/test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing] - C3: ...`.
- No XPASS, no FAILED, no ERROR.

**PASS.**

## Rubric 2: the pins pass on base

- The base suite passes in full (1211). Every pin rewritten in step 6 is a base test in its old
  form, so it passes on base.
- **Branch test files over base code.** I copied `test_dual_walk.py` and
  `test_symbolic_dynamics.py` into `basecopy/`. The patches that must hit base code were pointed
  back at `SD`: `_drive`'s `_image_chain`/`trellis_itinerary`, the two `registered_evidence` unit
  tests, and the logger names. I ran from `basecopy/` (its own `src` is first on pyproject's
  `pythonpath`, and `tanglepack.__file__` is under `basecopy/src`). Results:
  - `test_evidence_disagreement_names_a_stable_end[stable end]` and `[walk]`, with the new
    `a a` fixture: **2 passed on base**. The fixture change keeps position 1 and the same
    occurrences on the old code.
  - The base-failing set is listed under rubric 7-9.

**PASS.**

## Rubric 3: diff stat

`git diff --stat 63477e3..refactor/word-finished`:
- at `3486039`: 10 files, 598 insertions(+), 435 deletions(-);
- at `65edec4`: 10 files, 658 insertions(+), 435 deletions(-). My commit touches only
  `tests/unit/topology/test_symbolic_dynamics.py`, which is already in the list.

These are exactly plan §3's ten files.

**PASS.**

## Rubric 4: base-vs-branch equality on 63 items (re-run)

**How I ran it.**
- I copied the coder's `probe.py`, `run.sh` and `compare.py` into `te_b6b/` and re-pointed
  `run.sh` at a fresh `te_b6b/out/`.
- I read the probe before trusting it:
  - it asserts `tanglepack.__file__` lies under the given repo;
  - it imports the `SymbolicDynamics` module via importlib;
  - it blasts each series step as `cases.blast_series` does.
- I added a second independent reduction oracle (`naive_reduce`: delete the first adjacent
  inverse pair, repeat) and a `class_rawlen_eq_word` counter. The counter checks that the
  unreduced translation has the same length as `cd.symbols`.
- Run: `run.sh base <base worktree>` and `run.sh branch <branch worktree>`, both against
  `3486039`'s code. My later commit only adds tests.

**Checks on the run itself.**
- 63 items on each side: 16 builds plus 47 series steps (inversion 5, k28 9, k32 11,
  nested_outer 12, p3 10).
- All 63 base lines carry `"branch": false` and all 63 branch lines `"branch": true`.
- Every `.err` file is empty.

**Results** (`te_b6b/compare.txt`):

```
items: 63 63 same keys: True
equality mismatches: none
walk-status mapping (base (status, truncated) -> branch (status, reason)): count
  ('trivial', False)  -> ('resolved', None)  402
  ('unique', False)   -> ('resolved', None)  151
  None                -> None                349
mapping violations: none
```

- **Equal on all 63:**
  - letter-free `rules` and `refined_rules`;
  - `transition_matrix()` at both levels;
  - `is_reliable`, every `verified` and every `itinerary_verified`;
  - every `ClassDynamics.status` and every `unresolved_reason`;
  - closed-table size, inert set, and every evidence word.
- **Walk statuses.**
  - 553 searches, all mapped per the plan, with an equal walk count per search.
  - The corpus has no `ambiguous`, no failure value and no `truncated=True`. The failure,
    ambiguous and truncated rows of the mapping are therefore exercised only by unit tests
    (rubric 7 and my ambiguous test below).

**PASS.**

## Rubric 5: the move is complete

Re-run on the branch:

```
$ grep -n "EvidenceDiagnosis\|ImageEvidence\|_order_flip\|_diagnose\|_compare_evidence\|registered_evidence" src/tanglepack/topology/SymbolicDynamics.py
188:    from .diagnostics import ImageEvidence          (TYPE_CHECKING)
311:    evidence: list[ImageEvidence] = field(default_factory=list)
$ grep -n "^from .SymbolicDynamics" src/tanglepack/topology/diagnostics.py     -> nothing (rc 1)
```

- `diagnostics.py` imports `SymbolicDynamics` only under `TYPE_CHECKING` (`:57`) and inside
  functions (`:827` and `:932`).
- It imports `_image_chain` (`:45`, from `BridgeClass`) and `trellis_itinerary` (`:46`) at
  module level.
- The single `verify_words(...)` call is at `SymbolicDynamics.py:1301`, with the module-level
  import at `:179`.
- `def _member_pair` is still in `SymbolicDynamics.py:638`.
- In a fresh interpreter (cwd `/tmp`, `PYTHONPATH=<branch>/src`):
  - `ImageEvidence.__module__`, `EvidenceDiagnosis.__module__` and
    `registered_evidence.__module__` are all `tanglepack.topology.diagnostics`;
  - `from tanglepack.topology import ImageEvidence, EvidenceDiagnosis, registered_evidence, MembersDisagree`
    succeeds;
  - `import tanglepack` succeeds, and so do separate fresh imports of `topology.diagnostics`,
    `topology.SymbolicDynamics`, `topology.DualWalk` and `tanglepack.loom`.
- **Letter deviation:** the grep hits the annotation at `:311`, which the plan's "import line(s)
  and the one `verify_words` call" does not list. The plan's own step 4 asks for that import
  "for its annotation". The grep pattern cannot match the `verify_words` call.

**PASS** (letter deviation noted).

## Rubric 6: line counts (`wc -l`)

| File | Branch | Cap | Base |
|---|---|---|---|
| `SymbolicDynamics.py` | **1312** | ≤ 1330 | 1537 |
| `DualWalk.py` | **1383** | ≤ 1390 | 1370 |
| `diagnostics.py` | **990** | ≤ 990, exactly at the cap | 713 |

**PASS.** `diagnostics.py` has no headroom left.

## Rubric 7: three statuses

- `grep -rn '"unique"\|"trivial"' src/tanglepack/topology/DualWalk.py` finds nothing (rc 1).
- `grep -rn "search\.truncated\|\.truncated and\|not .*\.truncated" src tests` finds one line,
  `DualWalk.py:1214  if not found_walks and budget.truncated:`.
  - That is `_Budget.truncated`, which the rubric exempts.
  - The same grep finds 9 hits on base.
  - No other `.truncated` reader is left in `src`, `tests` or `scripts`.
- `test_every_failure_is_unresolved_with_its_reason`, six ids on the `tie` graph: passes on the
  branch. On base all six FAIL, each reading the old status (`'no_start_node' == 'unresolved'`,
  …, `'unique' == 'unresolved'` for `truncated`).
- `test_a_cap_cutting_a_tie_short_keeps_its_candidate_unresolved`: passes on the branch and
  FAILS on base with `assert ('unique' == 'unresolved'`.
  - I confirmed the base value directly. The tie at `max_walks=1` gives
    `status='unique', truncated=True, 1 walk`; the full tie is `ambiguous, truncated=False,
    2 walks`.
  - **The rubric's "on base the status is `ambiguous`" is wrong in its letter, as the coder
    said.** The test still fails on base, which is what the item requires.

**PASS** (letter deviation in the rubric text noted).

## Rubrics 8-9: free reduction

All tests below were run against base code through `basecopy/`.

| Test | Branch | Base |
|---|---|---|
| `test_free_reduction_cancels_adjacent_inverses_only` (`a b b⁻¹ a⁻¹`→`[]`, `a b a⁻¹`→`[0,1,2]`, `a⁻¹ a b`→`[2]`) | pass | FAIL (`no attribute 'free_reduce'`) |
| `test_a_walk_reading_c_cinv_is_the_empty_word` (`cd.symbols == []`, `status == "trivial"`) | pass | FAIL (`[Symbol(...), Symbol(...)] == []`) |
| `test_a_class_named_only_inside_a_cancelled_pair_is_not_added` | pass | FAIL (`[a, u] == [a]`, "Left contains one more item") |
| `test_loops_dropped_from_word_but_recorded` (now also asserts `free_reduce(...) == []`) | pass | FAIL (`no attribute 'free_reduce'`) |

**PASS.**

## Rubric 10: the evidence keeps its teeth

- On `brcopy/`, I reverted ONLY `_drive`'s two patch targets from `DG` to `SD`.
- Both `test_evidence_disagreement_names_a_stable_end` parametrizations FAIL, with
  `AttributeError: '_Crossings' object has no attribute 'iterate'`. The real `diagnostics`
  lookups run, so the patches demonstrably reach the evidence path.
- My new evidence test fails the same way. I restored with `git checkout -- .` (clean).
- The test still asserts:
  - position 1;
  - `walked == (L_1^2, L_3)` and `registered == (R_1^3, R_3^1)`;
  - one WARNING (`walk`) or none (`stable end`), counted on `DG.__name__`
    (`tanglepack.topology.diagnostics`).

**PASS.**

## Rubric 11: reduction probe (re-run)

Branch totals over all 63 items:

```
class_words 899  class_lib 0  class_oracle 0  class_occ 0  class_naive 0  class_rawlen_eq_word 899
evidence   1830  ev_lib 0     ev_oracle 0     ev_occ 0     ev_naive 0
base: class_words 899, evidence 1830, every oracle 0
CANCELLATION lines: none
```

- Every per-item column is 0 on all 63 items (table in `compare.txt`).
- 899 = the 902 classes minus the 3 unresolved singletons (`k28_series@0`, `k32_series@0`,
  `nested_outer_series@0`).
- **The oracles are not toothless.** On synthetic words, `stack_reduce` (class level and same
  occurrence) and `naive_reduce` give:
  - `a b b⁻¹ a⁻¹`: 2 / 2 / 2;
  - `a b a⁻¹`: 0 / 0 / 0;
  - `a⁻¹ a b`: 1 / 1 / 1;
  - `a a⁻¹` at different occurrences: 1 / **0** / 1.

**PASS. No stop condition.**

## Break-it checks

Each break is one edit on `brcopy/`. I ran `tests/unit/topology`,
`tests/regression/test_evidence_diagnosis.py` and `tests/plotting` (441 green before any break),
then restored with `git checkout -- .`. The tree was clean after every run (`breaks.py`).

| # | Break | Result |
|---|---|---|
| B1 | `free_reduce` cancels nothing | 4 FAIL: free_reduce unit, loops, `c c^-1` driver, cancelled-pair closure |
| B2 | closure/inertness words unreduced | 1 FAIL: cancelled-pair closure (this reader is guarded on its own) |
| B3 | `cd.symbols` unreduced | 2 FAIL: `c c^-1` driver, cancelled-pair closure |
| B4 | registered evidence word unreduced (`diagnostics:850`) | **0 FAIL on the coder's tests**; 1 FAIL after my test |
| B5 | `_drive` patches reverted to `SD` (rubric 10) | 7 FAIL, both stable-end cases among them |
| B6 | capped search with candidates reads `resolved`/`ambiguous` again | 5 FAIL: cap test, cap-on-tie, reason[truncated], truncated-no-wall, heteroclinic budget[2] |
| B7 | `unreachable` loses its reason | 3 FAIL |
| B8 | zero-crossing walk reads `ambiguous` | 3 FAIL (incl. `test_an_order_flip_gives_no_verdict`) |
| B9 | wall gate admits an ambiguous search (`!= "unresolved"`) | **0 FAIL on the coder's tests**; 1 FAIL after my test |
| B10 | the truncated reason gets the generic sentence | 3 FAIL |
| B11 | `cd.ambiguous` never set | **0 FAIL on the coder's tests**; 1 FAIL after my test |
| B12 | `plotting._class_status` returns `"-"` | **0 FAIL** (left untested, see below) |

## Tests I added (commit `65edec4`, trailer `Co-Authored-By: Claude Fable 5.1`)

Both are in `tests/unit/topology/test_symbolic_dynamics.py` and drive the real
`symbolic_dynamics` through `_drive`.

1. **`test_the_registered_word_is_reduced_before_it_is_compared`** guards plan step 4 /
   question 3 (B4).
   - The registered image reads `a u^-1 a^-1 a a^-1` and the walk reads `a u^-1 a^-1`.
   - On the branch they agree: `word_agrees is True`, `verified is True`,
     `itinerary_verified is False`, and the raw itinerary is kept.
   - It FAILS on base (`'a u^-1 a^-1 a a^-1' == 'a u^-1 a^-1'`). That failure is the intended
     behaviour change.
2. **`test_an_ambiguous_search_marks_its_class_and_is_no_wall`** guards plan step 3's wall gate
   and `cd.ambiguous` (B9, B11).
   - It uses the B6a tie searched in full: the class is `ambiguous`, `cd.ambiguous` is set,
     `is_reliable` is False, and the later class is walked with no wall.
   - It is a pin. In its old form (`not search.truncated` for `reason is None`) it passes on
     base. With base's gate broken the same way, it fails on base.

## Deviations from the plan's letter

The coder listed these, and I confirmed each one:
1. **Loops test.** It still pins `translate_itinerary` → `["a", "a^-1"]`, and adds
   `free_reduce(...) == []`. The plan said "now expects the empty word", which assumed the
   reduction lived in `translate_itinerary`. The plan's own step 2 keeps that function
   unreduced, so the coder's reading is consistent. The driver-level empty word is pinned
   separately.
2. **Extra exports** `WalkReason` and `free_reduce`, not in the plan's list.
3. **`verify_words(records, trellis, naming, table, letters)`.** Five arguments, which the plan
   allows ("may trim"). It reads `naming.homotopy`/`naming.iterated`, so `_drive`'s fake naming
   gained those two attributes.
4. **Dropped log line.** The per-class DEBUG line `"class %s: %s, word %r, verified %s"` is gone
   (log-only).
5. **`max_walks` WARNING wording.** It now ends "the search decides nothing". A capped multi-walk
   search no longer also logs the "ambiguous" WARNING.
6. **`diagnostics.py` is exactly at its 990 cap.**
7. **Rubric 5's grep** hits the `:311` annotation.
8. **Rubric 7's grep** hits `_Budget.truncated` (exempt by the rubric).

Found by me:
- **Rubric 7's text.** It says the capped tie is `ambiguous` on base; it is `unique` with
  `truncated=True`. The test fails on base either way.
- **Base `_compare_evidence` versus `verify_words`.** For a class with no evidence, base
  explicitly reset `verified`/`itinerary_verified` to None; `verify_words` `continue`s and
  leaves the dataclass defaults, which are also None. No observable difference, and the probe's
  equality covers it.
- **Coverage gaps.** Before my commit, three plan-step guards had no test: B4, B9, B11. They are
  now covered.
- **B12, `plotting._class_status`, is still untested.** I left it on purpose: it is a display
  label, no plotting test reads it, and adding one would pin wording.

## Not tested

- The failure, ambiguous and truncated rows of the walk-status mapping on real data: the corpus
  has none. They are covered by synthetic unit tests only.
- Behaviour change 3, a constrained out-of-budget search with no walk getting the budget
  sentence, at the driver level. It is pinned only through the heteroclinic budget test's
  now-unconditional `"budget" in reason`, and B10 shows that test bites.
