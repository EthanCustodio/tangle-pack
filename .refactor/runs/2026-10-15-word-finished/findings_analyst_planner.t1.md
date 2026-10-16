# Findings — B6b `word-finished` (analyst, iteration 1)

**The working tree is dirty.** `HEAD` = `63477e3` (B7b merged). There are uncommitted changes to
`.refactor/decisions.md`, `.refactor/knowledge/*.md`, `CLAUDE.md`, `docs/alignment_plan_2026-10-09.md`
and `docs/code_fidelity_audit_2026-10-09.md`, plus untracked run folders, figures and two sweep scripts.
The change branches from `HEAD` and must not touch any of them. This run updated
`.refactor/knowledge/topology.md` and `.refactor/knowledge/testing.md` (B6b sections at the end).

**Planner:** you read this, and the person running the session then grades your plan "approved",
"revise" or "answered". Every item under "What differs from the request" below is a decision or a
correction your plan has to deal with.

## Kind of request and target

This is a refactor: vocabulary, a move and one reduction rule. The targets are:

- `DualWalk.WalkStatus` / `WalkSearch`;
- `SymbolicDynamics._resolve_itinerary`, `translate_itinerary` and the evidence block;
- `topology/diagnostics.py` as the destination of the move.

I read the request the way it is written: re-derive the plan's numbers against `HEAD`, count the
consumers, check the import graph, run the probe under both rules, and do not re-plan.

## Short answer

1. **The word is A.8c's word except for the reduction, and the reduction is a no-op on the corpus
   under either rule** (verified). Neither rule cancels anything on any of the 63 builds and series
   steps:
   - (a) class-level adjacent `c c^-1`: 0 cancellations;
   - (b) same-occurrence, as in `laws.py` `reduced()`: 0 cancellations.

   Both were checked over all 902 class words (active and inert, walked and trellis-sourced) and
   all 1,830 registered evidence readings. No refined child changes and no inert verdict changes.
   The two rules DO differ on hand-built unit fixtures, though: three unit tests fail under
   class-level reduction (details below).
2. **The algorithm module is not only the algorithm** (verified). The plan's six evidence names take
   223 lines of definitions inside `SymbolicDynamics.py`. The algorithm never reads them: they are
   filled at `:1518-1519` after the words.
3. **`_member_pair` is NOT evidence and must not move** (verified). The request adds it to the plan's
   list, but it is called by `SymbolicDynamics.match_members` (`:852`, the refinement), and the test
   `_drive` helper monkeypatches it (`test_symbolic_dynamics.py:567`).

## 1. Line numbers at HEAD (verified by reading)

| Cited | At HEAD | Note |
|---|---|---|
| `DualWalk.py:557-567` `WalkStatus` | 557-567 | ✓ nine values |
| `WalkSearch :575-634` | 575-634 | ✓ no `reason` field; `truncated` is a separate flag |
| status returns `:1109,1118,1123,1153,1164,1224,1234,1280-1294` | same | ✓ no_start_node, no_goal_node, span_meets_no_node, trivial, unreachable, truncated, no_admissible_walk, trivial/unique/ambiguous |
| `MembersDisagree :790`, `trellis_itinerary :1301` | same | ✓ |
| `SymbolicDynamics.py` docstring `:1-130` | `:1-154` | the evidence bullet is `:80-99` |
| `EvidenceDiagnosis :268`, `ImageEvidence :300` | 267-296, 299-333 (decorator lines first) | ✓ |
| `ClassDynamics :337-414`, `.status :405-414` | 336-413, 404-413 | ✓ (±1) |
| `translate_itinerary :540-612` | 540-610 | ✓ |
| `_signature :613`, `refine :621-670` | 613, 621-668 | ✓ |
| `registered_evidence :687`, `_member_pair :739` | 687-736, 739-745 | ✓ |
| `is_reliable :900`, `describe :1070` | 899, 1070 | ✓ |
| `_resolve_itinerary :1171-1266`, mapping `:1249-1266` | 1171-1265, mapping 1250-1253, `cd.ambiguous` 1265 | ✓ |
| `_diagnose :1268`, `_order_flip :1299`, `_compare_evidence :1323` | same (end 1296, 1320, 1379) | ✓ |
| `symbolic_dynamics :1382`; "the driver's lines 137-138"; "translates at `:134`" | 1382-1537; the calls are at ABSOLUTE 1518-1519; the translate is at absolute 1515 | the request's `:134` is a driver-relative line, not a file line |
| `diagnostics.py` 713 lines | 713 | ✓ |
| `__init__.py:108-116`, `:179-187` | same | ✓ |
| `plotting.py:1810-1817` | `_class_status` 1807-1817 | ✓ |
| `laws.py:2775, 2817, 2932`; `:2781`; `reduced() :3038`; `:78` | same | ✓ |
| `p1_itinerary_disagreements :3009-3040` / `:2998-3080` | def 3000-3082 | small drift |
| `test_heteroclinic.py:76,94` | same | ✓ |
| `test_symbolic_dynamics.py:443-445` | same | these are `ClassDynamics.status`, which already has the four values; a WalkSearch collapse does not touch them |
| `test_dual_walk.py :88-530` | 88-530 | 29 assertion lines (see §2) |

## 2. Consumers and pins (verified by grep)

### The nine WalkSearch statuses and `search.status`

**src**

- `DualWalk.py`: the literal (`:557-567`), the docstrings (`:68`, `:94-104`, `:582-595`,
  `:1075-1098`), the eight return sites, and `__repr__` `:631`.
- `SymbolicDynamics.py`:
  - `:1250-1251` maps the statuses into `unresolved_reason` (`f"dual-graph walk {search.status}"`).
  - `:1265` sets `cd.ambiguous`.
  - **`:1477` is the wall gate `search.status in ("unique", "trivial")`, which the request does not
    list.**
  - docstring `:112-117`.
- `plotting.py:1817` `_class_status` (its docstring at `:1808` names the values).
- `topology/__init__.py` exports `WalkStatus` (`:97`, `:170`).
- `describe()` does **not** print `search.status`: it prints `cd.status` (`:1091`), which is already
  four-valued.

**tests (pins)**

- `test_dual_walk.py`: 29 assertion lines, made up of `unique` 12, `ambiguous` 9, `unreachable` 2,
  and one each of `trivial`, `no_start_node`, `no_goal_node`, `no_admissible_walk`, `truncated` and
  `span_meets_no_node`.
  - The truncated-flag pins are `:254-260` (status `ambiguous` WITH `truncated=True`) and `:450`.
- `test_heteroclinic.py:76`, `:94`.
- `laws.py:2775`, `:2817`, `:2932`.
- `test_symbolic_dynamics.py:754` pins `first.search.truncated and first.search.walk is not None`
  (the candidate is kept under a truncated budget).

**Outside the six files, print-only**

- `scripts/symbolic_figures.py:268` puts `search.status` into its health report. It keeps working
  and prints the new values.
- `stress/numerics/report/make_numbers.py:894-906` matches `"unreachable"` inside ARCHIVED
  `unresolved_reason` strings. Archived data, out of scope.

### The seven evidence names (src, tests, scripts, `__init__`)

- `EvidenceDiagnosis`, `ImageEvidence`, `registered_evidence`:
  - `__init__.py` `:108-116`, `:179-187`;
  - SD internals;
  - `test_symbolic_dynamics.py` `:470` and `:496` call `SD.registered_evidence`, after
    monkeypatching `SD._image_chain` / `SD.trellis_itinerary` (`:466-467`, `:492-493`);
  - `test_evidence_diagnosis.py:5` mentions the old path in its docstring only.
- `_member_pair`: `SD.match_members :852` and the `_drive` monkeypatch `:567`.
- `_diagnose`, `_order_flip`, `_compare_evidence`: SD only.
  - `test_evidence_diagnosis.py:53` defines its OWN `_order_flip` (an independent oracle), which is
    not a consumer.
- No script, no `stress/` file and nothing in `loom` imports any of the seven.

### Honest file list for a status collapse plus the evidence move

1. `src/tanglepack/topology/DualWalk.py`
2. `src/tanglepack/topology/SymbolicDynamics.py`
3. `src/tanglepack/topology/diagnostics.py`
4. `src/tanglepack/topology/__init__.py`
5. `src/tanglepack/topology/plotting.py`: `_class_status` and its docstring. It is required only
   if the label should read the new form; `str(status)` already works.
6. `tests/helpers/laws.py`: 3 pins. The T law also needs a decision; see §4.
7. `tests/unit/topology/test_dual_walk.py`: 29 pins, plus the new reachability test.
8. `tests/unit/topology/test_symbolic_dynamics.py`:
   - the evidence tests `:448` and `:483` (retarget to `diagnostics`);
   - `_drive` `:567-569`: its `SD._image_chain` / `SD.trellis_itinerary` patches will no longer
     reach a `registered_evidence` that lives in `diagnostics`, which would silently stop driving
     the evidence of `test_evidence_disagreement_names_a_stable_end`;
   - the logger name in `:673`;
   - the two fixtures in §4;
   - the new tests.
9. `tests/regression/test_heteroclinic.py`: 2 pins.
10. Optional: `tests/regression/test_evidence_diagnosis.py`. Its docstring (`:5`) and its caplog
    logger name (`:88`, `:120`, `"tanglepack.topology.SymbolicDynamics"`) go stale. I *believe* it
    still passes, because it counts WARNINGs by message text, and WARNING records from any logger
    reach caplog. I did not verify this with the move applied.

That is nine files, ten with the optional one. This matches the request's "nine or ten".

## 3. Import graph for the move (verified on a scratch copy of `src/`)

At call time, `diagnostics` needs from the algorithm side:

- `BridgeClass._image_chain` and `DualWalk.trellis_itinerary`;
- `SymbolicDynamics.translate_itinerary` and `_signature`.

It needs these only as types: `Symbol`, `ClassDynamics`, `ElementNaming`, `BridgeClass*`,
`ElementRef`.

I tested on a copy (`scratchpad/b6b/imp/`): `diagnostics` imports `_image_chain` and
`trellis_itinerary` at module level, and SD imports from `diagnostics` at module level. `import`
of `tanglepack`, `.topology`, `.diagnostics`, `.SymbolicDynamics`, `.StablePartition`, `.DualWalk`,
`.Trellis` and `tanglepack.loom` all succeed. Adding `from .SymbolicDynamics import
translate_itinerary, _signature` at module level in `diagnostics` gives an
`ImportError (partially initialized module)` from every entry point, including `import tanglepack`.

Cleanest route:

- SD imports `verify_words` (and `ImageEvidence` for the annotation) from `diagnostics` at module
  level.
- `diagnostics` imports `_image_chain` and `trellis_itinerary` at module level, and imports
  `translate_itinerary` / `_signature` INSIDE `verify_words` / `registered_evidence`.
  `diagnostics.measured_sides` already does this for `plotting.hole_point` (`diagnostics.py:458`).
- Type names go under `TYPE_CHECKING`.
- `ClassDynamics.evidence: list[ImageEvidence]` is a string annotation (`from __future__ import
  annotations`), so it needs only a `TYPE_CHECKING` import.
- Passing `translate` in as an argument also works, but it adds a parameter for no gain.

## 4. Free-reduction probe at HEAD (verified)

**Scripts:** `scratchpad/b6b/probe.py`, `run.sh`, `summ.py`, output in `out/*.jsonl`.

**What it builds:** the 16 builds (the 11 law cases, `p3_15`, `p3_16`, `p3_6_blasts`, `k28` at 8
blasts, and nested with 6 inner blasts) plus every step of the five `cases.SERIES`
(9 + 10 + 5 + 11 + 12 = 47). Each series step is blasted exactly as `cases.blast_series` does.

**What it does per item:**

- It takes every `ClassDynamics` with an itinerary and freely reduces `cd.symbols` (loops already
  dropped) with a stack, which reaches the fixed point.
- Rule (a) cancels a token followed by the same class in the opposite direction. Rule (b)
  additionally requires the same `occurrence`, which is B7's `reduced()`.
- It does the same for every `ImageEvidence.symbols`.
- It recomputes `refine` and `inert_classes` on the class-level-reduced words and compares them
  with the unreduced ones.

**Result:** all 63 items read 0 cancellations under both rules, for class words and evidence alike.

- `refine` is unchanged on 63/63 and `inert_classes` is unchanged on 63/63.
- The recomputed inert set equals `cd.kind` on 63/63, which checks the harness.

**Cross-check through the real test suite.** The scratch plugin `b6b_reduce_plugin.py` patches a
reduction into `SD.translate_itinerary` before the test modules import it, so the driver, the
evidence and the T law (`p1_itinerary_disagreements`) all go through it.

| Run | Result | Reduction hits |
|---|---|---|
| HEAD | 1211 passed, 47 skipped, 5 deselected, 1 xfailed (C3), 214 s | — |
| same-occurrence rule | identical | 0 |
| class-level rule | 3 failed, 1208 passed | 3, every one a hand-built `a a^-1` |

The three class-level failures:

- `test_loops_dropped_from_word_but_recorded` (`test_symbolic_dynamics.py:214-218`) pins
  `["a", "a^-1"]` after dropping a loop. The occurrences are `(R_1^1, R_3^3)` and `(R_1^3, R_3^1)`:
  different. **This existing pin directly contradicts A.8c's class-level rule**, and the plan's new
  "reduces to the empty word" test would replace it.
- `test_evidence_disagreement_names_a_stable_end[signs0-stable end-run0]` and `[signs1-walk-run1]`
  (`:645-675`): the registered fixture is `a a^-1` (`:656`), so the diagnosis position goes from
  1 to 0. To keep its teeth the fixture needs a non-cancelling registered word (for example
  `a u`), or the evidence path must read unreduced symbols.

**Per-item counts**, columns as printed. Every class-level and same-occurrence column is 0.
WalkSearch statuses are given for context:

| build / step | classes | class-level cancellations (words) | same-occurrence (words) | evidence readings | class-level (evidence) | same-occurrence (evidence) | WalkSearch statuses |
|---|---|---|---|---|---|---|---|
| heteroclinic | 75 | 0 | 0 | 70 | 0 | 0 | trivial 46, unique 6 |
| inversion | 4 | 0 | 0 | 4 | 0 | 0 | trivial 1 |
| inversion_one_blast | 7 | 0 | 0 | 7 | 0 | 0 | trivial 5, unique 2 |
| inversion_period2 | 6 | 0 | 0 | 6 | 0 | 0 | trivial 1 |
| inversion_two_blasts | 8 | 0 | 0 | 10 | 0 | 0 | trivial 6, unique 2 |
| k10 | 2 | 0 | 0 | 3 | 0 | 0 | trivial 1, unique 1 |
| k28_8 | 4 | 0 | 0 | 59 | 0 | 0 | trivial 2, unique 2 |
| k28_one_blast | 3 | 0 | 0 | 3 | 0 | 0 | trivial 1 |
| k28_two_blasts | 4 | 0 | 0 | 4 | 0 | 0 | trivial 2, unique 2 |
| nested | 9 | 0 | 0 | 9 | 0 | 0 | trivial 3, unique 2 |
| nested_inner6 | 47 | 0 | 0 | 25 | 0 | 0 | trivial 23, unique 4 |
| nested_inner_blasts | 9 | 0 | 0 | 9 | 0 | 0 | trivial 3, unique 2 |
| p3 | 5 | 0 | 0 | 5 | 0 | 0 | trivial 1 |
| p3_15 | 43 | 0 | 0 | 21 | 0 | 0 | trivial 21, unique 2 |
| p3_16 | 43 | 0 | 0 | 30 | 0 | 0 | trivial 21, unique 2 |
| p3_6_blasts | 43 | 0 | 0 | 21 | 0 | 0 | trivial 21, unique 2 |
| inversion_series@0 | 4 | 0 | 0 | 4 | 0 | 0 | trivial 1 |
| inversion_series@1 | 7 | 0 | 0 | 7 | 0 | 0 | trivial 5, unique 2 |
| inversion_series@2 | 8 | 0 | 0 | 10 | 0 | 0 | trivial 6, unique 2 |
| inversion_series@3 | 8 | 0 | 0 | 13 | 0 | 0 | trivial 6, unique 2 |
| inversion_series@4 | 8 | 0 | 0 | 19 | 0 | 0 | trivial 6, unique 2 |
| k28_series@0 | 3 | 0 | 0 | 1 | 0 | 0 | trivial 1 |
| k28_series@1 | 3 | 0 | 0 | 3 | 0 | 0 | trivial 1 |
| k28_series@2 | 4 | 0 | 0 | 4 | 0 | 0 | trivial 2, unique 2 |
| k28_series@3 | 4 | 0 | 0 | 7 | 0 | 0 | trivial 2, unique 2 |
| k28_series@4 | 4 | 0 | 0 | 13 | 0 | 0 | trivial 2, unique 2 |
| k28_series@5 | 4 | 0 | 0 | 21 | 0 | 0 | trivial 2, unique 2 |
| k28_series@6 | 4 | 0 | 0 | 35 | 0 | 0 | trivial 2, unique 2 |
| k28_series@7 | 4 | 0 | 0 | 59 | 0 | 0 | trivial 2, unique 2 |
| k28_series@8 | 4 | 0 | 0 | 59 | 0 | 0 | trivial 2, unique 2 |
| k32_series@0 | 3 | 0 | 0 | 1 | 0 | 0 | trivial 1 |
| k32_series@1 | 3 | 0 | 0 | 3 | 0 | 0 | trivial 1 |
| k32_series@2 | 9 | 0 | 0 | 4 | 0 | 0 | trivial 3, unique 3 |
| k32_series@3 | 9 | 0 | 0 | 9 | 0 | 0 | trivial 3, unique 3 |
| k32_series@4 | 9 | 0 | 0 | 15 | 0 | 0 | trivial 3, unique 3 |
| k32_series@5 | 9 | 0 | 0 | 23 | 0 | 0 | trivial 3, unique 3 |
| k32_series@6 | 18 | 0 | 0 | 40 | 0 | 0 | trivial 9, unique 9 |
| k32_series@7 | 18 | 0 | 0 | 71 | 0 | 0 | trivial 9, unique 9 |
| k32_series@8 | 45 | 0 | 0 | 124 | 0 | 0 | trivial 16, unique 14 |
| k32_series@9 | 45 | 0 | 0 | 124 | 0 | 0 | trivial 16, unique 14 |
| k32_series@10 | 45 | 0 | 0 | 124 | 0 | 0 | trivial 16, unique 14 |
| nested_outer_series@0 | 8 | 0 | 0 | 6 | 0 | 0 | trivial 2 |
| nested_outer_series@1 | 8 | 0 | 0 | 8 | 0 | 0 | trivial 2 |
| nested_outer_series@2 | 9 | 0 | 0 | 9 | 0 | 0 | trivial 3, unique 2 |
| nested_outer_series@3 | 9 | 0 | 0 | 12 | 0 | 0 | trivial 3, unique 2 |
| nested_outer_series@4 | 9 | 0 | 0 | 18 | 0 | 0 | trivial 3, unique 2 |
| nested_outer_series@5 | 9 | 0 | 0 | 26 | 0 | 0 | trivial 3, unique 2 |
| nested_outer_series@6 | 9 | 0 | 0 | 40 | 0 | 0 | trivial 3, unique 2 |
| nested_outer_series@7 | 9 | 0 | 0 | 66 | 0 | 0 | trivial 3, unique 2 |
| nested_outer_series@8 | 9 | 0 | 0 | 108 | 0 | 0 | trivial 3, unique 2 |
| nested_outer_series@9 | 9 | 0 | 0 | 108 | 0 | 0 | trivial 3, unique 2 |
| nested_outer_series@10 | 9 | 0 | 0 | 108 | 0 | 0 | trivial 3, unique 2 |
| nested_outer_series@11 | 9 | 0 | 0 | 108 | 0 | 0 | trivial 3, unique 2 |
| p3_series@0 | 5 | 0 | 0 | 5 | 0 | 0 | trivial 1 |
| p3_series@1 | 5 | 0 | 0 | 5 | 0 | 0 | trivial 1 |
| p3_series@2 | 5 | 0 | 0 | 5 | 0 | 0 | trivial 1 |
| p3_series@3 | 5 | 0 | 0 | 5 | 0 | 0 | trivial 1 |
| p3_series@4 | 5 | 0 | 0 | 5 | 0 | 0 | trivial 1 |
| p3_series@5 | 5 | 0 | 0 | 5 | 0 | 0 | trivial 1 |
| p3_series@6 | 43 | 0 | 0 | 21 | 0 | 0 | trivial 21, unique 2 |
| p3_series@7 | 43 | 0 | 0 | 26 | 0 | 0 | trivial 21, unique 2 |
| p3_series@8 | 43 | 0 | 0 | 28 | 0 | 0 | trivial 21, unique 2 |
| p3_series@9 | 43 | 0 | 0 | 29 | 0 | 0 | trivial 21, unique 2 |

Totals: 902 class words, 1,830 evidence readings, 0 / 0 / 0 / 0 cancellations. 553 WalkSearches:
`unique` 151 and `trivial` 402. No other status occurs, and the `truncated` flag is never set at
`max_walks=64`. The 3 unresolved classes (`k28_series@0`, `k32_series@0`, `nested_outer_series@0`)
are singleton landings with no registered chain (`cd.search` None; reason checked on `k28_series@0`).
Every fine failure status is therefore reachable only in hand-built unit tests.

**The 2026-10-09 probe.** The only record is the plan's sentence (`alignment_plan_2026-10-09.md:1074`:
"0 adjacent `c c^-1` pairs in any class word on 13 builds"). I found no script and no output in
`docs/`, `.refactor/` or the scratchpad. B7's analyst record
(`.refactor/runs/2026-10-14-refined-level-defined/findings_analyst_planner.t1.md:66`) reports that
the same-occurrence reduction "never fired" on member readings, over 15 builds and 47 steps. Today's
result is consistent with both, and it extends them to 63 items, both rules, and class words as
well as evidence.

## 5. Line counts and caps (measured)

| File | HEAD | After a pure move (measured) | Honest cap |
|---|---|---|---|
| `SymbolicDynamics.py` | 1,537 | **1,297** | about 1,330 |
| `DualWalk.py` | 1,370 | nothing moves out | about 1,390 |
| `diagnostics.py` | 713 | +about 240 | not set here |

How the numbers come out:

- **SD after a pure move.** I deleted the six definitions (223 lines) with their separator lines,
  the `# ── registered evidence ──` header (`:684-686`) and the two driver calls (`:1518-1519`),
  keeping `_member_pair`. That leaves 1,297 lines. The cap adds the import, the `verify_words`
  call, the reduction loop and the A.8c docstring, minus the 20-line evidence bullet moved to
  `diagnostics`.
- **DualWalk.** The literal shrinks by about 5 lines; the `reason` field and the A.8c paragraph add
  some. The cap is 1,370 plus the A.8c paragraph.
- **The plan's caps (SD ≤ 1,100, DualWalk ≤ 1,250) are unreachable by this slice.** The evidence
  block is about 240 lines, not the plan's "~300", and `_member_pair` stays.

## 6. What differs from the request (for the planner to resolve)

1. **`_member_pair` stays.** It is refinement code (`match_members :852`). The plan's own step 3
   list (`alignment_plan:1076-1077`) does not include it; the request added it.
2. **WalkSearch "trivial" is not A.8c's "trivial".** `WalkSearch.status == "trivial"` means a
   zero-crossing walk (start and goal share a face), whose itinerary `(start, goal)` is usually ONE
   token. On the corpus, 279 of the 402 trivial walks give a non-empty word. `ClassDynamics.status
   == "trivial"` means an EMPTY word (151 classes). Collapsing to `{resolved, trivial, ambiguous,
   unresolved}` with today's "trivial" puts two meanings on one word. The wall gate `:1477` and
   `laws.py:2932` read "trivial" as "decided, multiplicity 1", not as "empty word".
3. **"All nine fine values reachable as `reason`" (rubric 4) cannot hold.** `unique`, `trivial` and
   `ambiguous` become statuses, not reasons. The reasons are the six failure values (`unreachable`,
   `no_start_node`, `no_goal_node`, `no_admissible_walk`, `truncated`, `span_meets_no_node`), plus
   `members_disagree`.
4. **`members_disagree` is not a WalkSearch outcome today.** `MembersDisagree` is raised before
   `shortest_walks` runs, so `cd.search` stays None and no WalkSearch exists to carry the reason.
   Either `_resolve_itinerary` builds an unresolved WalkSearch for it, or the reason lives only on
   the class.
5. **There are two kinds of "truncated".** In the first, the status is `"truncated"` with no walk
   (`:1224`, constrained search). In the second, the `truncated=True` flag sits on a `unique` or
   `ambiguous` search that keeps its candidate (`:1293`; B6a turns it unresolved in
   `_resolve_itinerary :1252`). Pins `test_dual_walk.py:254-260` and `test_symbolic_dynamics.py:754`
   rely on the candidate surviving. The plan must say whether the second kind becomes status
   `unresolved` with reason `truncated` (with or without keeping `walks`).
6. **The reduction has two readers, not one.** The class-level `words` behind inertness (B5) and
   the closure are built from RAW pairs (`:1481-1487`), not from `translate_itinerary`. Reducing only
   in `translate_itinerary` would make `cd.symbols` and the inertness/closure graph disagree as soon
   as a cancellation fires. On the corpus nothing fires, so no number moves either way.
7. **The default reconciliation needs an unreduced path.** `registered_evidence` and the T law
   (`laws.py:3070`) both call `translate_itinerary`. If the library cancels at class level
   INSIDE `translate_itinerary`, the T law receives class-reduced symbols. Its same-occurrence
   `reduced()` can then never fire, and T silently becomes the COARSER class-level comparison. The
   default "the T law keeps its `P1`-level reduction" therefore needs one of two things: the
   reduction as a separate function the driver calls, or a flag on `translate_itinerary`. The
   evidence words face the same question: today's evidence-diagnosis fixture is `a a^-1`.
8. **The two rules differ on the hand-built fixtures (3 tests) but never on the corpus.** The
   request's "they cannot differ while there are 0 class-level cancellations" holds for the corpus
   only.
9. **`describe()` needs no change for the statuses.** It prints `cd.status`, which is already
   four-valued. The visible change is `plotting._class_status` (and the script's health line).
10. **The plan section says "eight" statuses.** That predates B6a's `span_meets_no_node`; there
    are nine at HEAD, as the request says.
11. **Logger name.** The evidence WARNING/INFO move from `tanglepack.topology.SymbolicDynamics` to
    the `diagnostics` logger, if `diagnostics`' own logger is used. `test_symbolic_dynamics.py:673`
    counts by that logger name.

## Recommendation (smallest change; the planner owns the plan)

- Collapse the WalkSearch statuses. Resolve items 2-5 explicitly. The smallest consistent choice
  that keeps the meaning: keep a decided zero-crossing walk out of the word "trivial", or document
  that WalkSearch's `trivial` is zero crossings while the class's `trivial` is an empty word.
- Move the six evidence definitions, and NOT `_member_pair`, behind `diagnostics.verify_words`,
  with the call-time import of `translate_itinerary` / `_signature`.
- Add free reduction as ONE function over `(class, direction)` tokens, used by the driver for
  `cd.symbols` AND for the inertness/closure `words`.
- Leave `translate_itinerary` unreduced so the T law and the evidence keep their own reading, or
  give it an explicit flag.
- Caps: SD ≤ about 1,330, DualWalk ≤ about 1,390.
- Corpus behaviour: no number changes (verified for the reduction; the move and the status rename
  are mechanical).

## Out of scope (noticed)

- `CLAUDE.md`'s symbolic-dynamics bullet still lists the old status set and "falls back to all of
  them" (B12).
- `MembersDisagree` is still not exported (B6a's open point; it is a one-liner if the status work
  touches it).
- `stress/numerics/report/make_numbers.py` matches `"unreachable"` in archived reason strings.
  Fine as long as the new reason text still contains the fine value.
