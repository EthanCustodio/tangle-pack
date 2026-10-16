# Plan — B6b `word-finished`

## 1. The answer

Mostly yes. The word the code computes is A.8c's word in every respect but one: adjacent inverse
pairs `c c^-1` are never cancelled. On the corpus this changes nothing. The analyst re-ran the probe
on all 63 builds and series steps under both rules (class-level, and B7's same-occurrence rule) and
found 0 cancellations in 902 class words and 1,830 evidence readings, so no refined child and no
inert verdict moves.

The algorithm module is not only the algorithm. About 240 lines of evidence comparison
(`ImageEvidence`, `EvidenceDiagnosis`, `registered_evidence`, `_diagnose`, `_order_flip`,
`_compare_evidence`) sit in `SymbolicDynamics.py`, and they are filled only after every word is
known (`:1518-1519`).

The vocabulary is the real work, because `WalkSearch.status == "trivial"` does not mean what A.8c's
"trivial" means:
- `WalkSearch` "trivial" is a walk with zero crossings, and 279 of its 402 corpus uses still give a
  one-token word.
- A.8c's (and `ClassDynamics`') "trivial" is an empty word.

## 2. Recommendation

Change it, as one slice, with four corrections to the request:
- **`trivial` keeps one meaning.** A `WalkSearch` has no "trivial". A zero-crossing walk is
  `resolved`.
- **`_member_pair` stays put.** It is refinement code (`match_members`).
- **The reduction is one function with two readers.** It runs on the class word and on the
  inertness/closure words. `translate_itinerary` stays unreduced, so the T law keeps its finer
  same-occurrence reading.
- **The `reason` field holds the six search failures.** The request's "all nine values as reasons"
  cannot hold: `unique`, `trivial` and `ambiguous` become statuses. `members_disagree` stays on the
  class (question 4).

## 3. Scope

**Slug:** `word-finished`. Branch `refactor/word-finished`, worktree `../tangle-pack-word-finished`,
base `63477e3`. The uncommitted `docs/`, `.refactor/` and `CLAUDE.md` are never touched.

**File list (ten):**

| # | File | Why |
|---|---|---|
| 1 | `src/tanglepack/topology/DualWalk.py` | statuses, `reason`, docstring |
| 2 | `src/tanglepack/topology/SymbolicDynamics.py` | `free_reduce`, the two readers, `_resolve_itinerary`, wall gate `:1477`, evidence moved out, docstring |
| 3 | `src/tanglepack/topology/diagnostics.py` | the evidence block plus `verify_words` |
| 4 | `src/tanglepack/topology/__init__.py` | re-export the moved names from `diagnostics`; export `MembersDisagree` |
| 5 | `src/tanglepack/topology/plotting.py` | `_class_status` prints the class's four-valued status |
| 6 | `tests/helpers/laws.py` | 3 status pins (`:2775, 2817, 2932`) |
| 7 | `tests/unit/topology/test_dual_walk.py` | 29 status pins and 6 `truncated` pins rewritten; the reason test |
| 8 | `tests/unit/topology/test_symbolic_dynamics.py` | evidence tests retargeted, `_drive` patches, logger name `:673`, the `a a^-1` pins `:214`/`:656`, `.truncated` `:754, 831`, new tests |
| 9 | `tests/regression/test_heteroclinic.py` | 2 status pins (`:76, 94`) |
| 10 | `tests/regression/test_evidence_diagnosis.py` | its docstring and caplog logger name go stale with the move |

**Size.** This slice is above the usual "about five files / 300 lines" bound:
- The move is roughly 240 lines out of one file and 250 into another.
- Files 6-10 are mechanical pin rewrites.
- The new logic is small: `free_reduce` (about 15 lines), the `reason` field, and `verify_words`
  (about 20 lines).

Question 6 offers a split.

**Caps (measured from `HEAD` minus the move):**
- `SymbolicDynamics.py` ≤ 1,330 (1,537 today)
- `DualWalk.py` ≤ 1,390 (1,370 today)
- `diagnostics.py` ≤ 990 (713 today)

The plan's old caps (1,100 / 1,250) cannot be reached: the block is 240 lines, not 300, and
`_member_pair` stays.

## 4. Steps

**1. Statuses (`DualWalk.py`).**

```python
WalkStatus = Literal["resolved", "ambiguous", "unresolved"]
WalkReason = Literal["unreachable", "no_start_node", "no_goal_node",
                     "no_admissible_walk", "truncated", "span_meets_no_node"]

@dataclass
class WalkSearch:
    start: ElementRef
    goal: ElementRef
    status: WalkStatus
    walks: list[Walk] = field(default_factory=list)
    start_faces: list["FaceNode"] = field(default_factory=list)
    reason: Optional[WalkReason] = None      # set iff status == "unresolved"
    constrained: bool = False
    rejected: int = 0
```

- The `truncated` field goes; `reason == "truncated"` says the same thing.
- Return sites change as follows:
  - the six failure sites return `WalkSearch(..., "unresolved", reason=...)`;
  - the zero-crossing return and `unique` become `"resolved"`;
  - `ambiguous` stays `"ambiguous"`;
  - a search whose budget ran out but kept candidates (`:1293`) returns
    `"unresolved", reason="truncated"` and keeps its `walks` for inspection.
- Callers that passed `truncated` positionally (`:1224, 1234, 1294`) switch to keyword arguments.
- `__repr__` prints the reason.

**2. One reduction, two readers (`SymbolicDynamics.py`).**

```python
def free_reduce(tokens: Sequence[tuple[BridgeClass, int]]) -> list[int]:
    """Positions surviving free reduction (A.8c): a token next to its inverse
    cancels; one stack pass reaches the fixed point."""
    kept: list[int] = []
    for i, (cls, direction) in enumerate(tokens):
        if kept and tokens[kept[-1]] == (cls, -direction):
            logger.debug("free reduction cancels %s^%+d against its inverse", cls.label, direction)
            kept.pop()
        else:
            kept.append(i)
    return kept
```

The driver applies it in two places:
- **Closure and inertness `words`** (`:1481-1487`): they read the reduced tokens. A class named only
  inside a cancelled pair is therefore not named.
- **`cd.symbols`**: the driver sets
  `cd.symbols = [s for i, s in enumerate(symbols) if i in kept]`, with
  `kept = free_reduce([(s.bridge_class, s.direction) for s in symbols])`.

`translate_itinerary` stays unreduced (loops dropped only). The T law
(`p1_itinerary_disagreements`) keeps its own same-occurrence `reduced()`, which is finer at the
`P1` level.

**3. `_resolve_itinerary` and the wall gate read the new form.**

```python
if search.status == "unresolved":
    cd.unresolved_reason = (
        "search budget reached before uniqueness was established"
        if search.reason == "truncated" else f"dual-graph walk {search.reason}")
```

- `cd.ambiguous = search.status == "ambiguous"`.
- The wall gate becomes `search.status == "resolved" and search.walk.multiplicity == 1`.
- `MembersDisagree` still sets `cd.unresolved_reason` with no `WalkSearch` (question 4).

**4. The move (`SymbolicDynamics.py` → `diagnostics.py`).**

- **What moves:** `ImageEvidence`, `EvidenceDiagnosis`, `registered_evidence`, `_diagnose`,
  `_order_flip` and `_compare_evidence`, plus the evidence bullet of the module docstring.
- **The new entry point:**
  `verify_words(dynamics_records, trellis, iterated, homotopy, naming, table, letters)`. The coder
  may trim the arguments to what the code reads. The driver calls it once, after all symbols are
  set, and it fills `evidence` / `verified` / `itinerary_verified`.
- **The registered word is compared as a word.** `verify_words` reduces it with the same
  `free_reduce` (question 3). `ImageEvidence.itinerary` stays the raw record.
- **Imports, as the analyst tested:**
  - `diagnostics` imports `_image_chain` and `trellis_itinerary` at module level.
  - It imports `translate_itinerary`, `_signature` and `free_reduce` inside the function, as
    `measured_sides` already does for `hole_point`.
  - Type names go under `TYPE_CHECKING`.
  - `SymbolicDynamics` imports `verify_words` (and `ImageEvidence` for its annotation) from
    `diagnostics`.
- **What does not change:** `_member_pair` stays in `SymbolicDynamics.py`, and `is_reliable` keeps
  `verified is not False`.
- **Exports:** `__init__.py` re-exports the moved names from `.diagnostics` and adds
  `MembersDisagree`.

**5. Docstrings.**
- The `DualWalk.py` and `SymbolicDynamics.py` module docstrings open with A.8c's two definitions:
  - the object: the reduced itinerary of `f(b)`, a class invariant;
  - the algorithm (D3): the constrained shortest walk.
- They then list what makes a class unresolved: the six reasons, `members_disagree`, and a
  singleton with no registered chain.
- The `truncated` paragraph (`DualWalk.py:94-104`) is rewritten to match.
- `plotting._class_status` returns `cd.status` (plus its `unresolved: <reason>` and `trellis`
  forms) in place of the walk status.

**6. Tests (pins rewritten, never weakened).**

Rewrites:
- `status == "unique"` and `"trivial"` → `status == "resolved"`.
- `not search.truncated` → `search.reason is None`.
- `:450` → `status == "unresolved" and reason == "truncated" and walk is None`.
- `:254-260` and `:754` → `reason == "truncated"` with the candidate kept.

The evidence tests import from `diagnostics`:
- `_drive` patches `diagnostics._image_chain` / `diagnostics.trellis_itinerary`, so it still drives
  the evidence.
- The logger name changes at `:673` and in `test_evidence_diagnosis.py`.
- Two fixtures change:
  - The loops test `:214` now expects the empty word.
  - The evidence fixture `:656` reads `a a` (swap the last two refs) in place of `a a^-1`. It keeps
    position 1 and the same `walked`/`registered` occurrences.

New tests are listed in rubric 7-10.

## 5. Behaviour

**No number changes** on any law case, deep run or series step. That covers:
- letter-free `rules` and `refined_rules`;
- both matrices;
- `is_reliable` and `verified`;
- the four class statuses;
- the closed-table size;
- `unresolved_reason` (the corpus has 3 unresolved classes, all singletons with no chain).

**What prints or reads differently:**
1. **`WalkSearch.status` values.**
   - `unique` and zero-crossing walks read `resolved`.
   - Failures read `unresolved` with a `reason`.
   - A search that was cut short but kept candidates reads `unresolved`/`truncated` (today it reads
     `unique` or `ambiguous` with `truncated=True`).
   - `WalkSearch.truncated` no longer exists.
2. **The plot's class label** shows the class status. A one-token zero-crossing class that showed
   "trivial" now shows "resolved", and "trivial" appears only for an empty word.
   `scripts/symbolic_figures.py:268` prints the new walk values (no edit).
3. **A constrained search that ran out with no walk** now gets the budget sentence as its
   `unresolved_reason`, where it read "dual-graph walk truncated". It never occurs on the corpus.
4. **Free reduction on hand-built data.** `a a^-1` is now the empty word, and a class named only
   inside a cancelled pair is no longer added to the closure. Corpus: 0 cancellations.
5. **Logger.** Evidence WARNING/INFO lines come from `tanglepack.topology.diagnostics`.
6. **Exports.** `MembersDisagree` is importable from `tanglepack.topology`.

## 6. Not doing

- W-b / W-c: the registered chain never decides a tie and never is the word.
- The `MAX_DETOUR` and `max_walks` values, and the marker and wall rules (D13, C17).
- No change to `land_element`, `shortest_walks`' search, `refine` or `match_members`.
- Moving `_member_pair`.
- B6a's second-tier span tolerance (`SCALING_RTOL`).
- Making the T law call the library reduction. It stays at the `P1` level.
- `stress/numerics/report/make_numbers.py`. It reads archived reasons, and "dual-graph walk
  unreachable" keeps its text.
- The stale `CLAUDE.md` symbolic-dynamics bullet (B12).

## 7. Rubric

1. `MPLBACKEND=Agg env/bin/python -m pytest -q` passes on the branch. `-rxX` lists exactly one
   xfail (C3) and no XPASS.
2. Every pinning test the change relies on passes on the base commit `63477e3` and on the branch.
   Pins rewritten in step 6 count as passing on base in their OLD form.
3. `git diff --stat 63477e3..refactor/word-finished` touches only the ten files of section 3.
4. **Base-vs-branch equality on all 63 items:** the 11 law cases, `p3_15`, `p3_16`,
   `p3_6_blasts`, k28 at 8 blasts, nested with 6 inner blasts, and every step of `cases.SERIES`.
   - **Equal everywhere:** letter-free `rules`, `refined_rules`, `transition_matrix()` at both
     levels, `is_reliable`, every `verified`, every `ClassDynamics.status`, every
     `unresolved_reason`, closed-table size, `inert` set.
   - **Walk statuses:** a per-class table maps the base `(status, truncated)` to the branch
     `(status, reason)` and follows exactly this mapping:
     - `unique` → `resolved`/None;
     - `trivial` → `resolved`/None;
     - `ambiguous` → `ambiguous`/None;
     - any status with `truncated=True` → `unresolved`/`truncated`;
     - each failure value → `unresolved`/that value.
   - The test engineer reruns this, starting from the `b6b/` probe scripts in the scratchpad but
     never trusting their output.
5. **The move is complete.**
   - `grep -n "EvidenceDiagnosis\|ImageEvidence\|_order_flip\|_diagnose\|_compare_evidence\|registered_evidence" src/tanglepack/topology/SymbolicDynamics.py`
     finds only the import line(s) and the one `verify_words` call.
   - `grep -n "^from .SymbolicDynamics" src/tanglepack/topology/diagnostics.py` finds nothing.
   - `_member_pair` is still defined in `SymbolicDynamics.py`.
   - `ImageEvidence.__module__ == "tanglepack.topology.diagnostics"`.
   - `from tanglepack.topology import ImageEvidence, EvidenceDiagnosis, registered_evidence, MembersDisagree`
     succeeds.
   - `import tanglepack` succeeds in a fresh interpreter.
6. **Line counts (`wc -l`):** `SymbolicDynamics.py` ≤ 1,330, `DualWalk.py` ≤ 1,390,
   `diagnostics.py` ≤ 990. The report states all three numbers.
7. **Statuses are three.**
   - `grep -rn '"unique"\|"trivial"' src/tanglepack/topology/DualWalk.py` finds no status literal.
   - `grep -rn "search\.truncated\|\.truncated and\|not .*\.truncated" src tests` finds nothing
     (`_Budget.truncated` is internal and stays).
   - A parametrized test over the hand-built dual graphs in `test_dual_walk.py` reaches each of the
     six reasons with `status == "unresolved"`.
   - A test on the `:254` tie at `max_walks=1` reads `status == "unresolved"`,
     `reason == "truncated"` and `walk is not None`. It fails on base, where the status is
     `ambiguous`.
8. **Free reduction.**
   - A unit test shows `free_reduce` over `a, b, b^-1, a^-1` keeps no position, and over
     `a, b, a^-1` keeps all three.
   - A driver test with a walked itinerary reading `c c^-1` gives `cd.symbols == []` and
     `cd.status == "trivial"`. It fails on base.
9. **Both readers reduce.** A driver test where a class's only walked tokens are `c c^-1`, with `c`
   outside the alphabet, does NOT add `c` to the closed table. It fails on base, where the closure
   adds it.
10. **The evidence keeps its teeth.** `test_evidence_disagreement_names_a_stable_end` still asserts
    position 1, the same `walked`/`registered` occurrences, and one WARNING (counted on the
    `diagnostics` logger) for the "walk" cause. Its `_drive` patches demonstrably reach the evidence
    path: the test fails if the patch target is reverted to `SD._image_chain`.
11. **The free-reduction probe** is re-run on the branch with the library's `free_reduce` on all 63
    items and shows 0 cancellations in class words and in evidence. Same-occurrence counts stay 0.
    It is attached to the coder's report and re-run by the test engineer. Any non-zero count STOPS
    the slice: report the class, the word before and after, and whether a refined child changes.

## 8. Questions

1. **How should a decided zero-crossing walk be named?** Default: `WalkSearch.status` has three
   values, `resolved` / `ambiguous` / `unresolved`. "trivial" belongs only to a class, meaning an
   empty word. (The alternative is to keep a fourth walk status and document the two meanings.)
2. **What happens to a search cut short with candidates?** Default: it becomes `unresolved` with
   reason `truncated`, its `walks` are kept for inspection, and the `truncated` field goes.
3. **Which reduction rule, and who uses it?** Default:
   - Class-level per A.8c, as one `free_reduce` applied to `cd.symbols` and to the
     closure/inertness words.
   - `translate_itinerary` stays unreduced; the T law keeps its same-occurrence `reduced()`.
   - The registered evidence word is reduced by the same function before comparison, so a word is
     compared with a word. The evidence fixture becomes `a a`.

   Alternative: the evidence stays unreduced, as the orchestrator suggested. That keeps the fixture,
   but a cancellation in a walk word would then read as a "walk" disagreement and make the table
   unreliable.
4. **Where does `members_disagree` live?** Default: on the class only (`unresolved_reason`, as
   today). No `WalkSearch` is built for a search that never ran, so the walk reasons are six.
5. **Export `MembersDisagree` from `tanglepack.topology`?** Default: yes. It is one line, and the
   file is in the list anyway.
6. **Keep this as one slice?** Default: yes. The alternative is two runs: statuses first, then the
   move and the reduction.
