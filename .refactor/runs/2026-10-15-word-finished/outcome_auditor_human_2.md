# Outcome — B6b `word-finished` (auditor)

**Author:** please grade this "confirmed", "touchup" or "figures". The branch is `refactor/word-finished`
@ `65edec4` (worktree `../tangle-pack-word-finished`), base `63477e3`. It is NOT merged; the merge is yours.
Both worktrees are clean. On main, `.refactor/decisions.md` and `.refactor/knowledge/{topology,testing}.md`
are updated and left uncommitted, as you set it up. `CLAUDE.md` and `docs/` were not touched.

## Your question, answered now

**Is the word of a class what A.8c defines?** Yes.
- The word is now the freely reduced itinerary of `f(b)`: loops are dropped and an adjacent `c c^-1`
  cancels at class level (`SymbolicDynamics.free_reduce`, `:545`).
- The same function reduces all three words: the class word (`cd.symbols`), the words the closure and
  inertness read, and the registered evidence word.
- `translate_itinerary` stays unreduced, so the T law keeps its finer same-occurrence reading.
- On the 63-item corpus nothing cancels, under either rule. So no word, refined child, inert verdict
  or matrix moved.

**Is the algorithm module only the algorithm?** Yes, for the evidence.
- `ImageEvidence`, `EvidenceDiagnosis`, `registered_evidence`, `_diagnose` and `_order_flip` now live
  in `topology/diagnostics.py`. `_compare_evidence` became `verify_words` there.
- The algorithm calls `verify_words` once, after every word is final (`SymbolicDynamics.py:1301`).
- `_member_pair` stays in `SymbolicDynamics.py`, because it is refinement code.

**The vocabulary is A.8c's.**
- A walk search is `resolved`, `ambiguous` or `unresolved`, with one of six reasons.
- "trivial" now means only one thing: a class whose word is empty.
- Before this change, `trivial` was also a walk status meaning "zero crossings". On k=2.8 two blasts,
  class `b` (word `c`) had walk status `trivial` and class status `resolved`. That double meaning is gone.

## What changed

There are 10 files, exactly the plan's list (`git diff --stat`: +658 / −435), in four commits:
three by the coder and one test commit by the test engineer.

| File | Before | After |
|---|---|---|
| `SymbolicDynamics.py` | 1,537 | **1,312** (cap 1,330) |
| `DualWalk.py` | 1,370 | **1,383** (cap 1,390; grew by the A.8c docstring and the `reason` field) |
| `diagnostics.py` | 713 | **990** (cap 990: exactly at it) |

The most important hunk is how a search becomes a class verdict, in `_resolve_itinerary`.

Before (`63477e3`):
```python
walk = search.walk
if walk is None or search.status not in ("unique", "trivial", "ambiguous"):
    cd.unresolved_reason = f"dual-graph walk {search.status}"
elif search.truncated:  # a budget cut short decides nothing, whatever the status
    cd.unresolved_reason = "search budget reached before uniqueness was established"
```
After (`65edec4`):
```python
if search.status == "unresolved":
    cd.unresolved_reason = (
        "search budget reached before uniqueness was established"
        if search.reason == "truncated"
        else f"dual-graph walk {search.reason}"
    )
```

What the search itself now returns (`shortest_walks`):
- zero crossings or one itinerary: `resolved`;
- several itineraries: `ambiguous`;
- the six failures: `unresolved` plus `reason`;
- a budget cut that still has candidates: `unresolved`/`truncated`, with the walks kept for
  inspection. The `WalkSearch.truncated` field is gone.

The wall gate now reads `status == "resolved" and multiplicity == 1`. Both module docstrings open with
A.8c's two definitions (the object, and the D3 algorithm) and the list of what leaves a class unresolved.

## Every behaviour change

No number changed on any of the 63 items. What reads or acts differently:

1. **`WalkSearch.status` values.**
   - Old `unique` and `trivial` both read `resolved`.
   - Failures read `unresolved` plus a `reason`.
   - A capped search reads `unresolved`/`truncated` and keeps its walks. On base it read `unique` or
     `ambiguous` with `truncated=True`.
   - `scripts/symbolic_figures.py:268` prints the new values without being edited.
2. **The itinerary-table status label (`plotting._class_status`).**
   - It now shows the class status, so a one-token zero-crossing class that used to show "trivial"
     shows "resolved".
   - The `-` fallback is gone.
3. **A constrained search that runs out of budget with no walk** now gets the budget sentence as its
   `unresolved_reason` (base: "dual-graph walk truncated"). This never happens on the corpus.
4. **Free reduction.**
   - A walk reading `a a^-1` is the empty word, so its class is `trivial`.
   - A class named only inside a cancelled pair is not added to the closure.
   - A registered word is reduced before it is compared.
   - The corpus has none of these.
5. **Logging.**
   - Evidence WARNING/INFO lines now come from `tanglepack.topology.diagnostics`.
   - The per-class DEBUG line (`class %s: %s, word %r, verified %s`) is gone.
   - A capped search with several walks logs only the max_walks WARNING (now ending "the search decides
     nothing"), not also the ambiguous one.
6. **Exports.** `MembersDisagree` is exported as you approved. The coder also exported `WalkReason` and
   `free_reduce`, which you did not ask for.

## What was verified, and how

**What I ran myself:**
- **Full suite on `65edec4`:** 1223 passed, 47 skipped, 5 deselected, 1 xfailed (C3), no XPASS, 204 s.
  Base is 1211; the 12 new tests account for the difference.
- **The diff.** I read the whole `src` diff. The `tests` diff I read for removed assertions: every
  removed pin has a same-strength or stricter replacement (`in ("unique","trivial")` became
  `== "resolved"`, which now also fails a capped search).
- **The rubric greps.**
  - No `"unique"`/`"trivial"` literal is left in `DualWalk.py`.
  - The only `.truncated` left in src/tests/scripts is `_Budget.truncated`, which is internal.
  - `SymbolicDynamics.py` names the evidence types only in the `TYPE_CHECKING` import (`:188`) and
    the `ClassDynamics.evidence` annotation (`:311`).
  - The line counts are the ones in the table above.
- **Imports in a fresh interpreter (cwd `/tmp`):**
  - `from tanglepack.topology import ImageEvidence, EvidenceDiagnosis, registered_evidence, MembersDisagree, free_reduce, WalkReason`
    works.
  - The moved names report module `tanglepack.topology.diagnostics`.
- **My own base-vs-branch spot check** (`scratchpad/auditor/spot.py`) on k10 and k28 two blasts:
  - Identical words (`a -> a u^-1 b^-1`, `b -> c`, `c -> a u^-1 a^-1`; k10 `a -> a u^-1 a^-1`),
    class statuses, `verified` and `is_reliable`.
  - The only diff is the walk status (`unique`/`trivial` became `resolved`).

**Re-run by the test engineer and the reviewer** (I read their scripts and outputs; I did not re-run all 63):
- **Base vs branch on all 63 items** (16 builds plus every step of the five series), equal on every
  field: letter-free `rules` and `refined_rules`, both matrices, `is_reliable`, every `verified` and
  `itinerary_verified`, every class status and `unresolved_reason`, the closed-table size, the inert
  set and the evidence words.
- **The walk-status mapping** covers 553 searches: 402 `trivial` and 151 `unique`, all now `resolved`.
  There are no violations, and 349 classes have no search.
- **Free-reduction probe:** 0 cancellations in 899 class words and 1,830 evidence readings, by the
  library function and by three independent oracles (class level, same occurrence, naive). The
  oracles were shown to fire on synthetic words.
- **Twelve break-its.** B1–B11 each make at least one test fail. Before the TE's commit, three plan
  steps had no guard:
  - the evidence word left unreduced (B4);
  - the wall gate admitting an ambiguous search (B9);
  - `cd.ambiguous` never set (B11).

  The TE's two tests close all three. B12 (the plot label) has no test.

## What did not go to plan

- **The plan's rubric 7 was wrong about base.** The capped tie at `max_walks=1` reads `unique` with
  `truncated=True` on base, not `ambiguous`. The new test still fails on base, which is what the
  item required.
- **The loops unit test.** The plan said it would "now expect the empty word". That would contradict
  its own step 2, which keeps `translate_itinerary` unreduced. So the test still pins
  `translate_itinerary -> ["a","a^-1"]` and adds `free_reduce(...) == []`. A separate driver test
  pins the empty class word.
- **Smaller deviations.**
  - The two extra exports, the dropped DEBUG line and the changed WARNING (see the behaviour list).
  - `_compare_evidence` was folded into `verify_words`, which takes five arguments rather than the
    plan's seven (the plan allowed the trim).
  - The rubric-5 grep also hits the `ClassDynamics.evidence` annotation, which plan step 4 itself
    asked for.
- **`diagnostics.py` reached its cap only after squashing.** It landed at exactly 990 lines. To get
  there the coder folded `_compare_evidence` and compacted two log calls (one logger call puts three
  arguments on a line, and `verdicts` uses the name `e`). The reviewer flagged both as off the
  module's idiom but readable. Any later addition to that file must remove a line.

## Left open for you

- **`CLAUDE.md`.** The symbolic-dynamics bullet and the `DualWalk`/`SymbolicDynamics` table rows still
  describe the nine statuses, `WalkSearch.truncated`, and the evidence living in `SymbolicDynamics`
  (B12, deliberately untouched).
- **`free_reduce` runs twice per class.** Once on the closure tokens (built from `oriented_class` of
  the homotopy parents) and once on `cd.symbols` (built from `translate_itinerary`). If a cancellation
  ever fires, its DEBUG line prints twice. These two derivations of one word predate B6b.
- **`plotting._class_status` is untested.** That is deliberate: a test would pin wording.
- **The ambiguous, failure and truncated rows of the status mapping are exercised only by hand-built
  unit tests.** The same holds for any real cancellation: the corpus has none of these.
- **Commit trailers.** All four commits end `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`,
  as B7b's did. The attribution in force now names `Claude Opus 5.5`. Amend before merging if that
  matters to you.

## Where I am least sure this is what you meant

1. **A capped search keeps its walks but is `unresolved`.** This is the default you approved. But
   `search.walks` is now non-empty on an unresolved search. Any future code that reads `search.walk`
   without first checking `status` will pick up an undecided candidate.
   - Today only `_resolve_itinerary`, which checks the status first, and the wall gate (which needs
     `resolved`) read it.
2. **The registered evidence word is reduced before it is compared.** This is question 3's default.
   The consequence: a registered chain reading `w c c^-1` now VERIFIES a walked `w`, and only
   `itinerary_verified` (element level, INFO) records the difference. The TE's test pins exactly this
   (`verified True`, `itinerary_verified False`). If you wanted the evidence to flag such a
   back-and-forth, that is the alternative you declined.
3. **The plot label changed.** A zero-crossing class with a one-token word used to show "trivial" in
   the itinerary table and now shows "resolved". This is the vocabulary decision working as intended.
   But it changes what readers of the existing figures see, and it is the one user-visible change no
   test guards.
